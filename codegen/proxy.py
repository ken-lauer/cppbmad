from __future__ import annotations

import logging
import pathlib
from dataclasses import dataclass
from string import Template
from typing import TYPE_CHECKING, Any

from .paths import CODEGEN_ROOT, CPPBMAD_ROOT, CPPBMAD_SRC
from .types import STANDARD_TYPES, ArgumentType, FullType, PointerType, native_type_containers
from .util import struct_to_proxy_class_name

if TYPE_CHECKING:
    from .arg import Argument
    from .gen import CodegenStructure

logger = logging.getLogger(__name__)


@dataclass
class TemplateEntry:
    fortran_getter: str
    fortran_setter: str | None
    cpp_get_decl: str
    cpp_get_accessors: list[str]
    cpp_set_decl: str | None
    cpp_set_accessors: list[str]


# Helper map for native allocatable types (e.g. 'real' -> 'RealAlloc1D')
# matches nt.name.lower() to the argument type key
_native_alloc_map = {nt.name.lower(): nt.cpp_container_name for nt in native_type_containers}


def get_condition(cpptype: str, logic_type: str, attr_name: str) -> str:
    target = f"struct_obj%{attr_name}"

    if logic_type == "NOT":
        base = ".true."
    elif logic_type == "ALLOC":
        base = f"allocated({target})"
    elif logic_type == "PTR":
        base = f"associated({target})"
    else:
        raise ValueError(f"Unknown logic_type: {logic_type}")

    # Arrays passed to C/C++ via c_loc must be contiguous otherwise the
    # calculated strides in C++ won't match memory gaps.
    # Note: 'is_contiguous' is standard F2008.
    if cpptype == "string":
        return base
    return f"{base} .and. is_contiguous({target})"


# ---------------------------------------------------------------------------
# Dynamic Fortran Code Generator for Arrays
# ---------------------------------------------------------------------------
def generate_fortran_array_routine(
    cpp_type: str, dim: int, is_derived_type: bool = False, logic_type: str = "NOT"
) -> str:
    """
    Generates the Fortran implementation string dynamically for 1D, 2D, 3D arrays,
    returns valid views (data + bounds + strides).
    """

    args = [
        "struct_obj_ptr",
        "data_ptr",
        "bounds",
        "strides" if dim > 1 else None,
        "is_allocated",
        "el_size" if is_derived_type else None,
    ]
    args = [a for a in args if a]  # Filter None
    arg_str = ", ".join(args)
    n_bounds = dim * 2

    decls = [
        "type(c_ptr), intent(in), value :: struct_obj_ptr",
        "type(c_ptr), intent(out) :: data_ptr",
        f"integer(c_int), dimension({n_bounds}), intent(out) :: bounds",
        "logical(c_bool), intent(out) :: is_allocated",
        "type(STRUCTNAME), pointer :: struct_obj",
    ]

    if dim > 1:
        decls.append(f"integer(c_int), dimension({dim}), intent(out) :: strides")
        # Helper ints for stride calc
        decls.append("integer :: " + ", ".join(f"d{i}" for i in range(1, dim)))

    if is_derived_type:
        decls.append("integer(c_size_t), intent(out) :: el_size")

    decls_str = "\n    ".join(decls)

    # Condition Block
    condition = get_condition(cpp_type, logic_type, "FATTRNAME")

    # Access Pattern (c_loc args)
    loc_args = ", ".join(f"lbound(struct_obj%FATTRNAME, {i})" for i in range(1, dim + 1))
    ref_args = ", ".join(f"bounds({(i - 1) * 2 + 1})" for i in range(1, dim + 1))

    bound_assigns = []
    for i in range(1, dim + 1):
        idx_l = (i - 1) * 2 + 1
        idx_u = (i - 1) * 2 + 2
        bound_assigns.append(f"bounds({idx_l}) = int(lbound(struct_obj%FATTRNAME, {i}), c_int)")
        bound_assigns.append(f"bounds({idx_u}) = int(ubound(struct_obj%FATTRNAME, {i}), c_int)")

    bound_str = "\n      ".join(bound_assigns)

    stride_assigns = []
    if dim > 1:
        stride_assigns.append("strides(1) = 1_c_int")
        calc_dim = []
        calc_dim.append("d1 = bounds(2) - bounds(1) + 1")
        if dim > 2:
            calc_dim.append("d2 = bounds(4) - bounds(3) + 1")
        stride_assigns.extend(calc_dim)
        stride_assigns.append("strides(2) = d1")
        if dim > 2:
            stride_assigns.append("strides(3) = d1 * d2")

    stride_str = "\n      ".join(stride_assigns)

    el_size_logic = ""
    if is_derived_type:
        el_size_logic = f"el_size = int(storage_size(struct_obj%FATTRNAME({ref_args})) / 8, c_size_t)"

    zero_assigns = ["data_ptr = c_null_ptr", "bounds = 0_c_int"]
    if dim > 1:
        zero_assigns.append("strides = 0_c_int")
    if is_derived_type:
        zero_assigns.append("el_size = 0")
    zero_assigns.append("is_allocated = .false.")
    zero_str = "\n      ".join(zero_assigns)

    return f"""
  subroutine STRUCTNAME_get_FATTRNAME_info({arg_str}) &
        bind(c, name='STRUCTNAME_get_FATTRNAME_info')
    {decls_str}

    call c_f_pointer(struct_obj_ptr, struct_obj)

    if ({condition}) then
      data_ptr = c_loc(struct_obj%FATTRNAME({loc_args}))
      {bound_str}
      {stride_str}
      {el_size_logic}
      is_allocated = .true.
    else
      {zero_str}
    endif
  end subroutine
"""


def generate_fortran_alloc_routine(
    cpp_type: str, dim: int, is_derived_type: bool = False, logic_type: str = "ALLOC"
) -> str:
    """
    Generates a simplified Fortran declaration for Allocatable arrays (specifically 1D)
    that returns just the data pointer and allocation status (and element size for types),
    skipping bounds and stride calculations, allowing the C++ Proxy class to manage logic.
    """
    if dim != 1:
        # Fallback to standard routine if multi-dim alloc support isn't simplified yet
        return generate_fortran_array_routine(cpp_type, dim, is_derived_type, logic_type)

    # Generate standard getter (now that we unified interfaces)
    dummy = generate_fortran_array_routine(cpp_type, dim, is_derived_type, logic_type)

    # Append Realloc routine
    realloc = """
  subroutine STRUCTNAME_reallocate_FATTRNAME(struct_obj_ptr, lbound_, n) &
        bind(c, name='STRUCTNAME_reallocate_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    integer(c_int), value :: lbound_
    integer(c_size_t), value :: n
    type(STRUCTNAME), pointer :: struct_obj

    call c_f_pointer(struct_obj_ptr, struct_obj)
    
    if (n == 0) then
      if (allocated(struct_obj%FATTRNAME)) deallocate(struct_obj%FATTRNAME)
    else
      if (allocated(struct_obj%FATTRNAME)) deallocate(struct_obj%FATTRNAME)
      allocate(struct_obj%FATTRNAME(lbound_:lbound_ + n - 1))
    endif
  end subroutine
"""
    return dummy + realloc


# ---------------------------------------------------------------------------
# Dynamic Fortran Code Generator for Array setters
# ---------------------------------------------------------------------------
def generate_fortran_array_setter(cpp_type: str, dim: int, logic_type: str = "NOT") -> str:
    """
    Generates Fortran implementation for setting array attributes from a flat C buffer.
    Supports 1D, 2D, 3D.
    """
    is_boolean = cpp_type == "bool"
    fortran_val_type = "integer(c_int)" if is_boolean else "FORTRANTYPE"

    decls = [
        "type(c_ptr), intent(in), value :: struct_obj_ptr",
        "type(c_ptr), intent(in), value :: val_ptr",
        f"integer(c_int), dimension({dim}), intent(in) :: shape",
        "type(STRUCTNAME), pointer :: struct_obj",
    ]

    colons = ",".join([":"] * dim)
    decls.append(f"{fortran_val_type}, pointer :: val({colons})")

    decls_str = "\n    ".join(decls)

    body = ["call c_f_pointer(struct_obj_ptr, struct_obj)"]
    body.append("if (c_associated(val_ptr)) then")
    body.append("  call c_f_pointer(val_ptr, val, shape)")

    target = "struct_obj%FATTRNAME"

    if is_boolean:
        rhs = "(val .ne. 0)"
    else:
        rhs = "val"

    if logic_type == "ALLOC":
        realloc_check = " .or. ".join(f"(size({target}, {d}) /= shape({d}))" for d in range(1, dim + 1))
        body.append(f"  if (allocated({target})) then")
        body.append(f"     if ({realloc_check}) deallocate({target})")
        body.append("  endif")
        body.append(
            f"  if (.not. allocated({target})) allocate({target}({', '.join([f'shape({d})' for d in range(1, dim + 1)])}))"
        )
        body.append(f"  {target} = {rhs}")

    elif logic_type == "PTR":
        size_check = " .and. ".join(f"(size({target}, {d}) == shape({d}))" for d in range(1, dim + 1))
        body.append(f"  if (associated({target})) then")
        body.append(f"    if ({size_check}) then")
        body.append(f"       {target} = {rhs}")
        body.append("    endif")
        body.append("  endif")

    else:
        body.append(f"  {target} = {rhs}")

    body.append("endif")

    body_str = "\n    ".join(body)

    return f"""
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, val_ptr, shape) &
      bind(c, name='STRUCTNAME_set_FATTRNAME')
    {decls_str}

    {body_str}
  end subroutine
"""


# ---------------------------------------------------------------------------
# Fortran pattern fragments - SCALARS
# ---------------------------------------------------------------------------
FORTRAN_SCALAR_GETTER = """
  subroutine STRUCTNAME_get_FATTRNAME(struct_obj_ptr, value_out) bind(c, name='STRUCTNAME_get_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    FORTRANTYPE, intent(out) :: value_out
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    value_out = struct_obj%FATTRNAME
  end subroutine
"""

# Gets a value from a pointer safely (returns value + validity boolean)
# Used for Optional Returns in C++
FORTRAN_SCALAR_PTR_GETTER_OPTIONAL = """
  subroutine STRUCTNAME_get_FATTRNAME(struct_obj_ptr, value_out, is_valid) bind(c, name='STRUCTNAME_get_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    FORTRANTYPE, intent(out) :: value_out
    logical(c_bool), intent(out) :: is_valid
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)

    if (associated(struct_obj%FATTRNAME)) then
      value_out = struct_obj%FATTRNAME ! Implicit cast/copy
      is_valid = .true.
    else
      is_valid = .false.
    endif
  end subroutine
"""

FORTRAN_POINTER_GETTER = """
  subroutine STRUCTNAME_get_FATTRNAME(struct_obj_ptr, ptr_out) bind(c, name='STRUCTNAME_get_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(out) :: ptr_out
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    if (associated(struct_obj%FATTRNAME)) then
      ptr_out = c_loc(struct_obj%FATTRNAME)
    else
      ptr_out = c_null_ptr
    endif
  end subroutine
"""

# ---------------------------------------------------------------------------
# C++ pattern fragments - GETTERS
# ---------------------------------------------------------------------------

CPP_SCALAR_DECL = "    void STRUCTNAME_get_FATTRNAME(const void* struct_obj, CTYPE* value_out);"
CPP_SCALAR_ACCESSOR = """
    CTYPE CATTRNAME() const {
        CTYPE value;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &value);
        return value;
    }
"""

CPP_SCALAR_PTR_OPTIONAL_DECL = (
    "    void STRUCTNAME_get_FATTRNAME(const void* struct_obj, CTYPE* value_out, bool* is_valid);"
)
CPP_SCALAR_PTR_OPTIONAL_ACCESSOR = """
    std::optional<CTYPE> CATTRNAME() const {
        CTYPE value;
        bool is_valid;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &value, &is_valid);
        if (is_valid) return value;
        return std::nullopt;
    }
"""

CPP_POINTER_DECL = "    void STRUCTNAME_get_FATTRNAME(const void* struct_obj, CTYPE** ptr_out);"
CPP_POINTER_ACCESSOR = """
    CTYPE* CATTRNAME() const {
        CTYPE* ptr;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &ptr);
        return ptr;
    }
"""

CPP_TYPE_SCALAR_DECL = "    void STRUCTNAME_get_FATTRNAME(const void* struct_obj, void** ptr_out);"

# ---------------------------------------------------------------------------
# Fortran/C++ pattern fragments - SETTERS
# ---------------------------------------------------------------------------
FORTRAN_SCALAR_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, value_in) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    FORTRANTYPE, intent(in), value :: value_in
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    struct_obj%FATTRNAME = value_in
  end subroutine
"""

FORTRAN_POINTER_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, value_in) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    FORTRANTYPE, intent(in), value :: value_in
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    if (associated(struct_obj%FATTRNAME)) then
      struct_obj%FATTRNAME = value_in
    endif
  end subroutine
"""
FORTRAN_TYPE_SCALAR_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, src_ptr) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(in), value :: src_ptr
    type(STRUCTNAME), pointer :: struct_obj
    type(ATTRTYPE), pointer :: src_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    call c_f_pointer(src_ptr, src_obj)
    struct_obj%FATTRNAME = src_obj
  end subroutine
"""
FORTRAN_TYPE_POINTER_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, src_ptr) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(in), value :: src_ptr
    type(STRUCTNAME), pointer :: struct_obj
    type(ATTRTYPE), pointer :: src_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    if (associated(struct_obj%FATTRNAME)) then
      call c_f_pointer(src_ptr, src_obj)
      struct_obj%FATTRNAME = src_obj
    endif
  end subroutine
"""

CPP_SCALAR_SET_DECL = "    void STRUCTNAME_set_FATTRNAME(void* struct_obj, CTYPE value_in);"
CPP_SCALAR_SET_ACCESSOR = """
    void set_CATTRNAME(CTYPE value) {
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, value);
    }
"""
CPP_POINTER_SET_DECL = CPP_SCALAR_SET_DECL
CPP_POINTER_SET_ACCESSOR = CPP_SCALAR_SET_ACCESSOR

CPP_TYPE_SCALAR_SET_DECL = "    void STRUCTNAME_set_FATTRNAME(void* struct_obj, const void* src_ptr);"
CPP_TYPE_SCALAR_SET_ACCESSOR = """
    void set_CATTRNAME(const ${return_proxy_name}& src) {
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, src.get_fortran_ptr());
    }
"""
CPP_TYPE_POINTER_SET_DECL = CPP_TYPE_SCALAR_SET_DECL
CPP_TYPE_POINTER_SET_ACCESSOR = CPP_TYPE_SCALAR_SET_ACCESSOR

# ---------------------------------------------------------------------------
# Character scalars / Array 1D
# ---------------------------------------------------------------------------
FORTRAN_CHAR_SCALAR_DYN_GETTER = """
  subroutine STRUCTNAME_get_FATTRNAME_info(struct_obj_ptr, data_ptr, str_len, is_allocated) &
    bind(c, name='STRUCTNAME_get_FATTRNAME_info')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(out) :: data_ptr
    integer(c_int), intent(out) :: str_len
    logical(c_bool), intent(out) :: is_allocated
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)

    if (CONDITION) then
      data_ptr = c_loc(struct_obj%FATTRNAME)
      str_len = int(len(struct_obj%FATTRNAME), c_int)
      is_allocated = .true.
    else
      data_ptr = c_null_ptr
      str_len = 0
      is_allocated = .false.
    endif
  end subroutine
"""

FORTRAN_CHAR_ALLOC_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, str_ptr, str_len) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(in), value :: str_ptr
    integer(c_int), intent(in), value :: str_len
    type(STRUCTNAME), pointer :: struct_obj
    character(len=str_len), pointer :: temp_str

    call c_f_pointer(struct_obj_ptr, struct_obj)

    if (allocated(struct_obj%FATTRNAME)) deallocate(struct_obj%FATTRNAME)

    if (str_len > 0) then
       call c_f_pointer(str_ptr, temp_str)
       allocate(struct_obj%FATTRNAME, source=temp_str)
       struct_obj%FATTRNAME = temp_str(1:str_len)
    endif
  end subroutine
"""
FORTRAN_CHAR_PTR_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, str_ptr, str_len) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(in), value :: str_ptr
    integer(c_int), intent(in), value :: str_len
    type(STRUCTNAME), pointer :: struct_obj
    character(len=str_len), pointer :: temp_str

    call c_f_pointer(struct_obj_ptr, struct_obj)

    if (associated(struct_obj%FATTRNAME)) deallocate(struct_obj%FATTRNAME)

    if (str_len > 0) then
        call c_f_pointer(str_ptr, temp_str)
        allocate(struct_obj%FATTRNAME, source=temp_str)
    else
        nullify(struct_obj%FATTRNAME)
    endif
  end subroutine
"""

# Character Arrays are distinct because they return length AND bounds
FORTRAN_CHAR_ARRAY_1D_ALL = """
  subroutine STRUCTNAME_get_FATTRNAME_info(struct_obj_ptr, data_ptr, bounds, str_len, is_allocated) &
      bind(c, name='STRUCTNAME_get_FATTRNAME_info')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(out) :: data_ptr
    integer(c_int), dimension(2), intent(out) :: bounds
    integer(c_int), intent(out) :: str_len
    logical(c_bool), intent(out) :: is_allocated
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)

    if (CONDITION) then
      data_ptr = c_loc(struct_obj%FATTRNAME(lbound(struct_obj%FATTRNAME, 1)))
      bounds(1) = int(lbound(struct_obj%FATTRNAME, 1), c_int)
      bounds(2) = int(ubound(struct_obj%FATTRNAME, 1), c_int)
      str_len = int(len(struct_obj%FATTRNAME), c_int)
      is_allocated = .true.
    else
      data_ptr = c_null_ptr
      bounds = 0
      str_len = 0
      is_allocated = .false.
    endif
  end subroutine
"""

CPP_CHAR_SCALAR_DYN_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(
        const void* s,
        char** d,
        int* len,
        bool* is_alloc
    );
"""

CPP_CHAR_SCALAR_DYN_ACCESSOR = """
    std::string CATTRNAME() const {
        return ProxyHelpers::get_string(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
    }
"""
CPP_CHAR_ARRAY_1D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(
        const void* s,
        char** d,
        int* bounds,    // [lower, upper]
        int* str_len,
        bool* is_alloc
    );
"""

CPP_CHAR_ARRAY_1D_ACCESSOR = """
    FCharArray1D CATTRNAME() const {
        return ProxyHelpers::get_char_array_1d(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
    }
"""


def subst(s: str, **kw) -> str:
    out = s
    for k, v in kw.items():
        out = out.replace(k.upper(), v)
    return out


# ---------------------------------------------------------------------------
# Template Makers
# ---------------------------------------------------------------------------
def make_scalar(fortran_type: str, cpp_type: str) -> TemplateEntry:
    return TemplateEntry(
        fortran_getter=subst(FORTRAN_SCALAR_GETTER, fortrantype=fortran_type),
        fortran_setter=subst(FORTRAN_SCALAR_SETTER, fortrantype=fortran_type),
        cpp_get_decl=subst(CPP_SCALAR_DECL, ctype=cpp_type),
        cpp_get_accessors=[subst(CPP_SCALAR_ACCESSOR, ctype=cpp_type)],
        cpp_set_decl=subst(CPP_SCALAR_SET_DECL, ctype=cpp_type),
        cpp_set_accessors=[subst(CPP_SCALAR_SET_ACCESSOR, ctype=cpp_type)],
    )


def make_scalar_pointer_optional(fortran_type: str, cpp_type: str) -> TemplateEntry:
    return TemplateEntry(
        fortran_getter=subst(FORTRAN_SCALAR_PTR_GETTER_OPTIONAL, fortrantype=fortran_type),
        fortran_setter=subst(FORTRAN_POINTER_SETTER, fortrantype=fortran_type),
        cpp_get_decl=subst(CPP_SCALAR_PTR_OPTIONAL_DECL, ctype=cpp_type),
        cpp_get_accessors=[subst(CPP_SCALAR_PTR_OPTIONAL_ACCESSOR, ctype=cpp_type)],
        cpp_set_decl=subst(CPP_POINTER_SET_DECL, ctype=cpp_type),
        cpp_set_accessors=[subst(CPP_POINTER_SET_ACCESSOR, ctype=cpp_type)],
    )


def make_char_array_1d(ptr_type: PointerType) -> TemplateEntry:
    return TemplateEntry(
        fortran_getter=subst(
            FORTRAN_CHAR_ARRAY_1D_ALL, condition=get_condition("string", ptr_type, "FATTRNAME")
        ),
        fortran_setter=None,
        cpp_get_decl=CPP_CHAR_ARRAY_1D_DECL,
        cpp_get_accessors=[CPP_CHAR_ARRAY_1D_ACCESSOR],
        cpp_set_decl=None,
        cpp_set_accessors=[],
    )


def make_array_template(
    dim: int,
    ptr_type: str,
    cpp_type: str,
    cpp_decl_tmpl: str,
    cpp_acc_tmpl: str,
    cpp_set_decl_tmpl: str = "",
    cpp_set_acc_tmpl: str = "",
    is_derived_type: bool = False,
    fortran_gen_func: Any = None,
) -> TemplateEntry:
    """
    Generic generator for 1D, 2D, 3D arrays of both Simple and Derived types.
    """
    fortran_s = None
    cpp_s_decl = None
    cpp_s_acc = []

    if cpp_set_decl_tmpl and cpp_set_acc_tmpl:
        fortran_s = generate_fortran_array_setter(cpp_type, dim, logic_type=ptr_type)
        cpp_s_decl = subst(cpp_set_decl_tmpl, ctype=cpp_type)
        cpp_s_acc = [subst(cpp_set_acc_tmpl, ctype=cpp_type)]

    gen_func = fortran_gen_func or generate_fortran_array_routine

    return TemplateEntry(
        fortran_getter=gen_func(cpp_type, dim, is_derived_type, ptr_type),
        fortran_setter=fortran_s,
        cpp_get_decl=subst(cpp_decl_tmpl, ctype=cpp_type),
        cpp_get_accessors=[subst(cpp_acc_tmpl, ctype=cpp_type)],
        cpp_set_decl=cpp_s_decl,
        cpp_set_accessors=cpp_s_acc,
    )


def make_char_scalar_dyn(ptr_type: str, setter_pattern: str) -> TemplateEntry:
    return TemplateEntry(
        fortran_getter=subst(
            FORTRAN_CHAR_SCALAR_DYN_GETTER, condition=get_condition("string", ptr_type, "FATTRNAME")
        ),
        fortran_setter=setter_pattern,
        cpp_get_decl=CPP_CHAR_SCALAR_DYN_DECL,
        cpp_get_accessors=[CPP_CHAR_SCALAR_DYN_ACCESSOR],
        cpp_set_decl=CPP_CHARACTER_SET_DECL,
        cpp_set_accessors=[CPP_CHARACTER_SET_ACCESSOR],
    )


# ---------------------------------------------------------------------------
# Build templates dictionary
# ---------------------------------------------------------------------------
templates: dict[FullType, TemplateEntry] = {}

# Scalar simple types
_simple_types: list[ArgumentType] = ["real", "real16", "integer", "integer8", "logical"]
for tname in _simple_types:
    info = STANDARD_TYPES[tname]
    templates[FullType(tname, 0, "NOT")] = make_scalar(info.fortran_type, info.c_type)
    templates[FullType(tname, 0, "PTR")] = make_scalar_pointer_optional(info.fortran_type, info.c_type)

# Complex Scalar (custom accessors needed)
info = STANDARD_TYPES["complex"]
templates[FullType("complex", 0, "NOT")] = TemplateEntry(
    fortran_getter=subst(FORTRAN_SCALAR_GETTER, fortrantype=info.fortran_type),
    fortran_setter=subst(FORTRAN_SCALAR_SETTER, fortrantype=info.fortran_type),
    cpp_get_decl=subst(CPP_SCALAR_DECL, ctype=info.c_type),
    cpp_get_accessors=[
        """
    std::complex<double> CATTRNAME() const {
        std::complex<double> c_value;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &c_value);
        return c_value;
    }
"""
    ],
    cpp_set_decl=subst(CPP_SCALAR_SET_DECL, ctype=info.c_type),
    cpp_set_accessors=[
        """
    void set_CATTRNAME(std::complex<double> value) {
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, value);
    }
"""
    ],
)

# Complex Pointer scalar
templates[FullType("complex", 0, "PTR")] = TemplateEntry(
    fortran_getter=subst(FORTRAN_SCALAR_PTR_GETTER_OPTIONAL, fortrantype=info.fortran_type),
    fortran_setter=subst(FORTRAN_POINTER_SETTER, fortrantype=info.fortran_type),
    cpp_get_decl="    void STRUCTNAME_get_FATTRNAME(const void* struct_obj, double _Complex* val_out, bool* is_valid);",
    cpp_get_accessors=[
        """
    std::optional<std::complex<double>> CATTRNAME() const {
        std::complex<double> val;
        bool is_valid;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, reinterpret_cast<double _Complex*>(&val), &is_valid);
        if(is_valid) return val;
        return std::nullopt;
    }
"""
    ],
    cpp_set_decl="    void STRUCTNAME_set_FATTRNAME(void* struct_obj, std::complex<double> value_in);",
    cpp_set_accessors=[
        """
    void set_CATTRNAME(std::complex<double> value) {
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, value);
    }
"""
    ],
)

# ---------------------------------------------------------------------------
# Simplified Character Handling
# ---------------------------------------------------------------------------
FORTRAN_CHARACTER_SIMPLE_SETTER = """
  subroutine STRUCTNAME_set_FATTRNAME(struct_obj_ptr, str_ptr, str_len) bind(c, name='STRUCTNAME_set_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(in), value :: str_ptr
    integer(c_int), intent(in), value :: str_len
    type(STRUCTNAME), pointer :: struct_obj
    character(len=str_len), pointer :: str_in
    call c_f_pointer(struct_obj_ptr, struct_obj)
    call c_f_pointer(str_ptr, str_in)
    struct_obj%FATTRNAME = str_in ! implicitly handles padding
  end subroutine
"""
CPP_CHARACTER_SET_DECL = (
    "    void STRUCTNAME_set_FATTRNAME(void* struct_obj, const char* str_ptr, int str_len);"
)
CPP_CHARACTER_SET_ACCESSOR = """
    void set_CATTRNAME(const std::string& value) {
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
    }
"""

# Character scalar (NOT)
templates[FullType("character", 0, "NOT")] = TemplateEntry(
    fortran_getter="""
  subroutine STRUCTNAME_get_FATTRNAME_info(struct_obj_ptr, data_ptr, bounds, is_allocated) &
    bind(c, name='STRUCTNAME_get_FATTRNAME_info')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(out) :: data_ptr
    integer(c_int), dimension(2), intent(out) :: bounds
    logical(c_bool), intent(out) :: is_allocated
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    data_ptr = c_loc(struct_obj%FATTRNAME)
    bounds(1) = 1_c_int
    bounds(2) = int(len_trim(struct_obj%FATTRNAME), c_int)
    is_allocated = .true.
  end subroutine
""",
    fortran_setter=FORTRAN_CHARACTER_SIMPLE_SETTER,
    cpp_get_decl="void STRUCTNAME_get_FATTRNAME_info(const void* s, char** d, int* bounds, bool* a);",
    cpp_get_accessors=[
        """
    std::string CATTRNAME() const {
        FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
        return std::string(arr.data(), arr.size());
    }
"""
    ],
    cpp_set_decl=CPP_CHARACTER_SET_DECL,
    cpp_set_accessors=[CPP_CHARACTER_SET_ACCESSOR],
)

templates[FullType("character", 0, "ALLOC")] = make_char_scalar_dyn("ALLOC", FORTRAN_CHAR_ALLOC_SETTER)
templates[FullType("character", 0, "PTR")] = make_char_scalar_dyn("PTR", FORTRAN_CHAR_PTR_SETTER)

# ---------------------------------------------------------------------------
# Arrays 1D / 2D / 3D (Real + Integer + Complex)
# ---------------------------------------------------------------------------

# Templates for C++ side logic
CPP_ARRAY_1D_DECL = (
    "    void STRUCTNAME_get_FATTRNAME_info(const void* s, CTYPE** d, int* bounds, bool* is_alloc);"
)
CPP_ARRAY_1D_ACCESSOR = """
    FArray1D<CTYPE> CATTRNAME() const {
        return ProxyHelpers::get_array_1d<CTYPE>(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
    }
"""

# Alloc 1D - Primitives
CPP_ALLOC_1D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(const void* s, void** d, int* bounds, bool* is_alloc);
    void STRUCTNAME_reallocate_FATTRNAME(void* s, int lb, size_t n);
"""
# Accessor returns standard FAlloc1D<T, ...> with struct-specific function pointers
CPP_ARRAY_1D_ALLOC_ACCESSOR = """
    ALLOC_TYPE CATTRNAME() const {
        return ALLOC_TYPE(
            const_cast<void*>(fortran_ptr_),
            STRUCTNAME_reallocate_FATTRNAME,
            STRUCTNAME_get_FATTRNAME_info
        );
    }
"""

CPP_ARRAY_SET_DECL = "    void STRUCTNAME_set_FATTRNAME(void* s, const void* d, const int* shape);"

# Generic 1D
CPP_ARRAY_1D_SET_ACCESSOR_GENERIC = """
    void set_CATTRNAME(const std::vector<CTYPE>& v) {
        int shape[] = { static_cast<int>(v.size()) };
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, v.data(), shape);
    }
"""
# Bool 1D
CPP_ARRAY_1D_SET_ACCESSOR_BOOL = """
    void set_CATTRNAME(const std::vector<bool>& v) {
        int shape[] = { static_cast<int>(v.size()) };
        std::vector<int> bv(v.size());
        for(size_t i=0; i<v.size(); ++i) bv[i] = v[i] ? 1 : 0;
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, bv.data(), shape);
    }
"""

CPP_ARRAY_2D_DECL = "    void STRUCTNAME_get_FATTRNAME_info(const void* s, CTYPE** d, int* bounds, int* strides, bool* is_alloc);"
CPP_ARRAY_2D_ACCESSOR = """
    FArray2D<CTYPE> CATTRNAME() const {
        return ProxyHelpers::get_array_2d<CTYPE>(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
    }
"""

# Generic 2D
CPP_ARRAY_2D_SET_ACCESSOR_GENERIC = """
    void set_CATTRNAME(const std::vector<std::vector<CTYPE>>& v) {
        ProxyHelpers::set_array_2d<CTYPE>(fortran_ptr_, STRUCTNAME_set_FATTRNAME, v);
    }
"""
# Bool 2D
CPP_ARRAY_2D_SET_ACCESSOR_BOOL = """
    void set_CATTRNAME(const std::vector<std::vector<bool>>& v) {
        int rows = static_cast<int>(v.size());
        int cols = rows > 0 ? static_cast<int>(v[0].size()) : 0;
        int shape[] = { cols, rows }; 
        
        std::vector<int> flat;
        flat.reserve(rows * cols);
        for (int j = 0; j < cols; ++j) {
            for (int i = 0; i < rows; ++i) {
                flat.push_back(v[i][j] ? 1 : 0);
            }
        }
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, flat.data(), shape);
    }
"""

CPP_ARRAY_3D_DECL = "    void STRUCTNAME_get_FATTRNAME_info(const void* s, CTYPE** d, int* bounds, int* strides, bool* is_alloc);"
CPP_ARRAY_3D_ACCESSOR = """
    FArray3D<CTYPE> CATTRNAME() const {
        return ProxyHelpers::get_array_3d<CTYPE>(fortran_ptr_, STRUCTNAME_get_FATTRNAME_info);
    }
"""

# Generic 3D
CPP_ARRAY_3D_SET_ACCESSOR_GENERIC = """
    void set_CATTRNAME(const std::vector<std::vector<std::vector<CTYPE>>>& v) {
        ProxyHelpers::set_array_3d<CTYPE>(fortran_ptr_, STRUCTNAME_set_FATTRNAME, v);
    }
"""

# Bool 3D
CPP_ARRAY_3D_SET_ACCESSOR_BOOL = """
    void set_CATTRNAME(const std::vector<std::vector<std::vector<bool>>>& v) {
        int n3 = static_cast<int>(v.size());
        int n2 = n3 > 0 ? static_cast<int>(v[0].size()) : 0;
        int n1 = n2 > 0 ? static_cast<int>(v[0][0].size()) : 0;
        int shape[] = { n1, n2, n3 }; 
        
        std::vector<int> flat; 
        flat.reserve(n1*n2*n3);
        for(int k=0; k<n3; ++k) {
          for(int j=0; j<n2; ++j) {
            for(int i=0; i<n1; ++i) {
               flat.push_back(v[k][j][i] ? 1 : 0);
            }
          }
        }
        STRUCTNAME_set_FATTRNAME(fortran_ptr_, flat.data(), shape);
    }
"""

std_types = ["real", "real16", "integer", "complex", "integer8", "logical"]

for tname in std_types:
    info = STANDARD_TYPES[tname]

    if tname == "logical":
        acc_1d = CPP_ARRAY_1D_SET_ACCESSOR_BOOL
        acc_2d = CPP_ARRAY_2D_SET_ACCESSOR_BOOL
        acc_3d = CPP_ARRAY_3D_SET_ACCESSOR_BOOL
    else:
        acc_1d = CPP_ARRAY_1D_SET_ACCESSOR_GENERIC
        acc_2d = CPP_ARRAY_2D_SET_ACCESSOR_GENERIC
        acc_3d = CPP_ARRAY_3D_SET_ACCESSOR_GENERIC

    for ptr_type in ["NOT", "ALLOC", "PTR"]:
        # 1D
        if ptr_type == "ALLOC":
            templates[FullType(tname, 1, ptr_type)] = make_array_template(
                1,
                ptr_type,
                info.c_type,
                CPP_ALLOC_1D_DECL,
                CPP_ARRAY_1D_ALLOC_ACCESSOR,
                CPP_ARRAY_SET_DECL,
                acc_1d,
                fortran_gen_func=generate_fortran_alloc_routine,
            )
        else:
            templates[FullType(tname, 1, ptr_type)] = make_array_template(
                1, ptr_type, info.c_type, CPP_ARRAY_1D_DECL, CPP_ARRAY_1D_ACCESSOR, CPP_ARRAY_SET_DECL, acc_1d
            )
        # 2D
        templates[FullType(tname, 2, ptr_type)] = make_array_template(
            2, ptr_type, info.c_type, CPP_ARRAY_2D_DECL, CPP_ARRAY_2D_ACCESSOR, CPP_ARRAY_SET_DECL, acc_2d
        )
        # 3D
        templates[FullType(tname, 3, ptr_type)] = make_array_template(
            3, ptr_type, info.c_type, CPP_ARRAY_3D_DECL, CPP_ARRAY_3D_ACCESSOR, CPP_ARRAY_SET_DECL, acc_3d
        )

# ---------------------------------------------------------------------------
# Derived Types
# ---------------------------------------------------------------------------

# Type Scalar (NOT)
templates[FullType("type", 0, "NOT")] = TemplateEntry(
    fortran_getter="""
  subroutine STRUCTNAME_get_FATTRNAME(struct_obj_ptr, ptr_out) bind(c, name='STRUCTNAME_get_FATTRNAME')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    type(c_ptr), intent(out) :: ptr_out
    type(STRUCTNAME), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    ptr_out = c_loc(struct_obj%FATTRNAME)
  end subroutine
""",
    fortran_setter=FORTRAN_TYPE_SCALAR_SETTER,
    cpp_get_decl=CPP_TYPE_SCALAR_DECL,
    cpp_get_accessors=[
        """
    ${return_proxy_name} CATTRNAME() const {
        void* ptr;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &ptr);
        return ${return_proxy_name}(ptr);
    }
"""
    ],
    cpp_set_decl=CPP_TYPE_SCALAR_SET_DECL,
    cpp_set_accessors=[CPP_TYPE_SCALAR_SET_ACCESSOR],
)

# Type Scalar (PTR) - Returns Optional
templates[FullType("type", 0, "PTR")] = TemplateEntry(
    fortran_getter=FORTRAN_POINTER_GETTER,
    fortran_setter=FORTRAN_TYPE_POINTER_SETTER,
    cpp_get_decl=CPP_TYPE_SCALAR_DECL,
    cpp_get_accessors=[
        """
    std::optional<${return_proxy_name}> CATTRNAME() const {
        void* ptr;
        STRUCTNAME_get_FATTRNAME(fortran_ptr_, &ptr);
        if (!ptr) return std::nullopt;
        return ${return_proxy_name}(ptr);
    }
"""
    ],
    cpp_set_decl=CPP_TYPE_POINTER_SET_DECL,
    cpp_set_accessors=[CPP_TYPE_POINTER_SET_ACCESSOR],
)

# Derived Type C++ definitions
CPP_TYPE_ARRAY_1D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(
        const void* s,
        void** d,
        int* bounds,
        bool* is_alloc,
        size_t* el_size
    );
"""

# Alloc 1D - Derived Types
CPP_TYPE_ALLOC_1D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(const void* s, void** d, int* bounds, bool* is_alloc, size_t* el_size);
    void STRUCTNAME_reallocate_FATTRNAME(void* s, int lb, size_t n);
"""

CPP_TYPE_ARRAY_1D_ACCESSOR = """
    ${return_proxy_name}Array1D CATTRNAME() const {
        return ProxyHelpers::get_type_array_1d<${return_proxy_name}Array1D>(
            fortran_ptr_,
            STRUCTNAME_get_FATTRNAME_info
        );
    }
"""

# Derived Type Alloc Accessor
CPP_TYPE_ARRAY_1D_ALLOC = """
    ${return_proxy_name}Alloc1D CATTRNAME() const {
        return ${return_proxy_name}Alloc1D(
            const_cast<void*>(fortran_ptr_),
            STRUCTNAME_reallocate_FATTRNAME,
            STRUCTNAME_get_FATTRNAME_info
        );
    }
"""

CPP_TYPE_ARRAY_2D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(
        const void* s, void** d, int* bounds, int* strides, bool* a, size_t* es
    );
"""

CPP_TYPE_ARRAY_2D_ACCESSOR = """
    ${return_proxy_name}Array2D CATTRNAME() const {
        return ProxyHelpers::get_type_array_2d<${return_proxy_name}Array2D>(
            fortran_ptr_,
            STRUCTNAME_get_FATTRNAME_info
        );
    }
"""
CPP_TYPE_ARRAY_3D_DECL = """
    void STRUCTNAME_get_FATTRNAME_info(
         const void* s, void** d, int* bounds, int* strides, bool* a, size_t* es
    );
"""
CPP_TYPE_ARRAY_3D_ACCESSOR = """
    ${return_proxy_name}Array3D CATTRNAME() const {
        return ProxyHelpers::get_type_array_3d<${return_proxy_name}Array3D>(
            fortran_ptr_,
            STRUCTNAME_get_FATTRNAME_info
        );
    }
"""

# Register Type Arrays
for ptr_type in ["NOT", "ALLOC", "PTR"]:
    # 1D Derived
    if ptr_type == "ALLOC":
        templates[FullType("type", 1, ptr_type)] = make_array_template(
            1,
            ptr_type,
            "void*",
            CPP_TYPE_ALLOC_1D_DECL,
            CPP_TYPE_ARRAY_1D_ALLOC,
            is_derived_type=True,
            fortran_gen_func=generate_fortran_alloc_routine,
        )
    else:
        templates[FullType("type", 1, ptr_type)] = make_array_template(
            1, ptr_type, "void*", CPP_TYPE_ARRAY_1D_DECL, CPP_TYPE_ARRAY_1D_ACCESSOR, is_derived_type=True
        )

    # 2D Derived
    templates[FullType("type", 2, ptr_type)] = make_array_template(
        2, ptr_type, "void*", CPP_TYPE_ARRAY_2D_DECL, CPP_TYPE_ARRAY_2D_ACCESSOR, is_derived_type=True
    )
    # 3D Derived
    templates[FullType("type", 3, ptr_type)] = make_array_template(
        3, ptr_type, "void*", CPP_TYPE_ARRAY_3D_DECL, CPP_TYPE_ARRAY_3D_ACCESSOR, is_derived_type=True
    )

# Character 1D (Standard fixed size array via generic char handler)
templates[FullType("character", 1, "NOT")] = make_char_array_1d("NOT")
templates[FullType("character", 1, "ALLOC")] = make_char_array_1d("ALLOC")


# ---------------------------------------------------------------------------
# Dispatch-based scalar access
# ---------------------------------------------------------------------------
# Types eligible for dispatch: simple scalars with ptr=NOT, dim=0
DISPATCHABLE_SCALAR_TYPES: set[ArgumentType] = {"real", "real16", "integer", "integer8", "logical", "complex"}


def is_dispatchable_scalar(arg: Argument) -> bool:
    ft = arg.full_type
    return ft.dim == 0 and ft.ptr == "NOT" and ft.type in DISPATCHABLE_SCALAR_TYPES and arg.is_component


def group_dispatch_fields(
    struct: CodegenStructure,
) -> dict[ArgumentType, list[tuple[int, Argument]]]:
    """Group dispatchable scalar fields by type, assigning sequential field IDs."""
    groups: dict[ArgumentType, list[tuple[int, Argument]]] = {}
    for arg in struct.arg:
        if not is_dispatchable_scalar(arg):
            continue
        type_name = arg.full_type.type
        if type_name not in groups:
            groups[type_name] = []
        groups[type_name].append((len(groups[type_name]), arg))
    return groups


def generate_fortran_dispatch_getter(
    struct_name: str, type_name: ArgumentType, fields: list[tuple[int, Argument]]
) -> str:
    info = STANDARD_TYPES[type_name]
    cases = "\n".join(f"    case({fid}); value_out = struct_obj%{arg.f_name}" for fid, arg in fields)
    return f"""
  subroutine {struct_name}_get_{type_name}(struct_obj_ptr, field_id, value_out) &
      bind(c, name='{struct_name}_get_{type_name}')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    integer(c_int), intent(in), value :: field_id
    {info.fortran_type}, intent(out) :: value_out
    type({struct_name}), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    select case(field_id)
{cases}
    end select
  end subroutine
"""


def generate_fortran_dispatch_setter(
    struct_name: str, type_name: ArgumentType, fields: list[tuple[int, Argument]]
) -> str:
    info = STANDARD_TYPES[type_name]
    cases = "\n".join(f"    case({fid}); struct_obj%{arg.f_name} = value_in" for fid, arg in fields)
    return f"""
  subroutine {struct_name}_set_{type_name}(struct_obj_ptr, field_id, value_in) &
      bind(c, name='{struct_name}_set_{type_name}')
    type(c_ptr), intent(in), value :: struct_obj_ptr
    integer(c_int), intent(in), value :: field_id
    {info.fortran_type}, intent(in), value :: value_in
    type({struct_name}), pointer :: struct_obj
    call c_f_pointer(struct_obj_ptr, struct_obj)
    select case(field_id)
{cases}
    end select
  end subroutine
"""


# ---------------------------------------------------------------------------
# Remaining Code Generation Functions
# ---------------------------------------------------------------------------
def split_signature(cpp_template: str, class_name: str) -> tuple[str, str]:
    clean_template = cpp_template.strip()
    assert "{" in clean_template
    signature = clean_template[: clean_template.find("{")].strip()
    header_declaration = signature + ";"

    if "(" not in signature:
        raise ValueError(
            f"Signature '{signature}' missing '(' in class '{class_name}' from template:\n{cpp_template}"
        )

    ret_type_and_method, args = signature.split("(", 1)
    ret_type, method = ret_type_and_method.rsplit(" ", 1)
    impl_sig = f"{ret_type} {class_name}::{method}({args}"
    implementation = clean_template.replace(signature, impl_sig)
    return header_declaration, implementation


def generate_accessor_code(
    struct_name: str,
    arg: Argument,
    full_type: FullType,
    attr_kind: str = "",
):
    try:
        tpl = templates[full_type]
    except KeyError as ex:
        raise ValueError(f"Unsupported type: {full_type}") from ex

    to_replace = {
        "structname": struct_name.lower(),
        "fattrname": arg.f_name,
        "cattrname": arg.c_name,
        "fortrantype": STANDARD_TYPES[full_type.type].fortran_type,
    }

    if full_type.type in _native_alloc_map:
        to_replace["alloc_type"] = _native_alloc_map[full_type.type]

    if attr_kind:
        to_replace["attrtype"] = attr_kind

    def replace_all(s: str) -> str:
        return subst(s, **to_replace)

    return {
        "fortran_getter": replace_all(tpl.fortran_getter),
        "fortran_setter": replace_all(tpl.fortran_setter) if tpl.fortran_setter else None,
        "cpp_get_decl": replace_all(tpl.cpp_get_decl),
        "cpp_get_accessors": [replace_all(acc) for acc in tpl.cpp_get_accessors],
        "cpp_set_decl": replace_all(tpl.cpp_set_decl) if tpl.cpp_set_decl else None,
        "cpp_set_accessors": [replace_all(acc) for acc in tpl.cpp_set_accessors]
        if tpl.cpp_set_accessors
        else [],
    }


def create_fortran_proxy_code(structs: list[CodegenStructure]) -> str:
    container_types = []
    # TODO: only generate containers if they're used
    for struct in structs:
        container_types.append(
            f"""\

  type :: {struct.container_alloc_name}
    type({struct.f_name}), allocatable :: data(:)
  end type {struct.container_alloc_name}
            """.rstrip()
        )

    for typ in native_type_containers:
        container_types.append(
            f"""\

  type :: {typ.name}_container_alloc
    {typ.fortran_type}, allocatable :: data(:)
  end type {typ.name}_container_alloc
            """.rstrip()
        )

    container_types.append(
        """\

  type :: character_container_alloc
    character(:), allocatable :: data(:)
  end type character_container_alloc
        """.rstrip()
    )

    imports_map = {}
    for struct in structs:
        module = struct.module
        imports_map.setdefault(module, set())
        imports_map[module].add(struct.f_name)

    import_lines = []
    for module, imps in imports_map.items():
        import_lines.append(f"  use {module}, only: " + ", ".join(sorted(imps)))

    imports_str = "\n".join(import_lines)
    containers_str = "\n".join(container_types)

    output = []
    output.append(
        f"""\
module bmad_struct_proxy_mod
  use bmad_struct
  use tao_struct
  use test_struct_defs
  use, intrinsic :: iso_c_binding

  {imports_str}
  {containers_str}

contains
"""
    )

    for nt in native_type_containers:
        struct_name = nt.name.lower()
        container_name = nt.fortran_container_struct
        output.append(
            f"""
  function allocate_{struct_name}_container() result(ptr) bind(c)
    implicit none
    type(c_ptr) :: ptr
    type({container_name}), pointer :: ctr
    allocate(ctr)
    ptr = c_loc(ctr)
  end function

  subroutine deallocate_{struct_name}_container(ptr) bind(c)
    implicit none
    type(c_ptr), value :: ptr
    type({container_name}), pointer :: ctr
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, ctr)
      deallocate(ctr)
    end if
  end subroutine

  subroutine reallocate_{struct_name}_container_data(container_ptr, lbound_, n) bind(c)
    implicit none
    type(c_ptr), value :: container_ptr
    integer(c_int), value :: lbound_
    integer(c_size_t), value :: n
    type({container_name}), pointer :: ctr

    if (.not. c_associated(container_ptr)) return
    call c_f_pointer(container_ptr, ctr)

    if (n == 0) then
      if (allocated(ctr%data)) deallocate(ctr%data)
    else
      if (allocated(ctr%data)) deallocate(ctr%data)
      allocate(ctr%data(lbound_:lbound_ + n - 1))
    end if
  end subroutine

  subroutine access_{struct_name}_container(container_ptr, d_ptr, bounds, is_allocated) bind(c)
    use iso_c_binding
    implicit none
    type(c_ptr), value :: container_ptr
    type(c_ptr), intent(out) :: d_ptr
    integer(c_int), dimension(2), intent(out) :: bounds
    logical(c_bool), intent(out) :: is_allocated

    type({container_name}), pointer :: ctr

    if (.not. c_associated(container_ptr)) then
       is_allocated = .false.
       return
    endif

    call c_f_pointer(container_ptr, ctr)

    if (allocated(ctr%data)) then
      is_allocated = .true.
      bounds(1) = int(lbound(ctr%data, 1), c_int)
      bounds(2) = int(ubound(ctr%data, 1), c_int)
      d_ptr = c_loc(ctr%data(bounds(1)))
    else
      is_allocated = .false.
      d_ptr = c_null_ptr
      bounds = 0
    endif
  end subroutine
            """
        )

    # Character container (special: deferred-length, needs str_len in realloc/access)
    output.append(
        """
  function allocate_character_container() result(ptr) bind(c)
    implicit none
    type(c_ptr) :: ptr
    type(character_container_alloc), pointer :: ctr
    allocate(ctr)
    ptr = c_loc(ctr)
  end function

  subroutine deallocate_character_container(ptr) bind(c)
    implicit none
    type(c_ptr), value :: ptr
    type(character_container_alloc), pointer :: ctr
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, ctr)
      deallocate(ctr)
    end if
  end subroutine

  subroutine reallocate_character_container_data(container_ptr, lbound_, n, str_len) bind(c)
    implicit none
    type(c_ptr), value :: container_ptr
    integer(c_int), value :: lbound_
    integer(c_size_t), value :: n
    integer(c_int), value :: str_len
    type(character_container_alloc), pointer :: ctr

    if (.not. c_associated(container_ptr)) return
    call c_f_pointer(container_ptr, ctr)

    if (n == 0) then
      if (allocated(ctr%data)) deallocate(ctr%data)
    else
      if (allocated(ctr%data)) deallocate(ctr%data)
      allocate(character(str_len) :: ctr%data(lbound_:lbound_ + n - 1))
    end if
  end subroutine

  subroutine access_character_container(container_ptr, d_ptr, bounds, str_len, is_allocated) bind(c)
    use iso_c_binding
    implicit none
    type(c_ptr), value :: container_ptr
    type(c_ptr), intent(out) :: d_ptr
    integer(c_int), dimension(2), intent(out) :: bounds
    integer(c_int), intent(out) :: str_len
    logical(c_bool), intent(out) :: is_allocated

    type(character_container_alloc), pointer :: ctr

    if (.not. c_associated(container_ptr)) then
       is_allocated = .false.
       return
    endif

    call c_f_pointer(container_ptr, ctr)

    if (allocated(ctr%data)) then
      is_allocated = .true.
      bounds(1) = int(lbound(ctr%data, 1), c_int)
      bounds(2) = int(ubound(ctr%data, 1), c_int)
      str_len = int(len(ctr%data), c_int)
      d_ptr = c_loc(ctr%data(bounds(1)))
    else
      is_allocated = .false.
      d_ptr = c_null_ptr
      bounds = 0
      str_len = 0
    endif
  end subroutine
        """
    )

    for struct in structs:
        struct_name = struct.f_name.lower()
        output.append(f"  !! {struct_name}")
        output.append(
            f"""
    function allocate_fortran_{struct_name}(n, element_size) result(ptr) bind(c)
      implicit none
      integer(c_int), value :: n
      integer(c_size_t), intent(out) :: element_size
      type(c_ptr) :: ptr
      type({struct_name}), pointer :: fptr
      type({struct_name}), pointer :: fptr_array(:)

      if (n <= 0) then
        allocate(fptr)
        ptr = c_loc(fptr)
        element_size = int(storage_size(fptr) / 8, c_size_t)
      else
        allocate(fptr_array(n))
        ptr = c_loc(fptr_array)
        element_size = int(storage_size(fptr_array(1)) / 8, c_size_t)
      end if
    end function

    subroutine deallocate_fortran_{struct_name}(ptr, n) bind(c)
      implicit none
      type(c_ptr), value :: ptr
      integer(c_int), value :: n
      type({struct_name}), pointer :: fptr
      type({struct_name}), pointer :: fptr_array(:)

      if (c_associated(ptr)) then
        if (n <= 0) then
          call c_f_pointer(ptr, fptr)
          deallocate(fptr)
        else
          call c_f_pointer(ptr, fptr_array, [n])
          deallocate(fptr_array)
        end if
      end if
    end subroutine

  subroutine copy_fortran_{struct_name}(src_ptr, dst_ptr) bind(c)
    implicit none
    type(c_ptr), value :: src_ptr, dst_ptr
    type({struct_name}), pointer :: src, dst

    if (c_associated(src_ptr) .and. c_associated(dst_ptr)) then
      call c_f_pointer(src_ptr, src)
      call c_f_pointer(dst_ptr, dst)
      dst = src  ! Fortran derived type assignment
    end if
  end subroutine

  function allocate_{struct_name}_container() result(ptr) bind(c)
    implicit none
    type(c_ptr) :: ptr
    type({struct.container_alloc_name}), pointer :: ctr
    allocate(ctr)
    ptr = c_loc(ctr)
  end function

  subroutine deallocate_{struct_name}_container(ptr) bind(c)
    implicit none
    type(c_ptr), value :: ptr
    type({struct.container_alloc_name}), pointer :: ctr
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, ctr)
      deallocate(ctr)
    end if
  end subroutine

  subroutine reallocate_{struct_name}_container_data(container_ptr, lbound_, n) bind(c)
    implicit none
    type(c_ptr), value :: container_ptr
    integer(c_int), value :: lbound_
    integer(c_size_t), value :: n
    type({struct.container_alloc_name}), pointer :: ctr

    if (.not. c_associated(container_ptr)) return
    call c_f_pointer(container_ptr, ctr)

    if (n == 0) then
      if (allocated(ctr%data)) deallocate(ctr%data)
    else
      if (allocated(ctr%data)) deallocate(ctr%data)
      allocate(ctr%data(lbound_:lbound_ + n - 1))
    end if
  end subroutine

  subroutine access_{struct_name}_container(container_ptr, d_ptr, bounds, is_allocated, elem_size) bind(c)
    use iso_c_binding
    implicit none
    type(c_ptr), value :: container_ptr
    type(c_ptr), intent(out) :: d_ptr
    integer(c_int), dimension(2), intent(out) :: bounds
    logical(c_bool), intent(out) :: is_allocated
    integer(c_size_t), intent(out) :: elem_size

    type({struct.container_alloc_name}), pointer :: ctr

    if (.not. c_associated(container_ptr)) then
       is_allocated = .false.
       return
    endif

    call c_f_pointer(container_ptr, ctr)

    if (allocated(ctr%data)) then
      is_allocated = .true.
      bounds(1) = int(lbound(ctr%data, 1), c_int)
      bounds(2) = int(ubound(ctr%data, 1), c_int)
      ! Use intrinsic storage_size (returns bits) divided by 8 for bytes
      elem_size = storage_size(ctr%data(bounds(1))) / 8
      d_ptr = c_loc(ctr%data(bounds(1)))
    else
      is_allocated = .false.
      d_ptr = c_null_ptr
      bounds = 0
      elem_size = 0
    endif
  end subroutine
    """
        )

        # Dispatch-based scalar accessors
        dispatch_groups = group_dispatch_fields(struct)
        dispatched_fields: set[str] = set()
        for type_name, fields in sorted(dispatch_groups.items()):
            output.append(f"  ! dispatch: {struct_name}%{type_name} ({len(fields)} fields)")
            output.append(generate_fortran_dispatch_getter(struct_name, type_name, fields))
            output.append(generate_fortran_dispatch_setter(struct_name, type_name, fields))
            for _, arg in fields:
                dispatched_fields.add(arg.f_name)

        for arg in struct.arg:
            if not arg.is_component:
                continue
            if arg.f_name in dispatched_fields:
                continue
            try:
                acc = generate_accessor_code(struct.f_name, arg, arg.full_type, arg.kind)
            except ValueError as ex:
                output.append(f"  ! skipped {struct.f_name}%{arg.f_name}: {ex}")
                continue

            output.append(f"  ! {struct.f_name}%{arg.f_name}: {arg.full_type}")
            output.append(acc["fortran_getter"])
            if acc["fortran_setter"]:
                output.append(acc["fortran_setter"])

    output.append("end module")
    return "\n".join(output)


def infer_cpp_type(arg) -> str | None:
    ft = arg.full_type
    if ft.dim > 0 and (ft.type in {"character", "type"}):
        return None

    if ft.type == "type":
        base = struct_to_proxy_class_name(arg.kind)
    elif ft.type == "character":
        base = "std::string"
    elif ft.type == "complex":
        base = "std::complex<double>"
    elif ft.type == "logical":
        base = "bool"
    else:
        if ft.type not in STANDARD_TYPES:
            return None
        base = STANDARD_TYPES[ft.type].c_type

    res = base
    if ft.dim > 0:
        for _ in range(ft.dim):
            res = f"std::vector<{res}>"

    return res


def _generate_proxy_constructor_arg(
    struct: CodegenStructure, arg: Argument
) -> tuple[str, str] | tuple[None, None]:
    if not arg.is_component:
        return None, None

    try:
        acc = generate_accessor_code(struct.f_name, arg, arg.full_type, arg.kind)
    except ValueError:
        return None, None

    if not acc["cpp_set_decl"]:
        return None, None

    setter_name = f"set_{arg.c_name}"
    cpp_type = infer_cpp_type(arg)

    if not cpp_type:
        return None, None

    if arg.full_type.type in {"type"}:
        ctor_arg = f"optional_ref<const {cpp_type}>"
        ctor_body = f"    if ({arg.c_name}) {setter_name}({arg.c_name}->get());"
    else:
        ctor_arg = f"std::optional<{cpp_type}>"
        ctor_body = f"    if ({arg.c_name}) {setter_name}(*{arg.c_name});"

    return (ctor_arg, ctor_body)


def _generate_proxy_constructor_args(struct: CodegenStructure) -> dict[str, tuple[str, str]]:
    res: dict[str, tuple[str, str]] = {}
    for arg in struct.arg:
        ctor_arg, ctor_body = _generate_proxy_constructor_arg(struct, arg)
        if ctor_arg is not None and ctor_body is not None:
            res[arg.f_name] = (ctor_arg, ctor_body)
    return res


def _generate_proxy_constructor(
    struct: CodegenStructure,
    proxy_class_name: str,
) -> str | None:
    ctor_args = _generate_proxy_constructor_args(struct)
    if not ctor_args:
        return None
    args = [f"{ctor_type} {name} = std::nullopt" for name, (ctor_type, _) in ctor_args.items()]
    init_body = [ctor_init for _, ctor_init in ctor_args.values()]
    ctor_args_str = ",\n        ".join(args)
    ctor_inits_str = "\n".join(init_body)
    return f"""
    explicit {proxy_class_name}(
        {ctor_args_str}
    ) : FortranProxy() {{
{ctor_inits_str}
    }}
"""


def get_proxy_header_and_code(
    header_template_src: str,
    cpp_template_src: str,
    structs: list[CodegenStructure],
) -> tuple[str, str]:
    c_forward_declarations: list[str] = []
    subs: dict[str, str] = {}

    class_template = Template(
        """

template <>
struct FortranTraits<${class_name}> {
  static void* allocate() {
    size_t sz;
    return allocate_fortran_${struct_name}(0, &sz);
  }
  static void deallocate(void* ptr) noexcept {
    deallocate_fortran_${struct_name}(ptr, 0);
  }
  static void copy(const void* src, void* dst) {
    copy_fortran_${struct_name}(src, dst);
  }
  static constexpr std::string_view type_name() {
    return "${struct_name}";
  }
};

class ${class_name} : public FortranProxy<${class_name}> {
 public:
  using FortranProxy::FortranProxy;
  using FortranProxy::operator=;

  ${class_body}
};

"""
    )

    classes = {}
    all_impl = []

    for struct in structs:
        proxy_class_name = struct_to_proxy_class_name(struct.f_name)
        class_body: list[str] = []

        ctor_code = _generate_proxy_constructor(struct, proxy_class_name)
        if ctor_code:
            class_body.append(ctor_code)

        # Dispatch-based scalar accessors
        dispatch_groups = group_dispatch_fields(struct)
        dispatched_fields: dict[str, tuple[ArgumentType, int]] = {}
        dispatch_struct_name = struct.f_name.lower()
        for type_name, fields in sorted(dispatch_groups.items()):
            info = STANDARD_TYPES[type_name]
            c_forward_declarations.append(
                f"    void {dispatch_struct_name}_get_{type_name}"
                f"(const void* struct_obj, int field_id, {info.c_type}* value_out);"
            )
            c_forward_declarations.append(
                f"    void {dispatch_struct_name}_set_{type_name}"
                f"(void* struct_obj, int field_id, {info.c_type} value_in);"
            )
            for field_id, arg in fields:
                dispatched_fields[arg.f_name] = (type_name, field_id)

        for arg in struct.arg:
            if not arg.is_component:
                continue

            if arg.f_name in dispatched_fields:
                type_name, field_id = dispatched_fields[arg.f_name]
                info = STANDARD_TYPES[type_name]
                c_type = info.c_type
                c_name = arg.c_name
                struct_name = dispatch_struct_name

                getter_body = (
                    f"    {c_type} {c_name}() const {{\n"
                    f"        {c_type} value;\n"
                    f"        {struct_name}_get_{type_name}(fortran_ptr_, {field_id}, &value);\n"
                    f"        return value;\n"
                    f"    }}"
                )
                sig, impl = split_signature(getter_body, proxy_class_name)
                all_impl.append(impl)
                class_body.append(f"{sig} // {arg.full_type} [dispatch:{field_id}]")

                setter_body = (
                    f"    void set_{c_name}({c_type} value) {{\n"
                    f"        {struct_name}_set_{type_name}(fortran_ptr_, {field_id}, value);\n"
                    f"    }}"
                )
                sig, impl = split_signature(setter_body, proxy_class_name)
                all_impl.append(impl)
                class_body.append(sig)
                continue

            try:
                acc = generate_accessor_code(struct.f_name, arg, arg.full_type, arg.kind)
            except ValueError as ex:
                logger.warning(f"Proxy class {struct.f_name}%{arg.f_name} skipped: {ex}")
                continue

            c_forward_declarations.append(acc["cpp_get_decl"])

            for accessor_body in acc["cpp_get_accessors"]:
                if arg.full_type.type == "type":
                    accessor_body = Template(accessor_body).substitute(
                        return_proxy_name=struct_to_proxy_class_name(arg.kind)
                    )
                sig, impl = split_signature(accessor_body, proxy_class_name)
                all_impl.append(impl)
                class_body.append(f"{sig} // {arg.full_type}")

            if acc["cpp_set_decl"]:
                c_forward_declarations.append(acc["cpp_set_decl"])
                for accessor_body in acc["cpp_set_accessors"]:
                    if arg.full_type.type == "type":
                        accessor_body = Template(accessor_body).substitute(
                            return_proxy_name=struct_to_proxy_class_name(arg.kind)
                        )
                    sig, impl = split_signature(accessor_body, proxy_class_name)
                    all_impl.append(impl)
                    class_body.append(sig)

        subs[f"{struct.f_name}_class_body"] = "\n".join(class_body)
        classes[struct.f_name.lower()] = class_body

    class_forward_declarations = []
    proxy_classes = []

    class_forward_declarations.append('extern "C" {')

    for nt in native_type_containers:
        name = nt.name
        class_forward_declarations.append(f"""
  void* allocate_{name}_container();
  void reallocate_{name}_container_data(void *, int, size_t) noexcept;
  void deallocate_{name}_container(void *) noexcept;
  void access_{name}_container(const void* handle, void** data, int* bounds, bool* alloc);
""")

    # Character container (special: deferred-length character arrays)
    class_forward_declarations.append("""
  void* allocate_character_container();
  void reallocate_character_container_data(void *, int, size_t, int) noexcept;
  void deallocate_character_container(void *) noexcept;
  void access_character_container(const void* handle, void** data, int* bounds, int* str_len, bool* alloc);
""")

    for struct_name in classes:
        class_forward_declarations.append(f"""
  void* allocate_fortran_{struct_name}(int n, size_t *element_size);
  void deallocate_fortran_{struct_name}(void* ptr, int n) noexcept;
  void copy_fortran_{struct_name}(const void* src, void* dst);

  void* allocate_{struct_name}_container();
  void reallocate_{struct_name}_container_data(void *, int, size_t) noexcept;
  void deallocate_{struct_name}_container(void *) noexcept;
  void access_{struct_name}_container(const void* handle, void** data, int* bounds, bool* alloc, size_t* elem_size);
  """)

    class_forward_declarations.append("}")

    class_forward_declarations.append("""
template <typename T>
using optional_ref = std::optional<std::reference_wrapper<T>>;
""")

    for nt in native_type_containers:
        name = nt.name
        class_forward_declarations.append(f"""
struct {nt.cpp_container_name} : public FAlloc1D<{nt.cpp_type}> {{
    using Base = FAlloc1D<{nt.cpp_type}>;
    using Base::Base;
    {nt.cpp_container_name}() : Base(
        allocate_{name}_container,
        deallocate_{name}_container,
        reallocate_{name}_container_data,
        access_{name}_container
    ) {{}}
    {nt.cpp_container_name}(int n) : Base(
        n,
        allocate_{name}_container,
        deallocate_{name}_container,
        reallocate_{name}_container_data,
        access_{name}_container
    ) {{}}
    {nt.cpp_container_name}(void* handle, ReallocFuncPtr realloc, PrimAccessFuncPtr access)
       : Base(handle, realloc, access) {{}}
}};
""")

    # CharacterAlloc1D container
    class_forward_declarations.append("""
struct CharacterAlloc1D : public FCharAlloc1D {
    using Base = FCharAlloc1D;
    using Base::Base;
    CharacterAlloc1D() : Base(
        allocate_character_container,
        deallocate_character_container,
        reallocate_character_container_data,
        access_character_container
    ) {}
    CharacterAlloc1D(int n, int str_len = 200) : Base(
        n, str_len,
        allocate_character_container,
        deallocate_character_container,
        reallocate_character_container_data,
        access_character_container
    ) {}
    CharacterAlloc1D(void* handle, CharReallocFuncPtr realloc, CharAccessFuncPtr access)
       : Base(handle, realloc, access) {}
};
""")

    for struct_name, class_body in classes.items():
        struct_name = struct_name.lower()
        class_name = struct_to_proxy_class_name(struct_name)
        class_forward_declarations.append(f"class {class_name};")
        class_forward_declarations.append(f"""
using {class_name}Array1D = FTypeArray1D<
    {class_name},
    allocate_fortran_{struct_name},
    deallocate_fortran_{struct_name}
>;
using {class_name}Array2D = FTypeArray2D<{class_name}>;
using {class_name}Array3D = FTypeArray3D<{class_name}>;

struct {class_name}Alloc1D : public FTypeAlloc1D<{class_name}Array1D> {{
    using Base = FTypeAlloc1D<{class_name}Array1D>;
    using Base::Base;
    {class_name}Alloc1D() : Base(
        allocate_{struct_name}_container,
        deallocate_{struct_name}_container,
        reallocate_{struct_name}_container_data,
        access_{struct_name}_container
    ) {{}}
    {class_name}Alloc1D(int n) : Base(
        n,
        allocate_{struct_name}_container,
        deallocate_{struct_name}_container,
        reallocate_{struct_name}_container_data,
        access_{struct_name}_container
    ) {{}}
    {class_name}Alloc1D(void* handle, ReallocFuncPtr realloc, TypeAccessFuncPtr access) 
       : Base(handle, realloc, access) {{}}
}};
        """)
        proxy_classes.append(
            class_template.substitute(
                struct_name=struct_name,
                class_name=struct_to_proxy_class_name(struct_name),
                class_body="\n".join(class_body),
            )
        )

    subs["c_forward_declarations"] = "\n".join(c_forward_declarations)
    subs["class_forward_declarations"] = "\n".join(class_forward_declarations)
    subs["proxy_classes"] = "\n".join(proxy_classes)

    header = Template(header_template_src.replace("// ${", "${")).substitute(subs)
    impl = cpp_template_src + "\n".join(all_impl)
    return header, impl


def get_proxy_classes(structs: list[CodegenStructure]) -> dict[pathlib.Path, str]:
    file_to_contents = {}
    generated = CPPBMAD_SRC / "generated"

    file_to_contents[generated / "proxy_mod.f90"] = create_fortran_proxy_code(structs)

    cpp_proxy_header_template = (CODEGEN_ROOT / "proxy.tpl.hpp").read_text()
    cpp_proxy_cpp_template = (CODEGEN_ROOT / "proxy.tpl.cpp").read_text()

    proxy_header, proxy_impl = get_proxy_header_and_code(
        cpp_proxy_header_template,
        cpp_proxy_cpp_template,
        structs,
    )
    file_to_contents[CPPBMAD_ROOT / "include" / "bmad" / "generated" / "proxy.hpp"] = proxy_header
    file_to_contents[generated / "proxy.cpp"] = proxy_impl
    return file_to_contents
