from __future__ import annotations

import os
import pathlib
import string
import subprocess
import textwrap

from .arg import Argument, CodegenStructure
from .context import SUPPORTED_ARRAY_DIMS, CodegenConfig
from .cpp import CppWrapperArgument
from .enums import EnumValue, get_ele_attributes, get_ele_keys
from .paths import CODEGEN_ROOT, PYBMAD_INCLUDE, PYBMAD_LIB, PYBMAD_ROOT, PYBMAD_SRC
from .proxy import _generate_proxy_constructor_arg
from .proxy import templates as proxy_templates
from .routines import FortranRoutine, RoutineArg
from .types import remove_optional
from .util import snake_to_camel, sorted_routines


def struct_array_usage_dimensions(
    routines_by_name: dict[str, FortranRoutine],
    structs: list[CodegenStructure],
) -> dict[str, set[int]]:
    by_kind = {}

    args = []
    for routine in routines_by_name.values():
        if routine.usable:
            for arg in routine.args:
                args.append(arg)

    for st in structs:
        for arg in st.arg:
            args.append(arg)

    for arg in args:
        if arg.full_type.type == "type":
            kind = arg.kind.lower()
            by_kind.setdefault(kind, set())
            by_kind[kind].add(arg.full_type.dim)

    return by_kind


def generate_enum_wrapper(clsname: str, enums: list[EnumValue], *, offset: int | None = 0) -> str:
    """
    Generate a Python enum.IntEnum class definition string.
    """
    code = ["", "", f"class {clsname}(enum.IntEnum):"]

    if not enums:
        code.append("    pass")
        return "\n".join(code)

    if not offset:
        offset_suffix = ""
    elif offset < 0:
        offset_suffix = f" - {abs(offset)}"
    else:
        offset_suffix = f" + {offset}"

    for attr in enums:
        comment = f"  # {attr.comment}" if attr.comment else ""
        code.append(f"    {attr.name} = {attr.value}{offset_suffix}{comment}")

    if offset:
        f_offset = -offset
        code.append(f"""
    @property
    def fortran_value(self) -> int:
        return self.value + {f_offset}
    """)
    return "\n".join(code)


def generate_enum_constants(
    all_enums: dict[str, dict[str, EnumValue]],
    custom_enums: dict[str, str],
) -> str:
    """
    Generate Python module-level constants.
    """
    result = []
    for fn, enums in all_enums.items():
        result.append("")
        result.append(f"# Constants from {fn}")
        for enum in enums.values():
            if enum.comment:
                result.append(f"# {enum.comment}")

            if enum.name in custom_enums:
                enum_cls = custom_enums[enum.name]
                result.append(f"{enum.name} = {enum_cls}.{enum.name}")
            else:
                result.append(f"{enum.name} = {enum.value}")

    return "\n".join(result)


def generate_enum_wrapper_code(enums: dict[str, dict[str, EnumValue]]) -> str:
    ele_attrs = get_ele_attributes(enums["bmad_struct.f90"])
    ele_keys = get_ele_keys(enums["bmad_struct.f90"])
    custom_enums = {}
    for enum in ele_attrs:
        custom_enums[enum.name] = "EleAttribute"
    for enum in ele_keys:
        custom_enums[enum.name] = "EleKey"
    return "\n".join(
        (
            "import enum",
            "",
            generate_enum_wrapper("EleAttribute", ele_attrs, offset=-1),
            generate_enum_wrapper("EleKey", ele_keys),
            generate_enum_constants(enums, custom_enums),
        )
    )


def generate_routine_return_value_wrapper(routine: FortranRoutine) -> list[str]:
    immut_args = [arg for arg in routine.args if arg.is_python_immutable]
    outputs = [*routine.outputs, *immut_args]

    if len(outputs) <= 1 and not immut_args:  # TODO: immut args even for 1 -> struct
        return []

    clsname, full_clsname = routine.python_class_return_type

    lines = []
    lines.append(f'    nb::class_<{full_clsname}>(m, "{clsname}", "{routine.name} return type")')
    for arg in outputs:
        lines.append(f'        .def_ro("{arg.python_name}", &{full_clsname}::{arg.c_name})')

    lines.append(f'        .def("__len__", [](const {full_clsname} &) {{ return {len(outputs)}; }})')
    lines.append(f'        .def("__getitem__", [](const {full_clsname} &s, int i) -> nb::object {{')
    lines.append(f"            if (i < 0) i += {len(outputs)};")
    for i, arg in enumerate(outputs):
        lines.append(f"            if (i == {i}) return nb::cast(s.{arg.c_name});")
    lines.append("            throw nb::index_error();")
    lines.append("        })")

    # lines.append(f'      .def("__repr__", [](const {full_clsname} &self){{ return to_string(self); }})')
    lines.append("        ;")
    return lines


def generate_py_routine_return_value_struct(routine: FortranRoutine) -> list[str]:
    args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]
    outputs = [arg for arg in args if arg.arg.intent == "out"]

    name, py_name = routine.python_class_return_type
    orig_struct = routine.cpp_return_type

    immut_args = [arg for arg in args if arg.arg.is_python_immutable]
    arg_to_decl: list[tuple[CppWrapperArgument, tuple[str, str]]] = [
        (arg, arg.struct_decl(ignore_intent=True)) for arg in immut_args
    ]

    lines = []
    if not outputs:
        header = [f"struct {py_name} {{"]
    elif len(outputs) == 1:
        header = [f"struct {py_name} {{"]
        (output,) = outputs
        lines.append("")
        decl = output.struct_decl()
        assert decl is not None
        type_, name = decl
        header.append(f"  {type_} {name};")
    else:
        header = [f"struct {py_name} : public {orig_struct} {{"]

    lines.extend(header)

    for arg, (type_, name) in arg_to_decl:
        if arg.arg.is_optional:
            type_ = f"std::optional<{type_}>"
        lines.append(f"  {type_} {name};")

    if len(outputs) in {0, 1}:
        lines.append("};")
    else:
        lines.append(f"{py_name}({orig_struct} _base,")
        for idx, (arg, (type_, name)) in enumerate(arg_to_decl):
            if arg.arg.is_optional:
                type_ = f"std::optional<{type_}>"
            lines.append(f"  {type_} {name}")
            if idx != len(arg_to_decl) - 1:
                lines[-1] += ","
        lines.append(f") : {orig_struct}(std::move(_base)),")
        for idx, (arg, _decl) in enumerate(arg_to_decl):
            name = arg.decl_arg_name
            lines.append(f"{name}({name})")
            if idx != len(arg_to_decl) - 1:
                lines[-1] += ","
        lines.append("  {}")

        lines.append("};")
    return lines


def _get_py_routine_arg_type(arg: RoutineArg) -> str:
    cpp_type = arg.transform.cpp_type
    if arg.is_python_immutable:
        if cpp_type.startswith("optional_ref"):
            # optional_ref<T> -> std::optional<T>
            cpp_type = remove_optional(cpp_type)
            cpp_type = f"std::optional<{cpp_type}>"
        cpp_type = cpp_type.replace("&", "")
    elif cpp_type.startswith("optional_ref"):
        # nanobind cannot bind optional_ref<T> (std::optional<std::reference_wrapper<T>>),
        # so expose as T* at the Python boundary and convert back in the wrapper.
        inner = remove_optional(cpp_type)
        cpp_type = f"{inner} *"
    return cpp_type


def _get_py_routine_decl_spec(routine: FortranRoutine, allow_defaults: bool) -> list[str]:
    """py_ routine (wrapped due to immutable inout has its own definition)"""
    specs = []

    for arg in reversed(routine.args):
        if not arg.is_input:
            continue

        if allow_defaults and arg.transform.cpp_default:
            default_str = arg.transform.cpp_default
        else:
            default_str = ""
            allow_defaults = False

        cpp_type = _get_py_routine_arg_type(arg)

        # When optional_ref<T> is converted to T*, fix the default from std::nullopt to nullptr
        if cpp_type.endswith("*") and default_str == " = std::nullopt":
            default_str = " = nullptr"

        specs.append(f"{cpp_type} {arg.c_name}{default_str}")

    return list(reversed(specs))


def get_py_routine_decl(
    routine: FortranRoutine, return_type: str, python_name: str, defaults: bool, namespace: bool
):
    """Declaration of special py_ routine (wrapped due to immutable inout has its own definition)"""

    decl_args = ", ".join(_get_py_routine_decl_spec(routine, defaults))

    routine_and_args = f"{python_name}({decl_args})"

    if routine.cpp_namespace and namespace:
        routine_and_args = "::".join((routine.cpp_namespace, routine_and_args))

    return f"{return_type} {routine_and_args}"


def generate_py_routine_wrapper(routine: FortranRoutine) -> list[str]:
    assert routine.docstring is not None
    lines = []

    name = snake_to_camel(routine.name)
    py_name = f"Py{name}"

    args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]
    lines.append(
        get_py_routine_decl(
            routine,
            return_type=py_name,
            python_name=f"python_{routine.name}",
            defaults=True,
            namespace=False,
        )
    )
    lines.append("{")

    immut_args = [arg for arg in args if arg.arg.is_python_immutable]
    outputs = [arg for arg in args if arg.arg.intent == "out"]

    if outputs:
        res = "auto _result = "
    else:
        res = ""

    lines.append(f"  {res}{routine.cpp_namespace}::{routine.overloaded_name}(")

    def get_call_arg(arg: CppWrapperArgument):
        if arg.arg.transform.is_optional_ref and arg.arg.is_python_immutable:
            return f"make_opt_ref({arg.arg.c_name})"
        if arg.arg.transform.is_optional_ref:
            return f"ptr_to_opt_ref({arg.arg.c_name})"
        return arg.arg.c_name

    call_args = [get_call_arg(arg) for arg in args if arg.arg.intent in ("in", "inout")]
    lines.append(", ".join(call_args))
    lines.append(");")

    def get_output(arg: CppWrapperArgument):
        name = arg.decl_arg_name
        if arg.arg.is_optional:
            return f"{name}"
        return name

    local_outputs = ", ".join(get_output(arg) for arg in immut_args)

    if not outputs:
        lines.append(f"  auto py_result {{ {py_name} {{ {local_outputs} }} }};")
    else:
        lines.append(f"  auto py_result {{ {py_name}{{_result, {local_outputs} }} }};")

    lines.append("  return py_result;")
    lines.append("}")

    return lines


def _has_optional_ref_args(routine: FortranRoutine) -> bool:
    """Check if routine has optional_ref args that need lambda wrapping for nanobind."""
    return any(arg.transform.is_optional_ref for arg in routine.args if arg.is_input)


def _generate_optional_ref_lambda(routine: FortranRoutine) -> str:
    """Generate an inline lambda that converts T* to optional_ref<T> for nanobind."""
    args = routine.wrapper_args
    lambda_params = []
    call_args = []
    for arg in args:
        if not arg.is_input:
            continue
        cpp_type = arg.transform.cpp_type
        if arg.transform.is_optional_ref:
            inner = remove_optional(cpp_type)
            lambda_params.append(f"{inner} *{arg.c_name}")
            call_args.append(f"ptr_to_opt_ref({arg.c_name})")
        else:
            lambda_params.append(f"{cpp_type} {arg.c_name}")
            call_args.append(arg.c_name)

    param_str = ", ".join(lambda_params)
    call_str = ", ".join(call_args)
    fn = f"{routine.cpp_namespace}::{routine.overloaded_name}"

    # For overloaded functions, disambiguate with a typed function pointer
    input_arg_types = [arg.transform.cpp_type for arg in args if arg.is_input]
    fptr_args = ", ".join(input_arg_types)
    ret_type = routine.cpp_return_type or "void"
    # Assign the correct overload to a local function pointer, then call it
    fptr_type = f"{ret_type} (*)({fptr_args})"
    return (
        f"[]({param_str}) {{\n"
        f"        auto fn = static_cast<{fptr_type}>(&{fn});\n"
        f"        return fn({call_str});\n"
        f"    }}"
    )


def generate_routine_pybind_def(routine: FortranRoutine, overloads: list[FortranRoutine]) -> list[str]:
    assert routine.docstring is not None
    lines = []

    lines.append("m.def(")
    lines.append(f'  "{routine.overloaded_name}",')

    args = routine.wrapper_args

    # nanobind cannot bind optional_ref<T> (std::optional<std::reference_wrapper<T>>).
    # For routines not already wrapped, generate a lambda that takes T* and converts.
    needs_lambda = not routine.needs_python_wrapper and _has_optional_ref_args(routine)

    if needs_lambda:
        lines.append(f"{_generate_optional_ref_lambda(routine)},")
    elif routine.needs_python_wrapper:
        routine_ref = f"&python_{routine.name}"
        if overloads:
            arg_types = [_get_py_routine_arg_type(arg) for arg in args if arg.is_input]
            overload_args = ", ".join(arg_types)
            lines.append(f"nb::overload_cast<{overload_args}>({routine_ref}),")
        else:
            lines.append(f"{routine_ref},")
    else:
        routine_ref = f"&{routine.cpp_namespace}::{routine.overloaded_name}"
        if overloads:
            arg_types = [arg.transform.cpp_type for arg in args if arg.is_input]
            overload_args = ", ".join(arg_types)
            lines.append(f"nb::overload_cast<{overload_args}>({routine_ref}),")
        else:
            lines.append(f"{routine_ref},")

    for arg in args:
        if arg.is_input:
            if arg.is_optional:
                lines.append(f'nb::arg("{arg.python_name}") = nb::none(),')
            else:
                lines.append(f'nb::arg("{arg.python_name}"),')

    doc = routine.docstring.to_numpy_docstring(args)
    lines.append(rf'R"""({doc})"""')
    lines.append(");")

    return lines


def generate_py_routine_wrappers(routines: dict[str, FortranRoutine]) -> tuple[str, str]:
    structures = []
    code = []
    to_wrap = [routine for routine in routines.values() if routine.usable and routine.needs_python_wrapper]
    for routine in sorted_routines(to_wrap):
        structures.extend(generate_py_routine_return_value_struct(routine))
        code.extend(generate_py_routine_wrapper(routine))

    return "\n".join(structures), "\n".join(code)


def generate_py_routine_defs(routines: dict[str, FortranRoutine]):
    code = []

    to_wrap = [routine for routine in routines.values() if routine.usable]
    for routine in sorted_routines(to_wrap):
        overloads = [
            rt for rt in to_wrap if rt.overloaded_name == routine.overloaded_name and rt is not routine
        ]
        code.extend(generate_routine_return_value_wrapper(routine))
        code.extend(generate_routine_pybind_def(routine, overloads))

    return "\n".join(code)


def generate_pybmad_header(template: str) -> str:
    forward_decls = ["namespace Pybmad {"]

    # for struct in structs:
    #     forward_decls.append(f"void init_{struct.f_name}(py::module &, py::class_<{struct.cpp_class}> &);")

    # for struct in structs:
    #     forward_decls.append(f"std::string to_string(const {struct.cpp_class}& self);")
    forward_decls.append("} // namespace Pybmad")
    tpl = string.Template(template.replace("// ${", "${"))
    return tpl.substitute(forward_declarations="\n".join(forward_decls))


def generate_member_property_binding(cpp_class: str, arg: Argument, tpl) -> str:
    """
    Generate the nanobind property binding for a single struct member.

    Parameters
    ----------
    cpp_class : str
        Name of the C++ proxy class the member belongs to.
    arg : Argument
        The struct member to bind.
    tpl : ProxyTemplate
        Proxy template for ``arg.full_type``; used to decide whether a Fortran
        setter (and thus a read-write property) exists.
    """
    comment = arg.comment.replace('"', "'") if arg.comment else ""
    comment = comment.replace("\n", "\\n")
    docstring = f', "{comment}"' if comment else ""

    getter = f"&{cpp_class}::{arg.c_name}"

    # The getter returns a container/proxy by value that points into the
    # parent's Fortran memory, so the result must keep the parent (self) alive.
    # reference_internal (nanobind's default getter policy) drops its keep-alive
    # when a by-value return is promoted to move, so attach an explicit
    # keep_alive<0, 1> post-call hook instead.
    if tpl.fortran_setter:
        keepalive = ", nb::for_getter(nb::keep_alive<0, 1>())" if arg.needs_python_keepalive else ""
        setter = f"&{cpp_class}::set_{arg.c_name}"
        return f'.def_prop_rw("{arg.python_name}", {getter}, {setter}{keepalive}{docstring})'

    keepalive = ", nb::keep_alive<0, 1>()" if arg.needs_python_keepalive else ""
    return f'.def_prop_ro("{arg.python_name}", {getter}{keepalive}{docstring})'


def generate_pybmad_struct_code(struct: CodegenStructure, used_array_dims: set[int]) -> list[str]:
    code_lines = [""]
    code_lines.append("// =============================================================================")
    code_lines.append(f"// {struct.f_name}")
    code_lines.append(f"void init_{struct.f_name}(nb::module_ &m, nb::class_<{struct.cpp_class}> &cls) {{")

    ctor_entries: list[tuple[Argument, str]] = []

    for arg in struct.arg:
        (ctor_type, _ctor_body) = _generate_proxy_constructor_arg(struct, arg)
        if ctor_type is not None:
            ctor_entries.append((arg, ctor_type))

    has_optional_ref = any(ct.startswith("optional_ref") for _, ct in ctor_entries)

    if ctor_entries:
        if has_optional_ref:
            # nanobind cannot bind optional_ref<T> (std::optional<std::reference_wrapper<T>>)
            # in nb::init<>(), so use a placement-new __init__ lambda that takes T* instead.
            lambda_params = [f"{struct.cpp_class} *self"]
            call_args = []
            nb_args = []
            for arg, ctor_type in ctor_entries:
                if ctor_type.startswith("optional_ref"):
                    inner = remove_optional(ctor_type)
                    lambda_params.append(f"{inner} *{arg.c_name}")
                    call_args.append(f"ptr_to_opt_ref({arg.c_name})")
                else:
                    lambda_params.append(f"{ctor_type} {arg.c_name}")
                    call_args.append(arg.c_name)
                nb_args.append(f'nb::arg("{arg.python_name}") = nb::none()')

            params_str = ", ".join(lambda_params)
            call_str = ", ".join(call_args)
            args_str = ", ".join(nb_args)

            code_lines.append('    cls.def("__init__",')
            code_lines.append(f"        []({params_str}) {{")
            code_lines.append(f"            new (self) {struct.cpp_class}({call_str});")
            code_lines.append("        },")
            code_lines.append(f"        {args_str})")
        else:
            types_str = ", ".join(ct for _, ct in ctor_entries)
            args_str = ", ".join(f'nb::arg("{arg.python_name}") = nb::none()' for arg, _ in ctor_entries)

            # cls.def(nb::init<T1, T2>(), nb::arg("x")=none, nb::arg("y")=none);
            code_lines.append(f"    cls.def(nb::init<{types_str}>(),")
            code_lines.append(f"        {args_str})")
    else:
        code_lines.append("    cls.def(nb::init<>())")

    for arg in struct.arg:
        if not arg.is_component:
            continue

        try:
            tpl = proxy_templates[arg.full_type]
        except KeyError:
            code_lines.append(f"        // {arg.full_type} {arg.c_name} proxy support missing")
            continue

        code_lines.append(f"        {generate_member_property_binding(struct.cpp_class, arg, tpl)}")

    if 1 in used_array_dims:
        container_cls = f"{struct.cpp_class}Alloc1D"
        code_lines.append(
            f'      .def_static("new_array1d", [](int sz) {{ return {container_cls}(sz); }}, '
            f'nb::arg("sz") = 0)'
        )
        code_lines.append(
            f'      .def_static("new_array1d_bounds", [](int lbound, int ubound) {{ auto cnt = {container_cls}(); cnt.resize_bounds(lbound, ubound); return cnt; }}, '
            f'nb::arg("lbound"), nb::arg("ubound"))'
        )

    # TODO json
    # code.append(f'.def("to_json", &instance_to_json<{struct.cpp_class}>)')
    # code.append(
    #     f'.def("to_bson", [](const {struct.cpp_class} &self){{ json j; std::vector<std::uint8_t> v = json::to_bson(self); return v; }})'
    # )
    # code.append(
    #     f'.def("to_msgpack", [](const {struct.cpp_class} &self){{ json j; std::vector<std::uint8_t> v = json::to_msgpack(self); return v; }})'
    # )

    code_lines.append(
        textwrap.dedent(f"""
            .def("__repr__", [](const {struct.cpp_class} &self){{
                return to_string(self);
            }})
            """)
    )
    code_lines.append(
        textwrap.dedent(f"""
            .def("__copy__", [](const {struct.cpp_class} &self){{
                return {struct.cpp_class}(self);  // under-the-hood fortran copy
            }})
            .def("__deepcopy__", [](const {struct.cpp_class} &self, nb::dict& memo){{
                return {struct.cpp_class}(self);
            }})
            .def("__eq__", [](const {struct.cpp_class} &self, const {struct.cpp_class} &other){{
                return self.get_fortran_ptr() == other.get_fortran_ptr();
            }}, nb::is_operator())
            .def("__hash__", [](const {struct.cpp_class} &self){{
                return std::hash<std::uintptr_t>{{}}(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
            }})
            """)
    )
    code_lines.append("        ;")
    code_lines.append("")

    for n in SUPPORTED_ARRAY_DIMS:
        if n in used_array_dims:
            t_array = f"{struct.cpp_class}Array{n}D"
            t_alloc = f"{struct.cpp_class}Alloc{n}D"
            t_python_array = f"{struct.python_class_name}Array{n}D"
            t_python_alloc = f"{struct.python_class_name}Alloc{n}D"
            if n == 1:
                code_lines.append(
                    f'    bind_1d_type_array_pair<{t_array}, {t_alloc}>(m, "{t_python_array}", "{t_python_alloc}");'
                )
            else:
                code_lines.append(f'    bind_FTypeArrayND<{t_array}>(m, "{t_python_array}");')
        else:
            code_lines.append(f"    // {n}D {struct.cpp_class} arrays are not used in structs/routines")

    code_lines.append("}")

    return code_lines


def _group_structures_by_char(
    structs: list[CodegenStructure],
) -> dict[str, list[CodegenStructure]]:
    """Group structures by the first letter of their Fortran name."""
    structs_by_char = {}
    for st in structs:
        ch = st.f_name[0].lower()
        structs_by_char.setdefault(ch, []).append(st)
    return structs_by_char


def _submodule_var(submodule: str) -> str:
    """C++ variable name holding the pybind submodule for ``submodule``."""
    return f"m_{submodule}"


def _group_routines_by_source_and_char(
    routines_by_name: dict[str, FortranRoutine],
) -> dict[str, dict[str, list[FortranRoutine]]]:
    """Group routines by source namespace and first letter."""
    routines_map = {}

    for routine in routines_by_name.values():
        if not routine.usable:
            continue

        src = routine.cpp_namespace
        char = routine.name[0].lower()

        routines_map.setdefault(src, {}).setdefault(char, [])
        routines_map[src][char].append(routine)

    return routines_map


def _generate_structure_files(
    files: dict[pathlib.Path, str],
    structs_by_char: dict[str, list[CodegenStructure]],
    array_usage: dict[str, set],
) -> list[str]:
    """
    Generate the split C++ files (structs_{letter}) for structures.
    Returns the list of generated headers and initialization function calls.

    Parameters
    ----------
    files : dict[pathlib.Path, str]
    structs_by_char : dict[str, list[CodegenStructure]]
        Structures grouped by their first character.
    all_structs : list[CodegenStructure]
        List of all available structures.
    routines_by_name : dict[str, FortranRoutine]
        Dictionary of Fortran routines keyed by name.
    array_usage : dict[str, set]
        Dictionary tracking array dimension usage per structure.

    Returns
    -------
    list[str]
        Header includes
    """
    headers: list[str] = ["#pragma once"]

    for char, char_structs in sorted(structs_by_char.items()):
        char_structs.sort(key=lambda st: (st.module, st.f_name))

        header_name = f"structs_{char}.hpp"
        header_path = PYBMAD_INCLUDE / "pybmad" / "generated" / header_name

        header_decls = []
        for st in char_structs:
            header_decls.append(f"void init_{st.f_name}(nb::module_ &m, nb::class_<{st.cpp_class}> &class_);")

        newline = "\n"
        files[header_path] = textwrap.dedent(f"""\
            #pragma once
            #include <nanobind/nanobind.h>
            #include "bmad/generated/proxy.hpp"
            #include "pybmad/generated/structs.hpp"
            namespace nb = nanobind;

            using namespace Bmad;

            // Per-struct init functions
            {newline.join(header_decls)}
        """)

        headers.append(f'#include "pybmad/generated/{header_name}"')

        src_lines = [
            f'#include "pybmad/generated/{header_name}"',
            '#include "bmad/generated/proxy.hpp"',
            '#include "bmad/generated/to_string.hpp"',
            '#include "bmad/to_string.hpp"',
            '#include "pybmad/arrays.hpp"',
            "#include <cstdint>",
            "#include <functional>",
            '#include "pybmad/util.hpp"',
            "",
            "using namespace Pybmad;",
            "namespace nb = nanobind;",
            "",
        ]

        for st in char_structs:
            lines = generate_pybmad_struct_code(st, array_usage.get(st.f_name.lower(), set()))
            src_lines.extend(lines)

        src_path = PYBMAD_SRC / "generated" / f"structs_{char}.cpp"
        files[src_path] = "\n".join(src_lines)

    return headers


def _generate_routine_files(
    files: dict[pathlib.Path, str],
    routines_map: dict[str, dict[str, list[FortranRoutine]]],
    submodule_map: dict[str, str],
) -> tuple[list[str], list[str]]:
    """
    Generate the split C++ files ({source}_routines_{letter}) for routines.
    Returns the list of generated headers and initialization function calls.

    Routines whose C++ namespace has a ``submodule_map`` entry are bound into
    that pybind submodule; the rest bind onto the top-level module ``m``.
    """
    headers: list[str] = ["#pragma once"]
    init_calls: list[str] = []

    for src, chars_dict in sorted(routines_map.items()):
        submodule = submodule_map.get(src, "")
        module_var = _submodule_var(submodule) if submodule else "m"
        for char, routine_list in sorted(chars_dict.items()):
            if not routine_list:
                continue

            subset_routines = {r.name: r for r in routine_list}

            base_name = f"{src}_routines_{char}"
            header_name = f"{base_name}.hpp"
            src_name = f"{base_name}.cpp"
            init_fn_name = f"init_{src}_routines_{char}"

            wrapper_structs, wrappers_code = generate_py_routine_wrappers(subset_routines)

            header_path = PYBMAD_INCLUDE / "pybmad" / "generated" / header_name
            files[header_path] = textwrap.dedent(f"""\
                #pragma once
                #include <nanobind/nanobind.h>
                #include <nanobind/stl/complex.h>
                #include <nanobind/stl/optional.h>
                #include <nanobind/stl/vector.h>
                #include <nanobind/stl/string.h>

                #include "pybmad/arrays.hpp"
                #include "pybmad/util.hpp"

                namespace nb = nanobind;

                void {init_fn_name}(nb::module_ &m);

                {wrapper_structs}
            """)

            headers.append(f'#include "pybmad/generated/{header_name}"')
            init_calls.append(f"{init_fn_name}({module_var});")

            defs_block = generate_py_routine_defs(subset_routines)
            # defs_block_indented = textwrap.indent(defs_block, "  ")

            cpp_content = textwrap.dedent(f"""\
                #include "pybmad/generated/{header_name}"

                namespace nb = nanobind;
                using namespace nanobind::literals;
                using namespace Pybmad;

                {wrappers_code}

                void {init_fn_name}(nb::module_ &m) {{
                {defs_block}
                }}
            """)

            files[PYBMAD_SRC / "generated" / src_name] = cpp_content

    return headers, init_calls


def _generate_struct_init(structs: list[CodegenStructure]):
    src_lines = []
    for st in structs:
        src_lines.append(
            f"    auto py_{st.python_class_name} = nb::class_<{st.cpp_class}>(m, "
            f'"{st.python_class_name}", "Fortran struct: {st.f_name}");'
        )

    for st in structs:
        src_lines.append(f"    init_{st.f_name}(m, py_{st.python_class_name});")
    return src_lines


def _generate_main_module_file(
    files: dict[pathlib.Path, str],
    enums: dict[str, dict[str, EnumValue]],
    struct_inits: list[str],
    routine_inits: list[str],
    submodule_decls: list[str],
) -> None:
    """Generate the main _pybmad.cpp file aggregating all bindings."""
    template_text = (CODEGEN_ROOT / "pybind_mod.tpl.cpp").read_text()
    template = string.Template(template_text.replace("// ${", "${"))

    newline = "\n"
    core_body = textwrap.dedent(f"""\
        m.doc() = "pybmad";

        bind_standard_arrays(m);

        // Structures
        {newline.join(struct_inits)}

        // Hand-written bindings
        init_common_structs(m);

        // Routine submodules (one per C++ namespace)
        {newline.join(f"    {s}" for s in submodule_decls)}

        // Routine initializers
        {newline.join(f"    {s}" for s in routine_inits)}
    """).strip()

    inclusions = textwrap.dedent("""
        #include "pybmad/generated/structs.hpp"
        #include "pybmad/generated/routines.hpp"
    """)

    substituted = template.substitute(
        pybind11_routine_wrappers=inclusions,
        pybind11_definitions=core_body,
    )

    files[PYBMAD_SRC / "generated" / "_pybmad.cpp"] = substituted
    files[PYBMAD_LIB / "_enums.py"] = generate_enum_wrapper_code(enums)


def _pybmad_version() -> str:
    """Version for the generated package.

    Defaults to the ``UPSTREAM_TAG`` env var (set by the track-upstream
    workflow to the Bmad release being wrapped); otherwise falls back to
    ``git describe`` of this repo.
    """
    tag = os.environ.get("UPSTREAM_TAG")
    if tag:
        return tag
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--no-dirty"],
            cwd=CODEGEN_ROOT,
            text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "0.0.0"


def generate_init_dot_py(
    config: CodegenConfig,
    structs: list[CodegenStructure],
    routines_by_name: dict[str, FortranRoutine],
    enums: dict[str, dict[str, EnumValue]],
    array_usage: dict[str, set[int]],
) -> str:
    imports = []
    all_ = []

    # Map each C++ namespace that has usable routines to its Python submodule
    # name, preserving config order for the submodule import block.
    present_namespaces = {rt.cpp_namespace for rt in routines_by_name.values() if rt.usable}
    namespace_to_submodule: dict[str, str] = {}
    namespace_top_level: dict[str, bool] = {}
    submodule_names: list[str] = []
    for r in config.routines:
        if r.python_submodule and r.cpp_namespace in present_namespaces:
            namespace_to_submodule[r.cpp_namespace] = r.python_submodule
            namespace_top_level[r.cpp_namespace] = r.python_top_level
            if r.python_submodule not in submodule_names:
                submodule_names.append(r.python_submodule)

    def add_name(name: str, mod: str = config.python_module_name):
        import_line = f"from .{mod} import {name}"
        all_line = f'    "{name}",'
        if import_line not in imports:
            imports.append(import_line)
            all_.append(all_line)

    imports.append("# Globals")
    all_.append("    # Globals")
    for name in config.python_imports:
        add_name(name)

    imports.append("# Classes")
    all_.append("    # Classes")

    for struct in structs:
        add_name(struct.python_class_name)
        for n in sorted(array_usage.get(struct.f_name.lower(), [])):
            if n > 0:
                add_name(f"{struct.python_class_name}Array{n}D")
                if n == 1:
                    add_name(f"{struct.python_class_name}Alloc{n}D")

    submodule_lines = []
    for sub in submodule_names:
        add_name(sub)
        submodule_lines.append(f'_sys.modules[f"{{__name__}}.{sub}"] = {sub}')

    imports.append("")
    imports.append("# Functions")
    all_.append("")
    all_.append("    # Functions")
    for _, rt in sorted(routines_by_name.items(), key=lambda item: item[0]):
        if not rt.usable:
            continue
        sub = namespace_to_submodule.get(rt.cpp_namespace)
        if sub and not namespace_top_level.get(rt.cpp_namespace, True):
            # Reachable only as pybmad.<sub>.<name>, not at the top level.
            continue
        if sub:
            # Top-level alias preserving the historic flat API; canonical
            # binding lives in the submodule.
            alias_line = f"{rt.overloaded_name} = {sub}.{rt.overloaded_name}"
            if alias_line not in imports:
                imports.append(alias_line)
                all_.append(f'    "{rt.overloaded_name}",')
        else:
            add_name(rt.overloaded_name)

    imports.append("")
    imports.append("# Enums")
    all_.append("")
    all_.append("    # Enums")
    for per_file_enums in enums.values():
        for enum in per_file_enums.values():
            add_name(enum.name, "_enums")

    add_name("EleAttribute", "_enums")
    add_name("EleKey", "_enums")

    nl = "\n"
    version = _pybmad_version()
    return f"""
from __future__ import annotations
import sys as _sys

__version__ = "{version}"

{nl.join(imports)}

# Submodules
{nl.join(submodule_lines)}

del _sys

__all__ = [
{nl.join(all_)}
]
    """


def expected_stub_files(
    config: CodegenConfig,
    routines_by_name: dict[str, FortranRoutine],
) -> list[pathlib.Path]:
    """
    Paths of the ``.pyi`` stubs that ``nanobind.stubgen`` will emit at build time.
    """
    namespaces = {rt.cpp_namespace for rt in routines_by_name.values() if rt.usable}
    submodule_names: set[str] = {
        r.python_submodule for r in config.routines if r.python_submodule and r.cpp_namespace in namespaces
    }

    package_dir = PYBMAD_LIB / config.python_module_name
    return [
        package_dir / "__init__.pyi",
        *[package_dir / f"{sub}.pyi" for sub in sorted(submodule_names)],
    ]


def generate_pybmad(
    config: CodegenConfig,
    structs: list[CodegenStructure],
    routines_by_name: dict[str, FortranRoutine],
    enums: dict[str, dict[str, EnumValue]],
) -> dict[pathlib.Path, str]:
    """
    Generate split pybind11 bindings for Bmad.

    Parameters
    ----------
    structs : list[CodegenStructure]
        List of structure definitions to bind.
    routines_by_name : dict[str, FortranRoutine]
        Dictionary mapping routine names to their definitions.
    enums : dict[str, dict[str, EnumValue]] | None
        Optional parsed enums. If None, they are parsed.

    Returns
    -------
    dict[pathlib.Path, str]
        A dictionary mapping file output paths to their generated string content.
    """
    files: dict[pathlib.Path, str] = {}
    array_usage = struct_array_usage_dimensions(routines_by_name, structs)

    init_dot_py = generate_init_dot_py(config, structs, routines_by_name, enums, array_usage)

    structs_by_char = _group_structures_by_char(structs)
    routines_map = _group_routines_by_source_and_char(routines_by_name)

    submodule_map = {r.cpp_namespace: r.python_submodule for r in config.routines if r.python_submodule}

    struct_headers = _generate_structure_files(files, structs_by_char, array_usage)
    struct_inits = _generate_struct_init(structs)

    routine_headers, routine_inits = _generate_routine_files(files, routines_map, submodule_map)

    # Declare one submodule per namespace that actually has routines.
    submodule_decls = [
        f"auto {_submodule_var(submodule_map[src])} = "
        f'm.def_submodule("{submodule_map[src]}", "{src} routines");'
        for src in sorted(routines_map)
        if src in submodule_map
    ]

    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "structs.hpp"] = "\n".join(struct_headers)
    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "routines.hpp"] = "\n".join(routine_headers)
    files[PYBMAD_ROOT / "pybmad" / "__init__.py"] = init_dot_py

    _generate_main_module_file(files, enums, struct_inits, routine_inits, submodule_decls)

    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "init.hpp"] = generate_pybmad_header(
        template=(CODEGEN_ROOT / "pybind.tpl.hpp").read_text(),
    )

    return files
