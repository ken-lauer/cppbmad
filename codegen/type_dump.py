"""
Generate all binding-layer code for a single Fortran type specification.

Helper for debugging all of the layers here.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import textwrap
from string import Template

from .arg import Argument
from .context import CodegenConfig, ConfigContext, config_context
from .cpp import CppWrapperArgument
from .fortran import FortranWrapperArgument
from .paths import CODEGEN_ROOT
from .proxy import (
    generate_accessor_code,
    split_signature,
)
from .proxy import (
    templates as proxy_templates,
)
from .routines import FortranRoutine, RoutineArg
from .struct_parser.parser import StructureMember, TypeInformation
from .struct_parser.util import FileLine
from .types import (
    ALLOC,
    NOT,
    PTR,
    STRUCT,
    ArgumentType,
    FullType,
    Intent,
    PointerType,
    get_type_transform,
)
from .util import struct_to_proxy_class_name

EXAMPLE_STRUCT = "example_struct"
EXAMPLE_CLASS = "ExampleStruct"
EXAMPLE_MEMBER = "member"
EXAMPLE_DERIVED_KIND = "branch_struct"


def _ensure_context() -> None:
    """Ensure a ConfigContext is available for codegen helpers that need it."""
    try:
        config_context.get()
    except LookupError:
        params = CodegenConfig.from_file(CODEGEN_ROOT / "default.toml")
        config_context.set(ConfigContext(params=params))


def _make_type_info(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
) -> TypeInformation:
    """Build a TypeInformation for a synthetic struct member."""
    dimension: str | None = None
    if dim > 0:
        if shape:
            dimension = ",".join(shape)
        elif ptr in (ALLOC, PTR):
            dimension = ",".join([":"] * dim)
        else:
            dimension = ",".join(["6"] * dim)

    actual_type = base_type
    if base_type == "integer8":
        actual_type = "integer"
    elif base_type == "real16":
        actual_type = "real"

    return TypeInformation(
        type=actual_type,
        kind=kind or None,
        pointer=ptr == PTR,
        allocatable=ptr == ALLOC,
        dimension=dimension,
    )


def _make_member(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
) -> StructureMember:
    """Build a synthetic StructureMember."""
    type_info = _make_type_info(base_type, kind, ptr, dim, shape)
    definition = type_info.to_fortran_declaration() + f" :: {EXAMPLE_MEMBER}"
    return StructureMember(
        definition=definition,
        type_info=type_info,
        name=EXAMPLE_MEMBER,
    )


def _make_argument(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
) -> Argument:
    """Build a synthetic Argument for struct member use."""
    member = _make_member(base_type, kind, ptr, dim, shape)

    if dim > 0:
        if shape:
            array = list(shape)
        elif ptr in (ALLOC, PTR):
            array = [":"] * dim
        else:
            array = ["6"] * dim
    else:
        array = []

    actual_kind = kind
    if base_type == "integer8":
        actual_kind = ""
    elif base_type == "real16":
        actual_kind = "qp"

    return Argument(
        is_component=True,
        f_name=EXAMPLE_MEMBER,
        c_name=EXAMPLE_MEMBER,
        python_name=EXAMPLE_MEMBER,
        type=base_type,
        kind=actual_kind,
        pointer_type=ptr,
        array=array,
        member=member,
    )


def _make_routine_arg(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
    intent: Intent,
    is_optional: bool,
):
    """Build a synthetic RoutineArg for routine wrapper generation."""
    member = _make_member(base_type, kind, ptr, dim, shape)
    if is_optional:
        member.type_info = member.type_info.replace(optional=True)

    if dim > 0:
        if shape:
            array = list(shape)
        elif ptr in (ALLOC, PTR):
            array = [":"] * dim
        else:
            array = ["6"] * dim
    else:
        array = []

    actual_kind = kind
    if base_type == "integer8":
        actual_kind = ""
    elif base_type == "real16":
        actual_kind = "qp"

    return RoutineArg(
        is_component=True,
        f_name=f"f_{EXAMPLE_MEMBER}",
        c_name=EXAMPLE_MEMBER,
        python_name=EXAMPLE_MEMBER,
        type=base_type,
        kind=actual_kind,
        pointer_type=ptr,
        array=array,
        member=member,
        intent=intent,
        doc_is_optional=is_optional,
    )


def _make_routine():
    """Build a minimal synthetic FortranRoutine."""
    return FortranRoutine(
        filename=pathlib.Path("synthetic.f90"),
        name="example_routine",
        proc_type="subroutine",
        start_line=FileLine(pathlib.Path("synthetic.f90"), 1, "subroutine example_routine()"),
    )


def generate_struct_member_code(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
) -> dict[str, str]:
    """Generate all layers of code for a struct member.

    All code is produced by the real codegen machinery in proxy.py.
    Returns a dict mapping descriptive labels to code strings.
    """
    ft = FullType(base_type, dim, ptr)
    arg = _make_argument(base_type, kind, ptr, dim, shape)
    result: dict[str, str] = {}

    if ft not in proxy_templates:
        result["error"] = (
            f"FullType {ft} has no struct member proxy template.\n"
            f"This combination is not currently supported for struct members."
        )
        return result

    acc = generate_accessor_code(EXAMPLE_STRUCT, arg, ft, kind)
    result["Fortran getter (proxy_mod.f90)"] = acc["fortran_getter"].strip()
    if acc["fortran_setter"]:
        result["Fortran setter (proxy_mod.f90)"] = acc["fortran_setter"].strip()

    proxy_cls = EXAMPLE_CLASS
    attr_proxy_cls = struct_to_proxy_class_name(kind) if base_type == STRUCT and kind else ""

    # C++ extern "C" declarations
    extern_c_lines = [acc["cpp_get_decl"]]
    if acc["cpp_set_decl"]:
        extern_c_lines.append(acc["cpp_set_decl"])
    result['C++ extern "C" declarations (proxy.hpp)'] = "\n".join(extern_c_lines).strip()

    header_lines = []
    impl_lines = []
    for accessor_body in acc["cpp_get_accessors"]:
        if base_type == STRUCT and kind:
            accessor_body = Template(accessor_body).substitute(
                return_proxy_name=attr_proxy_cls,
            )
        sig, impl = split_signature(accessor_body, proxy_cls)
        header_lines.append(sig.strip())
        impl_lines.append(impl.strip())

    if acc["cpp_set_accessors"]:
        for accessor_body in acc["cpp_set_accessors"]:
            if base_type == STRUCT and kind:
                accessor_body = Template(accessor_body).substitute(
                    return_proxy_name=attr_proxy_cls,
                )
            sig, impl = split_signature(accessor_body, proxy_cls)
            header_lines.append(sig.strip())
            impl_lines.append(impl.strip())

    result["C++ class declarations (proxy.hpp)"] = "\n".join(header_lines)
    result["C++ implementation (proxy.cpp)"] = "\n\n".join(impl_lines)

    # Python binding (TODO duplicated)
    tpl = proxy_templates[ft]
    getter = f"&{EXAMPLE_CLASS}::{EXAMPLE_MEMBER}"
    if arg.needs_python_keepalive:
        getter_fn = f"nb::cpp_function({getter}, nb::keep_alive<0, 1>())"
    else:
        getter_fn = getter

    if tpl.fortran_setter:
        setter = f"&{EXAMPLE_CLASS}::set_{EXAMPLE_MEMBER}"
        result["Python binding (nanobind)"] = f'    .def_prop_rw("{EXAMPLE_MEMBER}", {getter_fn}, {setter})'
    else:
        result["Python binding (nanobind)"] = f'    .def_prop_ro("{EXAMPLE_MEMBER}", {getter_fn})'

    return result


def generate_routine_arg_code(
    base_type: ArgumentType,
    kind: str,
    ptr: PointerType,
    dim: int,
    shape: list[str] | None,
    intent: Intent,
    is_optional: bool,
) -> dict[str, str]:
    """Generate all layers of wrapper code for a routine argument.

    Uses the real FortranWrapperArgument and CppWrapperArgument classes.
    Returns a dict mapping descriptive labels to code strings.
    """
    _ensure_context()

    ft = FullType(base_type, dim, ptr)
    result: dict[str, str] = {}
    try:
        transform = get_type_transform(
            ft,
            intent=intent,
            is_optional=is_optional,
            kind=kind,
            is_dynamic_array=ptr in (ALLOC, PTR) and dim > 0,
        )
    except Exception as ex:
        result["error"] = f"get_type_transform failed: {ex}"
        return result

    cpp_type = transform.cpp_type
    cpp_return = transform.cpp_return_type
    cpp_call = transform.cpp_call_fortran_type
    fortran_type = transform.fortran_type
    fortran_native = transform.fortran_native_type

    if base_type == STRUCT and kind:
        proxy_name = struct_to_proxy_class_name(kind)
        cpp_type = cpp_type.replace("PROXYCLS", proxy_name)
        cpp_return = cpp_return.replace("PROXYCLS", proxy_name)

    result["Transform summary"] = textwrap.dedent(f"""\
        Fortran native type:  {fortran_native}
        Fortran iso_c type:   {fortran_type}
        C++ argument type:    {cpp_type}
        C++ declare type:     {transform.cpp_declare_type}
        C++ extern "C" type:  {cpp_call}
        C++ return type:      {cpp_return}
        C++ default value:    {transform.cpp_default or "(none)"}""")

    try:
        routine_arg = _make_routine_arg(base_type, kind, ptr, dim, shape, intent, is_optional)
        routine = _make_routine()
    except Exception as ex:
        result["note"] = f"Could not construct synthetic routine arg: {ex}"
        return result

    try:
        f_lines: list[str] = []
        f_handler = FortranWrapperArgument.from_arg(routine_arg, routine, f_lines, have_err_flag=False)

        f_sections: list[str] = []
        decls = f_handler.get_declarations()
        if decls:
            f_sections.append("! Declarations")
            f_sections.extend(decls)
        input_conv = f_handler.get_input_conversion()
        if input_conv:
            f_sections.append("")
            f_sections.append("! Input conversion (C -> Fortran)")
            f_sections.extend(input_conv)
        output_conv = f_handler.get_output_conversion()
        if output_conv:
            f_sections.append("")
            f_sections.append("! Output conversion (Fortran -> C)")
            f_sections.extend(output_conv)
        f_sections.append("")
        f_sections.append(f"! Call argument name: {f_handler.get_call_arg_name()}")

        result["Fortran wrapper (routines_mod.f90)"] = "\n".join(f_sections)
    except Exception as ex:
        result["Fortran wrapper error"] = str(ex)

    try:
        cpp_handler = CppWrapperArgument.from_arg(routine_arg)

        cpp_sections: list[str] = []
        cpp_sections.append(f"// Declaration: {cpp_handler.cpp_decl}")
        pre = cpp_handler.pre_call_lines()
        if pre:
            cpp_sections.append("")
            cpp_sections.append("// Pre-call setup")
            cpp_sections.extend(pre)
        cpp_sections.append("")
        cpp_sections.append(f"// Fortran call argument: {cpp_handler.call_argument()}")
        post = cpp_handler.post_call_lines()
        if post:
            cpp_sections.append("")
            cpp_sections.append("// Post-call")
            cpp_sections.extend(post)

        struct_decl = cpp_handler.struct_decl()
        if struct_decl:
            type_, name = struct_decl
            cpp_sections.append("")
            cpp_sections.append(f"// Output struct field: {type_} {name}")

        result["C++ wrapper (routines.cpp)"] = "\n".join(cpp_sections)
    except Exception as ex:
        result["C++ wrapper error"] = str(ex)

    return result


def render_markdown(
    spec_desc: str,
    struct_code: dict[str, str] | None,
    routine_codes: list[tuple[str, dict[str, str]]],
) -> str:
    """Render generated code as a Markdown document."""
    lines: list[str] = []

    def add(text: str = "") -> None:
        lines.append(text)

    add(f"# Fortran/C++/Python wrappers for `{spec_desc}`")
    add()

    if struct_code is not None:
        add("## Struct Member")
        add()
        if "error" in struct_code:
            add(f"**Not supported:** {struct_code['error']}")
            add()
        else:
            for label, code in struct_code.items():
                lang = _lang_for_label(label)
                add(f"### {label}")
                add()
                add(f"```{lang}")
                add(code)
                add("```")
                add()

    if routine_codes:
        add("## Routine Arguments")
        add()
        for heading, code in routine_codes:
            add(f"### {heading}")
            add()
            if "error" in code:
                add(f"**Not supported:** {code['error']}")
                add()
            else:
                for label, content in code.items():
                    lang = _lang_for_label(label)
                    add(f"#### {label}")
                    add()
                    add(f"```{lang}")
                    add(content)
                    add("```")
                    add()

    return "\n".join(lines)


def _lang_for_label(label: str) -> str:
    ll = label.lower()
    if "fortran" in ll:
        return "fortran"
    if "c++" in ll or "cpp" in ll or "extern" in ll:
        return "cpp"
    if "python" in ll or "pybind" in ll:
        return "cpp"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate all binding-layer code for a single Fortran type specification.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              %(prog)s real                          # scalar real
              %(prog)s integer --dim 1 --alloc       # 1D allocatable integer
              %(prog)s type --kind branch_struct     # derived type scalar
              %(prog)s real --dim 2 --shape 6,6      # 2D fixed-size real
              %(prog)s real --optional               # optional inout scalar
              %(prog)s real --dim 1 --alloc
        """),
    )

    parser.add_argument(
        "base_type",
        choices=["real", "real16", "integer", "integer8", "logical", "complex", "character", "type"],
        help="Fortran base type",
    )
    parser.add_argument("--dim", type=int, default=0, help="Array dimensionality (0-3, default: 0)")
    parser.add_argument(
        "--shape",
        default=None,
        help="Comma-separated fixed dimensions (e.g. 6 or 6,6). Implies --dim if not set.",
    )

    ptr_group = parser.add_mutually_exclusive_group()
    ptr_group.add_argument("--ptr", action="store_true", help="Pointer type")
    ptr_group.add_argument("--alloc", action="store_true", help="Allocatable type")

    parser.add_argument(
        "--kind",
        default="",
        help="Type kind / derived type name (required for type, e.g. branch_struct)",
    )
    parser.add_argument("--optional", action="store_true", help="Mark routine argument as optional")
    parser.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")

    args = parser.parse_args()

    base_type: ArgumentType = args.base_type
    dim: int = args.dim
    shape: list[str] | None = None

    if args.shape:
        shape = args.shape.split(",")
        if dim == 0:
            dim = len(shape)

    ptr: PointerType = NOT
    if args.ptr:
        ptr = PTR
    elif args.alloc:
        ptr = ALLOC

    kind = args.kind
    if base_type == STRUCT and not kind:
        kind = EXAMPLE_DERIVED_KIND

    spec_parts = [base_type]
    if kind:
        spec_parts[-1] = f"type({kind})"
    if dim > 0:
        if shape:
            spec_parts.append(f"dim={','.join(shape)}")
        else:
            spec_parts.append(f"{dim}D")
    if ptr != NOT:
        spec_parts.append(ptr.lower())
    spec_desc = " ".join(spec_parts)

    struct_code: dict[str, str] | None = None
    struct_code = generate_struct_member_code(base_type, kind, ptr, dim, shape)

    routine_codes: list[tuple[str, dict[str, str]]] = []
    intents_to_show: list[Intent] = ["in", "out", "inout"]
    optionals_to_show = [False, True] if args.optional else [False]

    for intent in intents_to_show:
        for opt in optionals_to_show:
            label = f"intent({intent})"
            if opt:
                label += ", optional"
            code = generate_routine_arg_code(base_type, kind, ptr, dim, shape, intent, opt)
            routine_codes.append((label, code))

    output = render_markdown(spec_desc, struct_code, routine_codes)

    if args.output == "-":
        sys.stdout.write(output)
        if not output.endswith("\n"):
            sys.stdout.write("\n")
    else:
        pathlib.Path(args.output).write_text(output if output.endswith("\n") else output + "\n")


if __name__ == "__main__":
    main()
