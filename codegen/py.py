from __future__ import annotations

import pathlib
import string
import textwrap

from .arg import CodegenStructure
from .config import SUPPORTED_ARRAY_DIMS
from .cpp import CppWrapperArgument
from .enums import EnumValue, get_ele_attributes, get_ele_keys, parse_all_enums
from .paths import CODEGEN_ROOT, PYBMAD_INCLUDE, PYBMAD_SRC
from .proxy import templates as proxy_templates
from .routines import FortranRoutine, RoutineArg, is_python_immutable
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


def generate_enum_wrapper(clsname: str, enums: list[EnumValue]) -> str:
    code = [f'    py::native_enum<{clsname}>(m, "{clsname}", "enum.IntEnum") ']
    for attr in enums:
        doc = f', "{attr.comment}"' if attr.comment else ""
        code.append(f'        .value("{attr.name}", {clsname}::{attr.name}{doc})')
    code.append("        .export_values()")
    code.append("        .finalize();")
    return "\n".join(code)


def generate_enum_constants(all_enums: dict[str, dict[str, EnumValue]]) -> str:
    result = []
    for fn, enums in all_enums.items():
        result.append("")
        result.append(f"// Enums from {fn}")
        for enum in enums.values():
            if enum.comment:
                result.append(f"// {enum.comment}")
            pytype = "int_" if enum.type == "int" else "float_"
            result.append(f'm.attr("{enum.name}") = py::{pytype}(Bmad::{enum.name});')
    return "\n".join(result)


def generate_enum_wrapper_code(enums: dict[str, dict[str, EnumValue]]) -> str:
    return "\n".join(
        (
            generate_enum_wrapper("EleAttribute", get_ele_attributes(enums["bmad_struct.f90"])),
            generate_enum_wrapper("EleKey", get_ele_keys(enums["bmad_struct.f90"])),
            generate_enum_constants(enums),
        )
    )


def generate_routine_return_value_wrapper(routine: FortranRoutine) -> list[str]:
    immut_args = [arg for arg in routine.args if is_python_immutable(arg)]
    outputs = [*routine.outputs, *immut_args]

    if len(outputs) <= 1 and not immut_args:  # TODO: immut args even for 1 -> struct
        return []

    clsname = snake_to_camel(routine.name)
    if immut_args:
        full_clsname = f"Py{clsname}"
    else:
        full_clsname = routine.cpp_return_type
        clsname = full_clsname.split("::")[1]

    lines = []
    lines.append(
        f'    py::class_<{full_clsname}, std::unique_ptr<{full_clsname}>>(m, "{clsname}", "Fortran routine {routine.name} return value")'
    )
    for arg in outputs:
        lines.append(f'        .def_readonly("{arg.python_name}", &{full_clsname}::{arg.c_name})')

    lines.append(f'        .def("__len__", [](const {full_clsname} &) {{ return {len(outputs)}; }})')
    lines.append(f'        .def("__getitem__", [](const {full_clsname} &s, size_t i) -> py::object {{')
    lines.append(f"            if (i >= {len(outputs)}) throw py::index_error();")
    for i, arg in enumerate(outputs):
        lines.append(f"            if (i == {i}) return py::cast(s.{arg.c_name});")
    lines.append("            return py::none();")
    lines.append("        })")

    # lines.append(f'      .def("__repr__", [](const {full_clsname} &self){{ return to_string(self); }})')
    lines.append("        ;")
    return lines


def generate_py_routine_return_value_struct(routine: FortranRoutine) -> list[str]:
    args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]
    outputs = [arg for arg in args if arg.arg.intent == "out"]

    name = snake_to_camel(routine.name)
    py_name = f"Py{name}"
    orig_struct = routine.cpp_return_type

    immut_args = [arg for arg in args if is_python_immutable(arg.arg)]
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
    if is_python_immutable(arg):
        if cpp_type.startswith("optional_ref"):
            # optional_ref<T> -> std::optional<T>
            cpp_type = remove_optional(cpp_type)
            cpp_type = f"std::optional<{cpp_type}>"
        cpp_type = cpp_type.replace("&", "")
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

    immut_args = [arg for arg in args if is_python_immutable(arg.arg)]
    outputs = [arg for arg in args if arg.arg.intent == "out"]

    if outputs:
        res = "auto _result = "
    else:
        res = ""

    lines.append(f"  {res}{routine.cpp_namespace}::{routine.overloaded_name}(")

    def get_call_arg(arg: CppWrapperArgument):
        if arg.arg.transform.is_optional_ref and is_python_immutable(arg.arg):
            return f"make_opt_ref({arg.arg.c_name})"
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


def generate_routine_pybind_def(routine: FortranRoutine, overloads: list[FortranRoutine]) -> list[str]:
    assert routine.docstring is not None
    lines = []

    lines.append("m.def(")
    lines.append(f'  "{routine.overloaded_name}",')

    if routine.needs_python_wrapper:
        routine_ref = f"&python_{routine.name}"
    else:
        routine_ref = f"&{routine.cpp_namespace}::{routine.overloaded_name}"

    if overloads:
        if routine.needs_python_wrapper:
            arg_types = [_get_py_routine_arg_type(arg) for arg in routine.args if arg.is_input]
        else:
            arg_types = [arg.transform.cpp_type for arg in routine.args if arg.is_input]
        overload_args = ", ".join(arg_types)
        lines.append(f"py::overload_cast<{overload_args}>({routine_ref}),")

    else:
        lines.append(f"{routine_ref},")

    for arg in routine.args:
        if arg.is_input:
            if arg.is_optional:
                lines.append(f'py::arg("{arg.c_name}") = py::none(),')
            else:
                lines.append(f'py::arg("{arg.c_name}"),')

    doc = routine.docstring.to_numpy_docstring()
    lines.append(rf'R"""({doc})"""')
    lines.append(");")

    return lines


def generate_py_routine_wrappers(routines: dict[str, FortranRoutine]):
    code = []

    to_wrap = [routine for routine in routines.values() if routine.usable and routine.needs_python_wrapper]
    for routine in sorted_routines(to_wrap):
        code.extend(generate_py_routine_return_value_struct(routine))
        code.extend(generate_py_routine_wrapper(routine))

    return "\n".join(code)


def generate_py_routine_defs(routines: dict[str, FortranRoutine]):
    code = []

    to_wrap = [routine for routine in routines.values() if routine.usable]
    for routine in sorted_routines(to_wrap):
        overloads = [
            rt for rt in to_wrap if rt.overloaded_name == routine.overloaded_name and rt is not routine
        ]
        code.extend(generate_routine_pybind_def(routine, overloads))
        code.extend(generate_routine_return_value_wrapper(routine))

    return "\n".join(code)


def generate_pybmad_header(
    template: str,
    routines: dict[str, FortranRoutine],  # TODO # noqa: ARG001
    structs: list[CodegenStructure],
) -> str:
    forward_decls = ["namespace Pybmad {"]

    for struct in structs:
        forward_decls.append(f"void init_{struct.f_name}(py::module &, py::class_<{struct.cpp_class}> &);")

    # for struct in structs:
    #     forward_decls.append(f"std::string to_string(const {struct.cpp_class}& self);")
    forward_decls.append("} // namespace Pybmad")
    tpl = string.Template(template.replace("// ${", "${"))
    return tpl.substitute(forward_declarations="\n".join(forward_decls))


def generate_pybmad_struct_code(
    structs: list[CodegenStructure],  # noqa: ARG001
    routines_by_name: dict[str, FortranRoutine],  # noqa: ARG001
    struct: CodegenStructure,
    used_array_dims: set[int],
) -> list[str]:
    code_lines = [""]
    code_lines.append("// =============================================================================")
    code_lines.append(f"// {struct.f_name}")
    code_lines.append(f"void init_{struct.f_name}(py::module &m, py::class_<{struct.cpp_class}> &cls) {{")

    # if struct.c_constructor_arg_list:
    #     code_lines.append("        // TODO: add proper constructor with arguments")
    #     code_lines.append("        .def(py::init<>())")
    # else:
    code_lines.append("        cls.def(py::init<>())")

    for arg in struct.arg:
        if not arg.is_component:
            continue

        try:
            tpl = proxy_templates[arg.full_type]
        except KeyError:
            code_lines.append(f"        // {arg.full_type} {arg.c_name} proxy support missing")
            continue

        comment = arg.comment.replace('"', "'") if arg.comment else ""
        code_lines.append(f"        // {struct.cpp_class}.{arg.c_name} ({arg.full_type} - {comment}")
        if tpl.fortran_setter:
            code_lines.append(
                f'        .def_property("{arg.python_name}", &{struct.cpp_class}::{arg.c_name}, &{struct.cpp_class}::set_{arg.c_name})'
            )
        else:
            code_lines.append(
                f'        .def_property_readonly("{arg.python_name}", &{struct.cpp_class}::{arg.c_name})'
            )

    if 1 in used_array_dims:
        container_cls = f"{struct.cpp_class}Alloc1D"
        code_lines.append(
            f'      .def_static("new_array1d", [](int sz, int lbound) {{ return {container_cls}(lbound, sz); }}, '
            f'py::arg("sz"), py::arg("lbound") = 1)'
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
    code_lines.append("        ;")
    code_lines.append("")

    for n in SUPPORTED_ARRAY_DIMS:
        if n in used_array_dims:
            code_lines.append(
                f'    bind_FTypeArrayND<{struct.cpp_class}Array{n}D>(m, "{struct.python_class_name}Array{n}D");'
            )
            if n == 1:
                code_lines.append(
                    f'    bind_FTypeAlloc1D<{struct.cpp_class}Alloc1D>(m, "{struct.python_class_name}Alloc1D");'
                )
        else:
            code_lines.append(f"    // {n}D {struct.cpp_class} arrays are not used in structs/routines")

    code_lines.append("}")

    return code_lines


def generate_pybmad(
    structs: list[CodegenStructure],
    routines_by_name: dict[str, FortranRoutine],
    enums: dict[str, dict[str, EnumValue]] | None = None,
) -> dict[pathlib.Path, str]:
    """
    Generate a split pybind11 binding:
      - _core.cpp: defines PYBIND11_MODULE and calls per-structure init functions
      - One .cpp per structure: defines void init_<struct_name>(py::module &)
      - A header declaring init functions (included by _core.cpp and each per-struct file)
    """

    if enums is None:
        enums = parse_all_enums()
    files: dict[pathlib.Path, str] = {}

    code_lines = []
    code_lines.append("#include <pybind11/pybind11.h>")
    code_lines.append('#include "bmad/generated/proxy.hpp"')
    code_lines.append('#include "bmad/generated/to_string.hpp"')
    code_lines.append('#include "bmad/to_string.hpp"')
    code_lines.append('#include "pybmad/arrays.hpp"')
    code_lines.append("")
    code_lines.append("using namespace Pybmad;")
    code_lines.append("namespace py = pybind11;")
    code_lines.append("")
    code_lines.append("namespace Pybmad {")

    def by_module_and_name(struct: CodegenStructure):
        return (struct.module, struct.f_name)

    array_usage = struct_array_usage_dimensions(routines_by_name, structs)
    for struct in sorted(structs, key=by_module_and_name):
        code_lines.extend(
            generate_pybmad_struct_code(
                structs, routines_by_name, struct, array_usage.get(struct.f_name.lower(), set())
            )
        )

    code_lines.append("} // namespace Pybmad")
    filename = PYBMAD_SRC / "generated" / "structs.cpp"
    files[filename] = "\n".join(code_lines)

    enum_code = generate_enum_wrapper_code(enums)
    hand_written_bindings = "\n".join(
        (
            "bind_standard_arrays(m);",
            "init_common_structs(m);",
        )
    )

    template_text = (CODEGEN_ROOT / "pybind.tpl.cpp").read_text()
    # Ensure ${ substitution does not break comment placeholders
    template_text = template_text.replace("// ${", "${")
    template = string.Template(template_text)

    struct_init_lines = []
    for struct in structs:
        struct_init_lines.append(
            f'    auto py_{struct.python_class_name} = py::class_<{struct.cpp_class}>(m, "{struct.python_class_name}", "Fortran struct: {struct.f_name}");'
        )

    for struct in structs:
        struct_init_lines.append(f"    init_{struct.f_name}(m, py_{struct.python_class_name});")

    init_calls = "\n".join(struct_init_lines)

    routine_block = generate_py_routine_defs(routines_by_name)
    routine_wrappers = generate_py_routine_wrappers(routines_by_name)
    enum_code = enum_code.replace("\\n", "\\n    ")
    core_body = textwrap.dedent(
        f"""
        m.doc() = "pybmad";

        // Per-structure bindings
{init_calls}

        // Hand-written bindings
{hand_written_bindings}

        // Enums
{enum_code or "// (No enums)"}

        // Routines
{routine_block}
        """
    ).rstrip()

    substituted = template.substitute(
        pybind11_routine_wrappers=routine_wrappers,
        pybind11_definitions=core_body,
    )

    files[PYBMAD_SRC / "generated" / "_pybmad.cpp"] = substituted

    pybmad_header = generate_pybmad_header(
        template=(CODEGEN_ROOT / "pybind.tpl.hpp").read_text(),
        routines=routines_by_name,
        structs=structs,
    )
    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "init.hpp"] = pybmad_header
    return files
