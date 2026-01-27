from __future__ import annotations

import pathlib
import string
import textwrap

from .arg import CodegenStructure
from .context import SUPPORTED_ARRAY_DIMS
from .cpp import CppWrapperArgument
from .enums import EnumValue, get_ele_attributes, get_ele_keys
from .paths import CODEGEN_ROOT, PYBMAD_INCLUDE, PYBMAD_LIB, PYBMAD_SRC
from .proxy import _generate_proxy_constructor_arg
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
    """
    Generate a Python enum.IntEnum class definition string.
    """
    code = ["", "", f"class {clsname}(enum.IntEnum):"]

    if not enums:
        code.append("    pass")
        return "\n".join(code)

    for attr in enums:
        comment = f"  # {attr.comment}" if attr.comment else ""
        code.append(f"    {attr.name} = {attr.value}{comment}")

    return "\n".join(code)


def generate_enum_constants(all_enums: dict[str, dict[str, EnumValue]]) -> str:
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

            result.append(f"{enum.name} = {enum.value}")

    return "\n".join(result)


def generate_enum_wrapper_code(enums: dict[str, dict[str, EnumValue]]) -> str:
    return "\n".join(
        (
            "import enum",
            "",
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

    clsname, full_clsname = routine.python_class_return_type

    lines = []
    lines.append(
        f'    py::class_<{full_clsname}, std::unique_ptr<{full_clsname}>>(m, "{clsname}", "{routine.name} return type")'
    )
    for arg in outputs:
        lines.append(f'        .def_readonly("{arg.python_name}", &{full_clsname}::{arg.c_name})')

    lines.append(f'        .def("__len__", [](const {full_clsname} &) {{ return {len(outputs)}; }})')
    lines.append(f'        .def("__getitem__", [](const {full_clsname} &s, int i) -> py::object {{')
    lines.append(f"            if (i < 0) i += {len(outputs)};")
    for i, arg in enumerate(outputs):
        lines.append(f"            if (i == {i}) return py::cast(s.{arg.c_name});")
    lines.append("            throw py::index_error();")
    lines.append("        })")

    # lines.append(f'      .def("__repr__", [](const {full_clsname} &self){{ return to_string(self); }})')
    lines.append("        ;")
    return lines


def generate_py_routine_return_value_struct(routine: FortranRoutine) -> list[str]:
    args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]
    outputs = [arg for arg in args if arg.arg.intent == "out"]

    name, py_name = routine.python_class_return_type
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

    doc = routine.docstring.to_numpy_docstring(routine.arg_names_with_result)
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


def generate_pybmad_struct_code(struct: CodegenStructure, used_array_dims: set[int]) -> list[str]:
    code_lines = [""]
    code_lines.append("// =============================================================================")
    code_lines.append(f"// {struct.f_name}")
    code_lines.append(f"void init_{struct.f_name}(py::module &m, py::class_<{struct.cpp_class}> &cls) {{")

    ctor_args: list[str] = []
    ctor_types: list[str] = []

    for arg in struct.arg:
        (ctor_type, _ctor_body) = _generate_proxy_constructor_arg(struct, arg)
        if ctor_type is not None:
            ctor_args.append(f'py::arg("{arg.python_name}") = py::none()')
            ctor_types.append(ctor_type)

    if ctor_args:
        types_str = ", ".join(ctor_types)
        args_str = ", ".join(ctor_args)

        # cls.def(py::init<T1, T2>(), py::arg("x")=none, py::arg("y")=none);
        code_lines.append(f"    cls.def(py::init<{types_str}>(),")
        code_lines.append(f"        {args_str})")
    else:
        code_lines.append("    cls.def(py::init<>())")

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
    code_lines.append(
        textwrap.dedent(f"""
            .def("__copy__", [](const {struct.cpp_class} &self){{
                return {struct.cpp_class}(self);  // under-the-hood fortran copy
            }})
            .def("__deepcopy__", [](const {struct.cpp_class} &self, py::dict& memo){{
                return {struct.cpp_class}(self);
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


def _group_structures_by_char(
    structs: list[CodegenStructure],
) -> dict[str, list[CodegenStructure]]:
    """Group structures by the first letter of their Fortran name."""
    structs_by_char = {}
    for st in structs:
        ch = st.f_name[0].lower()
        structs_by_char.setdefault(ch, []).append(st)
    return structs_by_char


def _group_routines_by_source_and_char(
    routines_by_name: dict[str, FortranRoutine],
) -> dict[str, dict[str, list[FortranRoutine]]]:
    """Group routines by source namespace and first letter."""
    routines_map = {}

    for routine in routines_by_name.values():
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
            header_decls.append(f"void init_{st.f_name}(py::module &m, py::class_<{st.cpp_class}> &class_);")

        newline = "\n"
        files[header_path] = textwrap.dedent(f"""\
            #pragma once
            #include <pybind11/pybind11.h>
            #include "bmad/generated/proxy.hpp"
            #include "pybmad/generated/structs.hpp"
            namespace py = pybind11;

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
            "",
            "using namespace Pybmad;",
            "namespace py = pybind11;",
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
) -> tuple[list[str], list[str]]:
    """
    Generate the split C++ files ({source}_routines_{letter}) for routines.
    Returns the list of generated headers and initialization function calls.
    """
    headers: list[str] = ["#pragma once"]
    init_calls: list[str] = []

    for src, chars_dict in sorted(routines_map.items()):
        for char, routine_list in sorted(chars_dict.items()):
            # Reconstruct the dictionary using r.name as the key
            subset_routines = {r.name: r for r in routine_list}

            base_name = f"{src}_routines_{char}"
            header_name = f"{base_name}.hpp"
            src_name = f"{base_name}.cpp"
            init_fn_name = f"init_{src}_routines_{char}"

            wrapper_structs, wrappers_code = generate_py_routine_wrappers(subset_routines)

            header_path = PYBMAD_INCLUDE / "pybmad" / "generated" / header_name
            files[header_path] = textwrap.dedent(f"""\
                #pragma once
                #include <pybind11/complex.h>
                #include <pybind11/numpy.h>
                #include <pybind11/pybind11.h>
                #include <pybind11/stl.h>

                #include "pybmad/arrays.hpp"
                #include "pybmad/util.hpp"

                namespace py = pybind11;

                void {init_fn_name}(py::module &m);

                {wrapper_structs}
            """)

            headers.append(f'#include "pybmad/generated/{header_name}"')
            init_calls.append(f"{init_fn_name}(m);")

            defs_block = generate_py_routine_defs(subset_routines)
            # defs_block_indented = textwrap.indent(defs_block, "  ")

            cpp_content = textwrap.dedent(f"""\
                #include "pybmad/generated/{header_name}"

                namespace py = pybind11;
                using namespace pybind11::literals;
                using namespace Pybmad;

                {wrappers_code}
                
                void {init_fn_name}(py::module &m) {{
                {defs_block}
                }}
            """)

            files[PYBMAD_SRC / "generated" / src_name] = cpp_content

    return headers, init_calls


def _generate_struct_init(structs: list[CodegenStructure]):
    src_lines = []
    for st in structs:
        src_lines.append(
            f"    auto py_{st.python_class_name} = py::class_<{st.cpp_class}>(m, "
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


def generate_pybmad(
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

    structs_by_char = _group_structures_by_char(structs)
    routines_map = _group_routines_by_source_and_char(routines_by_name)

    struct_headers = _generate_structure_files(files, structs_by_char, array_usage)
    struct_inits = _generate_struct_init(structs)

    routine_headers, routine_inits = _generate_routine_files(files, routines_map)

    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "structs.hpp"] = "\n".join(struct_headers)
    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "routines.hpp"] = "\n".join(routine_headers)

    _generate_main_module_file(files, enums, struct_inits, routine_inits)

    files[PYBMAD_INCLUDE / "pybmad" / "generated" / "init.hpp"] = generate_pybmad_header(
        template=(CODEGEN_ROOT / "pybind.tpl.hpp").read_text(),
    )

    return files
