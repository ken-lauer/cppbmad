"""
Generate per-project API documentation pages for mkdocs.
"""

from __future__ import annotations

import logging
import pathlib

from . import paths
from .arg import Argument, CodegenStructure
from .context import CodegenConfig, ProjectSettings, UpstreamInfo
from .docstring import type_information_to_python_type
from .enums import EnumValue
from .routines import FortranRoutine
from .util import struct_to_proxy_class_name

logger = logging.getLogger(__name__)

DOCS_SOURCE = pathlib.Path(__file__).resolve().parent.parent / "python" / "docs" / "source"


def _source_link(
    filepath: pathlib.Path,
    lineno: int,
    source_dir: pathlib.Path,
    upstream: UpstreamInfo,
    upstream_url: str,
) -> str:
    """Markdown link to the source file on GitHub."""
    expanded = paths.normalize(filepath)
    try:
        rel = expanded.relative_to(source_dir.resolve())
    except ValueError:
        return f"`{filepath.name}:{lineno}`"

    if not upstream.commit_hash or not upstream_url:
        return f"`{rel}:{lineno}`"

    url = f"{upstream_url.rstrip('/')}/blob/{upstream.commit_hash}/{rel}#L{lineno}"
    return f"[`{rel}`]({url})"


def _python_type_for_arg(arg: Argument, struct_links: dict[str, str] | None = None) -> str:
    """Python type string for a struct attribute."""

    python_type = type_information_to_python_type(arg.member.type_info)
    if arg.type == "type":
        cls = struct_to_proxy_class_name(arg.kind)
        if struct_links and cls in struct_links:
            return f"[{python_type}]({struct_links[cls]})"
    return python_type


def _generate_struct_section(
    struct: CodegenStructure,
    struct_links: dict[str, str],
    source_dir: pathlib.Path,
    upstream: UpstreamInfo,
    upstream_url: str,
) -> list[str]:
    """Generate markdown for a single struct with an attribute table."""
    assert struct.parsed is not None
    src = _source_link(struct.parsed.filename, struct.parsed.line, source_dir, upstream, upstream_url)
    lines = [
        # Register identifier for autorefs cross-references
        f"::: pybmad.{struct.python_class_name}",
        "    options:",
        "      heading_level: 0",
        "      show_root_heading: false",
        "      members: false",
        "      show_signature: false",
        "      show_bases: false",
        "      show_docstring_description: false",
        "",
        f"### {struct.python_class_name}",
        "",
        f"Fortran struct: `{struct.f_name}` ({src})",
        "",
        "All attributes may be passed to the initializer as arguments:",
        "",
    ]

    attrs = [a for a in struct.arg if a.is_component]
    if attrs:
        lines.append("| Attribute | Type | Description |")
        lines.append("|-----------|------|-------------|")
        for attr in attrs:
            py_type = _python_type_for_arg(attr, struct_links)
            comment = attr.comment.replace("|", "\\|") if attr.comment else ""
            lines.append(f"| `{attr.python_name}` | {py_type} | {comment} |")
        lines.append("")

    return lines


def _generate_project_page(
    project: str,
    structs: list[CodegenStructure],
    routines: list[FortranRoutine],
    project_info: dict[str, ProjectSettings],
    struct_links: dict[str, str],
    source_dir: pathlib.Path,
    upstream: UpstreamInfo,
    upstream_url: str,
) -> str:
    """Generate a markdown page for a single project."""
    info = project_info.get(project)
    title = info.title if info else project.title()
    description = info.description if info else ""
    lines = [
        f"# {title}",
        "",
        description,
        "",
    ]

    if structs:
        lines.append("## Classes (Fortran Structures)")
        lines.append("")
        for struct in sorted(structs, key=lambda s: s.python_class_name):
            lines.extend(_generate_struct_section(struct, struct_links, source_dir, upstream, upstream_url))

    if routines:
        lines.append("## Procedures")
        lines.append("")
        # Group by overloaded_name to avoid duplicate mkdocstrings directives
        grouped: dict[str, list[FortranRoutine]] = {}
        for routine in routines:
            grouped.setdefault(routine.overloaded_name, []).append(routine)

        for name in sorted(grouped):
            variants = grouped[name]
            src_links = [
                _source_link(r.start_line.filename, r.start_line.lineno, source_dir, upstream, upstream_url)
                for r in variants
            ]
            lines.append(f"### {name}")
            lines.append("")
            if len(variants) == 1:
                lines.append(f"Fortran source: {src_links[0]}")
            else:
                lines.append("Fortran sources (overloaded):")
                lines.append("")
                for r, src in zip(variants, src_links, strict=True):
                    lines.append(f"- `{r.name}`: {src}")
            lines.append("")
            lines.append(f"::: pybmad.{name}")
            lines.append("    options:")
            lines.append("      show_root_heading: false")
            lines.append("      show_root_toc_entry: false")
            lines.append("")

    return "\n".join(lines)


def _generate_index_page(
    project_structs: dict[str, list[CodegenStructure]],
    project_routines: dict[str, list[FortranRoutine]],
    enums: dict[str, dict[str, EnumValue]],
    project_info: dict[str, ProjectSettings],
) -> str:
    """Generate an alphabetical index page linking to per-project pages."""
    lines = [
        "# API Index",
        "",
        "| Name | Type | Project |",
        "|------|------|---------|",
    ]

    entries: list[tuple[str, str, str, str]] = []  # (sort_key, name, type, project)

    for project, structs in project_structs.items():
        info = project_info.get(project)
        title = info.title if info else project.title()
        for struct in structs:
            name = struct.python_class_name
            anchor = name.lower()
            entries.append(
                (
                    name.lower(),
                    f"[{name}]({project}.md#{anchor})",
                    "Struct",
                    f"[{title}]({project}.md)",
                )
            )

    for project, routines in project_routines.items():
        info = project_info.get(project)
        title = info.title if info else project.title()
        seen = set()
        for routine in routines:
            name = routine.overloaded_name

            if name in seen:
                continue

            seen.add(name)
            anchor = name.lower()
            entries.append(
                (
                    name.lower(),
                    f"[`{name}`]({project}.md#{anchor})",
                    "Routine",
                    f"[{title}]({project}.md)",
                )
            )

    for per_file_enums in enums.values():
        for enum in per_file_enums.values():
            entries.append(
                (
                    enum.name.lower(),
                    f"`{enum.name}`",
                    "Enum",
                    "[Enums](enums.md)",
                )
            )

    entries.sort(key=lambda e: e[0])

    for _, name, type_, project in entries:
        lines.append(f"| {name} | {type_} | {project} |")

    lines.append("")
    return "\n".join(lines)


def _generate_enums_page(enums: dict[str, dict[str, EnumValue]]) -> str:
    """Generate a page listing all enums/constants."""
    lines = [
        "# Enums / Constants",
        "",
        "Fortran integer and real parameters exposed as Python constants.",
        "",
    ]

    for source_file, per_file_enums in sorted(enums.items()):
        if not per_file_enums:
            continue
        lines.append(f"## From `{source_file}`")
        lines.append("")
        lines.append("| Name | Value | Description |")
        lines.append("|------|-------|-------------|")
        for enum in per_file_enums.values():
            comment = enum.comment or ""
            lines.append(f"| `{enum.name}` | `{enum.value}` | {comment} |")
        lines.append("")

    return "\n".join(lines)


def generate_docs(
    structs: list[CodegenStructure],
    routines_by_name: dict[str, FortranRoutine],
    enums: dict[str, dict[str, EnumValue]],
    config: CodegenConfig,
    upstream: UpstreamInfo | None = None,
) -> dict[pathlib.Path, str]:
    """
    Generate per-project API documentation markdown files.

    Returns
    -------
    dict[pathlib.Path, str]
        Mapping of output file paths to their generated content.
    """
    api_dir = DOCS_SOURCE / "api"
    namespace_to_project = {p.cpp_namespace: p.name for p in config.projects if p.cpp_namespace}
    project_info = {p.name: p for p in config.projects}
    known_projects = set(project_info.keys())
    if upstream is None:
        upstream = UpstreamInfo()
    source_dir = config.source_dir
    upstream_url = config.upstream_source_url

    project_structs: dict[str, list[CodegenStructure]] = {}
    for struct in structs:
        project = struct.project
        assert project is not None
        if project not in known_projects:
            continue
        project_structs.setdefault(project, []).append(struct)

    project_routines: dict[str, list[FortranRoutine]] = {}
    for routine in routines_by_name.values():
        if not routine.usable:
            continue
        if routine.cpp_namespace == "CppBmadTest":
            continue
        project = namespace_to_project.get(routine.cpp_namespace, "bmad")
        project_routines.setdefault(project, []).append(routine)

    # Build cross-reference map: ClassName -> "project.md#classname"
    struct_links: dict[str, str] = {}
    for project, pstructs in project_structs.items():
        for s in pstructs:
            anchor = s.python_class_name.lower()
            struct_links[s.python_class_name] = f"{project}.md#{anchor}"

    files: dict[pathlib.Path, str] = {}

    all_projects = sorted(set(list(project_structs.keys()) + list(project_routines.keys())))
    for project in all_projects:
        page = _generate_project_page(
            project,
            project_structs.get(project, []),
            project_routines.get(project, []),
            project_info,
            struct_links,
            source_dir,
            upstream,
            upstream_url,
        )
        files[api_dir / f"{project}.md"] = page
        logger.info(
            "Docs: %s - %d structs, %d routines",
            project,
            len(project_structs.get(project, [])),
            len(project_routines.get(project, [])),
        )

    files[api_dir / "index.md"] = _generate_index_page(project_structs, project_routines, enums, project_info)
    files[api_dir / "enums.md"] = _generate_enums_page(enums)
    return files
