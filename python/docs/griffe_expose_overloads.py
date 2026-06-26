"""Griffe extension: expose overload-only functions from stub files.

Griffe stashes ``@overload`` defs in ``module.overloads`` and only promotes them
to real members when a non-overload implementation follows. Pure ``.pyi`` stubs
(e.g. nanobind output) have no implementation, so overload-only functions are
never exposed and mkdocstrings fails with "Could not collect". This extension
synthesizes a member for each such function, using the last overload as the
canonical definition and attaching all overloads.

It also pre-computes, for the docs template, which docstring sections are shared
across every overload versus unique to a single overload, so the template can
render the common prose/returns once and only the differing parameters per
overload (see ``overrides/python/readthedocs/function.html.jinja``).
"""

from __future__ import annotations

import griffe

# Namespace used for data stashed in ``Object.extra`` for the template.
EXTRA_NS = "expose_overloads"

# Common sections of these kinds render *before* the per-overload blocks
# (descriptive prose); any other common section renders *after* (returns, etc.).
_TOP_KINDS = frozenset({"text", "admonition", "examples", "deprecated"})


def _element_key(element: object) -> tuple | str:
    """A value-based identity for a single docstring element."""
    if isinstance(element, str):
        return element
    parts = []
    for attr in ("name", "title", "annotation", "description", "value", "contents", "kind"):
        if hasattr(element, attr):
            parts.append((attr, str(getattr(element, attr))))
    return tuple(parts) if parts else str(element)


def _section_key(section: griffe.DocstringSection) -> tuple:
    """A hashable, value-based identity for a parsed docstring section."""
    kind = section.kind.value
    value = section.value
    if isinstance(value, str):
        return (kind, value)
    if isinstance(value, (list, tuple)):
        return (kind, tuple(_element_key(item) for item in value))
    return (kind, _element_key(value))


def _split_sections(overloads: list[griffe.Function]) -> dict | None:
    """Partition overload docstring sections into common vs per-overload-unique.

    Returns None when there is nothing useful to dedupe (fewer than two
    overloads, or none carry a docstring).
    """
    if len(overloads) < 2:
        return None

    parsed = [o.docstring.parsed if o.docstring else [] for o in overloads]
    if not any(parsed):
        return None

    # A section is common when its value-key appears in every overload.
    key_sets = [{_section_key(s) for s in sections} for sections in parsed]
    common_keys = set.intersection(*key_sets) if key_sets else set()

    # Preserve document order using the first overload that has the section.
    common_seen: set = set()
    common_top: list = []
    common_bottom: list = []
    for sections in parsed:
        for s in sections:
            key = _section_key(s)
            if key in common_keys and key not in common_seen:
                common_seen.add(key)
                (common_top if s.kind.value in _TOP_KINDS else common_bottom).append(s)

    per_overload = []
    for overload, sections in zip(overloads, parsed, strict=True):
        unique = [s for s in sections if _section_key(s) not in common_keys]
        per_overload.append({"func": overload, "unique": unique})

    return {
        "enabled": True,
        "common_top": common_top,
        "common_bottom": common_bottom,
        "per_overload": per_overload,
    }


class ExposeOverloads(griffe.Extension):
    def on_module_members(self, *, mod: griffe.Module, **kwargs) -> None:  # noqa: ARG002
        for name, overloads in mod.overloads.items():
            if not overloads or name in mod.members:
                continue
            canonical = overloads[-1]
            canonical.overloads = list(overloads)
            mod[name] = canonical

            split = _split_sections(list(overloads))
            if split is not None:
                canonical.extra[EXTRA_NS] = split
