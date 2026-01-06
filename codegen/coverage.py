from __future__ import annotations

import collections
import html

from .arg import CodegenStructure
from .routines import FortranRoutine


def _tag(idx: str, content: str, cls: str = "", style: str = "", extra_attrs: str = "") -> str:
    """Helper to generate HTML tags with less verbosity."""
    attr_cls = f' class="{cls}"' if cls else ""
    attr_sty = f' style="{style}"' if style else ""
    return f"<{idx}{attr_cls}{attr_sty}{extra_attrs}>{content}</{idx}>"


def _classify_failure_reason(reason: str) -> str:
    """
    Groups specific detailed error messages into broader categories
    for the summary statistics table.
    """
    lower_r = reason.lower()
    if "argument not defined" in lower_r:
        if "have:" in lower_r:
            return "Arg string defined in docstring/interface but missing in struct"
        return "Argument definition missing"
    if "untranslated type" in lower_r:
        if lower_r.split()[-2] == "c_ptr":
            return "C pointer usage"
        return "Untranslated Fortran structure"
    # if "variable" in lower_r and "array" in lower_r:
    #     return "Unsupported Variable/Pointer Array dimension"
    if "docstring" in lower_r:
        return "Missing Docstring"
    if "configuration skip list" in lower_r:
        return "Skipped by Config"
    if "test_build" in lower_r:
        return "Test Build Skip"
    if "module name unset" in lower_r:
        return "Routine not in a module"

    return reason


CSS = """
    :root { 
        --p-color: #007bff; 
        --ok-color: #28a745; 
        --warn-color: #ffc107; 
        --err-color: #ab2e3e; 
        --bg-gray: #f8f9fa; 
        --border-color: #dee2e6;
    }
    body { 
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        padding: 20px; color: #333; max-width: 1400px; margin: 0 auto; line-height: 1.5; background-color: #fff;
    }

    /* Dashboard Grid */
    .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 40px; }
    .card { background: white; padding: 20px; border: 1px solid var(--border-color); border-radius: 8px; text-align: center; }
    .card h3 { margin: 0 0 10px 0; font-size: 13px; text-transform: uppercase; color: #6c757d; font-weight: 700; letter-spacing: 0.5px; } 
    .card .val { font-size: 28px; font-weight: 700; color: #212529; }
    .card small { display: block; margin-top: 5px; color: #6c757d; font-size: 0.9em; }

    /* Layout & Tables */
    h2 { border-bottom: 2px solid var(--bg-gray); padding-bottom: 10px; margin-top: 50px; scroll-margin-top: 60px; color: #495057; }
    p.intro { color: #6c757d; margin-bottom: 20px; }

    .table-container { overflow-x: auto; border: 1px solid var(--border-color); border-radius: 6px; box-shadow: 0 1px 3px rgba(0,0,0,0.02); }
    table { border-collapse: collapse; width: 100%; font-size: 13px; }
    th { background-color: var(--bg-gray); position: sticky; top: 0; font-weight: 600; text-align: left; padding: 12px; border-bottom: 2px solid var(--border-color); color: #495057; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: top; }
    tbody tr:last-child td { border-bottom: none; }
    tbody tr:hover { background-color: #f8f9fa; }
    
    /* Specific Columns */
    .col-mono { font-family: SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; color: #E83E8C; font-size: 12px; }
    .col-loc { font-family: monospace; color: #666; font-size: 11px; white-space: nowrap; }
    .col-stat { width: 100px; }

    /* Badges */
    .badge { display: inline-flex; align-items: center; justify-content: center; padding: 3px 8px; border-radius: 4px; font-size: 11px; font-weight: 700; min-width: 50px; }
    .badge-ok { background-color: #d4edda; color: #155724; }
    .badge-warn { background-color: #fff3cd; color: #856404; }
    .badge-err { background-color: #f8d7da; color: #721c24; }
    .badge-subtle { background-color: #e2e3e5; color: #383d41; }
    .badge-wrap { background-color: #d1ecf1; color: #0c5460; margin-left: 5px; }

    /* Lists */
    .reason-list { margin: 0; padding-left: 15px; color: #666; font-size: 0.95em; }
    .reason-list li { margin-bottom: 2px; }
    .f-sub { font-size: 0.85em; color: #999; margin-left: 5px; }
    
    /* Row Grouping for Duplicates */
    .group-start td { border-top: 2px solid #ced4da; }

    /* Search Input */
    .table-search {
        width: 100%;
        max-width: 300px;
        padding: 8px 12px;
        margin-bottom: 10px;
        border: 1px solid var(--border-color);
        border-radius: 4px;
        font-size: 13px;
    }
    .table-search:focus { outline: none; border-color: var(--p-color); box-shadow: 0 0 0 2px rgba(0,123,255,.1); }

    /* Table of Contents */
    .toc-container { background: var(--bg-gray); padding: 15px 20px; border-radius: 8px; margin-bottom: 30px; border: 1px solid var(--border-color); }
    .toc-title { font-weight: 700; margin-bottom: 10px; font-size: 0.9em; text-transform: uppercase; color: #6c757d; }
    .toc-list { list-style: none; padding: 0; margin: 0; display: flex; flex-wrap: wrap; gap: 15px; }
    .toc-list li a { text-decoration: none; color: var(--p-color); font-weight: 500; font-size: 14px; }
    .toc-list li a:hover { text-decoration: underline; }
"""

SCRIPT_JS = """
<script>
document.addEventListener('DOMContentLoaded', function() {
    const searchInputs = document.querySelectorAll('.table-search');
    
    searchInputs.forEach(input => {
        input.addEventListener('keyup', function() {
            const getCellValue = (tr, idx) => tr.children[idx].innerText || tr.children[idx].textContent;
            
            const tableId = this.dataset.target;
            const table = document.getElementById(tableId);
            const body = table.querySelector('tbody');
            const trs = body.querySelectorAll('tr');
            const filter = this.value.toUpperCase();
            
            trs.forEach(tr => {
                // Check both displayed text and the data-search attribute (for hidden terms)
                const visibleText = tr.textContent || tr.innerText;
                const hiddenText = tr.getAttribute('data-search') || '';
                
                // Allow row valid if either matches
                const content = (visibleText + " " + hiddenText).toUpperCase();
                
                if (content.indexOf(filter) > -1) {
                    tr.style.display = "";
                } else {
                    tr.style.display = "none";
                }
            });
        });
    });
});
</script>
"""


def generate_coverage_report(
    routines: list[FortranRoutine],
    structs: list[CodegenStructure],
    missing_struct_attrs: dict[str, dict[str, str]],
) -> str:
    """Generates a standalone HTML report."""

    routines_by_name: dict[str, list[FortranRoutine]] = collections.defaultdict(list)
    for r in routines:
        if "Skipped by Config" in [_classify_failure_reason(reason) for reason in r.unusable_reason]:
            continue

        routines_by_name[r.name].append(r)

    # Unique Coverage
    total_unique_routines = len(routines_by_name)
    usable_unique_routines = 0
    totally_skipped_routines = 0

    # Track stats for things that are completely broken (no usable version found)
    failure_counts = collections.Counter()
    missing_type_counts = collections.Counter()

    for variants in routines_by_name.values():
        is_covered = any(r.usable for r in variants)

        if is_covered:
            usable_unique_routines += 1
        else:
            totally_skipped_routines += 1
            seen_reasons = set()
            for r in variants:
                for reason in r.unusable_reason:
                    if reason.startswith("Untranslated type:"):
                        _, type_and_array = reason.split(":", 1)
                        type_ = type_and_array.strip().split()[0]
                        missing_type_counts[type_.strip()] += 1

                    cat = _classify_failure_reason(reason)
                    if cat not in seen_reasons:
                        failure_counts[cat] += 1
                        seen_reasons.add(cat)

    routine_health = (usable_unique_routines / total_unique_routines * 100) if total_unique_routines else 0

    # Struct Stats
    total_structs = len(structs)
    structs_partial = len(missing_struct_attrs)
    structs_full = total_structs - structs_partial
    struct_health = (structs_full / total_structs * 100) if total_structs else 0

    # --- Body ---

    # -- TOC & Header --
    toc_html = """
    <div class="toc-container">
        <div class="toc-title">Table of Contents</div>
        <ul class="toc-list">
            <li><a href="#sec-structs">Structure Definitions</a></li>
            <li><a href="#sec-types">Untranslated Types</a></li>
            <li><a href="#sec-blockers">Implementation Blockers</a></li>
            <li><a href="#sec-routines">Routine Definitions</a></li>
        </ul>
    </div>
    """

    summary_html = f"""
    <div class="summary-grid">
        <div class="card">
            <h3>Structures</h3>
            <div class="val" style="color:{"#28a745" if struct_health > 90 else "#dc3545"}">{struct_health:.1f}%</div>
            <small>Full<br>{structs_full}/{total_structs}</small>
        </div>
        <div class="card">
            <h3>Routines</h3>
            <div class="val" style="color:{"#28a745" if routine_health > 50 else "#dc3545"}">{routine_health:.1f}%</div>
            <small>Usable<br>{usable_unique_routines}/{total_unique_routines}</small>
        </div>
    </div>
    """

    # -- Structs Table --
    rows_structs = []
    for s in sorted(structs, key=lambda x: x.f_name):
        removed_map = missing_struct_attrs.get(s.f_name, {})
        loc_str = html.escape(s.module or "N/A")

        if not removed_map:
            status_html = _tag("span", "FULL", "badge badge-ok")
            dropped_html = ""
        else:
            status_html = _tag("span", "PARTIAL", "badge badge-warn")
            items = []
            for f_name, raw_type in removed_map.items():
                items.append(
                    f"<li><strong style='color:#333'>{html.escape(f_name)}</strong> <span class='f-sub'>({html.escape(raw_type)})</span></li>"
                )
            dropped_html = _tag("ul", "".join(items), "reason-list")

        # Create search text for data-search to allow easy filtering
        search_txt = f"{s.f_name} {loc_str}".lower()

        row = [
            _tag("td", html.escape(s.f_name), "col-mono"),
            _tag("td", loc_str, "col-loc"),
            _tag("td", str(len(s.arg))),
            _tag("td", status_html, "col-stat"),
            _tag("td", dropped_html),
        ]
        rows_structs.append(_tag("tr", "".join(row), extra_attrs=f' data-search="{html.escape(search_txt)}"'))

    # -- Routines Table --
    rows_routines = []

    # Sort by Name first, then filename to keep duplicates adjacent
    sorted_names = sorted(routines_by_name.keys())

    for name in sorted_names:
        variants = routines_by_name[name]

        # Sort variants so "Usable" ones usually appear on top if mixed
        variants_sorted = sorted(variants, key=lambda x: (not x.usable, str(x.module)))

        # Determine if the group (Name) is effectively covered
        group_covered = any(r.usable for r in variants)

        for i, r in enumerate(variants_sorted):
            # Styling for start of a new group
            tr_class = "group-start" if i == 0 else ""

            # Location Render
            fname = r.filename.name if r.filename else ""
            mod = r.module if r.module else ""
            if fname and mod:
                loc_txt = f"{mod}<br><span style='color:#999'>{fname}</span>"
            else:
                loc_txt = mod or fname or "Unknown"

            # Status Render
            if r.usable:
                status_html = _tag("span", "USABLE", "badge badge-ok")
                reason_html = ""
                if r.needs_python_wrapper:
                    status_html += _tag("span", "WRAP", "badge badge-wrap")
            else:
                # If the routine is skipped, but a sibling variant is valid,
                # mark it as "Redundant Skip" (grey) instead of "Error Skip" (Red)
                # to help the user focus on actual problems.
                if group_covered:
                    status_html = _tag("span", "SKIP", "badge badge-subtle")
                    reason_style = "color:#999"
                else:
                    status_html = _tag("span", "SKIP", "badge badge-err")
                    reason_style = ""

                items = [f"<li>{html.escape(reas)}</li>" for reas in r.unusable_reason]
                reason_html = _tag("ul", "".join(items), "reason-list", reason_style)

            # NOTE for Filtering:
            # We visually hide the name in subsequent rows (using empty string).
            # To make filtering work (e.g. searching the name finds all rows),
            # we embed the name in the data-search attribute.
            search_terms = f"{r.name} {mod} {fname} {loc_txt}"

            row = [
                # Only show Name in the first row of the group
                _tag("td", html.escape(r.name) if i == 0 else "", "col-mono"),
                _tag("td", loc_txt, "col-loc"),
                _tag("td", status_html, "col-stat"),
                _tag("td", reason_html),
            ]
            rows_routines.append(
                _tag("tr", "".join(row), tr_class, extra_attrs=f' data-search="{html.escape(search_terms)}"')
            )

    # -- MISSING TYPES TABLE --
    rows_types = []
    if missing_type_counts:
        for m_type, count in missing_type_counts.most_common():
            row = [_tag("td", html.escape(m_type), "col-mono"), _tag("td", str(count))]
            rows_types.append(_tag("tr", "".join(row)))
    else:
        rows_types.append(
            "<tr><td colspan='2' style='text-align:center;color:#999'>No specific Type errors found.</td></tr>"
        )

    # -- Failure Stats Table --
    rows_stats = []
    if totally_skipped_routines > 0:
        for reason, count in failure_counts.most_common():
            pct = count / totally_skipped_routines * 100
            row = [
                _tag("td", html.escape(reason)),
                _tag("td", str(count)),
                _tag("td", f"{pct:.1f}%"),
            ]
            rows_stats.append(_tag("tr", "".join(row)))
    else:
        rows_stats.append(
            "<tr><td colspan='3' style='text-align:center;color:#28a745'>100% Coverage reached!</td></tr>"
        )

    # --- Final Final Assembly ---
    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Codegen Coverage Report</title>
    <style>{CSS}</style>
</head>
<body>
    <h1>cppbmad / pybmad Wrapper Status</h1>
    
    {toc_html}
    {summary_html}

    <h2 id="sec-structs">Structure Definitions</h2>
    <p class="intro">Fields removed from structures due to unmapped types.</p>
    <input type="text" class="table-search" placeholder="Filter structures by name or module..." data-target="table-structs">
    <div class="table-container">
        <table id="table-structs">
            <thead><tr><th>Name</th><th>Module</th><th>Fields</th><th>Status</th><th>Dropped Attributes</th></tr></thead>
            <tbody>{"".join(rows_structs)}</tbody>
        </table>
    </div>

    <h2 id="sec-types">Untranslated Types</h2>
    <p class="intro">
      These missing Type definitions block entire routines from being wrapped. 
      (Counts reflect only routines heavily blocked by this type).
    </p>
    <div class="table-container" style="max-width: 600px;">
        <table id="table-types">
            <thead><tr><th>Missing C++ Type</th><th>Blocking Impact (Unique Routines)</th></tr></thead>
            <tbody>{"".join(rows_types)}</tbody>
        </table>
    </div>

    <h2 id="sec-blockers">Implementation Blockers</h2>
    <p class="intro">Why are routines completely skipped? (Excludes routines that have at least one working definition).</p>
    <div class="table-container" style="max-width: 800px;">
        <table id="table-blockers">
            <thead><tr><th>Reason</th><th>Unique Routines Affected</th><th>% of Missing</th></tr></thead>
            <tbody>{"".join(rows_stats)}</tbody>
        </table>
    </div>

    <h2 id="sec-routines">Routine Definitions</h2>
    <p class="intro">
      List of all scanned Fortran routines. Grouped by name. <br>
      <span class="badge badge-subtle">SKIP</span> indicates a definition was skipped, but another valid definition for this name exists (so it's safe to ignore).
      <span class="badge badge-err">SKIP</span> indicates the routine is completely unavailable in Python.
    </p>
    <input type="text" class="table-search" placeholder="Search routines ("my_func", "bmad_mod")..." data-target="table-routines">
    <div class="table-container">
        <table id="table-routines">
            <thead><tr><th>Routine Name</th><th>Location</th><th>Status</th><th>Issues</th></tr></thead>
            <tbody>{"".join(rows_routines)}</tbody>
        </table>
    </div>
    
    {SCRIPT_JS}
</body>
</html>
    """
