from __future__ import annotations

from .gen import main

ctx = main()

structs, parsed_structs, routines, routines_by_name = (
    ctx.codegen_structs,
    ctx.parsed_structs,
    ctx.routines,
    ctx.routines_by_name,
)
