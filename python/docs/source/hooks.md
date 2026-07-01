# Hooks

Bmad and Tao expose a number of _hook_ procedure pointers that let user code run
at specific points during tracking, lattice calculation, or optimization. pybmad
surfaces the supported ones as **read/write properties** on a per-submodule
registry object:

```python
import pybmad
from pybmad import bmad, tao

def on_postprocess(start_orb, ele, param, end_orb):
    # end_orb is a live proxy — mutate it to change the exit coordinates
    print(ele.name, end_orb.vec[0])

bmad.hooks.track1_postprocess = on_postprocess   # install (fires after every track1)
current = bmad.hooks.track1_postprocess          # read back the current callback
bmad.hooks.track1_postprocess = None             # clear
```

`pybmad.hooks` is an alias for `pybmad.bmad.hooks`; Tao hooks live under
`pybmad.tao.hooks`.

## Semantics

- **Assign** a callable to install a hook, assign `None` to clear it, and **read**
  the property to get the current callback (or `None`).
- **Arguments** are live, non-owning proxy/array views of the underlying Fortran
  objects. Mutating them (e.g. `end_orb.vec[...]`) affects the ongoing
  calculation. They are valid **only for the duration of the call** — do not stash
  them for later use.
- **Optional** arguments that are absent on the Fortran side are passed as `None`.
- **Out-parameters** are returned from the callback. Return `None` to leave them
  unchanged; otherwise return them in signature order (a tuple when there is more
  than one). Each property's docstring gives the exact signature and return value.
- **Exceptions** raised inside a callback are reported to stderr and swallowed —
  they never propagate back into Fortran (which would crash the process).
- The callback fires on the thread driving Bmad/Tao while the GIL is held; set and
  clear hooks from that same thread.

## Bmad tracking hooks

`pybmad.bmad.hooks` — nine tracking / custom hooks from Bmad's
`bmad_routine_interface`.

::: pybmad.bmad.BmadHooks
  options:
    show_root_heading: false
    show_root_toc_entry: false
    members_order: source
    show_if_no_docstring: false

## Tao hooks

`pybmad.tao.hooks` — the lattice-calculation, optimizer, and merit-function hooks.

::: pybmad.tao.TaoHooks
  options:
    show_root_heading: false
    show_root_toc_entry: false
    members_order: source
    show_if_no_docstring: false

!!! note
These bindings are hand-written and maintained (not code-generated). Not all are supported.
