# Hooks

Bmad and Tao expose a number of _hook_ procedure pointers that let user code run
at specific points during tracking, lattice calculation, or optimization. pybmad
surfaces the supported ones as **read/write properties** on a per-submodule
registry object:

```python
import pybmad as pb

def postprocess(
    start_orb: pb.CoordStruct,
    ele: pb.EleStruct,
    param: pb.LatParamStruct,
    end_orb: pb.CoordStruct,
):
    # Observe the exit orbit after each element. end_orb is a live proxy.
    print(f"Ele: {ele.name} Starting orbit: {start_orb.vec} Ending orbit: {end_orb.vec}")

pb.hooks.track1_postprocess = postprocess
result = pb.track_all(lat, orbit, ix_branch)
pb.hooks.track1_postprocess = None  # clear the hook
```

`pybmad.hooks` is an alias for `pybmad.bmad.hooks`.
Tao hooks live under `pybmad.tao.hooks`.

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

<!-- prettier-ignore-start -->

::: pybmad.bmad.BmadHooks
  options:
    show_root_heading: false
    show_root_toc_entry: false
    members_order: source
    show_if_no_docstring: false
<!-- prettier-ignore-end -->

## Tao hooks

`pybmad.tao.hooks` — the lattice-calculation, optimizer, and merit-function hooks.

<!-- prettier-ignore-start -->
::: pybmad.tao.TaoHooks
  options:
    show_root_heading: false
    show_root_toc_entry: false
    members_order: source
    show_if_no_docstring: false
<!-- prettier-ignore-end -->

!!! note
These bindings are hand-written and maintained (only function prototypes are code-generated). Not all are supported.
