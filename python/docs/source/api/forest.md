# Forest / PTC

PTC (Polymorphic Tracking Code) and Forest library.

## Classes (Fortran Structures)

::: pybmad.Fibre
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Fibre

Fortran struct: `fibre` ([`forest/include/a_def_element_fibre_layout.inc`](https://github.com/bmad-sim/bmad-ecosystem/blob/5132fc8eb2553dd5976821cea40e62ba5bab94e2/forest/include/a_def_element_fibre_layout.inc#L385))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `DIR` | int |  |
| `PREVIOUS` | [Fibre](forest.md#fibre) |  |
| `NEXT` | [Fibre](forest.md#fibre) | POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA |
| `PARENT_LAYOUT` | [Layout](forest.md#layout) |  |
| `pos` | int | POSITION IN LAYOUT NEW STUFF.... |
| `BETA0` | float | ,P0C |
| `GAMMA0I` | float | ,P0C |
| `GAMBET` | float | ,P0C |
| `MASS` | float | ,P0C |
| `CHARGE` | float |  |
| `AG` | float | spin g-2 TO TIE LAYOUTS |
| `P` | [Fibre](forest.md#fibre) | tying them in the so-called database universe M_u |
| `N` | [Fibre](forest.md#fibre) |  |
| `loc` | int |  |

::: pybmad.Layout
    options:
      heading_level: 0
      show_root_heading: false
      members: false
      show_signature: false
      show_bases: false
      show_docstring_description: false

### Layout

Fortran struct: `layout` ([`forest/include/a_def_element_fibre_layout.inc`](https://github.com/bmad-sim/bmad-ecosystem/blob/5132fc8eb2553dd5976821cea40e62ba5bab94e2/forest/include/a_def_element_fibre_layout.inc#L414))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `NAME` | str | IDENTIFICATION |
| `INDEX` | int | IDENTIFICATION, CHARGE SIGN |
| `HARMONIC_NUMBER` | float |  |
| `CLOSED` | bool |  |
| `N` | int | TOTAL ELEMENT IN THE CHAIN |
| `NTHIN` | int | NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES) |
| `THIN` | float | PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS POINTERS OF LINK LAYOUT |
| `LASTPOS` | int | POSITION OF LAST VISITED |
| `LAST` | [Fibre](forest.md#fibre) | LAST VISITED |
| `END` | [Fibre](forest.md#fibre) |  |
| `START` | [Fibre](forest.md#fibre) |  |
| `START_GROUND` | [Fibre](forest.md#fibre) | STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING |
| `END_GROUND` | [Fibre](forest.md#fibre) | STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING |
| `NEXT` | [Layout](forest.md#layout) |  |
| `PREVIOUS` | [Layout](forest.md#layout) |  |
