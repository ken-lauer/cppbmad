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

Fortran struct: `fibre` ([`forest/include/a_def_element_fibre_layout.inc`](https://github.com/bmad-sim/bmad-ecosystem/blob/af2f5465f9fc1f180df62221d54ba1e8c7182fd2/forest/include/a_def_element_fibre_layout.inc#L387))

All attributes may be passed to the initializer as arguments:

| Attribute | Type | Description |
|-----------|------|-------------|
| `DIR` | int |  |
| `pos` | int | POSITION IN LAYOUT NEW STUFF.... |
| `BETA0` | float | ,P0C |
| `GAMMA0I` | float | ,P0C |
| `GAMBET` | float | ,P0C |
| `MASS` | float | ,P0C |
| `CHARGE` | float |  |
| `AG` | float | spin g-2 TO TIE LAYOUTS |
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

Fortran struct: `layout` ([`forest/include/a_def_element_fibre_layout.inc`](https://github.com/bmad-sim/bmad-ecosystem/blob/af2f5465f9fc1f180df62221d54ba1e8c7182fd2/forest/include/a_def_element_fibre_layout.inc#L416))

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
