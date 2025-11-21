## cppbmad

C++ interface to lower level [Bmad](https://github.com/bmad-sim/bmad-ecosystem) and Tao structures.

Wrapped for Python in [pybmad](https://github.com/ken-lauer/pybmad).

## Code Generation Overview

`codegen` is a Python package which generates:

1. Proxy class C++ headers and source code
2. Fortran-side accessors (and setters)
3. C++ enum equivalents of Bmad constants. For example, the `proton$` parameter
   on the Fortran side is translated to `PROTON` on the C++ side, while element
   attributes get special treatment with the `enum class` `EleAttribute`.

It relies on a not-yet-merged branch upstream in Bmad, including Fortran JSON
support and `struct` parsing utilities.
