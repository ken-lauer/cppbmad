## cppbmad and pybmad

`cppbmad`: C++ interface to lower level [Bmad](https://github.com/bmad-sim/bmad-ecosystem) and Tao structures and routines.

`pybmad`: cppbmad wrapped for Python using pybind11.

## Code Generation Overview

`codegen` is a Python package which generates:

1. Proxy class C++ headers and source code
2. Fortran-side accessors (and setters)
3. C++ enum equivalents of Bmad constants. For example, the `proton$` parameter
   on the Fortran side is translated to `PROTON` on the C++ side, while element
   attributes get special treatment with the `enum class` `EleAttribute`.
4. Pybind11 wrappers for proxy classes, enums, routines, and Python-only
   return-value structures (due to immutability of Python booleans and such).

It relies on a not-yet-merged branch upstream in [Bmad](), including Fortran JSON
support and `struct` parsing utilities.

## Wrapping other codebases?

This approach could potentially be generalized for other Fortran codebases.
It would be neat one day to do that.

- I'm admittedly no Fortran expert to see where the limitations are. I've built
  this entirely around how Bmad uses Fortran - the parser, the routine
  handling, and everything else.
- I don't have any other projects in mind at the moment.
- It's doubtful any non-Bmad people would stumble across this.
- The codebase is a mess (not just humble developer speak)

So in all likelihood this will stay as a Bmad-only C++/Python wrapper
generator.

## Building

1. Clone [Bmad](https://github.com/bmad-sim/bmad-ecosystem) (currently: my JSON branch)
2. Set `ACC_ROOT_DIR` to that directory
3. Run `bash regenerate.sh`

This script:

1. Runs `codegen` to update all the bindings
2. Runs `cmake`
3. Builds the project
4. Runs a few automated tests

Use `-DBUILD_PYBMAD=OFF` to disable building the Python bindings.

### pybmad installation

`scikit-build` is used under the hood to compile and install pybmad for you.

```
python -m codegen
CMAKE_BUILD_PARALLEL_LEVEL=8 python -m pip install .
```

### pybmad development

For pybmad development, in-tree building is possible.

I use the following locally:

```
cd python/pybmad
ln -s ../../debug/python/_pybmad.so
```
