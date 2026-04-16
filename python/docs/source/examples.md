# Examples

Example scripts are in [`python/examples/`](https://github.com/ken-lauer/cppbmad/tree/master/python/examples).

## Lattice Basics

### simple.py

Port of the Fortran `simple_bmad_program` to pybmad. Reads a lattice, modifies
a solenoid's strength, recalculates transfer matrices, propagates Twiss
parameters, and prints element info.

```python
--8<-- "python/examples/simple.py"
```

### parse.py

Parse a Bmad lattice file and iterate over elements, printing names and properties.

```python
--8<-- "python/examples/parse.py"
```

### lat_ele_locator.py

Select elements within a lattice by name or wildcard pattern using `lat_ele_locator`.

```python
--8<-- "python/examples/lat_ele_locator.py"
```

## Beam Tracking

### csr.py

Track a particle distribution through a lattice with coherent synchrotron radiation
(CSR) and space charge effects. Shows how to configure `bmad_com`,
`space_charge_com`, initialize a beam, and call `track_beam`.

```python
--8<-- "python/examples/csr.py"
```

### bbu.py

Beam breakup (BBU) instability analysis: parse a lattice, compute Twiss
parameters, hybridize the lattice, and track with `bbu_track_all`.

```python
--8<-- "python/examples/bbu.py"
```

## Geometry

### wall_generator.py

Extract wall geometry from `wall3d` lattice sections. Computes wall radius and
surface normals at angular sample points and writes output files.

```python
--8<-- "python/examples/wall_generator.py"
```

## Tao

pybmad is also able to interface with Tao's data structures directly.
Loading a lattice with PyTao, we can poke into all of the Bmad and Tao
structures directly from Python.

Interaction with `pytao.Tao` instances is supported natively, but
`pytao.SubprocessTao` requires more effort to work with - as Tao's data is no
longer in the same process. See the `tao-subproc.py` example below for how to
work around this.

### tao.py

Initialize PyTao and access universe/lattice data structures via `get_super_universe`.

```python
--8<-- "python/examples/tao.py"
```

### tao-plot.py

Extract graph and curve data from a Tao session and plot with matplotlib.

```python
--8<-- "python/examples/tao-plot.py"
```

### tao-subproc.py

Run Tao in a subprocess, interact with Bmad/Tao structures using pybmad in that
subprocess, and return data to the parent process.

This example is in 2 parts `tao-subproc.py` and `tao_subproc_helper.py`. This
is due to a limitation in PyTao - it requires that custom functions be in an
importable module (which `__main__` is not).

```python
--8<-- "python/examples/tao-subproc.py"
```

**tao_subproc_helper.py**

Helper module for `tao-subproc.py` -- provides a callable that accesses
universe and lattice data from within the subprocess.

```python
--8<-- "python/examples/tao_subproc_helper.py"
```
