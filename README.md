
# SASMOL
# ======

import sasmol.system as system
This is the unstable developer branch of sasmol

Library for defining molecular objects to create simulation and analysis programs

To install:

python setup.py install

dependencies:

- Python
- NumPy

***

# [`basic usage`](doc_files/basic_usage.md)

# [`library documentation`](https://madscatt.github.io/zazmol/index.html)

# [`developer notes`](doc_files/developer_notes.md)

***


```python
```


```python
molecule = system.Molecule('hiv1_gag.pdb')
```

reading filename:  hiv1_gag.pdb
num_atoms = 6730
>> > found  1  model(s) or frame(s)
finished reading frame = 1


```python
molecule.natoms()
```

6730


```python
molecule.number_of_frames()
```

1


```python
frame = 0
```


```python
molecule.calculate_center_of_mass(frame)
```

array([-6.79114736, -23.71577133,   8.06558513])


```python
molecule.translate(frame, [88.3, 19.6, 14.7], point=True)
```


```python
molecule.calculate_center_of_mass(frame)
```

array([88.3,  19.6,  14.7])


```python
molecule.rotate(frame, 'x', 45*3.1515927/180.0)
```


```python
molecule.calculate_center_of_mass(frame)
```

array([88.3,   3.40417778,  24.26234889])


```python
molecule.calculate_principal_moments_of_inertia(frame)
```

(array([1.30834716e+07,   1.91993314e+08,   1.85015201e+08]),
 array([[-0.08711655, -0.97104917,  0.22242802],
        [-0.9670572,  0.13604234,  0.21515775],
        [0.23918838,  0.19635682,  0.95091162]]),
 array([[1.90290278e+08,  -1.54065145e+07,   2.25205595e+06],
        [-1.54065145e+07,   2.43538600e+07,   3.99557354e+07],
        [2.25205595e+06,   3.99557354e+07,   1.75447849e+08]]))


```python
molecule.write_pdb('rotated_hiv1_gag.pdb', frame, 'w')
```

## Installed build identity

Normal `setup.py build` / installation completes the existing native build and
includes `sasmol/_build_info.json`. The owner helper is `build_support.py`;
component ID is `zazmol`, canonical package version remains `2.0`. The schema-1
stamp carries the source revision (or null), tracked-source modified flag and
UTC build time. Eggs carry it as a package resource. Source distributions carry
the helper and `.build-info.json`, so a source archive can be built without Git.
A failed build does not create a success stamp. Missing evidence stays unknown.

Read installed evidence through `importlib.resources` in the interpreter that
runs the jobs, rather than reading a nearby checkout. Normal
`~/anaconda3/bin/python setup.py install` now verifies the installed
stamp/distribution, native extensions, modern API, headers and static library
in an isolated child interpreter. After verification it automatically updates
only ZazMol in an explicitly bound website's stack. No separate publish command
is required. The optional one-time binding lives at
`<selected Python prefix>/etc/sassie/version-stack.json`; its schema is owned by
`../genapp_zazzie/docs/source_revision_metadata.md`. No application dependency
is required for standalone installations without that binding.

Plain builds, eggs/archives and staged or alternate-target installs (`--root`,
`--home`, `--user`, alternate prefix/library directory) do not update a running
instance. A failed install or native verification cannot publish. A binding or
publication error fails the command with `Installed successfully; version-stack
update failed`; the verified new package remains installed. This does not add
rollback to the normal setuptools installer. Deployment and binding setup are
separate approved operations. Publication preserves
independent release dates/tags and reports modified or mixed builds explicitly.
Git hooks do not publish installed versions. Live hook migration and deployment
are separately authorized operations. These changes do not alter scientific
calculations, native build options, or the existing release version.
