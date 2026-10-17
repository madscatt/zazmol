
### Basic Usage

#### This section highlights basic usage from reading in atomic information from a protein data bank (PDB) file, extracting information from the molecule, performing simple operations on the molecule and how to save the new coordinates to a new PDB file
***

#### Import sasmol library


```python
import sasmol.system as system
```

#### Define "molecule" an instance of the Molecule class.  This is the "object" that will hold all information for the molecule


```python
molecule = system.Molecule()
```

#### Read  the contents of a PDB file into the molecule


```python
molecule.read_pdb('hiv1_gag.pdb')
```

    reading filename:  hiv1_gag.pdb
    num_atoms =  6730
    >>> found  1  model(s) or frame(s)
    finished reading frame =  1


#### Query the number of atoms


```python
molecule.natoms()
```




    6730



####  Determine the center of mass of the molecule for frame = 0


```python
molecule.calculate_center_of_mass(0) 
```




    array([ -6.79114736, -23.71577133,   8.06558513])



#### Set the center of mass to [0, 0, 0] for frame = 0


```python
molecule.center(0)
```

#### Check that the center of mass is indeed now at [0, 0, 0]


```python
molecule.calculate_center_of_mass(0)
```




    array([  7.11544707e-13,   2.48159571e-12,  -8.45832820e-13])



#### Save the coordinates to a new PDB file


```python
molecule.write_pdb("centered_hiv1_gag.pdb", 0, "w");
```


### Prepared atom selections

The separate `get_subset_mask_prepared()` method prepares a selection expression
once per call before evaluating it for each atom. Use it explicitly when repeated
expression translation is a measured cost:

```python
import sasmol.system as system

molecule = system.Molecule_Maker(3, name=['N', 'CA', 'C'])
error, mask = molecule.get_subset_mask_prepared('name[i] == "CA"')
if not error:
    indices = molecule.get_indices_from_mask(mask)  # array([1])
```

`get_subset_mask()` and existing callers are unchanged. Both methods read current
molecular descriptors and retain atom order. The prepared method stores no masks,
indices, molecular state or expressions between calls. It still traverses every
atom on each selection; it does not make a repeated full-molecule search disappear.
Preparation is a Python operation and needs no native compiler or new dependency.

As with the existing method, developers supply suitable expression strings and
handle returned results or exceptions. Successful masks are integer NumPy arrays;
existing no-match and partial-failure results are lists with selection messages.
No new error channel, input-validation layer or automatic fallback is introduced.
The prepared method is intended for molecular descriptor expressions, not exact
parity of arbitrary Python frame/local introspection or translation audit events.
A subclass that customizes `get_subset_mask()` must review the separate method
before using it: the prepared method does not call that custom override. Existing
`get_named_subset_mask()` calls continue using the original evaluator.
