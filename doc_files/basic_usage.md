
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
