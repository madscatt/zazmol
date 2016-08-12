
#### Import sasmol library


```python
import sasmol.sasmol as sasmol
```

#### Define "molecule" an instance of the SasMol class.  This is the "object" that will hold all information for the molecule


```python
molecule = sasmol.SasMol()
```

#### Read  the contents of a PDB file into the molecule


```python
molecule.read_pdb('hiv1_gag.pdb')
```

#### Query the number of atoms


```python
molecule.natoms()
```

####  Determine the center of mass of the molecule for frame = 0


```python
molecule.calculate_center_of_mass(0) 
```

#### Set the center of mass to [0, 0, 0] for frame = 0


```python
molecule.center(0)
```

#### Check that the center of mass is indeed now at [0, 0, 0]


```python
molecule.calculate_center_of_mass(0)
```
