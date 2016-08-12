## ZAZMOL
## ======

This is the unstable developer branch of sasmol

Library for defining molecular objects to create simulation and analysis programs

To install:

python setup.py install

dependencies:

numpy,
mocker

***

### [`basic usage`](doc_files/basic_usage.md)

### [`developer notes`](development_tools/notes.md)

[]([file i/o notes](development_tools/file_io_experiments.md))

***


```python
import sasmol.sasmol as sasmol
```


```python
molecule = sasmol.SasMol()
```


```python
molecule.read_pdb('hiv1_gag.pdb')
```

    reading filename:  hiv1_gag.pdb
    num_atoms =  6730
    >>> found  1  model(s) or frame(s)
    finished reading frame =  1



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




    array([ 88.3,  19.6,  14.7])




```python
molecule.moveto(frame,[88.3, 19.6, 14.7])
```


```python
molecule.calculate_center_of_mass(frame)
```




    array([ 88.3,  19.6,  14.7])




```python
molecule.rotate(frame,'x',45)
```


```python
molecule.calculate_center_of_mass(frame)
```




    array([ 88.3       , -21.92399383,  10.93565245])




```python
molecule.calculate_principle_moments_of_inertia(frame)
```




    (array([  1.30834716e+07,   1.91993314e+08,   1.85015201e+08]),
     array([[-0.08711655, -0.97104917,  0.22242802],
            [-0.53401862, -0.14296585, -0.8332976 ],
            [-0.84097255,  0.19137472,  0.50610364]]),
     array([[  1.90290278e+08,  -7.02983463e+06,  -1.38929432e+07],
            [ -7.02983463e+06,   1.36127046e+08,  -7.74046010e+07],
            [ -1.38929432e+07,  -7.74046010e+07,   6.36746633e+07]]))




```python
molecule.write_pdb('rotated_hiv1_gag.pdb',frame,'w');
```
