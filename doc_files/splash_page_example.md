

```python
import sasmol.sasmol as sasmol
```


```python
molecule = sasmol.SasMol('hiv1_gag.pdb')
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




    array([ -6.79114736, -23.71577133,   8.06558513])




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




    array([ 88.3       ,  -2.21197083,  24.39994232])




```python
molecule.calculate_principle_moments_of_inertia(frame)
```




    (array([  1.30834716e+07,   1.91993314e+08,   1.85015201e+08]),
     array([[-0.08711655, -0.97104917,  0.22242802],
            [-0.99611823,  0.08773832, -0.00710418],
            [ 0.01261695,  0.2221835 ,  0.97492323]]),
     array([[  1.90290278e+08,  -1.55144810e+07,  -1.31655748e+06],
            [ -1.55144810e+07,   1.44693990e+07,   2.29686530e+06],
            [ -1.31655748e+06,   2.29686530e+06,   1.85332310e+08]]))




```python
molecule.write_pdb('rotated_hiv1_gag.pdb',frame,'w');
```
