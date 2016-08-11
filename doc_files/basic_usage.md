

```python
import sasmol.sasmol as s
```


```python
m = s.SasMol(0)
```


```python
m.read_pdb('hiv1_gag.pdb')
```

    reading filename:  hiv1_gag.pdb
    num_atoms =  6730
    >>> found  1  model(s) or frame(s)
    finished reading frame =  1



```python
m
```




    sasmol object




```python
m.names()[:3]
```




    ['N', 'HT1', 'HT2']




```python
m.calculate_center_of_mass(0)
```




    array([ -6.79114736, -23.71577133,   8.06558513])




```python
m.center(0)
```


```python
m.calculate_center_of_mass(0)
```




    array([  7.11544707e-13,   2.48159571e-12,  -8.45832820e-13])


