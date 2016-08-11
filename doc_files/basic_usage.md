

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
import os ; os.chdir('/Users/curtisj/subversion_working_copies/svn_utk/sassie_2.0/trunk/sassie/simulate/monte_carlo')
```


```python
m.calculate_center_of_mass(0)
```




    array([  7.11544707e-13,   2.48159571e-12,  -8.45832820e-13])




```python
import gui_mimic_ten_mer
```


    ---------------------------------------------------------------------------

    ImportError                               Traceback (most recent call last)

    <ipython-input-12-7a9b78798219> in <module>()
    ----> 1 import gui_mimic_ten_mer
    

    /Users/curtisj/subversion_working_copies/svn_utk/sassie_2.0/trunk/sassie/simulate/monte_carlo/gui_mimic_ten_mer.py in <module>()
          6 import sys
          7 
    ----> 8 import sassie.simulate.monte_carlo.monte_carlo as monte_carlo
          9 import sassie.interface.input_filter as input_filter
         10 import multiprocessing


    ImportError: No module named sassie.simulate.monte_carlo.monte_carlo



```python

```
