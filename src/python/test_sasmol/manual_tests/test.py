import os
import sasmol.sasmol as sasmol
import sasmol.calculate as calculate


pdb_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

m = sasmol.SasMol()

try:
    m.read_pdb("hiv1_gag.pdb")
except:
    m.read_pdb(pdb_file_path+"hiv1_gag.pdb")

calc = calculate.Calculate

print 'calc.calculate_center_of_mass = ',calc.calculate_center_of_mass(m,0)

calc_com = calc.calculate_center_of_mass

calc_rg = calc.calculate_radius_of_gyration

print 'calc_com = ',calc_com(m,0)
     
print 'calc_rg = ',calc_rg(m,0)

print 'm.calculate_mass() = ', m.calculate_mass()


'''

    You just put them in __init__.py.

    So with test/classes.py being:

    class A(object): pass
    class B(object): pass
    ... and test/__init__.py being:

    from classes import *

    class Helper(object): pass
    You can import test and have access to A, B and Helper

    >>> import test
    >>> test.A
    <class 'test.classes.A'>
    >>> test.B
    <class 'test.classes.B'>
    >>> test.Helper
    <class 'test.Helper'>

'''
