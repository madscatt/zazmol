import os
import sasmol.sasmol as sasmol
import timeit

pdb_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

m = sasmol.SasMol(0)

try:
    m.read_pdb("hiv1_gag.pdb")
except:
    m.read_pdb(pdb_file_path+"hiv1_gag.pdb")

#dcdfile = 'big.dcd'

start_time = timeit.default_timer()
min_max_list = m.calculate_minimum_and_maximum_one_frame(0)
#min_max_list = m.calculate_minimum_and_maximum_all_frames(dcdfile)
final_time = timeit.default_timer()
without_time = final_time - start_time


start_time = timeit.default_timer()
try:
    min_max_list = m.calculate_minimum_and_maximum_one_frame(0)
#    min_max_list = m.calculate_minimum_and_maximum_all_frames(dcdfile)
except:
    pass

final_time = timeit.default_timer()
with_time = final_time - start_time

print 'without try/except: time = ', without_time
print 'with try/except: time = ', with_time



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
