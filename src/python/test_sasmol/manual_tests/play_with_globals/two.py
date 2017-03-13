import sys

import sasmol.sasconfig as sasconfig

print 'in test: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__
print 'in test: sasconfig.__logging__ = ', sasconfig.__logging__

class two_class():
    def __init__(self):
        self._id = sasconfig.__sasmol_id__ + 1
        sasconfig.__sasmol_id__ += 1

    def __del__(self):
        sasconfig.__sasmol_id__ -= 1
        del self
    
print ; print 'IN TWO' ; print

sasconfig.__logging__ = True

print 'creating instance a'
a = two_class()
#with two_class() as a
print 'in two: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__

print 'creating instance b'
b = two_class()
print 'in two: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__

del a,b


