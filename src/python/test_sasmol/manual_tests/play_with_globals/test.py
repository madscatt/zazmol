import sys
import weakref

sys.path.append('./')

import sasmol.sasmol as sasmol
import sasmol.sasconfig as sasconfig

sasconfig.__sasmol_id__ += 1

print 'in test: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__
print 'in test: sasconfig.__logging__ = ', sasconfig.__logging__
print 

for i in xrange(10):
    sasconfig.__sasmol_id__ += 1

print 'in test: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__
print 

print
print 'calling two'
print

import two as two

print
print 'back in test' ; print

print 'in test: sasconfig.__sasmol_id__ = ', sasconfig.__sasmol_id__
print 'in test: sasconfig.__logging__ = ', sasconfig.__logging__

print

