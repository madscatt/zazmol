import math
import os


def determine_float_type():
    '''
       determine whether it will be a 32/64-bit number for the float type
    '''
    a = 1.e40
    if math.isinf(a):
        return 'float32'
    # set up a pseudo-32bit machine for testing
    elif a == 1.e40:
        return 'float64'
    # Need a reasonable default
    else:
        return 'float'


os.environ.setdefault('SASMOL_LARGETEST', 'n')
os.environ.setdefault('SASMOL_HUGETEST', 'n')


# Genrate huge dcd files if SASMOL_HUGETEST is specified
# IMPORTANT, SASMOL_HUGETEST only implies the huge rna dcd files generated below will be used
if os.environ['SASMOL_HUGETEST'] == 'y':
    from sasmol.test_sasmol.utilities import generate_huge_dcd_onthefly
    generate_huge_dcd_onthefly.generate_huge_dcd()


os.environ.setdefault('SASMOL_FLOATTYPE', determine_float_type())
