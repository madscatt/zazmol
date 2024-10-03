import numpy
import dcdio_module as dcdio    
import time
import sys

import sasmol.system as system

A = system.Molecule(0)
A.read_pdb('hiv1_gag.pdb')

natoms = A.natoms()

x = A.coor()[:, 0]
y = A.coor()[:, 1]
z = A.coor()[:, 2]

#x=[1.0,2.0,3.0] ; y=[1.0,2.0,3.0] ; z=[1.0,2.0,3.0]

x = numpy.array(x, numpy.float32)
y = numpy.array(y, numpy.float32)
z = numpy.array(z, numpy.float32)

# Open the DCD file using dcdio.open_dcd_read
filename = 'h.dcd'
fp = dcdio.open_dcd_read(filename)
print('fp = ',fp)

nnatoms = 0
nset = 0
istart = 0
nsavc = 0
delta = 0.0
namnf = 0
freeindexes = []
reverseEndian = 0
charmm = 0

print('nnatoms = ', nnatoms)
print('nset = ', nset)
print('freeindexes = ', freeindexes)

# Read the header to get the necessary information
print('fp = ',fp)
result = dcdio.read_dcdheader(fp, 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0)

print('fp = ',fp)
# Extract the values from the result
nnatoms = int(result[0])
nset = int(result[1])
istart = int(result[2])
nsavc = int(result[3])
delta = result[4]
namnf = int(result[5])
reverseEndian = int(result[6])
charmm = int(result[7])
num_fixed = int(result[8])

# Initialize arrays to hold coordinates for all frames
x = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
y = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
z = numpy.zeros((nset, nnatoms), dtype=numpy.float32)

# Print the results
print("N:", nnatoms)
print("NSET:", nset)
print("ISTART:", istart)

i = 0
first = 0  # Initialize the 'first' variable

print('reading dcd file')
start_time = time.time()
sum = 0.0
for i in range(nset):
    print('.', end=' ')
    sys.stdout.flush()
    read_start_time = time.time()
    tx = numpy.zeros(nnatoms, dtype=numpy.float32)
    ty = numpy.zeros(nnatoms, dtype=numpy.float32)
    tz = numpy.zeros(nnatoms, dtype=numpy.float32)

    print('fp = ',fp)
    #final_result = dcdio.read_dcdstep(fp, nnatoms, num_fixed, first, reverseEndian, charmm, tx, ty, tz)
    try:
        final_result = dcdio.read_dcdstep(fp, nnatoms, num_fixed, first, reverseEndian, charmm, tx, ty, tz)
        if not final_result:
            print("Reached the end of the file or failed to read DCD step")
            break
    except RuntimeError as e:
        print(f"Error reading DCD step: {e}")
        print('step = ', i)
        pass


    #print("Final Result:", final_result)
    read_end_time = time.time()

    sum += read_end_time - read_start_time

    x[i][:] = tx
    y[i][:] = ty
    z[i][:] = tz

# Close the file
#dcdio.close_dcd_read(fp)
