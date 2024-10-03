import numpy
import dcdio_module as dcdio    
import time
import sys
import sasmol.system as system

# Initialize the molecule and read the PDB file
A = system.Molecule(0)
A.read_pdb('hiv1_gag.pdb')

# Get the number of atoms
natoms = A.natoms()

# Extract initial coordinates
x = A.coor()[:, 0]
y = A.coor()[:, 1]
z = A.coor()[:, 2]

# Convert coordinates to numpy arrays
x = numpy.array(x, numpy.float32)
y = numpy.array(y, numpy.float32)
z = numpy.array(z, numpy.float32)

# Open the DCD file using dcdio.open_dcd_read
filename = 'h.dcd'
fp = dcdio.open_dcd_read(filename)
print('fp = ', fp)

# Initialize variables
nnatoms = 0
nset = 0
istart = 0
nsavc = 0
delta = 0.0
namnf = 0
freeindexes = []
reverseEndian = 0
charmm = 0

# Read the header to get the necessary information
print('Reading DCD header...')
header_result = dcdio.read_dcdheader(fp, nnatoms, nset, istart, nsavc, namnf, delta, delta, reverseEndian, charmm, freeindexes)

# Extract the values from the result
nnatoms = int(header_result[0])
nset = int(header_result[1])
istart = int(header_result[2])
nsavc = int(header_result[3])
delta = header_result[4]
namnf = int(header_result[5])
reverseEndian = int(header_result[6])
charmm = int(header_result[7])
num_fixed = int(header_result[8])

# Print the results
print("N:", nnatoms)
print("NSET:", nset)
print("ISTART:", istart)

# Initialize arrays to hold coordinates for all frames
x = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
y = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
z = numpy.zeros((nset, nnatoms), dtype=numpy.float32)

# Read the DCD file
print('Reading DCD file...')
start_time = time.time()
sum = 0.0
first = 1  # Initialize the 'first' variable

for i in range(nset):
    print('.', end=' ')
    sys.stdout.flush()
    read_start_time = time.time()
    tx = numpy.zeros(nnatoms, dtype=numpy.float32)
    ty = numpy.zeros(nnatoms, dtype=numpy.float32)
    tz = numpy.zeros(nnatoms, dtype=numpy.float32)

    try:
        final_result = dcdio.read_dcdstep(fp, nnatoms, num_fixed, first, reverseEndian, charmm, tx, ty, tz)
        if not final_result:
            print("Reached the end of the file or failed to read DCD step")
            break
    except RuntimeError as e:
        print(f"Error reading DCD step: {e}")
        print('step = ', i)
        break

    read_end_time = time.time()
    sum += read_end_time - read_start_time

    x[i][:] = tx
    y[i][:] = ty
    z[i][:] = tz

    first = 0  # Only the first frame is considered the first

# Close the file
dcdio.close_dcd_read(fp)

# Print the total time taken
end_time = time.time()
print(f"\nTotal time taken: {end_time - start_time} seconds")