

import dcdio_module as dcdio
#import ctypes

filename = 'hiv1_gag_200_frames.dcd'

filepointer_capsule = dcdio.open_dcd_file(filename)


if filepointer_capsule:
    print(f"Successfully opened {filename}")
else:
    print(f"Failed to open {filename}")

N = 0
NSET = 0
ISTART = 0
NSAVC = 0
NAMNF = 0
DELTA = 0.0
data = 0.0
extra_arg = 0
reverseEndian = 0
charmm = 0

# Call the read_dcdheader function
#result = dcdio_module.read_dcdheader(file_ptr, N, NSET, ISTART, NSAVC, NAMNF, DELTA, data, extra_arg, reverseEndian, charmm)
readheaderresult, filepointer_capsule, nnatoms, nset, istart, nsavc, delta, namnf, reverseEndian, charmm = dcdio.read_dcdheader(filepointer_capsule)

# Print the results
print("readheaderresult, filepointer, nnatoms, nset, istart, nsavc, delta, namnf, reverseEndian, charmm")
print(readheaderresult, filepointer_capsule, nnatoms, nset, istart, nsavc, delta, namnf, reverseEndian, charmm)
