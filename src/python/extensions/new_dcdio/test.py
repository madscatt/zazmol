

import dcdio_module

filename = 'hiv1_gag_200_frames.dcd'

file_ptr = dcdio_module.open_dcd_read(filename)
if file_ptr:
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
result = dcdio_module.read_dcdheader(file_ptr, N, NSET, ISTART, NSAVC, NAMNF, DELTA, data, extra_arg, reverseEndian, charmm)

# Print the results
print("Result:", result[0])
print("N:", result[1])
print("NSET:", result[2])
print("ISTART:", result[3])
print("NSAVC:", result[4])
print("NAMNF:", result[5])
print("DELTA:", result[6])
print("Data:", result[7])
print("Extra Arg:", result[8])
print("Reverse Endian:", result[9])
print("Charmm:", result[10])
