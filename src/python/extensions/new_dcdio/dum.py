import dcdio_module

# Open the DCD file
fd = dcdio_module.open_dcd_read("hiv1_gag_200_frames.dcd")

# Call the read_dcdheader method
result = dcdio_module.read_dcdheader(fd, 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0)

# Print the results
print("N:", result[0])
print("NSET:", result[1])
print("ISTART:", result[2])
print("NSAVC:", result[3])
print("NAMNF:", result[4])
print("DELTA:", result[5])
print("Data:", result[6])
print("Extra Arg:", result[7])
print("Reverse Endian:", result[8])
print("Charmm:", result[9])