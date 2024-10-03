'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import sasmol.system as system
#from . import dcdio
import dcdio_module as dcdio
import sys
import numpy
import time
#sys.path.append('./')

A = system.Molecule(0)
A.read_pdb('min3.pdb')

natoms = A.natoms()

x = A.coor()[:, 0]
y = A.coor()[:, 1]
z = A.coor()[:, 2]

#x=[1.0,2.0,3.0] ; y=[1.0,2.0,3.0] ; z=[1.0,2.0,3.0]

x = numpy.array(x, numpy.float32)
y = numpy.array(y, numpy.float32)
z = numpy.array(z, numpy.float32)

filename = 'hiv1_gag_200_frames.dcd'
'''
fp = dcdio.open_dcd_write(filename)

nset = 200
istart = 1
nsavc = 1
delta = 1.0

headerresult = dcdio.write_dcdheader(
    fp, filename, natoms, nset, istart, nsavc, delta)

print('writing '+str(nset)+' to disk')
start_time = time.time()

for blah in range(nset):
    print(".", end=' ')
    sys.stdout.flush()

    x = x+5.0
    y = y+5.0
    z = z+5.0
    stepresult = dcdio.write_dcdstep(fp, x, y, z, blah)

end_time = time.time()

dt = end_time-start_time

print('\ntotal time = ', dt, ' time per structure = ', dt/nset)

dcdio.close_dcd_write(fp)


filename = '200c.dcd'
filename = 'c7.dcd'

'''
filename = 'hiv1_gag_200_frames.dcd'
fd = dcdio.open_dcd_read(filename)

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

#readheaderresult, nnatoms, nset, istart, nsavc, delta, namnf, reverseEndian, charmm = dcdio.read_dcdheader(fd, 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0)
result = dcdio.read_dcdheader(fd, 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0)

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

'''
print('read header result = ', readheaderresult)
print('nnatoms = ', nnatoms)
print('nset = ', nset)
print('istart = ', istart)
print('nsavc = ', nsavc)
print('delta = ', delta)
print('namnf = ', namnf)
print('reverseEndian = ', reverseEndian)
print('charmm = ', charmm)
'''
nnatoms = result[0] 
nset = result[1]
istart = result[2]
nsavc = result[3]
namnf = result[4]
delta = result[5]   
data = result[6]
reverseEndian = result[8]
charmm = result[9]


x = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
y = numpy.zeros((nset, nnatoms), dtype=numpy.float32)
z = numpy.zeros((nset, nnatoms), dtype=numpy.float32)

num_fixed = 0
first = 1
result = 1

i = 0

# try:
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

    #result = dcdio.read_dcdstep(
    #final_result = dcdio.read_dcdstep(
    #    fd, tx, ty, tz, num_fixed, i, reverseEndian, charmm)
    
    final_result = dcdio.read_dcdstep(fd, nnatoms, num_fixed, first, reverseEndian, charmm, tx, ty, tz)

    print("Final Result:", final_result)
    read_end_time = time.time()

    sum += read_end_time-read_start_time

    x[i][:] = tx
    y[i][:] = ty
    z[i][:] = tz

end_time = time.time()
dt = end_time-start_time

print('\nread total_time = ', sum, ' time per structure = ', sum/nset)
print('total_time = ', dt, ' time per structure = ', dt/nset)
print('ratio(total time) = ', dt/sum)
print('ratio(per structure) = ', (dt/nset)/(sum/nset))


# except:
#
#	print " I failed          :(     "

dcdio.close_dcd_read(ifp)
