import numpy as np
import struct

useHelioCentric = 1

if(useHelioCentric == 1):
	filename = 'All_h.dat'
	bfilename = 'All_h.bin'
else:
	filename = 'All_b.dat'
	bfilename = 'All_b.bin'



time, index, x, y, z, vx, vy, vz = np.loadtxt(filename, dtype='double', unpack=True)

f = open(bfilename,"wb")

for i in range(len(time)):

	if(i < 27 * 3):
		print(i, time[i], x[i], y[i], z[i])
	s = struct.pack('d', time[i])
	f.write(s)
	s = struct.pack('d', x[i])
	f.write(s)
	s = struct.pack('d', y[i])
	f.write(s)
	s = struct.pack('d', z[i])
	f.write(s)

f.close()
	
