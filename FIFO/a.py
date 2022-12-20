import numpy as np
import os
import struct

try:
	os.mkfifo('myfifo', 0o666)
except:
	pass

try:
	os.fifoCheck('fifoCheck', 0o666)
except:
	pass

N = 1

#################################
# send N
#################################
fd = os.open('myfifo', os.O_WRONLY)

line = N.to_bytes(2, 'little')

os.write(fd, line)

os.close(fd)
print('Sent', N)

##############################
# receive N to check
##############################

fd = open('fifoCheck', 'r')

line2 = fd.read()

line3= str.encode(line2)

N2 = int.from_bytes(line3, byteorder='little', signed=False)
print(N2)

fd.close()

##############################
# read initial conditions
##############################

time1, index, x, y, z, vx, vy, vz = np.loadtxt('initial1.dat', dtype='double', unpack=True)

fd = os.open('myfifo', os.O_WRONLY)

time = 2396770.5

print(time, index, x, y, z, vx, vy, vz)

#line = data.to_bytes(4, 'little')
line = bytearray(struct.pack("d", time1))
line += bytearray(struct.pack("i", int(index)))
line += bytearray(struct.pack("d", x))
line += bytearray(struct.pack("d", y))
line += bytearray(struct.pack("d", z))
line += bytearray(struct.pack("d", vx))
line += bytearray(struct.pack("d", vy))
line += bytearray(struct.pack("d", vz))

os.write(fd, line)

os.close(fd)


