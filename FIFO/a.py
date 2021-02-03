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

data = np.ones(2)
data[0] = 5.0
data[1] = 6.0

OK = np.zeros(1)

#fifo = open(path, "w")

fd = os.open('myfifo', os.O_WRONLY)

line = data.tobytes()

os.write(fd, line)

os.close(fd)


##############################

fd = open('fifoCheck', 'r')

line2 = fd.read()

line3= str.encode(line2)

d = int.from_bytes(line3, byteorder='little', signed=False)
print(d)

fd.close()

