import numpy as np
import math


t1, i1, x1, y1, z1 = np.loadtxt("outGR2_001B", usecols=(0,1,4,5,6), skiprows=16800, unpack=True)
t2, i2, x2, y2, z2 = np.loadtxt("out10", usecols=(0,1,2,3,4), skiprows=0, unpack=True)

print(t1[0], i1[0], x1[0], y1[0], z1[0]);
print(t2[0], i2[0], x2[0], y2[0], z2[0]);

for i in range(len(t2)):
	dx = x1[i] - x2[i]
	dy = y1[i] - y2[i]
	dz = z1[i] - z2[i]

	d = math.sqrt(dx * dx + dy * dy + dz * dz)

	print(t1[i], i1[i], d * 149597870700.0)
