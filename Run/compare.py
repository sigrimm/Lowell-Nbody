import numpy as np
import math

#startTime = 2433282.5
#startTime = 2396758.5
#startTime = 2459200.5
startTime = 2450800.5

#name ='105140'
#name ='Hidalgo'
#name = '2016RB1'
name = '367943'
#name = '397326'

#t2, i2, x2, y2, z2 = np.loadtxt("out%s_b.dat" % name, usecols=(0,1,3,4,5), skiprows=0, unpack=True)
#skr = int(t2[0] - startTime);
#t1, x1, y1, z1 = np.loadtxt("../TestSet/%s_d.dat" % name, usecols=(0,1,2,3), skiprows=skr, unpack=True)

t2, i2, x2, y2, z2 = np.loadtxt("out%s_h.dat" % name, usecols=(0,1,3,4,5), skiprows=0, unpack=True)
#skr = int(t2[0] - startTime);
skr = int(startTime - 2396758.5)
t1, x1, y1, z1 = np.loadtxt("../TestSet1850/%s_h.dat" % name, usecols=(0,1,2,3), skiprows=skr, unpack=True)

#ii = 27
ii = 0

tt2 = t2[i2==ii]
xx2 = x2[i2==ii]
yy2 = y2[i2==ii]
zz2 = z2[i2==ii]

print('A', t1[0], x1[0], y1[0], z1[0]);
print('B', tt2[0], xx2[0], yy2[0], zz2[0]);

for i in range(len(tt2)):
	dx = x1[i] - xx2[i]
	dy = y1[i] - yy2[i]
	dz = z1[i] - zz2[i]

	d = math.sqrt(dx * dx + dy * dy + dz * dz)

	print((t1[i] - startTime)/365.25, d * 149597870700.0)
