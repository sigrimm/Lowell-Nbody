import numpy as np
import math
import os


os.system("cp initial.dat initial_1.dat")


S = 16

n = np.zeros(S)
t1 = np.zeros(S)
t2 = np.zeros(S)
t3 = np.zeros(S)


N = 1
n[0] = N

for i in range(S-1):
	#print(2*N)

	os.system("cat initial_%d.dat initial_%d.dat > initial_%d.dat" % (N, N, 2*N))
	N *= 2
	n[i+1] = int(N)


N = 1
for i in range(S - 3):
	os.system("./integrator -in initial_%d.dat -outInterval 1000000 > test" % N)

	with open('test', "r") as f:
		lines = f.readlines()

		#last line
		time = lines[-1]

		time = time.split()

		print(int(n[i]), time[-1])

		t1[i] = time[-1]

	N *= 2

N = 1
for i in range(S - 1):
	os.system("./integratorGPU -in initial_%d.dat -mode 1 -outInterval 1000000 > test" % N)

	with open('test', "r") as f:
		lines = f.readlines()

		#last line
		time = lines[-1]

		time = time.split()

		print(int(n[i]), time[-1])

		t2[i] = time[-1]

	N *= 2
N = 1
for i in range(S):
	os.system("./integratorGPU -in initial_%d.dat -mode 0 -outInterval 1000000 > test" % N)

	with open('test', "r") as f:
		lines = f.readlines()

		#last line
		time = lines[-1]

		time = time.split()

		print(int(n[i]), time[-1])

		t3[i] = time[-1]

	N *= 2


outfile = 'Time.dat'

with open(outfile, "w") as f:
	print("#N	CPU	GPU	GPU2", file = f)
	for i in range(S):
		print("%d\t%g\t%g\t%g" % (n[i], t1[i], t2[i], t3[i]), file = f)

print("timings are written into file", outfile)

#Run time in seconds: 
