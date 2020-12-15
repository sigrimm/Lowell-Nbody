from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB
import sys

epoch_start = '1950-01-01'  # no earlier than 1950 due to Horizons limits
epoch_stop = '1960-12-31'


useHelioCentric = 1
###small bodies###
#single time step
#obj = Horizons(id='Ceres',epochs=2458133.33546,location='@0')


# multiple time steps
#barycentric:
#obj = Horizons(id='Ceres',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@0')

#heliocentric
#obj = Horizons(id='Ceres',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')

name = str(sys.argv[1])
print(name)


if(useHelioCentric == 1):
	obj = Horizons(id=name,epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')
else:
	obj = Horizons(id=name,epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@0')

#retrieve state vector
vec=obj.vectors(refplane='earth', cache=False)
#vec=obj.vectors(refplane='ecliptic', cache=False)

#print(vec['targetname','datetime_jd','x','y','z','vx','vy','vz'])

targetname = vec['targetname']
datetime_jd = vec['datetime_jd']
x = vec['x']
y = vec['y']
z = vec['z']
vx = vec['vx']
vy = vec['vy']
vz = vec['vz']


## NON Gravitational constants ##
a = SBDB.query(name, full_precision=True)
#print(SBDB.schematic(a['orbit']))

A1 = 0.0
A2 = 0.0
A3 = 0.0
ALN = 0.1112620426
NK = 4.6142
NM = 2.15
NN = 5.093
R0 = 2.808

try:
	A = str(a['orbit']['model_pars']['A1'])
	
except:
	A = '0.0 AU / d2'
AA1 = A.split(' AU / d2')
if(len(AA1) != 2):
	print("wrong unit in A1")
	quit()
A1 = float(AA1[0])

try:
	A = str(a['orbit']['model_pars']['A2'])
except:
	A = '0.0 AU / d2'
AA2 = A.split(' AU / d2')
if(len(AA2) != 2):
	print("wrong unit in A1")
	quit()
A2 = float(AA2[0])

try:
	A = str(a['orbit']['model_pars']['A3'])
except:
	A = '0.0 AU / d2'
AA3 = A.split(' AU / d2')
if(len(AA3) != 2):
	print("wrong unit in A1")
	quit()
A13 = float(AA3[0])

try:
	ALN = float(a['orbit']['model_pars']['ALN'])
except:
	ALN = 0.1112620426
try:
	NK = float(a['orbit']['model_pars']['NK'])
except:
	NK = 4.6142
try:
	NM = float(a['orbit']['model_pars']['NM'])
except:
	NM = 2.15
try:
	NN = float(a['orbit']['model_pars']['NN'])
except:
	NN = 5.093
try:
	R0 = a['orbit']['model_pars']['R0']
	if(R0 == 'AU'):
		R0 = 1.0
	else:
		R0 = float(R0)
except:
	R0 = 2.808


print("A1 %s" % A1)
print("A2 %s" % A2)
print("A3 %s" % A3)
print("ALN %s" % ALN)
print("NK %s" % NK)
print("NM %s" % NM)
print("NN %s" % NN)
print("R0 %s" % R0)

if(useHelioCentric == 1):
	f = open("%s_h.dat" % name, "w")
else:
	f = open("%s_d.dat" % name, "w")

for i in range(len(vec)):
	print("%.20g %.30g %.30g %.30g %.30g %.30g %.30g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g" % (datetime_jd[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], A1, A2, A3, ALN, NK, NM, NN, R0), file = f)
f.close()

