from astroquery.jplhorizons import Horizons

epoch_start = '1850-01-01'  # no earlier than 1950 due to Horizons limits
epoch_stop = '2026-12-31'


useHelioCentric = 1
###small bodies###
#single time step
#obj = Horizons(id='Ceres',epochs=2458133.33546,location='@0')


# multiple time steps
#barycentric:
#obj = Horizons(id='Ceres',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@0')

#heliocentric
#obj = Horizons(id='Ceres',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')

#obj = Horizons(id='Flora',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')

###large bodies###

#Earth 399
#Moon 301
#Earth-Moon system 3

#Moon :
#obj = Horizons(id_type='majorbody', id='301',epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')

ii = ([0, 1, 2, 399, 4, 5, 6, 7, 8, 9, 301, 10])
#name = (['Barycenter', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon', 'Sun', 'Ceres', 'Pallas', 'Juno', 'Vesta', 'Hygiea', 'Eunomia', 'Euphrosyne', 'Europa', 'Davida', 'Interamnia', 'Psyche', 'Cybele', 'Thisbe', 'Doris', 'Patientia', 'Sylvia'])
name = (['Barycenter', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon', 'Sun', 'Ceres', 'Pallas', 'Juno', 'Vesta', 'Hygiea', 'Eunomia', 'Euphrosyne', 'Europa', 'Davida', 'Interamnia', 'Psyche', 'Cybele', 'Thisbe', 'Iris', 'Camilla', 'Sylvia'])

for i in range(len(name)):

	if(i < 12):
		if(useHelioCentric == 1):
			obj = Horizons(id_type='majorbody', id=ii[i],epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')
		else:
			obj = Horizons(id_type='majorbody', id=ii[i],epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@0')
	else:
		if(useHelioCentric == 1):
			obj = Horizons(id=name[i],epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@sun')
		else:
			obj = Horizons(id=name[i],epochs={'start':epoch_start,'stop':epoch_stop,'step':'1d'},location='@0')

	#retrieve state vector
	#vec=obj.vectors(refplane='ecliptic', cache=False)
	vec=obj.vectors(refplane='earth', cache=False)

	#print(vec['targetname','datetime_jd','x','y','z','vx','vy','vz'])

	targetname = vec['targetname']
	datetime_jd = vec['datetime_jd']
	x = vec['x']
	y = vec['y']
	z = vec['z']
	vx = vec['vx']
	vy = vec['vy']
	vz = vec['vz']

	print("****** %s ********" % name[i])

	if(useHelioCentric == 1):
		f = open("%s_h.dat" % name[i], "w")
	else:
		f = open("%s_b.dat" % name[i], "w")

	for i in range(len(vec)):
		print("%.20g %.30g %.30g %.30g %.30g %.30g %.30g" % (datetime_jd[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]), file = f)
	f.close()
