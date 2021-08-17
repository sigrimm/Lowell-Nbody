import numpy as np

useHelioCentric = 1


if(useHelioCentric == 1):
	tS, xS, yS, zS, vxS, vyS, vzS = np.loadtxt('Sun_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t0, x0, y0, z0, vx0, vy0, vz0 = np.loadtxt('Mercury_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t1, x1, y1, z1, vx1, vy1, vz1 = np.loadtxt('Venus_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t2, x2, y2, z2, vx2, vy2, vz2 = np.loadtxt('Earth_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t3, x3, y3, z3, vx3, vy3, vz3 = np.loadtxt('Mars_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t4, x4, y4, z4, vx4, vy4, vz4 = np.loadtxt('Jupiter_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t5, x5, y5, z5, vx5, vy5, vz5 = np.loadtxt('Saturn_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t6, x6, y6, z6, vx6, vy6, vz6 = np.loadtxt('Uranus_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t7, x7, y7, z7, vx7, vy7, vz7 = np.loadtxt('Neptune_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t8, x8, y8, z8, vx8, vy8, vz8 = np.loadtxt('Pluto_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t9, x9, y9, z9, vx9, vy9, vz9 = np.loadtxt('Ceres_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t10, x10, y10, z10, vx10, vy10, vz10 = np.loadtxt('Pallas_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t11, x11, y11, z11, vx11, vy11, vz11 = np.loadtxt('Juno_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t12, x12, y12, z12, vx12, vy12, vz12 = np.loadtxt('Vesta_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t13, x13, y13, z13, vx13, vy13, vz13 = np.loadtxt('Hygiea_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t14, x14, y14, z14, vx14, vy14, vz14 = np.loadtxt('Eunomia_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t15, x15, y15, z15, vx15, vy15, vz15 = np.loadtxt('Euphrosyne_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t16, x16, y16, z16, vx16, vy16, vz16 = np.loadtxt('Europa_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t17, x17, y17, z17, vx17, vy17, vz17 = np.loadtxt('Davida_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t18, x18, y18, z18, vx18, vy18, vz18 = np.loadtxt('Interamnia_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t19, x19, y19, z19, vx19, vy19, vz19 = np.loadtxt('Moon_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t20, x20, y20, z20, vx20, vy20, vz20 = np.loadtxt('Psyche_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t21, x21, y21, z21, vx21, vy21, vz21 = np.loadtxt('Cybele_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t22, x22, y22, z22, vx22, vy22, vz22 = np.loadtxt('Thisbe_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	#t23, x23, y23, z23, vx23, vy23, vz23 = np.loadtxt('Doris_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t23, x23, y23, z23, vx23, vy23, vz23 = np.loadtxt('Iris_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	#t24, x24, y24, z24, vx24, vy24, vz24 = np.loadtxt('Patientia_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t24, x24, y24, z24, vx24, vy24, vz24 = np.loadtxt('Camilla_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t25, x25, y25, z25, vx25, vy25, vz25 = np.loadtxt('Sylvia_h.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
else:
	tS, xS, yS, zS, vxS, vyS, vzS = np.loadtxt('Sun_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t0, x0, y0, z0, vx0, vy0, vz0 = np.loadtxt('Mercury_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t1, x1, y1, z1, vx1, vy1, vz1 = np.loadtxt('Venus_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t2, x2, y2, z2, vx2, vy2, vz2 = np.loadtxt('Earth_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t3, x3, y3, z3, vx3, vy3, vz3 = np.loadtxt('Mars_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t4, x4, y4, z4, vx4, vy4, vz4 = np.loadtxt('Jupiter_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t5, x5, y5, z5, vx5, vy5, vz5 = np.loadtxt('Saturn_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t6, x6, y6, z6, vx6, vy6, vz6 = np.loadtxt('Uranus_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t7, x7, y7, z7, vx7, vy7, vz7 = np.loadtxt('Neptune_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t8, x8, y8, z8, vx8, vy8, vz8 = np.loadtxt('Pluto_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t9, x9, y9, z9, vx9, vy9, vz9 = np.loadtxt('Ceres_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t10, x10, y10, z10, vx10, vy10, vz10 = np.loadtxt('Pallas_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t11, x11, y11, z11, vx11, vy11, vz11 = np.loadtxt('Juno_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t12, x12, y12, z12, vx12, vy12, vz12 = np.loadtxt('Vesta_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t13, x13, y13, z13, vx13, vy13, vz13 = np.loadtxt('Hygiea_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t14, x14, y14, z14, vx14, vy14, vz14 = np.loadtxt('Eunomia_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t15, x15, y15, z15, vx15, vy15, vz15 = np.loadtxt('Euphrosyne_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t16, x16, y16, z16, vx16, vy16, vz16 = np.loadtxt('Europa_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t17, x17, y17, z17, vx17, vy17, vz17 = np.loadtxt('Davida_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t18, x18, y18, z18, vx18, vy18, vz18 = np.loadtxt('Interamnia_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t19, x19, y19, z19, vx19, vy19, vz19 = np.loadtxt('Moon_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t20, x20, y20, z20, vx20, vy20, vz20 = np.loadtxt('Psyche_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t21, x21, y21, z21, vx21, vy21, vz21 = np.loadtxt('Cybele_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t22, x22, y22, z22, vx22, vy22, vz22 = np.loadtxt('Thisbe_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	#t23, x23, y23, z23, vx23, vy23, vz23 = np.loadtxt('Doris_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t23, x23, y23, z23, vx23, vy23, vz23 = np.loadtxt('Iris_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	#t24, x24, y24, z24, vx24, vy24, vz24 = np.loadtxt('Patientia_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t24, x24, y24, z24, vx24, vy24, vz24 = np.loadtxt('Camilla_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)
	t25, x25, y25, z25, vx25, vy25, vz25 = np.loadtxt('Sylvia_b.dat', usecols=(0,1,2,3,4,5,6), unpack=True)

if(useHelioCentric == 1):
	f = open("All_h.dat", "w")
else:
	f = open("All_b.dat", "w")


for i in range(len(t0)):
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (tS[i], 0, xS[i], yS[i], zS[i], vxS[i], vyS[i], vzS[i]), file=f)		#Sun
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t0[i], 1, x0[i], y0[i], z0[i], vx0[i], vy0[i], vz0[i]), file=f)		#Mercury
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t1[i], 2, x1[i], y1[i], z1[i], vx1[i], vy1[i], vz1[i]), file=f)		#Venus
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t2[i], 3, x2[i], y2[i], z2[i], vx2[i], vy2[i], vz2[i]), file=f)		#Earth
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t3[i], 4, x3[i], y3[i], z3[i], vx3[i], vy3[i], vz3[i]), file=f)		#Mars
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t4[i], 5, x4[i], y4[i], z4[i], vx4[i], vy4[i], vz4[i]), file=f)		#Jupiter
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t5[i], 6, x5[i], y5[i], z5[i], vx5[i], vy5[i], vz5[i]), file=f)		#Saturn
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t6[i], 7, x6[i], y6[i], z6[i], vx6[i], vy6[i], vz6[i]), file=f)		#Uranus
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t7[i], 8, x7[i], y7[i], z7[i], vx7[i], vy7[i], vz7[i]), file=f)		#Neptune
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t8[i], 9, x8[i], y8[i], z8[i], vx8[i], vy8[i], vz8[i]), file=f)		#Pluto
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t19[i], 20, x19[i], y19[i], z19[i], vx19[i], vy19[i], vz19[i]), file=f)	#Moon
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t9[i], 10, x9[i], y9[i], z9[i], vx9[i], vy9[i], vz9[i]), file=f)	#Ceres
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t10[i], 11, x10[i], y10[i], z10[i], vx10[i], vy10[i], vz10[i]), file=f)	#Pallas
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t11[i], 12, x11[i], y11[i], z11[i], vx11[i], vy11[i], vz11[i]), file=f)	#Juno
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t12[i], 13, x12[i], y12[i], z12[i], vx12[i], vy12[i], vz12[i]), file=f)	#Vesta
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t23[i], 24, x23[i], y23[i], z23[i], vx23[i], vy23[i], vz23[i]), file=f)	#Iris
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t13[i], 14, x13[i], y13[i], z13[i], vx13[i], vy13[i], vz13[i]), file=f)	#Hygiea
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t14[i], 15, x14[i], y14[i], z14[i], vx14[i], vy14[i], vz14[i]), file=f)	#Eunomia
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t20[i], 21, x20[i], y20[i], z20[i], vx20[i], vy20[i], vz20[i]), file=f)	#Psyche
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t15[i], 16, x15[i], y15[i], z15[i], vx15[i], vy15[i], vz15[i]), file=f) #Euphrosyne
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t16[i], 17, x16[i], y16[i], z16[i], vx16[i], vy16[i], vz16[i]), file=f) #Europa
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t21[i], 22, x21[i], y21[i], z21[i], vx21[i], vy21[i], vz21[i]), file=f)	#Cybele
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t25[i], 26, x25[i], y25[i], z25[i], vx25[i], vy25[i], vz25[i]), file=f)	#Sylvia
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t22[i], 23, x22[i], y22[i], z22[i], vx22[i], vy22[i], vz22[i]), file=f)	#Thisbe
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t24[i], 25, x24[i], y24[i], z24[i], vx24[i], vy24[i], vz24[i]), file=f)	#Camilla
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t17[i], 18, x17[i], y17[i], z17[i], vx17[i], vy17[i], vz17[i]), file=f)	#Davida
	print("%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g" % (t18[i], 19, x18[i], y18[i], z18[i], vx18[i], vy18[i], vz18[i]), file=f) #Interamnia
f.close()
