
read: read.o asteroid.o
	g++ -g -o read read.o asteroid.o

read.o: read2.cpp define.h force.h integrator.h Chebyshev.h
	g++ -g -c -O3 -o read.o read2.cpp

asteroid.o: asteroid.cpp asteroid.h define.h
	g++ -g -c -O3 -o asteroid.o asteroid.cpp

