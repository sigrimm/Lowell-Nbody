
integrator: main.o asteroid.o
	g++ -g -o integrator main.o asteroid.o

main.o: main.cpp define.h force.h integrator.h Chebyshev.h
	g++ -g -c -O3 -o main.o main.cpp

asteroid.o: asteroid.cpp asteroid.h define.h
	g++ -g -c -O3 -o asteroid.o asteroid.cpp

