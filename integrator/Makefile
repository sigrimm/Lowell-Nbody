all: integrator integratorGPU

integrator: main.o asteroid.o RKF.o integrator.o read.o
	g++ -g -o integrator main.o asteroid.o RKF.o integrator.o read.o

main.o: main.cpp asteroid.h 
	g++ -g -c -O3 -DUSEGPU=0 -o main.o main.cpp

integrator.o: integrator.cpp asteroid.h Chebyshev.h force.h
	g++ -g -c -O3 -o integrator.o integrator.cpp

RKF.o: RKF.cpp asteroid.h
	g++ -g -c -O3 -o RKF.o RKF.cpp

asteroid.o: asteroid.cpp asteroid.h
	g++ -g -c -O3 -DUSEGPU=0 -o asteroid.o asteroid.cpp

read.o: read.cpp asteroid.h
	g++ -g -c -O3 -o read.o read.cpp



integratorGPU: mainGPU.o asteroid1.o asteroidGPU.o RKF.o integratorGPU.o read.o
	nvcc -o integratorGPU mainGPU.o asteroid1.o asteroidGPU.o RKF.o integratorGPU.o read.o

asteroid1.o: asteroid.cpp asteroid.h
	nvcc -g -c -O3 -DUSEGPU=1 -o asteroid1.o asteroid.cpp

mainGPU.o: main.cpp asteroid.h
	nvcc -g -c -O3 -DUSEGPU=1 -o mainGPU.o main.cpp

integratorGPU.o: integrator.cu asteroid.h ChebyshevGPU.h forceGPU.h
	nvcc -g -c -O3 -o integratorGPU.o integrator.cu

asteroidGPU.o: asteroid.cu asteroid.h
	nvcc -c -O3 -o asteroidGPU.o asteroid.cu

clean:
	rm *.o
	rm integrator
	rm integratorGPU
