all: integratorGPU integrator

integrator: main.o asteroid.o
	g++ -g -o integrator main.o asteroid.o

integratorGPU: mainGPU.o asteroid1.o asteroidGPU.o
	nvcc -o integratorGPU mainGPU.o asteroid1.o asteroidGPU.o

main.o: main.cpp define.h RKF.h force.h integrator.h Chebyshev.h
	g++ -g -c -O3 -DUSEGPU=0 -o main.o main.cpp

asteroid.o: asteroid.cpp asteroid.h define.h
	g++ -g -c -O3 -DUSEGPU=0 -o asteroid.o asteroid.cpp


asteroid1.o: asteroid.cpp asteroid.h define.h
	nvcc -g -c -O3 -DUSEGPU=1 -o asteroid1.o asteroid.cpp

mainGPU.o: main.cpp define.h RKF.h force.h integrator.h Chebyshev.h
	nvcc -g -c -O3 -DUSEGPU=1 -o mainGPU.o main.cpp

asteroidGPU.o: asteroid.cu asteroid.h define.h
	nvcc -c -O3 -o asteroidGPU.o asteroid.cu

clean:
	rm *.o
	rm integrator
	rm integratorGPU
