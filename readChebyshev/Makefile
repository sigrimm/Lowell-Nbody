readChebyshev: readChebyshev.o planets.o perturbers.o
	g++ planets.o perturbers.o readChebyshev.o -o readChebyshev

planets.o: planets.cpp planets.h define.h
	g++ -c -o planets.o planets.cpp

perturbers.o: perturbers.cpp perturbers.h define.h
	g++ -c -o perturbers.o perturbers.cpp

readChebyshev.o: readChebyshev.cpp planets.h perturbers.h define.h
	g++ -c -o readChebyshev.o readChebyshev.cpp

