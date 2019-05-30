main: main.o bagse_mixture.o GenEM_mixture.o  bagse_logistic.o GenEM_logistic.o logistic.o GenEM.o
	g++ -O3  main.o bagse_mixture.o GenEM_mixture.o  bagse_logistic.o GenEM_logistic.o logistic.o GenEM.o -lgsl -lgslcblas -lm -o bagse  
main.o: main.cc 
	g++ -c  main.cc
bagse_mixture.o: bagse_mixture.cc bagse_mixture.h
	g++    -c  bagse_mixture.cc 
GenEM_mixture.o: GenEM_mixture.cc GenEM_mixture.h
	g++ -c GenEM_mixture.cc 
bagse_logistic.o: bagse_logistic.cc bagse_logistic.h
	g++ -c bagse_logistic.cc
GenEM_logistic.o: GenEM_logistic.cc GenEM_logistic.h
	g++ -c GenEM_logistic.cc
logistic.o: logistic.cc logistic.h
	g++ -c logistic.cc
GenEM.o: GenEM.cc GenEM.h
	g++ -c GenEM.cc
clean:
	rm *.o bagse	
