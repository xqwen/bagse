main: main.o controller.o logistic.o 
	g++  -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas  -o bagse  
static: main.o controller.o logistic.o
	g++  -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas -static -o bagse.static
main.o: main.cc 
	g++ -c  main.cc	
controller.o: controller.cc classdef.h
	g++    -c  controller.cc 
logistic.o: logistic.cc logistic.h
	g++ -c logistic.cc -fpermissive
clean:
	rm *.o bagse	
