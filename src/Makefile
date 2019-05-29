main: bagse_main.o bagse.o GenEM_mixture.o 
	g++  -O3 bagse_main.o bagse.o GenEM_mixture.o  -lm -o bagse  
bagse_main.o: bagse_main.cc 
	g++ -c  bagse_main.cc
bagse.o: bagse.cc bagse.h
	g++    -c  bagse.cc 
GenEM.o: GenEM_mixture.cc GenEM_mixture.h
	g++ -c GenEM_mixture.cc 
clean:
	rm *.o bagse	
