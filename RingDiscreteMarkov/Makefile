.PHONY: clean all debug

all:  main.o update.o Ring.exe

debug: Ring-g.exe

update.o : update.cpp
	g++ -O3 -Wall  -c  update.cpp 

main.o: main.cpp
	g++ -O3 -Wall  -c  main.cpp 


Ring.exe: main.o update.o
	g++ -O3 -Wall  -o Ring.exe  main.o update.o -lm -lgsl -lgslcblas


Ring-g.exe:  main.cpp update.cpp
	g++ -O0 -Wall -g -o Ring-g.exe   main.cpp update.cpp -lm -lgsl -lgslcblas



clean: 
	rm -rf *exe
	rm -rf *.o

