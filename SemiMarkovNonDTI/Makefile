.PHONY: clean all debug

all: time_jumps.o main.o update.o OneSiteASEPCloning.exe

debug: OneSiteASEPCloning-g.exe

main.o: main.cpp type_defs.h
	g++ -O3 -Wall  -c  main.cpp 

update.o: update.cpp type_defs.h
	g++ -O3 -Wall  -c  update.cpp 

time_jumps.o: time_jumps.cpp type_defs.h
	g++ -O3 -Wall  -c time_jumps.cpp 

OneSiteASEPCloning.exe: time_jumps.o main.o update.o
	g++ -O3 -Wall  -o OneSiteASEPCloning.exe main.o  time_jumps.o   update.o -lm -lgsl -lgslcblas

OneSiteASEPCloning-g.exe:
	g++ -O0 -Wall -g -o OneSiteASEPCloning-g.exe  time_jumps.cpp  main.cpp update.cpp -lm -lgsl -lgslcblas

clean:
	rm -rf *exe
	rm -rf *.o

