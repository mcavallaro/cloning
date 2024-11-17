.PHONY: clean all debug

all: functions.o main.o update.o memTASEPCloning.exe time_increments_f.o memTASEPCloningDensity.exe

debug: memTASEPCloning-g.exe

main.o: main.cpp type_defs.h functions.h
	g++ -O3 -Wall  -c  main.cpp 

main-density.o: main-density.cpp type_defs.h functions.h
	g++ -O3 -Wall  -c  main-density.cpp


update.o: update.cpp type_defs.h functions.h
	g++ -O3 -Wall  -c  update.cpp 

functions.o: functions.c type_defs.h functions.h
	g++ -O3 -Wall  -c functions.c  

time_increments_f.o: time_increments_f.cpp type_defs.h functions.h
	g++ -O3 -Wall  -c time_increments_f.cpp  


memTASEPCloning.exe: functions.o main.o update.o time_increments_f.o
	g++ -O3 -Wall  -o memTASEPCloning.exe  functions.o  main.o update.o time_increments_f.o -lm -lgsl -lgslcblas

memTASEPCloningDensity.exe: functions.o main-density.o update.o time_increments_f.o
	g++ -O3 -Wall  -o memTASEPCloningDensity.exe  functions.o  main-density.o update.o time_increments_f.o -lm -lgsl -lgslcblas


memTASEPCloning-g.exe: *.h *cpp *.c
	g++ -O0 -Wall -g -o memTASEPCloning-g.exe  functions.c  time_increments_f.cpp main.cpp update.cpp -lm -lgsl -lgslcblas

clean:
	rm -rf *exe
	rm -rf *.o

