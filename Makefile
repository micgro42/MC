### Compile settings
CC=g++
CXXWARNINGS=-Wall -Wextra -pedantic #-Wshadow -Wlogical-op
CFLAGS= -Ofast -fopenmp -std=c++11 $(CXXWARNINGS) -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/

### Linker settings
LDFLAGS= -fopenmp
LDLIBS= -L/users/stud/micgro42/boost/lib -lpthread -lboost_log -lboost_unit_test_framework -lboost_thread -lboost_filesystem -lboost_date_time -lboost_chrono -lboost_system

### Debug Settings
DBGDIR = Debug
DBGCFLAGS = -g3 -DDEBUG

### Release Settings
RLSDIR = Release
RELCFLAGS = -DNDEBUG

SRCDIR=src/

all: unittest main

debug: CFLAGS += $(DBGCFLAGS)
debug: OUTDIR = $(DBGDIR)
debug: all

main: Mc.o main.o
	$(CC) $(LDFLAGS) Mc.o main.o $(LDLIBS) -o main

Mc.o: $(SRCDIR)Mc.cpp $(SRCDIR)Mc.h
	$(CC) $(CFLAGS) -c $(SRCDIR)Mc.cpp
	
main.o: $(SRCDIR)main.cpp 
	$(CC) $(CFLAGS) -c $(SRCDIR)main.cpp
	
unittest: Mc_test.o Mc.o
	$(CC) $(LDFLAGS) Mc.o Mc_test.o $(LDLIBS) -o unittest

Mc_test.o: $(SRCDIR)Mc_test.cpp Mc.o
	$(CC) $(CFLAGS) -c $(SRCDIR)Mc_test.cpp
	
geom_pbc.o: $(SRCDIR)geom_pbc.c $(SRCDIR)global.h
	 $(CC) $(CFLAGS) -c $(SRCDIR)geom_pbc.c
	 
stat5.o: $(SRCDIR)stat5.c $(SRCDIR)stat5.h
	$(CC) $(CFLAGS) -c $(SRCDIR)stat5.c

clean:
	rm -f *.o main unittest
	
.PHONY: all clean debug release 
