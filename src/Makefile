CC=g++
LDFLAGS=
DEPS = 
SRC = main.cc
PLATFORM := $(shell uname)
CPPFLAGS=-g -fopenmp

all: btfc

btfc: $(SRC) $(DEPS)
	$(CC) -std=c++11 $(CPPFLAGS) $(LDFLAGS) -O3 -pthread -o $@ $(SRC)
	
format:
	clang-format -i *.cpp *.h

clean :
	rm -rf *.o btfc