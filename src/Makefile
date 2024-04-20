CC=g++
LDFLAGS=
DEPS = btf-count.h utility.h
SRC = main.cc btf-count.cc
PLATFORM := $(shell uname)
CPPFLAGS=-g -fopenmp

all: btfc

btfc: $(SRC) $(DEPS)
	$(CC) -std=c++11 $(CPPFLAGS) $(LDFLAGS) -O3 -pthread -o $@ $(SRC)
	
format:
	clang-format -i *.cpp *.h

clean:
	rm -rf *.o btfc