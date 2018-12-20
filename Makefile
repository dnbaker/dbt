.PHONY=all tests clean obj update
CXX=g++
CC=gcc

GMATCH=$(findstring g++,$(CXX))

CLHASH_CHECKOUT = "&& git checkout master"
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -Wunused-variable -Wno-class-memaccess

ifndef EXTRA
	EXTRA:=
endif
ifndef INCPLUS
	INCPLUS:=
endif
ifndef EXTRA_LD
	EXTRA_LD:=
endif
DBG:=
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 \
	  -pipe -fno-strict-aliasing -march=native $(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++17 $(WARNINGS)
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS) $(INCLUDE)
LIB=-lz
LD=-L. $(EXTRA_LD)
INCLUDE += -Iinclude

ifneq (,$(findstring g++,$(CXX)))
	ifeq ($(shell uname),Darwin)
		ifeq (,$(findstring clang,$(CXX)))
			POPCNT_CXX:=clang
		else
			POPCNT_CXX:=$(CXX)
		endif
	else
		POPCNT_CXX:=$(CXX)
	endif
endif

EX=$(patsubst src/%.cpp,%,$(wildcard src/*.cpp))


all: pfparse

d: $(D_EX)

HEADERS=$(wildcard include/*.h) fpwrap/fpwrap.h

clhash.o:
	cd clhash && make $@ && cp $@ ..

SANITIZERS= # address undefined leak # pointer-compare # pointer-substract

%: src/%.cpp $(HEADERS) clhash.o
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o  $(OBJ) $(patsubst %,-fsanitize=%,$(SANITIZERS)) $< -o $@ $(LIB)
mpitest: src/mpitest.cpp $(HEADERS)
	mpiCC src/mpitest.cpp -lstdc++ -o mpitest
clean:
	rm -f $(EX)
mostlyclean: clean
