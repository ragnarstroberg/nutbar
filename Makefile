
CPP = g++

FLAGS = -std=c++11 -fopenmp -O3 -march=native
INCLUDE = -I$(HOME)/include/armadillo
LIBS = -lgsl -lopenblas

# This doesn't work everywhere yet.
#LIBS += -L/opt/boost/1.58.0/lib  -lboost_system -lboost_filesystem 
LIBS +=   -lboost_system -lboost_filesystem 
FLAGS += -DNOBOOST


ifeq ($(DEBUG),on)
  FLAGS += -g
endif

OBJ = NuVec.o NuBasis.o NuProj.o JMState.o JBasis.o TransitionDensity.o
EXE = nutbar
#EXE = nutbar_test

all: $(EXE)

%.o : %.cc
	$(CPP) -c $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

nutbar : nutbar.cc $(OBJ)
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

nutbar_test : nutbar.cc $(OBJ)
	@echo TESTING PURPOSES ONLY
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	rm -f $(EXE) $(OBJ)
