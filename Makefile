
CPP = g++

FLAGS = -std=c++11 -fopenmp -O3 -march=native
INCLUDE = -I$(HOME)/include/armadillo
LIBS = -lgsl -lopenblas

# This doesn't work everywhere yet.
#LIBS += -lboost_system -lboost_filesystem
FLAGS += -DNOBOOST


ifeq ($(DEBUG),on)
  FLAGS += -g
endif

OBJ = NuVec.o NuBasis.o NuProj.o JMState.o JBasis.o TransitionDensity.o
EXE = nutbar

all: $(EXE)

%.o : %.cc
	$(CPP) -c $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

nutbar : nutbar.cc $(OBJ)
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	rm -f $(EXE) $(OBJ)
