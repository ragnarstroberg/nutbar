
CPP = g++

FLAGS = -std=c++11 -fopenmp -O3 -march=native -Wall
#INCLUDE = -I$(HOME)/include/armadillo
INCLUDE = -Iarmadillo
LIBS = -lgsl -lopenblas

# This doesn't work everywhere yet.
#LIBS += -L/opt/boost/1.58.0/lib  -lboost_system -lboost_filesystem 
LIBS +=   -lboost_system -lboost_filesystem 
FLAGS += -DNOBOOST


# I assume we're running on linux
OS = LINUX
# But in case we're crazy enough to run on MacOS, might as well check...
ifneq (,$(findstring arwin,$(shell uname))) 
  OS = MACOS
endif



ifeq ($(OS),MACOS)
#  FLAGS     = -Xpreprocessor -fopenmp -O3 -march=native -std=c++11 -fPIC 
  FLAGS     = -Xpreprocessor -fopenmp -O3  -std=c++11  -Wall
  LIBS += -lomp
endif



ifeq ($(DEBUG),on)
  FLAGS += -g -DVERBOSE
endif

ifeq ($(VERBOSE),on)
  FLAGS += -DVERBOSE
endif

OBJ = NuVec.o NuBasis.o NuProj.o JMState.o JBasis.o TransitionDensity.o Profiler.o ReadWrite.o Operators.o
EXE = nutbar
#EXE = nutbar_test

all: $(EXE)

%.o : %.cc
	$(CPP) -c $(FLAGS) $(INCLUDE) $^ -o $@ 

nutbar : nutbar.cc $(OBJ)
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

nutbar_test : nutbar.cc $(OBJ)
	@echo TESTING PURPOSES ONLY
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	rm -f $(EXE) $(OBJ)
