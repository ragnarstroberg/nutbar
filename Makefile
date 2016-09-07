
CPP = g++

FLAGS = -std=c++11
INCLUDE = 
LIBS = -lgsl -lopenblas

ifeq ($(DEBUG),on)
  FLAGS += -g
endif

OBJ = NuVec.o NuBasis.o NuProj.o JMState.o JBasis.o TransitionDensity.o
EXE = NuShelltoMBPT

all: $(EXE)

%.o : %.cc
	$(CPP) -c $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

NuShelltoMBPT : NuShelltoMBPT.cc $(OBJ)
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	rm -f $(EXE) $(OBJ)
