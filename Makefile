
CPP = g++

FLAGS = -std=c++11
INCLUDE = 
LIBS = -lgsl -lopenblas

ifeq ($(DEBUG),on)
  FLAGS += -g
endif

OBJ = NuVec.o NuBasis.o NuProj.o JMState.o JBasis.o
EXE = NuShelltoMBPT

all: $(EXE)

%.o : %.cc
	$(CPP) -c $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

NuShelltoMBPT : NuSHelltoMBPT.cc $(OBJ)
	$(CPP) $(FLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	@rm $(EXE) $(OBJ)
