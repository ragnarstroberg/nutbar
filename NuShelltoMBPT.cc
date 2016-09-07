#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"
#include "TransitionDensity.hh"


using namespace std;






int main(int argc, char** argv)
{

TransitionDensity trans;

if (argc>1)
{
  trans.ReadInputFromFile( argv[1]);
}
else
{
  trans.ReadInputInteractive();
}


trans.ReadFiles();

trans.CalculateMschemeAmplitudes();

trans.WriteEGV("mbpt.egv");
trans.WriteTRDENS_input("trdens.in");

  return 0;
}

