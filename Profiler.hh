#ifndef Profiler_hh
#define Profiler_hh 1

#include <map>
#include <iostream>
#include <iomanip>

using namespace std;

class Profiler
{
  public:
   static map<string, double> timer;
   static map<string, int> counter;
   static float start_time;

   Profiler();
   map<string,float> GetTimes();
   void PrintTimes();
   void PrintCounters();
   void PrintAll();


};





#endif
