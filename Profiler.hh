#ifndef Profiler_hh
#define Profiler_hh 1

#include <map>
#include <string>

class Profiler
{
  public:
   static std::map<std::string, double> timer;
   static std::map<std::string, int> counter;
   static float start_time;

   Profiler();
   std::map<std::string,float> GetTimes();
   void PrintTimes();
   void PrintCounters();
   void PrintAll();


};





#endif
