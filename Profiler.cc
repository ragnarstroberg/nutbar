
#include "Profiler.hh"
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>


std::map<std::string, double> Profiler::timer;
std::map<std::string, int> Profiler::counter;
float Profiler::start_time = -1;

Profiler::Profiler()
{
  if (start_time < 0)
  {
    start_time = omp_get_wtime();
    counter["N_Threads"] = omp_get_max_threads();
  }
}



std::map<std::string,float> Profiler::GetTimes()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  std::map<std::string,float> times;
  times["user"] = ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec;
  times["system"] = ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
  times["real"] = omp_get_wtime() - start_time;
  return times;
}


void Profiler::PrintTimes()
{
  auto time_tot = GetTimes();
   
   std::cout << "====================== TIMES (s) ====================" << std::endl;
   std::cout.setf(std::ios::fixed);
   for ( auto it : timer )
   {
     int nfill = (int) (20 * it.second / time_tot["real"]);
     std::cout << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(5) << std::right << it.second;
     std::cout << " (" << std::setw(4) << std::setprecision(1) << 100*it.second / time_tot["real"] << "%) |";
     for (int ifill=0; ifill<nfill; ifill++) std::cout << "*";
     for (int ifill=nfill; ifill<20; ifill++) std::cout << " ";
     std::cout << "|";
     std::cout  << std::endl;
     
   }

}

void Profiler::PrintCounters()
{
   std::cout << "===================== COUNTERS =====================" << std::endl;
   std::cout.setf(std::ios::fixed);
   for ( auto it : counter )
     std::cout << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(0) << std::right << it.second  << std::endl;
}


void Profiler::PrintAll()
{
  PrintCounters();
  PrintTimes();
}
