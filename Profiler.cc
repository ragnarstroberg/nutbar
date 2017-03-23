
#include "Profiler.hh"
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>


map<string, double> Profiler::timer;
map<string, int> Profiler::counter;
float Profiler::start_time = -1;

Profiler::Profiler()
{
  if (start_time < 0)
  {
    start_time = omp_get_wtime();
    counter["N_Threads"] = omp_get_max_threads();
  }
}



map<string,float> Profiler::GetTimes()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  map<string,float> times;
  times["user"] = ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec;
  times["system"] = ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
  times["real"] = omp_get_wtime() - start_time;
  return times;
}


void Profiler::PrintTimes()
{
  auto time_tot = GetTimes();
   
   cout << "====================== TIMES (s) ====================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : timer )
   {
     int nfill = (int) (20 * it.second / time_tot["real"]);
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(5) << std::right << it.second;
     cout << " (" << setw(4) << setprecision(1) << 100*it.second / time_tot["real"] << "%) |";
     for (int ifill=0; ifill<nfill; ifill++) cout << "*";
     for (int ifill=nfill; ifill<20; ifill++) cout << " ";
     cout << "|";
     cout  << endl;
     
   }

}

void Profiler::PrintCounters()
{
   cout << "===================== COUNTERS =====================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : counter )
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(0) << std::right << it.second  << endl;
}


void Profiler::PrintAll()
{
  PrintCounters();
  PrintTimes();
}
