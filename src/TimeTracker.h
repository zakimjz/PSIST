/*
* TimeTracker.h 
*
*/


#include <sys/time.h>
#include <unistd.h>

#define microsec 1000000.0

class TimeTracker {
 private:
  struct timeval start_time;
  struct timeval stop_time;
  bool  running;
 public:
  TimeTracker() {
    running=false;
  }

  void start() {
    gettimeofday(&start_time, (struct timezone *)0);
    running=true;
  }

  void stop(){
    if (!running) return;
    else{
      gettimeofday(&stop_time, (struct timezone *)0);      
      running=false;
    }
  }
  

  double getTotalTime(){
    double st, en;
    //if (!running) return(-1.0);
    //else {
    st = start_time.tv_sec + (start_time.tv_usec/microsec);
    en = stop_time.tv_sec + (stop_time.tv_usec/microsec);
    //running=false;
    return (double)(en - st);
    //}
  }
};
