#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "various.h"




double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

double clock_ticks_to_seconds(clock_t nticks){
  return (double)nticks/(double)CLOCKS_PER_SEC;
}
