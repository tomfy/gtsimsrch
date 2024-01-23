#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "various.h"


double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

double clock_ticks_to_seconds(clock_t nticks){
  return (double)nticks/(double)CLOCKS_PER_SEC;
}

void chomp(char* str){ // remove any trailing newlines from str
  long len = strlen(str);
  while(str[len-1] == '\n'){
    str[len-1] = '\0';
    len--;
  }
}
