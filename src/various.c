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

double n_over_d(ND the_nd){ // for positive n and d. return -1 if denom is zero.
  return (the_nd.n > 0)? (double)the_nd.n/(double)the_nd.d : -1;
}
/* double n_over_d(ND nd){ */
/*   return (nd.d > 0)? (double)nd.n/(double)nd.d : NAN; */
/* } */
