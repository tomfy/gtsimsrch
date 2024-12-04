#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
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
  return (the_nd.d > 0)? (double)the_nd.n/(double)the_nd.d : NAN;
}
/* double n_over_d(ND nd){ */
/*   return (nd.d > 0)? (double)nd.n/(double)nd.d : NAN; */
/* } */

void print_d_r(FILE* fh, ND nd){
  if(nd.d == 0){
    fprintf(fh, "    0     -  ");
  }else{
    fprintf(fh, "%5ld %6.5lf  ", nd.d, (double)nd.n/(double)nd.d);
  }
}

void print_n_over_d(FILE* fh, ND nd){
  if(nd.d == 0){
    fprintf(fh, " -  ");
  }else{
    fprintf(fh, "%6.5lf  ", (double)nd.n/(double)nd.d);
  }
}

bool NDs_equal(ND nd1, ND nd2){
  return ((nd1.n == nd2.n) && (nd1.d == nd2.d));
}
