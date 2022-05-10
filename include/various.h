#include <time.h>

typedef struct{
  long n; // numerator
  long d; // denominator
}ND; // numerator and denominator

typedef struct{
  long l1;
  long l2;
  long l3;
  long l4;
}four_longs;

double hi_res_time(void);
double clock_ticks_to_seconds(clock_t nticks);
