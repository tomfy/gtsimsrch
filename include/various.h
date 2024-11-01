#include <time.h>
#include <stdbool.h>

#define ID_NA_INDEX -1 //

typedef struct{
  long n; // numerator
  long d; // denominator
}ND; // numerator and denominator

typedef struct{
  long l1;
  long l2;
}two_longs;

typedef struct{
  long l1;
  long l2;
  long l3;
} three_longs;

typedef struct{
  long l1;
  long l2;
  long l3;
  long l4;
}four_longs;

typedef struct{
  double x1;
  double x2;
}two_doubles;

typedef struct{
  char ch1;
  char ch2;
}two_chars;

typedef struct{
  long XFmin;
  long XFmax;
  long NFhet;
  long XMmin;
  long XMmax;
  long NMhet;
  long XFmin_triple;
  long XFmax_triple;
  long XMmin_triple;
  long XMmax_triple;
}Xover_info;

//double n_over_d(ND nd);
double hi_res_time(void);
double clock_ticks_to_seconds(clock_t nticks);
void chomp(char* str); // remove any trailing newlines from str
double n_over_d(ND the_nd);
void print_d_r(FILE* fh, ND nd);
void print_n_over_d(FILE* fh, ND nd);
bool NDs_equal(ND nd1, ND nd2);
