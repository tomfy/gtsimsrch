#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // needed for getopt
#include <sys/sysinfo.h> // needed for get_nprocs
#include <pthread.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>

// #include <math.h>
// #include <ctype.h>
// #include <limits.h>

#include "vect.h"
#include "various.h"

#define INIT_N_ACCESSIONS 10000
#define INIT_N_MARKERS 10000

// There will be one of these structs for each thread,
// each thread which will process the markers in the range from first_marker to last_marker
typedef struct{
  long n_accessions;
  Vstr* marker_lines; 
  long first_marker;
  long last_marker;
  bool use_alt_marker_id;
  double minGP;
 
  Vstr* marker_ids;
  Vchar** genotypes; // genotypes[i]->a[2*j+1] is genotype of ith accession, jth marker
} TD; // thread data

void* process_marker_range(void* x);

char token_to_genotype(char* token, long gtidx, long gpidx, double minGP);
char GTstr_to_dosage(char* tkn);
void get_GT_GQ_GP_indices(char* format, long* GTidx, long* GQidx, long* GPidx);
bool GP_to_quality_ok(char* token, double minGP);

double clock_time(clockid_t the_clock){
  struct timespec tspec;
  clock_gettime(the_clock, &tspec);
  return (double)(tspec.tv_sec + 1.0e-9*tspec.tv_nsec);
}

extern int errno; 

int main(int argc, char *argv[]){
  errno = 0;

  char* input_filename = NULL;
  FILE *in_stream = NULL;
  Vchar* output_filename = construct_vchar_from_str("vcftogts");
  FILE* out_stream = NULL; 

  // double minGQ = 0; // not implemented, but should be.
  double minGP = 0;

  long nprocs = (long)get_nprocs(); // returns 2*number of cores if hyperthreading.
  long Nthreads = (nprocs > 2)? nprocs/2 : 1; // default number of threads

  bool use_alt_marker_id = false; // default is to use marker ids in col 3 of vcf file.
  // (but if these are absent you can construct and use alternative marker ids from cols 1 and 2.)
  bool shuffle_accessions = false;
  long rand_seed = (unsigned)time(0);
     
 
  while (1) {
    int an_option;
    int option_index = 0;
    static struct option long_options[] = {
      {"input",   required_argument, 0,  'i'}, // vcf filename
      {"output",  required_argument, 0,  'o'}, // output filename
      {"pmin",  required_argument,  0,  'p'}, // min. 'estimated genotype probability'
      {"threads", required_argument, 0,  't'}, // number of threads to use
      {"alternate_marker_ids",  no_argument, 0, 'a'}, // construct marker ids from cols 1 and 2 (in case garbage in col 3)
      {"randomize",    no_argument, 0,  'r' }, // shuffle the order of the accessions in output
      {"seed", required_argument, 0, 's'}, // rng seed. Only relevant if shuffling.
      {0,         0,                 0,  0 }
    };
    if(an_option == -1) break;
    an_option = getopt_long_only(argc, argv, "", long_options, &option_index);
    //   fprintf(stderr, "# getopt_long retval: %d %c \n", c, c);
    switch(an_option){
    case 'i':
      input_filename = optarg;
      in_stream = fopen(input_filename, "r");
      if(in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'o':
      free_vchar(output_filename);
      output_filename = construct_vchar_from_str(optarg);
      break;
    case 'p' :
      if(sscanf(optarg, "%lf ", &minGP) != 1  ||  errno != 0){
	fprintf(stderr, "# minGP; conversion of argument %s to double failed.\n", optarg);
	exit(EXIT_FAILURE);
      }else if(minGP > 1){
	fprintf(stderr, "# minGP was set to %8.4lf , must be <= 1\n", minGP);
	exit(EXIT_FAILURE);
      }
      fprintf(stderr, "# minGP set to: %8.5lf\n", minGP);
      break;
    case 't' :
      if(sscanf(optarg, "%ld", &Nthreads) != 1  ||  errno != 0){
	fprintf(stderr, "# Nthreads; conversion of argument %s to long failed.\n", optarg);
	exit(EXIT_FAILURE);
      }else if(Nthreads > nprocs){
	fprintf(stderr, "# Setting Nthreads to max allowed value of %ld.\n", nprocs);
	Nthreads = nprocs;
      }else if(Nthreads < 0){
	fprintf(stderr, "# Setting Nthreads to min allowed value of 1.\n");
	Nthreads = 1;
      }
      break;
    case 'a' :
      use_alt_marker_id = true;
      break;
    case 'r' :
      shuffle_accessions = true;
      break;
  case 's' :
    if(sscanf(optarg, "%ld", &rand_seed) != 1  ||  errno != 0){
      fprintf(stderr, "# rand_seed; conversion of argument %s to long failed.\n", optarg);
      exit(EXIT_FAILURE);
    }
    break;
    } // end of command line processing loop
  }

  out_stream = fopen(output_filename->a, "w");
  if(out_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename->a);
    exit(EXIT_FAILURE);
  }

  clockid_t the_clock = CLOCK_MONOTONIC;
  struct timespec tspec;
  if(clock_getres(the_clock, &tspec) == 0){
    double t_resolution = tspec.tv_sec + 1.0e-9 * tspec.tv_nsec;
    if(t_resolution > 1.0e-3) fprintf(stderr, "# timing resolution is %8.5lf\n", t_resolution);
  }else{
    exit(EXIT_FAILURE);
  }
  double t0 = clock_time(the_clock);
  srand(rand_seed);
  
  // ****************************************************
  // *****   Read first line; store accession ids.  *****
  // ****************************************************

  char* line = NULL;
  size_t len = 0;
  ssize_t nread;
  
  long accid_count = 0;
  Vstr* accession_ids = construct_vstr(INIT_N_ACCESSIONS);
  while((nread = getline(&line, &len, in_stream)) != -1){
    char* saveptr = NULL;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    
    if((token == NULL) || (token[0] == '#' && token[1] == '#')) continue; // skip comments (starting with ##) and any empty lines
    if(token[0] == '#'){ // this the line with accession ids
      for(long ii = 1; ii <= 8; ii++){ // read in cols 1 through 8 "POS ID REF ..." but don't store them.
	token = strtok_r(NULL, "\t \n\r", &saveptr);
      }
      while(1){
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	if(token == NULL) break;
	char* acc_id = (char*)malloc((strlen(token)+1)*sizeof(char));
	push_to_vstr(accession_ids, strcpy(acc_id, token)); // store
	accid_count++;
      }
      break;
    }else{ 
      fprintf(stderr, "token: %s (should start with #)\n", token);
      exit(EXIT_FAILURE);
    }
  }
  long n_accessions = accid_count;

  // *********************************************************
  // *****   Read the rest of the lines, one per marker  *****
  // *****   Each line has genotypes for all accessions  *****
  // *********************************************************
 
  Vstr* marker_lines = construct_vstr(INIT_N_MARKERS);
  while((nread = getline(&line, &len, in_stream)) != -1){
    char* line_copy = (char*)malloc((nread+1)*sizeof(char));
    push_to_vstr(marker_lines, strcpy(line_copy, line));
  }
  long n_markers = marker_lines->size;
  free(line);
  fclose(in_stream);

  double t1 = clock_time(the_clock);
  fprintf(stderr, "# input done\n");
   
  // ********************************************************
  // *****  Extract genotypes, and quality information  *****
  // *****  Filter if requested and store genotypes     *****
  // ********************************************************
  
  Vstr* marker_ids = construct_vstr_empties(n_markers);

  char* str_of_spaces = (char*)malloc((2*n_markers + 1)*sizeof(char));
  for(long i=0; i<2*n_markers; i++){
    str_of_spaces[i] = ' ';
  }
  str_of_spaces[2*n_markers] = '\0'; // now a null terminated string, 2*n_markers spaces, followed by null.    
  Vchar** genotypes = (Vchar**)malloc(n_accessions*sizeof(Vchar*));
  for(long i=0; i<n_accessions; i++){
    genotypes[i] = construct_vchar_from_str(str_of_spaces);
  }
  free(str_of_spaces);

  TD* td;   
    
  if(Nthreads == 0){ // process without creating any new threads
    td = (TD*)malloc(sizeof(TD));
    fprintf(stderr, "# sizeof TD: %ld \n", (long)sizeof(TD));
    td->n_accessions = n_accessions;
    td->marker_lines = marker_lines;
    td->first_marker = 0;
    td->last_marker = n_markers-1;
    td->use_alt_marker_id = use_alt_marker_id;
    td->minGP = minGP;
    td->marker_ids = marker_ids;
    td->genotypes = genotypes;

    process_marker_range((void*)td);
    free(td);
 
  }else{ // 1 or more pthreads   
    td = (TD*)malloc(Nthreads*sizeof(TD));
    for(long i_thread = 0; i_thread<Nthreads; i_thread++){
      td[i_thread].n_accessions = n_accessions;
      td[i_thread].marker_lines = marker_lines;
      td[i_thread].first_marker = (i_thread == 0)? 0 : td[i_thread-1].last_marker + 1;
      td[i_thread].last_marker = (long)((double)(i_thread+1)*n_markers/Nthreads - 1);
      td[i_thread].use_alt_marker_id = use_alt_marker_id;
      td[i_thread].minGP = minGP;
      td[i_thread].marker_ids = marker_ids;
      td[i_thread].genotypes = genotypes;
    }
    td[Nthreads-1].last_marker = n_markers-1;
    
    pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
    for(long i=0; i<Nthreads; i++){
      int iret = pthread_create( thrids+i, NULL, process_marker_range, (void*) (td+i));
      if(iret > 0) fprintf(stderr, "# warning. pthread_create returned non-zero value. Thread %ld \n", (long)thrids[i]);
    }

    for(long i_thread=0; i_thread<Nthreads; i_thread++){
      pthread_join(thrids[i_thread], NULL);
    }
    free(td);
    free(thrids);
  }

  double t2 = clock_time(the_clock);
  fprintf(stderr, "# time for threads to run %8.4lf\n", t2 - t1);


  // **********************
  // *****  output  *******
  // **********************
  Vlong* accession_indices = construct_vlong_whole_numbers(accession_ids->size);
  if(shuffle_accessions) shuffle_vlong(accession_indices); 
    
  fprintf(out_stream, "MARKER");
  for(long i_marker=0; i_marker<marker_ids->size; i_marker++){
    fprintf(out_stream, " %s", marker_ids->a[i_marker]);
  }fprintf(out_stream, "\n");
  for(long i=0; i<accid_count; i++){
    long i_accession = accession_indices->a[i];
    fprintf(out_stream , "%s%s\n", accession_ids->a[i_accession], genotypes[i_accession]->a);
  }
  fclose(out_stream);

  double t3 = clock_time(the_clock);
  fprintf(stderr, "# Done. Times for input: %8.4lf, processing: %8.4lf, output: %8.4lf, total: %8.4lf\n",
	  t1-t0, t2-t1, t3-t2, t3-t0);

  // *****  cleanup  *****
  free_vchar(output_filename);
  free_vstr(marker_ids);
  free_vstr(accession_ids);
  free_vstr(marker_lines);
  for(long i=0; i<n_accessions; i++){
    free_vchar(genotypes[i]);
  }
  free(genotypes);
} // end of main


  //////////////////////////////////////////////////
  //         subroutine definitions               //  
  //////////////////////////////////////////////////

void* process_marker_range(void* x){

  TD* td = (TD*)x;
  Vstr* marker_lines = td->marker_lines;
  long first_marker = td->first_marker;
  long last_marker = td->last_marker;
  
  Vstr* marker_ids = td->marker_ids;
  Vchar** genotypes = td->genotypes;

  long marker_count = 0;
  char* line;
  for(long i_marker=first_marker; i_marker<=last_marker; i_marker++){
    line = marker_lines->a[i_marker];
      
    char* saveptr = line;
    Vchar* chromosome = construct_vchar_from_str(strtok_r(line, "\t \n\r", &saveptr)); // first token found is the chromosome number
    Vchar* position =  construct_vchar_from_str(strtok_r(NULL, "\t \n\r", &saveptr)); // next token is position within chromosome
    Vchar* marker_id;
    char* token =  strtok_r(NULL, "\t \n\r", &saveptr);    
    if(td->use_alt_marker_id){  
      marker_id = construct_vchar_from_str(chromosome->a);
      append_char_to_vchar(marker_id, '_');
      append_str_to_vchar(marker_id, position->a);
    }else{
      marker_id = construct_vchar_from_str(token); // strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    }
    marker_ids->a[i_marker] = marker_id->a; 

    free(marker_id); // free marker_id struct, but not its array of chars, which are now pointed to by marker_ids->a[i_marker]
    free_vchar(chromosome);
    free_vchar(position);
    char* ref_allele =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* alt_allele =  strtok_r(NULL, "\t \n\r", &saveptr);    

    for(long i = 1; i <= 3; i++){ // read the next 3 cols -
      token =  strtok_r(NULL, "\t \n\r", &saveptr);
    }

    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* format = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // format string, e.g. GT:DS:GP
    long GTidx, GQidx, GPidx;
    get_GT_GQ_GP_indices(format, &GTidx, &GQidx, &GPidx);
    free(format);
    long acc_index = 0;
    while(1){ // read genotypes from one line, i.e. one marker
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL)	break; // end of line has been reached.
      char genotype = token_to_genotype(token, GTidx, GPidx, td->minGP);
      genotypes[acc_index]->a[2*i_marker+1] = genotype;
      acc_index++;   
    } // done reading genotypes for all accessions of this marker
    assert(acc_index == td->n_accessions); // check that this line has number of accessions = number of accession ids.
    marker_count++;
  } // done reading all lines (markers)
  fprintf(stderr, "# A thread is done processing markers %ld through %ld; %ld markers x %ld accessions\n",
	  first_marker, last_marker, marker_count, td->n_accessions);
}

char token_to_genotype(char* token, long gtidx, long gpidx, double minGP){
  char result;
  char* saveptr;
  bool quality_ok = true; // (gpidx >= 0  &&  minGP > 0)? false : true;
  long idx = 0;
  char* tkn = strtok_r(token, ":", &saveptr);
  if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
    result = GTstr_to_dosage(tkn); // return result;
  }else if(idx == gpidx){
    quality_ok =
      GP_to_quality_ok(tkn, minGP);
  }
  idx++;
  
  while(1){ // get more subtokens
    tkn = strtok_r(NULL, ":", &saveptr);
    if(tkn == NULL)	break; // end of line has been reached.
    if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      result = GTstr_to_dosage(tkn);
    }else if(idx == gpidx){
      quality_ok =
	GP_to_quality_ok(tkn, minGP);
    }
    idx++;   
  }

  if(! quality_ok) result = 'X';
  return result;
}

bool GP_to_quality_ok(char* token, double minGP){
  if(minGP > 0.0){
    bool quality_ok = false;
    float p0, p1, p2;
    if(sscanf(token, "%f,%f,%f", &p0, &p1, &p2) == 3){
      quality_ok = (p0 >= minGP  ||  p1 >= minGP  || p2 >= minGP);
    }
    return quality_ok;
  }else{
    return true; 
  }
}

char GTstr_to_dosage(char* tkn){
  long d = 0;
  char a1 = tkn[0];
  char a2 = tkn[2];
  if(a1 == '1'){
    d++;
  }else if(a1 == '.'){
    return 'X';
  }
  if(a2 == '1'){
    d++;
  }else if(a2 == '.'){
    return 'X';
  }
  return (char)(d + 48);
}

void get_GT_GQ_GP_indices(char* format, long* GTp, long* GQp, long* GPp){
  *GQp = -1;
  *GPp = -1;
  long index = 0;
  char* saveptr = format;
  char* token = strtok_r(format, ":", &saveptr);
  if(strcmp(token, "GT") == 0) *GTp = index;
  if(strcmp(token, "GQ") == 0) *GQp = index;
  if(strcmp(token, "GP") == 0) *GPp = index;
  while(1){ // read in cols 1 through 8 "POS ID REF ..."
    index++;
    token = strtok_r(NULL, ":", &saveptr);
    if(token == NULL) break; 
    if(strcmp(token, "GT") == 0) *GTp = index;
    if(strcmp(token, "GQ") == 0) *GQp = index;
    if(strcmp(token, "GP") == 0) *GPp = index;   
  }
  return;
}


// unused

/* bool GP_to_quality_ok_xxx(char* token, double minGP){ */
/*    if(minGP > 0.0){ */
/*      char* saveptr; */
/*      char* tkn = strtok_r(token, ",", &saveptr); */
/*      if(atof(tkn) >= minGP) return true; */
/*      tkn = strtok_r(NULL, ",", &saveptr); */
/*      if(atof(tkn) >= minGP) return true; */
/*      tkn = strtok_r(NULL, ",", &saveptr); */
/*      if(atof(tkn) >= minGP) return true; */
/*      return false; */
/*    }else{ */
/*      return true; */
/*    } */
/*  } */

/* Vchar* token_to_plink_genotype(char* token, Vchar** alleles, long gtidx, long gpidx, double minGP){ */
/*   Vchar* plnkgt; */
/*   char* saveptr; */
/*   long idx = 0; */
/*   bool quality_ok = true; */
/*   char* tkn = strtok_r(token, ":", &saveptr); */
/*   if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1 */
/*     plnkgt = GT_to_plnkgt(tkn, alleles); */
/*   }else if(idx == gpidx){ */
/*     quality_ok = GP_to_quality_ok(tkn, minGP); */
/*   } */
/*   idx++;  */
/*   while(1){ // read genotypes from one line, i.e. one marker */
/*     tkn = strtok_r(NULL, ":", &saveptr); */
/*     if(tkn == NULL)	break; // end of line has been reached. */
/*     if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1 */
/* 	plnkgt = GT_to_plnkgt(tkn, alleles); */
/*     }else if(idx == gpidx){ */
/* 	quality_ok = GP_to_quality_ok(tkn, minGP); */
/*     } */
/*     idx++; */
/*   } */
/*   if(! quality_ok){ */
/*     free_vchar(plnkgt); */
/*     plnkgt = construct_vchar_from_str("\t0\t0"); */
/*   } */
/*   return plnkgt; */
/* } */

/* Vchar* GT_to_plnkgt(char* token, Vchar** alleles){ */
/*   // token would be for example "0/1", returns plnkgt, whose string part plnkgt->a would be e.g. "\tA\tACTA" */
/*   Vchar* plnkgt = construct_vchar(16); */
/*   append_char_to_vchar(plnkgt, '\t'); */
/*   if(token[0] == '0'){ */
/*     append_str_to_vchar(plnkgt, alleles[0]->a); */
/*   }else if(token[0] == '1'){ */
/*     append_str_to_vchar(plnkgt, alleles[1]->a); */
/*   }else{ */
/*     append_char_to_vchar(plnkgt, '0'); */
/*   } */
/*   append_char_to_vchar(plnkgt, '\t'); */
/*   if(token[2] == '0'){ */
/*     append_str_to_vchar(plnkgt, alleles[0]->a); */
/*   }else if(token[2] == '1'){ */
/*     append_str_to_vchar(plnkgt, alleles[1]->a); */
/*   }else{ */
/*     append_char_to_vchar(plnkgt, '0'); */
/*   } */
/*   return plnkgt; */
/* } */
