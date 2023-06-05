#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <sys/sysinfo.h>
#include <pthread.h>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>

#include "vect.h"
#include "various.h"

#define INIT_NMARKERS 6000

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
  Vchar** genotypes; // genotypes[i]->a[j] is genotypes of ith marker, jth accession
} TD; // thread data

void* process_marker_range(void* x);


char token_to_genotype(char* token, long gtidx, long gpidx, double minGP);
// Vchar* token_to_plink_genotype(char* token, Vchar** alleles, long gtidx, long gpidx, double minGP);
char GTstr_to_dosage(char* tkn);
void get_GT_GQ_GP_indices(char* format, long* GTidx, long* GQidx, long* GPidx);
// Vchar* GT_to_plnkgt(char* token, Vchar** alleles);
bool GP_to_quality_ok(char* token, double minGP);
bool GP_to_quality_ok_xxx(char* token, double minGP);

extern int errno; 

int main(int argc, char *argv[]){
  errno = 0;
  double t_start = hi_res_time();
  
  double minGQ = 0;
  double minGP = 0;
  // double max_marker_md_fraction = 0.25;

  char* input_filename = NULL;
  FILE *in_stream = NULL;
  Vchar* output_filename = construct_vchar_from_str("vcftogts");
  Vchar* marker_ids_filename = construct_vchar_from_str(""); 
  FILE* out_stream = NULL;

  //  bool plink = false;
  bool use_alt_marker_id = false;

  long nprocs = (long)get_nprocs();
  long Nthreads = (nprocs >= 3)? nprocs-2 : 1;

  bool oldway = false; // true;
    
  int c;
  while((c = getopt(argc, argv, "i:o:p:kaut:")) != -1){
    // i: input file name (required).
    // o: output file name. Default: "vcftogts"
    // p: minGP, there must be an gt with est. genotype prob. >= minGP, or it is considered missing data
    // a: -a to use marker ids constructed from separate chromosome and position entries in vcf file.
    // t: number of threads to use.
    switch(c){
    case 'i':
      input_filename = optarg;
      in_stream = fopen(input_filename, "r");
      if(in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
      //   fclose(in_stream);
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
      }else if(Nthreads < 1){
	fprintf(stderr, "# Setting Nthreads to min allowed value of 1.\n");
	Nthreads = 1;
      }
      break;
    case 'a' :
      use_alt_marker_id = true;
      break;
    }
  }
    /* if(plink){ // construct the two filenames with extensions .ped and .map */
    /*   append_str_to_vchar(marker_ids_filename, output_filename->a); */
    /*   append_str_to_vchar(output_filename, ".ped");     */
    /*   append_str_to_vchar(marker_ids_filename, ".map");     */
    /* }else{ // duplicatesearch */
    /*   // leave it alone */
    /* } */

    out_stream = fopen(output_filename->a, "w");
    if(out_stream == NULL){
      fprintf(stderr, "Failed to open %s for writing.\n", output_filename->a);
      exit(EXIT_FAILURE);
    }
  
    char* line = NULL;
    size_t len = 0;
    ssize_t nread;
  
    // ****************************************************
    // *****   Read first line; store accession ids.  *****
    // ****************************************************
  
    long accid_count = 0;
    Vstr* accession_ids = construct_vstr(1000);
    char* saveptr = NULL;
    while((nread = getline(&line, &len, in_stream)) != -1){
      saveptr = line;
      char* token = strtok_r(line, "\t \n\r", &saveptr);
    
      if((token == NULL) || (token[0] == '#' && token[1] == '#')) continue; // skip comments (starting with ##) and empty lines
      //  fprintf(stderr, "%s  ", token);
      if(token[0] == '#'){ // this the line with accession ids
	for(long ii = 1; ii <= 8; ii++){ // read in cols 1 through 8 "POS ID REF ..." but don't store them.
	  token = strtok_r(NULL, "\t \n\r", &saveptr);
	  //	fprintf(stderr, "%s ", token);
	}// fprintf(stderr, "\n");
	while(1){
	  token = strtok_r(NULL, "\t \n\r", &saveptr);
	  if(token == NULL) break;
	  // fprintf(stderr, "%s\n", token); // getchar();
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

    // *********************************************************
    // *****   Read the rest of the lines, one per marker  *****
    // *****   Each line has genotypes for all accessions  *****
    // *********************************************************
 
    Vstr* marker_lines = construct_vstr(1000);
    while((nread = getline(&line, &len, in_stream)) != -1){
      char* line_copy = (char*)malloc((nread+1)*sizeof(char));
      push_to_vstr(marker_lines, strcpy(line_copy, line));
    }

    
    long n_accessions = accid_count;
    long n_markers = marker_lines->size;

    clockid_t the_clock = CLOCK_MONOTONIC;
    struct timespec tspec;
    if(clock_getres(the_clock, &tspec) == 0){
      fprintf(stderr, "# clock resolution: %ld sec,  %ld nsec \n", (long)tspec.tv_sec, (long)tspec.tv_nsec);
    }else{
      exit(1);
    }

    TD* td;
    Vstr* marker_ids = construct_vstr(n_markers);
    marker_ids->size = n_markers;
    //fprintf(stderr, "#### %ld \n", n_markers); // getchar();

    char* str_of_spaces = (char*)malloc((2*n_markers + 1)*sizeof(char));
    for(long i=0; i<2*n_markers; i++){
      str_of_spaces[i] = ' ';
    }
    str_of_spaces[2*n_markers] = '\0'; // now a null terminated string, 2*n_markers spaces, followed by null.
    
    Vchar** genotypes = (Vchar**)malloc(n_accessions*sizeof(Vchar*));
    for(long i=0; i<n_accessions; i++){
      // genotypes[i] = construct_vchar(2*n_markers+1); // make cap 1 larger to accomodate terminating \0
      genotypes[i] = construct_vchar_from_str(str_of_spaces);
    }

    struct timespec tsp1;
    clock_gettime(the_clock, &tsp1);
    double starttime = tsp1.tv_sec + 1.0e-9*tsp1.tv_nsec;
  
    if(Nthreads == 0){    
      td = (TD*)malloc(sizeof(TD));
      td->n_accessions = n_accessions;
      td->marker_lines = marker_lines;
      td->first_marker = 0;
      td->last_marker = n_markers-1;
      td->use_alt_marker_id = use_alt_marker_id;
      td->minGP = minGP;
      td->marker_ids = marker_ids;
      td->genotypes = genotypes;

      process_marker_range((void*)td);
 
    }else{ // 1 or more pthreads   
      td = (TD*)malloc(Nthreads*sizeof(TD));
   
      td[0].n_accessions = n_accessions;
      td[0].marker_lines = marker_lines;
      td[0].first_marker = 0;
      td[0].last_marker = (long)((double)n_markers/Nthreads - 1);
      td[0].use_alt_marker_id = use_alt_marker_id;
      td[0].minGP = minGP;
      td[0].marker_ids = marker_ids;
      td[0].genotypes = genotypes;
    
      for(long i_thread = 1; i_thread<Nthreads; i_thread++){
	td[i_thread].n_accessions = n_accessions;
	td[i_thread].marker_lines = marker_lines;
	td[i_thread].first_marker = td[i_thread-1].last_marker + 1;
	td[i_thread].last_marker = (long)((double)(i_thread+1)*n_markers/Nthreads - 1);
	td[i_thread].use_alt_marker_id = use_alt_marker_id;
	td[i_thread].minGP = minGP;
	td[i_thread].marker_ids = marker_ids;
	td[i_thread].genotypes = genotypes;
      }
      td[Nthreads-1].last_marker = n_markers-1;
    
      int iret = 0;
      pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
      for(long i=0; i<Nthreads; i++){
	iret += pthread_create( thrids+i, NULL, process_marker_range, (void*) (td+i));
      }
      fprintf(stderr, "# iret: %d\n", iret);

      for(long i_thread=0; i_thread<Nthreads; i_thread++){
	pthread_join(thrids[i_thread], NULL);
      }
    }
    struct timespec tsp2;
    clock_gettime(the_clock, &tsp2);
    double endtime = tsp2.tv_sec + 1.0e-9*tsp2.tv_nsec;
    fprintf(stderr, "# time for threads to run %12.6lf\n", endtime-starttime); 
  
    //  *****  output  *****
    fprintf(out_stream, "MARKER");
    for(long i_marker=0; i_marker<td->marker_ids->size; i_marker++){
      fprintf(out_stream, " %s", td->marker_ids->a[i_marker]);
    }fprintf(out_stream, "\n");
    //  fprintf(stderr, "# # # n_accessions: %ld %ld \n", accid_count, n_accessions);
    for(long i_accession=0; i_accession<accid_count; i_accession++){      
      fprintf(out_stream , "%s%s\n", accession_ids->a[i_accession], genotypes[i_accession]->a);
    }
  
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

    //  fprintf(stderr, "# in pmr. A\n");
    // Vstr* chromosome_ids = construct_vstr(1000);
    // Vstr* positions = construct_vstr(1000);
    // Vlong* marker_md_counts = construct_vlong(1000);
    // long marker_md_count = 0;
    long marker_count = 0;
    //Vchar* alleles[2]; // an array of 2 Vchar*s

    //while((nread = getline(&line, &len, in_stream)) != -1){ //
    //fprintf(stderr, "# first, last markers: %ld %ld\n", first_marker, last_marker);
    char* line;
    for(long i_marker=first_marker; i_marker<=last_marker; i_marker++){
      line = marker_lines->a[i_marker];
      
      char* saveptr = line;
      // marker_md_count = 0;
      Vchar* chromosome = construct_vchar_from_str(strtok_r(line, "\t \n\r", &saveptr)); // first token found is the chromosome number
      // push_to_vstr(chromosome_ids, chromosome->a);
      Vchar* position =  construct_vchar_from_str(strtok_r(NULL, "\t \n\r", &saveptr)); // next token is position within chromosome
      // push_to_vstr(positions, position->a);
      Vchar* marker_id;
      char* token =  strtok_r(NULL, "\t \n\r", &saveptr);    
      if(td->use_alt_marker_id){  
	marker_id = construct_vchar_from_str(chromosome->a);
	append_char_to_vchar(marker_id, '_');
	append_str_to_vchar(marker_id, position->a);
      }else{
	marker_id = construct_vchar_from_str(token); // strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
      }

      // fprintf(stderr, "marker count, id: %ld %s \n", marker_count, marker_id->a);
      // push_to_vstr(marker_ids, marker_id->a);
      //  fprintf(stderr, "#marker id:  %s\n", marker_id->a); 
      marker_ids->a[i_marker] = marker_id->a; 
      //   fprintf(stderr, "# in pmr. B\n");

      free(marker_id); // free marker_id struct, but not its array of chars.
      free(chromosome);
      free(position);
      char* ref_allele =  strtok_r(NULL, "\t \n\r", &saveptr);
      //   alleles[0] = construct_vchar_from_str(ref_allele);
      char* alt_allele =  strtok_r(NULL, "\t \n\r", &saveptr);
      //   alleles[1] = construct_vchar_from_str(alt_allele);
    

      for(long i = 1; i <= 3; i++){ // read the next 3 cols -
	token =  strtok_r(NULL, "\t \n\r", &saveptr);
      }

      token =  strtok_r(NULL, "\t \n\r", &saveptr);
      char* format = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // format string, e.g. GT:DS:GP
      long GTidx, GQidx, GPidx;
      get_GT_GQ_GP_indices(format, &GTidx, &GQidx, &GPidx);

      long acc_index = 0;
      while(1){ // read genotypes from one line, i.e. one marker
	//  fprintf(stderr, "# in pmr. C %ld\n", acc_index);

	token = strtok_r(NULL, "\t \n\r", &saveptr);
	// token = strtok_r(NULL, "\t", &saveptr); 
	if(token == NULL)	break; // end of line has been reached.

	/* if(plink){ // gt will be something like "\tA\tC" or "\tAGC\tA"; alleles can be multi-character */
	/* 	Vchar* plink_gt = token_to_plink_genotype(token, /\*format,*\/ alleles, 0, GPidx, minGP); */
	/* 	append_str_to_vchar(accession_genotypes[acc_index], plink_gt->a); */
	/* 	free_vchar(plink_gt); */
	/* }else{ // */
	char genotype = token_to_genotype(token, /*format,*/ GTidx, GPidx, td->minGP); // i.e. GT:DS:GP
	//	char s[3] = "  "; // genotype is one char, 0, 1, 2, or X
	//	s[1] = genotype;
	
	// append_str_to_vchar(accession_genotypes[acc_index], s);
	//   genotypes[acc_index]->a[2*i_marker] = ' ';
	genotypes[acc_index]->a[2*i_marker+1] = genotype;
	//  }
	acc_index++;   
      } // done reading genotypes for all accessions of this marker
      assert(acc_index == td->n_accessions);
      // fprintf(stderr, "# n accessions: %ld %ld \n", acc_index, genotypes[i_marker]->capacity);
      /* if(acc_index == (genotypes[i_marker]->capacity-1)){ */
      /*         genotypes[i_marker]->length = genotypes[i_marker]->capacity; */
      /* }else{ */
      /*         exit(EXIT_FAILURE); */
      /* } */
      //  free_vchar(alleles[0]);
      //  free_vchar(alleles[1]);
      free(format);
      marker_count++;
      // if(marker_count % 500 == 0) fprintf(stderr, "# markers read so far: %ld %ld\n", marker_count, hi_res_time() - t_start);
    } // done reading all lines (markers)
    // fprintf(stderr, "### %ld %ld \n", marker_count, marker_ids->capacity);
    /* if(marker_count == marker_ids->capacity){ */
    /*   marker_ids->size = marker_count; */
    /* }else{ */
    /*   exit(EXIT_FAILURE); */
    /* } */
    // fprintf(stderr, "# in pmr. D\n"); //getchar();

    free(line);
    //fclose(in_stream);
    fprintf(stderr, "# A thread is done reading marker %ld through %ld; %ld markers x %ld accessions\n",
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
      //  GP_to_quality_ok_xxx(tkn, minGP);
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
	// GP_to_quality_ok_xxx(tkn, minGP);

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

  bool GP_to_quality_ok_xxx(char* token, double minGP){
    if(minGP > 0.0){
      char* saveptr;
      char* tkn = strtok_r(token, ",", &saveptr);
      if(atof(tkn) >= minGP) return true;
      tkn = strtok_r(NULL, ",", &saveptr);
      if(atof(tkn) >= minGP) return true;
      tkn = strtok_r(NULL, ",", &saveptr);
      if(atof(tkn) >= minGP) return true;
      return false;
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
