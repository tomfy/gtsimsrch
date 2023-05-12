#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>

#include "vect.h"
#include "various.h"

char token_to_genotype(char* token, char* format, long gtidx, long gpidx, double minGP);
char GTstr_to_dosage(char* tkn);
void get_GQ_GP_indices(char* format, long* GQidx, long* GPidx);


int main(int argc, char *argv[]){
  double t_start = hi_res_time();
  
  double minGQ = 0;
  double minGP = 0;

  char* input_filename = NULL;
  FILE *in_stream = NULL;
  char* output_filename = "vcftogts.out";
  FILE* out_stream = NULL;
    
  int c;
  while((c = getopt(argc, argv, "i:o:")) != -1){
    // i: input file name (required).
    // o: output file name. Default: "duplicatesearch.out"
  
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
      output_filename = optarg;
      break;
    }
  }

  out_stream = fopen(output_filename, "w");
  if(out_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
    exit(EXIT_FAILURE);
  }
  
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  // *****   Read first line; store marker ids.  *****
  long accid_count = 0;
  Vstr* accession_ids = construct_vstr(1000);
  char* saveptr = NULL;
  while((nread = getline(&line, &len, in_stream)) != -1){
    saveptr = line;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    
    if((token == NULL) || (token[0] == '#' && token[1] == '#')) continue; // skip comments (starting with ##) and empty lines
    //  fprintf(stderr, "%s  ", token);
    if(token[0] == '#'){ // this the line with accession ids
      for(long ii = 1; ii <= 8; ii++){ // read in cols 1 through 8 "POS ID REF ..."
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	//	fprintf(stderr, "%s ", token);
      }fprintf(stderr, "\n");
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

  Vchar** accession_genotypes = (Vchar**)malloc(accid_count * sizeof(Vchar*));
  // initialize the Vstr's:
  for(long i = 0; i < accid_count; i++){
    accession_genotypes[i] = construct_vchar(12000);
    // fprintf(stderr, "%ld  %ld %ld \n", i, accession_genotypes[i]->capacity, accession_genotypes[i]->length);
  }
  
  /* fprintf(stderr, "MARKER"); */
  /* for(long i=0; i<accession_ids->size; i++){ */
  /*   fprintf(stderr, " %s", accession_ids->a[i]); */
  /* } */
  /* fprintf(stderr, "\n"); */
  
  Vstr* marker_ids = construct_vstr(1000);
  // Read in the rest of the lines (markers)
  long marker_count = 0;
  while((nread = getline(&line, &len, in_stream)) != -1){ // 
    saveptr = line;
    
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    char* chromosome = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* position = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // might want to construct a marker id from chromosome and position    
    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    
    char* marker_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    push_to_vstr(marker_ids, marker_id);
       
    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* ref_allele = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // needed for plink format
    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* alt_allele = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // needed for plink format

    for(long i = 1; i <= 3; i++){ // read the next 3 cols -
      token =  strtok_r(NULL, "\t \n\r", &saveptr);
    }

    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* format = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // format string, e.g. GT:DS:GP
    long GQidx, GPidx;
    get_GQ_GP_indices(format, &GQidx, &GPidx);

    //  fprintf(stderr, "%ld %s %s %s %s %s %s \n", marker_count, chromosome, position, marker_id, ref_allele, alt_allele, format);
    long accession_count = 0;
    //  long marker_missing_data_count = 0;
  

    long acc_index = 0;
    while(1){ // read genotypes from one line, i.e. one marker
      
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      
      if(token == NULL)	break; // end of line has been reached.

      char genotype = token_to_genotype(token, format, 0, GPidx, minGP); // i.e. GT:DS:GP
      if(1){
	char s[3] = "  "; // genotype is one char (will need to do something more for plink)
	s[0] = genotype;
	append_str_to_vchar(accession_genotypes[acc_index], s);
      }else{
	append_char_to_vchar(accession_genotypes[acc_index], genotype);
	append_char_to_vchar(accession_genotypes[acc_index], ' ');
      }
      //     fprintf(stderr, "%ld  [%s]\n", accession_genotypes[acc_index]->length, accession_genotypes[acc_index]->a);
      accession_count++;
      //  fprintf(stderr, "mrkr, accession counts: %ld %ld\n", marker_count, accession_count);
      acc_index++;
    } // done reading genotypes for all accessions of this marker
   
    marker_count++;
    if(marker_count % 200 == 0) fprintf(stderr, "markers read: %ld %10.4f\n", marker_count, hi_res_time() - t_start);
  } // done reading all lines (markers)
  fprintf(stderr, "# Done reading all %ld markers x %ld accessions\n", marker_count, accid_count);
  fprintf(out_stream, "MARKER");
  for(long i=0; i<marker_ids->size; i++){
    fprintf(out_stream, " %s", marker_ids->a[i]);
  }fprintf(out_stream, "\n");
  //  fprintf(stderr, "### %ld %ld  %s\n", accid_count, accession_ids->size, accession_ids->a[0]);
  for(long i=0; i<accid_count; i++){
    // fprintf(out_stream, "# %ld \n", i);
    char* s = accession_genotypes[i]->a;
    //  fprintf(stderr, "%ld  %ld  \n", i, accession_genotypes[i]->length);
    // s[50] = '\0';
    fprintf(out_stream , "%s %s\n", accession_ids->a[i], s);
  }
  fclose(in_stream);
  fclose(out_stream);
} // end of main


// subroutine defintions
char token_to_genotype(char* token, char* format, long gtidx, long gpidx, double minGP){
  char result;
  char* saveptr;
  bool quality_ok = (gpidx >= 0  &&  minGP > 0)? false : true;
  long idx = 0;
  //  fprintf(stderr, "%s   ", token);
  assert(format[0] == 'G' && format[1] == 'T');
  char* tkn = strtok_r(token, ":", &saveptr);
  // fprintf(stderr, "A %s | %s | %s\n", token, tkn, saveptr);
  if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
    result = GTstr_to_dosage(tkn);
  }else if(idx == gpidx){
    if(minGP > 0){ // filter on genotype probability
      double p0, p1, p2;
      if(sscanf(tkn, "%f,%f,%f", &p0, &p1, &p2) == 3){
	quality_ok = (p0 >= minGP  ||  p1 >= minGP  || p2 >= minGP);
      }
    }
  }
  idx++; 
  while(1){ // read genotypes from one line, i.e. one marker
    tkn = strtok_r(NULL, "\t \n\r", &saveptr);
    if(tkn == NULL)	break; // end of line has been reached.
    if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      result = GTstr_to_dosage(tkn);
    }else if(idx == gpidx){
      if(minGP > 0){ // filter on genotype probability
	double p0, p1, p2;
	if(sscanf(tkn, "%f,%f,%f", &p0, &p1, &p2) == 3){
	  quality_ok = (p0 >= minGP  ||  p1 >= minGP  || p2 >= minGP);
	}
      }
    }
    idx++; 
    if(! quality_ok) result = '0';
  }    
  return result;
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

void get_GQ_GP_indices(char* format, long* GQp, long* GPp){
  *GQp = -1;
  *GPp = -1;
  //fprintf(stderr, "# # # %ld %ld \n", *GQp, *GPp);
  long index = 0;
  char* saveptr = format;
  char* token = strtok_r(format, ":", &saveptr);
  //fprintf(stderr, "#A ## ## %s\n", token);
  if(strcmp(token, "GQ") == 0) *GQp = index;
  if(strcmp(token, "GP") == 0) *GPp = index;
  while(1){ // read in cols 1 through 8 "POS ID REF ..."
    index++;
    token = strtok_r(NULL, ":", &saveptr);
    if(token == NULL) break; 
    //fprintf(stderr, "#B ## ## %s\n", token);
    if(strcmp(token, "GQ") == 0) *GQp = index;
    if(strcmp(token, "GP") == 0) *GPp = index;
    //fprintf(stderr, "ZZZZZZZZZZZZ: %ld %ld \n", *GQp, *GPp);

  }
  //fprintf(stderr, "ZZZ: %ld %ld \n", *GQp, *GPp);
  return;
}
