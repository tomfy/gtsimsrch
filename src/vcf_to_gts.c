#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
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

#define INIT_NMARKERS 6000

char token_to_genotype(char* token, long gtidx, long gpidx, double minGP);
Vchar* token_to_plink_genotype(char* token, Vchar** alleles, long gtidx, long gpidx, double minGP);
char GTstr_to_dosage(char* tkn);
void get_GT_GQ_GP_indices(char* format, long* GTidx, long* GQidx, long* GPidx);
Vchar* GT_to_plnkgt(char* token, Vchar** alleles);
bool GP_to_quality_ok(char* token, double minGP);

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

  bool plink = false;
  bool use_alt_marker_id = false;
    
  int c;
  while((c = getopt(argc, argv, "i:o:p:ka")) != -1){
    // i: input file name (required).
    // o: output file name. Default: "vcftogts"
    // p: minGP, there must be an gt with est. genotype prob. >= minGP, or it is considered missing data
    //   // x: max_marker_md_fraction. If a marker has > this proportion of missing data, omit from output. 
    // k: -k gives output suitable for plink input, default is output suitable for input to duplicatesearch
    // a: -a to use marker ids constructed from separate chromosome and position entries in vcf file.
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
	fprintf(stderr, "# conversion of argument %s to double (minGP) failed.\n", optarg);
	exit(EXIT_FAILURE);
      }else if(minGP > 1){
	fprintf(stderr, "# minGP was set to %8.4lf , must be <= 1\n", minGP);
	exit(EXIT_FAILURE);
      }
      fprintf(stderr, "# minGP set to: %8.5lf\n", minGP);
      break;
    /* case 'x': */
    /*    if(sscanf(optarg, "%lf ", &max_marker_md_fraction) != 1  ||  errno != 0){ */
    /* 	fprintf(stderr, "# conversion of argument %s to double (max_marker_md_fraction) failed.\n", optarg); */
    /* 	exit(EXIT_FAILURE); */
    /*   }else if(max_marker_md_fraction <= 0){ */
    /* 	fprintf(stderr, "# max_marker_md_fraction was set to %8.4lf , must be > 0\n", max_marker_md_fraction); */
    /* 	exit(EXIT_FAILURE); */
    /*   } */
    /*   fprintf(stderr, "# max_marker_md_fraction set to: %8.5lf\n", max_marker_md_fraction); */
    /*   break; */
    case 'k' :
      plink = true;
      break;
    case 'a' :
      use_alt_marker_id = true;
      break;
      
    }
  }
  if(plink){ // construct the two filenames with extensions .ped and .map
    append_str_to_vchar(marker_ids_filename, output_filename->a);
    append_str_to_vchar(output_filename, ".ped");    
    append_str_to_vchar(marker_ids_filename, ".map");    
  }else{
    // leave it alone
  }

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
      for(long ii = 1; ii <= 8; ii++){ // read in cols 1 through 8 "POS ID REF ..."
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

  Vchar* accession_genotypes[accid_count]; // We know how many accession ids (and therefore how many accessions) already.
  for(long i = 0; i < accid_count; i++){
    accession_genotypes[i] = construct_vchar(2*INIT_NMARKERS); 
  }  
  
  Vstr* marker_ids = construct_vstr(1000);
  Vstr* chromosome_ids = construct_vstr(1000);
  Vstr* positions = construct_vstr(1000);
  Vlong* marker_md_counts = construct_vlong(1000);
  long marker_md_count = 0;
  long marker_count = 0;
  Vchar* alleles[2]; // an array of 2 Vchar*s

  while((nread = getline(&line, &len, in_stream)) != -1){ // 
    saveptr = line;

    marker_md_count = 0;
    Vchar* chromosome = construct_vchar_from_str(strtok_r(line, "\t \n\r", &saveptr)); // first token found is the chromosome number
    push_to_vstr(chromosome_ids, chromosome->a);
    Vchar* position =  construct_vchar_from_str(strtok_r(NULL, "\t \n\r", &saveptr)); // next token is position within chromosome
    push_to_vstr(positions, position->a);
    Vchar* marker_id;
    char* token =  strtok_r(NULL, "\t \n\r", &saveptr);    
    if(use_alt_marker_id){
    
    marker_id = construct_vchar_from_str(chromosome->a);
    append_char_to_vchar(marker_id, '_');
    append_str_to_vchar(marker_id, position->a);
    }else{
      marker_id = construct_vchar_from_str(token); // strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    }
    push_to_vstr(marker_ids, marker_id->a);

    free(marker_id); // free marker_id struct, but not it's array of chars.
      free(chromosome);
      free(position);
      char* ref_allele =  strtok_r(NULL, "\t \n\r", &saveptr);
    alleles[0] = construct_vchar_from_str(ref_allele);
    char* alt_allele =  strtok_r(NULL, "\t \n\r", &saveptr);
    alleles[1] = construct_vchar_from_str(alt_allele);
    

    for(long i = 1; i <= 3; i++){ // read the next 3 cols -
      token =  strtok_r(NULL, "\t \n\r", &saveptr);
    }

    token =  strtok_r(NULL, "\t \n\r", &saveptr);
    char* format = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); // format string, e.g. GT:DS:GP
    long GTidx, GQidx, GPidx;
    get_GT_GQ_GP_indices(format, &GTidx, &GQidx, &GPidx);

    long acc_index = 0;
    while(1){ // read genotypes from one line, i.e. one marker
      
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      
      if(token == NULL)	break; // end of line has been reached.

      if(plink){ // gt will be something like "\tA\tC", but alleles can be multi-character
	Vchar* plink_gt = token_to_plink_genotype(token, /*format,*/ alleles, 0, GPidx, minGP);
	//	if(strcmp(plink_gt->a, "0") == 0) marker_md_count++;
	//	fprintf(stderr, "# # # : %s\n", plink_gt);
	append_str_to_vchar(accession_genotypes[acc_index], plink_gt->a);
	free_vchar(plink_gt);
      }else{ //
	//	fprintf(stderr, "#AAA: %ld %lf\n", GPidx, minGP);
	char genotype = token_to_genotype(token, /*format,*/ GTidx, GPidx, minGP); // i.e. GT:DS:GP
	//	if(genotype == 'X') marker_md_count++;
	   //	if(1){
	  char s[3] = "  "; // genotype is one char
	  s[1] = genotype;
	  append_str_to_vchar(accession_genotypes[acc_index], s);
	/* }else{ */
	/*   append_char_to_vchar(accession_genotypes[acc_index], ' '); // quicker to add spaces here, rather than when printing. */
	/*     append_char_to_vchar(accession_genotypes[acc_index], genotype); */
	
	/* } */
      }
     acc_index++;   
    } // done reading genotypes for all accessions of this marker
    if(acc_index != accession_ids->size) exit(EXIT_FAILURE);
    // push_to_vlong(marker_md_counts, marker_md_count);
    // fprintf(stderr, "# %ld  %ld \n", marker_md_counts->size, marker_md_count);

    free_vchar(alleles[0]);
    free_vchar(alleles[1]);
    free(format);
    
    marker_count++;
    if(marker_count % 500 == 0) fprintf(stderr, "markers read: %ld %10.4f\n", marker_count, hi_res_time() - t_start);
 
  } // done reading all lines (markers)
  free(line);
  fclose(in_stream);
  fprintf(stderr, "# Done reading all %ld markers x %ld accessions\n", marker_count, accid_count);


  //////////////////////////////////////////////////
  //                 output                       //
  //////////////////////////////////////////////////
  //  long max_marker_md_count = (long)(max_marker_md_fraction * accession_ids->size + 0.5);
  if(plink){
     for(long i=0; i<accid_count; i++){
      char* s = accession_genotypes[i]->a;
      fprintf(out_stream , "%s\t%s\t%s%s\n", accession_ids->a[i], accession_ids->a[i], "0\t0\t0\t0", s);
    }
     // and output marker ids in a separate file:
       FILE* out_stream_map = fopen(marker_ids_filename->a, "w");
       for(long i=0; i<marker_ids->size; i++){
	 fprintf(out_stream_map, "%s\t%s\t0\t%s\n", chromosome_ids->a[i], marker_ids->a[i], positions->a[i]);
       }
       fclose(out_stream_map);
     
  }else{
    fprintf(out_stream, "MARKER");
    for(long i=0; i<marker_ids->size; i++){
      // if(marker_md_counts[i] > max_marker_md_count)
	fprintf(out_stream, " %s", marker_ids->a[i]);
    }fprintf(out_stream, "\n");
    //  fprintf(stderr, "### %ld %ld  %s\n", accid_count, accession_ids->size, accession_ids->a[0]);
    for(long i=0; i<accid_count; i++){
      
      fprintf(out_stream , "%s%s\n", accession_ids->a[i], accession_genotypes[i]->a);
    }
  }
  fclose(out_stream);
  

  // clean up
  free_vstr(accession_ids);
  for(long i=0; i<accid_count; i++){
    free_vchar(accession_genotypes[i]);
  }
  free_vstr(marker_ids);
  free_vstr(chromosome_ids);
  free_vstr(positions);
  free_vchar(output_filename);
  free_vchar(marker_ids_filename);
  
  } // end of main


//////////////////////////////////////////////////
//         subroutine defintions                //  
//////////////////////////////////////////////////

char token_to_genotype(char* token, long gtidx, long gpidx, double minGP){
  char result;
  char* saveptr;
  bool quality_ok = true; // (gpidx >= 0  &&  minGP > 0)? false : true;
  long idx = 0;
  //  assert(format[0] == 'G' && format[1] == 'T');
  char* tkn = strtok_r(token, ":", &saveptr);
  if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      result = GTstr_to_dosage(tkn); // return result;
  }else if(idx == gpidx){
    quality_ok = GP_to_quality_ok(tkn, minGP);
  }
  idx++;
  
  while(1){ // get more subtokens
    tkn = strtok_r(NULL, ":", &saveptr);
    if(tkn == NULL)	break; // end of line has been reached.
    //  fprintf(stderr, "tkn: %s , idx: %ld %ld\n", tkn, idx, gtidx);
    if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      result = GTstr_to_dosage(tkn);
    }else if(idx == gpidx){
	quality_ok = GP_to_quality_ok(tkn, minGP);
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
      // fprintf(stderr, "%lf %lf %lf %lf  %ld\n", p0, p1, p2, minGP, (long)quality_ok);
    }
    return quality_ok;
  }else{
    return true; 
  }
}

Vchar* token_to_plink_genotype(char* token, Vchar** alleles, long gtidx, long gpidx, double minGP){
  Vchar* plnkgt;
  char* saveptr;
  long idx = 0;
  bool quality_ok = true; // (gpidx >= 0  &&  minGP > 0)? false : true;
  // fprintf(stderr, "%s  %ld %lf \n", token, gpidx, minGP);
  //  assert(format[0] == 'G' && format[1] == 'T');
  char* tkn = strtok_r(token, ":", &saveptr);
  //  fprintf(stderr, "A %s | %s | %s ;  %ld %ld ; %s %s\n", token, tkn, saveptr, idx, gtidx, alleles[0]->a, alleles[1]->a);
  if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      plnkgt = GT_to_plnkgt(tkn, alleles);
  }else if(idx == gpidx){
    quality_ok = GP_to_quality_ok(tkn, minGP);
  }
  idx++; 
  while(1){ // read genotypes from one line, i.e. one marker
    tkn = strtok_r(NULL, ":", &saveptr);
    if(tkn == NULL)	break; // end of line has been reached.
    //   fprintf(stderr, "tkn: %s , idx: %ld %ld\n", tkn, idx, gpidx);

    if(idx == gtidx){ // tkn should be e.g. 0|1 or 0/1 or 1/1
      plnkgt = GT_to_plnkgt(tkn, alleles);
    }else if(idx == gpidx){
      quality_ok = GP_to_quality_ok(tkn, minGP);
    }
    idx++;
  }
  if(! quality_ok){
    free_vchar(plnkgt);
    plnkgt = construct_vchar_from_str("\t0\t0");
  }
  //  fprintf(stderr, "#$@#: %s\n", plnkgt->a);
  return plnkgt;
}

Vchar* GT_to_plnkgt(char* token, Vchar** alleles){
  // token would be for example "0/1", returns plnkgt, whose string part plnkgt->a would be e.g. "\tA\tACTA"
  Vchar* plnkgt = construct_vchar(16);
  append_char_to_vchar(plnkgt, '\t');
    if(token[0] == '0'){
      append_str_to_vchar(plnkgt, alleles[0]->a);
    }else if(token[0] == '1'){
      append_str_to_vchar(plnkgt, alleles[1]->a);
    }else{
      append_char_to_vchar(plnkgt, '0');
    }
    append_char_to_vchar(plnkgt, '\t');
    if(token[2] == '0'){
      append_str_to_vchar(plnkgt, alleles[0]->a);
    }else if(token[2] == '1'){
      append_str_to_vchar(plnkgt, alleles[1]->a);
    }else{
      append_char_to_vchar(plnkgt, '0');
    }
    return plnkgt;
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
  //fprintf(stderr, "#A ## ## %s\n", token);
  if(strcmp(token, "GT") == 0) *GTp = index;
  if(strcmp(token, "GQ") == 0) *GQp = index;
  if(strcmp(token, "GP") == 0) *GPp = index;
  while(1){ // read in cols 1 through 8 "POS ID REF ..."
    index++;
    token = strtok_r(NULL, ":", &saveptr);
    if(token == NULL) break; 
    //fprintf(stderr, "#B ## ## %s\n", token);
    if(strcmp(token, "GT") == 0) *GTp = index;
    if(strcmp(token, "GQ") == 0) *GQp = index;
    if(strcmp(token, "GP") == 0) *GPp = index;
   
  }
  return;
}
