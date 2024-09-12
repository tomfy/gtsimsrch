// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <getopt.h>
#include <unistd.h>
#include <assert.h>
#include <stdbool.h>
#include "gtset.h"
#include "pedigree.h"

#define DO_ALL (0)

int do_checks = 0; // option -c sets this to 1 to do some checks.

//double hi_res_time(void);

void print_usage_info(FILE* stream);
Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents, double max_xhgmr, bool quick_xhgmr);
void sort_and_output_pedigrees(Vpedigree* the_pedigrees, long max_solns_out, FILE* o_stream, bool long_output_format);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();

  long alternative_pedigrees_level = 0; // find alt. pedigrees for  0: none, 1,2: accessions with bad pedigrees, 3,4: all accessions
  // 1,3: just accessions which are parents in pedigree file are considered as possible parents
  // 2,4: all accessions are considered are considered as possible parents.
  double max_marker_missing_data_fraction = 0.25; // default; control this with -x command line option.
  double max_accession_missing_data_fraction = 0.5;
  double min_minor_allele_frequency = 0.08; // 
  char* output_filename = "find_parents.out";
  double max_xhgmr = 0.18;
  double max_hgmr = 0.03;
  long max_candidate_parents = 80;
  bool quick_xhgmr = true;
  bool bitwise = false;
  long max_solns_out = 8;
  bool multiple_solns_on_one_line = true;
  // bool sort_by_z = false; // if false, sort by d
  bool use_xhgmr = true;
  // the following used in categorizing pedigrees from file as good or bad:
  double max_self_agmr = 0.04;
  double max_ok_hgmr = 0.04;
  double max_self_R = 0.04;
  double max_ok_d = 0.04;
  double max_ok_z = 0.04;
  bool long_output_format = false;
  
  double ploidy = 2;
  long Nthreads = 0;
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -in <dosages_file> [-out <output_filename> -xhgmr_max <max_xhgmr>] \n", argv[0]);
    print_usage_info(stderr);
    fprintf(stderr, "%d\n", (int)EXIT_FAILURE);
    exit(EXIT_FAILURE);
  }

  char* genotypes_filename = NULL;
  FILE *g_stream = NULL; // for reading in genotypes
  char* pedigrees_filename = NULL;
  FILE *p_stream = NULL; // for reading in pedigrees (if specified)

  int c;
  while(1){
    int option_index = 0;
    static struct option long_options[] = {
      {"input",   required_argument, 0,  'i'}, // filename of genotype data set
      {"pedigree", required_argument, 0, 'p'}, // filename of pedigree list
      {"output",  required_argument, 0,  'o'}, // output filename
      {"marker_max_missing_data", required_argument, 0, 'm'}, // markers with > this fraction missing data will not be used.
      {"maf_min", required_argument, 0, 'f'}, //
      {"accession_max_missing_data", required_argument, 0, 'a'},
      {"candidate_parents_max", required_argument, 0, 'c'},
      {"xhgmr_max", required_argument, 0, 'x'},
      {"help", no_argument, 0, 'h'},
      {"bitwise", no_argument, 0, 'b'},
      //  {"sort_by_z", no_argument, 0, 's'},
      {"check", no_argument, 0, 'k'}, // not implemented
      {"alternative_pedigrees_level", required_argument, 0, 'd'},
      {"hgmr_max", required_argument, 0, 'H'},
      {"agmr_self_max", required_argument, 0, 'S'},
      {"R_self_max", required_argument, 0, 'R'},
      {"d_max", required_argument, 0, 'D'},
      {"z_max", required_argument, 0, 'Z'},
      {0,         0,                 0,  0 }
    };
     
    c = getopt_long_only(argc, argv, "", long_options, &option_index);
    if(c == -1) break;
    switch(c){

  
      // while((c = getopt(argc, argv, "ci:x:o:m:h:r:D:p:")) != -1){
      // switch(c){

    case 'k':
      do_checks = 1;
      break;
    case 'i':
      genotypes_filename = optarg;
      g_stream = fopen(genotypes_filename, "r");
      if(g_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	exit(EXIT_FAILURE);
      }
      fclose(g_stream);
      break;
    case 'p':
      pedigrees_filename = optarg;
      p_stream = fopen(pedigrees_filename, "r");
      if(p_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", pedigrees_filename);
	exit(EXIT_FAILURE);
      }
      // fclose(p_stream);
      break;
    case 'o':
      output_filename = optarg;
      break;
    case 'c':
      max_candidate_parents = atoi(optarg);
      break;
    case 'm':
      if(optarg == 0){
	fprintf(stderr, "Option m requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_marker_missing_data_fraction = atof(optarg);
	if (max_marker_missing_data_fraction < 0){
	  fprintf(stderr, "# Max missing data fraction in markers will be set automatically.\n");
	}}
      break;
    case 'a':
      if(optarg == 0){
	fprintf(stderr, "Option a requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_accession_missing_data_fraction = atof(optarg);
	if (max_accession_missing_data_fraction < 0){
	  fprintf(stderr, "# max missing data fraction in accessions must be >= 0.\n");
	  exit(EXIT_FAILURE);
	}}
      break;
    case 'f': 
      min_minor_allele_frequency = (double)atof(optarg);
      if(min_minor_allele_frequency < 0){
	fprintf(stderr, "option f (min_minor_allele_frequency) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
     
    case 'x': 
      max_xhgmr = (double)atof(optarg);
      if(max_xhgmr < 0){
	fprintf(stderr, "option _xhgmr_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
 case 'H': 
      max_ok_hgmr = (double)atof(optarg);
      if(max_ok_hgmr < 0){
	fprintf(stderr, "option hgmr_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
 case 'S': 
      max_self_agmr = (double)atof(optarg);
      if(max_self_agmr < 0){
	fprintf(stderr, "option agmr_self_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
       case 'R': 
      max_self_R = (double)atof(optarg);
      if(max_self_R < 0){
	fprintf(stderr, "option R_self_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
       case 'D': 
      max_ok_d = (double)atof(optarg);
      if(max_ok_d < 0){
	fprintf(stderr, "option d_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
       case 'Z': 
      max_ok_z = (double)atof(optarg);
      if(max_ok_z < 0){
	fprintf(stderr, "option z_max requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;


      
    case 'd':
      alternative_pedigrees_level = (long)atoi(optarg);
      if(alternative_pedigrees_level < 0  ||  alternative_pedigrees_level > 4){
	fprintf(stderr, "option alternative_pedigrees requires one of 0, 1, 2, 3, or 4 as argument.\n");
	exit(EXIT_FAILURE);
      }
    case 'b':
      bitwise = true;
      break;
    /* case 's': */
    /*   sort_by_z = true; */
    /*   break; */
    case 'h':
      print_usage_info(stderr);
      exit(EXIT_FAILURE);
      break;
    case '?':
      fprintf(stderr, "? case in command line processing switch.\n");
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr, "default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
  if(optind < argc){
    fprintf(stderr, "Non-option arguments. Exiting.\n");
    exit(EXIT_FAILURE);
  }
  if(genotypes_filename == NULL){
    fprintf(stderr, "Must specify genotype filename: -i <filename>");

    fprintf(stderr, "%ld\n", (long)EXIT_FAILURE);
    exit(EXIT_FAILURE);
  }
  
  FILE *o_stream = NULL;
  o_stream = fopen(output_filename, "w");
  if(o_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
    exit(EXIT_FAILURE);
  }

  fprintf(o_stream, "# ");
  for(long i=0; i<argc; i++){
    fprintf(o_stream, "%s ", argv[i]);
  }fprintf(o_stream, "\n");
  // fprintf(stderr, "# genotypes file type: %d\n", genotype_file_type);
  //  char* geno_file_type = (genotype_file_type == DOSAGES)? "dosages" : "genotypes";
  fprintf(stdout, "# Genotypes filename: %s\n", genotypes_filename);
  // fprintf(o_stream, "# Genotypes filename: %s\n", genotypes_filename);
   
  // *****  done processing command line  *****

  // *************************************************************
  // ***  Read the genotypes file  *******************************
  // *************************************************************
  double t_a = hi_res_time();
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
  add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set, max_accession_missing_data_fraction, Nthreads); // load the new set of accessions
  print_vchar(stdout, the_genotypes_set->acc_filter_info);
  print_vchar(o_stream, the_genotypes_set->acc_filter_info);
  double t_b = hi_res_time();
  fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
  
  filter_genotypesset(the_genotypes_set, o_stream);
  print_vchar(stdout, the_genotypes_set->marker_filter_info);
  print_vchar(o_stream, the_genotypes_set->marker_filter_info);
  double t_c = hi_res_time();
  fprintf(stdout, "# Time to filter genotype data: %6.3f sec.\n", t_c - t_b);
  rectify_markers(the_genotypes_set); // needed for xhgmr to be fast.
  store_homozygs(the_genotypes_set);
  populate_marker_dosage_counts(the_genotypes_set);
  set_n_00_1_22_1s(the_genotypes_set);
  check_genotypesset(the_genotypes_set);
  set_Abits_Bbits(the_genotypes_set, Nthreads);
  Vidxid* the_gt_vidxid = construct_sorted_vidxid(the_genotypes_set->accessions); // ids and indexes of the_genotypes_set, sorted by id
  double t_d = hi_res_time();
  fprintf(stdout, "# Time for other initial processing of genotype data: %6.3f sec.\n", t_d - t_c);
  fflush(stdout);
  long n_gt_accessions = the_genotypes_set->accessions->size;
  

    // ****************************************************************
    // ***  If no pedigree file, or -alt 4 option  ********************
    // ***  Calculate pairwise quantities (agmr, xhgmr/hgmr, R?)  *****
    // ****************************************************************
  
   Viaxh** pairwise_info = (p_stream == NULL  ||  (alternative_pedigrees_level == 4))?
     calculate_pairwise_info(the_genotypes_set, max_candidate_parents, max_xhgmr, quick_xhgmr) : NULL;
  
    // ***************************************
    // ***  end of pairwise calculation  *****
    // ***************************************
 
   
   // ********************************************************* 
   // ***  Read the pedigrees file  ***************************
   // *********************************************************
   if(p_stream != NULL){ // have pedigree file
     fprintf(stderr, "# before read_and_store_pedigrees_3col\n");   
     const Vpedigree* pedigrees = read_and_store_pedigrees_3col(p_stream, the_gt_vidxid, the_genotypes_set); // acc id and then female, male parent ids in first 3 cols.
     Vpedigree* alt_pedigrees = construct_vpedigree(1000);
     fprintf(stderr, "# after read_and_store_pedigrees_3col\n");
     fclose(p_stream);
     fprintf(stdout, "# Done reading pedigree file. Time to read pedigree file: %6.3f\n", hi_res_time() - t_d);
     fprintf(stdout, "# Stored genotypes of %ld accessions, and  %ld pedigrees. \n", the_genotypes_set->accessions->size, pedigrees->size);

     /* */
     const Vlong* parent_idxs = accessions_with_offspring(pedigrees); // , the_genotypes_set->n_gt_accessions);
     fprintf(stdout, "# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);
     fprintf(stdout, "# Cumulative time so far: %6.3f sec.\n", hi_res_time() - t_begin_main);
    
     double t_start = hi_res_time();
     // char ploidy_char = (char)(the_genotypes_set->ploidy + 48);
     fprintf(stderr, "# Analyzing %ld pedigrees from table. \n", pedigrees->size);
     for(long i=0; i<pedigrees->size; i++){
       if(i % 100  == 0) fprintf(stdout, "# Done testing %ld pedigrees.\n", i);
       Pedigree* the_pedigree = pedigrees->a[i];
       Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(the_pedigree, the_genotypes_set); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
       the_pedigree->pedigree_stats = the_pedigree_stats;
       // assert(strcmp(the_pedigree->Accession->id, the_pedigree->A->id->a) == 0);
       Accession* A = the_pedigree->A;
       Accession* F = the_pedigree->F;
       Accession* M = the_pedigree->M;
       long A_gtset_idx = index_of_id_in_vidxid(the_gt_vidxid, A->id->a);
       if(0 || (F != NULL  ||  M != NULL)){ // at least one parent specified in pedigree
    
	 /* fprintf(o_stream, "%20s %5ld %20s %20s %ld  ", */
	 /* 	 A->id->a, A->missing_data_count, */
	 /* 	 (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA", the_pedigree_stats->all_good_count); */
	 fprintf(o_stream, "%s  P  ", A->id->a); // progeny accession and 'P' to indicate these are the parents from the pedigree file.
	 print_pedigree(o_stream, the_pedigree, long_output_format);
	 //print_pedigree_stats(o_stream, the_pedigree_stats, long_output_format); // print the stats for this pedigree

	 if(alternative_pedigrees_level == 1){
	   fprintf(stderr, "doing alt pedigrees for %s\n", A->id->a);
	   if(pedigree_ok(the_pedigree_stats, max_self_agmr, max_ok_hgmr, max_self_R, max_ok_d) == 0){	    
	     alt_pedigrees = pedigree_alternatives(the_pedigree, the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	     print_pedigree_alternatives(o_stream, alt_pedigrees, 4, long_output_format);
	     fprintf(o_stream, "\n");
	   }
	 }else if(alternative_pedigrees_level == 2){
	   fprintf(stderr, "alternative_pedigrees level 2 not implemented. Bye.\n"); exit(0);
	 }else if((alternative_pedigrees_level == 3)){
	   alt_pedigrees = pedigree_alternatives(the_pedigree, the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	   print_pedigree_alternatives(o_stream, alt_pedigrees, 4, long_output_format);
	   fprintf(o_stream, "\n");
	 }else if(alternative_pedigrees_level == 4){
	
	   // calculate d, z, etc. for triples involving the candidate parents in pairwise_info[i]	
	   Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pairwise_info[A_gtset_idx], max_candidate_parents);
	   //  fprintf(o_stream, "%s  %ld   ", prog->id->a, pairwise_info[i]->size); // output progeny id
	   sort_and_output_pedigrees(alt_pedigrees, max_solns_out, o_stream, long_output_format);
	   if(alt_pedigrees->size == 0) fprintf(o_stream, " this accession has no candidate parents.");
	   fprintf(o_stream, "\n");
	   free_vpedigree(alt_pedigrees);
	 }
	 alt_pedigrees->size = 0;
	 free(the_pedigree_stats);
       } //    
     } // end loop over pedigrees
     free_vpedigree(alt_pedigrees);
     /**/
    
  }else{ // no pedigree file, consider all accessions as possible parents
 
    // **************************************************
    // ***  evaluate parent1-parent2-progeny triples  ***
    // **************************************************

    double t_x = hi_res_time();
    long count_accs_w_no_cand_parents = 0;
    long count_accs_w_too_many_cand_parents = 0;
    
    for(long i=0; i<n_gt_accessions; i++){ 
      Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents.
	
      // calculate d, z, etc. for triples involving the candidate parents in pairwise_info[i]	
      Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(prog, the_genotypes_set, pairwise_info[i], max_candidate_parents);
      fprintf(o_stream, "%s  %ld   ", prog->id->a, pairwise_info[i]->size); // output progeny id
      sort_and_output_pedigrees(alt_pedigrees, max_solns_out, o_stream, long_output_format);
      if(alt_pedigrees->size == 0) fprintf(o_stream, " this accession has no candidate parents.");
      fprintf(o_stream, "\n");
      free_vpedigree(alt_pedigrees);
    } // end loop over offspring accessions

    fprintf(stderr, "# time for triple calculation: %8.3f\n", hi_res_time() - t_x);
    fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr);
    fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
    fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	    max_candidate_parents, count_accs_w_too_many_cand_parents);
  }
  
  // ********************  cleanup  **************************
  fclose(o_stream);
  free_genotypesset(the_genotypes_set);
  //  fprintf(stderr, "after free_genotypesset\n");
  free_vidxid(the_gt_vidxid);
  //free_vpedigree(alt_pedigrees);
  //  fprintf(stderr, "after free_vpedigree\n");
  //   free_vlong(parent_idxs);
  // fprintf(stderr, "# Done with cleanup\n");
  // getchar();
  fprintf(stdout, "# Total time: %10.4lf sec.\n", hi_res_time() - t_begin_main);
}
// **********************************************************
// ********************  end of main  ***********************
// **********************************************************




// **********************************************************
// ******************  functions  ***************************
// **********************************************************

/* double hi_res_time(void){ */
/*   return (double)clock()/(double)CLOCKS_PER_SEC; */
/* } */

void print_usage_info(FILE* stream){
  // c: do checks,
  // i: dosages input filename 
  // x: max fraction of missing data for markers,
  // o: output filename.
  // m: min minor allele frequency
  // h: max xhgmr for candidate parents
  // p: max_candidate_parents
  fprintf(stream, " -input <filename>               input filename (dosages matrix)\n");
  fprintf(stream, " -output <filename>              output filename (default: find_parents.out) \n");
  fprintf(stream, " -marker_max_missing_data <f>    don't use markers with missing data fraction > f\n");
  fprintf(stream, " -accession_max_missing_data <f> accessions with missing data fraction >f are ignored > f\n");
  fprintf(stream, " -maf_min <f>                    don't use markers with minor allele frequency < f\n");
  fprintf(stream, " -xhgmr_max <f>                  candidate parents must have xhgmr <= f\n");
  fprintf(stream, " -candidate_parents_max <n>      sort candidate parents by xhgmr and use only best n\n");
  fprintf(stream, " -help                           print this message and exit.\n");
}

Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents, double max_xhgmr, bool quick_xhgmr){
  long n_gt_accessions = the_genotypes_set->accessions->size;
  Viaxh** pairwise_info = (Viaxh**)malloc(n_gt_accessions*sizeof(Viaxh*)); // array of Viaxh*. pairwise_info[i]->a[j].xhgmr 
  for(long ii=0; ii<n_gt_accessions; ii++){
    pairwise_info[ii] = construct_viaxh(2*max_candidate_parents); // vector to hold candidate parents of accession with index ii
  }
    
  // if(use_xhgmr){ // Calculate xhgmr
  calculate_xhgmrs(the_genotypes_set, pairwise_info, quick_xhgmr, max_xhgmr);
  
  /* } else{  // hgmr - bitwise */
  /*   calculate_hgmrs(the_genotypes_set, pairwise_info, max_hgmr); */
  /*   fprintf(stdout, "# time for hgmrs: %10.3f \n", hi_res_time() - t0); */
  /* } */
  return pairwise_info;
}

void sort_and_output_pedigrees(Vpedigree* the_pedigrees, long max_solns_out, FILE* o_stream, bool long_output_format){
  // *******************************************************************************
  // ***  sort alt_pedigrees_array, output best parent pairs for this accession  ***
  // *******************************************************************************
  if(the_pedigrees->size > 0){
    if(the_pedigrees->size > 1)
      // sort_vpedigree_by_maxh1h2z(the_pedigrees);
      sort_vpedigree_by_maxdz(the_pedigrees);
      long n_out = (max_solns_out < the_pedigrees->size)? max_solns_out : the_pedigrees->size;
      for(long iii=0; iii < n_out; iii++){
	/* Accession* par1 = the_pedigrees->a[iii]->F; */
	/* Accession* par2 = the_pedigrees->a[iii]->M; */
	/* Pedigree_stats* the_ps = the_pedigrees->a[iii]->pedigree_stats; */
	/* // if((!multiple_solns_on_one_line) || (iii == 0)) fprintf(o_stream, "%s  %ld   ", prog->id->a, n_ok_parents); // output progeny id */
	/* fprintf(o_stream, "%s %s %ld  ", par1->id->a, par2->id->a, the_ps->all_good_count); // output candidate pair of parents */
	/* print_pedigree_stats(o_stream, the_ps, long_output_format); */
	fprintf(o_stream, "A  "); // to indicate alternative parents (not those from pedigree file)
	print_pedigree(o_stream, the_pedigrees->a[iii], long_output_format);
	// if((!multiple_solns_on_one_line) || (iii == n_out-1)) fprintf(o_stream, "\n");
      } // loop over solutions for one progeny accession
      //fprintf(o_stream, "\n");
  }
} // end of sort_and_output_pedigrees
