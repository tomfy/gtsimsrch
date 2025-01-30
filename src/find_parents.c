// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/sysinfo.h> // needed for get_nprocs
#include <assert.h>
#include <stdbool.h>
#include "gtset.h"
#include "pedigree.h"

#define DO_ALL (0)

#define DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION  0.2
#define DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION  0.5
#define DEFAULT_MIN_MAF  0.1
#define DEFAULT_MAX_XHGMR 0.18
#define DEFAULT_MAX_CANDIDATE_PARENTS 80

int do_checks = 0; // option -c sets this to 1 to do some checks.

//double hi_res_time(void);

void print_usage_info(FILE* stream);
Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents, double max_xhgmr, bool quick_xhgmr);
void sort_and_output_pedigrees(GenotypesSet* the_gtsset, Vpedigree* the_pedigrees, long max_solns_out, FILE* o_stream); //, bool long_output_format, double d_scale_factor);
void set_scaled_d_in_one_pedigree(Pedigree_stats* the_ps, double d_scale_factor);
void set_scaled_d_in_vpedigree(Vpedigree* the_pedigrees, double d_scale_factor);
void print_Xover2_rates(FILE* fh, Xcounts_2 X2);
void print_X3_info(FILE* o_stream, Xcounts_3 X3);
void print_Xover_rates(FILE* fh, Xcounts_3 X3);
// void d_z_of_random_set(GenotypesSet* gtset, long sample_size, double* mean_d, double* mean_z);
void random_set_means(GenotypesSet* gtset, long sample_size);
// double* mean_hgmr, double* mean_R, double* mean_d, double* mean_z);
// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();

  long alternative_pedigrees_level = 0; // find alt. pedigrees for  0: none, 1: accessions with bad pedigrees, 2: all accessions
  // 1,2: all accessions are considered are considered as possible parents.
  // 3,4: same as 1,2, respectively, but only accessions which are parents in pedigree file are considered as possible parents (not implemented)
  double max_marker_missing_data_fraction = DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION; // default; control this with -x command line option.
  double max_accession_missing_data_fraction = DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION;
  double min_minor_allele_frequency = DEFAULT_MIN_MAF; // 
  char* output_filename = "find_parents.out";
  double max_xhgmr = DEFAULT_MAX_XHGMR;
  // double max_hgmr = 0.03;
  long max_candidate_parents = 80;
  bool quick_xhgmr = true;
  //bool bitwise = false;
  long max_solns_out = 4;
  bool multiple_solns_on_one_line = true;
  // bool sort_by_z = false; // if false, sort by d
  //bool use_xhgmr = true;
  
  // the following used in categorizing pedigrees from file as good or bad:
  double max_self_agmr = 0.03;
  double max_ok_hgmr = 0.16;
  double max_self_R = 0.1;
  double max_ok_d = 0.15;
  double max_ok_z = 0.15;

  bool long_output_format = false; // true; true -> output denominators as well as 
  double d_scale_factor = 1.17; // sort on max(d_scale_factor*d, z)
  
  double ploidy = 2;
  long nprocs = (long)get_nprocs(); // returns 2*number of cores if hyperthreading.
  long Nthreads = -1; // (nprocs > 2)? nprocs/2 : 1; // default number of threads

  unsigned rand_seed = time(0); // (unsigned)tspec.tv_nsec;
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
      //  {"bitwise", no_argument, 0, 'b'},
      //  {"sort_by_z", no_argument, 0, 's'},
      {"check", no_argument, 0, 'k'}, // not implemented
      {"alternative_pedigrees_level", required_argument, 0, 'd'},
      {"seed", required_argument, 0, 'r'},
      {"hgmr_max", required_argument, 0, 'H'},
      {"agmr_self_max", required_argument, 0, 'S'},
      {"R_self_max", required_argument, 0, 'R'},
      {"d_max", required_argument, 0, 'D'},
      {"z_max", required_argument, 0, 'Z'},
      {"scale_factor", required_argument, 0, 'F'}, 
      //   {"threads", required_argument, 0, 't'},
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
    /* case 't': */
    /*   Nthreads = atoi(optarg); */
    /*   break; */
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
     case 'r': 
      rand_seed  = (unsigned)atoi(optarg);
      if(max_xhgmr < 0){
	fprintf(stderr, "option _xhgmr_max requires an real argument >= 0\n");
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
    case 'F': 
      d_scale_factor = (double)atof(optarg);
      if( d_scale_factor <= 0){
	fprintf(stderr, "option d_scale_factor requires an real argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;

      
    case 'd':
      alternative_pedigrees_level = (long)atoi(optarg);
      if(alternative_pedigrees_level < 0  ||  alternative_pedigrees_level > 4){
	fprintf(stderr, "option alternative_pedigrees requires one of 0, 1, 2, 3, or 4 as argument.\n");
	exit(EXIT_FAILURE);
      }
      break;
    /* case 'b': */
    /*   bitwise = true; */
    /*   break; */
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

  srand(rand_seed);
  // *****  done processing command line  *****

  // *************************************************************
  // ***  Read the genotypes file  *******************************
  // *************************************************************
  fprintf(stdout, "# Using %ld threads.\n", (long)Nthreads);
  fprintf(o_stream, "# Rng seed: %ld\n", (long)rand_seed);
  double t_a = hi_res_time();
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
  add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set, max_accession_missing_data_fraction, Nthreads); // load the new set of accessions
  print_vchar(stdout, the_genotypes_set->acc_filter_info);
  print_vchar(o_stream, the_genotypes_set->acc_filter_info);
  double t_b = hi_res_time();
  fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
  
  filter_genotypesset(the_genotypes_set, o_stream);
  fprintf(stderr, "after filter_genotypesset\n");
  print_vchar(stdout, the_genotypes_set->marker_filter_info);
  print_vchar(o_stream, the_genotypes_set->marker_filter_info);

  double t_c = hi_res_time();
  fprintf(stdout, "# Time to filter genotype data: %6.3f sec.\n", t_c - t_b);
  // rectify_markers(the_genotypes_set); // needed for xhgmr to be fast.
  store_homozygs(the_genotypes_set); // needed?
  fprintf(stderr, "before set_chromosome_start_indices.\n");
  if(the_genotypes_set->phased) set_chromosome_start_indices(the_genotypes_set);
  /*
  // populate_marker_dosage_counts(the_genotypes_set); // needed?
  for(long i=0; i<the_genotypes_set->accessions->size; i++){
    set_n2exp0s(the_genotypes_set, i);
  }
  // set_n_00_1_22_1s(the_genotypes_set);
  /* */
  check_genotypesset(the_genotypes_set);
  set_Abits_Bbits(the_genotypes_set, Nthreads);
  Vidxid* the_gt_vidxid = construct_sorted_vidxid(the_genotypes_set->accessions); // ids and indexes of the_genotypes_set, sorted by id
  double t_d = hi_res_time();
  fprintf(stdout, "# Time for other initial processing of genotype data: %6.3f sec.\n", t_d - t_c);

  // double mean_hgmr, mean_R, mean_D, mean_Z;
  long sample_size = 1000;
  random_set_means(the_genotypes_set, sample_size); // get and store mean hgmr, R, etc.
  the_genotypes_set->d_scale_factor = d_scale_factor;
  fprintf(stderr, "mean hgmr, R, d, z: %8.5f  %8.5f  %8.5f  %8.5f\n",
	  the_genotypes_set->mean_hgmr, the_genotypes_set->mean_R, the_genotypes_set->mean_d, the_genotypes_set->mean_z); 
  fflush(stdout);
  long n_gt_accessions = the_genotypes_set->accessions->size;
  // print_genotypesset(stderr, the_genotypes_set);

  // ****************************************************************
  // ***  If no pedigree file, or -alt 3 option  ********************
  // ***  Calculate pairwise quantities (agmr, xhgmr/hgmr, R?)  *****
  // ****************************************************************
  
  Viaxh** pairwise_info = (p_stream == NULL  ||  (alternative_pedigrees_level == 3))?
    calculate_pairwise_info(the_genotypes_set, max_candidate_parents, max_xhgmr, quick_xhgmr) : NULL;

  fprintf(stdout, "# time for calculate_pairwise_info (xhgmrs): %.5f\n", hi_res_time() - t_d);
  // ***************************************
  // ***  end of pairwise calculation  *****
  // ***************************************

   double initialization_time; 
 
  if(p_stream != NULL){ // have pedigree file
    // ********************************************************* 
    // ***  Read the pedigrees file  ***************************
    // *********************************************************
    fprintf(stderr, "About to read_and_store_pedigrees...\n");
    const Vpedigree* pedigrees = read_and_store_pedigrees_3col(p_stream, the_gt_vidxid, the_genotypes_set); // acc id and then female, male parent ids in first 3 cols.
    fclose(p_stream); 
    fprintf(stdout, "# Done reading pedigree file. Time to read pedigree file: %6.3f\n", hi_res_time() - t_d);
    fprintf(stdout, "# Stored genotypes of %ld accessions, and  %ld pedigrees. \n", the_genotypes_set->accessions->size, pedigrees->size);
     
    const Vlong* parent_idxs = accessions_with_offspring(pedigrees, the_genotypes_set->accessions->size); // , the_genotypes_set->n_gt_accessions);
    fprintf(stdout, "# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);
  
    initialization_time = hi_res_time(); 
    fprintf(stdout, "# Cumulative time so far: %6.3f sec.\n", initialization_time - t_begin_main);
     fprintf(stderr, "# Analyzing %ld pedigrees from table. \n", pedigrees->size);
    
     for(long i=0; i<pedigrees->size; i++){
       if(i % 100  == 0) fprintf(stdout, "# Done testing %ld pedigrees.\n", i);
       Pedigree* the_pedigree = pedigrees->a[i];
       Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(the_pedigree, the_genotypes_set); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
       //set_scaled_d_in_one_pedigree(the_pedigree_stats, d_scale_factor);
       the_pedigree->pedigree_stats = the_pedigree_stats;
      
       Accession* A = the_pedigree->A;
       Accession* F = the_pedigree->F;
       Accession* M = the_pedigree->M;
       
       long A_gtset_idx = index_of_id_in_vidxid(the_gt_vidxid, A->id->a);
       if(0 || (F != NULL  ||  M != NULL)){ // at least one parent specified in pedigree
    
	 fprintf(o_stream, "%s  x  P  ", A->id->a); // progeny accession and 'P' to indicate these are the parents from the pedigree file.
	 print_pedigree_normalized(o_stream, the_pedigree); //, the_genotypes_set);
	 two_doubles hratios = heterozyg_ratios(A, F);
	 
	 if(the_genotypes_set->phased){
	   Xcounts_3 X3 = count_crossovers(the_genotypes_set, F, M, A);
	   print_Xover_rates(o_stream, X3);
	   fprintf(o_stream, "  %7.5f  %7.5f ", hratios.x1, hratios.x2);
	 } // end of if phased branch

	 if(alternative_pedigrees_level == 1){ // iff pedigrees bad, do search for parents, considering only accessions in parent_idxs as possible parents
	   long pedigree_ok_return_value = pedigree_ok(the_pedigree_stats, max_self_agmr, max_self_R, max_ok_d, max_ok_z);
	   if(pedigree_ok_return_value <= 0){ // if pedigree is 'bad' look for alternatives 
	     Vpedigree* alt_pedigrees = pedigree_alternatives(the_pedigree, the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	     print_pedigree_alternatives(o_stream, alt_pedigrees, max_solns_out, long_output_format);
	     free_vpedigree(alt_pedigrees);
	   }
	 }else if(alternative_pedigrees_level == 2){ // iff pedigrees bad, do search for parents. Consider all accessions as possible parents
	   long pedigree_ok_return_value = pedigree_ok(the_pedigree_stats, max_self_agmr, max_self_R, max_ok_d, max_ok_z);
	   if(pedigree_ok_return_value <= 0){ // if pedigree is 'bad', look for alternatives
	     Viaxh* pwi_wrt_A = construct_viaxh(the_genotypes_set->accessions->size); // pairwise info for just A
	     long n_gt_accessions = the_genotypes_set->accessions->size;
    
	     // Calculate xhgmr
	     for(long ii=0; ii<n_gt_accessions; ii++){
	       Accession* A2 = the_genotypes_set->accessions->a[ii];
	       if(ii != A_gtset_idx){
		 ND the_xhgmr = xhgmr(the_genotypes_set, A, A2, quick_xhgmr);
		 //n_xhgmrs_calculated++;
		 double dbl_xhgmr = (the_xhgmr.d > 0)? (double)the_xhgmr.n/the_xhgmr.d : NAN;
		 if(dbl_xhgmr <= max_xhgmr){ // xhgmr is small enough 
		   //n_xhgmrs_le_max++;
		   push_to_viaxh(pwi_wrt_A, ii, dbl_xhgmr);
		 }
	       }
	     }
	     
	     Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pwi_wrt_A, max_candidate_parents);
	     sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, o_stream); //, long_output_format, d_scale_factor);
	     //if(alt_pedigrees->size == 0) 
	     free_vpedigree(alt_pedigrees);
	     free_viaxh(pwi_wrt_A);
	   }
	 }else if(alternative_pedigrees_level == 3){  // calculate d, z, etc. for triples involving the candidate parents in pairwise_info[i]	
	   Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pairwise_info[A_gtset_idx], max_candidate_parents);
	   sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, o_stream); //, long_output_format, d_scale_factor);
	   //if(alt_pedigrees->size == 0) fprintf(o_stream, " this accession has no candidate parents.");
	   free_vpedigree(alt_pedigrees);
	 }
	 
	 fprintf(o_stream, "\n");
      
	 free(the_pedigree_stats);
       } //    
     } // end loop over pedigrees
   
     //free_vpedigree(alt_pedigrees);
     /**/
     // NOW DO THE accessions without pedigrees if requested

     //
    
   }else{ // no pedigree file, consider all accessions as possible parents
 
    // **************************************************
    // ***  evaluate parent1-parent2-progeny triples  ***
    // **************************************************
    initialization_time = hi_res_time();
    
    long count_accs_w_no_cand_parents = 0;
    long count_accs_w_too_many_cand_parents = 0;
    
    for(long i=0; i<n_gt_accessions; i++){ 
      Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents.
	
      // calculate d, z, etc. for triples involving the candidate parents in pairwise_info[i]	
      Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(prog, the_genotypes_set, pairwise_info[i], max_candidate_parents);
      fprintf(o_stream, "%s  %ld   ", prog->id->a, pairwise_info[i]->size); // output the progeny id
      sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, o_stream); //, long_output_format, d_scale_factor);
      //if(alt_pedigrees->size == 0) fprintf(o_stream, " this accession has no candidate parents.");
      fprintf(o_stream, "\n");
      free_vpedigree(alt_pedigrees);
    } // end loop over offspring accessions

    fprintf(stderr, "# time for triple calculation: %8.3f\n", hi_res_time() - initialization_time);
    fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr);
    fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
    fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	    max_candidate_parents, count_accs_w_too_many_cand_parents);
  }
   fprintf(stdout, "time after initialization: %6.3f\n", hi_res_time() - initialization_time);
  
  // ********************  cleanup  **************************
  fclose(o_stream);
  free_genotypesset(the_genotypes_set);
  free_vidxid(the_gt_vidxid);
  fprintf(stdout, "# Total time: %10.4lf sec.\n", hi_res_time() - t_begin_main);
}
// **********************************************************
// ********************  end of main  ***********************
// **********************************************************




// **********************************************************
// ******************  functions  ***************************
// **********************************************************

void print_usage_info(FILE* stream){
  fprintf(stream, "-i   -input <filename>               input filename (dosages matrix, required)\n");
  fprintf(stream, "-o   -output <filename>              output filename (default: find_parents.out) \n");
  fprintf(stream, "-p   -pedigree <filename>            pedigree filename (optional)\n");
  fprintf(stream, "-m   -marker_max_missing_data <f>    don't use markers with missing data fraction > f. (default: %4.2f)\n", DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION);
  fprintf(stream, "-a   -accession_max_missing_data <f> accessions with missing data fraction >f are ignored > f (default: %4.2f)\n", DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION);
  fprintf(stream, "-f   -maf_min <f>                    don't use markers with minor allele frequency < f (default: %4.2f)\n", DEFAULT_MIN_MAF);
  fprintf(stream, "-x   -xhgmr_max <f>                  candidate parents must have xhgmr <= f (default: %4.2f)\n", DEFAULT_MAX_XHGMR);
  fprintf(stream, "-c   -candidate_parents_max <n>      sort candidate parents by xhgmr and use only best n (default: %4d)\n", DEFAULT_MAX_CANDIDATE_PARENTS);
  fprintf(stream, "-h   -help                           print this message and exit.\n");
}

Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents,
				double max_xhgmr, bool quick_xhgmr){
  long n_gt_accessions = the_genotypes_set->accessions->size;
  Viaxh** pairwise_info = (Viaxh**)malloc(n_gt_accessions*sizeof(Viaxh*)); // array of Viaxh*. pairwise_info[i]->a[j].xhgmr 
  for(long ii=0; ii<n_gt_accessions; ii++){
    pairwise_info[ii] = construct_viaxh(2*max_candidate_parents); // vector to hold candidate parents of accession with index ii
  }

  double t0 = hi_res_time();
  if(0){ // Calculate xhgmr
    calculate_xhgmrs(the_genotypes_set, pairwise_info, quick_xhgmr, max_xhgmr);
     fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0);
  } else{  // hgmr - bitwise
    calculate_hgmrs(the_genotypes_set, pairwise_info, max_xhgmr);
    fprintf(stdout, "# time for hgmrs: %10.3f \n", hi_res_time() - t0);
  }
  return pairwise_info;
}

void sort_and_output_pedigrees(GenotypesSet* the_gtsset, Vpedigree* the_pedigrees,
			       long max_solns_out, FILE* o_stream){ //, bool long_output_format, double d_scale_factor){
  // *******************************************************************************
  // ***  sort alt_pedigrees_array, output best parent pairs for this accession  ***
  // *******************************************************************************
  if(the_pedigrees->size > 0){
    if(the_pedigrees->size > 1) sort_vpedigree_by_d(the_pedigrees); // sort_vpedigree_by_maxdz(the_pedigrees);
      long n_out = (max_solns_out < the_pedigrees->size)? max_solns_out : the_pedigrees->size;
      for(long iii=0; iii < n_out; iii++){ // output the best n_out solutions
	Pedigree* a_pedigree = the_pedigrees->a[iii];
	
	a_pedigree->pedigree_stats->X3 = count_crossovers(the_gtsset, a_pedigree->F, a_pedigree->M, a_pedigree->A);
	fprintf(o_stream, "A  "); // to indicate alternative parents (not those from pedigree file)
	 print_pedigree_normalized(o_stream, the_pedigrees->a[iii]);
	 print_Xover_rates(o_stream, a_pedigree->pedigree_stats->X3);
      } // loop over solutions for one progeny accession
  }else{
    fprintf(o_stream, " this accession has no candidate parents.");
  }
} // end of sort_and_output_pedigrees

/* void set_scaled_d_in_vpedigree(Vpedigree* the_pedigrees, double d_scale_factor){ */
/*   for(long i = 0; i < the_pedigrees->size; i++){ */
/*     Pedigree_stats* the_ps = the_pedigrees->a[i]->pedigree_stats; */
/*     set_scaled_d_in_one_pedigree(the_ps, d_scale_factor); */
/*   } */
/* } */

/* void set_scaled_d_in_one_pedigree(Pedigree_stats* the_ps, double d_scale_factor){ */
/*   if(the_ps->d.d == 0){ */
/*     the_ps->scaled_d = NAN; */
/*   }else{ */
/*     double d = n_over_d(the_ps->d); */
/*     the_ps->scaled_d = d_scale_factor*d; */
/*   } */
/* } */

void print_X3_info(FILE* o_stream, Xcounts_3 X3){
  if(X3.XFA.Xa < X3.XFA.Xb){
    fprintf(o_stream, " %ld %ld %ld ", X3.XFA.Xa, X3.XFA.Xb, X3.XFA.Nhet);
  }else{
    fprintf(o_stream, " %ld %ld %ld ", X3.XFA.Xb, X3.XFA.Xa, X3.XFA.Nhet);
  }
  if(X3.XMA.Xa < X3.XMA.Xb){
    fprintf(o_stream, " %ld %ld %ld ", X3.XMA.Xa, X3.XMA.Xb, X3.XMA.Nhet);
  }else{
    fprintf(o_stream, " %ld %ld %ld ", X3.XMA.Xb, X3.XMA.Xa, X3.XMA.Nhet);
  }
  fprintf(o_stream, " %ld %ld %ld %ld ", X3.XFmin_3, X3.XFmax_3, X3.XMmin_3, X3.XMmax_3);
}

void random_set_means(GenotypesSet* gtset, long sample_size){
  // double* mean_hgmr, double* mean_R, double* mean_d, double* mean_z){
  
  double mean_hgmr = 0;
  double mean_R = 0;
  double mean_d = 0;
  double mean_z = 0;

  for(long j=0; j<sample_size; ){
    long i = (long)(rand()*(double)gtset->accessions->size/((double)RAND_MAX+1) );
    Accession* O = gtset->accessions->a[i];
    long ip1, ip2;
    while(1){
      ip1 = (long)(rand()*(double)gtset->accessions->size/((double)RAND_MAX+1) );
      if(ip1 != i) break;
    }
    Accession* P1 = gtset->accessions->a[ip1];
    while(1){
    ip2 = (long)(rand()*(double)gtset->accessions->size/((double)RAND_MAX+1) );
    if(ip2 != i) break;
    }
    Accession* P2 = gtset->accessions->a[ip2];
    
    Pedigree_stats* ps = bitwise_triple_counts(P1, P2, O);
    if(ps->d.d > 0  &&  ps->z.d > 0){
      mean_d += n_over_d(ps->d);
      mean_z += n_over_d(ps->z);
      mean_hgmr += n_over_d(ps->par1_hgmr) + n_over_d(ps->par2_hgmr);
      mean_R += n_over_d(ps->par1_R) + n_over_d(ps->par2_R);
      j++;
  }
  }
  mean_hgmr /= (double)(2*sample_size);
  mean_R /= (double)(2*sample_size);
  mean_d /= (double)sample_size;
  mean_z /= (double)sample_size;
  gtset->mean_hgmr = mean_hgmr;
  gtset->mean_R = mean_R;
  gtset->mean_d = mean_d;
  gtset->mean_z = mean_z;

  /* Vaccessions* accs = gtset->accessions; */
  /* for(long i=0; i< accs->size; i++){ */
  /*   Accession* A = accs->a[i]; */
  /*   A->hgmr1_n = A->par1_hgmr/mean_hgmr; */
  /*   A->hgmr2_n = A->par2_hgmr/mean_hgmr; */
  /*   A->R1_n = A->par1_R/mean_R; */
  /*   A->R2_n = A->par2_R/mean_R; */
  /*   A->d_n = A->d/mean_d; */
  /*   A->z_n = A->z/mean_z; */
  /* } */
}

void print_Xover2_rates(FILE* fh, Xcounts_2 X2){
  print_double_nan_as_hyphen(fh, X2.Xa/(double)X2.Nhet);
}

void print_Xover_rates(FILE* fh, Xcounts_3 X3){
  double Fnhet = X3.XFA.Nhet;
  double Mnhet = X3.XMA.Nhet;
  double XFA_rate = X3.XFA.Xa/Fnhet;
  double XMA_rate = X3.XMA.Xa/Mnhet;
  print_double_nan_as_hyphen(fh, XFA_rate);
  print_double_nan_as_hyphen(fh, XMA_rate);

  if(isnan(XFA_rate)  ||  isnan(XMA_rate)){ // if either parent is unknown, one of these will be nan.
    fprintf(fh, " - - - ");
  }else{
    print_double_nan_as_hyphen(fh, X3.XFmin_3/Fnhet);
    print_double_nan_as_hyphen(fh, X3.XMmin_3/Mnhet);
    print_double_nan_as_hyphen(fh, (X3.XFmin_3+X3.XMmin_3)/(Fnhet+Mnhet));
  }
}
