// C program to test pedigrees, search for alternatives using genotype data.
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
#include <sys/types.h>
#include <sys/sysinfo.h>

#include "gtset.h"
#include "pedigree.h"

#define DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION  0.25
#define DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION  1.0
#define DEFAULT_MIN_MAF  0.05
#define DEFAULT_MAX_NHGMR 0.2
#define DEFAULT_MAX_CANDIDATE_PARENTS 1000
#define DEFAULT_RANDOM_SAMPLE_SIZE 40000
#define DEFAULT_MAX_SOLNS_OUT 3  // in addition to pedigree from file.
#define ALL_ALTS 1 // search for alt. pedigrees for all accessions. all accessions are considered as parents.
#define NO_ALTS 0 // just check the pedigrees
#define BADPED_ALTS 2 // only search for alternatives if pedigree looks bad
#define PEDPAR_ALTS 3  // only accessions which appear as parents in the pedigree file are considered as parents


int do_checks = 0; // option -c sets this to 1 to do some checks.
struct sysinfo memInfo;

//double hi_res_time(void);


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVMP(){ //Note: this value is in KB! virtual mem this proc.
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPMP(){ //Note: this value is in KB! physical mem this proc
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}
void print_mem_info(FILE* stream);

void print_usage_info(FILE* stream);
Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents, double max_nhgmr, bool quick_xhgmr);
void sort_and_output_pedigrees(GenotypesSet* the_gtsset, Vpedigree* the_pedigrees, long max_solns_out, bool do_phased, FILE* o_stream); //, bool long_output_format, double d_scale_factor);
void set_scaled_d_in_one_pedigree(Pedigree_stats* the_ps, double d_scale_factor);
void set_scaled_d_in_vpedigree(Vpedigree* the_pedigrees, double d_scale_factor);
// void print_Xover2_rates(FILE* fh, Xcounts_2 X2);
void print_X3_info(FILE* o_stream, Xcounts_3 X3);
void print_crossover_info(FILE* o_stream, Xcounts_3 X3);
void print_Xover_rates(FILE* fh, Xcounts_3 X3);
// void d_z_of_random_set(GenotypesSet* gtset, long sample_size, double* mean_d, double* mean_z);
void random_set_means(GenotypesSet* gtset, long sample_size);
Viaxh* hgmrs_wrt_one_accession(GenotypesSet* the_genotypes_set, Accession* A, double max_nhgmr);
void print_dummy_pedigree_output(FILE* stream, bool phased, bool do_phased); // for when an accession has no pedigree

// double* mean_hgmr, double* mean_R, double* mean_d, double* mean_z);
// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();

  long alternative_pedigrees_level = NO_ALTS; // find alt. pedigrees for  0: none, 1: accessions with bad pedigrees, 2: all accessions
  // 1,2: all accessions are considered are considered as possible parents.
  // 3,4: same as 1,2, respectively, but only accessions which are parents in pedigree file are considered as possible parents (not implemented)
  double max_marker_missing_data_fraction = DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION; // default; control this with -x command line option.
  double max_accession_missing_data_fraction = DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION;
  double min_minor_allele_frequency = DEFAULT_MIN_MAF; // 
  char* output_filename = "find_parents.out";
  double max_nhgmr = DEFAULT_MAX_NHGMR;
  long max_candidate_parents = DEFAULT_MAX_CANDIDATE_PARENTS;
  bool do_phased = false;
  bool quick_xhgmr = true;
  long max_solns_out = DEFAULT_MAX_SOLNS_OUT; // in addition to pedigree
  bool multiple_solns_on_one_line = true;
  
  // the following used in categorizing pedigrees from file as good or bad:
  double max_self_agmr = 0.03;
  double max_ok_hgmr = 0.16;
  double max_self_R = 0.1;
  double max_ok_d = 0.15; // normalized

  bool long_output_format = false; // true; true -> output denominators as well as 
  
  double ploidy = 2;
  
  long nprocs = (long)get_nprocs(); // returns 2*number of cores if hyperthreading.
  long Nthreads = -1; // multithreading not implemented properly yet.  Fix threaded_input in gtset.c // (nprocs > 2)? nprocs/2 : 1; // default number of threads

  unsigned rand_seed = time(0); // (unsigned)tspec.tv_nsec;
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -in <dosages_file> [-out <output_filename> -nhgmr_max <max_nhgmr>] \n", argv[0]);
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
      {"nhgmr_max", required_argument, 0, 'n'},
      {"help", no_argument, 0, 'h'},
      {"check", no_argument, 0, 'k'}, // not implemented
      {"alternative_pedigrees_level", required_argument, 0, 'l'},
      {"seed", required_argument, 0, 'r'},
      {"hgmr_max", required_argument, 0, 'H'},
      {"agmr_self_max", required_argument, 0, 'S'},
      {"R_self_max", required_argument, 0, 'R'},
      {"d_max", required_argument, 0, 'D'},
      {"threads", required_argument, 0, 't'},
      {"do_phased", no_argument, 0, 'P'},
      {0,         0,                 0,  0 }
    };
     
    c = getopt_long_only(argc, argv, "", long_options, &option_index);
    if(c == -1) break;
    switch(c){
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
    case 't':
      Nthreads = atoi(optarg);
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
     case 'r': 
      rand_seed  = (unsigned)atoi(optarg);
      break;
    case 'n': 
      max_nhgmr = (double)atof(optarg);
      if(max_nhgmr < 0){
	fprintf(stderr, "option nhgmr_max requires an real argument >= 0\n");
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
    case 'l':
      alternative_pedigrees_level = (long)atoi(optarg);
      
      if(alternative_pedigrees_level < 0  ||  alternative_pedigrees_level > 3){
	fprintf(stderr, "option alternative_pedigrees requires one of 0, 1, 2, or 3 as argument.\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'P':
      do_phased = true;
      break;
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
  fprintf(stdout, "# Genotypes filename: %s\n", genotypes_filename);

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
  fprintf(stderr, "After add_accessions_...\n");
  print_vchar(stdout, the_genotypes_set->acc_filter_info);
  print_vchar(o_stream, the_genotypes_set->acc_filter_info);
  double t_b = hi_res_time();
  fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
  
  filter_genotypesset(the_genotypes_set);
  // fprintf(stderr, "after filter_genotypesset\n");
  print_vchar(stdout, the_genotypes_set->marker_filter_info);
  print_vchar(o_stream, the_genotypes_set->marker_filter_info);
  fprintf(stdout, "# Data is phased? %s\n", the_genotypes_set->phased? "true" : "false");
  fprintf(stdout, "# Phased analysis and output requested? %s\n", do_phased? "true" : "false");

  double t_c = hi_res_time();
  fprintf(stdout, "# Time to filter genotype data: %6.3f sec.\n", t_c - t_b);
  // rectify_markers(the_genotypes_set); // only needed for xhgmr to be fast - not using xhgmr
  // store_homozygs(the_genotypes_set); // needed?
  // fprintf(stderr, "before set_chromosome_start_indices.  \n");
  if(the_genotypes_set->phased)    set_chromosome_start_indices(the_genotypes_set);
  // fprintf(stderr, "phased: %c \n", the_genotypes_set->phased? 't' : 'f');
  // populate_marker_dosage_counts(the_genotypes_set); // needed?

  check_genotypesset(the_genotypes_set);
  set_Abits_Bbits(the_genotypes_set, Nthreads);
  Vidxid* the_gt_vidxid = construct_sorted_vidxid(the_genotypes_set->accessions); // ids and indexes of the_genotypes_set, sorted by id

  // print_vidxid(stderr, the_gt_vidxid);
  double t_d = hi_res_time();
  fprintf(stdout, "# Time for other initial processing of genotype data: %6.3f sec.\n", t_d - t_c);

  long sample_size = DEFAULT_RANDOM_SAMPLE_SIZE;
  random_set_means(the_genotypes_set, sample_size); // get and store mean hgmr, R, etc.
  fprintf(stderr, "# mean hgmr, R, d, z: %8.5f  %8.5f  %8.5f  %8.5f\n",
	  the_genotypes_set->mean_hgmr, the_genotypes_set->mean_R, the_genotypes_set->mean_d, the_genotypes_set->mean_z);
  fflush(stdout);

  double t_e = hi_res_time();
  long n_gt_accessions = the_genotypes_set->accessions->size;
  fprintf(stdout, "# Time to calculate means for %ld  random pairs.: %6.3f sec.\n", sample_size, t_e - t_d);
 
  if(p_stream != NULL){ // have pedigree file
    
    // ********************************************************* 
    // ***  Read the pedigrees file  ***************************
    // *********************************************************
    // fprintf(stderr, "About to read_and_store_pedigrees...\n");
    const Vpedigree* the_pedigrees = read_and_store_pedigrees_3col(p_stream, the_gt_vidxid, the_genotypes_set); // acc id and then female, male parent ids in first 3 cols.
    fclose(p_stream); 
    fprintf(stdout, "# Done reading pedigree file. Time to read pedigree file: %6.3f sec.\n", hi_res_time() - t_d);
    fprintf(stdout, "# Stored genotypes of %ld accessions, and  %ld pedigrees. \n", the_genotypes_set->accessions->size, the_pedigrees->size); 
    // *********************************************************
    // ***  Done reading in pedigrees from file  ***************
    // ********************************************************* 
  double initialization_time = hi_res_time();
  fprintf(stdout, "# total initialization time: %6.3f sec.\n", initialization_time - t_begin_main);

    // ****************************************************************
    // ***  If -alt 1 (i.e. ALL_ALTS)  option  ********************
    // ***  Calculate hgmrs for all pairs  ****************************
    // **************************************************************** 
    Viaxh** pairwise_info = NULL;
    fprintf(stderr, "# alternative pedigrees level: %ld\n", alternative_pedigrees_level);
    if(alternative_pedigrees_level == ALL_ALTS){ // calculate hgmrs for all pairs
      fprintf(stderr, "# calculating hgmrs for all pairs\n");
      Vlong* all_acc_idxs = construct_vlong_whole_numbers(the_genotypes_set->accessions->size);
      pairwise_info = calculate_hgmrs(the_genotypes_set, all_acc_idxs, max_candidate_parents, max_nhgmr);
      fprintf(stdout, "# time to calculate hgmrs: %.5f\n", hi_res_time() - initialization_time);
    }else if(alternative_pedigrees_level == PEDPAR_ALTS){ // only calculate pairs with at least one acc appearing as parent in pedigree file.
      const Vlong* parent_idxs = accessions_with_offspring(the_pedigrees, the_genotypes_set->accessions->size); // , the_genotypes_set->n_gt_accessions);
      fprintf(stdout, "# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);
      pairwise_info = calculate_hgmrs(the_genotypes_set, parent_idxs, max_candidate_parents, max_nhgmr);
      fprintf(stdout, "# time to calculate hgmrs: %.5f\n", hi_res_time() - initialization_time);
    }
    // ***************************************
    // ***  end of pairwise calculation  *****
    // *************************************** 
  
    double t_f = hi_res_time(); 
    fprintf(stdout, "# Cumulative time so far: %6.3f sec.\n", t_f - t_begin_main);
    fprintf(stderr, "# Analyzing %ld pedigrees from file. \n", the_pedigrees->size);
    //fprintf(stderr, "# Phased? %s\n", the_genotypes_set->phased? "true" : "false");
    Vaccession* the_accessions = the_genotypes_set->accessions;
    for(long i=0; i<the_accessions->size; i++){ 
      if(i % 500  == 0) fprintf(stdout, "# Done testing %ld pedigrees.\n", i);
      Accession* A = the_accessions->a[i]; // the_pedigree->A;
      Accession* F = (A->Fpar_idx != ID_NA_INDEX)? the_accessions->a[A->Fpar_idx]: NULL; //the_pedigree->F;
      Accession* M = (A->Mpar_idx != ID_NA_INDEX)? the_accessions->a[A->Mpar_idx] : NULL; //the_pedigree->M;     
     
      if(0 || (F != NULL  ||  M != NULL)){ // at least one parent is specified in pedigree
	Pedigree* the_pedigree = construct_pedigree(A, F, M); // the_pedigrees->a[i];
	Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(the_pedigree, the_genotypes_set); 
	the_pedigree->pedigree_stats = the_pedigree_stats;
	
	fprintf(o_stream, "%s\tP\t", A->id->a); // progeny accession and 'P' to indicate these are parents from the pedigree file.
	print_pedigree_normalized(o_stream, the_pedigree); //, the_genotypes_set);
	// two_doubles hratios = heterozyg_ratios(A, F);
	if(the_genotypes_set->phased && do_phased){
	  Xcounts_3 X3 = count_crossovers(the_genotypes_set, F, M, A);
	  //	  print_Xover_rates(o_stream, X3);
	  //      print_X3_info(o_stream, X3);
	  print_crossover_info(o_stream, X3);
	} // end of if phased branch
	if(alternative_pedigrees_level == PEDPAR_ALTS){ // search for parents, considering only accessions in parent_idxs as possible parents
	    Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pairwise_info[A->index], max_candidate_parents);
	    sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, do_phased, o_stream);
	    free_viaxh(pairwise_info[A->index]);
	    free_vpedigree(alt_pedigrees);	  
	}else if(alternative_pedigrees_level == BADPED_ALTS){ // iff pedigres bad, search for parents. Consider all accessions as possible parents
	  bool pedigree_ok =  d_ok(the_pedigree_stats, max_ok_d*the_genotypes_set->mean_hgmr);
	  // pedigree_ok(the_pedigree_stats, max_self_agmr, max_self_R, max_ok_d);
	  if(! pedigree_ok){ // if pedigree is 'bad', look for alternatives
	    Viaxh* hgmrs_wrt_A = hgmrs_wrt_one_accession(the_genotypes_set, A, max_nhgmr);	     
	    Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, hgmrs_wrt_A, max_candidate_parents);
	    sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, do_phased, o_stream); //, long_output_format, d_scale_factor);
	    free_vpedigree(alt_pedigrees);

	    free_viaxh(hgmrs_wrt_A);
	  }
	}else if(alternative_pedigrees_level == ALL_ALTS){  // calculate d etc. for triples involving the candidate parents in pairwise_info[i]	
	  Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pairwise_info[A->index], max_candidate_parents);
	  if(i%100 == 0){ fprintf(stdout, "%ld  ", i); print_mem_info(stdout); }
	  sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, do_phased, o_stream);
	  free_viaxh(pairwise_info[A->index]);
	  free_vpedigree(alt_pedigrees);
	}
	fprintf(o_stream, "\n");
	fflush(o_stream);
	// free(the_pedigree_stats);
	free(the_pedigree); // _stats); 
	//    end of accession has pedigree branch
      }else{  // accession has no pedigree
	if( (alternative_pedigrees_level == PEDPAR_ALTS)  ||  (alternative_pedigrees_level == ALL_ALTS ) ){
	  // Accession* prog =  A; // the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents
      	  //  fprintf(o_stream, "%s  P - -  - - - - - - -  - - - - - - - ", A->id->a); // dummy output for pedigree
	  
	  fprintf(o_stream, "%s\t", A->id->a);
	  print_dummy_pedigree_output(o_stream, the_genotypes_set->phased, do_phased);
	  Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(A, the_genotypes_set, pairwise_info[i], max_candidate_parents);
	  sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, do_phased, o_stream);
	  fprintf(o_stream, "\n");
	  free_vpedigree(alt_pedigrees);
	   fprintf(stdout, "%ld  ", i); print_mem_info(stdout);
	}
      }
      
    } // end loop over accessions
      fprintf(stdout, "# time for triple (FTR) calculation: %8.3f\n", hi_res_time() - t_f);
  }else{ //  *******  no pedigree file, consider all accessions as possible parents  ********
    Viaxh** pairwise_info = calculate_hgmrs(the_genotypes_set,
					    construct_vlong_whole_numbers(the_genotypes_set->accessions->size),
					    max_candidate_parents, max_nhgmr);
    fprintf(stdout, "# time to calculate hgmrs: %.5f\n", hi_res_time() - t_e);
    
    double triples_to_calculate = 0;
    for(long i=0; i<the_genotypes_set->accessions->size; i++){
      long n_parents = pairwise_info[i]->size;
      if(n_parents > max_candidate_parents) n_parents = max_candidate_parents;
      triples_to_calculate += n_parents*(n_parents-1)/2;
    }
    fprintf(stderr, "Number of triples to calculate: %8.4g\n", triples_to_calculate);  
  
    // **************************************************
    // ***  evaluate parent1-parent2-progeny triples  ***
    // **************************************************
    double t_g = hi_res_time();
    
    long count_accs_w_no_cand_parents = 0;
    long count_accs_w_too_many_cand_parents = 0;

    double triples_calculated = 0;
    for(long i=0; i<n_gt_accessions; i++){
      long n_parents = pairwise_info[i]->size;
      if(n_parents == 0){
	count_accs_w_no_cand_parents++;
      }else if(n_parents > max_candidate_parents){
	count_accs_w_too_many_cand_parents++;
	n_parents = max_candidate_parents;
      }
      Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents.	
      // calculate d, z, etc. for triples involving the candidate parents in pairwise_info[i]
      Vpedigree* alt_pedigrees = calculate_triples_for_one_accession(prog, the_genotypes_set, pairwise_info[i], max_candidate_parents);
      triples_calculated += n_parents*(n_parents-1)/2;
       if(i%100 == 0){
	// fprintf(stdout, "%ld  ", i); print_mem_info(stdout);
	 
	 fprintf(stderr, "#%lld  %8.5g %% done\n",(long long)triples_calculated, 100*triples_calculated/triples_to_calculate);
      }
      // output
      fprintf(o_stream, "%s\t", prog->id->a); // output the progeny id
      print_dummy_pedigree_output(o_stream, the_genotypes_set->phased, do_phased);
      sort_and_output_pedigrees(the_genotypes_set, alt_pedigrees, max_solns_out, do_phased, o_stream);
      fprintf(o_stream, "\n");
      
      free_viaxh(pairwise_info[i]);
      free_vpedigree(alt_pedigrees);
    } // end loop over offspring accessions

    fprintf(stderr, "# time for triple (FTR) calculation: %8.3f\n", hi_res_time() - t_g);
    fprintf(o_stream, "# candidate parents have nhgmr <= %8.5f\n", max_nhgmr);
    fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
    fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	    max_candidate_parents, count_accs_w_too_many_cand_parents);
  }
  
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

Viaxh* hgmrs_wrt_one_accession(GenotypesSet* the_genotypes_set, Accession* A, double max_hgmr_norm){
  Viaxh* hgmrs_wrt_A = construct_viaxh(the_genotypes_set->accessions->size); 
  for(long ii=0; ii < the_genotypes_set->accessions->size; ii++){
    Accession* A2 = the_genotypes_set->accessions->a[ii];
    if(A2->index != A->index){
      ND hgmr_nd = bitwise_hgmr(A, A2);
      double dbl_hgmr_norm = n_over_d(hgmr_nd)/the_genotypes_set->mean_hgmr;
      if(dbl_hgmr_norm <= max_hgmr_norm){ // xhgmr is small enough 
	push_to_viaxh(hgmrs_wrt_A, ii, dbl_hgmr_norm);
      }
    }
  }
  return hgmrs_wrt_A;
}

void print_usage_info(FILE* stream){
  fprintf(stream, "-i   -input <filename>               input filename (dosages matrix, required)\n");
  fprintf(stream, "-o   -output <filename>              output filename (default: find_parents.out) \n");
  fprintf(stream, "-p   -pedigree <filename>            pedigree filename (optional)\n");
  fprintf(stream, "-m   -marker_max_missing_data <f>    don't use markers with missing data fraction > f. (default: %4.2f)\n", DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION);
  fprintf(stream, "-a   -accession_max_missing_data <f> accessions with missing data fraction >f are ignored > f (default: %4.2f)\n", DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION);
  fprintf(stream, "-f   -maf_min <f>                    don't use markers with minor allele frequency < f (default: %4.2f)\n", DEFAULT_MIN_MAF);
  fprintf(stream, "-l    alternative_pedigrees_level    controls search for alternative pedigrees  (default: 0 )\n");
  fprintf(stream, "                                            0: just check pedigrees, 1: search all accessions for alternative pedigrees\n");
  fprintf(stream, "-n   -nhgmr_max <f>                  candidate parents must have normalized hgmr <= f (default: %4.2f)\n", DEFAULT_MAX_NHGMR);
  fprintf(stream, "-c   -candidate_parents_max <n>      sort candidate parents by nhgmr and use only best n (default: %4d)\n", DEFAULT_MAX_CANDIDATE_PARENTS);
  fprintf(stream, "-r   -seed                           random number generator seed. (default set from clock).\n");
  fprintf(stream, "-h   -help                           print this message and exit.\n");
}



void sort_and_output_pedigrees(GenotypesSet* the_gtsset, Vpedigree* the_pedigrees,
			       long max_solns_out, bool do_phased, FILE* o_stream){ //, bool long_output_format, double d_scale_factor){
  // *******************************************************************************
  // ***  sort alt_pedigrees_array, output best parent pairs for this accession  ***
  // *******************************************************************************
  //fprintf(stderr, "#xx  %ld %ld\n", the_pedigrees->size, max_solns_out);
  if(the_pedigrees->size > 0){
    if(the_pedigrees->size > 1) sort_vpedigree_by_d(the_pedigrees); // 
    long n_out = (max_solns_out < the_pedigrees->size)? max_solns_out : the_pedigrees->size;
    //fprintf(stderr, "#  %ld  %ld\n", the_pedigrees->size, n_out);
    for(long iii=0; iii < n_out; iii++){ // output the best n_out solutions
      Pedigree* a_pedigree = the_pedigrees->a[iii];
      /* fprintf(stderr, "before count_crossovers.\n");
       fprintf(stderr, "A %s\n", a_pedigree->A->id->a);
      fprintf(stderr, "F %s\n", a_pedigree->F->id->a);
      fprintf(stderr, "M %s\n", a_pedigree->M->id->a); /* */
    

      fprintf(o_stream, "A\t"); // to indicate alternative parents (not those from pedigree file)
      print_pedigree_normalized(o_stream, the_pedigrees->a[iii]);
      if(the_gtsset->phased && do_phased){
	a_pedigree->pedigree_stats->X3 = count_crossovers(the_gtsset, a_pedigree->F, a_pedigree->M, a_pedigree->A);
	//fprintf(stderr, "after count_crossovers.\n");
	// print_Xover_rates(o_stream, a_pedigree->pedigree_stats->X3);
	// print_X3_info(o_stream, a_pedigree->pedigree_stats->X3);
	print_crossover_info(o_stream, a_pedigree->pedigree_stats->X3);

      }
      // fprintf(o_stream, " XXX "); 
    } // loop over solutions for one progeny accession
  }else{
    fprintf(o_stream, " this accession has no candidate parents.");
  }
} // end of sort_and_output_pedigrees

void print_dummy_pedigree_output(FILE* stream, bool phased, bool do_phased){ // for when an accession has no pedigree
  if(phased && do_phased){
    fprintf(stream, "P\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"); // phased  17 fields
  }else{
    fprintf(stream, "P\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"); // unphased  10 fields
  }
}

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

void print_Xover_rates(FILE* fh, Xcounts_3 X3){
  double Fnhet = X3.XFA.Nhet;
  double Mnhet = X3.XMA.Nhet;
  double XFA_rate = X3.XFA.Xmin/Fnhet; // N crossovers implied if F is parent of A
  double XMA_rate = X3.XMA.Xmin/Mnhet; // N crossovers implied if M is parent of A
  print_double_nan_as_hyphen(fh, XFA_rate); 
  print_double_nan_as_hyphen(fh, XMA_rate);

  // now output numbers of crossovers implied by F and M being the parents of A, and
  // each pair of chromosomes in A must derive, one from F and the other from M.
  if(isnan(XFA_rate)  ||  isnan(XMA_rate)){ // if either parent is unknown, one of these will be nan.
    fprintf(fh, "\t-\t-\t-");
  }else{
    print_double_nan_as_hyphen(fh, X3.XFmin_3/Fnhet);
    print_double_nan_as_hyphen(fh, X3.XMmin_3/Mnhet);
    print_double_nan_as_hyphen(fh, (X3.XFmin_3+X3.XMmin_3)/(Fnhet+Mnhet));
  }
}

void print_X3_info(FILE* o_stream, Xcounts_3 X3){
  /* if(X3.XFA.Xa < X3.XFA.Xb){
    fprintf(o_stream, " %ld %ld %ld ", X3.XFA.Xa, X3.XFA.Xb, X3.XFA.Nhet);
  }else{
    fprintf(o_stream, " %ld %ld %ld ", X3.XFA.Xb, X3.XFA.Xa, X3.XFA.Nhet);
  }
  if(X3.XMA.Xa < X3.XMA.Xb){
    fprintf(o_stream, " %ld %ld %ld ", X3.XMA.Xa, X3.XMA.Xb, X3.XMA.Nhet);
  }else{
    fprintf(o_stream, " %ld %ld %ld ", X3.XMA.Xb, X3.XMA.Xa, X3.XMA.Nhet);
  } /* */
  fprintf(o_stream, "%ld\t%ld\t%ld\t", X3.XFA.Xmin, X3.XFA.Xmax, X3.XFA.Nhet);
  fprintf(o_stream, "%ld\t%ld\t%ld\t", X3.XMA.Xmin, X3.XMA.Xmax, X3.XMA.Nhet);
  fprintf(o_stream, "%ld\t%ld\t%ld\t%ld\t", X3.XFmin_3, X3.XFmax_3, X3.XMmin_3, X3.XMmax_3);
}

void print_crossover_info(FILE* fh, Xcounts_3 X3){
  long Fnhet = X3.XFA.Nhet;
  long Mnhet = X3.XMA.Nhet;
  double XFA_rate = X3.XFA.Xmin/Fnhet; // N crossovers implied if F is parent of A
  double XMA_rate = X3.XMA.Xmin/Mnhet; // N crossovers implied if M is parent of A
  // fprintf(fh, " %ld %ld %ld ", X3.XFA.Xmin, X3.XF);

  fprintf(fh, "%ld\t%ld\t%ld\t", X3.XFA.Xmin, X3.XFmin_3, Fnhet);
  fprintf(fh, "%ld\t%ld\t%ld\t", X3.XMA.Xmin, X3.XMmin_3, Mnhet);

  /*  fprintf(stderr, " %ld\t%ld\t%ld\t%ld\t\n", X3.XFmin_3, Fnhet, X3.XMmin_3, Mnhet);
    long n =  (X3.XFmin_3+X3.XMmin_3);
    long d = Fnhet+Mnhet;
    fprintf(stderr, " %ld\t%ld\t\n", n, d);
    print_double_nan_as_hyphen(stderr,(double)n/d); /* */
  //print_double_nan_as_hyphen(fh, XFA_rate); 
  //print_double_nan_as_hyphen(fh, XMA_rate);

  // now output numbers of crossovers implied by F and M being the parents of A, and
  // each pair of chromosomes in A much derive, one from F and the other from M.
  if(Fnhet == 0  ||  Mnhet == 0){ // if either parent is unknown, one of these will be nan.
    fprintf(fh, "\t- ");
  }else{
    print_double_nan_as_hyphen(fh, (double)(X3.XFmin_3+X3.XMmin_3)/(Fnhet+Mnhet));
    //print_double_nan_as_hyphen(fh, (X3.XFmin_3+X3.XMmin_3)/(Fnhet+Mnhet));
  } 
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
}

void print_mem_info(FILE* stream){
	sysinfo(&memInfo);
	long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
	virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
	virtualMemUsed *= memInfo.mem_unit;
	long vmem_this_proc = getVMP();
	long pmem_this_proc = getPMP();
	
	fprintf(stream, "  %lld   %ld  %ld\n", virtualMemUsed, vmem_this_proc, pmem_this_proc);
}

//  **************  unused  *********************
/* Viaxh** calculate_pairwise_info(GenotypesSet* the_genotypes_set, long max_candidate_parents, */
/* 				double max_xhgmr, bool quick_xhgmr){ */
 
/*   Viaxh** pairwise_info = (Viaxh**)malloc(the_genotypes_set->accessions->size*sizeof(Viaxh*)); // array of Viaxh*. pairwise_info[i]->a[j].xhgmr  */
/*   for(long ii=0; ii< the_genotypes_set->accessions->size; ii++){ */
/*     pairwise_info[ii] = construct_viaxh(2*max_candidate_parents); // vector to hold candidate parents of accession with index ii */
/*   } */

/*   double t0 = hi_res_time(); */
/*   if(0){ // Calculate xhgmr */
/*     calculate_xhgmrs(the_genotypes_set, pairwise_info, quick_xhgmr, max_xhgmr); */
/*      fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0); */
/*   } else{  // hgmr - bitwise */
/*     calculate_hgmrs(the_genotypes_set, pairwise_info, max_xhgmr); */
/*     fprintf(stdout, "# time for hgmrs: %10.3f \n", hi_res_time() - t0); */
/*   } */
/*   return pairwise_info; */
/* } */

/* void print_Xover2_rates(FILE* fh, Xcounts_2 X2){ */
/*   print_double_nan_as_hyphen(fh, X2.Xa/(double)X2.Nhet); */
/* } */
