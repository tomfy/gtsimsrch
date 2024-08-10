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



// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();

  
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
  bool sort_by_z = false; // if false, sort by d
  bool use_xhgmr = true;
    
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
  FILE *g_stream = NULL;

  // c: do checks,
  // i: dosages input filename
  // o: output filename.
  // x: max fraction of missing data for markers,
  // f: max fraction of missing data for accessions
  // a: min minor allele frequency
  // h: max xhgmr for candidate parents
  // p: max_candidate_parents

  int c;
  while(1){
    int option_index = 0;
    static struct option long_options[] = {
      {"input",   required_argument, 0,  'i'}, // filename of new data set
      {"output",  required_argument, 0,  'o'}, // output filename
      {"marker_max_missing_data", required_argument, 0, 'm'}, // markers with > this fraction missing data will not be used.
      {"maf_min", required_argument, 0, 'f'}, //
      {"accession_max_missing_data", required_argument, 0, 'a'},
      {"candidate_parents_max", required_argument, 0, 'p'},
      {"xhgmr_max", required_argument, 0, 'x'},
      {"help", no_argument, 0, 'h'},
      {"bitwise", no_argument, 0, 'b'},
      {"sort_by_z", no_argument, 0, 's'},
      {0,         0,                 0,  0 }
    };
     
    c = getopt_long_only(argc, argv, "", long_options, &option_index);
    if(c == -1) break;
    switch(c){

  
      // while((c = getopt(argc, argv, "ci:x:o:m:h:r:D:p:")) != -1){
      // switch(c){

    case 'c':
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
    case 'o':
      output_filename = optarg;
      break;
    case 'p':
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
	fprintf(stderr, "option x (max_xhgmr) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'b':
      bitwise = true;
      break;
       case 's':
      sort_by_z = true;
      break;
    case 'h':
      print_usage_info(stderr);
      exit(EXIT_FAILURE);
      break;
     
      /* case 'D': // d > this argument means this a is poor candidate triple of parents and offspring */
      /* 	if(optarg == 0){ */
      /* 	  perror("option x requires a numerical argument > 0\n"); */
      /* 	  exit(EXIT_FAILURE); */
      /* 	}else{ */
      /* 	  max_ok_d = atof(optarg); */
      /* 	  if (max_ok_d < 0) exit(EXIT_FAILURE); */
      /* 	} */
      /* 	break; */
    case '?':
      fprintf(stderr, "? case in command line processing switch.\n");
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr, "default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
    // printf("optind: %d argc: %d\n", optind, argc);
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

  
  // ***************  read the genotypes file  *******************************
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
  double t_d = hi_res_time();
  fprintf(stdout, "# Time for other initial processing of genotype data: %6.3f sec.\n", t_d - t_c);
  fflush(stdout);
  long n_accessions = the_genotypes_set->accessions->size;

  // ****************************************************************
  // ***  Calculate pairwise quantities (agmr, hgmr, xhgmr?, R?)  ***
  // ****************************************************************
  double t0 = hi_res_time();
  
  Viaxh** pairwise_info = (Viaxh**)malloc(n_accessions*sizeof(Viaxh*)); // array of Viaxh*. pairwise_info[i]->a[j].xhgmr 
  for(long ii=0; ii<n_accessions; ii++){
    pairwise_info[ii] = construct_viaxh(2*max_candidate_parents); // vector to hold candidate parents of accession with index ii
  }
    
  if(use_xhgmr){ // Calculate xhgmr
    calculate_xhgmrs(the_genotypes_set, pairwise_info, quick_xhgmr, max_xhgmr);
    fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0);
  } else{  // hgmr - bitwise
    calculate_hgmrs(the_genotypes_set, pairwise_info, max_hgmr);
    fprintf(stdout, "# time for hgmrs: %10.3f \n", hi_res_time() - t0);
  }
  // ***************************************
  // ***  end of pairwise calculation  *****
  // ***************************************
 

  // **************************************************
  // ***  evaluate parent1-parent2-progeny triples  ***
  // **************************************************

  double t_x = hi_res_time();

  long count_accs_w_no_cand_parents = 0;
  long count_accs_w_too_many_cand_parents = 0;

  Vpedigree* pedigrees = construct_vpedigree(1000);
  for(long i=0; i<n_accessions; i++){ // progeny accessions
    Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents.
    Viaxh* cppps = pairwise_info[i]; // cppps: candidate parent progeny pairs. These are the indices and xhgmrs (and agmrs, hgmrs?) of candidate parents to accession 'prog',. 
    long ncandpairs = cppps->size;
    if(ncandpairs == 0){ 
      count_accs_w_no_cand_parents++;
      fprintf(o_stream, "# %s has no candidate parents (xhmgr <= %8.5f)\n", prog->id->a, max_xhgmr);
    }else if(ncandpairs > max_candidate_parents){ // if too many candidates, just take the max_candidate_parents best ones
      sort_viaxh_by_xhgmr(cppps);
      cppps->size = max_candidate_parents; // limit candidate to max_candidate_parents
      count_accs_w_too_many_cand_parents++;
      // fprintf(stderr, "# %ld %7.4f  %ld %7.4f \n", cppps->a[0]->l, cppps->a[0]->d, cppps->a[1]->l, cppps->a[1]->d); 
    }

    calculate_triples_for_one_accession(prog, the_genotypes_set, cppps, pedigrees);

    // ***************************************************************************
    // ***  sort pedigrees_array, output best parent pairs for this accession  ***
    // ***************************************************************************
    if(pedigrees->size > 1) sort_vpedigree_by_maxdz(pedigrees);
    long n_out = (max_solns_out < pedigrees->size)? max_solns_out : pedigrees->size;
    //  fprintf(stderr, "n_out: %ld \n", n_out);
    for(long iii=0; iii < n_out; iii++){
      Accession* par1 = pedigrees->a[iii]->F;
      Accession* par2 = pedigrees->a[iii]->M;
      Pedigree_stats* the_ps = pedigrees->a[iii]->pedigree_stats;
      if((!multiple_solns_on_one_line) || (iii == 0)) fprintf(o_stream, "%s  %ld   ", prog->id->a, ncandpairs); // output progeny id
      fprintf(o_stream, "%s %s  ", par1->id->a, par2->id->a); // output candidate pair of parents
      print_pedigree_stats(o_stream, the_ps);
      if((!multiple_solns_on_one_line) || (iii == n_out-1)) fprintf(o_stream, "\n");
    } // loop over solutions for one progeny accession
    
    pedigrees->size = 0;
  } // end loop over offspring accessions
  
  fprintf(stderr, "# time for triple calculation: %8.3f\n", hi_res_time() - t_x);
  fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr);
  fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
  fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	  max_candidate_parents, count_accs_w_too_many_cand_parents);
  
  // ********************  cleanup  **************************
  fclose(o_stream);
  free_genotypesset(the_genotypes_set);
  //  fprintf(stderr, "after free_genotypesset\n");
  //   free_vidxid(the_vidxid);
  free_vpedigree(pedigrees);
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



