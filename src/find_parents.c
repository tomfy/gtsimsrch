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

int do_checks = 0; // option -c sets this to 1 to do some checks.

double hi_res_time(void);

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
  double max_xhgmr = 0.15;
  long max_candidate_parents = 80;
  bool quick_xhgmr = true;
  bool oldway = false;
    
  double ploidy = 2;
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

  // exit(EXIT_FAILURE);
  
  FILE *o_stream = NULL;
  o_stream = fopen(output_filename, "w");
  if(o_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
    exit(EXIT_FAILURE);
  }
  // fprintf(stderr, "# genotypes file type: %d\n", genotype_file_type);
  //  char* geno_file_type = (genotype_file_type == DOSAGES)? "dosages" : "genotypes";
  fprintf(stdout, "# Genotypes filename: %s max marker missing data: %5.2lf%c\n", genotypes_filename, max_marker_missing_data_fraction*100, '%');
  fprintf(o_stream, "# Genotypes filename: %s max marker missing data: %5.2lf%c\n", genotypes_filename, max_marker_missing_data_fraction*100, '%');
   
  
  // *****  done processing command line  *****

  // ***************  read the genotypes file  *******************************
  double t_a = hi_res_time();
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
  add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set, max_accession_missing_data_fraction); // load the new set of accessions
  double t_b = hi_res_time();
  fprintf(stdout, "# Done reading genotypes file. %ld accessions will be analyzed.\n", the_genotypes_set->n_accessions);
  fprintf(stdout, "# %ld accessions were excluded from analysis due to > %5.2lf percent missing data.\n",
	  the_genotypes_set->n_bad_accessions, max_accession_missing_data_fraction*100);
  fprintf(o_stream, "# %ld accessions were excluded from analysis due to > %5.2lf percent missing data.\n",
	  the_genotypes_set->n_bad_accessions, max_accession_missing_data_fraction*100);
  fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
  filter_genotypesset(the_genotypes_set, o_stream);
  double t_c = hi_res_time();
  fprintf(stdout, "# Time to clean genotype data: %6.3f sec.\n", t_c - t_b);
  rectify_markers(the_genotypes_set);
  store_homozygs(the_genotypes_set);
  populate_marker_dosage_counts(the_genotypes_set);
  check_genotypesset(the_genotypes_set);
  set_Abits_Bbits(the_genotypes_set);
  double t_d = hi_res_time();
  fprintf(stdout, "# Time to rectify genotype data: %6.3f sec.\n", t_d - t_c);
  fflush(stdout);

  /* if(0){ // old way */
  /*   // ***********  get xhgmr for all pairs:  ****************** */
  /*   double t0 = hi_res_time(); */
  /*   long n_xhgmrs_le_max = 0; */
  /*   long n_acc = the_genotypes_set->accessions->size; */

  /*   Vlong** cand_pppairs = (Vlong**)malloc(n_acc*sizeof(Vlong*)); */
  /*   for(long ii=0; ii<n_acc; ii++){ */
  /*     cand_pppairs[ii] = construct_vlong(20); */
  /*   } */
  /*   for(long ii=0; ii<n_acc; ii++){ */
  /*     if(ii % 100  == 0) fprintf(stderr, "# ii: %ld\n", ii); */
  /*     Accession* A1 = the_genotypes_set->accessions->a[ii]; */
  /*     for(long jj=ii+1; jj<n_acc; jj++){ */
  /* 	Accession* A2 = the_genotypes_set->accessions->a[jj]; */
  /* 	ND the_xhgmr = xhgmr(the_genotypes_set, A1, A2, quick_xhgmr); */
  /* 	if(the_xhgmr.d > 0){ */
  /* 	  double dbl_xhgmr = (double)the_xhgmr.n/the_xhgmr.d; */
  /* 	  if(dbl_xhgmr <= max_xhgmr){ */
  /* 	    n_xhgmrs_le_max++; */
  /* 	    push_to_vlong(cand_pppairs[ii], jj); */
  /* 	    push_to_vlong(cand_pppairs[jj], ii); */
  /* 	  } */
  /* 	}  */
  /*     } */
  /*   } */
  /*   fprintf(stdout, "# n xhgmrs <= %8.5f :  %ld\n", max_xhgmr, n_xhgmrs_le_max); */
  /*   fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0);  */
  /*   //  ******************************************************** */

  /*   //  ****************   */
  /*   long count_accs_w_no_cand_parents = 0; */
  /*   long count_accs_w_too_many_cand_parents = 0; */
  /*   for(long i=0; i<n_acc; i++){ */
  /*     Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents. */
  /*     Vlong* cppps = cand_pppairs[i]; // these are the indices of candidate parents to accession 'prog'.  */
  /*     long ncandpairs = cppps->size; */
  /*     //   fprintf(stdout, "%ld  %20s  %ld \n", i, prog->id->a, ncandpairs); */
  /*     if(ncandpairs == 0){ */
  /* 	count_accs_w_no_cand_parents++; */
  /* 	fprintf(o_stream, "### %s has no candidate parents (xhmgr <= %8.5f)\n", prog->id->a, max_xhgmr); */
  /*     }else if(ncandpairs <= max_candidate_parents){	   */
  /* 	for(long ii=0; ii<cppps->size; ii++){ */
  /* 	  long par1idx = cppps->a[ii]; */
  /* 	  Accession* par1 = the_genotypes_set->accessions->a[par1idx]; */
  /* 	  for(long jj=ii; jj<cppps->size; jj++){ */
  /* 	    long par2idx = cppps->a[jj]; */
  /* 	    Accession* par2 = the_genotypes_set->accessions->a[par2idx]; */
  /* 	    Pedigree_stats* the_ps = triple_counts( par1->genotypes->a,  par2->genotypes->a, prog->genotypes->a, ploidy ); */
  /* 	    the_ps->xhgmr1 = xhgmr(the_genotypes_set, par1, prog, 0); // do full (not 'quick') xhgmr */
  /* 	    the_ps->xhgmr2 = xhgmr(the_genotypes_set, par2, prog, 0); // do full (not 'quick') xhgmr */
  /* 	    fprintf(o_stream, "%s %s %s  %ld  ", prog->id->a, par1->id->a, par2->id->a, ncandpairs); */
  /* 	    print_pedigree_stats(o_stream, the_ps); */
  /* 	    long N_22 = the_ps->d_22.n; // n_00_2 + n_22_0 */
  /* 	    long D_22 = the_ps->d_22.d; // n_00_x + n_22_x */
  /* 	    long N_21 = the_ps->d_21.n; // n_02_0 + n_02_2 + n_20_0 + n_20_2 + n_00_1 + n_22_1 */
  /* 	    long D_21 = the_ps->d_21.d; // n_00_x + n_02_x + n_20_x + n_22 */
  /* 	    // 'apparent' style dissimilarity: */
  /* 	    ND nd_app = {N_22 + 0.5*N_21, D_21}; */
  /* 	    print_d_r(o_stream, nd_app); // col 24 */
  /* 	    print_d_r(o_stream, the_ps->d_22); // col 26 */
  /* 	    print_d_r(o_stream, the_ps->d_21); // col 28 */
  /* 	    print_d_r(o_stream, the_ps->d_11); // col 30 */
  /* 	    ND xxx = {N_22 + N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d}; */
  /* 	    print_d_r(o_stream, xxx); */
  /* 	    ND xxxx = {N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d}; */
  /* 	    print_d_r(o_stream, xxxx); */
  /* 	    fprintf(o_stream, "\n"); */

  /* 	  } */
  /* 	} */
  /*     }else{ */
  /* 	count_accs_w_too_many_cand_parents++; */
  /*     } */
  /*   } */
  /*   fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr); */
  /*   fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents); */
  /*   fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n", */
  /* 	    max_candidate_parents, count_accs_w_too_many_cand_parents); */
  /* }else{ */
  
  // new way   
  double t0 = hi_res_time();
  long n_xhgmrs_calculated = 0;
  long n_xhgmrs_le_max = 0;
  long n_acc = the_genotypes_set->accessions->size;

  Vld** progeny_cplds = (Vld**)malloc(n_acc*sizeof(Vld*)); // candidate parent
  for(long ii=0; ii<n_acc; ii++){
    progeny_cplds[ii] = construct_vld(2*max_candidate_parents);
  }
  for(long ii=0; ii<n_acc; ii++){
    if(ii % 100  == 0) fprintf(stderr, "# ii: %ld\n", ii);
    Accession* A1 = the_genotypes_set->accessions->a[ii];
    for(long jj=ii+1; jj<n_acc; jj++){
      Accession* A2 = the_genotypes_set->accessions->a[jj];
      ND the_xhgmr = {0, 0};
      if(0){
	the_xhgmr = xhgmr(the_genotypes_set, A1, A2, quick_xhgmr);
	n_xhgmrs_calculated++;
	if(the_xhgmr.d > 0){
	  double dbl_xhgmr = (double)the_xhgmr.n/the_xhgmr.d;
	  if(dbl_xhgmr <= max_xhgmr){
	    n_xhgmrs_le_max++;
	    // push_to_vlong(cand_pppairs[ii], jj);
	    //  push_to_vlong(cand_pppairs[jj], ii);
	    push_to_vld(progeny_cplds[ii], jj, dbl_xhgmr);
	    push_to_vld(progeny_cplds[jj], ii, dbl_xhgmr);
	  }
	}
      }else{ // new way - bitwise
	four_longs bfcs = bitwise_agmr_hgmr(A1, A2);
	long b_agmr_num = bfcs.l1 + bfcs.l3;
	long b_agmr_denom = b_agmr_num + bfcs.l2 + bfcs.l4;
	long b_hgmr_num = bfcs.l1;
	long b_hgmr_denom = bfcs.l1 + bfcs.l2;
	double true_dist = (b_agmr_denom > 0)? (double)(b_agmr_num + bfcs.l1)/(double)b_agmr_denom : -1;
	double agmr = (b_agmr_denom > 0)? (double)b_agmr_num / (double)b_agmr_denom : -1;
	double hgmr = (b_hgmr_denom > 0)? (double)b_hgmr_num / (double)b_hgmr_denom : -1;
	n_xhgmrs_calculated++;

	if(hgmr < max_xhgmr){
	  n_xhgmrs_le_max++;
	  push_to_vld(progeny_cplds[ii], jj, hgmr);
	  push_to_vld(progeny_cplds[jj], ii, hgmr);
	}
      }
    }
  }
  fprintf(stdout, "# n xhgmrs calculated: %ld ;  <= %8.5f :  %ld\n", n_xhgmrs_calculated, max_xhgmr, n_xhgmrs_le_max);
  fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0);

  long count_accs_w_no_cand_parents = 0;
  long count_accs_w_too_many_cand_parents = 0;
  exit(0);
  // long capacity = 10000;
  //  Pedigree_stats** pedstats_array = (Pedigree_stats**)malloc(capacity*sizeof(Pedigree_stats*));
  // long n_ppairs_stored = 0;

  Vpedigree* pedigrees = construct_vpedigree(1000);
    
  for(long i=0; i<n_acc; i++){ // progeny accessions
    Accession* prog = the_genotypes_set->accessions->a[i]; // the progeny accession, for which we seek parents.
    Vld* cppps = progeny_cplds[i]; // these are the indices of candidate parents to accession 'prog'. 
    long ncandpairs = cppps->size;
    //   fprintf(stdout, "%ld  %20s  %ld \n", i, prog->id->a, ncandpairs);
    if(ncandpairs == 0){
      count_accs_w_no_cand_parents++;
      fprintf(o_stream, "# %s has no candidate parents (xhmgr <= %8.5f)\n", prog->id->a, max_xhgmr);
    }else if(ncandpairs > max_candidate_parents){ // if too many candidates, just take the max_candidate_parents best ones
      sort_vld_by_d(cppps);
      cppps->size = max_candidate_parents;
      count_accs_w_too_many_cand_parents++;
      fprintf(stderr, "# %ld %7.4f  %ld %7.4f \n", cppps->a[0]->l, cppps->a[0]->d, cppps->a[1]->l, cppps->a[1]->d); 
    }
      
    for(long ii=0; ii<cppps->size; ii++){
      long par1idx = cppps->a[ii]->l;
      Accession* par1 = the_genotypes_set->accessions->a[par1idx];
      for(long jj=ii; jj<cppps->size; jj++){
	long par2idx = cppps->a[jj]->l;
	Accession* par2 = the_genotypes_set->accessions->a[par2idx];
	Pedigree_stats* the_ps = triple_counts( par1->genotypes->a,  par2->genotypes->a, prog->genotypes->a, ploidy );
	Pedigree* the_pedigree = construct_pedigree(prog, par1, par2);
	the_pedigree->pedigree_stats = the_ps;
	the_ps->xhgmr1 = xhgmr(the_genotypes_set, par1, prog, 0); // do full (not 'quick') xhgmr
	the_ps->xhgmr2 = xhgmr(the_genotypes_set, par2, prog, 0); // do full (not 'quick') xhgmr
	if(oldway){
	  
	  fprintf(o_stream, "%s %s %s  %ld  ", prog->id->a, par1->id->a, par2->id->a, ncandpairs);
	  print_pedigree_stats(o_stream, the_ps);
	  long N_22 = the_ps->d_22.n; // n_00_2 + n_22_0
	  long D_22 = the_ps->d_22.d; // n_00_x + n_22_x
	  long N_21 = the_ps->d_21.n; // n_02_0 + n_02_2 + n_20_0 + n_20_2 + n_00_1 + n_22_1
	  long D_21 = the_ps->d_21.d; // n_00_x + n_02_x + n_20_x + n_22
	  // 'apparent' style dissimilarity:
	  ND nd_app = {N_22 + 0.5*N_21, D_21};
	  print_d_r(o_stream, nd_app); // col 24
	  print_d_r(o_stream, the_ps->d_22); // col 26
	  print_d_r(o_stream, the_ps->d_21); // col 28
	  print_d_r(o_stream, the_ps->d_11); // col 30
	  ND xxx = {N_22 + N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d};
	  print_d_r(o_stream, xxx);
	  ND xxxx = {N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d};
	  print_d_r(o_stream, xxxx);
	  fprintf(o_stream, "\n");
	}else{ // new way
	  /* if(n_ppairs_stored >= capacity){ */
	  /*   capacity *= 2; */
	  /*   pedstats_array = (Pedigree_stats**)realloc(pedstats_array, capacity*sizeof(Pedigree_stats*)); */
	  /* } */
	  //   pedstats_array[n_ppairs_stored] = the_ps; // the pedigree_stats for this pair of parents
	  //   n_ppairs_stored++;
	  push_to_vpedigree(pedigrees, the_pedigree);
	}
      } // end loop over parent 2
    } // end loop over parent 1
    if(!oldway){ // sort pedstats_array, output best parent pairs for this offspring accession
      //	fprintf(stderr, "### n_ppairs_stored: %ld\n", n_ppairs_stored);
      if(pedigrees->size > 0){
	//	sort_pedigree_stats_by_d(pedstats_array, n_ppairs_stored);
	sort_vpedigree_by_d(pedigrees);
	/* if(pedigrees->size >= 2){ */
	/* Pedigree* p0 = pedigrees->a[0]; */
	/* Pedigree_stats* ps0 = pedigrees->a[0]->pedigree_stats; */
	/* double d0 = (double)ps0->d.n/ps0->d.d; */
	/* Pedigree* p1 = pedigrees->a[1]; */
	/* Pedigree_stats* ps1 = pedigrees->a[1]->pedigree_stats; */
	/* double d1 = (double)ps1->d.n/ps1->d.d; */
	/* fprintf(stderr, "# %s  %s %s  %8.5f  %s  %s %s  %8.5f\n", */
	/* 	p0->A->id->a, p0->F->id->a, p0->M->id->a, d0, */
	/* 	p1->A->id->a, p1->F->id->a, p1->M->id->a, d1); */
	/* } */
	for(long iii=0; iii<pedigrees->size; iii++){
	  Accession* par1 = pedigrees->a[iii]->F;
	  Accession* par2 = pedigrees->a[iii]->M;
	  Pedigree_stats* the_ps = pedigrees->a[iii]->pedigree_stats;
	  if(iii == 0) fprintf(o_stream, "%s ", prog->id->a);
	  if(iii<2){
	    fprintf(o_stream, "%s %s  %ld  ", par1->id->a, par2->id->a, ncandpairs);
	    print_pedigree_stats(o_stream, the_ps);
	    long N_22 = the_ps->d_22.n; // n_00_2 + n_22_0
	    long D_22 = the_ps->d_22.d; // n_00_x + n_22_x
	    long N_21 = the_ps->d_21.n; // n_02_0 + n_02_2 + n_20_0 + n_20_2 + n_00_1 + n_22_1
	    long D_21 = the_ps->d_21.d; // n_00_x + n_02_x + n_20_x + n_22
	    // 'apparent' style dissimilarity:
	    ND nd_app = {N_22 + 0.5*N_21, D_21};
	    print_d_r(o_stream, nd_app); // col 24
	    print_d_r(o_stream, the_ps->d_22); // col 26
	    print_d_r(o_stream, the_ps->d_21); // col 28
	    print_d_r(o_stream, the_ps->d_11); // col 30
	    ND xxx = {N_22 + N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d};
	    print_d_r(o_stream, xxx);
	    ND xxxx = {N_21 + the_ps->d_11.n, D_21 + the_ps->d_11.d};
	    print_d_r(o_stream, xxxx);
	    fprintf(o_stream, "  ");
	  } 
	  // free(pedstats_array[iii]);
	}
	fprintf(o_stream, "\n");
	//n_ppairs_stored = 0;
      }
    }
    pedigrees->size = 0;
  } // end loop over offspring accessions
  fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr);
  fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
  fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	  max_candidate_parents, count_accs_w_too_many_cand_parents);
  //}
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

double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

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
