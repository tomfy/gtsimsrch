// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "gtset.h"
#include "pedigree.h"

int do_checks = 0; // option -c sets this to 1 to do some checks.

double hi_res_time(void);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();

  
  double max_marker_missing_data_fraction = 0.2; // default; control this with -x command line option.
  double max_accession_missing_data_fraction = 0.5;
  double min_minor_allele_frequency = 0; // 
  char* output_filename = "find_parents.out";
  double max_xhgmr = 0.15;
  long max_candidate_parents = 100;
    
  double ploidy = 2;
  //double epsilon = 0.01;
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -i <dosages_file> [-o <output_filename> -h <max_xhgmr>] \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* genotypes_filename = NULL;
  FILE *g_stream = NULL;

  // c: do checks,
  // i: dosages input filename 
  // x: max fraction of missing data for markers,
  // o: output filename.
  // m: min minor allele frequency
  // h: max xhgmr for candidate parents

  int c;
  while((c = getopt(argc, argv, "ci:x:o:m:h:r:D:")) != -1){
    switch(c){

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
    case 'x':
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_marker_missing_data_fraction = atof(optarg);
	if (max_marker_missing_data_fraction < 0){
	  printf("# max missing data fraction in markers will be set automatically.\n");
	}}
      break;
    case 'm': 
      min_minor_allele_frequency = (double)atof(optarg);
      if(min_minor_allele_frequency < 0){
	fprintf(stderr, "option m (min_minor_allele_frequency) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'h': // xhgmr > this -> poor parent-offspring candidate
      if(optarg == 0){
	perror("option h requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_xhgmr = atof(optarg);
	if (max_xhgmr < 0) exit(EXIT_FAILURE);
	fprintf(stderr, "# max xhgmr: %7.4f\n", max_xhgmr);
      }
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
      printf("? case in command line processing switch.\n");
      if ((optopt == 'g') || (optopt == 'x') || (optopt == 'o'))
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf(stderr, "Unknown option character: %d\n", optopt);
      exit(EXIT_FAILURE);
    default:
      perror("default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
    // printf("optind: %d argc: %d\n", optind, argc);
  if(optind < argc){
    perror("Non-option arguments. Bye.\n");
    exit(EXIT_FAILURE);
  }
  if(genotypes_filename == NULL){
    perror("must specify genotype filename: -i <filename>");
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
  clean_genotypesset(the_genotypes_set);
  double t_c = hi_res_time();
  fprintf(stdout, "# Time to clean genotype data: %6.3f sec.\n", t_c - t_b);
  rectify_markers(the_genotypes_set);
  store_homozygs(the_genotypes_set);
  populate_marker_dosage_counts(the_genotypes_set);
  double t_d = hi_res_time();
  fprintf(stdout, "# Time to rectify genotype data: %6.3f sec.\n", t_d - t_c);
  fflush(stdout);
  
  // ***********  get xhgmr for all pairs:  ******************
  double t0 = hi_res_time();
  long n_xhgmrs_le_max = 0;
  long n_acc = the_genotypes_set->accessions->size;

  Vlong** cand_pppairs = (Vlong**)malloc(n_acc*sizeof(Vlong*));
  for(long ii=0; ii<n_acc; ii++){
    cand_pppairs[ii] = construct_vlong(20);
  }
  for(long ii=0; ii<n_acc; ii++){
    if(ii % 100  == 0) fprintf(stderr, "# ii: %ld\n", ii);
    Accession* A1 = the_genotypes_set->accessions->a[ii];
    for(long jj=ii+1; jj<n_acc; jj++){
      Accession* A2 = the_genotypes_set->accessions->a[jj];
      ND the_xhgmr = xhgmr(the_genotypes_set, A1, A2, 1); // quick version
      if(the_xhgmr.d > 0){
	double dbl_xhgmr = (double)the_xhgmr.n/the_xhgmr.d;
	if(dbl_xhgmr <= max_xhgmr){
	  n_xhgmrs_le_max++;
	  push_vlong(cand_pppairs[ii], jj);
	  push_vlong(cand_pppairs[jj], ii);
	}
      } 
    }
  }

  fprintf(stdout, "# n xhgmrs <= %8.5f :  %ld\n", max_xhgmr, n_xhgmrs_le_max);
  fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0); 
  //  ********************************************************
  
  long count_accs_w_no_cand_parents = 0;
  long count_accs_w_too_many_cand_parents = 0;
  for(long i=0; i<n_acc; i++){
    Accession* prog = the_genotypes_set->accessions->a[i];
    Vlong* cppps = cand_pppairs[i];
    long ncandpairs = cppps->size;
    //   fprintf(stdout, "%ld  %20s  %ld \n", i, prog->id->a, ncandpairs);
    if(ncandpairs == 0){
      count_accs_w_no_cand_parents++;
    }else if(ncandpairs <= max_candidate_parents){	  
      for(long ii=0; ii<cppps->size; ii++){
	long par1idx = cppps->a[ii];
	Accession* par1 = the_genotypes_set->accessions->a[par1idx];
	for(long jj=ii; jj<cppps->size; jj++){
	  long par2idx = cppps->a[jj];
	  Accession* par2 = the_genotypes_set->accessions->a[par2idx];
	  Pedigree_stats* the_ps = triple_counts( par1->genotypes->a,  par2->genotypes->a, prog->genotypes->a, ploidy );
	  the_ps->xhgmr1 = xhgmr(the_genotypes_set, par1, prog, 0); // do full (not 'quick') xhgmr
	  the_ps->xhgmr2 = xhgmr(the_genotypes_set, par2, prog, 0); // do full (not 'quick') xhgmr
	  fprintf(o_stream, "%s %s %s  %ld  ", prog->id->a, par1->id->a, par2->id->a, ncandpairs);
	  print_pedigree_stats(o_stream, the_ps); fprintf(o_stream, "\n");

	}
      }
    }else{
      count_accs_w_too_many_cand_parents++;
    }
  }
  fprintf(o_stream, "# candidate parents have xghmr <= %8.5f\n", max_xhgmr);
  fprintf(o_stream, "# number of accessions with no candidate parents found: %ld\n", count_accs_w_no_cand_parents);
  fprintf(o_stream, "# number of accessions with > %ld candidate parents found: %ld\n",
	  max_candidate_parents, count_accs_w_too_many_cand_parents);
      
  // ********************  cleanup  **************************
  fclose(o_stream);
  free_genotypesset(the_genotypes_set);
  //  fprintf(stderr, "after free_genotypesset\n");
  //   free_vidxid(the_vidxid);
  //   free_vpedigree(pedigrees);
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

