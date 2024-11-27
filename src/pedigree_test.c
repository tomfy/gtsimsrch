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
//#include "various.h"
//#define UNKNOWN -1
//#define DOSAGES 0
//#define GENOTYPES 1

int do_checks = 0; // option -c sets this to 1 to do some checks.

//double hi_res_time(void);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();
  
  int do_alternative_pedigrees = 0; // 0: none, 1: only when given pedigree is 'bad', 2: all
  // double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data_fraction = 0.25; // default; control this with -x command line option.
  double max_accession_missing_data_fraction = 1; //0.5;
  double min_minor_allele_frequency = 0; // 
    char* pedigree_test_output_filename = "pedigree_test_info";
    char* genotypes_matrix_output_filename = "genotype_matrix_out";
    
    double max_self_agmr12 = 1; // need to specify if doing alternative pedigrees 
    double max_ok_hgmr = 1; // accept everything as ok
    double max_self_r = 1; // need to specify if doing alternative pedigrees
    double max_ok_z = 1;
    double max_ok_d = 1; // accept everything as ok
    bool verbose = true;
    
    double ploidy = 2;
    double epsilon = 0.01;
    long Nthreads = 0;

    char* lpabpba_outfile = "logpabpba";
  
    // ***** process command line *****
    if (argc < 2) {
      fprintf(stderr, "Usage:  %s -g <genotypes_file>  -p <pedigree_file> \n", argv[0]);
      exit(EXIT_FAILURE);
    }
    //  char* dosages_filename = NULL;
    char* genotypes_filename = NULL;
    FILE *g_stream = NULL;
    char* pedigrees_filename = NULL;
    FILE *p_stream = NULL;

    // A: pedigree alternatives 0, 1, or 2
    // c: do checks,
    // d: dosages filename, g: genotypes filename (must have either d or g)
    // p: pedigree filename,
    // x: max fraction of missing data for markers,
    // o: output filename.
    // a: max 'self' agmr, h: max ok hgmr,  r: max 'self' r, D: max ok d;
    int c;
    //  int genotype_file_type = UNKNOWN;
    while((c = getopt(argc, argv, "A:cd:g:p:x:o:a:h:r:D:")) != -1){
      switch(c){
      case 'A':
	if(optarg == 0){
	  perror("option A requires an integer argument; 0, 1, or >=2\n");
	  exit(EXIT_FAILURE);
	}else{
	  do_alternative_pedigrees = atoi(optarg);
	  if(do_alternative_pedigrees < 0) exit(EXIT_FAILURE);
	  if(do_alternative_pedigrees > 2) do_alternative_pedigrees = 2;
	}
	break;
      case 'c':
	do_checks = 1;
	break;
      case 'd':
	genotypes_filename = optarg;
	g_stream = fopen(genotypes_filename, "r");
	if(g_stream == NULL){
	  fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	  exit(EXIT_FAILURE);
	}
	fclose(g_stream);
	//  genotype_file_type = DOSAGES;
	break;
      case 'g':
	genotypes_filename = optarg;
	g_stream = fopen(genotypes_filename, "r");
	if(g_stream == NULL){
	  fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	  exit(EXIT_FAILURE);
	}
	fclose(g_stream);
	//  genotype_file_type = GENOTYPES;
	break;
      case 'p':
	pedigrees_filename = optarg;
	p_stream = fopen(pedigrees_filename, "r");
	if(p_stream == NULL){
	  fprintf(stderr, "Failed to open %s for reading.\n", pedigrees_filename);
	  exit(EXIT_FAILURE);
	}
	break;
      case 'o':
	pedigree_test_output_filename = optarg;
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
      case 'a': // agmr12 < this -> looks like self
	if(optarg == 0){
	  perror("option a requires a numerical argument > 0\n");
	  exit(EXIT_FAILURE);
	}else{
	  max_self_agmr12 = atof(optarg);
	  if (max_self_agmr12 < 0) exit(EXIT_FAILURE);
	}
      case 'h': // hgmr > this -> poor parent-offspring candidate
	if(optarg == 0){
	  perror("option h requires a numerical argument > 0\n");
	  exit(EXIT_FAILURE);
	}else{
	  max_ok_hgmr = atof(optarg);
	  if (max_ok_hgmr < 0) exit(EXIT_FAILURE);
	}
	break;
      case 'r': // r > this -> probably biparental
	if(optarg == 0){
	  perror("option x requires a numerical argument > 0\n");
	  exit(EXIT_FAILURE);
	}else{
	  max_self_r = atof(optarg);
	  if (max_self_r < 0) exit(EXIT_FAILURE);
	}
	break;
      case 'D': // d > this argument means this a is poor candidate triple of parents and offspring
	if(optarg == 0){
	  perror("option x requires a numerical argument > 0\n");
	  exit(EXIT_FAILURE);
	}else{
	  max_ok_d = atof(optarg);
	  if (max_ok_d < 0) exit(EXIT_FAILURE);
	}
	break;
      case '?':
	printf("? case in command line processing switch.\n");
	if ((optopt == 'g') || (optopt == 'p') || (optopt == 'd') || (optopt == 'x') || (optopt == 'o'))
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
    if(pedigrees_filename == NULL){
      perror("must specify pedigrees filename: -i <filename>");
      exit(EXIT_FAILURE);
    }

    // exit(EXIT_FAILURE);

    FILE* lpabpba_stream = fopen(lpabpba_outfile, "w");
  
    FILE *o_stream = NULL;
    o_stream = fopen(pedigree_test_output_filename, "w");
    if(o_stream == NULL){
      fprintf(stderr, "Failed to open %s for writing.\n", pedigree_test_output_filename);
      exit(EXIT_FAILURE);
    }
    // fprintf(stderr, "# genotypes file type: %d\n", genotype_file_type);
    //  char* geno_file_type = (genotype_file_type == DOSAGES)? "dosages" : "genotypes";
    fprintf(stdout, "# Genotypes filename: %s max marker missing data: %5.2lf%c\n", genotypes_filename, max_marker_missing_data_fraction*100, '%');
    fprintf(stdout, "# Pedigrees filename: %s, output filename: %s \n", pedigrees_filename, pedigree_test_output_filename);
    fprintf(o_stream, "# Genotypes filename: %s max marker missing data: %5.2lf%c\n", genotypes_filename, max_marker_missing_data_fraction*100, '%');
    fprintf(o_stream, "# Pedigrees filename: %s, output filename: %s \n", pedigrees_filename, pedigree_test_output_filename);
  
    // *****  done processing command line  *****

    // ***************  read the genotypes file  *******************************
    double t_a = hi_res_time();
    GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
    add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set, max_accession_missing_data_fraction, Nthreads); // load the new set of accessions
    double t_b = hi_res_time();
       fprintf(stdout, "# Done reading genotypes file. %ld accessions:\n", the_genotypes_set->n_accessions);
     fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
     filter_genotypesset(the_genotypes_set, o_stream);
    double t_c = hi_res_time();
     fprintf(stdout, "# Time to clean genotype data: %6.3f sec.\n", t_c - t_b);
    rectify_markers(the_genotypes_set);
    store_homozygs(the_genotypes_set);
    populate_marker_dosage_counts(the_genotypes_set);
    set_n_00_1_22_1s(the_genotypes_set);
  check_genotypesset(the_genotypes_set);
  set_Abits_Bbits(the_genotypes_set, Nthreads);
    double t_d = hi_res_time();
    fprintf(stdout, "# Time to rectify genotype data: %6.3f sec.\n", t_d - t_c);
  
    /* if(0){ */
    /*   double t0 = hi_res_time(); */
    /*   long n_acc = the_genotypes_set->accessions->size; */
    /*   for(long ii=0; ii<n_acc; ii++){ */
    /* 	Accession* A1 = the_genotypes_set->accessions->a[ii]; */
    /* 	for(long jj=ii+1; jj<n_acc; jj++){ */
    /* 	  Accession* A2 = the_genotypes_set->accessions->a[jj]; */
    /* 	  ND the_xhgmr = */
    /* 	    xhgmr(the_genotypes_set, A1, A2, 0); */
    /* 	  //	  ghgmr(the_genotypes_set, A1, A2); */
    /* 	  if(1 && the_xhgmr.d > 0){ */
    /* 	    double dbl_xhgmr = (double)the_xhgmr.n/the_xhgmr.d; */
    /* 	    if(dbl_xhgmr < 0.5){ */
    /* 	      fprintf(stderr, "%ld %ld   %ld %ld\n", */
    /* 		      // A1->id->a, A2->id->a, */
    /* 		      ii, jj,  */
    /* 		      the_xhgmr.n, the_xhgmr.d); // , (double)the_xhgmr.n/the_xhgmr.d); */
    /* 	    } */
    /* 	  }  */
    /* 	} */
    /*   } */
    /*   /\* quick_and_dirty_hgmrs(the_genotypes_set);   *\/ */
    /*   fprintf(stdout, "# time for xhgmrs: %10.3f \n", hi_res_time() - t0);  */
    /*   exit(0); */
    /* } */
    double t_start = hi_res_time();

    

    // ***************  read the pedigrees file  ***************************
    fprintf(stderr, "# before read_and_store_pedigrees_3col\n"); 
    Vidxid* the_vidxid = construct_sorted_vidxid(the_genotypes_set->accessions); // ids and indexes of the_genotypes_set, sorted by id
    const Vpedigree* pedigrees = read_and_store_pedigrees_3col(p_stream, the_vidxid, the_genotypes_set);  
    fprintf(stderr, "# after read_and_store_pedigrees_3col\n");
    fclose(p_stream);
    fprintf(stdout, "# Done reading pedigree file. Time to read pedigree file: %6.3f\n", hi_res_time() - t_start);
    fprintf(stdout, "# Stored genotypes of %ld accessions, and  %ld pedigrees. \n", the_genotypes_set->accessions->size, pedigrees->size);

    // ***************  Done reading input files  ******************************
  
    const Vlong* parent_idxs = accessions_with_offspring(pedigrees, the_genotypes_set->accessions->size); // , the_genotypes_set->n_accessions);
    fprintf(stdout, "# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);
    fprintf(stdout, "# Cumulative time so far: %6.3f sec.\n", hi_res_time() - t_begin_main);
    
    t_start = hi_res_time();
    char ploidy_char = (char)(the_genotypes_set->ploidy + 48);
    fprintf(stderr, "# Analyzing %ld pedigrees from table. \n", pedigrees->size);
    for(long i=0; i<pedigrees->size; i++){
      if(i % 100  == 0) fprintf(stdout, "# Done testing %ld pedigrees.\n", i);
      
      Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(pedigrees->a[i], the_genotypes_set, 1.0); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
      // assert(strcmp(pedigrees->a[i]->Accession->id, pedigrees->a[i]->A->id->a) == 0);
      Accession* A = pedigrees->a[i]->A;
      Accession* F = pedigrees->a[i]->F;
      Accession* M = pedigrees->a[i]->M;
   
      if(0 || (F != NULL  ||  M != NULL)){ // at least one parent specified in pedigree

	/* if(A != NULL  &&  F != NULL  && M != NULL){ */
	/*   two_doubles lpabpba_AF = logPABPBA(the_genotypes_set, A, F); */
	/*   two_doubles lpabpba_AM = logPABPBA(the_genotypes_set, A, M); */
	/*   fprintf(lpabpba_stream, "A: %s  F: %s  %8.5f  %8.5f   M: %s  %8.5f  %8.5f\n", */
	/* 	  A->id->a, F->id->a, lpabpba_AF.x1, lpabpba_AF.x2, M->id->a, lpabpba_AM.x1, lpabpba_AM.x2); */
	/* } */
	
	fprintf(o_stream, "%20s %5ld %20s %20s %ld  ",
		A->id->a, A->missing_data_count,
		(F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA", the_pedigree_stats->all_good_count);

	print_pedigree_stats(o_stream, the_pedigree_stats, verbose); // print the stats for this pedigree

	if(do_alternative_pedigrees == 1){
	  fprintf(stderr, "doing alt pedigrees for %s\n", A->id->a);
	  if(pedigree_ok(the_pedigree_stats, max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d) == 0){	    
	    Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	    print_pedigree_alternatives(o_stream, alt_pedigrees, 4, verbose);
	  }
	}else if((do_alternative_pedigrees == 2)){
	  // || (pedigree_ok(the_pedigree_stats, max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d) == 0)){
	  //   fprintf(stderr, "About to get alt. pedigrees. max_ok_hgmr: %8.4lf,  max_ok_d: %8.4lf \n", max_ok_hgmr, max_ok_d);
	  Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	  print_pedigree_alternatives(o_stream, alt_pedigrees, 4, verbose);
	}
	free(the_pedigree_stats);
	fprintf(o_stream, "\n");   
      } //
    
    } // end loop over pedigrees
  
    fprintf(stdout, "# Done testing %ld pedigrees. Time to test pedigrees %6.3f sec.\n", pedigrees->size, hi_res_time() - t_start);

    // ********************  cleanup  **************************
    t_start = hi_res_time();
    fclose(o_stream);
    free_genotypesset(the_genotypes_set);
    //  fprintf(stderr, "after free_genotypesset\n");
    free_vidxid(the_vidxid);
    free_vpedigree(pedigrees);
    //  fprintf(stderr, "after free_vpedigree\n");
    free_vlong(parent_idxs);
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


