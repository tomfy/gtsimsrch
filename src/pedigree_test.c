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
//#define UNKNOWN -1
//#define DOSAGES 0
//#define GENOTYPES 1

int do_checks_flag = 0; // option -c sets this to 1 to do some checks.

double hi_res_time(void);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();
  
  int do_alternative_pedigrees = 0; // 0: none, 1: only when given pedigree is 'bad', 2: all
  // double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data_fraction = 0.2; // default; control this with -x command line option.
  double min_minor_allele_frequency = 0; // 
    char* pedigree_test_output_filename = "pedigree_test_info";
    char* genotypes_matrix_output_filename = "genotype_matrix_out";
    double max_self_agmr12 = 1; // need to specify if doing alternative pedigrees 
    double max_ok_hgmr = 1; // accept everything as ok
    double max_self_r = 1; // need to specify if doing alternative pedigrees
    double max_ok_z = 1;
    double max_ok_d = 1; // accept everything as ok
    
    double ploidy = 2;
    double epsilon = 0.01;
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
	do_checks_flag = 1;
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
    add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set); // load the new set of accessions
    double t_b = hi_res_time();
       fprintf(stdout, "# Done reading genotypes file. %ld accessions:\n", the_genotypes_set->n_accessions);
     fprintf(stdout, "# Time to read genotype data: %6.3f sec.\n", t_b - t_a);
    clean_genotypesset(the_genotypes_set);
    double t_c = hi_res_time();
     fprintf(stdout, "# Time to clean genotype data: %6.3f sec.\n", t_c - t_b);
    rectify_markers(the_genotypes_set);
    store_homozygs(the_genotypes_set);
    populate_marker_dosage_counts(the_genotypes_set);
    double t_d = hi_res_time();
    fprintf(stdout, "# Time to rectify genotype data: %6.3f sec.\n", t_d - t_c);
  
    if(0){
      double t0 = hi_res_time();
      quick_and_dirty_hgmrs(the_genotypes_set);  
      fprintf(stdout, "# time for hgmrs: %10.3f \n", hi_res_time() - t0);
      exit(0);
    }
    double t_start = hi_res_time();
 
    Vidxid* the_vidxid = construct_sorted_vidxid(the_genotypes_set);
    // fprintf(stderr, "# size of the_vidxid: %ld \n", the_vidxid->size);
    for(long i=0; i<the_vidxid->size; i++){
      long index = the_vidxid->a[i]->index;
      //   fprintf(stderr, "# i idx id: %ld %ld  %s\n", i, index, the_vidxid->a[i]->id); 
      Accession* a = the_genotypes_set->accessions->a[index];  
      //  fprintf(stderr, "# i, index, id, index: %ld %ld  %s  %ld\n", i, index, a->id->a, index_of_id_in_vidxid(the_vidxid, a->id->a));
    }
    // exit(0);
    // ***************  read the pedigrees file  ***************************
      const Vpedigree* pedigrees = read_and_store_pedigrees_3col(p_stream, the_vidxid, the_genotypes_set);  

      //  const Vpedigree* pedigrees = read_the_pedigrees_file_and_store(p_stream, the_vidxid, the_genotypes_set);
    fclose(p_stream);
    fprintf(stdout, "# Done reading pedigree file. Time to read pedigree file: %6.3f\n", hi_res_time() - t_start);
    fprintf(stdout, "# Stored genotypes of %ld accessions, and  %ld pedigrees. \n", the_genotypes_set->accessions->size, pedigrees->size);

    // ***************  Done reading input files  ******************************

  
    const Vlong* parent_idxs = accessions_with_offspring(pedigrees, the_genotypes_set->n_accessions);
    fprintf(stdout, "# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);
    fprintf(stdout, "# Cumulative time so far: %6.3f sec.\n", hi_res_time() - t_begin_main);
    
    t_start = hi_res_time();
    char ploidy_char = (char)(the_genotypes_set->ploidy + 48);
    for(long i=0; i<pedigrees->size; i++){
      if(i % 1000  == 0){
	fprintf(stdout, "# Done testing %ld pedigrees.\n", i);
      }
      
      Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(pedigrees->a[i], the_genotypes_set); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
      // assert(strcmp(pedigrees->a[i]->Accession->id, pedigrees->a[i]->A->id->a) == 0);
      Accession* A = pedigrees->a[i]->A;
      Accession* F = pedigrees->a[i]->F;
      Accession* M = pedigrees->a[i]->M;
      if(1 || F != NULL  &&  M != NULL){
   
	fprintf(o_stream, "%20s %5ld %20s %20s  ",
		A->id->a, A->missing_data_count,
		(F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA");
	/* Pedigree_stats* ps = the_pedigree_stats; */
	/* fprintf(o_stream, " %ld %ld  %ld %ld  %ld %ld  %ld %ld  %ld %ld  %ld %ld     ", */
	/* 	  ps->agmr12.n, ps->agmr12.d, ps->par1_hgmr.n, ps->par1_hgmr.d, ps->par1_R.n, ps->par1_R.d, */
	/* 	  ps->par2_hgmr.n, ps->par2_hgmr.d, ps->par2_R.n, ps->par2_R.d, ps->d.n, ps->d.d); */
  
	print_pedigree_stats(o_stream, the_pedigree_stats);
	//double epsilon = 0.01;
	//   two_doubles llF = lls(the_genotypes_set, F, A, o_stream, epsilon);
	//   two_doubles llM = lls(the_genotypes_set, M, A, o_stream, epsilon);
	//  four_longs F_fd = forbidden(the_genotypes_set, F, A);
	//  fprintf(o_stream, "  %ld %ld  %7.5f ", F_fd.l1, F_fd.l2, (F_fd.l2 > 0)? (double)F_fd.l1/F_fd.l2 : 2);

	//   two_longs ZZZ = diploid_quick_and_dirty_triple_counts(F, M, A);
	//  long ZZZd = ZZZ.l1 + ZZZ.l2;
	if(0){
	  four_longs zeroes4 = {0, 0, 0, 0};
	  four_longs FhR = (F != NULL)? hgmr_R(F->genotypes->a, A->genotypes->a, (char)(the_genotypes_set->ploidy+48)) : zeroes4;
	  four_longs MhR = (M != NULL)? hgmr_R(M->genotypes->a, A->genotypes->a, (char)(the_genotypes_set->ploidy+48)) : zeroes4;
	  fprintf(o_stream, "    %ld %8.5f  %ld %8.5f     %ld %8.5f %ld %8.5f  ",
		  FhR.l2, (FhR.l2>0)? (double)FhR.l1/FhR.l2 : 2.0,
		  FhR.l4, (FhR.l4>0)? (double)FhR.l3/FhR.l4 : 2.0,
		  MhR.l2, (MhR.l2>0)? (double)MhR.l1/MhR.l2 : 2.0,
		  MhR.l4, (MhR.l4>0)? (double)MhR.l3/MhR.l4 : 2.0
		  );
	}
	if(0){
	  four_longs tfc = //{0, 0, 0, 0};
	    tfca(F->genotypes->a, M->genotypes->a, A->genotypes->a, ploidy);
	  long numer1 = tfc.l1;
	  long denom1 = tfc.l2;
	  long numer2 = tfc.l3;
	  long denom2 = tfc.l4;
	    fprintf(o_stream, "   %ld %ld %8.5f  %ld %ld %8.5f  %ld %ld %8.5f", //  %ld %ld %8.5f",
		  numer1+numer2, denom1, (denom1 > 0)? (double)(numer1+numer2)/denom1 : 2,
		  numer1, denom1, (denom1>0)? (double)numer1/denom1 : 2,
		  numer2, denom2, (denom2>0)? (double)numer2/denom2 : 2);
	    four_longs qd_n22x = q_and_d_n22x_diploid(F,M,A);
	    fprintf(o_stream, "  %ld %ld %ld %ld  %8.5f  ",
		    qd_n22x.l1, qd_n22x.l2, qd_n22x.l3, qd_n22x.l4, (qd_n22x.l4 > 0)? (double)qd_n22x.l2/qd_n22x.l4 : 2.0);
	}
	if(0){
	  four_longs tfc_FF = tfca(F->genotypes->a, F->genotypes->a, A->genotypes->a, ploidy);
	  long numer1FF = tfc_FF.l1;
	  long denom1FF = tfc_FF.l2;
	  long numer2FF = tfc_FF.l3;
	  long denom2FF = tfc_FF.l4;	
	  fprintf(o_stream, "   %ld %ld %8.5f  %ld %ld %8.5f  %ld %ld %8.5f",
		  numer1FF+numer2FF, denom1FF, (denom1FF > 0)? (double)(numer1FF+numer2FF)/denom1FF : 2,
		  numer1FF, denom1FF, (denom1FF>0)? (double)numer1FF/denom1FF : 2,
		  numer2FF, denom2FF, (denom2FF>0)? (double)numer2FF/denom2FF : 2);
	}
	if(0){
	  Accession* a = pedigrees->a[i]->A;
	  Accession* f = pedigrees->a[i]->F;
	  Accession* m = pedigrees->a[i]->M;
	  ND hf = quick_and_dirty_hgmr(a, f, ploidy_char);
	  fprintf(o_stream, "   %ld %ld %6.5f  ", hf.n, hf.d, (hf.d>0)? (double)hf.n/hf.d : 2);
	}
	if(do_alternative_pedigrees == 1){
	  Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	  print_pedigree_alternatives(o_stream, alt_pedigrees);
	}else if((do_alternative_pedigrees == 2)){
	  // || (pedigree_ok(the_pedigree_stats, max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d) == 0)){
	  //   fprintf(stderr, "About to get alt. pedigrees. max_ok_hgmr: %8.4lf,  max_ok_d: %8.4lf \n", max_ok_hgmr, max_ok_d);
	  Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_z, max_ok_d);
	  print_pedigree_alternatives(o_stream, alt_pedigrees);
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

double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}


