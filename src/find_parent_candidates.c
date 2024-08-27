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
void print_command_line(FILE* ostream, int argc, char** argv);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************



int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();
  
  int do_alternative_pedigrees = 0; // 0: none, 1: only when given pedigree is 'bad', 2: all
  double max_marker_missing_data_fraction = 0.2; // default; control this with -x command line option.
  double min_minor_allele_frequency = -1; //
  char* pedigree_test_output_filename = "xpedigree_test_info";
  char* genotypes_matrix_output_filename = "genotype_matrix_out";
  double max_self_agmr12 = 1; // need to specify if doing alternative pedigrees
  double max_ok_hgmr = 0.02; // 
  double max_self_r = 1; // need to specify if doing alternative pedigrees
  double max_ok_d = 0.03; // accept everything as ok
  long max_parent_candidates = 30; // give up if more parent candidates than this.
  long max_queries = 1000000000;
  double ploidy = 2;
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
  // a: max 'self' agmr,
  // h: max ok hgmr,
  // r: max 'self' r,
  // D: max ok d;
  // n: max queries (default: do all)
  int c;
  //  int genotype_file_type = UNKNOWN;
  while((c = getopt(argc, argv, "A:cd:g:p:w:x:o:a:h:r:D:m:n:")) != -1){
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
    case 'n':
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_queries = atoi(optarg);
	if (max_queries < 0) exit(EXIT_FAILURE);
      }
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
      fprintf(stderr, "# min_minor_allele_freq: %7.5f \n", min_minor_allele_frequency);
      if(min_minor_allele_frequency < 0){
	fprintf(stderr, "option m (min_minor_allele_frequency) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
      /* case 'a': // agmr12 < this -> looks like self */
      /* 	if(optarg == 0){ */
      /* 	  perror("option a requires a numerical argument > 0\n"); */
      /* 	  exit(EXIT_FAILURE); */
      /* 	}else{ */
      /* 	  max_self_agmr12 = atof(optarg); */
      /* 	  if (max_self_agmr12 < 0) exit(EXIT_FAILURE); */
      /* 	} */
    case 'h': // hgmr>this -> poor parent-offspring candidate
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
  /* if(pedigrees_filename == NULL){ */
  /*   perror("must specify pedigrees filename: -i <filename>"); */
  /*   exit(EXIT_FAILURE); */
  /* } */

  
  // exit(EXIT_FAILURE);
  
  FILE *o_stream = NULL;
  o_stream = fopen(pedigree_test_output_filename, "w");
  if(o_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", pedigree_test_output_filename);
    exit(EXIT_FAILURE);
  }
   
  print_command_line(stdout, argc, argv);
  fprintf(stdout, "# genotypes filename: %s max marker missing data: %5.3lf\n",
	  genotypes_filename, max_marker_missing_data_fraction);
  fprintf(stdout, "# max_ok_hgmr: %8.5f  max_ok_d: %8.5f\n", max_ok_hgmr, max_ok_d);
  // fprintf(stderr, "# pedigrees filename: %s, output filename: %s \n", pedigrees_filename, pedigree_test_output_filename);

  
  // *****  done processing command line  *****

  // ***************  read the genotypes file  *******************************
  double t_start = hi_res_time();
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
  add_accessions_to_genotypesset_from_file(genotypes_filename, the_genotypes_set); // load the new set of accessions
  clean_genotypesset(the_genotypes_set);
  double ttt = hi_res_time();
  rectify_markers(the_genotypes_set);
  store_homozygs(the_genotypes_set);
  fprintf(stdout, "# time to rectify and store homozygs: %6.4f\n", hi_res_time() - ttt);

  double pppairs_time = 0;
  double fmptriples_time = 0;
  
  long nn = the_genotypes_set->accessions->size;
  if(nn > max_queries) nn =  max_queries;   
  for(long i=0; i<nn; i++){
    ttt = hi_res_time();
    Accession* the_accession = the_genotypes_set->accessions->a[i];
    //   Vlong* alt_parent_idxs = alternative_parents(the_accession, the_genotypes_set, max_ok_hgmr);
    /* fprintf(stderr, "%s  %ld  ", the_accession->id->a, alt_parent_idxs->size); */
    /* for(long j=0; j<alt_parent_idxs->size; j++){ */
    /*   fprintf(stderr, "%s  ", the_accession->id->a); */
    /*   fprintf(stderr, "%s \n", the_genotypes_set->accessions->a[alt_parent_idxs->a[j]]->id->a); */
    /* } */
    //fprintf(stderr, "\n");
    Vaccession* parent_candidates = construct_vaccession(max_parent_candidates+1);
    for(long j=0; j<the_genotypes_set->accessions->size; j++){
      if(j == i) continue; // accession cannot be it's own parent.
      Accession* the_other_accession = the_genotypes_set->accessions->a[j];
      //	four_longs x = forbidden(the_genotypes_set, the_accession, the_other_accession);
      double hgmr; //forbidden_rate, forbidden_rate_xx, forbidden_rate_xxx;
      ND x;
      if(ploidy == 2){
	x = quick_and_dirty_hgmr(the_accession, the_other_accession, (char)(the_genotypes_set->ploidy+48)); //forbidden_rate = (x.d > 0)? (double)x.n/x.d : 2; //
      }else{
        x = ghgmr(the_genotypes_set, the_accession, the_other_accession);
      }

	  hgmr = (x.d > 0)? (double)x.n/x.d : 2; //
	  fprintf(stderr, "%ld %ld %8.6f \n", x.n, x.d, hgmr);
      // 	ND    xxx =  hgmr_nd(the_accession->genotypes->a, the_other_accession->genotypes->a, (char)(the_genotypes_set->ploidy+48)); forbidden_rate_xxx = (xxx.d > 0)? (double)xxx.n/xxx.d : 2; //
	  //	ND xxxx = ghgmr_old(the_genotypes_set, the_accession, the_other_accession); long d = xxxx.d; long n = xxxx.n;	double forbidden_rate_xxxx = (d>0)? (double)n/d : 2;

	  //	  fprintf(stderr, "%ld %ld  %8.5f \n", xx.n, xx.d, forbidden_rate_xx); //, forbidden_rate_xxxx);
	//  	fprintf(stderr, "%8.5f  %8.5f  %8.5f  %8.5f\n", forbidden_rate, forbidden_rate_xx, forbidden_rate_xxx, forbidden_rate_xxxx);
	//double forbidden_rate_x = (xx.d > 0)? (double)xx.n/xx.d : 2;
	//	fprintf(stderr, "%ld %ld %8.5f \n", x.n, x.d, forbidden_rate); //, hgmrnd.n, hgmrnd.d, ((hgmrnd.d > 0)? (double)hgmrnd.n/hgmrnd.d : 2), forbidden_rate_x);
	if(hgmr < max_ok_hgmr){ // good parent candidate
	  add_accession_to_vaccession(parent_candidates, the_other_accession);
	}
       
      if(parent_candidates->size > max_parent_candidates) break;
    }
    pppairs_time += hi_res_time()-ttt;
    fprintf(stderr, "## max_ok_hgmr: %8.5f  n parent candidates: %ld \n", max_ok_hgmr, parent_candidates->size);
    long n_parent_candidates = parent_candidates->size;      
     
    // 	fprintf(stdout, "1 1 1 1 number of parent candidates: %ld\n", parent_candidates->size);
    long good_triple_count = 0;
    if(1 && parent_candidates->size < max_parent_candidates){ // do forbidden triples
      ttt = hi_res_time();
      char* prog_gts = the_accession->genotypes->a;
      //	fprintf(stderr, "accession %s  ; testing %ld parentages \n", the_accession->id->a, n_parent_candidates*(n_parent_candidates+1)/2);
      for(long ii = 0; ii < parent_candidates->size; ii++){
	Accession* parent1 = parent_candidates->a[ii];
	//	fprintf(stderr, "# ii: %ld  parent1 idx: %ld \n", ii, parent1->index);
	for(long jj = ii; jj < parent_candidates->size; jj++){
	  Accession* parent2 = parent_candidates->a[jj];
	  //	  fprintf(stderr, "# jj: %ld  parent2 idx: %ld \n", jj, parent2->index);
	  four_longs TFCs = tfca(parent1->genotypes->a, parent2->genotypes->a, prog_gts, ploidy);
	  ND tnd;
	  if(ploidy == 2){
	    tnd = tfc_diploid(parent1->genotypes->a, parent2->genotypes->a, prog_gts);
	  }else if(ploidy == 4){
	    tnd = tfc_tetraploid(parent1->genotypes->a, parent2->genotypes->a, prog_gts);
	  }
	  four_longs ftcs = triple_forbidden_counts(parent1->genotypes->a, parent2->genotypes->a, prog_gts, ploidy);
	  //	  two_longs dtcs = diploid_quick_and_dirty_triple_counts(parent1, parent2, the_accession);
	  double D1 = (TFCs.l2 > 0)? (double)TFCs.l1/TFCs.l2 : 2;
	   double D2 = (TFCs.l4 > 0)? (double)TFCs.l3/TFCs.l4 : 2;
	   double D3 = (TFCs.l2 > 0)? (double)(TFCs.l1 + TFCs.l3)/TFCs.l2 : 2;
	   double D4 = (tnd.d > 0)? (double)tnd.n/tnd.d  : 2;
	  double ddd = (ftcs.l4>0)? (double)ftcs.l1/ftcs.l4 : 2;
	  if(ddd < max_ok_d  ||  D1 < max_ok_d  ||  D2 < max_ok_d){
	    fprintf(stdout, "%20s     %20s %20s  %ld %ld %ld %ld  %8.6f  %8.6f  %8.6f  %8.6f %ld %ld %ld %ld  %8.6f %ld %ld\n", 
  // %ld %ld\n",
		    the_accession->id->a, parent1->id->a, parent2->id->a,
		    ftcs.l1, ftcs.l2, ftcs.l3, ftcs.l4, ddd, D1, D2, D3,
		    TFCs.l1, TFCs.l2, TFCs.l3, TFCs.l4, D4, tnd.n, tnd.d); //, dtcs.l1, dtcs.l2);
	    good_triple_count++;
	  }
	}
      }
      fmptriples_time += hi_res_time()-ttt;
      fprintf(stderr, "accession: %s  number of candidate parents: %ld  triples checked: %ld  good triples: %ld \n",
	      the_accession->id->a, parent_candidates->size, n_parent_candidates*(n_parent_candidates+1)/2, good_triple_count);
    } // end of getting triple forbidden counts for this progeny accession
     
      // free_vaccession(parent_candidates);
  }
  fprintf(stdout, "# time to get candidate parents for %ld accs: %6.5f\n", nn, pppairs_time);

  fprintf(stdout, "# time to get candidate parent pairs for %ld accs: %6.5f\n", nn, fmptriples_time);

  exit(0);

  // // // // // // //
    
  if(0){
    double t0 = hi_res_time();
    quick_and_dirty_hgmrs(the_genotypes_set);  
    fprintf(stderr, "# time for hgmrs: %10.3f \n", hi_res_time() - t0);
    exit(0);
  }
  t_start = hi_res_time();
  fprintf(stderr, "# n accessions: %ld \n", the_genotypes_set->n_accessions);
  Vidxid* the_vidxid = construct_sorted_vidxid(the_genotypes_set);
  fprintf(stderr, "# size of the_vidxid: %ld \n", the_vidxid->size);
  for(long i=0; i<the_vidxid->size; i++){
    long index = the_vidxid->a[i]->index;
    //   fprintf(stderr, "# i idx id: %ld %ld  %s\n", i, index, the_vidxid->a[i]->id); 
    Accession* a = the_genotypes_set->accessions->a[index];  
    //  fprintf(stderr, "# i, index, id, index: %ld %ld  %s  %ld\n", i, index, a->id->a, index_of_id_in_vidxid(the_vidxid, a->id->a));
  }
  // exit(0);
  // ***************  read the pedigrees file  ***************************
 
  const Vpedigree* pedigrees = read_the_pedigrees_file_and_store(p_stream, the_vidxid, the_genotypes_set);
  fclose(p_stream);
  fprintf(stderr, "# Done reading genotypes file and pedigree file. \n");
  fprintf(stderr, "# Stored genotypes of %ld accessions, and  %ld pedigrees. Time: %6.3f\n",
	  the_genotypes_set->accessions->size, pedigrees->size, hi_res_time() - t_start);

  // ***************  Done reading input files  ******************************

  
  const Vlong* parent_idxs = accessions_with_offspring(pedigrees, the_genotypes_set->n_accessions);
  printf("# According to pedigree file there are %ld accessions with offspring.\n", parent_idxs->size);

  t_start = hi_res_time();
  char ploidy_char = (char)(the_genotypes_set->ploidy + 48);
  for(long i=0; i<pedigrees->size; i++){
    if(i % 1000  == 0){
      fprintf(stderr, "# Done testing %ld pedigrees.\n", i);
    }
    Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(pedigrees->a[i], the_genotypes_set->ploidy); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
    // assert(strcmp(pedigrees->a[i]->Accession->id, pedigrees->a[i]->A->id->a) == 0);
    Accession* F = pedigrees->a[i]->F;
    Accession* M = pedigrees->a[i]->M;
    fprintf(o_stream, "%20s %5ld %20s %20s  ",
	    pedigrees->a[i]->A->id->a, pedigrees->a[i]->A->missing_data_count,
	    (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA");
    /* Pedigree_stats* ps = the_pedigree_stats; */
    /* fprintf(o_stream, " %ld %ld  %ld %ld  %ld %ld  %ld %ld  %ld %ld  %ld %ld     ", */
    /* 	  ps->agmr12.n, ps->agmr12.d, ps->par1_hgmr.n, ps->par1_hgmr.d, ps->par1_R.n, ps->par1_R.d, */
    /* 	  ps->par2_hgmr.n, ps->par2_hgmr.d, ps->par2_R.n, ps->par2_R.d, ps->d.n, ps->d.d); */
  
    print_pedigree_stats(o_stream, the_pedigree_stats, true);
 
    if(0){
      Accession* a = pedigrees->a[i]->A;
      Accession* f = pedigrees->a[i]->F;
      Accession* m = pedigrees->a[i]->M;
      //   ND af = quick_and_dirty_hgmr(a, f);
      //   ND am = quick_and_dirty_hgmr(a, m);
      ND hf = quick_and_dirty_hgmr(a, f, ploidy_char);
      //   ND hm = quick_and_dirty_hgmr_a(a, m, the_genotypes_set->ploidy);

      //  Pedigree_stats* ps = the_pedigree_stats; // triple_counts(f->genotypes->a, m->genotypes->a, a->genotypes->a);
      //  long n = ps->par1_hgmr.n;
      // long d = ps->par1_hgmr.d;
      //   fprintf(o_stream, "  %ld %ld  %7.5f  ",  n, d, (d > 0)? (double)n/d : 2);
      // fprintf(o_stream, "   %ld %ld %ld %ld   ", af.n, af.d, am.n, am.d);
      //    ND hhh = hgmr_nd(a->genotypes->a, f->genotypes->a);
      //   fprintf(o_stream, "  %ld %ld  %6.5f  ", hhh.n, hhh.d, (hhh.d>0)? (double)hhh.n/hhh.d : 2);
      fprintf(o_stream, "   %ld %ld %6.5f  ", hf.n, hf.d, (hf.d>0)? (double)hf.n/hf.d : 2);
    }
  
    if(do_alternative_pedigrees == 1){
      Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_d);
      print_pedigree_alternatives(o_stream, alt_pedigrees);
    }else if((do_alternative_pedigrees == 2) || (pedigree_ok(the_pedigree_stats, max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d) == 0)){
      //   fprintf(stderr, "About to get alt. pedigrees. max_ok_hgmr: %8.4lf,  max_ok_d: %8.4lf \n", max_ok_hgmr, max_ok_d);
      Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_d);
      print_pedigree_alternatives(o_stream, alt_pedigrees);
    }
    fprintf(o_stream, "\n");
    free(the_pedigree_stats);
  }
  
  fprintf(stderr, "# Done testing %ld pedigrees.\n", pedigrees->size);

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
  fprintf(stderr, "# Total time: %10.4lf sec.\n", hi_res_time() - t_begin_main);
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

void print_command_line(FILE* ostream, int argc, char** argv){
  fprintf(ostream, "# command line:  ");
  for(int i=0; i<argc; i++){
    fprintf(ostream, "%s  ", argv[i]);
  }fprintf(ostream, "\n");
}
