// C version of 'k-mer' search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <assert.h>

#include "gtset.h"
//#include "various.h"

//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // one of these for each chunk
  long capacity;
  long size;
  Vlong** a; // array of Vlong*; indices are patterns, values are Vlong*s containing indices of accessions with that pattern
} Pattern_ids;

typedef struct{
  long capacity;
  long size; // number of chunks
  long chunk_size; // number of chunks used.
  long n_patterns; // (ploidy+1)^chunk_size
  Pattern_ids** a; // array of Pattern_ids, one element for each chunk
} Chunk_pattern_ids;

typedef struct{
  long query_index;
  long match_index;
  double usable_chunks; // estimate or actual count
  long n_matching_chunks;
  double est_agmr;
  double agmr;
  Vdouble* agmrs;
  //  double d1;
  //  double hgmr; 
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci;

int  do_checks = 0;


// *********************** function declarations ************************************************
void print_usage_info(FILE* ostream);
char* ipat_to_strpat(long len, long ipat); // unused
long strpat_to_ipat(long len, char* strpat); // unused
double agmr(Accession* gts1, Accession* gts2);
Vdouble* maf_range_agmrs(GenotypesSet* the_gtset, Accession* acc1, Accession* acc2, Vdouble* maf_threshholds);

// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, double n_usable_chunks, long n_matching_chunks,
		   double est_agmr, double agmr, Vdouble* agmrs); //, double d1,double hgmr);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void add_mci_to_vmci(Vmci* the_vmci, Mci* the_mci);
void sort_vmci_by_agmr(Vmci* the_vmci); // sort Vmci by agmr
void sort_vmci_by_index(Vmci* the_vmci); // sort Vmci by index
int cmpmci_a(const void* v1, const void* v2);
int cmpmci_i(const void* v1, const void* v2);
void free_vmci(Vmci* the_vmci);

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long n_patterns);
void free_pattern_ids(Pattern_ids*);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size, long ploidy);
// Vlong** get_all_match_counts(long n_accessions, Chunk_pattern_ids* the_cpi);
void populate_chunk_pattern_ids_from_vaccession(Vaccession* the_accessions, Chunk_pattern_ids* the_cpi);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream);
void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Accession* the_gts, Chunk_pattern_ids* the_cpi);
Vmci** find_matches(GenotypesSet* the_genotypes_set,
		    //long n_ref_accessions, Vaccession* the_accessions,
		    Chunk_pattern_ids* the_cpi, double max_est_agmr);
long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long out_format);
long print_results_a(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long out_format);
void print_command_line(FILE* ostream, int argc, char** argv);
// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  double start0 = hi_res_time();
  
  long ploidy = 2; //
  long chunk_size = -1; // default: choose automatically based on ploidy, etc.
  long n_passes = 1; // n_chunks = n_passes*(int)(n_markers/chunk_size)
  // long n_chunks = 1000000; // default number of chunks (large number -> use all markers)
  unsigned rand_seed = (unsigned)time(0);
  double max_marker_missing_data_fraction = -1; // default: set to 2.0/chunk_size after chunk_size is set.
  double max_accession_missing_data_fraction = 0.5;
  double min_minor_allele_frequency = 0.0; //
  double max_est_agmr = 0.2;
  long output_format = 1; // 1 ->  acc_id1 acc_id2  n_usable_chunks n_matching_chunks est_agmr agmr
  char default_output_filename[] = "duplicatesearch.out";

  char* rparam_buf;
  size_t rparam_len;
  FILE* rparam_stream = open_memstream(&rparam_buf, &rparam_len);

  // ***** process command line *****
  if (argc < 2) {
    print_usage_info(stdout);
    exit(EXIT_FAILURE);
  }

  char* input_filename = NULL;
  FILE *in_stream = NULL;
  char* reference_set_filename = NULL;
  FILE *ref_in_stream = NULL;
  char* output_filename = default_output_filename;
  FILE* out_stream = NULL;
    
  int c;
  while((c = getopt(argc, argv, "i:r:o:n:k:e:s:x:a:f:h")) != -1){
    // i: input file name (required).
    // r: reference set file name (optional).
    // o: output file name. Default:
    // s: random number seed. Default: get seed from clock.
    // x: marker max missing data fraction. Default: 2.0/chunk_size
    // f: accession max missing data fraction. Default: 0.5
    // a: min minor allele frequency. Default: 0
    // k: chunk size (number of markers per chunk). Defaults: diploid: 8, tetraploid: 5, hexaploid: 4.
    // e: max agmr. Default: 0.2 (Calculate agmr only if quick est. is < this value; output match only if agmr < this value.)
    // n: use n*n_markers/chunk_size chunks. (i.e. each marker gets used in ~n chunks). Default: 1
    // v: verbosity. Default: 1, 2 gives slightly more output.
    // h: help. print usage info
     
    switch(c){
    case 'i':
      input_filename = optarg;
      in_stream = fopen(input_filename, "r");
      if(in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
      fclose(in_stream);
      break;
    case 'r':
      reference_set_filename = optarg;
      ref_in_stream = fopen(reference_set_filename, "r");
      if(ref_in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", reference_set_filename);
	exit(EXIT_FAILURE);
      }
      fclose(ref_in_stream);
      break;
    case 'o':
      output_filename = optarg;
      break;
    /* case 'p': // keep each accession with probability p (for testing with random smaller data set) */
    /*   ploidy = (long)atoi(optarg); */
    /*   if(ploidy <= 0){ */
    /* 	fprintf(stderr, "# ploidy specified as %ld. Must be non-negative integer; exiting.\n", ploidy); */
    /* 	exit(EXIT_FAILURE); */
    /*   } */
    /*   break; */
    case 'n': 
      n_passes = (long)atoi(optarg);
      if(n_passes <= 0){
	fprintf(stderr, "option n (n_passes) requires an integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'k': // chunk size
      chunk_size = (long)atoi(optarg);
      if(chunk_size <= 0){
	fprintf(stderr, "option k (chunk_size) requires a integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break; 
    case 'e': 
      max_est_agmr = (double)atof(optarg);
      if(max_est_agmr <= 0  || max_est_agmr > 1.0){
	fprintf(stderr, "option e (max_est_agmr) requires a numerical argument 0<x<=1 \n");
	exit(EXIT_FAILURE);
      }
      break;
    case 's': // random number generator seed
      rand_seed = (unsigned)atoi(optarg);
      if(rand_seed <= 0){
	fprintf(stderr, "option s (rng seed) requires a integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'x':
      max_marker_missing_data_fraction = (double)atof(optarg);
      if(max_marker_missing_data_fraction <= 0){
	fprintf(stderr, "option x (max_marker_missing_data_fraction) requires an real argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'f':
      max_accession_missing_data_fraction = (double)atof(optarg);
      if(max_accession_missing_data_fraction <= 0){
	fprintf(stderr, "option f (max_accessions_missing_data_fraction) requires an real argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'a': 
      min_minor_allele_frequency = (double)atof(optarg);
      if(min_minor_allele_frequency < 0){
	fprintf(stderr, "option a (min_minor_allele_frequency) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'v': 
      output_format = (long)atoi(optarg);
      if(output_format < 1){
	fprintf(stderr, "option v (output_format) requires an integer argument >= 1\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'h':
      print_usage_info(stderr);
      exit(EXIT_FAILURE);
      break;
    case '?':
      fprintf(stderr, "? case in command line processing switch.\n");
      if ((optopt == 'i') || (optopt == 'r') || (optopt == 'o') ||
	  (optopt == 'p') || (optopt == 'n') || (optopt == 'k') ||
	  (optopt == 'e') || (optopt == 's') || (optopt == 'x')  || (optopt == 'a') )
	fprintf (stderr, "  Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf (stderr, "  Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr, "  Unknown option character: %d\n", optopt);
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr, "default case (abort)\n");
      exit(EXIT_FAILURE);
    } // end of switch block
  } // end of loop over c.l. arguments

  if(optind < argc){
    fprintf(stderr, "Non-option arguments. Bye.\n");
    exit(EXIT_FAILURE);
  }  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  srand(rand_seed);
  out_stream = fopen(output_filename, "w");
  if(out_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
    exit(EXIT_FAILURE);
  }

  print_command_line(out_stream, argc, argv);
  print_command_line(stdout, argc, argv);	 
  fprintf(rparam_stream, "# Rng seed: %ld\n", (long)rand_seed);
  fprintf(rparam_stream, "# Min. marker minor allele frequency: %5.3lf\n", min_minor_allele_frequency);
  if(max_marker_missing_data_fraction > 0) {
    fprintf(rparam_stream, "# Max. marker missing data fraction: %5.3lf\n", max_marker_missing_data_fraction);
  }else{
    fprintf(rparam_stream, "# Max. marker missing data fraction will be set to 2.0/chunk_size.\n");
  }
  fprintf(rparam_stream, "# Max. agmr: %5.3lf\n", max_est_agmr);
    fclose(rparam_stream);
  fprintf(stdout, "%s", rparam_buf);
  fprintf(out_stream, "%s", rparam_buf);
  free(rparam_buf);
  
  // *****  done processing command line  *****************************************

  long n_ref_accessions = 0;
  long n_accessions = 0;
  long n_markers = 0;

  // *****  read in genotype data (including optionally a reference data set) and create a genotypesset object.

  double t_start = hi_res_time();
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy); 
  Vaccession* the_accessions; 
  if(reference_set_filename != NULL){ // load the reference set, if one was specified.
   
    add_accessions_to_genotypesset_from_file(reference_set_filename, the_genotypes_set, max_accession_missing_data_fraction);
    n_ref_accessions = the_genotypes_set->accessions->size;
    the_genotypes_set->n_ref_accessions = n_ref_accessions;
    fprintf(stdout, "# Done reading reference data set dosages from file %s. %ld accessions and %ld markers.\n",
	    reference_set_filename, the_genotypes_set->n_accessions, the_genotypes_set->n_markers);
  }
  add_accessions_to_genotypesset_from_file(input_filename, the_genotypes_set, max_accession_missing_data_fraction); // load the new leset of accessions
  fprintf(stdout, "# Done reading dosages from file %s. %ld accessions and %ld markers.\n",
	  input_filename, the_genotypes_set->n_accessions, the_genotypes_set->n_markers);
  fprintf(stdout, "# Time to load dosage data: %6.3lf sec.\n", hi_res_time() - t_start);

    ploidy = the_genotypes_set->ploidy;
  //  find  chunk_size  if not specified on command line; 
    if(chunk_size <= 0){ // with MAX_PATTERNS of 10000, chunk_size is 8 for diploid, 5 for tetraploid, 4 for hexaploid.
    chunk_size = (long)(log(MAX_PATTERNS)/log((double)ploidy+1.0));
  }
  the_genotypes_set->max_marker_missing_data_fraction = (max_marker_missing_data_fraction <= 0)? 2.0/chunk_size : max_marker_missing_data_fraction;
  fprintf(out_stream, "# Max. marker missing data fraction: %5.3lf\n", the_genotypes_set->max_marker_missing_data_fraction);
  fprintf(stdout, "# Max. marker missing data fraction: %5.3lf\n", the_genotypes_set->max_marker_missing_data_fraction);

  // fprintf(stdout, "# pre-cleaning ragmr: %8.6f \n", ragmr(the_genotypes_set));

  rectify_markers(the_genotypes_set); // swap dosage 0 and 2 for markers with dosage more common, so afterward 0 more common that 2 for all markers.
  clean_genotypesset(the_genotypes_set);
  store_homozygs(the_genotypes_set);
  // fprintf(stdout, "# post-cleaning ragmr: %8.6f \n", ragmr(the_genotypes_set));

  the_accessions = the_genotypes_set->accessions;
   
  n_markers = the_genotypes_set->n_markers;

  //  if(n_chunks*chunk_size > n_markers){
  //    n_chunks = n_markers/chunk_size;
  //  }
  long n_chunks_per_pass = n_markers/chunk_size;
  long n_chunks = n_passes*n_chunks_per_pass;
  fprintf(out_stream, "# Cleaned data has %ld markers.\n", the_genotypes_set->n_markers);
  fprintf(stdout, "# Chunk size: %ld  n_chunks: %ld\n", chunk_size, n_chunks);
    fprintf(out_stream, "# chunk size: %ld  n_chunks: %ld\n", chunk_size, n_chunks);

  // *****  done reading and storing input  **********
  
 
  //  fprintf(stderr, "# n_markers: %ld\n", n_markers);
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices);
  marker_indices->size = n_chunks_per_pass*chunk_size; //
  fprintf(stderr, "i: 0  marker_indices->size: %ld \n", marker_indices->size);
  for(long i=1; i<n_passes; i++){
    Vlong* more_marker_indices = construct_vlong_whole_numbers(n_markers);
    shuffle_vlong(more_marker_indices);
    more_marker_indices->size = n_chunks_per_pass*chunk_size;
    append_vlong_to_vlong(marker_indices, more_marker_indices);
    // fprintf(stderr, "i: %ld  marker_indices->size: %ld \n", i, marker_indices->size);
  }
  //  exit(0);
   t_start = hi_res_time();
  set_vaccession_chunk_patterns(the_accessions, marker_indices, n_chunks, chunk_size, ploidy);
  double t_1 = hi_res_time();
  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, chunk_size, ploidy);
  populate_chunk_pattern_ids_from_vaccession(the_accessions, the_cpi);
  double t_2 = hi_res_time();
  //  fprintf(stdout, "# Time to construct chunk_pattern_ids. patterns: %8.4f  populate %8.4f  total %8.4f\n", t_1 - t_start, t_2 - t_1, t_2 - t_start);
    fprintf(stdout, "# Time to construct & populate chunk_pattern_ids: %6.3f\n", t_2 - t_start);

  
  t_start = hi_res_time(); 
  // Vmci** query_vmcis = find_matches(n_ref_accessions, the_accessions, the_cpi, max_est_agmr);
  Vmci** query_vmcis = find_matches(the_genotypes_set, the_cpi, max_est_agmr);
  long true_agmr_count = print_results(the_accessions, query_vmcis, out_stream, output_format);
  fprintf(stdout, "# Time to find candidate matches and %ld true agmrs: %6.3f\n", true_agmr_count, hi_res_time() - t_start);
  fclose(out_stream);

  long cume_s = 0;
  // *****  clean up  *****
  for(long i=0; i< the_accessions->size; i++){
    long s = query_vmcis[i]->size;
    cume_s += s;
    //  fprintf(stderr, "# i, size of query_vmcis[i]: %ld  %ld \n", i, s);
     free_vmci(query_vmcis[i]);
  }
  // fprintf(stderr, "# mcis freed: %ld\n", cume_s);
  free(query_vmcis);
  free_genotypesset(the_genotypes_set); 
  free_vlong(marker_indices);
  free_chunk_pattern_ids(the_cpi);
  fprintf(stdout, "# total duplicatesearch run time: %9.3f\n", hi_res_time() - start0);
  exit(EXIT_SUCCESS);
}

// **********************************************************************************************
// **********************  end of main  *********************************************************
// **********************************************************************************************


// *******************  function definitions  ***************************************************
// *****  Vaccession  ***********************************************************


void populate_chunk_pattern_ids_from_vaccession(Vaccession* the_accessions, Chunk_pattern_ids* the_cpi){
  long n_patterns = the_cpi->n_patterns;
 
  for(long i_gts=0; i_gts<the_accessions->size; i_gts++){
    Accession* the_gts = the_accessions->a[i_gts];
    if(DO_ASSERT) assert(i_gts == the_gts->index);
    Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts 
    long mdcount = 0;
    
    for(long i=0; i<the_chunk_patterns->size; i++){
      if(the_chunk_patterns->a[i] == n_patterns){ mdcount++; }
    }
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      if(DO_ASSERT) assert(the_pat >= 0);

      if(the_cpi->a[i_chunk]->a[the_pat] == NULL){
	the_cpi->a[i_chunk]->a[the_pat] = construct_vlong(1);
      }
      Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];	
      add_long_to_vlong(the_accidxs, the_gts->index);
    }
  }
  
  long total_mdchunk_count = 0;
  for(long i=0; i<the_cpi->size; i++){
    long chunk_md_count = (the_cpi->a[i]->a[n_patterns] == NULL)?
      0 :  // there are no accessions having missing data for this chunk
      the_cpi->a[i]->a[n_patterns]->size;
    total_mdchunk_count += chunk_md_count;
  } 
}

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accidxs having that pattern.

Pattern_ids* construct_pattern_ids(long n_patterns){
  Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids));
  pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=(ploidy+1)^k) good patterns, and index n_pattern is for the missing data case.
  pat_ids->size = n_patterns+1;
  pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(long ipat=0; ipat< pat_ids->size; ipat++){
    pat_ids->a[ipat] = NULL; // construct_vlong(1); // waste of memory? set to NULL until needed?
  }
  return pat_ids;
}

void free_pattern_ids(Pattern_ids* pat_ids){
  if(pat_ids == NULL) return;
  for(long i=0; i<pat_ids->size; i++){
    if(pat_ids->a[i] != NULL) free_vlong(pat_ids->a[i]);
  }
  free(pat_ids->a);
  free(pat_ids);
}

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*

Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size, long ploidy){ // needed size is known at construct time, so one param for both cap and size
  Chunk_pattern_ids* chunk_pat_ids = (Chunk_pattern_ids*)malloc(1*sizeof(Chunk_pattern_ids));
  chunk_pat_ids->capacity = n_chunks;
  chunk_pat_ids->size = n_chunks;
  chunk_pat_ids->chunk_size = chunk_size;
  long n_patterns = int_power(ploidy+1, chunk_size);
  chunk_pat_ids->n_patterns = n_patterns;
  chunk_pat_ids->a = (Pattern_ids**)malloc(n_chunks*sizeof(Pattern_ids*));
  for(int i=0; i< chunk_pat_ids->size; i++){
    chunk_pat_ids->a[i] = construct_pattern_ids(n_patterns);
  }
  return chunk_pat_ids;
}

void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi){
  if(the_cpi == NULL) return;
  for(long i=0; i< the_cpi->size; i++){
    free_pattern_ids(the_cpi->a[i]);
  }
  free(the_cpi->a);
  free(the_cpi);
}

void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream){
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    fprintf(ostream, "i_chunk: %ld\n", i_chunk);
    Pattern_ids* the_pi = the_cpi->a[i_chunk];
    for(long i_pat=0; i_pat<the_pi->size; i_pat++){
      fprintf(ostream, "  i_pat: %ld \n", i_pat);
      Vlong* the_idxs = the_pi->a[i_pat];
      if(the_idxs->size > 0){
	fprintf(ostream, "     matches: ");
	for(long ii=0; ii<the_idxs->size; ii++){
	  fprintf(ostream, "%ld  ", the_idxs->a[ii]);
	}
	fprintf(ostream, "\n");
      }
    }
  }
}

// *****  Accession and Chunk_pattern_ids  ***********

Vlong* find_chunk_match_counts(Accession* the_gts, Chunk_pattern_ids* the_cpi){ //, long n_accessions){ //, Vlong** accidx_hmatchcounts){
  // fprintf(stderr, "top of find_chunk_match_counts\n");
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_matchcounts = construct_vlong_zeroes(the_gts->index); //   n_accessions);
 
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];  
    Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat
    // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data )
    if(the_pat == n_patterns){ // missing data in this chunk 
    }else{ // the_pat = 0..n_patterns-1 (good data)   
      // just get the counts for matches with index < index of the_gts
      // since those are the only ones used in find_matches. (slightly faster)
      for(long i=0; i<chunk_match_idxs->size; i++){	  
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	if(accidx >= the_gts->index) break;
	accidx_matchcounts->a[accidx]++; // accidx here is < the_gts->index
      }
    }
  }
  // fprintf(stderr, "bottom of find_chunk_match_counts\n");
  return accidx_matchcounts; 
} 

// ***** Mci  *****

Mci* construct_mci(long qidx, long midx, double usable_chunks, long n_matching_chunks,
		   // double est_matching_chunk_fraction, double matching_chunk_fraction){
		   double est_agmr, double agmr, Vdouble* agmrs){ //, double d1, double hgmr){
  Mci* the_mci = (Mci*)calloc(1,sizeof(Mci));
  the_mci->query_index = qidx;
  the_mci->match_index = midx;
  the_mci->usable_chunks = usable_chunks;
  the_mci->n_matching_chunks = n_matching_chunks;
  the_mci->est_agmr = est_agmr;
  the_mci->agmr = agmr;
  the_mci->agmrs = agmrs;
  //  the_mci->d1 = d1;
  //  the_mci->hgmr = hgmr;
  
  return the_mci;
}

// *****  Vmci  *********************************************************************************

Vmci* construct_vmci(long init_cap){
  Vmci* the_vmci = (Vmci*)malloc(1*sizeof(Vmci));
  the_vmci->capacity = init_cap;
  the_vmci->size = 0;
  the_vmci->a = (Mci**)malloc(init_cap*sizeof(Mci*));
  return the_vmci;
}

void add_mci_to_vmci(Vmci* the_vmci, Mci* mci){
  long cap = the_vmci->capacity;
  long n = the_vmci->size;
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vmci->a = (Mci**)realloc(the_vmci->a, cap*sizeof(Mci*));
    the_vmci->capacity = cap;
  }
  the_vmci->a[n] = mci;
  the_vmci->size++;
}

void free_vmci(Vmci* the_vmci){
  if(the_vmci == NULL) return;
  for(long i=0; i< the_vmci->size; i++){
    free(the_vmci->a[i]);
  }
  free(the_vmci->a);
  free(the_vmci);
}


int cmpmci_i(const void* v1, const void* v2){  // sort by query index, agmr is tie breaker
  const Mci** s1 = (const Mci**)v1;
  const Mci** s2 = (const Mci**)v2;

  long idx1 = (*s1)->query_index;
  long idx2 = (*s2)->query_index;

  if(idx1 > idx2){
    return 1;
  }else if(idx1 < idx2){
    return -1;
  }else{
    double a1 = (*s1)->agmr;
    double a2 = (*s2)->agmr;
    if(a1 > a2){
      return 1;
    }else if(a1 < a2){
      return -1;
    }else{
      return 0;
    }
  }
}

int cmpmci_a(const void* v1, const void* v2){ // sort by agmr, query index is tie breaker
  const Mci** s1 = (const Mci**)v1;
  const Mci** s2 = (const Mci**)v2;
  double a1 = (*s1)->agmr;
  double a2 = (*s2)->agmr;

  if(a1 > a2){
    return 1;
  }else if(a1 < a2){
    return -1;
  }else{
    long idx1 = (*s1)->query_index;
    long idx2 = (*s2)->query_index;
    if(idx1 > idx2){
      return 1;
    }else if(idx1 <  idx2){
      return -1;
    }else{
      return 0;
    }
  }
}

void sort_vmci_by_agmr(Vmci* the_vmci){ // sort by agmr (low to high), query index as tie-breaker
  qsort(the_vmci->a, the_vmci->size, sizeof(Mci*), cmpmci_a);
}

void sort_vmci_by_index(Vmci* the_vmci){ // sort by agmr (low to high), query index as tie-breaker
  qsort(the_vmci->a, the_vmci->size, sizeof(Mci*), cmpmci_i);
}



// *********************************************
// *********************************************

double agmr(Accession* gtset1, Accession* gtset2){
  char* gts1 = gtset1->genotypes->a;
  char* gts2 = gtset2->genotypes->a;
  long usable_pair_count = 0; // = agmr_denom
  long mismatches = 0; // = agmr_numerator
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break; // end of 
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0');
    if(a1 != MISSING_DATA_CHAR){
      if(a2 != MISSING_DATA_CHAR){
	usable_pair_count++;
	if(a1 != a2) mismatches++;
      }
    }
  }
  return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1;
}

Vdouble* maf_range_agmrs(GenotypesSet* the_gtset, Accession* acc1, Accession* acc2, Vdouble* maf_threshholds){
  add_double_to_vdouble(maf_threshholds, 1.0);
  Vdouble* agmrs = construct_vdouble(maf_threshholds->size);
  long n_markers = acc1->genotypes->length;
  if(DO_ASSERT) assert(acc2->genotypes->length == n_markers);
  char* gts1 = acc1->genotypes->a;
  char* gts2 = acc2->genotypes->a;
  Vlong* numerators =  construct_vlong_zeroes(maf_threshholds->size);
  Vlong* denominators =  construct_vlong_zeroes(maf_threshholds->size);
  Vlong* marker_alt_allele_counts = the_gtset->marker_alt_allele_counts;
  Vlong* marker_md_counts = the_gtset->marker_missing_data_counts;
  // fprintf(stderr, "n genotypes: %ld  %ld \n", acc1->genotypes->length, acc2->genotypes->length);
  for(long i=0; i<n_markers; i++){ 
    char a1 = gts1[i];
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0'  &&  a1 != '\0');
    if(a1 != MISSING_DATA_CHAR){
      if(a2 != MISSING_DATA_CHAR){
	double marker_alt_allele_freq = marker_alt_allele_counts->a[i]/(double)(the_gtset->n_accessions - marker_md_counts->a[i]);
	for(long j=0; j< maf_threshholds->size; j++){
	  if(marker_alt_allele_freq < maf_threshholds->a[j]){
	    
	    denominators->a[j]++;
	    if(a1 != a2) numerators->a[j]++;
	    break;
	  }
	} // end loop over maf ranges
      }
    }
  } // end loop over markers
  agmrs->size = numerators->size;
  for(long j=0; j< agmrs->size; j++){
    agmrs->a[j] = (double)numerators->a[j]/(double)denominators->a[j];
  }
  return agmrs;
}

Vmci** find_matches(GenotypesSet* the_genotypes_set,
		    // long n_ref_accessions, Vaccession* the_accessions,
		    Chunk_pattern_ids* the_cpi,
		    //long min_usable_chunks,
		    double max_est_agmr)
{
  clock_t start = clock();
  clock_t fcmc_ticks = 0;

  Vaccession* the_accessions = the_genotypes_set->accessions;
  long n_ref_accessions = the_genotypes_set->n_ref_accessions;
  long n_markers = the_accessions->a[0]->genotypes->length;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;
  
  long true_agmr_count = 0;
  // long xcount = 0;
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);
  Vmci** query_vmcis = (Vmci**)malloc(the_accessions->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_accessions->size; i++){
    query_vmcis[i] = construct_vmci(4);
  }
  for(long i_query=n_ref_accessions; i_query< the_accessions->size; i_query++){
  
    Accession* q_gts = the_accessions->a[i_query];
    long q_md_chunk_count = q_gts->md_chunk_count;
    //clock_t ticks_before_fcmc = clock();
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi); //, n_ref_accessions);
    //fcmc_ticks += clock() - ticks_before_fcmc;
    
    for (long i_match = 0; i_match < i_query; i_match++){
      long matching_chunk_count = chunk_match_counts->a[i_match];
      long match_md_chunk_count = the_accessions->a[i_match]->md_chunk_count;
      // xxx
	
      double usable_chunk_count = (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks; // estimate
      //  fprintf(stderr, "# n md chunks, query: %ld  match: %ld  est number of usable chunk pairs: %8.3lf \n", q_md_chunk_count, match_md_chunk_count, usable_chunk_count);
      
      if( matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count ){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	// double true_hgmr;
	double true_agmr = agmr(q_gts, the_accessions->a[i_match]); //, &true_hgmr);
	// Three_ds dists = poly_agmr(q_gts, the_accessions->a[i_match]);
	//double true_agmr = dists.d1;
	double mafs[1] = {0.07};
	Vdouble* maf_threshholds = construct_vdouble_from_array(1, mafs);
	//	fprintf(stderr, "number of maf categories: %ld \n", maf_threshholds->size);
	if(true_agmr <= max_est_agmr){
	  Vdouble* agmrs = maf_range_agmrs(the_genotypes_set, q_gts, the_accessions->a[i_match], maf_threshholds);
	  //	  fprintf(stderr, "number of agmrs: %ld \n", agmrs->size);
	  true_agmr_count++;
	  add_mci_to_vmci(query_vmcis[i_query],
			  construct_mci(i_query, i_match, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, agmrs)); //, dists.d2, dists.d3)); //true_hgmr))	  //  fprintf(stderr, "# i_query i_match: %ld %ld \n", i_query, i_match);
	  if(i_match >= n_ref_accessions){ add_mci_to_vmci(query_vmcis[i_match],
							   construct_mci(i_match, i_query, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, agmrs)); //, dists.d2, dists.d3)); // true_hgmr));
	    //	    xcount++;
	  }
	} // end if(true_agmr < max_est_agmr)
      } // end if(enough matching chunks)
    } // end loop over potential matches to query
    free_vlong(chunk_match_counts);
  } // end loop over queries.
  //  fprintf(stderr, "# n_ref_accessions: %ld true_agmr_count: %ld  xcount: %ld\n", n_ref_accessions, true_agmr_count, xcount);
  //clock_t find_matches_ticks = clock() - start; 
  /* fprintf(stderr, "# time in: find_chunk_match_count: %8.3lf; rest of find_matches: %8.3lf; find_matches total: %8.3lf\n", */
  /* 	  clock_ticks_to_seconds(fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks - fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks)); */
  return query_vmcis;
}

long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long output_format){
  long true_agmr_count = 0;
  for(long i_q=0; i_q<the_accessions->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    sort_vmci_by_agmr(the_vmci);
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      Mci* the_mci = the_vmci->a[i_m];
      //  long match_idx = the_mci->match_index; //(the_mci->query_index == i_q)? the_mci->match_index : the_mci->query_index;
      Accession* q_acc = the_accessions->a[i_q];
      Accession* m_acc = the_accessions->a[the_mci->match_index];

      
      fprintf(ostream, "%26s  %26s  %5.2f  %3ld %7.4f  %8.6f",
	      // i_q,
	      q_acc->id->a,   m_acc->id->a,  
	      the_mci->usable_chunks,  the_mci->n_matching_chunks,
	      the_mci->est_agmr,  the_mci->agmr);
      Vdouble* the_agmrs = the_mci->agmrs;
      for(long iii=0; iii<the_agmrs->size; iii++){
	fprintf(ostream, "  %7.4f", the_agmrs->a[iii]);
      }
      if (output_format == 1){
	// leave as is
      }else if(output_format == 2){ // add a bit more info
	fprintf(ostream, "  %4ld  %4ld  %3ld  %3ld",
		q_acc->missing_data_count,  q_acc->md_chunk_count, m_acc->missing_data_count, m_acc->md_chunk_count);
      }else{
	fprintf(stderr, "# output_format %ld is unknown. using default output_format.\n", output_format);
      }
      fprintf(ostream, "\n");
      true_agmr_count++;
    } // end of loop over matches to query
  } // end of loop over queries
  return true_agmr_count;
}

long print_results_a(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long output_format){
  Vmci* all_mcis = construct_vmci(1000);
  for(long i_q=0; i_q<the_accessions->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      // Mci* the_mci = the_vmci->a[i_m];
      add_mci_to_vmci(all_mcis, the_vmci->a[i_m]);
    }
  }
  long true_agmr_count = 0;
  // sort_vmci_by_agmr(all_mcis);
  sort_vmci_by_index(all_mcis);
  for(long i=0; i< all_mcis->size; i++){
    Mci* the_mci = all_mcis->a[i];
    Accession* q_acc = the_accessions->a[the_mci->query_index];
    Accession* m_acc = the_accessions->a[the_mci->match_index];
      fprintf(ostream, "%4ld  %30s  %30s  %5.2f  %3ld  %5.3f  %5.3f",
	      the_mci->query_index,  q_acc->id->a,   m_acc->id->a,  
	      the_mci->usable_chunks,  the_mci->n_matching_chunks,
	      the_mci->est_agmr,  the_mci->agmr);
      if (output_format == 1){
	// leave as is
      }else if(output_format == 2){ // add a bit more info
	fprintf(ostream, "  %4ld  %4ld  %3ld  %3ld",
		q_acc->missing_data_count,  q_acc->md_chunk_count, m_acc->missing_data_count, m_acc->md_chunk_count);
      }else{
	fprintf(stderr, "# output_format %ld is unknown. using default output_format.\n", output_format);
      }
      fprintf(ostream, "\n");
      true_agmr_count++;
    } // end of loop over matches to query
  free_vmci(all_mcis);
  return true_agmr_count;
} // end print_results_a


void print_usage_info(FILE* ostream){
  // i: input file name (required).
  // r: reference set file name.
  // o: output file name. Default: duplicatesearch.out
  // e: max estimated agmr. Default: 0.2 (Calculate agmr only if quick est. is < this value.) 
  // x: marker max missing data fraction
  // a: min minor allele frequency
  // k: chunk size (number of markers per chunk). Default: 8
  // s: random number seed. Default: get seed from clock.
  // f: output format control. Default: 1; 2 gives additional info.
  // n: number of chunks to use. Default: use each marker ~once.  
  // h: help. print usage info
  fprintf(ostream, "Options: \n");
  fprintf(ostream, "  -i \t input file name (required).\n");
  fprintf(ostream, "  -r \t file name of reference data set (optional).\n");
  fprintf(ostream, "  -o \t output file name. Default: duplicatesearch.out\n");
  fprintf(ostream, "  -e \t maximum agmr; calculate agmr only if est. agmr is < this value. Default: 0.2\n");
  fprintf(ostream, "  -x \t maximum marker missing data fraction. Default: 2.0/chunk_size \n");
  fprintf(ostream, "  -a \t minimum minor allele frequency. Default: 0 \n");
  fprintf(ostream, "  -k \t chunk_size (number of markers per chunk). Default: set automatically depending on ploidy:\n");
  fprintf(ostream, "     \t\t ploidy=2 k=8; ploidy=4 k=5; ploidy=6 k=4.\n");
  fprintf(ostream, "  -s \t random number generator seed. Default: get seed from clock. \n");
   fprintf(ostream, "  -f \t maximum accession missing data fraction. Default: 0.5\n");
  fprintf(ostream, "  -v \t control output format. Default 1; 2 for more info.\n");
  // fprintf(ostream, "  -p \t ploidy. (default: 2)\n"); this now is automatically detected from the data.
  fprintf(ostream, "  -n \t number of chunks to use. Default: (int)n_markers/chunk_size \n");
  fprintf(ostream, "  -h \t print this usage information. \n");
}

void print_command_line(FILE* ostream, int argc, char** argv){
  fprintf(ostream, "# command line:  ");
  for(int i=0; i<argc; i++){
    fprintf(ostream, "%s  ", argv[i]);
  }fprintf(ostream, "\n");
}

// ************ unused ******************************

// given a vlong of matches to the gts so far in chunk (first j gts, say),
// take the next gt in query (i012) and add matching
// gts to get vlong of matches to first (j+1) gts.
// i.e. multiply each by 3, and add 0,1, or 2 as appropriate
// matching here <->  0 matched by 0, 2 by 2, 1 by anything (0,1,2) 
/* Vlong* matching_ipats(Vlong* ipats, long i012){ */
/*   if(i012 == 1){ */
/*     Vlong* new_vlong = construct_vlong(3*ipats->size); */
/*     for(long i=0; i<ipats->size; i++){ */
/*       add_long_to_vlong(new_vlong, 3*ipats->a[i] + 0); */
/*       add_long_to_vlong(new_vlong, 3*ipats->a[i] + 1); */
/*       add_long_to_vlong(new_vlong, 3*ipats->a[i] + 2); */
/*     } */
/*     free_vlong(ipats); */
/*     return new_vlong; */
/*   }else{ */
/*     if(DO_ASSERT) assert(i012 == 0  || i012 == 2); */
/*     for(long i=0; i<ipats->size; i++){ */
/*       ipats->a[i] = 3*ipats->a[i] + i012; */
/*     } */
/*     return ipats; */
/*   } */
/* } */

/* char* ipat_to_strpat(long k, long ipat){ // note: this generates a string which looks like the base-3 representation of ipat, EXCEPT in reverse order */
/*   // i.e. (for k=4) ipat=1 -> pat = '1000' */
/*   char* pat = (char*)malloc((k+1)*sizeof(char)); */
/*   if(ipat >= 0){ */
/*     for(long i=0; i<k; i++){ */
/*       pat[i] = 48 + ipat % 3; */
/*       ipat /= 3; */
/*     } */
/*   }else{ */
/*     for(long i=0; i<k; i++){ */
/*       pat[i] = 'X'; */
/*     } */
/*   } */
/*   pat[k] = '\0'; */
/*   // printf("ipat %ld   strpat: %s \n", ipat, pat); */
/*   return pat; */
/* } */

/* long strpat_to_ipat(long len, char* strpat){ // the inverse of ipat_to_strpat */
/*   long pat = 0; */
/*   long f = 1; */
/*   for(long j=0; j < len; j++){ */
/*     char a = strpat[j] - 48; */
/*     if((a>=0) && (a<=2)){ */
/*       pat += f*a; */
/*       f*=3; */
/*     }else{ */
/*       pat = -1; */
/*       break; */
/*     } */
/*   } */
/*   return pat; */
/* } */

/* Three_ds poly_agmr(Accession* gtset1, Accession* gtset2){ */
/*   char* gts1 = gtset1->genotypes->a; */
/*   char* gts2 = gtset2->genotypes->a; */
/*   long usable_pair_count = 0; // = agmr_denom */
/*   long mismatches = 0; // = agmr_numerator */
/*   long L1dist = 0; */
/*   double Lxdist = 0; */
/*   // fprintf(stderr, "strlen gts1, gts2: %ld %ld \n", strlen(gts1), strlen(gts2)); */
/*   for(long i=0; ;i++){ */
/*     char a1 = gts1[i]; */
/*     if(a1 == '\0') break; // end of  */
/*     char a2 = gts2[i]; */
/*     if(DO_ASSERT) assert(a2 != '\0'); */
/*     // fprintf(stderr, "chars: %c %c\n", a1, a2); */
/*     if(a1 != MISSING_DATA_CHAR){ */
/*       if(a2 != MISSING_DATA_CHAR){ */
/* 	usable_pair_count++; */
/* 	if(a1 != a2){ */
/* 	  mismatches++; */
/* 	  long dif = abs(a1 - a2); */
/* 	  L1dist += dif; */
/* 	  Lxdist += (dif > 2)? 2 : dif; */
/* 	} */
/*       } */
/*     } */
/*   } */
 
/*   Three_ds result; */
/*   if(usable_pair_count > 0){ */
/*     result = (Three_ds) */
/*       { .d1 = (double)mismatches/(double)usable_pair_count, // hamming dist normalized so max is 1 */
/* 	.d2 = (double)L1dist/(double)usable_pair_count, // L1 dist normalized so max is ploidy */
/* 	.d3 = (double)Lxdist/(double)usable_pair_count }; // sqrt((double)L2dist)/(double)usable_pair_count }; */
/*   }else{ */
/*     result = (Three_ds){-1, -1, -1}; */
/*   } */
/*   return result; */
/* } */




/* double agmr_hgmr(Accession* gtset1, Accession* gtset2, double* hgmr){ */
/*   char* gts1 = gtset1->genotypes->a; */
/*   char* gts2 = gtset2->genotypes->a; */
/*   long usable_pair_count = 0; // = agmr_denom */
/*   long mismatches = 0; // = agmr_numerator */
/*   long hgmr_denom = 0; */
/*   long hgmr_numerator = 0; */
/*   // fprintf(stderr, "strlen gts1, gts2: %ld %ld \n", strlen(gts1), strlen(gts2)); */
/*   for(long i=0; ;i++){ */
/*     char a1 = gts1[i]; */
/*     if(a1 == '\0') break; // end of  */
/*     char a2 = gts2[i]; */
/*     if(DO_ASSERT) assert(a2 != '\0'); */
/*     // fprintf(stderr, "chars: %c %c\n", a1, a2); */
/*     if(a1 != MISSING_DATA_CHAR){ */
/*       if(a2 != MISSING_DATA_CHAR){ */
/* 	usable_pair_count++; */
/* 	if(a1 != a2) mismatches++; */
/* 	if(a1 != '1' && a2 != '1'){ */
/* 	  hgmr_denom++; */
/* 	  if(a1 != a2) hgmr_numerator++; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   *hgmr = (hgmr_denom > 0)? (double)hgmr_numerator/hgmr_denom : -1; */
/*   return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1; */
/* } */


// *****  end of function definitions  *****
