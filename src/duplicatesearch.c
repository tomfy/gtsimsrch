// Chunk-wise search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <sys/sysinfo.h> // needed for get_nprocs
#include <pthread.h>
#include <assert.h>
#include <stdbool.h>
#include <errno.h>

#include "gtset.h"

#define DISTANCE_NORM_FACTOR (1.0) // if 0.5 max possible distance is 1 (if all dosage pairs are 0|2)
// defaults
#define DEFAULT_DIPLOID_CHUNK_SIZE  6
#define DEFAULT_MAX_DISTANCE  0.5
#define DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION  0.25
#define DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION  0.5
#define DEFAULT_MIN_MAF  0.1
// #define DEFAULT_D_RANDOM_SAMPLE_SIZE  20000
#define BITWISE true

//long nnnn = 0;
//double fcmc_time = 0;
//double fcmc_time2 = 0;
//double TXX = 0;
//double AGMR0;

//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // one of these for each chunk
  long capacity;
  long size;
  Vlong** a; // array of Vlong*; indices are patterns, values are Vlong*s containing indices of accessions with that pattern
} Pattern_idxs;

typedef struct{
  long capacity;
  long size; // number of chunks
  long chunk_size; // number of chunks used.
  long n_patterns; // (ploidy+1)^chunk_size
  Pattern_idxs** a; // array of Pattern_idxs, one element for each chunk
} Chunk_pattern_idxs;


typedef struct{
  long query_index;
  long match_index;
  double usable_chunks; // estimate or actual count
  long n_matching_chunks;
  double est_dist;
  //double dist;
  double agmr;
  double hgmr;
  double agmr0;
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci;

typedef struct{
  Vaccession* the_accessions;
  Chunk_pattern_idxs* the_cpi;
  long first_chunk;
  long last_chunk;
  double chunk_pair_md_count;
}threaded_pop_cpi_struct;

typedef struct{
  const GenotypesSet* the_genotypes_set;
  const Chunk_pattern_idxs* the_cpi;
  double max_est_dist;
  long Nthreads;
  long thread_number;

  Vmci** query_vmcis; // results stored here
  long distance_count; // the number of full distance calculations
  long n_matches_checked;
  double n_chunk_matches;
  double pqr_time;
  double agmr_time;
}TD;

int  do_checks = 0;

// *********************** function declarations ************************************************
void print_usage_info(FILE* ostream);
char* ipat_to_strpat(long len, long ipat); // unused
long strpat_to_ipat(long len, char* strpat); // unused
//ND distance(Accession* acc1, Accession* acc2);
four_longs agmr_hgmr_diploid(Accession* gtset1, Accession* gtset2);
//Vdouble* distances_random_sample(GenotypesSet* the_gtset, long n);


// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, double n_usable_chunks, long n_matching_chunks,
		   double est_dist, double agmr, double hgmr, double agmr0);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void push_to_vmci(Vmci* the_vmci, Mci* the_mci);
void sort_vmci_by_agmr(Vmci* the_vmci); // sort Vmci by distance
void sort_vmci_by_index(Vmci* the_vmci); // sort Vmci by index
int cmpmci_a(const void* v1, const void* v2);
int cmpmci_i(const void* v1, const void* v2);
void free_vmci(Vmci* the_vmci);


// *****  Pattern_idxs; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_idxs* construct_pattern_idxs(long n_patterns);
void free_pattern_idxs(Pattern_idxs*);

// *****  Chunk_pattern_idxs; indices are chunk numbers; elements are Pattern_idxs*
Chunk_pattern_idxs* construct_chunk_pattern_idxs(long n_chunks, long chunk_size, long ploidy);
// Vlong** get_all_match_counts(long n_accessions, Chunk_pattern_idxs* the_cpi);
double populate_chunk_pattern_idxs_from_vaccession(Vaccession* the_accessions, Chunk_pattern_idxs* the_cpi, long first_chunk, long last_chunk);
double count_chunk_matches(Chunk_pattern_idxs* the_cpi); // count the number of matching chunks among all accessions.
void print_chunk_pattern_idxs(Chunk_pattern_idxs* the_cpi, FILE* ostream);
void free_chunk_pattern_idxs(Chunk_pattern_idxs* the_cpi);

// *****  Gts and Chunk_pattern_idxs  ***********
long* find_chunk_match_counts(const Accession* the_gts, const Chunk_pattern_idxs* the_cpi); //, long first_match_idx, long last_match_idx);
Vmci** find_matches(const GenotypesSet* the_genotypes_set,
		    const Chunk_pattern_idxs* the_cpi, long Nthreads, double max_est_dist, long* distance_calculations_performed);

long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long out_format);
void print_command_line(FILE* ostream, int argc, char** argv);
// *****  functions to pass to pthread_create  ********
void* process_chunk_range(void* x);
void* process_query_range(void* x);

double clock_time(clockid_t a_clock){
  struct timespec tspec;
  clock_gettime(a_clock, &tspec);
  return (double)(tspec.tv_sec + 1.0e-9*tspec.tv_nsec);
}


// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  errno = 0;

  clockid_t clock1 = CLOCK_MONOTONIC;
  struct timespec tspec;
  double t_start = clock_time(clock1);
  unsigned rand_seed = time(0); // (unsigned)tspec.tv_nsec;

  
  long ploidy = 2; // Polyploid case not implemented
  long chunk_size = DEFAULT_DIPLOID_CHUNK_SIZE; // default: 6 - good for diploid.
  long n_passes = 1; // n_chunks = n_passes*(int)(n_markers/chunk_size) n_passes = 1 -> 
  
  double max_marker_missing_data_fraction = DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION; // if < 0, gets set to 2.0/chunk_size after chunk_size is set.
  double max_accession_missing_data_fraction = DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION; 
  double min_minor_allele_frequency = DEFAULT_MIN_MAF; //
  double max_est_dist = DEFAULT_MAX_DISTANCE; 
  long output_format = 1; // 1 ->  acc_id1 acc_id2 agmr hgmr; otherwise add 3 more columns: usable_chunks matching_chunks est_distance agmr0
  char default_output_filename[] = "duplicatesearch.out";
  bool print_filtered_gtset = false;
  bool shuffle_accessions = true; 

  long nprocs = (long)get_nprocs(); // returns 2*number of cores if hyperthreading.
  long Nthreads = (nprocs > 2)? nprocs/2 : 1; // default number of threads
  // long distance_random_sample_size = DEFAULT_D_RANDOM_SAMPLE_SIZE;
  
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
  while(1){
    int option_index = 0;
    static struct option long_options[] = {
      {"input",   required_argument, 0,  'i'}, // filename of new data set
      {"reference", required_argument, 0, 'r'}, // filename of reference data set (optional)
      {"output",  required_argument, 0,  'o'}, // output filename
      {"format", required_argument, 0, 'f'},
	
      {"marker_max_missing_data", required_argument, 0, 'm'}, // markers with > this fraction missing data will not be used.
      {"maf_min", required_argument, 0, 'l'}, //
      {"accession_max_missing_data", required_argument, 0, 'a'},
      {"distance_max",  required_argument,  0,  'd'}, // max. estimated distance (agmr)
      //	{"max_distance",  required_argument,  0,  'd'}, // min. 'estimated genotype probability'

      {"chunk_size", required_argument, 0, 'k'}, // number of markers per chunk. D
      {"passes", required_argument, 0, 'n'}, // use each marker in ~passes chunks		
      {"unshuffled",    no_argument, 0,  'u' }, // default is to shuffle the order of the accessions in output
      {"threads", required_argument, 0,  't'}, // number of threads to use
      {"seed", required_argument, 0, 's'}, // rng seed. Only relevant if shuffling.

      {"help", no_argument, 0, 'h'},
      {0,         0,                 0,  0 }
    };
     
    c = getopt_long_only(argc, argv, "", long_options, &option_index);
    if(c == -1) break;
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
    case 'd': 
      max_est_dist = (double)atof(optarg);
      if(max_est_dist <= 0){
	fprintf(stderr, "option d (distance_max) requires a numerical argument > 0 \n");
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
    case 'm':
      max_marker_missing_data_fraction = (double)atof(optarg);
      if(max_marker_missing_data_fraction <= 0){
	fprintf(stderr, "max_marker_missing_data_fraction < 0; will be set to 1.5/chunk_size\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'a':
      max_accession_missing_data_fraction = (double)atof(optarg);
      if(max_accession_missing_data_fraction <= 0){
	fprintf(stderr, "option f (max_accessions_missing_data_fraction) requires an real argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'l': 
      min_minor_allele_frequency = (double)atof(optarg);
      if(min_minor_allele_frequency < 0){
	fprintf(stderr, "option a (min_minor_allele_frequency) requires an real argument >= 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 't' :
      if(sscanf(optarg, "%ld", &Nthreads) != 1  ||  errno != 0){
	fprintf(stderr, "# Nthreads; conversion of argument %s to long failed.\n", optarg);
	exit(EXIT_FAILURE);
      }else if(Nthreads > nprocs){
	fprintf(stderr, "# Setting Nthreads to max allowed value of %ld.\n", nprocs);
	Nthreads = nprocs;
      }else if(Nthreads < 0){
	fprintf(stderr, "# Setting Nthreads to min allowed value of 1.\n");
	Nthreads = 1;
      }
      break;
    case 'f': 
      output_format = (long)atoi(optarg);
      if(output_format < 1){
	fprintf(stderr, "option v (output_format) requires an integer argument >= 1\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'u' :
      shuffle_accessions = false;
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
  fprintf(rparam_stream, "# Chunk size: %ld\n", chunk_size);
  fprintf(rparam_stream, "# Min. marker minor allele frequency: %5.3lf\n", min_minor_allele_frequency);
  if(max_marker_missing_data_fraction < 0) max_marker_missing_data_fraction = 1.5/chunk_size;
  fprintf(rparam_stream, "# Max. marker missing data fraction: %5.3lf\n", max_marker_missing_data_fraction);
  fprintf(rparam_stream, "# Max. accession missing data fraction: %5.3lf\n", max_accession_missing_data_fraction);
  if(max_est_dist > 0){
    fprintf(rparam_stream, "# Max. estimated distance: %5.3lf\n", max_est_dist);
  }else{
    fprintf(rparam_stream, "# Max. estimated distance will be set automatically.\n");
  }
  fclose(rparam_stream);
  fprintf(stdout, "%s", rparam_buf);
  fprintf(out_stream, "%s", rparam_buf);
  free(rparam_buf);
 
  // *****  done processing command line  *****************************************


  // *****  read in genotype data (including optionally a reference data set) and create a genotypesset object.
 
  GenotypesSet* the_genotypes_set = construct_empty_genotypesset(max_marker_missing_data_fraction, min_minor_allele_frequency, ploidy);
 
  if(reference_set_filename != NULL){ // load the reference set, if one was specified.
    fprintf(stdout, "# loading ref set.\n");
    add_accessions_to_genotypesset_from_file(reference_set_filename, the_genotypes_set, max_accession_missing_data_fraction);
    the_genotypes_set->n_ref_accessions = the_genotypes_set->accessions->size;
    fprintf(stdout, "# Done reading reference data set dosages from file %s. %ld accessions stored (%ld rejected); %ld markers.\n",
	    reference_set_filename, the_genotypes_set->n_ref_accessions, the_genotypes_set->n_bad_accessions, the_genotypes_set->n_markers);
  }
  fprintf(stdout, "# loading data set.\n");
  add_accessions_to_genotypesset_from_file(input_filename, the_genotypes_set, max_accession_missing_data_fraction); // load the new set of accessions
  fprintf(stdout, "# Done reading dosages from file %s. %ld accessions stored (%ld rejected); %ld markers.\n",
	  input_filename, the_genotypes_set->n_accessions, the_genotypes_set->n_bad_accessions, the_genotypes_set->n_markers);
  
  if(shuffle_accessions) shuffle_order_of_accessions(the_genotypes_set); // this helps to spread load evenly over threads 
    
  ploidy = the_genotypes_set->ploidy;
  if(ploidy != 2) { fprintf(stderr, "# Ploidy of %ld detected. Non-diploid ploidy not implemented. Exiting.\n", ploidy); exit(EXIT_FAILURE); }
  // double t_after_read = clock_time(clock1);
  // 'rectification' not needed.
  // rectify_markers(the_genotypes_set); // swap dosage 0 and 2 for markers with dosage 2 more common, so afterward 0 more common than 2 for all markers.
 
  filter_genotypesset(the_genotypes_set, out_stream);
  // double t_after_filter = clock_time(clock1);
  //  store_homozygs(the_genotypes_set); // homozygs are not needed for duplicatesearch
  //  double t_after_store_homozygs = clock_time(clock1);
  // fprintf(stderr, "# times. \n# to read: %8.5lf\n# to filter: %8.5lf\n", t_after_read-t_start, t_after_filter-t_after_read);
  // getchar();
  
  long n_accessions = the_genotypes_set->accessions->size;   
  long n_markers = the_genotypes_set->n_markers;

  long n_chunks_per_pass = n_markers/chunk_size;
  long n_chunks = n_passes*n_chunks_per_pass;
  fprintf(out_stream, "# Filtered data has %ld accessions, and %ld markers.\n", n_accessions, n_markers);
  fprintf(stdout, "# Chunk size: %ld  n_chunks: %ld\n", chunk_size, n_chunks);
  fprintf(out_stream, "# Chunk size: %ld  n_chunks: %ld\n", chunk_size, n_chunks);

  if(print_filtered_gtset){
    FILE* fh_gtsout = fopen("filtered_gtset.out", "w");
    print_genotypesset(fh_gtsout, the_genotypes_set);
    fclose(fh_gtsout);
  }

  double t_after_input = clock_time(clock1);
  fprintf(stdout, "# Time to load & filter dosage data: %6.3lf sec.\n", t_after_input - t_start);
  // *****  done reading, filtering, and storing input  **********

  check_genotypesset(the_genotypes_set);
  double t_after_chk = clock_time(clock1);
  fprintf(stdout, "# Time for check_genotypesset: %lf\n", t_after_chk - t_after_input);
  set_Abits_Bbits(the_genotypes_set);
  double t_after_set_ABbits = clock_time(clock1);
  fprintf(stdout, "# Time for set_Abits_Bbits: %lf\n", t_after_set_ABbits - t_after_chk);
  populate_marker_dosage_counts(the_genotypes_set);
  // set_agmr0s(the_genotypes_set);
  double t_after_pop_marker_dosage_counts = clock_time(clock1);
  fprintf(stdout, "# Time for populate_marker_dosage_counts: %lf\n", t_after_pop_marker_dosage_counts - t_after_set_ABbits);


  agmr0(the_genotypes_set); // calculate the overall agmr0 for the genotypes set.
  /* if(max_est_dist < 0){ // get max_est_dist so as to do approx. n_ds_to_get distance calculations */
  /*   // get random sample of distances */
  /*   distance_random_sample_size = 2*n_accessions; */
  /*   Vdouble* rand_distances = distances_random_sample(the_genotypes_set, distance_random_sample_size); */
  /*   Vdouble* sorted_rand_distances = sort_vdouble(rand_distances); */
  /*   double n_ds_all = 0.5*(n_accessions - the_genotypes_set->n_ref_accessions)*(n_accessions + the_genotypes_set->n_ref_accessions - 1); */
  /*   long max_dist_idx = (long)((n_ds_to_get/n_ds_all)*distance_random_sample_size); */
  /*   fprintf(stderr, "# %ld  %ld  %ld\n", n_ds_to_get, distance_random_sample_size, (long)n_ds_all); */
  /*   if (max_dist_idx >= sorted_rand_distances->size) max_dist_idx = sorted_rand_distances->size-1; */
  /*   fprintf(stderr, "# %ld %ld \n", max_dist_idx, sorted_rand_distances->size); */
  /*   max_est_dist = sorted_rand_distances->a[max_dist_idx]; */
  /*   fprintf(stdout, "# Max. estimated distance: %8.5f\n", max_est_dist); */
  /*   fprintf(out_stream, "# Max. estimated distance: %8.5f\n", max_est_dist); */
  /*   free_vdouble(sorted_rand_distances); */
  /* } */


  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices);
  marker_indices->size = n_chunks_per_pass*chunk_size; //
  for(long i=1; i<n_passes; i++){
    Vlong* more_marker_indices = construct_vlong_whole_numbers(n_markers);
    shuffle_vlong(more_marker_indices);
    more_marker_indices->size = n_chunks_per_pass*chunk_size;
    append_vlong_to_vlong(marker_indices, more_marker_indices);
  }

  Vaccession* the_accessions = the_genotypes_set->accessions;
  set_vaccession_chunk_patterns(the_accessions, marker_indices, n_chunks, chunk_size, ploidy);
  Chunk_pattern_idxs* the_cpi = construct_chunk_pattern_idxs(n_chunks, chunk_size, ploidy);
  double chunk_pair_md_count = 0;
 
  if(Nthreads == 0){ // unthreaded
    chunk_pair_md_count += populate_chunk_pattern_idxs_from_vaccession(the_accessions, the_cpi, 0, n_chunks-1);
  }else{ // use 1 or more pthreads
    long first_chunk = 0;
    pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
    threaded_pop_cpi_struct* tpcs = (threaded_pop_cpi_struct*)malloc(Nthreads*sizeof(threaded_pop_cpi_struct));
    for(long i_thread = 0; i_thread < Nthreads; i_thread++){
      tpcs[i_thread].the_accessions = the_accessions;
      tpcs[i_thread].the_cpi = the_cpi;
      tpcs[i_thread].first_chunk = first_chunk;
      tpcs[i_thread].last_chunk = (i_thread == Nthreads-1)? n_chunks-1 : first_chunk + n_chunks/Nthreads - 1;
      first_chunk += n_chunks/Nthreads;
    } 
 
    for(long i_thread = 0; i_thread < Nthreads; i_thread++){
      int iret = pthread_create(thrids+i_thread, NULL, process_chunk_range, (void*) (tpcs+i_thread));
      if(iret > 0) fprintf(stderr, "# warning. pthread_create returned non-zero value. Thread %ld \n", (long)thrids[i_thread]);
    }

    for(long i_thread=0; i_thread<Nthreads; i_thread++){
      pthread_join(thrids[i_thread], NULL);
      chunk_pair_md_count += tpcs[i_thread].chunk_pair_md_count;
    }
  }
  
  
  double t_after_cpi = clock_time(clock1);
  fprintf(stdout, "# Time to construct & populate chunk_pattern_idxs: %6.3f\n", t_after_cpi - t_after_pop_marker_dosage_counts);
  //  exit(0);
  
   double chunk_match_count = count_chunk_matches(the_cpi);
   double total_possible_chunk_matches = n_chunks * 0.5*(double)n_accessions*(n_accessions-1);
   /* fprintf(stdout, "# There are %5.3e chunk matches of possible %5.3e; ratio: %6.2f\n", */
   /* 	   chunk_match_count, total_possible_chunk_matches, total_possible_chunk_matches/chunk_match_count); */
   /* fprintf(stdout, "# Total chunk match count: %12.4f; time: %8.5f\n", chunk_match_count, clock_time(clock1)-t_after_cpi); */

   fprintf(stdout, "# There are %6.4e possible chunk matches, of which %6.4e (%6.3f %%) are lost due to missing data,\n",
	   total_possible_chunk_matches, chunk_pair_md_count, 100*(chunk_pair_md_count/total_possible_chunk_matches));
   fprintf(stdout, "#   leaving %6.4e of which %6.4e (%6.3f %% ) actually match.\n",
	   total_possible_chunk_matches-chunk_pair_md_count, chunk_match_count,
	   100*(chunk_match_count/(total_possible_chunk_matches-chunk_pair_md_count)));
  

   long distance_calculations_performed;
   Vmci** query_vmcis = find_matches(the_genotypes_set, the_cpi, Nthreads, max_est_dist, &distance_calculations_performed); //, n_maf_categories, maf_category_marker_indices);
  double t_after_find_matches = clock_time(clock1);
  fprintf(stdout, "# Time to find %ld candidate matches and distances: %6.3f\n",
	  distance_calculations_performed, t_after_find_matches - t_after_cpi);
  long output_pairs_count = print_results(the_accessions, query_vmcis, out_stream, output_format);
  fprintf(stdout, "# Number of accession pairs output: %ld\n", output_pairs_count);
  fclose(out_stream);

  // *****  clean up  *****
  long cume_s = 0; 
  for(long i=0; i< the_accessions->size; i++){
    long s = query_vmcis[i]->size;
    cume_s += s;
    free_vmci(query_vmcis[i]);
  }
  //fprintf(stderr, "# TXX: %8.6f\n", TXX);
 
  free(query_vmcis);
  free_genotypesset(the_genotypes_set); 
  free_vlong(marker_indices);
  free_chunk_pattern_idxs(the_cpi);
  fprintf(stderr, "# total duplicatesearch run time: %9.3f\n", clock_time(clock1) - t_start);
  exit(EXIT_SUCCESS);
}

// **********************************************************************************************
// **********************  end of main  *********************************************************
// **********************************************************************************************


// *******************  function definitions  ***************************************************

// *****  Vaccession  ***********************************************************

double populate_chunk_pattern_idxs_from_vaccession(Vaccession* the_accessions, Chunk_pattern_idxs* the_cpi, long first_chunk, long last_chunk){
  long n_patterns = the_cpi->n_patterns;
 
  for(long i_gts=0; i_gts<the_accessions->size; i_gts++){
    Accession* the_gts = the_accessions->a[i_gts];
    if(DO_ASSERT) assert(i_gts == the_gts->index);
    Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts 
    
    // for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
    for(long i_chunk = first_chunk; i_chunk <= last_chunk; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      if(DO_ASSERT) assert(the_pat >= 0);

      if(the_cpi->a[i_chunk]->a[the_pat] == NULL){
	the_cpi->a[i_chunk]->a[the_pat] = construct_vlong(1);
      }
      Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];	
      push_to_vlong(the_accidxs, the_gts->index);
    }
  }
  
  long total_mdchunk_count = 0;
  double total_md_chunk_matches = 0;
  for(long i=first_chunk; i<=last_chunk; i++){
    long chunk_md_count = (the_cpi->a[i]->a[n_patterns] == NULL)?
      0 :  // there are no accessions having missing data for this chunk
      the_cpi->a[i]->a[n_patterns]->size;
    total_mdchunk_count += chunk_md_count;
    total_md_chunk_matches += (double)chunk_md_count*(double)(the_accessions->size - 0.5*(chunk_md_count+1));
  }
  //fprintf(stderr, "total_md_chunk_matches: %12.0f\n", total_md_chunk_matches);
  return total_md_chunk_matches;
}

double count_chunk_matches(Chunk_pattern_idxs* the_cpi){ // count the number of matching chunks among all accessions.
  // this should a small fraction
  if(DO_ASSERT) assert(the_cpi->a[0]->size == the_cpi->n_patterns);
  double cm_count = 0;
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    for(long i_pat=0; i_pat<the_cpi->n_patterns; i_pat++){
      long n_cp_matches = (the_cpi->a[i_chunk]->a[i_pat] != NULL)? the_cpi->a[i_chunk]->a[i_pat]->size : 0;
      cm_count += (n_cp_matches*(n_cp_matches-1))/2;
    }
  }
  return cm_count; 
}

// *****  Pattern_idxs; indices are patterns; elements are Vlong* of accidxs having that pattern.

Pattern_idxs* construct_pattern_idxs(long n_patterns){
  Pattern_idxs* pat_ids = (Pattern_idxs*)malloc(1*sizeof(Pattern_idxs));
  pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=(ploidy+1)^k) good patterns, and index n_pattern is for the missing data case.
  pat_ids->size = n_patterns+1;
  pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(long ipat=0; ipat< pat_ids->size; ipat++){
    pat_ids->a[ipat] = NULL; // construct_vlong(1); // waste of memory? set to NULL until needed?
  }
  return pat_ids;
}

void free_pattern_idxs(Pattern_idxs* pat_ids){
  if(pat_ids == NULL) return;
  for(long i=0; i<pat_ids->size; i++){
    if(pat_ids->a[i] != NULL) free_vlong(pat_ids->a[i]);
  }
  free(pat_ids->a);
  free(pat_ids);
}

// *****  Chunk_pattern_idxs; indices are chunk numbers; elements are Pattern_idxs*

Chunk_pattern_idxs* construct_chunk_pattern_idxs(long n_chunks, long chunk_size, long ploidy){ // needed size is known at construct time, so one param for both cap and size
  Chunk_pattern_idxs* chunk_pat_ids = (Chunk_pattern_idxs*)malloc(1*sizeof(Chunk_pattern_idxs));
  chunk_pat_ids->capacity = n_chunks;
  chunk_pat_ids->size = n_chunks;
  chunk_pat_ids->chunk_size = chunk_size;
  long n_patterns = int_power(ploidy+1, chunk_size);
  chunk_pat_ids->n_patterns = n_patterns;
  chunk_pat_ids->a = (Pattern_idxs**)malloc(n_chunks*sizeof(Pattern_idxs*));
  for(int i=0; i< chunk_pat_ids->size; i++){
    chunk_pat_ids->a[i] = construct_pattern_idxs(n_patterns);
  }
  return chunk_pat_ids;
}

void free_chunk_pattern_idxs(Chunk_pattern_idxs* the_cpi){
  if(the_cpi == NULL) return;
  for(long i=0; i< the_cpi->size; i++){
    free_pattern_idxs(the_cpi->a[i]);
  }
  free(the_cpi->a);
  free(the_cpi);
}

void print_chunk_pattern_idxs(Chunk_pattern_idxs* the_cpi, FILE* ostream){
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    fprintf(ostream, "i_chunk: %ld\n", i_chunk);
    Pattern_idxs* the_pi = the_cpi->a[i_chunk];
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

// *****  Accession and Chunk_pattern_idxs  ***********

long* find_chunk_match_counts(const Accession* the_accession, const Chunk_pattern_idxs* the_cpi){ //, long first_match_idx, long last_match_idx){
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_patterns = the_accession->chunk_patterns; // chunk patterns for Accession the_accession (i.e. the query accession)
  // Vlong* accidx_matchcounts = construct_vlong_zeroes(the_accession->index);
  long* accidx_matchcounts = (long*)calloc(the_accession->index, sizeof(long));
 
  for(long i_chunk=0; i_chunk < chunk_patterns->size; i_chunk++){
    long the_pattern = chunk_patterns->a[i_chunk];  
    

    if(DO_ASSERT) assert(the_pattern >= 0  &&  the_pattern <= n_patterns);
    // (patterns 0..n_patterns-1 are good, the_pattern == n_patterns (== (ploidy+1)^chunk_size) indicates a chunk with missing data )
    if(the_pattern != n_patterns){ // no missing data in this chunk, i.e. the_pattern = 0..n_patterns-1 (good data)   
      // just get the counts for matches with index < index of the_accession      
      // since those are the only ones used in find_matches. (slightly faster)
      Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pattern]; // array of accession indices of the matches to this chunk & pat
      //   for(long i=0; i<chunk_match_idxs->size; i++){ //
      for(long i=0; chunk_match_idxs->a[i] < the_accession->index; i++){
	  //	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	  // if(accidx >= the_accession->index) break;
	//	accidx_matchcounts->a[accidx]++; // accidx here is < the_accession->index
	accidx_matchcounts[ chunk_match_idxs->a[i] ]++;
      }
    }
  }
  return accidx_matchcounts; // array containing the number of matching chunks for each accession (with lower index than the_accession->index)
}

// ***** Mci  *****

Mci* construct_mci(long qidx, long midx, double usable_chunks, long n_matching_chunks,
		   // double est_matching_chunk_fraction, double matching_chunk_fraction){
		   double est_dist, double agmr, double hgmr, double agmr0){
  Mci* the_mci = (Mci*)calloc(1,sizeof(Mci));
  the_mci->query_index = qidx;
  the_mci->match_index = midx;
  the_mci->usable_chunks = usable_chunks;
  the_mci->n_matching_chunks = n_matching_chunks;
  the_mci->est_dist = est_dist;
  //  the_mci->dist = dist;
  the_mci->agmr = agmr;
  the_mci->hgmr = hgmr;
  the_mci->agmr0 = agmr0;
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

void push_to_vmci(Vmci* the_vmci, Mci* mci){
  long cap = the_vmci->capacity;
  long n = the_vmci->size;
  if(n == cap){   // if necessary, resize w realloc
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


int cmpmci_i(const void* v1, const void* v2){  // sort by query; distance is tie breaker
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

int cmpmci_a(const void* v1, const void* v2){ // sort by distance, query index is tie breaker
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

void sort_vmci_by_agmr(Vmci* the_vmci){ // sort by distance (low to high), query index as tie-breaker
  qsort(the_vmci->a, the_vmci->size, sizeof(Mci*), cmpmci_a);
}

void sort_vmci_by_index(Vmci* the_vmci){ // sort by query index
  qsort(the_vmci->a, the_vmci->size, sizeof(Mci*), cmpmci_i);
}

// *********************************************
// *********************************************

Vmci** find_matches(const GenotypesSet* the_genotypes_set,
		    const Chunk_pattern_idxs* the_cpi,
		    long Nthreads,
		    double max_est_dist, long* distance_calculations_performed)
{
  Vaccession* the_accessions = the_genotypes_set->accessions;
  long n_ref_accessions = the_genotypes_set->n_ref_accessions;
  long n_queries = the_accessions->size - n_ref_accessions;
  long n_markers = the_accessions->a[0]->genotypes->length;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;

  double min_matching_chunk_fraction = pow(1.0 - 2*max_est_dist, chunk_size);
  Vmci** query_vmcis = (Vmci**)malloc(the_accessions->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_accessions->size; i++){
    query_vmcis[i] = construct_vmci(4);
  }

  if(Nthreads == 0){ // i.e. no pthreads
    TD* td = (TD*)malloc(1*sizeof(TD));
    td[0].thread_number = 0;
    td[0].the_genotypes_set = the_genotypes_set;
    td[0].the_cpi = the_cpi;
    td[0].max_est_dist = max_est_dist;
    td[0].Nthreads = 1;

    td[0].query_vmcis = query_vmcis;
    process_query_range((void*) td);

     fprintf(stdout, "# Time to check %ld potential matches, and calculate %ld distances: %6.4f sec.\n",
	      td[0].n_matches_checked, td[0].distance_count, td[0].pqr_time);
    *distance_calculations_performed = td[0].distance_count;
    free(td);
  }else{ // use one or more pthreads
    long full_distance_calculations_performed = 0;
    pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
    TD* td = (TD*)malloc(Nthreads*sizeof(TD));
    fprintf(stdout, "# Nthreads = %ld\n", Nthreads);
    for(long i_thread = 0; i_thread < Nthreads; i_thread++){
      td[i_thread].thread_number = i_thread;
      td[i_thread].the_genotypes_set = the_genotypes_set;
      td[i_thread].the_cpi = the_cpi;
      td[i_thread].max_est_dist = max_est_dist;
      td[i_thread].Nthreads = Nthreads;
	
      td[i_thread].query_vmcis = query_vmcis;
    }
 
    for(long i_thread = 0; i_thread < Nthreads; i_thread++){
      int iret = pthread_create(thrids+i_thread, NULL, process_query_range, (void*) (td+i_thread));
      if(iret > 0) fprintf(stderr, "# warning. pthread_create returned non-zero value. Thread %ld \n", (long)thrids[i_thread]);
    }

    for(long i_thread=0; i_thread<Nthreads; i_thread++){
      pthread_join(thrids[i_thread], NULL);
      /* fprintf(stdout, "# Thread %ld; potential matches checked: %ld in %6.4f sec; distance calculations: %ld in %6.4f sec.\n", */
      /* 	      i_thread, td[i_thread].n_matches_checked, td[i_thread].pqr_time-td[i_thread].agmr_time, */
      /* 	      td[i_thread].distance_count, td[i_thread].agmr_time); */

      fprintf(stdout, "# Thread: %ld. Time to check %ld potential matches, and calculate %ld distances: %6.4f sec.\n",
	      i_thread, td[i_thread].n_matches_checked, td[i_thread].distance_count, td[i_thread].pqr_time);
      full_distance_calculations_performed += td[i_thread].distance_count;
    }
    *distance_calculations_performed = full_distance_calculations_performed;
    // fprintf(stdout, "# Total distance calculations: %ld\n", full_distance_calculations_performed);
    free(td);
    free(thrids);
  } // end >= 1 pthreads
  return query_vmcis;
}

void* process_chunk_range(void* x){
  threaded_pop_cpi_struct* tpcs = (threaded_pop_cpi_struct*)x;  
  tpcs->chunk_pair_md_count = populate_chunk_pattern_idxs_from_vaccession(tpcs->the_accessions, tpcs->the_cpi, tpcs->first_chunk, tpcs->last_chunk);
  // fprintf(stderr, "thread chunk_pair_md_count: %8.4f\n", tpcs->chunk_pair_md_count);
}

void* process_query_range(void* x){
  clockid_t clock2 = CLOCK_MONOTONIC;
  double start_pqr_time = clock_time(clock2);
  double true_distance_time = 0;
  TD* td = (TD*)x;
  const GenotypesSet* the_genotypes_set = td->the_genotypes_set;
  const Chunk_pattern_idxs* the_cpi = td->the_cpi;
  double max_est_dist = td->max_est_dist; 
  Vmci** query_vmcis = td->query_vmcis;

  Vaccession* the_accessions = the_genotypes_set->accessions;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;

  td->n_matches_checked = 0;
  td->distance_count = 0;

  for(long i_query = td->the_genotypes_set->n_ref_accessions + td->thread_number; i_query < the_accessions->size; i_query += td->Nthreads){
      
    Accession* q_gts = the_accessions->a[i_query];
     double agmr_nought = the_genotypes_set->agmr0;
    double min_matching_chunk_fraction = (max_est_dist*agmr_nought < 1.0)? pow(1.0 - max_est_dist*agmr_nought, chunk_size) : 0.0;
    double q_ok_chunk_fraction_x_mmcf = min_matching_chunk_fraction*(double)q_gts->ok_chunk_count/(double)n_chunks;
   
    long i_match_max = i_query-1; // (td->triangle)? i_query-1 : td->last_match_idx;
    long* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi); //, td->first_match_idx, i_match_max);
    for (long i_match = 0; i_match <= i_match_max; i_match++){ // compare the query just to other accessions with smaller indices.
      // est. usable chunks: (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks;
      //  double min_matching_chunk_count = q_ok_chunk_fraction_x_mmcf*the_accessions->a[i_match]->ok_chunk_count;
      if( chunk_match_counts[i_match] >= q_ok_chunk_fraction_x_mmcf*the_accessions->a[i_match]->ok_chunk_count ){ // min_matching_chunk_count
	//	double predistance_time = clock_time(clock2);
	double agmr, hgmr;
	if(BITWISE){ // bitwise ~20x faster
	  four_longs bfcs = bitwise_agmr_hgmr(q_gts, the_accessions->a[i_match]);
	  // bfcs: {count_02, count_00_22, count_01_12, count_11}
	  long b_agmr_num = bfcs.l1 + bfcs.l3;
	  long b_agmr_denom = b_agmr_num + bfcs.l2 + bfcs.l4;
	  long b_hgmr_num = bfcs.l1;
	  long b_hgmr_denom = bfcs.l1 + bfcs.l2;
	  agmr = (b_agmr_denom > 0)? (double)b_agmr_num / (double)b_agmr_denom : -1;
	  hgmr = (b_hgmr_denom > 0)? (double)b_hgmr_num / (double)b_hgmr_denom : -1;	
	}else{ // straightforward (slow) way
	  four_longs four_counts = agmr_hgmr_diploid(q_gts, the_accessions->a[i_match]);
	  // four_counts: {count_02, count_00_22, count_01_12, count_11};
	  long agmr_numerator = four_counts.l1 + four_counts.l3;
	  //  long d_numerator = agmr_numerator + four_counts.l1;
	  long d_denominator = agmr_numerator + four_counts.l2 + four_counts.l4;
	  agmr = (d_denominator > 0)?  // DISTANCE_NORM_FACTOR*(double)agmr_numerator/(double)d_denominator : -1;
	     DISTANCE_NORM_FACTOR*(double)agmr_numerator/(double)d_denominator : -1;
	  long hgmr_numerator = four_counts.l1;
	  long hgmr_denominator = hgmr_numerator + four_counts.l2;
	  hgmr = (hgmr_denominator > 0)? (double)hgmr_numerator/(double)hgmr_denominator : -1;
	}
	//	fprintf(stderr, "# agmr: %8.5f  hgmr: %8.5f \n", agmr, hgmr);
	//	agmr_nought = pow(agmr_nought * agmr0_accvsallx(the_genotypes_set, the_accessions->a[i_match]), 0.5); // too slow
	//	double agmr_norm = agmr/agmr_nought; // agmr_nought is typical dist for random pair of this data set.
		td->distance_count++;

		//	double postdistance_time = clock_time(clock2);
		//	true_distance_time += postdistance_time - predistance_time;
	if(agmr <= max_est_dist*agmr_nought){
	   long matching_chunk_count = chunk_match_counts[i_match];
	   double usable_chunk_count = (double)q_gts->ok_chunk_count * the_accessions->a[i_match]->ok_chunk_count/(double)n_chunks;
	   //	     min_matching_chunk_count/min_matching_chunk_fraction;
	   double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	  double est_dist = DISTANCE_NORM_FACTOR*(1.0 - pow(matching_chunk_fraction, 1.0/chunk_size)); // /agmr_nought;
	  //	  double agmr_nought_b = agmr0_accvsall(the_genotypes_set, the_accessions->a[i_match]);
	  push_to_vmci(query_vmcis[i_query],
		       construct_mci(i_query, i_match, usable_chunk_count, matching_chunk_count, est_dist, agmr, hgmr, agmr_nought));
	
	} // end if(true_dist < max_est_dist)
      } // end if(enough matching chunks) - i.e. est dist is small enough
    } // end loop over potential matches to query
  
    td->n_matches_checked += i_query; // (blockwise)? ((td->triangle)? i_query - td->first_match_idx : (td->last_match_idx - td->first_match_idx + 1)) : i_query;
    // free_vlong(chunk_match_counts);
    free(chunk_match_counts);
  } // end loop over queries.
  td->pqr_time = clock_time(clock2)  - start_pqr_time;
  td->agmr_time = true_distance_time;
}

long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long output_format){
  long distance_count = 0;

  fprintf(ostream, "#id1  id2  agmr  hgmr  ");
  if(output_format != 1){
    fprintf(ostream, "usable_chunks matching_chunks est_agmr normalized_agmr");
  }
  fprintf(ostream, "\n");
    
  for(long i_q=0; i_q<the_accessions->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    sort_vmci_by_agmr(the_vmci); // sort in ascending agmr order.   
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      Mci* the_mci = the_vmci->a[i_m];      
      Accession* q_acc = the_accessions->a[i_q];
      Accession* m_acc = the_accessions->a[the_mci->match_index];
      fprintf(ostream, "%s  %s  %8.6f %8.6f", q_acc->id->a, m_acc->id->a, the_mci->agmr, the_mci->hgmr);
      if(output_format != 1){
	double agmr_norm = (the_mci->agmr0 > 0)? the_mci->agmr/the_mci->agmr0 : -1;
	fprintf(ostream, "  %6.2f %ld %7.5f %7.5f ",
		the_mci->usable_chunks, the_mci->n_matching_chunks, the_mci->est_dist, agmr_norm); //the_mci->agmr0);
      }
      fprintf(ostream, "\n");
      distance_count++;
    } // end of loop over matches to query
  } // end of loop over queries
  return distance_count;
}

void print_usage_info(FILE* ostream){
  // i: input file name (required).
  // r: reference set file name.
  // o: output file name. Default: duplicatesearch.out
  // e: max estimated distance. default: DEFAULT_MAX_DISTANCE (Calculate distance only if quick est. is < this value.) 
  // x: marker max missing data fraction
  // a: min minor allele frequency
  // k: chunk size (number of markers per chunk). Default: 8
  // s: random number seed. Default: get seed from clock.
  // f: output format control. Default: 1; 2 gives additional info.
  // n: number of chunks to use. Default: use each marker ~once.  
  // h: help. print usage info
  fprintf(ostream, "Options: \n");
  fprintf(ostream, "  -i  -input       <input file name> (required).\n");
  fprintf(ostream, "  -r  -reference   <file name of reference data set> (optional).\n");
  fprintf(ostream, "  -o  -output <output file name>  Default: duplicatesearch.out\n");
  fprintf(ostream, "  -f  -format  output format. Default: 1; 2 for more info.\n");
    
  fprintf(ostream, "  -d  -distance_max <number> calculate distance only if est. distance <= this value. Default: %5.3f\n", DEFAULT_MAX_DISTANCE);
  fprintf(ostream, "  -m  -marker_max_missing_data <number> maximum marker missing data fraction. Default: %5.3f\n", DEFAULT_MAX_MARKER_MISSING_DATA_FRACTION);
  fprintf(ostream, "  -l  -maf_min <number>  minimum minor allele frequency. Default: %5.3F \n", DEFAULT_MIN_MAF);
  fprintf(ostream, "  -a  -accession_max_missing_data <number> maximum accession missing data fraction. Default: %5.3f\n", DEFAULT_MAX_ACCESSION_MISSING_DATA_FRACTION);
  
  fprintf(ostream, "  -k  -chunk_size <integer> number of markers per chunk). Default: %d\n", DEFAULT_DIPLOID_CHUNK_SIZE);
  fprintf(ostream, "  -s  -seed <integer> random number generator seed. Default: get seed from clock. \n");
  fprintf(ostream, "  -n  -passes <integer> The number of chunks to use is  (int)n_markers/chunk_size \n");
  fprintf(ostream, "  -u  -unshuffle  Leave the order of accessions as in input file. Default is to shuffle the order.\n");
  fprintf(ostream, "  -t  -threads <integer> The number of threads to use. Default is to set automatically based on the number of cores.\n");
  fprintf(ostream, "  -h  -help  print this usage information. \n");
}

void print_command_line(FILE* ostream, int argc, char** argv){
  fprintf(ostream, "# command:  ");
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
/*       push_to_vlong(new_vlong, 3*ipats->a[i] + 0); */
/*       push_to_vlong(new_vlong, 3*ipats->a[i] + 1); */
/*       push_to_vlong(new_vlong, 3*ipats->a[i] + 2); */
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




four_longs agmr_hgmr_diploid(Accession* gtset1, Accession* gtset2){
  char* gts1 = gtset1->genotypes->a;
  char* gts2 = gtset2->genotypes->a;
  long count_00_22 = 0; // count of 00, 22
  long count_02 = 0; // count of 02, 20
  long count_11 = 0; // count of 11
  long count_01_12 = 0; // count of 01, 10, 12, 21
  // agmr_numerator = count_02 + count_01_12
  // hgmr_numerator = count_02;
  // agmr_denominator = agmr_numerator + count_00_22 + count_11
  // hgmr_denominator = hgmr_numerator + count_00_22
  
  // fprintf(stderr, "strlen gts1, gts2: %ld %ld \n", strlen(gts1), strlen(gts2));
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break; // end of
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0');
    // fprintf(stderr, "chars: %c %c\n", a1, a2);
    if(a1 != MISSING_DATA_CHAR){ // skip this marker
      if(a2 != MISSING_DATA_CHAR){ // skip this marker
	//	if(a1 != a2) mismatches++;
	if(a1 != '1' && a2 != '1'){ // both homozyg
	  if(a1 != a2) { // not equal: 02 or 20
	    count_02++;
	  }else{
	    count_00_22++;
	  }
	}else{ // at least one heterozyg
	  if(a1 == a2){ 
	    count_11++;
	  }else{
	    count_01_12++;
	  }
	}
      }
    }
  }
  four_longs result // = {count_00_22, count_02, count_11, count_01_12};
  = {count_02, count_00_22, count_01_12, count_11};
  return result;
}


/* long print_results_a(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream, long output_format){ */
/*   Vmci* all_mcis = construct_vmci(1000); */
/*   for(long i_q=0; i_q<the_accessions->size; i_q++){ */
/*     Vmci* the_vmci = query_vmcis[i_q]; */
/*     for(long i_m=0; i_m < the_vmci->size; i_m++){ */
/*       // Mci* the_mci = the_vmci->a[i_m]; */
/*       push_to_vmci(all_mcis, the_vmci->a[i_m]); */
/*     } */
/*   } */
/*   long true_agmr_count = 0; */
/*   // sort_vmci_by_agmr(all_mcis); */
/*   sort_vmci_by_index(all_mcis); */
/*   for(long i=0; i< all_mcis->size; i++){ */
/*     Mci* the_mci = all_mcis->a[i]; */
/*     Accession* q_acc = the_accessions->a[the_mci->query_index]; */
/*     Accession* m_acc = the_accessions->a[the_mci->match_index]; */
 
/*     fprintf(ostream, "%4ld  %30s  %30s  %5.2f  %3ld  %5.3f  %5.3f", */
/* 	    the_mci->query_index,  q_acc->id->a,   m_acc->id->a,   */
/* 	    the_mci->usable_chunks,  the_mci->n_matching_chunks, */
/* 	    the_mci->est_dist,  the_mci->agmr); */
/*     if (output_format == 1){ */
/*       // leave as is */
/*     }else if(output_format == 2){ // add a bit more info */
/*       fprintf(ostream, "  %4ld  %4ld  %3ld  %3ld", */
/* 	      q_acc->missing_data_count,  q_acc->md_chunk_count, m_acc->missing_data_count, m_acc->md_chunk_count); */
/*     }else{ */
/*       fprintf(stderr, "# output_format %ld is unknown. using default output_format.\n", output_format); */
/*     } */
/*     fprintf(ostream, "\n"); */
/*     true_agmr_count++; */
/*   } // end of loop over matches to query */
/*   free_vmci(all_mcis); */
/*   return true_agmr_count; */
/* } // end print_results_a */




/* Vlong** get_maf_cat_marker_indices(GenotypesSet* the_genotypes_set, long n_maf_categories){ */
  
/*   Vdouble* sorted_mafs = sort_vdouble(copy_vdouble(get_minor_allele_frequencies(the_genotypes_set))); */
/*   Vdouble* maf_threshholds = construct_vdouble(n_maf_categories); */
/*   for(long i=1; i<n_maf_categories; i++){ */
/*     long quantile_index = sorted_mafs->size*i/n_maf_categories; */
/*     double quantile_maf = 0.5*(sorted_mafs->a[quantile_index] + sorted_mafs->a[quantile_index+1]); */
/*     push_to_vdouble(maf_threshholds, quantile_maf); */
/*     fprintf(stderr, "quantile: %lf index: %ld  maf: %lf \n", (double)i/n_maf_categories, quantile_index, quantile_maf); */
/*   } */
/*   push_to_vdouble(maf_threshholds, 1.0); */
  
/*   FILE* fh = fopen("mafs.out", "w"); */
/*   for(long i=0; i<sorted_mafs->size; i++){ */
/*     fprintf(fh, "%lf\n", sorted_mafs->a[i]); */
/*   } */
/*   fclose(fh); */

/*   //  ******   Now get arrays of indices of markers in the maf categories: */
/*   Vlong** maf_category_marker_indices = (Vlong**)malloc(n_maf_categories*sizeof(Vlong*)); */
  
/*   for(long i=0; i<n_maf_categories; i++){ //  */
/*     maf_category_marker_indices[i] = construct_vlong(the_genotypes_set->n_markers); */
/*   }    */
/*   Vlong* missing_data_counts = the_genotypes_set->marker_missing_data_counts; */
/*   Vlong* minor_allele_counts = the_genotypes_set->marker_alt_allele_counts; */
/*   if(DO_ASSERT) assert(missing_data_counts->size == minor_allele_counts->size); */
/*   for(long i_marker=0; i_marker<minor_allele_counts->size; i_marker++){ */
/*     long ok_count = the_genotypes_set->n_accessions - missing_data_counts->a[i_marker]; */
/*     if(ok_count > 0){ */
/*       double marker_minor_allele_frequency = (double)minor_allele_counts->a[i_marker]/(double)(2*ok_count); */
/*       for(long j_maf=0; j_maf<n_maf_categories; j_maf++){ */
/* 	if(marker_minor_allele_frequency <= maf_threshholds->a[j_maf]){ */
/* 	  push_to_vlong(maf_category_marker_indices[j_maf], i_marker); */
/* 	  break; */
/* 	} */
/*       } // end loop over maf ranges */
/*     } */
/*   } */

/*   free_vdouble(sorted_mafs); */
/*   free_vdouble(maf_threshholds); */
/*   return maf_category_marker_indices; */
/* } */

/* Vlong** get_maf_cat_marker_indices_x(GenotypesSet* the_genotypes_set){ */
  
/*   Vdouble* sorted_mafs = sort_vdouble(copy_vdouble(get_minor_allele_frequencies(the_genotypes_set))); */
/*   Vdouble* maf_threshholds = construct_vdouble(10); */
/*   /\* maf_threshholds->a[0] = 0.005; *\/ */
/*   /\*  maf_threshholds->a[1] = 0.01; *\/ */
/*   /\*   maf_threshholds->a[2] = 0.02; *\/ */
/*   /\*    maf_threshholds->a[3] = 0.04; *\/ */
/*   /\*     maf_threshholds->a[4] = 0.08; *\/ */
/*   /\*     maf_threshholds->a[5] = 0.16; *\/ */
/*   /\*     maf_threshholds->a[6] = 0.32; *\/ */
/*   //  push_to_vdouble(maf_threshholds, 0.01); */
/*   push_to_vdouble(maf_threshholds, 0.075); */
/*   push_to_vdouble(maf_threshholds, 0.1); */
/*   push_to_vdouble(maf_threshholds, 0.125); */
/*   push_to_vdouble(maf_threshholds, 0.15); */
/*   push_to_vdouble(maf_threshholds, 1.0); */

/*   /\* for(long i=1; i<n_maf_categories; i++){ *\/ */
/*   /\*   long quantile_index = sorted_mafs->size*i/n_maf_categories; *\/ */
/*   /\*   double quantile_maf = 0.5*(sorted_mafs->a[quantile_index] + sorted_mafs->a[quantile_index+1]); *\/ */
/*   /\*   push_to_vdouble(maf_threshholds, quantile_maf); *\/ */
/*   /\*   fprintf(stderr, "quantile: %lf index: %ld  maf: %lf \n", (double)i/n_maf_categories, quantile_index, quantile_maf); *\/ */
/*   /\* } *\/ */
 
/*   long n_maf_categories = maf_threshholds->size; */
/*   FILE* fh = fopen("mafs.out", "w"); */
/*   fprintf(fh, "# n maf categories: %ld\n", n_maf_categories); */
/*   for(long i=0; i<sorted_mafs->size; i++){ */
/*     fprintf(fh, "%lf\n", sorted_mafs->a[i]); */
/*   } */
/*   fclose(fh); */

/*   //  ******   Now get arrays of indices of markers in the maf categories: */
/*   Vlong** maf_category_marker_indices = (Vlong**)malloc(n_maf_categories*sizeof(Vlong*)); */
  
/*   for(long i=0; i<maf_threshholds->size; i++){ //  */
/*     maf_category_marker_indices[i] = construct_vlong(the_genotypes_set->n_markers); */
/*   } */
/*   Vlong* missing_data_counts = the_genotypes_set->marker_missing_data_counts; */
/*   Vlong* minor_allele_counts = the_genotypes_set->marker_alt_allele_counts; */
/*   if(DO_ASSERT) assert(missing_data_counts->size == minor_allele_counts->size); */
/*   for(long i_marker=0; i_marker<minor_allele_counts->size; i_marker++){ */
/*     long ok_count = the_genotypes_set->n_accessions - missing_data_counts->a[i_marker]; */
/*     if(ok_count > 0){ */
/*       double marker_minor_allele_frequency = (double)minor_allele_counts->a[i_marker]/(double)(2*ok_count); */
/*       for(long j_maf=0; j_maf<n_maf_categories; j_maf++){ */
/* 	if(marker_minor_allele_frequency <= maf_threshholds->a[j_maf]){ */
/* 	  push_to_vlong(maf_category_marker_indices[j_maf], i_marker); */
/* 	  break; */
/* 	} */
/*       } // end loop over maf ranges */
/*     } */
/*   } */
/*   for(long i=0; i<n_maf_categories; i++){ */
    
/*     fprintf(stderr, "maf category: %ld; max maf in category: %7.4f; n markers in category: %ld \n", */
/* 	    i, maf_threshholds->a[i], maf_category_marker_indices[i]->size); */
/*   } */
/*   free_vdouble(sorted_mafs); */
/*   free_vdouble(maf_threshholds); */
/*   return maf_category_marker_indices; */
/* } */
  


/* Vagmri* maf_category_agmrs(GenotypesSet* the_gtset, Accession* acc1, Accession* acc2, */
/* 			   long n_maf_categories, Vlong** maf_cat_marker_indices){ */
/*   Vdouble* agmrs = construct_vdouble(n_maf_categories+1); */
/*   Vagmri* the_agmris = construct_vagmri(n_maf_categories+1); */
  
/*   //  fprintf(stderr, "n_maf_categories %ld  agmrs size: %ld \n", n_maf_categories, agmrs->size); */
/*   long numerator = 0; // for the calculation using all markers */
/*   long denominator = 0; // for the calculation using all markers */
/*   double nhagmr = 0; */
/*   for(long i_maf=0; i_maf<n_maf_categories; i_maf++){  */
/*     Vlong* marker_indices = maf_cat_marker_indices[i_maf]; */
/*     /\* for(long jjjj=0; jjjj<8; jjjj++){ *\/ */
/*     /\*   fprintf(stderr, "  %ld  %ld  %ld  ", i_maf, marker_indices->size, marker_indices->a[jjjj]);  *\/ */
/*     /\* }fprintf(stderr, "\n"); getchar(); *\/ */
/*     long maf_numerator = 0; */
/*     long maf_denominator = 0; */
/*     double maf_nhagmr = 0; */
/*     for(long i=0; i<marker_indices->size; i++){ */
/*       long i_marker = marker_indices->a[i]; */
/*       char a1 = acc1->genotypes->a[i_marker]; */
/*       if(a1 != MISSING_DATA_CHAR){ */
/* 	char a2 = acc2->genotypes->a[i_marker]; */
/* 	if(a2 != MISSING_DATA_CHAR){ */
/* 	  double maf = the_gtset->mafs->a[i_marker]; */
/* 	  maf_nhagmr += 4.0*maf*(1.0-maf)*(1.0 - 1.5*maf*(1.0-maf)); */
/* 	  maf_denominator++; // count all genotype pairs with neither being missing data. */
/* 	  if(a1 != a2) { // count genotype pairs with differing dosages. */
/* 	    maf_numerator += abs(a1 - a2); //  */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*     Agmri* maf_agmri = construct_agmri(maf_numerator, maf_denominator, maf_nhagmr); */
/*     push_to_vagmri(the_agmris, maf_agmri); */
/*     //   push_to_vdouble(agmrs, (maf_denominator > 0)? (double)maf_numerator/(double)maf_denominator : -1); */
/*     //   push_to_vdouble(agmrs, (maf_nhagmr > 0)? (double)maf_numerator/(double)maf_nhagmr : -1); */
/*     numerator += maf_numerator; */
/*     denominator += maf_denominator; */
/*     nhagmr += maf_nhagmr; */
  
/*     // agmrs->size = n_maf_categories; */
/*   } // end loop over maf categories */
/*   //  fprintf(stderr, "Asize of vagmri: %ld\n", the_agmris->size); */
/*   /\* for(long jjjj=0; jjjj<the_agmris->size; jjjj++){ *\/ */
/*   /\*   fprintf(stderr, " %ld ", the_agmris->a[jjjj]->n); *\/ */
/*   /\* }fprintf(stderr, "##\n"); *\/ */
/*   for(long ii = the_agmris->size-2; ii >= 0; ii--){  */
/*     the_agmris->a[ii]->n += the_agmris->a[ii+1]->n; */
/*     the_agmris->a[ii]->en += the_agmris->a[ii+1]->en; */
/*     the_agmris->a[ii]->d += the_agmris->a[ii+1]->d; */

/*   } */
/*   /\* fprintf(stderr, "###  %ld %ld  %ld %ld  %7.4f %7.4f\n", *\/ */
/*   /\* 	  numerator, the_agmris->a[0]->n, *\/ */
/*   /\* 	   denominator, the_agmris->a[0]->d, *\/ */
/*   /\* 	  nhagmr, the_agmris->a[0]->en); *\/ */

/*   //    Agmri* the_agmri = construct_agmri(numerator, denominator, nhagmr); */
/*   //   push_to_vagmri(the_agmris, the_agmri); */
/*   //push_to_vdouble(agmrs, (denominator > 0)? (double)numerator/(double)denominator : -1); */
/*   //push_to_vdouble(agmrs, (nhagmr > 0)? (double)numerator/(double)nhagmr : -1); */
/*   /\* for(long ii=0; ii<agmrs->size; ii++){ *\/ */
/*   /\*   fprintf(stderr, "%8.4f  ", agmrs->a[ii]); *\/ */
/*   /\* }fprintf(stderr, "\n"); *\/ */
/*   //  fprintf(stderr, "Bsize of vagmri: %ld\n", the_agmris->size); */
/*   return the_agmris;   // agmrs; */
/* } */
      

/* Vdouble* maf_range_agmrs(GenotypesSet* the_gtset, Accession* acc1, Accession* acc2, Vdouble* maf_threshholds){ //, long n_maf_categories){ //, Vlong** maf_cat_marker_indices){ */
/*   //  push_to_vdouble(maf_threshholds, 1.0); */
/*   Vdouble* agmrs = construct_vdouble(maf_threshholds->size); */
/*   long n_markers = acc1->genotypes->length; */
/*   if(DO_ASSERT) assert(acc2->genotypes->length == n_markers); */
/*   char* gts1 = acc1->genotypes->a; */
/*   char* gts2 = acc2->genotypes->a; */
/*   Vlong* numerators =  construct_vlong_zeroes(maf_threshholds->size); */
/*   Vlong* denominators =  construct_vlong_zeroes(maf_threshholds->size); */
/*   Vlong* marker_alt_allele_counts = the_gtset->marker_alt_allele_counts; */
/*   Vlong* marker_md_counts = the_gtset->marker_missing_data_counts; */
/*   Vlong* markers_used = construct_vlong_zeroes(maf_threshholds->size); */
/*   // fprintf(stderr, "n genotypes: %ld  %ld \n", acc1->genotypes->length, acc2->genotypes->length); */
/*   for(long i=0; i<n_markers; i++){  */
/*     char a1 = gts1[i]; */
/*     char a2 = gts2[i]; */
/*     if(DO_ASSERT) assert(a2 != '\0'  &&  a1 != '\0'); */
/*     if(a1 != MISSING_DATA_CHAR){ */
/*       if(a2 != MISSING_DATA_CHAR){ */
/* 	double marker_alt_allele_freq = 0.5*marker_alt_allele_counts->a[i]/(double)(the_gtset->n_accessions - marker_md_counts->a[i]); */
/* 	for(long j=0; j< maf_threshholds->size; j++){ */
/* 	  if(marker_alt_allele_freq <= maf_threshholds->a[j]){	     */
/* 	    denominators->a[j]++; */
/* 	    if(a1 != a2) numerators->a[j]++; */
/* 	    markers_used->a[j]++; */
/* 	    break; */
/* 	  } */
/* 	} // end loop over maf ranges */
/*       } */
/*     } */
/*     /\* for(long jjj=0; jjj<markers_used->size; jjj++){ *\/ */
/*     /\*   fprintf(stderr, "%ld  %8.5lf  ", markers_used->a[jjj], maf_threshholds->a[jjj]); *\/ */
/*     /\* }fprintf(stderr, "\n"); *\/ */
/*   } // end loop over markers */
/*   agmrs->size = numerators->size; */
/*   for(long j=0; j< agmrs->size; j++){ */
/*     agmrs->a[j] = (double)numerators->a[j]/(double)denominators->a[j]; */
/*   } */
/*   //  fprintf(stderr, "### %ld \n", agmrs->size); */
/*   return agmrs; */
/* } */



/* typedef struct{ */
/*   long n; */
/*   long d; */
/*   double en; */
/* } Agmri; */

/* typedef struct{ */
/*   long capacity; */
/*   long size; */
/*   Agmri** a; */
/* } Vagmri; */


// *****  Agmri   *********************
// Agmri* construct_agmri(long n, long d, double en);

// *****  Vagmri  ********************************
/* Vagmri* construct_vagmri(long cap); */
/* void push_to_vagmri(Vagmri* the_vagmri, Agmri* the_agmri); */
/* Agmri* pop_from_vagmri(Vagmri* the_vagmri); */
/* Agmri* get_ith_agmri_from_vagmri(Vagmri* the_vagmri, long i); */
/* void free_vagmri(Vagmri* the_vagmri); */

// *****  Agmri   *********************
/* Agmri* construct_agmri(long n, long d, double en){ */
/*   Agmri* the_agmri = (Agmri*)calloc(1, sizeof(Agmri)); */
/*   the_agmri->n = n; */
/*   the_agmri->d = d; */
/*   the_agmri->en = en; */
/*   return the_agmri; */
/* } */

// *****  Vagmri  ********************************
/* Vagmri* construct_vagmri(long cap){ */
/*   Vagmri* the_vagmri = (Vagmri*)malloc(1*sizeof(Vagmri)); */
/*   the_vagmri->capacity = cap; */
/*   the_vagmri->size = 0; */
/*   the_vagmri->a = (Agmri**)malloc(cap*sizeof(Agmri*)); */
/*   return the_vagmri; */
/* } */
/* void push_to_vagmri(Vagmri* the_vagmri, Agmri* the_agmri){ */
/*   long cap = the_vagmri->capacity; */
/*   long n = the_vagmri->size; */
/*   // if necessary, resize w realloc */
/*   if(n == cap){ */
/*     cap *= 2; */
/*     the_vagmri->a = (Agmri**)realloc(the_vagmri->a, cap*sizeof(Agmri*)); */
/*     the_vagmri->capacity = cap; */
/*   } */
/*   the_vagmri->a[n] = the_agmri; */
/*   the_vagmri->size++; */
/* } */

/* Agmri* pop_from_vagmri(Vagmri* the_vagmri){ */
/*   Agmri* the_agmri = the_vagmri->a[the_vagmri->size - 1]; */
/*   the_vagmri->size--; */
/*   return the_agmri; */
/* } */

/* Agmri* get_ith_agmri_from_vagmri(Vagmri* the_vagmri, long i){ // i=-1  ->  last element, etc. */
/*   // for(long j=0; j<the_vdouble->size; j++){ fprintf(stderr, "%f ", the_vdouble->a[j]);} fprintf(stderr, "\n"); */
/*   Agmri* x = (i >= 0)? */
/*     the_vagmri->a[i] : */
/*     the_vagmri->a[the_vagmri->size + i]; */
/*   // fprintf(stderr, "i: %ld  x: %f  vdouble->size: %ld \n", i, x, the_vdouble->size-i); */
/*   return x; */
/* } */

/* void free_vagmri(Vagmri* the_vagmri){ */
/*   if(the_vagmri == NULL) return; */
/*   for(long i=0; i< the_vagmri->size; i++){ */
/*     free(the_vagmri->a[i]); */
/*   } */
/*   free(the_vagmri->a); */
/*   free(the_vagmri); */
/* } */


/* double agmr(Accession* gtset1, Accession* gtset2){ */
/*   char* gts1 = gtset1->genotypes->a; */
/*   char* gts2 = gtset2->genotypes->a; */
/*   long usable_pair_count = 0; // = agmr_denom */
/*   long mismatches = 0; // = agmr_numerator */
/*   for(long i=0; ;i++){ */
/*     char a1 = gts1[i]; */
/*     if(a1 == '\0') break; // end of  */
/*     char a2 = gts2[i]; */
/*     if(DO_ASSERT) assert(a2 != '\0'); */
/*     if(a1 != MISSING_DATA_CHAR){ */
/*       if(a2 != MISSING_DATA_CHAR){ */
/* 	usable_pair_count++; */
/* 	if(a1 != a2) mismatches++; */
/*       } */
/*     } */
/*   } */
/*   return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1; */
/* } */


/* Vstr* read_dosages_file(char* input_filename){ */
/*   FILE* i_stream = fopen(input_filename, "r"); */
/*   if (i_stream == NULL) { */
/*     perror("fopen"); */
/*     exit(EXIT_FAILURE); */
/*   } */

/*   Vstr* the_lines = construct_vstr(10000); */
/*   char* line = NULL; */
/*   size_t len = 0; */
/*   ssize_t nread; */

/*   // *****   Read all lines; store non-comments.  ***** */
/*   long markerid_count = 0; */
/*   Vstr* marker_ids = construct_vstr(1000); */
/*   char* saveptr = NULL; */
/*   while((nread = getline(&line, &len, i_stream)) != -1){ */
/*     if(line[0] == '#') continue; */
/*     char* a_line = (char*)malloc((nread+1)*sizeof(char)); */
/*     push_to_vstr(the_lines, strcpy(a_line, line)); */
/*   } */
/*   return the_lines; */
/* } */

/* ND distance(Accession* acc1, Accession* acc2){ */
/*   char* gts1 = acc1->genotypes->a; */
/*   char* gts2 = acc2->genotypes->a; */
  
/*   long mdcount1 = 0; */
/*   long mdcount2 = 0; */
  
/*   if(DO_ASSERT) assert(acc1->genotypes->length == acc2->genotypes->length); */
/*   long denom = 0; */
/*   long numer = 0; */
/*   long h_denom = 0; */
/*   long h_numer = 0; */
/*   for(long i=0; ;i++){ */
/*     char a1 = gts1[i]; */
/*     if(a1 == '\0') break; // end of  */
/*     char a2 = gts2[i]; */
/*     if(DO_ASSERT) assert(a2 != '\0'); */
/*     if(a1 != MISSING_DATA_CHAR){ */
/*       if(a2 != MISSING_DATA_CHAR){ */
/* 	denom++; */
/* 	if(a1 != a2) */
/* 	  //  numer++; */
/* 	  numer += abs(a1 - a2); */
/*       } */
/*     } */
/*   } */
/*   ND result = {numer, denom}; */
/*   // result.n = numer; */
/*   // result.d = denom; */
/*   return result; */
/* } */


/* Vdouble* distances_random_sample(GenotypesSet* the_gtset, long n){ */
/*   Vdouble* distances = construct_vdouble(n); */
/*   Vaccession* accessions = the_gtset->accessions; */
/*   long n_ref = the_gtset->n_ref_accessions; */
/*   long n_query = accessions->size - n_ref; */
/*   for(long i=0; i<n; i++){ */
/*     while(1){ */
/*       long iq = (long)(rand()*(double)n_query/((double)RAND_MAX+1) ) + n_ref; */
/*       long ia = (long)(rand()*(double)(accessions->size)/((double)RAND_MAX+1)); */
/*       if(ia != iq){ */
/* 	ND dnd = distance(accessions->a[iq], accessions->a[ia]); */
/* 	if(dnd.d > 0){ */
/* 	  push_to_vdouble(distances, DISTANCE_NORM_FACTOR*(double)dnd.n/dnd.d); */
/* 	  //	  fprintf(stderr, "%ld %ld %8.5f\n", iq, ia, DISTANCE_NORM_FACTOR*(double)dnd.n/dnd.d); */
/* 	  break; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   return distances; */
/* } */

// *****  end of function definitions  *****
