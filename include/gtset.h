#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "various.h"
#include "vect.h"
//#include "pedigree.h"

#define DO_ASSERT 0
#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1
#define MAX_PLOIDY 2
#define MAX_PATTERNS 10000
#define MISSING_DATA_CHAR 126 // value which will be stored in char for missing data; bigger than ploidy likely to be.
#define INIT_VACC_CAPACITY 4000
#define STRTOL_FAIL 1  // set errno to this in str_to_long if strtol fails (in some other way besides out of range)

typedef struct{
  Vchar* id;
  long index; // the index in the accessions array of GenotypesSet
  Vchar* genotypes;
  Vlong* chunk_patterns;
  //  long md_chunk_count; // the number of chunks with missing data
  long ok_chunk_count; // the number of chunks with no missing data = n_chunks - md_chunk_count
  long missing_data_count;
  Vlong* ref_homozygs; // indices of the markers for which this acc is homozyg (ref allele)
  Vlong* alt_homozygs; //
  Vull* Abits; // Abit and Bbist encode the gts, 64 gts to each pair of unsigend long longs
  Vull* Bbits; //
  double agmr0;
  double n_exp_00_1_22_1;
  double  n_exp_00_1_22_1_self;
}Accession;

typedef struct{
  long capacity;
  long size;
  Accession** a;
}Vaccession;

typedef struct{
  //  long capacity; // needed?
  double max_marker_missing_data_fraction;
  double min_minor_allele_frequency;
  long n_accessions; // accessions stored (ref + new but bad ones not counted)
  long n_bad_accessions; // accessions rejected due to excessive missing data
  long n_ref_accessions; // ref accessions stored
  long n_markers; // redundant.
  long ploidy; //

  Vaccession* accessions;
  
  Vstr* marker_ids; // vector of marker_ids
  //  Vmarker* markers; // vector of markers
  Vlong* marker_missing_data_counts; //
  Vlong* marker_alt_allele_counts; //
  Vdouble* mafs; 
  Vlong** marker_dosage_counts; // counts of dosages for each marker
  Vlong* dosage_counts; // counts of dosages for whole
  double agmr0;
}GenotypesSet;

typedef struct{
  Vstr* accession_lines;
  long first_line;
  long last_line;
  long markerid_count;
  long ploidy;
  double max_acc_missing_data_fraction;
  
  Vlong* marker_missing_data_counts;
  Vlong* marker_alt_allele_counts;
  long n_bad_accessions;
  Vaccession* accessions;
}threaded_input_struct;

typedef struct{
  GenotypesSet* gtss;
  long first; // first acc idx
  long last;
}threaded_setAB_struct;

// *****  functions  *****
long int_power(long base, long power);
long str_to_long(char* str);
char token_to_dosage(char* token, long* ploidy);
//long determine_file_format(char* filename);

// *****  Accession  *****
Accession* construct_accession(char* id, long idx, char* genotypes, long accession_md_count);
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count);
long set_accession_chunk_patterns(Accession* the_gts, Vlong* m_indices, long n_chunks, long k, long ploidy);
char* print_accession(Accession* the_gts, FILE* ostream);

void free_accession(Accession* the_accession);
void free_accession_innards(Accession* the_accession);

double agmr0(GenotypesSet* the_gtsset);
double agmr0_qvsall(const GenotypesSet* the_gtsset, Accession* A);
double agmr0_accvsall(const GenotypesSet* the_gtsset, Accession* A);
double pair_agmr0(Accession* A, Accession* B);
void n_00_1_22_1_accvsall(const GenotypesSet* the_gtsset, Accession* A );


// *****  Vaccession  *****
Vaccession* construct_vaccession(long cap);
void push_to_vaccession(Vaccession* the_vacc, Accession* the_acc);
void shuffle_order_of_accessions(GenotypesSet* the_genotypes_set);
void set_vaccession_chunk_patterns(Vaccession* the_accessions, Vlong* m_indices, long n_chunks, long k, long ploidy);
void print_vaccession(Vaccession* the_accessions, FILE* ostream);
void check_accession_indices(Vaccession* the_accessions);
void free_vaccession(Vaccession* the_vacc);

/* // *****  Marker  ***** */
/* Marker* construct_marker(char* id, double alt_allele_freq); */

/* // *****  Vmarker  ***** */
/* Vmarker* construct_vmarker(long cap); // construct empty Vmarker with capacity cap */
/* add_marker_to_vmarker(Vmarker* the_vmarker, Marker* the_marker); */
/* free_vmarker(Vmarker* the_vmarker); */

// *****  GenotypesSet  *****
GenotypesSet* construct_empty_genotypesset(double max_marker_md_fraction, double min_min_allele_freq, long ploidy);
//GenotypesSet* read_dosages_file_and_store(char* input_filename, double delta);
//GenotypesSet* read_genotypes_file_and_store(char* input_filename);

void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set, double max_acc_missing_data_fraction, long Nthreads);
void threaded_input(FILE* in_stream, long n_lines_in_chunk, double max_acc_md_fraction, long Nthreads, Vstr* marker_ids, GenotypesSet* the_genotypes_set);
void* input_lines_1thread(void* x); // for threaded processing of input lines.
// void read_dosages_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);
void populate_marker_dosage_counts(GenotypesSet* the_gtsset);
char token_to_dosage(char* token, long* ploidy);
// void read_genotypes_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);

void check_gtsset(GenotypesSet* gtsset);
//GenotypesSet* construct_genotypesset(Vaccession* accessions, Vstr* marker_ids, Vlong* md_counts, double delta, double max_marker_md_fraction);
void print_genotypesset_stats(GenotypesSet* gtss);
void check_genotypesset(GenotypesSet* gtss);
// GenotypesSet* construct_filtered_genotypesset(const GenotypesSet* the_gtsset, double max_md_fraction);
void filter_genotypesset(GenotypesSet* the_genotypes_set, FILE* ostream);
void rectify_markers(GenotypesSet* the_gtsset);
void set_Abits_Bbits(GenotypesSet* the_genotypesset, long Nthreads); // diploid only
void* set_Abits_Bbits_1thread(void* x);
void store_homozygs(GenotypesSet* the_gtsset);
void set_agmr0s(GenotypesSet* the_gtsset);
void set_n_00_1_22_1s(GenotypesSet* the_gtsset);
Vdouble* get_minor_allele_frequencies(GenotypesSet* the_gtset);

four_longs bitwise_agmr_hgmr(Accession* acc1, Accession* acc2);
void quick_and_dirty_hgmrs(GenotypesSet* the_gtsset);
ND quick_hgmr(Accession* acc1, Accession* acc2, char ploidy_char);
four_longs quick_hgmr_R(Accession* acc1, Accession* acc2, char ploidy_char);
//ND quick_and_dirty_hgmr_a(Accession* acc1, Accession* acc2);
double hgmr(char* gts1, char* gts2);
four_longs hgmr_R(char* par_gts, char* prog_gts, char ploidy_char);
ND xhgmr(GenotypesSet* gtset, Accession* a1, Accession* a2, int quick);
ND quick_and_dirty_hgmr(Accession* acc1, Accession* acc2, char ploidy_char); // get quick 'hgmr', and then if not large get true hgmr.
two_doubles lls(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny, FILE* stream, double epsilon);
ND ghgmr_old(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny);
ND ghgmr(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny);
double ragmr(GenotypesSet* the_gtsset);
void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);

Vidxid* construct_vidxid(const GenotypesSet* the_gtsset);
Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset);
long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset);
