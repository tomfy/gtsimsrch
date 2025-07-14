#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdbool.h>
#include "various.h"
#include "vect.h"
// #include "pedigree.h"

#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1
#define MAX_PLOIDY 2
#define MAX_PATTERNS 10000
#define MISSING_DATA_CHAR 'X'   // 126 // value which will be stored in char for missing data; bigger than ploidy likely to be.
#define INIT_VACC_CAPACITY 4000
#define STRTOL_FAIL 1  // set errno to this in str_to_long if strtol fails (in some other way besides out of range)

typedef struct{
  Vchar* id;
  long index; // the index in the accessions array of GenotypesSet
  Vchar* genotypes;
  Vchar* phases; // 'p' or 'm' for +, -, 
  Vlong* chunk_patterns;
  //  long md_chunk_count; // the number of chunks with missing data
  long ok_chunk_count; // the number of chunks with no missing data = n_chunks - md_chunk_count
  long missing_data_count;
  Vlong* ref_homozygs; // indices of the markers for which this acc is homozyg (ref allele)
  Vlong* alt_homozygs; //
  Vull* Abits; // Abit and Bbist encode the gts, 64 gts to each pair of unsigned long longs
  Vull* Bbits; //
  double agmr0;
  double n_exp_00_1_22_1;
  double  n_exp_00_1_22_1_self;
  long n2exp0s; // considering markers with 2 in this accession, expected number of 0's in a random other accession (for xhgmr denom).
  bool has_pedigree; // true iff pedigree file gives at least one parent for this accession.
  long Fpar_idx; // index of female parent according to pedigree file (or ID_NA_INDEX if not in file)
  long Mpar_idx; // index of male parent according to pedigree file (or ID_NA_INDEX if not in file)
 
  bool search_done; // true iff search for parents has been done
  //Accession* pedFpar;
  //Accession* pedMpar;
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
  long n_raw_accessions;
  long n_bad_accessions; // accessions rejected due to excessive missing data
   long n_accessions; // accessions stored (ref + new but bad ones not counted)
  long n_ref_accessions; // ref accessions stored
  long n_markers; // redundant.
  long ploidy; //
  bool phased; // true if input is phased data (as indicated by presence of CHROMOSOME line)

  Vaccession* accessions;
  
  Vstr* marker_ids; // vector of marker_ids
  Vlong* chromosomes; // the chromosome numbers for each marker.
  Vlong* chromosome_start_indices; // the marker indices at which new chromosomes start.
  //  Vmarker* markers; // vector of markers
  Vlong* marker_missing_data_counts; //
  Vlong* marker_alt_allele_counts; //
  Vdouble* mafs; 
  Vlong** marker_dosage_counts; // counts of dosages for each marker. marker_dosage_counts->[i]->a[j] is count of dosage i for marker j
  Vlong* dosage_counts; // counts of dosages for whole
  double agmr0;

  Vchar* acc_filter_info;
  Vchar* marker_filter_info;
  //double d_scale_factor;
  double mean_hgmr;
  double mean_R;
  double mean_d;
  double mean_z;
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
two_chars token_to_dosage(char* token, long* ploidy);
//long determine_file_format(char* filename);

// *****  Accession  *****
Accession* construct_accession(char* id, long idx, char* genotypes, char* phases, long accession_md_count);
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count);
long set_accession_chunk_patterns(Accession* the_gts, Vlong* m_indices, long n_chunks, long k, long ploidy);
char* print_accession(Accession* the_acc, FILE* ostream);
void free_accession(Accession* the_accession);

// *****  Vaccession  *****
Vaccession* construct_vaccession(long cap);
void push_to_vaccession(Vaccession* the_vacc, Accession* the_acc);
void set_vaccession_chunk_patterns(Vaccession* the_accessions, Vlong* m_indices, long n_chunks, long k, long ploidy);
void print_vaccession(Vaccession* the_accessions, FILE* ostream);
void check_accession_indices(Vaccession* the_accessions);
void free_vaccession(Vaccession* the_vacc);

// *****  GenotypesSet  *****
GenotypesSet* construct_empty_genotypesset(double max_marker_md_fraction, double min_min_allele_freq, long ploidy);
void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set, double max_acc_missing_data_fraction, long Nthreads);
void threaded_input(FILE* in_stream, long n_lines_in_chunk, double max_acc_md_fraction, long Nthreads, Vstr* marker_ids, GenotypesSet* the_genotypes_set);
void* input_lines_1thread(void* x); // for threaded processing of input lines.
void populate_marker_dosage_counts(GenotypesSet* the_gtsset);

void print_genotypesset_stats(GenotypesSet* gtss);
void check_genotypesset(GenotypesSet* gtss);
void filter_genotypesset(GenotypesSet* the_genotypes_set);
void rectify_markers(GenotypesSet* the_gtsset);
void set_Abits_Bbits(GenotypesSet* the_genotypesset, long Nthreads); // diploid only
void* set_Abits_Bbits_1thread(void* x);
void store_homozygs(GenotypesSet* the_gtsset);
void set_chromosome_start_indices(GenotypesSet* the_gtsset);
// void set_agmr0s(GenotypesSet* the_gtsset);
// void set_n2exp0s(GenotypesSet* gtsset, long i);
// Vdouble* get_minor_allele_frequencies(GenotypesSet* the_gtset);

four_longs bitwise_agmr_hgmr(Accession* acc1, Accession* acc2);
ND bitwise_hgmr(Accession* acc1, Accession* acc2);
// ND bitwise_R(Accession* parent, Accession* offspring);
// Viaxh** calculate_hgmrs_old(GenotypesSet* the_genotypes_set, long max_candidate_parents, double max_hgmr);
Viaxh** calculate_hgmrs(GenotypesSet* the_genotypes_set, const Vlong* cand_parent_idxs, long max_candidate_parents, double max_hgmr);
// void quick_and_dirty_hgmrs(GenotypesSet* the_gtsset);
// ND quick_hgmr(Accession* acc1, Accession* acc2, char ploidy_char);
// four_longs quick_hgmr_R(Accession* acc1, Accession* acc2, char ploidy_char);
//ND quick_and_dirty_hgmr_a(Accession* acc1, Accession* acc2);
double hgmr(char* gts1, char* gts2);
four_longs hgmr_R(char* par_gts, char* prog_gts, char ploidy_char);
// ND xhgmr(const GenotypesSet* gtset, Accession* a1, Accession* a2, bool quick);
// void calculate_xhgmrs(GenotypesSet* the_genotypes_set, Viaxh** progeny_cplds, bool quick_xhgmr, double max_xhgmr);
// ND quick_and_dirty_hgmr(Accession* acc1, Accession* acc2, char ploidy_char); // get quick 'hgmr', and then if not large get true hgmr.
// two_doubles lls(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny, FILE* stream, double epsilon);
//ND ghgmr_old(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny);
//ND ghgmr(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny);
//ND psr(Accession* acc1, Accession* acc2, Vlong* chroms); // phase mismatch rate
void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);

Vidxid* construct_vidxid(const Vaccession* accessions);
Vidxid* construct_sorted_vidxid(const Vaccession* accessions);
long check_idxid_map(Vidxid* vidxid, const Vaccession* accessions);

ND phase_switches_one_chrom(Vchar* p1s, Vchar* p2s, Vlong* chroms, long* start);
ND phase_switches(Accession* acc1, Accession* acc2, Vlong* chroms);

three_longs heterozyg_ratios(Accession* acc1, Accession* acc2);
void read_gts_line_add_accession_to_gtset(GenotypesSet* the_genotypes_set, char* acc_id, long markerid_count, char* saveptr, double max_acc_missing_data_fraction);
// ##### unused #####

double agmr0(GenotypesSet* the_gtsset);
double agmr0_qvsall(const GenotypesSet* the_gtsset, Accession* A);
double agmr0_accvsall(const GenotypesSet* the_gtsset, Accession* A);
double pair_agmr0(Accession* A, Accession* B);
void n_00_1_22_1_accvsall(const GenotypesSet* the_gtsset, Accession* A );

two_doubles logPABlogPBA(GenotypesSet* the_gtsset, Accession* A, Accession* B);
