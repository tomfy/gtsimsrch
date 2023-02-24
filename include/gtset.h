#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "various.h"
#include "vect.h"
//#include "pedigree.h"

#define DO_ASSERT 1
#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1
#define MAX_PLOIDY 60
#define MAX_PATTERNS 10000
#define MISSING_DATA_CHAR 126 // value which will be stored in char for missing data; bigger than ploidy likely to be.
#define INIT_VACC_CAPACITY 2000
#define STRTOL_FAIL 1  // set errno to this in str_to_long if strtol fails (in some other way besides out of range)

typedef struct{
  Vchar* id;
  long index; // the index in the accessions array of GenotypesSet
  // long n_markers;
  Vchar* genotypes;
  Vlong* chunk_patterns;
  long md_chunk_count;
  long missing_data_count;
  Vlong* ref_homozygs; // indices of the markers for which this acc is homozyg (ref allele)
  // Vlong* heterozygs;
  Vlong* alt_homozygs; //
  //  Vlong* high_dosages; // indices of the markers for which this acc has dosage > ploidy/2
}Accession;

typedef struct{
  long capacity;
  long size;
  Accession** a;
}Vaccession;

typedef struct{
  long capacity; // needed?
  double max_marker_missing_data_fraction;
  double min_minor_allele_frequency;
  long n_accessions; // redundant.
  long n_bad_accessions; // accessions rejected due to excessive missing data
  long n_ref_accessions; 
  long n_markers; // redundant.
  long ploidy; //

  Vaccession* accessions;
  
  Vstr* marker_ids; // vector of marker_ids
  //  Vmarker* markers; // vector of markers
  Vlong* marker_missing_data_counts; //
  Vlong* marker_alt_allele_counts; //
  Vdouble* mafs; 
  Vlong** marker_dose_counts; // counts of dosages.
}GenotypesSet;

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

// *****  Vaccession  *****
Vaccession* construct_vaccession(long cap);
void push_to_vaccession(Vaccession* the_vacc, Accession* the_acc);
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

void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set, double max_acc_missing_data_fraction);
// void read_dosages_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);
void populate_marker_dosage_counts(GenotypesSet* the_gtsset);
char token_to_dosage(char* token, long* ploidy);
// void read_genotypes_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);

void check_gtsset(GenotypesSet* gtsset);
GenotypesSet* construct_genotypesset(Vaccession* accessions, Vstr* marker_ids, Vlong* md_counts, double delta, double max_marker_md_fraction);
void check_genotypesset(GenotypesSet* gtss);
GenotypesSet* construct_cleaned_genotypesset(const GenotypesSet* the_gtsset, double max_md_fraction);
void clean_genotypesset(GenotypesSet* the_genotypes_set);
void rectify_markers(GenotypesSet* the_gtsset);
void store_homozygs(GenotypesSet* the_gtsset);
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
