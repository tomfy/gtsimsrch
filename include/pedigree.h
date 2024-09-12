// *****  typedefs for Pedigree, Vpedigree *****
typedef struct{ // 0: progeny, 1, 2: parents
  ND agmr01;
  ND agmr02;
  ND agmr12;
  ND par1_hgmr;
  ND par1_xhgmr;
  ND par1_R;
  ND par2_hgmr;
  ND par2_xhgmr; 
  ND par2_R;
  ND z; // (n00_1 + n22_1)/(n00_x + n22_x)
  ND d; // (n1+n2)/(n0+n1+n2);
  ND d_2; // 
  ND d_22; // both parents homozyg, delta = 2  i.e. 00_2 + 22_0
  ND d_21; // both parents homozyg, delta = 1
  ND d_11; // one parent homozyg, one heterozyg, delta = 1
  double max_dz;
  double hgmr1;
  double hgmr2;
  double xhgmr1;
  double xhgmr2;
  long all_good_count; // number of markers with genotypes for all 3 (or both in case of only one parent) genotypes known (i.e. no missing gts)
  long n_01or10_1; // alternative denom for z
}Pedigree_stats; // 

typedef struct{
  Accession* F; // female parent
  Accession* M; // male parent
  Accession* A; // accession
  Pedigree_stats* pedigree_stats; 
}Pedigree;

typedef struct{
  long capacity;
  long size;
  Pedigree** a;
}Vpedigree;

typedef struct{
  long idx;
  double hgmr;
}Idxhgmr;

// *****  function declarations  *****

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent);
//double hgmr(char* gts1, char* gts2);
Pedigree_stats* construct_pedigree_stats(void); // just initializing to 0's
Pedigree_stats* bitwise_triple_counts(Accession* par1, Accession* par2, Accession* prog);
Vpedigree* calculate_triples_for_one_accession(Accession* prog, GenotypesSet* the_genotypes_set, Viaxh* cppps, long max_candidate_parents);
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset);
//, long* d0counts, long* d1counts, long* d2counts); // , GenotypesSet* the_gtsset);
long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_ok_hgmr, double max_self_r, double max_ok_d);
void free_pedigree(const Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* read_and_store_pedigrees_3col(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset);
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset); 
Vpedigree* construct_vpedigree(long cap);
const Vlong* accessions_with_offspring(const Vpedigree* the_Vped); //, long n_accessions);

//Vlong* alternative_parents(Accession* the_acc, const GenotypesSet* const the_gtsset, double max_ok_hgmr);
Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_z, double max_ok_d);
Vpedigree* alternative_pedigrees(Accession* the_acc, const GenotypesSet* the_gtsset, Vlong* best_parent_candidate_idxs, long ub, double max_ok_d);
void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees, long max_to_print, bool verbose);
void push_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped);

/* void sort_vpedigree_by_d(Vpedigree* the_vped);
void sort_vpedigree_by_z(Vpedigree* the_vped);
int compare_pedigree_d(const void* a, const void* b);
int compare_pedigree_z(const void* a, const void* b); /* */

void sort_vpedigree_by_maxdz(Vpedigree* the_vped);
int compare_pedigree_maxdz(const void* a, const void* b);
void sort_vpedigree_by_maxh1h2z(Vpedigree* the_vped);
int compare_pedigrees_h1h2z(const void* a, const void* b);
void free_vpedigree(const Vpedigree* the_vped);

// *****  array of Idxhgmr  *****
int cmpidxhgmr(const void* v1, const void* v2);
void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array);

// *****  miscellaneous  *****
long long_min(long a, long b);
long long_max(long a, long b);

// two_longs gamete_dosage_range(long d, long ploidy);
// two_longs diploid_quick_and_dirty_triple_counts(Accession* acc1, Accession* acc2, Accession* progacc);
// four_longs q_and_d_n22x_diploid(Accession* acc1, Accession* acc2, Accession* progacc);
// ND tfc_tetraploid(char* gts1, char* gts2, char* proggt);
// ND tfc_diploid(char* gts1, char* gts2, char* proggts);
// ND TFC(char* gts1, char* gts2, char* proggts, long ploidy);
// four_longs tfca(char* gts1, char* gts2, char* proggts, long ploidy);
// four_longs triple_forbidden_counts(char* gts1, char* gts2, char* proggts, long ploidy);
// Pedigree_stats* triple_counts_x(char* gts1, char* gts2, char* proggts, long* d0counts, long* d1counts, long* d2counts);
Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts, long ploidy);

// long marker_d_counts(Pedigree* the_pedigree, long* d0counts, long* d1counts, long* d2counts);
void print_pedigree(FILE* fh, Pedigree* the_pedigree, bool verbose);
void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats, bool verbose);

//void sort_pedigree_stats_by_d(Pedigree_stats** the_pss, long size); // sort in place
//int pscmp(const void* v1, const void* v2);

void print_d_r(FILE* fh, ND nd);
// double n_over_d(ND nd);

//long check_idxid_map(Vidxid* vidxid, const Vaccession* accessions);
