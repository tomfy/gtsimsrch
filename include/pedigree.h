// *****  typedefs for Pedigree, Vpedigree *****
typedef struct{
  long Xa; 
  long Xb;
  long Nhet;
}Xcounts_2;

typedef struct{ // same as Xcounts_2, except Xmin (Xmax)
  // is the min (max) the 2 
  long Xmin; 
  long Xmax;
  long Nhet;
}Xcounts_2mmn;

typedef struct{
  Xcounts_2mmn XFA;
  Xcounts_2mmn XMA;
  long XFmin_3;
  long XFmax_3;
  long XMmin_3;
  long XMmax_3;
}Xcounts_3;

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
  ND d; // (n1+n2)/(n0+n1+n2);
  ND z; // (n00_1 + n22_1)/(n00_x + n22_x)

  // ND xz;
  // ND d_2; // 
  // ND d_22; // both parents homozyg, delta = 2  i.e. 00_2 + 22_0
  // ND d_21; // both parents homozyg, delta = 1
  // ND d_11; // one parent homozyg, one heterozyg, delta = 1

  double hgmr1_n;
  double R1_n;
  double hgmr2_n;
  double R2_n;
  double d_n;
  double z_n;
  double ftc_n;

  double hgmr1;
  double hgmr2;
  double xhgmr1;
  double xhgmr2;

  Xcounts_3 X3;

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

/* typedef struct{ */
/*   long idx; */
/*   double hgmr; */
/* }Idxhgmr; */

typedef struct{
  long XFmin;
  long XFmax;
  long NFhet;
  long XMmin;
  long XMmax;
  long NMhet;
  long XFmin_triple;
  long XFmax_triple;
  long XMmin_triple;
  long XMmax_triple;
}Xover_info;




// *****  function declarations  *****

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent);
//double hgmr(char* gts1, char* gts2);
Pedigree_stats* construct_pedigree_stats(void); // just initializing to 0's
Pedigree_stats* bitwise_triple_counts(Accession* par1, Accession* par2, Accession* prog);
Vpedigree* calculate_triples_for_one_accession(Accession* prog, const GenotypesSet* the_genotypes_set, Viaxh* cppps, long max_candidate_parents);
double calculate_xFTR(Pedigree* the_pedigree, const GenotypesSet* the_gtsset);
double calculate_xxFTR(Pedigree* the_pedigree, const GenotypesSet* the_gtsset, double alpha);
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, const GenotypesSet* the_gtsset);

// long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_self_r, double max_ok_d);
bool d_ok(Pedigree_stats* p, double max_ok_d); // 1: d looks good; 0: d too large.

void free_pedigree(const Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* construct_vpedigree(long cap);
Vpedigree* read_and_store_pedigrees_3col(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset);
const Vlong* accessions_with_offspring(const Vpedigree* the_Vped, long n_accessions);
// Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_d);
// void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees, long max_to_print, bool verbose);
void push_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped);
void sort_vpedigree_by_d(Vpedigree* the_vped);
int compare_pedigree_d(const void* a, const void* b);

void free_vpedigree(const Vpedigree* the_vped);

// *****  array of Idxhgmr  *****
//int cmpidxhgmr(const void* v1, const void* v2);
//void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array);

//Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts, long ploidy);

Xcounts_3 count_crossovers(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring); 
Xcounts_2mmn count_crossovers_one_parent(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring);
Xcounts_2 count_crossovers_one_chromosome(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring, long first, long last);
Xcounts_3 count_crossovers_two_parents(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring);
two_longs get_1marker_phases_wrt_1parent(char p_phase, char o_gt, char o_phase);

void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats, bool verbose);
void print_pedigree_normalized(FILE* fh, Pedigree* the_pedigree);
void print_normalized_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats);
void print_double_nan_as_hyphen(FILE* fh, double x);
