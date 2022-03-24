// *****  typedefs for Pedigree, Vpedigree *****
typedef struct{
  long n; // numerator
  long d; // denominator
}ND; // numerator and denominator

typedef struct{
  ND agmr12;
  ND par1_hgmr;
  //  ND par1_r;
  ND par1_R;
  ND par2_hgmr;
  //  ND par2_r;
  ND par2_R;
  //  ND d1;
  //  ND d2;
  ND d;
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
double hgmr(char* gts1, char* gts2);
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree);
//, long* d0counts, long* d1counts, long* d2counts); // , GenotypesSet* the_gtsset);
long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_ok_hgmr, double max_self_r, double max_ok_d);
void free_pedigree(const Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset); 
Vpedigree* construct_vpedigree(long cap);
const Vlong* accessions_with_offspring(const Vpedigree* the_Vped, long n_accessions);
const Vaccession* accessions_with_offspring_x(const Vpedigree* the_vped, const GenotypesSet* the_gtsset);
Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_d);
void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees);
void add_pedigree_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped);
void free_vpedigree(const Vpedigree* the_vped);

// *****  array of Idxhgmr  *****
int cmpidxhgmr(const void* v1, const void* v2);
void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array);

// *****  miscellaneous  *****
long long_min(long a, long b);
long long_max(long a, long b);


Pedigree_stats* triple_counts_x(char* gts1, char* gts2, char* proggts,
			   long* d0counts, long* d1counts, long* d2counts);
Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts);

long marker_d_counts(Pedigree* the_pedigree,
		     // char* gts1, char* gts2, char* proggts,
		     long* d0counts, long* d1counts, long* d2counts);
void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats);
// void print_pedigree_stats_x(FILE* fh, Pedigree_stats* the_pedigree_stats);
double get_agmr12(Pedigree_stats* p);
double get_hgmr1(Pedigree_stats* p);
//double get_r1(Pedigree_stats* p);
double get_R1(Pedigree_stats* p);
double get_hgmr2(Pedigree_stats* p);
//double get_r2(Pedigree_stats* p);
double get_R2(Pedigree_stats* p);
//double get_d1(Pedigree_stats* p);
//double get_d2(Pedigree_stats* p);
double get_d(Pedigree_stats* p);

long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset);
