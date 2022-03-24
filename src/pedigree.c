#include <stdlib.h>
#include <stdio.h>
#include "gtset.h"
// #include "vect.h"
#include "pedigree.h"

#define PEDIGREE_FIELDS 7 // number of whitespace separated fields in pedigree file, with ids in last 3 fields.

// extern int do_checks_flag; // option -c sets this to 1 to do some checks.

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent){ //IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid){
  Pedigree* the_pedigree = (Pedigree*)malloc(sizeof(Pedigree));
  the_pedigree->F = Fparent;
  the_pedigree->M = Mparent;
  the_pedigree->A = Acc;
  the_pedigree->pedigree_stats = NULL;
  return the_pedigree;
}

double hgmr(char* gts1, char* gts2){
  char c1, c2;
  long n_numer = 0;
  long n_denom = 0;
  long i=0;
  while((c1 = gts1[i]) != '\0'){
    if((c1 == '0') || (c1 == '2')){
      c2 = gts2[i];
      if((c2 == '0') || (c2 == '2')){
	n_denom++;
	if(c1 != c2) n_numer++;
      }
    }
    i++;
  }
  return (n_denom > 0)? (double)n_numer/(double)n_denom : 2.0;  
}


long marker_d_counts(Pedigree* the_pedigree,
		     // char* gts1, char* gts2, char* proggts,
		     long* d0counts, long* d1counts, long* d2counts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){
  char* gts1 = the_pedigree->F->genotypes->a;
  char* gts2 = the_pedigree->M->genotypes->a;
  char* proggts = the_pedigree->A->genotypes->a;
 
  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0;
  
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0;
  
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
  long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0;

  long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0;
  long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0;
  long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0;
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){
    if(c3 != '3'){
      c1 = gts1[i];
      c2 = gts2[i];
      if(c1 == '0'){ // 0xy
	if(c2 == '0'){ // 00y
	  if(c3 == '0'){ 
	    n_00_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_00_1++; d1counts[i]++;
	  }else if(c3 == '2'){
	    n_00_2++; d2counts[i]++;
	  }
	}else if(c2 == '1'){ // 01y
	  if(c3 == '0'){ 
	    n_01_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_01_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_01_2++; d1counts[i]++;
	  }
	}else if(c2 == '2'){ // 02y
	  if(c3 == '0'){
	    n_02_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_02_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_02_2++; d1counts[i]++;
	  }
	}
      }else if(c1 == '1'){ // 1xy
	if(c2 == '0'){ // 10y
	  if(c3 == '0'){
	    n_10_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_10_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_10_2++; d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 11y
	  if(c3 == '0'){
	    n_11_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_11_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_11_2++; d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 12y
	  if(c3 == '0'){
	    n_12_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_12_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_12_2++; d0counts[i]++;
	  }
	}
      }else if(c1 == '2'){ // 2xy
	if(c2 == '0'){ // 20y
	  if(c3 == '0'){
	    n_20_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_20_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_20_2++; d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 21y
	  if(c3 == '0'){
	    n_21_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_21_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_21_2++; d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 22y
	  if(c3 == '0'){
	    n_22_0++; d2counts[i]++;
	  }else if(c3 == '1'){
	    n_22_1++; d1counts[i]++;
	  }else if(c3 == '2'){
	    n_22_2++; d0counts[i]++;
	  }
	}
      }
    }
    i++;
  }
  return i;
}


Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts){ // 
  //  long* d0counts, long* d1counts, long* d2counts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){

  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0;
  
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0;
  
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
  long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0;

  long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0;
  long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0;
  long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0;
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){
    if(c3 != '3'){
      c1 = gts1[i];
      c2 = gts2[i];
      if(c1 == '0'){ // 0xy
	if(c2 == '0'){ // 00y
	  if(c3 == '0'){ 
	    n_00_0++; //d0counts[i]++;
	  }else if(c3 == '1'){
	    n_00_1++; //d1counts[i]++;
	  }else if(c3 == '2'){
	    n_00_2++; //d2counts[i]++;
	  }
	}else if(c2 == '1'){ // 01y
	  if(c3 == '0'){ 
	    n_01_0++; //d0counts[i]++;
	  }else if(c3 == '1'){
	    n_01_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_01_2++; //d1counts[i]++;
	  }
	}else if(c2 == '2'){ // 02y
	  if(c3 == '0'){
	    n_02_0++; //d1counts[i]++;
	  }else if(c3 == '1'){
	    n_02_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_02_2++; //d1counts[i]++;
	  }
	}
	else if(c2 == '3'){ // 03y
	  if(c3 == '0'){
	    n_03_0++;
	  }else if(c3 == '1'){
	    n_03_1++;
	  }else if(c3 == '2'){
	    n_03_2++;
	  }
	}
      }else if(c1 == '1'){ // 1xy
	if(c2 == '0'){ // 10y
	  if(c3 == '0'){
	    n_10_0++; //d0counts[i]++;
	  }else if(c3 == '1'){
	    n_10_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_10_2++; //d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 11y
	  if(c3 == '0'){
	    n_11_0++; //d0counts[i]++;
	  }else if(c3 == '1'){
	    n_11_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_11_2++; //d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 12y
	  if(c3 == '0'){
	    n_12_0++; //d1counts[i]++;
	  }else if(c3 == '1'){
	    n_12_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_12_2++; //d0counts[i]++;
	  }
	}else if(c2 == '3'){ // 13y
	  if(c3 == '0'){
	    n_13_0++;
	  }else if(c3 == '1'){
	    n_13_1++;
	  }else if(c3 == '2'){
	    n_13_2++;
	  }
	}
      }else if(c1 == '2'){ // 2xy
	if(c2 == '0'){ // 20y
	  if(c3 == '0'){
	    n_20_0++; //d1counts[i]++;
	  }else if(c3 == '1'){
	    n_20_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_20_2++; //d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 21y
	  if(c3 == '0'){
	    n_21_0++; //d1counts[i]++;
	  }else if(c3 == '1'){
	    n_21_1++; //d0counts[i]++;
	  }else if(c3 == '2'){
	    n_21_2++; //d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 22y
	  if(c3 == '0'){
	    n_22_0++; //d2counts[i]++;
	  }else if(c3 == '1'){
	    n_22_1++; //d1counts[i]++;
	  }else if(c3 == '2'){
	    n_22_2++; //d0counts[i]++;
	  }
	}else if(c2 == '3'){ // 23y
	  if(c3 == '0'){
	    n_23_0++;
	  }else if(c3 == '1'){
	    n_23_1++;
	  }else if(c3 == '2'){
	    n_23_2++;
	  }
	}
      }else if(c1 == '3'){ // 3xy
	if(c2 == '0'){ // 30y
	  if(c3 == '0'){
	    n_30_0++;
	  }else if(c3 == '1'){
	    n_30_1++;
	  }else if(c3 == '2'){
	    n_30_2++;
	  }
	}else if(c2 == '1'){ // 31y
	  if(c3 == '0'){
	    n_31_0++;
	  }else if(c3 == '1'){
	    n_31_1++;
	  }else if(c3 == '2'){
	    n_31_2++;
	  }
	}else if(c2 == '2'){ // 32y
	  if(c3 == '0'){
	    n_32_0++;
	  }else if(c3 == '1'){
	    n_32_1++;
	  }else if(c3 == '2'){
	    n_32_2++;
	  }
	}else{ // 33y
	  // do nothing
	}
      }
    }
    i++;
  }
  long n_0 = // 15 triples consistent with true parents-offspring relationship,
    // with no genotyping errors
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 + n_11_0 + n_11_1 + n_11_2 + n_12_1 + n_12_2 +
    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 = // 10 triple requiring one (0<->1 or 1<->2) error for consistency 
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 of 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0;
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1;
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2;
  
  //  long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0;
  //  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1;
  //  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2;
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0;
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1;
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2;

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0;
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1;
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2;
  
  //  long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0;
  //  long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1;
  //  long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2;
  
  long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0;
  long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1;
  long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2;

  long hgmr1_numer = n_0x_2 + n_2x_0;
  long hgmr2_numer = n_x0_2 + n_x2_0;
  long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer;
  long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer;

  long r0x1or2_numer = n_0x_1 + n_0x_2 + n_2x_1 + n_2x_0;
  long r0x1or2_denom = r0x1or2_numer + n_0x_0 + n_2x_2;
  long rx01or2_numer = n_x0_1 + n_x0_2 + n_x2_1 + n_x2_0;
  long rx01or2_denom = rx01or2_numer + n_x0_0 + n_x2_2;

  long agmr12_numer =
    n_01_0 + n_01_1 + n_01_2 + //n_01_3 +
    n_02_0 + n_02_1 + n_02_2 + //n_02_3 +
    n_10_0 + n_10_1 + n_10_2 + //n_10_3 +
    n_12_0 + n_12_1 + n_12_2 + //n_12_3 +
    n_20_0 + n_20_1 + n_20_2 + //n_20_3 +
    n_21_0 + n_21_1 + n_21_2; //n_21_3;
  
  long agmr12_denom = agmr12_numer +
    n_00_0 + n_00_1 + n_00_2 + //n_00_3 +
    n_11_0 + n_11_1 + n_11_2 + //n_11_3 +
    n_22_0 + n_22_1 + n_22_2; //n_22_3;

  Pedigree_stats* pedigree_stats = (Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  ND agmr12_nd = {agmr12_numer, agmr12_denom};
  pedigree_stats->agmr12 = agmr12_nd;
  
  ND hgmr1_nd = {hgmr1_numer, hgmr1_denom};
  pedigree_stats->par1_hgmr = hgmr1_nd;
  ND R1_nd = {r0x1or2_numer, r0x1or2_denom};
  pedigree_stats->par1_R = R1_nd;
  
  ND hgmr2_nd = {hgmr2_numer, hgmr2_denom};
  pedigree_stats->par2_hgmr = hgmr2_nd; 
  ND R2_nd = {rx01or2_numer, rx01or2_denom};
  pedigree_stats->par2_R = R2_nd;
  
  ND d_nd = {n_1 + n_2, n_0 + n_1 + n_2};
  pedigree_stats->d = d_nd;

  return pedigree_stats;
}

Pedigree_stats* triple_counts_x(char* gts1, char* gts2, char* proggts, // 
    long* d0counts, long* d1counts, long* d2counts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){

  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0;
  
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0;
  
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
  long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0;

  long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0;
  long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0;
  long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0;
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){
    if(c3 != '3'){
      c1 = gts1[i];
      c2 = gts2[i];
      if(c1 == '0'){ // 0xy
	if(c2 == '0'){ // 00y
	  if(c3 == '0'){ 
	    n_00_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_00_1++; d1counts[i]++;
	  }else if(c3 == '2'){
	    n_00_2++; d2counts[i]++;
	  }
	}else if(c2 == '1'){ // 01y
	  if(c3 == '0'){ 
	    n_01_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_01_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_01_2++; d1counts[i]++;
	  }
	}else if(c2 == '2'){ // 02y
	  if(c3 == '0'){
	    n_02_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_02_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_02_2++; d1counts[i]++;
	  }
	}
	else if(c2 == '3'){ // 03y
	  if(c3 == '0'){
	    n_03_0++;
	  }else if(c3 == '1'){
	    n_03_1++;
	  }else if(c3 == '2'){
	    n_03_2++;
	  }
	}
      }else if(c1 == '1'){ // 1xy
	if(c2 == '0'){ // 10y
	  if(c3 == '0'){
	    n_10_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_10_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_10_2++; d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 11y
	  if(c3 == '0'){
	    n_11_0++; d0counts[i]++;
	  }else if(c3 == '1'){
	    n_11_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_11_2++; d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 12y
	  if(c3 == '0'){
	    n_12_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_12_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_12_2++; d0counts[i]++;
	  }
	}else if(c2 == '3'){ // 13y
	  if(c3 == '0'){
	    n_13_0++;
	  }else if(c3 == '1'){
	    n_13_1++;
	  }else if(c3 == '2'){
	    n_13_2++;
	  }
	}
      }else if(c1 == '2'){ // 2xy
	if(c2 == '0'){ // 20y
	  if(c3 == '0'){
	    n_20_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_20_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_20_2++; d1counts[i]++;
	  }
	}else if(c2 == '1'){ // 21y
	  if(c3 == '0'){
	    n_21_0++; d1counts[i]++;
	  }else if(c3 == '1'){
	    n_21_1++; d0counts[i]++;
	  }else if(c3 == '2'){
	    n_21_2++; d0counts[i]++;
	  }
	}else if(c2 == '2'){ // 22y
	  if(c3 == '0'){
	    n_22_0++; d2counts[i]++;
	  }else if(c3 == '1'){
	    n_22_1++; d1counts[i]++;
	  }else if(c3 == '2'){
	    n_22_2++; d0counts[i]++;
	  }
	}else if(c2 == '3'){ // 23y
	  if(c3 == '0'){
	    n_23_0++;
	  }else if(c3 == '1'){
	    n_23_1++;
	  }else if(c3 == '2'){
	    n_23_2++;
	  }
	}
      }else if(c1 == '3'){ // 3xy
	if(c2 == '0'){ // 30y
	  if(c3 == '0'){
	    n_30_0++;
	  }else if(c3 == '1'){
	    n_30_1++;
	  }else if(c3 == '2'){
	    n_30_2++;
	  }
	}else if(c2 == '1'){ // 31y
	  if(c3 == '0'){
	    n_31_0++;
	  }else if(c3 == '1'){
	    n_31_1++;
	  }else if(c3 == '2'){
	    n_31_2++;
	  }
	}else if(c2 == '2'){ // 32y
	  if(c3 == '0'){
	    n_32_0++;
	  }else if(c3 == '1'){
	    n_32_1++;
	  }else if(c3 == '2'){
	    n_32_2++;
	  }
	}else{ // 33y
	  // do nothing
	}
      }
    }
    i++;
  }
  long n_0 =
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 + n_11_0 + n_11_1 + n_11_2 + n_12_1 + n_12_2 +

    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 =
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 of 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.

  //  double d1 = (n_0 > 0)? (double)n_1/(double)(n_0 + n_1) : 2;
  //  double d2 = (n_0 > 0)? (double)n_2/(double)(n_0 + n_2) : 2;
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0;
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1;
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2;
  
  long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0;
  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1;
  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2;
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0;
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1;
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2;

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0;
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1;
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2;
  
  long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0;
  long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1;
  long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2;
  
  long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0;
  long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1;
  long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2;

  long hgmr1_numer = n_0x_2 + n_2x_0;
  long hgmr2_numer = n_x0_2 + n_x2_0;
  long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer;
  long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer;
  // double hgmr1 = (hgmr1_denom > 0)? (double)hgmr1_numer/(double)hgmr1_denom : 2;
  // double hgmr2 = (hgmr2_denom > 0)? (double)hgmr2_numer/(double)hgmr2_denom : 2;

  //  long r0x1_numer = n_0x_1 + n_2x_1;
  //  long r0x1_denom = r0x1_numer + n_0x_0 + n_2x_2;
  //  long rx01_numer = n_x0_1 + n_x2_1;
  //  long rx01_denom = rx01_numer + n_x0_0 + n_x2_2;
  long r0x1or2_numer = n_0x_1 + n_0x_2 + n_2x_1 + n_2x_0;
  long r0x1or2_denom = r0x1or2_numer + n_0x_0 + n_2x_2;
  long rx01or2_numer = n_x0_1 + n_x0_2 + n_x2_1 + n_x2_0;
  long rx01or2_denom = rx01or2_numer + n_x0_0 + n_x2_2;
  // double r0x1 = (r0x1_denom > 0)? (double)r0x1_numer/(double)r0x1_denom : 2;
  // double rx01 = (rx01_denom > 0)? (double)rx01_numer/(double)rx01_denom : 2;

  long agmr12_numer =
    n_01_0 + n_01_1 + n_01_2 + //n_01_3 +
    n_02_0 + n_02_1 + n_02_2 + //n_02_3 +
    n_10_0 + n_10_1 + n_10_2 + //n_10_3 +
    n_12_0 + n_12_1 + n_12_2 + //n_12_3 +
    n_20_0 + n_20_1 + n_20_2 + //n_20_3 +
    n_21_0 + n_21_1 + n_21_2; //n_21_3;
  
  long agmr12_denom = agmr12_numer +
    n_00_0 + n_00_1 + n_00_2 + //n_00_3 +
    n_11_0 + n_11_1 + n_11_2 + //n_11_3 +
    n_22_0 + n_22_1 + n_22_2; //n_22_3;
    
  /* long agmr1_numer = hgmr1_numer + n_0x_1 + n_2x_1 + n_1x_0 + n_1x_2; */
  /* long agmr2_numer = hgmr2_numer + n_x0_1 + n_x2_1 + n_x1_0 + n_x1_2; */
  /* long agmr1_denom = agmr1_numer + n_0x_0 + n_1x_1 + n_2x_2; */
  /* long agmr2_denom = agmr2_numer + n_x0_0 + n_x1_1 + n_x2_2; */
  /* double agmr1 = (agmr1_denom > 0)? (double)agmr1_numer/(double)agmr1_denom : 2; */
  /* double agmr2 = (agmr2_denom > 0)? (double)agmr2_numer/(double)agmr2_denom : 2; */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr1,  hgmr1, r0x1); */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr2,  hgmr2, rx01); */
  /* fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); */

  Pedigree_stats* pedigree_stats = (Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  ND agmr12_nd = {agmr12_numer, agmr12_denom};
  pedigree_stats->agmr12 = agmr12_nd;
  ND hgmr1_nd = {hgmr1_numer, hgmr1_denom};
  pedigree_stats->par1_hgmr = hgmr1_nd;
  //  ND r1_nd = {r0x1_numer, r0x1_denom};
  //  pedigree_stats->par1_r = r1_nd;
  ND R1_nd = {r0x1or2_numer, r0x1or2_denom};
  pedigree_stats->par1_R = R1_nd;
  
  ND hgmr2_nd = {hgmr2_numer, hgmr2_denom};
  pedigree_stats->par2_hgmr = hgmr2_nd; 
  //  ND r2_nd = {rx01_numer, rx01_denom};
  //  pedigree_stats->par2_r = r2_nd;
 ND R2_nd = {rx01or2_numer, rx01or2_denom};
  pedigree_stats->par2_R = R2_nd;
  
  /* ND d1_nd = {n_1, n_0 + n_1}; */
  /* pedigree_stats->d1 = d1_nd; */
  /* ND d2_nd = {n_2, n_0 + n_2}; */
  /* pedigree_stats->d2 = d2_nd; */

  ND d_nd = {n_1 + n_2, n_0 + n_1 + n_2};
  pedigree_stats->d = d_nd;
  
  /* if(hgmr1 <= max_hgmr  &&  hgmr2 <= max_hgmr  &&  d1 <= max_d1){ */
  /*   fprintf(stderr, "  %-30s %-30s  ", id1, id2); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr1, r0x1); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr2, rx01); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); */
  /*   fprintf(stderr, "\n"); */
  /* } */

  return pedigree_stats;
}



/* void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset){ */

/*   char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index]; */
/*   char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index]; */
/*   char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index]; */
    
/*   Ahr fp_ahr, mp_ahr, fm_ahr; */
/*   agmr_hgmr_r(fempar_gts, acc_gts, &(the_pedigree->fp_ahr)); */
/*   agmr_hgmr_r(malpar_gts, acc_gts, &(the_pedigree->mp_ahr)); */
/*   agmr_hgmr_r(fempar_gts, malpar_gts, &(the_pedigree->fm_ahr)); */
/* } */

Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree){ //, long* d0counts, long* d1counts, long* d2counts){ //, GenotypesSet* the_gtsset){
  return triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a );			 
}


const Vlong* accessions_with_offspring(const Vpedigree* the_vped, long n_accessions){
  Vlong* offspring_counts = construct_vlong_zeroes(n_accessions);
  for(long i=0; i<the_vped->size; i++){
    const Pedigree* the_ped = the_vped->a[i];
    // fprintf(stderr, "%ld %ld \n", the_ped->F->index, the_ped->Fparent->index);
    long fp_index = the_ped->F->index; // Fparent->index;
    offspring_counts->a[fp_index]++;
    long mp_index = the_ped->M->index; // parent->index;
    offspring_counts->a[mp_index]++; 
  }
  Vlong* accidxs_with_offspring = construct_vlong(100);
  //  Vaccession* accessions = construct_vaccession(100);
  for(long i=0; i<offspring_counts->size; i++){
    if(offspring_counts->a[i] > 0){
      add_long_to_vlong(accidxs_with_offspring, i);
      //    add_accession_to_vaccession(accessions, the_gtsset->accessions->a[i];
    }
  }
  free_vlong(offspring_counts);
  return accidxs_with_offspring;
}

const Vaccession* accessions_with_offspring_x(const Vpedigree* the_vped, const GenotypesSet* the_gtsset){
  Vlong* offspring_counts = construct_vlong_zeroes(the_gtsset->n_accessions);
  for(long i=0; i<the_vped->size; i++){
    const Pedigree* the_ped = the_vped->a[i];
    // fprintf(stderr, "%ld %ld \n", the_ped->F->index, the_ped->Fparent->index);
    long fp_index = the_ped->F->index; // Fparent->index;
    offspring_counts->a[fp_index]++;
    long mp_index = the_ped->M->index; // parent->index;
    offspring_counts->a[mp_index]++; 
  }
  //  Vlong* accidxs_with_offspring = construct_vlong(100);
  Vaccession* accessions_with_offspring = construct_vaccession(100);
  for(long i=0; i<offspring_counts->size; i++){
    if(offspring_counts->a[i] > 0){
      // add_long_to_vlong(accidxs_with_offspring, i);
      add_accession_to_vaccession(accessions_with_offspring, the_gtsset->accessions->a[i]);
    }
  }
  free_vlong(offspring_counts);
  return accessions_with_offspring;
}

/* void print_pedigree_stats_x(FILE* fh, Pedigree_stats* the_pedigree_stats){ */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->agmr12.d > 0)? (double)the_pedigree_stats->agmr12.n/(double)the_pedigree_stats->agmr12.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->par1_hgmr.d > 0)? (double)the_pedigree_stats->par1_hgmr.n/(double)the_pedigree_stats->par1_hgmr.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->par1_r.d > 0)? (double)the_pedigree_stats->par1_r.n/(double)the_pedigree_stats->par1_r.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->par2_hgmr.d > 0)? (double)the_pedigree_stats->par2_hgmr.n/(double)the_pedigree_stats->par2_hgmr.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->par2_r.d > 0)? (double)the_pedigree_stats->par2_r.n/(double)the_pedigree_stats->par2_r.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->d1.d > 0)? (double)the_pedigree_stats->d1.n/(double)the_pedigree_stats->d1.d : 2); */
/*   fprintf(fh, "%6.5lf ", (the_pedigree_stats->d2.d > 0)? (double)the_pedigree_stats->d2.n/(double)the_pedigree_stats->d2.d : 2); */
/* } */

/* void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){ */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->agmr12.d, (the_pedigree_stats->agmr12.d > 0)? */
/* 	  (double)(the_pedigree_stats->agmr12.n+1)/(double)(the_pedigree_stats->agmr12.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->par1_hgmr.d, (the_pedigree_stats->par1_hgmr.d > 0)? */
/* 	  (double)(the_pedigree_stats->par1_hgmr.n+1)/(double)(the_pedigree_stats->par1_hgmr.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->par1_r.d, (the_pedigree_stats->par1_r.d > 0)? */
/* 	  (double)(the_pedigree_stats->par1_r.n+1)/(double)(the_pedigree_stats->par1_r.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->par2_hgmr.d, (the_pedigree_stats->par2_hgmr.d > 0)? */
/* 	  (double)(the_pedigree_stats->par2_hgmr.n+1)/(double)(the_pedigree_stats->par2_hgmr.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->par2_r.d, (the_pedigree_stats->par2_r.d > 0)? */
/* 	  (double)(the_pedigree_stats->par2_r.n+1)/(double)(the_pedigree_stats->par2_r.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->d1.d, (the_pedigree_stats->d1.d > 0)? */
/* 	  (double)(the_pedigree_stats->d1.n+1)/(double)(the_pedigree_stats->d1.d+1) : 2); */
/*   fprintf(fh, "%4ld %6.5lf ", the_pedigree_stats->d2.d, (the_pedigree_stats->d2.d > 0)? */
/* 	  (double)(the_pedigree_stats->d2.n+1)/(double)(the_pedigree_stats->d2.d+1) : 2); */
/* } */

void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->agmr12.d, (double)(the_pedigree_stats->agmr12.n+1)/(double)(the_pedigree_stats->agmr12.d+1));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par1_hgmr.d, (double)(the_pedigree_stats->par1_hgmr.n+1)/(double)(the_pedigree_stats->par1_hgmr.d+1));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par1_R.d, (double)(the_pedigree_stats->par1_R.n+1)/(double)(the_pedigree_stats->par1_R.d+1));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par2_hgmr.d, (double)(the_pedigree_stats->par2_hgmr.n+1)/(double)(the_pedigree_stats->par2_hgmr.d+1));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par2_R.d, (double)(the_pedigree_stats->par2_R.n+1)/(double)(the_pedigree_stats->par2_R.d+1));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->d.d, (double)(the_pedigree_stats->d.n+1)/(double)(the_pedigree_stats->d.d+1));
}

double get_agmr12(Pedigree_stats* p){
  return (p->agmr12.d > 0)? (double)(p->agmr12.n+1)/(double)(p->agmr12.d+1) : 2;
}
double get_hgmr1(Pedigree_stats* p){
  return (p->par1_hgmr.d > 0)? (double)(p->par1_hgmr.n+1)/(double)(p->par1_hgmr.d+1) : 2;
}
double get_R1(Pedigree_stats* p){
  return (p->par1_R.d > 0)? (double)(p->par1_R.n+1)/(double)(p->par1_R.d+1) : 2;
}
double get_hgmr2(Pedigree_stats* p){
  return (p->par2_hgmr.d > 0)? (double)(p->par2_hgmr.n+1)/(double)(p->par2_hgmr.d+1) : 2;
}
double get_R2(Pedigree_stats* p){
  return (p->par2_R.d > 0)? (double)(p->par2_R.n+1)/(double)(p->par2_R.d+1) : 2;
}
/* double get_d1(Pedigree_stats* p){ */
/*   return (p->d1.d > 0)? (double)(p->d1.n+1)/(double)(p->d1.d+1) : 2; */
/* } */
/* double get_d2(Pedigree_stats* p){ */
/*   return (p->d2.d > 0)? (double)(p->d2.n+1)/(double)(p->d2.d+1) : 2; */
/* } */
double get_d(Pedigree_stats* p){
  return (p->d.d > 0)? (double)(p->d.n+1)/(double)(p->d.d+1) : 2;
}

void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees){
  fprintf(fh, " %3ld  ", alt_pedigrees->size);
  for(long i=0; i<alt_pedigrees->size; i++){
    Pedigree* alt_pedigree = alt_pedigrees->a[i];
    fprintf(fh, "%20s %20s ", alt_pedigree->F->id->a, alt_pedigree->M->id->a);
    print_pedigree_stats(fh, alt_pedigree->pedigree_stats);
  } 
}

long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_ok_hgmr, double max_self_r, double max_ok_d){ // returns 2 for ok biparental, 1 for ok self, 0 for bad
  double agmr12 = get_agmr12(p);
    double hgmr1 = get_hgmr1(p);
  double r1 = get_R1(p);
  double hgmr2 = get_hgmr2(p);
  double r2 = get_R2(p);
  double d = get_d(p);
  long result = 0;
  /* fprintf(stderr, "%7.4lf %7.4lf %7.4lf %7.4lf    %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", */
  /* 	  max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d1, */
  /* 	  agmr12, hgmr1, r1, hgmr2, r2, d1); */
  if(agmr12 <= max_self_agmr12){ // pedigree says self (or parents very similar)
    if( /*(agmr12 <= max_self_agmr12) && */ (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 <= max_self_r) && (r2 <= max_self_r) && (d <= max_ok_d) ){
      result = 1;
    }
  }else{ // pedigree says biparental
    if( (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 > max_self_r) && (r2 > max_self_r) && (d <= max_ok_d) ){
      result = 2;
    }
  }
  //  fprintf(stderr, "pedigree_ok?: %ld\n", result);
  return result;
}

void free_pedigree(const Pedigree* the_pedigree){
    if(the_pedigree == NULL) return;
  free(the_pedigree->pedigree_stats);
  free((Pedigree*) the_pedigree);
}

// *****  Vpedigree  *****

Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset){

  char* line = NULL;
  size_t len = 0;
  ssize_t nread;
  char* saveptr = NULL;
  if((nread = getline(&line, &len, p_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    if((token == NULL)  || (strcmp(token, "Accession") != 0)){
      exit(EXIT_FAILURE);
    }
  }
  Vpedigree* pedigrees = construct_vpedigree(1000);
  while((nread = getline(&line, &len, p_stream)) != -1){
    Vstr* fields = construct_vstr(PEDIGREE_FIELDS);
    char* token = strtok_r(line, "\t \n\r", &saveptr);
 
    add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    }
    // construct a Pedigree struct from last 3 fields  
    char* acc_id = ith_str_from_vstr(fields, -3); // ith_str ... copies the string, i.e. allocates more memory
    char* fempar_id = ith_str_from_vstr(fields, -2);
    char* malpar_id = ith_str_from_vstr(fields, -1);   
  
    long acc_idx, fempar_idx, malpar_idx;
    if(
       (strcmp(acc_id, "NA") != 0) &&
       (strcmp(fempar_id, "NA") != 0) &&
       (strcmp(malpar_id, "NA") != 0) &&
       ((acc_idx = index_of_id_in_vidxid(the_vidxid, acc_id)) != -1) &&
       ((fempar_idx = index_of_id_in_vidxid(the_vidxid, fempar_id)) != -1) &&
       ((malpar_idx = index_of_id_in_vidxid(the_vidxid, malpar_id)) != -1)
       ){
     
      Accession* Acc = the_gtsset->accessions->a[acc_idx]; // id->index;
      Accession* Fpar = the_gtsset->accessions->a[fempar_idx]; // id->index;
      Accession* Mpar = the_gtsset->accessions->a[malpar_idx]; // id->index;
     
      Pedigree* a_pedigree = construct_pedigree(Acc, Fpar, Mpar);
      add_pedigree_to_vpedigree(pedigrees, a_pedigree);
    }
    free_vstr(fields);
  } // done reading all lines
  free(line); // only needs to be freed once.
  // fprintf(stderr, "# size of Vpedigree pedigrees: %ld \n", pedigrees->size);
  return pedigrees;
}

Vpedigree* construct_vpedigree(long cap){
  Vpedigree* the_vped = (Vpedigree*)malloc(sizeof(Vpedigree));
  the_vped->capacity = cap;
  the_vped->size = 0;
  the_vped->a = (Pedigree**)malloc(the_vped->capacity*sizeof(Pedigree*));
}
  
void add_pedigree_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped){
  long cap = the_vped->capacity;
  long n = the_vped->size;
  if(n == cap){
    cap *= 2;
    the_vped->a = (Pedigree**)realloc(the_vped->a, cap*sizeof(Pedigree*));
    the_vped->capacity = cap;
  }
  the_vped->a[n] = the_ped;
  the_vped->size++;
}

Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_d){
  long n_parents = parent_idxs->size;

  // double max_ok_hgmr = 0.05;

  char* acc_id = the_pedigree->A->id->a; //Accession->id;
  long acc_idx = the_pedigree->A->index; // Accession->index;
  char* acc_gts = the_pedigree->A->genotypes->a; // the_gtsset->genotype_sets->a[the_pedigree->Accession->index];
  long fparent_idx = the_pedigree->F->index; //Fparent->index;
  char* fparent_gts = the_pedigree->F->genotypes->a; // the_gtsset->genotype_sets->a[fparent_idx];
  long mparent_idx = the_pedigree->M->index; //  Mparent->index;
  char* mparent_gts = the_pedigree->M->genotypes->a; // the_gtsset->genotype_sets->a[mparent_idx];


  // get best candidate parents on basis of hgmr (plus those in pedigree)
  // the_pedigree parent_idxs the_gtsset
  Vlong* best_parent_candidate_idxs = construct_vlong(10); 
  add_long_to_vlong(best_parent_candidate_idxs, fparent_idx); // add female parent (from pedigree)
  if(mparent_idx != fparent_idx) add_long_to_vlong(best_parent_candidate_idxs, mparent_idx); // add male parent (from pedigree) if distinct

  // sort accession indices by hgmr   
  Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_parents*sizeof(Idxhgmr));
  for(long i=0; i<parent_idxs->size; i++){
    long idx = parent_idxs->a[i];
    char* pgts = the_gtsset->accessions->a[idx]->genotypes->a; // genotype_sets->a[idx];
    the_idxhgmrs[i].idx = idx;
    the_idxhgmrs[i].hgmr = hgmr(acc_gts, pgts); //the_hgmr; //the_idxhgmr = {idx, the_hgmr);
  }
  sort_idxhgmr_by_hgmr(n_parents, the_idxhgmrs);
  //fprintf(stderr, "n_parents: %ld \n", n_parents);
  for(long i=0; i<n_parents; i++){ // store 'good' hgmr indices 
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != acc_idx){
      double the_hgmr = the_idxhgmrs[i].hgmr;
      if(the_hgmr >= max_ok_hgmr) break; 
      //   if(the_hgmr >= 0){
      if(the_idx != fparent_idx  &&  the_idx != mparent_idx  &&  the_idx != acc_idx){
	add_long_to_vlong(best_parent_candidate_idxs, the_idx);
      }     
      // }
    }
  }
  free(the_idxhgmrs);

  long ub = long_min(best_parent_candidate_idxs->size, 8); // set the number of possible parents to consider.
  // fprintf(stderr, "XXX: %8.4lf %8.4lf  %ld \n", max_ok_hgmr, max_ok_d1, ub); 
  Vpedigree* alt_pedigrees = construct_vpedigree(10);
  for(long i=0; i<ub; i++){
    long idx1 = best_parent_candidate_idxs->a[i];
    Accession* acc1 = the_gtsset->accessions->a[idx1];
    char* id1 = acc1->id->a; // accqession_ids->a[idx1];
    char* gts1 = acc1->genotypes->a; // genotype_sets->a[idx1];
      for(long j=i; j<ub; j++){
	long idx2 = best_parent_candidate_idxs->a[j];
	Accession* acc2 = the_gtsset->accessions->a[idx2];
	char* id2 = acc2->id->a; // _ids->a[idx2];
	if(! ((idx1 == fparent_idx && idx2 == mparent_idx) || (idx1 == mparent_idx && idx2 == fparent_idx))){  
	  char* gts2 = acc2->genotypes->a; // the_gtsset->genotype_sets->a[idx2];
	  Pedigree* alt_pedigree = construct_pedigree(the_pedigree->A, acc1, acc2); // arbitrarily put acc1 as Female parent, acc2 as male
	  Pedigree_stats* alt_pedigree_stats = triple_counts(gts1, gts2, acc_gts);
	  //if(get_hgmr1(alt_pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree_stats) <= 0.05  &&
	  if(get_d(alt_pedigree_stats) <= max_ok_d){   
	    alt_pedigree->pedigree_stats = alt_pedigree_stats;
	    add_pedigree_to_vpedigree(alt_pedigrees, alt_pedigree);
	  }
	}
      }
  }
  free_vlong(best_parent_candidate_idxs);
  return alt_pedigrees;
}


void free_vpedigree(const Vpedigree* the_vped){
    if(the_vped == NULL) return;
  for(long i=0; i<the_vped->size; i++){
    free_pedigree(the_vped->a[i]);
  }
  free(the_vped->a);
  free((Vpedigree*)the_vped);
}

// *****  sorting an array of Idxhgmr  *****
int cmpidxhgmr(const void* v1, const void* v2){
  const Idxhgmr* s1 = (const Idxhgmr*)v1;
  const Idxhgmr* s2 = (const Idxhgmr*)v2;
  return (s1->hgmr < s2->hgmr)? -1 : (s1->hgmr > s2->hgmr)? 1 : 0;
}

void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array){ // sort in place
  qsort(array, size, sizeof(Idxhgmr), cmpidxhgmr);
}

long long_min(long a, long b){
  return (a <= b)? a : b;
}

long long_max(long a, long b){
  return (a >= b)? a : b;
}

long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset){
  for(long i=0; i<the_gtsset->n_accessions; i++){
    char* id = the_gtsset->accessions->a[i]->id->a;
    long idx = index_of_id_in_vidxid(vidxid, id);
    if(idx != i) return 0;
  }
  return 1;
}


/* void print_pedigree_alternatives(FILE* fh, const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs){ */
/*   long n_parents = parent_idxs->size; */

/*   double max_ok_hgmr = 0.05; */

/*   char* acc_id = the_pedigree->A->id->a; //Accession->id; */
/*   long acc_idx = the_pedigree->A->index; // Accession->index; */
/*   char* acc_gts = the_pedigree->A->genotypes->a; // the_gtsset->genotype_sets->a[the_pedigree->Accession->index]; */
/*   long fparent_idx = the_pedigree->F->index; //Fparent->index; */
/*   char* fparent_gts = the_pedigree->F->genotypes->a; // the_gtsset->genotype_sets->a[fparent_idx]; */
/*   long mparent_idx = the_pedigree->M->index; //  Mparent->index; */
/*   char* mparent_gts = the_pedigree->M->genotypes->a; // the_gtsset->genotype_sets->a[mparent_idx]; */


/*   // get best candidate parents on basis of hgmr (plus those in pedigree) */
/*   // the_pedigree parent_idxs the_gtsset */
/*   Vlong* best_parent_candidate_idxs = construct_vlong(10);  */
/*   add_long_to_vlong(best_parent_candidate_idxs, fparent_idx); // add female parent (from pedigree) */
/*   if(mparent_idx != fparent_idx) add_long_to_vlong(best_parent_candidate_idxs, mparent_idx); // add male parent (from pedigree) if distinct */

/*   // sort accession indices by hgmr    */
/*   Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_parents*sizeof(Idxhgmr)); */
/*   for(long i=0; i<parent_idxs->size; i++){ */
/*     long idx = parent_idxs->a[i]; */
/*     char* pgts = the_gtsset->accessions->a[idx]->genotypes->a; // genotype_sets->a[idx]; */
/*     the_idxhgmrs[i].idx = idx; */
/*     the_idxhgmrs[i].hgmr = hgmr(acc_gts, pgts); //the_hgmr; //the_idxhgmr = {idx, the_hgmr); */
/*   } */
/*   sort_idxhgmr_by_hgmr(n_parents, the_idxhgmrs); */
 
/*   for(long i=0; i<n_parents; i++){ // store 'good' hgmr indices  */
/*     long the_idx = the_idxhgmrs[i].idx; */
/*     if(the_idx != acc_idx){ */
/*       double the_hgmr = the_idxhgmrs[i].hgmr; */
/*       if(the_hgmr >= 0.06) break;  */
/*       if(the_hgmr >= 0){ */
/* 	if(the_idx != fparent_idx  &&  the_idx != mparent_idx  &&  the_idx != acc_idx){ */
/* 	  add_long_to_vlong(best_parent_candidate_idxs, the_idx); */
/* 	}      */
/*       } */
/*     } */
/*   } */
/*   free(the_idxhgmrs); */

/*   long ub = long_min(best_parent_candidate_idxs->size, 8); // set the number of possible parents to consider. */
/*   Vpedigree* alt_pedigrees = construct_vpedigree(10); */
/*   for(long i=0; i<ub; i++){ */
/*     long idx1 = best_parent_candidate_idxs->a[i]; */
/*     Accession* acc1 = the_gtsset->accessions->a[idx1]; */
/*     char* id1 = acc1->id->a; // accession_ids->a[idx1]; */
/*     char* gts1 = acc1->genotypes->a; // genotype_sets->a[idx1]; */
/*     for(long j=i; j<ub; j++){ */
/*       long idx2 = best_parent_candidate_idxs->a[j]; */
/*       Accession* acc2 = the_gtsset->accessions->a[idx2]; */
/*       char* id2 = acc2->id->a; // _ids->a[idx2]; */
/*       char* gts2 = acc2->genotypes->a; // the_gtsset->genotype_sets->a[idx2]; */
/*       Pedigree* alt_pedigree = construct_pedigree(the_pedigree->A, acc1, acc2); // arbitrarily put acc1 as Female parent, acc2 as male */
/*       Pedigree_stats* alt_pedigree_stats = triple_counts(gts1, gts2, acc_gts); */
/*       alt_pedigree->pedigree_stats = alt_pedigree_stats; */
/*       add_pedigree_to_vpedigree(alt_pedigrees, alt_pedigree);    */
/*     } */
/*   } */
/*  free_vlong(best_parent_candidate_idxs); */
 
/*   // print */
/*   for(long i=0; i<alt_pedigrees->size; i++){ */
/*     Pedigree* alt_pedigree = alt_pedigrees->a[i]; */
/*     if(get_hgmr1(alt_pedigree->pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree->pedigree_stats) <= 0.05  &&  get_d1(alt_pedigree->pedigree_stats) <= 0.015){    */
/*       fprintf(fh, "%s %s ", alt_pedigree->F->id->a, alt_pedigree->M->id->a); */
/*       print_pedigree_stats(fh, alt_pedigree->pedigree_stats); */
/*     } */
/*   }  */
/* } */
