#include <stdlib.h>
#include <stdio.h>
#include "gtset.h"
// #include "vect.h"
#include "pedigree.h"
// #include "various.h"

#define PEDIGREE_FIELDS 7 // number of whitespace separated fields in pedigree file, with ids in last 3 fields.

// extern int do_checks_flag; // option -c sets this to 1 to do some checks.

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent){ //IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid){
  Pedigree* the_pedigree = (Pedigree*)malloc(sizeof(Pedigree));
  the_pedigree->F = Fparent;
  the_pedigree->M = Mparent;
  the_pedigree->A = Acc;
  the_pedigree->pedigree_stats = NULL;
  //  fprintf(stderr, "# in construct_pedigree, just before return.\n");
  return the_pedigree;
}

/* double hgmr(char* gts1, char* gts2){ */
/*   char c1, c2; */
/*   long n_numer = 0; */
/*   long n_denom = 0; */
/*   long i=0; */
/*   while((c1 = gts1[i]) != '\0'){ */
/*     if((c1 == '0') || (c1 == '2')){ */
/*       c2 = gts2[i]; */
/*       if((c2 == '0') || (c2 == '2')){ */
/* 	n_denom++; */
/* 	if(c1 != c2) n_numer++; */
/*       } */
/*     } */
/*     i++; */
/*   } */
/*   fprintf(stderr, "hgmr   n,d: %ld %ld\n", n_numer, n_denom);  */
/*   return (n_denom > 0)? (double)n_numer/(double)n_denom : 2.0;   */
/* } */


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


two_longs gamete_dosage_range(long d, long ploidy){
  long half_ploidy = ploidy/2;
  long min_dosage, max_dosage;
  if(d <= half_ploidy){
    min_dosage = 0;
    max_dosage = d;
  }else{ // d > half_ploidy
    min_dosage = d - half_ploidy;
    max_dosage = half_ploidy;
  }
  two_longs result = {min_dosage, max_dosage};
  return result;
}

two_longs diploid_quick_and_dirty_triple_counts(Accession* acc1, Accession* acc2, Accession* progacc){
  long allowed_count = 0;
  long forbidden_count = 0;
  long count22_2a = 0;
  long count22_01a = 0;
  long n01_2 = 0;
  long n10_2 = 0;
  
  for(long i=0; i<acc1->alt_homozygs->size; i++){ // just look at markers with gt1 = 2
    long index = acc1->alt_homozygs->a[i];
    assert(acc1->genotypes->a[index] == '2');
    char gt2 = acc2->genotypes->a[index];
    char gtprog = progacc->genotypes->a[index];
    if(gt2 == MISSING_DATA_CHAR  ||  gtprog == MISSING_DATA_CHAR) continue;
    if(gt2 == '0'){
      if(gtprog == '1'){
	allowed_count++;
      }else{
	forbidden_count++;
      }
    }else if(gt2 == '1'){
      if(gtprog == '0'){
	forbidden_count++;
      }else{
	allowed_count++;
      }
    }else if(gt2 == '2'){
      if(gtprog == '2'){
	count22_2a++;
      }else{
	count22_01a++;
      }

    }
  }

  long count22_2b = 0;
  long count22_01b = 0;
  for(long i=0; i<acc2->alt_homozygs->size; i++){
    long index = acc2->alt_homozygs->a[i];
    char gt1 = acc1->genotypes->a[index];
    char gtprog = progacc->genotypes->a[index];
    if(gt1 == MISSING_DATA_CHAR  ||  gtprog == MISSING_DATA_CHAR) continue;
    if(gt1 == '0'){
      if(gtprog == '1'){
	allowed_count++;
      }else{
	forbidden_count++;
      }
    }else if(gt1 == '1'){
      if(gtprog == '0'){
	forbidden_count++;
      }else{
	allowed_count++;
      }
    }else if(gt1 == '2'){
      if(gtprog == '2'){
	count22_2b++;
      }else{
	count22_01b++;
      }
    }
  }
  assert(count22_2a == count22_2b);
  assert(count22_01a == count22_01b);
  forbidden_count += count22_01a;
  allowed_count += count22_2a;
  // fprintf(stderr, "%ld %ld  %ld %ld   %ld %ld \n", count22_2a, count22_01a, count22_2b, count22_01b, allowed_count, forbidden_count);
  two_longs result = {forbidden_count, allowed_count};
  return result;
}

four_longs q_and_d_n22x_diploid(Accession* acc1, Accession* acc2, Accession* progacc){
  long n22_0 = 0;
  long n22_1 = 0;
  long n22_2 = 0;
  long n00_0 = 0;
  long n00_1 = 0;
  long n00_2 = 0;
  long n00 = 0;
  long n22 = 0;
  
  for(long i=0; i<acc1->alt_homozygs->size; i++){ // just look at markers with gt1 = 2
    long index = acc1->alt_homozygs->a[i];
    assert(acc1->genotypes->a[index] == '2');
    char gt2 = acc2->genotypes->a[index];
    if(gt2 == MISSING_DATA_CHAR) continue;
    char gtprog = progacc->genotypes->a[index];
    if(gtprog == MISSING_DATA_CHAR) continue;
    if(gt2 == '2'){
      n22++;
      if(gtprog == '2'){
	n22_2++;
      }else if(gtprog == '1'){
	n22_1++;
      }else if(gtprog == '0'){
	n22_0++;
      }
    }
  }
 
  for(long i=0; i<acc1->ref_homozygs->size; i++){ // just look at markers with gt1 = 2
    long index = acc1->ref_homozygs->a[i];
    assert(acc1->genotypes->a[index] == '0');
    char gt2 = acc2->genotypes->a[index];
    if(gt2 == MISSING_DATA_CHAR) continue;
    char gtprog = progacc->genotypes->a[index];
    if(gtprog == MISSING_DATA_CHAR) continue;
    if(gt2 == '0'){
      n00++;
      if(gtprog == '0'){
	n00_0++;
      }else if(gtprog == '1'){
	n00_1++;
      }else if(gtprog == '2'){
	n00_2++;
      }
    }
    if(n00 >= 400) break;
  }

  //fprintf(stderr, "%ld %ld %ld  %ld    %ld %ld %ld  %ld \n", n00_0, n00_1, n00_2, n00, n22_0, n22_1, n22_2, n22);
  //   fprintf(stderr, "%ld %ld %ld  %ld \n", n00_0 + n22_2, n00_1 + n22_1, n00_2 + n22_0, n00 + n22);
  n22_0 += n00_2;
  n22_1 += n00_1;
  n22_2 += n00_0;
 
  four_longs result = {n22_0, n22_1, n22_2, n22_0 + n22_1 + n22_2};
  return result;
}
    
ND tfc_tetraploid(char* gts1, char* gts2, char* proggts){ // 
  char c1, c2, c3;
  long numer = 0; // small if gts1, gts2 parents of proggts
  long denom = 0;
 
  long i = 0;
  while((c3 = proggts[i]) != '\0'){ // go until hit null termination of proggts
    if(c3 == MISSING_DATA_CHAR) {i++; continue;}
    if((c1 = gts1[i]) == MISSING_DATA_CHAR) {i++; continue;}
    if((c2 = gts2[i]) == MISSING_DATA_CHAR) {i++; continue;}
    long progd = c3 - 48;
    long p1d = c1 - 48;
    long p2d = c2 - 48;
    if(p1d == 0){
      if(p2d == 0){ // x,0,0
	denom++;
	if(progd == 1 || progd == 2) numer++;
      }else if(p2d == 1){ // x,0,1
	denom++;
	if(progd == 2) numer++;
      }
    }else if(p1d == 1){
      if(p2d == 0){ // x,1,0
	denom++;
	if(progd == 2) numer++;
      }else if(p2d == 1){ // x,1,1
	denom++;
	if(progd == 3) numer++;
      }
    }else if(p1d + p2d  <= 5){
      // do nothing
    }else if(p1d == 3){
      if(p2d == 3){ // x,3,3
	denom++;
	if(progd == 1) numer++;
      }else if(p2d == 4){ // x,3,4
	denom++;
	if(progd == 2) numer++;
      }
    }else if(p1d == 4){
      if(p2d == 3){ // x,4,3
	denom++;
	if(progd == 2) numer++;
      }else if(p2d == 4){ // x,1,1
	denom++;
	if(progd == 2 || progd == 3) numer++;
      }
    }
    i++;
  }
  ND result = {numer, denom};
  return result;
}


ND tfc_diploid(char* gts1, char* gts2, char* proggts){ // just (n100 + n122)/(nx00 + nx22)
  char c1, c2, c3;
  long numer = 0; // small if gts1, gts2 parents of proggts
  long denom = 0;
 
  long i = 0;
  while((c3 = proggts[i]) != '\0'){ // go until hit null termination of proggts
    if(c3 == MISSING_DATA_CHAR) {i++; continue;}
    if((c1 = gts1[i]) == MISSING_DATA_CHAR) {i++; continue;}
    if((c2 = gts2[i]) == MISSING_DATA_CHAR) {i++; continue;}
    long progd = c3 - 48;
    long p1d = c1 - 48;
    long p2d = c2 - 48;
    if(p1d == 0  &&  p2d == 0){
      denom++;
      if(progd == 1) numer++;
    }else if(p1d == 2  &&  p2d == 2){
      denom++;
      if(progd == 1) numer++;
    }
    i++;
  }
  ND result = {numer, denom};
  return result;
}

ND TFC(char* gts1, char* gts2, char* proggts, long ploidy){
  ND result = {-1, -1};
  if(ploidy == 2){
    result = tfc_diploid(gts1, gts2, proggts);
  }else if(ploidy == 4){
    result = tfc_tetraploid(gts1, gts2, proggts);
  }else{
    fprintf(stderr, "# TFC not implemented for ploidy = %ld \n", ploidy);
  }
  return result;
}

four_longs tfca(char* gts1, char* gts2, char* proggts, long ploidy){
  char c1, c2, c3;
  long f1_count = 0; // n002 + n220
  long f2_count = 0; // small if gts1, gts2 parents of proggts
  long f1_denom = 0;
  long f2_denom = 0;
  long i = 0;
  while((c3 = proggts[i]) != '\0'){ // go until hit null termination of proggts
    if(c3 == MISSING_DATA_CHAR) {i++; continue;}
    if((c1 = gts1[i]) == MISSING_DATA_CHAR) {i++; continue;}
    if((c2 = gts2[i]) == MISSING_DATA_CHAR) {i++; continue;}
    long progd = c3 - 48;
    long p1d = c1 - 48;
    long p2d = c2 - 48;
   
    if(p1d == 1  &&  p2d == 1) {i++; continue;}
    f1_denom++;
    if(p1d == 0  && p2d == 0){
      f2_denom++;
      if(progd == 2){
	f1_count++; // n002
      }else if(progd == 1){
	f2_count++;
      }
	  
    }else if(p1d == 2  &&  p2d == 2){
      f2_denom++;
      if(progd == 0){
	f1_count++; // n220
	
      }else if(progd == 1){
	f2_count++;
      }
    }else{
      if(p1d > p2d){ long t = p1d; p1d = p2d; p2d = t;}
      if(p1d == 0){
	if(p2d == 1){
	  if(progd == 2){
	    f1_count++;
	  }
	}else if(p2d == 2){
	  if(progd != 1){
	    f1_count++;
	  }
	}
      }else if(p1d == 1){
	if(p2d == 2){
	  if(progd == 0){
	    f1_count++;
	  }
	}
      }
    }
    i++;
  } // end of loop
  four_longs result = {f1_count, f1_denom, f2_count, f2_denom};
  return result;
}


four_longs triple_forbidden_counts(char* gts1, char* gts2, char* proggts, long ploidy){ 
  // version for any (even) ploidy, counts all triples that shouldn't happen if
  // gts1, gts2 are parents of proggts.
  char c1, c2, c3;
  long i = 0;
  long forbidden_count = 0;
  long forbidden1_count = 0;
  long forbidden2_count = 0;
  long no_md_count = 0;
  long not_po2po2_count = 0;
  while((c3 = proggts[i]) != '\0'){ // go until hit null termination of proggts
    if(c3 == MISSING_DATA_CHAR) {i++; continue;}
    if((c1 = gts1[i]) == MISSING_DATA_CHAR) {i++; continue;}
    if((c2 = gts2[i]) == MISSING_DATA_CHAR) {i++; continue;}
    long progeny_dosage = (long)(c3-48);
    two_longs fgam_d_range = gamete_dosage_range(c1-48, ploidy);
    two_longs mgam_d_range = gamete_dosage_range(c2-48, ploidy);
    //  fprintf(stderr, "# %ld %ld  gamete ranges: f: %ld %ld  m: %ld %ld \n", (long)(c1-48), (long)(c2-48), fgam_d_range.l1, fgam_d_range.l2, mgam_d_range.l1, mgam_d_range.l2);
    long min_prog_d = fgam_d_range.l1 + mgam_d_range.l1;
    long max_prog_d = fgam_d_range.l2 + mgam_d_range.l2;
    if(progeny_dosage < min_prog_d  ||  progeny_dosage > max_prog_d){
      forbidden_count++;
      if((min_prog_d - progeny_dosage == 1) || (progeny_dosage - max_prog_d == 1)) forbidden1_count++;
      if((min_prog_d - progeny_dosage == 2) || (progeny_dosage - max_prog_d == 2)) forbidden2_count++;
    }
    no_md_count++;
    if( (2*(c1-48) != ploidy) || (2*(c2-48) != ploidy) ) not_po2po2_count++;
    i++;
    //  fprintf(stderr, "# %ld %ld \n", i, no_md_count);
  }
  // fprintf(stderr, "### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
  four_longs result = {forbidden_count, forbidden1_count,
		       //forbidden2_count,
		       no_md_count, not_po2po2_count};
  return result;
}

Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts, long ploidy){ // 
  //  long* d0counts, long* d1counts, long* d2counts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){

  char c1, c2, c3;
  char alt_char =(char)(ploidy + 48); // character corresponding to alt allele homozygous
  // distinguish 46 types of triples:
  // here 0 means homozygous (ref allele)
  // 2 means homozygous (alt allele) i.e. dosage == ploidy
  // 1 means heterozygous, i.e. any dosage from 1 to ploidy-1
  // 3 in, e.g. n_03_1 means missing data
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0; // n_01_0 <-> parent1 has 0, parent2 has 1, and progeny has 0
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0; // n_03_0 <-> parent 1 has 0, parent 2 has missing data, progeny has 0
  
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

  long n_33_x = 0; // both parents have md, progeny anything except md.
  long n_xy_3 = 0; // progeny md, parents anything.
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){ // go until hit null termination of proggts
    // fprintf(stderr, "#zzz: %c %c %c \n", c3, gts1[i], gts2[i]);
    if(c3 != MISSING_DATA_CHAR){
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
	else if(c2 == MISSING_DATA_CHAR){ // 03y
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
	}else if(c2 == MISSING_DATA_CHAR){ // 13y
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
	}else if(c2 == MISSING_DATA_CHAR){ // 23y
	  if(c3 == '0'){
	    n_23_0++;
	  }else if(c3 == '1'){
	    n_23_1++;
	  }else if(c3 == '2'){
	    n_23_2++;
	  }
	}
      }else if(c1 == MISSING_DATA_CHAR){ // 3xy
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
	  n_33_x++;
	}
      }
    }else{ // c3 == MISSING_DATA_CHAR
      n_xy_3++;
    }
    i++;
  }
  long n_0 = // 15 triples consistent with true parents-offspring relationship,
    // with no genotyping errors
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 +
    n_11_0 + n_11_1 + n_11_2 +
    n_12_1 + n_12_2 +
    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 = // 10 triple requiring one (0<->1 or 1<->2) error for consistency 
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 or 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0;
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1;
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2;
  
  long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0;
  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1;
  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2;
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0;
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1;
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2;

  long n_3x_0 = n_30_0 + n_31_0 + n_32_0; // + n_03_0;
  long n_3x_1 = n_30_1 + n_31_1 + n_32_1; // + n_03_1;
  long n_3x_2 = n_30_2 + n_31_2 + n_32_2; // + n_03_2;

  long z_numer = n_00_1 + n_22_1;
  long z_denom = z_numer + n_00_0 + n_22_2 + n_00_2 + n_22_0;

  long total =
    n_0x_0 + n_0x_1 + n_0x_2
    + n_1x_0 + n_1x_1 + n_1x_2
    + n_2x_0 + n_2x_1 + n_2x_2
    + n_3x_0 + n_3x_1 + n_3x_2
    + n_33_x + n_xy_3;

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

  ND z_nd = {z_numer, z_denom};
  pedigree_stats->z = z_nd; 

  return pedigree_stats;
}

/* Pedigree_stats* triple_counts_x(char* gts1, char* gts2, char* proggts, //  */
/* 				long* d0counts, long* d1counts, long* d2counts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){ */

/*   char c1, c2, c3; */
/*   long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0; */
/*   long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0; */
/*   long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0; */
/*   long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0; */
  
/*   long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0; */
/*   long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0; */
/*   long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0; */
/*   long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0; */
  
/*   long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0; */
/*   long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0; */
/*   long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0; */
/*   long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0; */

/*   long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0; */
/*   long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0; */
/*   long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0; */
  
/*   long i=0; */
/*   while((c3 = proggts[i]) != '\0'){ */
/*     if(c3 != '3'){ */
/*       c1 = gts1[i]; */
/*       c2 = gts2[i]; */
/*       if(c1 == '0'){ // 0xy */
/* 	if(c2 == '0'){ // 00y */
/* 	  if(c3 == '0'){  */
/* 	    n_00_0++; d0counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_00_1++; d1counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_00_2++; d2counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '1'){ // 01y */
/* 	  if(c3 == '0'){  */
/* 	    n_01_0++; d0counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_01_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_01_2++; d1counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '2'){ // 02y */
/* 	  if(c3 == '0'){ */
/* 	    n_02_0++; d1counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_02_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_02_2++; d1counts[i]++; */
/* 	  } */
/* 	} */
/* 	else if(c2 == '3'){ // 03y */
/* 	  if(c3 == '0'){ */
/* 	    n_03_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_03_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_03_2++; */
/* 	  } */
/* 	} */
/*       }else if(c1 == '1'){ // 1xy */
/* 	if(c2 == '0'){ // 10y */
/* 	  if(c3 == '0'){ */
/* 	    n_10_0++; d0counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_10_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_10_2++; d1counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '1'){ // 11y */
/* 	  if(c3 == '0'){ */
/* 	    n_11_0++; d0counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_11_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_11_2++; d0counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '2'){ // 12y */
/* 	  if(c3 == '0'){ */
/* 	    n_12_0++; d1counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_12_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_12_2++; d0counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '3'){ // 13y */
/* 	  if(c3 == '0'){ */
/* 	    n_13_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_13_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_13_2++; */
/* 	  } */
/* 	} */
/*       }else if(c1 == '2'){ // 2xy */
/* 	if(c2 == '0'){ // 20y */
/* 	  if(c3 == '0'){ */
/* 	    n_20_0++; d1counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_20_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_20_2++; d1counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '1'){ // 21y */
/* 	  if(c3 == '0'){ */
/* 	    n_21_0++; d1counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_21_1++; d0counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_21_2++; d0counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '2'){ // 22y */
/* 	  if(c3 == '0'){ */
/* 	    n_22_0++; d2counts[i]++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_22_1++; d1counts[i]++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_22_2++; d0counts[i]++; */
/* 	  } */
/* 	}else if(c2 == '3'){ // 23y */
/* 	  if(c3 == '0'){ */
/* 	    n_23_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_23_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_23_2++; */
/* 	  } */
/* 	} */
/*       }else if(c1 == '3'){ // 3xy */
/* 	if(c2 == '0'){ // 30y */
/* 	  if(c3 == '0'){ */
/* 	    n_30_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_30_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_30_2++; */
/* 	  } */
/* 	}else if(c2 == '1'){ // 31y */
/* 	  if(c3 == '0'){ */
/* 	    n_31_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_31_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_31_2++; */
/* 	  } */
/* 	}else if(c2 == '2'){ // 32y */
/* 	  if(c3 == '0'){ */
/* 	    n_32_0++; */
/* 	  }else if(c3 == '1'){ */
/* 	    n_32_1++; */
/* 	  }else if(c3 == '2'){ */
/* 	    n_32_2++; */
/* 	  } */
/* 	}else{ // 33y */
/* 	  // do nothing */
/* 	} */
/*       } */
/*     } */
/*     i++; */
/*   } */
/*   long n_0 = */
/*     n_00_0 + n_01_0 + n_01_1 + n_02_1 + */
/*     n_10_0 + n_10_1 + n_11_0 + n_11_1 + n_11_2 + n_12_1 + n_12_2 + */

/*     n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case */
/*   long n_1 = */
/*     n_00_1 + n_01_2 + n_02_0 + n_02_2 + */
/*     n_10_2 + n_12_0 + */
/*     n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 of 1<->2 type. */
/*   long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type. */

/*   //  double d1 = (n_0 > 0)? (double)n_1/(double)(n_0 + n_1) : 2; */
/*   //  double d2 = (n_0 > 0)? (double)n_2/(double)(n_0 + n_2) : 2; */
  
/*   long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0; */
/*   long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1; */
/*   long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2; */
  
/*   long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0; */
/*   long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1; */
/*   long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2; */
  
/*   long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0; */
/*   long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1; */
/*   long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2; */

/*   long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0; */
/*   long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1; */
/*   long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2; */
  
/*   long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0; */
/*   long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1; */
/*   long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2; */
  
/*   long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0; */
/*   long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1; */
/*   long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2; */

/*   long hgmr1_numer = n_0x_2 + n_2x_0; */
/*   long hgmr2_numer = n_x0_2 + n_x2_0; */
/*   long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer; */
/*   long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer; */
/*   // double hgmr1 = (hgmr1_denom > 0)? (double)hgmr1_numer/(double)hgmr1_denom : 2; */
/*   // double hgmr2 = (hgmr2_denom > 0)? (double)hgmr2_numer/(double)hgmr2_denom : 2; */

/*   //  long r0x1_numer = n_0x_1 + n_2x_1; */
/*   //  long r0x1_denom = r0x1_numer + n_0x_0 + n_2x_2; */
/*   //  long rx01_numer = n_x0_1 + n_x2_1; */
/*   //  long rx01_denom = rx01_numer + n_x0_0 + n_x2_2; */
/*   long r0x1or2_numer = n_0x_1 + n_0x_2 + n_2x_1 + n_2x_0; */
/*   long r0x1or2_denom = r0x1or2_numer + n_0x_0 + n_2x_2; */
/*   long rx01or2_numer = n_x0_1 + n_x0_2 + n_x2_1 + n_x2_0; */
/*   long rx01or2_denom = rx01or2_numer + n_x0_0 + n_x2_2; */
/*   // double r0x1 = (r0x1_denom > 0)? (double)r0x1_numer/(double)r0x1_denom : 2; */
/*   // double rx01 = (rx01_denom > 0)? (double)rx01_numer/(double)rx01_denom : 2; */

/*   long agmr12_numer = */
/*     n_01_0 + n_01_1 + n_01_2 + //n_01_3 + */
/*     n_02_0 + n_02_1 + n_02_2 + //n_02_3 + */
/*     n_10_0 + n_10_1 + n_10_2 + //n_10_3 + */
/*     n_12_0 + n_12_1 + n_12_2 + //n_12_3 + */
/*     n_20_0 + n_20_1 + n_20_2 + //n_20_3 + */
/*     n_21_0 + n_21_1 + n_21_2; //n_21_3; */
  
/*   long agmr12_denom = agmr12_numer + */
/*     n_00_0 + n_00_1 + n_00_2 + //n_00_3 + */
/*     n_11_0 + n_11_1 + n_11_2 + //n_11_3 + */
/*     n_22_0 + n_22_1 + n_22_2; //n_22_3; */
    
/*   /\* long agmr1_numer = hgmr1_numer + n_0x_1 + n_2x_1 + n_1x_0 + n_1x_2; *\/ */
/*   /\* long agmr2_numer = hgmr2_numer + n_x0_1 + n_x2_1 + n_x1_0 + n_x1_2; *\/ */
/*   /\* long agmr1_denom = agmr1_numer + n_0x_0 + n_1x_1 + n_2x_2; *\/ */
/*   /\* long agmr2_denom = agmr2_numer + n_x0_0 + n_x1_1 + n_x2_2; *\/ */
/*   /\* double agmr1 = (agmr1_denom > 0)? (double)agmr1_numer/(double)agmr1_denom : 2; *\/ */
/*   /\* double agmr2 = (agmr2_denom > 0)? (double)agmr2_numer/(double)agmr2_denom : 2; *\/ */
/*   /\* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr1,  hgmr1, r0x1); *\/ */
/*   /\* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr2,  hgmr2, rx01); *\/ */
/*   /\* fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); *\/ */

/*   Pedigree_stats* pedigree_stats = (Pedigree_stats*)malloc(sizeof(Pedigree_stats)); */
/*   ND agmr12_nd = {agmr12_numer, agmr12_denom}; */
/*   pedigree_stats->agmr12 = agmr12_nd; */
/*   ND hgmr1_nd = {hgmr1_numer, hgmr1_denom}; */
/*   pedigree_stats->par1_hgmr = hgmr1_nd; */
/*   // fprintf(stderr, "  %7.5f ", (hgmr1_denom > 0)? (double)hgmr1_numer/hgmr1_denom : 2); */
/*   //  ND r1_nd = {r0x1_numer, r0x1_denom}; */
/*   //  pedigree_stats->par1_r = r1_nd; */
/*   ND R1_nd = {r0x1or2_numer, r0x1or2_denom}; */
/*   pedigree_stats->par1_R = R1_nd; */
  
/*   ND hgmr2_nd = {hgmr2_numer, hgmr2_denom}; */
/*   pedigree_stats->par2_hgmr = hgmr2_nd;  */
/*   //  ND r2_nd = {rx01_numer, rx01_denom}; */
/*   //  pedigree_stats->par2_r = r2_nd; */
/*   ND R2_nd = {rx01or2_numer, rx01or2_denom}; */
/*   pedigree_stats->par2_R = R2_nd; */
  
/*   /\* ND d1_nd = {n_1, n_0 + n_1}; *\/ */
/*   /\* pedigree_stats->d1 = d1_nd; *\/ */
/*   /\* ND d2_nd = {n_2, n_0 + n_2}; *\/ */
/*   /\* pedigree_stats->d2 = d2_nd; *\/ */

/*   ND d_nd = {n_1 + n_2, n_0 + n_1 + n_2}; */
/*   pedigree_stats->d = d_nd; */
  
/*   /\* if(hgmr1 <= max_hgmr  &&  hgmr2 <= max_hgmr  &&  d1 <= max_d1){ *\/ */
/*   /\*   fprintf(stderr, "  %-30s %-30s  ", id1, id2); *\/ */
/*   /\*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr1, r0x1); *\/ */
/*   /\*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr2, rx01); *\/ */
/*   /\*   fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); *\/ */
/*   /\*   fprintf(stderr, "\n"); *\/ */
/*   /\* } *\/ */

/*   return pedigree_stats; */
/* } */



/* void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset){ */

/*   char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index]; */
/*   char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index]; */
/*   char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index]; */
    
/*   Ahr fp_ahr, mp_ahr, fm_ahr; */
/*   agmr_hgmr_r(fempar_gts, acc_gts, &(the_pedigree->fp_ahr)); */
/*   agmr_hgmr_r(malpar_gts, acc_gts, &(the_pedigree->mp_ahr)); */
/*   agmr_hgmr_r(fempar_gts, malpar_gts, &(the_pedigree->fm_ahr)); */
/* } */
Pedigree_stats* construct_pedigree_stats(Pedigree* the_pedigree, long ploidy){
  Pedigree_stats* the_ps = (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
  the_ps->agmr12.n = 0;
  the_ps->agmr12.d = 0;
  the_ps->par1_hgmr.n = 0;
  the_ps->par1_hgmr.d = 0;
  the_ps->par1_R.n = 0;
  the_ps->par1_R.d = 0;
  the_ps->par2_hgmr.n = 0;
  the_ps->par2_hgmr.d = 0;
  the_ps->par2_R.n = 0;
  the_ps->par2_R.d = 0;
  //  the_ps->d.n = 0;
  //  the_ps->d.d = 0;
  the_ps->pseudo_hgmr.n = 0;
  the_ps->pseudo_hgmr.d = 0;  
}
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, long ploidy){ //, long* d0counts, long* d1counts, long* d2counts){ //, GenotypesSet* the_gtsset){

  assert(the_pedigree->F != NULL  ||  the_pedigree->M != NULL); // shouldn't have both parents NULL //
  if(the_pedigree->F != NULL  &&  the_pedigree->M != NULL){
    return triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a, ploidy );
  }else{ // one of the parents is NULL 
    Pedigree_stats* the_ps = (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
    the_ps->agmr12.n = 0;
    the_ps->agmr12.d = 0;
    the_ps->z.n = 0;
    the_ps->z.d = 0;
    if(the_pedigree->F != NULL){ // we have female parent id, no male parent id
      four_longs hgmrR = hgmr_R(the_pedigree->F->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
      the_ps->par1_hgmr.n = hgmrR.l1;
      the_ps->par1_hgmr.d = hgmrR.l2;
      the_ps->par1_R.n = hgmrR.l3;
      the_ps->par1_R.d = hgmrR.l4;
      the_ps->par2_hgmr.n = 0;
      the_ps->par2_hgmr.d = 0;
      the_ps->par2_R.n = 0;
      the_ps->par2_R.d = 0;

    }else{ // we have male parent id, no female parent id
      if(DO_ASSERT) assert(the_pedigree->M != NULL);
      the_ps->par1_hgmr.n = 0;
      the_ps->par1_hgmr.d = 0;
      the_ps->par1_R.n = 0;
      the_ps->par1_R.d = 0;
      four_longs hgmrR = hgmr_R(the_pedigree->M->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
      the_ps->par2_hgmr.n = hgmrR.l1;
      the_ps->par2_hgmr.d = hgmrR.l2;
        the_ps->par2_R.n = hgmrR.l3;
      the_ps->par2_R.d = hgmrR.l4;
    }
    //construct_pedigree_stats(the_pedigree, ploidy); //the_ps = (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
    // the_ps->agmr12.n = 0;
    //  fprintf(stderr, "######  %ld %ld \n", the_ps->par1_hgmr.n, the_ps->par2_hgmr.d);
    return the_ps;
  }
}
const Vlong* accessions_with_offspring(const Vpedigree* the_vped, long n_accessions){
  Vlong* offspring_counts = construct_vlong_zeroes(n_accessions);
  for(long i=0; i<the_vped->size; i++){
    const Pedigree* the_ped = the_vped->a[i];
    // fprintf(stderr, "%ld %ld \n", the_ped->F->index, the_ped->Fparent->index);
    if(the_ped->F != NULL) offspring_counts->a[the_ped->F->index]++;
    if(the_ped->M != NULL) offspring_counts->a[the_ped->M->index]++; 
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

void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->agmr12.d, n_over_d(the_pedigree_stats->agmr12));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par1_hgmr.d, n_over_d(the_pedigree_stats->par1_hgmr));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par1_R.d, n_over_d(the_pedigree_stats->par1_R));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par2_hgmr.d, n_over_d(the_pedigree_stats->par2_hgmr));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->par2_R.d, n_over_d(the_pedigree_stats->par2_R));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->z.d, n_over_d(the_pedigree_stats->z));
  fprintf(fh, "%5ld %6.5lf  ", the_pedigree_stats->d.d, n_over_d(the_pedigree_stats->d));
 
}

double n_over_d(ND nd){
  return (nd.d > 0)? (double)nd.n/(double)nd.d : 2;
}

void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees){
  fprintf(fh, " %3ld  ", alt_pedigrees->size);
  for(long i=0; i<alt_pedigrees->size; i++){
    Pedigree* alt_pedigree = alt_pedigrees->a[i];
    fprintf(fh, "%20s %20s ", alt_pedigree->F->id->a, alt_pedigree->M->id->a);
    print_pedigree_stats(fh, alt_pedigree->pedigree_stats);
  } 
}

long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_ok_hgmr, double max_self_r, double max_ok_z){ // returns 2 for ok biparental, 1 for ok self, 0 for bad
  double agmr12 = n_over_d(p->agmr12);
  double hgmr1 = n_over_d(p->par1_hgmr);
  double r1 = n_over_d(p->par1_R);
  double hgmr2 = n_over_d(p->par2_hgmr);
  double r2 = n_over_d(p->par2_R);
  // double d = n_over_d(p->d);
  double z = n_over_d(p->z);
  long result = 0;
  /* fprintf(stderr, "%7.4lf %7.4lf %7.4lf %7.4lf    %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", */
  /* 	  max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d1, */
  /* 	  agmr12, hgmr1, r1, hgmr2, r2, d1); */
  if(agmr12 <= max_self_agmr12){ // pedigree says self (or parents very similar)
    if( /*(agmr12 <= max_self_agmr12) && */ (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 <= max_self_r) && (r2 <= max_self_r) && (z <= max_ok_z) ){
      result = 1;
    }
  }else{ // pedigree says biparental
    if( (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 > max_self_r) && (r2 > max_self_r) && (z <= max_ok_z) ){
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
// read pedigree from file in 3 whitespace-separated cols format
Vpedigree* read_and_store_pedigrees_3col(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset){ 

  char* line = NULL;
  size_t len = 0;
  ssize_t nread;
  char* saveptr = NULL;
  /* if((nread = getline(&line, &len, p_stream)) != -1){ */
  /*   char* token = strtok_r(line, "\t \n\r", &saveptr); */
  /*   if((token == NULL)  || (strcmp(token, "Accession") != 0)){ */
  /*     exit(EXIT_FAILURE); */
  /*   } */
  /* } */
  Vpedigree* pedigrees = construct_vpedigree(1000);
  while((nread = getline(&line, &len, p_stream)) != -1){
    Vstr* fields = construct_vstr(3);
    char* token = strtok_r(line, "\t \n\r", &saveptr);
 
    add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    }
    // construct a Pedigree struct from first 3 fields  
    char* acc_id = ith_str_from_vstr(fields, 0); // ith_str ... copies the string, i.e. allocates more memory
    char* fempar_id = ith_str_from_vstr(fields, 1);
    char* malpar_id = ith_str_from_vstr(fields, 2);
    // fprintf(stderr, "### %s %s %s \n", acc_id, fempar_id, malpar_id);
  
    long acc_idx, fempar_idx, malpar_idx;
    if( // progeny id is not "NA" and is present in genotypes data set
	( strcmp(acc_id, "NA") != 0) &&
	((acc_idx = index_of_id_in_vidxid(the_vidxid, acc_id)) != ID_NA_INDEX)
	){ 
      fempar_idx = index_of_id_in_vidxid(the_vidxid, fempar_id);
      malpar_idx = index_of_id_in_vidxid(the_vidxid, malpar_id);
      
      if( (fempar_idx != ID_NA_INDEX) || (malpar_idx != ID_NA_INDEX) ){ // pedigree file has a valid id for at least 1 parent
	Accession* Acc = the_gtsset->accessions->a[acc_idx]; // id->index;
       	Accession* Fpar = (fempar_idx != ID_NA_INDEX)? the_gtsset->accessions->a[fempar_idx] : NULL; // id->index;
	Accession* Mpar = (malpar_idx != ID_NA_INDEX)? the_gtsset->accessions->a[malpar_idx] : NULL; // id->index;
	//	fprintf(stderr, "# %p %p %p \n", Acc, Fpar, Mpar);
	Pedigree* a_pedigree = construct_pedigree(Acc, Fpar, Mpar);
	//	fprintf(stderr, "## added %s %s %s  to pedigree\n", Acc->id->a, Fpar->id->a, Mpar->id->a);
	add_pedigree_to_vpedigree(pedigrees, a_pedigree);
      }
    }
    free_vstr(fields);
  } // done reading all lines
  free(line); // only needs to be freed once.
  // fprintf(stderr, "# size of Vpedigree pedigrees: %ld \n", pedigrees->size);
  return pedigrees;
} // end of read pedigrees (3 whitespace-separated column format)

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
    // fprintf(stderr, "### %s %s %s \n", acc_id, fempar_id, malpar_id);
  
    long acc_idx, fempar_idx, malpar_idx;
    if( // none of the ids are "NA" and id is present in genotypes data set
	strcmp(acc_id, "NA") != 0){
      acc_idx = index_of_id_in_vidxid(the_vidxid, acc_id);
      fempar_idx = index_of_id_in_vidxid(the_vidxid, fempar_id);
      malpar_idx = index_of_id_in_vidxid(the_vidxid, malpar_id);
      /* if(acc_idx != ID_NA_INDEX  && (fempar_idx != ID_NA_INDEX  ||  malpar_idx != ID_NA_INDEX)){ */
      /*   fprintf(stderr, "## %s %s %s %ld %ld %ld \n", acc_id, fempar_id, malpar_id, */
      /* 	    index_of_id_in_vidxid(the_vidxid, acc_id), */
      /* 	    index_of_id_in_vidxid(the_vidxid, fempar_id), */
      /* 	    index_of_id_in_vidxid(the_vidxid, malpar_id)); */
      /* } */
      if(
	 ((strcmp(fempar_id, "NA") != 0) || (strcmp(malpar_id, "NA") != 0)) &&
	 ((acc_idx = index_of_id_in_vidxid(the_vidxid, acc_id)) != ID_NA_INDEX) &&
	 (((fempar_idx = index_of_id_in_vidxid(the_vidxid, fempar_id)) != ID_NA_INDEX) ||
	  ((malpar_idx = index_of_id_in_vidxid(the_vidxid, malpar_id)) != ID_NA_INDEX))
	 ){
     
	Accession* Acc = the_gtsset->accessions->a[acc_idx]; // id->index;
	Accession* Fpar = (fempar_idx != ID_NA_INDEX)? the_gtsset->accessions->a[fempar_idx] : NULL; // id->index;
	Accession* Mpar = (malpar_idx != ID_NA_INDEX)? the_gtsset->accessions->a[malpar_idx] : NULL; // id->index;
	//	fprintf(stderr, "# %p %p %p \n", Acc, Fpar, Mpar);
	Pedigree* a_pedigree = construct_pedigree(Acc, Fpar, Mpar);
	//	fprintf(stderr, "## added %s %s %s  to pedigree\n", Acc->id->a, Fpar->id->a, Mpar->id->a);
	add_pedigree_to_vpedigree(pedigrees, a_pedigree);
      }
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

Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_z, double max_ok_d){
  long n_parents = parent_idxs->size;

  char* acc_id = the_pedigree->A->id->a; //Accession->id;
  long acc_idx = the_pedigree->A->index; // Accession->index;
  char* acc_gts = the_pedigree->A->genotypes->a; // the_gtsset->genotype_sets->a[the_pedigree->Accession->index];
  Accession* fparent = the_pedigree->F;
  long fparent_idx = (fparent != NULL)? fparent->index : ID_NA_INDEX; //Fparent->index;
  char* fparent_gts = (fparent != NULL)? fparent->genotypes->a : ""; // the_gtsset->genotype_sets->a[fparent_idx];
  Accession* mparent = the_pedigree->M;
  long mparent_idx = (mparent != NULL)? mparent->index : ID_NA_INDEX; //  Mparent->index;
  char* mparent_gts = (mparent != NULL)? mparent->genotypes->a : ""; // the_gtsset->genotype_sets->a[mparent_idx];
  assert(fparent != NULL  ||  mparent != NULL);

  // get best candidate parents on basis of hgmr (plus those in pedigree)
  Vlong* best_parent_candidate_idxs = construct_vlong(10); 
  if(fparent != NULL) add_long_to_vlong(best_parent_candidate_idxs, fparent_idx); // add female parent (from pedigree)
  if(mparent != NULL  &&  mparent_idx != fparent_idx) add_long_to_vlong(best_parent_candidate_idxs, mparent_idx); // add male parent (from pedigree) if distinct

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
      if(the_hgmr >= max_ok_hgmr) break; // all the rest are worse, so skip them.
      //   if(the_hgmr >= 0){
      if(the_idx != fparent_idx  &&  the_idx != mparent_idx){
	add_long_to_vlong(best_parent_candidate_idxs, the_idx);
      }     
      // }
    }
  }
  free(the_idxhgmrs);

  long ub = long_min(best_parent_candidate_idxs->size, 10000); // set the number of possible parents to consider.
  // fprintf(stderr, "XXX: %8.4lf %8.4lf  %ld \n", max_ok_hgmr, max_ok_d1, ub); 
  Vpedigree* alt_pedigrees = construct_vpedigree(10);
  for(long i=0; i<ub; i++){
    long idx1 = best_parent_candidate_idxs->a[i];
    Accession* acc1 = the_gtsset->accessions->a[idx1];
    char* id1 = acc1->id->a; // accession_ids->a[idx1];
    char* gts1 = acc1->genotypes->a; // genotype_sets->a[idx1];
    for(long j=i; j<ub; j++){
      long idx2 = best_parent_candidate_idxs->a[j];
      Accession* acc2 = the_gtsset->accessions->a[idx2];
      char* id2 = acc2->id->a; // _ids->a[idx2];
      if(! ((idx1 == fparent_idx && idx2 == mparent_idx) || (idx1 == mparent_idx && idx2 == fparent_idx))){  
	char* gts2 = acc2->genotypes->a; // the_gtsset->genotype_sets->a[idx2];
	Pedigree* alt_pedigree = construct_pedigree(the_pedigree->A, acc1, acc2); // arbitrarily put acc1 as Female parent, acc2 as male
	Pedigree_stats* alt_pedigree_stats = triple_counts(gts1, gts2, acc_gts, the_gtsset->ploidy);
	four_longs fls = q_and_d_n22x_diploid(acc1, acc2, the_pedigree->A);
	//if(get_hgmr1(alt_pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree_stats) <= 0.05  &&

	/* */
	if(n_over_d(alt_pedigree_stats->z) <= max_ok_z
	   &&
	   n_over_d(alt_pedigree_stats->d) <= max_ok_d){
	  alt_pedigree->pedigree_stats = alt_pedigree_stats;
	  add_pedigree_to_vpedigree(alt_pedigrees, alt_pedigree);
	}  /* */
      }
    }
  }
  free_vlong(best_parent_candidate_idxs);
  return alt_pedigrees;
}

Vlong* alternative_parents(Accession* the_acc, const GenotypesSet* const the_gtsset, double max_ok_hgmr){
  Vlong* best_parent_candidate_idxs = construct_vlong(10);
  long n_accessions = the_gtsset->accessions->size;
  // sort accession indices by hgmr   
  Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_accessions*sizeof(Idxhgmr));
  for(long i=0; i<n_accessions; i++){
    if(i == the_acc->index) continue;
    //long idx = parent_idxs->a[i];
    char* pgts = the_gtsset->accessions->a[i]->genotypes->a; // genotype_sets->a[idx];
    the_idxhgmrs[i].idx = i;
    //  fprintf(stderr, "# max_ok_hgmr: %6.5f\n", max_ok_hgmr);
    if(1){
      four_longs x = quick_hgmr_R(the_gtsset->accessions->a[i], the_acc, (char)(the_gtsset->ploidy+48));
      the_idxhgmrs[i].hgmr = (x.l2 > 0)? (double)x.l1/x.l2 : 2;
      if(the_idxhgmrs[i].hgmr <= max_ok_hgmr){
	fprintf(stderr, "%s  %s  %ld %ld %ld %6.5f\n",
		the_acc->id->a, the_gtsset->accessions->a[i]->id->a,
		i, x.l1, x.l2, the_idxhgmrs[i].hgmr);
      }
    }else{
      four_longs hR =  hgmr_R(pgts, the_acc->genotypes->a, (char)(the_gtsset->ploidy+48));
      the_idxhgmrs[i].hgmr = (hR.l2 > 0)? (double)hR.l1/hR.l2 : 2; //hgmr_nd(the_acc->genotypes->a, pgts);
      /* fprintf(stderr, "%s  %s  %ld %ld %ld %6.5f\n", */
      /* 	      the_acc->id->a, the_gtsset->accessions->a[i]->id->a, */
      /* 	      i, hnd.n, hnd.d, the_idxhgmrs[i].hgmr); */
      // fprintf(stderr, "%s %ld %ld %ld \n", the_acc->id, i, x.l1, x.l2);
    }
  }
  sort_idxhgmr_by_hgmr(n_accessions, the_idxhgmrs);
  //fprintf(stderr, "n_parents: %ld \n", n_parents);
  for(long i=0; i<n_accessions; i++){ // store 'good' hgmr indices 
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != the_acc->index){ 
      double the_hgmr = the_idxhgmrs[i].hgmr;
      // fprintf(stderr, "%s %ld  %8.5f\n", the_acc->id->a, the_idx, the_hgmr);
      if(the_hgmr >= max_ok_hgmr) break; // all the rest are worse, so skip them.
      //   if(the_hgmr >= 0){
      add_long_to_vlong(best_parent_candidate_idxs, the_idx);
          
      // }
    }
  }
  free(the_idxhgmrs);
  return best_parent_candidate_idxs;
} // end of alternative_parents

Vpedigree* alternative_pedigrees(Accession* the_acc, const GenotypesSet* the_gtsset, Vlong* best_parent_candidate_idxs, long ub, double max_ok_z){
  ub = long_min(best_parent_candidate_idxs->size, ub); // set the number of possible parents to consider.
  // fprintf(stderr, "XXX: %8.4lf %8.4lf  %ld \n", max_ok_hgmr, max_ok_d1, ub); 
  Vpedigree* alt_pedigrees = construct_vpedigree(10);
  for(long i=0; i<ub; i++){
    long idx1 = best_parent_candidate_idxs->a[i];
    Accession* acc1 = the_gtsset->accessions->a[idx1];
    char* id1 = acc1->id->a; // accession_ids->a[idx1];
    char* gts1 = acc1->genotypes->a; // genotype_sets->a[idx1];
    for(long j=i; j<ub; j++){
      long idx2 = best_parent_candidate_idxs->a[j];
      Accession* acc2 = the_gtsset->accessions->a[idx2];
      char* id2 = acc2->id->a; // _ids->a[idx2];
      //	if(! ((idx1 == fparent_idx && idx2 == mparent_idx) || (idx1 == mparent_idx && idx2 == fparent_idx))){  
      char* gts2 = acc2->genotypes->a; // the_gtsset->genotype_sets->a[idx2];
      Pedigree* alt_pedigree = construct_pedigree(the_acc, acc1, acc2); // arbitrarily put acc1 as Female parent, acc2 as male
      Pedigree_stats* alt_pedigree_stats = triple_counts(gts1, gts2, the_acc->genotypes->a, the_gtsset->ploidy);
      //if(get_hgmr1(alt_pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree_stats) <= 0.05  &&
      if(n_over_d(alt_pedigree_stats->z) <= max_ok_z){   
	alt_pedigree->pedigree_stats = alt_pedigree_stats;
	add_pedigree_to_vpedigree(alt_pedigrees, alt_pedigree);
      }
      //	}
    }
  }
  free_vlong(best_parent_candidate_idxs);
  return alt_pedigrees;
} // end of alternative_pedigrees


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
    // fprintf(stderr, "%ld %ld %s\n", i, idx, id);
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
