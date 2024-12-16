#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gtset.h"
#include <stdbool.h>
#include "pedigree.h"

// extern int do_checks; // option -c sets this to 1 to do some checks.

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent){ 
  Pedigree* the_pedigree = (Pedigree*)malloc(sizeof(Pedigree));
  the_pedigree->F = Fparent;
  the_pedigree->M = Mparent;
  the_pedigree->A = Acc;
  the_pedigree->pedigree_stats = construct_pedigree_stats(); // NULL;
  return the_pedigree;
}

Xcounts_3 count_crossovers(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  if(offspring == NULL){
    fprintf(stderr, "# in count_crossovers offspring Accession* is NULL. Bye.\n");
    exit(EXIT_FAILURE);
  }
  Xcounts_3 X3;
    if(Fparent == NULL){ // only have male parent in pedigree
      X3 = (Xcounts_3){(Xcounts_2){0, 0, 0}, count_crossovers_one_parent(the_gtsset, Mparent, offspring), 0, 0, 0, 0};
    }else if(Mparent == NULL){ // only have female parent in pedigree
      X3 = (Xcounts_3){count_crossovers_one_parent(the_gtsset, Fparent, offspring), (Xcounts_2){0, 0, 0}, 0, 0, 0, 0};
    }else{ // both F and M are non-NULL
      X3 = count_crossovers_two_parents(the_gtsset, Fparent, Mparent, offspring);   
    }
  return X3;
}

Xcounts_2 count_crossovers_one_parent(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  if(parent == NULL  ||  offspring == NULL) {
    return (Xcounts_2){-1, -1, -1};
  }
  long Xmin = 0, Xmax = 0, Nhet = 0; // Nhet = number of heterozyg gts in parent
  long Xa = 0, Xb = 0;
  long prev_chrom_number = -1, prev_phase_a = -1, prev_phase_b = -1;
  long phase_a = -1, phase_b = -1, chrom_number;
 
  for(long i=0; i < parent->genotypes->length; i++){

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    char o_phase = offspring->phases->a[i];

       chrom_number = the_gtsset->chromosomes->a[i];
    if(chrom_number != prev_chrom_number){ // now on next chromosome
      if(Xa < Xb){ // add chromosome Xmin, Xmax to totals
	Xmin += Xa; Xmax += Xb;
      }else{
	Xmin += Xb; Xmax += Xa;
      }
      // reset for new chromosome:
      Xa = 0; Xb = 0;
      prev_chrom_number = chrom_number;
      phase_a = -1;
      phase_b = -1;
      prev_phase_a = -1; // needed
      prev_phase_b = -1;
    }

    // ########################################
    char p_gt = parent->genotypes->a[i];
    char p_phase = parent->phases->a[i];

    if(p_gt == '1'){
      Nhet++; // counts the number of markers which are heterozyg in the parent, and non-missing in the offspring

      two_longs phases_ab = get_1marker_phases_wrt_1parent(p_phase, o_gt, o_phase);
      phase_a = phases_ab.l1;
      phase_b = phases_ab.l2;
    }
    // ######################################

    // compare current phases with previous values,
    // update crossover counts, and  update prev_phase_a, prev_phase_b
    if(prev_phase_a >= 0  &&  phase_a != prev_phase_a) Xa++; // phase has changed - crossover
    prev_phase_a = phase_a;
    if(prev_phase_b >= 0  &&  phase_b != prev_phase_b) Xb++; // phase has changed - crossover
    prev_phase_b = phase_b;
    
  } // end loop over markers
  if(Xa < Xb){ // add the crossovers from the last chromosome.
    Xmin += Xa; Xmax += Xb;
  }else{
    Xmin += Xb; Xmax += Xa;
  }
 
  return (Xcounts_2){Xmin, Xmax, Nhet};
} // end of count_crossovers

Xcounts_2 count_crossovers_one_chromosome(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring, long first, long next){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  if(parent == NULL  ||  offspring == NULL) {
    return (Xcounts_2){-1, -1, -1};
  }
 
  long Xa = 0, Xb = 0, Nhet = 0;
  long prev_chrom_number = -1, chrom_number;
  long prev_phase_a = -1, prev_phase_b = -1, phase_a = -1, phase_b = -1;
 
  for(long i=first; i < next; i++){

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    char o_phase = offspring->phases->a[i];

    // ########################################
    char p_gt = parent->genotypes->a[i];
    char p_phase = parent->phases->a[i];

    if(p_gt == '1'){
      Nhet++; // counts the number of markers which are heterozyg in the parent, and non-missing in the offspring
      two_longs phases_ab = get_1marker_phases_wrt_1parent(p_phase, o_gt, o_phase);
      phase_a = phases_ab.l1;
      phase_b = phases_ab.l2;
    }
    // ######################################

    // compare current phases with previous values, update crossover counts, 
    // and  update prev_phase_a, prev_phase_b
    if(prev_phase_a >= 0  &&  phase_a != prev_phase_a) Xa++; // phase has changed - crossover
    prev_phase_a = phase_a;
    if(prev_phase_b >= 0  &&  phase_b != prev_phase_b) Xb++; // phase has changed - crossover
    prev_phase_b = phase_b;
    
  } // end loop over markers
  return (Xcounts_2){Xa, Xb, Nhet};
} // end of count_crossovers_one_chromosome

Xcounts_3 count_crossovers_two_parents(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  long NhetF = 0, XFmin_2 = 0, XFmax_2 = 0, XFmin_3 = 0, XFmax_3 = 0;
  long NhetM = 0, XMmin_2 = 0, XMmax_2 = 0, XMmin_3 = 0, XMmax_3 = 0;
 
  long n_chroms = the_gtsset->chromosome_start_indices->size - 1;
  for(long i=0; i < n_chroms; i++){
   
    long start_index = the_gtsset->chromosome_start_indices->a[i];
    long next_start_index = the_gtsset->chromosome_start_indices->a[i+1];
    
    Xcounts_2 FX = count_crossovers_one_chromosome(the_gtsset, Fparent, offspring, start_index, next_start_index);
   
    NhetF += FX.Nhet;
    if(FX.Xa < FX.Xb){
      XFmin_2 += FX.Xa; XFmax_2 += FX.Xb;
    }else{
      XFmin_2 += FX.Xb; XFmax_2 += FX.Xa;
    }

    Xcounts_2 MX = count_crossovers_one_chromosome(the_gtsset, Mparent, offspring, start_index, next_start_index);
    NhetM += MX.Nhet;
    if(MX.Xa < MX.Xb){
      XMmin_2 += MX.Xa; XMmax_2 += MX.Xb;
    }else{
      XMmin_2 += MX.Xb; XMmax_2 += MX.Xa;
    }

    long X_Fa_Mb = FX.Xa + MX.Xb; // crossovers if F is parent of a, M is parent of b
    long X_Fb_Ma = FX.Xb + MX.Xa; // crossovers if F is parent of b, M is parent of a
    if(X_Fa_Mb < X_Fb_Ma){
      XFmin_3 += FX.Xa;
      XFmax_3 += FX.Xb;
      XMmin_3 += MX.Xb;
      XMmax_3 += MX.Xa;
    }else{
      XFmin_3 += FX.Xb;
      XFmax_3 += FX.Xa;
      XMmin_3 += MX.Xa;
      XMmax_3 += MX.Xb;
    }
  } // end loop over chromosomes
  Xcounts_3 X3 = {(Xcounts_2){XFmin_2, XFmax_2, NhetF}, (Xcounts_2){XMmin_2, XMmax_2, NhetM}, XFmin_3, XFmax_3, XMmin_3, XMmax_3};
  return X3;
}

two_longs get_1marker_phases_wrt_1parent(char p_phase, char o_gt, char o_phase){
  long phase_a, phase_b;
  long offA = (o_gt == '0'  ||  (o_gt == '1'  &&  o_phase == 'p'))? 0 : 1; // 0 <-> ref allele, 1<->alt.
  long offB = (o_gt == '0'  ||  (o_gt == '1'  &&  o_phase == 'm'))? 0 : 1; // 0 <-> ref allele, 1<->alt.

  if(p_phase == 'p'){ // i.e. 0|1 in vcf    
    phase_a = (offA == 0)? 0 : 1;
    phase_b = (offB == 0)? 0 : 1; 
  }else if(p_phase == 'm'){ // i.e. 1|0 in vcf
    phase_a = (offA == 0)? 1 : 0;
    phase_b = (offB == 0)? 1 : 0;	
  }else{
    fprintf(stderr, "p_phase is neither p nor m: %c\n", p_phase);
    exit(0);
  }
  return (two_longs){phase_a, phase_b};
}

Pedigree_stats* triple_counts(char* gts1, char* gts2, char* proggts, long ploidy){ // 

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
  Pedigree_stats* pedigree_stats = construct_pedigree_stats(); // (Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  long n_0 = // 15 triples consistent with true parents-offspring relationship,
    // with no genotyping errors
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 +
    n_11_0 + n_11_1 + n_11_2 + // include these in denom of d or not?
    n_12_1 + n_12_2 +
    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 = // 10 triples requiring one (0<->1 or 1<->2) error for consistency 
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // these 10 can happen if just one error of 0<->1 or 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // these 2 can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.
 
  long n_11_x =  n_11_0 + n_11_1 + n_11_2;
  // ************************************
  long n_a2 = n_00_2 + n_22_0; // delta = 2 triples
  long n_a1 =   n_00_1 + n_02_0 + n_02_2 +
    n_20_0 + n_20_2 + n_22_1; // delta = 1 triples used by apparent
  long n_o1 =   n_01_2 + n_10_2 + n_12_0 + n_21_0; // other delta = 1 triples.

  long n_00 = n_00_0 + n_00_1 + n_00_2;
  long n_22 = n_22_0 + n_22_1 + n_22_2;

  long n_02 = n_02_0 + n_02_1 + n_02_2;
  long n_20 = n_20_0 + n_20_1 + n_20_2;

  long n_01 = n_01_0 + n_01_1 + n_01_2;
  long n_10 = n_10_0 + n_10_1 + n_10_2;
	
  long n_12 = n_12_0 + n_12_1 + n_12_2;
  long n_21 = n_21_0 + n_21_1 + n_21_2;

  // ************************************
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0; // par2 gt can be missing
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1; // par2 gt can be missing
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2; // par2 gt can be missing
  
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0; // par2 gt can be missing
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1; // par2 gt can be missing
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2; // par2 gt can be missing

  long n_3x_0 = n_30_0 + n_31_0 + n_32_0; // + n_03_0;
  long n_3x_1 = n_30_1 + n_31_1 + n_32_1; // + n_03_1;
  long n_3x_2 = n_30_2 + n_31_2 + n_32_2; // + n_03_2;

  long z_numer = n_00_1 + n_22_1;
  long z_denom = z_numer + n_00_0 + n_22_2 + n_00_2 + n_22_0; // = n_00_x + n_22_x

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0; // par1 gt can be missing
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1; // par1 gt can be missing
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2; // par1 gt can be missing
  
  long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0; // par1 gt can be missing
  long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1; // par1 gt can be missing
  long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2; // par1 gt can be missing

  long hgmr1_numer = n_0x_2 + n_2x_0;
  long hgmr2_numer = n_x0_2 + n_x2_0;
  long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer;
  long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer;

  long r1_numer = n_0x_1 + n_2x_1; // + n_0x_2 + n_2x_1 + n_2x_0;
  long r1_denom = r1_numer + n_0x_0 + n_2x_2;
  long r2_numer = n_x0_1 + n_x2_1; // n_x0_2 + n_x2_1 + n_x2_0;
  long r2_denom = r2_numer + n_x0_0 + n_x2_2;
 
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
  
  pedigree_stats->agmr12 = (ND) {agmr12_numer, agmr12_denom};

  pedigree_stats->par1_hgmr = (ND) {hgmr1_numer, hgmr1_denom};
  pedigree_stats->par1_R = (ND) {r1_numer, r1_denom};
 
  pedigree_stats->par2_hgmr = (ND) {hgmr2_numer, hgmr2_denom};
  pedigree_stats->par2_R = (ND)  {r2_numer, r2_denom};

  pedigree_stats->d = (ND){n_1 + n_2, n_0 + n_1 + n_2};

  pedigree_stats->z = (ND) {z_numer, z_denom};

  pedigree_stats->all_good_count = n_0 + n_1 + n_2;

  return pedigree_stats;
} // end of triple_counts

Pedigree_stats* bitwise_triple_counts(Accession* par1, Accession* par2, Accession* prog){
  
  long n_00_1_22_1 = 0, n_total_no_md = 0, n_00_22 = 0, n_total_x_11 = 0;
  long n_00_2_22_0 = 0, n_00_2 = 0, n_22_0 = 0;
  long n_0xorx0_2_nomd = 0; // prog is 2 and at least one parent is 0.
  long n_2xorx2_0_nomd = 0; // prog is 0 and at least one parent is 2.
  long n_01or10_1 = 0; // parents are 0 and 1, prog is 1;
  long n_0x_1_2x_1 = 0, n_x0_1_x2_1 = 0;
  long n_0x_0_2x_2 = 0, n_x0_0_x2_2 = 0;
  long ndiff12 = 0, ndiff01 = 0, ndiff02 = 0;
  long hgmr1_numerator = 0, hgmr1_denominator = 0;
  long hgmr2_numerator = 0, hgmr2_denominator = 0;
  unsigned long long i0, j0, k0,  i1, j1, k1,  i2, j2, k2, missing;

  // A,B->dosage 0,0->0, 0,1->1, 1,0->missing, 1,1->2 

  for(long i_long = 0; i_long < par1->Abits->size; i_long++){
    unsigned long long iA = par1->Abits->a[i_long]; // i : parent 1
    unsigned long long iB = par1->Bbits->a[i_long];
    unsigned long long jA = par2->Abits->a[i_long]; // j : parent 2
    unsigned long long jB = par2->Bbits->a[i_long];
    unsigned long long kA = prog->Abits->a[i_long]; // k : progeny
    unsigned long long kB = prog->Bbits->a[i_long];
    
    missing = (iA & ~iB) | (jA & ~jB) | (kA & ~kB); // one of the three is missing.
    //  i02 = ~(iA ^ iB); // is i homozyg?
    //  j02 = ~(jA ^ jB); // is j homozyg?
    i0 = ~(iA | iB); // is i ref homozyg?
    j0 = ~(jA | jB); // is j ref homozyg?
    k0 = ~(kA | kB); // is k ref homozyg?
    i1 = ~iA & iB; // is i heterozyg?
    j1 = ~jA & jB; // is j heterozyg?
    k1 = ~kA & kB; // is k heterozyg?
    i2 = iA & iB; // is i alt homozyg?
    j2 = jA & jB; // is j alt homozyg?
    k2 = kA & kB; // is k alt homozyg?

    //  unsigned long long is2_0x = k2 & (i0 | j0) & ~missing; // 2_0x, 2_x0
    //  unsigned long long is0_2x = k0 & (i2 | j2) & ~missing; // 0_2x, 0_x2
    
    unsigned long long is2_0x_or_0_2x = (k2 & i0) | (k0 & i2); // for hgmr1 numerator (j can be missing)
    unsigned long long is2_x0_or_0_x2 = (k2 & j0) | (k0 & j2); // for hgmr2 numerator (i can be missing)
    unsigned long long ik_not_1 = (i0|i2) & (k0|k2); // neither i nor k is 1, i.e. 00 | 22 | 02 | 20  for hgmr1 denom (j can be missing)
    unsigned long long jk_not_1 = (j0|j2) & (k0|k2); // neither j nor k is 1, i.e. 00 | 22 | 02 | 20  for hgmr2 denom (i can be missing)

    unsigned long long is2_00 = (i0 & j0 & k2);
    unsigned long long is0_22 = (i2 & j2 & k0);
   
    unsigned long long is2_0xorx0_nomd = k2 & (i0 | j0) & ~missing; // 2_0x, 2_x0
    unsigned long long is0_2xorx2_nomd = k0 & (i2 | j2) & ~missing; // 0_2x, 0_x2
    unsigned long long is00_22 = ((i0 & j0) | (i2 & j2)) & ~missing; // k is 0, 1, or 2
   
    //   unsigned long long is00_2_22_0 = (i0 & j0 & k2) | (i2 & j2 & k0) & ~missing;
    unsigned long long is1_00_1_22 = k1 & is00_22 & ~missing;
    unsigned long long isx_11 = i1 & j1 & ~missing; //
    
    unsigned long long is_i0or2_k1 = (i0 | i2) & k1; // & ~missing;    
    unsigned long long is0x_0_or_2x_2 = (i0 & k0) | (i2 & k2);
     
    unsigned long long is_j0or2_k1 = (j0 | j2) & k1; // & ~missing;
    unsigned long long isx0_0_or_x2_2 = (j0 & k0) | (j2 & k2);
     
    /* unsigned long long isx_01or10 = ((i0 & j1) | (i1 & j0)) &  ~missing; */
    /* unsigned long long is1_01or10 = k1 & isx_01or10; */
    /* unsigned long long is1_diag = k1 & (~(iA^jA)) & (~(iB^jB)) & ~missing; //prog is 1, parents equal */
    /* unsigned long long is1_xx = k1 & ~missing; */
    // fprintf(stderr, "%llu  %llu\n", isx_01or10, is1_01or10);

    unsigned long long i_ne_j = ((iA ^ jA) | (iB ^ jB)) & ~missing;
    unsigned long long i_ne_k = ((iA ^ kA) | (iB ^ kB)) & ~missing;
    unsigned long long j_ne_k = ((jA ^ kA) | (jB ^ kB)) & ~missing;
  
    n_0xorx0_2_nomd += __builtin_popcountll(is2_0xorx0_nomd);
    n_2xorx2_0_nomd += __builtin_popcountll(is0_2xorx2_nomd);
    n_00_1_22_1 +=  __builtin_popcountll(is1_00_1_22); // either 00_1 or 22_1
    n_00_2 += __builtin_popcountll(is2_00);
    n_22_0 += __builtin_popcountll(is0_22);
    n_00_2_22_0 = n_00_2 + n_22_0;

    n_00_22 += __builtin_popcountll(is00_22);
    
    n_0x_1_2x_1 += __builtin_popcountll(is_i0or2_k1); // R1 numerator
    n_0x_0_2x_2 += __builtin_popcountll(is0x_0_or_2x_2); // R1 denominator is this + R1 numerator
    n_x0_1_x2_1 += __builtin_popcountll(is_j0or2_k1); // R2 numerator
    n_x0_0_x2_2 += __builtin_popcountll(isx0_0_or_x2_2); // R2 denominator is this + R2 numerator
    //   n_01or10_1 += __builtin_popcountll(is1_01or10);

    hgmr1_numerator += __builtin_popcountll(is2_0x_or_0_2x);
    hgmr1_denominator += __builtin_popcountll(ik_not_1);
    hgmr2_numerator += __builtin_popcountll(is2_x0_or_0_x2);
    hgmr2_denominator += __builtin_popcountll(jk_not_1);
  
    n_total_no_md += 64 - __builtin_popcountll(missing);
    n_total_x_11 += __builtin_popcountll(isx_11); // 

    ndiff12 += __builtin_popcountll(i_ne_j);
    ndiff01 += __builtin_popcountll(i_ne_k);
    ndiff02 += __builtin_popcountll(j_ne_k);

  } // end of loop over 64 bit chunks
  
  //fprintf(stderr, "ZZZAAA: %ld %ld %ld %ld \n", n_0x_1_2x_1, n_0x_0_2x_2, n_x0_1_x2_1, n_x0_0_x2_2);
	
  Pedigree_stats* pedigree_stats = construct_pedigree_stats(); //(Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  pedigree_stats->agmr12 = (ND) {ndiff12, n_total_no_md}; // agmr between parents 
  pedigree_stats->agmr01 = (ND) {ndiff01, n_total_no_md}; 
  pedigree_stats->agmr02 = (ND) {ndiff02, n_total_no_md};

  long F = 1;
  pedigree_stats->d = (ND) {n_0xorx0_2_nomd + n_2xorx2_0_nomd
			    + F*n_00_1_22_1
			    , n_total_no_md};
  pedigree_stats->z = (ND) {n_00_1_22_1, n_00_22};
  pedigree_stats->par1_hgmr = (ND) {hgmr1_numerator, hgmr1_denominator};
  pedigree_stats->par2_hgmr = (ND) {hgmr2_numerator, hgmr2_denominator};
  pedigree_stats->par1_R = (ND) {n_0x_1_2x_1, n_0x_1_2x_1 + n_0x_0_2x_2};
  pedigree_stats->par2_R = (ND) {n_x0_1_x2_1, n_x0_1_x2_1 + n_x0_0_x2_2};
  pedigree_stats->n_01or10_1 = n_01or10_1;
  pedigree_stats->all_good_count = n_total_no_md;
  //assert(pedigree_stats->par1_R.n == Rnd.n);
  //assert(pedigree_stats->par1_R.d == Rnd.d);
  return pedigree_stats;
} // end of bitwise_triple_counts

Vpedigree*  calculate_triples_for_one_accession(Accession* prog, GenotypesSet* the_genotypes_set, Viaxh* cppps, long max_candidate_parents){

  // sort the parent candidates and keep the best ones
  long limited_n_candpairs = cppps->size;
  if(cppps->size == 0){ 
  }else if(cppps->size > max_candidate_parents){ // if too many parent candidates, just take the max_candidate_parents best ones
    sort_viaxh_by_xhgmr(cppps);
    limited_n_candpairs = max_candidate_parents; // limit candidate to max_candidate_parents
  }
 
  Vpedigree* alt_pedigrees = construct_vpedigree(1000);
  for(long ii=0; ii<limited_n_candpairs; ii++){
    long par1idx = cppps->a[ii]->idx;
    Accession* par1 = the_genotypes_set->accessions->a[par1idx];
    for(long jj=ii; jj<limited_n_candpairs; jj++){	
      long par2idx = cppps->a[jj]->idx;
      Accession* par2 = the_genotypes_set->accessions->a[par2idx];
      Pedigree_stats* the_ps;
      Pedigree* the_pedigree = construct_pedigree(prog, par1, par2);
    
      if(0){
	the_ps = bitwise_triple_counts(par1, par2, prog);
	the_ps->xhgmr1 = cppps->a[ii]->xhgmr;
	the_ps->xhgmr2 = cppps->a[jj]->xhgmr;
	//	the_ps->hgmr1 = cppps->a[ii]->hgmr;
	//	the_ps->hgmr2 = cppps->a[jj]->hgmr;
      }else{
	the_ps = calculate_pedigree_stats(the_pedigree, the_genotypes_set);
      }
      the_pedigree->pedigree_stats = the_ps;
      push_to_vpedigree(alt_pedigrees, the_pedigree);
    
    } // end loop over parent 2
  } // end loop over parent 1
  return alt_pedigrees;
} //

Pedigree_stats* construct_pedigree_stats(void){
  Pedigree_stats* the_ps = (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
  the_ps->agmr12 = (ND) {0, 0};
  the_ps->agmr01 = (ND) {0, 0};
  the_ps->agmr02 = (ND) {0, 0};
    
  the_ps->par1_hgmr = (ND) {0, 0};
  the_ps->par1_xhgmr = (ND) {0, 0}; 
  the_ps->par1_R = (ND) {0, 0};
  
  the_ps->par2_hgmr = (ND) {0, 0}; 
  the_ps->par2_xhgmr = (ND) {0, 0};
  the_ps->par2_R = (ND) {0, 0};

  the_ps->d = (ND) {0, 0};
  the_ps->z = (ND) {0, 0};

  the_ps->all_good_count = 0;

  the_ps->hgmr1_n = NAN;
  the_ps->hgmr2_n = NAN;
  the_ps->R1_n = NAN;
  the_ps->R2_n = NAN;
  the_ps->d_n = NAN;
  the_ps->z_n = NAN;
  //the_ps->scaled_d = NAN;
  //the_ps->max_scaleddz = NAN;
  the_ps->xhgmr1 = NAN;
  the_ps->xhgmr2 = NAN;
  return the_ps;
}
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset){ //, long* d0counts, long* d1counts, long* d2counts){ //, GenotypesSet* the_gtsset){
  long ploidy = the_gtsset->ploidy;
  double d_scale_factor = the_gtsset->d_scale_factor;
  Pedigree_stats* the_ps; //  = construct_pedigree_stats(); // (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
  assert(the_pedigree->F != NULL  ||  the_pedigree->M != NULL); // shouldn't have both parents NULL //
  if(the_pedigree->F != NULL  &&  the_pedigree->M != NULL){
   
    the_ps = bitwise_triple_counts(the_pedigree->F,  the_pedigree->M,  the_pedigree->A);

    if(0){ // compare bitwise, nonbitwise calculations as check - slow
      Pedigree_stats* nobw_ps = triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a, ploidy );
      assert
	(NDs_equal(the_ps->agmr12, nobw_ps->agmr12));
      assert(NDs_equal(the_ps->par1_hgmr, nobw_ps->par1_hgmr));
      assert(NDs_equal(the_ps->par1_R, nobw_ps->par1_R));
      assert(NDs_equal(the_ps->par2_hgmr, nobw_ps->par2_hgmr));
      assert(NDs_equal(the_ps->par2_R, nobw_ps->par2_R));
      assert(NDs_equal(the_ps->d, nobw_ps->d));
      assert(NDs_equal(the_ps->z, nobw_ps->z));
    }
       
     the_ps->par1_xhgmr = xhgmr(the_gtsset, the_pedigree->F, the_pedigree->A, false);
     the_ps->par2_xhgmr = xhgmr(the_gtsset, the_pedigree->M, the_pedigree->A, false);
  }else{ // one of the parents is NULL
    the_ps = construct_pedigree_stats();
    the_ps->agmr12 = (ND) {0, 0};
    the_ps->z = (ND) {0, 0};
    // the_ps->xz = (ND) {0, 0};
    if(the_pedigree->F != NULL){ // we have female parent id, no male parent id
      four_longs hgmrR = hgmr_R(the_pedigree->F->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
   
      the_ps->par1_hgmr.n = hgmrR.l1;
      the_ps->par1_hgmr.d = hgmrR.l2;
      the_ps->par1_R.n = hgmrR.l3;
      the_ps->par1_R.d = hgmrR.l4;
      the_ps->par2_hgmr = (ND) {0, 0};
      the_ps->par2_R = (ND) {0, 0};
        the_ps->par1_xhgmr = xhgmr(the_gtsset, the_pedigree->F, the_pedigree->A, false);
        the_ps->par2_xhgmr = (ND) {0, 0};
    }else{ // we have male parent id, no female parent id
      if(DO_ASSERT) assert(the_pedigree->M != NULL);
      //        fprintf(stderr, "pedigree with male parent only.\n");
      the_ps->par1_hgmr = (ND) {0, 0};
      the_ps->par1_R = (ND) {0, 0};
      four_longs hgmrR = hgmr_R(the_pedigree->M->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
      the_ps->par2_hgmr.n = hgmrR.l1;
      the_ps->par2_hgmr.d = hgmrR.l2;
      the_ps->par2_R.n = hgmrR.l3;
      the_ps->par2_R.d = hgmrR.l4;
        the_ps->par1_xhgmr = (ND) {0, 0};
        the_ps->par2_xhgmr = xhgmr(the_gtsset, the_pedigree->M, the_pedigree->A, false);  
    }   
  }
  the_ps->hgmr1_n = n_over_d(the_ps->par1_hgmr)/the_gtsset->mean_hgmr;
  the_ps->hgmr2_n = n_over_d(the_ps->par2_hgmr)/the_gtsset->mean_hgmr;
  the_ps->R1_n = n_over_d(the_ps->par1_R)/the_gtsset->mean_R;
  the_ps->R2_n = n_over_d(the_ps->par2_R)/the_gtsset->mean_R;
  the_ps->d_n = n_over_d(the_ps->d)/the_gtsset->mean_d;
  the_ps->z_n = n_over_d(the_ps->z)/the_gtsset->mean_z;
  
  the_ps->xhgmr1 = n_over_d(the_ps->par1_xhgmr);
  the_ps->xhgmr2 = n_over_d(the_ps->par2_xhgmr);
  return the_ps;
}

// get the indices of all the accessions which, according to the_vped, have offspring.
const Vlong* accessions_with_offspring(const Vpedigree* the_vped, long n_accessions){ // , long n_accessions){
  Vlong* offspring_counts = construct_vlong_zeroes(n_accessions);
  for(long i=0; i<the_vped->size; i++){
    const Pedigree* the_ped = the_vped->a[i];
    if(the_ped->F != NULL) {
   
      offspring_counts->a[the_ped->F->index]++;
    }
    if(the_ped->M != NULL) {
      offspring_counts->a[the_ped->M->index]++;
    }
  }

  Vlong* accidxs_with_offspring = construct_vlong(100);
  for(long i=0; i<offspring_counts->size; i++){
    if(offspring_counts->a[i] > 0){
      push_to_vlong(accidxs_with_offspring, i);
    }
  }
  free_vlong(offspring_counts);
  return accidxs_with_offspring;
}


void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats, bool verbose){
  if(verbose){ // for each of 9 quantities of interest, print denominator and ratio (numerator/denominator)
    print_d_r(fh, the_pedigree_stats->agmr12);
    print_d_r(fh, the_pedigree_stats->par1_hgmr);
    // print_d_r(fh, the_pedigree_stats->par1_xhgmr);
    print_d_r(fh, the_pedigree_stats->par1_R);
    print_d_r(fh, the_pedigree_stats->par2_hgmr);
    //  print_d_r(fh, the_pedigree_stats->par2_xhgmr);
    print_d_r(fh, the_pedigree_stats->par2_R);
    print_d_r(fh, the_pedigree_stats->d);
    print_d_r(fh, the_pedigree_stats->z);
    //print_d_r(fh, the_pedigree_stats->d_old);

  }else{ // print ratios but not denominators
    print_n_over_d(fh, the_pedigree_stats->agmr12);
    print_n_over_d(fh, the_pedigree_stats->par1_hgmr);
    // print_n_over_d(fh, the_pedigree_stats->par1_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par1_R);
    print_n_over_d(fh, the_pedigree_stats->par2_hgmr);
    //  print_n_over_d(fh, the_pedigree_stats->par2_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par2_R);
    print_n_over_d(fh, the_pedigree_stats->d);
    print_n_over_d(fh, the_pedigree_stats->z);
    // print_n_over_d(fh, the_pedigree_stats->d_old);
  }	     	   
}

void print_pedigree_normalized(FILE* fh, Pedigree* the_pedigree){
  //double mean_hgmr, double mean_R, double mean_d, double mean_z){
  Accession* F = the_pedigree->F;
  Accession* M = the_pedigree->M;
  fprintf(fh, "%s  %s  %ld  ", (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA",  the_pedigree->pedigree_stats->all_good_count);
  print_normalized_pedigree_stats(fh, the_pedigree->pedigree_stats);
}

void print_normalized_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){ 
  print_n_over_d(fh, the_pedigree_stats->agmr12);

  print_double_nan_as_hyphen(fh, the_pedigree_stats->hgmr1_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->R1_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->hgmr2_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->R2_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->d_n);
}

void print_double_nan_as_hyphen(FILE* fh, double x){
  if(isnan(x)){
    fprintf(fh, "-  ");
  }else{
    fprintf(fh, "%7.5lf  ", x);
  }
}

void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees, long max_to_print, bool verbose){
  fprintf(fh, " %3ld  ", alt_pedigrees->size);
  long n_to_print = (max_to_print < alt_pedigrees->size)? max_to_print : alt_pedigrees->size;
  for(long i=0; i<n_to_print; i++){
    Pedigree* alt_pedigree = alt_pedigrees->a[i];
    fprintf(fh, "%20s %20s %ld ", alt_pedigree->F->id->a, alt_pedigree->M->id->a, alt_pedigree->pedigree_stats->all_good_count);
    print_pedigree_stats(fh, alt_pedigree->pedigree_stats, verbose);
  } 
}


long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_self_r, double max_ok_d, double max_ok_z){
  //  > 0 pedigree looks good
  //  <= 0 pedigree looks bad
  //  0  pedigree bad, not self
  //  1  pedigree ok, self
  //  2  pedigree ok, distinct parents
  // -1  pedigree bad, but looks like self (female parent correct)
  // -2  pedigree bad, but looks like self (male parent correct)
 
  double agmr12 = n_over_d(p->agmr12);
  //double hgmr1 = n_over_d(p->par1_hgmr);
  double r1 = n_over_d(p->par1_R);
  //double hgmr2 = n_over_d(p->par2_hgmr);
  double r2 = n_over_d(p->par2_R);
  double d = n_over_d(p->d);
  double z = n_over_d(p->z);
  long result = 0;

  if(agmr12 <= max_self_agmr12){ // parents in pedigree are identical or very similar)
    if((r1 <= max_self_r) && (r2 <= max_self_r) && (d <= max_ok_d) && (z <= max_ok_z) ){
      result = 1;
    }
  }else{ // parents in pedigree are not very similar
    if((r1 > max_self_r) && (r2 > max_self_r) && (d <= max_ok_d) && (z <= max_ok_z) ){
      result = 2;
    }else if( (r1 <= max_self_r) && ((d > max_ok_d) || (z > max_ok_z)) ){
      result = -1;
    }else if( (r2 <= max_self_r) && ((d > max_ok_d) || (z > max_ok_z)) ){
      result = -2;
    }
  }
  return result;
}

void free_pedigree(const Pedigree* the_pedigree){
  if(the_pedigree == NULL) return;
  free(the_pedigree->pedigree_stats);
  free((Pedigree*) the_pedigree);
}

// *****  Vpedigree  *****
// read pedigree from file in 3 whitespace-separated cols format
Vpedigree* read_and_store_pedigrees_3col(FILE* p_stream, Vidxid* the_gt_vidxid, GenotypesSet* the_gtsset){ 

  char* line = NULL;
  size_t len = 0;
  ssize_t nread;
  char* saveptr = NULL;
 
  long prog_id_NA_count = 0;
  long not_in_genotypes_set_count = 0;
  long no_parent_gts_count = 0; // counts the accessions with neither parent having genotypes
  Vpedigree* pedigrees = construct_vpedigree(1000);
  while((nread = getline(&line, &len, p_stream)) != -1){
    Vstr* fields = construct_vstr(3);
    char* token = strtok_r(line, "\t \n\r", &saveptr);
 
    push_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      push_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of parent ids
    }
    // construct a Pedigree struct from first 3 fields  
    char* acc_id = ith_str_from_vstr(fields, 0); // ith_str ... copies the string, i.e. allocates more memory
    char* fempar_id = ith_str_from_vstr(fields, 1);
    char* malpar_id = ith_str_from_vstr(fields, 2);  
    long acc_idx, fempar_idx, malpar_idx;
  
    if(strcmp(acc_id, "NA") != 0){  // id of this accession is not "NA" 
      acc_idx = index_of_id_in_vidxid(the_gt_vidxid, acc_id);
      Accession* Acc = the_gtsset->accessions->a[acc_idx]; // id->index;
      if(acc_idx == ID_NA_INDEX){ // no genotypes for this accession
	not_in_genotypes_set_count++;
      }else{ // have genotypes for this accession
	fempar_idx = index_of_id_in_vidxid(the_gt_vidxid, fempar_id);
	malpar_idx = index_of_id_in_vidxid(the_gt_vidxid, malpar_id);	 
	if( (fempar_idx != ID_NA_INDEX) || (malpar_idx != ID_NA_INDEX) ){ // pedigree file gives at least 1 parent for which there are genotypes
	   
	  Acc->has_pedigree = true;
	  Acc->Fpar_idx = fempar_idx;
	  Acc->Mpar_idx = malpar_idx;
	  Accession* Fpar = (fempar_idx != ID_NA_INDEX)?
	    the_gtsset->accessions->a[fempar_idx] : NULL; // id->index;
	  Accession* Mpar = (malpar_idx != ID_NA_INDEX)?
	    the_gtsset->accessions->a[malpar_idx] : NULL; // id->index;
	  Pedigree* a_pedigree = construct_pedigree(Acc, Fpar, Mpar);
	  push_to_vpedigree(pedigrees, a_pedigree); // store pedigrees (with genotypes for accession and at least 1 of the parents)
	}else{
	  no_parent_gts_count++;
	  fprintf(stdout, "Invalid pedigree for accession: %s ; both parents are NA or lack genotypes.\n", Acc->id->a);  
	}
      }
    }else{
      prog_id_NA_count++;
    }
    free_vstr(fields);
  } // done reading all lines
  fprintf(stderr, "# pedigrees with NA for progeny id: %ld\n", prog_id_NA_count);
  fprintf(stderr, "# pedigrees with progeny not in genotypes set: %ld\n", not_in_genotypes_set_count);
  fprintf(stderr, "# pedigrees with both parent ids NA or not in genotypes set: %ld\n", no_parent_gts_count);
  fprintf(stderr, "# pedigrees with progeny and at least one parent in genotypes set: %ld\n", pedigrees->size);
  free(line); // only needs to be freed once.
  // fprintf(stderr, "# size of Vpedigree pedigrees: %ld \n", pedigrees->size);
  return pedigrees;
} // end of read pedigrees (3 whitespace-separated column format)

Vpedigree* construct_vpedigree(long cap){
  Vpedigree* the_vped = (Vpedigree*)malloc(sizeof(Vpedigree));
  the_vped->capacity = cap;
  the_vped->size = 0;
  the_vped->a = (Pedigree**)malloc(the_vped->capacity*sizeof(Pedigree*));
  return the_vped;
}
  
void push_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped){
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

void sort_vpedigree_by_d(Vpedigree* the_vped){ // sort by max(scaled_d, z) (increasing) 
  qsort(the_vped->a, the_vped->size, sizeof(Pedigree*), compare_pedigree_d);
}

int compare_pedigree_d(const void* a, const void* b){
  ND dnd1 = (*((Pedigree**)a))->pedigree_stats->d;
  ND dnd2 = (*((Pedigree**)b))->pedigree_stats->d;
  
   double d1 = (dnd1.d > 0)? (double)dnd1.n/dnd1.d : 100.0;
  double d2 = (dnd2.d > 0)? (double)dnd2.n/dnd2.d : 100.0;

  if(d1 > d2){
     return 1;
  }else if(d1 < d2){
     return -1;
  }else{
     return 0;
  }
}



Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_z, double max_ok_d){
  
  long n_parents = parent_idxs->size;

  // get the accession ids and indices in the pedigree:
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
  if(fparent != NULL) push_to_vlong(best_parent_candidate_idxs, fparent_idx); // add female parent (from pedigree)
  if(mparent != NULL  &&  mparent_idx != fparent_idx) push_to_vlong(best_parent_candidate_idxs, mparent_idx); // add male parent (from pedigree) if distinct

  // sort accession indices by hgmr   
  Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_parents*sizeof(Idxhgmr));
  for(long i=0; i<parent_idxs->size; i++){
    long idx = parent_idxs->a[i];
    char* pgts = the_gtsset->accessions->a[idx]->genotypes->a; // genotype_sets->a[idx];
    the_idxhgmrs[i].idx = idx;
    the_idxhgmrs[i].hgmr = hgmr(acc_gts, pgts); //the_hgmr; //the_idxhgmr = {idx, the_hgmr);
  }
  sort_idxhgmr_by_hgmr(n_parents, the_idxhgmrs);
  for(long i=0; i<n_parents; i++){ // store 'good' hgmr indices 
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != acc_idx){
      double the_hgmr = the_idxhgmrs[i].hgmr;
      if(the_hgmr >= max_ok_hgmr) break; // all the rest are worse, so skip them.
      if(the_idx != fparent_idx  &&  the_idx != mparent_idx){
	push_to_vlong(best_parent_candidate_idxs, the_idx);
      }     
    }
  }
  free(the_idxhgmrs);

  long ub = long_min(best_parent_candidate_idxs->size, 10000); // set the number of possible parents to consider.
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

	if(n_over_d(alt_pedigree_stats->z) <= max_ok_z
	   &&
	   n_over_d(alt_pedigree_stats->d) <= max_ok_d){
	  alt_pedigree->pedigree_stats = alt_pedigree_stats;
	  push_to_vpedigree(alt_pedigrees, alt_pedigree);
	}  /* */
      }
    }
  }
  the_pedigree->A->search_done = true;
  sort_vpedigree_by_d(alt_pedigrees);
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





