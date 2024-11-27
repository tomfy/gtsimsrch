#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gtset.h"
#include <stdbool.h>
#include "pedigree.h"

#define PEDIGREE_FIELDS 7 // number of whitespace-separated fields in pedigree file, with ids in last 3 fields.
// extern int do_checks; // option -c sets this to 1 to do some checks.

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent){ //IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid){
  Pedigree* the_pedigree = (Pedigree*)malloc(sizeof(Pedigree));
  the_pedigree->F = Fparent;
  the_pedigree->M = Mparent;
  the_pedigree->A = Acc;
  the_pedigree->pedigree_stats = construct_pedigree_stats(); // NULL;
  return the_pedigree;
}


three_longs count_crossovers(GenotypesSet* the_gtsset, Accession* parent, Accession* offspring){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  if(parent == NULL  ||  offspring == NULL) {
    return (three_longs){-1, -1, -1};
  }
  long Xmin = 0, Xmax = 0, Nhet = 0; // number of heterozyg gts in parent
  long Xa = 0, Xb = 0;
  long prev_chrom_number = -1, prev_phase_a = -1, prev_phase_b = -1;
  long phase_a = -1, phase_b = -1, chrom_number;
 
  for(long i=0; i < parent->genotypes->length; i++){

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    char o_phase = offspring->phases->a[i];

       chrom_number = the_gtsset->chromosomes->a[i];
       // fprintf(stderr, "i:  %ld  prevchr, chr: %ld %ld   Xa, Xb:  %ld %ld   %ld %ld  %ld %ld\n", i, prev_chrom_number, chrom_number, Xa, Xb, prev_phase_a, prev_phase_b, phase_a, phase_b);
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
    //if(p_gt != '1') continue; // skip if parent gt not heterozyg, i.e. if homozyg or missing.

    if(p_gt == '1'){
      Nhet++; // counts the number of markers which are heterozyg in the parent, and non-missing in the offspring

      two_longs phases_ab = get_1marker_phases_wrt_1parent(p_phase, o_gt, o_phase);
      phase_a = phases_ab.l1;
      phase_b = phases_ab.l2;
    }
    // ######################################

 

    // compare current phases with previous values,
    // update crossover counts, and
    // and  update prev_phase_a, prev_phase_b
    if(prev_phase_a >= 0  &&  phase_a != prev_phase_a) Xa++; // phase has changed - crossover
    prev_phase_a = phase_a;
    if(prev_phase_b >= 0  &&  phase_b != prev_phase_b) Xb++; // phase has changed - crossover
    prev_phase_b = phase_b;
    
  } // end loop over markers
  // fprintf(stderr, "XXprevchr, chr: %ld %ld   Xa, Xb:  %ld %ld \n",  prev_chrom_number, chrom_number, Xa, Xb);
  if(Xa < Xb){ // add the crossovers from the last chromosome.
    Xmin += Xa; Xmax += Xb;
  }else{
    Xmin += Xb; Xmax += Xa;
  }
 
  return (three_longs){Xmin, Xmax, Nhet};
} // end of count_crossovers

Xcounts_2 count_crossovers_one_chromosome(GenotypesSet* the_gtsset, Accession* parent, Accession* offspring, long first, long next){
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

    // compare current phases with previous values,
    // update crossover counts, and
    // and  update prev_phase_a, prev_phase_b
    if(prev_phase_a >= 0  &&  phase_a != prev_phase_a) Xa++; // phase has changed - crossover
    prev_phase_a = phase_a;
    if(prev_phase_b >= 0  &&  phase_b != prev_phase_b) Xb++; // phase has changed - crossover
    prev_phase_b = phase_b;
    
  } // end loop over markers
  return (Xcounts_2){Xa, Xb, Nhet};
} // end of count_crossovers_one_chromosome

Xcounts_3 count_crossovers_two_parents(GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  Xcounts_3 result = (Xcounts_3){(Xcounts_2){-1,-1,-1}, (Xcounts_2){-1,-1,-1}, -1, -1, -1, -1};
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
    
    // long chrnumber = the_gtsset->chromosomes->a[start_index];
    // fprintf(stderr, "M  %ld %ld   %ld %ld  %ld\n", i, chrnumber, MX.Xa, MX.Xb, MX.Nhet);
    // fprintf(stderr, "F  %ld %ld   %ld %ld  %ld\n", i, chrnumber, FX.Xa, FX.Xb, FX.Nhet);

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
  Xcounts_2 FX2 = (Xcounts_2){XFmin_2, XFmax_2, NhetF};
  Xcounts_2 MX2 = (Xcounts_2){XMmin_2, XMmax_2, NhetM};
  result.FA = FX2;
  result.MA = MX2;
  result.XFmin_3 = XFmin_3;
  result.XFmax_3 = XFmax_3;
  result.XMmin_3 = XMmin_3;
  result.XMmax_3 = XMmax_3;
  return result;
}

Xover_info count_crossovers_two_parents_old(GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  Xover_info result = (Xover_info){-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  if( offspring == NULL
      || Fparent == NULL
      || Mparent == NULL) return result;
  
  long prev_chrom_number = -1, chrom_number;
  
  long prev_phase_Fa = -1, prev_phase_Fb = -1;
  long phase_Fa, phase_Fb;
  long XFa = 0, XFb = 0;

   long prev_phase_Ma = -1, prev_phase_Mb = -1;
  long phase_Ma, phase_Mb;
  long XMa = 0, XMb = 0;

  long XFmin = 0, XFmax = 0, NFhet = 0;
  long XMmin = 0, XMmax = 0, NMhet = 0;

  long XFmin_triple = 0, XFmax_triple = 0;
  long XMmin_triple = 0, XMmax_triple = 0;
  //long Xmin = 0, Xmax = 0; // number of heterozyg gts in parent
  
  for(long i=0; i < offspring->genotypes->length; i++){
    
    char o_gt = offspring->genotypes->a[i];
    char o_phase = offspring->phases->a[i];  
    if(o_gt == MISSING_DATA_CHAR) continue;

    chrom_number = the_gtsset->chromosomes->a[i];
    if(chrom_number != prev_chrom_number){ // now on next chromosome
      if(XFa < XFb){ // add chromosome Xmin, Xmax to totals	
	XFmin += XFa; XFmax += XFb;
      }else{
	XFmin += XFb; XFmax += XFa;
      }
      if(XMa < XMb){ // add the crossovers from the last chromosome.
	XMmin += XMa; XMmax += XMb;
      }else{
	XMmin += XMb; XMmax += XMa;
      }
      long X_Fa_Mb = XFa + XMb; // crossovers for F parent of a, M parent of b.
      long X_Fb_Ma = XFb + XMa; // crossovers for F parent of b, M parent of a.
      if(X_Fa_Mb < X_Fb_Ma){
	XFmin_triple += XFa;
	XFmax_triple += XFb;
	XMmin_triple += XMb;
	XMmax_triple += XMa;
      }else{
	XFmin_triple += XFb;
	XFmax_triple += XFa;
	XMmin_triple += XMa;
	XMmax_triple += XMb;
      }
      
      // reset for new chromosome:
         
      prev_chrom_number = chrom_number;
     
      phase_Fa = -1; phase_Fb = -1;
      prev_phase_Fa = -1; // needed
      prev_phase_Fb = -1;
      XFa = 0; XFb = 0;

      phase_Ma = -1; phase_Mb = -1;
      prev_phase_Ma = -1; // needed
      prev_phase_Mb = -1;
      XMa = 0; XMb = 0;
    }
    
    // ##############  Female parent  ######################
    if(1  ||  Fparent != NULL){
    char Fp_gt = Fparent->genotypes->a[i];
    if(Fp_gt == '1'){ // skip if parent gt not heterozyg, i.e. if homozyg or missing.
      char Fp_phase = Fparent->phases->a[i]; 
  
      NFhet++; // counts the number of markers which are heterozyg in the parent, and non-missing in the offspring

      two_longs phases_ab = get_1marker_phases_wrt_1parent(Fp_phase, o_gt, o_phase);
      phase_Fa = phases_ab.l1;
      phase_Fb = phases_ab.l2;
    }
    // compare current phases with previous values,
    // update crossover counts, and
    // and  update prev_phase_a, prev_phase_b
    if(prev_phase_Fa >= 0  &&  phase_Fa != prev_phase_Fa) XFa++; // phase has changed - crossover
    prev_phase_Fa = phase_Fa;
    if(prev_phase_Fb >= 0  &&  phase_Fb != prev_phase_Fb) XFb++; // phase has changed - crossover
    prev_phase_Fb = phase_Fb;
    }
    // ####################################################

    
    // ##############  Male parent  #######################
    if(1  ||  Mparent != NULL){
      char Mp_gt = Mparent->genotypes->a[i];
      if(Mp_gt == '1'){ // skip if parent gt not heterozyg, i.e. if homozyg or missing.
	char Mp_phase = Mparent->phases->a[i]; 
  
	NMhet++; // counts the number of markers which are heterozyg in the parent, and non-missing in the offspring

	two_longs phases_ab = get_1marker_phases_wrt_1parent(Mp_phase, o_gt, o_phase);
	phase_Ma = phases_ab.l1;
	phase_Mb = phases_ab.l2;

      }
      // compare current phases with previous values,
      // update crossover counts, and
      // and  update prev_phase_a, prev_phase_b
      if(prev_phase_Ma >= 0  &&  phase_Ma != prev_phase_Ma) XMa++; // phase has changed - crossover
      prev_phase_Ma = phase_Ma;
      if(prev_phase_Mb >= 0  &&  phase_Mb != prev_phase_Mb) XMb++; // phase has changed - crossover
      prev_phase_Mb = phase_Mb;
    }
    // ####################################################

    
  } // ##########  end loop over markers  ##################
  if(XFa < XFb){ // add the crossovers from the last chromosome.
    XFmin += XFa; XFmax += XFb;
  }else{
    XFmin += XFb; XFmax += XFa;
  }
   if(XMa < XMb){ // add the crossovers from the last chromosome.
    XMmin += XMa; XMmax += XMb;
  }else{
    XMmin += XMb; XMmax += XMa;
  }
      long X_Fa_Mb = XFa + XMb; // crossovers for F parent of a, M parent of b.
      long X_Fb_Ma = XFb + XMa; // crossovers for F parent of b, M parent of a.
      if(X_Fa_Mb < X_Fb_Ma){
	XFmin_triple += XFa;
	XFmax_triple += XFb;
	XMmin_triple += XMb;
	XMmax_triple += XMa;
      }else{
	XFmin_triple += XFb;
	XFmax_triple += XFa;
	XMmin_triple += XMa;
	XMmax_triple += XMb;
      }
      
      result = (Xover_info){XFmin, XFmax, NFhet, XMmin, XMmax, NMhet,
			    XFmin_triple, XFmax_triple, XMmin_triple, XMmax_triple};
  return result;
} // end of count_chromosome_crossovers

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
  //  pedigree_stats->agmr12 = {0, -1};
  //  pedigree_stats->par1_hgmr = {0, -1};
  /* ND par1_R; */
  /* ND par2_hgmr; */
  /* ND par2_R; */
  /* ND z; // (n00_1 + n22_1)/(n00_x + n22_x) */
  /* ND d; */
  /* // ND pseudo_hgmr; */
  /* ND xhgmr1; */
  /* ND xhgmr2; */
  /* ND d_22; // both parents homozyg, delta = 2  */
  /* ND d_21; // both parents homozyg, delta = 1 */
  /* ND d_11; // one  */
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
  // fprintf(stderr, "n00_2, n22_0, n_2: %ld %ld %ld \n", n_00_2, n_22_0, n_2);
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

  // 'apparent': (n_a2 + 0.5*n_a1)/(n_00 + n_22 + n_02 + n_20)
  // find_parents: (n_a2 + n_a1 + n_o1)/(n_00 + n_02 + n_20 + n_22 + n_01 + n_10 + n_12 + n_21);
  // to try:  (n_a2 + wa*n_a1 + wo*n_o1))/(n_00 + n_02 + n_20 + n_22 + wx*(n_01 + n_10 + n_12 + n_21));
  // so wa=1/2, wo=wx=0
  // also just look at n_o1/(n_01 + n_10 + n_02 + n_20);
  // ************************************
  ND d_22 = {n_a2, n_00 + n_22}; // parents both homozyg, delta = 2
  pedigree_stats->d_22 = d_22;
  ND d_21 = {n_a1, n_00 + n_22 + n_02 + n_20}; // parents both homozyg, delta = 1
  pedigree_stats->d_21 = d_21;
  
  ND d_11 = {n_o1, n_01 + n_10 + n_12 + n_21}; // one parents homozyg, one heterozyg, delta = 1;
  //  fprintf(stderr, "### %ld %ld \n", d_11.n, d_11.d);
  pedigree_stats->d_11 = d_11; 
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0; // par2 gt can be missing
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1; // par2 gt can be missing
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2; // par2 gt can be missing
  
  /* long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0; // par2 gt can be missing
  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1; // par2 gt can be missing
  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2; // par2 gt can be missing
  /* */
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0; // par2 gt can be missing
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1; // par2 gt can be missing
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2; // par2 gt can be missing

  long n_3x_0 = n_30_0 + n_31_0 + n_32_0; // + n_03_0;
  long n_3x_1 = n_30_1 + n_31_1 + n_32_1; // + n_03_1;
  long n_3x_2 = n_30_2 + n_31_2 + n_32_2; // + n_03_2;

  long z_numer = n_00_1 + n_22_1;
  long z_denom = z_numer + n_00_0 + n_22_2 + n_00_2 + n_22_0; // = n_00_x + n_22_x
  // fprintf(stderr, "%ld %ld  %ld %ld %ld   %ld %ld\n", z_numer, z_denom, n_0, n_1, n_2, n_1 + n_2, n_0 _1 + n_2 );
  /* long total =
    n_0x_0 + n_0x_1 + n_0x_2
    + n_1x_0 + n_1x_1 + n_1x_2
    + n_2x_0 + n_2x_1 + n_2x_2
    + n_3x_0 + n_3x_1 + n_3x_2
    + n_33_x + n_xy_3; /* */

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0; // par1 gt can be missing
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1; // par1 gt can be missing
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2; // par1 gt can be missing
  
  //  long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0;
  //  long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1;
  //  long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2;
  
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
  // fprintf(stderr, "rgtc: %ld %ld %ld %ld \n", r1_numer, r1_denom, r2_numer, r2_denom);
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

  pedigree_stats->d_2 = (ND) {n_1 + n_2 - n_00_1 - n_22_1, n_0 + n_1 + n_2 - n_11_x}; // d_nd;

  pedigree_stats->d = (ND){n_1 + n_2, n_0 + n_1 + n_2};
  // fprintf(stderr, "rg d_old. N, D, N/D:  %ld  %ld  %8.5f  n1 n2: %ld %ld\n", d_old_numer, d_old_denom, (d_old_denom > 0)? (double)d_old_numer/(double)d_old_denom : -1, n_1, n_2);
   // fprintf(stderr, "n1, n2, n0,  num, denom: %ld  %ld  %ld   %ld %ld  n1 n2: %ld %ld\n", n_1, n_2, n_0, n_1+n_2, n_1+n_2+n_0, n_1, n_2);
  pedigree_stats->z = (ND) {z_numer, z_denom};

  /* double scaled_d = pedigree_stats->scaled_d; */
  /* double z = n_over_d(pedigree_stats->z); */
  /* pedigree_stats->max_scaleddz = (scaled_d>z)? scaled_d : z; */

  pedigree_stats->all_good_count = n_0 + n_1 + n_2;

  return pedigree_stats;
} // end of triple_counts

Pedigree_stats* bitwise_triple_counts(Accession* par1, Accession* par2, Accession* prog){ // , GenotypesSet* gtset){
  
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
    
    n_0x_1_2x_1 += __builtin_popcountll(is_i0or2_k1);
    n_0x_0_2x_2 += __builtin_popcountll(is0x_0_or_2x_2);   
    n_x0_1_x2_1 += __builtin_popcountll(is_j0or2_k1);
    n_x0_0_x2_2 += __builtin_popcountll(isx0_0_or_x2_2);   
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
	
  Pedigree_stats* pedigree_stats = construct_pedigree_stats(); //(Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  // if(ndiff12 > n_total_no_md) fprintf(stderr, "%ld  %ld\n", ndiff12, n_total_no_md);
  pedigree_stats->agmr12 = (ND) {ndiff12, n_total_no_md}; // agmr between parents
  pedigree_stats->agmr01 = (ND) {ndiff01, n_total_no_md}; 
  pedigree_stats->agmr02 = (ND) {ndiff02, n_total_no_md};

  pedigree_stats->d_2 = (ND) {n_0xorx0_2_nomd + n_2xorx2_0_nomd, n_total_no_md - n_total_x_11};
  /* long n_1 = n_0xorx0_2_nomd + n_2xorx2_0_nomd + n_00_1_22_1 - n_00_2_22_0;; */
  /* long n_2 = n_00_2_22_0; */
  /* long d_old_numer = n_1 + n_2; */
  /* long d_old_denom = n_total_no_md; */
  pedigree_stats->d = (ND) {n_0xorx0_2_nomd + n_2xorx2_0_nomd + n_00_1_22_1, n_total_no_md};
  // fprintf(stderr, "d d d %8.5f\n", n_over_d(pedigree_stats->d));
  // fprintf(stdout, "d num, denom:  %ld %ld\n",  n_0xorx0_2_nomd + n_2xorx2_0_nomd, n_total_no_md - n_total_x_11);
  pedigree_stats->z = (ND) {n_00_1_22_1, n_00_22};
  pedigree_stats->par1_hgmr = (ND) {hgmr1_numerator, hgmr1_denominator};
  pedigree_stats->par2_hgmr = (ND) {hgmr2_numerator, hgmr2_denominator};
  pedigree_stats->par1_R = (ND) {n_0x_1_2x_1, n_0x_1_2x_1 + n_0x_0_2x_2};
  pedigree_stats->par2_R = (ND) {n_x0_1_x2_1, n_x0_1_x2_1 + n_x0_0_x2_2};
  /* double scaled_d = pedigree_stats->scaled_d; */
  /* double z = n_over_d(pedigree_stats->z); */
  /* pedigree_stats->max_scaleddz = (scaled_d>z)? scaled_d : z; */
  //fprintf(stderr, "bwtc: %ld %ld %ld %ld \n", n_0x_1_2x_1, n_0x_1_2x_1 + n_0x_0_2x_2, n_x0_1_x2_1, n_x0_1_x2_1 + n_x0_0_x2_2);
  //fprintf(stderr, "bw d_old. N, D, N/D:  %ld  %ld  %8.5f n1 n2: %ld %ld\n", d_old_numer, d_old_denom, (d_old_denom > 0)? (double)d_old_numer/(double)d_old_denom : -1, n_1, n_2 );
  // fprintf(stderr, "xxx: %ld %ld %ld %ld  %ld %ld\n", n_0xorx0_2_nomd, n_2xorx2_0_nomd, n_00_1_22_1, n_00_2_22_0, n_00_2, n_22_0);
  pedigree_stats->n_01or10_1 = n_01or10_1;
  pedigree_stats->all_good_count = n_total_no_md;
  
  return pedigree_stats;
} // end of bitwise_triple_counts

Vpedigree*  calculate_triples_for_one_accession(Accession* prog, GenotypesSet* the_genotypes_set, Viaxh* cppps, long max_candidate_parents){

  // sort the parent candidates and keep the best ones
  long limited_n_candpairs = cppps->size;
  if(cppps->size == 0){ 
    // count_accs_w_no_cand_parents++;
    // fprintf(o_stream, "# %s has no candidate parents (xhmgr <= %8.5f)\n", prog->id->a, max_xhgmr);-
  }else if(cppps->size > max_candidate_parents){ // if too many parent candidates, just take the max_candidate_parents best ones
    sort_viaxh_by_xhgmr(cppps);
    limited_n_candpairs = max_candidate_parents; // limit candidate to max_candidate_parents
    // count_accs_w_too_many_cand_parents++;
  }
  //   fprintf(stderr, "cppps->size: %ld limited_n_candpairs: %ld\n", cppps->size, limited_n_candpairs);


  
  Vpedigree* alt_pedigrees = construct_vpedigree(1000);
  for(long ii=0; ii<limited_n_candpairs; ii++){
    long par1idx = cppps->a[ii]->idx;
    //  Iaxh* iaxh = cppps->a[ii];
    //   fprintf(stderr, "ii: %ld   %ld  %8.5f %8.5f %8.5f\n", ii, iaxh->idx, iaxh->agmr, iaxh->hgmr, iaxh->xhgmr);
    Accession* par1 = the_genotypes_set->accessions->a[par1idx];
    for(long jj=ii; jj<limited_n_candpairs; jj++){	
      long par2idx = cppps->a[jj]->idx;
      // Iaxh* iaxhjj = cppps->a[jj];
      //	fprintf(stderr, "jj: %ld   %ld  %8.5f %8.5f %8.5f\n", jj, iaxhjj->idx, iaxhjj->agmr, iaxhjj->hgmr, iaxhjj->xhgmr);
      Accession* par2 = the_genotypes_set->accessions->a[par2idx];
      // fprintf(stderr, "# %s %s %s \n", prog->id->a, par1->id->a, par2->id->a);
      Pedigree_stats* the_ps = bitwise_triple_counts(par1, par2, prog);

      the_ps->xhgmr1 = cppps->a[ii]->xhgmr;
      the_ps->xhgmr2 = cppps->a[jj]->xhgmr;

      the_ps->hgmr1 = cppps->a[ii]->hgmr;
      the_ps->hgmr2 = cppps->a[jj]->hgmr;
	
      Pedigree* the_pedigree = construct_pedigree(prog, par1, par2);
      the_pedigree->pedigree_stats = the_ps;
      //	the_ps->xhgmr1 = xhgmr(the_genotypes_set, par1, prog, 0); // do full (not 'quick') xhgmr
      //	the_ps->xhgmr2 = xhgmr(the_genotypes_set, par2, prog, 0); // do full (not 'quick') xhgmr
	
      push_to_vpedigree(alt_pedigrees, the_pedigree);
      //	fprintf(stderr, "AFM: %s    %s %s   %ld \n\n", the_pedigree->A->id->a, the_pedigree->F->id->a, the_pedigree->M->id->a, alt_pedigrees->size); 
    } // end loop over parent 2
  } // end loop over parent 1
  return alt_pedigrees;
} //

ND xz(GenotypesSet* gtset, Accession* O, Accession* P1, Accession* P2){
  //fprintf(stderr, "top of Zn\n");
  long numerator = 0;
  double denominator = 0;
  for(long i = 0; i<P1->genotypes->length; i++){
    char gt1 = P1->genotypes->a[i];
    if(gt1 == '0'){
      char gt2 = P2->genotypes->a[i];
      char gtO = O->genotypes->a[i];
      if(gt2 == '0'){
	if(gtO == '1') numerator++;
	double Ohet_fraction = (double)gtset->marker_dosage_counts[1]->a[i]/(double)gtset->accessions->size;
	denominator += Ohet_fraction;
      }
    }else if(gt1 == '2'){
      char gt2 = P2->genotypes->a[i];
      char gtO = O->genotypes->a[i];
      if(gt2 == '2'){
	if(gtO == '1') numerator++;
	double Ohet_fraction = (double)gtset->marker_dosage_counts[1]->a[i]/(double)gtset->accessions->size;
	denominator += Ohet_fraction;
      }
    }
  }
  //fprintf(stderr, "bottom of xz\n");
  // fprintf(stderr, "xz, numerator: %ld  denominator: %7.5f\n", numerator, denominator);
  return (ND) {numerator, (long)(denominator+0.5)};
}

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

  the_ps->scaled_d = NAN;
  the_ps->max_scaleddz = NAN;
  the_ps->xhgmr1 = NAN;
  the_ps->xhgmr2 = NAN;
  return the_ps;
}
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset, double d_scale_factor){ //, long* d0counts, long* d1counts, long* d2counts){ //, GenotypesSet* the_gtsset){
  long ploidy = the_gtsset->ploidy;
  Pedigree_stats* the_ps; //  = construct_pedigree_stats(); // (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
  assert(the_pedigree->F != NULL  ||  the_pedigree->M != NULL); // shouldn't have both parents NULL //
  if(the_pedigree->F != NULL  &&  the_pedigree->M != NULL){
   
    // the_ps = triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a, ploidy );
     the_ps = bitwise_triple_counts(the_pedigree->F,  the_pedigree->M,  the_pedigree->A);

    /* if(0){ // compare bitwise, nonbitwise calculations as check */
    /*   Pedigree_stats* nobw_ps = triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a, ploidy );  */
    /*   assert(NDs_equal(the_ps->agmr12, nobw_ps->agmr12)); */
    /*   assert(NDs_equal(the_ps->par1_hgmr, nobw_ps->par1_hgmr)); */
    /*   assert(NDs_equal(the_ps->par1_R, nobw_ps->par1_R)); */
    /*   assert(NDs_equal(the_ps->par2_hgmr, nobw_ps->par2_hgmr)); */
    /*   assert(NDs_equal(the_ps->par2_R, nobw_ps->par2_R)); */
    /*   assert(NDs_equal(the_ps->d, nobw_ps->d)); */
    /*   assert(NDs_equal(the_ps->z, nobw_ps->z)); */
    /* } */
     /* double d = n_over_d(the_ps->d); */
     /* fprintf(stderr, "d: %8.5lf  ", d); */
     /* d /= the_gtsset->mean_d; // NAN if mean_d is 0 */
     /* //pedigree_stats->s */
     /* fprintf(stderr, " %8.5f   ", d); */
     double scaled_d = n_over_d(the_ps->d)*d_scale_factor/the_gtsset->mean_d;
     the_ps->scaled_d = scaled_d;
      
     // double scaled_d = pedigree_stats->scaled_d;
     double zn = n_over_d(the_ps->z)/the_gtsset->mean_z;
     if(! isnan(zn) && !isnan(scaled_d)){
       the_ps->max_scaleddz = (scaled_d>zn)? scaled_d : zn;
     }
     // fprintf(stderr, " %8.5lf  %8.5lf  %8.5lf\n", d_scale_factor, scaled_d, zn);
       
     the_ps->par1_xhgmr = xhgmr(the_gtsset, the_pedigree->F, the_pedigree->A, false);
     the_ps->par2_xhgmr = xhgmr(the_gtsset, the_pedigree->M, the_pedigree->A, false);
     the_ps->xz = xz(the_gtsset, the_pedigree->A, the_pedigree->F, the_pedigree->M);
  }else{ // one of the parents is NULL
    the_ps = construct_pedigree_stats();
    the_ps->agmr12 = (ND) {0, 0};
    the_ps->z = (ND) {0, 0};
    the_ps->xz = (ND) {0, 0};
    // return the_ps;
    if(the_pedigree->F != NULL){ // we have female parent id, no male parent id
      //       fprintf(stderr, "pedigree with female parent only.\n");
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
    //construct_pedigree_stats(the_pedigree); //the_ps = (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
    // the_ps->agmr12.n = 0;
    //   fprintf(stderr, "######  %ld %ld   %ld %ld\n", the_ps->xhgmr1.n, the_ps->xhgmr1.d, the_ps->xhgmr2.n, the_ps->xhgmr2.d);
  }
  the_ps->xhgmr1 = n_over_d(the_ps->par1_xhgmr);
  the_ps->xhgmr2 = n_over_d(the_ps->par2_xhgmr);
  return the_ps;
}

// get the indices of all the accessions which, according to the_vped, have offspring.
const Vlong* accessions_with_offspring(const Vpedigree* the_vped, long n_accessions){ // , long n_accessions){
  //fprintf(stderr, "n peds: %ld\n", the_vped->size);
  Vlong* offspring_counts = construct_vlong_zeroes(n_accessions);
  for(long i=0; i<the_vped->size; i++){
    const Pedigree* the_ped = the_vped->a[i];
    //fprintf(stderr, "i:  %ld   A index: %ld\n", i, the_vped->a[i]->A->index);
    // fprintf(stderr, "%ld %ld \n", the_ped->F->index, the_ped->Fparent->index);
    if(the_ped->F != NULL) {
      // fprintf(stderr, "n peds: %ld\n", the_vped->size);
      //fprintf(stderr, "F index: %ld\n", the_ped->F->index);
      offspring_counts->a[the_ped->F->index]++;
      //fprintf(stderr, "i: %ld  F index: %ld\n", i, the_ped->F->index);
    }
    //else{fprintf(stderr, "F is NULL\n");}
    if(the_ped->M != NULL) {
      // fprintf(stderr, "n peds: %ld\n", the_vped->size);
      // fprintf(stderr, "M index: %ld\n", the_ped->M->index);
      offspring_counts->a[the_ped->M->index]++;
      // fprintf(stderr, "i: %ld  M index: %ld\n", i, the_ped->M->index);
    }
    //else{fprintf(stderr, "M is NULL\n");
  }
  //fprintf(stderr, "XYZ\n");

  //fprintf(stderr, "XXXXX %ld \n", offspring_counts->size);
  Vlong* accidxs_with_offspring = construct_vlong(100);
  //  Vaccession* accessions = construct_vaccession(100);
  for(long i=0; i<offspring_counts->size; i++){
    if(offspring_counts->a[i] > 0){
      push_to_vlong(accidxs_with_offspring, i);
      //    add_accession_to_vaccession(accessions, the_gtsset->accessions->a[i];
    }
  }
  free_vlong(offspring_counts);
  return accidxs_with_offspring;
}

void print_pedigree(FILE* fh, Pedigree* the_pedigree, bool verbose){
  Accession* F = the_pedigree->F;
  Accession* M = the_pedigree->M;
  fprintf(fh, "%s  %s  %ld  ", (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA",  the_pedigree->pedigree_stats->all_good_count);
  print_pedigree_stats(fh, the_pedigree->pedigree_stats, verbose);
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
    print_n_over_d(fh, the_pedigree_stats->agmr12, 1.0);
    print_n_over_d(fh, the_pedigree_stats->par1_hgmr, 1.0);
    // print_n_over_d(fh, the_pedigree_stats->par1_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par1_R, 1.0);
    print_n_over_d(fh, the_pedigree_stats->par2_hgmr, 1.0);
    //  print_n_over_d(fh, the_pedigree_stats->par2_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par2_R, 1.0);
    print_n_over_d(fh, the_pedigree_stats->d, 1.0);
    print_n_over_d(fh, the_pedigree_stats->z, 1.0);
    // print_n_over_d(fh, the_pedigree_stats->d_old);
  }
  fprintf(fh, "%7.5f  ", the_pedigree_stats->scaled_d);
  fprintf(fh, "%7.5lf  ", the_pedigree_stats->max_scaleddz);
  //   fprintf(fh, "%7.5lf  ", the_pedigree_stats->hgmr1);
  fprintf(fh, "%7.5lf  ", the_pedigree_stats->xhgmr1);
  //  fprintf(fh, "%7.5lf  ", the_pedigree_stats->hgmr2);
  fprintf(fh, "%7.5lf  ", the_pedigree_stats->xhgmr2);	     	   
}

void print_pedigree_normalized(FILE* fh, Pedigree* the_pedigree, GenotypesSet* gtset){
  //double mean_hgmr, double mean_R, double mean_d, double mean_z){
  Accession* F = the_pedigree->F;
  Accession* M = the_pedigree->M;
  fprintf(fh, "%s  %s  %ld  ", (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA",  the_pedigree->pedigree_stats->all_good_count);
  print_normalized_pedigree_stats(fh, the_pedigree->pedigree_stats, gtset);
}

void print_normalized_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats, GenotypesSet* gtset){ 
  double mean_hgmr = gtset->mean_hgmr;
  double mean_R = gtset->mean_R;
  double mean_d = gtset->mean_d;
  double mean_z = gtset->mean_z;
  // print ratios but not denominators
  print_n_over_d(fh, the_pedigree_stats->agmr12, 1.0);
  print_n_over_d(fh, the_pedigree_stats->par1_hgmr, mean_hgmr);
  // print_n_over_d(fh, the_pedigree_stats->par1_xhgmr);
  print_n_over_d(fh, the_pedigree_stats->par1_R, mean_R);
  print_n_over_d(fh, the_pedigree_stats->par2_hgmr, mean_hgmr);
  //  print_n_over_d(fh, the_pedigree_stats->par2_xhgmr);
  print_n_over_d(fh, the_pedigree_stats->par2_R, mean_R);
  print_n_over_d(fh, the_pedigree_stats->d, mean_d);
  print_n_over_d(fh, the_pedigree_stats->z, mean_z);
  // print_n_over_d(fh, the_pedigree_stats->d_old);
  
  // fprintf(fh, "%7.5f  ", the_pedigree_stats->scaled_d);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->scaled_d);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->max_scaleddz);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->xhgmr1);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->xhgmr2);
  // fprintf(fh, "%7.5lf  ", the_pedigree_stats->max_scaleddz);
  // fprintf(fh, "%7.5lf  ", the_pedigree_stats->hgmr1);
  // fprintf(fh, "%7.5lf  ", the_pedigree_stats->xhgmr1);
  // fprintf(fh, "%7.5lf  ", the_pedigree_stats->hgmr2);
  // fprintf(fh, "%7.5lf  ", the_pedigree_stats->xhgmr2);	
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

void print_Xcounts_2(FILE* fh, Xcounts_2 X){
  if(X.Xa < X.Xb){
    fprintf(fh, " %ld %ld %ld ", X.Xa, X.Xb, X.Nhet);
  }else{
        fprintf(fh, " %ld %ld %ld ", X.Xb, X.Xa, X.Nhet);
  }
}

long pedigree_ok_x(Pedigree_stats* p, double max_self_agmr12, double max_ok_hgmr, double max_self_r, double max_ok_z){
  // returns 2 for ok biparental, 1 for ok self, 0 for bad
  double agmr12 = n_over_d(p->agmr12);
  double hgmr1 = n_over_d(p->par1_hgmr);
  double r1 = n_over_d(p->par1_R);
  double hgmr2 = n_over_d(p->par2_hgmr);
  double r2 = n_over_d(p->par2_R);
  double z = n_over_d(p->z);
  long result = 0;
  if(agmr12 <= max_self_agmr12){ // pedigree says self (or parents very similar)
    if( /*(agmr12 <= max_self_agmr12) && */ (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 <= max_self_r) && (r2 <= max_self_r) && (z <= max_ok_z) ){
      result = 1;
    }
  }else{ // pedigree says biparental
    if( (hgmr1 <= max_ok_hgmr) && (hgmr2 <= max_ok_hgmr) && (r1 > max_self_r) && (r2 > max_self_r) && (z <= max_ok_z) ){
      result = 2;
    }
  }
  return result;
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
      //fprintf(stderr, "# acc_idx: %ld\n", acc_idx);
      if(acc_idx == ID_NA_INDEX){ // no genotypes for this accession
	not_in_genotypes_set_count++;
      }else{ // have genotypes for this accession
	fempar_idx = index_of_id_in_vidxid(the_gt_vidxid, fempar_id);
	//fprintf(stderr, "# fem_idx: %ld\n", fempar_idx);
	malpar_idx = index_of_id_in_vidxid(the_gt_vidxid, malpar_id);	 
	//fprintf(stderr, "# mal_idx: %ld\n", malpar_idx);
	//	if(1){
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
 
    push_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      push_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
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
	push_to_vpedigree(pedigrees, a_pedigree);
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
  return the_vped;
}
  
void push_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped){
  long cap = the_vped->capacity;
  long n = the_vped->size;
  //fprintf(stderr, "### cap: %ld  size: %ld\n", cap, n);
  if(n == cap){
    cap *= 2;
    the_vped->a = (Pedigree**)realloc(the_vped->a, cap*sizeof(Pedigree*));
    the_vped->capacity = cap;
  }
  //fprintf(stderr, "# the_ped %p\n", the_ped);
  the_vped->a[n] = the_ped;
  the_vped->size++;
  //fprintf(stderr, "return from push_to_vped\n");
}



void sort_vpedigree_by_maxdz(Vpedigree* the_vped){ // sort by max(scaled_d, z) (increasing) 
  qsort(the_vped->a, the_vped->size, sizeof(Pedigree*), compare_pedigree_maxdz);
}

int compare_pedigree_maxdz(const void* a, const void* b){
  //int retval;
  double d1 = (*((Pedigree**)a))->pedigree_stats->scaled_d;
  double d2 = (*((Pedigree**)b))->pedigree_stats->scaled_d;
  
  //  double d1 = (dnd1.d > 0)? (double)dnd1.n/dnd1.d : 1.0;
  // double d2 = (dnd2.d > 0)? (double)dnd2.n/dnd2.d : 1.0;
  ND znd1 = (*((Pedigree**)a))->pedigree_stats->z;
  ND znd2 = (*((Pedigree**)b))->pedigree_stats->z;
  double z1 = (znd1.d > 0)? (double)znd1.n/znd1.d : 1.0;
  double z2 = (znd2.d > 0)? (double)znd2.n/znd2.d : 1.0;
  
  double x1 = (d1 > z1)? d1 : z1;
  double x2 = (d2 > z2)? d2 : z2;
  if(x1 > x2){
    // retval = 1;
     return 1;
  }else if(x1 < x2){
    //  retval = -1;
     return -1;
  }else{
    //   retval = 0;
     return 0;
  }
  // fprintf(stderr, "in compare... %7.4f %7.4f %7.4f\n", d1, z1, x1);
  //  return retval;
}

void sort_vpedigree_by_maxh1h2z(Vpedigree* the_vped){
  qsort(the_vped->a, the_vped->size, sizeof(Pedigree*), compare_pedigrees_h1h2z);
}

int compare_pedigrees_h1h2z(const void* a, const void* b){
  ND nd_h1_a = (*((Pedigree**)a))->pedigree_stats->par1_hgmr;
  ND nd_h2_a = (*((Pedigree**)a))->pedigree_stats->par2_hgmr;
  ND nd_z_a = (*((Pedigree**)a))->pedigree_stats->z;

  double h1_a = (nd_h1_a.d > 0)? (double)nd_h1_a.n/nd_h1_a.d : 1.0;
  double h2_a = (nd_h2_a.d > 0)? (double)nd_h2_a.n/nd_h2_a.d : 1.0;
  double z_a = (nd_z_a.d > 0)? (double)nd_z_a.n/nd_z_a.d : 1.0;
  
  ND nd_h1_b = (*((Pedigree**)b))->pedigree_stats->par1_hgmr;
  ND nd_h2_b = (*((Pedigree**)b))->pedigree_stats->par2_hgmr;
  ND nd_z_b = (*((Pedigree**)b))->pedigree_stats->z;

  double h1_b = (nd_h1_b.d > 0)? (double)nd_h1_b.n/nd_h1_b.d : 1.0;
  double h2_b = (nd_h2_b.d > 0)? (double)nd_h2_b.n/nd_h2_b.d : 1.0;
  double z_b = (nd_z_b.d > 0)? (double)nd_z_b.n/nd_z_b.d : 1.0;
 
  double xa = (h1_a > h2_a)? h1_a : h2_a;
  if(z_a > xa) xa = z_a;

  double xb = (h1_b > h2_b)? h1_b : h2_b;
  if(z_b > xb) xb = z_b;
 
  if(xa > xb){
    return 1;
  }else if(xa < xb){
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
  // fprintf(stderr, "n_parents: %ld \n", n_parents);
  for(long i=0; i<n_parents; i++){ // store 'good' hgmr indices 
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != acc_idx){
      double the_hgmr = the_idxhgmrs[i].hgmr;
      // fprintf(stderr, "%s %ld %8.5f\n", acc_id, the_idx, the_hgmr);
      if(the_hgmr >= max_ok_hgmr) break; // all the rest are worse, so skip them.
      //   if(the_hgmr >= 0){
      if(the_idx != fparent_idx  &&  the_idx != mparent_idx){
	push_to_vlong(best_parent_candidate_idxs, the_idx);
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
	//if(get_hgmr1(alt_pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree_stats) <= 0.05  &&

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
  sort_vpedigree_by_maxdz(alt_pedigrees);

  free_vlong(best_parent_candidate_idxs);
  return alt_pedigrees;
}


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
	push_to_vpedigree(alt_pedigrees, alt_pedigree);
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



// ******************************************************************************
// *****  unused ****************************************************************
// ******************************************************************************

/* long check_idxid_map(Vidxid* vidxid, const Vaccession* accessions){ */
/*   for(long i=0; i<accessions->size; i++){ */
/*     char* id = accessions->a[i]->id->a; */
/*     long idx = index_of_id_in_vidxid(vidxid, id); */
/*     // fprintf(stderr, "%ld %ld %s\n", i, idx, id); */
/*     if(idx != i) return 0; */
/*   } */
/*   return 1; */
/* } */


/* void sort_vpedigree_by_d(Vpedigree* the_vped){
  qsort(the_vped->a, the_vped->size, sizeof(Pedigree*), compare_pedigree_d);
}
int compare_pedigree_d(const void* a, const void* b){
  ND nd1 = (*((Pedigree**)a))->pedigree_stats->d;
  ND nd2 = (*((Pedigree**)b))->pedigree_stats->d;
  double d1 = (nd1.d > 0)? (double)nd1.n/nd1.d : 1.0;
  double d2 = (nd2.d > 0)? (double)nd2.n/nd2.d : 1.0;
  if(d1 > d2){
    return 1;
  }else if(d1 < d2){
    return -1;
  }else{
    return 0;
  }
} /* */

/* void sort_vpedigree_by_z(Vpedigree* the_vped){
  qsort(the_vped->a, the_vped->size, sizeof(Pedigree*), compare_pedigree_z);
}
int compare_pedigree_z(const void* a, const void* b){
  ND nd1 = (*((Pedigree**)a))->pedigree_stats->z;
  ND nd2 = (*((Pedigree**)b))->pedigree_stats->z;
  double d1 = (nd1.d > 0)? (double)nd1.n/nd1.d : 1.0;
  double d2 = (nd2.d > 0)? (double)nd2.n/nd2.d : 1.0;
  if(d1 > d2){
    return 1;
  }else if(d1 < d2){
    return -1;
  }else{
    return 0;
  }
} /* */


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

/* 
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
/* */

/*
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
      push_to_vlong(best_parent_candidate_idxs, the_idx);
          
      // }
    }
  }
  free(the_idxhgmrs);
  return best_parent_candidate_idxs;
} // end of alternative_parents
/* */

/*
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
} /* */

/* 
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
} /* */

/*
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
}/* */

/*
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
}/* */

/*
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
} /* */

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



/* two_longs diploid_quick_and_dirty_triple_counts(Accession* acc1, Accession* acc2, Accession* progacc){
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
  for(long i=0; i<acc2->alt_homozygs->size; i++){ // gt2 is '2'
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
} /* */

/*
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
 
  for(long i=0; i<acc1->ref_homozygs->size; i++){ // just look at markers with gt1 = 0
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
} /* */

/* void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset){ */

/*   char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index]; */
/*   char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index]; */
/*   char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index]; */
    
/*   Ahr fp_ahr, mp_ahr, fm_ahr; */
/*   agmr_hgmr_r(fempar_gts, acc_gts, &(the_pedigree->fp_ahr)); */
/*   agmr_hgmr_r(malpar_gts, acc_gts, &(the_pedigree->mp_ahr)); */
/*   agmr_hgmr_r(fempar_gts, malpar_gts, &(the_pedigree->fm_ahr)); */
/* } */


/*  // sorting an array of Pedigree_stats*
void sort_pedigree_stats_by_d(Pedigree_stats** the_pss, long size){ // sort in place
  qsort(the_pss, size, sizeof(Pedigree_stats*), pscmp);
}

int pscmp(const void* v1, const void* v2){
  const Pedigree_stats** s1 = (const Pedigree_stats**)v1;
  const Pedigree_stats** s2 = (const Pedigree_stats**)v2;
  long denom1 = ((*s1)->d).d;
  double d1 = (denom1 > 0)? (double)((*s1)->d).n/denom1 : 2.0;
  long denom2 = ((*s2)->d).d;
  double d2 = (denom2 > 0)? (double)((*s2)->d).n/denom2 : 2.0;
  if(d1 < d2){
    return -1;
  }else if(d1 > d2){
    return 1;
  }else{
    return 0;
  }
} /* */
