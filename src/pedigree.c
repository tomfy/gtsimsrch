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
  the_pedigree->pedigree_stats = NULL; // construct_pedigree_stats(); 
 
  return the_pedigree;
}

four_longs count_crossovers(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  if(offspring == NULL){
    fprintf(stderr, "# in count_crossovers offspring Accession* is NULL. Bye.\n");
    exit(EXIT_FAILURE);
  }
  /* Xcounts_2mmn FXcounts = {0, 0, 0}; */
  /* Xcounts_2mmn MXcounts = {0, 0, 0}; */
  /* if(Fparent != NULL) FXcounts = count_crossovers_one_parent(the_gtsset, Fparent, offspring); */
  /* if(Mparent != NULL) MXcounts = count_crossovers_one_parent(the_gtsset, Mparent, offspring); */
  /* Xcounts_3 X3 = (Xcounts_3){FXcounts, MXcounts, 0, 0, 0, 0}; */

  ND FXcounts = {0, 0};
  ND MXcounts = {0, 0};
  if(Fparent != NULL) FXcounts = count_crossovers_one_parent(the_gtsset, Fparent, offspring);
  if(Mparent != NULL) MXcounts = count_crossovers_one_parent(the_gtsset, Mparent, offspring);
  four_longs XX = (four_longs){FXcounts.n, FXcounts.d, MXcounts.n, MXcounts.d};
  
  return XX;
}


ND count_crossovers_one_parent(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  if(parent == NULL  ||  offspring == NULL) {
    return (ND){-1, -1};
  }
  long NX = 0, Nhet = 0;
  long prev_chrom_number = -1, chrom_number;
  long prev_phase = -1, phase = -1;
 
  for(long i=0; i < parent->genotypes->length; i++){

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    if(o_gt == '1') continue; // don't use markers heterozygous in offspring
    char o_phase = offspring->phases->a[i];
       chrom_number = the_gtsset->chromosomes->a[i];
    if(chrom_number != prev_chrom_number){ // now on next chromosome
      // reset for new chromosome:
      prev_chrom_number = chrom_number;
      phase = -1;
      prev_phase = -1; // needed
    }

    // ########################################
    char p_gt = parent->genotypes->a[i];
    char p_phase = parent->phases->a[i];

    if(p_gt == '1'){
      Nhet++; // counts the number of markers which are heterozyg in the parent, and homozygous in the offspring
      two_longs phases_ab = get_1marker_phases_wrt_1parent(p_phase, o_gt, o_phase);
      phase = phases_ab.l1;
    }
    // ######################################

    // compare current phases with previous values,
    // update crossover counts, and  update prev_phase
    if(prev_phase >= 0  &&  phase != prev_phase) NX++; // phase has changed - crossover
    prev_phase = phase;
    
  } // end loop over markers

  return (ND){NX, Nhet};
} // end of count_crossovers




two_longs get_1marker_phases_wrt_1parent(char p_phase, char o_gt, char o_phase){
  long phase_a, phase_b;
  // get 
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


Pedigree_stats* bitwise_triple_counts(Accession* par1, Accession* par2, Accession* prog){ //, GenotypesSet* the_gtset){
  
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
    unsigned long long isx_11 = i1 & j1 & ~missing; // both parents het, none missing
    
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
    //fprintf(stderr, "AAAAAAAAA: %ld  %ld \n", n_00_2, n_22_0);

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
    n_total_x_11 += __builtin_popcountll(isx_11); // both parents are het, none missing.

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
  pedigree_stats->d = (ND) { n_0xorx0_2_nomd + n_2xorx2_0_nomd
			     + n_00_2_22_0 // to double count 00_2 and 22_0
			     + F*n_00_1_22_1,
			     n_total_no_md 
			     //  - n_total_x_11   // count these in denominator ??
  };
  pedigree_stats->z = (ND) {n_00_1_22_1, n_00_22};
  pedigree_stats->par1_hgmr = (ND) {hgmr1_numerator, hgmr1_denominator};
  pedigree_stats->par2_hgmr = (ND) {hgmr2_numerator, hgmr2_denominator};
  pedigree_stats->par1_R = (ND) {n_0x_1_2x_1, n_0x_1_2x_1 + n_0x_0_2x_2};
  pedigree_stats->par2_R = (ND) {n_x0_1_x2_1, n_x0_1_x2_1 + n_x0_0_x2_2};
  pedigree_stats->n_01or10_1 = n_01or10_1;
  pedigree_stats->all_good_count = n_total_no_md;

  //fprintf(stderr, "XXXX: %ld  %ld   %ld\n", n_total_no_md, n_total_x_11, n_total_no_md-n_total_x_11);
  // pedigree_stats->d_n = n_over_d(pedigree_stats->d)/the_gtset->mean_d;
  
  //assert(pedigree_stats->par1_R.n == Rnd.n);
  //assert(pedigree_stats->par1_R.d == Rnd.d);
  return pedigree_stats;
} // end of bitwise_triple_counts



Vpedigree*  calculate_triples_for_one_accession(Accession* prog, const GenotypesSet* the_genotypes_set, Viaxh* cppps, long max_candidate_parents){

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
      // free(the_pedigree->pedigree_stats); 
      the_ps = calculate_pedigree_stats(the_pedigree, the_genotypes_set);
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
  //the_ps->par1_xhgmr = (ND) {0, 0}; 
  the_ps->par1_R = (ND) {0, 0};
  
  the_ps->par2_hgmr = (ND) {0, 0}; 
  //the_ps->par2_xhgmr = (ND) {0, 0};
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

double calculate_xFTR(Pedigree* the_pedigree, const GenotypesSet* the_gtsset){
  long N = 0; // counts forbidden triples
  double D = 0; //  
  if(the_pedigree->A == NULL  || the_pedigree->F == NULL  ||  the_pedigree->M == NULL){
    return -1;
  }else{
    // Vlong** marker_dosage_counts; // counts of dosages for each marker. marker_dosage_counts->[i]->a[j] is count of dosage i for marker j
    Vchar* Agts = the_pedigree->A->genotypes;
    Vchar* Fgts = the_pedigree->F->genotypes;
    Vchar* Mgts = the_pedigree->M->genotypes;
    for(long i=0; i<Agts->length; i++){
      double ok_count = (double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
      char Fgt = Fgts->a[i];
      char Mgt = Mgts->a[i];
      char Agt = Agts->a[i];
      if(Fgt == '0'){
	if(Mgt == '0'){	  
	  D += (the_gtsset->marker_dosage_counts[1]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/ok_count;
	  if(Agt == '1' || Agt == '2') N++;
	}else if(Mgt == '1'){
	  D += (the_gtsset->marker_dosage_counts[2]->a[i])/ok_count;
	  if(Agt == '2') N++;
	}else if(Mgt == '2'){
	  D += (the_gtsset->marker_dosage_counts[0]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/ok_count;
	  if(Agt == '0' || Agt == '2') N++;
	}
      }else if(Fgt == '1'){
	if(Mgt == '0'){	  
	  D += (the_gtsset->marker_dosage_counts[2]->a[i])/ok_count;
	  if(Agt == '2') N++;
	}else if(Mgt == '2'){
	  D += (the_gtsset->marker_dosage_counts[0]->a[i])/ok_count;
	  if(Agt == '0') N++;
	}
      }else if(Fgt == '2'){
	if(Mgt == '0'){	  
	  D += (the_gtsset->marker_dosage_counts[0]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/ok_count;
	  if(Agt == '0' || Agt == '2') N++;
	}else if(Mgt == '1'){
	  D += (the_gtsset->marker_dosage_counts[0]->a[i])/ok_count;
	  if(Agt == '0') N++;
	}else if(Mgt == '2'){
	  if(Agt == '1' || Agt == '0') N++;
	  D += (the_gtsset->marker_dosage_counts[1]->a[i] +  the_gtsset->marker_dosage_counts[0]->a[i])/ok_count;
	}  
      }
    }
    return (D>0)? N/D : -1;
  }
}

double calculate_xxFTR(Pedigree* the_pedigree, const GenotypesSet* the_gtsset, double alpha){
  long N = 0; // counts forbidden triples
  double d = 0; // expected number of forbidden triples based on proportions of 0, 1, 2 in this accession
  double D = 0; // expected number of forbidden triples based proportions of 0, 1, 2 in each marker.
  if(the_pedigree->A == NULL  || the_pedigree->F == NULL  ||  the_pedigree->M == NULL){
    return -1;
  }else{
    // Vlong** marker_dosage_counts; // counts of dosages for each marker. marker_dosage_counts->[i]->a[j] is count of dosage i for marker j
    Accession* A = the_pedigree->A;
    Accession* F = the_pedigree->F;
    Accession* M = the_pedigree->M;
    Vchar* Agts = A->genotypes;
    Vchar* Fgts = F->genotypes;
    Vchar* Mgts = M->genotypes;
    double acc_ok_count_x = (double)(A->dosage_counts[0] + A->dosage_counts[1] + A->dosage_counts[2]); // count of missing gts in this accessions,
    long A0 = 0; long A1 = 0; long A2 = 0;
    for(long i=0; i<Agts->length; i++){
      char Fgt = Fgts->a[i];
      char Mgt = Mgts->a[i];
      char Agt = Agts->a[i];
      if(Fgt == 'X' || Mgt == 'X' || Agt == 'X') continue;
      if(Fgt == '1'  &&  Mgt == '1') continue;
      if(Agt == '0'){ A0++; }
      else if(Agt == '1'){ A1++; }
      else if(Agt == '2'){ A2++; }
    }
    double acc_ok_count = A0 + A1 + A2;
    long ok_triple_count = 0;
    fprintf(stderr, "fdsa: %7.5f %7.5f \n", acc_ok_count_x, acc_ok_count);
    // but should exclude missing gts in parents also, i.e. only count markers with valid gts in all 3.
    for(long i=0; i<Agts->length; i++){
      double marker_ok_count =   // the number of accessions with non-missing genotypes for this marker
	(double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
      char Fgt = Fgts->a[i];
      char Mgt = Mgts->a[i];
      char Agt = Agts->a[i];
  
      if(Fgt == '0'){
	if(Mgt == '0'){ // Agt = 1, 2 forbidden
	  // fprintf(stderr, "ZZZ: %ld  %ld  %ld \n", acc_ok_count, A->dosage_counts[1], A->dosage_counts[2]);
	  d += A1 + A2; // /acc_ok_count;
	  D += (the_gtsset->marker_dosage_counts[1]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '1' || Agt == '2') N++; // 
	}else if(Mgt == '1'){ // 2 forbidden
	  d += A2; // /acc_ok_count;
	  D += (the_gtsset->marker_dosage_counts[2]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '2') N++;
	}else if(Mgt == '2'){ // 0, 2 forbidden
	  d += A0 + A2; // /acc_ok_count;
	  D += (the_gtsset->marker_dosage_counts[0]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '0' || Agt == '2') N++;
	}
      }else if(Fgt == '1'){
	if(Mgt == '0'){	// 2 forbidden
	  d += A2;
	  D += (the_gtsset->marker_dosage_counts[2]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '2') N++;
	}else if(Mgt == '1'){
	  
	}else if(Mgt == '2'){ // 0 forbidden
	  d += A0;
	  D += (the_gtsset->marker_dosage_counts[0]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '0') N++;
	}
      }else if(Fgt == '2'){
	if(Mgt == '0'){	// 0, 2 forbidden  
	  d += A0 + A2;
	  D += (the_gtsset->marker_dosage_counts[0]->a[i] + the_gtsset->marker_dosage_counts[2]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '0' || Agt == '2') N++;
	}else if(Mgt == '1'){ // 0 forbidden
	  d += A0;
	  D += (the_gtsset->marker_dosage_counts[0]->a[i])/marker_ok_count;
	  ok_triple_count++;
	  if(Agt == '0') N++;
	}else if(Mgt == '2'){ // 0, 1 forbidden
	  if(Agt == '1' || Agt == '0') N++;
	  d += A0 + A1;
	  D += (the_gtsset->marker_dosage_counts[1]->a[i] +  the_gtsset->marker_dosage_counts[0]->a[i])/marker_ok_count;
	  ok_triple_count++;
	}  
      }
    }
    d /= acc_ok_count;
    fprintf(stderr, "asdf  %ld  %7.5f  %7.5f  %ld\n", N, d, D, ok_triple_count);
    return (D>0 && d>0)? alpha*(N/d) + (1.0-alpha)*(N/D) : -1;
  }
}




Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, const GenotypesSet* the_gtsset){
  long ploidy = the_gtsset->ploidy;
  Pedigree_stats* the_ps; //  = construct_pedigree_stats(); // (Pedigree_stats*)calloc(1, sizeof(Pedigree_stats));
  assert(the_pedigree->F != NULL  ||  the_pedigree->M != NULL); // shouldn't have both parents NULL //
  if(the_pedigree->F != NULL  &&  the_pedigree->M != NULL){
   
    the_ps = bitwise_triple_counts(the_pedigree->F,  the_pedigree->M,  the_pedigree->A);

    /*   if(0){ // compare bitwise, nonbitwise calculations as check - slow
      Pedigree_stats* nobw_ps = triple_counts( the_pedigree->F->genotypes->a,  the_pedigree->M->genotypes->a,  the_pedigree->A->genotypes->a, ploidy );
      assert
	(NDs_equal(the_ps->agmr12, nobw_ps->agmr12));
      assert(NDs_equal(the_ps->par1_hgmr, nobw_ps->par1_hgmr));
      assert(NDs_equal(the_ps->par1_R, nobw_ps->par1_R));
      assert(NDs_equal(the_ps->par2_hgmr, nobw_ps->par2_hgmr));
      assert(NDs_equal(the_ps->par2_R, nobw_ps->par2_R));
      assert(NDs_equal(the_ps->d, nobw_ps->d));
      assert(NDs_equal(the_ps->z, nobw_ps->z));
    } /* */
       
    //  the_ps->par1_xhgmr = xhgmr(the_gtsset, the_pedigree->F, the_pedigree->A, false);
    //  the_ps->par2_xhgmr = xhgmr(the_gtsset, the_pedigree->M, the_pedigree->A, false);
  }else{ // one of the parents is NULL
    the_ps = construct_pedigree_stats();
    the_ps->agmr12 = (ND) {0, 0};
    the_ps->z = (ND) {0, 0};
    if(the_pedigree->F != NULL){ // we have female parent id, no male parent id
      four_longs hgmrR = hgmr_R(the_pedigree->F->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
   
      the_ps->par1_hgmr.n = hgmrR.l1;
      the_ps->par1_hgmr.d = hgmrR.l2;
      the_ps->par1_R.n = hgmrR.l3;
      the_ps->par1_R.d = hgmrR.l4;
      the_ps->par2_hgmr = (ND) {0, 0};
      the_ps->par2_R = (ND) {0, 0};
      //  the_ps->par1_xhgmr = xhgmr(the_gtsset, the_pedigree->F, the_pedigree->A, false);
      //  the_ps->par2_xhgmr = (ND) {0, 0};
    }else{ // we have male parent id, no female parent id
      if(DO_ASSERT) assert(the_pedigree->M != NULL);
      the_ps->par1_hgmr = (ND) {0, 0};
      the_ps->par1_R = (ND) {0, 0};
      four_longs hgmrR = hgmr_R(the_pedigree->M->genotypes->a, the_pedigree->A->genotypes->a, (char)(ploidy + 48));
      the_ps->par2_hgmr.n = hgmrR.l1;
      the_ps->par2_hgmr.d = hgmrR.l2;
      the_ps->par2_R.n = hgmrR.l3;
      the_ps->par2_R.d = hgmrR.l4;
      //  the_ps->par1_xhgmr = (ND) {0, 0};
      //   the_ps->par2_xhgmr = xhgmr(the_gtsset, the_pedigree->M, the_pedigree->A, false);  
    }   
  }
  the_ps->hgmr1_n = n_over_d(the_ps->par1_hgmr)/the_gtsset->mean_hgmr;
  the_ps->hgmr2_n = n_over_d(the_ps->par2_hgmr)/the_gtsset->mean_hgmr;
  the_ps->R1_n = n_over_d(the_ps->par1_R)/the_gtsset->mean_R;
  the_ps->R2_n = n_over_d(the_ps->par2_R)/the_gtsset->mean_R;
  the_ps->d_n = n_over_d(the_ps->d)/the_gtsset->mean_d;
  the_ps->z_n = n_over_d(the_ps->z)/the_gtsset->mean_z;
  the_ps->ftc_n = the_ps->d.n/the_gtsset->mean_ftc;
  
  //  the_ps->xhgmr1 = n_over_d(the_ps->par1_xhgmr);
  //  the_ps->xhgmr2 = n_over_d(the_ps->par2_xhgmr);
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
    //fprintf(fh, " %ld %ld ", the_pedigree_stats->d.n, the_pedigree_stats->d.d);
    //print_d_r(fh, the_pedigree_stats->d_old);

  }else{ // print ratios but not numerators and denominators
    print_n_over_d(fh, the_pedigree_stats->agmr12);
    print_n_over_d(fh, the_pedigree_stats->par1_hgmr);
    // print_n_over_d(fh, the_pedigree_stats->par1_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par1_R);
    print_n_over_d(fh, the_pedigree_stats->par2_hgmr);
    //  print_n_over_d(fh, the_pedigree_stats->par2_xhgmr);
    print_n_over_d(fh, the_pedigree_stats->par2_R);
    print_n_over_d(fh, the_pedigree_stats->d);
    // print_n_over_d(fh, the_pedigree_stats->z);
    // print_n_over_d(fh, the_pedigree_stats->d_old);
  }	     	   
}

void print_pedigree_normalized(FILE* fh, Pedigree* the_pedigree){
  //double mean_hgmr, double mean_R, double mean_d, double mean_z){
  Accession* F = the_pedigree->F;
  Accession* M = the_pedigree->M;
  fprintf(fh, "\t%s\t%s\t%ld", (F != NULL)? F->id->a : "NA", (M != NULL)? M->id->a : "NA",
	  the_pedigree->pedigree_stats->all_good_count);
  print_normalized_pedigree_stats(fh, the_pedigree->pedigree_stats);
}

void print_normalized_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){
  print_n_over_d(fh, the_pedigree_stats->agmr12);

  print_double_nan_as_hyphen(fh, the_pedigree_stats->hgmr1_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->R1_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->hgmr2_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->R2_n);
  print_double_nan_as_hyphen(fh, the_pedigree_stats->d_n);
  // print_double_nan_as_hyphen(fh, the_pedigree_stats->z_n);
}

void print_double_nan_as_hyphen(FILE* fh, double x){
  if(isnan(x)){
    fprintf(fh, "\t-");
  }else{
    fprintf(fh, "\t%7.5lf", x);
  }
}

bool d_ok(Pedigree_stats* p, double max_ok_d){ // true: d looks good; false: d too large.
  double d = n_over_d(p->d);
  return (d <= max_ok_d)? true : false;
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
 
    if(strcmp(acc_id, "NA") != 0){ // id of this accession is not "NA" 
      acc_idx = index_of_id_in_vidxid(the_gt_vidxid, acc_id);
      Accession* Acc = the_gtsset->accessions->a[acc_idx]; // id->index;
      if(acc_idx == ID_NA_INDEX){ // no genotypes for this accession
	not_in_genotypes_set_count++;
      }else{ // have genotypes for this accession

	fempar_idx = index_of_id_in_vidxid(the_gt_vidxid, fempar_id);
	malpar_idx = index_of_id_in_vidxid(the_gt_vidxid, malpar_id);	 
	if( (fempar_idx != ID_NA_INDEX) || (malpar_idx != ID_NA_INDEX) ){ // pedigree file has at least 1 parent for which there are genotypes

	  /* bool gt1pedigree = false; */
	  /* if(Acc->has_pedigree){ */
	  /*   gt1pedigree = true; */
	  /*   long fidx = Acc->Fpar_idx;	     */
	  /*   char* fid = (fidx == ID_NA_INDEX)? "NA": the_gtsset->accessions->a[fidx]->id->a; */
	  /*   long midx = Acc->Mpar_idx;	     */
	  /*   char* mid = (midx == ID_NA_INDEX)? "NA" : the_gtsset->accessions->a[midx]->id->a; */
	  /*   fprintf(stderr, "offspring %s  Fpar %s  Mpar %s \n", Acc->id->a, fid, mid);  */
	  /* } */
	  Acc->has_pedigree = true;
	  Acc->Fpar_idx = fempar_idx;
	  Acc->Mpar_idx = malpar_idx;
	  /* if(gt1pedigree){ */
	  /*   long fidx = Acc->Fpar_idx;	     */
	  /*   char* fid = (fidx == ID_NA_INDEX)? "NA": the_gtsset->accessions->a[fidx]->id->a; */
	  /*   long midx = Acc->Mpar_idx;	     */
	  /*   char* mid = (midx == ID_NA_INDEX)? "NA" : the_gtsset->accessions->a[midx]->id->a; */
	  /*   fprintf(stderr, "   offspring %s  Fpar %s  Mpar %s \n\n", Acc->id->a, fid, mid); */
	  /* } */
	  Accession* Fpar = (fempar_idx != ID_NA_INDEX)?
	    the_gtsset->accessions->a[fempar_idx] : NULL; // id->index;
	  Accession* Mpar = (malpar_idx != ID_NA_INDEX)?
	    the_gtsset->accessions->a[malpar_idx] : NULL; // id->index;
	  Pedigree* a_pedigree = construct_pedigree(Acc, Fpar, Mpar);
	  // Acc->pedigree = a_pedigree;
	  push_to_vpedigree(pedigrees, a_pedigree); // store pedigrees (with genotypes for accession and at least 1 of the parents)
	}else{
	  no_parent_gts_count++;
	  // fprintf(stdout, "Invalid pedigree for accession: %s ; both parents are NA or lack genotypes.\n", Acc->id->a);  
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

void sort_vpedigree_by_d(Vpedigree* the_vped){ 
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

void free_vpedigree(const Vpedigree* the_vped){
  if(the_vped == NULL) return;
  for(long i=0; i<the_vped->size; i++){
    free_pedigree(the_vped->a[i]);
  }
  free(the_vped->a);
  free((Vpedigree*)the_vped);
}

//  ************ unused *************************

/* // *****  sorting an array of Idxhgmr  ***** 
int cmpidxhgmr(const void* v1, const void* v2){
  const Idxhgmr* s1 = (const Idxhgmr*)v1;
  const Idxhgmr* s2 = (const Idxhgmr*)v2;
  return (s1->hgmr < s2->hgmr)? -1 : (s1->hgmr > s2->hgmr)? 1 : 0;
}

void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array){ // sort in place
  qsort(array, size, sizeof(Idxhgmr), cmpidxhgmr);
} /* */


/*
void print_pedigree_alternatives(FILE* fh, const Vpedigree* alt_pedigrees, long max_to_print, bool verbose){
  fprintf(fh, " %3ld  ", alt_pedigrees->size);
  long n_to_print = (max_to_print < alt_pedigrees->size)? max_to_print : alt_pedigrees->size;
  for(long i=0; i<n_to_print; i++){
    Pedigree* alt_pedigree = alt_pedigrees->a[i];
    fprintf(fh, "%20s %20s %ld ", alt_pedigree->F->id->a, alt_pedigree->M->id->a, alt_pedigree->pedigree_stats->all_good_count);
    print_pedigree_stats(fh, alt_pedigree->pedigree_stats, verbose);
  } 
} /* */


/* 
Vpedigree* pedigree_alternatives(const Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, const Vlong* parent_idxs, double max_ok_hgmr, double max_ok_d){
  
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
  // first, include the accessions in the pedigree
  if(fparent != NULL) push_to_vlong(best_parent_candidate_idxs, fparent_idx); // add female parent (from pedigree)
  if(mparent != NULL  &&  mparent_idx != fparent_idx) push_to_vlong(best_parent_candidate_idxs, mparent_idx); // add male parent (from pedigree) if distinct

  //fprintf(stderr, "size of best_parent_candidate_idxs: %ld \n", best_parent_candidate_idxs->size);

  // get hgmrs w.r.t. the other accessions (i.e. besides fparent and mparent in pedigree) 
  Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_parents*sizeof(Idxhgmr));
  for(long i=0; i<parent_idxs->size; i++){
    long idx = parent_idxs->a[i];
    //fprintf(stderr, "Z %ld %ld %ld\n", idx, fparent_idx, mparent_idx);
    //  if(idx > 10000) {fprintf(stderr, "bad idx: %ld\n", idx); getchar(); }
    the_idxhgmrs[i].idx = idx;
    if(idx == fparent_idx  ||  idx == mparent_idx){
      the_idxhgmrs[i].hgmr = 1000.0; // big - fparent_idx, mparent_idx are already in best_candidate_idxs
    }else{
      char* pgts = the_gtsset->accessions->a[idx]->genotypes->a; // genotype_sets->a[idx];
      //the_idxhgmrs[i].idx = idx;
      four_longs bwah = bitwise_agmr_hgmr(the_pedigree->A, the_gtsset->accessions->a[idx]);
      //  long b_agmr_num = bfcs.l1 + bfcs.l3;
      //	  long b_agmr_denom = b_agmr_num + bfcs.l2 + bfcs.l4;
      long b_hgmr_num = bwah.l1;
      long b_hgmr_denom = bwah.l1 + bwah.l2;
      //	  agmr = (b_agmr_denom > 0)? (double)b_agmr_num / (double)b_agmr_denom : -1;
      double b_hgmr = (b_hgmr_denom > 0)? (double)b_hgmr_num / (double)b_hgmr_denom : -1;
      //fprintf(stderr, "EDFG: %7.4f  %7.4f  %7.4f  %7.4f \n", b_hgmr, the_gtsset->mean_hgmr, max_ok_hgmr,  b_hgmr/the_gtsset->mean_hgmr);
      the_idxhgmrs[i].hgmr = b_hgmr/the_gtsset->mean_hgmr; // (acc_gts, pgts); //the_hgmr; //the_idxhgmr = {idx, the_hgmr);
    
      //fprintf(stderr, "hgmrs: %7.5f  \n", the_idxhgmrs[i].hgmr);
    }
  }
  // sort accession indices by hgmr
  sort_idxhgmr_by_hgmr(n_parents, the_idxhgmrs);
  for(long i=0; i<n_parents; i++){ // store 'good' hgmr indices 
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != acc_idx){
      double the_hgmr = the_idxhgmrs[i].hgmr;
      fprintf(stderr, " DFDF: %7.4f  %7.4f\n", the_hgmr, max_ok_hgmr);
      if(the_hgmr >= max_ok_hgmr) break; // all the rest are worse, so skip them.
      if(the_idx != fparent_idx  &&  the_idx != mparent_idx){
	push_to_vlong(best_parent_candidate_idxs, the_idx);
      }     
    }
  }
  free(the_idxhgmrs);
  
  long ub = long_min(best_parent_candidate_idxs->size, 80); // set the number of possible parents to consider.
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
	Pedigree* alt_pedigree = construct_pedigree(the_pedigree->A, acc1, acc2); // arbitrarily acc1 is Fem parent, acc2 is male
	Pedigree_stats* alt_pedigree_stats = calculate_pedigree_stats(alt_pedigree, the_gtsset);
	  alt_pedigree->pedigree_stats = alt_pedigree_stats;
	  push_to_vpedigree(alt_pedigrees, alt_pedigree);
      }
    } // j loop
  } // i loop
  the_pedigree->A->search_done = true;
  sort_vpedigree_by_d(alt_pedigrees);
  free_vlong(best_parent_candidate_idxs);
  return alt_pedigrees;
} /* */


/*
three_longs count_FT2_one_chromosome(const GenotypesSet* the_gtsset, Accession* parent1, Accession* parent2,
				  Accession* offspring, long first, long next){
  // count FT2 triples, i.e. those  have dosage 1 in the offspring
  // and may be forbidden, depending on the phase of the offspring (e.g. 01_1, 21_1, 02_1)
  // these are very sensitive to mistakes in the phase in the offspring
  // if there is a switch of phase somewhere in the middle of the offspring
  // (e.g. phase should be + everywhere but is + for first half, and - for second half)
  // then instead of comparing parent 1 with the two chromosomes (call them a and b) of each homologous pair
  // (and similarly for parent 2), we are, in effect, comparing parent 1 with something that is
  // partially a and partially b. I suspect this often happens because the two chromosome arms
  // (separated by centromeric region with low density of markers) do not have the correct
  // relative phase. It would probably be good to analyze the chromosome arms separately
  // As an alternative to that this function counts FT2 for 1->a,2->b  and 1->b,2->a as a function of
  // the index of the triples being considered, and assumes a rel phase switch position in the offspring
  // and minimizes FT2 over possible phase switch positions, and over 1->a,2->b,  1->b,2->a.


  if(parent1 == NULL  ||  parent2 == NULL  ||  offspring == NULL) {
    return (three_longs){-1, -1, -1};
  }

  Vlong* p1ap2b_FT2_counts = construct_vlong(1000);
  Vlong* FT2_indices = construct_vlong(1000);
  push_to_vlong(p1ap2b_FT2_counts, 0);
  push_to_vlong(FT2_indices, -1);
  Vlong* p1bp2a_FT2_counts = construct_vlong(1000);
  Vlong* p1bp2a_FT2_indices = construct_vlong(1000);
  push_to_vlong(p1bp2a_FT2_counts, 0);
  // push_to_vlong(p1bp2a_FT2_indices, -1);
  long FT2_index = 0;
  long nophase_count = 0;
   
  long chrom_number = the_gtsset->chromosomes->a[first];
  // e.g. Fa count markers with o_gt = 0 (ref,ref) and a has alt allele, and o_gt = 2(alt,alt) and a has ref allele. 
  for(long i=first; i < next; i++){
   

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt != '1') continue;
    char o_phase = offspring->phases->a[i];
    
    char p1_gt = parent1->genotypes->a[i];
    if(p1_gt == MISSING_DATA_CHAR) continue;
    //char p_phase = parent1->phases->a[i];
  
    char p2_gt = parent2->genotypes->a[i];
    if(p2_gt == MISSING_DATA_CHAR) continue; 

    char offA = 'x', offB = 'x';
    if(o_phase == 'p'){ // p -> a has ref allele, b has alt allele.
      offA = 'r'; offB = 'a';
    }else if(o_phase == 'm'){
      offA = 'a'; offB = 'r';
    }else{
      nophase_count++;
      continue;
    }
    long ab = 0, ba = 0;
    if(p1_gt == '0'){
      if(p2_gt == '0'){ // 00_1
	continue;
      }else if(p2_gt == '1'){ // 01_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 1; // p1->b: 1, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  ab = 1; // p1->a: 1, p2->b: 0
	  ba = 0; // p1->b: 0, p2->a: 0
	}
      }else{ // p2_gt == '2'   02_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 2; // p1->b: 1, p2->a: 1
	}else{ // offA == 'a', offB = 'r'
	  ab = 2; // p1->a: 1, p2->b: 1
	  ba = 0; // p1->b: 0, p2->a: 0
	}
      }
    }else if(p1_gt == '1'){
      if(p2_gt == '0'){ // 10_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 1; // p1->a: 0, p2->b: 1
	  ba = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 1; // p1->b: 0, p2->a: 1
	}
      }else if(p2_gt == '2'){ // 12_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 1; // p1->b: 0, p2->a: 1
	}else{ // offA == 'a', offB = 'r'
	  ab = 1; // p1->a: 0, p2->b: 1
	  ba = 0; // p1->b: 0, p2->a: 0
	}
      }else{ // 11_1
	continue; 
      }
    }else if(p1_gt == '2'){
      if(p2_gt == '1'){ // 21_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 1; // p1->a: 1, p2->b: 0
	  ba = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 1; // p1->b: 1, p2->a: 0
	}
      }else if(p2_gt == '0'){ // 20_1
	if(offA == 'r'){ // and offB == 'a'
	  ab = 2; // p1->a: 1, p2->b: 1
	  ba = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  ab = 0; // p1->a: 0, p2->b: 0
	  ba = 2; // p1->b: 1, p2->a: 1
	}
      }else{ // 22_1
	continue;
      }
    }
           
    push_to_vlong(p1ap2b_FT2_counts, p1ap2b_FT2_counts->a[FT2_index] + ab);
    push_to_vlong(FT2_indices, i);
    push_to_vlong(p1bp2a_FT2_counts, p1bp2a_FT2_counts->a[FT2_index] + ba);
    // push_to_vlong(p1bp2a_FT2_indices, i);
    FT2_index++;
  } // end of loop over markers on chromosome


   long ab_ba_min_FT2 = 10000000; // will be best (min) FT2 when switching from p1bp2a to p1ap2b
  long ab_ba_min_idx = 0; // range: 0 to number of FT2 markers; best place to switch from p1bp2a to p1ap2b
  long ba_ab_min_FT2 = 10000000; // will be best (min) FT2 when switching from p1bp2a to p1ap2b
  long ba_ab_min_idx = 0; // range: 0 to number of FT2 markers; best place to switch from p1bp2a to p1ap2b
  long N_FT2_markers = p1ap2b_FT2_counts->size-1; // number of FT2 ('type 2' forbidden triples) markers on chromosome
  long N_markers = next - first; // number of markers on chromosome
  long ab_total = p1ap2b_FT2_counts->a[N_FT2_markers]; // FT2 count when p1ap2b for whole chromosome
  long ba_total = p1bp2a_FT2_counts->a[N_FT2_markers]; // FT2 count when p1bp2a for whole chromosome
  fprintf(stderr, "sizes: %ld  %ld  %ld\n", FT2_index, N_FT2_markers, N_markers);
  fprintf(stderr, "ab, ba FT2 totals: %ld %ld \n", ab_total, ba_total);

  long Gab_ba, Gba_ab;
  for(long i=0; i<p1ap2b_FT2_counts->size; i++){
    // fprintf(stderr, "chrom: %ld  k, Gab(k), Gba(k): %ld %ld %ld \n",
    //   chrom_number, i, p1ap2b_FT2_counts->a[i], p1bp2a_FT2_counts->a[i]); 
    Gab_ba = p1ap2b_FT2_counts->a[i] + (ba_total - p1bp2a_FT2_counts->a[i]); // first part 1a2b, 2nd part 1b2a
    if(Gab_ba < ab_ba_min_FT2){
      ab_ba_min_FT2 = Gab_ba; ab_ba_min_idx = i;
    }
    Gba_ab = p1bp2a_FT2_counts->a[i] + (ab_total - p1ap2b_FT2_counts->a[i]); // first part 1b2a, 2nd part 1a2b
    if(Gba_ab < ba_ab_min_FT2){
      ba_ab_min_FT2 = Gba_ab; ba_ab_min_idx = i;
    }
    // fprintf(stderr, "ABVF: chrom: %ld   %ld   %ld %ld \n", chrom_number, Gab_ba, Gba_ab, ab_total+ba_total - Gab_ba); 
  }
  fprintf(stderr, "ZBAVL chrom: %ld   %ld %ld  FT2s un:  %ld %ld  %ld corrected: %ld %ld   FT2s  %ld %ld  %ld  1st,nxt: %ld %ld    %ld  %ld \n", //  %s %s %s\n",
	  chrom_number,
	  N_markers, FT2_index,

	  ab_total, ba_total,
	  (ab_total < ba_total)? ab_total : ba_total,
	  
	  ab_ba_min_idx, ba_ab_min_idx,
	  ab_ba_min_FT2, ba_ab_min_FT2,
	  (ab_ba_min_FT2 < ba_ab_min_FT2)? ab_ba_min_FT2 : ba_ab_min_FT2,
	  
	  first, next, 
	  FT2_indices->a[ab_ba_min_idx], FT2_indices->a[ba_ab_min_idx]
	
	  // parent1->id->a, parent2->id->a, offspring->id->a
	  );

  long FT2_count_uncorrected = (ab_total < ba_total)? ab_total : ba_total;
  long FT2_count_corrected = (ab_ba_min_FT2 < ba_ab_min_FT2)? ab_ba_min_FT2 : ba_ab_min_FT2;
  long jlo, max_min_idx;
  if(ab_ba_min_FT2 < ba_ab_min_FT2){
    fprintf(stderr, "ab_ba < ba_ab: %ld %ld \n", ab_ba_min_FT2, ba_ab_min_FT2);
    jlo = FT2_indices->a[ab_ba_min_idx];
    fprintf(stderr, "ab_ba < ba_ab: %ld %ld  %ld\n", ab_ba_min_FT2, ba_ab_min_FT2, jlo);
  }else if(ba_ab_min_FT2 < ab_ba_min_FT2){
    fprintf(stderr, "ab_ba > ba_ab: %ld %ld \n", ab_ba_min_FT2, ba_ab_min_FT2);
    jlo = FT2_indices->a[ba_ab_min_idx];
    fprintf(stderr, "ab_ba > ba_ab: %ld %ld  %ld\n", ab_ba_min_FT2, ba_ab_min_FT2, jlo);
  }else{ // equal
    long jlo_abba = FT2_indices->a[ab_ba_min_idx];
    long jlo_baab = FT2_indices->a[ba_ab_min_idx];
    if(ab_ba_min_idx >= ba_ab_min_idx){
      jlo = jlo_abba; max_min_idx = ab_ba_min_idx;
    }else{
      jlo = jlo_baab; max_min_idx = ba_ab_min_idx;
    }
    fprintf(stderr, "jlo_abba, jlo_baab, jlo: %ld %ld %ld   %ld %ld\n",
	    jlo_abba, jlo_baab, jlo, 
	    max_min_idx, FT2_indices->size);
  }
  fprintf(stderr, "max_min_idx, FT2_indices->size: %ld %ld \n", max_min_idx, FT2_indices->size);
  if(max_min_idx == FT2_indices->size - 1){ // no correction needed
    for(long ii = first; ii < next; ii++){
      offspring->corrected_phases->a[ii] = offspring->phases->a[ii];
      // fprintf(stderr, "Aph, corrph: %c  %c \n", offspring->phases->a[ii], offspring->corrected_phases->a[ii]);
    }
  }else{ // flip the phases for later part of chromosome
    long jhi = FT2_indices->a[max_min_idx + 1];

    long switch_index = (long)(0.5*(jlo + jhi)); // half way between
    fprintf(stderr, "switch_index: %ld\n", switch_index);
    for(long ii = first; ii < next; ii++){
      char ph = offspring->phases->a[ii];
      if(ii >= switch_index){
	offspring->corrected_phases->a[ii] = (ph == 'p')? 'm' : ((ph == 'm')? 'p' : ph); // p<->m, x unchanged
      }else{
	offspring->corrected_phases->a[ii] = offspring->phases->a[ii];
      }
      // fprintf(stderr, "Bph, corrph: %c  %c \n", ph, offspring->corrected_phases->a[ii]);
    }
  }
  fprintf(stderr, "end of count_FT2_one_chrom\n");
  return (three_longs){FT2_count_uncorrected, FT2_count_corrected, N_FT2_markers};
  } // end of count_FT2_one_chromosome  /* */

/* 
long pedigree_ok(Pedigree_stats* p, double max_self_agmr12, double max_self_r, double max_ok_d){
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
    if((r1 <= max_self_r) && (r2 <= max_self_r) && (d <= max_ok_d) ){
      result = 1;
    }
  }else{ // parents in pedigree are not very similar
    if((r1 > max_self_r) && (r2 > max_self_r) && (d <= max_ok_d)){
      result = 2;
    }else if( (r1 <= max_self_r) && (d > max_ok_d) ){
      result = -1;
    }else if( (r2 <= max_self_r) && (d > max_ok_d) ){
      result = -2;
    }
  }
  return result;
} /* */

/* void FT2_phase_correction(const GenotypesSet* the_gtsset, Accession* parent1, Accession* parent2,
				  Accession* offspring, long first, long next){

  // do error correction on the phases of the offspring accession
  // by considering 'type 2' forbidden triples
  // these are those which are forbidden or not depending on the phase of the offspring
  // i.e.: 01_1, 10_1, 21_1, 12_1, 02_1, 20_1
  // the idea is to see if reversing the phase of one end of the chromosome
  // (i.e. beyond some point to be determined)
  // will reduce the number of these type 2 forbidden triples (FT2's)
  // to a small value.
  // k indexes the FT2 markers
  // i indexes all the markers
  
  if(parent1 == NULL  ||  parent2 == NULL  ||  offspring == NULL) {
    return; //  (three_longs){-1, -1, -1};
  }
 
  Vlong* ft2_p1a_k = construct_vlong(1000); // count of p1->a forbidden ft2 markers 1 through k
  Vlong* ft2_p1b_k = construct_vlong(1000);
  Vlong* ft2_p2a_k = construct_vlong(1000);
  Vlong* ft2_p2b_k = construct_vlong(1000);
  push_to_vlong(ft2_p1a_k, 0);
  push_to_vlong(ft2_p1b_k, 0);
  push_to_vlong(ft2_p2a_k, 0);
  push_to_vlong(ft2_p2b_k, 0);
  Vlong* marker_index_of_ft2_k = construct_vlong(1000);  
  push_to_vlong(marker_index_of_ft2_k, -1);

  long nophase_count = 0;

  long chrom_number = the_gtsset->chromosomes->a[first];
  // e.g. Fa count markers with o_gt = 0 (ref,ref) and a has alt allele, and o_gt = 2(alt,alt) and a has ref allele.
  long k = 0; // index of ft2 markers on chromosome
  for(long i=first; i < next; i++){
    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt != '1') continue;
    char o_phase = offspring->phases->a[i];
    
    char p1_gt = parent1->genotypes->a[i];
    if(p1_gt == MISSING_DATA_CHAR) continue;
    //char p_phase = parent1->phases->a[i];
  
    char p2_gt = parent2->genotypes->a[i];
    if(p2_gt == MISSING_DATA_CHAR) continue; 

    char offA = 'x', offB = 'x';
    if(o_phase == 'p'){ // a has ref allele, b has alt allele.
      offA = 'r'; offB = 'a';
    }else if(o_phase == 'm'){
      offA = 'a'; offB = 'r';
    }else{ // no phase for offspring 
      nophase_count++;
      continue;
    }

    long p1a = 0, p1b = 0, p2a = 0, p2b = 0;
    if(p1_gt == '0'){
      if(p2_gt == '0'){ // 00_1
	continue;
      }else if(p2_gt == '1'){ // 01_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 0; p2b = 0;  // p1->a: 0, p2->b: 0
	  p1b = 1; p2a = 0; // p1->b: 1, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  p1a = 1; p2b = 0; // p1->a: 1, p2->b: 0
	  p1b = 0; p2a = 0;  // p1->b: 0, p2->a: 0
	}
      }else{ // p2_gt == '2'   02_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 0; p2b = 0; // p1->a: 0, p2->b: 0
	  p1b = 1; p2a = 1; // p1->b: 1, p2->a: 1
	}else{ // offA == 'a', offB = 'r'
	  p1a = 1; p2b = 1; // p1->a: 1, p2->b: 1
	  p1b = 0; p2a = 0; // p1->b: 0, p2->a: 0
	}
      }
    }else if(p1_gt == '1'){
      if(p2_gt == '0'){ // 10_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 0; p2b = 1; // p1->a: 0, p2->b: 1
	  p1b = 0; p2a = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  p1a = 0; p2b = 0; // p1->a: 0, p2->b: 0
	  p1b = 0; p2a = 1; // p1->b: 0, p2->a: 1
	}
      }else if(p2_gt == '2'){ // 12_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 0; p2b = 0; // ab = 0; // p1->a: 0, p2->b: 0
	  p1b = 0; p2a = 1; // ba = 1; // p1->b: 0, p2->a: 1
	}else{ // offA == 'a', offB = 'r'
	  p1a = 0; p2b = 1; // ab = 1; // p1->a: 0, p2->b: 1
	  p1b = 0; p2a = 0; // ba = 0; // p1->b: 0, p2->a: 0
	}
      }else{ // 11_1
	continue; 
      }
    }else if(p1_gt == '2'){
      if(p2_gt == '1'){ // 21_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 1; p2b = 0; // ab = 1; // p1->a: 1, p2->b: 0
	  p1b = 0; p2a = 0; // ba = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  p1a = 0; p2b = 0; // ab = 0; // p1->a: 0, p2->b: 0
	  p1b = 1; p2a = 0; // ba = 1; // p1->b: 1, p2->a: 0
	}
      }else if(p2_gt == '0'){ // 20_1
	if(offA == 'r'){ // and offB == 'a'
	  p1a = 1; p2b = 1; // ab = 2; // p1->a: 1, p2->b: 1
	  p1b = 0; p2a = 0; // ba = 0; // p1->b: 0, p2->a: 0
	}else{ // offA == 'a', offB = 'r'
	  p1a = 0; p2b = 0; // ab = 0; // p1->a: 0, p2->b: 0
	  p1b = 1; p2a = 1; // ba = 2; // p1->b: 1, p2->a: 1
	}
      }else{ // 22_1
	continue;
      }
    }

    push_to_vlong(ft2_p1a_k, ft2_p1a_k->a[k] + p1a);
    push_to_vlong(ft2_p1b_k, ft2_p1b_k->a[k] + p1b);
    push_to_vlong(ft2_p2a_k, ft2_p2a_k->a[k] + p2a);
    push_to_vlong(ft2_p2b_k, ft2_p2b_k->a[k] + p2b);
    push_to_vlong(marker_index_of_ft2_k, i);
    k++;
  } // end of loop over markers on chromosome.

    long p1ab_min_ft2 = 10000000; // will be best (min) FT2 when switching from p1a to p1b
    long p1ab_min_k = 0; // range: 0 to number of ft2 markers; best place to switch from p1a to p1b
    long p1ba_min_ft2 = 10000000; // will be best (min) FT2 when switching from p1a to p1b
    long p1ba_min_k = 0; // range: 0 to number of ft2 markers; best place to switch from p1a to p1b
    long p2ab_min_ft2 = 10000000; // will be best (min) FT2 when switching from p1a to p1b
    long p2ab_min_k = 0; // range: 0 to number of ft2 markers; best place to switch from p1a to p1b
    long p2ba_min_ft2 = 10000000; // will be best (min) FT2 when switching from p1a to p1b
    long p2ba_min_k = 0; // range: 0 to number of ft2 markers; best place to switch from p1a to p1b

  long K = ft2_p1a_k->size-1; // number of FT2 ('type 2' forbidden triples) markers on chromosome
  long N_markers = next - first; // number of markers on chromosome
  long ft2_p1a_K = ft2_p1a_k->a[K]; // FT2 count when p1a for whole chromosome
  long ft2_p1b_K = ft2_p1b_k->a[K]; // FT2 count when p1b for whole chromosome
  long ft2_p2a_K = ft2_p2a_k->a[K]; // FT2 count when p2a for whole chromosome
  long ft2_p2b_K = ft2_p2b_k->a[K]; // FT2 count when p2b for whole chromosome
  
  fprintf(stderr, "sizes: %ld  %ld  %ld\n", k, K, N_markers);
  fprintf(stderr, "total ft2s. p1a, p1b: %ld %ld   p2a, p2b: %ld %ld   %ld %ld\n",
	  ft2_p1a_K, ft2_p1b_K, ft2_p2a_K, ft2_p2b_K,
	  ft2_p1a_K + ft2_p2b_K,  ft2_p1b_K + ft2_p2a_K);

  long p1ab_k, p1ba_k, p2ab_k, p2ba_k;
  
  for(long k=0; k<ft2_p1a_k->size; k++){
   
    p1ab_k = ft2_p1a_k->a[k] + (ft2_p1b_K - ft2_p1b_k->a[k]);
    p1ba_k = ft2_p1b_k->a[k] + (ft2_p1a_K - ft2_p1a_k->a[k]);
    p2ab_k = ft2_p2a_k->a[k] + (ft2_p2b_K - ft2_p2b_k->a[k]);
    p2ba_k = ft2_p2b_k->a[k] + (ft2_p2a_K - ft2_p2a_k->a[k]);
     fprintf(stderr, "xchrom: %ld  k:  %ld     %ld %ld     %ld %ld   %ld %ld   %ld %ld  %ld %ld\n",
	    chrom_number, k, ft2_p1a_k->a[k], ft2_p1b_k->a[k], ft2_p2a_k->a[k], ft2_p2b_k->a[k],
	     ft2_p1a_k->a[k] + ft2_p2b_k->a[k], ft2_p1b_k->a[k] + ft2_p2a_k->a[k],
	     p1ab_k, p1ba_k, p2ab_k, p2ba_k
	     );
    if(p1ab_k < p1ab_min_ft2){
      p1ab_min_ft2 = p1ab_k;
      p1ab_min_k = k;
    }
    if(p1ba_k < p1ba_min_ft2){
      p1ba_min_ft2 = p1ba_k;
      p1ba_min_k = k;
    }
    if(p2ab_k < p2ab_min_ft2){
      p2ab_min_ft2 = p2ab_k;
      p2ab_min_k = k;
    }
    if(p2ba_k < p2ba_min_ft2){
      p2ba_min_ft2 = p2ba_k;
      p2ba_min_k = k;
    }
  }
    fprintf(stderr, "ABVF: %s   chrom: %ld    %ld %ld   %ld %ld    %ld %ld   %ld %ld \n",
	    offspring->id->a, chrom_number,
	    ft2_p1a_K, ft2_p1b_K, ft2_p2a_K, ft2_p2b_K, 
	    p1ab_min_ft2, p1ba_min_ft2, p2ab_min_ft2, p2ba_min_ft2); 
  return;
} // end of FT2_phase_correction
// **********************************************************************************************************
/* */

/* XFcounts count_XF_one_chromosome(const GenotypesSet* the_gtsset, Accession* parent, Accession* other_parent,
				 Accession* offspring, long first, long next){
  // Assuming that parent is indeed a parent of offspring,
  // and calling the two homologous chromosomes in the offspring a and b
  // count the number of crossovers:
  //    for a to be derived from the parent, and
  //    for b to be derived from the parent,
  // and count the number of forbidden combinations:
  //    if a is derived from the parent, and
  //    if b is derived from the parent
  // i.e. if parent is actually the parent of offspring then
  // either a or b is derived from the parent, i.e. from the
  // some combination of the pair of homologous chromosomes in the parent
  // For each of the hypotheses 1):'a is derived from parent', and 2):'b is derived from parent'
  // we want to count both the number of crossovers required, and the number of
  // forbidden combinations
  // e.g. if the parent has dosage 0 (i.e. homozyg, ref allele) and  a = ref, b = alt (dosage +1)
  // then this counts as a forbidden combination under hypothesis 2), but is allowed under 1)

  if(parent == NULL  ||  offspring == NULL) {
    return (XFcounts){-1, -1, -1, -1, -1, -1};
  }

  //bool corrected = false;
  //Vchar* ophases = (corrected)? offspring->corrected_phases : offspring->phases;
  Vchar* ophases = offspring->phases;  //  offspring->corrected_phases;
 
  long Xa = 0, Xb = 0, Nhet = 0;
  long prev_chrom_number = -1, chrom_number;
  long prev_phase_a = -1, prev_phase_b = -1, phase_a = -1, phase_b = -1;
  long Fa = 0, Fb = 0, Nhom = 0; // counts of forbidden combinations
  long Faa = 0, Fbb = 0;

  // e.g. Fa count markers with o_gt = 0 (ref,ref) and a has alt allele, and o_gt = 2(alt,alt) and a has ref allele. 
  for(long i=first; i < next; i++){
    //long j = i-first
    chrom_number = the_gtsset->chromosomes->a[i];
    char p_gt = parent->genotypes->a[i];
    if(p_gt == MISSING_DATA_CHAR) continue;
    char p_phase = parent->phases->a[i];
  
    // do next two lines to only use markers with ok gts in all three.
    char other_parent_gt = other_parent->genotypes->a[i];
    if(other_parent_gt == MISSING_DATA_CHAR) continue;

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    char o_phase = ophases->a[i];
    long offA = (o_gt == '0'  ||  (o_gt == '1'  &&  o_phase == 'p'))? 0 : 1; // 0 <-> ref allele, 1<->alt.
    long offB = (o_gt == '0'  ||  (o_gt == '1'  &&  o_phase == 'm'))? 0 : 1; // 0 <-> ref allele, 1<->alt.
    
    // ########################################
  
    // we need to know the convention: 1+ means a is ref, and b is alt.

    if(p_gt == '1'){
      if(o_gt == '1') continue; // skip these as there errors in phases to make these introduce many spurious crossovers
      Nhet++; // counts the number of markers which are heterozyg in the parent, and homozyg in the offspring
      two_longs phases_ab = get_1marker_phases_wrt_1parent(p_phase, o_gt, o_phase);
      phase_a = phases_ab.l1;
      phase_b = phases_ab.l2;

      // compare current phases with previous values, update crossover counts, 
      // and  update prev_phase_a, prev_phase_b
      if(prev_phase_a >= 0  &&  phase_a != prev_phase_a) Xa++; // phase has changed - crossover
      prev_phase_a = phase_a;
      if(prev_phase_b >= 0  &&  phase_b != prev_phase_b) Xb++; // phase has changed - crossover
      prev_phase_b = phase_b;
    
    }else if(p_gt == '0'){
      Nhom++;
      if(o_gt == '1'  &&  (other_parent_gt == '1'  ||  other_parent_gt == '2') ){ // i.e. 01_1 or 02_1
	Faa += (offA == 1)? 1 : 0;
	Fbb += (offB == 1)? 1 : 0;
	//fprintf(stderr, "chrom,i,Faa,Fbb: %ld  %ld  %ld  %ld\n", chrom_number, i, Faa, Fbb);
      }else {
	Fa += (offA == 1)? 1 : 0;
	Fb += (offB == 1)? 1 : 0;
      }
    }else if(p_gt == '2'){
      Nhom++;
      if(o_gt == '1'  &&  (other_parent_gt == '1'  ||  other_parent_gt == '0') ){ // i.e. 21_1 or 20_1
	Faa += (offA == 1)? 0 : 1;
	Fbb += (offB == 1)? 0 : 1;
	//fprintf(stderr, "chrom,i,Faa,Fbb: %ld  %ld  %ld  %ld\n", chrom_number, i, Faa, Fbb);
      }else{
	Fa += (offA == 1)? 0 : 1;
	Fb += (offB == 1)? 0 : 1;
      }
    } // else p_gt == 'X' skip this marker
      // ######################################   
  } // end loop over markers
  return (XFcounts){Xa, Xb, Nhet, Fa, Fb, Nhom};
  } // end of count_XF_one_chromosome /* */

/* Xcounts_3 count_XF_two_parents(const GenotypesSet* the_gtsset, Accession* Fparent, Accession* Mparent, Accession* offspring){
  bool arms = true;
  long NhetF = 0, XFmin_2 = 0, XFmax_2 = 0, XFmin_3 = 0, XFmax_3 = 0;
  long NhetM = 0, XMmin_2 = 0, XMmax_2 = 0, XMmin_3 = 0, XMmax_3 = 0;
  long n_chroms = the_gtsset->chromosome_start_indices->size - 1;

  four_longs XFM_3 = {0, 0, 0, 0};
  for(long i=0; i < n_chroms; i++){
    long start_index = the_gtsset->chromosome_start_indices->a[i];
    long next_chromosome_start_index = the_gtsset->chromosome_start_indices->a[i+1];

    // the following is experimental  correction of offspring phases
    // get 'corrected' phases for offspring
    //  FT2_phase_correction(the_gtsset, Fparent, Mparent, offspring, start_index, next_chromosome_start_index);
    // three_longs FT2_counts_one_chrom =     
    //   count_FT2_one_chromosome(the_gtsset, Fparent, Mparent, offspring, start_index, next_chromosome_start_index); 

    if(ANALYZE_CHROMOSOME_ARMS){ // analyze chromosome 'arms' separately
        long next_arm_start_index = (start_index + next_chromosome_start_index)/2; // simplistic assumption: centeromere is in middle (in terms of numbers of markers)
       XFM_3 = count_XF_F_and_M_one_chromosome( the_gtsset,
					   Fparent, Mparent, offspring,
					   start_index, next_arm_start_index,
					   &XFmin_2, &XFmax_2, &NhetF,
					   &XMmin_2, &XMmax_2, &NhetM,
					   XFM_3);
       
        XFM_3 = count_XF_F_and_M_one_chromosome( the_gtsset,
					   Fparent, Mparent, offspring,
					   next_arm_start_index, next_chromosome_start_index,
					   &XFmin_2, &XFmax_2, &NhetF,
					   &XMmin_2, &XMmax_2, &NhetM,
					   XFM_3);
      XFmin_3 = XFM_3.l1;
      XFmax_3 = XFM_3.l2;
      XMmin_3 = XFM_3.l3;
      XMmax_3 = XFM_3.l4;

    }else{ // 
    
      XFM_3 = count_XF_F_and_M_one_chromosome( the_gtsset,
					   Fparent, Mparent, offspring,
					   start_index, next_chromosome_start_index,
					   &XFmin_2, &XFmax_2, &NhetF,
					   &XMmin_2, &XMmax_2, &NhetM,
					   XFM_3);
      XFmin_3 = XFM_3.l1;
      XFmax_3 = XFM_3.l2;
      XMmin_3 = XFM_3.l3;
      XMmax_3 = XFM_3.l4;
    }
    
  } // end loop over chromosomes
  Xcounts_3 X3 = {(Xcounts_2mmn){XFmin_2, XFmax_2, NhetF}, (Xcounts_2mmn){XMmin_2, XMmax_2, NhetM}, XFmin_3, XFmax_3, XMmin_3, XMmax_3};
  return X3;
} /* */

// *********************************************************************************
/* four_longs count_XF_F_and_M_one_chromosome(const GenotypesSet* the_gtsset,
					   Accession* Fparent, Accession* Mparent, Accession* offspring,
					   long start_index, long next_start_index,
					   long* XFmin_2, long* XFmax_2, long* NhetF,
					   long* XMmin_2, long* XMmax_2, long* NhetM,
					   four_longs XFM_3){
  long XFmin_3 = XFM_3.l1;
  long XFmax_3 = XFM_3.l2;
  long XMmin_3 = XFM_3.l3;
  long XMmax_3 = XFM_3.l4;
  
  XFcounts FX = count_XF_one_chromosome(the_gtsset,
					Fparent, Mparent, offspring,
					start_index, next_start_index);
  *NhetF += FX.Nhet;
   
  if(FX.Xa < FX.Xb){
    *XFmin_2 += FX.Xa; *XFmax_2 += FX.Xb;
  }else{
    *XFmin_2 += FX.Xb; *XFmax_2 += FX.Xa;
  }

  XFcounts MX = count_XF_one_chromosome(the_gtsset,
					Mparent, Fparent, offspring,
					start_index, next_start_index);
  *NhetM += MX.Nhet;
  if(MX.Xa < MX.Xb){
    *XMmin_2 += MX.Xa; *XMmax_2 += MX.Xb;
  }else{
    *XMmin_2 += MX.Xb; *XMmax_2 += MX.Xa;
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
  return (four_longs){XFmin_3, XFmax_3, XMmin_3, XMmax_3};
} /* */


// *********************************************************************************

/* XFcounts choose_P1aP2b_or_P1bP2a(XFcounts* P1ab, XFcounts* P2ab){
  double alpha = 0; // 1->use crossover rate, 0->use forbidden rate

  double P1aP2b_Xr = ((P1ab->Xa + P2ab->Xb) > 0)? (P1ab->Xa + P2ab->Xb)/(double)(P1ab->Nhet-1 + P2ab->Nhet-1) : 0;
  double P1bP2a_Xr = ((P1ab->Xb + P2ab->Xa) > 0)? (P1ab->Xb + P2ab->Xa)/(double)(P1ab->Nhet-1 + P2ab->Nhet-1) : 0;

  double P1aP2b_Fr = ((P1ab->Fa + P2ab->Fb) > 0)? (P1ab->Fa + P2ab->Fb)/(double)(P1ab->Nhom + P2ab->Nhom) : 0;
  double P1bP2a_Fr = ((P1ab->Fb + P2ab->Fa) > 0)? (P1ab->Fb + P2ab->Fa)/(double)(P1ab->Nhom + P2ab->Nhom) : 0;

  double P1aP2b_fom = alpha*P1aP2b_Xr + (1.0-alpha)*P1aP2b_Fr;
  double P1bP2a_fom = alpha*P1bP2a_Xr + (1.0-alpha)*P1bP2a_Fr;

  if(P1aP2b_fom < P1bP2a_fom){ //  P1->a,P2->b  is better
    return (XFcounts){
      P1ab->Xa + P2ab->Xb, P1ab->Xb + P2ab->Xa, P1ab->Nhet + P2ab->Nhet,
      P1ab->Fa + P2ab->Fb, P1ab->Fb + P2ab->Fa, P1ab->Nhom + P2ab->Nhom
    };
  }else{
    return (XFcounts){ //  P1->b,P2->a  is better
      P1ab->Xb + P2ab->Xa, P1ab->Xa + P2ab->Xb, P1ab->Nhet + P2ab->Nhet,
      P1ab->Fb + P2ab->Fa, P1ab->Fa + P2ab->Fb, P1ab->Nhom + P2ab->Nhom
    };
  } 
} /* */

/* void add_to_XFcounts(XFcounts* T, XFcounts I){
  T->Xa += I.Xa;
  T->Xb += I.Xb;
  T->Nhet += I.Nhet;
  T->Fa += I.Fa;
  T->Fb += I.Fb;
  T->Nhom += I.Nhom;
} /* */


/*
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
/* */


/* Vpedigree*  calculate_triples_for_one_accession_x(Accession* prog, const GenotypesSet* the_genotypes_set, Vlong* parent_idxs, long max_candidate_parents){

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
      the_ps = calculate_pedigree_stats(the_pedigree, the_genotypes_set);
      the_pedigree->pedigree_stats = the_ps;
      push_to_vpedigree(alt_pedigrees, the_pedigree);
    
    } // end loop over parent 2
  } // end loop over parent 1
  return alt_pedigrees;
} /* */

/* // this version keeps track of crossover required to get a or b separately, and chooses min for each chromosome
   // needed if want to use markers heterozygous in offspring
Xcounts_2mmn count_crossovers_one_parent_old(const GenotypesSet* the_gtsset, Accession* parent, Accession* offspring){
  // Assuming that parent is indeed a parent of offspring,
  // count the min number of crossovers needed to reconcile them

  if(parent == NULL  ||  offspring == NULL) {
    return (Xcounts_2mmn){-1, -1, -1};
  }
  long Xmin = 0, Xmax = 0, Nhet = 0; // Nhet = number of heterozyg gts in parent
  long Xa = 0, Xb = 0;
  long prev_chrom_number = -1, prev_phase_a = -1, prev_phase_b = -1;
  long phase_a = -1, phase_b = -1, chrom_number;
 
  for(long i=0; i < parent->genotypes->length; i++){

    char o_gt = offspring->genotypes->a[i]; 
    if(o_gt == MISSING_DATA_CHAR) continue;
    if(o_gt == '1') continue;
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
  fprintf(stderr, "AAAAAAABBBBBBB: %ld %ld %ld\n", Xmin, Xmax, Nhet);
  return (Xcounts_2mmn){Xmin, Xmax, Nhet};
 
} // end of count_crossovers_old 
/* */
