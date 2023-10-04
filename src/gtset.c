#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include "gtset.h"
//#include "various.h"
//#include "pedigree.h"


extern int do_checks; // option -c sets this to 1 to do some checks.
//unsigned long long bits[64];

long int_power(long base, long power){ // calculate base^power using integer math.
  long result = 1;
  for(int i=0; i<power; i++){
    result *= base;
  }
  return result;
}

// *****  Accession  implementation *****
Accession* construct_accession(char* id, long idx, char* genotypes, long accession_md_count){
  Accession* the_accession = (Accession*)malloc(1*sizeof(Accession));
  the_accession->id = construct_vchar_from_str(id);
  the_accession->index = idx;
  the_accession->genotypes = construct_vchar_from_str(genotypes);
  the_accession->chunk_patterns = NULL; // set to NULL to avoid containing garbage which causes crash when freeing accession
  the_accession->missing_data_count = accession_md_count;
  the_accession->ref_homozygs = NULL; // construct_vlong(10); // can construct in store_homozygs (so no mem used if not calling store_homozygs)
  //  the_accession->heterozygs = construct_vlong(1000);
  the_accession->alt_homozygs = NULL; // construct_vlong(10);
  the_accession->Abits = NULL;
  the_accession->Bbits = NULL;
  return the_accession;
}
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count){
  the_accession->missing_data_count = missing_data_count;
}

// for one accession's set of genotypes, loop over chunks and find the gt patterns. Store in the_gts->chunk_patterns
long set_accession_chunk_patterns(Accession* the_gts, Vlong* m_indices, long n_chunks, long k, long ploidy){
  long gts_mdchunk_count = 0;
  long n_patterns = int_power(ploidy+1, k); // 3^k, the number of valid patterns, also there is a 'pattern' for missing data, making 3^k + 1 in all
  //  the_gts->chunk_homozyg_counts = construct_vlong_zeroes(n_chunks);
  Vlong* chunk_pats = construct_vlong(n_chunks); // (Vlong*)malloc(n_chunks*sizeof(Vlong));
  //  fprintf(stderr, "m_indices->size %ld  n_chunks: %ld \n", m_indices->size, n_chunks);
  long n_h_usable_chunks = 0;
  long n_usable_homozygs = 0; // the number of homozygous gts in usable chunks.
  for(long i_chunk=0; i_chunk < n_chunks; i_chunk++){
    long n_homozygs = 0;
    long i_chunkstart = k*i_chunk;
    long i_pat = 0;
    long f = 1;

    // loop over characters in the chunk and construct a corresponding long index, in range [0..3^k] (3^k is the index for a chunk with any missing data)
    for(long j=0; j < k; j++){ 
      long idx = m_indices->a[i_chunkstart + j]; // 
      char a = the_gts->genotypes->a[idx];
      long l = (long)a - 48;
      if((l>=0) && (l<=ploidy)){ // this char is ok (0, 1, or 2, not missing data)
	i_pat += f*l;
	f *=3;
	if(l == 0  ||  l == 2){
	  n_homozygs++;
	}
      }else{ // missing data in (at least) one of the markers in the chunk
	i_pat = n_patterns;
	gts_mdchunk_count++;
	n_homozygs = -1;
	break;
      }
    } // end of loop over the k chars in a chunk.
    push_to_vlong(chunk_pats, i_pat);
    
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  the_gts->md_chunk_count = gts_mdchunk_count;
  return gts_mdchunk_count;
}


char* print_accession(Accession* the_gts, FILE* ostream){
  // fprintf(ostream, "Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->genotypes);
  // fprintf(ostream, "%s  %s\n", the_gts->id, the_gts->genotypes);
  fprintf(ostream, "Accession id: %s index: %ld length: %ld missing data count:  %ld\n",
	  the_gts->id->a, the_gts->index, the_gts->genotypes->length, the_gts->missing_data_count);
}

void free_accession(Accession* the_accession){
  if(the_accession == NULL) return;
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
  free_vlong(the_accession->chunk_patterns);
  free_vlong(the_accession->alt_homozygs);
  free_vlong(the_accession->ref_homozygs);
  free_vull(the_accession->Abits);
  free_vull(the_accession->Bbits);
  free(the_accession);
}

double agmr0(GenotypesSet* the_gtsset){  
  /* Vlong* dosage_counts = the_gtsset->dosage_counts; */
  /* double Meff = (double)the_gtsset->n_markers*the_gtsset->n_accessions - dosage_counts->a[3]; */
  /* double f0 = (double)dosage_counts->a[0]/Meff; */
  /* double f1 = (double)dosage_counts->a[1]/Meff; */
  /* double f2 = (double)dosage_counts->a[2]/Meff; */
  
  /* fprintf(stderr, "# Meff: %8.5f  %ld %ld %ld \n", Meff, dosage_counts->a[0], dosage_counts->a[1], dosage_counts->a[2]); */
  /* fprintf(stderr, "# %8.5f  %8.5f  %8.5f\n", f0, f1, f2); */
  /* return 2*(f1*(f0+f2) + f0*f2); */
  double n_exp_different = 0;
  double n_exp_ok = 0; 
  
  for(long i=0; i<the_gtsset->marker_ids->size; i++){
    double ok_count = (double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
    double f0 = the_gtsset->marker_dosage_counts[0]->a[i]/ok_count;
    double f1 = the_gtsset->marker_dosage_counts[1]->a[i]/ok_count;
    double f2 = the_gtsset->marker_dosage_counts[2]->a[i]/ok_count;;
    n_exp_different += 2*(f0*(f1+f2) + f1*f2);
    n_exp_ok += ok_count/(double)the_gtsset->accessions->size;
  
    // fprintf(stderr, "## f0, f1, f2: %8.5f %8.5f %8.5f  ok_count: %8.5f  n_exp_diff: %8.5f nok %8.5f\n",
    //	  f0, f1, f2, ok_count, n_exp_different, n_exp_ok);
}
  return (n_exp_ok > 0)? n_exp_different/n_exp_ok : -1;
}
double agmr0_accvsall(const GenotypesSet* the_gtsset, Accession* A){
   double n_exp_different = 0;
   double n_exp_ok = 0;
   
   for(long i=0; i<the_gtsset->marker_ids->size; i++){
    double ok_count = (double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
    double f0 = the_gtsset->marker_dosage_counts[0]->a[i]/ok_count;
    double f1 = the_gtsset->marker_dosage_counts[1]->a[i]/ok_count;
    double f2 = the_gtsset->marker_dosage_counts[2]->a[i]/ok_count;;

    if(A->genotypes->a[i] - 48  == 0){
      n_exp_different += (f1+f2);
      n_exp_ok += ok_count/(double)the_gtsset->accessions->size;
    }else if(A->genotypes->a[i] - 48 == 1){
      n_exp_different += (f0+f2);
      n_exp_ok += ok_count/(double)the_gtsset->accessions->size;
    }else if(A->genotypes->a[i] - 48 == 2){
      n_exp_different += (f0+f1);
      n_exp_ok += ok_count/(double)the_gtsset->accessions->size;
    }else{ // missing data
      
    }  
  /* fprintf(stderr, "## f0, f1, f2: %8.5f %8.5f %8.5f  ok_count: %8.5f  n_exp_diff: %8.5f nok %8.5f\n", */
  /* 	  f0, f1, f2, ok_count, n_exp_different, n_exp_ok); */
   }
   return (n_exp_ok > 0)? n_exp_different/n_exp_ok : -1;
}


double pair_agmr0(Accession* A, Accession* B){ // find expected agmr for two accs with totally uncorrelated gts.
  long M = A->genotypes->length;
  long A0 = A->ref_homozygs->size;
  long A2 = A->alt_homozygs->size;
  long AX = A->missing_data_count;
  long A1 = M - (A0 + A2 + AX); // heterozygs
  
  long B0 = B->ref_homozygs->size;
  long B2 = B->alt_homozygs->size;
  long BX = B->missing_data_count;
  long B1 = M - (B0 + B2 + BX); // heterozygs
    
  if( (M - AX) > 0  &&  (M - BX) > 0){
    return (double)(A0*(B1+B2) + A1*(B0+B2) + A2*(B0+B1))/(double)((M-AX)*(M-BX));
  }else{
    return -1;
  }
}

double agmr0_qvsall(const GenotypesSet* the_gtsset, Accession* A){
  Vlong* dosage_counts = the_gtsset->dosage_counts;
  double Meff = (double)the_gtsset->n_markers*the_gtsset->n_accessions - dosage_counts->a[3];
  double f0 = (double)dosage_counts->a[0]/Meff;
  double f1 = (double)dosage_counts->a[1]/Meff;
  double f2 = (double)dosage_counts->a[2]/Meff;

  long M = A->genotypes->length;
  long A0 = A->ref_homozygs->size;
  long A2 = A->alt_homozygs->size;
  long AX = A->missing_data_count;
  long A1 = M - (A0 + A2 + AX); // heterozygs

  double agmr0 = (double)(f0*(A1+A2) + f1*(A0+A2) + f2*(A0 + A1)) / (double)(M - AX);
  return agmr0;
}

/* void free_accession_innards(Accession* the_accession){ // doesn't free the_accession itself */
/*   // use to free each element of an array of Accession (not an array of Accession*) */
/*   if(the_accession == NULL) return; */
/*   free_vchar(the_accession->id); */
/*   free_vchar(the_accession->genotypes); */
/* } */

// *****  Vaccession implementation  *****
Vaccession* construct_vaccession(long cap){ // construct Vaccession with capacity of cap, but empty.
  Vaccession* the_vacc = (Vaccession*)malloc(sizeof(Vaccession));
  the_vacc->capacity = cap;
  the_vacc->size = 0;
  the_vacc->a = (Accession**)malloc(cap*sizeof(Accession*));
  return the_vacc;
}

void push_to_vaccession(Vaccession* the_vacc, Accession* the_acc){
  long cap = the_vacc->capacity;
  long n = the_vacc->size;
   if(n == cap){   // if necessary, resize w realloc
    cap *= 2;
    the_vacc->a = (Accession**)realloc(the_vacc->a, cap*sizeof(Accession*));
    the_vacc->capacity = cap;
  }
  the_vacc->a[n] = the_acc;
   the_vacc->size++;
}


void shuffle_order_of_accessions(GenotypesSet* the_genotypes_set){
  long n_ref = the_genotypes_set->n_ref_accessions;
  long n_new = the_genotypes_set->accessions->size - n_ref;
  Vaccession* accessions = the_genotypes_set->accessions; // construct_vaccession(n_ref + n_new);

  Accession* temp_acc;
   for(long i=0; i<n_ref-1; i++){
    long j = i+1 + (long)((n_ref-1-i)*(double)rand()/RAND_MAX); // get an integer in the range [i+1,n_ref-1]
    temp_acc = accessions->a[j];
    accessions->a[j] = accessions->a[i];
    //   accessions->a[j]->index = j;
    accessions->a[i] = temp_acc;
    //   accessions->a[i]->index = i;
  }

  for(long i=0; i<n_new-1; i++){
    long j = i+1 + (long)((n_new-1-i)*(double)rand()/RAND_MAX); // get an integer in the range [i+1,n_new-1]
    temp_acc = accessions->a[j + n_ref];
    accessions->a[j + n_ref] = accessions->a[i + n_ref];
    //  accessions->a[j + n_ref]->index = j + n_ref;
    accessions->a[i + n_ref] = temp_acc;
    //  accessions->a[i + n_ref]->index = i + n_ref;
  }
  for(long i=0; i<the_genotypes_set->accessions->size; i++){
    accessions->a[i]->index = i;
  }
 
}

void set_vaccession_chunk_patterns(Vaccession* the_accessions, Vlong* m_indices, long n_chunks, long k, long ploidy){
  long total_mdchunk_count = 0;
  for(long i=0; i < the_accessions->size; i++){
    long mdchcount = set_accession_chunk_patterns(the_accessions->a[i], m_indices, n_chunks, k, ploidy);
    total_mdchunk_count += mdchcount;
  }
}

void print_vaccession(Vaccession* the_accessions, FILE* ostream){
  for(int i=0; i<the_accessions->size; i++){
    print_accession(the_accessions->a[i], ostream);
  }

}

void check_accession_indices(Vaccession* the_accessions){
  for(long i=0; i<the_accessions->size; i++){
    Accession* a_gts = the_accessions->a[i];
    if(a_gts->index != i){
      fprintf(stderr, "In check_gts_indices. i: %ld  index: %ld\n", i, a_gts->index);
      exit(EXIT_FAILURE);
    }
  }
}

void free_vaccession(Vaccession* the_vacc){
  if(the_vacc == NULL) return;
  //  fprintf(stderr, "# in free_vaccessions. XXXXXXXXXXXXXXXXXX\n");
  for(long i=0; i<the_vacc->size; i++){
    free_accession(the_vacc->a[i]);
  }
  free(the_vacc->a);
  free(the_vacc);
}

// *****  GenotypesSet implementation *****
GenotypesSet* construct_empty_genotypesset(double max_marker_md_fraction, double min_min_allele_freq, long ploidy){
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->max_marker_missing_data_fraction = max_marker_md_fraction;
  the_gtsset->min_minor_allele_frequency = min_min_allele_freq;
  the_gtsset->n_accessions = -1;
  the_gtsset->n_bad_accessions = 0;
  the_gtsset->n_ref_accessions = 0;
  the_gtsset->n_markers = -1;
  the_gtsset->ploidy = ploidy;
  the_gtsset->accessions = construct_vaccession(INIT_VACC_CAPACITY);
  the_gtsset->marker_ids = NULL;
  the_gtsset->marker_missing_data_counts = NULL;
  the_gtsset->marker_alt_allele_counts = NULL;
  the_gtsset->mafs = NULL;
  the_gtsset->marker_dosage_counts = NULL;
  return the_gtsset;
}

void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set, double max_acc_missing_data_fraction){
  FILE* g_stream = fopen(input_filename, "r");
  if (g_stream == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  // *****   Read first line; store marker ids.  *****
  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  char* saveptr = NULL;
  while((nread = getline(&line, &len, g_stream)) != -1){
    saveptr = line;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    
    if((token == NULL) || (token[0] == '#')) continue; // skip comments, empty lines
    if((strcmp(token, "MARKER") == 0)){
      while(1){
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	if(token == NULL) break;
	char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
	push_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
	markerid_count++;
      }
      break;
    }else{ 
      fprintf(stderr, "token: %s (should be MARKER)\n", token);
      exit(EXIT_FAILURE);
    }
  }
  // *****  done reading first line (with marker ids)  *****
  if(the_genotypes_set->marker_missing_data_counts == NULL){    
    the_genotypes_set->marker_missing_data_counts = construct_vlong_zeroes(marker_ids->size);
    the_genotypes_set->marker_alt_allele_counts = construct_vlong_zeroes(marker_ids->size);
    
    the_genotypes_set->marker_dosage_counts = (Vlong**)malloc((MAX_PLOIDY+1)*sizeof(Vlong*));
    for(long i=0; i<=MAX_PLOIDY; i++){
      the_genotypes_set->marker_dosage_counts[i] = construct_vlong_zeroes(marker_ids->size);
    }
  }
  
  // if the_genotypes_set already has marker_ids set (from a reference set), check for agreement with marker_ids
  // i.e. in case of using reference set, make sure the set of markers is the same in both data sets.
  if(the_genotypes_set->marker_ids  == NULL){
    the_genotypes_set->marker_ids = marker_ids;
  }else{
    if(the_genotypes_set->marker_ids->size != marker_ids->size){
      fprintf(stderr, "# Data sets have different numbers of markers: %5ld  %5ld. Exiting.\n", the_genotypes_set->marker_ids->size, marker_ids->size);
      exit(EXIT_FAILURE);
    }else{
      long idx = compare_vstrs(marker_ids, the_genotypes_set->marker_ids);
      assert(idx >= 0);
      if(idx < marker_ids->size){ 
	fprintf(stderr, "# Warning: Set of marker ids does not agree with reference.\n");
	fprintf(stderr, "# First pair of non-matching ids: %s %s %ld\n", the_genotypes_set->marker_ids->a[idx], marker_ids->a[idx], idx);
      }else{
	fprintf(stdout, "# Marker ids agree with reference set.\n");
      }
    }
    free_vstr(marker_ids);
  }

  // Read in the rest of the lines, construct Accession for each line
  long accession_count = 0;
  while((nread = getline(&line, &len, g_stream)) != -1){
    // fprintf(stderr, "# reading accession %ld\n", accession_count);
    saveptr = line;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    long marker_count = 0;
    long accession_missing_data_count = 0;
    char* genotypes = (char*)calloc((markerid_count+1), sizeof(char));    
    genotypes[markerid_count] = '\0'; // terminate with null.

    while(1){ // read dosages from one line.   
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL)	break;    
      genotypes[marker_count] = token_to_dosage(token, &(the_genotypes_set->ploidy));
      if(genotypes[marker_count] == MISSING_DATA_CHAR) accession_missing_data_count++;
      marker_count++;
    } // done reading dosages for all markers of this accession
    if(marker_count != markerid_count) {
      fprintf(stderr, "# marker_count, markerid_count: %ld %ld \n", marker_count, markerid_count);
      exit(EXIT_FAILURE); 
    }
    double fraction_to_keep = 0.75;
    // if accession does not have too much missing data, construct Accession and store in the_genotypes_set
    // if(0 || (double)rand()/((double)(RAND_MAX)+1) < fraction_to_keep){
     if(1 || (double)rand()/((double)(RAND_MAX)+1) < fraction_to_keep){   
    if(accession_missing_data_count <= max_acc_missing_data_fraction * the_genotypes_set->marker_ids->size){
      Accession* the_accession = construct_accession(acc_id, accession_count, genotypes, accession_missing_data_count);
      for(long jjj=0; jjj<markerid_count; jjj++){ //
	if(genotypes[jjj] == MISSING_DATA_CHAR){
	  the_genotypes_set->marker_missing_data_counts->a[jjj]++;
	  accession_missing_data_count++;
	}else{
	  the_genotypes_set->marker_alt_allele_counts->a[jjj] += (long)(genotypes[jjj]-48); // 48->+=0, 49->+=1, 50->+=2, etc.
	}
      }   
      push_to_vaccession(the_genotypes_set->accessions, the_accession);
      accession_count++;
    }else{
      fprintf(stdout, "# Accession: %s rejected due to missing data at %ld out of %ld markers.\n",
	      acc_id, accession_missing_data_count, the_genotypes_set->marker_ids->size);
      the_genotypes_set->n_bad_accessions++;
    }
    }
     free(acc_id); // or cut out the middleman (acc_id)?
      free(genotypes);
  } // done reading all lines
  fclose(g_stream);
  // fprintf(stderr, "# %ld accessions removed for excessive missing data; %ld accessions kept.\n",
  //	  the_genotypes_set->n_bad_accessions, accession_count);
  free(line); // only needs to be freed once.
  the_genotypes_set->n_accessions = the_genotypes_set->accessions->size;
  the_genotypes_set->n_markers = the_genotypes_set->marker_ids->size;
  if(DBUG && do_checks) check_genotypesset(the_genotypes_set);
}

long str_to_long(char* str){ // using strtol and checking for various problems.
  char *end;
  errno = STRTOL_FAIL; // I guess this is global declared in errno.h ?
  const long sl = strtol(str, &end, 10);
  if (end == str) {
    //  fprintf(stderr, "%s: not a decimal number\n", str);
  } else if ('\0' != *end) {
    //  fprintf(stderr, "%s: extra characters at end of input: %s\n", str, end);
  } else if ((LONG_MIN == sl || LONG_MAX == sl) && ERANGE == errno) {
    //  fprintf(stderr, "%s out of range of type long\n", str);
  } else if (sl > INT_MAX) {
    //  fprintf(stderr, "%ld greater than INT_MAX\n", sl);
  } else if (sl < INT_MIN) {
    //  fprintf(stderr, "%ld less than INT_MIN\n", sl);
  }else{ // success
    errno = 0;
  }
  return sl;
}

char token_to_dosage(char* token, long* ploidy){
  if(strcmp(token, "X") == 0){ // missing data
    return MISSING_DATA_CHAR;
  }else{
    long l = str_to_long(token);
    if(errno != 0){
      exit(EXIT_FAILURE);
    }else if(l < 0  ||  l > MAX_PLOIDY){ // if outside the range of possible dosages.
      return MISSING_DATA_CHAR;
    }else{
      if(l > *ploidy) *ploidy = l; // so ploidy ends up as the max of all dosages.
      return (char)l + 48;
    }
  }
}

void populate_marker_dosage_counts(GenotypesSet* the_gtsset){
  long nX=0;
  the_gtsset->dosage_counts = construct_vlong_zeroes(4); // (long*)malloc(4*sizeof(long));
  for(long j=0; j<the_gtsset->n_accessions; j++){
    Accession* the_acc = the_gtsset->accessions->a[j];
    for(long i=0; i<the_gtsset->n_markers; i++){
      long dosage = (long)the_acc->genotypes->a[i];
      if(dosage != MISSING_DATA_CHAR){
	the_gtsset->marker_dosage_counts[dosage-48]->a[i]++;
      }else{
	nX++;
      }
    }
  }
   for(long i=0; i<the_gtsset->n_markers; i++){
     for(long j=0; j<=2; j++){
       the_gtsset->dosage_counts->a[j] += the_gtsset->marker_dosage_counts[j]->a[i];
     }
   }
   the_gtsset->dosage_counts->a[3] += nX;
}

void set_agmr0s(GenotypesSet* the_gtsset){
  for(long i=0; i< the_gtsset->accessions->size; i++){
    Accession* A = the_gtsset->accessions->a[i];
    A->agmr0 = agmr0_accvsall(the_gtsset, A);
  }
}

double ragmr(GenotypesSet* the_gtsset){
  for(long j=0; j<the_gtsset->n_accessions; j++){
    Accession* the_acc = the_gtsset->accessions->a[j];
    for(long i=0; i<the_gtsset->n_markers; i++){
      long dosage = (long)the_acc->genotypes->a[i];
      if(dosage != MISSING_DATA_CHAR){
	the_gtsset->marker_dosage_counts[dosage-48]->a[i]++;
      }
    }
  }

  double ragmr = 0;
  for(long i=0; i<the_gtsset->n_markers; i++){
    double sum_nsq = 0;
    double sum_n = 0;
    for(long k=0; k<=the_gtsset->ploidy; k++){
      double n = (double)the_gtsset->marker_dosage_counts[k]->a[i];
      sum_n += n;
      sum_nsq += n*n;
    }
    double x = (1.0 - sum_nsq/(sum_n*sum_n));
    // fprintf(stderr, "# sum_n: %8.0f  sum_nsq  %12.0f  x: %
    ragmr += (1.0 - sum_nsq/(sum_n*sum_n)); 
  }
  ragmr /= the_gtsset->n_markers;
  //   fprintf(stderr, "#  ragmr: 
  // fprintf(stderr, "bottom of ragmr\n");

  return ragmr;
}

void check_genotypesset(GenotypesSet* gtss){
  // recalculates the maf and missing data of each marker,
  // checks that these are same as values stored in gtss
  // do this after filtering to check that filtering
  // correctly adjusted the maf and missing data for each marker.
  assert(gtss->marker_ids->size == gtss->n_markers);
  assert(gtss->accessions->size == gtss->n_accessions);
  long* marker_md_counts = (long*)calloc(gtss->n_markers, sizeof(long));
  long* marker_alt_allele_counts = (long*)calloc(gtss->n_markers, sizeof(long));
  for(long i=0; i<gtss->n_accessions; i++){
    Accession* an_acc = gtss->accessions->a[i];
    assert(strlen(an_acc->genotypes->a) == gtss->n_markers);
    long accmdcount = 0;
    for(long j=0; j<gtss->n_markers; j++){
      if(an_acc->genotypes->a[j] == MISSING_DATA_CHAR){
	marker_md_counts[j]++;
	accmdcount++;
      }else{
	marker_alt_allele_counts[j] += (long)(an_acc->genotypes->a[j] - 48);
      }
    }
    assert(an_acc->missing_data_count == accmdcount);
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(marker_md_counts[j] == gtss->marker_missing_data_counts->a[j]);
    assert(marker_alt_allele_counts[j] == gtss->marker_alt_allele_counts->a[j]);
  }
  free(marker_md_counts);
  fprintf(stderr, "# Successfully completed check_genotypesset\n");
}

void filter_genotypesset(GenotypesSet* the_gtsset){ // construct a new set of 'filtered' accession genotypes, which replace the raw ones.
  double max_marker_md_fraction = the_gtsset->max_marker_missing_data_fraction;
  Vlong* marker_md_counts = the_gtsset->marker_missing_data_counts; // the number of missing data for each marker
  Vlong* alt_allele_counts = the_gtsset->marker_alt_allele_counts;  
  long n_accs = the_gtsset->n_accessions;
  
  if(max_marker_md_fraction < 0){ // set max_marker_md_fraction to some multiple of median md fraction.
    // (but max_marker_md_fraction is set by default in duplicatesearch.c to 2/chunk_size;
    // to use this way of defining max_marker_md_fraction, specify, e.g. -x -5 to set max_marker_md_fraction to 5 time median.)
    long* mdcount_histogram = (long*)calloc(n_accs+1, sizeof(long));
    double factor = -1*max_marker_md_fraction;
    for(long i=0; i< marker_md_counts->size; i++){
      mdcount_histogram[marker_md_counts->a[i]]++;
    }
    long mrkrs_so_far = 0;
    long median_md_count;
    for(long i=0; i<=n_accs; i++){
      mrkrs_so_far += mdcount_histogram[i];
      if(mrkrs_so_far > 0.5*marker_md_counts->size){
	median_md_count = i;
	break;
      }
    }
    max_marker_md_fraction = factor*(double)median_md_count/(double)n_accs;
  } // end of optionally setting max_marker_md_fraction to multiple of median
  
  // identify the markers to keep:
  long n_markers_to_keep = 0;
  Vlong* md_ok = construct_vlong_zeroes(marker_md_counts->size); // set to 1 if number of missing data gts is small enough
  Vstr* filtered_marker_ids = construct_vstr(1000);
  Vlong* filtered_marker_md_counts = construct_vlong(1000);
  Vlong* filtered_alt_allele_counts = construct_vlong(1000);
  long mdsum_all = 0;
  long mdsum_kept = 0;
  long altallelesum_all = 0;
  long altallelesum_kept = 0;
  long only_one_allele_count = 0;
  long too_much_missing_data_count = 0;
  long maf_too_low_count = 0;
  for(long i=0; i<marker_md_counts->size; i++){
    mdsum_all += marker_md_counts->a[i];
    altallelesum_all += alt_allele_counts->a[i];
  
    if (marker_md_counts->a[i] <= max_marker_md_fraction*the_gtsset->n_accessions){  // not too much missing data
      long max_possible_alt_alleles = the_gtsset->ploidy * (the_gtsset->n_accessions - marker_md_counts->a[i]);
      double min_min_allele_count = fmax(max_possible_alt_alleles * the_gtsset->min_minor_allele_frequency, 1);
      double max_min_allele_count = fmin(max_possible_alt_alleles * (1.0 - the_gtsset->min_minor_allele_frequency), max_possible_alt_alleles-1);
      bool alt_allele_freq_ok = (alt_allele_counts->a[i] >= min_min_allele_count)  && // alternative allele frequency not too small,
	(alt_allele_counts->a[i] <= max_min_allele_count);     // and not too large
      if ( alt_allele_freq_ok ){ // the alt_allele_frequency is in the right range
	md_ok->a[i] = 1;
	n_markers_to_keep++;
	mdsum_kept += marker_md_counts->a[i];
	altallelesum_kept += alt_allele_counts->a[i];
	long marker_id_length = strlen(the_gtsset->marker_ids->a[i]);
	char* marker_id_to_keep = strcpy((char*)malloc((marker_id_length+1)*sizeof(char)), the_gtsset->marker_ids->a[i]);
	push_to_vstr(filtered_marker_ids, marker_id_to_keep); // store the kept marker ids.
	push_to_vlong(filtered_marker_md_counts, marker_md_counts->a[i]); // store the md counts for the kept markers.
	push_to_vlong(filtered_alt_allele_counts, alt_allele_counts->a[i]);
      }else{
	maf_too_low_count++;
      }
    }else{
      too_much_missing_data_count++;
    }
  } // end of loop over markers
  // fprintf(stderr, "# mdsum_all: %ld  mdsum_kept: %ld  n_markers all: %ld  n_markers_to_keep: %ld\n", mdsum_all, mdsum_kept, marker_md_counts->size, n_markers_to_keep);
  double raw_md_fraction = (double)mdsum_all/(double)(marker_md_counts->size*n_accs);
  double filtered_md_fraction = (double)mdsum_kept/(double)(n_markers_to_keep*n_accs);
  double raw_minor_allele_freq = (double)altallelesum_all/(marker_md_counts->size*n_accs*the_gtsset->ploidy);
  double filtered_minor_allele_freq = (double)altallelesum_kept/(double)(n_markers_to_keep*n_accs*the_gtsset->ploidy);

  /* fprintf(stderr, "# Removing markers with missing data fraction > %5.3lf or minor allele frequency < %5.3f\n", */
  /* 	  max_marker_md_fraction, the_gtsset->min_minor_allele_frequency); */
  /* fprintf(stderr, "#   removed %ld markers for excessive missing data.\n", too_much_missing_data_count); */
  /* fprintf(stderr, "#   removed an additional %ld markers for too small maf.\n", maf_too_low_count); */
   fprintf(stdout, "# Removing markers with missing data fraction > %5.3lf or minor allele frequency < %5.3f\n",
  	  max_marker_md_fraction, the_gtsset->min_minor_allele_frequency);
  fprintf(stdout, "#   removed %ld markers for excessive missing data.\n", too_much_missing_data_count);
  fprintf(stdout, "#   removed an additional %ld markers for too small maf.\n", maf_too_low_count);
  
  fprintf(stdout, "# Raw data has %ld markers, missing data fraction = %5.3lf, minor allele frequencey = %5.3f\n",
	  marker_md_counts->size, raw_md_fraction, raw_minor_allele_freq);
  fprintf(stdout, "# Filtered data has %ld markers, missing data fraction = %5.3lf, minor allele frequency = %5.3lf\n",
	  n_markers_to_keep, filtered_md_fraction, filtered_minor_allele_freq);

  Vaccession* the_accessions = construct_vaccession(the_gtsset->n_accessions); //(Accession*)malloc(the_gtsset->n_accessions*sizeof(Accession)); 
  for(long i=0; i<the_gtsset->n_accessions; i++){ // loop over accessions
    char* raw_gts = the_gtsset->accessions->a[i]->genotypes->a; // the string with all the genotypes for accession i
    char* filtered_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    long acc_md_count = 0;
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	filtered_gts[k] = raw_gts[j];				       
	k++;
	if(raw_gts[j] == MISSING_DATA_CHAR) acc_md_count++;
      }
    }
    filtered_gts[k] = '\0'; // terminate with null.
    if(DBUG && do_checks) assert(k == n_markers_to_keep);
    Accession* the_accession = construct_accession( the_gtsset->accessions->a[i]->id->a, i, filtered_gts, acc_md_count ); //(Accession*)malloc(sizeof(Accession));
    free(filtered_gts); // this str was copied into newly allocated memory in the_accession->
    push_to_vaccession(the_accessions, the_accession);
  }
  free_vlong(md_ok);

  free_vaccession(the_gtsset->accessions);
  the_gtsset->accessions = the_accessions;
  
  free_vstr(the_gtsset->marker_ids);
  the_gtsset->marker_ids = filtered_marker_ids;

  free_vlong(the_gtsset->marker_missing_data_counts);
  the_gtsset->marker_missing_data_counts = filtered_marker_md_counts;

  free_vlong(the_gtsset->marker_alt_allele_counts);
  the_gtsset->marker_alt_allele_counts = filtered_alt_allele_counts;

  the_gtsset->n_markers = the_gtsset->accessions->a[0]->genotypes->length;
} // end filter_genotypesset

void rectify_markers(GenotypesSet* the_gtsset){ // if alt allele has frequency > 0.5, dosage -> ploidy-dosage
  long n_accessions = the_gtsset->accessions->size;
  long n_markers = the_gtsset->marker_alt_allele_counts->size;
  long ploidy = the_gtsset->ploidy;
  long delta_dosage;
  long n_markers_rectified = 0;
  for(long i=0; i<n_markers; i++){
    delta_dosage = 0;
    //  fprintf(stderr, "altallcounts, okmrkrs: %ld  %ld\n", the_gtsset->marker_alt_allele_counts->a[i], (n_accessions - the_gtsset->marker_missing_data_counts->a[i]) );
    if(2*the_gtsset->marker_alt_allele_counts->a[i] > ploidy*(n_accessions - the_gtsset->marker_missing_data_counts->a[i])){ //
      n_markers_rectified++;
      for(long j=0; j<n_accessions; j++){
	char dosage = the_gtsset->accessions->a[j]->genotypes->a[i];
	if(dosage != MISSING_DATA_CHAR){
	  long ldosage = dosage-48;
	  if(2*ldosage != ploidy){
	    the_gtsset->accessions->a[j]->genotypes->a[i] = (char)(ploidy - ldosage) + 48;
	    delta_dosage += ploidy - 2*ldosage;
	  }
	}
      }
    }
    // so far genotypes of each accession have (if necessary) been 'rectified', but now marker_alt_allele_counts needs to be updated
    //  fprintf(stderr, "%ld  %ld   ", the_gtsset->marker_alt_allele_counts->a[i], delta_dosage);
    the_gtsset->marker_alt_allele_counts->a[i] += delta_dosage;
  }
  fprintf(stderr, "##### n_markers_rectified: %ld\n", n_markers_rectified);
}

void set_Abits_Bbits(GenotypesSet* the_genotypesset){ // diploid only
  unsigned long long x = 1, bits[64];
  for(long j=0; j<64; j++, x = x<<1){
    bits[j] = x; // bits[i] is an unsigned long long with the ith bit set (i.e. 1, with all others 0)
  }

  for(long i_acc =0; i_acc<the_genotypesset->accessions->size; i_acc++){
    Accession* the_acc = the_genotypesset->accessions->a[i_acc];
    long n_markers = the_acc->genotypes->length;
    long n_ulls = n_markers/64; // Abits, Bbits will each have this many longs
    n_ulls++; // add another ull which will have only some bits used (because n markers may not be multiple of chunk size)
    Vull* Abits = construct_vull_zeroes(n_ulls);
    Vull* Bbits = construct_vull_zeroes(n_ulls);
    Abits->size = n_ulls;
    Bbits->size = n_ulls;
    for(long i_ull=0; i_ull<n_ulls; i_ull++){
      unsigned long long A=0, B=0;
      for(long j=0; j<64; j++){ // loop over the bits 
	long i_gt = 64*i_ull + j;
	char gt = (i_gt < n_markers)? the_acc->genotypes->a[i_gt] : 'X';
	if(gt == '0'){ // 00
	  // leave the bit as 0 in both A and B.
	}else if(gt == '1'){ // 01
	  B |= bits[j];
	 
	}else if(gt == '2'){ // 11
	  A |= bits[j];
	  B |= bits[j];
	}else{ // missing data; 10
	  A |= bits[j];
	}
	// fprintf(stderr, "# acc: %s  %ld  %c  %llu  %llu\n", the_acc->id->a, i_gt, gt, A, B);
      }
      Abits->a[i_ull] = A;
      Bbits->a[i_ull] = B;
      /* fprintf(stderr, "# acc id: %s\n", the_acc->id->a); */
      /* fprintf(stderr, "# A: %llu\n", A); */
      /* fprintf(stderr, "# B: %llu\n", B); */
    }
    the_acc->Abits = Abits;
    the_acc->Bbits = Bbits;
  }
}


void store_homozygs(GenotypesSet* the_gtsset){ // for each accession,
  // store indices of alt homozygous genotypes, and (separately) ref homozyg gts.
  for(long i=0; i<the_gtsset->accessions->size; i++){
    Accession* acc = the_gtsset->accessions->a[i];
     acc->alt_homozygs = construct_vlong(10000);
     acc->ref_homozygs = construct_vlong(10000);
    for(long j=0; j<the_gtsset->marker_alt_allele_counts->size; j++){
      if(acc->genotypes->a[j] == MISSING_DATA_CHAR) continue;
      long dosage = (long)(acc->genotypes->a[j] - 48); // 0, 1, ..., ploidy
      if( dosage == the_gtsset->ploidy){ // if dosage == ploidy
	push_to_vlong(acc->alt_homozygs, j);	
      }else if(dosage == 0){
	push_to_vlong(acc->ref_homozygs, j);
      }
    }
  }
}

Vdouble* get_minor_allele_frequencies(GenotypesSet* the_gtset){
  Vlong* missing_data_counts = the_gtset->marker_missing_data_counts;
  Vlong* minor_allele_counts = the_gtset->marker_alt_allele_counts;
  if(DO_ASSERT) assert(missing_data_counts->size == minor_allele_counts->size);
  Vdouble* marker_mafs = construct_vdouble(missing_data_counts->size);
  for(long i=0; i<missing_data_counts->size; i++){
    long ok_count = the_gtset->n_accessions - missing_data_counts->a[i];
    if(ok_count > 0){
      double minor_allele_frequency = (double)minor_allele_counts->a[i]/(double)(2.0*ok_count);
      push_to_vdouble(marker_mafs, minor_allele_frequency);
    }
  }
  the_gtset->mafs = marker_mafs;
  return marker_mafs;
}

four_longs bitwise_agmr_hgmr(Accession* acc1, Accession* acc2){
  /* fprintf(stderr, "Abits->size: %ld %ld\n", acc1->Abits->size, acc1->Bbits->size); */
  /*  */
  //long Ndo = 0, Nso = 0, Ndx = 0, Nsx = 0;
  four_longs rval = {0, 0, 0, 0};
  unsigned long long isOi, isXi, isOj, isXj;
  for(long i_long = 0; i_long < acc1->Abits->size; i_long++){
    unsigned long long iA = acc1->Abits->a[i_long];
    unsigned long long iB = acc1->Bbits->a[i_long];
    unsigned long long jA = acc2->Abits->a[i_long];
    unsigned long long jB = acc2->Bbits->a[i_long];
    isOi = ~(iA ^ iB);
    isOj = ~(jA ^ jB);
    isXi = ~iA & iB;
    isXj = ~jA & jB;
    unsigned long long isDo = (iA ^ jB) & isOi & isOj;
    unsigned long long isSo = ~(iA ^ jB) & isOi & isOj;
    unsigned long long isDx = (isOi & isXj) | (isOj & isXi);
    unsigned long long isSx = isXi & isXj;
    //   Ndo  i.e. number of markers with both homozyg, but different (i.e. 02 and 20)
      rval.l1 += __builtin_popcountll(isDo);
      //   Nso   00, 22
      rval.l2 += __builtin_popcountll(isSo);
      //    Ndx  01, 10, 12, 21
      rval.l3 += __builtin_popcountll(isDx);
      //    Nsx   11
      rval.l4 += __builtin_popcountll(isSx);
  }
  return rval;
}

void quick_and_dirty_hgmrs(GenotypesSet* the_gtsset){ // get q and d 'hgmr' for all accession pairs
  long good_count = 0;
  long bad_count = 0;
  char ploidy_char = (char)(the_gtsset->ploidy + 48);
  //  fprintf(stderr, "# number of accessions: %ld\n", the_gtsset->accessions->size);
  for(long i=0; i<the_gtsset->accessions->size; i++){
    Accession* acc1 = the_gtsset->accessions->a[i];
    Vlong* alt_homozygs = acc1->alt_homozygs;
    for(long j=i+1; j<the_gtsset->accessions->size; j++){
      Accession* acc2 = the_gtsset->accessions->a[j];
      ND n_d12 = quick_hgmr(acc1, acc2, ploidy_char);
      if(n_d12.d == 1) {bad_count++;}else{good_count++;}
    }
  }
  fprintf(stderr, "# good_count: %ld  bad_count: %ld \n", good_count, bad_count);
}

ND quick_hgmr(Accession* acc1, Accession* acc2, char ploidy_char){
  // get hgmr by considering just the markers with dosage == ploidy or dosage = 0 in acc1
  long n_p0_0p = 0;
  long n_00_pp = 0;
  for(long k=0; k<acc1->alt_homozygs->size; k++){ // dosage = ploidy in acc1
    long idx = acc1->alt_homozygs->a[k];
    char gt2 = acc2->genotypes->a[idx];
    if(gt2 == '0'){ // p0
      n_p0_0p++;
    }else if(gt2 == ploidy_char){ // pp
      n_00_pp++;
    }
  }
  for(long k=0; k<acc1->ref_homozygs->size; k++){ // dosage = 0 in acc1
    long idx = acc1->ref_homozygs->a[k];
    char gt2 = acc2->genotypes->a[idx];
    if(gt2 == '0'){ // 00
      n_00_pp++;
    }else if(gt2 == ploidy_char){ // 0p
      n_p0_0p++;
    }
  }
  ND result = {n_p0_0p, n_00_pp + n_p0_0p}; //, r_numer, r_denom};
  return result;
}

four_longs quick_hgmr_R(Accession* acc1, Accession* acc2, char ploidy_char){
  // get hgmr by considering just the markers with dosage == ploidy or dosage = 0 in acc1
  long n_p0_0p = 0;
  long n_01_p1 = 0;
  long n_00_pp = 0;
  for(long k=0; k<acc1->alt_homozygs->size; k++){ // dosage == ploidy in acc1
    //long idx = acc1->alt_homozygs->a[k];
    char gt2 = acc2->genotypes->a[acc1->alt_homozygs->a[k]]; //  idx];
    if(gt2 == '0'){ // p0
      n_p0_0p++;
    }else if(gt2 == ploidy_char){ // pp
      n_00_pp++;
    }else if(gt2 != MISSING_DATA_CHAR){ // acc2 heterozygous (1 ... ploidy-1)
      n_01_p1++;
    }
  }
  for(long k=0; k<acc1->ref_homozygs->size; k++){ // dosage == 0 in acc1
    //long idx = acc1->ref_homozygs->a[k];
    char gt2 = acc2->genotypes->a[acc1->ref_homozygs->a[k]]; //   idx];
    if(gt2 == '0'){ // 00
      n_00_pp++;
    }else if(gt2 == ploidy_char){ // 0p
      n_p0_0p++;
    }else if(gt2 != MISSING_DATA_CHAR){ // acc2 heterozygous (1 ... ploidy-1)
      n_01_p1++;
    }
  }
  long hgmr_denom = n_00_pp + n_p0_0p;
  long r_numer = n_p0_0p + n_01_p1; // r_numer/r_denom should be small if acc1 is self parent of acc2
  long r_denom = n_00_pp + r_numer;
  four_longs result = {n_p0_0p, hgmr_denom, r_numer, r_denom};
  return result;
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

ND xhgmr(GenotypesSet* gtset, Accession* a1, Accession* a2, int quick){
  // numerator is same as hgmr, but denominator is based on expected numbers of
  // dosage = 0 markers in random accession
  Vlong* a1d2s = a1->alt_homozygs;
  double expected_refds = 0;
  long counted_refds = 0;
  long n0s = 0;
  
  for(long i=0; i<a1d2s->size; i++){ // consider markers with dosage==2 in accession 1
    long idx = a1d2s->a[i];
    char a2_dosage = a2->genotypes->a[idx];
    if(a2_dosage == '0') counted_refds++;
    n0s += gtset->marker_dosage_counts[0]->a[idx]; // n0s_this_marker;
  }
  expected_refds = (double)n0s/(double)gtset->accessions->size;
  if(quick && counted_refds/expected_refds > 0.5  && counted_refds > 100) return (ND){counted_refds, expected_refds};
  
  Vlong* a2d2s = a2->alt_homozygs;
  for(long i=0; i<a2d2s->size; i++){ // consider markers with dosage==2 in accession 2
    long idx = a2d2s->a[i];
    char a1_dosage = a1->genotypes->a[idx];
    if(a1_dosage == '0') counted_refds++;
    n0s += gtset->marker_dosage_counts[0]->a[idx]; // n0s_this_marker;
  }
  expected_refds = (double)n0s/(double)gtset->accessions->size;
  return (ND){counted_refds, expected_refds};
}

four_longs hgmr_R(char* par_gts, char* prog_gts, char ploidy_char){ // return hgmr numerator and denominator
  char c1, c2;
  long i=0;
  long n02n20 = 0;
  long n00n22 = 0;
  long n0xn2x = 0;
  long n01n21 = 0;
  while((c1 = par_gts[i]) != '\0'){
    if((c1 == '0') || (c1 == ploidy_char)){ // c1 homozygous 
      c2 = prog_gts[i];
      if(c2 != MISSING_DATA_CHAR){
	n0xn2x++;
	if((c2 == '0') || (c2 == ploidy_char)){ // c2 homozygous
	  if(c1 != c2) {
	    n02n20++;
	  }else{
	    n00n22++;
	  }
	}else{
	  n01n21++;
	}
      }
    }
    i++;
  }
  four_longs result = {n02n20, n02n20+n00n22, n02n20+n01n21, n0xn2x};
  return result; 
}

ND quick_and_dirty_hgmr(Accession* acc1, Accession* acc2, char ploidy_char){ // get quick 'hgmr', and then if not large get true hgmr.
  long n_p0 = 0; // dosages = ploidy, 0
  long n_phet = 0; // dosages = ploidy, heterozyg
  long n_pp = 0; // both dosages = ploidy
  long n_md = 0; // missing data
  long denom = 0;
  long start = 0;
  long end = 100;

 
  while(1){
    if(end > acc1->alt_homozygs->size) end = acc1->alt_homozygs->size;
    for(long k=start; k<end; k++){
      char gt2 = acc2->genotypes->a[acc1->alt_homozygs->a[k]];
      if(gt2 == '0'){ // ploidy-0
	n_p0++;
      }else if(gt2 == MISSING_DATA_CHAR){
	n_md++; // (missing data)
      }else if(gt2 == ploidy_char){
	n_pp++;
      }else{
	n_phet++;
      }
    }
    denom = n_p0 + n_phet + n_pp;
    double sigma = sqrt((double)n_p0*(n_phet+n_pp)/denom);
    if(denom > 0 &&(n_p0 - 3*sigma > 0.1*denom )) { // quick and dirty 'hgmr' is too big, reject.
      ND result = {1, 1};
      return result; //
    }
    start += 100;
    if(start >= acc1->alt_homozygs->size) break; // have used up all the alt_homozygs
    end += 100;
  }

  // if not rejected, look at acc1 dosage=0 markers to get hgmr
  long n_00 = 0;
  long n_0p = 0;
  for(long k=0; k<acc1->ref_homozygs->size; k++){
    long idx = acc1->ref_homozygs->a[k];
    char gt2 = acc2->genotypes->a[idx];
    if(gt2 == '0'){ // 00
      n_00++;
    }else if(gt2 == ploidy_char){
      n_0p++; // 0p
    }
  }
  
  denom = n_00 + n_0p + n_p0 + n_pp;
  // double h = hgmr(acc1->genotypes->a, acc2->genotypes->a);
   
  //  if(denom == 0) fprintf(stderr, "# B %ld %ld %ld %ld\n", n_00, n_0p, n_p0, n_pp);
  double hgmr = (denom > 0)? (double)(n_0p + n_p0)/(double)denom : -1;
  // fprintf(stderr, "%s %s  %7.5f      %ld %ld %7.5f  \n", acc1->id->a, acc2->id->a,  h, n_p0 + n_0p, denom, hgmr);
  
  //  if(numer < 30) fprintf(stderr, "ijnd: %ld %ld %ld %ld\n", i, j, numer, denom);
  ND result = {n_0p+n_p0, denom};
  return result; // hgmr;
}

ND ghgmr(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny){ // generalized hgmr (ploidy can be > 2; reduces to hgmr for diploid)
  if(parent1 == NULL  ||  progeny == NULL){
    ND result = {0, 0};
    return result;
  }
  long forbidden_count = 0;
  long denom = 0;
  long ploidy = the_gtsset->ploidy;
  for(long i=0; i<the_gtsset->n_markers; i++){
    long a_dosage = parent1->genotypes->a[i];
    long b_dosage = progeny->genotypes->a[i];
    if(a_dosage == MISSING_DATA_CHAR  ||  b_dosage == MISSING_DATA_CHAR) continue;
    a_dosage -= 48; b_dosage -= 48;
    long delta_dosage = labs(a_dosage - b_dosage);
    if(2*delta_dosage > ploidy){
      forbidden_count++;
    }
    else if(2*a_dosage != ploidy  &&  2*b_dosage != ploidy){
      denom++;
    }
  }
  denom += forbidden_count;
  ND result = {forbidden_count, denom};
  return result;
}


ND ghgmr_old(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny){
  if(parent1 == NULL  ||  progeny == NULL){
    ND result = {0, 0};
    return result;
  }
  long forbidden1_count = 0;
  long denom = 0;
  long forbidden2_count = 0;
  // long denom2 = 0;
  // long denom2 = 0;
  long ploidy = the_gtsset->ploidy;
  if(ploidy == 2){
    for(long i=0; i<the_gtsset->n_markers; i++){
      long a_dosage = parent1->genotypes->a[i];
      long b_dosage = progeny->genotypes->a[i];
      if(a_dosage == MISSING_DATA_CHAR  ||  b_dosage == MISSING_DATA_CHAR) continue;
      //   if(ploidy == 2){
      if(a_dosage == '0'){
	if(b_dosage == '2'){
	  forbidden1_count++;
	} else if(b_dosage == '0'){
	  denom++;
	}	
      }else if(a_dosage == '2'){
	if(b_dosage == '0'){
	  forbidden1_count++;
	}else if (b_dosage == '2'){
	  denom++;
	}
      }
    }
  }

 
  if(ploidy == 4){
    for(long i=0; i<the_gtsset->n_markers; i++){
      long a_dosage = parent1->genotypes->a[i];
      long b_dosage = progeny->genotypes->a[i];
      if(a_dosage == MISSING_DATA_CHAR  ||  b_dosage == MISSING_DATA_CHAR) continue;
      if(a_dosage != '2'){
	if(a_dosage == '0'){ // 03, 04 forbidden
	  if(b_dosage == '2'){ // don't count
	  }else if(b_dosage == '3'){
	    forbidden1_count++;
	  } else if(b_dosage == '4'){
	    forbidden2_count++;
	  }else{ // b_dosage is '0' or '1'
	    denom++;
	  }
	}else if(a_dosage == '1'){ // 14 forbidden
	  if(b_dosage == '2'){  // don't count
	  }else if(b_dosage == '4'){
	    forbidden1_count++;
	  }else{ // b_dosage is '0', '1', or '3'
	    denom++;
	  }
	}else if(a_dosage == '3'){ // 30 forbidden
	  if(b_dosage == '0'){
	    forbidden1_count++;
	  }else if(b_dosage != '2'){
	    denom++;	}
	}else if(a_dosage == '4'){ // 40, 41 forbidden
	  if(b_dosage == '1'){
	    forbidden1_count++;
	  }else if(b_dosage == '0'){
	    forbidden2_count++;
	  }else if(b_dosage != '2'){
	    denom++;
	  }
	}
      }
      // *******************************
    }
  }

  long numer = forbidden1_count+forbidden2_count;
  denom += numer;
  //denom2 += forbidden2_count;
  ND result = {numer, denom};
  return result;
} 

void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# max_marker_missing_data_fraction: %8.4lf \n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]);
  }fprintf(fh, "\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    Accession* acc = the_gtsset->accessions->a[i];
    fprintf(fh, "%s  ", acc->id->a);
    for(long j=0; j < acc->genotypes->length; j++){
      char gt = acc->genotypes->a[j];
      if(gt == '~') gt = 'X';
      fprintf(fh, "%c ", gt); // acc->genotypes->a[j]);
    }
    fprintf(fh, "\n");
  }
}

void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(fh, "# max marker missing data fraction: %6.4f\n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  if(the_gtsset == NULL) return;
  free_vstr(the_gtsset->marker_ids);
  free_vlong(the_gtsset->marker_missing_data_counts);
  free_vlong(the_gtsset->marker_alt_allele_counts);
  if(the_gtsset->marker_dosage_counts != NULL){
    for(long i=0; i<=MAX_PLOIDY; i++){
      free_vlong(the_gtsset->marker_dosage_counts[i]);
    }
    free(the_gtsset->marker_dosage_counts);
  }
  free_vaccession(the_gtsset->accessions);
  free(the_gtsset);
}

Vidxid* construct_vidxid(const GenotypesSet* the_gtsset){
  Vidxid* the_vidxid = (Vidxid*)malloc(sizeof(Vidxid));
  the_vidxid->capacity = the_gtsset->n_accessions;
  the_vidxid->size = the_vidxid->capacity;
  the_vidxid->a = (IndexId**)malloc(the_vidxid->size*sizeof(IndexId*));
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId));
    the_idxid->index = i;
    char* the_id =  the_gtsset->accessions->a[i]->id->a;
    the_idxid->id = strcpy((char*)malloc((strlen(the_id)+1)*sizeof(char)), the_id);
    the_vidxid->a[i] = the_idxid;
  }
  return the_vidxid;  
}

Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset){
  Vidxid* the_vidxid = construct_vidxid(the_gtsset);
  sort_vidxid_by_id(the_vidxid);
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* x = the_vidxid->a[i];  }
  if(DBUG) assert(check_idxid_map(the_vidxid, the_gtsset) == 1);
  return the_vidxid;
}


 
/* two_doubles lls(GenotypesSet* the_gtsset, Accession* parent1, Accession* progeny, FILE* stream, double epsilon){ */
/*   // */
/*   if(parent1 == NULL  ||  progeny == NULL){ */
/*     two_doubles result = {0, 0}; */
/*     return result; */
/*   } */
/*   double ll_aob = 0; */
/*   long forbidden_count = 0; */
/*   long count_00_22 = 0; */
/*   double ll_oob = 0; */
/*   double ll_forbidden_oob = 0; */
/*   long ploidy = the_gtsset->ploidy; */
/*   // fprintf(stderr, "## ploidy: %ld \n", ploidy); exit(0); */
/*   long n_markers = the_gtsset->marker_alt_allele_counts->size; */
/*   for(long i=0; i<n_markers; i++){ */
/*     long a_dosage = parent1->genotypes->a[i]; */
/*     long b_dosage = progeny->genotypes->a[i]; */
/*     if(a_dosage == MISSING_DATA_CHAR  ||  b_dosage == MISSING_DATA_CHAR) continue; */
/*     double f = (double)the_gtsset->marker_alt_allele_counts->a[i]/((the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i])*ploidy); */
/*     //   fprintf(stderr, "#  %ld %ld %10.5g \n", n_markers, the_gtsset->marker_missing_data_counts->a[i], f); */
/*     if(ploidy == 2){ */
/*       double pr_oob = pow(f, b_dosage) * pow(1.0-f, ploidy-b_dosage) * ((b_dosage == '1')? 2 : 1); */
/*       if(a_dosage == '0'){ */
/* 	if(b_dosage == '2'){ */
/* 	  forbidden_count++; */
/* 	  ll_forbidden_oob += log(pr_oob); */
/* 	  ll_oob += log(pr_oob); */
/* 	  ll_aob += log(epsilon); */
/* 	} else if(b_dosage == '1'){ */
/* 	  ll_aob += log((f)*(1-epsilon)); */
/* 	  //	    fprintf(stderr, "# f, log(f), ll_aob: %7.5f  %9.5f  %9.5g\n", f, log(f), ll_aob); */
/* 	  ll_oob += log(pr_oob); */
/* 	} else if(b_dosage == '0'){ */
/* 	  ll_aob += log((1-f)*(1-epsilon)); */
/* 	  ll_oob += log(pr_oob); */
/* 	  count_00_22++; */
/* 	}	 */
/*       }else if(a_dosage == '1'){ */
/* 	ll_oob += log(pr_oob); */
/* 	if(b_dosage == '0'){ */
/* 	  ll_aob += log(0.5*(1-f)); */
/* 	}else if(b_dosage == '1'){ */
/* 	  ll_aob += log(0.5); */
/* 	}else if(b_dosage == '2'){ */
/* 	  ll_aob += log(0.5*f); */
/* 	} */
/*       }else if(a_dosage == '2'){ */
/* 	if(b_dosage == '0'){ */
/* 	  forbidden_count++; */
/* 	  ll_forbidden_oob += log(pr_oob); */
/* 	  ll_aob += log(epsilon); */
/* 	}else if(b_dosage == '1'){ */
/* 	  ll_oob += log(pr_oob); */
/* 	  ll_aob += log((1-f)*(1-epsilon)); */
/* 	}else if (b_dosage == '2'){ */
/* 	  ll_oob += log(pr_oob); */
/* 	  ll_aob += log(f*(1-epsilon)); */
/* 	  count_00_22++; */
/* 	} */
/*       } */
/*     }else{ */
/*       fprintf(stderr, "# lls not implemented for ploidy > 2"); */
/*     } // ploidy == 2 */
/*   } // loop over markers */


/*   two_doubles result = {ll_aob, ll_oob}; */
/*   fprintf(stream, " %10.8g %10.8g  ", ll_aob, ll_oob); */
/*   fprintf(stream, " %ld %ld  ", forbidden_count, forbidden_count+count_00_22); */
/*   return result; */
/* } */

