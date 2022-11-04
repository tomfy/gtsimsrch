#include <errno.h>
#include <limits.h>
#include <math.h>
#include "gtset.h"
//#include "various.h"
//#include "pedigree.h"


extern int do_checks_flag; // option -c sets this to 1 to do some checks.

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
  the_accession->ref_homozygs = construct_vlong(1000);
  //  the_accession->heterozygs = construct_vlong(1000);
  the_accession->alt_homozygs = construct_vlong(1000);
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
    add_long_to_vlong(chunk_pats, i_pat);
    
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  the_gts->md_chunk_count = gts_mdchunk_count;
  return gts_mdchunk_count;
}

char* print_accession(Accession* the_gts, FILE* ostream){
  // fprintf(ostream, "Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->genotypes);
  // fprintf(ostream, "%s  %s\n", the_gts->id, the_gts->genotypes);
  fprintf(ostream, "XXX  %s %ld  %ld %ld %ld\n",
	  the_gts->id->a, the_gts->index, the_gts->genotypes->length, the_gts->genotypes->length, the_gts->missing_data_count);
}

void free_accession(Accession* the_accession){
  // print_accession(the_accession, stderr);
  if(the_accession == NULL) return;
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
  free_vlong(the_accession->chunk_patterns);
  free_vlong(the_accession->alt_homozygs);
  free_vlong(the_accession->ref_homozygs);
  free(the_accession);
}
void free_accession_innards(Accession* the_accession){ // doesn't free the_accession itself
  // use to free each element of an array of Accession (not an array of Accession*)
  if(the_accession == NULL) return;
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
}

// *****  Vaccession implementation  *****
Vaccession* construct_vaccession(long cap){ // construct Vaccession with capacity of cap, but empty.
  Vaccession* the_vacc = (Vaccession*)malloc(sizeof(Vaccession));
  the_vacc->capacity = cap;
  the_vacc->size = 0;
  the_vacc->a = (Accession**)malloc(cap*sizeof(Accession*));
}

void add_accession_to_vaccession(Vaccession* the_vacc, Accession* the_acc){
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
  for(long i=0; i<the_vacc->size; i++){
    free_accession(the_vacc->a[i]);
  }
  free(the_vacc->a);
  free(the_vacc);
}

// *****  GenotypesSet implementation *****
GenotypesSet* construct_empty_genotypesset(double max_marker_md_fraction, double min_min_allele_freq, long ploidy){
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  //  the_gtsset->delta = delta;
  the_gtsset->max_marker_missing_data_fraction = max_marker_md_fraction;
  the_gtsset->min_minor_allele_frequency = min_min_allele_freq;
  the_gtsset->n_accessions = -1; // accessions->size;
  the_gtsset->n_bad_accessions = 0;
  the_gtsset->n_markers = -1; // marker_ids->size;
  the_gtsset->ploidy = ploidy;
  the_gtsset->accessions = construct_vaccession(INIT_VACC_CAPACITY);
  the_gtsset->marker_ids = NULL; //marker_ids;
  the_gtsset->marker_missing_data_counts = NULL; // md_counts;
  the_gtsset->marker_alt_allele_counts = NULL; // md_counts;
  the_gtsset->marker_dose_counts = NULL;
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
    /* if(strcmp(token, "ploidy:") == 0){ // can specify ploidy in input this way, but not needed as will be inferred from set of dosages. */
    /*   token = strtok_r(NULL, "\t \n\r", &saveptr); */
    /*   long pldy = str_to_long(token); */
    /*   if(errno != 0){ */
    /* 	fprintf(stderr, "# ploidy not specified in input file %s, bye.\n", input_filename); */
    /* 	exit(EXIT_FAILURE); */
    /*   } */
    /*   if(pldy > the_genotypes_set->ploidy) the_genotypes_set->ploidy = pldy; */
    /*   fprintf(stderr, "# ploidy: %ld %ld\n", pldy, the_genotypes_set->ploidy); */
    /*   continue; */
    /* }else */
    if((strcmp(token, "MARKER") == 0)){
      while(1){
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	if(token == NULL) break;
	char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
	add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
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
    
    the_genotypes_set->marker_dose_counts = (Vlong**)malloc((MAX_PLOIDY+1)*sizeof(Vlong*));
    for(long i=0; i<=MAX_PLOIDY; i++){
      the_genotypes_set->marker_dose_counts[i] = construct_vlong_zeroes(marker_ids->size);
    }
  }
  
  // if the_genotypes_set already has marker_ids set, check for agreement with marker_ids
  // i.e. in case of using reference set, make sure the number of markers is the same in both data sets.
  if(the_genotypes_set->marker_ids  == NULL){
    the_genotypes_set->marker_ids = marker_ids;
  }else{
    if(the_genotypes_set->marker_ids->size != marker_ids->size){
      fprintf(stderr, "# data sets have different numbers of markers: %5ld  %5ld. Exiting.\n", the_genotypes_set->marker_ids->size, marker_ids->size);
      exit(EXIT_FAILURE);
    }
    else if(compare_vstrs(marker_ids, the_genotypes_set->marker_ids) != 0){
      fprintf(stderr, "# sets of marker ids do not agree. Exiting.\n");
      exit(EXIT_FAILURE);
    }
    free_vstr(marker_ids);
  }

  // Read in the rest of the lines, construct Accession for each line
  long accession_count = 0;
  while((nread = getline(&line, &len, g_stream)) != -1){
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
      
      //  fprintf(stdout, "token, ploidy: %s  %ld\n", token, the_genotypes_set->ploidy); // marker_count, genotypes[marker_count]);
      if(genotypes[marker_count] == MISSING_DATA_CHAR){
	the_genotypes_set->marker_missing_data_counts->a[marker_count]++;
	accession_missing_data_count++;
      }else{
	the_genotypes_set->marker_alt_allele_counts->a[marker_count] += (long)(genotypes[marker_count]-48);
      }

      marker_count++;
    } // done reading dosages for all markers of this accession
   
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    Accession* the_accession = construct_accession(acc_id, accession_count, genotypes, accession_missing_data_count);
    free(acc_id); // or cut out the middleman (acc_id)?
    free(genotypes);
    //  the_accession->missing_data_count = accession_missing_data_count;
    if(accession_missing_data_count <= max_acc_missing_data_fraction * the_genotypes_set->marker_ids->size){
      add_accession_to_vaccession(the_genotypes_set->accessions, the_accession);
      accession_count++;
    }else{
      fprintf(stderr, "# Accession: %s rejected due to missing data at %ld out of %ld markers.\n",
	      the_accession->id->a, accession_missing_data_count, the_genotypes_set->marker_ids->size);
      the_genotypes_set->n_bad_accessions++;
    }
  
   
  } // done reading all lines
  fclose(g_stream);
  free(line); // only needs to be freed once.

  the_genotypes_set->n_accessions = the_genotypes_set->accessions->size;
  the_genotypes_set->n_markers = the_genotypes_set->marker_ids->size;
  if(DBUG && do_checks_flag) check_genotypesset(the_genotypes_set); 
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
  if(strcmp(token,"NA") == 0){ // missing data
    return MISSING_DATA_CHAR;
  }else{
    long l = str_to_long(token);
    //  fprintf(stderr, "# l: %ld  errno: %ld \n", l, (long)errno);
    //   fprintf(stderr, "#l: %ld   ploidy: %ld  errno: %ld\n", l, ploidy, (long)errno);
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
  for(long j=0; j<the_gtsset->n_accessions; j++){
    Accession* the_acc = the_gtsset->accessions->a[j];
    for(long i=0; i<the_gtsset->n_markers; i++){
      long dosage = (long)the_acc->genotypes->a[i];
      if(dosage != MISSING_DATA_CHAR){
	the_gtsset->marker_dose_counts[dosage-48]->a[i]++;
      }
    }
  }
}

double ragmr(GenotypesSet* the_gtsset){
  for(long j=0; j<the_gtsset->n_accessions; j++){
    Accession* the_acc = the_gtsset->accessions->a[j];
    for(long i=0; i<the_gtsset->n_markers; i++){
      long dosage = (long)the_acc->genotypes->a[i];
      if(dosage != MISSING_DATA_CHAR){
	the_gtsset->marker_dose_counts[dosage-48]->a[i]++;
      }
    }
  }

  double ragmr = 0;
  for(long i=0; i<the_gtsset->n_markers; i++){
    double sum_nsq = 0;
    double sum_n = 0;
    for(long k=0; k<=the_gtsset->ploidy; k++){
      double n = (double)the_gtsset->marker_dose_counts[k]->a[i];
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


/* void read_genotypes_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set){ // */

/*   FILE* g_stream = fopen(input_filename, "r"); */
/*   if (g_stream == NULL) { */
/*     perror("fopen"); */
/*     exit(EXIT_FAILURE); */
/*   } */
  
/*   char* line = NULL; */
/*   size_t len = 0; */
/*   ssize_t nread; */

/*   long accessions_capacity = 2000; */
/*   Vaccession* the_accessions = construct_vaccession(accessions_capacity); //(Vaccession*)malloc(accessions_capacity*sizeof(ccession)); */

/*   // read delta from first (comment) line */
/*   double delta; */
/*   double max_marker_missing_data_fraction; */
/*   char* saveptr = NULL; */
/*   if((nread = getline(&line, &len, g_stream)) != -1){ */
/*     char* token = strtok_r(line, "\t \n\r", &saveptr); */
/*     if((token == NULL)  || (strcmp(token, "#") != 0)){ */
/*       exit(EXIT_FAILURE); */
/*     } */
/*     long field_count = 1; */
/*     while(1){ */
/*       token = strtok_r(NULL, "\t \n\r", &saveptr); */
/*       if(token == NULL) break; */
/*       if(field_count == 2){ */
/* 	delta = atof(token); */
/*       }else if(field_count == 4){ */
/* 	max_marker_missing_data_fraction = atof(token); */
/*       } */
/*     } */
/*   } */
/*   // ********************************************************* */
/*   // read first non-comment line, with MARKER and marker ids */
/*   long markerid_count = 0; */
/*   Vstr* marker_ids = construct_vstr(1000); */
/*   saveptr = NULL; */
/*   if((nread = getline(&line, &len, g_stream)) != -1){ */
/*     char* token = strtok_r(line, "\t \n\r", &saveptr); */
/*     if((token == NULL)  || (strcmp(token, "MARKER") != 0)){ */
/*       exit(EXIT_FAILURE); */
/*     } */
/*     while(1){ */
/*       token = strtok_r(NULL, "\t \n\r", &saveptr); */
/*       if(token == NULL) break; */
/*       char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char)); */
/*       add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store */
/*       markerid_count++; */
/*     } */
/*   } */
/*   // ********************************************************* */
/*   if(the_genotypes_set->marker_missing_data_counts == NULL){ */
/*     the_genotypes_set->marker_missing_data_counts = construct_vlong_zeroes(marker_ids->size); */
/*   } */

/*   long accession_count = 0; */
/*   // long* marker_missing_data_counts = (long*)calloc(markerid_count, sizeof(long)); */
/*   // Vlong* the_md_vlong = construct_vlong_zeroes(markerid_count); // construct_vlong_from_array(markerid_count, marker_missing_data_counts); */
/*   saveptr = NULL; */
/*   while((nread = getline(&line, &len, g_stream)) != -1){ */
  
/*     char* token = strtok_r(line, "\t \n\r", &saveptr); */
/*     char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token); */
/*     long token_count = 1; */
/*     long accession_missing_data_count = 0; */
/*     char* genotypes = (char*)malloc((markerid_count+1) * sizeof(char)); */
/*     genotypes[markerid_count] = '\0'; // terminate with null. */
/*     while(1){ */
/*       token = strtok_r(NULL, "\t \n\r", &saveptr); */
/*       if(token == NULL) break; */
/*       token_count++; */
/*       genotypes = strcpy(genotypes, token); */
/*     } */
/*     assert(token_count == 2); // because line should be: accession_id 00102011000201020....\n  just 2 whitespace-separated tokens. */
/*     assert(strlen(genotypes) == markerid_count); */
/*     for(long igt = 0; igt < markerid_count; igt++){ // count missing data      */
/*       if(genotypes[igt] == MISSING_DATA_CHARACTER){ */
/* 	the_genotypes_set->marker_missing_data_counts->a[igt]++; */
/* 	accession_missing_data_count++; */
/*       } */
/*     } // done reading genotypes for all markers */
  
/*     Accession* the_accession = construct_accession( acc_id, accession_count, genotypes, accession_missing_data_count); */
/*     free(acc_id); // or cut out the middleman (acc_id)? */
/*     free(genotypes); */
/*     //   the_accession->missing_data_count = accession_missing_data_count; */
/*     add_accession_to_vaccession(the_accessions, the_accession); */
/*     accession_count++; */
/*   } // done reading all lines */
/*   fclose(g_stream); */
/*   free(line); // only needs to be freed once. */
/*   the_genotypes_set->accessions = the_accessions; */
/*   //return the_genotypes_set; */
/* } */

void check_genotypesset(GenotypesSet* gtss){
  assert(gtss->marker_ids->size == gtss->n_markers);
  long* md_counts = (long*)calloc(gtss->n_markers, sizeof(long)); 
  for(long i=0; i<gtss->n_accessions; i++){
    Accession* an_acc = gtss->accessions->a[i];
    assert(strlen(an_acc->genotypes->a) == gtss->n_markers);
    for(long j=0; j<gtss->n_markers; j++){
      if(an_acc->genotypes->a[j] == MISSING_DATA_CHAR)md_counts[j]++;
    }
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(md_counts[j] == gtss->marker_missing_data_counts->a[j]);
  }
  free(md_counts);
  fprintf(stderr, "# Successfully completed check_genotypesset\n");
}

void clean_genotypesset(GenotypesSet* the_gtsset){ // construct a new set of 'cleaned' accessions, which replace the raw accessions
  double max_marker_md_fraction = the_gtsset->max_marker_missing_data_fraction;
  Vlong* md_counts = the_gtsset->marker_missing_data_counts; // the number of missing data for each marker
  Vlong* alt_allele_counts = the_gtsset->marker_alt_allele_counts;
  
  long n_accs = the_gtsset->n_accessions;
  if(max_marker_md_fraction < 0){ //set max_marker_md_fraction to some multiple of median md fraction.
    long* mdcount_histogram = (long*)calloc(n_accs+1, sizeof(long));
    double factor = -1*max_marker_md_fraction;
    for(long i=0; i< md_counts->size; i++){
      mdcount_histogram[md_counts->a[i]]++;
    }
    long mrkrs_so_far = 0;
    long median_md_count;
    for(long i=0; i<=n_accs; i++){
      mrkrs_so_far += mdcount_histogram[i];
      if(mrkrs_so_far > 0.5*md_counts->size){
	median_md_count = i;
	break;
      }
    }
    max_marker_md_fraction = factor*(double)median_md_count/(double)n_accs;
  }
  // identify the markers to keep:
  long n_markers_to_keep = 0;
  Vlong* md_ok = construct_vlong_zeroes(md_counts->size); // set to 1 if number of missing data gts is small enough
  Vstr* cleaned_marker_ids = construct_vstr(1000);
  Vlong* cleaned_md_counts = construct_vlong(1000);
  Vlong* cleaned_alt_allele_counts = construct_vlong(1000);
  long mdsum_all = 0;
  long mdsum_kept = 0;
  long altallelesum_all = 0;
  long altallelesum_kept = 0;
  
  for(long i=0; i<md_counts->size; i++){
    mdsum_all += md_counts->a[i];
    altallelesum_all += alt_allele_counts->a[i];
    /*  fprintf(stderr, "mmaf: %10.5f  %ld %ld  %ld  %8.4f  %ld  %8.4f\n",
	the_gtsset->min_minor_allele_frequency, the_gtsset->ploidy,
	the_gtsset->n_accessions, alt_allele_counts->a[i],
	the_gtsset->ploidy * the_gtsset->n_accessions * the_gtsset->min_minor_allele_frequency, md_counts->a[i],
	max_marker_md_fraction*the_gtsset->n_accessions );  /* */
    double min_min_allele_count = the_gtsset->ploidy * (the_gtsset->n_accessions - md_counts->a[i]) * the_gtsset->min_minor_allele_frequency;
    double max_min_allele_count = the_gtsset->ploidy * (the_gtsset->n_accessions - md_counts->a[i]) * (1.0 - the_gtsset->min_minor_allele_frequency);
    if(
       (md_counts->a[i] <= max_marker_md_fraction*the_gtsset->n_accessions) // not too much missing data
       &&
       ( (alt_allele_counts->a[i] >= min_min_allele_count)  && // minor allele frequency not too small,
	 (alt_allele_counts->a[i] <= max_min_allele_count) ) // (or too large)
       ) {
      md_ok->a[i] = 1;
      n_markers_to_keep++;
      mdsum_kept += md_counts->a[i];
      altallelesum_kept += alt_allele_counts->a[i];
      long marker_id_length = strlen(the_gtsset->marker_ids->a[i]);
      char* marker_id_to_keep = strcpy((char*)malloc((marker_id_length+1)*sizeof(char)), the_gtsset->marker_ids->a[i]);
      add_string_to_vstr(cleaned_marker_ids, marker_id_to_keep); // store the kept marker ids.
      add_long_to_vlong(cleaned_md_counts, md_counts->a[i]); // store the md counts for the kept markers.
      add_long_to_vlong(cleaned_alt_allele_counts, alt_allele_counts->a[i]);
    }
  }
  double raw_md_fraction = (double)mdsum_all/(double)(md_counts->size*n_accs);
  double cleaned_md_fraction = (double)mdsum_kept/(double)(n_markers_to_keep*n_accs);
  double raw_minor_allele_freq = (double)altallelesum_all/(md_counts->size*n_accs*the_gtsset->ploidy);
  double cleaned_minor_allele_freq = (double)altallelesum_kept/(double)(n_markers_to_keep*n_accs*the_gtsset->ploidy);
  fprintf(stdout, "# Removing markers with missing data fraction > %5.3lf or minor allele frequency < %5.3f\n",
	  max_marker_md_fraction, the_gtsset->min_minor_allele_frequency);
  fprintf(stdout, "# Raw data has %ld markers, missing data fraction = %5.3lf, minor allele frequencey = %5.3f\n",
	  md_counts->size, raw_md_fraction, raw_minor_allele_freq);
  fprintf(stdout, "# Cleaned data has %ld markers, missing data fraction = %5.3lf, minor allele frequency = %5.3lf\n",
	  n_markers_to_keep, cleaned_md_fraction, cleaned_minor_allele_freq);

  Vaccession* the_accessions = construct_vaccession(the_gtsset->n_accessions); //(Accession*)malloc(the_gtsset->n_accessions*sizeof(Accession)); 
  for(long i=0; i<the_gtsset->n_accessions; i++){ // loop over accessions
    char* raw_gts = the_gtsset->accessions->a[i]->genotypes->a; // the string with all the genotypes for accession i
    char* cleaned_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    long acc_md_count = 0;
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      // fprintf(stderr, "j: %ld  gt: %d \n", j, raw_gts[j]);
      if(md_ok->a[j] == 1){
	cleaned_gts[k] = raw_gts[j];
	k++;
	if(raw_gts[j] == MISSING_DATA_CHAR) acc_md_count++;
      }
    }
    cleaned_gts[k] = '\0'; // terminate with null.
    if(DBUG && do_checks_flag) assert(k == n_markers_to_keep);
    // fprintf(stderr, "cleaned genotypes: %s \n", cleaned_gts);
    Accession* the_accession = construct_accession( the_gtsset->accessions->a[i]->id->a, i, cleaned_gts, acc_md_count ); //(Accession*)malloc(sizeof(Accession));
    free(cleaned_gts);
    //   the_accession->missing_data_count = acc_md_count;
    add_accession_to_vaccession(the_accessions, the_accession);
    //   fprintf(stderr, "xx: %ld %s\n", the_accession->index, the_accession->id->a);
  }
  free_vlong(md_ok);

  free_vaccession(the_gtsset->accessions);
  the_gtsset->accessions = the_accessions;
  free_vstr(the_gtsset->marker_ids);
  the_gtsset->marker_ids = cleaned_marker_ids;
  free_vlong(the_gtsset->marker_missing_data_counts);
  the_gtsset->marker_missing_data_counts = cleaned_md_counts;
  free_vlong(the_gtsset->marker_alt_allele_counts);
  the_gtsset->marker_alt_allele_counts = cleaned_alt_allele_counts;
  the_gtsset->n_markers = the_gtsset->accessions->a[0]->genotypes->length;
} // end clean_genotypesset

void rectify_markers(GenotypesSet* the_gtsset){ // if alt allele has frequency > 0.5, dosage -> ploidy-dosage
  long n_accessions = the_gtsset->accessions->size;
  long n_markers = the_gtsset->marker_alt_allele_counts->size;
  long ploidy = the_gtsset->ploidy;
  // fprintf(stderr, "n_acc: %ld n_markers: %ld ploidy: %ld \n", n_accessions, n_markers, ploidy);
  for(long i=0; i<n_markers; i++){
    //  fprintf(stderr, "i: %ld altcount  %ld   md count: %ld \n", i, the_gtsset->marker_alt_allele_counts->a[i], the_gtsset->marker_missing_data_counts->a[i]);
    if(the_gtsset->marker_alt_allele_counts->a[i] > (n_accessions - the_gtsset->marker_missing_data_counts->a[i])){
      for(long j=0; j<n_accessions; j++){
	char dosage = the_gtsset->accessions->a[j]->genotypes->a[i];
	if(dosage != MISSING_DATA_CHAR){
	  the_gtsset->accessions->a[j]->genotypes->a[i] = (ploidy - (dosage - 48)) + 48;
	}
      }
    }
  }
}

void store_homozygs(GenotypesSet* the_gtsset){ // 
  for(long i=0; i<the_gtsset->accessions->size; i++){
    Accession* acc = the_gtsset->accessions->a[i];
   
    for(long j=0; j<the_gtsset->marker_alt_allele_counts->size; j++){
      if(acc->genotypes->a[j] == MISSING_DATA_CHAR) continue;
      long dosage = (long)(acc->genotypes->a[j] - 48); // 0, 1, ..., ploidy
      if( dosage == the_gtsset->ploidy){ // if dosage == ploidy
	add_long_to_vlong(acc->alt_homozygs, j);
	
	//	 fprintf(stderr, "i %ld  dosage: %ld  ploidy: %ld \n", i, (long)acc->genotypes->a[j] - 48, the_gtsset->ploidy);
      }else if(dosage == 0){
	add_long_to_vlong(acc->ref_homozygs, j);
      }
      /* else{ */
      /* 	add_long_to_vlong(acc->heterozygs, j); */
      /* } */
    }
  }
}

void quick_and_dirty_hgmrs(GenotypesSet* the_gtsset){ // get q and d 'hgmr' for all accession pairs
  long good_count = 0;
  long bad_count = 0;
  char ploidy_char = (char)(the_gtsset->ploidy + 48);
  fprintf(stderr, "# number of accessions: %ld\n", the_gtsset->accessions->size);
  for(long i=0; i<the_gtsset->accessions->size; i++){
    Accession* acc1 = the_gtsset->accessions->a[i];
    Vlong* alt_homozygs = acc1->alt_homozygs;
    for(long j=i+1; j<the_gtsset->accessions->size; j++){
      Accession* acc2 = the_gtsset->accessions->a[j];
      //   ND n_d12 = quick_and_dirty_hgmr(acc1, acc2, ploidy_char);
      ND n_d12 = quick_hgmr(acc1, acc2, ploidy_char);
       
      //	four_longs llll = quick_hgmr_R(acc1, acc2, ploidy_char);
      //	if(llll.l2 == 1) {bad_count++;}else{good_count++;}
      //	ND n_d12 = hgmr_nd(acc1->genotypes->a, acc2->genotypes->a, ploidy_char);
      if(n_d12.d == 1) {bad_count++;}else{good_count++;}
      /* long n = n_d12.l1; */
      /* long d = n_d12.l2; */
      /* long rn = n_d12.l3; */
      /* long rd = n_d12.l4; */
      /* fprintf(stderr, " %ld %ld  %ld %ld  %7.5f ", i, j, n, d, (d>0)? (double)n/d : 2.0); */
      /* fprintf(stderr, " %ld %ld  %7.5f", rn, rd, (rd>0)? (double)rn/rd : 2.0); */
      /* fprintf(stderr, "\n"); */
    }
  }
  fprintf(stderr, "# good_count: %ld  bad_count: %ld \n", good_count, bad_count);
}

ND quick_hgmr(Accession* acc1, Accession* acc2, char ploidy_char){
  // get hgmr by considering just the markers with dosage == ploidy or dosage = 0 in acc1
  long n_p0_0p = 0;
  //  long n_01_p1 = 0;
  long n_00_pp = 0;
  for(long k=0; k<acc1->alt_homozygs->size; k++){ // dosage = ploidy in acc1
    long idx = acc1->alt_homozygs->a[k];
    char gt2 = acc2->genotypes->a[idx];
    if(gt2 == '0'){ // p0
      n_p0_0p++;
    }else if(gt2 == ploidy_char){ // pp
      n_00_pp++;
    }
    /* else if(gt2 != MISSING_DATA_CHAR){ // acc2 heterozygous (1 ... ploidy-1) */
    /*   n_01_p1++; */
    /* } */
  }
  for(long k=0; k<acc1->ref_homozygs->size; k++){ // dosage = 0 in acc1
    long idx = acc1->ref_homozygs->a[k];
    char gt2 = acc2->genotypes->a[idx];
    if(gt2 == '0'){ // 00
      n_00_pp++;
    }else if(gt2 == ploidy_char){ // 0p
      n_p0_0p++;
    }
    /* else if(gt2 != MISSING_DATA_CHAR){ // acc2 heterozygous (1 ... ploidy-1) */
    /*   n_01_p1++; */
    /* } */
  }
  //long hgmr_denom = ;
  //  long r_numer = n_p0_0p + n_01_p1;
  //long r_denom = n_00_pp + r_numer;
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
  long r_numer = n_p0_0p + n_01_p1;
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
  //  ND result = {n_numer, n_denom};
  //return result;
  //  fprintf(stderr, "hgmr   n,d: %ld %ld  ", n_numer, n_denom); 
  return (n_denom > 0)? (double)n_numer/(double)n_denom : 2.0;  
}

ND xhgmr(GenotypesSet* gtset, Accession* a1, Accession* a2){
  // numerator is same as hgmr, but denominator is based on expected numbers of
  // dosage = 0 markers in random accession
  Vlong* a1d2s = a1->alt_homozygs;
  // fprintf(stderr, "a1 acc id: %s   alt_homozygs: %ld \n", a1->id->a, a1d2s->size); // getchar();
  double expected_refds = 0;
  long counted_refds = 0;
  
  for(long i=0; i<a1d2s->size; i++){
    long idx = a1d2s->a[i];
    char a2_dosage = a2->genotypes->a[idx];
    if(a2_dosage == '0') counted_refds++;
    long n0s_this_marker = gtset->marker_dose_counts[0]->a[idx];
    long n012s_this_marker = gtset->accessions->size - gtset->marker_missing_data_counts->a[idx];
    expected_refds += (double)n0s_this_marker / (double)n012s_this_marker;
  }
  
  Vlong* a2d2s = a2->alt_homozygs;
  for(long i=0; i<a2d2s->size; i++){
    long idx = a2d2s->a[i];
    char a1_dosage = a1->genotypes->a[idx];
    if(a1_dosage == '0') counted_refds++;
    long n0s_this_marker = gtset->marker_dose_counts[0]->a[idx];
    long n012s_this_marker = gtset->accessions->size - gtset->marker_missing_data_counts->a[idx];
    expected_refds += (double)n0s_this_marker / (double)n012s_this_marker;
  }
  
  // fprintf(stderr, "%ld %ld  %ld %8.4lf  %ld %ld\n", a1->index, a2->index, counted_refds, expected_refds, a1d2s->size, a2d2s->size);
  ND result = {counted_refds, expected_refds};
  return result;
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

  /*   else if(ploidy == 8){ */
  /*     if(a_dosage == '0'){ // 05,06,07,08 forbidden */
  /* 	if(b_dosage-48 >= 5){ */
  /* 	  forbidden1_count++; */
  /* 	}else{ */
  /* 	  denom++; */
  /* 	} */
  /*     }else if(a_dosage == '1'){ // 16, 17, 18 forbidden */
  /* 	if(b_dosage-48 >= 6){ */
  /* 	  forbidden1_count++; */
  /* 	}else{ */
  /* 	  denom++; */
  /* 	} */
  /*     }else if(a_dosage == '2'){ // 27, 28 forbidden */
  /* 	if(b_dosage-48 >=7){ */
  /* 	  forbidden1_count++; */
  /* 	}else{ */
  /* 	  denom++; */
  /* 	} */
  /*     }else if(a_dosage == '3'){ // 38 forbidden */
  /* 	if(b_dosage == '8'){ */
  /* 	  forbidden1_count++; */
  /* 	}else{ */
  /* 	  denom++; */
  /* 	} */
  /*     } */
  /*   }else{ */
  /*     fprintf(stderr, "# lls not implemented for ploidy == %ld", ploidy); */
  /*   }  */
  /* } // loop over markers */
  long numer = forbidden1_count+forbidden2_count;
  denom += numer;
  //denom2 += forbidden2_count;
  ND result = {numer, denom};
  return result;
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

void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# max_marker_missing_data_fraction: %8.4lf \n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]);
  }fprintf(fh, "\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    fprintf(fh, "%s  %s\n", the_gtsset->accessions->a[i]->id->a, the_gtsset->accessions->a[i]->genotypes->a);
  }
}

void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  // fprintf(fh, "# delta: %6.4lf\n", the_gtsset->delta);
  fprintf(fh, "# max marker missing data fraction: %6.4f\n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "# n_markers: %ld\n", the_gtsset->n_markers);
  /* for(long i=0; i<the_gtsset->accessions->size; i++){ */
  /*   fprintf(stderr, "acc: %ld dosage=ploidy_count: %ld \n", i, the_gtsset->accessions->a[i]->alt_homozygs->size); */
  /* } */
}

void free_genotypesset(GenotypesSet* the_gtsset){
  if(the_gtsset == NULL) return;
  free_vstr(the_gtsset->marker_ids);
  free_vlong(the_gtsset->marker_missing_data_counts);
  free_vlong(the_gtsset->marker_alt_allele_counts);
  for(long i=0; i<=MAX_PLOIDY; i++){
    free_vlong(the_gtsset->marker_dose_counts[i]);
  }
  free(the_gtsset->marker_dose_counts);
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
    IndexId* x = the_vidxid->a[i];
    // fprintf(stderr, "%ld %ld %s\n", i, x->index, x->id);
  }
  if(DBUG) assert(check_idxid_map(the_vidxid, the_gtsset) == 1);
  return the_vidxid;
}

