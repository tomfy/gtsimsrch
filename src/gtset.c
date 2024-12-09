#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include "gtset.h"
//#include "various.h"
//#include "pedigree.h"

#define MIN_ALT_ALLELES 0
#define NO_PHASE_CHAR 'x' // used for homozygous, missing data, or phase unknown

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
Accession* construct_accession(char* id, long idx, char* genotypes, char* phases, long accession_md_count){
  Accession* the_accession = (Accession*)malloc(1*sizeof(Accession));
  the_accession->id = construct_vchar_from_str(id);
  the_accession->index = idx;  
  the_accession->genotypes = construct_vchar_from_str(genotypes);
  //fprintf(stderr, "genotypes->length: %ld    ", the_accession->genotypes->length);
   the_accession->phases = construct_vchar_from_str(phases);
   //fprintf(stderr, "phases->length: %ld\n", the_accession->phases->length);
  the_accession->chunk_patterns = NULL; // set to NULL to avoid containing garbage which causes crash when freeing accession
  the_accession->missing_data_count = accession_md_count;
  the_accession->ref_homozygs = NULL; // construct_vlong(10); // can construct in store_homozygs (so no mem used if not calling store_homozygs)
  //  the_accession->heterozygs = construct_vlong(1000);
  the_accession->alt_homozygs = NULL; // construct_vlong(10);
  the_accession->Abits = NULL;
  the_accession->Bbits = NULL;
  the_accession->Fpar_idx = ID_NA_INDEX;
  the_accession->Mpar_idx = ID_NA_INDEX;
  the_accession->has_pedigree = false;
  the_accession->search_done = false;
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
  // the_gts->md_chunk_count = gts_mdchunk_count;
  the_gts->ok_chunk_count = n_chunks - gts_mdchunk_count;
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
  the_gtsset->phased = false;
  the_gtsset->accessions = construct_vaccession(INIT_VACC_CAPACITY);
  the_gtsset->marker_ids = NULL;
  the_gtsset->chromosomes = NULL;
  the_gtsset->marker_missing_data_counts = NULL;
  the_gtsset->marker_alt_allele_counts = NULL;
  the_gtsset->mafs = NULL;
  the_gtsset->marker_dosage_counts = NULL;
  the_gtsset->acc_filter_info = construct_vchar(2000);
  the_gtsset->marker_filter_info = construct_vchar(2000);
  return the_gtsset;
}

void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set, double max_acc_missing_data_fraction, long Nthreads){
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
	// fprintf(stderr, "### %30s  %ld %ld \n", mrkr_id, markerid_count, marker_ids->size); 
      }
      break;
    }else{ 
      fprintf(stderr, "token: %s (should be MARKER)\n", token);
      exit(EXIT_FAILURE);
    }
  }

   while((nread = getline(&line, &len, g_stream)) != -1){

 // fprintf(stderr, "# reading accession %ld\n", accession_count);
    saveptr = line; // these are char*
    char* token = strtok_r(saveptr, "\t \n\r", &saveptr); // read first token in line, either "CHROMOSOME" or and accession id
    if((token == NULL) || (token[0] == '#')) continue; // skip comments, empty lines
    //fprintf(stderr, "token: %s\n", token); getc(stdin);
    if((strcmp(token, "CHROMOSOME") == 0)){ // CHROMOSOME line is present, store chromosome info (needed for phased analysis)
      Vlong* chromosome_numbers = construct_vlong(1000);
      while(1){
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	if(token == NULL) break;
	long i_chrom = str_to_long(token);
	push_to_vlong(chromosome_numbers, i_chrom);
	//	char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
	//	push_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
	//markerid_count++;
	// fprintf(stderr, "### %30s  %ld %ld \n", mrkr_id, markerid_count, marker_ids->size); 
      }
      the_genotypes_set->chromosomes = chromosome_numbers;
      fprintf(stderr, "in add_accession_... size of chromosomes vlong. %ld \n", the_genotypes_set->chromosomes->size);
      break;
    }else{ 
      fprintf(stderr, "token: %s (should be CHROMOSOME)\n", token);
      exit(EXIT_FAILURE);
    }
  }

  
  // *****  done reading first two lines (with marker ids, chromosome numbers)  *****
  if(the_genotypes_set->marker_missing_data_counts == NULL){    
    the_genotypes_set->marker_missing_data_counts = construct_vlong_zeroes(marker_ids->size);
    the_genotypes_set->marker_alt_allele_counts = construct_vlong_zeroes(marker_ids->size);
    //  fprintf(stderr, "XXXXXXXXXXXXXX: %ld \n", the_genotypes_set->marker_missing_data_counts->a[0]);
    
    the_genotypes_set->marker_dosage_counts = (Vlong**)malloc((MAX_PLOIDY+1)*sizeof(Vlong*));
    for(long i=0; i<=MAX_PLOIDY; i++){
      the_genotypes_set->marker_dosage_counts[i] = construct_vlong_zeroes(marker_ids->size);
    }
  }
  
  // if the_genotypes_set already has marker_ids set (from a reference set), check for agreement with marker_ids
  // i.e. in case of using reference set, make sure the set of markers is the same in both data sets.
  if(the_genotypes_set->marker_ids  == NULL){
    the_genotypes_set->marker_ids = marker_ids;
    // the_genotypes_set->chromosomes = chromosome_numbers;
    // fprintf(stderr, "the_genotypes_set->marker_ids->size: %ld \n", the_genotypes_set->marker_ids->size);
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

  double t1 = hi_res_time();
  if(Nthreads == -1){ // old, unthreaded way

    /*
    // read through next non-comment line
    // if CHROMOSOME line, store chromosomes, set phased to true
    // else read & store one line of dosages

      while((nread = getline(&line, &len, g_stream)) != -1){
    // fprintf(stderr, "# reading accession %ld\n", accession_count);
  saveptr = line;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    
    if((token == NULL) || (token[0] == '#')) continue; // skip comments, empty lines
    if((strcmp(token, "CHROMOSOME") == 0)){ // CHROMOSOME line is present; data is phased
        Vlong* chromosome_numbers = construct_vlong(1000);
       while(1){
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	if(token == NULL) break;
	long i_chrom = str_to_long(token);
	push_to_vlong(chromosome_numbers, i_chrom);
	//	char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
	//	push_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
	//markerid_count++;
	// fprintf(stderr, "### %30s  %ld %ld \n", mrkr_id, markerid_count, marker_ids->size); 
      }
       the_genotypes_set->chromosomes = chromosome_numbers;
       the_genotypes_set->phased = true;
  }
      } /* */

    
  // Read in the rest of the lines, construct an Accession for each line
  long accession_count = 0;
  while((nread = getline(&line, &len, g_stream)) != -1){
    // fprintf(stderr, "# reading accession %ld\n", accession_count);
    saveptr = line; // these are char*
    char* token = strtok_r(saveptr, "\t \n\r", &saveptr); // read first token in line, either "CHROMOSOME" or and accession id
    if((token == NULL) || (token[0] == '#')) continue; // skip comments, empty lines
    if(0){
    /* fprintf(stderr, "token: %s\n", token); getc(stdin); */
    /* if((strcmp(token, "CHROMOSOME") == 0)){ // CHROMOSOME line is present, store chromosome info (needed for phased analysis) */
    /*   Vlong* chromosome_numbers = construct_vlong(1000); */
    /*   while(1){ */
    /* 	token = strtok_r(NULL, "\t \n\r", &saveptr); */
    /* 	if(token == NULL) break; */
    /* 	long i_chrom = str_to_long(token); */
    /* 	push_to_vlong(chromosome_numbers, i_chrom); */
    /* 	//	char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char)); */
    /* 	//	push_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store */
    /* 	//markerid_count++; */
    /* 	// fprintf(stderr, "### %30s  %ld %ld \n", mrkr_id, markerid_count, marker_ids->size);  */
    /*   } */
    /*   the_genotypes_set->chromosomes = chromosome_numbers; */
    /*   fprintf(stderr, "in add_accession_... size of chromosomes vlong. %ld \n", the_genotypes_set->chromosomes->size); */
    /*   // the_genotypes_set->phased = true; // now allow for presence of chromosome numbers in unphased data. */
     }else{ // not CHROMOSOME line, should be acc id, followed by dosages */
     
    //char* token = strtok_r(saveptr, "\t \n\r", &saveptr);
      char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
      long marker_count = 0;
      long accession_missing_data_count = 0;
      char* genotypes = (char*)calloc((markerid_count+1), sizeof(char));
      genotypes[markerid_count] = '\0'; // terminate with null.
      char* phases =  (char*)calloc((markerid_count+1), sizeof(char));
      phases[markerid_count] = '\0'; // indices 0 through markerid_count-1 are phases ('x', 'p', or 'm'), index markerid_count has terminating null.

      while(1){ // read dosages from one line.
	//fprintf(stderr, "before saveptr: [%s]\n", saveptr);
	token = strtok_r(NULL, "\t \n\r", &saveptr);
	//fprintf(stderr, "saveptr: [%s]  [%s]\n", saveptr, token);
	if(token == NULL)	break;
      two_chars dsg_ph = token_to_dosage(token, &(the_genotypes_set->ploidy));
      genotypes[marker_count] = dsg_ph.ch1;
      //fprintf(stderr, "dsg_ph.ch2: %c\n", dsg_ph.ch2);
      if(dsg_ph.ch2 != NO_PHASE_CHAR) the_genotypes_set->phased = true; 
      phases[marker_count] = dsg_ph.ch2;
      // fprintf(stderr, "XXX %c %c \n", dsg_ph.ch1, dsg_ph.ch2);
      if(genotypes[marker_count] == MISSING_DATA_CHAR) accession_missing_data_count++;
      marker_count++;
    } // done reading dosages for all markers of this accession
    if(marker_count != markerid_count) {
      fprintf(stderr, "# marker_count, markerid_count: %ld %ld \n", marker_count, markerid_count);
      exit(EXIT_FAILURE); 
    }
    // double fraction_to_keep = 0.75;  
    // if((double)rand()/((double)(RAND_MAX)+1) < fraction_to_keep){
    if(1){
      // if accession does not have too much missing data, construct Accession and store in the_genotypes_set
      if(accession_missing_data_count <= max_acc_missing_data_fraction * the_genotypes_set->marker_ids->size){
	
	Accession* the_accession = construct_accession(acc_id, accession_count, genotypes, phases, accession_missing_data_count);
	for(long jjj=0; jjj<markerid_count; jjj++){ //
	  if(genotypes[jjj] == MISSING_DATA_CHAR){
	    the_genotypes_set->marker_missing_data_counts->a[jjj]++;
	    // accession_missing_data_count++; // not needed; has no effect?
	  }else{
	    the_genotypes_set->marker_alt_allele_counts->a[jjj] += (long)(genotypes[jjj]-48); // 48->+=0, 49->+=1, 50->+=2, etc.
	  }
	}   
	push_to_vaccession(the_genotypes_set->accessions, the_accession);
	accession_count++;
      }else{
	fprintf(stderr, "# Accession: %s rejected due to missing data at %ld out of %ld markers.\n",
		acc_id, accession_missing_data_count, the_genotypes_set->marker_ids->size);
	the_genotypes_set->n_bad_accessions++;
      }
    }
    free(acc_id); // or cut out the middleman (acc_id)?
    free(genotypes);
    }
  } // done reading all lines
  }else{ // threaded!
    threaded_input(g_stream, 720, max_acc_missing_data_fraction, Nthreads, marker_ids, the_genotypes_set);
  }
  double t2 = hi_res_time();
  fprintf(stderr, "t2-t1: %7.5f\n", t2-t1);
  //exit(0);
  fclose(g_stream);
  // fprintf(stderr, "# %ld accessions removed for excessive missing data; %ld accessions kept.\n",
  //	  the_genotypes_set->n_bad_accessions, accession_count);


  
  free(line); // only needs to be freed once.
  the_genotypes_set->n_accessions = the_genotypes_set->accessions->size;
  the_genotypes_set->n_markers = the_genotypes_set->marker_ids->size;
  //fprintf(stderr, "the_genotypes_set->marker_ids->size: %ld   %ld\n", the_genotypes_set->marker_ids->size, the_genotypes_set->n_markers);
  if(DBUG && do_checks) check_genotypesset(the_genotypes_set);

  char buffer[2000];
  sprintf(buffer, "# Number of accessions with genotypes, before filtering: %ld\n", the_genotypes_set->n_accessions + the_genotypes_set->n_bad_accessions);
  append_str_to_vchar(the_genotypes_set->acc_filter_info, buffer);
  sprintf(buffer, "# Removed %ld accessions with missing data fraction > %6.4lf\n", the_genotypes_set->n_bad_accessions, max_acc_missing_data_fraction);
  append_str_to_vchar(the_genotypes_set->acc_filter_info, buffer);
  sprintf(buffer, "# Number of accessions with genotypes, after filtering: %ld\n", the_genotypes_set->n_accessions);
  append_str_to_vchar(the_genotypes_set->acc_filter_info, buffer);
  
}

void threaded_input(FILE* in_stream, long n_lines_in_chunk, double max_acc_md_fraction, long Nthreads, Vstr* marker_ids, GenotypesSet* the_genotypes_set){
  char* line = NULL;
  size_t len = 0;
  long nread = 0;
  long total_lines_read = 0;
  Vstr* accession_lines = construct_vstr(1000);
  Vaccession* all_used_accessions = construct_vaccession(1000);
  while(nread >= 0){ // loop over chunks
    long line_count = 0;
    // read n_markers_in_chunk lines (or up to eof)
    while(
	  (line_count < n_lines_in_chunk) &&
	  ((nread = getline(&line, &len, in_stream)) != -1)
	  ){
      char* line_copy = strcpy( (char*)malloc((nread+1)*sizeof(char)), line);
      chomp(line_copy);
      push_to_vstr(accession_lines, line_copy);
      line_count++;
    }
    total_lines_read += line_count;
    // fprintf(stdout, "line_count: %ld   total lines read: %ld\n", line_count, total_lines_read);

    // ********************************************************
    // *****  Extract genotypes, and quality information  *****
    // *****  Filter if requested and store genotypes     *****
    // ********************************************************
    
    if(Nthreads == 0){ // process without creating any new threads
      threaded_input_struct tis;
      tis.accession_lines = accession_lines; // Vstr*
      tis.first_line = 0;
      tis.last_line = line_count - 1;
      tis.markerid_count = marker_ids->size;
      tis.max_acc_missing_data_fraction = max_acc_md_fraction;
      tis.n_bad_accessions = 0;
      
      // fprintf(stderr, "Nthreads=0; before process_input_lines.\n");
      input_lines_1thread((void*)(&tis));
      //  fprintf(stderr, "Nthreads=0; after process_input_lines.\n");

      for(long im=0; im<marker_ids->size; im++){ // loop over stored markers
	the_genotypes_set->marker_alt_allele_counts->a[im] += tis.marker_alt_allele_counts->a[im];
	the_genotypes_set->marker_missing_data_counts->a[im] += tis.marker_missing_data_counts->a[im];
      }
      for(long ia=0; ia<tis.accessions->size; ia++){
	push_to_vaccession(the_genotypes_set->accessions, tis.accessions->a[ia]);
      }
      	the_genotypes_set->n_bad_accessions += tis.n_bad_accessions;
      free(tis.accessions); // but don't free the c strings containing the actual genotypes (dosages), which are stored in all_used_accessions.
      free(tis.accessions->a);
    }else{ // 1 or more pthreads
      
      threaded_input_struct* tis = (threaded_input_struct*)malloc(Nthreads*sizeof(threaded_input_struct));
      for(long i_thread = 0; i_thread<Nthreads; i_thread++){
	tis[i_thread].accession_lines = accession_lines; // Vstr*
	tis[i_thread].first_line = (i_thread == 0)? 0 : tis[i_thread-1].last_line + 1;;
	tis[i_thread].last_line = (long)((double)(i_thread+1)*(line_count)/Nthreads - 1); //n_lines_in_chunk - 1;
	tis[i_thread].markerid_count = marker_ids->size;
	tis[i_thread].max_acc_missing_data_fraction = max_acc_md_fraction;
	tis[i_thread].n_bad_accessions = 0;
      }
      tis[Nthreads-1].last_line = line_count - 1;

      pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
      for(long i=0; i<Nthreads; i++){ // run the threads
	int iret = pthread_create( thrids+i, NULL, input_lines_1thread, (void*) (tis+i));
	if(iret > 0) fprintf(stderr, "# warning. pthread_create returned non-zero value. Thread %ld \n", (long)thrids[i]);
      }
      for(long i_thread=0; i_thread<Nthreads; i_thread++){ // wait for threads to terminate.
	pthread_join(thrids[i_thread], NULL);
      }
    
      // store results from this chunk
      for(long ith=0; ith<Nthreads; ith++){ // loop over threads
	for(long im=0; im<marker_ids->size; im++){ // loop over stored markers
	  the_genotypes_set->marker_alt_allele_counts->a[im] += tis[ith].marker_alt_allele_counts->a[im];
	  the_genotypes_set->marker_missing_data_counts->a[im] += tis[ith].marker_missing_data_counts->a[im];
	}
	for(long ia=0; ia<tis[ith].accessions->size; ia++){ // loop over markers stored by thread ith
	  push_to_vaccession(the_genotypes_set->accessions, tis[ith].accessions->a[ia]); //chunk_genos[ith]->a[im]);
	}
	the_genotypes_set->n_bad_accessions += tis[ith].n_bad_accessions;
	free(tis[ith].accessions); // but don't free the c strings containing the actual genotypes (dosages), which are stored in all_used_genos.
	free(tis[ith].accessions->a);
      }    
      free(tis);
      free(thrids);
    } // end >=1 pthreads branch

    for(long im=0; im < accession_lines->size; im++){
      free(accession_lines->a[im]); // free the c-strings containing the lines of this chunk.
    }
    accession_lines->size = 0; // done with this chunk, let next chunk overwrite marker_lines->a
  } // end of loop over chunks
}


void* input_lines_1thread(void* x){ // process 1 thread's set of lines
  threaded_input_struct* tis = (threaded_input_struct*)x;
  tis->marker_missing_data_counts = construct_vlong_zeroes(tis->markerid_count);
  tis->marker_alt_allele_counts = construct_vlong_zeroes(tis->markerid_count);
  tis->n_bad_accessions = 0;
  tis->accessions = construct_vaccession(tis->last_line - tis->first_line + 1);
  
  long accession_count = 0;
  for(long i=tis->first_line; i<=tis->last_line; i++){
    char* line = tis->accession_lines->a[i];
    // fprintf(stderr, "line: %s\n", line);
  
    // while((nread = getline(&line, &len, g_stream)) != -1){
    // fprintf(stderr, "# reading accession %ld\n", accession_count);
    char* saveptr = line;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    long marker_count = 0;
    long accession_missing_data_count = 0;
    char* genotypes = (char*)calloc((tis->markerid_count+1), sizeof(char));    
    genotypes[tis->markerid_count] = '\0'; // terminate with null. NEEDED?
    char* phases = (char*)calloc((tis->markerid_count+1), sizeof(char));
    phases[tis->markerid_count] = '\0'; // NEEDED?
    while(1){ // read dosages from one line.   
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL)	break;
      // fprintf(stderr, "token: %s\n", token);
         two_chars dsg_ph = token_to_dosage(token, &(tis->ploidy));
      genotypes[marker_count] = dsg_ph.ch1;
      phases[marker_count] = dsg_ph.ch2;
     
      if(genotypes[marker_count] == MISSING_DATA_CHAR) accession_missing_data_count++;
      marker_count++;
    } // done reading dosages for all markers of this accession
    //  fprintf(stderr, "# line number: %ld ; marker_count, markerid_count: %ld %ld \n", i, marker_count, tis->markerid_count);
    if(marker_count != tis->markerid_count) {
      fprintf(stderr, "# marker_count, markerid_count: %ld %ld \n", marker_count, tis->markerid_count);
      exit(EXIT_FAILURE); 
    }
    // double fraction_to_keep = 0.75;  
    // if((double)rand()/((double)(RAND_MAX)+1) < fraction_to_keep){
    if(1){
      // if accession does not have too much missing data, construct Accession and store in the_genotypes_set
      if(accession_missing_data_count <= tis->max_acc_missing_data_fraction * tis->markerid_count){
	Accession* the_accession = construct_accession(acc_id, accession_count, genotypes, phases, accession_missing_data_count);
	for(long jjj=0; jjj<tis->markerid_count; jjj++){ //
	  if(genotypes[jjj] == MISSING_DATA_CHAR){
	    tis->marker_missing_data_counts->a[jjj]++;
	    accession_missing_data_count++;
	  }else{
	    tis->marker_alt_allele_counts->a[jjj] += (long)(genotypes[jjj]-48); // 48->+=0, 49->+=1, 50->+=2, etc.
	  }
	}   
	push_to_vaccession(tis->accessions, the_accession);
	accession_count++;
      }else{
	fprintf(stderr, "# Accession: %s rejected due to missing data at %ld out of %ld markers.\n",
		acc_id, accession_missing_data_count, marker_count);
	tis->n_bad_accessions++;
      }
    }
    free(acc_id); // or cut out the middleman (acc_id)?
    free(genotypes);
  } // loop over lines
  /* for(long j=0; j<tis->markerid_count; j++){ */
  /*   fprintf(stderr, "%ld  %ld\n", j, tis->marker_missing_data_counts->a[j]); */
  /* }fprintf(stderr, "\n"); */
} // end of  input_lines_1thread

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

two_chars token_to_dosage(char* token, long* ploidy){
  char dsg;
  char phase = NO_PHASE_CHAR; // for 'plus'; corresponds to 1 in pdsgs file, phase = 'm' for -1 in pdsgs file.
  if(strcmp(token, "X") == 0){ // missing data
    dsg = MISSING_DATA_CHAR;
  }else{
    if(token[0] == '-'){ // set phase to 'm', remove the minus, process the rest of token.
      phase = 'm';
      token++; 
    }else if(token[0] == '+'){
      phase = 'p';
      token++;
    }
    long l = str_to_long(token);
    if(errno != 0){
      exit(EXIT_FAILURE);
    }else if(l < 0  ||  l > MAX_PLOIDY){ // if outside the range of possible dosages.
      dsg = MISSING_DATA_CHAR;
    }else{
      if(l > *ploidy) *ploidy = l; // so ploidy ends up as the max of all dosages.
      dsg =  (char)l + 48;
    }
  }
  return (two_chars){dsg, phase};
}

void populate_marker_dosage_counts(GenotypesSet* the_gtsset){
  long nX=0;
  the_gtsset->dosage_counts = construct_vlong_zeroes(4); // (long*)malloc(4*sizeof(long));
  for(long j=0; j<the_gtsset->n_accessions; j++){
    Accession* the_acc = the_gtsset->accessions->a[j];
    for(long i=0; i<the_gtsset->n_markers; i++){
      long dosage = (long)the_acc->genotypes->a[i];
      if(dosage != MISSING_DATA_CHAR){
	the_gtsset->marker_dosage_counts[dosage-48]->a[i]++; // just store the ok dosages (not missing)
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



void set_n_00_1_22_1s(GenotypesSet* the_gtsset){
  for(long i=0; i< the_gtsset->accessions->size; i++){
    Accession* A = the_gtsset->accessions->a[i];
    n_00_1_22_1_accvsall(the_gtsset, A);
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

ND psr(Accession* acc1, Accession* acc2, Vlong* chroms){ // phase mismatch rate
  // considering just markers with both accessions heterozygous, and both having phase info,
  // count how many agree and how many disagree about the phase.

  long n_agree = 0;
  long n_disagree = 0;
  long n_agree_chrom = 0;
  long n_disagree_chrom = 0;
  long n_agree_b = 0;
  long n_disagree_b = 0;
  long n_switches = 0;
  Vchar* p1s = acc1->phases;
  Vchar* p2s = acc2->phases;
  long chrom = chroms->a[0];
  long prev_chrom = -1;
  long prev_switch_position = 0;
  long spacing = -1;
  if(p1s->length != p2s->length  || p1s->length != chroms->size){
    fprintf(stderr, "phases lengths disagree: %ld %ld  chroms->size  %ld\n", p1s->length, p2s->length, chroms->size);
    exit(EXIT_FAILURE);
  }

  bool agree;
  bool prev_agree;
  for(long i=0; i<p1s->length; i++){
    chrom = chroms->a[i];

    if(chrom != prev_chrom){
      fprintf(stderr, "Z na, nda chrom: %ld %ld   %ld %ld \n", prev_chrom, chrom, n_agree_chrom, n_disagree_chrom);
      fprintf(stderr, "BBB: %s  %s   %ld %ld   %ld %ld   %lf\n", acc1->id->a, acc2->id->a, prev_chrom, chrom, n_agree_chrom, n_disagree_chrom, n_disagree_chrom/(double)(n_agree_chrom + n_disagree_chrom));
      if(n_disagree_chrom > n_agree_chrom){
	long tmp = n_disagree_chrom;
	n_disagree_chrom = n_agree_chrom;
	n_agree_chrom = tmp;
      }
      n_agree += n_agree_chrom;
      n_disagree += n_disagree_chrom;
      n_agree_chrom = 0;
      n_disagree_chrom = 0;
      if(prev_chrom > 0){
      spacing = i - prev_switch_position;
      fprintf(stderr, "Aspacing: %ld \n", spacing);
      }
      prev_switch_position = i;
    }
    
    char p1 = p1s->a[i];
    char p2 = p2s->a[i];
    //long chrom = 
    if(p1 == 'p'){
      if(p2 == 'p'){
	agree = true;
	n_agree_chrom++;
	n_agree_b++;
      }else if(p2 == 'm'){
	agree = false;
	n_disagree_chrom++;
	n_disagree_b++;
      }
    }else if(p1 == 'm'){
      if(p2 == 'm'){
	agree = true;
	n_agree_chrom++;
	n_agree_b++;
      }else if(p2 == 'p'){
	agree = false;
	n_disagree_chrom++;
	n_disagree_b++;
      }
    }
    
    if((chrom == prev_chrom)  &&  (agree != prev_agree)){
      n_switches++;
      spacing = i - prev_switch_position;
      fprintf(stderr, "Bspacing: %ld \n", spacing);
      prev_switch_position = i;
    }
    prev_agree = agree;
    prev_chrom = chrom;
  } // end loop over markers
  spacing = p1s->length - prev_switch_position;
  fprintf(stderr, "Cspacing: %ld \n", spacing);
  if(n_disagree_chrom > n_agree_chrom){
    long tmp = n_disagree_chrom;
    n_disagree_chrom = n_agree_chrom;
    n_agree_chrom = tmp;
  }
  n_agree += n_agree_chrom;
  n_disagree += n_disagree_chrom;

  fprintf(stderr, "Z na, nda chrom: %ld %ld \n", n_agree_chrom, n_disagree_chrom);
  fprintf(stderr, "ZZZ: %ld %ld  %ld %ld  %ld  %7.5f\n", n_agree, n_disagree, n_agree_b, n_disagree_b, n_switches, n_over_d((ND){n_switches, n_agree+n_disagree-1}));
																				   
  ND result = {n_switches, n_agree + n_disagree - 1};
  return result;
}

void print_genotypesset_stats(GenotypesSet* gtss){
  fprintf(stderr, "max_marker_missing_data_fraction: %8.4f\n", gtss->max_marker_missing_data_fraction);
  fprintf(stderr, "min_minor_allele_frequency: %8.4f\n", gtss->min_minor_allele_frequency);
  fprintf(stderr, "n_accessions: %ld  %ld\n", gtss->n_accessions, gtss->accessions->size);
  fprintf(stderr, "n bad, ref accessions: %ld %ld \n", gtss->n_bad_accessions, gtss->n_ref_accessions);
  fprintf(stderr, "n_markers: %ld %ld \n", gtss->n_markers, gtss->marker_ids->size);
  fprintf(stderr, "sizes, marker_md_counts %ld,  marker_alt_allele_counts %ld\n",
	  gtss->marker_missing_data_counts->size, gtss->marker_alt_allele_counts->size);
  if(gtss->mafs != NULL) fprintf(stderr, "sizes, mafs %ld \n", gtss->mafs->size);
									 //, marker_dosage_counts: %ld\n", gtss->mafs->size, gtss->marker_dosage_counts);
  if(gtss->dosage_counts != NULL) fprintf(stderr, "size of dosage_counts %ld \n", gtss->dosage_counts->size);
  fprintf(stderr, "agmr0: %8.5f\n", gtss->agmr0);
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
    //fprintf(stderr, "j, %ld;  mdcounts: %ld %ld \n", j, marker_md_counts[j], gtss->marker_missing_data_counts->a[j]);
    assert(marker_md_counts[j] == gtss->marker_missing_data_counts->a[j]);
    assert(marker_alt_allele_counts[j] == gtss->marker_alt_allele_counts->a[j]);
  }
  free(marker_md_counts);
  // fprintf(stderr, "# Successfully completed check_genotypesset\n");
}

void filter_genotypesset(GenotypesSet* the_gtsset, FILE* ostream){ // construct a new set of 'filtered' accession genotypes, which replace the raw ones.
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
      double min_min_allele_count = fmax(max_possible_alt_alleles * the_gtsset->min_minor_allele_frequency, MIN_ALT_ALLELES);
      double max_min_allele_count = fmin(max_possible_alt_alleles * (1.0 - the_gtsset->min_minor_allele_frequency), max_possible_alt_alleles-MIN_ALT_ALLELES);
      bool alt_allele_freq_ok = (alt_allele_counts->a[i] >= min_min_allele_count)  && // alternative allele frequency not too small,
	(alt_allele_counts->a[i] <= max_min_allele_count);     // and not too large
      // fprintf(stderr, "# alt allele info: %ld  %ld  %f  %f \n", max_possible_alt_alleles, alt_allele_counts->a[i], min_min_allele_count, max_min_allele_count);
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
	// fprintf(stderr, "min/max allele_counts: %f %f  alt_allele_count: %ld\n", min_min_allele_count, max_min_allele_count, alt_allele_counts->a[i]);
      }
    }else{
      too_much_missing_data_count++;
    }
  } // end of loop over markers

  double raw_md_fraction = (double)mdsum_all/(double)(marker_md_counts->size*n_accs);
  double filtered_md_fraction = (double)mdsum_kept/(double)(n_markers_to_keep*n_accs);
  double raw_minor_allele_freq = (double)altallelesum_all/(marker_md_counts->size*n_accs*the_gtsset->ploidy);
  double filtered_minor_allele_freq = (double)altallelesum_kept/(double)(n_markers_to_keep*n_accs*the_gtsset->ploidy);

  // store information on results of filtering
  char buffer[2000];
  sprintf(buffer, "# There are %ld markers before filtering, missing data fraction = %5.3lf, minor allele frequency = %5.3f\n",
	  marker_md_counts->size, raw_md_fraction, raw_minor_allele_freq);
  append_str_to_vchar(the_gtsset->marker_filter_info, buffer);
  sprintf(buffer, "# Removed %ld markers with missing data fraction > %6.4f\n", too_much_missing_data_count, the_gtsset->max_marker_missing_data_fraction);
  append_str_to_vchar(the_gtsset->marker_filter_info, buffer);
  sprintf(buffer, "# Removed an additional %ld markers with maf < %6.4f\n", maf_too_low_count, the_gtsset->min_minor_allele_frequency);
  append_str_to_vchar(the_gtsset->marker_filter_info, buffer);
  sprintf(buffer, "# Filtered data has %ld markers, missing data fraction = %6.4lf, minor allele frequency = %5.3lf\n",
	  n_markers_to_keep, filtered_md_fraction, filtered_minor_allele_freq);
  append_str_to_vchar(the_gtsset->marker_filter_info, buffer);

  // ***********************************
  // *****  filter the genotypes  ******
  // ***********************************
  Vaccession* the_accessions = construct_vaccession(the_gtsset->n_accessions); //(Accession*)malloc(the_gtsset->n_accessions*sizeof(Accession)); 
  for(long i=0; i<the_gtsset->n_accessions; i++){ // loop over accessions
    char* raw_gts = the_gtsset->accessions->a[i]->genotypes->a; // the string with all the genotypes for accession i
    char* filtered_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    char* raw_phases = the_gtsset->accessions->a[i]->phases->a;
    //fprintf(stderr, "raw gts, phases lengths: %ld %ld    ", the_gtsset->accessions->a[i]->genotypes->length, the_gtsset->accessions->a[i]->phases->length);
    char* filtered_phases = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    long acc_md_count = 0;
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	filtered_gts[k] = raw_gts[j];
	filtered_phases[k] = raw_phases[j]; // if data is unphased, these will be meaningless and unused.
	// if(i==0  &&  the_gtsset->phased) filtered_chromosomes->a[k] = raw_chromosomes->a[j];
	k++;
	if(raw_gts[j] == MISSING_DATA_CHAR) acc_md_count++;
      }
    }
    filtered_gts[k] = '\0'; // terminate with null.
    filtered_phases[k] = '\0';
    if(DBUG && do_checks) assert(k == n_markers_to_keep);
    Accession* the_accession = construct_accession( the_gtsset->accessions->a[i]->id->a, i, filtered_gts, filtered_phases, acc_md_count ); //(Accession*)malloc(sizeof(Accession));
    free(filtered_gts); // this str was copied into newly allocated memory in the_accession->
    free(filtered_phases);
    push_to_vaccession(the_accessions, the_accession);
  }
  // *****  end of filtering of genotypes  ******
  
  // *******************************************************
  // *****  if phased, filter the chromosomes array.  ******
  // *******************************************************
  //if(the_gtsset->phased){ 
  Vlong* raw_chromosomes = the_gtsset->chromosomes;
  Vlong* filtered_chromosomes = construct_vlong(n_markers_to_keep+1);
  long k=0;
  for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
    if(md_ok->a[j] == 1){
      //	filtered_gts[k] = raw_gts[j];
      //	filtered_phases[k] = raw_phases[j];
      filtered_chromosomes->a[k] = raw_chromosomes->a[j];
      k++;
    }
  }
    
  filtered_chromosomes->size = n_markers_to_keep;
  free_vlong(the_gtsset->chromosomes);
  the_gtsset->chromosomes = filtered_chromosomes;
  fprintf(stderr, "in filter_genotypesset. filtered chroms size: %ld \n", the_gtsset->chromosomes->size);
  
  // *****  end of filtering of chromosomes array (phased data only)

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
    if(2*the_gtsset->marker_alt_allele_counts->a[i] > ploidy*(n_accessions - the_gtsset->marker_missing_data_counts->a[i])){ // for this marker switch 0 and 2, and 1, -1
      n_markers_rectified++;
      for(long j=0; j<n_accessions; j++){
	char dosage = the_gtsset->accessions->a[j]->genotypes->a[i];

	if(dosage != MISSING_DATA_CHAR){
	  if(0){
	    long ldosage = dosage-48;
	    if(2*ldosage != ploidy){ 
	      the_gtsset->accessions->a[j]->genotypes->a[i] = (char)(ploidy - ldosage) + 48;
	      delta_dosage += ploidy - 2*ldosage;
	    }
	  }else{ // version for ploidy == 2 only
	    if(dosage == '0'){
	      dosage = '2';
	      delta_dosage += 2;
	    }else if(dosage == '2'){
	      dosage = '0';
	      delta_dosage -= 2;
	    }
	    the_gtsset->accessions->a[j]->genotypes->a[i] = dosage;
	  }
	
	  // reverse phases of heterozygs p->m, and m->p.  homozygs unaffected
	  char phase = the_gtsset->accessions->a[j]->phases->a[i];
	  if(phase == 'p'){
	    phase = 'm';
	  }else if(phase == 'm'){
	    phase = 'p';
	  }
	  the_gtsset->accessions->a[j]->phases->a[i] = phase;
	}
      }
    }
    // so far genotypes of each accession have (if necessary) been 'rectified', but now marker_alt_allele_counts needs to be updated
    //  fprintf(stderr, "%ld  %ld   ", the_gtsset->marker_alt_allele_counts->a[i], delta_dosage);
    the_gtsset->marker_alt_allele_counts->a[i] += delta_dosage;
  }
  fprintf(stderr, "##### n_markers_rectified: %ld\n", n_markers_rectified);
}

void set_Abits_Bbits(GenotypesSet* the_genotypesset, long Nthreads){ // diploid only
  if(Nthreads < 0){ // unthreaded
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
	  if(gt == '0'){ // 00 i.e. A=0, B=0
	    // leave the bit as 0 in both A and B.
	  }else if(gt == '1'){ // 01  A=0, B=1
	    B |= bits[j]; // set jth bit of B
	 
	  }else if(gt == '2'){ // 11  A=1, B=1
	    A |= bits[j];
	    B |= bits[j];
	  }else{ // missing data; 10  A=1, B=0
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
  }else if(Nthreads == 0){
    threaded_setAB_struct tsAB;
    tsAB.gtss = the_genotypesset;
    tsAB.first = 0;
    tsAB.last = the_genotypesset->accessions->size - 1;
    //  fprintf(stderr, "XXXXXXXX calling set_Abits_Bbits_1thread...\n"); //getchar();
    set_Abits_Bbits_1thread((void*) &tsAB);
  }else{ // Nthreads >= 1
    threaded_setAB_struct* tsAB = (threaded_setAB_struct*)malloc(Nthreads*sizeof(threaded_setAB_struct));
    pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
    for(long ith=0; ith < Nthreads; ith++){
      tsAB[ith].gtss = the_genotypesset;
      tsAB[ith].first = (ith == 0)? 0 : tsAB[ith-1].last + 1;    
      tsAB[ith].last = (ith == Nthreads-1)? the_genotypesset->accessions->size - 1 : tsAB[ith].first + (long)(the_genotypesset->accessions->size/Nthreads) - 1;
    
      
      int iret = pthread_create( thrids+ith, NULL, set_Abits_Bbits_1thread, (void*) (tsAB+ith));
      if(iret > 0) fprintf(stderr, "# warning. pthread_create returned non-zero value. Thread %ld \n", (long)thrids[ith]);
    }
    for(long i_thread=0; i_thread<Nthreads; i_thread++){ // wait for threads to terminate.
      pthread_join(thrids[i_thread], NULL);
    }
  }
}

void* set_Abits_Bbits_1thread(void* x){
  threaded_setAB_struct* tsAB =  (threaded_setAB_struct*)x;
   unsigned long long xx = 1, bits[64];
  for(long j=0; j<64; j++, xx = xx<<1){
    bits[j] = xx; // bits[i] is an unsigned long long with the ith bit set (i.e. 1, with all others 0)
  }

  for(long i_acc = tsAB->first; i_acc <= tsAB->last; i_acc++){
    Accession* the_acc = tsAB->gtss->accessions->a[i_acc];
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
} // end of one_thread_set_Abits_Bbits


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
	push_to_vlong(acc->alt_homozygs, j); // store index of alt homozyg marker	
      }else if(dosage == 0){
	push_to_vlong(acc->ref_homozygs, j);
      }
      
    }
  }
}

void set_chromosome_start_indices(GenotypesSet* the_gtsset){
  Vlong* chrom_start_indices = construct_vlong(50);
  long prev_chrom_number = -1; // the_gtsset->chromosomes->a[0];
  //fprintf(stderr, "chroms size: %ld \n", the_gtsset->chromosomes->size);
  for(long i=0; i<the_gtsset->chromosomes->size; i++){
    long chrom_number = the_gtsset->chromosomes->a[i];
    //fprintf(stderr, "XXX: %ld %ld \n", prev_chrom_number, chrom_number);
    if(chrom_number != prev_chrom_number) push_to_vlong(chrom_start_indices, i); 
    prev_chrom_number = chrom_number;
  }
  push_to_vlong(chrom_start_indices, the_gtsset->chromosomes->size); 
  the_gtsset->chromosome_start_indices = chrom_start_indices;
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

ND bitwise_R(Accession* parent, Accession* offspring){
  unsigned long long i0, i1, i2, j0, j1, j2;
  long N = 0;
  long D = 0;
  for(long k=0; k<parent->Abits->size; k++){ // loop over 64 bit chunks
    unsigned long long iA = parent->Abits->a[k];
    unsigned long long iB = parent->Bbits->a[k];
    unsigned long long jA = offspring->Abits->a[k];
    unsigned long long jB = offspring->Bbits->a[k];
    i0 = ~(iA | iB);
    // i1 = ~iA & iB;
    i2 = iA & iB;
    j0 = ~(jA | jB);
    j1 = ~jA & jB;
    j2 = jA & jB;
    unsigned long long is01_21 = ~(iA ^ iB) & j1;
    unsigned long long is00_22 = (i0 & j0) | (i2 & j2); // ~( (iA ^ iB) | (jA ^ jB) ); // = ~(iA ^ iB) & ~(jA ^ jB) 
    N += __builtin_popcountll(is01_21);
    D += __builtin_popcountll(is00_22);
  }
  D += N;
  return (ND){N, D};
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
    // D: difference, S: same ; o: homozyg, x: (at least 1) heterozyg.
    //   Ndo  i.e. number of markers with both homozyg, but different
    //     (i.e. 02 and 20)
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

ND bitwise_hgmr(Accession* acc1, Accession* acc2){
  ND rval = {0, 0};
  unsigned long long isOi, isXi, isOj, isXj;
  for(long i_long = 0; i_long < acc1->Abits->size; i_long++){
    unsigned long long iA = acc1->Abits->a[i_long];
    unsigned long long iB = acc1->Bbits->a[i_long];
    unsigned long long jA = acc2->Abits->a[i_long];
    unsigned long long jB = acc2->Bbits->a[i_long];
    isOi = ~(iA ^ iB);
    isOj = ~(jA ^ jB);
    //isXi = ~iA & iB;
    //isXj = ~jA & jB;
    unsigned long long isDo = (iA ^ jB) & isOi & isOj;
    unsigned long long isSo = ~(iA ^ jB) & isOi & isOj;
    //unsigned long long isDx = (isOi & isXj) | (isOj & isXi);
    //unsigned long long isSx = isXi & isXj;
    // D: difference, S: same, o: homozyg, x: (at least 1) heterozyg.
    //   Ndo  i.e. number of markers with both homozyg, but different
    //     (i.e. 02 and 20)
      rval.n += __builtin_popcountll(isDo); 
      //   Nso   00, 22
      rval.d += __builtin_popcountll(isSo);
  }
  rval.d += rval.n; 
  return rval;
}

void calculate_hgmrs(GenotypesSet* the_genotypes_set, Viaxh** pairwise_info, double max_hgmr){
  long n_hgmrs_calculated = 0;
  long n_hgmrs_le_max = 0;
  long n_acc = the_genotypes_set->accessions->size;
  for(long ii=0; ii<n_acc; ii++){
    if(ii % 100  == 0) fprintf(stderr, "# ii: %ld\n", ii);
    Accession* A1 = the_genotypes_set->accessions->a[ii];
    for(long jj=ii+1; jj<n_acc; jj++){
      Accession* A2 = the_genotypes_set->accessions->a[jj];
	   
      ND hgmr_nd = bitwise_hgmr(A1, A2);
      // fprintf(stderr, "hgmr: %7.5f  %7.5f\n", n_over_d(hgmr_nd), the_genotypes_set->mean_hgmr);
      double hgmr = n_over_d(hgmr_nd)/the_genotypes_set->mean_hgmr;	      
      n_hgmrs_calculated++;
	    
      if(hgmr <= max_hgmr){	      
	n_hgmrs_le_max++;
	push_to_viaxh(pairwise_info[ii], jj, hgmr);
	push_to_viaxh(pairwise_info[jj], ii, hgmr);
      }
    }
  }
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

void set_n2exp0s(GenotypesSet* gtset, long i){
  Accession* A = gtset->accessions->a[i];
  Vlong* Ad2s = A->alt_homozygs;
  long n0s = 0;
 for(long i=0; i<Ad2s->size; i++){ // consider markers with dosage==2 in accession 1
    long idx = Ad2s->a[i];
    n0s += gtset->marker_dosage_counts[0]->a[idx]; // n0s_this_marker;
  }
  A->n2exp0s = n0s;
}

ND xhgmr(GenotypesSet* gtset, Accession* a1, Accession* a2, bool quick){
  // if(quick) count markers with  dosage 2 in a1, dosage 0 in a2 (for numerator)
  // denominator is number expected of a2 dosage = 0 for those markers,
  // just based on fraction of accessions having 0 for each of those markers.
  // if(!quick) also do in the other direction (2's in a2, 0's in a1).
  Vlong* a1d2s = a1->alt_homozygs;
  double expected_refds = 0;
  long counted_refds = 0;
  long n0s = a1->n2exp0s; // 0
  
  for(long i=0; i<a1d2s->size; i++){ // consider markers with dosage==2 in accession 1
    long idx = a1d2s->a[i];
    char a2_dosage = a2->genotypes->a[idx];
    if(a2_dosage == '0') counted_refds++;
    //n0s += gtset->marker_dosage_counts[0]->a[idx]; // n0s_this_marker;
  }
  /* if(n0s != n0sx){ */
  /*   fprintf(stderr, "AAA: %ld  %ld \n", n0s, n0sx); */
  /*   exit(0); */
  /* } */
  expected_refds = (double)n0s/(double)gtset->accessions->size;
  if(quick && ( counted_refds > 100 )  && ( counted_refds/expected_refds > 0.5 )) return (ND){counted_refds, expected_refds};

  n0s += a2->n2exp0s;
  Vlong* a2d2s = a2->alt_homozygs;
  for(long i=0; i<a2d2s->size; i++){ // consider markers with dosage==2 in accession 2
    long idx = a2d2s->a[i];
    char a1_dosage = a1->genotypes->a[idx];
    if(a1_dosage == '0') counted_refds++;
    //n0s += gtset->marker_dosage_counts[0]->a[idx]; // n0s_this_marker;
  }
  /* if(n0s != n0sx){ */
  /*   fprintf(stderr, "ABC: %ld  %ld \n", n0s, n0sx); */
  /*   exit(0); */
  /* } */
  expected_refds = (double)n0s/(double)gtset->accessions->size;
  return (ND){counted_refds, expected_refds};
}

void calculate_xhgmrs(GenotypesSet* the_genotypes_set, Viaxh** pairwise_info, bool quick_xhgmr, double max_xhgmr){
  // calculate xhgmr for all pairs of accessions in the_gtset
  long n_xhgmrs_calculated = 0;
  long n_xhgmrs_le_max = 0;
  long n_acc = the_genotypes_set->accessions->size;
  for(long ii=0; ii<n_acc; ii++){
    if(ii % 100  == 0) fprintf(stdout, "# ii: %ld\n", ii);
    Accession* A1 = the_genotypes_set->accessions->a[ii];
    for(long jj=ii+1; jj<n_acc; jj++){
      Accession* A2 = the_genotypes_set->accessions->a[jj];
      ND the_xhgmr = xhgmr(the_genotypes_set, A1, A2, quick_xhgmr);
      n_xhgmrs_calculated++;
      double dbl_xhgmr = (the_xhgmr.d > 0)? (double)the_xhgmr.n/the_xhgmr.d : NAN;
      if(dbl_xhgmr <= max_xhgmr){
	n_xhgmrs_le_max++;
	push_to_viaxh(pairwise_info[ii], jj, dbl_xhgmr);
	push_to_viaxh(pairwise_info[jj], ii, dbl_xhgmr);
      }
    }
  }
  fprintf(stdout, "# n xhgmrs calculated: %ld ;  <= %8.5f :  %ld\n", n_xhgmrs_calculated, max_xhgmr, n_xhgmrs_le_max);
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


void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){ 
  fprintf(fh, "# max_marker_missing_data_fraction: %8.4lf \n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "# min_minor_allele_frequency: %8.4lf \n", the_gtsset->min_minor_allele_frequency);
  fprintf(fh, "MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]);
  }fprintf(fh, "\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    Accession* acc = the_gtsset->accessions->a[i];
    fprintf(fh, "%24s  ", acc->id->a);
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

Vidxid* construct_vidxid(const Vaccession* accessions){
  Vidxid* the_vidxid = (Vidxid*)malloc(sizeof(Vidxid));
  the_vidxid->capacity = accessions->size;
  the_vidxid->size = the_vidxid->capacity;
  the_vidxid->a = (IndexId**)malloc(the_vidxid->size*sizeof(IndexId*));
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId));
    the_idxid->index = i;
    char* the_id =  accessions->a[i]->id->a;
    the_idxid->id = strcpy((char*)malloc((strlen(the_id)+1)*sizeof(char)), the_id);
    the_vidxid->a[i] = the_idxid;
  }
  return the_vidxid;  
}

Vidxid* construct_sorted_vidxid(const Vaccession* accessions){
  Vidxid* the_vidxid = construct_vidxid(accessions);
  sort_vidxid_by_id(the_vidxid);
  // for(long i=0; i<the_vidxid->size; i++){ IndexId* x = the_vidxid->a[i];  }
  //  if(DBUG) assert(check_idxid_map(the_vidxid, accessions) == 1);
  return the_vidxid;
}

long check_idxid_map(Vidxid* vidxid, const Vaccession* accessions){
  for(long i=0; i<accessions->size; i++){
    char* id = accessions->a[i]->id->a;
    long idx = index_of_id_in_vidxid(vidxid, id);
    // fprintf(stderr, "%ld %ld %s\n", i, idx, id);
    if(idx != i) return 0;
  }
  return 1;
}



//  #####  unused  ############################
 
double agmr0(GenotypesSet* the_gtsset){  
  double n_exp_different = 0;
  double n_exp_ok = 0; 
  
  for(long i=0; i<the_gtsset->marker_ids->size; i++){
    double ok_count = (double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
    double f0 = the_gtsset->marker_dosage_counts[0]->a[i]/ok_count;
    double f1 = the_gtsset->marker_dosage_counts[1]->a[i]/ok_count;
    double f2 = the_gtsset->marker_dosage_counts[2]->a[i]/ok_count;;
    n_exp_different += 2*(f0*(f1+f2) + f1*f2);
    n_exp_ok += ok_count/(double)the_gtsset->accessions->size;
  }
 
  double agmr0 = (n_exp_ok > 0)? n_exp_different/n_exp_ok : -1;
  the_gtsset->agmr0 = agmr0;
  return agmr0;
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

void  n_00_1_22_1_accvsall(const GenotypesSet* the_gtsset, Accession* A ){
  double n_exp_00_1_22_1 = 0;
  double n_exp_00_1_22_1_self = 0;
   
  for(long i=0; i<the_gtsset->marker_ids->size; i++){
    if(A->genotypes->a[i] - 48  ==  1){ // sum over markers with A having dosage = 1
      double ok_count = (double)(the_gtsset->accessions->size - the_gtsset->marker_missing_data_counts->a[i]);
      double f0 = the_gtsset->marker_dosage_counts[0]->a[i]/ok_count;
      double f2 = the_gtsset->marker_dosage_counts[2]->a[i]/ok_count;
      n_exp_00_1_22_1 += f0*f0 + f2*f2;
      n_exp_00_1_22_1_self += f0 + f2;
	}
  }
  A->n_exp_00_1_22_1 = n_exp_00_1_22_1;
  A->n_exp_00_1_22_1_self = n_exp_00_1_22_1_self;
  // fprintf(stderr, "%8.5f  %8.5f \n", n_exp_00_1_22_1, n_exp_00_1_22_1_self);
  // return (two_longs) {n_exp_00_1_22_1, n_exp_00_1_22_1_self};
}

/*
void set_agmr0s(GenotypesSet* the_gtsset){
  for(long i=0; i< the_gtsset->accessions->size; i++){
    Accession* A = the_gtsset->accessions->a[i];
    A->agmr0 = agmr0_accvsall(the_gtsset, A);
  }
  }/* */

/*
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
}/* */

/*
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
} /* */

void print_info_re_filtering(GenotypesSet* the_gtsset, FILE* fh){
  /* fprintf(fh, "# Number of accessions with genotypes, before filtering: %ld\n", the_gtsset->n_accessions + the_gtsset->n_bad_accessions); */
  /* fprintf(fh, "# Removing accessions with missing data fraction > %6.4lf\n"); */
  /* fprintf(fh, "# There are %ld markers before filtering, missing data fraction = %5.3lf, minor allele frequency = %5.3f\n", */
  /* 	  the_gtsset->marker_md_counts->size, the_gtsset->raw_md_fraction, the_gtsset->raw_minor_allele_freq);    */
  /* fprintf(fh, "# Removing markers with missing data fraction > %6.4lf\n", the_gtsset->max_marker_md_fraction); */
  /* fprintf(fh, "#   removed %ld markers for missing data fraction > %f5.2\n", the_gtsset->too_much_missing_data_count, the_gtsset->max_marker_missing_data_fraction); */
  /* fprintf(fh, "# Removing markers with minor allele frequency < %5.3f\n", the_gtsset->min_minor_allele_frequency); */
  /* fprintf(fh, "#   removed an additional %ld markers for maf < %5.2f\n", the_gtsset->maf_too_low_count, the_gtsset->min_minor_allele_frequency);  */
  /* fprintf(fh, "# Filtered data has %ld markers, missing data fraction = %5.3lf, minor allele frequency = %5.3lf\n", */
  /* 	  the_gtsset->n_markers_to_keep, the_gtsset->filtered_md_fraction, the_gtsset->filtered_minor_allele_freq); */

}

/*
two_doubles logPABPBA(GenotypesSet* the_gtsset, Accession* A, Accession* B){ // generalized hgmr (ploidy can be > 2; reduces to hgmr for diploid)
  if(A == NULL  ||  B == NULL){
    two_doubles result = {0, 0};
    return result;
  }
  long forbidden_count = 0;
  long denom = 0;
  long ploidy = the_gtsset->ploidy;
  double logPAB = 0;
  double logPBA = 0;
   Vlong* mdcs0 = the_gtsset->marker_dosage_counts[0];
    Vlong* mdcs1 = the_gtsset->marker_dosage_counts[1];
     Vlong* mdcs2 = the_gtsset->marker_dosage_counts[2];
  for(long i=0; i<the_gtsset->n_markers; i++){
   
    long a_dosage = A->genotypes->a[i];
    long b_dosage = B->genotypes->a[i];
    if(a_dosage == MISSING_DATA_CHAR  ||  b_dosage == MISSING_DATA_CHAR) continue;
    a_dosage -= 48; b_dosage -= 48;
    
    long delta_dosage = labs(a_dosage - b_dosage);
    if(2*delta_dosage > ploidy){
      forbidden_count++;
    }
    else if(2*a_dosage != ploidy  &&  2*b_dosage != ploidy){
      denom++;
    }
    long n0 = mdcs0->a[i];
    long n1 = mdcs1->a[i];
    long n2 = mdcs2->a[i];
    long ntot = n0 + n1 + n2;
    double f0 = (double)n0/ntot;
    double f1 = (double)n1/ntot;
    double f2 = (double)n2/ntot;
    if(a_dosage == 0){
      if(b_dosage == 1){
	double pab = 0.5*f1 + f2; // a parent of b
	double pba = 0.5*(f0 + 0.5*f1);
	logPAB += log(pab);
	logPBA += log(pba);
      }
    }else if(a_dosage == 1){
      if(b_dosage == 0){
	double pab = 0.5*(f0 + 0.5*f1);
	  double pba = 0.5*f1 + f2; // a parent of b
	logPAB += log(pab);
	logPBA += log(pba);
      }else if(b_dosage == 2){
	double pab = 0.5*(0.5*f1 + f2);
	double pba = f0 + 0.5*f1;
	logPAB += log(pab);
	logPBA += log(pba);
      }
    }else if(a_dosage == 2){
      if(b_dosage == 1){
	double pab = 0.5*f1 + f0;
	double pba = 0.5*(f2 + 0.5*f1);
	logPAB += log(pab);
	logPBA += log(pba);
      }
    }
  }
  denom += forbidden_count;
  
  two_doubles result = {logPAB, logPBA};
  return result;
} /* */




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

/*
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
} /* */
/*
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
} /* */

/*
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
}   /* */
