#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <pthread.h>
//#include <assert.h>
//#include <stdbool.h>
#include "vect.h"

// #include "gtset.h"

typedef struct{
  long N; // allocated size
  long first;
  long last;
  Vdouble** vds;
  Vdouble* dists;
} Nvds; // an array of N Vdouble*  

double d(Vdouble* x1, Vdouble* x2){
  double d = 0;
  for(long i=0; i< x1->size; i++){
    d += fabs(x1->a[i] - x2->a[2]);
  }
}

double dist(Vdouble* v1, Vdouble* v2);
//Vdouble* all_distances(long N, Vdouble** nvs);
void* all_distances(void* x); 

int
main(int argc, char *argv[])
{
  clockid_t clkid = CLOCK_MONOTONIC;
  struct timespec tspec;
  if(clock_getres(clkid, &tspec) == 0){
    fprintf(stderr, "# clock resolution: %ld %ld \n", (long)tspec.tv_sec, (long)tspec.tv_nsec);
  }else{
    exit(1);
  }

  long Nthreads = 6;
  
  long Nelem = 10000; // number of doubles in each vector
  long Nvectors = 1000; // number of vectors
  unsigned rand_seed = 1234567; // (unsigned)time(0);
  srand(rand_seed);

  fprintf(stdout, "# rng seed: %ld\n", (long)rand_seed);

  // *****  initialize  ***** //
  Vdouble** number_vectors = (Vdouble**)malloc(Nvectors*sizeof(Vdouble*));
  //  fprintf(stderr, "# Nelem %ld \n", Nelem);
  for(long i=0; i<Nvectors; i++){
    //   fprintf(stderr, "#i: %ld \n", i);
    Vdouble* nv = construct_vdouble(Nelem);
    for(long j=0; j<Nelem; j++){
      double x = (double)rand()/RAND_MAX;
      //    fprintf(stderr, "# %ld %ld  %8.5f\n", i, j, x);
      //     fprintf(stderr, "# cap, size: %ld %ld \n", nv->capacity, nv->size);
      push_to_vdouble(nv, x);
      
    }
    number_vectors[i] = nv;
  }

  
 
 
  
  if(0){
    struct timespec tsp1;
    clock_gettime(clkid, &tsp1);
    double starttime = tsp1.tv_sec + 1.0e-9*tsp1.tv_nsec;

     Nvds nvds1;
  nvds1.N = Nvectors;
  nvds1.vds = number_vectors;
  nvds1.dists = construct_vdouble(1000);
    nvds1.first = 0;
    nvds1.last = Nvectors-1;
    all_distances((void*) &nvds1);
   
     struct timespec tsp2;
    clock_gettime(clkid, &tsp2);
    double endtime = tsp2.tv_sec + 1.0e-9*tsp2.tv_nsec;
    fprintf(stdout, "# time for dist calcs %12.6lf seconds.\n", endtime-starttime);
    
    for(long i=0; i<nvds1.dists->size; i++){
      fprintf(stdout, "d1: %8.3f\n", nvds1.dists->a[i]);
    }
  }else{
    int iret = 0;
    Nvds* vdses = (Nvds*)malloc(Nthreads*sizeof(Nvds));
    pthread_t* thrids = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
    struct timespec tsp1;
    clock_gettime(clkid, &tsp1);
    double starttime = tsp1.tv_sec + 1.0e-9*tsp1.tv_nsec;
    
    for(long i=0; i<Nthreads; i++){
   
      vdses[i].N = Nvectors;
      vdses[i].vds = number_vectors;
      vdses[i].dists = construct_vdouble(1000);
      if(i == 0){
      vdses[i].first = 0;
      vdses[i].last = (long)(Nvectors * (1.0 - sqrt((Nthreads-1.0)/Nthreads)));
      }else{
      vdses[i].first = vdses[i-1].last + 1;
      vdses[i].last = (long)(Nvectors * (1.0 - sqrt((double)(Nthreads-i-1)/Nthreads)));
      }
      if(i == (Nthreads-1)) vdses[i].last = Nvectors-1;

      iret += pthread_create( thrids+i, NULL, all_distances, (void*) (vdses+i)); 
    }

    for(long i=0; i<Nthreads; i++){
    pthread_join(thrids[i], NULL);
    }
    

    struct timespec tsp2;
    clock_gettime(clkid, &tsp2);
    double endtime = tsp2.tv_sec + 1.0e-9*tsp2.tv_nsec;
    fprintf(stdout, "# time for dist calcs %12.6lf\n", endtime-starttime);
    
    // fprintf(stdout, "# nvds1.dists->size: %ld \n", nvds1.dists->size);
    // fprintf(stdout, "# nvds2.dists->size: %ld \n", nvds2.dists->size);
    for(long ii=0; ii<Nthreads; ii++){
      fprintf(stderr, "# N vectors analyzed by thread %ld, is: %ld\n", ii, vdses[ii].dists->size);
    for(long i=0; i<vdses[ii].dists->size; i++){
      fprintf(stdout, "d: %8.3f thread: %ld\n", vdses[ii].dists->a[i], ii);
    }   
    }
  }


 
}
// end of main

double dist(Vdouble* v1, Vdouble* v2){
  double d = 0;
  if(v1->size != v2->size) exit(1);
  for(long i=0; i<v1->size; i++){
    d += pow(v1->a[i] - v2->a[i], 2);
  }
  return sqrt(d);
}

void* all_distances(void* x){
  Nvds* nvds = (Nvds*)x;
  long N = nvds->N;
  long f = nvds->first;
  long l = nvds->last;
  Vdouble** vds = nvds->vds;
  Vdouble* dists = construct_vdouble(10);
  for(long i=f; i<=l; i++){
    Vdouble* v1 = vds[i];
    for(long j=i; j<N; j++){
      Vdouble* v2 = vds[j];
      double d = dist(v1, v2);
      push_to_vdouble(nvds->dists, d);
    }
  }
}




  
    
