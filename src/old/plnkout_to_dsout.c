#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#define INIT_ARRAY_SIZE 10000  // initial size of array for storing ids.

// read in two files, one with N accession ids, and the other with an NxN matrix of distances
// and output in format similar to duplicate_search output, i.e. on each line
// id1 id2 distance
// only distances <= max_distance are output. 

// usage:  plinkout_to_dsout -i xxx.mdist.id -d xxx.mdist -m 0.16  -o xxx.dists

char *id_file = NULL;
char *distmatrix_file = NULL;
char *output_file = NULL;
FILE *idstream;
FILE *distmatrixstream;
FILE *outputstream;
double max_distance = 0.2; // default max distance

typedef struct // Define a struct to store id1 and id2. a struct is a collection of variables of different types under a single name.
{
  char id1[256]; // id1 and id2 are strings of length 256
  char id2[256];
} IdPair;


int main(int argc, char *argv[])
{
  errno = 0;
  // reading in the args
  int c;
  while ((c = getopt(argc, argv, "i:d:o:m:")) != -1)
    {
      switch (c)
        {
        case 'i':
	  idstream = fopen(optarg, "r");
	  if (idstream == NULL)
            {
	      printf("Error opening file %s for reading.\n", optarg);
	      exit(EXIT_FAILURE);
            }
	  break;
        case 'd':
	  distmatrixstream = fopen(optarg, "r");
	  if (distmatrixstream == NULL)
            {
	      printf("Error opening file %s for reading.\n", optarg);
	      exit(EXIT_FAILURE);
            }
	  break;
        case 'o':
	  outputstream = fopen(optarg, "w");
	  if (outputstream == NULL)
            {
	      printf("Error opening file %s for writing.\n", optarg);
	      exit(EXIT_FAILURE);
            }
	  break;
        case 'm':
	  max_distance = atof(optarg);
	  if(errno != 0){
	    printf("Error converting argument %s to double.\n", optarg);
	  }
	  break;
        default:
	  printf("Invalid argument(s): %c\n", c);
	  exit(EXIT_FAILURE);
        }
    }

  fprintf(outputstream, "# command:  ");
  for(int i=0; i<argc; i++){
    fprintf(outputstream, "%s  ", argv[i]);
  }fprintf(outputstream, "\n");

  // **********  read in the id file  ************
  long array_size = INIT_ARRAY_SIZE;
  // malloc() is used to dynamically allocate a single large block of memory with the specified size.
  IdPair *id_pairs = malloc(array_size * sizeof(IdPair));
  long num_ids = 0;                                        // Keep track of the number of ids read in
  
  //  read each line from the file with getline
  char *line = NULL;
  size_t len = 0;
  __ssize_t nread;
  while (nread = getline(&line, &len, idstream) != -1) // Read each line from the file; -1 is returned when there is no more line to read
    {                                                    // getline reads an entire line from stream, storing the address of the buffer containing the text into *lineptr.
      if (num_ids >= array_size)                       // If the number of pairs exceeds the array size, double the array size
        {
	  array_size *= 2;
	  id_pairs = realloc(id_pairs, array_size * sizeof(IdPair)); // realloc() is used to dynamically change the memory allocation of a previously allocated memory.
        }
      sscanf(line, "%s %s", id_pairs[num_ids].id1, id_pairs[num_ids].id2); // Read id1 and id2 from the line
      num_ids++;
    }
  fprintf(stderr, "# Number of ids read from id file: %ld\n", num_ids);
  fclose(idstream);
  // **********  done reading id file  *********************
    
  // ***********  read in the distance matrix  *********************
  long num_cols = 0;
  long num_rows = 0;
  long row = 0;
    
  // *****  read first line and count the number of columns  *****
  nread = getline (&line, &len, distmatrixstream);
  num_rows++;
  char *token = strtok(line, " \t\n");
  long col = 0; // 0-based col number of this token

  while (token != NULL) // Read each token from the line
    {
      token = strtok(NULL, " \t\n"); // Get the next token
      col++;
      //	fprintf(stdout, "# token: [%s] col: %ld\n", token, col);
    }
  row++;
  fprintf(stdout, "# Done reading first line of matrix. Number of columns is: %ld\n", col);
  num_cols = col;
  if(num_cols != num_ids){ // check that the number of columns matches the number of ids.
    fprintf(stderr, "# Number of cols in matrix (%ld) not equal to number of ids (%ld). Bye.\n", num_cols, num_ids);
    exit(EXIT_FAILURE);
  }
  
  // ***********  read the rest of the lines  ************
  long num_distances_out = 0;
  while (nread = getline(&line, &len, distmatrixstream) != -1){
      char *token = strtok(line, " \t"); // Split the line by space or tab
      col = 0;   
      while (token != NULL) // Read each token from the line
        {
	  if (col >= row) // Only read the upper triangle of the distance matrix
  	    break;
	  // printf("row = %d, col = %d, token = %s\n", row, col, token);
	  // distances[row][col] = atof(token); // Convert the token to double and store it in the matrix. atof() is used to convert a string to a double.
	  double distance = atof(token);
	  if (distance <= max_distance)
            {
	      // printf("id1 = %s, id2 = %s, distance = %f, %ld %ld \n", id_pairs[row].id1, id_pairs[col].id1, distance, row, col);
	      fprintf(outputstream, "%s %s %f\n", id_pairs[row].id1, id_pairs[col].id1, distance);
	      num_distances_out++;
            }
	  col++;
	  token = strtok(NULL, " \t"); // Get the next token
        }
      row++;
    }
  num_rows = row;
  fprintf(stdout, "# Done reading distance matrix, %ld rows.\n", num_rows);
  free(line);

  // check if the number of ids matches the number of rows in the distance matrix
  if (num_ids != num_rows){
      fprintf(stderr, "Error: Number of rows in the distance matrix (%ld) not equal to the number of ids (%ld)\n", num_rows, num_ids);
      exit(EXIT_FAILURE);
  }else{
    fprintf(stdout, "# Read %ld accession ids and distance matrix with %ld columns and %ld rows.\n", num_ids, num_cols, num_rows);
    fprintf(stdout, "# Output %ld accession pairs with distance <= %7.5lf\n", num_distances_out, max_distance);
  }
  
  // free the memory
  free(id_pairs);
  fclose(distmatrixstream);
  fclose(outputstream);
  return 0;
};
