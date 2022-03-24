// ********  various 'vectors' arrays which know their size ****
// ********  and (some of them) allow adding elements (a la 'push')
// ********  with realloc to increase capacity if necessary.

#define DBUG 1

// *********  typedefs  ********

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  long* a; // array
} Vlong;

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  char** a; // array of strings
} Vstr;

typedef struct{
  long capacity;
  long length; // length of the string not including term. null
  char* a; // a regular null-terminated string
}Vchar; 

typedef struct{
  long index;
  char* id;
} IndexId;

typedef struct{
  long capacity;
  long size;
  IndexId**a;
} Vidxid;


// *****  Function declarations  ************************************************************

// ***** Vlong ******************************************************************************
Vlong* construct_vlong(long cap); // set capacity = cap, size = 0
Vlong* construct_vlong_zeroes(long size);
Vlong* construct_vlong_from_array(long size, long* array); // initialize with array of known size
Vlong* construct_vlong_whole_numbers(long size); // initialize to 0,1,2,3,...size-1
void add_long_to_vlong(Vlong* the_vlong, long x); // push, realloc if necessary
void shuffle_vlong(Vlong* the_vlong); // randomize order of array elements
void free_vlong(const Vlong* the_vlong); // free memory

// *****  Vstr  *****************************************************************************
Vstr* construct_vstr(long cap); // set capacity = cap, size = 0
Vstr* construct_vstr_copy(Vstr* the_vstr);
void add_string_to_vstr(Vstr* the_vstr, char* str); // push, realloc if necessary
char* ith_str_from_vstr(Vstr* the_vstr, long i); // perl-like: index -1 -> last element, etc.
char* copy_ith_str_from_vstr(Vstr* the_vstr, long i); // perl-like: index -1 -> last element, etc.
void print_vstr(FILE* fh, Vstr* the_vstr);
void free_vstr(const Vstr* the_vstr); // free memory

// *****  Vchar  *****
Vchar* construct_vchar(long cap);
Vchar* construct_vchar_from_str(char* str); // str is null-terminated str
Vchar* copy_vchar(Vchar* avchar);
Vchar* append_str_to_vchar(Vchar* the_vchar, char* str);
void print_vchar(FILE* fh, Vchar* the_vchar);
void free_vchar(const Vchar* the_vchar);

// *****  IndexId  *****
IndexId* construct_indexid(long idx, char* id);
void free_indexid(const IndexId* the_idxid);

// *****  Vidxid  *****
/* Vidxid* construct_vidxid_from_vstr(Vstr* ids); */
/* Vidxid* construct_sorted_vidxid_from_vstr(Vstr* ids); */
int strcmpx(const void* v1, const void* v2);
void sort_vidxid_by_id(Vidxid* the_vidxid);
long index_of_id_in_vidxid(Vidxid* the_vidxid, char* id);
void print_vidxid(FILE* fh, Vidxid* the_vidxid);
void free_vidxid(const Vidxid* the_vidxid);


