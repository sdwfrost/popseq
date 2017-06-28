//header file for all functions in the package
#include <R.h>
#include <Rinternals.h>






//internal function to check for duplicated combinations
void check_duplicates(int start, int ncomb, int ncol, int *comb, int *ind);

//log-factorial function
double lfactorial(int n);

//factorial function
int factorial(int n);

//choose function
int choose_c(int n, int k);


//function to generate all possible models based on an intermediate set of structures
void genfullmodels(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num);


/*This code to calculate combinations was adapted from a post by scvalex at
 http://compprog.wordpress.com/2007/10/17/generating-combinations-1/.
It's very elegant and saved me a headache---many thanks!*/
int next_comb(int *comb, int k, int n);

/*function to generate all model structures in an efficient manner
(returns an R vector of a given size, though it is not necessary to
specify the size in advance)*/
SEXP genmodels(SEXP nsamp1);


//recursive function to be called from 'genmodels'
void genmodels_recur(int nsamp, int nstart, int *indexl, int *indexl_ind, int repeatind, int *comb_next, int ncomb_next, int ntotsamp, int *currstruct, int *q, int **structure, int *maxmodels, int incmodels);


//recursive function to be called from 'genfullmodels'
void genfullmodels_recur(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num);


//indexing function for matrices
inline int index2(int i, int j, int nrow)
{
	return (nrow*j)+i;
}

//indexing function for matrices
inline int index2_col(int i, int j, int ncol)
{
	return (ncol*i)+j;
}

//indexing function for matrices
inline int index3(int i, int j, int k, int nrow, int ncol)
{
	return (k*nrow*ncol)+(nrow*j)+i;
}




//function for calculating integer powers
int pow_int(int base, int exp);

//function to perform memory reallocation of intermediate structures if required
void realloc_maxmod(int nsamp, int *maxmodels, int incmodels, int **structure);

//function to perform memory reallocation of full model set if required
void realloc_maxmod_full(int nsamp, int *maxmodels, int incmodels, int **models_num);


/*bubble sort algorithm (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec(double *values, int *ind, int n);

/*bubble sort algorithm for INTEGERS (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec_int(int *values, int *ind, int n);

/*bubble sort algorithm for INTEGERS (sorts into increasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_inc_int(int *values, int *ind, int n);
