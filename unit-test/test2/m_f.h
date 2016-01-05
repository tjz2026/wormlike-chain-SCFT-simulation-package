#ifndef _M_F_H_
#define _M_F_H_
#include "constants.h"
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*standard error handler */
void nrerror(char error_text[]);

/* allocate an int vector with subscript range v[nl..nh] */
int *ivector_NR(long nl, long nh);

/* free an int vector allocated with ivector() */
void free_ivector_NR(int *v, long nl, long nh);

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector_NR(long nl, long nh);

/* free a double vector allocated with dvector() */
void free_dvector_NR(double *v, long nl, long nh);

int **imatrix_NR(long nrl, long nrh, long ncl, long nch);
void free_imatrix_NR(int **m, long nrl, long nrh, long ncl, long nch);
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix_NR(long nrl, long nrh, long ncl, long nch);

/* free a double matrix allocated by dmatrix() */
void free_dmatrix_NR(double **m, long nrl, long nrh, long ncl, long nch);

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(int nl, int nh) ;  

/* free a double vector allocated with dvector() */
void free_dvector(double *v, int nl, int nh);

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl, int nrh, int ncl, int nch);   

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
double ***f3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh); 

/* free a double f3tensor allocated by f3tensor() */ 
void free_f3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
	int ndl, int ndh);
/* allocate a int f3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
int ***f3tensor_int(int nrl, int nrh, int ncl, int nch, int ndl, int ndh); 
/* free a int f3tensor allocated by f3tensor() */ 
void free_f3tensor_int(int ***t, int nrl, int nrh, int ncl, int nch,
					   int ndl, int ndh) ; 
/* allocate a double f3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
double ****f4tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh, int nsl, int nsh);

/* free a double f4tensor allocated by f3tensor() */ 
void free_f4tensor(double ****tt, int nrl, int nrh, int ncl, int nch,
	int ndl, int ndh, int nsl, int nsh);

/* define a structure array */
struct array_struct_int_int{
	int l;
	int m;
};

struct array_struct_int_int_double{
	int i;
	int j;
	double x;
};

struct array_struct_int_int_int_double{
	int i;
	int j;
	int k;
	double x;
};

/* define structure chain for some spatial symmetry in order to keep the basis function */
typedef struct symmetry{
	double normal_coeff;
	int nx; // positive represents cos function; minus represents sin function; 0 resprsents function 1 i.e. cos(0)=1
	int ny; // positive represents cos function; minus represents sin function; 0 resprsents function 1 i.e. cos(0)=1
	int nz; // positive represents cos function; minus represents sin function; 0 resprsents function 1 i.e. cos(0)=1
	struct symmetry *next;
}SymmetryType;

#endif

