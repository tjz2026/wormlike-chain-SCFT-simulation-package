#ifndef	MATRIX_H
#define MATRIX_H

#include "constants.h"
#include "m_f.h"
#include "Matrix2.h"



void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void sprsin(double **a, int n, double thresh, int nmax, double sa[], int ija[]); 
void sprsax(double sa[], int ija[], double x[], double b[], int n);
void compx_sparse_matrx_mul_vector(int n, double *saR, int *ijaR, double *saI, int *ijaI, double *bR, double *bI, double *cR, double *cI);
void compx_matrx_mul_matrx(int n, double **aR, double **aI, double **bR, double **bI, double **cR, double **cI);
void compx_matrx_mul_vector(int n, double **aR, double **aI, double *bR, double *bI, double *cR, double *cI);
double norm_compx_vector_mul_vector(int n, double *aR, double *aI);
void real_matrx_mul_matrx(int n, double **a, double **b, double **c);
void inverse_matrix_real(int n, double **a_origin, double **a_origin_inverse);
void Convert(Matrix& Mat, double** DMax, int RowNum, int ColNum);
void ConvertInverse(double** DMax, const Matrix& Mat, int RowNum, int ColNum);
void inverse_matrix_complex(int n, double **a_originR, double **a_originI, double **a_originR_inverse, double **a_originI_inverse);
void tridag(double a[], double b[], double c[], double r[], double u[], unsigned long n);
void spline(double x[],double y[],int n,double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
double fact (int l);
double fact_con (int m1, int m2);
int maximum (int a1, int a2, int a3);
int minimum (int a1, int a2, int a3);
int sign (int l);
double three_j (int l1, int l2, int l3, int m1, int m2, int m3);
double Triple_product (int l1, int l2, int l3, int m1, int m2, int m3);
double simposon_1D (int n1, double delta1, double *a1);
/////////////////////////////////////////////////////////////////////////////////////////
double simposon_1D_NR(int n0, int n1, double delta1, double *a1); 
/////////////////////////////////////////////////////////////////////////////////////////
double simposon_2D (int n1,int n2,double deltalx,double deltaly,double **a2);
double Sum_sparse_2D(int ndim_i, vector<array_struct_int_int_double> Jij_1D, double *ai, double *aj);
/////////////////////////////////////////////////////////////////////////////////////////
double simposon_2D_mpi (int n1, int n2, double deltalx, double deltaly, double **a2);
double simposon_2D_1D_mpi (int n1, int n2, double deltalx, double deltaly, double *a21);
double simposon_3D_mpi (int n1, int n2, int n3, double deltalx, double deltaly, double deltalz, double ***a3);
double simposon_3D_1D_mpi (int n1, int n2, int n3, double deltalx, double deltaly, double deltalz, double *a31); 
/////////////////////////////////////////////////////////////////////////////////////////
double Sum_2D(int Ni, int Nj, double **a);
/////////////////////////////////////////////////////////////////////////////////////////
double double_dot_multi (double **a, double **b);
void inverse_matrix_NR(double **a_origin, double **a_origin_inverse,int N);
void multi_matrix2_NR(double **a, double *b, double *c, int N);
/////////////////////////////////////////////////////////////////////////////////////////
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void fourn (double data[], unsigned long nn[], int ndim, int isign);
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
#endif
