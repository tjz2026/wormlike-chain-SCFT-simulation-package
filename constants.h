#ifndef	_CONSTANTS_H_
#define _CONSTANTS_H_
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <vector>
//#include <fftw3.h>
//#include <fftw3_mpi.h>
#include <mpi.h>

using namespace std;
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#define CC 2e-4
#define CC_FINAL 1e-4

#define CC_M 1e-2  
#define CC_M_FINAL 1e-2
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#define LLS 1e-8  
#define TINY 1.0e-40;
#define NR_END 1
#define FREE_ARG char*
#define MAXITS 10000
#define EPS 1.0e-10
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6
#define ALF 1.0e-4
#define PI 3.141592653589792
#define N_dim_ddm 3 // the dimension for the tensor of double dot multiply (ddm) on the basis of Maier-Saupe potential
////////////////////////////////////////////////////
//////////////////////////////////////////////////////

#define L_Bar 4
#define M_Bar (L_Bar+1)*(L_Bar+1)
//#define L_Bar_A 6
//#define M_Bar_A (L_Bar_A+1)*(L_Bar_A+1)
//#define L_Bar_B 6
//#define M_Bar_B (L_Bar_B+1)*(L_Bar_B+1)

#define Thresh_sprase_matrix 1.0e-10

#define Dim_ordering  (M_Bar*M_Bar+L_Bar)
/////////////////////////////////////////////////////
const int NDIM=1;
///////////////////////////////////////////////////
// chain parameters

// Cell parameters



// grid number

// model parameters


// SCFT control parameters 
const int DIMENSION=2;
const double NN=5.0;  // polymerization
const int SIDEx=1;
const int SIDEy=64;
const int SIDEz=32;
const double M_grid=1.0*SIDEx*SIDEy*SIDEz; // the number of grids for the spacial division
extern int LOCAL_SIZE;
const int NMAX=100;
const double ds=1.0/NMAX;
const double kapa_A=0.5;
const double kapa_B=0.5;  
extern double dx;
extern double dy;
////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/////////////////// Anderson Mixing ////////////////////////////////////////
const int n_r_WAB=10; // the steps including the preceding n_r_WAB historic steps
const int Num_step_simple_mixing_for_WAB=30;
extern double lambda_WAB_anderson;
const double lambda_WAB_anderson_const=0.1;


const int Num_step_simple_mixing_for_M=500000;
const int Num_step_iterae_M=4;  // update the orientational field "M" after every Num_step_iterae_M steps for anderson mixing
const int n_r_M=3;
extern double lambda_M_anderson;
const double lambda_M_anderson_const=0.1;
extern int Num_iteration_step_M; // counting the practical number of iterated M field
const char method='B';
// NOTE, uniaxial phase i.e. M_xx  = M_yy = -(1/2)M_zz; xx-->00, yy-->11, zz-->22
// NOTE, biaxial phase i.e. M_xx != M_yy ; xx-->00, yy-->11, zz-->22
// A represents uniaxial phase and only elements on the diagonal are nonezero (i.e. Smectic A); 
// B represents general uniaxial phase (i.e. Smectic C);
// not A or not B represents biaxial phase;
///////////////////////////////////////////////////
const double fa=0.3;
const double fb=1.0-fa;
const int NA=int(NMAX*fa);
const double lanbtWA=0.2;
const double lanbtWB=0.2;
const double lanbtM=0.0;
extern int Num_iteration_step_WM; // counting the number of total iteration steps
const double NXab=16.00;
const double NXac=0.00;
const double NXbc=0.00;
const double NMu_NXab=1.0; // NMu/NXab
const double NMu=0.0*(0.6666666)*NMu_NXab;
const double M_initial=36.974695*(0.6666666)*NMu_NXab;
///////////////////////////////////////////////////
#endif
