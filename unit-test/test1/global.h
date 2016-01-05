#ifndef	MDE_H
#define MDE_H
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
#include "matrix.h"
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
extern vector<array_struct_int_int> basis_function_SHF; // basis function for spherical harmonics function
extern vector<array_struct_int_int_int_double> GAMA_nonzero_1D;
extern vector<array_struct_int_int_double> J11ij_nonzero_1D;
extern vector<array_struct_int_int_double> J22ij_nonzero_1D;
extern vector<array_struct_int_int_double> J12ij_nonzero_1D;
extern vector<array_struct_int_int_double> J13ij_nonzero_1D;
extern vector<array_struct_int_int_double> J23ij_nonzero_1D;
extern vector<vector<array_struct_int_int_double> > THETA_nonzero_2D;
extern double **matrix_Rx,**matrix_Ry;
extern double ***GAMA;
struct grid {    
        int dim;
        int n_grid[3];
        int grid_num;   
        int start[3];         
        int end[3];
        double length[3];
        double dx,dy,dz; 
        void set_up_global_grid(int n, double lx, double ly, double lz);            
         };

struct chain {
   int n_blk; // block number per chain, maximum number of block is 4.
   int block_begin[4];  //starting index of each block
   int block_end[4];  // ending index of each block
   int block_len[4];  // length of each block
   int block_spe[4];  // chemical species of each block
   int N_s;  // discretization of total contour length
  // int N_sp;  // total grid number in real space r
  // int N_Mbar;  // expansion rank number of sphere harmonics function
   double ***q_f; // forward propagator
   double ***q_b; // backward propagator
   void set_up_chain( int nb, double *fb, int *sp);
   void clean_up_chain();
};

struct chemical {
    int n_spe;  // number of different chemical molecualr
    const int max_spe=3;  // three different molecualr at current version at most.
    double **W_sp;   // chemical potential of each chemical species
    double **R_sp;   // density of each chemical species
    double ****M_OP; // M order parameter potential
    double ****S_OP; // S order parameter
    void set_up_chemical( int n);
    void clean_up_chemical();
};

extern chain diblock;
extern chemical AB_melt;

extern double ***THETAij;
extern double **THETAij_M11_M22;
extern double **THETAij_M33;
extern double **THETAij_M12;
extern double **THETAij_M13;
extern double **THETAij_M23;
extern double **J11ij, **J22ij, **J12ij, **J13ij, **J23ij;

extern double   ****G_R_inverse,    ****G_I_inverse;
extern double ***sa_G_R_inverse,  ***sa_G_I_inverse;
extern int   ***ija_G_R_inverse, ***ija_G_I_inverse;

extern double pff_global,FE_global;
extern double t_diff_global[3],tm_diff_global,M_v;//M_v;
// M_v: integrate on z, that means the volume; 
//////////////////  For  MPI  //////////////////////
extern int myid, numprocs;
extern int local_nx, local_x_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;
//extern fftwnd_mpi_plan plan3d;
//extern fftwnd_mpi_plan plan3d_bak;
//extern fftwnd_mpi_plan plan2d;
//extern fftwnd_mpi_plan plan2d_bak;
//extern fftw_complex *local_data;
void create_global_variable();
//extern fftw_complex *local_data_inv;
/////////////////////////////////////////////////////////////////////////////////////////
#endif
