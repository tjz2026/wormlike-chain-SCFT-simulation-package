/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
#include "global.h"
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
vector<array_struct_int_int> basis_function_SHF; // basis function for spherical harmonics function
vector<array_struct_int_int_int_double> GAMA_nonzero_1D;
vector<array_struct_int_int_double> J11ij_nonzero_1D;
vector<array_struct_int_int_double> J22ij_nonzero_1D;
vector<array_struct_int_int_double> J12ij_nonzero_1D;
vector<array_struct_int_int_double> J13ij_nonzero_1D;
vector<array_struct_int_int_double> J23ij_nonzero_1D;
vector<vector<array_struct_int_int_double> > THETA_nonzero_2D;

    void grid::set_up_global_grid(int n, double lx, double ly, double lz)
{        
      dim=n;   
      assert(dim >= 1 && dim <= 3); 
      if (dim==1) {
          n_grid[0]=SIDEx;
          n_grid[1]=1;
          n_grid[2]=1;
          dx=lx/SIDEx;
          dy=0.0;
          dz=0.0;
                  }      
 
      if (dim==2) {
          n_grid[0]=SIDEx;
          n_grid[1]=SIDEy;
          n_grid[2]=1;
          dx=lx/SIDEx;
          dy=ly/SIDEy;
          dz=0.0;
                  }  
    
      if (dim==3) {
          n_grid[0]=SIDEx;
          n_grid[1]=SIDEy;
          n_grid[2]=SIDEz;
          dx=lx/SIDEx;
          dy=ly/SIDEy;
          dz=lz/SIDEz;
                  } 
     
}

grid global_grid;
grid local_grid;

   void chain::set_up_chain(int nb, double *fb, int *sp)
{
      n_blk=nb;
      N_s=NMAX;  
      assert(n_blk > 1 && n_blk <= 4);  // diblcok, triblock, and four blocks. 
       int K;
       int start_indx=0;    
       for(K=0;K<n_blk-1;K++){
         block_spe[K]=sp[K]; 
         block_begin[K]=start_indx; 
         block_len[K]=int(N_s*fb[K]);
         start_indx=start_indx+block_len[K];    
         block_end[K]=start_indx; 
                              }
         block_spe[n_blk-1]=sp[n_blk-1]; 
         block_begin[n_blk-1]=block_end[n_blk-2]; 
         block_len[n_blk-1]=N_s-block_begin[n_blk-1];
         block_end[K]=N_s; 
	 q_f=f3tensor(0,N_s,0,LOCAL_SIZE-1,0,M_Bar-1);
	 q_b=f3tensor(0,N_s,0,LOCAL_SIZE-1,0,M_Bar-1);
         printf("sizeof(a): %d\n", sizeof(q_f));



}


   void chain::clean_up_chain()
{
    free_f3tensor(q_f,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
    free_f3tensor(q_b,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
    cout<<" q_f & q_b are deleted "<<endl;
}


    void chemical::set_up_chemical(int n)   
{   
     n_spe=n;
     assert(n_spe > 1 && n_spe <= 3);
     W_sp=dmatrix(0,n_spe-1,0,LOCAL_SIZE-1);
     R_sp=dmatrix(0,n_spe,0,LOCAL_SIZE-1);
     M_OP=f4tensor(0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
     S_OP=f4tensor(0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);

}

    void chemical::clean_up_chemical()   
{
     free_dmatrix(W_sp,0,n_spe-1,0,LOCAL_SIZE-1);
     free_dmatrix(R_sp,0,n_spe,0,LOCAL_SIZE-1);
     free_f4tensor(M_OP,0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
     free_f4tensor(S_OP,0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
      
}

chain diblock; // set up a diblock chain;
chemical AB_melt; // set up a AB melt.


double **matrix_Rx,**matrix_Ry;
double ***GAMA;
double ***THETAij;
double **THETAij_M11_M22;
double **THETAij_M33;
double **THETAij_M12;
double **THETAij_M13;
double **THETAij_M23;
double **J11ij, **J22ij, **J12ij, **J13ij, **J23ij;
    
double   **G_R_inverse,    **G_I_inverse;
double ***sa_G_R_inverse,  ***sa_G_I_inverse;
int   ***ija_G_R_inverse, ***ija_G_I_inverse;

double pff_global,FE_global;
double t_diff_global[3];
double tm_diff_global,M_v;//M_v;

//double ta_diff_global,tb_diff_global,tm_diff_global,M_v;//M_v;

double lambda_WAB_anderson;
int Num_iteration_step_WM;
double lambda_M_anderson;
int Num_iteration_step_M;

//double dx;
//double dy;

//////////////////  For  MPI  //////////////////////
int myid, numprocs;
int local_nx, local_x_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;
char name1[80],name2[80],name3[80],name4[80],name5[80],name6[80],name7[80],name8[80],name9[80],name10[80];
//fftwnd_mpi_plan plan3d;
//fftwnd_mpi_plan plan3d_bak;
//fftwnd_mpi_plan plan2d; // fftw forward transform
//fftwnd_mpi_plan plan2d_bak; // fftw backward transform
//fftw_complex *local_data;  // calculate the propagator from the positive direction
//fftw_complex *local_data_inv;  // calculate the propagator from the another direction
/////////////////////////////////////////////////////////////////////////////////////////////

  void create_global_variable()
{
   
  global_grid.set_up_global_grid( 2, 1.0,1.0,0.0);
  
  cout<<"grid dim="<<global_grid.dim<<endl;
  cout<<"grid nx="<<global_grid.n_grid[0]<<endl;
  cout<<"local_size"<<" "<<LOCAL_SIZE<<endl;
  AB_melt.set_up_chemical(2);

  double fb[2];
  fb[0]=0.3; 
  fb[1]=0.7;
  int sp[2];
  sp[0]=0; 
  sp[1]=1;

  diblock.set_up_chain(2, fb, sp);
  cout<<"chain Ns="<<diblock.N_s<<endl;

}


