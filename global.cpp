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

          
    void grid::set_up_global_grid(int n )
{        
      dim=n;   
      assert(dim >= 1 && dim <= 3); 
// for 1d case, nx=1,ny=1,nz=SIDEz, for 2d case, nx=1,ny=SIDEy,nz=SIDEz.
      if (dim==1) {
          n_grid[2]=SIDEz;
          n_grid[1]=1;
          n_grid[0]=1;
          //dz=lz/SIDEz;
          //dy=0.0;
          //dx=0.0;
          //vol=SIDEz*dz;
          grid_num=SIDEz;  

                  }      
 
      if (dim==2) {
          n_grid[1]=SIDEy;
          n_grid[2]=SIDEz;
          n_grid[0]=1;
          //dy=ly/SIDEy;
          //dz=lz/SIDEz;
          //dx=0.0;
          //vol=SIDEy*dy*SIDEz*dz;  
          grid_num=SIDEy*SIDEz;  
                  }  
    
      if (dim==3) {
          n_grid[0]=SIDEx;
          n_grid[1]=SIDEy;
          n_grid[2]=SIDEz;
          //dx=lx/SIDEx;
          //dy=ly/SIDEy;
          //dz=lz/SIDEz;
          //vol=SIDEx*dx*SIDEy*dy*SIDEz*dz;  
          grid_num=SIDEx*SIDEy*SIDEz;  
                  } 
     
}

    void grid::set_up_1d_indx()  // do not use this function for global_grid!!
{
    int K_i,K_j,K_k; 
    array_struct_int_int_int test;

    for(K_i=start[0];K_i<=end[0];K_i++){
        for(K_j=start[1];K_j<=end[1];K_j++){
             for(K_k=start[2];K_k<=end[2];K_k++){
                test.i=K_i;    
                test.j=K_j;    
                test.k=K_k;
                indx_1d.push_back(test);
                                                 }   
                                            }
                                        }
//    cout<<"loc grid_num="<<grid_num<<endl; 
//    cout<<"indx_1d.size"<<indx_1d.size()<<endl; 
//    cout<<"start[0],end[0]"<<start[0]<<" "<<end[0]<<endl; 
//    cout<<"start[1],end[1]"<<start[1]<<" "<<end[1]<<endl; 
//    cout<<"start[2],end[2]"<<start[2]<<" "<<end[2]<<endl; 

    assert(indx_1d.size()==grid_num);                    

}

grid global_grid;
grid local_r_grid;
grid local_k_grid;

   void chain::set_up_chain(int nn, int nb, double *fb, int *sp, chemical &Chemical )
{

      Loa=nn;  
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
         block_kapa[K]=Chemical.spe_kapa[block_spe[K]]; 
                              }
// the last block is treated specially.
         block_spe[n_blk-1]=sp[n_blk-1]; 
         block_begin[n_blk-1]=block_end[n_blk-2]; 
         block_len[n_blk-1]=N_s-block_begin[n_blk-1];
         block_end[n_blk-1]=N_s; 
         block_kapa[n_blk-1]=Chemical.spe_kapa[block_spe[n_blk-1]]; 

	 q_f=f3tensor(0,N_s,0,LOCAL_SIZE-1,0,M_Bar-1);
	 q_b=f3tensor(0,N_s,0,LOCAL_SIZE-1,0,M_Bar-1);
         printf("sizeof(a): %d\n", sizeof(q_f));
         for (K=0; K<n_blk;K++) {
         cout<<"blk_begin"<<block_begin[K]<<" for"<<K<<" th block"<<endl; 
         cout<<"blk_end"<<block_end[K]<<" for"<<K<<" th block"<<endl; }


}


   void chain::clean_up_chain()
{
    free_f3tensor(q_f,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
    free_f3tensor(q_b,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
    cout<<" q_f & q_b are deleted "<<endl;
}


    void chemical::set_up_chemical(int n, double *kapa,double **X_N)   
{   
     n_spe=n;
     assert(n_spe >= 2 && n_spe <= 3);
     int k,j;
     for (k=0;k<n_spe;k++) {
     spe_kapa[k]=kapa[k];
        for (j=0;j<n_spe;j++) {
         XN[k][j]=X_N[k][j];
                              } 
                           } 

     cout<<"chemical::set_up:: local_size="<<LOCAL_SIZE<<endl;
     cout<<"chemical::set_up::n_spe="<<n_spe<<endl;
     W_sp=dmatrix(0,n_spe-1,0,LOCAL_SIZE-1);
     R_sp=dmatrix(0,n_spe,0,LOCAL_SIZE-1);
     M_OP=f3tensor(0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
     S_OP=f4tensor(0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);

}

    void chemical::clean_up_chemical()   
{
     free_dmatrix(W_sp,0,n_spe-1,0,LOCAL_SIZE-1);
     free_dmatrix(R_sp,0,n_spe,0,LOCAL_SIZE-1);
     free_f3tensor(M_OP,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
     free_f4tensor(S_OP,0,n_spe,0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
      
}

chain diblock; // set up a diblock chain;
chemical AB_melt; // set up a AB melt.


int LOCAL_SIZE;
int LOCAL_SIZE_K;

double **matrix_Rx,**matrix_Ry,**matrix_Rz;
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

  void create_global_grid()
{
   
  global_grid.set_up_global_grid( DIMENSION );
  
  cout<<"grid dim="<<global_grid.dim<<endl;
  cout<<"grid nx="<<global_grid.n_grid[0]<<endl;
  cout<<"grid ny="<<global_grid.n_grid[1]<<endl;
  cout<<"grid nz="<<global_grid.n_grid[2]<<endl;

}


  void create_chain_chemistry()
{
   
  double kapa[2];
  kapa[0]=0.5; 
  kapa[1]=0.5;
  double **X_N;
  X_N=dmatrix(0,2,0,2);
  X_N[0][0]=0.0;  
  X_N[0][1]=NXab;  
  X_N[0][2]=NXac;  
  X_N[1][0]=NXab;  
  X_N[1][1]=0.0;  
  X_N[1][2]=NXbc;  
  X_N[2][0]=NXac;  
  X_N[2][1]=NXbc;  
  X_N[2][2]=0.0;  
  AB_melt.set_up_chemical(2, kapa, X_N);
  free_dmatrix(X_N,0,2,0,2);
  double fb[2];
  fb[0]=0.3; 
  fb[1]=0.7;
  int sp[2];
  sp[0]=0; 
  sp[1]=1;

  diblock.set_up_chain(NN, 2, fb, sp, AB_melt);
  cout<<"chain Ns="<<diblock.N_s<<endl;

}


