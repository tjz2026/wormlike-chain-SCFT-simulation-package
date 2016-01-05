#include "global.h"
#include "fft.h"

fft_job fft_pll;

void fft_job::create_mpi_fft_plan()  // for mpi fftw 2D & 3D
{
/* get local data size and allocate */
      int dim; 
      dim=global_grid.dim;
      cout<<"dim of fft"<<dim<<endl;
      assert(dim >= 2 && dim <= 3); 
       int K,L,M,N;  
       for(K=0;K<3;K++){
         n_fft[K]=global_grid.n_grid[K];
                         }  
         L=n_fft[0]; 
         M=n_fft[1]; 
         N=n_fft[2];
         cout<<"L-M-N"<<L<<" "<<M<<" "<<N<<endl; 
// Note that for r2c fft, the input real-data array must be padded with extra values to accommodate the size of the complex data
// See Ref: http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
// Do not forget this when doing assignment!!! 
         if(dim==2) {
         alloc_local = fftw_mpi_local_size_2d(M, N/2+1, MPI_COMM_WORLD,
                                               &local_n0, &local_n0_start);
         r_in = fftw_alloc_real(2 * alloc_local);
         c_out = fftw_alloc_complex(alloc_local);
         /* create plan for out-of-place r2c DFT */
         plan_f = fftw_mpi_plan_dft_r2c_2d(M, N, r_in, c_out, MPI_COMM_WORLD,
                                         FFTW_MEASURE);
         plan_b = fftw_mpi_plan_dft_c2r_2d(M, N, c_out, r_in, MPI_COMM_WORLD,
                                         FFTW_MEASURE);}
         else  {  
         alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, MPI_COMM_WORLD,
                                               &local_n0, &local_n0_start);
         r_in = fftw_alloc_real(2 * alloc_local);
         c_out = fftw_alloc_complex(alloc_local);
         plan_f = fftw_mpi_plan_dft_r2c_3d(L, M, N, r_in, c_out, MPI_COMM_WORLD,
                                         FFTW_MEASURE);
         plan_b = fftw_mpi_plan_dft_c2r_3d(L, M, N, c_out, r_in, MPI_COMM_WORLD,
                                         FFTW_MEASURE);}

// SLAB decomposition, sliced in X axis.
// Note that the physical size of r_in is not equal to LOCAL_SIZE, but the size of c_out does equal to LOCAL_SIZE_K,
// make sure you understand this.
// set up local grid structure in real and fourier space.
         if(dim==2) { LOCAL_SIZE=local_n0*n_fft[2];
                    LOCAL_SIZE_K=local_n0*(n_fft[2]/2+1);
             local_r_grid.start[1]=local_n0_start;      
             local_r_grid.start[0]=0;      
             local_r_grid.start[2]=0;      
             local_r_grid.end[1]=local_n0_start+local_n0-1;      
             local_r_grid.end[0]=n_fft[0]-1;      
             local_r_grid.end[2]=n_fft[2]-1; 

             local_r_grid.n_grid[0]=n_fft[0]; 
             local_r_grid.n_grid[1]=local_n0; 
             local_r_grid.n_grid[2]=n_fft[2]; 

                      
        }
 
         else  { LOCAL_SIZE=local_n0*n_fft[1]*n_fft[2];
                    LOCAL_SIZE_K=local_n0*(n_fft[2]/2+1)*n_fft[1];
             local_r_grid.start[0]=local_n0_start;      
             local_r_grid.start[1]=0;      
             local_r_grid.start[2]=0;      
             local_r_grid.end[0]=local_n0_start+local_n0-1;      
             local_r_grid.end[1]=n_fft[1]-1;      
             local_r_grid.end[2]=n_fft[2]-1;      

             local_r_grid.n_grid[1]=n_fft[1]; 
             local_r_grid.n_grid[0]=local_n0; 
             local_r_grid.n_grid[2]=n_fft[2]; 

               }
             local_r_grid.dim=dim;
             local_r_grid.grid_num=LOCAL_SIZE;

             local_r_grid.set_up_1d_indx();

}


void fft_job::set_up_fft_grid_indx()

{

//That is, we only store the lower half (non-negative frequencies), plus one element, of the last dimension of the data from the ordinary complex transform.
// this means that we only store roughly half the size of last dimension.

             local_k_grid.dim=local_r_grid.dim;
             local_k_grid.grid_num=LOCAL_SIZE_K;

            cout<<"Local_size of k"<<LOCAL_SIZE_K<<"on"<<myid<<endl; 
             if (local_k_grid.dim==2) {       
             local_k_grid.start[1]=local_n0_start;      
             local_k_grid.start[0]=0;      
             local_k_grid.start[2]=0;      
             local_k_grid.end[1]=local_n0_start+local_n0-1; 
             local_k_grid.end[2]=n_fft[2]/2;      
             local_k_grid.end[0]=n_fft[0]-1;
              
             local_k_grid.n_grid[0]=n_fft[0];
             local_k_grid.n_grid[1]=local_n0;
             local_k_grid.n_grid[2]=n_fft[2]/2+1;
                                      }
             else  {
             local_k_grid.start[0]=local_n0_start;      
             local_k_grid.start[1]=0;      
             local_k_grid.start[2]=0;      
             local_k_grid.end[0]=local_n0_start+local_n0-1; 
             local_k_grid.end[1]=n_fft[1]-1;      
             local_k_grid.end[2]=n_fft[2]/2;
             local_k_grid.n_grid[1]=n_fft[0];
             local_k_grid.n_grid[0]=local_n0;
             local_k_grid.n_grid[2]=n_fft[2]/2+1;
                  } 

        local_k_grid.set_up_1d_indx();

        int K;
        vector<array_struct_int_int_int>::iterator iter;
    for( iter = local_k_grid.indx_1d.begin(); iter != local_k_grid.indx_1d.end(); iter++ )
        {
             if(n_fft[0]!=1 && (*iter).i>(n_fft[0]/2)) {
             (*iter).i=((*iter).i)%(n_fft[0]/2) -(n_fft[0]/2);}
             if(n_fft[1]!=1 && (*iter).j>(n_fft[1]/2)) {
             (*iter).j=((*iter).j)%(n_fft[1]/2) -(n_fft[1]/2);}
             if(n_fft[2]!=1 && (*iter).k>(n_fft[2]/2)) {
             (*iter).k=((*iter).k%(n_fft[2]/2) -(n_fft[2]/2));}


     //       cout<<"indx_1d.i="<<(*iter).i<<endl;
     //       cout<<"indx_1d.j="<<(*iter).j<<endl;
     //       cout<<"indx_1d.k="<<(*iter).k<<endl;
        }


}


void fft_job::create_serial_fft_plan()  // for 1D,2D,3D fft on single cpu processor
{
/* get local data size and allocate */

      dim=global_grid.dim;
      assert(dim >= 2 && dim <= 3); 
       int K,L,M,N;  
       for(K=0;K<dim;K++){
         n_fft[K]=global_grid.n_grid[K];
                         }  
         L=n_fft[0]; 
         M=n_fft[1]; 
         N=n_fft[2]; 
/*
 to be implemented soon!
 */
};

  void fft_job::clean_fft_plan()
{
fftw_destroy_plan(plan_f);
fftw_destroy_plan(plan_b);
}


void set_up_pll_fft()

{ 
  fft_pll.create_mpi_fft_plan();
  fft_pll.set_up_fft_grid_indx();
 
}


void test_fft_r2c()

{ 
  
 int i,j,k;
 int i_i,j_j,k_k,kk;
 int n1,n2,n3;
    n1=local_r_grid.n_grid[0];
    n2=local_r_grid.n_grid[1];
    n3=local_r_grid.n_grid[2];

 for (i=local_r_grid.start[0],kk=0;i<=local_r_grid.end[0];i++) {
    for (j=local_r_grid.start[1];j<=local_r_grid.end[1];j++)      {
// why k<=local_r_grid.end[2]+2 ??, because we assume SIDEz to be even.
       for (k=local_r_grid.start[2];k<=local_r_grid.end[2]+2;k++,kk++) {
            i_i=i-local_r_grid.start[0];            
            j_j=j-local_r_grid.start[1];            
            k_k=k-local_r_grid.start[2];
            fft_pll.r_in[i_i*(n2*(n3+2))+j_j*(n3+2)+k_k]=(k+2*j);
                                                                     }    
                                                                   }
                                                                 } 


      fftw_execute(fft_pll.plan_f);  

      fftw_execute(fft_pll.plan_b);  
      for (kk=0;kk<n1*n2*(n3+2);kk++) {
         fft_pll.r_in[kk]=fft_pll.r_in[kk]/(1.0*local_r_grid.grid_num); } 
  
 
}


