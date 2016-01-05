#include <mpi.h>
#include <fftw3-mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include <fstream>
 using namespace std;    
     int main(int argc, char **argv)
     {
         const ptrdiff_t N0 = 8, N1 = 4;
         fftw_plan plan_f,plan_b;
         fftw_complex *data;
         fftw_complex *c_out;
         double *r_in;
         ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
     
         MPI_Init(&argc, &argv);
         fftw_mpi_init();
     
         /* get local data size and allocate */
  //       alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
  //                                            &local_n0, &local_0_start);
  //       data = fftw_alloc_complex(alloc_local);
     
         /* create plan for in-place forward DFT */
  //       plan_f = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
  //                                   FFTW_FORWARD, FFTW_ESTIMATE);
     
 //        plan_b = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
   //                                  FFTW_BACKWARD, FFTW_ESTIMATE);


         alloc_local = fftw_mpi_local_size_2d(N0, N1/2+1, MPI_COMM_WORLD,
                                               &local_n0, &local_0_start);
         cout<<"alloc_local"<<alloc_local<<endl;
         cout<<"local_n0"<<local_n0<<endl;
         r_in = fftw_alloc_real(2 * alloc_local);
         c_out = fftw_alloc_complex(alloc_local);
         /* create plan for out-of-place r2c DFT */
         plan_f = fftw_mpi_plan_dft_r2c_2d(N0, N1, r_in, c_out, MPI_COMM_WORLD,
                                         FFTW_MEASURE);
         plan_b = fftw_mpi_plan_dft_c2r_2d(N0, N1, c_out, r_in, MPI_COMM_WORLD,
                                         FFTW_MEASURE);

         /* initialize data to some function my_function(x,y) */
         //for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
         //   data[i*N1 + j][0] = i+2*j;
         //   data[i*N1 + j][1] = 0.0;
     
         for (i = 0; i < local_n0; ++i) for (j = 0; j < N1+2; ++j) {
            r_in[i*(N1+2) + j] = 2*i+j;}

         /* compute transforms, in-place, as many times as desired */
         fftw_execute(plan_f);
        // for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
        //    cout<<"data r"<<data[i*N1 + j][0]<<endl;

        // for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
        //    cout<<"data i"<<data[i*N1 + j][1]<<endl;

     
      /*   for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
           cout<<"c_out r"<<c_out[i*N1 + j][0]<<endl;

         for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
            cout<<"c_out i"<<c_out[i*N1 + j][1]<<endl;
      */
            fftw_execute(plan_b);
         for (i = 0; i < local_n0; ++i) for (j = 0; j < N1+2; ++j) {
              r_in[i*(N1+2)+j]=r_in[i*(N1+2)+j]/16.0; }
         fftw_destroy_plan(plan_f);
         fftw_destroy_plan(plan_b);


     
         MPI_Finalize();
     }
