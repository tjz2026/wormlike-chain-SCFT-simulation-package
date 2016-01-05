#ifndef	FFT_H
#define FFT_H
#include <fftw3.h>
#include <fftw3-mpi.h>
struct fft_job {
       int dim;
       int n_fft[3];
       double *r_in;
       fftw_complex *c_out;
       ptrdiff_t alloc_local, local_n0, local_n0_start;
       fftw_plan plan_f,plan_b;  
       void create_mpi_fft_plan();
       void create_serial_fft_plan();
       void set_up_fft_grid_indx();
       void clean_fft_plan();  
       
};

extern fft_job fft_pll;
void set_up_pll_fft();
void test_fft_r2c();
void rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,unsigned long nn3, int isign);
#endif
