////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
#include "global.h"
#include "cell.h"
#include "initial.h"
#include "fft.h"
#include "G_matrix.h"
#include "minimal.h"
#include <ctime>
#include "AB_diblock_driver.h"
#define tol 1.0e-2
//int LOCAL_SIZE;
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){


    clock_t start;
	start=clock();

	FILE *fp;
	
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	
        MPI_Barrier(MPI_COMM_WORLD);

//	plan2d = fftw2d_mpi_create_plan (MPI_COMM_WORLD, SIDEx, SIDEy, FFTW_FORWARD, FFTW_MEASURE);
//	plan2d_bak = fftw2d_mpi_create_plan (MPI_COMM_WORLD, SIDEx, SIDEy, FFTW_BACKWARD, FFTW_MEASURE);
//	fftwnd_mpi_local_sizes (plan2d, &local_nx, &local_x_start, 
//		&local_ny_after_transpose, &local_y_start_after_transpose, &total_local_size);

//	local_data = (fftw_complex*) malloc(sizeof(fftw_complex)*total_local_size);

/////////////////////////////////////////////////////////////////////////////////////////
/*
	if(fp=fopen("mpi.dat","a"))
	{
		fprintf(fp,"local_nx = %d\n local_x_start = %d\n local_ny_after_transpose = %d local_y_start_after_transpose = %d\n total_local_size = %d\n",local_nx, local_x_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size );
		fclose(fp);
	}
*/	


//	initial();
//	initw();
	double ax,bx,cx;
//      readw();
	

	ax=1.0/double(SIDEx);
	bx=1.3/double(SIDEx);
	cx=1.8/double(SIDEx);
	double minFE;
	double minx;
  double **G_R_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
  double **G_I_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
       //minFE=brent(ax,bx,cx,cal_scft,tol,&minx);
	
	create_global_grid();
        HEX.init_cell( true, false, 4, 0.0, 1.040654, 0.60082 );
       // set up parallel fft job and build local r and k grid.
        set_up_pll_fft();
        test_fft_r2c();
        create_chain_chemistry();
        initial(AB_melt); 

        matrix_G_inverse( 1, 0, 0.5, matrix_Rx, matrix_Ry, matrix_Rz, \
                G_R_inv,G_I_inv,global_grid, diblock, HEX );

       // cout<<"matrix_Rx"<<"0,0 "<<matrix_Rx[0][0]<<" 1,1"<<matrix_Rx[1][1]<<endl;
       // cout<<"G_R_inv"<<"0,0x "<<G_R_inv[0][0]<<endl;
       // cout<<"G_R_inv"<<"1,1x "<<G_R_inv[1][1]<<endl;
        //cout<<"G_R_inv"<<"2,2x "<<G_R_inv[2][2]<<endl;
         
        //exit(1);  
        compress_G_matrix( global_grid, AB_melt, diblock, sa_G_R_inverse, sa_G_I_inverse, \ 
           ija_G_R_inverse, ija_G_I_inverse, HEX);
// scft loop
        AB_diblock_scft(4, 3000, CC);

//prepare to finalize

	//fftwnd_mpi_destroy_plan (plan2d);
	//fftwnd_mpi_destroy_plan (plan2d_bak);
	MPI_Finalize ();

        double speed=(double)(clock()-start)/CLOCKS_PER_SEC;
	cout<<speed<<"seconds"<<endl;		
	free_dmatrix(G_R_inv,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(G_I_inv,0,M_Bar-1,0,M_Bar-1);
	return 0;
}
