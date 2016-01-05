////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
#include "global.h"
#include "minimal.h"
#include <ctime>

#define tol 1.0e-2
int LOCAL_SIZE;
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
	LOCAL_SIZE=SIDEx*SIDEy;


//	initial();
//	initw();
	double ax,bx,cx;
//      readw();
	

	ax=1.0/double(SIDEx);
	bx=1.3/double(SIDEx);
	cx=1.8/double(SIDEx);
	double minFE;
	double minx;
       //minFE=brent(ax,bx,cx,cal_scft,tol,&minx);
	
	create_global_variable();

	//fftwnd_mpi_destroy_plan (plan2d);
	//fftwnd_mpi_destroy_plan (plan2d_bak);
	MPI_Finalize ();

        double speed=(double)(clock()-start)/CLOCKS_PER_SEC;
	cout<<speed<<"seconds"<<endl;		
	return 0;
}
