/////////////////////////////////////////////////////////////////////////
#include "global.h"
#include "cell.h"
void density( chain *Chain, chemical *Chemical, cell *Cell )
{
	int i,j,s;
	int K,K_sp;
        int nx,ny,nz;
        double dx,dy,dz,ds;
	double pff_temp_global;
	double *ar=dvector(0,LOCAL_SIZE-1);

        nx=local_r_grid.n_grid[0];
        ny=local_r_grid.n_grid[1];
        nz=local_r_grid.n_grid[2];
        dx=Cell->dsize[0];  
        dy=Cell->dsize[1];  
        dz=Cell->dsize[2];  

	for(K=0;K<LOCAL_SIZE;K++){
		ar[K]=Chain->q_f[Chain->N_s][K][0];
	}

       double pff_temp;
       if(Cell->dim==3) 
             pff_temp=simposon_3D_1D_mpi (nx,ny,nz, dx, dy, dz, ar)/(Cell->vol);
       else if(Cell->dim==2) 
             pff_temp=simposon_2D_1D_mpi ( ny,nz, dy, dz, ar)/(Cell->vol);
       else
             pff_temp=simposon_1D_NR_pbc(0,nz,dz,ar)/(Cell->vol);

        pff_temp_global=0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&pff_temp,&pff_temp_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&pff_temp_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (myid == 0){
		pff_global=log(pff_temp_global);
	}
	
	free_dvector(ar,0,LOCAL_SIZE-1);
	
 	double totden_global;
        int NS,nch;
        NS=Chain->N_s;
        ds=1.0/NS;
	double *c_sp=dvector(0,NS);
        double  sum_sp[Chain->n_blk]; 
	//double *c_A=dvector(0,NA);
	//double *c_B=dvector(0,NB);
	double *denz=dvector(0,LOCAL_SIZE-1);

	for(K=0;K<LOCAL_SIZE;K++){
                denz[K]=0.0;
             for (j=0;j<Chemical->n_spe;j++) {
                 Chemical->R_sp[j][K]=0.0;   }                
 
             for(K_sp=0;K_sp<Chain->n_blk;K_sp++) {            
		 sum_sp[K_sp]=0.0;
		 for(i=0;i<M_Bar;i++){
			for(s=Chain->block_begin[K_sp];s<=Chain->block_end[K_sp];s++){
			c_sp[s]=Chain->q_f[s][K][i]*Chain->q_b[s][K][i];				
			                                                             }
		  sum_sp[K_sp]=sum_sp[K_sp]+simposon_1D_NR(Chain->block_begin[K_sp],Chain->block_end[K_sp],ds,c_sp);
		                     }
                 nch=Chain->block_spe[K_sp]; 
                 Chemical->R_sp[nch][K]=Chemical->R_sp[nch][K]+sum_sp[K_sp];
		 denz[K]=denz[K]+sum_sp[K_sp];
		                                  }
	                         }

        double totden;
       if(Cell->dim==3) 
	     totden=simposon_3D_1D_mpi (nx, ny, nz, dx, dy, dz, denz);
       else if(Cell->dim==2) 
	     totden=simposon_2D_1D_mpi ( ny, nz, dy, dz, denz);
       else
	     totden=simposon_1D_NR_pbc ( 0, nz, dz, denz);

        totden_global=0.0; 
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&totden,&totden_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&totden_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	
	for(K=0;K<LOCAL_SIZE;K++){
              Chemical->R_sp[Chemical->n_spe][K]=0.0;                   
           for(K_sp=0;K_sp<Chemical->n_spe;K_sp++) {            
              Chemical->R_sp[K_sp][K]=Chemical->R_sp[K_sp][K]*(Cell->vol/totden_global);
              Chemical->R_sp[Chemical->n_spe][K]=Chemical->R_sp[Chemical->n_spe][K]+Chemical->R_sp[K_sp][K];
		                                   }
	                         }



	free_dvector(c_sp,0,NS);
	free_dvector(denz,0,LOCAL_SIZE-1);



}

double free_energy_diblock( chemical *Chemical, cell *Cell  )
{
	int K,K_sp,k_sp;
	int j,k;
	double tempEE,tempEE_global;
	double *ar=dvector(0,LOCAL_SIZE-1);
        int nx,ny,nz;
        double dx,dy,dz;
        double temp1,temp2; 
	//double **b=dmatrix(0,N_dim_ddm-1,0,N_dim_ddm-1);

        nx=local_r_grid.n_grid[0];
        ny=local_r_grid.n_grid[1];
        nz=local_r_grid.n_grid[2];
        dx=Cell->dsize[0];  
        dy=Cell->dsize[1];  
        dz=Cell->dsize[2];  

         cout<<"XN[0][1]"<<Chemical->XN[0][1]<<endl; 

	for(K=0;K<LOCAL_SIZE;K++){
            ar[K]=0.0; 
	    for(K_sp=0;K_sp<Chemical->n_spe-1;K_sp++){
	       for(k_sp=K_sp+1;k_sp<Chemical->n_spe;k_sp++){
                ar[K]=ar[K]+Chemical->XN[K_sp][k_sp]*Chemical->R_sp[K_sp][K]*Chemical->R_sp[k_sp][K];
                                                           }
                                                  }

                 temp1=0.0;
                 temp2=0.0;
	    for(K_sp=0;K_sp<=Chemical->n_spe-1;K_sp++){
                ar[K]=ar[K]-Chemical->W_sp[K_sp][K]*Chemical->R_sp[K_sp][K];
                 temp1=temp1+Chemical->W_sp[K_sp][K];
                 temp2=temp2+Chemical->R_sp[K_sp][K];
                                                     }
             ar[K]=ar[K]+0.5*temp1*(temp2-1.0);   // compressibility yita
 
                                 }  

	// ?? 	ar[K]=NXab*RA[K]*RB[K]-WA[K]*RA[K]-WB[K]*RB[K]+0.5*(WA[K]+WB[K])*(RA[K]+RB[K]-1.0);	
//		if (fabs(NMu-0.0)>1.0e-4) 
//		{
//			for(j=0;j<N_dim_ddm;j++)
//				for(k=0;k<N_dim_ddm;k++)
//				{
//					b[j][k]=M_OP[K][j][k];
//				}
//				ar[K]=ar[K]+0.5*(double_dot_multi(b, b)/NMu);
//		}
	

       if(Cell->dim==3) 
	     tempEE=simposon_3D_1D_mpi ( nx, ny, nz, dx,dy, dz, ar);
       else if(Cell->dim==2) 
	     tempEE=simposon_2D_1D_mpi ( ny, nz, dy, dz, ar);
       else
	     tempEE=simposon_1D_NR_pbc(0, nz, dz, ar);

        tempEE_global=0.0;
	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&tempEE,&tempEE_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
          
        
	if (myid == 0)	FE_global=-pff_global+tempEE_global/Cell->vol;
	if (myid == 0)	cout<<"FE_tempEE"<<tempEE_global/Cell->vol<<endl;
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&FE_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	free_dvector(ar,0,LOCAL_SIZE-1);
	//free_dmatrix(b,0,N_dim_ddm-1,0,N_dim_ddm-1);
	
	return FE_global;
}

/////////////////////////////////////////////////////////////////////////
