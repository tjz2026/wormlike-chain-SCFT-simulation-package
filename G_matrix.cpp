#include "global.h"
#include "G_matrix.h"

void matrix_G_inverse( int blk_step, int K_th, double kapa_temp, double **matrix_Rxx, \ 
     double **matrix_Ryy,double **matrix_Rzz, double **G_R_inv, \
     double **G_I_inv, grid &Grid, chain &Chain, cell &Cell ) 
// blk_step=+-(1,2,3).  blk_step =+-1 or +-2 represents the first two steps along the chain. "+/-" for forward/backward propagator.
// K_th : the 1d indx of local k space grid.

//Grid is the class references for global_grid
{// coeff_Rij=1.0 for G; coeff_Rij=-1.0 for G_star
	int i,j;
	int K_i,K_ii;
	int K_j,K_jj;
	double **G_R=dmatrix(0,M_Bar-1,0,M_Bar-1);
	double **G_I=dmatrix(0,M_Bar-1,0,M_Bar-1);
	//double **G_R_inverse=dmatrix(0,M_Bar-1,0,M_Bar-1);
	//double **G_I_inverse=dmatrix(0,M_Bar-1,0,M_Bar-1);
	const double eps=1.0e-10;
        int K, NN;
        double Kx,Ky,Kz;
        double d_s,dx,dy,dz; 
        double coeff_deltaij, coeff_Rij;

      
      assert(abs(blk_step)<=3 && blk_step!=0);
//coefficients for  BDF3 algorithm.
        if(blk_step>0) 
           coeff_Rij=1.0;
        else 
           coeff_Rij=-1.0;

        if(abs(blk_step)==1)
          coeff_deltaij=1.0;
        else if (abs(blk_step)==2) 
          coeff_deltaij=1.5;
        else
          coeff_deltaij=11.0/6.0;

        d_s=1.0/(Chain.N_s);
        NN=Chain.Loa;
        //cout<<"kapa_temp="<<kapa_temp<<endl;  

         array_struct_int_int_int iter;
         iter.i = local_k_grid.indx_1d[K_th].i;
         iter.j = local_k_grid.indx_1d[K_th].j;
         iter.k = local_k_grid.indx_1d[K_th].k;

         dx=Cell.dsize[0];
         dy=Cell.dsize[1];
         dz=Cell.dsize[2];

         if(fabs(dx)<eps) 
             Kx=0.0;
         else  
             Kx=iter.i*2.0*PI/(Grid.n_grid[0]*dx);
         if(fabs(dy)<eps) 
             Ky=0.0;
         else  
             Ky=iter.j*2.0*PI/(Grid.n_grid[1]*dy);
         if(fabs(dz)<eps) 
             Kz=0.0;
         else  
             Kz=iter.k*2.0*PI/(Grid.n_grid[2]*dz);

          
           //cout<<"coeff_deltaij"<<coeff_deltaij<<endl;
           //cout<<"ds "<<ds<<" NN "<<NN<<endl;
           //cout<<"SHF 1x"<<basis_function_SHF[1].l<<endl;


	  for(i=0;i<M_Bar;i++){
	     for(j=0;j<M_Bar;j++){
	         G_R[i][j]=0.0;
// Warning: this is wrong :: G_I[i][j]=coeff_Rij*d_s*(matrix_Rxx[i][j]*Kx +  matrix_Ryy[i][j]*Ky +  matrix_Rzz[i][j]*Kz  );
		 G_I[i][j]=coeff_Rij*d_s*(matrix_Rxx[i][j]*Kx +  matrix_Rxx[i][j]*Ky +  matrix_Ryy[i][j]*Kz  );
		if (i==j) G_R[i][j]=coeff_deltaij + d_s*NN*basis_function_SHF[j].l*(basis_function_SHF[j].l+1)/(2.0*kapa_temp);
				 }
			}

           
//	  for(i=0;i<M_Bar;i++){
//	     for(j=0;j<M_Bar;j++){
//	         G_R_inv[i][j]=G_R[i][j];
//	         G_I_inv[i][j]=G_I[i][j];
//                                 }
//                              }   


         inverse_matrix_complex(M_Bar, G_R, G_I, G_R_inv, G_I_inv);

		//	for(i=0;i<M_Bar;i++){
		//		for(j=0;j<M_Bar;j++){
		//			GAB_R_inverse[K][i][j]=G_R_inverse[i][j];
		//			GAB_I_inverse[K][i][j]=G_I_inverse[i][j];
		//		                    }
		//                            }


///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////	
	free_dmatrix(G_R,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(G_I,0,M_Bar-1,0,M_Bar-1);
	//free_dmatrix(G_R_inverse,0,M_Bar-1,0,M_Bar-1);
	//free_dmatrix(G_I_inverse,0,M_Bar-1,0,M_Bar-1);

}



void compress_G_matrix( grid &Grid, chemical &Chemical, chain &Chain, \
                double ***sa_G_R_inv, double ***sa_G_I_inv, int ***ija_G_R_inv, int ***ija_G_I_inv, cell &Cell)
{
	double **atemp1=dmatrix_NR(1,M_Bar,1,M_Bar);
	double **atemp2=dmatrix_NR(1,M_Bar,1,M_Bar);
	double **G_R_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
	double **G_I_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
         int K;
         int i,j,i_sp,n_sp,blk_step;
         n_sp=Chemical.n_spe;

    for(i_sp=0;i_sp<2*n_sp;i_sp++){
         if(i_sp<n_sp)
           blk_step=3;
         else       
           blk_step=-3;
                        
	for(K=0;K<=LOCAL_SIZE_K-1;K++){
 // store only G matrixs for s>=3. 
    matrix_G_inverse( blk_step, K, Chemical.spe_kapa[i_sp%n_sp], matrix_Rx, matrix_Ry, matrix_Rz, \
                     G_R_inv,G_I_inv,Grid, Chain,Cell ); // Does this pass by reference really work? Check out!!!             
		for(i=1;i<=M_Bar;i++)
			for(j=1;j<=M_Bar;j++){
				atemp1[i][j]=G_R_inv[i-1][j-1]; // note that the starting index changes from 0 to 1
				atemp2[i][j]=G_I_inv[i-1][j-1];
	 	                             }
		
	sprsin(atemp1, M_Bar, Thresh_sprase_matrix, Dim_ordering, sa_G_R_inverse[i_sp][K],ija_G_R_inverse[i_sp][K]);
	sprsin(atemp2, M_Bar, Thresh_sprase_matrix, Dim_ordering, sa_G_I_inverse[i_sp][K],ija_G_I_inverse[i_sp][K]);
                                       }
                                   }

	free_dmatrix_NR(atemp1,1,M_Bar,1,M_Bar);
	free_dmatrix_NR(atemp2,1,M_Bar,1,M_Bar);
	free_dmatrix(G_R_inv,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(G_I_inv,0,M_Bar-1,0,M_Bar-1);

}

