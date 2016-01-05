#include "global.h"
#include "initial.h"
#include "fft.h"
#include "G_matrix.h"
#include "MDE_linear.h"
# define Maier_Saupe 0
/////////////////////////////////////////////////////////////////////////////////////////////

void MDE_linear_driver ( chain *Chain, chemical *Chemical )

{
 int s;
 const int forward=1;
 const int backward=-1;
 double **Pr=dmatrix(0,LOCAL_SIZE-1,0,M_Bar-1);
  
    Propagator_init( forward,  Chain );
 // cout<<"WA(0)"<<Chemical->W_sp[0][0]<<endl;
 // cout<<"WB(0)"<<Chemical->W_sp[1][0]<<endl;
  for (s=1;s<=Chain->N_s;s++) {
  Pr_stepping(s, forward, Chemical, Chain, Pr);
  MDE_one_step( s,forward, local_r_grid,  Chain, Chemical, Pr );
 // if( myid==0) {
 // cout<<"0q_f[s][0][0] s= "<<s<<" "<<Chain->q_f[s][0][0]<<endl;
 //              } 
  assert(Chain->q_f[s][0][0]<10);

                            }


  Propagator_init( backward, Chain );
  for (s=Chain->N_s-1;s>=0;s--) {
  Pr_stepping(s, backward,  Chemical, Chain, Pr);
  MDE_one_step( s,backward, local_r_grid,  Chain, Chemical, Pr );
//  if(s==Chain->N_s-1 && myid==0 ) {
 //      cout<<"Pr(0,0) s=NS-1"<<Pr[0][0]<<endl;
 //                                 }

 // if( myid==0 ) {
 // cout<<"0q_b[s][0][0] s= "<<s<<" "<<Chain->q_b[s][0][0]<<endl;
 //              } 
                             }

  free_dmatrix(Pr,0,LOCAL_SIZE-1,0,M_Bar-1);
}


void calc_THETAij( chemical &Chemical    )
{
	int i,j;
	int K;
	for(K=0;K<LOCAL_SIZE;K++){
		for(i=0;i<M_Bar;i++){
			for(j=0;j<M_Bar;j++){
		THETAij[K][i][j]=(Chemical.M_OP[K][0][0]-Chemical.M_OP[K][1][1])*THETAij_M11_M22[i][j] 
			+ Chemical.M_OP[K][2][2]*THETAij_M33[i][j]
			+ Chemical.M_OP[K][0][1]*THETAij_M12[i][j] 
			+ Chemical.M_OP[K][0][2]*THETAij_M13[i][j] 
			+ Chemical.M_OP[K][1][2]*THETAij_M23[i][j];
			}
		}
	}

	THETA_nonzero_2D.clear();   // remove the existing elements in "THETA_nonzero_2D"
	

	array_struct_int_int_double test1;
	
	for(K=0;K<LOCAL_SIZE;K++){
		vector<array_struct_int_int_double> test_THETA;
		for(i=0;i<M_Bar;i++){
			for(j=0;j<M_Bar;j++){
				if (fabs(THETAij[K][i][j]) >= Thresh_sprase_matrix)
				{
					test1.i=i;
					test1.j=j;
					test1.x=THETAij[K][i][j];
					test_THETA.push_back(test1);
				}
			}
		}
		THETA_nonzero_2D.push_back(test_THETA);
	}
	
}
/////////////////////////////////////////////////////////////////////////////////////////////

// initial condition of propagator q_f and q_b
void Propagator_init( int orient, chain *Chain )
{
       int K,i;
     if(orient==1) {

	for(K=0;K<LOCAL_SIZE;K++){
		for(i=0;i<M_Bar;i++){
			Chain->q_f[0][K][i]=0.0;
			if (i==0) Chain->q_f[0][K][i]=1.0;
	                            }
                        	 }
                    }

      else if (orient==-1) {
	for(K=0;K<LOCAL_SIZE;K++){
		for(i=0;i<M_Bar;i++){
			Chain->q_b[Chain->N_s][K][i]=0.0;
			if (i==0) Chain->q_b[Chain->N_s][K][i]=1.0;
	                            }
                        	 }

                           }  

}


void Pr_stepping(int s_step, int orient, chemical *Chemical, chain *Chain, double **Pr)

{
     int i,j,l,n,n1,nblk,nch,K,index_nonzero;
     double ds;
     double q_coef1[3][3];
     double q_coef2[3][3];
     double temp1,temp2,Wtemp,tempx;  
// set up coefficient for s_step=1,2,s(s>=3)
         q_coef1[0][0]=1.0; 
         q_coef1[0][1]=0.0; 
         q_coef1[0][2]=0.0; 

         q_coef1[1][0]=2.0; 
         q_coef1[1][1]=-0.5; 
         q_coef1[1][2]=0.0; 

         q_coef1[2][0]=3.0; 
         q_coef1[2][1]=-1.5; 
         q_coef1[2][2]=1.0/3.0; 

         q_coef2[0][0]=1.0; 
         q_coef2[0][1]=0.0; 
         q_coef2[0][2]=0.0; 

         q_coef2[1][0]=2.0; 
         q_coef2[1][1]=-1.0; 
         q_coef2[1][2]=0.0; 

         q_coef2[2][0]=3.0; 
         q_coef2[2][1]=-3.0; 
         q_coef2[2][2]=1.0; 

        ds=1.0/Chain->N_s;
	for(K=0;K<LOCAL_SIZE;K++){
		for(i=0;i<M_Bar;i++){
			Pr[K][i]=0.0;
				    }
			         }

//  if(s_step==31) {
//  cout<<"block num"<<Chain->n_blk<<endl;
//  cout<<"block begin[0]"<<Chain->block_begin[0]<<endl;
//  cout<<"block end[0]"<<Chain->block_end[0]<<endl;
//  cout<<"block begin[1]"<<Chain->block_begin[1]<<endl;
//  cout<<"block end[1]"<<Chain->block_end[1]<<endl;
//                 }
 
  assert(orient == 1 || orient==-1);  // two directions of a chain
  if(orient==1) {
      n = 3 < s_step ? 3 : s_step;
      assert(s_step > 0 );}  // s_step should not be the starting indx of a chain.
  else {
      assert(s_step < Chain->N_s ); 
      n = 3 < (Chain->N_s-s_step) ? 3 : (Chain->N_s-s_step);}

//get block indx
      for(n1=0;n1<Chain->n_blk;n1++){
          if(s_step>Chain->block_begin[n1] && s_step<=Chain->block_end[n1]) {
                nblk=n1;
                break; }
                                }
         if(s_step==0) nblk=0; 

         nch=Chain->block_spe[nblk];  
//check if n, nch and nblk are all correct
 //    if(myid==0 && orient==-1 && s_step==Chain->N_s-1) {
  //       cout<<"s_step="<<s_step<<endl;
   //        cout<<"n="<<n<<endl;
   //        cout<<"nblk="<<nblk<<endl;
   //        cout<<"nch="<<nch<<endl;
   //         }

     MPI_Barrier(MPI_COMM_WORLD);
      
	for(K=0;K<LOCAL_SIZE;K++){
		Wtemp=Chemical->W_sp[nch][K];
# if Maier_Saupe
		for(index_nonzero=0;index_nonzero<THETA_nonzero_2D[K].size();index_nonzero++){
			i=THETA_nonzero_2D[K][index_nonzero].i;
			j=THETA_nonzero_2D[K][index_nonzero].j;
			tempx=THETA_nonzero_2D[K][index_nonzero].x;
                        temp2=0.0;
                     if(orient==1)  
                       for(l=1;l<=n;l++){
                         temp2=temp2+ q_coef2[n-1][l-1]*Chain->q_f[s_step-l][K][j];}
                     else
                       for(l=1;l<=n;l++){
                         temp2=temp2+ q_coef2[n-1][l-1]*Chain->q_b[s_step+l][K][j];}
			//Pr[K][i]=Pr[K][i] + ds*tempx*q[0][K][j];
			Pr[K][i]=Pr[K][i] + ds*tempx*temp2;
		}
# else
# endif
		for(i=0;i<M_Bar;i++){
                        temp1=0.0;
                        temp2=0.0;
                     if(orient==1)  
                       for(l=1;l<=n;l++){
                        assert(l<=3 && n<=3 && s_step-l>=0);
                         temp1=temp1+ q_coef1[n-1][l-1]*Chain->q_f[s_step-l][K][i];
                         temp2=temp2+ q_coef2[n-1][l-1]*Chain->q_f[s_step-l][K][i];}
                     else
                       for(l=1;l<=n;l++){
                         temp1=temp1+ q_coef1[n-1][l-1]*Chain->q_b[s_step+l][K][i];
                         temp2=temp2+ q_coef2[n-1][l-1]*Chain->q_b[s_step+l][K][i];}
  
			//Pr[K][i]=Pr[K][i] + (1.0 - ds*Wtemp)*q[0][K][i];
                        //if(K==0 && s_step==Chain->N_s-1 && myid==0 && i==0) cout<<"temp1,2"<<temp1<<"  "<<temp2<<endl;
                        //if(K==0 && s_step==Chain->N_s-1 && myid==0 && i==0) cout<<"q_b[s_step+l][K][i]"<<Chain->q_b[100][K][i]<<endl;
			Pr[K][i]=Pr[K][i] + temp1 - ds*Wtemp*temp2;
		}
	}
           
}


/////////////////////////////////////////////////////////////////////////////////////////////
void MDE_one_step(int s_step,int orient,grid &Grid,chain *Chain, chemical *Chemical, double **Pr )
{
	int i,j,i1,j1,k1,n,nblk,nch,nG,s,blk_step;
	int i_i,j_j,k_k,K;
        int n1,n2,n3; 

	double **Pk_Real=dmatrix(0,LOCAL_SIZE_K-1,0,M_Bar-1);
	double **Pk_Imag=dmatrix(0,LOCAL_SIZE_K-1,0,M_Bar-1);
	double **q_temp_Real=dmatrix(0,LOCAL_SIZE_K-1,0,M_Bar-1);
	double **q_temp_Imag=dmatrix(0,LOCAL_SIZE_K-1,0,M_Bar-1);
	double **G_R_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
	double **G_I_inv=dmatrix(0,M_Bar-1,0,M_Bar-1);
	////////////////////////////////////////////////////////

         n1=Grid.n_grid[0];
         n2=Grid.n_grid[1];
         n3=Grid.n_grid[2];
//get block indx
      for(n=0;n<Chain->n_blk;n++){
          if(s_step>Chain->block_begin[n] && s_step<=Chain->block_end[n]) {
                nblk=n;
     //    cout<<"nsp="<<nsp<<" s_step= "<<s_step<<endl;
                break; }
                                } 
         if(s_step==0) nblk=0;

         nch=Chain->block_spe[nblk];

         if(orient==1) { nG=n;}
         else {nG=n+ Chemical->n_spe;}
         if(s_step==0) nG=0;

        // cout<<"nG=,nch"<<nG<<" "<<nch<<endl;
         assert(nG >=0 && nG<2*Chemical->n_spe); 
         

 for(i=0;i<M_Bar;i++){
    for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
       for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
// why k<=local_r_grid.end[2]+2 ??, because we assume SIDEz to be even.
// Warning: it is wrong!! Make sure you understand this. 
//for (k1=Grid.start[2];k1<=Grid.end[2]+2;k1++,K++) {
          for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {
            i_i=i1-Grid.start[0];            
            j_j=j1-Grid.start[1];            
            k_k=k1-Grid.start[2];
            fft_pll.r_in[i_i*(n2*(n3+2))+j_j*(n3+2)+k_k]=Pr[K][i]; // the compilcated subindex could be inefficient, try to 
// optimize in future.
                                                            }    
                                                         }
                                                     } 
// this actually ineffcient when placing the FFT as a inner loop, try to use fftw MANYPLAN in futrure 
        fftw_execute(fft_pll.plan_f);  

       for (K=0;K<LOCAL_SIZE_K;K++) {
         Pk_Real[K][i]=fft_pll.c_out[K][0];
         Pk_Imag[K][i]=fft_pll.c_out[K][1]; 
   //  if(K==0) cout<<"Pk_R(0,0)"<<Pk_Real[0][0]<<endl;  
                                    }
                    } 

 if(s_step<=2 && orient==1) {
      //if(s_step==1) blk_step=1; 
      //if(s_step==2) blk_step=2;
      blk_step=s_step; 
  for(K=0;K<LOCAL_SIZE_K;K++)   {
    matrix_G_inverse( blk_step, K, 0.5, matrix_Rx, matrix_Ry, matrix_Rz, \
                G_R_inv,G_I_inv,global_grid, diblock, HEX ); // Does this pass by reference really work? Check out!!!             
    compx_matrx_mul_vector(M_Bar, G_R_inv, G_I_inv, Pk_Real[K], Pk_Imag[K], q_temp_Real[K], q_temp_Imag[K]);

   //  if(K==0) cout<<"q_temp_R(0,0)"<<q_temp_Real[0][0]<<endl;  
                                }             
                            }
 else if(s_step>=Chain->N_s-2 && orient==-1) {
      if(s_step==Chain->N_s-1) blk_step=-1; 
      if(s_step==Chain->N_s-2) blk_step=-2;
  for(K=0;K<LOCAL_SIZE_K;K++){
    matrix_G_inverse( blk_step, K, 0.5, matrix_Rx, matrix_Ry, matrix_Rz, \
                G_R_inv,G_I_inv,global_grid, diblock, HEX ); // Does this pass by reference really work? Check out!!!             
    compx_matrx_mul_vector(M_Bar, G_R_inv, G_I_inv, Pk_Real[K], Pk_Imag[K], q_temp_Real[K], q_temp_Imag[K]);
                             }             
                                 }
 else {
 for(K=0;K<LOCAL_SIZE_K;K++){
	//	compx_matrx_mul_vector(M_Bar, GA1_R_inverse[K],GA1_I_inverse[K],Pk_Real[K],
	//		Pk_Imag[K],q_temp_Real[K],q_temp_Imag[K]);
	compx_sparse_matrx_mul_vector(M_Bar, sa_G_R_inverse[nG][K], ija_G_R_inverse[nG][K], \
                                      sa_G_I_inverse[nG][K], ija_G_I_inverse[nG][K], \
					Pk_Real[K], Pk_Imag[K],q_temp_Real[K],q_temp_Imag[K]);
		  	    }
     }
     

	for(i=0;i<M_Bar;i++){
		for(K=0;K<LOCAL_SIZE_K;K++){
			fft_pll.c_out[K][0]=q_temp_Real[K][i];
			fft_pll.c_out[K][1]=q_temp_Imag[K][i];
	                                   }

        fftw_execute(fft_pll.plan_b);  

    for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
       for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
          for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {
            i_i=i1-Grid.start[0];            
            j_j=j1-Grid.start[1];            
            k_k=k1-Grid.start[2];
// find a better way to deal with orient in future. How to pass Chain.q_f as a pointer??
            switch(orient) 
            {
              case 1 :  
              {
              Chain->q_f[s_step][K][i]=fft_pll.r_in[i_i*(n2*(n3+2))+j_j*(n3+2)+k_k]/global_grid.grid_num;
              break;
              }
              case -1 :  
              {
              Chain->q_b[s_step][K][i]=fft_pll.r_in[i_i*(n2*(n3+2))+j_j*(n3+2)+k_k]/global_grid.grid_num;
              break;
              }
              default: {}
             }
                                                            }
                                                         } 
			                              }
                 } // end loop i for M_Bar
	free_dmatrix(Pk_Real,0,LOCAL_SIZE_K-1,0,M_Bar-1);
	free_dmatrix(Pk_Imag,0,LOCAL_SIZE_K-1,0,M_Bar-1);
	free_dmatrix(q_temp_Real,0,LOCAL_SIZE_K-1,0,M_Bar-1);
	free_dmatrix(q_temp_Imag,0,LOCAL_SIZE_K-1,0,M_Bar-1);
	free_dmatrix(G_R_inv,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(G_I_inv,0,M_Bar-1,0,M_Bar-1);

}












