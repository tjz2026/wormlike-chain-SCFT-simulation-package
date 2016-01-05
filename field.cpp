#include "global.h"
#include "cell.h"
#include "field.h"
// initialize the field with certain specification.
int andsn_dim;
int simple_num;
int iter_counter;
double lambda_andsn; 
double ***andsn_W_sp;
double ***andsn_dW_sp;
double *global_t_diff;
double *yita;  // pressure field 

void field_init(int field_type, grid &Grid, chemical &Chemical, chain &Chain )
{

   int i1,j1,k1;
   int K;
   double temp_coff=0.8;
   double fa;
// init the pressure field
    yita=dvector(0,LOCAL_SIZE-1);

// select the initial field type
switch (field_type) {
// random noise 
   case 0 : {
     break;
         }
// FCC
   case 1 : {
             

  
     break;
         }
// BCC
   case 2 : {
             

  
     break;
         }
// GYR
   case 3 : {
             

  
     break;
         }
// HEX
   case 4 : {
    assert(Grid.dim==2 || Grid.dim==3);
     //  fa=Chain.block_len[0]/(1.0*Chain.N_s);
     //  cout<<"init field: fa = "<<fa<<endl;
       cout<<"Grid.start[0]:"<<Grid.start[0]<<endl;
       cout<<"Grid.start[1]:"<<Grid.start[1]<<endl;
       cout<<"Grid.start[2]:"<<Grid.start[2]<<endl;
       cout<<"Grid.end[0]:"<<Grid.end[0]<<endl;
       cout<<"Grid.end[2]:"<<Grid.end[1]<<endl;
       cout<<"Grid.end[2]:"<<Grid.end[2]<<endl;
       cout<<"Grid.n_grid[0]:"<<Grid.n_grid[0]<<endl;
       cout<<"Grid.n_grid[1]:"<<Grid.n_grid[1]<<endl;
       cout<<"Grid.n_grid[2]:"<<Grid.n_grid[2]<<endl;


// the hexagonal structure is kept translation constant along the X axis, and the ration of size on Y and Z axis is sqrt(3)         
      for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
         for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
            for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {

          //  if(K==0) {
          //   cout<<"fa="<<(Chain.block_len[0]/(1.0*Chain.N_s))<<endl;
          //   cout<<"coscos"<<cos(2.0*PI*(j1+1)/Grid.n_grid[1])*cos(2.0*PI*(k1+1)/Grid.n_grid[2])<<endl;
          //           }


                Chemical.R_sp[0][K]=(Chain.block_len[0]/(1.0*Chain.N_s))*(1.0+ \ 
                 temp_coff*cos(2.0*PI*(j1+1)/global_grid.n_grid[1])*cos(2.0*PI*(k1+1)/global_grid.n_grid[2]));
            if(K==0) {
             cout<<"RA(0) in case 4"<<Chemical.R_sp[0][K]<<endl;
                     }
                                                              }
                                                           }

                                                       }    
     break;
         }
// LAM
   case 5 : {
    assert(Grid.dim==2 || Grid.dim==3);
// the sine-like lam is  along the Z axis
      for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
         for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
            for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {
                Chemical.R_sp[0][K]=(Chain.block_len[0]/(1.0*Chain.N_s))*(1.0+ \ 
                 temp_coff*sin(2.0*PI*(k1+1)/global_grid.n_grid[2]));
                                                              }
                                                           }
                                                       }    

  
     break;
         }

   default : {}
     break;
 } 

// currently, only AB melts is supported, field initialization of more molecular species is to be implemented.

      for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
         for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
            for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {
                assert(K <LOCAL_SIZE);  
            if(K==0) {
             cout<<"RA(0)"<<Chemical.R_sp[0][K]<<endl;
                     }
                Chemical.R_sp[1][K]=1.0-Chemical.R_sp[0][K];
                Chemical.W_sp[0][K]=NXab*Chemical.R_sp[1][K];
                Chemical.W_sp[1][K]=NXab*Chemical.R_sp[0][K];
                                                              }
                                                           }
                                                       }    


}

void andsn_init( int andsn_dim,int num, double lambda, grid &Grid, chemical &Chemical, chain &Chain )
{
    int nch,nb;
    double *f_sp;
    double f_tmp;
    f_sp=dvector(0,Chemical.n_spe-1);

    andsn_dim=andsn_dim;
    cout<<"set andsn_dim="<<andsn_dim<<endl;
    simple_num=num;
    lambda_andsn=lambda; 
    andsn_W_sp=f3tensor(0,Chemical.n_spe-1,0,andsn_dim,0,LOCAL_SIZE-1);
    andsn_dW_sp=f3tensor(0,Chemical.n_spe-1,0,andsn_dim,0,LOCAL_SIZE-1);
    global_t_diff=dvector(0,Chemical.n_spe-1);

           for(nb=0;nb<Chain.n_blk;nb++){
                f_tmp=Chain.block_len[nb]/(1.0*Chain.N_s);
                nch=Chain.block_spe[nb];
                f_sp[nch]=f_sp[nch]+f_tmp;
                                        }
          cout<<"fa,fb= in andsn"<<f_sp[0]<<" "<<f_sp[1]<<endl;
    free_dvector(f_sp,0,Chemical.n_spe-1);
}


void field_clean(chemical &Chemical )
{
    free_dvector(yita,0,LOCAL_SIZE-1);

    free_dvector(global_t_diff,0,Chemical.n_spe-1);
    free_f3tensor(andsn_W_sp,0,Chemical.n_spe-1,0,andsn_dim,0,LOCAL_SIZE-1);
    free_f3tensor(andsn_dW_sp,0,Chemical.n_spe-1,0,andsn_dim,0,LOCAL_SIZE-1);
}


int indx_update_andsn(int index_i, int n_r)
{
	// this function is designed for locating the index of the preceding steps e.g. k-m  or k-n in the notes
	// if n_r=andsn_dim=5
	// 0  --  0  6   12
	// 1  --  1  7   13
	// 2  --  2  8   14
	// 3  --  3  9   15
	// 4  --  4  10  16
	// 5  --  5  11  17
	///////////////////////
	
	int result_i=index_i;
	
	if (index_i<0)
		result_i = index_i + (n_r+1);
	
	return result_i;
}
/************************************************************************/


void andsn_iterate_diblock( int andsn_dim, grid &Grid, chemical &Chemical, chain &Chain, cell &Cell  )  
{

        int K,nch,nb;
	double temp1;
	double temp2;
	double f_tmp;
	double pressure_coeff;
	double *t_diff;
	double  ta_diff,tb_diff,t_diff_global;
	double *f_sp;
        int nx,ny,nz;
        double dx,dy,dz;  

        nx=Grid.n_grid[0];
        ny=Grid.n_grid[1];
        nz=Grid.n_grid[2];

        dx=Cell.dsize[0];  
        dy=Cell.dsize[1];  
        dz=Cell.dsize[2];  

	t_diff=dvector(0,Chemical.n_spe-1);
	f_sp=dvector(0,Chemical.n_spe-1);

           for(nch=0;nch<Chemical.n_spe;nch++){
               f_sp[nch]=0.0;
                                           }
 
           for(nb=0;nb<Chain.n_blk;nb++){
                f_tmp=Chain.block_len[nb]/(1.0*Chain.N_s);
                nch=Chain.block_spe[nb];
                f_sp[nch]=f_sp[nch]+f_tmp;
                                        }
	
	int nr_temp=(iter_counter - simple_num < andsn_dim)? (iter_counter - simple_num):andsn_dim; // get the smaller one

	if (iter_counter<simple_num){
		t_diff[0]=0.0;
		t_diff[1]=0.0;
		for(K=0;K<LOCAL_SIZE;K++){
           //  density_to_field();
// !!!Warning: in future, a more general function to produce the
// new field from density profile according to the saddle point approximation needs to be implemented.
// currently, only AB diblock is supported.
			//yita[K]=0.5*(Chemical.W_sp[0][K]+Chemical.W_sp[1][K] -Chemical.XN[0][1] );
			yita[K]=0.5*(Chemical.W_sp[0][K]+Chemical.W_sp[1][K] );
			temp1=Chemical.XN[0][1]*(Chemical.R_sp[1][K]-f_sp[1])+yita[K] - Chemical.W_sp[0][K];
			temp2=Chemical.XN[1][0]*(Chemical.R_sp[0][K]-f_sp[0])+yita[K] - Chemical.W_sp[1][K];
                        Chemical.W_sp[0][K]=Chemical.W_sp[0][K]+lambda_andsn*(temp1);
                        Chemical.W_sp[1][K]=Chemical.W_sp[1][K]+lambda_andsn*(temp2);
			if(fabs(temp1)>t_diff[0]) t_diff[0]=fabs(temp1);
			if(fabs(temp2)>t_diff[1]) t_diff[1]=fabs(temp2);
		}
                ta_diff=t_diff[0];
                tb_diff=t_diff[0];
                t_diff_global=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&ta_diff,&t_diff_global,1,MPI_DOUBLE, MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&t_diff_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
                global_t_diff[0]=t_diff_global;

                t_diff_global=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&tb_diff,&t_diff_global,1,MPI_DOUBLE, MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&t_diff_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
                global_t_diff[1]=t_diff_global;

	}

	else if (iter_counter == simple_num){
		t_diff[0]=0.0;
		t_diff[1]=0.0;
		for(K=0;K<LOCAL_SIZE;K++){
	       yita[K]=0.5*(Chemical.W_sp[0][K]+Chemical.W_sp[1][K] );
               temp1=Chemical.XN[0][1]*(Chemical.R_sp[1][K]-f_sp[1])+yita[K] - Chemical.W_sp[0][K];
               temp2=Chemical.XN[1][0]*(Chemical.R_sp[0][K]-f_sp[0])+yita[K] - Chemical.W_sp[1][K];
               Chemical.W_sp[0][K]=Chemical.W_sp[0][K]+lambda_andsn*(temp1);
               Chemical.W_sp[1][K]=Chemical.W_sp[1][K]+lambda_andsn*(temp2);

	       yita[K]=0.5*(Chemical.W_sp[0][K]+Chemical.W_sp[1][K] - Chemical.XN[0][1] );
               andsn_W_sp[0][0][K]=Chemical.XN[0][1]*Chemical.R_sp[1][K]+yita[K];
               andsn_W_sp[1][0][K]=Chemical.XN[0][1]*Chemical.R_sp[0][K]+yita[K];
               andsn_dW_sp[0][0][K]=andsn_W_sp[0][0][K]-Chemical.W_sp[0][K];
               andsn_dW_sp[1][0][K]=andsn_W_sp[1][0][K]-Chemical.W_sp[1][K];

	       if(fabs(temp1)>t_diff[0]) t_diff[0]=fabs(temp1);
	       if(fabs(temp2)>t_diff[1]) t_diff[1]=fabs(temp2);
		}

                ta_diff=t_diff[0];
                tb_diff=t_diff[0];
                t_diff_global=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&ta_diff,&t_diff_global,1,MPI_DOUBLE, MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&t_diff_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
                global_t_diff[0]=t_diff_global;

                t_diff_global=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&tb_diff,&t_diff_global,1,MPI_DOUBLE, MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&t_diff_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
                global_t_diff[1]=t_diff_global;

		}

	else {
		double *C_WAB_andsn=dvector_NR(1,nr_temp);
		double *V_WAB_andsn=dvector_NR(1,nr_temp);
		double **U_WAB_andsn=dmatrix_NR(1,nr_temp,1,nr_temp);
		double **inverse_U_WAB_andsn=dmatrix_NR(1,nr_temp,1,nr_temp);
		double *g=dvector(0,LOCAL_SIZE-1);
		double *h=dvector(0,LOCAL_SIZE-1);
		double *gh=dvector(0,LOCAL_SIZE-1);
		double U_temp[2], V_temp[2], U_temp_global[2], V_temp_global[2];
		double U_tmp, V_tmp, U_tmp_global, V_tmp_global;

		double WA_temp, WB_temp, DA_temp, DB_temp;
		int m_andsn, n_andsn;
		int k_andsn=(iter_counter - simple_num)%(andsn_dim+1);

		for(K=0;K<LOCAL_SIZE;K++){
	           yita[K]=0.5*(Chemical.W_sp[0][K]+Chemical.W_sp[1][K] - Chemical.XN[0][1] );
                   andsn_W_sp[0][k_andsn][K]=Chemical.XN[0][1]*Chemical.R_sp[1][K]+yita[K];
                   andsn_W_sp[1][k_andsn][K]=Chemical.XN[0][1]*Chemical.R_sp[0][K]+yita[K];
                   andsn_dW_sp[0][k_andsn][K]=andsn_W_sp[0][k_andsn][K]-Chemical.W_sp[0][K];
                   andsn_dW_sp[1][k_andsn][K]=andsn_W_sp[1][k_andsn][K]-Chemical.W_sp[1][K];
		                         }



                cout<<"Vol cell:"<<Cell.vol<<endl;   
                cout<<"nr_temp"<<nr_temp<<endl;   
                cout<<"andsn_dim"<<andsn_dim<<endl;   
                cout<<"simple num"<<simple_num<<endl;   
                cout<<"iter"<<iter_counter<<endl;   
	
// for AB diblock only, thus nch<2, removed and rewrite for more general purpose in future. 
  for(m_andsn=1;m_andsn<=nr_temp;m_andsn++){
          V_WAB_andsn[m_andsn]=0.0;
      for(n_andsn=m_andsn;n_andsn<=nr_temp;n_andsn++){
          U_WAB_andsn[m_andsn][n_andsn]=0.0;
          for (nch=0;nch<2;nch++) { 
             for(K=0;K<LOCAL_SIZE;K++) {
	g[K]= andsn_dW_sp[nch][k_andsn][K] - andsn_dW_sp[nch][indx_update_andsn(k_andsn-m_andsn,andsn_dim)][K];
	h[K]= andsn_dW_sp[nch][k_andsn][K] - andsn_dW_sp[nch][indx_update_andsn(k_andsn-n_andsn,andsn_dim)][K];
	gh[K]=g[K]*h[K];
				       }
       if(Cell.dim==3) 
	U_tmp=simposon_3D_1D_mpi (nx, ny,nz, dx,dy, dz, gh)/Cell.vol;
       else if(Cell.dim==2) 
	U_tmp=simposon_2D_1D_mpi (ny, nz, dy, dz, gh)/Cell.vol;
       else
	U_tmp=simposon_1D_NR_pbc ( 0, nz, dz, gh)/Cell.vol;

                   U_tmp_global=0.0;
		   MPI_Barrier(MPI_COMM_WORLD);
		   MPI_Reduce(&U_tmp,&U_tmp_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
	           MPI_Bcast(&U_tmp_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		   U_WAB_andsn[m_andsn][n_andsn]=U_WAB_andsn[m_andsn][n_andsn]+U_tmp_global;
                                 }
                                                     }
		
         for (nch=0;nch<2;nch++) { 
	   for(K=0;K<LOCAL_SIZE;K++){
		 g[K]= andsn_dW_sp[nch][k_andsn][K] - andsn_dW_sp[nch][indx_update_andsn(k_andsn-m_andsn,andsn_dim)][K];
		 h[K]= andsn_dW_sp[nch][k_andsn][K];
		 gh[K]=g[K]*h[K];
			            }
       if(Cell.dim==3) 
	V_tmp=simposon_3D_1D_mpi (nx, ny,nz, dx,dy, dz, gh)/Cell.vol;
       else if(Cell.dim==2) 
	V_tmp=simposon_2D_1D_mpi (ny, nz, dy, dz, gh)/Cell.vol;
       else
	V_tmp=simposon_1D_NR_pbc ( 0, nz, dz, gh)/Cell.vol;
                        V_tmp_global=0.0;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&V_tmp,&V_tmp_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Bcast(&V_tmp_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		        V_WAB_andsn[m_andsn]=V_WAB_andsn[m_andsn]+V_tmp_global;
		                 }
                    } // for(m_andsn=1;m_andsn<=nr_temp;m_andsn++)


	for(m_andsn=2;m_andsn<=nr_temp;m_andsn++){
	     for(n_andsn=1;n_andsn<m_andsn;n_andsn++){
				U_WAB_andsn[m_andsn][n_andsn]=U_WAB_andsn[n_andsn][m_andsn];
			                             }
	                                  	 } 

		inverse_matrix_NR(U_WAB_andsn, inverse_U_WAB_andsn, nr_temp);

		multi_matrix2_NR(inverse_U_WAB_andsn, V_WAB_andsn, C_WAB_andsn, nr_temp);

	for(K=0;K<LOCAL_SIZE;K++){
		WA_temp=andsn_W_sp[0][k_andsn][K];
		WB_temp=andsn_W_sp[1][k_andsn][K];
		DA_temp=andsn_dW_sp[0][k_andsn][K];
		DB_temp=andsn_dW_sp[1][k_andsn][K];
		for(n_andsn=1;n_andsn<=nr_temp;n_andsn++){
	WA_temp= WA_temp + C_WAB_andsn[n_andsn]*(andsn_W_sp[0][indx_update_andsn(k_andsn-n_andsn,andsn_dim)][K] - andsn_W_sp[0][k_andsn][K]);
        WB_temp= WB_temp + C_WAB_andsn[n_andsn]*(andsn_W_sp[1][indx_update_andsn(k_andsn-n_andsn,andsn_dim)][K] - andsn_W_sp[1][k_andsn][K]);
	DA_temp= DA_temp + C_WAB_andsn[n_andsn]*(andsn_dW_sp[0][indx_update_andsn(k_andsn-n_andsn,andsn_dim)][K] - andsn_dW_sp[0][k_andsn][K]);
        DB_temp= DB_temp + C_WAB_andsn[n_andsn]*(andsn_dW_sp[1][indx_update_andsn(k_andsn-n_andsn,andsn_dim)][K] - andsn_dW_sp[1][k_andsn][K]);
			}

			Chemical.W_sp[0][K]=WA_temp+lambda_andsn*DA_temp;
			Chemical.W_sp[1][K]=WB_temp+lambda_andsn*DA_temp;
		}


//////////////// convergency criterion ///////
	////////	For A   /////////
		double tempd, tempw;
		double tempd_global, tempw_global;

  for (nch=0;nch<2;nch++) { 
		for(K=0;K<LOCAL_SIZE;K++){
			g[K]= andsn_dW_sp[nch][k_andsn][K];
			h[K]= andsn_dW_sp[nch][k_andsn][K];
			gh[K]=g[K]*h[K];
		                         }
       if(Cell.dim==3) 
	tempd=simposon_3D_1D_mpi (nx, ny,nz, dx,dy, dz, gh)/Cell.vol;
       else if(Cell.dim==2) 
	tempd=simposon_2D_1D_mpi (ny, nz, dy, dz, gh)/Cell.vol;
       else
	tempd=simposon_1D_NR_pbc ( 0, nz, dz, gh)/Cell.vol;
		for(K=0;K<LOCAL_SIZE;K++){
			g[K]= andsn_W_sp[nch][k_andsn][K];
			h[K]= andsn_W_sp[nch][k_andsn][K];
			gh[K]=g[K]*h[K];
		                         }
       if(Cell.dim==3) 
	tempw=simposon_3D_1D_mpi (nx, ny,nz, dx,dy, dz, gh)/Cell.vol;
       else if(Cell.dim==2) 
	tempw=simposon_2D_1D_mpi (ny, nz, dy, dz, gh)/Cell.vol;
       else
	tempw=simposon_1D_NR_pbc ( 0, nz, dz, gh)/Cell.vol;
	        tempd_global=0.0;	
	        tempw_global=0.0;	
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&tempd,&tempd_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&tempw,&tempw_global,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(&tempd_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&tempw_global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		global_t_diff[nch]=sqrt(fabs(tempd_global/tempw_global));
                        }

		free_dvector_NR(C_WAB_andsn,1,nr_temp);
		free_dvector_NR(V_WAB_andsn,1,nr_temp);
		free_dmatrix_NR(U_WAB_andsn,1,nr_temp,1,nr_temp);
		free_dmatrix_NR(inverse_U_WAB_andsn,1,nr_temp,1,nr_temp);
		free_dvector(g,0,LOCAL_SIZE-1);
		free_dvector(h,0,LOCAL_SIZE-1);
		free_dvector(gh,0,LOCAL_SIZE-1);
	        free_dvector(t_diff,0,Chemical.n_spe-1);
	        free_dvector(f_sp,0,Chemical.n_spe-1);
    }

 }

	
