#include "global.h"
void BF() // constructing the basis function for spherical harmonics
{
	int l,m;
	array_struct_int_int test;

	for(l=0;l<=L_Bar;l++){
		for(m=-l;m<=l;m++){
			test.l=l;
			test.m=m;
			basis_function_SHF.push_back(test);
		}
	}

/*
	int i;
	for(i=0;i<M_Bar;i++){
		printf("%d %d\n",basis_function_SHF[i].l, basis_function_SHF[i].m);
	}
*/
}

void initial( chemical &Chemical)
{
	matrix_Rx=dmatrix(0,M_Bar-1,0,M_Bar-1);
	matrix_Ry=dmatrix(0,M_Bar-1,0,M_Bar-1);
	matrix_Rz=dmatrix(0,M_Bar-1,0,M_Bar-1);
	GAMA=f3tensor(0,M_Bar-1,0,M_Bar-1,0,M_Bar-1);
         
	//q=f3tensor(0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
	//qstar=f3tensor(0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
	//WA=dvector(0,LOCAL_SIZE-1);
	//WB=dvector(0,LOCAL_SIZE-1);
	//WA_out=dmatrix(0,n_r_WAB,0,LOCAL_SIZE-1);
	//WB_out=dmatrix(0,n_r_WAB,0,LOCAL_SIZE-1);
	//dA_anderson=dmatrix(0,n_r_WAB,0,LOCAL_SIZE-1);
	//dB_anderson=dmatrix(0,n_r_WAB,0,LOCAL_SIZE-1);
	
//	M11_out=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	M33_out=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	M12_out=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	M13_out=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	M23_out=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	d11_anderson=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	d33_anderson=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	d12_anderson=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	d13_anderson=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//	d23_anderson=dmatrix(0,n_r_M,0,LOCAL_SIZE-1);
//
//	RA=dvector(0,LOCAL_SIZE-1);
//	RB=dvector(0,LOCAL_SIZE-1);
//	RHOA=dmatrix(0,LOCAL_SIZE-1,0,M_Bar-1);
//	RHOB=dmatrix(0,LOCAL_SIZE-1,0,M_Bar-1);

//	M_OP=f3tensor(0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
//	S_OP=f3tensor(0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
//	SA_OP=f3tensor(0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);
//	SB_OP=f3tensor(0, LOCAL_SIZE-1, 0, N_dim_ddm-1, 0, N_dim_ddm-1);

        
	THETAij=f3tensor(0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	THETAij_M11_M22=dmatrix(0,M_Bar-1,0,M_Bar-1);
	THETAij_M33=dmatrix(0,M_Bar-1,0,M_Bar-1);
	THETAij_M12=dmatrix(0,M_Bar-1,0,M_Bar-1);
	THETAij_M13=dmatrix(0,M_Bar-1,0,M_Bar-1);
	THETAij_M23=dmatrix(0,M_Bar-1,0,M_Bar-1);

	J11ij=dmatrix(0,M_Bar-1,0,M_Bar-1);
	J22ij=dmatrix(0,M_Bar-1,0,M_Bar-1);
	J12ij=dmatrix(0,M_Bar-1,0,M_Bar-1);
	J13ij=dmatrix(0,M_Bar-1,0,M_Bar-1);
	J23ij=dmatrix(0,M_Bar-1,0,M_Bar-1);
             int n_sp;
             //n_sp=AB_melt.n_spe;
             n_sp=Chemical.n_spe;

	    G_R_inverse=dmatrix(0,M_Bar-1,0,M_Bar-1);
	 sa_G_R_inverse=f3tensor(0,2*n_sp-1,0,LOCAL_SIZE_K-1,1,Dim_ordering);
	ija_G_R_inverse=f3tensor_int(0,2*n_sp-1,0,LOCAL_SIZE_K-1,1,Dim_ordering);
	
	    G_I_inverse=dmatrix(0,M_Bar-1,0,M_Bar-1);
	 sa_G_I_inverse=f3tensor(0,2*n_sp-1,0,LOCAL_SIZE_K-1,1,Dim_ordering);
	ija_G_I_inverse=f3tensor_int(0,2*n_sp-1,0,LOCAL_SIZE_K-1,1,Dim_ordering);

    /////////////////////////////////////for M///////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////


	int i,j,k;
	double temp;
	BF();

	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			for(k=0;k<M_Bar;k++){
				GAMA[i][j][k]=Triple_product(basis_function_SHF[i].l, 
					basis_function_SHF[j].l, basis_function_SHF[k].l, basis_function_SHF[i].m, basis_function_SHF[j].m, basis_function_SHF[k].m);
			}
		}
	}

	
	temp=-1.0/sqrt(3.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			matrix_Rx[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 1, basis_function_SHF[i].m, basis_function_SHF[j].m, 1);
		}
	}

	temp=-1.0/sqrt(3.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			matrix_Ry[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 1, basis_function_SHF[i].m, basis_function_SHF[j].m, -1);
		}
	}


       temp=1.0/sqrt(3.0);

	for(i=0;i<M_Bar;i++){

		for(j=0;j<M_Bar;j++){

			matrix_Rz[i][j]=temp*Triple_product(basis_function_SHF[i].l, 

				basis_function_SHF[j].l, 1, basis_function_SHF[i].m, basis_function_SHF[j].m, 0);

		}
	}

	
	temp=1.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			THETAij_M11_M22[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 2);
		}
	}
	
	temp=1.0/sqrt(5.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			THETAij_M33[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 0);
		}
	}
	
	temp=2.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			THETAij_M12[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, -2);
		}
	}
	
	temp=-2.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			THETAij_M13[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 1);
		}
	}
	
	temp=-2.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			THETAij_M23[i][j]=temp*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, -1);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////for J///////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	double temp1, temp2;

	temp1= 1.0/sqrt(15.0);
	temp2=-1.0/(3.0*sqrt(5.0));
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			J11ij[i][j]=temp1*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 2)
				+ temp2*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 0);
		}
	}

	temp1=-1.0/sqrt(15.0);
	temp2=-1.0/(3.0*sqrt(5.0));
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			J22ij[i][j]=temp1*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 2)
				+ temp2*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 0);
		}
	}

	temp1= 1.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			J12ij[i][j]=temp1*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, -2);
		}
	}

	temp1=-1.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			J13ij[i][j]=temp1*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, 1);
		}
	}

	temp1=-1.0/sqrt(15.0);
	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			J23ij[i][j]=temp1*Triple_product(basis_function_SHF[i].l, 
				basis_function_SHF[j].l, 2, basis_function_SHF[i].m, basis_function_SHF[j].m, -1);
		}
	}



//	int n_anderson;
//	int K;
//	
//	for(n_anderson=0;n_anderson<=n_r_WAB;n_anderson++){
//		for(K=0;K<LOCAL_SIZE;K++){
//			WA_out[n_anderson][K]=0.0;
//			WB_out[n_anderson][K]=0.0;
//			dA_anderson[n_anderson][K]=0.0;
//			dB_anderson[n_anderson][K]=0.0;
//		}		
//	}
//
//
//	for(n_anderson=0;n_anderson<=n_r_M;n_anderson++){
//		for(K=0;K<LOCAL_SIZE;K++){
//			M11_out[n_anderson][K]=0.0;
//			M33_out[n_anderson][K]=0.0;
//			M12_out[n_anderson][K]=0.0;
//			M13_out[n_anderson][K]=0.0;
//			M23_out[n_anderson][K]=0.0;
//			d11_anderson[n_anderson][K]=0.0;
//			d33_anderson[n_anderson][K]=0.0;
//			d12_anderson[n_anderson][K]=0.0;
//			d13_anderson[n_anderson][K]=0.0;
//			d23_anderson[n_anderson][K]=0.0;
//		}
//	}
//
//		/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////initialize THETA_nonzero_2D/////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	array_struct_int_int_double test1;
        int K;    
	for(K=0;K<LOCAL_SIZE;K++){
	vector<array_struct_int_int_double> test_THETA;
		for(i=0;i<M_Bar;i++){
			for(j=0;j<M_Bar;j++){
				test1.i=i;
				test1.j=j;
				test1.x=0.0;
				test_THETA.push_back(test1);
			}
		}
		THETA_nonzero_2D.push_back(test_THETA);
	}


	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////input J??ij_nonzero_1D////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			if (fabs(J11ij[i][j]) >= Thresh_sprase_matrix) 
			{
				test1.i=i;
				test1.j=j;
				test1.x=J11ij[i][j];
				J11ij_nonzero_1D.push_back(test1);
			}
		}
	}


	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			if (fabs(J22ij[i][j]) >= Thresh_sprase_matrix) 
			{
				test1.i=i;
				test1.j=j;
				test1.x=J22ij[i][j];
				J22ij_nonzero_1D.push_back(test1);
			}
		}
	}


	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			if (fabs(J12ij[i][j]) >= Thresh_sprase_matrix) 
			{
				test1.i=i;
				test1.j=j;
				test1.x=J12ij[i][j];
				J12ij_nonzero_1D.push_back(test1);
			}
		}
	}


	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			if (fabs(J13ij[i][j]) >= Thresh_sprase_matrix) 
			{
				test1.i=i;
				test1.j=j;
				test1.x=J13ij[i][j];
				J13ij_nonzero_1D.push_back(test1);
			}
		}
	}


	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			if (fabs(J23ij[i][j]) >= Thresh_sprase_matrix) 
			{
				test1.i=i;
				test1.j=j;
				test1.x=J23ij[i][j];
				J23ij_nonzero_1D.push_back(test1);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////input GAMA_nonzero_1D   (ijk) ////////////////
	///////////////////////////////////////////////////////////////////////////////////

	array_struct_int_int_int_double test2;

	for(i=0;i<M_Bar;i++){
		for(j=0;j<M_Bar;j++){
			for(k=0;k<M_Bar;k++){
				if (fabs(GAMA[i][j][k]) >= Thresh_sprase_matrix)
				{
					test2.i=i;
					test2.j=j;
					test2.k=k;
					test2.x=GAMA[i][j][k];
					GAMA_nonzero_1D.push_back(test2);
				}
			}
		}
	}


	printf("GAMA_nonzero_1D.size()=%d\n",GAMA_nonzero_1D.size());

        cout<<"done initialization"<<endl;
	///////////////////////////////////////////////////////////////////////////////////

	

}
/////////////////////////////////////////////////////////////////////////////////////////
