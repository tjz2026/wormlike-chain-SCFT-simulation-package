/////////////////////////////////////////////////////////////////////////////////////////
void clean()
{
    free_dmatrix(matrix_Rx,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(matrix_Ry,0,M_Bar-1,0,M_Bar-1);
	free_f3tensor(GAMA,0,M_Bar-1,0,M_Bar-1,0,M_Bar-1);
	free_f3tensor(q,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
	free_f3tensor(qstar,0,NMAX,0,LOCAL_SIZE-1,0,M_Bar-1);
	free_dvector(WA,0,LOCAL_SIZE-1);
	free_dvector(WB,0,LOCAL_SIZE-1);
	free_dmatrix(WA_out,0,n_r_WAB,0,LOCAL_SIZE-1);
	free_dmatrix(WB_out,0,n_r_WAB,0,LOCAL_SIZE-1);
	free_dmatrix(dA_anderson,0,n_r_WAB,0,LOCAL_SIZE-1);
	free_dmatrix(dB_anderson,0,n_r_WAB,0,LOCAL_SIZE-1);

	free_dmatrix(M11_out,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(M33_out,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(M12_out,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(M13_out,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(M23_out,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(d11_anderson,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(d33_anderson,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(d12_anderson,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(d13_anderson,0,n_r_M,0,LOCAL_SIZE-1);
	free_dmatrix(d23_anderson,0,n_r_M,0,LOCAL_SIZE-1);

	free_dvector(RA,0,LOCAL_SIZE-1);
	free_dvector(RB,0,LOCAL_SIZE-1);
	free_dmatrix(RHOA,0,LOCAL_SIZE-1,0,M_Bar-1);
	free_dmatrix(RHOB,0,LOCAL_SIZE-1,0,M_Bar-1);


	free_f3tensor(M_OP,0,LOCAL_SIZE-1,0,N_dim_ddm-1,0,N_dim_ddm-1);
	free_f3tensor(S_OP,0,LOCAL_SIZE-1,0,N_dim_ddm-1,0,N_dim_ddm-1);
	free_f3tensor(SA_OP,0,LOCAL_SIZE-1,0,N_dim_ddm-1,0,N_dim_ddm-1);
	free_f3tensor(SB_OP,0,LOCAL_SIZE-1,0,N_dim_ddm-1,0,N_dim_ddm-1);


	free_f3tensor(THETAij,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(THETAij_M11_M22,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(THETAij_M33,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(THETAij_M12,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(THETAij_M13,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(THETAij_M23,0,M_Bar-1,0,M_Bar-1);

	free_dmatrix(J11ij,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(J22ij,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(J12ij,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(J13ij,0,M_Bar-1,0,M_Bar-1);
	free_dmatrix(J23ij,0,M_Bar-1,0,M_Bar-1);

	      free_f3tensor(GA1_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA1_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA1_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA1_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA1_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA1_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA2_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA2_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA2_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA2_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA2_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA2_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA3_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA3_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA3_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA3_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA3_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA3_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB3_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB3_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB3_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	
	      free_f3tensor(GB3_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB3_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB3_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA3_star_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA3_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA3_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GA3_star_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GA3_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GA3_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB1_star_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB1_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB1_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB1_star_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB1_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB1_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB2_star_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB2_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB2_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB2_star_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB2_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB2_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB3_star_R_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB3_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB3_star_R_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

	      free_f3tensor(GB3_star_I_inverse,0,LOCAL_SIZE-1,0,M_Bar-1,0,M_Bar-1);
	 free_dmatrix_NR(sa_GB3_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);
	free_imatrix_NR(ija_GB3_star_I_inverse,1,LOCAL_SIZE,1,M_Bar*M_Bar);

}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
