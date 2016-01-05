#include "matrix.h"
#include "Matrix2.h"
#include "NumMeth.h"

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	
	vv=dvector_NR(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
			if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
			vv[i]=1.0/big;
	}
	
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			
			*d = -(*d);
			
			vv[imax]=vv[j];
			
		}
		
		indx[j]=imax;
		
		if (a[j][j] == 0.0) a[j][j]=TINY;
		
		if (j != n) {
			
			dum=1.0/(a[j][j]);
			
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
			
		}
		
	}
	
	free_dvector_NR(vv,1,n);
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void lubksb(double **a, int n, int *indx, double b[])
{
	
	int i,ii=0,ip,j;
	double sum;
	
	for (i=1;i<=n;i++) {
		
		ip=indx[i];
		
		sum=b[ip];
		
		b[ip]=b[i];
		
		if (ii)
			
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
			
		else if (sum) ii=i;
			
		b[i]=sum;
			
	}
	
	for (i=n;i>=1;i--) {
		
		sum=b[i];
		
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		
		b[i]=sum/a[i][i];
		
	}
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void sprsin(double **a, int n, double thresh, int nmax, double sa[],
			int ija[])//  ordering for the sparse matrix
{
	void nrerror(char error_text[]);
	int i,j;
	int k;
	
	
	
	for (j=1;j<=n;j++) sa[j]=a[j][j];
	ija[1]=n+2;
	k=n+1;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) {
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) nrerror("sprsin: nmax too small");
				sa[k]=a[i][j];
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;
	}
	
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void sprsax(double sa[], int ija[], double x[], double b[],
			int n)// sparse matrix multiplies a vector x
{
	void nrerror(char error_text[]);
	int i,k;
	
	if (ija[1] != n+2) nrerror("sprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++)
			b[i] += sa[k]*x[ija[k]];
		
	}
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void compx_sparse_matrx_mul_vector(int n, double *saR, int *ijaR, double *saI, int *ijaI, double *bR, double *bI, double *cR, double *cI) 
{ // the multiply between a complex sparse matrix and a complex vector, cR for real part, cI for imaginary part
	int i;
	double *bRR=dvector_NR(1,n);  //  Real * Real
	double *bRI=dvector_NR(1,n);  //  Real * Imag
	double *bIR=dvector_NR(1,n);  //  Imag * Real
	double *bII=dvector_NR(1,n);  //  Imag * Imag
	
	double *xR=dvector_NR(1,n);  
	double *xI=dvector_NR(1,n);  
	
	for(i=0;i<n;i++){
		xR[i+1]=bR[i];
		xI[i+1]=bI[i];
	}
	
	
	sprsax(saR, ijaR, xR, bRR, n);
	sprsax(saR, ijaR, xI, bRI, n);
	sprsax(saI, ijaI, xR, bIR, n);
	sprsax(saI, ijaI, xI, bII, n);
	
	
	for(i=0;i<n;i++){
		cR[i]=bRR[i+1] - bII[i+1];
		cI[i]=bRI[i+1] + bIR[i+1];
	}
	
	free_dvector_NR(bRR,1,n);
	free_dvector_NR(bRI,1,n);
	free_dvector_NR(bIR,1,n);
	free_dvector_NR(bII,1,n);
	free_dvector_NR(xR,1,n);
	free_dvector_NR(xI,1,n);
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void compx_matrx_mul_matrx(int n, double **aR, double **aI, double **bR, double **bI, double **cR, double **cI) 
{ // the multiply between two complex matrices, cR for real part, cI for imaginary part
	int i,j,k;
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
		{
			cR[i][j]=0.0;
			cI[i][j]=0.0;
		}
	}
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			for(k=0;k<n;k++){
				cR[i][j]=cR[i][j]+(aR[i][k]*bR[k][j] - aI[i][k]*bI[k][j]);
				cI[i][j]=cI[i][j]+(aR[i][k]*bI[k][j] + aI[i][k]*bR[k][j]);
			}
		}
	}
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void compx_matrx_mul_vector(int n, double **aR, double **aI, double *bR, double *bI, double *cR, double *cI) 
{ // the multiply between a complex matrix and a complex vector, cR for real part, cI for imaginary part
	int i,j;
	
	for(i=0;i<n;i++)
	{
		cR[i]=0.0;
		cI[i]=0.0;
	}
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			cR[i]=cR[i]+(aR[i][j]*bR[j] - aI[i][j]*bI[j]);
			cI[i]=cI[i]+(aR[i][j]*bI[j] + aI[i][j]*bR[j]);
		}
	}
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
double norm_compx_vector_mul_vector(int n, double *aR, double *aI) 
{ // norm of a complex vector
	int i;
	double result=0.0;
	
	for(i=0;i<n;i++){
		result=result+(aR[i]*aR[i] + aI[i]*aI[i]);
	}
	
	return sqrt(result);	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void real_matrx_mul_matrx(int n, double **a, double **b, double **c) 
{ // the multiply between two real matrices
	int i,j,k;
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
		{
			c[i][j]=0.0;
		}
	}
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			for(k=0;k<n;k++){
				c[i][j]=c[i][j]+a[i][k]*b[k][j];
			}
		}
	}
	
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void inverse_matrix_real(int n, double **a_origin, double **a_origin_inverse) ///calculate the inverse of a real matrix
{
	double **a,**y,d,*col;
	int i,j,*indx;
	a=dmatrix_NR(1,n,1,n);
	col=dvector_NR(1,n);
	y=dmatrix_NR(1,n,1,n);
	indx=ivector_NR(1,n);
	
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
		{			
			a[i][j]=a_origin[i-1][j-1];			
		}
		
		ludcmp(a,n,indx,&d);
		for(j=1;j<=n;j++)
		{			
			for(i=1;i<=n;i++) col[i]=0.0;
			col[j]=1.0;
			lubksb(a,n,indx,col);
			for(i=1;i<=n;i++) y[i][j]=col[i];
		}
		
		for(i=1;i<=n;i++)
			for(j=1;j<=n;j++)
			{
				a_origin_inverse[i-1][j-1]=y[i][j];
			}		
			
			free_dmatrix_NR(a,1,n,1,n);
			free_dvector_NR(col,1,n);
			free_dmatrix_NR(y,1,n,1,n);
			free_ivector_NR(indx,1,n);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// convert from double** to Matrix
void Convert(Matrix& Mat, double** DMax, int RowNum, int ColNum)
{
	for(int i = 0; i < RowNum; i++)
	{
		for(int j = 0; j < ColNum; j++)
		{
			Mat(i + 1, j + 1) = DMax[i][j];
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// convert from Matrix to double**
void ConvertInverse(double** DMax, const Matrix& Mat, int RowNum, int ColNum)
{
	for(int i = 0; i < RowNum; i++){
		for(int j = 0; j < ColNum; j++){
			DMax[i][j] = Mat(i + 1, j + 1);
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inverse_matrix_complex(int n, double **a_originR, double **a_originI, double **a_originR_inverse, double **a_originI_inverse) ///calculate the inverse of a complex matrix
{
	Matrix matR(n,n), matI(n,n), matInvR(n,n), matInvI(n,n);
	Convert(matR, a_originR, n, n);
	Convert(matI, a_originI, n, n);
	cinv(matR, matI, matInvR, matInvI);
	ConvertInverse(a_originR_inverse, matInvR, n, n);
	ConvertInverse(a_originI_inverse, matInvI, n, n);
}
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void tridag(double a[], double b[], double c[], double r[], double u[], unsigned long n)
{

	unsigned long j;

	double bet,*gam;

	gam=dvector_NR(1,n);

	if (b[1] == 0.0) nrerror("Error 1 in tridag");

	u[1]=r[1]/(bet=b[1]);

	for (j=1;j<=n;j++) {

		gam[j]=c[j-1]/bet;

		bet=b[j]-a[j]*gam[j];

		if (bet == 0.0)	nrerror("Error 2 in tridag");

		u[j]=(r[j]-a[j]*u[j-1])/bet;

	}

	for (j=(n-1);j>=1;j--){
		
		u[j] -= gam[j+1]*u[j+1];
	}

	free_dvector_NR(gam,1,n);

}
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//////////////////       Cubic Spline Interpolation  /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void spline(double x[],double y[],int n,double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;
	
	u=dvector_NR(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector_NR(u,1,n-1);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;
	
	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double fact (int l) // factorial
{
	int i;
	double result=1.0;
	if (l==0 || l==1)
		result=1.0;
	else
		for(i=2;i<=l;i++)
			result=result*i;

	return result;
}
//////////////////////////////////////////////////////
double fact_con (int m1, int m2) // factorial continue: calculate the factorial from m1*(m1+1)*(m1+2)...*(m2)
{
	int i;
	double result=1.0;
	for(i=m1;i<=m2;i++)
	{
		result=result*i;
	}
		
	return result;
}
//////////////////////////////////////////////////////
int maximum (int a1, int a2, int a3)
{
	int tempi,tempj;
	tempi=(a1>a2)? a1:a2;
	tempj=(tempi>a3)? tempi:a3;
	return tempj;

}
//////////////////////////////////////////////////////

int minimum (int a1, int a2, int a3)
{
	int tempi,tempj;
	tempi=(a1<a2)? a1:a2;
	tempj=(tempi<a3)? tempi:a3;
	return tempj;
	
}
////////////////////////////////////////////////////////
int sign (int l)
{
	int result;
	if ( int(fabs(l)) %2==0)
		result=1;
	else
		result=-1;

	return result;

}
////////////////////////////////////////////////////////
double three_j (int l1, int l2, int l3, int m1, int m2, int m3) // 3j symbol
{
	int FACTVALUE=20;  // if the m>FACTVALUE, factorial m! will be calculate according to (FACTVALUE)!*(FACTVALUE+1)*...*m
	int k,k_min,k_max;
	double temp,tempk;
	double temp1,temp2; // denominator
	double sum;
	double result;
	
	if (m1+m2+m3 !=0 || l1+l2-l3<0 || l1-l2+l3<0 || -l1+l2+l3<0 || fabs(l1)<fabs(m1) || fabs(l2)<fabs(m2) || fabs(l3)<fabs(m3))
	{
		result=0.0;
	}
	else
	{
		if (l1+l2+l3+1>FACTVALUE)
		{
			temp1=1.0/fact(FACTVALUE);
			temp2=1.0/fact_con(FACTVALUE+1,l1+l2+l3+1);
		}
		else
		{
			temp1=1.0/fact(l1+l2+l3+1);
			temp2=1.0;
		}
		
		temp=fact(l1+l2-l3)*fact(l1-l2+l3)*fact(-l1+l2+l3)*temp1;
		temp=temp*fact(l1-m1)*fact(l1+m1)*fact(l2-m2)*temp2*fact(l2+m2)*fact(l3-m3)*fact(l3+m3);
		temp=sign(l1-l2-m3)*sqrt(temp);

		k_min=maximum(0,l2-l3-m1,l1-l3+m2);
		k_max=minimum(l1+l2-l3,l1-m1,l2+m2);
		

		sum=0.0;
		for(k=k_min;k<=k_max;k++)
		{
			tempk=temp/fact(k);
			tempk=tempk/fact(l1+l2-l3-k);
			tempk=tempk/fact(l1-m1-k);
			tempk=tempk/fact(l2+m2-k);
			tempk=tempk/fact(l3-l2+m1+k);
			tempk=tempk/fact(l3-l1-m2+k);
			sum=sum+sign(k)*tempk;
		}
		result=sum;
		
	}
	
	return result;
}
////////////////////////////////////////////////////////
double Triple_product (int l1, int l2, int l3, int m1, int m2, int m3) // triple product of real spherical harmonics functions
{
	double temp;
	double result=0.0;
	double sum;

	if (sign(l1+l2+l3)==-1)
	{
		result=0.0;
	}
	else
	{
		if (m1==0 && m2==0 && m3==0)  // 1
		{
			temp=three_j(l1,l2,l3,0,0,0);
			result=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1))*temp*temp;
		}
		else if (m1>0 && m2>0 && m3>0) // 2
		{
			temp=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1)/2.0);
			temp=temp*three_j(l1,l2,l3,0,0,0);
			sum=sign(m1)*three_j(l1,l2,l3,-m1,m2,m3)+sign(m2)*three_j(l1,l2,l3,m1,-m2,m3)+sign(m3)*three_j(l1,l2,l3,m1,m2,-m3);
			result=temp*sum;
		}
		else if (m1==0 && m2>0 && m3>0) // 3
		{
			temp=sign(m2)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,0,-m2,m3);
		}
		else if (m1>0 && m2==0 && m3>0) // 4
		{
			temp=sign(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,-m1,0,m3);
		}
		else if (m1>0 && m2>0 && m3==0) // 5
		{
			temp=sign(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,-m1,m2,0);
		}
		else if (m1==0 && m2<0 && m3<0) // 6
		{
			temp=sign(m2)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,0,m2,-m3);
		}
		else if (m1<0 && m2==0 && m3<0) // 7
		{
			temp=sign(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,m1,0,-m3);
		}
		else if (m1<0 && m2<0 && m3==0) // 8
		{
			temp=sign(m1)*sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1));
			result=temp*three_j(l1,l2,l3,0,0,0)*three_j(l1,l2,l3,m1,-m2,0);
		}
		else if (m1>0 && m2<0 && m3<0) // 9
		{
			temp=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1)/2.0);
			temp=temp*three_j(l1,l2,l3,0,0,0);
			sum=sign(m2)*three_j(l1,l2,l3,m1,m2,-m3)+sign(m3)*three_j(l1,l2,l3,m1,-m2,m3)-sign(m1)*three_j(l1,l2,l3,m1,m2,m3);
			result=temp*sum;
		}
		else if (m1<0 && m2>0 && m3<0) // 10
		{
			temp=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1)/2.0);
			temp=temp*three_j(l1,l2,l3,0,0,0);
			sum=sign(m1)*three_j(l1,l2,l3,m1,m2,-m3)-sign(m2)*three_j(l1,l2,l3,m1,m2,m3)+sign(m3)*three_j(l1,l2,l3,-m1,m2,m3);
			result=temp*sum;
		}
		else if (m1<0 && m2<0 && m3>0) // 11
		{
			temp=sqrt(1.0*(2*l1+1)*(2*l2+1)*(2*l3+1)/2.0);
			temp=temp*three_j(l1,l2,l3,0,0,0);
			sum=sign(m1)*three_j(l1,l2,l3,m1,-m2,m3)+sign(m2)*three_j(l1,l2,l3,-m1,m2,m3)-sign(m3)*three_j(l1,l2,l3,m1,m2,m3);
			result=temp*sum;
		}
		else
		{
			result=0.0;
		}
	}

	return result;

}
/////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_1D (int n1, double delta1, double *a1)   
{   
	int s1,flag;
	double FUN_A,FUN_B,FUN2,FUN4; 
    static double result_sum;
	FUN_A=a1[0];
	FUN_B=a1[n1];
	FUN2=0.0;
	FUN4=a1[1];
	for (s1=1,flag=1;flag<(n1/2);flag++)
	{
		s1=s1+1;
		FUN2=FUN2+a1[s1];
		s1=s1+1;
		FUN4=FUN4+a1[s1];
	}
	result_sum=(delta1/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2);
	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_1D_NR(int n0, int n1, double delta1, double *a1)   
{  
	
	int s1,flag;
	double FUN_A,FUN_B,FUN2,FUN4;
    double result_sum;
	FUN_A=a1[n0];
	FUN_B=a1[n1];
	FUN2=0.0;
	FUN4=a1[n0+1];
	for (s1=n0+1,flag=1;flag<((n1-n0)/2);flag++)
	{
		s1=s1+1;
		FUN2=FUN2+a1[s1];
		s1=s1+1;
		FUN4=FUN4+a1[s1];
	}
	result_sum=(delta1/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2);
	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_2D (int n1,int n2,double deltalx,double deltaly,double **a2)
{
	int s1,flag,i;
	double FUN_A,FUN_B,FUN2,FUN4; 
	double result_sum;
	double *a1;
	
	a1=new double[n1+1];
	
	for(i=0;i<=n1;i++){
		FUN_A=a2[i][0];
		FUN_B=a2[i][n2];
		FUN2=0.0;
		FUN4=a2[i][1];
		for (s1=1,flag=1;flag<(n2/2);flag++){
			s1=s1+1;
			FUN2=FUN2+a2[i][s1];
			s1=s1+1;
			FUN4=FUN4+a2[i][s1];
		}
		a1[i]=(deltaly/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2); 
	}
	
	FUN_A=a1[0];
	FUN_B=a1[n1];
	FUN2=0.0;
	FUN4=a1[1];
	for (s1=1,flag=1;flag<(n1/2);flag++){
		s1=s1+1;
		FUN2=FUN2+a1[s1];
		s1=s1+1;
		FUN4=FUN4+a1[s1];
	}
	result_sum=(deltalx/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2); 
	
	delete []a1;
	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_2D_mpi (int n1, int n2, double deltalx, double deltaly, double **a2)  
{   // n1, n2 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
	// no volume average.  
	// double a2[n1][n2] --- 0 ~ n1-1 ,  0 ~ n2-1  ;  it assumes the periodic boundary condition here.
	
	int s1,flag,i;
    static double result_sum;
	double FUN_A,FUN_B,FUN2,FUN4; 
    double *a1;
	//a1 = dvector(0,n1-1);
	a1=new double[n1];
	
	for (i=0;i<n1;i++)
	{	
		FUN_A=a2[i][0];
		FUN_B=a2[i][0];
		FUN2=0.0;
		FUN4=a2[i][1];
		for (s1=1,flag=1;flag<(n2/2);flag++)
		{
			s1=s1+1;
			FUN2=FUN2+a2[i][s1];
			s1=s1+1;
			FUN4=FUN4+a2[i][s1];
		}
		a1[i]=(deltaly/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2);  // deltay=delta1
	}
    
	
	FUN2=0.0;
	FUN4=0.0;
	for (s1=0;s1<n1;s1++)
	{
		FUN2=FUN2+a1[s1];
		s1=s1+1;
		FUN4=FUN4+a1[s1];
	}
	result_sum=(deltalx/3.0)*(4.0*FUN4+2.0*FUN2);  // deltax=delta1
	   
	
	//free_dvector(a1,0,n1-1);
	delete []a1;
	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_2D_1D_mpi (int n1, int n2, double deltalx, double deltaly, double *a21)  
{   // it assumes the periodic boundary condition here.
	// calculate the integration through transforming 2D matrix to 1D matrix
	int i,j;
	static double result_sum;
	double **a2;
	a2 = dmatrix(0,n1-1,0,n2-1);

	for(i=0;i<n1;i++)
		for(j=0;j<n2;j++)
			a2[i][j]=a21[i*n2+j];

	result_sum=	simposon_2D_mpi(n1, n2, deltalx, deltaly,a2);

	free_dmatrix(a2,0,n1-1,0,n2-1);
	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double simposon_3D_mpi (int n1, int n2, int n3, double deltalx, double deltaly, double deltalz, double ***a3)
{   // n1, n2, n3 should be even numbers i.e. 2*N (N is an arbitrary integer larger than zero);
	// no volume average.  
	// double a2[n1][n2][n3] --- 0 ~ n1-1 ,  0 ~ n2-1, 0 ~ n3-1  ;  it assumes the periodic boundary condition here.
	
	int s1,flag,i,j;
    static double result_sum;
	double FUN_A,FUN_B,FUN2,FUN4;
    double *a1, **a2;
	a1 = dvector(0,n1-1);
	a2 = dmatrix(0,n1-1,0,n2-1);
	
	for (i=0;i<n1;i++)
		for (j=0;j<n2;j++) 
		{
			FUN_A=a3[i][j][0];
			FUN_B=a3[i][j][0];
			FUN2=0.0;
			FUN4=a3[i][j][1];
			for (s1=1,flag=1;flag<(n3/2);flag++)
			{
				s1=s1+1;
				FUN2=FUN2+a3[i][j][s1];
				s1=s1+1;
				FUN4=FUN4+a3[i][j][s1];
			}
			a2[i][j]=(deltalz/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2);  // deltaz=delta1
		}
		
		for (i=0;i<n1;i++)
		{	
			FUN_A=a2[i][0];
			FUN_B=a2[i][0];
			FUN2=0.0;
			FUN4=a2[i][1];
			for (s1=1,flag=1;flag<(n2/2);flag++)
			{
				s1=s1+1;
				FUN2=FUN2+a2[i][s1];
				s1=s1+1;
				FUN4=FUN4+a2[i][s1];
			}
			a1[i]=(deltaly/3.0)*(FUN_A+FUN_B+4.0*FUN4+2.0*FUN2);    // deltay=delta1
		}
		
		
		FUN2=0.0;
		FUN4=0.0;
		for (s1=0;s1<n1;s1++)
		{
			FUN2=FUN2+a1[s1];
			s1=s1+1;
			FUN4=FUN4+a1[s1];
		}
		result_sum=(deltalx/3.0)*(4.0*FUN4+2.0*FUN2);  // deltax=delta1
		
		free_dvector(a1,0,n1-1);
		free_dmatrix(a2,0,n1-1,0,n2-1);
		return result_sum;
}
/*************************************************************************/
///////////////////////////////////////////////////////////////////////////
double simposon_3D_1D_mpi (int n1, int n2, int n3, double deltalx, double deltaly, double deltalz, double *a31)
{   // it assumes the periodic boundary condition here.
	// calculate the integration through transforming 3D matrix to 1D matrix

	int i,j,k;
	static double result_sum;
	double ***a3=f3tensor(0,n1-1,0,n2-1,0,n3-1);

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			for(k=0;k<n3;k++){
				a3[i][j][k]=a31[i*n2*n3+j*n3+k];
			}
		}
	}

	result_sum=simposon_3D_mpi (n1, n2, n3, deltalx, deltaly, deltalz, a3);

	free_f3tensor(a3,0,n1-1,0,n2-1,0,n3-1);

	return result_sum;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double Sum_2D(int Ni, int Nj, double **a)
{// doing summation for a 2D array
	int i,j;
	double sum_i, sum_j;
	double result;

	sum_i=0.0;
	for(i=0;i<Ni;i++){
		sum_j=0.0;
		for(j=0;j<Nj;j++){
			sum_j=sum_j+a[i][j];
		}
		sum_i=sum_i+sum_j;
	}
	result=sum_i;

	return result;
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
double Sum_sparse_2D(int ndim_i, vector<array_struct_int_int_double> Jij_1D, double *ai, double *aj)
// ndim_i: the dimension for "i" direction
{// doing summation for a 2D array (i,j)
	int i,j;
	int index_nonzero;
	double tempx;
	double result;
	double *temp_i=dvector(0,ndim_i-1);
	
	for(i=0;i<ndim_i;i++){
		temp_i[i]=0.0;
	}
	
	for(index_nonzero=0;index_nonzero<Jij_1D.size();index_nonzero++){
		i=Jij_1D[index_nonzero].i;
		j=Jij_1D[index_nonzero].j;
		tempx=Jij_1D[index_nonzero].x;
		temp_i[i]=temp_i[i] + tempx*aj[j];
	}
	
	result=0.0;
	for(i=0;i<ndim_i;i++){
		result=result + temp_i[i]*ai[i];
	}
	
	free_dvector(temp_i,0,ndim_i-1);
	return result;
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
double double_dot_multi (double **a, double **b) // double dot multiply between two tensors
{
	int i,j;
	double result;
	
	result=0.0;
	
	for(i=0;i<N_dim_ddm;i++)
		for(j=0;j<N_dim_ddm;j++)
		{
			result=result+a[i][j]*b[j][i];
		}
		
		return result;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void inverse_matrix_NR(double **a_origin, double **a_origin_inverse,int N) ///ÇóÄæ¾ØÕó
{
	
	double **a,**y,d,*col;
	int i,j,*indx;
	a=dmatrix_NR(1,N,1,N);
	col=dvector_NR(1,N);
	y=dmatrix_NR(1,N,1,N);
	indx=ivector_NR(1,N);
	
	for(i=1;i<=N;i++)
		for(j=1;j<=N;j++)
		{
			
			a[i][j]=a_origin[i][j];
			
		}
		
		ludcmp(a,N,indx,&d);
		for(j=1;j<=N;j++)
		{
			
			for(i=1;i<=N;i++) col[i]=0.0;
			col[j]=1.0;
			lubksb(a,N,indx,col);
			for(i=1;i<=N;i++) y[i][j]=col[i];
		}
		
		for(i=1;i<=N;i++)
			for(j=1;j<=N;j++)
			{
				a_origin_inverse[i][j]=y[i][j];
			}
			
			
			
			free_dmatrix_NR(a,1,N,1,N);
			free_dvector_NR(col,1,N);
			free_dmatrix_NR(y,1,N,1,N);
			free_ivector_NR(indx,1,N);
}
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void multi_matrix2_NR(double **a, double *b, double *c, int N)
{
	int i,j;
	
	for(i=1;i<=N;i++)
	{
		c[i]=0.0;
	}
	
	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			c[i]=c[i]+a[i][j]*b[j];
		}
	}
	
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void fourn (double data[], unsigned long nn[], int ndim, int isign){//fourn transformation
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;
	
	for(ntot=1,idim=1;idim<=ndim;idim++){
		ntot *= nn[idim];
	}
	nprev=1;
	
	for(idim=ndim;idim>=1;idim--){
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		
		for(i2=1;i2<=ip2;i2+=ip1){
			if(i2 < i2rev){
				for(i1=i2;i1<=i2+ip1-2;i1+=2){
					for(i3=i1;i3<=ip3;i3+=ip2){
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while(ibit >= ip1 && i2rev > ibit){
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		
		ifp1=ip1;
		while(ifp1 < ip2){
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for(i3=1;i3<=ifp1;i3+=ip1){
				for(i1=i3;i1<=i3+ip1-2;i1+=2){
					for(i2=i1;i2<=ip3;i2+=ifp2){
						k1=i2;
						k2=k1+ifp1;
						tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
						tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
/************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////////////
