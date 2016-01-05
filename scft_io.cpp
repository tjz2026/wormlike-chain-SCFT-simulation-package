#include "global.h"
#include "cell.h"
#include "field.h"
/////////////////////////////////////////////////////////////////////////
void output1( chemical &Chemical, cell &Cell)
{

	FILE *pf;
	if (myid == 0){
		if (pf=fopen("results.dat","a")) {
             fprintf(pf,"dx*SIDEx=%lf  dy*SIDEy=%lf  dz*SIDEz=%lf  FE=%lf  pff=%lf  %e %e\n",\
             Cell.dsize[0]*SIDEx,Cell.dsize[1]*SIDEy,Cell.dsize[2]*SIDEz,FE_global,pff_global,global_t_diff[0],global_t_diff[1]);
		    fclose(pf);	

		                                 }
	             }

}

void output( grid &Grid, chemical &Chemical, cell &Cell )
{
char name1[80],name2[80],name3[80];


	register int i;
        int K,K_i, K_j, K_k,i1,j1,k1;
        double dx,dy,dz;
	FILE *pf;
        
        dx=Cell.dsize[0];
        dy=Cell.dsize[1];
        dz=Cell.dsize[2];
///////////////////////////////////////////////////////////////////////////////
	if (myid == 0){
		if (pf=fopen("free_energy.dat","a"))
		{	

             fprintf(pf,"dx*SIDEx=%lf  dy*SIDEy=%lf  dz*SIDEz=%lf  FE=%lf  pff=%lf  %e %e\n",\
             Cell.dsize[0]*SIDEx,Cell.dsize[1]*SIDEy,Cell.dsize[2]*SIDEz,FE_global,pff_global,global_t_diff[0],global_t_diff[1]);

			fclose(pf);	
		}
	}

	sprintf(name1,"RHO_total%d.dat",myid);
	if (pf=fopen(name1,"wb"))
	{
    for (i1=Grid.start[0],K=0;i1<=Grid.end[0];i1++) {
       for (j1=Grid.start[1];j1<=Grid.end[1];j1++)      {
          for (k1=Grid.start[2];k1<=Grid.end[2];k1++,K++) {
	fprintf(pf,"%lf %lf %lf %lf %lf %lf\n",i1*dx, j1*dy, k1*dz, Chemical.R_sp[0][K], \
               Chemical.R_sp[1][K], (1.0-Chemical.R_sp[0][K]-Chemical.R_sp[1][K]));
				                          }
			                                }
		                                    }
		fclose(pf);	
	}

/////////////////////////////////////////////////////////////////////


}

/************************************************************************/
