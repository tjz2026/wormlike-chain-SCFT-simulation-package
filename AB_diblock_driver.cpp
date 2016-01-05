#include "global.h"
#include "cell.h"
#include "field.h"
#include "MDE_linear.h"
#include "density_energy.h"
#include "AB_diblock_driver.h"
#include "scft_io.h"

 int iter_max;
 double scft_cc;

void AB_diblock_scft(int field_type, int max_iter_num, double cc )
{
 
 bool converged;
 iter_max=max_iter_num;
 scft_cc=cc;
 int simple_mixing_steps=30; 
// init the extern field
 field_init( field_type, local_r_grid, AB_melt, diblock );
 andsn_init(  n_r_WAB,  simple_mixing_steps, lambda_WAB_anderson_const, \
            local_r_grid, AB_melt, diblock );

 converged=true;
 for (iter_counter=1; iter_counter<=iter_max;iter_counter++) {
 MDE_linear_driver ( &diblock, &AB_melt );
 density( &diblock, &AB_melt, &HEX );

  if(myid==0) {
  cout<<"RA(0)"<<AB_melt.R_sp[0][0]<<endl;
  cout<<"RB(0)"<<AB_melt.R_sp[1][0]<<endl;
  cout<<"WA(0)"<<AB_melt.W_sp[0][0]<<endl;
  cout<<"WB(0)"<<AB_melt.W_sp[1][0]<<endl;
              } 

 free_energy_diblock( &AB_melt, &HEX  );

 andsn_iterate_diblock(n_r_WAB, local_r_grid, AB_melt, diblock, HEX  );

 if(myid==0) {
  cout<<"field err: "<<global_t_diff[0]<<"  "<<global_t_diff[1]<<" on"<<iter_counter<<endl;
             } 
 
 converged=converge_or_not();
 if(converged) {
 output1( AB_melt, HEX);
 output( local_r_grid, AB_melt, HEX );
 if(myid==0) cout<<"scft converged,exit"<<endl;
    break; }
                                                             }

}

bool converge_or_not()
{
 bool convg;
if(global_t_diff[0]<scft_cc && global_t_diff[1]<scft_cc) 
   convg=true;
else
  convg=false;
                                                         
 return convg;

}


