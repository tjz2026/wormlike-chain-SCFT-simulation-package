#include "global.h"
#include "cell.h"
cell FCC;
cell BCC;
cell GYR;
cell HEX;
cell LAM;

void cell::init_cell( bool Cell_symmetry, bool Cell_fixed, int Cell_type, double lx, double ly, double lz )
{

      dim=global_grid.dim;
      cell_sym=Cell_symmetry;
      cell_fixed=Cell_fixed;
      cell_type=Cell_type; 
      size_change=false;
      size_change_counter=0;
  
      assert(dim >= 1 && dim <= 3); 
// for 1d case, nx=1,ny=1,nz=SIDEz, for 2d case, nx=1,ny=SIDEy,nz=SIDEz.
      if (dim==1) {
          dsize[0]=0.0;
          dsize[1]=0.0;
          dsize[2]=lz/global_grid.n_grid[2];
          vol=global_grid.n_grid[2]*dsize[2];
                  }      
       if (dim==2) {
          dsize[0]=0.0;
          dsize[1]=ly/global_grid.n_grid[1];
          dsize[2]=lz/global_grid.n_grid[2];
          vol=global_grid.n_grid[1]*dsize[1]*global_grid.n_grid[2]*dsize[2];
                  }  
      if (dim==3) {
          dsize[0]=lx/global_grid.n_grid[0];
          dsize[1]=ly/global_grid.n_grid[1];
          dsize[2]=lz/global_grid.n_grid[2];
          vol=global_grid.n_grid[0]*dsize[0]*global_grid.n_grid[1]*dsize[1]*global_grid.n_grid[2]*dsize[2];
                  }  


}





// check whether the cell size matches the group symmetry.

bool cell::check_cell(int field_type)
{
 const double eps=1.0e-5;
 bool ret;
if (dim==1) {
        return true;}
else  if (dim==2) {
    if(abs(field_type)==4) {
       if(fabs((global_grid.n_grid[1]*dsize[1])/(global_grid.n_grid[2]*dsize[2]) - sqrt(3.0))>eps ) 
                 ret=false;
       else 
                 ret= true;}
    else {
       if(fabs((global_grid.n_grid[1]*dsize[1])-(global_grid.n_grid[2]*dsize[2]))>eps) 
                 ret=false;
       else 
                 ret=true;}
    }

else {
    if(abs(field_type)==4) {
       if(fabs((global_grid.n_grid[1]*dsize[1])/(global_grid.n_grid[2]*dsize[2]) - sqrt(3.0))>eps ) 
                 ret=false;
       else  
                 ret=true; }
    else {
       if(fabs((global_grid.n_grid[1]*dsize[1])-(global_grid.n_grid[2]*dsize[2]))>eps || 
          fabs((global_grid.n_grid[0]*dsize[0])-(global_grid.n_grid[1]*dsize[1]))>eps           
              ) 
                 ret=false;
       else 
                 ret=true;
         }
      }
            
 
    return ret;

}

void cell::update_cell( bool Size_change)
{

   size_change=Size_change;
   if(size_change) {
   size_change_counter=size_change_counter+1;
                   }


}











