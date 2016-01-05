#ifndef	G_MATRIX_H
#define G_MATRIX_H
#include "cell.h"
void matrix_G_inverse( int blk_step, int K_th, double kapa_temp, double **matrix_Rxx, \ 
     double **matrix_Ryy,double **matrix_Rzz, double **G_R_inv, \
     double **G_I_inv, grid &Grid, chain &Chain, cell &Cell ); 
void compress_G_matrix( grid &Grid, chemical &Chemical, chain &Chain, double ***sa_G_R_inv, \ 
                      double ***sa_G_I_inv, int ***ija_G_R_inv, int ***ija_G_I_inv, cell &Cell);
#endif
