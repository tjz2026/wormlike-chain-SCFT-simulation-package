#ifndef	MDE_LINEAR_H
#define MDE_LINEAR_H

void MDE_linear_driver ( chain *Chain, chemical *Chemical );
void calc_THETAij( chemical &Chemical    );
void Propagator_init( int orient, chain *Chain );
void Pr_stepping(int s_step, int orient, chemical *Chemical, chain *Chain, double **Pr);
void MDE_one_step(int s_step,int orient,grid &Grid,chain *Chain, chemical *Chemical, double **Pr );
#endif
