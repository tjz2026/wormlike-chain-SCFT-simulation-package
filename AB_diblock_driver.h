#ifndef	AB_DIBLOCK_DRIVER_H
#define AB_DIBLOCK_DRIVER_H
extern int iter_max;
extern double scft_cc;
void AB_diblock_scft(int field_type, int max_iter_num, double cc);
bool converge_or_not();
#endif
