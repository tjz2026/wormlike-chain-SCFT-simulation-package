#ifndef	CELL_H
#define CELL_H

struct cell {
       bool cell_sym;  // whether cell is of a certain symmetry or not, for unit cell calculation, cell_sym==true.
// otherwise cell_sym==false, which means cell size on each axis is user defined.
       bool cell_fixed; // cell size can not be altered if cell_fixed == true.
       int cell_type;
       bool size_change;
       int  size_change_counter;
       int  dim;
       double dsize[3];
       double vol; 
void init_cell( bool Cell_symmetry, bool Cell_fixed, int Cell_type, double lx, double ly, double lz );
bool check_cell(int field_type );
void update_cell( bool Size_change);
};

extern cell FCC;
extern cell BCC;
extern cell HEX;
extern cell GYR;
extern cell LAM;

#endif
