// General header file for C++ programs 
// in "Numerical Methods for Physics" 

#ifndef NUM_METH_H
#define NUM_METH_H

#include <iostream>
#include <fstream>
#include <assert.h>  
#include <cmath>
#include "Matrix2.h"

void cinv( Matrix RealA, Matrix ImagA, 
			 Matrix& RealAinv, Matrix& ImagAinv );
#endif
