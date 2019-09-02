#ifndef XY_GRID_HPP
#define XY_GRID_HPP


#include "rng.hpp"


struct vec_2d {
	double x, y;
};



struct xy_model_grid
{

	enum DIRECTIONS { TOP = 0,
	                  BOTTOM,
	                  RIGHT,
	                  LEFT,
	                  FRONT,
	                  BACK
	};

	xy_model_grid();
	xy_model_grid(int Nx, int Ny, int Nz);
	xy_model_grid(const xy_model_grid &o);
	xy_model_grid &operator=(const xy_model_grid &o);
	
	~xy_model_grid();
	
	int Nx, Ny, Nz;  // Dimensions
	
	double *s;       // spins, encoded as a simple angle.

};



void swap(xy_model_grid &f, xy_model_grid &s);




#endif // XY_GRID_HPP
