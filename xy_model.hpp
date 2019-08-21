#ifndef XY_MODEL_HPP
#define XY_MODEL_HPP


#include "rng.hpp"


struct vec_2d {
	double x, y;
};


constexpr const double pi = 3.141592653589793238462643383279502884L;
	

struct xy_model_grid
{

	enum DIRECTIONS { TOP = 0,
	                  BOTTOM,
	                  RIGHT,
	                  LEFT,
	                  FRONT,
	                  BACK
	};
	
	xy_model_grid(int Nx, int Ny, int Nz);
	~xy_model_grid();
	
	int Nx, Ny, Nz;  // Dimensions
	
	
	// The boundary_cycle indicates which domain boundary is currently
	// static. The communication works as follows: At the start of each
	// cycle, each processor communicates the changing boundary to
	// its neighbors. Then each process performs simple MC moves on all
	// but the passive boundary. At the end of the cycle, the active
	// boundaries are communicated and a new passive boundary is chosen.
	// in 2D, boundary cycle goes top (0), right(2), bottom(1), left(3)
	// 
	// in 3D, boundary cycle goes top (0), right face(2), front face(4),
	//                       bottom(1),     left face(3), back face(5)
	//
	// So in a 2D simulation, during cycle 0, on an 8x8 grid, the grid
	// points marked with a 0 can be changed, but the blank ones cannot
	// because their direct neighbors are being changed by another process.
	// 
	// +--------+
	// | 000000 |
	// | 000000 |
	// | 000000 |
	// | 000000 |
	// | 000000 |
	// | 000000 |
	// | 000000 |
	// |        |
	// +--------+
	//
	// Then, in cycle 1, the picture changes to
	//
	// +--------+
	// |        |
	// | 0000000|
	// | 0000000|
	// | 0000000|
	// | 0000000|
	// | 0000000|
	// | 0000000|
	// |        |
	// +--------+	
	//
	int boundary_cycle;
	double *s;          // spins, encoded as a simple angle.

};





xy_model_grid rand_init_grid(int Nx, int Ny, int Nz, xoshiro256 &rng);



xy_model_grid rand_init_grid_ising(int Nx, int Ny, int Nz, xoshiro256 &rng);



void ascii_pain_ising_grid(const xy_model_grid &g);


double average_spin(const xy_model_grid &g);
double average_plus_spin(const xy_model_grid &g);

void print_xy_grid_ising(const xy_model_grid &g, std::ostream &out);




#endif // XY_MODEL_HPP
