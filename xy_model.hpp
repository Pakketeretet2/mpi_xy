#ifndef XY_MODEL_HPP
#define XY_MODEL_HPP


#include "xy_grid.hpp"
#include "rng.hpp"


constexpr const double pi = 3.141592653589793238462643383279502884L;
	

class xy {
public:
	xy(int Nx, int Ny, int Nz, xoshiro256 &random_generator,
	   std::ostream &log);

	void rand_init_grid();
	void rand_init_grid_ising();
	void set_grid(const xy_model_grid &new_grid);

	double average_spin() const;
	double average_up_spin() const;

	void print_xy_grid_ising(std::ostream &out);

	int metropolis_move();

	double spin_energy(int ix, int iy, int iz = 0) const;
	double total_spin_energy() const;

	void set_temp(double T);
	void temp() const;

	
	// Some MPI-related functions:
	void exchange_boundaries();
	void set_my_neighbors(int neighs[6]);

	const xy_model_grid &grid() const
	{ return g; }

	const char *boundary2name[6] = { "TOP", "BOTTOM", "LEFT", "RIGHT",
	                                 "FRONT", "BACK" };
	
private:
	xy_model_grid g;
	xoshiro256 rng;

	// Physics-related:
	double T, beta;

	// Grid-related:
	// These regulate which points are the grid you can change.
	// In MPI runs, 0 and N-1 are usually off-limits because they
	// belong to your neighboring procs. Also, 1 and N-2 are off-
	// limits unless you can guarantee that its neighbors on the other
	// proc are not changed. See also boundary_cycle.
	int xlo, xhi, ylo, yhi, zlo, zhi;

	
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
	// |        |
	// |  0000  |
	// |  0000  |
	// |  0000  |
	// |  0000  |
	// |  0000  |
	// |        |
	// |        |
	// +--------+
	//
	// Then, in cycle 1, the picture changes to
	//
	// +--------+
	// |        |
	// |        |
	// |  00000 |
	// |  00000 |
	// |  00000 |
	// |  00000 |
	// |        |
	// |        |
	// +--------+	
	//
	int boundary_cycle;
	

	// MPI-related:
	int my_rank;
	int comm_size;
	int my_neighbors[6];
	std::ostream &log;
};





#endif // XY_MODEL_HPP
