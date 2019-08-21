#include "xy_model.hpp"
#include "rng.hpp"

#include <cassert>




xy_model_grid::xy_model_grid(int Nx, int Ny, int Nz)
	: Nx(Nx), Ny(Ny), Nz(Nz), boundary_cycle(0), s(nullptr)
{
	long N = Nx*Ny*Nz;
	s = new double[N];
	assert(s && "Failed to allocate grid of given size!");
}

xy_model_grid::~xy_model_grid()
{
	delete [] s;
}


xy_model_grid rand_init_grid(int Nx, int Ny, int Nz, xoshiro256 &rng)
{
	xy_model_grid xy_grid(Nx, Ny, Nz);
	long N = Nx*Ny*Nz;
	
	auto u = std::uniform_real_distribution<double>(0, 2*pi);
	for (long i = 0; i < N; ++i) {
		xy_grid.s[i] = u(rng);
	}
	
	return xy_grid;
}

// "fake" the model to be simple ising.
xy_model_grid rand_init_grid_ising(int Nx, int Ny, int Nz, xoshiro256 &rng)
{
	xy_model_grid xy_grid(Nx, Ny, Nz);
	long N = Nx*Ny*Nz;

	double ts = 0.0;
	
	auto u = std::uniform_real_distribution<double>(0, 1);
	for (long i = 0; i < N; ++i) {
		if (u(rng) > 0.5) {
			xy_grid.s[i] = 1;
		} else {
			xy_grid.s[i] = -1;
		}
	}

	return xy_grid;
}


double average_spin(const xy_model_grid &g)
{
	long N = g.Nx*g.Ny*g.Nz;
	double c = 1.0 / static_cast<double>(N);
	double ts = 0.0;
	for (long i = 0; i < N; ++i) {
		ts += g.s[i]*c;
	}

	return ts;
}

double average_plus_spin(const xy_model_grid &g)
{
	long N = g.Nx*g.Ny*g.Nz;
	double c = 1.0 / static_cast<double>(N);
	double tp = 0.0;

	for (long i = 0; i < N; ++i) {
		if (g.s[i] > 0) {
			tp += c;
		}
	}
	return tp;
}


void print_xy_grid_ising(const xy_model_grid &g, std::ostream &out)
{
	long N = g.Nx*g.Ny*g.Nz;

	for (long iy = 0; iy < g.Ny; ++iy) {
		for (long ix = 0; ix < g.Nx; ++ix) {
			long i = ix + g.Nx*iy;
			if (g.s[i] > 0) {
				out << "+";
			} else {
				out << " ";
			}
		}
		out << "\n";
	}
}
