#include "xy_grid.hpp"
#include "rng.hpp"
#include <cassert>




xy_model_grid::xy_model_grid(int Nx, int Ny, int Nz)
	: Nx(Nx), Ny(Ny), Nz(Nz), s(nullptr)
{
	long N = Nx*Ny*Nz;
	s = new double[N];
	assert(s && "Failed to allocate grid of given size!");
}


xy_model_grid::xy_model_grid()
	: Nx(0), Ny(0), Nz(0), s(nullptr)
{}


xy_model_grid::xy_model_grid(const xy_model_grid &o)
	: Nx(o.Nx), Ny(o.Ny), Nz(o.Nz)
{
	long N = Nx*Ny*Nz;
	s = new double[N];
	assert(s && "Failed to allocate grid of given size in copy!");

	for (long i = 0; i < N; ++i) {
		s[i] = o.s[i];
	}
}


xy_model_grid &xy_model_grid::operator=(const xy_model_grid &o)
{
	xy_model_grid cpy(o);
	using std::swap;
	swap(*this, cpy);
	return *this;
}




xy_model_grid::~xy_model_grid()
{
	delete [] s;
}


void swap(xy_model_grid &f, xy_model_grid &s)
{
	using std::swap;
	swap(f.Nx, s.Nx);
	swap(f.Ny, s.Ny);
	swap(f.Nz, s.Nz);
	swap(f.s, s.s);
	
}
