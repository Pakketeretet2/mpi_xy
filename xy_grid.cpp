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





template <typename T>
void binary_write(std::ostream &out, const T &t)
{
	const void *vp = static_cast<const void*>(&t);
	const char *p = static_cast<const char *>(vp);
	out.write(p, sizeof(T));
}


template <typename T>
void binary_write(std::ostream &out, const T* t, std::size_t n_elements)
{
	const void *vp = static_cast<const void*>(t);
	const char *p  = static_cast<const char *>(vp);
	out.write(p, sizeof(T)*n_elements);
}



template <typename T>
void binary_read(std::istream &in, T &t)
{
	void *vp = static_cast<void*>(&t);
	char *p = static_cast<char *>(vp);
	in.read(p, sizeof(T));
}


template <typename T>
void binary_read(std::istream &in, T* t, std::size_t n_elements)
{
	void *vp = static_cast<void*>(t);
	char *p  = static_cast<char *>(vp);
	in.read(p, sizeof(T)*n_elements);
}





// Save grid to binary file:
void xy_model_grid::save_grid(std::ostream &out) const
{
	binary_write(out, Nx);
	binary_write(out, Ny);
	binary_write(out, Nz);

	binary_write(out, s, Nx*Ny*Nz);
}



xy_model_grid::xy_model_grid(std::istream &in)
	: Nx(0), Ny(0), Nz(0), s(nullptr)
{
	binary_read(in, Nx);
	binary_read(in, Ny);
	binary_read(in, Nz);

	s = new double[Nx*Ny*Nz];
	assert(s && "Failed to allocate grid of given size!");

	binary_read(in, s, Nx*Ny*Nz);
		
}
