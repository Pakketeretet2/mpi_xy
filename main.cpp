#include <lyra/lyra.hpp>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <memory>


#include "rng.hpp"
#include "xy_model.hpp"

struct cmd_opts
{
	cmd_opts() : Nx(32), Ny(32), Nz(1), seed(2019){}
	int Nx, Ny, Nz;
	uint_fast64_t seed;
};



// A little RAII around MPI:
struct mpi_main
{
	mpi_main(int argc, char **argv)
	{
		MPI_Init(&argc, &argv);
	}
	
	~mpi_main()
	{
		MPI_Finalize();
	}

	int run(int argc, char **argv);
	
};


void create_grid_layout(std::vector<int> &rank_layout, int my_neighbors[6],
                        int my_rank, int n_tiles_x, int n_tiles_y,
                        int n_tiles_z)
{
	int my_x_tile = my_rank % n_tiles_x;
	int my_y_tile = my_rank / n_tiles_x;

	int rank = 0;
	for (int iy = 0; iy < n_tiles_y; ++iy) {
		for (int ix = 0; ix < n_tiles_x; ++ix) {
			rank_layout[ix + iy*n_tiles_x] = rank;
			++rank;
		}
	}
	
	// top, so who is above me:
	std::vector<int> neigh_x_tiles(6, my_x_tile);
	std::vector<int> neigh_y_tiles(6, my_y_tile);
	
	neigh_y_tiles[xy_model_grid::TOP] = my_y_tile - 1;
	if (neigh_y_tiles[xy_model_grid::TOP] < 0) {
		neigh_y_tiles[xy_model_grid::TOP] += n_tiles_y;
	}

	// bottom, so who is below me:
	neigh_y_tiles[xy_model_grid::BOTTOM] = my_y_tile + 1;
	if (neigh_y_tiles[xy_model_grid::BOTTOM] > n_tiles_y) {
		neigh_y_tiles[xy_model_grid::BOTTOM] -= n_tiles_y;
	}

	neigh_x_tiles[xy_model_grid::RIGHT] = my_x_tile + 1;
	if (neigh_x_tiles[xy_model_grid::RIGHT] > n_tiles_x) {
		neigh_x_tiles[xy_model_grid::RIGHT] -= n_tiles_x;
	}
	
        neigh_x_tiles[xy_model_grid::LEFT] = my_x_tile - 1;
	if (neigh_x_tiles[xy_model_grid::LEFT] < 0) {
		neigh_x_tiles[xy_model_grid::LEFT] += n_tiles_x;
	}

	// Leave the last two for other time.
	my_neighbors[0] = neigh_x_tiles[0] + neigh_y_tiles[0]*n_tiles_x;
	my_neighbors[1] = neigh_x_tiles[1] + neigh_y_tiles[1]*n_tiles_x;
	my_neighbors[2] = neigh_x_tiles[2] + neigh_y_tiles[2]*n_tiles_x;
	my_neighbors[3] = neigh_x_tiles[3] + neigh_y_tiles[3]*n_tiles_x;
	my_neighbors[4] = neigh_x_tiles[4] + neigh_y_tiles[4]*n_tiles_x;
	my_neighbors[5] = neigh_x_tiles[5] + neigh_y_tiles[5]*n_tiles_x;

}


void print_total_grid(const xy_model_grid &g, int my_rank,
                      int n_tiles_x, int n_tiles_y, std::ostream &out)
{
	int my_x_tile = my_rank % n_tiles_x;
	int my_y_tile = my_rank / n_tiles_x;
	long n_all = n_tiles_x * n_tiles_y * g.Nx * g.Ny;
	
	std::vector<double> all_thetas(n_all,0);
	std::vector<double> gather_thetas(n_all,0);

	std::string log_name = "proc." + std::to_string(my_rank) + ".log";
	std::ofstream log(log_name);

	long x_offset = my_x_tile * g.Nx;
	long y_offset = my_y_tile * g.Ny;

	log << "Writing grid of proc " << my_rank << " to tile ("
	    << my_x_tile << ", " << my_y_tile << "), global offsets ("
	    << x_offset << ", " << y_offset << ").\n";

	for (long iy = 0; iy < g.Ny; ++iy) {
		for (long ix = 0; ix < g.Nx; ++ix) {
			long i = ix + g.Nx*iy; // index in my local theta vector

			// Construct corresponding index in global theta vector:
			int idx_x = x_offset + ix;
			int idx_y = y_offset + iy;
			// (idx_x, idy_x) is double index in the global vector.
			// The stride in y is g.Ny*n_tiles_y.
			int idx = idx_x + idx_y * g.Ny*n_tiles_y;
			all_thetas[idx] = g.s[i];
		}
	}

	
	MPI_Allreduce(all_thetas.data(), gather_thetas.data(),
	              n_all, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	double total_theta = 0.0;
	double total_plus  = 0.0;
	double N = 0.0;


	if (my_rank == 0) {
		long n_total_y = n_tiles_y*g.Ny;
		long n_total_x = n_tiles_x*g.Nx;

		std::cerr << "Total grid is " << n_total_x << "x" << n_total_y << "\n";
		for (long i = 0; i < n_all; ++i) {
			if (gather_thetas[i] > 0) {
				total_theta += 1.0;
				total_plus += 1.0;
				out << "+";
			} else {
				total_theta -= 1.0;
				out << " ";
			}
			N += 1.0;

			if (i % n_total_x == n_total_x - 1) {
				out << "\n";
			}
		}
	}
}



void warm_up(xoshiro256 &rng)
{
	std::uniform_real_distribution<double> u(0.0,1.0);
	for (int i = 0; i < 50; ++i) {
		u(rng);
	}
}




int mpi_main::run(int argc, char **argv)
{
	cmd_opts opts;
	bool user_help = false;

	int my_rank = -1;
	int comm_size = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	
	// Parse command line args:
	auto cli = lyra::opt(opts.seed, "seed")["-s"]["--seed"]
		("Random seed to use")
		| lyra::opt(opts.Nx, "Nx")["-Nx"]["--Nx"]
		("Size of grid in x-direction.")
		| lyra::opt(opts.Ny, "Ny")["-Ny"]["--Ny"]
		("Size of grid in y-direction.")
		| lyra::opt(opts.Nz, "Nz")["-Nz"]["--Nz"]
		("Size of grid in z-direction.")
		| lyra::help(user_help);
	
	if (!cli.parse({argc, argv})) {
		if (my_rank == 0) {
			std::cerr << "Error parsing options!";
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
		return -1;
	}

	// Determine the tiling of space:
	int n_tiles_x = 0;
	int n_tiles_y = 0;
	int n_tiles_z = 0;
	
	if (opts.Nz == 1) {

		double sqrt_size = std::sqrt(comm_size);
		n_tiles_x = std::ceil(sqrt_size);
		n_tiles_y = comm_size / n_tiles_x;
		n_tiles_z = 1;

		if (my_rank == 0) {
			std::cerr << "Using the following tiling:\n";
			std::cerr << "tiles in x: " << n_tiles_x << "\n"
			          << "tiles in y: " << n_tiles_y << "\n";
		}

		// Also up-scale Nx and Ny with 2 to have room for the boundaries.
		opts.Nx += 2;
		opts.Ny += 2;
		
	} else {
		if (my_rank == 0) {
			std::cerr << "3D case not implemented yet!\n";
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
		return -2;
	}
	
	
	int my_seed = opts.seed + my_rank;
	xoshiro256 rng(my_seed);
	warm_up(rng);
	xy_model_grid xy_grid = rand_init_grid_ising(opts.Nx, opts.Ny,
	                                             opts.Nz, rng);
	std::cerr << "Grid " << my_rank << ": avg(s) = "
	          << average_spin(xy_grid) << ", n(+) = "
	          << average_plus_spin(xy_grid) << "\n";
	

	// This defines the layout:
	int n_tiles = n_tiles_x*n_tiles_y;
	std::vector<int> rank_layout(n_tiles);
	int my_neighbors[6];

	create_grid_layout(rank_layout, my_neighbors, my_rank,
	                   n_tiles_x, n_tiles_y, n_tiles_z);
	
	
	if (my_rank == 0) {
		std::cerr << "rank layout:\n";
		for (int iy = 0; iy < n_tiles_y; ++iy) {
			for (int ix = 0; ix < n_tiles_x; ++ix) {
				std::cerr << std::setw(8)
				          << rank_layout[ix + iy*n_tiles_x];
			}
			std::cerr << "\n";
		}

	}
	
	std::ofstream total_grid_out;
	if (my_rank == 0) {
		total_grid_out.open("total_grid.dat");
	}
	print_total_grid(xy_grid, my_rank, n_tiles_x, n_tiles_y, total_grid_out);
	

	std::string grid_out = "grid." + std::to_string(my_rank) + ".dat";
	std::ofstream my_out(grid_out);
	print_xy_grid_ising(xy_grid, my_out);
	
	
	return 0;
}




int main( int argc, char **argv )
{
	mpi_main m(argc, argv);

	return m.run(argc, argv);
}

