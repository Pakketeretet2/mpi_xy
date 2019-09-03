#include <lyra/lyra.hpp>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <memory>

#include "my_timer.hpp"
#include "rng.hpp"
#include "xy_model.hpp"
#include "output_png.hpp"


struct cmd_opts
{
	cmd_opts() : Nx(32), Ny(32), Nz(1),
	             n_tiles_x(3), n_tiles_y(3), n_tiles_z(1),
	             seed(2019){}
	int Nx, Ny, Nz;
	int n_tiles_x, n_tiles_y, n_tiles_z;
	uint_fast64_t seed;
};



// A little RAII around MPI:
struct mpi_main
{
	mpi_main(int argc, char **argv) : my_rank(-1), comm_size(-1)
	{
		if (USE_MPI) {
			MPI_Init(&argc, &argv);
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
			MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
		} else {
			my_rank = 0;
			comm_size = 1;
		}
		
		log.open("proc." + std::to_string(my_rank) + ".log");

	}
	
	~mpi_main()
	{
		if (USE_MPI) {
			MPI_Finalize();
		}
	}

	int run(int argc, char **argv);

	// MPI-based member variables:
	std::ofstream log;

	int my_rank, comm_size;
};



void create_grid_layout(std::vector<int> &rank_layout, int my_neighbors[6],
                        int my_rank, int n_tiles_x, int n_tiles_y,
                        int n_tiles_z, std::ostream &log)
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
	if (neigh_y_tiles[xy_model_grid::BOTTOM] >= n_tiles_y) {
		neigh_y_tiles[xy_model_grid::BOTTOM] -= n_tiles_y;
	}

	neigh_x_tiles[xy_model_grid::RIGHT] = my_x_tile + 1;
	if (neigh_x_tiles[xy_model_grid::RIGHT] >= n_tiles_x) {
		neigh_x_tiles[xy_model_grid::RIGHT] -= n_tiles_x;
	}
	
        neigh_x_tiles[xy_model_grid::LEFT] = my_x_tile - 1;
	if (neigh_x_tiles[xy_model_grid::LEFT] < 0) {
		neigh_x_tiles[xy_model_grid::LEFT] += n_tiles_x;
	}
	log << "My neigh tiles: (^, >, v, <): \n  ("
	    << neigh_x_tiles[xy_model_grid::TOP] << ", "
	    << neigh_y_tiles[xy_model_grid::TOP] << ")\n  ("
	    << neigh_x_tiles[xy_model_grid::RIGHT] << ", "
	    << neigh_y_tiles[xy_model_grid::RIGHT] << ")\n  ("
	    << neigh_x_tiles[xy_model_grid::BOTTOM] << ", "
	    << neigh_y_tiles[xy_model_grid::BOTTOM] << ")\n  ("
	    << neigh_x_tiles[xy_model_grid::LEFT] << ", "
	    << neigh_y_tiles[xy_model_grid::LEFT] << ")\n";	

	// Leave the last two for other time.
	my_neighbors[0] = neigh_x_tiles[0] + neigh_y_tiles[0]*n_tiles_x;
	my_neighbors[1] = neigh_x_tiles[1] + neigh_y_tiles[1]*n_tiles_x;
	my_neighbors[2] = neigh_x_tiles[2] + neigh_y_tiles[2]*n_tiles_x;
	my_neighbors[3] = neigh_x_tiles[3] + neigh_y_tiles[3]*n_tiles_x;
	my_neighbors[4] = neigh_x_tiles[4] + neigh_y_tiles[4]*n_tiles_x;
	my_neighbors[5] = neigh_x_tiles[5] + neigh_y_tiles[5]*n_tiles_x;

	log << "my neighbors: (top, bottom, right, left, front, back): ("
	    << my_neighbors[0] << ", " << my_neighbors[1] << ", "
	    << my_neighbors[2] << ", " << my_neighbors[3] << ", "
	    << my_neighbors[4] << ", " << my_neighbors[5] << "\n";
}


// Reduces all grids into one and returns it. It throws out the boundary
// layers that belong to other processors, i.e., there are no duplicate
// rows or columns in the resulting grid.
xy_model_grid all_reduce_grid(const xy_model_grid &g,
                              int my_rank, int n_tiles_x, int n_tiles_y,
                              std::ostream &log)
{
	int my_x_tile = my_rank % n_tiles_x;
	int my_y_tile = my_rank / n_tiles_x;

	long Nx_all = n_tiles_x * (g.Nx-2);
	long Ny_all = n_tiles_y * (g.Ny-2);
	
	long n_all = Nx_all * Ny_all;

	xy_model_grid all_grid(Nx_all, Ny_all, 1);
		
	std::vector<double> all_thetas(n_all, 0);
	std::vector<double> gather_thetas(n_all, 0);
	
	long x_offset = my_x_tile * (g.Nx-2);
	long y_offset = my_y_tile * (g.Ny-2);

	log << "Writing grid of proc " << my_rank << " to tile ("
	    << my_x_tile << ", " << my_y_tile << "), global offsets ("
	    << x_offset << ", " << y_offset << ").\n";

	// Suppose we have the following case:
	// 2x2, inner grid of 2x2. Then, we want to map the following
	//
	// [  0,  1,  2,  3  ] [ 16, 17, 18, 19 ]
	// [  4,  5,  6,  7  ] [ 20, 21, 22, 23 ]
	// [  8,  9, 10, 11  ] [ 24, 25, 26, 27 ]
	// [ 12, 13, 14, 15  ] [ 28, 29, 30, 31 ]
	//
	// [ 32, 33, 34, 35  ] [ 48, 49, 50, 51 ]
	// [ 36, 37, 38, 39  ] [ 52, 53, 54, 55 ]
	// [ 40, 41, 42, 43  ] [ 56, 57, 58, 59 ]
	// [ 44, 45, 46, 47  ] [ 60, 61, 62, 63 ]
	//
	// to
	// [  5   6 21 22 ]
	// [  9  10 25 26 ]
	// [  37 38 53 54 ]
	// [  41 42 57 58 ]
	// indexed as
	// [  0  1  2  3 ]
	// [  4  5  6  7 ]
	// [  8  9 10 11 ]
	// [ 12 13 14 15 ]
	//

	for (long iy = 1; iy < g.Ny-1; ++iy) {
		for (long ix = 1; ix < g.Nx-1; ++ix) {
			long i = ix + g.Nx*iy; // index in my local theta vector

			// Construct corresponding index in global theta vector:
			int idx_x = x_offset + ix-1;
			int idx_y = y_offset + iy-1;
			// (idx_x, idy_x) is double index in the global vector.
			// The stride in y is Nx_all
			int idx = idx_x + idx_y * Nx_all;
			
			all_thetas[idx] = g.s[i];

		}
	}

	
	MPI_Allreduce(all_thetas.data(), gather_thetas.data(),
	              n_all, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// Global theta index:
	for (long i = 0; i < n_all; ++i) {
		all_grid.s[i] = gather_thetas[i];
	}

	return all_grid;
}




void print_total_grid(const xy_model_grid &g, int my_rank,
                      int n_tiles_x, int n_tiles_y, std::ostream &out,
                      std::ostream &log, bool print_pad = false)
{

	xy_model_grid all_grid = all_reduce_grid(g, my_rank,
	                                         n_tiles_x, n_tiles_y, log);
	
	double total_theta = 0.0;
	double total_plus  = 0.0;
	double N = 0.0;

	// If we do not print padding, the loop is more tricky.
	// If you do not print pad, then....
	// every (i%g.Nx)th column should be skipped.
	// Also, every (

	std::string y_line = std::string(g.Nx, '-');
	for (long i = 1; i < n_tiles_x; ++i) {
		y_line += "0";
		y_line += std::string(g.Nx,'-');
	}
	
	long Nx_all = n_tiles_x * g.Nx;
	long Ny_all = n_tiles_y * g.Ny;
	long n_all = Nx_all * Ny_all;

	if (my_rank == 0) {
		long n_total_x = n_tiles_x*g.Nx;

		long ix = 0, iy = 0;
		std::ofstream grid_coords("grid_coords_test.dat");
		for (long i = 0; i < n_all; ++i) {
			// Just convert i back to an ix and an iy:
			// i = idx_x + idx_y * g.Ny*n_tiles_y;
			if (all_grid.s[i] > 0) {
				total_theta += 1.0;
				total_plus += 1.0;
				out << "+";
			} else {
				total_theta -= 1.0;
				out << " ";
			}

			++ix;
		        
			if (ix == n_total_x) {
				ix = 0;
				++iy;
				out << "\n";
			} else if (ix % g.Nx == 0) {
				out << "|";
			}

			if ( (ix % n_total_x == 0) &&
			     (iy % g.Ny == 0) ) {
				out << y_line << "\n";
			}
			
			N += 1.0;

			
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



template <typename T>
std::string to_string_w_width(const T &t, std::size_t width, char fill)
{
	std::string res = std::to_string(t);
	std::size_t added = 0;
	std::string front;
	while (added + res.size() < width) {
		front += fill;
		++added;
	}
	return front + res;
}



int mpi_main::run(int argc, char **argv)
{
	cmd_opts opts;
	bool user_help = false;

	
	long tot_sweeps = 10000;
	double T = 2.0;
	// Parse command line args:
	auto cli = lyra::opt(opts.seed, "seed")["-s"]["--seed"]
		("Random seed to use")
		| lyra::opt(opts.Nx, "Nx")["-Nx"]["--Nx"]
		("Size of grid in x-direction.")
		| lyra::opt(opts.Ny, "Ny")["-Ny"]["--Ny"]
		("Size of grid in y-direction.")
		| lyra::opt(opts.Nz, "Nz")["-Nz"]["--Nz"]
		("Size of grid in z-direction.")
		| lyra::opt(opts.n_tiles_x, "n_tiles_x")["-tx"]["--ntiles-x"]
		("Number of processes in x direction.")
		| lyra::opt(opts.n_tiles_y, "n_tiles_y")["-ty"]["--ntiles-y"]
		("Number of processes in y direction.")
		| lyra::opt(opts.n_tiles_z, "n_tiles_z")["-tz"]["--ntiles-z"]
		("Number of processes in z direction.")
		| lyra::opt(T, "T")["-T"]["--temperature"]
		("The temperature of the simulation.")
		| lyra::opt(tot_sweeps, "tot_sweeps")["-n"]["--nsweeps"]
		("Number of Monte Carlo sweeps to perform.")
		| lyra::help(user_help);

	if (!cli.parse({argc, argv})) {
		if (my_rank == 0) {
			std::cerr << "Error parsing options!";
		}
		if (USE_MPI) {
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		return -1;
	}

	int n_tiles_x = opts.n_tiles_x;
	int n_tiles_y = opts.n_tiles_y;
	int n_tiles_z = opts.n_tiles_z;

	
	if (opts.Nz == 1) {

		int n_tiles_all = n_tiles_x * n_tiles_y;
		if (n_tiles_all != comm_size) {
			if (my_rank == 0) {
				std::cerr << "Tiling " << n_tiles_x << "x"
				          << n_tiles_y
				          << " inconsistent with comm size "
				          << comm_size << "!\n";
			}
			return -3;
		}
		
		
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
		if (USE_MPI) {
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		return -2;
	}
	
	
	int my_seed = opts.seed + my_rank;
	xoshiro256 rng(my_seed);
	warm_up(rng);
	xy spins(opts.Nx, opts.Ny, opts.Nz, rng, log);
	
	// This defines the layout:
	int n_tiles = n_tiles_x*n_tiles_y;
	std::vector<int> rank_layout(n_tiles);
	int my_neighbors[6];

	create_grid_layout(rank_layout, my_neighbors, my_rank,
	                   n_tiles_x, n_tiles_y, n_tiles_z, log);
	spins.set_my_neighbors(my_neighbors);
	// spins.rand_init_grid_ising();
	xy_model_grid grid(opts.Nx, opts.Ny, 1);
	for (int i = 0; i < opts.Ny; ++i) {
		for (int j = 0; j < opts.Nx; ++j) {
			int idx = j + i*opts.Nx;
			if (j <= i) {
				grid.s[idx] = 1;
			} else {
				grid.s[idx] = -1;
			}
		}
	}

	spins.set_grid(grid);
	std::string init_grid_name_png = "init_grid_" + std::to_string(my_rank);
	init_grid_name_png += ".png";
	xy_grid_to_png(spins.grid(), init_grid_name_png.c_str());

	std::string init_all_grid_name_png = "init_grid_all.png";
	xy_model_grid total_grid_init = all_reduce_grid(spins.grid(), my_rank,
	                                                n_tiles_x, n_tiles_y,
	                                                log);
	if (my_rank == 0) {
		xy_grid_to_png(total_grid_init, init_all_grid_name_png.c_str());
	}
	
	
	spins.set_temp(T);
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
	
	// Exchange all boundaries first:
	for (int i = 0; i < 4; ++i) {
		spins.exchange_boundaries();
	}

	if (my_rank == 0) {
		std::cerr << "Waiting for all MPI procs to reach pre-run.\n";
	}
	if (USE_MPI) MPI_Barrier(MPI_COMM_WORLD);
	
	if (my_rank == 0) {
		std::cerr << "\n\n****  Starting run!  ****\n";
	}
	my_timer timer;
	if (my_rank == 0) {
		timer.enable_output(std::cerr);
	}
	
	double accepts = 0.0;
	double total = 0.0;
	long long steps_taken = 0;
	
	long tot_moves = 10*grid.Nx*grid.Ny;
	long n_steps_to_exchange = grid.Nx*grid.Ny;

	std::vector<std::string> boundary_names = {"TOP", "BOTTOM", "RIGHT", "LEFT"};
	double total_U = spins.total_spin_energy();
	if (my_rank == 0) {
		std::cerr << "Total energy = " << total_U << "\n";
	}

	double c = 1.0 / static_cast<double>(opts.Nx*opts.Ny);
	
	int png_grids_written = 0;
	std::string out_name = "grids_out/total_grid_";
	out_name += to_string_w_width(png_grids_written, 6, '0');
	out_name += ".png";
	xy_model_grid total_grid = all_reduce_grid(spins.grid(), my_rank,
	                                           n_tiles_x, n_tiles_y,
	                                           log);
	if (my_rank == 0) {
		xy_grid_to_png(total_grid, out_name.c_str());
	}
	++png_grids_written;

	int sweeps_per_write = 50;

	long total_steps = tot_sweeps * tot_moves;
	long n_out = 10000000;
	long f_out = n_out / 50;
	for (long n_sweeps = 0; n_sweeps < tot_sweeps; ++n_sweeps) {
		for (long n_moves = 0; n_moves < tot_moves; ++n_moves) {
			accepts += spins.metropolis_move();
			++total;
			++steps_taken;

			if (n_moves % n_steps_to_exchange == 0) {
				spins.exchange_boundaries();
			}

			if (my_rank == 0 && (steps_taken % n_out == 0) ) {
				double prog = steps_taken;
				prog /= total_steps;
				std::cerr << steps_taken << " / " << total_steps
				          << " (" << prog*100 << " %)\n";
			}
			if (my_rank == 0 && (steps_taken % f_out == 0) ) {
				std::cout << steps_taken << " " << accepts/total << " "
				          << total_U*c << " " << spins.average_up_spin()
				          << " " << spins.average_spin() << "\n";
			}
		
		}
		total_U = spins.total_spin_energy();
		
	        

		if (n_sweeps % sweeps_per_write == 0) {
			total_grid = all_reduce_grid(spins.grid(), my_rank,
			                             n_tiles_x, n_tiles_y,
			                             log);
		
			out_name = "grids_out/total_grid_";
			out_name += to_string_w_width(png_grids_written, 6, '0');
			out_name += ".png";
			
			if (my_rank == 0) {
				xy_grid_to_png(total_grid, out_name.c_str());
			}
			++png_grids_written;
		}
	}
	timer.toc("MPI ising model");
	
	
	return 0;
}




int main( int argc, char **argv )
{
	mpi_main m(argc, argv);

	return m.run(argc, argv);
}

