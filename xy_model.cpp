#include "xy_grid.hpp"
#include "xy_model.hpp"
#include "rng.hpp"

#include <cassert>
#include <mpi.h>



xy::xy(int Nx, int Ny, int Nz, xoshiro256 &random_generator,
       std::ostream &log)
	: g(Nx, Ny, Nz), rng(random_generator), T(1.0), beta(1.0),
	  xlo(1), xhi(Nx-2), ylo(1), yhi(Ny-2), zlo(1), zhi(Nz-2),
	  boundary_cycle(0), my_rank(-1), my_neighbors{0,0,0,0,0,0}, log(log)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if (comm_size == 1) {
		std::cerr << "Serial run!\n";
	} else {
		if (my_rank == 0) {
			std::cerr << "Parallel run over "
			          << comm_size << " processes!\n";
		}
	}	
}



void xy::rand_init_grid()
{
	long N = g.Nx*g.Ny*g.Nz;
	
	auto u = std::uniform_real_distribution<double>(0, 2*pi);
	for (long i = 0; i < N; ++i) {
		g.s[i] = u(rng);
	}
}

// "fake" the model to be simple ising.
void xy::rand_init_grid_ising()
{
	long N = g.Nx*g.Ny*g.Nz;
	auto u = std::uniform_real_distribution<double>(0, 1);
	for (long i = 0; i < N; ++i) {
		if (u(rng) > 0.5) {
		        g.s[i] = 1;
		} else {
			g.s[i] = -1;
		}
	}
}


double xy::average_spin() const
{
	long N = g.Nx*g.Ny*g.Nz;
	double c = 1.0 / static_cast<double>(N);
	double ts = 0.0;
	for (long i = 0; i < N; ++i) {
		ts += g.s[i]*c;
	}

	return ts;
}

double xy::average_up_spin() const
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


void xy::print_xy_grid_ising(std::ostream &out)
{
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


void xy::exchange_boundaries()
{
	log << "Starting boundary exchange cycle.\n";
	// The size depends on the current cycle.
	// in 2D, boundary cycle goes top (0), right(2), bottom(1), left(3)
	// so the size cycle goes Ny, Nx, Ny, Nx.

        	
	int boundary_size = -1;
	int direction_send = boundary_cycle;
	int direction_recv = -1;
	switch(direction_send) {
	case xy_model_grid::TOP:
		direction_recv = xy_model_grid::BOTTOM;
		break;
	case xy_model_grid::BOTTOM:
		direction_recv = xy_model_grid::TOP;
		break;
	case xy_model_grid::RIGHT:
		direction_recv = xy_model_grid::LEFT;
		break;
	case xy_model_grid::LEFT:
		direction_recv = xy_model_grid::RIGHT;
		break;
	default:
		if (my_rank == 0) {
			std::cerr << "Something tragic just happened.\n";
		}
		
		if (USE_MPI) {
			MPI_Abort(MPI_COMM_WORLD,1);
		} else {
			assert(false);
		}
	}

	log << "Receiving from " << boundary2name[direction_recv]
	    << ", sending to " << boundary2name[direction_send] << "\n";
	
	int send_to   = my_neighbors[direction_send];
	int recv_from = my_neighbors[direction_recv];

	log << "Exchanging boundary " << boundary2name[direction_send]
	    << "; I send to " << send_to << " and receive from "
	    << recv_from << "\n";
	if (direction_send == xy_model_grid::TOP ||
	    direction_send == xy_model_grid::BOTTOM) {
		boundary_size = g.Nx-2;
	} else {
		boundary_size = g.Ny-2;
	}
	
	std::vector<double> recv_boundary(boundary_size);
	std::vector<double> send_boundary(boundary_size);
	
	int top_y_offset = 1, bottom_y_offset = g.Ny-2;
	int right_x_offset = g.Nx-2, left_x_offset = 1;
	long max_idx = g.Nx*g.Ny;
	
	// Fill boundary. The boundary you send is the second-last
	// row or column, and you store it in the last row or column
	// on the other side.
	switch(direction_send) {
	case xy_model_grid::TOP: {
		int y_offset = top_y_offset;
		log << "  Sending boundary at (x,y) = (x," << y_offset << ")\n";
		for (int ix = 0; ix < boundary_size; ++ix) {
			int idx = y_offset*g.Nx + (ix+1);
			send_boundary[ix] = g.s[idx];
		}
		break;
	}
	case xy_model_grid::BOTTOM: {
		int y_offset = bottom_y_offset;
		log << "  Sending boundary at (x,y) = (x," << y_offset << ")\n";
		for (int ix = 0; ix < boundary_size; ++ix) {
			int idx = y_offset*g.Nx + (ix+1);
			send_boundary[ix] = g.s[idx];
		}
		break;
	}
	case xy_model_grid::RIGHT: {
		int x_offset = right_x_offset;
		log << "  Sending boundary at (x,y) = (" << x_offset << ",y)\n";
		for (int iy = 0; iy < boundary_size; ++iy) {
			int idx = x_offset + (iy+1)*g.Nx;
			send_boundary[iy] = g.s[idx];
		}
		break;
	}
	case xy_model_grid::LEFT: {
		int x_offset = left_x_offset;
		log << "  Sending boundary at (x,y) = (" << x_offset << ",y)\n";
		for (int iy = 0; iy < boundary_size; ++iy) {
			int idx = x_offset + (iy+1)*g.Nx;
			send_boundary[iy] = g.s[idx];
		}
		break;
	}
		
	default:
		break;
	};

	
	log << "Filled send_boundary.\n";

	
	int send_tag = 1000 + my_rank;
	int recv_tag = 1000 + recv_from;
	MPI_Status mpi_status;
	int mpi_send_status = MPI_SUCCESS;

	log << "Calling Sendrecv; receiving from " << recv_from
	    << ", sending to " << send_to << "...\n";

	log.flush();

	if (USE_MPI) {
		mpi_send_status = MPI_Sendrecv(send_boundary.data(),
		                               boundary_size, MPI_DOUBLE,
		                               send_to, send_tag,
		                               recv_boundary.data(),
		                               boundary_size, MPI_DOUBLE,
		                               recv_from, recv_tag,
		                               MPI_COMM_WORLD,
		                               &mpi_status);
	} else {
		mpi_status.MPI_ERROR = 0;
		mpi_send_status = 0;
	}
	
	log << "send/receive status: " << mpi_send_status << "/"
	    << mpi_status.MPI_ERROR << ".\n";
	
	log << "Sent following boundary: (";
	for (int i = 0; i < boundary_size; ++i) {
		if (send_boundary[i] > 0) log << "+";
		else log << " ";
	}
	log << ")\n";
	
	log << "Received following boundary: (";
	for (int i = 0; i < boundary_size; ++i) {
		if (recv_boundary[i] > 0) log << "+";
		else log << " ";
	}
	log << ")\n";

	
	// Now set your own boundary accordingly:
	switch(direction_send) {
	case xy_model_grid::TOP: {
		// If you sent to TOP, you received from BOTTOM:
		int y_offset = bottom_y_offset+1;
		log << "  Assigning boundary to y = " << y_offset << "\n";
		for (int ix = 0; ix < boundary_size; ++ix) {
			int idx = y_offset*g.Nx + (ix+1);
			if (idx < 0 || idx >= max_idx) {
				log << "Index " << idx << " either negative or "
				    << ">= " << max_idx << "!\n";
				assert(false && "Index out of bounds!");
			}
			log << "    idx = " << idx << "; (ix,iy) = (" << ix+1
			    << ", " << y_offset << ").\n";
			log << "    s[" << idx << "] was " << g.s[idx]
			    << ", is now ";
			g.s[idx] = recv_boundary[ix];
			log << g.s[idx] << "\n";
		}
		break;
	}
	case xy_model_grid::BOTTOM: {
		int y_offset = top_y_offset-1;
		for (int ix = 0; ix < boundary_size; ++ix) {
			int idx = y_offset*g.Nx + (ix+1);
			if (idx < 0 || idx >= max_idx) {
				log << "Index " << idx << " either negative or "
				    << ">= " << max_idx << "!\n";
				assert(false && "Index out of bounds!");
			}
			g.s[idx] = recv_boundary[ix];
		}
		break;
	}
	case xy_model_grid::RIGHT: {
		int x_offset = left_x_offset-1;
		for (int iy = 0; iy < boundary_size; ++iy) {
			int idx = x_offset + (iy+1)*g.Nx;
			if (idx < 0 || idx >= max_idx) {
				log << "Index " << idx << " either negative or "
				    << ">= " << max_idx << "!\n";
				assert(false && "Index out of bounds!");
			}
			g.s[idx] = recv_boundary[iy];
		}
		break;
	}
	case xy_model_grid::LEFT: {
		int x_offset = right_x_offset+1;
		for (int iy = 0; iy < boundary_size; ++iy) {
			int idx = x_offset + (iy+1)*g.Nx;
			if (idx < 0 || idx >= max_idx) {
				log << "Index " << idx << " either negative or "
				    << ">= " << max_idx << "!\n";
				assert(false && "Index out of bounds!");
			}
			g.s[idx] = recv_boundary[iy];
		}
		break;
	}
		
	default:
		break;
	};

	//                           0  1  2  3
	int boundary_rotate_2D[4] = {3, 2, 0, 1};
	// Update to the next boundary cycle.
	boundary_cycle = boundary_rotate_2D[boundary_cycle];
	log << "Set up next cycle to " << boundary2name[boundary_cycle]
	    << ".\n\n";
	
	// TODO: Change xlo, ylo, xhi and yhi to match the new exchange dir
	switch(boundary_cycle) {
	case xy_model_grid::TOP:
		// You send to top, receive from bottom. In this case, you
		// cannot touch your own bottom row or sides.
		xlo = 2;
		xhi = g.Nx-3;
		ylo = 1;
		yhi = g.Ny-3;
		break;
	case xy_model_grid::BOTTOM:
		xlo = 2;
		xhi = g.Nx-3;
		ylo = 2;
		yhi = g.Ny-2;
		break;
	case xy_model_grid::RIGHT:
		xlo = 2;
		xhi = g.Nx-2;
		ylo = 2;
		yhi = g.Ny-3;
		break;
	case xy_model_grid::LEFT:
		xlo = 1;
		xhi = g.Nx-3;
		ylo = 2;
		yhi = g.Ny-3;
		break;
	default:
		assert(false && "Reached unreachable case!");
	}
	
}



double xy::spin_energy(int ix, int iy, int iz) const
{	
	int idx = ix + g.Nx * iy;
	int ixp = (ix + 1) % g.Nx;
	int ixm = ix - 1;
	if (ixm < 0) ixm += g.Nx;
	
	int iyp = (iy + 1) % g.Ny;
	int iym = iy - 1;
	if (iym < 0) iym += g.Ny;

	int i_top   = g.Nx*iyp + ix;
	int i_bot   = g.Nx*iym + ix;
	int i_left  = g.Nx*iy + ixm;
	int i_right = g.Nx*iy + ixp;

	double U[4];
	double s_mine = -g.s[idx];
	U[0] = s_mine*g.s[i_top];
	U[1] = s_mine*g.s[i_bot];
	U[2] = s_mine*g.s[i_left];
	U[3] = s_mine*g.s[i_right];

	double Ut = U[0] + U[1] + U[2] + U[3];
	return Ut;
}




int xy::metropolis_move()
{
	std::uniform_int_distribution<int> ux(xlo, xhi);
	std::uniform_int_distribution<int> uy(ylo, yhi);
	std::uniform_real_distribution<double> u(0.0, 1.0);
	int rand_x = ux(rng);
	int rand_y = uy(rng);

	assert(rand_x > 0 && rand_x < g.Nx-1 && "rand_x on boundary!");
	assert(rand_y > 0 && rand_y < g.Ny-1 && "rand_y on boundary!");

	int idx = rand_x + g.Nx * rand_y;
	assert(idx >= 0 && idx <= g.Nx*g.Ny && "idx out of bounds!");
	
	double old_E = spin_energy(rand_x, rand_y);
	g.s[idx] *= -1.0;
	double new_E = spin_energy(rand_x, rand_y);

	double dU = new_E - old_E;
	bool accept = true;

	// If new_E is more negative than old_E, dU < 0.
	if (dU < 0) {
		// accept
	} else if (std::exp(-beta*dU) > u(rng)) {
		// accept as well.
	} else {
		accept = false;
	}

	if (!accept) {
		g.s[rand_x + g.Nx*rand_y] *= -1.0;
		return 0;
	} else {
		return 1;
	}
}


double xy::total_spin_energy() const
{
	long N = g.Nx*g.Ny;
	std::vector<double> spin_energies(N, 0.0);
	for (long iy = 1; iy < g.Ny - 1; ++iy) {
		for (long ix = 1; ix < g.Nx-1; ++ix) {
			long ii = ix + iy*g.Nx;
			spin_energies[ii] = spin_energy(ix, iy);
		}
	}

	std::vector<double> all_spin_energies(N,0.0);
	MPI_Allreduce(spin_energies.data(), all_spin_energies.data(), N,
	              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	double U_total = 0.0;
	for (long i = 0; i < N; ++i) {
		U_total += all_spin_energies[i];
	}
	return U_total;
}


void xy::set_my_neighbors(int neighs[6])
{
	my_neighbors[0] = neighs[0];
	my_neighbors[1] = neighs[1];
	my_neighbors[2] = neighs[2];
	my_neighbors[3] = neighs[3];
	my_neighbors[4] = neighs[4];
	my_neighbors[5] = neighs[5];
}

void xy::set_temp(double new_T)
{
	T = new_T;
	beta = 1.0 / T;
}


void xy::set_grid(const xy_model_grid &new_grid)
{
	assert(new_grid.Nx == g.Nx &&
	       new_grid.Ny == g.Ny &&
	       new_grid.Nz == g.Nz &&
	       "Grid sizes inconsistent, cannot assign!");

	for (int j = 0; j < g.Ny; ++j) {
		for (int i = 0; i < g.Nx; ++i) {
			g.s[i + j*g.Nx] = new_grid.s[i + j*g.Nx];
		}
	}
	
}
