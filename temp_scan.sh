#!/bin/zsh
#
# Temperature scan, to see if everything sort-of works:
#

NX=128
NY=256

for SEED in 1 2 3 4 5;
do
	for T in 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5
	do
		DIR="T_"$T
		mkdir -p $DIR
		cd $DIR
		mkdir -p grids_out
		OUT="ising_S_"$SEED".dat"
		mpiexec -np 8 --use-hwthread-cpus ../mpi_xy --ntiles-x 4 --ntiles-y 2 --Nx $NX --Ny $NY -T $T -s 1000 > $OUT

		cd ../
	done
done

