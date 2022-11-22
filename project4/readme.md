# FYS4150
# Project 4
#### Ising model, Markov Chain Monte Carlo sampling and parallelisation

This repository contains the following files:

###### src/Lattice.cpp
Defines a Lattice class containing spins. Has methods to flip spin, calculate E and M.

###### main.cpp
Main file. Contains the metropolis algorithm. To compile, run the command
	g++ main.cpp src/Lattice.cpp -fopenmp -o main -std=c++11 (optimisation -O1 or -O2 or -O3) (link to armadillo)
To run the file, type
	main L n_cycles ordered_lattice save_results temperature
OR
	main L n_cycles ordered_lattice save_results minimum_temperature maximum_temperature n_temperatures use_parallelisation
and hit enter. Here,
* ordered_lattice is a bool = 0, 1. If true, initial configuration is all parallel spins. If false, randomly initialised spin lattice.
* save_results is a bool = 0, 1. If true, saves cycle number, energy and magnetisation to a .csv file
* n_temperatures is no. of temperatures to scan, including minimum and maximum.
* use_parallelisation is a bool = 0, 1. If true, uses OpenMP to run parallel threads in loop over temperature values.

(Running main without arguments prints a message about required arguments.)

###### analytical.py
Calculates the analytical values for L=2.

###### problem4.py
Plots heat capacity and magnetic susceptibility as function of no. of MC cycles, with comparison to analytical results.
The .csv data files were created using
	main 2 3000000 0 1 1.0

###### problem5.py
Plots mean energy per spin and mean magnetisation per spin as function of no. of MC cycles, by reading appropriate .csv datafiles.
The .csv files were created using
	main 20 3000000 ordered_lattice 1 T
for ordered_lattice = 0, 1 and T = 1.0, 2.4.

###### problem6.py
Plots the distribution of energy per spipn for T=1.0 and T=2.4.
Reads .csv files gotten by running
	main 20 3000000 1 1 T
for T = 1.0, 2.4.

###### problem7.py
Plots results from a timing test on parallelised code.

###### problem8.py
Plots mean energy per spin, mean magnetisation per spin, heat capactiy and magnetic susceptibility as function of temperature, for lattice sizes L = 40,60,80,100. Temperature ranges from 2.1 to 2.4. 
The .csv files it reads were gotten by running
	main L 1000000 0 1 2.1 2.4 9 1
and
	main L 1000000 0 1 2.27 2.35 9 1
for L = 40,60,80,100.
####### warning: this creates a lot of .csv files: 4*20 = 80 of them...
