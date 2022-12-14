# FYS3150
# Project 5

#### Solving the Schr. Eq.
This repository contains the following files:

##### main.cpp

Main file. 
Contains the relevant methods to compute the A and B matrices as well as initialising the Gaussian wave packet.
Carries out the Crank-Nicolson method in solving the SE.

To compile, run the command
```
	g++ main.cpp -o main -std=c++11 -lblas -llapack
```
To run the file, type
```
	main <potential filename> <timestep> <total time> <x_c> <sigma_x> <p_x> <y_c> <sigma_y> <p_y>
```
The potential file is ''\<potential filename\>.dat'' and is generated by the python script potential_generator.py. The six final arguments describe the parameters of the initial Gaussian wave packet.
Running main with incorrect no. of arguments prints a help message.
The program saves as binary files the probability function p and the real and imaginary components of the wavefunction u for all timesteps.

##### problem7.py
Plots deviation of total probability from 1 for no slit and for double slit as function of time.

##### problem8.py
Plots colourmaps of probability function and real & imaginary parts of wavefunction for t=0, 0.001, 0.002.

##### problem9.py
Plots the conditional probability along y for x=0.8 and t=0.002, for single, double and triple slit setups.

##### potential_generator.py
Generates a .dat file containing the potential v_ij for no, single, double and triple slit setups.
