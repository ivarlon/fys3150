
________________.___. _________________  ____ ._______________   
\_   _____/\__  |   |/   _____/\_____  \/_   ||   ____/\   _  \  
 |    __)   /   |   |\_____  \   _(__  < |   ||____  \ /  /_\  \ 
 |     \    \____   |/        \ /       \|   |/       \\  \_/   \
 \___  /    / ______/_______  //______  /|___/______  / \_____  /
     \/     \/              \/        \/            \/        \/ 
__________                   __               __    ________     
\______   \_______  ____    |__| ____   _____/  |_  \_____  \    
 |     ___/\_  __ \/  _ \   |  |/ __ \_/ ___\   __\   _(__  <    
 |    |     |  | \(  <_> )  |  \  ___/\  \___|  |    /       \   
 |____|     |__|   \____/\__|  |\___  >\___  >__|   /______  /   
                        \______|    \/     \/              \/    
                                                                 
                                                                 

  ___               _             _                 
 | _ \___ _ _  _ _ (_)_ _  __ _  | |_ _ _ __ _ _ __ 
 |  _/ -_) ' \| ' \| | ' \/ _` | |  _| '_/ _` | '_ \
 |_| \___|_||_|_||_|_|_||_\__, |  \__|_| \__,_| .__/
                          |___/               |_|   


This repo contains basically the following files:
* main.cpp
Contains functions to run the various parts of the project.
Uncomment a line in the main function to run the desired function.
How to compile:
run the following command in the terminal:
g++ main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/ -o main.exe -O2 -std=c++11 -larmadillo
then run
./main.exe

* src/Particle.cpp
defines the Particle class, which holds particle parameters such as mass, charge, position and velocity.

* src/PenningTrap.cpp
defines the PenningTrap class to contain Particle instances and evolve system in time using RK4 or FE.

* test_run_one_particle.py
Plots results for one particle stored in appropriate .csv files, gotten by running the method test_run_one_particle in main.cpp.

* simulate_two_particles.py
Plots results for two particles stored in appropriate .csv files, gotten by running the method simulate_two_particles in main.cpp

* performance_tests.py
Plots numerical errors for RK4 and FE stored in appropriate .csv files, gotten by running the method performance_tests in main.cpp

* resonance_exploration.py
Plots resonance spectrum based on data in appropriate .csv files, gotten by running the method resonance_exploration in main.cpp

* resonance_finestructure.py
Plots resonance spectrum based on data in appropriate .csv files, gotten by running the method resonance_finestructure in main.cpp
