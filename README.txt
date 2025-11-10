RNS

This code computes the structure of rapidly rotating neutron star. This version computes only one stellar model, 
but it is easily altered to calculate a series of models. 

Original code written by Nik Stergioulas. This version has been modified by: Sharon Morsink, Sheldon Campbell, and John Braga.

1. INPUT: Provide the code with an equation of state file, a central energy densities and either a spin frequency or a value of r_ratio.

2. OUTPUT: Mass, Radius

3. COMPILATION: At the command line type "make rns" to create the executable rns
> make rns
The warnings produced should not be problematic. It's possible that you may have to change the name of your compiler if you are not using the gnu compiler. 

4. RUN: To run, you will need to choose an equation of state file

EOS FILES: I keep my equation of state files in a directory called "eos" that has the relative path given below. The string "../eos/eosAPR" could be replaced with just "eosAPR" if your eos files are in the same directory as the executable. 

EXECUTION: To create a model star using eos APR spinning at 300 Hz type at the command line 

> rns -f "../eos/eosAPR" -e 1e15 -s 300
This command should take about 15 seconds to run.

Sample Output:
spin=300
spin=300
../eos/eosAPR,  MDIVxSDIV=151x301
e_c 	 Mass 	 Radius	 Spin  
e15 	 Msun 	 km    	 Hz    
1 	 1.453 	 11.48 	 300.0 

5. Code files:
rns.c = main() Edit the main() procedure to make rns do the calculations you would like it to do!
struct.h = structure definitions. The NeutronStar, Metric, and EOS structures are defined here, so take a look here to see what you're working with!
surface.c = routines that find the star's surface and compute various physical quantities on the surface.

6. Neutron Star Routines:
To compute a star:
(a) Input the desired EOS, central density, and spin frequency
(b) Call SetUpStar() -- This sets up the coordinates, grid, and metric.
(c) Compute a spherical star with the same central density as the desired star: MakeSphere()
(d) Call SetSpin() to get a star with the desired spin frequency
(e) Call Surface() to return a collection of vectors defining physical quantities on the star's surface.



