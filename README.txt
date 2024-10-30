RNS

This code computes the structure of rapidly rotating neutron star. This version computes only one stellar model, but it easily be altered to compute a series of models. Based on code written by Nik Stergioulas.

1. INPUT: Provide the code with an equation of state file, a central energy densities and either a spin frequency or a value of r_ratio.

2. OUTPUT: Mass, Radius

3. COMPILATION: At the command line type "make rns" to create the executable rns
> make rns
The warnings produced should not be problematic. It's possible that you may have to change the name of your compiler if you are not using the gnu compiler. 

4. RUN: To run, you will need to choose an equation of state file

EOS FILES: I keep my equation of state files in a directory called "eos" that has the relative path given below. The string "../eos/eosAPR" could be replace with just "eosAPR" if your eos files are in the same directory as the executable. 

EXECUTION: To create a models using eos APR spinning at 300 Hz type at the command line 

> rns -f "../eos/eosAPR" -e 1e15 -s 300
This command should take about 15 seconds to run.

Sample Output:
spin=300
spin=300
../eos/eosAPR,  MDIVxSDIV=151x301
e_c 	 Mass 	 Radius	 Spin  
e15 	 Msun 	 km    	 Hz    
1 	 1.453 	 11.48 	 300.0 




