/************************************************************************** 
*                         rns.c                                           *
*                                                                         *
* Computes the structure of a rapidly rotating neutron star               *
* Based on code written by Nik Stergioulas                                *
*                                                                         *
* Routines to compute the Kepler limit and other useful things can be 
* found in the file findmodel.c
*                                                                         *
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
      
  int i, ierr;
  double
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
    r_ratio=1.0,
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
   spin_freq=0.0,                  /* Spin Frequency */
    Gamma_P=1.0;                      /* Gamma for polytropic EOS */  
                

  
  int nummodels=1;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */



  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly" or "quark"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;  

      case 'b':
	sscanf(argv[i+1],"%lf",&B);
	B *= 1.602e33*KSCALE;
	break;       

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'd':
	/* CHOOSE THE NAME OF THE OUTPUT DIRECTORY */
	sscanf(argv[i+1],"%s",data_dir);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	break;

      case 'l':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_max);
	if(strcmp(eos_type,"poly")!=0)
	  e_max *= C*C*KSCALE;
	break;

      case 'n':
	sscanf(argv[i+1],"%d",&nummodels);
	break;

     case 's':
	/* CHOOSE THE SPIN FREQUENCY (HZ) */
	sscanf(argv[i+1],"%lf",&spin_freq);
	printf("spin=%g\n",spin_freq);
	break;

      case 'r':
	/* CHOOSE r_ratio (instead of spin) */
	sscanf(argv[i+1],"%lf",&r_ratio);
	//printf("spin=%g\n",r_ratio);
	break;

	

      case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        //fprintf(stderr,"     quark : simple quark model \n"); 
	//fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -d directory output goes to \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	//fprintf(stderr,"  -l largest central energy density \n");
	fprintf(stderr,"  Specify either s or r: \n");
	fprintf(stderr,"     -s spin frequency in Hz (0)\n");
	fprintf(stderr,"     -r r_ratio (1.0)\n");
	fprintf(stderr,"  -h this menu\n\n");
	exit(1);
	break;  
      }
    }


  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */

  ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,
		    &eos, &star);

  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;

  printf("e_c \t Mass \t Radius\t Spin  \n");
  printf("e15 \t Msun \t km    \t Hz    \n");

  // Compute a spherical star with the correct central energy density
    ierr = MakeSphere(&eos, &star, e_center);

    if (r_ratio != 1.0 && spin_freq == 0.0) // r_ratio was set (spin was not)
        ierr = rns( r_ratio,  e_center, &eos, &star);
    
    if (spin_freq != 0.0) // Spin frequency was set
      ierr = SetSpin(&eos, &star, e_center, spin_freq);

    //ierr = Surface(&eos,&star);

    //orbit(&eos, &star);
  
  
    printf("%g \t %4.3f \t %4.2f \t %4.1f \n",
	      star.e_center, star.Mass/MSUN, star.R_e*1e-5, star.Omega/(2.0*PI));
    
 

  return 0;
}









