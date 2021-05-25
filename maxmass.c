/************************************************************************** 
*                         maxmass.c                                       *
*                                                                         *
* Computes the maximum mass non-rotating neutron star.                    *
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


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
      
  int i, ierr;
  double
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
   r_ratio,                      /* axis ratio */
   spin_freq=0,                  /* Spin Frequency */
   Gamma_P,                      /* Gamma for polytropic EOS */  
   spin_kep,                     /* Kepler Frequency */
   spin_step=100,                /* Step size for frequency changes*/
   r_start=1.0,                  /* Equator-Pole radius ratio*/
   r_kep,                        /* Kepler r_ratio*/
   r_step=0.01,                  /* Step size for r_ratio changes*/  
   Mass,                         /* Gravitational mass */
   Mass_0,                       /* Baryon Mass */
   Omega,			 /* Angular Velocity */
   J,				 /* Angular Momentum, in cgs units */
   R_e,                          /* Circumferential radius at equator */
   Omega_K,                      /* Keplerian velocity of particle orbiting at equator */
   r_e                           /* coord. radius at equator 	*/
   ;
  
  int printflag=1;

  double xx,
    zeta,
    delta,
    Radius,
    dr,
    R_0,                            /* Radius of the zero spin star */
    Ratio,
    kepler_ratio,                   /* r_ratio at the Kepler limit  */
    kepler_freq,                    /* freq at Kepler limit */
    stuff;

  

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  char abbrev[80] = "abbrev";
  char gnu_file[80] = "data-zero/gplot-MRdata.txt";

  FILE *output;
  FILE *gplot;
  

  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'a':
	sscanf(argv[i+1],"%s",abbrev);
	break;

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



      case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        fprintf(stderr,"     quark : simple quark model \n"); 
	fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -d directory output goes to \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	fprintf(stderr,"  -l largest central energy density \n");
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

  //printf("Star set up \n");

  ierr = MaxMass(e_min, e_max, &eos, &star, data_dir);

  if (ierr !=1){

  output = fopen("data-zero/zero.dat","a");
  gplot = fopen(gnu_file,"a");

  xx = G*star.Mass/(star.R_e*C*C);

  e_max = star.e_center;

  fprintf(output,"#EOS e_c M R M/R alpha\n");
  fprintf(output,"%s %g %g %g %g %g \n",
	  eos_file,
	  e_max, star.Mass/MSUN, star.R_e*1e-5, xx,
	  e_max*4.0*PI*pow(star.R_e,3)/(3*star.Mass)*1e15);
   fclose(output);

   fprintf(gplot,"E%s=%g; M%s=%g; R%s=%g; Z%s=%g; A%s=%g\n",
	   abbrev,e_max,
	   abbrev,star.Mass/MSUN,abbrev,star.R_e*1e-5,
	   abbrev,xx,
	   abbrev,e_max*4.0*PI*pow(star.R_e,3)/(3*star.Mass)*1e15);
   close(gplot);

   output= fopen("data-zero/zero.tex","a");

   fprintf(output,"%s & %g & %g & %g & %g  \\\\  \n",
	   eos_file,
	   star.Mass/MSUN,
	   star.R_e*1e-5,
	   xx,
	   sqrt(star.Mass*G/pow(star.R_e,3)));


  }


  /*
  if (!ierr){ 
    output = fopen(data_dir,"w");

    fprintf(output,"Maximum Mass Spherical star for EOS:\n");
    fprintf(output,"%s, MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
    fprintf(output,"e_center = %g e15 g/cm^3\n",star.e_center);
    fprintf(output,"Mass   = %g Msun \n", star.Mass/MSUN);
    fprintf(output,"Mass_0 = %g Msun \n", star.Mass_0/MSUN);
    fprintf(output,"Radius = %g km \n", star.R_e*1e-5);
    fprintf(output,"M/R = %g \n", star.Mass*G/(star.R_e*C*C));
    fprintf(output,"(R^3/GM)^{1/2} = %g \n",
	    sqrt(pow(star.R_e,3)/(G*star.Mass)));

    fclose(output); 
  }
  */

  return 0;
}









