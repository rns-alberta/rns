/************************************************************************** 
*                            stableorbit.c                                        
*                                                                         
*       The routines in here are used to find and fit the inner most
*       stable orbit of the star.
*                                                                         
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
#include "findmodel.h"

#include "surface.h"
#include "stableorbit.h"

void orbit(EOS *eos, NeutronStar *star){
  double s_surf_eqa = star->s_surf[1];
  double s_used[4], potP[4], potN[4], velP[4], velN[4];
  double potenP, potenN, pottest;
  double potzeP[SDIV+1], potzeN[SDIV+1];
  double s_stable, zero = RMIN;
  int k, j, p;
  double velcheck;

  //printf("\t\t\t entered orbit, s_surf_eqa = %g\n",s_surf_eqa); 
  /* The positive case */
  /* j is the last grid point before the surface of the star is reached at the equator*/
  j = 0;
  while( star->metric.s_gp[j] <= s_surf_eqa){
    j++;
  }
  j--;
  //printf("surface of star at grid j=%d \n",j);

  
  /* Check if fabs(vel) < 1.0 at the surface of the star */
  p = j+1;
  //p=j;
  velcheck = 1.1;
  while(fabs(velcheck) >= 1){
    pottest = pot_dr2P(eos, star, p, 1, &velcheck);
    //  printf("\t\t\t p = %d \t\tpottest=%6.5e\tvelcheck=%6.5e\n",p, pottest, velcheck);
    p++;
  }
  p--;

  //printf("\t\t\t p=%d \t\tpottest=%6.5e\n",p, pottest);


  /* Separate and deal with Dpot being positive or negative */
  if(p == j + 1){
    if(pottest > 0.0){
      //printf("ISCO is outside the star's surface\n");
      /* Find where Dpot = 0.0 */
      while(pottest > 0.0){
	pottest = pot_dr2P(eos, star, p+1, 1, &velcheck);
	//	printf("\t\t\t p=%d \t\tpottest=%6.5e\n",p+1, pottest);
	p++;
      }
      p--;

      potP[0] = pot_dr2P(eos, star, p, 1, &velcheck);
      potP[1] = pot_dr2P(eos, star, p+1, 1, &velcheck);
      potP[2] = pot_dr2P(eos, star, p+2, 1, &velcheck);
      potP[3] = pot_dr2P(eos, star, p+3, 1, &velcheck);
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[3], s_used[3]);
      s_stable = polation(s_used, potP, zero);
      if(s_stable < s_used[0]){
	s_stable = s_used[0] - potP[0]*((s_used[1] - s_used[0])/(potP[1] - potP[0]));
      }
      if(s_stable > s_used[1]){
	s_stable = s_used[0] - potP[0]*((s_used[1] - s_used[0])/(potP[1] - potP[0]));
      }
      k = p + 1;
      //      printf("\t\t ISCO at s = \t%6.5e\n", s_stable);
      //printf("\t\t\t%6.5e\t%6.5e\t%6.5e\n", s_used[0], s_stable, s_used[1]);
      //printf("\t\t\tcase 1\n");
    }
    else{ // pottest < 0
      /* Simplest case, the surface is s_stable */
      s_used[0] = star->metric.s_gp[j+1];
      s_used[1] = star->metric.s_gp[j+2];
      s_used[2] = star->metric.s_gp[j+3];
      s_used[3] = star->metric.s_gp[j+4];
      potP[0] = pot_dr2P(eos, star, j+1, 1, &velcheck);
      potP[1] = pot_dr2P(eos, star, j+2, 1, &velcheck);
      potP[2] = pot_dr2P(eos, star, j+3, 1, &velcheck);
      potP[3] = pot_dr2P(eos, star, j+4, 1, &velcheck);
      
      /* Set it up as the surface, since all orbits are stable */
      s_stable = s_surf_eqa;
      k = j+1;
      //printf("\t\t\t%6.5e\n", s_stable);
      //printf("\t\t\tcase 2\n");
    }
  }
  else{
    if(pottest > 0.0){
      /* Find where Dpot = 0.0 */
      while(pottest > 0.0){
	pottest = pot_dr2P(eos, star, p+1, 1, &velcheck);
	//printf("\t\t\t\t\tpottest=%6.5e\n", pottest);
	p++;
      }
      p--;

      potP[0] = pot_dr2P(eos, star, p, 1, &velcheck);
      potP[1] = pot_dr2P(eos, star, p+1, 1, &velcheck);
      potP[2] = pot_dr2P(eos, star, p+2, 1, &velcheck);
      potP[3] = pot_dr2P(eos, star, p+3, 1, &velcheck);
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potP[3], s_used[3]);
      s_stable = polation(s_used, potP, zero);
      if(s_stable < s_used[0]){
	s_stable = s_used[0] - potP[0]*((s_used[1] - s_used[0])/(potP[1] - potP[0]));
      }
      if(s_stable > s_used[1]){
	s_stable = s_used[0] - potP[0]*((s_used[1] - s_used[0])/(potP[1] - potP[0]));
      }
      k = p + 1;
      //printf("\t\t\t%6.5e\n", s_stable);
      //printf("\t\t\t%6.5e\t%6.5e\t%6.5e\n", s_used[0], s_stable, s_used[1]);
      //printf("\t\t\tcase 3\n");
    }
    else{
      /* Find where vel = 1.0 */      
      potP[0] = pot_dr2P(eos, star, p, 1, &velcheck);
      velP[0] = velcheck;
      potP[1] = pot_dr2P(eos, star, p+1, 1, &velcheck);
      velP[1] = velcheck;
      potP[2] = pot_dr2P(eos, star, p+2, 1, &velcheck);
      velP[2] = velcheck;
      potP[3] = pot_dr2P(eos, star, p+3, 1, &velcheck);
      velP[3] = velcheck;
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", velP[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velP[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velP[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velP[3], s_used[3]);
      s_stable = polation(s_used, velP, 1.0 - RMIN);
      k= p + 1;
      //printf("\t\t\t%6.5e\n", s_stable);
      //printf("\t\t\tcase 4\n");
    }
  }

  //printf("s_stable = %g k=%d \n", s_stable,k);

  /* Determine OMEGA */
  star->orbitP = OMEGA(eos, star, k-1, s_stable, 1)*C/(sqrt(KAPPA))/(2.0*PI);
  //printf("Orbital Freq at ISCO = \torbitP=%6.5e Hz (co-rotating)\n", star->orbitP);

 

  /* Counter-rotating Orbits */
  /* j is the last grid point before the surface of the star is reached at the equator*/
  j = 0;
  while( star->metric.s_gp[j] <= s_surf_eqa){
    j++;
  }
  j--;

  /* Check that fabs(vel) < 1.0 at the surface of the star */
  p = j+1;
  velcheck = 1.1;
  while(fabs(velcheck) >= 1){
    pottest = pot_dr2P(eos, star, p, 0, &velcheck);
    p++;
  }
  p--;

  /* Separate and deal with Dpot being positive or negative */
  if(p == j + 1){
    if(pottest > 0.0){
      /* Find where Dpot = 0.0 */
      while(pottest > 0.0){
	pottest = pot_dr2P(eos, star, p+1, 0, &velcheck);
	p++;
      }
      p--;

      potN[0] = pot_dr2P(eos, star, p, 0, &velcheck);
      potN[1] = pot_dr2P(eos, star, p+1, 0, &velcheck);
      potN[2] = pot_dr2P(eos, star, p+2, 0, &velcheck);
      potN[3] = pot_dr2P(eos, star, p+3, 0, &velcheck);
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[3], s_used[3]);
      s_stable = polation(s_used, potN, zero);
      if(s_stable < s_used[0]){
	s_stable = s_used[0] - potN[0]*((s_used[1] - s_used[0])/(potN[1] - potN[0]));
      }
      if(s_stable > s_used[1]){
	s_stable = s_used[0] - potN[0]*((s_used[1] - s_used[0])/(potN[1] - potN[0]));
      }
      k = p + 1;
      //printf("\t\t\t%6.5e\t%6.5e\t%6.5e\n", s_used[0], s_stable, s_used[1]);
      //printf("\t\t\tcase 1\n");
    }
    else{
      /* Simplest case, the surface is s_stable */
      s_used[0] = star->metric.s_gp[j+1];
      s_used[1] = star->metric.s_gp[j+2];
      s_used[2] = star->metric.s_gp[j+3];
      s_used[3] = star->metric.s_gp[j+4];
      potN[0] = pot_dr2P(eos, star, j+1, 0, &velcheck);
      potN[1] = pot_dr2P(eos, star, j+2, 0, &velcheck);
      potN[2] = pot_dr2P(eos, star, j+3, 0, &velcheck);
      potN[3] = pot_dr2P(eos, star, j+4, 0, &velcheck);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[3], s_used[3]);
      potenN = polation(potN, s_used, s_surf_eqa);  
      /* Set it up as the surface, since all orbits are stable */
      s_stable = s_surf_eqa;
      k = j+1;
      //printf("\t\t\tcase 2\n");
    }
  }
  else{
    if(pottest > 0.0){
      /* Find where Dpot = 0.0 */
      while(pottest > 0.0){
	pottest = pot_dr2P(eos, star, p+1, 0, &velcheck);
	p++;
      }
      p--;

      potN[0] = pot_dr2P(eos, star, p, 0, &velcheck);
      potN[1] = pot_dr2P(eos, star, p+1, 0, &velcheck);
      potN[2] = pot_dr2P(eos, star, p+2, 0, &velcheck);
      potN[3] = pot_dr2P(eos, star, p+3, 0, &velcheck);
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", potN[3], s_used[3]);
      s_stable = polation(s_used, potN, zero);
      if(s_stable < s_used[0]){
	s_stable = s_used[0] - potN[0]*((s_used[1] - s_used[0])/(potN[1] - potN[0]));
      }
      if(s_stable > s_used[1]){
	s_stable = s_used[0] - potN[0]*((s_used[1] - s_used[0])/(potN[1] - potN[0]));
      }
      k = p + 1;
      //printf("\t\t\tcase 3\n");
    }
    else{
      /* Find where vel = 1.0 */      
      potN[0] = pot_dr2P(eos, star, p, 0, &velcheck);
      velN[0] = velcheck;
      potN[1] = pot_dr2P(eos, star, p+1, 0, &velcheck);
      velN[1] = velcheck;
      potN[2] = pot_dr2P(eos, star, p+2, 0, &velcheck);
      velN[2] = velcheck;
      potN[3] = pot_dr2P(eos, star, p+3, 0, &velcheck);
      velN[3] = velcheck;
      s_used[0] = star->metric.s_gp[p];
      s_used[1] = star->metric.s_gp[p+1];
      s_used[2] = star->metric.s_gp[p+2];
      s_used[3] = star->metric.s_gp[p+3];
      //printf("\t\t\t%6.5e\t%6.5e\n", velN[0], s_used[0]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velN[1], s_used[1]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velN[2], s_used[2]);
      //printf("\t\t\t%6.5e\t%6.5e\n", velN[3], s_used[3]);
      s_stable = polation(s_used, velN, -1.0 + RMIN);
      k= p + 1;
      //printf("\t\t\tcase 4\n");
    }
  }

  /* Determine OMEGA */
  star->orbitN = OMEGA(eos, star, k-1, s_stable, 0)*C/(sqrt(KAPPA))/(2.0*PI);
  //printf("Orbital Freq at ISCO = \torbitN=%6.5e Hz (counter-rot)\n", star->orbitN);
}

double OMEGA(EOS *eos, NeutronStar *star, int k, double s_stable, int sign){
  double r0s, r1s, rho0r, rho1r, ome0r, ome1r, gam0r, gam1r;
  double OMEGA_P, omega_P, vel_P, r_is_P, rho_P;
  double OMEGAP[4], omega[4], velP[4], r_is[4], rho[4];
  double one[4], two[4], three[4];
  double first[4], second[4], third[4];
  double s_used[4];
  int i, q;
  char filename[80]; 
  FILE *fp;
  
  //printf("\t\t\t\t OMEGA \t\ts_stable=%6.5e\n", s_stable);


  for(i=0; i<4; i++){
    q = k+i;

    /* Convert the derivatives from ds to dr values */
  r0s    = star->r_e * star->metric.s_gp[q] / (1.0 - star->metric.s_gp[q]);
  r1s    = star->r_e * pow(1.0 - star->metric.s_gp[q], -2.0);
  rho0r  = star->metric.rho[q][1];
  rho1r  = star->metric.rho_s[q][1] / r1s;
  ome0r  = star->metric.omega[q][1];
  ome1r  = star->metric.omega_s[q][1] / r1s;
  gam0r  = star->metric.gama[q][1];
  gam1r  = star->metric.gama_s[q][1] / r1s;

    r_is[i]   = r0s;
    rho[i]    = rho0r;
    omega[i]  = ome0r;
    s_used[i] = star->metric.s_gp[q];

    one[i]    = SQ(ome1r*r0s*exp(-rho0r));
    two[i]    = (2.0/r0s)*(rho1r + gam1r);
    three[i]  = -(SQ(rho1r) - SQ(gam1r)); 

    first[i]  = r0s / (2.0 - r0s*(rho1r - gam1r));
    second[i] = ome1r * r0s * exp(-rho0r);
    third[i]  = one[i] + two[i] + three[i];

    if(sign == 1){
      /* The positive case */
      velP[i]   = first[i]*(second[i] + pow(third[i], 0.5));
    }
    if(sign == 0){
      /* The negative case */
      velP[i]   = first[i]*(second[i] - pow(third[i], 0.5));
    }

    OMEGAP[i] = omega[i] + (velP[i] / r_is[i])*exp(rho[i]);
    //printf("\t\t\t\t s=%6.5e \t\tOmega=%6.5e\n", s_used[i],OMEGAP[i]);
    //printf("\t\t\t\t\t\ts_used=%6.5e\n", s_used[i]);


  }
  
  /*  vel_P = polation(velP, s_used, s_stable);
  r_is_P = polation(r_is, s_used, s_stable);
  rho_P = polation(rho, s_used, s_stable);
  omega_P = polation(omega, s_used, s_stable);*/

  //sprintf(filename,"%s/surface/archive/IMSCO.out", eos->data_dir);
  //fp = fopen(filename, "a");
  /*
  OMEGA_P = omega_P + (vel_P / r_is_P)*exp(rho_P);
  printf("\t\t\t\tinterpolate then combine: OMEGA=%6.5e\n", OMEGA_P);
  */
  //fprintf(fp, "%6.5e\t", OMEGA_P);

  OMEGA_P = polation(OMEGAP, s_used, s_stable);
  //printf("\t\t\t\tcombine then interpolate: OMEGA=%6.5e\n", OMEGA_P);  

  //fprintf(fp, "%6.5e\t%d\n", OMEGA_P, sign);
  //fclose(fp);

  return OMEGA_P;
}

double polation(double *yg, double *xg, double xs){
  double one, two, three, four;
  double ys;
  
  one   = (xs-xg[1])*(xs-xg[2])*(xs-xg[3])*yg[0]/
    ((xg[0]-xg[1])*(xg[0]-xg[2])*(xg[0]-xg[3]));
  two   = (xs-xg[0])*(xs-xg[2])*(xs-xg[3])*yg[1]/
    ((xg[1]-xg[0])*(xg[1]-xg[2])*(xg[1]-xg[3]));
  three = (xs-xg[0])*(xs-xg[1])*(xs-xg[3])*yg[2]/
    ((xg[2]-xg[0])*(xg[2]-xg[1])*(xg[2]-xg[3]));
  four  = (xs-xg[0])*(xs-xg[1])*(xs-xg[2])*yg[3]/
    ((xg[3]-xg[0])*(xg[3]-xg[1])*(xg[3]-xg[2]));
  ys = one + two + three + four;
  //printf("\t\t\t\tpolation %6.5e\t%6.5e\t%6.5e\t%6.5e\n", one, two, three, four);
  return ys;
}

double pot_dr2P(EOS *eos, NeutronStar *star, int s, int sign, double *velcheck){
  double dpot_dr2; 
  double one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen;
  double ell, emol, velsq, vel;
  double aay, bee, delsq;
  double first, second, third;
  double rho2r, rho1r, rho0r, ome2r, ome1r, ome0r, gam2r, gam1r, gam0r;
  double rho_ss, ome_ss, gam_ss, r2s, r1s, r0s;
  char filename[80]; 
  FILE *fp;

  double dsdrr;

  //printf("\t\t\t\ts = %d\n", s);

  /* Remember that we are in the plane of the equator, so theta = pi/2 and mu = 0 */

  /* Covert the derivatives from ds to dr values */
  r0s    = star->r_e * star->metric.s_gp[s] / (1.0 - star->metric.s_gp[s]);
  r1s    = star->r_e * pow(1.0 - star->metric.s_gp[s], -2.0);
  r2s    = 2.0 * star->r_e * pow(1.0 - star->metric.s_gp[s], -3.0);
  dsdrr = -2.0*pow(1.0 - star->metric.s_gp[s], 3.0)/(pow(star->r_e,2));

  rho_ss = (star->metric.rho_s[s+1][1] - star->metric.rho_s[s-1][1]) / (star->metric.s_gp[s+1] - star->metric.s_gp[s-1]);
  rho0r  = star->metric.rho[s][1];
  rho1r  = star->metric.rho_s[s][1] / r1s;
  //  rho2r  = star->metric.rho_s[s][1] / r2s + rho_ss / SQ(r1s);
  rho2r  = star->metric.rho_s[s][1] * dsdrr + rho_ss / SQ(r1s);

  ome_ss = (star->metric.omega_s[s+1][1] - star->metric.omega_s[s-1][1]) / (star->metric.s_gp[s+1] - star->metric.s_gp[s-1]);
  ome0r  = star->metric.omega[s][1];
  ome1r  = star->metric.omega_s[s][1] / r1s;
  //ome2r  = star->metric.omega_s[s][1]/ r2s + ome_ss / SQ(r1s);
  ome2r  = star->metric.omega_s[s][1] * dsdrr + ome_ss / SQ(r1s);


  gam_ss = (star->metric.gama_s[s+1][1] - star->metric.gama_s[s-1][1]) / (star->metric.s_gp[s+1] - star->metric.s_gp[s-1]);
  gam0r  = star->metric.gama[s][1];
  gam1r  = star->metric.gama_s[s][1] / r1s;
  //gam2r  = star->metric.gama_s[s][1] / r2s + gam_ss / SQ(r1s);
  gam2r  = star->metric.gama_s[s][1] * dsdrr + gam_ss / SQ(r1s);

  /* Calculate v, L, and E-Om*L*/
  first  = SQ(ome1r * r0s * exp(-rho0r));
  second = 2.0 *( rho1r + gam1r) / r0s;
  third  = -(SQ(rho1r) - SQ(gam1r));
  delsq  = first + second + third;
  aay    = r0s / (2.0 - r0s*(rho1r - gam1r));
  bee    = ome1r * r0s * exp(-rho0r);
  if(sign == 1){
    vel   = aay * (bee + pow(delsq, 0.5)); /* Positive case */
  }
  if(sign == 0){
    vel   = aay * (bee - pow(delsq, 0.5)); /* Negative case */
  }
  //printf("\t\t\t\tvel = %6.5e\n", vel);
  *velcheck = vel;

  velsq = pow(1.0 - SQ(vel), -0.5);
  ell   = vel * velsq * r0s * exp(0.5*(gam0r - rho0r));
  emol  = velsq * exp(0.5*(rho0r + gam0r));

  /* Calculate dpot_dr_dr */
  one      = -rho2r * exp(-rho0r) * SQ(emol);
  two      = SQ(rho1r * emol) * exp(-rho0r);
  three    = 2.0 * rho1r * exp(-rho0r) * emol * ell * ome1r;
  four     = -2.0 * ome2r * ell * emol * exp(-rho0r);
  five     = 2.0 * ome1r * ell * emol * rho1r * exp(-rho0r);
  six      = 2.0 * SQ(ome1r * ell) * exp(-rho0r);
  seven    = -rho2r * exp(rho0r) * SQ(ell / r0s);
  eight    = -SQ(rho1r) * exp(rho0r) * SQ(ell / r0s);
  nine     = 2.0 * rho1r * exp(rho0r) * SQ(ell) * pow(r0s, -3.0);
  ten      = nine;
  eleven   = -6.0 * exp(rho0r) * SQ(ell / SQ(r0s));
  twelve   = -gam2r * exp(gam0r);
  thirteen = -SQ(gam1r) * exp(gam0r);

  dpot_dr2 = one + two + three + four + five + six + seven + eight + nine + ten + eleven + twelve + thirteen;

  //printf("\t\t\t\t*velcheck = %6.5e\n", *velcheck);
  /*if(*velcheck == 0){
     if(sign == 0){
       sprintf(filename,"%s/surface/archive/%1.2lf,%1.2lfminus.out", eos->data_dir, star->e_center, star->r_ratio);
     }
     if(sign == 1){
       sprintf(filename,"%s/surface/archive/%1.2lf,%1.2lfplus.out", eos->data_dir, star->e_center, star->r_ratio);
     }
     fp = fopen(filename, "a");
     fprintf(fp, "%d\t", s);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t", r0s, r1s, r2s);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t%6.5e\t", rho_ss, rho0r, rho1r, rho2r);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t%6.5e\t", ome_ss, ome0r, ome1r, ome2r);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t%6.5e\t", gam_ss, gam0r, gam1r, gam2r);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t", first, second, third);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t", aay, bee, delsq);
     if(sign == 0){
       fprintf(fp, "%6.5e\t%6.5e\t", pow(delsq, 0.5), bee - pow(delsq, 0.5));
     }
     if(sign == 1){
       fprintf(fp, "%6.5e\t%6.5e\t", pow(delsq, 0.5), bee + pow(delsq, 0.5));
     }
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t%6.5e\t", vel, velsq, ell, emol);
     fprintf(fp, "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t", one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen);
     fprintf(fp, "%6.5e\n", dpot_dr2);
     fclose(fp);
  }*/

  // Negative value of dpot_dr2 means stable orbit

  return (dpot_dr2);
}

