/************************************************************************** 
 * surface.c                                                                                                                 
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
#include "equil_util.h"
#include "findmodel.h"

#include "surface.h"




/* Interpolates to find the value of a metric function 
   at the star's surface */
double SurfaceMetric(double *s_vec, double **f, double s_surf, int s, int m)
{
  
  double value;
  
  double metric[4];

  metric[0] = f[s][m];
  metric[1] = f[s+1][m];
  metric[2] = f[s+2][m];
  metric[3] = f[s+3][m];

  value = interpolate(s_vec, metric, s_surf);
  
  return value;
}

/****************************************************************************
 * Surface finds the tabulated values of s, r_is (isotropic radius), and r 
 *       (Schwartzchild radius) as a function of mu.
 ***************************************************************************/

int Surface( EOS *eos, NeutronStar *star){
  int i,j,k;  
  double 
    **enthalpy,      /*Enthalpy*/
    gama,          /*Gamma potential*/
    rho,           /*Rho potential*/
    alpha,
    omega,
    gama_s,
    gama_mu,
    rho_s,
    rho_mu,
    omega_s,
    omega_mu,
    *s_gp,           /*The value of s*/
    *mu,             /*The value of mu*/
    *s_surf,         /*Surface, with s*/
    *r_is_surf,      /*Surface, with isotropic r*/
    *r_surf,         /*Surface, with Schwartzchild r*/
    *gravity_surf,   /*Surface gravity (acceleration) */
    enthalpy_min,    /*The minimum enthalpy, used to find surface of the star*/
    r_e,             /*Used to convert s to r_is, r*/
    s_edge[4],
    enthal_edge[4],
    s_str[4];

  double Omega, drds, r, sintheta, vel, vel_s, vel_mu,
    partial_r, partial_th, grav;

  char filename[80]; 
  FILE *fp;
  
  enthalpy = star->enthalpy;

  s_gp = star->metric.s_gp;
  mu = star->metric.mu;
  s_surf = star->s_surf;
  r_is_surf = star->r_is_surf;
  r_surf = star->r_surf;
  gravity_surf = star->gravity_surf;
  r_e = star->r_e;
  enthalpy_min = star->enthalpy_min;

  Omega = star->Omega * sqrt(KAPPA)/C;

  for(i=1; i<= MDIV; i++){
    s_surf[i] = 0.0;
    r_is_surf[i] = 0.0;
    r_surf[i] = 0.0;
  }

  /*Find s(mu) by discovering where enthalpy goes to min, and extrapolation*/
  /*Find r_is(mu) by a simple conversion of s(mu)*/
  /*Find r(mu) by gama, rho from interpolation with s(mu)*/

  for(i=1; i<= MDIV; i++){
    j=1;
    while(enthalpy[j][i] >= enthalpy_min){
      j++;
      if(enthalpy[j][i] < enthalpy[j+1][i]){
	break;
      }
    }
    j--;
   
    /* At this point j is the last grid point before the surface is reached */
    if(enthalpy[j+1][i] > enthalpy[j+2][i]){
      //printf("\t\t\tcase 1\n");
    }
    else{
      //printf("\t\t\tcase 2\n");
      j++;
    }

    s_edge[0] = s_gp[j-3];
    s_edge[1] = s_gp[j-2];
    s_edge[2] = s_gp[j-1];
    s_edge[3] = s_gp[j];
    enthal_edge[0] = enthalpy[j-3][i];
    enthal_edge[1] = enthalpy[j-2][i];
    enthal_edge[2] = enthalpy[j-1][i];
    enthal_edge[3] = enthalpy[j][i];

      s_surf[i] = extrapolate(s_edge, enthal_edge, enthalpy_min );  
      r_is_surf[i] = (r_e * s_surf[i]) / (1 - s_surf[i]);

	s_str[0] = s_gp[j-1];
	s_str[1] = s_gp[j];
	s_str[2] = s_gp[j+1];
	s_str[3] = s_gp[j+2];

	/* Interpolate to find the values of the metric functions */
	
	gama = SurfaceMetric(s_str, star->metric.gama, s_surf[i], j-1, i);
	rho = SurfaceMetric(s_str, star->metric.rho, s_surf[i], j-1, i);
	alpha = SurfaceMetric(s_str, star->metric.alpha, s_surf[i], j-1, i);
	omega = SurfaceMetric(s_str, star->metric.omega, s_surf[i], j-1, i);

	gama_s = SurfaceMetric(s_str,star->metric.gama_s, s_surf[i],j-1,i);
	rho_s = SurfaceMetric(s_str,star->metric.rho_s, s_surf[i],j-1,i);
	omega_s = SurfaceMetric(s_str,star->metric.omega_s, s_surf[i],j-1,i);

	gama_mu = SurfaceMetric(s_str,star->metric.gama_mu, s_surf[i],j-1,i);
	rho_mu = SurfaceMetric(s_str,star->metric.rho_mu, s_surf[i],j-1,i);
	omega_mu = SurfaceMetric(s_str,star->metric.omega_mu, s_surf[i],j-1,i);

	if(i == 1){
	  star->frame = omega;
	}
	
	r_surf[i] = pow(KAPPA,0.5) * r_is_surf[i] * exp((gama-rho)*0.5);

	drds = r_e/pow(1.0-s_surf[i],2);
	sintheta = sqrt( (1+mu[i])*(1-mu[i]) );

	r = r_is_surf[i];
	vel = (Omega - omega) * exp(-rho) * r * sintheta;

	vel_s = exp(-rho)*r*sintheta *
	  ( (Omega-omega)*(-rho_s + drds/r) - omega_s);

	partial_r = (1-vel*vel)*(gama_s+rho_s) - 2.0*vel*vel_s;
	partial_r *= 1.0/drds;

	if ( sintheta !=0 ){
	  vel_mu =  exp(-rho)*r*sintheta * 
	    ( (Omega-omega)*(-rho_mu - mu[i]/pow(sintheta,2)) - omega_mu);

	  partial_th = (1-vel*vel)*(gama_mu+rho_mu) - 2.0*vel*vel_mu;
	  partial_th *= -1.0*sintheta;
	}
	else partial_th = 0.0;
	
	gravity_surf[i] = sqrt( pow(partial_r,2) + pow(partial_th/r,2));

	gravity_surf[i] *= 0.5*exp(-alpha)/(1.0-vel*vel) * G * C*C/G* 1.0/sqrt(KAPPA);

	grav = G*pow(star->R_e,-2)*star->Mass * 1.0/sqrt(1.0-2.0*G*star->Mass/(C*C*star->R_e));

	star->gravity_surf[i] = gravity_surf[i]/grav; 

	/*
	printf("i=%d sintheta=%lf grav=%6.5e = %6.5e\n",
	       i,sintheta, 
	       gravity_surf[i],
	       G*pow(star->R_e,-2)*star->Mass * 1.0/sqrt(1.0-2.0*G*star->Mass/(C*C*star->R_e))
	       );
	*/
  }

  star->s_surf = s_surf;
  star->r_is_surf = r_is_surf;
  star->r_surf = r_surf;

  return 0;
}



/****************************************************************************
 * LegnFit finds the Legndre Coefficients that fit the 
 *      r_is expansion. The A_0, A_2, and A_4 are given. 
 ***************************************************************************/

int LegnFit( EOS *eos, NeutronStar *star){
  double
    *mu,             /*The value of mu*/
    P_zero[MDIV+1], /*The product of r(mu) and P_k, integrated to get coefficients*/
    P_two[MDIV+1],
    P_four[MDIV+1],
    P_six[MDIV+1],
    A_zero,          /*Legendre Coefficients*/
    A_two,
    A_four,
    A_six;

  double g_zero, g_two, g_four, g_six;
  
  int i;
 
  mu = star->metric.mu;

  for(i=1; i<=MDIV; i++){
    P_zero[i] = 1.0;
    P_two[i]  = 0.5 * ((3.0 * pow(mu[i],2)) - 1.0);
    P_four[i] = 0.125 * ((35.0 * pow(mu[i],4)) - (30.0 * pow(mu[i],2)) + 3.0);
    P_six[i]  = 0.0625 * ((231.0*pow(mu[i],6)) - (315.0*pow(mu[i],4)) 
			  + (105.0*pow(mu[i],2)) - 5.0);
  }

  A_zero = BodeInt(star->r_surf,P_zero);
  A_two  = BodeInt(star->r_surf,P_two);
  A_four = BodeInt(star->r_surf,P_four);
  A_six  = BodeInt(star->r_surf,P_six);

 
  /* Make the curve fit dimensionless by dividing by the Schwarschild Equatorial Radius */
  A_zero = A_zero / star->R_e - 1.0;
  A_two = A_two / star->R_e;
  A_four = A_four / star->R_e;
  A_six = A_six / star->R_e;

  /*
  printf("a0=%lf a2=%lf a4=%lf a6=%lf,   0=%6.5e \n",
	 A_zero, A_two, A_four, A_six,
	 A_zero -0.5*A_two + 3.0/8.0*A_four - 5.0/16.0*A_six);
  */

  star->R_zero = A_zero;
  star->R_two = A_two;
  star->R_four = A_four;
  star->R_six = A_six;


  g_zero = BodeInt(star->gravity_surf,P_zero);
  g_two  = BodeInt(star->gravity_surf,P_two);
  g_four  = BodeInt(star->gravity_surf,P_four);
  g_six  = BodeInt(star->gravity_surf,P_six);

  /*
  printf("g0=%lf g2=%lf g4=%lf g6=%lf,   0=%6.5e \n",
	 g_zero, g_two, g_four, g_six,
	 g_zero -0.5*g_two + 3.0/8.0*g_four - 5.0/16.0*g_six - star->gravity_surf[1]);
  */

  star->g_zero = g_zero;
  star->g_two = g_two;
  star->g_four = g_four;
  star->g_six = g_six;

  return 0;
}


/****************************************************************************
 * BodeInt convolves two functions and integrates using Bode's Law.
 ***************************************************************************/

double  BodeInt(double *f1, double *f2){
 
  int i;
  double integrand[MDIV+1];
  double answer=0.0;

  for (i=1; i<=MDIV; i++){
    integrand[i] = f1[i] * f2[i];
  }

  
  /*Bodes's Rule*/
  if( ( (MDIV-1)%4 ) == 0 ){
    answer = 14.0 * (integrand[1] + integrand[MDIV]);
    for(i=2; i<MDIV; i++){
      if(i%2 == 0){
	answer += 64.0 * integrand[i];
      }
      if((i-1)%4 == 0){
	answer += 28.0 * integrand[i];
      }
      if((i-3)%4 == 0){
	answer += 24.0 * integrand[i];
      }
    }
  } else {  /*Includes Simpson's Rule for the last three points*/
    answer = 14.0 * (integrand[1] + integrand[MDIV-2]) + 
      (45.0/3.0)*(integrand[MDIV-2] + 4.0*integrand[MDIV-1] + integrand[MDIV]);

    for(i=2; i<(MDIV-2); i++){
      if(i%2 == 0){
      	answer += 64.0 * integrand[i];
      }
      if((i-1)%4 == 0){
    	answer += 28.0 * integrand[i];
      }
      if((i-3)%4 == 0){
	answer += 24.0 * integrand[i];
      }
    }
  }
 
  answer *= (1.0/45.0)  * (1.0 / (MDIV - 1.0));

  return answer;
}

