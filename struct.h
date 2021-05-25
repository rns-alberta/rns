#ifndef STRUCT_H

#define STRUCT_H
/* Development of structures done by Sheldon Campbell */

typedef struct{               
  char eos_file[80];          /*The name of the EOS file as stored on the computer*/
  char eos_type[10];          /*Tabulated, Polytropic, or Quark*/
  char data_dir[80];          /*The location the output legendre files are sent to*/
  double log_e_tab[2001];      /*?*/
  double log_p_tab[2001];      /*?*/
  double log_h_tab[2001];      /*?*/
  double log_n0_tab[2001];     /*?*/
  int n_tab;                  /*?*/
  double Gamma_P;             /*?*/
  double B;                   /*?*/
  double K;                   /*?*/
} EOS;


typedef struct {
  double s_gp[SDIV+1], mu[MDIV+1],          /* The value of s and mu at the grid points */
    **alpha, **gama, **rho, **omega,        /* The four gravitational potentials and derivatives */
    **alpha_s, **gama_s, **rho_s, **omega_s,        
    **alpha_mu, **gama_mu, **rho_mu, **omega_mu,    
    **alpha_ms, **gama_ms, **rho_ms, **omega_ms;    
} Metric;

typedef struct {
  Metric metric;               /*All the geometrical information*/
  double *s_surf;       /*The surface of the star, s(mu)*/
  double *r_is_surf;    /*The surface of the star, r_is(mu), isotropic radius*/
  double *r_surf;       /*The surface of the star, r(mu), Schwarzchild radius*/
  double *surf_der;     /*The surface derivative dr/dtheta */
  double *gravity_surf; /*The acceleration due to gravity at the surface */
  double Mass;                 /*In Solar Mass Units*/
  double Mass_0;               /*In Solar Mass Units*/
  double r_ratio;              /*Ratio of polar to equatoral radius*/
  double e_center;             /*Energy Density at the center of the star*/
  double p_center;             /*Pressure at the center of the star*/
  double h_center;             /*Enthalpy at the center of the star*/
  double e_surface;            /*Energy Density at the surface of the star*/
  double enthalpy_min;         /*The minimum enthalpy, roughly zero, finds surface of the star*/
  double r_e;                  /*Needed to convert from s to isotropic radius (unitless)*/
  double R_e;                  /*Schwarzschild radial coordinate (unitless) */
  double frame;                /*The frame dragging term at the surface of the equator*/
  double Omega;                /*The angular velocity of the star (radians per second)*/
  double Omega_K;              /**/
  double I;                    /*Moment of inertia of the star*/
  double Q;                    /*Quadrupole moment of the star*/
  double q;
  double zeta;                 /*Used in calculating the invariant quadrupole moment*/
  double m_2;                  /*Coordinate invariant quadrupole moment*/
  double ang_mom;              /*Angular momentum of the star, in units of g*cm^2/s*/
  double **energy, **pressure, **enthalpy, **velocity_sq;    /*s by mu grid of Energy Density, Pressure, Enthalpy */
  double *v_plus,*v_minus,*dpot_dr_dr;                       /**/
  double orbitP, orbitN;       /* Inermost stable co and counter circular orbits (unitless?) */
  double xval, yval;
  double R_zero, R_two, R_four, R_six; /* Legendre Coefficients of Schwarzschild Radial Coord */
  double g_zero, g_two, g_four, g_six; /* Legendre Coefficients of the accel. due to gravity */
  int    signal;
} NeutronStar;





#endif
