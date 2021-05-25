/***************************************************************************
 * consts.h 
 *
 * This file contains definitions of constants used in the procedures.
 *
 **************************************************************************/

#define DM (1.0/(MDIV-1.0))          /* spacing in mu direction */ 
#define RDIV 900                     /* grid point in RK integration */ 
#define SMAX 0.9999                  /* maximum value of s-coordinate */  
#define DS (SMAX/(SDIV-1.0))         /* spacing in s-direction */
#define LMAX 10                      /* max. term in Legendre poly. */
#define C 2.99792458e10              /* speed of light in vacuum */
#define G 6.673e-8                   /* gravitational constant */ 
#define KAPPA (1.0e-15*C*C/G)        /* scaling factor */
#define KSCALE (KAPPA*G/(C*C*C*C))   /* another scaling factor */
#define MSUN 1.989e33                /* Mass of Sun */
#define SQ(x) ((x)*(x))              /* square macro */
#define MB 1.66e-24                  /* baryon mass */
#define RMIN 1.0e-15                 /* use approximate TOV equations when
                                        computing spherical star and r<RMIN */
#define UNUSED (-1.11e30)            /* Used in the Ridder zero-finder */ 
#ifndef PI
#define PI 3.141592653589793238462643              /* what else */
#endif

#define M_MAX 6         /* Maximum l=m mode considered */
#define BOUND 1         /* Excludes last grid point, since we set everything
                           to zero there */
#define MAXIT 30        /* Maximum number of iterations in rtsec_G */ 
#include <float.h>
#include <limits.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif



