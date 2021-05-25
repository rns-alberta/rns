								    
void print(double r_ratio,
	   double e_center, double Mass, double Mass_0, double R_e,
	   double Omega, double Omega_K, double J
	   );

int MakeSphere(EOS *eos,
	       NeutronStar *star,
	       double e_center);

double Kepler(EOS *eos, NeutronStar *star,
	       double e_center);

int rns(double r_ratio, 
	double e_center,
	EOS *eos, 
	NeutronStar *star
	);

int SetUpStar( char eos_file[80],
	       char eos_type[10],
	       char data_dir[80],
	       double Gamma_P,
	       double B,
	       double K,
	       EOS *eos,
	       NeutronStar *star);

int SetSpin( EOS *eos, NeutronStar *star,
	     double e_center, double spinfreq);

int SetJ( EOS *eos, NeutronStar *star,
	  double e_center, double J);

int MaxMass(double emin, double emax,
	    EOS *eos, NeutronStar *star, char name[80]);

int SetMassRatio(double Ratio,
                double Mass,
                double e_center,
                EOS *eos,
                NeutronStar *star,
		 int printflag);


double printpolint(double *xp, double *yp, int order, double xb,
                   double *err, int printflag);
