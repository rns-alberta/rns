void surfacefit(EOS *eos, NeutronStar *star, double e_center, double ep_step, double r_step);

int SurfPrint(EOS *eos, NeutronStar *star);

void SurfGraph(EOS *eos, NeutronStar *star);

int Surface( EOS *eos, NeutronStar *star);

int LegnFit( EOS *eos, NeutronStar *star);

double SurfaceMetric(double *s_vec, double **f, double s_surf, int s, int m);

double  BodeInt(double *f1, double *f2);
