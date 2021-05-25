void orbit(EOS *eos, NeutronStar *star);

double OMEGA(EOS *eos, NeutronStar *star, int k, double s_stable, int sign);

double polation(double *yg, double *xg, double xs);

double pot_dr2P(EOS *eos, NeutronStar *star, int s, int sign, double *velcheck);
