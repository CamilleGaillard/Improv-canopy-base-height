
//-----------------------------------------------------------------------------

// psychromatic constant
double FOAgama(double P);


// slope of curve, expressing the relation between vapor saturation and temp.
double FOAs(double T);

// air density
double FOArho(double P, double T);

double FOAe0(double T);

double FOAeS(double Tmax, double Tmin);

double FOAVPD(double eA, double eS);

// Compute the evapotranspiration
double FOAPenMan(double Rn, double ga, double gn, double s, double cp, double gama, double rho, double VPD);

// Compute the reference evapotranspiration
double FOAPenManRef( double temp_actual, double temp_last, double delta, double radiation, 
					 double gamma, double wind, double e_mean, double e_actual );

