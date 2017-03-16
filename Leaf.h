

double GetAmax( double sn, double sc );

// Collatz 1992 equation A12
// temperature in celcius
double Q10Func( double K25, double Q10, double T );

// Collatz 1991 equation A3
double FGammastar( double Oi, double tau );

// eq a2 collatz 1991
double FJeC3( double a, double alpha, double Q0, double Ci, double gammastar );

double FJeC4( double Q0, double a, double alpharf );

//eq a5 collatz 1991
double FJcC3( double Vm, double Ci, double gammastar, double Kc, double Ko, double Oi );

// Collatz 1992 eq5
double FJcC4( double Vm );

// eq A6 Collatz 1991
double FJsC3( double Vm );

// collatz 1992 eq4
double FJpC4( double Ci, double P, double k );


// used for both c3 and c4.
// follows woodward Ci is 0.7 of atm CO2 for C3
// for C4 we assume Ci is 5 times atm CO2 ???
// based on Collatz et al. 1991 equations but woodward concept
// for calculating Vmax
// double FVm( double An, double Rd, double Ca, double CiFrac, double gammastar, 
// 		   double Kc, double Ko, double Oi, double T );

double FVm( double An, double T );

// day respiration
// rate should be 0.015 for C3 and 0.025 for C4
// many sources, not sure who was first, see Arora
double RmL( double Vm, double rate );

// takes gs in mol/m2/s gives Ad in micromols/m2/s
// re-arrangement for equation 6B from collatz 1992
double CollatzAd( double gs, double Ci, double Ca, double P );

// collatz 1991 equation 1
// hs rh at leaf surface
// gives conductance in mol/m2/s
double CollatzBallBerryLeaf( double A, double Ca, double P, double m, double b, 
							double hs, double gb );

double GetAd( double m, double b, double Ab, double Rh, double P, 
			 double Ca, double Ci, double gb );

// gives aerodynamic conductance in mol/m2/s
// based on Jones 1992 empirical function
// coeff 0.271 assumes standard atm press and 25 deg C
double GetgbLEAF( double u, double cld );

double FindC3Ci( double Ci, double *aa );

double FindC4Ci( double Ci, double *aa );

//var lVm,lO2conc,ltau,lKc,lKo,laC4,lalpha,lQp,lCit,lP,lT,lws,lRh:real
double FindC3Min( double Vm, double Oi, double gammastar, double Kc, double Ko, double aC3, 
				  double mtree, double btree, double GammaC3, double alpha, double Q0, 
				  double P, double Ca, double T, double Rh, double gb );
 
//var lVm,lO2conc,ltau,lKc,lKo,laC4,lalpha,lQp,lCit,lP,lT,lws,lRh:real
double FindC4Min( double Vm, double aC4, double k, double mgrass, double bgrass, 
				  double GammaC4, double alpharf, double Q0, double P, double Ca, 
				  double T, double Rh, double gb);
  
// compute A0 and R for C3-plants
void GetC3A( double T, double sn, double sc, double Ca, double Oi, double Q0, 
			 double GammaC3, double aC3, double alpha,
			 double Rh, double P, double wind, double CLDtree,
			 double mtree, double btree, double tauC3K25, double tauC3Q10, 
			 double KoC3K25, double KoC3Q10, double KcC3K25, double KcC3Q10,
			 double *A, double *RmLv );

// compute A0 and R for C4-plants
void GetC4A( double T, double sn, double sc, double Ca, double k, 
			 double Q0, double GammaC4, double aC4,
			 double alpharf, double Rh, double P,
			 double wind, double CLDgrass, double mgrass, double bgrass,
 			 double *A, double *RmLv );

