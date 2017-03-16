#ifndef LeafGlobals_h___
#define LeafGlobals_h___
#include <cmath> 

// Constants for leaf photosynthesis

const double ABS_PHOTONS_C3			= 0.86;			// absorbtance to incident flux of photons, Collatz C3
const double ALPHA_C3				= 0.08;			// intrinsic quantum efficiency for C02 uptake, Collatz C3
const double ALPHAR_F_C4			= 0.067;		// product of alpha_r (intrinsic quantum yield of C3
const double ABS_PHOTONS_C4			= 0.80;			// absorbtance to incident flux of photons, Collatz C4
													// photosynthesis) and f (fraction of absorbed photons used
													// by C3 reactions), FJiC4, Collatz C4 units mol/mol
const double KAPPA_C4				= 0.7*1e6;		// initial slope of photosynthetic CO2 resp for C4, 0.7 mumol/m2/s
													// FC4Jc, Collatz C4
const double OI_PAR_PREASSURE		= 21.0;			// O2 Partial pressure KPa, Collatz

const double B_TREE_C3				= 0.01;			// ball berry constant conductance units (take care with reporting units here)
const double M_TREE_C3				= 9.0;			// ball berry constant dimensionless
const double B_GRASS_C4				= 0.04;			// ball berry constant conductance units
const double M_GRASS_C4				= 4.0;			// ball berry constant dimensionless

const double KC_K25_C3				= 30.0;			// c3 value of Kc at 25degC, Pa, Collatz C3
const double KC_Q10_C3				= 2.1;			// c3 rate of change of Kc with temp, Collatz C3
const double KO_K25_C3				= 30.0;			// c3 value of Ko at 25degC, KPa, Collatz C3
const double KO_Q10_C3				= 1.2;			// c3 rate of change of Ko with temp, Collatz C3
const double TAU_K25_C3				= 2600.0;		// c3 value of tau at 25degC,  Collatz C3
const double TAU_Q10_C3				= 0.57;			// c3 rate of change of tau with temp, Collatz C3

const double CS_FACTOR				= 1.4;

const double CLD_TREE				= 0.02;			// charac leaf dimension tree (assumes 2 cm leaf width)
const double CLD_GRASS				= 0.005;		// charac leaf dimension tree (assumes 2 cm leaf width)

const double ZD_CONST				= 0.86;			// constant for estimating displacement height, Jones 1992
const double Z0_CONST				= 0.06;			// constant for estimating roughness length, Jones 1992
const double VEGETATION_HEIGHT		= 1.5;			// H_bar average vegetation height in m (needed to calculate wind)
const double ROUGHNESS_LENGTH		= Z0_CONST*VEGETATION_HEIGHT;
const double DISPLACEMENT_HEIGHT	= ZD_CONST*VEGETATION_HEIGHT;
const double REF_HEIGHT_Z			= 10.;			// Height, at which wind is meassured
const double WIND_AT_HEIGHT_HELPER	= 1./log( (REF_HEIGHT_Z-DISPLACEMENT_HEIGHT)/ROUGHNESS_LENGTH );

const double KARMAN_CONST			= 0.41;			// von Karman constant, for canopy boundary layer conductance, no dimension
const double KARMAN_CONST_QUAD		= pow(KARMAN_CONST,2.);	// quadrat of von Karman constant

const double R_MAINT_RESP_C4		= 0.025;		// C3 RmL leaf  maint respiration as proportion of Vm, value from Collatz  original
const double R_MAINT_RESP_C3		= 0.015;		// C4 RmL leaf  maint respiration as proportion of Vm, value from Collatz  original

const double SGC					= 0.287;		// specific gas constant, KJ/kg/K, FOArho}
const double LAMBDA					= 2.45;			// MJ/kg from FOA latent heat of air at 20degC}
const double GSC					= 0.082;		// solar constant - extrateresstrial solar radiation MJ/m2/minute
const double ANGSTRONG_A			= 0.25;			// constant defining prop of extrat. radiation reaching earth on overcast days
const double ANGSTRONG_B			= 0.50;			// constant defining the additional prop of extrat. radiation reaching earth on clear days
const double ALBEDO					= 0.23;			// canopy reflection coefficient estimate for a grass reference crop (dimensionless)
const double SBC					= 4.903e-09;	// stephan bolzman constant MJ.k^-4.day^-1  

#endif























