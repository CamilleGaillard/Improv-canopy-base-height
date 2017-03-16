#ifndef PlantClassConstants_h___
#define PlantClassConstants_h___
#include <cmath>
#include "MyMath.h"

// double const CANOPY_AREA_0_INIT			= 0.0452;
// double const CANOPY_AREA_Z_INIT			= 0.0452;
// double const ROOT_AREA_0_INIT			= 0.0452;
// double const ROOT_AREA_D_INIT			= 0.0452;

double const MMS_TO_KGD_HELPER			= 44.*(12./44.)*1e-9*1.544;  // required to transform mumol/m^2/s carbon into kg/m^2/s biomass

double const MIN_ROOT_RADIUS			= 0.015;		// minimum radius of roots (1.5cm)
// double const MIN_ROOT_RADIUS			= 0.05;			// minimum radius of roots (5cm)

double const BETA_HELPER				= MIN_ROOT_RADIUS*MIN_ROOT_RADIUS*M_PI;

double const EVERGREEN_THR				= 0.5;			// at which value of "evergreen_" do we assume that plant is evergreen?
double const RAIN_LIGHT_THR				= 0.5;			// at which value of "rain_light_" do we assume that plants are summer green?

double const MAX_LAI					= 20.;

const int    NUM_LIGHT_NGB      = 8; // number of neighbours for light competition

//const double LIGHT_COMP_PAR     = 0.1;		// BEFORE LIAM CHANGE
// const double LIGHT_COMP_PAR     = 0.05;			// NEW LIAM
const double LIGHT_COMP_PAR     = 0.5;			// NEW LIAM

const double LIGHT_COMP_PAR_HLP = 1.-LIGHT_COMP_PAR;

const double LIGHT_COMP_SELF_DEAD = 0.1;            // SIMON self shading grasses

const int    N_PHEN             = 13;

//const int    MAX_SEED_NUM       = 25;

#endif





