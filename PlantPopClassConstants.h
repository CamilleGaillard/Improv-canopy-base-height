/// -P------------------------
// asdf
#include "GridCellClassConstants.h"
#ifndef PlantPopClassConstants_h___
#define PlantPopClassConstants_h___

int const    MAX_POP_SIZE_SQRT  = 60;   // square root of pop size
int const    MAX_POP_SIZE       = MAX_POP_SIZE_SQRT*MAX_POP_SIZE_SQRT;


double const PLOT_SIZE          = 100.*100.;   // in m^2, required for z* calculations

double const SIZE_FACTOR        = 10000./PLOT_SIZE; // for calculations of (unit)/ha

double const IND_AREA_MULT 		= PLOT_SIZE / (double(MAX_POP_SIZE) * 0.9971082773)	;  // in  m^2, subgrid area available per individual based on current max. canopy_area_0 of grasses

const double SP_HEAT            = 1.013e-3;   	// MJ/kg/degC = J/g/degC from FOA specific heat of moist air}


#endif









