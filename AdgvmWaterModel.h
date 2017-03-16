#ifndef AdgvmWaterModel_h___
#define AdgvmWaterModel_h___
#include "GridCellClass.h"
/*void BucketIn(  double drain, double *Theta, double ThetaFC, double lai, double ctr, double ksat, 
double soil_texture, double sat_water_cont, double perc_exponent );*/
void BucketIn( GridCellSoilData *soil, GridCellClimateData *clim,GridCellPlantData *plantdata );


// Input  : Et, Theta, ThetaWP
// Output : BOut
void BucketOut( GridCellSoilData *soil, GridCellPlantData *plantdata );

#endif
