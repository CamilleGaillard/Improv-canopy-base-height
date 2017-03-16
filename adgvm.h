#ifndef adgvm_h___
#define adgvm_h___
#include <libconfig.h++> 
// global variables
extern libconfig::Config cfg;
// global constants

//const int YEARS_TO_RUN = 501;

const int    DAYS_IN_YEAR  = 365;
const int    MONTH_IN_YEAR = 12;

const unsigned int END_OF_MONTH[MONTH_IN_YEAR] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};  // CAVE: these are the real days, count-days are -1 due to start at 0!!!

const unsigned int N_BUFF_YEARS = 25;
const unsigned int BUFF_SIZE = N_BUFF_YEARS*MONTH_IN_YEAR;
//const unsigned int trait_years = N_BUFF_YEARS / N_BUFF_YEARS;		// in this case once every 25 years
//const unsigned int POP_VEC_LENGTH = 31;							// has become obsolete since using map to collect output variables
//const unsigned int TRAIT_VEC_LENGTH = 48;							// has become obsolete since using map to collect output variables


const int IMISSING = -9999;
const float FMISSING = -9999.;
const int N_STEXTURE = 14;		// number of soil texture categories

// NOTE: RhoneAGG - table 4.3 -- very high ksat for sandy texture classes which are dubious
const double R[N_STEXTURE] = 			{0.,  12.13, 	10.21, 	10.21, 	8.16, 	8.16, 	3.63, 	5.42,  	9.39, 	6.10, 	6.77, 	4.73, 	4.32, 	3.91};
const double W_COND[N_STEXTURE] = 		{0., 0.004428, 0.003924,  0.003924, 0.00648, 0.007992, 0.00396, 0.03042, 0.008496, 0.027792, 0.022716, 0.11124, 0.381816, 0.871452};
const double RESID_CONTENT[N_STEXTURE] = {0., 0.020, 0.03, 0.041, 0.027, 0.015, 0.068, 0.075, 0.040, 0.109, 0.056, 0.090, 0.090, 0.090};
const double PHI_S[N_STEXTURE] = 		{0.0, 0.405, 0.49, 0.49, 0.356, 0.63, 0.786, 0.63, 0.153, 0.478, 0.356, 0.218, 0.09, 0.121};

/*

	r[j] = R[(int)soil[j].soil_texture]

*/

#endif
