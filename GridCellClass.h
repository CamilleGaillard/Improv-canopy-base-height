#ifndef GridCellClass_h___
#define GridCellClass_h___

#include "adgvm.h"
#include "GridCellClassStructures.h"
#include "GridCellClassConstants.h"
#include "PlantPopClass.h"
#include "SoilClass.h"		// SIMON soil carbon

class clGridCell
{
	private:
		int    runid_;
		int    fire_;
		double lat_;
		double lon_;
		double rd0_[MONTH_IN_YEAR];
		double pwet_[MONTH_IN_YEAR];
		double ralpha_[MONTH_IN_YEAR];
		double rbeta_[MONTH_IN_YEAR];


		//double wilt_point_;
		//double field_cap_;

		//other state variables
		GridCellSoilData soil[N_SOIL_LAYERS];				// array with soil properties
		GridCellClimateData clim[DAYS_IN_YEAR];				// array with climate properties
		GridCellPlantData plantdata;
		clPlantPop mypop_;
		
		clSoil soilcarbon_; 	// SIMON soil carbon		
		
		void getLeafPhotosyn();
		double WindAtHeight( double height, double ref_wind );
		// MP: new annual 365day storage of all relevant cell datas
		
	public:
		unsigned int latcount;
		unsigned int loncount;

		clGridCell();
		void initialize(int runid, int fire, MyInData *IData, unsigned int latcount, unsigned int loncount);
		void initializeyear( MyInData *IData);
		void RainFallYear();
		int runDailyProcesses( int year, int day /*, double par, double rad, double tmp, double hum, double wnd,
								double apr, double co2, double pre, double sun, double soc, double son, double *smo, double *root_bm */);
		int runAnnualProcesses( int year );
		void finalize();
		double getlat()			{ return lat_; }
		double getlon()			{ return lon_; }
		unsigned int getlatc()		{ return latcount; }
		unsigned int getlonc()		{ return loncount; }
		int    getrunid()		{ return runid_; }
		int    getfire()		{ return fire_; }
	
		GridCellSoilData* getsoil()	{ return soil; }
		GridCellClimateData* getclim(unsigned int day)	{ return &clim[day]; }
		GridCellPlantData* getplant()	{ return &plantdata; }
		clPlantPop* getplantpop() {return &mypop_;}
};
#endif

