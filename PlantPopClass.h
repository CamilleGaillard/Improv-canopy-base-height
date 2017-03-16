#ifndef PlantPopClass_h___
#define PlantPopClass_h___
#include <netcdf.h>
//#include <fstream>   // Simon
#include <map>
#include "PlantPopClassConstants.h"
#include "PlantClass.h"

#include "NcInputClasses.h"
class clGridCell;
#include "SeedBankClass.h"
class clPlantPop
{
	private:
		clPlant plant_pop_[MAX_POP_SIZE];
		clSeedBank seed_bank_[VTYPE_NUM];
		float aux_[3]; //note this is just a temporary storage not to be used for anything else;
		
		int    runid_;
		int    fire_;
		int    firecount_;
		int YEARS_TO_RUN ;
		int year_ ;
		int dct_ ;
		
//		int    pop_size_;
		int pop_size_[3];		// [0]: size of entire population; [1]: size of tree population; [2]; size of C4-grass population; (nunber of alive individuals at a given point in time)
		float lat_;
		float lon_;
				
//		float npp_annual_;
		
		float soilcarbon_nwl_;  // SIMON non woody litter for soil C model
		float soilcarbon_fwl_;  // SIMON fine woody litter for soil C model
		float soilcarbon_cwl_;  // SIMON coarses woody litter for soil C model		
		
		float root_bm_[N_SOIL_LAYERS];
// 		float n_total_;			// total N in vegetation, g/m^2
// 		float n_uptake_annual_;	// N uptake, g/m^2/year
		
		//float b_leaf_dead_;     // dead leaf biomass pool
		
		double crown_area_0_ ;
		
		float getCanopyArea();
		
		void   runSeedProd();
		void   runMortalityDaily();
		void   runMortalityAnnual();
		void   runRecruitment(clGridCell *gridcell);
		
		//float getBLeafDead() { return b_leaf_dead_*SIZE_FACTOR/1000.; }
		
		float getMeanStemDiam(float *tmp);
		float getMeanWiltPoint(float *tmp);
		float getMeanActive(float *tmp);
		float getMeanVegType( int type );
		float getMeanCanopyArea0(float *tmp);
		float getMeanCGain(float *tmp);
		float getMeanCGain2(float *tmp); //STEVE
		float getCGain_Plot(float *tmp); //STEVE
		float getMeanRainLight(float *tmp);
		float getMeanAlive(float *tmp);
		
		float getSumBLeaf(float *tmp);
//         float getSumBLeafDead();
//         float getSumBLeafDead( return b_leaf_dead_ );
		float getSumBStem(float *tmp);
		float getSumBWood(float *tmp);
		float getSumBBranch(float *tmp);				
		float getSumBRoot(float *tmp);
		float getSumBBark(float *tmp);
		float getSumBRepr(float *tmp);
		float getSumBStor(float *tmp);
		
		float getDryFuel( float *tmp);
		float getWetFuel();
		
// 		float getSumN()  { return n_total_; }
		
//         float getSumLai();
//         float getSumHeight();
		float getSumBasalArea(float *tmp);
//         float getSumWiltPoint();
		float getSumActive(float *tmp);
		float getSumCGain(float *tmp);
		float getSumCGainCum(float *tmp);
		float getSumRainLight(float *tmp);
		float getSumAlive(float *tmp);


	// grazing-related local variables		FLAG: CHECK AFTERWARDS WHICH OF THOSE ARE ACTUALLY STILL USED!
		int with_grazing_ ;
		float meristem_ ;
		float cows_per_ha_ ;
		float demand_cow_ ;
		float demand_ ;
		float deficit_ ;	// negative: cows' demand exceeds supply
		float gross_supply_ ;
		float net_supply_ ;
		float mean_avail_ind_ ;
		float balance_ ;
		float frequency_; 
    	double bm_slope_ ;
	 	double bm_intercept_ ;
	 	double a_ ;
		double b_ ;
		double c_ ;
		double weight_coeff_[MAX_POP_SIZE] ;		    
	float* data_pointer(std::string varname, std::string longname, std::string units, unsigned int dimlength, nc_type type); //this is for pop_data	
	float* data_pointer(std::string varname, std::string longname, std::string units, nc_type type); //this is for ind_data
    public:
#ifdef NC_OUTPUT
		int    buffidx;
		std::map<std::string,vector<float>> ind_data;
		std::map<std::string,vector<float>> pop_data;

		/*MLL_array* trait_matrix;
		MLL_array* pop_matrix;													
		MLL_array* soil_moist;	*/
#endif

		clPlantPop();
//		clPlantPop( float thetawp, float thetafc );
		//void initialize( int runid, int fire, float lat, float lon,  GridCellSoilData* soil, unsigned int latcount_in, unsigned int loncount_in );
		void initialize( clGridCell *gridcell_);
		void finalize();
		void writePlantSmall( int i );   // write data from plant i in population
		void write(GridCellSoilData* soil);
		int runDailyProcesses( clGridCell *gridcell, int year, int day, float A0C3, float A0C4,
								float wnd, float pre, float tmp, float apr,
								float co2, float hum, float sun_dur, float rad );
		int runAnnualProcesses( int year, clGridCell *gridcell );
		void runGrazing();
		void bm_weight_coeffs();
		void sla_weight_coeffs();		
		void runFireModel( int year, int day, float wnd, float hum, float pre, GridCellSoilData* soil );
		float getTreeAlive( int num )  { return plant_pop_[num].getAlive(); }
		float getMeanLai(float *tmp);
// 		float getMeanLaiTree();
// 		float getMeanLaiNew();
		float getSumStomCond(float *tmp);
		float getMeanStomCond(float *tmp);
		float getMeanAge(float *tmp);				// NEW MIRJAM
		float getTreeFraction();
		float getC4GrassFraction();
		float getC3GrassFraction();
		float getMeanHeight(float *tmp);
		float getSumCanopyArea0(float *tmp);
		float getSumEvapotr(float *tmp);
		float getMeanWaterAvail(float *tmp);		// NEW LIAM	
		float getMeanGi(float *tmp);				// NEW LIAM
		float getNAliveByType(float *tmp);			// NEW MIRJAM
		float getShoot2Root(float *tmp);			// NEW MIRJAM
		float getMeanSla(float *tmp);				// NEW MIRJAM
		float getMeancanopy_h_base(float *tmp);		// NEW CAMILLE		
		float getMeanStem (float *tmp);				// NEW CAMILLE

		float getSoilCarbonNWL() { float tmp = soilcarbon_nwl_; soilcarbon_nwl_ = 0.; return tmp; } // SIMON non woody litter for soil C model
		float getSoilCarbonFWL() { float tmp = soilcarbon_fwl_; soilcarbon_fwl_ = 0.; return tmp; } // SIMON fine woody litter for soil C model
		float getSoilCarbonCWL() { float tmp = soilcarbon_cwl_; soilcarbon_cwl_ = 0.; return tmp; } // SIMON coarses woody litter for soil C model		
		
};




#endif
