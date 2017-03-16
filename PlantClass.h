#ifndef PlantClass_h___
#define PlantClass_h___
#include "SeedClass.h"
#include "PlantClassConstants.h"
#include "GridCellClassConstants.h"
#include "GridCellClassStructures.h"
class clGridCell;
class clPlant
{
	private:
		int    alive_;         // individual alive? required for reproduction/mortality
		int    age_;           // age of plant (years)
		int    dpos_;
		int    dneg_;
// 		int    c_def_counter_;
		int    dleaf_alloc_;   // counts number of days for which carbon is only allocated to leaves
		int    evg_flush_flag_;
		int    c_gain_cum_flag_;
		int    day_; 		//counts the number of days mortC
		int    id_;			// individual ID-number, allows to keep track of specific individuals
		int    dct_;		// day counter
		int    year_ ;		// year counter
// 		int    day1_; 		// counts the number of days Water Transport				
		
		double b_leaf_;        // state var: leaf biomass, live
		double b_leaf_d_st_;   // state var: leaf biomass, dead standing/hanging
		double b_leaf_d_ly_;   // state var: leaf biomass, dead lying
		double b_wood_;		   // state var: aboveground woody heartwood biomass
		#ifdef W_BRANCHES
			double b_branch_;	   // state var: branch heartwood biomass
		#endif	
		double bstem_old_ ; 	// FLAG: for testing, remove later
		double b_stem_;        // state var: stem heartwood biomass
		double b_root_;        // state var: coarse root biomass
		double b_bark_;        // biomass in bark/defence
		double b_repr_;        // biomass for reproduction
		double b_stor_;        // biomass in storage compartments for resprouting
		double light_ext_;	// SIMON
		double wood_density_;
		double sla_;
		double ldmc_;
		
		float soilcarbon_nwl_;  // SIMON non woody litter for soil C model
		float soilcarbon_fwl_;  // SIMON fine woody litter for soil C model
		float soilcarbon_cwl_;  // SIMON coarses woody litter for soil C model
				
		double leaf_thickness_grass_ = 0.2 ;
		double leaf_thickness_tree_	= 0.1 ;

		double root_fracs_[N_SOIL_LAYERS];  // fraction of root bimoass in different soil layers
//		double compheight_[NUM_LIGHT_NGB];   // heights of competitors	
//		double compheight_lai_[NUM_LIGHT_NGB];   // lai of competitors		// NEW LIAM
		double root_bm_[N_SOIL_LAYERS];		// root biomass for each soil layer
		
		double lai_;         //state var: leaf area index az z=0
		
		double height_;			// state var: plant height
		double stem_diam_int_;	// state var: stem diameter inner part, only b_stem_
		double stem_diam_tot_;	// state var: stem diameter total= stem_diam_int+bark_thick
		double bark_thick_;		// state var: bark thickness
		double wilt_point_;		// state var: wilt point
//		double wilt_avg_;
		double active_;			// state var: photosynthetically active
		double sum_active_ ;	// cumulative sum of active days within a year
		double sum_water_ ; 	// cumulative water availability 
		double sum_light_ ; 	// cumulative light availability
		double length_leaf_total_grass_;		// NEW LIAM
		
		double leaf_longevity_;
		double root_turnover_;		// NEW LIAM	
		double root_turnover_tree_ ;	// NEW LIAM	
		
		double crown_area_0_;    // crown diameter at z=0
		double crown_h_base_;    // SIMON crown base height, this is height where crown "starts"

		
		#ifdef W_BRANCHES
		
			void calCrownVolume() ;
			
			double rho_crown_;
			double crown_volume_;
			double crown_r_;		 // crown radius at base of the crown
			double crown_l_; 		 // length of the crown, from basis to top
			double par_d_ ; 		 // opening width of crown at top, for a given crown shape 	
			
			double b_leaf_old_ ;
			double b_branch_old_ ;
			double b_stem_old_ ;
			double b_bark_old_ ;
			double b_stor_old_ ;
			double b_repr_old_ ;
			double b_root_old_ ;
			
		#endif	

		double c_gain_;         // carbon gain
		double c_gain_mmolm2s_; // STEVE carbon gain in mmol per m2 per s
		double c_balance_;		// total carbon balance
		double c_gain_cum_;		// cumulative carbon gain since last leaf flush
		
		double stom_cond_;		// crown stomatal conductance
		double stom_cond_max_;		// crown stomatal conductance
		double boul_cond_;		// crown boundary layer conductance
		double evapotr_;		// evapotranspiration
		double evapotr_max_;		// evapotranspiration
		double leaf_temp_;		// leaf temperature
		
		double light_avail_;    // helper for light stressed photosynthesis
		double water_avail_;	// helper for water stressed photosynthesis
//		double theta_helper_;
		double water_t_max_;	//Liam--maximum daily water transport based on stem diameter Meinzer et al.(2005)
		double Gi_weighted_; 	//scaling plant water conductivity			// NEW LIAM
		double soil_mat_pot_weighted_;		// NEW LIAM
		

		double leaf_resp_helper_;	// required for leaf respiration calculations
		double wood_resp_helper_;	// required for wood respiration calculations
		double root_resp_helper_;	// required for root respiration calculations
		
		double can_frm_helper_;		// stores value 0.5^(1/CanfrmPar2)
		
		double phen_rain_[N_PHEN];	// stores raingreen trigger variable
		double phen_lght_[N_PHEN];	// stores summergreen trigger variable
		int    phen_loop_;
		int    stem_count_;		// NEW CAMILLE
		double area_sapwood_;	// NEW CAMILLE		
		double sap_frac_;   // fraction of stem area that is sapwood
// 		double RmSum_;

		double weight_ ;	// to create a palatability factor that will influence which grasses will get eaten preferentially
		double ind_consumption_ ; // biomass removed through grazing, per grass individual
		double loss_frac_ ;		// fraction of leaf biomass that got eaten
		double pot_realloc_ ;	// non-water limited reallocation from storage to leaves	
			
		
// 		double site_factor_;		// plant specific site factor, scales light availability
		
		clSeed sd_;
		
		void   calHeight();		// calculate stem height from biomass
		void   calStemDiam();	// calculate stem diameter from tree height
		void   calLai();
		void   calBoulCond( double wnd, double mms_helper );  // calculate stomatal conductance
		void   calStomCond( double tmp, double apr, double co2, double hum, double As, double As_max );  // calculate stomatal conductance
		void   calWaterAvailability( GridCellSoilData* soil, int year );
		void   calEvapotr(  double et_helper_0, double et_helper_1, double et_helper_2, double et_helper_3, double mms_helper  );
		void   calLeafTemp( double hum, double eS, double tmp, double apr, double mms_helper );
		void   calLightAvail();
		void   calWater_trans_max();	//Liam--maximum daily water transport based on stem diameter Meinzer et al.(2005)
		void   calCrownArea();    // calculate crown area at height z above ground
		void   calCrownBHeight();  // calculate crown base height from trait
		void   calLDMC();
		void   calSLA();
		void   calWoodDensity();
		void   calLeafLongevity();
		void   calStemCount(); // NEW CAMILLE
		void   calAreaSapwood(); // NEW CAMILLE		

		
	public:
		clPlant();
		void initialize( clGridCell *gridcell_, int id_in); 
		void   getRootBmVector( double *vec );
// 		void   runAnnualProcesses( double z_star, double thetaWP, double thetaFC );
// 		void   runAnnualProcesses( double ca_frac, double z_star );
		void   runAnnualProcesses(/* double *compheight */);
		void   runDailyProcesses( int year, int dct,GridCellSoilData* soil, double A0C3, double A0C4, double wnd, double pre,
								double tmp, double apr, double co2, double hum, double sun_dur, double rad, double resp_temp_fac,
								double et_helper_0, double et_helper_1, double et_helper_2, double et_helper_3, double eS,
								double mms_helper, double *compheight, double *compheight_lai, double *compheight_crown_area);
		int    runMortalityDaily();
		int    runMortalityAnnual();
		
		double runGrazing(double balance_, double demand_, double meristem_, int nind_grass_, double weight_coeff_) ; 
		double calWeightCoeff(double bm_slope_, double bm_intercept_, double a_, double b_, double c_) ;
		double cropInd(double demand_,double meristem_, int nind_grass_) ;
				
		void   runFireModel( double intensity, double flame_height );
		int    getSeedProd();
		void   calWiltPoint(GridCellSoilData* soil ); // calculate wilt point
		void   runDecomposition();

		void   setDead();
		double getBLeaf()           { return b_leaf_; }   // biomass (not C!)
		double getSla()           { return sla_; } 
		double getWoodDensity()		{ return wood_density_; }
		double getLightExt()		{ return light_ext_; }


		double getBWood()			{ return b_wood_; }
		#ifdef W_BRANCHES
			double getBBranch()			{ return b_branch_; }
		#endif	
		double getBStem()           { return b_stem_; }
		double getBRoot()           { return b_root_; }
		double getBBark()           { return b_bark_; }
		double getBRepr()           { return b_repr_; }
		double getBStor()           { return b_stor_; }
		double getBLeafDead()       { return b_leaf_d_st_ + b_leaf_d_ly_; }
		double getBLeafDeadSt()     { return b_leaf_d_st_; }
		double getBLeafDeadLy()     { return b_leaf_d_ly_; }
		
		float getSoilCarbonNWL(); // SIMON non woody litter for soil C model
		float getSoilCarbonFWL(); // SIMON fine woody litter for soil C model
		float getSoilCarbonCWL(); // SIMON coarses woody litter for soil C model		
	
		double getNLeaf()           { return b_leaf_*0.44/sd_.getCNRatioLeaf(); }   // nitrogen content, 0.44 for conversion to C
		double getNWood()           { return b_wood_*0.44/sd_.getCNRatioWood(); }
		#ifdef W_BRANCHES
			double getNBranch()			{ return b_branch_ *0.44/sd_.getCNRatioWood(); }
		#endif	
		double getNStem()			{ return b_stem_ *0.44/sd_.getCNRatioWood(); }
		double getNRoot()           { return b_root_*0.44/sd_.getCNRatioRoot(); }
		double getNBark()           { return b_bark_*0.44/sd_.getCNRatioWood(); }
		double getNRepr()           { return b_repr_*0.44/sd_.getCNRatioLeaf(); }
		double getNStor()           { return b_stor_*0.44/sd_.getCNRatioRoot(); }
		
		double getLai()             { return lai_; }
		double getHeight()          { return height_; }
		double getStemDiam()        { return stem_diam_tot_; }
		double getBarkThick()       { return bark_thick_; }
		double getWiltPoint()       { return wilt_point_; }
		double getActive()          { return active_; }
		double getCrownArea()      	{ return crown_area_0_; }
		double getCrown_h_base()   	{ return crown_h_base_; }
		double getStem_count()		{ return stem_count_; }
		double getCGain()           { return c_gain_; }
		double getCGain2()          { return c_gain_mmolm2s_; } //STEVE
		double getCBalance()        { return c_balance_; }
		double getCGainCum()        { return c_gain_cum_; }
		double getStomCond()        { return stom_cond_; }
		double getStomCond_max()    { return stom_cond_max_; }
		double getBoulCond()        { return boul_cond_; }
		double getEvapotr()         { return evapotr_; }
		double getEvapotr_max()     { return evapotr_max_; }
		double getLeafTemp()        { return leaf_temp_; }
		double getWater_trans_max() { return water_t_max_; }
		double getWaterAvail()      { return water_avail_; }		// NEW LIAM
		double getGi_avg()          { return Gi_weighted_; }		// NEW LIAM	

		double getCanFrmHelper()    { return can_frm_helper_; }
		double getRainLight()       { return sd_.getRainLight(); }
		int    getVegType()         { return sd_.getVegType(); }
		int    getAlive()           { return alive_; }
		int    getEvergreen()       { return sd_.getEvergreen()>EVERGREEN_THR ? 1. : 0.; }
		int    getAge()             { return age_;}							// NEW MIRJAM
		double getSeedWeight()      { return sd_.getSeedWeight(); }
		//double getSLA()             { return sd_.getSla(); }            // SIMON
		double getSLA()             { return sla_; }                      // SIMON
		double getLDMC()            { return ldmc_; }                     // SIMON
		double getDLeafAlloc()      { return dleaf_alloc_; }
   	 	int    getSpeciesIndex()    { return sd_.getSpeciesIndex(); }
		double getDay()             { return day_; }
		clSeed getSeed()            { return sd_; }
		clSeed* getsd()            { return &sd_; }

//		void   setSeed( clSeed sd );
		void   setInitialValues( clSeed sd, GridCellSoilData* soil );
};
#endif




