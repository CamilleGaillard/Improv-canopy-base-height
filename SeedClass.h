#ifndef SeedClass_h___
#define SeedClass_h___
#include "SeedClassConstants.h"
class clGridCell;
class clSeed
{
    private:
		int    species_index_;
		int    veg_type_;      // woody plant or grass: 0=tree, 1=C4 grass
		int	   w_shrubs_ = 1;
		
// 		these traits are optimized
		double sla_;           // specific leaf area
		double alloc_root_;    // allocation to roots
		double alloc_leaf_;    // allocation to leaves
		double alloc_wood_;    // allocation to aboveground woody biomass; alloc_wood_ = alloc_stem_ + alloc_branch_ , when simulating with branches; otherwise it's equal to alloc_stem_!
		#ifdef W_BRANCHES
			double alloc_branch_;  // allocation to branches
		#endif
		double alloc_stem_;    // allocation to stem; 
		double alloc_bark_;    // allocation to bark
		double alloc_repr_;    // allocation to reproduction
		double alloc_stor_;    // allocation to storage
		
		double rain_light_;    // phenology control
		double evergreen_;     // evergreen plant?
		double rain_thr_on_;   // rain threshold on
		double rain_thr_off_;  // rain threshold off
		double light_thr_on_;  // light threshold on
		double light_thr_off_; // light threshold off
		
		double biomas_par_1_;       // height to biomass parameter
		double biomas_par_2_;       // height to biomass parameter
		double rotfrm_par_1_;       // root form parameter
		double rotfrm_par_2_;       // root form parameter
		double rotfrm_max_d_;       // maximum rooting depth
		
		double light_ext_;          // crown light extinction parameter
		double wood_density_;       // wood density, kg/m^3
		double seed_weight_;        // seed weight, in kg
		
		double canfrm_par_1_;       // crown form parameter
		double canfrm_par_2_;       // crown form parameter
		double canfrm_h_base_;      // crown height parameter - Cam

		double p_50_;             // p_50 value
		
		double stor_to_wood_;		// fraction of storage allocated to wood after fire
		double stor_to_stem_ ;		// fraction of storage allocated to stem (sub-fraction of stor_to_wood_; stor_to_stem_ + stor_to_stem_ = stor_to_wood_)
		double stor_to_branch_;		// fraction of storage allocated to branches	
		double stor_to_leaf_;		// fraction of storage allocated to leaf after fire
		
		double stem_count_;			// number of stems - NEW CAMILLE
		
// 		these traits are not optimized
		double mort_carbon_;        // mortality due to negative carbon balance
		double mort_background_;    // background mortaliy
		double mort_biomass_;       // mortality due to low biomass
		double mort_mechanic_1_;	// mortality due to mechanic instability
		double mort_mechanic_2_;	// mortality due to mechanic instability
		
		double topkill_const_;      // topkill parameter
		double topkill_h_;          // topkill parameter
		double topkill_i_;          // topkill parameter
		
		double ball_berry_m_;       // Parameter for Ball Berry equation (stomatal conductance)
		double ball_berry_b_;       // Parameter for Ball Berry equation (stomatal conductance)
		
		double r_maint_resp_;       // Maintainance respiration (fraction of A), required for stom. cond.
		double r_growth_resp_;      // Growth respiration (fraction of A)
		
		double beta_leaf_resp_;		// fraction of leaf biomass that respires
		double beta_wood_resp_;		// fraction of wood biomass that respires
		double beta_root_resp_;		// fraction of root biomass that respires
		double beta_n_resp_;		// parameter for respiration model

		double cn_ratio_leaf_;		// C:N raio of leaves
		double cn_ratio_wood_;		// C:N raio of woody biomass
		double cn_ratio_root_;		// C:N raio of roots
		
    public:
        clSeed();                       // standard constructor
        clSeed &operator = (const clSeed&orig_seed );
        void initialize( clGridCell *gridcell_); 
        int    getSpeciesIndex()    { return species_index_; }
        int    getVegType()         { return veg_type_; }
        double getEvergreen()       { return evergreen_; }
        double getSla()             { return sla_; }
        double getAllocRoot()       { return alloc_root_; }
        double getAllocLeaf()       { return alloc_leaf_; }
        double getAllocWood()       { return alloc_wood_; }
        double getAllocStem()		{ return alloc_stem_; }
        #ifdef W_BRANCHES
       		double getAllocBranch()		{ return alloc_branch_; }
       	#endif	
        double getAllocBark()       { return alloc_bark_; }
        double getAllocRepr()       { return alloc_repr_; }
        double getAllocStor()       { return alloc_stor_; }

        double getRainLight()       { return rain_light_; }
        double getRainThrOn()       { return rain_thr_on_; }
        double getRainThrOff()      { return rain_thr_off_; }
        double getLightThrOn()      { return light_thr_on_; }
        double getLightThrOff()     { return light_thr_off_; }
        double getCanfrmPar1()      { return canfrm_par_1_; }
        double getCanfrmPar2()      { return canfrm_par_2_; }
        double getCanfrmHBase()     { return canfrm_h_base_; } // crown height trait - Cam
        double getP50()		    { return p_50_; }
        double getRotfrmPar1()      { return rotfrm_par_1_; }
        double getRotfrmPar2()      { return rotfrm_par_2_; }
        double getRotfrmMaxD()      { return rotfrm_max_d_; }
        double getBiomasPar1()      { return biomas_par_1_; }
        double getBiomasPar2()      { return biomas_par_2_; }
        double getLightExt()        { return light_ext_; }
        double getWoodDensity()     { return wood_density_; }
        double getSeedWeight()      { return seed_weight_; }
        double getMortCarbon()      { return mort_carbon_; }
        double getMortBackground()  { return mort_background_; }
        double getMortBiomass()     { return mort_biomass_; }
        double getMortMechanic1()   { return mort_mechanic_1_; }
        double getMortMechanic2()   { return mort_mechanic_2_; }
        double getTopkillConst()    { return topkill_const_; }
        double getTopkillH()        { return topkill_h_; }
        double getTopkillI()        { return topkill_i_; }
	double getBallBerryM()		{ return ball_berry_m_; }
	double getBallBerryB()		{ return ball_berry_b_; }
	double getRMaintResp()		{ return r_maint_resp_; }
	double getRGrowthResp()		{ return r_growth_resp_; }
	double getBetaLeafResp()	{ return beta_leaf_resp_; }
	double getBetaWoodResp()	{ return beta_wood_resp_; }
	double getBetaRootResp()	{ return beta_root_resp_; }
	double getBetaNResp()		{ return beta_n_resp_; }
	double getCNRatioLeaf()		{ return cn_ratio_leaf_; }
	double getCNRatioWood()		{ return cn_ratio_wood_; }
	double getCNRatioRoot()		{ return cn_ratio_root_; }
		
	double getStorToWood()		{ return stor_to_wood_; }	// MP
	double getStorToStem()		{ return stor_to_stem_; }	// MP
	double getStorToBranch()	{ return stor_to_branch_; }	// MP
	double getStorToLeaf()		{ return stor_to_leaf_; }
	double getStemCount()		{ return stem_count_;	} // NEW CAMILLE
        
        void   setSpeciesIndex(   int    val ) { species_index_=val; }
        void   setCheckEvergreen( double val );
        void   setCheckSla(       double val );
        void   setCheckAllocRoot( double val );
        void   setCheckAllocLeaf( double val );
        void   setCheckAllocWood( double val );
        void   setCheckAllocBranch( double val, double val2 );
        void   setCheckAllocStem( double woody, double branch );
        void   setCheckAllocBark( double val );
        void   setCheckAllocRepr( double val );
        void   setCheckAllocStor( double val );
        void   setCheckRainLight( double val );
        void   setCheckRainThrOn( double val );
        void   setCheckRainThrOff( double val );
        void   setCheckLightThrOn( double val );
        void   setCheckLightThrOff( double val );
        void   setCheckBiomasPar1( double val );
        void   setCheckBiomasPar2( double val );
        void   setCheckRotfrmPar1( double val );
        void   setCheckRotfrmPar2( double val );
        void   setCheckRotfrmMaxD( double val );
        void   setCheckCanfrmPar1( double val );
        void   setCheckCanfrmPar2( double val );
        void   setCheckCanfrmHBase( double val );
        void   setCheckP50( double val );
        void   setCheckWoodDensity( double val );
        void   setCheckLightExt( double val );
	void   setCheckSeedWeight( double val );
	void   setCheckStorToWood( double val );
	void   setCheckStorToLeaf( double val );
	void   setCheckAllocParameters();   // this scales allocation parameters to sum=1
	void   setCheckStemCount( double val ); // NEW CAMILLE	
};

#endif


