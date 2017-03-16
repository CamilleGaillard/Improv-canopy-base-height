/// Seed class
/// ----- Seed class

#include <iostream>
#include <iomanip>
#include <cmath>

#include "SeedClass.h"

#ifdef D_ONE_TRAIT_COMB
	#include "SeedClassConstants_one_comb_only.h"
#else
	#include "SeedClassConstants.h"
#endif
#include "MyMath.h"
#include "GridCellClass.h"

extern double global_max_root_depth[2];   // FLAG GLOBAL VARIABLE WORKAROUND SHOULD BE REMOVED IN FUTURE!!!

using namespace std;

//------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

clSeed::clSeed()//constructor
{
}

void clSeed::initialize( clGridCell * gridcell_)
{
	w_shrubs_			= cfg.lookup("simulation.W_SHRUBS");
// 	cout << "Random traits" << "  species_index_  " << species_index_ << endl;
	veg_type_           = MyRound( RUnif() );		// BEFORE LIAM CHANGE
//	veg_type_			= GR_4 ;					// allow only grasses to establish
//	veg_type_			= TREE ;					// allow only trees to establish
	species_index_      = (int)floor( RUnif()*SPECIES_NUM); // SIMON set to SPECIES_NUM again
	p_50_               = RUnif( P_50_MIN[veg_type_],  P_50_MAX[veg_type_] );	// NEW LIAM, NOW UP HERE
//	sla_                = RUnif( SLA_MIN[veg_type_],         SLA_MAX[veg_type_] );		// BEFORE LIAM CHANGE
	sla_                = 0.; //1/((0.224+(-0.41*p_50_))*0.1);  // NEW LIAM  SIMON set to 0 because function of p50 now
	alloc_wood_         = RUnif( ALLOC_WOOD_MIN[veg_type_],  ALLOC_WOOD_MAX[veg_type_] );
	#ifdef W_BRANCHES
		alloc_branch_		= RUnif( ALLOC_BRANCH_MIN[veg_type_], ALLOC_BRANCH_MAX[veg_type_] ) * alloc_wood_;
		alloc_stem_			= alloc_wood_ - alloc_branch_ ;
	#else
		alloc_stem_			= alloc_wood_ ;
	#endif	
//	wood_density_       = RUnif( WOOD_DENSITY_MIN[veg_type_], WOOD_DENSITY_MAX[veg_type_] );		// BEFORE LIAM CHANGE
	wood_density_       = 0.; // ((0.259+(-0.05921*p_50_))*1000); // NEW LIAM   SIMON set to 0 because function of p50 now
	seed_weight_        = RUnif( SEED_WEIGHT_MIN[veg_type_], SEED_WEIGHT_MAX[veg_type_] );
	rain_light_         = RUnif( RAIN_LIGHT_MIN[veg_type_],  RAIN_LIGHT_MAX[veg_type_] );
	evergreen_          = RUnif( EVERGREEN_MIN[veg_type_],     EVERGREEN_MAX[veg_type_] );
	alloc_root_         = RUnif( ALLOC_ROOT_MIN[veg_type_],    ALLOC_ROOT_MAX[veg_type_] );
	alloc_leaf_         = RUnif( ALLOC_LEAF_MIN[veg_type_],    ALLOC_LEAF_MAX[veg_type_] );
	alloc_bark_         = RUnif( ALLOC_BARK_MIN[veg_type_],    ALLOC_BARK_MAX[veg_type_] );
	alloc_repr_         = RUnif( ALLOC_REPR_MIN[veg_type_],    ALLOC_REPR_MAX[veg_type_] );
	alloc_stor_         = RUnif( ALLOC_STOR_MIN[veg_type_],    ALLOC_STOR_MAX[veg_type_] );
//	rain_thr_on_        = RUnif( THR_RAIN_ON_MIN[veg_type_],   THR_RAIN_ON_MAX[veg_type_] );
//	rain_thr_off_       = RUnif( THR_RAIN_OFF_MIN[veg_type_],  THR_RAIN_OFF_MAX[veg_type_] );
	rain_thr_on_		= RUnif( p_50_, 0.) ;
	rain_thr_off_       = RUnif( THR_RAIN_OFF_MIN[veg_type_], THR_RAIN_OFF_MAX[veg_type_] );    // MP: this is the value that needs to be subtracted from the ON-threshold;
	rain_thr_off_       = MyMax( rain_thr_on_ - rain_thr_off_, P_50_MIN[veg_type_] ) ;     // MP: now calculate the actual OFF-threshold's value, and make sure it's never lower than the lowest p50 of -3.0!
//	rain_thr_off_		= MyMax( rain_thr_on_ - rain_thr_off_, THR_RAIN_ON_MIN[veg_type_] ); 	// MIRJAM: OFF-threshold always <= ON-threshold, but never lower than the lowest allowed ON-threshold
	light_thr_on_       = RUnif( THR_LIGHT_ON_MIN[veg_type_],  THR_LIGHT_ON_MAX[veg_type_] );
	light_thr_off_      = RUnif( THR_LIGHT_OFF_MIN[veg_type_], THR_LIGHT_OFF_MAX[veg_type_] );
	light_thr_off_		= MyMax( light_thr_on_ - light_thr_off_, THR_LIGHT_ON_MIN[veg_type_]); 	// MIRJAM: OFF-threshold always <= ON-threshold, but never lower than the lowest allowed ON-threshold
	biomas_par_1_       = RUnif( BIOMAS_PAR_1_MIN[veg_type_],  BIOMAS_PAR_1_MAX[veg_type_] );
	biomas_par_2_       = RUnif( BIOMAS_PAR_2_MIN[veg_type_],  BIOMAS_PAR_2_MAX[veg_type_] );
	rotfrm_par_1_       = RUnif( ROTFRM_PAR_1_MIN[veg_type_],  ROTFRM_PAR_1_MAX[veg_type_] );
	rotfrm_par_2_       = RUnif( ROTFRM_PAR_2_MIN[veg_type_],  ROTFRM_PAR_2_MAX[veg_type_] );
//	rotfrm_max_d_       = RUnif( ROTFRM_MAX_D_MIN[veg_type_],  ROTFRM_MAX_D_MAX[veg_type_] );	// BEFORE LIAM CHANGE
// 	rotfrm_max_d_       = RUnif( global_max_root_depth[veg_type_],  global_max_root_depth[veg_type_] );

	//cout << "global_max_root_depth: "<< gridcell_->getplant()->max_root_depth[0]<<" "<<gridcell_->getplant()->max_root_depth[1]<<endl;
//	cout << "woody allocation: " << setw(14) << alloc_wood_ << setw(14) << alloc_branch_ << setw(14) << alloc_stem_ << setw(14) << alloc_branch_ + alloc_stem_ << endl ;
	
	if (gridcell_->getplant()->max_root_depth[veg_type_] > ROTFRM_MAX_D_MIN[veg_type_])
	{
		rotfrm_max_d_       = RUnif( ROTFRM_MAX_D_MIN[veg_type_],   MyMin( ROTFRM_MAX_D_MAX[veg_type_], gridcell_->getplant()->max_root_depth[veg_type_]) );		// NEW LIAM, FLAG!!!! NEED TO SOMEHOW PARSE THIS TRHOUGH ALL THE WAY FROM NcInput!!!
	}
	else 
	{
		rotfrm_max_d_ = ROTFRM_MAX_D_MIN[veg_type_];
	//	cout << "warning global_max_root_depth[veg_type_] <=ROTFRM_MAX_D_MIN[veg_type_]"<<endl;
	}
	
	if(veg_type_==GR_4 && rotfrm_max_d_ > 0.6) cout << "SeedClass: rotfrm_max_d GR_4 too high: " << rotfrm_max_d_ << endl;			// FLAG: only temporary to check, remove later! MP

	light_ext_          = RUnif( LIGHT_EXT_MIN[veg_type_],     LIGHT_EXT_MAX[veg_type_] );
//	wood_density_       = RUnif( WOOD_DENSITY_MIN[veg_type_],  WOOD_DENSITY_MAX[veg_type_] );	// BEFORE LIAM CHANGE
	canfrm_par_1_       = RUnif( CANFRM_PAR_1_MIN[veg_type_],  CANFRM_PAR_1_MAX[veg_type_] );
	canfrm_par_2_       = RUnif( CANFRM_PAR_2_MIN[veg_type_],  CANFRM_PAR_2_MAX[veg_type_] );
	canfrm_h_base_      = RUnif( CANFRM_H_BASE_MIN[veg_type_], CANFRM_H_BASE_MAX[veg_type_] );
	stor_to_wood_       = RUnif( STOR_TO_WOOD_MIN[veg_type_],  STOR_TO_WOOD_MAX[veg_type_] );
	stor_to_leaf_       = RUnif( STOR_TO_LEAF_MIN[veg_type_],  STOR_TO_LEAF_MAX[veg_type_] );
	
	#ifdef W_BRANCHES
		if(alloc_wood_ > 0.)	// avoid div0 here for grasses, which have no woody allocation!
		{
			stor_to_stem_		= stor_to_wood_ * alloc_stem_/alloc_wood_ ;		// MP, for simpliticty for now let's assume the ratio for reallocation to stem and branches is the same as the regular allocation ratio
			stor_to_branch_		= stor_to_wood_ - stor_to_stem_ ;
		}
		else
		{
			stor_to_stem_ = 0. ;
			stor_to_branch_ = 0. ;
		}	
	#else
		stor_to_stem_		= stor_to_wood_ ;
	#endif	
	
	if (w_shrubs_==1) stem_count_ = RUnif( STEM_COUNT_MIN[veg_type_],    STEM_COUNT_MAX[veg_type_] );	// NEW CAMILLE
	else stem_count_ = 1.;
				
// 	cout << "Random traits" << "  species_index_  " << species_index_ << endl;
	if ( species_index_<0 ) cout << "random" << endl;

	
	mort_carbon_        = MORT_CARBON[veg_type_];
	mort_background_    = MORT_BACKGROUND[veg_type_];
	mort_biomass_       = MORT_BIOMASS[veg_type_];
	mort_mechanic_1_    = MORT_MECHANIC_1[veg_type_];
	mort_mechanic_2_    = MORT_MECHANIC_2[veg_type_];
	
	topkill_const_      = TOPKILL_CONST[veg_type_];
	topkill_h_          = TOPKILL_H[veg_type_];
	topkill_i_          = TOPKILL_I[veg_type_];
	ball_berry_m_       = BALL_BERRY_M[veg_type_];
	ball_berry_b_       = BALL_BERRY_B[veg_type_];
	r_maint_resp_       = R_MAINT_RESP[veg_type_];
	r_growth_resp_      = R_GROWTH_RESP[veg_type_];
	beta_leaf_resp_     = BETA_LEAF_RESP[veg_type_];
	beta_wood_resp_     = BETA_WOOD_RESP[veg_type_];
	beta_root_resp_     = BETA_ROOT_RESP[veg_type_];
	beta_n_resp_        = BETA_N_RESP[veg_type_];
	cn_ratio_leaf_      = CN_RATIO_LEAF[veg_type_];
	cn_ratio_wood_      = CN_RATIO_WOOD[veg_type_];
	cn_ratio_root_      = CN_RATIO_ROOT[veg_type_];
}


//------------------------------------------------------------------------------------------------------------

clSeed& clSeed::operator = (const clSeed& orig_seed )
{
//     cout << " = operator " << endl;
    this -> species_index_      = orig_seed.species_index_;
    this -> veg_type_           = orig_seed.veg_type_;
    this -> evergreen_          = orig_seed.evergreen_;
    this -> sla_                = orig_seed.sla_;
    this -> alloc_root_         = orig_seed.alloc_root_;
    this -> alloc_leaf_         = orig_seed.alloc_leaf_;
    this -> alloc_wood_         = orig_seed.alloc_wood_;
    # ifdef W_BRANCHES
    	this -> alloc_branch_		= orig_seed.alloc_branch_;
    #endif	
    this -> alloc_stem_			= orig_seed.alloc_stem_;
    this -> alloc_bark_         = orig_seed.alloc_bark_;
    this -> alloc_repr_         = orig_seed.alloc_repr_;
    this -> alloc_stor_         = orig_seed.alloc_stor_;
    this -> rain_light_         = orig_seed.rain_light_;
    this -> rain_thr_on_        = orig_seed.rain_thr_on_;
    this -> rain_thr_off_       = orig_seed.rain_thr_off_;
    this -> light_thr_on_       = orig_seed.light_thr_on_;
    this -> light_thr_off_      = orig_seed.light_thr_off_;
    this -> canfrm_par_1_       = orig_seed.canfrm_par_1_;
    this -> canfrm_par_2_       = orig_seed.canfrm_par_2_;
    this -> canfrm_h_base_      = orig_seed.canfrm_h_base_;
    this -> p_50_               = orig_seed.p_50_;
    this -> rotfrm_par_1_       = orig_seed.rotfrm_par_1_;
    this -> rotfrm_par_2_       = orig_seed.rotfrm_par_2_;
    this -> rotfrm_max_d_       = orig_seed.rotfrm_max_d_;
    this -> biomas_par_1_       = orig_seed.biomas_par_1_;
    this -> biomas_par_2_       = orig_seed.biomas_par_2_;
    this -> light_ext_          = orig_seed.light_ext_;
    this -> wood_density_       = orig_seed.wood_density_;
    this -> seed_weight_        = orig_seed.seed_weight_;
    this -> mort_carbon_        = orig_seed.mort_carbon_;
    this -> mort_background_    = orig_seed.mort_background_;
    this -> mort_biomass_       = orig_seed.mort_biomass_;
    this -> mort_mechanic_1_    = orig_seed.mort_mechanic_1_;
    this -> mort_mechanic_2_    = orig_seed.mort_mechanic_2_;
    this -> topkill_const_      = orig_seed.topkill_const_;
    this -> topkill_h_          = orig_seed.topkill_h_;
    this -> topkill_i_          = orig_seed.topkill_i_;
    this -> ball_berry_m_       = orig_seed.ball_berry_m_;
    this -> ball_berry_b_       = orig_seed.ball_berry_b_;
    this -> r_maint_resp_       = orig_seed.r_maint_resp_;
    this -> r_growth_resp_      = orig_seed.r_growth_resp_;
    this -> beta_leaf_resp_     = orig_seed.beta_leaf_resp_;
    this -> beta_wood_resp_     = orig_seed.beta_wood_resp_;
    this -> beta_root_resp_     = orig_seed.beta_root_resp_;
    this -> beta_n_resp_        = orig_seed.beta_n_resp_;
    this -> cn_ratio_leaf_      = orig_seed.cn_ratio_leaf_;
    this -> cn_ratio_wood_      = orig_seed.cn_ratio_wood_;
    this -> cn_ratio_root_      = orig_seed.cn_ratio_root_;
    this -> stor_to_wood_       = orig_seed.stor_to_wood_;		// MP
    this -> stor_to_stem_		= orig_seed.stor_to_stem_;		// MP
    #ifdef W_BRANCHES
	    this -> stor_to_branch_		= orig_seed.stor_to_branch_;	// MP
	#endif   
    this -> stor_to_leaf_       = orig_seed.stor_to_leaf_;
    this -> stem_count_			= orig_seed.stem_count_; 	// NEW CAMILLE    
// 	cout << "SPINI " << 1000+species_index_ << endl;
	
    return *this;
}

//---------------------------------------------------------------------------

// setCheck functions ensure that traits remain between the minimum and maximum
// values as defined in SeedClassConstants.h after the optimization step.

// functions are called in SeedBankClass runMutation, and runCrossover, when seeds with new trait combinations are created

void clSeed::setCheckEvergreen( double val ) {
    evergreen_     = MyMax( EVERGREEN_MIN[veg_type_], MyMin( EVERGREEN_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckSla( double val ) {
    sla_           = MyMax( SLA_MIN[veg_type_], MyMin( SLA_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckAllocRoot( double val ) {
    alloc_root_    = MyMax( ALLOC_ROOT_MIN[veg_type_], MyMin( ALLOC_ROOT_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckAllocLeaf( double val ) {
    alloc_leaf_    = MyMax( ALLOC_LEAF_MIN[veg_type_], MyMin( ALLOC_LEAF_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckAllocWood( double val ) {
    alloc_wood_    = MyMax( ALLOC_WOOD_MIN[veg_type_], MyMin( ALLOC_WOOD_MAX[veg_type_], val ) );
	return;
}

#ifdef W_BRANCHES
	// FLAG MP: not sure about this function's correctness; multiplication-factor alloc_wood_? alloc_branch_ is supposed to be a trait-defined fraction of alloc_wood_ and
	// alloc_stem_ the remainder of the woody allocation...

	void clSeed::setCheckAllocBranch (double val, double val2 ) {
   	 alloc_branch_    = MyMax( ALLOC_BRANCH_MIN[veg_type_], MyMin( ALLOC_BRANCH_MAX[veg_type_], val ) ) * val2 ;
		return;
	}

	// FLAG MP: not sure I need this, but think so; double-check!

	void clSeed::setCheckAllocStem(double woody, double branch) {
		alloc_stem_		= woody - branch ;
		return;
	}
#endif	

void clSeed::setCheckAllocBark( double val ) {
    alloc_bark_    = MyMax( ALLOC_BARK_MIN[veg_type_], MyMin( ALLOC_BARK_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckAllocStor( double val ) {
    alloc_stor_    = MyMax( ALLOC_STOR_MIN[veg_type_], MyMin( ALLOC_STOR_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckAllocRepr( double val ) {
    alloc_repr_    = MyMax( ALLOC_REPR_MIN[veg_type_], MyMin( ALLOC_REPR_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRainLight( double val ) {
    rain_light_    = MyMax( RAIN_LIGHT_MIN[veg_type_], MyMin( RAIN_LIGHT_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRainThrOn( double val ) {
    rain_thr_on_   = MyMax( THR_RAIN_ON_MIN[veg_type_], MyMin( THR_RAIN_ON_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRainThrOff( double val ) {
   rain_thr_off_ = MyMin( MyMax( val, P_50_MIN[veg_type_]), MyMax(rain_thr_on_ - 0.1, P_50_MIN[veg_type_]) ) ;        // MP: constrain rain OFF between the absolute Minimum ON-value and the actual ON-value
   return;
}

void clSeed::setCheckLightThrOn( double val ) {
    light_thr_on_  = MyMax( THR_LIGHT_ON_MIN[veg_type_], MyMin( THR_LIGHT_ON_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckLightThrOff( double val ) {
   light_thr_off_ = MyMin( MyMax( val, THR_LIGHT_ON_MIN[veg_type_] ), MyMax( light_thr_on_ - 0.1, THR_LIGHT_ON_MIN[veg_type_]) ) ; // MP: constrain light OFF between the absolute Minimum ON-value and the actual ON-value
   return;
}

void clSeed::setCheckBiomasPar1( double val ) {
    biomas_par_1_  = MyMax( BIOMAS_PAR_1_MIN[veg_type_], MyMin( BIOMAS_PAR_1_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckBiomasPar2( double val ) {
    biomas_par_2_  = MyMax( BIOMAS_PAR_2_MIN[veg_type_], MyMin( BIOMAS_PAR_2_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRotfrmPar1( double val ) {
    rotfrm_par_1_  = MyMax( ROTFRM_PAR_1_MIN[veg_type_], MyMin( ROTFRM_PAR_1_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRotfrmPar2( double val ) {
    rotfrm_par_2_  = MyMax( ROTFRM_PAR_2_MIN[veg_type_], MyMin( ROTFRM_PAR_2_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckCanfrmPar1( double val ) {
    canfrm_par_1_  = MyMax( CANFRM_PAR_1_MIN[veg_type_], MyMin( CANFRM_PAR_1_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckCanfrmPar2( double val ) {
    canfrm_par_2_  = MyMax( CANFRM_PAR_2_MIN[veg_type_], MyMin( CANFRM_PAR_2_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckCanfrmHBase( double val ) {
    canfrm_h_base_  = MyMax( CANFRM_H_BASE_MIN[veg_type_], MyMin( CANFRM_H_BASE_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckP50( double val ) {
    p_50_  = MyMax( P_50_MIN[veg_type_], MyMin( P_50_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckRotfrmMaxD( double val ) {
    rotfrm_max_d_  = MyMax( ROTFRM_MAX_D_MIN[veg_type_], MyMin( ROTFRM_MAX_D_MAX[veg_type_], val ) ); 
	return;
}

void clSeed::setCheckWoodDensity( double val ) {
    wood_density_  = MyMax( WOOD_DENSITY_MIN[veg_type_], MyMin( WOOD_DENSITY_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckLightExt( double val ) {
    light_ext_     = MyMax( LIGHT_EXT_MIN[veg_type_], MyMin( LIGHT_EXT_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckSeedWeight( double val ) {
    seed_weight_   = MyMax( SEED_WEIGHT_MIN[veg_type_], MyMin( SEED_WEIGHT_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckStorToWood( double val ) {
    stor_to_wood_  = MyMax( STOR_TO_WOOD_MIN[veg_type_], MyMin( STOR_TO_WOOD_MAX[veg_type_], val ) );
	return;
}

void clSeed::setCheckStorToLeaf( double val ) {
    stor_to_leaf_  = MyMax( STOR_TO_LEAF_MIN[veg_type_], MyMin( STOR_TO_LEAF_MAX[veg_type_], val ) );
	return;
}

//----------------------------------------------------------

// called in PlantClass initialize (beginning of run); in PlantClass setInitialValues (runRecruitment); 

void clSeed::setCheckAllocParameters()   // this scales allocation parameters such that they sum to 1
{
	double sum = alloc_root_+alloc_leaf_+alloc_wood_+alloc_bark_+alloc_repr_+alloc_stor_;
	alloc_root_ /= sum;
	alloc_leaf_ /= sum;
	alloc_wood_ /= sum;
	alloc_bark_ /= sum;
	alloc_repr_ /= sum;
	alloc_stor_ /= sum;
	
	sum = stor_to_wood_+stor_to_leaf_;
	if ( sum>1 )
	{
		stor_to_wood_ /= sum;
		stor_to_leaf_ /= sum;
	}
	
	#ifdef W_BRANCHES
		if(alloc_wood_ > 0.)	// dealing with trees
		{
			sum = alloc_branch_ + alloc_stem_ ;  // should equal alloc_wood; in case they don't, scale both to assure they will while keeping their ratio towards one another
			alloc_branch_ = alloc_branch_ * alloc_wood_ / sum ;
			alloc_stem_ = alloc_stem_ * alloc_wood_ / sum ;
	
			stor_to_stem_ = stor_to_wood_ * alloc_stem_/alloc_wood_ ;		// MP, for simpliticty for now let's assume the ratio for reallocation to stem and branches is the same as the regular allocation ratio
			stor_to_branch_ = stor_to_wood_ - stor_to_stem_ ; 
		}	
		else	// dealing with grasses, alloc_wood should never be 0 for trees
		{
			alloc_branch_ = 0. ;
			alloc_stem_ = 0. ;
			stor_to_stem_ = 0. ;
			stor_to_branch_ = 0. ;
		}
	#else
		alloc_stem_	  = alloc_wood_ ;		// check for correctness of alloc_wood_ already done a couple of lines above
		stor_to_stem_ = stor_to_wood_ ;	
	#endif
	
	return;
}

//----------------------------------------------------------

void clSeed::setCheckStemCount( double val ) {		// NEW CAMILLE
    if(w_shrubs_==1) stem_count_  = MyMax( STEM_COUNT_MIN[veg_type_], MyMin( STEM_COUNT_MAX[veg_type_], val ) );
    else stem_count_ = 1.;
	return;
}






















