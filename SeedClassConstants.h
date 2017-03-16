#ifndef SeedClassConstants_h___
#define SeedClassConstants_h___

using namespace::std;


int const    VTYPE_NUM                      = 2;    // number of vegetation types

int const    TREE                           = 0;
int const    GR_4                           = 1;
int const    GR_3                           = 2;

int const    SPECIES_NUM                    = 50;
// int const    SPECIES_NUM                    = 2;

// double mrd = IData.root_depth_;
// -------------------------------------------------------------------------------
// --- TRAITS THAT ARE MODIFIED BY THE GENETIC ALGORITHM -------------------------
// -------------------------------------------------------------------------------

// double const SLA_MIN[VTYPE_NUM]             = {        2.,        2. };	// BEFORE LIAM CHANGES
// double const SLA_MAX[VTYPE_NUM]             = {       30.,       25. };	// BEFORE LIAM CHANGES
double const SLA_MIN[VTYPE_NUM]             = {        0.,        0. };	// NEW LIAM
double const SLA_MAX[VTYPE_NUM]             = {        0.,        0. };	// NEW LIAM

double const ALLOC_ROOT_MIN[VTYPE_NUM]      = {      0.2,      0.2 };		// carbon allocation to roots
double const ALLOC_ROOT_MAX[VTYPE_NUM]      = {       0.4,       0.8 };

double const ALLOC_LEAF_MIN[VTYPE_NUM]      = {      0.35,       0.25 };		// carbon allocation to leaves
double const ALLOC_LEAF_MAX[VTYPE_NUM]      = {       0.5,       0.5 };		// NEW LIAM

double const ALLOC_WOOD_MIN[VTYPE_NUM]      = {      0.25,       0.0 };		// carbon allocation to stem and branches combined
double const ALLOC_WOOD_MAX[VTYPE_NUM]      = {       0.35,       0.0 };

#ifdef W_BRANCHES
	double const ALLOC_BRANCH_MIN[VTYPE_NUM] 	= {		  /*0.05*/0.15, 		 0.0 };		// fraction of alloc_wood that goes to branches; autecological, genetically determined fraction; MP
	double const ALLOC_BRANCH_MAX[VTYPE_NUM] 	= {		  /*0.1*/0.25, 		 0.0 };	
#endif	

double const ALLOC_BARK_MIN[VTYPE_NUM]      = {       0.001,       0.0 };		// carbon allocation to bark/defence
double const ALLOC_BARK_MAX[VTYPE_NUM]      = {       0.05,       0.0 };		

double const ALLOC_REPR_MIN[VTYPE_NUM]      = {       0.05,     0.05 };		// carbon allocation to reproduction
double const ALLOC_REPR_MAX[VTYPE_NUM]      = {       0.2,       0.2 };


double const ALLOC_STOR_MIN[VTYPE_NUM]      = {       0.1,       0.0 };		// carbon allocation to storage
double const ALLOC_STOR_MAX[VTYPE_NUM]      = {       0.4,       0.4 };

double const RAIN_LIGHT_MIN[VTYPE_NUM]      = {        0,         0  };		// pheneology control by rain or light
double const RAIN_LIGHT_MAX[VTYPE_NUM]      = {        1.,        1. };		// 0=rain, 1=light

double const EVERGREEN_MIN[VTYPE_NUM]       = {        0.,        1. };		// evergreen
double const EVERGREEN_MAX[VTYPE_NUM]       = {        1.,        1. };		// 0=no/deciduous, 1=yes/evergreen

double const THR_RAIN_ON_MIN[VTYPE_NUM]     = {       -3.0,       -3.0 };	// NEW LIAM
double const THR_RAIN_ON_MAX[VTYPE_NUM]     = {        -0.2,        -0.2};	// NEW LIAM

//double const THR_RAIN_OFF_MIN[VTYPE_NUM]    = {       -12.0,       -12.0 };	// NEW LIAM	
//double const THR_RAIN_OFF_MAX[VTYPE_NUM]    = {        -0.5,        -0.5};	// NEW LIAM

double const THR_RAIN_OFF_MIN[VTYPE_NUM]    = {       0.1,       0.1 };	// NEW MIRJAM; to subtract from ON-threshold to obtain OFF-threshold; eyeballed this, can be fine-tuned if necessary
double const THR_RAIN_OFF_MAX[VTYPE_NUM]    = {       1.5,       1.5 };		// NEW MIRJAM

double const THR_LIGHT_ON_MIN[VTYPE_NUM]    = {        6.,        6. };
double const THR_LIGHT_ON_MAX[VTYPE_NUM]    = {       14.,       14. };

double const THR_LIGHT_OFF_MIN[VTYPE_NUM]   = {        0.1,       0.1 };		// NEW MIRJAM; to subtract from ON-threshold to obtain OFF-threshold; eyeballed this, can be fine-tuned if necessary
double const THR_LIGHT_OFF_MAX[VTYPE_NUM]   = {        2.,        2. };

double const BIOMAS_PAR_1_MIN[VTYPE_NUM]    = {       0.1,       2.3 };
double const BIOMAS_PAR_1_MAX[VTYPE_NUM]    = {       0.8,       2.7 };

double const BIOMAS_PAR_2_MIN[VTYPE_NUM]    = {       2.3,       1.8 };
double const BIOMAS_PAR_2_MAX[VTYPE_NUM]    = {       3.6,       2.2 };

double const ROTFRM_PAR_1_MIN[VTYPE_NUM]    = {      0.01,      0.01 };
double const ROTFRM_PAR_1_MAX[VTYPE_NUM]    = {       10.,       10. };

// double const ROTFRM_PAR_2_MIN[VTYPE_NUM]    = {      -20.,      -20. };	// BEFORE LIAM CHANGES
double const ROTFRM_PAR_2_MIN[VTYPE_NUM]    = {       -1.,        1. };		// NEW LIAM
double const ROTFRM_PAR_2_MAX[VTYPE_NUM]    = {       20.,       20. };

double const ROTFRM_MAX_D_MIN[VTYPE_NUM]    = {       0.3,       0.3 };		// meters
double const ROTFRM_MAX_D_MAX[VTYPE_NUM]    = {       3.6,       2.4 };		// meters
// double const ROTFRM_MAX_D_MIN[VTYPE_NUM]    = {       0.3,       0.3 };
// double const ROTFRM_MAX_D_MAX[VTYPE_NUM]    = {       0.3,       0.3 };

// double const ROTFRM_MAX_D_MIN[VTYPE_NUM]    = {       1.0,       1.0 };
// double const ROTFRM_MAX_D_MAX[VTYPE_NUM]    = {       13.0,       2.0 };

double const LIGHT_EXT_MIN[VTYPE_NUM]       = {       0.4,       0.4 };		// crown light extinction parameter
double const LIGHT_EXT_MAX[VTYPE_NUM]       = {       0.5,       0.5 };		// crown light extinction parameter

// double const WOOD_DENSITY_MIN[VTYPE_NUM]    = {      500.,      200. };	// wood density, kg/m^3
// double const WOOD_DENSITY_MAX[VTYPE_NUM]    = {      900.,      900. };	// wood density, kg/m^3
// double const WOOD_DENSITY_MIN[VTYPE_NUM]    = {      200.,      200. };	// wood density, kg/m^3
// double const WOOD_DENSITY_MAX[VTYPE_NUM]    = {     1000.,      200. };	// wood density, kg/m^3
double const WOOD_DENSITY_MIN[VTYPE_NUM]    = {      100.,      100. };		// wood density, kg/m^3
double const WOOD_DENSITY_MAX[VTYPE_NUM]    = {     1200.,      1200. };	// wood density, kg/m^3

double const SEED_WEIGHT_MIN[VTYPE_NUM]     = {     0.001,   0.0001 };		// seed weight, kg
double const SEED_WEIGHT_MAX[VTYPE_NUM]     = {      0.05,     0.05 };		// seed weight, kg

// double const CANFRM_PAR_1_MIN[VTYPE_NUM]    = {      0.01,      0.01 };	// maximum crown radius
// double const CANFRM_PAR_1_MAX[VTYPE_NUM]    = {       0.4,       0.4 };	// (as fraction of height)
// double const CANFRM_PAR_1_MIN[VTYPE_NUM]    = {       0.2,      0.15 };	// maximum crown radius
// double const CANFRM_PAR_1_MAX[VTYPE_NUM]    = {       0.4,      0.25 };	// (as fraction of height)
// double const CANFRM_PAR_1_MIN[VTYPE_NUM]    = {       0.3,      0.3 };	// maximum crown radius
// double const CANFRM_PAR_1_MAX[VTYPE_NUM]    = {       0.3,      0.3 };     	// (as fraction of height)
double const CANFRM_PAR_1_MIN[VTYPE_NUM]    = {       0.2,      1.0 };     // maximum crown radius
double const CANFRM_PAR_1_MAX[VTYPE_NUM]    = {       0.4,      1.0};     // (as fraction of height)

double const CANFRM_PAR_2_MIN[VTYPE_NUM]    = {        1.,       20. };     // slope of crown: 1=linear
double const CANFRM_PAR_2_MAX[VTYPE_NUM]    = {       25.,       60. };     // >1=more cylindrical, <1=more copped

double const CANFRM_H_BASE_MIN[VTYPE_NUM]   = {       0.01,       0.01 };     // Cam - crown lower limit (0 = up to the ground)
double const CANFRM_H_BASE_MAX[VTYPE_NUM]   = {        .5,        .95 };     // Cam - crown upper limit (1 = up to max height) (grass assumed to always have crown up to base)


// double const CANFRM_PAR_2_MIN[VTYPE_NUM]    = {        1.,       1. };     // slope of crown: 1=linear
// double const CANFRM_PAR_2_MAX[VTYPE_NUM]    = {       1.,       1. };     // >1=more cylindrical, <1=more copped
// double const CANFRM_PAR_2_MIN[VTYPE_NUM]    = {       0.1,       0.1 };     // slope of crown: 1=linear
// double const CANFRM_PAR_2_MAX[VTYPE_NUM]    = {       20.,       30. };     // >1=more cylindrical, <1=more copped

double const P_50_MIN[VTYPE_NUM]            = {      -3.0,     -3.0 };		//  p50
double const P_50_MAX[VTYPE_NUM]            = {       -0.2,      -0.2 };

const double STOR_TO_WOOD_MIN[VTYPE_NUM]    = {        0.,        0. };		// allocation of storage to wood after fire
const double STOR_TO_WOOD_MAX[VTYPE_NUM]    = {       0.4,        0. };		// allocation of storage to wood after fire

const double STOR_TO_LEAF_MIN[VTYPE_NUM]    = {        0.6,        0. };		// allocation of storage to leaves after fire
const double STOR_TO_LEAF_MAX[VTYPE_NUM]    = {       0.9,       0.7 };		// allocation of storage to leaves after fire

double const STEM_COUNT_MIN[VTYPE_NUM]      = {      1.,     1. };		//  number of stems - NEW CAMILLE
double const STEM_COUNT_MAX[VTYPE_NUM]      = {      10.,    1. };		//  number of stems - NEW CAMILLE  // SIMON 1 for grasses


// -------------------------------------------------------------------------------
// --- TRAITS THAT ARE RANDOMLY SELECTED BUT NOT OPTIMIZED -----------------------
// -------------------------------------------------------------------------------


// Value for                                        TREES   C4 GRASS



// -------------------------------------------------------------------------------
// --- TRAITS THAT ARE CONSTANT --------------------------------------------------
// -------------------------------------------------------------------------------

// Value for                                        TREES   C4 GRASS
double const MORT_CARBON[VTYPE_NUM]		= {      0.3,        0.3 };     // mortality due to negative carbon balance
// double const MORT_BACKGROUND[VTYPE_NUM]		= {     0.00,      0.05 };     	// background mortality
double const MORT_BACKGROUND[VTYPE_NUM]		= {     0.00,      0.01 };     	// background mortality

double const MORT_BIOMASS[VTYPE_NUM]		= {      0.05,      0.05 };     // mortality due to low height
double const MORT_MECHANIC_1[VTYPE_NUM]		= {       10.,        5. };     // mortality due to mechanic instability
double const MORT_MECHANIC_2[VTYPE_NUM]		= {        6.,        6. };     // mortality due to mechanic instability
//double const MORT_MECHANIC_1[VTYPE_NUM]		= {       75.,       10. };     // mortality due to mechanic instability
//double const MORT_MECHANIC_2[VTYPE_NUM]		= {        6.,        6. };     // mortality due to mechanic instability

const double TOPKILL_CONST[VTYPE_NUM]		= {       4.3,        0. };     // topkill constants (Higgins et al. 2000)
const double TOPKILL_H[VTYPE_NUM]		= {     5.003,        0. };     // topkill constants (Higgins et al. 2000)
const double TOPKILL_I[VTYPE_NUM]		= {  0.004408,        0. };     // topkill constants (Higgins et al. 2000)

const double BALL_BERRY_M[VTYPE_NUM]		= {       9.0,       4.0 };     // ball berry equation parameter
const double BALL_BERRY_B[VTYPE_NUM]		= {      0.01,      0.04 };     // ball berry equation parameter

const double R_MAINT_RESP[VTYPE_NUM]		= {     0.015,     0.025 };     // maintenance respiration parameter, stomatal conductance
const double R_GROWTH_RESP[VTYPE_NUM]		= {      0.35,      0.35 };     // growth respiration parameter

const double BETA_LEAF_RESP[VTYPE_NUM]		= {      0.01,      0.01 };		// fraction of leaf biomass that respires
// const double BETA_LEAF_RESP[VTYPE_NUM]		= {      0.04,      0.02 };		// fraction of leaf biomass that respires--from some reference I found
const double BETA_WOOD_RESP[VTYPE_NUM]		= {      0.01,      0.01 };		// fraction of wood biomass that respires
const double BETA_ROOT_RESP[VTYPE_NUM]		= {      0.01,      0.01 };		// fraction of root biomass that respires
// const double BETA_ROOT_RESP[VTYPE_NUM]		= {      0.04,      0.02 };		// fraction of root biomass that respires

const double BETA_N_RESP[VTYPE_NUM]			= {     0.281,     0.281 };		// parameter for respiration model

const double CN_RATIO_LEAF[VTYPE_NUM]		= {      120.,      120. };		// C:N ratio of leaves
const double CN_RATIO_WOOD[VTYPE_NUM]		= {      150.,      120. };		// C:N ratio of woody biomass
const double CN_RATIO_ROOT[VTYPE_NUM]		= {       60.,      120. };		// C:N ratio of roots


#endif











