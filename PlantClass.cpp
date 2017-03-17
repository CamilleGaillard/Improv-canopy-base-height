#include "PlantClass.h"
#include "LeafGlobals.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "GridCellClass.h"
using namespace std;

//-----------------------------------------------------------------------------------

/// Standard constructor
/// Initialization with random seed
/// Only called for initial plant population

clPlant::clPlant()
{
}

void clPlant::initialize( clGridCell * gridcell_, int id_in)
{	
	id_ = id_in;		// individual number index; to keep track of specific individuals
	sd_.initialize(gridcell_);
	active_ = 1.;		// everybody active at beginning; turn-off-switches should switch them off if they shouldn't be active; otherwise, we ought to initialize leaf biomass to zero if we initialize active to 0!
	c_gain_ = 0.;
	c_gain_cum_ = 0.;
	c_gain_mmolm2s_ =0.;
	c_gain_cum_flag_ = 0;
	
	double tmp_bm   = RUnif( 0.05, .1 );  // changed minimum value from 0.001 to 0.05, to make sure the initial height does not go below 10 cm (this only works for single-stem at the moment)
	if ( sd_.getVegType()== TREE )
	{
		b_leaf_			= tmp_bm*0.07;		// default 0.1
		b_wood_			= tmp_bm*0.43;
		#ifdef W_BRANCHES
			b_stem_			= tmp_bm*0.42;	// MP
			b_branch_		= b_wood_-b_stem_;			
		#else
			b_stem_			= b_wood_ ;
		#endif		
		b_root_			= tmp_bm*0.3;
		b_bark_			= tmp_bm*0.05;
		b_repr_			= tmp_bm*0.05;
		b_stor_			= tmp_bm*0.1;
	}
	else
	{
		b_leaf_			= tmp_bm*0.2;
		b_wood_			= 0.;
		#ifdef W_BRANCHES
			b_stem_			= b_wood_;
			b_branch_		= 0.;
		#else
			b_stem_			= b_wood_;
		#endif		
		b_root_			= tmp_bm*0.6;
		b_bark_			= 0.;
		b_repr_			= tmp_bm*0.1;
		b_stor_			= tmp_bm*0.1;
	}
	
//	if(sd_.getEvergreen() < EVERGREEN_THR) cout << id_ << setw(14) << tmp_bm << setw(14) << b_leaf_ << setw(14) << b_wood_ << setw(14) << b_stem_ << setw(14) << b_bark_ << setw(14) << b_root_ << setw(14) << b_stor_ << setw(14) << b_repr_ << endl;
	
	b_leaf_d_st_ = 0.;
	b_leaf_d_ly_ = 0.;
	
	age_			= MyRound( RUnif()*15 );
	dpos_			= 0;
	dneg_			= 0;
// 	c_def_counter_	= 0;
	dleaf_alloc_	= 0;
	evg_flush_flag_ = 0;
	day_ = 0;
// 	day1_ = 0;
// 	Rm_Sum_ = 0;
	
	soilcarbon_nwl_ = 0.;  // SIMON non woody litter for soil C model
	soilcarbon_fwl_ = 0.;  // SIMON fine woody litter for soil C model
	soilcarbon_cwl_ = 0.;  // SIMON coarses woody litter for soil C model
	
	if ( sd_.getEvergreen() > EVERGREEN_THR ) active_ = 1.;

	calWoodDensity() ;						// needed to calculate stem diameter
	#ifdef W_BRANCHES
		calStemCount() ;		// needed in height calucation
		height_ = exp(sd_.getBiomasPar1()/sd_.getBiomasPar2()) * Mypow(b_stem_/stem_count_, 1./sd_.getBiomasPar2()) ; // starting height
		calCrownBHeight();
		calCrownVolume();
	#else		
		calStemCount(); 		// NEW CAMILLE								// MIRJAM: moved this up above calHeight(), otherwise stem_count_ won't be initialized in calHeight() !!!
		calHeight();  			// update tree height
		calCrownBHeight();
		calStemDiam(); 			// update stem diameter
	#endif	
	calCrownArea();
	calLai();
	stom_cond_ = 0.;
	stom_cond_max_ = 0.;
	boul_cond_ = 0.;
	evapotr_   = 0.;
	evapotr_max_   = 0.;
	leaf_temp_ = 0.;
// 	calWiltPoint( thetawp, thetafc );
	alive_ = 1;
	
//	cout << crown_l_*100. << setw(14) << crown_volume_/crown_l_ << setw(14) << sd_.getCanfrmPar2() << setw(14) << 2./sd_.getCanfrmPar2()+1. << setw(14) << crown_r_*100. << setw(14) << par_d_ << setw(14) << crown_area_0_ << setw(14) << crown_r_/crown_l_ 
//	<< setw(14) << b_stem_ << setw(14) << (b_branch_+b_leaf_)/b_stem_ << setw(14) << rho_crown_ << setw(14) << (b_leaf_+b_branch_)*(2./sd_.getCanfrmPar2()+1.)/(M_PI*0.16*pow(crown_l_,3.)) << endl ;

//	calAreaSapwood(); // NEW CAMILLE	

	//light_ext_ = 0.3+0.3*sd_.getCanfrmHBase();	// SIMON
	light_ext_ = 0.4+0.3*sd_.getCanfrmHBase();	// SIMON

	calLightAvail();
	
	calLDMC();	// leaf dry matter content fn. of P50				// NEW LIAM
	calSLA();    // 0.1 instead of 0.01 to transform from g/cm^2 to kg/m^2				// NEW LIAM

	calLeafLongevity();  // factor 10 to transform units // using p50
	//root_turnover_ = 0.00137;  // NEW LIAM; ~70 % grass root turnover per year
	root_turnover_ = 0.00245; // ca. 90% turnover  SIMON
	root_turnover_tree_ = 0.0005; // ca. 20% turnover

	can_frm_helper_ = pow( 0.5, 1./sd_.getCanfrmPar2() );
	
	double sum = 0;
	double C1 = sd_.getRotfrmPar1();
	double C2 = sd_.getRotfrmPar2();
	double MD = sd_.getRotfrmMaxD();
	
	for ( int i=0; i<N_SOIL_LAYERS; i++ )
	{
		if ( (DEPTH[i]-THICK[i]*0.5) < MD )
			root_fracs_[i] =      pow( ( 1. - pow( (DEPTH[i]-THICK[i]*0.5)/MD, C1 )  ), C2 )    *THICK[i];
		else
			root_fracs_[i] = 0;
		
		sum += root_fracs_[i];
		root_bm_[i]=0;
	}
	if (sum>0)
	{
		for ( int i=0; i<N_SOIL_LAYERS; i++ )
		{
			root_fracs_[i] /= sum;
		}
	}

	phen_loop_=0;
	for ( int i=0; i<N_PHEN; i++ ) phen_rain_[i] = 0.;
	for ( int i=0; i<N_PHEN; i++ ) phen_lght_[i] = 0.;
		
	// minimum avaids that respiration exceeds biomass, resp=helper*biomass*temperature_factor where temperature_factor <2.05

	leaf_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaLeafResp()/sd_.getCNRatioLeaf() );  // *MAIN_RESP_CARB_PROPORT(=1)
	wood_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaWoodResp()/sd_.getCNRatioWood() );  // *MAIN_RESP_CARB_PROPORT(=1)
	root_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaRootResp()/sd_.getCNRatioRoot() );  // *MAIN_RESP_CARB_PROPORT(=1)
	//cout<<"leafresp: "<<leaf_resp_helper_<<"woodresp: "<<wood_resp_helper_<<"root resp: "<<root_resp_helper_<<endl;
	sd_.setCheckAllocParameters();
	
//	if(id_==1287) cout << "PLANT CLASS INIT: " <<setw(14)<< sd_.getVegType() <<setw(14)<< sd_.getEvergreen() <<setw(14)<< active_ <<setw(14)<< b_leaf_ <<setw(14)<< b_root_ <<setw(14)<< b_repr_ <<setw(14)<< b_stor_ <<setw(14)<< age_ << endl;
}

//-----------------------------------------------------------------------------------


/// Fill a vector with trait data for output

void clPlant::getRootBmVector( double *vec )
{
	for (int i=0; i<N_SOIL_LAYERS; i++ )
	  {
	    vec[i]=root_bm_[i];  
	  }
    return;
}

//-----------------------------------------------------------------------------------


void clPlant::calLDMC()
{
	ldmc_ = 0.224+(-0.41*sd_.getP50());	// leaf dry matter content fn. of P50
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calSLA()
{
//	double leaf_thickness_grass_ = 0.3;
	if (sd_.getVegType()== GR_4) sla_  = 1./(ldmc_*leaf_thickness_grass_);   
	if (sd_.getVegType()== TREE)  sla_  = 1./(ldmc_*leaf_thickness_tree_); 
	//if (RUnif()<0.05) cout << endl << setw(40) << sd_.getVegType() << setw(14) << sla_ << setw(14) << ldmc_ << setw(14)<< sd_.getP50() ;
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calWoodDensity()
{
	// 1.54 to account for moisture content of measured dry wood density at 15% moisture content to 50% moisture content
	wood_density_ = 1.54*(0.259+(-0.05921*sd_.getP50()))*1000.; // 1000 converts density to old units
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calLeafLongevity()
{
	//leaf_longevity_ = 1./(29664.*pow( 10.*sd_.getSla(), -0.909 ))*0.1;  // factor 10 to transform units		// BEFORE LIAM CHANGES; MP: the factor *0.1 is definitively not correct!
	leaf_longevity_ = 1./(29664.*pow( 10.*sla_, -0.909 ));  // factor 10 to transform units // using p50
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calStemCount()  // NEW CAMILLE
{

	//if(sd_.getVegType()==TREE)
	stem_count_ = sd_.getStemCount();

	return;

}

//-----------------------------------------------------------------------------------


void clPlant::calHeight()
{

	length_leaf_total_grass_ = 0.0;		// NEW LIAM

	
	if(sd_.getVegType()==GR_4)			// NEW LIAM
	{
		//double leaf_thickness_grass_ = 0.1;	//mm		// NEW LIAM
		//length_leaf_total_grass_ = b_leaf_/((leaf_thickness_grass_*leaf_width_grass)/1000);		// NEW LIAM
		//height_ = MyMin(1.5, (length_leaf_total_grass_)); //m 

		double leaf_width_grass = (1.25*sd_.getP50())+22.5; //mm		// NEW LIAM
		double LTD = ((-0.025*sd_.getP50()) + 0.294) *0.001;		// MP: REF??? => this delivers values ranging between 0,000469 and 0,000319 for leaf tissue density leaf tissue density g/cm^3 
		double VOL_leaf = 0.7*1000.*(b_leaf_/LTD);// mm^3
	
		length_leaf_total_grass_ = (VOL_leaf / (leaf_thickness_grass_ * leaf_width_grass)) / 1000;	// MP: works if VOL_leaf is indeed in mm^3, 0.1 equals leaf thickness in mm, and leaf_width is in mm; then lenght_leaf is in [m]
		height_ = MyMin(1.5, (length_leaf_total_grass_)); //m 
		
//		cout << "LTD, p50, VOL_leaf, VOL_leaf^(1/3), b_leaf, height, No. : " << LTD << ", " << sd_.getP50() << ", " << VOL_leaf/1000. << ", " << Mypow((VOL_leaf/1000.),(1./3.)) << ", " << b_leaf_ << ", " << height_ << ", " << length_leaf_total_grass_/1.5 << endl;		
	}
	else  // we're dealing with trees
	{
	
		#ifdef W_BRANCHES
		{	
			//First step: take the newly updated stem biomass to calculate the potential height and crown lenght that the tree could have based only on stem biomass. Keep the old values, in case we need them again. 

			double crown_l_old = crown_l_ ;			// need to initialize crown_l_ in initialize() and setInitialValues(), where height_ and crown_h_base_ is initialized
			double crown_r_old = crown_r_ ;			// initialize crown_r_ in initialize() and setInitialValues()
			double crown_v_old = crown_volume_ ;	// initialize crown_volume_ in initialize() and setInitialValues()
			double height_old = height_;

			height_ = MyMax(exp( sd_.getBiomasPar1()/sd_.getBiomasPar2()) * pow((b_stem_)/stem_count_, 1./sd_.getBiomasPar2()), 0.001 );		// this is a mathematically equivalent simplification of the preceeding equation  // avoids that height=0 which leads to problems in stem diameter calculations;

			if(height_ >= height_old)	// stem biomass increment must have been positive
			{
				crown_l_ += (height_ - height_old) ;		// try to increase crown length by the increment in height 
				calCrownVolume() ;
				
//				if(crown_volume_>= crown_v_old && crown_l_ >= crown_l_old && crown_r_ >= crown_r_old) cout << "Case 1 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ >= crown_l_old && crown_r_ < crown_r_old)  cout << "Case 2 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ < crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 3 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ < crown_l_old && crown_r_ < crown_r_old)	  cout << "Case 4 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(year_==0 && id_==6 && crown_volume_< crown_v_old && crown_l_ >= crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 5 "<< dct_ <<setw(14)<< height_ <<setw(14)<< crown_h_base_ <<setw(14)<< b_leaf_ <<setw(14)<< b_branch_ <<setw(14)<< b_leaf_+b_branch_ <<setw(14)<< crown_volume_ <<setw(14)<< crown_v_old <<setw(14)<< crown_l_ <<setw(14)<< crown_l_old <<setw(14)<< crown_r_ <<setw(14)<< crown_r_old << endl;
//				if(year_==0 && id_==6 && crown_volume_< crown_v_old && crown_l_ >= crown_l_old && crown_r_ < crown_r_old)  cout << "Case 6 "<< dct_  <<setw(14)<< height_ <<setw(14)<< crown_h_base_ <<setw(14)<< b_leaf_ <<setw(14)<< b_branch_ <<setw(14)<< b_leaf_+b_branch_ <<setw(14)<< crown_volume_ <<setw(14)<< crown_v_old <<setw(14)<< crown_l_ <<setw(14)<< crown_l_old <<setw(14)<< crown_r_ <<setw(14)<< crown_r_old << endl;
//				if(crown_volume_< crown_v_old && crown_l_ < crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 7 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_< crown_v_old && crown_l_ < crown_l_old && crown_r_ < crown_r_old)	cout << "Case 8 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
				
				if(crown_r_ < crown_r_old)	// crown bimass increment not sufficient to keep up with height increment; 
				{
					crown_l_ = crown_volume_ * (2./sd_.getCanfrmPar2() + 1.) / (M_PI * crown_r_old * crown_r_old) ;  // keep radius constant to old value, shorten crown
					crown_h_base_ = height_ - crown_l_ ;	// crown base height migrates upwards 					 // push crown base upward
				}	
					
			}
			else	// stem biomass increment must have been negative (e.g., because resp. > net C-gain). 
			{
				crown_l_ *= height_ / height_old ;		// try to shorten crown length to avoid that it can get longer than the height; can this get zero or negative?
				
				calCrownVolume() ;		// if stem biomass decreased, crown biomass ought to have decreased as well (c_balance_ negative)
				
//				if(crown_volume_>= crown_v_old && crown_l_ >= crown_l_old && crown_r_ >= crown_r_old) cout << "Case 9 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ >= crown_l_old && crown_r_ < crown_r_old)  cout << "Case 10 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ < crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 11 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_>= crown_v_old && crown_l_ < crown_l_old && crown_r_ < crown_r_old)	  cout << "Case 12 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_< crown_v_old && crown_l_ >= crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 13 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_< crown_v_old && crown_l_ >= crown_l_old && crown_r_ < crown_r_old)  cout << "Case 14 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_< crown_v_old && crown_l_ < crown_l_old && crown_r_ >= crown_r_old)  cout << "Case 15 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
//				if(crown_volume_< crown_v_old && crown_l_ < crown_l_old && crown_r_ < crown_r_old)	cout << "Case 16 "<< id_ <<setw(14)<< year_ <<setw(14)<< dct_ << endl;
				
				if(crown_l_ <= crown_l_old && crown_r_ > crown_r_old)	// crown biomass didn't decrease as much as crown length
				{
					crown_l_ = crown_l_ = crown_volume_ * (2./sd_.getCanfrmPar2() + 1.) / (M_PI * crown_r_old * crown_r_old) ;	// do not let the crown radius get larger albeit crown volume decreased!
					crown_h_base_ = height_ - crown_l_ ;	// crown base height will go down
					if(crown_h_base_ <= 0.)		// bad-ass case, rescale completely to get out of that situation (although not very realistic, but then trees shrinking in such a way isn't either...)
					{
						calCrownBHeight() ;  	// re-intializes crown base height and crown length
						calCrownVolume() ;		// re-initalize crown_r_ and par_d_ based on new crown length delivered from calCrownBHeight						
					}
				}
				else if(crown_volume_ > crown_v_old)	// should be unlikely; maybe if phenology just got switched on?; radius will increase, but that should be okay in this case
				{
					crown_h_base_ = height_ - crown_l_ ;	// crown base height will go down
					if(crown_h_base_ <= 0.)		// bad-ass case, rescale completely to get out of that situation (although not very realistic, but then trees shrinking in such a way isn't either...)
					{	
						calCrownBHeight() ;  	// re-intializes crown base height and crown length
						calCrownVolume() ;		// re-initalize crown_r_ and par_d_ based on new crown length delivered from calCrownBHeight						
					}					
				}
			}	
		}
		#else
		
//			height_ = exp( (log(b_stem_/stem_count_)+sd_.getBiomasPar1())/sd_.getBiomasPar2() ); // NEW CAMILLE 
			double height_old = height_ ;		// FLAG: for debugging-purposes 
			height_ = MyMax(exp( sd_.getBiomasPar1()/sd_.getBiomasPar2()) * pow((b_stem_)/stem_count_, 1./sd_.getBiomasPar2()), 0.001 );		// this is a mathematically equivalent simplification of the preceeding equation  // avoids that height=0 which leads to problems in stem diameter calculations
			calStemDiam() ;
//			if (height_<height_old) cout << "PROBLEM, height just decreased! " << year_ << setw(14) << dct_ << setw(14) << id_ << setw(14) << active_ << setw(14) << round(sd_.getEvergreen()) << setw(14) << b_stem_ << setw(14) << b_leaf_ 
//													<< setw(14) << height_old << setw(14) << height_ << setw(14) << c_balance_ << endl;

			//double heightoldequiv = exp( (log((double)b_stem_+(double)b_leaf_)+(double)sd_.getBiomasPar1())/(double)sd_.getBiomasPar2() );
			//if (abs((height_/heightoldequiv)-1.)>0.1) cout << "Height/oldcomp: "<<heightoldequiv<<" "<<height_<<" "<<exp((double) sd_.getBiomasPar1()/(double)sd_.getBiomasPar2())<<" "<<pow((double)b_stem_/(double)stem_count_, 1./(double)sd_.getBiomasPar2())<<endl;
		#endif	 			 				
	}		
#ifdef NANCHECK
		if(std::isnan(height_)) cout << "soil moisture NaN at calHeight" <<endl;
		if(std::isnan(crown_r_)) cout << "crown radius NaN at calHeight" <<endl;
#endif


	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calCrownBHeight()
{
	crown_h_base_ = height_*sd_.getCanfrmHBase();   // SIMON, CAMILLE
	#ifdef W_BRANCHES
		crown_l_ = height_ - crown_h_base_ ;
	#endif	
	//cout << sd_.getVegType() << setw(14) << sd_.getCanfrmHBase() << setw(14) << crown_h_base_ << setw(14) << height_ << endl;
	return;
}


//-----------------------------------------------------------------------------------

void clPlant::calStemDiam()
{
	
	double stem_diam_int_one_ ; // NEW CAMILLE
	double stem_diam_tot_one_ ; // NEW CAMILLE
	
	// stem diameter is now calculated by using plant height and wood density
	// it is assumed that the stem is a cylinder
	// here, the bark biomass is included for stem diameter calculations,
	// we assume that bark density is 1/2 of wood density  NOTE need to find parameter!
  //----------------------------------------------	
  // P50 method
//	if(sd_.getVegType()==GR_4) wood_density_ = 200.;	// FLAG: MP: This is how it used to be for grasses in the old version; breaks connection between rooting depth and SLA
//	stem_diam_int_ = 2.*Mypow( (b_stem_            ) / ( M_PI*height_*wood_density_ ), 0.5 );						// NEW LIAM
//	stem_diam_tot_ = 2.*Mypow( (b_stem_+2.*b_bark_ ) / ( M_PI*height_*wood_density_ ), 0.5 );						// NEW LIAM
	stem_diam_int_ = 2.*Mypow( (b_stem_/stem_count_            ) / ( M_PI*height_*wood_density_ ), 0.5 );	// NEW LIAM + NEW CAMILLE
	stem_diam_int_one_ = 2.*Mypow( (b_stem_) / ( M_PI*height_*wood_density_ ), 0.5 );	// NEW CAMILLE
	
	stem_diam_tot_ = 2.*Mypow( ((b_stem_/stem_count_)+2.*(b_bark_/stem_count_) ) / ( M_PI*height_*wood_density_ ), 0.5 );		// NEW LIAM + NEW CAMILLE
	stem_diam_tot_one_ = 2.*Mypow( ((b_stem_)+2.*(b_bark_) ) / ( M_PI*height_*wood_density_ ), 0.5 );  // NEW CAMILLE
	
	bark_thick_    = (stem_diam_tot_ - stem_diam_int_)*0.5;
 	//cout << height_<<endl;
	
//	stem_diam_int_ = 2.*Mypow( (b_stem_            ) / ( M_PI*height_*sd_.getWoodDensity() ), 0.5 );		// BEFORE LIAM CHANGE
//	stem_diam_tot_ = 2.*Mypow( (b_stem_+2.*b_bark_ ) / ( M_PI*height_*sd_.getWoodDensity() ), 0.5 );		// BEFORE LIAM CHANGE
//	bark_thick_    = (stem_diam_tot_ - stem_diam_int_)*0.5;													// BEFORE LIAM CHANGE
//	#ifdef W_BRANCHES
//		calCrownVolume(); //FLAG! refresh crown Volume calculation based on new stem thickness and therefore new crown density. This currently might reduce the crown diameter in case the crown density has increased.
//	#endif	
	return;
}

//-----------------------------------------------------------------------------------

//void clPlant::calAreaSapwood()  // NEW CAMILLE 
//{
	////if(sd_.getVegType()==TREE)
	//area_sapwood_ = 1.582*pow(stem_diam_tot_ , 1.764) ;	//ref = Meinzer et al. - Regulation of water flux through tropical forest canopy trees- Do universal rules apply - Tree Physiol-2001 - fig 1
	
	//return;

//}


//-----------------------------------------------------------------------------------

// calculates crown volume based on crown biomass and crown bulk density

#ifdef W_BRANCHES

// calculates crown volume based on crown biomass and crown bulk density, and crown radius based on volume, crown shape, and crown length

void clPlant::calCrownVolume()
{

			// crown volume, based on crown biomass and crown bulk density; crown bulk density may turn out to be a parameter that needs further fine tuning, since values in the literature for crown bulk density
			// are mostly related to crown fire modelling, i.e., only taking into account leaves and fine twigs, whereas I only have leaves and total branch biomass available at this point to deduce crown volume;
			// most helpful information I could find is from "Utilization of Residual Forest Biomass" (Pentti Hakkila), with values for Picea abies (pp. 55). Nonlinear least squares regression of crown bulk density against 
			// stem diameter based on values in table on p. 55). Parameters may need further adjustment, or could become trait values. 
			
			double par_a = 2.6	;			// original: 2.5486465 ;
			double par_b = -3.	;			// original: -6.6365189 ;
			double par_c = 0.5	;			// original: 0.6865899 ;

			calStemDiam() ;							// needed to calculate rho_crown 

			rho_crown_ = MyMax(par_a + par_b / (pow( stem_diam_int_*100., par_c)), 1.0) ;  // stem diameter needs to be in cm for this regression
//			rho_crown_ = 0.9 ;
			crown_volume_ = (b_leaf_ + b_branch_) / rho_crown_ / stem_count_ ;
			
			// convert crown volume to crown radius, including crown shape parameter CanfrmPar2; CanfrmPar2 = 1 equals conical crown shape, CanfrmPar2 = 25 equals approx. cylindrical crown shape;
			// use the function describing radius at a given distance l downwards from the tip as a function of l and CanfrmPar2: r(l) = a * l^(1/CanfrmPar2)
			// the crown volume is defined as the solid of revolution of this function around the x-axis: V = M_PI * a^2 * Mypow(l,2/CanfrmPar2+1)/(2/CanfrmPar2+1)
			// combining both equations allows to calculate crown radius as a function of crown volume, crown length, and CanfrmPar2
			if (height_ < crown_h_base_)  {calCrownBHeight(); cout<<"calCrownBHeight() workaround in calCrownVolume() "<< id_ << setw(14)<< sd_.getVegType() <<setw(14)<< height_ <<setw(14)<< crown_h_base_ <<endl;}
//			crown_l_ = height_ - crown_h_base_ ;	// crown lenght from top to lower edge	
			
			crown_r_ = sqrt(crown_volume_*(2./sd_.getCanfrmPar2()+1.)/(M_PI*crown_l_)) ; 	// for multistem trees, this is the radius of one sub-crown, as it is calculated based on sub-crown volume	
			
			par_d_ = crown_r_ / Mypow(crown_l_, 1./sd_.getCanfrmPar2()) ; 		// opening width of crown at top for a given crown shape; required for fire routine to calculate radius at any given length of the crown
//	cout << "calcrownvolume: "<<crown_r_<<" "<<crown_volume_<<" "<<stem_diam_int_<<" "<<rho_crown_<<endl;
#ifdef NANCHECK
		if(std::isnan(crown_r_)) cout << "crown radius NaN at calCrownVolume() crown_l_: "<< crown_l_<<"rho_crown_: "<<rho_crown_ <<"crown_volume_ : "<<crown_volume_<<"b_leaf_: "<<b_leaf_<<"b_branch_: "<<b_branch_ <<"height_ "<<height_<<endl;
#endif
	
		return ;
}
 #endif
 
//-----------------------------------------------------------------------------------

/// Returns crown area at height z and stores height in member variables,
/// stores ca at height 0 and at height zstar

void clPlant::calCrownArea()
{
	double crown_area_0_one_ ; // NEW CAMILLE
	
	if(sd_.getVegType()==GR_4)
	{
		double diam_gr = MyMin(0.5+((1.5*b_root_)/(b_root_+0.250)), 1.128); // check max numbers; now correct to obtain a max crown area of 1 m2 for grasses; before was 1.28; REF?
		crown_area_0_ = M_PI*pow((diam_gr/2.), 2.);// less killing of trees? //300315
	}

	else   //if(sd_.getVegType()==TREE)
	{
		#ifdef W_BRANCHES
			
			crown_area_0_ = stem_count_ * crown_r_ * crown_r_ * M_PI ;  		// estimated crown area at the base of the crown; making multiple stems behave fully additive is probably overestimating as there will be 
																				// overlaps between the different sub-crowns, need to come up with something better here
				
		#else
			//	crown_area_0_ = Mypow( sd_.getCanfrmPar1()*height_*( 1.-Mypow( crown_h_base_/height_, sd_.getCanfrmPar2() ) ), 2 )*M_PI;   // SIMON
			// crown_area_0_ = stem_count_*Mypow( sd_.getCanfrmPar1()*height_*( 1.-Mypow( crown_h_base_/height_, sd_.getCanfrmPar2() ) ), 2 )*M_PI;   // SIMON added stem_count_; based on Strigul et al., 2008, eq. 6
			crown_area_0_ = MyMax(0.2, (Mypow( 10.*stem_diam_int_*( 1.-Mypow( crown_h_base_/height_, sd_.getCanfrmPar2() ) ), 2)*M_PI)*stem_count_); // Liam (Strigul-ish) - check with Camille & Simon about calculating can area from stem diam  with this implimentation 
			crown_area_0_one_ = MyMax(0.2, (Mypow( 10.*stem_diam_int_*( 1.-Mypow( crown_h_base_/height_, sd_.getCanfrmPar2() ) ), 2)*M_PI));// NEW CAMILLE
			
			if (crown_area_0_ > crown_area_0_one_) // NEW CAMILLE
			{
				crown_area_0_ = crown_area_0_one_ ;// NEW CAMILLE
			}
			else
			{
				crown_area_0_ = crown_area_0_ ;// NEW CAMILLE
			}
			
		#endif	
	
		
		
	}
	
	//if (RUnif()<0.01) cout << setw(4) << sd_.getVegType() << setw(14) << height_ << setw(14) << crown_h_base_ << setw(14) << sd_.getCanfrmHBase() << setw(14) << crown_area_0_ << endl;

	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calWiltPoint(GridCellSoilData* soil)		// THIS FUNCTION SERIOUSLY CHANGED IN LIAM'S NEW VERSION; AS FAR AS I CAN SEE, THERE IS JUST ONE USEFUL LINE IN THERE RIGHT NOW (wilt_point_=sd._getCanfrmPar())
{														// REST COULD BE COMMENTED OUT OR DELETED JUST AS WELL WITHOUT MAKING ANY DIFFERENCE (?)
//note: SWC = Saturated Water Content
//	float sla  = sd_.getSla();				// FLAG: THIS COMMENTED OUT SINCE THE FUNCTION getSla HAS BEEN BROKEN WITH LIAM INTRODUCING THE P50 METHOD
//	float ar   = sd_.getAllocRoot();		//	min(1,(root biomass/leaf biomass)*scalingfactor(0,1) )  // NEW LIAM COMMENTED OUT
	
	// new function, function is always field capacity for ar=0 and wiltpoint for ar=1 
// 	wilt_point_ = pow( (1.-ar), 0.5*(SLA_MAX[sd_.getVegType()]-sla) ) * (SWC-thetawp)*0.99 + thetawp;		
// 	theta_helper_ = 1./(SWC-thetawp);

// 	wilt_point_ = -4.4610 + 0.1932*sla;	 // Ackerly 2004 with 1 outlier removed	
//	double water_potential_wilt = -5 + 0.1166*sla;		// chosen to give wilt point of -1.5 MPa at SLA of 30m^2/kg, i.e. approx. that of sunflower			// WAS COMMENTED IN IN LIAM'S NEW VERSION, BUT DOESN'T GET USED ANYWHERE
//	double height_adjust_potential = 0.01*height_; 		// wilt point less negative as tree gets taller -- 0.1 MPa per 10 m									// WAS COMMENTED IN IN LIAM'S NEW VERSION, BUT DOESN'T GET USED ANYWHERE
//	water_potential_wilt += height_adjust_potential; 																										// WAS COMMENTED IN IN LIAM'S NEW VERSION, BUT DOESN'T GET USED ANYWHERE		
// 	wilt_point_ = water_potential_wilt;						//NOTE: This is in MPa
// 	wilt_point_ = (SWC)*Mypow(phi_s/(-(water_potential_wilt*100)), (1/r));		// convert back to soil moisture content NOTE: need to change wilt point in water model then
// 	wilt_point_ = (0.534*-2.7)+0.047; //proxy for minimum leaf water potential
	wilt_point_ = sd_.getP50();	//decided against min leaf wat pot and will just use p50 instead
// 	cout << "  smp_WP  " << soil_mat_pot_weighted_ << "  wilt_point_  " << wilt_point_ << endl;	
	
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calLai()
{
// 	lai_ = sd_.getSla()*b_leaf_/MyMax(crown_area_0_,0.1); // 0.1 avoids huge LAIs, particularly after fire
//	double C1 = sd_.getCanfrmPar1();
	
	if(sd_.getVegType()==GR_4)			// NEW LIAM
	{
		lai_ = (length_leaf_total_grass_*(((1.25*sd_.getP50())+22.5)/1000.))/crown_area_0_; 		// NEW LIAM
		lai_ = MyMin( lai_, 10. );
	}
	else
	{
		lai_ = MyMax(0, sla_*b_leaf_/crown_area_0_); // 0.1 avoids huge LAIs, particularly after fire		// NEW LIAM
		lai_ = MyMin( lai_, 10. );
	}
	return;
}

//-----------------------------------------------------------------------------------

// calculates boundary layer conductance in mmol/m^2/s
void clPlant::calBoulCond( double wnd, double mms_helper  )
{
	float Zd;
	float Z0;
	float hbar = height_;
	
	hbar = MyMax(MyMin(hbar,11.5),0.1) ;	// constrain between [0.1, 11.5]; necessary to avoid negative arguments in log
	
	float log_tmp = 1./log(((float)REF_HEIGHT_Z-(float)ZD_CONST*hbar)/(float)Z0_CONST*hbar);
	
	boul_cond_ = KARMAN_CONST_QUAD*wnd*log_tmp*log_tmp;   // in m/s (from wnd), KARMAN and log_tmp dimensionless
	boul_cond_ = boul_cond_/mms_helper;  // now in mumol/m^2/s
	return;
}

//-----------------------------------------------------------------------------------

// calculates stomatal conductance in mmol/m^2/s

void clPlant::calStomCond( double tmp, double apr, double co2, double hum, double As, double As_max )
{
// 	tmp    deg C
// 	apr    Pa
//	co2    Pa
//	hum    %
//	As     mumol/m^2/s
	double cs = 0;
	double cs_max_ = 0;				// used to calculate potential evapotranspiration
	
	if ( active_ )
	{
		cs = co2-CS_FACTOR*As*apr/boul_cond_;	// in Pa, this line follows Arora
		cs = MyMax( cs, 0.1*co2 );  // NOTE cs must be positive and if it is too small, we get problems with too high stom. cond as we divide by cs
									// NOTE with the factor 0.1 I observed max stom. cond. of 1 mol/m^2/s at tapajos (eucalypts: up to 0.8 mol/m^2/s)

		cs_max_ = co2-CS_FACTOR*As_max*apr/boul_cond_;	// in Pa, this line follows Arora
		cs_max_ = MyMax( cs_max_, 0.1*co2 );  // NOTE cs must be positive and if it is too small, we get problems with too high stom. cond as we divide by cs
			
		stom_cond_ = sd_.getBallBerryM()*(1.-sd_.getRMaintResp())*As*hum*apr/cs+sd_.getBallBerryB(); // in mumol/m^2/s

		stom_cond_max_ = sd_.getBallBerryM()*(1.-sd_.getRMaintResp())*As_max*hum*apr/cs_max_+sd_.getBallBerryB(); // in mumol/m^2/s

	}
	else
	{
		stom_cond_ = 0.;
		stom_cond_max_ = 0.;
	}
	
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calEvapotr( double et_helper_0, double et_helper_1, double et_helper_2, double et_helper_3, double mms_helper )
{
// 	The following helpers are defined in PlantPopClass.h
// 	et_helper_0 = (s*rad)
// 	et_helper_1 = (86400.*rho*SP_HEAT*vpd);
// 	et_helper_2 = LAMBDA*s;
// 	et_helper_3 = LAMBDA*gama;
// 	mms_helper  = 1e-6*0.0224*((tmp+273.15)/273.16)*(apr/1.013e5);  // transformation mmol/m^2/s to m/s

	double stom_cond_tmp = stom_cond_*mms_helper; // transform to non-molar units (m/s)
	double boul_cond_tmp = boul_cond_*mms_helper; // transform to non-molar units (m/s)

	double stom_cond_max_tmp = stom_cond_max_*mms_helper; // transform to non-molar units (m/s)
	
	evapotr_  = (et_helper_0 + et_helper_1*boul_cond_tmp)/( et_helper_2 + et_helper_3 + et_helper_3*boul_cond_tmp/stom_cond_tmp );  // in mm/day
	evapotr_ = MyMax(0, evapotr_*crown_area_0_*light_avail_) ;

	evapotr_max_  = (et_helper_0 + et_helper_1*boul_cond_tmp)/( et_helper_2 + et_helper_3 + et_helper_3*boul_cond_tmp/stom_cond_max_tmp );  // in mm/day
	evapotr_max_ = evapotr_max_ * crown_area_0_ * light_avail_ ;	// FLAG MP: somehwere in-between of using crown area (underestimate!) and total leaf area (overestimate); not sure how good this really is...

	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calLeafTemp( double hum, double eS, double tmp, double apr, double mms_helper )
{
	double leaf_helper = 0;
	
	double stom_cond_tmp = stom_cond_*mms_helper; // transform to non-molar units (m/s)
	double boul_cond_tmp = boul_cond_*mms_helper; // transform to non-molar units (m/s)
	
	if ( active_==0 ) leaf_temp_ = tmp;
	else
	{
		leaf_helper =  log( (evapotr_*(1./(stom_cond_tmp*1000.)+1./(boul_cond_tmp*1000.) ) + hum*eS)/0.6108 );
		leaf_temp_ = 237.3*leaf_helper/(17.27-leaf_helper);
	}
// 	if ( RUnif()<0.001 ) cout << setw(14) << tmp << setw(14) << leaf_temp_ << setw(14) << leaf_temp_-tmp << endl;
	
	return;
}

//-----------------------------------------------------------------------------------

void clPlant::calLightAvail()
{
	// light ext small -> has light availability even at small LAIs
	// light ext high  -> requires minimum LAI but light_avail_ much higher at higher LAI
	//light_avail_    = (1.-Myexp(-sd_.getLightExt()*lai_))/sd_.getLightExt();
	
	light_avail_    = ((1.-Myexp(-light_ext_*lai_))/light_ext_);  // SIMON

	//light_avail_    = ((1.-Myexp(-sd_.getLightExt()*lai_))/sd_.getLightExt())*(1.-(0.5*crown_h_base_/height_));  // CAM: "*(1-(0.5*can8a_/height_));" new part of the equation; 0.5 is an arbitrary; SIMON second term assumes that in a small crown (high crown_h_b_) leaves are packed densly and light availability goes down due to self shading, in contrast to a large crown (crown_h_base_ close to zero).

	//if (RUnif()<0.01) cout << setw(4) << sd_.getVegType() << setw(14) << (1.-Myexp(-sd_.getLightExt()*lai_))/sd_.getLightExt() << setw(14) << light_avail_ << setw(14) << crown_h_base_ << setw(14) << sd_.getCanfrmHBase() << setw(14) << crown_area_0_ << endl;


	return;
}

//-----------------------------------------------------------------------------------

// void clPlant::runAnnualProcesses( double z_star, double thetaWP, double thetaFC )
// void clPlant::runAnnualProcesses( double ca_frac, double z_star )

void clPlant::runAnnualProcesses(/* double *compheight */)
{
// 	if ( c_gain_cum_<0 && active_==1) c_gain_cum_flag_= 1;
	if ( c_gain_cum_<0 ) c_gain_cum_flag_= 1;	// if(id_==99) cout << "C-flag: " << c_gain_cum_flag_ << endl;
	c_gain_cum_ = 0.;
	day_ = 0.;		// NEW LIAM, WAS COMMENTED IN BEFORE
	
	evg_flush_flag_ = 0; // NEW LIAM 121115	
	
	age_++;
	
	return;
}

//-----------------------------------------------------------------------------------


void clPlant::runDecomposition()
{
	// transition from standing/hanging dead leaves to dead leaves lying on ground
//	b_leaf_d_ly_ += 0.00001*b_leaf_d_st_;
//	b_leaf_d_st_ -= 0.00001*b_leaf_d_st_;

	if ( sd_.getVegType()==TREE )
	{
		b_leaf_d_ly_ += 0.15*b_leaf_d_st_;  // 0.15: reduction to ca 10% after 2 weeks
		b_leaf_d_st_ -= 0.15*b_leaf_d_st_;
	}
	else
	{
		b_leaf_d_ly_ += 0.005*b_leaf_d_st_;  // 0.005: reduction to ca 16% after one year
		b_leaf_d_st_ -= 0.005*b_leaf_d_st_;
	}

	// transition from dead leaves lying on the ground to soil pool
	soilcarbon_nwl_ += 0.001*b_leaf_d_ly_;
	b_leaf_d_ly_    -= 0.001*b_leaf_d_ly_;	

	// decomposition of dead material
	soilcarbon_nwl_ *= 0.9999;  // SIMON fine roots and leaves
	soilcarbon_fwl_ *= 0.9999;  // SIMON coarse roots, branches
	soilcarbon_cwl_ *= 0.9999;  // SIMON coars stem 

	return;
}




//-----------------------------------------------------------------------------------

void clPlant::runDailyProcesses( int year, int dct, GridCellSoilData* soil, double A0C3, double A0C4,  double wnd, double pre,
						double tmp, double apr, double co2, double hum, double sun_dur, double rad, double resp_temp_fac,
						double et_helper_0, double et_helper_1, double et_helper_2, double et_helper_3, double eS, double mms_helper,
						double *compheight, double *compheight_lai, double *compheight_crown_area )					
								 
{   // A0C3 and A0C4 given in mumol/m^2/s CO2
	
	double As;	 	// limited crown photosynthesis
	double As_max;	// non-limited crown photosynthesis
	double As_net;	// limited crown photosynthesis - growth respiration
	double RmLeaf;	// leaf maintainance respiration
	double RmWood;	// aboveground woody biomass respiration (stem + branches)
	#ifdef W_BRANCHES
		double RmBranch; // branch maintenance respiration
	#endif	
	double RmStem;	// stem maintainance respiration
	double RmRoot;	// root maintainance respiration
	double RmStor;	// storage maintainance respiration
	double RmBark;	// bark maintenance respiration
	double RmRepr;
	
	double Rm_Sum_ = 0;
	
	double height_comp = crown_h_base_*can_frm_helper_;
	
	if(dct_==0) {sum_active_ = 0 ; sum_light_ = 0. ; sum_water_ = 0.; }
	if(active_==1) sum_active_++ ;
	
	dct_ = dct ;
	year_ = year ;

//	if(id_==1287 && sd_.getVegType()==GR_4 && (year_==5 || year_==6)) cout << "Begin daily: " << year_ <<setw(14)<< dct_ <<setw(14)<< active_ <<setw(14)<< sd_.getEvergreen() <<setw(14)<< b_leaf_ <<setw(14)<< b_root_ <<setw(14)<< b_stor_ <<setw(14)<< b_repr_<<setw(14)<< c_balance_ <<setw(14)<< c_gain_cum_ <<setw(14)<< age_ << endl;
//	if(b_repr_<0) cout << "Begin daily: " << year_ <<setw(14)<< dct_ <<setw(14)<< active_ <<setw(14)<< sd_.getEvergreen() <<setw(14)<< b_leaf_ <<setw(14)<< b_root_ <<setw(14)<< b_stor_ <<setw(14)<< b_repr_<<setw(14)<< c_balance_ <<setw(14)<< c_gain_cum_ <<setw(14)<< age_ << endl;

		
//	if (year==1 && id_==1) cout << "POS1: " << setw(14) << dct_ << setw(14) << b_stem_ << setw(14) << b_leaf_ << setw(14) << height_ << setw(14) << alive_ << setw(14) << active_ << endl; 
//	if(id_==99) cout << year_ <<setw(14)<< dct_ <<setw(14)<< active_ <<setw(14)<< sd_.getEvergreen() <<setw(14)<< b_stem_ <<setw(14)<< b_leaf_ <<setw(14)<< b_bark_ <<setw(14)<< b_root_ <<setw(14)<< b_repr_ <<setw(14)<< b_stor_
//					 <<setw(14)<< c_balance_ <<setw(14)<< c_gain_cum_ <<setw(14)<< sum_active_ << endl;
//	if(dct_==364 && c_gain_cum_ <= 0. && year_ > 90) cout << year_ <<setw(14)<< id_ <<setw(14)<< sd_.getEvergreen() <<setw(14)<< c_gain_cum_ <<setw(14)<< sum_active_ <<setw(14)<< height_ << endl;
//	if(dct_==364 && year_==999) cout << id_ <<setw(14)<< sd_.getEvergreen() <<setw(14)<< sum_active_ <<setw(14)<< height_ <<setw(14)<< age_ <<setw(14)<< crown_area_0_ <<setw(14)<< lai_ <<setw(14)<< sum_water_/(dct_+1.) <<setw(14)<< sum_light_/(dct_+1.) << endl;
 
	// only for debugging: previous day's biomass values
	#ifdef W_BRANCHES
		b_leaf_old_ = b_leaf_ ;
		b_branch_old_ = b_branch_ ;
		b_stem_old_ = b_stem_ ;
		b_bark_old_ = b_bark_ ;
		b_root_old_ = b_root_ ;
		b_stor_old_ = b_stor_ ;
		b_repr_old_ = b_repr_ ;
	#endif	 
 	
	if ( sd_.getVegType()==GR_4 )  As = A0C4; // mumol/m^2/s
	else                /* TREE */ As = A0C3; // mumol/m^2/s
	
	calWiltPoint( soil );
//	if(year_==16 && id_> 1597) cout << "Before calWaterAvailability: " << year_ <<setw(14)<< dct_ <<setw(14)<< id_ << endl; 
	calWaterAvailability( soil, year );   // calculate water availability index water_avail_
//	if(year_==16 && id_>1597) cout << "After calWaterAvailability: " << year_ <<setw(14)<< dct_ <<setw(14)<< id_ << endl;
	calLightAvail();
	calLai();	// cout << "Height 1: " << year << setw(14) << id_ << setw(14) << height_*100. << endl;
	calHeight();	// cout << "Height 2: " << year << setw(14) << id_ << setw(14) << height_*100. << endl;
	#ifndef W_BRANCHES
		calCrownBHeight(); 		// with branches, the upward shift of crown base height happens due to adjustment of crown length in height calculation
	#endif	
//	calStemDiam();		// now called in calHeight!
	calCrownArea();	
	
	for ( int i=0; i<NUM_LIGHT_NGB; i++ )
	{
//		if ( height_comp < compheight[i] ) light_avail_ *= (LIGHT_COMP_PAR*(height_comp/compheight[i])+(1.-LIGHT_COMP_PAR) );		// BEFORE LIAM CHANGE

// 		if ( height_comp < compheight[i] && compheight_crown_area[i] > 6.25 && compheight_lai[i] > 1.0) // Liam: the 6.25 based on size of individual cell so fn. of plot size and no. individuals, the lai is also changable
// 		if ( height_comp < compheight[i] && compheight[i] > 2.0 && compheight_lai[i] > 0.0 )
// 		{
// 			light_avail_ *= Mypow((1.-LIGHT_COMP_PAR), compheight_lai[i] );
// 			//cout << setw(3) << sd_.getVegType() << setw(12) << Mypow((1.-LIGHT_COMP_PAR), compheight_lai[i] ) << endl;
// // 			light_avail_ *= Mypow((1.-LIGHT_COMP_PAR)*(1-(height_comp/compheight[i])), (compheight_lai[i])*MyMin(1.,(compheight_crown_area[i]/60.25)) );	// Liam: I'm working on this at the moment so will updata again
// 			
// 		}
//NOTE: I will need to come back to this. 

		if ( height_comp < compheight[i] && compheight_crown_area[i] > 2.77 && compheight_lai[i] > 0.1 /*&& sd_.getVegType()==0*/) 
		{		  
			light_avail_ *= Mypow((1.-LIGHT_COMP_PAR)*(1-(height_comp/compheight[i])), (compheight_lai[i]*1.)*MyMin(1.,(compheight_crown_area[i]/25.0)) );	// 151015	  			
		}

//		if ( height_comp < compheight[i] && compheight[i] > 2.0 && compheight_lai[i] > 0.0 )
//			light_avail_ *= (LIGHT_COMP_PAR*(height_comp/compheight[i])+Mypow((1.-LIGHT_COMP_PAR), compheight_lai[i]) );	// NEW LIAM		
	}


	if (sd_.getVegType ()== GR_4)                        // SIMON CAMILLE grass self shading
	{                                                    // SIMON CAMILLE grass self shading
		if ( b_leaf_ < b_leaf_d_st_ )                // SIMON CAMILLE grass self shading
		{
			light_avail_ *= (LIGHT_COMP_SELF_DEAD*b_leaf_/b_leaf_d_st_ + (1.-LIGHT_COMP_SELF_DEAD));  // SIMON CAMILLE grass self shading
		}
	}
	
	As_max = As * light_avail_ * active_;
	As = As * water_avail_ * light_avail_ * active_;   // water and light stressed photosynthesis  // mumol/m^2/s
	
//	cout << "water_avail_: " << water_avail_ << endl;

	calBoulCond( wnd, mms_helper );
	calStomCond( tmp, apr, co2, hum, As, As_max );  // requires water stressed A in mumol C/m^2/s
	calEvapotr( et_helper_0, et_helper_1, et_helper_2, et_helper_3, mms_helper );
// 	calLeafTemp( hum, eS, tmp, apr, mms_helper );

	As  = As * crown_area_0_;	// water&light stressed crown photos., mumol/s per plant, only when active 			// FLAG MP: crown area vs. total leaf area of plant?
	As_max  = As_max * crown_area_0_;	// water&light stressed crown photos., mumol/s per plant, only when active 
	
	c_gain_mmolm2s_ = As; // gpp in mumol C/plant/s     NOTE: per plant or per m^2?
	
	As      = As * sun_dur * MMS_TO_KGD_HELPER;	// convert As from micromol CO2/m^2/s to kg C/day/m^2 and into kg biomass/day/m^2
	c_gain_ = As;  // in kg biomass/day/plant
	
	// maintenance respiration
	RmLeaf = ( leaf_resp_helper_ * resp_temp_fac * b_leaf_ );
	RmRepr = ( leaf_resp_helper_ * resp_temp_fac * b_repr_ );
	#ifdef W_BRANCHES
		RmBranch = ( wood_resp_helper_ * resp_temp_fac * b_branch_ );					// MP: don't know what to do about sapwood respiration of branches, cause branch diameters unknown?
		RmStem = ( wood_resp_helper_ * resp_temp_fac * b_stem_ /* sap_frac_*/ );		// MP: only sapwood respires, heartwood is dead and does not respire
		RmWood = RmBranch + RmStem ;
	#else	
		RmStem = (wood_resp_helper_ * resp_temp_fac * b_wood_ /* sap_frac_*/ ); 
		RmWood = RmStem ;
	#endif	
	RmBark = ( wood_resp_helper_ * resp_temp_fac * b_bark_ );	// MP: not sure if bark should respire, or if it is all dead? Cambium would be alive and respire, but isn't that included in RmStem?
	RmRoot = ( root_resp_helper_ * resp_temp_fac * b_root_ );
	RmStor = ( root_resp_helper_ * resp_temp_fac * b_stor_ );
	
	double RmTot = RmLeaf + RmWood + RmRepr + RmBark + RmRoot + RmStor ; 	// total overall  maintenance respiration 
	
//	if(sd_.getVegType()==GR_4) cout << "Maintenance respiration: " << id_ <<setw(14)<< RmTot <<setw(14)<< RmLeaf <<setw(14)<< b_stem_ <<setw(14)<< b_branch_ <<setw(14)<< RmRepr <<setw(14)<< RmBark << setw(14)<< RmRoot <<setw(14)<< RmStor << endl;

// MP: instead of taking maintenance respiration from the biomass pools, now taking it from net assimilation; using c_balance as the total remaining C-gain that will be allocated to the different biomass compartments
// if C-gain is negative, i.e., total respiratory expenditure exceeds C-gain through assimilation, take away the remainder from the living biomass pools to get to zero.
/*	b_leaf_ -= RmLeaf;
	b_repr_ -= RmRepr;
	b_wood_ -= RmWood;
	#ifdef W_BRANCHES
		b_branch_ -= RmBranch;
	#endif	
	b_stem_ -= RmStem;
	b_bark_ -= RmBark;
	b_root_ -= RmRoot;
	b_stor_ -= RmStor;
*/	
	
// 	Calculate carbon gain, given in kg biomass/day/plant
	As_net     = As * ( 1.-sd_.getRGrowthResp() );   // this is always > 0 as Rg= const*As

	//Calculate the net carbon gain of plant, that is As - growth resp - maintenance resp - turnover, this can be <0
	c_balance_ = As_net - RmTot - leaf_longevity_*b_leaf_;  // MP: this is what remains for the plant to be allocated to the respective biomass pools, if it's > 0

	// leaf turnover
	double loss = leaf_longevity_*b_leaf_;
	b_leaf_d_st_ += leaf_longevity_*b_leaf_;
	b_leaf_      -= leaf_longevity_*b_leaf_;	// if(id_==6) cout << "Leaf longevity: " << b_leaf_*leaf_longevity_ << endl;	
	
	// root turnover
	if(sd_.getVegType()==GR_4) b_root_ -= root_turnover_*b_root_;   // SIMON also for trees? Only fine roots = certain proportion of total roots
	if(sd_.getVegType()==GR_4) c_balance_ -= root_turnover_*b_root_;		// MP: should we not also have fine root turnover for trees in some way?
// 	if(sd_.getVegType()==TREE) b_root_ -= root_turnover_tree_*b_root_;   // SIMON also for trees? Only fine roots = certain proportion of total roots	

	if ( soil_mat_pot_weighted_ <= wilt_point_ && active_==1) 	// NEW LIAM; plants are experiencing water stress
	{
  		day_++;
		
		if (day_ > 30)
		  {
		    c_balance_   -= 0.02*b_leaf_;
		    c_balance_   -= 0.02*b_root_;
		    b_leaf_d_st_ += 0.02*b_leaf_;		// NEW LIAM
		    b_leaf_       = 0.98*b_leaf_;		// NEW LIAM
		    b_root_       = 0.98*b_root_;		// NEW LIAM
		}
	}

	if ( soil_mat_pot_weighted_ > wilt_point_ ) day_--;
	if ( day_<0 ) day_=0;

//	if (year==1 && id_==1) cout << "POS6: " << setw(14) << dct_ << setw(14) << b_stem_ << setw(14) << b_leaf_ << setw(14) << height_ << setw(14) << alive_ << setw(14) << active_ << endl;

	c_gain_cum_ += c_balance_;
	
	// allocate carbon to biomass pools
	if ( c_balance_ > 0 )		// MP: changed this from As_net to c_balance	FLAG: if negative, need to take carbon away from biomass pools to get to zero
	{
	
//		double alloc_leaf_dynamic = sd_.getAllocLeaf()/(1.+Myexp(3.*(lai_-7.)));
//		if (lai_ > 6.) alloc_leaf_dynamic = 0.0; 
//		double alloc_leaf_dynamic = sd_.getAllocLeaf()/(1.+Myexp(1.5*(lai_-6.)));
		double alloc_leaf_dynamic = sd_.getAllocLeaf()/(1.+Myexp(3*(lai_-4.)));
		double alloc_root_dynamic = sd_.getAllocRoot();
		double alloc_wood_dynamic = sd_.getAllocWood();
		double alloc_bark_dynamic = sd_.getAllocBark();
		double alloc_repr_dynamic = sd_.getAllocRepr();
		double alloc_stor_dynamic = sd_.getAllocStor();
		double alloc_stem_dynamic = sd_.getAllocStem();		// if not with branches, this is equivalent to alloc_wood_dynamic (taken care of in SeedClass)
		#ifdef W_BRANCHES
			double alloc_branch_dynamic = sd_.getAllocBranch();
		#endif

		//redistribute excess leaf allocation to stem
		//--------------
		
		double diff_leaf_alloc = sd_.getAllocLeaf()-alloc_leaf_dynamic; 		// NEW LIAM: what does not go to leaves, goes to stem instead		
		if(sd_.getVegType()==TREE)
		{
			alloc_stem_dynamic += diff_leaf_alloc; 			// only stem benefits from excess allocation
			#ifdef W_BRANCHES
				alloc_wood_dynamic = alloc_stem_dynamic + alloc_branch_dynamic; 	// need to update wood dynamic as it contains stems!
			#else
				alloc_wood_dynamic = alloc_stem_dynamic; 		// no branches, therefore wood = stem, excess goes to wood=stem
			#endif	
		}	
		//--------------
		
		
		double max_b_stor=50.;
		switch(sd_.getVegType())
		{
		case TREE: 
			max_b_stor=50.; //Trees
			break;
		case GR_4:
			max_b_stor=0.5; //C4 Grasses
			break;
		default:
			max_b_stor=50.;
			cout<<"Warning in PlantClass rundailyprocesses allocation undefined vegtype occurred"<<endl;
			break;
		}

		// Stop all allocation to storage when b_stor_ exceeds 50% of root biomass or max_b_stor value
		if((b_stor_ > b_root_*0.5 || b_stor_ > max_b_stor))	alloc_stor_dynamic = 0.;	
		
		//normalize all alloc_xxxx_dynamic variables (works for grasses as well as trees; for grasses, wood and bark are zero, so still fine with the sum to include them)
		double sum_dynamic = alloc_root_dynamic + alloc_leaf_dynamic + alloc_wood_dynamic + alloc_bark_dynamic + alloc_repr_dynamic + alloc_stor_dynamic; 
		alloc_root_dynamic /= sum_dynamic;
		alloc_leaf_dynamic /= sum_dynamic;
		alloc_wood_dynamic /= sum_dynamic;
		alloc_bark_dynamic /= sum_dynamic;
		alloc_repr_dynamic /= sum_dynamic;
		alloc_stor_dynamic /= sum_dynamic;
		alloc_stem_dynamic /= sum_dynamic;
		#ifdef W_BRANCHES
			alloc_branch_dynamic /= sum_dynamic;  //Notice: sum_dynamic contains already stem+branch in form of wood
		#endif

		//--------------- 
				
		b_leaf_ += c_balance_ * alloc_leaf_dynamic; 	// if(id_==6) cout << "Leaf allocation: " << c_balance_ * alloc_leaf_dynamic << endl;
		b_wood_ += c_balance_ * alloc_wood_dynamic;		// MP
		#ifdef W_BRANCHES
			b_branch_ += c_balance_ * alloc_branch_dynamic;	// MP
		#endif
		b_stem_ += c_balance_ * alloc_stem_dynamic;		// MP
		b_root_ += c_balance_ * alloc_root_dynamic;		//cout << "Root allocation: " << b_root_ <<setw(14)<< c_balance_ <<setw(14)<< alloc_root_dynamic << endl;
		b_bark_ += c_balance_ * alloc_bark_dynamic;
		b_repr_ = MyMax(0.,b_repr_+c_balance_ * alloc_repr_dynamic);
		b_stor_ += c_balance_ * alloc_stor_dynamic;
		
//	if(id_==460) cout << year_ <<setw(14)<< dct_ <<setw(14)<< As_net <<setw(14)<< As_net*alloc_stem_dynamic-RmStem <<setw(14)<< As_net*alloc_leaf_dynamic-RmLeaf-loss <<setw(14)<< As_net*alloc_bark_dynamic-RmBark <<setw(14)<< As_net*alloc_root_dynamic-RmRoot
//					<<setw(14)<< As_net*alloc_repr_dynamic-RmRepr <<setw(14)<< As_net*alloc_stor_dynamic-RmStor <<setw(14)<< c_balance_ << endl;
			
	}
	else // MP: c_balance_ <= 0, no allocation, biomass removal instead	
	{
		b_leaf_ = MyMax(b_leaf_ + RmLeaf/RmTot * c_balance_, 0.) ;	// if no leaves, RmLeaf is zero
		b_repr_ = MyMax(b_repr_ + RmRepr/RmTot	* c_balance_, 0.) ;	// + because c_balance is negative or zero when we end up in this case!
		#ifdef W_BRANCHES
			b_branch_ = MyMax(b_branch_ + RmBranch/RmTot * c_balance_, 0.) ;
		#endif	
		b_stem_ += RmStem/RmTot * c_balance_ ;
		b_wood_ += RmWood/RmTot * c_balance_ ;
		b_bark_ += RmBark/RmTot * c_balance_ ;
		b_root_ += RmRoot/RmTot * c_balance_ ;		// cout << "Root allocation: " << b_root_ <<setw(14)<< c_balance_ <<setw(14)<< RmRoot <<setw(14)<< RmTot << endl;
		b_stor_ = MyMax(b_stor_ + RmStor/RmTot * c_balance_, 0.) ; 

		// MP: check biomass pools crucial for immediate survival don't get <= 0; if any does, kill that individual

//		if(veg_type_==TREE && (b_stem_ <= 0. || b_wood_ <= 0. || b_bark_ <= 0. || b_root_ <= 0.)) setDead();		// doesn't seem to happen	FLAG: decrease pop_size in PlantPopClass!!! doing it that way is not correct!
//		if((veg_type_==GR_C4) && b_root <= 0.) setDead() ;	
	}

//	if ( water_avail_ < wilt_avg_ && active_==1) day_++;		// BEFORE LIAM CHANGE	
//	calLai();													// BEFORE LIAM CHANGE

	// FLAG: this moving window average is TERRIBLY inefficient
	for ( int i=1; i<N_PHEN; i++ )
	{
		phen_rain_[i-1] = phen_rain_[i];
		phen_lght_[i-1] = phen_lght_[i];
	}
//	phen_rain_[N_PHEN-1] = pre;
	phen_rain_[N_PHEN-1] = soil_mat_pot_weighted_;		// NEW LIAM		
	phen_lght_[N_PHEN-1] = rad;
	
	double phen_rain_mean = 0.;
	double phen_lght_mean = 0.;
	
	// FLAG: Bug: here there is a sum over 13 elements only, but the mean is calculated by dividing by 14; fixed N_PHEN to 13...
	for ( int i=1; i<N_PHEN; i++ ) 
	{
		phen_rain_mean += phen_rain_[i];
		phen_lght_mean += phen_lght_[i];
	}
	phen_rain_mean /= (double)N_PHEN;
	phen_lght_mean /= (double)N_PHEN;
	
	// phenology control
	
//	if(id_==1 && year==1) cout << "Phen stats: " << sd_.getEvergreen() << setw(14) << sd_.getRainLight() << setw(14) << phen_rain_mean << setw(14) << sd_.getRainThrOn() << setw(14) << sd_.getRainThrOff() << setw(14) << phen_lght_mean << setw(14) << sd_.getLightThrOn() << setw(14) 
//			<< sd_.getLightThrOff() << setw(14) << endl;
	
	if ( sd_.getEvergreen() <= EVERGREEN_THR )   // only deciduous vegetation switches
	{
		if( sd_.getRainLight() <= RAIN_LIGHT_THR ) // rain control
		{
			if( active_==0 && phen_rain_mean > sd_.getRainThrOn() )
			{   // move to metabolic state
//				b_leaf_      += sd_.getStorToLeaf()*b_stor_;		// BEFORE LIAM CHANGE
//				b_stor_      -= sd_.getStorToLeaf()*b_stor_;		// BEFORE LIAM CHANGE
//				b_leaf_      += 0.65*b_stor_;		// NEW LIAM
			  	if(lai_ < 7 && b_leaf_ < 75.)
				{
					b_leaf_      += sd_.getStorToLeaf()*0.65*b_stor_;		// NEW LIAM
// 					b_stor_       = 0.;					// NEW LIAM
					b_stor_      -= sd_.getStorToLeaf()*b_stor_;	
				}
				active_       = 1.;
				dpos_         = 0;
				dneg_         = 0;
				dleaf_alloc_  = 0;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;		// MP: crown volume and radius increase due to adding leaf biomass
				#endif	
//				day_ = 0.;				// NEW LIAM
			}
			if( active_==1 && phen_rain_mean < sd_.getRainThrOff() )	// MP: made sure in SeedClass that RainThrOff is never > RainThrOn
			{   // move to dormant state
//				b_leaf_d_st_ += 1.*b_leaf_;		// BEFORE LIAM CHANGE
//				b_leaf_       = 0.*b_leaf_;		// BEFORE LIAM CHANGE
				if(lai_ < 7 && b_leaf_ < 75.)
				{			  
					b_leaf_d_st_ += 0.7*b_leaf_;	// NEW LIAM
					b_stor_      += 0.3*b_leaf_;	// NEW LIAM
					b_leaf_       = 0.*b_leaf_;		// NEW LIAM
				}
				active_       = 0.;
				dpos_         = 0;
				dneg_         = 0;
				dleaf_alloc_  = 0;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;		// MP: crown volume and radius decrease due to leaf loss
				#endif	
//				day_ = 0.;				// NEW LIAM
			}
		}
	else  // getRainLight>RAIN_LIGHT_THR, light control
		{
			if( active_==0 && phen_lght_mean  > sd_.getLightThrOn() )
			{   // move to metabolic state
//				b_leaf_      += sd_.getStorToLeaf()*b_stor_;		// BEFORE LIAM CHANGE
//				b_stor_      -= sd_.getStorToLeaf()*b_stor_;		// BEFORE LIAM CHANGE
				if(lai_ < 7 && b_leaf_ < 75.)
				{			  
					b_leaf_      += sd_.getStorToLeaf()*0.65*b_stor_;		// NEW LIAM
					b_stor_      -= sd_.getStorToLeaf()*b_stor_;
				}				
				active_       = 1.;
				dpos_         = 0;
				dneg_         = 0;
				dleaf_alloc_  = 0;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;	// MP: crown volume and radius increase due to leaf biomass increase
				#endif	
// 				day_ = 0.;				// NEW LIAM
			}
			if( active_==1 && phen_lght_mean < sd_.getLightThrOff() )
			{   // move to dormant state
//				b_leaf_d_st_ += 1.*b_leaf_;		// BFORE LIAM CHANGE
//				b_leaf_       = 0.*b_leaf_;		// BEFORE LIAM CHANGE
				b_leaf_d_st_ += 0.7*b_leaf_;	// NEW LIAM
				b_stor_      += 0.3*b_leaf_;	// NEW LIAM
				b_leaf_       = 0.*b_leaf_;		// NEW LIAM
				active_       = 0.;
				dpos_         = 0;
				dneg_         = 0;
				dleaf_alloc_  = 0;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;	// MP: crown volume and radius decrease due to loss of leaf biomass
				#endif	
// 				day_ = 0.;				// NEW LIAM
			}
		}
	}
	else  // evergreen
	{
		if( sd_.getRainLight() <= RAIN_LIGHT_THR ) // rain control
		{
			if( evg_flush_flag_==0 && phen_rain_mean > sd_.getRainThrOn() )
			{   // move to metabolic state
// 				b_leaf_        += sd_.getStorToLeaf()*b_stor_;
// 				b_stor_        -= sd_.getStorToLeaf()*b_stor_;
//				b_leaf_      += 0.4*b_stor_;	// NEW LIAM
				if(lai_ < 7 && b_leaf_ < 75)
				{			  
			  		b_leaf_      += sd_.getStorToLeaf()*0.65*b_stor_;				
					b_stor_      -= sd_.getStorToLeaf()*b_stor_;		
				}
				dleaf_alloc_    = 0;
				evg_flush_flag_ = 1;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;	// MP: crown gets larger as new leafs are added
				#endif	
			}
			
			if( evg_flush_flag_==1 && phen_rain_mean < sd_.getRainThrOff() )
			{   // move to dormant state
				evg_flush_flag_ = 0;					
			}
		}
		else  // getRainLight>RAIN_LIGHT_THR, light control
		{
			if( evg_flush_flag_==0 && phen_lght_mean  > sd_.getLightThrOn() )
			{   // move to metabolic state
// 				b_leaf_        += sd_.getStorToLeaf()*b_stor_;
// 				b_stor_        -= sd_.getStorToLeaf()*b_stor_;
//				b_leaf_      += 0.4*b_stor_;		// NEW LIAM
				if(lai_ < 7 && b_leaf_ < 75)
				{			  
					b_leaf_        += sd_.getStorToLeaf()*b_stor_*0.65;				
					b_stor_        -= sd_.getStorToLeaf()*b_stor_;		
				}
				dleaf_alloc_    = 0;
				evg_flush_flag_ = 1;
				#ifdef W_BRANCHES
					if(sd_.getVegType()==TREE) calCrownVolume() ;	// MP: crown gets larger as new leafs are added
				#endif					
			}
			if( evg_flush_flag_==1 && phen_lght_mean < sd_.getLightThrOff() )
			{   // move to dormant state
				evg_flush_flag_ = 0;					
			}
		}
		active_ = 1.;
	}
	
//	if(active_==0.) cout << id_ << setw(14) << dct_ << setw(14) << sd_.getEvergreen() << setw(14) << sd_.getLightThrOn() << setw(14) << sd_.getLightThrOff() << setw(14) << sd_.getRainThrOn() << setw(14) << sd_.getRainThrOff() << endl; 
//	if (id_==28) cout << "End daily: " << dct_ << setw(14) << active_ << setw(14) << sd_.getEvergreen() << setw(14) << As << setw(14) << As_net << setw(14) << c_balance_ << setw(14) << c_gain_cum_ << setw(14) << height_ << setw(14) 
//				<< b_stem_ << setw(14) << b_leaf_ << endl; 
//	if (id_==28) cout << "End daily:   " << year_ << setw(14) << dct_<< setw(14) << id_ << setw(14) << active_ << setw(14) << sd_.getEvergreen() << setw(14) << height_ << setw(14) << b_stem_ << setw(14) << b_leaf_ << setw(14) << c_balance_ << endl; 
	
	//-------- 
	
	for (int i=0; i<N_SOIL_LAYERS; i++ )
	  {
	    root_bm_[i] = b_root_*root_fracs_[i];
	  }
	  
//	if (year==1 && id_==1) cout << "POS8: " << setw(14) << dct_ << setw(14) << b_stem_ << setw(14) << b_leaf_ << setw(14) << height_ << setw(14) << alive_ << setw(14) << active_ << endl;	

	sum_water_ += water_avail_ ;
	sum_light_ += light_avail_ ;
		
	return;
}

//-----------------------------------------------------------------------------------

int clPlant::getSeedProd()
{
	int seed_num = 0;

	if(sd_.getVegType()==TREE) 
		{
			seed_num = floor( ( b_repr_/sd_.getSeedWeight()) );
			b_repr_ -= ( seed_num*sd_.getSeedWeight() );
		}
	else 
		{
			seed_num = floor( ( b_repr_/sd_.getSeedWeight()) )* IND_AREA_MULT;  // MP: * IND_AREA_MULT to get from one clone to total of that grass patch!
			b_repr_ -= ( seed_num/IND_AREA_MULT*sd_.getSeedWeight() );			// b_repr_ is for one clone therefore need to divide by IND_AREA_MULT !!! This was wrong before and caused b_repr to go negative!!!
		}
	if (b_repr_<0.) 	
		{
			cout<<"Warning getSeedProd(): b_repr_<0"<<endl;
			b_repr_=0.;
		}
// 	cout << "SEEDS " << "    Height    " << height_ << "   Seed Num   "  << seed_num << "   b_repr   " << b_repr_ << "   c_gain_cum   " <<  c_gain_cum_ << "  dleaf_alloc   " << dleaf_alloc_ << "    c_gain_cum_flag   " << c_gain_cum_flag_ << "   Active   " << active_  << endl;
	
	
// 	seed_num     = MyMin( seed_num, MAX_SEED_NUM );
// 	b_repr_ -= (seed_num*sd_.getSeedWeight());
	
	
	return seed_num;
//	return 0;
}

//-----------------------------------------------------------------------------------

int clPlant::runMortalityDaily()		// FLAG: this function currently out of use
{
//  if ( day_>90 )		// NEW LIAM, WAS COMMENTED IN BEFORE
//	{
		//if ( RUnif(0,1) < sd_.getMortCarbon() )
//		if ( RUnif() < sd_.getMortCarbon() )		// NEW LIAM, WAS COMMENTED IN BEFORE
//		{
//			setDead();					// NEW LIAM, WAS COMMENTED IN BEFORE
// 			cout << "MORTALITY_WATER_AVAIL" << "  height  " << height_ << "  age  " << age_ << "  veg_type  " << sd_.getVegType() << "  evergr  " << MyRound(sd_.getEvergreen()) << " ccum  "<< c_gain_cum_ <<endl;
// 			cout << "MORTALITY_NEG_C" << "  day  " << day_ <<"  height  " << height_ << "  age  " << age_ << "  veg_type  " << sd_.getVegType() << "  evergr  " << MyRound(sd_.getEvergreen()) << " ccum  "<< c_gain_cum_ <<endl;
//			return 1;		// NEW LIAM, WAS COMMENTED IN BEFORE
//		}
//		day_= 0;  // NEW LIAM, WAS COMMENTED IN BEFORE
//	}

	return 0;
}

//-----------------------------------------------------------------------------------


void clPlant::setDead()
{
	alive_  = 0;
	active_ = 0;

	// Need to keep standing and lying leaves (rather than directly adding it to soilcarbon_nwl_)
        // because standing leaves influence light competition within grass individuals (self-shading)
        // and standing+lying leaves influence fire intensity. NWL also includes fine roots, storage, ...,
        // that is material that does not contribute to fire.

	b_leaf_d_st_    += b_leaf_;                                           // SIMON leaves added to dead standing/hanging leaves
	soilcarbon_nwl_ += (0.2*b_root_               + b_repr_ + b_stor_ );  // SIMON fine roots (20% of root biomass)
	soilcarbon_fwl_ += (0.8*b_root_ + 0.2*b_stem_ + b_bark_ );            // SIMON coarse roots (80% of root biomass), branches (20% of stem)
	soilcarbon_cwl_ += (            + 0.8*b_stem_);                       // SIMON coars stem (80% of stem)

	b_leaf_ = 0.;
	b_stem_ = 0.;
	b_root_ = 0.;
	b_repr_ = 0.;
	b_stor_ = 0.;
	b_bark_ = 0.;
	age_	= 0.;
	height_	= 0.;
	
	//cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx PlantClass    " << setw(14) << soilcarbon_nwl_ << setw(14) << soilcarbon_fwl_ << setw(14) << soilcarbon_cwl_ << endl;

	return;
}


//-----------------------------------------------------------------------------------

int clPlant::runMortalityAnnual()
{
//	if ( c_gain_cum_flag_==1 && active_==1 )		// BEFORE LIAM CHANGE
	if ( c_gain_cum_flag_==1 ) 	// NEW LIAM; plant dies regardless of its activity status
	{	double C_mort = RUnif(0,1) ;	
		if (  C_mort < sd_.getMortCarbon() )
		{
			setDead();
 			//cout << "MORTALITY_CARBON" << "  height " << setw(10) << height_ << ", age " << setw(4) << age_ << ", veg_type " << setw(2) << sd_.getVegType() << endl;
			return 1;
		}
		c_gain_cum_flag_= 0;
	}
	
	if ( RUnif(0,1) < sd_.getMortBackground() )          // m4, background
	{
		setDead();
 		//cout << "MORTALITY_BACKGR" << "  height " << setw(10) << height_ << ", age " << setw(4) << age_ << ", veg_type " << setw(2) << sd_.getVegType() << endl;
		return 1;
	}
	
//	if ( age_ > 10 && height_<sd_.getMortBiomass() ) // m5, low height
//	{
//		setDead();
// 		//cout << "MORTALITY_GROWTH" << "  height " << setw(10) << height_ << ", age " << setw(4) << age_ << ", veg_type " << setw(2) << sd_.getVegType() << endl;
//		return 1;
//	}

	
	if ( height_<1e-7 ) // m5, low height
	{
//		cout << "MORTALITY_GROWTH" << "  height " << setw(10) << height_ << ", age " << setw(4) << age_ << ", veg_type " << setw(2) << sd_.getVegType() << endl;
		setDead(); //MHE-MIRJAM #added setDead();
		return 1;
	}
		
		
		
	// Critical vegetation height, Niklas and Spatz (2010) Am J Botany
	// Use only internal stem diameter, this assumes that bark does not contribute to stability
	double Hcrit = 3.;
//	if ( sd_.getVegType()==TREE ) Hcrit  = 0.79 * pow( (11.852*sd_.getWoodDensity()+37.)/9.81*sd_.getWoodDensity(), 0.33333 ) * pow( stem_diam_int_, 0.66666 );		// BEFORE LIAM CHANGE

/*	#ifdef W_BRANCHES
		// mechanical stability based on stem diameter, wood density, and the ratio between crown biomass and stem biomass;
		
		if ( sd_.getVegType()==TREE )
		{
			double King_C = 60.511 * (b_leaf_ + b_branch_)/b_stem_ + 5.794	;	// based on King 1981, Oecologia, eq. 5, approximated by simplified linear regression, which is hardly different from original
			double Young_E = 13.321*wood_density_ + 3061.317; 		// Young's modulus of elasticity; linear regression against wood density based on ca. 225 samples, chapter 5 Wood Handbook (USDA, 2010)
			Hcrit = Mypow(0.5 * stem_diam_int_ / sqrt(wood_density_ / (King_C * Young_E)), 2./3.) ;  	// this based on O'Brien et al., 1995, Ecology, eq. 1, assuming bark does not contribute to mechanical stability
		}
	#else*/
		if ( sd_.getVegType()==TREE ) Hcrit  = 0.79 * pow( (11.852*wood_density_+37.)/9.81*wood_density_, 0.33333 ) * pow( stem_diam_int_, 0.66666 );			// NEW LIAM
		// 	if ( RUnif(0,1) < pow( height_, sd_.getMortMechanic2() )/( pow( height_, sd_.getMortMechanic2() )+ pow( 2.*Hcrit, sd_.getMortMechanic2() )))
		// 	if ( RUnif(0,1) < 0.1*pow( height_, sd_.getMortMechanic2() )/( pow( height_, sd_.getMortMechanic2() )+ pow( 2.*Hcrit, sd_.getMortMechanic2() )))
	//#endif
		
	if ( RUnif(0,1) < sd_.getMortMechanic1()*(height_/Hcrit-1.) )
	{
		setDead();
		// cout << "MORTALITY_INSTAB" << "  height " << setw(10) << height_ << ", age " << setw(4) << age_ << ", veg_type " << setw(2) << sd_.getVegType() << endl;
		return 1;
	}
	
	return 0;
}

//-----------------------------------------------------------------------------------

double clPlant::runGrazing(double balance_, double demand_, double meristem_, int nind_grass_, double weight_coeff_)
{
//	cout << "Demand in PlantClass: " << demand_ << " kg" << endl;
//	cout << "Weight_coeff: " << weight_coeff_ << setw(14) << b_leaf_ << setw(14) << (b_leaf_ - meristem_) * weight_coeff_ << endl;
	
	if(balance_ <= 0.) 	// there is not enough net grass biomass to satisfy the demand of the animals; all available grass gets cropped
	{
		ind_consumption_ = MyMax((b_leaf_ - meristem_), 0.) ; // MyMax to make sure ind only gets cropped when leaf biomass exceeds meristem biomass	
		if(b_leaf_ > 0.) loss_frac_ = ind_consumption_ / b_leaf_ ;	// determine what fraction of the original leaf biomass got consumed; if-statement to avoid potential div-0
		if(b_leaf_ == 0.) loss_frac_ = 0.;  // no leaf biomass means nothing gets consumed anyway
		
		// update plant biomass
		b_leaf_ = MyMin(meristem_, b_leaf_) ;	// all grass inds get cropped down to the meristem; or stay below limit if they already were below due to other reasons
	}
	
	if(balance_ > 0. && b_leaf_ > meristem_)	// there is more grass biomass available than the animals need to satisfy their demand
	{		
//		ind_consumption_ = cropInd(demand_,meristem_,nind_grass_) ;		// determine the amount of biomass removed from a specific individual, based on its SLA (leaning on Lloyd et al., 2010)
		ind_consumption_ = (b_leaf_ - meristem_) * weight_coeff_ + ((b_leaf_ -meristem_) -(b_leaf_ -meristem_) * weight_coeff_)*RUnif(0.,0.1); 		// FLAG: quick surrogate for Lloyd et al. approach
		if(b_leaf_ > 0.) loss_frac_ = ind_consumption_ / b_leaf_ ;		// determine what fraction of the original leaf biomass got consumed; if-statement to avoid potential div-0	
		if(b_leaf_ == 0.) loss_frac_ = 0.; 							    // no leaf biomass means nothing gets consumed anyway	
	
		//update plant biomass
		if(b_leaf_ > 0.) b_leaf_ = b_leaf_ - ind_consumption_ ;		
	}
	
	// re-allocate carbon from storage to leaf biomass (not sure if one should do this here, and about temporal lag between re-allocation and grazing. Think this shouldn't be instantaneous the way it is now, need to
	// come up with a better solution! Suspect this might also be responsible for bm_weight_coeffs > 1 when not doing the recalculation of bm_weight before calling calWeightCoeffs the second time in PlantPopClass?)
	// Liam had the idea to tie reallocation to water_avail, since during drought it does not make sense for a plant to invest into a lot of leaves anyway	
	
	pot_realloc_ = sd_.getStorToLeaf() * loss_frac_ * b_stor_ ;			// non-water limited reallocation from storage to leaves
	b_leaf_+= pot_realloc_ ;
	b_stor_-= pot_realloc_ ;
		
	// update plant architecture
		//cout << "Height 3: " << height_ << endl;
		calHeight();	//cout << "Height 4: " << height_ << endl;
		calCrownArea();
		calLai();
		
	if(ind_consumption_ < 0.) cout << "Something wrong with calculation of ind_consumption!!!" << endl;	
			
	return ind_consumption_ ;
}

//-----------------------------------------------------------------------------------

double clPlant::cropInd(double demand_,double meristem_, int nind_grass_)		// FLAG: This seems to be crap!!! What gets removed is way too low to get anywhere within a decent number of loop iterations in PlantPopClass!
{	
	// A and B regression coefficients derived from  Lloyd et al. (2010) regression, Fig. 2a, scaling their regression to our min and max SLA-range
	const double A = 0.0030644961 ;
	const double B = 0.0083199784 ; 
	double cropped_biomass_ ;
		
	cropped_biomass_ = MyMin(MyMax(0., demand_ * (A * sqrt(sla_) - B + 1./ ((double) nind_grass_))), (MyMax((b_leaf_ - meristem_),0.))) ;	// function 1 in Lloyd et al. (2010)
	// first MyMin to make sure that the calculated cropped biomas based on the Lloyd et al function does not exceed the biomass that is actually available from that individual
	// first MyMax to make sure that consumption does not get negative (as the Lloyd et al function could become negative
	// second MyMax to make sure that biomass only gets removed from that individual if there is more than the meristem on the plant
	
	cropped_biomass_ = MyMax(cropped_biomass_,0.) ;		// FLAG: write-statement in PlantPopClass revealed this to be < 0. sometimes albeit the above Min/Max checks ?!?
		
	return cropped_biomass_ ;
}

//-----------------------------------------------------------------------------------

double clPlant::calWeightCoeff(double bm_slope_, double bm_intercept_, double a_, double b_, double c_)
{
	double bm_weight_ = bm_slope_ * b_leaf_ + bm_intercept_ ;
	double sla_weight_ = a_ * Mypow(sla_,2.) + b_ * sla_ + c_ ;
	
	if(bm_weight_ < -0.0001 || bm_weight_ > 1.0001) 	// that's more than mathemactical slop can be made responsible for!
	{
		cout << "bm_weight_ out of range! " << bm_weight_ << setw(14) << b_leaf_ << setw(14) << bm_slope_ << setw(14) << bm_intercept_ << endl;		// technically this should never happen since line fit between (bm_min,0) and (bm_max,1)
//		exit(1); 	
	}
	
	bm_weight_ = MyMax(MyMin(bm_weight_,1.),0.) ; 	
	
	if(sla_weight_ < -0.0001 || sla_weight_ > 1.0001) // that's more than mathemactical slop can be made responsible for!
	{
		cout << "sla_weight_ out of range! " << sla_weight_ << setw(14) << sla_ << a_ << setw(14) << b_ << setw(14) << c_ << endl;	// technically this should never happen since parabola fit between (sla_min,0) and (sla_max,1)
//		exit(1);
	}
	
	sla_weight_ = MyMax(MyMin(sla_weight_,1.),0.) ;
	
	weight_ = (bm_weight_ + sla_weight_) / 2.;			// for now, give both factors equal weight
	
	return weight_ ;
}

//-----------------------------------------------------------------------------------

void clPlant::runFireModel( double intensity, double flame_height )
{
/*
 * NEW CAMILLE
Fire doesn't kill the tree directly at any stage, it would require to kill the rooting system, and thus be intense and long enough, but we do not simulate that.
As such, death due to fire happens only "mechanically" when plant tries to resprout.
* This needs to be checked, to see how it happens actually.
*/
	// remove dead biomass
	b_leaf_d_st_ = 0. ;
	b_leaf_d_ly_ = 0. ;
	
	// topkill only for live individuals
	if ( alive_==1 )
	{
						
		if ( sd_.getVegType()==TREE )
		{
		//1 - Calc how much bark is burned; cf ref adgvm1: Williams 2003, Russell Smith 2009
			// moved to PlantPopClass.cpp : flame_height=20.1*(1+ exp(-((0.41*intensity)/1000))) ;
			float BurnC = (0.884/(1+exp(2.7184-(1.0601*(log(flame_height))))))-0.00884 ;
			b_bark_ = b_bark_*(1-BurnC) ;
		//1' - Calculate new bark thickness
			stem_diam_tot_ = 2.*Mypow( (b_stem_+2.*b_bark_ ) / ( M_PI*height_*wood_density_ ), 0.5 ) ; // unit : meters
			bark_thick_ = (stem_diam_tot_ - stem_diam_int_)*0.5 ;
// NEW CAMILLE - Note: maybe should call bark_thick calculation equation and add some of those things in it instead?
			// can it be updated in a more efficient way?
		//2 - Probability of survival depending on bark thickness remaining; based on Hoffmann 2012
   
			if  ( bark_thick_ <= 0.02 && RUnif()*0.02 >= bark_thick_ ) // ==> Set Dead due to fire Intensity kill 
			// 20 mm is "TWEAK" threshold from Hoffmann 2012, but in model, diameter, and thus bark_thick_, are in m, so 20mm->0.02m
			// Or should it be > instead of >= ? 
			{
				if (active_==0)   // inactive  - NEW CAMILLE : leaves will resprout next season, but biomass is still there ==> need to try with resprouting some leaves here now too
				{
				
				b_wood_  = 0.05*b_wood_ + (sd_.getStorToWood()*b_stor_*0.65); 			
				b_stem_  = 0.05*b_stem_ + (sd_.getStorToStem()*b_stor_*0.65);				
				#ifdef W_BRANCHES
				b_branch_ = 0.05*b_branch_ + (sd_.getStorToBranch()*b_stor_*0.65); 	
				#endif	
				b_stor_ -=                 sd_.getStorToWood()*b_stor_; // QUESTION : What about the reallocation to stem and branches? check Mirjam
				
				b_leaf_  = 0;
				
				b_bark_  = 0.05*b_bark_ + (sd_.getAllocBark()*b_stor_*0.65);
				b_stor_ -=                 sd_.getStorToStem()*b_stor_;

				calCrownBHeight();
				}


				else // active - resprout now . Was : evergreen 
				{
				// NEW CAMILLE : Order of calculation changed, used to be leaf first, then bark then stem. Think it makes more sense this way. Is it fine?
							
				b_wood_  = 0.05*b_wood_ + (sd_.getStorToWood()*b_stor_*0.65); 			
				b_stem_  = 0.05*b_stem_ + (sd_.getStorToStem()*b_stor_*0.65);				
				#ifdef W_BRANCHES
				b_branch_ = 0.05*b_branch_ + (sd_.getStorToBranch()*b_stor_*0.65); 	
				#endif	
				b_stor_ -=                 sd_.getStorToWood()*b_stor_; // QUESTION : What about the reallocation to stem and branches? check Mirjam
				
				b_leaf_  = 0.05*b_leaf_ + (sd_.getStorToLeaf()*b_stor_*0.65);
				b_stor_ -=                 sd_.getStorToLeaf()*b_stor_;
				
				b_bark_  = 0.05*b_bark_ + (sd_.getAllocBark()*b_stor_*0.65);
				b_stor_ -=                 sd_.getStorToStem()*b_stor_;

				calCrownBHeight();
				}
			}
				
				
			else //Fire intensity doesn't kill
			{
				if ( height_<flame_height ) // Set Dead due to Flame Height kill 
				{
					if (active_==0)   // inactive  - NEW CAMILLE : leaves will resprout next season, but biomass is still there ==> need to try with resprouting some leaves here now too
					{	
					b_wood_  = 0.05*b_wood_ + (sd_.getStorToWood()*b_stor_*0.65); 			
					b_stem_  = 0.05*b_stem_ + (sd_.getStorToStem()*b_stor_*0.65);				
					#ifdef W_BRANCHES
					b_branch_ = 0.05*b_branch_ + (sd_.getStorToBranch()*b_stor_*0.65); 	
					#endif	
					b_stor_ -=                 sd_.getStorToWood()*b_stor_; // QUESTION : What about the reallocation to stem and branches? check Mirjam
				
					b_leaf_  = 0;
					
					b_bark_  = 0.05*b_bark_ + (sd_.getAllocBark()*b_stor_*0.65);
					b_stor_ -=                 sd_.getStorToStem()*b_stor_;
					
					calCrownBHeight();
					}
					else // active - resprout now . Was : evergreen 
					{
					// NEW CAMILLE : Order of calculation changed, used to be leaf first, then bark then stem. Think it makes more sense this way. Is it fine?
					b_wood_  = 0.05*b_wood_ + (sd_.getStorToWood()*b_stor_*0.65); 			
					b_stem_  = 0.05*b_stem_ + (sd_.getStorToStem()*b_stor_*0.65);				
					#ifdef W_BRANCHES
					b_branch_ = 0.05*b_branch_ + (sd_.getStorToBranch()*b_stor_*0.65); 	
					#endif	
					b_stor_ -=                 sd_.getStorToWood()*b_stor_; // QUESTION : What about the reallocation to stem and branches? check Mirjam
				
					b_leaf_  = 0.05*b_leaf_ + (sd_.getStorToLeaf()*b_stor_*0.65);
					b_stor_ -=                 sd_.getStorToLeaf()*b_stor_;
				
					b_bark_  = 0.05*b_bark_ + (sd_.getAllocBark()*b_stor_*0.65);
					b_stor_ -=                 sd_.getStorToStem()*b_stor_;
					

					calCrownBHeight();
					}
				}
				else if (flame_height>=crown_h_base_ && flame_height<height_ )	// height_>flame_height ==> plant survives	
				// need to estimate how much of the crown gets consumed in the fire; update leaf biomass (and branch biomass?), as well as crown_h_base_?
				{
					if (active_==0)   // inactive  - resprout next season. Was : ( sd_.getEvergreen() <= EVERGREEN_THR ) : only deciduous vegetation switches 
        			{
					b_leaf_ = 0; // NEW Camille - consider that when inactive leaves are dry and thus burn (if they didn't already fell), but without killing stem or branches
					crown_h_base_ = flame_height ; // NEW Camille - but flames coming from fuel are still damaging the lower branches
					}
					else // active - resprout now . Was : evergreen 
					{					
        			#ifdef W_BRANCHES
        			
        				// crown lenght changes due to flames consuming the lower part of the crown => crown radius and crown area also changed, leaf and branch biomass is reduced
						// re-set crown length

						crown_l_ = height_ - flame_height ;

						// calculate the remnant crown volume and consumed crown volume

						double remn_crown_vol = M_PI * par_d_ * par_d_ * Mypow(crown_l_, 2./sd_.getCanfrmPar2()+1.)/(2./sd_.getCanfrmPar2()+1.) ;

						// calculate the remaining leaf and branch biomass, assuming it is proportional to the fraction of crown volume lost in the fire

						b_leaf_ *= remn_crown_vol / crown_volume_ ;
						b_branch_ *= remn_crown_vol / crown_volume_ ;

						// recalculate the radius at the new canopy base; par_d_ remains constant, as the outward shape of the crowd won't change by shortening it from the bottom

						crown_r_ = par_d_ * Mypow(crown_l_, 1./sd_.getCanfrmPar2()) ;
						crown_area_0_ = stem_count_ * crown_r_ * crown_r_ * M_PI ; 

						// update crown volume and crown base height

						crown_volume_ = remn_crown_vol ;
						crown_h_base_ = flame_height ;
        				        				      		
        			#else
        				double Re = ((height_+flame_height)/2);			// height for the remnant part
						double Bu = ((flame_height+crown_h_base_)/2);	// height for the burned part
        		
            			double tmp_ca_h_base = Mypow( sd_.getCanfrmPar1()*height_*( 1.-Mypow( flame_height/height_, sd_.getCanfrmPar2() ) ), 2 )*M_PI;   // SIMON 
    
						double tmp_ca_h_base_Re = Mypow( sd_.getCanfrmPar1()*height_*( 1.-Mypow( Re/height_, sd_.getCanfrmPar2() ) ), 2 )*M_PI;   // area at height for the remnant part
						double tmp_ca_h_base_Bu = Mypow( sd_.getCanfrmPar1()*height_*( 1.-Mypow( Bu/height_, sd_.getCanfrmPar2() ) ), 2 )*M_PI;   // area at height for the burned part

						tmp_ca_h_base_Re *= (height_-flame_height); // volume remnant
           				tmp_ca_h_base_Bu *= (flame_height-crown_h_base_); // volume burned
 
            			//double BCF = tmp_ca_h_base_Bu/(tmp_ca_h_base_Re+tmp_ca_h_base_Bu);
 
            			double RCF = tmp_ca_h_base_Re/(tmp_ca_h_base_Re+tmp_ca_h_base_Bu); // remnant crown fraction

						b_leaf_  = RCF*b_leaf_;    // simply remove a proportion of leaf biomass
						crown_h_base_ = flame_height ;
					#endif
					}
				}
// (canopy_h_base_>flame_height) -> nothing else happens
			
			}    //########### can be deleted?
		} // if ( sd_.getVegType()==TREE )
						
		else // grasses
		{
			// always remove leavses
			b_leaf_  = (sd_.getStorToLeaf()*b_stor_*0.65); // used to be : b_leaf_  = 0.05*b_leaf_ + (sd_.getStorToLeaf()*b_stor_*0.65);
			b_stor_ -=                 sd_.getStorToLeaf()*b_stor_;
		}



	calHeight();     // QUESTION : Mirjam does it, but I used not to. Should it be?
	//calCanopyBHeight();                      ############
	//calStemDiam();   // update stem diameter ############
	calStemDiam();   // update stem diameter
	calCrownArea();
	calLai();
	//age_ = 0;   // commented this out

	//if ( sd_.getEvergreen() <= EVERGREEN_THR ) active_ = 0;  // QUESTION : What's the point? - CAMILLE
	//if ( sd_.getEvergreen() >  EVERGREEN_THR ) evg_flush_flag_=0;  // QUESTION : What's the point? - CAMILLE
// cout << setw(14) << crown_area_0_ << setw(4) << sd_.getVegType() << endl;
    
   	}  // end if ( alive_==1 ) statement

	return;
}



//-----------------------------------------------------------------------------------


void clPlant::calWaterAvailability(GridCellSoilData* soil, int year )
{
	float beta_tmp = 0;
	float gtheta   = 0;
	float root_fracs_tmp[N_SOIL_LAYERS];
	
	float stem_temp = 1.;
	float root_temp = 1.;
	// Calculate rooting depth, assume that roots form a cylinder with radius MIN_ROOT_RADIUS.
	// This depth is only relevent when it is less than the maximum rooting depth.
//	if(sd_.getVegType()==GR_4) wood_density_ = 200.;	// FLAG: MP: This is how it used to be for grasses in the old version; breaks connection between rooting depth and SLA
	
//	float rdepth   = b_root_/( BETA_HELPER*sd_.getWoodDensity() );					// BEFORE LIAM CHANGE
	float rdepth   = b_root_/( BETA_HELPER*wood_density_ ); // p50 method	// NEW LIAM FLAG MP: Check that there is no division by zero here for grasses (wood density proxy for grasses?)...
	
	//if(sd_.getVegType()==GR_4 && b_root_ > 10. && year==2000) cout << "wood density  GR_4: " << wood_density_ <<"  "<< rdepth << "   " << b_root_ << "   " << b_leaf_ << "   " << b_leaf_/b_root_<< endl; 
	
	// When the rooting depth is less than the maximum rooting depth, then the roots form cylinder
	// with constant radius MIN_ROOT_RADIUS. Water upake in all soil layers where the plant has
	// roots is the same. When the rooting depth reaches the maximum rooting depth, then water
	// uptake is calculated by using the fractions defined by the root form.
//-----------------------------------------------------------
	// root fracs original
	if ( rdepth<sd_.getRotfrmMaxD() )
	{
		int l=0;
		while( DEPTH[l]<rdepth ) l++;
		
		for ( int k=0; k<(l+1); k++ ) 
		{
		  root_fracs_tmp[k]=1./(float)(l+1);
// 		  cout << "  root_fraction_top1_" << k << "  " << root_fracs_tmp[k] << endl; 
		}
		for ( int k=(l+1); k<N_SOIL_LAYERS; k++ ) 
		{
		  root_fracs_tmp[k]=0;
// 		  cout << "  root_fraction_top2_" << k << "  " << root_fracs_tmp[k] << endl; 
		}  
	}
	else
	{
		for ( int k=0; k<N_SOIL_LAYERS; k++ ) 
		{
		  root_fracs_tmp[k]=root_fracs_[k];
// 		  cout << "  root_fraction_bottom_" << k << "  " << root_fracs_tmp[k] << endl; 
		}
	}

// 		for ( int k=0; k<N_SOIL_LAYERS; k++ ) 
// 		{
// 		  root_fracs_tmp[k]=root_fracs_[k];
// // 		  cout << "  root_fraction_bottom_" << k << "  " << root_fracs_tmp[k] << endl; 
// 		}
// end root fracs original


// SOIL TEXTURE--------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// NOTE: redoing water availability based the cohesion-tension theory 
//---------------------------------------------------------------------------------------------
	double effective_sat[N_SOIL_LAYERS]; // Saturation of soil layer_i for Matric potential method		// NEW LIAM
	double soil_mat_pot[N_SOIL_LAYERS]; // Turn saturation into matric potential - MPa					// NEW LIAM

	//---------------------------------------------------------------------------------------------
	// Cohesion tension parameters		// NEW LIAM DOWN TO *********
	double p50 = sd_.getP50(); //50% loss of conductivity -- make trait later --- now trait --seedclass constants
	double a_shape = 5; // shape parameter -- make trait later? 
// 	double ks_max = 15; // grass maybe//#kg m –1 s –1 MPa –1 from McDowell Tree Physiology 22, 763–774 NOT SURE THAT THIS IS CORRECT
	double pw = 999.97; //#g/m^3 specific gravity of water I think
	double grav = 9.8; // m/s^2
	double nu = 1.002; // #ast 20degC pascal per sec = 10^-3 kg/m/s NOTE: It's really easy to make this a fn. of temperature 
// 	double p50_root = MyMin(p50+0.2+((height_*pw*grav)/1000000), -0.0001);		// FLAG: using this one now instead of the one below, 26.03.2015
	double p50_root = MyMin(p50+((height_*pw*grav)/1000000), -0.0001);	// Liam: 171115

	// CAMILLE, in your version grass sapwood area was always 0 because stem_diam_tot_ is zero, this is why you did not get grasses
	double area_sapwood_ = 0;                                             // SIMON
// 	if(sd_.getVegType()==TREE) area_sapwood_ = stem_count_*1.582*pow(stem_diam_tot_ , 1.764); // SIMON I think multiplication must be done outside the power function
	if(sd_.getVegType()==TREE) area_sapwood_ = stem_count_*1.582*Mypow((stem_diam_tot_ * 100), 1.764)/10000; //m^2 
	//if(sd_.getVegType()==TREE) area_sapwood_ = 1.582*pow(stem_diam_tot_*stem_count_ , 1.764); // SIMON
	// 	else area_sapwood_ = (14.79*Mypow(b_root_, 0.8164))/10000; //m^2 //sapwood no longer necessary for grasses

	double area_leaf_tot = crown_area_0_*lai_; // m^2
//-------------------------------------------------------------------------------------------------------
// leaf conductance //NOTE: watch out for infinities from 0 leaf area
	double k_leaf = ((61.922 + (7.758*p50))*(18*0.001))/1000;		// Markenstein 2011 digitised ((18*0.001) mmol to g)(1000 g to kg) (kg m^-1 s^-1 MPa^-1)
	double k_leaf_normalised = 0.0;
	if(sd_.getVegType()==TREE) k_leaf_normalised = (k_leaf * area_leaf_tot)/0.2;	//kg.s^-1.MPa^-1 // 0.2 assumes 20cm leaf length 
	if(sd_.getVegType()==GR_4) k_leaf_normalised = (k_leaf * area_leaf_tot)/MyMax(1.00, height_);	//kg.s^-1.MPa^-1 // I'm assuming grass leaf length is a minimum of 1m and a max of height(1.5m) 
	double leaf_resistance = 1/k_leaf_normalised;
// -------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// sapwood conductivity //NOTE: watch out for infinities from 0 leaf area
	double k_sw = ((581.85 + (90.07*p50))*18)/1000;		// Markenstein 2011 digitised (18 mol to g)(1000 g to kg) (kg m^-1 s^-1 MPa^-1)
	double k_sw_normalised = ((k_sw * (area_sapwood_))/MyMax(0.05, height_));	//kg.s^-1.MPa^-1
	double sw_resistivity = 1/k_sw_normalised;
//-------------------------------------------------------------------------------------------------------
	double min_plant_wat_pot = (p50*-1); // We assusme this is p50 for now
	double E_canopy[N_SOIL_LAYERS] = {0.};
	double ks_max[N_SOIL_LAYERS] = {0.}; 	
	double ks_soil[N_SOIL_LAYERS] = {0.};
	double R_soil_root[N_SOIL_LAYERS] = {0.};
	
	Gi_weighted_ = 0.;
	soil_mat_pot_weighted_ = 0.; // this will be the value for the previous days water conditions in daily processes
	double len_root = 0.; // root length calculated following Fisher et al. 2007
	double R_soil_root_tot = 0.;
	double ks_root = 0.; 
	
// ********* END NEW LIAM

//---------------------------------------------------------------------------------------------


	for ( int j=0; j<N_SOIL_LAYERS; j++ ) 
	{
				
//----------------------------------------------------------------------------------------------------
// Root length
		if ( sd_.getVegType()==TREE ) len_root = (b_root_*0.2*1000*root_fracs_tmp[j])/(((0.259+(-0.05921*sd_.getP50()))*1000*1.548462*1000)*3.14*Mypow(0.005, 2)); 
		if ( sd_.getVegType()==GR_4 ) len_root = (b_root_*0.9*1000*root_fracs_tmp[j])/(((0.259+(-0.05921*sd_.getP50()))*1000*1.548462*1000)*3.14*Mypow(0.005, 2)); 
//----------------------------------------------------------------------------------------------------		
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------		
// Soil hydraulic conductivity		
		ks_soil[j] = soil[j].ksat*Mypow((soil[j].smo/soil[j].sat_water_cont), ((2*(R[(int)soil[j].soil_texture]))+2));	
//----------------------------------------------------------------------------------------------------		
//----------------------------------------------------------------------------------------------------
// Soil to root resistance 
		if (root_fracs_tmp[j] > 0 && sd_.getVegType()==TREE) R_soil_root[j] = MyMax(0.0, log((sqrt((1/(len_root*3.14))/0.0005))/(2*3.14*(len_root)*(((ks_soil[j]*0.0098)*1000)/3600)*9.81)));
		if (root_fracs_tmp[j] > 0 && sd_.getVegType()==GR_4) R_soil_root[j] = MyMax(0.0, log((sqrt((1/(len_root*3.14))/0.0005))/(2*3.14*(len_root)*(((ks_soil[j]*0.0098)*1000)/3600)*9.81)));
		if (root_fracs_tmp[j] == 0) R_soil_root[j] = 0.0;
		R_soil_root_tot += R_soil_root[j];
// 		cout << "  R_soil_root_" << j << "    " << N_SOIL_LAYERS << "    " << R_soil_root[j] << endl;
//----------------------------------------------------------------------------------------------------
// Root resistance
		if(sd_.getVegType()==0) ks_max[j] = 0.0000050*(b_root_*0.2)*1000*root_fracs_tmp[j];		// Hickler 2006 (10-7 m3 kg-1 s-1 MPa-1) here (kg kg-1 s-1 MPa-1) no. for grasses used for both trees and grasses
		if(sd_.getVegType()==1) ks_max[j] = 0.0000050*(b_root_*0.9)*1000*root_fracs_tmp[j];
		ks_root += ks_max[j];		
//----------------------------------------------------------------------------------------------------
		
//Cohesion tension start------------------------------------------------------------------------------
		// Matric potential method		// NEW LIAM, ENDS AT ***
		// NOTE: Check wilt_point calculation is in correct units	
		if(soil[j].smo != soil[j].smo) {cout << "soil moisture NaN for layer " << j << " in calWaterAvailability at pos 1! " << soil[j].twp << "   " << soil[j].tfc << "   " << soil[j].soc << "   " << soil[j].son << "   " << soil[j].ksat << " " << soil[j].smo << "   " << soil[j].root_bm << "   " << soil[j].sat_water_cont << "   " << soil[j].soil_texture << endl; exit(1);}
		effective_sat[j] = (soil[j].smo-RESID_CONTENT[(int)soil[j].soil_texture])/(soil[j].sat_water_cont-RESID_CONTENT[(int)soil[j].soil_texture]);	// FLAG: CHECK! This may sometimes lead to the production of NaN-values/negative values

		if( effective_sat[j] < 0. || effective_sat[j] != effective_sat[j] ) effective_sat[j] = 0.;
		soil_mat_pot[j] = (-(PHI_S[(int)soil[j].soil_texture]*(Mypow(effective_sat[j], -R[(int)soil[j].soil_texture]))))/100.; //Psi soil
		if(soil_mat_pot[j] != soil_mat_pot[j]) cout << "soil_mat_pot of soil layer " << j << " is NaN!!!!" << endl;
		soil_mat_pot_weighted_+=(soil_mat_pot[j]*root_fracs_tmp[j]);
//----------------------------------------------------------------------------------------------------		
//----------------------------------------------------------------------------------------------------
// Soil water content method
// NOTE: Check wilt_point calculation is in correct units		

		if( soil[j].smo != soil[j].smo)
		{
			if(soil[j].smo != soil[j].smo) cout << "soil moisture NaN for layer " << j << " in calWaterAvailability at pos 2! " << soil[j].twp << " " << soil[j].tfc << " " << soil[j].soc << " " << soil[j].son << " " << soil[j].ksat << " " << soil[j].smo << " " << soil[j].root_bm << " " << soil[j].sat_water_cont << " " << soil[j].soil_texture << endl;
			beta_tmp = 0.0; 
		}
	}

//----------------------------------------------------------------------------------------------------
// Liam 201115
	double delta_psi = MyMin(((min_plant_wat_pot) - ((height_*pw*grav)/1000000)) + soil_mat_pot_weighted_, 1.0); 
	double test_E_canopy = 0.0; 
	if(sd_.getVegType()==TREE && delta_psi > 0) test_E_canopy = ((delta_psi)/(leaf_resistance + sw_resistivity + (R_soil_root_tot) + (1/ks_root)))*86400; //NOTE: Is 86400 reasonalble? 24hrs
	if(sd_.getVegType()==GR_4 && delta_psi > 0) test_E_canopy = ((delta_psi)/(leaf_resistance + (R_soil_root_tot) + (1/ks_root)))*86400; //NOTE: Is 86400 reasonalble? 24hrs
	double Gi_other = (1-(1/(1+Myexp(a_shape*(soil_mat_pot_weighted_-p50_root)))));
	Gi_weighted_ = Gi_other;
	double WA = 0.0;
	if (test_E_canopy > 0 && evapotr_ > 0) WA = MyMin(1, (test_E_canopy/evapotr_)*Gi_other);
	beta_tmp = WA;
//----------------------------------------------------------------------------------------------------
	if (beta_tmp > 1.0)
	{
	  beta_tmp = 1.0;
	}
	if (beta_tmp != beta_tmp) //flag and safety for nans 		// NEW LIAM
	{	
	  beta_tmp =0.0; 		// NEW LIAM
	}
	water_avail_ = beta_tmp;
	
	return;
}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
// Proportion of nwl added to soil model, this avoids that all biomass is
// immediately added to soil carbon model when a plant dies.
float clPlant::getSoilCarbonNWL()
{
	float tmp = 0.01*soilcarbon_nwl_;
	soilcarbon_nwl_ -= tmp;

	return tmp;
}

//-----------------------------------------------------------------------------------
// Proportion of fwl added to soil model, this avoids that all biomass is
// immediately added to soil carbon model when a plant dies.
float clPlant::getSoilCarbonFWL()
{
	float tmp = 0.001*soilcarbon_fwl_;
	soilcarbon_fwl_ -= tmp;

	return tmp;
}

//-----------------------------------------------------------------------------------
// Proportion of cwl added to soil model, this avoids that all biomass is
// immediately added to soil carbon model when a plant dies.
float clPlant::getSoilCarbonCWL()
{
	float tmp = 0.0001*soilcarbon_cwl_;
	soilcarbon_cwl_ -= tmp;

	return tmp;
}
//-----------------------------------------------------------------------------------

void clPlant::setInitialValues(clSeed sd, GridCellSoilData* soil)
{
	sd_ = sd;
	active_          = 1.;		// MP changed to active instead of inactive, as otherwise we shouldn't initialize leaf biomass either!
	c_gain_          = 0.;
	c_gain_cum_      = 0.;
	c_gain_mmolm2s_  = 0.;
	c_gain_cum_flag_ = 0;
	
	if ( sd_.getVegType()==TREE )
	{
		b_leaf_			= sd_.getSeedWeight()*0.07;
		b_wood_			= sd_.getSeedWeight()*0.43;
		#ifdef W_BRANCHES
			b_stem_			= sd_.getSeedWeight()*0.43;	// MP
			b_branch_		= b_wood_-b_stem_;	// MPx	
		#else
			b_stem_			= b_wood_;
		#endif	
		b_root_			= sd_.getSeedWeight()*0.3;
		b_bark_			= sd_.getSeedWeight()*0.05;
		b_repr_			= sd_.getSeedWeight()*0.05;
		b_stor_			= sd_.getSeedWeight()*0.1;
	}
	else
	{
		b_leaf_			= sd_.getSeedWeight()*0.4;
		b_wood_			= 0.;
		#ifdef W_BRANCHES
			b_branch_		= 0.;
			b_stem_			= 0.;
		#else
			b_stem_			= b_wood_;
		#endif		
		b_root_			= sd_.getSeedWeight()*0.4;
		b_bark_			= sd_.getSeedWeight()*0.0;
		b_repr_			= sd_.getSeedWeight()*0.1;
		b_stor_			= sd_.getSeedWeight()*0.1;
	}
	//b_leaf_d_st_ = 0.;  // SIMON not set to zero, we keep dead biomas! This influences light competition in establishment phase.
	//b_leaf_d_ly_ = 0.;  // SIMON not set to zero, we keep dead biomas! This influences light competition in establishment phase.
	
	age_			= 0;
	dpos_			= 0;
	dneg_			= 0;
// 	c_def_counter_	= 0;
	dleaf_alloc_	= 0;
	evg_flush_flag_ = 0;
	day_ = 0;
	
	if ( sd_.getEvergreen() > EVERGREEN_THR ) active_ = 1.;
	
	calWoodDensity() ;						// needed to calculate stem diameter
	#ifdef W_BRANCHES
		calStemCount() ;		// needed in height calucation
		height_ = exp( sd_.getBiomasPar1()/sd_.getBiomasPar2()) * pow((b_stem_)/stem_count_, 1./sd_.getBiomasPar2()) ; //height for starting
		calCrownBHeight();
		calCrownVolume() ;
	#else		
		calStemCount(); 		// NEW CAMILLE								// MIRJAM: moved this up above calHeight(), otherwise stem_count_ won't be initialized in calHeight() !!!
		calHeight();  			// update tree height
		calCrownBHeight();
		calStemDiam(); 			// update stem diameter
	#endif	
	
	calCrownArea();
	calLai();
	stom_cond_ = 0.;
	stom_cond_max_ = 0.;
	boul_cond_ = 0.;
	evapotr_   = 0.;
	evapotr_max_   = 0.;
	leaf_temp_ = 0.;

	calWiltPoint( soil );
	alive_ = 1;
	
	//light_ext_ = 0.3+0.3*sd_.getCanfrmHBase();	// SIMON
	light_ext_ = 0.4+0.3*sd_.getCanfrmHBase();	// SIMON

	calLightAvail();
	
	calLDMC();	// leaf dry matter content fn. of P50				// NEW LIAM
	calSLA();    // 0.1 instead of 0.01 to transform from g/cm^2 to kg/m^2				// NEW LIAM

	calLeafLongevity();  // factor 10 to transform units // using p50
// 	root_turnover_ = 0.00137;  // NEW LIAM; ~70 % grass root turnover per year
	root_turnover_ = 0.00245; // ca. 90% turnover  SIMON
	root_turnover_tree_ = 0.0005;
	
//	cout << crown_l_*100. << setw(14) << sd_.getSeedWeight() << setw(14) << sd_.getCanfrmPar2() << setw(14) << 2./sd_.getCanfrmPar2()+1. << setw(14) << crown_r_*100. << setw(14) << par_d_ << setw(14) << crown_area_0_ << setw(14) << crown_r_/crown_l_ 
//	<< setw(14) << b_stem_ << setw(14) << (b_branch_+b_leaf_)/b_stem_ << setw(14) << rho_crown_ << setw(14) << (b_leaf_+b_branch_)*(2./sd_.getCanfrmPar2()+1.)/(M_PI*0.16*pow(crown_l_,3.)) << endl ;
 	
 	double r_conventional = sd_.getCanfrmPar1()*height_*pow(1.-crown_h_base_/height_ ,sd_.getCanfrmPar2()) ;		// based on old crown area calculation
// 	double r_volume	= sqrt((b_leaf_+b_branch_)*(2./sd_.getCanfrmPar2()+1.)/(M_PI*rho_crown_*(height_-crown_h_base_))) ;		// based on crown volume calculation
	
//	cout << (height_-crown_h_base_) * 100. << setw(14) << r_conventional*100. << setw(14) << r_volume*100. << setw(14) << rho_crown_ << setw(14) << r_conventional/(height_-crown_h_base_) << setw(14) << r_volume/(height_-crown_h_base_) << endl ;

	can_frm_helper_ = pow( 0.5, 1./sd_.getCanfrmPar2() );
	
	soilcarbon_nwl_ = 0.;  // SIMON non woody litter for soil C model
	soilcarbon_fwl_ = 0.;  // SIMON fine woody litter for soil C model
	soilcarbon_cwl_ = 0.;  // SIMON coarses woody litter for soil C model
	
	double sum = 0;
	double C1 = sd_.getRotfrmPar1();
	double C2 = sd_.getRotfrmPar2();
	double MD = sd_.getRotfrmMaxD();
	
	for ( int i=0; i<N_SOIL_LAYERS; i++ )
	{
		if ( (DEPTH[i]-THICK[i]*0.5) < MD )
			root_fracs_[i] =      pow( ( 1. - pow( (DEPTH[i]-THICK[i]*0.5)/MD, C1 )  ), C2 )    *THICK[i];
		else
			root_fracs_[i] = 0;
		
		sum += root_fracs_[i];
		root_bm_[i] = 0;
	}
	if (sum>0)
	{
		for ( int i=0; i<N_SOIL_LAYERS; i++ )
		{
			root_fracs_[i] /= sum;
		}
	}
	
	phen_loop_=0;
	for ( int i=0; i<N_PHEN; i++ ) phen_rain_[i] = 0.;
	for ( int i=0; i<N_PHEN; i++ ) phen_lght_[i] = 0.;

	
	// minimum avoids that respiration exceeds biomass, resp=helper*biomass*temperature_factor where temperature_factor <2.05
	leaf_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaLeafResp()/sd_.getCNRatioLeaf() );  // *MAIN_RESP_CARB_PROPORT(=1)
	wood_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaWoodResp()/sd_.getCNRatioWood() );  // *MAIN_RESP_CARB_PROPORT(=1)
	root_resp_helper_ = MyMin( 0.45, sd_.getBetaNResp()*sd_.getBetaRootResp()/sd_.getCNRatioRoot() );  // *MAIN_RESP_CARB_PROPORT(=1)
	sd_.setCheckAllocParameters();
	
	return;
}

//-----------------------------------------------------------------------------------





















