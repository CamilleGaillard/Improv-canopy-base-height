#include "GridCellClassConstants.h"
#include "LeafGlobals.h"
#include "Leaf.h"

#include "MyMath.h"
#include "GridCellClass.h"
#include "Radiation.h"

#include <iostream>
#include <iomanip>



using namespace::std;


// This constructor does not initialize the wilting points,
// need to call calWiltPoints.

clGridCell::clGridCell()
{
	latcount = 0;
	loncount = 0;
	plantdata.lai  = 0;
	plantdata.gsc  = 0;
	plantdata.cc3  = 0.;
	plantdata.cc4  = 0.;
	plantdata.ctr  = 1.-plantdata.cc3-plantdata.cc4;	// MP FLAG: This does not allow the existence of bare ground!
	plantdata.A0C3 = 0;
	plantdata.A0C4 = 0;
	plantdata.RLC4 = 0;
	plantdata.RLC3 = 0;
	plantdata.etr = 0;
for(int i =0; i<N_SOIL_LAYERS;i++)
	{
	soil[i].twp = 0;
	soil[i].tfc = 0;
	soil[i].root_bm=0;
	}

	soilcarbon_.PrintCarbonPools(); // SIMON soil carbon
		
}

//---------------------------------------------------------------------------------------

void clGridCell::initialize(int runid, int fire, MyInData *IData, unsigned int latcount_in, unsigned int loncount_in)
{
	runid_ = runid;
	fire_  = fire;
	latcount = latcount_in;
	loncount = loncount_in;
	lat_   = IData->climlats(latcount);
	lon_   = IData->climlons(loncount);

	#ifdef REPRODUCIBLE	
		unsigned int gridseed = (lon_ * 18000. + lat_) * 1000 + runid_ * 100 ;
		gen.seed(gridseed) ; 	// to produce a unique reproducible random seed for each grid cell coordinate
	
		cout << "using unique reproducible random seed: " << lon_ << setw(14) << lat_ << setw(14) << runid_ << setw(14) << dec << gridseed << endl; 
		
		double test = 0.; 
		for(int i=0;i<10000000;i++) test +=RUnif();
//		cout << "average of RUnif()="<< test/10000000. <<endl;
	
	#endif
	
	//cout<<lat_<<" " <<lon_<<endl;

	//initialize all local soil and climate datas from the IData structure.
	for ( int j=0; j<N_SOIL_LAYERS; j++ )		// loop over soil layers
			{
				soil[j].smo = IData->theta_fc(j,latcount,loncount);			
				soil[j].ksat = IData->k_sat(j,latcount,loncount);	 //used to be a constant in old model
				soil[j].sat_water_cont = IData->sat_water_cont(j,latcount,loncount);  //used to be a constant in old model
				soil[j].twp = IData->theta_wp(j,latcount,loncount);	 // used to be reset through old IData.getAnnualData
				soil[j].tfc = IData->theta_fc(j,latcount,loncount);	 // used to be reset through old IData.getAnnualData		
				soil[j].soc = IData->soil_C(j,latcount,loncount);	 // used to be reset through old IData.getAnnualData	
				soil[j].son = IData->soil_N(j,latcount,loncount);        // used to be reset through old IData.getAnnualData
				soil[j].soil_texture = IData->soil_texture(j,latcount,loncount);  //used to be a constant in old model
				//soil[j].perc_exponent = IData.perc_exponent(j,latcount,loncount);		
						
				if(soil[j].smo != soil[j].smo) cout << "soil moisture NaN for layer " << j << " in GridCell initialize!" << endl;
				
			}
			
	// initialize soil carbon model with soc from input data, convert from g to kg         // SIMON soil carbon
	cout << endl << "Initialize soil carbon model with soc from data: " << soil[0].soc << endl;
	soilcarbon_.Initialize(soil[0].soc/1000.);   // SIMON soil carbon
	soilcarbon_.PrintCarbonPools();   // SIMON soil carbon			
			
	//initialize climate on a 365 day basis (this used to be in clInData::getAnnualData)	
	int month_in_year;
	double radnet_day;
	double SunHrs;
	
	for (int i=0;i<MONTH_IN_YEAR;i++)
	{
		rd0_[i]=IData->rd0(i,latcount,loncount);
		pwet_[i]=IData->pwet(i,latcount,loncount);
		ralpha_[i]=IData->ralpha(i,latcount,loncount);
		rbeta_[i]=IData->rbeta(i,latcount,loncount);
	}
	RainFallYear();	// "pre" used to be reset through old IData.getAnnualData
	
//	cout << "GridCell initialize, after RainFallYear(): " << RUnif() << endl;
	
	for ( int day=0; day<DAYS_IN_YEAR; day++ )
	{
		month_in_year = (int) floor( ((double) (day+1))/30.42 );
		
		GetNetRadiation( lat_, IData->sunp(month_in_year,latcount,loncount), IData->tmp_max(month_in_year,latcount,loncount),
						 IData->tmp_min(month_in_year,latcount,loncount), IData->eA12(month_in_year,latcount,loncount),
						 day, &radnet_day ,&SunHrs);
		
		clim[day].par = IData->Q012(month_in_year,latcount,loncount)/10.;  
		clim[day].rad = radnet_day; 					   
		clim[day].tmp = IData->tmp_day(month_in_year,latcount,loncount);   
		clim[day].hum = IData->reh(month_in_year,latcount,loncount);       
		clim[day].wnd = IData->wnd(month_in_year,latcount,loncount);	   
		clim[day].apr = IData->atm_press(latcount,loncount);		   
		clim[day].co2 = 387.;  //NOTE CO2 is global parameter here!	   
		clim[day].sun = IData->sunp(month_in_year,latcount,loncount);	   
	}
		
	plantdata.max_root_depth[0] = IData->root_depth(0,latcount,loncount);	// NEW LIAM, FROM adgvm.cpp line 209; this is the max root depth allowed for trees
	plantdata.max_root_depth[1] = MyMin(0.6,IData->root_depth(0,latcount,loncount));	// NEW LIAM, FROM adgvm.cpp line 210; this is the max root depth allowed for grasses; constrain to a maximum depth of 60 cm
//	cout << "global_max_root_depth: "<< plantdata.max_root_depth[0]<<" "<<plantdata.max_root_depth[1]<<endl;
	//initialize plant population
	
//	mypop_.initialize( runid, fire, lat_, lon_,soil, latcount_in, loncount_in);  //FLAG! for now the population only gets the "top" soil layer
//	cout << "before mypop.initialize(): " << RUnif() << endl; 
	mypop_.initialize( this );  //FLAG! for now the population only gets the "top" soil layer

//	cout << " ... aDGVM gridCell ("<<loncount<<", "<<latcount<< ") initialized, coords: "<<lon_<<", "<<lat_ <<" degree, start simulations. " << RUnif() << endl;
}

//-----------------------------------------------------------------------------------------
void clGridCell::initializeyear(MyInData *IData)
{
/*  // some extra tests to see if all the constants are still as they are supposed to be

	//initialize all local soil and climate datas from the IData structure.
	for ( int j=0; j<N_SOIL_LAYERS; j++ )		// loop over soil layers
			{
				//soil[j].smo = IData->theta_wp(j,latcount,loncount);			
				 if (soil[j].ksat != IData->k_sat(j,latcount,loncount)) cout << "Warning ksat unexpectedly changed.."<<endl;			
				 if (soil[j].sat_water_cont != IData->sat_water_cont(j,latcount,loncount)) cout << "Warning sat_water_cont unexpectedly changed.."<<endl;	
				 if (soil[j].twp != IData->theta_wp(j,latcount,loncount)) cout << "Warning twp unexpectedly changed.."<<endl;			 
				 if (soil[j].tfc != IData->theta_fc(j,latcount,loncount)) cout << "Warning tfc unexpectedly changed.."<<endl;			 
				 if (soil[j].soc != IData->soil_C(j,latcount,loncount)) cout << "Warning soc unexpectedly changed.."<<endl;			 
				 if (soil[j].son != IData->soil_N(j,latcount,loncount)) cout << "Warning son unexpectedly changed.."<<endl;			 
				 if (soil[j].soil_texture != IData->soil_texture(j,latcount,loncount)) cout << "Warning soil_texture unexpectedly changed.."<<endl;	
				//soil[j].perc_exponent = IData.perc_exponent(j,latcount,loncount);				

			}
	//initialize climate on a 365 day basis (this used to be in clInData::getAnnualData)	
	int month_in_year;
	double radnet_day;
	double SunHrs;
	
	for (int i=0;i<MONTH_IN_YEAR;i++)
	{
		 if (rd0_[i]!=IData->rd0(i,latcount,loncount)) cout << "Warning rd0_ unexpectedly changed.."<<endl;			
		 if (pwet_[i]!=IData->pwet(i,latcount,loncount)) cout << "Warning pwet_ unexpectedly changed.."<<endl;			
		 if (ralpha_[i]!=IData->ralpha(i,latcount,loncount)) cout << "Warning ralpha_ unexpectedly changed.."<<endl;			
		 if (rbeta_[i]!=IData->rbeta(i,latcount,loncount)) cout << "Warning rbeta_ unexpectedly changed.."<<endl;			
	}
	
	for ( int day=0; day<DAYS_IN_YEAR; day++ )
	{
		month_in_year = (int) floor( ((double) (day+1))/30.42 );
		
		GetNetRadiation( lat_, IData->sunp(month_in_year,latcount,loncount), IData->tmp_max(month_in_year,latcount,loncount),
						 IData->tmp_min(month_in_year,latcount,loncount), IData->eA12(month_in_year,latcount,loncount),
						 day, &radnet_day ,&SunHrs);
		
		if (clim[day].par != IData->Q012(month_in_year,latcount,loncount)/10.) cout << "Warning par unexpectedly changed.."<<endl;	
		if (clim[day].rad != radnet_day) cout << "Warning rad unexpectedly changed.."<<endl;						
		if (clim[day].tmp != IData->tmp_day(month_in_year,latcount,loncount)) cout << "Warning tmp unexpectedly changed.."<<endl;	
		if (clim[day].hum != IData->reh(month_in_year,latcount,loncount)) cout << "Warning hum unexpectedly changed.."<<endl;	
		if (clim[day].wnd != IData->wnd(month_in_year,latcount,loncount)) cout << "Warning wnd unexpectedly changed.."<<endl;	
		if (clim[day].apr != IData->atm_press(latcount,loncount)) cout << "Warning apr unexpectedly changed.."<<endl;		
		if (clim[day].co2 != 387.) cout << "Warning co2 unexpectedly changed.."<<endl;	  //NOTE CO2 is global parameter here!	
		if (clim[day].sun != IData->sunp(month_in_year,latcount,loncount)) cout << "Warning sun unexpectedly changed.."<<endl;	
	}

*/
	RainFallYear();
}

//-----------------------------------------------------------------------------------------

// This procedure generates a vector RYear, which contains the rain/day

void clGridCell::RainFallYear()
{
	int month;
	int lastmonth;
	int d;
	double monthlyrain;
	double eventsize = 0;
	double rain_sum=0;
//	double rain_sum_test = 0;
//	double rain_sum_test_b = 0;
	
	lastmonth = -1;

	// for all days in the year
	for ( d=0; d<DAYS_IN_YEAR; d++ )
	{
		// compute the month
		month = (int) floor( (((double) d)+1)/30.42 );
		
		// If a new month begins, we must compute the total rain for it.
		// Rain is assumed to be gamma-distributed over a month
		// with the parameters ralpha_ and rbeta_ from climate files.
		if ( month != lastmonth )
		{
		//	monthlyrain =tgamma(ralpha_[month])*rbeta_[month]; 	// use the stdc++11 tgamma function which is faster; NOTE: this is the gamma function, but we need the distribution!
			monthlyrain =RGamma(ralpha_[month], rbeta_[month]);	// this generates a gamma-distributed random number
			
			lastmonth = month;
			// if the number of rainy days in the current month is 0, we
			// will have no rain in this month.
			if ( rd0_[month]<=0 )
				eventsize = 0; 
			else
				// else we compute the average amount of rain on a
				// rainy day.
				eventsize = monthlyrain/rd0_[month];
			//rain_sum_test_b += monthlyrain;
		}
		
		if ( RUnif()  <= pwet_[month] )	{clim[d].pre = RExp(eventsize); /* cout << "RUnif() 1: " << RUnif() << endl; */}
		else clim[d].pre = 0;
		rain_sum += clim[d].pre;
		//rain_sum_test +=RGamma(ralpha_[month],rbeta_[month]);
		//cout <<RGamma(ralpha_[month],rbeta_[month])<< " ";
	}
	//cout << rain_sum <<endl;// << " " <<rain_sum_test<< "  " <<rain_sum_test_b <<endl;
//	cout << "End of RainFallYear(): " << RUnif() << endl; 
	return;
}

//-----------------------------------------------------------------------------------------

int clGridCell::runDailyProcesses( int year, int day)
{
	int rval;
	double avg_rain_tmp     = 0.;   // these averages are required for the phenology model
	double avg_rain_min_tmp = 0.;
	double avg_rain_max_tmp = 0.;
	
	double avg_rad_tmp      = 0;   // these averages are required for the phenology model
	double avg_rad_min_tmp  = 0.;
	double avg_rad_max_tmp  = 0.;
	
	double psi     = lat_*M_PI/180.0;
	double dr      = 1.0+0.033*cos(day*2.0*M_PI/365.0);
	double delta   = 0.409*sin(day*2.0*M_PI/365.0-1.39);
	double omega   = acos(-tan(psi)*tan(delta));
	double hrs     = omega*24.0/M_PI;   // daylength defined by latitude (h)
	double sun_dur = clim[day].sun*hrs*3600.;     // sunshine hours, (in s)
	
	double co2_Pa = clim[day].co2*clim[day].apr*1e-6; // convert co2 from ppm into partial pressure (Pa)

// 	NOTE par given in mumol/m^2/s, we don't need these calculations
//STEVE	par = par*5.*1e6/(86400.*sun);  // convert from J/m^2/s to mumol/m^2/s  NOTE need to know sunshine duration!
// 	par = par*5.*1e6/(86400.*sun);  // tapajos data is  mumol/m^2/s 
	
// 	calculate leaf level photosynthetic rate of the grid cell
	GetC3A( clim[day].tmp, soil[0].son, soil[0].soc, co2_Pa, OI_PAR_PREASSURE, clim[day].par, R_MAINT_RESP_C3, ABS_PHOTONS_C3, ALPHA_C3,
			clim[day].hum, clim[day].apr, WindAtHeight(VEGETATION_HEIGHT, clim[day].wnd), CLD_TREE, M_TREE_C3, B_TREE_C3,
			TAU_K25_C3, TAU_Q10_C3, KO_K25_C3, KO_Q10_C3, KC_K25_C3, KC_Q10_C3,
			&plantdata.A0C3, &plantdata.RLC3 );
	
	GetC4A( clim[day].tmp, soil[0].son, soil[0].soc, co2_Pa, KAPPA_C4, clim[day].par, R_MAINT_RESP_C4, ABS_PHOTONS_C4, ALPHAR_F_C4,
			clim[day].hum, clim[day].apr, WindAtHeight(VEGETATION_HEIGHT, clim[day].wnd), CLD_GRASS, M_GRASS_C4, B_GRASS_C4,
			&plantdata.A0C4, &plantdata.RLC4 );	
	
	rval= mypop_.runDailyProcesses(this, year, day, plantdata.A0C3, plantdata.A0C4,  clim[day].wnd, clim[day].pre, clim[day].tmp, clim[day].apr, co2_Pa, 
							  clim[day].hum, sun_dur, clim[day].rad*RUnif( 0.9, 1.1 ) ); //FLAG! for now the population only gets the "top" soil layer
		
	// update parameters
	float aux_[3];
	plantdata.lai = mypop_.getMeanLai(aux_);
	plantdata.gsc = mypop_.getSumStomCond(aux_);
// 	plantdata.gsc = mypop_.getMeanStomCond();   // NOTE do we need mean or sum here? which units? now it is mmol/m^2/s
	plantdata.cc3 = mypop_.getC3GrassFraction();
	plantdata.cc4 = mypop_.getC4GrassFraction();
	plantdata.ctr = mypop_.getTreeFraction();
	plantdata.etr = mypop_.getSumEvapotr(aux_);
	
// check moisture (last) argument, in aDGVM1: (1.-MyTreePop.getDormant(0))*(Rain[day_in_year]-EtSiteRef ) );   // SIMON 

//	double tmp_nwl = mypop_.getSoilCarbonNWL();
//	double tmp_fwl = mypop_.getSoilCarbonFWL();
//	double tmp_cwl = mypop_.getSoilCarbonCWL();
//	cout << "YYY" << setw(14) << tmp_nwl << setw(14) << tmp_fwl << setw(14) << tmp_cwl << endl;
//	soilcarbon_.UpdateCarbonPools( tmp_nwl, tmp_fwl, tmp_cwl, clim[day].tmp, clim[day].pre-aux_[0] );   // SIMON soil carbon


	soilcarbon_.UpdateCarbonPools( mypop_.getSoilCarbonNWL(), mypop_.getSoilCarbonFWL(), mypop_.getSoilCarbonCWL(), clim[day].tmp, clim[day].pre-aux_[0] );   // SIMON soil carbon
//	soilcarbon_.PrintCarbonPools();   // SIMON soil carbon

	return rval;
}

//-----------------------------------------------------------------------------------------

int clGridCell::runAnnualProcesses( int year )
{	
	return mypop_.runAnnualProcesses(year, this); 
}

//-----------------------------------------------------------------------------------------

void clGridCell::finalize()
{
	cout << "   Finalize GridCell" << endl;
	mypop_.finalize();

	cout << "   ... done" << endl;
	
	return;
}

//-----------------------------------------------------------------------------------------

double clGridCell::WindAtHeight( double height, double ref_wind )
{
	double ws = 0.1;
	
	if ( height > 1.3 )
	{
		ws = ref_wind * log( (height-DISPLACEMENT_HEIGHT)/ROUGHNESS_LENGTH ) * WIND_AT_HEIGHT_HELPER;
	}
	
	return ws;
}

//-----------------------------------------------------------------------------------------

