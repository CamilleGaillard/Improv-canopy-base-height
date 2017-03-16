/// Plant population class
#include <iostream>
#include <iomanip>
#include <fstream>

#include "PlantPopClass.h"
#include "SeedBankClass.h"
#include "LeafGlobals.h"
#include "adgvm.h"
#include "GridCellClass.h"



using namespace::std;


//this used to reside in functions.h
// Start constants for fire model, Table 14
// constants for fire intensity (Higgins 2008)
const float FIRE_H                     = 16890.;	// heat yield of fuel (J g-1)
const float FIRE_C                     = 301.;		
const float FIRE_AW                    = 119.7;
const float FIRE_QM                    = 2.6e6;	// energy required to vaporise water in the fuel 
const float FIRE_QV                    = 160749.;  // heat of pre-ignition for the vegetation 


//global variables for netcdf header definition!!
std::map<std::string,std::string> netcdf_longname, netcdf_units;
std::map<std::string,nc_type> netcdf_datatype;

//NOTICE: the varnames for POP and IND should not be identical, especially if their units/nctype/longname differs.

//--------------------------------------------------------------------------------------

#ifdef NC_OUTPUT
	inline float* clPlantPop::data_pointer(std::string varname, std::string longname, std::string units, unsigned int dimlength, nc_type type)
	{//this is for pop_data

		if (buffidx==0)
		{
			pop_data.emplace(make_pair(varname,vector <float>(BUFF_SIZE*dimlength)));//prepare vector
			float* ptr = pop_data[varname].data();
			for (int i=0;i<BUFF_SIZE*dimlength;++i) ptr[i] = FMISSING; //initialize all values to FMISSING!
			netcdf_units[varname]=units;
			netcdf_longname[varname]=longname;
			netcdf_datatype[varname]=type;
		}
		return &((pop_data[varname])[buffidx*dimlength]); 
	}


	inline float* clPlantPop::data_pointer(std::string varname, std::string longname, std::string units, nc_type type)
	{//this is for ind_data
		int dimlength=1;
		if (varname =="time") dimlength = 1;
		else dimlength = MAX_POP_SIZE;
		ind_data.emplace(make_pair(varname,vector <float>(dimlength)));//prepare vector
		netcdf_units[varname]=units;
		netcdf_longname[varname]=longname;
		netcdf_datatype[varname]=type;
		float* ptr = &((ind_data[varname])[0]);
		for (int i=0;i<dimlength;++i) ptr[i] = FMISSING; //initialize all values to FMISSING!
		return ptr;
	}
#endif

//--------------------------------------------------------------------------------------
clPlantPop::clPlantPop()
{
	return;
}
//--------------------------------------------------------------------------------------

float FireIntensity( float fuel,float fuel_moisture, float wind_speed )
{
	return (FIRE_H*1000.*fuel*atan(wind_speed)*FIRE_C*1000.*fuel/(1000.*fuel+FIRE_AW))/(FIRE_QM*fuel_moisture+FIRE_QV*(1.-fuel_moisture));
}
//end of historic functions.h


// This constructor does not initialize the wilting points,
// need to call calWiltPoints.

//--------------------------------------------------------------------------------------

void clPlantPop::finalize()
{
	cout << "     Finalize PlantPop" << endl;

	return;
}
//--------------------------------------------------------------------------------------

void clPlantPop::initialize( clGridCell *gridcell_) //(  int runid, int fire, float lat, float lon,  GridCellSoilData* soil, unsigned int latcount_in, unsigned int loncount_in )
{
//	cout << "Initialize PlantPop " << RUnif() << endl;
	runid_ = gridcell_->getrunid();
	fire_  = gridcell_->getfire();
	lat_   = gridcell_->getlat();
	lon_   = gridcell_->getlon();
	for ( int i=0; i<VTYPE_NUM; i++ ) { seed_bank_[i].initialize(gridcell_);}

	for ( int i=0; i<MAX_POP_SIZE; i++ ) {
		plant_pop_[i].initialize(gridcell_,i);
		plant_pop_[i].calWiltPoint( gridcell_->getsoil() ); // calculate plant wilt points
	}
	
//	npp_annual_ = 0.;
	
	YEARS_TO_RUN = cfg.lookup("simulation.YEARS_TO_RUN") ;
//	cout << "YEARS_TO_RUN: " << YEARS_TO_RUN << endl ;

// 	n_total_=0;
// 	for ( int i=0; i<MAX_POP_SIZE; i++ )
// 	{
// 		if (plant_pop_[i].getAlive()==1 ) 
// 			n_total_ += plant_pop_[i].getNLeaf() + plant_pop_[i].getNWood() + plant_pop_[i].getNRoot()
// 					 +  plant_pop_[i].getNBark() + plant_pop_[i].getNRepr() + plant_pop_[i].getNStor();
// 	}
	
// 	n_total_ *= (SIZE_FACTOR*0.1);   // SIZE_FACTOR->kg/ha, *0.1->g/m^2
// 	cout << n_total_ << endl;
	
// 	n_uptake_annual_ = 0.;		// N uptake, g/m^2/year


	crown_area_0_   = 0.;
	
	firecount_ = 0;		// re-set to zero at beginning of each new year!
	
	for (int i=0; i<N_SOIL_LAYERS; i++ )
	{
		root_bm_[i]= 0.;  
	} 
		
	
	pop_size_[0] = MAX_POP_SIZE;
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{  // not the whole population is alive
		if ( RUnif()<0.8)
//		if (RUnif() < 0.2)
		{
			plant_pop_[i].setDead();
			pop_size_[0]--;
		}
	}
	
	soilcarbon_nwl_ = 0.;  // SIMON non woody litter for soil C model
	soilcarbon_fwl_ = 0.;  // SIMON fine woody litter for soil C model
	soilcarbon_cwl_ = 0.;  // SIMON coarses woody litter for soil C model	
	

#ifdef NC_OUTPUT
	buffidx=0;
	// write to netcdf file
 /*	//this part new to collect output data in the format required to write output to netcdf-files
	// create trait matrix mll_array to collect the trait vector over all individuals up to MAX_POP_SIZE

	size_t DimLengths[2]={MAX_POP_SIZE, TRAIT_VEC_LENGTH}; 
	trait_matrix = new MLL_array(2,DimLengths);
	
	// create the two two population matrix mll_arrays needed for the pop-data, so that the population matrix will hold BUFF_SIZE time steps
	
	size_t PopLen[2]  = {BUFF_SIZE, POP_VEC_LENGTH+1};	
	size_t SoilLen[2] = {BUFF_SIZE, N_SOIL_LAYERS};
	
	pop_matrix = new MLL_array(2,PopLen);													
	soil_moist = new MLL_array(2,SoilLen);*/

	//trait_data= new float[MAX_POP_SIZE*TRAIT_VEC_LENGTH+2];
	//pop_data = new float[BUFF_SIZE*POP_DATA_SIZE+2];
//	trait_data.emplace(make_pair("info",vector <float>(2)));
	//pop_data.emplace(make_pair("info",vector <float>(2)));
	float* p = data_pointer("info","internal info","deg",2,NC_DOUBLE);
	float* p2 = data_pointer("info","internal info","deg",NC_DOUBLE);
//	(trait_data["info"])[0]=gridcell_->getlatc();
//	(trait_data["info"])[1]=gridcell_->getlonc();
	p[0]=gridcell_->getlatc();
	p[1]=gridcell_->getlonc();
	p2[0]=gridcell_->getlatc();
	p2[1]=gridcell_->getlonc();

	

#endif
		
	#ifdef FIXED_TRAITS
	// Take seed from first plant and use it for all plants.
	// Trait combination of this seed is randomly generated.
	clSeed tmp_seed = plant_pop_[0].getSeed();
	for ( int cnt_pop=1; cnt_pop<pop_size_[0]; cnt_pop++ )
	{
		plant_pop_[cnt_pop].setInitialValues(tmp_seed, gridcell_->getsoil() );     // set traits from seed // state variables to initial val.
	}
	#endif
	
	// initializations for grazing 
	
	with_grazing_ = cfg.lookup("grazing.WITH_GRAZING") ;
	if(with_grazing_ == 1)
	{
		meristem_ = cfg.lookup("grazing.MERISTEM") ;			// bottom part of grass ind which is not accessible for animal mouth
		demand_cow_ = cfg.lookup("grazing.DEMAND_COW") ;  		// put this in initialize for now as long as it stays constant
		cows_per_ha_ = cfg.lookup("grazing.COWS_PER_HA") ;		// put this in initialize for now as long as it stays constant
		frequency_ = cfg.lookup("grazing.FREQUENCY") ; 
		demand_ = demand_cow_ * cows_per_ha_ ; 				    // total grass biomass required for all animals per day 
		deficit_ = 0. ;											// grass biomass that lacks to meet the animals' demand
		
		cout << "FREQUENCY: " << frequency_ << " (" << frequency_*365. << " days per year)" << ", daily demand: " << demand_ << " kg, annual demand: " << (frequency_*365.)*demand_ << " kg" << endl; 
	}
		
	//----------------------------------------------------------------------------------
	
	
//	cout << " ... PlantPop initialized." << RUnif() << endl;
	
	return;
}

//--------------------------------------------------------------------------------------

// calculate canopy area of population

float clPlantPop::getCanopyArea()
{
    float canopy_area=0;
    for ( int i=0; i<MAX_POP_SIZE; i++ )
        if (plant_pop_[i].getAlive()==1 ) canopy_area += plant_pop_[i].getCrownArea();
    return canopy_area;
}

//--------------------------------------------------------------------------------------

float clPlantPop::getDryFuel(float *tmp)
{	
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBLeafDeadLy() ;		// no standing dead leaf biomass for trees
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += (plant_pop_[i].getBLeafDeadLy() + plant_pop_[i].getBLeafDeadSt()) * IND_AREA_MULT ;
	}
	
	tmp[0] = tmp[1] + tmp[2] ;
		
	for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.;    	// convert from kg to tons/ha    
       
    return tmp[0]*0.1;		// total dry fuel in units kg/m2
}    

//--------------------------------------------------------------------------------------

float clPlantPop::getWetFuel()		// FLAG: Do this separately for trees and grasses as well, or leave the way it is?
{
	float tmp=0;
	
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		// live grass leaf biomass
		if ( plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 ) tmp += plant_pop_[i].getBLeaf() ;		// in kg
	}
	
	return tmp/PLOT_SIZE ;  // units kg/m2
}

//--------------------------------------------------------------------------------------

void clPlantPop::runGrazing()
{	
	getSumBLeaf(aux_) ; 
	aux_[2] = aux_[2] * 1000. ;			// convert from t/ha to kg/ha !!!
	gross_supply_ = aux_[2] ;							// total grass leaf biomass per ha stand, in kg (getSumBLeaf returns t/ha!!!!!)
	
	net_supply_ = gross_supply_ - (IND_AREA_MULT * meristem_ * pop_size_[2]) ;	// leaf biomass that is actually available for grazers
//	mean_avail_ind_ = net_supply_ / pop_size_[2] ;		// average available leaf biomass per individual
//	min_ind_req_ = ceil(demand_ / mean_avail_ind_) ;			// minimum number of individuals that will statistically be needed to satisfy the demand
	balance_ = net_supply_ - demand_ ;  
//	eat_frac_ = demand_ / net_supply_ ; 	// fraction of net biomass that will be eaten

//	cout << "SumBLeaf before grazing: " << aux_[2] << ", demand: " << demand_ << ", balance: " << balance_ << endl; 	
			
		if(balance_ > 0.) 	// create a preference rank of available grass individuals based on preference weighting coeffecients
		{	
			bm_weight_coeffs() ;	// delivers slope and intercept based on min and max biomass range of grass individuals
			sla_weight_coeffs() ;	// delivers coeffs a, b, and c based on min and max sla range of grass individuals
			
//			cout << "weight coeffs: " << bm_slope_ << ", " << bm_intercept_ << ", " << a_ << ", " << b_ << ", " << c_ << endl;
			
			for(int i=0; i<MAX_POP_SIZE; i++)
			{
				if(plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getAlive()==1 && plant_pop_[i].getBLeaf() > meristem_) 	// only do this for living inds with enough biomass to be grazed; we're not eating dead biomass for now
				{
					weight_coeff_[i] = plant_pop_[i].calWeightCoeff(bm_slope_,bm_intercept_,a_,b_,c_);
//					cout << "aggregated weight coeff: " << i << ", " << weight_coeff_[i] << endl;
				}		
				else
				{
					 weight_coeff_[i] = FMISSING ;
				}	 
			}
			// now draw a random grass individual (random number between 0 and MAX_POP_SIZE, masked by getVegType==GR_4, and for that individual, draw a second random number between 0 and 1
			// compare that second random number against the weight_coeff of that individual. If the second random number is <= weight_coeff, the individual will be grazed (call runGrazing from
			// PlantClass (plant_pop_[i].runGrazing()). Keep track (make a running sum here!) of how much biomass gets removed per individual through runGrazing(). Keep drawing random grass individuals 
			// until the running sum >= demand. Ascertain that no grass individual will get grazed twice/several times per day (? Don't know if this will be necessary, as the reduced biomass of an
			// already-grazed ind will decrease the likelyhood of being grazed again anyway through the recalculation of the weight coefficient). 
			
			double checksum = 0. ;
			int rand_ind ;			// selection of random individual number
			int nind_affected = 0. ;
			int loop_cycles = 0 ;
			do
			{	
				rand_ind = (int) floor(RUnif(0.,(double)MAX_POP_SIZE)) ;									
				if(plant_pop_[rand_ind].getVegType()==GR_4 && plant_pop_[rand_ind].getAlive()==1 && plant_pop_[rand_ind].getBLeaf() > meristem_ && RUnif(0.,1.)<=weight_coeff_[rand_ind]) 
				{	
					double ind_cons = plant_pop_[rand_ind].runGrazing(balance_, (demand_/ IND_AREA_MULT), meristem_, pop_size_[2], weight_coeff_[rand_ind]) ;  // runGrazing needs to return how much biomass got removed, so that the checksum can be updated
//					cout << "IndConsumption: " << ind_cons << setw(14) << ind_cons/( plant_pop_[rand_ind].getBLeaf()-meristem_)*100. << " %" << setw(14) << plant_pop_[rand_ind].getBLeaf()-meristem_ << setw(14) << plant_pop_[rand_ind].getSLA() << endl; 
					checksum+= IND_AREA_MULT * ind_cons ;	
					
					bm_weight_coeffs();		// not sure if I have to do that again here
					sla_weight_coeffs();	// not sure if I have to do that again here
					
					weight_coeff_[rand_ind] = plant_pop_[rand_ind].calWeightCoeff(bm_slope_,bm_intercept_,a_,b_,c_);	// update the grazed ind's weight coeff, to make sure it will be less likely consumed if its number gets drawn 
																														// again during that day's round of grazing (correct?)
					nind_affected++ ;																									
																														
//					cout << "aggregated weight coeff: " << rand_ind << ", " << weight_coeff_[rand_ind] << endl;	
				}
				#ifndef SANITY_CHECKS_OFF  // compile with -DSANITY_CHECKS_OFF to deactivate the following test
    				if (loop_cycles > MAX_POP_SIZE*10000.)
    				{	
    					getSumBLeaf(aux_);
        				cout << "ERROR Too many loop_cycles in run_grazing. Skipping loop" << setw(14) << nind_affected << setw(14) << checksum << setw(14) << aux_[2]*1000.-meristem_*pop_size_[2] << endl;
        				break;  //break if in average every individual was tried to be grazed 10000 times
    				}
    				loop_cycles++;
				#endif	
			}
			while(checksum < demand_);
			
//			cout << "Nind_affected: " << nind_affected << "pop_size_[2]: " << pop_size_[2] << endl;
//			cout << "Supply > Demand, BLeaf removed: " << checksum << setw(14) << nind_affected << setw(14) << demand_ << endl; 
		}
		
		
		else	// demand equals/exceeds net supply, all available grass biomass will be eaten
		{
			for(int i=0; i<MAX_POP_SIZE; i++)
			{
			 	if(plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getAlive()==1) plant_pop_[i].runGrazing(balance_, demand_/ IND_AREA_MULT, meristem_, pop_size_[2], weight_coeff_[i]) ;	
			} 
			
//			cout << "ALL grass inds affected: " << pop_size_[2] << ", Deficit (kg): " << balance_ << endl; 
		}		
	
//	getSumBLeaf(aux_) ;	
//	cout << "SumBLeaf after grazing: " << aux_[2]*1000. << "  SumBLeaf-meristem: " << aux_[2]*1000. - IND_AREA_MULT*meristem_ << endl; 
}

//--------------------------------------------------------------------------------------

// fire potentially gets started if the run has progressed for more than 50 years, and if the soil moisture of the top soil level drops to a level between 60 % field capacity and 100 % field capacity (60 % only if field capacity gets very close to 1 and 
// wilting point gets very close to 0; in all other cases it's higher than 60 % and gets percentually closer to the value of field capacity in cases when field capacity and wilting point are closing in on each other), and at the same time
// precipitation is less than 10 mm, and a drawn random number is less than  

void clPlantPop::runFireModel( int year, int day, float wnd, float hum, float pre, GridCellSoilData* soil )
{
	getSumCanopyArea0(aux_);
//	cout << setw(14) << aux_[0] << setw(14) << aux_[1] << setw(14) << aux_[2] << endl;
	int ignition=0;
	//if (  pre < 5. && RUnif() < (-0.083333*getSumCanopyArea0(aux_)/10000. + 0.2) ) {ignition = 1.;} // cout << setw(4) << day << setw(13) << pre << endl;}
	//if (  pre < 5. && RUnif() < (-0.083333*aux_[1]/10000. + 0.2) ) {ignition = 1.; } // cout << setw(4) << day << setw(13) << pre << endl;}
	
	
	
	
	if (  pre < 5. && RUnif() < ( 0.09*(1. - pow(aux_[1]/10000.,6)/(pow(aux_[1]/10000.,6)+pow(0.25,6))) ) ) {ignition = 1.;}
	
	
//	float tmp_smo = (soil[0].smo+soil[1].smo+soil[2].smo)/3.;  // SIMON
	float tmp_smo = (soil[0].smo) ;		// Simon

	if(soil[0].smo != soil[0].smo) cout << "soil moisture NaN for layer 0 in runFireModel!" << endl; 

	//if ( year>15 && RUnif()< 0.005*(0.5-getSumCanopyArea0(aux_)/10000.) && pre < 5.)
	if ( year>50 && ignition==1 && RUnif()< 0.005 )
	{
		float fuel_biomass  = getDryFuel(aux_)+getWetFuel();
		// float fuel_moisture = (tmp_smo*getDryFuel(aux_)+hum*getWetFuel())/fuel_biomass;
		// float fuel_moisture = (tmp_smo*getDryFuel(aux_)+0.85*getWetFuel())/fuel_biomass;  // SIMON
		//float fuel_moisture = (hum*getDryFuel(aux_)+tmp_smo*getWetFuel())/fuel_biomass; // SIMON

		// first expression for dry fuel is based on Chaney book
		float fuel_moisture = (exp(0.46 + 0.03297*hum*100.)/100*getDryFuel(aux_)+tmp_smo*getWetFuel())/fuel_biomass; // SIMON
		
		float intensity = FireIntensity( fuel_biomass, fuel_moisture, wnd );

		float patchiness   = 0.989-0.595*exp(-0.000837*intensity);
		//float flame_height = 0.061*pow(intensity, 0.667);  //
		float flame_height = 20.1*(1. - exp(-((0.41*intensity)/1000))) ;  // new equation
		cout << "---------" << setw(13) << aux_[1]/10000 << setw(13) << intensity << endl;
		if ( intensity > 300. )
		{
			cout << endl << "FIRE0 " << setw(5) << year << setw(5) << day << setw(14) << getDryFuel(aux_) << setw(14) << getWetFuel() << setw(14) << fuel_biomass << setw(14) << fuel_moisture << setw(14) << intensity << setw(14) << patchiness << setw(14) << flame_height << setw(7) << pop_size_[1] << setw(7) << pop_size_[2] << setw(14) << getSumCanopyArea0(aux_)/10000 << endl;
			
			
			for ( int i=0; i<MAX_POP_SIZE; i++ )
				if ( RUnif()<patchiness ) plant_pop_[i].runFireModel( intensity, flame_height );

			firecount_++ ;
			cout <<         "FIRE1 " << setw(5) << year << setw(5) << day << setw(14) << getDryFuel(aux_) << setw(14) << getWetFuel() << setw(14) << fuel_biomass << setw(14) << fuel_moisture << setw(14) << intensity << setw(14) << patchiness << setw(14) << flame_height << setw(7) << pop_size_[1] << setw(7) << pop_size_[2] << setw(14)<< getSumCanopyArea0(aux_)/10000 << setw(5)  << endl;
		}
	}
	
	return;
}

//--------------------------------------------------------------------------------------

int clPlantPop::runDailyProcesses( clGridCell *gridcell, int year, int day, float A0C3, float A0C4,
						float wnd, float pre, float tmp, float apr,
						float co2, float hum, float sun_dur, float rad )
{
	year_ = year ;
	dct_ = day ;
	GridCellSoilData* soil =gridcell->getsoil();

	int bufferoverflow = 0;

	if (gridcell->getplant()->max_root_depth[0]>=0 && soil[0].twp >=0)
	{ //we are running an actual pixel (NOT ocean pixel)

		if(day==0) firecount_=0	;	// re-set the fire counter at the beginning of each new year

		for (int i=0; i<N_SOIL_LAYERS; i++ )
	 	{
	    	soil[i].root_bm= 0.;  
	    } 
	
		pop_size_[1] = 0. ;
		pop_size_[2] = 0. ;
	  
		for (int i=0; i<MAX_POP_SIZE; i++)
		{
			if(plant_pop_[i].getVegType()==TREE && plant_pop_[i].getAlive()==1) pop_size_[1]++ ;
			if(plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getAlive()==1) pop_size_[2]++ ;
		}  
	
//		cout << "after pop_size calculation: " << pop_size_[0] << "   " << pop_size_[1] << "   " << pop_size_[2] << endl; 
	
		double tmp_root_bm[N_SOIL_LAYERS]; 
// 		float tmp_n_total = n_total_;
	
	
// 		The following calculates helpers required by and constant for all plants. These helpers are passed as arguments
// 		to the plants in the plant population.
	
// 		Helper for temperature effect on maintenance respiration
		float resp_temp_fac = pow( (3.22-0.046*tmp), (tmp-20.)/10. ); 		// MP: this has a maxium at approx. 37 °C and decreases again for higher temperatures
	
		runMortalityDaily();
	
// 		Helpers for evapotranspiration, most variables for Et are constant for all plants, only stom. cond. gs differs
		float s           = 2504.0*exp((17.27*tmp)/(tmp+237.2))/pow(tmp+237.3,2);
		float gama        = apr*0.001*0.000665;
		float rho         = (apr*1e-3)/(0.287*(1.01*(tmp+273.3)));
		float eS          = 0.6108*exp((17.27*tmp)/(tmp+237.3));
		float eA          = hum/100.*eS;
		float vpd         = eS-eA;
	
		float et_helper_0 = (s*rad);  // here rad is in MJ/m^2/day
		float et_helper_1 = (86400.*rho*SP_HEAT*vpd); // 
		float et_helper_2 = (LAMBDA*s);
		float et_helper_3 = (LAMBDA*gama);
	
// 		Helper for transformation mmol/m^2/s to m/s
		float mms_helper  = 1e-6*0.0224*((tmp+273.15)/273.16)*(1.013e5/apr);
	

		// set soilcarbon variables to zero, these variables are re-calculated in the for-loop
		// i.e. a proportion of biomass of plants moves into these pools
		soilcarbon_nwl_ = 0.;  // SIMON non woody litter for soil C model
		soilcarbon_fwl_ = 0.;  // SIMON fine woody litter for soil C model
		soilcarbon_cwl_ = 0.;  // SIMON coarses woody litter for soil C model	
	
		double compheight[NUM_LIGHT_NGB];	// height of competing individuals
		double compheight_lai[NUM_LIGHT_NGB];	// LAI of competing individuals
		double compheight_crown_area[NUM_LIGHT_NGB];	// crown_area_0_ of neighbours	
		int    compinds[NUM_LIGHT_NGB];		// indexes of competing individuals

		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			if (plant_pop_[i].getAlive()==1 )
			{
				compinds[0] = (i-MAX_POP_SIZE_SQRT+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[1] = (i+MAX_POP_SIZE_SQRT+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[2] = (i-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[3] = (i+1+MAX_POP_SIZE)%MAX_POP_SIZE;
			
				compinds[4] = (i-MAX_POP_SIZE_SQRT-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[5] = (i+MAX_POP_SIZE_SQRT-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[6] = (i-MAX_POP_SIZE_SQRT+1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds[7] = (i+MAX_POP_SIZE_SQRT+1+MAX_POP_SIZE)%MAX_POP_SIZE;
			
				for ( int cmp=0; cmp<NUM_LIGHT_NGB; cmp++ )
				{
					// get heigh and lai of competing individuals, only if they are active
					compheight[cmp]     = plant_pop_[compinds[cmp]].getHeight()*(float)plant_pop_[compinds[cmp]].getActive();
					compheight_lai[cmp] = plant_pop_[compinds[cmp]].getLai()   *(float)plant_pop_[compinds[cmp]].getActive();
					compheight_crown_area[cmp] = plant_pop_[compinds[cmp]].getCrownArea()*(double)plant_pop_[compinds[cmp]].getActive();

					// this weights the competitor height to account for the crown form 
					if ( compheight[cmp]>0 ) compheight[cmp] *= plant_pop_[compinds[cmp]].getCanFrmHelper();

					//cout << setw(3) << plant_pop_[compinds[cmp]].getActive() << setw(3) << plant_pop_[compinds[cmp]].getActive() << setw(3) << plant_pop_[compinds[cmp]].getVegType() << setw(12) << compheight[cmp] << setw(12) << compheight_lai[cmp] << endl;
				}
			
			}
		}

// 		Run daily processes of all plants
		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			plant_pop_[i].runDecomposition();   // SIMON this is done for all plants, also dead

			if ( plant_pop_[i].getAlive()==1 )
			{
				plant_pop_[i].runDailyProcesses(year, day, soil, A0C3, A0C4,  wnd, pre, tmp, apr, co2, hum/100., sun_dur, rad,
							resp_temp_fac, et_helper_0, et_helper_1, et_helper_2, et_helper_3, eS, mms_helper, compheight, compheight_lai, compheight_crown_area);
											 
				plant_pop_[i].getRootBmVector( tmp_root_bm ); 
				for (int j=0; j<N_SOIL_LAYERS; j++ )
				{	
					soil[j].root_bm += tmp_root_bm[j];  
				}
			}

			soilcarbon_nwl_ += plant_pop_[i].getSoilCarbonNWL();  // SIMON non woody litter for soil C model
			soilcarbon_fwl_ += plant_pop_[i].getSoilCarbonFWL();  // SIMON fine woody litter for soil C model
			soilcarbon_cwl_ += plant_pop_[i].getSoilCarbonCWL();  // SIMON coarses woody litter for soil C model				
		}

		float tmp_sumBRoot=getSumBRoot(aux_)*1000.;
		for (int j=0; j<N_SOIL_LAYERS; j++ )
		{	
//			cout << "PlantPop: " << soil[j].root_bm <<setw(14)<< tmp_sumBRoot << endl;
		if ( tmp_sumBRoot<=0 ) soil[j].root_bm = 0;
		else                   soil[j].root_bm /=tmp_sumBRoot;       // MP: risk to get NaN here if tmp_sumBRoot==0...
		}
// 		cout << plant_pop_[3].getLeafTemp() << endl;

	
		if ((with_grazing_ == 1) && (RUnif(0., 1.) <= frequency_) && (year > 800) && (year < 1501))  // grazing will occur on that specific day; prioritize grazing over fire, if both happen to occur on the same day
		{	
			runGrazing() ;
//			cout << "GRAZING! " << year << ",  " << day << endl;
		}	
			
		if ( fire_==1 )	runFireModel( year, day, wnd, hum/100., pre, soil );
	
/*		if (day==364)  // last day of the year
		{
//			cout << "No. of fires in year " << year << ": " << firecount_ << endl;		
// 			npp_annual_ += (tmp_biomass_1-tmp_biomass_0)/1000.*0.44*SIZE_FACTOR;  // in tC/ha/year
// 			npp_annual_ += (tmp_biomass_1-tmp_biomass_0)/10.*0.65*SIZE_FACTOR;  // in g C/m^2/year
			npp_annual_ = (getSumCGainCum(aux_))/10.*0.65*SIZE_FACTOR;  // in g C/ha/year
		}*/
	
		//SIMON  moved this to to it in spring of northern or souhern hemisphere
		if ( year>5 && ((gridcell->getlat()>0 && day==30) || (gridcell->getlat()<=0 && day==210)) )
		{
			for ( int i=0; i<MAX_POP_SIZE; i++ )
			{
				if (plant_pop_[i].getAlive()==1 )
				{
					plant_pop_[i].runAnnualProcesses();
				}
			}
			
			runMortalityAnnual();
			runSeedProd();
			seed_bank_[0].runMutation();	// tree seed bank    NOTE mutation now done in crossover
			seed_bank_[0].runCrossover();	// tree seed bank
			seed_bank_[1].runMutation();	// grass seed bank   NOTE mutation now done in crossover
			seed_bank_[1].runCrossover();	// grass seed bank
			runRecruitment(gridcell);
		
		}	
	
// 		n_total_ = 0;
// 		for ( int i=0; i<MAX_POP_SIZE; i++ )
// 		{
// 			if (plant_pop_[i].getAlive()==1 ) 
// 				n_total_ += plant_pop_[i].getNLeaf() + plant_pop_[i].getNWood() + plant_pop_[i].getNRoot()
// 				+  plant_pop_[i].getNBark() + plant_pop_[i].getNRepr() + plant_pop_[i].getNStor();
// 		}
	
// 		n_total_ *= (SIZE_FACTOR*0.1);   // SIZE_FACTOR->kg/ha, *0.1->g/m^2
// 		n_uptake_annual_ += n_total_ - tmp_n_total;		// N uptake, g/m^2/year
	
	}


	#ifdef NC_OUTPUT	

	//-----------------------------------------------------------------------------------

	// create a population vector that holds all POP_VEC_LENGTH population variables
	// NOTE: all monthly output values are not really monthly values, instead they seem to be daily values on the last day of each month, or other chosen
	//            days within a month; 

	//------------------------------


	int m = -1;
	for (int ml=0; ml<MONTH_IN_YEAR; ml++)
		if ((END_OF_MONTH[ml]-1) == day) {m =ml; break;}
	if (m >=0)
	{	
		data_pointer("time","year fraction","a",1,NC_FLOAT)[0]   = (float) year+((float)day)/365.;		// year plus day in decimal presentation; allows any choice of output frequency
		if (gridcell->getplant()->max_root_depth[0]>=0 && soil[0].twp>=0)
		{ //we are running a normal valid pixel. Output real data!

			getSumBLeaf(data_pointer("SumBLeaf","plot leaf biomass","t/ha",3,NC_FLOAT)); 
			getSumBWood(data_pointer("SumBWood","plot wood biomass","t/ha",3,NC_FLOAT));
			#ifdef W_BRANCHES
				getSumBBranch(data_pointer("SumBBranch","plot branch biomass","t/ha",3,NC_FLOAT));
			#endif	
			getSumBStem(data_pointer("SumBStem","plot stem biomass","t/ha",3,NC_FLOAT));
			getSumBRoot(data_pointer("SumBRoot","plot coarse root biomass","t/ha",3,NC_FLOAT));
			getSumBBark(data_pointer("SumBBark","plot bark biomass","t/ha",3,NC_FLOAT));
			getSumBRepr(data_pointer("SumBRepr","plot biomass allocated to reproduction","t/ha",3,NC_FLOAT));
			getSumBStor(data_pointer("SumBStor","plot level biomass in storage compartments","t/ha",3,NC_FLOAT));
			getMeanLai(data_pointer("MeanLai","plot mean leaf area index","m2/m2",3,NC_FLOAT));
			getMeanHeight(data_pointer("MeanHeight","plot mean vegetation height","m",3,NC_FLOAT));
			getSumBasalArea(data_pointer("SumBasalArea","plot total basal area","m2/ha",3,NC_FLOAT));
			getSumCanopyArea0(data_pointer("SumCanopyArea0","canopy area at base","m2/ha",3,NC_FLOAT));
			getMeanWiltPoint(data_pointer("MeanWiltPoint","mean of ind. plant wilting point","volume/volume",3,NC_FLOAT));
			getMeanActive(data_pointer("MeanActive","proportion of plants in active stage","fraction",3,NC_FLOAT));
			getMeanCGain(data_pointer("MeanCGain","mean carbon gain of plot","CHECK",3,NC_FLOAT));
			getMeanRainLight(data_pointer("MeanRainLight","average of phenology control by rain or light","[-]",3,NC_FLOAT));
			getMeanAlive(data_pointer("meanAlive","proportion of ind. alive","fraction",3,NC_FLOAT));
			getSumStomCond(data_pointer("SumStomCond","plot level stomatal conductivity","CHECK",3,NC_FLOAT));
			data_pointer("MeanVegType","fraction of alive individuals that are trees","fraction",1,NC_FLOAT)[0] = getMeanVegType( TREE ) ;	
			//data_pointer("MeanWoodDens", "average wood density of all trees","kg m-3", 1, NC_FLOAT) [0] = getMeanWoodDensity() ;
			getMeanCGain2(data_pointer("MeanCGain2","mean C gain of plot","mmol/plant/s (CHECK)",3,NC_FLOAT));
			getSumEvapotr(data_pointer("PlotEvapotr","avg. plot level evapotr.","mm/(day*m-2)",3,NC_FLOAT));
			getCGain_Plot(data_pointer("CGain_Plot","total C-gain pers stand","mmol/ha/s (CHECK)",3,NC_FLOAT));
			getDryFuel(data_pointer("BLeafDead","plot biomass stored in dead leaves","t/ha",3,NC_FLOAT));	 //getBLeafDead() ;   // SIMON
			getSumCGainCum(data_pointer("SumCGainCum","cumulative total C gain of plot","tC/ha",3,NC_FLOAT));
			getNAliveByType(data_pointer("nind_alive","number of alive individuals","ind/plot",3,NC_INT));
			getMeanAge(data_pointer("meanAge","mean age of alive individuals","years",3,NC_FLOAT));
			getShoot2Root(data_pointer("shoot2root","aboveground to belowground biomass ratio","fraction",3,NC_FLOAT));
			getMeanSla(data_pointer("meanSla","average specific leaf area index","m2/kg",3,NC_FLOAT));
			data_pointer("firecount","No. of fires per year","per year",1,NC_INT)[0] = firecount_ ;		// NEW MIRJAM, 16.03.2015
			getMeancanopy_h_base(data_pointer("Meancanopy_h_base","Canopy height","m",3,NC_FLOAT));
			getMeanStem(data_pointer("MeanStem","Stem number","integer",3,NC_FLOAT)); // NEW CAMILLE

			float* psmo = data_pointer("soilm","soil moisture","fraction",N_SOIL_LAYERS,NC_FLOAT);
			for (int l=0; l<N_SOIL_LAYERS; l++) psmo[l] = soil[l].smo;

		}
		else
		{ //we are running an invalid (ocean) pixel. set all output variables to FMISSING!	
		//for (int k=1;k<POP_DATA_SIZE;k++) pop_data[buffidx*POP_DATA_SIZE+2+k]=FMISSING;
		}
		buffidx++;  //increment buff idx
	}	
		
	if (buffidx>=BUFF_SIZE)	// last year of N_BUFF_YEARS has been reached, time to flush the buffer
	{
	//	cout << "time to flush the output buffer!" << endl;
		buffidx=0; //reset buffer
		// flush the pop_matrix and soil_moist-cube into the overall mll_array from NetcdfOutputClass, where it will be written to the population nc-file,
		// and then reset the two mll_arrays to "missing" (or just overwrite with next block of N-BUFF_YEARS worth of data)
		// This needs to be done externally from the caller routine. We return "bufferoverflow = 1" to indicate this condition
		bufferoverflow=1;

	}			
	//-----------------------------------------------------------------------------------


	#endif
return bufferoverflow;
}

//--------------------------------------------------------------------------------------

int clPlantPop::runAnnualProcesses( int year, clGridCell *gridcell )
{
int bufferoverflow = 0;	
GridCellSoilData* soil =gridcell->getsoil();

if (gridcell->getplant()->max_root_depth[0]>=0 && soil[0].twp >=0)
{
//	npp_annual_ = 0;
// 	n_uptake_annual_ = 0;		// N uptake, g/m^2/year

/* --------------------------------------------
	SIMON moved this to annual process to be able to do this on different days on northern and southern hemisphere
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{
		if (plant_pop_[i].getAlive()==1 )
		{
			plant_pop_[i].runAnnualProcesses();
		}
	}
	
	runMortalityAnnual();
	runSeedProd();
	seed_bank_[0].runMutation();	// tree seed bank    NOTE mutation now done in crossover
	seed_bank_[0].runCrossover();	// tree seed bank
	seed_bank_[1].runMutation();	// grass seed bank   NOTE mutation now done in crossover
	seed_bank_[1].runCrossover();	// grass seed bank
	runRecruitment(soil );
----------------------------------------------- */
}
//-----------------------------------------------

#ifdef NC_OUTPUT

// trait matrix mll_array to collect the trait vector over all individuals up to MAX_POP_SIZE is created in function initialize in PlantPopClass, 
// at the beginning of each new grid cell loop where this function gets called in main program.
// destruction of trait matrix mll_array in function finalize() in PlantPopClass, which is called at end of grid cell loop in main program.


// loop over all individuals up to MAX_POP_SIZE and put the trait vectors into the trait matrix


if ((year%N_BUFF_YEARS==0))		// output trait data every N_BUFF_YEARS years during the run
{
	data_pointer("time","year fraction","a",NC_FLOAT)[0]   = (float) year;
	//define here the data pointers for the variables we want to write out
	float* alive_ptr = data_pointer("alive","plant life status", "unitless",NC_INT);
	float* VegType_ptr = data_pointer("VegType","vegetation type (0: trees; 1: grass/herbaceous", "unitless",NC_INT);
	float* active_ptr = data_pointer("active","activity status (0 or 1)", "unitless",  NC_INT);
	float* RainLight_ptr = data_pointer("RainLight","phenology control by rain or light (0: rain; 1: light)", "unitless",  NC_FLOAT);
	float* Evergreen_ptr = data_pointer("Evergreen","0: deciduous; 1: evergreen", "unitless",NC_FLOAT);
	float* Sla_ptr = data_pointer("Sla","specific leaf area", "m2/kg",NC_FLOAT);
	float* AllocRoot_ptr = data_pointer("AllocRoot","allocation to root compartment", "fraction",NC_FLOAT);
	float* AllocLeaf_ptr = data_pointer("AllocLeaf","allocation to leaf compartment", "fraction",NC_FLOAT);
	float* RainThrOn_ptr = data_pointer("RainThrOn","	water threshold for leaf-on", "CHECK",NC_FLOAT);
	float* RainThrOff_ptr = data_pointer("RainThrOff","water threshold for leaf-off", "CHECK",NC_FLOAT);
	float* LightThrOn_ptr = data_pointer("LightThrOn","light threshold for leaf-on", "MJ/m2/day",NC_FLOAT);
	float* LightThrOff_ptr = data_pointer("LightThrOff","light threshold for leaf-off", "MJ/m2/day",NC_FLOAT);
	float* BiomassPar1_ptr = data_pointer("BiomassPar1","ratio of height to diameter growth", "CHECK",NC_FLOAT);
	float* BiomassPar2_ptr = data_pointer("BiomassPar2","ratio of height to diameter growth", "CHECK",NC_FLOAT);
	float* RotfrmPar1_ptr = data_pointer("RotfrmPar1","root form, par1 (slope)", "CHECK",NC_FLOAT);
	float* RotfrmPar2_ptr = data_pointer("RotfrmPar2","root form, par2 (top/deep roots)", "CHECK",NC_FLOAT);
	float* RotfrmMaxD_ptr = data_pointer("RotfrmMaxD","root form, maximum depth", "m",NC_FLOAT);
	float* WoodDensity_ptr = data_pointer("WoodDensity","wood density at 15%% moisture content", "kg/m3",NC_FLOAT);
	float* LightExt_ptr = data_pointer("LightExt","canopy light extinction parameter", "CHECK",NC_FLOAT);
	float* AllocWood_ptr = data_pointer("AllocWood","allocation to stem compartment", "fraction",NC_FLOAT);
	float* AllocBark_ptr = data_pointer("AllocBark","allocation to bark compartment", "fraction",NC_FLOAT);
	float* AllocRepr_ptr = data_pointer("AllocRepr","allocation to reproduction", "fraction",NC_FLOAT);
	float* AllocStor_ptr = data_pointer("AllocStor","allocation to storage compartments", "fraction",NC_FLOAT);
	float* SeedWeight_ptr = data_pointer("SeedWeight","seed weight", "kg",NC_FLOAT);
	float* StorToWood_ptr = data_pointer("StorToWood","allocation from storage to stem after fire", "fraction",NC_FLOAT);
	float* StorToLeaf_ptr = data_pointer("StorToLeaf","allocation from storage to leaves afer fire", "fraction",NC_FLOAT);
	float* StemDiamTot_ptr = data_pointer("StemDiamTot","stem diameter", "m",NC_FLOAT);
	float* Height_ptr = data_pointer("Height","plant height", "m",NC_FLOAT);
	float* BLeaf_ptr = data_pointer("BLeaf","leaf biomass", "kg/individual",NC_FLOAT);
	float* BWood_ptr = data_pointer("BWood","stem heartwood biomass", "kg/individual",NC_FLOAT);
	float* BRoot_ptr = data_pointer("BRoot","root biomass", "kg/individual",NC_FLOAT);
	float* CanfrmPar1_ptr = data_pointer("CanfrmPar1","canopy shape parameter (radius/height)", "unitless",NC_FLOAT);
	float* CanfrmPar2_ptr = data_pointer("CanfrmPar2","canopy shape parameter (canopy slope)", "CHECK",NC_FLOAT);
	float* BBark_ptr = data_pointer("BBark","biomass stored in bark/defense", "kg/individual",NC_FLOAT);
	float* BRepr_ptr = data_pointer("BRepr","biomass stored in reproduction compartments", "kg/individual",NC_FLOAT);
	float* BStor_ptr = data_pointer("BStor","biomass in storage compartments for resprouting", "kg/individual",NC_FLOAT);
	float* CBalance_ptr = data_pointer("CBalance","total C balance of individual", "kg/individual",NC_FLOAT);
	float* CrownArea_ptr = data_pointer("CrownArea","canopy area 0", "m2",NC_FLOAT);
	float* Lai_ptr = data_pointer("Lai","leaf area index", "m2/m2",NC_FLOAT);
	float* SpeciesIndex_ptr = data_pointer("SpeciesIndex","species index", "unitless",NC_INT);
	float* IndWaterAvail_ptr = data_pointer("IndWaterAvail","Individual water availability", "fraction",NC_FLOAT);
	float* age_ptr = data_pointer("age","age of individual, in years", "years",NC_FLOAT);
	float* CanfrmHBase_ptr = data_pointer("CanfrmHBase","Canopy base height", "fraction",NC_FLOAT);
	float* StemCount_ptr = data_pointer("StemCount","Number of stems", "unitless",NC_FLOAT);
	float* BStem_ptr = data_pointer("BStem","heartwood biomass stored in stem", "kg/individual",NC_FLOAT);
	float* BBranch_ptr = data_pointer("BBranch","woody biomass stored in branches", "kg/individual",NC_FLOAT);
	float* StorToStem_ptr = data_pointer("StorToStem","allocation of storage to stem after fire", "fraction",NC_FLOAT);
	float* StorToBranch_ptr = data_pointer("StorToBranch","allocation of storage to branches after fire", "fraction",NC_FLOAT);


	for (int i=0; i<MAX_POP_SIZE; i++)
		{
		if (!(gridcell->getplant()->max_root_depth[0]>=0 && soil[0].twp>=0)) continue;
		//if we get until here, we are in a real gridcell and now write the data
		alive_ptr[i]=plant_pop_[i].getAlive();

		if (plant_pop_[i].getAlive() != 1) continue;

		//if we get until here, we are in a real gridcell with an alive individual now write the trait data
		VegType_ptr[i]=		plant_pop_[i].getVegType();
		active_ptr[i]=		plant_pop_[i].getActive();
		RainLight_ptr[i]=	plant_pop_[i].getRainLight();
		Evergreen_ptr[i]=	plant_pop_[i].getEvergreen();
		Sla_ptr[i]=		plant_pop_[i].getSla();
		AllocRoot_ptr[i]=	plant_pop_[i].getsd()->getAllocRoot();
		AllocLeaf_ptr[i]=	plant_pop_[i].getsd()->getAllocLeaf();
		RainThrOn_ptr[i]=	plant_pop_[i].getsd()->getRainThrOn();
		RainThrOff_ptr[i]=	plant_pop_[i].getsd()->getRainThrOff();
		LightThrOn_ptr[i]=	plant_pop_[i].getsd()->getLightThrOn();
		LightThrOff_ptr[i]=	plant_pop_[i].getsd()->getLightThrOff();
		BiomassPar1_ptr[i]=	plant_pop_[i].getsd()->getBiomasPar1();
		BiomassPar2_ptr[i]=	plant_pop_[i].getsd()->getBiomasPar2();
		RotfrmPar1_ptr[i]=	plant_pop_[i].getsd()->getRotfrmPar1();
		RotfrmPar2_ptr[i]=	plant_pop_[i].getsd()->getRotfrmPar2();
		RotfrmMaxD_ptr[i]=	plant_pop_[i].getsd()->getRotfrmMaxD();
		WoodDensity_ptr[i]=	plant_pop_[i].getWoodDensity()/1.54; // double, wood density, kg m-3 at 15% moisture content (->1.54)
		LightExt_ptr[i]=	plant_pop_[i].getLightExt();
		AllocWood_ptr[i]=	plant_pop_[i].getsd()->getAllocWood();
		AllocBark_ptr[i]=	plant_pop_[i].getsd()->getAllocBark();
		AllocRepr_ptr[i]=	plant_pop_[i].getsd()->getAllocRepr();
		AllocStor_ptr[i]=	plant_pop_[i].getsd()->getAllocStor();
		SeedWeight_ptr[i]=	plant_pop_[i].getsd()->getSeedWeight();
		StorToWood_ptr[i]=	plant_pop_[i].getsd()->getStorToWood();
		StorToLeaf_ptr[i]=	plant_pop_[i].getsd()->getStorToLeaf();
		StemDiamTot_ptr[i]=	plant_pop_[i].getStemDiam();
		Height_ptr[i]=		plant_pop_[i].getHeight();
		BLeaf_ptr[i]=		plant_pop_[i].getBLeaf();
		BWood_ptr[i]=		plant_pop_[i].getBWood();
		BRoot_ptr[i]=		plant_pop_[i].getBRoot();
		CanfrmPar1_ptr[i]=	plant_pop_[i].getsd()->getCanfrmPar1();
		CanfrmPar2_ptr[i]=	plant_pop_[i].getsd()->getCanfrmPar2();
		BBark_ptr[i]=		plant_pop_[i].getBBark();
		BRepr_ptr[i]=		plant_pop_[i].getBRepr();
		BStor_ptr[i]=		plant_pop_[i].getBStor();
		CBalance_ptr[i]=	plant_pop_[i].getCBalance();
		CrownArea_ptr[i]=	plant_pop_[i].getCrownArea();
		Lai_ptr[i]=		plant_pop_[i].getLai();
		SpeciesIndex_ptr[i]=	plant_pop_[i].getSpeciesIndex();
		IndWaterAvail_ptr[i]=	plant_pop_[i].getWaterAvail();
		age_ptr[i]=		plant_pop_[i].getAge();
		CanfrmHBase_ptr[i]=	plant_pop_[i].getsd()->getCanfrmHBase();
		StemCount_ptr[i]=	plant_pop_[i].getsd()->getStemCount();
		BStem_ptr[i]=		plant_pop_[i].getBStem();
		#ifdef W_BRANCHES
			BBranch_ptr[i]=		plant_pop_[i].getBBranch();
		#endif 	
		StorToStem_ptr[i]=	plant_pop_[i].getsd()->getStorToStem();
		StorToBranch_ptr[i]=	plant_pop_[i].getsd()->getStorToBranch();

		
		}
bufferoverflow=1;
}

#endif //NC_OUTPUT
 
//-----------------------------------------------
      
return bufferoverflow;
}

//--------------------------------------------------------------------------------------

void clPlantPop::runSeedProd()
{	
	seed_bank_[0].emptySeedBank();  // annual seed production, old seeds die
	seed_bank_[1].emptySeedBank();  // annual seed production, old seeds die
	
	// 	cout << "--0------ " << seed_bank_[0].getSeedNumber() << " " << seed_bank_[1].getSeedNumber() << endl;
	
	float seed_nums_0[MAX_POP_SIZE];
	float seed_nums_1[MAX_POP_SIZE];
	for ( int i=0; i<MAX_POP_SIZE; i++ ) seed_nums_0[i] = 0;
	for ( int i=0; i<MAX_POP_SIZE; i++ ) seed_nums_1[i] = 0;

	float seed_prod[VTYPE_NUM];
	for ( int i=0; i<VTYPE_NUM; i++ ) seed_prod[i] = 0;
	
	
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{
		if ( plant_pop_[i].getAlive()==1 )
		{
			if ( plant_pop_[i].getVegType()==TREE )
			{
				seed_nums_0[i] = plant_pop_[i].getSeedProd();
				seed_prod[ plant_pop_[i].getVegType() ] += seed_nums_0[i];
			}
			else if ( plant_pop_[i].getVegType()==GR_4 )
			{
				seed_nums_1[i] = plant_pop_[i].getSeedProd();
				seed_prod[ plant_pop_[i].getVegType() ] += seed_nums_1[i];
			}
			//seed_prod[ plant_pop_[i].getVegType() ] += plant_pop_[i].getBRepr();
		}
	}
	
	//cout << setw(15) << seed_prod[0] << setw(15) << seed_prod[1] << endl;
	
	if ( seed_prod[0] >= MAX_SEED_BANK_SIZE )
	{
		double tmp_fac = (double)MAX_SEED_BANK_SIZE/(double)seed_prod[0];
		seed_prod[0] = 0;
		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			seed_nums_0[i] = (int)floor( seed_nums_0[i]*tmp_fac );
			seed_prod[0] += seed_nums_0[i];
		}
		//cout << "New tree seed bank size  " << setw(8) << seed_prod[0] << endl;
	}

	if ( seed_prod[1] >= MAX_SEED_BANK_SIZE )
	{
		double tmp_fac = (double)MAX_SEED_BANK_SIZE/(double)seed_prod[1];
		seed_prod[1] = 0;
		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			seed_nums_1[i] = (int)floor( seed_nums_1[i]*tmp_fac );
			seed_prod[1] += seed_nums_1[i];
		}
		//cout << "New grass seed bank size " << setw(8) << seed_prod[1] << endl;
	}

	//for ( int i=0; i<MAX_POP_SIZE; i++ ) cout << setw(4) << seed_nums_0[i] << setw(4) << seed_nums_1[i] << setw(4) << plant_pop_[i].getVegType() << endl;

	for ( int sp=0; sp<SPECIES_NUM; sp++ )
	{
		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			if ( plant_pop_[i].getSpeciesIndex()==sp && seed_nums_0[i]>0 )
			{
				
				for ( int j=0; j<seed_nums_0[i]; j++ ) seed_bank_[0].addSeed( plant_pop_[i].getSeed() );
			}
		}
	}

	for ( int sp=0; sp<SPECIES_NUM; sp++ )
	{
		for ( int i=0; i<MAX_POP_SIZE; i++ )
		{
			if ( plant_pop_[i].getSpeciesIndex()==sp && seed_nums_1[i]>0 )
			{
				for ( int j=0; j<seed_nums_1[i]; j++ ) seed_bank_[1].addSeed( plant_pop_[i].getSeed() );
			}
		}
	}
	
	return;
}


//--------------------------------------------------------------------------------------

void clPlantPop::runMortalityDaily()
{
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{
		if (plant_pop_[i].getAlive()==1 && pop_size_[0]>1 ) // one individual remains alive to avoid nan
		{
			if ( plant_pop_[i].runMortalityDaily()==1 )
			{
				pop_size_[0]--; 
			}
		}
	}
	
	return;
}

//--------------------------------------------------------------------------------------

void clPlantPop::runMortalityAnnual()
{
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{
		if (plant_pop_[i].getAlive()==1 && pop_size_[0]>1 ) // one individual remains alive to avoid nan
		{
			if ( plant_pop_[i].runMortalityAnnual()==1 )
			{
				pop_size_[0]--;
			}
		}
	}
	
	return;
}

//--------------------------------------------------------------------------------------

void clPlantPop::runRecruitment(clGridCell *gridcell)
{
    GridCellSoilData* soil =gridcell->getsoil();
	int veg_type = TREE ;
	int seed_num = seed_bank_[TREE].getSeedNumber()+seed_bank_[GR_4].getSeedNumber(); // total seed number
// 	cout << "SEEDS " << setw(10) << seed_bank_[0].getSeedNumber() << setw(10) << seed_bank_[1].getSeedNumber() << endl;
// 	float type_frac = (float)seed_bank_[0].getSeedNumber()/(float)seed_num;  // probability for veg type 0
	
// 	float type_frac = getMeanVegType( 0 );
// 	float type_frac = MyMax( MyMin( getTreeFraction()/(getTreeFraction()+getC4GrassFraction()), 0.9 ), 0.1 );
	float type_frac = 0.5;  // NOTE simply used 0.5 here; 50 % of the recruited seedlings are trees
//	float type_frac = 0. ;	// 0 % of the recruits are trees (?)
// 	float type_frac = getTreeFraction();
	int compinds_gap[NUM_LIGHT_NGB];

	
	// MAX_POP_SIZE-pop_size_ gives the number of dead individuals/empty spots in plant pop vector
	// These empty spots are distributed to grass and tree recruits
	int tree__recruits = type_frac*(MAX_POP_SIZE-pop_size_[0]);
	int grass_recruits = (MAX_POP_SIZE-pop_size_[0])-tree__recruits;
	
	if ( tree__recruits > seed_bank_[TREE].getSeedNumber() )
	{
		int diff = tree__recruits - seed_bank_[TREE].getSeedNumber();
		tree__recruits -= diff;
		grass_recruits += diff;
	}
	
	if ( grass_recruits > seed_bank_[GR_4].getSeedNumber() )
	{
		int diff = grass_recruits - seed_bank_[GR_4].getSeedNumber();
		tree__recruits += diff;
		grass_recruits -= diff;
	}
	
// 	cout << "SEEDS " << setw(10) << seed_bank_[TREE].getSeedNumber() << setw(10) << seed_bank_[GR_4].getSeedNumber() << setw(10) << seed_num << setw(10) << pop_size_[0] << setw(10) << (MAX_POP_SIZE-pop_size_[0]) << setw(14) << type_frac << setw(14) << tree__recruits << setw(14) << grass_recruits << endl;
//      FLAG: currently all seeds get planted in the empty cellc starting in the "upper left" corner of stand "square"	
//	FLAG: Rethink this planting concept
	if ( seed_num > 0)
	{
		int    cnt_pop   = 0;
		int	   cnt_germ_success = 0; 	// counter to keep track how many seeds have germinated
		int	   empty_spots = tree__recruits+grass_recruits;
		while ( cnt_pop<MAX_POP_SIZE && tree__recruits+grass_recruits>0 )  // go once through the entire pop-grid; stop either once the end has been reached, or if there are no more seeds that can be planted
		{
			if( plant_pop_[cnt_pop].getAlive()==0 )  // look for empty spots in the pop-grid that can be planted with new seeds
			{			
  			  	double compheight_crown_area_gap = 0.0;	
// 				cout <<  " compheight_crown_area_gap  " << compheight_crown_area_gap << endl;
    
				compinds_gap[0] = (cnt_pop-MAX_POP_SIZE_SQRT+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[1] = (cnt_pop+MAX_POP_SIZE_SQRT+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[2] = (cnt_pop-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[3] = (cnt_pop+1+MAX_POP_SIZE)%MAX_POP_SIZE;
			
				compinds_gap[4] = (cnt_pop-MAX_POP_SIZE_SQRT-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[5] = (cnt_pop+MAX_POP_SIZE_SQRT-1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[6] = (cnt_pop-MAX_POP_SIZE_SQRT+1+MAX_POP_SIZE)%MAX_POP_SIZE;
				compinds_gap[7] = (cnt_pop+MAX_POP_SIZE_SQRT+1+MAX_POP_SIZE)%MAX_POP_SIZE;

			  for ( int cmp=0; cmp<NUM_LIGHT_NGB; cmp++ )
			  {
				if ( plant_pop_[compinds_gap[cmp]].getVegType()==0 ) compheight_crown_area_gap += plant_pop_[compinds_gap[cmp]].getCrownArea()*(double)plant_pop_[compinds_gap[cmp]].getActive();
// 					cout << " compheight_crown_area_gap_" << cmp << "  " << compheight_crown_area_gap << "  plant_pop_[compinds_gap[cmp]].getCrownArea0()  " << plant_pop_[compinds_gap[cmp]].getCrownArea0() << "  (double)plant_pop_[compinds_gap[cmp]].getActive()  " << (double)plant_pop_[compinds_gap[cmp]].getActive() << 
//                  "  plant_pop_[compinds_gap[cmp]].getCrownArea0()*(double)plant_pop_[compinds_gap[cmp]].getActive()  " << plant_pop_[compinds_gap[cmp]].getCrownArea0()*(double)plant_pop_[compinds_gap[cmp]].getActive()<< endl;
			  }
			  
			  
				if ( RUnif() < type_frac ) veg_type = TREE;	// decide if the seed to be planted will be a tree seed or a grass seed, based on the fraction of the pop constituted by trees (not sure this is the best strategy...)
				else                       veg_type = GR_4;
				
				if ( veg_type==TREE && tree__recruits==0 ) veg_type = GR_4;  // if we run out of tree seeds (grass seeds), switch to grass seeds (tree seeds); this is NOT circular, since the while-loop prevents both seed pools
				if ( veg_type==GR_4 && grass_recruits==0 ) veg_type = TREE;		// from being ZERO at the same time!
				
				if ( veg_type==TREE ) tree__recruits--;		// decrease the respective seed pool for the seed type that got selected to be planted
				if ( veg_type==GR_4 ) grass_recruits--;
				
//				cout << setw(14) << tree__recruits << setw(14) << grass_recruits << endl;
								
				clSeed tmp_seed = seed_bank_[veg_type].getRandomSeed();		// draw a random seed from the seed pool		
				
// 				if ( RUnif()<(0.05+25.*tmp_seed.getSeedWeight()) ) // germination probability of the drawn random seed,based on its weight;
// 				if ( RUnif()<(0.2+1.*tmp_seed.getSeedWeight()) ) // germination probability of the drawn random seed,based on its weight; 	
				if ( RUnif()<(0.04+1.*tmp_seed.getSeedWeight()) - (0.15*MyMin(1., compheight_crown_area_gap/14.4)) ) // germination probability of the drawn random seed,based on its weight; 
				{
					plant_pop_[cnt_pop].setInitialValues( tmp_seed, soil );     // set traits from seed // state variables to initial val.
					pop_size_[0]++;
					cnt_germ_success++;
				//	cout << "add " << veg_type << setw(14) << tmp_seed.getSeedWeight() << setw(14) << 10.*tmp_seed.getSeedWeight() << endl;
				}
			}
			cnt_pop++;
		}
//		cout << "Germination rate (%): " << ((double) cnt_germ_success/(double) empty_spots) * 100. << endl;
	}
	
	return;
}


// ---------------------------------------------------------------------
// --- Calculate means of plant population -----------------------------
// ---------------------------------------------------------------------

//--------------------------------------------------------------------------------------

// float clPlantPop::getMeanLai()		// NEW LIAM, THIS ENTIRE FUNCTION WAS COMMENTED IN BEFORE!
// {
//   float tmp=0;
//   for ( int i=0; i<MAX_POP_SIZE; i++ ) 
//		tmp += ( plant_pop_[i].getLai()*plant_pop_[i].getAlive() );
//    return tmp/(float)pop_size_[0];
// }
//--------------------------------------------------------------------------------------

//float clPlantPop::getMeanWoodDensity()
//{
	//float tmp=0.;
	//for (int i=0; i<MAX_POP_SIZE; i++)
	//{
		//if(plant_pop_[i].getVegType()==TREE && plant_pop_[i].getAlive()==1) tmp += plant_pop_[i].getWoodDensity() ;
	//}
	//tmp = tmp/(pop_size_[1] * 1.54) ;	// 1.54 => convert from wood density @ 50% moisture to wood density @ 15 % moisture
	
////	cout << "Average wood density: " << year_ <<setw(14)<< tmp <<setw(14)<< pop_size_[1] <<setw(14)<< pop_size_[2] <<setw(14)<< pop_size_[0] << endl;
	
	//return tmp ;
//}
//--------------------------------------------------------------------------------------

float clPlantPop::getMeancanopy_h_base(float *tmp)	// NEW CAMILLE	
{
   
   //double tmp=0;
    //for ( int i=0; i<MAX_POP_SIZE; i++ ) 
		//tmp += ( plant_pop_[i].getCAN8A()*plant_pop_[i].getAlive() );
	////cout<<tmp<<" " <<plant_pop_[0].getCAN8A()<< " "<<pop_size_<<endl;
    //return tmp/(double)pop_size_;
   
        
    for (int i=0; i<3; i++) tmp[i]=0;	
    for (int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()==1 ) tmp[0] += ( plant_pop_[i].getCrown_h_base()*plant_pop_[i].getAlive() );
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getCrown_h_base()*plant_pop_[i].getAlive() );
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getCrown_h_base()*plant_pop_[i].getAlive() );
	//cout<<tmp<<" " <<plant_pop_[0].getCanopy_h_base()<< " "<<pop_size_<<endl;
	}
  // return tmp[i]/(float)pop_size_[i];
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
	return tmp[0] ;
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanStem(float *tmp) // NEW CAMILLE
{
          
    for (int i=0; i<3; i++) tmp[i]=0;	
    for (int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()==1 ) tmp[0] += ( plant_pop_[i].getStem_count()*plant_pop_[i].getAlive() );
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getStem_count()*plant_pop_[i].getAlive() );
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getStem_count()*plant_pop_[i].getAlive() );
	//cout<<tmp<<" " <<plant_pop_[0].getStem_count()<< " "<<pop_size_<<endl;
	}
  // return tmp[i]/(float)pop_size_[i];
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
	return tmp[0] ;
}


//--------------------------------------------------------------------------------------

float clPlantPop::getMeanWaterAvail(float *tmp)
{

	for (int i=0; i<3; i++) tmp[i]=0;
   
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		tmp[0] += ( plant_pop_[i].getWaterAvail()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );		// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getWaterAvail()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );	// position 1 holds the average accross all trees of the population	
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getWaterAvail()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );	// position 2 holds the average accross all C4-grasses of the population	
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	   
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanGi(float *tmp)
{

	for (int i=0; i<3; i++) tmp[i]=0;
    
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{ 
		tmp[0] += ( plant_pop_[i].getGi_avg()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );											// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getGi_avg()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );	// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getGi_avg()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );	// position 2 holds the average accross all C4-grasses of the population
	}	
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanLai(float *tmp)		// NEW LIAM
{

	for (int i=0; i<3; i++) tmp[i]=0;
    
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		tmp[0] += ( plant_pop_[i].getBLeaf()*plant_pop_[i].getSLA()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );		// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getBLeaf()*plant_pop_[i].getSLA()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );		// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getBLeaf()*plant_pop_[i].getSLA()*plant_pop_[i].getAlive()*plant_pop_[i].getActive() );		// position 2 holds the average accross all C4-grasses of the population
	}	
	
	for (int i=0; i<3; i++) tmp[i] = tmp[i] /(((float)pop_size_[i]/(float)MAX_POP_SIZE)*PLOT_SIZE); 	 // MP: for now scale LAI to the fraction of total plot size assigned to the specific vegtype according to its number of individuals
	
	return tmp[0];	
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanHeight(float *tmp)
{

	for (int i=0; i<3; i++) tmp[i]=0;
    
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		if (plant_pop_[i].getAlive()!=1) continue;
		tmp[0] += ( plant_pop_[i].getHeight());					// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) {tmp[1] += ( plant_pop_[i].getHeight());}// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) {tmp[2] += ( plant_pop_[i].getHeight());}	// position 2 holds the average accross all C4-grasses of the population	
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0)	tmp[i] = tmp[i]/(float) pop_size_[i] ; 
		else			tmp[i] = FMISSING ; 	
	}
	
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanStemDiam(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
    
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		tmp[0] += ( plant_pop_[i].getStemDiam()*plant_pop_[i].getAlive() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  ( plant_pop_[i].getStemDiam()*plant_pop_[i].getAlive() );		// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] +=  ( plant_pop_[i].getStemDiam()*plant_pop_[i].getAlive() );		// position 2 holds the average accross all C4-grasses of the population		
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
		
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanWiltPoint(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
   
    for ( int i=0; i<MAX_POP_SIZE; i++ )
    { 
		tmp[0] += ( plant_pop_[i].getWiltPoint()*plant_pop_[i].getAlive() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  ( plant_pop_[i].getWiltPoint()*plant_pop_[i].getAlive() );		// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += 	( plant_pop_[i].getWiltPoint()*plant_pop_[i].getAlive() );		// position 2 holds the average accross all C4-grasses of the population	
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
		
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanActive(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
   
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		tmp[0] += ( plant_pop_[i].getActive() *plant_pop_[i].getAlive() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getActive() *plant_pop_[i].getAlive() ); 		// position 1 holds the average accross all trees of the population
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getActive() *plant_pop_[i].getAlive() );			// position 2 holds the average accross all C4-grasses of the population	
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
		
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanVegType( int type )		
{
    float tmp=0;
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
        if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==type ) tmp += 1.;
        
        if(pop_size_[0] != 0)
        	{tmp = tmp/(float)pop_size_[0]; }
        else
        	{tmp = FMISSING; }	
//	if (tmp> 1){ cout << " Error meanvegtype unplausible "<< tmp<<" popsize: " <<pop_size_[0]<<" "<<pop_size_[1]<<" "<<pop_size_[2]<<" runid:"<<YEARS_TO_RUN<<lat_<<endl;}
    
    return tmp ;
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanCanopyArea0(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
   
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
//     	if (plant_pop_[i].getAlive()==1 ) tmp[0] += plant_pop_[i].getCrownArea();		// BEFORE LIAM CHANGES			// position zero holds the average accross the entire population
        if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getCrownArea(); 	//NOTE: This is for trees only like this			// NEW LIAM
        if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getCrownArea() * IND_AREA_MULT ;  		// position 2 holds the average accross all C4-grasses of the population      
//         tmp += ( plant_pop_[i].getCrownArea()*plant_pop_[i].getAlive() );
	}
	
	tmp[0] = tmp[1] + tmp[2] ;	
    	
	for (int i=0; i<3; i++)
	{
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	} 
	
    return tmp[1]/(float) pop_size_[0];  //FLAG2 //NEWLIAM this is for trees only but divided with TREE+GRASS complete popsize
}

//--------------------------------------------------------------------------------------


float clPlantPop::getMeanCGain(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
   
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
//		tmp[0] += ( plant_pop_[i].getCGain()*plant_pop_[i].getAlive() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getCGain()); 		// position 1 holds the average accross all trees of the population   
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getCGain()) * IND_AREA_MULT ;   	// position 2 holds the average accross all C4-grasses of the population 				
	}

	tmp[0] = tmp[1] + tmp[2] ;
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanCGain2(float *tmp) //STEVE this is per plant
{
 
    for (int i=0; i<3; i++) tmp[i]=0;
   
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
//		tmp[0] += ( plant_pop_[i].getCGain2()*plant_pop_[i].getAlive() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( plant_pop_[i].getCGain2()); 		// position 1 holds the average accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( plant_pop_[i].getCGain2()) * IND_AREA_MULT ;   	// position 2 holds the average accross all C4-grasses of the population		
	}

	tmp[0] = tmp[1] + tmp[2] ;
		
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getCGain_Plot(float *tmp) //STEVE this is per stand ie CGainplant[i]*canopy_area_plant[i]/total_area_of_plot
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
//      if (plant_pop_[i].getAlive()==1 ) tmp[0] += (plant_pop_[i].getCGain2() * plant_pop_[i].getCrownArea());		// SIMON replaced Z by 0; position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += (plant_pop_[i].getCGain2() * plant_pop_[i].getCrownArea()); 	// SIMON replaced Z by 0; position 1 holds the average accross all trees of the population 
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += (plant_pop_[i].getCGain2() * plant_pop_[i].getCrownArea()) * IND_AREA_MULT ;    // SIMON replaced Z by 0; position 2 holds the average accross all C4-grasses of the population      
    }

	tmp[0] = tmp[1] + tmp[2] ;
	    
    for (int i=0; i<3; i++) tmp[i] = tmp[i] /(float)PLOT_SIZE;
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanRainLight(float *tmp)
{
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		tmp[0] += ( MyRound(plant_pop_[i].getRainLight()) );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += ( MyRound(plant_pop_[i].getRainLight()) ); 		// position 1 holds the average accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += ( MyRound(plant_pop_[i].getRainLight()) );    	// position 2 holds the average accross all C4-grasses of the population 		
	}	
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
   return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanAlive(float *tmp)
{
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		tmp[0] +=  plant_pop_[i].getAlive() ;			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  plant_pop_[i].getAlive() ; 			// position 1 holds the average accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += 	 plant_pop_[i].getAlive() ;	    	// position 2 holds the average accross all C4-grasses of the population 		
	}	
	
	for (int i=0; i<3; i++) tmp[i] = tmp[i] /(float)MAX_POP_SIZE;
	
    return tmp[0];		// FLAG: division by MAX_POP_SIZE, desired like that?
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanStomCond(float *tmp)
{
    for (int i=0; i<3; i++) tmp[i]=0;
  
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
// 		tmp += ( plant_pop_[i].getStomCond()*plant_pop_[i].getAlive() );
		tmp[0] += ( plant_pop_[i].getStomCond()*plant_pop_[i].getCrownArea() );			// position zero holds the average accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  ( plant_pop_[i].getStomCond()*plant_pop_[i].getCrownArea() ); 		// position 1 holds the average accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] +=  ( plant_pop_[i].getStomCond()*plant_pop_[i].getCrownArea() );    	// position 2 holds the average accross all C4-grasses of the population 	
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}

// 	return tmp/getSumCanopyArea0();
	return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getMeanAge(float *tmp)		// NEW MIRJAM; 04.02.2015
{
    for (int i=0; i<3; i++) tmp[i]=0;	

	for (int i=0; i<MAX_POP_SIZE; i++)
	{
	if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
        tmp[0] += plant_pop_[i].getAge();	
        if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getAge();
        if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getAge();        		
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
	return tmp[0] ;
}

//--------------------------------------------------------------------------------------

float clPlantPop::getShoot2Root(float *tmp)		// NEW MIRJAM; 05.02.2015
{
	float temp[3];	// this to temporarily store the root biomass
	
    for (int i=0; i<3; i++) { tmp[i]=0;	temp[i]=0;}

	for (int i=0; i<MAX_POP_SIZE; i++)
	{
	if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
        {tmp[0] += plant_pop_[i].getBLeaf() + plant_pop_[i].getBWood() + plant_pop_[i].getBBark(); temp[0] += plant_pop_[i].getBRoot();}
        if (plant_pop_[i].getVegType()==TREE) { tmp[1] += plant_pop_[i].getBLeaf() + plant_pop_[i].getBWood() + plant_pop_[i].getBBark(); temp[1] += plant_pop_[i].getBRoot(); }
        if (plant_pop_[i].getVegType()==GR_4) { tmp[2] += plant_pop_[i].getBLeaf() ; temp[2] += plant_pop_[i].getBRoot(); }       		
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) temp[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
	return tmp[0] ;	
}	
	
//--------------------------------------------------------------------------------------

float clPlantPop::getMeanSla(float *tmp)		// NEW MIRJAM; 16.02.2015	
{
    for (int i=0; i<3; i++) tmp[i]=0;	

	for (int i=0; i<MAX_POP_SIZE; i++)
	{
	if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
        tmp[0] += plant_pop_[i].getSLA();	
        if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getSLA();
        if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getSLA();        		
	}
	
	for (int i=0; i<3; i++) 
	{	
		if(pop_size_[i] != 0) 
			{ tmp[i] = tmp[i] /(float) pop_size_[i] ; }
		else	
		    { tmp[i] = FMISSING ; }	
	}
	
	return tmp[0] ;		
}


//--------------------------------------------------------------------------------------

// ---------------------------------------------------------------------
// --- Calculate sums of plant population ------------------------------
// ---------------------------------------------------------------------


float clPlantPop::getSumEvapotr(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
	for ( int i=0; i<MAX_POP_SIZE; i++ )
	{
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  ( plant_pop_[i].getEvapotr());	 		// position 1 holds the sum accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] +=  ( plant_pop_[i].getEvapotr()) * IND_AREA_MULT ;	    	// position 2 holds the sum accross all C4-grasses of the population		
	}
	
	tmp[0] = tmp[1] + tmp[2] ;
	
//	for (int i=0; i<3; i++) tmp[i] = tmp[i] /(float)PLOT_SIZE;			// mm ha-1?  

//	for (int i=0; i<3; i++) tmp[i] = tmp[i] / (((float)pop_size_[i]/(float)MAX_POP_SIZE)*PLOT_SIZE);   	// MP: for now scale each vegtype's ET to the fraction of total plot size assigned to the specific vegtype
	for (int i=0; i<3; i++) tmp[i] = tmp[i] / PLOT_SIZE ; 	// mm m-2 (?); plot-level average for each veg-type; sum of all veg-types = total mean plot level evapotr., in unit mm m-2 day-1																								
		
	return tmp[0];		// FLAG: Check if this is actually correct! I 
// 	return tmp*SIZE_FACTOR*0.001; // from m to mm
// 	return tmp/getSumCanopyArea0();
// 	return tmp/pop_size_[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBLeaf(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBLeaf();	 		// position 1 holds the sum accross all trees of the population 
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBLeaf() * IND_AREA_MULT ;    	 	// position 2 holds the sum accross all grasses of the population   
    }
    
    tmp[0] = tmp[1] + tmp[2] ;
    
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBStem(float *tmp)
{
  
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBStem();	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBStem() * IND_AREA_MULT ;    		// position 2 holds the sum accross all C4-grasses of the population		
	} 
    
    tmp[0] = tmp[1] + tmp[2] ;
    	
	for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.;        
       
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBWood(float *tmp)
{
  
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBWood();	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBWood() * IND_AREA_MULT ;    		// position 2 holds the sum accross all C4-grasses of the population		
	} 
    
    tmp[0] = tmp[1] + tmp[2] ;
    	
	for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.;        
       
    return tmp[0];
}

//--------------------------------------------------------------------------------------
#ifdef W_BRANCHES

	float clPlantPop::getSumBBranch(float *tmp)
	{
  
	    for (int i=0; i<3; i++) tmp[i]=0;
  
	    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	    {
			if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
			if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBBranch();	 		// position 1 holds the sum accross all trees of the population  
			if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBBranch() * IND_AREA_MULT ;    		// position 2 holds the sum accross all C4-grasses of the population		
		} 
    
	    tmp[0] = tmp[1] + tmp[2] ;
    	
		for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.;        
       
	    return tmp[0];
}

#endif

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBRoot(float *tmp)
{
 
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] +=  plant_pop_[i].getBRoot();	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] +=  plant_pop_[i].getBRoot() * IND_AREA_MULT ;    		// position 2 holds the sum accross all C4-grasses of the population	   
    }
    
    tmp[0] = tmp[1] + tmp[2] ;
        
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBBark(float *tmp)
{
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBBark();	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBBark() * IND_AREA_MULT ;     		// position 2 holds the sum accross all C4-grasses of the population      
    }
    
    tmp[0] = tmp[1] + tmp[2] ;
        
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBRepr(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getBRepr();	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getBRepr() * IND_AREA_MULT ;     		// position 2 holds the sum accross all C4-grasses of the population        
    }    
    
    tmp[0] = tmp[1] + tmp[2] ;
        
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBStor(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] +=  plant_pop_[i].getBStor();  // position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] +=  plant_pop_[i].getBStor() * IND_AREA_MULT ;  // position 2 holds the sum accross all C4-grasses of the population          
    }
    
    tmp[0] = tmp[1] + tmp[2] ;
        
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR/1000.; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumStomCond(float *tmp)
{
   
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  (plant_pop_[i].getStomCond()*plant_pop_[i].getCrownArea());	 		// position 1 holds trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] +=  (plant_pop_[i].getStomCond()*plant_pop_[i].getCrownArea()) * IND_AREA_MULT ;      	// position 2 holds all C4-grasses of the population   	
	}
    
    tmp[0] = tmp[1] + tmp[2] ;
    	
	for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
		
//     return tmp*SIZE_FACTOR;
	return tmp[0];
// 	return tmp/getCrownArea(0.);
// 		return tmp;   // NOTE no multiplication because stom. cond. given in mmol/m^2/s
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumBasalArea(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += pow( (plant_pop_[i].getStemDiam()*0.5), 2 )*M_PI;	 		// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += pow( (plant_pop_[i].getStemDiam()*0.5), 2 )*M_PI * IND_AREA_MULT ;       	// position 2 holds the sum accross all C4-grasses of the population       
    }    
    
    tmp[0] = tmp[1] + tmp[2] ;
        
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
    
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumActive(float *tmp)
{
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		tmp[0] += plant_pop_[i].getActive();						// position zero holds the sum accross the entire population
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getActive();	// position 1 holds the sum accross all trees of the population  
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getActive();       // position 2 holds the sum accross all C4-grasses of the population        
    }
    
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
    
   return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumCanopyArea0(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
	        tmp[0] += plant_pop_[i].getCrownArea();					// position zero holds the sum accross the entire population
        	if (plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getCrownArea();	// position 1 holds the sum accross all trees of the population  
	        if (plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getCrownArea();	       // position 2 holds the sum accross all C4-grasses of the population       
    }
    
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR;  //FLAG2
    
    return tmp[1];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumCGain(float *tmp)
{
  
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
                tmp[0] += plant_pop_[i].getCGain();				// position zero holds the sum accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] +=  plant_pop_[i].getCGain(); 		// position 1 holds the sum accross all trees of the population  
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] +=  plant_pop_[i].getCGain(); 	    // position 2 holds the sum accross all C4-grasses of the population      
    }
    
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
    
	return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumCGainCum(float *tmp)
{

    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		tmp[0] += plant_pop_[i].getCGainCum();				// position zero holds the sum accross the entire population
		if(plant_pop_[i].getVegType()==TREE) tmp[1] += plant_pop_[i].getCGainCum(); 		// position 1 holds the sum accross all trees of the population   
		if(plant_pop_[i].getVegType()==GR_4) tmp[2] += plant_pop_[i].getCGainCum();	 	    // position 2 holds the sum accross all C4-grasses of the population 
	}
	
	for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
	
    return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumRainLight(float *tmp)
{
   
    for (int i=0; i<3; i++) tmp[i]=0;
  
    for ( int i=0; i<MAX_POP_SIZE; i++ ) 
    {
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		tmp[0] += MyRound(plant_pop_[i].getRainLight());					// position zero holds the sum accross the entire population
		if (plant_pop_[i].getVegType()==TREE) tmp[1] += MyRound(plant_pop_[i].getRainLight());	// position 1 holds the sum accross all trees of the population   
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] += MyRound(plant_pop_[i].getRainLight());	// position 2 holds the sum accross all C4-grasses of the population        
    }
    
    for (int i=0; i<3; i++) tmp[i] = tmp[i]*SIZE_FACTOR; 
        
	return tmp[0];
}

//--------------------------------------------------------------------------------------

float clPlantPop::getSumAlive(float *tmp)
{

	for (int i=0; i<3; i++) tmp[i]=0;
  
	for ( int i=0; i<MAX_POP_SIZE; i++ ) 
	{
		if (plant_pop_[i].getAlive()!=1) continue; //only sum over alive individuals
		tmp[0] += plant_pop_[i].getAlive();						// position zero holds the sum accross the entire population
		if (plant_pop_[i].getVegType()==TREE) tmp[1] +=  plant_pop_[i].getAlive();	// position 1 holds the sum accross all trees of the population   
		if (plant_pop_[i].getVegType()==GR_4) tmp[2] +=  plant_pop_[i].getAlive();	// position 2 holds the sum accross all C4-grasses of the population         
	}
    
	return tmp[0];		
}

//--------------------------------------------------------------------------------------

float clPlantPop::getTreeFraction()		// FLAG: could combine this function and the following one into one by using vector. Leave the way it is, or alter?
{
	float sum=0;
	for ( int i=0; i<MAX_POP_SIZE; i++ )
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==TREE ) sum += plant_pop_[i].getCrownArea();  // SIMON Replaced Z by 0
		
	return sum/PLOT_SIZE;
}

//--------------------------------------------------------------------------------------

float clPlantPop::getC4GrassFraction()
{
	float sum=0;
	for ( int i=0; i<MAX_POP_SIZE; i++ )
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 ) sum += plant_pop_[i].getCrownArea();  // SIMON Replaced Z by 0
		
	return sum/PLOT_SIZE;
}

//--------------------------------------------------------------------------------------

// C3 grasses not implemented yet
float clPlantPop::getC3GrassFraction()
{
	float sum=0;
	for ( int i=0; i<MAX_POP_SIZE; i++ )
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_3 ) sum += plant_pop_[i].getCrownArea();  // SIMON Replaced Z by 0
		
	return sum/PLOT_SIZE;
}

//--------------------------------------------------------------------------------------

// number of alive individuals at a given time, total, and by vtype

float clPlantPop::getNAliveByType(float *tmp)	
{
	
	for (int i=0; i<3; i++) tmp[i] = (float) pop_size_[i] ;  	
	    
    return tmp[0];

}





//--------------------------------------------------------------------------------------

// function to determine how much a grass individual's leaf biomass will influence grazing selection using a linear function fit to min and max grass
// leaf biomass range

void clPlantPop::bm_weight_coeffs()
{
	double bm_max_ = 0.;
	double bm_min_ = 3.4E38;
	
	for (int i=0; i<MAX_POP_SIZE; i++)		// put this in runDailyProcesses, since it doesn't have to be called for each individual separately!
	{	
//		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4) cout << "BLeaf ind: " << plant_pop_[i].getBLeaf() << endl; 
//		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getBLeaf()<bm_min_) bm_min_ = plant_pop_[i].getBLeaf();		// find the grass ind with the lowest biomass
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getBLeaf()> bm_max_) bm_max_ = plant_pop_[i].getBLeaf();		// find the grass ind with the highest biomass
	}
	
	bm_min_ = meristem_ ; 	
	
//	cout << "bm_weight_coeffs(), bm_min, bm_max: " << bm_min_ << setw(14) << bm_max_ << endl;
	
	bm_slope_ = 1. / (bm_max_ - bm_min_);			// slope of a line between (bm_min,0) and (bm_max,1)
	bm_intercept_ = 1. / (1. - bm_max_/bm_min_);	// intercept of a line between (bm_min,0) and (bm_max,1)
}

//--------------------------------------------------------------------------------------

// function to determine how much a grass individual's SLA will influence grazing selection using a quadratic function fit to min and max SLA range

void clPlantPop::sla_weight_coeffs()
{
	double sla_max_ = 0.; 
	double sla_min_ = 3.4E38 ;
	
	for (int i=0; i<MAX_POP_SIZE; i++)		// put this in runDailyProcesses, since it doesn't have to be called for each individual separately!
	{
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getSLA()< sla_min_) sla_min_ = plant_pop_[i].getSLA();		// find the grass ind with the lowest SLA
		if (plant_pop_[i].getAlive()==1 && plant_pop_[i].getVegType()==GR_4 && plant_pop_[i].getSLA()> sla_max_) sla_max_ = plant_pop_[i].getSLA();		// find the grass ind with the highest SLA
	}
	
//	cout << "sla_weight_coeffs(), sla_min, sla_max: " << sla_min_ << setw(14) << sla_max_ << endl; 

	// coefficients for a downwardly open parabola, with coordinates (sla_min_,0) and (sla_max_,1), with (sla_max_,1) being the vertex of the parabola
		
	a_ = -1. / (Mypow((sla_max_-sla_min_),2.)) ;
	b_ = 2. * sla_max_ / (Mypow((sla_min_-sla_max_),2.)) ;
	c_ = (Mypow(sla_min_,3.) + 2.*sla_min_*Mypow(sla_max_,2.) - 3.*Mypow(sla_min_,2.)*sla_max_) / (Mypow((sla_min_ - sla_max_),3.)) ;

	// coefficients for an upwardly open parabola, with coordinates (sla_min_,0) and (sla_max_,1), with (sla_min_,0) being the vertex of the parabola
/*	a = 1. / (Mypow(sla_max_,2.) - 2*sla_min_*sla_max_ + Mypow(sla_min_,2.)) ;
	b = (-2.*sla_min_)/(Mypow(sla_max_,2.) - 2*sla_min_*sla_max_ + Mypow(sla_min_,2.)) ;
	c = (Mypow(sla_min_,2.))/(Mypow(sla_max_,2.) - 2*sla_min_*sla_max_ + Mypow(sla_min_,2.)) ;  */
	

}












