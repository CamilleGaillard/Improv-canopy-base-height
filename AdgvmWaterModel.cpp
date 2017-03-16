
// #include <fstream>
// #include <iostream>
#include "AdgvmWaterModel.h"
#include "GridCellClassConstants.h"
#include "MyMath.h"
#include "LeafGlobals.h"
#include <iostream>
#include <iomanip>

using namespace std;


// const double *THICK;
// const double *DEPTH;


// Update water content in soil layers
// Input  : drain, theta, soil[k].tfc
// Output : BIn

// ----------------------------------------------------------------------------------------------------------

void BucketIn( GridCellSoilData *soil, GridCellClimateData *clim,GridCellPlantData *plantdata )
{
	double soilin;			   	//rain input (m/day)
	double drain = clim->pre;
	double sumSD;
	double *ThetaTmp = new double[N_SOIL_LAYERS];
	double *SD       = new double[N_SOIL_LAYERS];	//this is how much water can take up by soil SoilDry (m)
// 	const double DataFC[N_SOIL_LAYERS] = { 0.47, 0.5, 0.46, 0.5, 0.5, 0.5, 0.5 }; //STEVE
//	double DataFC[N_SOIL_LAYERS];
//	for ( int k=0; k<N_SOIL_LAYERS; k++ ) DataFC[k] = soil[k].tfc;

	double DataSWC[N_SOIL_LAYERS];
	for ( int k=0; k<N_SOIL_LAYERS; k++ ) DataSWC[k] = soil[k].sat_water_cont;
	
	double frac_cov = plantdata->ctr;
	if ( frac_cov > 1) frac_cov = 1;
	double biome = 0.02;				// Gerten et al. 2004 Table 1
	
	if ( drain > 0 )
	{
	    double can_stor = drain * plantdata->lai * biome * frac_cov;
	    
	      if ( drain > can_stor )
	      {
		drain -= can_stor;
	      }
	      else
	      {
		drain=0;
	      }
	      
	  can_stor = 0;
	}
	
	sumSD    = 0;
	soilin   = drain*1e-3;
	
//	cout << "Water status: " << DataSWC[0] <<setw(14)<< drain <<setw(14)<< soil[0].smo << setw(14)<<  soil[1].smo << setw(14)<< soil[2].smo << setw(14)<< soil[3].smo << setw(14)<< soil[4].smo << setw(14)<< soil[5].smo << setw(14)<< soil[6].smo <<endl; 

	for ( int i=0; i<N_SOIL_LAYERS; i++ )
		{
			if(soil[i].smo != soil[i].smo) {cout << "soil moisture NaN for layer " << i << " in BucketIn at pos 1!  " << endl; exit(1); }
			ThetaTmp[i] = soil[i].smo;			
		}
//------------------------------------------------------------------------------------------------
// original tipping bucket soilin version
// 	for ( int i=0; i<N_SOIL_LAYERS; i++) cout << " Theta_top_"<< i << "  " << ThetaTmp[i] << "  soilin  " << soilin << endl;  	
		

	for ( int i=0; i<N_SOIL_LAYERS; i++ )
	{
//		SD[i] = (soil[i].tfc-ThetaTmp[i])*THICK[i];
		SD[i] = (DataSWC[i]-ThetaTmp[i])*THICK[i]; //STEVE
		sumSD = sumSD+SD[i];			// how much water can be taken up by all soil layers
	}

	if ( soilin > sumSD )				// more rain than can be taken up
	{
		for ( int i=0; i<N_SOIL_LAYERS; i++ )
// 			ThetaTmp[i] = soil[k].tfc;	// set water contents to fc
			ThetaTmp[i] = DataSWC[i];	// STEVE
		
		soilin = 0;
	}
	
	if ( soilin > 0 )
	{
		for ( int i=0; i<N_SOIL_LAYERS; i++ )
		{
			if ( soilin>SD[i] )
			{
// 				ThetaTmp[i] = soil[i].tfc;
				ThetaTmp[i] = DataSWC[i]; //STEVE
				soilin -= SD[i];
				soilin  = MyMax( 0, soilin );
			}
			else
			{
				ThetaTmp[i] = ThetaTmp[i]+soilin/THICK[i];
				soilin = 0;
			}
		}
	}
// 	for (int i=0; i<N_SOIL_LAYERS; i++) cout << " Theta_bottom_"<< i << "  " << ThetaTmp[i] << endl;  	
// end original tipping bucket soilin version
//----------------------------------------------------------------------------------------------------
 	int time_step=24;
//	int time_step=1; // Liam 04.04.2016
	double perc;
//--------------------------------------------------------------------------------------------------------
// percolation calculated at sub-daily timesteps -- avoids overestimation of percolation
// ksat(input data) in meters per hour
// can replace ksat with w_cond to quickly try new parameters	
//--------------------------------------------------------------------------------------------------------	
	for ( int j=0; j<time_step; j++ )
	{
		for ( int i=1; i<N_SOIL_LAYERS; i++ )
		{
					// exponent(r)---values from RhoneAGG data set (www.iges.org/gswp/input.html)

// 		  perc = soil[i-1].ksat*Mypow((ThetaTmp[i-1]/DataSWC[i-1]), ((2*R[(int) soil[i].soil_texture])+2));	//FLAG3	// Percolation in m/day
		  perc = (soil[i-1].ksat)*Mypow((ThetaTmp[i-1]/DataSWC[i-1]), ((2*R[(int) soil[i].soil_texture])+2));	//FLAG3	// I have reduced the number of time steps for efficiency
		  perc = MyMin(ThetaTmp[i-1],perc / THICK[i]);	// FLAG4: temporary workaround for cased when drainage exceeds water content in soil layer!!!
//		  perc = perc / THICK[i];		// LIAM 02.03.2015
		  
		  if ( perc>0.1) perc=0.1;	//NOTE:BUGFIX(Liam) - if percolation exceeds the amount of water in a bucket then it produces Nans. This is problem with discretisation of continuous process			  
		if ( ThetaTmp[i]<perc ) perc=0.;
		  
		  ThetaTmp[i-1] -= perc;						
		  ThetaTmp[i] += (perc*(THICK[i-1]/THICK[i]));		// N.B. add (Fc[i-1]/Fc[i]) when soil layers are diff. 
		  if ( ThetaTmp[i]>DataSWC[i]) ThetaTmp[i]=DataSWC[i];	//Liam: if water exceeds holding capacity set to holding capacity-4 when we have different soil layers!	
		  
	          if (i == (N_SOIL_LAYERS-1) )
	          {
        	      	perc = (soil[i].ksat)*Mypow((ThetaTmp[i]/DataSWC[i]), ((2*R[(int) soil[i].soil_texture])+2));		// Percolation in m/day --ksat m/hr at the moment
              		perc = perc/THICK[i];
              		if ( perc>0.1) perc=0.1;
              		if ( ThetaTmp[i]<perc ) perc=0.;
              		ThetaTmp[i] -= perc;
              		if ( ThetaTmp[i]>DataSWC[i]) ThetaTmp[i]=DataSWC[i];
        	  }
	  
		  if(ThetaTmp[i] != ThetaTmp[i]) {cout << "ThetaTmp got NaN for layer " << i << ": " << perc << "  " << soil[i-1].ksat << "  " << soil[i-1].ksat << "  " << ThetaTmp[i-1] << "  " << DataSWC[i-1] << endl; exit(1); }
		  perc = 0.0;		  
		}		
	}
     
	for ( int i=0; i<N_SOIL_LAYERS; i++ ) 
	{
		soil[i].smo = ThetaTmp[i];
		if(soil[i].smo != soil[i].smo) {cout << "soil moisture NaN for layer " << i << " in BucketIn pos 2!" << endl; exit(1); }
	}	

	delete [] ThetaTmp;
	delete [] SD;
	return;
}

//--------------------------------------------------------------------------------------------------------	

void BucketOut( GridCellSoilData *soil, GridCellPlantData *plantdata )
{

	double soilout;
	double sumSL;
	double *ThetaTmp = new double[N_SOIL_LAYERS];
	double *SL       = new double[N_SOIL_LAYERS];
	double DataWP[N_SOIL_LAYERS];
	//int rank_SL[N_SOIL_LAYERS];
	double s_out_weight[N_SOIL_LAYERS];		// soilout for each soil layer
	double s_out_remainder[N_SOIL_LAYERS];		// difference between soilout for each soil layer and moisture that is removed
	for ( int i=0; i<N_SOIL_LAYERS; i++) s_out_remainder[i]=0.;		// difference between soilout for each soil layer and moisture that is removed

// 	double sum_s_out_remainder;		// total difference between soilout and removed soil moisture - ?set to zero at every time-step or keep track of the cumulative deficit? 
	double sum_s_out_remainder = 0.;
		
//	double water_potential_wilt = (sat_water_cont)*Mypow(phi_s/(12*100), (1/r)); // wilt point from matric potential method		// NEW LIAM		FLAG: ADJUSTED THIS BELOW TO MAKE IT BE BY SOIL LAYER!
	double water_potential_wilt[N_SOIL_LAYERS];
	
//	for ( int k=0; k<N_SOIL_LAYERS; k++ ) DataWP[k] = soil[k].twp;		//All plant WP are for sunflowers or cabbagges or something similar!
	
	for ( int i=0; i<N_SOIL_LAYERS; i++ )
	{
// 		water_potential_wilt[i]	= soil[i].sat_water_cont * Mypow(PHI_S[(int)soil[i].soil_texture]/(12*100), (1/R[(int)soil[i].soil_texture]));		// NEW LIAM, ADJUSTED TO BE BY SOIL LAYER	
//		water_potential_wilt[i]	= soil[i].sat_water_cont * Mypow(PHI_S[(int)soil[i].soil_texture]/(3*100), (1/R[(int)soil[i].soil_texture])) + RESID_CONTENT[(int)soil[i].soil_texture];		// NEW LIAM, ADJUSTED TO BE BY SOIL LAYER	
		water_potential_wilt[i]	= (soil[i].sat_water_cont - RESID_CONTENT[(int)soil[i].soil_texture]) * Mypow(PHI_S[(int)soil[i].soil_texture]/(3*100), (1/R[(int)soil[i].soil_texture])) + RESID_CONTENT[(int)soil[i].soil_texture];		// NEW LIAM, ADJUSTED TO BE BY SOIL LAYER	
		
		DataWP[i] = water_potential_wilt[i]; // Minimum WP (-3 MPa) from Plantclass.h based on p50		// NEW LIAM		
	  	ThetaTmp[i] = soil[i].smo;
	  	if(soil[i].smo != soil[i].smo) {cout << "soil moisture NaN for layer " << i << " in BucketOut pos1!" << endl; exit(1);}
	}
	
	sumSL   = 0.;
	soilout = plantdata->etr*1e-3; 							//convert mm/day into m/day => updated in GridCellClass on a daily basis, based on sum of Evapotranspiration calculated in PlantPopClass
	
	for ( int i=0; i<N_SOIL_LAYERS; i++ )
	{
		SL[i] = (ThetaTmp[i]-DataWP[i])*THICK[i]; 			//how much water can be lost from soil (m) STEVE
		if (SL[i]>0) sumSL += SL[i];
		s_out_weight[i] = soilout*soil[i].root_bm;					//calc water to be taken from layers based on total root biomass in soil layer; 
	}
	
	if ( soilout>0 )
	{
		for ( int i=0; i<N_SOIL_LAYERS; i++ )
		{
		  if (SL[i]>0)
		  {
		      if ( s_out_weight[i]>SL[i] )
			{
				ThetaTmp[i] = DataWP[i]; 	//STEVE
				s_out_weight[i] -= SL[i];
				s_out_weight[i]  = MyMax( 0, s_out_weight[i] );
				s_out_remainder[i] = s_out_weight[i];				// store the difference between what should and can be taken out and re-iterate below
			}
			else
			{	//cout << "THETA-ELSE 1: " << ThetaTmp[i] <<setw(14)<< s_out_weight[i] <<setw(14)<< s_out_weight[i]/THICK[i]<< endl;
				ThetaTmp[i] = ThetaTmp[i]-s_out_weight[i]/THICK[i];	//cout << "THETA-ELSE 2: " << ThetaTmp[i] << endl;
				s_out_weight[i] = 0.;
			}
		  sum_s_out_remainder += s_out_remainder[i];
		  }
		}
	}	
;
//-------------------------------------------------------------------------------------------
// soilout original
// 	if ( soilout>0 )
// 	{
// 		for ( i=0; i<N_SOIL_LAYERS; i++ )
// 		{
// 			if ( soilout>SL[i] )
// 			{
// //				ThetaTmp[i] =  soil[i].twp;
// 				ThetaTmp[i] = DataWP[i]; //STEVE
// 				soilout -= SL[i];
// 				soilout  = MyMax( 0, soilout );
// 			}
// 			else
// 			{
// 				ThetaTmp[i] = ThetaTmp[i]-soilout/THICK[i];
// 				soilout = 0;
// 			}
// 		}
// 	}
// soilout orininal end	
//-------------------------------------------------------------------------------------------
// NOTE: The intention here was to make sure transpired water is taken from the soil somehow
// 	 if the previous soil layers had no water left 
//	 Commented out as it shouldn't make a difference - will make a difference if we try to calculate ground water recharge

// 	if ( sum_s_out_remainder>0 )		// Taking this out for now as it allows water in deep soil to be stolen by shallow rooting plants
// 	{
// 		for ( int i=0; i<N_SOIL_LAYERS; i++ )
// 		{
// 		  SL[i] = (ThetaTmp[i]-DataWP[i])*THICK[i]; 			//STEVE
// 		  if ( SL[i]>0 ) sumSL += SL[i];
//  		    if ( SL[i]>0 )
// 		    {
// 
// 		      if ( sum_s_out_remainder>SL[i] )
// 			  {
// 				ThetaTmp[i] = DataWP[i];
// 				sum_s_out_remainder -= SL[i];
// 				sum_s_out_remainder  = MyMax( 0, sum_s_out_remainder );
// 				
// 			}
// 			else
// 			{
// 				ThetaTmp[i] = ThetaTmp[i]-sum_s_out_remainder/THICK[i];
// 				sum_s_out_remainder = 0.;
// 			}
// 		    }
// 		}
// 	}

	for ( int i=0; i<N_SOIL_LAYERS; i++ ) 
	{
		if(soil[i].smo != soil[i].smo) {cout << "soil moisture NaN for layer " << i << " in BucketOut pos 2!" << endl; exit(1); }
		soil[i].smo = ThetaTmp[i];
	}	
	
//	cout << "ThetaTmp bottom BucketOut: " << ThetaTmp[0] << endl;
	
	delete [] ThetaTmp;
	delete [] SL;	
	return;
}

//---------------------------------------------------------------------------------------------------------------------------------------------

/* double getET( double rad, double gs, double tmp, double apr, double hum, double wnd, double vhe, double car )		// MP: THIS FUNCTION SEEMS TO BE LEGACY, NOT USED ANYWHERE; EVAPOTR. CURRENTLY CALCULATED IN PlantClass.cpp
{
	double Zd   = ZD_CONST*vhe;
	double Z0   = Z0_CONST*vhe;
	double log_tmp = 1./log((REF_HEIGHT_Z-Zd)/Z0);
	
	double ga = KARMAN_CONST_QUAD*wnd*log_tmp*log_tmp;
	double Rn = rad;
	
// 	cout << setw(14) << ga << setw(14) << gn << endl;
	
	double s = 2504.0*exp((17.27*tmp)/(tmp+237.2))/pow(tmp+237.3,2);
	double cp = 1.013e-3;
	double gama = apr*0.001*0.000665;
	double rho = (apr*1e-3)/(0.287*(1.01*(tmp+273.3)));
	double eS = 0.6108*exp((17.27*tmp)/(tmp+237.3));
	double eA = hum/100.*eS;
	double VPD = eS-eA;
	
	double et_tmp = 0;
	
// 	cout << "GETET" << setw(14) << rad << setw(14) << gs << setw(14) << tmp << setw(14) << apr << setw(14) << hum << setw(14) << wnd << setw(14) << vhe << setw(14) << Rn << setw(14) << ga << setw(14) << s << setw(14) << cp << setw(14) << gama << setw(14) << rho << setw(14) <<  VPD << setw(14) << (s*Rn+86400.*rho*cp*VPD*ga)/(LAMBDA*(s+gama*(1.0+ga/gs))) << endl;
	
	if(gs>0.0001)
	{
		et_tmp = (s*Rn+86400.*rho*cp*VPD*ga)/(LAMBDA*(s+gama*(1.0+ga/gs)))*car/10000.;
// 		et_tmp = 0.;
	}
	
	double leaf_helper =  log( (et_tmp*(1./gs+1./gs) + hum/100.*eS)/0.6108 );
	
	double leaf_temp   = 237.3*leaf_helper/(17.27-leaf_helper);
	
// 	cout << setw(14) << tmp << setw(14) << leaf_temp << setw(14) << et_tmp << setw(14) << et_tmp*(2./gs) << setw(14) << hum/100.*eS << setw(14) << et_tmp*(2./gs) + hum/100.*eS << setw(14) << hum << setw(14) << eS << setw(14) << gs << endl;
	
	return et_tmp;
}
*/




















