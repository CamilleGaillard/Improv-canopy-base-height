#include <iostream>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <iomanip>
#include "adgvm.h"

#include "MyMath.h"
#include "NcInputClasses.h"
#include "LeafGlobals.h"
#include "Radiation.h"
#include "PenmanMonteith.h"


/* To compile in C++ with objects:
g++ -c Ncread2.cpp -I/usr/local/include -L/usr/local/lib -lnetcdf
g++ -c test.cpp -I/usr/local/include -L/usr/local/lib -lnetcdf
g++ -o test test.o Ncread2.o  (links the compiled object files)
*/

/* These are the names of the data files we will read. */
//#define CLIMFILE "/home/sscheiter/NcInputFiles/climvars_10min.nc"	// normally such stuff should be in a joboptions-file rather than 
//#define SOILFILE "/home/sscheiter/NcInputFiles/WISE_soilvars_10min.nc"		// being compiled into the code

using namespace std;

//================================

static void handle_error(int status, int ident)
{
	cout<<"Netcdf input file error, status = "<< status << ", ident = " << ident <<endl;
}
//======================================================================


NcInfile::NcInfile(int ncid, int count1, double minlon, double maxlon, double minlat, double maxlat)	// constructor; ncid identifies the infile on which the
																																									// following actions are supposed to be performed
	{   
	
	
	
		// inquire the dimid and dimlength for longitude, and read in the latitude values of the entire infile
	
		status = nc_inq_dimid(ncid, "lon", &dimid);
		if(status != NC_NOERR) handle_error(status, 1); 
	
		status = nc_inq_dimlen(ncid, dimid, &xlen);
		if(status != NC_NOERR) handle_error(status, 2);	
	
		double *lons_tot;	
		lons_tot = new double[xlen];
	
		status = nc_get_var(ncid, dimid, lons_tot);
		if(status != NC_NOERR) handle_error(status, 3);

		// inquire the dimid and dimlength for longitude, and read in the longitude values of the entire infile
	
		status = nc_inq_dimid(ncid, "lat", &dimid);
		if(status != NC_NOERR) handle_error(status, 4); 
	
		status = nc_inq_dimlen(ncid,dimid, &ylen);
		if(status != NC_NOERR) handle_error(status, 5);
	
		double *lats_tot;
		lats_tot = new double[ylen];
		
		status = nc_get_var(ncid, dimid, lats_tot);
		if(status != NC_NOERR) handle_error(status, 6);
	
		//--------------------------------------------------------
				
		// check that specified coordinates are within the range of coordinates available from input file
	
		if((minlon < lons_tot[0]) || (maxlon > lons_tot[xlen-1]) || (minlat < lats_tot[0]) || (maxlat > lats_tot[ylen-1])) 
		{
			cout << "specified coordinate range outside the range of input file, setting range to maximum available range" << endl;
	        	minlon=lons_tot[0];
			maxlon=lons_tot[xlen-1];
			minlat=lats_tot[0];
			maxlat=lats_tot[ylen-1];
		}
	
		// establish relationship between lon and lat values and array indices
	
		double lon_a = double(xlen-1) / (lons_tot[xlen-1] - lons_tot[0]) ;	// slope for longitude	
		double lon_b = lons_tot[0] * (-1.) * lon_a ;									// intercept for longitude
	
		double lat_a = double(ylen-1) / (lats_tot[ylen-1] - lats_tot[0]) ;		// slope for latitude
		double lat_b = lats_tot[0] * (-1.) * lat_a ;										// intercept for latitude
	
		delete[] lons_tot;
		delete[] lats_tot;
	
		double lon_ix_1 = lon_a * minlon + lon_b ;
		double lon_ix_2 = lon_a * maxlon + lon_b ;
		double lat_ix_1 = lat_a * minlat + lat_b ;
		double lat_ix_2 = lat_a * maxlat + lat_b ; 

		// round calculated indexes to nearest integer to access the closest lon/lat coordinate 

		lon_ix_1 = round(lon_ix_1);
		lon_ix_2 = round(lon_ix_2);
		lat_ix_1 = round(lat_ix_1);
		lat_ix_2 = round(lat_ix_2);
	
		//define number of data points to be read (important to know the count value)
	
		xread = abs(lon_ix_2 - lon_ix_1)+1 ;
		yread = abs(lat_ix_2 - lat_ix_1)+1 ;
	
		//---------------------------------------------------------------------------------------
	
		// starting points and counts for reading the longitude and latitude values from nc-infile
			
		xsrt = new size_t[1];
		xsrt[0]=min(lon_ix_1,lon_ix_2);
		xcount = new size_t[1];
		xcount[0]=xread;

		ysrt = new size_t[1];
		ysrt[0]=min(lat_ix_1,lat_ix_2);
		ycount = new size_t[1];
		ycount[0]=yread;

		//------------------------------------------------

		// starting points and counts for reading 3D-variables from the input file
		
		srt = new size_t[3];	// starting point for reading the 3D-variables from the input file
		srt[0]=0;
		srt[1]=lat_ix_1;
		srt[2]=lon_ix_1;
		count = new size_t[3];	
		count[0]=count1; 
		count[1]=yread;
		count[2]=xread;
		
		// starting points and counts for reading 2D-variables from the input file
		
		srt2d = new size_t[2];
		srt2d[0] = lat_ix_1;
		srt2d[1] = lon_ix_1;
		count2d = new size_t[2];
		count2d[0] = yread;
		count2d[1] = xread;
		
	}	// end of constructor

NcInfile::~NcInfile()   			//destructor
	{
		delete srt;
		delete count;
		delete ysrt;
		delete ycount;
		delete xsrt;
		delete xcount;
		delete srt2d;
		delete count2d;
	}

//======================================================================

MLL_array::MLL_array()  //empty MLL_array constructor; to create MLL-type variables from template variables read from Nc-file
	{
		totsize=0;
		dimensions=0;
	}
MLL_array::MLL_array(unsigned int newdimensions, size_t *newdimsizes)	// constructor 1: creates an empty MLL_array with specified dimensions
{																														// newdimension: number of new dimensions to be created

    for(int i=0;i<5;i++) n[i]=1;		
	dimensions=newdimensions;
	for(int i=0;i<dimensions;i++) n[i]=newdimsizes[i];
	    totsize=1;
	    for(int i=0;i<5;i++) totsize=totsize*n[i];
	    storage = new float[totsize];	// allocate sufficient memory
	for(int i=0;i<totsize;i++) storage[i] = FMISSING;

}	
	
MLL_array::MLL_array(int ncid, const char *name,size_t *start, size_t *count)		// constructor 2; does actual netcdf-handling and allocates
																																// sufficient memory for the respective input variables
    {																															// for now up to 5D-variables possible
        for(int i=0;i<5;i++) n[i]=1;		
		
		int varid;
		int status;
			
		status = nc_inq_varid(ncid, name, &varid);					// link name of variable in nc-file with an ID
	   	if(status != NC_NOERR) handle_error(status, 7);
	        
	    status = nc_inq_varndims(ncid,varid,&dimensions);	// inquire the numer of dimensions of the input variable
	    if(status != NC_NOERR) handle_error(status, 8);
	        
	    for(int i=0;i<dimensions;i++) n[i]=count[i];
	    totsize=1;
	    for(int i=0;i<5;i++) totsize=totsize*n[i];
	        
		storage = new float[totsize];	// allocate sufficient memory
			
	    status = nc_get_vara_float(ncid, varid, start, count, storage);		//directly read into MLL_array structure
	    if(status != NC_NOERR) handle_error(status, 9);
	} 
		
MLL_array::~MLL_array()		// free memory again using the destructor
	{       
		delete [] storage;
	}
		
void MLL_array::EmptyFromTemplate(MLL_array * templatemllarray)
	{
		for(int i=0;i<5;i++) n[i]=1;
		dimensions=templatemllarray->dimensions;
		for(int i=0;i<dimensions;i++) n[i]=templatemllarray->n[i];
		totsize=1;
	    for(int i=0;i<5;i++) totsize=totsize*n[i];
		storage = new float[totsize];		
	}
		
float& MLL_array::operator () (unsigned int p0, unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4)		// for 5D-variables
	{	// p4=last dimension, varies fastest; p0=first dimension, varies slowest
		
		if ((dimensions==5) && (p0 < n[0]) && (p1 < n[1]) && (p2 < n[2]) && (p3 < n[3]) && (p4 < n[4])) return storage[(n[4]*n[3]*n[2]*n[1]*p0+n[4]*n[3]*n[2]*p1+n[4]*n[3]*p2+n[4]*p3+p4)]; 
		invalid = std::numeric_limits<float>::signaling_NaN();
		return invalid;
	}
		
float& MLL_array::operator () (unsigned int p0, unsigned int p1, unsigned int p2, unsigned int p3)		// for 4D-variables
	{
		if ((dimensions==4) && (p0 < n[0]) && (p1 < n[1]) && (p2 < n[2]) && (p3 < n[3])) return storage[(n[3]*n[2]*n[1]*p0+n[3]*n[2]*p1+n[3]*p2+p3)];
		invalid = std::numeric_limits<float>::signaling_NaN();
		return invalid;
	}	
		
float& MLL_array::operator () (unsigned int p0, unsigned int p1, unsigned int p2)		// for 3D-variables
	{
		if ((dimensions==3) && (p0 < n[0]) && (p1 < n[1]) && (p2 < n[2])) return storage[(n[2]*n[1]*p0+n[2]*p1+p2)];
		invalid = std::numeric_limits<float>::signaling_NaN();
		return invalid;
	}	
			
float& MLL_array::operator () (unsigned int p0, unsigned int p1)		// for 2D-variables
	{
		if ((dimensions==2) && (p0 < n[0]) && (p1 < n[1])) return storage[(n[1]*p0+p1)];
		invalid = std::numeric_limits<float>::signaling_NaN();
		return invalid;
	}	
			
float& MLL_array::operator () (unsigned int p0)		// for 1D-variables
	{
		if ((dimensions==1) && (p0 < n[0])) return storage[(p0)];
		invalid = std::numeric_limits<float>::signaling_NaN();
		return invalid;
	}	

//======================================================================

MyInData::MyInData(double minlon, double maxlon, double minlat, double maxlat)		// constructor
	{	
		string CLIMFILE = cfg.lookup("infile.climfile").c_str();
		string SOILFILE = cfg.lookup("infile.soilfile").c_str();

		minlon_ = minlon;
		maxlon_ = maxlon;
		minlat_ = minlat;
		maxlat_ = maxlat;
	
		//===================================================================
	
		// open the climate data input file
			
		status = nc_open(CLIMFILE.c_str(), NC_NOWRITE, &ncid);
		if(status != NC_NOERR) handle_error(status, 10);
	
		NcInfile climdata(ncid, 12, minlon, maxlon, minlat, maxlat)	;	// constructor of NcInfile class initializes xsrt, xcount,
																											// ysrt, ycount, srt, and count (12 = number of months)
		xread=climdata.xread;
		yread=climdata.yread;
		
		// get actual variables
	
		new (&climlons)  MLL_array(ncid,"lon",climdata.xsrt,climdata.xcount);	// 1D-array holding longitude coordinate values
		new (&climlats)  MLL_array(ncid,"lat",climdata.ysrt,climdata.ycount);		// 1D-array holding latitude coordinate values
			
		new (&tmp)   MLL_array(ncid,"tmp",climdata.srt,climdata.count);  
		new (&dtr)   MLL_array(ncid,"dtr",climdata.srt,climdata.count); 
		new (&sunp)  MLL_array(ncid,"sunp",climdata.srt,climdata.count);
		new (&reh)   MLL_array(ncid,"reh",climdata.srt,climdata.count); 
		new (&rd0)   MLL_array(ncid,"rd0",climdata.srt,climdata.count); 
		new (&wnd)   MLL_array(ncid,"wnd",climdata.srt,climdata.count); 
		new (&frost) MLL_array(ncid,"frs",climdata.srt,climdata.count);
		new (&pre)   MLL_array(ncid,"pre",climdata.srt,climdata.count);
		new (&preCV) MLL_array(ncid,"preCV",climdata.srt,climdata.count);
	
		new (&elv)   MLL_array(ncid, "elv", climdata.srt2d,climdata.count2d);

		// convert elevation from kilometers to meters
		for (int i=0; i<elv.totsize; i++) elv.storage[i]=elv.storage[i]*1000.;
	
	
		//--------------------------------------
	
		//close climate input file
	
		status = nc_close(ncid);
		if(status != NC_NOERR) handle_error(status, 11);		
	
		//===================================================================
		
		// number of soil layers used in the model
		
		soil_layers = 12;
	
		// open the soil data input file
	
		status = nc_open(SOILFILE.c_str(), NC_NOWRITE, &ncid);
		if(status != NC_NOERR) handle_error(status, 12);	
		
		NcInfile soildata(ncid, 2, minlon, maxlon, minlat, maxlat)	;	// constructor of NcInfile class initializes xsrt, xcount, 
																										// ysrt, ycount, srt, and count (2 = number of soil layers; [0] = topsoil, [1] = subsoil)
																										
		if (climdata.xread != soildata.xread || climdata.yread != soildata.yread)
			{
				cout << "WARNING: GRID MISMATCH! grid resolution of soil data file differs from climate data grid resolution! " << endl;
			}
																											
		// get actual variables and expand from 2 soil layers (toptoil, subsoil) to 12 model layers
	
		new (&soillons) MLL_array(ncid,"lon",soildata.xsrt,soildata.xcount);	// 1D-array holding longitude coordinate values
		new (&soillats)  MLL_array(ncid,"lat",soildata.ysrt,soildata.ycount);		// 1D-array holding latitude coordinate values
	    
	    // wilting point    
	    MLL_array theta_wp_in(ncid,"WP",soildata.srt,soildata.count);				// get wilting point from nc-infile
	    size_t ModelLayerDimsizes[3]={soil_layers,soildata.count[1],soildata.count[2]};		// create an empty array with 12 instead of 2 soil layers
		new (&theta_wp)   MLL_array(3,ModelLayerDimsizes);  
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (theta_wp_in(1,j,k)!=-9999)) theta_wp(i,j,k)=theta_wp_in(1,j,k)/100. ;	// assign topsoil values to model layer 1, and subsoil values (where available)
		      else theta_wp(i,j,k)=theta_wp_in(0,j,k)/100. ;       				// to layer 2-12; if no subsoil data are available, topsoil data are used for all 12 layers
		    }  
		
		// field capacity
		MLL_array theta_fc_in(ncid,"FC",soildata.srt,soildata.count);	
		new (&theta_fc)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (theta_fc_in(1,j,k)!=-9999)) theta_fc(i,j,k)=theta_fc_in(1,j,k)/100. ;
			  else theta_fc(i,j,k)=theta_fc_in(0,j,k)/100. ; 
			}
		
		// saturated soil hydraulic conductivity
		MLL_array k_sat_in(ncid,"Ksat",soildata.srt,soildata.count);	  
		new (&k_sat)  MLL_array(3,ModelLayerDimsizes);
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (k_sat_in(1,j,k)!=-9999)) k_sat(i,j,k)=k_sat_in(1,j,k);
			  else k_sat(i,j,k)=k_sat_in(0,j,k); 
			}		
		
		// soil texture
		MLL_array soil_texture_in(ncid,"texture",soildata.srt,soildata.count);	  
		new (&soil_texture)  MLL_array(3,ModelLayerDimsizes);
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (soil_texture_in(1,j,k)!=-9999)) soil_texture(i,j,k)=soil_texture_in(1,j,k);
			  else soil_texture(i,j,k)=soil_texture_in(0,j,k); 
			}		
				
		// water content at saturation
		MLL_array sat_water_cont_in(ncid,"Wsat",soildata.srt,soildata.count);			
		new (&sat_water_cont)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (sat_water_cont_in(1,j,k)!=-9999)) sat_water_cont(i,j,k)=sat_water_cont_in(1,j,k)/100. ;
			  else sat_water_cont(i,j,k)=sat_water_cont_in(0,j,k)/100. ; 
			}				
		
		// percolation exponent
		MLL_array perc_exponent_in(ncid,"exponent",soildata.srt,soildata.count);				
		new (&perc_exponent)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (perc_exponent_in(1,j,k)!=-9999)) perc_exponent(i,j,k)=perc_exponent_in(1,j,k);
			  else perc_exponent(i,j,k)=perc_exponent_in(0,j,k); 
			}						
		
		// soil bulk density
		MLL_array bulk_dens_in(ncid,"bulkdensity",soildata.srt,soildata.count);			
		new (&bulk_dens)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (bulk_dens_in(1,j,k)!=-9999)) bulk_dens(i,j,k)=bulk_dens_in(1,j,k);
			  else bulk_dens(i,j,k)=bulk_dens_in(0,j,k); 
			}		

		// soil nitrogen content
		MLL_array soil_N_in(ncid,"Ntot",soildata.srt,soildata.count);					
		new (&soil_N)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (soil_N_in(1,j,k)!=-9999)) soil_N(i,j,k)=soil_N_in(1,j,k);
			  else soil_N(i,j,k)=soil_N_in(0,j,k); 
			}		
		
		// soil organic carbon content
		MLL_array soil_C_in(ncid,"soil_C",soildata.srt,soildata.count);							
		new (&soil_C)  MLL_array(3,ModelLayerDimsizes); 
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (soil_C_in(1,j,k)!=-9999)) soil_C(i,j,k)=soil_C_in(1,j,k);
			  else soil_C(i,j,k)=soil_C_in(0,j,k); 
			}			
		
		// rooting depth
		MLL_array root_depth_in(ncid,"root_depth",soildata.srt,soildata.count);		
		new (&root_depth)  MLL_array(3,ModelLayerDimsizes); 	
		for (int i=0;i<ModelLayerDimsizes[0];i++)
		  for (int j=0;j<ModelLayerDimsizes[1];j++)
		    for (int k=0;k<ModelLayerDimsizes[2];k++)
		    {
		      if ((i>0) && (root_depth_in(1,j,k)!=-9999)) root_depth(i,j,k)=root_depth_in(1,j,k);
			  else root_depth(i,j,k)=root_depth_in(0,j,k); 
			}					
		//--------------------------------------
	
		//close soil input file
	
		status = nc_close(ncid);
		if(status != NC_NOERR) handle_error(status, 13);
		
		//--------------------------------------	

		// Bring input data into the format required by the model
		
		// create variables that have the same structure as the ones read in from the Nc-files
		
	tmp_min.EmptyFromTemplate(&tmp);		
       	tmp_max.EmptyFromTemplate(&tmp);
       	ralpha.EmptyFromTemplate(&tmp);
       	rbeta.EmptyFromTemplate(&tmp);  
       	pwet.EmptyFromTemplate(&tmp); 
       	tmp_day.EmptyFromTemplate(&tmp);
       	atm_press.EmptyFromTemplate(&elv); 
       	
       	for( int i=0;i<tmp.totsize;i++)
       	{
       		 tmp_min.storage[i] = tmp.storage[i] - 0.5*dtr.storage[i];
      		 tmp_max.storage[i] = tmp.storage[i] + 0.5*dtr.storage[i];
      		 if (pre.storage[i] ==0.)
      		 {
      		      ralpha.storage[i] = 1.;
      		      rbeta.storage[i] = 0.;
      		 }
      		 else
      		 {
      		      ralpha.storage[i] = 1./pow(preCV.storage[i]/100.,2);
      		      rbeta.storage[i] = pre.storage[i]/pow(preCV.storage[i]/100.,2);
      		}

			sunp.storage[i] = sunp.storage[i]/100.;
			tmp_day.storage[i] = tmp.storage[i] + (tmp_max.storage[i] - tmp.storage[i])/2.; // this from InDataReaderClass l. 89. Not sure this makes sense!!!			
		}
		
		// convert from number of days per month to fraction of month
		
		for (int i=0; i<12; i++)
			for (int j=0; j<yread; j++)
				for (int k=0; k<xread; k++)
					{
						pwet(i,j,k) = rd0(i,j,k) / NDAYS_IN_MONTH[i];	
						frost(i,j,k) = frost(i,j,k) / NDAYS_IN_MONTH[i];
					}
		
		// calculate altitude-dependent site-specific atmospheric pressure
		
		for (int i=0; i<elv.totsize; i++)
		{
			atm_press.storage[i] = 101.325*pow((293.0-0.0065*elv.storage[i])/293.0,5.26)*1000.;
		}
		
		// apply constraints to minimum and maximum soil-C content and soil-N content (l. 162-180 InDataReaderClass)
		
		for (int i=0; i<(soil_layers*elv.totsize); i++)
		 {
		 	if (soil_C.storage[i] < 0.01 && soil_C.storage[i] > -9999.) soil_C.storage[i] = 0.01;  // FLAG: CHECK ON SOIL-C DATA AT HWSD, LOOKS SPURIOUS WAY IT IS!
		 	if (soil_N.storage[i] < 0.01 && soil_N.storage[i] > -9999.) soil_N.storage[i] = 0.01;
		 	if (soil_C.storage[i] > 30000.) soil_C.storage[i] = 30000;
		 }	
		 
		 // define soil layer thickness (array of 12 values, identical for all grid cells)
		
		 for (int i=0; i<soil_layers; i++) thickness[i] = 0.3;
		 		 
		 // calculate depth of lower soil layer boundaries
		 
		 depth[0] = thickness[0];
		 for (int i=1; i<soil_layers; i++) depth[i] = depth[i-1]+thickness[i];
		 
		//-------------------------------------------------------
		
		//derived variables from calcatmospheric
		eS12.EmptyFromTemplate(&tmp);
		eA12.EmptyFromTemplate(&tmp);
		s12.EmptyFromTemplate(&tmp);
		gama12.EmptyFromTemplate(&tmp);
		VPD12.EmptyFromTemplate(&tmp);
		rho12.EmptyFromTemplate(&tmp);
		Q012.EmptyFromTemplate(&tmp);

        // used to be calcAtmospheric
		int midmonthday;
	
		// compute some atmospheric values of the current coordinates
		for ( int i=0; i<12; i++ )
		for ( int j=0; j<yread; j++ )
		for ( int k=0; k<xread; k++ )
		{
			midmonthday		= (int) floor((double) (i+1.)*30.42-15.21);

			eS12(i,j,k)		= FOAeS(tmp_max(i,j,k),tmp_min(i,j,k));		// average saturation vapor pressure, KPa; from PenmanMonteith.h
			// FOAe0 = saturation vapor pressure, reh = rel. hum in % => reh/100 in KPa
			// actual actual vapor pressure
			eA12(i,j,k)		= reh(i,j,k)/100.*eS12(i,j,k);  // ((FOAe0(tmax[i])+FOAe0(tmin[i]))/2.);
			// Radiation in MJ/m^2/day
			s12(i,j,k)			= FOAs(tmp(i,j,k)); 							// KPa/degC
			// psychochromatic constant 
			gama12(i,j,k)		= FOAgama(atm_press(j,k)/1000.);  				// KPa/degC
			// vapor pressure difference
			VPD12(i,j,k)		= FOAVPD(eA12(i,j,k),eS12(i,j,k)); 				// KPa
			// density of air
			rho12(i,j,k)		= FOArho(atm_press(j,k),tmp(i,j,k));				// g/m3 //takes atm_press in Pascals!
			// PAR
			Q012(i,j,k)		= GetPARRadiation(climlats(j),sunp(i,j,k), midmonthday);		// from Radiation.h
		}

	// end calcAtmospheric


	}





