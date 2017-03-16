
#include <cstdlib>
#include <cstring>
#include "adgvm.h"
#include <iostream>
#include <string>
#include "NcInputClasses.h"
#include "GridCellClass.h"
#include "GridCellClassConstants.h"
#include "PlantPopClassConstants.h"
#include "NcOutputClass.h"
using namespace std;
															 

extern std::map<std::string,std::string> netcdf_longname, netcdf_units;
extern std::map<std::string,nc_type> netcdf_datatype;


static void handle_error(int status, int ident)
{
	cout<<"Netcdf output error, status = "<< nc_strerror(status) << ", ident = " << ident <<endl;
}

//=================================================================================
	
NcOutputClass::NcOutputClass(const char* filename, unsigned int ctype, MyInData *IData)			// constructor
{


	libconfig::Setting& pop_outvars = cfg.lookup("output.outvars");
	pop_written.clear();
	pop_written["time"]=1;
	for(int i=0;i<pop_outvars.getLength();i++) 
		{
		const string s = pop_outvars[i];
		pop_written[s]=1;
		}
	
	contenttype=ctype;

	
	if (contenttype == POP_TYPE) BUFF_SIZEx=BUFF_SIZE;  //for the POP_TYPE we use a buffer of "BUFF_SIZE" time stamps in one write
	else BUFF_SIZEx = 1;


	
	latcount=IData->yread;
	loncount=IData->xread;
	recnum = new size_t[latcount*loncount];
//	valid_data = new int[latcount*loncount];
	for(int i=0;i<latcount*loncount;i++) {
		recnum[i]=0;
//		valid_data=0;
		}
	// create the netcdf output file for the output data 
	 
	status = nc_create(filename,NC_NETCDF4, &ncid);
	if (status != NC_NOERR) handle_error(status, 1);

	// define the dimensions. The time dimension for timestamps/iterations will be of unlimited length

	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &timestamp_dimid); //always present
	if (status != NC_NOERR) handle_error(status, 2);
	status = nc_def_var(ncid, "time", NC_FLOAT, 1, &timestamp_dimid, &varid);
	data_varid["time"] = varid;
	if (status != NC_NOERR) handle_error(status, 3);	
	
	status = nc_def_dim(ncid, "lat", latcount, &lat_dimid); 	
	if (status != NC_NOERR) handle_error(status, 4);
	status = nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid);		// 1: the dimension of latitude
	if (status != NC_NOERR) handle_error(status, 5);

	
	status = nc_def_dim(ncid, "lon", loncount, &lon_dimid);
	if (status != NC_NOERR) handle_error(status, 6);
	status = nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid);
	if (status != NC_NOERR) handle_error(status, 7);	

	
	if (contenttype == POP_TYPE)
	{
		//only in popfile
	
		status = nc_def_dim(ncid, "layer", N_SOIL_LAYERS, &layer_dimid);
		if (status != NC_NOERR) handle_error(status, 8);	
		status = nc_def_var(ncid, "layer", NC_INT, 1, &layer_dimid, &layer_varid);
		if (status != NC_NOERR) handle_error(status, 9);	
		status = nc_def_dim(ncid, "vtype", N_VTYPE, &vtype_dimid);
		if (status != NC_NOERR) handle_error(status, 8);	
		status = nc_def_var(ncid, "vtype", NC_INT, 1, &vtype_dimid, &vtype_varid);
		if (status != NC_NOERR) handle_error(status, 9);	

	}
	if (contenttype == TRAIT_TYPE)
	{
		//only in traitfile
		status = nc_def_dim(ncid, "individual", MAX_POP_SIZE, &ind_dimid); 
		if (status != NC_NOERR) handle_error(status, 10);
		status = nc_def_var(ncid, "individual", NC_INT, 1, &ind_dimid, &ind_varid);
		if (status != NC_NOERR) handle_error(status, 11);	
	}

	// assign attributes to the dimension variables where applicable
	
	
	status = nc_put_att_text(ncid,time_varid, "long_name", strlen("time"), "time");
	if (status != NC_NOERR) handle_error(status, 110);
	
	status = nc_put_att_text(ncid,time_varid, "units", strlen("years since 0000-00-00 "), "years since 0000-00-00 ");
	if (status != NC_NOERR) handle_error(status, 111);
	
	status = nc_put_att_text(ncid,time_varid, "calendar", strlen("noleap"), "noleap");
	if (status != NC_NOERR) handle_error(status, 112);	
	
	
	
	status = nc_put_att_text(ncid, lat_varid, "long_name", strlen("latitude"), "latitude");
	if (status != NC_NOERR) handle_error(status, 12);		
		
	status = nc_put_att_text(ncid, lat_varid, "units", strlen("degrees north"), "degrees north");
	if (status != NC_NOERR) handle_error(status, 13);	

	double latrange[2]={IData->minlat_,IData->maxlat_};
	status = nc_put_att_double(ncid, lat_varid, "actual_range", NC_DOUBLE, 2, latrange);  
	if (status != NC_NOERR) handle_error(status, 14);	
	
	status = nc_put_att_text(ncid, lat_varid, "long_name", strlen("longitude"), "longitude");
	if (status != NC_NOERR) handle_error(status, 15);			
	
	status = nc_put_att_text(ncid, lon_varid, "units", strlen("degrees east"), "degrees east");
	if (status != NC_NOERR) handle_error(status,16);	

	double lonrange[2]={IData->minlon_,IData->maxlon_};
	status = nc_put_att_double(ncid, lon_varid, "actual_range", NC_DOUBLE, 2, lonrange);  
	if (status != NC_NOERR) handle_error(status, 17);	
	
	if (contenttype == POP_TYPE)
	{
		status = nc_put_att_text(ncid, layer_varid, "long_name", strlen("soil layer number"), "soil layer number");
		if (status != NC_NOERR) handle_error(status, 18);	
		int l_range[2] = {1, N_SOIL_LAYERS};
		status = nc_put_att_int(ncid, layer_varid, "actual_range", NC_INT, 2, l_range);	
		if (status != NC_NOERR) handle_error(status, 19);	

		status = nc_put_att_text(ncid, vtype_varid, "long_name", strlen("vegetation type (ALL,TREE,GRASS)"), "vegetation type (ALL,TREE,GRASS)");
		if (status != NC_NOERR) handle_error(status, 18);	
		int l_range2[2] = {1, N_VTYPE};
		status = nc_put_att_int(ncid, vtype_varid, "actual_range", NC_INT, 2, l_range2);	
		if (status != NC_NOERR) handle_error(status, 19);	

	}
	
	if (contenttype == TRAIT_TYPE)
	{

		status = nc_put_att_text(ncid, ind_varid, "long_name", strlen("individual number"), "individual number");
		if (status != NC_NOERR) handle_error(status, 20);	
	
		int l_range[2] = {1, MAX_POP_SIZE};
		status = nc_put_att_int(ncid, ind_varid, "actual_range", NC_INT, 2,  l_range);	
		if (status != NC_NOERR) handle_error(status, 21);	
	}


	
	// write the dimension variable data, i.e., put the latitudes and longitudes of our data grid into the nc-file.
	
	status = nc_put_var_float(ncid, lat_varid, IData->climlats.storage);
	if (status != NC_NOERR) handle_error(status, 38);	
	
	status = nc_put_var_float(ncid, lon_varid, IData->climlons.storage);
	if (status != NC_NOERR) handle_error(status, 39);			
	
	// also write soil layer and individual vector out to ncfiles
	
	if (contenttype == POP_TYPE)
	{
		int layer_vector[N_SOIL_LAYERS];	
		for (int i=0; i<N_SOIL_LAYERS; i++) layer_vector[i] = i + 1;
		status = nc_put_var_int(ncid, layer_varid, layer_vector);
		if (status != NC_NOERR) handle_error(status, 40);	
		int vtype_vector[N_VTYPE];	
		for (int i=0; i<N_VTYPE; i++) vtype_vector[i] = i;
		status = nc_put_var_int(ncid, vtype_varid, vtype_vector);
		if (status != NC_NOERR) handle_error(status, 40);	

	}	
	
	if (contenttype == TRAIT_TYPE)
	{
		int ind_vector[MAX_POP_SIZE];
		for (int i=0; i<MAX_POP_SIZE;i++) ind_vector[i] = i +1;
		status = nc_put_var_int(ncid, ind_varid, ind_vector);
		if (status != NC_NOERR) handle_error(status, 41);	
	}			




	firstrun=true;
}
//-------------------------------------------------------------------------------------------------------- end of constructor

void NcOutputClass::WriteHeader(std::map <std::string,vector<float>> data)
/*
WriteHeader is responsible for writing out the Variable specific header to the Netcdf file. The data originates from the global variables
netcdf_longname, netcdf_units and netcdf_datatype These Variables are dynamically filled with information in PlantPopClass.cpp
*/
{	
	firstrun=false;
	status = nc_redef(ncid);  //put nc file back in define mode in order to define new variables!
	if (status != NC_NOERR) handle_error(status, 36);	
	if (contenttype == POP_TYPE) cout <<endl<< "NOTICE: Following variables not written to pop.nc file: ";
	if (contenttype == TRAIT_TYPE) cout <<endl<< "NOTICE: Following variables not written to trait.nc file: ";

// now define all the data variables:
 		
		for  (std::map<std::string,vector<float>>::iterator it=data.begin(); it!=data.end(); ++it)
		 {
			if ((it->first == "info") || (it->first == "time")) continue; //info/time contains internal data that we do not want to put in header so skip it
			
			//create the corresponding variable in the netcdf outputfile
			
			if ((pop_written.count(it->first)==1)) 
			{ // cout<<"Creating Netcdf variable: "<<pop_outvarname[i-1]<<endl;
				size_t DimLengths[5]={1,1,1,1,1};
				int dimids[5]={timestamp_dimid,1,1,1,1};
				unsigned int dims=1;
				DimLengths[0]=BUFF_SIZEx;
				switch (it->second.size()/BUFF_SIZEx)
				{
				case 1: //OT_SCALAR
				//create the mll_array for the scalar type
				DimLengths[1]=1;	
				DimLengths[2]=1;
				dimids[1]=lat_dimid;
				dimids[2]=lon_dimid;				
				dims=3;
				break;
				case N_SOIL_LAYERS:// OT_SOIL
			// create an empty MLL_array for vector population output variable, e.g., smo, for BUFF_SIZE, N_SOIL_LAYERS, latcount, loncount)
				DimLengths[1]=N_SOIL_LAYERS;
				DimLengths[2]=1;
				DimLengths[3]=1;	
				dimids[1]=layer_dimid;
				dimids[2]=lat_dimid;
				dimids[3]=lon_dimid;				
				dims=4;	
			 	break;
			        case N_VTYPE: //OT_VTYPE
				//create the mll_array for the scalar type
				DimLengths[1]=N_VTYPE;
				DimLengths[2]=1;
				DimLengths[3]=1;	
				dimids[1]=vtype_dimid;
				dimids[2]=lat_dimid;
				dimids[3]=lon_dimid;				
				dims=4;	
				break;
			        case MAX_POP_SIZE: //OT_VTYPE
				//create the mll_array for the scalar type
				DimLengths[1]=MAX_POP_SIZE;
				DimLengths[2]=1;
				DimLengths[3]=1;	
				dimids[1]=ind_dimid;
				dimids[2]=lat_dimid;
				dimids[3]=lon_dimid;				
				dims=4;	
				break;
				default:
				cout<<"warning undefined Netcdf outputtype dimension OT"<<endl;
				break;
				}

	
				status = nc_def_var(ncid, it->first.c_str(),netcdf_datatype[it->first], dims, dimids, &varid);		
				data_varid[it->first]=varid;
				if (status != NC_NOERR) handle_error(status, 22);
				status = nc_def_var_chunking(ncid, varid,NC_CHUNKED,DimLengths);		
				if (status != NC_NOERR) handle_error(status, 22);
			
				// PUT DOWN THE LONG_NAME AND UNIT OF EACH VARIABLE, AS WELL AS _FillValue AND missing_value
				float missing_=FMISSING;
				status = nc_put_att_text(ncid, varid, "long_name", netcdf_longname[it->first].size(), netcdf_longname[it->first].c_str());
				if (status != NC_NOERR) handle_error(status, 23);
			
				status = nc_put_att_text(ncid,varid, "units", netcdf_units[it->first].size(), netcdf_units[it->first].c_str());
				if (status != NC_NOERR) handle_error(status, 24);
			
				status = nc_put_att_float(ncid,varid, "_FillValue", netcdf_datatype[it->first], 1, &missing_);		
				if (status != NC_NOERR) handle_error(status, 25);
			
				status = nc_put_att_float(ncid, varid, "missing_value", netcdf_datatype[it->first], 1, &missing_);		
				if (status != NC_NOERR) handle_error(status, 26);								

			 }
			 else {cout<<"\""<<it->first <<"\",";}
			
		}
		
	// end define mode
	status = nc_enddef(ncid);
	if (status != NC_NOERR) handle_error(status, 37);	
	cout<<endl;
}	


int NcOutputClass::WriteBuffer(std::map <std::string,vector<float>> data)
{
if (firstrun) WriteHeader(data);
// get the lat/lon values from the info fields of the "data" map!
unsigned int lati = MyRound((data["info"])[0]);
unsigned int loni = MyRound((data["info"])[1]);

 for (std::map<std::string,vector<float>>::iterator it=data.begin(); it!=data.end(); ++it)
   	{	//cout<<"writebuffer b"<<endl;
		if (it->first == "info") continue; //info contains internal data that we do not want to write out, so skip it
		// write out the actual data; do this using a loop over all variables with the same structure, to avoid endless code duplication
		// the count and start settings tell netcdf to write BUFF_SIZE of data
		
		if (pop_written.count(it->first)==1) 
		{  // cout<<"Writing out: "<<pop_outvarname[i]<<endl;
				size_t start[5]={1,1,1,1,1};
				size_t count[5]={1,1,1,1,1};
				unsigned int dims=1;
				start[0]=recnum[lati*loncount+loni];
				count[0]=BUFF_SIZEx;

				if(it->second.size() < BUFF_SIZEx) {cout<<"ERROR undefined Netcdf data buffer size OT"<<endl; continue;} //this should never happen....
				if(it->second.size()==BUFF_SIZEx)
				{
				//for the scalar type
				start[1]=lati;	
				start[2]=loni;
				count[1]=1;	
				count[2]=1;
				dims=3;
				}
				else
				{
			// for vector population output variable, e.g., smo, for BUFF_SIZEx, N_SOIL_LAYERS, latcount, loncount)
				start[1]=0;
				start[2]=lati;
				start[3]=loni;	
				count[1]=it->second.size()/BUFF_SIZEx;
				count[2]=1;
				count[3]=1;	
				dims=4;	
				}	
				//This is where the actual writeout happens!
				if (it->first == "time") varid=time_varid;  //time varid is stored locally as it is different for different concurrent netcdffiles
				else varid = data_varid[it->first];

				status = nc_put_vara_float(ncid,varid, start, count,it->second.data());
//				cout<<"written variable " << it->first.c_str()<<" "<<count[1]<<" "<<it->second.data()[BUFF_SIZEx-1]<<" " <<it->second.data()[0]<<endl;
				if (status != NC_NOERR)  { handle_error(status, 42);}
		}
	

	}
recnum[lati*loncount+loni] += BUFF_SIZEx;	//increase recnum to prepare for next iteration	
	//---------------------
	

return int(recnum[lati*loncount+loni]);
}

//-----------------------------------------------------------------------------------------------------

int NcOutputClass::BufferFull()
{
int example = 0;
int i;
for(i=0;i<latcount*loncount;i++) if (recnum[i] !=recnum[0]) break;
//cout<<"bufferdebug:" <<i<<endl;
return (i == latcount*loncount);
}

//-----------------------------------------------------------------------------------------------------
	
NcOutputClass::~NcOutputClass()
{	
	// to close the popdata nc-file and traitdata nc-file at the end of the run, and to free dynamically allocated arrays and objects
	
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status, 46);	
		data_size=0;
}




