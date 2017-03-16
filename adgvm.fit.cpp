#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include <stdlib.h>
#include "adgvm.h"
#include "MyMath.h"

using namespace std;


//----------------------------------------------------------------------------------
// global constants and variables

//#include "adgvm_interface.cpp"
#include "GridCellClassConstants.h"
#include "GridCellClass.h"

// use bucket model and evapotranspiration from adgvm
#ifdef USE_ADGVM_WATER_MODEL
#include "AdgvmWaterModel.h"
#endif
#ifdef NC_OUTPUT
#include "NcOutputClass.h"
#endif

const string IN_DATA_HOME  = "InputData/";
int nrank,nprocs;
libconfig::Config cfg;

//--------------------------

unsigned int cellindex(unsigned int cidx)
{
	return cidx*nprocs+nrank;
}
//--------------------------


double LIGHT_COMP_PAR = 0.05;
double LIGHT_COMP_PAR_HLP = 1.-LIGHT_COMP_PAR;
double MORT_BIOMASS[VTYPE_NUM] = {  0.05,      0.05 };


int main( int argc, char **argv )
{try {
	
	#ifdef WITH_MPI
	MPI::Init();	
	nrank=MPI::COMM_WORLD.Get_rank();
	nprocs=MPI::COMM_WORLD.Get_size();
	#else
	nrank=0;
	nprocs=1;
	#endif

	time_t start;		// to measure time to run the program
	time_t end;
	time (&start);	
	
	//cout << "RSEED     " << time(NULL) << endl;
	/*srand48(time(NULL));
	srand(time(NULL));
	{	
		double test=0;
		int i;
		for(i=0;i<10000000;i++) test +=RUnif();
		cout << "average of RUnif()="<< test/i<<endl;
	}*/
	
	int file_counter = 0;
	
	int    RUNID = -999;
	int    FIRE  =  0;       // no fire by default
	
	cout << endl;
	
//------------------------------------------------------------------------	

	// transform lat and lon strings from arguments to double values
	
	double minlon;
	double maxlon;
	double minlat;
	double maxlat;

	if(argc == 7)	// 4 coordinates, need to read data for a grid box
	{	
		RUNID = atoi( argv[1] );	// convert string to int-value using atoi
		FIRE  = atoi( argv[2] );
		minlon = atof(argv[3]) ;		
		maxlon = atof(argv[4]);
		minlat = atof(argv[5]);
		maxlat = atof(argv[6]);
		
		double tmp;
		if(maxlon < minlon) {tmp = maxlon; maxlon = minlon; minlon = tmp;}
		if(maxlat < minlat) {tmp = maxlat; maxlat = minlat; minlat = tmp;}
	}
	else if(argc == 5)	// 2 coordinates, need to read data for a point
	{
		RUNID = atoi( argv[1] );  
		FIRE  = atoi( argv[2] );  
		minlon = atof(argv[3]) ; 
		maxlon = minlon ;         
		minlat = atof(argv[4]);  
		maxlat = minlat ;         

	}	
	else if(argc ==3)
	{
		RUNID = atoi( argv[1] );
		FIRE  = atoi( argv[2] );
		minlon = 31.77;  // skukuza by default
		maxlon = 31.77;
		minlat = -24.39;
		maxlat = -24.39;
		
		cout << "WARNING: Wrong number of arguments, using standard arguments for Skukuza site." << endl;
	}

	else if(argc ==8 )
	{
		RUNID = atoi( argv[1] );  
		FIRE  = atoi( argv[2] );  
		minlon = atof(argv[3]) ; 
		maxlon = minlon ;         
		minlat = atof(argv[4]);  
		maxlat = minlat ;         

		LIGHT_COMP_PAR = atof(argv[5]);
		LIGHT_COMP_PAR_HLP = 1.-LIGHT_COMP_PAR;
		MORT_BIOMASS[TREE] = argv[6];     // mortality due to low height
		MORT_BIOMASS[GR_C4] = argv[7];     // mortality due to low height
		
		cout << "WARNING: Wrong number of arguments, using standard arguments for Skukuza site." << endl;
	}



	else
	{
		RUNID = 1;
		FIRE  = 0;
		minlon = 31.75;  // skukuza by default
		maxlon = 31.75;
		minlat = -24.25;
		maxlat = -24.25;
		
		cout << "WARNING: Wrong number of arguments, using standard arguments for Skukuza site." << endl;
	}
	
	int    runid;
	int    fire;
	
	
	runid = RUNID;
	fire  = FIRE;
	if ((fire<0) || (fire >1)) 
	{
	cout << "parameter (fire) is implausible: "<<fire<<endl;
	return 1;
	}
	if (nrank==0){
	cout << "N_SOIL_LAYERS      " << N_SOIL_LAYERS << endl;
	cout << "RUN_ID    " << runid << endl;
	cout << "FIRE      " << fire << endl;}
	
	//-------------------------------------------------------
	// get meta-information on the run from the run configuration file (aka, namelist/joboptions file)
	
	// Read the configuration file. If there is an error, report it and exit.
	try
		{
		cfg.readFile("runconfig.cfg");
		}
	catch(const libconfig::FileIOException &fioex)
		{
		std::cerr << "I/O error while reading runconfig.cfg file." << std::endl;
		return(EXIT_FAILURE);
		}
	catch(const libconfig::ParseException &pex)
		{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
		      << " - " << pex.getError() << std::endl;
		return(EXIT_FAILURE);
		}
	int YEARS_TO_RUN = 1;
	try
		{
		YEARS_TO_RUN = cfg.lookup("simulation.YEARS_TO_RUN");
		}
	catch(const libconfig::SettingNotFoundException &nfex)
		{
		cerr << "Error: Setting '" <<nfex.getPath()<< "' not found in configuration file." << endl;
		return(EXIT_FAILURE);
		}
	cout<<"Years to run: "<<YEARS_TO_RUN<<endl;
	
	//-------------------------------------------------------
	



	// get input data from netcdf file and calculate some atmospheric values (one year, for all specified grid cells)
	
	MyInData IData(minlon,maxlon,minlat,maxlat);
#ifdef NC_OUTPUT
 	NcOutputClass* trait_out;
	NcOutputClass* pop_out;
	if (nrank==0){
	char trait_filename_[250];
	char pop_filename_[250];
	sprintf( trait_filename_, "trait_%i_%i.nc", runid, fire );
	sprintf( pop_filename_, "pop_%i_%i.nc", runid, fire );
	trait_out= new NcOutputClass(trait_filename_,TRAIT_TYPE,&IData);
	pop_out= new NcOutputClass(pop_filename_,POP_TYPE,&IData);
	}
#endif
	
//=========================================================	
	
	unsigned int ncells = IData.yread*IData.xread;
	unsigned int ncells_local = ncells/nprocs+(nrank<(ncells%nprocs));		// number of grid cells run by each core. If the devision has a remainder then these get distributed until the last is distributed
	unsigned int year;							// each individual grid cell keeps track of the simulation year it is at; different grid cells may be at different years at a given time, depending on work load!
	unsigned int day;			// each individual grid cell keeps track of the simulation day it is at
	
	if (nrank==0) cout << "Running aDGVM2 for "<<ncells<<" gridcells on a ("<<IData.xread<<"," <<IData.yread << ") grid for "<<YEARS_TO_RUN<<" years" <<endl;	
	cout<<"Node " << nrank<<" running "<< ncells_local<<" gridcells"<<endl;
	
	clGridCell mycell;		// allocate a grid cell object for all grid cells that one core will have to run
	

	


	
	//--- break the infinite loop once all years are done for all local gridcells
	for (int cidx = 0; cidx < ncells_local;cidx++)
		{
			int loncount = cellindex(cidx) % IData.xread;
			int latcount = cellindex(cidx) / IData.xread;
			mycell.initialize(runid,fire,&IData,latcount,loncount);		
			// 	run aDGVM simulations	

			for (year=0 ; year<YEARS_TO_RUN; year++ )
			  {			
				if (nrank==0) {printf("\r Working on year %i of gridcell %i/%i...", year,cellindex(cidx)+1,ncells);
				fflush(stdout);}
				
				
				if (day==0) mycell.RainFallYear();
				for (day=0 ; day<365; day++ )
				{
					#ifdef USE_ADGVM_WATER_MODEL
					if (mycell.getplant()->max_root_depth[0]>=0) BucketIn(mycell.getsoil(),mycell.getclim(day), mycell.getplant() );
					#endif
/*					for (int i=0;i<N_SOIL_LAYERS;i++) cout << mycell.getsoil()[i].tfc << " ";
					cout <<endl;*/
					if(mycell.runDailyProcesses(year,day))
						if (nrank==0)		
						{
						#ifdef NC_OUTPUT
							pop_out->FillOutputBuffer(mycell.getplantpop()->pop_data);
						}
						else
						{ 
						#ifdef WITH_MPI
							MPI::COMM_WORLD.Send(mycell.getplantpop()->pop_data, BUFF_SIZE*POP_DATA_SIZE+2, MPI::FLOAT, 0, 2);
						#endif //WITH_MPI
						#endif //NC_OUTPUT
						}

					#ifdef USE_ADGVM_WATER_MODEL
					if (mycell.getplant()->max_root_depth[0]>=0) BucketOut( mycell.getsoil(), mycell.getplant() );	
					#endif  		
				}
				if(mycell.runAnnualProcesses(year))  // run annual processes
					if (nrank==0)		
					{
					#ifdef NC_OUTPUT
						trait_out->FillOutputBuffer(mycell.getplantpop()->trait_data);
					}
					else
					{ 
					#ifdef WITH_MPI
						MPI::COMM_WORLD.Send(mycell.getplantpop()->trait_data, MAX_POP_SIZE*TRAIT_VEC_LENGTH+2, MPI::FLOAT, 0, 1);
					#endif //WITH_MPI
					#endif //NC_OUTPUT
					}
		#ifdef NC_OUTPUT		
		#ifdef WITH_MPI
		// receive data from slaves until all data is received
			while(MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE,1))
			{	
				float* trait_data= new float[MAX_POP_SIZE*TRAIT_VEC_LENGTH+2];
				MPI::COMM_WORLD.Recv(trait_data, MAX_POP_SIZE*TRAIT_VEC_LENGTH+2, MPI::FLOAT, MPI_ANY_SOURCE, 1);
				trait_out->FillOutputBuffer(trait_data);
				delete[] trait_data;
			}

			while(MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE,2))
			{	
				float* pop_data = new float[BUFF_SIZE*POP_DATA_SIZE+2];
				MPI::COMM_WORLD.Recv(pop_data, BUFF_SIZE*POP_DATA_SIZE+2, MPI::FLOAT, MPI_ANY_SOURCE, 2);
				pop_out->FillOutputBuffer(pop_data);
				delete[] pop_data;

			}
		#endif
		#endif


			   } //end of year loop

		mycell.finalize();
		}
if (nrank==0)	
#ifdef NC_OUTPUT
while(true)
{
		#ifdef WITH_MPI
		// receive data from slaves until all data is received
			while(MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE,1))
			{	
				float* trait_data= new float[MAX_POP_SIZE*TRAIT_VEC_LENGTH+2];
				MPI::COMM_WORLD.Recv(trait_data, MAX_POP_SIZE*TRAIT_VEC_LENGTH+2, MPI::FLOAT, MPI_ANY_SOURCE, 1);
				trait_out->FillOutputBuffer(trait_data);
				delete[] trait_data;
			}

			while(MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE,2))
			{	
				float* pop_data = new float[BUFF_SIZE*POP_DATA_SIZE+2];
				MPI::COMM_WORLD.Recv(pop_data, BUFF_SIZE*POP_DATA_SIZE+2, MPI::FLOAT, MPI_ANY_SOURCE, 2);
				pop_out->FillOutputBuffer(pop_data);
				delete[] pop_data;

			}
		#endif
	if (pop_out->BufferFull() && trait_out->BufferFull()) break;
}// end "while(true) loop"
#endif

	
	time(&end);
	double timedif = difftime(end,start);
	cout<<"Node " << nrank<<" finished calculations"<<endl;
	
	#ifdef NC_OUTPUT                                                           
	if (nrank==0){ delete trait_out; delete pop_out; }                  
	#endif                                                                     
	cout << "Time to run program was " << timedif/60. << " minutes." << endl;  
	
	#ifdef WITH_MPI
	MPI::Finalize(); 		
	#endif
	return 0;
}
catch(const libconfig::SettingNotFoundException &nfex)
	{
	//this captures all libconfig "SettingNotFound" errors and prints a human readable error message.
	cerr << "Error: Setting '" <<nfex.getPath()<< "' not found in configuration file." << endl;
	return(EXIT_FAILURE);
	}
}




