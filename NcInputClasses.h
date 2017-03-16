#ifndef NcInputClasses_h___
#define NcInputClasses_h___

/* This class does the basic netcdf file handling. It inquires dimension IDs and dimension lenghts
for longitude and latitude in the infile, checks that the user-specified coordinates are within the 
coordinate range available from the input file, and establishes the relationship between lon and lat values
and array indices. It uses the calculated lon and lat array indices to find the lon and lat coordinates in the 
nc-infile that are closest to the user-specified coordinates and defines the number of grid pixels (count) that
need to be read in as well as the starting indices in the ncfile.
*/

class NcInfile	
{	
	private:
		int status;
		int dimid;
		size_t xlen;		// lenght of dimension longitude
		size_t ylen;		// lenght of dimension latitude
		int varid;
		
	public:
		size_t* xsrt;
		size_t* xcount;
		size_t* ysrt;
		size_t* ycount;
		size_t* srt;
		size_t* count;
		size_t* srt2d;
		size_t* count2d;
		int xread;		// number of values to be read in x-direction (lon)
		int yread;		// number of values to be read in y-direction (lat)
		
	NcInfile(int ncid, int count1, double minlon, double maxlon, double minlat, double maxlat);	// constructor prototype
	
	// Destructor to free the location and counter variables again
	
	~NcInfile();   			//destructor prototype
};	

//------------------------------------------------------------------------------------------------------------
/* This class does the actual netcdf-handling required to read variables from the nc-infiles into program variables.
It allocates sufficient memory for the respective input variables and creates a user-friendly access-scheme for handling multi-dimensional
input variables. In addition, new multi-dimensional variables with the same user-friendly scheme can be created using the input variables'
scheme as empty template for these new variables. 
*/

class MLL_array
{
	private:
		int n[5];
		int dimensions;
	public:
		float *storage;
		unsigned int totsize;
		float invalid;
		
		MLL_array();  //empty constructor; create variables of type MLL_array from template variables read in from netcdf-file
		MLL_array(int ncid, const char *name,size_t *start, size_t *count);		// constructor; does actual netcdf-handling and allocates 
																														// sufficient memory for the respective variables
		~MLL_array();		// destructor;
		
		void EmptyFromTemplate(MLL_array * templatemllarray); // make this object identical to the template MLL_array object
		MLL_array(unsigned int newdimensions, size_t *newdimsizes);
		float& operator () (unsigned int p0, unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4);		// for 5D-variables
		float& operator () (unsigned int p0, unsigned int p1, unsigned int p2, unsigned int p3);		// for 4D-variables
		float& operator () (unsigned int p0, unsigned int p1, unsigned int p2);		// for 3D-variables
		float& operator () (unsigned int p0, unsigned int p1);		// for 2D-variables
		float& operator () (unsigned int p0);		// for 1D-variables
};	

//------------------------------------------------------------------------------------------------------------
/* This class uses the NcInfile class to initialize the start and count values required to read data from a specific 
nc-infile. It then creates the actual variables to store the input data and gets the input data by using the MLL_array class. It also handles the 
opening and closing of nc-infiles.
*/
class MyInData		
{	
	private:
		int ncid;
		int status;
	public:
		double minlon_;
		double maxlon_;
		double minlat_;
		double maxlat_;
		
	
		int xread;
		int yread;
		const unsigned int NDAYS_IN_MONTH[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
		unsigned int soil_layers; 
		double thickness[12];
		double depth[12];
		MLL_array climlons;
		MLL_array climlats;
		MLL_array tmp;
		MLL_array dtr;
		MLL_array sunp;
		MLL_array reh;
		MLL_array rd0;
		MLL_array wnd;
		MLL_array frost;
		MLL_array pre;
		MLL_array preCV;
		MLL_array elv;
		MLL_array soillons;
		MLL_array soillats;
		MLL_array theta_wp;
		MLL_array theta_fc;
		MLL_array k_sat;
		MLL_array soil_texture;
		MLL_array sat_water_cont;
		MLL_array perc_exponent;
		MLL_array bulk_dens;
		MLL_array soil_N;
		MLL_array soil_C;
		MLL_array root_depth;
		MLL_array tmp_min;
		MLL_array tmp_max;
		MLL_array ralpha;
		MLL_array rbeta;
		MLL_array pwet;
		MLL_array tmp_day;
		MLL_array atm_press;

		//derived variables from calcatmospheric
		MLL_array eS12;
		MLL_array eA12;
		MLL_array s12;
		MLL_array gama12;
		MLL_array VPD12;
		MLL_array rho12;
		MLL_array Q012;


		
		MyInData(double minlon, double maxlon, double minlat, double maxlat);		// constructor
};
#endif

