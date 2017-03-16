#ifndef NcOutputClass_h___
#define NcOutputClass_h___
#include "adgvm.h"
#include "PlantPopClass.h"


const unsigned int POP_TYPE=1;
const unsigned int TRAIT_TYPE=2;

																				
class NcOutputClass
{
	
	private:
		std::map<std::string,int> pop_written;
		std::map<std::string,int> trait_written;
		int BUFF_SIZEx;
		int ncid;
		int status;
		int dimid;
		int timestamp_dimid;
		int layer_dimid;
		int vtype_dimid;
		int ind_dimid;
		int lat_dimid;
		int lon_dimid;
		int layer_varid;
		int vtype_varid;
		int lat_varid;
		int lon_varid;
		int ind_varid;
		int time_varid;
		int varid;
		std::map <std::string,int> data_varid;
		unsigned int contenttype;
		unsigned int latcount;
		unsigned int loncount;
		size_t *recnum;
		int data_size; //number of recods in the data_varid and data_mll arrays
		void WriteHeader(std::map <std::string,vector<float>> data);
		bool firstrun;
	
	public:		
		NcOutputClass(const char* filename, unsigned int ctype, MyInData *IData);			// constructor prototype
		~NcOutputClass();			// destructor
		int BufferFull();
		int WriteBuffer(std::map <std::string,vector<float>> data);				// not sure if this should be a member function of NcOutputClass or just a separate function? Tend to think it should be member function
//		int FillOutputBuffer(MLL_array* data, MLL_array* data2, int latcount, int loncount);
//		int FillOutputBuffer(std::string, vector<float>);
						
					
};				

#endif															
