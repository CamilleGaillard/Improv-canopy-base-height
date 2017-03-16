#ifndef GridCellClassStructures_h___
#define GridCellClassStructures_h___
struct GridCellClimateData {
		double par;
		double rad;
		double tmp;
		double hum;
		double wnd;
		double apr;
		double co2;
		double pre;
		double sun;
		};

struct GridCellSoilData {
		double twp;	// theta wilting point
		double tfc;	// theta field capacity
		double soc;
		double son;
		double ksat;
		double smo;
		double root_bm;
		double sat_water_cont;
		double soil_texture;
//		double perc_exponent;
		};
struct GridCellPlantData {
		//locally calculated values
		double lai;
		double gsc;
		double cc3;
		double cc4;
		double ctr;
		double etr;
		double A0C3;
		double A0C4;
		double RLC3;
		double RLC4;
		double max_root_depth[2]; 	// NEW LIAM; 0 for trees, 1 for grasses
};
#endif
