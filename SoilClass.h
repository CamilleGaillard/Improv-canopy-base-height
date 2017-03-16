#ifndef SoilClass_h___
#define SoilClass_h___

// This class implements the Yasso soil model
// following Liski et al 2005 Ecological Modelling

//#include "adgvm.h"
#include "SoilClassGlobals.h"



class clSoil
{
    private:
        double x_fwl_;      // fine woody litter
        double x_cwl_;      // coarse woody litter
        double x_ext_;      // extractives
        double x_cel_;      // celluloses
        double x_lig_;      // lignin-like compounds
        double x_hum1_;     // humus
        double x_hum2_;     // more recalcitrant humus

        double r_ext_;      // CO2 release from extractives
        double r_cel_;      // CO2 release from cellulose
        double r_lig_;      // CO2 release from lignin-like
        double r_hum1_;     // CO2 release from humus1
        double r_hum2_;     // CO2 release from humus2

        double u_nwl_;      // non-woody litter input
        double u_fwl_;      // fine woodly litter input
        double u_cwl_;      // coarse woody litter input

        double sens_hum_;   // tmp and drought sensitivity of k for humus
        double sens_oth_;   // tmp and drought sensitivity of k and a for other compartments

    public:
        clSoil();
        clSoil( double soil_carbon);
        ~clSoil() {};

	void   Initialize( double soil_carbon );

        double UpdateCarbonPools( double u_nwl, double u_fwl, double u_cwl, double T, double D );
        double GetCarbonInput();
        double GetCarbonStored();
        double GetCarbonRelease();
        double GetCarbonInHumus();
        void   PrintCarbonPools();
//         void   UpdateCarbonPoolsAfterFire( double above_below_fraction );

        double GetXfwl() { return x_fwl_;  }
        double GetXcwl() { return x_cwl_;  }
        double GetXext() { return x_ext_;  }
        double GetXcel() { return x_cel_;  }
        double GetXlig() { return x_lig_;  }
        double GetXhu1() { return x_hum1_; }
        double GetXhu2() { return x_hum2_; }
        double GetUnwl() { return u_nwl_;  }
        double GetUfwl() { return u_fwl_;  }
        double GetUcwl() { return u_cwl_;  }
        double GetRext() { return r_ext_;  }
        double GetRcel() { return r_cel_;  }
        double GetRlig() { return r_lig_;  }
        double GetRhu1() { return r_hum1_; }
        double GetRhu2() { return r_hum2_; }

    private:
        void getSensHum( double T, double D );
        void getSensOth( double T, double D );
};



#endif














