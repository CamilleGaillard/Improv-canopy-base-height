#include "SoilClassGlobals.h"
#include "SoilClass.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace::std;


clSoil::clSoil()
{
    x_fwl_  = 0.;    // fine woody litter
    x_cwl_  = 0.;    // coarse woody litter
    x_ext_  = 0.;    // extractives
    x_cel_  = 0.;    // celluloses
    x_lig_  = 0.;    // lignin-like compounds
    x_hum1_ = 1.;    // humus
    x_hum2_ = 1.;    // more recalcitrant humus

    r_ext_  = 0.;     // CO2 release from extractives
    r_cel_  = 0.;     // CO2 release from cellulose
    r_lig_  = 0.;     // CO2 release from lignin-like
    r_hum1_ = 0.;     // CO2 release from humus1
    r_hum2_ = 0.;     // CO2 release from humus2

    u_nwl_  = 0.;     // non-woody litter input
    u_fwl_  = 0.;     // fine woodly litter input
    u_cwl_  = 0.;     // coarse woody litter input
}

clSoil::clSoil(double soil_carbon)
{
//     x_fwl_  = 0.;    // fine woody litter
//     x_cwl_  = 0.;    // coarse woody litter
//     x_ext_  = 0.;    // extractives
//     x_cel_  = 0.;    // celluloses
//     x_lig_  = 0.;    // lignin-like compounds
    x_fwl_  = soil_carbon*0.01152994;    // fine woody litter
    x_cwl_  = soil_carbon*0.06949435;    // coarse woody litter
    x_ext_  = soil_carbon*0.00851789;    // extractives
    x_cel_  = soil_carbon*0.06624854;    // celluloses
    x_lig_  = soil_carbon*0.06360159;    // lignin-like compounds
    x_hum1_ = soil_carbon*0.2377162;     // humus
    x_hum2_ = soil_carbon*0.5428915;     // more recalcitrant humus

    r_ext_  = 0.;     // CO2 release from extractives
    r_cel_  = 0.;     // CO2 release from cellulose
    r_lig_  = 0.;     // CO2 release from lignin-like
    r_hum1_ = 0.;     // CO2 release from humus1
    r_hum2_ = 0.;     // CO2 release from humus2

    u_nwl_  = 0.;     // non-woody litter input
    u_fwl_  = 0.;     // fine woodly litter input
    u_cwl_  = 0.;     // coarse woody litter input
}




// ---------------------------------------------------------------
// Initialize soil carbon pools with soil_carbon (units??)
void clSoil::Initialize(double soil_carbon)
{
//     x_fwl_  = 0.;    // fine woody litter
//     x_cwl_  = 0.;    // coarse woody litter
//     x_ext_  = 0.;    // extractives
//     x_cel_  = 0.;    // celluloses
//     x_lig_  = 0.;    // lignin-like compounds
    x_fwl_  = soil_carbon*0.01152994;    // fine woody litter
    x_cwl_  = soil_carbon*0.06949435;    // coarse woody litter
    x_ext_  = soil_carbon*0.00851789;    // extractives
    x_cel_  = soil_carbon*0.06624854;    // celluloses
    x_lig_  = soil_carbon*0.06360159;    // lignin-like compounds
    x_hum1_ = soil_carbon*0.2377162;     // humus
    x_hum2_ = soil_carbon*0.5428915;     // more recalcitrant humus

    r_ext_  = 0.;     // CO2 release from extractives
    r_cel_  = 0.;     // CO2 release from cellulose
    r_lig_  = 0.;     // CO2 release from lignin-like
    r_hum1_ = 0.;     // CO2 release from humus1
    r_hum2_ = 0.;     // CO2 release from humus2

    u_nwl_  = 0.;     // non-woody litter input
    u_fwl_  = 0.;     // fine woodly litter input
    u_cwl_  = 0.;     // coarse woody litter input

}




// ---------------------------------------------------------------
// update carbon pools (run yasso)
// function takes input values in biomass (kg/m^2) and
// transforms it to C (kg/m^2)
double clSoil::UpdateCarbonPools( double u_nwl, double u_fwl, double u_cwl, double T, double D )
{
    getSensHum( T, D );
    getSensOth( T, D );

    u_nwl_  = 0.44*u_nwl;     // non-woody litter input
    u_fwl_  = 0.44*u_fwl;     // fine woodly litter input
    u_cwl_  = 0.44*u_cwl;     // coarse woody litter input


    r_ext_  = (1.-YASSO_P_EXT)  * YASSO_K_EXT * x_ext_;
    r_cel_  = (1.-YASSO_P_CEL)  * YASSO_K_CEL * x_cel_;
    r_lig_  = (1.-YASSO_P_LIG)  * YASSO_K_LIG * x_lig_;
    r_hum1_ = (1.-YASSO_P_HUM1) * YASSO_K_HUM1* x_hum1_;
    r_hum2_ =                     YASSO_K_HUM2* x_hum2_;


    x_fwl_  += (u_fwl_ - sens_oth_*YASSO_A_FWL*x_fwl_);

    x_cwl_  += (u_cwl_ - sens_oth_*YASSO_A_CWL*x_cwl_);

    x_ext_  += (u_nwl_*YASSO_C_NWL_EXT + YASSO_C_FWL_EXT*sens_oth_*YASSO_A_FWL*x_fwl_
            + YASSO_C_CWL_EXT*sens_oth_*YASSO_A_CWL*x_cwl_ - sens_oth_*YASSO_K_EXT*x_ext_);

    x_cel_  += (u_nwl_*YASSO_C_NWL_CEL + YASSO_C_FWL_CEL*sens_oth_*YASSO_A_FWL*x_fwl_
            + YASSO_C_CWL_CEL*sens_oth_*YASSO_A_CWL*x_cwl_ - sens_oth_*YASSO_K_CEL*x_cel_);

    x_lig_  += (u_nwl_*YASSO_C_NWL_LIG + YASSO_C_FWL_LIG*sens_oth_*YASSO_A_FWL*x_fwl_
            + YASSO_C_CWL_LIG*sens_oth_*YASSO_A_CWL*x_cwl_ - sens_oth_*YASSO_K_LIG*x_lig_
            + YASSO_P_EXT*sens_oth_*YASSO_K_EXT*x_ext_ + YASSO_P_CEL*sens_oth_*YASSO_K_CEL*x_cel_);

    x_hum1_ += (YASSO_P_LIG*sens_oth_*YASSO_K_LIG*x_lig_ - sens_hum_*YASSO_K_HUM1*x_hum1_);

    x_hum2_ += (YASSO_P_HUM1*sens_hum_*YASSO_K_HUM1*x_hum1_ - sens_hum_*YASSO_K_HUM2*x_hum2_);

    return GetCarbonRelease();
}



double clSoil::GetCarbonInput()
{
    return u_nwl_ + u_fwl_ + u_cwl_;
}


double clSoil::GetCarbonStored()
{
    return x_fwl_ + x_cwl_ + x_ext_ + x_cel_ + x_lig_ + x_hum1_ + x_hum2_;
}


double clSoil::GetCarbonRelease()
{
    return r_ext_ + r_cel_ + r_lig_ + r_hum1_ + r_hum2_;
}

// ------------------------------------------------------------------------
// returns soil carbon stroed in humus, in kgC/m^2
double clSoil::GetCarbonInHumus()
{
    return r_hum1_ + r_hum2_;
}

void clSoil::PrintCarbonPools()
{
    cout << setw( 7) << "SOILC" <<
            setw(14) << x_fwl_  <<
            setw(14) << x_cwl_  <<
            setw(14) << x_ext_  <<
            setw(14) << x_cel_  <<
            setw(14) << x_lig_  <<
            setw(14) << x_hum1_ <<
            setw(14) << x_hum2_ <<
            setw(14) << r_ext_  <<
            setw(14) << r_cel_  <<
            setw(14) << r_lig_  <<
            setw(14) << r_hum1_ <<
            setw(14) << r_hum2_ <<
            setw(14) << u_nwl_  <<
            setw(14) << u_fwl_  <<
            setw(14) << u_cwl_  <<
            endl;

    return;
}


void clSoil::getSensHum( double T, double D )
{
//     sens_hum_ = 1.+ YASSO_S_HUM*YASSO_BETA*(T-YASSO_T0) + YASSO_GAMMA*(D-YASSO_D0);
    sens_hum_ = 1.+ YASSO_S_HUM*YASSO_BETA*(T-YASSO_T0) + YASSO_GAMMA*D;
}

void clSoil::getSensOth( double T, double D )
{
//     sens_oth_ = 1.+             YASSO_BETA*(T-YASSO_T0) + YASSO_GAMMA*(D-YASSO_D0);
    sens_oth_ = 1.+             YASSO_BETA*(T-YASSO_T0) + YASSO_GAMMA*D;
}





