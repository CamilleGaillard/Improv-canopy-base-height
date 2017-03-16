#include <fstream>
//#include <stdlib.h>
#include <cmath>
#include "MyMath.h"
#include "SeedBankClass.h"
#include "SeedClass.h"
#include "PlantClassConstants.h"

using namespace::std;

//---------------------------------------------------------------------------------------

clSeedBank::clSeedBank()
{
    sb_ = new clSeed[MAX_SEED_BANK_SIZE];
    seed_number_ = 0;   // empty seedbank
    seed_weight_ = 0.;    
    for ( int i=0; i<SPECIES_NUM; i++ ) species_seed_num_[i] = 0;
}
//---------------------------------------------------------------------------------------

clSeedBank::~clSeedBank()
{
    delete[] sb_;
}
//---------------------------------------------------------------------------------------


void clSeedBank::initialize( clGridCell * gridcell_)
{
    seed_number_ = 0;   // empty seedbank
    seed_weight_ = 0.;   
    for ( int i=0; i<MAX_SEED_BANK_SIZE; i++ ) sb_[i].initialize(gridcell_); //FLAGFLAG
}
//---------------------------------------------------------------------------------------

void clSeedBank::addSeed( clSeed sd )
{
    if ( seed_number_<MAX_SEED_BANK_SIZE )
    {
        sb_[seed_number_] = sd;
        seed_number_++;
        seed_weight_ += sd.getSeedWeight();
        species_seed_num_[sd.getSpeciesIndex()]++;
    }
        
    return;
}

//---------------------------------------------------------------------------------------

void clSeedBank::emptySeedBank()
{
    for ( int i=0; i<SPECIES_NUM; i++ ) species_seed_num_[i] = 0;
	seed_number_ = 0;
	seed_weight_ = 0.;	
	return;
}

//---------------------------------------------------------------------------------------

clSeed clSeedBank::getRandomSeed()
{
	return sb_[(int) RUnif(0, seed_number_)];
}

//---------------------------------------------------------------------------------------

void clSeedBank::runMutation()
{
	for ( int i=0; i<seed_number_; i++ )
	{
		
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckEvergreen(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRainLight() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckSla(         RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getSla() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocRoot(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocRoot() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocLeaf(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocLeaf() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRainThrOn(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRainThrOn() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRainThrOff(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRainThrOff() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckLightThrOn(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getLightThrOn() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckLightThrOff( RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getLightThrOff() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRainLight(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRainLight() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckBiomasPar1(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getBiomasPar1() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckBiomasPar2(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getBiomasPar2() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRotfrmPar1(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRotfrmPar1() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRotfrmPar2(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRotfrmPar2() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckRotfrmMaxD(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getRotfrmMaxD() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckWoodDensity( RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getWoodDensity() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckLightExt(    RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getLightExt() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocWood(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocWood() );
		#ifdef W_BRANCHES
			if ( RUnif()<MUT_PROBABILITY ) {
				sb_[i].setCheckAllocBranch( RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocBranch()/sb_[i].getAllocWood(), sb_[i].getAllocWood());	// MP; modification of the "raw" branch allocation trait
				sb_[i].setCheckAllocStem( sb_[i].getAllocWood(), sb_[i].getAllocBranch());																// MP; think I need to adjust branch and stem at the same time, make sure they add up to alloc_wood again in the end
			}	
		#endif	
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocBark(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocBark() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocRepr(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocRepr() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckAllocStor(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getAllocStor() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckSeedWeight(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getSeedWeight() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckCanfrmPar1(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getCanfrmPar1() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckCanfrmPar2(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getCanfrmPar2() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckCanfrmHBase( RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getCanfrmHBase() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckP50(         RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getP50() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckStorToWood(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getStorToWood() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckStorToLeaf(  RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getStorToLeaf() );
		if ( RUnif()<MUT_PROBABILITY )sb_[i].setCheckStemCount(   RUnif( MUT_RATE_LO, MUT_RATE_UP )  *sb_[i].getStemCount()  ); // NEW CAMILLE		
		sb_[i].setCheckAllocParameters();
		
	}
	
	return;
}

// -----------------------------------------------------------------------
// --- crossover ---------------------------------------------------------
// -----------------------------------------------------------------------

void clSeedBank::runCrossover()
{
	int cum_seed_num[SPECIES_NUM+1];
	for ( int i=0; i<SPECIES_NUM; i++ ) cum_seed_num[i] = -1;
	cum_seed_num[0] = 0;
	int ind = 0;
	int xs1= 0;
	int xs2= 0;
	int xs3 = 0;
	
	for ( int i=0; i<SPECIES_NUM; i++ )
	{
		if ( species_seed_num_[i]>0 )
		{
			cum_seed_num[ind+1] = cum_seed_num[ind]+species_seed_num_[i];
			ind++;
		}
	}
	
	
	int sp_ind_cur = sb_[0].getSpeciesIndex();
	int sp_ind_old = sp_ind_cur;
	int cnt_cum    = 0;
	
	for ( int i=0; i<seed_number_; i++ )
	{
		sp_ind_cur = sb_[i].getSpeciesIndex();
		if ( sp_ind_cur != sp_ind_old ) cnt_cum++;
		
		xs1 = RUnifInt(seed_number_);
		//xs2 = (int)MyRound( RUnif( cum_seed_num[cnt_cum]-0.49, cum_seed_num[cnt_cum+1]-1.+0.49 ) ); // Useless now
		//xs3 = (int)MyRound( RUnif( cum_seed_num[cnt_cum]-0.49, cum_seed_num[cnt_cum+1]-1.+0.49 ) ); // Useless now
		
		//cout << setw(5) << i << setw(5) << xs1 << setw(5) << sb_[i].getSpeciesIndex() << setw(5) << sb_[xs1].getSpeciesIndex() << endl;
		
		sp_ind_old = sp_ind_cur;
		
		// species with a low number of seeds can crossover with any individual
		if ( RUnif() < 1.-(species_seed_num_[sp_ind_cur]/(double)MAX_SEED_NUM*5.) )
		{
			// 			cout << species_seed_num_[sp_ind_cur] << "  " << 1.-(species_seed_num_[sp_ind_cur]/(double)MAX_SEED_NUM*5.) << endl;
			xs1 = (int)MyRound( RUnif( 0.49, seed_number_+0.49 ) );
			//xs1 = RUnifInt(seed_number_); //FLAG MP proposed equivalent solution using integer random number generator directly
		}
		
// 		cout << setw(10) << cum_seed_num[cnt_cum] << setw(10) << cum_seed_num[cnt_cum+1]-1 << setw(10) << xs1 << setw(10) << xs2 << setw(10) << xs3 << setw(10) << sb_[xs1].getSpeciesIndex() << setw(10) << sb_[xs2].getSpeciesIndex() << setw(10) << sb_[xs3].getSpeciesIndex() << setw(10) << sb_[i].getSpeciesIndex() << endl;
		
		// 		cout << setw(5) << sb_[i].getSpeciesIndex() << setw(5) << sb_[i].getVegType() << setw(15) << sb_[i].getSla()  << setw(15) << sb_[xs1].getSla()  << setw(15) << sb_[xs2].getSla()  << setw(15) << sb_[xs3].getSla() << setw(15) << tv << endl;
		
		
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckEvergreen(   sb_[i].getEvergreen()   + RUnif()*(sb_[xs1].getEvergreen()  -sb_[i].getEvergreen())  );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckSla(         sb_[i].getSla()         + RUnif()*(sb_[xs1].getSla()        -sb_[i].getSla())  );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocRoot(   sb_[i].getAllocRoot()   + RUnif()*(sb_[xs1].getAllocRoot()  -sb_[i].getAllocRoot()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocLeaf(   sb_[i].getAllocLeaf()   + RUnif()*(sb_[xs1].getAllocLeaf()  -sb_[i].getAllocLeaf()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRainThrOn(   sb_[i].getRainThrOn()   + RUnif()*(sb_[xs1].getRainThrOn()  -sb_[i].getRainThrOn()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRainThrOff(  sb_[i].getRainThrOff()  + RUnif()*(sb_[xs1].getRainThrOff() -sb_[i].getRainThrOff()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckLightThrOn(  sb_[i].getLightThrOn()  + RUnif()*(sb_[xs1].getLightThrOn() -sb_[i].getLightThrOn()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckLightThrOff( sb_[i].getLightThrOff() + RUnif()*(sb_[xs1].getLightThrOff()-sb_[i].getLightThrOff()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRainLight(   sb_[i].getRainLight()   + RUnif()*(sb_[xs1].getRainLight()  -sb_[i].getRainLight()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckBiomasPar1(  sb_[i].getBiomasPar1()  + RUnif()*(sb_[xs1].getBiomasPar1() -sb_[i].getBiomasPar1()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckBiomasPar2(  sb_[i].getBiomasPar2()  + RUnif()*(sb_[xs1].getBiomasPar2() -sb_[i].getBiomasPar2()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRotfrmPar1(  sb_[i].getRotfrmPar1()  + RUnif()*(sb_[xs1].getRotfrmPar1() -sb_[i].getRotfrmPar1()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRotfrmPar2(  sb_[i].getRotfrmPar2()  + RUnif()*(sb_[xs1].getRotfrmPar2() -sb_[i].getRotfrmPar2()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckRotfrmMaxD(  sb_[i].getRotfrmMaxD()  + RUnif()*(sb_[xs1].getRotfrmMaxD() -sb_[i].getRotfrmMaxD()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckWoodDensity( sb_[i].getWoodDensity() + RUnif()*(sb_[xs1].getWoodDensity()-sb_[i].getWoodDensity()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckLightExt(    sb_[i].getLightExt()    + RUnif()*(sb_[xs1].getLightExt()   -sb_[i].getLightExt()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocWood(   sb_[i].getAllocWood()   + RUnif()*(sb_[xs1].getAllocWood()  -sb_[i].getAllocWood()) );
		#ifdef W_BRANCHES
			if (RUnif()<CROSSOVER_CR) {
				sb_[i].setCheckAllocBranch( sb_[i].getAllocBranch()/sb_[i].getAllocWood()   + RUnif()*(sb_[xs1].getAllocBranch()/sb_[xs1].getAllocWood()  -sb_[i].getAllocBranch()/sb_[i].getAllocWood()), sb_[i].getAllocWood() );	// MP; modification of the "raw" branch allocation trait
				sb_[i].setCheckAllocStem( sb_[i].getAllocWood(), sb_[i].getAllocBranch());																// MP; think I need to adjust branch and stem at the same time, make sure they add up to alloc_wood again in the end
			}
		#endif
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocBark(   sb_[i].getAllocBark()   + RUnif()*(sb_[xs1].getAllocBark()  -sb_[i].getAllocBark()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocRepr(   sb_[i].getAllocRepr()   + RUnif()*(sb_[xs1].getAllocRepr()  -sb_[i].getAllocRepr()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckAllocStor(   sb_[i].getAllocStor()   + RUnif()*(sb_[xs1].getAllocStor()  -sb_[i].getAllocStor()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckSeedWeight(  sb_[i].getSeedWeight()  + RUnif()*(sb_[xs1].getSeedWeight() -sb_[i].getSeedWeight()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckCanfrmPar1(  sb_[i].getCanfrmPar1()  + RUnif()*(sb_[xs1].getCanfrmPar1() -sb_[i].getCanfrmPar1()));
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckCanfrmPar2(  sb_[i].getCanfrmPar2()  + RUnif()*(sb_[xs1].getCanfrmPar2() -sb_[i].getCanfrmPar2()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckCanfrmHBase( sb_[i].getCanfrmHBase() + RUnif()*(sb_[xs1].getCanfrmHBase()-sb_[i].getCanfrmHBase()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckP50(         sb_[i].getP50()         + RUnif()*(sb_[xs1].getP50()        -sb_[i].getP50()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckStorToWood(  sb_[i].getStorToWood()  + RUnif()*(sb_[xs1].getStorToWood() -sb_[i].getStorToWood()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckStorToLeaf(  sb_[i].getStorToLeaf()  + RUnif()*(sb_[xs1].getStorToLeaf() -sb_[i].getStorToLeaf()) );
		if (RUnif()<CROSSOVER_CR) sb_[i].setCheckStemCount(   sb_[i].getStemCount()   + RUnif()*(sb_[xs1].getStemCount()  -sb_[i].getStemCount())  ); // NEW CAMILLE
		sb_[i].setCheckAllocParameters();		
		
	}
	
	return;
}











