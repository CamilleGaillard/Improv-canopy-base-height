#ifndef SeedBankClass_h___
#define SeedBankClass_h___
#include "SeedBankClassConstants.h"
#include "SeedClass.h"
class clGridCell;
class clSeedBank
{
    private:
        clSeed* sb_;   // stores seeds produces by trees
        int seed_number_;         // indicates how full seedbank is
    	int species_seed_num_[SPECIES_NUM];      // number of seeds of a species
    	double seed_weight_;
        
    public:
        clSeedBank();
	~clSeedBank();
	void initialize( clGridCell *gridcell_);
        void addSeed( clSeed sd );
        void emptySeedBank();
        clSeed getRandomSeed();
        int  getSeedNumber()          { return seed_number_; }
        double getSeedWeight()		{ return seed_weight_; }
        void runOptim();
        void runMutation();
        void runCrossover();

};

#endif
