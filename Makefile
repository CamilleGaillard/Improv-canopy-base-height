CC=gcc
CXX=mpic++
FC = gfortran
RM=rm -f
OPTIMFLAGS= -O3
DEBUGFLAGS= -Warray-bounds -DWITH_MPI #-Wunused-variable #-Wmaybe-uninitialized #-Wunused-variable #-Wall #-pg
CPPFLAGS=-I/usr/include -std=gnu++11 $(OPTIMFLAGS) -DUSE_ADGVM_WATER_MODEL -march=native $(DEBUGFLAGS) -DREPRODUCIBLE -DNC_OUTPUT #-DW_BRANCHES #-DNANCHECK
LDFLAGS= -L/usr/include $(OPTIMFLAGS) $(DEBUGFLAGS) 
LDLIBS=-lnetcdf -lconfig++ -L/usr/lib 

OBJS= adgvm.o NcOutputClass.o NcInputClasses.o GridCellClass.o PlantPopClass.o \
PlantClass.o Leaf.o SeedBankClass.o SeedClass.o Radiation.o PenmanMonteith.o \
AdgvmWaterModel.o SoilClass.o

TOBJS= adgvm.o NcOutputClass.o NcInputClasses.o GridCellClass.o PlantPopClass.o \
PlantClass.o Leaf.o SeedBankClass.o SeedClass.o Radiation.o PenmanMonteith.o \
AdgvmWaterModel.o SoilClass.o

#TOBJS= NcOutputClass.o Ncread2.o GridCellClass.o PlantPopClass.o \
#PlantClass.o Leaf.o SeedBankClass.o SeedClass.o Radiation.o PenmanMonteith.o \
#AdgvmWaterModel.o SoilClass.o


all: aDGVM

aDGVM: $(OBJS)
	$(CXX) $(LDFLAGS) -o aDGVM $(OBJS) $(LDLIBS) 

aTEST: $(TOBJS)
	$(CXX) $(LDFLAGS) -o aTEST $(TOBJS) $(LDLIBS) 


%.o: %.cpp 
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $<

							
clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) aDGVM

# the following ifdef options are available:
# -DUSE_ADGVM_WATER_MODEL	=>tells the model to use Liam's water model
# -DNC_OUTPUT			=> writes output to two separate nc-files whose file names are specified in adgvm.cpp (1 file for pop-data, 1 file for trait data)
# -DWITH_MPI				=> will run parallel version of model; to use this option, CHANGE COMPILER from g++ to mpi++ 
# -DDEBUG_ON			=> provides additional output to stdout from PlantClass.cpp to check on photosynthesis
# -DNANCHECK			=> activates internal checks/debug output for "NaN" values obtained throughout program
# -DREPRODUCIBLE	=> will seed the random number generator with a combination of longitude, latitude, and runid, which allows to create a reproducible but coordinate-specific unique random number sequence
# -DW_BRANCHES		=> will simulate trees with explicit crowns, i.e., branch biomass and trait-defined crown shape
