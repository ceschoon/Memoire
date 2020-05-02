#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <stdlib.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"
#include "DFT.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

#include "../SolidPhase.h"
#include "../UniformPhases.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the free energies at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDir = "data";
	double kT = 1;
	double mu = -1;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	////////////////// Compute free energy for the solid ///////////////////
	
	double freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min;
	
	int statusSolid = fixedkTMuSolid(kT, mu, argc, argv, log, 
		freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min);
	
	// Print results in log file
	
	if (statusSolid==0)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
		log << "Solid Density = " << densitySolid << endl;
		
		log << "Min Ngrid = " << Ngrid_min << endl;
		log << "Min Cvac  = " << Cvac_min  << endl;
		log << "Min alpha = " << alpha_min << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Solid Computation FAILED" << endl;
	}
	
	////////////////// Compute free energy for the fluid ///////////////////
	
	double freeEnergyVapour, freeEnergyLiquid;
	double densityVapour, densityLiquid;
	bool superCritical;
	
	int statusFluid = fixedkTMuFluid(kT, mu, argc, argv, log, 
		freeEnergyVapour, densityVapour, freeEnergyLiquid, densityLiquid,
		superCritical);
	
	// Print results in log file
	
	if (statusFluid==0 && superCritical)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Supercritical Fluid Free Energy = " << freeEnergyVapour << endl;
		log << "Supercritical Fluid Density = " << densityVapour << endl;
	}
	else if (statusFluid==0 && !superCritical)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Vapour Free Energy = " << freeEnergyVapour << endl;
		log << "Vapour Density = " << densityVapour << endl;
		log << endl;
		log << "Liquid Free Energy = " << freeEnergyLiquid << endl;
		log << "Liquid Density = " << densityLiquid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Fluid Computation FAILED" << endl;
	}
	
	return 1;
}
