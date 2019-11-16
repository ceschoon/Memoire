#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "SolidPhase.h"
#include "SolidPhase.cpp"
#include "../LJUniformFluid/UniformFluid.h"
#include "../LJUniformFluid/UniformFluid.cpp"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the free energy of the FCC Solid at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT = 1;
	double mu = -1;
	
	Options options;
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////// Run DFT computations for the solid ///////////////////
	
	double aLattice = 0;
	double densitySolid = 0;
	double numParticlesInCell = 0;
	double freeEnergySolid = 0;
	double aVdWSolid = 0;
	double hsdSolid = 0;
	bool successSolid = true;
	
	minOverNpoints(argc, argv, log, kT, mu, aLattice, densitySolid, 
	               numParticlesInCell, freeEnergySolid, aVdWSolid, 
	               hsdSolid, successSolid);
	
	// Print results in log file
	
	if (successSolid)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Lattice Constant = " << aLattice << endl;
		log << "Number of Particles in FCC Cell = " << numParticlesInCell << endl;
		log << "Solid Density = " << densitySolid << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Minimisation Npoints FAILED" << endl;
	}
	
	///////////////// Run DFT computations for the fluid ///////////////////
	
	// TODO: use bisection method?? 
	//       current problem if we are in vapor phase and start in the liquid
	
	bool successFluid = true;
	
	double densityFluid = findRootdOmegadRho(1.5, kT, mu, aVdWSolid, hsdSolid, successFluid);
	double freeEnergyFluid = uniformFluidOmega(kT, mu, aVdWSolid, hsdSolid, densityFluid);
	
	// Print results in log file
	
	if (successFluid)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Fluid Density = " <<densityFluid << endl;
		log << "Fluid Free Energy = " << freeEnergyFluid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Fluid computation FAILED" << endl;
	}
	
	return 1;
}

