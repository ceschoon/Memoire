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
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

#include "../SolidFCC.h"
#include "../SolidPhase.h"
#include "../UniformPhases.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the free energy of the FCC Solid at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT = 1;
	double mu = -1;
	int Npoints = 125;
	double rhoFluidMax = 1.5;
	double rhoFluidMin = 0.001;
	
	Options options;
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	options.addOption("Npoints", &Npoints);
	options.addOption("RhoFluidMax", &rhoFluidMax);
	options.addOption("RhoFluidMin", &rhoFluidMin);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////// Run DFT computations for the solid ///////////////////
	
	double densitySolid = 0;
	double numParticlesInCell = 0;
	double freeEnergySolid = 0;
	double aVdWSolid = 0;
	double hsdSolid = 0;
	bool useSnapshot = false;
	bool successSolid = true;
	
	DFTcomputation(argc, argv, log, kT, mu, Npoints, densitySolid, 
	               numParticlesInCell, freeEnergySolid, aVdWSolid, 
	               hsdSolid, useSnapshot, successSolid);
	
	// Print results in log file
	
	if (successSolid)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Number of Particles in FCC Cell = " << numParticlesInCell << endl;
		log << "Solid Density = " << densitySolid << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Solid Computation FAILED" << endl;
	}
	
	///////////////// Run DFT computations for the fluid ///////////////////
	
	bool successFluid = true;
	double densityFluid = 0;
	double freeEnergyFluid = 0;
	
	try
	{
		vector<double> roots;
		findRootsdOmegadRhoSpinodal(kT, mu, aVdWSolid, hsdSolid, 
			rhoFluidMin, rhoFluidMax, roots);
		
		if (roots.size()==1)
		{
			densityFluid = roots[0];
			freeEnergyFluid = uniformOmega(kT, mu, aVdWSolid, hsdSolid, densityFluid);
		}
		else if (roots.size()==2)
		{
			double freeEnergyVapor  = uniformOmega(
				kT, mu, aVdWSolid, hsdSolid, roots[0]);
			double freeEnergyLiquid = uniformOmega(
				kT, mu, aVdWSolid, hsdSolid, roots[1]);
			
			if (freeEnergyVapor<freeEnergyLiquid)
			{
				densityFluid = roots[0];
				freeEnergyFluid = freeEnergyVapor;
			}
			else
			{
				densityFluid = roots[1];
				freeEnergyFluid = freeEnergyLiquid;
			}
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Anormal number of roots for the fluid dOmega/drho" << endl;
			successFluid = false;
		}
	}
	catch (...)
	{
		successFluid = false;
	}
	
	if (successFluid)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Fluid Density = " << densityFluid << endl;
		log << "Fluid Free Energy = " << freeEnergyFluid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Fluid computation FAILED" << endl;
	}
	
	return 1;
}
