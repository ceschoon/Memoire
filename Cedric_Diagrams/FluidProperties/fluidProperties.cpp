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

#include "../UniformPhases.h"
#include "../utilities.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program that computes the properties of the fluid phase corresponding
// (same kT,mu) to solid dft computations

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double maxFileIndex = 100;
	double rhoFluidMin = 0.001;
	double rhoFluidMax = 1.5;
	
	Options options;
	options.addOption("MaxFileIndex", &maxFileIndex);
	options.addOption("RhoFluidMin", &rhoFluidMin);
	options.addOption("RhoFluidMax", &rhoFluidMax);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	string dataDir1 = "Data-kt-mu-Solid";
	string dataDir2 = "Data-kt-mu-Solid-Fluid";
	
	int sysresult = system(("mkdir -p " + dataDir2).c_str());
	
	
	////////// Loop over all possible files in Solid directory  ////////////
	
	// not optimal, but fast enough
	
	for (int kTCounter=0; kTCounter<maxFileIndex; kTCounter++)
	{
		for (int muCounter=0; muCounter<maxFileIndex; muCounter++)
		{
			////////////////// Check if file exists ////////////////////
			
			ifstream dataFile(dataDir1+"/results_"+to_string(kTCounter)+
			                                   "_"+to_string(muCounter)+".dat");
			if (!dataFile) continue;
			
			//////////// Read solid computation parameters /////////////
			
			double kT = 0;
			double mu = 0;
			double aVdW = 0;
			double hsd  = 0;
			double densitySolid = 0;
			double numParticlesSolid = 0;
			double freeEnergySolid = 0;
			bool successSolid = false;
			int successSolid2 = 0;
			
			readDataFromFile(dataFile, "kT", kT);
			readDataFromFile(dataFile, "mu", mu);
			readDataFromFile(dataFile, "aVdW", aVdW);
			readDataFromFile(dataFile, "hsd", hsd);
			readDataFromFile(dataFile, "densitySolid", densitySolid);
			readDataFromFile(dataFile, "numParticlesSolid", numParticlesSolid);
			readDataFromFile(dataFile, "freeEnergySolid", freeEnergySolid);
			readDataFromFile(dataFile, "successSolid", successSolid2);
			successSolid = (successSolid2==1);
			
			dataFile.close();
			
			///////////////// Compute fluid properties /////////////////
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "kT = " << kT
			    << "mu = " << mu << endl;
			
			// properties of the most stable fluid phase
			double densityFluid = 0;
			double freeEnergyFluid = 0;
			
			// properties of the vapor and liquid phases independently
			double densityVapor = 0;
			double freeEnergyVapor = 0;
			double densityLiquid = 0;
			double freeEnergyLiquid = 0;
			
			bool successFluid = true;
			
			try
			{
				vector<double> roots;
				findRootsdOmegadRhoSpinodal(kT, mu, aVdW, hsd, 
					rhoFluidMin, rhoFluidMax, roots);
				
				if (roots.size()==1)
				{
					densityFluid = roots[0];
					freeEnergyFluid = uniformOmega(kT, mu, aVdW, 
						hsd, densityFluid);
					
					densityVapor = densityFluid;
					densityLiquid = densityFluid;
					freeEnergyVapor = freeEnergyFluid;
					freeEnergyLiquid = freeEnergyFluid;
				}
				else if (roots.size()==2)
				{
					densityVapor = roots[0];
					densityLiquid = roots[1];
					
					freeEnergyVapor  = uniformOmega(
						kT, mu, aVdW, hsd, roots[0]);
					freeEnergyLiquid = uniformOmega(
						kT, mu, aVdW, hsd, roots[1]);
					
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
				log << "Vapor Density = " << densityVapor << endl;
				log << "Vapor Free Energy = " << freeEnergyVapor << endl;
				log << "Liquid Density = " << densityLiquid << endl;
				log << "Liquid Free Energy = " << freeEnergyLiquid << endl;
			}
			else
			{
				log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
				log << "ERROR: Fluid computation FAILED"
				    << " at kT = " << kT << " mu = " << mu << endl;
			}
			
			/////////////////// Save results in file ///////////////////
			
			ofstream dataFile2(dataDir2+"/results_"+to_string(kTCounter)+
			                                    "_"+to_string(muCounter)+".dat");
			
			dataFile2 << "kT = " << kT << endl;
			dataFile2 << "mu = " << mu << endl;
			dataFile2 << endl;
			dataFile2 << "aVdW = " << aVdW << endl;
			dataFile2 << "hsd = " << hsd << endl;
			dataFile2 << endl;
			dataFile2 << "densitySolid = " << densitySolid << endl;
			dataFile2 << "numParticlesSolid = " << numParticlesSolid << endl;
			dataFile2 << "freeEnergySolid = " << freeEnergySolid << endl;
			dataFile2 << endl;
			dataFile2 << "densityFluid = " << densityFluid << endl;
			dataFile2 << "freeEnergyFluid = " << freeEnergyFluid << endl;
			dataFile2 << endl;
			dataFile2 << "densityVapor = " << densityVapor << endl;
			dataFile2 << "freeEnergyVapor = " << freeEnergyVapor << endl;
			dataFile2 << endl;
			dataFile2 << "densityLiquid = " << densityLiquid << endl;
			dataFile2 << "freeEnergyLiquid = " << freeEnergyLiquid << endl;
			dataFile2 << endl;
			dataFile2 << "successSolid = " << successSolid << endl;
			dataFile2 << "successFluid = " << successFluid << endl;
			dataFile2 << "success      = " << (successSolid && successFluid) << endl;
			
			dataFile2.close();
		}
	}
	
	return 1;
}

