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

// Computation at given (kT,mu,Ngrid,Cvac,alpha)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDir = "data";
	double kT = 1;
	double mu = -1;
	int    Ngrid = 125;
	double Cvac = 1e-4;
	double alpha = 100;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	options.addOption("Ngrid", &Ngrid);
	options.addOption("Cvac", &Cvac);
	options.addOption("alpha", &alpha);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	////////////////// Compute free energy for the solid ///////////////////
	
	double freeEnergy, freeEnergyErr, density;
	
	int status = DFTgaussian( argc, argv, log, kT, mu, Ngrid, Cvac, alpha,
	                          freeEnergy, freeEnergyErr, density );
	
	// Print results in log file
	
	if (status==0)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Free Energy = " << freeEnergy << endl;
		log << "Free Energy (Error) = " << freeEnergyErr << endl;
		log << "Density = " << density << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Solid Computation FAILED" << endl;
	}
	
	return 1;
}
