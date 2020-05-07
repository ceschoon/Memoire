#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace std;

#include "options.h"

#include "Potential1.h"
#include "Log.h"
#include "myColor.h"
#include "../UniformFluid.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the density minima of the uniform fluid free energy
// as well the corresponding free energies by finding the roots of the 
// free energy derivative

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// thermodynamics
	double kT = 1;
	double mu = -3;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	// root initial guesses
	double rhoLowInit = 0.01;
	double rhoHighInit = 1;
	
	/////////// Read options
	
	Options options;
	
	// thermodynamics
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	// root initial guesses
	options.addOption("RhoLowInit", &rhoLowInit);
	options.addOption("RhoHighInit",  &rhoHighInit);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	///////////
	
	// Compute Hard sphere diameter and VdW parameter
	
	LJ potential1(sigma1, eps1, rcut1);
	double hsd = potential1.getHSD(kT);
	double aVdW = potential1.getVDW_Parameter(kT); // TODO: what about the eps=1 ?
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "hsd = " << hsd << endl;
	log << "aVdW = " << aVdW << endl;
	
	//////////////////////// Find the density minimum //////////////////////
	// we do this by finding root of the potential derivative
	
	bool success = true;
	
	double rhoLow  = findRootdOmegadRhoNewton(rhoLowInit , kT, mu, aVdW, hsd, success);
	double rhoHigh = findRootdOmegadRhoNewton(rhoHighInit, kT, mu, aVdW, hsd, success);
	
	double OmegaLow  = uniformFluidOmega(kT, mu, aVdW, hsd, rhoLow );
	double OmegaHigh = uniformFluidOmega(kT, mu, aVdW, hsd, rhoHigh);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "rhoLow  = " << rhoLow << endl;
	log << "rhoHigh = " << rhoHigh << endl;
	log << "OmegaLow = " << OmegaLow << endl;
	log << "OmegaHigh = " << OmegaHigh << endl;
	log << "success = " << success << endl;
	
	return 1;
}
