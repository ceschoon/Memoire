#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace std;

#include "options.h"

//#include "Potential1.h"
#include "Log.h"
#include "myColor.h"
#include "../UniformPhases.h"



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
	
	// root min and max boundaries
	double rhoMin = 0.001;
	double rhoMax = 1.5;
	
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
	options.addOption("RhoMin", &rhoMin);
	options.addOption("RhoMax",  &rhoMax);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	///////////
	
	// Compute Hard sphere diameter and VdW parameter
	
	/* // TODO: there is currently a problem in the potential class
	LJ potential1(sigma1, eps1, rcut1);
	double hsd = potential1.getHSD(kT);
	double aVdW = potential1.getVDW_Parameter(kT); 
	*/
	
	// manually set the parameters
	double hsd = 1.01561;
	double aVdW = -7.27846;
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "hsd = " << hsd << endl;
	log << "aVdW = " << aVdW << endl;
	
	///////////////////// Find the roots of dOmega/drho ////////////////////
	
	vector<double> roots;
	vector<double> omegaAtroots;
	
	findRootsdOmegadRhoSpinodal(kT, mu, aVdW, hsd, rhoMin, rhoMax, roots);
	
	for (int i=0; i<roots.size(); i++)
	{
		double omega = uniformOmega(kT, mu, aVdW, hsd, roots[i] );
		omegaAtroots.push_back(omega);
	}
	
	// Report
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "number of roots  = " << roots.size() << endl;
	log << "roots = ";
	for (int i=0; i<roots.size(); i++)
		log << roots[i] << " ";
	log << endl;
	log << "omega at roots = ";
	for (int i=0; i<omegaAtroots.size(); i++)
		log << omegaAtroots[i] << " ";
	log << endl;
	
	return 1;
}
