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

#include <gsl/gsl_sf_lambert.h>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "UniformFluid.h"
#include "UniformFluid.cpp"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to plot the polynomial function used to find the coexistence curve

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// thermodynamics
	double kT = 1;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	// first guess at density
	double rho1 = 0.04;
	
	/////////// Read options
	
	Options options;
	
	// thermodynamics
	options.addOption("kT", &kT);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	// first guess at density
	options.addOption("rho1",  &rho1);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	////////////////// Compute polynomial parameters ///////////////////////
	
	// Compute Hard sphere diameter and VdW parameter
	
	LJ potential1(sigma1, eps1, rcut1);
	double hsd = potential1.getHSD(kT);
	double aVdW = potential1.getVDW_Parameter(kT); // what about the eps=1 ?
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "hsd = " << hsd << endl;
	log << "aVdW = " << aVdW << endl;
	
	// Compute polynomial parameters for the coexistence equation
	// omega - rho*domega/drho = omega1
	
	double mu1 = muFromCoexDensity(rho1, kT, aVdW, hsd);
	double omega1 = uniformFluidOmega(kT, mu1, aVdW, hsd, rho1);
	double omegaTilde = omega1 / kT * 4*M_PI/3*pow(hsd/2,3);
	double aVdWTilde = aVdW / (4*M_PI/3*pow(hsd/2,3));
	
	log << "omegaTilde = " << omegaTilde << endl;
	log << "aVdWTilde= " << aVdWTilde << endl;
	
	/////////////////////////// Evaluate and save //////////////////////////
	
	ofstream dataFile("polyEta.dat");
	
	double rhoMin = 0;
	double rhoMax = 1.3;
	double rhoStep = 0.01;
	
	double rho = rhoMin;
	
	while (rho < rhoMax)
	{
		double etaFMT = 4*M_PI/3*pow(hsd/2,3) * rho;
		
		double polyEta =                     aVdWTilde * pow(etaFMT,5)
		                 +             (1-3*aVdWTilde) * pow(etaFMT,4)
		                 +  (omegaTilde-1+3*aVdWTilde) * pow(etaFMT,3)
		                 + (-1-aVdWTilde-3*omegaTilde) * pow(etaFMT,2)
		                 +           (-1+3*omegaTilde) * etaFMT
		                 -                 omegaTilde  * 1;
		
		dataFile << rho << " " << polyEta << endl;
		
		rho += rhoStep;
	}
	
	//////////////////////// Create gnuplot script /////////////////////////
	
	ofstream outPlotFile("gnuplot_plot");
	
	outPlotFile << "########## Plot coexistence polynomial in etaFMT #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'polyEta.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence polynomial\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"poly\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key off" << endl;
	outPlotFile << "plot 'polyEta.dat'" << endl;
	outPlotFile << endl;
	
	return 1;
}

