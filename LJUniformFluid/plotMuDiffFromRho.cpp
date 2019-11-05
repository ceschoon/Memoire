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
	
	/////////// Read options
	
	Options options;
	
	// thermodynamics
	options.addOption("kT", &kT);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// Compute Hard sphere diameter and VdW parameter
	
	LJ potential1(sigma1, eps1, rcut1);
	double hsd = potential1.getHSD(kT);
	double aVdW = potential1.getVDW_Parameter(kT); // what about the eps=1 ?
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "hsd = " << hsd << endl;
	log << "aVdW = " << aVdW << endl;
	
	/////////////////////////// Evaluate and save //////////////////////////
	
	ofstream dataFile("muDiff.dat");
	
	double rhoMin = 0.01;
	double rhoMax = 0.8;
	double rhoStep = 0.01;
	
	double rho = rhoMin;
	
	while (rho < rhoMax)
	{
		cout << "rho = " << rho << endl;
		
		double rho2 = rho2FromRho1(rho, kT, aVdW, hsd);
		double mu1 = muFromCoexDensity(rho , kT, aVdW, hsd);
		double mu2 = muFromCoexDensity(rho2, kT, aVdW, hsd);
		
		dataFile << rho << " " << mu2-mu1 << endl;
		
		rho += rhoStep;
	}
	
	//////////////////////// Create gnuplot script /////////////////////////
	
	ofstream outPlotFile("gnuplot_plot");
	
	outPlotFile << "########## Plot mu difference #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'polyEta.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Mu Difference for coexistence densities candidates\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho1\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"mu2-mu1\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key off" << endl;
	outPlotFile << "plot 'muDiff.dat'" << endl;
	outPlotFile << endl;
	
	return 1;
}

