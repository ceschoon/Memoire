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
#include "UniformFluid.h"
#include "UniformFluid.cpp"




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to find the coexistence curve between the vapor and liquid phases

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// thermodynamics
	double kTMin = 1;
	double kTMax = 1.5;
	double kTStep = 0.05;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	/////////// Read options
	
	Options options;
	
	// thermodynamics
	options.addOption("kTMin", &kTMin);
	options.addOption("kTMax", &kTMax);
	options.addOption("kTStep", &kTStep);
	
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
	
	////////////////// Coexistence properties for each kT //////////////////
	
	ofstream outDataFile("coexistence_uniform.dat");
	outDataFile << "#kT #muCoex #rho1 #rho2" << endl;
	
	int sizekT = int( (kTMax-kTMin)/kTStep ) +1;
	vector<double> kT(sizekT);
	for (int i=0; i<sizekT; i++) kT[i] = kTMin + i*kTStep;
	
	for (int i=0; i<sizekT; i++)
	{
		// Compute Hard sphere diameter and VdW parameter
		
		LJ potential1(sigma1, eps1, rcut1);
		double hsd = potential1.getHSD(kT[i]);
		double aVdW = potential1.getVDW_Parameter(kT[i]); // what about the eps=1 ?
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "kT = " << kT[i] << endl;
		log << "hsd = " << hsd << endl;
		log << "aVdW = " << aVdW << endl;
		
		// Find coexistence densities and mu
		
		double rho1,rho2,muCoex;
		bool critical = false;
		bool success = true;
		
		coexistenceDensities(argc, argv, log, kT[i], aVdW, hsd, rho1, rho2,
		                     muCoex, critical, success);
		
		if (!critical && success)
			outDataFile << kT[i] << " " << rho1 << " " << rho2 << " " 
			            << muCoex << endl;
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Coexistence properties at kT = " << kT[i] << endl;
		log << "#" << endl;
		log << "vapor density  = " << rho1 << endl;
		log << "liquid density = " << rho2 << endl;
		log << "mu = " << muCoex << endl;
		log << "critical = " << critical << endl;
		log << "success = " << success << endl;
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("gnuplot_plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexistence_uniform.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Vapor-Liquid\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexistence_uniform.dat' using 2:1, \\" << endl;
	outPlotFile << "     'coexistence_uniform.dat' using 3:1" << endl;
	
	return 1;
}






