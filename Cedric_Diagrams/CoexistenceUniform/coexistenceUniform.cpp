#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <stdlib.h>

using namespace std;

#include "options.h"

#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "../UniformPhases.h"
#include "../SolidFCC.h"




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to find the coexistence curve between the vapor and liquid phases

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// geometry
	double dx = 0.0125;
	string pointsFile;
	
	// thermodynamics
	double kTMin = 1;
	double kTMax = 1.5;
	double kTStepMax = 0.1;
	double kTStepMin = 0.01;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	/////////// Read options
	
	Options options;
	
	// geometry
	options.addOption("dx", &dx);
	options.addOption("IntegrationPointsFile", &pointsFile);
	
	// thermodynamics
	options.addOption("kTMin", &kTMin);
	options.addOption("kTMax", &kTMax);
	options.addOption("kTStepMax", &kTStepMax);
	options.addOption("kTStepMin", &kTStepMin);
	
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
	outDataFile << "#kT #rhoV #rhoL #muCoex" << endl;
	
	double kT = kTMin;
	double kTStep = kTStepMax;
	
	while (kT<kTMax && kTStep>kTStepMin)
	{
		// Compute Hard sphere diameter and VdW parameter
		
		double hsd, aVdW;
		
		LJ potential1(sigma1, eps1, rcut1);
		hsd = potential1.getHSD(kT);
		
		if (dx>0) // if we want to use the VdW parameter of a grid of spacing dx
		{
			// values doesn't matter here
			int Npoints = 127;
			double L[3] = {Npoints*dx, Npoints*dx, Npoints*dx};
			double mu = -6;
			
			SolidFCC theDensity1(dx, L, hsd);
			FMT_Species_Numeric species1(theDensity1,hsd,pointsFile, mu, Npoints);
			Interaction i1(species1,species1,potential1,kT,log,pointsFile);
			i1.initialize();
			aVdW = i1.getVDWParameter()/2; // divide by 2 as not the same definition
		}
		else
		{
			aVdW = potential1.getVDW_Parameter(kT);
		}
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "kT = " << kT << endl;
		log << "hsd = " << hsd << endl;
		log << "aVdW = " << aVdW << endl;
		
		// Find coexistence densities and mu
		
		double rhoV,rhoL,muCoex;
		bool superCritical = false;
		bool success = true;
		
		coexistenceDensities(argc, argv, log, kT, aVdW, hsd, rhoV, rhoL,
		                     muCoex, superCritical, success);
		
		// Save coexistence densities under critical point
		
		if (!superCritical && success)
			outDataFile << kT << " " << rhoV << " " << rhoL << " " 
			            << muCoex << endl;
		
		// Report
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Coexistence properties at kT = " << kT << endl;
		log << "#" << endl;
		log << "vapor density  = " << rhoV << endl;
		log << "liquid density = " << rhoL << endl;
		log << "mu = " << muCoex << endl;
		log << "superCritical = " << superCritical << endl;
		log << "success = " << success << endl;
		
		// Prepare next step
		
		if (superCritical)
		{
			kT -= kTStep/2;
			kTStep = kTStep/2;
			
			log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Current best estimate for the critical temperature is " 
			    << kT << endl;
			log << "Current best estimate for the critical density is     " 
			    << (rhoV+rhoL)/2 << endl;
		}
		else
		{
			kT += kTStep;
		}
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("plot");
	
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
	outPlotFile << "set key bottom left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexistence_uniform.dat' using 2:1, \\" << endl;
	outPlotFile << "     'coexistence_uniform.dat' using 3:1" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}






