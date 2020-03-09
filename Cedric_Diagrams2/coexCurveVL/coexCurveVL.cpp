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
	
	ofstream outDataFile("coexCurveVL.dat");
	outDataFile << "#kT				rhoV			rhoL			"
	            << "muCoex				freeEnergy		"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	double kT = kTMin;
	double kTStep = kTStepMax;
	
	while (kT<kTMax && kTStep>kTStepMin)
	{
		// Compute Hard sphere diameter and VdW parameter
		
		double hsd, aVdW;
		
		#ifdef POTENTIAL_WHDF
		log << "Using WHDF Potential" << endl;
		WHDF potential1(sigma1, eps1, rcut1);
		
		#else
		log << "Using LJ Potential" << endl;
		LJ potential1(sigma1, eps1, rcut1); // default
		#endif
		
		hsd = potential1.getHSD(kT);
		
		if (dx>0) // if we want to use the VdW parameter of a grid of spacing dx
		{
			// values doesn't matter here
			int Ngrid = 127;
			double L[3] = {Ngrid*dx, Ngrid*dx, Ngrid*dx};
			double mu = -6;
			
			SolidFCC theDensity1(dx, L, hsd);
			
			#ifdef ANALYTIC_WEIGHTS
			log << "Using analytic evaluation of weights" << endl;
			FMT_Species_Analytic species1(theDensity1,hsd, mu, Ngrid);
			Interaction_Linear_Interpolation i1(species1,species1,potential1,kT,log);
			i1.initialize();
			
			#else
			log << "Using numeric evaluation of weights" << endl;
			FMT_Species_Numeric species1(theDensity1,hsd,pointsFile, mu, Ngrid);
			Interaction i1(species1,species1,potential1,kT,log,pointsFile);
			i1.initialize();
			#endif
			
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
		
		// Also compute free energy
		
		double freeEnergyCoex;
		double densityCoex;
		
		Log log_freeEnergyCalc("log_freeEnergyCalc.dat");
		
		int statusFreeEnergyCalc = fixedkTMuFluid(kT, muCoex, argc, argv,
			log_freeEnergyCalc, freeEnergyCoex, densityCoex);
		
		// Save coexistence densities under critical point
		
		if (!superCritical && success && statusFreeEnergyCalc==0)
		{
			outDataFile << kT << " 	" << rhoV << " 	" << rhoL << " 	"
			            << muCoex << " 	" << freeEnergyCoex << " 	"
			            << endl;
		}
		else if (!superCritical && success)
		{
			outDataFile << kT << " 	" << rhoV << " 	" << rhoL << " 	"
			            << muCoex << " 	" << "failed" << " 	"
			            << endl;
		}
		
		// Report
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Coexistence properties at kT = " << kT << endl;
		log << "#" << endl;
		if (!superCritical && success)
		{
			log << "vapor density  = " << rhoV << endl;
			log << "liquid density = " << rhoL << endl;
			log << "mu = " << muCoex << endl;
		}
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
	outPlotFile << "set output 'coexCurveVL.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Vapor-Liquid\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key bottom left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveVL.dat' using 2:1, \\" << endl;
	outPlotFile << "     'coexCurveVL.dat' using 3:1" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}






