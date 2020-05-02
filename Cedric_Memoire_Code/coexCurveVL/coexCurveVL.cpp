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




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to find the coexistence curve between the vapor and liquid phases

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDir = "data";
	double kTMin = 1;
	double kTMax = 1.5;
	double kTStepMax = 0.1;
	double kTStepMin = 0.01;
	
	/////////// Read options
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kTMin", &kTMin);
	options.addOption("kTMax", &kTMax);
	options.addOption("kTStepMax", &kTStepMax);
	options.addOption("kTStepMin", &kTStepMin);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	////////////////// Coexistence properties for each kT //////////////////
	
	ofstream dataFile(dataDir+"/coexCurveVL.dat");
	dataFile << "#kT				rhoV			rhoL			"
	         << "muCoex				freeEnergy		"
	         << endl;
	dataFile << scientific << setprecision(8);
	
	double kT = kTMin;
	double kTStep = kTStepMax;
	
	while (kT<kTMax && kTStep>kTStepMin)
	{
		// Find coexistence densities and mu
		
		double rhoV,rhoL,mu,omega;
		bool superCritical = false;
		
		int status = vapourLiquidCoexistence(argc, argv, log, kT, rhoV, rhoL, 
			mu, omega, superCritical);
		
		// Save coexistence densities under critical point
		
		if (!superCritical && status==0)
		{
			dataFile << kT << " 	" << rhoV << " 	" << rhoL << " 	"
			         << mu << " 	" << omega << " 	"
			         << endl;
		}
		
		// Report
		
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Coexistence properties at kT = " << kT << endl;
		log << "#" << endl;
		if (!superCritical && status==0)
		{
			log << "vapor density  = " << rhoV << endl;
			log << "liquid density = " << rhoL << endl;
			log << "mu = " << mu << endl;
		}
		log << "superCritical = " << superCritical << endl;
		log << "success = " << (status==0) << endl;
		
		// Prepare next step
		
		if (superCritical)
		{
			kT -= kTStep/2;
			kTStep = kTStep/2;
			
			log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Current best estimate for the critical temperature is " 
			    << kT << endl;
		}
		else
		{
			kT += kTStep;
		}
	}
	
	dataFile.close();
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream plotFile(dataDir+"/coexCurveVL.gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output 'coexCurveVL.tex'" << endl;
	plotFile << endl;
	plotFile << "set title 'Coexistence Vapour-Liquid'" << endl;
	plotFile << "set xlabel '$\\rho$'" << endl;
	plotFile << "set ylabel '$kT$'" << endl;
	plotFile << "set style data lines" << endl;
	plotFile << "set grid" << endl;
	plotFile << endl;
	plotFile << "set key bottom left" << endl;
	plotFile << endl;
	plotFile << "plot 'coexCurveVL.dat' using 2:1 title 'vapour', \\" << endl;
	plotFile << "     'coexCurveVL.dat' using 3:1 title 'liquid'    " << endl;
	plotFile << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	ofstream shellFile(dataDir+"/plot_everything.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   coexCurveVL.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  coexCurveVL.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   coexCurveVL.aux   coexCurveVL.log   coexCurveVL.tex   coexCurveVL-inc.eps   coexCurveVL-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	sysresult = system(("chmod +x "+dataDir+"/plot_everything.sh; cd "+dataDir+"; ./plot_everything.sh").c_str());
	
	
	return 0;
}






