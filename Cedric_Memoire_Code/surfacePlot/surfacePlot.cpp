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

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"
#include "DFT.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

#include "../SolidPhase.h"
#include "../UniformPhases.h"


int surfacePlot( double kT, double mu, int Ngrid,
                 int argc, char** argv, Log &log);



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the free energies at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT_rangeMin = 1;
	double kT_rangeMax = 1;
	double kT_rangeStep = 1;
	
	double mu_rangeMin = -1;
	double mu_rangeMax = -1;
	double mu_rangeStep = 1;
	
	int Ngrid_rangeMin = 125;
	int Ngrid_rangeMax = 125;
	int Ngrid_rangeStep = 1;
	
	Options options;
	
	options.addOption("kT_rangeMin", &kT_rangeMin);
	options.addOption("kT_rangeMax", &kT_rangeMax);
	options.addOption("kT_rangeStep", &kT_rangeStep);
	
	options.addOption("mu_rangeMin", &mu_rangeMin);
	options.addOption("mu_rangeMax", &mu_rangeMax);
	options.addOption("mu_rangeStep", &mu_rangeStep);
	
	options.addOption("Ngrid_rangeMin", &Ngrid_rangeMin);
	options.addOption("Ngrid_rangeMax", &Ngrid_rangeMax);
	options.addOption("Ngrid_rangeStep", &Ngrid_rangeStep);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////////////// Loop over variables //////////////////////////
	
	for (double kT=kT_rangeMin; kT<=kT_rangeMax; kT+=kT_rangeStep)
		for (double mu=mu_rangeMin; mu<=mu_rangeMax; mu+=mu_rangeStep)
			for (int Ngrid=Ngrid_rangeMin; Ngrid<=Ngrid_rangeMax; Ngrid+=Ngrid_rangeStep)
				surfacePlot( kT, mu, Ngrid, argc, argv, log);
	
	return 0;
}











////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



// Function to compute the freeEnergy on the (Cvac,alpha) plane in order
// to plot a heatmap.

int surfacePlot( double kT, double mu, int Ngrid,
                 int argc, char** argv, Log &log)
{
	////////////////////////////////// Log /////////////////////////////////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Surface plot on (Cvac,alpha) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// Seperate log file for the DFT computation
	// by default, we don't want the non-interesting log from the DFT to 
	// pollute the current log
	
	Log logDFT("log_DFT.dat");
	
	///////////////////////////// Read options /////////////////////////////
	
	string dataDir = "data";
	double Cvac_rangeMin;
	double Cvac_rangeMax;
	double Cvac_logStep;
	double alpha_rangeMin;
	double alpha_rangeMax;
	double alpha_logStep;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("Cvac_rangeMin", &Cvac_rangeMin);
	options.addOption("Cvac_rangeMax", &Cvac_rangeMax);
	options.addOption("Cvac_logStep", &Cvac_logStep);
	options.addOption("alpha_rangeMin", &alpha_rangeMin);
	options.addOption("alpha_rangeMax", &alpha_rangeMax);
	options.addOption("alpha_logStep", &alpha_logStep);
	
	options.read(argc, argv);
	
	/////////////////////////////// Init file //////////////////////////////
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	stringstream ssNgrid; ssNgrid << Ngrid;
	
	string surfacePlotFileName = "";
	surfacePlotFileName = surfacePlotFileName + "surfacePlot_"
	                    + "kT="+sskT.str() + "_mu="+ssMu.str()
	                    + "_Ngrid="+ssNgrid.str();
	ofstream surfacePlotFile(dataDir+"/"+surfacePlotFileName+".dat");
	
	surfacePlotFile << "# Parameters:" << endl;
	surfacePlotFile << "# kT = " << kT << endl;
	surfacePlotFile << "# mu = " << mu << endl;
	surfacePlotFile << "# Ngrid = " << Ngrid << endl;
	surfacePlotFile << "# " << endl;
	surfacePlotFile << "# Plot:" << endl;
	surfacePlotFile << "# log10Cvac		log10alpha		freeEnergy	"
	                << endl;
	surfacePlotFile << "# " << endl;
	surfacePlotFile << scientific;
	
	//////////////// Compute free energy on (Cvac,alpha) plane /////////////
	
	for (double Cvac=Cvac_rangeMin; Cvac<Cvac_rangeMax; Cvac*=exp(Cvac_logStep))
	{
		for (double alpha=alpha_rangeMin; alpha<alpha_rangeMax; alpha*=exp(alpha_logStep))
		{
			log << endl;
			log << "Computing for Cvac = " << Cvac << " and alpha = " 
			    << alpha << endl;
			
			double freeEnergy, density;
			int statusDFT = DFTgaussian( argc, argv, logDFT, kT, mu, Ngrid, 
				Cvac, alpha, freeEnergy, density);
			
			if (statusDFT==0)
			{
				log << "freeEnergy = " << freeEnergy << endl;
				surfacePlotFile << std::log(Cvac)  / std::log(10) << " 	"
				                << std::log(alpha) / std::log(10) << " 	"
				                << freeEnergy << endl;
			}
			else
			{
				log << "freeEnergy computation failed" << endl;
				surfacePlotFile << std::log(Cvac)  / std::log(10) << " 	"
				                << std::log(alpha) / std::log(10) << " 	"
				                << "failed" << endl;
			}
		}
		surfacePlotFile << endl; // seperate data blocks
	}
	
	/////////////////////////// Plot with gnuplot //////////////////////////
	
	double log10Cvac_rangeMin_plot = std::log( Cvac_rangeMin )/std::log(10);
	double log10Cvac_rangeMax_plot = std::log( Cvac_rangeMax * exp(-Cvac_logStep) )/std::log(10);
	double log10alpha_rangeMin_plot = std::log( alpha_rangeMin )/std::log(10);
	double log10alpha_rangeMax_plot = std::log( alpha_rangeMax * exp(-alpha_logStep) )/std::log(10);
	
	ofstream plotFile(dataDir+"/surfacePlot.gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color" << endl;
	plotFile << "set output 'surfacePlot.tex'" << endl;
	plotFile << endl;
	plotFile << "set title '$kT = " << kT << " \\ \\beta\\mu = " << mu << " \\ N_{grid} = " << Ngrid << "$'" << endl;
	plotFile << "set xlabel '$\\log_{10} c_{vac}$'" << endl;
	plotFile << "set ylabel '$\\log_{10}\\alpha$'" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << "set xrange[" << log10Cvac_rangeMin_plot  << ":" << log10Cvac_rangeMax_plot  << "]" << endl;
	plotFile << "set yrange[" << log10alpha_rangeMin_plot << ":" << log10alpha_rangeMax_plot << "]" << endl;
	plotFile << endl;
	plotFile << "set datafile missing 'failed'" << endl;
	//plotFile << "set view map" << endl;
	//plotFile << "set dgrid3d" << endl;
	//plotFile << "set pm3d interpolate 10,10" << endl;
	//plotFile << "splot \"" << surfacePlotFileName << ".dat\" using 1:2:3 with pm3d" << endl;
	plotFile << "plot '" << surfacePlotFileName << ".dat' using 1:2:($3/1.0) with image" << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	
	ofstream shellFile(dataDir+"/plot_everything.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   surfacePlot.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  surfacePlot.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   surfacePlot.aux   surfacePlot.log   surfacePlot.tex   surfacePlot-inc.eps   surfacePlot-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	sysresult = system(("chmod +x "+dataDir+"/plot_everything.sh; cd "+dataDir+"; ./plot_everything.sh").c_str());
	
	// save the plot image with a specific name
	sysresult = system(("cd "+dataDir+"; mv surfacePlot.pdf "+surfacePlotFileName+".pdf").c_str());
	
	return 0;
}

