////////////////////////////////////////////////////////////////////////////
//                                                                        //
//      Program to compute physical properties at the critical point      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>

using namespace std;

#include "options.h"

#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "../UniformPhases.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to plot typical free energy profiles for the uniform phases

int main(int argc, char** argv)
{
	/////////////////////////////// Options ////////////////////////////////
	
	string dataDir = "data";
	double kT = 1.0;
	double mu = -3.0;
	double rhoMin = 0.01;
	double rhoMax = 0.6;
	double rhoStep = 0.001;
	
	Options options;
	
	options.addOption("dataDir", &dataDir);
	options.addOption("kT", &kT);
	options.addOption("mu", &mu);
	options.addOption("RhoMinPlot", &rhoMin);
	options.addOption("RhoMaxPlot", &rhoMax);
	options.addOption("RhoStepPlot", &rhoStep);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	ofstream dataFile(dataDir+"/omega_uniform.dat");
	dataFile << scientific << setprecision(8);
	dataFile << "# Uniform free energy profile and its derivatives" << endl;
	dataFile << "# kT = " << kT << " mu = " << mu << endl;
	dataFile << "# rho         omega         domega        d2omega" << endl;
	
	
	////////////////// Compute hsd and aVdW parameters /////////////////////
	
	log << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Computing hsd and VdW parameter" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	double hsd, aVdW;
	compute_hsd_aVdW(argc, argv, log, kT, hsd, aVdW);
	
	
	///////////// Compute free energy profile and derivatives //////////////
	
	log << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Computing free energy profile and derivatives" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << endl;
	
	for (double rho=rhoMin; rho<rhoMax; rho+=rhoStep)
	{
		log << "currently at rho = " << rho << endl;
		
		double omega = uniformOmega(kT, mu, aVdW, hsd, rho);
		double domega = uniformOmegaDerivative(kT, mu, aVdW, hsd, rho);
		double d2omega = uniformOmegaDerivative2(kT, mu, aVdW, hsd, rho);
		
		dataFile << rho << " " << omega << " " << domega << " " 
		         << d2omega << endl; 
	}
	
	
	////////////////////// Plot free energy profile ////////////////////////
	
	
	log << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Making plot scripts" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	
	string fileName = "omega";
	ofstream plotFile(dataDir+"/"+fileName+".gp");
	
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output '" << fileName << ".tex'" << endl;
	plotFile << endl;
	plotFile << "#set title 'Free energy profile'" << endl;
	plotFile << "set title '\\Large$kT = " << kT << "$, \\Large$\\beta\\mu = " << mu << "$'" << endl;
	plotFile << "set xlabel '$\\rho$'" << endl;
	plotFile << "set ylabel '$\\beta\\omega$'" << endl;
	plotFile << "set grid linecolor rgb '#B0B0B0'" << endl;
	plotFile << "set key off" << endl;
	plotFile << "set xrange [0:0.6]" << endl;
	plotFile << "set yrange [-0.12:0]" << endl;
	plotFile << endl;
	plotFile << "plot 'omega_uniform.dat' using 1:2 with lines" << endl;
	
	plotFile.close();
	
	/*
	sysresult = system(("cd "+dataDir+"; gnuplot "+fileName+".gp").c_str());
	sysresult = system(("cd "+dataDir+"; pdflatex "+fileName+".tex").c_str());
	sysresult = system(("cd "+dataDir+"; rm "+fileName+".aux "+fileName+
	                    ".log "+fileName+"-inc.eps "+fileName+
	                    "-inc-eps-converted-to.pdf "+fileName+".tex").c_str());
	*/
	
	////////////////////// Plot first derivative ///////////////////////////
	
	
	fileName = "domega";
	plotFile.open(dataDir+"/"+fileName+".gp");
	
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output '" << fileName << ".tex'" << endl;
	plotFile << endl;
	plotFile << "#set title 'First derivative'" << endl;
	plotFile << "set title '\\Large$kT = " << kT << "$, \\Large$\\beta\\mu = " << mu << "$'" << endl;
	plotFile << "set xlabel '\\Large$\\rho$'" << endl;
	plotFile << "set ylabel '\\Large$\\beta\\partial\\omega/\\partial\\rho$'" << endl;
	plotFile << "set grid linecolor rgb '#B0B0B0'" << endl;
	plotFile << "set key off" << endl;
	plotFile << "set xrange [0:0.6]" << endl;
	plotFile << "set yrange [-0.25:0.25]" << endl;
	plotFile << endl;
	plotFile << "plot 'omega_uniform.dat' using 1:3 with lines,\\" << endl;
	plotFile << "     'omega_uniform.dat' using 1:(0) with lines linecolor 'black'" << endl;
	
	plotFile.close();
	
	/*
	sysresult = system(("cd "+dataDir+"; gnuplot "+fileName+".gp").c_str());
	sysresult = system(("cd "+dataDir+"; pdflatex "+fileName+".tex").c_str());
	sysresult = system(("cd "+dataDir+"; rm "+fileName+".aux "+fileName+
	                    ".log "+fileName+"-inc.eps "+fileName+
	                    "-inc-eps-converted-to.pdf "+fileName+".tex").c_str());
	*/
	
	////////////////////// Plot second derivative //////////////////////////
	
	
	fileName = "d2omega";
	plotFile.open(dataDir+"/"+fileName+".gp");
	
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output '" << fileName << ".tex'" << endl;
	plotFile << endl;
	plotFile << "#set title 'Second derivative'" << endl;
	plotFile << "set title '\\Large$kT = " << kT << "$, \\Large$\\beta\\mu = " << mu << "$'" << endl;
	plotFile << "set xlabel '\\Large$\\rho$'" << endl;
	plotFile << "set ylabel '\\Large$\\beta\\partial^2\\omega/\\partial\\rho^2$'" << endl;
	plotFile << "set grid linecolor rgb '#B0B0B0'" << endl;
	plotFile << "set key off" << endl;
	plotFile << "set xrange [0:0.6]" << endl;
	plotFile << "set yrange [-2:3]" << endl;
	plotFile << endl;
	plotFile << "plot 'omega_uniform.dat' using 1:4 with lines,\\" << endl;
	plotFile << "     'omega_uniform.dat' using 1:(0) with lines linecolor 'black'" << endl;
	
	plotFile.close();
	
	/*
	sysresult = system(("cd "+dataDir+"; gnuplot "+fileName+".gp").c_str());
	sysresult = system(("cd "+dataDir+"; pdflatex "+fileName+".tex").c_str());
	sysresult = system(("cd "+dataDir+"; rm "+fileName+".aux "+fileName+
	                    ".log "+fileName+"-inc.eps "+fileName+
	                    "-inc-eps-converted-to.pdf "+fileName+".tex").c_str());
	*/
	
	//////////////////// Bash script to plot everything ////////////////////
	
	
	ofstream shellFile(dataDir+"/plot_everything.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   omega.gp" << endl;
	shellFile << "gnuplot  domega.gp" << endl;
	shellFile << "gnuplot d2omega.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex   omega.tex" << endl;
	shellFile << "pdflatex  domega.tex" << endl;
	shellFile << "pdflatex d2omega.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   omega.aux   omega.log   omega.tex   omega-inc.eps   omega-inc-eps-converted-to.pdf" << endl;
	shellFile << "rm  domega.aux  domega.log  domega.tex  domega-inc.eps  domega-inc-eps-converted-to.pdf" << endl;
	shellFile << "rm d2omega.aux d2omega.log d2omega.tex d2omega-inc.eps d2omega-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	shellFile << "if [ $# -gt 0 ]" << endl;
	shellFile << "then" << endl;
	shellFile << "mv   omega.pdf   omega_$1.pdf" << endl;
	shellFile << "mv  domega.pdf  domega_$1.pdf" << endl;
	shellFile << "mv d2omega.pdf d2omega_$1.pdf" << endl;
	shellFile << "fi" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	sysresult = system(("chmod +x "+dataDir+"/plot_everything.sh; cd "+dataDir+"; ./plot_everything.sh").c_str());
	
	return 0;
}






