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
#include <armadillo>

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

// Program to plot the Omega(rho) curve for the uniform liquid

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// thermodynamics
	double kT = 1;
	double muMin = -1.5;
	double muMax = -0.5;
	double muStep = 0.5;
	double rhoMin = 0.01;
	double rhoMax = 0.80;
	double rhoStep = 0.01;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	// plot details
	double yRangeMin = 0;
	double yRangeMax = 0;
	
	Options options;
	
	// thermodynamics
	options.addOption("kT", &kT);
	options.addOption("MuMin", &muMin);
	options.addOption("MuMax", &muMax);
	options.addOption("MuStep", &muStep);
	options.addOption("RhoMin", &rhoMin);
	options.addOption("RhoMax", &rhoMax);
	options.addOption("RhoStep", &rhoStep);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	// plot details
	options.addOption("yRangeMin", &yRangeMin);
	options.addOption("yRangeMax", &yRangeMax);
	
	options.read(argc, argv);
	
	/////////////////////// Multiple plots, loop over mu ///////////////////
	
	// Compute Hard sphere diameter and VdW parameter
	
	LJ potential1(sigma1, eps1, rcut1);
	double hsd1 = potential1.getHSD(kT);
	double aVdW = potential1.getVDW_Parameter(kT); // TODO: what about the eps=1 ?
	
	// create gnuplot script
	
	ofstream outPlotFile("gnuplot_plot");
	
	outPlotFile << "########## Plot Uniform Free Energy #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'Omega.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Uniform Free Energy\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"Omega\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	if (yRangeMin!=0 && yRangeMin!=0) 
		outPlotFile << "set yrange ["+to_string(yRangeMin)+":"+to_string(yRangeMax)+"]" << endl;
	outPlotFile << endl;
		
	// loop over mu
	
	int sizeMu = int( (muMax-muMin)/muStep +1 );
	vector<double> mu(sizeMu);
	
	for (int i=0; i<sizeMu; i++)
		mu[i] = muMin + muStep*i;
	
	for (int j=0; j<sizeMu; j++)
	{
		// Compute the function
		
		int sizeRho = int( (rhoMax-rhoMin)/rhoStep +1 );
		vector<double> rho(sizeRho);
		vector<double> Omega(sizeRho);
		
		for (int i=0; i<sizeRho; i++) 
		{
			rho[i] = rhoMin + rhoStep*i;
			Omega[i] = uniformFluidOmega(kT,mu[j],aVdW,hsd1,rho[i]);
		}
		
		// Save data to plot
		
		ofstream outDataFile("Omega_"+to_string(j)+".dat");
		
		for (int i=0; i<sizeRho; i++)
			outDataFile << rho[i] << " " << Omega[i] << endl;
		
		if (sizeMu==1)
			outPlotFile << "plot 'Omega_"+to_string(j)+".dat' title 'kT="+to_string(kT)+",mu="+to_string(mu[j])+"'" << endl;
		else if (j==0)
			outPlotFile << "plot 'Omega_"+to_string(j)+".dat' title 'kT="+to_string(kT)+",mu="+to_string(mu[j])+"', \\" << endl;
		else if (j==sizeMu-1)
			outPlotFile << "     'Omega_"+to_string(j)+".dat' title 'kT="+to_string(kT)+",mu="+to_string(mu[j])+"'" << endl;
		else 
			outPlotFile << "     'Omega_"+to_string(j)+".dat' title 'kT="+to_string(kT)+",mu="+to_string(mu[j])+"', \\" << endl;
		
	}
	
	return 1;
}


