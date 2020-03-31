#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_multifit.h>

using namespace std;

#include "options.h"

#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "../UniformPhases.h"
#include "../SolidFCC.h"
#include "../utilities.h"



int coexVL(int argc, char** argv, double &kT_c, double &rho_c);

int myPolyFit(vector<double> data_x, vector<double> data_y,
              vector<double> &coefficients, int orderPolyFit);



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the critical exponent beta defined by
// rhoL-rhoV ~ (kTc-kT)^beta

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// none
	Log log("log_exponent.dat");
	
	///////////////////// Compute VL Coexistence Curve /////////////////////
	
	double kT_c,rho_c;
	int status_coex = coexVL(argc, argv, kT_c, rho_c);
	
	if (status_coex != 0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Coexistence computation not successful" << endl;
		return 1;
	}
	
	/////////////////// Get and tranform data for the fit //////////////////
	
	ifstream coexDataFile("coexCurveVL.dat");
	if (!coexDataFile) 
	{
		log << "ERROR: Could not open file \"coexCurveVL.dat\"" << endl;
		return 1;
	}
	
	vector<double> kT;
	vector<double> rhoV;
	vector<double> rhoL;
	
	readColumnVectorFromFile(coexDataFile, 0, kT);
	readColumnVectorFromFile(coexDataFile, 1, rhoV);
	readColumnVectorFromFile(coexDataFile, 2, rhoL);
	
	double N = kT.size();
	if (rhoV.size()!=N || rhoL.size()!=N)
	{
		log << "ERROR: mismatch in data length between kT,rhoV,rhoL" << endl;
		return 1;
	}
	
	int M = 5; // do not take the last M temperatures as they are too close to our kT_c estimate
	vector<double> data_x(N-M,0);
	vector<double> data_y(N-M,0);
	
	for (int i=0; i<N-M; i++)
	{
		data_x[i] = std::log(kT_c - kT[N-1-M-i]);
		data_y[i] = std::log(rhoL[N-1-M-i] - rhoV[N-1-M-i]);
	}
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "data_x = ";
	for (int i=0; i<data_x.size(); i++) log << data_x[i] << " ";
	log << endl;
	log << "data_y = ";
	for (int i=0; i<data_y.size(); i++) log << data_y[i] << " ";
	log << endl;
	
	/////////// Linear fit to compute the critical exponent beta ///////////
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << endl;
	
	vector<double> coefficients;
	int status_fit = myPolyFit(data_x, data_y, coefficients, 1);
	
	if (status_fit != 0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: fit not successful" << endl;
		return 1;
	}
	
	double beta = coefficients[1];
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "critical exponent beta = " << beta << endl;
	
	////////////////////////////// data to plot ////////////////////////////
	
	ofstream plotDataFile("dataFit.dat");
	
	plotDataFile << "# x y_data y_fit" << endl;
	
	for (int i=0; i<data_x.size(); i++)
	{
		double x = data_x[i];
		double y_data = data_y[i];
		double y_fit = coefficients[1]*x+coefficients[0];
		
		plotDataFile << x << " " << y_data << " " << y_fit << endl;
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("plot_fit");
	
	outPlotFile << "########## Plot Fit #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'fit.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Fit to extract exponent beta\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"log(kT_c-kT)\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"log(rhoL-rhoV)\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key bottom right" << endl;
	outPlotFile << "set label \"beta = " << beta << "\" font \",20\" at -6.3,-4.4" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'dataFit.dat' using 1:2 with points title 'data', \\" << endl;
	outPlotFile << "     'dataFit.dat' using 1:3 with lines  title 'fit'" << endl;
	
	int sysresult = system("gnuplot plot_fit");
	
	return 0;
}










////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



int coexVL(int argc, char** argv, double &kT_c, double &rho_c)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// geometry
	double dx = 0.0125;
	string pointsFile;
	
	// thermodynamics
	double kTMin = 1.1;
	double kTMax = 1.5;
	double kTStepMax = 0.01;
	double kTStepMin = 1e-6;
	
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
	
	Log log("log_coexVL.dat");
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
			
			// update critical point estimates
			kT_c = kT;
			rho_c = (rhoV+rhoL)/2;
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
	
	ofstream outPlotFile("plot_CoexVL");
	
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
	
	int sysresult = system("gnuplot plot_CoexVL");
	
	return 0;
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// polynomial fit (wrapper for the gsl code)


int myPolyFit(vector<double> data_x, vector<double> data_y,
              vector<double> &coefficients, int orderPolyFit)
{
	cout << "=================================" << endl;
	cout << "== Entering myPolyFit function ==" << endl;
	cout << "=================================" << endl;
	
	// Change units in order to have quantities of order 1
	
	double data_x_scale = data_x[0];
	double data_y_scale = data_y[0];
	
	cout << endl;
	cout << "data_x_scale = " << data_x_scale << endl;
	cout << "data_y_scale = " << data_y_scale << endl;
	
	for (int i=0; i<data_x.size(); i++) data_x[i] /= data_x_scale;
	for (int i=0; i<data_y.size(); i++) data_y[i] /= data_y_scale;
	
	// Convert input to GSL formats
	
	int n = data_x.size();
	int p = orderPolyFit+1;
	
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;

	X = gsl_matrix_alloc (n, p);
	y = gsl_vector_alloc (n);
	w = gsl_vector_alloc (n);
	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc (p, p);
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			gsl_matrix_set (X, i, j, pow(data_x[i],j));
			gsl_vector_set (y, i, data_y[i]);
			gsl_vector_set (w, i, 1.0);
		}
	}
	
	// Solve
	
	double chisq;
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);
	
	// Convert output from GSL format
	// don't forget to return to the original units
	
	coefficients.clear();
	for (int i=0; i<=orderPolyFit; i++)
		coefficients.push_back(gsl_vector_get(c,(i)));
	
	for (int i=0; i<=orderPolyFit; i++)
		coefficients[i] *= data_y_scale / pow(data_x_scale,i);
	
	// Deallocate memory
	
	gsl_matrix_free (X);
	gsl_vector_free (y);
	gsl_vector_free (w);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
	
	return 0;
}



