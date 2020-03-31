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



int isotherm(int argc, char** argv, double kT);

int myPolyFit(vector<double> data_x, vector<double> data_y,
              vector<double> &coefficients, int orderPolyFit);



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the critical exponent delta defined by
// omega_c-omega ~ (rho-rho_c)^delta

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// estimation of the critical point properties
	double kT_c = 1.13781055 + 1e-6;
	double rho_c = (0.2383+0.2404)/2;
	double mu_c = -2.8311;
	double omega_c = -0.085925;
	
	Options options;
	options.addOption("kT_c", &kT_c);
	options.addOption("rho_c", &rho_c);
	options.addOption("mu_c", &mu_c);
	options.addOption("omega_c", &omega_c);
	options.read(argc, argv);
	
	Log log("log_exponent.dat");
	options.write(log);
	
	//////////////////////// Compute Isotherm Curve ////////////////////////
	
	int status_isotherm = isotherm(argc, argv, kT_c);
	
	if (status_isotherm != 0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Isotherm computation not successful" << endl;
		return 1;
	}
	
	/////////////////// Get and tranform data for the fit //////////////////
	
	ifstream isothermDataFile("isotherm.dat");
	if (!isothermDataFile) 
	{
		log << "ERROR: Could not open file \"isotherm.dat\"" << endl;
		return 1;
	}
	
	vector<double> rho;
	vector<double> omega;
	
	readColumnVectorFromFile(isothermDataFile, 0, rho);
	readColumnVectorFromFile(isothermDataFile, 1, omega);
	
	double N = rho.size();
	if (omega.size()!=N)
	{
		log << "ERROR: mismatch in data length between rho,omega" << endl;
		return 1;
	}
	
	int M = 0; // do not take the last M densities as they are too close to our rho_c estimate
	vector<double> data_x(N-M,0);
	vector<double> data_y(N-M,0);
	
	for (int i=0; i<N-M; i++)
	{
		data_x[i] = std::log(rho[i] - rho_c);
		data_y[i] = std::log(omega_c - omega[i]);
	}
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "data_x = ";
	for (int i=0; i<data_x.size(); i++) log << data_x[i] << " ";
	log << endl;
	log << "data_y = ";
	for (int i=0; i<data_y.size(); i++) log << data_y[i] << " ";
	log << endl;
	
	/////////// Linear fit to compute the critical exponent delta //////////
	
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
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "fit coefficient 0 = " << coefficients[0] << endl;
	log << "fit coefficient 1 = " << coefficients[1] << endl;
	
	double delta = coefficients[1];
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "critical exponent delta = " << delta << endl;
	
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
	outPlotFile << "set title \"Fit to extract exponent delta\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"log(rho-rho_c)\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"log(omega_c-omega)\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key bottom right" << endl;
	outPlotFile << "#set label \"delta = " << delta << "\" font \",20\" at -3,-9" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'dataFit.dat' using 1:2 with points title 'data', \\" << endl;
	outPlotFile << "     'dataFit.dat' using 1:3 with lines  title 'fit (delta = "<<delta<<")'" << endl;
	
	int sysresult = system("gnuplot plot_fit");
	
	return 0;
}










////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



int isotherm(int argc, char** argv, double kT)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// geometry
	double dx = 0.0125;
	string pointsFile;
	
	// thermodynamics
	double mu_min = -2.8311+0.01;
	double mu_max = -2.8311+0.1;
	double mu_step = 0.001;
	
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
	options.addOption("mu_min", &mu_min);
	options.addOption("mu_max", &mu_max);
	options.addOption("mu_step", &mu_step);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	options.read(argc, argv);
	
	Log log("log_isotherm.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	////////////////// Coexistence properties for each kT //////////////////
	
	ofstream outDataFile("isotherm.dat");
	outDataFile << scientific << setprecision(8);
	outDataFile << "# Isotherm curve (fluid) for kT = " << kT << endl;
	outDataFile << "# rho          	omega          	"
	            << endl;
	
	for (double mu=mu_min; mu<=mu_max; mu+=mu_step)
	{
		// Compute the free energy
		
		double omega, rho;
		
		Log log_freeEnergyCalc("log_freeEnergyCalc.dat");
		
		int statusFreeEnergyCalc = fixedkTMuFluid(kT, mu, argc, argv,
			log_freeEnergyCalc, omega, rho);
		
		// Save densities and free energies
		
		if (statusFreeEnergyCalc==0)
		{
			outDataFile << rho << " 	" << omega << endl;
			log << "For mu = " << mu << " we have rho = " << rho 
			    << " and omega = " << omega << endl;
		}
		else
		{
			log << "Computation (fixedkTMuFluid) failed for mu = " << mu << endl;
		}
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("plot_isotherm");
	
	outPlotFile << "########## Plot Isotherm Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'isotherm.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Isotherm for kT=" << kT << "\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"pressure\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key bottom left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'isotherm.dat' using 1:(-($2)) notitle" << endl;
	
	int sysresult = system("gnuplot plot_isotherm");
	
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



