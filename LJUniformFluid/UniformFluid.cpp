#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <stdio.h>

using namespace std;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>

#include "options.h"

#include "Potential1.h"
#include "Log.h"
#include "myColor.h"
#include "UniformFluid.h"




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume

double uniformFluidOmega(double kT, double mu, double aVdW, 
                         double hsd, double rho)
{
	double ideal = kT*(rho*std::log(rho)-rho);
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCore = kT*rho*(  -std::log(1-etaFMT) + 3.0/2*etaFMT*(2-etaFMT)/pow(1-etaFMT,2)  );
	// Carnahan-Starling
	double hardCore = kT*rho*etaFMT*(4-3*etaFMT)/pow(1-etaFMT,2);
	double meanField = kT*aVdW*rho*rho;
	double muN = mu*rho;
	
	return ideal + hardCore + meanField - muN;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume

double uniformFluidOmegaDerivative(double kT, double mu, double aVdW, 
                                   double hsd, double rho)
{
	double idealDerivative = kT*std::log(rho);
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCoreDerivative = kT*(  -std::log(1-etaFMT) + etaFMT*(7-13.0/2*etaFMT+5.0/2*etaFMT*etaFMT)/pow(1-etaFMT,3)  );
	// Carnahan-Starling
	double hardCoreDerivative = kT*( 3*pow(etaFMT,3) - 9*pow(etaFMT,2) + 8*etaFMT) / pow(1-etaFMT,3);
	double meanFieldDerivative = kT*2*aVdW*rho;
	double muNDerivative = mu;
	
	return idealDerivative + hardCoreDerivative + meanFieldDerivative - muNDerivative;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume

double uniformFluidOmegaDerivative2(double kT, double mu, double aVdW, 
                                    double hsd, double rho)
{
	double idealDerivative2 = kT/rho;
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCoreDerivative2 = kT*4*M_PI/3*pow(hsd/2,3)*(  pow(1-etaFMT,-1)  +  3*pow(1-etaFMT,-4)*etaFMT*(7-13.0/2*etaFMT+5.0/2*etaFMT*etaFMT)  +  pow(1-etaFMT,-3)*(7-13*etaFMT+15.0/2*etaFMT*etaFMT)  );
	// Carnahan-Starling
	double hardCoreDerivative2 = kT*4*M_PI/3*pow(hsd/2,3)*( -2*etaFMT + 8 ) / pow(1-etaFMT,4);
	double meanFieldDerivative2 = kT*2*aVdW;
	double muNDerivative2 = 0;
	
	return idealDerivative2 + hardCoreDerivative2 + meanFieldDerivative2 - muNDerivative2;
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Returns the chemical potential from the density at coexistence

double muFromCoexDensity(double rho, double kT, double aVdW, double hsd)
{
	double ideal = kT * std::log(rho);
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	double hardCore = kT * etaFMT * (3*etaFMT*etaFMT-9*etaFMT+8) * pow(1-etaFMT,-3);
	double meanField = kT * 2*aVdW*rho;
	
	return ideal + hardCore + meanField;
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Returns the associated density given by the coexistence polynomial

double rho2FromRho1(double rho1, double kT, double aVdW, double hsd)
{
	// Compute polynomial parameters for the coexistence equation
	// omega - rho*domega/drho = omega1
	
	double mu1 = muFromCoexDensity(rho1, kT, aVdW, hsd);
	double omega1 = uniformFluidOmega(kT, mu1, aVdW, hsd, rho1);
	double omegaTilde = omega1 / kT * 4*M_PI/3*pow(hsd/2,3);
	double aVdWTilde = aVdW / (4*M_PI/3*pow(hsd/2,3));
	
	// Compute polynomial coefficients
	
	double a5 = aVdWTilde;
	double a4 = 1-3*aVdWTilde;
	double a3 = -1+omegaTilde+3*aVdWTilde;
	double a2 = -1-aVdWTilde-3*omegaTilde;
	double a1 = -1+3*omegaTilde;
	double a0 = -omegaTilde;
	
	// Solve polynomial equation
	
	double a[6] = {a0, a1, a2, a3, a4, a5};
	double z[10];

	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (6);
	gsl_poly_complex_solve (a, 6, w, z);
	gsl_poly_complex_workspace_free (w);

	for (int i = 0; i < 5; i++)
	{
		//printf ("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]);
	}
	
	// Find the biggest real solution
	
	double small_value = 1e-10;
	double eta2 = 0;
	
	for (int i=0; i<5; i++)
	{
		if (abs(z[2*i+1]) < small_value) // is a real root
			if (z[2*i] > eta2)
				eta2 = z[2*i];
	}
	
	return eta2 / (4*M_PI/3 * pow(hsd/2,3));
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Locally defined structures and functions needed by the root finding
// algorithm in "coexistenceDensities"

struct MuDiff_params
{
	double kT, aVdW, hsd;
};

double muDiff_f(double rho1, void *params)
{
	struct MuDiff_params *p = (struct MuDiff_params *) params;
	
	double rho2 = rho2FromRho1(rho1, p->kT, p->aVdW, p->hsd);
	double mu1 = muFromCoexDensity(rho1, p->kT, p->aVdW, p->hsd);
	double mu2 = muFromCoexDensity(rho2, p->kT, p->aVdW, p->hsd);
	
	return mu2-mu1;
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Compute densities at coexistence 

void coexistenceDensities(int argc, char **argv, Log &log,
                          double kT, double aVdW, double hsd,
                          double &rho1, double &rho2, double &muCoex,
                          bool &critical, bool &success)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// rho1 search parameters
	double rho1Min = 0.001;
	double rho1Max = 0.8;
	double rho1FinestStep = 0.001;
	
	/////////// Read options
	
	Options options;
	
	options.addOption("Rho1Min",&rho1Min);
	options.addOption("Rho1Max",&rho1Max);
	options.addOption("Rho1FinestStep",&rho1FinestStep);
	
	options.read(argc, argv);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters (Coex Densities):" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	struct MuDiff_params params = {kT, aVdW, hsd};
	
	//////////////////////////// Search interval ///////////////////////////
	
	// To start the bracketing root finding algorithm, we need an initial
	// interval that contains only one zero, the one we are looking for.
	// See the plot of muDiff      TODO: this algorithm is specific to the
	//                                   shape of the function -> bad idea
	// We use the fact that the zero we are looking for is between the only
	// positive and negative parts of muDiff
	
	double step = rho1Max-rho1Min;
	double a = rho1Min; // current best lower boundary
	double b = rho1Max; // current best upper boundary
	double small_value = 1e-10;
	bool foundLowerBoundary = false;
	bool foundUpperBoundary = false;
	
	while (!foundLowerBoundary || !foundUpperBoundary)
	{
		double c = a;
		
		while (c<b)
		{
			double muDiffAtc = muDiff_f(c, &params);
			
			if (c>a && muDiffAtc>small_value) 
			{
				foundLowerBoundary = true;
				a = c;
				
				log << "found a = " << a << "     muDiff = " << muDiffAtc << endl;
			}
			if (c<b && muDiffAtc<-small_value) 
			{
				foundUpperBoundary = true;
				b = c;
				
				log << "found b = " << b << "     muDiff = " << muDiffAtc << endl;
			}
			if (foundLowerBoundary && foundUpperBoundary) break;
			
			c += step;
		}
		
		step = step/2;
		if (step < rho1FinestStep) break;
	}
	
	if (!foundLowerBoundary || !foundUpperBoundary)
	{
		critical = true;
		
		log << "--------------" << endl;
		log << "The system is above or too close to the critical temperature." << endl;
		log << "mu2-mu1 is zero for all rho1 candidates." << endl;
		log << "--------------" << endl;
		
		return ;
	}
	
	log << "a = " << a << "     muDiff = " << muDiff_f(a,&params) << endl;
	log << "b = " << b << "     muDiff = " << muDiff_f(b,&params) << endl;
	
	//////////////////////// Root finding algorithm ////////////////////////
	
	// local customisations
	int max_iter = 100;
	double precisionRel = 1e-5;
	
	int status;
	int iter = 0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = a;
	double x_lo = a, x_hi = b;
	gsl_function F;

	F.function = &muDiff_f;
	F.params = &params;

	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	printf ("Using %s method\n", gsl_root_fsolver_name (s));
	printf ("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err(est)");

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, precisionRel);

		if (status == GSL_SUCCESS) printf ("Converged:\n");

		printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	
	//////////////////////// Coexistence properties ////////////////////////
	
	rho1 = r;
	rho2 = rho2FromRho1(rho1, kT, aVdW, hsd);
	muCoex = muFromCoexDensity(rho1, kT, aVdW, hsd);
	critical = false;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Locally defined structures and functions needed by findRootdOmegadRho

struct Omega_params
{
	double kT, mu, hsd, aVdW;
};

double uniformFluidOmegaDerivative_f(double rho, void *params)
{
	struct Omega_params *p = (struct Omega_params *) params;
	return uniformFluidOmegaDerivative(p->kT, p->mu, p->aVdW, p->hsd, rho);
}

double uniformFluidOmegaDerivative_df(double rho, void *params)
{
	struct Omega_params *p = (struct Omega_params *) params;
	return uniformFluidOmegaDerivative2(p->kT, p->mu, p->aVdW, p->hsd, rho);
}

void uniformFluidOmegaDerivative_fdf(double rho, void *params, double *y, double *dy)
{
	struct Omega_params *p = (struct Omega_params *) params;
	*y  = uniformFluidOmegaDerivative (p->kT, p->mu, p->aVdW, p->hsd, rho);
	*dy = uniformFluidOmegaDerivative2(p->kT, p->mu, p->aVdW, p->hsd, rho);
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Note: The code for the findRootdOmegadRho function is adapted from the 
//       gsl reference
// Returns the root

double findRootdOmegadRho(double rhoInit, double kT, double mu, double aVdW,
                          double hsd , bool &success)
{
	// customisations
	int max_iter = 100;
	double precisionRel = 1e-5;
	
	int status;
	int iter = 0;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = rhoInit;
	gsl_function_fdf FDF;
	struct Omega_params params = {kT, mu, hsd, aVdW};
	
	FDF.f = &uniformFluidOmegaDerivative_f;
	FDF.df = &uniformFluidOmegaDerivative_df;
	FDF.fdf = &uniformFluidOmegaDerivative_fdf;
	FDF.params = &params;
	
	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, x);
	
	printf ("Using %s method\n",gsl_root_fdfsolver_name (s));
	printf ("%-5s %10s %10s\n","iter", "root", "err (est)");
	
	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (x, x0, 0, precisionRel);
		
		if (status == GSL_SUCCESS) printf ("Converged:\n");
		
		printf ("%5d %10.7f %10.7f\n",iter, x, x - x0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	
	gsl_root_fdfsolver_free (s);
	
	if (status!=GSL_SUCCESS) success = false;
	
	return x;
}

