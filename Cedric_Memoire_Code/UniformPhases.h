////////////////////////////////////////////////////////////////////////////
//                                                                        //
//               DFT Computations for the Uniform Phases                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


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
#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "SolidDensity.h"

#include "utilities.h"

// Note on units: 
//   the chemical potential mu is expressed in kT (it is beta*mu)
//   the free energies (given per unit volume) are expressed in kT too

// free energy and derivatives of the uniform phases

double uniformOmega(double kT, double mu, double aVdW, 
                    double hsd, double rho);

double uniformOmegaDerivative(double kT, double mu, double aVdW, 
                              double hsd, double rho);

double uniformOmegaDerivative2(double kT, double mu, double aVdW, 
                               double hsd, double rho);

// functions used to find the coexistence properties between the uniform phases

double muFromCoexDensity(double rho, double kT, double aVdW, double hsd);

int vapourLiquidCoexistence (int argc, char **argv, Log &log,
                             double kT, double &rhoV, double &rhoL, 
                             double &muCoex, double &freeEnergyCoex,
                             bool &superCritical);

// find roots of dOmegadRho to extract densities

void spinodalDensities(int argc, char **argv, Log &log, double kT,
                       double &rho1, double &rho2, bool &spinodalExists);

void spinodalDensities(double aVdW, double hsd, double &rho1, double &rho2,
                       bool &spinodalExists);

void findRootsdOmegadRhoSpinodal(double kT, double mu, double aVdW, 
                                 double hsd, double rhoMin, double rhoMax,
                                 vector<double> &roots);

// fluid properties for given (kT,mu)

int fixedkTMuFluid(double kT, double mu, int argc, char** argv, Log &log, 
                   double &freeEnergyVapour, double &densityVapour,
                   double &freeEnergyLiquid, double &densityLiquid,
                   bool &supercritical);

int compute_hsd_aVdW(int argc, char **argv, Log &log,
                     double kT, double &hsd, double &aVdW,
                     bool hsdOnly = false);




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume (-pressure)

double uniformOmega(double kT, double mu, double aVdW, 
                    double hsd, double rho)
{
	double ideal = rho*std::log(rho)-rho;
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCore = rho*(  -std::log(1-etaFMT) + 3.0/2*etaFMT*(2-etaFMT)/pow(1-etaFMT,2)  );
	// Carnahan-Starling
	double hardCore = rho*etaFMT*(4-3*etaFMT)/pow(1-etaFMT,2);
	double meanField = aVdW*rho*rho;
	double muN = mu*rho;
	
	return ideal + hardCore + meanField - muN;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume

double uniformOmegaDerivative(double kT, double mu, double aVdW, 
                              double hsd, double rho)
{
	double idealDerivative = std::log(rho);
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCoreDerivative = (  -std::log(1-etaFMT) + etaFMT*(7-13.0/2*etaFMT+5.0/2*etaFMT*etaFMT)/pow(1-etaFMT,3)  );
	// Carnahan-Starling
	double hardCoreDerivative = ( 3*pow(etaFMT,3) - 9*pow(etaFMT,2) + 8*etaFMT) / pow(1-etaFMT,3);
	double meanFieldDerivative = 2*aVdW*rho;
	double muNDerivative = mu;
	
	return idealDerivative + hardCoreDerivative + meanFieldDerivative - muNDerivative;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remark: free energy per unit volume

double uniformOmegaDerivative2(double kT, double mu, double aVdW, 
                               double hsd, double rho)
{
	double idealDerivative2 = 1/rho;
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	// Percus Yevick
	//double hardCoreDerivative2 = 4*M_PI/3*pow(hsd/2,3)*(  pow(1-etaFMT,-1)  +  3*pow(1-etaFMT,-4)*etaFMT*(7-13.0/2*etaFMT+5.0/2*etaFMT*etaFMT)  +  pow(1-etaFMT,-3)*(7-13*etaFMT+15.0/2*etaFMT*etaFMT)  );
	// Carnahan-Starling
	double hardCoreDerivative2 = 4*M_PI/3*pow(hsd/2,3)*( -2*etaFMT + 8 ) / pow(1-etaFMT,4);
	double meanFieldDerivative2 = 2*aVdW;
	double muNDerivative2 = 0;
	
	return idealDerivative2 + hardCoreDerivative2 + meanFieldDerivative2 - muNDerivative2;
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Returns the chemical potential from the density at coexistence

double muFromCoexDensity(double rho, double kT, double aVdW, double hsd)
{
	double ideal = std::log(rho);
	double etaFMT = 4*M_PI/3*pow(hsd/2,3)*rho;
	double hardCore = etaFMT * (3*etaFMT*etaFMT-9*etaFMT+8) * pow(1-etaFMT,-3);
	double meanField = 2*aVdW*rho;
	
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
	double omega1 = uniformOmega(kT, mu1, aVdW, hsd, rho1);
	double omegaTilde = omega1 * 4*M_PI/3*pow(hsd/2,3);
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


int vapourLiquidCoexistence (int argc, char **argv, Log &log,
                             double kT, double &rhoV, double &rhoL, 
                             double &muCoex, double &freeEnergyCoex,
                             bool &superCritical)
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
	
	//////////////////// Compute hsd and VdW parameters ////////////////////
	
	double hsd,aVdW;
	int statusHSD = compute_hsd_aVdW(argc, argv, log, kT, hsd, aVdW);
	if (statusHSD != 0) 
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Failed to compute hsd and aVdW" << endl;
		
		return 1;
	}
	
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
		superCritical = true;
		
		log << "--------------" << endl;
		log << "The system is above or too close to the critical temperature." << endl;
		log << "mu2-mu1 is zero for all rho1 candidates." << endl;
		log << "--------------" << endl;
		
		return 1;
	}
	
	log << "a = " << a << "     muDiff = " << muDiff_f(a,&params) << endl;
	log << "b = " << b << "     muDiff = " << muDiff_f(b,&params) << endl;
	
	//////////////////////// Root finding algorithm ////////////////////////
	
	// local customisations
	int max_iter = 100;
	double precisionRel = 1e-8;
	
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
	
	if (iter >= max_iter) 
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Exceeded max number of iterations allowed" << endl;
		
		return 1;
	}

	gsl_root_fsolver_free (s);
	
	//////////////////////// Coexistence properties ////////////////////////
	
	rhoV = r;
	rhoL = rho2FromRho1(r, kT, aVdW, hsd);
	
	if (rhoV>rhoL) 
	{
		double temp = rhoV;
		rhoV = rhoL;
		rhoL = temp;
	}
	
	muCoex = muFromCoexDensity(rhoV, kT, aVdW, hsd);
	freeEnergyCoex = uniformOmega(kT, muCoex, aVdW, hsd, rhoV);
	superCritical = false;
	
	return 0;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Locally defined structures and functions needed by findRootdOmegadRho

struct Omega_params
{
	double kT, mu, hsd, aVdW;
};

double uniformOmegaDerivative_f(double rho, void *params)
{
	struct Omega_params *p = (struct Omega_params *) params;
	return uniformOmegaDerivative(p->kT, p->mu, p->aVdW, p->hsd, rho);
}

double uniformOmegaDerivative_df(double rho, void *params)
{
	struct Omega_params *p = (struct Omega_params *) params;
	return uniformOmegaDerivative2(p->kT, p->mu, p->aVdW, p->hsd, rho);
}

void uniformOmegaDerivative_fdf(double rho, void *params, double *y, double *dy)
{
	struct Omega_params *p = (struct Omega_params *) params;
	*y  = uniformOmegaDerivative (p->kT, p->mu, p->aVdW, p->hsd, rho);
	*dy = uniformOmegaDerivative2(p->kT, p->mu, p->aVdW, p->hsd, rho);
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Find the spinodal densities by solving the polynomial equation
// 0 = d2Omega/drho2
// the solutions are returned such as rho1 < rho2

void spinodalDensities(int argc, char **argv, Log &log, double kT,
                       double &rho1, double &rho2, bool &spinodalExists)
{
	double hsd, aVdW;
	compute_hsd_aVdW(argc, argv, log, kT, hsd, aVdW);
	
	spinodalDensities(aVdW, hsd, rho1, rho2, spinodalExists);
}

void spinodalDensities(double aVdW, double hsd, double &rho1, double &rho2,
                       bool &spinodalExists)
{
	// Compute polynomial parameters
	
	double aVdWTilde = aVdW / (4*M_PI/3*pow(hsd/2,3));
	
	// Compute polynomial coefficients
	
	double a5 = 2*aVdWTilde;
	double a4 = 1-8*aVdWTilde;
	double a3 = -4+12*aVdWTilde;
	double a2 = 4-8*aVdWTilde;
	double a1 = 4+2*aVdWTilde;
	double a0 = 1;
	
	// Solve polynomial equation
	
	double a[6] = {a0, a1, a2, a3, a4, a5};
	double z[10];

	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (6);
	int status = gsl_poly_complex_solve (a, 6, w, z);
	if (status) //GSL_EFAILED
	{
		cout << "ERROR: In spinodalDensities: Root solving qr method failed to converge" << endl;
		return ;
	}
	gsl_poly_complex_workspace_free (w);

	cout << "Spinodal roots are" << endl;
	for (int i = 0; i < 5; i++)
	{
		printf ("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]);
	}
	
	// Find the valid real solutions (between 0 and 1)
	
	double small_value = 1e-10;
	vector<double> realRoots;
	
	for (int i=0; i<5; i++)
	{
		if (abs(z[2*i+1]) < small_value) // is a real root
			if (small_value < z[2*i] && z[2*i] < 1) // is a valid packing fraction
				realRoots.push_back(z[2*i]);
	}
	
	// There should be 0 or two valid roots
	
	if (realRoots.size()!=0 && realRoots.size()!=2)
	{
		cout << "ERROR: In spinodalDensities: Incorrect number of valid "
		     << "real solutions for the spinodal equation" << endl;
		return;
	}
	
	if (realRoots.size()==0)
		spinodalExists = false;
	
	if (realRoots.size()==2)
	{
		spinodalExists = true;
		
		if (realRoots[0]<realRoots[1])
		{
			rho1 = realRoots[0]/(4*M_PI/3*pow(hsd/2,3));
			rho2 = realRoots[1]/(4*M_PI/3*pow(hsd/2,3));
		}
		else
		{
			rho1 = realRoots[1]/(4*M_PI/3*pow(hsd/2,3));
			rho2 = realRoots[0]/(4*M_PI/3*pow(hsd/2,3));
		}
	}
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Locally defined structures and functions needed by the root finding
// algorithm in "findRootsdOmegadRhoSpinodal"

struct dOmegadrho_params
{
	double kT, mu, aVdW, hsd;
};

double dOmegadrho_f(double rho, void *params)
{
	struct dOmegadrho_params *p = (struct dOmegadrho_params *) params;
	
	return uniformOmegaDerivative(p->kT, p->mu, p->aVdW, p->hsd, rho);
}

void zeroBisection(double xMin, double xMax, double &root, 
                   dOmegadrho_params params)
{
	// local customisations
	int max_iter = 100;
	double precisionRel = 1e-8;
	
	////////////// 
	int status;
	int iter = 0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = xMin;
	double x_lo = xMin, x_hi = xMax;
	gsl_function F;

	F.function = &dOmegadrho_f;
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
	//////////////
	
	// save the root
	root = r;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Find the roots of dOmega/drho 

// We use the spinodal solutions to limit the search interval for the
// bisection algorithm

// Returns the physical roots only

void findRootsdOmegadRhoSpinodal(double kT, double mu, double aVdW, 
                                 double hsd, double rhoMin, double rhoMax,
                                 vector<double> &roots)
{
	roots.clear();
	struct dOmegadrho_params params = {kT, mu, aVdW, hsd};
	
	// Find the spinodal densities
	
	double rho1 = 0;
	double rho2 = 0;
	bool spinodalExists = false;
	
	spinodalDensities(aVdW, hsd, rho1, rho2, spinodalExists);
	
	cout << "spinodalExists = " << spinodalExists << endl;
	cout << "rho1 = " << rho1 << endl;
	cout << "rho2 = " << rho2 << endl;
	
	// Check if a spinodal exists
	// If not, the is only one solution that we search between rhoMin and 
	// rhoMax
	
	if (!spinodalExists) 
	{
		double root;
		zeroBisection(rhoMin, rhoMax, root, params);
		roots.push_back(root);
		return;
	}
	
	// If there is a spinodal, we need to know how many roots there are
	// and in which interval to look for them
	// We do this by evaluating dOmega/drho at the spinodal densities
	
	double dOmegadrho_at_rho1 = uniformOmegaDerivative(kT, mu, aVdW, 
		hsd, rho1);
	double dOmegadrho_at_rho2 = uniformOmegaDerivative(kT, mu, aVdW, 
		hsd, rho2);
	
	cout << "dOmegadrho_at_rho1 = " << dOmegadrho_at_rho1 << endl;
	cout << "dOmegadrho_at_rho2 = " << dOmegadrho_at_rho2 << endl;
	
	// In the case where there is only the vapor solution,
	// we need to search the root between rhoMin and the first spinodal
	
	if (dOmegadrho_at_rho1>0 && dOmegadrho_at_rho2>0)
	{
		double root;
		zeroBisection(rhoMin, rho1, root, params);
		roots.push_back(root);
		return;
	}
	
	// In the case where there is only the liquid solution,
	// we need to search the root between the second spinodal and rhoMax
	
	if (dOmegadrho_at_rho1<0 && dOmegadrho_at_rho2<0)
	{
		double root;
		zeroBisection(rho2, rhoMax, root, params);
		roots.push_back(root);
		return;
	}
	
	// In the case where there are two solutions,
	// we need to search the first root between rhoMin and the first spinodal,
	// and the third between the second spinodal and rhoMax.
	// The second root between the two spinodals has no physical meaning
	
	if (dOmegadrho_at_rho1>0 && dOmegadrho_at_rho2<0)
	{
		double root;
		zeroBisection(rhoMin, rho1, root, params); // vapor solution
		roots.push_back(root);
		zeroBisection(rho2, rhoMax, root, params); // liquid solution
		roots.push_back(root);
		return;
	}
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Function to compute properties of the fluid at fixed (kT,mu)


int fixedkTMuFluid(double kT, double mu, int argc, char** argv, Log &log, 
                   double &freeEnergyVapour, double &densityVapour,
                   double &freeEnergyLiquid, double &densityLiquid,
                   bool &supercritical)
{
	/////////////////////////////// Options ////////////////////////////////
	
	// data directory
	string dataDir = "data";
	
	// control for the limits
	double rhoFluidMax = 1.5;
	double rhoFluidMin = 1e-100;
	
	Options options;
	
	options.addOption("dataDir", &dataDir);
	options.addOption("RhoFluidMax", &rhoFluidMax);
	options.addOption("RhoFluidMin", &rhoFluidMin);
	
	options.read(argc, argv);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysres = system(("mkdir -p "+dataDir).c_str());
	
	//////////// Check if computation has not already been done ////////////
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	
	ifstream dataInFile(dataDir+"/kTMuFluid_"+"kT="+sskT.str()+"_mu="
	                    +ssMu.str()+".dat");
	
	// if file exists, get the result from there
	if (dataInFile)
	{
		// check if it was a successful computation
		int success;
		readDataFromFile(dataInFile, "success", success);
		
		if (success==1)  // bool "true" 
		{
			log << endl;
			log << "Returning existing computation" << endl;
			
			// read results in file
			int supercritical_temp;
			readDataFromFile(dataInFile, "freeEnergyVapour", freeEnergyVapour);
			readDataFromFile(dataInFile, "freeEnergyLiquid", freeEnergyLiquid);
			readDataFromFile(dataInFile, "densityVapour", densityVapour);
			readDataFromFile(dataInFile, "densityLiquid", densityLiquid);
			readDataFromFile(dataInFile, "supercritical", supercritical_temp);
			supercritical = (supercritical_temp==1); // convert to bool
			
			log << endl;
			log << "freeEnergyVapour = " << freeEnergyVapour << endl;
			log << "freeEnergyLiquid = " << freeEnergyLiquid << endl;
			log << "densityVapour = " << densityVapour << endl;
			log << "densityLiquid = " << densityLiquid << endl;
			log << "supercritical = " << supercritical << endl;
			
			return 0;
		}
		
		if (success==0)  // bool "false" 
		{
			log << endl;
			log << "Returning existing computation (is a failed one)" << endl;
			
			return 1;
		}
	}
	
	////////////// By default, store the result as a failure ///////////////
	// it will be overwritten in the case of a success
	
	ofstream dataFile(dataDir+"/kTMuFluid_"+"kT="+sskT.str()+"_mu="
	                  +ssMu.str()+".dat");
	
	dataFile << "# Result for fluid computation at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << endl;
	dataFile << "success = " << false << endl;
	
	
	////////////////// Compute aVdW and hsd parameters /////////////////////
	
	double hsd,aVdW;
	
	int statusHSD = compute_hsd_aVdW(argc, argv, log, kT, hsd, aVdW);
	if (statusHSD != 0) 
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Failed to compute hsd and aVdW" << endl;
		
		return 1;
	}
	
	// Report
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "kT = " << kT << endl;
	log << "hsd = " << hsd << endl;
	log << "aVdW = " << aVdW << endl;
	
	
	///////////////// Compute free energy for the fluid ////////////////////
	
	bool successFluid = true;
	
	try
	{
		vector<double> roots;
		findRootsdOmegadRhoSpinodal(kT, mu, aVdW, hsd, 
			rhoFluidMin, rhoFluidMax, roots);
		
		if (roots.size()==1)
		{
			densityVapour = roots[0];
			freeEnergyVapour = uniformOmega(kT, mu, aVdW, hsd, roots[0]);
			densityLiquid = densityVapour;
			freeEnergyLiquid = freeEnergyVapour;
			supercritical = true;
		}
		else if (roots.size()==2)
		{
			densityVapour = (roots[0]<roots[1])?roots[0]:roots[1];
			densityLiquid = (roots[0]<roots[1])?roots[1]:roots[0];
			
			freeEnergyVapour  = uniformOmega(
				kT, mu, aVdW, hsd, densityVapour);
			freeEnergyLiquid = uniformOmega(
				kT, mu, aVdW, hsd, densityLiquid);
			
			supercritical = false;
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Anormal number of roots for the fluid dOmega/drho" << endl;
			successFluid = false;
		}
	}
	catch (...)
	{
		successFluid = false;
	}
	
	// Check status
	
	if (successFluid)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Vapour Density = " << densityVapour << endl;
		log << "Vapour Free Energy = " << freeEnergyVapour << endl;
		log << "Liquid Density = " << densityLiquid << endl;
		log << "Liquid Free Energy = " << freeEnergyLiquid << endl;
		log << "supercritical = " << supercritical << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Fluid computation FAILED" << endl;
		
		return 1;
	}
	
	//////////////////////////// Store the result //////////////////////////
	
	dataFile.open(dataDir+"/kTMuFluid_"+"kT="+sskT.str()+"_mu="
	              +ssMu.str()+".dat");
	
	dataFile << "# Result for fluid computation at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << endl;
	dataFile << "freeEnergyVapour = " << scientific << setprecision(8) << freeEnergyVapour << endl;
	dataFile << "freeEnergyLiquid = " << scientific << setprecision(8) << freeEnergyLiquid << endl;
	dataFile << "densityVapour    = " << scientific << setprecision(8) << densityVapour << endl;
	dataFile << "densityLiquid    = " << scientific << setprecision(8) << densityLiquid << endl;
	dataFile << "supercritical = " << supercritical << endl;
	dataFile << endl;
	dataFile << "success = " << true << endl;
	
	return 0;
}








////////////////////////////////////////////////////////////////////////////



int compute_hsd_aVdW(int argc, char **argv, Log &log,
                     double kT, double &hsd, double &aVdW,
                     bool hsdOnly)
{
	///////////////////////////////// Log //////////////////////////////////
	
	log << myColor::GREEN << "==========================================" << myColor::RESET << endl << "#" << endl;
	log << "== Computing hsd and aVdW for kT=" << kT << " ==" << endl;
	
	/////////////////////////////// Options ////////////////////////////////
	
	// geometry
	double dx = -1; // -1 means do not use the grid aVdW
	string pointsFile;
	string useAnalyticWeights = "false";
	
	// potential
	string potentialType = "LJ";
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	
	Options options;
	
	// geometry
	options.addOption("dx", &dx);
	options.addOption("IntegrationPointsFile", &pointsFile);
	options.addOption("UseAnalyticWeights", &useAnalyticWeights);
	
	// potential
	options.addOption("PotentialType", &potentialType);
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	options.read(argc, argv);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	///////////////////////// Initialising Potentials //////////////////////
	
	LJ potentialLJ(sigma1, eps1, rcut1);
	WHDF potentialWHDF(sigma1, eps1, rcut1);
	
	if (potentialType == "LJ") log << "Using the LJ potential" << endl;
	else if (potentialType == "WHDF") log << "Using the WHDF potential" << endl;
	else
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Potential type not identified" << endl;
		log << "       Unable to compute hsd and aVdW for fluid computations" << endl;
		
		return 1;
	}
	
	////////////////////////////////////////////////////////////////////////
	
	if (potentialType == "LJ")   hsd = potentialLJ.getHSD(kT);
	if (potentialType == "WHDF") hsd = potentialWHDF.getHSD(kT);
	
	if (!hsdOnly && potentialType == "LJ")   aVdW = potentialLJ.getVDW_Parameter(kT);;
	if (!hsdOnly && potentialType == "WHDF") aVdW = potentialWHDF.getVDW_Parameter(kT);;
	
	// if we want to use the VdW parameter of a grid of spacing dx
	if (!hsdOnly && dx>0) 
	{
		// values doesn't matter here
		int Ngrid = 127;
		double L[3] = {Ngrid*dx, Ngrid*dx, Ngrid*dx};
		double mu = -6;
		
		SolidDensity theDensity1(dx, L, hsd);
		
		if (useAnalyticWeights == "true")
		{
			log << "Using analytic evaluation of weights" << endl;
			FMT_Species_Analytic species1(theDensity1,hsd, mu, Ngrid);
			
			Interaction_Interpolation_QF interacLJ(species1,species1,potentialLJ,kT,log);
			Interaction_Interpolation_QF interacWHDF(species1,species1,potentialWHDF,kT,log);
			
			if (potentialType == "LJ")
			{
				interacLJ.initialize();
				aVdW = interacLJ.getVDWParameter()/2; // divide by 2 as not the same definition
			}
			if (potentialType == "WHDF")
			{
				interacWHDF.initialize();
				aVdW = interacWHDF.getVDWParameter()/2; // divide by 2 as not the same definition
			}
		}
		else if (useAnalyticWeights == "false")
		{
			log << "Using numeric evaluation of weights" << endl;
			FMT_Species_Numeric species1(theDensity1,hsd,pointsFile, mu, Ngrid);
			
			Interaction interacLJ(species1,species1,potentialLJ,kT,log, pointsFile);
			Interaction interacWHDF(species1,species1,potentialWHDF,kT,log, pointsFile);
			
			if (potentialType == "LJ")
			{
				interacLJ.initialize();
				aVdW = interacLJ.getVDWParameter()/2; // divide by 2 as not the same definition
			}
			if (potentialType == "WHDF")
			{
				interacWHDF.initialize();
				aVdW = interacWHDF.getVDWParameter()/2; // divide by 2 as not the same definition
			}
		}
		else
		{
			log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Unrecognised option for UseAnalyticWeights" << endl;
			log << "       Unable to compute hsd and aVdW for fluid computations" << endl;
			
			return 1;
		}
	}
	
	return 0;
}



////////////////////////////////////////////////////////////////////////////