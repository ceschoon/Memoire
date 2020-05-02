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


int criticalEta(double &eta_c);
int criticalTemperature(int argc, char **argv, double eta_c, 
                        double &kT_c, double &hsd_c, double &aVdW_c);
int criticalDensity(double eta_c, double hsd_c, double &rho_c);
int criticalChemicalPotential(double eta_c, double rho_c, double &mu_c);
int criticalFreeEnergy(double eta_c, double rho_c, double &omega_c);


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to find the critical properties of the uniform phases

int main(int argc, char** argv)
{
	ofstream dataFile("critical_point.dat");
	dataFile << "# Physical Properties at the Critical Point" << endl;
	dataFile << scientific << setprecision(8);
	
	double eta_c, kT_c, hsd_c, aVdW_c, rho_c, mu_c, omega_c;
	
	criticalEta(eta_c);
	dataFile << "eta = " << eta_c << endl;
	
	criticalTemperature(argc, argv, eta_c, kT_c, hsd_c, aVdW_c);
	dataFile << "hsd = " << hsd_c << endl;
	dataFile << "aVdW = " << aVdW_c << endl;
	dataFile << "kT = " << kT_c << endl;
	
	criticalDensity(eta_c, hsd_c, rho_c);
	dataFile << "rho = " << rho_c << endl;
	
	criticalChemicalPotential(eta_c, rho_c, mu_c);
	dataFile << "mu = " << mu_c << endl;
	
	criticalFreeEnergy(eta_c, rho_c, omega_c);
	dataFile << "omega = " << omega_c << endl;
	
	return 0;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


int criticalEta(double &eta_c)
{
	// Solve polynomial equation
	
	printf ("Computing critical eta\n");
	
	double a[6] = {-1, 5, 20, 4, -5, 1};
	double z[10];
	
	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (6);
	gsl_poly_complex_solve (a, 6, w, z);
	gsl_poly_complex_workspace_free (w);
	
	for (int i = 0; i < 5; i++)
	{
		printf ("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]);
	}
	
	// Find the only positive real solution (about 0.13044...)
	
	double small_value = 1e-10;
	
	for (int i=0; i<5; i++)
		if (abs(z[2*i+1]) < small_value) // is a real root
			if (z[2*i] > 0) // is positive
				eta_c = z[2*i];
	
	return 0;
}






int criticalTemperature(int argc, char **argv, double eta_c, 
                        double &kT_c, double &hsd_c, double &aVdW_c)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT_init = 1;
	double rel_precision = 1e-8;
	
	/////////// Read options
	
	Options options;
	
	options.addOption("kT_init", &kT_init);
	options.addOption("rel_precision", &rel_precision);
	
	options.read(argc, argv);
	
	Log log("log_kTc.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	//////////////////////// Compute aVdW for kT=1 /////////////////////////
	
	double hsd0, aVdW0;
	compute_hsd_aVdW(argc, argv, log, 1.0, hsd0, aVdW0, false);
	
	cout << "Computed aVdW0 = " << aVdW0 << endl;
	log << endl;
	
	/////////////////////// Compute kT_c Iteratively ///////////////////////
	
	// We need to find the temperature at which
	//   kT*hsd^3 = -12/PI * aVdW0 * eta*(1-eta)^4 / (eta^4-4*eta^3+4*eta^2+4*eta+1)
	// This is equivalent to finding the zero of the following expression
	//   kT*hsd^3 + 12/PI * aVdW0 * eta*(1-eta)^4 / (eta^4-4*eta^3+4*eta^2+4*eta+1)
	
	// Initialisation
	
	double kT = kT_init;
	double kT_step = 0.05;
	
	double hsd, aVdW_dummy;
	compute_hsd_aVdW(argc, argv, log, kT, hsd, aVdW_dummy, true);
	
	double hsd3 = hsd*hsd*hsd;
	double eta2 = eta_c*eta_c;
	double eta3 = eta2*eta_c;
	double eta4 = eta3*eta_c;
	
	double expr = kT*hsd3 + 12/PI*aVdW0 * eta_c * pow(1-eta_c,4) 
	              / (eta4 -4*eta3 +4*eta2 +4*eta_c +1);
	
	log << "Initialising with kT = " << kT << endl;
	log << "                 hsd = " << hsd << endl;
	log << "                expr = " << expr << endl;
	
	// Loop
	
	while (kT_step > kT*rel_precision)
	{
		double kT_new = kT + kT_step;
		double hsd_new;
		compute_hsd_aVdW(argc, argv, log, kT, hsd_new, aVdW_dummy, true);
		
		log << "Testing kT_new = " << kT_new << endl;
		log << "       hsd_new = " << hsd_new << endl;
		
		double hsd3_new = hsd_new*hsd_new*hsd_new;
		double expr_new = kT_new*hsd3_new + 12/PI*aVdW0 * eta_c * pow(1-eta_c,4)
		                  / (eta4 -4*eta3 +4*eta2 +4*eta_c +1);
		
		log << "      expr_new = " << expr_new << endl;
		
		if (expr*expr_new<0) // if changed sign
		{
			log << "Dividing kT_step by two" << endl;
			kT_step /= 2;
		}
		else
		{
			kT = kT_new;
			hsd = hsd_new;
			expr = expr_new;
		}
	}
	
	// Save final kT as the critical temperature
	
	kT_c = kT;
	hsd_c = hsd;
	aVdW_c = aVdW0/kT_c;
	
	return 0;
}






int criticalDensity(double eta_c, double hsd_c, double &rho_c)
{
	double hsd3 = hsd_c*hsd_c*hsd_c;
	rho_c = eta_c / ( 4*PI/3*hsd3/8 );
	
	return 0;
}






int criticalChemicalPotential(double eta_c, double rho_c, double &mu_c)
{
	double mu_ideal = std::log(rho_c);
	
	double eta4 = pow(eta_c,4);
	double eta3 = pow(eta_c,3);
	double eta2 = pow(eta_c,2);
	
	double mu_excess = - (eta4 - 16*eta3 + 21*eta2 - 4*eta_c +1) 
	                     / pow(1-eta_c,4);
	
	mu_c = mu_ideal + mu_excess;
	
	return 0;
}






int criticalFreeEnergy(double eta_c, double rho_c, double &omega_c)
{
	double eta3 = pow(eta_c,3);
	double eta2 = pow(eta_c,2);
	
	omega_c = rho_c * (7*eta3 + 7*eta2 + 3*eta_c -1) /2 /pow(1-eta_c,3);
	
	return 0;
}




