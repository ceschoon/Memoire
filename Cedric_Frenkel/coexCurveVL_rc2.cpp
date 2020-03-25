
// Program to compute the vapour-liquid coexistence curve from Frenkel's
// fits results from https://arxiv.org/pdf/1910.05746.pdf




#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Frenkel_WHDF_rc2.h"

using namespace std;


// For the VL coexistence, the initial guess for the densities need to
// be more sophisticated.

// It is very hard to find a suitable guess for the densities near the
// critical point. I think this is due to a lack of simulation data close
// to the critical point that makes the fit inaccurate in this region.
// We thus restrict ourselves to a temperature kT<0.95.

double guessV(double kT)
{
	double rho = 0.3;
	
	if (kT < 0.99)  rho = 0.15;
	if (kT < 0.95)  rho = 0.08;
	if (kT < 0.85)  rho = 0.03;
	if (kT < 0.75) rho = 0.01;
	if (kT < 0.60) rho = 0.001;
	
	return rho;
}

double guessL(double kT)
{
	double kTc = 1.0;
	double rhoc = 0.3;
	double A = 10;
	
	return rhoc + pow((kTc-kT)/A, 0.25);
}




double liquidDensityAtSamePressure(double kT, double rhoV)
{
	double rhoL = guessL(kT); // initial density
	double pressureL = pressureRc2Liquid(rhoL, kT);
	double pressureV = pressureRc2Vapour(rhoV, kT);
	
	double rho_step = 0.01; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and pressure
		double rhoL_new = rhoL + rho_step;
		double pressureL_new = pressureRc2Liquid(rhoL_new, kT);
		
		// keep new density if it decreases the pressure difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(pressureL_new-pressureV) < abs(pressureL-pressureV))
		{
			rhoL = rhoL_new;
			pressureL = pressureL_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	return rhoL;
}




double liquidDensityAtSameMu(double kT, double rhoV)
{
	double rhoL = guessL(kT); // initial density
	double muL = muRc2Liquid(rhoL, kT);
	double muV = muRc2Vapour(rhoV, kT);
	
	double rho_step = 0.01; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and mu
		double rhoL_new = rhoL + rho_step;
		double muL_new = muRc2Liquid(rhoL_new, kT);
		
		// keep new density if it decreases the mu difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(muL_new-muV) < abs(muL-muV))
		{
			rhoL = rhoL_new;
			muL = muL_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	return rhoL;
}




int coexistenceVL(double kT, double &rhoV, double &rhoL)
{
	rhoV = guessV(kT); // initial density
	double rhoL_p  = liquidDensityAtSamePressure(kT, rhoV);
	double rhoL_mu = liquidDensityAtSameMu(kT, rhoV);
	
	double rho_step = 0.03*guessV(kT); // initial step
	double rho_relprecision = 1e-5;
	
	while ( abs(rho_step/rhoV)           > rho_relprecision )//||
	       // abs((rhoL_p-rhoL_mu)/rhoL_p) > rho_relprecision    )
	{
		// new densities
		double rhoV_new = rhoV + rho_step;
		double rhoL_p_new  = liquidDensityAtSamePressure(kT, rhoV_new);
		double rhoL_mu_new = liquidDensityAtSameMu(kT, rhoV_new);
		
		// keep new density if it decreases the difference between the 
		// corresponding liquid densities, otherwise continue to search in 
		// the other direction with a smaller step
		if (abs(rhoL_p_new-rhoL_mu_new) < abs(rhoL_p-rhoL_mu))
		{
			rhoV = rhoV_new;
			rhoL_p  = rhoL_p_new;
			rhoL_mu = rhoL_mu_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	// at the end, return the average of the two liquid densities
	rhoL = (rhoL_p + rhoL_mu)/2;
	
	// For a manual check
	//rhoV = guessL(kT);
	//rhoL = guessV(kT);
	
	return 0; // assume it went right, no check currently implemented
}




int main()
{
	// Output file
	
	ofstream outDataFile("coexCurveVL_Frenkel.dat");
	outDataFile << "#kT				rhoV			rhoL			"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	// Compute coexistence curve
	
	for (double kT=0.5; kT<=0.95; kT+=0.001)
	{
		double rhoV,rhoL;
		int status = coexistenceVL(kT, rhoV, rhoL);
		
		if (status==0) 
		{
			outDataFile << kT << " 	" << rhoV << " 	" << rhoL << " 	"
			            << endl;
		}
		else
		{
			cout << "ERROR: coexistenceVL failed for kT = " << kT << endl;
		}
	}
	
	// Gnuplot script
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexCurveVL_Frenkel.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Vapour-Liquid (Frenkel)\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveVL_Frenkel.dat' using 2:1 title 'vapour', \\" << endl;
	outPlotFile << "     'coexCurveVL_Frenkel.dat' using 3:1 title 'liquid'" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}




/*
// For manual checks 

int main()
{
	double kT = 0.8;
	double rhoV = 0.5;
	
	double rhoL_p = liquidDensityAtSamePressure(kT, rhoV);
	double rhoL_mu = liquidDensityAtSameMu(kT, rhoV);
	
	cout << "kT = " << kT << endl;
	cout << "rhoV = " << rhoV << endl;
	cout << "rhoL_p  = " << rhoL_p << endl;
	cout << "rhoL_mu = " << rhoL_mu << endl;
	
	return 0;
}
*/




