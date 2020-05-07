
// Program to compute the vapour-solid coexistence curve from Frenkel's
// fits results from https://arxiv.org/pdf/1910.05746.pdf




#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Frenkel_WHDF_rc2.h"

using namespace std;



// For the VS coexistence, the initial guess for the densities needs to
// be more sophisticated.

// It is very hard to find a suitable guess for the densities at low 
// temperatures. We thus restrict ourselves to a temperature kT>0.1 .

double guessV(double kT)
{
	double rho = 1e-4;
	
	if (kT < 0.35)  rho = 1e-7;
	if (kT < 0.25)  rho = 1e-10;
	if (kT < 0.17)  rho = 1e-15;
	if (kT < 0.13)  rho = 1e-20;
	
	return rho;
}



double solidDensityAtSamePressure(double kT, double rhoV)
{
	double rhoS = 1; // initial guess
	double pressureS = pressureRc2Solid(rhoS, kT);
	double pressureV = pressureRc2Vapour(rhoV, kT);
	
	double rho_step = 0.1; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and pressure
		double rhoS_new = rhoS + rho_step;
		double pressureS_new = pressureRc2Solid(rhoS_new, kT);
		
		// keep new density if it decreases the pressure difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(pressureS_new-pressureV) < abs(pressureS-pressureV))
		{
			rhoS = rhoS_new;
			pressureS = pressureS_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	return rhoS;
}




double solidDensityAtSameMu(double kT, double rhoV)
{
	double rhoS = 1; // initial guess
	double muS = muRc2Solid(rhoS, kT);
	double muV = muRc2Vapour(rhoV, kT);
	
	double rho_step = 0.1; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and mu
		double rhoS_new = rhoS + rho_step;
		double muS_new = muRc2Solid(rhoS_new, kT);
		
		// keep new density if it decreases the mu difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(muS_new-muV) < abs(muS-muV))
		{
			rhoS = rhoS_new;
			muS = muS_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	return rhoS;
}




int coexistenceVS(double kT, double &rhoV, double &rhoS)
{
	rhoV = guessV(kT); // initial density
	double rhoS_p  = solidDensityAtSamePressure(kT, rhoV);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoV);
	
	double rho_step = 0.03*guessV(kT); // initial step
	double rho_relprecision = 1e-5;
	
	while ( abs(rho_step/rhoV)           > rho_relprecision )//||
	        //abs((rhoS_p-rhoS_mu)/rhoS_p) > rho_relprecision    )
	{
		// new densities
		double rhoV_new = rhoV + rho_step;
		double rhoS_p_new  = solidDensityAtSamePressure(kT, rhoV_new);
		double rhoS_mu_new = solidDensityAtSameMu(kT, rhoV_new);
		
		// keep new density if it decreases the difference between the 
		// corresponding solid densities, otherwise continue to search in 
		// the other direction with a smaller step
		if (abs(rhoS_p_new-rhoS_mu_new) < abs(rhoS_p-rhoS_mu))
		{
			rhoV = rhoV_new;
			rhoS_p  = rhoS_p_new;
			rhoS_mu = rhoS_mu_new;
		}
		else
		{
			rho_step *= -1;
			rho_step /=  2;
		}
	}
	
	// at the end, return the average of the two solid densities
	rhoS = (rhoS_p + rhoS_mu)/2;
	
	return 0; // assume it went right, no check currently implemented
}




int main()
{
	// Output file
	
	ofstream outDataFile("coexCurveVS_Frenkel.dat");
	outDataFile << "#kT				rhoV			rhoS			"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	// Compute coexistence curve
	
	for (double kT=0.1; kT<=0.5; kT+=0.001)
	{
		double rhoV,rhoS;
		int status = coexistenceVS(kT, rhoV, rhoS);
		
		if (status==0) 
		{
			outDataFile << kT << " 	" << rhoV << " 	" << rhoS << " 	"
			            << endl;
		}
		else
		{
			cout << "ERROR: coexistenceVS failed for kT = " << kT << endl;
		}
	}
	
	// Gnuplot script
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexCurveVS_Frenkel.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Vapour-Solid (Frenkel)\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveVS_Frenkel.dat' using 2:1 title 'vapour', \\" << endl;
	outPlotFile << "     'coexCurveVS_Frenkel.dat' using 3:1 title 'solid'" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}




/*
// For manual checks 

int main()
{
	double kT = 0.8;
	double rhoV = 0.5;
	
	double rhoS_p = solidDensityAtSamePressure(kT, rhoV);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoV);
	
	cout << "kT = " << kT << endl;
	cout << "rhoV = " << rhoV << endl;
	cout << "rhoS_p  = " << rhoS_p << endl;
	cout << "rhoS_mu = " << rhoS_mu << endl;
	
	return 0;
}
*/




