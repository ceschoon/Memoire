
// Program to compute the liquid-solid coexistence curve from Frenkel's
// fits results from https://arxiv.org/pdf/1910.05746.pdf




#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Frenkel_WHDF_rc2.h"

using namespace std;




double solidDensityAtSamePressure(double kT, double rhoL)
{
	double rhoS = 1; // initial guess
	double pressureS = pressureRc2Solid(rhoS, kT);
	double pressureL = pressureRc2Liquid(rhoL, kT);
	
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
		if (abs(pressureS_new-pressureL) < abs(pressureS-pressureL))
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




double solidDensityAtSameMu(double kT, double rhoL)
{
	double rhoS = 1; // initial guess
	double muS = muRc2Solid(rhoS, kT);
	double muL = muRc2Liquid(rhoL, kT);
	
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
		if (abs(muS_new-muL) < abs(muS-muL))
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




int coexistenceLS(double kT, double &rhoL, double &rhoS)
{
	rhoL = 0.8; // initial density
	double rhoS_p  = solidDensityAtSamePressure(kT, rhoL);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoL);
	
	double rho_step = 0.1; // initial step
	double rho_relprecision = 1e-5;
	
	while ( abs(rho_step/rhoL)           > rho_relprecision )//||
	       // abs((rhoS_p-rhoS_mu)/rhoS_p) > rho_relprecision    )
	{
		// new densities
		double rhoL_new = rhoL + rho_step;
		double rhoS_p_new  = solidDensityAtSamePressure(kT, rhoL_new);
		double rhoS_mu_new = solidDensityAtSameMu(kT, rhoL_new);
		
		// keep new density if it decreases the difference between the 
		// corresponding solid densities, otherwise continue to search in 
		// the other direction with a smaller step
		if (abs(rhoS_p_new-rhoS_mu_new) < abs(rhoS_p-rhoS_mu))
		{
			rhoL = rhoL_new;
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
	
	ofstream outDataFile("coexCurveLS_Frenkel.dat");
	outDataFile << "#kT				rhoL			rhoS			"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	// Compute coexistence curve
	
	for (double kT=0.5; kT<=1.4; kT+=0.01)
	{
		double rhoL,rhoS;
		int status = coexistenceLS(kT, rhoL, rhoS);
		
		if (status==0) 
		{
			outDataFile << kT << " 	" << rhoL << " 	" << rhoS << " 	"
			            << endl;
		}
		else
		{
			cout << "ERROR: coexistenceLS failed for kT = " << kT << endl;
		}
	}
	
	// Gnuplot script
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexCurveLS_Frenkel.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Liquid-Solid (Frenkel)\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveLS_Frenkel.dat' using 2:1 title 'liquid', \\" << endl;
	outPlotFile << "     'coexCurveLS_Frenkel.dat' using 3:1 title 'solid'" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}




/*
// For manual checks 

int main()
{
	double kT = 0.8;
	double rhoL = 0.5;
	
	double rhoS_p = solidDensityAtSamePressure(kT, rhoL);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoL);
	
	cout << "kT = " << kT << endl;
	cout << "rhoL = " << rhoL << endl;
	cout << "rhoS_p  = " << rhoS_p << endl;
	cout << "rhoS_mu = " << rhoS_mu << endl;
	
	return 0;
}
*/




