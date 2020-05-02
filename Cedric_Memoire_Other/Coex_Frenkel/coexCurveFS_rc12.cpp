
// Program to compute the fluid-solid coexistence curve from Frenkel's
// fits results from https://arxiv.org/pdf/1910.05746.pdf




#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Frenkel_WHDF_rc12.h"

using namespace std;




// For the rc=1.2 FS coexistence, the initial guess for the densities needs
// to be more sophisticated.

double guessF(double kT)
{
	double rho = 0.9;
	
	if (kT < 0.70) rho = 0.5;
	if (kT < 0.50) rho = 0.01;
	if (kT < 0.30) rho = 1e-6;
	if (kT < 0.20) rho = 1e-10;
	if (kT < 0.15) rho = 1e-15;
	if (kT < 0.11) rho = 1e-20;
	
	return rho;
}




double solidDensityAtSamePressure(double kT, double rhoF)
{
	double rhoS = 1.2; // initial guess
	double pressureS = pressureRc12Solid(rhoS, kT);
	double pressureF = pressureRc12Fluid(rhoF, kT);
	
	double rho_step = 0.01; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and pressure
		double rhoS_new = rhoS + rho_step;
		double pressureS_new = pressureRc12Solid(rhoS_new, kT);
		
		// keep new density if it decreases the pressure difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(pressureS_new-pressureF) < abs(pressureS-pressureF))
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




double solidDensityAtSameMu(double kT, double rhoF)
{
	double rhoS = 1.2; // initial guess
	double muS = muRc12Solid(rhoS, kT);
	double muF = muRc12Fluid(rhoF, kT);
	
	double rho_step = 0.01; // initial step
	double rho_precision = 1e-6;
	
	while ( abs(rho_step) > rho_precision )
	{
		// new density and mu
		double rhoS_new = rhoS + rho_step;
		double muS_new = muRc12Solid(rhoS_new, kT);
		
		// keep new density if it decreases the mu difference,
		// otherwise continue to search in the other direction with a 
		// smaller step
		if (abs(muS_new-muF) < abs(muS-muF))
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




int coexistenceFS(double kT, double &rhoF, double &rhoS)
{
	rhoF = guessF(kT); // initial density
	double rhoS_p  = solidDensityAtSamePressure(kT, rhoF);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoF);
	
	double rho_step = 0.03*guessF(kT); // initial step
	double rho_relprecision = 1e-5;
	
	while ( abs(rho_step/rhoF)           > rho_relprecision )//||
	        //abs((rhoS_p-rhoS_mu)/rhoS_p) > rho_relprecision    )
	{
		// new densities
		double rhoF_new = rhoF + rho_step;
		double rhoS_p_new  = solidDensityAtSamePressure(kT, rhoF_new);
		double rhoS_mu_new = solidDensityAtSameMu(kT, rhoF_new);
		
		// keep new density if it decreases the difference between the 
		// corresponding solid densities, otherwise continue to search in 
		// the other direction with a smaller step
		if (abs(rhoS_p_new-rhoS_mu_new) < abs(rhoS_p-rhoS_mu))
		{
			rhoF = rhoF_new;
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
	
	ofstream outDataFile("coexCurveFS_Frenkel.dat");
	outDataFile << "#kT				rhoF			rhoS			"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	// Compute coexistence curve
	
	for (double kT=0.1; kT<=1.4; kT+=0.001)
	{
		double rhoF,rhoS;
		int status = coexistenceFS(kT, rhoF, rhoS);
		
		if (status==0) 
		{
			outDataFile << kT << " 	" << rhoF << " 	" << rhoS << " 	"
			            << endl;
		}
		else
		{
			cout << "ERROR: coexistenceFS failed for kT = " << kT << endl;
		}
	}
	
	// Gnuplot script
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexCurveFS_Frenkel.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Fluid-Solid (Frenkel)\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveFS_Frenkel.dat' using 2:1 title 'fluid', \\" << endl;
	outPlotFile << "     'coexCurveFS_Frenkel.dat' using 3:1 title 'solid'" << endl;
	
	int sysresult = system("gnuplot plot");
	
	return 0;
}




/*
// For manual checks 

int main()
{
	double kT = 0.8;
	double rhoF = 0.5;
	
	double rhoS_p = solidDensityAtSamePressure(kT, rhoF);
	double rhoS_mu = solidDensityAtSameMu(kT, rhoF);
	
	cout << "kT = " << kT << endl;
	cout << "rhoF = " << rhoF << endl;
	cout << "rhoS_p  = " << rhoS_p << endl;
	cout << "rhoS_mu = " << rhoS_mu << endl;
	
	return 0;
}
*/




