
// Header file for functions that return thermodynamic quantities from
// fit coefficients of the excess free energy associated to Frenkel's 
// potential from https://arxiv.org/pdf/1910.05746.pdf


#include <cmath>
#include <iostream>

using namespace std;



// function that returns Frenkel's fit coefficients for the liquid state
// the potential is WHDF rc=2

double coefRc2Liquid(int n, int m)
{
	if (n==2)
	{
		if (m==-3) return  0.76902;
		if (m==-2) return -3.15391;
		if (m==-1) return  2.88640;
		if (m== 0) return  5.17074;
		if (m== 1) return -10.4462;
		if (m== 2) return 0.865638;
	}
	if (n==3)
	{
		if (m==-3) return  3.26974;
		if (m==-2) return -77.7391;
		if (m==-1) return  368.594;
		if (m== 0) return -723.031;
		if (m== 1) return  646.553;
		if (m== 2) return -219.3018;
	}
	if (n==4)
	{
		if (m==-3) return -8.29352;
		if (m==-2) return  377.872;
		if (m==-1) return -1929.64;
		if (m== 0) return  3888.37;
		if (m== 1) return -3522.43;
		if (m== 2) return 1219.511;
	}
	if (n==5)
	{
		if (m==-3) return -31.8858;
		if (m==-2) return -597.874;
		if (m==-1) return  3874.72;
		if (m== 0) return -8352.09;
		if (m== 1) return  7764.53;
		if (m== 2) return -2720.454;
	}
	if (n==6)
	{
		if (m==-3) return  102.257;
		if (m==-2) return  349.298;
		if (m==-1) return -3783.50;
		if (m== 0) return  8937.04;
		if (m== 1) return -8548.36;
		if (m== 2) return 3013.089;
	}
	if (n==7)
	{
		if (m==-3) return -97.1753;
		if (m==-2) return -4.10970;
		if (m==-1) return  1783.38;
		if (m== 0) return -4738.69;
		if (m== 1) return  4672.47;
		if (m== 2) return -1647.794;
	}
	if (n==8)
	{
		if (m==-3) return  31.1298;
		if (m==-2) return -44.7404;
		if (m==-1) return -315.575;
		if (m== 0) return  985.658;
		if (m== 1) return -1004.95;
		if (m== 2) return 353.2964;
	}
	
	cout << "ERROR: fit coefficient n=" << n << " m=" << m 
	     << " for the WHDF rc=2 liquid was not found." << endl;
	
	return 0;
}










// function that returns Frenkel's fit coefficients for the vapour state
// the potential is WHDF rc=2

double coefRc2Vapour(int n, int m)
{
	if (n==2)
	{
		if (m==-3) return  19.1395;
		if (m==-2) return -92.8140;
		if (m==-1) return  151.915;
		if (m== 0) return -99.3440;
		if (m== 1) return  16.8835;
	}
	if (n==3)
	{
		if (m==-3) return -21265.4;
		if (m==-2) return  75376.4;
		if (m==-1) return  -100238;
		if (m== 0) return  59293.4;
		if (m== 1) return -13166.0;
	}
	
	cout << "ERROR: fit coefficient n=" << n << " m=" << m 
	     << " for the WHDF rc=2 vapour was not found." << endl;
	
	return 0;
}










// function that returns Frenkel's fit coefficients for the solid state
// the potential is WHDF rc=2 and the solid has a fcc lattice structure

double coefRc2Solid(int n, int m)
{
	if (n==0)
	{
		if (m==-2) return  59.9947;
		if (m==-1) return -111.081;
		if (m== 0) return  23.4348;
		if (m== 1) return  19.8311;
	}
	if (n==1)
	{
		if (m==-2) return -254.708;
		if (m==-1) return  458.412;
		if (m== 0) return -85.2926;
		if (m== 1) return -88.2034;
	}
	if (n==2)
	{
		if (m==-2) return  427.727;
		if (m==-1) return -749.243;
		if (m== 0) return  124.562;
		if (m== 1) return  154.892;
	}
	if (n==3)
	{
		if (m==-2) return -356.028;
		if (m==-1) return  607.871;
		if (m== 0) return -75.9448;
		if (m== 1) return -157.340;
	}
	if (n==4)
	{
		if (m==-2) return  147.127;
		if (m==-1) return -245.199;
		if (m== 0) return  21.1821;
		if (m== 1) return  76.6613;
	}
	if (n==5)
	{
		if (m==-2) return -24.1761;
		if (m==-1) return  39.3881;
		if (m== 0) return -1.96241;
		if (m== 1) return -12.3529;
	}
	
	cout << "ERROR: fit coefficient n=" << n << " m=" << m 
	     << " for the WHDF rc=2 solid was not found." << endl;
	
	return 0;
}









// function that returns the pressure of the WHDF rc=2 liquid
// computed from the liquid fit coefficients
// expressed in units of kT and potential's sigma


double pressureRc2Liquid(double rho, double kT)
{
	double pressure_ideal = rho;
	double pressure_excess = 0;
	
	for (int n=2; n<=8; n++)
		for (int m=-3; m<=2; m++)
			pressure_excess += (n-1)*coefRc2Liquid(n,m)*pow(rho,n)*pow(1.0/kT,m);
	
	return pressure_ideal + pressure_excess;
}



// function that returns the pressure of the WHDF rc=2 vapour
// computed from the vapour fit coefficients
// expressed in units of kT and potential's sigma


double pressureRc2Vapour(double rho, double kT)
{
	double pressure_ideal = rho;
	double pressure_excess = 0;
	
	for (int n=2; n<=3; n++)
		for (int m=-3; m<=1; m++)
			pressure_excess += (n-1)*coefRc2Vapour(n,m)*pow(rho,n)*pow(1.0/kT,m);
	
	return pressure_ideal + pressure_excess;
}



// function that returns the pressure of the WHDF rc=2 solid
// computed from the solid fit coefficients
// expressed in units of kT and potential's sigma


double pressureRc2Solid(double rho, double kT)
{
	double pressure_ideal = rho;
	double pressure_excess = 0;
	
	for (int n=0; n<=5; n++)
		for (int m=-2; m<=1; m++)
			pressure_excess += (n-1)*coefRc2Solid(n,m)*pow(rho,n)*pow(1.0/kT,m);
	
	return pressure_ideal + pressure_excess;
}











// function that returns the chemical potential of the WHDF rc=2 liquid
// computed from the liquid fit coefficients
// expressed in units of kT and potential's sigma


double muRc2Liquid(double rho, double kT)
{
	double mu_ideal = std::log(rho);
	double mu_excess = 0;
	
	for (int n=2; n<=8; n++)
		for (int m=-3; m<=2; m++)
			mu_excess += n*coefRc2Liquid(n,m)*pow(rho,n-1)*pow(1.0/kT,m);
	
	return mu_ideal + mu_excess;
}



// function that returns the chemical potential of the WHDF rc=2 vapour
// computed from the vapour fit coefficients
// expressed in units of kT and potential's sigma


double muRc2Vapour(double rho, double kT)
{
	double mu_ideal = std::log(rho);
	double mu_excess = 0;
	
	for (int n=2; n<=3; n++)
		for (int m=-3; m<=1; m++)
			mu_excess += n*coefRc2Vapour(n,m)*pow(rho,n-1)*pow(1.0/kT,m);
	
	return mu_ideal + mu_excess;
}



// function that returns the chemical potential of the WHDF rc=2 solid
// computed from the solid fit coefficients
// expressed in units of kT and potential's sigma


double muRc2Solid(double rho, double kT)
{
	double mu_ideal = std::log(rho);
	double mu_excess = 3.0/2*std::log(1.0/kT);
	
	for (int n=0; n<=5; n++)
		for (int m=-2; m<=1; m++)
			mu_excess += n*coefRc2Solid(n,m)*pow(rho,n-1)*pow(1.0/kT,m);
	
	return mu_ideal + mu_excess;
}








