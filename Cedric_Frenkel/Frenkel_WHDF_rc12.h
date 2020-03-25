
// Header file for functions that return thermodynamic quantities from
// fit coefficients of the excess free energy associated to Frenkel's 
// potential from https://arxiv.org/pdf/1910.05746.pdf


#include <cmath>
#include <iostream>

using namespace std;



// function that returns Frenkel's fit coefficients for the fluid state
// the potential is WHDF rc=1.2

double coefRc12Fluid(int n, int m)
{
	if (n==2)
	{
		if (m==-3) return -1.10009;
		if (m==-2) return  6.27642;
		if (m==-1) return -13.5539;
		if (m== 0) return  15.5901;
		if (m== 1) return -7.05497;
		if (m== 2) return  0.63096;
	}
	if (n==3)
	{
		if (m==-3) return  4.79833;
		if (m==-2) return -31.7208;
		if (m==-1) return  76.0381;
		if (m== 0) return -84.1657;
		if (m== 1) return  46.4876;
		if (m== 2) return -10.4596;
	}
	if (n==4)
	{
		if (m==-3) return -0.47723;
		if (m==-2) return  19.8875;
		if (m==-1) return -62.2031;
		if (m== 0) return  78.4142;
		if (m== 1) return -55.7024;
		if (m== 2) return  18.4785;
	}
	if (n==5)
	{
		if (m==-3) return -7.60729;
		if (m==-2) return  39.2723;
		if (m==-1) return -131.249;
		if (m== 0) return  205.539;
		if (m== 1) return -110.462;
		if (m== 2) return  9.10103;
	}
	if (n==6)
	{
		if (m==-3) return -22.4105;
		if (m==-2) return  65.7709;
		if (m==-1) return  66.3596;
		if (m== 0) return -315.674;
		if (m== 1) return  244.739;
		if (m== 2) return -43.4275;
	}
	if (n==7)
	{
		if (m==-3) return  56.4536;
		if (m==-2) return -228.002;
		if (m==-1) return  243.158;
		if (m== 0) return  30.6737;
		if (m== 1) return -133.904;
		if (m== 2) return  33.7645;
	}
	if (n==8)
	{
		if (m==-3) return -29.5309;
		if (m==-2) return  127.279;
		if (m==-1) return -174.625;
		if (m== 0) return  68.3272;
		if (m== 1) return  18.5640;
		if (m== 2) return -9.79081;
	}
	
	cout << "ERROR: fit coefficient n=" << n << " m=" << m 
	     << " for the WHDF rc=1.2 fluid was not found." << endl;
	
	return 0;
}











// function that returns Frenkel's fit coefficients for the solid state
// the potential is WHDF rc=1.2 and the solid has a fcc lattice structure

double coefRc12Solid(int n, int m)
{
	if (n==0)
	{
		if (m==-2) return  1585.85;
		if (m==-1) return -2137.72;
		if (m== 0) return -1114.08;
		if (m== 1) return -0.10780;
	}
	if (n==1)
	{
		if (m==-2) return  -6114.51;
		if (m==-1) return  8239.19;
		if (m== 0) return  4047.61;
		if (m== 1) return  -195.525;
	}
	if (n==2)
	{
		if (m==-2) return  9422.37;
		if (m==-1) return -12692.5;
		if (m== 0) return -5907.11;
		if (m== 1) return  682.374;
	}
	if (n==3)
	{
		if (m==-2) return -7253.89;
		if (m==-1) return  9769.46;
		if (m== 0) return  4342.62;
		if (m== 1) return -807.773;
	}
	if (n==4)
	{
		if (m==-2) return  2789.93;
		if (m==-1) return -3757.12;
		if (m== 0) return -1600.17;
		if (m== 1) return  364.210;
	}
	if (n==5)
	{
		if (m==-2) return -428.864;
		if (m==-1) return  577.560;
		if (m== 0) return  236.283;
		if (m== 1) return -46.0042;
	}
	
	cout << "ERROR: fit coefficient n=" << n << " m=" << m 
	     << " for the WHDF rc=1.2 solid was not found." << endl;
	
	return 0;
}









// function that returns the pressure of the WHDF rc=1.2 fluid
// computed from the fluid fit coefficients
// expressed in units of kT and potential's sigma


double pressureRc12Fluid(double rho, double kT)
{
	double pressure_ideal = rho;
	double pressure_excess = 0;
	
	for (int n=2; n<=8; n++)
		for (int m=-3; m<=2; m++)
			pressure_excess += (n-1)*coefRc12Fluid(n,m)*pow(rho,n)*pow(1.0/kT,m);
	
	return pressure_ideal + pressure_excess;
}



// function that returns the pressure of the WHDF rc=1.2 solid
// computed from the solid fit coefficients
// expressed in units of kT and potential's sigma


double pressureRc12Solid(double rho, double kT)
{
	double pressure_ideal = rho;
	double pressure_excess = 0;
	
	for (int n=0; n<=5; n++)
		for (int m=-2; m<=1; m++)
			pressure_excess += (n-1)*coefRc12Solid(n,m)*pow(rho,n)*pow(1.0/kT,m);
	
	return pressure_ideal + pressure_excess;
}











// function that returns the chemical potential of the WHDF rc=1.2 fluid
// computed from the fluid fit coefficients
// expressed in units of kT and potential's sigma


double muRc12Fluid(double rho, double kT)
{
	double mu_ideal = std::log(rho);
	double mu_excess = 0;
	
	for (int n=2; n<=8; n++)
		for (int m=-3; m<=2; m++)
			mu_excess += n*coefRc12Fluid(n,m)*pow(rho,n-1)*pow(1.0/kT,m);
	
	return mu_ideal + mu_excess;
}



// function that returns the chemical potential of the WHDF rc=1.2 solid
// computed from the solid fit coefficients
// expressed in units of kT and potential's sigma


double muRc12Solid(double rho, double kT)
{
	double mu_ideal = std::log(rho);
	double mu_excess = 3.0/2*std::log(1.0/kT);
	
	for (int n=0; n<=5; n++)
		for (int m=-2; m<=1; m++)
			mu_excess += n*coefRc12Solid(n,m)*pow(rho,n-1)*pow(1.0/kT,m);
	
	return mu_ideal + mu_excess;
}








