#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>

#include "../utilities.h"
#include "options.h"
#include "Log.h"
#include "myColor.h"


using namespace std; 

int max(int a, int b) {return (a>b)?a:b;}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the triple point (kT and other properties) from
// the intersection of the polynomial fits of the solid density.

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string file_VS_rhoV = "You must specify the data files in the options";
	string file_VS_rhoS = "You must specify the data files in the options";
	string file_VS_mu = "You must specify the data files in the options";
	string file_VS_omega = "You must specify the data files in the options";
	string file_LS_rhoL = "You must specify the data files in the options";
	string file_LS_rhoS = "You must specify the data files in the options";
	string file_output = "triple_point.dat";
	double estimated_kT_triple = 1;
	
	Options options;
	options.addOption("file_VS_rhoV", &file_VS_rhoV);
	options.addOption("file_VS_rhoS", &file_VS_rhoS);
	options.addOption("file_VS_mu", &file_VS_mu);
	options.addOption("file_VS_omega", &file_VS_omega);
	options.addOption("file_LS_rhoL", &file_LS_rhoL);
	options.addOption("file_LS_rhoS", &file_LS_rhoS);
	options.addOption("estimated_kT_triple", &estimated_kT_triple);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////////// Compute kT at triple point ///////////////////////
	
	vector<double> coefs_VS_rhoS;
	vector<double> coefs_LS_rhoS;
	int orderPoly1;
	int orderPoly2;
	
	ifstream dataFile(file_VS_rhoS);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_VS_rhoS.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_VS_rhoS << endl;
	}
	dataFile.close();
	
	dataFile.open(file_LS_rhoS);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_LS_rhoS.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_LS_rhoS << endl;
	}
	dataFile.close();
	
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "coefs_VS_rhoS = ";
	for (int i=0; i<coefs_VS_rhoS.size(); i++) log << coefs_VS_rhoS[i] << " ";
	log << endl;
	log << "coefs_LS_rhoS = ";
	for (int i=0; i<coefs_LS_rhoS.size(); i++) log << coefs_LS_rhoS[i] << " ";
	log << endl;
	
	// Solve polynomial equation "poly1-poly2 = 0"
	
	#define N max(coefs_VS_rhoS.size(), coefs_LS_rhoS.size())
	double a[N];
	double z[2*(N-1)];
	
	for (int i=0; i<N; i++) 
	{
		double ai = 0;
		if (i<coefs_VS_rhoS.size()) ai += coefs_VS_rhoS[i];
		if (i<coefs_LS_rhoS.size()) ai -= coefs_LS_rhoS[i];
		
		a[i] = ai;
	}

	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (N);
	gsl_poly_complex_solve (a, N, w, z);
	gsl_poly_complex_workspace_free (w);
	
	// Find the real solutions
	
	double small_value = 1e-10;
	vector<double> real_roots;
	
	cout << "real roots: ";
	for (int i=0; i<N-1; i++)
	{
		if (abs(z[2*i+1]) < small_value) // is a real root
		{
			real_roots.push_back(z[2*i]);
			cout << z[2*i] << " " << endl;
		}
	}
	cout << endl;
	
	if (real_roots.size()==0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: no real root found "<< endl;
	}
	
	// Select the real root that is the closest to the given estimation
	
	double kT_triple;
	
	for (int i=0; i<real_roots.size(); i++)
	{
		double root = real_roots[i];
		if (abs(root-estimated_kT_triple) < abs(kT_triple-estimated_kT_triple))
			kT_triple = root;
	}
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Best estimation for triple point kT:     " << kT_triple << endl;
	
	
	//////////////////// Compute rhoS at triple point //////////////////////
	
	double rhoS_triple = 0;
	
	for (int i=0; i<=coefs_VS_rhoS.size(); i++)
		rhoS_triple += coefs_VS_rhoS[i] * pow(kT_triple,i);
	
	log << "Best estimation for triple point rhoS:   " << rhoS_triple << endl;
	
	
	//////////////////// Compute rhoV at triple point //////////////////////
	
	vector<double> coefs_VS_rhoV;
	
	dataFile.open(file_VS_rhoV);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_VS_rhoV.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_VS_rhoV << endl;
	}
	dataFile.close();
	
	double rhoV_triple = 0;
	
	for (int i=0; i<=coefs_VS_rhoV.size(); i++)
		rhoV_triple += coefs_VS_rhoV[i] * pow(kT_triple,i);
	
	log << "Best estimation for triple point rhoV:   " << rhoV_triple << endl;
	
	
	//////////////////// Compute rhoL at triple point //////////////////////
	
	vector<double> coefs_LS_rhoL;
	
	dataFile.open(file_LS_rhoL);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_LS_rhoL.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_LS_rhoL << endl;
	}
	dataFile.close();
	
	double rhoL_triple = 0;
	
	for (int i=0; i<=coefs_LS_rhoL.size(); i++)
		rhoL_triple += coefs_LS_rhoL[i] * pow(kT_triple,i);
	
	log << "Best estimation for triple point rhoL:   " << rhoL_triple << endl;
	
	///////////////////// Compute mu at triple point ///////////////////////
	
	vector<double> coefs_VS_mu;
	
	dataFile.open(file_VS_mu);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_VS_mu.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_VS_mu << endl;
	}
	dataFile.close();
	
	double mu_triple = 0;
	
	for (int i=0; i<=coefs_VS_mu.size(); i++)
		mu_triple += coefs_VS_mu[i] * pow(kT_triple,i);
	
	log << "Best estimation for triple point mu:     " << mu_triple << endl;
	
	//////////////////// Compute omega at triple point //////////////////////
	
	vector<double> coefs_VS_omega;
	
	dataFile.open(file_VS_omega);
	if (dataFile)
	{
		int i=0;
		while (true)
		{
			double ai;
			int status = readDataFromFile(dataFile, "a"+to_string(i), ai);
			if (status==0) coefs_VS_omega.push_back(ai);
			else break;
			i += 1;
		}
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << file_VS_omega << endl;
	}
	dataFile.close();
	
	double omega_triple = 0;
	
	for (int i=0; i<=coefs_VS_omega.size(); i++)
		omega_triple += coefs_VS_omega[i] * pow(kT_triple,i);
	
	log << "Best estimation for triple point omega:  " << omega_triple << endl;
	
	
	/////////////////////// Save results in a file /////////////////////////
	
	ofstream ofile(file_output);
	ofile << "# Properties at the triple point" << endl;
	ofile << scientific << setprecision(8);
	
	ofile << "kT    = " << kT_triple    << endl;
	ofile << "rhoV  = " << rhoV_triple  << endl;
	ofile << "rhoL  = " << rhoL_triple  << endl;
	ofile << "rhoS  = " << rhoS_triple  << endl;
	ofile << "mu    = " << mu_triple    << endl;
	ofile << "omega = " << omega_triple << endl;
	
	return 0;
}








