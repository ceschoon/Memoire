#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <stdlib.h>

using namespace std; 

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "utilities.h"


void DFTcomputation( int argc, char** argv, Log &log, 
                     double kT, double mu, int Npoints, 
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, double &aVdW, double &hsd, 
                     bool useSnapshot, bool &success );

void DFTcomputation2( int argc, char** argv, Log &log, 
                      double kT, double mu, int Npoints,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, double &aVdW, double &hsd, 
                      bool useSnapshot, bool &success );

void DFTGaussianProfiles( Log &log, DFT &dft, SolidFCC &theDensity1,
                          double *L, int ncopy, double alatt, 
                          double NvacMin, double NvacMax, double NvacStep,
                          double alphaMin, double alphaMax, double alphaStep,
                          double &density_final, double &numParticles_final,
                          double &freeEnergy_final, bool &success);

void DFTGaussianProfilesFixedNvac( Log &log, DFT &dft, SolidFCC &theDensity1,
                          double *L, int ncopy, double alatt, double Nvac,
                          double alphaMin, double alphaMax, double alphaStep,
                          double &density_final, double &numParticles_final,
                          double &freeEnergy_final, bool &success);

void DFTGaussianProfilesMin2D( Log &log, DFT &dft, SolidFCC &theDensity1,
                               double *L, int ncopy, double alatt, string inFileName, 
                               double NvacMin, double NvacMax, 
                               double NvacLogStepMin, double NvacLogStepMax,
                               double alphaMin, double alphaMax, 
                               double alphaLogStepMin, double alphaLogStepMax,
                               double &density_final, double &numParticles_final,
                               double &freeEnergy_final, bool &success);




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// This function is just an intermediate to catch EtaTooLarge exceptions

void DFTcomputation( int argc, char** argv, Log &log, 
                     double kT, double mu, int Npoints, 
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, double &aVdW, double &hsd, 
                     bool useSnapshot, bool &success )
{
	try
	{
		DFTcomputation2( argc, argv, log, kT, mu, Npoints, density_final,
		                 numParticles_final, freeEnergy_final, 
		                 aVdW, hsd, useSnapshot, success);
	}
	catch( Eta_Too_Large_Exception &e)
	{
		success = false;
		log << myColor::RED << "=================================" << myColor::RESET << endl;
		log << "ERROR: Eta too large in DFTcomputation" << endl;
	}
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// Remark: the free energy is given per unit volume

void DFTcomputation2( int argc, char** argv, Log &log, 
                      double kT, double mu, int Npoints,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, double &aVdW, double &hsd, 
                      bool useSnapshot, bool &success )
{
	////////////////////////////  Parameters ///////////////////////////////

	// control
	int nCores = 6;
	bool showGraphics = true;
	string pointsFile("..//SS31-Mar-2016//ss109.05998");
	string outfile("dump.dat");
	string infile;
	string infileGaussian;

	// geometry
	double dx = 0.1;

	// thermodynamics -> as input

	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	double hsd1 = -1;

	// density initialisation
	int ncopy = 1;
	int ncopy_range = 2;
	double Nvac = 0.0025;
	double NvacMin = 0.001;
	double NvacMax = 0.01;
	double NvacLogStep = 0.3;
	double NvacLogStepMin = 0.003;
	double NvacLogStepMax = 0.3;
	double alpha = 300;
	double alphaMin = 50;
	double alphaMax = 5000;
	double alphaLogStep = 0.3;
	double alphaLogStepMin = 0.003;
	double alphaLogStepMax = 0.3;

	// minimisation
	double forceLimit = 1.0;
	double dt = 1e-3;
	double dtMax = 1;
	double alpha_start = 0.01;
	double alphaFac = 1;
	int maxSteps = -1;

	////////////////////////////////////////////

	Options options;

	// control
	options.addOption("nCores", &nCores);
	options.addOption("ShowGraphics", &showGraphics);
	options.addOption("OutputFile", &outfile);
	options.addOption("IntegrationPointsFile", &pointsFile);
	if (useSnapshot) options.addOption("InputFile", &infile);
	if (useSnapshot) options.addOption("InputFileGaussian", &infileGaussian);

	// geometry
	options.addOption("Dx", &dx);

	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);

	// density initialisation
	options.addOption("ncopy", &ncopy);
	options.addOption("ncopy_range", &ncopy_range);
	options.addOption("Nvac", &Nvac);
	options.addOption("NvacMin", &NvacMin);
	options.addOption("NvacMax", &NvacMax);
	options.addOption("NvacLogStep", &NvacLogStep);
	options.addOption("NvacLogStepMin", &NvacLogStepMin);
	options.addOption("NvacLogStepMax", &NvacLogStepMax);
	options.addOption("GaussianAlpha", &alpha);
	options.addOption("GaussianAlphaMin", &alphaMin);
	options.addOption("GaussianAlphaMax", &alphaMax);
	options.addOption("GaussianAlphaLogStep", &alphaLogStep);
	options.addOption("GaussianAlphaLogStepMin", &alphaLogStepMin);
	options.addOption("GaussianAlphaLogStepMax", &alphaLogStepMax);

	// minimisation
	options.addOption("ForceTerminationCriterion",&forceLimit);
	options.addOption("TimeStep", &dt);
	options.addOption("TimeStepMax", &dtMax);
	options.addOption("AlphaStart", &alpha_start);
	options.addOption("AlphaFac", &alphaFac);
	options.addOption("MaxSteps", &maxSteps);

	options.read(argc, argv);

	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- DFT computation at fixed kT, mu, Npoints and alpha ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

	log << myColor::RED << myColor::BOLD << "Input parameters (arguments):" << myColor::RESET << endl <<  "#" << endl;

	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Npoints = " << Npoints << endl;

	log << myColor::GREEN << "=================================" << myColor::RESET << endl;

	log << myColor::RED << myColor::BOLD << "Input parameters (file):" << myColor::RESET << endl <<  "#" << endl;

	options.write(log);

	#ifdef USE_OMP    
	omp_set_dynamic(0);
	omp_set_num_threads(nCores);

	int fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	log << "omp max threads = " << omp_get_max_threads() << endl;
	#endif

	log << myColor::GREEN << "=================================" << myColor::RESET << endl;

	////////////////////////////////////

	success = true;
	
	int ncopy3= ncopy*ncopy*ncopy;

	double L[3] = {10,10,10};
	L[0] = L[1] = L[2] = Npoints*dx;
	
	double alatt = L[0]/ncopy;



	///////////////////////  Initialise the System /////////////////////////


	//////////////////////////////////////
	////// Create potential && effective hsd

	LJ potential1(sigma1, eps1, rcut1);

	if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
	hsd = hsd1;

	/////////////////////////////////////
	// Create density objects


	SolidFCC theDensity1(dx, L, hsd1);


	/////////////////////////////////////
	// Create the species objects

	FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile, mu, Npoints);

	Interaction i1(species1,species1,potential1,kT,log,pointsFile);
	i1.initialize();
	aVdW = i1.getVDWParameter()/2; // divide by 2 as not the same definition

	/////////////////////////////////////
	// Create the hard-sphere object
	RSLT fmt;

	/////////////////////////////////////
	// DFT object

	DFT dft(&species1);

	dft.addHardCoreContribution(&fmt);  
	dft.addInteraction(&i1);

	/////////////////////////////////////////////////////
	// Report
	log << "Lx = " << L[0] << endl;
	log << "HSD 1 = " << hsd1 << endl;


	///////////////////////////////////////////////
	// Fix the chemical potential of the surfactant species.

	species1.setChemPotential(mu);

	//  check(theDensity1, dft,i1);

	log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;





	////////////////////// Run with gaussian profiles //////////////////////

#ifdef RESTRICT_GAUSSIAN_PROFILE

	// We first do multiple runs with different variations of the gaussian 
	// in order to select the best profile as the initial condition of the 
	// full density minimisation. 
	// --> NOT ANYMORE
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "--- Restriction to gaussian profiles ---" << endl <<  "#" << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	double density_gaussian = 0;
	double numParticles_gaussian = 0;
	double freeEnergy_gaussian = 0;
	
	bool successGaussian = false;
	
	// Min the two parameters (alpha, Nvac) at the same time -> faster
	DFTGaussianProfilesMin2D( log, dft, theDensity1,
	                          L, ncopy, alatt, infileGaussian, 
	                          NvacMin, NvacMax, 
	                          NvacLogStepMin, NvacLogStepMax,
	                          alphaMin, alphaMax, 
	                          alphaLogStepMin, alphaLogStepMax,
	                          density_gaussian, numParticles_gaussian,
	                          freeEnergy_gaussian, successGaussian);
	
	/*
	DFTGaussianProfiles( log, dft, theDensity1,
	                     L, ncopy, alatt,
	                     NvacMin, NvacMax, NvacLogStep,
	                     alphaMin, alphaMax, alphaLogStep,
	                     density_gaussian, numParticles_gaussian,
	                     freeEnergy_gaussian, successGaussian);
	*/
	
	if (!successGaussian) success = false;
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Gaussian profile) Free Energy = " << freeEnergy_gaussian << endl;
	log << "(Gaussian profile) Density = " << density_gaussian << endl;
	log << "(Gaussian profile) Number of Particles = " << numParticles_gaussian << endl;

#endif




	//////////////////////////  Full Minimisation //////////////////////////

#ifndef RESTRICT_GAUSSIAN_PROFILE
	
	if (useSnapshot && !infile.empty())  theDensity1.initialize_from_file(infile.c_str());
	else  theDensity1.initializeGaussian(alpha, alatt, ncopy, 1-Nvac);
	
	fireMinimizer2 minimizer(dft,log);
	minimizer.setForceTerminationCriterion(forceLimit);
	minimizer.setTimeStep(dt);
	minimizer.setTimeStepMax(dtMax);
	minimizer.setAlphaStart(alpha_start);
	minimizer.setAlphaFac(alphaFac);
	minimizer.run(maxSteps);
	
	// Sanity check

	//if (minimizer.getCalls() >= maxSteps-1 && maxSteps > 0)
	// TODO: check if F is stable with err>1
	// success = false;

#endif


	////////////////////////////  Final Results ////////////////////////////

#ifdef RESTRICT_GAUSSIAN_PROFILE
	
	density_final = density_gaussian;
	numParticles_final = numParticles_gaussian;
	freeEnergy_final = freeEnergy_gaussian;
	
#endif

#ifndef RESTRICT_GAUSSIAN_PROFILE
	
	double Natoms = theDensity1.getNumberAtoms();
	double Omega = minimizer.getF();
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Final Omega/V: " << Omega/(L[0]*L[1]*L[2]) << endl;
	
	density_final = Natoms/(L[0]*L[1]*L[2]);
	numParticles_final = Natoms;
	freeEnergy_final = Omega/(L[0]*L[1]*L[2]);
	
#endif

}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Minimise over Nvac

void DFTGaussianProfiles( Log &log, DFT &dft, SolidFCC &theDensity1,
                          double *L, int ncopy, double alatt, 
                          double NvacMin, double NvacMax, double NvacLogStep,
                          double alphaMin, double alphaMax, double alphaLogStep,
                          double &density_final, double &numParticles_final,
                          double &freeEnergy_final, bool &success)
{
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "--- Minimisation over Nvac ---" << endl <<  "#" << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	///////// initialise variables
	
	success = true;
	
	int NumNvacs = std::log(NvacMax/NvacMin)/NvacLogStep+1;
	vector<double> NvacVec(NumNvacs);
	
	for(int i=0; i<NumNvacs; i++)
		NvacVec[i] = NvacMin * exp(NvacLogStep*i);
	
	vector<bool> successPerNvac(NumNvacs, true);
	vector<double> density(NumNvacs);
	vector<double> numParticles(NumNvacs);
	vector<double> freeEnergy(NumNvacs);
	
	///////// DFT computation for each Nvac
	
	for(int i=0; i<NumNvacs; i++)
	{
		try
		{
			double density_temp = 0;
			double numParticles_temp = 0;
			double freeEnergy_temp = 0;
			bool success_temp = false;
			
			DFTGaussianProfilesFixedNvac(log, dft, theDensity1, 
				L, ncopy, alatt, NvacVec[i],
				alphaMin, alphaMax, alphaLogStep, density_temp,
				numParticles_temp, freeEnergy_temp, success_temp);
			
			numParticles[i] = numParticles_temp;
			density[i] = density_temp;
			freeEnergy[i] = freeEnergy_temp;
			successPerNvac[i] = success_temp;
		}
		catch(...)
		{
			successPerNvac[i] = false;
		}
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "NvacVec = ";
	for (int i=0;i<NumNvacs; i++) log << NvacVec[i] << " ";
	log << endl;
	
	log << "freeEnergy = ";
	for (int i=0;i<NumNvacs; i++) log << freeEnergy[i] << " ";
	log << endl;
	
	log << "density = ";
	for (int i=0;i<NumNvacs; i++) log << density[i] << " ";
	log << endl;
	
	log << "numParticles = ";
	for (int i=0;i<NumNvacs; i++) log << numParticles[i] << " ";
	log << endl;
	
	log << "successPerNvac = ";
	for (int i=0;i<NumNvacs; i++) log << successPerNvac[i] << " ";
	log << endl;
	
	
	//////// We find the best Nvac by minimizing the free energy 
	// We do the minimisation of a parabolic fit near the min data point
	// First, we clean the dataset of values from failed computations
	
	vector<double> NvacVec_clean;
	vector<double> freeEnergy_clean;
	vector<double> density_clean;
	vector<double> numParticles_clean;
	
	for (int i=0; i<NumNvacs; i++)
	{
		if (successPerNvac[i]) 
		{
			NvacVec_clean.push_back(NvacVec[i]);
			freeEnergy_clean.push_back(freeEnergy[i]);
			density_clean.push_back(density[i]);
			numParticles_clean.push_back(numParticles[i]);
		}
	}
	
	NvacVec = NvacVec_clean;
	freeEnergy = freeEnergy_clean;
	density = density_clean;
	numParticles = numParticles_clean;
	
	double Nvac_min = 0;
	double freeEnergy_min = 0;
	double density_min = 0;
	double numParticles_min = 0;
	
	int statusMinNvac = minFromDataParabola(NvacVec, freeEnergy, Nvac_min, 
	                                       freeEnergy_min);
	
	if (statusMinNvac!=0) success = false;
	
	/////////
	
	if (statusMinNvac==0)
	{
		evalFromDataInterpolation(NvacVec, density, Nvac_min, density_min);
		evalFromDataInterpolation(NvacVec, numParticles, Nvac_min, numParticles_min);
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min Nvac) Min Nvac = " << Nvac_min << endl;
	log << "(Min Nvac) Min density = " << density_min << endl;
	log << "(Min Nvac) Min numParticles = " << numParticles_min << endl;
	log << "(Min Nvac) Min freeEnergy = " << freeEnergy_min << endl;
	
	// Return values
	
	density_final = density_min;
	numParticles_final = numParticles_min;
	freeEnergy_final = freeEnergy_min;
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Minimise over alpha

void DFTGaussianProfilesFixedNvac( Log &log, DFT &dft, SolidFCC &theDensity1,
                          double *L, int ncopy, double alatt, double Nvac,
                          double alphaMin, double alphaMax, double alphaLogStep,
                          double &density_final, double &numParticles_final,
                          double &freeEnergy_final, bool &success)
{
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "--- Minimisation over alpha ---" << endl <<  "#" << endl;
	log << "Nvac = " << Nvac << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	///////// initialise variables
	
	success = true;
	
	int NumAlphas = std::log(alphaMax/alphaMin)/alphaLogStep+1;
	vector<double> alphaVec(NumAlphas);
	
	for(int i=0; i<NumAlphas; i++)
		alphaVec[i] = alphaMin * exp(alphaLogStep*i);
	
	vector<bool> successPerAlpha(NumAlphas, true);
	vector<double> density(NumAlphas);
	vector<double> numParticles(NumAlphas);
	vector<double> freeEnergy(NumAlphas);
	
	///////// DFT computation for each alpha
	
	for(int i=0; i<NumAlphas; i++)
	{
		try
		{
			theDensity1.initializeGaussian(alphaVec[i], alatt, ncopy, 1-Nvac);
			
			numParticles[i] = theDensity1.getNumberAtoms();
			density[i] = numParticles[i]/(L[0]*L[1]*L[2]);
			freeEnergy[i] = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
		}
		catch(...)
		{
			successPerAlpha[i] = false;
		}
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "alphaVec = ";
	for (int i=0;i<NumAlphas; i++) log << alphaVec[i] << " ";
	log << endl;
	
	log << "freeEnergy = ";
	for (int i=0;i<NumAlphas; i++) log << freeEnergy[i] << " ";
	log << endl;
	
	log << "density = ";
	for (int i=0;i<NumAlphas; i++) log << density[i] << " ";
	log << endl;
	
	log << "numParticles = ";
	for (int i=0;i<NumAlphas; i++) log << numParticles[i] << " ";
	log << endl;
	
	log << "successPerAlpha = ";
	for (int i=0;i<NumAlphas; i++) log << successPerAlpha[i] << " ";
	log << endl;
	
	
	//////// We find the best alpha by minimizing the free energy 
	// We do the minimisation of a parabolic fit near the min data point
	// First, we clean the dataset of values from failed computations
	
	vector<double> alphaVec_clean;
	vector<double> freeEnergy_clean;
	vector<double> density_clean;
	vector<double> numParticles_clean;
	
	for (int i=0; i<NumAlphas; i++)
	{
		if (successPerAlpha[i]) 
		{
			alphaVec_clean.push_back(alphaVec[i]);
			freeEnergy_clean.push_back(freeEnergy[i]);
			density_clean.push_back(density[i]);
			numParticles_clean.push_back(numParticles[i]);
		}
	}
	
	alphaVec = alphaVec_clean;
	freeEnergy = freeEnergy_clean;
	density = density_clean;
	numParticles = numParticles_clean;
	
	double alpha_min = 0;
	double freeEnergy_min = 0;
	double density_min = 0;
	double numParticles_min = 0;
	
	int statusMinAlpha = minFromDataParabola(alphaVec, freeEnergy, alpha_min, 
	                                         freeEnergy_min);
	
	if (statusMinAlpha!=0) success = false;
	
	/////////
	
	if (statusMinAlpha==0)
	{
		evalFromDataInterpolation(alphaVec, density, alpha_min, density_min);
		evalFromDataInterpolation(alphaVec, numParticles, alpha_min, numParticles_min);
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min alpha) Min alpha = " << alpha_min << endl;
	log << "(Min alpha) Min density = " << density_min << endl;
	log << "(Min alpha) Min numParticles = " << numParticles_min << endl;
	log << "(Min alpha) Min freeEnergy = " << freeEnergy_min << endl;
	
	// Return values
	
	density_final = density_min;
	numParticles_final = numParticles_min;
	freeEnergy_final = freeEnergy_min;
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Minimise over alpha and Nvac at the same time

void DFTGaussianProfilesMin2D( Log &log, DFT &dft, SolidFCC &theDensity1,
                               double *L, int ncopy, double alatt, string inFileName, 
                               double NvacMin, double NvacMax, 
                               double NvacLogStepMin, double NvacLogStepMax,
                               double alphaMin, double alphaMax, 
                               double alphaLogStepMin, double alphaLogStepMax,
                               double &density_final, double &numParticles_final,
                               double &freeEnergy_final, bool &success)
{
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "--- Minimisation over alpha and Nvac ---" << endl <<  "#" << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	///////////////////////// Initialise variables /////////////////////////
	
	success = true;
	bool foundMinimum = false;
	
	double alphaLogStep = alphaLogStepMax;
	double NvacLogStep = NvacLogStepMax;
	
	double alpha_current;
	double Nvac_current;
	double freeEnergy_current;
	
	// max number of dft computations allowed in this function 
	// above this limit, I don't consider that this function is useful
	
	int maxDFTcomputations = int (std::log(alphaMax/alphaMin)/alphaLogStepMax)
	                       * int (std::log( NvacMax/ NvacMin)/NvacLogStepMax);
	int numDFTcomputations = 0;
	
	// Random number engine to choose starting position in (alpha,Nvac)
	
	int seed = chrono::system_clock::now().time_since_epoch().count();
	log << "seed = " << seed << endl;
	
	default_random_engine gen(seed);
	uniform_real_distribution<double> dist01(0,1);
	
	bool initialisePosition = true;
	
	// initialise from file if inFile non empty
	
	if(! inFileName.empty())
	{
		ifstream inFile(inFileName);
		
		if (inFile)
		{
			readDataFromFile(inFile, "alpha", alpha_current);
			readDataFromFile(inFile, "Nvac", Nvac_current);
			
			try
			{
				theDensity1.initializeGaussian(alpha_current, alatt, ncopy, 1-Nvac_current);
				freeEnergy_current = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
				numDFTcomputations ++;
				
				initialisePosition = false;
			}
			catch (...)
			{
				log << "DFT computation failed for inFile initial condition" << endl;
			}
		}
		else
		{
			log << "ERROR: Unable to read inFile for initial condition" << endl;
		}
	}
	
	///////////////////////////// Minimisation /////////////////////////////
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Starting minimisation" << endl;
	
	int indexPrev = 0; // previous direction with lower free energy
	
	while (!foundMinimum && numDFTcomputations<=maxDFTcomputations)
	{
		// Starting position in (alpha,Nvac) plane initialised randomly 
		
		while (initialisePosition)
		{
			try
			{
				double ran1 = dist01(gen);
				double ran2 = dist01(gen);
				
				alpha_current = exp(   std::log(alphaMin)*(1-ran1)
				                     + std::log(alphaMax)*ran1      );
				Nvac_current  = exp(   std::log(NvacMin )*(1-ran1)
				                     + std::log(NvacMax )*ran1      );
				
				log << "Trying initial condition alpha = " << alpha_current
					<< " and Nvac = " << Nvac_current << endl;
				
				// initial free energy
				theDensity1.initializeGaussian(alpha_current, alatt, ncopy, 1-Nvac_current);
				freeEnergy_current = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
				numDFTcomputations ++;
				
				initialisePosition = false;
				
				log << "Found successful initial condition" << endl;
			}
			catch(...) 
			{
				numDFTcomputations ++;
			}
		}
		
		// compute the free energy in every direction
		
		vector<double> alpha_candidates(4);
		vector<double> Nvac_candidates(4);
		vector<double> freeEnergy_candidates(4);
		vector<bool> valid_candidates(4);
		
		alpha_candidates[0] = alpha_current*exp(alphaLogStep);
		Nvac_candidates[0] = Nvac_current;
		alpha_candidates[1] = alpha_current;
		Nvac_candidates[1] = Nvac_current*exp(NvacLogStep);
		alpha_candidates[2] = alpha_current*exp(-alphaLogStep);
		Nvac_candidates[2] = Nvac_current;
		alpha_candidates[3] = alpha_current;
		Nvac_candidates[3] = Nvac_current*exp(-NvacLogStep);
		
		int indexMin = -1; // -1 means current position
		double freeEnergyMin = freeEnergy_current;
		
		for (int i=0; i<4; i++)
		{
			// begin searching a lower free energy in the previous working
			// direction
			
			int j = indexPrev + i;
			while (j>=4) j-=4;
			
			try
			{
				theDensity1.initializeGaussian(alpha_candidates[j],
					alatt, ncopy, 1-Nvac_candidates[j]);
				freeEnergy_candidates[j] = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
				numDFTcomputations ++;
				valid_candidates[j] = true;
				
				if(freeEnergy_candidates[j] < freeEnergyMin)
				{
					indexMin = j;
					freeEnergyMin = freeEnergy_candidates[j];
					indexPrev = j;
					break; // for speed, go in the first direction better than the current position
				}
			}
			catch(...) 
			{
				numDFTcomputations ++;
				valid_candidates[j] = false;
			}
		}
		
		// Report
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		
		log << "alpha_current = " << alpha_current << endl;
		log << "Nvac_current = " << Nvac_current << endl;
		log << "freeEnergy_current = " << freeEnergy_current << endl;
		
		// Check if current position is inside of the limits
		// If not, restart from somewhere else
		
		if (alpha_current<alphaMin || alpha_current>alphaMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Went outside of the limits, retry with a new starting point" << endl;
			
			// try with a smaller step (less likely to make the same mistake)
			if (alphaLogStep>alphaLogStepMin) 
			{
				alphaLogStep = alphaLogStep/2;
				log << "dividing alphaLogStep by 2" << endl;
			}
			
			initialisePosition = true;
			continue;
		}
		if (Nvac_current<NvacMin || Nvac_current>NvacMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Went outside of the limits, retry with a new starting point" << endl;
			
			// try with a smaller step (less likely to make the same mistake)
			if (NvacLogStep>NvacLogStepMin)
			{
				NvacLogStep  = NvacLogStep/2;
				log << "dividing alphaLogStep by 2" << endl;
			}
			
			initialisePosition = true;
			continue;
		}
		
		// If min freeEnergy is the current position, stay there
		// and try again with a smaller step
		
		if (indexMin==-1)
		{
			if (alphaLogStep>alphaLogStepMin) 
			{
				alphaLogStep = alphaLogStep/4;
				log << "dividing alphaLogStep by 4" << endl;
			}
			if (NvacLogStep>NvacLogStepMin)
			{
				NvacLogStep  = NvacLogStep/4;
				log << "dividing alphaLogStep by 4" << endl;
			}
			
			// If we reach the min steps allowed, we return the current 
			// location as the minimum
			
			if (alphaLogStep<alphaLogStepMin && NvacLogStep<NvacLogStepMin) 
				foundMinimum = true;
		}
		
		// Otherwise, move to the position with lowest freeEnergy
		
		else
		{
			alpha_current = alpha_candidates[indexMin];
			Nvac_current = Nvac_candidates[indexMin];
			freeEnergy_current = freeEnergy_candidates[indexMin];
		}
	}
	
	/////////////////////////////// Finalise ///////////////////////////////
	
	// Exit if failure
	
	if (numDFTcomputations>maxDFTcomputations)
	{
		log << "ERROR: Exceeded max number of dft computations allowed" << endl;
		log << "maxDFTcomputations = " << maxDFTcomputations << endl;
		success = false;
		return;
	}
	else
	{
		log << "numDFTcomputations = " << numDFTcomputations << endl;
		log << "maxDFTcomputations = " << maxDFTcomputations << endl;
	}
	
	// Last dft computation to get the values of interest at the minimum
	
	try
	{
		theDensity1.initializeGaussian(alpha_current, 
			alatt, ncopy, 1-Nvac_current);
		
		freeEnergy_final = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
		numParticles_final = theDensity1.getNumberAtoms();
		density_final = numParticles_final/(L[0]*L[1]*L[2]);
	}
	catch(...) 
	{
		log << "Failed last computation" << endl;
	}
	
	// save in file for reuse
	
	ofstream outFile("gaussian_profile_min.dat");
	if (outFile)
	{
		outFile << "alpha = " << alpha_current << endl;
		outFile << "Nvac = " << Nvac_current << endl;
		outFile << "density = " << density_final << endl;
		outFile << "numParticles = " << numParticles_final << endl;
		outFile << "freeEnergy = " << freeEnergy_final << endl;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min alpha,Nvac) Min alpha = " << alpha_current << endl;
	log << "(Min alpha,Nvac) Min Nvac = " << Nvac_current << endl;
	log << "(Min alpha,Nvac) Min density = " << density_final << endl;
	log << "(Min alpha,Nvac) Min numParticles = " << numParticles_final << endl;
	log << "(Min alpha,Nvac) Min freeEnergy = " << freeEnergy_final << endl;
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/*
void check(Density &density, DFT &dft, Interaction &i1)
{

  int Nx = density.Nx();
  int Ny = density.Ny();
  int Nz = density.Nz();
  
  // Find largest density
  int jx = 0;
  int jy = 0;
  int jz = 0;
  double dmax = 0;

  for(int ix = 0; ix < Nx; ix++)
    for(int iy = 0; iy < Ny; iy++)
      for(int iz = 0; iz < Nz; iz++)
        if(density.getDensity(ix,iy,iz) > dmax)
          {
            dmax = density.getDensity(ix,iy,iz);
            jx = ix;
            jy = iy;
            jz = iz;
          }


  double  F = dft.calculateFreeEnergyAndDerivatives(false);
  DFT_Vec dF;
  dF.zeros(dft.getDF(0).size());
  dF.set(dft.getDF(0));

    
  for(int iz = jz; iz < Nz; iz++)
    {
      double eps = 1e-6;

      double d = density.getDensity(jx,jy,iz);

      double test = dF.get(density.pos(jx,jy,iz));


      double test1 = F+eps*d*test+0.5*eps*eps*d*d*i1.getW(0);
      double test2 = F-eps*d*test+0.5*eps*eps*d*d*i1.getW(0);
      
      

      double d1 = d+eps; //*d;
      density.set_Density_Elem(jx,jy,iz,d1);
      double fp = dft.calculateFreeEnergyAndDerivatives(false);
      cout << "Test1 = " << test1 << endl;

      double d2 = d-eps; //*d;
      density.set_Density_Elem(jx,jy,iz,d2);
      double fm = dft.calculateFreeEnergyAndDerivatives(false);
      cout << "Test2 = " << test2 << endl;
      
      density.set_Density_Elem(jx,jy,iz,d);

      double direct = i1.checkCalc(jx,jy,iz);


      
      cout << jx << " " << jy << " " << iz << " " << d << " " << ((fp-fm)/(d1-d2)) << " " << dF.get(density.pos(jx,jy,iz)) << " direct = " << direct << endl;
      
    }
  exit(0);  

}
*/

