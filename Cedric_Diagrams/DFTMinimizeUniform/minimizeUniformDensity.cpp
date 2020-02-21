#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <gsl/gsl_sf_lambert.h>
#include <armadillo>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "../SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the uniform free energy and density at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// control
	int nCores = 6;
	bool showGraphics = true;
	string pointsFile("..//SS31-Mar-2016//ss109.05998");
	string outfile("dump.dat");
	string infile;
	
	// geometry
	double dx = 0.1;
	double L[3] = {10,10,10};
	int Npoints = 100;
	
	// thermodynamics
	double kT = 1;
	double mu = -1;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	double hsd1 = -1;

	// density initialisation
	double rho = 0.3;

	// minimisation
	double forceLimit = 1e-4;
	double dt = 1e-3;
	double dtMax = 1;
	double alpha_start = 0.01;
	double alphaFac = 1;
	int maxSteps = -1;
	
	Options options;
	
	// control
	options.addOption("nCores", &nCores);
	options.addOption("ShowGraphics", &showGraphics);
	options.addOption("OutputFile", &outfile);
	options.addOption("IntegrationPointsFile", &pointsFile);
	options.addOption("InputFile", &infile);

	// geometry
	options.addOption("Dx", &dx);
	options.addOption("Npoints", &Npoints);
	
	// thermodynamics
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	
	// density initialisation
	options.addOption("Rho", &rho);

	// minimisation
	options.addOption("ForceTerminationCriterion",&forceLimit);
	options.addOption("TimeStep", &dt);
	options.addOption("TimeStepMax", &dtMax);
	options.addOption("AlphaStart", &alpha_start);
	options.addOption("AlphaFac", &alphaFac);
	options.addOption("MaxSteps", &maxSteps);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	#ifdef USE_OMP    
	omp_set_dynamic(0);
	omp_set_num_threads(nCores);

	int fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	log << "omp max threads = " << omp_get_max_threads() << endl;
	#endif
	
	L[0] = L[1] = L[2] = Npoints*dx;


	////////////////////////  Initialise the System //////////////////////////


	//////////////////////////////////////
	////// Create potential && effective hsd

	LJ potential1(sigma1, eps1, rcut1);

	if(hsd1 < 0) hsd1 = potential1.getHSD(kT);

	/////////////////////////////////////
	// Create density objects

	SolidFCC theDensity1(dx, L, hsd1);
	
	theDensity1.initializeUniform(rho);

	/////////////////////////////////////
	// Create the species objects

	FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile);

	Interaction i1(species1,species1,potential1,kT,log,pointsFile);
	i1.initialize();

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


	if(! infile.empty())
	theDensity1.readDensity(infile.c_str());


	///////////////////////////////////////////////
	// Fix the chemical potential

	//    species1.setFixedMass(N);
	//    log << "N fixed at " << N << endl;
	//    species1.setChemPotential(0);

	species1.setChemPotential(mu);


	log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;

	/////////////////////////////  Minimisation //////////////////////////////

	fireMinimizer2 minimizer(dft,log);
	minimizer.setForceTerminationCriterion(forceLimit);
	minimizer.setTimeStep(dt);
	minimizer.setTimeStepMax(dtMax);
	minimizer.setAlphaStart(alpha_start);
	minimizer.setAlphaFac(alphaFac);
	minimizer.run(maxSteps);

	/////////////////////////////  Final Results /////////////////////////////

	double Natoms = theDensity1.getNumberAtoms();
	double Omega = minimizer.getF();

	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Final Density: " << Natoms/(L[0]*L[1]*L[2]) << endl;
	log << "Final Free Energy: " << Omega/(L[0]*L[1]*L[2]) << endl;
	
	return 1;
}

