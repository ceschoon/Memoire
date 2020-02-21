#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <chrono>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"
#include "DFT.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

#include "../SolidFCC.h"
#include "../SolidPhase.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to analyse the initial condition

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
	int Npoints = 125;
	
	// thermodynamics
	double kT = 1.0;
	double mu = -1;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	double hsd1 = -1;
	
	// density initialisation
	int ncopy = 1;
	int ncopy_range = 2;
	double Nvac = 0.0025;
	double alpha = 100;
	double alphaMin = 50;
	double alphaMax = 5000;
	double alphaStep = 0.3;
	
	////////////////////////////////////////////
	
	Options options;
	
	// control
	options.addOption("nCores", &nCores);
	options.addOption("ShowGraphics", &showGraphics);
	options.addOption("OutputFile", &outfile);
	options.addOption("IntegrationPointsFile", &pointsFile);
	options.addOption("InputFile", &infile);
	
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
	options.addOption("GaussianAlpha", &alpha);
	options.addOption("GaussianAlphaMin", &alphaMin);
	options.addOption("GaussianAlphaMax", &alphaMax);
	options.addOption("GaussianAlphaLogStep", &alphaStep);
	
	options.read(argc, argv);
	
	////////////////////////////////////////////
	
	Log log("log.dat");
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	
	options.write(log);
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Npoints = " << Npoints << endl;
	
	#ifdef USE_OMP    
	omp_set_dynamic(0);
	omp_set_num_threads(nCores);

	int fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	log << "omp max threads = " << omp_get_max_threads() << endl;
	#endif
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	////////////////////////////////////
	
	int ncopy3= ncopy*ncopy*ncopy;

	double Natoms  = 4*(1-Nvac)*ncopy3;
	double alatt   = 0;
	double Density = 0;

	if(dx > 0)
	{
		L[0] = L[1] = L[2] = Npoints*dx;
		alatt = L[0]/ncopy;
		Density = 4*(1-Nvac)/(alatt*alatt*alatt);
	} else {
		dx = L[0]/Npoints;
		alatt = pow(4*(1-Nvac)/Density,1.0/3.0);
		L[0] = L[1] = L[2] = ncopy*alatt;
	}
	
	chrono::system_clock::time_point start;
	chrono::system_clock::time_point stop;
	
	
	///////////////////////  Initialise the System /////////////////////////


	// Create potential && effective hsd

	LJ potential1(sigma1, eps1, rcut1);

	if(hsd1 < 0) hsd1 = potential1.getHSD(kT);

	// Create density objects


	SolidFCC theDensity1(dx, L, hsd1);
	
	// time the initialisation
	start = chrono::system_clock::now();
	
	theDensity1.initializeGaussian(alpha, alatt, ncopy, 1-Nvac);
	//theDensity1.initializeGaussian2(alpha, alatt, ncopy_range, 1-Nvac);
	
	stop = chrono::system_clock::now();
	chrono::duration<double> time_init = stop-start;

	// Create the species objects

	FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile, mu, Npoints);

	Interaction i1(species1,species1,potential1,kT,log,pointsFile);
	i1.initialize();

	// Create the hard-sphere object
	RSLT fmt;

	// DFT object

	DFT dft(&species1);

	dft.addHardCoreContribution(&fmt);  
	dft.addInteraction(&i1);

	// Report
	log << "Lx = " << L[0] << endl;
	log << "HSD 1 = " << hsd1 << endl;
	
	// Fix the chemical potential of the surfactant species.

	species1.setChemPotential(mu);
	
	
	
	/////////////////////// Analyse initial condition //////////////////////
	
	
	double expectedNumParticles = 4*(1-Nvac);
	
	// effect of the small constant added to avoid EtaTooLargeException
	double numParticlesSmallConstant = 1e-6*dx*dx*dx*Npoints*Npoints*Npoints*1.0/ncopy3;
	
	double finalNumParticles = theDensity1.getNumberAtoms()/ncopy3;
	
	// number of particles lost due to the gaussian extending to infinity
	double numParticlesCutOff = expectedNumParticles - 
		(finalNumParticles - numParticlesSmallConstant);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << "alpha = " << alpha << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << "Results normalised for one cell: " << endl;
	log << "expectedNumParticles = " << expectedNumParticles << endl;
	log << "numParticlesSmallConstant = " << numParticlesSmallConstant << endl;
	log << "finalNumParticles = " << finalNumParticles << endl;
	log << "numParticlesCutOff = " << numParticlesCutOff << endl;
	
	
	//////////////////////// Test a DFT computation ////////////////////////
	
	start = chrono::system_clock::now();
	
	double freeEnergy = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
	
	stop = chrono::system_clock::now();
	chrono::duration<double> time_dft = stop-start;
	
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << "freeEnergy = " << freeEnergy << endl;
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << "time_init = " << time_init.count() << endl;
	log << "time_dft = " << time_dft.count() << endl;
	
	return 0;
}
