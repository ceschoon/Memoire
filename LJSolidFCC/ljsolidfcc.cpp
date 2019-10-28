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
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"


void minOverNpoints( int argc, char** argv, double kT, double mu, 
                     double &aLattice_min, double &density_min, 
                     double &numParticles_min, double &freeEnergy_min,
                     bool &success);

void minOverAlpha( int argc, char** argv, double kT, double mu, int Npoints,
                   double &alpha_min, double &density_min, 
                   double &numParticles_min, double &freeEnergy_min,
                   bool &success);

void DFTcomputation( int argc, char** argv, double kT, double mu, 
                     int Npoints, double alpha, bool runMinimizer,
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, bool &success            );

void DFTcomputation2( int argc, char** argv, double kT, double mu, 
                      int Npoints, double alpha, bool runMinimizer,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, bool &success            );

void check(Density &density, DFT &dft, Interaction &i);






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT = 1;
	double mu = -1;
	
	Options options;
	options.addOption("kT", &kT);
	options.addOption("Mu", &mu);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	//////////////////////// Run DFT computations //////////////////////////
	
	double aLattice = 0;
	double densitySolid = 0;
	double numParticlesInCell = 0;
	double freeEnergySolid = 0;
	bool success = false;
	
	minOverNpoints(argc, argv, kT, mu, aLattice, densitySolid, 
	               numParticlesInCell, freeEnergySolid, success);
	
	// Report
	
	if (success)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Lattice Constant = " << aLattice << endl;
		log << "Number of Particles in FCC Cell = " << numParticlesInCell << endl;
		log << "Solid Density = " << densitySolid << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Minimisation FAILED = " << endl;
	}
	
	return 1;
}










////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


void minOverNpoints( int argc, char** argv, double kT, double mu, 
                     double &aLattice_min, double &density_min, 
                     double &numParticles_min, double &freeEnergy_min,
                     bool &success)
{
	////////////////////////////// Parameters //////////////////////////////
	
	int NpointsMin = 1;
	int NpointsMax = 1;
	int NpointsStep = 1;
	int ncopy = 1;
	double dx = 0.1;
	
	Options options;
	options.addOption("NpointsMin", &NpointsMin);
	options.addOption("NpointsMax", &NpointsMax);
	options.addOption("NpointsStep", &NpointsStep);
	options.addOption("ncopy", &ncopy);
	options.addOption("Dx", &dx);
	options.read(argc, argv);
	
	Log log("log_min_Npoints.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	/////////////////////////// DFT Computations ///////////////////////////
	
	
	// initialise variables
	
	success = true;
	
	int NumLatticeSizes = (NpointsMax-NpointsMin)/NpointsStep+1;
	vector<int> Npoints(NumLatticeSizes);
	vector<double> aLattice(NumLatticeSizes);
	
	for(int i=0; i<NumLatticeSizes; i++)
	{
		Npoints[i] = NpointsMin+NpointsStep*i;
		aLattice[i] = Npoints[i]*dx;
	}
	
	vector<bool> successPerNpoint(NumLatticeSizes, true);
	vector<double> density(NumLatticeSizes);
	vector<double> numParticles(NumLatticeSizes);
	vector<double> freeEnergy(NumLatticeSizes);
	
	// perform DFT computation for each Npoints
	
	for(int i=0; i<NumLatticeSizes; i++)
	{
		bool success_final;
		double alpha_final;
		double density_final;
		double numParticles_final;
		double freeEnergy_final;
		
#ifdef RESTRICT_GAUSSIAN_PROFILE
		minOverAlpha( argc, argv, kT, mu, Npoints[i], alpha_final, density_final, 
                      numParticles_final, freeEnergy_final, success_final);
#endif

#ifndef RESTRICT_GAUSSIAN_PROFILE
		DFTcomputation( argc, argv, kT, mu, Npoints[i], -1, false, density_final,
		                numParticles_final, freeEnergy_final, success_final);
#endif

		successPerNpoint[i] = success_final;
		density[i] = density_final;
		numParticles[i] = numParticles_final;
		freeEnergy[i] = freeEnergy_final;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Npoints = ";
	for (int i=0;i<NumLatticeSizes; i++) log << Npoints[i] << " ";
	log << endl;
	
	log << "freeEnergy = ";
	for (int i=0;i<NumLatticeSizes; i++) log << freeEnergy[i] << " ";
	log << endl;
	
	log << "density = ";
	for (int i=0;i<NumLatticeSizes; i++) log << density[i] << " ";
	log << endl;
	
	log << "numParticles = ";
	for (int i=0;i<NumLatticeSizes; i++) log << numParticles[i] << " ";
	log << endl;
	
	log << "successPerNpoint = ";
	for (int i=0;i<NumLatticeSizes; i++) log << successPerNpoint[i] << " ";
	log << endl;
	
	
	/////////////////////// Plot Free Energy Profile ///////////////////////
	
	/*
	string outFile = "free_energy_profile_Npoints";
	//string outFile =  "free_energy_profile_Npoints_mu="+mu+"_kT="+kT
	
	Grace graphFreeEnergy(800,600,1,false);
	graphFreeEnergy.addDataSet(aLattice, freeEnergy);
	//graphFreeEnergy.setTitle("Free Energy profile for mu="+mu+" and kT="+kT);
	graphFreeEnergy.setTitle("Free Energy profile");
	graphFreeEnergy.setXAxisLabel("Lattice constant");
	graphFreeEnergy.setYAxisLabel("Free Energy");
	graphFreeEnergy.redraw();
	graphFreeEnergy.printToFile(outFile,4);
	*/
	
	
	//////////////////////////// Parabolic Fit /////////////////////////////
	
	// find the data point that minimise the free energy
	
	int index_min = 0;
	for (int i=0; i<NumLatticeSizes; i++)
		if (freeEnergy[i]<freeEnergy[index_min] && successPerNpoint[i]) 
			index_min = i;
	
	// select a few points in the vicinity of the min data point
	
	vector<int> indicesToFit;
	for (int i=-2; i<3; i++)
	{
		if (index_min+i>=0 && index_min+i<NumLatticeSizes && successPerNpoint[index_min+i])
			indicesToFit.push_back(index_min+i);
		else 
		{
			success = false;
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: The set of lattice sizes is inappropriate, the minimum is either out of the set or too close to its boundaries" << endl;
			return;
		}
	}
	int numIndicesToFit = indicesToFit.size();
	
	// find the Npoints that minimises the free energy and the minimum value
	
	arma::vec x(numIndicesToFit);
	arma::vec y(numIndicesToFit);
	
	for (int i=0; i<numIndicesToFit; i++) x(i) = Npoints[indicesToFit[i]];
	for (int i=0; i<numIndicesToFit; i++) y(i) = freeEnergy[indicesToFit[i]];
	
	arma::vec polyFreeEnergy = arma::polyfit(x,y,2);
	
	double Npoints_min = -1.0/2*polyFreeEnergy(1)/polyFreeEnergy(0);
	aLattice_min = Npoints_min*dx/ncopy;
	freeEnergy_min = polyFreeEnergy(0)*Npoints_min*Npoints_min +
	                 polyFreeEnergy(1)*Npoints_min + polyFreeEnergy(2);
	
	// find the corresponding values of the density and number of particles
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = density[indicesToFit[i]];
	arma::vec polyDensity = arma::polyfit(x,y,2);
	
	density_min = polyDensity(0)*Npoints_min*Npoints_min +
	              polyDensity(1)*Npoints_min + polyDensity(2);
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = numParticles[indicesToFit[i]];
	arma::vec polyNumParticles = arma::polyfit(x,y,2);
	
	numParticles_min = polyNumParticles(0)*Npoints_min*Npoints_min +
	                   polyNumParticles(1)*Npoints_min + polyNumParticles(2);
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Npoints) Min Npoints = " << Npoints_min << endl;
	log << "(Npoints) Min Free Energy = " << freeEnergy_min << endl;
	log << "(Npoints) Min Density = " << density_min << endl;
	log << "(Npoints) Min Number of Particles = " << numParticles_min << endl;
}










////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


void minOverAlpha( int argc, char** argv, double kT, double mu, int Npoints,
                   double &alpha_min, double &density_min, 
                   double &numParticles_min, double &freeEnergy_min,
                   bool &success)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double alphaMin = 100;
	double alphaMax = 100;
	double alphaStep = 10;
	
	Options options;
	options.addOption("GaussianAlphaMin", &alphaMin);
	options.addOption("GaussianAlphaMax", &alphaMax);
	options.addOption("GaussianAlphaStep", &alphaStep);
	options.read(argc, argv);
	
	Log log(("log_min_alpha_Npoints="+to_string(Npoints)+".dat").c_str());
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	/////////////////////////// DFT Computations ///////////////////////////
	
	
	// initialise variables
	
	success = true;
	
	//int NumAlphas = (alphaMax-alphaMin)/alphaStep+1;
	int NumAlphas = std::log(alphaMax/alphaMin)/alphaStep+1;
	vector<double> alpha(NumAlphas);
	
	for(int i=0; i<NumAlphas; i++)
		//alpha[i] = alphaMin+alphaStep*i;
		alpha[i] = alphaMin * exp(alphaStep*i);
	
	vector<bool> successPerAlpha(NumAlphas, true);
	vector<double> density(NumAlphas);
	vector<double> numParticles(NumAlphas);
	vector<double> freeEnergy(NumAlphas);
	
	// perform DFT computation for each Npoints
	
	for(int i=0; i<NumAlphas; i++)
	{
		bool success_final;
		double density_final;
		double numParticles_final;
		double freeEnergy_final;
		
		DFTcomputation( argc, argv, kT, mu, Npoints, alpha[i], false, density_final,
		                numParticles_final, freeEnergy_final, success_final);
		
		successPerAlpha[i] = success_final;
		density[i] = density_final;
		numParticles[i] = numParticles_final;
		freeEnergy[i] = freeEnergy_final;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "alpha = ";
	for (int i=0;i<NumAlphas; i++) log << alpha[i] << " ";
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
	
	
	/////////////////////// Plot Free Energy Profile ///////////////////////
	
	/*
	string outFile = "free_energy_profile_alpha";
	//string outFile =  "free_energy_profile_alpha_Npoints="+Npoints+"_mu="+mu+"_kT="+kT
	
	Grace graphFreeEnergy(800,600,1,false);
	graphFreeEnergy.addDataSet(alpha, freeEnergy);
	//graphFreeEnergy.setTitle("Free Energy profile for Npoints="+Npoints+", mu="+mu+" and kT="+kT);
	graphFreeEnergy.setTitle("Free Energy profile");
	graphFreeEnergy.setXAxisLabel("Gaussian alpha parameter");
	graphFreeEnergy.setYAxisLabel("Free Energy");
	graphFreeEnergy.redraw();
	graphFreeEnergy.printToFile(outFile,4);
	*/
	
	
	//////////////////////////// Parabolic Fit /////////////////////////////
	
	// find the data point that minimise the free energy
	
	int index_min = 0;
	for (int i=0; i<NumAlphas; i++)
		if (freeEnergy[i]<freeEnergy[index_min] && successPerAlpha[i]) 
			index_min = i;
	
	// select a few points in the vicinity of the min data point
	
	vector<int> indicesToFit;
	for (int i=-2; i<3; i++)
	{
		if (index_min+i>=0 && index_min+i<NumAlphas && successPerAlpha[index_min+i])
			indicesToFit.push_back(index_min+i);
		else 
		{
			success = false;
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: The set of gaussian alpha parameters is inappropriate, the minimum is either out of the set or too close to its boundaries" << endl;
			return;
		}
	}
	int numIndicesToFit = indicesToFit.size();
	
	// find the Npoints that minimises the free energy and the minimum value
	
	arma::vec x(numIndicesToFit);
	arma::vec y(numIndicesToFit);
	
	for (int i=0; i<numIndicesToFit; i++) x(i) = alpha[indicesToFit[i]];
	for (int i=0; i<numIndicesToFit; i++) y(i) = freeEnergy[indicesToFit[i]];
	
	arma::vec polyFreeEnergy = arma::polyfit(x,y,2);
	
	alpha_min = -1.0/2*polyFreeEnergy(1)/polyFreeEnergy(0);
	freeEnergy_min = polyFreeEnergy(0)*alpha_min*alpha_min +
	                 polyFreeEnergy(1)*alpha_min + polyFreeEnergy(2);
	
	// find the corresponding values of the density and number of particles
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = density[indicesToFit[i]];
	arma::vec polyDensity = arma::polyfit(x,y,2);
	
	density_min = polyDensity(0)*alpha_min*alpha_min +
	              polyDensity(1)*alpha_min + polyDensity(2);
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = numParticles[indicesToFit[i]];
	arma::vec polyNumParticles = arma::polyfit(x,y,2);
	
	numParticles_min = polyNumParticles(0)*alpha_min*alpha_min +
	                   polyNumParticles(1)*alpha_min + polyNumParticles(2);
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Gaussian alpha) Min Npoints = " << alpha_min << endl;
	log << "(Gaussian alpha) Min Free Energy = " << freeEnergy_min << endl;
	log << "(Gaussian alpha) Min Density = " << density_min << endl;
	log << "(Gaussian alpha) Min Number of Particles = " << numParticles_min << endl;
}











////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// This function is just an intermediate to catch EtaTooLarge exceptions

void DFTcomputation( int argc, char** argv, double kT, double mu, 
                     int Npoints, double alpha, bool runMinimizer,
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, bool &success            )
{
	try
	{
		DFTcomputation2( argc, argv, kT, mu, Npoints, alpha, runMinimizer,
		                 density_final, numParticles_final, freeEnergy_final,
		                 success);
	}
	catch( Eta_Too_Large_Exception &e)
	{
		success = false;
	}
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


  /**
   *   @brief  Performs DFT computation at fixed FCC cell size
   */

void DFTcomputation2( int argc, char** argv, double kT, double mu, 
                      int Npoints, double alpha, bool runMinimizer,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, bool &success            )
{
  //////////////////////////////  Parameters ///////////////////////////////
  
  // control
  int nCores = 6;
  bool showGraphics = true;
  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  
  // geometry
  double dx = 0.1;
  double L[3] = {10,10,10};

  // thermodynamics -> as input

  // potential
  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;
  double hsd1 = -1;
  
  // density initialisation
  int ncopy = 1;
  double prefac = 1;
  double Nvac = 0;
  
  // minimisation
  double forceLimit = 1e-4;
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
  options.addOption("InputFile", &infile);

  // geometry
  options.addOption("Dx", &dx);
  
  // potential
  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);
  
  // density initialisation
  options.addOption("prefac", &prefac);
  options.addOption("ncopy", &ncopy);
  options.addOption("Nvac", &Nvac);
  if (alpha<0) options.addOption("GaussianAlpha", &alpha);

  // minimisation
  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);
  options.addOption("AlphaFac", &alphaFac);
  options.addOption("MaxSteps", &maxSteps);
  
  options.read(argc, argv);

  Log log("log_DFT.dat");
  log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

  log << myColor::RED << myColor::BOLD << "Input parameters (arguments):" << myColor::RESET << endl <<  "#" << endl;

  log << "kT = " << kT << endl;
  log << "mu = " << mu << endl;
  log << "Npoints = " << Npoints << endl;
  log << "alpha = " << alpha << endl;
  
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

  ////////////////////////  Initialise the System //////////////////////////
  
  
  //////////////////////////////////////
  ////// Create potential && effective hsd

  LJ potential1(sigma1, eps1, rcut1);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
  
  /////////////////////////////////////
  // Create density objects


  SolidFCC theDensity1(dx, L, hsd1);

  //theDensity1.initialize2(alpha, alatt, ncopy, prefac);
  theDensity1.initialize(alpha, alatt, ncopy, prefac);

  
  /////////////////////////////////////
  // Create the species objects
  
  FMT_Species species1(theDensity1,hsd1,pointsFile);

  Interaction i1(species1,species1,potential1,kT,log);

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
  // Fix the chemical potential of the surfactant species.

  double N = theDensity1.getNumberAtoms();
  //    species1.setFixedMass(N);
  //    log << "N fixed at " << N << endl;
  //    species1.setChemPotential(0);

  species1.setChemPotential(mu);
  
  
  /////////////////  Thermodynamics of the Liquid state ////////////////////
  
  vector<double> densities;
  densities.push_back(N/(L[0]*L[1]*L[2]));

  double mu1 = dft.Mu(densities,0);
  
  log << "Omega/(V kT) = " << dft.Omega(densities) << endl;  
  log << "mu liq = " << mu1 << endl;
  log << "F liq = " << dft.Fhelmholtz(densities)*L[0]*L[1]*L[2] << endl;
  

  //  check(theDensity1, dft,i1);


  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
  
  /////////////////////////////  Minimisation //////////////////////////////
  
  fireMinimizer2 minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  if (runMinimizer) minimizer.run(maxSteps);

  /////////////////////////////  Final Results /////////////////////////////

  Natoms = theDensity1.getNumberAtoms();
  double Omega;
  if (runMinimizer)
    Omega = minimizer.getF();
  else
    Omega = dft.calculateFreeEnergyAndDerivatives(false);
  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log << "dft.Omega = " << dft.Omega(densities) << endl;
  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Omega: " << Omega << endl;
  log << "Excess Omega = " << dOmega << endl;
  
  density_final = Natoms/(L[0]*L[1]*L[2]);
  numParticles_final = Natoms;
  freeEnergy_final = Omega;
  
  // Sanity check
  
  if (minimizer.getCalls() >= maxSteps-1 && maxSteps > 0)
    // TODO: check if F is stable with err>1
    success = false;
}









////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


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
