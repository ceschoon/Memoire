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

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "SolidPhase.h"


// Things to be done

// TODO: More significant digits in number of particles
//       (print vacancy concentration instead)


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


void minOverNpoints( int argc, char** argv, Log &log,  double kT, double mu, 
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
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Minimisation over Npoints ---" << myColor::RESET << endl <<  "#" << endl;
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
		double density_final;
		double numParticles_final;
		double freeEnergy_final;

		DFTcomputation( argc, argv, log, kT, mu, Npoints[i], true, density_final,
		                numParticles_final, freeEnergy_final, success_final);

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

// This function is just an intermediate to catch EtaTooLarge exceptions

void DFTcomputation( int argc, char** argv, Log &log, double kT, double mu, 
                     int Npoints, bool runMinimizer,
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, bool &success            )
{
	try
	{
		DFTcomputation2( argc, argv, log, kT, mu, Npoints, runMinimizer,
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


// Remark: the free energy is given per unit volume

void DFTcomputation2( int argc, char** argv, Log &log, double kT, double mu, 
                      int Npoints, bool runMinimizer,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, bool &success            )
{
	////////////////////////////  Parameters ///////////////////////////////

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
	double alpha = 100;
	double alphaMin = 50;
	double alphaMax = 5000;
	double alphaStep = 0.3;

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
	options.addOption("GaussianAlpha", &alpha);
	options.addOption("GaussianAlphaMin", &alphaMin);
	options.addOption("GaussianAlphaMax", &alphaMax);
	options.addOption("GaussianAlphaLogStep", &alphaStep);

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





	///////////////////////  Initialise the System /////////////////////////


	//////////////////////////////////////
	////// Create potential && effective hsd

	LJ potential1(sigma1, eps1, rcut1);

	if(hsd1 < 0) hsd1 = potential1.getHSD(kT);

	/////////////////////////////////////
	// Create density objects


	SolidFCC theDensity1(dx, L, hsd1);

	//theDensity1.initialize2(alpha, alatt, ncopy, prefac);
	//theDensity1.initialize(alpha, alatt, ncopy, prefac);


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

	species1.setChemPotential(mu);

	//  check(theDensity1, dft,i1);

	log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;





	//////////////////// Pre-run with gaussian profiles ////////////////////

#ifndef NO_GAUSSIAN_PRE_RUN

	// We first do multiple runs with different widths of the gaussian 
	// in order to select the best profile as the initial condition of the 
	// full computation. 
	
	// initialise variables
	
	int NumAlphas = std::log(alphaMax/alphaMin)/alphaStep+1;
	vector<double> alphaVec(NumAlphas);
	
	for(int i=0; i<NumAlphas; i++)
		alphaVec[i] = alphaMin * exp(alphaStep*i);
	
	vector<bool> successPerAlpha(NumAlphas, true); // useless in the current code
	vector<double> density(NumAlphas);
	vector<double> numParticles(NumAlphas);
	vector<double> freeEnergy(NumAlphas);
	
	// DFT computation for each alpha
	
	for(int i=0; i<NumAlphas; i++)
	{
		bool success_temp = true; 
		double density_temp;
		double numParticles_temp;
		double freeEnergy_temp; // per unit volume 
		
		theDensity1.initialize(alphaVec[i], alatt, ncopy, prefac);
		
		numParticles_temp = theDensity1.getNumberAtoms();
		density_temp = numParticles_temp/(L[0]*L[1]*L[2]);
		freeEnergy_temp = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
		
		successPerAlpha[i] = success_temp;
		density[i] = density_temp;
		numParticles[i] = numParticles_temp;
		freeEnergy[i] = freeEnergy_temp;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "--- Pre-run with gaussian profiles ---" << endl <<  "#" << endl;
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






	//////////////// Find the best alpha for the gaussians /////////////////
	
	// We find the best alpha by minimizing the free energy 
	// We do the minimisation of a parabolic fit near the min data point
	
	// find the data point that minimises the free energy
	
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
	
	for (int i=0; i<numIndicesToFit; i++) x(i) = alphaVec[indicesToFit[i]];
	for (int i=0; i<numIndicesToFit; i++) y(i) = freeEnergy[indicesToFit[i]];
	
	arma::vec polyFreeEnergy = arma::polyfit(x,y,2);
	
	double alpha_min = -1.0/2*polyFreeEnergy(1)/polyFreeEnergy(0);
	double freeEnergy_min = polyFreeEnergy(0)*alpha_min*alpha_min +
	                        polyFreeEnergy(1)*alpha_min + polyFreeEnergy(2);
	
	// find the corresponding values of the density and number of particles
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = density[indicesToFit[i]];
	arma::vec polyDensity = arma::polyfit(x,y,2);
	
	double density_min = polyDensity(0)*alpha_min*alpha_min +
	                     polyDensity(1)*alpha_min + polyDensity(2);
	
	for (int i=0; i<numIndicesToFit; i++) y(i) = numParticles[indicesToFit[i]];
	arma::vec polyNumParticles = arma::polyfit(x,y,2);
	
	double numParticles_min = polyNumParticles(0)*alpha_min*alpha_min +
	                          polyNumParticles(1)*alpha_min + polyNumParticles(2);
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Gaussian alpha) Min alpha = " << alpha_min << endl;
	log << "(Gaussian alpha) Min Free Energy = " << freeEnergy_min << endl;
	log << "(Gaussian alpha) Min Density = " << density_min << endl;
	log << "(Gaussian alpha) Min Number of Particles = " << numParticles_min << endl;

	// choose the minimum as the starting point for the full (unrestricted) minimisation 

	alpha = alpha_min;

#endif




	//////////////  Full Minimisation with the best profile ////////////////

#ifndef RESTRICT_GAUSSIAN_PROFILE

	theDensity1.initialize(alpha, alatt, ncopy, prefac);
	
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

#ifndef NO_GAUSSIAN_PRE_RUN
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Final Omega/V: " << freeEnergy_min << endl;
	
	density_final = density_min;
	numParticles_final = numParticles_min;
	freeEnergy_final = freeEnergy_min;
	
#endif

#ifndef RESTRICT_GAUSSIAN_PROFILE
	
	Natoms = theDensity1.getNumberAtoms();
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

