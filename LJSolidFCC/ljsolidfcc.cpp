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

void check(Density &density, DFT &dft, Interaction &i);
void DFTcomputation( int argc, char** argv, int Npoints,
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final                          );


// TODO
// Subtility: Do not know where the minima relative to Npoints is.
//            and there can be more than one minimum.
// => Plot the free energy landscape (gnuplot) (grace? if in file)
// => Fit only the points near the global minimum


int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	int NpointsMin = 1;
	int NpointsMax = 1;
	int NpointsStep = 1;
	
	Options options;
	options.addOption("NpointsMin", &NpointsMin);
	options.addOption("NpointsMax", &NpointsMax);
	options.addOption("NpointsStep", &NpointsStep);
	options.read(argc, argv);
	
	Log log("log_min_cell_size.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	/////////////////////////// DFT Computations ///////////////////////////
	
	
	// initialise variables
	
	int NumLatticeSizes = (NpointsMax-NpointsMin)/NpointsStep+1;
	int Npoints[NumLatticeSizes];
	
	for(int i=0; i<NumLatticeSizes; i++)
		Npoints[i] = NpointsMin+NpointsStep*i;
	
	double density[NumLatticeSizes];
	double numParticles[NumLatticeSizes];
	double freeEnergy[NumLatticeSizes];
	
	// perform DFT computation for each Npoints
	
	for(int i=0; i<NumLatticeSizes; i++)
	{
		double density_final;
		double numParticles_final;
		double freeEnergy_final;
		
		DFTcomputation( argc, argv, Npoints[i], density_final, 
		                numParticles_final, freeEnergy_final   );
		
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
	
	
	//////////////////////////// Parabolic Fit /////////////////////////////
	
	
	// find the Npoints that minimises the free energy and the minimum value
	
	arma::vec x(NumLatticeSizes);
	arma::vec y(NumLatticeSizes);
	
	for (int i=0; i<NumLatticeSizes; i++) x(i) = Npoints[i];
	for (int i=0; i<NumLatticeSizes; i++) y(i) = freeEnergy[i];
	
	arma::vec polyFreeEnergy = arma::polyfit(x,y,2);
	
	double Npoints_min = -1.0/2*polyFreeEnergy(1)/polyFreeEnergy(0);
	double freeEnergy_min = polyFreeEnergy(0)*Npoints_min*Npoints_min +
	                        polyFreeEnergy(1)*Npoints_min + polyFreeEnergy(2);
	
	// find the corresponding values of the density and number of particles
	
	for (int i=0; i<NumLatticeSizes; i++) y(i) = density[i];
	arma::vec polyDensity = arma::polyfit(x,y,2);
	
	double density_min = polyDensity(0)*Npoints_min*Npoints_min +
	                     polyDensity(1)*Npoints_min + polyDensity(2);
	
	for (int i=0; i<NumLatticeSizes; i++) y(i) = numParticles[i];
	arma::vec polyNumParticles = arma::polyfit(x,y,2);
	
	double numParticles_min = polyNumParticles(0)*Npoints_min*Npoints_min +
	                          polyNumParticles(1)*Npoints_min + polyNumParticles(2);
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Min Npoints = " << Npoints_min << endl;
	log << "Min Free Energy = " << freeEnergy_min << endl;
	log << "Min Density = " << density_min << endl;
	log << "Min Number of Particles = " << numParticles_min << endl;
	
	return 1;
}



  /**
   *   @brief  Performs DFT computation at fixed FCC cell size
   */

void DFTcomputation( int argc, char** argv, int Npoints, 
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final                          )
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

  // thermodynamics
  double kT = 1;
  double Mu = -1.0;

  // potential
  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;
  double hsd1 = -1;
  
  // density initialisation
  int ncopy = 1;
  double prefac = 1;
  double Nvac = 0;
  double alpha = 50;
  double alphaFac = 1;
  
  // minimisation
  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;
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
  
  // thermodynamics
  options.addOption("kT", &kT);
  options.addOption("Mu", &Mu);
  
  // potential
  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);
  
  // density initialisation
  options.addOption("prefac", &prefac);
  options.addOption("ncopy", &ncopy);
  options.addOption("Nvac", &Nvac);
  options.addOption("alpha", &alpha);

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

  log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

  options.write(log);
  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;


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
  
#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  log << "omp max threads = " << omp_get_max_threads() << endl;
#endif



  ////////////////////////  Initialise the System //////////////////////////
  
  
  //////////////////////////////////////
  ////// Create potential && effective hsd

  LJ potential1(sigma1, eps1, rcut1);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
  
  /////////////////////////////////////
  // Create density objects


  SolidFCC theDensity1(dx, L, hsd1);

  theDensity1.initialize2(alpha, alatt, ncopy, prefac, Mu/kT);

  
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

  species1.setChemPotential(Mu);
  
  
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
  minimizer.run(maxSteps); 

  /////////////////////////////  Final Results /////////////////////////////

  Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log << "dft.Omega = " << dft.Omega(densities) << endl;
  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Omega: " << Omega << endl;
  log << "Excess Omega = " << dOmega << endl;
  
  density_final = Natoms/(L[0]*L[1]*L[2]);
  numParticles_final = Natoms;
  freeEnergy_final = Omega;
}



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
