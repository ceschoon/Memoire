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

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

void check(Density &density, DFT &dft, Interaction &i);

int main(int argc, char** argv)
{
  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 6;

  double zPos = 0;

  double gaussian_std = 0.05;
  double n_atoms_in_cell = 3.9999;
  double eta2 = 0.2;

  double hsd1 = -1;
  double hsd2 = -1;
  
  double epsWall = 1;
  double sigWall = 1;

  double kT = 1;

  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;

  double eps2   = 1;
  double sigma2 = 1;
  double rcut2  = 3;

  double eps_surf   = 1;
  double sigma_surf = 1;
  double rcut_surf  = 3;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  
  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;

  bool showGraphics = true;

  int maxSteps = -1;

  double Asurf = 0;
  double rho_surf = -1;

  double Mu = -1.0;


  int ncopy = 1;
  double prefac = 1;
  double Nvac = 0;
  int Npoints = 1;
  double alpha = 50;
  double alphaFac = 1;
  double dx = 0.1;
  
  Options options;

  options.addOption("prefac", &prefac);
  options.addOption("ncopy", &ncopy);
  options.addOption("Nval", &Nvac);
  options.addOption("Npoints", &Npoints);
  options.addOption("alpha", &alpha);
  options.addOption("Dx", &dx);
  
  options.addOption("nCores", &nCores);
  options.addOption("GaussianStandardDeviation", &gaussian_std);
  options.addOption("NAtomsInCell", &n_atoms_in_cell);

  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("kT", &kT);
  options.addOption("Mu",   &Mu);

  
  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);

  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);
  options.addOption("AlphaFac", &alphaFac);
  options.addOption("MaxSteps", &maxSteps);

  options.addOption("InputFile", &infile);
  
  options.addOption("ShowGraphics", &showGraphics);
  
  options.read(argc, argv);

  Log log("log.dat");
  log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

  log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

  options.write(log);
  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;

  //  double dx = 1.0/PointsPerHardSphere;

  //  if(L[0] < 0) L[0] = dx;
  //  if(L[1] < 0) L[1] = dx;
  //  if(L[2] < 0) L[2] = dx;  


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
  log << "Num Atoms In Unit Cell = " << n_atoms_in_cell << endl;
  log << "Gaussian Std = " << gaussian_std << endl;
  log << "Lx = " << L[0] << endl;
  log << "HSD 1 = " << hsd1 << endl;



  if(! infile.empty())
    theDensity1.readDensity(infile.c_str());

  
  ///////////////////////////////////////////////
  // Fix the mass of the surfactant species.

      double N = theDensity1.getNumberAtoms();
  //    species1.setFixedMass(N);
  //    log << "N fixed at " << N << endl;
  //    species1.setChemPotential(0);

  species1.setChemPotential(Mu);
  
  
  /////////////////////////////////////////////////////
  // The thermodynamics
  // TODO ??
  
  vector<double> densities;
  densities.push_back(N/(L[0]*L[1]*L[2]));

  double mu1 = dft.Mu(densities,0);
  
  log << "Omega/(V kT) = " << dft.Omega(densities) << endl;  
  log << "mu liq = " << mu1 << endl;
  log << "F liq = " << dft.Fhelmholtz(densities)*L[0]*L[1]*L[2] << endl;
  

  //  check(theDensity1, dft,i1);


  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
  
  //  fireMinimizer_Mu minimizer(dft,log);
  fireMinimizer2 minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps); // TODO


  Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log << "dft.Omega = " << dft.Omega(densities) << endl;
  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Omega: " << Omega << endl;
  log << "Excess Omega = " << dOmega << endl;

    
  return 1;
}


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
