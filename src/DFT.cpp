#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>

using namespace std;


#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"

//static double xmax = 1.4;
//static int    Nmax = 200; 


double DFT::F_IdealGas(Density& density, DFT_Vec& dF)
{
  double Fideal = 0.0;
  double dV = density.dV();
  
  bool calcForces = (dF.size() >= density.Ntot());
  
  for(long i=0;i<density.Ntot();i++)
    {
      double d0 = density.getDensity(i);
      if(d0 > -SMALL_VALUE)
	{
	  Fideal += (d0*log(d0)-d0)*dV;
	  if(calcForces) dF.addTo(i,log(d0)*dV);
	} else {
	if(calcForces) dF.addTo(i,log(SMALL_VALUE)*dV);
      }
    }
  return Fideal;
}

double DFT::F_External(Density& density, double mu, DFT_Vec& dF)
{
  double dV = density.dV();
  double Fx = density.getExternalFieldEnergy()*dV - density.getNumberAtoms()*mu;
  
  if(dF.size() >= density.Ntot())
    dF.Increment_Shift_And_Scale(density.getExternalField(),dV,mu);

  return Fx;
}


template <class T> DFT_FMT<T>::DFT_FMT(Lattice &lattice, FMT_Species &species)
  : fmt_(lattice) { fmt_.addSpecies(species);}

template <class T>
double DFT_FMT<T>::calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex)
{
  double F = 0;
  dF.zeros(density.Ntot());

  try{
    F = fmt_.calculateFreeEnergyAndDerivatives(density,dF);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  cout << "F fmt = " << F << endl;
  
  if(!onlyFex) // add in ideal gas and external field and chem potential
    {  
      F += F_IdealGas(density, dF);       // Ideal gas 
      F += F_External(density, mu, dF);   // External field
      cout << "F ideal = " << F_IdealGas(density, dF)  << endl;
      cout << "F ext = " << F_External(density, mu, dF) << endl;
    }
  return F;
}


template <class T>
double DFT_FMT<T>::Xliq_From_Mu(double mu) const 
{ 
  double xa = 0.5;
  double xb = 0.5;

  while(Enskog(xa).chemPotentialCS() > mu) { xa /= 2;}

  double xmax = 6/M_PI;

  while(Enskog(xb).chemPotentialCS() < mu) { xb = (xmax+xb)/2; }

  while(fabs(xa-xb) > 1e-8)
    {
      double x = (xa+xb)/2;
      if(Enskog(x).chemPotentialCS() > mu) xb = x;
      else xa = x;
    }
  return (xa+xb)/2;
}

#include "Potential1.h"

template <class T>
DFT_VDW<T>::DFT_VDW(Lattice& lattice,   Potential1& potential,  string& pointsFile, double kT)
  : DFT(),  vdw_(0,0), species_(NULL)
{
  // The lattice
  long Nx = lattice.Nx();
  long Ny = lattice.Ny();
  long Nz = lattice.Nz();

  double dx = lattice.getDX();
  double dy = lattice.getDY();
  double dz = lattice.getDZ();
  
  // Set up the mean field potential

  w_att_.initialize(Nx,Ny,Nz);

  for(int nx = 0;nx<Nx;nx++)
    for(int ny = 0;ny<Ny;ny++)
      for(int nz = 0;nz<Nz;nz++)
	{
	  long pos = nz+Nz*(ny+Ny*nx);
	  
	  double x = nx*dx;
	  double y = ny*dy;
	  double z = nz*dz;

	  if(nx > Nx/2) x -= Nx*dx;
	  if(ny > Ny/2) y -= Ny*dy;
	  if(nz > Nz/2) z -= Nz*dz;

	  double r2 = x*x+y*y+z*z;
	  w_att_.Real().addTo(pos,potential.Watt(sqrt(r2))/kT); 
 	}
  
  w_att_.do_real_2_fourier();

  // Set the parameters in the VDW object  
  double hsd   = potential.getHSD(kT);
  double a_vdw = 0.5*w_att_.Real().accu()*dx*dy*dz;

  vdw_.set_VDW_Parameter(a_vdw);
  vdw_.set_HardSphere_Diameter(hsd);

  species_ = new FMT_Species(hsd,lattice,pointsFile);
  
  // Create hard-sphere object
  dft_fmt_ = new DFT_FMT<T>(lattice, *species_); 

  // initialize working space for calculations ...
  v_mean_field_.initialize(Nx,Ny,Nz);

}

template <class T>
double DFT_VDW<T>::calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex)
{
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  long Ntot = density.Ntot();
  long Nout = density.Nout();

  double dV = density.dV();

  dF.zeros(Ntot);
  double F = 0;

  // Hard sphere contributions to F and dF
  try {
    F = dft_fmt_->calculateFreeEnergyAndDerivatives(density,mu,dF,onlyFex);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  // Mean field contribution to F and dF
  // Divide by Ntot because of normalization of fft

  v_mean_field_.Four().Schur(density.getDK(),w_att_.Four());
  v_mean_field_.Four().multBy(dV*dV/Ntot);
  v_mean_field_.do_fourier_2_real(); 

  F  += 0.5*density.getInteractionEnergy(v_mean_field_.Real());
  dF.IncrementBy(v_mean_field_.Real());

  return F;
}



// These lines are needed to force instantiation of the relevant classes. 


template class DFT_FMT<WhiteBearI>;
template class DFT_FMT<WhiteBearII>;
template class DFT_FMT<RSLT>;
template class DFT_FMT<RSLT2>;

template class DFT_VDW<WhiteBearI>;
template class DFT_VDW<WhiteBearII>;
template class DFT_VDW<RSLT>;
template class DFT_VDW<RSLT2>;
