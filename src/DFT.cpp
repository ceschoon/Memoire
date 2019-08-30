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


double DFT::Mu(const vector<double> &x, int species) const
{
  double mu = log(x[species]);

  if(fmt_)
    mu += fmt_->BulkMuex(x, allSpecies_, species);
  
  for(auto &interaction: Interactions_)
    mu += interaction->Mu(x,species);
  return mu;
}

double DFT::Omega(const vector<double> &x) const
{
  double omega = Fhelmholtz(x);
  for(int i=0;i<allSpecies_.size();i++)
    if(fabs(x[i]) > SMALL_VALUE) omega -= x[i]*Mu(x,i);
  return omega;
}

double DFT::Fhelmholtz(const vector<double> &x) const
{
  double F = 0.0;
  for(auto &y: x)
    if(fabs(y) > SMALL_VALUE)
      F += y*log(y)-y;

  if(fmt_)
    F += fmt_->BulkFex(x, allSpecies_);

  for(auto &interaction: Interactions_)
    F += interaction->Fhelmholtz(x);  

  return F;  
}  


double DFT::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  for(auto &species : allSpecies_)  
    species->zeroForce();
  
  for(auto &s: allSpecies_) s->beginForceCalculation();

  double F = calculateFreeEnergyAndDerivatives_internal_(onlyFex);

  for(auto &s: allSpecies_)
    s->endForceCalculation();

  return F;
}


double DFT::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  double F = 0.0;
  
  if(onlyFex) return F; // do nothing.

  for(auto &species : allSpecies_)
    {
      const Density& density = species->getDensity();
      double dV = density.dV();
      
      for(long i=0;i<density.Ntot();i++)
	{
	  double d0 = density.getDensity(i);
	  if(d0 > -SMALL_VALUE)
	    {
	      F += (d0*log(d0)-d0)*dV;
	      species->addToForce(i,log(d0)*dV);
	    } else {
	    species->addToForce(i,log(SMALL_VALUE)*dV);
	  }
	}
    }

  for(auto &species : allSpecies_)
    F += species->externalField(true); // bCalcForces = true: obsolete?  

  if(fmt_)
    {    
      try{
	F += fmt_->calculateFreeEnergyAndDerivatives(allSpecies_);
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    }
  
  return F;
  
}
/*
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
*/



#include "Potential1.h"
/*
template <class T>
double DFT_VDW<T>::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  double F = 0;

  // Hard sphere contributions to F and dF
  try {
    F = DFT_FMT<T>::calculateFreeEnergyAndDerivatives_internal_(onlyFex);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  // Mean field contribution to F and dF
  for(auto &s: DFT::allSpecies_)
      F += ((VDW_Species *) s)->getInteractionEnergyAndForces();

  // Mean field contribution to F and dF
  for(auto &interaction: DFT::Interactions_)
    F += interaction->getInteractionEnergyAndForces();  

  return F;
}

*/
