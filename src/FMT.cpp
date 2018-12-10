#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "FMT.h"
#include "Enskog.h"


ostringstream Eta_Too_Large_Exception::cnvt;

double sigma2 = 1e-4;


FMT::FMT(Lattice &lattice)
  : dx_(lattice.getDX()), dy_(lattice.getDY()), dz_(lattice.getDZ()), etaMax_(1e30), dPhi_(lattice.Nx(),lattice.Ny(),lattice.Nz()){}

double FMT::dPHI(long i)
{
  // the weighted densities for the lattice site under consideration (i).
  // The basic densities are the scalars eta and s, the vector v and the tensor T.
  // Note that s0,s1 etc only differ in constant prefactors if there is only one species present.
  double eta = 0.0;
  double s0  = 0.0;
  double s1  = 0.0;
  double s2  = 0.0;

  double v1[3] = {0.0,0.0,0.0};
  double v2[3] = {0.0,0.0,0.0};
  double T0[3][3] = {{0.0,0.0,0.0},
		 {0.0,0.0,0.0},
		 {0.0,0.0,0.0}};

  // Collect the contributions to the various weighted densities at lattice position i
  // hsd1_ has to do with mixtures (this has not been fully implementd. hsd1_ > 0 is a flag for whether or not there is a second species).
  for(FMT_Species &species : AllSpecies_)
    {
        double hsd = species.getHSD();
  
	eta += species.getEta(i); //Eta(dd).r(i);

	s0 += species.getS(i)/(hsd*hsd) ;
	s1 += species.getS(i)/hsd;
	s2 += species.getS(i);

	for(int j=0;j<3;j++)
	  {
	    v1[j] += species.getV(j,i)/hsd;
	    v2[j] += species.getV(j,i);
	    for(int k=0;k<3;k++)
	      T0[j][k] += species.getT(j,k,i);
	  }
    }

  // touch up the normalization of the tensor density to account for any numerical errors
  double ss = T0[0][0]+T0[1][1]+T0[2][2];

  T0[0][0] +=  (s2-ss)/3;
  T0[1][1] +=  (s2-ss)/3;
  T0[2][2] +=  (s2-ss)/3;

  // Calculate some useful quantities: v dot T, T dot T etc
  double vT[3] = { T0[0][0]*v2[0] +T0[0][1]*v2[1] +T0[0][2]*v2[2],
		   T0[1][0]*v2[0] +T0[1][1]*v2[1] +T0[1][2]*v2[2],
		   T0[2][0]*v2[0] +T0[2][1]*v2[1] +T0[2][2]*v2[2]};
  double TT[3][3] =  { {T0[0][0]*T0[0][0] +T0[0][1]*T0[1][0] +T0[0][2]*T0[2][0],
		     T0[1][0]*T0[0][0] +T0[1][1]*T0[1][0] +T0[1][2]*T0[2][0],
		     T0[2][0]*T0[0][0] +T0[2][1]*T0[1][0] +T0[2][2]*T0[2][0]},

		    {T0[0][0]*T0[0][1] +T0[0][1]*T0[1][1] +T0[0][2]*T0[2][1],
		     T0[1][0]*T0[0][1] +T0[1][1]*T0[1][1] +T0[1][2]*T0[2][1],
		     T0[2][0]*T0[0][1] +T0[2][1]*T0[1][1] +T0[2][2]*T0[2][1]},

		    {T0[0][0]*T0[0][2] +T0[0][1]*T0[1][2] +T0[0][2]*T0[2][2],
		     T0[1][0]*T0[0][2] +T0[1][1]*T0[1][2] +T0[1][2]*T0[2][2],
		     T0[2][0]*T0[0][2] +T0[2][1]*T0[1][2] +T0[2][2]*T0[2][2]}};

  double v1_v2  = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; //dot(v1,v2);
  double v2_v2  = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]; //dot(v1,v2);
  double vTv    = v2[0]*vT[0]+v2[1]*vT[1]+v2[2]*vT[2]; //dot(v1,v2);
  double T2     = TT[0][0]+TT[1][1]+TT[2][2]; //trace(TT);
  double T3     = TT[0][0]*T0[0][0]+TT[0][1]*T0[1][0]+TT[0][2]*T0[2][0]
    +TT[1][0]*T0[0][1]+TT[1][1]*T0[1][1]+TT[1][2]*T0[2][1]
    +TT[2][0]*T0[0][2]+TT[2][1]*T0[1][2]+TT[2][2]*T0[2][2];

  // These are the eta-dependendent cofactors that lie at the heart of FMT
  double f1 = f1_(eta);
  double f2 = f2_(eta);
  double f3 = f3_(eta);

  // Now, construct the local value of phi.
  double phi = 0;
  phi -= (1/M_PI)*s0*f1; //log(1-eta);
  phi += (1/(2*M_PI))*(s1*s2-v1_v2)*f2;
  phi += Phi3(s2,v2_v2,vTv,T2,T3)*f3; 

  if(etaMax_ > 1.0)
    if(eta > 0.5 && 1-eta < 0.0)
      throw Eta_Too_Large_Exception();

  // Also add in the contributions to the derivative of phi (at lattice site i) wrt the various weighted densities
  // (part of the chain-rule evaluation of dPhi/drho(j) = dPhi/deta(i) * deta(i)/drho(j) + ...)
  for(FMT_Species &species : AllSpecies_)
    {
      double hsd = species.getHSD();
      
        double f1 = f1_(eta);
  
	double f1p = f1p_(eta);
	double f2p = f2p_(eta);
	double f3p = f3p_(eta);

	double dPhi_dEta = 0;
	dPhi_dEta -= (1/M_PI)*s0*f1p; ///(1-eta);
	dPhi_dEta += (1/(2*M_PI))*(s1*s2-v1_v2)*f2p;
	dPhi_dEta += Phi3(s2,v2_v2,vTv,T2,T3)*f3p; 

	double dPhi_dS2 = 0;
	dPhi_dS2 += -(1/M_PI)*f1/(hsd*hsd);
	dPhi_dS2 += (1/(2*M_PI))*((s2/hsd)+s1)*f2;
	dPhi_dS2 += dPhi3_dS2(s2,v2_v2,vTv,T2,T3)*f3; //(1.0/(8*M_PI))*(-v2_v2+T2)*f3; 

	double dPhi_dV0[3];
	
	for(int k=0;k<3;k++)
	  {
	    dPhi_dV0[k] = 0.0;
	    dPhi_dV0[k] += (1/(2*M_PI))*(-v1[k]-v2[k]/hsd)*f2;
	    dPhi_dV0[k] += dPhi3_dV2(k, s2, v2_v2, v2, vT)*f3; //(1.0/(8*M_PI))*(2*vT[k]-2*s2*v2[k])*f3;
	  }
	
	double dPhi_dT[3][3];

	for(int k=0;k<3;k++)
	  for(int j=0;j<3;j++)
	    dPhi_dT[j][k] = dPhi3_dT(j,k,s2,v2,T0, TT)*f3; // (1.0/(8*M_PI))*(v2[j]*v2[k]-3*TT(j,k)+2*s2*T0(j,k))*f3; 
		
	species.Set_dPhi_Eta(i,dPhi_dEta);
	species.Set_dPhi_S(i,dPhi_dS2);

	for(int j=0;j<3;j++)
	  {
	    species.Set_dPhi_V(j,i,dPhi_dV0[j]);
	    for(int k=j;k<3;k++)	   
	      species.Set_dPhi_T(j,k,i,(j ==k ? 1 : 2)*dPhi_dT[j][k]); // taking account that we only use half the entries
	  }
    }

  return phi;
}


// MIXTURES NOT IMPLEMENTED
double FMT::calculateFreeEnergy(Density &density)
{


  /*
  // Compute the FFT of density 
  density.doFFT();

  // reference to Fourier-space array of density
  const DFT_Vec_Complex &rho_k = density.getDK();

  // This does the convolution of the density and the weight for each weighted density after which it converts back to real space 
  // ( so this computes the weighted densities n(r) = int w(r-r')rho(r')dr'). The results are all stored in parts of FMT_Weighted_Density
  for(FMT_Weighted_Density &d: d0_)
    d.convolute(rho_k);

  */

  cout << "Density max in real space: " << density.getDensity().max() << endl;
  AllSpecies_.front().convolute(density);

  long   Ntot = density.Ntot();  
  
  double emax = 0.0;  
  for(long i=0;i<Ntot;i++)
    emax = max(emax, AllSpecies_.front().getEta(i));

  cout << "Eta = " << emax << endl;
  


  

  double F = doFreeEnergyLoop(Ntot);

  return F;
}

// MIXTURES NOT IMPLEMENTED
double FMT::doFreeEnergyLoop(long Ntot)
{
  // Now compute the free energy. Here we loop over all lattice sites and compute Phi(r_i) for each one. This presupposes that we did the convolution above. 
  
  double F = 0;
  int chunk = Ntot/20;
  long i;

  // There is a problem throwing exceptions from an OMP loop - I think that only the affected thread stops and the others continue.
  // So we eat the exception, remember it, and rethrow it when all loops are finished.
  bool hadCatch = false;

  #pragma omp parallel for		\
    shared( chunk    )				\
    private(i)					\
    schedule(static,chunk)			\
    reduction(+:F)
  for(i=0;i<Ntot;i++)
    {
      try {
	F += dPHI(i);
      } catch( Eta_Too_Large_Exception &e) {
	hadCatch = true;
      }
    }
  // rethrow exception if it occurred: this messiness is do to the parallel evaluation. 
  if(hadCatch) 
    throw   Eta_Too_Large_Exception();

  cout << "From FMT: F = " << F << " dPHI(0) = " << dPHI(0) << endl;
  return F;
}

// MIXTURES NOT IMPLEMENTED
//HERE: this cannot be right for a mixture ... TBD

// Calculate dF[i] = dPhi/drho(i)
//                 = sum_j dV * dPhi(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  dn_{alpha}(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  w_{alpha}(j,i)
// This is done with convolutions: FMT_Weighted_Density is an array with the index (alpha)
// and holds the weights, w_{alpha}(j,i) and their FFT's AND dPhi(j)/dn_{alpha}(j).
// It therefore FFT's both of these and adds them to dPhi_.Four.
// Once this is done of all alpha, dPhi_ FFT's back to real space and the result is put into dF (with a factor of dV thrown in).

void FMT::calculateFreeEnergyDerivatives(double dV, DFT_Vec& dF)
{
  dPhi_.Four().zeros();

  AllSpecies_.front().AccumulatedPhi(dPhi_.Four());
  
  //  for(FMT_Weighted_Density &d: dd)
  //    d.add_to_dPhi(dPhi_.Four());

  dPhi_.do_fourier_2_real();

  dF.set(dPhi_.cReal());
  dF.multBy(dV);
}


double FMT::calculateFreeEnergyAndDerivatives(Density& density, DFT_Vec& dF0) 
{
  double F = 0;

  try {
    F = calculateFreeEnergy(density);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  double dV   = density.dV();
  calculateFreeEnergyDerivatives(dV, dF0);
  return F*dV;
};

