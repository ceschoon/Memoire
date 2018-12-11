#ifndef __LUTSKO__SPECIES__
#define __LUTSKO__SPECIES__

#include "FMT_Weighted_Density.h"
#include "Lattice.h"
#include "Density.h"
/**
  *  @brief Species Class: hard sphere diameters, etc.
  *
  *    This class holds an 11-dimentional array called d_ . This in turn holds the weighted densities as
  *             d_ = {eta(N),s(N),V1(N),...,VD(N), T11(N), T12(N), ... T1D(N), T22(N), ..., TDD(N)}
  *             where for D=3 there are 1+1+3+(3+2+1) = 11 entries.
  *             Note that each of these, e.g. eta(N), is an N-dimensional vector holding the value of the weighted density at each point.
  */

class FMT_Species
{
 public:
  /**
  *   @brief  Default  constructor for FMT_Species 
  *  
  *   @param  hsd is the hard-sphere diameter
  *   @param  lattice describes the mesh
  *   @param  pointsFile contains the points for spherical integration
  *   @return nothing 
  */    
 FMT_Species(double hsd, string &pointsFile): hsd_(hsd), d_(11), bInitialized_(false), pointsFile_(pointsFile){}

  void Initialize(Lattice &lattice);

  FMT_Species(const FMT_Species &) = delete;
  
  ~FMT_Species(){}

/**
  *   @brief  Accessor for hard sphere diameter
  *  
  *   @param  none
  *   @return hsd_
  */      
  double getHSD() const { return hsd_;}

/**
  *   @brief  get value of Eta at pos
  *  
  *   @param  pos is the mesh position
  *   @return value of Eta at pos
  */    
  double getEta(long pos)           const { return d_[EI()].r(pos);}

/**
  *   @brief  get value of S at pos
  *  
  *   @param  pos is the mesh position
  *   @return value of S at pos
  */      
  double getS(long pos)             const { return d_[SI()].r(pos);}

/**
  *   @brief  get value of component j of V at pos
  *  
  *   @param  j is index of V
  *   @param  pos is the mesh position
  *   @return value of V(j) at pos
  */      
  double getV(int j, long pos)      const { return d_[VI(j)].r(pos);}

/**
  *   @brief  get value of component j,k of T at pos
  *  
  *   @param  j is first index of T
  *   @param  k is second index of T
  *   @param  pos is the mesh position
  *   @return value of T(j,k) at pos
  */      
  double getT(int j,int k,long pos) const { return d_[TI(j,k)].r(pos);}
  
  const DFT_Vec_Complex& getWEK() const { return d_[0].wk();}

/**
  *   @brief  set value of dPhi_dEta at pos
  *  
  *   @param  pos is the mesh position
  *   @return none
  */        
  void Set_dPhi_Eta(long i,double val)             {d_[EI()].Set_dPhi(i,val);}

/**
  *   @brief  set value of dPhi_dS at pos
  *  
  *   @param  pos is the mesh position
  *   @return none
  */        
  void Set_dPhi_S(long i,double val)               {d_[SI()].Set_dPhi(i,val);}

/**
  *   @brief  set value of dPhi_dV_j at pos
  *  
  *   @param  j is index of V
  *   @param  pos is the mesh position
  *   @return none
  */        
  void Set_dPhi_V(int j, long i,double val)        {d_[VI(j)].Set_dPhi(i,val);}

/**
  *   @brief  set value of dPhi_dT(j,k) at pos
  *  
  *   @param  j is first index of T
  *   @param  k is second index of T
  *   @param  pos is the mesh position
  *   @return none
  */        
  void Set_dPhi_T(int j, int k, long i,double val) {d_[TI(j,k)].Set_dPhi(i,val);}

  void convolute(Density &density)
  {
    if(!bInitialized_) Initialize(density);
    
    // reference to Fourier-space array of density
    const DFT_Vec_Complex &rho_k = density.getDK();

    // This does the convolution of the density and the weight for each weighted density after which it converts back to real space 
    // ( so this computes the weighted densities n(r) = int w(r-r')rho(r')dr'). The results are all stored in parts of FMT_Weighted_Density
    for(FMT_Weighted_Density &d: d_)
      d.convolute(rho_k);
  }

  void Accumulate_dPhi(DFT_Vec_Complex& dPhi)
  {
    for(FMT_Weighted_Density &d: d_)
      d.add_to_dPhi(dPhi);
  }


  const DFT_Vec &getV_Real(int J) const { return d_[2+J].Real();}
  const DFT_Vec_Complex& getVweight_Four(int J) const { return d_[2+J].wk();}  
  
protected:
 /**
  *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
  *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). Their evaluation involves real-space integrations for which the 
  *           integration points are given in the file pointsFile. Most of the work occurs via a call to initializeWeightedDensities@
  *
  *   @param  densities: the array of weighted densities
  *   @param  hsd is the hard-sphere diameter
  *   @param  pointsFile is the file holding the integration points for integrating a spherical shell
  *   @param  Nx is the number of lattice points in the x-direction
  *   @param  Ny is the number of lattice points in the y-direction
  *   @param  Nz is the number of lattice points in the z-direction
  */        
 void generateWeights(Lattice &lattice);

 int EI() const {return 0;}
 int SI() const {return 1;}
 int VI(int j) const {return 2+j;}
 int TI(int j, int k) const
 {
   if(j > k) swap(j,k);
   if(j == 0) return 5+k;
   else if (j == 1) return 7+k;
   return 10;
 }

 protected:
  double hsd_ = 0.0; ///< hard sphere diameter 
 
  vector<FMT_Weighted_Density>  d_; ///< all weighted densities in real & fourier space

  bool bInitialized_ = false;
  string pointsFile_;
};


#endif // __LUTSKO__SPECIES__
