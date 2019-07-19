#ifndef __LUTSKO_DFT__
#define __LUTSKO_DFT__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>




#include "Enskog.h"

#include "Potential1.h"
#include "VDW1.h"
#include "Density.h"
#include "FMT.h"


/*! \mainpage classicalDFT: Finite Temperature Density Functional Theory in 3 dimensions
 *
 * \section intro_sec Introduction
 *
 * This code is a library that can be used to build practical applications of ftDFT. It supports the most common Fundamental Measure Theory DFT models for hard spheres (including White Bear) and can be used for crystallization. At present, it only handles single-component systems: generalization to multicomponent systems is in the works. 
 *
 * \section install_sec Installation
 *
 *  The library just needs to be deposited into a directory. Applications should be built in sub-directories of this main directory.
 *
 * etc...
 */


/**
  *  @brief DFT base Class. This class enapsulates the DFT models (e.g. Fundamental measure Theory - FMT; Mean Field Theories; etc.)
  *
  *  @detailed This is the base class for all of the DFT objects.
  */  

class DFT
{
 public:
  /**
   *   @brief  Default  constructor for DFT 
   *  
   *   @return nothing 
   */  
  DFT(){}
  /**
   *   @brief  Default  destructor for DFT
   *  
   *   @return nothing 
   */  
  ~DFT(){}


  /**
   *   @brief  Tells the DFT that there is another species
   *  
   *   @param  s is the Species object
   */ 
  void addSpecies(Species* s) {allSpecies_.push_back(s);}

  int getNumberOfSpecies() const {return allSpecies_.size();}

  double getNumberAtoms(int i) const {return allSpecies_[i]->getDensity().getNumberAtoms();}

  const Lattice& lattice() const {return allSpecies_.front()->getLattice();}

  double get_convergence_monitor() const { double d = 0; for(auto& s: allSpecies_) d += s->get_convergence_monitor(); return d;}
  
  const Density& getDensity(int i) const {return allSpecies_[i]->getDensity();}

  DFT_Vec &getDF(int i) {return allSpecies_[i]->getDF();}
  
  void doDisplay(string &title, string &file)
  {
    for(auto &x: allSpecies_) x->doDisplay(title,file);
  }

  /**
   *   @brief  Sets the densities based on amplitudes by passing to the species object
   *  
   *   @param  i is the species
   *   @param x is the amplitde
   */  
  void set_density_from_amplitude(int i,DFT_Vec &x) {allSpecies_[i]->set_density_from_amplitude(x);}

  void writeDensity(int i, string &of) const {allSpecies_[i]->getDensity().writeDensity(of);}

  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions. 
   *            This is just a wrapper that calls calculateFreeEnergyAndDerivatives_internal_ which does the real work. 
   *  
   *   @param  onlyFex : if true, the ideal and external contributions are not added: IS THIS PARAMETER NEEDED???
   *   @return total free energy of system
   */  
  double calculateFreeEnergyAndDerivatives(bool onlyFex);

  /**
   *   @brief  Compute chemical potential/kT for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @param  species is the species for which we calculate the chemical potential
   *   @return mu/kT
   */   
  virtual double Mu(const vector<double> &x, int species) const {return log(x[species]);}

  /**
   *   @Brief  Compute grand potenial/kT/V for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @return Omega/(kT * V)
   */   
  virtual double Omega(const vector<double> &x) const {double omega = Fhelmholtz(x); for(int i=0;i<allSpecies_.size();i++) if(fabs(x[i]) > SMALL_VALUE) omega -= x[i]*Mu(x,i); return omega;}

  /**
   *   @brief  Compute Helmholtz free energy/kT/V for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @return F/(kT * V)
   */   
  virtual double Fhelmholtz(const vector<double> &x) const {double F = 0.0; for(auto &y: x) if(fabs(y) > SMALL_VALUE) F += y*log(y)-y; return F;}  

  /**
   *   @brief  Returns the name of the DFT object - i.e. identifies the model being used.
   *   @return the name as a string object.
   */     
  virtual string Name() const = 0;

  /**
   *   @brief  This conveys information to the FMT object for models that use it.
   */     
  virtual void setEtaMax(double etaMax) {}

  /**
   *   @brief  Only implemented for DFT_FMT object.
   */     
  virtual double Xliq_From_Mu(double mu) const {throw std::runtime_error("Not implemented");}

 protected:
  
  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions. Here, this is just the ideal gas contribution.
   *  
   *   @param  onlyFex : if true, the ideal and external contributions are not added
   *   @return total free energy of system
   */  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex);


 protected:
  vector<Species*> allSpecies_;

  
};


/**
  *  @brief DFT_FMT Class: A hard-sphere object.
  *
  *   @detailed A single-species hard-sphere object.
  *  
  */  

template <class T> class DFT_FMT : public DFT
{
 public:
  /**
  *   @brief  Default  constructor for DFT_FMT
  *
  *   @return nothing 
  */  
  
  // External field
  DFT_FMT(int Nx, int Ny, int Nz)
    : fmt_(Nx,Ny,Nz) {}


  /**
  *   @brief  Default  destructor for DFT_FMT
  *  
  *   @return nothing 
  */  
  ~DFT_FMT(){}


/**
  *   @brief  Compute chemical potential/kT for given density
  *  
  *   @param  x is the array of densities
  *   @param  species is the species for which we calculate the chemcical potential
  *   @return mu/kT
  */   
  virtual double Mu(const vector<double> &x, int species) const {return DFT::Mu(x,species) + fmt_.BulkMuex(x, allSpecies_, species);}

/**
  *   @brief  Compute Helmhltz Free Energy/kT/V for given density
  *  
  *   @param  x is the array of densities
  *   @return F/(kT * V)
  */   
  virtual double Fhelmholtz(const vector<double> &x) const {return  DFT::Fhelmholtz(x) + fmt_.BulkFex(x, allSpecies_);}


/**
  *   @brief  Find liquid with given chem potential/kT
  *  
  *   @param  mu is chem potential
  *   @return mu/kT
  */   
  virtual double Xliq_From_Mu(double mu) const;

  
/**
  *   @brief  The name of this DFT object
  *  
  *   @return name as string
  */   
  virtual string Name() const { return string("DFT_FMT : ") + fmt_.Name();}

 /**
  *   @brief  This conveys information to the FMT object for models that use it.
  */     
 virtual void setEtaMax(double etaMax) {fmt_.setEtaMax(etaMax);}

 protected:
 /**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
 virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex = false);
 
 protected:
  T         fmt_;   ///< Hard-sphere FMT object
};



/**
  *  @brief FMT plus van der Waals mean field tail.
  *
  *   @detailed A single-species mean-field DFT model. It holds a DFT_FMT object to do the hard-sphere stuff.
  *  
  */  

template <class T> class DFT_VDW : public DFT_FMT<T>
{
 public:
  /**
   *   @brief  Default  constructor for DFT 
   *  
   *   @param  Nx is number of points in the x-direction
   *   @param  Ny is number of points in the y-direction
   *   @param  Nz is number of points in the z-direction
   */  
  DFT_VDW(int Nx, int Ny, int Nz);

  /**
   *   @brief  Default  destructor for DFT 
   *  
   *   @return nothing 
   */  
  ~DFT_VDW(){} 

  /**
   *   @brief  Compute chemical potential/kT for given density: note that interspecies interaction fields are not supported.
   *  
   *   @param  x is the  array of densities
   *   @param  species is the species for which we calculate the chemical potential
   *   @return mu/kT
   */   
  virtual double Mu(const vector<double> &x, int species) const;

  /**
   *   @brief  Compute helmholtz free energy/kT/V for given density: note that interspecies interaction fields are not supported.
   *  
   *   @param  x is the  array of densities
   *   @return F/(kT * V)
   */     
  virtual double Fhelmholtz(const vector<double> &x) const;

  
  virtual string Name() const { return string("DFT_VDW : ") + DFT_FMT<T>::Name();} 

 protected:
  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
   *  
   *   @param  density is current density
   *   @param  mu is the chemical potential
   *   @return total free energy of system
   */  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex = false);    
  
 protected:
  DFT_FFT v_mean_field_;

};

/**
  *  @brief DFT_VDW_Surfactant class
  *
  *   @detailed This extends DFT_VDW with functionality to model a surfactant. The additional contributions depend only on gradients
  *              of the other species, so the bulk thermodynamics is just the same as that inherited from DFT_VDW. 
  *  
  */  

template <class T> class DFT_VDW_Surfactant : public DFT_VDW<T>
{
 public:
  /**
  *   @brief  Default  constructor for DFT 
  *  
  *   @param  Nx is number of points in x direction
  *   @param  Ny is number of points in y direction
  *   @param  Nz is number of points in z direction
  *   @return nothing 
  */  
  DFT_VDW_Surfactant(int Nx, int Ny, int Nz);

  /**
  *   @brief  Default  destructor
  */  
  ~DFT_VDW_Surfactant(){}
  
  void addSurfactantPotential(Potential1& pot, double kT);

  virtual string Name() const { return string("DFT_VDW_Surfactant : ") +  DFT_VDW<T>::Name();} 

 protected:
  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
   *  
   *   @param  density is current density
   *   @param  mu is the chemical potential
   *   @return total free energy of system
   */  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex = false);
  
 protected:
  DFT_FFT surfactant_potential_; //< Arrays holding surfactant assymetric potential
};

#endif // __LUTSKO_DFT__
