#ifndef __LUTSKO_DROPLET__
#define __LUTSKO_DROPLET__

#include <iostream>
#include <vector>
#include <cmath>

#include "Density.h"

#ifdef USE_MGL
#include "Display.h" 
#endif

#ifdef USE_GRACE
#include "Grace.h"
#endif

double const PI = 3.14159265358979;

/**
  *  @brief Class that specializes Density to a FCC cell
  */  

class SolidFCC : public Density
{
 public:
  /**
   *   @brief  Default  constructor for a FCC cell in a periodic box 
   *  
   *   @param  dx is lattice spacing: assumed to be the same in all directions
   *   @param  L[] are the dimensions of the physical region (FCC cube side)
   *   @param  hsd is the hard sphere effective diameter
   *   @return nothing 
   */  
 SolidFCC(double dx, double L[], double hsd)
   : Density(dx,L),  hsd_(hsd)
   {
#ifdef USE_GRACE
      grace_ = new Grace();
#endif
   }

   ~SolidFCC()
   {
#ifdef USE_GRACE
     if(grace_  != NULL) delete grace_;
#endif
   }
   
  /**
   *   @brief  Generates an initial guess at the density
   *
   *   @detailed This guess is the sum of multiple gaussian profiles located at each FCC atom position for one cell
   *  
   *   @param  cell dimensions
   *   @param  standard deviation of gaussian profile
   *   @param  number of atoms in the FCC unit cell (taking vacancies into account)
   *   @return  
   */  
/*
  virtual void initialize(double *L, double gaussian_std, double n_atoms_in_cell)
  {    
    // list of atoms intersecting FCC unit cell
    vector<vector<double>> atomsPositions;
    
    atomsPositions.push_back({0,0,0});
    atomsPositions.push_back({1,0,0});
    atomsPositions.push_back({0,1,0});
    atomsPositions.push_back({0,0,1});
    atomsPositions.push_back({0,1,1});
    atomsPositions.push_back({1,0,1});
    atomsPositions.push_back({1,1,0});
    atomsPositions.push_back({1,1,1});
    atomsPositions.push_back({0,0.5,0.5});
    atomsPositions.push_back({0.5,0,0.5});
    atomsPositions.push_back({0.5,0.5,0});
    atomsPositions.push_back({1,0.5,0.5});
    atomsPositions.push_back({0.5,1,0.5});
    atomsPositions.push_back({0.5,0.5,1});
    
    double Ntot = 0;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        for(int k=0;k<Nz_; k++)
          {
            double x = getX(i);
            double y = getY(j);
            double z = getZ(k);
            
            double dd = 0;
            
            for (int l=0; l<atomsPositions.size(); l++)
            {
                double xl = -L[0]/2+L[0]*atomsPositions[l][0];
                double yl = -L[1]/2+L[1]*atomsPositions[l][1];
                double zl = -L[2]/2+L[2]*atomsPositions[l][2];
                
                // add gaussian contribution from atom l
                double r2 = (x-xl)*(x-xl)+(y-yl)*(y-yl)+(z-zl)*(z-zl);
                double alpha = 1.0/(2*gaussian_std*gaussian_std);
                double rho0 = n_atoms_in_cell/4.0;
                dd += rho0 * pow(alpha/PI,3.0/2) * exp(-alpha*r2);
            }
            
            set_Density_Elem(i,j,k,dd);
          }
    
#ifdef USE_MGL
    // This is an object that writes a png snapshot whenever doDisplay gets called.
    display_ = new Display(Nx_,Ny_);
#endif
  }
*/

  virtual void initialize(double alpha, double a_latt, int ncells, double prefac)
  {
    double atoms[4][3];

    atoms[0][0] = 0.0;
    atoms[0][1] = 0.0;
    atoms[0][2] = 0.0;

    atoms[1][0] = a_latt/2;
    atoms[1][1] = a_latt/2;
    atoms[1][2] = 0.0;

    atoms[2][0] = a_latt/2;
    atoms[2][1] = 0.0;
    atoms[2][2] = a_latt/2;

    atoms[3][0] = 0.0;
    atoms[3][1] = a_latt/2;
    atoms[3][2] = a_latt/2;

    cout << "Nx = " << Nx_ << " Ny = " << Ny_ << " Nz = " << Nz_ << endl;
    cout << "Lx = " << L_[0] << " Ly = " << L_[1] << " Lz = " << L_[2] << endl;
    cout << "atoms[3]: " << atoms[3][0] << " " << atoms[3][1] << " " << atoms[3][2] << endl;
    
    cout << "alpha = " << alpha << " prefac = " << prefac << endl;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        for(int k=0;k<Nz_; k++)
          {
            double x = getX(i);
            double y = getY(j);
            double z = getZ(k);

            double dsum = 0;

            for(int icell=0;icell < ncells; icell++)
              for(int jcell=0;jcell < ncells; jcell++)
                for(int kcell=0;kcell < ncells; kcell++)
                  for(int l=0;l<4;l++)         
                    {
                      double dx = fabs(x-atoms[l][0]-icell*a_latt); if(dx > L_[0]/2) dx -= L_[0];
                      double dy = fabs(y-atoms[l][1]-jcell*a_latt); if(dy > L_[1]/2) dy -= L_[1];
                      double dz = fabs(z-atoms[l][2]-kcell*a_latt); if(dz > L_[2]/2) dz -= L_[2];

                      double r2 = dx*dx+dy*dy+dz*dz;
                      
                      // wide gaussian to avoid EtaTooLargeException
                      double mixing = 0.01;
                      double alpha_wide = 1/(2*a_latt/2*a_latt/2);
                      dsum += mixing*prefac*pow(alpha_wide/M_PI,1.5)*exp(-alpha_wide*r2);
                      
                      dsum += (1-mixing)*prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
                    }
            
            if(dsum < 1e-4) dsum = 1e-4;
            set_Density_Elem(i,j,k,dsum);
          }
    cout << "Nx = " << Nx_ << " Ny = " << Ny_ << " Nz = " << Nz_ << endl;
    cout << "Lx = " << L_[0] << " Ly = " << L_[1] << " Lz = " << L_[2] << endl;
    cout << "a = " << a_latt << endl;
    cout << "Density = " << 4/(a_latt*a_latt*a_latt) << endl;
    
#ifdef USE_MGL
    // This is an object that writes a png snapshot whenever doDisplay gets called.
    display_ = new Display(Nx_,Ny_);
#endif
  }
  
  
  
  
  virtual void initialize2(double alpha, double a_latt, int ncells, double prefac)
  {
    double atoms[4][3];

    atoms[0][0] = 0.0;
    atoms[0][1] = 0.0;
    atoms[0][2] = 0.0;

    atoms[1][0] = a_latt/2;
    atoms[1][1] = a_latt/2;
    atoms[1][2] = 0.0;

    atoms[2][0] = a_latt/2;
    atoms[2][1] = 0.0;
    atoms[2][2] = a_latt/2;

    atoms[3][0] = 0.0;
    atoms[3][1] = a_latt/2;
    atoms[3][2] = a_latt/2;

    cout << "Nx = " << Nx_ << " Ny = " << Ny_ << " Nz = " << Nz_ << endl;
    cout << "Lx = " << L_[0] << " Ly = " << L_[1] << " Lz = " << L_[2] << endl;
    cout << "atoms[3]: " << atoms[3][0] << " " << atoms[3][1] << " " << atoms[3][2] << endl;
    
    cout << "alpha = " << alpha << " prefac = " << prefac << endl;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        for(int k=0;k<Nz_; k++)
          {
            double x = getX(i);
            double y = getY(j);
            double z = getZ(k);

            double dsum = 0;

            for(int icell=0;icell < ncells; icell++)
              for(int jcell=0;jcell < ncells; jcell++)
                for(int kcell=0;kcell < ncells; kcell++)
                  for(int l=0;l<4;l++)         
                    {
                      double dx = fabs(x-atoms[l][0]-icell*a_latt); if(dx > L_[0]/2) dx -= L_[0];
                      double dy = fabs(y-atoms[l][1]-jcell*a_latt); if(dy > L_[1]/2) dy -= L_[1];
                      double dz = fabs(z-atoms[l][2]-kcell*a_latt); if(dz > L_[2]/2) dz -= L_[2];

                      double r2 = dx*dx+dy*dy+dz*dz;
                      dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
                    }
            
            if(dsum < 1e-4) dsum = 1e-4;
            set_Density_Elem(i,j,k,dsum);
          }
    cout << "Nx = " << Nx_ << " Ny = " << Ny_ << " Nz = " << Nz_ << endl;
    cout << "Lx = " << L_[0] << " Ly = " << L_[1] << " Lz = " << L_[2] << endl;
    cout << "a = " << a_latt << endl;
    cout << "Density = " << 4/(a_latt*a_latt*a_latt) << endl;
    
#ifdef USE_MGL
    // This is an object that writes a png snapshot whenever doDisplay gets called.
    display_ = new Display(Nx_,Ny_);
#endif
  }

  // This gets called after every update and for each species: seq is the species number. 
  
  virtual void doDisplay(string &title, string &file, int seq) const
  {
    // find max density point
    
    int jx = 0;
    int jy = 0;
    int jz = 0;
    double dmax = 0;

    for(int ix = 0; ix < Nx(); ix++)
      for(int iy = 0; iy < Ny(); iy++)
        for(int iz = 0; iz < Nz(); iz++)
            if(getDensity(ix,iy,iz) > dmax)
              {
                dmax = getDensity(ix,iy,iz);
                jx = ix;
                jy = iy;
                jz = iz;
              }
    
#ifdef USE_GRACE    
    // Write to a grace window
    if(grace_ == NULL) return;

    grace_->deleteDataSet(seq);

    
    for(int i= 0; i <= Nz(); i++)
      {
        int iz = i-jz;
        while(iz > Nz()/2) iz -= Nz();
        while(iz < -Nz()/2) iz += Nz();
    
        double x = getZ(iz);
        grace_->addPoint(x,getDensity(jx,jy,i),seq);
      }
    grace_->setTitle(title.c_str());
    grace_->redraw();
#endif
#ifdef USE_MGL
    // find max density on plane
    
    double max_density = 0;
    
    for(int i=0;i<Nx_;i++)
    {
        for(int j=0;j<Ny_;j++)
        {
            double d = getDensity(i,j,jz);
            if (d>max_density) max_density = d;
        }
    }
    
    // make a png snapshot --> Adjusted to max density
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        display_->setData(i+Nx_*j, min(getDensity(i,j,jz)/max_density,1.0)); // TODO original line doesn't work when 1.0 is replaced by a bigger number in min(,) 

    stringstream title1; title1 << "Species " << seq << " " << title;
    stringstream file1;  file1 << seq << "_" << file;

    string t1 = title1.str();
    string f1 = file1.str();
    
    display_->doDisplay(t1, f1);
#endif       
  }

 protected:
  int sequence_;
  double hsd_;

#ifdef USE_MGL
  Display *display_;
#endif  
  

#ifdef USE_GRACE  
  Grace *grace_ = NULL;
#endif  

};



#endif // __LUTSKO_DROPLET__
