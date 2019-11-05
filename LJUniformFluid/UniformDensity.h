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
  *  @brief Class that specializes Density to a uniform fluid
  */  

class UniformDensity : public Density
{
 public:
 
 UniformDensity(double dx, double L[], double hsd)
   : Density(dx,L),  hsd_(hsd)
   {
#ifdef USE_GRACE
      grace_ = new Grace();
#endif
   }

   ~UniformDensity()
   {
#ifdef USE_GRACE
     if(grace_  != NULL) delete grace_;
#endif
   }







  virtual void initialize(double rho)
  {
    cout << "Nx = " << Nx_ << " Ny = " << Ny_ << " Nz = " << Nz_ << endl;
    cout << "Lx = " << L_[0] << " Ly = " << L_[1] << " Lz = " << L_[2] << endl;
    cout << " rho = " << rho << endl;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        for(int k=0;k<Nz_; k++)
          {
            set_Density_Elem(i,j,k,rho);
          }
    
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
