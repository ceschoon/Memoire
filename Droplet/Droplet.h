#ifndef __LUTSKO_DROPLET__
#define __LUTSKO_DROPLET__

#include "Density.h"
#include "Display.h" 
 
/**
  *  @brief Class that specializes Density to a slit-pore geometry with a fixed spherical particle
  */  

class Droplet : public Density
{
 public:
  /**
   *   @brief  Default  constructor for a droplet in a periodic box 
   *  
   *   @param  dx is lattice spacing: assumed to be the same in all directions
   *   @param  L[] are the dimensions of the physical region
   *   @param  sigV is the length scale of the wall potential
   *   @param  epsV is the energy scale of the wall potential
   *   @return nothing 
   */  
 Droplet(double dx, double L[], Grace *g, double hsd)
   : Density(dx,L),  g_(g), hsd_(hsd){}
  
  /**
   *   @brief  Generates an intial guess at the density
   *  
   *   @param  density_inside is the interior density
   *   @param  density_outside is exterior density
   *   @param  R is the radius
   *   @return  
   */  
  virtual void initialize(double density_inside, double density_outside, double R)
  {
    Density::initialize(density_inside,density_outside);

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k);
	    
	    double r2 = x*x+y*y+z*z;

	    double dd = (r2 < R*R ? density_inside : density_outside);

	    set_Density_Elem(i,j,k,dd);
	  }
    // This is an object that writes a png snapshot whenever doDisplay gets called.
    display_ = new Display(Nx_,Ny_);
  }

  // This gets called after every update and for each species: seq is the species number. 
  
  virtual void doDisplay(string &title, string &file, int seq) const
  {
    // Write to a grace window
    if(g_ == NULL) return;

    g_->deleteDataSet(seq);

    for(int i= 0; i < Nz(); i++)
      {
	double x = getZ(i);
	g_->addPoint(x,getDensity(Nx()/2,Ny()/2,i),seq);
      }
    g_->setTitle(title.c_str());
    g_->redraw();

    // make a png snapshot
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	display_->setData(i+Nx_*j, min(getDensity(i,j,Nz_/2),1.0));

    stringstream title1; title1 << "Species " << seq << " " << title;
    stringstream file1;  file1 << seq << "_" << file;

    string t1 = title1.str();
    string f1 = file1.str();
    
    display_->doDisplay(t1, f1);
  }  

 protected:
  Grace *g_;
  int sequence_;
  double hsd_;
  Display *display_;
};



#endif // SLIT_PORE2__
