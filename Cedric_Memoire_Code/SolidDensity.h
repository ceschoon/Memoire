////////////////////////////////////////////////////////////////////////////
//                                                                        //
//      Header file for methods related to the generation of solid        //
//      density profiles.                                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef SOLID_DENSITY_H
#define SOLID_DENSITY_H

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


class SolidDensity : public Density
{
	public:
	
	SolidDensity(double dx, double L[], double hsd) : Density(dx,L),  hsd_(hsd)	
	{
		#ifdef USE_GRACE
			grace_ = new Grace();
		#endif
		#ifdef USE_MGL
			// This is an object that writes a png snapshot whenever doDisplay gets called.
			display_ = new Display(Nx_,Ny_);
		#endif
	}

	~SolidDensity()
	{
		#ifdef USE_GRACE
			if(grace_ != NULL) delete grace_;
		#endif
		#ifdef USE_MGL
			if(display_ != NULL) delete display_;
		#endif
	}
	
	
	virtual void initialiseUniform(double rho)
	{
		for(int i=0;i<Nx_;i++)
			for(int j=0;j<Ny_;j++)
				for(int k=0;k<Nz_; k++)
					set_Density_Elem(i,j,k,rho);
	}
	
	
	virtual void initialiseFCCGaussianBasis(double alpha, double a_latt, int ncells, double prefac)
	{
		// basis for one unit cell (conventional cubic cell)
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
		
		// loop over grid points
		for(int i=0;i<Nx_;i++)
		for(int j=0;j<Ny_;j++)
		for(int k=0;k<Nz_;k++)
		{
			double x = getX(i);
			double y = getY(j);
			double z = getZ(k);
			
			double dsum = 0;
			
			// sum the contributions from all lattice sites
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
			
			// Uncomment these lines to help avoid eta>1 exceptions
			double small_value = 1e-10;
			if(dsum < small_value) dsum = small_value;
			
			set_Density_Elem(i,j,k,dsum);
		}
	}
	
	
	virtual void initialiseFCCGaussianFromSC(double alpha, double a_latt, int ncells, double prefac)
	{
		// loop over grid points
		for(int i=0;i<Nx_;i++)
		for(int j=0;j<Ny_;j++)
		for(int k=0;k<Nz_;k++)
		{
			double x = getX(i);
			double y = getY(j);
			double z = getZ(k);
			
			double dsum = 0;
			
			// sum the contributions from all lattice sites
			// we use a simple cubic with half the FCC lattice constant
			// and keep the indices whose sum is a multiple of two
			for(int icell=0;icell < 2*ncells; icell++) 
			for(int jcell=0;jcell < 2*ncells; jcell++)
			for(int kcell=0;kcell < 2*ncells; kcell++)
			{
				// condition that the SC lattice site must fullfill in
				// order to be part of the FCC lattice
				if ((icell+jcell+kcell)%2!=0) continue;
				
				double dx = fabs(x-icell*a_latt/2); if(dx > L_[0]/2) dx -= L_[0];
				double dy = fabs(y-jcell*a_latt/2); if(dy > L_[1]/2) dy -= L_[1];
				double dz = fabs(z-kcell*a_latt/2); if(dz > L_[2]/2) dz -= L_[2];
				
				double r2 = dx*dx+dy*dy+dz*dz;
				dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
			}
			
			// Uncomment these lines to help avoiding eta>1 exceptions
			double small_value = 1e-10;
			if(dsum < small_value) dsum = small_value;
			
			set_Density_Elem(i,j,k,dsum);
		}
	}
	
	// WRONG LATTICE: THIS IS NOT HCP (THE ARTICLE MUST BE INCORRECT)
	// For details on the generation of the HCP lattice, see
	// DOI: 10.1063/1.1997138
	
	virtual void initialiseHCPGaussianFromSC(double alpha, double a_latt, int ncells, double prefac)
	{
		// loop over grid points
		for(int i=0;i<Nx_;i++)
		for(int j=0;j<Ny_;j++)
		for(int k=0;k<Nz_;k++)
		{
			double x = getX(i);
			double y = getY(j);
			double z = getZ(k);
			
			double dsum = 0;
			
			// sum the contributions from all lattice sites
			// we use a simple cubic with half the FCC lattice constant
			// and keep the indices whose sum is a multiple of two
			for(int icell=0;icell < 6*ncells; icell++) 
			for(int jcell=0;jcell < 6*ncells; jcell++)
			for(int kcell=0;kcell < 6*ncells; kcell++)
			{
				// conditions that the SC lattice site must fullfill in
				// order to be part of the HCP lattice
				if ((icell+jcell+kcell)%6!=0) continue;
				if (abs(icell-jcell)%3!=0) continue;
				if (abs(icell-kcell)%3!=0) continue;
				if (abs(jcell-kcell)%3!=0) continue;
				if (((icell+jcell+kcell)/6 - icell)%3!=0) continue;
				if (((icell+jcell+kcell)/6 - jcell)%3!=0) continue;
				if (((icell+jcell+kcell)/6 - kcell)%3!=0) continue;
				
				if (i==0 && j==0 && k==0)
				{
					cout << "icell = " << icell << " "
					     << "jcell = " << jcell << " "
					     << "kcell = " << kcell << " "
					     << endl;
				}
				
				double dx = fabs(x-icell*a_latt/6); if(dx > L_[0]/2) dx -= L_[0];
				double dy = fabs(y-jcell*a_latt/6); if(dy > L_[1]/2) dy -= L_[1];
				double dz = fabs(z-kcell*a_latt/6); if(dz > L_[2]/2) dz -= L_[2];
				
				double r2 = dx*dx+dy*dy+dz*dz;
				dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
			}
			
			// Uncomment these lines to help avoid eta>1 exceptions
			double small_value = 1e-10;
			if(dsum < small_value) dsum = small_value;
			
			set_Density_Elem(i,j,k,dsum);
		}
	}
	
	
	
	
	
	
	virtual void initialiseHexagonal(double alpha, double a_latt, 
				double b_latt, double c_latt, int ncells, double prefac)
	{
		cout << "Initialisation of Solid Density (Hexagonal Lattice)" << endl;
		
		// basis for one unit cell (custom rectangular cell)
		double atoms[4][3];
		
		atoms[0][0] = 0.0;
		atoms[0][1] = 0.0;
		atoms[0][2] = 0.0;
		
		atoms[1][0] = a_latt/2;
		atoms[1][1] = b_latt/2;
		atoms[1][2] = 0.0;
		
		atoms[2][0] = a_latt/2;
		atoms[2][1] = b_latt/6;
		atoms[2][2] = c_latt/2;
		
		atoms[3][0] = 0.0;
		atoms[3][1] = 4*b_latt/6;
		atoms[3][2] = c_latt/2;
		
		// loop over grid points
		for(int i=0;i<Nx_;i++)
		for(int j=0;j<Ny_;j++)
		for(int k=0;k<Nz_;k++)
		{
			double x = getX(i);
			double y = getY(j);
			double z = getZ(k);
			
			double dsum = 0;
			
			// sum the contributions from all lattice sites
			for(int icell=0;icell < ncells; icell++) 
			for(int jcell=0;jcell < ncells; jcell++)
			for(int kcell=0;kcell < ncells; kcell++)
			for(int l=0;l<4;l++)
			{
				double dx = fabs(x-atoms[l][0]-icell*a_latt); if(dx > L_[0]/2) dx -= L_[0];
				double dy = fabs(y-atoms[l][1]-jcell*b_latt); if(dy > L_[1]/2) dy -= L_[1];
				double dz = fabs(z-atoms[l][2]-kcell*c_latt); if(dz > L_[2]/2) dz -= L_[2];
				
				double r2 = dx*dx+dy*dy+dz*dz;
				dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
			}
			
			// Uncomment these lines to help avoid eta>1 exceptions
			double small_value = 1e-10;
			if(dsum < small_value) dsum = small_value;
			
			set_Density_Elem(i,j,k,dsum);
		}
		
		cout << "Finished Initialisation of Solid Density" << endl;
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
		{
			if(getDensity(ix,iy,iz) > dmax)
			{
				dmax = getDensity(ix,iy,iz);
				jx = ix;
				jy = iy;
				jz = iz;
			}
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
	
	#ifdef USE_GRACE
		Grace *grace_ = NULL;
	#endif
	#ifdef USE_MGL
		Display *display_ = NULL;
	#endif
	
};

#endif // SOLID_DENSITY_H
