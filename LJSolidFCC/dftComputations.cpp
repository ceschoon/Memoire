#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <stdlib.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "SolidPhase.h"
#include "SolidPhase.cpp"
#include "../LJUniformFluid/UniformFluid.h"
#include "../LJUniformFluid/UniformFluid.cpp"
#include "utilities.cpp"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to perform dft computations for a given range of the parameters
// (kT,mu,Npoints).

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kTMin = 1;
	double kTMax = 1.5;
	double kTStep = 0.1;
	double muMin = -1.5;
	double muMax = -0.5;
	double muStep = 0.2;
	int NpointsMin = 120;
	int NpointsMax = 135;
	int NpointsStep = 1;
	
	Options options;
	options.addOption("kTMin", &kTMin);
	options.addOption("kTMax", &kTMax);
	options.addOption("kTStep", &kTStep);
	options.addOption("MuMin", &muMin);
	options.addOption("MuMax", &muMax);
	options.addOption("MuStep", &muStep);
	options.addOption("NpointsMin", &NpointsMin);
	options.addOption("NpointsMax", &NpointsMax);
	options.addOption("NpointsStep", &NpointsStep);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////////////// Prepare computations /////////////////////////
	
	string dataDir1 = "Data-kt-mu-Npoints";
	string dataDir2 = "Data-kt-mu";
	string plotDir  = "Plots-kt-mu"; 
	
	int sysresult = system(("mkdir " + dataDir1).c_str());
	    sysresult = system(("mkdir " + dataDir2).c_str());
	    sysresult = system(("mkdir " + plotDir ).c_str());
	
	int kTSize = int( (kTMax-kTMin)/kTStep +1 );
	vector<double> kT(kTSize);
	for (int i=0; i<kTSize; i++) kT[i] = kTMin + kTStep*i;
	
	int muSize = int( (muMax-muMin)/muStep +1 );
	vector<double> mu(muSize);
	for (int i=0; i<muSize; i++) mu[i] = muMin + muStep*i;
	
	int NpointsSize = int( (NpointsMax-NpointsMin)/NpointsStep +1 );
	vector<int> Npoints(NpointsSize);
	for (int i=0; i<NpointsSize; i++) Npoints[i] = NpointsMin + NpointsStep*i;
	
	//////////////////////// Run DFT computations //////////////////////////
	
	double aVdW = 0;
	double hsd = 0;
	
	for (int i=0; i<NpointsSize; i++)
	{
		for (int j=0; j<kTSize; j++)
		{
			for (int k=0; k<muSize; k++)
			{
				// dft computation
				
				double densitySolid = 0;
				double numParticlesSolid = 0;
				double freeEnergySolid = 0;
				bool successSolid = false;
				
				DFTcomputation( argc, argv, log, kT[j], mu[k], Npoints[i], 
				                densitySolid, numParticlesSolid, 
				                freeEnergySolid, aVdW, hsd, successSolid);
				
				// save in file
				
				ofstream dataFile(dataDir1+"/results_"+to_string(i)+
				                                   "_"+to_string(j)+
				                                   "_"+to_string(k)+".dat");
				
				dataFile << "Npoints = " << Npoints[i] << endl;
				dataFile << "kT = " << kT[j] << endl;
				dataFile << "mu = " << mu[k] << endl;
				dataFile << endl;
				dataFile << "densitySolid = " << densitySolid << endl;
				dataFile << "numParticlesSolid = " << numParticlesSolid << endl;
				dataFile << "freeEnergySolid = " << freeEnergySolid << endl;
				dataFile << endl;
				dataFile << "successSolid = " << successSolid << endl;
				
				dataFile.close();
			}
		}
		
		// clean weights
		
		int sysresult = system("rm weights*");
	}
	
	//////////////////// Eliminate the parameter Npoints ///////////////////
	
	// minimise over N points and compute associated fluid free energy
	
	for (int j=0; j<kTSize; j++)
	{
		for (int k=0; k<muSize; k++)
		{
			// read data in file
			
			vector<double> density(NpointsSize,0);
			vector<double> numParticles(NpointsSize,0);
			vector<double> freeEnergy(NpointsSize,0);
			vector<bool> successPerNpoint(NpointsSize,0);
			
			for (int i=0; i<NpointsSize; i++)
			{
				ifstream dataFile(dataDir1+"/results_"+to_string(i)+
				                                   "_"+to_string(j)+
				                                   "_"+to_string(k)+".dat");
				
				double density2 = 0;
				double numParticles2 = 0;
				double freeEnergy2 = 0;
				int success = 0;
				
				readDataFromFile(dataFile, "densitySolid", density2);
				readDataFromFile(dataFile, "numParticlesSolid", numParticles2);
				readDataFromFile(dataFile, "freeEnergySolid", freeEnergy2);
				readDataFromFile(dataFile, "successSolid", success);
				
				density[i] = density2;
				numParticles[i] = numParticles2;
				freeEnergy[i] = freeEnergy2;
				successPerNpoint[i] = (bool) success;
				
				dataFile.close();
			}
			
			// find the data point that minimise the free energy
			
			int index_min = 0;
			for (int i=0; i<NpointsSize; i++)
				if (freeEnergy[i]<freeEnergy[index_min] && successPerNpoint[i]) 
					index_min = i;
			
			// select a few points in the vicinity of the min data point
			
			bool successSolid = true;
			
			vector<int> indicesToFit;
			for (int i=-1; i<2; i++)
			{
				if (index_min+i>=0 && index_min+i<NpointsSize && successPerNpoint[index_min+i])
					indicesToFit.push_back(index_min+i);
				else 
				{
					successSolid = false;
					log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
					log << "ERROR: The set of lattice sizes is inappropriate, the minimum is either out of the set or too close to its boundaries" << endl;
					//log << "DEBUG: index_min                     is " << index_min << endl;
					//log << "DEBUG: i                             is " << i << endl;
					//log << "DEBUG: index_min+i>=0                is " << (index_min+i>=0) << endl;
					//log << "DEBUG: index_min+i<NpointsSize       is " << (index_min+i<NpointsSize) << endl;
					//log << "DEBUG: successPerNpoint[index_min+i] is " << successPerNpoint[index_min+i] << endl;
				}
			}
			int numIndicesToFit = indicesToFit.size();
			
			double Npoints_min, density_min, numParticles_min, freeEnergy_min;
			
			if (successSolid)
			{
				// find the Npoints that minimises the free energy and the minimum value
				
				arma::vec x(numIndicesToFit);
				arma::vec y(numIndicesToFit);
				
				for (int i=0; i<numIndicesToFit; i++) x(i) = Npoints[indicesToFit[i]];
				for (int i=0; i<numIndicesToFit; i++) y(i) = freeEnergy[indicesToFit[i]];
				
				arma::vec polyFreeEnergy = arma::polyfit(x,y,2);
				
				double Npoints_min = -1.0/2*polyFreeEnergy(1)/polyFreeEnergy(0);
				
				freeEnergy_min = polyFreeEnergy(0)*Npoints_min*Npoints_min +
				                 polyFreeEnergy(1)*Npoints_min + polyFreeEnergy(2);
				
				// find the corresponding values of the density and number of particles
				
				for (int i=0; i<numIndicesToFit; i++) y(i) = density[indicesToFit[i]];
				arma::vec polyDensity = arma::polyfit(x,y,2);
				
				density_min = polyDensity(0)*Npoints_min*Npoints_min +
				              polyDensity(1)*Npoints_min + polyDensity(2);
				
				for (int i=0; i<numIndicesToFit; i++) y(i) = numParticles[indicesToFit[i]];
				arma::vec polyNumParticles = arma::polyfit(x,y,2);
				
				numParticles_min = polyNumParticles(0)*Npoints_min*Npoints_min +
				                   polyNumParticles(1)*Npoints_min + polyNumParticles(2);
			}
			
			// Report
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "(Npoints) Min Npoints = " << Npoints_min << endl;
			log << "(Npoints) Min Free Energy = " << freeEnergy_min << endl;
			log << "(Npoints) Min Density = " << density_min << endl;
			log << "(Npoints) Min Number of Particles = " << numParticles_min << endl;
			
			// Compute associated fluid properties
			
			// TODO: use bisection method?? 
			//       current problem if we are in vapor phase and start in the liquid
			
			bool successFluid = true;
			
			double densityFluid = findRootdOmegadRho(1.5, kT[j], mu[k], aVdW, hsd, successFluid);
			double freeEnergyFluid = uniformFluidOmega(kT[j], mu[k], aVdW, hsd, densityFluid);
			
			if (successFluid)
			{
				log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
				log << "Fluid Density = " << densityFluid << endl;
				log << "Fluid Free Energy = " << freeEnergyFluid << endl;
			}
			else
			{
				log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
				log << "ERROR: Fluid computation FAILED" << endl;
			}
			
			// Save in File
			
			ofstream dataFile(dataDir2+"/results_"+to_string(j)+
			                                   "_"+to_string(k)+".dat");
			
			dataFile << "kT = " << kT[j] << endl;
			dataFile << "mu = " << mu[k] << endl;
			dataFile << endl;
			dataFile << "densitySolid = " << density_min << endl;
			dataFile << "numParticlesSolid = " << numParticles_min << endl;
			dataFile << "freeEnergySolid = " << freeEnergy_min << endl;
			dataFile << endl;
			dataFile << "densityFluid = " << densityFluid << endl;
			dataFile << "freeEnergyFluid = " << freeEnergyFluid << endl;
			dataFile << endl;
			dataFile << "successSolid = " << successSolid << endl;
			dataFile << "successFluid = " << successFluid << endl;
			dataFile << "success      = " << (successSolid && successFluid) << endl;
			
			dataFile.close();
		}
	}
	
	/////////////// Plot free energy profiles with Npoints /////////////////
	
	// make data file with plot data
	
	for (int j=0; j<kTSize; j++)
	{
		for (int k=0; k<muSize; k++)
		{
			ofstream plotFile(plotDir+"/freeEnergy_"+to_string(j)+
			                                     "_"+to_string(k)+".dat");
			
			plotFile << "# kT=" << kT[j] << " mu=" << mu[k] << endl;
			plotFile << "# Npoints freeEnergy" << endl;
			plotFile << endl;
			
			for (int i=0; i<NpointsSize; i++)
			{
				ifstream dataFile(dataDir1+"/results_"+to_string(i)+
				                                   "_"+to_string(j)+
				                                   "_"+to_string(k)+".dat");
				
				double freeEnergy = 0;
				readDataFromFile(dataFile, "freeEnergySolid", freeEnergy);
				
				dataFile.close();
				
				plotFile << Npoints[i] << " " << freeEnergy << endl;
			}
			
			plotFile.close();
		}
	}
	
	// create the gnuplot script 
	
	ofstream plotFile(plotDir+"/plot");
	
	for (int j=0; j<kTSize; j++)
	{
		for (int k=0; k<muSize; k++)
		{
			plotFile << "#----------------------------------------" << endl;
			plotFile << "set term svg enhanced mouse #size 600,500" << endl;
			plotFile << "set output 'freeEnergy_" << j << "_" << k << ".svg'" << endl;
			plotFile << endl;
			plotFile << "set title \"Free Energy Profile with Npoints\" font \",20\"" << endl;
			plotFile << "set label \"kT = " << kT[j] << " mu = " << mu[k] << "\" at graph 0.1,0.9 font \",16\"" << endl; //0.29,1.05 is under title
			plotFile << "set xlabel \"Npoints\" font \",20\"" << endl;
			plotFile << "set ylabel \"Omega/(kT*V)\" font \",20\"" << endl;
			plotFile << endl;
			plotFile << "set key off" << endl;
			plotFile << endl;
			plotFile << "plot \"freeEnergy_" << to_string(j) << "_" 
			         << to_string(k) << ".dat\" with lines" << endl;
			plotFile << endl;
			
		}
	}
	
	plotFile.close();
	
	return 1;
}

