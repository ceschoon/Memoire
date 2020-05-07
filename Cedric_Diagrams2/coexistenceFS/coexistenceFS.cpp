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
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

#include "../SolidFCC.h"
#include "../SolidPhase.h"
#include "../UniformPhases.h"
#include "../utilities.h"


// computation task to perform for each (kT,mu)

int muComputation( double kT, double mu, 
                   int argc, char **argv, Log &log,
                   double &freeEnergyDiff, 
                   vector<double> &mu_vec, 
                   vector<double> &freeEnergySolid_vec,
                   vector<double> &densitySolid_vec,
                   vector<double> &Ngrid_vec,
                   vector<double> &Cvac_vec,
                   vector<double> &alpha_vec,
                   vector<double> &freeEnergyFluid_vec,
                   vector<double> &densityFluid_vec)
{
	// Compute free energy for the solid
	
	double freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min;
	
	int statusSolid = fixedkTMuSolid(kT, mu, argc, argv, log, 
		freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min);
	
	// Compute free energy for the fluid
	
	double freeEnergyFluid, densityFluid;
	
	int statusFluid = fixedkTMuFluid(kT, mu, argc, argv, log, 
		freeEnergyFluid, densityFluid);
	
	// Save results
	
	freeEnergyDiff = freeEnergyFluid-freeEnergySolid;
	
	if (statusSolid==0 && statusFluid==0)
	{
		mu_vec.push_back(mu);
		freeEnergySolid_vec.push_back(freeEnergySolid);
		densitySolid_vec.push_back(densitySolid);
		Ngrid_vec.push_back(Ngrid_min);
		Cvac_vec.push_back(Cvac_min);
		alpha_vec.push_back(alpha_min);
		freeEnergyFluid_vec.push_back(freeEnergyFluid);
		densityFluid_vec.push_back(densityFluid);
	}
	
	// Print results in log file
	
	if (statusSolid==0 && statusFluid==0)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
		log << "Solid Density     = " << densitySolid << endl;
		log << "Ngrid = " << Ngrid_min << endl;
		log << "Cvac  = " << Cvac_min  << endl;
		log << "alpha = " << alpha_min << endl;
		log << "Fluid Free Energy = " << freeEnergyFluid << endl;
		log << "Fluid Density     = " << densityFluid << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		if (statusSolid!=0)
			log << "ERROR: Solid Computation FAILED" << endl;
		if (statusFluid!=0)
			log << "ERROR: Fluid Computation FAILED" << endl;
	}
	
	// return status
	
	if (statusSolid==0 && statusFluid==0) return 0;
	return 1;
}







int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDir = "data";
	double kTStart;
	double kTStop;
	double kTStep;
	double muCoexGuess;
	double muRange;
	double muStep;
	
	Options options;
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kTStart", &kTStart);
	options.addOption("kTStop", &kTStop);
	options.addOption("kTStep", &kTStep);
	options.addOption("muCoexGuess", &muCoexGuess);
	options.addOption("muRange", &muRange);
	options.addOption("muStep", &muStep);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	////////////////////////////// Loop over kT ////////////////////////////
	
	double kT = kTStart;
	
	while ( (kTStep>0 && kT<=kTStop) || // stop condition if kTStep>0
	        (kTStep<0 && kT>=kTStop) )  // stop condition if kTStep<0
	{
		////////// Check if the computation has been done already //////////
		
		log << endl;
		log << "Checking if the computation has been done already" << endl;
		
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		ifstream dataInFile(dataDir+"/coexistenceFS_"+"kT="+sskT.str()+".dat");
		
		// if file exists, get the result from there
		if (dataInFile)
		{
			// check if it was a successful computation
			int success;
			readDataFromFile(dataInFile, "success", success);
			
			if (success==1)  // bool "true" 
			{
				log << endl;
				log << "Skipping existing computation" << endl;
				
				// update muCoexGuess based on file's result
				readDataFromFile(dataInFile, "muCoex", muCoexGuess);
				
				log << endl;
				log << "read muCoex = " << muCoexGuess << endl;
				
				kT += kTStep;
				
				continue;
			}
		}
		
		
		////////////// Free energies computation for this kT ///////////////
		
		bool crossedCoex = false;
		vector<double> mu_vec;
		vector<double> freeEnergySolid_vec;
		vector<double> densitySolid_vec;
		vector<double> Ngrid_vec;
		vector<double> Cvac_vec;
		vector<double> alpha_vec;
		vector<double> freeEnergyFluid_vec;
		vector<double> densityFluid_vec;
		
		//// start with our current guess and see in which direction we need
		//// to follow to reach the coexistence
		
		double mu_0 = muCoexGuess;
		double mu_1 = muCoexGuess + muStep; 
		
		double freeEnergyDiff_0; // fluid - solid
		double freeEnergyDiff_1;
		
		// TODO: currently not robust to failures !!!
		int status_0 = muComputation(kT, mu_0, argc, argv, log, 
			freeEnergyDiff_0, mu_vec, freeEnergySolid_vec, densitySolid_vec,
			Ngrid_vec, Cvac_vec,alpha_vec,
			freeEnergyFluid_vec, densityFluid_vec);
		int status_1 = muComputation(kT, mu_1, argc, argv, log, 
			freeEnergyDiff_1, mu_vec, freeEnergySolid_vec, densitySolid_vec,
			Ngrid_vec, Cvac_vec,alpha_vec,
			freeEnergyFluid_vec, densityFluid_vec);
		 
		double mu_current,freeEnergyDiff_current;
		
		log << endl;
		log << "mu_0 = " << mu_0 << " freeEnergyDiff_0 = "
		    << freeEnergyDiff_0 << " status_0 = " << status_0 << endl;
		log << "mu_1 = " << mu_1 << " freeEnergyDiff_1 = " 
		    << freeEnergyDiff_1 << " status_1 = " << status_1 << endl;
		
		// if wrong direction
		//   that is if we go away from freeEnergyDiff_0 = 0 (coex)
		//   or if only computation "1" fails (in that case, it is likely 
		//   that the comptations will work in direction of "0"
		if ( (status_0==0 && status_1==0 && freeEnergyDiff_0 < 0 && freeEnergyDiff_1 < freeEnergyDiff_0) ||
		     (status_0==0 && status_1==0 && freeEnergyDiff_0 > 0 && freeEnergyDiff_1 > freeEnergyDiff_0) ||
		     (status_0==0 && status_1!=0) ) 
		{
			log << "Going in the wrong direction, flipping muStep" << endl;
			
			muStep *= -1;
			mu_current = mu_0;
			freeEnergyDiff_current = freeEnergyDiff_0;
		}
		// if right direction
		//   note that if both computation fails, we continue to search in 
		//   the default direction
		else
		{
			log << "Going in the right direction" << endl;
			
			mu_current = mu_1;
			freeEnergyDiff_current = freeEnergyDiff_1;
			
			// if we already crossed the coexistence point, do one more 
			// computation as only two point will not be enough to do a 
			// precise interpolation 
			if (freeEnergyDiff_1*freeEnergyDiff_0<0)
			{
				log << "Already crossed the coexistence point" << endl;
				
				mu_current += muStep;
				muComputation(kT, mu_current, argc, argv, log, freeEnergyDiff_current, 
					mu_vec, freeEnergySolid_vec, densitySolid_vec,
					Ngrid_vec, Cvac_vec,alpha_vec,
					freeEnergyFluid_vec, densityFluid_vec);
				
				crossedCoex = true;
			}
		}
		
		//// continue in the right direction untill we cross 
		//// freeEnergyDiff = 0 or go too far relative to a given range
		
		while (abs(mu_current-muCoexGuess)<muRange && !crossedCoex)
		{
			double mu_old = mu_current;
			double freeEnergyDiff_old = freeEnergyDiff_current;
			mu_current += muStep;
			
			// TODO: currently not robust to failures !!!
			muComputation(kT, mu_current, argc, argv, log, 
			              freeEnergyDiff_current, mu_vec, 
			              freeEnergySolid_vec, densitySolid_vec,
			              Ngrid_vec, Cvac_vec,alpha_vec,
			              freeEnergyFluid_vec, densityFluid_vec);
			
			// stop if we cross coexistence point (change in sign)
			// and if we have at least three valid computations
			if ( (freeEnergyDiff_old*freeEnergyDiff_current<0) &&
			     (mu_vec.size()>=3) )
				crossedCoex = true;
		}
		
		
		///////////////////// Extract coexistence mu ///////////////////////
		
		// sort by increasing order in mu
		
		for (int i=0; i<mu_vec.size(); i++)
		{
			// find index of the min in the remaining j=i,i+1,... indices
			int j_min = i; 
			for (int j=i; j<mu_vec.size(); j++)
				if (mu_vec[j]<mu_vec[j_min]) j_min = j;
			
			// swap index i and j_min
			double temp = mu_vec[i]; 
			mu_vec[i]     = mu_vec[j_min]; 
			mu_vec[j_min] = temp;
			
			temp = freeEnergySolid_vec[i]; 
			freeEnergySolid_vec[i]     = freeEnergySolid_vec[j_min]; 
			freeEnergySolid_vec[j_min] = temp;
			
			temp = densitySolid_vec[i]; 
			densitySolid_vec[i]     = densitySolid_vec[j_min]; 
			densitySolid_vec[j_min] = temp;
			
			temp = Ngrid_vec[i]; 
			Ngrid_vec[i]     = Ngrid_vec[j_min]; 
			Ngrid_vec[j_min] = temp;
			
			temp = Cvac_vec[i]; 
			Cvac_vec[i]     = Cvac_vec[j_min]; 
			Cvac_vec[j_min] = temp;
			
			temp = alpha_vec[i]; 
			alpha_vec[i]     = alpha_vec[j_min]; 
			alpha_vec[j_min] = temp;
			
			temp = freeEnergyFluid_vec[i]; 
			freeEnergyFluid_vec[i]     = freeEnergyFluid_vec[j_min]; 
			freeEnergyFluid_vec[j_min] = temp;
			
			temp = densityFluid_vec[i]; 
			densityFluid_vec[i]     = densityFluid_vec[j_min]; 
			densityFluid_vec[j_min] = temp;
		}
		
		// report computations
		
		log << endl;
		log << "mu_vec = ";
		for (int i=0; i<mu_vec.size(); i++)
			log << mu_vec[i] << " ";
		
		log << endl;
		log << "freeEnergySolid_vec = ";
		for (int i=0; i<freeEnergySolid_vec.size(); i++)
			log << freeEnergySolid_vec[i] << " ";
		
		log << endl;
		log << "densitySolid_vec = ";
		for (int i=0; i<densitySolid_vec.size(); i++)
			log << densitySolid_vec[i] << " ";
		
		log << endl;
		log << "Ngrid_vec = ";
		for (int i=0; i<Ngrid_vec.size(); i++)
			log << Ngrid_vec[i] << " ";
		
		log << endl;
		log << "Cvac_vec = ";
		for (int i=0; i<Cvac_vec.size(); i++)
			log << Cvac_vec[i] << " ";
		
		log << endl;
		log << "alpha_vec = ";
		for (int i=0; i<alpha_vec.size(); i++)
			log << alpha_vec[i] << " ";
		
		log << endl;
		log << "freeEnergyFluid_vec = ";
		for (int i=0; i<freeEnergyFluid_vec.size(); i++)
			log << freeEnergyFluid_vec[i] << " ";
		
		log << endl;
		log << "densityFluid_vec = ";
		for (int i=0; i<densityFluid_vec.size(); i++)
			log << densityFluid_vec[i] << " ";
		
		// find the coexistence value for mu by searching the zero of
		// freeEnergyFluid - freeEnergySolid (using interpolation)
		
		vector<double> freeEnergyDiff_vec(mu_vec.size()); 
		for (int i=0; i<mu_vec.size(); i++)
			freeEnergyDiff_vec[i] = freeEnergyFluid_vec[i] - freeEnergySolid_vec[i];
		
		bool foundCoex = true;
		double muCoex = 0;
		double freeEnergyCoex = 0;
		double densitySolidCoex = 0;
		double NgridCoex = 0;
		double CvacCoex = 0;
		double alphaCoex = 0;
		double densityFluidCoex = 0;
		
		// check if has the same sign before searching the zero
		if (freeEnergyDiff_vec[0] * freeEnergyDiff_vec[freeEnergyDiff_vec.size()-1] < 0)
		{
			// find min mu
			zeroFromDataInterpolation(mu_vec, freeEnergyDiff_vec, muCoex);
			
			// find min free energy
			evalFromDataInterpolation(mu_vec, freeEnergySolid_vec, muCoex, freeEnergyCoex);
			
			// find min densities
			evalFromDataInterpolation(mu_vec, densitySolid_vec, muCoex, densitySolidCoex);
			evalFromDataInterpolation(mu_vec, Ngrid_vec, muCoex, NgridCoex);
			evalFromDataInterpolation(mu_vec, Cvac_vec, muCoex, CvacCoex);
			evalFromDataInterpolation(mu_vec, alpha_vec, muCoex, alphaCoex);
			evalFromDataInterpolation(mu_vec, densityFluid_vec, muCoex, densityFluidCoex);
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Range invalid for kT = " << kT 
			    << ". It does not contain the coexistence mu" << endl;
			
			foundCoex = false;
		}
		
		/////////////////////// Store the result ///////////////////////////
		
		ofstream dataFile(dataDir+"/coexistenceFS_"+"kT="+sskT.str()+".dat");
		
		dataFile << "# Result for fluid-solid coexistence at " << endl;
		dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
		dataFile << endl;
		dataFile << "muCoex = " << scientific << setprecision(8) << muCoex << endl;
		dataFile << "freeEnergyCoex = " << scientific << setprecision(8) << freeEnergyCoex << endl;
		dataFile << "densitySolidCoex = " << scientific << setprecision(8) << densitySolidCoex << endl;
		dataFile << "NgridCoex = " << scientific << setprecision(8) << NgridCoex << endl;
		dataFile << "CvacCoex = " << scientific << setprecision(8) << CvacCoex << endl;
		dataFile << "alphaCoex = " << scientific << setprecision(8) << alphaCoex << endl;
		dataFile << "densityFluidCoex = " << scientific << setprecision(8) << densityFluidCoex << endl;
		dataFile << endl;
		dataFile << "success = " << foundCoex << endl;
		
		dataFile.close();
		
		////////////////////////////////////////////////////////////////////
		
		// increment
		kT += kTStep;
		
		// update muCoexGuess
		if (foundCoex) muCoexGuess = muCoex;
	}
	
	
	
	
	
	return 0;
}