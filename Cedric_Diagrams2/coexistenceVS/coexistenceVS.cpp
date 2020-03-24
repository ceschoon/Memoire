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
                   vector<double> &freeEnergyVapour_vec,
                   vector<double> &densityVapour_vec)
{
	// Compute free energy for the solid
	
	double freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min;
	
	int statusSolid = fixedkTMuSolid(kT, mu, argc, argv, log, 
		freeEnergySolid, densitySolid, Ngrid_min, Cvac_min, alpha_min);
	
	// Compute free energy for the vapour
	
	double freeEnergyVapour, densityVapour;
	double freeEnergyLiquid, densityLiquid;
	bool supercritical;
	
	int statusVapour = fixedkTMuFluid2(kT, mu, argc, argv, log, 
		freeEnergyVapour, densityVapour, freeEnergyLiquid, densityLiquid,
		supercritical);
	
	// Save results
	
	freeEnergyDiff = freeEnergyVapour-freeEnergySolid;
	
	if (statusSolid==0 && statusVapour==0)
	{
		mu_vec.push_back(mu);
		freeEnergySolid_vec.push_back(freeEnergySolid);
		densitySolid_vec.push_back(densitySolid);
		freeEnergyVapour_vec.push_back(freeEnergyVapour);
		densityVapour_vec.push_back(densityVapour);
	}
	
	// Print results in log file
	
	if (statusSolid==0 && statusVapour==0)
	{
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Solid Free Energy = " << freeEnergySolid << endl;
		log << "Solid Density     = " << densitySolid << endl;
		log << "Vapour Free Energy = " << freeEnergyVapour << endl;
		log << "Vapour Density     = " << densityVapour << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		if (statusSolid!=0)
			log << "ERROR: Solid Computation FAILED" << endl;
		if (statusVapour!=0)
			log << "ERROR: Vapour Computation FAILED" << endl;
	}
	
	// return status
	
	if (statusSolid==0 && statusVapour==0) return 0;
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
		ifstream dataInFile(dataDir+"/coexistenceVS_"+"kT="+sskT.str()+".dat");
		
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
		vector<double> freeEnergyVapour_vec;
		vector<double> densityVapour_vec;
		
		//// start with our current guess and see in which direction we need
		//// to follow to reach the coexistence
		
		double mu_0 = muCoexGuess;
		double mu_1 = muCoexGuess + muStep; 
		
		double freeEnergyDiff_0; // vapour - solid
		double freeEnergyDiff_1;
		
		// TODO: currently not robust to failures !!!
		int status_0 = muComputation(kT, mu_0, argc, argv, log, 
			freeEnergyDiff_0, mu_vec, freeEnergySolid_vec, densitySolid_vec,
			freeEnergyVapour_vec, densityVapour_vec);
		int status_1 = muComputation(kT, mu_1, argc, argv, log, 
			freeEnergyDiff_1, mu_vec, freeEnergySolid_vec, densitySolid_vec,
			freeEnergyVapour_vec, densityVapour_vec);
		 
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
					freeEnergyVapour_vec, densityVapour_vec);
				
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
			              freeEnergyVapour_vec, densityVapour_vec);
			
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
			
			temp = freeEnergyVapour_vec[i]; 
			freeEnergyVapour_vec[i]     = freeEnergyVapour_vec[j_min]; 
			freeEnergyVapour_vec[j_min] = temp;
			
			temp = densityVapour_vec[i]; 
			densityVapour_vec[i]     = densityVapour_vec[j_min]; 
			densityVapour_vec[j_min] = temp;
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
		log << "freeEnergyVapour_vec = ";
		for (int i=0; i<freeEnergyVapour_vec.size(); i++)
			log << freeEnergyVapour_vec[i] << " ";
		
		log << endl;
		log << "densityVapour_vec = ";
		for (int i=0; i<densityVapour_vec.size(); i++)
			log << densityVapour_vec[i] << " ";
		
		// find the coexistence value for mu by searching the zero of
		// freeEnergyVapour - freeEnergySolid (using interpolation)
		
		vector<double> freeEnergyDiff_vec(mu_vec.size()); 
		for (int i=0; i<mu_vec.size(); i++)
			freeEnergyDiff_vec[i] = freeEnergyVapour_vec[i] - freeEnergySolid_vec[i];
		
		bool foundCoex = true;
		double muCoex = 0;
		double freeEnergyCoex = 0;
		double densitySolidCoex = 0;
		double densityVapourCoex = 0;
		
		// check if has the same sign before searching the zero
		if (freeEnergyDiff_vec[0] * freeEnergyDiff_vec[freeEnergyDiff_vec.size()-1] < 0)
		{
			// find min mu
			zeroFromDataInterpolation(mu_vec, freeEnergyDiff_vec, muCoex);
			
			// find min free energy
			evalFromDataInterpolation(mu_vec, freeEnergySolid_vec, muCoex, freeEnergyCoex);
			
			// find min densities
			evalFromDataInterpolation(mu_vec, densitySolid_vec, muCoex, densitySolidCoex);
			evalFromDataInterpolation(mu_vec, densityVapour_vec, muCoex, densityVapourCoex);
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Range invalid for kT = " << kT 
			    << ". It does not contain the coexistence mu" << endl;
			
			foundCoex = false;
		}
		
		/////////////////////// Store the result ///////////////////////////
		
		ofstream dataFile(dataDir+"/coexistenceVS_"+"kT="+sskT.str()+".dat");
		
		dataFile << "# Result for vapour-solid coexistence at " << endl;
		dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
		dataFile << endl;
		dataFile << "muCoex = " << scientific << setprecision(8) << muCoex << endl;
		dataFile << "freeEnergyCoex = " << scientific << setprecision(8) << freeEnergyCoex << endl;
		dataFile << "densitySolidCoex = " << scientific << setprecision(8) << densitySolidCoex << endl;
		dataFile << "densityVapourCoex = " << scientific << setprecision(8) << densityVapourCoex << endl;
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