////////////////////////////////////////////////////////////////////////////
//                                                                        //
//                 DFT Computations for the Solid Phase                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <dirent.h>

using namespace std; 

#ifdef USE_OMP
#include <omp.h>
#endif

#include <armadillo>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_multimin.h>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidDensity.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"
#include "utilities.h"


// Parameter space is (kT,mu,Ngrid,Cvac,alpha),
// where kT     is the temperature
//       mu     is the chemical potential (unit kT)
//       Ngrid  is the number of grid points on the cell side
//       Cvac   is the concentration of vacancies
//       alpha  is the width parameters of the gaussian density profile
// (kT,mu)            thermodynamics parameters
// (Ngrid,Cvac,alpha) parameters of the density profile with respect to
//                    which we minimise the free energy functional


int fixedkTMuSolid(double kT, double mu,
                   int argc, char** argv, Log &log,
                   double &freeEnergy_min, double &density_min, 
                   double &Ngrid_min, double &Cvac_min, double &alpha_min);


int minOverAlphaOnly(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &alpha_min);


int minOverCvacAlpha(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &Cvac_min, double &alpha_min);


int DFTgaussian( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density);




////////////////////////////////////////////////////////////////////////////

// Function to find minimum of free energy in (Ngrid,Cvac,alpha)-space.
// The plan is:
// o read options
//    get ranges for Ngrid
// o try to use prev result for (kT,mu) that is close
//    if none available, find valid Ngrid and then search for Ngrid_min 
// o calc min (Cvac,alpha) for each Ngrid selected
// o find minimum
// o store result somewhere

// TODO: ugly but works, could be re-written

int fixedkTMuSolid(
	double kT, double mu,                        // input
	int argc, char** argv, Log &log,             // options for algo behaviour
	double &freeEnergy_min, double &density_min, // result 
	double &Ngrid_min, double &Cvac_min, double &alpha_min)
{
	///// log /////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Computation of free energy at given (kT,mu) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	
	///// read options /////
	
	string dataDir = "data";
	double kT_rangePrevResult = 0.03;
	double mu_rangePrevResult = 0.15;
	int Ngrid_rangeMin;
	int Ngrid_rangeMax;
	int Ngrid_rangeStep;
	int Ngrid_maxNumTrials = 4;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kT_rangePrevResult", &kT_rangePrevResult);
	options.addOption("mu_rangePrevResult", &mu_rangePrevResult);
	options.addOption("Ngrid_rangeMin", &Ngrid_rangeMin);
	options.addOption("Ngrid_rangeMax", &Ngrid_rangeMax);
	options.addOption("Ngrid_rangeStep", &Ngrid_rangeStep);
	options.addOption("Ngrid_maxNumTrials", &Ngrid_maxNumTrials);
	
	options.read(argc, argv);
	
	int sysres = system(("mkdir -p "+dataDir).c_str());
	
	///// Check if computation has not already been done /////
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	
	ifstream dataInFile(dataDir+"/kTMuSolid_"+"kT="+sskT.str()+"_mu="
	                    +ssMu.str()+".dat");
	
	// if file exists, get the result from there
	if (dataInFile)
	{
		// check if it was a successful computation
		int success;
		readDataFromFile(dataInFile, "success", success);
		
		if (success==1)  // bool "true" 
		{
			log << endl;
			log << "Returning existing computation" << endl;
			
			// read results in file
			readDataFromFile(dataInFile, "freeEnergy_min", freeEnergy_min);
			readDataFromFile(dataInFile, "density_min", density_min);
			readDataFromFile(dataInFile, "Ngrid_min", Ngrid_min);
			readDataFromFile(dataInFile, "Cvac_min", Cvac_min);
			readDataFromFile(dataInFile, "alpha_min", alpha_min);
			
			log << endl;
			log << "freeEnergy_min = " << freeEnergy_min << endl;
			log << "density_min = " << density_min << endl;
			log << "Ngrid_min = " << Ngrid_min << endl;
			log << "Cvac_min = " << Cvac_min << endl;
			log << "alpha_min = " << alpha_min << endl;
			
			return 0;
		}
		
		if (success==0)  // bool "false" 
		{
			log << endl;
			log << "Returning existing computation (is a failed one)" << endl;
			
			return 1;
		}
	}
	
	
	///// By default, store the result as a failure /////
	// It will be overwritten in the case of a successful computation
	
	sskT.str(""); sskT << scientific << setprecision(4) << kT;
	ssMu.str(""); ssMu << scientific << setprecision(4) << mu;
	
	ofstream dataFile(dataDir+"/kTMuSolid_"+"kT="+sskT.str()+"_mu="
	                  +ssMu.str()+".dat");
	
	dataFile << "# Result for solid computation at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << endl;
	dataFile << "success = " << false << endl;
	
	dataFile.close();
	
	///// Initialisations /////
	
	struct CompResult 
	{
		int Ngrid;
		double freeEnergy;
		double density;
		double Cvac;
		double alpha;
	};
	
	vector<CompResult> results;
	
	int Ngrid_intMin; // Ngrid point with lowest freeEnergy (integer)
	                  // Ngrid_min will be the minimum of an interpolation
	                  // (double float)
	
	
	///////// First, find a value of Ngrid that gives a successful /////////
	/////////// computation of the free energy -> starting point ///////////
	
	///// try to use prev result for (kT,mu) that is close /////
	
	bool usePreviousResult = false;
	
	log << endl;
	log << "Trying to find a value of Ngrid that leads to a valid computation" << endl;
	log << "Start looking in previous results." << endl;
	
	DIR *dir;
	struct dirent *ent;

	if ((dir = opendir (dataDir.c_str())) != NULL) 
	{
		// print all the candidates files within directory 
		
		while ((ent = readdir (dir)) != NULL) 
		{
			string ent_name(ent->d_name);
			
			if (ent_name.size()>=9 && ent_name.substr(0,9) == "kTMuSolid")
			{
				log << "Checking file: " << ent_name << endl;
				
				// get kT and mu from file name
				
				int pos1_=-1; // position of underscore in name (_kT=)
				int pos2_=-1; // position of underscore in name (_mu=)
				
				for (int i=0; i<ent_name.size(); i++)
				{
					if (ent_name[i]=='_' && pos1_<0) pos1_ = i;
					else if (ent_name[i]=='_' && pos2_<0) pos2_ = i;
				}
				
				int lenkTstr = pos2_-pos1_-4;
				int lenmustr = ent_name.size()-4-pos2_-4;
				
				string kT_str = ent_name.substr(pos1_+4,lenkTstr);
				string mu_str = ent_name.substr(pos2_+4,lenmustr);
				
				double kT_read = stod(kT_str);
				double mu_read = stod(mu_str);
				
				// check if (kT,mu) in file are close enough to be used
				
				if ( (abs(kT-kT_read)<kT_rangePrevResult) ||
				     (abs(mu-mu_read)<mu_rangePrevResult)   )
				{
					ifstream dataInFile(dataDir+"/"+ent_name);
					if (dataInFile)
					{
						int success_read;
						double Ngrid_read;
						
						readDataFromFile(dataInFile, "success", success_read);
						
						if (success_read==1) // here 1 means true, thus success
						{
							readDataFromFile(dataInFile, "Ngrid_min", Ngrid_read);
							int Ngrid_candidate = int(Ngrid_read+0.5);
							
							log << "Trying file Ngrid=" << Ngrid_candidate << endl;
							
							double Cvac, alpha, freeEnergy, density;
							int status_minOverCvacAlpha = 
								minOverCvacAlpha(kT, mu, Ngrid_candidate,
								argc, argv, log, freeEnergy, density, Cvac, alpha);
							
							if (status_minOverCvacAlpha==0) // success
							{
								log << "success! freeEnergy=" << freeEnergy << endl;
								
								struct CompResult res;
								res.Ngrid = Ngrid_candidate;
								res.freeEnergy = freeEnergy;
								res.density = density;
								res.Cvac = Cvac;
								res.alpha = alpha;
								
								results.push_back(res);
								log << "File Ngrid=" << res.Ngrid << " accepted as starting point" << endl;
								usePreviousResult = true;
								break;
							}
							else
							{
								log << "computation not successful" << endl;
							}
						}
					}
					else
					{
						// could not open directory 
						log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
						log << "ERROR: could not open data file \"" << dataDir+"/"+ent_name << "\"" << endl;
					}
				}
			}
		}
		
		if (!usePreviousResult)
			log << "Did not found useful data in previous results." << endl;
		
		closedir (dir);
	} 
	else 
	{
		// could not open directory 
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: could not open data directory (fixedkTMu)" << endl;
	}
	
	
	///// if no previous result is available, use a random search to find a 
	///// value of Ngrid that leads to a successful computation.
	///// Try only a few ones, if these few ones fail, we stop trying.
	///// This saves a lot of computation time instead of spending them on 
	///// useless computations.
		
	if (!usePreviousResult)
	{
		log << endl;
		log << "Trying to find a value of Ngrid that leads to a valid computation" << endl;
		
		vector<int> Ngrid_toTest;
		for (int Ngrid=Ngrid_rangeMin; Ngrid<Ngrid_rangeMax; Ngrid+=Ngrid_rangeStep)
			Ngrid_toTest.push_back(Ngrid);
		random_shuffle(Ngrid_toTest.begin(), Ngrid_toTest.end());
		
		int trialsCounter = 0;
		
		bool noSuccessfulComp = true;
		for (int Ngrid : Ngrid_toTest)
		{
			log << "Trying Ngrid=" << Ngrid << endl;
			
			trialsCounter++;
			if (trialsCounter>Ngrid_maxNumTrials) break;
			
			double Cvac, alpha, freeEnergy, density;
			int status_minOverCvacAlpha = minOverCvacAlpha(
				kT, mu, Ngrid, argc, argv, log, freeEnergy, density, 
				Cvac, alpha);
			
			if (status_minOverCvacAlpha==0) // success
			{
				log << "success! freeEnergy=" << freeEnergy << endl;
				
				struct CompResult res;
				res.Ngrid = Ngrid;
				res.freeEnergy = freeEnergy;
				res.density = density;
				res.Cvac = Cvac;
				res.alpha = alpha;
				
				results.push_back(res);
				Ngrid_intMin = Ngrid; // update current best Min
				noSuccessfulComp = false;
				break;
			}
		}
		
		if (noSuccessfulComp)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: No successful computation found for kT=" << kT << " and mu=" << mu << endl;
			
			return 1;
		}
	}
	
	///////////// Follow low free energies to reach a minimum //////////////
	
	// descent from the valid Ngrid that we found untill we
	// fall in the minimum (we assume only one exists)
	
	log << endl;
	log << "Trying to find the minimum over Ngrid from our first valid computation" << endl;
	
	bool foundMin = false;
	
	// here, there should be only one result in the vector
	struct CompResult lastRes = results[results.size()-1];
	
	// try neighbour +
	struct CompResult resPlus; resPlus.Ngrid = lastRes.Ngrid + Ngrid_rangeStep;
	log << "Trying Ngrid=" << resPlus.Ngrid << endl;
	
	double CvacPlus, alphaPlus, freeEnergyPlus, densityPlus;
	int status_Plus = minOverCvacAlpha(
		kT, mu, resPlus.Ngrid, argc, argv, log, freeEnergyPlus, densityPlus, 
		CvacPlus, alphaPlus);
	
	if (status_Plus==0)
	{
		log << "success! freeEnergy=" << freeEnergyPlus << endl;
		
		resPlus.freeEnergy = freeEnergyPlus;
		resPlus.density = densityPlus;
		resPlus.Cvac = CvacPlus;
		resPlus.alpha = alphaPlus;
		
		results.push_back(resPlus);
		
		if (freeEnergyPlus<lastRes.freeEnergy) 
			Ngrid_intMin = resPlus.Ngrid; // update current best Min
	}
	
	// try neighbour -
	struct CompResult resMinus; resMinus.Ngrid = lastRes.Ngrid - Ngrid_rangeStep;
	log << "Trying Ngrid=" << resMinus.Ngrid << endl;
	
	double CvacMinus, alphaMinus, freeEnergyMinus, densityMinus;
	int status_Minus = minOverCvacAlpha(
		kT, mu, resMinus.Ngrid, argc, argv, log, freeEnergyMinus, 
		densityMinus, CvacMinus, alphaMinus);
	
	if (status_Minus==0)
	{
		log << "success! freeEnergy=" << freeEnergyMinus << endl;
		
		resMinus.freeEnergy = freeEnergyMinus;
		resMinus.density = densityMinus;
		resMinus.Cvac = CvacMinus;
		resMinus.alpha = alphaMinus;
		
		results.push_back(resMinus);
		
		if (freeEnergyMinus<lastRes.freeEnergy) 
			Ngrid_intMin = resMinus.Ngrid; // update current best Min
	}
	
	// check if not already at the minimum
	if (freeEnergyPlus  > lastRes.freeEnergy && 
	    freeEnergyMinus > lastRes.freeEnergy &&
	    status_Plus ==0 && status_Minus == 0)
	{
		log << "Already at the minimum (Ngrid)" << endl;
		foundMin = true;
		Ngrid_intMin = lastRes.Ngrid; // update current best Min
	}
	
	// if not at the minimum, continue to search in the min direction
	int NgridStep = Ngrid_rangeStep;
	if (freeEnergyMinus<lastRes.freeEnergy && status_Minus==0) 
		NgridStep *= -1;
	
	while (!foundMin)
	{
		// try new Ngrid point, one step further than the current min
		int resMinIndex = 0;
		for (int i=0; i<results.size(); i++) 
			if (results[i].Ngrid == Ngrid_intMin) resMinIndex = i;
		
		lastRes = results[resMinIndex];
		struct CompResult resNew; resNew.Ngrid = lastRes.Ngrid + NgridStep;
		
		// check if we don't cross the range limits
		if (resNew.Ngrid > Ngrid_rangeMax || resNew.Ngrid < Ngrid_rangeMin)
			break;
		
		double CvacNew, alphaNew, freeEnergyNew, densityNew;
		int status_New = minOverCvacAlpha(
			kT, mu, resNew.Ngrid, argc, argv, log, freeEnergyNew, densityNew,
			CvacNew, alphaNew);
		
		if (status_New==0)
		{
			log << "success! freeEnergy=" << freeEnergyNew << endl;
			
			resNew.freeEnergy = freeEnergyNew;
			resNew.density = densityNew;
			resNew.Cvac = CvacNew;
			resNew.alpha = alphaNew;
			
			results.push_back(resNew);
			
			if (freeEnergyNew<lastRes.freeEnergy)
				Ngrid_intMin = resNew.Ngrid;  // update current best Min
			else
				foundMin = true; // the min is the previous one
		}
		else
		{
			log << "failure !" << endl;
			break;
		}
	}
	
	// At this point, we should have found a minimum unless the minimum
	// lies outside of the given range
	
	if (!foundMin) 
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "Failure: No minimum over Ngrid found for kT=" << kT << " and mu=" << mu << endl;
		log << "        (No minimum at the end of the descent)";
		
		return 1;
	}
	
	///// find minimum of interpolation around Ngrid_intMin /////
	
	// we make sure that we have at least 3 data points centered around 
	// Ngrid_intMin
	
	log << endl;
	log << "Preparing interpolation for Ngrid_intMin = " << Ngrid_intMin << endl;
	
	struct CompResult resMin;
	struct CompResult resMinPlus;
	struct CompResult resMinMinus;
	
	bool foundMinMinus = false;
	bool foundMinPlus  = false;
	
	for (int i=0; i<results.size(); i++)
	{
		struct CompResult res = results[i];
		if (res.Ngrid==Ngrid_intMin) 
		{
			resMin = res;
		}
		if (res.Ngrid==Ngrid_intMin+1) 
		{
			foundMinPlus = true;
			resMinPlus = res;
		}
		if (res.Ngrid==Ngrid_intMin-1) 
		{
			foundMinMinus = true;
			resMinMinus = res;
		}
	}
	
	// do the calculation if not found
	bool calcPlusFailed = true;
	if (!foundMinPlus)
	{
		resMinPlus.Ngrid = Ngrid_intMin + Ngrid_rangeStep;
		
		log << "Computing missing Ngrid=" << resMinPlus.Ngrid << endl;
		
		// check that we don't cross the range limits
		if (resMinPlus.Ngrid <= Ngrid_rangeMax &&
		    resMinPlus.Ngrid >= Ngrid_rangeMin)
		{
			double Cvac, alpha, freeEnergy, density;
			int status_minOverCvacAlpha = minOverCvacAlpha(
				kT, mu, resMinPlus.Ngrid, argc, argv, log, freeEnergy, 
				density, Cvac, alpha);
			
			if (status_minOverCvacAlpha==0)
			{
				log << "success! freeEnergy=" << freeEnergy << endl;
				
				resMinPlus.freeEnergy = freeEnergy;
				resMinPlus.density = density;
				resMinPlus.Cvac = Cvac;
				resMinPlus.alpha = alpha;
				
				calcPlusFailed = false;
			}
		}
	}
	
	if (!foundMinPlus && calcPlusFailed) // in that case, we cannot proceed to the minimisation
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: No minimum found for kT=" << kT << " and mu=" << mu << endl;
		log << "       (Failure in neighbour+ calc for interpolation)" << endl;
		
		return 1;
	}
	
	
	// do the calculation if not found
	bool calcMinusFailed = true;
	if (!foundMinMinus)
	{
		resMinMinus.Ngrid = Ngrid_intMin + Ngrid_rangeStep;
		
		log << "Computing missing Ngrid=" << resMinMinus.Ngrid << endl;
		
		// check that we don't cross the range limits
		if (resMinMinus.Ngrid <= Ngrid_rangeMax &&
		    resMinMinus.Ngrid >= Ngrid_rangeMin)
		{
			double Cvac, alpha, freeEnergy, density;
			int status_minOverCvacAlpha = minOverCvacAlpha(
				kT, mu, resMinMinus.Ngrid, argc, argv, log, freeEnergy, 
				density, Cvac, alpha);
			
			if (status_minOverCvacAlpha==0)
			{
				log << "success! freeEnergy=" << freeEnergy << endl;
				
				resMinMinus.freeEnergy = freeEnergy;
				resMinMinus.density = density;
				resMinMinus.Cvac = Cvac;
				resMinMinus.alpha = alpha;
				
				calcMinusFailed = false;
			}
		}
	}
	
	if (!foundMinMinus && calcMinusFailed) // in that case, we cannot proceed to the minimisation
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: No minimum found for kT=" << kT << " and mu=" << mu << endl;
		log << "       (Failure in neighbour- calc for interpolation)" << endl;
		
		return 1;
	}
	
	// store in vectors for minimisation
	
	vector<int> Ngrid_vec = {resMinMinus.Ngrid, 
	                         resMin.Ngrid, 
	                         resMinPlus.Ngrid};
	
	vector<double> freeEnergy_vec = {resMinMinus.freeEnergy, 
	                                 resMin.freeEnergy, 
	                                 resMinPlus.freeEnergy};
	
	vector<double> density_vec = {resMinMinus.density, 
	                              resMin.density, 
	                              resMinPlus.density};
	
	vector<double> Cvac_vec = {resMinMinus.Cvac, 
	                           resMin.Cvac, 
	                           resMinPlus.Cvac};
	
	vector<double> alpha_vec = {resMinMinus.alpha, 
	                            resMin.alpha, 
	                            resMinPlus.alpha};
	
	
	// find (double float) minimum from interpolation
	
	int statusMin = minFromDataCSpline(Ngrid_vec, freeEnergy_vec, 
		Ngrid_min, freeEnergy_min);
	//int statusMin = minFromDataParabola(Ngrid_vec, freeEnergy_vec, 
	//	Ngrid_min, freeEnergy_min);
	
	// Evaluate other quantities at the minimum
	
	statusMin += evalFromDataInterpolation(Ngrid_vec, density_vec, Ngrid_min, density_min);
	statusMin += evalFromDataInterpolation(Ngrid_vec, Cvac_vec, Ngrid_min, Cvac_min);
	statusMin += evalFromDataInterpolation(Ngrid_vec, alpha_vec, Ngrid_min, alpha_min);
	
	// Report
	
	log << endl;
	log << "used Ngrid = ";
	for (int Ngrid : Ngrid_vec) log << Ngrid << " ";
	log << endl;
	log << "free energies = ";
	for (double freeEnergy : freeEnergy_vec) log << freeEnergy << " ";
	log << endl;
	log << "densities = ";
	for (double density : density_vec) log << density << " ";
	log << endl;
	log << "Cvac = ";
	for (double Cvac : Cvac_vec) log << Cvac << " ";
	log << endl;
	log << "alpha = ";
	for (double alpha : alpha_vec) log << alpha << " ";
	log << endl;
	log << "minimum Ngrid = " << Ngrid_min << endl;
	log << "minimum free energy = " << freeEnergy_min << endl;
	log << "minimum density = " << density_min << endl;
	log << "minimum Cvac = " << Cvac_min << endl;
	log << "minimum alpha = " << alpha_min << endl;
	log << "status = " << statusMin << endl;
	
	// Check validity of the result
	
	if (statusMin!=0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: No minimum found for kT=" << kT << " and mu=" << mu << endl;
		log << "       (Unable to perform minimisation with interpolation)" << endl;
		
		return 1;
	}
	
	///// Store the result /////
	
	sskT.str(""); sskT << scientific << setprecision(4) << kT;
	ssMu.str(""); ssMu << scientific << setprecision(4) << mu;
	
	dataFile.open(dataDir+"/kTMuSolid_"+"kT="+sskT.str()+"_mu="
	              +ssMu.str()+".dat");
	
	dataFile << "# Result for solid computation at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << endl;
	dataFile << "Ngrid_min = " << scientific << setprecision(8) << Ngrid_min << endl;
	dataFile << "Cvac_min = " << scientific << setprecision(8) << Cvac_min << endl;
	dataFile << "alpha_min = " << scientific << setprecision(8) << alpha_min << endl;
	dataFile << "freeEnergy_min = " << scientific << setprecision(8) << freeEnergy_min << endl;
	dataFile << "density_min = " << scientific << setprecision(8) << density_min << endl;
	dataFile << endl;
	dataFile << "success = " << true << endl;
	
	dataFile.close();
	
	return 0; // everything was fine
}











////////////////////////////////////////////////////////////////////////////


// Only minimise over the gaussian parameter alpha

int minOverAlphaOnly(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &alpha_min)
{
	///////////////////////////////// log //////////////////////////////////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Minimisation over alpha (only) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	
	// Seperate log file for the DFT computation
	// by default, we don't want the non-interesting log from the DFT to 
	// pollute the current log
	
	Log logDFT("log_DFT.dat");
	
	//////////////////////////// read options //////////////////////////////
	
	string dataDir = "data";
	double kT_rangePrevResult = 0.03;
	double mu_rangePrevResult = 0.15;
	int Ngrid_rangePrevResult = 2;
	double Cvac_min = 0; // default value
	double Cvac_step = 0.0001;
	double Cvac_max = 0.01;
	double alpha_rangeMin;
	double alpha_rangeMax;
	double alpha_logStepMin;
	double alpha_logStepMax;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kT_rangePrevResult", &kT_rangePrevResult);
	options.addOption("mu_rangePrevResult", &mu_rangePrevResult);
	options.addOption("Ngrid_rangePrevResult", &Ngrid_rangePrevResult);
	options.addOption("alpha_rangeMin", &alpha_rangeMin);
	options.addOption("alpha_rangeMax", &alpha_rangeMax);
	options.addOption("alpha_logStepMin", &alpha_logStepMin);
	options.addOption("alpha_logStepMax", &alpha_logStepMax);
	
	options.read(argc, argv);
	
	int sysres = system(("mkdir -p "+dataDir).c_str());
	
	///// Check if computation has not already been done /////
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	stringstream ssNgrid; ssNgrid << Ngrid;
	
	ifstream dataInFile(dataDir+"/minAlphaOnly_"+"kT="+sskT.str()+"_mu="
	                    +ssMu.str()+"_Ngrid="+ssNgrid.str()+".dat");
	
	// if file exists, get the result from there
	if (dataInFile)
	{
		// check if it was a successful computation
		int success;
		readDataFromFile(dataInFile, "success", success);
		
		if (success==1)  // bool "true" 
		{
			log << endl;
			log << "Returning existing computation" << endl;
			
			// read results in file
			readDataFromFile(dataInFile, "freeEnergy_min", freeEnergy_min);
			readDataFromFile(dataInFile, "density_min", density_min);
			readDataFromFile(dataInFile, "alpha_min", alpha_min);
			
			log << endl;
			log << "freeEnergy_min = " << freeEnergy_min << endl;
			log << "density_min = " << density_min << endl;
			log << "alpha_min = " << alpha_min << endl;
			
			return 0;
		}
		
		if (success==0)  // bool "false" 
		{
			log << endl;
			log << "Returning existing computation (is a failed one)" << endl;
			
			return 1;
		}
	}
	
	
	//////////////// By default, save result as a failure //////////////////
	
	sskT.str(""); sskT << scientific << setprecision(4) << kT;
	ssMu.str(""); ssMu << scientific << setprecision(4) << mu;
	ssNgrid.str(""); ssNgrid << Ngrid;
	
	ofstream dataFile(dataDir+"/minAlphaOnly_"+"kT="+sskT.str()+"_mu="+ssMu.str()
	                  +"_Ngrid="+ssNgrid.str()+".dat");
	
	dataFile << "# Result for min over alpha at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << "Ngrid = " << Ngrid << endl;
	dataFile << endl;
	dataFile << "success = " << false << endl;
	
	dataFile.close();
	
	//////////////////////////// Minimisation //////////////////////////////
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Starting minimisation" << endl;
	
	struct CompResult 
	{
		double freeEnergy;
		double density;
		double alpha;
	};
	
	vector<CompResult> results;
	bool foundMinimum = false;
	bool init = true;
	double alpha_current = alpha_rangeMax;
	double alpha_logStep = -alpha_logStepMax;
		
	// loop to repeat calculations for different alpha until we find 
	// the minimum or go outside the limits
	while (!foundMinimum)
	{
		// increment except the first time
		
		if (init) init = false;
		else alpha_current *= exp(alpha_logStep);
		
		// free energy computation
		
		int status_current;
		double freeEnergy_current, density_current;
		double Cvac_current = Cvac_min;
		
		do
		{
			log << "Trying alpha = " << alpha_current 
			    << " and Cvac = " << Cvac_current << endl;
			
			status_current = DFTgaussian( argc, argv, logDFT, kT, mu, Ngrid, 
				Cvac_current, alpha_current, freeEnergy_current, density_current);
			
			Cvac_current += Cvac_step;
		}
		while (status_current!=0 && Cvac_current<Cvac_max);
		
		if (status_current==0)
		{
			// report
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "alpha_current = " << alpha_current << endl;
			log << "freeEnergy_current = " << freeEnergy_current << endl;
			log << "density_current = " << density_current << endl;
			
			// add to the results
			
			struct CompResult res;
			res.freeEnergy = freeEnergy_current;
			res.density = density_current;
			res.alpha = alpha_current;
			results.push_back(res);
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Free energy computation failed for alpha_current = "
			    << alpha_current << endl;
		}
		
		
		// checking if alpha still in the authorised boundaries
		
		if (alpha_current < alpha_rangeMin || alpha_current > alpha_rangeMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Went outside the authorised limits for alpha";
			return 1;
		}
		
		// go back and lower step if we detect an increase in free energy
		// that means we crossed the minimum
		
		if (results.size()>=2)
		{
			double freeEnergy0 = results[results.size()-1].freeEnergy;
			double freeEnergy1 = results[results.size()-2].freeEnergy;
			
			if (freeEnergy0 > freeEnergy1)
			{
				// reverse the search direction and lower the step
				alpha_logStep *= -1;
				alpha_logStep /=  2;
				
				// checking if we reached the max precision
				// if we did, return values at alpha1 as the minimum
				if (abs(alpha_logStep) < alpha_logStepMin)
				{
					alpha_min = alpha_current;
					freeEnergy_min = freeEnergy_current;
					density_min = density_current;
					foundMinimum = true;
				}
			}
		}
	}
	
	
	/////////////////////////////// Finalise ///////////////////////////////
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min alpha) Min alpha = " << alpha_min << endl;
	log << "(Min alpha) Min freeEnergy = " << freeEnergy_min << endl;
	log << "(Min alpha) Min density = " << density_min << endl;
	
	
	////////////////////////// Save result in file /////////////////////////
	
	sskT.str(""); sskT << scientific << setprecision(4) << kT;
	ssMu.str(""); ssMu << scientific << setprecision(4) << mu;
	ssNgrid.str(""); ssNgrid << Ngrid;
	
	dataFile.open(dataDir+"/minAlphaOnly_"+"kT="+sskT.str()+"_mu="+ssMu.str()
	                  +"_Ngrid="+ssNgrid.str()+".dat");
	
	dataFile << "# Result for min over alpha at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << "Ngrid = " << Ngrid << endl;
	dataFile << endl;
	dataFile << "alpha_min = " << scientific << setprecision(8) << alpha_min << endl;
	dataFile << "freeEnergy_min = " << scientific << setprecision(8) << freeEnergy_min << endl;
	dataFile << "density_min = " << scientific << setprecision(8) << density_min << endl;
	dataFile << endl;
	dataFile << "success = " << true << endl;
	
	dataFile.close();
	
	return 0;
}










////////////////////////////////////////////////////////////////////////////

// Structure of parameters for the DFT function for the GSL algorithm

struct DFT_params
{
	double kT, mu;
	int Ngrid;
	int argc;
	char **argv;
};

// Function to minimise using the GSL minimisation algorithm (see below)

double myfun_DFT (const gsl_vector *v, void *params)
{
	double Cvac  = exp( gsl_vector_get(v, 0) );
	double alpha = exp( gsl_vector_get(v, 1) );
	
	struct DFT_params *p = (struct DFT_params *)params;
	double kT = p->kT;
	double mu = p->mu;
	int Ngrid = p->Ngrid;
	int argc = p->argc;
	char **argv = p->argv;
	
	Log logDFT("log_DFT.dat");
	
	double freeEnergy,density;
	int status = DFTgaussian( argc, argv, logDFT,
	             kT, mu, Ngrid, Cvac, alpha, freeEnergy, density);
	
	if (status==0) return freeEnergy;
	//else return GSL_NAN;
	else return 1e8*(1-Cvac); // trick to avoid NAN error and lead the minimiser to larger Cvac
}



// Function to find the minimum of the free energy in the (Cvac,alpha) plane
// using the Nelder-Mead algorithm (from the GSL library)

int minOverCvacAlphaNM_noCatch(double kT, double mu, int Ngrid,
                       int argc, char** argv, Log &log,
                       double &freeEnergy_min, double &density_min,
                       double &Cvac_min, double &alpha_min)
{
	///////////////////////////////// log //////////////////////////////////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Minimisation over alpha and Cvac (Nelder-Mead) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	
	
	//////////////////////////// read options //////////////////////////////
	
	string dataDir = "data";
	
	double Cvac_init  = 1e-3;
	double alpha_init = 300;
	double logStepInit = 0.25;
	double logStepEnd = 0.0025;
	
	double Cvac_rangeMin = 1e-6;
	double Cvac_rangeMax = 1e-1;
	double alpha_rangeMin = 30;
	double alpha_rangeMax = 3000;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	
	options.addOption("Cvac_init", &Cvac_init);
	options.addOption("alpha_init", &alpha_init);
	options.addOption("logStepInit", &logStepInit);
	options.addOption("logStepEnd", &logStepEnd);
	
	options.addOption("Cvac_rangeMin", &Cvac_rangeMin);
	options.addOption("Cvac_rangeMax", &Cvac_rangeMax);
	options.addOption("alpha_rangeMin", &alpha_rangeMin);
	options.addOption("alpha_rangeMax", &alpha_rangeMax);
	
	options.read(argc, argv);
	options.write(log);
	
	int sysres = system(("mkdir -p "+dataDir).c_str());
	
	double logCvac_init  = std::log(Cvac_init);
	double logalpha_init = std::log(alpha_init);
	
	//////////////// By default, save result as a failure //////////////////
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	stringstream ssNgrid; ssNgrid << Ngrid;
	
	ofstream dataFile(dataDir+"/minCvacAlpha_"+"kT="+sskT.str()+"_mu="+ssMu.str()
	                  +"_Ngrid="+ssNgrid.str()+".dat");
	
	dataFile << "# Result for min over (Cvac,alpha) at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << "Ngrid = " << Ngrid << endl;
	dataFile << endl;
	dataFile << "success = " << false << endl;
	
	dataFile.close();
	
	//////////////////////////// minimisation //////////////////////////////
	
	// parameters
	struct DFT_params par;
	par.kT = kT;
	par.mu = mu;
	par.Ngrid = Ngrid;
	par.argc = argc;
	par.argv = argv;
	
	// init variables
	const gsl_multimin_fminimizer_type *T = 
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;
	
	// Starting point 
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, logCvac_init);
	gsl_vector_set (x, 1, logalpha_init);
	
	// Set initial step sizes to 1
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, logStepInit);
	
	// Initialize method and iterate
	minex_func.n = 2;
	minex_func.f = myfun_DFT;
	minex_func.params = &par;
	
	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		
		if (status) break;
		
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, logStepEnd);
		
		if (status == GSL_SUCCESS)
			log << "Converged to a minimum" << endl;
		
		log << "Iteration " << iter
		    << "  Cvac  = " << exp(gsl_vector_get (s->x, 0))
		    << "  alpha = " << exp(gsl_vector_get (s->x, 1))
		    << "  freeEnergy = " << s->fval
		    << "  simplex size = " << size
		    << endl;
		
		// check if exceeded given range
		if (exp(gsl_vector_get (s->x, 0))<Cvac_rangeMin || 
		    exp(gsl_vector_get (s->x, 0))>Cvac_rangeMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Went outside of the limits" << endl;
			break;
		}
		if (exp(gsl_vector_get (s->x, 1))<alpha_rangeMin || 
		    exp(gsl_vector_get (s->x, 1))>alpha_rangeMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Went outside of the limits" << endl;
			break;
		}
	}
	while (status == GSL_CONTINUE && iter < 100);
	
	// save minimum
	Cvac_min  = exp(gsl_vector_get (s->x, 0));
	alpha_min = exp(gsl_vector_get (s->x, 1));
	
	// Deallocate memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	
	if (status!=GSL_SUCCESS)	return 1;
	
	/////////////////////////////// Finalise ///////////////////////////////
	
	// Run at minimum to get the density
	
	Log logDFT("log_DFT.dat");
	int statusMin = DFTgaussian( argc, argv, logDFT,
		kT, mu, Ngrid, Cvac_min, alpha_min, freeEnergy_min, density_min);
	
	if (statusMin!=0)
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Unable to perform the last computation at the minimum" << endl;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min Cvac,alpha) Min Cvac = " << Cvac_min << endl;
	log << "(Min Cvac,alpha) Min alpha = " << alpha_min << endl;
	log << "(Min Cvac,alpha) Min freeEnergy = " << freeEnergy_min << endl;
	log << "(Min Cvac,alpha) Min density = " << density_min << endl;
	
	
	////////////////////////// Save result in file /////////////////////////
	
	sskT.str(""); sskT << scientific << setprecision(4) << kT;
	ssMu.str(""); ssMu << scientific << setprecision(4) << mu;
	ssNgrid.str(""); ssNgrid << Ngrid;
	
	dataFile.open(dataDir+"/minCvacAlpha_"+"kT="+sskT.str()+"_mu="+ssMu.str()
	                  +"_Ngrid="+ssNgrid.str()+".dat");
	
	dataFile << "# Result for min over (Cvac,alpha) at " << endl;
	dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
	dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
	dataFile << "Ngrid = " << Ngrid << endl;
	dataFile << endl;
	dataFile << "Cvac_min = " << scientific << setprecision(8) << Cvac_min << endl;
	dataFile << "alpha_min = " << scientific << setprecision(8) << alpha_min << endl;
	dataFile << "freeEnergy_min = " << scientific << setprecision(8) << freeEnergy_min << endl;
	dataFile << "density_min = " << scientific << setprecision(8) << density_min << endl;
	dataFile << endl;
	dataFile << "success = " << true << endl;
	
	dataFile.close();
	
	return 0;
}







int minOverCvacAlpha(double kT, double mu, int Ngrid,
                       int argc, char** argv, Log &log,
                       double &freeEnergy_min, double &density_min,
                       double &Cvac_min, double &alpha_min)
{
	int status;
	
	try
	{
		status = minOverCvacAlphaNM_noCatch(kT, mu, Ngrid, argc, argv, log,
			freeEnergy_min, density_min, Cvac_min, alpha_min);
	}
	catch (...) {return 1;}
	
	return status;
}






////////////////////////////////////////////////////////////////////////////









int DFTgaussianNoCatch( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density)
{
	////////////////////////////  Parameters ///////////////////////////////
	
	// control
	int nCores = 6;
	
	// geometry
	double dx = 0.1;
	string pointsFile("..//SS31-Mar-2016//ss109.05998");
	string useAnalyticWeights = "false";
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	string potentialType = "LJ";
	
	// density initialisation
	int ncopy = 1;
	int ncopy_range = 2;
	string solidType = "FCC";
	
	////////////////////////////////////////////
	
	Options options;
	
	// control
	options.addOption("nCores", &nCores);
	
	// geometry
	options.addOption("dx", &dx);
	options.addOption("IntegrationPointsFile", &pointsFile);
	options.addOption("UseAnalyticWeights", &useAnalyticWeights);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	options.addOption("PotentialType",  &potentialType);
	
	// density initialisation
	options.addOption("ncopy", &ncopy);
	options.addOption("ncopy_range", &ncopy_range);
	options.addOption("SolidType",  &solidType);
	
	options.read(argc, argv);
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- DFT gaussian computation at fixed (kT,mu,Ngrid,Cvac,alpha) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	log << "Cvac = " << Cvac << endl;
	log << "alpha = " << alpha << endl;
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	options.write(log);
	
	#ifdef USE_OMP    
	omp_set_dynamic(0);
	omp_set_num_threads(nCores);
	
	int fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	log << "omp max threads = " << omp_get_max_threads() << endl;
	#endif
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// Check the potential option
	
	if (potentialType == "LJ") log << "Using the LJ potential" << endl;
	else if (potentialType == "WHDF") log << "Using the WHDF potential" << endl;
	else
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Potential type not identified" << endl;
		log << "       Unable to perform Solid DFT computations" << endl;
		
		return 1;
	}
	
	// Check the weights option
	
	if (useAnalyticWeights == "true") log << "Using analytic evaluation of weights" << endl;
	else if (useAnalyticWeights == "false") log << "Using numeric evaluation of weights" << endl;
	else
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Unrecognised option for UseAnalyticWeights" << endl;
		log << "       Unable to perform Solid DFT computations" << endl;
		
		return 1;
	}
	
	///////////////////////  Initialise the System /////////////////////////
	
	
	int ncopy3= ncopy*ncopy*ncopy;
	double L[3] = {Ngrid*dx,Ngrid*dx,Ngrid*dx};
	double alatt = L[0]/ncopy;
	
	////////////////////////////////////
	// Create potential && effective hsd
	
	LJ potentialLJ(sigma1, eps1, rcut1);
	WHDF potentialWHDF(sigma1, eps1, rcut1);
	
	double hsd1;
	if (potentialType == "LJ")   hsd1 = potentialLJ.getHSD(kT);
	if (potentialType == "WHDF") hsd1 = potentialWHDF.getHSD(kT);
	
	/////////////////////////
	// Create density objects
	
	SolidDensity theDensity1(dx, L, hsd1);
	
	if (solidType=="FCC")
	{
		log << "Using FCC lattice to initialise the Solid" << endl;
		theDensity1.initialiseFCCGaussianFromSC(alpha, alatt, ncopy, 1-Cvac);
	}
	// IT TURNS OUT THE HCP CONSTRUCTION IS INCORRECT (see SolidDensity.h)
	//else if (solidType=="HCP")
	//{
	//	log << "Using HCP lattice to initialise the Solid" << endl;
	//	theDensity1.initialiseHCPGaussianFromSC(alpha, alatt, ncopy, 1-Cvac);
	//}
	else
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Solid type not identified" << endl;
		log << "       Unable to initialise solid density" << endl;
		
		return 1;
	}
	
	log << "Natoms = " << theDensity1.getNumberAtoms() << endl;
	
	////////////////////////////////////////////////////////////////////////
	
	// Create DFT and related objects 
	// Then compute physical properties
	
	// Duplicate code for both analytic and numeric weights
	// (how to do it better?)
	
	
	//////// Analytic weights version ////////
	
	if (useAnalyticWeights == "true")
	{
		// Create DFT and related objects 
		
		FMT_Species_Analytic speciesA(theDensity1,hsd1, mu, Ngrid);
		
		DFT dftA(&speciesA);
		
		RSLT fmtA;
		dftA.addHardCoreContribution(&fmtA);
		
		Interaction_Interpolation_QF iA_LJ(speciesA,speciesA,potentialLJ,kT,log);
		Interaction_Interpolation_QF iA_WHDF(speciesA,speciesA,potentialWHDF,kT,log);
		
		if (potentialType == "LJ") 
		{
			iA_LJ.initialize();
			dftA.addInteraction(&iA_LJ);
		}
		if (potentialType == "WHDF")
		{
			iA_WHDF.initialize();
			dftA.addInteraction(&iA_WHDF);
		}
		
		speciesA.setChemPotential(mu);
		
		// Compute physical properties
		
		freeEnergy = dftA.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
		density = theDensity1.getNumberAtoms()/(L[0]*L[1]*L[2]);
	}
	
	
	
	//////// Numeric weights version ////////
	
	if (useAnalyticWeights == "false")
	{
		// Create DFT and related objects 
		
		FMT_Species_Numeric speciesN(theDensity1,hsd1,pointsFile, mu, Ngrid);
		
		DFT dftN(&speciesN);
		
		RSLT fmtN;
		dftN.addHardCoreContribution(&fmtN);
		
		Interaction iN_LJ(speciesN,speciesN,potentialLJ,kT,log,pointsFile);
		Interaction iN_WHDF(speciesN,speciesN,potentialWHDF,kT,log,pointsFile);
		
		if (potentialType == "LJ") 
		{
			iN_LJ.initialize();
			dftN.addInteraction(&iN_LJ);
		}
		if (potentialType == "WHDF")
		{
			iN_WHDF.initialize();
			dftN.addInteraction(&iN_WHDF);
		}
		
		speciesN.setChemPotential(mu);
		
		// Compute physical properties
		
		freeEnergy = dftN.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
		density = theDensity1.getNumberAtoms()/(L[0]*L[1]*L[2]);
	}
	
	
	
	/////// Report ///////
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(DFT Gaussian profile) Free Energy = " << freeEnergy << endl;
	log << "(DFT Gaussian profile) Density = " << density << endl;
	
	return 0;
}








int DFTgaussian( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density)
{
	try
	{
		int status = DFTgaussianNoCatch( argc, argv, log, kT, mu, Ngrid, 
		                                 Cvac, alpha, freeEnergy, density);
		return status;
	}
	catch (...) 
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl;
		log << "ERROR in DFTgaussian (error type not tested)" << endl;
		
		return 1;
	}
	
	return 1;
}





