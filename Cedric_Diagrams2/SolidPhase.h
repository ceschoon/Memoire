////////////////////////////////////////////////////////////////////////////
//                                                                        //
//                  DFT Computation for the Solid Phase                   //
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
#include "SolidFCC.h"
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
// (kT,mu)            are the control parameters that we will play with 
//                    when we build the phase diagram (thermodynamics 
//                    parameters, they describe the macroscopic state)
// (Ngrid,Cvac,alpha) are parameters that will be fixed for each point 
//                    (kT,mu) in thermodynamic space; we find the value of
//                    these parameters by minimising the free energy in the
//                    3d-space given by (Ngrid,Cvac,alpha).



int fixedkTMuSolid(double kT, double mu,
                   int argc, char** argv, Log &log,
                   double &freeEnergy_min, double &density_min, 
                   double &Ngrid_min, double &Cvac_min, double &alpha_min);


int minOverCvacAlpha2D(double kT, double mu, int Ngrid,
                       int argc, char** argv, Log &log,
                       double &freeEnergy_min, double &density_min,
                       double &Cvac_min, double &alpha_min);


int minOverAlphaOnly(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &alpha_min);

// generic function that points to the selected method for the minimisation
int minOverCvacAlpha(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &Cvac_min, double &alpha_min);


int DFTgaussian( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density);


int surfacePlot( double kT, double mu, int Ngrid,
                 int argc, char** argv, Log &log);



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
	
	string dataDir;
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
	
	
	///// Check if computation has not already been done /////
	
	if (!dataDir.empty())
	{
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
	}
	
	
	///// By default, store the result as a failure /////
	// It will be overwritten in the case of a successful computation
	
	if (!dataDir.empty())
	{
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
		
		ofstream dataFile(dataDir+"/kTMuSolid_"+"kT="+sskT.str()+"_mu="
		                  +ssMu.str()+".dat");
		
		dataFile << "# Result for solid computation at " << endl;
		dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
		dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
		dataFile << endl;
		dataFile << "success = " << false << endl;
	}
	
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
	
	if (!dataDir.empty())
	{
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
	
	int statusMin = minFromDataParabola(Ngrid_vec, freeEnergy_vec, 
		Ngrid_min, freeEnergy_min);
	
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
	
	if (!dataDir.empty())
	{
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
		
		ofstream dataFile(dataDir+"/kTMuSolid_"+"kT="+sskT.str()+"_mu="
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
	}
	
	return 0; // everything was fine
}










////////////////////////////////////////////////////////////////////////////


// Function to find minimum of free energy in (Cvac,alpha)-space.
// The plan is:
// o read options
//    get ranges for Cvac and alpha
// o try to use prev result for (kT,mu,Ngrid) that is close
// o calc min (Cvac,alpha)
// o store result somewhere


int minOverCvacAlpha2D(double kT, double mu, int Ngrid,
                       int argc, char** argv, Log &log,
                       double &freeEnergy_min, double &density_min,
                       double &Cvac_min, double &alpha_min)
{
	///////////////////////////////// log //////////////////////////////////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Minimisation over alpha and Cvac ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	
	// Seperate log file for the DFT computation
	// by default, we don't want the non-interesting log from the DFT to 
	// pollute the current log
	
	Log logDFT("log_DFT.dat");
	
	//////////////////////////// read options //////////////////////////////
	
	string dataDir;
	double kT_rangePrevResult = 0.03;
	double mu_rangePrevResult = 0.15;
	int Ngrid_rangePrevResult = 2;
	double Cvac_rangeMin;
	double Cvac_rangeMax;
	double Cvac_logStepMin;
	double Cvac_logStepMax;
	double alpha_rangeMin;
	double alpha_rangeMax;
	double alpha_logStepMin;
	double alpha_logStepMax;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("kT_rangePrevResult", &kT_rangePrevResult);
	options.addOption("mu_rangePrevResult", &mu_rangePrevResult);
	options.addOption("Ngrid_rangePrevResult", &Ngrid_rangePrevResult);
	options.addOption("Cvac_rangeMin", &Cvac_rangeMin);
	options.addOption("Cvac_rangeMax", &Cvac_rangeMax);
	options.addOption("Cvac_logStepMin", &Cvac_logStepMin);
	options.addOption("Cvac_logStepMax", &Cvac_logStepMax);
	options.addOption("alpha_rangeMin", &alpha_rangeMin);
	options.addOption("alpha_rangeMax", &alpha_rangeMax);
	options.addOption("alpha_logStepMin", &alpha_logStepMin);
	options.addOption("alpha_logStepMax", &alpha_logStepMax);
	
	options.read(argc, argv);
	
	
	///// Check if computation has not already been done /////
	
	if (!dataDir.empty())
	{
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
		stringstream ssNgrid; ssNgrid << Ngrid;
		
		ifstream dataInFile(dataDir+"/minCvacAlpha_"+"kT="+sskT.str()+"_mu="
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
				readDataFromFile(dataInFile, "Cvac_min", Cvac_min);
				readDataFromFile(dataInFile, "alpha_min", alpha_min);
				
				log << endl;
				log << "freeEnergy_min = " << freeEnergy_min << endl;
				log << "density_min = " << density_min << endl;
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
	}
	
	
	//////////////// By default, save result as a failure //////////////////
	
	if (!dataDir.empty())
	{
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
	}
	
	///////////////////////// init minimisation ////////////////////////////
	
	// minimisation variables
	
	bool foundMinimum = false;
	
	double alpha_current;
	double Cvac_current;
	double freeEnergy_current;
	double density_current;
	
	double alpha_logStep = alpha_logStepMax;
	double Cvac_logStep = Cvac_logStepMax;
	
	// max number of dft computations allowed in this function 
	
	int expectedNumComp = int (std::log(alpha_rangeMax/alpha_rangeMin)/alpha_logStepMax)+
	                    + int (std::log(Cvac_rangeMax/Cvac_rangeMin)/Cvac_logStepMax);
	int maxDFTcomputations = 8*expectedNumComp;
	int numDFTcomputations = 0;
	
	// Random number engine to choose starting position in (Cvac,alpha)
	
	int seed = chrono::system_clock::now().time_since_epoch().count();
	log << "seed = " << seed << endl;
	
	default_random_engine gen(seed);
	uniform_real_distribution<double> dist01(0,1);
	
	bool initialisePosition = true;
	
	
	/////// Initialise (Cvac,alpha) from previous result if possible ///////
	
	if (!dataDir.empty())
	{
		log << endl;
		log << "Trying to initialise position in (Cvac,alpha) from previous results" << endl;
		
		DIR *dir;
		struct dirent *ent;

		if ((dir = opendir (dataDir.c_str())) != NULL) 
		{
			while ((ent = readdir (dir)) != NULL) 
			{
				string ent_name(ent->d_name);
				
				if (ent_name.size()>=12 && ent_name.substr(0,12) == "minCvacAlpha")
				{
					log << "Checking file: " << ent_name << endl;
					
					// get kT, mu and Ngrid from file name
					
					int pos1_=-1; // position of underscore in name (_kT=)
					int pos2_=-1; // position of underscore in name (_mu=)
					int pos3_=-1; // position of underscore in name (_Ngrid=)
					
					for (int i=0; i<ent_name.size(); i++)
					{
						if (ent_name[i]=='_' && pos1_<0) pos1_ = i;
						else if (ent_name[i]=='_' && pos2_<0) pos2_ = i;
						else if (ent_name[i]=='_' && pos3_<0) pos3_ = i;
					}
					
					int lenkTstr = pos2_-pos1_-4;
					int lenmustr = pos3_-pos2_-4;
					int lenNgridstr = ent_name.size()-4-pos3_-7;
					
					string kT_str = ent_name.substr(pos1_+4,lenkTstr);
					string mu_str = ent_name.substr(pos2_+4,lenmustr);
					string Ngrid_str = ent_name.substr(pos3_+7,lenNgridstr);
					
					double kT_read = stod(kT_str);
					double mu_read = stod(mu_str);
					int Ngrid_read = stoi(Ngrid_str);
					
					// check if (kT,mu,Ngrid) in file are close enough to be used
					
					if ( (abs(kT-kT_read)<kT_rangePrevResult) ||
					     (abs(mu-mu_read)<mu_rangePrevResult) ||
					     (abs(Ngrid-Ngrid_read)<=Ngrid_rangePrevResult)   )
					{
						ifstream dataInFile(dataDir+"/"+ent_name);
						if (dataInFile)
						{
							int success_read;
							double Cvac_read;
							double alpha_read;
							
							readDataFromFile(dataInFile, "success", success_read);
							
							if (success_read==1) // here 1 means true, thus success
							{
								readDataFromFile(dataInFile, "Cvac_min", Cvac_read);
								readDataFromFile(dataInFile, "alpha_min", alpha_read);
								
								log << "Trying initial condition alpha = " << alpha_read
								    << " and Cvac = " << Cvac_read << endl;
								
								double freeEnergy_temp, density_temp;
								int success_comp = 
									DFTgaussian( argc, argv, logDFT, 
									kT, mu, Ngrid, Cvac_read, 
									alpha_read, freeEnergy_temp, density_temp);
								
								if (success_comp==0)
								{
									Cvac_current = Cvac_read;
									alpha_current = alpha_read;
									freeEnergy_current = freeEnergy_temp;
									density_current = density_temp;
									
									log << "Initialise from file to Cvac=" 
										<< Cvac_read << " and alpha=" 
										<< alpha_read << endl;
									
									initialisePosition = false;
									break;
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
			
			closedir (dir);
		} 
		else 
		{
			// could not open directory 
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: could not open data directory (minOverCvacAlpha2D)" << endl;
		}
	}
	
	
	//////////////////////////// Minimisation //////////////////////////////
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "Starting minimisation" << endl;
	
	int indexPrev = 0; // previous direction with lower free energy
	
	while (!foundMinimum && numDFTcomputations<=maxDFTcomputations)
	{
		// Starting position in (Cvac,alpha) plane initialised randomly 
		
		while (initialisePosition)
		{
			double ran1 = dist01(gen);
			double ran2 = dist01(gen);
			
			alpha_current = exp(   std::log(alpha_rangeMin)*(1-ran1)
			                     + std::log(alpha_rangeMax)*ran1      );
			Cvac_current  = exp(   std::log(Cvac_rangeMin )*(1-ran1)
			                     + std::log(Cvac_rangeMax )*ran1      );
			
			log << "Trying initial condition alpha = " << alpha_current
				<< " and Cvac = " << Cvac_current << endl;
			
			// initial free energy
			int success_comp = DFTgaussian( argc, argv, logDFT, kT, mu, Ngrid, 
				Cvac_current, alpha_current, freeEnergy_current, density_current);
			
			numDFTcomputations ++;
			
			if (success_comp==0)
			{
				initialisePosition = false;
				log << "Found successful initial condition" << endl;
			}
		}
		
		// compute the free energy in every direction
		
		vector<double> alpha_candidates(4);
		vector<double> Cvac_candidates(4);
		vector<double> freeEnergy_candidates(4);
		vector<double> density_candidates(4);
		vector<bool> valid_candidates(4);
		
		alpha_candidates[0] = alpha_current*exp(alpha_logStep);
		Cvac_candidates[0] = Cvac_current;
		alpha_candidates[1] = alpha_current;
		Cvac_candidates[1] = Cvac_current*exp(Cvac_logStep);
		alpha_candidates[2] = alpha_current*exp(-alpha_logStep);
		Cvac_candidates[2] = Cvac_current;
		alpha_candidates[3] = alpha_current;
		Cvac_candidates[3] = Cvac_current*exp(-Cvac_logStep);
		
		int indexMin = -1; // -1 means current position
		alpha_min = alpha_current;
		Cvac_min = Cvac_current;
		freeEnergy_min = freeEnergy_current;
		density_min = density_current;
		
		for (int i=0; i<4; i++)
		{
			// begin searching a lower free energy in the previous working
			// direction
			
			int j = indexPrev + i;
			while (j>=4) j-=4;
			
			double freeEnergy_temp, density_temp;
			int status_comp = DFTgaussian( argc, argv, logDFT, kT, mu, Ngrid, 
				Cvac_candidates[j], alpha_candidates[j], freeEnergy_temp, density_temp);
			
			numDFTcomputations ++;
			
			log << endl;
			log << "Testing candidate number" << j << endl;
			log << "  alpha = " << alpha_candidates[j] << endl;
			log << "  Cvac = " << Cvac_candidates[j] << endl;
			log << "  freeEnergy = " << freeEnergy_temp << endl;
			log << "  status = " << status_comp << endl;
			
			if (status_comp==0)
			{
				freeEnergy_candidates[j] = freeEnergy_temp;
				density_candidates[j] = density_temp;
				valid_candidates[j] = true;
				
				if(freeEnergy_candidates[j] < freeEnergy_min)
				{
					indexMin = j;
					alpha_min = alpha_candidates[j];
					Cvac_min = Cvac_candidates[j];
					freeEnergy_min = freeEnergy_candidates[j];
					density_min = density_candidates[j];
					indexPrev = j;
					break; // for speed, go in the first direction better than the current position
				}
			}
			else
			{
				valid_candidates[j] = false;
			}
		}
		
		// Report
		
		log << endl;
		log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		
		log << "alpha_current = " << alpha_current << endl;
		log << "Cvac_current = " << Cvac_current << endl;
		log << "freeEnergy_current = " << freeEnergy_current << endl;
		log << "density_current = " << density_current << endl;
		
		// Check if current position is inside of the limits
		// If not, restart from somewhere else
		
		if (alpha_current<alpha_rangeMin || alpha_current>alpha_rangeMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Went outside of the limits, retry with a new starting point" << endl;
			
			// try with a smaller step (less likely to make the same mistake)
			if (alpha_logStep>alpha_logStepMin) 
			{
				alpha_logStep = alpha_logStep/2;
				log << "dividing alpha_logStep by 2" << endl;
			}
			
			initialisePosition = true;
			continue;
		}
		if (Cvac_current<Cvac_rangeMin || Cvac_current>Cvac_rangeMax)
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Went outside of the limits, retry with a new starting point" << endl;
			
			// try with a smaller step (less likely to make the same mistake)
			if (Cvac_logStep>Cvac_logStepMin)
			{
				Cvac_logStep  = Cvac_logStep/2;
				log << "dividing Cvac_logStep by 2" << endl;
			}
			
			initialisePosition = true;
			continue;
		}
		
		// If min freeEnergy is the current position, stay there
		// and try again with a smaller step
		
		if (indexMin==-1)
		{
			if (alpha_logStep>alpha_logStepMin) 
			{
				alpha_logStep = alpha_logStep/4;
				log << "dividing alpha_logStep by 4" << endl;
			}
			if (Cvac_logStep>Cvac_logStepMin)
			{
				Cvac_logStep  = Cvac_logStep/4;
				log << "dividing Cvac_logStep by 4" << endl;
			}
			
			// If we reach the min steps allowed, we return the current 
			// location as the minimum
			
			if (alpha_logStep<alpha_logStepMin && Cvac_logStep<Cvac_logStepMin) 
				foundMinimum = true;
		}
		
		// Otherwise, move to the position with lowest freeEnergy
		
		else
		{
			alpha_current = alpha_candidates[indexMin];
			Cvac_current = Cvac_candidates[indexMin];
			freeEnergy_current = freeEnergy_candidates[indexMin];
			density_current = density_candidates[indexMin];
		}
	}
	
	
	/////////////////////////////// Finalise ///////////////////////////////
	
	// Exit if failure
	
	if (numDFTcomputations>maxDFTcomputations)
	{
		log << "ERROR: Exceeded max number of dft computations allowed" << endl;
		log << "maxDFTcomputations = " << maxDFTcomputations << endl;
		return 1;
	}
	else
	{
		log << "numDFTcomputations = " << numDFTcomputations << endl;
		log << "maxDFTcomputations = " << maxDFTcomputations << endl;
	}
	
	// Report
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(Min Cvac,alpha) Min Cvac = " << Cvac_min << endl;
	log << "(Min Cvac,alpha) Min alpha = " << alpha_min << endl;
	log << "(Min Cvac,alpha) Min freeEnergy = " << freeEnergy_min << endl;
	log << "(Min Cvac,alpha) Min density = " << density_min << endl;
	
	
	////////////////////////// Save result in file /////////////////////////
	
	if (!dataDir.empty())
	{
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
		dataFile << "Cvac_min = " << scientific << setprecision(8) << Cvac_min << endl;
		dataFile << "alpha_min = " << scientific << setprecision(8) << alpha_min << endl;
		dataFile << "freeEnergy_min = " << scientific << setprecision(8) << freeEnergy_min << endl;
		dataFile << "density_min = " << scientific << setprecision(8) << density_min << endl;
		dataFile << endl;
		dataFile << "success = " << true << endl;
	}
	
	return 0;
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
	
	string dataDir;
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
	
	///// Check if computation has not already been done /////
	
	if (!dataDir.empty())
	{
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
	}
	
	
	//////////////// By default, save result as a failure //////////////////
	
	if (!dataDir.empty())
	{
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
		stringstream ssNgrid; ssNgrid << Ngrid;
		
		ofstream dataFile(dataDir+"/minAlphaOnly_"+"kT="+sskT.str()+"_mu="+ssMu.str()
		                  +"_Ngrid="+ssNgrid.str()+".dat");
		
		dataFile << "# Result for min over alpha at " << endl;
		dataFile << "kT = " << scientific << setprecision(4) << kT << endl;
		dataFile << "mu = " << scientific << setprecision(4) << mu << endl;
		dataFile << "Ngrid = " << Ngrid << endl;
		dataFile << endl;
		dataFile << "success = " << false << endl;
	}
	
	
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
	
	if (!dataDir.empty())
	{
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
		stringstream ssNgrid; ssNgrid << Ngrid;
		
		ofstream dataFile(dataDir+"/minAlphaOnly_"+"kT="+sskT.str()+"_mu="+ssMu.str()
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
	}
	
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
	
	string dataDir;
	
	double Cvac_init  = 1e-3;
	double alpha_init = 300;
	double logStepInit = 0.25;
	double logStepEnd = 0.0025;
	
	double Cvac_rangeMin = 1e-6;
	double Cvac_rangeMax = 1e-1;
	double alpha_rangeMin = 30;
	double alpha_rangeMax = 3000;
	
	Options options;
	
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
	
	double logCvac_init  = std::log(Cvac_init);
	double logalpha_init = std::log(alpha_init);
	
	//////////////// By default, save result as a failure //////////////////
	
	if (!dataDir.empty())
	{
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
	}
	
	
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
	
	if (!dataDir.empty())
	{
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
		dataFile << "Cvac_min = " << scientific << setprecision(8) << Cvac_min << endl;
		dataFile << "alpha_min = " << scientific << setprecision(8) << alpha_min << endl;
		dataFile << "freeEnergy_min = " << scientific << setprecision(8) << freeEnergy_min << endl;
		dataFile << "density_min = " << scientific << setprecision(8) << density_min << endl;
		dataFile << endl;
		dataFile << "success = " << true << endl;
	}
	
	return 0;
}


int minOverCvacAlphaNM(double kT, double mu, int Ngrid,
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



int minOverCvacAlpha(double kT, double mu, int Ngrid,
                     int argc, char** argv, Log &log,
                     double &freeEnergy_min, double &density_min,
                     double &Cvac_min, double &alpha_min)
{
	#ifdef MIN_OVER_ALPHA_ONLY
	return minOverAlphaOnly(kT,mu,Ngrid,argc,argv,log,
		freeEnergy_min,density_min,alpha_min);
	
	#else
	return minOverCvacAlphaNM(kT,mu,Ngrid,argc,argv,log,
		freeEnergy_min,density_min,Cvac_min,alpha_min);
	//return minOverCvacAlpha2D(kT,mu,Ngrid,argc,argv,log,
	//	freeEnergy_min,density_min,Cvac_min,alpha_min);
	
	#endif
}



////////////////////////////////////////////////////////////////////////////






void DFTgaussianNoCatch( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density)
{
	////////////////////////////  Parameters ///////////////////////////////
	
	// control
	int nCores = 6;
	string pointsFile("..//SS31-Mar-2016//ss109.05998");
	
	// geometry
	double dx = 0.1;
	
	// potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	double hsd1 = -1;
	
	// density initialisation
	int ncopy = 1;
	int ncopy_range = 2;
	
	////////////////////////////////////////////
	
	Options options;
	
	// control
	options.addOption("nCores", &nCores);
	options.addOption("IntegrationPointsFile", &pointsFile);
	
	// geometry
	options.addOption("dx", &dx);
	
	// potential
	options.addOption("eps1",   &eps1);
	options.addOption("sigma1", &sigma1);
	options.addOption("rcut1",  &rcut1);
	
	// density initialisation
	options.addOption("ncopy", &ncopy);
	options.addOption("ncopy_range", &ncopy_range);
	
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
	
	////////////////////////////////////
	
	int ncopy3= ncopy*ncopy*ncopy;
	double L[3] = {Ngrid*dx,Ngrid*dx,Ngrid*dx};
	double alatt = L[0]/ncopy;
	
	
	///////////////////////  Initialise the System /////////////////////////
	
	//////////////////////////////////////
	////// Create potential && effective hsd
	
	#ifdef POTENTIAL_WHDF
	log << "Using WHDF Potential" << endl;
	WHDF potential1(sigma1, eps1, rcut1);
	#else
	log << "Using LJ Potential" << endl;
	LJ potential1(sigma1, eps1, rcut1); // default
	#endif
	
	if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
	
	/////////////////////////////////////
	// Create density objects
	
	SolidFCC theDensity1(dx, L, hsd1);
	
	/////////////////////////////////////
	// Create the species objects
	
	#ifdef ANALYTIC_WEIGHTS
	log << "Using analytic evaluation of weights" << endl;
	FMT_Species_Analytic species1(theDensity1,hsd1, mu, Ngrid);
	Interaction_Linear_Interpolation i1(species1,species1,potential1,kT,log);
	i1.initialize();
	
	#else
	log << "Using numeric evaluation of weights" << endl;
	FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile, mu, Ngrid);
	Interaction i1(species1,species1,potential1,kT,log,pointsFile);
	i1.initialize();
	#endif
	
	/////////////////////////////////////
	// Create the hard-sphere object
	RSLT fmt;
	
	/////////////////////////////////////
	// DFT object
	
	DFT dft(&species1);
	
	dft.addHardCoreContribution(&fmt);  
	dft.addInteraction(&i1);
	
	/////////////////////////////////////////////////////
	// Report
	log << "Lx = " << L[0] << endl;
	log << "HSD 1 = " << hsd1 << endl;
	
	
	///////////////////////////////////////////////
	// Fix the chemical potential of the surfactant species.
	
	species1.setChemPotential(mu);
	
	//  check(theDensity1, dft,i1);
	
	log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
	
	
	////////////////////// Run with gaussian profiles //////////////////////
	
	theDensity1.initializeGaussian(alpha, alatt, ncopy, 1-Cvac);
	freeEnergy = dft.calculateFreeEnergyAndDerivatives(false)/(L[0]*L[1]*L[2]);
	density = theDensity1.getNumberAtoms()/(L[0]*L[1]*L[2]);
	
	// Report
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "(DFT Gaussian profile) Free Energy = " << freeEnergy << endl;
	log << "(DFT Gaussian profile) Density = " << density << endl;
}


int DFTgaussian( int argc, char** argv, Log &log,
                 double kT, double mu, int Ngrid, double Cvac, double alpha,
                 double &freeEnergy, double &density)
{
	try
	{
		DFTgaussianNoCatch( argc, argv, log, kT, mu, Ngrid, Cvac, alpha, 
		                    freeEnergy, density);
	}
	catch (...) //any
	{
		log << myColor::RED << "=================================" << myColor::RESET << endl;
		//log << "ERROR: Eta too large in DFTgaussian" << endl;
		log << "ERROR in DFTgaussian (error type not identified)" << endl;
		
		return 1;
	}
	
	return 0;
}












////////////////////////////////////////////////////////////////////////////

// Function to compute the freeEnergy on the (Cvac,alpha) plane in order
// to plot a heatmap.

int surfacePlot( double kT, double mu, int Ngrid,
                 int argc, char** argv, Log &log)
{
	////////////////////////////////// Log /////////////////////////////////
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	log << myColor::RED << myColor::BOLD << "--- Surface plot on (Cvac,alpha) ---" << myColor::RESET << endl <<  "#" << endl;
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	
	log << "kT = " << kT << endl;
	log << "mu = " << mu << endl;
	log << "Ngrid = " << Ngrid << endl;
	
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// Seperate log file for the DFT computation
	// by default, we don't want the non-interesting log from the DFT to 
	// pollute the current log
	
	Log logDFT("log_DFT.dat");
	
	///////////////////////////// Read options /////////////////////////////
	
	string dataDir = "surface_plots";
	double Cvac_rangeMin;
	double Cvac_rangeMax;
	double Cvac_logStep;
	double alpha_rangeMin;
	double alpha_rangeMax;
	double alpha_logStep;
	
	Options options;
	
	options.addOption("DataDirectory", &dataDir);
	options.addOption("Cvac_rangeMin", &Cvac_rangeMin);
	options.addOption("Cvac_rangeMax", &Cvac_rangeMax);
	options.addOption("Cvac_logStep", &Cvac_logStep);
	options.addOption("alpha_rangeMin", &alpha_rangeMin);
	options.addOption("alpha_rangeMax", &alpha_rangeMax);
	options.addOption("alpha_logStep", &alpha_logStep);
	
	options.read(argc, argv);
	
	/////////////////////////////// Init file //////////////////////////////
	
	int sysresult = system(("mkdir -p "+dataDir).c_str());
	
	stringstream sskT; sskT << scientific << setprecision(4) << kT;
	stringstream ssMu; ssMu << scientific << setprecision(4) << mu;
	stringstream ssNgrid; ssNgrid << Ngrid;
	
	string surfacePlotFileName = "";
	surfacePlotFileName = surfacePlotFileName + "surfacePlot_"
	                    + "kT="+sskT.str() + "_mu="+ssMu.str()
	                    + "_Ngrid="+ssNgrid.str();
	ofstream surfacePlotFile(dataDir+"/"+surfacePlotFileName+".dat");
	
	surfacePlotFile << "# Parameters:" << endl;
	surfacePlotFile << "# kT = " << kT << endl;
	surfacePlotFile << "# mu = " << mu << endl;
	surfacePlotFile << "# Ngrid = " << Ngrid << endl;
	surfacePlotFile << "# " << endl;
	surfacePlotFile << "# Plot:" << endl;
	surfacePlotFile << "# logCvac		logalpha		freeEnergy	"
	                << endl;
	surfacePlotFile << "# " << endl;
	surfacePlotFile << scientific;
	
	//////////////// Compute free energy on (Cvac,alpha) plane /////////////
	
	for (double Cvac=Cvac_rangeMin; Cvac<Cvac_rangeMax; Cvac*=exp(Cvac_logStep))
	{
		for (double alpha=alpha_rangeMin; alpha<alpha_rangeMax; alpha*=exp(alpha_logStep))
		{
			log << endl;
			log << "Computing for Cvac = " << Cvac << " and alpha = " 
			    << alpha << endl;
			
			double freeEnergy, density;
			int statusDFT = DFTgaussian( argc, argv, logDFT, kT, mu, Ngrid, 
				Cvac, alpha, freeEnergy, density);
			
			if (statusDFT==0)
			{
				log << "freeEnergy = " << freeEnergy << endl;
				surfacePlotFile << std::log(Cvac)  << " 	"
				                << std::log(alpha) << " 	"
				                << freeEnergy << endl;
			}
			else
			{
				log << "freeEnergy computation failed" << endl;
				surfacePlotFile << std::log(Cvac)  << " 	"
				                << std::log(alpha) << " 	"
				                << "failed" << endl;
			}
		}
		surfacePlotFile << endl; // seperate data blocks
	}
	
	/////////////////////////// Plot with gnuplot //////////////////////////
	
	ofstream plotFile(dataDir+"/plot_surface");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output '" << surfacePlotFileName << ".svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Solid Free Energy Surface\" font \",20\"" << endl;
	plotFile << "set label \"kT = " << kT << " mu = " << mu << " Ngrid = " << Ngrid << "\" at graph 0.03,1.03 font \",16\"" << endl;
	plotFile << "set xlabel \"log(Cvac)\" font \",20\"" << endl;
	plotFile << "set ylabel \"log(alpha)\" font \",20\"" << endl;
	//plotFile << "set zlabel \"Omega/(kT*V)\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "set datafile missing \"failed\"" << endl;
	//plotFile << "set view map" << endl;
	//plotFile << "set dgrid3d" << endl;
	//plotFile << "set pm3d interpolate 10,10" << endl;
	//plotFile << "splot \"" << surfacePlotFileName << ".dat\" using 1:2:3 with pm3d" << endl;
	plotFile << "plot \"" << surfacePlotFileName << ".dat\" using 1:2:($3/1.0) with image" << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+dataDir+"; gnuplot plot_surface").c_str());
	
	return 0;
}



