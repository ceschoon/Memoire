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



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Estimation of the coexistence mu from kT
// Used to limit the range of mu for the dft computations
// It is assumed that the temperatures are sorted from the smallest to the
// largest

// This is a linear interpolation of coexistence data from previous calculations

void estimateMuCoexFromkT(vector<double> kTEsti, vector<double> muEsti, 
                          double kT, double &mu, bool &success)
{
	double small_value = 1e-10;
	
	success = true;
	if (kT<kTEsti[0]-small_value || kT>kTEsti[kTEsti.size()-1]+small_value) 
		success = false;
	
	for (int i=0; i<kTEsti.size()-1; i++)
	{
		if (kT>kTEsti[i]-small_value && kT<kTEsti[i+1]+small_value)
		{
			double x = (kT-kTEsti[i])/(kTEsti[i+1]-kTEsti[i]);
			mu = muEsti[i]*(1-x)+muEsti[i+1]*x;
		}
	}
}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to perform dft computations for a given range of the parameters
// (kT,mu,Npoints).

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kTMin = 1;
	double kTMax = 1;
	double kTStep = 0.1;
	int kTStartIndex = 0;
	double muMin = -3;
	double muMax = 1;
	double muStep = 0.5;
	int muStartIndex = 0;
	int NpointsMin = 122;
	int NpointsMax = 135;
	int NpointsStep = 1;
	int NpointsStartIndex = 0;
	
	string muCoexEstimationFileName;
	double muCoexEstimationRange = 2.0;
	
	Options options;
	options.addOption("kTMin", &kTMin);
	options.addOption("kTMax", &kTMax);
	options.addOption("kTStep", &kTStep);
	options.addOption("kTStartIndex", &kTStartIndex);
	options.addOption("MuMin", &muMin);
	options.addOption("MuMax", &muMax);
	options.addOption("MuStep", &muStep);
	options.addOption("MuStartIndex", &muStartIndex);
	options.addOption("NpointsMin", &NpointsMin);
	options.addOption("NpointsMax", &NpointsMax);
	options.addOption("NpointsStep", &NpointsStep);
	options.addOption("NpointsStartIndex", &NpointsStartIndex);
	options.addOption("MuCoexEstimationFile", &muCoexEstimationFileName);
	options.addOption("MuCoexEstimationRange", &muCoexEstimationRange);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	
	
	
	///////////////////////// Prepare computations /////////////////////////
	
	string dataDir1 = "Data-kt-mu-Npoints-Solid";
	
	int sysresult = system(("mkdir " + dataDir1).c_str());
	
	int kTSize = int( (kTMax-kTMin)/kTStep +1 );
	vector<double> kT(kTSize);
	for (int i=0; i<kTSize; i++) kT[i] = kTMin + kTStep*i;
	
	int muSize = int( (muMax-muMin)/muStep +1 );
	vector<double> mu(muSize);
	for (int i=0; i<muSize; i++) mu[i] = muMin + muStep*i;
	
	int NpointsSize = int( (NpointsMax-NpointsMin)/NpointsStep +1 );
	vector<int> Npoints(NpointsSize);
	for (int i=0; i<NpointsSize; i++) Npoints[i] = NpointsMin + NpointsStep*i;
	
	
	
	/////////////////// Estimation of the coexistence mu ///////////////////
	
	vector<double> kTEsti;
	vector<double> muEsti;
	
	// Read data file
	
	if (!muCoexEstimationFileName.empty())
	{
		ifstream dataFile(muCoexEstimationFileName);
		
		if (dataFile)
		{
			readColumnVectorFromFile(dataFile, 3, kTEsti);
			readColumnVectorFromFile(dataFile, 2, muEsti);
			
			// Report 
			
			log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "kTEsti = ";
			for (int i=0; i<kTEsti.size(); i++)
				log << kTEsti[i] << " ";
			log << endl;
			log << "muEsti = ";
			for (int i=0; i<muEsti.size(); i++)
				log << muEsti[i] << " ";
			log << endl;
		}
		else
		{
			log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Unable to open coexistence estimation file" << endl;
		}
	}
	
	//////////////////////// Run DFT computations //////////////////////////
	
	// TODO: parallelise loop ?
	for (int i=0; i<NpointsSize; i++)
	{
		bool useSnapshot = false;
		
		for (int j=0; j<kTSize; j++)
		{
			// limit the range of mu with an estimation of where the 
			// coexistence is
			
			double muCoexEstimated;
			bool successEstimation = false;
			
			if (!muCoexEstimationFileName.empty())
			{
				estimateMuCoexFromkT(kTEsti,muEsti,kT[j],muCoexEstimated,
				                     successEstimation);
				if (!successEstimation) 
				{
					log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
					log << "ERROR: Estimation of mu not successful for kT = " << kT[j] << endl;
				}
				else
				{
					log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
					log << "Estimated coexistence mu = " << muCoexEstimated
					    << " for kT = " << kT[j] << endl;
				}
			}
			
			for (int k=0; k<muSize; k++)
			{
				// option to limit range with estimation of the coexistence mu
				if (!muCoexEstimationFileName.empty() && successEstimation)
					if (mu[k]<muCoexEstimated-muCoexEstimationRange/2 || 
					    mu[k]>muCoexEstimated+muCoexEstimationRange/2) 
						continue; // skip this mu
				
				// dft computation
				
				double densitySolid = 0;
				double numParticlesSolid = 0;
				double freeEnergySolid = 0;
				double aVdW = 0;
				double hsd = 0;
				bool successSolid = false;
				
				try
				{
					DFTcomputation( argc, argv, log, kT[j], mu[k], Npoints[i], 
					                densitySolid, numParticlesSolid, 
					                freeEnergySolid, aVdW, hsd, 
					                useSnapshot, successSolid);
				}
				catch (...)
				{
					log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
					log << "ERROR: DFT computation failed for "
					    << "kT = " << kT[j] << " " 
					    << "mu = " << mu[k] << " " 
					    << "Npoints = " << Npoints[i] << " " 
					    << endl;
					log << myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
					
					successSolid = false;
				}
				
				if (successSolid) useSnapshot = true;
				else useSnapshot = false;
				
				// save in file
				
				ofstream dataFile(dataDir1+"/results_"+to_string(NpointsStartIndex+i)+
				                                   "_"+to_string(kTStartIndex+j)+
				                                   "_"+to_string(muStartIndex+k)+".dat");
				
				dataFile << "Npoints = " << Npoints[i] << endl;
				dataFile << "kT = " << kT[j] << endl;
				dataFile << "mu = " << mu[k] << endl;
				dataFile << endl;
				dataFile << "aVdW = " << aVdW << endl;
				dataFile << "hsd  = " << hsd  << endl;
				dataFile << endl;
				dataFile << "densitySolid = " << densitySolid << endl;
				dataFile << "numParticlesSolid = " << numParticlesSolid << endl;
				dataFile << "freeEnergySolid = " << freeEnergySolid << endl;
				dataFile << endl;
				dataFile << "successSolid = " << successSolid << endl;
				
				dataFile.close();
				
			}
			
			// clean weights
			
			//int sysresult = system("rm weights*");
		}
	}
	
	return 1;
}

