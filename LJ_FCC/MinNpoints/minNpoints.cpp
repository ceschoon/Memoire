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

#include "../utilities.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program thats performs minimisation over Npoints of the DFT data

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double maxFileIndex = 100;
	double densityLimitSolid = 0.8; // lowest solid density tolerated
	                                // used to discard computations that
	                                // converged to a fluid phase
	
	// options 
	
	Options options;
	options.addOption("MaxFileIndex", &maxFileIndex);
	options.addOption("DensityLimitSolid", &densityLimitSolid);
	options.read(argc, argv);
	
	// log
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// files and directories
	
	string dataDir1 = "Data-kt-mu-Npoints-Solid";
	string dataDir2 = "Data-kt-mu-Solid";
	string plotDir  = "Plots-kt-mu-Solid"; 
	
	int sysresult = system(("mkdir " + dataDir2).c_str());
	    sysresult = system(("mkdir " + plotDir ).c_str());
	
	ofstream plotFileErase(plotDir+"/plot");
	plotFileErase << endl;
	plotFileErase.close();
	
	////////// Loop over all possible files in Solid directory  ////////////
	
	// not optimal, but fast enough
	
	for (int kTCounter=0; kTCounter<maxFileIndex; kTCounter++)
	{
		for (int muCounter=0; muCounter<maxFileIndex; muCounter++)
		{
			bool existsFileForThatkTmu = false;
			
			double kT = 0;
			double mu = 0;
			
			double aVdW;
			double hsd;
			
			vector<int> Npoints;
			vector<double> densitySolid;
			vector<double> numParticlesSolid;
			vector<double> freeEnergySolid;
			
			for (int NpointsCounter=0; NpointsCounter<maxFileIndex; NpointsCounter++)
			{
				////////////////// Check if file exists ////////////////////
				
				ifstream dataFile(dataDir1+"/results_"+to_string(NpointsCounter)+
				                                   "_"+to_string(kTCounter)+
				                                   "_"+to_string(muCounter)+".dat");
				if (!dataFile) continue;
				else existsFileForThatkTmu = true;
				
				//////////////// Read computation results //////////////////
				
				int Npoints2;
				double densitySolid2;
				double numParticlesSolid2;
				double freeEnergySolid2;
				int success2;
				
				readDataFromFile(dataFile, "Npoints", Npoints2);
				readDataFromFile(dataFile, "kT", kT);
				readDataFromFile(dataFile, "mu", mu);
				readDataFromFile(dataFile, "aVdW", aVdW); // independant of Npoints
				readDataFromFile(dataFile, "hsd", hsd);   // independant of Npoints
				readDataFromFile(dataFile, "densitySolid", densitySolid2);
				readDataFromFile(dataFile, "numParticlesSolid", numParticlesSolid2);
				readDataFromFile(dataFile, "freeEnergySolid", freeEnergySolid2);
				readDataFromFile(dataFile, "successSolid", success2);
				
				dataFile.close();
				
				//////////////// Keep result if successful /////////////////
				
				if (success2==1 && densitySolid2>densityLimitSolid)
				{
					Npoints.push_back(Npoints2);
					densitySolid.push_back(densitySolid2);
					numParticlesSolid.push_back(numParticlesSolid2);
					freeEnergySolid.push_back(freeEnergySolid2);
				}
				
			} // end loop Npoints
			
			if (!existsFileForThatkTmu) continue;
			
			///////////////////////////// Report ///////////////////////////
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "kT = " << kT << " " << "mu = " << mu << endl;
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "Npoints = ";
			for (int i=0; i<Npoints.size(); i++) 
				log << Npoints[i] << " ";
			log << endl;
			log << "densitySolid = ";
			for (int i=0; i<densitySolid.size(); i++) 
				log << densitySolid[i] << " ";
			log << endl;
			log << "numParticlesSolid = ";
			for (int i=0; i<numParticlesSolid.size(); i++) 
				log << numParticlesSolid[i] << " ";
			log << endl;
			log << "freeEnergySolid = ";
			for (int i=0; i<freeEnergySolid.size(); i++) 
				log << freeEnergySolid[i] << " ";
			log << endl;
			
			////////////////////////// Minimisation ////////////////////////
			
			bool successMinimisation = true;
			double Npoints_min;
			double densitySolid_min;
			double numParticlesSolid_min;
			double freeEnergySolid_min;
			
			// We minimize over Npoints
			
			int statusMinNpoints = minFromDataParabola(Npoints, 
				freeEnergySolid, Npoints_min, freeEnergySolid_min);
			
			if (statusMinNpoints!=0)
			{
				successMinimisation = false;
				log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
				log << "ERROR: The set of lattice sizes is inappropriate, the minimum is either out of the set or too close to its boundaries" << endl;
			}
			
			// We evaluate other quantities at the minimum
			
			if (successMinimisation)
			{
				evalFromDataInterpolation(Npoints, densitySolid, 
					Npoints_min, densitySolid_min);
				evalFromDataInterpolation(Npoints, numParticlesSolid, 
					Npoints_min, numParticlesSolid_min);
			}
			
			// Report
			
			log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "(Npoints) Min Npoints = " << Npoints_min << endl;
			log << "(Npoints) Min Density Solid = " << densitySolid_min << endl;
			log << "(Npoints) Min Number of Particles Solid = " << numParticlesSolid_min << endl;
			log << "(Npoints) Min Free Energy Solid = " << freeEnergySolid_min << endl;
			
			////////////// Save minimisation results in file ///////////////
			
			ofstream dataFile2(dataDir2+"/results_"+to_string(kTCounter)+
			                                    "_"+to_string(muCounter)+".dat");
			
			dataFile2 << "kT = " << kT << endl;
			dataFile2 << "mu = " << mu << endl;
			dataFile2 << endl;
			dataFile2 << "aVdW = " << aVdW << endl;
			dataFile2 << "hsd = " << hsd << endl;
			dataFile2 << endl;
			dataFile2 << "densitySolid = " << densitySolid_min << endl;
			dataFile2 << "numParticlesSolid = " << numParticlesSolid_min << endl;
			dataFile2 << "freeEnergySolid = " << freeEnergySolid_min << endl;
			dataFile2 << endl;
			dataFile2 << "successSolid = " << successMinimisation << endl;
			
			dataFile2.close();
			
			
			/////////// Plot free energy profiles with Npoints /////////////
			
			// save plot data in file
			
			ofstream plotFile(plotDir+"/freeEnergySolid_"+to_string(kTCounter)+
			                                          "_"+to_string(muCounter)+".dat");
			
			plotFile << "# kT=" << kT << " mu=" << mu << endl;
			plotFile << "# Npoints freeEnergySolid" << endl;
			plotFile << endl;
			
			for (int i=0; i<Npoints.size(); i++)
				plotFile << Npoints[i] << " " << freeEnergySolid[i] << endl;
			
			plotFile.close();
			
			// create the gnuplot script 
			
			plotFile.open(plotDir+"/plot");
			
			plotFile << "#----------------------------------------" << endl;
			plotFile << "set term svg enhanced mouse #size 600,500" << endl;
			plotFile << "set output 'freeEnergySolid_" << kTCounter << "_" << muCounter << ".svg'" << endl;
			plotFile << endl;
			plotFile << "set title \"Solid Free Energy Profile with Npoints\" font \",20\"" << endl;
			plotFile << "set label \"kT = " << kT << " mu = " << mu << "\" at graph 0.1,0.9 font \",16\"" << endl; //0.29,1.05 is under title
			plotFile << "set xlabel \"Npoints\" font \",20\"" << endl;
			plotFile << "set ylabel \"Omega/(kT*V)\" font \",20\"" << endl;
			plotFile << endl;
			plotFile << "set key off" << endl;
			plotFile << endl;
			plotFile << "plot \"freeEnergySolid_" << to_string(kTCounter) << "_" 
			         << to_string(muCounter) << ".dat\" with linespoints" << endl;
			plotFile << endl;
			
			plotFile.close();
			
			sysresult = system(("cd "+plotDir+"; gnuplot plot").c_str());
			
		} // end loop mu
		
	} // end loop kT
	
	return 0;
}

