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
#include "../UniformPhases.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute coexistence properties of the solid and fluid phases.

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double maxFileIndex = 100;
	double rhoFluidMin = 0.001;
	double rhoFluidMax = 1.5;
	
	//TODO: Change the way the vapor-solid and liquid-solid curves are
	//      seperated -> not very clean
	double densityLimitVL = 0.3;
	
	// options
	
	Options options;
	options.addOption("MaxFileIndex", &maxFileIndex);
	options.addOption("RhoFluidMin", &rhoFluidMin);
	options.addOption("RhoFluidMax", &rhoFluidMax);
	options.addOption("DensityLimitVL", &densityLimitVL);
	options.read(argc, argv);
	
	// log
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	// files and directories
	
	string dataDir  = "Data-kt-mu-Solid";
	string plotDir1 = "Plots-omega-diff";
	string plotDir2 = "Plots-coexistence-curve";
	
	int sysresult = system(("mkdir -p " + plotDir1).c_str());
	    sysresult = system(("mkdir -p " + plotDir2).c_str());
	
	
	////////// Loop over all possible files in Solid directory  ////////////
	
	vector<double> kT_vec;
	vector<double> muCoex;
	vector<double> densityFluidCoex; 
	vector<double> densitySolidCoex;
	vector<double> freeEnergyCoex;
	vector<double> numParticlesSolidCoex;
	vector<double> success_vec;
	
	for (int kTCounter=0; kTCounter<maxFileIndex; kTCounter++)
	{
		bool existsFileForThatkT = false;
		
		double kT = 0;
		double aVdW, hsd;
		vector<double> mu;
		vector<double> densitySolid;
		vector<double> freeEnergyFluid;
		vector<double> freeEnergySolid;
		vector<double> numParticlesSolid;
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "kTCounter = " << kTCounter << endl;
		
		for (int muCounter=0; muCounter<maxFileIndex; muCounter++)
		{
			////////////////// Check if file exists ////////////////////
			
			string fileName = dataDir+"/results_"+to_string(kTCounter)+
			                                  "_"+to_string(muCounter)+".dat";
			ifstream dataFile(fileName);
			
			if (!dataFile)
			{
				log << "INFO: Unable to open file \"" << fileName << "\"" << endl;
				continue;
			}
			else 
			{
				log << "INFO: Successfully opened file \"" << fileName << "\"" << endl;
				existsFileForThatkT = true;
			}
			
			//////////////// Read computation results //////////////////
			
			bool successReading = true;
			
			double mu2;
			double densitySolid2;
			double numParticlesSolid2;
			double freeEnergySolid2;
			int success2;
			
			try
			{
				readDataFromFile(dataFile, "kT", kT);
				readDataFromFile(dataFile, "mu", mu2);
				readDataFromFile(dataFile, "aVdW", aVdW);
				readDataFromFile(dataFile, "hsd", hsd);
				readDataFromFile(dataFile, "densitySolid", densitySolid2);
				readDataFromFile(dataFile, "numParticlesSolid", numParticlesSolid2);
				readDataFromFile(dataFile, "freeEnergySolid", freeEnergySolid2);
				readDataFromFile(dataFile, "successSolid", success2);
			}
			catch(...)
			{
				successReading = false;
				log << "ERROR: Failed to read computation results" << endl;
			}
			
			dataFile.close();
			
			//////////////// Keep result if successful /////////////////
			
			if (successReading && success2==1)
			{
				mu.push_back(mu2);
				densitySolid.push_back(densitySolid2);
				numParticlesSolid.push_back(numParticlesSolid2);
				freeEnergySolid.push_back(freeEnergySolid2);
			}
			
		} // end loop mu
		
		if (!existsFileForThatkT) continue;
		
		//////////////// Compute free energy of the fluid //////////////
		
		for (int i=0; i<mu.size(); i++)
		{
			vector<double> roots;
			findRootsdOmegadRhoSpinodal(kT, mu[i], aVdW, hsd, 
				rhoFluidMin, rhoFluidMax, roots);
			
			if (roots.size()==1)
			{
				freeEnergyFluid.push_back( uniformOmega(kT, mu[i], aVdW, 
					hsd, roots[0]) );
			}
			else if (roots.size()==2)
			{
				double freeEnergyVapor  = uniformOmega(
					kT, mu[i], aVdW, hsd, roots[0]);
				double freeEnergyLiquid = uniformOmega(
					kT, mu[i], aVdW, hsd, roots[1]);
				
				if (freeEnergyVapor<freeEnergyLiquid)
					freeEnergyFluid.push_back( freeEnergyVapor );
				else
					freeEnergyFluid.push_back( freeEnergyLiquid );
			}
			else
			{
				log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
				log << "ERROR: Anormal number of roots for the fluid dOmega/drho" << endl;
			}
		}
		
		///////////////////////////// Report ///////////////////////////
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "kT = " << kT << endl;
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "mu = ";
		for (int i=0; i<mu.size(); i++) 
			log << mu[i] << " ";
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
		log << "freeEnergyFluid = ";
		for (int i=0; i<freeEnergyFluid.size(); i++) 
			log << freeEnergyFluid[i] << " ";
		log << endl;
		
		///////////////////// Free Energy Differences //////////////////////
		
		// Compute free energy difference and save in file for plotting
		// omega fluid - omega solid
		
		vector<double> freeEnergyDiff(mu.size()); 
		
		for (int i=0; i<mu.size(); i++)
			freeEnergyDiff[i] = freeEnergyFluid[i] - freeEnergySolid[i];
		
		ofstream plotFile(plotDir1+"/freeEnergyDiff_"+to_string(kT_vec.size())+".dat");
		
		plotFile << "# kT = " << kT << endl;
		plotFile << "# mu freeEnergyDiff" << endl;
		
		for (int i=0; i<mu.size(); i++)
			plotFile << mu[i] << " " << freeEnergyDiff[i] << endl;
		
		// Report
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "freeEnergyDiff = ";
		for (int i=0; i<freeEnergyDiff.size(); i++) 
			log << freeEnergyDiff[i] << " ";
		log << endl;
		
		///////////////////////// Find Coexistence /////////////////////////
		
		// find the coexistence value for mu by searching the zero of
		// omega fluid - omega solid
		
		bool success = true;
		
		double muCoex2 = 0;
		double densitySolidCoex2 = 0;
		double freeEnergyCoex2 = 0;
		double numParticlesSolidCoex2 = 0;
		
		// check if has the same sign before searching the zero
		if (freeEnergyDiff[0] * freeEnergyDiff[freeEnergyDiff.size()-1] < 0)
		{
			zeroFromDataInterpolation(mu, freeEnergyDiff, muCoex2);
			
			evalFromDataInterpolation(mu, densitySolid, muCoex2, densitySolidCoex2);
			evalFromDataInterpolation(mu, freeEnergySolid, muCoex2, freeEnergyCoex2);
			evalFromDataInterpolation(mu, numParticlesSolid, muCoex2, numParticlesSolidCoex2);
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Range invalid for kT = " << kT 
			    << ", does not contain the coexistence mu" << endl;
			success = false;
		}
		
		///////////////// Fluid density at coexistence /////////////////////
		
		double densityFluidCoex2;
		
		vector<double> roots;
		findRootsdOmegadRhoSpinodal(kT, muCoex2, aVdW, hsd, 
			rhoFluidMin, rhoFluidMax, roots);
		
		if (roots.size()==1)
		{
			densityFluidCoex2 = roots[0];
		}
		else if (roots.size()==2)
		{
			double freeEnergyVapor  = uniformOmega(
				kT, muCoex2, aVdW, hsd, roots[0]);
			double freeEnergyLiquid = uniformOmega(
				kT, muCoex2, aVdW, hsd, roots[1]);
			
			if (freeEnergyVapor<freeEnergyLiquid)
				densityFluidCoex2 = roots[0];
			else
				densityFluidCoex2 = roots[1];
		}
		else
		{
			log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
			log << "ERROR: Anormal number of roots for the fluid dOmega/drho" << endl;
		}
		
		// Report 
		
		log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "muCoex = " << muCoex2 << endl;
		log << "densityFluidCoex = " << densityFluidCoex2 << endl;
		log << "densitySolidCoex = " << densitySolidCoex2 << endl;
		log << "freeEnergyCoex = " << freeEnergyCoex2 << endl;
		log << "numParticlesSolidCoex = " << numParticlesSolidCoex2 << endl;
		log << "success = " << success << endl;
		
		// Save
		
		kT_vec.push_back(kT);
		muCoex.push_back(muCoex2);
		densityFluidCoex.push_back(densityFluidCoex2);
		densitySolidCoex.push_back(densitySolidCoex2);
		freeEnergyCoex.push_back(freeEnergyCoex2);
		numParticlesSolidCoex.push_back(numParticlesSolidCoex2);
		success_vec.push_back(success);
		
	} // end loop kT
	
	
	///////////////////// Save coexistence properties //////////////////////
	
	// save in file in a format that can be plotted with gnuplot
	// sort according to the temperatures
	
	ofstream plotFile(plotDir2+"/coexistenceFS.dat");
	
	plotFile << "# rhoF rhoS muCoex kT" << endl;
	
	for (int i=0; i<kT_vec.size(); i++)
		if (success_vec[i])
			plotFile << densityFluidCoex[i] << " " << densitySolidCoex[i] 
			         << " " << muCoex[i] << " " << kT_vec[i] << endl;
	
	plotFile.close();
	
	sysresult = system(("sort -k 4 "+plotDir2+"/coexistenceFS.dat "
	                    +"> "+plotDir2+"/coexistenceFS_temp.dat").c_str());
	sysresult = system(("mv "+plotDir2+"/coexistenceFS_temp.dat "
	                    +plotDir2+"/coexistenceFS.dat").c_str());
	                    
	// do also one with seperation of the vapor-solid and liquid-solid curves
	
	ofstream plotFileVS(plotDir2+"/coexistenceVS.dat");
	ofstream plotFileLS(plotDir2+"/coexistenceLS.dat");
	
	plotFileVS << "# rhoV rhoS muCoex kT" << endl;
	plotFileLS << "# rhoL rhoS muCoex kT" << endl;
	
	for (int i=0; i<kT_vec.size(); i++)
		if (success_vec[i] && densityFluidCoex[i]<densityLimitVL)
			plotFileVS << densityFluidCoex[i] << " " << densitySolidCoex[i] 
			           << " " << muCoex[i] << " " << kT_vec[i] << endl;
		else if (success_vec[i])
			plotFileLS << densityFluidCoex[i] << " " << densitySolidCoex[i] 
			           << " " << muCoex[i] << " " << kT_vec[i] << endl;
	
	plotFileVS.close();
	plotFileLS.close();
	
	sysresult = system(("sort -k 4 "+plotDir2+"/coexistenceVS.dat "
	                    +"> "+plotDir2+"/coexistenceVS_temp.dat").c_str());
	sysresult = system(("sort -k 4  "+plotDir2+"/coexistenceLS.dat "
	                    +"> "+plotDir2+"/coexistenceLS_temp.dat").c_str());
	sysresult = system(("mv "+plotDir2+"/coexistenceVS_temp.dat "
	                    +plotDir2+"/coexistenceVS.dat").c_str());
	sysresult = system(("mv "+plotDir2+"/coexistenceLS_temp.dat "
	                    +plotDir2+"/coexistenceLS.dat").c_str());
	
	/////////////////////////// Gnuplot scripts ////////////////////////////
	
	// free energy differences
	
	for (int j=0; j<kT_vec.size(); j++)
	{
		plotFile.open(plotDir1+"/plot");
		
		plotFile << "#----------------------------------------" << endl;
		plotFile << "set term svg enhanced mouse #size 600,500" << endl;
		plotFile << "set output 'freeEnergyDiff_" << j << ".svg'" << endl;
		plotFile << endl;
		plotFile << "set title \"Free Energy Difference with mu\" font \",20\"" << endl;
		plotFile << "set label \"kT = " << kT_vec[j] << "\" at graph 0.1,0.9 font \",16\"" << endl; //0.29,1.05 is under title
		plotFile << "set xlabel \"mu\" font \",20\"" << endl;
		plotFile << "set ylabel \"(OmegaF-OmegaS)/(kT*V)\" font \",20\"" << endl;
		plotFile << endl;
		plotFile << "set key off" << endl;
		plotFile << endl;
		plotFile << "plot \"freeEnergyDiff_" << to_string(j) << ".dat\" with linespoints" << endl;
		plotFile << endl;
		
		plotFile.close();
		
		sysresult = system(("cd "+plotDir1+"; gnuplot plot").c_str());
	}
	
	// coexistence curve with densities (fluid-solid)
	
	plotFile.open(plotDir2+"/plot_rhoFS");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'coexistenceFS.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Fluid-Solid coexistence curve\" font \",20\"" << endl;
	plotFile << "set xlabel \"rho\" font \",20\"" << endl;
	plotFile << "set ylabel \"kT\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key top left" << endl;
	plotFile << endl;
	plotFile << "plot \"coexistenceFS.dat\" using 1:4 with linespoints title 'fluid', \\" << endl;
	plotFile << "     \"coexistenceFS.dat\" using 2:4 with linespoints title 'solid'" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+plotDir2+"; gnuplot plot_rhoFS").c_str());
	
	// coexistence curve with chemical potentials
	
	plotFile.open(plotDir2+"/plot_muFS");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'coexistenceFS_mu.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Fluid-Solid coexistence curve\" font \",20\"" << endl;
	plotFile << "set xlabel \"mu\" font \",20\"" << endl;
	plotFile << "set ylabel \"kT\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "plot \"coexistenceFS.dat\" using 3:4 with linespoints" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+plotDir2+"; gnuplot plot_muFS").c_str());
	
	// vapor-solid coexistence curve (densities)
	
	plotFile.open(plotDir2+"/plot_rhoVS");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'coexistenceVS.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Vapor-Solid coexistence curve\" font \",20\"" << endl;
	plotFile << "set xlabel \"rho\" font \",20\"" << endl;
	plotFile << "set ylabel \"kT\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key top left" << endl;
	plotFile << endl;
	plotFile << "plot \"coexistenceVS.dat\" using 1:4 with linespoints title 'vapor', \\" << endl;
	plotFile << "     \"coexistenceVS.dat\" using 2:4 with linespoints title 'solid'" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+plotDir2+"; gnuplot plot_rhoVS").c_str());
	
	// liquid-solid coexistence curve (densities)
	
	plotFile.open(plotDir2+"/plot_rhoLS");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'coexistenceLS.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Liquid-Solid coexistence curve\" font \",20\"" << endl;
	plotFile << "set xlabel \"rho\" font \",20\"" << endl;
	plotFile << "set ylabel \"kT\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key top left" << endl;
	plotFile << endl;
	plotFile << "plot \"coexistenceLS.dat\" using 1:4 with linespoints title 'liquid', \\" << endl;
	plotFile << "     \"coexistenceLS.dat\" using 2:4 with linespoints title 'solid'" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+plotDir2+"; gnuplot plot_rhoLS").c_str());
	
	return 0;
}

