#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "../utilities.h"
#include "options.h"
#include "Log.h"
#include "myColor.h"


using namespace std; 

double min(double a, double b) {return (a<b)?a:b;}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program that puts the coexistence curves together in a phase diagram

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDir = "data";
	
	Options options;
	options.addOption("DataDirectory", &dataDir);
	options.read(argc, argv);
	
	// no need to add to options
	double kTMin = 0;
	double kTMax = 2;
	double kTStep = 1e-4; // precision of data saves
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	///////////////////////////// Read files ///////////////////////////////
	
	ofstream outDataFile("coexCurveVS.dat");
	outDataFile << "#kT				rhoV			rhoS			"
	            << "muCoex			freeEnergy		Ngrid			"
	            << "Cvac			alpha			"
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	for (double kT=kTMin; kT<=kTMax; kT+=kTStep)
	{
		//log << endl;
		//log << "Checking for a valid computation at kT=" << kT << endl;
		
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		ifstream dataInFile(dataDir+"/coexistenceVS_"+"kT="+sskT.str()+".dat");
		
		if (dataInFile)
		{
			// check if it was a successful computation
			int success;
			readDataFromFile(dataInFile, "success", success);
			
			if (success==1)  // bool "true" 
			{
				log << endl;
				log << "Found valid computation for kT=" << kT << endl;
				
				double rhoV, rhoS, muCoex, freeEnergyCoex, Ngrid, Cvac, alpha;
				readDataFromFile(dataInFile, "muCoex", muCoex);
				readDataFromFile(dataInFile, "freeEnergyCoex", freeEnergyCoex);
				readDataFromFile(dataInFile, "densityVapourCoex", rhoV);
				readDataFromFile(dataInFile, "densitySolidCoex", rhoS);
				readDataFromFile(dataInFile, "NgridCoex", Ngrid);
				readDataFromFile(dataInFile, "CvacCoex", Cvac);
				readDataFromFile(dataInFile, "alphaCoex", alpha);
				
				outDataFile << kT << " 	" << rhoV << " 	" << rhoS << " 	"
				            << muCoex << " 	" << freeEnergyCoex << " 	"
				            << Ngrid << " 	" << Cvac << " 	" << alpha << " 	"
				            << endl;
			}
		}
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Coexistence Curve #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output 'coexCurveVS.svg'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Coexistence Vapour-Solid\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"rho\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"kT\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot 'coexCurveVS.dat' using 2:1 title 'vapour', \\" << endl;
	outPlotFile << "     'coexCurveVS.dat' using 3:1 title 'solid' " << endl;
	
	int sysresult = system("gnuplot plot");
	
	////////////////////////////////////////////////////////////////////////
	
	
	return 0;
}

