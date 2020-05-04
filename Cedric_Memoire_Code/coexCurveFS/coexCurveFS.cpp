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

// Program that collects the coexistence data in a file and draws
// the coexistence curve

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataDirIn = "dataIn";
	string dataDirOut = "dataOut";
	string fluidType = "fluid";
	
	Options options;
	options.addOption("DataDirectoryIn", &dataDirIn);
	options.addOption("DataDirectoryOut", &dataDirOut);
	options.addOption("FluidType", &fluidType);
	options.read(argc, argv);
	
	// no need to add to options
	double kTMin = 0;
	double kTMax = 5;
	double kTStep = 1e-4; // precision of data saves
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	int sysresult = system(("mkdir -p "+dataDirIn ).c_str());
	    sysresult = system(("mkdir -p "+dataDirOut).c_str());
	
	///////////////////////////// Fluid type ///////////////////////////////
	
	string fluid_str = "";
	string Fluid_str = "";
	string rhoF_str = "";
	string fileName = "";
	
	if (fluidType == "fluid")  
	{
		fluid_str = "fluid"; 
		Fluid_str = "Fluid";
		rhoF_str = "rhoF";
		fileName = "coexCurveFS";
	}
	else if (fluidType == "vapour")  
	{
		fluid_str = "vapour"; 
		Fluid_str = "Vapour";
		rhoF_str = "rhoV";
		fileName = "coexCurveVS";
	}
	else if (fluidType == "liquid")  
	{
		fluid_str = "liquid"; 
		Fluid_str = "Liquid";
		rhoF_str = "rhoL";
		fileName = "coexCurveLS";
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: Fluid Type Unrecognised" << endl;
	}
	
	///////////////////////////// Read files ///////////////////////////////
	
	ofstream outDataFile(dataDirOut+"/"+fileName+".dat");
	outDataFile << "#kT				"+rhoF_str+"			rhoS			"
	            << "muCoex			freeEnergy		Ngrid			"
	            << "Cvac			alpha			" 
	            << endl;
	outDataFile << scientific << setprecision(8);
	
	for (double kT=kTMin; kT<=kTMax; kT+=kTStep)
	{
		//log << endl;
		//log << "Checking for a valid computation at kT=" << kT << endl;
		
		stringstream sskT; sskT << scientific << setprecision(4) << kT;
		ifstream dataInFile(dataDirIn+"/coexistenceFS_"+"kT="+sskT.str()+".dat");
		
		if (dataInFile)
		{
			// check if it was a successful computation
			int success;
			readDataFromFile(dataInFile, "success", success);
			
			if (success==1)  // bool "true" 
			{
				log << endl;
				log << "Found valid computation for kT=" << kT << endl;
				
				double rhoF, rhoS, muCoex, freeEnergyCoex, Ngrid, Cvac, alpha;
				readDataFromFile(dataInFile, "muCoex", muCoex);
				readDataFromFile(dataInFile, "freeEnergyCoex", freeEnergyCoex);
				readDataFromFile(dataInFile, "densityFluidCoex", rhoF);
				readDataFromFile(dataInFile, "densitySolidCoex", rhoS);
				readDataFromFile(dataInFile, "NgridCoex", Ngrid);
				readDataFromFile(dataInFile, "CvacCoex", Cvac);
				readDataFromFile(dataInFile, "alphaCoex", alpha);
				
				outDataFile << kT << " 	" << rhoF << " 	" << rhoS << " 	"
				            << muCoex << " 	" << freeEnergyCoex << " 	"
				            << Ngrid << " 	" << Cvac << " 	" << alpha << " 	"
				            << endl;
			}
		}
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream plotFile(dataDirOut+"/"+fileName+".gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output '"+fileName+".tex'" << endl;
	plotFile << endl;
	plotFile << "set title 'Coexistence "+Fluid_str+"-Solid'" << endl;
	plotFile << "set xlabel '$\\rho$'" << endl;
	plotFile << "set ylabel '$kT$'" << endl;
	plotFile << "set style data lines" << endl;
	plotFile << "set grid" << endl;
	plotFile << endl;
	plotFile << "set key bottom center" << endl;
	plotFile << endl;
	plotFile << "plot '"+fileName+".dat' using 2:1 title '"+fluid_str+"', \\" << endl;
	plotFile << "     '"+fileName+".dat' using 3:1 title 'solid'    " << endl;
	plotFile << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	ofstream shellFile(dataDirOut+"/plot_everything.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   "+fileName+".gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  "+fileName+".tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   "+fileName+".aux   "+fileName+".log   "+fileName+".tex   "+fileName+"-inc.eps   "+fileName+"-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	sysresult = system(("chmod +x "+dataDirOut+"/plot_everything.sh; cd "+dataDirOut+"; ./plot_everything.sh").c_str());
	
	
	return 0;
}

