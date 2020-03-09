#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

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

#include "../SolidPhase.h"
#include "../UniformPhases.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to compute the free energies at fixed (kT,mu)

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	double kT_rangeMin = 1;
	double kT_rangeMax = 1;
	double kT_rangeStep = 1;
	
	double mu_rangeMin = 1;
	double mu_rangeMax = 1;
	double mu_rangeStep = 1;
	
	int Ngrid_rangeMin = 120;
	int Ngrid_rangeMax = 150;
	int Ngrid_rangeStep = 5;
	
	Options options;
	
	options.addOption("kT_rangeMin", &kT_rangeMin);
	options.addOption("kT_rangeMax", &kT_rangeMax);
	options.addOption("kT_rangeStep", &kT_rangeStep);
	
	options.addOption("mu_rangeMin", &mu_rangeMin);
	options.addOption("mu_rangeMax", &mu_rangeMax);
	options.addOption("mu_rangeStep", &mu_rangeStep);
	
	options.addOption("Ngrid_rangeMin", &Ngrid_rangeMin);
	options.addOption("Ngrid_rangeMax", &Ngrid_rangeMax);
	options.addOption("Ngrid_rangeStep", &Ngrid_rangeStep);
	
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	///////////////////////// Loop over variables //////////////////////////
	
	for (double kT=kT_rangeMin; kT<=kT_rangeMax; kT+=kT_rangeStep)
		for (double mu=mu_rangeMin; mu<=mu_rangeMax; mu+=mu_rangeStep)
			for (int Ngrid=Ngrid_rangeMin; Ngrid<=Ngrid_rangeMax; Ngrid+=Ngrid_rangeStep)
				surfacePlot( kT, mu, Ngrid, argc, argv, log);
	
	return 0;
}
