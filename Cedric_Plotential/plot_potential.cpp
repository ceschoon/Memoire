#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <stdlib.h>

#include "Potential1.h"

using namespace std;


// Program to plot the potential

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// files
	int sysresult = system("mkdir -p data");
	ofstream dataFile("data/potentialDataPlot.dat");
	dataFile << "r				LJ (rc=3)		WHDF (rc=2)		"
	         << "WHDF (rc=1.2)	tWF (rc=3)		"
	         << scientific << endl;
	
	// LJ potential
	double eps1   = 1;
	double sigma1 = 1;
	double rcut1  = 3;
	double kT = 1;
	
	// WHDF potential r=2
	double eps2   = 1;
	double sigma2 = 1;
	double rcut2  = 2;
	
	// WHDF potential r=1.2
	double eps3   = 1;
	double sigma3 = 1;
	double rcut3  = 1.2;
	
	// tWF potential
	double eps4   = 1;
	double sigma4 = 1;
	double rcut4  = 3;
	
	// range
	double rMin = 0.9;
	double rMax = 2.1;
	double rStep = 0.001;
	
	//////////////////////////// Data to plot //////////////////////////////
	
	// define potentials
	LJ   potential1(sigma1, eps1, rcut1);
	WHDF potential2(sigma2, eps2, rcut2);
	WHDF potential3(sigma3, eps3, rcut3);
	tWF  potential4(sigma4, eps4, rcut4);
	
	// print potential minima
	cout << endl;
	cout << "--------------------------" << endl;
	cout << "LJ:             rmin = " << potential1.getRmin() << endl;
	cout << "WHDF (rc=2):    rmin = " << potential2.getRmin() << endl;
	cout << "WHDF (rc=1.2):  rmin = " << potential3.getRmin() << endl;
	cout << "tWF:            rmin = " << potential4.getRmin() << endl;
	
	// print potential hsd
	cout << endl;
	cout << "--------------------------" << endl;
	cout << "hsd for kT = " << kT << endl;
	cout << "LJ:             hsd = " << potential1.getHSD(kT) << endl;
	cout << "WHDF (rc=2):    hsd = " << potential2.getHSD(kT) << endl;
	cout << "WHDF (rc=1.2):  hsd = " << potential3.getHSD(kT) << endl;
	cout << "tWF:            hsd = " << potential4.getHSD(kT) << endl;
	
	// print potential aVdW
	cout << endl;
	cout << "--------------------------" << endl;
	cout << "aVdW for kT = " << kT << endl;
	cout << "LJ:             aVdW = " << potential1.getVDW_Parameter(kT) << endl;
	cout << "WHDF (rc=2):    aVdW = " << potential2.getVDW_Parameter(kT) << endl;
	cout << "WHDF (rc=1.2):  aVdW = " << potential3.getVDW_Parameter(kT) << endl;
	cout << "tWF:            aVdW = " << potential4.getVDW_Parameter(kT) << endl;
	
	// plot v(r)
	for (double r=rMin; r<rMax; r+=rStep)
	{
		dataFile << r << " 	"
		         << potential1.V(r) << " 	"
		         << potential2.V(r) << " 	"
		         << potential3.V(r) << " 	"
		         << potential4.V(r) << " 	"
		         << endl;
	}
	
	/////////////////////////// gnuplot script /////////////////////////////
	
	ofstream plotFile;
	plotFile.open("data/plot_potential");
	
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'potential.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Potential shapes\" font \",20\"" << endl;
	plotFile << "set xlabel \"r/sigma\" font \",20\"" << endl;
	plotFile << "set ylabel \"V(r)/epsilon\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set xrange [" << rMin << ":" << rMax << "]" << endl;
	plotFile << "set yrange [-1.2:1.2]" << endl;
	plotFile << endl;
	plotFile << "set key top right" << endl;
	plotFile << endl;
	
	plotFile << "plot \"potentialDataPlot.dat\" using 1:2 with lines "
	         << "title \"LJ rc=3\" ,\\" << endl;
	plotFile << "     \"potentialDataPlot.dat\" using 1:3 with lines "
	         << "title \"WHDF rc=2\" ,\\" << endl;
	plotFile << "     \"potentialDataPlot.dat\" using 1:4 with lines "
	         << "title \"WHDF rc=1.2\" ,\\" << endl;
	//plotFile << "     \"potentialDataPlot.dat\" using 1:5 with lines "
	//         << "title \"tWF rc=3\"" << endl;
	
	plotFile.close();
	sysresult = system("cd data; gnuplot plot_potential");
	
	return 0;
}

