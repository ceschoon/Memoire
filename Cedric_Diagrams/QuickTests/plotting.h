#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

using namespace std;

void plotInFile(vector<double> x, vector<double> y, string filename,
                string title="", string xlabel="", string ylabel="")
{
	// save data in file
	
	ofstream dataFile(filename+".dat");
	
	dataFile << "# x y" << endl;
	
	for (int i=0; i<x.size(); i++)
		dataFile << x[i] << " " << y[i] << endl;
	
	dataFile.close();
	
	// create gnuplot script
	
	ofstream plotFile(filename+".gp");
	
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output '" << filename << ".svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"" << title << "\" font \",20\"" << endl;
	plotFile << "set xlabel \"" << xlabel << "\" font \",20\"" << endl;
	plotFile << "set ylabel \"" << ylabel << "\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	
	plotFile << "plot \"" << filename << ".dat\" with lines" << endl;
	
	plotFile.close();
	
	// execute script and clean afterwards
	
	int sysres = system(("gnuplot "+filename+".gp").c_str());
	    sysres = system(("rm "+filename+".dat").c_str());
	    sysres = system(("rm "+filename+".gp").c_str());
}