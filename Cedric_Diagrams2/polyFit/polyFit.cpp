#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_multifit.h>

#include "../utilities.h"
#include "options.h"
#include "Log.h"
#include "myColor.h"


using namespace std; 


int myPolyFit(vector<double> data_x, vector<double> data_y,
              vector<double> &coefficients, int orderPolyFit);


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program to fit data with polynomials
// We use the gsl non linear fit functions

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	string dataFileName = "You must specify the data file in the options";
	string coefFileName = "fitCoefs.dat";
	string polyFileName = "fitCurve.dat";
	string plotFileName = "fitCurve.svg";
	int orderPolyFit = 4;
	int colNumDataX = 0;
	int colNumDataY = 1;
	
	Options options;
	options.addOption("dataFileName", &dataFileName);
	options.addOption("coefFileName", &coefFileName);
	options.addOption("polyFileName", &polyFileName);
	options.addOption("plotFileName", &plotFileName);
	options.addOption("orderPolyFit", &orderPolyFit);
	options.addOption("colNumDataX", &colNumDataX);
	options.addOption("colNumDataY", &colNumDataY);
	options.read(argc, argv);
	
	Log log("log.dat");
	log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;
	options.write(log);
	log << myColor::GREEN << "=================================" << myColor::RESET << endl;
	
	
	//////////////// Loop over all files in data directory /////////////////
	/*
	DIR *dir;
	struct dirent *ent;

	if ((dir = opendir (dataDir.c_str())) != NULL) 
	{
		while ((ent = readdir (dir)) != NULL) 
		{
			string ent_name(ent->d_name);
			int ent_size = ent_name.size();
			if (ent_size<=4) continue;
			
			int posLastDot;
			for (int i=0; i<ent_size; i++)
				if (ent_name[i]=='.') posLastDot = i;
			
			string ent_namebody = ent_name.substr(0,ent_size-posLastDot);
			string ent_extension = ent_name.substr(ent_size-posLastDot,ent_size);
			
			log << "Checking file:    " << ent_name      << endl;
			log << "file name (body): " << ent_namebody  << endl;
			log << "file extension:   " << ent_extension << endl;
			
			if (ent_extension == ".dat")
			{
				// Files
				// Read data from file
				// Fit
			}
		}
		
		closedir (dir);
	} 
	else 
	{
		// could not open directory 
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: could not open data directory (polyFit)" << endl;
	}
	*/
	
	
	///////////////////////// Read data from file //////////////////////////
	
	vector<double> data_x;
	vector<double> data_y;
	
	ifstream dataFile(dataFileName);
	if (dataFile)
	{
		readColumnVectorFromFile(dataFile, colNumDataX, data_x);
		readColumnVectorFromFile(dataFile, colNumDataY, data_y);
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: cannot read file: " << dataFileName << endl;
	}
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "data_x = ";
	for (int i=0; i<data_x.size(); i++) log << data_x[i] << " ";
	log << endl;
	log << "data_y = ";
	for (int i=0; i<data_y.size(); i++) log << data_y[i] << " ";
	log << endl;
	
	//////////////////////////////// poly Fit //////////////////////////////
	
	log << endl;
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << endl;
	
	ofstream coefFile(coefFileName);
	coefFile << scientific << setprecision(8);
	
	vector<double> coefficients;
	int statusFit = myPolyFit(data_x, data_y, coefficients, orderPolyFit);
	
	if (statusFit == 0)
	{
		// Write polynomial fit coefficients in file
		for (int i=0; i<coefficients.size(); i++)
			coefFile << "a" << i << " = " << coefficients[i] << endl;
	}
	else
	{
		log <<  myColor::RED << "=================================" << myColor::RESET << endl << "#" << endl;
		log << "ERROR: polynomial fit not successful for file: " 
		    << coefFileName << endl;
	}
	
	log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
	log << "coefficients = ";
	for (int i=0; i<coefficients.size(); i++) log << coefficients[i] << " ";
	log << endl;
	
	/////////////////////// save smooth polynomial curve ///////////////////
	
	ofstream polyFile(polyFileName);
	coefFile << scientific << setprecision(8);
	
	polyFile << "# x             " << " 	" << "y             " << endl;
	
	int NfitCurve = 1000;
	for (int i=0; i<NfitCurve; i++)
	{
		double x = data_x[0] * (1-i/double(NfitCurve-1)) + data_x[data_x.size()-1] * i/double(NfitCurve-1);
		
		double y = 0;
		for (int j=0; j<coefficients.size(); j++)
			y += coefficients[j] * pow(x,j);
		
		polyFile << x << " 	" << y << endl;
	}
	
	///////////////////////////// gnuplot script ///////////////////////////
	
	ofstream outPlotFile("plot");
	
	outPlotFile << "########## Plot Data Versus Fit #########" << endl;
	outPlotFile << endl;
	outPlotFile << "set term svg enhanced mouse #size 800,600" << endl;
	outPlotFile << "set output '" << plotFileName << "'" << endl;
	outPlotFile << endl;
	outPlotFile << "set title \"Data Versus Fit\" font \",20\"" << endl;
	outPlotFile << "set xlabel \"y\" font \",20\"" << endl;
	outPlotFile << "set ylabel \"x\" font \",20\"" << endl;
	outPlotFile << "set style data lines" << endl;
	outPlotFile << "set grid" << endl;
	outPlotFile << "set key top left" << endl;
	outPlotFile << endl;
	outPlotFile << "plot '" << dataFileName << "' using " << colNumDataX+1 
	            << ":" << colNumDataY+1 << " with points title 'data', \\" << endl;
	outPlotFile << "     '" << polyFileName << "' using 1:2 title 'fit'" << endl;
	
	int sysresult = system("gnuplot plot");
	
	////////////////////////////////////////////////////////////////////////
	
	return 0;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// polynomial fit (wrapper for the gsl code)


int myPolyFit(vector<double> data_x, vector<double> data_y,
              vector<double> &coefficients, int orderPolyFit)
{
	cout << "=================================" << endl;
	cout << "== Entering myPolyFit function ==" << endl;
	cout << "=================================" << endl;
	
	// Change units in order to have quantities of order 1
	
	double data_x_scale = data_x[0];
	double data_y_scale = data_y[0];
	
	cout << endl;
	cout << "data_x_scale = " << data_x_scale << endl;
	cout << "data_y_scale = " << data_y_scale << endl;
	
	for (int i=0; i<data_x.size(); i++) data_x[i] /= data_x_scale;
	for (int i=0; i<data_y.size(); i++) data_y[i] /= data_y_scale;
	
	// Convert input to GSL formats
	
	int n = data_x.size();
	int p = orderPolyFit+1;
	
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;

	X = gsl_matrix_alloc (n, p);
	y = gsl_vector_alloc (n);
	w = gsl_vector_alloc (n);
	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc (p, p);
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			gsl_matrix_set (X, i, j, pow(data_x[i],j));
			gsl_vector_set (y, i, data_y[i]);
			gsl_vector_set (w, i, 1.0);
		}
	}
	
	// Solve
	
	double chisq;
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);
	
	// Convert output from GSL format
	// don't forget to return to the original units
	
	coefficients.clear();
	for (int i=0; i<=orderPolyFit; i++)
		coefficients.push_back(gsl_vector_get(c,(i)));
	
	for (int i=0; i<=orderPolyFit; i++)
		coefficients[i] *= data_y_scale / pow(data_x_scale,i);
	
	// Deallocate memory
	
	gsl_matrix_free (X);
	gsl_vector_free (y);
	gsl_vector_free (w);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
	
	return 0;
}





