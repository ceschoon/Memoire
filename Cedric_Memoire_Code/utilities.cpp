#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>

#include "utilities.h"

using namespace std; 

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remove tabs and spaces from string

// Function used by readDataFromFile

void removeTabsAndSpaces(string &str)
{
	int N = str.size();
	string str_new = "";
	
	for (int i=0; i<N; i++)
		if (str[i]!=' ' && str[i]!='	')
			str_new = str_new + str[i];
	
	str = str_new;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Read data from file

// Expected format: parameter = value, all spaces and tabs are ignored

int readDataFromFile(ifstream &file, string parameter, double &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	bool foundParameter = false;
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter && line[b]=='=') 
		{
			value = stod(line.substr(b+1, a-b-1));
			foundParameter = true;
		}
	}
	
	if (!foundParameter) return 1;
	return 0;
}

int readDataFromFile(ifstream &file, string parameter, int &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	bool foundParameter = false;
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter && line[b]=='=')
		{ 
			value = stoi(line.substr(b+1, a-b-1));
			foundParameter = true;
		}
	}
	
	if (!foundParameter) return 1;
	return 0;
}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Read column vector

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<double> &vec)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	vec.empty();
	
	string line;
	
	while (getline(file, line))
	{
		char separator = ' ';
		char comment = '#';
		int indexStart = 0;
		int indexStop = line.size();
		int separatorCount = 0;
		
		if (line.size()>0) if(line[0]==comment) continue; // skip comments
		
		for (int i=0; i<line.size(); i++) 
		{
			if (line[i]==separator) 
			{
				separatorCount ++;
				if (separatorCount==columnNumber)   indexStart = i+1;
				if (separatorCount==columnNumber+1) indexStop = i-1;
			}
		}
		
		string value_str = line.substr(indexStart,indexStop-indexStart+1);
		vec.push_back(stod(value_str));
	}
}

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<int> &vec)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	vec.empty();
	
	string line;
	
	while (getline(file, line))
	{
		char separator = ' ';
		char comment = '#';
		int indexStart = 0;
		int indexStop;
		int separatorCount = 0;
		
		if (line.size()>0) if(line[0]==comment) continue; // skip comments
		
		for (int i=0; i<line.size(); i++) 
		{
			if (line[i]==separator) 
			{
				separatorCount ++;
				if (separatorCount==columnNumber)   indexStart = i+1;
				if (separatorCount==columnNumber+1) indexStop = i-1;
			}
		}
		
		string value_str = line.substr(indexStart,indexStop);
		vec.push_back(stoi(value_str));
	}
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Evaluate function sampled by a few data points
// Cubic Spline interpolation 

int evalFromDataInterpolation(vector<double> x, vector<double> y, double xIn,
                              double &yOut)
{
	// convert to double array
	
	int N = x.size();
	double x2[N];
	double y2[N];
	
	for (int i=0; i<N; i++)
	{
		x2[i] = x[i];
		y2[i] = y[i];
	}
	
	// cubic spline interpolation
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
	
	gsl_spline_init (spline, x2, y2, N);
	
	// evaluate function
	
	yOut = gsl_spline_eval (spline, xIn, acc);
	
	// free pointers
	
	gsl_interp_accel_free (acc);
	gsl_spline_free (spline);
	
	return 0;
}


int evalFromDataInterpolation(vector<int> x, vector<double> y, double xIn,
                              double &yOut)
{
	vector<double> xd(x.size(),0);
	for (int i=0; i<x.size(); i++) xd[i] = x[i];
	
	return evalFromDataInterpolation(xd,y,xIn,yOut);
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// functions and structures needed by the gsl algorithms

struct spline_param
{
	gsl_spline *spline;
};

double function_spline (double x, void *params)
{
	struct spline_param *p = (struct spline_param *) params;
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	
	double y = gsl_spline_eval (p->spline, x, acc);
	
	gsl_interp_accel_free (acc);
	
	return y;
}

double function_spline_deriv (double x, void *params)
{
	struct spline_param *p = (struct spline_param *) params;
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	
	double y = gsl_spline_eval_deriv (p->spline, x, acc);
	
	gsl_interp_accel_free (acc);
	
	return y;
}








////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Find zero of function sampled by a few data points

// First, fit the data set with a cubic spline
// Second, find the zero of the fit by bisection

// The function assumes that a zero exists between the min and max values
// of the x vector.

int zeroFromDataInterpolation(vector<double> x, vector<double> y, double &root)
{
	///////////// Find x-axis boundaries for the root finding //////////////
	
	int N = x.size();
	
	double xMin = x[0];
	double xMax = x[0];
	
	for (int i=0; i<N; i++)
	{
		if (xMin>x[i]) xMin = x[i];
		if (xMax<x[i]) xMax = x[i];
	}
	
	/////////////////////// Cubic Spline interpolation /////////////////////
	
	double x2[N];
	double y2[N];
	
	for (int i=0; i<N; i++)
	{
		x2[i] = x[i];
		y2[i] = y[i];
	}
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
	
	gsl_spline_init (spline, x2, y2, N);
	
	//////////////////////// Find zero of the spline ///////////////////////
	
	// use gsl bracketing algorithm
	
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0;
	double x_lo = xMin, x_hi = xMax;
	gsl_function F;
	struct spline_param params = {spline};

	F.function = &function_spline;
	F.params = &params;

	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	printf ("using %s method\n", 
	        gsl_root_fsolver_name (s));

	printf ("%5s [%9s, %9s] %9s %9s\n",
	        "iter", "lower", "upper", "root", "err(est)");

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-8);

		if (status == GSL_SUCCESS)
			printf ("Converged:\n");

		printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
		        iter, x_lo, x_hi, r, x_hi - x_lo);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	
	// save result
	root = r;
	
	// free pointers
	gsl_interp_accel_free (acc);
	gsl_spline_free (spline);
	
	return 0;
}






////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// TODO: Can fail because the cubic can have multiple extrema !!??

// Find minimum of function sampled by a few data points

// First, fit the data set with a cubic spline
// Second, find the minimum of the fit with a bisection root finding on
// the derivative

// returns 0 if all went well

int minFromDataCSpline(vector<double> x, vector<double> y, double &xMin,
                       double &yMin)
{
	//////////////////////////// Sanity Check //////////////////////////////
	
	// Check that the min point is not on the edge of the dataset
	
	int N = x.size();
	
	int iMin = 0;
	for (int i=0; i<N; i++)  if (y[i]<y[iMin])  iMin = i;
	
	if (iMin==0 || iMin==N-1)
	{
		cout << "ERROR: In minFromData: Minimum on the edge of the dataset" << endl;
		return 1;
	}
	
	/////////////////////// Cubic Spline interpolation /////////////////////
	
	double x2[N];
	double y2[N];
	
	for (int i=0; i<N; i++)
	{
		x2[i] = x[i];
		y2[i] = y[i];
	}
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
	
	gsl_spline_init (spline, x2, y2, N);
	
	////////////////// Find zero of the spline derivative //////////////////
	
	// use gsl bracketing algorithm
	// search the minimum between the two points surrounding the min data point
	
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0;
	double x_lo = x[iMin-1], x_hi = x[iMin+1];
	gsl_function F;
	struct spline_param params = {spline};

	F.function = &function_spline_deriv;
	F.params = &params;

	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	printf ("using %s method\n", 
	        gsl_root_fsolver_name (s));

	printf ("%5s [%9s, %9s] %9s %9s\n",
	        "iter", "lower", "upper", "root", "err(est)");

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-8);

		if (status == GSL_SUCCESS)
			printf ("Converged:\n");

		printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
		        iter, x_lo, x_hi, r, x_hi - x_lo);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	
	// save result
	xMin = r;
	yMin = gsl_spline_eval (spline, xMin, acc);
	
	// free pointers
	gsl_interp_accel_free (acc);
	gsl_spline_free (spline);
	
	return 0;
}

int minFromDataCSpline(vector<int> x, vector<double> y, double &xMin, 
                        double &yMin)
{
	vector<double> xd(x.size(),0);
	for (int i=0; i<x.size(); i++) xd[i] = x[i];
	
	return minFromDataCSpline(xd,y,xMin,yMin);
}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Find minimum of function sampled by a few data points

// First, we find the min data point
// Second, we interpolate the min point and its neighbours using a parabola
// Third, we select the min of the parabola

// returns 0 if all went well

int minFromDataParabola(vector<double> x, vector<double> y, double &xMin, 
                        double &yMin)
{
	/////////////////////// Find the min data point ////////////////////////
	
	int N = x.size();
	
	int iMin = 0;
	for (int i=0; i<N; i++)  if (y[i]<y[iMin])  iMin = i;
	
	if (iMin==0 || iMin==N-1)
	{
		cout << "ERROR: In minFromData: Minimum on the edge of the dataset" << endl;
		return 1;
	}
	
	/////////////////////// Parabolic Interpolation ////////////////////////
	
	// We want to find y = ax^2 + bx + c such as y(xi)=xi
	// This is solving the matrix problem
	// y1     x1^2 x1 1    a
	// y2  =  x2^2 x2 1    b
	// y3     x3^2 x3 1    c
	
	double x1 = x[iMin-1];
	double x2 = x[iMin];
	double x3 = x[iMin+1];
	
	double y1 = y[iMin-1];
	double y2 = y[iMin];
	double y3 = y[iMin+1];
	
	double det = (x1-x3)*(x3-x2)*(x2-x1);
	
	double a = - 1/det * (        (x3-x2)*y1 +         (x1-x3)*y2 +         (x2-x1)*y3);
	double b =   1/det * ((x3+x2)*(x3-x2)*y1 + (x1+x3)*(x1-x3)*y2 + (x2+x1)*(x2-x1)*y3);
	double c = - 1/det * (  x3*x2*(x3-x2)*y1 +   x1*x3*(x1-x3)*y2 +   x2*x1*(x2-x1)*y3);
	
	xMin = -b/(2*a);
	yMin = -b*b/(4*a) + c;
	
	return 0;
}


int minFromDataParabola(vector<int> x, vector<double> y, double &xMin, 
                        double &yMin)
{
	vector<double> xd(x.size(),0);
	for (int i=0; i<x.size(); i++) xd[i] = x[i];
	
	return minFromDataParabola(xd,y,xMin,yMin);
}
