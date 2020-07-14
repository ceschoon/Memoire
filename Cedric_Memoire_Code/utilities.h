#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std; 

// Read data from file
// Expected format: parameter = value, all spaces and tabs are ignored

int readDataFromFile(ifstream &file, string parameter, double &value);
int readDataFromFile(ifstream &file, string parameter, int &value);

// Read column vector from file
// Expected format:

//# comments
//value11 value21
//value12 value22
//value13 value23
//value14 value24
//value15 value25
//value16 value26
//value17 value27

// columnNumber starts at 0

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<double> &vec);
void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<int> &vec);

// Evaluate function sampled by a few data points

int evalFromDataInterpolation(vector<double> x, vector<double> y, double xIn,
                              double &yOut);
int evalFromDataInterpolation(vector<int> x   , vector<double> y, double xIn,
                              double &yOut);

// Find zero of function sampled by a few data points
// The function assumes that a zero exists between the min and max values
// of the x vector.

int zeroFromDataInterpolation(vector<double> x, vector<double> y, double &root);

// TODO: Can fail because the cubic can have multiple minima !!??
// Find minimum of function sampled by a few data points

int minFromDataCSpline(vector<double> x, vector<double> y, double &xMin,
                       double &yMin);
int minFromDataCSpline(vector<int> x   , vector<double> y, double &xMin,
                       double &yMin);

// Find minimum of function sampled by a few data points

int minFromDataParabola(vector<double> x, vector<double> y, double &xMin, 
                        double &yMin);
int minFromDataParabola(vector<int> x   , vector<double> y, double &xMin, 
                        double &yMin);

// Evaluate function sampled by a few data points (using Parabola)

int evalFromDataParabola(vector<double> x, vector<double> y, double xIn,
                              double &yOut);
int evalFromDataParabola(vector<int> x, vector<double> y, double xIn,
                              double &yOut);

// Estimation od the error in parabolic interpolation

int errorInDataParabola(vector<double> x, vector<double> y, double xIn,
                        double &yErr);
int errorInDataParabola(vector<int> x, vector<double> y, double xIn,
                        double &yErr);



