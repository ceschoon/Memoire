#include <iostream>
#include <vector>
#include "../utilities.h"
#include "plotting.h"

using namespace std;

int main()
{
	// create data set
	/*
	int N = 10;
	
	vector<double> x(N,0);
	vector<double> y(N,0);
	
	x = {    0,    1,    2,    3,    4,    5,    6,    7,    8,    9};
	y = { 24.2, 15.9, 8.55, 3.89, 1.09, 0.03, 1.22, 4.05, 9.35, 16.8};
	*/
	
	int N=3;
	
	vector<double> x(N,0);
	vector<double> y(N,0);
	
	x = {-2.2,-1.1, 0.3};
	y = {1.3,-0.05,1.2};
	
	// plot data set
	
	plotInFile(x,y,"dataset");
	
	// find zero
	
	double xMin = 0;
	double yMin = 0;
	
	int status = minFromDataParabola(x,y,xMin,yMin);
	
	cout << "xMin = " << xMin << endl;
	cout << "yMin = " << yMin << endl;
	cout << "status = " << status << endl;
	
	return 0;
}