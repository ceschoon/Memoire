#include <iostream>
#include <vector>
#include "../utilities.h"
#include "plotting.h"

using namespace std;

int main()
{
	// create data set
	
	int N = 10;
	
	vector<double> x(N,0);
	vector<double> y(N,0);
	
	x = {    0,    1,    2,    3,    4,    5,    6,    7,    8,    9};
	y = {-2.32,-1.92,-1.55,-1.09,-0.49, 0.03, 0.59, 0.99, 1.35, 1.78};
	
	// plot data set
	
	plotInFile(x,y,"dataset");
	
	// find zero
	
	double root = 0;
	
	int status = zeroFromDataInterpolation(x,y,root);
	
	cout << "status = " << status << endl;
	
	return 0;
}