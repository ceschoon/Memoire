#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;

const double PI = 3.14159265358979;

int main()
{
	////////////////////////// Generate HCP Cell ///////////////////////////
	
	// Output file
	
	ofstream dataFile;
	dataFile.open("FCC_cell.dat");
	
	// Cell parameters
	
	double a = 1.0;
	
	// Colour code
	
	int colourback = 8;
	int colourhighlight = 2;
	
	// Generate lower square
	
	double xc = a/2;
	double yc = a/2;
	double zc = 0;
	
	dataFile << endl;
	dataFile << xc << " " << yc << " " << zc << " " << colourhighlight << endl;
	
	double angle = 5*PI/4;
	double radius = a*sqrt(2)/2;
	
	for (int i=0; i<4; i++)
	{
		int colour = colourback;
		if (i==0) colour = colourhighlight;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += PI/2;
	}
	
	// Generate middle square
	
	xc = a/2;
	yc = a/2;
	zc = a/2;
	
	dataFile << endl;
	
	angle = PI;
	radius = a/2;
	
	for (int i=0; i<4; i++)
	{
		int colour = colourback;
		if (i==0 || i==1) colour = colourhighlight;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += PI/2;
	}
	
	// Generate upper square
	
	xc = a/2;
	yc = a/2;
	zc = a;
	
	dataFile << endl;
	dataFile << xc << " " << yc << " " << zc << " " << colourback << endl;
	
	angle = 5*PI/4;
	radius = a*sqrt(2)/2;
	
	for (int i=0; i<4; i++)
	{
		int colour = colourback;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += PI/2;
	}
	
	return 0;
}