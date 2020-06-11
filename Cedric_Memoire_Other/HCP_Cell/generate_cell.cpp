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
	dataFile.open("HCP_cell.dat");
	
	// Cell parameters
	
	double a = 1;
	double b = a*sqrt(3);
	double c = a*sqrt(8.0/3);
	
	// Colour code
	
	int colourback = 8;
	int colourhighlight = 2;
	
	// Generate lower hexagon
	
	double xc = a/2;
	double yc = b/2;
	double zc = 0;
	
	dataFile << endl;
	dataFile << xc << " " << yc << " " << zc << " " << colourhighlight << endl;
	
	double angle = PI/3*4;
	double radius = a;
	
	for (int i=0; i<6; i++)
	{
		int colour = colourback;
		if (i==0) colour = colourhighlight;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += PI/3;
	}
	
	// Generate upper hexagon
	
	xc = a/2;
	yc = b/2;
	zc = c;
	
	dataFile << endl;
	dataFile << xc << " " << yc << " " << zc << " " << colourback << endl;
	
	angle = PI/3*4;
	radius = a;
	
	for (int i=0; i<6; i++)
	{
		int colour = colourback;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += PI/3;
	}
	
	// Generate middle triangle
	
	xc = a/2;
	yc = b/2;
	zc = c/2;
	
	dataFile << endl;
	
	angle = PI/6*7;
	radius = a*(sqrt(3)/2-1.0/4);
	
	for (int i=0; i<3; i++)
	{
		int colour = colourback;
		if (i==0 || i==2) colour = colourhighlight;
		
		dataFile << xc+radius*cos(angle) << " " << yc+radius*sin(angle) 
		         << " " << zc << " " << colour << endl;
		
		angle += 2*PI/3;
	}
	
	return 0;
}