#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;

const double PI = 3.14159265358979;

void addUnitCell(double posX, double posY, double aCell, int colour, int pointtype, ofstream &dataFile)
{
	double x1 = posX + 0;
	double y1 = posY + 0;
	
	double x2 = posX + aCell/2;
	double y2 = posY + aCell/2;
	
	dataFile << endl;
	dataFile << x1 << " " << y1 << " " << colour << " " << pointtype << endl;
	dataFile << x2 << " " << y2 << " " << colour << " " << pointtype << endl;
}

int main()
{
	////////////////////////// Generate HCP Cell ///////////////////////////
	
	// Output file
	
	ofstream dataFile;
	dataFile.open("lattice.dat");
	
	// Cell parameters
	
	double aCell = 1.0;
	
	// Colour code
	
	int colourback = 8;
	int colourhighlight = 2;
	int pointtypeback = 7;
	int pointtypehighlight = 65;
	
	// Generate 
	
	for (int i=-1; i<5; i++)
	{
		for (int j=-1; j<3; j++)
		{
			double posX = i*aCell;
			double posY = j*aCell;
			
			
			// particular cell containing the vacancy defect
			if (i==2 && j==1)
			{
				double x1 = posX + 0;
				double y1 = posY + 0;
				
				double x2 = posX + aCell/2;
				double y2 = posY + aCell/2;
				
				dataFile << endl;
				dataFile << x1 << " " << y1 << " " << colourhighlight << " " << pointtypehighlight << endl;
				dataFile << x2 << " " << y2 << " " << colourback << " " << pointtypeback << endl;
			}
			// regular cells
			else
			{
				addUnitCell(posX, posY, aCell, colourback, pointtypeback, dataFile);
			}
		}
	}
	
	return 0;
}