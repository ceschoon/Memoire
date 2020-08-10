#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;




int main(int argc, char** argv)
{
	// Check for correct number of arguments
	if (argc != 3)
	{
		cout << "Adds a given number of spaces between to consecutive tabs" << endl;
		cout << "Usage: ./align_correctly <file> <numSpaces>" << endl;
		return 1;
	}
	
	// Read Arguments
	
	ifstream inFile(argv[1]);
	ofstream outFile(string(argv[1])+".out");
	int numSpaces = stoi(argv[2]);
	
	// Output
	
	string lineIn;
	string lineOut;
	
	while (getline(inFile, lineIn))
	{
		lineOut = "";
		
		// Add spaces if line starts with a tab
		
		if (lineIn[0]=='	') 
		{
			for (int j=0; j<numSpaces; j++) lineOut += ' ';
		}
		
		// Scan the line
		
		for (int i=0; i<lineIn.size()-1; i++)
		{
			lineOut += lineIn[i];
			
			// Add spaces if the current character is a tab and so is the 
			// next one
			
			if (lineIn[i]=='	' && lineIn[i+1]=='	')
			{
				for (int j=0; j<numSpaces; j++) lineOut += ' ';
			}
		}
		
		lineOut += lineIn[lineIn.size()-1];
		
		outFile << lineOut << endl;
	}
	
	return 0;
}