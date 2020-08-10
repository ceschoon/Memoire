#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;




string getFieldXFromLine(string line, int fieldTarget)
{
	string field_str;
	
	int fieldCounter = 0;
	bool newField = false;
	
	// skip comment line
	if (line[0]=='#') return "";
	
	for (int i=0; i<line.size(); i++)
	{
		// delimiter
		if (line[i]==' ' || line[i] == '	') newField = true;
		
		// field character
		else
		{
			// increment field counter if necessary
			if (newField) 
			{
				fieldCounter++;
				newField = false;
			}
			else ;
			
			// copy field value if we are at the target field
			if (fieldCounter == fieldTarget) field_str += line[i];
		}
	}
	
	return field_str;
}



double mean(vector<double> values)
{
	double sum = 0;
	for (int i=0; i<values.size(); i++) sum += values[i];
	
	return sum/values.size();
}



int main(int argc, char** argv)
{
	// Check for correct number of arguments
	if (argc != 2)
	{
		cout << "Usage: ./compute_mean <file>" << endl;
		return 1;
	}
	
	// Read Arguments
	
	ifstream inFile(argv[1]);
	
	// Read values in file
	
	string line;
	vector<double> values;
	
	while (getline(inFile, line))
	{
		// skip empty lines
		if (line.size()==0) continue;
		
		string value_str = getFieldXFromLine(line, 2);
		double value = stod(value_str);
		
		values.push_back(value);
	}
	
	// Compute and print mean value
	
	cout << scientific << setprecision(2) << endl;
	cout << "mean = " << mean(values) << endl;
	
	return 0;
}