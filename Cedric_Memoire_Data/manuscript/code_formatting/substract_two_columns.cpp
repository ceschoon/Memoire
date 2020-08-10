#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

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





int main(int argc, char** argv)
{
	// Check for correct number of arguments
	if (argc != 3)
	{
		cout << "Usage: ./extract_one_column <file1> <file2>" << endl;
		return 1;
	}
	
	// Read Arguments
	
	ifstream inFile1(argv[1]);
	ifstream inFile2(argv[2]);
	
	// Output
	
	string line1;
	string line2;
	
	while (getline(inFile1, line1))
	{
		// Note: It is assumed that both files have the same number of lines
		getline(inFile2, line2);
		
		// skip but reproduce empty lines
		if (line1.size()==0) 
		{
			cout << endl;
			continue;
		}
		
		string value1_str = getFieldXFromLine(line1, 0);
		string value2_str = getFieldXFromLine(line2, 0);
		
		double value1 = stod(value1_str);
		double value2 = stod(value2_str);
		
		cout << value1-value2 << endl;
		
	}
	
	return 0;
}