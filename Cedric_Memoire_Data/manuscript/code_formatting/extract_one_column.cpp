#include <iostream>
#include <fstream>
#include <iomanip>

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
		cout << "Usage: ./extract_one_column <in file> <column number>" << endl;
		return 1;
	}
	
	// Read Arguments
	
	ifstream inFile(argv[1]);
	int colNum = stoi(argv[2]);
	//string colTitle = argv[3];
	
	// Output
	
	string line;
	
	while (getline(inFile, line))
	{
		// skip comment line
		if (line[0]=='#') continue;
		
		string value_str = getFieldXFromLine(line, colNum);
		cout << value_str << endl;
		
	}
	
	return 0;
}