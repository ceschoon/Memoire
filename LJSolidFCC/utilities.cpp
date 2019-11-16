#include <fstream>
#include <string>

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remove tabs and spaces from string

// Function used by readDataFromFile

void removeTabsAndSpaces(string &str)
{
	int N = str.size();
	string str_new = "";
	
	for (int i=0; i<N; i++)
		if (str[i]!=' ' && str[i]!='	')
			str_new = str_new + str[i];
	
	str = str_new;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Read data from file

// Expected format: parameter = value, all spaces and tabs are ignored

void readDataFromFile(ifstream &file, string parameter, double &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter) 
			value = stod(line.substr(b+1, a-b-1));
	}
}

void readDataFromFile(ifstream &file, string parameter, int &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter) 
			value = stoi(line.substr(b+1, a-b-1));
	}
}