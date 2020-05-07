#include <iostream>
#include <fstream>
#include <vector>
#include "../utilities.h"

using namespace std;

int main()
{
	ifstream file("coexistenceFS.dat");
	vector<double> muCoex;
	
	readColumnVectorFromFile(file, 2, muCoex);
	
	cout << "muCoex = ";
	for (int i=0; i<muCoex.size(); i++)
		cout << muCoex[i] << " ";
	cout << endl;
	
	return 0;
}