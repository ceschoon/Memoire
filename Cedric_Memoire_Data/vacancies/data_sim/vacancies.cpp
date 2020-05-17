#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
	// Check for correct usage
	
	if (argc != 2)
	{
		cout << "Usage: ./vacancies dataFile" << endl;
		return 1;
	}
	
	// Containers for extracted data
	
	vector<double> Y;
	vector<double> Ylnphi;
	vector<double> v2;
	vector<double> phi;
	vector<double> errabs_phi;
	vector<double> T;
	vector<double> errabs_T;
	
	// Hardcoded errors on Y and Ylnphi
	
	double errabs_Y = 4e-3;
	double errabs_Ylnphi = 5e-2;
	
	// Read data from file
	
	string fileName = argv[1];
	ifstream dataFile(fileName);
	if (!dataFile) {cout << "Cannot read file!" << endl; return 1;}
	else {cout << "Reading file \""<< fileName << "\"" << endl;}
	
	string line;
	while (getline(dataFile, line))
	{
		if (line.size()>0 && line[0]!='#')
		{
			int pos_delimitor1 = 0; bool foundDel1 = false;
			int pos_delimitor2 = 0; bool foundDel2 = false;
			
			for (int i=0; i<line.size(); i++)
			{
				if (line[i]==' ' && !foundDel1) {pos_delimitor1 = i; foundDel1 = true;}
				else if (line[i]==' ' && !foundDel2) {pos_delimitor2 = i; foundDel2 = true;}
			}
			
			string Y_str = line.substr(0,pos_delimitor1);
			string Ylnphi_str = line.substr(pos_delimitor1+1, 
			                                pos_delimitor2-pos_delimitor1-1);
			string v2_str = line.substr(pos_delimitor2+1, 
			                            line.size()-pos_delimitor2-1);
			cout << "Read Y_str = \"" << Y_str 
			     << "\" Ylnphi_str = \"" << Ylnphi_str 
			     << "\" v2_str = \"" << v2_str 
			     << "\"" << endl;
			
			Y.push_back(stod(Y_str));
			Ylnphi.push_back(stod(Ylnphi_str));
			v2.push_back(stod(v2_str));
		}
	}
	
	// Compute phi and error
	
	phi = vector<double>(Y.size(),0);
	errabs_phi = vector<double>(Y.size(),0);
	T = vector<double>(Y.size(),0);
	errabs_T = vector<double>(Y.size(),0);
	
	for (int i=0; i<Y.size(); i++)
	{
		double errrel_Y = abs(errabs_Y / Y[i]);
		double errrel_Ylnphi = abs(errabs_Ylnphi / Ylnphi[i]);
		double ratio = Ylnphi[i]/Y[i];
		double errrel_ratio = sqrt(errrel_Y * errrel_Y
		                          +errrel_Ylnphi * errrel_Ylnphi);
		
		phi[i] = exp(Ylnphi[i]/Y[i]);
		errabs_phi[i] = phi[i] * abs(ratio) * errrel_ratio; // propagation of std errors
		//errabs_phi[i] = exp((Ylnphi[i]+errabs_Ylnphi)/(Y[i]+errabs_Y)) // max deviation possible
		//               -exp(Ylnphi[i]/Y[i]);
		
		T[i] = 4*Y[i]/(v2[i]*v2[i]);
		errabs_T[i] = 4*errabs_Y/(v2[i]*v2[i]);
		
		cout << scientific << setprecision(2);
		cout << "Computed phi = " << phi[i] << " +- " << errabs_phi[i] << endl;
		cout << "           T = " << T[i] << " +- " << errabs_T[i] << endl;
	}
	
	// Store computed quantities in a new file
	
	ofstream outFile("out_"+fileName);
	
	// copy the comment header from the original file
	dataFile.clear();
	dataFile.seekg(0, ios_base::beg);
	while (getline(dataFile, line))
	{
		if (line.size()>0 && line[0]=='#') {outFile << line << endl;}
		else {break;}
	}
	
	// store new data
	outFile << "#Y        Ylnphi   v2       T        errabs_T phi      errabs_phi" << endl;
	outFile << scientific << setprecision(2);
	for (int i=0; i<Y.size(); i++)
	{
		outFile << Y[i] << " "
		        << Ylnphi[i] << " "
		        << v2[i] << " "
		        << T[i] << " "
		        << errabs_T[i] << " "
		        << phi[i] << " "
		        << errabs_phi[i] << " "
		        << endl;
	}
	
	return 0;
}