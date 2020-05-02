#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int planeColorFCC(int plane)
{
	while (plane < 0) plane += 3;
	
	if (plane%3==0) return 6;
	if (plane%3==1) return 7;
	if (plane%3==2) return 2;
	
	return 1;
}

int planeColorHCP(int plane)
{
	while (plane < 0) plane += 2;
	
	if (plane%2==0) return 6;
	if (plane%2==1) return 7;
	
	return 1;
}

int main()
{
	// Generate FCC lattice
	
	ofstream dataFile;
	dataFile.open("FCC.dat");
	
	vector<int> normal = {1,1,1};
	
	vector<vector<int>> basis;
	basis.push_back({0,0,0});
	basis.push_back({1,1,0});
	basis.push_back({1,0,1});
	basis.push_back({0,1,1});
	
	for (int ix=-1; ix<=1; ix++)
	for (int iy=-1; iy<=1; iy++)
	for (int iz=-1; iz<=1; iz++)
	{
		dataFile << endl;
		
		for (vector<int> point : basis)
		{
			point[0] += ix*2;
			point[1] += iy*2;
			point[2] += iz*2;
			
			int plane = point[0]*normal[0] + point[1]*normal[1] + point[2]*normal[2];
			plane /= 2;
			
			dataFile << double(point[0])/2 << " "
			         << double(point[1])/2 << " "
			         << double(point[2])/2 << " "
			         << planeColorFCC(plane)
			         << endl;
		}
	}
	
	dataFile.close();
	
	
	
	// Generate HCP lattice
	
	dataFile.open("HCP.dat");
	
	normal = {1,-1,0};
	
	basis.clear();
	basis.push_back({0,0,0});
	basis.push_back({1,1,4});
	basis.push_back({1,4,1});
	basis.push_back({4,1,1});
	basis.push_back({5,5,2});
	basis.push_back({5,2,5});
	basis.push_back({2,5,5});
	
	for (int ix=-1; ix<=1; ix++)
	for (int iy=-1; iy<=1; iy++)
	for (int iz=-1; iz<=1; iz++)
	{
		dataFile << endl;
		
		for (vector<int> point : basis)
		{
			point[0] += ix*6;
			point[1] += iy*6;
			point[2] += iz*6;
			
			int plane = point[0]*normal[0] + point[1]*normal[1] + point[2]*normal[2];
			plane /= 3;
			
			dataFile << double(point[0])/6 << " "
			         << double(point[1])/6 << " "
			         << double(point[2])/6 << " "
			         << planeColorHCP(plane)
			         << endl;
		}
	}
	
	dataFile.close();
	
	return 0;
}