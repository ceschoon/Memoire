#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "../utilities.h"

using namespace std; 

double min(double a, double b) {return (a<b)?a:b;}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Program that puts the coexistence curves together in a phase diagram

int main(int argc, char** argv)
{
	////////////////////////////// Parameters //////////////////////////////
	
	// files and directories
	
	string dataUniformCoexDir = "CoexistenceUniform";
	string dataSolidCoexDir = "Plots-coexistence-curve";
	string dataCurvesDir = "Plots-diagram";
	
	int sysresult = system(("mkdir -p "+dataCurvesDir).c_str());
	
	
	////////////////////////////////////////////////////////////////////////
	
	vector<double> rhoV_coexVL;
	vector<double> rhoL_coexVL;
	vector<double> kT_coexVL;
	
	vector<double> rhoV_coexVS;
	vector<double> rhoS_coexVS;
	vector<double> kT_coexVS;
	
	vector<double> rhoL_coexLS;
	vector<double> rhoS_coexLS;
	vector<double> kT_coexLS;
	
	
	///////////////////////////// Read files ///////////////////////////////
	
	
	// read VL coexistence file
	
	ifstream fileVL(dataUniformCoexDir+"/coexistence_uniform.dat");
	
	readColumnVectorFromFile(fileVL, 1, rhoV_coexVL);
	readColumnVectorFromFile(fileVL, 2, rhoL_coexVL);
	readColumnVectorFromFile(fileVL, 0, kT_coexVL);
	
	// read VS coexistence file
	
	ifstream fileVS(dataSolidCoexDir+"/coexistenceVS.dat");
	
	readColumnVectorFromFile(fileVS, 0, rhoV_coexVS);
	readColumnVectorFromFile(fileVS, 1, rhoS_coexVS);
	readColumnVectorFromFile(fileVS, 3, kT_coexVS);
	
	// read LS coexistence file
	
	ifstream fileLS(dataSolidCoexDir+"/coexistenceLS.dat");
	
	readColumnVectorFromFile(fileLS, 0, rhoL_coexLS);
	readColumnVectorFromFile(fileLS, 1, rhoS_coexLS);
	readColumnVectorFromFile(fileLS, 3, kT_coexLS);
	
	//////////////// Find triple and critical properties ///////////////////
	
	// triple point is between the smallest kT in the LS coexistence curve
	// and the largest kT in the VS coex curve
	// assume that the coex curves are sorted in kT 
	
	double kT_triple = (kT_coexLS[0]+kT_coexVS[kT_coexVS.size()-1])/2;
	
	// find index of kT just below triple point
	int index_below_triple = 0;
	for (int i=0; i<kT_coexVL.size()-2; i++)
		if (kT_coexVL[i]<kT_triple && kT_coexVL[i+1]>kT_triple)
			index_below_triple = i;
	
	double rhoV_triple = ( rhoV_coexVL[index_below_triple] +
	                       rhoV_coexVL[index_below_triple+1] )/2;
	
	double rhoL_triple = ( rhoL_coexVL[index_below_triple] +
	                       rhoL_coexVL[index_below_triple+1] )/2;
	
	double rhoS_triple = min( rhoS_coexLS[0] ,
	                          rhoS_coexVS[rhoS_coexVS.size()-1] );
	
	// the critical point is just above the largest kT in the VL coex curve
	
	double kT_critical = kT_coexVL[kT_coexVL.size()-1];
	
	double rho_critical = ( rhoV_coexVL[kT_coexVL.size()-1] +
	                        rhoL_coexVL[kT_coexVL.size()-1]   )/2;
	
	// save in file
	
	ofstream fileProp(dataCurvesDir+"/triple_and_critical_properties.dat");
	
	fileProp << "kT_triple   = " << kT_triple   << endl;
	fileProp << "rhoV_triple = " << rhoV_triple << endl;
	fileProp << "rhoL_triple = " << rhoL_triple << endl;
	fileProp << "rhoS_triple = " << rhoS_triple << endl;
	fileProp << " " << endl;
	fileProp << "kT_critical  = " << kT_critical  << endl;
	fileProp << "rho_critical = " << rho_critical << endl;
	
	fileProp.close();
	
	
	////////////////////// Curve for uniform phases ////////////////////////
	
	
	ofstream fileCurveUniform(dataCurvesDir+"/curve_uniform.dat");
	
	fileCurveUniform << "# rho kT" << endl;
	
	// rhoV from VS coex curve
	for (int i=0; i<rhoV_coexVS.size(); i++)
		fileCurveUniform << rhoV_coexVS[i] << " " << kT_coexVS[i] << endl;
	
	// rhoV from VL coex curve
	// only above the triple point
	for (int i=0; i<rhoV_coexVL.size(); i++)
		if (kT_coexVL[i] > kT_triple)
			fileCurveUniform << rhoV_coexVL[i] << " " << kT_coexVL[i] << endl;
	
	// rhoL from VL coex curve
	// only above the triple point, reversed ordering
	int n = rhoL_coexVL.size();
	for (int i=0; i<rhoL_coexVL.size(); i++)
		if (kT_coexVL[n-i-1] > kT_triple)
			fileCurveUniform << rhoL_coexVL[n-i-1] << " " << kT_coexVL[n-i-1] << endl;
	
	// rhoL from LS coex curve
	for (int i=0; i<rhoL_coexLS.size(); i++)
		fileCurveUniform << rhoL_coexLS[i] << " " << kT_coexLS[i] << endl;
	
	fileCurveUniform.close();
	
	
	//////////////////////// Curve for solid phase /////////////////////////
	
	
	ofstream fileCurveSolid(dataCurvesDir+"/curve_solid.dat");
	
	fileCurveSolid << "# rho kT" << endl;
	
	// rhoS from VS coex curve
	for (int i=0; i<rhoS_coexVS.size(); i++)
		fileCurveSolid << rhoS_coexVS[i] << " " << kT_coexVS[i] << endl;
	
	// rhoS from LS coex curve
	for (int i=0; i<rhoS_coexLS.size(); i++)
		fileCurveSolid << rhoS_coexLS[i] << " " << kT_coexLS[i] << endl;
	
	fileCurveSolid.close();
	
	
	//////////////////////////// Gnuplot script ////////////////////////////
	
	
	ofstream plotFile(dataCurvesDir+"/plot_diagram");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'diagram.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"LJ Phase Diagram\" font \",20\"" << endl;
	plotFile << "set xlabel \"rho\" font \",20\"" << endl;
	plotFile << "set ylabel \"kT\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "plot \"curve_uniform.dat\" with lines linetype 1, \\" << endl;
	plotFile << "     \"curve_solid.dat\"   with lines linetype 1" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	sysresult = system(("cd "+dataCurvesDir+"; gnuplot plot_diagram").c_str());
	
	
	return 0;
}

