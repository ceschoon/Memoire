#include "Log.h"


using namespace std;


void minOverNpoints( int argc, char** argv, Log &log, double kT, double mu, 
                     double &aLattice_min, double &density_min, 
                     double &numParticles_min, double &freeEnergy_min,
                     double &aVdW, double &hsd, bool &success );

void DFTcomputation( int argc, char** argv, Log &log, 
                     double kT, double mu, int Npoints, 
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, double &aVdW, double &hsd, 
                     bool &success );

void DFTcomputation2( int argc, char** argv, Log &log, 
                      double kT, double mu, int Npoints,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, double &aVdW, double &hsd, 
                      bool &success );