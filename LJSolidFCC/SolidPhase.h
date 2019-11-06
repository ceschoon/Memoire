#include "Log.h"


using namespace std;


void minOverNpoints( int argc, char** argv, Log &log, double kT, double mu, 
                     double &aLattice_min, double &density_min, 
                     double &numParticles_min, double &freeEnergy_min,
                     bool &success);

/*
void minOverAlpha( int argc, char** argv, Log &log,  double kT, double mu, 
                   int Npoints, double &alpha_min, double &density_min, 
                   double &numParticles_min, double &freeEnergy_min,
                   bool &success);
*/

void DFTcomputation( int argc, char** argv, Log &log,  double kT, double mu, 
                     int Npoints, bool runMinimizer,
                     double &density_final, double &numParticles_final,
                     double &freeEnergy_final, bool &success            );

void DFTcomputation2( int argc, char** argv, Log &log, double kT, double mu, 
                      int Npoints, bool runMinimizer,
                      double &density_final, double &numParticles_final,
                      double &freeEnergy_final, bool &success            );