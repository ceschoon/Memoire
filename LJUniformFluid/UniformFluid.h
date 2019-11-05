// Functions related to the uniform fluid

using namespace std;


double uniformFluidOmega(double kT, double mu, double aVdW, 
                         double hsd, double rho);

double uniformFluidOmegaDerivative(double kT, double mu, double aVdW, 
                                   double hsd, double rho);

double uniformFluidOmegaDerivative2(double kT, double mu, double aVdW, 
                                    double hsd, double rho);

double muFromCoexDensity(double rho, double kT, double aVdW, double hsd);


void coexistenceDensities(int argc, char **argv, Log &log,
                          double kT, double aVdW, double hsd,
                          double &rho1, double &rho2, double &muCoex,
                          bool &critical, bool &success);

/*
double findRoot(double rhoInit, double kT, double mu, double aVdW, 
                double hsd , bool &success);

double coexistenceMu(int argc, char **argv, double kT, double aVdW, 
                     double hsd , bool &success);
*/