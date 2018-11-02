#ifndef __LUTSKO_CONJUG2__
#define __LUTSKO_CONJUG2__

#include <mgl2/mgl.h>

using namespace std;

#include "Grace.h"

#include "DFT.h"

/**
  *  @brief Minimizer base class
  *
  */  

class Minimizer
{
 public:
 Minimizer(DFT &dft, Density &density, double mu) : dft_(dft), density_(density), mu_(mu), forceLimit_(0.1), err_(0.0), bFrozenBoundary_(false)
  {
    x_.resize(density_.Ntot());
    dF_.resize(x_.size());
  }

  void setFrozenBoundaryFlag(bool f) {bFrozenBoundary_ = f;}
  
  virtual void initialize();

  void run(string& logfile, long maxSteps = -1);

  virtual double step() = 0;
  virtual void finish(const char *) = 0;

  virtual int draw_before() // Display something before the next step
  {
    stringstream ts;
    ts << "calls = " << calls_ << " F = " << F_ << " N = " << density_.getNumberAtoms();

    string title = ts.str();
    string file("image_current.png");
    density_.doDisplay(title, file);

  }
  virtual int draw_during() = 0;  // Display something during the minimization
  virtual void draw_after()  = 0;  // Display something after the minimization

  const Density & getDensity() const { return density_;}
  int getCalls() const { return calls_;}
  double getF() const { return F_;}
  double getErr() const { return err_;}

  double getForceTerminationCriterion() const {return forceLimit_;}
  void   setForceTerminationCriterion(double v) {forceLimit_ = v;}

  Density & getDensity(){return density_;}
  
  virtual double getDF_DX();

  virtual double get_convergence_monitor() const { return dF_.inf_norm()/density_.dV();}

  
 protected:
  DFT &dft_;
  Density &density_;
  double mu_;
  double err_;

  // A hook to allow for graphical displays
  void (*display_)(Minimizer &minimizer);
  
  // Working space for the minimization

  DFT_Vec dF_;
  DFT_Vec x_;

  int calls_ = 0;
  int step_counter_ = 0;
  double F_ = 0;

  double forceLimit_;
  double f_abs_max_; // max absolute value of dF_
  bool bFrozenBoundary_;
};

/**
  *  @brief Minimizer_Fixed_N: specializes to the case of fixed number of particles
  *
  */  

class Minimizer_Fixed_N : public Minimizer
{
 public:
 Minimizer_Fixed_N(DFT &dft, Density &density, double Nfixed) : Minimizer(dft, density, 0.0), N_fixed_target_(Nfixed)
  {}

  virtual void initialize();
  
  virtual double getDF_DX();

 protected:
  double mu_eff_;
  double err_;

  double N_fixed_target_;
};


/**
  *  @brief Conjugate Gradients Class: a hand-coded CG minimizer with backtracking.
  */  

class ConjugateGradients2 : public Minimizer
{
 public:
 ConjugateGradients2(DFT &dft, Density &density, double mu) : Minimizer(dft, density, mu) , Mixing_(1.0){}
  
  virtual void initialize();
  
  virtual double step();
  double tryStep(double alpha, double &dalf, DFT_Vec &y, double &df);

  void check(bool write);

  virtual double getDF_DX();

  virtual int draw_during();  // Display something during the minimization
  virtual void draw_after();  // Display something after the minimization

  double get_sigma_conj_grad_secant() const { return sigma_conj_grad_secant_;}
  void   set_sigma_conj_grad_secant(double v) { sigma_conj_grad_secant_ = v;}

  double get_df_limit_conj_grad() const { return df_limit_conj_grad_;}
  void   set_df_limit_conj_grad(double v) { df_limit_conj_grad_= v;}

  int  get_n_reset_conj_grad() const { return n_reset_conj_grad_;}
  void set_n_reset_conj_grad(double v) { n_reset_conj_grad_ = v;}

  double get_mixing_parameter() const {return Mixing_;}
  void   set_mixing_parameter(double v) {Mixing_ = v;}

 protected:

  // Working space for the minimization

  DFT_Vec r_;
  DFT_Vec d_;

  double delta_new_ = 0.0;
  double delta_0_   = 0.0;

  // constants that control the minimization

  double sigma_conj_grad_secant_ = 1;

  int jmax_conj_grad_secant_ = 40;
  int n_reset_conj_grad_ = 100;
  double df_limit_conj_grad_ = 1e-10;

  double Mixing_;
};


class ConjugateGradients2_variable : public ConjugateGradients2
{
 public:
 ConjugateGradients2_variable(DFT &dft, Density &density, double mu, bool showGraphics = true) : ConjugateGradients2(dft, density, mu)
  {}

  //  virtual int draw_during();
  //  virtual void draw_after();

  virtual void finish(const char *c)
  {
  }

 protected:
};

// Here we minimize  F[rho0 + x*x*a] over x with a chosen to fix the total mass.
// There are a couple of problems.
// Since the total mass is fixed, there are really only Ntot-1 independent variables. We should 
// really reduce the size of the variable array by one to account for this ... 

class ConjugateGradients2_Fixed : public ConjugateGradients2
{
 public:
 ConjugateGradients2_Fixed(DFT &dft, Density &density, double Ntarget, bool showGraphics = true) : ConjugateGradients2(dft, density, 0.0), N_fixed_target_(Ntarget)
  {}

  virtual int draw_during();
  virtual void draw_after();

  virtual double getDF_DX();

  virtual void finish(const char *c)
  {
  }

 protected:
  double N_fixed_target_;
  double mu_eff_;
};


/**
  *  @brief Picard minimizer Class
  *
  */  

class Picard : public Minimizer
{
 public:
 Picard(DFT &dft, Density &density, double mu) : Minimizer(dft, density, mu) , Mixing_(1.0){}
  
  virtual void initialize();
  
  virtual double step();

  virtual int draw_during(){ return 1;}  // Display something during the minimization
  virtual void draw_after()  // Display something after the minimization
  {
    cout << "After picard step " << step_counter_ 
	 << " F-mu*N = " << F_ 
	 << " max(dF) = " << f_abs_max_
	 << " and N = " << density_.getNumberAtoms() 
	 << endl;
  }

  virtual void finish(const char *){}

  double get_mixing_parameter() const {return Mixing_;}
  void   set_mixing_parameter(double v) {Mixing_ = v;}

 protected:

  double Mixing_;
};




/**
  *  @brief DDFT minimizer Class 
  *
  */  

class DDFT : public Minimizer
{
 public:
 DDFT(DFT &dft, Density &density, Grace *g = NULL, bool showGraphics = true)
   : Minimizer(dft, density, 0.0), show_(showGraphics) ,grace_(g), tolerence_fixed_point_(1e-4), successes_(0), fixedBorder_(false), modified_(false)
  {
    dt_ = 10*0.1*density_.getDX() * density_.getDX();
    dt_ = 0.0001*density_.getDX() * density_.getDX();
    dtMax_ = density_.getDX() * density_.getDX();
  }
  ~DDFT() {}

  virtual void initialize();
  
  //  virtual int draw_before(); // Display something before the next step
  virtual int draw_during(){}  // Display something during the minimization
  virtual void draw_after()  // Display something after the minimization
  {
    cout << "After DDFT step " << step_counter_ 
	 << " F = " << F_ 
	 << " max(dF) = " << f_abs_max_
	 << " and N = " << density_.getNumberAtoms() 
	 << endl;
  }
  virtual void finish(const char *c){};

  void Display(double F, double dFmin, double dFmax, double N);
    
  void set_tolerence_fixed_point(double e) { tolerence_fixed_point_ = e;}
  void set_max_time_step(double t) { dtMax_ = t;}
  void setTimeStep(double dt) { dt_ = dt;}

  void setFixedBoundary() {fixedBorder_ = true;}
  void setModified() {modified_ = true;}
    
  double getTimeStep() const { return dt_;}
  
  double F_string(Density &d, double *fmax = NULL);
  void reverseForce(DFT_Vec *tangent);

  virtual double step() = 0;
  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true) = 0;

 protected:

  bool show_;
  
  Grace *grace_;
  double dt_;
  DFT_Vec oldF_;
  double tolerence_fixed_point_;

  // control of adaptive time step
  int successes_;
  double dtMax_;

  bool fixedBorder_;
  bool modified_; // for certain derived classes
};


/**
  *  @brief DDFT minimizer Class using discretized equations
  *
  */  

class DDFT_Discrete : public DDFT
{
 public:
 DDFT_Discrete(DFT &dft, Density &density, bool bFixedBoundaries = false, Grace *g = NULL, bool showGraphics = true)
   : DDFT(dft, density, g, showGraphics), bFixedBoundaries_(bFixedBoundaries)
  {}
  ~DDFT_Discrete() {}

  virtual void initialize();
  
  void sub_step_x(DFT_Vec &new_density, const Density &original_density, bool bFixedBoundaries = false);
  void sub_step_y(DFT_Vec &new_density, const Density &original_density, bool bFixedBoundaries = false);
  bool sub_step_z(DFT_Vec &new_density, const Density &origina_density,  bool bSelfConsistent = false, bool bFixedBoundaries = false);

  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  void   solv_tridiag(const DFT_Vec &b, DFT_Vec &RHS, double D, bool bFixedBoundaries = false);
  void   solv_periodic_tridiag(DFT_Vec &RHS, double D);
  void   solv_periodic_tridiag_2(DFT_Vec &b, DFT_Vec &RHS, double D);

  void test_solv_tridiag();
  
 protected:
  bool bFixedBoundaries_;
};

/**
  *  @brief DDFT minimizer Class using integrating factor
  *
  */  

class DDFT_IF : public DDFT
{
 public:
 DDFT_IF(DFT &dft, Density &density, Grace *g = NULL, bool showGraphics = true)
   : DDFT(dft, density, g, showGraphics)
    {}
  ~DDFT_IF() {}

  virtual void initialize();

  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  double fftDiffusion(const Density &density, DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1);
  void calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1);
  void restore_values_on_border(DFT_Vec& d1, const DFT_Vec &d0);

  virtual double get_convergence_monitor() const { return fabs(dF_.max() - dF_.min())/density_.dV();}
 protected:

  double largest_change_on_border_; // for reporting effect of restore_values_on_border()
};

/**
  *  @brief DDFT minimizer Class using integrating factor
  *
  */  

class DDFT_IF_Open : public DDFT
{
 public:
 DDFT_IF_Open(DFT &dft, Density &density, double background, Grace *g = NULL, bool showGraphics = true)
   : DDFT(dft, density, g, showGraphics), background_(background), sin_in_(NULL), sin_out_(NULL)
    {
    }
  ~DDFT_IF_Open() {if(sin_in_) delete sin_in_; if(sin_out_) delete sin_out_;}

  virtual void initialize();

  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  double fftDiffusion(const Density &density, DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1) {throw std::runtime_error("Need to adapt fftDiffusion for non-string application");}
  double fftDiffusion(DFT_Vec &d1, const double *RHS0_sin_transform, const double *RHS1_sin_transform);
  void calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1);

  void pack_for_sin_transform(const double *x, double val);
  
  void unpack_after_transform(double *x, double val);

  
 protected:
  double background_;
  
  fftw_plan sin_plan_;
  unsigned sin_Ntot_;
  unsigned sin_Norm_;
  double *sin_in_;
  double *sin_out_;

  vector<double> Lamx;
  vector<double> Lamy;
  vector<double> Lamz;
};

class DDFT_Open : public DDFT
{
 public:
 DDFT_Open(DFT &dft, Density &density, Grace *g = NULL, bool showGraphics = true)
   : DDFT(dft, density, g, showGraphics)
    {}
  ~DDFT_Open() {}

  virtual void initialize();
  
  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  double fftDiffusion(const Density &density,DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1);
  void calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1);

 protected:
};

/**
  *  @brief Minimizer using nlOpt library
  *
  */
/*

class nlOptMinimizer : public Minimizer
{
 public:
 nlOptMinimizer(DFT &dft, Density &density, double Ntarget) :  Minimizer(dft, density, 0.0), N_fixed_target_(Ntarget) 
  {
  }
  
  //  virtual void initialize();

  void run(string& logfile, long maxSteps = -1);

  virtual double step() { return 1;}
  virtual void finish(const char *) {}

  virtual int draw_during() {return 1;}  // Display something during the minimization
  virtual void draw_after() {}  // Display something after the minimization

  double objectiveFunction(unsigned n, const double *x, double *grad);
  
 protected:
  double N_fixed_target_;
};
*/

/**
  *  @brief Minimizer using FIRE algorithm
  *
  */

class fireMinimizer : public Minimizer_Fixed_N
{
 public:
 fireMinimizer(DFT &dft, Density &density, double Ntarget) :  Minimizer_Fixed_N(dft, density, Ntarget)
  { 
    v_.resize(density_.Ntot());
    //dt_ = 1e-3; // my guess
    dt_ = 1e-3; // my guess
    dt_max_ = 10*dt_;
    
    N_min_ = 5;
  
    alpha_start_ = 0.1;
    f_dec_ = 0.5;
    f_inc_ = 1.1;
    f_alf_ = 0.99;
  }
  
  virtual void initialize();

  virtual double step();
  virtual void finish(const char *) {}

  virtual int draw_during();
  virtual void draw_after();

  void verlet();

  void setTimeStep(double dt) { dt_ = dt;}
  void setTimeStepMax(double dt) { dt_max_ = dt;}
  
  void setAlphaStart(double a) { alpha_start_ = a;}
  void setAlphaFac(double a) { f_alf_ = a;}

 protected:
  DFT_Vec v_;


  double alpha_start_;
  double f_dec_;
  double f_inc_;
  double f_alf_;
  double dt_max_;
  double dt_;

  unsigned it_;
  unsigned cut_;
  int N_min_;
  double alpha_;

};

/**
  *  @brief Minimizer using FIRE algorithm
  *
  */

class fireMinimizer_Mu : public Minimizer
{
 public:
 fireMinimizer_Mu(DFT &dft, Density &density, double mu) :  Minimizer(dft, density, mu)
  { 
    v_.resize(density_.Ntot());
    //dt_ = 1e-3; // my guess
    dt_ = 1e-3; // my guess
    dt_max_ = 10*dt_;
    
    N_min_ = 5;
  
    alpha_start_ = 0.1;
    f_dec_ = 0.5;
    f_inc_ = 1.1;
    f_alf_ = 0.99;
  }
  
  virtual void initialize();

  virtual double step();
  virtual void finish(const char *) {}

  virtual int draw_during();
  virtual void draw_after();

  void verlet();

  void setTimeStep(double dt) { dt_ = dt;}
  void setTimeStepMax(double dt) { dt_max_ = dt;}
  
  void setAlphaStart(double a) { alpha_start_ = a;}
  void setAlphaFac(double a) { f_alf_ = a;}

 protected:
  DFT_Vec v_;


  double alpha_start_;
  double f_dec_;
  double f_inc_;
  double f_alf_;
  double dt_max_;
  double dt_;

  unsigned it_;
  unsigned cut_;
  int N_min_;
  double alpha_;

};



/**
  *  @brief Minimizer using Ceres library
  *
  */  
/*
class ceresOpt : public FirstOrderFunction
{
 public:
 ceresOpt(DFT &dft, Density &density) :  FirstOrderFunction(), dft_(dft), density_(density) {}
  
  virtual ~ceresOpt() {}
  virtual bool Evaluate(const double* const parameters,
                         double* cost,
                         double* gradient) const;
  virtual int NumParameters() const {return density_.Ntot();}

  
 protected:
  DFT &dft_;
  Density &density_;
};


*/
#endif // sentinal
