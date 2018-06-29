#ifndef __LUTSKO_STRING_METHOD__
#define __LUTSKO_STRING_METHOD__

#include "Minimizer.h"
#include "Grace.h"

class StringMethod_MPI
{
 public:
 StringMethod_MPI(double mu, bool freeEnd) : mu_(mu), freeEnd_(freeEnd), step_counter_(0){};
  ~StringMethod_MPI(){};

  void setMu(double m) { mu_ = m;}
  virtual void run(string& logfile) = 0;

  
 protected:
  double mu_;
  bool freeEnd_;
  long step_counter_;
};


class StringMethod_MPI_Master : public StringMethod_MPI
{
 public:
 StringMethod_MPI_Master(int Nimages, Density &finalDensity, double bav, double F_final, double mu, Grace *g = NULL, bool freeEnd = false) 
   : StringMethod_MPI(mu, freeEnd), finalDensity_(finalDensity), bav_(bav), grace_(g)
    {
      Ntot_ = finalDensity.Ntot();
      gr_ = new mglGraph;
      Images_.resize(Nimages);
      dF_.resize(Nimages);
      dF_[Nimages-1] = F_final;

      
      //We just need one density to set up the array which is then shared by all
      finalDensity.initialize_2D_data(data_2D_);
  }
  ~StringMethod_MPI_Master(){if(gr_) delete gr_;}

  virtual void run(string& logfile);
  void interpolate();
  void processImages();
  void report(string &logfile);
  
  void Display(int dataSet, double dFmax = 0.0, double dFav = 0.0);
  void Draw(vector<double> &data, int image_number, double F);
  
  void addTask(int images) { taskList.push_back(images);}
  void archive(string &filename) const;
  
 private:
    Grace *grace_;
    vector<int> taskList;

    Density &finalDensity_;
    double bav_; // background density so we can reconstruct initial state
    
    vector< vector<double> > Images_;
    vector<double> dF_;
    long Ntot_;

    mglGraph *gr_;
    mglData data_2D_;

};


class StringMethod_MPI_Slave : public StringMethod_MPI
{
 public:
 StringMethod_MPI_Slave(DDFT &ddft, vector<Density*> string, double mu, int id, int offset) 
   : StringMethod_MPI(mu, false), ddft_(ddft), string_(string), id_(id), offset_(offset)
  {
    int N = string_.size();  

    DT_.resize(N,0.0);
    oldF_.resize(N,0.0);
    N_.resize(N,0.0);
    distances_.resize(N,0.0);

    // Initialize free energies  
    for(int J=0;J<N;J++)
      {
	oldF_[J] = ddft_.F_string(*(string_[J])) - mu_*string_[J]->getNumberAtoms();
	DT_[J] = 0.01; 
      }

    // this holds a copy of the string
    string_copy_.resize(N);
    for(int J=0;J<N;J++)
      string_copy_[J].zeros(string_[0]->Ntot());
    
  }
  
  ~StringMethod_MPI_Slave(){}

  virtual void run(string& logfile);
  virtual void report();
  
  void archive(string &filename) const;
  void read(string &filename);
  
 private:
  void rescale(double Ntarget);
  void rescale_linear();
  void Draw(Density& density, int image_number, double F);
  
 private:
  DDFT &ddft_;
  vector<Density*> string_;
  vector<DFT_Vec> string_copy_;   // this holds a copy of the string

  vector<double> DT_;
  vector<double> oldF_;
  vector<double> N_;
  vector<double> distances_;

  double dFmax_;
  double dFav_;
  double delta_max_;
  
  int id_;
  int offset_; // this tells us which are the real image indexes
};



#endif //#ifndef __LUTSKO_STRING_METHOD__
