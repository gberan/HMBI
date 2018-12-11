#ifndef _quasiharmonic_h
#define _quasiharmonic_h

#include <iostream>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "monomer.h"
#include "dimer.h"
#include "params.h"
#include "matrix.h"
#include "vector.h" // my vector class, "Vector"
#include <vector> // C STL library vector, "vector"
using std::vector;
#include <time.h>

class Quasiharmonic {

 private:

  Quasiharmonic();
  ~Quasiharmonic();

  //references for reference cell
  double volume_ref;
  double a_ref;
  double b_ref;
  double c_ref;
  double alpha_ref;
  double beta_ref;
  double gamma_ref;
  int total_k_point;
  Vector freq_ref;
  Vector Coords_ref;
  Matrix Modes_ref;
  bool ref_init;
  Vector freq_ref_k;//reference frequency for a particular k point
  
  //while getting the hessian frequencies, this string is the type of the hessian,
  //(e.g. ref, Plus or Minus freq)
  string hess_type;

  Vector Gruneisen_parameters;
  //Gruneisen_parameter for lattice lengths
  Vector Grun_a;
  Vector Grun_b;
  Vector Grun_c;
  

  //This bool now handled by Are_QHF_Available in the params class
  //bool Gruneisen_init;

 public:
  static Quasiharmonic& quasiharmonic(){
    static Quasiharmonic quasi;
    return quasi;
  }

  void Initialize();

  //Determining frequencies and volumes
  void DetermineReferenceFrequencies(ifstream &infile,int &imag_num);
  void DetermineIsotropicFrequencies(ifstream &infile,double &VolumePlus,double &VolumeMinus,Vector &FreqPlus,Vector &FreqMinus,int &imag_num);
  void DetermineLatticeFrequencies(ifstream &infile,double &LatticePlus,double &LatticeMinus,Vector &FreqPlus,Vector &FreqMinus,int &imag_num,int params);

  Vector GetReferenceFrequencies() {return freq_ref;};
  Vector GetGruneisenParameters() {return Gruneisen_parameters;};
  double GetReferenceVolume() {return volume_ref;};
  void SetReferenceMode() {ref_init = true;};
  void SetReferenceMode(Matrix Modes) { Modes_ref = Modes; ref_init = true;};
  Vector GetReferenceMode(int i) {return Modes_ref.GetColumnVector(i);};
  bool IsReferenceCellInitialized(){return ref_init;};

  //Setting and Reading Gruneisen parameters
  void ReadIsotropicFrequencies(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus);
 void ReadAnisotropicFrequencies(Vector& plus_a_freq,double& plus_a,Vector& minus_a_freq,double& minus_a,Vector& plus_b_freq,double& plus_b,
				Vector& minus_b_freq,double& minus_b,Vector& plus_c_freq,double& plus_c,Vector& minus_c_freq,double& minus_c);
 void SetGruneisenParameters(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus,int Param = 0);
 //bool IsGruneisenInitialized() {return Gruneisen_init;};
 //void InitializingGruneisen(bool init){Gruneisen_init = init;};
 
 //Calculations
 Vector ComputeAnisotropicFrequencies();
 Vector ComputeAnisotropicGradient();

 void ReadMoldenFrequencies();
 void ReadMoldenFrequencies(int k_point);
 Vector OrderFrequencies(Vector OriginalFreq,Matrix Modes);
 void MangageModeVectors(Matrix& Modes,Vector& OriginalFreq,Vector& Freq,vector<int> entered_entries);

 string GetHessianType(){return hess_type;};
 Vector ReadGeometryFromInput(double& volume);

};

#endif
