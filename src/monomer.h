#ifndef _monomer_h
#define _monomer_h

#include <iostream>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream.h>
#include <math.h>// by shuhao
#include <algorithm>// by shuhao
using namespace std;// by shuhao
#include <cassert>
#include <iomanip.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h> //for getcwd
using std::string;
#include "atom.h"
#include "params.h"
#include "vector.h"

/*
  A Monomer class, storing a collection of atoms, the total energy
(and possibly gradients), and other relevant parameters.

Unit conventions:
  Coordinates in Angstroms
  Energies in hartrees
  Gradients in hartrees/bohr  

GJB 01/09

 */


class Monomer {

  int Monomer_index; // assigns a unique number to each monomer
  int ref_mon; // for PBC image monomers: associate this monomer with
	       // its central unit cell version
  string mN_; // monomer label, e.g. m3 if Monomer_index = 3;
  int spin, charge;
  int Natoms;
  Atom *Atom_List; // Coordinates in angstroms

  // If doing AIFF, need IP for long-range DFT correction
  double IonizationPotential; 

  double MonomerMass;
  Vector CenterOfMass;

  // Flags to determine whether or not the jobs have been run
  bool QM_Job_Complete; 
  bool MM_Job_Complete;

  double Energy_QM, Energy_MM; // stored in hartrees
  Vector Grad_QM, Grad_MM;  // stored in hartrees/bohr  

  int QM_Grad_Init, MM_Grad_Init; // whether or not gradients have been read

  // Define relationship between monomer local coordinates and global
  // coordinates
  double RotationAngle; // in radians by Ali
  Vector RotationVector; // by Ali

 public:
  Monomer();  // The real work is done by Initialize
  Monomer(const Monomer& other);  // The real work is done by Initialize
  ~Monomer();
  // The function to populate the objects:

  // Q-Chem style geometry init
  void Initialize(int ind, int charge_in, int spin_in, string atoms[], 
		  double *xyz, int natoms);

  // Q-Chem style geometry init, with point charges
  void Initialize(int ind, int charge_in, int spin_in, string atoms[], 
		  double *xyz, int natoms, Vector charges);
  
  // Tinker-style geometry init
  void Initialize(int ind, int charge_in, int spin_in, string atoms[], 
		  double *xyz, int natoms, int *atom_types,
		  int *Nconnected, int* connectivity);

  // Tinker-style init, with point charges
  void Initialize(int ind, int charge_in, int spin_in, string atoms[], 
		  double *xyz, int natoms, int *atom_types,
		  int *Nconnected, int* connectivity, Vector charges);

  // Reset all the data
  void ResetEnergiesAndGradients() {
    Energy_QM = 0.0;
    Energy_MM = 0.0;
    Grad_QM.Set();
    Grad_MM.Set();
  };

  // Get some basic properties
  void FindCenterOfMass();
  double GetCenterOfMass(int dim) {return CenterOfMass[dim];};
  Vector GetCenterOfMass() {return CenterOfMass;};
  void SetIndex(int i) {Monomer_index = i;};
  int GetIndex() {return Monomer_index;};
  void SetReferenceMonomerIndex(int i) {ref_mon = i;};
  int GetReferenceMonomerIndex() {return ref_mon;};
  void SetLabel();
  string GetLabel() {return mN_;};
  int GetSpinState() {return spin;};
  int GetChargeState() {return charge;};
  int GetNumberOfAtoms() {return Natoms;};
  double GetIonizationPotential() {return IonizationPotential;};
  void SetIonizationPotential(double IP) {IonizationPotential = IP;};
  double GetMonomerMass() {return MonomerMass;};
  double GetRotationAngle() {return RotationAngle;};  // by Ali
  double GetRotationVector(int dim) {return RotationVector[dim];};  // by Ali
  Vector GetRotationVector() {return RotationVector;};  

  Atom& GetAtom(int i) {return Atom_List[i];};
  string GetSymbol(int i) {return Atom_List[i].GetSymbol();};
  int GetAtomicNumber(int i) {return Atom_List[i].GetAtomicNumber();};

  // Create and Run jobs
  void CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job=false);
  void CreateMMJob(Monomer Monomers[], int NMon);
  void CreateTinkerJob(Monomer Monomers[], int NMon);
  void CreateCamCaspJob();  // by Ali
  string RunQChemJob(bool MM_job=false);
  string RunMMJob();
  string RunTinkerJob();
  string RunCamCaspJob();  

  // Read Monomer energies from disk
  void ReadQMResults();
  void ReadMMResults();
  double ReadQChemEnergy(bool MM_job=false);
  void SetQMEnergy(double E) {Energy_QM = E;};
  void SetMMEnergy();
  void SetMMEnergy(double E) {Energy_MM = E;};
  double ReadTinkerEnergy();
  void ReadMultipoleMoments();
  void ReadPolarizabilities();
  void ReadIsotropicDipolePolarizability();
  void ReadDispersionCoefficients();
  void ReadFreqPolarizabilities(); // by shuhao
  void ReadDiagonalAnisotropicFreqPolarizabilities();
  void SetEmpiricalAtomDispersionType();


  void SetRotationAngle(double RotAng){RotationAngle = RotAng;};
  void SetRotationVector(int dim, double RotVec){RotationVector[dim] = RotVec;};  // by Ali
  
  double GetQMEnergy() {return Energy_QM;};
  double GetMMEnergy() {return Energy_MM;};

  // Read Monomer gradients from disk
  void SetQMGradient();
  void SetMMGradient();
  // or from another monomer
  void SetQMGradient(Monomer& other);
  void SetMMGradient(Monomer& other);
  Vector ReadGradient(string path, int type); 
  // or just from a vector
  void SetQMGradient(Vector grad) {Grad_QM = grad;};
  void SetMMGradient(Vector grad) {Grad_MM = grad;};

  Vector GetQMGradient() {return Grad_QM;};
  Vector GetMMGradient() {return Grad_MM;};


  Vector FindDistance(const Monomer& Mon);

  // Printing
  void PrintMonomerCartesian(FILE *outfile=stdout);
  void PrintGhostMonomerCartesian(FILE *outfile=stdout);
  // this one lets you substitute an alternative atomic symbol for printing,
  // i.e. to change color of key fragments when visualizing
  void PrintMonomerCartesian(FILE *outfile, string symbol); 
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintTinkerCartesian(int shift=0, bool monomer=true, 
			    FILE *outfile=stdout);
  void PrintQChemEmbeddingCharges(FILE* outfile=stdout);
  void PrintTinkerEmbeddingCharges(int shift=0, FILE* outfile=stdout);

  void PrintQMGradient(string title);
  void PrintMMGradient(string title);
  void PrintGradient(string title, Vector grad);
  void PrintAll();

  // Modify monomer geometry
  void Translate(Vector new_com);
  void SetAtomPosition(int i, double xyz[3]) {Atom_List[i].SetPosition(xyz);};

  // 
  void ComputeIntramolecularDistances();

  // Overload operators
  Monomer& operator=(const Monomer& other);

  void FindLocalCoord(bool use_old_axes=false); // Ali

  
  /*** AIFF functions ***/
  
  // Really, the AIFF has no intramolecular terms.  However, we do
  // have some routines that look at interatomic 3-body dispersion or
  // to correct MP2 2-body dispersion via Delta(vdw).

  /* 3-Body Dispersion */
  // We currently do not include this sort of contribution, since it
  // is hard to get the short-range damping right.
  double ComputeThreeBodyDispersion(); 


  // Old empirical Tkatchenko-Scheffler 3-body dispersion functions
  double ComputeInteratomicThreeBodyDispersion(); 

  // Delta(vdw) style dispersion correction between pairs of atoms on
  // a single monomer.
  double ComputeInteratomicMP2TwoBodyDispersionCorrection(); 




  string StringToUpper(string str) {
  for (int j=0; j<str.length(); ++j) 
    str[j]=toupper(str[j]);
  return str;
  }


};

#endif
