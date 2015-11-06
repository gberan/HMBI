#ifndef _dimer_h
#define _dimer_h

#include <iostream>
using std::istringstream; // part of iostream

#include "monomer.h"
#include "params.h"

/*

Dimer class. 

Stores 2 monomers, total & interaction energies, forces, etc.

When doing Periodic boundary conditions, MonA is in primary unit cell
and MonB is a periodic image.

Unit conventions:
  Coordinates in Angstroms
  Energies in hartrees
  Gradients in hartrees/bohr  


GJB 12/08
*/

class Dimer {

  int indexA, indexB; // labels for each monomer
  Monomer MonA, MonB; // individual monomer objects

  int spin, charge; // spin multiplicity, net charge
  int Natoms; // # of atoms in the dimer

  // for PBC
  bool is_image; // Mon B is a periodic image?  
  int K_vector[3]; // describes which cell MonB comes from.
  int reference_MonB; // monomer with same geom in primary unit cell

  // Symmetry
  int sym_fac; // number of times this dimer occurs due to symmetry

  // Energies, gradients, etc.
  double Energy_QM, Energy_MM; // stored in hartrees
  double dEint_QM, dEint_MM; // stored in hartrees
  double E_Electrostatic_MM, E_Induction_MM, E_2b_disp_MM; // stored in hartrees
  Vector Grad_QM, Grad_MM; // Gradients in hartrees/bohr  
  Vector Grad_Electrostatic; // 2-body electrostatic vector

  double Separation; // separation between the 2 monomers
  int minA, minB; // indecies for the two closest atoms

  int QM_Grad_Init, MM_Grad_Init; // whether or not gradients have been read

  double OrientEnergy; // by Ali

  Matrix *Tabs; // List of matrices, storing one electrostatic interaction
  Matrix *DampedTabs; // version of Tabs with damping applied to prevent
  // the polarization catastrophe at short ranges.  These are lists of
  // Tab matrices, one matrix for each pair of atoms between the two
  // monomers.

  Matrix *TabsGrad;
  Matrix *DampedTabsGrad;


  // Private copy constructors/assignment function declarations to
  // prevent inadvertant implicit copying.  If you need copying, make
  // these public and define them.
  Dimer(const Dimer &other);
  Dimer& operator=(const Dimer &other);
  


 public:
  Dimer(); // The real work is done by Initialize

  
  ~Dimer();
  // The function to populate the objects:
  void Initialize(Monomer& m1, Monomer& m2);
  // The function updates the object, in case monomers/params have changed
  void UpdateObjects(Monomer& m1, Monomer& m2);
  // Updates for PBC, uses reference monomer for monomer energies/gradients
  void UpdateObjects(Monomer& m1, Monomer& m2, Monomer& ref_m2);

  // Reset all the parameters
  void ResetEnergiesAndGradients() {
    Energy_QM = 0.0; Energy_MM = 0.0;
    dEint_QM = 0.0; dEint_MM = 0.0;
    E_Electrostatic_MM = 0.0; 
    E_Induction_MM = 0.0;
    E_2b_disp_MM = 0.0;
    Grad_QM.Set();
    Grad_MM.Set();
  };

  // Basic functions
  int GetNumberOfAtoms() {return Natoms;};
  Monomer GetMonomerA() {return MonA;};
  Monomer GetMonomerB() {return MonB;};
  int GetIndexA() {return indexA;};
  int GetIndexB() {return indexB;};
  void SetIsImage(bool image) {is_image = image;};
  bool IsPeriodicImage() {return is_image;};

  // Symmetry
  int GetSymmetryFactor() {return sym_fac;};
  void SetSymmetryFactor(int fac) {sym_fac = fac;};

  // Periodic Boundary Conditions
  void SetImageCell(int x, int y, int z) {K_vector[0]=x; K_vector[1]=y; 
    K_vector[2]=z; is_image = true;};
  int* GetImageCell() {return K_vector;};
  
  void SetReferenceMonomerIndex(int index) {reference_MonB = index; };
  int GetReferenceMonomerIndex() {return reference_MonB;};

  // Create and Run jobs
  void CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job = false);
  void CreateMMJob(Monomer Monomers[], int NMon);
  void CreateTinkerJob(Monomer Monomers[], int NMon);
  string RunQChemJob(bool MM_job=false);
  string RunMMJob();
  string RunTinkerJob();

  // Read Dimer energies from disk
  void ReadQMResults();
  void ReadMMResults();
  void ReadQMResults(Monomer& otherB); // use alternate MonB as reference
  void ReadMMResults(Monomer& otherB); // use alternate MonB as reference
  double ReadQChemEnergy(bool MM_job=false);
  double ReadQChemCounterpoiseInteractionEnergy();
  void SetMMEnergy();
  double ReadTinkerEnergy();
  void SetOrientEnergy(double energy){OrientEnergy=energy;};  // by Ali
  double GetOrientEnergy() {return OrientEnergy;};
  
  double GetQMEnergy() {return Energy_QM;};
  double GetMMEnergy() {return Energy_MM;};
  double GetQMIntEnergy() {return dEint_QM;};
  double GetMMIntEnergy() {return dEint_MM;};


  // Read Dimer gradients from disk
  void SetQMGradient();
  void SetMMGradient();
  Vector ReadGradient(string path, int type); 
  // or set them from a vector
  //  Vector ReadQMCPGradient(string path);
  void SetQMGradient(Vector grad) {Grad_QM = grad;};
  void SetMMGradient(Vector grad) {Grad_MM = grad;};

  Vector GetQMGradient() {return Grad_QM;};
  Vector GetMMGradient() {return Grad_MM;};
  Vector GetMM2BodyElectrostaticGradient() {return Grad_Electrostatic;};

  // For Local truncation
  double GetDampingFactor(double c0, double c1);
  int GetNzero(double c0, double c1);
  int GetNfull(double c0, double c1);
  int GetNdamp(double c0, double c1);


  double GetDimerSeparation() {return Separation;};
  int GetMinimumDistanceAtomA() {return minA;};
  int GetMinimumDistanceAtomB() {return minB;};


  // Printing
  void PrintDimerCartesian();
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintTinkerCartesian(FILE *outfile=stdout);
  void PrintQMGradient(string title);
  void PrintMMGradient(string title);
  void PrintGradient(string title, Vector grad);
  void PrintAll();

  void ComputeIntermolecularDistances();



  /*** AIFF routines ***/

  double ComputeAIFFEnergy(); // Main driver

  /* Electrostatics & Induction */
  double ComputeMultipoleInteractions(); // uses old, iterative solver for polarization.  Slow.
  double ComputeAIFFElectrostatics(); // Just electrostatics
  double ComputeAIFFInduction(); // Faster matrix-based approach for induction.
  Vector ComputeMultipoleGradient(); 

  void BuildDampedTabInteractionMatrices();
  void BuildTabInteractionMatrices();

  Matrix GetTabInteractionMatrix(int iatomA, int iatomB) 
  {return Tabs[iatomB*MonA.GetNumberOfAtoms()+iatomA];};
  
  Matrix GetDampedTabInteractionMatrix(int iatomA, int iatomB) 
  {return DampedTabs[iatomB*MonA.GetNumberOfAtoms()+iatomA];};

  Matrix GetTabGradient(int iatomA, int iatomB) 
  {return TabsGrad[iatomB*MonA.GetNumberOfAtoms()+iatomA];};
  
  Matrix GetDampedTabGradient(int iatomA, int iatomB) 
  {return DampedTabsGrad[iatomB*MonA.GetNumberOfAtoms()+iatomA];};
  
  // Grab components for AIFF:
  double GetMMElectrostaticEnergy() {return E_Electrostatic_MM;};
  double GetMMInductionEnergy() {return E_Induction_MM;};
  double GetMMDispersionEnergy() {return E_2b_disp_MM;};

  /* 2-Body Dispersion */
  double ComputeTwoBodyDispersion(); // using AIFF C6 & C8 coefficientes
  
  // The next 3 routines calculate just a piece of the 2-body
  // dispersion.  They are mostly for analysis purposes.  Use
  // ComputeTwoBodyDispersion() for practical applications.
  double ComputeTwoBodyDispersionC6(); // C6 piece only
  double ComputeTwoBodyDispersionC8(); // C8 piece only
  double ComputeTwoBodyDispersionC10(); // C10 piece only


  /* 3-Body Dispersion */
  // This is only used if you want to look at the interatomic
  // contribution involving 2 monomers.
  double ComputeThreeBodyDispersion();


  // Old empirical Tkatchenko-Scheffler dispersion functions
  double EstimateTwoBodyDispersion(); 
  double EstimateThreeBodyDispersion(); // (old)

  // For Delta(vdw)-type MP2 corrections
  double ComputeInteratomicMP2TwoBodyDispersionCorrection();

};

#endif
