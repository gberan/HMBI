#ifndef _dimer_h
#define _dimer_h

#include <limits>
#include <iostream>
#include <cmath> 
using std::istringstream; // part of iostream

#include "monomer.h"
#include "params.h"
#include <vector> // C STL library vector, "vector"

/*

Dimer class. 

Stores 2 monomers, total & interaction energies, forces, etc.

When doing Periodic boundary conditions, MonA is in primary unit cell
and MonB is a periodic image.

Unit conventions:
  Coordinates in Angstroms
  Energies in hartrees
  Gradients in hartrees/bohr  
  Hessians in hartrees/bohr/bohr

GJB 12/08
*/

class Dimer {

  int indexA, indexB; // labels for each monomer
  string typeA, typeB; //monomer type
  Monomer MonA, MonB; // individual monomer objects

  int spin, charge; // spin multiplicity, net charge
  int Natoms; // # of atoms in the dimer

  // for PBC
  bool is_image; // Mon B is a periodic image?  
  int K_vector[3]; // describes which cell MonB comes from.
  int reference_MonB; // monomer with same geom in primary unit cell

  // Symmetry
  int sym_fac; // number of times this dimer occurs due to symmetry
  int sym_fac_period; //number of times this dimer appears in as a periodic
                        //image. For the periodic image dimers, this value is 
                        //is stored in sym_fac
  vector<int> List_Sym_Dim;//list of symmetrical dimers
  vector<int> List_Sym_Periodic;//list of symmetrical periodic dimers
  vector<Matrix> Rotation; //rotation matrix of all symmetrical equivalent dimer
  vector<Matrix> Rotation_Periodic;
  vector<Vector> Atom_Equivalency;//equivalency of atoms on List_Sym_Dim
  vector<Vector> Atom_Equivalency_Periodic;
  vector<int> k_list;//Compliments List_Sym_Dim (for image-dimers) and List_Sym_Periodic (for non-image-dimers)
                     //List of which cell MonB on the list is in


  //Complimentry list for List_Sym_Dim and List_Sym_Periodic 
  vector<bool> MonB_List;//Compliments List_Sym_Dim (for image-dimers) and List_Sym_Periodic (for non-image-dimers)
                        //List whether the even or odd entree monomer on List_Sym_Dimer or List_Sym_Periodic is 
                        //outside the unit cell.
                        //If an entry is 0, an odd entry is outside. If 1 an even entry is.
                        //Needed for finding gradient for lattice parameters when symmetry is imposed.


  // Energies, gradients, etc.
  double Energy_QM, Energy_MM; // stored in hartrees
  double dEint_QM, dEint_MM; // stored in hartrees
  double E_Electrostatic_MM, E_Induction_MM, E_2b_disp_MM; // stored in hartrees
  Vector Grad_QM, Grad_MM; // Gradients in hartrees/bohr  
  Matrix Hess_QM, Hess_MM; // Hessians in hartrees/bohr/bohr
  Vector Grad_Electrostatic; // 2-body electrostatic vector

  double Separation; // separation between the 2 monomers
  int minA, minB; // indecies for the two closest atoms

  int QM_Grad_Init, MM_Grad_Init,Grad_Electrostatic_Init; // whether or not gradients have been read
  int QM_Hess_Init, MM_Hess_Init; //whether or not the Hessian has been read

  double OrientEnergy; // by Ali

  int Tabs_Init, DampedTabs_Init;// whether electrostatic interactions are stored
  Matrix *Tabs; // List of matrices, storing one electrostatic interaction
  Matrix *DampedTabs; // version of Tabs with damping applied to prevent
  // the polarization catastrophe at short ranges.  These are lists of
  // Tab matrices, one matrix for each pair of atoms between the two
  // monomers.

  int TabsGrad_Init, DampedTabsGrad_Init;//whether TabsGrad and DampedTabsGrad are initialized
  Matrix *TabsGrad;
  Matrix *DampedTabsGrad;

  // NMR stuff... JDH
  bool UseInTwoBody_; 


  // Watit - IR & Raman Intensities
	Matrix HessQM;
	Matrix HessMM;
	Matrix DipDA;
	Matrix DipDB;
	Matrix PolDA;
	Matrix PolDB;
  //

 public:
  Dimer(); // The real work is done by Initialize

  Dimer(const Dimer &other);
  
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
    if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
      Grad_QM.Set();
      Grad_MM.Set();
      QM_Grad_Init = 0;
      MM_Grad_Init = 0;
    }
  };
  void ResetHessians() {
    Hess_QM.Set();
    Hess_MM.Set();
    
    QM_Hess_Init = 0;
    MM_Hess_Init = 0;
  };

  // Basic functions
  int GetNumberOfAtoms() {return Natoms;};
  Monomer GetMonomerA() {return MonA;};
  Monomer GetMonomerB() {return MonB;};
  int GetIndexA() {return indexA;};
  int GetIndexB() {return indexB;};
  string GetTypeA() {return typeA;};
  string GetTypeB() {return typeB;};
  Vector GetCurrentCoordinates();

  void SetIsImage(bool image) {is_image = image;};
  bool IsPeriodicImage() {return is_image;};


  // Symmetry
  bool SymmetryCheck(Dimer Dimers[],int index,bool ImageToNonImage = false);
  bool DetermineSymmetry(Dimer Dimers[],Matrix LocalCoor, Matrix ThisInerita,string AtomList[], int index,bool ImageToNonImage = false);
  int GetSymmetryFactor() {return sym_fac;};
  int GetPeriodicSymmetryFactor() {return sym_fac_period;};
  void SetSymmetryFactor(int fac) {sym_fac = fac;};
  void SetPeriodicSymmetryFactor(int fac) {sym_fac_period = fac;};
  void AddToSymmetryList(Dimer& other,Matrix Rot,bool ImageToNonImage = false);
  void MergeEquivalencyList(vector<Vector> OtherList,Vector ThisList,bool ImageToNonImage = false);
//  void MergeSymmetryList(vector<int> OtherSymList,vector<bool> OtherMonBList,vector<int> OtherKList,bool ImageToNonImage,bool InverseList = false);
  void MergeSymmetryList(vector<int> OtherSymList,vector<int> OtherKList,bool ImageToNonImage);
  void AddToRotationList(Matrix Rot,vector<Matrix> OtherMatrixList,bool ImageToNonImage,bool InverseList = false);
  void ResetSymmetryList();
  void AlterSymmetryList(Dimer& SymDimer,int ListNumb,Matrix Rot,bool ImageToNonImage);
  vector<int> GetSymmetryList() {return List_Sym_Dim;};
  vector<int> GetPeriodicSymmetryList() {return List_Sym_Periodic;};
  vector<bool> GetMonBList() {return MonB_List;};
  vector<int> GetSymmetricalImageCell() {return k_list;}
  void AlterRotationList(Matrix Rot,int i){Rotation[i] = Rot;};
  void AlterPeriodicRotationList(Matrix Rot,int i){Rotation_Periodic[i] = Rot;};
  bool IsDimerPlanarOrLinear();

  vector<Matrix> GetRotationList() {return Rotation;};
  vector<Matrix> GetPeriodicRotationList() {return Rotation_Periodic;};
  vector<Vector> GetAtomEquivalency() {return Atom_Equivalency;};
  vector<Vector> GetPeriodicAtomEquivalency() {return Atom_Equivalency_Periodic;};

  // Periodic Boundary Conditions
  void SetImageCell(int x, int y, int z) {K_vector[0]=x; K_vector[1]=y; 
    K_vector[2]=z; is_image = true;};
  int* GetImageCell() {return K_vector;};
  
  void SetReferenceMonomerIndex(int index);// {reference_MonB = index; };
  int GetReferenceMonomerIndex() {return reference_MonB;};

  // Create and Run jobs
  void CreateQMJob(Monomer Monomers[], int NMon);
  void CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job = false);
  void CreateMolProJob(Monomer Monomers[], int NMon);
  void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images); //JDH NMR code
  void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon); //JDH NMR code
  void CreateMolProCCSDTJob(Monomer Monomers[0], int NMon);
  void CreatePSI4Job(Monomer Monomers[], int NMon);
  void CreateFiniteDifferenceMolProJob(Monomer Monomers[], int NMon,Vector Coord, bool CCSDT);
  void CreateFiniteDifferenceMolProJob(Monomer Monomers[], int NMon,bool CCSDT = true);
  void CreateEnergyFiniteDifferenceHessianJobs(Monomer Monomers[], int NMon,bool CCSDT);
  void CreateMMJob(Monomer Monomers[], int NMon);
  void CreateTinkerJob(Monomer Monomers[], int NMon);
  string RunQMJob(bool CCSDT=false);
  string RunQChemJob(bool MM_job=false);
  string RunMolProJob(bool CCSDT=false);
  string RunPSI4Job(bool MM_job=false);
  string RunFiniteDifferenceMolProJob();
  string RunFiniteDifferenceMolProJob(int i,bool Plus,bool CCSDT);
  string RunFiniteDifferenceMolProJob(int i,bool Plus1,int j,bool Plus2, bool CCSDT);
  string RunMMJob();
  string RunTinkerJob();

  // Read Dimer energies from disk
  void ReadQMResults();
  void ReadMMResults();
  void ReadQMResults(Monomer& otherB); // use alternate MonB as reference
  void ReadMMResults(Monomer& otherB); // use alternate MonB as reference
  double ReadQChemEnergy(bool MM_job=false);
  double ReadMolProEnergy();
  double ReadPSI4Energy(); //CSG
  double ReadG09Energy(); //Watit
  double ReadQChemCounterpoiseInteractionEnergy();
  void SetMMEnergy();
  double ReadTinkerEnergy();
  void SetOrientEnergy(double energy){OrientEnergy=energy;};  // by Ali
  double GetOrientEnergy() {return OrientEnergy;};
  
  double GetQMEnergy() {return Energy_QM;};
  double GetMMEnergy() {return Energy_MM;};
  double GetQMIntEnergy() {return dEint_QM;};
  double GetMMIntEnergy() {return dEint_MM;};


  Vector GetSpatialDampingFunctionGradient(double c0, double c1);
  Matrix GetSpatialDampingFunctionHessian(double c0, double c1);

  // Read Dimer gradients from disk
  void SetQMGradient();
  void SetMMGradient();
  Vector ReadGradient(string path, int type); 
  Vector ReadFiniteDifferenceMolProGradient(string path,bool CCSDT);
  // or set them from a vector
  //  Vector ReadQMCPGradient(string path);
  void SetQMGradient(Vector grad) {Grad_QM = grad; QM_Grad_Init = 1;};
  void SetMMGradient(Vector grad) {Grad_MM = grad; MM_Grad_Init = 1;};
  void SetQMHessian(Matrix hess) {Hess_QM = hess; QM_Hess_Init = 1;};
  void SetMMHessian(Matrix hess) {Hess_MM = hess; QM_Hess_Init = 1;};

  Vector GetQMGradient() {return Grad_QM;};
  Vector GetMMGradient() {return Grad_MM;};
  Vector GetMM2BodyElectrostaticGradient() {return Grad_Electrostatic;};

  Matrix GetQMHessian() {return Hess_QM;};
  double GetQMHess_value(int i,int j){return Hess_QM(i,j);};      //Watit
  Matrix GetMMHessian() {return Hess_MM;};

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
  void PrintDimerCartesian(FILE *outfile=stdout);
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintQChemGhostCartesian(FILE *outfile=stdout,bool printGhostMonA=false,bool printGhostMonB=false);
  void PrintTinkerCartesian(FILE *outfile=stdout);
  void PrintQMGradient(string title);
  void PrintMMGradient(string title);
  void PrintGradient(string title, Vector grad);
  void PrintAll();

  void ComputeIntermolecularDistances();
  Matrix ReadHessian(string path, int type);
  Matrix ReadFiniteDifferenceMolProHessian(string path);
  Matrix ReadEnergyFiniteDifferenceHessian(string path,bool CCSDT);
  void SetQMHessian();
  void SetMMHessian();


//yoni :added from Kaushik's code: used for multipole damping factors.
  double GetAIFFDampingFactor();




  /*** AIFF routines ***/

  double ComputeAIFFEnergy(); // Main driver

  /* Electrostatics & Induction */
  double ComputeMultipoleInteractions();
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

  //overload operators
  Dimer& operator=(const Dimer &other);

  // JDH NMR stuff...
  void CreateQChemJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );

  void CreateOrcaJob(Monomer Monomers[], int NMon );
  void CreateOrcaJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ); 
  void CreateOrcaJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );

  void CreateDaltonJob(Monomer Monomers[], int NMon );
  void CreateDaltonJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );
  void CreateDaltonJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );

  void CreateG09Job(Monomer Monomers[], int NMon, bool MM_job = false );
  void CreateG09Job(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );
  void CreateG09Job(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );

  string RunG09Job(bool MM_job=false);
  string RunDaltonJob();
  string RunOrcaJob();
  
  void ReadMolProNMRdata();
  void ReadDaltonNMRdata();
  void ReadG09NMRdata();
  void ReadG09EFGData();
  void ReadQChemNMRdata();
  void ReadOrcaNMRdata();

  bool GetUseInTwoBodyCalculation() { return UseInTwoBody_;}; //JDH
  void SetUseInTwoBodyCalculation(bool INPUT){ UseInTwoBody_ = INPUT;}; //JDH

  

  void ReadChelpGCharges(); // BETA JDH
  void ReadHirshfeldCharges(); // BETA JDH

  void CreateG09ChelpGJob(); // BETA JDHA
  void CreateG09ChelpGJob(int NMon, Monomer Monomers[]); // BETA JDH
  void CreateG09ChelpGJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ); // BETA JDH

  void CreateG09HirshfeldJob(); // BETA JDH
  void CreateG09HirshfeldJob(int NMon, Monomer Monomers[]); // BETA JDH
  void CreateG09HirshfeldJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ); // BETA JDH

  string RunG09ChelpGJob(); //BETA JDH
  string RunG09HirshfeldJob(); //BETA JDH

  // Watit - IR & Raman Intensities
	void ReadHess();
	void SetHessQM(Matrix in) {HessQM = in;};
	void SetHessMM(Matrix in) {HessMM = in;};
	Matrix GetHessQM_Matrix(){return HessQM;};
	Matrix GetHessMM_Matrix(){return HessMM;};
	double GetHessQM_value(int i,int j){return HessQM(i,j);};
	double GetHessMM_value(int i,int j){return HessMM(i,j);};

	void ReadIntensity();
        void SetDipDA(Matrix inA) {DipDA = inA;};
        void SetDipDB(Matrix inB) {DipDB = inB;};
        Matrix GetDipDA_Matrix(){return DipDA;};
        Matrix GetDipDB_Matrix(){return DipDB;};  
        double GetDipDA_value(int i,int j){return DipDA(i,j);};
        double GetDipDB_value(int i,int j){return DipDB(i,j);};
	void SetPolDA(Matrix inA) {PolDA = inA;};
	void SetPolDB(Matrix inB) {PolDB = inB;};
	Matrix GetPolDA_Matrix(){return PolDA;};
	Matrix GetPolDB_Matrix(){return PolDB;};  
	double GetPolDA_value(int i,int j){return PolDA(i,j);};
	double GetPolDB_value(int i,int j){return PolDB(i,j);};
  //
};

#endif
