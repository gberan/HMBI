#ifndef _monomer_h
#define _monomer_h

#include <cstring>
#include <iostream>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream>
#include <math.h>// by shuhao
#include <algorithm>// by shuhao
using namespace std;// by shuhao
#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h> //for getcwd
using std::string;
#include "atom.h"
#include "params.h"
#include "vector.h"
#include "matrix.h"

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
  int ref_mon; // for PBC and supercell image monomers: associate this monomer with
	       // its central unit cell version
  string mN_; // monomer label, e.g. m3 if Monomer_index = 3;
  int spin, charge;
  int Natoms;
  int unique_atoms;// # of symmetrical unique atoms in the molecule
  Atom *Atom_List; // Coordinates in angstroms
  string Monomer_type; 

  // Symmetry
  int sym_fac; // number of times this monomers occurs due to symmetry  
  int SymMon;//the index number of the monomer this monomer is symmetrical to
  int *Sym_Atom;//indices of the atoms in the symmetrical unique monomer that the atom in this monomer is equivalent. Assumes a one-to-one atom equivalency
  Matrix Symmetry_Rotation;//for gradiates and AIFF, tell us how to rotate to molecule to match the unique monomer.
  Matrix Fractional_Rotation;//the symmetry rotation to the unique monomer using fractional coordinates
  Vector Fractional_Translate;//the translation after the rotation matrix is applied to get the coordinates of the unique monomer
  //Matrix Local_Rotation;//rotation to the local coordinates during optimization. Necessary for some point groups
  //Vector Opt_Center;//The center when atoms are opimized around if they are in local coordinates


  //string SymMon;//definds which monomer the nonunique monomer is symmetrical too. 
  //vector<int> SymmetricalMonomers;//list the monomers that are symmetrical for the unique monomer only

  
                                        
  // If doing AIFF, need IP for long-range DFT correction
  double IonizationPotential; 

  double MonomerMass;
  Vector CenterOfMass;

  // Flags to determine whether or not the jobs have been run
  bool QM_Job_Complete; 
  bool MM_Job_Complete;

  double Energy_QM, Energy_MM; // stored in hartrees
  Vector Grad_QM, Grad_MM;  // stored in hartrees/bohr  
  Matrix Hess_QM, Hess_MM;  // stored in hartress/bohr/bohr

  int QM_Grad_Init, MM_Grad_Init, QM_Hess_Init, MM_Hess_Init; // whether or not gradients have been read


  //No longer using, Rotation Matrix used instead.
  // Define relationship between monomer local coordinates and global
  // coordinates
  double RotationAngle; // in radians by Ali
  Vector RotationVector; // by Ali

  // NMR stuff...
  bool UseInEwald_;
  //bool UseInEwaldBuffer_;
  bool UseInEmbedding_;
  bool UseInTwoBodyCalculation_; 
  bool UseInClusterCalculation_;

  // Watit - IR & Raman Intensities
	Matrix HessQM;
	Matrix HessMM;
	Matrix DipD;
	Matrix PolD;
  //

 public:
  Monomer();  // The real work is done by Initialize
  Monomer(const Monomer& other);  // The real work is done by Initialize
  
  ~Monomer();
  // The function to populate the objects:

 // Q-Chem style geometry init
  void Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, string atoms[], 
		  double *xyz, int natoms);

  // Q-Chem style geometry init, with point charges
  void Initialize(int ind, int natoms_previous,int charge_in, int spin_in, string type, string atoms[], 
		  double *xyz, int natoms, Vector charges);
  
  // Tinker-style geometry init
  void Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, string atoms[], 
		  double *xyz, int natoms, int *atom_types,
		  int *Nconnected, int* connectivity);


  // Tinker-style init, with point charges
  void Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, string atoms[], 
		  double *xyz, int natoms, int *atom_types,
		  int *Nconnected, int* connectivity, Vector charges);


  // Reset all the data
  void ResetEnergiesAndGradients() {
    Energy_QM = 0.0;
    Energy_MM = 0.0;
    if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
      Grad_QM.Set();
      Grad_MM.Set();
    }
  };

  void ResetHessians() {
    Hess_QM.Set();
    Hess_MM.Set();
  };

  // Get some basic properties
  void FindCenterOfMass();
  double GetCenterOfMass(int dim) {return CenterOfMass[dim];};
  Vector GetCenterOfMass() {return CenterOfMass;};
  void SetIndex(int i) {Monomer_index = i;};
  int GetIndex() {return Monomer_index;};
  void SetReferenceMonomerIndex(int i) {ref_mon = i;};
  int GetReferenceMonomerIndex() {return ref_mon;};
  string GetType() {return Monomer_type;};  
  void SetLabel();
  string GetLabel() {return mN_;};
  int GetSpinState() {return spin;};
  int GetChargeState() {return charge;};
  int GetNumberOfAtoms() {return Natoms;};
  int GetNumberOfUniqueAtoms() {return unique_atoms;};
  double GetIonizationPotential() {return IonizationPotential;};
  void SetIonizationPotential(double IP) {IonizationPotential = IP;};
  double GetMonomerMass() {return MonomerMass;};
  double GetRotationAngle() {return RotationAngle;};  // by Ali
  double GetRotationVector(int dim) {return RotationVector[dim];};  // by Ali
  Vector GetRotationVector() {return RotationVector;};  
  Atom& GetAtom(int i) {return Atom_List[i];};
  string GetSymbol(int i) {return Atom_List[i].GetSymbol();};
  int GetAtomicNumber(int i) {return Atom_List[i].GetAtomicNumber();};
  double GetAtomicMass(int i) {return Atom_List[i].GetAtomicMass();};
  Vector GetCoordinates();

  // Symmetry
  bool SymmetryCheck(Monomer Monomers[]);
  bool SymmetryForAsymmetricTop(Monomer Monomers[],Matrix LocalCoor,Matrix Inertia);
  bool SymmetryForLinearTop(Monomer Monomers[],Matrix LocalCoor,Matrix Rot);
  bool SymmetryForSymmetricTop(Monomer Monomers[],Matrix LocalCoor, Matrix Inertia);
  bool SymmetryForSphericalTop(Monomer Monomers[],Matrix LocalCoor, Matrix Rot);
  int GetSymmetryFactor() {return sym_fac;};
  void SetSymmetryFactor(int sym) {sym_fac = sym;};
  void IncreaseSymmetryFactor() {sym_fac++;};
  Matrix OrientateSphericalTop(Matrix LocalCoor, Matrix &Rotation);
  Matrix OrientateLinearTop(Matrix LocalCoor, Matrix &Rotation);
  Matrix OrientateSymmetricalTop(Matrix LocalCoor,Matrix &PrevRot);
  void FindSymmetricalEquivalentAtoms(Monomer& UniqueMon);
  void DeterminePointGroup(Matrix LocalCoor,Matrix RotToLocal,int molecule_type);
  bool MatchCoordinates(Matrix Coor1,Matrix Coor2);
  //void DefindSymmetricalMonomer(int i) {SymmetricalMonomers = i;};
  int GetSymmetricalMonomer() {return SymMon;};
  int GetSymmetricalAtom(int i) {return Sym_Atom[i];};
  //void AddToSymmetryList(int i) {SymmetricalMonomers.push_back(i);};
  void SetRotationMatrix(Matrix Rotation) {Symmetry_Rotation = Rotation;};
  Matrix GetRotationMatrix() {return Symmetry_Rotation;};
  //Matrix GetLocalMatrix() {return Local_Rotation;};
  //Vector GetLocalCenter() {return Opt_Center;};
  void SetFractRotationMatrix(Matrix FractRot) {Fractional_Rotation = FractRot;};
  Matrix GetFractionalRotationMatrix() {return Fractional_Rotation;};
  void SetFractTranslationVector(Monomer& UniqueMon);
  Vector GetFractTranslationVector(bool FromUniqueToThis); //{return Fractional_Translate;};
  void ResetSymmetry();

 //this functions are out of date. Do not use.
  void FindMapping(Monomer& Monomer);
  bool MappingMatrixForNonPlanar(Monomer& UniqueMon, Matrix COM, Matrix AtomMatrix, string Symbol[]);
  bool MappingMatrixForPlanar (Monomer& UniqueMon, Matrix COM, Matrix AtomMatrix, string Symbol[]);
  bool MappingMatrixForLinear (Monomer& UniqueMon, Matrix COM);


  //////////////////////////////////////////////////////
  // JDH's functions for Magnetic property calculations:
  //////////////////////////////////////////////////////

  // CREATE JOBS
  void CreateGDMAJob();
  void ModifyGDMAMomFile();
  
  // CHELPG
  void CreateG09ChelpGJob();
  void CreateG09ChelpGJob(Monomer Monomers[], int NMon);
  void CreateG09ChelpGJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ); 
  void CreateG09ChelpGJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );

  void CreateOrcaChelpGJob();


  // HIRSHFELD
  void CreateG09HirshfeldJob();
  void CreateG09HirshfeldJob(Monomer Monomers[], int NMon);
  void CreateG09HirshfeldJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );


  void CreateG09Job( Monomer Monomers[], int NMon, bool MM_job=false );
  void CreateG09Job( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );
  void CreateG09Job( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );

  // Dalton and QChem not yet compitible with advanced charge embedding schemes
  void CreateDaltonJob( Monomer Monomers[], int NMon);
  void CreateDaltonJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );
  void CreateDaltonJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );  

  
  void CreateOrcaJob( Monomer Monomers[], int NMon );
  void CreateOrcaJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images );
  void CreateOrcaJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges );


  void CreateQChemJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ); 

  void CreateMolProJob(Monomer Monomers[],int NMon);
  void CreateMolProCCSDTJob(Monomer Monomers[],int NMon);
  void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon);
  void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images);

  // RUN JOBS
  string RunGDMAJob();
  string RunG09GDMAJob();
  string RunFormChk();
  
  string RunG09ChelpGJob();
  string RunG09HirshfeldJob();
  string RunOrcaChelpGJob();
 
  string RunG09Job(bool MM_job=false); //JDH + Watit
  string RunOrcaJob();
  string RunDaltonJob();

  // LABEL MONOMERS
  bool GetUseInTwoBodyCalculation() { return UseInTwoBodyCalculation_;}; //JDH
  void SetUseInTwoBodyCalculation(bool INPUT){ UseInTwoBodyCalculation_ = INPUT;}; //JDH
  bool GetUseInClusterCalculation() { return UseInClusterCalculation_;}; //JDH
  void SetUseInClusterCalculation(bool INPUT){ UseInClusterCalculation_ = INPUT;};  


  
  // ELECTROSTATIC EMBEDDING
  bool GetUseInEmbedding() { return UseInEmbedding_;};
  void SetUseInEmbedding(bool INPUT){ UseInEmbedding_ = INPUT;};
  void PrintEmbeddingCharges(FILE* outfile); // JDH NMR code


  // EWALD
  bool GetUseInEwald() { return UseInEwald_;};
  void SetUseInEwald(bool INPUT){ UseInEwald_ = INPUT;};
  //bool GetUseInEwaldBuffer() { return UseInEwaldBuffer_;};
  //void SetUseInEwaldBuffer(bool INPUT){ UseInEwaldBuffer_ = INPUT;};
  

  // READ DATA
  void ReadDaltonNMRdata();
  void ReadDaltonClusterNMRData();
  void ReadOrcaNMRdata();

  double ReadG09Energy();
  void ReadG09NMRdata();
  void ReadG09ClusterNMRData();
  void ReadG09EFGData();
  void ReadG09ClusterEFGData();
  void ReadG09ChelpGCharges();
  void ReadOrcaChelpGCharges();
  void ReadHirshfeldCharges();

  void ReadMolProNMRdata();
  
  void ReadQChemNMRdata();


  int GetUnitCellIndex( int i);
  void SetUnitCellIndex(int i, int index); 
  int na_; // Inidialize these all to zero here? nope... ISO C++ "forbids" it... *sigh*
  int nb_;
  int nc_;
  //////////////////////////////////////////////////////
  // END JDH's functions for Magnetic property calculations:
  //////////////////////////////////////////////////////


  
  // Create and Run jobs
  //void CreateGDMAJob(); //JDH NMR code
  //void CreateG09Job( Monomer Monomers[], int NMon ); //JDH NMR code
  //void CreateG09Job( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ); // JDH NMR code 
  //void CreateQChemJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ); // JDH NMR code 
  //void CreateDaltonJob( Monomer Monomers[], int NMon); // JDH NMR code
  //void ModifyGDMAMomFile(); // JDH NMR code
  //string RunGDMAJob(); // JDH NMR code
  //string RunG09Job(); // JDH NMR code
  //string RunG09GDMAJob(); // JDH NMR code
  //string RunFormChk(); // JDH NMR code
  //bool GetUseInEmbedding() { return UseInEmbedding_;}; //JDH NMR code
  //void SetUseInEmbedding(bool INPUT){ UseInEmbedding_ = INPUT;}; //JDH
  //bool GetUseInTwoBodyNMR() { return UseInTwoBodyNMR_;}; //JDH
  //void SetUseInTwoBodyNMR(bool INPUT){ UseInTwoBodyNMR_ = INPUT;}; //JDH
  //bool GetUseInClusterNMRJob() { return UseInClusterNMR_;}; //JDH
  //void SetUseInClusterNMRJob(bool INPUT){ UseInClusterNMR_ = INPUT;}; //JDH

  
  //void ReadDaltonNMRdata(); // JDH NMR code
  //void ReadG09NMRdata(); // JDH NMR code
  //double ReadG09Energy(); // JDH
  //void ReadMolProNMRdata(); //JDH
  //void ReadQChemNMRdata(); // JDH

  void CreateQMJob(Monomer Monomers[], int NMon);
  void CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job=false);
  //void CreateMolProJob(Monomer Monomers[],int NMon);
  //void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon);
  //void CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images); // JDH NMR code
  void CreateFiniteDifferenceMolProJob(Monomer Monomers[],int NMon);
  void CreateFiniteDifferenceMolProJob(Monomer Monomers[],int NMon,bool CCSDT=false);
  void CreateMMJob(Monomer Monomers[], int NMon);
  void CreateTinkerJob(Monomer Monomers[], int NMon);
  void CreateCamCaspJob();  // by Ali
  string RunQMJob(bool CCSDT=false);
  string RunQChemJob(bool MM_job=false);
  string RunFiniteDifferenceMolProJob();
  string RunFiniteDifferenceMolProJobs(int i,bool plus,bool CCSDT);
  string RunMolProJob(bool CCSDT=false);
  string RunMMJob();
  string RunTinkerJob();
  string RunCamCaspJob();  

  //Create and Run PSI4 Jobs // CSG
  void CreatePSI4Job(Monomer Monomers[], int NMon);
  string RunPSI4Job(bool MM_job=false);

  // Read Monomer energies from disk
  void ReadQMResults(Monomer& SymMonomer);
  void ReadMMResults(Monomer& SymMonomer);
  double ReadQChemEnergy(bool MM_job=false);
  double ReadMolProEnergy();
  double ReadPSI4Energy(); // CSG
  void SetQMEnergy(double E) {Energy_QM = E;};
  void SetMMEnergy();
  void SetMMEnergy(double E) {Energy_MM = E;};
  double ReadTinkerEnergy();
  void ReadMultipoleMoments(Monomer& SymMonomer);
  void ReadMultipoleMoments(); //JDH
  void ReadPolarizabilities(Monomer& SymMonomer);
  void ReadIsotropicDipolePolarizability();
  void ReadDispersionCoefficients();
  void ReadFreqPolarizabilities(); // by shuhao
  void ReadDiagonalAnisotropicFreqPolarizabilities(Monomer& SymMonomer);
  void SetEmpiricalAtomDispersionType();

  //don't use Rotation Angles or vectors anymore. Use Rotation Matrix instend
  void SetRotationAngle(double RotAng){RotationAngle = RotAng;};
  void SetRotationVector(int dim, double RotVec){RotationVector[dim] = RotVec;};  // by Ali
  
  double GetQMEnergy() {return Energy_QM;};
  double GetMMEnergy() {return Energy_MM;};

  // Read Monomer gradients and hessian from disk
  void SetQMGradient();
  void SetMMGradient();
  void SetQMHessian();
  double GetQMHess_value(int i,int j){return Hess_QM(i,j);};	//Watit
  void SetMMHessian();
  // or from another monomer
  void SetQMGradient(Monomer& other);
  void SetMMGradient(Monomer& other);
  Vector ReadGradient(string path, int type); 
  void SetMMHessian(Monomer& other);
  Matrix SetSymmetricalHessian(int type,Monomer& Sym_Mon);
  Matrix ReadHessian(string path,int type);
  Matrix ReadFiniteDifferenceMolProHessian(string path);
  Matrix ReadFiniteDifferenceCCSDTHessian(string path);
  // or just from a vector or matrix
  void SetQMGradient(Vector grad) {Grad_QM = grad; QM_Grad_Init=1;};
  void SetMMGradient(Vector grad) {Grad_MM = grad; MM_Grad_Init=1;};
  void SetQMHessian(Matrix hess) {Hess_QM = hess; QM_Hess_Init=1;};
  void SetMMHessian(Matrix hess) {Hess_MM = hess; MM_Hess_Init=1;};
  void SetQMHessian(Monomer& other);



  Vector GetQMGradient() {return Grad_QM;};
  Vector GetMMGradient() {return Grad_MM;};
  Matrix GetQMHessian() {return Hess_QM;};
  Matrix GetMMHessian() {return Hess_MM;};

  Vector FindDistance(const Monomer& Mon);

  // Printing
  void PrintMonomerCartesian(FILE *outfile=stdout);
  void PrintGhostMonomerCartesian(FILE *outfile=stdout);
  // this one lets you substitute an alternative atomic symbol for printing,
  // i.e. to change color of key fragments when visualizing
  void PrintMonomerCartesian(FILE *outfile, string symbol); 
  void PrintMolProEmbeddingCharges(FILE* outfile); //JDH NMR code
  void PrintMolProMonomerCastesian(FILE *outfile,int start);
  void PrintMolProMonomerCastesian(FILE *outfile,Vector Coords,int start);
  void PrintPSI4MonomerCartesian(FILE *outfile,int start); //CSG
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintTinkerCartesian(int shift=0, bool monomer=true, 
			    FILE *outfile=stdout);
  void PrintTinkerCartesian(Vector coords, int shift=0, bool monomer=true, 
  			    FILE *outfile=stdout);
  void PrintQChemEmbeddingCharges(FILE* outfile=stdout);
  void PrintTinkerEmbeddingCharges(int shift=0, FILE* outfile=stdout);

  void PrintQMGradient(string title);
  void PrintMMGradient(string title);
  void PrintGradient(string title, Vector grad);
  void PrintAll();

  // Modify monomer geometry
  void Translate(Vector new_com);
  void ShiftCoords(Vector MonomerRModeVector); //Watit
  void FractionalTranslate(int x,int y, int z);
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

  // Watit - IR & Raman Intensities
	void ReadHess();
	void SetHessQM(Matrix in) {HessQM = in;};
	void SetHessMM(Matrix in) {HessMM = in;};
	Matrix GetHessQM_Matrix(){return HessQM;};
	Matrix GetHessMM_Matrix(){return HessMM;};
	double GetHessQM_value(int i, int j){return HessQM(i,j);};
	double GetHessMM_value(int i, int j){return HessMM(i,j);};

	void ReadIntensity();
	void SetDipD(Matrix in) {DipD = in;};
	Matrix GetDipD_Matrix(){return DipD;};
	double GetDipD_value(int i,int j){return DipD(i,j);};
	void SetPolD(Matrix in) {PolD = in;}; 
	Matrix GetPolD_Matrix(){return PolD;};
	double GetPolD_value(int i,int j){return PolD(i,j);};
  //

};

#endif
