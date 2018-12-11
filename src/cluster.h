#ifndef _cluster_h
#define _cluster_h
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
#include "nmr.h" // JDH
#include "ee.h" // JDH
#include "input.h" // JDH
#include "vector.h" // my vector class, "Vector"
#include <vector> // C STL library vector, "vector"
using std::vector;
#include <time.h>

/*
  A cluster class.  A cluster is composed of NMon monomers, and it has
NDim = NMon*(NMon-1) dimers.  Other relevant properties are contained, as
well.  

This class also handles most of the HMBI calculations.

Unit conventions:
  Coordinates in Angstroms
  Energies in hartrees
  Gradients in hartrees/bohr
  Hessians in hartrees/bohr^2


GJB 01/09

Made into a singleton class: only a single instance of it can exist.
This is a hack to work with the DL-FIND libraries.  

 */


class Cluster {
 private:

  Cluster();
  ~Cluster();

  string Title; //job name, from $comment section
  string InputFilename;
  string *InputFile;
  int NumberOfLinesInInput;
  int NMon, NDim,NDim_trunc,NDimImages_trunc; // number of monomers, dimers
  // Pointers for arrays to store info about each monomer.  Size NMon+1.
  // Element 0 = full cluster, 1->NMon = monomer info
  int *spins, *charges, *Natoms_per_monomer; 
  string *types;

  int UniqueAtoms,UniqueMon,UniqueDim;//# of symmetry unique defining Atom, 
                                      //Monomers and Dimers

  //  Arrays for Monomers and Dimers.  In each case, the arrays
  // count from 1->NMon/NDim
  Monomer *Monomers;
  Dimer *Dimers;

  // for periodic systems:
  int NMon_images;
  Monomer *MonomerImages;
  int NDim_images;
  Dimer *DimerImages;
  int TotalDim_nosym;//# of dimers with no symmetry
  
  double TotalMass;
  Vector CenterOfMass;
  Vector AtomicCoordinates; // Vector of length 3*GetTotalNumberOfAtoms() 
  Vector AcceptedAtomicCoordinates; // Vector of length 3*UniqueAtoms(); is the last accepted step in an optimization
  Vector SymmetryUniqueCoordinates;// Cartesian coordinates of the symmetry unique atoms
  Vector FractionalCoordinates;// Fractional Coordinates of the atoms in the unit cell
  Vector SymmetricalFractCoord;// Fractional Coordinates of the symmetry unique atoms 

  string *AtomicSymbols;
  int *AtomicNumbers;
  double *AtomicMasses;
  string *DispersionAtomTypes;
  Vector AtomTypes;

  Vector* unit_cell; // stores unit cell.  First dim loops over vectors,
                          // and second over xyz coordinates for each vector.
  Vector* reciprocal_cell;// stores reciprocal cell.

  // Also store angle/axis form;
  Vector UnitCellAngles;
  Vector UnitCellAxes;
  Vector AcceptedUnitCellAngles;
  Vector AcceptedUnitCellAxes;  
  // store the crystal volume
  double cell_volume;

  //lattice parameters locked under symmetry
  bool b_locked_by_a;//b and a will be the same length during optimization if true 
  bool c_locked_by_a;//c and a will be the same length during optimization if true 
  bool c_locked_by_b;//c and b will be the same length during optimization if true 
  bool beta_locked_by_alpha;//beta and alpha will always be the same during optimization if true 
  bool gamma_locked_by_alpha;//gamma and alpha will always be the same during optimization if true 
  bool gamma_locked_by_beta;//gamma and beta will always be the same during optimization if true
  bool lock_alpha;//alpha is fixed, used if alpha is 90 or 120 degrees
  bool lock_beta;//beta is fixed, used if beta is 90 or 120 degres
  bool lock_gamma;//gamma is fixed, used if gamma is 90 or 120 degrees


  double Energy_QM; // only used if doing benchmark full cluster QM calc.
  double E_Electrostatic_MM, E_Induction_MM; // stored in hartrees
  double Lattice_E_Electrostatic_MM, Lattice_E_Induction_MM; // stored in hartrees
  double Energy_MM;
  double Energy_Vib; //vibrational energy
  double Entropy;//vibrational entropy
  double Cv;//constant volume heat capacity
  double gruneisen_parameter;//macroscopic gruneisen parameter
  double Energy_HMBI; // stored in hartrees
  double Energy_HMBI_old; // stored in hartrees
  // Note: Gradients are stored in hartrees/bohr
  Vector Grad_QM; //JLM 
  Vector Grad_MM; 
  Matrix Hess_QM; //JLM
  Matrix Hess_MM;
  Vector Grad_HMBI; 
  Matrix Hess_HMBI; 
  Matrix UMWHess_HMBI; // Watit 
  Vector Freqs;//in cm^-1
  Matrix Vib; // Watit
  int* Grad_key; // list of indexing for each monomer in full gradient

  bool Do_Forces; // whether or not to load force files
  int QM_Grad_Init, MM_Grad_Init, QM_Hess_Init, MM_Hess_Init;

  //Quasiarmonic Approximation, these are now handled by the Quasiharmonic class
  //double reference_volume;
  //Vector reference_frequencies;
  //Vector Gruneisen_parameters;
  //bool Gruneisen_init;


  double HartreeToKJpMM; // conversion factor, hartree -> kJ/mol per monomer

  time_t start_time; //used to determine the time that elapsed

  // Watit - IR & Raman Intensities
  	Vector IntIR;
  	Vector IntRaman;
  	int Nzero;
  //

 public:
  // instantiate the single instance of this class
  static Cluster& cluster() {
    static Cluster thecluster;
    return thecluster;
  }

 void Initialize(ifstream& infile, int Nproc, string input_filename);

  // Print program introduction
  void PrintHeader();

  // Reset all the parameters
  void ResetEnergiesAndGradients() {
    Energy_HMBI = 0.0; Energy_MM = 0.0; Energy_QM = 0.0;
    Grad_MM.Set(); Grad_QM.Set();
    Grad_HMBI.Set();
  };


  // Input File processing routines
  void StoreInputFilename(string input_file);
  void CreateWorkingDirectories();
  void CreateWorkingDirectories(string Dir_name,bool HessJob);
  int FindNumberOfMonomers(ifstream& infile);
  string GetJobTitle(ifstream& infile);
  void PrintInputFile(ifstream& infile);
  void StoreInputFile(ifstream& infile);
  void CountNumberOfAtoms(ifstream& infile);
  void ReadHMBIParameters(ifstream& infile);
  void ReadAIFFParameters(ifstream& infile);  // by Ali
  void PrintExternalPrograms(); // print out paths of external programs
  string GetInputFilename() {return InputFilename;};

  void ReadGeometry(ifstream& infile);
  void ReadTinkerGeometry(ifstream& infile);
  void ReadSimpleXYZGeometry(ifstream& infile);
  Vector ReadEmbeddingCharges(ifstream& infile);
  void ReadUnitCellDimensions(ifstream& infile);
  void OldReadUnitCellDimensions(ifstream& infile);
  double GetUnitCellParameter(string type);
  void SetUnitCellParameter(string type, double value);
  void SetUnitCellVectors(Vector v1, Vector v2, Vector v3);
  Vector GetUnitCellVector(int index) {return unit_cell[index];};
  Vector GetReciprocalUnitCellVector(int index) {return reciprocal_cell[index];};
  Vector GetUnitCellAxes() {return UnitCellAxes;};
  Vector GetUnitCellAngles() {return UnitCellAngles;};
  Vector GetLastAcceptedUnitCellAxes() {return AcceptedUnitCellAxes;};
  Vector GetLastAcceptedUnitCellAngles() {return AcceptedUnitCellAngles;};
  Vector ComputeCellParametersFromLatticeVectors(Matrix latticeVec);
  Matrix ComputeLatticeVectorsFromCellParameters(Vector cellParams);

  void FindCenterOfMass();

  string ReadCIFSection(ifstream& infile); // Watit
  void CreateInputfromCIF(string ciffile);
  void RemoveRedundanceCoord(int NAtom, string Element[], double Position[]);

  string ReadQChemRemSection(ifstream& infile, int type=1);
  
  string ReadQChemBasisSection(ifstream& infile); // by shuhao

  // PSI4 Stuff by CSG
  string ReadPSI4Section(ifstream& infile); //CSG
  //PSI4 rem section //CSG
  string ReadPSI4RemSection(ifstream& infile); 

  // NMR Stuff by JDH
  string ReadDaltonSection(ifstream& infile); //JDH
  string ReadG09Section(ifstream& infile); //JDH + Watit
  string ReadOrcaSection(ifstream& infile); //JDH
  void PE_potential(int k, int type); //JDH
  //void PE_molecule(int k, bool dimer); //JDH
  //void PE_potential(int k, int type); //JDH

  // JDH Custom Basis Stuff:
  string ReadHBasis(ifstream& infile); //JDH
  string ReadCBasis(ifstream& infile); //JDH
  string ReadNBasis(ifstream& infile); //JDH
  string ReadOBasis(ifstream& infile); //JDH
  string ReadSBasis(ifstream& infile); //JDH
  string ReadClBasis(ifstream& infile); //JDH
  string ReadIBasis(ifstream& infile); //JDH
  string ReadSnBasis(ifstream& infile); //JDH
  string ReadPBasis(ifstream& infile); //JDH
  string ReadKBasis(ifstream& infile); //JDH
  string ReadNaBasis(ifstream& infile); //JDH
  string ReadBrBasis(ifstream& infile); //JDH
  string ReadFBasis(ifstream& infile); // JDH
  string ReadWBasis(ifstream& infile); // JDH
  string ReadVBasis(ifstream& infile); // JDH

  //MolPro rem sections
  string ReadMolProSection(ifstream& infile); // by yoni
  string ReadMolProInstSection(ifstream& infile); // by yoni
  string ReadMolProCBSBasisSection(ifstream& infile); // by yoni
  string ReadMolProInstHFSection(ifstream& infile); // by yoni
  string ReadMolProCCSDTBasisSection(ifstream& infile);//by yoni
  string ReadMolProCCSDTMP2Section(ifstream& infile);//by yoni
  string ReadMolProCCSDTInstSection(ifstream& infile);//by yoni

  string ReadQuantumEspressoSpeciesSection(ifstream& infile); // JLM
  string ReadQuantumEspressoSupercellSection(ifstream& infile); // JLM
  string ReadKPoints(ifstream& infile); // JLM
  string ReadSlaterKoster(ifstream& infile); // JLM
  string ReadHubbardDerivatives(ifstream& infile); // JLM
  string ReadMaxAngularMomentum(ifstream& infile); // JLM
  string ReadTinkerRemSection(ifstream& infile);
  string ReadAIFFRemSection(ifstream& infile);  // by Ali

  //Crystal rem sections
  string ReadCrystalHeadingSection(ifstream& infile);//by yoni
  string ReadCrystalBasisSection(ifstream& infile);//by yoni
  string ReadCrystalEndingSection(ifstream& infile);//by yoni

  void Rewind(ifstream& infile);

  // Other routines
  int GetNumberOfMonomers() {return NMon;};
  int GetNumberOfMonomerImages() {return NMon_images;}; //JDH
  int GetNumberOfDimers() {return NDim;};
  int GetNumberOfDimerImages() {return NDim_images;};
  int GetNumberofSymmetryUniqueDimers() {return UniqueDim;};
  int GetTotalNumberOfAtoms() {return Natoms_per_monomer[0];};
  int GetTypesOfAtoms(); //JLM
  int GetNumberOfAtoms(int iMon) {return Natoms_per_monomer[iMon];};
  int GetClusterSpin() {return spins[0];}; //by yoni
  int GetMonomerSpin(int iMon) {return spins[iMon];}; //by yoni
  int GetClusterCharge() {return charges[0];}; //by yoni
  int GetMonomerChange(int iMon) {return charges[iMon];}; //by yoni
  int GetNumberOfUniqueAtoms() {return UniqueAtoms;};
  int GetMonomerSymmetryFactor(int iMon) {return Monomers[iMon].GetSymmetryFactor();};
  void SetAtomicSymbols();
  string GetAtomicSymbol(int i) {return AtomicSymbols[i];};
  int GetAtomicNumber(int i) {return AtomicNumbers[i];};
  string GetDFTBGeomInput();
  double GetTotalMass() {return TotalMass;};
  
  void ReadCurrentCoordinates();
  Vector GetCurrentCoordinates(bool GetLatticeParams=true);// {return AtomicCoordinates;};
  Vector GetLastAcceptedCoordinates(bool GetLatticeParams=true);

  void CorrectingDimerRotations(int iDim,bool image, bool NonImagetoImage,int ListNumb);
  void RedetermineDimerSymmetry(int thisDim,Matrix Rot2,bool image,bool NonImagetoImage,int ListNumb);
  void CorrectingMonomerRotations(int iMon);
  bool CheckMonomerRotationOperator(int iMon, Matrix Rot);
  bool CheckDimerRotationOperator(int iDimer,Matrix Rot2,Matrix WrongMat,bool image);
  void SetFractionalCoordinates();
  Vector GetFractionalCoordinates(bool IncludeLatParams=true); //{return FractionalCoordinates;};
  void PrintFractionalCoordinates(FILE *outfile=stdout, bool includeLatParams=true);
  void ReadSymmetryUniqueCoordinates();
  Vector GetSymmetryUniqueCoordinates(bool IncludeParams=true);// {return SymmetryUniqueCoordinates;};
  Matrix GetChangeInRotationFromLatticeVectors(int Mon,int ChangeValue,bool FromUniqueToThis );
  void DetermineMonomerFractionalSymmetry();
  void ReadSymmetryUniqueFractionalCoordinates(); 
  Vector GetSymmetryUniqueFractionalCoordinates(bool IncludeLatticeParams=true);//{return SymmetricalFractCoord;};
  void SetNewCoordinates(Vector NewCoords, bool ResetEnergies=true);
  Vector GetSymmetryImposedCoordinates(Vector reduced_coord,bool usefract,bool GradIncludeLattice=true);

  Vector RearrangeCoordsForConstraintOpt(Vector coords);
  Vector RearrangeBackCoordsForConstraintOpt(Vector rearranged_coords);
 
  void MakeOptimizationStepFile();
  Vector EnergyOfFiniteDifference(double delta,bool IsPositive);
  double GetCellvolume() {return  cell_volume;};
  void SetCellvolume(double CellVol){cell_volume = CellVol;};
  // Return index for dimer composed of two monomers
  int DimerLookup(int imon1, int imon2); 


  // JDH MAGNETIC PROPERTIES CODE
  void PredictMagneticProperties();
  void DetermineSelfConsistentElectrostaticEmbeddingEnvironment();
  void ComputeSelfConsistentEmbeddingCharges();
  void AssignMonomerLists(); // FOR use with 'create_jobs_only'
  void CreateFragmentJobs();
  void RunJobsAndComputeEnergy(); // JDH NMR code


  
  //Creating and Running QM and MM jobs
  /* void RunGDMAJobs(); // JDH NMR code */
  /* void RunJobsAndComputeEnergy(); // JDH NMR code */
  /* void RunJobsAndComputeNMRShieldingTensors(); // JDH NMR code */
  /* void ComputeNMRShieldingTensors(); // JDH NMR code */
  /* void CreateQMNMRJobs(); // JDH NMR code */
  /* void RunQMNMRJobs(); // JDH NMR code */
  /* string* GetNMRJobList(int Njobs); // JDH NMR code */
  /* void ReadNMRData(); // JDH NMR code */
  void RunJobsAndComputeEnergy(int argc, char **argv);
  void CreateQMJobs();
  void CreateMMJobs();
  void CreateFiniteDifferenceDimers(bool CCSDT);
  void CreateFullClusterQChemJob(bool MM_job=false);
  void CreateFullClusterTinkerJob(int path, bool DiffName, string Name);
  void CreateQuantumEspressoJob();
  void CreateDFTBJob();
  void CreateFiniteDifferenceTinkerJobs(int path);
  void CreateCrystalJob();
  void CreateFiniteDifferenceCrystalJob();
  void CreateOrientJob();  // by Ali

  void RunJobs();
  string* RunQMJobs(int Njobs);
  string* RunMMJobs();
  string RunQuantumEspressoJob();
  string RunDFTBJob();
  string RunFullClusterQChemJob(bool MM_job=false);
  string RunFullClusterTinkerJob();
  string RunFiniteDifferenceTinkerJob(int CoordEntry,bool IsPlus);
  string RunCrystalJob();
  string RunFiniteDifferenceCrystal(int num,bool Plus);
  string RunOrientJob();

  // Read the Results
  void ReadQMResults();
  void ReadMMResults();
  void SetQMEnergy(); // for full cluster qm energy JLM
  void SetMMEnergy(); // for full cluster mm energy  

  // Energy routines
  double ReadQChemEnergy(bool MM_job=false);
  double ReadTinkerEnergy(string filename = "full.out");
  double ReadOrientEnergy();
  double ReadCrystalEnergy();
  double ReadQuantumEspressoEnergy();
  double ReadDFTBEnergy();

  double GetMMEnergy() {return Energy_MM;};
  double GetQMEnergy() {return Energy_QM;};
  double GetHMBIEnergy() {return Energy_HMBI;};
  double GetPreviousHMBIEnergy() {return Energy_HMBI_old;};
  void SetHelmholtzEnergy(double energy) {Energy_Vib = energy;};
  double GetHelmholtzEnergy() {return Energy_Vib;};
  void SetEntropy(double S) { Entropy = S;}
  double GetEntropy() {return Entropy;}
  void SetCv(double heat_capacity){Cv = heat_capacity;};
  double GetCv(){return Cv;};
  void SetGruneisenParameter(double gamma){gruneisen_parameter = gamma;};
  double GetGruneisenParameter(){return gruneisen_parameter;};
  double ComputeHMBIEnergy();
  double GetTotalOneBodyEnergy(string type);
  double GetTotalTwoBodyEnergy(string type);

  // add a 2-body contribution to the gradient
  void AddGradientContribution(int indA, int indB, Vector dimerGrad);

  // Gradient routines
  Vector ComputeHMBIGradient();
  Vector GetHMBIGradient() {return Grad_HMBI;};
  void SetQMGradient(); //JLM
  void SetMMGradient();
  Vector ReadGradient(int type); 
  Vector ReadQuantumEspressoGradient(); //JLM
  Vector ReadDFTBGradient(); //JLM
  Vector ReadCrystalGradient();
  Vector ReadFiniteDifferenceCrystalGradient();
  void PrintGradient(string title, const Vector& Grad);
  

  // Frequency routines
  void ComputeHarmonicFrequencies(Matrix& Hessian);
  Matrix ComputeMassWeightedHessian(Matrix& Hessian);
  Matrix UnMassWeightedHessian(Matrix& Hessian);
  Matrix UnMassWeightedEigenVectors(Matrix& Hessian);
  Matrix ProjectOutTranslationAndRotation(Matrix& Hessian);
  Matrix ComputeInertiaTensor();
  double ComputeZeroPointEnergy(Vector frequncies);
  void ComputeVibrationalContribution(Vector frequencies);
  void PrintFrequencies();


  void PrintNormalModes(Vector freqs, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index);
  void PrintNormalModes(Vector freqs, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index,string Name);
  Vector GetFrequency() {return Freqs;};

  //QuasiHarmonic Approximation routines
  void PrintVibrationalFrequencies(Vector freq,string Header, bool RemoveOldFile);
  //void ReadQuasiHarmonicFrequencies(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus);
  //void SetReferenceVolume(double volume) {reference_volume = volume;};
  //double GetReferenceVolume() {return reference_volume;};
  //void SetReferenceFrequencies(Vector frequencies) {reference_frequencies = frequencies;};
  //Vector GetReferenceFrequencies() {return reference_frequencies;};
  //void SetGruneisenParameters(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus);
  //Vector GetGruneisenParameters() {return Gruneisen_parameters;};
  //void InitializeQuasiHarmonicForSuperCell();
  void ComputeQuasiHarmonicFrequencies();
  Matrix ComputeQuasiHarmonicGradient();

  Matrix ComputeHMBIHessian();
  void PrintHessian(string title, Matrix hess);
  void SetQMHessian(); //JLM
  void SetMMHessian();
  void ReadQuantumEspressoHessian(); //JLM
  void ReadDFTBHessian(); //JLM
  Matrix ReadHessian(int type);
  Matrix ReadCrystalFrequencies();
  double GetEnthalpyPV();
  Matrix ReadFiniteDifferenceHessian();
  Matrix CreateAndRunFiniteDifferenceHessian();
  Matrix GetHMBIHessian() {return Hess_HMBI;};
  Vector CompareSortedOrigVecs(Vector orig, Vector sorted);
  Matrix SortColumns(Matrix orig,Vector sorted);
  Matrix GetFiniteDifferenceTinkerHessian(Vector Original_coords);

  Monomer GetMonomer(int iMon) {return Monomers[iMon];};
  Dimer GetDimer(int i) {return Dimers[i];};
  Dimer GetDimerImage(int i) {return DimerImages[i];};
  double GetAtomicMasses(int i) {return AtomicMasses[i];};

  // Watit - IR & Raman Intensities
  	Matrix VPolD;
	Matrix IRIntensity;
  	Matrix RamanIntensity;
  	Matrix RMode;
  //

  // Printing
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintCartesianCoordinates(FILE *outfile=stdout);
  void PrintTinkerCartesian(FILE *outfile=stdout);
  void PrintTinkerCartesian(Vector coords, FILE *outfile=stdout);
  void PrintXYZCartesian(FILE *outfile=stdout);
  void PrintInputFile(FILE *outfile=stdout);
  void OldPrintInputFile(FILE *outfile=stdout);

  // Write out crystal xyz coords for N1xN2xN3 cells
  void WriteCrystalXYZFile(int N1, int N2, int N3); 

  // Modify Cluster geometry
  void AdjustIntermolecularSpacing(double factor);
  void AdjustCoord(int mode, Matrix RModeVector); // Watit

  // Track geometry with xyz file
  void UpdateTrajectoryFile(int step, bool new_trajectory=false);

  void ComputeDistanceMatrix();
 
  void CreatePeriodicImageDimerList();
  void CreatePeriodicImageDimerList(double r_cutoff_2bd); //JDH
  void CreatePeriodicImageMonomerList(double r_cutoff);

  void UpdateJobStatus(int ijob);
  void OldUpdateJobStatus(int ijob);

  void DebugLatticeGradient();
  Vector ComputeMMCrystalGradient();
  Matrix ComputeStressTensor(Matrix LatticeGradient);

  void ComputeLatticeParamGrad(Vector& Grad);
  void ComputeLatticeParamGradFullMM();
  //Vector CartFromFractAndLatticeParams(Vector FractCoords,bool reverse,double a,double  b,double c,
  //				       double alpha,double beta,double gamma);
  Vector CartFromFractAndLatticeParams(Vector FractCoords,bool reverse);
  Vector ReadLatticeParamGradFullMM(Vector grad_lat_full_MM);
  double UpdateLatticeParamsAndComputeEnergy(int i, Vector coords, double delta);
  void UpdateLatticeParams(Vector coords); 
  void UpdateGeometryAndLatticeParams();
  void UpdateQuantumEspressoLattice();//JLM
  void UpdateQuantumEspressoGeometry();//JLM
  Vector LatticeParamsGradFromVectorsGrad(Vector VectorsGrad);
  void MaintainCartesianSymmetry(Vector& Coord,bool LatticeIncluded);
  Vector UniqueMonomerCoord(Vector& Coord,bool LatticeIncluded);
  Vector ConvertBetweenFractionAndCartesianCoordinates(Vector Coord, bool InFractional, bool GradIncludeLattice,bool Transpose = false);

  void AlterLatticeVectors(int VectorEntry,double ChangedBy);
  void ChangeCellVolume();

  void TranslateUsingLatticeParams(Vector OldUnitCellAngles,Vector OldUnitCellAxes);
  double New_ComputeLatticeParamGradFullMM(int i, double delta);
  double AssignGradContributionSignAlpha(int i, int j);
  double AssignGradContributionSignGamma(int i, int j);
  double AssignGradContributionSignBeta(int i, int j);
  double AssignGradContributionSignGammaExtra(int i, int j);
  void PrintUnitCellParams();

  double UnitCellVolume();

  //lattice parameters locked under symmetry
  bool IsBLocked() {return b_locked_by_a;};
  bool IsCLocked() {
    if(c_locked_by_a ||c_locked_by_b) return 1;
    else return 0;};
  bool IsAlphaLocked() {return lock_alpha;};
  bool IsBetaLocked() {
    if (lock_beta ||beta_locked_by_alpha) return 1;
    else return 0;}
  bool IsGammaLocked() {
    if(lock_gamma ||gamma_locked_by_alpha ||gamma_locked_by_beta)
      return 1;
    else return 0;};

  void ReadFrozenUnitCellParams(ifstream& infile);
  //JLM being lazy here...
  bool IsFrozenUnitCellParam(int i); 
  bool FreezeLatParams[6];

  Matrix ComputeInverseUnitCellMatrix();
  Matrix FormExternalStressTensor();
  Matrix FormGradientExternalPressureTerm();
  Vector ComputeVibrationalGradient();

  /*** AIFF functions ***/

  double ComputeAIFFEnergy(); // Main driver, both regular & periodic systems

  /* Electrostatics & Induction */
  double ComputeClusterMultipoleInteractions();  // Electrostatics + Induction
  Vector ComputeClusterMultipoleGradient();  // ES + Ind, energy & gradient
  // Used to preconverge the non-periodic multipole moments before
  // doing periodic case:
  void PreconvergeInducedMultipoleMomentsInUnitCell(Multipole dQ[], Multipole old_dQ[], int multipole_key[]); 
  double ComputePeriodicMultipoleInteractions();  // Periodic ES + Ind

  /* 2-Body Dispersion */
  // note: The periodic version calls the normal one to get interactions
  // within the unit cell.
  double ComputeTwoBodyDispersion(); // Normal
  double ComputePeriodicTwoBodyDispersion(); // Periodic systems

  /* 3-Body Dispersion */
  // Main callable routine for 3-body dispersion:
  double ComputeManyBodyDispersion(); 
  // Handle the calculating of dispersion coefficients, looping:
  double ComputeThreeBodyDispersion(string type = "default"); 
  // Evaluate the actual dispersion contribution:
  double ComputeAxilrodTellerMutoDispersionTerm(Atom I, Atom J, Atom K, string type = "default");



  // Old empirical Tkatchenko-Scheffler dispersion functions
  void ReadDispersionAtomTypes(ifstream& infile);
  void ReadDispersionCoefficients(ifstream& infile);
  double ComputeEmpiricalTwoBodyDispersion(); 
  // 3-body dispersion in this model is done via an altered call to
  // ComputeThreeBodyDispersion, using type = "Tkatchenko"



  //Used to determine the time that elapsed
  void SetStartTime(time_t time) {start_time = time;}
  time_t GetStartTime() {return start_time;};

  void ReadHessianData();

  // Watit - IR & Raman Intensities
	void ReadIntensityData();
  	void SetRMode(Matrix in) {RMode = in;};
  	Matrix GetRMode() {return RMode;};
  	double GetRMode_value(int i,int j) {return RMode(i,j);};
  	//void SetVPolD(Matrix in) {VPolD = in;};
  	//Matrix GetVPolD_Matrix(){return VPolD;};
  	//double GetVPolD_value(int i,int j){return VPolD(i,j);};
	void SetIntIR(Vector in) {IntIR = in;};
	//Vector GetIntIR(){return IntIR;};
  	void SetIntRaman(Vector in) {IntRaman = in;};
  	//Vector GetIntRaman(){return IntRaman;};
  	void SetVib(Matrix in) {Vib = in;};
	void PlotIntensityData(int NVib);
  //

};

#endif
