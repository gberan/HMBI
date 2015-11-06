#ifndef _cluster_h
#define _cluster_h
#include <iostream>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream.h>
#include <cassert>
#include <iomanip.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "monomer.h"
#include "dimer.h"
#include "params.h"
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



GJB 01/09

Made into a singleton class: only a single instance of it can exist.
This is a hack to work with the DL-FIND libraries.  

 */


class Cluster {
 private:

  Cluster();
  ~Cluster();

  string Title; //job name, from $comment section

  int NMon, NDim; // number of monomers, dimers
  // Pointers for arrays to store info about each monomer.  Size NMon+1.
  // Element 0 = full cluster, 1->NMon = monomer info
  int *spins, *charges, *Natoms_per_monomer; 

  //  Arrays for Monomers and Dimers.  In each case, the arrays
  // count from 1->NMon/NDim
  Monomer *Monomers;
  Dimer *Dimers;

  // for periodic systems:
  int NMon_images;
  Monomer *MonomerImages;
  int NDim_images;
  Dimer *DimerImages;

  double TotalMass;
  Vector CenterOfMass;
  Vector AtomicCoordinates; // Vector of length 3*GetTotalNumberOfAtoms() 

  string *AtomicSymbols;
  int *AtomicNumbers;
  double *AtomicMasses;
  string *DispersionAtomTypes;


  Vector* unit_cell; // stores unit cell.  First dim loops over vectors,
                          // and second over xyz coordinates for each vector.
  Vector* reciprocal_cell;// stores reciprocal cell.


  // Also store angle/axis form;
  Vector UnitCellAngles; // in degrees
  Vector UnitCellAxes; // in Angstroms
  // store the crystal volume
  double cell_volume;

  double Energy_QM; // only used if doing benchmark full cluster QM calc.
  double E_Electrostatic_MM, E_Induction_MM; // stored in hartrees
  double Lattice_E_Electrostatic_MM, Lattice_E_Induction_MM; // stored in hartrees
  Multipole *InducedMultipoles;
  // Used in Ewald sum
  Matrix *ReciprocalTabs;
  Matrix *DirectTabs;
  Matrix *DeltaDampTabs;

  double Energy_MM;
  double Energy_HMBI; // stored in hartrees
  // Note: Gradients are stored in hartrees/bohr
  Vector Grad_MM; 
  Vector Grad_HMBI; 
  int* Grad_key; // list of indexing for each monomer in full gradient

  bool Do_Forces; // whether or not to load force files
  int QM_Grad_Init, MM_Grad_Init;

  double HartreeToKJpMM; // conversion factor, hartree -> kJ/mol per monomer

 public:
  // instantiate the single instance of this class
  static Cluster& cluster() {
    static Cluster thecluster;
    return thecluster;
  }

  void Initialize(ifstream& infile, int Nproc);

  // Print program introduction
  void PrintHeader();

  // Reset all the parameters
  void ResetEnergiesAndGradients() {
    Energy_HMBI = 0.0; Energy_MM = 0.0;
    Grad_MM.Set();
    Grad_HMBI.Set();
  };


  // Input File processing routines
  void CreateWorkingDirectories();
  int FindNumberOfMonomers(ifstream& infile);
  string GetJobTitle(ifstream& infile);
  void PrintInputFile(ifstream& infile);
  void CountNumberOfAtoms(ifstream& infile);
  void ReadHMBIParameters(ifstream& infile);
  void ReadAIFFParameters(ifstream& infile);  // by Ali
  void PrintExternalPrograms(); // print out paths of external programs

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
  void ConvertFractionalToCartesianCoordinates();


  void FindCenterOfMass();

  string ReadQChemRemSection(ifstream& infile, int type=1);
  string ReadQChemBasisSection(ifstream& infile); // by shuhao
  string ReadTinkerRemSection(ifstream& infile);
  string ReadAIFFRemSection(ifstream& infile);  // by Ali

  void Rewind(ifstream& infile);

  // Other routines
  int GetNumberOfMonomers() {return NMon;};
  int GetTotalNumberOfAtoms() {return Natoms_per_monomer[0];};
  int GetNumberOfAtoms(int iMon) {return Natoms_per_monomer[iMon];};
  void SetAtomicSymbols();
  string GetAtomicSymbol(int i) {return AtomicSymbols[i];};
  int GetAtomicNumber(int i) {return AtomicNumbers[i];};

  double GetTotalMass() {return TotalMass;};

  void ReadCurrentCoordinates();
  Vector GetCurrentCoordinates() {return AtomicCoordinates;};
  void SetNewCoordinates(Vector NewCoords, bool ResetEnergies=true);

  double GetCellvolume() {return  cell_volume;};
  void SetCellvolume(double CellVol){cell_volume = CellVol;};
  // Return index for dimer composed of two monomers
  int DimerLookup(int imon1, int imon2); 


  //Creating and Running QM and MM jobs
  void RunJobsAndComputeEnergy();
  void RunJobsAndComputeEnergy(int argc, char **argv);
  void CreateQMJobs();
  void CreateMMJobs();
  void CreateFullClusterQChemJob(bool MM_job=false);
  void CreateFullClusterTinkerJob();
  void CreateOrientJob();  // by Ali

  void RunJobs();
  string* RunQMJobs();
  string* RunMMJobs();
  string RunFullClusterQChemJob(bool MM_job=false);
  string RunFullClusterTinkerJob();
  string RunOrientJob();

  // Read the Results
  void ReadQMResults();
  void ReadMMResults();
  void SetMMEnergy(); // for full cluster energy  

  // Energy routines
  double ReadQChemEnergy(bool MM_job=false);
  double ReadTinkerEnergy(string filename = "full.out");
  double ReadOrientEnergy();

  double GetMMEnergy() {return Energy_MM;};
  double GetHMBIEnergy() {return Energy_HMBI;};
  double ComputeHMBIEnergy();
  double GetTotalOneBodyEnergy(string type);
  double GetTotalTwoBodyEnergy(string type);

  // add a 2-body contribution to the gradient
  void AddGradientContribution(int indA, int indB, Vector dimerGrad);

  // Gradient routines
  Vector ComputeHMBIGradient();
  Vector GetHMBIGradient() {return Grad_HMBI;};
  void SetMMGradient();
  Vector ReadGradient(int type); 
  void PrintGradient(string title, const Vector& Grad);

  // Frequency routines
  void ComputeHarmonicFrequencies(Matrix& Hessian);
  Matrix ComputeMassWeightedHessian(Matrix& Hessian);
  Matrix ProjectOutTranslationAndRotation(Matrix& Hessian);
  Matrix ComputeInertiaTensor();
  double ComputeZeroPointEnergy(Vector frequencies);

  // Printing
  void PrintQChemCartesian(FILE *outfile=stdout);
  void PrintTinkerCartesian(FILE *outfile=stdout);
  void PrintXYZCartesian(FILE *outfile=stdout);
  void PrintInputFile(FILE *outfile=stdout);

  // Write out crystal xyz coords for N1xN2xN3 cells
  void WriteCrystalXYZFile(int N1, int N2, int N3); 

  // Modify Cluster geometry
  void AdjustIntermolecularSpacing(double factor);

  // Track geometry with xyz file
  void UpdateTrajectoryFile(int step, bool new_trajectory=false);

  void ComputeDistanceMatrix();
 
  void CreatePeriodicImageDimerList();
  void CreatePeriodicImageMonomerList(double r_cutoff);

  void UpdateJobStatus(int ijob);
  void OldUpdateJobStatus(int ijob);

  void DebugLatticeGradient();
  Vector ComputeMMCrystalGradient();
  Matrix ComputeStressTensor(Matrix LatticeGradient);


  /*** AIFF functions ***/

  double ComputeAIFFEnergy(); // Main driver, both regular & periodic systems

  /* Electrostatics & Induction */
  double ComputeClusterMultipoleInteractions();  // Electrostatics + Induction
  double ComputeClusterAIFFElectrostatics(); // Just electrostatics
  double ComputeClusterAIFFInduction(); // Faster matrix-based approach for induction.
  Vector ComputeClusterMultipoleGradient();  // ES + Ind, energy & gradient

  // Used to preconverge the non-periodic multipole moments before
  // doing periodic case:
  void PreconvergeInducedMultipoleMomentsInUnitCell(Multipole dQ[], Multipole old_dQ[], int multipole_key[]); 
  double ComputePeriodicMultipoleInteractions();  // Periodic ES + Ind
  void ComputePeriodicAIFFInducedMultipoles(); // Faster matrix-based approach for induction.
  double ComputeEwaldDampingCorrection(); // computed difference between damped and undamped induction in periodic cluster.


  // Find the direct/reciprocal space cutoffs and kappa parameter used in the Ewald summation
  double DetermineDirectAndRecipSpaceCutoffs(int &nX, int &nY, int &nZ, int &kX, int &kY, int &kZ); 
  double ComputePeriodicAIFFElectrostatics(); // GJB cleaned-up multipolar Ewald
  void ComputePeriodicAIFFElectrostaticsAndInduction(); 
  void ComputePeriodicAIFFInducedMultipolesMadelung(); // Uses Madelung Potential instead of finite cluster for multipoles.
  Vector ComputeMadelungPotentialForInduction(bool include_induced_moments); // Compute crystal Madelung potential
  Vector SlowComputeMadelungPotentialForInduction(bool include_induced_moments); // older, slower version.
  void BuildEwaldInteractionMatrices(); // Builds Tab-type matrices for PBC Madelung potential.
  int EwaldLookupAtomPair(int imonA, int iatomA, int imonB, int iatomB); // Look up index for Tab-type matrices for Madelung routines



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



  // GJB exploratory induction damping code
  void ReadAtomicInductionDampingFactors(ifstream& infile);

  // Old empirical Tkatchenko-Scheffler dispersion functions
  void ReadDispersionAtomTypes(ifstream& infile);
  void ReadDispersionCoefficients(ifstream& infile);
  double ComputeEmpiricalTwoBodyDispersion(); 
  // 3-body dispersion in this model is done via an altered call to
  // ComputeThreeBodyDispersion, using type = "Tkatchenko"


};

#endif
