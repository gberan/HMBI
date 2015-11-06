#ifndef _atom_h
#define _atom_h

#include <iostream>
#include <math.h>
#include <algorithm>
using namespace std;
#include "multipole.h"
#include "polarizability.h"

/* 
   A class for individual atoms.  

   GJB 12/08
*/

class Atom {

  string AtSym; // Atomic symbol
  int AtNum; // Atomic number
  double mass; // Atomic mass

  string DispersionAtomType_; // atom type for empirical disperion corrections

  // Dispersion coefficients - Note: we are phasing these out in favor
  // of ones calculated on the fly from the frequency-dependent
  // polarizabilities.  Do not use these.
  double C6, C8, C10; // C8/C10 are generally zero for H atoms
  // For correcting MP2 dispersion
  double C6_UCHF, C8_UCHF, C10_UCHF; // UCHF C6 coefficent
  // note: haven't yet implemented code for C8_UCHF or C10_UCHF

  // (Uncoupled) isotropic frequency-dependent polarizabilities
  Vector freq_pol_dipole; // dipole-dipole
  Vector freq_pol_quad; // quadrupole-quadrupole
  

  Vector xyz; // Cartesian Position in global coordinates, Angstroms
  Vector local_xyz; // Position in local coordinates, Angstroms

  double point_charge; // point charge, used for charge embedding

  // For specifying Tinker parameters:
  int atom_index; // counter inside the given fragment
  int MM_atom_type; 
  int Nconnected; // number of connected atoms
  int *connectivity; // list of connected atoms numbers

  Multipole MultipoleMoments; // Stores the permanent multipole moments
  // Induced moments are stored in the Dimer or Cluster, since they depend on
  // the environment.  They are not purely an Atomic property.  Exception,
  // we store the induced multipole moments in a periodic system here:
  Multipole InduceMultipoleMoments; // by shuhao; the converged
				    // induced multipole moments for
				    // the central unit cell atoms 

  // Store the isotropic_dipole_polarizability for this atom. This
  // gets used in combination rules for estimated C6 coefficients for
  // a pair of atoms based on the individual atomic C6 coefficients.
  // It also can be used to estimate C9 3-body dispersion coefficients
  // from the C6 ones.  Obsolete now.  Use freq_pol_* instead.
  double isotropic_dipole_polarizability;
  
  // Store the anisotropic distributed polarizabilities in spherical
  // tensor form.  These are used for calculating the induction
  // energy.
  Polarizability Pols; 

  // Logical flags to say if the Multipoles/Polarizabilities have been
  // initialized
  bool init_multipoles;
  bool init_InduceMultipoles; // by shuhao; for periodic systems
  bool init_pols;


 public:
  // Constructors & Destructor
  Atom(); // The real work is done by Initialize
  Atom(const Atom &other); // copy constructor
  ~Atom();
  // The function to populate the objects:
  void Initialize(int index, string symbol, double x, double y, double z, 
		  double charge=0.0); //q-chem style xyz
  void Initialize(int index, string symbol, double coord[3], int type, int
		  Nconnect, int *connect, double charge = 0.0); // tinker-style xyz

  void SetPosition(double new_xyz[3]); // for moving atoms
  void SetPosition(Vector new_xyz);  // for moving atoms
  double GetPosition(int dim) {return xyz[dim];}; // get one of xyz (dim =0,1,2)
  Vector GetPosition() {return xyz;}; // return global coords vector
  
  void SetLocalPosition(double new_loc_xyz[3]); // for local positions by Ali
  Vector GetLocalPosition() {return local_xyz;}; // return local coords vector
  double GetLocalPosition(int dim) {return local_xyz[dim];}; // get one element of local_xyz (dim =0,1,2) by Ali
 
  // For working with isotropic frequency-dependent polarizabilities - by shuhao
  // dipole-dipole polarizabilities:
  void SetFreq_Pol_Dipole(Vector new_freq_pol_dipole);
  Vector GetFreq_Pol_Dipole() {return freq_pol_dipole;};
  // quadrupole-quadrupole polarizabilities:
  void SetFreq_Pol_Quad(Vector new_freq_pol_quad);
  Vector GetFreq_Pol_Quad() {return freq_pol_quad;};

  
  // Set multipole moments/polarizabilities
  void SetMultipoleMoments(const Multipole& moments) {MultipoleMoments = moments; init_multipoles = true;};
  Multipole GetMultipoleMoments() const {return MultipoleMoments;};
  void SetPolarizability(const Polarizability& other) {Pols = other; init_pols = true;};
  Polarizability GetPolarizability() const {return Pols;};
  void SetIsotropicDipolePolarizability(double value) {isotropic_dipole_polarizability = value;};
  double GetIsotropicDipolePolarizability() {return isotropic_dipole_polarizability;};


  void SetInduceMultipoleMoments(const Multipole& Inducemoments) {InduceMultipoleMoments = Inducemoments;  init_InduceMultipoles = true;};// by shuhao
  Multipole GetInduceMultipoleMoments() {return InduceMultipoleMoments;};// by shuhao
  


  bool MultipolesAreInitialized() {return init_multipoles;};  
  bool InduceMultipolesAreInitialized() {return init_InduceMultipoles;}; // by shuhao
  bool PolarizabilityIsInitialized() {return init_pols;};

  double GetAtomicInductionDampingFactor(); // GJB exploratory

  int GetAtomIndex() {return atom_index;}; // get atom index by Ali
  int GetAtomicNumber();
  void SetAtomicMass();
  double GetAtomicMass() {return mass;};
  string GetSymbol() {return AtSym;};
  int GetMMAtomType() {return MM_atom_type;};
  int GetNumberOfConnectedAtoms() {return Nconnected;};
  int GetConnectivity(int i) {return connectivity[i];};

  void SetDispersionAtomType(string type) {DispersionAtomType_ = type;};
  string GetDispersionAtomType() {return DispersionAtomType_;};

  void SetDispersionCoefficients(double C6_in, double C8_in=0.0, double C10_in=0.0) {C6 = C6_in; C8=C8_in; C10=C10_in;};
  void SetUCHFDispersionCoefficients(double C6_in, double C8_in=0.0, double C10_in=0.0) {C6_UCHF = C6_in; C8_UCHF = C8_in; C10_UCHF = C10_in;};
  double GetC6Coefficient( string type = "normal") {
    if (type=="UCHF") 
      return C6_UCHF; 
    else 
      return C6;
  };
  double GetC8Coefficient( string type = "normal") {
    if (type=="UCHF") 
      return C8_UCHF; 
    else 
      return C8;
  };
  double GetC10Coefficient( string type = "normal") {
    if (type=="UCHF") 
      return C10_UCHF; 
    else 
      return C10;
  };

  void PrintDispersionCoefficients() {printf("%s: C6 = %f, C8 = %f, C10 = %f\n",AtSym.c_str(), C6,C8,C10);};
  

  double GetCoordinate(int dim) { return xyz[dim]; };
  double GetInterAtomicDistance(const Atom& atomB);

  // Routines to evaluate classical interactions
  // Construct Tab matrix
  Matrix BuildInteractionMatrix(Vector thisRotVec, double thisRotAng, 
				Atom& other, Vector otherRotVec,
				double otherRotAng, double beta_damp);
  Matrix BuildInteractionMatrixGradient(Matrix Tab, Vector thisRotVec, 
					double thisRotAng, 
					Atom& other, Vector otherRotVec,
					double otherRotAng, double beta_damp);
  // Routines to evaluate classical interactions in periodic systems
  // Construct reciprocal space Tab matrix for
  // Aatom-with-imageatoms----Batom-with-imageatoms // by shuhao
   Matrix BuildRecipInteractionMatrix(Vector thisRotVec, double thisRotAng,
				      Atom& other, Vector otherRotVec, 
				      double otherRotAng, 
				      int kx, int ky, int kz, double CellV,
				      Vector RecipCellx, Vector RecipCelly, 
				      Vector RecipCellz, double beta_damp); 
   
  // Construct direct space Tab matrix for
  // Aatom-with-imageatoms----Batom-with-imageatoms // by shuhao
   Matrix BuildDirecInteractionMatrix(Vector thisRotVec, double thisRotAng,
				      Atom& other, Vector otherRotVec, 
				      double otherRotAng, int nx, int ny, 
				      int nz, double CellV, Vector UnitCellx, 
				      Vector UnitCelly, Vector UnitCellz, 
				      double beta_damp); 
   // for the case of |kn| =0 in the reciprocal space; by shuhao
   Matrix OldBuildRecipInteractionMatrix_kn0(Vector thisRotVec, double thisRotAng,
					  Atom& other, Vector otherRotVec,
					  double otherRotAng, int kx, int ky, 
					  int kz, double CellV, 
					  Vector RecipCellx, Vector RecipCelly,
					  Vector RecipCellz, double beta_damp);

   // for the case of |kn| =0 in the reciprocal space; GJB version
   // that is slightly more efficient than original.  Could still be made
   // better by reducing the size of Tab.
   Matrix BuildRecipInteractionMatrix_kn0(Vector thisRotVec, double thisRotAng,
					  Atom& other, Vector otherRotVec,
					  double otherRotAng, double CellV, 
					  Vector RecipCellx, Vector RecipCelly,
					  Vector RecipCellz);
   Matrix BuildRecipInteractionMatrix_kn0_2(Vector thisRotVec, double thisRotAng,
					  Atom& other, Vector otherRotVec,
					  double otherRotAng, int kx, int ky, 
					  int kz, double CellV, 
					  Vector RecipCellx, Vector RecipCelly,
					  Vector RecipCellz, double beta_damp);

   // for the case of intramolecular correction in the direct space; by shuhao 
   Matrix BuildDirecInteractionMatrix_intramol(Vector thisRotVec, 
					       double thisRotAng, Atom& other, 
					       Vector otherRotVec,
					       double otherRotAng, int nx, 
					       int ny, int nz, double CellV,
					       Vector UnitCellx, 
					       Vector UnitCelly, 
					       Vector UnitCellz, 
					       double beta_damp);

   // Damping factor used in determining induction contribution
   double TangToenniesDampingFactor(int n, double beta, double R, 
				    int ideriv=0);
  


  // Printing commands
  void PrintQChemCartesian(FILE *outfile=stdout) {
    fprintf(outfile, "%-2s  %10.6f  %10.6f  %10.6f\n", 
	    AtSym.c_str(),xyz[0],xyz[1],xyz[2]);};
  void PrintQChemCartesian(FILE *outfile, string symbol) {
    fprintf(outfile, "%-2s  %10.6f  %10.6f  %10.6f\n", 
	    symbol.c_str(),xyz[0],xyz[1],xyz[2]);};
  void PrintQChemEmbeddingCharge(FILE *outfile=stdout) {
    fprintf(outfile, "%10.6f  %10.6f  %10.6f  %10.6f\n", xyz[0],xyz[1],
	    xyz[2], point_charge);};
  void PrintTinkerCartesian(int shift=0, FILE *outfile=stdout);
  void PrintTinkerEmbeddingCharge(int shift=0, FILE *outfile=stdout);

  // Operator overloading
  Atom& operator=(const Atom& other);


  /* Empirical dispersion coefficients */

  // These functions are used to estimate dispersion coefficients
  // based on Tkatchenko & Scheffler's tables.  The AIFF functions
  // below should be used instead.

  // Get the parameters: property is one of "Pol", "C6", "C9" or "Rvdw".
  double LookupAtomicDispersionParameter(string property);
  // Compute the C6 coefficient between this atom and another
  double EstimateC6Coefficient(Atom other);
  // Compute the C9 coefficient between this atom and 2 others
  double EstimateC9Coefficient(Atom other1, Atom other2, string type = "default");


  /* AIFF dispersion coefficients */

  // Obtained via Casimir-Polder integration over the
  // frequency-dependent polarizabilities.

  // Compute the C6, C8, or C10 coefficient between this atom and
  // another
  double CasimirC6Coefficient(Atom other);
  double CasimirC8Coefficient(Atom other);
  double CasimirC10Coefficient(Atom other);

  // Compute the C9 coefficient between this atom and 2 others
  double CasimirC9Coefficient(Atom other1, Atom other2);

  // GaussLegendre method for the Casimir-Polder integration
  Vector GaussLegendreWeights(double x1, double x2, int n);
  


};

#endif
