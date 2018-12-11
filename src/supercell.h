#ifndef _supercell_h
#define _supercell_h

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
#include "vector.h" // my vector class, "Vector"               
#include <vector> // C STL library vector, "vector"
#include "cluster.h"
using std::vector;  
#include <time.h>

/*
Class description: To handle supercell calculations
As of 7/18/12, it should not run any QM jobs
*/


class Supercell {

 private:
  
  Supercell();
  ~Supercell();

  int Ncells;//Number of cells
  int Natoms;
  int NMon;
  int Supercell_Natoms;
  int Supercell_NMon;
  int Natoms_lp;  // # atoms in unitcell minus the explicit lone pairs treated as atoms (in AMOEBA parameters)
  Vector Supercell_Size;
  Vector* Supercell_AtomicCoordinates;
  Vector* Supercell_FractionalCoordinates;

  double *AtomicMasses, *AtomicMasses_lp;


//  vector<vector<vector<vector<double>>>> Order_4_Tensor;
//  Order_4_Tensor Supercell_AtomicCoordinates;

  Vector* unit_cell; // stores primitive unit cell.  First dim loops over vectors,
                          // and second over xyz coordinates for each vector.

  Vector* recip_boundary_pts; // stores reciprocal space zone boundary points 
                              // which the user specifies in the $reciprocal_space_points section in the inputfile

  Vector* recip_grid_pts; // stores reciprocal space points based on the grid formed using the user-specified stepsize

  Monomer* Supercell_Monomers; // stores the monomers within the supercell

//  Dimer Dimers;

  Matrix SupercellHess_MM;
  Matrix SupercellHess_QM;
  Matrix hess;
  bool MM_SupercellHess_Init;
  bool QM_SupercellHess_Init;
  bool supercell_init;
  int primitive_cell_index;
  int tot_recip_boundary_pts;
  double recip_stepsize;
  int total_k_points;
  bool *plot_break;
  Vector phonopy_mapping;
  Vector* supercell_unit_cell;

  Matrix phonon_mat;

//  double C_one[3][3][3][Natoms][Natoms];
//  double C_two[3][3][3][3][Natoms][Natoms];

 public:
  // instantiate the single instance of this class
  static Supercell& supercell() {
    static Supercell thesupercell;
    return thesupercell;
  }                  

  void Initialize(ifstream &infile, int Nproc, string input_filename);
  void UpdateSupercellCoordinates();
  void CreateSupercellCoordinates();
  void CreateSupercellUnitCell();
  void SetSupercellInit(bool boolean) {supercell_init = boolean;};
  bool GetSupercellInit() {return supercell_init;};
  void CreateSupercellMonomerList();
  Vector CreatePhonopyMapping();
  Matrix ShiftOddMatrix(Matrix supercell);
  void CreateSupercellTinkerJob();
  void GetSuperCellLatticeParameters(double& a,double& b,double& c,
				     double& alpha,double& beta,double& gamma);
  void PrintTinkerCartesian(FILE *outfile);
  void RunHMBISupercellTinkerHessianJob();
  void SetSupercellHessian();
  Matrix ReadSupercellHessian(int type);
  void ComputeHMBISupercellHessian();
  int GetCellIndex(int na, int nb, int nc, int na1, int nb1, int nc1);
  Vector GetUnitCellVector(int index) {return supercell_unit_cell[index];};
  int GetTotalKPoint(){return total_k_points;};
  int GetSupercellNAtoms(){return Supercell_Natoms;};
  Matrix GetPhononMatrix(){return phonon_mat;};

  void ReadReciprocalSpaceBoundaryPointsAndCreateGrid(ifstream& infile) ;
  void CreateSupercellFractionalCoordinates();
  string GetDFTBGeomInput();
  void CreateReciprocalSpaceGrid();
  void ComputePhonons();
  void SetAtomicSymbols();

  void CreateReciprocalSpaceSamplingGrid(ifstream& infile);
  void ComputeThermalPropertiesUsingPhonons(ifstream& infile);
  void ComputeVibrationalEnergy();
  void UnMassWeightedEigenVectors(Matrix& Hessian);

  void Compute_C_one_matrix();
  void Compute_C_two_matrix();       
  void Compute_Elastic_Square_Brackets();
  void Compute_Elastic_Round_Brackets();

  void RunHMBISupercellTinkerFDHessianJob();
  Matrix ReadSupercellFDHessian(int type);

  void ComputePDOS_old(Matrix phonon_matrix);
  void ComputePDOS(Matrix phonon_matrix);       

  void PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index, int k_index);
  void PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index, int k_index, string file_name);
  
  void CreateSupercellInputFile();

  Vector GetSupercellCoordinates(int i){return Supercell_AtomicCoordinates[i];};  Matrix GetOnlyRealMatrixComponent(Matrix FullMatrix,Vector  eigvals);

};


#endif

