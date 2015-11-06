#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <unistd.h> //for getcwd
// following, plus iostream, were all I had before
#include <string>
using std::string;
#include "params.h"
#include "constants.h"
using namespace hmbi_constants;

// Constructor sets default parameters
Params::Params() {

  Nproc = 1; // Default is serial (1 processor)

  IPRINT = 0;
  Warnings = 0;


  // Set path prefixes from the current path
  int MAX_PATH_LEN = 500;
  char cwd[MAX_PATH_LEN];
  getcwd(cwd, MAX_PATH_LEN); 
  base_path = cwd;
  base_path += "/";
  path_QM = base_path + "/qchem";
  path_MM = base_path + "/aiff";

  Analyze_Only = false;
  CreateJobsOnly_ = false;
  Use_Embedding_Charges = false;
  
  jobtype = 1; // single point energy
  jobtype_string = "energy";

  MM_Type = 2; // default is AIFF
  Periodic = false; // not periodic
  read_lattice_vectors = false; // default to angles & axis lengths
  read_fractional_coordinates = false; // default to Cartesian coords (Angstroms)
  crystal_symmetry = true; // default to use symmetry;
  tin_foil_boundary_conds = true; // default to use tin-foil boundary conditions

  Neglect_Many_Body = false; // do MB terms by default
  Counterpoise = false; // no CP correction
  DoQMBenchmark_ = false; // no full cluster QM benchmark by default
  
  // Q-Chem $rem section - must be read from input file
  // ought to have flag to check if it has been set
  QC_rem = "";
  Tinker_rem = "";
  HMBI_rem = "";
  AIFF_rem = "";  // by Ali
  QC_rem2 = "";
  QC_basis = ""; // by shuhao

  // Local truncation parameter defaults - fairly conservative
  do_local_two_body_truncation = false;  // no local truncation
  cutoff1 = 9.0; 
  cutoff0 = cutoff1 + 1.0;
  c0_set_explicitly = false;
  // The next two are used in toy routines to scan over possible cutoffs
  scan_local_cutoffs = false;
  min_cutoff1 = 2.0;
  N_cutoff_steps = 20;


  // AIFF parameters  --- by Ali
  LocalCoordinatesInitialized_ = false;
  AIFFBasisSet = "Sadlej";  // default
  IonizationPotential = 0.4638; // default is for water in au
  DampingFactor = 1.45; // default
  DoAtomicInductionDamping_ = false;
  DoDispersionDamping_ = true;
  MaxPolarizationCycles = 100; //default
  OrientDebug_ = false;
  UseGlobalCoords_ = false;
  BuildForceFieldOnly_ = false;
  RunJobsOnly_ = false;

  for (int i=0;i<93;i++)
    AtomicInductionDampingFactor[i]=DampingFactor;



  // Read default CamCasp version from Unix environment
  char *camcasp;
  camcasp = getenv("CAMCASP");
  if (camcasp == NULL) {
    printf("Parameters(): Error: Shell environment variable CAMCASP must be set.\n");
    exit(1);
  }
  CamCaspHome = camcasp;
  printf("CamCaspHome = %s\n",CamCaspHome.c_str());

  //CamCaspHome = "/home/software/camcasp-5.6";  // default

  // AIFF cutoffs & convergence parameters
  polarization_cutoff = 15.0; // default (angstroms)
  induction_convergence = 5; // default (10^-5 kJ/mol)
  induction_gradconvergence = 5; // default (10^-5 kJ/mol)
  induction_iter_scaling = 1.0; // default = no scaling (1.0).  
  preconverge_unit_cell_induced_moments = false; // by default, don't bother
  two_body_dispersion_cutoff = 20.0; // default (angstroms)
  three_body_dispersion_cutoff = 10.0; // default (angstroms)

  // This parameter is stored in hartrees, but input as kJ/mol:
  three_body_dispersion_print_threshold = 100000.0; // huge, no
						    // printing by
						    // default

  // Turn on/off different AIFF parts.  All true by default.
  do_aiff_electrostatics = true;
  do_aiff_induction = true;
  do_aiff_2body_dispersion = true;
  do_aiff_3body_dispersion = true;


  // Ewald parameters: 

  // We now have an automatic algorithm that determines nearly optimal
  // Ewald summation cutoffs and kappa.  You just need to specify
  // Ewald_accuracy, and the code does the rest.  Bigger values of
  // Ewald accuracy => better results.  Recommended values:
  // 10.0: decent, probably shouldn't go smaller if you want kJ/mol accuracy
  // 15.0: Usually sufficiently converged (default)
  // 20.0: Very good.
  // Note: we have found one case, acetamide, where 25 was needed to converge.
  // This might have to do with its really large unit cell.
  Ewald_accuracy = 15; 

  // other Ewald parameters, usually auto-determined from Ewald_accuracy.
  Ewald_kappa = -1.0; // With a large-enough Ewald_accuracy parameter,
		      // the resulting energy should be relatively
		      // invariant to this parameter, though the
		      // computational efficiency can vary a lot.
		      // Just let the code determine this parameter
		      // unless you have a good reason not to.

  Ewald_induction_type = 1; // Use finite cluster by default.

  // Ewald cutoffs.. Again, don't play with these unless you have a
  // good reason to.
  Recip_cutoffx = -1;
  Recip_cutoffy = -1;
  Recip_cutoffz = -1;
  Direc_cutoffx = -1;
  Direc_cutoffy = -1;
  Direc_cutoffz = -1;
  
  // Type of polarizabilities used to compute dispersion coefficients
  use_diagonal_anisotropic_pols = true;



  // Geometry optimization
  Do_Forces = false; // no forces
  UseFiniteDifferenceGradients_ = false; // Analytical gradients
  DoFiniteDifferenceFreqs_ = false; // FD hessian only right now
  MaxOptCycles = 100; // default number of geometry opt cycles
  OptCycle_ = 0; // cycle counter
  // The following 3 optimizer options seem to work well:
  CoordType_ = 2; // totally-connected highly deloc internal coords
  OptType_ = 3; // L-BFGS
  LineSearchType_ = 1; // Trust radius based on energy
  // Use DL-FIND optimization library for geom opts, etc?  Should be
  // TRUE unless you have a good reason not to (e.g. simplified
  // debugging)
  UseDLFind_ = true;

  // By default, compute the full hessian. Set to a monomer index
  // to do it for a given monomer.
  monomer_frequency_ = 0;

  //do_raman_intensities_ = 0;

  EstimateThreeBodyDispersion_ = false;
  MP2DispersionCorrection_ = false;

  tinker_debug_ = false;

}

Params::~Params() {

}

void Params::SetParameter(string param, string value) {
  if (IPRINT > 1)
    printf("Setting parameter <%s> to '%s'\n",param.c_str(),value.c_str());

  param = StringToUpper(param);

  // Job type
  if (param=="JOBTYPE") {
    value = StringToUpper(value);
    if (value=="SP" || value=="SINGLEPOINT" || value=="ENERGY") 
      {jobtype = 1; jobtype_string = "energy";}
    else if (value=="FORCE" || value =="OPT") 
      {jobtype = 2; Do_Forces = true; jobtype_string = "force";
	if (value=="FORCE") {
	  MaxOptCycles = 0;
	  UseDLFind_ = false;
	}
      }
    else if (value=="HESSIAN" || value=="FREQ" || value=="FREQUENCY") 
      {jobtype=3; jobtype_string = "hessian"; DoFiniteDifferenceFreqs_ = true;}
    else if (value=="EXPAND")
      {jobtype=101; jobtype_string = "expand";}
    else 
      {printf("Unknown Jobtype: %s\n",value.c_str()); exit(1);}
  }

  else if (param=="ANALYZE_ONLY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Analyze_Only = StringToInt(value); 
  }
  else if (param=="CREATE_JOBS_ONLY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    CreateJobsOnly_ = StringToInt(value); 
  }
  else if (param=="EMBEDDING_CHARGES") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Use_Embedding_Charges = StringToInt(value); 

  }
  else if (param=="NEGLECT_MANY_BODY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Neglect_Many_Body = StringToInt(value); 
  }

  else if (param=="COUNTERPOISE") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Counterpoise = StringToInt(value); 
  }

  else if (param=="QM_BENCHMARK") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    DoQMBenchmark_ = StringToInt(value); 
  }

  // Periodic boundary conditions
  else if (param=="PERIODIC") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Periodic = StringToInt(value);
  }

  else if (param=="READ_LATTICE_VECTORS") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    read_lattice_vectors = StringToInt(value);
  }

  else if (param=="READ_FRACTIONAL_COORDINATES") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    read_fractional_coordinates = StringToInt(value);
  }

  // Crystal symmetry
  else if (param=="SYMMETRY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    crystal_symmetry = StringToInt(value);
  }
  // Ewald boundary conditions for PBC calcs
  else if (param=="EWALD_TINFOIL_BOUNDARY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    tin_foil_boundary_conds = StringToInt(value);
  }

  // Paths
  else if (param=="QM_PATH") path_QM = base_path + value;
  else if (param=="MM_PATH") path_MM = base_path + value;

  // Print Level
  else if (param=="IPRINT") IPRINT = StringToInt(value);

  // MM type
  else if (param=="MM_CODE") {
    value = StringToUpper(value);
    if (value=="TINKER") MM_Type = 1;
    else if (value=="ORIENT" || value=="AIFF") MM_Type = 2;
    else if (value=="QCHEM") MM_Type = 3;
    else if (value=="EE-PA") MM_Type = 4;
    else {printf("Unknown MM_Type: %s\n",value.c_str()); exit(1);}
  }

  // Geometry optimization
  else if (param=="MAX_OPT_CYCLES") MaxOptCycles = StringToInt(value);
  else if (param=="OPTIMIZER") {
    value = StringToUpper(value);
    if (value=="DLFIND") {
      UseDLFind_ = true;
    }
    else if (value=="INTERNAL") {
      UseDLFind_ = false;

    }
    else {
      printf("Error: Unknown Optimizer: %s\n",value.c_str()); exit(1);
    }
  }
  else if (param=="OPT_TYPE") OptType_ = StringToInt(value);
  else if (param=="LINE_SEARCH_TYPE") LineSearchType_ = StringToInt(value);
  else if (param=="COORD_TYPE") CoordType_ = StringToInt(value);
  else if (param=="FDIFF_GRAD") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    UseFiniteDifferenceGradients_ = StringToInt(value); 
  }
  // QChem $rem section
  else if (param=="QC_REM") QC_rem = value;
  // QChem $basis section
  else if (param=="QC_BASIS") QC_basis = value;
  // Tinker $rem section - for keyfile
  else if (param=="TINKER_REM") Tinker_rem = value;
  // Orient $rem section --- by Ali
  else if (param=="AIFF_REM") AIFF_rem = value;
  // Secondary QChem $rem section
  else if (param=="QC_REM2") QC_rem2 = value;
  // Full parameter string - saved, just in case we want it later
  else if (param=="HMBI_PARAMS") HMBI_rem = value; 
    
  // Local truncation parameters
  else if (param=="LOCAL_2_BODY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
      do_local_two_body_truncation = StringToInt(value);
  }
  else if (param=="CUTOFF1") {
    cutoff1 = StringToDouble(value); 
    if (!c0_set_explicitly)
      cutoff0 = cutoff1 + 1.0; // c0 = c1 + 1.0 by default
  }
  else if (param=="CUTOFF0") {
    cutoff0 = StringToDouble(value); 
    c0_set_explicitly = true;
  }
  // and for scanning across possible local cutoff values:
  else if (param=="SCAN_LOCAL_CUTOFFS") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
      scan_local_cutoffs = StringToInt(value);
  }
  else if (param=="MIN_CUTOFF1") min_cutoff1 = StringToDouble(value); 
  else if (param=="N_CUTOFF_STEPS") N_cutoff_steps = StringToInt(value);



  else if (param=="EXPANSION_FACTOR") ExpansionFactor = StringToDouble(value); 

  // AIFF parameters
  else if (param=="DAMPINGFACTOR") {
    DampingFactor = StringToDouble(value); 
  }
  else if (param=="POLARIZATION_CUTOFF") {
    polarization_cutoff = StringToDouble(value); 
  }
  else if (param=="MAX_POLARIZATION_CYCLES") {
    MaxPolarizationCycles = StringToInt(value);
  }
  else if (param=="INDUCTION_CONVERGENCE") {
    induction_convergence = StringToDouble(value);
  }
  else if (param=="INDUCTION_GRADCONVERGENCE") {
    induction_gradconvergence = StringToDouble(value);
  }
  else if (param=="INDUCTION_ITER_SCALING") {
    induction_iter_scaling = StringToDouble(value);
    if (induction_iter_scaling > 1.0 || induction_iter_scaling <= 0.0) {
      printf("ERROR: Parameter INDUCTION_ITER_SCALING should be greater than zero and less than or equalt to 1.0.\n");
      exit(1);
    }
  }
    else if (param=="ATOMIC_INDUCTION_DAMPING") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    DoAtomicInductionDamping_ = StringToInt(value);
  }
    else if (param=="DISPERSION_DAMPING") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    DoDispersionDamping_ = StringToInt(value);
  }

  else if (param=="PRECONVERGE_UNIT_CELL_MOMENTS") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    preconverge_unit_cell_induced_moments = StringToInt(value);
  }
  else if (param=="EWALD_KAPPA" || param=="CONVERG_KAPPA") { // second-case for backward compatibility
    Ewald_kappa = StringToDouble(value);
  } 
  else if (param=="EWALD_ACCURACY") {
    Ewald_accuracy = StringToDouble(value);
  } 
  else if (param=="EWALD_INDUCTION_TYPE") {
    Ewald_induction_type = StringToInt(value);
  } 
  else if (param=="RECIP_CUTOFFX") {
    Recip_cutoffx = StringToInt(value);
  } 
  else if (param=="RECIP_CUTOFFY") {
    Recip_cutoffy = StringToInt(value);
  }
  else if (param=="RECIP_CUTOFFZ") {
    Recip_cutoffz = StringToInt(value);
  }
  else if (param=="DIREC_CUTOFFX") {
    Direc_cutoffx = StringToInt(value);
  }
  else if (param=="DIREC_CUTOFFY") {
    Direc_cutoffy = StringToInt(value);
  }
  else if (param=="DIREC_CUTOFFZ") {
    Direc_cutoffz = StringToInt(value);
  } 
  else if (param=="TWO_BODY_DISPERSION_CUTOFF") {
    two_body_dispersion_cutoff = StringToDouble(value); 
  }
  else if (param=="THREE_BODY_DISPERSION_CUTOFF") {
    three_body_dispersion_cutoff = StringToDouble(value); 
  }
  else if (param=="THREE_BODY_DISPERSION_PRINT_THRESHOLD") {
    // Input in kJ/mol, but store in hartrees
    three_body_dispersion_print_threshold = StringToDouble(value)/HartreesToKJpermole; 
  }
  else if (param=="DO_AIFF_ELECTROSTATICS") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    do_aiff_electrostatics = StringToInt(value);
    // if no electrostatics, turn off induction too.
    if (!do_aiff_electrostatics)
      do_aiff_induction = false;
  }
  else if (param=="DO_AIFF_INDUCTION") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    do_aiff_induction = StringToInt(value);
    if (do_aiff_induction && !do_aiff_electrostatics) {
      printf("Error: Cannot do AIFF induction without AIFF electrostatics.\n");
      printf("       Check your settings for DO_AIFF_ELECTROSTATICS and\n");
      printf("       DO_AIFF_INDUCTION.\n");
      exit(1);
    }
  }
  else if (param=="DO_AIFF_2BODY_DISPERSION") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    do_aiff_2body_dispersion = StringToInt(value);
  }
  else if (param=="DO_AIFF_3BODY_DISPERSION") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    do_aiff_3body_dispersion = StringToInt(value);
  }
  else if (param=="USE_DIAGONAL_ANISOTROPIC_POLARIZABILITIES") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    use_diagonal_anisotropic_pols = StringToInt(value);
  }

  else if (param=="AIFF_BASIS" || param=="ORIENTBASISSET") {
    AIFFBasisSet = StringToUpper(value); 
  }
  else if (param=="CAMCASPHOME") {
    printf("WARNING: Keyword CAMCASPHOME is no longer used.  Set the CAMCASP\n          environment variable instead.\n");
    //CamCaspHome = value; 

  }
  else if (param=="USE_GLOBAL_COORDINATES") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    UseGlobalCoords_ = StringToInt(value);
  }
  else if (param=="BUILD_FORCE_FIELD_ONLY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    BuildForceFieldOnly_ = StringToInt(value);
  }
  else if (param=="RUN_JOBS_ONLY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    RunJobsOnly_ = StringToInt(value);
  }
  else if (param=="INCLUDE_3BODY_DISPERSION") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    EstimateThreeBodyDispersion_ = StringToInt(value);
  }
  else if (param=="CORRECT_MP2_DISPERSION") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    MP2DispersionCorrection_ = StringToInt(value);
  }
  else if (param=="MONOMER_FREQUENCY") {
    monomer_frequency_ = StringToInt(value);
  }
  /*
  else if (param=="DO_RAMAN") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    do_raman_intensities_ = StringToInt(value);
  }
  */

  // Debugging
  else if (param=="TINKER_DEBUG") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    tinker_debug_ = StringToInt(value);
  }

  else {
    printf("Error: Params::SetParameter(): Unknown parameter %s\n",
	   param.c_str());
    exit(1);
  }
  
}

void Params::Print() {

  printf("\n************************\n");
  printf("   HMBI Parameters\n");
  printf("************************\n");
  printf("  Number of Processors = %d\n",Nproc);
  printf("  Main Process ID: %d\n",PID_);
  printf("  Job type = %d (%s)\n",jobtype,jobtype_string.c_str());
  printf("  Analyze_Only = %d\n",Analyze_Only);
  printf("  Do_Forces = %d\n",Do_Forces);
  printf("  Counterpoise = %d\n",Counterpoise);
  printf("  Periodic = %d\n",Periodic);
  printf("  Read_Lattice_Vectors = %d\n",read_lattice_vectors);
  printf("  Neglect_Many_Body = %d\n",Neglect_Many_Body);
  printf("\n");
  printf("  path_QM = %s\n",path_QM.c_str());
  printf("  path_MM = %s\n",path_MM.c_str());
  printf("\n");
  printf("  Print_Level = %d\n",IPRINT);
  printf("  MM_type = %d\n",MM_Type);
  printf("\n");
  printf("  Local_2_body = %d\n",do_local_two_body_truncation);
  printf("  cutoff1 = %.3f   cutoff0 = %.3f\n",cutoff1,cutoff0);
  printf("\n");

  printf("  Q-Chem $rem section:\n");
  printf("%s\n",QC_rem.c_str());
  printf("\n");

  printf("  Q-Chem $basis section:\n"); // by shuhao
  printf("%s\n",QC_basis.c_str());
  printf("\n");

  if (MM_Type==1) {
    printf("  Tinker $rem section:\n");
    printf("%s\n",Tinker_rem.c_str());
    printf("\n");
  }
  else if (MM_Type==2) {
    printf("  AIFF $rem section:\n");  // by Ali
    printf("%s\n",AIFF_rem.c_str());
    printf("\n");
  }
  else if (MM_Type==3) {
  printf("  Secondary QChem $rem section:\n");  // by Ali
  printf("%s\n",QC_rem2.c_str());
  printf("\n");
  }


  printf("************************\n\n");

}


// Some string functions, ought to be in a separate string class

// Convert C++ string to integer. 
int Params::StringToInt(string str) {

  int i;
  std::istringstream ss( str );
  ss >> i;
  /*  Not sure if this works:
  if (! ss.good()) {
    printf("Error: Params::StringToInt: conversion failed!\n");
    exit(1);
  }
  */
  return i;
}

// Convert C++ string to double.
double Params::StringToDouble(string str) {

  double value;
  std::istringstream ss( str );
  ss.clear();
  ss >> value;

  return value;
}

string Params::StringToUpper(string str) {
  for (int j=0; j<str.length(); ++j) 
    str[j]=toupper(str[j]);
  return str;
}
