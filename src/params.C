#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <unistd.h> //for getcwd
// following, plus iostream, were all I had before
#include <string>
using std::string;
using namespace std;
#include "params.h"

// Constructor sets default parameters
Params::Params() {

  Nproc = 1; // Default is serial (1 processor)

  IPRINT = 0;
  Warnings = 0;

  input_filename = "";

  // Set path prefixes from the current path
  int MAX_PATH_LEN = 500;
  char cwd[MAX_PATH_LEN];
  getcwd(cwd, MAX_PATH_LEN); 
  base_path = cwd;
  base_path += "/";
  path_QM = base_path + "/qm";
  path_MM = base_path + "/mm";
  path_tmp_files = base_path + "/file_storage";
  path_qm_tmp_files = path_tmp_files + "/qm";
  path_mm_tmp_files = path_tmp_files + "/mm";

  path_hessian_files = base_path + "/hessian_files";
  path_hessian_qm = path_hessian_files + "/qm";
  path_hessian_mm = path_hessian_files + "/mm";

  Is_Opt_Complete = false;
  Are_Forces_Available = false;
  Are_Hessians_Available = false;
  Initialize_Hessians = false;
  Analyze_Only = false;
  CreateJobsOnly_ = false;
  Use_Embedding_Charges = false;
  Read_Cif_File = false; //Watit 
  
  jobtype = 1; // single point energy
  jobtype_string = "energy";
  NumberOfOptConstraints = 1; //One optimization constraints by default
  Find_Enthalpy = false;
  ExternalPressure = 0.000101325; //external pressure (user input in GPa but a.u. is hartrees/bohr^3) is 1 atm by default
  Temperature = 298 ; // in kelvin
  tolerance_max_g = 4.5e-5;
  tolerance_rms_g = 3.0e-5;
  tolerance_e = 100e-8;

  //Quasiharmonic calculations
  Do_QuasiHarmonic = false;//No quasi-harmonic calculates be default
  Are_QHF_Available = false;//Quasi-harmonic frequecuies are not availble by default
  Volume_Gruneisen = true;//volume based gruneisen parameters
  Save_QuasiHarmonic_Freq = true;// values the frequencies calculations are 
                                  // saved rather an being written over.
  QuasiStepsDone = 0;//If Quasiharmonic steps are partially done.
                      //If gruneisen is isotropic, there are 3 steps, otherwise there are 6 steps (not including steps for lattice angles)
  ReadQuasiGeometries = false;//If Plus or minus geometries in quasiarharmonic already existing

  QM_Type = 1; // default is Qchem
  MM_Type = 2; // default is AIFF
  Periodic = false; // not periodic
  read_lattice_vectors = false; // default to angles & axis lengths

  crystal_symmetry = true; // default is to use symmetry
  space_symmetry = false; // default to not use symmety for non-periodic systems and to use symmetry for periodic crystals;
  Symmetry_Tolerance = 4;// default tolerance to symmetry
  monomer_symmetry = false;//default for determining monomer symmetry
  mm_symmetry = true;//default for using symmetry for the AIFF terms
                     //AIFF only
  lattice_symmetry = false;//use symmetry for the lattice when space symmetry
                           //is false. If periodic is true, lattice symmetry                         
                          // will be set to true unless otherwise specified
  print_symmetry_info = false;
  symmetry_set = false;//bool allows for different default symmetry settings for periodic and non-periodic systems settings 
  old_dimer_symm_names_ = false;

  tin_foil_boundary_conds = true; // default to use tin-foil boundary conditions
  //freeze_unitcellparams = true;//true until the symmetry is working into the lattice
  freeze_unitcellparams = false;
  freeze_atom = false;
  freeze_params_until = -1;
  do_change_volume = false;
  change_volume = 0.0;
  change_a = 0.0;
  change_b = 0.0;
  change_c = 0.0;
  change_alpha = 0.0;
  change_beta = 0.0;
  change_gamma = 0.0;

  Neglect_Many_Body = false; // do MB terms by default
  MM_Only = false;// MM only
  QM_Only = false;// QM only
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
  MolPro_rem = "";//by yoni
  MolPro_inst = "";//by yoni
  MolPro_CBS = "";//by yoni
  MolPro_HF = "";//by yoni
  MolPro_CCSDT_Basis = "";//by yoni
  MolPro_CCSDT_MP2 = "";//by yoni
  MolPro_CCSDT_Inst = "";//by yoni
  PSI4_rem = ""; //CSG
  crystal_heading = "";//by yoni
  crystal_basis = "";//by yoni
  crystal_ending = "";//by yoni
  GaussianHeader_ = ""; //JDH
  DaltonSection_ = ""; //JDH
  OrcaHeader_ = ""; //JDH
  Molpro1_ = ""; //JDH
  Molpro2_ = ""; //JDH
  PSI4Header_ = ""; //CSG
  //GaussianHeader2_ = ""; // JDH
  QESupercellSection_ = "1 1 1"; // JLM
  QESpeciesSection_ = ""; // JLM
  dispersion_correction_ = true;
  dispersion_type_ = "XDM"; //JLM
  force_xc_set_ = false; //JLM
  force_xc_ = ""; //JLM
  qe_basis_ = 80; //JLM
  qe_basis_multiplier_ = 10; //JLM
  molecule_ = false; //JLM
  no_freq_disp_ = false; //JLM
  use_phonopy_freq_ = true; //JLM
  use_phonon_freq_ = false; //JLM
  use_dftb_freq_ = false; //JLM

  //////////////////////////////////////////////////////
  // JDH's functions for Magnetic property calculations:
  ////////////////////////////////////////////////////// 
  // set the default behavior of a bunch of NMR stuff...

  // Job control
  NMRJob_ = false; 
  EFGJob_ = false; //xxxxx
  Num_Asymmetric_Monomers_ = 1;//xxxxx
  TwoBodyCutoff_ = 4;//xxxxx
  ClusterCutoff_ = 0.0;//xxxxx
  UseScaledEFGTensors_ = true; 
  QM_Package_ = "QCHEM"; // Default to using qchem JLM
  
  // Electrostatic embedding
  UseElectrostaticEmbedding_ = false; //xxxxx
  ElectrostaticEmbeddingCutoff_ = 30;//xxxxx
  UseTwoBodyCharge_ = false;

  AdvancedEmbeddingScheme_ = "NONE";//xxxxx

  UseEwald_ = false;//xxxxx
  Ewald_Cutoff_ = 5.0; // JDH (DEBUG: CHECK THIS)//xxxxx
  Rcond_ = 0.000001; // Default values used throughout 2017 Charge embedding paper: 0.000001 //xxxxx
  
  UseSelfConsistentEmbedding_ = false;//xxxxx
  SelfConsistentMaxThreshold_ = 0.001;//xxxxx
  SelfConsistentRMSThreshold_ = 0.005;//xxxxx
  SeedSCE_ = false;

  
  charge_embedding_rank_ = 0;  // Default to only using charges (dalton code for higher-order multiploar embedding is unfinished)//xxxxx

  // Custom basis functions:
  Custom_Basis_ = 0;//xxxxx
  MixedBasisCutOff_ = 0; //xxxxx
  MixedBasisCutOff2_ = 0; //xxxxx
  MixedBasisCutOff3_ = 0; //xxxxx
  MixedBasisLevel1_ = ""; //xxxxx
  MixedBasisLevel2_ = ""; //xxxxx
  MixedBasisLevel3_ = ""; //xxxxx
  H_Basis_ = "";
  C_Basis_ = "";
  N_Basis_ = "";
  O_Basis_ = "";
  S_Basis_ = "";
  Cl_Basis_ = "";
  I_Basis_ = "";
  Sn_Basis_ = "";
  P_Basis_ = "";
  K_Basis_ = "";
  Na_Basis_ = "";
  Br_Basis_ = "";
  F_Basis_ = "";
  W_Basis_ = "";
  V_Basis_ = "";
  
  
  
  // QM Package interfaces:


  // CBS extrapolation
  CBS = false;
  CBS_basis1 = 3;//tz or aTZ by default
  CBS_basis2 = 4;//qz or aQZ by default

  //CCSD(T) correction
  CCSDT_correction = false;

  
  single_file_monomer_hess = false;

  memory = 200;//currently only implimated for molpro

  // Local truncation parameter defaults - fairly conservative
  do_local_two_body_truncation = false;  // no local truncation
  cutoff1 = 9.0; 
  cutoff0 = cutoff1 + 1.0;
  c0_set_explicitly = false;
  // The next two are used in toy routines to scan over possible cutoffs
  scan_local_cutoffs = false;
  min_cutoff1 = 2.0;
  N_cutoff_steps = 20;

  find_spatial_damping_gradient = true;


  // AIFF parameters  --- by Ali
  LocalCoordinatesInitialized_ = false;
  AIFFBasisSet = "Sadlej";  // default
  IonizationPotential = 0.4638; // default is for water in au
  DampingFactor = 1.5; // default
  MaxPolarizationCycles = 100; //default
  OrientDebug_ = false;
  UseGlobalCoords_ = false;
  BuildForceFieldOnly_ = false;
  RunJobsOnly_ = false;
  NumberMonomerTypes = 1; // default only 1 type of monomer
  montype = new string[NumberMonomerTypes];

  // Read default CamCasp version from Unix environment
  /*char *camcasp;
  camcasp = getenv("CAMCASP");
  if (camcasp == NULL) {
    printf("Parameters(): Error: Shell environment variable CAMCASP must be set.\n");
    exit(1);
  }
  CamCaspHome = camcasp;
  printf("CamCaspHome = %s\n",CamCaspHome.c_str());*/

  // AIFF cutoffs & convergence parameters
  polarization_cutoff = 15.0; // default (angstroms)
  induction_convergence = 5; // default (10^-5 kJ/mol)
  induction_gradconvergence = 5; // default (10^-5 kJ/mol)
  induction_iter_scaling = 1.0; // default = no scaling (1.0).  
  preconverge_unit_cell_induced_moments = false; // by default, don't bother
  two_body_dispersion_cutoff = 20.0; // default (angstroms)
  three_body_dispersion_cutoff = 10.0; // default (angstroms)

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
  UseFiniteDifferenceGradients_ = false; // FD gradient
  DoFiniteDifferenceFreqs_ = false; // FD hessian 
  DoEnergyFiniteDifferenceFreq_ = false; //Finite difference hessians using energies instead of gradients
  couple_gradients_ = false;

  //if you don't have Dr. Berans correction to the  hessian for
  //polarization in tinker set Do_fdTinkerHessian_ to true
  Do_fdTinkerHessian_ = false;// FD hessian for tinker on full system 

  LatticeVectorFiniteDifference = false;// finite difference on lattice vectors
  do_freq_after_opt = false;
  MaxOptCycles = 150; // default number of geometry opt cycles
  OptCycle_ = 0; // cycle counter
  // The following 3 optimizer options seem to work well:
  CoordType_ = 2; // totally-connected highly deloc internal coords
  OptType_ = 3; // L-BFGS
  LineSearchType_ = 1; // Trust radius based on energy

  use_approx_hessian = false; //by default we use identity matrix but using approx hessian at lower basis might improve convergence
  

  // Use DL-FIND optimization library for geom opts, etc?  Should be
  // TRUE unless you have a good reason not to (e.g. simplified
  // debugging)
  UseDLFind_ = true;
  UseKNitro_ = false;
  UseConjugate_ = false;
  UseSteepestDescent_ = false;
  UseLBFGS_ = false;
  use_local_optimizer_ = false;
  step_size_init_ = 0.65; //JLM

  // By default, compute the full hessian. Set to a monomer index
  // to do it for a given monomer.
  monomer_frequency_ = 0;


  EstimateThreeBodyDispersion_ = false;
  MP2DispersionCorrection_ = false;
  Convergence = 1;

  //Frequecy calc
  Do_Freq = false; //no HMBI freq (hessian) calc 
  Do_Raman = false; //Watit
  tinker_debug_ = false;

  // supercell params
  supercell_job = false;
  supercell_job_exist = false;
  supercell_size.Initialize(3);
  supercell_size[0] = 3;
  supercell_size[1] = 3;
  supercell_size[2] = 3;

  IsThermalUsingPhonons = true;
  ThermalForMultiTemp = false;
  create_supercell_input = false;
  TotalKPoint = 1;

  deuterated = false;

  stupid_var_ = 1.0;
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
      {jobtype = 1; jobtype_string = "energy";Do_Freq = false; Do_Forces = false;}
    else if (value=="FORCE" || value =="OPT") 
      {jobtype = 2; Do_Forces = true; jobtype_string = "force";
	Do_Freq = false;
	if (value=="FORCE") {
	  //MaxOptCycles = 0;
	  UseDLFind_ = false;
	  UseKNitro_ = false;
          UseConjugate_ = false;
          UseSteepestDescent_ = false;
          UseLBFGS_ = false;
	}
      }
    else if (value=="HESSIAN" || value=="FREQ" || value=="FREQUENCY") 
      {jobtype=3; Do_Freq = true; Do_Forces = false;jobtype_string = "hessian";}
    else if (value=="NMR")
      {jobtype=5; jobtype_string = "nmr";} // JDH
    else if (value=="EXPAND")
      {jobtype=101; jobtype_string = "expand";}
    else if (value=="PHONON")
      {jobtype=4; jobtype_string = "phonon";}
    else if (value=="RAMAN") // Watit
      {jobtype=6; Do_Raman = true; Do_Freq = true; Do_Forces = false;jobtype_string = "raman";}
    else 
      {printf("Unknown Jobtype: %s\n",value.c_str()); exit(1);}
  }
  
  else if (param=="QUASIHARMONIC"){
    value = StringToUpper(value);
    if (value == "FALSE" || value == "0")
      Do_QuasiHarmonic = false;
    else{
      Do_QuasiHarmonic = true;
      do_freq_after_opt = true;
    }
  }
  else if(param== "VOLUME_GRUNEISEN" || param=="VOLUME_GRUN"){
    value = StringToUpper(value);
    if(value == "FALSE" || value == "0")
      value = "0";
    else
     value = "1";
    Volume_Gruneisen = StringToInt(value);
  }
  else if(param== "SAVE_QHF"){
    value = StringToUpper(value);
    if(value == "FALSE" || value == "0")
      value = "0";
    else
      value = "1";
    Save_QuasiHarmonic_Freq = StringToInt(value);
  }
  else if(param=="NUMBER_OF_QUASI_STEPS" ||param== "NUMBER_OF_QUASIHARMONIC_STEPS"){
    QuasiStepsDone =  StringToInt(value);
  }
  else if(param=="READ_QUASI_GEOMETRIES" || param=="READ_QUASIHARMONIC_GEOMETRIES" ){
    value = StringToUpper(value);
    if(value == "FALSE" || value == "0")
      value = "0";
    else
     value = "1";
   ReadQuasiGeometries = StringToInt(value);
  } 
  else if(param=="FIND_ENTHALPY"){
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Find_Enthalpy = StringToInt(value);  
  }
  else if (param=="NUMBER_CONSTRAINTS") {
    NumberOfOptConstraints = StringToInt(value);
  }
  else if (param=="CONSTANT_VOLUME") {
    constant_volume = StringToDouble(value);
  }
  
  else if (param=="EXTERNAL_PRESSURE" || param=="PRESSURE") {               
    ExternalPressure = StringToDouble(value);
    Find_Enthalpy = true;
    
  }
  else if(param=="TEMPERATURE"){
    Temperature = StringToDouble(value);
  }
  else if(param=="TOLERANCE_MAX_G"){
    tolerance_max_g = StringToDouble(value);
  }
  else if(param=="TOLERANCE_RMS_G"){
    tolerance_rms_g = StringToDouble(value);
  }
  else if(param=="TOLERANCE_E"){
    tolerance_e = StringToDouble(value);
  }
  else if (param=="USE_APPROX_HESSIAN") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    use_approx_hessian = StringToInt(value);
  }
  else if (param=="IS_OPT_COMPLETE") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";               
    Is_Opt_Complete = StringToInt(value);
  }
  else if (param=="ARE_FORCES_AVAILABLE") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";  
    Are_Forces_Available = StringToInt(value);
  }
  else if (param=="ARE_HESSIANS_AVAILABLE") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";           
    Are_Hessians_Available = StringToInt(value);
  }
  else if(param=="ARE_QHF_AVAILABLE" || param=="ARE_QHA_AVAILABLE" ){
    value = StringToUpper(value);
    if(value=="FALSE" ||value=="0") value = "0";
   else value = "1";
    Are_QHF_Available = StringToInt(value); 
  }
  else if (param=="INITIALIZE_HESSIANS") {               
    value = StringToUpper(value);        
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Initialize_Hessians = StringToInt(value);               
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
  else if (param=="READ_CIF_FILE") {  //Watit
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Read_Cif_File = StringToInt(value);
  }
  else if (param=="CONVERGENCE") {  //Watit
    value = StringToUpper(value);
    if (value=="LOOSE" || value=="0") value = "0";
    else if (value=="TIGHT" || value=="2") value = "2";
    else if (value=="VTIGHT" || value=="VERYTIGHT" || value=="3") value = "3";
    else value = "1";
    Convergence = StringToInt(value);
  }
  else if (param=="NEGLECT_MANY_BODY") {
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    Neglect_Many_Body = StringToInt(value); 
  }
  else if (param=="FULL_MM_ONLY"){
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    MM_Only = StringToInt(value);  
  }
  else if (param=="FULL_QM_ONLY"){
    value = StringToUpper(value);
    if (value=="FALSE" || value=="0") value = "0";
    else value = "1";
    QM_Only = StringToInt(value);  
    couple_gradients_ = StringToInt(value);
  }
  // CBS exptrapolation
  else if(param=="CBS"){
    value = StringToUpper(value);
    if(value=="FALSE" || value=="0") value = "0";
    else value = "1";
    CBS = StringToInt(value);
  }
  else if(param=="CBS_BASIS1"){
   value = StringToUpper(value);
   if(value == "2" ||value == "DZ" || value == "ADZ" || value == "CC-PVDZ" || value == "AUG-CC-PVDZ")
     value = "2";
   else if(value == "3" || value == "TZ" || value == "ATZ" || value == "CC-PVDZ" || value == "AUG-CC-PVTZ")
     value = "3";
   else if(value == "4" || value == "QZ" || value == "AQZ" || value == "CC-PVQZ" || value == "AUG-CC_PVQZ")
     value = "4";
   else{
    printf("Basis set not recognize CBS_BASIS1 = %s\n",value.c_str());
    exit(0);
   }
   CBS_basis1 = StringToInt(value);
 }
 else if(param=="CBS_BASIS2"){
   value = StringToUpper(value);
   if(value == "2" ||value == "DZ" || value == "ADZ" || value == "CC-PVDZ" || value == "AUG-CC-PVDZ")
     value = "2";
   else if(value == "3" || value == "TZ" || value == "ATZ" || value == "CC-PVDZ" || value == "AUG-CC-PVTZ")
     value = "3";
   else if(value == "4" || value == "QZ" || value == "AQZ" || value == "CC-PVQZ" || value == "AUG-CC_PVQZ")
     value = "4";
   else{
    printf("Basis set not recognize CBS_BASIS2 = %s\n",value.c_str());
    exit(0);
   }
   CBS_basis2 = StringToInt(value);
 }
 else if(param=="CCSDT_CORRECTION"){
   value = StringToUpper(value);
   if(value=="FALSE" || value=="0") value = "0";
   else value = "1";
   CCSDT_correction = StringToInt(value);
 }
 else if(param=="SINGLE_FILE_MONOMER_HESS"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   single_file_monomer_hess = StringToInt(value);
 }else if (param=="MEMORY"){
   memory = StringToInt(value);
 }
  //If using MolPro, Equations for counterpoise correction needs to be defined in its rem section
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
   if(value == "1")
     do_local_two_body_truncation = StringToInt(value); //If periodic, local truncation is implied else lookout for "local_2_body" parameters
   //if space symmetry is not already set
   if(!symmetry_set){
     space_symmetry  = StringToInt(value);
     lattice_symmetry  = StringToInt(value);
   }
 }
  
 else if (param=="READ_LATTICE_VECTORS") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   read_lattice_vectors = StringToInt(value);
 }
  
 else if (param=="FREEZE_UNITCELLPARAMS") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   freeze_unitcellparams = StringToInt(value);
 }
  
 else if(param=="FREEZE_ATOM_COORDINATES"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   freeze_atom = StringToInt(value);
 }
  
 else if(param=="UNFREEZE_LATTICE_PARAMS_AFTER"){
   freeze_params_until = StringToInt(value);
   //printf("freeze_params_until = %i\n",freeze_params_until);
   if(freeze_params_until >= 0)
     freeze_unitcellparams = true;
   
  }
 else if(param=="CHANGE_VOLUME"){
   change_volume = StringToDouble(value);
   do_change_volume = true;
 }
 else if(param=="CHANGE_A"){
   change_a = StringToDouble(value);
   do_change_volume = true;
 }
 else if(param=="CHANGE_B"){
   change_b = StringToDouble(value);
   do_change_volume = true;
 }
 else if(param=="CHANGE_C"){
    change_c = StringToDouble(value);
    do_change_volume = true;
 }
 else if(param=="CHANGE_ALPHA"){
   change_alpha = StringToDouble(value);
   do_change_volume = true;
 }
 else if(param=="CHANGE_BETA"){
   change_beta = StringToDouble(value);
   do_change_volume = true;
 }
 else if(param=="CHANGE_GAMMA"){
   change_gamma = StringToDouble(value);
   do_change_volume = true;
 }
  
  // Crystal symmetry
 else if (param=="SYMMETRY") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   crystal_symmetry = StringToInt(value);

   //if space symmetry is not already set
   if(!symmetry_set)
     space_symmetry = StringToInt(value);
 }
  // space symmetry
 else if(param=="SPACE_SYMMETRY"){
   value = StringToUpper(value);
   if(value=="FALSE" || value=="0") value = "0";
   else value = "1";
   space_symmetry = StringToInt(value);
   symmetry_set = true;
 }
  //for symmetry
 else if (param=="SYMMETRY_TOLERANCE"){
   Symmetry_Tolerance = StringToInt(value);
 }
  //Lattice Params Symmetry
 else if (param=="LATTICE_SYMMETRY"){
   value = StringToUpper(value);
   if(value=="FALSE" || value =="0") value = "0";
   else value = "1";
   lattice_symmetry = StringToInt(value);
 }
  //Print symmetry info
 else if(param=="PRINT_SYMMETRY_INFO"){
   value = StringToUpper(value);
   if(value=="FALSE" || value =="0") value = "0";
   else value = "1";
   print_symmetry_info = StringToInt(value);
 }
  //monomer symmetry
 else if (param=="MONOMER_SYMMETRY"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";   
   monomer_symmetry = StringToInt(value);
 }
  //MM symmetery
 else if (param == "AIFF_SYMMETRY"||param=="MM_SYMMETRY"){
   value = StringToUpper(value);
   if(value=="FALSE" || value =="0") value = "0";
   else value = "1";
   mm_symmetry = StringToInt(value);
 }
  //Old Symmetry Names for Dimers
 else if (param == "OLDDIMERSYMM"||param=="OLD_DIMER_SYMM"){
   value = StringToUpper(value);
   if(value=="FALSE" || value =="0") value = "0";
   else value = "1";
   old_dimer_symm_names_ = StringToInt(value);
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
 else if (param=="SCRATCH_PATH") path_scratch = value;
 else if (param=="TMP_FILES_PATH") {
   path_tmp_files = base_path + value;
   path_qm_tmp_files = path_tmp_files + "/qm";
   path_mm_tmp_files = path_tmp_files + "/mm";
 }
 else if (param=="HESSIAN_FILES_PATH") {
   path_hessian_files = base_path + value;
   path_hessian_qm = path_hessian_files + "/qm";
   path_hessian_mm = path_hessian_files + "/mm";
 }
  
  // Print Level
 else if (param=="IPRINT" || param=="PRINT_LEVEL") IPRINT = StringToInt(value);
  
  // QM type
 else if (param=="QM_CODE"){
   value = StringToUpper(value);
   QM_Package_ = value;
   if (value== "QCHEM" ||value== "Q-CHEM") {
	QM_Type = 1;
   }	
   else if (value== "MOLPRO") {
     QM_Type = 2;
   }
   else if (value== "G09") {
     QM_Type = 3;
   }
   else if (value== "DALTON") {
     QM_Type = 4;
   }
   else if (value== "QE"|| value== "QUANTUMESPRESSO" || value=="QUANTUM_ESPRESSO") { //JLM
     QM_Type = 5;
   }
   else if (value== "ORCA") {
     QM_Type = 6;
   }
   else if (value== "PSI4") { //CSG
     QM_Type = 7;
   }
   else if (value== "DFTB" || value== "DFTB+" || value== "DFTB_PLUS" || value== "DFTBPLUS") {
     QM_Type = 8;
     use_dftb_freq_ = 1;
   }
   else {printf("Unknown QM_Type: %s\n",value.c_str()); exit(1);}
 }

  // MM type
 else if (param=="MM_CODE") {
   value = StringToUpper(value);
   if (value=="TINKER") MM_Type = 1;
   else if (value=="ORIENT" || value=="AIFF") MM_Type = 2;
   else if (value=="QCHEM") MM_Type = 3;
   else if (value=="EE-PA") MM_Type = 4;
   else if (value=="GDMA") MM_Type = 99; //JDH
   else if (value=="CHELPG") MM_Type = 98; //JDH
   else if (value=="HIRSHFELD") MM_Type = 97; // JDH
   else if (value=="CRYSTAL09"||value=="CRY2K9"){
     MM_Type = 5;
     MM_Only = true;
   }
   else {printf("Unknown MM_Type: %s\n",value.c_str()); exit(1);}
 }
  
  // Geometry optimization
 else if (param=="MAX_OPT_CYCLES") MaxOptCycles = StringToInt(value);
 else if (param=="OPTIMIZER") {
   value = StringToUpper(value);
   if (value=="DLFIND") {
     UseDLFind_ = true;
     UseKNitro_ = false;
     UseConjugate_ = false;
     UseSteepestDescent_ = false;
     UseLBFGS_ = false;
   }
   else if (value=="KNITRO") {
     UseDLFind_ = false;
     UseKNitro_ = true;
     UseConjugate_ = false;
     UseSteepestDescent_ = false;
     UseLBFGS_ = false;
   }
   else if (value=="CONJUGATEGRADIENT" || value == "CG") {
     UseDLFind_ = false;
     UseKNitro_ = false;
     UseConjugate_  = true;
     UseSteepestDescent_ = false;
     UseLBFGS_ = false;
     use_local_optimizer_ = true;
   }
   else if (value=="STEEPESTDESCENT" || value == "SD" || value == "STEEPEST-DESCENT") {
     UseDLFind_ = false;
     UseKNitro_ = false;
     UseConjugate_  = false;
     UseSteepestDescent_ = true;
     UseLBFGS_ = false;
     use_local_optimizer_ = true;
   }
   else if (value=="L-BFGS" || value == "LBFGS") {
     UseDLFind_ = false;
     UseKNitro_ = false;
     UseConjugate_  = false;
     UseSteepestDescent_ = false;
     UseLBFGS_ = true;
     use_local_optimizer_ = true;
   }
   else {
     printf("Error: Unknown Optimizer: %s\n",value.c_str()); exit(1);
   }
 }
 else if (param=="STEPSIZE") {
   step_size_init_ = StringToDouble(value);
 }
 else if (param=="COUPLE_GRADIENTS") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   couple_gradients_ = StringToInt(value);
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
  //finite difference on lattice vectors
 else if (param=="FDIFF_VECTOR"){
   value = StringToUpper(value);
   if(value=="FALSE" || value=="0") value = "0";
   else value = "1";
   LatticeVectorFiniteDifference = StringToInt(value);
 }
 else if (param=="FD_TINKERHESSIAN") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   Do_fdTinkerHessian_ = StringToInt(value);
 }
 else if (param=="DO_FREQ_AFTER_OPT"){
   value = StringToUpper(value);
   if (value=="FALSE" ||value=="0") value = "0";
   else value = "1";
   do_freq_after_opt = StringToInt(value);
 }
 else if (param=="FDIFF_HESS"){
   value = StringToUpper(value);
   if(value=="FALSE" || value=="0") value = "0";
   else value = "1";
   DoFiniteDifferenceFreqs_ = StringToInt(value);
 }
 else if (param=="ENERGY_FDIFF_HESS"){
   value = StringToUpper(value);
   if(value=="FALSE" || value=="0") value = "0";
   else value = "1";
   DoEnergyFiniteDifferenceFreq_ = StringToInt(value);
 }
  
 else if (param=="FINDSPATIALDAMPINGGRADIENT" || param=="FIND_SPATIAL_DAMPING_GRADIENT" ){
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   find_spatial_damping_gradient = StringToInt(value);
 }
 else if (param=="INPUT_FILENAME") {
   input_filename = value;
 }
  // Custom Basis Stuff JDH
 else if (param=="HBASIS" ) H_Basis_ = value;
 else if (param=="CBASIS" ) C_Basis_ = value;
 else if (param=="NBASIS" ) N_Basis_ = value;
 else if (param=="OBASIS" ) O_Basis_ = value;
 else if (param=="SBASIS" ) S_Basis_ = value;
 else if (param=="CLBASIS" ) Cl_Basis_ = value;
 else if (param=="IBASIS" ) I_Basis_ = value;
 else if (param=="SNBASIS" ) Sn_Basis_ = value;
 else if (param=="PBASIS" ) P_Basis_ = value;
 else if (param=="KBASIS" ) K_Basis_ = value;
 else if (param=="NABASIS" ) Na_Basis_ = value;
 else if (param=="BRBASIS" ) Br_Basis_ = value;
 else if (param=="FBASIS" ) F_Basis_ = value;
 else if (param=="WBASIS" ) W_Basis_ = value;
 else if (param=="VBASIS" ) V_Basis_ = value;
  // Gaussian stuff //
 else if (param=="G09") GaussianHeader_ = value;
 //else if (param=="G092") GaussianHeader2_ = value; 
  // Dalton stuff //
 else if (param=="DALTON") DaltonSection_ = value;
  // Orca stuff //
 else if (param=="ORCA") OrcaHeader_ = value;
  // PSI4 stuff //
 else if (param=="PSI4") PSI4Header_ = value;
 else if (param=="PSI4_REM") PSI4_rem = value;
  // QChem $rem section
 else if (param=="QC_REM") QC_rem = value;
  // QChem $basis section
 else if (param=="QC_BASIS") QC_basis = value;
  // Molpro $basis
 else if (param=="MOLPRO_REM") MolPro_rem = value;
  // Molpro calculations method
 else if (param=="MOLPRO_INST") MolPro_inst = value;
  // Molpro second $basis for CBS extrapolation
 else if (param=="MOLPRO_CBS") MolPro_CBS = value;
  // Instructions for HF in CBS extrapolation
 else if (param=="MOLPRO_HF") MolPro_HF = value;
  // CCSD(T) correction basis 
 else if (param=="MOLPRO_BASIS_CCSDT") MolPro_CCSDT_Basis = value;
  // MP2 for the CCSD(T) correction
  else if (param=="MOLPRO_CCSDT_MP2") MolPro_CCSDT_MP2 = value;
  else if (param=="MOLPRO_CCSDT_INST") MolPro_CCSDT_Inst = value;
  // Tinker $rem section - for keyfile
 else if (param=="TINKER_REM") Tinker_rem = value;
  // Orient $rem section --- by Ali
 else if (param=="AIFF_REM") AIFF_rem = value;
  // Secondary QChem $rem section
 else if (param=="QC_REM2") QC_rem2 = value;
  // Full parameter string - saved, just in case we want it later
 else if (param=="CRYSTAL_HEADING") crystal_heading = value;
 else if (param=="CRYSTAL_BASIS") crystal_basis = value;
 else if (param=="CRYSTAL_ENDING") crystal_ending = value;
 else if (param=="HMBI_PARAMS") HMBI_rem = value; 
 else if (param=="QE_SUPERCELL") QESupercellSection_ = value;
 else if (param=="QE_SPECIES") QESpeciesSection_ = value;
 else if (param=="QE_BASIS") qe_basis_ = StringToDouble(value);
 else if (param=="QE_BASIS_MULT" || param=="QE_BASIS_MULTIPLIER")  qe_basis_multiplier_ = StringToDouble(value);
 else if (param=="DISP_CORRECTION" || param=="DISPERSION_CORRECTION"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   dispersion_correction_ = StringToInt(value);
 }
 else if (param=="DISP_TYPE" || param=="DISPERSION_TYPE"){
   value = StringToUpper(value);
   dispersion_type_ = value;
 }
 else if (param=="FORCE_XC"){
   value = StringToUpper(value);
   force_xc_set_ = true;
   force_xc_ = value;
 }
 else if (param=="KPOINTS") kpoints_ = value;
 else if (param=="HUBBARD_DERIVS" || param=="HUBBARD_DERIVATIVES") hubbard_derivs_ = value;
 else if (param=="MAX_ANGULAR_MOMENTUM" | param=="MAXANGULARMOMENTUM" | param=="MAX_ANG_MOM") max_angular_momentum_ = value;
 else if (param=="SLATER_KOSTER") slater_koster_ = value;
 else if (param=="MOLECULE"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   molecule_ = StringToInt(value);
 }
 else if (param=="NO_FREQ_DISP" || param=="NO_FREQ_DISPLACEMENT"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   no_freq_disp_ = StringToInt(value);
 }
 else if (param=="PHONOPY" || param=="USE_PHONOPY"|| param=="USE_PHONOPY_FREQ"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   use_phonopy_freq_ = StringToInt(value);
 }
 else if (param=="PHONON" || param=="USE_PHONON"|| param=="USE_PHONON_FREQ"){
   value = StringToUpper(value);
   string val2;
   if (value=="FALSE" || value=="0"){ 
    value = "0";
    val2 = "1";
   }
   else{
     value = "1";
     val2 = "0";
   }
   use_phonon_freq_ = StringToInt(value);
   use_phonopy_freq_ = StringToInt(val2);
 }
 else if (param=="DFTB" || param=="USE_DFTB"|| param=="USE_DFTB_FREQ"){
   value = StringToUpper(value);
   string val2;
   if (value=="FALSE" || value=="0"){ 
    value = "0";
    val2 = "1";
   }
   else{
     value = "1";
     val2 = "0";
   }
   use_dftb_freq_ = StringToInt(value);
   use_phonopy_freq_ = StringToInt(val2);
 }
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
  
  // option for multiple damping factors in AIFF (added on 6th june, 2010)
  
 else if (param=="N_MONOMER_TYPES" ||
          param=="NUMBERMONOMERTYPES") {
   NumberMonomerTypes = StringToInt(value);            
   delete [] montype;
   montype = new string[NumberMonomerTypes];
   //montype.Initialize(1,NumberMonomerTypes);
   DimerDampingFactor.Initialize(NumberMonomerTypes,NumberMonomerTypes);
   for (int i=0;i<NumberMonomerTypes; i++) {
     for (int j=0;j<NumberMonomerTypes; j++) {
       DimerDampingFactor.Element(i,j) = 1.4; //default
     }
   }
   //montype.Print("montype 1");
   //DimerDampingFactor.Print("ddf 1");
   //fflush(stdout);
   
   ifstream infile1;
   infile1.open( input_filename.c_str() );
   if ( !infile1.is_open() ) {  
     printf("Params::SetParameters(): Cannot open file '%s'\n",
	    input_filename.c_str());
     exit(1);
   }                  
   string line;
   // Rewind the file, just in case
   infile1.clear();
   infile1.seekg(0,ios::beg);
   
   bool found = false; // flag for if aiff damping section found
   //stringmatrix montype has to be defined previously
   //int NoMonTypes = Params::Parameters().GetNumberMonomerTypes();
   
   // Start reading the file
   while ( !infile1.eof() ) {
     getline(infile1,line);
     int count = 0;
     if (line.substr(0,21)=="$damping_factors_aiff" ) {
       found = true;

       getline(infile1,line);
       
       // get the monomer type strings    
       istringstream iss(line);
       string mon_el;
       for (int j=0;j<NumberMonomerTypes;j++) {
	 iss >> mon_el;
	 montype[j] = mon_el;
	 fflush(stdout);
       }
       for (int k=0; k<( NumberMonomerTypes*(NumberMonomerTypes+1)/2 ); k++) {
	 fflush(stdout);
	 getline(infile1,line);
	 if (line.substr(0,4)!="$end") {
	   // get individual damping factors
	   istringstream iss(line);
	   string  montypeA,  montypeB;
	   int montypeA_index, montypeB_index;
	   string df;                
	   iss >> montypeA;
	   iss >> montypeB;
	   fflush(stdout);
	   
	   for (int i=0;i<NumberMonomerTypes;i++) {
	     if (montypeA.compare( montype[i] ) == 0) {
	       //              if (montypeA.compare( montype.Element(0,i) ) == 0) { 
	       montypeA_index = i;
	     }
	   }
	   for (int i=0;i<NumberMonomerTypes;i++) {
	     if (montypeB.compare( montype[i] ) == 0) {
	       //              if (montypeB.compare( montype.Element(0,i) ) == 0) {
	       montypeB_index = i;
	     } 
	   }
	   iss >> df;
	   DimerDampingFactor.Element( montypeA_index , montypeB_index ) =  StringToDouble(df);
	   DimerDampingFactor.Element( montypeB_index , montypeA_index ) =  StringToDouble(df);
	   
	 }
	 else {
	   goto Point1;           
	 }
	 //DimerDampingFactor.Print("ddf final");
       }
     }               
   }
 Point1:
   infile1.clear();
   //for (int i=0;i<NumberMonomerTypes;i++) {
   //  cout << montype[i] << "\t";
   //}
   //printf("\n\n");
 }
  
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
 else if (param=="SUPERCELL_JOB") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   supercell_job = StringToInt(value);
 }
 else if (param=="CORRECT_MP2_DISPERSION") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   MP2DispersionCorrection_ = StringToInt(value);
 }
 else if (param=="SUPERCELL_ANALYZE_ONLY") {
   value = StringToUpper(value);                 
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   supercell_job_exist = StringToInt(value);
 }
 else if (param=="MONOMER_FREQUENCY") {
   monomer_frequency_ = StringToInt(value);
 }
  
 else if (param=="SUPERCELL_SIZE_A") {
   supercell_size[0] = StringToInt(value);
 }
 else if (param=="SUPERCELL_SIZE_B") {
   supercell_size[1] = StringToInt(value);
 }        
 else if (param=="SUPERCELL_SIZE_C") {
   supercell_size[2] = StringToInt(value);
 }        
 else if (param=="THERMAL_USING_PHONONS") {
   value = StringToUpper(value);                 
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   IsThermalUsingPhonons = StringToInt(value);
 }
 else if (param=="THERMAL_FOR_MULTIPLE_TEMPERATURES"||
          param=="THERMAL_FOR_MULTIPLE_TEMPERATURE"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   ThermalForMultiTemp = StringToInt(value);
 }
 else if(param=="CREATE_SUPERCELL_INPUT"){
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   create_supercell_input = StringToInt(value);
 }
 else if(param=="DEUTERATED"||param=="DEUTERATE"){
   value = StringToUpper(value);
   if(value =="FALSE"||value=="0") value = "0";
   else value = "1";
   deuterated = StringToInt(value);
 }
  
  // Debugging
 else if (param=="TINKER_DEBUG") {
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   tinker_debug_ = StringToInt(value);
 }
 else if (param=="ADVANCED_EMBEDDING_SCHEME") { //JDH NEW
   AdvancedEmbeddingScheme_ = StringToUpper(value);
 }
 else if ( param=="MIXED_BASIS_LEVEL_1") { // JDH NEW 1/2017
   MixedBasisLevel1_ = value;
 }
 else if ( param=="MIXED_BASIS_LEVEL_2") { // JDH NEW 1/2017
   MixedBasisLevel2_ = value;
 }
 else if ( param=="MIXED_BASIS_LEVEL_3") { // JDH NEW  1/2017
   MixedBasisLevel3_ = value;
 }
 else if (param=="USE_ELECTROSTATIC_EMBEDDING") { //JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   UseElectrostaticEmbedding_ = StringToInt(value);
 }
 else if (param=="USE_TWO_BODY_CHARGE") { //JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   UseTwoBodyCharge_ = StringToInt(value);
 }
 else if (param=="CALCULATE_EFG") { //JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   EFGJob_ = StringToInt(value);
 }
 else if (param=="USE_SCALED_EFG_TENSORS") { //JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   UseScaledEFGTensors_ = StringToInt(value);
 } 
 else if (param=="USE_EWALD") { //JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   UseEwald_ = StringToInt(value);
 }
 else if (param=="EWALD_CUTOFF") { //JDH  (DEBUG Do i really need this?)
    Ewald_Cutoff_ = StringToDouble(value); //JDH
 }
 else if (param=="RCOND") { //JDH
   Rcond_ = StringToDouble(value); //JDH
 } 
 else if ( param=="CHARGE_EMBEDDING_RANK") {  //JDH 
   charge_embedding_rank_ = StringToInt(value);
 }
 else if (param=="NUMBER_OF_ASYMMETRIC_MONOMERS") { //JDH
   Num_Asymmetric_Monomers_ = StringToInt(value); //JDH
 }
 else if (param=="CUSTOM_BASIS") { // JDH
   Custom_Basis_ = StringToInt(value);
 } 
 else if (param=="MIXED_BASIS_CUT_OFF_1") {  // JDH 1/2017
   MixedBasisCutOff_ = StringToDouble(value);
 }
 else if (param=="MIXED_BASIS_CUT_OFF_2") {  // JDH 1/2017
   MixedBasisCutOff2_ = StringToDouble(value);
 }
 else if (param=="MIXED_BASIS_CUT_OFF_3") {  // JDH 1/2017
   MixedBasisCutOff3_ = StringToDouble(value);
 }
 else if (param=="TWO_BODY_CUTOFF") { //JDH
    TwoBodyCutoff_ = StringToDouble(value);
 }
 else if (param=="CLUSTER_CUTOFF") { //JDH 1/2017
   ClusterCutoff_ = StringToInt(value);
 }
 else if (param=="SELF_CONSISTENT_EMBEDDING") { // JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   UseSelfConsistentEmbedding_ = StringToInt(value);
 } 
 else if (param=="SELF_CONSISTENT_MAX_THRESHOLD") { //JDH
   SelfConsistentMaxThreshold_ = StringToDouble(value);
 }
 else if (param=="SELF_CONSISTENT_RMS_THRESHOLD") { //JDH
   SelfConsistentRMSThreshold_ = StringToDouble(value);
 }
 else if (param=="SEED_SCE") { // JDH
   value = StringToUpper(value);
   if (value=="FALSE" || value=="0") value = "0";
   else value = "1";
   SeedSCE_ = StringToInt(value);
 }
 else if (param=="ELECTROSTATIC_EMBEDDING_CUTOFF") { //JDH NEW 1/2017
   ElectrostaticEmbeddingCutoff_ = StringToDouble(value);
 }
  else if (param=="STUPIDVAR" || param=="STUPID_VARIABLE" || param=="STUPIDVARIABLE") { //JLM
    //value = StringToUpper(value);
    stupid_var_ = StringToDouble(value);
    /*if (value == "FALSE" || value == "0")
      stupid_var_ = false;
    else{
      stupid_var_ = true;
    }*/
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
