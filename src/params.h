#ifndef _params_h
#define _params_h
#include <cmath> //yoni : for symmetry tolerance
#include <string>
#include "matrix.h" 
using std::string;

/*
  A class that store the HMBI job parameters.
  
  GJB 01/09from HMBI

  Transformed to a singleton class, 11/09, GJB
 */


/* Instructions for Adding Keywords to the $hmbi input section
GJB 11/09

Step 1: Create a variable for the keyword in the class Params
        (params.h).  Make the name descriptive, and end it with an
        underscore.

Step 2: Create a member function to retreive the value of this
        variable inside params.h (usually inline functions are
        sufficient).

Step 3: (optional) Create a member function to modify this variable.
        Only necessary if the value can change during a calculation.

Step 4: Set a default value for the variable in the Constructor
        Params::Params() in params.C

Step 5: Add logic to Params::SetParameter (in params.C) to process the
        input file string.  Note, the input keyword does not
        necessarily have to be identical to the internal variable
        name, but make sure the names are meaningful and specific.

*/

class Params {
 private:
  Params();  // The real work is done by ReadInputFile
  ~Params(); 
  
  int Nproc; // # of processors used to run the job
  int PID_; // Unix Process ID, used to create unique scratch directories

  int IPRINT; // print level, higher = more output

  int Warnings; // A counter to track warnings that pop up and to notify
  // the user at the end.

  string base_path; // current working directory
  string path_QM, path_MM, path_scratch; // QM, MM, and Scratch paths
  string path_tmp_files; // store hessian, ykmat, skmat at this path
  string path_qm_tmp_files, path_mm_tmp_files; //to store approx hessian qm/mm files in tmp folder 
  string path_hessian_files; // store hmbi hessian files during freq calculation
  string path_hessian_qm, path_hessian_mm; // qm and mm folders for freq calculation
  
  string constraint;
  int NumberOfOptConstraints; // number of optimization constraints
  double constant_volume; //constant volume for constant volume optimization in (Angstroms)^3
  bool Find_Enthalpy; //perform enthalpy caluations. Currently for periodic systems only
  double ExternalPressure; //external pressure for constant pressure optimization in
                           //(user input in giga-pascals but atomic unit is hartrees/bohr^3)
  double Temperature; //in kelvin; Used to calculate vibrational energy contribution
  //Cameron
  double tolerance_max_g; //Optimization convergence tolerance for the maximum gradient 
  double tolerance_rms_g; //Optimization convergence tolerance for the RMS gradient
  double tolerance_e; //Optimization convergence tolerance for the change in energy

  bool Do_QuasiHarmonic;//approximate dependance of frequencies on values
  bool Volume_Gruneisen;//use volume gruneisen parameter or separator
                         //for each lattice parameters
  bool Save_QuasiHarmonic_Freq;// values the frequencies calculations are 
                          // saved rather an being written over.
  int QuasiStepsDone;//If Quasiharmonic steps are partially done.
                      //If gruneisen is isotropic, there are 3 steps, otherwise there are 6 steps (not including steps for lattice angles)
  bool ReadQuasiGeometries;//If Plus or minus geometries in quasiarharmonic already existing

  bool Analyze_Only; // whether to run jobs, or just read old results
  bool CreateJobsOnly_; // if true, exit after creating QM & MM jobs
  bool Neglect_Many_Body; // if true, do only 1+2 body at high level, no MB terms
  bool MM_Only; // Only the full MM is calulated.
  bool QM_Only; // Only the full QM is calulated. Currently intended only for planewave DFT in Quantum Espresso
  bool Use_Embedding_Charges; 
  bool Read_Cif_File; // Watit
  int Convergence; // Watit 
  bool Counterpoise; // whether or not doing BSSE correction for dimers
  bool DoQMBenchmark_; // whether or not to do full cluster calc at high level
                       // for benchmarking purposes

  /* 
     jobtype: 
     1 = single point energy
     2 = nuclear gradient 
     3 = nuclear hessian (not implemented)
     
     101 = expand geometry, run nothing
     
     In general: 1-99 normal jobs
                 100+ unusual jobs

  */
  int jobtype; 
  string jobtype_string; 

  bool Do_Forces; // whether or not to read force files
  bool Do_Freq; // whether or not to read hessians
  bool UseFiniteDifferenceGradients_; // Finite diff gradients?
  bool LatticeVectorFiniteDifference; //Finite difference on lattice vectors
  bool DoFiniteDifferenceFreqs_; // Finite diff freqs?
  bool DoEnergyFiniteDifferenceFreq_;
  bool Do_fdTinkerHessian_; //use finite difference
  bool do_freq_after_opt; // logic for freq after opt
  bool Is_Opt_Complete; //is optimization done??
  bool Are_Forces_Available; // bool for doing freq calc given forces
  bool Are_Hessians_Available; 
  bool Are_QHF_Available;// are QuasiHarmonic available
  bool Initialize_Hessians;
  bool get_lattice_hessian;
  bool get_elastic_constants;
  bool Do_Raman; // Watit
  bool couple_gradients_; //whether or not to couple lattice gradients with atomic gradients
  bool use_local_optimizer_;
  double stupid_var_;

  int MaxOptCycles; // maximum number of geometry optimization steps
  int OptCycle_; // counter for opt cycles
  int UseDLFind_; // Use DL-FIND optimizer
  int UseKNitro_; // Use K-NITRO optimizer
  int UseConjugate_; // Use internal Conjuated Gradient optimizer 
  int UseSteepestDescent_; // Use internal Steepest Descent optimizer 
  int UseLBFGS_; // Use internal L-BFGS optimizer 

  int monomer_frequency_; // if = 0, do full hessian.  if = n > 0, do hessian
                        // for only monomer n.

  int QM_Type;  /* QM_Type:  1=QChem, 2=Molpro 3=G09 4=DALTON 5=Quantum Espresso 6=ORCA 7=PSI4 8=DFTB+*/
  int MM_Type;  /* MM_Type:  1=TINKER, 2=AIFF, 3=QChem, 55=GDMA 56=ChelpG (JDH) */
  bool freeze_unitcellparams;
  bool freeze_atom;//for optimzation of the lattice parameters alone
  int freeze_params_until;//freezes the lattice parameters until a given number of cycles
  bool do_change_volume; //bool for changing volume, 
  double change_volume;//changes the volume of the unit cell isotropically while preserving symmetry
  double change_a;// change lattice parameter while preserving symmetry
  double change_b;// change lattice parameter while preserving symmetry
  double change_c;// change lattice parameter while preserving symmetry
  double change_alpha;// change lattice parameter while preserving symmetry
  double change_beta;// change lattice parameter while preserving symmetry
  double change_gamma;// change lattice parameter while preserving symmetry

  // for periodic boundary conditions (PBC)
  bool Periodic; // control flag
  bool press_set;
  bool temp_set;
  bool read_lattice_vectors; // whether to read angles or vectors for
			     // lattice
  bool crystal_symmetry; // flag for using symmetry with PBC
  bool space_symmetry; //flag for using space symmetry
  int Symmetry_Tolerance;//10^-int
  bool monomer_symmetry;// flag for using monomer symmetry
  bool mm_symmetry;//flag for using symmetry for the molecular mechanics
                   //currently only used of the AIFF
  bool lattice_symmetry;//symmetry is used for the lattice parameters
  bool print_symmetry_info;
  bool old_dimer_symm_names_;

  //This bool allows for different default symmetry settings for periodic and non-periodic systems.
  bool symmetry_set;

  int memory;//memory used by QM, only implemented for molpro

  bool tin_foil_boundary_conds; // if true, use tin-foil boundary
				// conds (infinite dielectric)
  // otherwise, use vacuum boundary conds (dielectric = 1).
  // Unimportant if unit cell lacks a dipole.  But for polar cells,
  // prefer tin-foil boundaries.

  int NumberMonomerTypes;//from Kaushiks code added by yoni
  string *montype;
//  StringMatrix montype;
//  DampFac DampFactor();
//  int NumberMonomerTypesSq = NumberMonomerTypes*NumberMonomerTypes;
  Matrix DimerDampingFactor;//(int NumberMonomerTypes);

  // for local 2-body truncations.  Cutoff distances are in Angstroms
  bool do_local_two_body_truncation; // logical flag to turn
				     // truncation on/off
  double cutoff0, cutoff1; // parameters for the truncation damping
			   // region
  string input_filename;  //yoni : added from Kaushik's code
  bool c0_set_explicitly; // was c0 set by the user explicitly?
  bool find_spatial_damping_gradient;

  // New variables to allow for scanning energies as a function of cutoff
  bool scan_local_cutoffs;
  double min_cutoff1; 
  int N_cutoff_steps;

  // AIFF cutoff parameters
  double polarization_cutoff; // max interaction radius for PBC
			      // induction			      
  double induction_convergence; // induction convergence threshold 
  double induction_gradconvergence; // induction gradient convergence threshold
  double induction_iter_scaling; // 0->1, slow down induction
				 // iterations.  Value = fraction of
				 // iteration step to mix with result
				 // from the previous cycle.
  bool preconverge_unit_cell_induced_moments; // if true, find induced
					      // multipole moments in
					      // non-periodic unit
					      // cell first, use them
					      // as a guess
  double two_body_dispersion_cutoff; // max interaction radius for PBC
				     // 2-body dispersion
  double three_body_dispersion_cutoff; // max interaction radius for
				       // PBC 3-body dispersion

  // Turn on/off different AIFF parts
  bool do_aiff_electrostatics; // if electrostatics is off, so is
			       // induction
  bool do_aiff_induction; 
  bool do_aiff_2body_dispersion;
  bool do_aiff_3body_dispersion;

  // AIFF Ewald parameters & cutoffs
  double Ewald_accuracy; // GJB, used with Ewald_kappa to determine Ewald cutoffs
  double Ewald_kappa; // partitioning factor for real & reciprocal space

  int Recip_cutoffx;    // by shuhao
  int Recip_cutoffy;   // by shuhao
  int Recip_cutoffz;  // by shuhao
  
  int Direc_cutoffx;    // by shuhao
  int Direc_cutoffy;   // by shuhao
  int Direc_cutoffz;  // by shuhao

  // Determine which frequency-dependent polarizabilities we use for
  // calculating dispersion coefficients.
  bool use_diagonal_anisotropic_pols; // true = use avg of diagonal
				      // aniso pols, else use
				      // isotropic pols.

  // Rem sections
  string HMBI_rem; //HMBI $hmbi section, stored just in case
  string QC_rem; // Q-Chem $rem section
  string MolPro_rem; //MolPro $rem section, 
  string MolPro_inst; //MolPro, instructions on how to to calucations
                     //Need to be seperate from rem for counterpoise correction
                     //In an CBS calculation, this is the MP2 section
  string MolPro_CBS; //second basis set for the complete basis set extrapolation
  string MolPro_HF; //$rem for the HF in a CBS calculation.
  string MolPro_CCSDT_Basis;//basis for the CCSD(T) correction
  string MolPro_CCSDT_MP2;//MP2 section for the CCSD(T) correction
  string MolPro_CCSDT_Inst;//instruction for the CCSD(T) part of the CCSD(T) correction
  string PSI4_rem; //PsI4 $rem section //CSG

  string QC_basis; // Q-Chem $basis section
  string QC_rem2; // Secondary Q-Chem $rem section for "MM" jobs
  string Tinker_rem; // Tinker keyfile parameters
  string AIFF_rem; // AIFF parameters --- by Ali
  string crystal_heading;//Heading for the crystal09 input file
  string crystal_basis;//path to the basis set for crystal09 input file
  string crystal_ending;//Ending for the crystal09 input file

  //JLM Quantum Espresso variables
  string QESupercellSection_;
  string QESpeciesSection_;
  double qe_basis_;
  double qe_basis_multiplier_;
  bool molecule_; //A Parameter needed for quantum Espresso when you intend to have 
                  //only a molecule in the cell; not the whole periodic crystal.
  bool no_freq_disp_; //Does not displpace molecule if negative frequencies exist 
  bool use_phonopy_freq_;
  bool use_phonon_freq_;
  bool use_dftb_freq_;

  //////////////////////////////////////////////////////
  // JDH's functions for Magnetic property calculations:
  //////////////////////////////////////////////////////

  // Job control
  bool NMRJob_; 
  bool EFGJob_; 
  int Num_Asymmetric_Monomers_; //JDH
  double TwoBodyCutoff_;
  double ClusterCutoff_; //JDH
  bool UseScaledEFGTensors_; // JDH
  string QM_Package_; //JDH NEW support for Q-Chem, MolPro, dalton and G09

  // Electrostatic embedding
  bool UseElectrostaticEmbedding_;
  double ElectrostaticEmbeddingCutoff_;
  bool UseTwoBodyCharge_;

  string AdvancedEmbeddingScheme_; // EWALD, SCE, SCRMP, NO
  bool UseEwald_; // EWALD
  double Ewald_Cutoff_; // EWALD 
  double Rcond_; // EWALD (Probably shouldn't ever have to adjust this, but just in case...)

  bool UseSelfConsistentEmbedding_; // SCE
  double SelfConsistentMaxThreshold_; // SCE
  double SelfConsistentRMSThreshold_; // SCE
  bool SeedSCE_; // SCE
  
  int charge_embedding_rank_; //JDH
  
  
  // Custom Basis functions:
  int Custom_Basis_; // JDH 0 = no custum basis (default), 1 = custom basis on asymmetric monomer, 2 = custom basis for all instances of the specified atom
  double MixedBasisCutOff_; // JDH
  double MixedBasisCutOff2_; // JDH
  double MixedBasisCutOff3_; // JDH
  string MixedBasisLevel1_; // JDH (this is the large one...)
  string MixedBasisLevel2_; // JDH
  string MixedBasisLevel3_; // JDH
  string H_Basis_;
  string C_Basis_;
  string N_Basis_;
  string O_Basis_;
  string S_Basis_;
  string Cl_Basis_;
  string I_Basis_;
  string Sn_Basis_;
  string P_Basis_;
  string K_Basis_;
  string Na_Basis_;
  string Br_Basis_;
  string F_Basis_;
  string W_Basis_;
  string V_Basis_;
  

  // QM Package interfaces:
  string GaussianHeader_; //JDH
  string DaltonSection_; //JDH
  string OrcaHeader_; // JDH
  string Molpro1_; //JDH
  string Molpro2_; //JDH
  string PSI4Header_; //CSG

  
  // Functions flagged for removal:
  // double NShellsEwaldTestPoints_; // JDH (controls the resolution of the test point sampling)
  // int NShellsEwaldFitPoints_;
  // int EwaldNFitScaleFactor_;
  // bool ReadEwaldPotentialFromFile_;
  // double EwaldShellThickness_;
  // double EwaldProbeShellBuffer_;

  

    
  //JDH's OLD functions DELETE NOW
  //bool Keep_NMR_Files_; //JDH
  //bool NMRChargeEmbedding_; //JDH
  //double Periodic_NMR_Cutoff_Two_Body_; //JDH
  //double Inner_Charge_Embedding_Cutoff_; //JDH
  //double Outer_Charge_Embedding_Cutoff_; //JDH
  //double NMRClusterFragment_; // JDH
  //bool NMRChargeEmbedding_; //JDH
  //double Periodic_NMR_Cutoff_Two_Body_; //JDH
  //double Inner_Charge_Embedding_Cutoff_; //JDH
  //double Outer_Charge_Embedding_Cutoff_; //JDH


  
  


  //CBS limit extrapolation
  bool CBS;
  int CBS_basis1;
  int CBS_basis2;

  //CCSD(T) correction
  bool CCSDT_correction;

  bool single_file_monomer_hess;

  double ExpansionFactor; // scale factor for adjusting intermolecular spacing
  // ExpansionFactor < 1 shrinks it, and > 1 expands it.

  string AIFFBasisSet; // by Ali
  string CamCaspHome; // by Ali
  double IonizationPotential;
  // Damping factor for the AIFF
  double DampingFactor;
  int MaxPolarizationCycles;
  bool OrientDebug_;

  bool UseGlobalCoords_;
  bool LocalCoordinatesInitialized_;
  bool BuildForceFieldOnly_;
  bool RunJobsOnly_;
  bool dispersion_correction_;
  string dispersion_type_;
  bool force_xc_set_;
  string force_xc_;

  int OptType_; // DL-FIND optimizer. 1 = Steepest Descent, 2 = Conj. Grad, 3 = L-BFGS
  int LineSearchType_; // DL-FIND line search type: 
  /*
                     0 = simple scaled
		     1 = trust radius based on energy
                     2 = trust radius based on gradient
                     3 = line search
  */

  bool use_approx_hessian; // do we calculate the approx hessian for lbfgs optimization??
  int CoordType_; // Coordinate representation, 0 = cartesian, 1-4 are various internal
  

  bool EstimateThreeBodyDispersion_; // Use simple ATM 3-body disp correction
  bool MP2DispersionCorrection_; // Correct MP2 dispersion

  //supercell info params
  bool supercell_job; //is this a supercell job??
  bool supercell_job_exist; //is this job only for analyzing?? 
  Vector supercell_size; //size of supercell
  bool create_supercell_input;//create an input of the supercell
  bool IsThermalUsingPhonons; // do we want thermal properties calculated using phonons
  bool ThermalForMultiTemp; // do we want thermal properties at multiple temperatures
  int TotalKPoint; 

  bool deuterated;


  bool tinker_debug_; // Do pure MM Tinker calculation

  //JLM opt variables
  double step_size_init_;
  string kpoints_;
  string hubbard_derivs_;
  string max_angular_momentum_;
  string slater_koster_;

 public:
  static Params& Parameters() {
    static Params params_;
    return params_;
  }
  
  //Basic functions to make life easy
  void SetParameter(string param, string value); 
  void Print();
  int StringToInt(string str);
  double StringToDouble(string str);
  string StringToUpper(string str);

  // Basic functions for accessing parameters
  int PrintLevel() {return IPRINT;};
  void Warning() {Warnings++;};
  int CheckWarnings() {return Warnings;};

  //MPI related 
  void SetNumberOfProcessors(int N) {Nproc = N;};
  int GetNumberOfProcessors() {return Nproc;};
  void SetPID(int pid) {PID_ = pid;};
  int GetPID() {return PID_;};

  //Job related
  int GetJobType() {return jobtype;};
  string GetJobTypeStr() {return jobtype_string;};
  bool RunJobs() {return !Analyze_Only;};
  bool CreateJobsOnly() {return CreateJobsOnly_;};
  bool RunJobsOnly() {return RunJobsOnly_;};
  bool BuildForceFieldOnly() {return BuildForceFieldOnly_;};
  bool FindEnthalpy() {return Find_Enthalpy;};
  double GetExternalPressure() {return ExternalPressure;};
  double GetTemperature() {return Temperature;};
  double GetMaxg() {return tolerance_max_g;};
  double GetMaxRms() {return tolerance_rms_g;};
  double GetETol() {return tolerance_e;};
  bool DoForces() {return Do_Forces;};
  bool DoFreq() {return Do_Freq;};
  bool UseFiniteDifferenceGradients() {return UseFiniteDifferenceGradients_;};
  bool DoLatticeVectorFiniteDifference() {return LatticeVectorFiniteDifference;};
  bool DoFiniteDifferenceFreqs() {return DoFiniteDifferenceFreqs_;};
  bool DoEnergyFiniteDifferenceFreqs() {return DoEnergyFiniteDifferenceFreq_;};
  bool Do_fdTinkerHessian() {return Do_fdTinkerHessian_;};
  bool DoFreqAfterOpt() { return do_freq_after_opt;};
  bool AreForcesAvailable() {return Are_Forces_Available;};
  bool AreHessiansAvailable() {return Are_Hessians_Available;};
  bool InitializeHessians() {return Initialize_Hessians;};
  bool EstimateThreeBodyDispersion() {return EstimateThreeBodyDispersion_;};
  bool DoMP2DispersionCorrection() {return MP2DispersionCorrection_;};
  bool DoCounterpoise() {return Counterpoise;};
  bool DoQMBenchmark() {return DoQMBenchmark_;};
  bool DispersionCorrection() {return dispersion_correction_;};
  string DispersionType() {return dispersion_type_;};
  bool ForceXCSet() {return force_xc_set_;};
  string ForceXCType() {return force_xc_;};
  bool IsMolecule() {return molecule_;};
  bool NoFreqDisplacement() {return no_freq_disp_;};
  bool UsePhonopyFreq() {return use_phonopy_freq_;};
  bool UsePhononFreq() {return use_phonon_freq_;};
  bool UseDFTBFreq() {return use_dftb_freq_;};
  bool DoRaman() {return Do_Raman;}; //Watit
  bool ReadCIFFile() {return Read_Cif_File;}; //Watit

  //Quasiharmonic approximation
  bool AreQHFAvailable(){return Are_QHF_Available;}//quasiharmonic frequencies
  bool DoQuasiHarmonic() {return Do_QuasiHarmonic;};
  bool UseVolumeGruneisen() {return Volume_Gruneisen;};
  bool SaveQuasiHarmonicCalculations() {return Save_QuasiHarmonic_Freq;};
  int NumberOfCompletedQuasiSteps() {return QuasiStepsDone;};
  bool QuasiReadQuasiGeometries() {return ReadQuasiGeometries;};
  bool GetThermalPropertiesUsingPhonons() {return IsThermalUsingPhonons;};
  bool GetThermalPropertiesAtMultipleTemperatures() {return ThermalForMultiTemp;};


  //////////////////////////////////////////////////////
  // JDH's functions for Magnetic property calculations:
  //////////////////////////////////////////////////////
  // Job control
  bool CalculateNMR() {return NMRJob_;};
  bool CalculateEFG() {return EFGJob_;};
  int GetNumAsymmetricMonomers() { return Num_Asymmetric_Monomers_; }; //JDH
  double GetTwoBodyCutoff() {return TwoBodyCutoff_;};
  double GetClusterCutoff() {return ClusterCutoff_;}; // JDH
  bool UseScaledEFGTensors() { return UseScaledEFGTensors_;}; // JDH
  string GetQMPackage() {return QM_Package_;}; // JDH

  // Electrostatic embedding
  bool   UseElectrostaticEmbedding() { return UseElectrostaticEmbedding_;};
  double GetElectrostaticEmbeddingCutoff() {return ElectrostaticEmbeddingCutoff_;};
  bool UseTwoBodyCharge() { return UseTwoBodyCharge_;};

  string GetAdvancedEmbeddingScheme(){ return AdvancedEmbeddingScheme_;};
  
  bool UseEwald() {return UseEwald_;}; // EWALD
  double GetEwaldCutoff() { return Ewald_Cutoff_; };  // EWALD
  double GetRcond(){ return Rcond_;};  // EWALD

  bool UseSelfConsistentEmbedding() { return UseSelfConsistentEmbedding_;}; // SCE
  double GetSelfConsistentMaxThreshold() { return SelfConsistentMaxThreshold_;}; // SCE
  double GetSelfConsistentRMSThreshold() { return SelfConsistentRMSThreshold_;}; // SCE
  bool SeedSCE() { return SeedSCE_; }; // SCE

  int GetChargeEmbeddingRank() {return charge_embedding_rank_;};
  

  // Custom Basis functions:
  int CustomBasis() {return Custom_Basis_;};
  double GetMixedBasisCutOff() { return MixedBasisCutOff_;};
  double GetMixedBasisCutOff2() { return MixedBasisCutOff2_;};
  double GetMixedBasisCutOff3() { return MixedBasisCutOff3_;};
  string GetMixedBasisLevel1() { return MixedBasisLevel1_;}; 
  string GetMixedBasisLevel2() { return MixedBasisLevel2_;}; 
  string GetMixedBasisLevel3() { return MixedBasisLevel3_;};
  string GetHBasis() {return H_Basis_;}; 
  string GetCBasis() {return C_Basis_;}; 
  string GetNBasis() {return N_Basis_;}; 
  string GetOBasis() {return O_Basis_;}; 
  string GetSBasis() {return S_Basis_;}; 
  string GetClBasis() {return Cl_Basis_;}; 
  string GetIBasis() {return I_Basis_;}; 
  string GetSnBasis() {return Sn_Basis_;}; 
  string GetPBasis() {return P_Basis_;};
  string GetKBasis() {return K_Basis_;};
  string GetNaBasis() {return Na_Basis_;};
  string GetBrBasis() {return Br_Basis_;};
  string GetFBasis() {return F_Basis_;};
  string GetWBasis() {return W_Basis_;};
  string GetVBasis() {return V_Basis_;};


  // QM Package interfaces:
  string GetGaussianHeader() {return GaussianHeader_;};
  string GetDaltonSection() {return DaltonSection_;};
  string GetOrcaHeader() {return OrcaHeader_;};
  string GetMolPro1() {return Molpro1_;}; 
  string GetMolPro2() {return Molpro2_;}; 
  string GetPSI4Header() {return PSI4Header_;};
  
  // JDH old  NMR stuff...
  /* double GetNMRClusterFragment() {return NMRClusterFragment_;}; */
  /* double GetNMRMixedBasisCutOff() { return NMRMixedBasisCutOff_;}; // JDH */
  /* double GetNMRMixedBasisCutOff2() { return NMRMixedBasisCutOff2_;}; // JDH */
  /* double GetNMRMixedBasisCutOff3() { return NMRMixedBasisCutOff3_;}; // JDH */
  /* string GetNMRMixedBasisLevel1() { return NMRMixedBasisLevel1_;}; // JDH */
  /* string GetNMRMixedBasisLevel2() { return NMRMixedBasisLevel2_;}; // JDH */
  /* string GetNMRMixedBasisLevel3() { return NMRMixedBasisLevel3_;}; // JDH */
  /* bool NMRuseChargeEmbedding() { return NMRChargeEmbedding_;}; // JDH */
  /* double GetPeriodicNMRCutoffTwoBody() {return Periodic_NMR_Cutoff_Two_Body_; }; //JDH */
  /* double GetClusterFragmentCutoff() {return Cluster_Fragment_Cutoff_;}; */
  /* double GetInnerEmbeddingCutoff() { return Inner_Charge_Embedding_Cutoff_; }; //JDH */
  /* double GetOuterEmbeddingCutoff() { return Outer_Charge_Embedding_Cutoff_; }; // JDH */
  /* int GetNumAsymmetricMonomers() { return Num_Asymmetric_Monomers_; }; //JDH */
  /* string GetQMPackage() {return QM_Package_;}; // JDH  */

  //Optimization related tags
  int GetMaxOptCycles() {return MaxOptCycles;};
  bool UseDLFind() {return UseDLFind_;};
  void SetUseDLFind(bool value) {UseDLFind_ = value;};
  bool UseApproxHessian() {return use_approx_hessian;};
  bool UseKNitro() {return UseKNitro_;};
  void SetKNitro(bool value) {UseKNitro_ = value;};
  bool UseConjugateGradient() {return UseConjugate_;};
  bool UseSteepestDescent() {return UseSteepestDescent_;};
  bool UseLBFGS() {return UseLBFGS_;};
  bool IsOptOver() {return Is_Opt_Complete;};
  int GetOptType() {return OptType_;};
  int GetLineSearchType() {return LineSearchType_;};
  int GetCoordType() {return CoordType_;};
  void IncrementOptCycle() {OptCycle_++;};
  int GetOptCycle() {return OptCycle_;};
  void ResetOptCycle() {OptCycle_ = 0;};
  double GetInitialStepSize() {return step_size_init_;};
  bool GetUseLocalOptimizer() {return use_local_optimizer_;};
  double GetStupidVariable() {return stupid_var_;};
  int GetConvergence() {return Convergence;}; //Watit


  //Symmetry
  bool UseCrystalSymmetry() {return crystal_symmetry;};
  double GetSymmetryTolerance() {return pow(10.0,-Symmetry_Tolerance);};
  bool UseSpaceSymmetry() {if(!crystal_symmetry) return 0;
                           else return space_symmetry;};
  bool UseMonomerSymmetry() {return monomer_symmetry;};
  bool UseMMSymmetry() {return mm_symmetry;};
  bool UseLatticeSymmetry() {if(!crystal_symmetry) return 0;
                              else if(space_symmetry) return 1;
			       else return lattice_symmetry;};
  bool PrintSymmetryInfo() {return print_symmetry_info;};
  bool OldSymmetryDimers() {return old_dimer_symm_names_;};

//Constraints
  string GetConstraintStr() {return constraint;};
  int GetNumberOfConstraints() { return NumberOfOptConstraints;};
  void SetConstantVolume(double volume) {constant_volume = volume;};
  double GetConstantVolume() { return constant_volume;};
  bool FreezeUnitCellParams() {return freeze_unitcellparams;};
  bool FreezeNuclearCoordinates() {return freeze_atom;};
  int UnFreezeParamsAfter() {return freeze_params_until;};
  bool ChangeVolume() {return do_change_volume;};
  double ChangeVolumeBy() {return change_volume;};
  double ChangeParameterA() {return change_a;};
  double ChangeParameterB()  {return change_b;};
  double ChangeParameterC()  {return change_c;};
  double ChangeParameterAlpha() {return change_alpha;};
  double ChangeParameterBeta() {return change_beta;};
  double ChangeParameterGamma() {return change_gamma;};

//Periodic or Long-range 
  bool IsPeriodic() {return Periodic;};
  bool ReadLatticeVectors() {return read_lattice_vectors;};
  bool UseEmbeddingCharges() {return Use_Embedding_Charges;};
  bool NeglectManyBody() {return Neglect_Many_Body;};
  bool UseFullMMOnly() {return MM_Only;};
  bool UseFullQMOnly() {return QM_Only;};
  bool DoLocal2BodyTruncation() {return do_local_two_body_truncation;};
  double GetLocalCutoff(int i) {if (i==0) return cutoff0; else return cutoff1;};
  bool FindSpatialDampingGradient(){return find_spatial_damping_gradient;};
  void SetLocalCutoff(int i, double value) {if (i==0) cutoff0=value; else cutoff1=value;};
  bool ScanLocalTruncationParameters() {return scan_local_cutoffs;};
  double GetMinLocalCutoff() {return min_cutoff1;};
  int GetNumberOfCutoffScanSteps() {return N_cutoff_steps;};
  bool TinFoilBoundaryConditions() {return tin_foil_boundary_conds;};
  bool IsSupercellJob() {return supercell_job;};
  bool Supercell_Analyze_Only() {return supercell_job_exist;};
  Vector GetSupercellSize() {return supercell_size;};
  void SetSupercellSize(Vector size) {supercell_size = size;};
  bool CreateSupercellInput(){return create_supercell_input;};
  //total k point for the super cell
  void SetSuperCellTotalKPoint(int total_k_point) {TotalKPoint = total_k_point;};
  int GetSuperCellTotalKPoint() {return TotalKPoint;};
  bool AreGradientsCoupled() {return couple_gradients_;};

//Other useful commands
  string GetBasePath() {return base_path;};
  string GetQMPath() {return path_QM;};
  string GetMMPath() {return path_MM;};
  string GetScratchPath() {return path_scratch;};
  void SetScratchPath(string path) {path_scratch = path;};
  string GetTmpFilesPath() {return path_tmp_files;};
  string GetTmpQMFilesPath() {return path_qm_tmp_files;};
  string GetTmpMMFilesPath() {return path_mm_tmp_files;};
  string GetHessianPath() {return path_hessian_files;};
  string GetHessianQMPath() {return path_hessian_qm;};
  string GetHessianMMPath() {return path_hessian_mm;};
  int GetQMType() {return QM_Type;};
  int GetMMType() {return MM_Type;};
  int MemoryUse() {return memory;};
  bool IsDeuterated() {return deuterated;};

  
  //Additional input information (usually copied from input file without alteration)
  string GetQChemRem() {return QC_rem;};
  string GetMolProRem() {return MolPro_rem;};
  string GetMolProInst() {return MolPro_inst;};
  //second basis set for a CBS extrapolation
  string GetMolProCBSRem() {return MolPro_CBS;}
  //calculations for the HF calculation in the CBS extrapolation
  string GetMolProHFRem() {return MolPro_HF;};
  string GetMolProCCSDTBasis(){return MolPro_CCSDT_Basis;};//basis for the CCSD(T) correction
  string GetMolProCCSDTMP2(){return MolPro_CCSDT_MP2;};//MP2 section for the CCSD(T) correction
  string GetMolProCCSDTInst(){return MolPro_CCSDT_Inst;};//instruction for the CCSD(T) part of the CCSD(T) correction
  string GetPSI4Rem() {return PSI4_rem;}; //CSG
  string GetQChemRem2() {return QC_rem2;};
  string GetQChemBasis() {return QC_basis;};
  string GetTinkerRem() {return Tinker_rem;};
  string GetAIFFRem() {return AIFF_rem;};
  string GetHMBIRem() {return HMBI_rem;};
  string GetCrystalHeading() {return crystal_heading;};
  string GetCrystalBasis() {return crystal_basis;};
  string GetCrystalEnding() {return crystal_ending;};
  string GetQuantumEspressoSupercellSection() {return QESupercellSection_;};  //JLM  
  string GetQuantumEspressoSpeciesSection() {return QESpeciesSection_;};  //JLM  
  double GetQEBasis() {return qe_basis_;};
  double GetQEBasisMultiplier() {return qe_basis_multiplier_;};
  string GetKPoints() {return kpoints_;};
  string GetHubbardDerivatives() {return hubbard_derivs_;};
  string GetMaxAngularMomentum() {return max_angular_momentum_;};
  string GetSlaterKoster() {return slater_koster_;};
  bool TinkerDebug() {return tinker_debug_;};

  //CBS extrapolation
  bool DoCBS() {return CBS;};
  int CBSBasis1() {return CBS_basis1;};
  int CBSBasis2() {return CBS_basis2;};

  //CCSD(T) correction
  bool DoCCSDTCorrection() {return CCSDT_correction;};
  bool SingleFileMonomerHess() {return single_file_monomer_hess;};

  double GetExpansionFactor() {return ExpansionFactor;};
  void SetMonomerForFrequency(int imon) {monomer_frequency_ = imon;};
  int GetFrequencyForMonomer() {return monomer_frequency_;};

  string GetAIFFBasisSet() {return AIFFBasisSet;};  // by Ali
  string GetCamCaspHome() {return CamCaspHome;};  // by Ali
  double GetIonizationPotential() {return IonizationPotential;};  // by Ali

  string* GetMonTypeArray() { return montype;};
  //StringMatrix GetMonTypeArray() { return montype;};
  Matrix GetDimerDampingFactor() {return DimerDampingFactor;};
  int GetNumberMonomerTypes() {return NumberMonomerTypes;};

  double GetDampingFactor() {return DampingFactor;};  // by Ali
  int GetMaxPolarizationCycles() {return MaxPolarizationCycles;};
  bool OrientDebug() {return OrientDebug_;};
  bool UseGlobalCoordinates() {return UseGlobalCoords_;};


  // Functions for working with AIFF force field parameters
  void SetMaxPolarizationRadius(double r) {polarization_cutoff = r;};
  double GetMaxPolarizationRadius() {return polarization_cutoff;};

  void SetMaxTwoBodyDispersionRadius(double r) {two_body_dispersion_cutoff = r;};
  double GetMaxTwoBodyDispersionRadius() {return two_body_dispersion_cutoff;};

  void SetMaxThreeBodyDispersionRadius(double r) {three_body_dispersion_cutoff = r;};
  double GetMaxThreeBodyDispersionRadius() {return three_body_dispersion_cutoff;};

  void SetInductionConvergence(double E) {induction_convergence = E;};
  double GetInductionConvergence() {return induction_convergence;};
  
  void SetInductionGradConvergence(double E) {induction_gradconvergence = E;};
  double GetInductionGradConvergence() {return induction_gradconvergence;};

  void SetInductionIterScaling(double s) {induction_iter_scaling = s;};
  double GetInductionIterScaling() {return induction_iter_scaling;};

  bool PreconvergeInducedMomentsInUnitCell() {return preconverge_unit_cell_induced_moments;};


  // Functions for checking which parts of AIFF to evaluate.  In general,
  // all should be true.
  bool DoAIFFElectrostatics() {return do_aiff_electrostatics;};
  bool DoAIFFInduction() {return do_aiff_induction;};
  bool DoAIFF2BodyDispersion() {return do_aiff_2body_dispersion;};
  bool DoAIFF3BodyDispersion() {return do_aiff_3body_dispersion;};

  void SetEwaldKappa(double value) {Ewald_kappa = value;};
  double GetEwaldKappa() {return Ewald_kappa;}; 

  void SetEwaldAccuracy(double value) {Ewald_accuracy = value;};
  double GetEwaldAccuracy() {return Ewald_accuracy;}; 

  void SetRecipLoopX(int Rx) {Recip_cutoffx = Rx;};
  int  GetRecipLoopX() {return Recip_cutoffx;};

  void SetRecipLoopY(int Ry) {Recip_cutoffy = Ry;};
  int  GetRecipLoopY() {return Recip_cutoffy;};
  
  void SetRecipLoopZ(int Rz) {Recip_cutoffz = Rz;};
  int  GetRecipLoopZ() {return Recip_cutoffz;};

  void SetDirecLoopX(int Dx) {Direc_cutoffx = Dx;};
  int  GetDirecLoopX() {return Direc_cutoffx;};

  void SetDirecLoopY(int Dy) {Direc_cutoffy = Dy;};
  int  GetDirecLoopY() {return Direc_cutoffy;};

  void SetDirecLoopZ(int Dz) {Direc_cutoffz = Dz;};
  int  GetDirecLoopZ() {return Direc_cutoffz;};


  bool UseDiagonalAnisotropicPolarizabilities() {return use_diagonal_anisotropic_pols;};


  /*
  bool LocalCoordinateAxesInitialized() {
    return LocalCoordinatesInitialized_;}
  void LocalCoordinateAxesInitialized(bool value) {LocalCoordinatesInitialized_ = value;};
  */


};

#endif
