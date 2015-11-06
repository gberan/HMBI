#ifndef _params_h
#define _params_h

#include <string>
using std::string;

/*
  A class that store the HMBI job parameters.
  
  GJB 01/09

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
  string path_QM, path_MM; // QM and MM paths

  bool Analyze_Only; // whether to run jobs, or just read old results
  bool CreateJobsOnly_; // if true, exit after creating QM & MM jobs
  bool Neglect_Many_Body; // if true, do only 1+2 body at high level, no MB terms
  bool Use_Embedding_Charges; 
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
  bool UseFiniteDifferenceGradients_; // Finite diff gradients?
  bool DoFiniteDifferenceFreqs_; // Finite diff freqs?
  int MaxOptCycles; // maximum number of geometry optimization steps
  int OptCycle_; // counter for opt cycles
  int UseDLFind_; // Use DL-FIND optimizer?
  int UseKNitro_; // Use K-NITRO optimizer?

  int monomer_frequency_; // if = 0, do full hessian.  if = n > 0, do hessian
                        // for only monomer n.

  //bool do_raman_intensities_; // if TRUE, compute Raman frequency intensities

  int MM_Type;  /* MM_Type:  1=TINKER, 2=AIFF, 3=QChem */

  // for periodic boundary conditions (PBC)
  bool Periodic; // control flag
  bool read_lattice_vectors; // whether to read angles or vectors for
			     // lattice
  bool read_fractional_coordinates; // input geom in fractional coords
  bool crystal_symmetry; // flag for using symmetry with PBC
  bool tin_foil_boundary_conds; // if true, use tin-foil boundary
				// conds (infinite dielectric)
  // otherwise, use vacuum boundary conds (dielectric = 1).
  // Unimportant if unit cell lacks a dipole.  But for polar cells,
  // prefer tin-foil boundaries.

  // for local 2-body truncations.  Cutoff distances are in Angstroms
  bool do_local_two_body_truncation; // logical flag to turn
				     // truncation on/off
  double cutoff0, cutoff1; // parameters for the truncation damping
			   // region
  bool c0_set_explicitly; // was c0 set by the user explicitly?

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

  double three_body_dispersion_print_threshold; // min size before
						// printing out a
						// 3-body contribution

  // Turn on/off different AIFF parts
  bool do_aiff_electrostatics; // if electrostatics is off, so is
			       // induction
  bool do_aiff_induction; 
  bool do_aiff_2body_dispersion;
  bool do_aiff_3body_dispersion;

  // AIFF Ewald parameters & cutoffs
  double Ewald_accuracy; // GJB, used with Ewald_kappa to determine Ewald cutoffs
  double Ewald_kappa; // partitioning factor for real & reciprocal space
  int Ewald_induction_type; // use large cluster (=1) or Madelung potential (=2) for Ewald?

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
  string QC_basis; // Q-Chem $basis section
  string QC_rem2; // Secondary Q-Chem $rem section for "MM" jobs
  string Tinker_rem; // Tinker keyfile parameters
  string AIFF_rem; // AIFF parameters --- by Ali

  double ExpansionFactor; // scale factor for adjusting intermolecular spacing
  // ExpansionFactor < 1 shrinks it, and > 1 expands it.

  string AIFFBasisSet; // by Ali
  string CamCaspHome; // by Ali
  double IonizationPotential;
  // Damping factor for the AIFF
  double DampingFactor;
  bool DoAtomicInductionDamping_; // use different values for each
				 // atom. Exploratory code.
  double AtomicInductionDampingFactor[93]; // store values for each element by atomic number
  bool DoDispersionDamping_;

  int MaxPolarizationCycles;
  bool OrientDebug_;

  bool UseGlobalCoords_;
  bool LocalCoordinatesInitialized_;
  bool BuildForceFieldOnly_;
  bool RunJobsOnly_;


  int OptType_; // DL-FIND optimizer. 1 = Steepest Descent, 2 = Conj. Grad, 3 = L-BFGS
  int LineSearchType_; // DL-FIND line search type: 
  /*
                     0 = simple scaled
		     1 = trust radius based on energy
                     2 = trust radius based on gradient
                     3 = line search
  */
  int CoordType_; // Coordinate representation, 0 = cartesian, 1-4 are various internal
  

  bool EstimateThreeBodyDispersion_; // Use simple ATM 3-body disp correction
  bool MP2DispersionCorrection_; // Correct MP2 dispersion



  bool tinker_debug_; // Do pure MM Tinker calculation

 public:

  // instantiate the single instance of this class
  static Params& Parameters() {
    static Params params_;
    return params_;
  }

  void SetParameter(string param, string value); 
  void Print();

  int StringToInt(string str);
  double StringToDouble(string str);
  string StringToUpper(string str);

  // Basic functions for accessing parameters
  int PrintLevel() {return IPRINT;};
  void Warning() {Warnings++;};
  int CheckWarnings() {return Warnings;};

  void SetNumberOfProcessors(int N) {Nproc = N;};
  int GetNumberOfProcessors() {return Nproc;};
  void SetPID(int pid) {PID_ = pid;};
  int GetPID() {return PID_;};
  int GetJobType() {return jobtype;};
  string GetJobTypeStr() {return jobtype_string;};
  bool RunJobs() {return !Analyze_Only;};
  bool CreateJobsOnly() {return CreateJobsOnly_;};
  bool DoForces() {return Do_Forces;};
  bool UseFiniteDifferenceGradients() {return UseFiniteDifferenceGradients_;};
  bool DoFiniteDifferenceFreqs() {return DoFiniteDifferenceFreqs_;};
  int GetMaxOptCycles() {return MaxOptCycles;};
  bool UseDLFind() {return UseDLFind_;};
  void SetUseDLFind(bool value) {UseDLFind_ = value;};
  bool UseKNitro() {return UseKNitro_;};
  void SetKNitro(bool value) {UseKNitro_ = value;};



  void SetMonomerForFrequency(int imon) {monomer_frequency_ = imon;};
  int GetFrequencyForMonomer() {return monomer_frequency_;};

  //bool DoRamanIntensties() {return do_raman_intensities_;};

  bool IsPeriodic() {return Periodic;};
  bool ReadLatticeVectors() {return read_lattice_vectors;};
  bool ReadFractionalCoordinates() {return read_fractional_coordinates;};
  bool UseCrystalSymmetry() {return crystal_symmetry;};
  bool TinFoilBoundaryConditions() {return tin_foil_boundary_conds;};

  bool UseEmbeddingCharges() {return Use_Embedding_Charges;};
  bool NeglectManyBody() {return Neglect_Many_Body;};
  bool DoCounterpoise() {return Counterpoise;};
  bool DoQMBenchmark() {return DoQMBenchmark_;};
  string GetBasePath() {return base_path;};
  string GetQMPath() {return path_QM;};
  string GetMMPath() {return path_MM;};
  
  string GetQChemRem() {return QC_rem;};
  string GetQChemRem2() {return QC_rem2;};
  string GetQChemBasis() {return QC_basis;};
  string GetTinkerRem() {return Tinker_rem;};
  string GetAIFFRem() {return AIFF_rem;};
  string GetHMBIRem() {return HMBI_rem;};
  
  int GetMMType() {return MM_Type;};

  bool DoLocal2BodyTruncation() {return do_local_two_body_truncation;};
  double GetLocalCutoff(int i) {if (i==0) return cutoff0; else return cutoff1;};
  void SetLocalCutoff(int i, double value) {if (i==0) cutoff0=value; else cutoff1=value;};
  bool ScanLocalTruncationParameters() {return scan_local_cutoffs;};
  double GetMinLocalCutoff() {return min_cutoff1;};
  int GetNumberOfCutoffScanSteps() {return N_cutoff_steps;};

  double GetExpansionFactor() {return ExpansionFactor;};

  string GetAIFFBasisSet() {return AIFFBasisSet;};  // by Ali
  string GetCamCaspHome() {return CamCaspHome;};  // by Ali
  double GetIonizationPotential() {return IonizationPotential;};  // by Ali
  double GetDampingFactor() {return DampingFactor;};  // by Ali
  bool DoAtomicInductionDamping() {return DoAtomicInductionDamping_;};
  double GetAtomicInductionDampingFactor(int atomic_number) {return AtomicInductionDampingFactor[atomic_number];}; // GJB exploratory
  void SetAtomicInductionDampingFactor(int atomic_number, double value) {AtomicInductionDampingFactor[atomic_number]=value;}; // GJB exploratory
  bool DoDispersionDamping() {return DoDispersionDamping_;}

  int GetMaxPolarizationCycles() {return MaxPolarizationCycles;};
  bool OrientDebug() {return OrientDebug_;};
  bool UseGlobalCoordinates() {return UseGlobalCoords_;};
  bool BuildForceFieldOnly() {return BuildForceFieldOnly_;};
  bool RunJobsOnly() {return RunJobsOnly_;};

  // Functions for working with AIFF force field parameters
  void SetMaxPolarizationRadius(double r) {polarization_cutoff = r;};
  double GetMaxPolarizationRadius() {return polarization_cutoff;};

  void SetMaxTwoBodyDispersionRadius(double r) {two_body_dispersion_cutoff = r;};
  double GetMaxTwoBodyDispersionRadius() {return two_body_dispersion_cutoff;};

  void SetMaxThreeBodyDispersionRadius(double r) {three_body_dispersion_cutoff = r;};
  double GetMaxThreeBodyDispersionRadius() {return three_body_dispersion_cutoff;};

  double GetThreeBodyDispersionPrintThreshold() {return three_body_dispersion_print_threshold;};



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

  void SetEwaldInductionType(int value) {Ewald_induction_type = value;};
  double GetEwaldInductionType() {return Ewald_induction_type;}; 

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

  int GetOptType() {return OptType_;};
  int GetLineSearchType() {return LineSearchType_;};
  int GetCoordType() {return CoordType_;};
  void IncrementOptCycle() {OptCycle_++;};
  int GetOptCycle() {return OptCycle_;};

  bool EstimateThreeBodyDispersion() {return EstimateThreeBodyDispersion_;};
  bool DoMP2DispersionCorrection() {return MP2DispersionCorrection_;};

  bool TinkerDebug() {return tinker_debug_;};

  


};

#endif
