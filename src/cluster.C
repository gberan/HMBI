// for getpid:
#include <sys/types.h>
#include <unistd.h>
// other:
#include "cluster.h"
#include "atom.h"  // by Ali
#include <sys/stat.h>
#include "opt.h"
#include "nmr.h" // JDH
#include "ee.h"  // JDH
#include "input.h" // JDH
#include "constants.h"
#include "quasiharmonic.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <vector>
#include <string>
#include <sstream>
#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */

using namespace hmbi_constants;

// Define some tags for parallel job types
#define QCHEM_TAG 1
#define TINKER_TAG 2
#define CAMCASP_TAG 3
#define QCHEMCP_TAG 5
#define MOLPRO_TAG 6
#define CRYSTAL_TAG 7
#define G09_TAG 8
#define GDMA_TAG 9
#define DALTON_TAG 11
#define QUANTUM_ESPRESSO_TAG 12
#define ORCA_TAG 13
#define PSI4_TAG 14
#define DFTB_TAG 15
#define DIETAG 10


Cluster::Cluster() : spins(NULL), charges(NULL), Natoms_per_monomer(NULL), 
		     Monomers(NULL), Dimers(NULL), DimerImages(NULL),
		     AtomicSymbols(NULL), AtomicNumbers(NULL), AtomicMasses(NULL),
		     reciprocal_cell(NULL), unit_cell(NULL), 
		     Grad_key(NULL), InputFile(NULL) {

  
  Energy_HMBI = 0.0;
  Energy_Vib = 0.0;
  Entropy = 0.0;
  Cv = 0.0;
  gruneisen_parameter = 0.0;
  Energy_MM = 0.0;
  NumberOfLinesInInput = 0;

}

void Cluster::Initialize(ifstream &infile, int Nproc, string input_filename) {

  // Print Header
  PrintHeader(); 

  StoreInputFilename(input_filename);
  // Search for job title string
  //Title = GetJobTitle(infile);
  //printf("Job: %s\n",Title.c_str() );
  // Print the input file back out to keep a record of what was run.
  PrintInputFile(infile);
  StoreInputFile(infile);
  Energy_HMBI_old = 0.0;
  

  // Read the HMBI job parameters
  ReadHMBIParameters(infile);  // Note, QC_rem set later on.



  //forces needed for HMBI frequency calculation
  if( Params::Parameters().DoFreq() && !Params::Parameters().CreateJobsOnly() && Params::Parameters().GetMMType() != 5 &&
      !(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable()) ){


    //At this point, if ANALYZE_ONLY = TRUE, assume that both the forces and hessians are availible
    if(!Params::Parameters().RunJobs()){
      Params::Parameters().SetParameter("ARE_FORCES_AVAILABLE","TRUE");
      Params::Parameters().SetParameter("ARE_HESSIANS_AVAILABLE","TRUE");
    }

    //are force calculations already available
    if(Params::Parameters().AreForcesAvailable() ) {
      // no need to perform force calculations...
      // so turn Do_Force flag on and analyze_only = true to read the gradients first
      int maxoptcycles = Params::Parameters().GetMaxOptCycles();
      stringstream ss;//create a stringstream
      ss << maxoptcycles;//add number to the stream
      Params::Parameters().SetParameter("ANALYZE_ONLY","TRUE");
        Params::Parameters().SetParameter("JOBTYPE","FORCE");
        Params::Parameters().SetParameter("DO_FREQ_AFTER_OPT","1");
        //remember that SetParameter("JOBTYPE","FORCE") makes max_opt_cycles =0, so we do the following reset
        //Params::Parameters().SetParameter("MAX_OPT_CYCLES",ss.str() ); //JLM got rid of this line to make it truly an analyze only calculation

      Params::Parameters().SetParameter("INITIALIZE_HESSIANS","1");
    }
    else {
      int maxoptcycles = Params::Parameters().GetMaxOptCycles();
      stringstream ss;//create a stringstream
      ss << maxoptcycles;//add number to the stream
      // now we need to perform force calculations...
      // so turn Do_Force flag on and analyze_only = false to run the gradients' jobs first
      //Params::Parameters().SetParameter("ANALYZE_ONLY","FALSE");
	  Params::Parameters().SetParameter("JOBTYPE","FORCE");
	  Params::Parameters().SetParameter("DO_FREQ_AFTER_OPT","1");
	  //remember that SetParameter("JOBTYPE","FORCE") makes max_opt_cycles =0, so we do the following reset
	  //Params::Parameters().SetParameter("MAX_OPT_CYCLES",ss.str() ); //JLM got rid of this line to make it truly an analyze only calculation

	Params::Parameters().SetParameter("INITIALIZE_HESSIANS","1");
    }

   }
   else if ( Params::Parameters().DoFreq() && Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable() ){
    Params::Parameters().SetParameter("JOBTYPE","ENERGY"); 
  }

  // Read the AIFF parameters  --- by Ali
  ReadAIFFParameters(infile); 

  // Create QM & MM directories, if needed
  CreateWorkingDirectories();
  
  // Print out information about external software packages in use
  PrintExternalPrograms();

  //check that cutoff0 is not less than cutoff1
  if(Params::Parameters().GetLocalCutoff(0) < Params::Parameters().GetLocalCutoff(1)){
    printf("ERROR::Cluster::Initialize(): cutoff0 is less than cutoff1. Alter cutoffs and rerun.\n");
    exit(0);
  }

// Watit - Read CIF File
  if (Params::Parameters().ReadCIFFile()) {
    printf("CIF code not working yet. Exiting HMBI\n");
    exit(0);
    //string cif = ReadCIFSection(infile);
    //printf("CIF_file: %s\n",cif.c_str());
    //CreateInputfromCIF(cif);	
/*	    	if (Params::Parameters().PrintLevel() > 0) {
	    		printf("Read CIF path from $cif section:\n%s",cif.c_str() );	
		}
*/
  }
//

  // Count how many monomers are present
  NMon = FindNumberOfMonomers(infile);
  NDim = NMon*(NMon-1)/2;
  TotalDim_nosym = NDim;//yoni: counter to determine saving with 
  //with symmetry. Will include image-dimers
  UniqueAtoms = 0;//yoni: # of Defining atoms that all atoms relate to by symmetry
  UniqueMon = 0;//yoni: number of monomers under symmetry
  UniqueDim = 0;//yoni: number of dimers under symmetry
  NDim_trunc =0;
  NDimImages_trunc =0;
  printf("The system contains %d monomers, %d dimers\n\n",NMon,NDim);

  //Gruneisen parameters for QuasiHarmonic approximation has not been set
  //Gruneisen_init = 0;

  // Define the conversion factor for hartrees -> kJ/mol/monomer
  HartreeToKJpMM = HartreesToKJpermole/(double) NMon;
  
  // Find number of atoms in each monomer and in the total cluster
  CountNumberOfAtoms(infile);

  //lock lattice parameters under symmetry
  b_locked_by_a=0; 
  c_locked_by_a=0; 
  c_locked_by_b=0; 
  beta_locked_by_alpha=0; 
  gamma_locked_by_alpha=0; 
  gamma_locked_by_beta=0;
  lock_alpha=0;
  lock_beta=0;
  lock_gamma=0;
  
  if(Params::Parameters().GetQMType()==1) {//QChem QM

    // Read in Q-Chem $rem section
    string qc_rem = ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem); 
    //if (Params::Parameters().PrintLevel() > 0) 
    printf("Q-Chem $rem section:\n%s",qc_rem.c_str() );
    
    // Read in Q-Chem $basis section
    string qc_basis = ReadQChemBasisSection(infile);
    Params::Parameters().SetParameter("QC_BASIS",qc_basis);
    if (Params::Parameters().PrintLevel() > 0)
      printf("Q-Chem $basis section:\n%s",qc_basis.c_str() );

    if(Params::Parameters().DoCBS()){
      printf("ERROR::Cluster.Initialize(): CBS extrapolation has not been implemented for Qchem.\n");
      exit(0);
    }

    if(Params::Parameters().Params::Parameters().DoCCSDTCorrection()){
      printf("ERROR::Cluster.Initialize(): CCSD(T) correction has not been implemented for Qchem.\n");
      exit(0);
    }
  }
  else if(Params::Parameters().GetQMType()==2){//Molpro
    string molpro_rem = ReadMolProSection(infile);
    Params::Parameters().SetParameter("MOLPRO_REM",molpro_rem);
    if (Params::Parameters().PrintLevel() > 0)
      printf("MolPro $rem section:\n%s",molpro_rem.c_str() );
    string molpro_inst = ReadMolProInstSection(infile);
    if (Params::Parameters().PrintLevel() > 0)
      printf("MolPro $inst section:\n%s",molpro_inst.c_str() );
    Params::Parameters().SetParameter("MOLPRO_INST",molpro_inst);

    //CBS limit
    if(Params::Parameters().DoCBS()){
      string molpro_basis_CBS = ReadMolProCBSBasisSection(infile);
      if (Params::Parameters().PrintLevel() > 0)
	printf("MolPro CBS $basis section:\n%s",molpro_basis_CBS.c_str() );
      Params::Parameters().SetParameter("MOLPRO_CBS",molpro_basis_CBS);
      string molpro_inst_HF = ReadMolProInstHFSection(infile);
       if (Params::Parameters().PrintLevel() > 0)
	printf("MolPro HF $inst section:\n%s",molpro_inst_HF.c_str() ); 
       Params::Parameters().SetParameter("MOLPRO_HF",molpro_inst_HF);
    }

    //CCSD(T) correction
    if(Params::Parameters().DoCCSDTCorrection()){
      string molpro_basis_CCSDT = ReadMolProCCSDTBasisSection(infile);
      if (Params::Parameters().PrintLevel() > 0)
	printf("MolPro CCSD(T) $basis section:\n%s",molpro_basis_CCSDT.c_str() );
      Params::Parameters().SetParameter("MOLPRO_BASIS_CCSDT",molpro_basis_CCSDT);
      string molpro_CCSDT_MP2 = ReadMolProCCSDTMP2Section(infile);
      if (Params::Parameters().PrintLevel() > 0)
	printf("MolPro CCSD(T) MP2 section:\n%s",molpro_CCSDT_MP2.c_str() );
      Params::Parameters().SetParameter("MOLPRO_CCSDT_MP2",molpro_CCSDT_MP2);
      string molpro_CCSDT_inst = ReadMolProCCSDTInstSection(infile);
      if (Params::Parameters().PrintLevel() > 0)
	printf("MolPro CCSD(T) Inst section:\n%s",molpro_CCSDT_inst.c_str() );
      Params::Parameters().SetParameter("MOLPRO_CCSDT_INST",molpro_CCSDT_inst);
      
    }

  }
  else if ( Params::Parameters().GetQMType()==3){//G09
    string g09 = ReadG09Section(infile);
    Params::Parameters().SetParameter("G09",g09);
  }
  else if ( Params::Parameters().GetQMType()==4){//DALTON
    string dalton = ReadDaltonSection(infile);
    Params::Parameters().SetParameter("DALTON",dalton);
    if ( Params::Parameters().UseElectrostaticEmbedding() ) {
      string orca = ReadOrcaSection(infile);
      Params::Parameters().SetParameter("ORCA",orca);
    }
  }
  else if ( Params::Parameters().GetQMType()==5){//Quantum Espresso
    string qe_species = ReadQuantumEspressoSpeciesSection(infile);
    string qe_supercell = ReadQuantumEspressoSupercellSection(infile);
    string kpoints = ReadKPoints(infile);
    if(!qe_species.empty()) Params::Parameters().SetParameter("QE_SPECIES",qe_species);
    if(!qe_supercell.empty()) Params::Parameters().SetParameter("QE_SUPERCELL",qe_supercell);
    if(!kpoints.empty()) Params::Parameters().SetParameter("KPOINTS",kpoints);
    //printf("Checking to make sure it set the parameter correctly: \n%s\n", Params::Parameters().GetQuantumEspressoControlSection().c_str() );
    //fflush(stdout);
    //printf("For comparison the Quantum Espresso section:\n%s",qe.c_str() );
    //fflush(stdout);
  }
  else if ( Params::Parameters().GetQMType()==6){//ORCA
    string orca = ReadOrcaSection(infile);
    Params::Parameters().SetParameter("ORCA",orca);
  }
  else if ( Params::Parameters().GetQMType()==7) { //PSI4
    //printf("DEBUG Read in Psi4 Section...\n"); fflush(stdout);
    string psi4 = ReadPSI4Section(infile);
    Params::Parameters().SetParameter("PSI4",psi4);
    string PSI4_rem = ReadPSI4RemSection(infile);
    Params::Parameters().SetParameter("PSI4_REM",PSI4_rem);
    if (Params::Parameters().PrintLevel() > 0)
      printf("PSI4 $rem section:\n%s",PSI4_rem.c_str() );
  }
  else if ( Params::Parameters().GetQMType()==8){//DFTB+
    //Note: Re-using a lot of the QE input sections. Consider making a more general name for these functions //JLM
    string qe_supercell = ReadQuantumEspressoSupercellSection(infile);
    string kpoints = ReadKPoints(infile);
    string slater_koster = ReadSlaterKoster(infile);
    string hubbard = ReadHubbardDerivatives(infile);
    string max_angular_momentum = ReadMaxAngularMomentum(infile);
    if(!qe_supercell.empty()) Params::Parameters().SetParameter("QE_SUPERCELL",qe_supercell);
    if(!kpoints.empty()) Params::Parameters().SetParameter("KPOINTS",kpoints);
    if(!slater_koster.empty()) Params::Parameters().SetParameter("SLATER_KOSTER",slater_koster);
    if(!hubbard.empty()) Params::Parameters().SetParameter("HUBBARD_DERIVS",hubbard);
    if(!max_angular_momentum.empty()) Params::Parameters().SetParameter("MAX_ANGULAR_MOMENTUM",max_angular_momentum);
    //printf("Checking to make sure it set the parameter correctly: \n%s\n", Params::Parameters().GetQuantumEspressoControlSection().c_str() );
    //fflush(stdout);
    //printf("For comparison the Quantum Espresso section:\n%s",qe.c_str() );
    //fflush(stdout);
  }
  else{
    printf("Cluster::Cluster() - Unknown QM Type = %d\n",
	   Params::Parameters().GetQMType());
    exit(1);

  }

  if (Params::Parameters().CustomBasis() != 0 ) {
    // Read in all of the custom basis info
    string basis = ReadHBasis(infile);
    Params::Parameters().SetParameter("HBasis",basis);
    
    
    basis = ReadCBasis(infile);
    Params::Parameters().SetParameter("CBasis",basis);
        
    
    basis = ReadNBasis(infile);
    Params::Parameters().SetParameter("NBasis",basis);

    
    basis = ReadOBasis(infile);
    Params::Parameters().SetParameter("OBasis",basis);

    
    basis = ReadSBasis(infile);
    Params::Parameters().SetParameter("SBasis",basis);

    
    basis = ReadClBasis(infile);
    Params::Parameters().SetParameter("ClBasis",basis);

    
    basis = ReadIBasis(infile);
    Params::Parameters().SetParameter("IBasis",basis);

    
    basis = ReadSnBasis(infile);
    Params::Parameters().SetParameter("SnBasis",basis);

    
    basis = ReadPBasis(infile);
    Params::Parameters().SetParameter("PBasis",basis);

    basis = ReadKBasis(infile);
    Params::Parameters().SetParameter("KBasis",basis);

    basis = ReadNaBasis(infile);
    Params::Parameters().SetParameter("NaBasis",basis);

    basis = ReadBrBasis(infile);
    Params::Parameters().SetParameter("BrBasis",basis);

    basis = ReadFBasis(infile);
    Params::Parameters().SetParameter("FBasis",basis);

    
    basis = ReadWBasis(infile);
    Params::Parameters().SetParameter("WBasis",basis);

    basis = ReadVBasis(infile);
    Params::Parameters().SetParameter("VBasis",basis);
    
    printf("All specified custom basis sets read in \n"); fflush(stdout);
  }



  if (Params::Parameters().GetMMType()==1) { // Tinker MM
    // Read in Tinker rem section 
    printf("Using Tinker to compute many-body interactions\n");
    string tinker_rem = ReadTinkerRemSection(infile);
    Params::Parameters().SetParameter("TINKER_REM",tinker_rem); 
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Tinker $rem section:\n%s",tinker_rem.c_str() );
  }
  else if (Params::Parameters().GetMMType()==2) { // AIFF
    // Read in AIFF rem section  --- by Ali
    printf("Using an ab initio force field to compute many-body interactions\n");
    string aiff_rem = ReadAIFFRemSection(infile);
    Params::Parameters().SetParameter("AIFF_REM",aiff_rem); 
    if (Params::Parameters().PrintLevel() > 0) 
      printf("AIFF keyword section:\n%s",aiff_rem.c_str() );
  }
  else if (Params::Parameters().GetMMType()==3) { // QChem for "MM"
    // Read in Secondary QChem rem section
    printf("Using Q-Chem to compute many-body interactions\n");
    string qchem2_rem = ReadQChemRemSection(infile,2);
    Params::Parameters().SetParameter("QC_REM2",qchem2_rem); 
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Secondary Q-Chem $rem section:\n%s",qchem2_rem.c_str() );
  }
  else if (Params::Parameters().GetMMType()==4) { // EE-PA
    printf("Performing an EE-PA calculation.\n");
    Params::Parameters().SetParameter("NEGLECT_MANY_BODY","true");
  }
  else if (Params::Parameters().GetMMType()==5) {
    printf("Using Crystal for calucations\n");
    //heading
    string crystal_heading = ReadCrystalHeadingSection(infile);
    Params::Parameters().SetParameter("CRYSTAL_HEADING",crystal_heading);

    //basis
    string crystal_basis = ReadCrystalBasisSection(infile);
    Params::Parameters().SetParameter("CRYSTAL_BASIS",crystal_basis);

    //ending
    string crystal_ending = ReadCrystalEndingSection(infile);
    Params::Parameters().SetParameter("CRYSTAL_ENDING",crystal_ending);
  }
  else if ( Params::Parameters().GetMMType()==99) {
    printf("Using GDMA for calculations\n");
    string g09 = ReadG09Section(infile);
    Params::Parameters().SetParameter("G09",g09);
  }
  else if (Params::Parameters().GetMMType()==55) { // GDMA
    if ( Params::Parameters().GetQMPackage() != "G09" ) {
      string g09 = ReadG09Section(infile);
      Params::Parameters().SetParameter("G09",g09);
    }
  }
  else if ( Params::Parameters().GetMMType()==98) {
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      printf("Using G09 - CHelpG charges for embedding environment\n"); fflush(stdout);
      string g09 = ReadG09Section(infile);
      Params::Parameters().SetParameter("G09",g09);
    } else {
      printf("Using ORCA to calculate chelpG charges for the monomers\n");
      string orca = ReadOrcaSection(infile);
      Params::Parameters().SetParameter("ORCA",orca);
    }
  }
  else if ( Params::Parameters().GetMMType()==97) {
    printf("Using G09 - Hirshfeld charges for embedding environment\n"); fflush(stdout);
    if ( Params::Parameters().GetQMPackage() != "G09" ) {
      string g09 = ReadG09Section(infile);
      Params::Parameters().SetParameter("G09",g09);
    }
  }
  else {
    printf("Cluster::Cluster() - Unknown MM Type = %d\n",
	   Params::Parameters().GetMMType());
    exit(1);
    
  }

  // Read in Monomer geometries, and initialize the monomers
  ReadGeometry(infile);
  //Params::Parameters().LocalCoordinateAxesInitialized(true);
  //printf("Local axes defined: %d\n",Params::Parameters().LocalCoordinateAxesInitialized());

  if ( Params::Parameters().GetJobTypeStr() == "expand" ) {
    AdjustIntermolecularSpacing( Params::Parameters().GetExpansionFactor());
    printf("Exiting.\n");
    exit(1);
  }

  // Create a list of the Atomic symbols
  SetAtomicSymbols();
  

  // Read current coordinates into the vector AtomicCoordinates
  // This gets used if we do geometry optimizations, for example
  // Initialize the vector first
  AtomicCoordinates.Initialize( 3 * GetTotalNumberOfAtoms() ); 
  ReadCurrentCoordinates();


  //Initialize the frequency vector
  if(Params::Parameters().DoFreq() || Params::Parameters().DoFiniteDifferenceFreqs() ||
     Params::Parameters().DoFreqAfterOpt() || Params::Parameters().DoQuasiHarmonic()){
    Freqs.Initialize(3*GetTotalNumberOfAtoms());
    IntIR.Initialize(3*GetTotalNumberOfAtoms());
    IntRaman.Initialize(3*GetTotalNumberOfAtoms());
  }

  //QuasiHarmonic Approximation
  //Gruneisen_init = 0;
  /*
  if(Params::Parameters().DoQuasiHarmonic()){
    reference_volume = 0.0;
    reference_frequencies.Initialize(3*GetTotalNumberOfAtoms());
    Gruneisen_parameters.Initialize(3*GetTotalNumberOfAtoms());
  }
  */


  //Finding symmetry for monomers
  for(int i=1;i<=NMon;i++){
    if(!Params::Parameters().UseSpaceSymmetry()){
      UniqueMon++;
      UniqueAtoms += Monomers[i].GetNumberOfAtoms();
    }else if(!Monomers[i].SymmetryCheck(Monomers)){
      UniqueMon++;
      //UniqueAtoms += Monomers[i].GetNumberOfAtoms();
       UniqueAtoms += Monomers[i].GetNumberOfUniqueAtoms();
    }
  }

  //if none of the monomers are symmetical. Turn off space symmetry
  if(UniqueMon == NMon && Params::Parameters().UseSpaceSymmetry()){
    printf("No symmetrical monomers. Turning off symmetry.\n");
    Params::Parameters().SetParameter("SPACE_SYMMETRY","FALSE");
    for(int i=1;i<=NMon;i++){
      Monomers[i].ResetSymmetry();
    }
    UniqueAtoms = GetTotalNumberOfAtoms();
  }
  
  //Reads current coordinates for the Symmetry Unique Atoms
  //Creates a reduced gradient for the geometry optimizations.
  //Using symmetry, the opimization for the entire system can be determined.
  SymmetryUniqueCoordinates.Initialize(3 * UniqueAtoms);
  ReadSymmetryUniqueCoordinates();

  if(Params::Parameters().UseSpaceSymmetry() && Params::Parameters().PrintSymmetryInfo()) {
    //gives the symmetry factor
    for (int i=1;i<=NMon;i++) {
      printf("The symmetry factor of m%i is %i\n",
	     Monomers[i].GetIndex(),Monomers[i].GetSymmetryFactor());

    }
    printf("Number of unique atoms is %i\n", UniqueAtoms);
    printf("Number of unique monomers is %i\n", UniqueMon);
  }

  // Generate the dimers.  Note, counting starts at 1;
  Dimers = new Dimer[NDim+1];
  int index = 1;
  for (int a=1;a<NMon;a++){
    for (int b=a+1;b<=NMon;b++) {
      //printf("index = %d, d(%d,%d)\n",index,a,b);
      Dimers[index].Initialize(Monomers[a], Monomers[b]);
      if(!Params::Parameters().UseSpaceSymmetry())
	UniqueDim++;
      else if(!Dimers[index].SymmetryCheck(Dimers, index))
	UniqueDim++;
      index++;
    }
  }


  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the monomers and the dimers
    printf("*** Print out the monomers\n");
    for (int i=1;i<=NMon;i++) {
      Monomers[i].PrintAll();
      printf("The number of unique monomers is %i\n",UniqueMon);
      printf("The symfac Monomer %i is %i\n",i,Monomers[i].GetSymmetryFactor());
    }
    printf("*** Print out the dimers\n");
    for (int i=1;i<=NDim;i++) {
      Dimers[i].PrintAll();
      printf("The symmetry factor of d(%i,%i) is %i\n", 
	     Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),Dimers[i].GetSymmetryFactor());
    }
  }

  // Read in the dispersion atom types if doing 3-body ATM dispersion
  if ( Params::Parameters().EstimateThreeBodyDispersion() ){
    ReadDispersionAtomTypes(infile);
  }
  if ( Params::Parameters().DoMP2DispersionCorrection() ) {
    ReadDispersionCoefficients(infile);
  }
 
  // If doing full QM we need unit cell parameters but no symmetrical dimers 
  if ( Params::Parameters().UseFullQMOnly() ) {
    if ( Params::Parameters().ReadLatticeVectors() ) //JLM
      OldReadUnitCellDimensions(infile);
    else{
      printf("Reading cell dimensions\n");
      ReadUnitCellDimensions(infile);
    }

    //Setting Fractional Coordinates for each of the atoms in the unit cell
    SetFractionalCoordinates();
    ReadSymmetryUniqueFractionalCoordinates();
    DetermineMonomerFractionalSymmetry();
  }

  // Handle periodic boundary conditions (PBC)
  if ( Params::Parameters().IsPeriodic() && !Params::Parameters().UseFullQMOnly() ) {
    // Read in the unit cell parameters
    fflush(stdout); // flush the output stream*/
    if ( Params::Parameters().ReadLatticeVectors() )
      OldReadUnitCellDimensions(infile);
    else{
      printf("Reading cell dimensions");
      ReadUnitCellDimensions(infile);
    }

    //Setting Fractional Coordinates for each of the atoms in the unit cell
    SetFractionalCoordinates();


    ReadSymmetryUniqueFractionalCoordinates();



    // Turn on local truncation
    Params::Parameters().SetParameter("LOCAL_2_BODY","TRUE");

    // Identify real space periodic images
    CreatePeriodicImageDimerList();

    //Sometimes, monomer and dimer rotation are created incorrectly
    //is has been noticed in linear molecules
    if(Params::Parameters().UseSpaceSymmetry()){

      //Check dimers first
      /*for(int i=1;i<=NDim;i++){
	for(int j=0;j<Dimers[i].GetSymmetryFactor();j++){
	  CorrectingDimerRotations(i,0,0,j);   
	} 
	for(int j=0;j<Dimers[i].GetPeriodicSymmetryFactor();j++){
	  CorrectingDimerRotations(i,0,1,j);   
	}
      }

      //Check image dimers
      for(int i=1;i<=NDim_images;i++){
	for(int j=0;j<DimerImages[i].GetSymmetryFactor();j++){
	  CorrectingDimerRotations(i,1,0,j); 
	}
      }*/


      //check monomers
      for(int i = 1; i<=NMon;i++)
	if((Monomers[i].GetSymmetryFactor() == 0))
	  CorrectingMonomerRotations(i);
    }


    //Setting Fractional Coordinates for each monomer
    DetermineMonomerFractionalSymmetry();

  }
  else
    NDim_images = 0;

  //Initializing Gradient if force calculations are being performed
  if(Params::Parameters().DoForces() || Params::Parameters().DoFreq()){
    if(Params::Parameters().UseFullQMOnly())
      Grad_HMBI.Initialize(3*UniqueAtoms+9);
    else if( (Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() == 1
	&& !(Params::Parameters().DoFiniteDifferenceFreqs())))
      Grad_HMBI.Initialize(3*UniqueAtoms+6);
    else
      Grad_HMBI.Initialize(3*UniqueAtoms);  
  }

    
  //for (int i=1;i<=NMon;i++) {//yoni:gives symmetry matrix
  //  printf(" m%i\n",i);
  //  Monomers[i].GetRotationMatrix().Print("Rotation");
 
  //}

  //printing savings due to symmetry
  if(Params::Parameters().UseSpaceSymmetry()){
    
    printf("\n");
    printf("Applying space group symmetry to reduce the numbers of monomer/dimer calculations required\n");
    printf("There are %i dimers with symmetry and %i without. (Approximate savings factor = %.1f)\n",
	   UniqueDim+NDim_images,TotalDim_nosym, (double) TotalDim_nosym/(UniqueDim+NDim_images));
    printf("\n");
    
    if(Params::Parameters().PrintSymmetryInfo()){
      for (int i=1;i<=NDim;i++) {
	printf("The symmetry factor and periodic symmetry factor of d(%i,%i) is %i and  %i\n", 
	       Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),Dimers[i].GetSymmetryFactor(),
	       Dimers[i].GetPeriodicSymmetryFactor());
	//if(Dimers[i].GetSymmetryFactor()){
	//printf("Symmetrical to ");
	//for(int j=0;j<Dimers[i].GetSymmetryList().size()/2;j++)
	//printf("d(%i,%i) ",
	//Dimers[i].GetSymmetryList()[2*j],Dimers[i].GetSymmetryList()[2*j+1]);
	// printf("\n");
	//for(int j=0;j<Dimers[i].GetAtomEquivalency().size();j++)
	// 	Dimers[i].GetAtomEquivalency()[j].Print("Atom Equivalency");
	//for(int j=0;j<Dimers[i].GetRotationList().size();j++)
	//	Dimers[i].GetRotationList()[j].Print("Dimer Rotation");
	//printf("and Symmetrical to Image Dimers ");
	//for(int j=0;j<Dimers[i].GetPeriodicSymmetryList().size()/2;j++)
	//	printf("d(%i,%i) ",
	//       Dimers[i].GetPeriodicSymmetryList()[2*j],Dimers[i].GetPeriodicSymmetryList()[2*j+1]);
	//printf("\n");
	//for(int j=0;j<Dimers[i].GetPeriodicAtomEquivalency().size();j++)
	//Dimers[i].GetPeriodicAtomEquivalency()[j].Print("Atom Equivalency");
	
	//for(int j=0;j<Dimers[i].GetPeriodicRotationList().size();j++)
	//  Dimers[i].GetPeriodicRotationList()[j].Print("Dimer Periodic  Rotation");
	// printf("\n");
	// printf("k_list = ");
	// for(int j=0;j<Dimers[i].GetSymmetricalImageCell().size();j++)
	//	 cout << Dimers[i].GetSymmetricalImageCell()[j] << " ";
	// printf("\n");
	//}
      } 
      for (int i=1;i<=NDim_images;i++) {  
	printf("The symmetry factor of d(%i,%i) is %i\n", 
	       DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB(),DimerImages[i].GetSymmetryFactor());
	//printf("Symmetrical to ");
	//for(int j=0;j<DimerImages[i].GetSymmetryList().size()/2;j++)
	//printf("d(%i,%i) ",
	//       DimerImages[i].GetSymmetryList()[2*j],DimerImages[i].GetSymmetryList()[2*j+1]);
	//printf("\n");
	//for(int j=0;j<DimerImages[i].GetAtomEquivalency().size();j++)
	//	  DimerImages[i].GetAtomEquivalency()[j].Print("Atom Equivalency");
	
	//for(int j=0;j<DimerImages[i].GetRotationList().size();j++)
	//  DimerImages[i].GetRotationList()[j].Print("Dimer Rotation");
	//printf("\n");
	//printf("k_list = ");
	//for(int j=0;j<DimerImages[i].GetSymmetricalImageCell().size();j++)
	//  cout << DimerImages[i].GetSymmetricalImageCell()[j] << " " ;
	//printf("\n");
      }
    }
  }
  fflush(stdout); // flush the output stream

  //exit(0);
  // Initialize Gradient memory, if needed
  if ( Params::Parameters().DoForces() ) {
    if ( !Params::Parameters().IsPeriodic() 
	 /*&& Params::Parameters().DoFiniteDifferenceFreqs()*/) {
      Grad_MM.Initialize(3*UniqueAtoms);
	//Grad_MM.Initialize( 3*GetTotalNumberOfAtoms );
    }
      else if ( Params::Parameters().IsPeriodic() ) {
	if ( (Params::Parameters().GetMMType()==2) ||  Params::Parameters().DoFiniteDifferenceFreqs()) {  // AIFF
	  Grad_MM.Initialize(3*UniqueAtoms );                
	} 
	else        
	  Grad_MM.Initialize( 3*UniqueAtoms + 6);
      }

    // Create a list of the starting index for each monomer in the full
    // gradient.  Like the Monomers, the key indexes count from 1->NMon.
    // On the other hand, in the gradient, the indexing starts at 0.
    Grad_key= new int[NMon+1];
    Grad_key[1] = 0; // set the first one by hand
    for (int i=2;i<=NMon;i++) {
      Grad_key[i] = Grad_key[i-1] + 3*Monomers[i-1].GetNumberOfAtoms();
    }
  }

  if ( Params::Parameters().UseFullQMOnly()) { //JLM
    Grad_QM.Initialize(3*UniqueAtoms+9);
    //Grad_QM.Initialize(3*GetTotalNumberOfAtoms()+9);
  }

  if ( Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())
    WriteCrystalXYZFile(2,2,2);

  // Done with all the initialization.  
  // Run RunJobsAndComputeEnergy() to start doing the hard work.

  //Change the volume isotropically if requested.
  if(Params::Parameters().ChangeVolume() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
    ChangeCellVolume();
    //ChangeCellVolumeIsotropically();
  }

  AcceptedAtomicCoordinates = SymmetryUniqueCoordinates;
  AcceptedUnitCellAxes = UnitCellAxes;
  AcceptedUnitCellAngles = UnitCellAngles;

}

// Destructor
Cluster::~Cluster() {

  delete [] Monomers;
  delete [] Dimers;
  
  delete [] spins;
  delete [] charges;
  delete [] Natoms_per_monomer;
  delete [] types;

  delete [] AtomicSymbols;
  delete [] AtomicNumbers;
  delete [] AtomicMasses;
  
  if ( Params::Parameters().IsPeriodic()) {
    delete [] MonomerImages;
    delete [] DimerImages;
    delete [] unit_cell;
    delete [] reciprocal_cell;
  }

  if ( Params::Parameters().UseFullQMOnly() && !Params::Parameters().IsPeriodic() ) {
    delete [] unit_cell;
    delete [] reciprocal_cell;
  }

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    delete [] Grad_key;
  }

  if ( Params::Parameters().EstimateThreeBodyDispersion() ) {
    delete [] DispersionAtomTypes;
  }


}


void Cluster::StoreInputFilename(string input_filename) {
  InputFilename = input_filename;
  string param = "INPUT_FILENAME";
  Params::Parameters().SetParameter(param,InputFilename);
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//////                     BEGIN JDH MAGNETIC PROPERTIES CODE
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////




void Cluster::PredictMagneticProperties() {

  printf("\n\n\n");
  printf("***************************************************\n");
  printf("*   Fragment-Based Magnetic Propery Prediction    *\n");
  printf("***************************************************\n\n");


  // ELECTROSTATIC EMBEDDING
  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    
    printf("\n\nUsing Electrostatic Embedding... \n");

    // BUILD JOBS
    if ( Params::Parameters().GetMMType()==99 ) {
      for (int i=1;i<=NMon;i++) {
	if (Params::Parameters().RunJobs() || Params::Parameters().BuildForceFieldOnly() ) {
	  Monomers[i].CreateGDMAJob(); 
	}
      }
    } else if ( Params::Parameters().GetMMType()==98 ) {
      for (int i=1; i<=NMon;i++) {
	if (Params::Parameters().RunJobs() || Params::Parameters().BuildForceFieldOnly() ) {
	  if ( Params::Parameters().GetQMPackage() == "G09" ) {
	    Monomers[i].CreateG09ChelpGJob();  
	  } else {
	    printf("Using Orca to compute chelpG charges for monomer %d\n",i);
	    Monomers[i].CreateOrcaChelpGJob();
	  }
	}
      }
    } else if ( Params::Parameters().GetMMType()==97 ) {
      for (int i=1; i<=NMon;i++) {
	if (Params::Parameters().RunJobs() || Params::Parameters().BuildForceFieldOnly() ) {
	  if ( Params::Parameters().GetQMPackage() == "G09" ) {
	    Monomers[i].CreateG09HirshfeldJob();  
	  } else {
	    printf("ERROR: cannot create HirshfeldJob for requested QM package\n");
	    exit(1);
	  }

	}
      }
    } else {
      printf("The NMR code only uses GDMA or ChelpG to save a lot of time! \n");
      exit(1);
    }

    // RUN JOBS
    if ( Params::Parameters().RunJobs() && !Params::Parameters().SeedSCE() ) {
      RunElectrostaticEmbeddingJobs(); // ee.C
    }

    // READ EE DATA
    if ( Params::Parameters().GetMMType()==99 ) {
      for (int imon=1;imon<=NMon;imon++) {
	Monomers[imon].ReadMultipoleMoments();
      }
    } else if (Params::Parameters().GetMMType()==98 ) {
      for (int imon=1;imon<=NMon;imon++) {
	if ( Params::Parameters().GetQMPackage() == "G09" ) {
	  Monomers[imon].ReadG09ChelpGCharges();
	} else if ( Params::Parameters().GetQMPackage() == "ORCA" || Params::Parameters().GetQMPackage() == "DALTON" ) {
	  Monomers[imon].ReadOrcaChelpGCharges();
	} else {
	  printf("ChelpG with %s isn't ready or isn't possible\n", Params::Parameters().GetQMPackage().c_str() );
	  exit(1);
	}
      }
    } else if (Params::Parameters().GetMMType()==97 ) {
      for (int imon=1;imon<=NMon;imon++) {
	Monomers[imon].ReadHirshfeldCharges();
      }
    }
    

    // BETA: Explore two-body contributions to charge embedding
    if ( Params::Parameters().UseTwoBodyCharge() ) {
      printf("Two Body charge calculations simply don't work. \n For example, consider a Cl- amino acid dimer. The dimer calculation can place up to 25 percent of the charge -1 charge on Cl- on the amino acid, which contributes an order of magnitude more two-body charge contribution than the rest of the dimers and if there are 3 amino acids close to the Cl- then the Cl- ends up wth a positive charge in the two-body calcualtion. \n");
      exit(1);
      
      // CREATE IMAGE MONOMER LIST
      // delete [] MonomerImages; // Clear this out just in case
      // delete [] DimerImages;
      // double MAX_CUTOFF = max(
      // 			      Params::Parameters().GetTwoBodyCutoff(),
      // 			      Params::Parameters().GetElectrostaticEmbeddingCutoff()
      // 			      );

      // // Create a new image monomer & dimer list (shouldn't be needed)
      // CreatePeriodicImageMonomerList(MAX_CUTOFF);

      // // assign basis sets using all monomers in the unit cell
      // AssignBasisSetsToAtomsTwoBodyCharge( NMon, Monomers, NMon_images, MonomerImages);
      // // assign monomer lists
      // AssignMonomerListsForTwoBodyElectrostaticEmbedding(NMon, Monomers, NMon_images, MonomerImages);
      
      // CreatePeriodicImageDimerList(Params::Parameters().GetTwoBodyCutoff());
      
      

      // // assign dimer lists
      // AssignDimerListsForFragmentTwoBodyChargeCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );

      

      
      // if (Params::Parameters().IsPeriodic() ) {
      // 	DetermineTwoBodyElectrostaticEmbeddingEnvironment(NMon, Monomers, NDim, Dimers, NMon_images, MonomerImages, NDim_images, DimerImages );
      // } else {
      // 	DetermineTwoBodyElectrostaticEmbeddingEnvironment(NMon, Monomers, NDim, Dimers );
      // }
      // //DetermineTwoBodyElectrostaticEmbeddingEnvironment();

      // // Now put the assignments back
      // delete [] MonomerImages; // Clear this out just in case
      // delete [] DimerImages;
      // MAX_CUTOFF = max(
      // 			      Params::Parameters().GetTwoBodyCutoff(),
      // 			      Params::Parameters().GetElectrostaticEmbeddingCutoff()
      // 			      );
      // CreatePeriodicImageMonomerList(MAX_CUTOFF);
      // AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers, NMon_images, MonomerImages);
      // AssignBasisSetsToAtoms( NMon, Monomers, NMon_images, MonomerImages);
      // CreatePeriodicImageDimerList(Params::Parameters().GetTwoBodyCutoff());
      // AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );
    }

    
    // Self-Consistent EE (optional)
    if ( Params::Parameters().UseSelfConsistentEmbedding() ) {
      if ( Params::Parameters().RunJobs() ) { 
	DetermineSelfConsistentElectrostaticEmbeddingEnvironment();
      }
    }

    

    
    // Electrostatic Embedding Environment Complete
    if ( Params::Parameters().BuildForceFieldOnly() && !Params::Parameters().CreateJobsOnly() ) {
      printf("Electrostatic Embedding Environment Built.  Exiting.\n");
      exit(0);
    } 
    
    // Clean up the directories
    if ( Params::Parameters().GetMMType()==99 ) {
      // Clean up the directories
      string cmd = "cd " + Params::Parameters().GetQMPath();
      cmd += "; rm -f *gdma*; rm -f *.mom tmp*";
      system(cmd.c_str());
    }

  } // END Electrostatic Embedding:

  if ( Params::Parameters().CreateJobsOnly() || Params::Parameters().RunJobs() ) {
    CreateFragmentJobs();
  } else {
    // Don't make the jobs or run the ewald, still need to assign Monomers
    //if (Params::Parameters().IsPeriodic() ) {
    AssignMonomerLists();
      //}
  }

  if ( Params::Parameters().RunJobs() && !Params::Parameters().CreateJobsOnly() ) {
    if (Params::Parameters().IsPeriodic() ) {
      RunFragmentJobs(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );
    } else {
      RunFragmentJobs(NMon, Monomers, NDim, Dimers );
    }
  }


    /* Let's do a quick clean up of the GMDA stuff, move this later
   once I decide on a charge embedding protocol*/
  //if ( !Params::Parameters().KeepNMRFiles()  ) {
  string cmd = "cd " + Params::Parameters().GetQMPath();
  cmd += "; rm -f *gdma*; rm -f *.mom tmp*";
  system(cmd.c_str());
    //}

  if ( Params::Parameters().CreateJobsOnly() ) {
    printf("Fragment jobs created, exiting... \n");
    exit(1);
  }


  if ( Params::Parameters().CalculateEFG() ) {
    printf("\n\n");
    printf("Calculating the EFG tensors and quadrupole coupling constant Cq\n");
    if (Params::Parameters().IsPeriodic() ) {
      ReadEFGData(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);
      ComputeEFGTensors(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);
    } else {
      ReadEFGData(NMon, Monomers, NDim, Dimers );
      ComputeEFGTensors(NMon, Monomers, NDim, Dimers);
    }
  }

  if (Params::Parameters().IsPeriodic() ) {
    ReadNMRData(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages); 
    ComputeNMRShieldingTensors(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);
  } else {
    printf("Reading NMR Data in cluster.C \n"); fflush(stdout);
    ReadNMRData(NMon, Monomers, NDim, Dimers );
    ComputeNMRShieldingTensors(NMon, Monomers, NDim, Dimers );
  }




}


void Cluster::DetermineSelfConsistentElectrostaticEmbeddingEnvironment() {

  printf("\n\nDetermining Self Consistent Electrostatic Embedding Environment\n"); fflush(stdout);

  if ( !(Params::Parameters().GetMMType()==98 || Params::Parameters().GetMMType()==97)  ) {
    printf("ERROR: Self Consistent Electrostatic Embedding code only works with G09/ChelpG or Hirshfeld Charges\n");
  }

  //Get Number of Atoms in Unit Cell
  int natoms_unitcell;
  natoms_unitcell = 0;
  for(int i=1; i<=NMon; i++) {
    natoms_unitcell += Monomers[i].GetNumberOfAtoms();
  }


  // Save Current Atomic Charges into array of doubles
  double *old_charges = new double[natoms_unitcell];
  int old_charges_index;
  old_charges_index = 0;
  //Loop through all monomers
  for(int j=1; j<=NMon; j++) {
    //Loop through all atom
    for(int k=0; k<Monomers[j].GetNumberOfAtoms(); k++) {
      old_charges[old_charges_index] = Cluster::cluster().GetMonomer(j).GetAtom(k).GetMultipoleMoments().GetMoments().Element(0);
      old_charges_index++;
    }
  }

  //First Pass for Converging Charges
  ComputeSelfConsistentEmbeddingCharges();

  //Save New Atomic Charges into array of doubles  
  double *new_charges = new double[natoms_unitcell];
  int new_charges_index;
  new_charges_index = 0;   
  //Loop through all monomers
  //cout << "NEW CHARGES \n";
  for(int j=1; j<=NMon; j++) {
    //Loop through all atoms
    for(int k=0; k<Monomers[j].GetNumberOfAtoms(); k++) {
      new_charges[new_charges_index] = Monomers[j].GetAtom(k).GetMultipoleMoments().GetMoments().Element(0);
      //new_charges[new_charges_index] = Cluster::cluster().GetMonomer(j).GetAtom(k).GetMultipoleMoments().GetMoments().Element(0);
      new_charges_index++;
    }
  }

  //Create Two Vectors of Old and New Charges
  Vector new_charges_vector(new_charges, natoms_unitcell);
  Vector old_charges_vector(old_charges, natoms_unitcell);

  //Create Vector containing the differences between each old charge and new charge    
  Vector charge_difference(new_charges, natoms_unitcell);
  charge_difference.operator-=(old_charges_vector);

  double max_abs_deviation;
  max_abs_deviation = fabs(charge_difference.Max(true));
  
  double rms_deviation;
  rms_deviation = charge_difference.RMS();

  //Count to keep track of number of iterations for convergence
  int count;
  count = 1;

  printf("\n\nCalculating the Self Consistent Electrostatic Embedding Environment: \n");
  printf("Iter. \t  Max(|diff|)  \t RMS \n");
  printf("%d    \t  %f           \t %f  \n", 0, max_abs_deviation, rms_deviation);

  // Make HMBI input parameters for both of these
  double rms_thresh = 0.00005; // too stringent
  double max_thresh = 0.00001; // too stringent
  int nBadSteps = 0;
  while( (rms_deviation > Params::Parameters().GetSelfConsistentRMSThreshold()  || max_abs_deviation > Params::Parameters().GetSelfConsistentMaxThreshold() )  && count < 30  ) {
       
    //Place "new charges" into old charges vector
    for(int i=0; i<natoms_unitcell; i++) {
      old_charges_vector[i] = new_charges_vector.Element(i);
    }
  
    ComputeSelfConsistentEmbeddingCharges();

    //Places New Charges calculated by ComputeEmbeddingCharges() into new charges vector
    new_charges_index = 0;
    for(int j=1; j<=NMon; j++) {
      //Loop through all atoms
      for(int k=0; k<Monomers[j].GetNumberOfAtoms(); k++) {
 
        new_charges_vector[new_charges_index] = Cluster::cluster().GetMonomer(j).GetAtom(k).GetMultipoleMoments().GetMoments().Element(0);
        new_charges_index++;
      }
    }

    if ( Params::Parameters().PrintLevel() >= 1 ) {
      int itr = 0;
      for(int j=1; j<=NMon; j++) {
	//Loop through all atoms
	for(int k=0; k<Monomers[j].GetNumberOfAtoms(); k++) {
	  printf("\t Monomer %d\t atom %d old charges %d\t new charges %d\n", j, k+1, old_charges_vector.Element(itr), new_charges_vector.Element(itr)  );
	  itr++;
	}
      }

    }

    //Updates Charge Difference Vector with new charge differences
    for(int i=0; i<natoms_unitcell; i++) {
      charge_difference[i] = new_charges_vector.Element(i);
    }   
    charge_difference.operator-=(old_charges_vector);


    if ( Params::Parameters().PrintLevel() >= 1 ) {
      int itr = 0;
      printf("\n Charge \t Old Charge \t New Charge \t Difference \n");
      for(int j=1; j<=NMon; j++) {
	//Loop through all atoms
	for(int k=0; k<Monomers[j].GetNumberOfAtoms(); k++) {
	  printf("Monomer %d, atom %d: \t %f \t %f \t %f  \n", j, k+1, old_charges_vector[itr], new_charges_vector.Element(itr), charge_difference.Element(itr));
	  itr++;
	}
      }
    }

    double old_max_abs_deviation = max_abs_deviation;
    double old_rms_deviation = rms_deviation;


    //Finds new max absolute dev. and rms dev.
    max_abs_deviation = fabs(charge_difference.Max(true));
    rms_deviation = charge_difference.RMS();
    printf("%d    \t  %f           \t %f  \n", count, max_abs_deviation, rms_deviation);


    if ( (old_max_abs_deviation < max_abs_deviation && old_rms_deviation < rms_deviation) && Params::Parameters().UseEwald() ) {
      nBadSteps++;
    }

    if ( nBadSteps >= 3) {
      printf("SCRMP convergence errors encountered: 3 bad steps taken\n");
      printf("Run SCE job without EWALD and use SCE_SEED = true to read in those resutls\n");
      exit(1);
    }


    

    count++;
  }


  if ( count > 20 ) {
    printf("Max iterations (20) achieved but convervence criteria not met, check the behavior of the system\n");
  }
  printf("\n\n Number of Iterations for Convergence: %d\n\n\n", count);


  delete [] old_charges;
  delete [] new_charges;
}


void Cluster::ComputeSelfConsistentEmbeddingCharges() {

  if ( Params::Parameters().IsPeriodic() ) {
    //PERIODIC
    CreatePeriodicImageMonomerList( Params::Parameters().GetElectrostaticEmbeddingCutoff() );
    
    AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers, NMon_images, MonomerImages);
    // for ( int imon=1; imon<= NMon_images; imon++ ) {
    //  MonomerImages[imon].SetUseInEmbedding(true);
    // }
    
    if ( Params::Parameters().UseEwald() ) {
      // Using Ewald Embedding:
      //Matrix EwaldCharges = GetEwaldEmbedding(NMon, Monomers, NMon_images, MonomerImages,  unit_cell, reciprocal_cell);
      Matrix EwaldCharges = ComputeEwaldCharges(NMon, Monomers, NMon_images, MonomerImages,  unit_cell);
      for ( int i=1; i<= NMon; i++ ) {
	if  ( Params::Parameters().GetMMType()==98  ) {
	  if ( Params::Parameters().GetQMPackage() == "G09" ) {
	    Monomers[i].CreateG09ChelpGJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
	  } else if ( Params::Parameters().GetQMPackage() == "ORCA" || Params::Parameters().GetQMPackage() == "DALTON" ) {
	    //Monomers[i].CreateOrcaChelpGJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
	    printf("ERROR SCE not compatible with any code other than G09\n");
	    exit(1); 
	  }
	} else if ( Params::Parameters().GetMMType()==97  ) {
	  printf("ERROR ComputeSelfConsistentEmbeddingCharges() not built yet...\n");
	  exit(1);
	} else {
	  printf("ERROR ComputeSelfConsistentEmbeddingCharges() only compatible with CHelpG\n");
	  exit(1);
	}
      }
      
    } else {
      for ( int i=1; i<= NMon; i++ ) {
	if  ( Params::Parameters().GetMMType()==98  ) {
	  if ( Params::Parameters().GetQMPackage() == "G09" ) {
	    Monomers[i].CreateG09ChelpGJob(Monomers, NMon, MonomerImages, NMon_images);
	  } else if ( Params::Parameters().GetQMPackage() == "ORCA" || Params::Parameters().GetQMPackage() == "DALTON" ) {
	    //Monomers[i].CreateOrcaChelpGJob(Monomers, NMon, MonomerImages, NMon_images);
	    printf("ERROR SCE not compatible with any code other than G09\n");
	    exit(1);
	  }
	} else if ( Params::Parameters().GetMMType()==97  ) {
	  Monomers[i].CreateG09HirshfeldJob(Monomers, NMon, MonomerImages, NMon_images);
	}
      }
    }
        
    
  } else{
    /// NON-PERIODIC
    AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers);

    for ( int i=1; i<= NMon; i++ ) {
      if  ( Params::Parameters().GetMMType()==98  ) {
	if ( Params::Parameters().GetQMPackage() == "G09" ) {
	  Monomers[i].CreateG09ChelpGJob(Monomers, NMon );
	} else if ( Params::Parameters().GetQMPackage() == "ORCA" || Params::Parameters().GetQMPackage() == "DALTON" ) {
	  Monomers[i].CreateOrcaChelpGJob();
	}
      } else if ( Params::Parameters().GetMMType()==97  ) {
	if ( Params::Parameters().GetQMPackage() == "G09" ) {
	  Monomers[i].CreateG09HirshfeldJob(Monomers, NMon );
	} else {
	  printf("ERROR ComputeSelfConsistentEmbeddingCharges() Hirshfeld only compatible with G09 right now...\n");
	  exit(1);
	}
      } else {
	printf("ERROR ComputeSelfConsistentEmbeddingCharges() only compatible with CHelpG or Hirshfeld\n");
	exit(1);
      }
    }

    
  }

  
  RunElectrostaticEmbeddingJobs();
  
  for (int imon=1;imon<=NMon;imon++) {
    if  ( Params::Parameters().GetMMType()==98  ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Monomers[imon].ReadG09ChelpGCharges();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" || Params::Parameters().GetQMPackage() == "DALTON" ) {
	Monomers[imon].ReadOrcaChelpGCharges();
      }
    } else if ( Params::Parameters().GetMMType()==97  ) {
      Monomers[imon].ReadHirshfeldCharges();
    }
  }



}

void Cluster::AssignMonomerLists() {
  double LARGE_NUMBER = 999999;
  // CREATE IMAGE MONOMER LIST
  delete [] MonomerImages; // Clear this out just in case
  delete [] DimerImages;
  double MAX_CUTOFF = max(
			  Params::Parameters().GetTwoBodyCutoff(),
			  Params::Parameters().GetElectrostaticEmbeddingCutoff()
			  );

  if (Params::Parameters().IsPeriodic() ) {
  CreatePeriodicImageMonomerList(MAX_CUTOFF);
  
  CreatePeriodicImageDimerList(Params::Parameters().GetTwoBodyCutoff());
  }
  
  // SET ELECTROSTATIC EMBEDDING MONOMERS
  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers, NMon_images, MonomerImages);
  }

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  // ASSIGN CLUSTER JOBS AS NEEDED:
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {

    
    for (int iclust=1; iclust<= NMon_jobs; iclust++ ) {
      int total_charge = 0;
      
      // MONOMER LOOP
      for ( int imon=1; imon<=NMon; imon++ ) {
	Monomers[imon].SetUseInClusterCalculation(false);
	
	if ( Monomers[iclust].FindDistance( Monomers[imon]).Element(0) <= Params::Parameters().GetClusterCutoff()) {
	  Monomers[imon].SetUseInClusterCalculation(true);
	  total_charge += Monomers[imon].GetChargeState();
	}
	
      }
      
      // IMAGE MOMONER LOOP
      if (Params::Parameters().IsPeriodic() ) {
	for (int imon=1; imon<=NMon_images;imon++) {
	  MonomerImages[imon].SetUseInClusterCalculation(false);
	  
	  if ( Monomers[iclust].FindDistance( MonomerImages[imon] ).Element(0) <= Params::Parameters().GetClusterCutoff() ) {
	    MonomerImages[imon].SetUseInClusterCalculation(true);
	    total_charge += MonomerImages[imon].GetChargeState();
	  }
	}
      }
    }
  }


  // ASSIGN DIMER JOBS:
  for ( int idim=1; idim<=NDim; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NMon_jobs; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
      USE_DIMER = false;
    }

    if ( USE_DIMER) {
      Dimers[idim].SetUseInTwoBodyCalculation(true);
    }

  } 

  if (Params::Parameters().IsPeriodic() ) {
  for ( int idim=1; idim<=NDim_images; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( DimerImages[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NMon_jobs; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 ) { // && separation <= Params::Parameters().GetClusterCutoff() ) {
      
      // Actually need to see if it's within the cluster cut off for all of the cluster
      // jobs, if its not then we still need to make the dimer job:
      bool inside = true;
      for (int imon=1; imon<=NMon_jobs; imon++ ) {
	if ( DimerImages[idim].GetMonomerA().GetIndex() == imon) {
	  //double tmp = Monomer[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0);
	  double distance = DimerImages[idim].GetDimerSeparation();
	  if ( distance > Params::Parameters().GetClusterCutoff() ) {
	    inside = false;
	  }
	}
      
      }
      
      if ( inside ){
	USE_DIMER = false;
      }
    }

    if ( USE_DIMER) {
      DimerImages[idim].SetUseInTwoBodyCalculation(true);
    }

  } 
  }

}

void Cluster::CreateFragmentJobs() {

  double LARGE_NUMBER = 999999;

  printf("\nCreating Fragment Jobs to predict magnetic properties\n");


  if ( Params::Parameters().IsPeriodic() ) {
    // Periodic Case:
    

    // if ( Params::Parameters().PrintLevel() >= 1 ) {
    //   printf("\nCreating Image monomer list with monomers out to %f Angstroms... \n",  MAX_CUTOFF);
    //   printf("Crateing Image dimer list with dimers out to %f  Angsroms...\n", Params::Parameters().GetTwoBodyCutoff() );
    // }


    // CREATE IMAGE MONOMER LIST
    delete [] MonomerImages; // Clear this out just in case
    delete [] DimerImages;
    double MAX_CUTOFF = max(
			    Params::Parameters().GetTwoBodyCutoff(),
			    Params::Parameters().GetElectrostaticEmbeddingCutoff()
			    );
    
    CreatePeriodicImageMonomerList(MAX_CUTOFF);
    
    CreatePeriodicImageDimerList(Params::Parameters().GetTwoBodyCutoff());

    AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );
      
    
    
    // SET ELECTROSTATIC EMBEDDING MONOMERS
    if ( Params::Parameters().UseElectrostaticEmbedding() ) {
      AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers, NMon_images, MonomerImages);
    }
    

    
    // ASSIGN MIXED BASIS REGIONS USING CUTOFFS DEFINED FROM AYSMMETRIC CELL
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 &&
	 Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      AssignBasisSetsToAtoms( NMon, Monomers, NMon_images, MonomerImages);
    }

    
    // DETERMINE EWALD POTENTIAL (UNFINISHED)
    if ( Params::Parameters().UseEwald() ) {
      if ( Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
	Matrix EwaldCharges = ComputeEwaldCharges(NMon, Monomers, NMon_images, MonomerImages,  unit_cell);
	//Matrix EwaldCharges = GetEwaldEmbedding(NMon, Monomers, NMon_images, MonomerImages,  unit_cell, reciprocal_cell); 
	
	// CLUSTER-BASED CALCULATION
	if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	  CreateClusterJobs(NMon, Monomers, NMon_images, MonomerImages, EwaldCharges );
	}
	
	
	// CREATE MONOMER AND DIMER JOBS
	CreateMonomerAndDimerJobs(NMon, Monomers, NMon_images, MonomerImages, NDim, Dimers, NDim_images, DimerImages, EwaldCharges); 
      }
      
    } else {

      // CLUSTER-BASED CALCULATION
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	CreateClusterJobs(NMon, Monomers, NMon_images, MonomerImages );
      }
      
      // CREATE MONOMER AND DIMER JOBS
      //if ( Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
      CreateMonomerAndDimerJobs(NMon, Monomers, NMon_images, MonomerImages, NDim, Dimers, NDim_images, DimerImages); 
      //}

    }
 

  } // END Periodic Case
  else {
    // Non- Periodic Case

    // CREATE MONOMER LIST
    double MAX_CUTOFF = max(
			    Params::Parameters().GetTwoBodyCutoff(),
			    Params::Parameters().GetElectrostaticEmbeddingCutoff()
			    );

    if ( Params::Parameters().UseElectrostaticEmbedding() ) {
      AssignMonomerListsForElectrostaticEmbedding(NMon, Monomers);
    }


    // ASSIGN MIXED BASIS REGIONS USING CUTOFFS DEFINED FROM AYSMMETRIC CELL
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 &&
	 Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      AssignBasisSetsToAtoms( NMon, Monomers );
    }


    // CLUSTER-BASED CALCULATION
    if ( Params::Parameters().GetClusterCutoff() > 0 ) {
      CreateClusterJobs(NMon, Monomers  );
    }
    
    // CREATE MONOMER AND DIMER JOBS
    CreateMonomerAndDimerJobs(NMon, Monomers, NDim, Dimers ); 
    

  } // END non-periodic case


}



// Creates & runs the jobs, computes the HMBI energy from the results
void Cluster::RunJobsAndComputeEnergy() {

  printf("\n");
  // Create and Run QM and MM calculations, if requested 
  if ( Params::Parameters().RunJobs() ) {
    // Create the jobs
    if ( !Params::Parameters().BuildForceFieldOnly() && !Params::Parameters().UseFullMMOnly())
      CreateQMJobs();

    if ( !Params::Parameters().NeglectManyBody() && !Params::Parameters().UseFullQMOnly() )
      CreateMMJobs();  

    if (!Params::Parameters().UseFullQMOnly())
      printf("Number of dimer jobs = %d\n",UniqueDim+NDim_images); //HACK

    //exit(1);
    if (Params::Parameters().CreateJobsOnly()) {
      // create the shell scripts for camcasp, if needed, then exit.

      if (Params::Parameters().GetMMType()==2 && !Params::Parameters().NeglectManyBody()) {
	for (int i=1;i<=NMon;i++) {
	  bool UseMMSym = Params::Parameters().UseMMSymmetry();
	  if(Monomers[i].GetSymmetryFactor() !=0 || !UseMMSym) 
	    Monomers[i].RunCamCaspJob(); // create the shell script
	  
	}
      } 
      printf("QM and MM jobs created.  Exiting.\n");
      exit(1);
    }
    // Now run the jobs
    RunJobs();
    // If we only wanted to run the QM & MM jobs, exit now.
    if ( Params::Parameters().RunJobsOnly()) {
      printf("QM and MM jobs complete.  Exiting.\n");
      exit(1);
    }

	/*
    // If we only wanted the force field, exit now
    if ( Params::Parameters().BuildForceFieldOnly()) {
    printf("Analyzing results from previously run QM jobs.\n");
      printf("Force field parameters have been calculated and saved.  Exiting.\n");
      exit(0);
    }
	*/
  } 
  else {
    printf("Analyzing results from previously run QM jobs.\n"); 
    if ( Params::Parameters().GetMMType()==2 ) {
      printf("  Note: Ab initio force field contribution will be re-evaluated.\n");
      if (Params::Parameters().OrientDebug() ) {
	system(RunOrientJob().c_str());;
	printf("Re-running Orient job as well.\n");
      }
    }
  } 


  if (!Params::Parameters().BuildForceFieldOnly() && !Params::Parameters().UseFullMMOnly()){
    ReadQMResults();
  }


  if (!Params::Parameters().NeglectManyBody() && !Params::Parameters().UseFullQMOnly()) 
    ReadMMResults();



  /* Process QM and MM results */

  //Optionally, Compute the analytical HMBI Hessian
  //Calculated before energy to determine the Helmholtz vibrational contribution.
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable()) {
    ComputeQuasiHarmonicFrequencies();//Do QuasiHarmonic Approximation
  }
  else if ( Params::Parameters().DoFreq() && 
	    !Params::Parameters().DoFiniteDifferenceFreqs()) {//Find frequency from the hessian

    Matrix mwHess_HMBI;
    if(Params::Parameters().GetMMType() == 5)//Crystal Hessian
      mwHess_HMBI = Hess_MM;
    else if (Params::Parameters().UseFullQMOnly()) {
      mwHess_HMBI = Hess_QM;
    }
    else{
      Hess_HMBI = ComputeHMBIHessian();   
	
	
	//PrintHessian("\nFull HMBI Hessian", Hess_HMBI);
	//printf("\n");
	//Since Hess_HMBI is not already mass-weighed
      mwHess_HMBI =  ComputeMassWeightedHessian(Hess_HMBI);  
    }
    
    //PrintHessian("\nFull Mass-Weighted HMBI Hessian", mwHess_HMBI);
    ComputeHarmonicFrequencies(mwHess_HMBI);
  }

  // A simple routine for playing with different possible local
  // truncation parameters,
  //printf("Params: ScanLocalTruncationParams = %d\n",Params::Parameters().ScanLocalTruncationParameters() );
  if ( Params::Parameters().ScanLocalTruncationParameters() ) {
    printf("Scanning over different local cutoff values");
    int Nsteps = Params::Parameters().GetNumberOfCutoffScanSteps();
    double c1 = Params::Parameters().GetLocalCutoff(1);
    double c0 = Params::Parameters().GetLocalCutoff(0);
    printf("Initial cutoffs: %.2f / %.2f\n",c1,c0);

    double c1_min = Params::Parameters().GetMinLocalCutoff();
    if (c1_min > c1) {
      printf("ScanLocalTruncationParameter: Error: c1_min = %.2f must be less than c1 = %f\n",c1_min,c1);
      exit(1);
    }

    double stepsize = (c1 - c1_min) / (double) Nsteps;
    int istep = 0;
    while (c1-stepsize >= c1_min) {
      istep++;
      c1 -= stepsize;
      c0 -= stepsize;
      Params::Parameters().SetLocalCutoff(1,c1);
      Params::Parameters().SetLocalCutoff(0,c0);
      printf("Step %d: c1 = %.2f, c0 = %.2f\n",istep,c1,c0);
      
      Energy_HMBI = ComputeHMBIEnergy();
      if ( Params::Parameters().GetMMType() == 4)
	printf("EE-PA Energy = %15.9f hartrees\n",Energy_HMBI);
      else
	printf("HMBI Energy = %15.9f hartrees\n",Energy_HMBI);



     //Vibrational Contribution
    if( (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly() )&& 
        (Params::Parameters().DoFreq() ||
         (Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) ){
      double Temperature = Params::Parameters().GetTemperature();
      printf("\n\n");
      printf("------Thermodynamic Properties-----------------------\n");
      printf("Entropy (at %.0f K)  = %15.9f J/(mol*K)\n",
	     Temperature,Entropy*1000*2625.5);
      printf("Enthalpy (at %.0f K)  = %15.9f Hartrees\n",
	     Temperature,(Energy_HMBI + Temperature*Entropy)); 
      printf("Electronic Enthalpy = %15.9f hartrees\n",Energy_HMBI - Energy_Vib);
      
      printf("Vibrational Energy (at %.0f K)  = %.3f kJ/mol\n",
	     Temperature,(Energy_Vib + Temperature*Entropy)*2625.5);
      printf("constant volume heat capacity (at %.0f K) = %.3f J/(mol*K)\n",Temperature,Cv);
      if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
        printf("Gruneisen parameter = %.3f\n", gruneisen_parameter);
      printf("Gibb's Free Energy (at %.0f K)  = %15.9f Hartrees\n",
	   Temperature,Energy_HMBI); 
     printf("-----------------------------------------------------\n");
    }

      
      printf("\n");
    }
  } // End ScanLocalTruncationParameters code


  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().OrientDebug() ) {
    printf("--------------------------------------------------\n");
    printf("Debug: Orient output:\n");
    double Edebug = ReadOrientEnergy();
    printf("--------------------------------------------------\n\n");
  }

  // Optionally, Compute the analytical HMBI gradient
  // Gradients need for finite difference 
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {

    Grad_HMBI = ComputeHMBIGradient();
    
    //Grad_HMBI.PrintGradient("\nFull HMBI Gradient");
    //printf("\n");
    
    //print time in optimization
    time_t current_time;
    current_time = time(NULL);
    double elapsed_time = difftime(current_time,start_time);
    printf("\n TIME UPDATE: %.0f seconds have elapsed\n",elapsed_time); 
    
    if (Params::Parameters().PrintLevel() > 0 )
      PrintGradient("\nFull HMBI Gradient",Grad_HMBI);
    }

    // Perform HMBI calculations
    Energy_HMBI = ComputeHMBIEnergy();

    if(Energy_HMBI < Energy_HMBI_old){
      //printf("accepting step\n");
      Energy_HMBI_old = Energy_HMBI;
      //AcceptedAtomicCoordinates = AtomicCoordinates; //JLM before symmetry
      AcceptedAtomicCoordinates = SymmetryUniqueCoordinates;
      AcceptedUnitCellAxes = UnitCellAxes;
      AcceptedUnitCellAngles = UnitCellAngles;
    }
    
    if ( Params::Parameters().GetMMType() == 4)
      printf("EE-PA Energy = %15.9f hartrees\n",Energy_HMBI);
    else
      printf("HMBI Energy = %15.9f hartrees\n",Energy_HMBI);
 

     //Vibrational Contribution
    if( (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()) && 
        (Params::Parameters().DoFreq() ||
         (Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) ){
      double Temperature = Params::Parameters().GetTemperature();
      printf("\n\n");
      printf("------Thermodynamic Properaties-----------------------\n");
      printf("Entropy (at %.0f K)  = %15.9f J/(mol*K)\n",
	     Temperature,Entropy*1000*2625.5);
      printf("Enthalpy (at %.0f K)  = %15.9f Hartrees\n",
	     Temperature,(Energy_HMBI + Temperature*Entropy)); 
      printf("Electronic Enthalpy = %15.9f hartrees\n",Energy_HMBI - Energy_Vib);
      
      printf("Vibrational Energy (at %.0f K)  = %.3f kJ/mol\n",
	     Temperature,(Energy_Vib + Temperature*Entropy)*2625.5);
      printf("constant volume heat capacity (at %.0f K) = %.3f J/(mol*K)\n",Temperature,Cv);
      if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
        printf("Gruneisen parameter = %.3f\n", gruneisen_parameter);
      printf("Gibb's Free Energy (at %.0f K)  = %15.9f Hartrees\n",
	   Temperature,Energy_HMBI); 
     printf("-----------------------------------------------------\n");
    }

}

// JDH overload it for the NMR calcs...
void Cluster::CreatePeriodicImageDimerList(double r_cutoff_2bd) {

  int NMonomers = min(NMon,Params::Parameters().GetNumAsymmetricMonomers() );


  NDim_images = 0;
  for (int imon=1; imon<=NMonomers; imon++) {
    for (int jmon=1;jmon<=NMon_images; jmon++ ) {
      double dist = Monomers[imon].FindDistance( MonomerImages[jmon] ).Element(0);
      if ( dist <= r_cutoff_2bd) {
	NDim_images++;
      }

    }
  }
  printf("\t Number of Image Dimers: %d\n", NDim_images);


  DimerImages = new Dimer[NDim_images+1];

  int itr=1;
  for (int imon=1; imon<=NMonomers; imon++) {
    for (int jmon=1;jmon<=NMon_images; jmon++ ) {
      double dist = Monomers[imon].FindDistance( MonomerImages[jmon] ).Element(0);
      if ( dist <= r_cutoff_2bd) {
	DimerImages[itr].Initialize( Monomers[imon], MonomerImages[jmon] );
	itr++;
      }

    }
  }

  
}

// JDH + Watit
string Cluster::ReadG09Section(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$G09") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  if(line.substr(0,1)=="#") {
	    rem += " nosymm Int=(Grid=Ultrafine)";  //symmetry=none Integral(Grid=99590,Acc2E=12) scf=( tight,xqc,maxcycle=256,maxconventionalcycles=100)
	  }
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "$end\n\n";

  infile.clear();
  return rem;
}

string Cluster::ReadDaltonSection(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$DALTON") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();
  return rem;


}

// Read in PSI4 section
string Cluster::ReadPSI4Section(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$PSI4") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();
  return rem;


}

string Cluster::ReadOrcaSection(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$ORCA") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();
  return rem;


}




/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//////                     END JDH MAGNETIC PROPERTIES CODE
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////



void Cluster::PrintHeader() {

    printf("\n");
    printf("\t\t\t    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    printf("\t\t\t  %%%% *                                         * %%%%\n");
    printf("\t\t\t %%%% *  Hybrid Many-Body Interaction QM/MM Code  * %%%%\n");
    printf("\t\t\t%%%% *                                             * %%%%\n");
    printf("\t\t\t%% *         by Greg Beran, Kaushik Nanda,         * %%\n");
    printf("\t\t\t%% *           Ali Sebetci, Kelly Theel,           * %%\n");
    printf("\t\t\t%% *           Yuanhang Huang, Shuhao Wen,         * %%\n");
    printf("\t\t\t%% *           Joshua Hartman, Yin Luo,            * %%\n");
    printf("\t\t\t%% *         Yonaton Heit, Dominique Nocito,       * %%\n");
    printf("\t\t\t%% *      Jessica McKinley, Watit Sontising,       * %%\n");
    printf("\t\t\t%% *             and Chandler Greenwell            * %%\n");
    printf("\t\t\t%% *                                               * %%\n");
    printf("\t\t\t%% *     University of California, Riverside       * %%\n");
    printf("\t\t\t%%%% *                                             * %%%%\n");
    printf("\t\t\t %%%% *                 Oct. 2017                 * %%%%\n");
    printf("\t\t\t  %%%% *                                         * %%%%\n"); 
    printf("\t\t\t    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
}


// Prints out the input file so that we can keep a copy in the output
void Cluster::PrintInputFile(ifstream& infile) {

  printf("----------------------------------------------------------\n");
  printf(" Input file:\n");
  printf("----------------------------------------------------------\n");
  // Rewind the file, just in case
  string line;
  Rewind(infile);
  while ( !infile.eof() ) {
    getline(infile,line);
    printf("%s\n",line.c_str());
  }
  printf("----------------------------------------------------------\n\n");

}

void Cluster::StoreInputFile(ifstream& infile){

  string line;
  Rewind(infile);
  while ( !infile.eof() ) {
   getline(infile,line);
   NumberOfLinesInInput++; 
  }
  //printf("NumberOfLinesInput = %i\n",NumberOfLinesInput);
  InputFile = new string[NumberOfLinesInInput];

  int i = 0;
  Rewind(infile);
  while ( !infile.eof() ) {
   getline(infile,line);
   InputFile[i] = line;
   //printf("%s\n",InputFile[i].c_str());
   i++;
  }  

  //exit(0);

}

void Cluster::CreateWorkingDirectories() {
    // Check if QM and MM directories exist.  If not, create them. 
    struct stat st;
    // QM directory
    if(stat(Params::Parameters().GetQMPath().c_str(),&st) != 0) {
      printf("Directory %s not found.  Creating it.\n",
	     Params::Parameters().GetQMPath().c_str());
      string cmd = "mkdir " + Params::Parameters().GetQMPath();
      system(cmd.c_str());
    }
    else 
      printf("Directory %s exists.\n",
	     Params::Parameters().GetQMPath().c_str());

    // MM directory
    if (! Params::Parameters().NeglectManyBody() ) {
      if(stat(Params::Parameters().GetMMPath().c_str(),&st) != 0) {
	printf("Directory %s not found.  Creating it.\n",
	       Params::Parameters().GetMMPath().c_str());
	string cmd = "mkdir " + Params::Parameters().GetMMPath();
	system(cmd.c_str());
      }
      else{ 
	printf("Directory %s exists.\n",
	       Params::Parameters().GetMMPath().c_str());
	// For Orient monomer directories --- by Ali
	if ( Params::Parameters().GetMMType()==2 && Params::Parameters().RunJobs() ) { 
	  string cmd = "rm -rf " + Params::Parameters().GetMMPath() + "/* " ;
	  system(cmd.c_str());
	}
	
      }
    }
    // directories for handling hessian QM/MM files plus storing current hessian/ykmat/skmat
    //if ( (Params::Parameters().DoFreqAfterOpt())  ) { // i.e if we are using optimizer created by KDN and GB for constrained opt or DLFIND
    if(Params::Parameters().DoFreqAfterOpt() || Params::Parameters().DoFreq() ||
       (Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable() )){
      if(stat(Params::Parameters().GetHessianPath().c_str(),&st) != 0) {
	printf("creating hessian Directors");
        printf("Directory %s not found.  Creating it.\n",
               Params::Parameters().GetHessianPath().c_str());
        string cmd = "mkdir " + Params::Parameters().GetHessianPath();
        system(cmd.c_str());
      }
      else {
        printf("Directory %s exists.\n",
               Params::Parameters().GetHessianPath().c_str());
      }

      // QM directory
      if(stat(Params::Parameters().GetHessianQMPath().c_str(),&st) != 0) {
        printf("Directory %s not found.  Creating it.\n",
               Params::Parameters().GetHessianQMPath().c_str());
        string cmd = "mkdir " + Params::Parameters().GetHessianQMPath();
        system(cmd.c_str());
      }
      else {
        printf("Directory %s exists.\n",
               Params::Parameters().GetHessianQMPath().c_str());
      }

      // MM directory
      if (! Params::Parameters().NeglectManyBody() ) {
        if(stat(Params::Parameters().GetHessianMMPath().c_str(),&st) != 0) {
          printf("Directory %s not found.  Creating it.\n",
                 Params::Parameters().GetHessianMMPath().c_str());
          string cmd = "mkdir " + Params::Parameters().GetHessianMMPath();
	  cmd += ";";
	  system(cmd.c_str());
	}
	else {
	  printf("Directory %s exists.\n",
		 Params::Parameters().GetHessianMMPath().c_str());
	}
	//Full MM tinker found by finite difference Hessian
	if(Params::Parameters().GetMMType()==1 && Params::Parameters().Do_fdTinkerHessian()){
	   //&& Params::Parameters().IsPeriodic()){
	  string full_dir = Params::Parameters().GetHessianMMPath() + "/full" ;
	  if( stat(full_dir.c_str(),&st) != 0){
	    printf("Directory %s not found.  Creating it.\n",
		   full_dir.c_str());
	    string cmd = "mkdir " + full_dir + ";";	 
	    //printf("%s\n",cmd.c_str());
	    system(cmd.c_str());
	  }
	  else
	    printf("Directory %s/full exists.\n",
		   Params::Parameters().GetHessianMMPath().c_str());
	}
	// For Orient monomer directories --- by Ali
	if ( Params::Parameters().GetMMType()==2 && Params::Parameters().RunJobs() ) {
	  string cmd = "rm -rf " + Params::Parameters().GetHessianMMPath() + "/* " ;
	  system(cmd.c_str());
	}
      }
    }
    // Tmp directories for handling approx_hessian QM/MM files plus storing current hessian/ykmat/skmat
    if ( (Params::Parameters().GetJobType()==2)  ) { // i.e if we are using optimizer created by KDN and GB for constrained opt or DLFIND
      if(stat(Params::Parameters().GetTmpFilesPath().c_str(),&st) != 0) {
        printf("Directory %s not found.  Creating it.\n",
               Params::Parameters().GetTmpFilesPath().c_str());
        string cmd = "mkdir " + Params::Parameters().GetTmpFilesPath();
        system(cmd.c_str());
      }
      else {
        printf("Directory %s exists.\n",
               Params::Parameters().GetTmpFilesPath().c_str());
      }
    }

    printf("\n");   
}

//create custom working direction
void Cluster::CreateWorkingDirectories(string Dir_name,bool hess_job){

    // Check if QM and MM directories exist.  If not, create them. 
    struct stat st;

    //Pathway for QM    
    string path;
   if(hess_job)
      path = Params::Parameters().GetHessianQMPath();
   else
    path = Params::Parameters().GetQMPath();
   path += "/" + Dir_name;
   if(stat(path.c_str(),&st) != 0) {
     printf("Directory %s not found.  Creating it.\n",path.c_str());
     string cmd = " mkdir " + path;
     system(cmd.c_str());
   }
   else {
        printf("Directory %s exists.\n",
               path.c_str());
   }

   // MM directory
   if (! Params::Parameters().NeglectManyBody() ) {
    if(hess_job)
      path = Params::Parameters().GetHessianMMPath();
    else
     path = Params::Parameters().GetMMPath();

   path += "/" + Dir_name;

     //printf("\nMM path = %s\n\n",path.c_str());
     if(stat(path.c_str(),&st) != 0) {
       printf("Directory %s not found.  Creating it.\n",path.c_str());
       string cmd = " mkdir " + path;
       system(cmd.c_str());
      }
      else {
        printf("Directory %s exists.\n",
               path.c_str());
      }

      //Full MM tinker found by finite difference Hessian found by finite difference
      if(Params::Parameters().GetMMType()==1 && Params::Parameters().Do_fdTinkerHessian()){
       //&& Params::Parameters().IsPeriodic()){
    

       string full_dir = path + "/full" ;
       if( stat(full_dir.c_str(),&st) != 0){
         printf("Directory %s not found.  Creating it.\n",
	       full_dir.c_str());
         string cmd = "mkdir " + full_dir;	  
         printf("%s\n",cmd.c_str());
         system(cmd.c_str());
       }
       else{
        printf("Directory %s exists.\n",
               full_dir.c_str());
      }
    }
   }

} 


// Scans the input file to determine the number of monomers
int Cluster::FindNumberOfMonomers(ifstream& infile) {
  int num_mon=0;
  string line;
  bool molec_sxn = false;

  // Rewind the file, just in case
  Rewind(infile);
  // Count the monomers
  while ( !infile.eof() ) {
    getline(infile,line);

    if (line.substr(0,9) == "$molecule") {
      molec_sxn = true;
    }

    if (line.substr(0,2) == "--" && molec_sxn==true) {
      num_mon++;
    }
    if (line.substr(0,4)=="$end" && molec_sxn==true) {
      break;
    }
  }
  infile.clear();
  return num_mon;
}

// Scans the input file for a $comment section, from which it gets a title\n"
string Cluster::GetJobTitle(ifstream& infile) {
  string title = "HMBI job\n"; // Generic default title
  string line;

  // Rewind the file, just in case
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,8)=="$comment") {
      title = "";
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  title += line;
	  title += "\n";
	}
	else
	  break;
      }
    }
  }
  infile.clear();
  return title;
}

void Cluster::ReadHMBIParameters(ifstream& infile) {
  string line;
  // Rewind the file, just in case
  Rewind(infile);
  
  string rem_hmbi;
  
  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,5)=="$hmbi") {
      rem_hmbi += line + "\n";
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) == "$end") {
	  rem_hmbi += line + "\n";
	  break;
	}
	else if ( !line.empty() ) {
	  
	  istringstream iss(line);
	  int i = 0;
	  string parameter, value, read;
	  iss >> parameter;
	  iss >> value;
	  if (value =="=")
	    iss >> value;
	  
	  //printf("Parameter = %s, value = %s\n",parameter.c_str(),
	  // value.c_str());
	  
	  // Test for success?
	  
	  // Set the parameter value
	  Params::Parameters().SetParameter(parameter,value);
	  
	  rem_hmbi += line + "\n";
	}
      }
    }
  }
  infile.clear();
  
  string param = "HMBI_PARAMS";
  Params::Parameters().SetParameter(param,rem_hmbi);
  
  // Optionally print out the parameters
  if (Params::Parameters().PrintLevel() > 1) Params::Parameters().Print();
  
}

// Reads unit cell in terms of unit cell side lengths and angles
// Expects to find 6 data points: 3 axis lengths a,b,c, on 1 line and
// three angles alpha, beta, gamma on next line.  Alpha is angle
// between axes b & c, Beta is between axes a & c, etc.
void Cluster::ReadUnitCellDimensions(ifstream& infile) {

  string line;
  // Rewind the file, just in case
  Rewind(infile);
  
  bool found = false; // flag for if unit cell section found

  // Initialize the unit cell vector list
  unit_cell = new Vector[3];
  reciprocal_cell = new Vector[3];
  for (int i=0;i<3;i++) {
    unit_cell[i].Initialize(3);
    reciprocal_cell[i].Initialize(3);
  }
  double a,b,c,alpha,beta,gamma;

  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,10)=="$unit_cell") {
      found = true;
      getline(infile,line);
      
      //printf("line = %s\n",line.c_str());
      { 
	// read in axis lengths
	istringstream iss(line);
	iss >> a;
	iss >> b;
	iss >> c;
      }
      getline(infile,line);
      { 
	// read in lattice angles
	istringstream iss(line);
	iss >> alpha;
	iss >> beta;
	iss >> gamma;
	
      }
      break;
    }
  }

  //lock lattice parameters under symmetry
  if(Params::Parameters().UseLatticeSymmetry()){
    printf("\n\nChecking symmetry of lattice parameters\n");
    
    double tolerance = Params::Parameters().GetSymmetryTolerance();
    //symmetry of the unit cell axes
    if(fabs(a-b) < tolerance){
      b_locked_by_a = 1;
      printf("length of a and b are the same\n");
    }
    if(fabs(a-c) < tolerance){
      c_locked_by_a = 1;
      printf("length of a and c are the same\n");
    }else if(fabs(b-c) < tolerance){
      c_locked_by_b = 1;
      printf("length of b and c are the same\n");
    }

    //symmetry of the unit cell angles. Each angle is first checked to see if they are fixed
    //JLM testing what happens if I decouple the locking and the checking which angles are the same
    
    //Is alpha locked
    if(fabs(alpha-90.0) < tolerance  ||fabs(alpha-120.0) < tolerance){
      lock_alpha = 1;
      printf("angle alpha locked to %f\n",alpha);
    }
    //Is beta locked?
    if(fabs(beta-90.0) < tolerance || fabs(beta-120.0) < tolerance){
      lock_beta = 1;
      printf("angle beta locked to %f\n",beta);
    }//Is beta the same angle as alpha? 
    else if(fabs(beta-alpha) < tolerance){
      beta_locked_by_alpha = 1;
      printf("angles alpha and beta are the same\n");
    }//Is gamma locked?
    if(fabs(gamma-90.0) < tolerance || fabs(gamma-120.0) < tolerance){
      lock_gamma = 1;
      printf("angles gamma locked to %f\n",gamma);
    }//Is gamma and alpha the same angle
    else if (fabs(gamma-alpha) < tolerance){
      gamma_locked_by_alpha = 1;
      printf("angles alpha and gamma are the same\n");
    }
    else if(fabs(gamma-beta)< tolerance){
      gamma_locked_by_beta = 1;
      printf("angles beta and gamma are the same\n");
    }

    
    //Is alpha locked
    /*if(fabs(alpha-90.0) < tolerance  ||fabs(alpha-120.0) < tolerance){
      lock_alpha = 1;
      printf("angle alpha locked to %f\n",alpha);
    }
    //Is beta locked?
    if(fabs(beta-90.0) < tolerance || fabs(beta-120.0) < tolerance){
      lock_beta = 1;
      printf("angle beta locked to %f\n",beta);
    }
    //Is gamma locked?
    if(fabs(gamma-90.0) < tolerance || fabs(gamma-120.0) < tolerance){
      lock_gamma = 1;
      printf("angle gamma locked to %f\n",gamma);
    }
    //Are gamma and alpha the same angle
    if (fabs(gamma-alpha) < tolerance){
      gamma_locked_by_alpha = 1;
      printf("angles alpha and gamma are the same\n");
    }
    //Are beta and alpha the same angle
    else if(fabs(beta-alpha) < tolerance){
      beta_locked_by_alpha = 1;
      printf("angles alpha and beta are the same\n");
    }
    //Are gamma and beta the same angle
    if(fabs(gamma-beta)< tolerance){
      gamma_locked_by_beta = 1;
      printf("angles beta and gamma are the same\n");
    }*/
  }
  else{
    printf("Not using lattice symmetry");
  }
  
  // Store the basic angles/axes data
  UnitCellAngles.Initialize(3);
  UnitCellAxes.Initialize(3);
  UnitCellAngles[0] = alpha;
  UnitCellAngles[1] = beta;
  UnitCellAngles[2] = gamma;
  UnitCellAxes[0] = a;
  UnitCellAxes[1] = b;
  UnitCellAxes[2] = c;
  
  
  
  // Now convert these to lattice vectors
  // first define a few helpful intermediates
  double alpha_rad = alpha*DegreesToRadians;
  double beta_rad = beta*DegreesToRadians;
  double gamma_rad = gamma*DegreesToRadians;
  
  double beta_term = 
    (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
  double gamma_term = 
    sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);
  
  // v1
  unit_cell[0][0] = a;
  unit_cell[0][1] = 0;
  unit_cell[0][2] = 0;
  
  // v2
  unit_cell[1][0] = b*cos(gamma_rad);
  unit_cell[1][1] = b*sin(gamma_rad);
  unit_cell[1][2] = 0;
  
  // v3
  unit_cell[2][0] = c*cos(beta_rad);
  unit_cell[2][1] = c*beta_term;
  unit_cell[2][2] = c*gamma_term;

  // shuhao: create the reciprocal cell vectors;
  
  // calculate the vector products of cell vectors
  double AxB_0 = unit_cell[0][1]*unit_cell[1][2]-unit_cell[0][2]*unit_cell[1][1];
  double AxB_1 = unit_cell[0][2]*unit_cell[1][0]-unit_cell[0][0]*unit_cell[1][2];
  double AxB_2 = unit_cell[0][0]*unit_cell[1][1]-unit_cell[0][1]*unit_cell[1][0];
  
  double BxC_0 = unit_cell[1][1]*unit_cell[2][2]-unit_cell[1][2]*unit_cell[2][1];
  double BxC_1 = unit_cell[1][2]*unit_cell[2][0]-unit_cell[1][0]*unit_cell[2][2];
  double BxC_2 = unit_cell[1][0]*unit_cell[2][1]-unit_cell[1][1]*unit_cell[2][0];
  
  double CxA_0 = unit_cell[2][1]*unit_cell[0][2]-unit_cell[2][2]*unit_cell[0][1];
  double CxA_1 = unit_cell[2][2]*unit_cell[0][0]-unit_cell[2][0]*unit_cell[0][2];
  double CxA_2 = unit_cell[2][0]*unit_cell[0][1]-unit_cell[2][1]*unit_cell[0][0];
  
  // crystal cell volume = Adot(BxC)
  //double CellVol = unit_cell[0][0]*BxC_0 + unit_cell[0][1]*BxC_1 + unit_cell[0][2]*BxC_2;

  //Volume of a column vectors that make triangular matrix is product of the main diagonal
  double CellVol = unit_cell[0][0]*unit_cell[1][1]*unit_cell[2][2];
  SetCellvolume(CellVol);
  
  // A* = (2pi/V)*(BxC); B* = (2pi/V)*(CxA); C* = (2pi/V)*(AxB);
  double dupi_v = 2*3.14159265359/CellVol;
  //A*
  reciprocal_cell[0][0]= dupi_v*BxC_0;
  reciprocal_cell[0][1]= dupi_v*BxC_1;
  reciprocal_cell[0][2]= dupi_v*BxC_2;
  
  //B*
  reciprocal_cell[1][0]= dupi_v*CxA_0;
  reciprocal_cell[1][1]= dupi_v*CxA_1;
  reciprocal_cell[1][2]= dupi_v*CxA_2;
  
  //C*
  reciprocal_cell[2][0]= dupi_v*AxB_0;
  reciprocal_cell[2][1]= dupi_v*AxB_1;
  reciprocal_cell[2][2]= dupi_v*AxB_2;

  
  if (found) {
    printf("\nUnit cell vectors: in Angstroms\n");
    printf("A: (%f,%f,%f)\n",unit_cell[0][0],unit_cell[0][1],unit_cell[0][2]);
    printf("B: (%f,%f,%f)\n",unit_cell[1][0],unit_cell[1][1],unit_cell[1][2]);
    printf("C: (%f,%f,%f)\n",unit_cell[2][0],unit_cell[2][1],unit_cell[2][2]);
    //shuhao: reciprocal cell vector
    printf("\nUnit reciprocal  cell vectors: in Angstroms^-1\n");
    printf("A*: (%f,%f,%f)\n",reciprocal_cell[0][0],reciprocal_cell[0][1],reciprocal_cell[0][2]);
    printf("B*: (%f,%f,%f)\n",reciprocal_cell[1][0],reciprocal_cell[1][1],reciprocal_cell[1][2]);
    printf("C*: (%f,%f,%f)\n",reciprocal_cell[2][0],reciprocal_cell[2][1],reciprocal_cell[2][2]);
    // cell volume
    printf("Unit cell volume: %f angstrom^3; CellVol: %f angstrom^3\n\n",cell_volume,CellVol);   
  }
  else {
    printf("Error: Cluster::ReadUnitCellDimensions(): $unit_cell section not found in input file!\n");
    exit(1);
  }
  
}


// Reads unit cell in terms of 3 lattice vectors, 1 vector per line 
// of the input file
void Cluster::OldReadUnitCellDimensions(ifstream& infile) {
  
  string line;
  // Rewind the file, just in case
  Rewind(infile);
  
  //double unit_cell[3][3]; // stores 3 unit cell vectors
  // first dim = i-th vector, second dim = elements of vector

  bool found = false; // flag for if unit cell section found

  unit_cell = new Vector[3];
  reciprocal_cell = new Vector[3];

  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,10)=="$unit_cell") {
      found = true;
      for (int i=0;i<3;i++) {
	unit_cell[i].Initialize(3);
        reciprocal_cell[i].Initialize(3);
	getline(infile,line);

	//printf("line = %s\n",line.c_str());
	istringstream iss(line);
	for (int j=0;j<3;j++) {
	  iss >> unit_cell[i][j];
	}
      }
    }
  }

  // shuhao: create the reciprocal cell vectors; 
  
  // calculate the vector products of cell vectors
  double AxB_0 = unit_cell[0][1]*unit_cell[1][2]-unit_cell[0][2]*unit_cell[1][1];
  double AxB_1 = unit_cell[0][2]*unit_cell[1][0]-unit_cell[0][0]*unit_cell[1][2];
  double AxB_2 = unit_cell[0][0]*unit_cell[1][1]-unit_cell[0][1]*unit_cell[1][0];
  
  double BxC_0 = unit_cell[1][1]*unit_cell[2][2]-unit_cell[1][2]*unit_cell[2][1];
  double BxC_1 = unit_cell[1][2]*unit_cell[2][0]-unit_cell[1][0]*unit_cell[2][2];
  double BxC_2 = unit_cell[1][0]*unit_cell[2][1]-unit_cell[1][1]*unit_cell[2][0];
  
  double CxA_0 = unit_cell[2][1]*unit_cell[0][2]-unit_cell[2][2]*unit_cell[0][1];
  double CxA_1 = unit_cell[2][2]*unit_cell[0][0]-unit_cell[2][0]*unit_cell[0][2];
  double CxA_2 = unit_cell[2][0]*unit_cell[0][1]-unit_cell[2][1]*unit_cell[0][0];
  
  
  // crystal cell volume = Adot(BxC)
  double CellVol = unit_cell[0][0]*BxC_0 + unit_cell[0][1]*BxC_1 + unit_cell[0][2]*BxC_2;
  SetCellvolume(CellVol);
  
  // A* = (2pi/V)*(BxC); B* = (2pi/V)*(CxA); C* = (2pi/V)*(AxB);
  double dupi_v = 2*3.14159265359/CellVol;
  
  //A*
  reciprocal_cell[0][0]= dupi_v*BxC_0;
  reciprocal_cell[0][1]= dupi_v*BxC_1;
  reciprocal_cell[0][2]= dupi_v*BxC_2;
  //B*
  reciprocal_cell[1][0]= dupi_v*CxA_0;
  reciprocal_cell[1][1]= dupi_v*CxA_1;
  reciprocal_cell[1][2]= dupi_v*CxA_2;
  //C*
  reciprocal_cell[2][0]= dupi_v*AxB_0;
  reciprocal_cell[2][1]= dupi_v*AxB_1;
  reciprocal_cell[2][2]= dupi_v*AxB_2;


  if (found) {
    printf("Unit cell vectors: in Angstroms\n");
    printf("A: (%f,%f,%f)\n",unit_cell[0][0],unit_cell[0][1],unit_cell[0][2]);
    printf("B: (%f,%f,%f)\n",unit_cell[1][0],unit_cell[1][1],unit_cell[1][2]);
    printf("C: (%f,%f,%f)\n",unit_cell[2][0],unit_cell[2][1],unit_cell[2][2]);
    //shuhao: reciprocal cell vector 
    printf("Unit reciprocal  cell vectors: in Angstroms^-1\n");
    printf("A*: (%f,%f,%f)\n",reciprocal_cell[0][0],reciprocal_cell[0][1],reciprocal_cell[0][2]);
    printf("B*: (%f,%f,%f)\n",reciprocal_cell[1][0],reciprocal_cell[1][1],reciprocal_cell[1][2]);
    printf("C*: (%f,%f,%f)\n",reciprocal_cell[2][0],reciprocal_cell[2][1],reciprocal_cell[2][2]);
    // cell volume
    printf("cell_volume: %f angstrom^3; CellVol: %f angstrom^3\n",cell_volume,CellVol);
  }
  else {
    printf("Error: Cluster::ReadUnitCellDimensions(): $unit_cell section not found in input file!\n");
    exit(1);
  }

}

Vector Cluster::ComputeCellParametersFromLatticeVectors(Matrix latticeVec) {

 Vector lattice_param_new(6);
 double a = latticeVec.GetRowVector(0).Norm();
 double b = latticeVec.GetRowVector(1).Norm();
 double c = latticeVec.GetRowVector(2).Norm();
 double alpha = acos(latticeVec.GetRowVector(1).DotProduct(latticeVec.GetRowVector(2))/(b*c))*RadiansToDegrees;
 double beta =  acos(latticeVec.GetRowVector(0).DotProduct(latticeVec.GetRowVector(2))/(a*c))*RadiansToDegrees;
 double gamma = acos(latticeVec.GetRowVector(0).DotProduct(latticeVec.GetRowVector(1))/(a*b))*RadiansToDegrees;
 //printf("ComputeCellParametersFromLatticeVectors: a = %f, b = %f, c = %f\n",a,b,c);
 //printf("ComputeCellParametersFromLatticeVectors: alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);

 lattice_param_new[0] = a;
 lattice_param_new[1] = b;
 lattice_param_new[2] = c;
 lattice_param_new[3] = alpha;
 lattice_param_new[4] = beta;
 lattice_param_new[5] = gamma;

 return lattice_param_new;

}

Matrix Cluster::ComputeLatticeVectorsFromCellParameters(Vector cellParams) {

    Matrix latticeVecs(3,3);
    double a = cellParams[0];
    double b = cellParams[1];
    double c = cellParams[2];
    double alpha = cellParams[3];
    double beta = cellParams[4];
    double gamma = cellParams[5];
    double alpha_rad = alpha*DegreesToRadians;
    double beta_rad = beta*DegreesToRadians;
    double gamma_rad = gamma*DegreesToRadians;
    double beta_term = (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
    double gamma_term = sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);
    
    /*
      cout << "after \n" <<"\n";
      cout << "a = " << a << "\n"; 
      cout << "b = " << b << "\n";
      cout << "c = " << c << "\n";
      cout << "alpha = " << alpha << "\n";
      cout << "beta = " << beta << "\n";
      cout << "gamma = " << gamma << "\n";
      fflush(stdout);
    */    
    
    // v1
    latticeVecs(0,0) = a;
    latticeVecs(0,1) = 0;
    latticeVecs(0,2) = 0;
    
    // v2
    latticeVecs(1,0) = b*cos(gamma_rad);
    latticeVecs(1,1) = b*sin(gamma_rad);             
    latticeVecs(1,2) = 0;
    
    // v3
    latticeVecs(2,0) = c*cos(beta_rad);
    latticeVecs(2,1) = c*beta_term;    
    latticeVecs(2,2) = c*gamma_term;
    
    //cell_volume = latticeVecs[0][0]*latticeVecs[1][1]*latticeVecs[2][2];
    return latticeVecs;

}

double Cluster::GetUnitCellParameter(string type) {

  // Case 1: we have the 3 lattice vectors, need to compute a, b, c,
  // alpha, beta, gamma
  if (Params::Parameters().ReadLatticeVectors() ) {
    double a = unit_cell[0].Norm();
    double b = unit_cell[1].Norm();
    double c = unit_cell[2].Norm();
    double alpha = acos(unit_cell[1].DotProduct(unit_cell[2])/(b*c))*RadiansToDegrees;
    double beta =  acos(unit_cell[0].DotProduct(unit_cell[2])/(a*c))*RadiansToDegrees;
    double gamma = acos(unit_cell[0].DotProduct(unit_cell[1])/(a*b))*RadiansToDegrees;
    printf("GetUnitCellParameter(): a = %f, b = %f, c = %f\n",a,b,c);
    printf("GetUnitCellParameter(): alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);

    if (type== "a" || type == "A") {return a;}
    else if (type == "b" || type == "B") {return b;}
    else if (type == "c" || type == "C") {return c;}
    else if (type == "alpha" || type == "ALPHA") {return alpha;}
    else if (type == "beta" || type == "BETA") {return beta;}
    else if (type == "gamma" || type == "GAMMA") {return gamma;}
    else {
      printf("ERROR: Cluster::GetUnitCellParameter(): Unknown parameter %s\n",type.c_str());
      exit(1);
    }
  }
  // Case 2: Just return the a, b, c, alpha, beta, gamma we already
  // have:
  else {
    if (type == "a" || type == "A") {return UnitCellAxes[0];}
    else if (type == "b" || type == "B") {return UnitCellAxes[1];}
    else if (type == "c" || type == "C") {return UnitCellAxes[2];}
    else if (type == "alpha" || type == "ALPHA") {return UnitCellAngles[0];}
    else if (type == "beta" || type == "BETA") {return UnitCellAngles[1];}
    else if (type == "gamma" || type == "GAMMA") {return UnitCellAngles[2];}
    else {
      printf("ERROR: Cluster::GetUnitCellParameter(): Unknown parameter %s\n",type.c_str());
      exit(1);
    }
  }
  
}

void Cluster::SetUnitCellParameter(string type, double value) {
  
  if (type == "a" || type == "A") { UnitCellAxes[0] = value;}
  else if (type == "b" || type == "B") { UnitCellAxes[1] = value;}
  else if (type == "c" || type == "C") { UnitCellAxes[2] = value;}
  else if (type == "alpha" || type == "ALPHA") { UnitCellAngles[0] = value;}
  else if (type == "beta" || type == "BETA") { UnitCellAngles[1] = value;}
  else if (type == "gamma" || type == "GAMMA") { UnitCellAngles[2] = value;}
  else {
    printf("ERROR: Cluster::SetUnitCellParameter(): Unknown parameter %s\n",type.c_str());
    exit(1);
  }
  
  
  double alpha_rad = UnitCellAngles[0]*DegreesToRadians;
  double beta_rad = UnitCellAngles[1]*DegreesToRadians;
  double gamma_rad = UnitCellAngles[2]*DegreesToRadians;
  
  double beta_term = 
    (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
  double gamma_term = 
    sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);
  
  // v1
  unit_cell[0][0] = UnitCellAxes[0];
  unit_cell[0][1] = 0;
  unit_cell[0][2] = 0;

  // v2
  unit_cell[1][0] = UnitCellAxes[1]*cos(gamma_rad);
  unit_cell[1][1] = UnitCellAxes[1]*sin(gamma_rad);
  unit_cell[1][2] = 0;
  
  // v3
  unit_cell[2][0] = UnitCellAxes[2]*cos(beta_rad);
  unit_cell[2][1] = UnitCellAxes[2]*beta_term;
  unit_cell[2][2] = UnitCellAxes[2]*gamma_term;
  
  printf("Lattice Params after moving\n");
  printf("UnitCellLength %f %f %f\n",
	 UnitCellAxes[0],UnitCellAxes[1],UnitCellAxes[2]);
  printf("UnitCellAngles %f %f %f\n",
	 UnitCellAngles[0],UnitCellAngles[1],UnitCellAngles[2]);
  
  
  
}

void Cluster::SetUnitCellVectors(Vector v1, Vector v2, Vector v3) {
  unit_cell[0] = v1;
  unit_cell[1] = v2;
  unit_cell[2] = v3;
}



// Generate a list of NMon+1 that contains the number of atoms in the full
// cluster (element 0) and in each monomer (elements 1->NMon).
void Cluster::CountNumberOfAtoms(ifstream& infile) {

  string line;
  // Rewind the file, just in case
  Rewind(infile);
  
  // Allocate array to store number of atoms per monomer
  // Element 0 is for full cluster, and monomers count from 1->NMon
  Natoms_per_monomer = new int[NMon+1];  
  
  // First create a list of how many atoms we have per monomer and total
  int index = 0; // monomer counter
  int natoms = 0; // atom counter within a monomer
  bool molec_sxn = false;
  
  // Note: Natoms_per_monomer[0] will be nonsense in this loop.  It gets fixed
  // after the loop is done.
  while ( !infile.eof() ) {
    getline(infile,line);

    // Identify molecule section, and read up to first monomer
    if (line.substr(0,9)=="$molecule") {
      molec_sxn = true;
    }
    
    if (line.substr(0,2)=="--" && molec_sxn==true) {
      Natoms_per_monomer[index] = natoms-1; // it overcounts by 1  
      natoms = 0; 
      getline(infile,line); // throw away charge/spin line
      index++; // increment monomer counter
    }
    
    // If we have reached the end of the section, break
    if (line.substr(0,4)=="$end" && molec_sxn==true) {
      Natoms_per_monomer[index] = natoms-1; // it overcounts by 1  
      break;
    }

    natoms++;
  }

  Natoms_per_monomer[0] = 0;  
  for (int i=1;i<=NMon;i++) {
    Natoms_per_monomer[0] += Natoms_per_monomer[i];
  }


}

int Cluster::GetTypesOfAtoms(){ //JLM 
  int tmp; 
  int atomicNum;
  int lengthArray;
  bool match = false;

  vector<int> check;
  int numUniqueAtoms;

  //count up total number of atom types
  for ( int i=1; i<=NMon; i++){
    tmp = Monomers[i].GetNumberOfAtoms();
    for ( int j=0; j<tmp; j++){
      atomicNum = Monomers[i].GetAtom(j).GetAtomicNumber();
      //printf("Monomer %d, atom %d = %d\n",i,j,atomicNum);
      //fflush(stdout);
      lengthArray=check.size();

      if(i==1 && j==0) {
        match=false;
      }
      else{
        for (int k=0; k<lengthArray; k++){
          //printf("check[k] = %d\natomicNum = %d\n",check[k],atomicNum);
          //fflush(stdout);
          if(check[k]==atomicNum){
            match=true;
          }
        }
      }
      if(!match){
        numUniqueAtoms++;
        check.push_back(atomicNum);
      }
      match=false;
    }
  }
  //printf("numUniqueAtoms = %d\n",numUniqueAtoms);

  return numUniqueAtoms;
}

void Cluster::ReadGeometry(ifstream& infile) {

  // Right now we have 2 options: Tinker or Q-Chem-style XYZ geometries
  if (Params::Parameters().GetMMType()==1) {
    ReadTinkerGeometry(infile);
  }
  else if (Params::Parameters().GetMMType()==2)  {
    ReadSimpleXYZGeometry(infile);
  } 
  else {
    ReadSimpleXYZGeometry(infile);
  }

 
  // Now set Center of Mass and Total Mass;
  TotalMass = 0.0;
  for (int i=1;i<=NMon;i++) {
    TotalMass += Monomers[i].GetMonomerMass();
  }
  
  FindCenterOfMass();
}

// Read in the geometries/charges/spin states for each monomer
// Use this info to initialize the Monomer objects
// Looks for tinker-format XYZ
void Cluster::ReadTinkerGeometry(ifstream& infile) {
  string line;
  
  // Allocate arrays to store charges, spin states, and # of atoms per monomer
  // Element 0 is for full cluster, and monomers count from 1->NMon
  charges = new int[NMon+1];
  spins = new int[NMon+1];
 
  //kdn: need mon type for aiff (date 12th july,2010)
  types = new string[NMon+1];

 // GJB - for AIFF, need ionization potentials
  double *IPs;
  IPs = new double[NMon+1];

  // Allocate arrays to store atom names/positions
  // Confusingly, these start counting at zero;
  int Ntot = GetTotalNumberOfAtoms();
  string symbols[Ntot]; // stores atomic symbol
  double *xyz;
  xyz = new double[3*Ntot]; // stores coordinates
  int *atom_types;
  atom_types = new int [Ntot]; // stores MM atom type

  // and to store connectivity
  int tmp[10]; // temporarily array, for scratch
  int *Nconnected; // list of # of connections for each atom
  Nconnected = new int [Ntot]; // stores MM atom type
  int *Connectivity;
  Connectivity = new int [6*Ntot]; // up to 6 connections

  for (int i=0;i<6*Ntot;i++)
    Connectivity[i] = 0;

  Rewind(infile);

  int index = 0, atom = 0;
  bool molec_sxn = false;

  while ( !infile.eof() ) {
    getline(infile,line);
    
    // Identify molecule section, and read up to first monomer
    if (line.substr(0,9)=="$molecule") {
      molec_sxn = true;
      infile >> charges[index]; infile >> spins[index]; // Full cluster state
    }

    if (line.substr(0,2)=="--" && molec_sxn==true) {
      index++; // increment monomer counter
      infile >> charges[index]; infile >> spins[index];
      // GJB: if need IP and mon type (kdn)
      if (Params::Parameters().GetMMType()==2) {
	infile >> IPs[index];
      }
      

      getline(infile,line); // advance line by one
      int tmp_int;
      for (int i=0;i<GetNumberOfAtoms(index);i++)  {
	getline(infile,line);
	//printf("Line = %s\n",line.c_str());
	istringstream iss(line);
	iss >> tmp_int; // atom number, discarded
	iss >> symbols[atom]; // Atomic symbol
	iss >> xyz[3*atom]; // x
	iss >> xyz[3*atom+1]; // y
	iss >> xyz[3*atom+2]; // z
	iss >> atom_types[atom]; // MM atom type
	
	// Read in connectivity into temporary array
	int item = 0;
	
	//printf("Atom %d is connected to:",tmp_int);	
	while (iss >> Connectivity[atom*6+item] ) {
	  //printf(" (i=%d) %d ",atom*6+item,Connectivity[atom*6+item]);
	  item++;
	}
	//printf("  Total of %d connections\n",item);
	Nconnected[atom] = item; // set number of connections
	atom++; // increment our atom counter
      }
    }
    
    // If we have reached the end of the section, break
    if (line.substr(0,4)=="$end" && molec_sxn==true) {
      break;
    }
  }

  // Read embedding charges, if requested
  Vector embedding_charges;
  if ( Params::Parameters().UseEmbeddingCharges() ) {
    embedding_charges.Initialize(Ntot);
    embedding_charges = ReadEmbeddingCharges(infile);
    //embedding_charges.Print("Embedding charges");
  }
  
  // Now initialize the Monomer objects
  Monomers = new Monomer[NMon+1];

  int offset = 0;// counter for extracing from symbols & xyz arrays
  
  for (int i=1;i<=NMon;i++) {
    int charge = charges[i];
    int spin = spins[i];
    string type = types[i]; 
    int natoms = GetNumberOfAtoms(i);
    
    // Grab subset of cartesian coordinates
    double *mon_xyz;
    mon_xyz = new double[3*natoms];
    for (int j=0;j<3*natoms;j++) {
      mon_xyz[j] = xyz[3*offset+j];
    }
    
    // Grab subset of atomic symbols, atom types, and Nconnected list
    string atoms[natoms];
    int *atypes;
    atypes = new int[natoms];
    int *nconn;
    nconn = new int[natoms];
    for (int j=0;j<natoms;j++) {
      atoms[j] = symbols[offset+j];
      atypes[j] = atom_types[offset+j];
      nconn[j] = Nconnected[offset+j];
    }
    
    // Grab subset of Connectivity list
    int *conn;
    conn = new int[6*natoms];
    for (int j=0;j<6*natoms;j++) {
      conn[j] = Connectivity[6*offset+j];
    }
    
    // Handle embedding charges, if appropriate
    Vector Monomer_charges;
    if ( Params::Parameters().UseEmbeddingCharges() ) {
      Monomer_charges.Initialize(natoms);
      for (int j=0;j<natoms;j++) {
	Monomer_charges[j] = embedding_charges[offset+j];
      }
      
      if ( Params::Parameters().PrintLevel() > 1 ) {
	printf("Monomer %d\n",i);
	Monomer_charges.Print("Charges\n");
      }
    }


    if ( Params::Parameters().UseEmbeddingCharges() ) {
      Monomers[i].Initialize( i, offset, charge, spin, type, atoms, mon_xyz, natoms, atypes,
			      nconn, conn, Monomer_charges);
    }
    
    else { // no embedding charges case 
      Monomers[i].Initialize( i, offset,charge, spin, type, atoms, mon_xyz, natoms, atypes,
			      nconn, conn);     

    }
  

    // Handle Ionization Potentials (AIFF), if necessary
    if (Params::Parameters().GetMMType()==2) {
      Monomers[i].SetIonizationPotential(IPs[i]);
      //if (Params::Parameters().PrintLevel() > 0) 
      printf("Setting IP for monomer %d to %f\n",i,IPs[i]); 
    }
    
    offset += natoms;
    
    delete [] mon_xyz;
    delete [] atypes;
    delete [] nconn;
    delete [] conn;
  }
  
  delete [] IPs;
  delete [] xyz;
  delete [] atom_types;
  delete [] Nconnected;
  delete [] Connectivity;
  
  
}

// Read in the geometries/charges/spin states for each monomer
// Use this info to initialize the Monomer objects
// This is old version, which reads q-chem-style cartesian data
void Cluster::ReadSimpleXYZGeometry(ifstream& infile) {
  string line;
  
  // Allocate arrays to store charges, spin states, and # of atoms per monomer
  // Element 0 is for full cluster, and monomers count from 1->NMon
  double *IPs;
  charges = new int[NMon+1];
  spins = new int[NMon+1];
  IPs = new double[NMon+1];
  types = new string[NMon+1];
  
  // Allocate arrays to store atom names/positions
  // Confusingly, these start counting at zero;
  int Ntot = GetTotalNumberOfAtoms();
  string symbols[Ntot];
  double *xyz;
  xyz = new double[3*Ntot];
  
  Rewind(infile);
  
  int index = 0, atom = 0;
  bool molec_sxn = false;
  while ( !infile.eof() ) {
    getline(infile,line);
    // Identify molecule section, and read up to first monomer
    if (line.substr(0,9)=="$molecule") {
      molec_sxn = true;
      infile >> charges[index]; infile >> spins[index]; // Full cluster state
      
    }
    
    if (line.substr(0,2)=="--" && molec_sxn==true) {

      index++; // increment monomer counter
      infile >> charges[index]; infile >> spins[index];
      // If doing AIFF, read in IP
      if (Params::Parameters().GetMMType()==2) {
	infile >> IPs[index]; infile >> types[index];
      } else if ( Params::Parameters().GetMMType()==99) { // GDMA
	infile >> IPs[index]; infile >> types[index];
      } else if ( Params::Parameters().GetMMType()==98) { // ChelpG
	infile >> IPs[index]; infile >> types[index];
      } else if ( Params::Parameters().GetMMType()==97) { // Hirshfeld
	infile >> IPs[index]; infile >> types[index];
      }

      
      
      for (int i=0;i<GetNumberOfAtoms(index);i++)  {

	infile >> symbols[atom]; // Atomic symbol
	infile >> xyz[3*atom]; // x
	infile >> xyz[3*atom+1]; // y
	infile >> xyz[3*atom+2]; // z
	atom++;

      }
    }

    // If we have reached the end of the section, break
    if (line.substr(0,4)=="$end" && molec_sxn==true) {
      break;

    }
  }


  // Read embedding charges, if requested
    Vector embedding_charges;
    if ( Params::Parameters().UseEmbeddingCharges() ) {
      embedding_charges.Initialize(Ntot);
      embedding_charges = ReadEmbeddingCharges(infile);
      //embedding_charges.Print("Embedding charges");
    }


  // Now initialize the Monomer objects
  Monomers = new Monomer[NMon+1];

  int offset = 0;// counter for extracing from symbols & xyz arrays

  for (int i=1;i<=NMon;i++) {
    int charge = charges[i];
    int spin = spins[i];
    string type = types[i];
    int natoms = GetNumberOfAtoms(i);
    // Grab subset of cartesian coordinates
    double *mon_xyz;
    mon_xyz = new double[3*natoms];
    for (int j=0;j<3*natoms;j++) {
      mon_xyz[j] = xyz[3*offset+j];
    }
    // Grab subset of atomic symbols
    string atoms[natoms];
    for (int j=0;j<natoms;j++) {
      atoms[j] = symbols[offset+j];
    }

    // Handle embedding charges, if appropriate
    Vector Monomer_charges;
    if ( Params::Parameters().UseEmbeddingCharges() ) {
      Monomer_charges.Initialize(natoms);
      for (int j=0;j<natoms;j++) {
	Monomer_charges[j] = embedding_charges[offset+j];
      }
      
      if ( Params::Parameters().PrintLevel() > 1 ) {
	printf("Monomer %d\n",i);
	Monomer_charges.Print("Charges\n");
      }
    }

    if ( Params::Parameters().UseEmbeddingCharges() ) {
      Monomers[i].Initialize(i,offset,charge,spin,type,atoms,mon_xyz,natoms,
			   Monomer_charges);
    }
    else { // no embedding charges case
      if ( Params::Parameters().PrintLevel() > 1 )
        printf("Monomer %i\n ",i);
      Monomers[i].Initialize(i,offset,charge,spin,type,atoms,mon_xyz,natoms);
    }
    // Handle Ionization Potentials (AIFF), if necessary
    if (Params::Parameters().GetMMType()==2) {
      Monomers[i].SetIonizationPotential(IPs[i]);
      //if (Params::Parameters().PrintLevel() > 0) 
	printf("Setting IP for monomer %d to %f\n",i,IPs[i]); 
    }

    offset += natoms;
    delete [] mon_xyz;


  }  
  delete [] IPs;
  delete [] xyz;

} 

// Uses monomer center of masses to find the cluster center of mass;
void Cluster::FindCenterOfMass() {
  /*     COM = Sum(m_i*r_i)/Sum(m_i)   */

  double totalmass = GetTotalMass();

  CenterOfMass.Initialize(3);

  for (int i=1;i<=NMon;i++) {
    Vector mon_com = Monomers[i].GetCenterOfMass(); // pos
    mon_com.Scale( Monomers[i].GetMonomerMass() ); // mass*pos
    CenterOfMass += mon_com;
  }
  CenterOfMass.Scale(1.0/totalmass); // mass*pos/totalmass


  if ( Params::Parameters().PrintLevel() > 2 ) 
    printf("Cluster Center of Mass = (%f, %f, %f)\n",
	 CenterOfMass[0],CenterOfMass[1],CenterOfMass[2]);

}

// Watit - Read $cif section
string Cluster::ReadCIFSection(ifstream& infile) {

	string cif_file;
	string line;
	Rewind(infile);
	while (!infile.eof()) {
		getline(infile,line);
		if (line.substr(0,4)=="$cif") {
			for (int i=0;;i++) {
				getline(infile,line);
				if (line.substr(0,4)!="$end") {
					if (line.substr(0,8)=="cif_file") {
						istringstream iss(line);
                                		string trash;
                                		iss >> trash;
						iss >> cif_file;
					}
				} else {
					break;
				}
			}
		}
	}
	infile.clear();
	return cif_file;
}

void Cluster::CreateInputfromCIF(string ciffile) { //Watit

	string cif_file;
	cif_file = Params::Parameters().GetBasePath() + ciffile;
	string line;
	ifstream infile;
	infile.open(cif_file.c_str());

	cout << cif_file.c_str() << endl;
	if(!infile.is_open()) {
		printf("Cluster::CreateInputfromCIF : Cannot open file '%s'\n",cif_file.c_str());
		exit(1);
	}

	//Get Info from cif file

	int NSym = 0;
	int NAtom = 0;
	int checksymcol = 0;
	int checkatomcol = 0;
	double Lattice[6];
	string SymOperater;
	string FracCoord;

	while (!infile.eof()) {
		getline(infile,line);
		if (line.substr(0,5)=="loop_") {
			for (int i=0;;i++) {
				getline(infile,line);
				if (line.substr(0,26)!="_symmetry_equiv_pos_as_xyz") {
					checksymcol++;
				} else {
					break;
				}
			}
		}
		if (line.substr(0,26)=="_symmetry_equiv_pos_as_xyz") {
			for (int i=0;;i++) {
				getline(infile,line);
				if (line.substr(0,14)!="_cell_length_a") {
					istringstream iss(line);
					string tmp;
					for (int j=0;j<checksymcol;j++) {
						iss >> tmp;
					}
					iss >> tmp;
					SymOperater += "\n";
					SymOperater += tmp;
					NSym++;
				} else {
					break;
				}
			}
			//cout << SymOperater.c_str() << endl;
		}
		if (line.substr(0,14)=="_cell_length_a") {
			for (int i=0;;i++) {
				int found_open = line.find("(");
				if (found_open != -1) {
					for (int i=0;i<line.size();i++) {
						found_open = line.find("(",found_open+1,1);
						if (found_open != -1) {
							line.erase(line.begin()+found_open,line.end());
						}
					}
				}
				for (int j=0;j<6;j++) {
					istringstream iss(line);
					string tmp;
					iss >> tmp;
					iss >> Lattice[j];
					getline(infile,line);
				}
				getline(infile,line);
				if (line.substr(0,5)=="loop_") {
					break;
				}
                        }
		}
		if (line.substr(0,5)=="loop_") {
                        for (int i=0;;i++) {
                                getline(infile,line);
                                if (line.substr(0,22)!="_atom_site_type_symbol") {
                                        checkatomcol++;
                                } else {
                                        break;
                                }
                        }
                }
		if (line.substr(0,18)=="_atom_site_fract_z") {
			for (int i=0;;i++) {
				getline(infile,line);
				if (!line.empty()) {
					istringstream iss(line);
					string tmp;
					for (int j=0;j<checkatomcol;j++) {
						iss >> tmp;
					}
					for (int j=0;j<4;j++) {
						iss >> tmp;
						FracCoord += tmp;
						FracCoord += " ";
					}
					FracCoord += "\n";
					NAtom++;
				} else {
					break;
				}
			}
			cout << FracCoord.c_str() << endl;
		}
	}

	cout << "There are " << NSym << " space group operator(s)" << endl;
	cout << "There are " << NAtom << " atom(s) contained in cif file" << endl; 


	// Get LatticeMatrix
	const double pi = 3.14159265359;
	double vol;
	double LatticeMatrix[9];
	vol = Lattice[0]*Lattice[1]*Lattice[2]*sqrt(1-(cos(Lattice[3]*pi/180))*(cos(Lattice[3]*pi/180))-(cos(Lattice[4]*pi/180))*(cos(Lattice[4]*pi/180))-(cos(Lattice[5]*pi/180))*(cos(Lattice[5]*pi/180))+2*cos(Lattice[3]*pi/180)*cos(Lattice[4]*pi/180)*cos(Lattice[5]*pi/180));

	LatticeMatrix[0] = Lattice[0];
    	LatticeMatrix[1] = 0;
    	LatticeMatrix[2] = 0;
    	LatticeMatrix[3] = Lattice[1]*cos(Lattice[5]*pi/180);
    	LatticeMatrix[4] = Lattice[1]*sin(Lattice[5]*pi/180);
    	LatticeMatrix[5] = 0;
    	LatticeMatrix[6] = Lattice[2]*cos(Lattice[4]*pi/180);
    	LatticeMatrix[7] = Lattice[2]*(cos(Lattice[3]*pi/180)-cos(Lattice[4]*pi/180)*cos(Lattice[5]*pi/180))/sin(Lattice[5]*pi/180);
    	LatticeMatrix[8] = vol/(Lattice[0]*Lattice[1]*sin(Lattice[5]*pi/180));

	cout << "Recieve Unit Cell parameters" << endl;

	// Get SpaceSymMatrix
	double SpaceSymTrans[NSym*3];
	double SpaceSymRot[NSym*9];
	int m = 0;
	int n = 0;
	int o = 0;

	string SpaceSym[NSym];
	istringstream iss(SymOperater);
	while (getline(iss,line)) {
		iss >> SpaceSym[m];
		int found_comma = SpaceSym[m].find(",");
		int last_comma;
		if (found_comma!=-1) {
			for (int j=0;j<SpaceSym[m].size();j++) {
				found_comma = SpaceSym[m].find(",",found_comma+1,1);
				if (found_comma!=-1) {
					SpaceSym[m].replace(found_comma,1," ");
					SpaceSym[m].append(" .");
				}
			}
		}
		m++;
	}

	m = 0;
	string tmp1[3];
	double tmp2[3];
	double tmp3[3];
	double tmp4;
	for (int i=0; i<NSym; i++) {
		istringstream iss(SpaceSym[i]);
		for (int j=0;j<3;j++) {
			iss >> tmp1[j];
			if (tmp1[j].length()==1&&tmp1[j].substr(0,1)!=".") {
				for (int k=0;k<3;k++) {
					tmp2[k] = 0;
					if (tmp1[j]=="x") { tmp2[0] = 1;}
					if (tmp1[j]=="y") { tmp2[1] = 1;}
					if (tmp1[j]=="z") { tmp2[2] = 1;}
					SpaceSymRot[n] = tmp2[k];
					n++;
				}
				tmp3[j] = 0;
				SpaceSymTrans[m] = 0;
				m++;
			} else if (tmp1[j].length()==2) {
				for (int k=0;k<3;k++) {
					tmp2[k] = 0;
					if (tmp1[j]=="-x") { tmp2[0] = -1;}
					if (tmp1[j]=="-y") { tmp2[1] = -1;}
					if (tmp1[j]=="-z") { tmp2[2] = -1;}
					SpaceSymRot[n] = tmp2[k];
					n++;
				}
				tmp3[j] = 0;
				SpaceSymTrans[m] = 0;
				m++;
			} else if (tmp1[j].length()>2) {
				for (int k=0;k<3;k++) {
					tmp2[k] = 0;
					tmp3[k] = 0;
					if (tmp1[j].find("+x")!=-1||tmp1[j].substr(0,1)=="x") { tmp2[0] = 1;}
					if (tmp1[j].find("-x") != -1 || tmp1[j].substr(0,2) == "-x") { tmp2[0] = -1;}
					if (tmp1[j].find("+y") != -1 || tmp1[j].substr(0,1) == "y") { tmp2[1] = 1;}
					if (tmp1[j].find("-y") != -1 || tmp1[j].substr(0,2) == "-y") { tmp2[1] = -1;}
					if (tmp1[j].find("+z") != -1 || tmp1[j].substr(0,1) == "z") { tmp2[2] = 1;}
					if (tmp1[j].find("-z") != -1 || tmp1[j].substr(0,2) == "-z") { tmp2[2] = -1;}
					SpaceSymRot[n] = tmp2[k];
					n++;
				}
				tmp3[j] = 0;
				if (atoi(tmp1[j].substr(0,1).c_str()) != 0 && atoi(tmp1[j].   substr(2,3).c_str()) != 0) {
					tmp4 = atof(tmp1[j].substr(0,1).c_str())/atof(tmp1[j].  substr(2,3).c_str());
					tmp3[j] = tmp4;}
				SpaceSymTrans[m] = tmp4;
				m++;
			} else {
				for (int k=0;k<3;k++) {
					tmp2[k] = 0;
					SpaceSymRot[n] = 0;
					n++;
				}
				tmp3[j] = 0;
				SpaceSymTrans[m] = 0;
				m++;
			}
		}		
	}

	double SpaceSymMatrix[NSym*16];
	for (int k=0;k<NSym;k++) {
		for (int j=0;j<3;j++) {
			for (int i=0;i<3;i++) {
				int l = i+j*3+k*3*3;
				SpaceSymMatrix[o] = SpaceSymRot[l];
				o++;
			}
			int m = j+k*3;
			SpaceSymMatrix[o] = SpaceSymTrans[m];
			o++;
		}
		for (int n=0;n<3;n++) {
			SpaceSymMatrix[o] = 0;
			o++;
		}
		SpaceSymMatrix[o] = 1;
		o++;
	}
/*
    	for (int k=0;k<NSym;k++) {
        	for (int j=0;j<4;j++) {
            		for (int i=0;i<4;i++) {
                		int l=j+i*4+k*4*4;
                		cout <<  SpaceSymMatrix[l] << " ";
			}	
			cout << endl;
		}
		cout << endl;
	}
*/
	cout << "Recieve Space Group Operators" << endl;

	// Get Coord
	m = 0;
	string Element[NAtom];
	double Position[NAtom*3];
        istringstream jss(FracCoord);
        while (getline(jss,line)) {
                int found_open = line.find("(");
		int found_close;
              	if (found_open != -1) {
                 	for (int i=0;i<line.size();i++) {
                              	found_open = line.find("(",found_open+1,1);
                             	if (found_open != -1) {
					found_close = line.find(")",found_open+1,1);
					int del = found_close-found_open+1;
                                      	line.erase(found_open,del);
                              	}
                       	}
		}
		istringstream ss(line);
        	ss >> Element[m];
		for (int i=0;i<3;i++) {
			//string tmpxx;
			//ss >> tmpxx;
			//cout << "Position " << tmpxx << endl;
			ss >> Position[3*m+i];
                }
		m++;
        }

	cout << "Before " << NAtom << endl;
	RemoveRedundanceCoord(NAtom,Element,Position);
	cout << "After " << NAtom << endl;
//////////////////
/*

	// Find Redundance Coord
	int checkcoord[NAtom];
	for (int i=0;i<NAtom;i++) {
		checkcoord[i] = 0;		
	}

	string refElement;
	double refPosition[3];
	for (int i=0;i<NAtom-1;i++) {
		if (checkcoord[i]==0) {
			refElement = Element[i];
			for (int j=0;j<3;j++) {
				refPosition[j] = Position[3*i+j];
			}
			for (int k=i+1;k<NAtom;k++) {
				m = 0;
				if(checkcoord[k]==0) {
					if (refElement==Element[k]) {m++;}
					for (int l=0;l<3;l++) {
						if (abs(refPosition[l]-Position[3*k+l])< 0.0001) {m++;} // threshold
					}
				}
				if (m==4) {checkcoord[k] += 1;}
			}
		}
	}

	cout << "Recieve initial Frac Coord" << endl;

	int tmpNAtom = 0;
	for (int i=0;i<NAtom;i++) {
		if (checkcoord[i]==0) {tmpNAtom++;}
	}
*/
	// Make 4x1 Position's Matrix
/*
	n = 0;
        string tmpElement[tmpNAtom];
        double tmpPosition[tmpNAtom*4];
	for (int i=0;i<NAtom;i++) {
		if (checkcoord[i]==0) {
			tmpElement[n] = Element[i];
			cout << tmpElement[n] << " " << endl;
			for (int j=0;j<4;j++) {
				if (j==3) {
					tmpPosition[3*n+j] = 1;
					cout << tmpPosition[3*n+j] << " " << endl;
				} else {
					tmpPosition[3*n+j] = Position[3*i+j];
					cout << tmpPosition[3*n+j] << " ";
				}
			}
			tmpNAtom++;
			n++;
		}
	}

*//*
	// Operate space sym
	m = 0;
	n = 0;
	double tmpSpaceSymMatrix[NSym*];
	for (int i=0; i<NSym; i++) {
		for (int j=0; j<NAtom; j++) {
			for (int k=0; k<4; k++) {
				for (int l=0; l<4; l++) {					
*/		


}

void Cluster::RemoveRedundanceCoord(int NAtom, string Element[], double Position[]) {  //Watit
	
	//int NAtom = natom;

	cout << "Remove::NAtom " << NAtom << endl;
	for (int i=0;i<NAtom;i++) {
		cout << Element[i] << " ";
		for (int j=0;j<3;j++) {
			cout << Position[3*i+j] << " ";
		}
		cout << endl;
	}

        int checkcoord[NAtom];
        for (int i=0;i<NAtom;i++) {
                checkcoord[i] = 0;
        }

        string refElement;
        double refPosition[3];
        for (int i=0;i<NAtom-1;i++) {
                if (checkcoord[i]==0) {
                        refElement = Element[i];
                        for (int j=0;j<3;j++) {
                                refPosition[j] = Position[3*i+j];
                        }
                        for (int k=i+1;k<NAtom;k++) {
                                int m = 0;
                                if(checkcoord[k]==0) {
                                        if (refElement==Element[k]) {m++;}
                                        for (int l=0;l<3;l++) {
                                                if (abs(refPosition[l]-Position[3*k+l])< 0.0001) {m++;} // threshold
                                        }
                                }
                                if (m==4) {checkcoord[k] += 1;}
                        }
                }
        }

        int tmpNAtom = 0;
        for (int i=0;i<NAtom;i++) {
		cout << "checkcoord " << checkcoord[i] << endl;
                if (checkcoord[i]==0) {tmpNAtom++;}
        }
	NAtom = tmpNAtom;
	cout << "Remove:: " << NAtom << endl;
}

// Read the Q-Chem $rem section
string Cluster::ReadQChemRemSection(ifstream& infile, int type) {

  // The $rem section gets stored as one long string, "rem"
  string rem = "\n"; 
  string line;

  Rewind(infile); // Rewind the file, just in case

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( ( (line.substr(0,6)=="$qchem" && line.substr(0,7)!="$qchem2") && type==1) || (line.substr(0,7)=="$qchem2" && type==2) ) {
       rem +="$rem\n";
      if ( (Params::Parameters().DoForces()) || Params::Parameters().DoFiniteDifferenceFreqs() )
        rem += "jobtype = force\n";
      if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
        rem += "jobtype = freq\n"; 
        rem += "vibman_print = 4\n"; //print out the mass-weighted hessian

      if (Params::Parameters().DoRaman()) //Watit
        rem += "doraman = true\n";
    }
	
    for (int i=0;;i++) {
      getline(infile,line);
      if (line.substr(0,4) != "$end") {
          rem += line;
	  rem += "\n";
      }
      else {
	// Insert keywords if doing gradients/hessians/Raman
	if ( Params::Parameters().GetJobType() > 1 )
	  rem += "no_reorient = true\nsym_ignore = true\n";
	  break;
	}
      }
    }
  }
  rem += "$end\n\n";

  infile.clear();
  return rem; 
  
}


// Read the Qchem $basis section
string Cluster::ReadQChemBasisSection(ifstream& infile) {

  // The $basis section gets stored as one long string, "basis"
  string basis = "\n";
  string line;
  
  Rewind(infile); // Rewind the file, just in case
  
  bool custom_basis = false;
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,7)=="$basis") {
      custom_basis = true;
      basis +="$basis\n";
      for (int i=0;;i++) {
        getline(infile,line);
        if (line.substr(0,4) != "$end") {
          basis += line;
          basis += "\n";
        }
        else
          break;
      }
    }
  }
  
  if (custom_basis)
    basis += "$end\n";
  
  infile.clear();
  return basis;
  
}

// Yoni: Read the basis of the molpro input file,
// This does not include the header
// If doing counterpoise coorection, 
// this section only contain basis
string Cluster::ReadMolProSection(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,7)=="$molpro" && line.substr(0,8) !="$molpro_" ) {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();
  return rem;
}



// Yoni: Read the instructions of the molpro input file,
// This does not include the header
// If doing counterpoise coorection, 
// this section only contain basis
string Cluster::ReadMolProInstSection(ifstream& infile) {
  string inst = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if( (line.find("$molpro_inst") != -1 ) && (line.find("$molpro_inst_hf") == -1  && line.find("$molpro_inst_HF") == -1) ) {
    //if ( line.substr(0,12)=="$molpro_inst") {
      for (int i=0;;i++) {
	getline(infile,line);
	//if (line.substr(0,4) != "$end") {
	if (line.find("$end") == -1) {
	  inst += line;
	  inst += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //inst += "\n\n";

  //printf("%s\n",inst.c_str());
  infile.clear();
  return inst;
}

// Read rem section from PSI4 input file
// This section does not include the header or cartesian coordinates
string Cluster::ReadPSI4RemSection(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$PSI4_rem") {
      for (int i=0;;i++) {
        getline(infile,line);
        if (line.substr(0,4) != "$end") {
          rem += line;
          rem += "\n";
        } else {
          break;
        }
      }
    }
  }
  infile.clear();
  return rem;
}  

// this section only contain basis
string Cluster::ReadMolProCBSBasisSection(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,11)=="$molpro_CBS"||line.substr(0,11)=="$molpro_cbs") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  //printf("%s\n",rem.c_str());

  infile.clear();
  return rem;
}

//Separate HF instructions for the CBS extrapolation
string Cluster::ReadMolProInstHFSection(ifstream& infile) {
  string inst = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,15)=="$molpro_inst_HF"||line.substr(0,15)=="$molpro_inst_hf") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  inst += line;
	  inst += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //inst += "\n\n";
  //printf("%s\n",inst.c_str());

  infile.clear();
  return inst;
}

// Read the basis for a CCSDT correction calculation
string Cluster::ReadMolProCCSDTBasisSection(ifstream& infile){
  string basis = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,19)=="$molpro_ccsdt_basis") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  basis += line;
	  basis += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //basis += "\n\n";
  //printf("%s\n",basis.c_str());

  infile.clear();
  return basis;

}

//Read MP2 section for the CCSD(T) correction
string Cluster::ReadMolProCCSDTMP2Section(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line=="$molpro_ccsdt_mp2") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  //printf("%s\n",rem.c_str());

  infile.clear();
  return rem;
}

// Read CCSD(T) section for CCSD(T) correction
string Cluster::ReadMolProCCSDTInstSection(ifstream& infile){
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line=="$molpro_ccsdt_inst") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  //printf("%s\n",rem.c_str());

  infile.clear();
  return rem;
}

string Cluster::ReadQuantumEspressoSpeciesSection(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,17).c_str(),line.substr(0,18).c_str());
    if ( line.substr(0,17)=="$quantum_espresso" || line.substr(0,16)=="$quantumEspresso" || line.substr(0,11)=="$qe_species") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();
  return rem;
}

string Cluster::ReadQuantumEspressoSupercellSection(ifstream& infile) {
  string rem = "";
  string line;
  int length_a, length_b, length_c;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,17).c_str(),line.substr(0,18).c_str());
    if ( line.substr(0,13)=="$qe_supercell") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
          istringstream iss(line);              
          iss >> length_a;
          iss >> length_b;
          iss >> length_c;
	  break;
	}
      }
    }
  }
  
  //printf("length_a = %d, length_b = %d, length_c = %d\n", length_a,length_b,length_c);

  if (length_a > 1 || length_b > 1 || length_c > 1){
    //printf("Initializing Supercell!\n\n");
    Vector size(3);
    size[0] = length_a;
    size[1] = length_b;
    size[2] = length_c;
    Params::Parameters().SetParameter("SUPERCELL_JOB","TRUE");
    Params::Parameters().SetParameter("SUPERCELL_ANALYZE_ONLY","TRUE");
    //size.Print("Supercell size before setting");
    Params::Parameters().SetSupercellSize(size);
    //size = Params::Parameters().GetSupercellSize();
    //size.Print("Supercell size");
  }
  fflush(stdout);
  //rem += "\n\n";
  infile.clear();
  return rem;
}

string Cluster::ReadKPoints(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,17).c_str(),line.substr(0,18).c_str());
    if ( line.substr(0,11)=="$kpoints") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();

  if(rem == "") {  
    printf("Kpoints not set. Assuming gamma point\n");
    rem = "1 1 1 1 1 1\n";
  }

  return rem;
}

string Cluster::ReadHubbardDerivatives(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,9).c_str(),line.substr(0,10).c_str());
    if ( line.substr(0,10)=="$hubbard") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();

  if(rem == "") {  
    printf("Hubbard Derivatives not set. Needed to run DFTB. Exiting...\n");
    exit(0);
    //rem = "1 1 1 1 1 1\n";
  }

  return rem;
}

string Cluster::ReadSlaterKoster(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,8).c_str(),line.substr(0,9).c_str());
    if ( line.substr(0,9)=="$slater") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();

  if(rem == "") {  
    printf("Slater Koster files not set. Needed to run DFTB. Exiting...\n");
    exit(0);
    //rem = "1 1 1 1 1 1\n";
  }

  return rem;
}

string Cluster::ReadMaxAngularMomentum(ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("this: %s, or this: %s\n",line.substr(0,13).c_str(),line.substr(0,14).c_str());
    if ( line.substr(0,12)=="$max_angular") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "\n\n";
  infile.clear();

  if(rem == "") {  
    printf("Max Angular Momentum not set. Needed to run DFTB. Exiting...\n");
    exit(0);
    //rem = "1 1 1 1 1 1\n";
  }

  return rem;
}

// Read the Tinker $tinker section
string Cluster::ReadTinkerRemSection(ifstream& infile) {

  // The $tinker section gets stored as one long string, "tinker"
  string rem = "\n"; 
  string line;
  
  Rewind(infile); // Rewind the file, just in case
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,7)=="$tinker") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	}
	else
	  break;
      }
    }
  }
  //rem += "\n";

  infile.clear();
  return rem; 
  
}

// Orient Parameters --- by Ali
void Cluster::ReadAIFFParameters(ifstream& infile) {
  string line;
  // Rewind the file, just in case
  Rewind(infile);
  
  string rem_aiff;
  
  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,5)=="$aiff" || line.substr(0,7)=="$orient") {
      rem_aiff += line + "\n";
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) == "$end") {
	  rem_aiff += line + "\n";
	  break;
	}
	else if ( !line.empty() ) {
	  
	  istringstream iss(line);
	  int i = 0;
	  string parameter, value, read;
	  iss >> parameter;
	  iss >> value;
	  if (value =="=")
	    iss >> value;
	  
	  //printf("Parameter = %s, value = %s\n",parameter.c_str(),
	  // value.c_str());
	  
	  // Test for success?
	  
	  // Set the parameter value
	  Params::Parameters().SetParameter(parameter,value);
	  
	  rem_aiff += line + "\n";
	}
      }
    }
  }
  infile.clear();
  
  string param = "AIFF_REM";
  Params::Parameters().SetParameter(param,rem_aiff);
  
  // Optionally print out the parameters
  if (Params::Parameters().PrintLevel() > 0) Params::Parameters().Print();
  
}

// Read the $aiff section  --- by Ali
string Cluster::ReadAIFFRemSection(ifstream& infile) {
  
  // The $aiff section gets stored as one long string, "rem"
  string rem = "\n"; 
  string line;
  
  Rewind(infile); // Rewind the file, just in case
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,5)=="$aiff" || line.substr(0,7)=="$orient") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	}
	else
	  break;
      }
    }
  }
  //rem += "\n";
  
  
  infile.clear();
  return rem; 
  
}

string Cluster::ReadCrystalHeadingSection(ifstream& infile){

  // The $crystal_head section gets stored as one long string, "hear"
  string heading = ""; 
  string line;

  Rewind(infile); // Rewind the file, just in case
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,16) == "$crystal_heading"){
      for(int i=0;;i++){
	getline(infile,line);
	if(line.substr(0,4) != "$end"){
	  heading += line;
	  heading += "\n";
	}
	else
	  break;
      }
    }
  }
  infile.clear();
  //printf("\n%s\n",heading.c_str());
  return heading;
}

string Cluster::ReadCrystalBasisSection(ifstream& infile){

  // The $crystal_head section gets stored as one long string, "hear"
  string basis = ""; 
  string line;

  Rewind(infile); // Rewind the file, just in case
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,14) == "$crystal_basis"){
      for(int i=0;;i++){
	getline(infile,line);
	if(line.substr(0,4) != "$end"){
	  basis += line;
	  basis += "\n";
	}
	else
	  break;
      }
    }
  }
  infile.clear();
  //printf("\n%s\n",basis.c_str());
  //exit(0);
  return basis;
}



string Cluster::ReadCrystalEndingSection(ifstream& infile){

  // The $crystal_head section gets stored as one long string, "hear"
  string ending = ""; 
  string line;

  Rewind(infile); // Rewind the file, just in case
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,15) == "$crystal_ending"){
 
      for(int i=0;;i++){
	getline(infile,line);
	if(line.substr(0,4) != "$end"){
	  ending += line;
	  ending += "\n";
	}
	else
	  break;
      }
    }
  }
  infile.clear();
  //printf("\n%s\n",ending.c_str());
  return ending;
}

// Print out the paths to any external programs we are using
void Cluster::PrintExternalPrograms() {

  printf("The following external software packages will be used:\n");
  
  // First QM
  if (Params::Parameters().GetQMType() == 1 ) { // Qchem
    printf("   Q-Chem: %s\n",getenv("QC"));
  }
  else if(Params::Parameters().GetQMType() == 2 ) { // MolPro
    string cmd = "which molpro";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {printf("ERROR reading MolPro path!\n"); exit(1);}
    char buffer[256];
    string MolPro_path = "";
    while (!feof(pipe)) {
      if (fgets(buffer, 256, pipe) != NULL)
	MolPro_path += buffer;
    }
    pclose(pipe);
    printf("   MolPro: %s\n",MolPro_path.c_str());
  }
  else if(Params::Parameters().GetQMType() == 7 ) { //PSI4
    string cmd = "which psi4";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {printf("ERROR reading PSI4 path!\n"); exit(1);}
    char buffer[256];
    string PSI4_path = "";
    while (!feof(pipe)) {
      if (fgets(buffer, 256, pipe) != NULL)
        PSI4_path += buffer;
      }
      pclose(pipe);
      printf("    Psi4: %s\n",PSI4_path.c_str());
  }

  // Now MM
  if (Params::Parameters().GetMMType() == 1 ) { // Tinker
    string cmd = "which analyze";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {printf("ERROR reading tinker path!\n"); exit(1);}
    char buffer[256];
    string tinker_path = "";
    while (!feof(pipe)) {
      if (fgets(buffer, 256, pipe) != NULL)
	tinker_path += buffer;
    }
    pclose(pipe);
    
    //string tinker_path = system(cmd.c_str());
    printf("   Tinker: %s\n",tinker_path.c_str());
  }
  
  if (Params::Parameters().GetMMType()==2) { // AIFF
    printf("  CamCASP: %s\n",Params::Parameters().GetCamCaspHome().c_str());
  }

  printf("\n");

}


// Read the embedding charges section
Vector Cluster::ReadEmbeddingCharges(ifstream& infile) {
  
  printf("Reading Embedding charges.\n");

  // do some initialization
  int atom = 0; // atom counter
  int index = 0; // monomer counter
  bool chg_sxn = false;
  string line;

  int Ntot = GetTotalNumberOfAtoms();
  Vector embed_chg(Ntot);
  embed_chg.Set();

  Rewind(infile); // Rewind the file, just in case
 

  // Begin reading the file
  while ( !infile.eof() ) {
    getline(infile,line);

    // Section is marked by "$embedding_charges"
    if (line.substr(0,18)=="$embedding_charges") {
      chg_sxn = true;
    }

    if (line.substr(0,2)=="--" && chg_sxn==true) {
      index++; // increment monomer counter
      
      for (int i=0;i<GetNumberOfAtoms(index);i++)  {
	getline(infile,line);
	istringstream iss(line);
	iss >> embed_chg[atom];
	//printf("Fragment %d, Line %d = %s\n",index, atom,line.c_str());

	atom++; // increment our atom counter
      }
    }

    // If we have reached the end of the section, break
    if (line.substr(0,4)=="$end" && chg_sxn==true) {
      break;
    }
  }

    if (atom != Ntot) {
      printf("Error: Cluster::ReadEmbeddingCharges(). Error parsing charges: Natoms = %d, Ntot = %d\n",atom,Ntot);
      exit(1);
    }

  infile.clear();

  return embed_chg;
}


// Rewinds the file
void Cluster::Rewind(ifstream& infile) {
  infile.clear();
  infile.seekg(0,ios::beg);
}


// Create arrays to store atomic numbers and symbols
void Cluster::SetAtomicSymbols() {
  AtomicSymbols = new string[GetTotalNumberOfAtoms()];
  AtomicNumbers = new int[GetTotalNumberOfAtoms()];
  AtomicMasses = new double[GetTotalNumberOfAtoms()];

  int index = 0;
  for (int imon=1;imon<=NMon;imon++) {
    int N = GetNumberOfAtoms(imon);
    for (int j=0;j<N;j++) {
      AtomicSymbols[index] = Monomers[imon].GetSymbol(j);
      AtomicNumbers[index] = Monomers[imon].GetAtomicNumber(j);
      AtomicMasses[index] = Monomers[imon].GetAtom(j).GetAtomicMass();
      index++;
    }
  }
}

// Reads xyz coordinates from monomers to create a single vector
void Cluster::ReadCurrentCoordinates() {

  int Ntot = GetTotalNumberOfAtoms();

  // Initialize the vector
  //AtomicCoordinates.Initialize( 3 * Ntot );
  // moved initialization to main routine

  int i = 0;

  // Loop over monomers
  for (int imon=1; imon<= NMon; imon++) {
    int N = GetNumberOfAtoms(imon);
    
    // Loop over atoms in each monomer
    for (int iatom=0; iatom<N; iatom++) {

      // Loop over xyz
      for (int dim=0; dim<3; dim++) {

	AtomicCoordinates[i] = Monomers[imon].GetAtom(iatom).GetPosition(dim);
	//printf("AtomicCoordinates[%d] = %f\n",i,AtomicCoordinates[i]);
	i++;
      }
    }
  }   

}


//check to make sure that the dimer rotation matrix is in the space group by checking it against the other dimer rotation matrix
void Cluster::CorrectingDimerRotations(int thisDim,bool image, bool NonImagetoImage,int ListNumb){

  //Dimer& ThisDim;
  Matrix Rot;
  if(image){
    //ThisDim = DimerImages[iDim];
    Rot = DimerImages[thisDim].GetRotationList()[ListNumb]; 
  }
  else{
    //ThisDim = Dimers[iDim];
    if(!NonImagetoImage)
      Rot = Dimers[thisDim].GetRotationList()[ListNumb]; 
    else
      Rot = Dimers[thisDim].GetPeriodicRotationList()[ListNumb]; 
  }
  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance

  bool RotCorrect = false;
  //printf("dimer = %i, ListNumb = %i\n",thisDim,ListNumb);
  //Rot.Print("Rot");

  //first check to see if Rot used in any of the dimers
  for(int iDim = 1; iDim <= NDim;iDim++){

    int ListSize = (Dimers[iDim].GetRotationList().size() + Dimers[iDim].GetPeriodicRotationList().size());

    bool skip = 0;
    if(!image && thisDim == iDim)
      skip = 1;
    //skip if dimer is planar or linear top because their rotation matrix may be wrong
    else if(Dimers[iDim].IsDimerPlanarOrLinear()){
      skip = 1;
      //printf("skip d(%i,%i)\n",Dimers[iDim].GetIndexA(),Dimers[iDim].GetIndexB());
    }


    if(!skip){
      //printf("thisDim = %i iDim = %i, image = %i no peri\n",thisDim,iDim,image);

      for(int iList = 0; iList<ListSize;iList++){
	
	Matrix diffRot;
	if(iList < Dimers[iDim].GetRotationList().size()){
	  diffRot = Dimers[iDim].GetRotationList()[iList];
	  
	}
	else{
	  diffRot = Dimers[iDim].GetPeriodicRotationList()[iList - Dimers[iDim].GetRotationList().size()];
	}      


	//diffRot.Print("Dimer Rot");
	
	//subtraction matrix
	double SumDiff = 0.0;
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++){
	    diffRot(i,j) -= Rot(i,j);
	    SumDiff += fabs(diffRot(i,j));
	  }
	//diffRot.Print("diffRot");	
	// printf("SumDiff= %f\n",SumDiff);
	//match found
	if(fabs(SumDiff) < tolerance){

	  RotCorrect = true;	  
	  return;
	}
      }//end of if statement 
    }//loop over symlist
  }//end loop over iDim

  if(!RotCorrect){
  //first check to see if Rot used in any of the image dimers
    for(int iDim = 1; iDim <= NDim_images;iDim++){
      
      int ListSize = (DimerImages[iDim].GetRotationList().size());
    

      
      bool skip = 0;
      if(image && thisDim == iDim)
	skip = 1;
      //skip if dimer is non an asymmetry top because their rotation matrix may be wrong
      else if(DimerImages[iDim].IsDimerPlanarOrLinear()){
	skip = 1;     
	//printf("skip d(%i,%i)\n",DimerImages[iDim].GetIndexA(),DimerImages[iDim].GetIndexB());
      }

      if(!skip){
	


	 //printf("thisDim = %i iDim = %i, image = %i\n",thisDim,iDim,image);

	for(int iList = 0; iList<ListSize;iList++){
	  
	  Matrix diffRot = DimerImages[iDim].GetRotationList()[iList];
	  // diffRot.Print("Dimer Rot");
	  
	  //subtraction matrix
	  double SumDiff = 0.0;
	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++){
	      diffRot(i,j) -= Rot(i,j);
	      SumDiff += fabs(diffRot(i,j));
	    }
	  //diffRot.Print("diffRot");	
	  // printf("SumDiff= %f\n",SumDiff);
	  //match found
	  if(fabs(SumDiff) < tolerance){
	    //printf("rot found\n\n");
	    RotCorrect = true;	    
	    return;
	  }
	}//end of if statement 
      }//loop over symlist
    }//end loop over iDim

  }
  //printf("RotCorrect = %i\n",RotCorrect);

  //printf("rot not found for d(%i,%i)\n\n",DimerImages[thisDim].GetIndexA(),DimerImages[thisDim].GetIndexB());
    //exit(0);



  //if(image)
  //  printf("d(%i,%i)",DimerImages[thisDim].GetIndexA(),DimerImages[thisDim].GetIndexB());
  //else
  //  printf("d(%i,%i)",Dimers[thisDim].GetIndexA(),Dimers[thisDim].GetIndexB());
  //Rot.Print("Rot");

  //no match found; replace rotation matrix
  if(!RotCorrect){
    
    
    //check against dimers first
    for(int iDim = 1; iDim <= NDim;iDim++){
      
      int ListSize = (Dimers[iDim].GetRotationList().size() + Dimers[iDim].GetPeriodicRotationList().size());
      
      //printf("image = %i thisDim = %i iDim = %i, NDim = %i\n",image,thisDim,iDim,NDim);
      
      bool skip = 0;
      if(!image && thisDim == iDim)
	skip = 1;
      //skip if dimer is planar or linear top because their rotation matrix may be wrong
      else if(Dimers[iDim].IsDimerPlanarOrLinear()){
	skip = 1;
      //printf("skip d(%i,%i)\n",Dimers[iDim].GetIndexA(),Dimers[iDim].GetIndexB());
      }
      
      if(!skip){
	for(int iList = 0; iList<ListSize;iList++){
	  Matrix Rot2;
          if(iList < Dimers[iDim].GetRotationList().size())
	    Rot2 = Dimers[iDim].GetRotationList()[iList];
          else{
            int j = iList - Dimers[iDim].GetRotationList().size();
            Rot2 = Dimers[iDim].GetPeriodicRotationList()[j];
          }
          //if(Dimers[iDim].GetIndexA() == 3 && Dimers[iDim].GetIndexB() == 6)
	  //printf("\n\ncompare d(%i,%i)\n",Dimers[iDim].GetIndexA(),Dimers[iDim].GetIndexB());
	  //Rot2.Print("Rot2");
	  if(CheckDimerRotationOperator(thisDim,Rot2,Rot,image)){
	      
	    //printf("MatchFound\n");
	    //Rot2.Print("Rot2");
            RedetermineDimerSymmetry(thisDim,Rot2,image,NonImagetoImage,ListNumb);     
            //exit(0);
	    return;
	  }
	}
      }  
    }  
  
    //first check to see if Rot used in any of the image dimers
    for(int iDim = 1; iDim <= NDim_images;iDim++){

      bool skip = 0;
      if(image && thisDim == iDim){
	skip = 1;
	//printf("skip d(%i,%i)\n",DimerImages[iDim].GetIndexA(),DimerImages[iDim].GetIndexB());
      //skip if dimer is non an asymmetry top because their rotation matrix may be wrong
      }else if(DimerImages[iDim].IsDimerPlanarOrLinear()){
	skip = 1;     
	//printf("skip d(%i,%i)\n",DimerImages[iDim].GetIndexA(),DimerImages[iDim].GetIndexB());
      }
      if(!skip){

	int ListSize = DimerImages[iDim].GetRotationList().size();
      
	for(int iList = 0; iList<ListSize;iList++){
	  Matrix Rot2 = DimerImages[iDim].GetRotationList()[iList];
          //printf("iList = %i",iList);

          //if(DimerImages[iDim].GetIndexA() == 3 && DimerImages[iDim].GetIndexB() == 6)
	  //printf("\n\ncompare d(%i,%i)\n",DimerImages[iDim].GetIndexA(),DimerImages[iDim].GetIndexB());
	  //Rot2.Print("Rot2");
	  if(CheckDimerRotationOperator(thisDim,Rot2,Rot,image)){
	    //printf("MatchFound\n");
	    //Rot2.Print("Rot2");
            RedetermineDimerSymmetry(thisDim,Rot2,image,NonImagetoImage,ListNumb);
            //exit(0);
	    return;
          }
	}	
      }
    }
  }
  if(image)
    printf("WARNING:CorrectingDimerRotations():: Could not verify d(%i,%i) rotation matrix ListNumb = %i\n",
	   DimerImages[thisDim].GetIndexA(),DimerImages[thisDim].GetIndexB(),ListNumb);
  else
    printf("WARNING:CorrectingDimerRotations():: Could not verify d(%i,%i) rotation matrix ListNumb = %i\n",
	   Dimers[thisDim].GetIndexA(),Dimers[thisDim].GetIndexB(),ListNumb);
  printf("If rotational operator is not in space group. Turn of symmetry with \"SPACE_SYMMETRY = false\"\n");
  Rot.Print("Rot");
  printf("\n\n");


}


//if dimer rotation matrix is not in space group, determine which rotation matrix of the other dimers works
bool Cluster::CheckDimerRotationOperator(int iDim, Matrix Rot,Matrix WrongMat,bool image){

  double tolerance = Params::Parameters().GetSymmetryTolerance();
  //WrongMat.Print("WrongMat");
  //Rot.Print("Rot");

  int N;

  int SymCOM;
  Monomer MonA;
  Monomer MonB;
  
  if(!image){
    
   N = Dimers[iDim].GetNumberOfAtoms();
   MonA = Dimers[iDim].GetMonomerA();
   MonB = Dimers[iDim].GetMonomerB();
  }else{

    
    N = DimerImages[iDim].GetNumberOfAtoms();
    MonA = DimerImages[iDim].GetMonomerA();
    MonB = DimerImages[iDim].GetMonomerB();

  }


  //printf("indexA = %i indexB = %i\n",MonA.GetIndex(),MonB.GetIndex());

  //center of mass coordinates
  MonA.FindCenterOfMass();
  MonB.FindCenterOfMass();
  //find the center of mass of dimer
  Vector COM(3);//center of mass vector of dimer
  COM.Set();
  COM[0] = (MonA.GetCenterOfMass(0)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(0)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[1] = (MonA.GetCenterOfMass(1)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(1)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[2] = (MonA.GetCenterOfMass(2)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(2)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  //printf("center of mass is %10.6f  %10.6f  %10.6f\n",COM[0],COM[1],COM[2]);;


  Matrix COMCoord(3,N);
  for(int i=0;i<N;i++) {

    if(i<MonA.GetNumberOfAtoms()){

      COMCoord(0,i) = MonA.GetAtom(i).GetPosition(0)-COM[0];// x center of mass
      COMCoord(1,i) = MonA.GetAtom(i).GetPosition(1)-COM[1];// y center of mass
      COMCoord(2,i) = MonA.GetAtom(i).GetPosition(2)-COM[2];// z center of mass
      
    }else{
      int j = i - MonA.GetNumberOfAtoms();
      COMCoord(0,i) = MonB.GetAtom(j).GetPosition(0)-COM[0];// x center of mass
      COMCoord(1,i) = MonB.GetAtom(j).GetPosition(1)-COM[1];// y center of mass
      COMCoord(2,i) = MonB.GetAtom(j).GetPosition(2)-COM[2];// z center of mass 
    }
  }


  //While the Rotation Matrix is not in the space group, it should reproduce the original dimer that this one is symmetrical to.
  WrongMat.Transpose();
  Matrix SymCoord = WrongMat.Multiply(COMCoord);


  //Applying the rotation matrix
  Matrix RotCoord = Rot.Multiply(COMCoord);

  
  //if(MonA.GetIndex() == 1 && MonB.GetIndex() == 25){
  //  printf("/n d(%i,%i)\n",MonA.GetIndex(),MonB.GetIndex());
  //  COMCoord.Print("COMCoord");
  //  SymCoord.Print("SymCoord");
  //  RotCoord.Print("RotCoord");
  //  Rot.Print("Rot");
  //  WrongMat.Print("Wrong");
  //}
  

  bool CoordMatched = true;
  int i = 0;//index for atoms on this monomer
  while(i < N && CoordMatched){
    int j=0;
   
    string thisAtom;
    if(i < MonA.GetNumberOfAtoms())
      thisAtom = MonA.GetAtom(i).GetSymbol();
    else
      thisAtom = MonB.GetAtom(i-MonA.GetNumberOfAtoms()).GetSymbol();
 
    bool AtomMatched = false;
    while(j < N && !AtomMatched){
      double x_diff = fabs(RotCoord(0,i)-SymCoord(0,j));
      double y_diff = fabs(RotCoord(1,i)-SymCoord(1,j));
      double z_diff = fabs(RotCoord(2,i)-SymCoord(2,j)); 
      bool SameAtomType = false;
 
      string SymAtom;
      if(j < MonA.GetNumberOfAtoms())
       SymAtom = MonA.GetAtom(j).GetSymbol();
      else
       SymAtom = MonB.GetAtom(j-MonA.GetNumberOfAtoms()).GetSymbol();

      if(thisAtom == SymAtom) {
	SameAtomType = true;
      }
      if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance && SameAtomType){
	AtomMatched = true;
	//printf("atom %i = %s atom %i = %s\n",i,thisAtom.c_str(),j,SymAtom.c_str());
	//printf(" i %f %f %f \n",COMCoord(0,i),COMCoord(1,i),COMCoord(2,i));
	//printf(" j %f %f %f \n\n",SymCoord(0,j),SymCoord(1,j),SymCoord(2,j));   
      }
      
      if(j == N-1 && !AtomMatched)
	CoordMatched = false;   
      j++;
    }
    i++;
  }

  //printf("\n\n");
  return CoordMatched;

}

//Update all dimer symmetry information with new Rotation matrix
void Cluster::RedetermineDimerSymmetry(int thisDim,Matrix Rot,bool image,bool NonImagetoImage,int ListNumb){


  /*
  if(image)
    DimerImages[thisDim].AlterRotationList(Rot,ListNumb);
  else if(NonImagetoImage)
    Dimers[thisDim].AlterPeriodicRotationList(Rot,ListNumb);
  else
    Dimers[thisDim].AlterRotationList(Rot,ListNumb);
  */


  //recreate old dimer to determine which atoms in the dimers are symmetrical equivalent
  int SymA;
  int SymB;
  int NA;
  int NB;
  int k_point[3];
  //int MonOutsideCell;//if -1, neither monomer is outside cell.
                      //if 0, monomer A is outside cell
                      //if 1, monomer B is outside cell
  Vector SymMasses;

  //Getting varibles needed to recreate old dimer
  if(image){
    SymA = DimerImages[thisDim].GetSymmetryList()[2*ListNumb];
    SymB = DimerImages[thisDim].GetSymmetryList()[2*ListNumb+1];

    k_point[0] = DimerImages[thisDim].GetSymmetricalImageCell()[3*ListNumb];
    k_point[1] = DimerImages[thisDim].GetSymmetricalImageCell()[3*ListNumb+1];
    k_point[2] = DimerImages[thisDim].GetSymmetricalImageCell()[3*ListNumb+2];

    NA = Monomers[SymA].GetNumberOfAtoms();
    NB = Monomers[SymB].GetNumberOfAtoms();

  }else if(NonImagetoImage){

    SymA = Dimers[thisDim].GetPeriodicSymmetryList()[2*ListNumb];
    SymB = Dimers[thisDim].GetPeriodicSymmetryList()[2*ListNumb+1];

    k_point[0] = Dimers[thisDim].GetSymmetricalImageCell()[3*ListNumb];
    k_point[1] = Dimers[thisDim].GetSymmetricalImageCell()[3*ListNumb+1];
    k_point[2] = Dimers[thisDim].GetSymmetricalImageCell()[3*ListNumb+2];

    NA = Monomers[SymA].GetNumberOfAtoms();
    NB = Monomers[SymB].GetNumberOfAtoms();

  }else{
    SymA = Dimers[thisDim].GetSymmetryList()[2*ListNumb];
    SymB = Dimers[thisDim].GetSymmetryList()[2*ListNumb+1];

    NA = Monomers[SymA].GetNumberOfAtoms();
    NB = Monomers[SymB].GetNumberOfAtoms();

    k_point[0] = 0;
    k_point[1] = 0;
    k_point[2] = 0;

  }

  // Determine the shift from the central cell to the image cell
  Vector shift(3);
  for (int i=0;i<3;i++) {
     shift[i] = k_point[0]*unit_cell[0][i] + k_point[1]*unit_cell[1][i] + 
        k_point[2]*unit_cell[2][i];
  }
 
  Monomer SymMonA = Monomers[SymA];
  Monomer SymMonB = Monomers[SymB];

  if(image || NonImagetoImage){
    Vector new_com(3);
    new_com = SymMonB.GetCenterOfMass();
    new_com += shift;
    SymMonB.Translate(new_com);
    SymMonB.FractionalTranslate(k_point[0],k_point[1],k_point[2]);
  }

  Dimer SymDimer;
  SymDimer.Initialize(SymMonA,SymMonB);

  //SymDimer.PrintDimerCartesian();
  fflush(stdout);

  //SymDimer.PrintDimerCartesian();
  if(image)
    DimerImages[thisDim].AlterSymmetryList(SymDimer,ListNumb,Rot,NonImagetoImage);
  else
    Dimers[thisDim].AlterSymmetryList(SymDimer,ListNumb,Rot,NonImagetoImage);
  //exit(0);
}


//checks to make sure that the monomer matrix is in the space group by matching it against the dimer rotation matrix 
void Cluster::CorrectingMonomerRotations(int iMon){
 
  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance

  bool RotCorrect = false;
  Matrix Rot = Monomers[iMon].GetRotationMatrix();
  

  //Rot.Print("\nRot");
  //first check to see if Rot used in any of the dimers
  for(int iDim = 1; iDim <= NDim;iDim++){
    
    int ListSize = (Dimers[iDim].GetRotationList().size() + Dimers[iDim].GetPeriodicRotationList().size());
    
    for(int iList = 1; iList<ListSize;iList++){

      Matrix diffRot;
      if(iList < Dimers[iDim].GetRotationList().size()){
	diffRot = Dimers[iDim].GetRotationList()[iList];
 
      }
      else{
	diffRot = Dimers[iDim].GetPeriodicRotationList()[iList - Dimers[iDim].GetRotationList().size()];
      }      

       //diffRot.Print("Dimer Rot");

      //subtraction matrix
      double SumDiff = 0.0;
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++){
	  diffRot(i,j) -= Rot(i,j);
	  SumDiff += fabs(diffRot(i,j));
	}
       //diffRot.Print("diffRot");	
      // printf("SumDiff= %f\n",SumDiff);
      //match found
      if(fabs(SumDiff) < tolerance){
        //printf("rot found\n\n");
	RotCorrect = true;
        break;
      }
   
    }//loop over symlist
  }//end loop over iDim

  //no match found replace rotate
  if(!RotCorrect){
    for(int iDim = 1; iDim <= NDim;iDim++){



      //int ListSize = Dimers[iDim].GetRotationList().size() + Dimers[iDim].GetPeriodicRotationList().size();
      int ListSize = Dimers[iDim].GetSymmetryFactor() + Dimers[iDim].GetPeriodicSymmetryFactor();

      for(int iList = 1; iList<ListSize;iList++){
	Matrix Rot2;
	if(iList < Dimers[iDim].GetRotationList().size()){
	  Rot2 = Dimers[iDim].GetRotationList()[iList];
	}
	else
	  Rot2 = Dimers[iDim].GetPeriodicRotationList()[iList - Dimers[iDim].GetRotationList().size()];
	
       //Rot2.Print("Rot2");

	//Check to make sure that Rot2, a rotation operator for a dimer, will rotate the monomer correctly
	if(CheckMonomerRotationOperator(iMon,Rot2)){
          //printf("Sym found %i\n",iMon);
	  Monomers[iMon].SetRotationMatrix(Rot2);
          int Sym = Monomers[iMon].GetSymmetricalMonomer();
          Monomers[iMon].FindSymmetricalEquivalentAtoms(Monomers[Sym]);
	  return;
	}
	
      }
    }
    printf("Warning::Cluster::CorrectingMonomerRotations() m%i rotation matrix not found among dimers\n",iMon);
    Rot.Print("Rot\n");
    //exit(0);

  }


     //if(Monomers[iMon].GetIndex() == 4)
     // exit(0);
  //exit(0);
  
}

//Checks if matrix Rot rotates monomer iMon into it's symmetrical equivalent monomer
bool Cluster::CheckMonomerRotationOperator(int iMon, Matrix Rot){
 
  int N = Monomers[iMon].GetNumberOfAtoms();
  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance
  
  Matrix COM(3,N);//center of mass coordinates

  int Sym = Monomers[iMon].GetSymmetricalMonomer();
  Matrix SymCOM(3,N);
  for(int i=0;i<N;i++) {
    COM(0,i) = Monomers[iMon].GetAtom(i).GetPosition(0)-Monomers[iMon].GetCenterOfMass(0);// x center of mass
    COM(1,i) = Monomers[iMon].GetAtom(i).GetPosition(1)-Monomers[iMon].GetCenterOfMass(1);// y center of mass
    COM(2,i) = Monomers[iMon].GetAtom(i).GetPosition(2)-Monomers[iMon].GetCenterOfMass(2);// z center of mass

    SymCOM(0,i) = Monomers[Sym].GetAtom(i).GetPosition(0)-Monomers[Sym].GetCenterOfMass(0);// x center of mass
    SymCOM(1,i) = Monomers[Sym].GetAtom(i).GetPosition(1)-Monomers[Sym].GetCenterOfMass(1);// y center of mass
    SymCOM(2,i) = Monomers[Sym].GetAtom(i).GetPosition(2)-Monomers[Sym].GetCenterOfMass(2);// z center of mass

  }

  Matrix RotCOM = Rot.Multiply(COM);
  
  //COM.Print("COM");
  //COM.Print("Rot");
  //COM.Print("RotCOM");
  //SymCOM.Print("SymCOM");

  bool CoordMatched = true;
  int i = 0;//index for atoms on this monomer
  while(i < N && CoordMatched){
    int j=0;
    bool AtomMatched = false;
    while(j < N && !AtomMatched){
      double x_diff = fabs(RotCOM(0,i)-SymCOM(0,j));
      double y_diff = fabs(RotCOM(1,i)-SymCOM(1,j));
      double z_diff = fabs(RotCOM(2,i)-SymCOM(2,j)); 
      bool SameAtomType = false;
      if(Monomers[iMon].GetSymbol(i) == Monomers[Sym].GetSymbol(j)) {
	SameAtomType = true;
      }
      if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance && SameAtomType){
	AtomMatched = true;
	//printf("atom %i = %s atom %i = %s\n",i,GetSymbol(i).c_str(),j,GetSymbol(j).c_str());
	//printf(" i %f %f %f \n",Coord1(0,i),Coord1(1,i),Coord1(2,i));
	//printf(" j %f %f %f \n\n",Coord2(0,j),Coord2(0,j),Coord2(2,j));   
      }
      
      if(j == N-1 && !AtomMatched)
	CoordMatched = false;   
      j++;
    }
    i++;
  }
  return CoordMatched;
}

//Set the fractional coordinate for each atom
void Cluster::SetFractionalCoordinates(){

  //getting Cartesian Coordinates
  Vector Coords = GetCurrentCoordinates();
  double tolerance = Params::Parameters().GetSymmetryTolerance();
  //Convert from cartesian to fractional coordinate
  FractionalCoordinates = ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);

  //FractionalCoordinates.PrintGradient("Complete Fractional Coordinates");

  //Setting Fractional Coordinates for each atom
  int i= 0;
  for(int imon=1;imon<=NMon;imon++)
    for(int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++){
	//Setting Fractional Coordinates for each atom
	Vector AtomCoords(3);
	AtomCoords[0] = FractionalCoordinates[3*i];
	AtomCoords[1] = FractionalCoordinates[3*i+1];
	AtomCoords[2] = FractionalCoordinates[3*i+2];
	Monomers[imon].GetAtom(iatom).SetFractionalPosition(AtomCoords);
	i++;
   
    }

}

Vector Cluster::GetFractionalCoordinates(bool IncludeLatParams){

  Vector ReturnedCoords;
    
  if(IncludeLatParams)//Include the lattice Params
    ReturnedCoords = FractionalCoordinates;
  else{//Do not include the lattice params
    ReturnedCoords.Initialize(3*GetTotalNumberOfAtoms());
    for(int i=0;i<3*GetTotalNumberOfAtoms();i++)
      ReturnedCoords[i] = FractionalCoordinates[i];
  }

  return ReturnedCoords;

}


//Reads the coordinates of the atoms on the symmetrical unique atoms or atoms in the asymmetric unit
//Reads the coordinates of the atoms on the symmetrical unique monomers
void Cluster::ReadSymmetryUniqueCoordinates(){ 


  //if( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2))
  //     SymmetryUniqueCoordinates.Initialize(3 * UniqueAtoms + 6); 
  //else
  SymmetryUniqueCoordinates.Initialize(3 * UniqueAtoms);
  int i = 0;

  //Loop over monomers
  for (int imon=1; imon<= NMon; imon++){

    //only count atoms in monomers that's symmetry factor is not zero
    //if(Monomers[imon].GetSymmetryFactor()!=0){
    int N = GetNumberOfAtoms(imon);
    //printf("Monomer %i\n",Monomers[imon].GetIndex());

    //Loop over atoms
    for (int iatom=0; iatom<N; iatom++){
       
      //check if atom is symmetrical unique.
      if(Monomers[imon].GetAtom(iatom).InAsymmetricUnit()){

        //Loop over xyz
        //printf("%s ",Monomers[imon].GetSymbol(iatom).c_str());
        for(int dim=0; dim<3; dim++){
	  SymmetryUniqueCoordinates[i] = Monomers[imon].GetAtom(iatom).GetPosition(dim);
	  //printf("imon = %i iatom = %i %f \n",imon,iatom,SymmetryUniqueCoordinates[i])
	  i++;
        }//end of loop over dim
      //printf("\n");
      }//end of if statement to get symmetrical uniqueness
    }//end of loop over atoms
  }//end of loop over monomers

  //for(int j=0;i<3;i++){
  //   SymmetryUniqueCoordinates[3*UniqueAtoms+i] = UnitCellAxes[i];
  //  SymmetryUniqueCoordinates[3*UniqueAtoms+3+i] = UnitCellAngles[i];
  // }
  //}

  //SymmetryUniqueCoordinates.PrintGradient("SymmetryUniqueCoordinates");

}

Vector Cluster::GetSymmetryUniqueCoordinates (bool IncludeLattice){
  
    if( (((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly()) && IncludeLattice){
      Vector SymmetryUniqueCoordinatesPlusLattice(3*UniqueAtoms+6);
      for(int i=0;i<3*UniqueAtoms;i++){
	SymmetryUniqueCoordinatesPlusLattice[i] = SymmetryUniqueCoordinates[i];
      }
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms]=UnitCellAxes[0];
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms+1]=UnitCellAxes[1];
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms+2]=UnitCellAxes[2];
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms+3]=UnitCellAngles[0];
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms+4]=UnitCellAngles[1];
      SymmetryUniqueCoordinatesPlusLattice[3*UniqueAtoms+5]=UnitCellAngles[2];

      //SymmetryUniqueCoordinatesPlusLattice.PrintGradient("SymmetryUniqueCoordinatesPlusLattice");

      return SymmetryUniqueCoordinatesPlusLattice;
    }
    else 
      return SymmetryUniqueCoordinates;
}


//This function determines the rotation matrix in fractional coordinates
void Cluster::DetermineMonomerFractionalSymmetry(){



  Matrix cell(3,3);
  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];
  double alpha = UnitCellAngles[0];
  double beta = UnitCellAngles[1];
  double gamma = UnitCellAngles[2];
  alpha *= DegreesToRadians;
  beta *= DegreesToRadians;
  gamma *= DegreesToRadians;
  double beta_term = ( cos(alpha) - cos(beta)*cos(gamma) ) / sin(gamma) ;
  double gamma_term = sqrt( 1 - cos(beta)*cos(beta) - beta_term*beta_term ) ;
  
  cell(0,0) = a;
  cell(0,1) = b*cos(gamma);
  cell(0,2) = c*cos(beta);  
  cell(1,0)= 0.0;
  cell(1,1) = b*sin(gamma);       
  cell(1,2) = c*beta_term;                        
  cell(2,0) = 0.0;
  cell(2,1) = 0.0;
  cell(2,2) = c*gamma_term;
  
  Matrix cell_inverse = cell;
  cell_inverse.Inverse();
  
  //cell.Print("cell");
  //cell_inverse.Print("cell_inverse");
  
  //printf("m%i\n",imon);
  
  
  //printf("Determine fractional symmetry\n");
  for(int imon=1;imon<=NMon;imon++){
    //if(Monomers[imon].GetSymmetryFactor()==0){
      
    //set fractional rotation matrix for monomers to the symmetrical unique monomer
    Matrix Rot = Monomers[imon].GetRotationMatrix();
    Matrix FractRot = Rot.Multiply(cell);
    FractRot = cell_inverse.Multiply(FractRot);
    Monomers[imon].SetFractRotationMatrix(FractRot);
    
      //set fractional rotation matrix for the atom to the atoms in the asymetrical unit
    for(int iatom = 0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++){
      Rot = Monomers[imon].GetAtom(iatom).GetRotationMatrix();
      FractRot = Rot.Multiply(cell);
      FractRot = cell_inverse.Multiply(FractRot);
      Monomers[imon].GetAtom(iatom).SetFractRotationMatrix(FractRot);



  
    }
    
    //Finding the translational component of the space group symmetry operator
    int Sym_Mon = Monomers[imon].GetSymmetricalMonomer();
    Monomers[imon].SetFractTranslationVector(Monomers[Sym_Mon]);
    
    //}
  }     

  

}

//Function never completed. Do not use
Matrix Cluster::GetChangeInRotationFromLatticeVectors(int Mon,int ChangeValue,bool FromUniqueToThis){
  
  double delta = 0.001;
  
  //if change value = 0, find dR/dx1
  //if change value = 1, find dR/dy1
  //if change value = 2, find dR/dz1
  //if change value = 3, find dR/dx2
  //if change value = 4, find dR/dy2
  //if change value = 5, find dR/dz2
  //if change value = 6, find dR/dx3
  //if change value = 7, find dR/dy3
  //if change value = 8, find dR/dz3

  Matrix Lattice(3,3);
  Lattice.SetColumnVector(unit_cell[0],0);
  Lattice.SetColumnVector(unit_cell[1],1);
  Lattice.SetColumnVector(unit_cell[2],2);


  //Lattice.Print("Lattice 1");
    
  //changing lattice
    if(ChangeValue==0){//find dR/dx1
      Lattice(0,0) -= delta;
    }
    else if(ChangeValue==1){//find dR/dy1
      Lattice(1,0) -= delta;
    }
    else if(ChangeValue==2){//find dR/dz1
      Lattice(2,0) -= delta;
    }else if(ChangeValue==3){//find dR/dx2
      Lattice(0,1) -= delta;
    }else if(ChangeValue==4){//find dR/dy2
      Lattice(1,1) -= delta;
    }else if(ChangeValue==5){//find dR/dz2
      Lattice(2,1) -= delta;
    }else if(ChangeValue==6){//find dR/dx3
      Lattice(0,2) -= delta;
    }else if(ChangeValue==7){//find dR/dy3
      Lattice(1,2) -= delta;
    }else if(ChangeValue==8){//find dR/dz3
      Lattice(2,2) -= delta;
    }else{
    printf("ERROR::Cluster::GetChangeInRotationFromLatticeVectors:: ChangeValue unknown\n");
    printf("ChangeValue = %i\n",ChangeValue);
    exit(0);
  }

  Matrix InverseLattice = Lattice;
  InverseLattice.Inverse();
  Matrix dRot = Monomers[Mon].GetFractionalRotationMatrix();
    if(FromUniqueToThis)
      dRot.Inverse();
    dRot = dRot.Multiply(InverseLattice);
    dRot = Lattice.Multiply(dRot);

    //Lattice.Print("Lattice 2");
    // InverseLattice.Print("InverseLattice");

    //dRot.Print("dRot");


  //changing lattice
    if(ChangeValue==0){//find dR/dx1
      Lattice(0,0) -= delta;
    }
    else if(ChangeValue==1){//find dR/dy1
      Lattice(1,0) += 2*delta;
    }
    else if(ChangeValue==2){//find dR/dz1
      Lattice(2,0) += 2*delta;
    }else if(ChangeValue==3){//find dR/dx2
      Lattice(0,1) += 2*delta;
    }else if(ChangeValue==4){//find dR/dy2
      Lattice(1,1) += 2*delta;
    }else if(ChangeValue==5){//find dR/dz2
      Lattice(2,1) += 2*delta;
    }else if(ChangeValue==6){//find dR/dx3
      Lattice(0,2) += 2*delta;
    }else if(ChangeValue==7){//find dR/dy3
      Lattice(1,2) += 2*delta;
    }else if(ChangeValue==8){//find dR/dz3
      Lattice(2,2) += 2*delta;
    }else{
    printf("ERROR::Cluster::GetChangeInRotationFromLatticeVectors:: ChangeValue unknown\n");
    printf("ChangeValue = %i\n",ChangeValue);
    exit(0);
  }

  InverseLattice = Lattice;
  InverseLattice.Inverse();
  Matrix dRot2 = Monomers[Mon].GetFractionalRotationMatrix();
    if(FromUniqueToThis)
      dRot2.Inverse();
    dRot2 = dRot2.Multiply(InverseLattice);
    dRot2 = Lattice.Multiply(dRot2);

    //Lattice.Print("Lattice 3");
    //InverseLattice.Print("InverseLattice 3");

    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++){
	//printf("dRot2 = %f \n dRot = %f\n",dRot2(i,j),dRot(i,j));
	dRot2(i,j) -= dRot(i,j);
      }
    dRot2.Scale(1/(2*delta));
    dRot2.Print("dRot2 2");

    return dRot2;


}

//Reads the fractional coordinates of the atoms on the symmetrical unique monomers
void Cluster::ReadSymmetryUniqueFractionalCoordinates(){

  //if( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2))
  //     SymmetricalFractCoord.Initialize(3 * UniqueAtoms + 6); 
  //else
  SymmetricalFractCoord.Initialize(3 * UniqueAtoms);
  // printf("SymmetricalFractCoord\n");
  int i = 0;



  //Loop over monomers
  for (int imon=1; imon<= NMon; imon++){
    //only count atoms in monomers that's symmetry factor is not zero
    //if(Monomers[imon].GetSymmetryFactor()!=0){
    int N = GetNumberOfAtoms(imon);
      //printf("Monomer %i\n",Monomers[imon].GetIndex());


    //Loop over atoms
    for (int iatom=0; iatom<N; iatom++){
      //check if atom is symmetrical unique.
      if(Monomers[imon].GetAtom(iatom).InAsymmetricUnit()){
	
       	//for(int dim=0; dim<3; dim++){
	//Loop over xyz
	//printf("%s ",Monomers[imon].GetSymbol(iatom).c_str());
	Vector AtomCoord = Monomers[imon].GetAtom(iatom).GetFractionalPosition();
	SymmetricalFractCoord[3*i] = AtomCoord[0];
	SymmetricalFractCoord[3*i+1] = AtomCoord[1];
	SymmetricalFractCoord[3*i+2] = AtomCoord[2];
	i++;
	//}//end of loop over dim
      }//end of if statement for symmetrical uniqueness
    }//end for iatom
  }//end for imon

  if ( Params::Parameters().PrintSymmetryInfo())
    SymmetricalFractCoord.PrintGradient("SymmetricalFractCoord");
  //exit(0);
}

Vector Cluster::GetSymmetryUniqueFractionalCoordinates(bool IncludeLatticeParams){
   
    if(IncludeLatticeParams){
      Vector SymmetryFractCoordPlusLattice(3*UniqueAtoms+6);
      for(int i=0;i<3*UniqueAtoms;i++){
	SymmetryFractCoordPlusLattice[i] = SymmetricalFractCoord[i];
      }
      SymmetryFractCoordPlusLattice[3*UniqueAtoms]=UnitCellAxes[0];
      SymmetryFractCoordPlusLattice[3*UniqueAtoms+1]=UnitCellAxes[1];
      SymmetryFractCoordPlusLattice[3*UniqueAtoms+2]=UnitCellAxes[2];
      SymmetryFractCoordPlusLattice[3*UniqueAtoms+3]=UnitCellAngles[0];
      SymmetryFractCoordPlusLattice[3*UniqueAtoms+4]=UnitCellAngles[1];
      SymmetryFractCoordPlusLattice[3*UniqueAtoms+5]=UnitCellAngles[2];
      
      return SymmetryFractCoordPlusLattice;
    }
    else
      return SymmetricalFractCoord;
}

//Gets the coordinates of the unit cell from the coordinates of the asymmetric unit 
Vector Cluster::GetSymmetryImposedCoordinates(Vector reduced_coord,bool usefract,bool GradIncludeLattice){

  int Natoms = GetTotalNumberOfAtoms();
  //Making key for the reduced_coordinates.
  int i=0;
  int Coord_key[UniqueAtoms];
  //Vector updated_coords = GetCurrentCoordinates();//update coordinates of the complete unit cell
  Vector updated_coords;
  //printf("grad include lattice parameter: %d\n",GradIncludeLattice);
  if(GradIncludeLattice){
    updated_coords.Initialize(3*Natoms+6);
  }
  else {
    updated_coords.Initialize(3*Natoms);
  }


  //Use fractional coordiantes
  if(usefract){
    //making the Coord_key for the Coordinates
    for(int iMon=1;iMon<=NMon;iMon++){
      //if(Monomers[iMon].GetSymmetryFactor()!=0)
      for(int iatom=0; iatom<GetNumberOfAtoms(iMon);iatom++){
	if(Monomers[iMon].GetAtom(iatom).InAsymmetricUnit()){
	  Coord_key[i] = Monomers[iMon].GetAtom(iatom).GetGlobalIndex();
	  i++;
	}
      }
    }
    
    int coord = 0;
    for(int iMon=1;iMon<=NMon;iMon++){
      for(int iatom=0; iatom<GetNumberOfAtoms(iMon);iatom++){
	

	//Get Space Group operator
	Matrix RotMat = Monomers[iMon].GetAtom(iatom).GetFractionalRotationMatrix();
	RotMat.Inverse();
	//RotMat.Print("RotMat");
	Vector Translation = Monomers[iMon].GetAtom(iatom).GetFractTranslationVector(true);

	//Degree atoms have to be shifted to preserve symmetry 
	//while degrees of freedom have the same cartesian coordinates
	//	Vector Shift = Monomers[iMon].GetAtom(iatom).GetShiftVector();

	//GlobalSym used to match entries in Coord_key
	int GlobalSym = Monomers[iMon].GetAtom(iatom).GetSymmetricalAtom();	
	for(int i=0;i<UniqueAtoms;i++){
	  if(Coord_key[i] == GlobalSym){
	    for(int xyz=0;xyz<3;xyz++){
	      updated_coords[coord] += RotMat(xyz,0)*reduced_coord[3*i]+RotMat(xyz,1)*reduced_coord[3*i+1]+RotMat(xyz,2)*reduced_coord[3*i+2];
	      updated_coords[coord] += Translation[xyz];
	      coord++;	      
	    } 

	  }
	}

	//RotMat.Print("RotMat");
	//Translation.Print("Translation");
	//Shift.Print("Shift");
	//reduced_coord.PrintGradient("reduced_coord");
	//updated_coords.PrintGradient("updated_coord");
	//exit(0);
      }


    }
    /*	   




    //find the change in the fractional coordinates for the asymetrical unit
    i = 0;
    Vector old_coord = GetSymmetryUniqueFractionalCoordinates(false);
    for(int iMon=1;iMon<=NMon;iMon++)
      //if(Monomers[iMon].GetSymmetryFactor()!=0)
      for(int iatom=0; iatom<GetNumberOfAtoms(iMon);iatom++)
	if(Monomers[iMon].GetAtom(iatom).InAsymmetricUnit()){
	  Coord_key[i] = Monomers[iMon].GetAtom(iatom).GetGlobalIndex();
	  for(int xyz=0;xyz<3;xyz++)
	    dCoord[3*i+xyz] = reduced_coord[3*i+xyz]-old_coord[3*i+xyz];
	  //printf("dCoord for atom %i = %f %f %f\n",
	  //	 i,dCoord[3*i],dCoord[3*i+1],dCoord[3*i+2]);
	  i++;
	}

    printf("\n");

    if(GradIncludeLattice)
      updated_coords = GetFractionalCoordinates(true);
    else
      updated_coords = GetFractionalCoordinates(false);

    //updating the coordinates of the comple unit cell;
    int coord=0;//counter for the coordinates in update_coords
    for(int iMon=1;iMon<=NMon;iMon++){

      //Matrix RotMat = Monomers[iMon].GetFractionalRotationMatrix();
      //RotMat.Inverse();
      //RotMat.Print("RotMat");
      
      // updated_coords.PrintGradient("update_coords");
      for(int iatom=0;iatom<GetNumberOfAtoms(iMon);iatom++){


	Matrix RotMat = Monomers[iMon].GetAtom(iatom).GetFractionalRotationMatrix();
	RotMat.Inverse();
	//RotMat.Print("RotMat");
	
	//GlobalSym used to match entries in Coord_key
	int GlobalSym = Monomers[iMon].GetAtom(iatom).GetSymmetricalAtom();
	for(int i=0;i<UniqueAtoms;i++){
	  if(Coord_key[i] == GlobalSym){
	    //printf("iMon = %i iatom = %i GlobalSym = %i Coord_key[i]= %i\n",
	    //	 iMon,iatom,GlobalSym,Coord_key[i]); 
	    //rotationg dCoord elements using the rotation matrix
	    for(int xyz=0;xyz<3;xyz++){
	      //printf("iMon = %i iatom = %i coord = %i i=%i\n",iMon,iatom,coord,i);  
	      updated_coords[coord] += RotMat(xyz,0)*dCoord[3*i]+RotMat(xyz,1)*dCoord[3*i+1]+RotMat(xyz,2)*dCoord[3*i+2];
	      coord++;
	    }//end of loop over xyz
	  }//end of if statement Coord_key[i] = GlobalSym
	}//end of loop over i
      }//end of loop over atoms in monomer
    }//end of loop over monomers

    */
    if(coord != 3*Natoms){
      printf("ERROR::Cluster::GetSymmetryImposedCoordinates() Unable to update atomic coordinates\n coord = %i\n",coord);
      exit(1); 
    }
    
  }else{
    //Vector dCoord(3*UniqueAtoms);
    
    //Finding the difference between the updated coordinates for the symmetrical unique coordinates and the old coordinates.
    //This difference will be used to find change in the rest of the system.
    Vector dCoord(3*UniqueAtoms);
    

      Vector old_coord = GetSymmetryUniqueCoordinates();
      if(GradIncludeLattice)
	old_coord = GetSymmetryUniqueCoordinates(false);


    for(int iMon=1;iMon<=NMon;iMon++)
	for(int iatom=0; iatom<GetNumberOfAtoms(iMon);iatom++)
	  if(Monomers[iMon].GetAtom(iatom).InAsymmetricUnit()){
	    Coord_key[i] = Monomers[iMon].GetAtom(iatom).GetGlobalIndex();
	    //printf("iMon =%i iatom =%i i=%i\n",iMon,iatom,i);
	    fflush(stdout);
	    for(int xyz=0;xyz<3;xyz++)
	      dCoord[3*i+xyz] = reduced_coord[3*i+xyz]-old_coord[3*i+xyz];
	    i++;
	  }
    //printf("Coord_key =[");
    //for(int i=0;i<UniqueAtoms;i++)
    //  printf(" %i",Coord_key[i]);
    //printf("]\n");
     
    //dCoord.PrintGradient("dCoord");

    if ( (((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly())
	  && GradIncludeLattice)
      updated_coords.Initialize(3*Natoms+6);
    else
      updated_coords.Initialize(3*Natoms);   
    for(int i=0;i<3*Natoms;i++)
      updated_coords[i] = GetCurrentCoordinates()[i];//update coordinates of the complete unit cell
    
    // updated_coords.PrintGradient("update_coords");
    int coord=0;//counter for the coordinates in update_coords
    
    //updating the coordinates of the complete unit cell;
    for(int iMon=1;iMon<=NMon;iMon++){
      //Matrix RotMat = Monomers[iMon].GetRotationMatrix();
      //RotMat.Transpose();
      for(int iatom=0;iatom<GetNumberOfAtoms(iMon);iatom++){
	Matrix RotMat = Monomers[iMon].GetAtom(iatom).GetRotationMatrix();
	RotMat.Inverse();
	//GlobalSym used to match entries in Coord_key
	int GlobalSym = Monomers[iMon].GetAtom(iatom).GetSymmetricalAtom();
	for(int i=0;i<UniqueAtoms;i++){
	  //printf("iMon = %i iatom = %i GlobalSym = %i Coord_key[i]= %i\n",
	  //       iMon,iatom,GlobalSym,Coord_key[i]); 
	  if(Coord_key[i] == GlobalSym){
	    //rotationg dCoord elements using the rotation matrix
	    for(int xyz=0;xyz<3;xyz++){
	      //printf("iMon = %i iatom = %i coord = %i i=%i\n",iMon,iatom,coord,i);  
	      updated_coords[coord] += RotMat(xyz,0)*dCoord[3*i]+RotMat(xyz,1)*dCoord[3*i+1]+RotMat(xyz,2)*dCoord[3*i+2];
	      coord++;
	    }//end of loop over xyz
	  }//end of if statement Coord_key[i] = GlobalSym
	}//end of loop over i
      }//end of loop over atoms in monomer
    }//end of loop over monomers


    if(coord != 3*Natoms){
      printf("ERROR::Cluster::GetSymmetryImposedCoordinates() Unable to update atomic coordinates\n coord = %i\n",coord);
      exit(1); 
    }
  } //End In Fractional logic test
  
  //including the lattice parameters
  if( ((Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2) || Params::Parameters().UseFullQMOnly()) && GradIncludeLattice ){
      /* && !Params::Parameters().UseCartesianCoordinates()*/
    for(int i = 0;i< 6;i++){
      updated_coords[3*Natoms+i] = reduced_coord[3*UniqueAtoms+i];
      //printf("the updated coord: %f\n",updated_coords[3*Natoms+i]);
    } 
  }
  //updated_coords.PrintGradient("updated_coords");
  //exit(0);
  return updated_coords;
  
}

Vector Cluster::EnergyOfFiniteDifference(double delta, bool IsPositive){

 
  Vector Energy(3*UniqueAtoms);
  Vector Original_Coord = GetCurrentCoordinates();
  Vector Original_SymCoord = GetSymmetryUniqueCoordinates(true);

  //Original_SymCoord.Print("Original_SymCoord");
  /*

  Vector Original_Coord = GetSymmetryUniqueCoordinates(false);
  //Original_Coord.PrintGradient("Original_Coords");
  int = 0;
  while (i<3*UniqueAtoms){
    Vector Unique_Coord =  GetSymmetryUniqueCoordinates();

    
    //changing coords 
    if(IsPositive)
      Unique_Coord[i] += delta;
    else
      Unique_Coord[i] -= delta;

    GetSymmetryUniqueFractionalCoordinates(false);
    

  }
}
  */


  int counter=0;//this is a counter for energy
  for(int iMon=1;iMon<=NMon;iMon++){

    //printf("Monomer %i has symfac of %i\n",iMon,sym_fac);
    int Natoms = Cluster::cluster().GetNumberOfAtoms(iMon);
    for(int iAtom=0;iAtom<Natoms;iAtom++){

      if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	
	for(int xyz=0;xyz<3;xyz++){
	  
	  //reseting coordinates
	  Vector Coord = Original_Coord;
	  Vector SymCoord = Original_SymCoord;

	  int GlobalIndex = Monomers[iMon].GetAtom(iAtom).GetGlobalIndex();  


	  //Skip due to a restriction from symmetry
	  bool skip = 0;
	  if(Monomers[iMon].GetAtom(iAtom).IsAtomFrozen())
	    skip = 1;
	  else if((Monomers[iMon].GetAtom(iAtom).IsXLocked() > 0 || Monomers[iMon].GetAtom(iAtom).FreezeX()) 
		  && xyz == 0 )
	    skip = 1;
	  else if((Monomers[iMon].GetAtom(iAtom).IsYLocked() || Monomers[iMon].GetAtom(iAtom).FreezeY())
		  && xyz == 1)
	    skip = 1;
	  //else if((Monomers[iMon].GetAtom(iAtom).IsZLocked() > 0 || Monomers[iMon].GetAtom(iAtom).FreezeZ())
		  //&& xyz == 2)
	  else if(Monomers[iMon].GetAtom(iAtom).FreezeZ() && xyz == 2)
	    skip = 1;

	  
	  printf("Mon = %i atom = %i xyz = %i\n",iMon,iAtom,xyz);
	  printf("counter = 0\n");
	  printf("skip = %i\n",skip);


	  if(!skip){

	    //plus
	    if(IsPositive){
	      SymCoord[counter] += delta;
	      printf("Shifting atom %i by +%f\n",
		     GlobalIndex,delta);
	    }else{
	      SymCoord[counter] -= delta;
	      printf("Shifting atom %i by -%f\n",
		     GlobalIndex,delta);
	    }	      
	    //Correcting coordiate locked under symmetry
	    Cluster::cluster().MaintainCartesianSymmetry(SymCoord,true);
	    Coord = Cluster::cluster().GetSymmetryImposedCoordinates(SymCoord,false); 
	    //Coord.PrintGradient("Coord");
	    //exit(0);

	    printf("counter = %i\n",counter);
	    Coord.PrintGradient("Coord");
	    
	    //Update the coordinates & get the energy
	    SetNewCoordinates(Coord);
	    RunJobsAndComputeEnergy();
	    Energy[counter]=GetHMBIEnergy();

	    //restoring the coordinates   
	    SetNewCoordinates(Original_Coord);

	  }//end of if statement checking if coord is skipped
	  counter++;
	}//end of xyz loop
      }//end of check if atom is in the asymmetric unit
    }//end of iAtom loop 
  }//end of iMon loop





  Energy.Print("Energy");
  
  return Energy;
}

// Updates all Monomers and Dimers with new coordinates Also erases
// any energies and gradients unless ResetEnergies flag is false.
// Careful: this flag should be true unless you have a good reason for
// it not to be.
void Cluster::SetNewCoordinates(Vector NewCoords, bool ResetEnergies) {


  /*  Deal with Monomers first: */
  int i = 0;
  double xyz[3];
  
  //Lattice Parameters are updated outside function.
  //if ( Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2)
  //  UpdateLatticeParams(NewCoords);  

  // Loop over monomers
  for (int imon=1; imon<= NMon; imon++) {
    int N = GetNumberOfAtoms(imon);

    //printf("Old Cartesian coords\n");
    //Monomers[imon].PrintMonomerCartesian();

    // Loop over atoms in each monomer
    for (int iatom=0; iatom<N; iatom++) {

      // Loop over xyz
      for (int dim=0; dim<3; dim++) {
	xyz[dim] = NewCoords[i];
	i++; 
      }
      //printf("Mon %d, Atom %d, (%f,%f,%f)\n",imon,iatom,xyz[0],xyz[1],xyz[2]);
      Monomers[imon].SetAtomPosition(iatom,xyz);
    }
    // Update local coordinates, but maintain original coordinate systems
    if (Params::Parameters().GetMMType()==2) 
      Monomers[imon].FindLocalCoord(true);

    //printf("New Cartesian coords\n\n");
    //Monomers[imon].PrintMonomerCartesian();
    
    // Update the center of mass of the monomer
    Monomers[imon].FindCenterOfMass();
    
    if (ResetEnergies && Params::Parameters().GetMMType() != 5) 
      Monomers[imon].ResetEnergiesAndGradients();
    
  }  

  /*
  //recheck monomer symmetry
  //should be the same
  int Mon_Number = 0;

  for(int i=1;i<=NMon;i++){
    Monomers[i].SetSymmetryFactor(1);
    if(!Params::Parameters().UseSpaceSymmetry()){
      Mon_Number++;
    }else if(!Monomers[i].SymmetryCheck(Monomers)){
      Mon_Number++;
    }
  }
  if(Mon_Number != UniqueMon){
    printf("ERROR::Cluster::SetNewCoordinates monomer symmetry broken, number of symmetry monomers = %i\n",Mon_Number);
    exit(0);
  }
  */

  /* Update Dimers */

  // Dimers
  for (int i=1;i<=NDim;i++) {
    int m1 = Dimers[i].GetIndexA();
    int m2 = Dimers[i].GetIndexB();
    Dimers[i].ResetSymmetryList();
    // Update dimers with modified monomer data
    //printf("m1 = %d, m2 = %d\n",m1,m2);
    Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);

    if (ResetEnergies) {
      Dimers[i].ResetEnergiesAndGradients();
      if(Params::Parameters().DoFreq())
	 Dimers[i].ResetHessians();
    }

  }

  //check dimer symmetry factor
  //should be the same but in case they are not
  UniqueDim=0;
  for(int i=1;i<=NDim;i++){
    Dimers[i].SetSymmetryFactor(1);
    Dimers[i].SetPeriodicSymmetryFactor(0);
    if(!Params::Parameters().UseSpaceSymmetry())
      UniqueDim++;
    else if(!Dimers[i].SymmetryCheck(Dimers, i))
      UniqueDim++;
  }
  //Reset the current and symmetry unique coordinates in the cluster object
  ReadCurrentCoordinates();
  ReadSymmetryUniqueCoordinates();

  if ( Params::Parameters().IsPeriodic() ) {
    delete [] DimerImages;

    // Identify real space periodic images - local fragments might have changed during optimization
    TotalDim_nosym=NDim;
    CreatePeriodicImageDimerList();
    
    //Reset current and symmetry  unique fractional coordinates in the cluster object
    SetFractionalCoordinates();
    ReadSymmetryUniqueFractionalCoordinates();
  }

  if (Params::Parameters().UseFullQMOnly()) {
    //Reset current and symmetry  unique fractional coordinates in the cluster object
    SetFractionalCoordinates();
    ReadSymmetryUniqueFractionalCoordinates();
  }

  //printing savings due to symmetry
  if(Params::Parameters().UseSpaceSymmetry()){
    printf("\n");
    printf("Applying space group symmetry to reduce the numbers of monomer/dimer calculations required\n");
    printf("There are %i dimers with symmetry and %i without. (Approximate savings factor = %.1f)\n",
	   UniqueDim+NDim_images,TotalDim_nosym, (double) TotalDim_nosym/(UniqueDim+NDim_images));
    printf("Expected savings factor = %.1f\n",
	   (double) TotalDim_nosym/(UniqueDim+NDim_images));
    printf("\n");

    if(Params::Parameters().PrintSymmetryInfo()){
      for (int i=1;i<=NDim;i++) {
	printf("The symmetry factor and periodic symmetry factor of d(%i,%i) is %i and  %i\n", 
	       Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),Dimers[i].GetSymmetryFactor(),
	       Dimers[i].GetPeriodicSymmetryFactor());
	//if(Dimers[i].GetSymmetryFactor()){
	//printf("Symmetrical to ");
	//for(int j=0;j<Dimers[i].GetSymmetryList().size()/2;j++)
	//	printf("d(%i,%i) ",
	//    Dimers[i].GetSymmetryList()[2*j],Dimers[i].GetSymmetryList()[2*j+1]);
	//printf("\n");
	//printf("and Symmetrical to Image Dimers ");
	//for(int j=0;j<Dimers[i].GetPeriodicSymmetryList().size()/2;j++)
	  //	printf("d(%i,%i) ",
		       //	       Dimers[i].GetPeriodicSymmetryList()[2*j],Dimers[i].GetPeriodicSymmetryList()[2*j+1]);
	//printf("\n");
	//}
      }  
    }
    if(Params::Parameters().PrintSymmetryInfo()) {
      for (int i=1;i<=NDim_images;i++) {  
	printf("The symmetry factor of d(%i,%i) is %i\n", 
	       DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB(),DimerImages[i].GetSymmetryFactor());
	//printf("Symmetrical to ");
	//for(int j=0;j<DimerImages[i].GetSymmetryList().size()/2;j++)
	  // printf("d(%i,%i) ",
		 //	     DimerImages[i].GetSymmetryList()[2*j],DimerImages[i].GetSymmetryList()[2*j+1]);
	//printf("\n");
      }
    }
  }
  /* Full Cluster */
  // geometry is read directly from monomers.  So just reset the energies.
  if(ResetEnergies) 
    ResetEnergiesAndGradients();

}

// Following should be used only for constraint volume optimization.
// In other words, check usability of this function before using it.
Vector Cluster::RearrangeCoordsForConstraintOpt(Vector coords) {

  Vector rearranged_coords(coords,false);

  for (int i = 0; i<6;i++) {       
    rearranged_coords[i] = coords[coords.GetLength()-6+i];
  }

  for (int i = 6; i<coords.GetLength();i++) {
    rearranged_coords[i] = coords[i-6];     
  }            
  return rearranged_coords;
}

// Following should be used only for constraint volume optimization.
// In other words, check usability of this function before using it.
Vector Cluster::RearrangeBackCoordsForConstraintOpt(Vector rearranged_coords) {

  Vector coords(rearranged_coords,false);

  for (int i = 0; i<(rearranged_coords.GetLength() - 6);i++) {
    coords[i] = rearranged_coords[i+6];
  }

  for (int i = (rearranged_coords.GetLength() - 6); i < rearranged_coords.GetLength();i++) {
    coords[i] = rearranged_coords[i - (rearranged_coords.GetLength() - 6) ];
  }

  return coords;
}

//Prints optimization step, due to using a reduced gradient, dlfind's opimization step output file does not output the entire crystal 
void Cluster::MakeOptimizationStepFile(){

  string filename = "Opt_Steps.xyz";


  FILE *Opt;
  int OptStep = Params::Parameters().GetOptCycle();
  //for the first step of the optimization, delete the Optimization step file
  if (OptStep==1){
    if ((Opt = fopen(filename.c_str(),"w"))==NULL){
      printf("Cluster::MakeOptimizationStepFile : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
  }
  //for every other step, append the file
  else
    if ((Opt = fopen(filename.c_str(),"a"))==NULL){
      printf("Cluster::MakeOptimizationStepFile : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
 
  //first line of step prints the number of atoms in the unit cell 
  int Natoms = GetTotalNumberOfAtoms();
  //Include Lattice Parameters
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())
    Natoms += 2; 

  fprintf(Opt,"%i\n",Natoms);
  //second line is the optimization step
  fprintf(Opt,"Opt step %i\n",OptStep);

  //Printing 
  for (int i=1;i<=NMon;i++) {
    Monomers[i].PrintMonomerCartesian(Opt);
  }

  //print lattice parameters
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    fprintf(Opt,"%-2s  %10.6f  %10.6f  %10.6f\n", 
	    "XX",UnitCellAxes[0],UnitCellAxes[1],UnitCellAxes[2]);
    fprintf(Opt,"%-2s  %10.6f  %10.6f  %10.6f\n", 
	    "XX",UnitCellAngles[0]*DegreesToRadians,
	    UnitCellAngles[1]*DegreesToRadians,UnitCellAngles[2]*DegreesToRadians);
  }
  
  
  fclose(Opt);
}

void Cluster::CreateQMJobs() {

 if(Params::Parameters().UseFullQMOnly()){
  //printf("Creating full QM job only...\n");
  if(Params::Parameters().GetQMType()==5){
    CreateQuantumEspressoJob();
  }
  else if (Params::Parameters().GetQMType()==8){
    CreateDFTBJob();
  }
 }
 else {
  // Monomers
  printf("Creating QM Monomer jobs...  ");
  fflush(stdout); // flush the output stream


  //use monomer symmetry
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();
    
  for (int i=1;i<=NMon;i++) {
    if (Params::Parameters().TinkerDebug()) {
      Monomers[i].CreateMMJob(Monomers, NMon);
    }
    else{
      //if monomer symmetry is being exploited, 
      //only do job if symmetry factor is not zero
      int sym_fac = Monomers[i].GetSymmetryFactor();
      if(!UseMonSym ||sym_fac != 0)
	Monomers[i].CreateQMJob(Monomers,NMon);
	//Monomers[i].CreateQChemJob(Monomers, NMon);
    }
  }

  double c0 = Params::Parameters().GetLocalCutoff(0); 
  double c1 = Params::Parameters().GetLocalCutoff(1); 

  // Dimers
  printf("Creating QM Dimer jobs.\n");
  fflush(stdout); // flush the output stream
  for (int i=1;i<=NDim;i++) {
    if(Dimers[i].GetSymmetryFactor() != 0){//if symmetry factor is zero, do not create a qchem input file
      double separation =  Dimers[i].GetDimerSeparation();
      if( Params::Parameters().DoLocal2BodyTruncation()  && (separation > c0) )
	NDim_trunc++;
      else{
	if (Params::Parameters().TinkerDebug()) 
	  Dimers[i].CreateMMJob(Monomers, NMon);
	else
	  Dimers[i].CreateQMJob(Monomers,NMon);
	  //Dimers[i].CreateQChemJob(Monomers, NMon);
      }
    }
  }
  
  // If doing periodic boundary conditions, create additional dimer jobs
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      double separation =  DimerImages[i].GetDimerSeparation();
      if( Params::Parameters().DoLocal2BodyTruncation()  && (separation > c0) ) { // check the greater than sign
	NDimImages_trunc++;
      	printf("discarding MM D(%i,%i)\n",
	       DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB());
      }else{
	if (Params::Parameters().UseEmbeddingCharges() ) {
	  printf("Cluster::CreateQMJobs() ERROR: Embedding charges not yet implemented for periodic systems\n");
	  exit(1);
	}
      if (Params::Parameters().TinkerDebug()) 
	DimerImages[i].CreateMMJob(Monomers, NMon);
      else
	DimerImages[i].CreateQMJob(Monomers,NMon);
      }
    }
  }

  string filename = "Opt_Steps.xyz";


  FILE *Opt;
  int OptStep = Params::Parameters().GetOptCycle();
  //for the first step of the optimization, delete the Optimization step file
  if (OptStep==1){
    if ((Opt = fopen(filename.c_str(),"w"))==NULL){
      printf("Cluster::MakeOptimizationStepFile : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
  }
  //for every other step, append the file
  else
    if ((Opt = fopen(filename.c_str(),"a"))==NULL){
      printf("Cluster::MakeOptimizationStepFile : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
 
  //first line of step prints the number of atoms in the unit cell 
  int Natoms = GetTotalNumberOfAtoms();
  //Include Lattice Parameters
  if(Params::Parameters().IsPeriodic())
    Natoms += 2; 

  fprintf(Opt,"%i\n",Natoms);
  //second line is the optimization step
  fprintf(Opt,"Opt step %i\n",OptStep);

  //Printing 
  for (int i=1;i<=NMon;i++) {
    Monomers[i].PrintMonomerCartesian(Opt);
  }

  //print lattice parameters
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    fprintf(Opt,"%-2s  %10.6f  %10.6f  %10.6f\n", 
	    "XX",UnitCellAxes[0],UnitCellAxes[1],UnitCellAxes[2]);
    fprintf(Opt,"%-2s  %10.6f  %10.6f  %10.6f\n", 
	    "XX",UnitCellAngles[0]*DegreesToRadians,
	    UnitCellAngles[1]*DegreesToRadians,UnitCellAngles[2]*DegreesToRadians);
  }
 }
  //yoni: not using for now
  /*
  CreateFiniteDifferenceDimers(true);
  exit(0);
  */

}

void Cluster::CreateMMJobs() {

  //MM Type
  int MM_Type = Params::Parameters().GetMMType();

  //use AIFF symmetry for monomers
  bool UseMMSym = Params::Parameters().UseMMSymmetry();
  //use monomer symmetry
  //used by tinker and QChem
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();


  //if true only create job for the full MM
  if(Params::Parameters().UseFullMMOnly() && MM_Type != 2){
    printf("Creating full MM job only...");
  }else{

    // Monomers
    printf("Creating MM Monomer jobs...  ");
    fflush(stdout); // flush the output stream
    
    for (int i=1;i<=NMon;i++) {
      int sym_fac = Monomers[i].GetSymmetryFactor();
      if(MM_Type==2 && (sym_fac!=0 ||!UseMMSym) )//for AIFF, check if using MM symmetry
	Monomers[i].CreateMMJob(Monomers, NMon);
      else if(sym_fac!=0 || !UseMonSym){
	Monomers[i].CreateMMJob(Monomers, NMon);
      }
    }
    
    double c0 = Params::Parameters().GetLocalCutoff(0); 
    double c1 = Params::Parameters().GetLocalCutoff(1);

    // Dimers
    printf("Creating MM Dimer jobs.\n\n");
    fflush(stdout); // flush the output stream
    for (int i=1;i<=NDim;i++) {
      double separation =  Dimers[i].GetDimerSeparation();
      //for a series of symmetrical dimers, only run one
      if(Dimers[i].GetSymmetryFactor()!=0)
	if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0)
	  Dimers[i].CreateMMJob(Monomers, NMon);
    }
    
    // If doing periodic boundary conditions, create additional dimer jobs
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NDim_images;i++) {
	double separation =  DimerImages[i].GetDimerSeparation();
	if (Params::Parameters().UseEmbeddingCharges() ) {
	  printf("Cluster::CreateQMJobs() ERROR: Embedding charges not yet implemented for periodic systems\n");
	  exit(1);
	}
	//for a series of symmetrical dimers, only run one
	if(separation <= c0)
	  DimerImages[i].CreateMMJob(Monomers,NMon);
	else
	  printf("discarding MM D(%i,%i)\n",
		 DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB());
      }
    }
  }
  
  // Full Cluster
  if ( Params::Parameters().DoQMBenchmark() ) { // full cluster benchmark 
    CreateFullClusterQChemJob(false);
  }
  
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
    if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() ) {
      //if(Params::Parameters().IsPeriodic()){
      
      //if you don't have Dr. Beran's correction for the polarization, use this option
      if(Params::Parameters().Do_fdTinkerHessian()){
	CreateFiniteDifferenceTinkerJobs(3);
      }
      else{
	CreateFullClusterTinkerJob(3,false,"");
      }
    }
    else{
      CreateFullClusterTinkerJob(1,false,"");
    }
  }
  else if ( Params::Parameters().GetMMType()==2 ) { // AIFF
    // The AIFF energy for the full cluster is evaluated later, when
    // we read the MM energies.
  }
  else if ( Params::Parameters().GetMMType()==3 ) { // Q-Chem
    CreateFullClusterQChemJob(true);
  }
  else if (Params::Parameters().GetMMType()==5 ) { //Crystal09
    CreateCrystalJob();
    if(Params::Parameters().DoForces())
      CreateFiniteDifferenceCrystalJob();
  }
  else {
    printf("Cluster::CreateMMJobs: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  } 

}


void Cluster::CreateFiniteDifferenceDimers(bool CCSDT){

  double delta = 0.01;//in Angstoms
  double deltaLength = 0.01;//in Angstroms
  double deltaAngle = 0.1;//in degrees


  int i=0;
  int iCoord = 0;
  for(int iMon=1;iMon<=NMon;iMon++){
    for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
      if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	for(int xyz = 0;xyz<3;xyz++){
	  bool Skip = false;
	  if(Monomers[iMon].GetAtom(iAtom).IsAtomFrozen())
	    Skip = true;
	  if(Monomers[iMon].GetAtom(iAtom).IsYLocked() && xyz==1)
	    Skip = true;
	  if(Monomers[iMon].GetAtom(iAtom).IsZLocked() && xyz==2)
	    Skip = true;
	  
	  if(!Skip){

	    

	    //Plus Geometry
	    Vector Shift = SymmetryUniqueCoordinates;
	    Shift[iCoord] += delta;
	    MaintainCartesianSymmetry(Shift,false);
	    //Shift.PrintGradient("Unique Coord");
	    //looping over dimers and shifting dimers
	    for(int iDim=1;iDim<=NDim;iDim++){
	      if(Dimers[iDim].GetSymmetryFactor()){
		int MonA = Dimers[iDim].GetIndexA();
		int MonB = Dimers[iDim].GetIndexB();

		if(Monomers[MonA].GetSymmetricalMonomer() == iMon){
		  printf("d(%i,%i) MonA symmetrical to iMon)\n",
			 MonA,MonB);
		}
		if(Monomers[MonB].GetSymmetricalMonomer() == iMon){
		  printf("d(%i,%i) MonB symmetrical to iMon)\n",
			 MonA,MonB);
		}

	      }
	    }
	  }
	  iCoord++;
	}
	i++;
      }
    }
    
  }
  exit(0);
}

// -----------------------------------------------------------------
// by Ali begin
//
// GJB: only used for backup debugging as of 9/09.  We now use our own
// code for computing the AIFF energy.
void Cluster::CreateOrientJob() {
  
// Set up the Orient input filename, with the full path.  
   string datafile = Params::Parameters().GetMMPath() + "/orient.in";

/* Create data file */
   FILE *data;
   if ((data = fopen(datafile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateOrientJob : Cannot open file '%s'\n",
        datafile.c_str());
        exit(1);
    }

   double DampingFactor = Params::Parameters().GetDampingFactor();
   int Natoms = GetNumberOfAtoms(1); // number of atoms in the first monomer
   
   // write the data file
   fprintf(data,"\nUnits Bohr kJ/mol\n\n");

   fprintf(data,"Parameters\n");
   fprintf(data,"  molecules %d \n", NMon);
   fprintf(data,"  Sites %d Polarizable %d \n", NMon*(Natoms+1), NMon*Natoms);
   fprintf(data,"  S_functions %d \n", NMon*20000);
   fprintf(data,"End\n\n");

   //fprintf(data,"#include orient.defn\n\n");

   // start what was old defn file
   fprintf(data,"Variables\n");

   for (int imon=1;imon<NMon+1;imon++) {
     if (Params::Parameters().UseGlobalCoordinates() ) {
       fprintf(data,"x%d 0.0 A \n", imon);
       fprintf(data,"y%d 0.0 A \n", imon);
       fprintf(data,"z%d 0.0 A \n", imon);
     }
     else {
       fprintf(data,"x%d %10.6f A \n", imon, Monomers[imon].GetAtom(0).GetCoordinate(0));
       fprintf(data,"y%d %10.6f A \n", imon, Monomers[imon].GetAtom(0).GetCoordinate(1));
       fprintf(data,"z%d %10.6f A \n", imon, Monomers[imon].GetAtom(0).GetCoordinate(2));
     }
      fprintf(data,"alpha%d %10.6f D \n", imon, Monomers[imon].GetRotationAngle());
      fprintf(data,"a%d %10.6f A \n", imon, Monomers[imon].GetRotationVector(0));
      fprintf(data,"b%d %10.6f A \n", imon, Monomers[imon].GetRotationVector(1));
      fprintf(data,"c%d %10.6f A \n\n", imon, Monomers[imon].GetRotationVector(2));
   }

   fprintf(data,"End\n\n");
   fprintf(data,"Types\n");

   string AtmInt[Natoms];
   char label[10]; // index ought to have fewer than 10 digits

   for (int i=0;i<Natoms;i++){
      fprintf(data,"   %s%d   Z   %d   \n", Monomers[1].GetAtom(i).GetSymbol().c_str(), Monomers[1].GetAtom(i).GetAtomIndex(), Monomers[1].GetAtom(i).GetAtomicNumber());
   
   // We need the following when the moments file is read 
      AtmInt[i] = Monomers[1].GetAtom(i).GetSymbol().c_str();
      sprintf(label,"%d",Monomers[1].GetAtom(i).GetAtomIndex());
      AtmInt[i] += label;
   }
   
   fprintf(data,"End\n\n");
   //   fprintf(data,"Units Angstrom\n\n");
   
   for (int imon=1;imon<NMon+1;imon++){
     
     fprintf(data,"Molecule monomer%d at x%d y%d z%d rotated by alpha%d about a%d b%d c%d \n", imon, imon, imon, imon, imon, imon, imon, imon);
     
     // Set up the filename, with the full path.  
     char label[10];
     string filename = Params::Parameters().GetMMPath() + "/m";
     sprintf(label,"%d",imon);
     filename += label;
     filename += ".mom";   
     
     ifstream infile;
     infile.open( filename.c_str() );
     if ( !infile.is_open() ) {
       printf("Cluster::CreateOrientJob: Cannot open file '%s'\n",
	      filename.c_str());
       exit(1);
     }
     
     // Read in moments file,  write to data file
     string line;
     int i=0;
     string entryAtm;
     while ( !infile.eof() ) {
       getline(infile,line);
       i++;
       if (i>2){
	 
	 istringstream iss(line);
	 string entry;
	 iss >> entry;
	 bool entryBool = false;
	 for(int n=0;n<Natoms;n++){
	   if (entry==AtmInt[n]){
	     entryBool = true;
	     entryAtm = entry;}
	 }
	 if (entryBool){
	   fprintf(data,"%s at ", entry.c_str());
	   iss >> entry;
	   fprintf(data,"%s ", entry.c_str());
	   iss >> entry;
	   fprintf(data,"%s ", entry.c_str());
	   iss >> entry;
	   fprintf(data,"%s ", entry.c_str());
	   fprintf(data,"Type %s ", entryAtm.c_str());
	   iss >> entry;
	   fprintf(data,"%s ", entry.c_str());
	   iss >> entry;
	   fprintf(data,"%s \n", entry.c_str());
	   getline(infile,line);
	 }
	 
         fprintf(data,"%s \n", line.c_str());
       }
     }
     
     fprintf(data,"End \n\n");
     
     infile.close();
     
   } // end for                                              
   
   for (int imon=1;imon<NMon+1;imon++){
     
     // Set up the filename, with the full path.  
     char label[10];
     string filename = Params::Parameters().GetMMPath() + "/m";
     sprintf(label,"%d",imon);
     filename += label;
     filename += ".pol";
     
     ifstream infile;
     infile.open( filename.c_str() );
     if ( !infile.is_open() ) {
       printf("Cluster::CreateOrientJob: Cannot open file '%s'\n",
	      filename.c_str());
       exit(1);
     }
     
     fprintf(data,"Polarizabilities for monomer%d \n\n", imon);
     
     // Read in polarizabilities file, write to data file
     string line;
     int i=0;
     string entryAtm;
     while ( !infile.eof() ) {
       getline(infile,line);
       i++;
       if (i>2){
	 istringstream iss(line);
	 string entry;
	 iss >> entry;
	 bool entryBool = false;
	 for (int n=0;n<Natoms;n++){
	   if (entry==AtmInt[n]){
	     entryBool = true;
	     entryAtm = entry;}
	 }
	 if (entryBool){
	   fprintf(data,"Read rank 2 site %s\n", entryAtm.c_str());
	   getline(infile,line);}
	 
	 fprintf(data,"%s \n", line.c_str());}
     }
     
     infile.close();
     
     fprintf(data,"End \n\n");
     
   } // end for                                              
   
   // end old data file
   
   
   
   fprintf(data,"plot xyz new file geom.xyz title \"Geometry\"\n\n");
   
   fprintf(data,"Switch iterate on\n");
   fprintf(data,"Pairs\n");
   fprintf(data,"   induction damping factor %f \n", DampingFactor); 
   fprintf(data,"End\n\n");
   
   // write the part to control energy calculations
   fprintf(data,"Note *** Total Cluster Energy\n");
   fprintf(data,"Energy\n\n");
   fprintf(data,"Note *** Compute Dimer energies\n\n");

   int store[NMon+1];
   for (int x=1; x<NMon+1;x++){
	   store[x]=1;}
   int a;
   for (int x=1; x<NMon+1;x++){
     if (store[x]==1)
       a = 1;
     else{
       fprintf(data,"Restore Molecule monomer%d \n", x);
       store[x]=1;}
     for (int y=x+1; y<NMon+1;y++){
       fprintf(data,"Note (%d,%d) Interaction \n", x, y);
       for (int z=1; z<NMon+1;z++){
	 if (store[z]==1){
	   if (z!=x){
	     fprintf(data,"Suppress Molecule monomer%d \n", z);
	     store[z]=0;}
	   else
	     a=1;	  
	 }
	 else
	   a=1;     
       }
       if (store[y]==1)
	 a=1;
       else{
	 fprintf(data,"Restore Molecule monomer%d \n", y);
	 store[y]=1;}
       fprintf(data,"Energy \n\n");
     }
   }
   // end energy part
   
   fprintf(data,"time\n\n");
   fprintf(data,"finish\n\n");
   
   
   fclose(data);
   
}
// end
// -----------------------------------------------------------------

void Cluster::CreateQuantumEspressoJob() {

  // Set up the filename, with the full path.  File is e.g. 'fullQE.in'
  string filename,filename2,filename3,filename4; 
  string path = Params::Parameters().GetQMPath();
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();
 
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable())
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  filename = path + "/fullQE.in";
  if(!Params::Parameters().UsePhonopyFreq()){
    filename2 = path + "/run.sh";
    filename3 = path + "/qe.ph.in";
  }
  else{
    filename2 = path + "/header";
    filename3 = path + "/run.sh";
  }
  filename4 = path + "/qpoints";
  //printf("%s\n",path.c_str());

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::CreateQuantumEspressoJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }  

  int numAtomTypes = GetTypesOfAtoms();
  //printf("Monomers[1].GetNumberOfAtomTypes()= %d\n",Monomers[1].GetNumberOfAtomTypes());
  
 // fprintf(job,"&Control\n%s\n/",Params::Parameters().GetQuantumEspressoControlSection().c_str());
  fprintf(job,"&CONTROL\n");
  fprintf(job,"  calculation = 'scf',\n");
  fprintf(job,"  restart_mode = 'from_scratch',\n");
  fprintf(job,"  prefix = \'fullQE%i\',\n",Params::Parameters().GetPID());
  fprintf(job,"  disk_io = \'high\',\n");
  fprintf(job,"  verbosity = \'high\',\n");
  fprintf(job,"  etot_conv_thr = 2.0D-6,\n");
  fprintf(job,"  forc_conv_thr = 6.0D-4,\n");
  if(Params::Parameters().DoForces()){
    fprintf(job,"  tprnfor = .TRUE.,\n");
    fprintf(job,"  tstress = .TRUE.,\n");
  }
  fprintf(job,"  wf_collect = .TRUE.,\n");
  fprintf(job,"/\n&SYSTEM\n");
  fprintf(job,"  ibrav = 0,\n");
  fprintf(job,"  nat = %i,\n",GetTotalNumberOfAtoms());
  fprintf(job,"  ntyp = %i,\n",numAtomTypes);
  fprintf(job,"  ecutwfc = %f,\n",Params::Parameters().GetQEBasis());
  fprintf(job,"  ecutrho = %f,\n",Params::Parameters().GetQEBasis()*Params::Parameters().GetQEBasisMultiplier());
  fprintf(job,"  nosym = .TRUE.,\n");
  if(Params::Parameters().ForceXCSet()){
    printf("WARNING!!! Enforcing exchange functional provided in input\n");
    fprintf(job,"  input_dft = '%s',\n",Params::Parameters().ForceXCType().c_str());
  }
  if(Params::Parameters().DispersionCorrection()){
    if(Params::Parameters().DispersionType()=="D2"){
      printf("Using Grimme's D2 dispersion correction\n");
      fprintf(job,"  vdw_corr = 'DFT-D',\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-PW86PBE"){
      printf("Using PW86PBE XDM dispersion correction\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.6836,\n");
      fprintf(job,"  xdm_a2 = 1.5045,\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-B86BPBE"){
      printf("Using B86BPBE XDM dispersion correction\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.6512,\n");
      fprintf(job,"  xdm_a2 = 1.4633,\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-BLYP"){
      printf("Using BLYP XDM dispersion correction\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.4502,\n");
      fprintf(job,"  xdm_a2 = 1.6210,\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-REVPBE"){
      printf("WARNING!!! Using REVPBE XDM dispersion correction. Might need reparameterization\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.3454,\n");
      fprintf(job,"  xdm_a2 = 1.9225,\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-PBESOL"){
      printf("WARNING!!! Using PBESOL XDM dispersion correction. Might need reparameterization\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.0000,\n");
      fprintf(job,"  xdm_a2 = 4.1503,\n/");
    }
    else if(Params::Parameters().DispersionType()=="XDM-PBE" || Params::Parameters().DispersionType()=="XDM"){
      printf("Using PBE XDM dispersion correction\n");
      fprintf(job,"  vdw_corr = 'XDM',\n");
      fprintf(job,"  xdm_a1 = 0.3275,\n");
      fprintf(job,"  xdm_a2 = 2.7673,\n/");
    }
    else{
      printf("WARNING!!! Unknown dispersion correction. Exiting HMBI\n");
      exit(0);
    }
    
  }
  else
    fprintf(job,"/");
  fprintf(job,"\n&ELECTRONS\n");
  fprintf(job,"  electron_maxstep = 1500,\n");
  fprintf(job,"  conv_thr = 1.D-8,\n");
  fprintf(job,"  mixing_beta = 0.5D0,\n/");
  //fprintf(job,"\n&IONS\n");
  //fprintf(job,"  ion_dynamics = \"bfgs\",\n/");
  //fprintf(job,"\n&CELL\n/ \n");
  fprintf(job,"\nATOMIC_SPECIES\n%s\n",Params::Parameters().GetQuantumEspressoSpeciesSection().c_str());
  fprintf(job,"ATOMIC_POSITIONS crystal\n");
  PrintFractionalCoordinates(job,false);
  //fprintf(job,"ATOMIC_POSITIONS {angstrom}\n");
  //PrintCartesianCoordinates(job);
  fprintf(job,"\nK_POINTS automatic\n");
  fprintf(job,"%s\n",Params::Parameters().GetKPoints().c_str());
  //fprintf(job,"1 1 1 1 1 1\n\n");
  fprintf(job,"CELL_PARAMETERS angstrom\n");
  Matrix unit_cell(3,3);
  unit_cell.SetRowVector(Cluster::cluster().GetUnitCellVector(0),0);
  unit_cell.SetRowVector(Cluster::cluster().GetUnitCellVector(1),1);
  unit_cell.SetRowVector(Cluster::cluster().GetUnitCellVector(2),2);
  //unit_cell.Print("unit_cell");
  for (int i=0;i<3;i++){
    fprintf(job,"%.12f   %.12f   %.12f\n",unit_cell(i,0),unit_cell(i,1),unit_cell(i,2));
  }

  fclose(job);

  if (Params::Parameters().DoFreq()) {

    int sizeA = 1;
    int sizeB = 1;
    int sizeC = 1;
    int supercell_Natoms;
    string qpoints = "";

    FILE *job4;
    if ((job4 = fopen(filename4.c_str(),"w"))==NULL) {
      printf("Cluster::CreateQuantumEspressoJob : Cannot open file '%s'\n",filename4.c_str());
      exit(1);
    }

    if(Params::Parameters().IsSupercellJob()){
      Vector Supercell_Size = Params::Parameters().GetSupercellSize();
      sizeA = Supercell_Size[0];
      sizeB = Supercell_Size[1];
      sizeC = Supercell_Size[2];
      int Ncells = int(sizeA * sizeB * sizeC);
      supercell_Natoms = Ncells * GetTotalNumberOfAtoms();

      double tmp;
      string tmpA, tmpB, tmpC;
      
      //Now make the uniformly spaced grid according to Monkhorst-Pack scheme (Phys. Rev. B, Volume 13, No. 12, Page 5188, 1976)
      //As of 31st October 2017, we do not consider any symmetry while creating this simple grid
      //This code was borrowed from Supercell.C (CreateReciprocalSpaceSamplingGrid)
      for (int i=1;i<=sizeA;i++) {
        for (int j=1;j<=sizeB;j++) {
          for (int k=1;k<=sizeC;k++) {
            tmp = (2*i-sizeA-1.0)/2/sizeA;
            stringstream stream;
            stream << tmp;
            tmpA = stream.str();
            stream.str("");
            //stream.clear();
            tmp = (2*j-sizeB-1.0)/2/sizeB;
            stream << tmp;
            tmpB = stream.str();
            stream.str("");
            //stream.clear();
            tmp = (2*k-sizeC-1.0)/2/sizeC;
            stream << tmp;
            tmpC = stream.str();
            stream.str("");
            //stream.clear();
            qpoints +=  tmpA + " " + tmpB + " " + tmpC + "  ";
            fprintf(job4,"%i%i%i\t%s\t%s\t%s\n",i-1,j-1,k-1,tmpA.c_str(),tmpB.c_str(),tmpC.c_str());
          }
        }
      }
    }
    else {
     qpoints = "0.00 0.00 0.00";
     fprintf(job4,"000\t%s\n",qpoints.c_str());
     supercell_Natoms = GetTotalNumberOfAtoms();
    }

    fclose(job4);

    // Open the input file for writing
    FILE *job2;
    if ((job2 = fopen(filename2.c_str(),"w"))==NULL) {
      printf("Cluster::CreateQuantumEspressoJob : Cannot open file '%s'\n",filename2.c_str());
      exit(1);
    } 

    if(!Params::Parameters().UsePhonopyFreq()) {
      fprintf(job2,"#!/bin/bash\n\n");
      fprintf(job2,"numProc=$1\n\n");
      fprintf(job2,"#First need to run the scf calculation\n");
      fprintf(job2,"sed -i \"/disk_io/d\" fullQE.in\n");
      fprintf(job2,"mpirun -np ${numProc} pw.x -i fullQE.in > fullQE.out\n\n");
      fprintf(job2,"front=\"qe.ph\"\n\n");
      fprintf(job2,"list=`cat qpoints |wc -l`\n");
      fprintf(job2,"for i in `seq 1 1 $list`\n");
      fprintf(job2,"do\n");
      fprintf(job2,"  line=`awk \"NR==$i\" qpoints`\n");
      fprintf(job2,"  back=`echo $line | awk '{print $1}'`\n");
      fprintf(job2,"  qpoints=`echo $line | awk '{print $2,$3,$4}'`\n");
      fprintf(job2,"  name=${front}${back}.in\n");
      fprintf(job2,"  out=${front}${back}.out\n");
      fprintf(job2,"  mkdir -p $back\n");
      fprintf(job2,"  cp qe.ph.in ${back}/${name}\n");
      fprintf(job2,"  cp -r fullQE* ${back}\n");
      fprintf(job2,"  cd $back\n");
      fprintf(job2,"  sed -i \"s/qpoints/${qpoints}/\" $name\n");
      fprintf(job2,"  mpirun -np ${numProc} ph.x -i $name > $out\n");
      fprintf(job2,"  cp dyn.G ../dyn${back}.G\n");
      fprintf(job2,"  cd ../\n");
      fprintf(job2,"done\n");
    }
    else{
      fprintf(job2,"&CONTROL\n");
      fprintf(job2,"  calculation = 'scf',\n");
      fprintf(job2,"  etot_conv_thr = 2.0D-6,\n");
      fprintf(job2,"  forc_conv_thr = 6.0D-4,\n");
      fprintf(job2,"  tprnfor = .TRUE.,\n");
      fprintf(job2,"  tstress = .TRUE.,\n");
      fprintf(job2,"  prefix = \'fullQEthis%i\',\n",Params::Parameters().GetPID());
      fprintf(job2,"  disk_io = \'high\',\n");
      fprintf(job2,"  wf_collect = .TRUE.,\n");
      fprintf(job2,"/\n&SYSTEM\n");
      fprintf(job2,"  ibrav = 0,\n");
      fprintf(job2,"  nat = %i,\n",supercell_Natoms);
      fprintf(job2,"  ntyp = %i,\n",numAtomTypes);
      fprintf(job2,"  ecutwfc = %f,\n",Params::Parameters().GetQEBasis());
      fprintf(job2,"  ecutrho = %f,\n",Params::Parameters().GetQEBasis()*Params::Parameters().GetQEBasisMultiplier());
      fprintf(job2,"  nosym = .TRUE.,\n");
      if(Params::Parameters().ForceXCSet()){
        printf("WARNING!!! Enforcing exchange functional provided in input\n");
        fprintf(job2,"  input_dft = '%s',\n",Params::Parameters().ForceXCType().c_str());
      }
      if(Params::Parameters().DispersionCorrection()){
        if(Params::Parameters().DispersionType()=="D2"){
          printf("Using Grimme's D2 dispersion correction\n");
          fprintf(job2,"  vdw_corr = 'DFT-D',\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-PW86PBE"){
          printf("Using PW86PBE XDM dispersion correction\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.6836,\n");
          fprintf(job2,"  xdm_a2 = 1.5045,\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-B86BPBE"){
          printf("Using B86BPBE XDM dispersion correction\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.6512,\n");
          fprintf(job2,"  xdm_a2 = 1.4633,\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-BLYP"){
          printf("Using BLYP XDM dispersion correction\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.4502,\n");
          fprintf(job2,"  xdm_a2 = 1.6210,\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-REVPBE"){
          printf("WARNING!!! Using REVPBE XDM dispersion correction. Might need reparameterization\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.3454,\n");
          fprintf(job2,"  xdm_a2 = 1.9225,\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-PBESOL"){
          printf("WARNING!!! Using PBESOL XDM dispersion correction. Might need reparameterization\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.0000,\n");
          fprintf(job2,"  xdm_a2 = 4.1503,\n/");
        }
        else if(Params::Parameters().DispersionType()=="XDM-PBE" || Params::Parameters().DispersionType()=="XDM"){
          printf("Using PBE XDM dispersion correction\n");
          fprintf(job2,"  vdw_corr = 'XDM',\n");
          fprintf(job2,"  xdm_a1 = 0.3275,\n");
          fprintf(job2,"  xdm_a2 = 2.7673,\n/");
        }
        else{
          printf("WARNING!!! Unknown dispersion correction. Exiting HMBI\n");
          exit(0);
        }
    
      }
      else
        fprintf(job2,"/");
      fprintf(job2,"\n&ELECTRONS\n");
      fprintf(job2,"  conv_thr = 1.D-8,\n");
      fprintf(job2,"  electron_maxstep = 1500,\n");
      fprintf(job2,"/\nK_POINTS automatic\n");
      if(sizeA > 1 || sizeB > 1 || sizeC > 1){
        printf("Doing Supercell so dropping KPOINTS to gamma point only\n");
        fprintf(job2,"1 1 1 1 1 1\n");
      }
      else
        fprintf(job2,"%s\n",Params::Parameters().GetKPoints().c_str());
    }

    fclose(job2);

    FILE *job3;
    if ((job3 = fopen(filename3.c_str(),"w"))==NULL) {
      printf("Cluster::CreateQuantumEspressoJob : Cannot open file '%s'\n",filename3.c_str());
      exit(1);
    } 

    if(!Params::Parameters().UsePhonopyFreq()) {
      //When using QE's ph.x to create frequencies
      //**currently cannot be done with dispersion correction
      fprintf(job3,"freq calc\n");
      fprintf(job3,"&INPUTPH\n");
      fprintf(job3,"  prefix = \'fullQE\',\n");
      fprintf(job3,"  verbosity = \'high\',\n");
      fprintf(job3,"  fildyn = \'dyn.G\',\n");
      fprintf(job3,"  tr2_ph = 1.0d-14,\n");
      fprintf(job3,"  outdir = \"./\",\n/");
      fprintf(job3,"\nqpoints\n");
    }
    else {
      fprintf(job3,"#!/bin/bash\n\n");
      fprintf(job3,"name=$1\n");
      fprintf(job3,"numproc=$2\n");
      fprintf(job3,"dim=$3\n\n");
      fprintf(job3,"#create displacements\n");
      fprintf(job3,"phonopy --pwscf -d --dim=\"$dim\" -c ${name}.in\n\n");
      fprintf(job3,"#make a list of all the jobs to run\n");
      fprintf(job3,"ls supercell-* | awk '{print $1}' | cut -d'.' -f 1 | sed \"s/^[^-]*-//\" > list\n");
      fprintf(job3,"num_atoms=`wc -l list | awk '{print $1}'`\n\n");
      fprintf(job3,"#modifies the number of atoms to account for the supercell\n");
      fprintf(job3,"sed -i \"s/tmp/`grep nat supercell.in | awk '{print $7}' | cut -d, -f1`/\" header\n");
      fprintf(job3,"#loop to run jobs\n");
      fprintf(job3,"for i in `seq 1 $num_atoms`;\n");
      fprintf(job3,"do\n\n");
      fprintf(job3,"#create jobs\n");
      fprintf(job3,"j=`awk \"NR==$i\" list`\n");
      fprintf(job3,"cat header supercell-$j.in > $name-$j.in\n");
      fprintf(job3,"mkdir -p $j\n");
      fprintf(job3,"mv $name-$j.in $j\n");
      fprintf(job3,"cd $j\n\n");
      fprintf(job3,"#run job\n");
      fprintf(job3,"sed -i \"s/this/${j}/g\" $name-$j.in\n");
      fprintf(job3,"mpirun -np $numproc pw.x -i $name-$j.in > $name-$j.out\n");
      fprintf(job3,"rm -rf pwscf.*\n");
      fprintf(job3,"cd ../\n\n");
      fprintf(job3,"done\n\n");
      fprintf(job3,"#analyze results\n");
      fprintf(job3,"phonopy --pwscf -f */*.out #-c $name.in\n");
      fprintf(job3,"#write dynamical matrix\n");
      fprintf(job3,"phonopy --pwscf --dim=\"$dim\" -c $name.in --qpoints=\"0.0 0.0 0.0\" --writedm\n");
      fprintf(job3,"mv qpoints.yaml qpoints_gamma.yaml\n");
      if(Params::Parameters().IsSupercellJob())
        fprintf(job3,"phonopy --pwscf --dim=\"$dim\" -c $name.in --qpoints=\"%s\" --writedm --writefc --fc-spg-symmetry --fc-symmetry\n",qpoints.c_str());
      fprintf(job3,"sed -i \"s/\\,//g\" qpoints*.yaml\n");
      fprintf(job3,"rm list\n");
    }

    fclose(job3);
  }
 
}

void Cluster::CreateDFTBJob() {

  // Set up the filename, with the full path.  File is e.g. 'dftb_in.hsd'
  string filename,filename2,filename3,filename4; 
  string path = Params::Parameters().GetQMPath();
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();
 
  printf("Correctly found the create DFTB job section\n");

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable())
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  filename = path + "/dftb_gamma_in.hsd";
  filename2 = path + "/modes_gamma_in.hsd";
  filename3 = path + "/dftb_phn_in.hsd";
  filename4 = path + "/modes_phn_in.hsd";
  //printf("%s\n",path.c_str());
  //fflush(stdout);

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::CreateDFTBJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }  

  fprintf(job,"Geometry {\n");
  //Maybe do an if statement for supercell?
  fprintf(job," %s\n",GetDFTBGeomInput().c_str());
  //fprintf(job," }\n");
  fprintf(job,"  Periodic = Yes\n");
  fprintf(job,"  LatticeVectors [Angstrom] = {\n");
  //Vector tmp_vec;
  for (int i=0; i<3; i++){
    Vector tmp_vec;
    tmp_vec=GetUnitCellVector(i);
    fprintf(job,"  %f\t%f\t%f\n",tmp_vec[0],tmp_vec[1],tmp_vec[2]);
  }
  fprintf(job,"  }\n");
  fprintf(job,"}\n");

  if(Params::Parameters().DoFreq()){
    //For Freq; Also set SCCTolerance = 1e-7 or better
    fprintf(job,"Driver = SecondDerivatives {\n");
    fprintf(job,"   Atoms = 1:-1\n");
    fprintf(job,"   Delta = 1e-4\n");
    fprintf(job,"}\n");
  }
  else
    fprintf(job,"Driver = {}\n");
    /*fprintf(job,"Driver = SteepestDescent {\n");
    fprintf(job,"   MovedAtoms = 1:-1\n");
    fprintf(job,"   MaxSteps = 0\n");
    fprintf(job,"   OutputPrefix = \"geo_end\"\n");
    //fprintf(job,"   ConvergentForcesOnly = Yes\n");
    //Add a tag here to c
    if(Params::Parameters().FreezeUnitCellParams())
      fprintf(job,"   LatticeOpt = No\n");
    fprintf(job,"   MaxForceComponent = 1e-6\n");
    fprintf(job,"}\n");*/
  fprintf(job,"Hamiltonian = DFTB {\n");
  fprintf(job,"  Filling = Fermi { \n");
  fprintf(job,"    Temperature [K] = 0  \n");
  fprintf(job,"  } \n");
  fprintf(job,"  Charge = 0 \n");
  fprintf(job,"  SCC = Yes \n");
  fprintf(job,"  SCCTolerance = 1.0E-08\n");
  fprintf(job,"  MaxSCCIterations = 300\n");
  fprintf(job,"  EwaldParameter = 0.0\n");
  fprintf(job,"  EwaldTolerance = 1e-9\n");
  //fprintf(job,"  OrbitalResolvedSCC = yes\n"); <--Check this parameter! JLM
  fprintf(job,"  DampXH = Yes\n");
  fprintf(job,"  DampXHExponent = 0.42000E+01\n");
  fprintf(job,"  ThirdOrderFull = Yes\n");
  fprintf(job,"  ReadInitialCharges = no \n");
  fprintf(job,"  Mixer = Broyden{} \n");
  fprintf(job,"  SpinPolarisation = {} \n");
  //fprintf(job," ForceEvaluation = \'traditional\'\n");

  //Need to add Hubbard Derivatives here
  fprintf(job,"%s\n",Params::Parameters().GetHubbardDerivatives().c_str());

  //Control if Disp. is on/off and if 3B Dispersion is on/off
  if(Params::Parameters().DispersionCorrection()){
    fprintf(job," Dispersion = DftD3 {\n");
    fprintf(job,"   Damping = BeckeJohnson {}\n");
    //JLM this is for adjusted BJ damping
    /*fprintf(job,"    Damping = BeckeJohnson {\n");
    fprintf(job,"      a1 = 0.5719\n");
    fprintf(job,"      a2 = 3.6017\n");
    fprintf(job,"    }\n");
    fprintf(job,"    s6 = 1.0\n");
    fprintf(job,"    s8 = 0.5883\n");*/ 
    if(Params::Parameters().EstimateThreeBodyDispersion())
      fprintf(job,"   threebody = yes\n");
    fprintf(job," }\n");
  }

  //Angular momentum data here
  fprintf(job,"%s\n",Params::Parameters().GetMaxAngularMomentum().c_str());

  //CHECK THIS VERY CAREFULLY -- Is the basis set?
  //Make this a read-in parameter
  fprintf(job,"%s\n",Params::Parameters().GetSlaterKoster().c_str());

  fprintf(job," KPointsAndWeights = { \n");
  //Get Kpoints and set = size{A,B,C}
  string kpoints=Params::Parameters().GetKPoints();
  //printf("kpoints = %s\n",kpoints.c_str());
  
  int sizeA, sizeB, sizeC, supercell_sizeA, supercell_sizeB, supercell_sizeC;

  istringstream iss(kpoints);
  iss >> sizeA;
  iss >> sizeB;
  iss >> sizeC;
  //Maybe add weight{A/B/C} as well?
  /*printf("sizeA = %i\tsizeB = %i\tsizeC = %i\n",sizeA,sizeB,sizeC);
  fflush(stdout);
  exit(0);*/

  double tmpA, tmpB, tmpC;
  //double tmp;
  //Add Kpoint grid here
  if(sizeA>1 || sizeB>1 ||sizeC > 1){
    //Now make the uniformly spaced grid according to Monkhorst-Pack scheme (Phys. Rev. B, Volume 13, No. 12, Page 5188, 1976)
    //As of 31st October 2017, we do not consider any symmetry while creating this simple grid
    //This code was borrowed from Supercell.C (CreateReciprocalSpaceSamplingGrid)
    for (int i=1;i<=sizeA;i++) {
      for (int j=1;j<=sizeB;j++) {
        for (int k=1;k<=sizeC;k++) {
          tmpA = (2*i-sizeA-1.0)/2/sizeA;
          tmpB = (2*j-sizeB-1.0)/2/sizeB;
          tmpC = (2*k-sizeC-1.0)/2/sizeC;
          fprintf(job,"   %f\t%f\t%f\t1.0\n",tmpA,tmpB,tmpC);
        }
      }
    }
  }
  else {
    fprintf(job," 0.00\t0.00\t0.00\t1.0\n");
  }
  fprintf(job," } \n");

  fprintf(job,"}\n");
  if(Params::Parameters().DoForces()){
    fprintf(job,"Analysis {\n");
    fprintf(job," CalculateForces = yes\n");
    fprintf(job,"}\n");
  }
  fprintf(job," \n");

  fclose(job);

  //Note: We don't actually need this section. The modes is used purely to ensure
  //Our answer is the same as DFTB+
  if (Params::Parameters().DoFreq()) {
    FILE *job2;
    if ((job2 = fopen(filename2.c_str(),"w"))==NULL) {
      printf("Cluster::CreateDFTBJob : Cannot open file '%s'\n",filename2.c_str());
      exit(1);
    }

    fprintf(job2,"Geometry {\n");
    //fprintf(job,"  %i\n",GetTotalNumberOfAtoms()); //<-- Don't need this line
    //Add in Molecular Geometry here...
    fprintf(job2," %s\n",GetDFTBGeomInput().c_str());
    //fprintf(job," }\n");
    fprintf(job2,"  Periodic = Yes\n");
    fprintf(job2,"  LatticeVectors [Angstrom] = {\n");
    //Vector tmp_vec;
     for (int i=0; i<3; i++){
      Vector tmp_vec=GetUnitCellVector(i);
      fprintf(job2,"  %f\t%f\t%f\n",tmp_vec[0],tmp_vec[1],tmp_vec[2]);
    }
    fprintf(job2,"  }\n");
    fprintf(job2,"}\n");
    fprintf(job2,"Hessian = {\n");
    fprintf(job2,"  <<< \"hessian.out\"\n");
    fprintf(job2,"}\n");
    fprintf(job2,"%s\n",Params::Parameters().GetSlaterKoster().c_str());

    fclose(job2);

    if(Params::Parameters().IsSupercellJob()) {
      /*if(!Supercell::supercell().GetSupercellInit()){
        printf("Attempting to initialize Supercell\n");
        fflush(stdout);
        
        Supercell::supercell().Initialize(Main::Main().GetInfileName().c_str(),Params::Parameters().GetNumberOfProcessors(),Main::Main().GetFilename().c_str());
        
      }*/

      FILE *job3;
      if ((job3 = fopen(filename3.c_str(),"w"))==NULL) {
        printf("Cluster::CreateDFTBJob : Cannot open file '%s'\n",filename3.c_str());
        exit(1);
      }  

      fprintf(job3,"Geometry {\n");
      fprintf(job3," %s\n",Supercell::supercell().GetDFTBGeomInput().c_str());
      fprintf(job3,"  Periodic = Yes\n");
      fprintf(job3,"  LatticeVectors [Angstrom] = {\n");
      for (int i=0; i<3; i++){
        Vector tmp_vec;
        tmp_vec=Supercell::supercell().GetUnitCellVector(i);
        fprintf(job3,"  %f\t%f\t%f\n",tmp_vec[0],tmp_vec[1],tmp_vec[2]);
      }
      fprintf(job3,"  }\n");
      fprintf(job3,"}\n");

      //For Freq; Also set SCCTolerance = 1e-7 or better
      fprintf(job3,"Driver = SecondDerivatives {\n");
      fprintf(job3,"   Atoms = 1:-1\n");
      fprintf(job3,"   Delta = 1e-4\n");
      fprintf(job3,"}\n");
      fprintf(job3,"Hamiltonian = DFTB {\n");
      fprintf(job3,"  Filling = Fermi { \n");
      fprintf(job3,"    Temperature [K] = 0  \n");
      fprintf(job3,"  } \n");
      fprintf(job3,"  Charge = 0 \n");
      fprintf(job3,"  SCC = Yes \n");
      fprintf(job3,"  SCCTolerance = 1.0E-08\n");
      fprintf(job3,"  MaxSCCIterations = 300\n");
      fprintf(job3,"  EwaldParameter = 0.0\n");
      fprintf(job3,"  EwaldTolerance = 1e-9\n");
      //fprintf(job3,"  OrbitalResolvedSCC = yes\n"); <--Check this parameter! JLM
      fprintf(job3,"  DampXH = Yes\n");
      fprintf(job3,"  DampXHExponent = 0.42000E+01\n");
      fprintf(job3,"  ThirdOrderFull = Yes\n");
      fprintf(job3,"  ReadInitialCharges = no \n");
      fprintf(job3,"  Mixer = Broyden{} \n");
      fprintf(job3,"  SpinPolarisation = {} \n");
      //fprintf(job3," ForceEvaluation = \'traditional\'\n");

      //Need to add Hubbard Derivatives here
      fprintf(job3,"%s\n",Params::Parameters().GetHubbardDerivatives().c_str());

      //Control if Disp. is on/off and if 3B Dispersion is on/off
      if(Params::Parameters().DispersionCorrection()){
        fprintf(job3," Dispersion = DftD3 {\n");
        fprintf(job3,"   Damping = BeckeJohnson {}\n");
        if(Params::Parameters().EstimateThreeBodyDispersion())
          fprintf(job3,"   threebody = yes\n");
        fprintf(job3," }\n");
      }

      //Angular momentum data here
      fprintf(job3,"%s\n",Params::Parameters().GetMaxAngularMomentum().c_str());

      //CHECK THIS VERY CAREFULLY -- Is the basis set?
      //Make this a read-in parameter
      fprintf(job3,"%s\n",Params::Parameters().GetSlaterKoster().c_str());

      fprintf(job3," KPointsAndWeights = { \n");
      //Get Kpoints and set = size{A,B,C}
      string kpoints=Params::Parameters().GetKPoints();
      //printf("kpoints = %s\n",kpoints.c_str());
  
      int sizeA, sizeB, sizeC, supercell_sizeA, supercell_sizeB, supercell_sizeC;

      istringstream iss(kpoints);
      iss >> sizeA;
      iss >> sizeB;
      iss >> sizeC;

      if(Params::Parameters().IsSupercellJob() && Params::Parameters().DoFreq()){
        Vector Supercell_Size = Params::Parameters().GetSupercellSize();
        supercell_sizeA = Supercell_Size[0];
        supercell_sizeB = Supercell_Size[1];
        supercell_sizeC = Supercell_Size[2];
        if(supercell_sizeA > 1)
          sizeA = 1;
        if(supercell_sizeB > 1)
          sizeB = 1;
        if(supercell_sizeC > 1)
          sizeC = 1;
        }

      double tmpA, tmpB, tmpC;
      //double tmp;
      //Add Kpoint grid here
      if(sizeA>1 || sizeB>1 ||sizeC > 1){
        //Now make the uniformly spaced grid according to Monkhorst-Pack scheme (Phys. Rev. B, Volume 13, No. 12, Page 5188, 1976)
        //As of 31st October 2017, we do not consider any symmetry while creating this simple grid
        //This code was borrowed from Supercell.C (CreateReciprocalSpaceSamplingGrid)
        for (int i=1;i<=sizeA;i++) {
          for (int j=1;j<=sizeB;j++) {
            for (int k=1;k<=sizeC;k++) {
              tmpA = (2*i-sizeA-1.0)/2/sizeA;
              tmpB = (2*j-sizeB-1.0)/2/sizeB;
              tmpC = (2*k-sizeC-1.0)/2/sizeC;
              fprintf(job3,"   %f\t%f\t%f\t1.0\n",tmpA,tmpB,tmpC);
            }
          }
        }
      }
      else {
        fprintf(job3," 0.00\t0.00\t0.00\t1.0\n");
      }
      fprintf(job3," } \n");

      fprintf(job3,"}\n");
      fprintf(job3," \n");

      fclose(job3);

      FILE *job4;
      if ((job4 = fopen(filename4.c_str(),"w"))==NULL) {
        printf("Cluster::CreateDFTBJob : Cannot open file '%s'\n",filename4.c_str());
        exit(1);
      }

      fprintf(job4,"Geometry {\n");
      //fprintf(job,"  %i\n",GetTotalNumberOfAtoms()); //<-- Don't need this line
      //Add in Molecular Geometry here...
      fprintf(job4," %s\n",Supercell::supercell().GetDFTBGeomInput().c_str());
      //fprintf(job," }\n");
      fprintf(job4,"  Periodic = Yes\n");
      fprintf(job4,"  LatticeVectors [Angstrom] = {\n");
      //Vector tmp_vec;
      for (int i=0; i<3; i++){
        Vector tmp_vec=Supercell::supercell().GetUnitCellVector(i);
       fprintf(job4,"  %f\t%f\t%f\n",tmp_vec[0],tmp_vec[1],tmp_vec[2]);
      }
      fprintf(job4,"  }\n");
      fprintf(job4,"}\n");
      fprintf(job4,"Hessian = {\n");
      fprintf(job4,"  <<< \"hessian.out\"\n");
      fprintf(job4,"}\n");
      fprintf(job4,"%s\n",Params::Parameters().GetSlaterKoster().c_str());

      fclose(job4);
    }

  }
  fflush(stdout);
  //exit(0);

}

string Cluster::GetDFTBGeomInput() {
  //printf("Attemping to create DFTB geometry input\n");
  int tmp; 
  int atomicNum;
  int lengthArray;
  bool match = false;
  string typeAtom;
  string geomInput = " TypeNames = {";
  
  //printf("geomInput = %s\n",geomInput.c_str());
  
  vector<int> check;
  vector<string> atomTypes;
  int numUniqueAtoms;

  //count up total number of atom types
  for ( int i=1; i<=NMon; i++){
    tmp = Monomers[i].GetNumberOfAtoms();
    for ( int j=0; j<tmp; j++){
      atomicNum = Monomers[i].GetAtom(j).GetAtomicNumber();
      typeAtom = Monomers[i].GetAtom(j).GetSymbol();
      //printf("Monomer %d, atom %d = %d\n",i,j,atomicNum);
      //fflush(stdout);
      lengthArray=check.size();

      if(i==1 && j==0) {
        match=false;
      }
      else{
        for (int k=0; k<lengthArray; k++){
          //printf("check[k] = %d\natomicNum = %d\n",check[k],atomicNum);
          //fflush(stdout);
          if(check[k]==atomicNum){
            match=true;
          }
        }
      }
      if(!match){
        numUniqueAtoms++;
        check.push_back(atomicNum);
        atomTypes.push_back(typeAtom);
      }
      match=false;
    }
  }
  
  //Make the beginning of the atom type list
  for (int k=0; k<atomTypes.size(); k++){
    //printf("k = %i\tatomTypes[k] = %s\n",k,atomTypes[k].c_str());
    geomInput = geomInput + " \"" + atomTypes[k] + "\"";
  }

  //Now construct the actual geometry input...
  geomInput = geomInput + " }\n  TypesAndCoordinates [Angstrom] = {\n";
  for ( int i=1; i<=NMon; i++){
    tmp = Monomers[i].GetNumberOfAtoms();
    for ( int j=0; j<tmp; j++){
      typeAtom = Monomers[i].GetAtom(j).GetSymbol();

      //Locate atom Type in the list of atomTypes...
      for (int k=0; k<atomTypes.size(); k++){
        //printf("check[k] = %d\natomicNum = %d\n",check[k],atomicNum);
        //fflush(stdout);
        if(atomTypes[k]==typeAtom){
          stringstream ss;
          ss << k+1;
          string str = ss.str();
          geomInput = geomInput + "    " + str + "\t";
        }
      }
      Vector tmp_xyz = Monomers[i].GetAtom(j).GetPosition();
      //construct string with current atom Type num, and coordinates
      for( int l=0; l< 3; l++){
        stringstream ss;
        ss << tmp_xyz[l];
        string str = ss.str();
        geomInput = geomInput + str + " ";
      }
      geomInput = geomInput + "\n";
    }
  }
  geomInput = geomInput + "  }\n";
  //printf("geomInput = \n%s\n",geomInput.c_str());
  //fflush(stdout);
  //exit(0);

  return geomInput;

}

void Cluster::CreateFullClusterTinkerJob(int path, bool DiffName, string Name) {
  // path = 1 --> normal mm path
  // path = 2 --> in tmp directory
  // path = 3 --> in hessian directory

/*
  // Full cluster
  string xyzfile, keyfile;
  if ( (path==2) && !DiffName) {
    xyzfile = Params::Parameters().GetTmpMMFilesPath() + "/full.xyz";
    keyfile = Params::Parameters().GetTmpMMFilesPath() + "/full.key";
  }
  if ( (path==2) && DiffName) {
    xyzfile = Params::Parameters().GetTmpMMFilesPath() + "/"+Name+".xyz";
    keyfile = Params::Parameters().GetTmpMMFilesPath() + "/"+Name+".key";
  }
  if (path==3) {
    xyzfile = Params::Parameters().GetHessianMMPath() + "/full.xyz";
    keyfile = Params::Parameters().GetHessianMMPath() + "/full.key";
  }
  else if (path==1) {
    xyzfile = Params::Parameters().GetMMPath() + "/full.xyz"; 
    keyfile = Params::Parameters().GetMMPath() + "/full.key"; 
  }
*/
  // Full cluster
  string pathstring;
  if (path==1) {
    pathstring = Params::Parameters().GetMMPath(); 
  }
  else if ( (path==2)) {
    pathstring = Params::Parameters().GetTmpMMFilesPath();

  }
  else if (path==3) {
    pathstring = Params::Parameters().GetHessianMMPath();
  }
  else{
    printf("Error::Cluster::CreateFullClusterTinkerJob() path number not recognized path = %i\n",path);
    exit(0);
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

 string xyzfile, keyfile;
 if(DiffName){
   xyzfile += pathstring + "/" + Name + ".xyz";
   keyfile += pathstring + "/" + Name + ".key";
 }
 else {
    xyzfile = pathstring + "/full.xyz";
    keyfile = pathstring + "/full.key";
 }

  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Cluster::CreateFullClusterTinkerJob : Cannot open file '%s'\n",
	   xyzfile.c_str());
    exit(1);
  }
  PrintTinkerCartesian(xyz);
  fclose(xyz);


  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
  if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
    printf("Cluster::CreateFullClusterTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }
  fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );

  // If Periodic, add periodicity info - unit cell, ewald summation
  if ( Params::Parameters().IsPeriodic() ) {
    double a = unit_cell[0].Norm();
    double b = unit_cell[1].Norm();
    double c = unit_cell[2].Norm();

    double alpha = RadiansToDegrees*
      acos(unit_cell[1].DotProduct( unit_cell[2] ) / (b*c));
    double beta  = RadiansToDegrees*
      acos(unit_cell[0].DotProduct( unit_cell[2] ) / (a*c));
    double gamma = RadiansToDegrees*
      acos(unit_cell[0].DotProduct( unit_cell[1] ) / (a*b));
    
    printf("Tinker lattice parameters:\n");
    printf("a = %f, b = %f, c = %f\n",a,b,c);
    printf("alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);
    
    fprintf(key,"# Periodic boundary conditions\n");
    fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",a,b,c);
    fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",alpha,beta,gamma);
    fprintf(key,"EWALD\t\tTRUE\n");
    // If we want vacuum boundary conditions for the ewald sum, set
    // EWALD-BOUNDARY flag.  For tinfoil (infinite dielectric) boundary
    // conditions, we omit this keyword.  For non-polar unit cells, it
    // doesn't matter.  But for polar ones, tinfoil boundary conditions
    // seem to work better.
    if (!Params::Parameters().TinFoilBoundaryConditions())
      fprintf(key,"EWALD-BOUNDARY\tTRUE\n");
    
  }

  fclose(key);
}

void Cluster::CreateCrystalJob() {
  // Full cluster
  string pathstring = Params::Parameters().GetMMPath(); ;
  if (Params::Parameters().DoFreq()) {
    pathstring = Params::Parameters().GetHessianMMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  string filename = pathstring + "/full.d12";
  
  // number of atom in the assymmetrical unit
  int Natom = SymmetricalFractCoord.GetLength()/3;

  /* Create the input file */
  FILE *input;
  if ((input = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::CreateCrystalJob : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  //print heading
  fprintf(input,"%s",Params::Parameters().GetCrystalHeading().c_str());

  //print lattice parameters
  fprintf(input,"%f ",UnitCellAxes[0]);
  if(!IsBLocked()) fprintf(input,"%12.10f ",UnitCellAxes[1]);
  if(!IsCLocked()) fprintf(input,"%12.10f ",UnitCellAxes[2]); 
  if(!IsAlphaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[0]);
  if(!IsBetaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[1]);
  if(!IsGammaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[2]);
  fprintf(input,"\n");

  //print atom fractional coordinates
  fprintf(input,"%i\n",Natom);
  for(int iMon = 1; iMon <= NMon; iMon++){
    for(int iAtom = 0; iAtom < Monomers[iMon].GetNumberOfAtoms();iAtom++){
      if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
        fprintf(input,"%-2i  %12.10f  %12.10f  %12.10f\n",
                 Monomers[iMon].GetAtom(iAtom).GetAtomicNumber(),
                 Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(0),
                 Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(1),
                 Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(2));
      }
    }
  } 
  //Force jobs created in CreateFiniteDifferenceCrystalJob
  /*
  if(Params::Parameters().DoForces()){
    fprintf(input,"OPTGEOM\n");
    fprintf(input,"FULLOPTG\n");
    fprintf(input,"MAXCYCLE\n");
    fprintf(input,"1\n");
    fprintf(input,"FINALRUN\n");
    fprintf(input,"4\n");
    fprintf(input,"ENDOPT\n");
  }
  */

  if(Params::Parameters().DoFreq()){
    fprintf(input,"FREQCALC\n");
    fprintf(input,"NOOPTGEOM\n");
    fprintf(input,"MODES\n");
    fprintf(input,"END\n");
  }

  fprintf(input,"END\n");

  //Including Basis set
  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalBasis().c_str());

  //Including Ending 
  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalEnding().c_str());

  fclose(input);

}

void Cluster::CreateFiniteDifferenceCrystalJob(){
  
  //Finite difference step size
  double delta = 0.0002;//in fractional coordinate
  double lengthdelta = 0.001;//in angstroms
  double angledelta = 0.01;//in degrees

  // Full cluster
  string pathstring = Params::Parameters().GetMMPath();
  if (Params::Parameters().DoFreq()) {
    pathstring = Params::Parameters().GetHessianMMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  // number of atom in the assymmetrical unit
  int Natom = SymmetricalFractCoord.GetLength()/3;

  //Nuclear finite difference step
  int i = 0;
  for(int iMon=0;iMon<NMon;iMon++){
    for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
      for(int xyz=0;xyz<3;xyz++){
	bool skip = 0;
	
	if(Params::Parameters().FreezeNuclearCoordinates())
	   skip = 1;
	if(!Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit())
	  skip = 1;
	if(Monomers[iMon].GetAtom(iAtom).IsAtomFrozen())
	  skip = 1;
	if(xyz == 1 && Monomers[iMon].GetAtom(iAtom).IsYLocked())
	  skip = 1;
	if(xyz == 2 && Monomers[iMon].GetAtom(iAtom).IsZLocked())
	  skip = 1;
	if(!skip){
	  //creating string from i for the filenames
	  char num[10];
	  sprintf(num,"%d",i);
	  

	  //Plus geometry
	  string filename = pathstring + "/full+" + num + ".d12";
	  //printf("filename = %s\n",filename.c_str());
	  
	  /* Create the input file */
	  FILE *input;
	  if ((input = fopen(filename.c_str(),"w"))==NULL) {
	    printf("Cluster::CreateCrystalJob : Cannot open file '%s'\n",
		   filename.c_str());
          exit(1);
	  }
	  
	  //print heading
	  fprintf(input,"%s",Params::Parameters().GetCrystalHeading().c_str());
	  
	  //print lattice parameters
	  fprintf(input,"%f ",UnitCellAxes[0]);
	  if(!IsBLocked()) fprintf(input,"%12.10f ",UnitCellAxes[1]);
	  if(!IsCLocked()) fprintf(input,"%12.10f ",UnitCellAxes[2]); 
	  if(!IsAlphaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[0]);
	  if(!IsBetaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[1]);
	  if(!IsGammaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[2]);
	  fprintf(input,"\n");      
	  
	  //Shifting atoms while maintain symmetry
	  //Vector Coord = SymmetryUniqueCoordinates;
	  Vector Coord = SymmetricalFractCoord;
	  Coord[i] += delta;
	  Coord = ConvertBetweenFractionAndCartesianCoordinates(Coord,1,0);
	  Cluster::cluster().MaintainCartesianSymmetry(Coord,false);
	  Coord = ConvertBetweenFractionAndCartesianCoordinates(Coord,0,0); 
	  
	  //print atom fractional coordinates
	  int j=0;
	  fprintf(input,"%i\n",Natom);
	  for(int jMon = 1; jMon <= NMon; jMon++){
	    for(int jAtom = 0; jAtom < Monomers[jMon].GetNumberOfAtoms();jAtom++){
	      if(Monomers[jMon].GetAtom(jAtom).InAsymmetricUnit()){
		fprintf(input,"%-2i  %12.10f  %12.10f  %12.10f\n",
			Monomers[jMon].GetAtom(jAtom).GetAtomicNumber(),
			Coord[3*j],Coord[3*j+1],Coord[3*j+2]);
		j++;
	      }
	    }
	
	  }
	  if(Params::Parameters().DoFreq()){
	    fprintf(input,"OPTGEOM\n");
	    fprintf(input,"FULLOPTG\n");
	    fprintf(input,"MAXCYCLE\n");
	    fprintf(input,"1\n");
	    fprintf(input,"FINALRUN\n");
	    fprintf(input,"4\n");
	    fprintf(input,"ENDOPT\n");
	  }
	  
	  fprintf(input,"END\n");
	  
	  //Including Basis set
	  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalBasis().c_str());
	  
	  //Including Ending 
	  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalEnding().c_str());
	  
	  fclose(input);

	  //Minus geometry
	  filename = pathstring + "/full-" + num + ".d12";
	  //printf("filename = %s\n",filename.c_str());
	  
	  /* Create the input file */
	  if ((input = fopen(filename.c_str(),"w"))==NULL) {
	    printf("Cluster::CreateCrystalJob : Cannot open file '%s'\n",
		   filename.c_str());
          exit(1);
	  }

	  //print heading
	  fprintf(input,"%s",Params::Parameters().GetCrystalHeading().c_str());
	  
	  //print lattice parameters
	  fprintf(input,"%f ",UnitCellAxes[0]);
	  if(!IsBLocked()) fprintf(input,"%12.10f ",UnitCellAxes[1]);
	  if(!IsCLocked()) fprintf(input,"%12.10f ",UnitCellAxes[2]); 
	  if(!IsAlphaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[0]);
	  if(!IsBetaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[1]);
	  if(!IsGammaLocked()) fprintf(input,"%12.10f ",UnitCellAngles[2]);
	  fprintf(input,"\n");      
	  
	  //Shifting atoms while maintain symmetry
	  Coord = SymmetricalFractCoord;
	  Coord[i] -= delta;
	  Coord = ConvertBetweenFractionAndCartesianCoordinates(Coord,1,0);
	  Cluster::cluster().MaintainCartesianSymmetry(Coord,false);
	  Coord = ConvertBetweenFractionAndCartesianCoordinates(Coord,0,0); 

	  //print atom fractional coordinates
	  j=0;
	  fprintf(input,"%i\n",Natom);
	  for(int jMon = 1; jMon <= NMon; jMon++){
	    for(int jAtom = 0; jAtom < Monomers[jMon].GetNumberOfAtoms();jAtom++){
	      if(Monomers[jMon].GetAtom(jAtom).InAsymmetricUnit()){
		fprintf(input,"%-2i  %12.10f  %12.10f  %12.10f\n",
			Monomers[jMon].GetAtom(jAtom).GetAtomicNumber(),
			Coord[3*j],Coord[3*j+1],Coord[3*j+2]);
		j++;
	      }
	    }
	
	  }
	  if(Params::Parameters().DoFreq()){
	    fprintf(input,"OPTGEOM\n");
	    fprintf(input,"FULLOPTG\n");
	    fprintf(input,"MAXCYCLE\n");
	    fprintf(input,"1\n");
	    fprintf(input,"FINALRUN\n");
	    fprintf(input,"4\n");
	    fprintf(input,"ENDOPT\n");
	  }
	  
	  fprintf(input,"END\n");
	  
	  //Including Basis set
	  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalBasis().c_str());
	  
	  //Including Ending 
	  fprintf(input,"%sEND\n",Params::Parameters().GetCrystalEnding().c_str());
	  
	  fclose(input);
	}
	i++;

      }

    }
  }

  //Lattice parameters finite stepsize
  for(int i=0;i<6;i++){
    bool skip = 0;
    
    if(Params::Parameters().FreezeUnitCellParams())
      skip = 1;
    if(i==1 && IsBLocked())
      skip = 1;
    if(i==2 && IsCLocked())
      skip = 1;
    if(i==3 && IsAlphaLocked())
      skip = 1;
    if(i==4 && IsBetaLocked())
      skip = 1;
    if(i==5 && IsGammaLocked())
      skip = 1;

    if(!skip) {
      
      //Length or angle string
      string LatticeType;
      if(i==0)
	LatticeType = "a";
      else if(i==1)
	LatticeType = "b";
      else if(i==2)
	LatticeType = "c";
      else if(i==3)
	LatticeType = "alpha";
      else if(i==4)
      LatticeType = "beta";
      else if(i==5)
	LatticeType = "gamma";
      
      //Plus geometry
      string filename = pathstring + "/full+" + LatticeType + ".d12";
      //printf("filename = %s\n",filename.c_str());
      
      /* Create the input file */
      FILE *input;
      if ((input = fopen(filename.c_str(),"w"))==NULL) {
	printf("Cluster::CreateCrystalJob : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }
      
      //print heading
      fprintf(input,"%s",Params::Parameters().GetCrystalHeading().c_str());
      
      //Altering Lattice Params
      Vector LatticeValues(6);
      for(int j=0;j<3;j++){
	LatticeValues[j] = UnitCellAxes[j];
	LatticeValues[j+3] = UnitCellAngles[j];
      }
      if(i < 3)
	LatticeValues[i] += lengthdelta;
      else
	LatticeValues[i] += angledelta;
      //print lattice parameters
      fprintf(input,"%f ",LatticeValues[0]);
      if(!IsBLocked()) fprintf(input,"%12.10f ",LatticeValues[1]);
      if(!IsCLocked()) fprintf(input,"%12.10f ",LatticeValues[2]); 
      if(!IsAlphaLocked()) fprintf(input,"%12.10f ",LatticeValues[3]);
      if(!IsBetaLocked()) fprintf(input,"%12.10f ",LatticeValues[4]);
      if(!IsGammaLocked()) fprintf(input,"%12.10f ",LatticeValues[5]);
      fprintf(input,"\n");

      //print atom fractional coordinates
      fprintf(input,"%i\n",Natom);
      for(int iMon = 1; iMon <= NMon; iMon++){
	for(int iAtom = 0; iAtom < Monomers[iMon].GetNumberOfAtoms();iAtom++){
	  if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	    fprintf(input,"%-2i  %12.10f  %12.10f  %12.10f\n",
		    Monomers[iMon].GetAtom(iAtom).GetAtomicNumber(),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(0),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(1),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(2));
	  }
	}
      } 
      
      if(Params::Parameters().DoFreq()){
	fprintf(input,"OPTGEOM\n");
	fprintf(input,"FULLOPTG\n");
	fprintf(input,"MAXCYCLE\n");
	fprintf(input,"1\n");
	fprintf(input,"FINALRUN\n");
	fprintf(input,"4\n");
	fprintf(input,"ENDOPT\n");
      }
      
      fprintf(input,"END\n");
      //Including Basis set
      fprintf(input,"%sEND\n",Params::Parameters().GetCrystalBasis().c_str());
      //Including Ending 
      fprintf(input,"%sEND\n",Params::Parameters().GetCrystalEnding().c_str());
      fclose(input);  

      //Minus geometry
      filename = pathstring + "/full-" + LatticeType + ".d12";
      //printf("filename = %s\n",filename.c_str());
      
      /* Create the input file */
      if ((input = fopen(filename.c_str(),"w"))==NULL) {
	printf("Cluster::CreateCrystalJob : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }
      
      //print heading
      fprintf(input,"%s",Params::Parameters().GetCrystalHeading().c_str());
      
      //Altering Lattice Params
      for(int j=0;j<3;j++){
	LatticeValues[j] = UnitCellAxes[j];
	LatticeValues[j+3] = UnitCellAngles[j];
      }
      if(i < 3)
	LatticeValues[i] -= lengthdelta;
      else
	LatticeValues[i] -= angledelta;
      //print lattice parameters
      fprintf(input,"%f ",LatticeValues[0]);
      if(!IsBLocked()) fprintf(input,"%12.10f ",LatticeValues[0]);
      if(!IsCLocked()) fprintf(input,"%12.10f ",LatticeValues[1]); 
      if(!IsAlphaLocked()) fprintf(input,"%12.10f ",LatticeValues[2]);
      if(!IsBetaLocked()) fprintf(input,"%12.10f ",LatticeValues[3]);
      if(!IsGammaLocked()) fprintf(input,"%12.10f ",LatticeValues[4]);
      fprintf(input,"\n");

      //print atom fractional coordinates
      fprintf(input,"%i\n",Natom);
      for(int iMon = 1; iMon <= NMon; iMon++){
	for(int iAtom = 0; iAtom < Monomers[iMon].GetNumberOfAtoms();iAtom++){
	  if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	    fprintf(input,"%-2i  %12.10f  %12.10f  %12.10f\n",
		    Monomers[iMon].GetAtom(iAtom).GetAtomicNumber(),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(0),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(1),
		    Monomers[iMon].GetAtom(iAtom).GetFractionalPosition(2));
	  }
	}
      } 
      
      if(Params::Parameters().DoFreq()){
	fprintf(input,"OPTGEOM\n");
	fprintf(input,"FULLOPTG\n");
	fprintf(input,"MAXCYCLE\n");
	fprintf(input,"1\n");
	fprintf(input,"FINALRUN\n");
	fprintf(input,"4\n");
	fprintf(input,"ENDOPT\n");
      }
      
      fprintf(input,"END\n");
      //Including Basis set
      fprintf(input,"%sEND\n",Params::Parameters().GetCrystalBasis().c_str());
      //Including Ending 
      fprintf(input,"%sEND\n",Params::Parameters().GetCrystalEnding().c_str());
      fclose(input);  
    } 

  }
}


double Cluster::GetEnthalpyPV() {

  double pressure = Params::Parameters().GetExternalPressure() * GigaPascalsToHartreesperBohrcube; //pressure in a.u.

  double volume = UnitCellVolume() * AngToBohr * AngToBohr * AngToBohr;

  double enthalpyPV = pressure*volume;    

  cout << "enthalpyPV = " << enthalpyPV*2625.5 << " KJ/mol \n";
  return enthalpyPV; 

}   

//Use if gradients or hessians are found with finite difference
void Cluster::CreateFiniteDifferenceTinkerJobs(int path) {
  // path = 1 --> normal mm path
  // path = 2 --> in tmp directory
  // path = 3 --> in hessian directory


  string pathstring;
  if (path==1) {
     pathstring = Params::Parameters().GetMMPath();
  }
  else if (path==2) {
    pathstring = Params::Parameters().GetTmpMMFilesPath();
  }
  else if (path==3) {
    pathstring = Params::Parameters().GetHessianMMPath();
  }
  else{
    printf("Error::Cluster::CreateFiniteDifferenceTinkerJobs() path number not recognized path = %i\n",path);
    exit(0);
  }



  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  int Natoms = GetTotalNumberOfAtoms();
  for(int i=0;i<3*Natoms;i++){

    //creating string from i for the filenames
    char num[10];
    sprintf(num,"%d",i);

  // The plus input file
    string xyzfile = pathstring + "/full/full" + num + "+.xyz";
    string keyfile = pathstring + "/full/full" + num + "+.key";;

/*
    // The plus input file
    string xyzfile, keyfile;
    if (path==2) {
      xyzfile = Params::Parameters().GetTmpMMFilesPath() + "/full/full" + num + "+.xyz";
      keyfile = Params::Parameters().GetTmpMMFilesPath() + "/full/full" + num + "+.key";
    }
    if (path==3) {
      xyzfile = Params::Parameters().GetHessianMMPath() + "/full/full" + num + "+.xyz";
      keyfile = Params::Parameters().GetHessianMMPath() + "/full/full" + num + "+.key";
    }
    else if (path==1) {CreateFiniteDifferenceCrystalJob
      xyzfile = Params::Parameters().GetMMPath() + "/full/full" + num + "+.xyz"; 
      keyfile = Params::Parameters().GetMMPath() + "/full/full" + num + "+.key"; 
    }
    //xyzfile += num + "+.xyz";
    //keyfile += num + "+.key";

*/



    //printf("xyz = %s\n",xyzfile.c_str());
    //printf("key = %s\n",keyfile.c_str());
    //exit(0);

    /* Create the xyz file */
    FILE *xyz;
    if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateFullClusterTinkerJobForFiniteDifference : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    
    //Set Coordiantes of file
    Vector Coords = GetCurrentCoordinates(false);
    Coords[i] -= 0.001;
    PrintTinkerCartesian(Coords,xyz);
    fclose(xyz);
    /* Create the keyfile */
    // Open the file for writing, write the Tinker rem section to it,
    // and close the file.
    FILE *key;
    if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateFullClusterTinkerJobForFiniteDifferenc : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
    
    // If Periodic, add periodicity info - unit cell, ewald summation
    if ( Params::Parameters().IsPeriodic() ) {
      double a = unit_cell[0].Norm();
      double b = unit_cell[1].Norm();
      double c = unit_cell[2].Norm();
      
      double alpha = RadiansToDegrees*
	acos(unit_cell[1].DotProduct( unit_cell[2] ) / (b*c));
      double beta  = RadiansToDegrees*
	acos(unit_cell[0].DotProduct( unit_cell[2] ) / (a*c));
      double gamma = RadiansToDegrees*
	acos(unit_cell[0].DotProduct( unit_cell[1] ) / (a*b));
      
      //printf("Tinker lattice parameters:\n");
      //printf("a = %f, b = %f, c = %f\n",a,b,c);
      //printf("alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);
    
      fprintf(key,"# Periodic boundary conditions\n");
      fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",a,b,c);
      fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",alpha,beta,gamma);
      fprintf(key,"EWALD\t\tTRUE\n");
      // If we want vacuum boundary conditions for the ewald sum, set
      // EWALD-BOUNDARY flag.  For tinfoil (infinite dielectric) boundary
      // conditions, we omit this keyword.  For non-polar unit cells, it
      // doesn't matter.  But for polar ones, tinfoil boundary conditions
      // seem to work better.
      if (!Params::Parameters().TinFoilBoundaryConditions())
	fprintf(key,"EWALD-BOUNDARY\tTRUE\n"); 
    }
    fclose(key);

    // The minus input file
    xyzfile = pathstring + "/full/full" + num + "-.xyz";
    keyfile = pathstring + "/full/full" + num + "-.key"; 

/*
    // The minus input file
    if (path==2) {
      xyzfile = Params::Parameters().GetTmpMMFilesPath() + "/full/full" + num + "-.xyz";
      keyfile = Params::Parameters().GetTmpMMFilesPath() + "/full/full" + num + "-.key"; 

    }
    if (path==3) {
      xyzfile = Params::Parameters().GetHessianMMPath() + "/full/full" + num + "-.xyz";
      keyfile = Params::Parameters().GetHessianMMPath() + "/full/full" + num + "-.key"; 
    }
    else if (path==1) {
      xyzfile = Params::Parameters().GetMMPath() + "/full/full" + num + "-.xyz"; 
      keyfile = Params::Parameters().GetMMPath() + "/full/full" + num + "-.key"; 
    }
    //xyzfile += num + "-.xyz";
    //keyfile += num + "-.key";
*/

    if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateFullClusterTinkerJobForFiniteDifference : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }

    //Set Coordinates of file
    Coords[i] += 0.002;
    PrintTinkerCartesian(Coords,xyz);
    fclose(xyz);    

   /* Create the keyfile */
    // Open the file for writing, write the Tinker rem section to it,
    // and close the file.
    if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateTinkerJob : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );

    // If Periodic, add periodicity info - unit cell, ewald summation
    if ( Params::Parameters().IsPeriodic() ) {
      double a = unit_cell[0].Norm();
      double b = unit_cell[1].Norm();
      double c = unit_cell[2].Norm();
      
      double alpha = RadiansToDegrees*
	acos(unit_cell[1].DotProduct( unit_cell[2] ) / (b*c));
      double beta  = RadiansToDegrees*
	acos(unit_cell[0].DotProduct( unit_cell[2] ) / (a*c));
      double gamma = RadiansToDegrees*
	acos(unit_cell[0].DotProduct( unit_cell[1] ) / (a*b));
      
      //printf("Tinker lattice parameters:\n");
      //printf("a = %f, b = %f, c = %f\n",a,b,c);
      //printf("alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);
    
      fprintf(key,"# Periodic boundary conditions\n");
      fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",a,b,c);
      fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",alpha,beta,gamma);
      fprintf(key,"EWALD\t\tTRUE\n");
      // If we want vacuum boundary conditions for the ewald sum, set
      // EWALD-BOUNDARY flag.  For tinfoil (infinite dielectric) boundary
      // conditions, we omit this keyword.  For non-polar unit cells, it
      // doesn't matter.  But for polar ones, tinfoil boundary conditions
      // seem to work better.
      if (!Params::Parameters().TinFoilBoundaryConditions())
	fprintf(key,"EWALD-BOUNDARY\tTRUE\n"); 
    }
    fclose(key);
  }
  
}

void Cluster::CreateFullClusterQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.force'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  


  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string filename = path + "/full.in"; 

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::CreateFullClusterQChemJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  // Print comment section
  fprintf(job,"$comment\nFull cluster\n$end\n\n");

  // Print $molecule section
  PrintQChemCartesian(job);
  
  // Print $rem section
  string rem;
  string basis;
  if (MM_job)
    rem = Params::Parameters().GetQChemRem2();
  else
    rem = Params::Parameters().GetQChemRem();
  basis = Params::Parameters().GetQChemBasis();

  // for CP correction, insert jobtype = BSSE at end of rem section.  
  // Note: This will override any previous setting for jobtype.
  if ( Params::Parameters().DoCounterpoise() ) {
    int spot = rem.find("$end");
    string bsse = "jobtype = bsse\n";
    rem.insert(spot,bsse);
    //printf("New rem:\n%s\n",rem.c_str());
  }
  
  fprintf(job,"%s\n",rem.c_str());
  fprintf(job,"%s\n",basis.c_str()); // by shuhao basis
  fclose(job);
}

// Manages the running of jobs, either in serial or parallel
void Cluster::RunJobs() {

  
  // NOTE: debug: need to test that MM force jobs work properly

  //Number of QM and MM Monomer Jobs depends whether monomer symmetry is being used
  int QM_Mon = NMon;
  if(Params::Parameters().UseMonomerSymmetry())
    QM_Mon = UniqueMon;

  /* Step 1: create list of jobs to run */
  // Count the jobs
  int N_QM_jobs = QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc;


  //MolPro Frequency calculations are determined entirely by 
  //if(Params::Parameters().DoCBS() && Params::Parameters().DoFreq() ){
  if((Params::Parameters().GetQMType() == 2) && Params::Parameters().DoFreq()){
    N_QM_jobs = 0;


    if(Params::Parameters().DoEnergyFiniteDifferenceFreqs()){
      for(int iMon = 1; iMon <= NMon; iMon++){
	if((Monomers[iMon].GetSymmetryFactor() != 0) || !Params::Parameters().UseMonomerSymmetry()){
	  N_QM_jobs += 6*Monomers[iMon].GetNumberOfAtoms();
	}
      }

      double c0 = Params::Parameters().GetLocalCutoff(0);
      double c1 = Params::Parameters().GetLocalCutoff(1);
      for(int i=1;i<=NDim;i++){
	double separation =  Dimers[i].GetDimerSeparation();  
	if(Dimers[i].GetSymmetryFactor() != 0 &&
	   (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	  N_QM_jobs += 2*(3*Dimers[i].GetNumberOfAtoms())*(3*Dimers[i].GetNumberOfAtoms()+1);
	}
      }

      for(int i=1;i<=NDim_images;i++){
	double separation =  DimerImages[i].GetDimerSeparation();
	if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	  N_QM_jobs += 2*(3*DimerImages[i].GetNumberOfAtoms())*(3*DimerImages[i].GetNumberOfAtoms()+1);
	}
      }

    }
    else{
      if(Params::Parameters().SingleFileMonomerHess()){
	N_QM_jobs += QM_Mon;

      }else{
	for(int iMon = 1; iMon <= NMon; iMon++){
	  if((Monomers[iMon].GetSymmetryFactor() != 0) || !Params::Parameters().UseMonomerSymmetry()){
	    N_QM_jobs += 6*Monomers[iMon].GetNumberOfAtoms();
	  }
	}
      } 
   
      double c0 = Params::Parameters().GetLocalCutoff(0);
      double c1 = Params::Parameters().GetLocalCutoff(1);
      for(int i=1;i<=NDim;i++){
	double separation =  Dimers[i].GetDimerSeparation();
	if(Dimers[i].GetSymmetryFactor() != 0 &&
	   (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	  N_QM_jobs += 6*Dimers[i].GetNumberOfAtoms();
	}
	
      }
      for(int i=1;i<=NDim_images;i++){
	double separation =  DimerImages[i].GetDimerSeparation();
	if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	  N_QM_jobs += 6*DimerImages[i].GetNumberOfAtoms();
	}
      }
    }
    
    if(Params::Parameters().DoCCSDTCorrection()){
      for(int iMon = 1; iMon <= NMon; iMon++){
	 if((Monomers[iMon].GetSymmetryFactor() != 0) || !Params::Parameters().UseMonomerSymmetry()){
	   N_QM_jobs += 6*Monomers[iMon].GetNumberOfAtoms();
	 }
       }
	
       double c0 = Params::Parameters().GetLocalCutoff(0);
       double c1 = Params::Parameters().GetLocalCutoff(1);
       for(int i=1;i<=NDim;i++){
	 double separation =  Dimers[i].GetDimerSeparation();  
	 if(Dimers[i].GetSymmetryFactor() != 0 &&
	   (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	    N_QM_jobs += 2*(3*Dimers[i].GetNumberOfAtoms())*(3*Dimers[i].GetNumberOfAtoms()+1);
	  }
	}
	
	for(int i=1;i<=NDim_images;i++){
	  double separation =  DimerImages[i].GetDimerSeparation();
	  if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	    N_QM_jobs += 2*(3*DimerImages[i].GetNumberOfAtoms())*(3*DimerImages[i].GetNumberOfAtoms()+1);
	  }
        }
      }
     }
     else{
      //N_QM_jobs *= 2;
       if(Params::Parameters().DoCCSDTCorrection()){
         N_QM_jobs += QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc;
	

         //CCSD(T) counterpoise dimer forces are handle in seperate finite difference jobs
         if(Params::Parameters().DoForces()){
           double c0 = Params::Parameters().GetLocalCutoff(0);
           double c1 = Params::Parameters().GetLocalCutoff(1);
           for(int i=1;i<=NDim;i++){
	     double separation =  Dimers[i].GetDimerSeparation();
	     if(Dimers[i].GetSymmetryFactor() != 0 &&
	       (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	         N_QM_jobs += 6*Dimers[i].GetNumberOfAtoms();
	      }
	    } 
	    for(int i=1;i<=NDim_images;i++){
	      double separation =  DimerImages[i].GetDimerSeparation();
	      if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	        N_QM_jobs += 6*DimerImages[i].GetNumberOfAtoms();
	      }
	    }
	 }
      }
   }
  //int N_MM_jobs = N_QM_jobs + 1; 
  // unless MM is AIFF, do all QM jobs 
  int N_MM_jobs = QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc +1;
  if ( Params::Parameters().GetMMType() == 2 )// for AIFF only do monomer QM jobs
    if(Params::Parameters().UseMMSymmetry())
      N_MM_jobs = UniqueMon; //AIFF exploiting symmetry
    else
      N_MM_jobs = NMon; // AIFF not exploiting symmetry

  //only performing the full MM
  if(Params::Parameters().UseFullMMOnly() && Params::Parameters().GetMMType() == 1 ){
    N_MM_jobs = 1;
  }

  //only performing the full QM
  if(Params::Parameters().UseFullQMOnly()){
    N_MM_jobs = 0;
    N_QM_jobs = 1;

  }
  
  //Crystal
  if(Params::Parameters().GetMMType() == 5){
    N_MM_jobs = 1;
    //Finite Difference for force calculations
    if(Params::Parameters().DoForces()){
 
      //Finite Difference for Nuclear coordinates  
      //Only for atoms in the asymetrical unit and coordinates
      //not fixed by symmetry
      if(!Params::Parameters().FreezeNuclearCoordinates()){
	for(int iMon=0;iMon<NMon;iMon++){
	  for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
	    if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){

	      if(!Monomers[iMon].GetAtom(iAtom).IsAtomFrozen()){
		if(Monomers[iMon].GetAtom(iAtom).IsYLocked() &&
		   Monomers[iMon].GetAtom(iAtom).IsZLocked() != 0)
		  N_MM_jobs += 2;
		else if(Monomers[iMon].GetAtom(iAtom).IsYLocked() ||
			Monomers[iMon].GetAtom(iAtom).IsZLocked() != 0)
		  N_MM_jobs += 4;
		else
		  N_MM_jobs += 6;
	      }   
	    }
	  }
	}
      }
      //Finite Difference for Lattice Parameters
      if(!Params::Parameters().FreezeUnitCellParams()){
	//Finite Diff for length A
	N_MM_jobs += 2;
	if(!IsBLocked())//Finite Diff for length B
	  N_MM_jobs += 2;
	if(!IsCLocked()) //Finite Diff for length C
	  N_MM_jobs += 2;
	if(!IsAlphaLocked())//Finite Diff for angle alpha
	  N_MM_jobs += 2;
	if(!IsBetaLocked())//Finite Diff for angle beta
	  N_MM_jobs += 2;
	if(!IsGammaLocked())//Finite Diff for angle gamma
	  N_MM_jobs += 2;
      }
    }
    //printf("N_MM_job = %i\n",N_MM_jobs);

  }


  //If full tinker hessian is found by finite difference, run them later on.
  if( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() 
      && Params::Parameters().Do_fdTinkerHessian() && Params::Parameters().GetMMType() == 1){
    //N_MM_jobs += 6*GetTotalNumberOfAtoms() - 1;
    N_MM_jobs += -1;
    //printf("N_MM_jobs = %i\n",N_MM_jobs);
  }

  
  if ( Params::Parameters().DoQMBenchmark() )
    N_QM_jobs++;

  int Njobs = N_QM_jobs + N_MM_jobs;
  printf("QM jobs = %d, MM = %d, total = %d\n",N_QM_jobs,N_MM_jobs,Njobs);

  // Create and Combine the job lists
  string *JobList = new string[Njobs];
  string *QMList = NULL;// = new string[N_QM_jobs];
  string *MMList = NULL;// = new string[N_MM_jobs];

  if (Params::Parameters().BuildForceFieldOnly() ||
      Params::Parameters().UseFullMMOnly()) {
    printf("Skipping QM jobs\n");
    Njobs -= N_QM_jobs;
    N_QM_jobs = 0;
  }
  else {
    QMList = RunQMJobs(N_QM_jobs);
    for (int i=0;i<N_QM_jobs;i++) {
      JobList[i] = QMList[i];    
    }
  }

  if (Params::Parameters().NeglectManyBody() || Params::Parameters().UseFullQMOnly()) {
    printf("Many-body terms being neglected - Skipping MM jobs.\n");
    Njobs -= N_MM_jobs;
  }
  else {
    MMList = RunMMJobs();
    //printf("MMList =%i\n", sizeof(MMList));
    for (int i=0;i<N_MM_jobs;i++) {
      JobList[N_QM_jobs+i] = MMList[i];
    }
  }

  delete [] QMList;
  if (! Params::Parameters().NeglectManyBody() )
    delete [] MMList;
  /*
  printf("List of Jobs to run:\n");
  for (int q=0;q<Njobs;q++) {
    printf("%d: %s\n",q,JobList[q].c_str());
  }
  */


  /* Step 2: Now actually run the jobs */

#ifdef PARALLEL

  int Nproc = Params::Parameters().GetNumberOfProcessors();
  if (Nproc > 1 ){
    /* Parallel  Version */
   printf("The jobs will be run in parallel on %d slave processors\n",Nproc-1);


    // may want Nproc + 1 nodes in code, to use all Nproc efficiently
    // for the actual calculations

    // load up MPI data locally
    MPI_Status status;
    int mynode,totalnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode

    // Verify that we don't have too many processors (rarely a problem!)
    int ExpensiveJobs = N_QM_jobs;
    if ( Params::Parameters().GetMMType() == 2) {
      //ExpensiveJobs += NMon; // AIFF monomer jobs
      ExpensiveJobs += N_MM_jobs;
    }
    if ( Params::Parameters().GetMMType() == 3) 
      ExpensiveJobs += N_MM_jobs; // Secondary QChem jobs are assumed to be expensive
    if ( totalnodes > ExpensiveJobs + 1) { // +1 accounts for master node
      printf("ERROR: Nodes assigned (%d) exceeds number of QM jobs (%d).",
	     totalnodes,ExpensiveJobs);
      printf("  Exiting to avoid wasting precious processors!\n");
      exit(1);
    }

    // Create a few variables
    int ijob = 0; // tracks job number
    int success = 0; // outcome of jobs;
    int tag, rank; // tag = type of job, rank = processor number
    int result; 


    // Seed a job to each processor initially
    for (rank=1;rank<totalnodes;rank++) {

      // Set the tag
      if (ijob < N_QM_jobs) {
	if(Params::Parameters().GetQMType() == 1){
	  //separate tags for QCHEM with and without counterpoise correction
	  if(Params::Parameters().DoCounterpoise())
	    tag = QCHEMCP_TAG;
	  else
	    tag = QCHEM_TAG;
	}
	else if(Params::Parameters().GetQMType() == 2)
	  tag = MOLPRO_TAG;
        else if(Params::Parameters().GetQMType() == 3)
	  tag = G09_TAG;
	else if(Params::Parameters().GetQMType() == 4)
	  tag = DALTON_TAG;
	else if(Params::Parameters().GetQMType() == 5)
	  tag = QUANTUM_ESPRESSO_TAG;
	else if(Params::Parameters().GetQMType() == 6)
	  tag = ORCA_TAG;
	else if(Params::Parameters().GetQMType() == 7)
	  tag = PSI4_TAG;
	else if(Params::Parameters().GetQMType() == 8)
	  tag = DFTB_TAG;
	else {
	  printf("ERROR: Cluster::RunJobs()- Unknown QM_Type = %d\n",
		 Params::Parameters().GetQMType());
	  exit(1);
	}
      }else {
	if (Params::Parameters().GetMMType() == 1 )
	  tag = TINKER_TAG;
	else if (Params::Parameters().GetMMType() == 2 )
	  tag = CAMCASP_TAG;
	else if (Params::Parameters().GetMMType() == 3 )
	  //separate tags for QCHEM with and without counterpoise correction
	  if(Params::Parameters().DoCounterpoise())
	    tag = QCHEMCP_TAG;
	  else
	    tag = QCHEM_TAG;
	else {
	printf("ERROR: Cluster::RunJobs() - Unknown MM_Type = %d\n",
	       Params::Parameters().GetMMType());
	exit(1);
	}
      }
      
      if (Params::Parameters().TinkerDebug())
	tag = TINKER_TAG;
      //UpdateJobStatus(ijob);

      char job[BUFFSIZE]; // 
      sprintf(job,"%s",JobList[ijob].c_str());


      MPI_Send(&job,BUFFSIZE,MPI_CHAR,rank,tag,MPI_COMM_WORLD);
      ijob++;
    }

    int jobs_run = 0;


    // Now continue running jobs until all are done
    while (ijob < Njobs) {
      MPI_Recv(&result,1,MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD,&status);
      
      success += result;

      // Set the tag for the next job
      if (ijob < N_QM_jobs){
	//separate tag for Qchem with and without counterpoise correction
        if(Params::Parameters().GetQMType() == 1){

	  if(Params::Parameters().DoCounterpoise())
	    tag = QCHEMCP_TAG;
	  else
	    tag = QCHEM_TAG;
        }
        else if(Params::Parameters().GetQMType() == 2)
 	  tag = MOLPRO_TAG;
	else if(Params::Parameters().GetQMType() == 3)
	  tag = G09_TAG;
	else if(Params::Parameters().GetQMType() == 4)
	  tag = DALTON_TAG;
	else if(Params::Parameters().GetQMType() == 5)
	  tag = QUANTUM_ESPRESSO_TAG;
	else if(Params::Parameters().GetQMType() == 6)
	  tag = ORCA_TAG;
        else if(Params::Parameters().GetQMType() == 7)
	  tag = PSI4_TAG;
        else if(Params::Parameters().GetQMType() == 8)
	  tag = DFTB_TAG;
	else {
	  printf("ERROR: Cluster::RunJobs() - Unknown QM_Type = %d\n",
		 Params::Parameters().GetQMType());
	  exit(1);        
        }
      }else {
	if (Params::Parameters().GetMMType() == 1 )
	  tag = TINKER_TAG;
	else if (Params::Parameters().GetMMType() == 2 )
	  tag = CAMCASP_TAG;
	else if (Params::Parameters().GetMMType() == 3 )
	  //separate tag for with and without counterpoise correction
	  if(Params::Parameters().DoCounterpoise())
	    tag = QCHEMCP_TAG;
	  else
	    tag = QCHEM_TAG;
	else if(Params::Parameters().GetMMType() == 5)
	  tag = CRYSTAL_TAG;
	else {
	  printf("ERROR: Cluster::RunJobs() - Unknown MM_Type = %d\n",
		 Params::Parameters().GetMMType());
	  exit(1);
	}
      }
      if (Params::Parameters().TinkerDebug())
	tag = TINKER_TAG;
      //UpdateJobStatus(ijob);

      UpdateJobStatus(jobs_run);
      jobs_run++;


      char job[BUFFSIZE];
      sprintf(job,"%s",JobList[ijob].c_str());

	
      MPI_Send(&job,BUFFSIZE,MPI_CHAR, status.MPI_SOURCE,tag,MPI_COMM_WORLD);
      ijob++;

    }

    // Collect all remaining results
    for (int rank = 1; rank < totalnodes; rank++) {
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);


      UpdateJobStatus(jobs_run);
      jobs_run++;

      success += result;
    }
    
    printf("Done running parallel jobs\n");
    if (success == Njobs)
      printf("All %d jobs completed successfully\n", success);
    else {
      printf("Error: Only %d of %d jobs completed successfully\n",success,Njobs);
      printf("Correct failed jobs by hand, and run the HMBI code again with option\n ANALYZE_ONLY = TRUE to continue processing the results\n");
      exit(1);
    }
  } 
  else {
#endif /* PARALLEL */

    /* Serial  Version */
    int ijob = 0;

    time_t QM_start_time, QM_stop_time;
    time_t MM_start_time, MM_stop_time;
    
    QM_start_time = time(NULL); // start the QM timer

    while (ijob < Njobs) {
      if (ijob==N_QM_jobs) {
	QM_stop_time = time(NULL);
	MM_start_time = time(NULL); 
      }
      //printf("execute job %d: %s \n",ijob,JobList[ijob].c_str());
      //fflush(stdout);
      system(JobList[ijob].c_str());
      UpdateJobStatus(ijob);
      ijob++;
    }

    // If we have no MM jobs, above timer stop/start never occurs.
    // So do it here:
    if (Params::Parameters().NeglectManyBody() ) {
      QM_stop_time = time(NULL);
      MM_start_time = time(NULL); 
    }

    MM_stop_time = time(NULL);
    double QM_time = difftime(QM_stop_time,QM_start_time);
    double MM_time = difftime(MM_stop_time,MM_start_time);

    printf("\nTime spent on QM jobs: %0.f sec.    MM jobs: %0.f sec\n",QM_time,MM_time);


#ifdef PARALLEL
  }
#endif /* PARALLEL */
  
  if ( Params::Parameters().GetMMType() == 2  && Params::Parameters().OrientDebug() ) { // Ali
    system(RunOrientJob().c_str());
  }
  NDim_trunc = 0;  
  NDimImages_trunc = 0;
  delete [] JobList;
}

// Returns a list of commands for running each job
string* Cluster::RunQMJobs(int Njobs) {

 //Number of QM Monomer Jobs depends whether monomer symmetry is being used
 bool UseMonSym = Params::Parameters().UseMonomerSymmetry();
 int QM_Mon = NMon;
 if(Params::Parameters().UseMonomerSymmetry()){
   QM_Mon = UniqueMon;
 }

 string *JobList = new string[Njobs];
  
 int ijob = 0;
 if(Params::Parameters().UseFullQMOnly()){
  if(Params::Parameters().GetQMType()==5)
    JobList[ijob] = RunQuantumEspressoJob();
  else if (Params::Parameters().GetQMType()==8)
    JobList[ijob] = RunDFTBJob();
  ijob++;   
 }
 else {
  // Monomers
  for (int i=1;i<=NMon;i++) {
    if (Params::Parameters().TinkerDebug()){
      JobList[ijob] = Monomers[i].RunMMJob();
      ijob++;
    }else {
	 //( Params::Parameters().DoCBS() || Params::Parameters().DoCCSDTCorrection() )){
      if(Params::Parameters().DoFreq() && ( Params::Parameters().GetQMType() == 2 )){
	if(!( UseMonSym && (Monomers[i].GetSymmetryFactor() == 0) )){//yoni: if monomer symmetry and  sym_fac = 0, do not run
	  if(Params::Parameters().SingleFileMonomerHess()){
	    JobList[ijob] = Monomers[i].RunQMJob();
	    ijob++;
	    for(int j=0;j<3*Monomers[i].GetNumberOfAtoms();j++){
	      if(Params::Parameters().DoCCSDTCorrection()){
		JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,true,true);
		ijob++;
		JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,false,true);
		ijob++;
	      }
	    }
	  }
	  else{
	    for(int j=0;j<3*Monomers[i].GetNumberOfAtoms();j++){
	      //if(Params::Parameters().DoCBS()){
	      JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,true,false);
	      ijob++;
	      JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,false,false);
	      ijob++;
	      //}
	      if(Params::Parameters().DoCCSDTCorrection()){
		JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,true,true);
		ijob++;
		JobList[ijob] = Monomers[i].RunFiniteDifferenceMolProJobs(j,false,true);
		ijob++;
	      }
	    } 
	  }
	}	
      }
      else{
	if(UseMonSym){//if monomer symmetry is being used
	  if(Monomers[i].GetSymmetryFactor() != 0){//yoni: if sym_fac = 0, do not run
	    JobList[ijob] = Monomers[i].RunQMJob();
	    ijob++;
	    if(Params::Parameters().DoCCSDTCorrection()){
	      JobList[ijob] = Monomers[i].RunQMJob(true);
	      ijob++;
	    }
	    //JobList[ijob] = Monomers[i].RunQChemJob();
	  }
	}else{//no monomer symmetry
	  JobList[ijob] = Monomers[i].RunQMJob();
	  ijob++;
	  if(Params::Parameters().DoCCSDTCorrection()){
	    JobList[ijob] = Monomers[i].RunQMJob(true);
	    ijob++;
	  }
	//JobList[ijob] = Monomers[i].RunQChemJob();
	}
      }
    }
  }

  double c0 = Params::Parameters().GetLocalCutoff(0);
  double c1 = Params::Parameters().GetLocalCutoff(1);

  // Dimers
  for (int i=1;i<=NDim;i++) {
    double separation =  Dimers[i].GetDimerSeparation();
   
    if(Dimers[i].GetSymmetryFactor()!=0 &&
       (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
      if(Params::Parameters().DoFreq() &&
	 //( Params::Parameters().DoCBS() || Params::Parameters().DoCCSDTCorrection() )){
	 ( Params::Parameters().GetQMType() == 2 )){
	for(int j=0;j<3*Dimers[i].GetNumberOfAtoms();j++){
	  if(Params::Parameters().DoEnergyFiniteDifferenceFreqs()){
	    for(int k=0;k<=j;k++){
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,k,true,false);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,k,false,false);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,k,true,false);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,k,false,false);
	      ijob++; 
	    }
	  }
	  else{
	    JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,false);
	    ijob++;
	    JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,false);
	    ijob++;
	  }	  
	  if(Params::Parameters().DoCCSDTCorrection()){
	    
	    for(int k=0;k<=j;k++){
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,k,true,true);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,k,false,true);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,k,true,true);
	      ijob++;
	      JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,k,false,true);
	      ijob++;
	      
	    }
	  }
	}
      }else{
	if (Params::Parameters().TinkerDebug()){
	  JobList[ijob] = Dimers[i].RunMMJob();
	  ijob++;
	}else{
	  //JobList[ijob] = Dimers[i].RunQChemJob();
	  JobList[ijob] = Dimers[i].RunQMJob();
	  ijob++;
	  if(Params::Parameters().DoCCSDTCorrection()){
	    JobList[ijob] = Dimers[i].RunQMJob(true);
	    ijob++;
	    
	    //counterpoise corrected CCSD(T) force calculations are numerical with each step in a different input file
	    if(Params::Parameters().DoForces() && Params::Parameters().DoCounterpoise() ){
	      for(int j=0;j<3*Dimers[i].GetNumberOfAtoms();j++){
		JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,true,true);
		ijob++;
		JobList[ijob] = Dimers[i].RunFiniteDifferenceMolProJob(j,false,true);
		ijob++;
		
	      }
	    }
	  }
	}
      }
    }
  }
  


  // If doing periodic boundary conditions, have additional dimer jobs
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      double separation =  DimerImages[i].GetDimerSeparation();
      
      if(DimerImages[i].GetSymmetryFactor()!=0 &&
	 (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	if(Params::Parameters().DoFreq() &&
	 ( Params::Parameters().GetQMType() == 2 )){
	   //( Params::Parameters().DoCBS() || Params::Parameters().DoCCSDTCorrection() )){
	  
	  for(int j=0;j<3*DimerImages[i].GetNumberOfAtoms();j++){
	    if(Params::Parameters().DoEnergyFiniteDifferenceFreqs()){
	      for(int k=0;k<=j;k++){
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,k,true,false);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,k,false,false);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,k,true,false);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,k,false,false);
		ijob++; 
	      }
	    }else{
	      JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,false);
	      ijob++;
	      JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,false);
	      ijob++;
	    }
	    if(Params::Parameters().DoCCSDTCorrection()){
	      
	      for(int k=0;k<=j;k++){
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,k,true,true);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,k,false,true);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,k,true,true);
		ijob++;
		JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,k,false,true);
		ijob++;
		
	      }
	    }
	  }
	  
	}else{
	  
	  if (Params::Parameters().TinkerDebug()){
	    JobList[ijob] = DimerImages[i].RunMMJob();
	    ijob++;
	  }else{
	    //JobList[ijob] = DimerImages[i].RunQChemJob();
	    JobList[ijob] = DimerImages[i].RunQMJob();
	    ijob++;
	    if(Params::Parameters().DoCCSDTCorrection()){
	      JobList[ijob] = DimerImages[i].RunQMJob(true);
	      ijob++;  
	      
	      //counterpoise corrected CCSD(T) force calculations are numerical with each step in a different input file
	      if(Params::Parameters().DoForces() && Params::Parameters().DoCounterpoise() ){
		for(int j=0;j<3*DimerImages[i].GetNumberOfAtoms();j++){
		  JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,true,true);
		  ijob++;
		  JobList[ijob] = DimerImages[i].RunFiniteDifferenceMolProJob(j,false,true);
		  ijob++;
		}
	      }
	      
	    }
	  }
	}
      }
    }
  }

  if (Params::Parameters().IsPeriodic() ){
    if(Params::Parameters().DoCCSDTCorrection())
      printf("Need to run %d QM monomer jobs, %d dimer jobs, %d dimer image jobs.\n",2*QM_Mon,2*UniqueDim,2*NDim_images);
    else
      printf("Need to run %d QM monomer jobs, %d dimer jobs, %d dimer image jobs.\n",QM_Mon,UniqueDim,NDim_images);
  }      
  else
    if(Params::Parameters().DoCCSDTCorrection())
      printf("Need to run %d QM monomer jobs and %d dimer jobs.\n",2*QM_Mon,2*UniqueDim);
    else
      printf("Need to run %d QM monomer jobs and %d dimer jobs.\n",QM_Mon,UniqueDim);

  // If benchmarking with full QM job
  if ( Params::Parameters().DoQMBenchmark() ) {
    printf("Also running full cluster QM job as a benchmark\n");
    JobList[ijob] = RunFullClusterQChemJob(false);
    ijob++;
  }
 }

 return JobList;
}

// Returns a list of commands for running each job
string* Cluster::RunMMJobs() {

  //flag for using symmetry for the MM; AIFF only
  bool UseMMSym = Params::Parameters().UseMMSymmetry();

  int Njobs = UniqueMon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc + 1;
  
  //if monomer symmetry is not being exploited
  if(!Params::Parameters().UseMonomerSymmetry())
    Njobs = NMon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc + 1;

  //only run Full MM
  if(Params::Parameters().UseFullMMOnly() && Params::Parameters().GetMMType() == 1)
    Njobs = 1;

  //full PBC tinker Hessian will be found by finite difference later on
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
      if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs()  
	   && Params::Parameters().Do_fdTinkerHessian()){
	   // && Params::Parameters().IsPeriodic() ){

	//Njobs += 6*GetTotalNumberOfAtoms() - 1; 
	printf("Removing job DoFreq = %i\n",Params::Parameters().DoFreq());
	Njobs += -1;
      }
    }

   if ( Params::Parameters().GetMMType() == 2 ) {// AIFF - monomers only for now
    if(UseMMSym)
      Njobs = UniqueMon;
    else
      Njobs = NMon;  
  } 
  

  //Crystal
  if(Params::Parameters().GetMMType() == 5){
    Njobs = 1;
    //Finite Difference for force calculations
    if(Params::Parameters().DoForces()){
 
      //Finite Difference for Nuclear coordinates  
      //Only for atoms in the asymetrical unit and coordinates
      //not fixed by symmetry
      if(!Params::Parameters().FreezeNuclearCoordinates()){
	for(int iMon=0;iMon<NMon;iMon++){
	  for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
	    if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){

	      if(!Monomers[iMon].GetAtom(iAtom).IsAtomFrozen()){
		if(Monomers[iMon].GetAtom(iAtom).IsYLocked() &&
		   Monomers[iMon].GetAtom(iAtom).IsZLocked() != 0)
		  Njobs += 2;
		else if(Monomers[iMon].GetAtom(iAtom).IsYLocked() ||
			Monomers[iMon].GetAtom(iAtom).IsZLocked() != 0)
		  Njobs += 4;
		else
		  Njobs += 6;
	      }   
	    } 
	  }
	}
      }
      //Finite Difference for Lattice Parameters
      if(!Params::Parameters().FreezeUnitCellParams()){
	//Finite Diff for length A
	Njobs += 2;
	if(!IsBLocked())//Finite Diff for length B
	  Njobs += 2;
	if(!IsCLocked()) //Finite Diff for length C
	  Njobs += 2;
	if(!IsAlphaLocked())//Finite Diff for angle alpha
	  Njobs += 2;
	if(!IsBetaLocked())//Finite Diff for angle beta
	  Njobs += 2;
	if(!IsGammaLocked())//Finite Diff for angle gamma
	  Njobs += 2;
      }
    }
  }

  string *JobList = new string[Njobs];

  int ijob = 0;


  //printf("Njobs = %i\n",Njobs);

  if(!Params::Parameters().UseFullMMOnly() /*&& Params::Parameters().GetMMType() == 1*/){
    // Monomers
    for (int i=1;i<=NMon;i++) {
      int sym_fac = Monomers[i].GetSymmetryFactor();
      if(Params::Parameters().GetMMType() == 2){//AIFF
	//run if symmetry factor is not zero, if the MM is exploiting symmetry
	if((sym_fac && UseMMSym) || !UseMMSym){
	  JobList[ijob] = Monomers[i].RunMMJob();
	  ijob++;
	}
      }else{//all other job types
	//run if symmetry factor is not zero, if the monomer is exploiting symmetry
	if(sym_fac || !Params::Parameters().UseMonomerSymmetry()){
	  JobList[ijob] = Monomers[i].RunMMJob();
	  ijob++;
	}
      }
    }
    // If we are not using the AIFF, we also evaluate dimers and full
    // cluster here.  If we are using the AIFF, this happens later, when
    // we read the energies.  For the AIFF, only the monomer properties
    // are determined here.
    if ( Params::Parameters().GetMMType() != 2 ) {// not the AIFF
      // Dimers
      double c0 = Params::Parameters().GetLocalCutoff(0);
      double c1 = Params::Parameters().GetLocalCutoff(1);    
      
      for (int i=1;i<=NDim;i++) {
	if(Dimers[i].GetSymmetryFactor()){// run if symmetry factor does not equal zero
	  double separation =  Dimers[i].GetDimerSeparation();
	  if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	    JobList[ijob] = Dimers[i].RunMMJob();
	    ijob++;
	  }
	}
      }
      
      // If doing periodic boundary conditions, have additional dimer jobs
      if ( Params::Parameters().IsPeriodic() ) {
	printf("Including MM Dimer Image jobs\n");
	for (int i=1;i<=NDim_images;i++) {
	  double separation =  DimerImages[i].GetDimerSeparation();
	  if(!Params::Parameters().DoLocal2BodyTruncation() ||separation <= c0){
	    JobList[ijob] = DimerImages[i].RunMMJob();
	    ijob++;
	  }
	}
	
      }
    }
  }

  // Full cluster
  if ( Params::Parameters().GetMMType() == 1 ) {// Tinker
    //if(!Params::Parameters().DoFreq() || !Params::Parameters().IsPeriodic()){
      if(!Params::Parameters().DoFreq() || !Params::Parameters().Do_fdTinkerHessian()){
      JobList[ijob] = RunFullClusterTinkerJob();
      ijob++;
    }
  }
  else if ( Params::Parameters().GetMMType() == 3 ) {// QChem
    JobList[ijob] = RunFullClusterQChemJob(true);
    ijob++;
  }
  else if( Params::Parameters().GetMMType() == 5 ) {// Crystal
    JobList[ijob] = RunCrystalJob();
    ijob++;
    if(Params::Parameters().DoForces()){

       int i=0;
       if(!Params::Parameters().FreezeNuclearCoordinates()){
         for(int iMon=0;iMon<NMon;iMon++){
	   for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
	     if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	       for(int xyz=0;xyz<3;xyz++){
		 if( !Monomers[iMon].GetAtom(iAtom).IsAtomFrozen() && //step atom that is frozen by symmetry
		    (xyz!=1 ||!Monomers[iMon].GetAtom(iAtom).IsYLocked()) && //step y coord if locked by symmetry
		    (xyz!=2 ||!Monomers[iMon].GetAtom(iAtom).IsZLocked()) ){ //step z coord if locked by symmetry
		  //plus jobs
		  JobList[ijob] = RunFiniteDifferenceCrystal(i,1);
		  ijob++;
		  //minus jobs
		  JobList[ijob] = RunFiniteDifferenceCrystal(i,0);
		  ijob++;
	        }
	        i++;
	      }
	    }
	  } 
        }
      }  

      if(!Params::Parameters().FreezeUnitCellParams()){
	int i = SymmetryUniqueCoordinates.GetLength();
        for(int j=0;j<6;j++){
          if((j!=1 || !IsBLocked()) && //Skip Length B if it's locked by symmetry
	     (j!=2 || !IsCLocked())  && //Skip Length C if it's locked by symmetry
	     (j!=3 || !IsAlphaLocked()) && //Skip Angle Alpha if it's locked by symmetry
	     (j!=4 || !IsBetaLocked()) && //Skip Angle Beta if it's locked by symmetry
	     (j!=5 || !IsGammaLocked()) ){ //Skip Angle Gamma if it's locked by symmetry
            JobList[ijob] = RunFiniteDifferenceCrystal(i+j,1); 
            ijob++;
            JobList[ijob] = RunFiniteDifferenceCrystal(i+j,0);  
            ijob++;
          }         
        }      
      }
    }
  }
  return JobList;
  
}

// ---------------------------------------------------------------
// begin    by Ali
//
// Returns command to run the Orient job
string Cluster::RunOrientJob() {
 
   string cmd;

   CreateOrientJob();

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  printf("%s \n", job_path.c_str());
  cmd += "pwd ;";
  cmd += "orient < orient.in > orient.out ";
  printf("%s \n", cmd.c_str());

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  cmd += "; ";

  printf("Here is RunOrientJob \n");

  return cmd;

}
// end by Ali
// ---------------------------------------------------------------

// Returns command string for running the qchem job
string Cluster::RunFullClusterQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();


  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  string infile = path + "/full.in";
  string outfile = path + "/full.out"; 

  // Execute Q-Chem
  // Q-Chem wants to be in the local directory, so we have to work in
  // the proper directory and with local filenames

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_outfile = "full.out";  
  string local_infile = "full.in";
 
  cmd += "qchem " + local_infile;
  cmd += " ";
  cmd += local_outfile;
  cmd += "; ";

  /*
  // Rename force file, if appropriate
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    string force_file = "full.force";
    cmd += "mv -f force.dat " + force_file;
    cmd += ";";
  }
  */

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// Returns command to run the full cluster Quantum Espresso job
string Cluster::RunQuantumEspressoJob() { //JLM
  // The job:
  string infile = "fullQE.in";
  string outfile = "fullQE.out";
  string infile2 = "phQE.in";
  string outfile2 = "phQE.out";  

  // First command, change to local directory
  string job_path = Params::Parameters().GetQMPath(); 
  if(Params::Parameters().DoFreq())
    job_path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += Quasiharmonic::quasiharmonic().GetHessianType();

  string cmd = "cd " + job_path;
  cmd += "; ";

  int Nproc = Params::Parameters().GetNumberOfProcessors();

  char NProc[10];
  sprintf(NProc,"%d",Nproc);
  // Second command, run the job
  if (Params::Parameters().DoFreq()){
    if(!Params::Parameters().UsePhonopyFreq()){
      //cmd += "mpirun -np " + std::string(NProc) + " ph.x " + " < " + infile2 + " > " + outfile2 + "; ";
      cmd += "bash run.sh "+ std::string(NProc) + "; ";
    }
    else{
    string dimension = Params::Parameters().GetQuantumEspressoSupercellSection().c_str();
    cmd += "bash run.sh fullQE "+ std::string(NProc) +" \"" + dimension + "\"; ";
    }
  
  }
  else {//if (Params::Parameters().DoForces()){
    cmd += "mpirun -np " + std::string(NProc) + " pw.x " + " -i " + infile + " > " + outfile + "; ";
  }

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("cmd = %s\n",cmd.c_str());
  fflush(stdout);
  //exit(0);
  return cmd;
}

// Returns command to run the full cluster DFTB job
string Cluster::RunDFTBJob() { //JLM
  // The job:
  string infile = "dftb_in.hsd";
  string outfile = "dftb_out.hsd";
  string modes_infile = "modes_in.hsd";
  string modes_outfile = "modes_out.hsd";
  //string infile2 = "phQE.in";
  //string outfile2 = "phQE.out";  

  // First command, change to local directory
  string job_path = Params::Parameters().GetQMPath(); 
  if(Params::Parameters().DoFreq())
    job_path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += Quasiharmonic::quasiharmonic().GetHessianType();

  string cmd = "cd " + job_path;
  cmd += "; ";

  if (Params::Parameters().DoFreq()){
    cmd += " mkdir gamma ; ";
    cmd += " mv *gamma*.hsd gamma ; ";
    cmd += " cd gamma ; ";
    cmd += " mv dftb_gamma_in.hsd " + infile + "; ";
    cmd += " mv modes_gamma_in.hsd " + modes_infile + "; ";
    cmd += " dftb+ " + infile + " > " + outfile + "; ";
    cmd += " modes " + modes_infile + " > " + modes_outfile + "; ";
    cmd += " mv hessian.out ../hessian_gamma.out; ";
    cmd += " cd ../ ; ";
    if(Supercell::supercell().GetSupercellInit()){  
      cmd += " mkdir phn ; ";
      cmd += " mv *phn*.hsd phn ; ";
      cmd += " cd phn ; ";
      cmd += " mv dftb_phn_in.hsd " + infile + "; ";
      cmd += " mv modes_phn_in.hsd " + modes_infile + "; ";
      cmd += " dftb+ " + infile + " > " + outfile + "; ";
      cmd += " modes " + modes_infile + " > " + modes_outfile + "; ";
      cmd += " mv hessian.out ../hessian_phn.out; ";
      cmd += " cd ../ ; ";
    }
  }
  else{
    cmd += " mv dftb_gamma_in.hsd " + infile + "; ";
    cmd += " dftb+ " + infile + " > " + outfile + "; ";
  }
  //Unfortunately DFTB+ is currently only compiled in serial. Compile in parallel later
  //int Nproc = Params::Parameters().GetNumberOfProcessors();

  //char NProc[10];
  //sprintf(NProc,"%d",Nproc);
  // Second command, run the job
  /*if (Params::Parameters().DoFreq()){
    if(!Params::Parameters().UsePhonopyFreq()){
      //cmd += "mpirun -np " + std::string(NProc) + " ph.x " + " < " + infile2 + " > " + outfile2 + "; ";
      cmd += "bash run.sh "+ std::string(NProc) + "; ";
    }
    else{
    string dimension = Params::Parameters().GetQuantumEspressoSupercellSection().c_str();
    cmd += "bash run.sh fullQE "+ std::string(NProc) +" \"" + dimension + "\"; ";
    }
  
  }
  else {//if (Params::Parameters().DoForces()){*/

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  printf("cmd = %s\n",cmd.c_str());
  //fflush(stdout);
  //exit(0);
  return cmd;
}


// Returns command to run the full tinker job
string Cluster::RunFullClusterTinkerJob() {

  // The job:
  string infile = "full.xyz";
  string outfile = "full.out";

  // To Execute Tinker:
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  if(Params::Parameters().GetJobTypeStr() == "hessian" 
     && Params::Parameters().DoFiniteDifferenceFreqs())
    job_path = Params::Parameters().GetHessianMMPath();


  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "analyze " + infile;
  cmd += " e > full.out; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  
  /* Actual running of jobs is now handled by main.C to simplify parallel code
  printf("Running Tinker energy calculation on %s\n",infile.c_str());
  fflush(stdout); // flush the output stream
  //printf("Executing: %s\n",cmd.c_str());
  system(cmd.c_str());
  */

    // If doing force job, compute the gradient
  if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    // Change to the local directory
    string cmd2 = "cd " + job_path;
    cmd2 += "; ";

    // Run the job
    cmd2 += "minimize full.xyz 100000000 > full.force;";
    
    // Rename the force file
    //cmd2 += "mv -f force.dat full.force;";

    // Remove tmp file and extra geom file created by minimize job
    cmd2 += "rm -f full.xyz_2; ";

    // Switch back to base directory
    cmd2 += "cd " + Params::Parameters().GetBasePath();

    // Actual running of jobs is handled by main.C to simplify parallel code
    //printf("Executing: %s\n",cmd2.c_str());
    //system(cmd2.c_str());

    // Combine the two commands
    cmd += ";" + cmd2;
  }

  if (Params::Parameters().GetJobTypeStr() == "hessian" &&
      !Params::Parameters().DoFiniteDifferenceFreqs() ) {

    // Change to the local directory 
    string job_path_hess = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path_hess += "/" + Quasiharmonic::quasiharmonic().GetHessianType();    

    string cmd2 = "cd " + job_path_hess;
    cmd2 += "; ";                 
    
    // Run the job
    cmd2 += "vibrate full.xyz 1 > full.freq;";           

    // Remove tmp file and extra geom file created by vibrate job
    cmd2 += "rm -f full.001; ";

    // Switch back to base directory
    cmd2 += "cd " + Params::Parameters().GetBasePath();

    // Actual running of jobs is handled by main.C to simplify parallel code
    //printf("Executing: %s\n",cmd2.c_str());
    //system(cmd2.c_str());

    // Combine the two commands
    cmd += ";" + cmd2;
  }

  return cmd;

}

string Cluster::RunFiniteDifferenceTinkerJob(int CoordEntry,bool IsPlus){

  //creating string from CoordEntry for the filenames
  char number[10];
  sprintf(number,"%d",CoordEntry);
  
  string PlusOrMinus;
  
  if(IsPlus) PlusOrMinus = "+";
  else PlusOrMinus = "-";
  
  string infile; 
  infile  = "full" + string(number) + string(PlusOrMinus) + ".xyz";
  string outfile;
  outfile  = "full" + string(number) + string(PlusOrMinus);
  if(Params::Parameters().DoForces())
    outfile += ".out";
  else
    outfile += ".force";
  
  // To Execute Tinker:
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.
  string cmd;
  //force by finite difference, untested
  if(Params::Parameters().DoForces()|| Params::Parameters().DoFiniteDifferenceFreqs() ) {
    
    // Change to the local directory
    string job_path = Params::Parameters().GetMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
        job_path +=  "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    job_path += "/full/"; 

    cmd += "cd " + job_path;
    cmd += "; ";

    // Second command, run the job
    cmd += "analyze " + infile;
    cmd += " e >" + outfile + "; ";
    
    // Third command, switch back to base directory
    cmd += "cd " + Params::Parameters().GetBasePath();
    
  }
  
  // If doing force job to find hessian by finite difference
  if (Params::Parameters().DoFreq()) {
	
    // Change to the local directory
    string job_path = Params::Parameters().GetHessianMMPath(); 

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable() ) 
        job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    job_path += "/full/";
    cmd = "cd " + job_path;
    cmd += "; ";

    // Run the job
    cmd += "minimize " + infile + " 100000000 > " + outfile + ";";
    
    // Rename the force file
    //cmd2 += "mv -f force.dat full.force;";

    // Remove tmp file and extra geom file created by minimize job
    cmd += "rm -f " + infile + "_2; ";
    //cmd += "rm -f full.xyz_2; ";

    // Switch back to base directory
    cmd += "cd " + Params::Parameters().GetBasePath();

    // Actual running of jobs is handled by main.C to simplify parallel code
    //printf("Executing: %s\n",cmd.c_str());
    //system(cmd.c_str());
  }

  return cmd;
    

}

string Cluster::RunCrystalJob() {
  
 // The job:
  string name = "full";

  // To Execute Tinker:
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  if(Params::Parameters().DoFreq() 
     && !Params::Parameters().DoFiniteDifferenceFreqs())
    job_path = Params::Parameters().GetHessianMMPath();



  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //First command, change directories
  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "runcry09 " + name;
  cmd += ";";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
    

}

string Cluster::RunFiniteDifferenceCrystal(int number,bool Plus){

  int NCoord = SymmetryUniqueCoordinates.GetLength();

  string PlusOrMinus;
  if(Plus) PlusOrMinus = "+";
  else PlusOrMinus = "-";

  string name;
  //Nuclear Coordiates
  if(number< NCoord){
    char num[10];
    sprintf(num,"%d",number);
    name = "full" + PlusOrMinus + num; 
  }//Length a
  else if(number == NCoord){
    name = "full" + PlusOrMinus + "a";
  }//Length b
  else if(number == NCoord+1){
    name = "full" + PlusOrMinus + "b";
  }//Length c
  else if(number == NCoord+2){
    name = "full" + PlusOrMinus + "c";
  }//angle alpha
  else if(number == NCoord+3){
    name = "full" + PlusOrMinus + "alpha";
  }//angle beta
  else if (number == NCoord+4){
    name = "full" + PlusOrMinus + "beta";
  }//angle gamma
  else if(number == NCoord+5){
    name = "full" + PlusOrMinus + "gamma";
  }else{
    printf("ERROR::Cluster::RunFiniteDifferenceCrystal(): varable number too high number = %i\n",
            number);
    exit(0);
  }

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  if(Params::Parameters().GetJobTypeStr() == "hessian" 
     && Params::Parameters().DoFiniteDifferenceFreqs())
    job_path = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //First command, change directories
  string cmd = "cd " + job_path;
  cmd += "; ";  
  
  //Second command, run crystal
  cmd += "runcry09 " + name + ";";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  //printf("cmd = %s\n",cmd.c_str());
  return cmd;

}

void Cluster::ReadQMResults() {

 //If we're doing a full cluster only skip all monomers and dimers
 if (!Params::Parameters().UseFullQMOnly()){
   // Monomers
   for (int i=1;i<=NMon;i++) {
     int Sym_Mon = Monomers[i].GetSymmetricalMonomer();
     //if(Monomers[i].GetSymmetryFactor()!=0){
     if ( Params::Parameters().TinkerDebug()){
       Monomers[i].ReadMMResults(Monomers[Sym_Mon]);
       Monomers[i].SetQMEnergy( Monomers[i].GetMMEnergy());
       Monomers[i].SetQMGradient( Monomers[i].GetMMGradient());
     }
     else
       Monomers[i].ReadQMResults(Monomers[Sym_Mon]);
     //}
   }

   double c0 = Params::Parameters().GetLocalCutoff(0);
   double c1 = Params::Parameters().GetLocalCutoff(1);  

   // Dimers
   for (int i=1;i<=NDim;i++) {
     double separation =  Dimers[i].GetDimerSeparation();
     int m1 = Dimers[i].GetIndexA();
     int m2 = Dimers[i].GetIndexB();
     // Update dimers with modified monomer data
     //printf("m1 = %d, m2 = %d\n",m1,m2);

     if(Dimers[i].GetSymmetryFactor() != 0 && (!Params::Parameters().DoLocal2BodyTruncation() ||separation <c0)){
       Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);
       if ( Params::Parameters().TinkerDebug() ) 
	 Dimers[i].ReadMMResults();
       else{
        Dimers[i].ReadQMResults();
       }
       if ( Params::Parameters().TinkerDebug() && Params::Parameters().DoForces()  ) {
	 Dimers[i].SetQMGradient( Dimers[i].GetMMGradient() );
       }
     }
   }

   // Additional dimers for periodic system
   if ( Params::Parameters().IsPeriodic() ) {
     //printf("Reading Dimer_image QM results\n");
     for (int i=1;i<=NDim_images;i++) {
       double separation =  DimerImages[i].GetDimerSeparation();
       int m1 = DimerImages[i].GetIndexA();
       int m2 = DimerImages[i].GetIndexB();
       int ref_mon = DimerImages[i].GetReferenceMonomerIndex();
      
       // Since we don't store a separate list of all the image monomers, we
       // grab it directly from the dimer object.
       Monomer MonB = DimerImages[i].GetMonomerB();
      
       //fflush(stdout);
       DimerImages[i].UpdateObjects(Monomers[m1],MonB,Monomers[ref_mon]);
       //DimerImages[i].UpdateObjects(Monomers[m1],MonB);
       //if ( Params::Parameters().TinkerDebug() ) {
       //Dimers[i].ReadMMResults();
       //DimerImages[i].SetQMGradient( DimerImages[i].GetMMGradient() );
       //}
       //else
       //printf("ReadQMResults for the DimerImages\n"); fflush(stdout);
       //if(DimerImages[i].GetSymmetryFactor() != 0)
       if(separation <=c0)
	 DimerImages[i].ReadQMResults( Monomers[ref_mon] );
      
       m1 = DimerImages[i].GetIndexA();
       m2 = DimerImages[i].GetIndexB();
       ref_mon = DimerImages[i].GetReferenceMonomerIndex();
       //printf("verify: m1 = %d, m2 = %d, ref_m2 = %d\n",m1,m2,ref_mon);

       //debug printing - GJB
       //DimerImages[i].PrintQMGradient("QM dimer gradient\n");
       //DimerImages[i].GetMonomerA().PrintQMGradient("QM Mon A gradient\n");
       //DimerImages[i].GetMonomerB().PrintQMGradient("QM Mon B gradient\n");
	
     }
     //printf("Done updating QM in PBC objects\n"); fflush(stdout);
   }
 }
 else { //Full Cluster; JLM
   SetQMEnergy();
   if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
     SetQMGradient();
     //PrintGradient("Full QM Gradient",Grad_QM);
   }
   if (Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() ) {
     SetQMHessian();
     //PrintHessian("Full QM Hessian",Hess_QM);
   }
 } 

  // If we did a benchmark full cluster QM calculation
  if ( Params::Parameters().DoQMBenchmark() ) {
    Energy_QM = ReadQChemEnergy(false);
  }

}

void Cluster::SetQMEnergy() {
  double energy;
  if (Params::Parameters().GetQMType()==5) // Quantum Espresso
    energy = ReadQuantumEspressoEnergy();
  else if (Params::Parameters().GetQMType()==8) // DFTB+
    energy = ReadDFTBEnergy();
  else {
    printf("Cluster::SetQMEnergy: Unknown QM_type: %d\n",
	   Params::Parameters().GetQMType() );
    exit(1);
  }
  Energy_QM = energy;
}

//Read Quantum Espresso Energy
double Cluster::ReadQuantumEspressoEnergy(){

  double energy = 0.0;
  //Path of QM
  string path = Params::Parameters().GetQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/fullQE"; 

  string out_filename = filename + ".out";

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadQuantumEspressoEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;

  infile.seekg(0,ios::beg);
  ifstream::pos_type posBeg = infile.tellg();
  infile.seekg(-1, ios::end);
  char buf;

  // Read in the energy from the bottom of the file
  while ( infile.tellg() != posBeg ) {
    buf = static_cast <char>(infile.peek());
    if (buf != '\n') {
      line += buf;
    }
    else {
      reverse(line.begin(), line.end());
      // Do something interesting with the line.
      string match = line.substr(0,17);
      // Search for final Quantum Espresso energy
      if ( match=="!    total energy" ) {
        istringstream iss(line);
        string tmp;
        for (int i=0;i<4;i++)
          iss >> tmp; // throw away text
        iss >> energy; // read energy
        break;
      }
      line.clear();
    }

    infile.seekg(-1, ios::cur); 
  }

  //Note that Quantum Espresso energy is in units of Rydberg. Need to convert to Hartrees.
  energy /= 2;
  if ( Params::Parameters().PrintLevel() > 0) printf("Quantum Espresso obtained QM Energy = %15.9f\n",energy);

  // Close the output file
  infile.close();

  return energy;
}

//Read DFTB Energy
double Cluster::ReadDFTBEnergy(){

  double energy = 0.0;
  //Path of QM
  string path = Params::Parameters().GetQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/dftb_out.hsd"; 

  //string out_filename = filename + ".out";

  // Open the energy file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadDFTBEnergy : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;

  infile.seekg(0,ios::beg);
  ifstream::pos_type posBeg = infile.tellg();
  infile.seekg(-1, ios::end);
  char buf;

  // Read in the energy from the bottom of the file
  while ( infile.tellg() != posBeg ) {
    buf = static_cast <char>(infile.peek());
    if (buf != '\n') {
      line += buf;
    }
    else {
      reverse(line.begin(), line.end());
      // Do something interesting with the line.
      string match = line.substr(0,13);
      // Search for final Quantum Espresso energy
      if ( match=="Total Energy:" ) {
        istringstream iss(line);
        string tmp;
        iss >> tmp; // throw away text
        iss >> tmp; // throw away text
        iss >> energy; // read energy
        break;
      }
      line.clear();
    }

    infile.seekg(-1, ios::cur); 
  }

  if ( Params::Parameters().PrintLevel() > 0) printf("DFTB obtained QM Energy = %15.9f\n",energy);

  // Close the output file
  infile.close();
  //Energy is in units of Hartrees
  return energy;
}

void Cluster::ReadMMResults() {


  //skip monomer and dimer if calculating only full MM
  if(!Params::Parameters().UseFullMMOnly())
    // Monomers
    for (int i=1;i<=NMon;i++) {
      int sym_Mon = Monomers[i].GetSymmetricalMonomer();
      Monomers[i].ReadMMResults(Monomers[sym_Mon]);
    }
    
    double c0 = Params::Parameters().GetLocalCutoff(0);
    double c1 = Params::Parameters().GetLocalCutoff(1);   

    if(!Params::Parameters().UseFullMMOnly()){
      // Dimers
      for (int i=1;i<=NDim;i++) {
	double separation =  Dimers[i].GetDimerSeparation();
	//do not get energy if the symmetry factor = 0 unless the MM type is AIFF
	if(Dimers[i].GetSymmetryFactor() || Params::Parameters().GetMMType() == 2){
	  int m1 = Dimers[i].GetIndexA();
	  int m2 = Dimers[i].GetIndexB();
	  
	  // Update dimers with modified monomer data
	  Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);
	  if(!Params::Parameters().DoLocal2BodyTruncation() || separation <= c0)
	    Dimers[i].ReadMMResults();
	}
      }
      
      // Additional dimers for periodic system
      if ( Params::Parameters().IsPeriodic() ) {
	/*
	  if ( Params::Parameters().GetMMType()==2) {
	  printf("Warning: Trying to do PBC AIFF electrostatics.  Not implemented\n");
	  //exit(1);
	  }
	*/
	for (int i=1;i<=NDim_images;i++) {
	  double separation =  DimerImages[i].GetDimerSeparation();
	  int m1 = DimerImages[i].GetIndexA();
	  int m2 = DimerImages[i].GetIndexB();
	  int ref_mon = DimerImages[i].GetReferenceMonomerIndex();
	  
	  Monomer MonB = DimerImages[i].GetMonomerB();
	  DimerImages[i].UpdateObjects(Monomers[m1],MonB,Monomers[ref_mon]);
	  if (separation <= c0)
	    DimerImages[i].ReadMMResults( Monomers[ref_mon] );
	}
	printf("Done updating MM in PBC objects\n"); fflush(stdout);
      }
    }
    // Full Cluster
    // Read in the energy of the full cluster at the MM level
    if ( Params::Parameters().GetMMType()==2 && Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) { 
      SetMMGradient(); // get energies & gradient simultaneously in this case
    //PrintGradient("Full MM Gradient",Grad_MM);
  }
  if(Params::Parameters().GetMMType()==2 && Params::Parameters().DoFreq() &&
     !Params::Parameters().DoFiniteDifferenceFreqs()){
    SetMMHessian(); //get hessian
    //PrintHessian("Full MM Hessian",Hess_MM)
  }
  else {
    SetMMEnergy();
    if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
      SetMMGradient();
     //PrintGradient("Full MM Gradient",Grad_MM);
    }
    if (Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() ) {
      SetMMHessian();
     //PrintHessian("Full MM Hessian",Hess_MM);
    }

  }
}

void Cluster::SetMMEnergy() {
  double energy;
  if (Params::Parameters().GetMMType()==1) // Tinker
    energy = ReadTinkerEnergy();
  else if (Params::Parameters().GetMMType()==2) //Orient
    energy = ComputeAIFFEnergy();
  else if (Params::Parameters().GetMMType()==3) // QChem
    energy = ReadQChemEnergy(true);
  else if (Params::Parameters().GetMMType()==5) //Crystal
    energy = ReadCrystalEnergy();
  else {
    printf("Cluster::SetMMEnergy: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  }

  Energy_MM = energy;
  //printf("MM Energy = %15.9f\n",Energy_MM);

}

/*
//added
void Cluster::GetFailedJobs(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  if (MM_job)
    path = Params::Parameters().GetMMPath();
  else
    path = Params::Parameters().GetQMPath();
  string out_filename = path + "/full.out";

  // Open the energy file       
  ifstream infile;
  infile.open( out_filename.c_str() ); 
  if ( !infile.is_open() ) {
    

//    printf("Cluster::ReadQChemEnergy : Cannot open file '%s'\n", out_filename.c_str());             





//added
*/




double Cluster::ReadQChemEnergy(bool MM_job) {
  double energy;

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  string out_filename = path + "/full.out"; 
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadQChemEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }
  
  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++)
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }
  }

  if ( Params::Parameters().PrintLevel() > 0) printf("QM Energy = %15.9f\n",energy);
  
  // Close the force file
  infile.close();
  return energy;
}

double Cluster::ReadTinkerEnergy(string filename) {

  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  double energy = 0.0;

  // Set up the filename with the full path.  Default file is "full.out".
  string out_filename = path + "/" + filename; 

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadTinkerEnergy : Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    /*
    string match = line.substr(0,13);
    if ( match==" Polarization" ) {
      istringstream iss(line);
      string tmp; 
      iss >> tmp; // throw away "Polarization" tag
      if (Params::Parameters().IsPeriodic() ) iss >> tmp; // get rid of "Ewald text"
      iss >> energy; // read the energy
    }
    */

    string match = line.substr(0,23);
    if ( match==" Total Potential Energy" ) {
      istringstream iss(line);
      string tmp; 
      // throw away "text tags"
      iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; 
      iss >> energy; // read the energy
    }


  }
  if ( Params::Parameters().PrintLevel() > 0 ) 
    printf("Tinker Energy = %f kcal/mol\n",energy);
    
  energy /= HartreesToKcalpermole; // convert to hartrees
 
  // Close the force file
  infile.close();
  return energy;
}

double Cluster::ReadCrystalEnergy(){ // by yoni

  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  double energy = 0.0;

  // Set up the filename with the full path.  Default file is "full.out".
  string out_filename = path + "/full.out"; 

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadCrystalEnergy : Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

 // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,20)==" TOTAL ENERGY + DISP"){
      istringstream iss(line);
      string tmp; 
      // throw away "text tags"
      iss >> tmp;
      iss >> tmp;
      iss >> tmp;
      iss>> tmp;
      iss >> tmp;
      iss >> energy; // read the energy
    }
    if(line.substr(0,18)==" TOTAL ENERGY(DFT)") {
      istringstream iss(line);
      string tmp; 
      // throw away "text tags"
      iss >> tmp; 
      iss >> tmp;
      iss >> tmp;
      iss >> energy; // read the energy
    }

  }

  //printf("Crystal energy = %f hartrees\n",energy);
  infile.close();
  return energy;
}

double Cluster::ReadOrientEnergy() {   // by Ali

  string path = Params::Parameters().GetMMPath();

  double energy = 0.0;
  double dim_energy = 0.0;
  double sum_dim_energy = 0.0;

  // Set up the filename 'full.out' with the full path
  string out_filename = path + "/orient.out"; 

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadAIFFEnergy : Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;
  int i=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    /*
    string match = line.substr(0,13);
    if ( match==" Polarization" ) {
      istringstream iss(line);
      string tmp; 
      iss >> tmp; // throw away "Polarization" tag
      if (Params::Parameters().IsPeriodic() ) iss >> tmp; // get rid of "Ewald text"
      iss >> energy; // read the energy
    }
    */
  

    string match = line.substr(0,23);
    
    if ( match=="Total                  " ) {
	 istringstream iss(line);
         string tmp; 
         // throw away "text tags"
         iss >> tmp; // iss >> tmp; iss >> tmp; iss>> tmp; 
	 if (i==0){
             iss >> energy; // read the energy
             printf("Cluster::OrientFullEnergy :  '%f'\n",energy);}
	 else{	 
             iss >> dim_energy; // read the energy
             Dimers[i].SetOrientEnergy(dim_energy);
	     sum_dim_energy += dim_energy;
 //            printf("Cluster::OrientDimerEnergy :  '%f'\n",dim_energy);
             printf("Cluster::OrientDimer-%d Energy :  %f \n",i,Dimers[i].GetOrientEnergy());}
         i++;
     }
    
  } // end while
  
  energy = energy - sum_dim_energy;
  
  if ( Params::Parameters().PrintLevel() > 0 ) 
    printf("Orient Energy = %f kJ/mol\n",energy);
  
  printf("Cluster::OrientManyBodyEnergy (kJ/mol):  '%f'\n",energy);
  
  energy /= HartreesToKJpermole; // convert to hartrees
   
  // Close the force file
  printf("Cluster::OrientManyBodyEnergy (hartrees):  '%f'\n",energy);
  infile.close();
  return energy;
}

// Get the QM Gradient, wrapper routine
void Cluster::SetQMGradient() {
  if (Params::Parameters().GetQMType() == 5){  // Quantum Espresso 
    Grad_QM = ReadQuantumEspressoGradient();
    //Grad_QM.PrintGradient("Current QM Gradient");
    QM_Grad_Init = 1;
  }
  else if (Params::Parameters().GetQMType() == 8){  // DFTB 
    Grad_QM = ReadDFTBGradient();
    //Grad_QM.PrintGradient("Current QM Gradient");
    QM_Grad_Init = 1;
  }
  else {
    printf("Unrecognized QM gradient tag %i. Exiting...",Params::Parameters().GetQMType());
    exit(1);
  }
  
}

Vector Cluster::ReadQuantumEspressoGradient(){ //JLM

  int Natoms=GetTotalNumberOfAtoms();
  Vector grad(3*Natoms+9);
  Matrix stress_tensor(3,3);

  //Path of QM
  string path = Params::Parameters().GetQMPath();
  bool skip_stress = false;
  double volume_scalar = GetCellvolume()*AngToBohr*AngToBohr*AngToBohr;
  Matrix Volume(3,3);
  Volume.SetRowVector(GetUnitCellVector(0),0);
  Volume.SetRowVector(GetUnitCellVector(1),1);
  Volume.SetRowVector(GetUnitCellVector(2),2);
  Vector tmp_vec = GetUnitCellAxes();
  double alat = tmp_vec[0]*AngToBohr;
  //printf("alat = %f   volume_scalar = %f \n",alat,volume_scalar);
  Volume.Scale(AngToBohr);
  //Volume.Print(" Normal volume ");  

  //Invert matrix
  Volume.Inverse();
  //Volume.Print(" Volume inverted before scaling ");

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/fullQE"; 
  string out_filename = filename + ".out"; 

  // Open the file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadQuantumEspressoGradient : Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

  // Read in the data from the bottom of the file
  string line;
  infile.seekg(0,ios::beg);
  ifstream::pos_type posBeg = infile.tellg();
  infile.seekg(-1, ios::end);
  char buf;

  // Read in the energy from the bottom of the file
  while ( infile.tellg() != posBeg ) {
    buf = static_cast <char>(infile.peek());
    if (buf != '\n') {
      line += buf;
    }
    else {
      reverse(line.begin(), line.end());
      //printf("%s\n",line.c_str());
      //Stress tensor always prints out last. So it will be the first one read
      if ( line.substr(0,38)=="          total   stress  (Ry/bohr**3)" && !skip_stress) {
        getline(infile,line); // throw away blank line
        getline(infile,line); // regrab if statement line
        for (int i=0;i<3;i++) {
          getline(infile,line); // Now we have our gradient
          istringstream iss(line);
          Vector tmp(3);
          for (int j=0;j<3;j++) {
            iss >> tmp[j]; 
            //Their Stress tensor is in units of Ry/bohr^3. We want Hartree/bohr (bohr=au)
            //This is a force which = -gradient. Correct by multiplying by -1
            tmp[j] /= -2.0;
          }
          stress_tensor.SetRowVector(tmp,i);
        }
        skip_stress = true;
              
      }

      // Start reading in the molecular gradient
      if ( line=="     Forces acting on atoms (cartesian axes, Ry/au):" || line=="     Forces acting on atoms (Ry/au):" ) {
        getline(infile,line); // regrab previous line
        getline(infile,line); // throw away header line
        getline(infile,line); // throw away blank line
        for (int i=0;i<Natoms;i++) {
          getline(infile,line); // Now we have our gradient
          istringstream iss(line);
          string tmp; 
          for (int k=0;k<6; k++)
            iss >> tmp; // throw away indexes

          for (int j=0;j<3;j++) {
            iss >> grad[3*i+j];//gradients will be rotated if need due to symmetry 
            //Gradients are in units of Ry/au. We want Hartree/au
            //This is a force which = -gradient. Correct by multiplying by -1
            grad[3*i+j] /= -2.0;
	    // Store the gradient elements
          }
        }
        break;
      }        
      line.clear();
    }
    infile.seekg(-1, ios::cur); 
  }

  // Close the output file
  infile.close(); 

  Matrix lattice_vector_grad = stress_tensor.Multiply(Volume);

  /*printf("\n");
  lattice_vector_grad.Print("My fcell minus omega");
  printf("\n");*/

  double press = Params::Parameters().GetExternalPressure();

  if(fabs(press) != 0.0){
    Volume.Scale(press*GigaPascalsToHartreesperBohrcube);
    lattice_vector_grad += Volume;
    Volume.Scale(1.0/(press*GigaPascalsToHartreesperBohrcube));
  }

  lattice_vector_grad.Scale(volume_scalar); //We should divide by the weight as well (if given) but ours is 1.0
  
  //lattice_vector_grad.Print("lattice vector gradient");
  //printf("\n");

  //Vibrational contribution to the gradient using QuasiHarmonic Approximation
  if((Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())){
    //ComputeQuasiHarmonicFrequencies();//Do QuasiHarmonic Approximation

    //Isotropic Gruneisen
    if(Params::Parameters().UseVolumeGruneisen()){
      Matrix VibGradMat;
      VibGradMat = ComputeQuasiHarmonicGradient();
      printf("\n");
      VibGradMat.Print("VibGradMat");
      printf("\n");

/*      dE_dx1 += VibGradMat.Element(0,0); //JLM add in once understood what this does
      dE_dx2 += VibGradMat.Element(1,0);
      dE_dx3 += VibGradMat.Element(2,0);
      dE_dy2 += VibGradMat.Element(1,1);
      dE_dy3 += VibGradMat.Element(2,1);
      dE_dz3 += VibGradMat.Element(2,2);*/
      lattice_vector_grad += VibGradMat;
      lattice_vector_grad.Print("lattice vector gradient after Quasiharmonic addition");
      printf("\n");

    }//Anisotropic Gruneisen
    else{ //JLM might not work with my current code... just needs adjusting
      printf("Warning!!! Anisotropic Grueneisen parameters do not currently work with the Quantum Espresso code.\nExiting HMBI\n");
      exit(1);
      Vector VibGrad;
      VibGrad = Quasiharmonic::quasiharmonic().ComputeAnisotropicGradient();

      printf("Helmholtz Gradient\n");
      printf("dF/da = %f\n",VibGrad[0]);
      printf("dF/db = %f\n",VibGrad[1]);
      printf("dF/dc = %f\n",VibGrad[2]);
      printf("dF/dalpha = %f\n",VibGrad[3]);
      printf("dF/dbeta = %f\n",VibGrad[4]);
      printf("dF/dgamma = %f\n\n",VibGrad[5]);

      //Adding Helmholtz Gradient to Full MM gradient since both are already expressed in terms of params 
      //rather than vectors 
      /*Grad_Lat_Full_MM[0] += VibGrad[0]; //JLM
      Grad_Lat_Full_MM[1] += VibGrad[1];
      Grad_Lat_Full_MM[2] += VibGrad[2];
      Grad_Lat_Full_MM[3] += VibGrad[3];
      Grad_Lat_Full_MM[4] += VibGrad[4];
      Grad_Lat_Full_MM[5] += VibGrad[5];*/
    }
  }

  lattice_vector_grad(0,1) = 0.0;
  lattice_vector_grad(0,2) = 0.0;
  lattice_vector_grad(1,2) = 0.0;

  //Account for symmetry in the lattice vectors
  /*if(Params::Parameters().UseLatticeSymmetry()){
    Vector previousAngles = GetLastAcceptedUnitCellAngles();
    //if gamma = 90 degrees: zero out [2,1]
    //if alpha also = 90 degrees: zero out [3,2]
    if(previousAngles[2]==90){
      lattice_vector_grad(1,0) = 0.0;
      if(previousAngles[0]==90){
        lattice_vector_grad(2,1) = 0.0;
      }
    }
    
    //if beta = 90 degrees: zero out [3,1]
    if(previousAngles[1]==90){
      lattice_vector_grad(2,0) = 0.0;
    }
  }*/

  /*grad.Print("total gradient");
  printf("\n");
  lattice_vector_grad.Print("final lattice vector gradient");
  printf("\n");*/

  Vector tmp(9);
  for (int i = 0; i<3; i++)
    for (int j = 0; j<3; j++)
      tmp[3*i+j] = lattice_vector_grad(i,j);

  for (int i=0; i<9; i++)
    grad[3*Natoms+i] = tmp[i];


  //grad.PrintGradient("Final Gradient");
  //fflush(stdout);
  return grad;
}


Vector Cluster::ReadDFTBGradient(){ //JLM

  printf("Correctly found DFTB+ gradient!\n");
  
  int Natoms=GetTotalNumberOfAtoms();
  Vector grad(3*Natoms+9);
  Matrix stress_tensor(3,3);
  Matrix lattice_derivs(3,3);

  //Path of QM
  string path = Params::Parameters().GetQMPath();
  bool skip_stress = false;
  bool skip_lattice_derivs = false;
  double volume_scalar = GetCellvolume()*AngToBohr*AngToBohr*AngToBohr;
  Matrix Volume(3,3);
  Volume.SetRowVector(GetUnitCellVector(0),0);
  Volume.SetRowVector(GetUnitCellVector(1),1);
  Volume.SetRowVector(GetUnitCellVector(2),2);
  Vector tmp_vec = GetUnitCellAxes();
  double alat = tmp_vec[0]*AngToBohr;
  //printf("alat = %f   volume_scalar = %f \n",alat,volume_scalar);
  Volume.Scale(AngToBohr);
  //Volume.Print(" Normal volume ");  

  //Invert matrix
  Volume.Inverse();
  //Volume.Print(" Volume inverted before scaling ");

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/detailed.out"; 

  // Open the file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadDFTBGradient : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  // Read in the data from the bottom of the file
  string line;
  infile.seekg(0,ios::beg);
  ifstream::pos_type posBeg = infile.tellg();
  infile.seekg(-1, ios::end);
  char buf;

  // Read in the energy from the bottom of the file
  while ( infile.tellg() != posBeg ) {
    buf = static_cast <char>(infile.peek());
    if (buf != '\n') {
      line += buf;
    }
    else {
      reverse(line.begin(), line.end());
      //printf("%s\n",line.c_str());

      //Total Lattice derivs always prints out last. So it will be the first one read
      if ( line=="Total lattice derivs" && !skip_lattice_derivs) {
        getline(infile,line); // throw away blank line
        getline(infile,line); // regrab if statement line
        for (int i=0;i<3;i++) {
          getline(infile,line); // Now we have our gradient
          istringstream iss(line);
          Vector tmp(3);
          for (int j=0;j<3;j++) {
            iss >> tmp[j]; 
            //Gradients are in units of Hartree/bohr
            //This is a force which = -gradient. Correct by multiplying by -1
            //tmp[j] *= -1.0;          
            //tmp[j] /= -2.0;
          }
          lattice_derivs.SetRowVector(tmp,i);
        }
        skip_lattice_derivs = true;
        //lattice_derivs.Print("Lattice Derivs");
        //printf("\n");
      }

      //if ( line.substr(0,20)=="Total stress tensor" && !skip_stress) {
      if ( line=="Total stress tensor" && !skip_stress) {
        getline(infile,line); // throw away blank line
        getline(infile,line); // regrab if statement line
        for (int i=0;i<3;i++) {
          getline(infile,line); // Now we have our gradient
          istringstream iss(line);
          Vector tmp(3);
          for (int j=0;j<3;j++) {
            iss >> tmp[j]; 
            //Gradients are in units of Hartree/bohr
            //This is a force which = -gradient. Correct by multiplying by -1
            tmp[j] *= -1.0;          
            //tmp[j] /= -2.0;
          }
          stress_tensor.SetRowVector(tmp,i);
        }
        skip_stress = true;
        //stress_tensor.Print("Stress Tensor");
        //printf("\n");
              
      }

      // Start reading in the molecular gradient
      if ( line=="Total Forces") {
        getline(infile,line); // regrab previous line
        getline(infile,line); // throw away header line
        //getline(infile,line); // throw away blank line
        for (int i=0;i<Natoms;i++) {
          getline(infile,line); // Now we have our gradient
          istringstream iss(line);
          //string tmp; 
          //for (int k=0;k<6; k++)
          //  iss >> tmp; // throw away indexes

          for (int j=0;j<3;j++) {
            iss >> grad[3*i+j]; 
            //Gradients are in units of Hartree/bohr
            //This is a force which = -gradient. Correct by multiplying by -1
            grad[3*i+j] *= -1.0;
            //grad[3*i+j] /= -2.0;
	    // Store the gradient elements
          }
        }
        break;
      }        
      line.clear();
    }
    infile.seekg(-1, ios::cur); 
  }

  // Close the output file
  infile.close(); 

  Matrix lattice_vector_grad = stress_tensor.Multiply(Volume);

  /*printf("\n");
  lattice_vector_grad.Print("My fcell minus omega");
  printf("\n");*/

  double press = Params::Parameters().GetExternalPressure();

  if(fabs(press) != 0.0){
    Volume.Scale(press*GigaPascalsToHartreesperBohrcube);
    lattice_vector_grad += Volume;
    Volume.Scale(1.0/(press*GigaPascalsToHartreesperBohrcube));
  }

  lattice_vector_grad.Scale(volume_scalar); //We should divide by the weight as well (if given) but ours is 1.0
  
  /*grad.Print("atomic gradient");
  printf("\n");

  lattice_vector_grad.Print("lattice vector gradient");
  printf("\n");

  lattice_derivs.Print("lattice vector derivatives (from DFTB)");
  printf("\n");*/

  //Vibrational contribution to the gradient using QuasiHarmonic Approximation
  if((Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())){
    //ComputeQuasiHarmonicFrequencies();//Do QuasiHarmonic Approximation

    //Isotropic Gruneisen
    if(Params::Parameters().UseVolumeGruneisen()){
      Matrix VibGradMat;
      VibGradMat = ComputeQuasiHarmonicGradient();
      printf("\n");
      VibGradMat.Print("VibGradMat");
      printf("\n");

/*      dE_dx1 += VibGradMat.Element(0,0); //JLM add in once understood what this does
      dE_dx2 += VibGradMat.Element(1,0);
      dE_dx3 += VibGradMat.Element(2,0);
      dE_dy2 += VibGradMat.Element(1,1);
      dE_dy3 += VibGradMat.Element(2,1);
      dE_dz3 += VibGradMat.Element(2,2);*/
      lattice_vector_grad += VibGradMat;
      lattice_vector_grad.Print("lattice vector gradient after Quasiharmonic addition");
      printf("\n");

    }//Anisotropic Gruneisen
    else{ //JLM might not work with my current code... just needs adjusting
      printf("Warning!!! Anisotropic Grueneisen parameters do not currently work with the Quantum Espresso code.\nExiting HMBI\n");
      exit(1);
      Vector VibGrad;
      VibGrad = Quasiharmonic::quasiharmonic().ComputeAnisotropicGradient();

      printf("Helmholtz Gradient\n");
      printf("dF/da = %f\n",VibGrad[0]);
      printf("dF/db = %f\n",VibGrad[1]);
      printf("dF/dc = %f\n",VibGrad[2]);
      printf("dF/dalpha = %f\n",VibGrad[3]);
      printf("dF/dbeta = %f\n",VibGrad[4]);
      printf("dF/dgamma = %f\n\n",VibGrad[5]);

      //Adding Helmholtz Gradient to Full MM gradient since both are already expressed in terms of params 
      //rather than vectors 
      /*Grad_Lat_Full_MM[0] += VibGrad[0]; //JLM
      Grad_Lat_Full_MM[1] += VibGrad[1];
      Grad_Lat_Full_MM[2] += VibGrad[2];
      Grad_Lat_Full_MM[3] += VibGrad[3];
      Grad_Lat_Full_MM[4] += VibGrad[4];
      Grad_Lat_Full_MM[5] += VibGrad[5];*/
    }
  }

  // FD hessians could be non-symmetric due to numerical noise, so symmetrize them
  /*for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) {
      lattice_vector_grad.Element(i,j) += lattice_vector_grad.Element(j,i); 
      lattice_vector_grad.Element(i,j) *= 0.5;
    }
  }

  for (int i=1;i<3;i++) {
    for (int j=0;j<i;j++) {
      lattice_vector_grad.Element(i,j) *= 0.0; 
      lattice_vector_grad.Element(i,j) += lattice_vector_grad.Element(j,i);
    }
  }

  lattice_vector_grad(0,1) = 0.0;
  lattice_vector_grad(0,2) = 0.0;
  lattice_vector_grad(1,2) = 0.0;*/

  /*grad.Print("total gradient");
  printf("\n");
  lattice_vector_grad.Print("final lattice vector gradient");
  printf("\n");*/

  Vector tmp(9);
  for (int i = 0; i<3; i++)
    for (int j = 0; j<3; j++)
      tmp[3*i+j] = lattice_vector_grad(i,j);

  for (int i=0; i<9; i++)
    grad[3*Natoms+i] = tmp[i];

  //grad.PrintGradient("Final Gradient");
  //fflush(stdout);
  //fflush(stdout);
  //exit(0);

  return grad;
}

// Get the MM Gradient, wrapper routine
void Cluster::SetMMGradient() {
  if (Params::Parameters().GetMMType() == 1)  // Tinker 
    Grad_MM = ReadGradient(2);
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    Grad_MM = ComputeClusterMultipoleGradient();
  }
  else if (Params::Parameters().GetMMType() == 3) { // QChem
    Grad_MM = ReadGradient(1);
  }
  else if (Params::Parameters().GetMMType() == 5){ // Crystal
    //Grad_MM = ReadCrystalGradient();
    Grad_MM = ReadFiniteDifferenceCrystalGradient();
  }
 
  MM_Grad_Init = 1;
}

// Main routine for reading in the MM gradient 
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
// type 1 = qchem style, type 2 = tinker style
Vector Cluster::ReadGradient(int type) {
  
  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  int Natoms = GetTotalNumberOfAtoms();
  Vector TempGrad(3*Natoms);//stores entrees from gradient before being rotated by symmetry
  Vector grad;

  if ( (Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2 ) && !Params::Parameters().DoFiniteDifferenceFreqs())
    grad.Initialize(3*UniqueAtoms+6);
  else
    grad.Initialize(3*UniqueAtoms);
  //KDN: conditionally (i.e. if periodic and jobtype=opt and not doing hessian calculations) change this to a vec of length 3n+6 to accomodate components arising due to the 6 lattice constants

  if ( Params::Parameters().NeglectManyBody() )
    grad.Set();
  else {

    string filename;
    // Set up the filename, with the full path.  File is 'full.force'
    if (type == 2) // Tinker MM job
      filename = path + "/full.force";
    else
      filename = path + "/full.out";
    
    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadGradient : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
    
    // Read in the data - search for the "Nuclear forces:" string
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      // Search for final q-chem energy
      if ( line==" Nuclear forces:" ) {
	getline(infile,line); // throw away header line
	
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  for (int j=0;j<3;j++) {
	    iss >> TempGrad[3*i+j]; // Store the gradient elements into a temporary gradient
	  }
	}
	break;
      }
    }
    
    infile.close();  // Close the file
  }
  

  //PrintGradient("New Full MM Gradient",grad);

  //Rotating elements of gradient and forming a gradient reduced by symmetry

  //First make a key will indicate which atoms are in the gradient
  int gradkey[UniqueAtoms];
  int gradindex=0;
  for(int i=1; i<=NMon;i++){
    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
      if(Monomers[i].GetAtom(j).InAsymmetricUnit()){
	gradkey[gradindex] = Monomers[i].GetAtom(j).GetGlobalIndex();
	gradindex++;
      }
    }
  }
  
  //Reducing the size of the gradient restricted  under symmetry.
  for(int iMon=1;iMon<=NMon;iMon++){
    //Matrix RotMat = Monomers[iMon].GetRotationMatrix();
    // RotMat.Print("Rotation Matrix");
    for(int iatom=0; iatom<Monomers[iMon].GetNumberOfAtoms();iatom++){
      Matrix RotMat = Monomers[iMon].GetAtom(iatom).GetRotationMatrix();
      int SymGlobalIndex = Monomers[iMon].GetAtom(iatom).GetSymmetricalAtom();//Global index of the atom that the iAtom is symmetrical too

      // printf("iMon = %i iatom = %i\n",iMon,iatom);
      //RotMat.Print("RotMat");

      for(int k=0; k<UniqueAtoms;k++){
	if(SymGlobalIndex==gradkey[k]){//matching Global index of the atom it is symmetric to the index in the key
	  int GlobIndex = Monomers[iMon].GetAtom(iatom).GetGlobalIndex() - 1;//Global index used to get entree from TempGrad
	  
	  //printf("GlobIndex = %i SymGlobalIndex = %i\n",GlobIndex,SymGlobalIndex);
		 
	  //rotating elements of the Gradient and placing it into the reduced gradient 
	  grad[3*k]   += RotMat(0,0)*TempGrad[3*GlobIndex] + RotMat(0,1)*TempGrad[3*GlobIndex+1] + RotMat(0,2)*TempGrad[3*GlobIndex+2];
	  grad[3*k+1] += RotMat(1,0)*TempGrad[3*GlobIndex] + RotMat(1,1)*TempGrad[3*GlobIndex+1] + RotMat(1,2)*TempGrad[3*GlobIndex+2];
	  grad[3*k+2] += RotMat(2,0)*TempGrad[3*GlobIndex] + RotMat(2,1)*TempGrad[3*GlobIndex+1] + RotMat(2,2)*TempGrad[3*GlobIndex+2];
	}
      }
    }
  }
  //TempGrad.PrintGradient("Full gradiant");
  //grad.PrintGradient("Reduced gradiant");
  return grad;
}


// Main routine for reading in gradient 
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
// type 1 = qchem style, type 2 = tinker style
Vector Cluster::ReadCrystalGradient() {
  
  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  // Set up the filename, with the full path.  File is 'full.out'
  string filename;
  filename = path + "/full.out";
    
    // Open the output file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadCrystalGradient : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }

    int Natoms = GetTotalNumberOfAtoms();

    //Determine placement of assymetrical atoms in crystals gradient
    int placement[UniqueAtoms];
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      // Search string listing whether atoms are in asymmetrical unit
      if ( line.substr(0,29)==" ATOMS IN THE ASYMMETRIC UNIT" ) {
	getline(infile,line); // throw away first line
	getline(infile,line); // throw away second line
	int j = 0;
	for (int i=0;i<Natoms ;i++) {
          string InAsymmetric;
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  iss >> InAsymmetric; 
	  fflush(stdout);
	  if(InAsymmetric == "T"){
	    placement[j] = i;
            j++;
	  }
	}
	break;
      }
    }
    //printf("placement = [");
    //for(int i = 0;i<UniqueAtoms;i++)
    //  printf(" %i",placement[i]);
    //printf("]\n");

    //read gradient
    Rewind(infile);
    Vector Grad;
  if ( (Params::Parameters().IsPeriodic() ) && !Params::Parameters().DoFiniteDifferenceFreqs())
    Grad.Initialize(3*UniqueAtoms+6);
  else
    Grad.Initialize(3*UniqueAtoms);

    while ( !infile.eof() ) {
      getline(infile,line);
      if(line == " CARTESIAN FORCES IN HARTREE/BOHR (ANALYTICAL)"){
        int k = 0;
        getline(infile,line);  //throw away first line     
	for (int i=0;i<Natoms ;i++) {
          getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
          iss >> tmp; // throw away atomic number
          if(i == placement[k]){
	    for (int j=0;j<3;j++) {
	      iss >> Grad[3*k+j]; // Store the gradient elements into a temporary gradient
	    }
            k++;
          }
        }  
	break;    
      }
    }
    //Grad.PrintGradient("Grad");

    //atom key for gradient
    int gradkey[UniqueAtoms];
    int gradindex=0;
    for(int i=1; i<=NMon;i++){
      if(Monomers[i].GetSymmetryFactor()!=0){
        for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
	  if(Monomers[i].GetAtom(j).InAsymmetricUnit()){
	    gradkey[gradindex] = Monomers[i].GetAtom(j).GetGlobalIndex();
	    gradindex++;
	  }
        }
      }
    }

    //Multiply Gradient values by number of symmetrical atoms.
    int multiplier[UniqueAtoms];
    for(int i=0;i<UniqueAtoms;i++){
      multiplier[i] = 0;
    }
    for(int iMon = 1;iMon <= NMon; iMon++){
      for(int iAtom = 0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
        for(int i = 0;i<UniqueAtoms;i++){
          if(Monomers[iMon].GetAtom(iAtom).GetSymmetricalAtom() == gradkey[i]){
            multiplier[i]++;
          }
        }
      }
    }
    /*
    printf("multiplier = [");
    for(int i = 0;i<UniqueAtoms;i++)
      printf(" %i",multiplier[i]);
    printf("]\n");
    */
    for(int i=0;i<UniqueAtoms;i++){
      for(int j=0;j<3;j++)
       Grad[3*i+j] *= double(multiplier[i]);
    }

    //Grad.PrintGradient("Grad after multiplier");
    infile.close();  // Close the file
    return Grad;
}

Vector Cluster::ReadFiniteDifferenceCrystalGradient(){
  
  //Finite difference step size, must be the same as in CreateFiniteDifferenceCrystalJob
  double delta = 0.0002;//in fractional coordinates
  double lengthdelta = 0.001;//in angstroms
  double angledelta = 0.01;//in degrees

  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //Read Gradient
  Vector Grad;
  if ( (Params::Parameters().IsPeriodic() ) && !Params::Parameters().DoFiniteDifferenceFreqs())
    Grad.Initialize(3*UniqueAtoms+6);
  else
    Grad.Initialize(3*UniqueAtoms);
  
  //Read Nuclear coordinates finite difference step
  int i=0;
  if(!Params::Parameters().FreezeNuclearCoordinates()){
    for(int iMon=0;iMon<NMon;iMon++){
      for(int iAtom=0;iAtom<Monomers[iMon].GetNumberOfAtoms();iAtom++){
	if(Monomers[iMon].GetAtom(iAtom).InAsymmetricUnit()){
	  for(int xyz=0;xyz<3;xyz++){
	    if( !Monomers[iMon].GetAtom(iAtom).IsAtomFrozen() && //step atom that is frozen by symmetry
		(xyz!=1 ||!Monomers[iMon].GetAtom(iAtom).IsYLocked()) && //step y coord if locked by symmetry
		(xyz!=2 ||!Monomers[iMon].GetAtom(iAtom).IsZLocked()) ){ //step z coord if locked by symmetry
	      
	      //Plus job
	      char num[10];
	      sprintf(num,"%d",i);
	      string filename = path + "/full+" + num  + ".out";
	      
	      // Open the force file
	      ifstream infile;
	      infile.open( filename.c_str() );
	      if ( !infile.is_open() ) {
		printf("Cluster::ReadFiniteDifferenceCrystalGradient : Cannot open file '%s'\n",
		       filename.c_str());
		exit(1);
	      }
	      
	      // Read in the data
	      string line;
	      double energy;
	      printf("reading file %s\n",filename.c_str());
	      while ( !infile.eof() ) {
		getline(infile,line);
		string match = line.substr(0,20);
		if ( match==" TOTAL ENERGY + DISP" ) {
		  printf("Reading energy\n");
		  istringstream iss(line);
		  string tmp; 
		  // throw away "text tags"
		  iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; iss>> tmp;
		  iss >> energy; // read the energy
		}
		
	      }
	      infile.close();  // Close the file
	      Grad[i] = energy;

	      //Minus job
	      filename = path + "/full-" + num  + ".out";
	      
	      // Open the force file
	      infile.open( filename.c_str() );
	      if ( !infile.is_open() ) {
		printf("Cluster::ReadFiniteDifferenceCrystalGradient : Cannot open file '%s'\n",
		       filename.c_str());
		exit(1);
	      }
	      
	      // Read in the data
	      while ( !infile.eof() ) {
		getline(infile,line);
		string match = line.substr(0,20);
		if ( match==" TOTAL ENERGY + DISP" ) {
		  istringstream iss(line);
		  string tmp; 
		  // throw away "text tags"
		  iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; iss>> tmp;
		  iss >> energy; // read the energy
		}
		
	      }
	      infile.close();  // Close the file
	      Grad[i] -= energy; 

	    }
	    i++;
	  }
	}
      } 
    }
  }  
 
  //Grad.Print("Energy diff");
  Grad.Scale(1/(2*delta));
  //Grad.Print("frac Grad in Hartrees\n");
  Grad = ConvertBetweenFractionAndCartesianCoordinates(Grad,1,1,1);
  //Grad.Print("Grad in Hartrees/Angs\n");

  Grad.Scale(1/AngToBohr);
  //Grad.Print("Grad in Hartrees/Bohr\n");

  if(!Params::Parameters().FreezeUnitCellParams()){
    for(int i=0;i<6;i++){
      if((i!=1 || !IsBLocked()) && //Skip Length B if it's locked by symmetry
	 (i!=2 || !IsCLocked())  && //Skip Length C if it's locked by symmetry
	 (i!=3 || !IsAlphaLocked()) && //Skip Angle Alpha if it's locked by symmetry
	 (i!=4 || !IsBetaLocked()) && //Skip Angle Beta if it's locked by symmetry
	 (i!=5 || !IsGammaLocked()) ){ //Skip Angle Gamma if it's locked by symmetry

	string LatticeString;
	if(i==0) LatticeString = "a";
	else if(i==1) LatticeString = "b";
	else if(i==2) LatticeString = "c";
	else if(i==3) LatticeString = "alpha";
	else if(i==4) LatticeString = "beta";
	else if(i==5) LatticeString = "gamma";
	
	//Plus job
	string filename = path + "/full+" + LatticeString  + ".out";
	
	// Open the force file
	ifstream infile;
	infile.open( filename.c_str() );
	if ( !infile.is_open() ) {
	  printf("Cluster::ReadFiniteDifferenceCrystalGradient : Cannot open file '%s'\n",
		 filename.c_str());
	  exit(1);
	}
	
	// Read in the data
	string line;
	double energy;
	printf("reading file %s\n",filename.c_str());
	while ( !infile.eof() ) {
	  getline(infile,line);
	  string match = line.substr(0,20);
	  if ( match==" TOTAL ENERGY + DISP" ) {
	    istringstream iss(line);
	    string tmp; 
	    // throw away "text tags"
	    iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; iss>> tmp;
	    iss >> energy; // read the energy
	  }
	  
	}
	infile.close();  // Close the file
	Grad[3*UniqueAtoms+i] = energy;
	
	//Minus job
	filename = path + "/full-" + LatticeString + ".out";
	
	// Open the force file
	infile.open( filename.c_str() );
	if ( !infile.is_open() ) {
	  printf("Cluster::ReadFiniteDifferenceCrystalGradient : Cannot open file '%s'\n",
		 filename.c_str());
	  exit(1);
	}
	
	// Read in the data
	while ( !infile.eof() ) {
	  getline(infile,line);
	  string match = line.substr(0,20);
	  if ( match==" TOTAL ENERGY + DISP" ) {
	    istringstream iss(line);
	    string tmp; 
	    // throw away "text tags"
	    iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; iss>> tmp;
	    iss >> energy; // read the energy
	  }
	  
	}
	infile.close();  // Close the file
	Grad[3*UniqueAtoms+i] -= energy; 	

	if(i<3)
	  Grad[3*UniqueAtoms+i] *= 1/(2*lengthdelta*AngToBohr);
	else
	  Grad[3*UniqueAtoms+i] *= 1/(2*angledelta*DegreesToRadians);
      }         
    }      
  }
  //Grad.Print("Final Grad");
  return Grad;
}

/*
Vector Cluster::ReadQMCPGradient(int type) {

{
      getline(infile,line);
      // Search for 
      if ( line==" Nuclear forces:" ) {
	getline(infile,line); // throw away header line
	
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  for (int j=0;j<3;j++) {
	    iss >> TempGrad[3*i+j]; // Store the gradient elements into a temporary gradient
	  }
	}
	break;
      }
    }  string path = Params::Parameters().GetMMPath();
  int Natoms = GetTotalNumberOfAtoms();
  Vector grad(3*Natoms);

  if ( Params::Parameters().NeglectManyBody() )
    grad.Set();
  else {

    string filename;
    // Set up the filename, with the full path.  File is 'full.force'
    if (type == 2) // Tinker MM job
      filename = path + "/full.force";
    else
      filename = path + "/full.out";

    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadGradient : Cannot open file '%s'\n",
             filename.c_str());
      exit(1);
    }

    // Read in the data - search for the "Nuclear forces:" string
    string line;
    while ( !infile.eof() ) {
     getline(infile,line);
      // Search for final q-chem energy
      if ( line==" Nuclear forces:" ) {
        getline(infile,line); // throw away header line

        for (int i=0;i<Natoms;i++) {
          getline(infile,line);
          istringstream iss(line);
          string tmp;
          iss >> tmp; // throw away the atom index
          for (int j=0;j<3;j++) {
            iss >> grad[3*i+j]; // Store the gradient elements
          }
        }
        break;
      }
    }

    infile.close();  // Close the file

  }
  //PrintGradient("New Full MM Gradient",grad);
  return grad;





}
*/


// Compute Sum(Emon)
double Cluster::GetTotalOneBodyEnergy(string type) {

  double energy = 0.0;
  // Sanity check
  if (type != "QM" && type != "MM") {
    printf("Error: Cluster::GetTotalOneBodyEnergy(): Type should be QM or MM\n");
    exit(1);
  }

  if (type=="MM" && Params::Parameters().NeglectManyBody() ) {
    energy = 0.0;
    return energy;
  }
  // 1-body contribution
  for (int i=1;i<=NMon;i++) {
    if (type=="QM") {
      energy += Monomers[i].GetQMEnergy();
    }
    else {
      energy += Monomers[i].GetMMEnergy();
      //printf("m%d = %15.9f\n",i,Monomers[i].GetMMEnergy()*2625.5);
    }
  }
  return energy;
}

// Compute Sum(dEdim_int)
double Cluster::GetTotalTwoBodyEnergy(string type) {
  
  double dE = 0.0, dEes = 0.0, dEind = 0.0, dEdisp = 0.0; 
  double dEdispc6 = 0.0, dEdispc8 = 0.0, dEdispc10 = 0.0;
  double energy = 0.0, energy_es = 0.0, energy_ind = 0.0;
  double energy_disp = 0.0, energy_dispc6 = 0.0, energy_dispc8 = 0.0, energy_dispc10 = 0.0;

  int Ndamp = 0, Nfull = 0, Nzero = 0;

  // Sanity check
  if (type != "QM" && type != "MM") {
    printf("Error: Cluster::GetTotalTwoBodyEnergy(): Type should be QM or MM\n");
    exit(1);
  }
  
  if (type=="MM" && Params::Parameters().NeglectManyBody() ) {
    energy = 0.0;
    return energy;
  }
  
  // 2-body interaction contribution
  double c0 = Params::Parameters().GetLocalCutoff(0);
  double c1 = Params::Parameters().GetLocalCutoff(1);
  
  printf("GetTotalTwoBodyEnergy: c1=%.2f, c0=%.2f\n",c1,c0);
  
  for (int i=1;i<=NDim;i++) {
    // Handle symmetry
    float TotalSym = Dimers[i].GetSymmetryFactor();
    TotalSym += 0.5*Dimers[i].GetPeriodicSymmetryFactor();
    
    if (type=="QM"){
      dE = TotalSym*Dimers[i].GetQMIntEnergy(); 
    }else {
      dE = TotalSym*Dimers[i].GetMMIntEnergy();
      if (Params::Parameters().GetMMType() == 2) {
	// divide the es and ind 
	dEes = Dimers[i].GetMMElectrostaticEnergy();
	dEind = Dimers[i].GetMMInductionEnergy();
	dEdisp = Dimers[i].ComputeTwoBodyDispersion();
	
	// divide the C6, C8 and C10 disp contribution
	// dEdispc6 = Dimers[i].ComputeTwoBodyDispersionC6();
	// dEdispc8 = Dimers[i].ComputeTwoBodyDispersionC8();
	// dEdispc10 = Dimers[i].ComputeTwoBodyDispersionC10();
      }
      
      /*
      int iA = Dimers[i].GetIndexA();
      int iB = Dimers[i].GetIndexB();
      printf("E(%3d %3d ) = %9.3f  E(%3d) = %9.3f E(%3d) = %9.3f  Eint = %9.3f\n",
	     iA,iB,Dimers[i].GetMMEnergy()*2625.5, 
	     iA, Dimers[i].GetMonomerA().GetMMEnergy()*2625.5,
	     iB, Dimers[i].GetMonomerB().GetMMEnergy()*2625.5,
	     Dimers[i].GetMMIntEnergy()*2625.5);
      */
    }
      
    if (Params::Parameters().GetMMType() == 2) {
      dEes *= TotalSym;
      dEind *= TotalSym;
      dEdisp *= TotalSym; 
      
      //dEdispc6 *= Dimers[i].GetSymmetryFactor();
      //dEdispc8 *= Dimers[i].GetSymmetryFactor();
      //dEdispc10 *= Dimers[i].GetSymmetryFactor();
    }
    if ( Params::Parameters().DoLocal2BodyTruncation() ) {
      double damp = Dimers[i].GetDampingFactor(c0,c1);
      
      
      //printf("d(%i,%i) damp %15.9f\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),damp);
      //printf("%s: E(raw) = %15.9f, E(damped) = %15.9f\n",type.c_str(),dE,dE*damp);
      Nzero += Dimers[i].GetNzero(c0,c1);
      Nfull += Dimers[i].GetNfull(c0,c1);
      Ndamp += Dimers[i].GetNdamp(c0,c1);
      
      /* the problem with following code is that when r approx. equals c1,
	 we can still get a damping function of 1.0 which has to be avoided.
	 Note that the energy does not get affected.
	 Nzero,Nfull and Ndamp should be computed separately to avoid this problem.
      */
      /*
	if (damp==0.0)
	Nzero++;
	else if (damp==1.0)
	Nfull++;
	else
	  Ndamp++;
      */
      
      dE *= damp;
      dEes *= damp;
      dEind *= damp;
      dEdisp *= damp;
      
      // dEdispc6 *= damp;
      // dEdispc8 *= damp;
      // dEdispc10 *= damp;
      
    }
    
    energy += dE;
    energy_es += dEes;
    energy_ind += dEind;
    energy_disp += dEdisp;
    
    
    
    
    // energy_dispc6 += dEdispc6;
    // energy_dispc8 += dEdispc8;
    // energy_dispc10 += dEdispc10;
  }
  // If periodic: additional 2-body interaction contributions
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      if (type=="QM"){
	dE = DimerImages[i].GetQMIntEnergy();
      }else {
	dE = DimerImages[i].GetMMIntEnergy();
	// divide the es and ind 
	if (Params::Parameters().GetMMType() == 2) {
	  dEes = DimerImages[i].GetMMElectrostaticEnergy();
	  dEind = DimerImages[i].GetMMInductionEnergy();
	  dEdisp = DimerImages[i].ComputeTwoBodyDispersion();
	  // divide the C6, C8 and C10 disp contribution
	  //  dEdispc6 = DimerImages[i].ComputeTwoBodyDispersionC6();
	  // dEdispc8 = DimerImages[i].ComputeTwoBodyDispersionC8();
	  // dEdispc10 = DimerImages[i].ComputeTwoBodyDispersionC10();
	}	  
	  //int iA = DimerImages[i].GetIndexA();
	  //int iB = DimerImages[i].GetIndexB();
	  //int *K_vector = NULL;
	  //K_vector = DimerImages[i].GetImageCell();
	
	  //printf("E(%3d %3d ) = %9.3f  E(%3d) = %9.3f E(%3d) = %9.3f  Eint = %9.3f  Ref = %d, k=(%d,%d,%d)  sym_fac = %d\n",
	  //iA,iB,DimerImages[i].GetMMEnergy()*2625.5, 
	  //iA, DimerImages[i].GetMonomerA().GetMMEnergy()*2625.5,
	  //iB, DimerImages[i].GetMonomerB().GetMMEnergy()*2625.5,
	  //DimerImages[i].GetMMIntEnergy()*2625.5,
	  //DimerImages[i].GetReferenceMonomerIndex(),K_vector[0], K_vector[1], K_vector[2], DimerImages[i].GetSymmetryFactor());
      }
      // Handle symmetry
      dE *= DimerImages[i].GetSymmetryFactor();
      if (Params::Parameters().GetMMType() == 2) {
	dEes *= DimerImages[i].GetSymmetryFactor();
	dEind *= DimerImages[i].GetSymmetryFactor();
	dEdisp *= DimerImages[i].GetSymmetryFactor();
	
	// dEdispc6 *= DimerImages[i].GetSymmetryFactor();
	// dEdispc8 *= DimerImages[i].GetSymmetryFactor();
	// dEdispc10 *= DimerImages[i].GetSymmetryFactor();
      }
      
      double damp = DimerImages[i].GetDampingFactor(c0,c1);
      //printf("d(%i,%i) damp %15.9f\n",DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB(),damp);


      //printf("%s: E(raw) = %15.9f, E(damped) = %15.9f\n",
      //       type.c_str(),0.5*dE,0.5*dE*damp);
      //printf("local damping factor = %f\n",damp);
      dE *= damp;
      energy += 0.5*dE;
      if (Params::Parameters().GetMMType() == 2) {
	dEes *= damp;
	dEind *= damp;
	dEdisp *= damp;
	
	// dEdispc6 *= damp;
	// dEdispc8 *= damp;
	// dEdispc10 *= damp;
	
	energy_es += 0.5*dEes;
	energy_ind += 0.5*dEind; 
	energy_disp += 0.5*dEdisp;
	
	// energy_dispc6 += 0.5*dEdispc6;
	// energy_dispc8 += 0.5*dEdispc8;
	// energy_dispc10 += 0.5*dEdispc10;

	
	/*
	int iA = DimerImages[i].GetIndexA();
	int iB = DimerImages[i].GetIndexB();
	printf("E(%3d %3d ) = %9.3f  E(%3d) = %9.3f E(%3d) = %9.3f  Eint = %9.3f\n",
	       iA,iB,DimerImages[i].GetMMEnergy()*2625.5, 
	       iA, DimerImages[i].GetMonomerA().GetMMEnergy()*2625.5,
	       iB, DimerImages[i].GetMonomerB().GetMMEnergy()*2625.5,
	       DimerImages[i].GetMMIntEnergy()*2625.5);
	*/
	
      }
    }
  }
  
  // output some summary data
  // check type to ensure this only gets printed once
  if ( Params::Parameters().DoLocal2BodyTruncation() && type=="QM" ) {
    printf("----------------------------------------------\n");
    if ( Params::Parameters().IsPeriodic() ) 
      // this printout is less meaningful usually for periodic calcs
      printf("Local 2-Body Truncation within central unit cell:\n\n");
    else
      printf("Local 2-Body Truncation\n\n");
    printf("c1 = %.3f Ang, c0 = %.3f Ang\n\n",Params::Parameters().GetLocalCutoff(1),
	   Params::Parameters().GetLocalCutoff(0) );
    printf("Number of terms treated completely: %d\n",Nfull);
    printf("            Number of terms damped: %d\n",Ndamp);
    printf("         Number of terms neglected: %d\n",Nzero);
    printf("----------------------------------------------\n");
    printf("      Total number of 2-body terms: %d\n\n",Nfull+Nzero+Ndamp);
  }
  
  if ( Params::Parameters().GetMMType()==2 ) {
    printf("\nPartition Total 2-Body AIFF contributions:\n");
    printf("  Electrostatics = %9.3f kJ/mol\n",energy_es*HartreesToKJpermole);
    printf("       Induction = %9.3f kJ/mol\n",energy_ind*HartreesToKJpermole);
    printf("      Dispersion = %9.3f kJ/mol\n",energy_disp*HartreesToKJpermole);
    printf("           Total = %9.3f kJ/mol\n\n",
	   (energy_es+energy_ind+energy_disp)*HartreesToKJpermole);

    //printf("dimer aiff dispc6 = %9.3f kJ/mol\n",energy_dispc6*2625.5);
    //printf("dimer aiff dispc8 = %9.3f kJ/mol\n",energy_dispc8*2625.5);
    //printf("dimer aiff dispc10 = %9.3f kJ/mol\n",energy_dispc10*2625.5);
    
    // printf(" dimer aiff disp percentage analyses: C6/tot = %f, C8/tot = %f, C10/tot =%f\n",
    //	 energy_dispc6/energy_disp, energy_dispc8/energy_disp, energy_dispc10/energy_disp);
  }

  return energy;
  }

// E(HMBI) = Sum(EMon) + Sum(Edim_int)
double Cluster::ComputeHMBIEnergy() {
  double energy = 0;
  double Efull_MM, scafac; // scafac = 1 unless neglecting many-body terms
  double MM_1b,MM_2b;
  double enthalpyPV = 0;//for enthalpy

  if ( Params::Parameters().DoQMBenchmark() ) {
    double E_mb_QM;
    printf("\n---------------------------------------------------------------------\n");
    printf("Benchmark QM results on full cluster:\n");
    E_mb_QM = Energy_QM - GetTotalOneBodyEnergy("QM") 
      - GetTotalTwoBodyEnergy("QM");
    printf("   Full cluster QM energy = %15.9f hartrees\n",
	   Energy_QM);
    printf("   QM many-body energy = %.3f kJ/mol (%.3f kJ/mol per monomer)\n",
	   E_mb_QM*HartreesToKJpermole, E_mb_QM*HartreeToKJpMM);
    printf("---------------------------------------------------------------------\n");
  }

  double QM_1b, QM_2b, QM_full;

  if (!Params::Parameters().BuildForceFieldOnly() && !Params::Parameters().UseFullQMOnly()) {

    QM_1b = GetTotalOneBodyEnergy("QM");
    QM_2b = GetTotalTwoBodyEnergy("QM");
  }
  else if(Params::Parameters().UseFullQMOnly()){ //JLM
    QM_1b = 0.0; QM_2b = 0.0;
    QM_full = GetQMEnergy();
  }
  else {
    QM_1b = 0.0; QM_2b = 0.0;
  }

  if ( Params::Parameters().NeglectManyBody() || Params::Parameters().TinkerDebug() || Params::Parameters().UseFullQMOnly()) {
    printf("Neglecting many-body terms\n");
    scafac = 0.0; 
    Efull_MM = 0.0;
    MM_1b = 0.0; MM_2b = 0.0;
  }
  else if(Params::Parameters().UseFullMMOnly()){
    printf("Only including full MM \n");
    Efull_MM = GetMMEnergy();
    QM_1b = 0.0;
    QM_2b = 0.0;
    MM_1b = 0.0;
    MM_2b = 0.0;
  }else {
    Efull_MM = GetMMEnergy();
    MM_1b = GetTotalOneBodyEnergy("MM");
    MM_2b = GetTotalTwoBodyEnergy("MM");  
    //Efull_MM = 0.0;
    scafac = 1.0;
  }

  double MM_mb = Efull_MM - scafac*MM_1b - scafac*MM_2b;
  
  if(Params::Parameters().UseFullQMOnly())
    printf("QM full      = %15.9f hartrees",QM_full);
  else if(Params::Parameters().UseFullMMOnly())
    printf("MM full      = %9.3f kJ/mol\n",Efull_MM*2625.5);  
  else{
    printf("QM 1-body    = %15.9f hartrees\n",QM_1b);
    printf("QM 2-body    = %15.9f hartrees\n",QM_2b);
    printf("MM full      = %9.3f kJ/mol\n",Efull_MM*2625.5); 
    printf("MM 1-body    = %9.3f kJ/mol\n",MM_1b*2625.5);
    printf("MM 2-body    = %9.3f kJ/mol\n",MM_2b*2625.5);
    printf("MM many-body = %9.3f kJ/mol\n",MM_mb*2625.5);
  }
  
  if(Params::Parameters().UseFullQMOnly()){
    // Use only QM energy to get HMBI energy
    // Note: This is intended ONLY for planewave DFT 
    energy = QM_full;
  }
  else{
    // Combine QM 1- & 2-body with MM many-body to get HMBI energy
    energy = QM_1b + QM_2b + MM_mb;
  }

  //double E12_QM = GetTotalOneBodyEnergy("QM") + GetTotalTwoBodyEnergy("QM");
  double E12_QM = QM_1b + QM_2b;
  printf("\nQM 1+2-body energy = %15.9f hartrees\n",E12_QM);

  //double E_mb_MM = Efull_MM - scafac*GetTotalOneBodyEnergy("MM")
  //  - scafac*GetTotalTwoBodyEnergy("MM");
  if ( Params::Parameters().GetMMType() != 4)
    printf("MM many-body energy = %.3f kJ/mol (%.3f kJ/mol per monomer)\n",
	   MM_mb*HartreesToKJpermole, MM_mb*HartreeToKJpMM);

  //total Electronic Energy in HMBI
  printf("HMBI Electronic Energy = %9.3f hartrees\n",energy);

  //PV term
  if( (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()) && Params::Parameters().FindEnthalpy()){
    enthalpyPV = GetEnthalpyPV();
    energy += enthalpyPV;
  }

  //Vibrational Contribution
  if( (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()) && 
      (Params::Parameters().DoFreq() ||
       (Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) ){
    double Temperature = Params::Parameters().GetTemperature();

    printf("Vibrational Contribution (at %.0f K)  = %9.3f kJ/mol\n",
	   Temperature,Energy_Vib*2625.5);


    energy += Energy_Vib;

    /*
    printf("\n\n");
    printf("Entropy (at %.0f K)  = %15.9f J/(mol*K)\n",
	   Temperature,Entropy*1000*2625.5);
    printf("Enthalpy (at %.0f K)  = %15.9f Hartrees\n",
	   Temperature,(energy + Temperature*Entropy)); 
    printf("constant volume heat capacity (at %.0f K) = %.3f J/(mol*K)\n",Temperature,Cv);
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
      printf("Gruneisen parameter = %.3f\n", gruneisen_parameter);
    printf("Gibb's Free Energy (at %.0f K)  = %15.9f Hartrees\n\n",
	 Temperature,energy);
    */ 
  }

  return energy;
}

void Cluster::PrintGradient(string title, const Vector& grad) {
  
  //int Natoms = GetTotalNumberOfAtoms();

  printf("%s\n",title.c_str());
  //for (int i=0;i<Natoms;i++) {
  //  printf("%2s %15.10f %15.10f %15.10f\n",GetAtomicSymbol(i).c_str(),
  //  	   grad[3*i],grad[3*i+1],grad[3*i+2]);
  //}
  int row = 0;//counter for the entries in the gradient
  for(int i = 1;i <=NMon; i++)
      for(int j = 0;j < Monomers[i].GetNumberOfAtoms();j++)
	if(Monomers[i].GetAtom(j).InAsymmetricUnit()){
	  printf("%2s %15.10f %15.10f %15.10f\n", Monomers[i].GetAtom(j).GetSymbol().c_str(),
		 grad[3*row],grad[3*row+1],grad[3*row+2]);
	  fflush(stdout);
	row++;
      }

  if ( (((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly()) &&
       !(Params::Parameters().DoFiniteDifferenceFreqs())  ) {
    //printf("UnitCellLength %15.10f %15.10f %15.10f\n",
    //       grad[3*Natoms],grad[3*Natoms+1],grad[3*Natoms+2]);
    //printf("UnitCellAngles %15.10f %15.10f %15.10f\n",
    //       grad[3*Natoms+3],grad[3*Natoms+4],grad[3*Natoms+5]);
    printf("UnitCellLength %15.10f %15.10f %15.10f\n",
           grad[3*UniqueAtoms],grad[3*UniqueAtoms+1],grad[3*UniqueAtoms+2]);
    printf("UnitCellAngles %15.10f %15.10f %15.10f\n",
           grad[3*UniqueAtoms+3],grad[3*UniqueAtoms+4],grad[3*UniqueAtoms+5]);
  }

}

// dE/dr(HMBI) = dEfull_MM/dr + Sum(dEMon_QM/dr - dEMon_MM/dr) 
//              + Sum(dEdim_int_QM/dr - dEdim_int_MM/dr )
Vector Cluster::ComputeHMBIGradient() {
  /* The messy part of this routine is taking the "sparse" 
     gradients from each monomer and mapping their contributions 
     to the full gradient.  For example, the monomer with Na 
     atoms has a gradient of size 3*Na, while the full gradient 
     is 3*Ntot in size.

     Using space symmetry, the size of the gradient is reduced to 3*UniqueAtoms
  */

  //yoni:found earlier in cluster
  //Finding gradient for only the atoms of monomers that are not found to be symmetrical to other monomers.
  //int UniqueAtoms=0;
  //for(int i=1; i<=NMon;i++){
  //  if(Monomers[i].GetSymmetryFactor()!=0){
  //    UniqueAtoms += Monomers[i].GetNumberOfAtoms();
  //  }
  // }
  
  //printf("The gradient has a length of %i\n",3*UniqueAtoms);

  int Natoms=GetTotalNumberOfAtoms();
  //this key will indicate which atoms are in the gradient
  int gradkey[UniqueAtoms];
  int gradindex=0;
  for(int i=1; i<=NMon;i++){
    if(Monomers[i].GetSymmetryFactor()!=0){
      for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
	if(Monomers[i].GetAtom(j).InAsymmetricUnit()){
	  gradkey[gradindex] = Monomers[i].GetAtom(j).GetGlobalIndex();
	  gradindex++;
	}
      }
    }
  }

  
  printf("gradkey =[");
  for(int i=0;i<UniqueAtoms;i++)
    printf(" %i",gradkey[i]);
  printf("]\n");

  //printf("Number of unique atoms = %i\n",UniqueAtoms);

  //int Ntot = GetTotalNumberOfAtoms();
  //Vector grad(3*UniqueAtoms);
 // Vector FullGrad(3*Ntot);
  Vector grad;
  if (Params::Parameters().UseFullQMOnly()) {
    //grad.Initialize(3*Natoms+9); //JLM before symmetry
    grad.Initialize(3*UniqueAtoms+9);
  }
  else if ( (Params::Parameters().IsPeriodic()  && Params::Parameters().GetMMType() != 2)
       && !Params::Parameters().DoFiniteDifferenceFreqs()  ) {
    grad.Initialize(3*UniqueAtoms+6);
  }
  else {
    grad.Initialize(3*UniqueAtoms);
  }
  double scafac=1.0, scafac2=1.0; 
 
  if ((Params::Parameters().NeglectManyBody() || Params::Parameters().TinkerDebug()) && !Params::Parameters().UseFullQMOnly()) {
    scafac = 0.0; // neglect many-body terms
  }
  else {
    scafac = 1.0; // treat them as usual
  }

  // Contribution from full cluster MM gradient
  if (!Params::Parameters().NeglectManyBody()){
    grad = Grad_MM;
    //grad.Set();
  }else {
    // No MM case.  Set all MM gradients to zero.
    // This wastes memory, but it was simple.
    grad.Set();
    for (int i=1;i<=NMon;i++) {
      Vector tmp(3*Monomers[i].GetNumberOfAtoms());
      tmp.Set();
      Monomers[i].SetMMGradient(tmp);
    }
    for (int i=1;i<=NDim;i++) {
      Vector tmp(3*Dimers[i].GetNumberOfAtoms());
      tmp.Set();
      Dimers[i].SetMMGradient(tmp);
    }
    for(int i=1;i<=NDim_images;i++){
      Vector tmp(3*DimerImages[i].GetNumberOfAtoms());
      tmp.Set();
      DimerImages[i].SetMMGradient(tmp);
    }
  }

  // Contribution from full cluster QM gradient JLM
  if (Params::Parameters().UseFullQMOnly()){
   grad.Set();
    for (int i = 0; i < UniqueAtoms; i++){
      for( int k = 0; k<3; k++){
        grad[3*i + k] = Grad_QM[3*(gradkey[i]-1)+k];
      }
    }

    for (int j = 9; j > 0; j--){
      grad[grad.GetLength()-j] = Grad_QM[Grad_QM.GetLength()-j];
    }

    //grad = Grad_QM; // before symmetry
    //grad.PrintGradient("Current QM gradient");
    //PrintGradient("Current QM gradient",grad);
    //exit(0);
  }
  else if (Params::Parameters().BuildForceFieldOnly()) {
    //No QM case. Set all QM gradients to zero
    for (int i=1;i<=NMon;i++) {
      Vector tmp(3*Monomers[i].GetNumberOfAtoms());
      tmp.Set();
      Monomers[i].SetQMGradient(tmp);
    }
    for (int i=1;i<=NDim;i++) {
      Vector tmp(3*Dimers[i].GetNumberOfAtoms());
      tmp.Set();
      Dimers[i].SetQMGradient(tmp);
    }
    for(int i=1;i<=NDim_images;i++){
      Vector tmp(3*DimerImages[i].GetNumberOfAtoms());
      tmp.Set();
      DimerImages[i].SetQMGradient(tmp);
    }
  }

  //Case where only the full MM or the full QM contributes to the gradient
  if(Params::Parameters().UseFullMMOnly() || Params::Parameters().UseFullQMOnly()){
    for (int i=1;i<=NMon;i++) {
      Vector tmp(3*Monomers[i].GetNumberOfAtoms());
      tmp.Set();
      Monomers[i].SetQMGradient(tmp);
      Monomers[i].SetMMGradient(tmp);
    }
    for (int i=1;i<=NDim;i++) {
      Vector tmp(3*Dimers[i].GetNumberOfAtoms());
      tmp.Set();
      Dimers[i].SetQMGradient(tmp);
      Dimers[i].SetMMGradient(tmp);
    }
    for(int i=1;i<=NDim_images;i++){
      Vector tmp(3*DimerImages[i].GetNumberOfAtoms());
      tmp.Set();
      DimerImages[i].SetQMGradient(tmp);
      DimerImages[i].SetMMGradient(tmp);
    }
  }  

  /*
  //Hack, zero 2 Body interactions
  printf("Hack no 2-body interactions\n");
  for (int i=1;i<=NMon;i++) {
	Vector tmp(3*Monomers[i].GetNumberOfAtoms());
	tmp.Set();
	Monomers[i].SetMMGradient(tmp
  }
  
  for (int i=1;i<=NDim;i++) {
    Vector tmp(3*Dimers[i].GetNumberOfAtoms());
    tmp.Set();
    Dimers[i].SetQMGradient(tmp);
    Dimers[i].SetMMGradient(tmp);
  }
  for(int i=1;i<=NDim_images;i++){
    Vector tmp(3*DimerImages[i].GetNumberOfAtoms());
    tmp.Set();
    DimerImages[i].SetQMGradient(tmp);
    DimerImages[i].SetMMGradient(tmp);
  }
  */
  grad.Scale(scafac);

  //grad.PrintGradient("QM full gradient");

 //If we're only using a cluster no need for monomer and dimer contributions 
 //if (!Params::Parameters().UseFullQMOnly() && !Params::Parameters().UseFullMMOnly()){
 if (!Params::Parameters().UseFullQMOnly()){
  //1 body contributions
  for(int i=1;i<=NMon;i++){

    //Matrix Rot = Monomers[i].GetRotationMatrix();

    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){

      Matrix Rot = Monomers[i].GetAtom(j).GetRotationMatrix();
      int SymGlobalIndex = Monomers[i].GetAtom(j).GetSymmetricalAtom();

      //int GlobalIndex = Monomers[i].GetAtom(j).GetSymmetricalAtom();
      for(int k=0; k<UniqueAtoms;k++){
	if(SymGlobalIndex==gradkey[k]){//matching Global index of the atom it is symmetric to the index in the key
	  for(int xyz=0;xyz<3;xyz++){
	    //printf("%i matches %i\n",
	    //	 Monomers[i].GetAtom(j).GetSymmetricalAtom(),gradkey[k]);
	    //printf("atom %s\n",Monomers[i].GetAtom(j).GetSymbol().c_str());
	    //printf("using grad %f %f %f\n",
	     //	 Monomers[i].GetQMGradient()[3*j],Monomers[i].GetQMGradient()[3*j+1],Monomers[i].GetQMGradient()[3*j+2]);
	    //QM contribution
	    grad[3*k] += Rot(0,xyz)*Monomers[i].GetQMGradient()[3*j+xyz];
	    grad[3*k+1] +=  Rot(1,xyz)*Monomers[i].GetQMGradient()[3*j+xyz];
	    grad[3*k+2] +=  Rot(2,xyz)*Monomers[i].GetQMGradient()[3*j+xyz];	  

	    //MM contribution
	    grad[3*k] -= scafac*Rot(0,xyz)*Monomers[i].GetMMGradient()[3*j+xyz];
	    grad[3*k+1] -=  scafac*Rot(1,xyz)*Monomers[i].GetMMGradient()[3*j+xyz];
	    grad[3*k+2] -=  scafac*Rot(2,xyz)*Monomers[i].GetMMGradient()[3*j+xyz];
         }//end of loop over xyz
	}//end of match global index to the gradkey
      }//end of loop over k
    }//end of loop over j
  }//end of loop over i

  //grad.PrintGradient(" Including One-Body Gradiant");  

  // grab the cutoffs for the QM -> MM damping
  double c0 = Params::Parameters().GetLocalCutoff(0);
  double c1 = Params::Parameters().GetLocalCutoff(1);

  //printf("finding two-body contribution\n");


  //using monomer symmetry
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();


  // if ( !(Params::Parameters().DoCounterpoise() ) ) {
  // 2-body interaction contributions
  for (int i=1;i<=NDim;i++) {

    float symfac = Dimers[i].GetSymmetryFactor();
    symfac += Dimers[i].GetPeriodicSymmetryFactor()/2;
    double separation = Dimers[i].GetDimerSeparation();
    if(symfac!=0 &&
       (!Params::Parameters().DoLocal2BodyTruncation() ||separation <= c0)){
      int Na = Dimers[i].GetMonomerA().GetNumberOfAtoms();
      int Nb = Dimers[i].GetMonomerB().GetNumberOfAtoms();
      int indexA = Dimers[i].GetIndexA();
      int indexB = Dimers[i].GetIndexB();
       


      //printf("D(%i,%i) symfac =%f\n",
      //	      indexA,indexB,symfac);

      
      //Find the local index for the symmetrical atom 
      int SymA = Monomers[indexA].GetSymmetricalMonomer();
      int SymB = Monomers[indexB].GetSymmetricalMonomer();
      
      //int startA = key[indexA];
      //int startB = key[indexB];
      
      //Find the damping factor but only if doing local 2 body truncation
      double damp;
      Vector DampGrad(3*(Na+Nb));  
      if(Params::Parameters().DoLocal2BodyTruncation()){
	damp = Dimers[i].GetDampingFactor(c0,c1);
	DampGrad = Dimers[i].GetSpatialDampingFunctionGradient(c0,c1);
      } 
      else{
	damp = 1.0;
	DampGrad.Set();
      }
      //printf("damp = %f\n",damp);
      
      //MonomerA piece of the gradient
      for(int j=0;j<Na;j++){
	
        //Matrix Rot = Monomers[indexA].GetRotationMatrix();
	Matrix Rot = Monomers[indexA].GetAtom(j).GetRotationMatrix();

	//use the Global Index find which atoms the atoms in a are symmetrical to
	int SymGlobalIndex = Monomers[indexA].GetAtom(j).GetSymmetricalAtom();

	//This is used to find which atom in the symmetrical monomer the atom is symmetrical to.
	//We need to know the local index of the atom to subtract the correct one body gradient
	bool MatchLocal =0;
	int SymAtom = 0;
	while(MatchLocal==0 && SymAtom<Na){

	  //printf("does %i that is sym to %i match %i\n",Monomers[indexA].GetAtom(j).GetGlobalIndex(),
	  //  Monomers[indexA].GetAtom(j).GetSymmetricalAtom(),Monomers[SymA].GetAtom(SymAtom).GetGlobalIndex());
	  
	  if(Monomers[indexA].GetAtom(j).GetSymmetricalAtom() == Monomers[SymA].GetAtom(SymAtom).GetGlobalIndex()){

	    //printf("Match Found\n");
	    MatchLocal=1;
	  }else if(!UseMonSym)//this process is unnecessary if monomer symmetry is not being use
	    MatchLocal=1;
	  else
	    SymAtom++;
	}
	if(MatchLocal==0){
	  printf("ERROR::Cluster::ComputeHMBIGradient() cannot find which atom %i is symmetrical\n",
		 Monomers[indexA].GetAtom(j).GetGlobalIndex());
	  exit(0);
	}
	
	for(int k=0;k<UniqueAtoms;k++){

	  if(gradkey[k]==SymGlobalIndex){

	    for(int xyz=0;xyz<3;xyz++){
	      //added the dimer contribution to the pair-wise interactions
	      grad[3*k] += damp*symfac*Rot(0,xyz)*(Dimers[i].GetQMGradient()[3*j+xyz]-scafac*Dimers[i].GetMMGradient()[3*j+xyz]);
	      grad[3*k+1] += damp*symfac*Rot(1,xyz)*(Dimers[i].GetQMGradient()[3*j+xyz]-scafac*Dimers[i].GetMMGradient()[3*j+xyz]);
	      grad[3*k+2] += damp*symfac*Rot(2,xyz)*(Dimers[i].GetQMGradient()[3*j+xyz]-scafac*Dimers[i].GetMMGradient()[3*j+xyz]);
	      
	      /*
	      if(3*k == 0 && fabs(Rot(0,xyz)) > 0.00001){
		printf("MonA\n");
		printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		printf("symfac = %f\n",symfac);
		printf("damp = %f\n",damp);
		printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		printf("3*j = %i\n",3*j);
		printf("QM Dimer Grad = %f\n",Dimers[i].GetQMGradient()[3*j+xyz]);
		printf("Sum Grad = %f\n\n",grad[3*k]);
	      }
	      */
	      
	      //Subtracting the monomer contribution.
	      //If counterpoise corrected, the monomer contribution is already subtracted for the QM.
	      //If monomer symmetry is in use, then the grad of the symmetry atoms will be used
	      if(!Params::Parameters().DoCounterpoise()){
		if(UseMonSym){
		  grad[3*k] -=  damp*symfac*Rot(0,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		  grad[3*k+1] -=  damp*symfac*Rot(1,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		  grad[3*k+2] -=  damp*symfac*Rot(2,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		}else{
		  grad[3*k] -=  damp*symfac*Rot(0,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		  grad[3*k+1] -=  damp*symfac*Rot(1,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		  grad[3*k+2] -=  damp*symfac*Rot(2,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		}
	      }
	      
	      /*
	      if(3*k == 0 && xyz == 0){
	      printf("MonA\n");
	      printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
	      printf("symfac = %f\n",symfac);
	      printf("damp = %f\n",damp);
	      printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
	      printf("3*j = %i\n",3*j);
	      printf("QM Monomer Grad = %f\n",Monomers[indexA].GetQMGradient()[3*j+xyz]);
	      printf("Sum Grad = %f\n\n",grad[3*k]);
	    }
	      */
	      
	      //For the monomer MM of the two-body interaction, a negative is being subtrated therefore it is positive.
	      //If the MM uses qchem and counterpoise corected, the two the monomer contribution is already subtrated
	      if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3){
		if(UseMonSym){
		  grad[3*k] +=  scafac*damp*symfac*Rot(0,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		  grad[3*k+1] +=  scafac*damp*symfac*Rot(1,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		  grad[3*k+2] += scafac*damp*symfac*Rot(2,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		}else{
		  grad[3*k] +=  scafac*damp*symfac*Rot(0,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		  grad[3*k+1] +=  scafac*damp*symfac*Rot(1,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		  grad[3*k+2] +=  scafac*damp*symfac*Rot(2,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		}
	      }
	      
	      //gradient for the damping function
	      grad[3*k] += Rot(0,xyz)*DampGrad[3*j+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());
	      grad[3*k+1] += Rot(1,xyz)*DampGrad[3*j+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());
	      grad[3*k+2] += Rot(2,xyz)*DampGrad[3*j+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());

	      /*
	      if(3*k == 0 && fabs(Rot(0,xyz)) > 0.00001){
		printf("MonA\n");
		printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		printf("DampGrad = %f\n",DampGrad[3*j+xyz]);
		printf("symfac = %f\n",symfac);
		printf("3*j = %i\n",3*j);
		printf("QM Int = %f\n",Dimers[i].GetQMIntEnergy());
		printf("Sum Grad = %f\n\n",grad[3*k]);
	      }
	      */

	    }
	  }
	}
	
      }
      //MonomerB piece of the gradient
      for(int j=0;j<Nb;j++){
	int k=Na+j;

        //Matrix Rot = Monomers[indexB].GetRotationMatrix();
	Matrix Rot = Monomers[indexB].GetAtom(j).GetRotationMatrix();


	//use the Global Index find which atoms the atoms in a are symmetrical to
	int SymGlobalIndex = Monomers[indexB].GetAtom(j).GetSymmetricalAtom();
	
	//These varibles are used to find which atom in the symmetrical monomer the atom is symmetrical to.
	//We need to know the local index of the atom to subtract the correct one body gradient
	bool MatchLocal =0;
	int SymAtom = 0;
	
	while(MatchLocal==0 && SymAtom<Nb){
	  //  printf("does %i that is sym to %i match %i\n",Monomers[indexB].GetAtom(j).GetGlobalIndex(),
	  // Monomers[indexB].GetAtom(j).GetSymmetricalAtom(),Monomers[SymB].GetAtom(SymAtom).GetGlobalIndex());
	  
	  if(Monomers[indexB].GetAtom(j).GetSymmetricalAtom() == Monomers[SymB].GetAtom(SymAtom).GetGlobalIndex()){
	    //printf("Match Found\n");
	    MatchLocal=1;
	  }else if(!UseMonSym)//this process is unnecessary if monomer symmetry is not being use
	    MatchLocal=1;
	  else
	    SymAtom++;
	}
	if(MatchLocal==0){
	  printf("ERROR::Cluster::ComputeHMBIGradient() cannot find which atom %i is symmetrical\n",
		 Monomers[indexB].GetAtom(j).GetAtomIndex());
	  exit(0);
	}
	
	for(int l=0;l<UniqueAtoms;l++){
	  if(gradkey[l]==SymGlobalIndex){

            for(int xyz=0;xyz<3;xyz++){
	      //added the dimer contribution to the pair-wise interactions
	      grad[3*l] += damp*symfac*Rot(0,xyz)*(Dimers[i].GetQMGradient()[3*k+xyz]-scafac*Dimers[i].GetMMGradient()[3*k+xyz]);
	      grad[3*l+1] += damp*symfac*Rot(1,xyz)*(Dimers[i].GetQMGradient()[3*k+xyz]-scafac*Dimers[i].GetMMGradient()[3*k+xyz]);
	      grad[3*l+2] += damp*symfac*Rot(2,xyz)*(Dimers[i].GetQMGradient()[3*k+xyz]-scafac*Dimers[i].GetMMGradient()[3*k+xyz]);
	      
	      /*
	      if(3*l == 0 && fabs(Rot(0,xyz)) > 0.00001){
		printf("MonB\n");
		printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		printf("symfac = %f\n",symfac);
		printf("damp = %f\n",damp);
		printf("3*l = %i\n",3*l);
		printf("QM Dimer Grad = %f\n",Dimers[i].GetQMGradient()[3*k+xyz]);
		printf("Sum Grad = %f\n\n",grad[3*l]);
	      }
	      */

	      //subtracting the monomer contribution
	      //if counterpoise corrected, the monomer contribution is already subtracted
	      //If monomer symmetry is in use, then the grad of the symmetry atoms will be used
	      if(!Params::Parameters().DoCounterpoise()){
	        if(UseMonSym){
	   	  grad[3*l] -=  damp*symfac*Rot(0,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
		  grad[3*l+1] -=  damp*symfac*Rot(1,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
		  grad[3*l+2] -=  damp*symfac*Rot(2,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
	        }else{
		  grad[3*l] -=  damp*symfac*Rot(0,xyz)*Monomers[indexB].GetQMGradient()[3*j+xyz];
		  grad[3*l+1] -=  damp*symfac*Rot(1,xyz)*Monomers[indexB].GetQMGradient()[3*j+xyz];
		  grad[3*l+2] -=  damp*symfac*Rot(2,xyz)*Monomers[indexB].GetQMGradient()[3*j+xyz];
	        }
	      }
	      /*
	    if(3*l == 0 && xyz == 0){
	      printf("MonB\n");
	      printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
	      printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
	      printf("symfac = %f\n",symfac);
	      printf("damp = %f\n",damp);
	      printf("3*k = %i\n",3*j);
	      printf("QM Monomer Grad = %f\n",Monomers[indexB].GetQMGradient()[3*j]);
	      printf("Sum Grad = %f\n\n",grad[3*l]);
	    }
	    */

	      //For the monomer MM of the two-body interaction, a negative is being subtracted therefore it is positive.
	      //If the MM uses qchem and counterpoise corrected, the two monomer contributions are already subtracted
	      if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3){
	        if(UseMonSym){
		  grad[3*l] +=  scafac*damp*symfac*Rot(0,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
		  grad[3*l+1] +=  scafac*damp*symfac*Rot(1,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
		  grad[3*l+2] += scafac*damp*symfac*Rot(2,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
	        }else{
		  grad[3*l] +=  scafac*damp*symfac*Rot(0,xyz)*Monomers[indexB].GetMMGradient()[3*j+xyz];
		  grad[3*l+1] +=  scafac*damp*symfac*Rot(1,xyz)*Monomers[indexB].GetMMGradient()[3*j+xyz];
		  grad[3*l+2] += scafac*damp*symfac*Rot(2,xyz)*Monomers[indexB].GetMMGradient()[3*j+xyz];
	        }
	      }

	      //gradient for the damping function
	      grad[3*l] += Rot(0,xyz)*DampGrad[3*k+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());
	      grad[3*l+1] += Rot(1,xyz)*DampGrad[3*k+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());
	      grad[3*l+2] += Rot(2,xyz)*DampGrad[3*k+xyz]*symfac*(Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());

	      /*
	      if(3*l == 0 && fabs(Rot(0,xyz)) > 0.00001){
		printf("MonB");
		printf("D(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		printf("symfac = %f\n",symfac);
		printf("DampGrad = %f\n",DampGrad[3*k]);
		printf("3*l = %i\n",3*k);
		printf("Int QM= %f\n",Dimers[i].GetQMIntEnergy());
		printf("Sum Grad = %f\n\n",grad[3*l]);
	      }
	      */

            }//end of loop over xyz
	  }//end of match global index to the gradkey
	}//end of loop over l
      }//end of loop over j
    }//end of if statement for symmetry factor != 0 and for truncation
  }//end of loop over i

  //grad.PrintGradient("Gradient after two body contribution");
  
  //printf("finding periodic two body gradient\n");
  
  //2-body contributions due to periodic boundary conditions
  if (Params::Parameters().IsPeriodic() ) {
       // 2-body interaction contributions
    for (int i=1;i<=NDim_images;i++) { 
      
      double separation = DimerImages[i].GetDimerSeparation();
      if(separation <= c0){           
	int indexA = DimerImages[i].GetIndexA();
	 int indexB = DimerImages[i].GetIndexB();
	 int indexB_ref = DimerImages[i].GetReferenceMonomerIndex();
	 int Na = DimerImages[i].GetMonomerA().GetNumberOfAtoms();
	 int Nb = DimerImages[i].GetMonomerB().GetNumberOfAtoms();   
	 // int startA = key[indexA];
	 //int startB = key[indexB];
	 
	 //printf("D(%i,%i) refB =%i\n",
	 //	      indexA,indexB, indexB_ref);
	 //fflush(stdout);
	 
	 //Find the local index for the symmetrical atom 
	 int SymA = Monomers[indexA].GetSymmetricalMonomer();
	 int SymB = Monomers[indexB_ref].GetSymmetricalMonomer();
	 
	 double damp = DimerImages[i].GetDampingFactor(c0,c1);
	 Vector DampGrad = DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1);
	 //   printf("damp = %f\n",damp);
	 //fflush(stdout);
	 
	 int symfac = DimerImages[i].GetSymmetryFactor();
	 
	 // piece arising from first monomer
	 for(int j=0;j<Na;j++){
	   
	   //use the Global Index find which atoms the atoms in monomer A are symmetrical to
	   int SymGlobalIndex = Monomers[indexA].GetAtom(j).GetSymmetricalAtom();
	   
           //Matrix Rot = Monomers[indexA].GetRotationMatrix();
	   Matrix Rot = Monomers[indexA].GetAtom(j).GetRotationMatrix();
	   //These varibles are used to find which atom in the symmetrical monomer the atom is symmetrical to.
	   //We need to know the local index of the atom to subtract the correct one body gradient
	   bool MatchLocal =0;
	   int SymAtom = 0;
	   
	   while(MatchLocal==0 && SymAtom<Na){
	   //printf("does %i that is sym to %i match %i\n",Monomers[indexA].GetAtom(j).GetGlobalIndex(),
	   //  Monomers[indexA].GetAtom(j).GetSymmetricalAtom(),Monomers[SymA].GetAtom(SymAtom).GetGlobalIndex());
	     
	     if(Monomers[indexA].GetAtom(j).GetSymmetricalAtom() == Monomers[SymA].GetAtom(SymAtom).GetGlobalIndex()){
	       //printf("Match Found\n");
	       MatchLocal=1;
	     }else if(!UseMonSym)//this process is unnecessary if monomer symmetry is not being use
	       MatchLocal=1;
	     else
	       SymAtom++; 
	   }
	   //untested check
	   if(MatchLocal==0){
	     printf("Error::Cluster::ComputeHMBIGradient() cannot find which atom %i is symmetrical\n",
		    Monomers[indexA].GetAtom(j).GetAtomIndex());
	   }
	   
	   for(int k=0;k<UniqueAtoms;k++){
	     if(gradkey[k]==SymGlobalIndex){
	       for(int xyz=0;xyz<3;xyz++){

	       //added the dimer contribution to the pair-wise interactions
	       grad[3*k] += 0.5*damp*symfac*Rot(0,xyz)*(DimerImages[i].GetQMGradient()[3*j+xyz]-scafac*DimerImages[i].GetMMGradient()[3*j+xyz]);
	       grad[3*k+1] += 0.5*damp*symfac*Rot(1,xyz)*(DimerImages[i].GetQMGradient()[3*j+xyz]-scafac*DimerImages[i].GetMMGradient()[3*j+xyz]);
	       grad[3*k+2] += 0.5*damp*symfac*Rot(2,xyz)*(DimerImages[i].GetQMGradient()[3*j+xyz]-scafac*DimerImages[i].GetMMGradient()[3*j+xyz]);
	       
	       /*
		 if(k==2){
		   printf("MonA: Added to atom %i from atom %i, %f %f %f \n",gradkey[k],Monomers[indexA].GetAtom(j).GetGlobalIndex(),
			  DimerImages[i].GetQMGradient()[3*j],DimerImages[i].GetQMGradient()[3*j+1],DimerImages[i].GetQMGradient()[3*j+2]);
		   printf("damp = %f\n",damp);
		   printf("symfac = %i\n",symfac);
		   printf("grad = %f\n\n",grad[3*k]);
		   fflush(stdout);
		 }
	       */

	       //Subtracting the monomer contribution.
	       //If counterpoise corrected, the monomer contribution is already subtracted for the QM.
	       //Since symmetry is in use, the grad of the symmetry atoms will be used
	       if(!Params::Parameters().DoCounterpoise()){
		 if(UseMonSym){
		   grad[3*k] -=  0.5*damp*symfac*Rot(0,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		   grad[3*k+1] -=  0.5*damp*symfac*Rot(1,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		   grad[3*k+2] -=  0.5*damp*symfac*Rot(2,xyz)*Monomers[SymA].GetQMGradient()[3*SymAtom+xyz];
		 }else{
		   grad[3*k] -=  0.5*damp*symfac*Rot(0,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		   grad[3*k+1] -=  0.5*damp*symfac*Rot(1,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		   grad[3*k+2] -=  0.5*damp*symfac*Rot(2,xyz)*Monomers[indexA].GetQMGradient()[3*j+xyz];
		 }
	       }

	       /*
		 if(k==2){
		   printf("     Monomer grad: %f %f %f\n",
			  Monomers[indexA].GetQMGradient()[3*j],Monomers[indexA].GetQMGradient()[3*j+1],Monomers[indexA].GetQMGradient()[3*j+2]);
		   printf("grad = %f\n\n",grad[3*k]);
		 }
	       */
	       
	       //For the monomer MM of the two-body interaction, a negative is being subtracted therefore it is positive.
	       //If the MM uses qchem and counterpoise corected, the two the monomer contribution is already subtrated
	       if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3){
		 if(UseMonSym){
		   grad[3*k] +=  0.5*scafac*damp*symfac*Rot(0,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		   grad[3*k+1] +=  0.5*scafac*damp*symfac*Rot(1,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		   grad[3*k+2] += 0.5*scafac*damp*symfac*Rot(2,xyz)*Monomers[SymA].GetMMGradient()[3*SymAtom+xyz];
		 }else{
		   grad[3*k] +=  0.5*scafac*damp*symfac*Rot(0,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		   grad[3*k+1] +=  0.5*scafac*damp*symfac*Rot(1,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		   grad[3*k+2] += 0.5*scafac*damp*symfac*Rot(2,xyz)*Monomers[indexA].GetMMGradient()[3*j+xyz];
		 }
	       }
	       
	       //gradient for the damping function
	       grad[3*k] += 0.5*Rot(0,xyz)*DampGrad[3*j+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       grad[3*k+1] += 0.5*Rot(1,xyz)*DampGrad[3*j+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       grad[3*k+2] += 0.5*Rot(2,xyz)*DampGrad[3*j+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       
	       /*
	       if(k==2){
		 printf("DampGrad = %f\n",DampGrad[3*j]);
		 printf("QMInt = %f\n",DimerImages[i].GetQMIntEnergy());
		 printf("grad = %f\n\n",grad[3*k]);
		 fflush(stdout);
	       }
	       */

	       //QM
	       //grad[3*k] += 0.5*damp*symfac*(DimerImages[i].GetQMGradient()[3*j]-Monomers[SymA].GetQMGradient()[3*SymAtom]);
	       //grad[3*k+1] += 0.5*damp*symfac*(DimerImages[i].GetQMGradient()[3*j+1]-Monomers[SymA].GetQMGradient()[3*SymAtom+1]);
	       //grad[3*k+2] += 0.5*damp*symfac*(DimerImages[i].GetQMGradient()[3*j+2]-Monomers[SymA].GetQMGradient()[3*SymAtom+2]);
	       
	       //MM
	       //grad[3*k] -= 0.5*damp*symfac*scafac*(DimerImages[i].GetMMGradient()[3*j]-Monomers[SymA].GetMMGradient()[3*SymAtom]);
	       //grad[3*k+1] -= 0.5*damp*symfac*scafac*(DimerImages[i].GetMMGradient()[3*j+1]-Monomers[SymA].GetMMGradient()[3*SymAtom+1]);
	       //grad[3*k+2] -= 0.5*damp*symfac*scafac*(DimerImages[i].GetMMGradient()[3*j+2]-Monomers[SymA].GetMMGradient()[3*SymAtom+2]);	
	       
               }//end of loop over xyz
	     }//end of if statement matching atom symmetry to the gradkey
	   }//end loop over k
	 }//end loop over j
	 
	 // Piece arising from second monomer 
	 for(int j=0;j<Nb;j++){
	   
	   int k=Na+j;
	   //use the Global Index find which atoms the atoms in a are symmetrical to
	   int SymGlobalIndex = Monomers[indexB_ref].GetAtom(j).GetSymmetricalAtom();
  	   
           Matrix Rot = Monomers[indexB_ref].GetAtom(j).GetRotationMatrix();
 
	   //These varibles are used to find which atom in the symmetrical monomer the atom is symmetrical to.
	   //We need to know the local index of the atom to subtract the correct one body gradient
	   bool MatchLocal =0;
	   int SymAtom = 0;
	   
	   while(MatchLocal==0 && SymAtom<Nb){
	     
	     // printf("does %i that is sym to %i match %i\n",Monomers[indexB_ref].GetAtom(j).GetGlobalIndex(),
	     //	 Monomers[indexB_ref].GetAtom(j).GetSymmetricalAtom(),Monomers[SymB].GetAtom(SymAtom).GetGlobalIndex());
	     
	     if(Monomers[indexB_ref].GetAtom(j).GetSymmetricalAtom() == Monomers[SymB].GetAtom(SymAtom).GetGlobalIndex()){
	       //printf("Match Found\n");
	       MatchLocal=1;
	     }else if(!UseMonSym)//this process is unnecessary if monomer symmetry is not being use
	       MatchLocal=1;
	     else
	       SymAtom++;
	   }
	   if(MatchLocal==0){//untested check
	     printf("Error::Cluster::ComputeHMBIGradient() cannot find which atom %i is symmetrical\n",
		    Monomers[indexB].GetAtom(j).GetAtomIndex());
	     exit(0);
	   }
	   
	   for(int l=0;l<UniqueAtoms;l++){
	     if(gradkey[l]==SymGlobalIndex){
               for(int xyz=0;xyz<3;xyz++){
	         // printf("MonB Added to atom %i from %i, %f %f %f \n",gradkey[l], SymAtom,
	         //  DimerImages[i].GetQMGradient()[3*k],DimerImages[i].GetQMGradient()[3*k+1],DimerImages[i].GetQMGradient()[3*k+2]);
	       
	         // printf("    Monomer grad: %f %f %f\n",
	         //   Monomers[SymB].GetQMGradient()[3*SymAtom],Monomers[SymB].GetQMGradient()[3*SymAtom+1],Monomers[SymB].GetQMGradient()[3*SymAtom+2]);
	         fflush(stdout); 
	         //added the dimer contribution to the pair-wise interactions
	         grad[3*l] += 0.5*damp*symfac*Rot(0,xyz)*(DimerImages[i].GetQMGradient()[3*k+xyz]-scafac*DimerImages[i].GetMMGradient()[3*k+xyz]);
	         grad[3*l+1] += 0.5*damp*symfac*Rot(1,xyz)*(DimerImages[i].GetQMGradient()[3*k+xyz]-scafac*DimerImages[i].GetMMGradient()[3*k+xyz]);
	         grad[3*l+2] += 0.5*damp*symfac*Rot(2,xyz)*(DimerImages[i].GetQMGradient()[3*k+xyz]-scafac*DimerImages[i].GetMMGradient()[3*k+xyz]);
	       
	       /*
		 if(l==2){
		   printf("MonB: Added to atom %i from atom %i, %f %f %f \n",gradkey[l],Monomers[indexB_ref].GetAtom(j).GetGlobalIndex(),
			  DimerImages[i].GetQMGradient()[3*k],DimerImages[i].GetQMGradient()[3*k+1],DimerImages[i].GetQMGradient()[3*k+2]);
		   printf("damp = %f\n",damp);
		   printf("symfac = %i\n",symfac);
		   printf("grad = %f\n\n",grad[3*l]);
		   fflush(stdout);
		 }
	       */

	       //Subracting the monomer contribution.
	       //If counterpoise corrected, the monomer contribution is already subtrated for the QM.
	       //Since symmetry is in use, the grad of the symmetry atoms will be used
	       if(!Params::Parameters().DoCounterpoise()){
		 if(UseMonSym){
		   grad[3*l] -=  0.5*damp*symfac*Rot(0,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
		   grad[3*l+1] -=  0.5*damp*symfac*Rot(1,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
		   grad[3*l+2] -=  0.5*damp*symfac*Rot(2,xyz)*Monomers[SymB].GetQMGradient()[3*SymAtom+xyz];
		 }else{
		   grad[3*l] -=  0.5*damp*symfac*Rot(0,xyz)*Monomers[indexB_ref].GetQMGradient()[3*j+xyz];
		   grad[3*l+1] -=  0.5*damp*symfac*Rot(1,xyz)*Monomers[indexB_ref].GetQMGradient()[3*j+xyz];
		   grad[3*l+2] -=  0.5*damp*symfac*Rot(2,xyz)*Monomers[indexB_ref].GetQMGradient()[3*j+xyz];	
		 }
	       }

	       /*
	       if(l==2){
		 printf("     Monomer grad: %f %f %f\n",
			Monomers[indexA].GetQMGradient()[3*j],Monomers[indexA].GetQMGradient()[3*j+1],Monomers[indexA].GetQMGradient()[3*j+2]);
		 printf("grad = %f\n\n",grad[3*l]);
	       }
	       */
	       
	       
	       //For the monomer MM of the two-body interaction, a negative is being subtrated therefore it is positive.
	       //If the MM uses qchem and counterpoise corected, the two the monomer contribution is already subtrated
	       if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3){
		 if(UseMonSym){
		   grad[3*l] +=  0.5*damp*symfac*scafac*Rot(0,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
		   grad[3*l+1] += 0.5*damp*symfac*scafac*Rot(1,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
		   grad[3*l+2] += 0.5*damp*symfac*scafac*Rot(2,xyz)*Monomers[SymB].GetMMGradient()[3*SymAtom+xyz];
		 }else{
		   grad[3*l] +=  0.5*damp*symfac*scafac*Rot(0,xyz)*Monomers[indexB_ref].GetMMGradient()[3*j+xyz];
		   grad[3*l+1] += 0.5*damp*symfac*scafac*Rot(1,xyz)*Monomers[indexB_ref].GetMMGradient()[3*j+xyz];
		   grad[3*l+2] += 0.5*damp*symfac*scafac*Rot(2,xyz)*Monomers[indexB_ref].GetMMGradient()[3*j+xyz];
		 }
	       }
	       
	       //gradient for the damping function
	       grad[3*l] += 0.5*Rot(0,xyz)*DampGrad[3*k+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       grad[3*l+1] += 0.5*Rot(1,xyz)*DampGrad[3*k+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       grad[3*l+2] += 0.5*Rot(2,xyz)*DampGrad[3*k+xyz]*symfac*(DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
	       
	       /*
	       if(l==2){
		 printf("DampGrad = %f\n",DampGrad[3*k]);
		 printf("QMInt = %f\n",DimerImages[i].GetQMIntEnergy());
		 printf("grad = %f\n\n",grad[3*l]);
		 fflush(stdout);
	       }
	       */
              }//end of loop over xyz
	     }
	   } 
	 }
      }
     }

     //grad.PrintGradient("periodic two body");
    
    if(Params::Parameters().FreezeNuclearCoordinates())
      grad.Set();

    if(Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() == 1
       && !Params::Parameters().DoFiniteDifferenceFreqs() ) {
      ComputeLatticeParamGrad(grad);
    }


   }


  //zeroing out any atoms that are frozen by symmetry and combining x,y,z gradients that are linked by symmetry
  for(int i=1;i<=NMon;i++){
    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
      for(int k = 0;k<UniqueAtoms;k++){
	if(gradkey[k] == Monomers[i].GetAtom(j).GetGlobalIndex()){
	  
	  //yoni print statement
	  //printf("atom %i Frozenatom = %i LockY = %i SetSignY = %i LockZ = %i SetSignZ = %i\n",
	  //	 Monomers[i].GetAtom(j).GetGlobalIndex(),Monomers[i].GetAtom(j).IsAtomFrozen(),
	  //	 Monomers[i].GetAtom(j).IsYLocked(),Monomers[i].GetAtom(j).ChangeYSign(),
	  //	 Monomers[i].GetAtom(j).IsZLocked(),Monomers[i].GetAtom(j).ChangeZSign());
    


	  //if atom is frozen zero out grad
	  if(Monomers[i].GetAtom(j).IsAtomFrozen()){
	    grad[3*k] = 0.0;
	    grad[3*k+1] = 0.0;
	    grad[3*k+2] = 0.0;
	  }

	  //change sign of X or y to preserve symmetry
	  
	  if(Monomers[i].GetAtom(j).ChangeXSign())
	    grad[3*k] = -grad[3*k];
	  if(Monomers[i].GetAtom(j).ChangeYSign())
	    grad[3*k+1] = -grad[3*k+1];



	  //grad.PrintGradient("grad before combinding\n");	  


          //if coordinates are locked under symmetry, make it zero
	  //or if coordinates of atoms are locked under symmetry, add them together
      

	  //SPECIAL CASE: Hexagonal cell with an atom is a fixed fractional x but with a cartesian coordinate with a y degree of freedom.
	  //The updated x cartesian coordinate (which is not a degree of freedom) is shifted so y can be optimized but the x fractional coordinate
	  // is a constant.
	  //Only tested for R-3c phase IV CO2
	  //delta_x = 1/tan(gamma) * delta_y where gamma = 120 degrees
	  double tolerance = Params::Parameters().GetSymmetryTolerance() ;
	  if((Monomers[i].GetAtom(j).FreezeX()||Monomers[i].GetAtom(j).IsXLocked()==1) && !(Monomers[i].GetAtom(j).FreezeY()||Monomers[i].GetAtom(j).IsYLocked())) {
	    if(Params::Parameters().IsPeriodic()) {
              if((UnitCellAngles[2]-120)<tolerance) {
	        grad[3*k+1] += grad[3*k]/tan(UnitCellAngles[2]*DegreesToRadians);
	        grad[3*k] = 0.0;
              }
            }
	  }
          //X frozen or locked
          else if(Monomers[i].GetAtom(j).FreezeX()){
           grad[3*k] = 0.0;
          }
	  else if(Monomers[i].GetAtom(j).IsXLocked() == 1){
	    grad[3*k+1] += grad[3*k];
	    grad[3*k] = 0.0;
	  }
	  else if(Monomers[i].GetAtom(j).IsXLocked() == 2){
	    grad[3*k+2] += grad[3*k];
	    grad[3*k] = 0.0;
	  }


	  
          //Y frozen or locked
          if(Monomers[i].GetAtom(j).FreezeY()){
	    grad[3*k+1] = 0.0;
          }
	  else if(Monomers[i].GetAtom(j).IsYLocked()){
	    grad[3*k+2] += grad[3*k+1];
	    grad[3*k+1] = 0.0;
	  }

	  /*
	  else if(Monomers[i].GetAtom(j).IsYLocked()){
	    grad[3*k] += grad[3*k+1];
	    grad[3*k+1] = 0.0;
	  }
	  */

          //Z frozen or locked
          if(Monomers[i].GetAtom(j).FreezeZ()){
	    grad[3*k+2] = 0.0;
	  }

	  /*else if(Monomers[i].GetAtom(j).IsZLocked() == 1){
	    grad[3*k] += grad[3*k+2];
	    grad[3*k+2] = 0.0;
	  }
	  else if(Monomers[i].GetAtom(j).IsZLocked() == 2){
	    grad[3*k+1] += grad[3*k+2];
	    grad[3*k+2] = 0.0;
	  }
	  */
	}
      }
    }
  }

 } //End If statement which skips the monomer and dimer contributions if only using a full cluster 

  //printf("Calling debug lattice gradient routine\n");
  //DebugLatticeGradient();
  //cell params are frozen
  //Note: have not yet coded in the check to see which unit cell params are frozen
  if(Params::Parameters().UseFullQMOnly() ){

   if(Params::Parameters().FreezeUnitCellParams()){
      for(int i = 0; i < 9; i++){
        grad[3*UniqueAtoms + i] = 0;
      }
    }
 
    if(Params::Parameters().FreezeNuclearCoordinates()){
      for(int i = 0; i < grad.GetLength() - 9; i++){
        grad[i] = 0;
      }
    }
  }

  PrintGradient("HMBI Gradient",grad);

  //grad.Print("Analytical");
   // delete grad;
  // delete [] key;
  //yoni::RMS and Max of parts of Gradient
  double RMS_nuclear = grad[0]*grad[0];
  int Max_nuclear = 0;
  double RMS_length;
  int Max_length;
  double RMS_angle;
  int Max_angle;

  //nuclear gradient

  for(int i=0;i<3*UniqueAtoms;i++){
    RMS_nuclear += grad[i]*grad[i];
    if(fabs(grad[i]) > fabs(grad[Max_nuclear]))
      Max_nuclear = i;
  }
  RMS_nuclear = sqrt(RMS_nuclear/UniqueAtoms);

  printf("nuclear grad RMS = %f Max = %f\n",RMS_nuclear, grad[Max_nuclear]);
  //lattice gradient
  if(((Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2) || Params::Parameters().UseFullQMOnly() )&& !Params::Parameters().DoFiniteDifferenceFreqs()) {
    RMS_length = grad[3*UniqueAtoms]*grad[3*UniqueAtoms];
    Max_length = 3*UniqueAtoms;
    RMS_angle = grad[3*UniqueAtoms+3]*grad[3*UniqueAtoms+3];
    Max_angle = 3*UniqueAtoms+3;

    for(int i=3*UniqueAtoms+1;i<3*UniqueAtoms+3;i++){
      RMS_length += grad[i]*grad[i]; 
      RMS_angle += grad[i+3]*grad[i+3];
      if(fabs(grad[i]) > fabs(grad[Max_length]))
	Max_length = i;
      if(fabs(grad[i+3]) > fabs(grad[Max_angle]))
	Max_angle = i+3;
    }
    RMS_length = sqrt(RMS_length/3);
    RMS_angle = sqrt(RMS_angle/3);

    printf("length grad RMS = %f Max = %f\n",RMS_length,grad[Max_length]);
    printf("angle grad RMS = %f Max = %f\n",RMS_angle,grad[Max_angle]);
   
  }

  //return FullGrad;

   return grad;
}

 // Prints Q-Chem-style $molecule section
 void Cluster::PrintQChemCartesian(FILE *outfile) {
   fprintf(outfile,"$molecule\n%d %d\n",charges[0],spins[0]);
   for (int i=1;i<=NMon;i++) {
     Monomers[i].PrintMonomerCartesian(outfile);
   }
   fprintf(outfile,"$end\n");

 }

 // Prints only the xyz coordinates
 void Cluster::PrintCartesianCoordinates(FILE *outfile) {
   for (int i=1;i<=NMon;i++) {
     Monomers[i].PrintMonomerCartesian(outfile);
   }

 }

 void Cluster::PrintTinkerCartesian(FILE *outfile) {

   // Print title line
   fprintf(outfile,"%d  Full Cluster\n",GetTotalNumberOfAtoms()); 
   int shift = 0;

   // apply shifts to ensure proper indexing of atom numbers/connectivity
   for (int i=1;i<=NMon;i++) {
     Monomers[i].PrintTinkerCartesian(shift,false,outfile);
     shift += Monomers[i].GetNumberOfAtoms();
   }
 }

 //Print coordinates that are not the same as the cluster
 void Cluster::PrintTinkerCartesian(Vector coords,FILE *outfile) {


   if(coords.GetLength() != 3*GetTotalNumberOfAtoms()){
     printf("Error::Cluster::PrintTinkerCartesian() Length of inputed coords is not the same size as the coordinates of the cluster\n coords = %i 3*Natoms = %i"
	    ,coords.GetLength(),3*GetTotalNumberOfAtoms());
     exit(0);
   }

   //Print title line
   fprintf(outfile, "%d Full Cluster\n",GetTotalNumberOfAtoms());

   //coordinates counter to ensure monomers have the correct coordinates
   int icoords = 0;

   //index to apply shift to ensure proper indexing of atom numbers/connectivity
   int shift = 0;

   //writing coordinates for each monomer with the proper indexes for atom  numbers/connectivity
   for (int i=1;i<=NMon;i++) {

     //Get the coordinates for each monomer
     Vector MonomerCoords(3*Monomers[i].GetNumberOfAtoms());
     for(int j=0;j<3*Monomers[i].GetNumberOfAtoms();j++){
       MonomerCoords[j] = coords[icoords];
       icoords++;
     }

     //MonomerCoords.PrintGradient("monomer coords");

     Monomers[i].PrintTinkerCartesian(MonomerCoords,shift,false,outfile);
     shift += Monomers[i].GetNumberOfAtoms();   
   }

 }


// Prints XYZ file cartesian coordinates, including header lines
void Cluster::PrintXYZCartesian(FILE *outfile) {
  // Print tile line
  fprintf(outfile,"%d\nFull Cluster\n",GetTotalNumberOfAtoms()); 
  for (int i=1;i<=NMon;i++) {
     Monomers[i].PrintMonomerCartesian(outfile);
  }

}

 /*
 void Cluster::SaveCurrentGeomToDisk(FILE *xyzfile) {

   if (

   fprintf(outfile,"%d\nCycle %d, Energy = %15.9f\n",GetTotalNumberOfAtoms(),
	   Parameters.GetOptCycle(),Energy_HMBI); 
   for (int i=1;i<=NMon;i++) {
     Monomers[i].PrintMonomerCartesian(outfile);
   }
 }
  */

// Prints an entirely new input file, e.g. in case you changed geometry
void Cluster::PrintInputFile(FILE *outfile) {
 
 bool printline = true;
 for(int j = 0;j<NumberOfLinesInInput;j++){

   if(InputFile[j].substr(0,9) == "$molecule"){
     printline = false;
     fprintf(outfile,"$molecule\n%d %d\n",charges[0], spins[0]);
     for (int i=1;i<=NMon;i++) {
       if (Params::Parameters().GetMMType() == 1) {
         fprintf(outfile,"--\n%d %d\n",Monomers[i].GetChargeState(), 
	       Monomers[i].GetSpinState() );
         Monomers[i].PrintTinkerCartesian(0,false,outfile);
       }
       else if (Params::Parameters().GetMMType() == 2) {
         fprintf(outfile,"--\n%d %d %f\n",Monomers[i].GetChargeState(), 
	         Monomers[i].GetSpinState(), 
	         Monomers[i].GetIonizationPotential());
         Monomers[i].PrintMonomerCartesian(outfile);
       }
       else {
         fprintf(outfile,"--\n%d %d\n",Monomers[i].GetChargeState(), 
	         Monomers[i].GetSpinState());
         Monomers[i].PrintMonomerCartesian(outfile);
       }

     }  

   }



   if(InputFile[j].substr(0,10) == "$unit_cell" && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
     printline = false;
     fprintf(outfile,"$unit_cell\n%f %f %f\n%f %f %f\n",UnitCellAxes[0],
       UnitCellAxes[1], UnitCellAxes[2], UnitCellAngles[0],
       UnitCellAngles[1], UnitCellAngles[2]);

       }
   if(InputFile[j].substr(0,4) == "$end" && printline == false){
      printline = true;
   }
   if(printline)
     fprintf(outfile,"%s\n",InputFile[j].c_str());

 }

 fflush(stdout);
}

// Prints an entirely new input file, e.g. in case you changed geometry
void Cluster::OldPrintInputFile(FILE *outfile) {

   // Print $comment section
   fprintf(outfile,"$comment\n%s$end\n\n",Title.c_str());

   // Print $hmbi section
   fprintf(outfile,"%s\n", Params::Parameters().GetHMBIRem().c_str() );

   // Print $molecule section
   fprintf(outfile,"$molecule\n%d %d\n",charges[0], spins[0]);

   for (int i=1;i<=NMon;i++) {
     if (Params::Parameters().GetMMType() == 1) {
       fprintf(outfile,"--\n%d %d\n",Monomers[i].GetChargeState(), 
	     Monomers[i].GetSpinState() );
       Monomers[i].PrintTinkerCartesian(0,false,outfile);
     }
     else if (Params::Parameters().GetMMType() == 2) {
       fprintf(outfile,"--\n%d %d %f\n",Monomers[i].GetChargeState(), 
	       Monomers[i].GetSpinState(), 
	       Monomers[i].GetIonizationPotential());
       Monomers[i].PrintMonomerCartesian(outfile);
     }
     else {
       fprintf(outfile,"--\n%d %d\n",Monomers[i].GetChargeState(), 
	       Monomers[i].GetSpinState());
       Monomers[i].PrintMonomerCartesian(outfile);
     }

   }
   fprintf(outfile,"$end\n\n");

   // Print $qchem section - erase "$rem" and replace it with "$qchem"
   string tmp = Params::Parameters().GetQChemRem();
   tmp.erase(0,5);
   fprintf(outfile,"$qchem%s", tmp.c_str() );

   // Print MM section
   if (Params::Parameters().GetMMType() == 1) // Tinker
     fprintf(outfile,"$tinker%s$end\n\n", Params::Parameters().GetTinkerRem().c_str() );
   else if (Params::Parameters().GetMMType() == 2) // AIFF
     fprintf(outfile,"$aiff%s$end\n\n", Params::Parameters().GetAIFFRem().c_str() );
   else if (Params::Parameters().GetMMType() == 3) {// Qchem2;
     string tmp = Params::Parameters().GetQChemRem2();
     tmp.erase(0,5); //erase "$rem" & replace it with "$qchem2"
     fprintf(outfile,"$qchem2%s", tmp.c_str() );
   }

   if ( Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()) {
     fprintf(outfile,"$unit_cell\n%f %f %f\n%f %f %f\n$end\n",UnitCellAxes[0],
	     UnitCellAxes[1], UnitCellAxes[2], UnitCellAngles[0],
	     UnitCellAngles[1], UnitCellAngles[2]);
   }

   fprintf(outfile,"\n");

 }

void Cluster::PrintFractionalCoordinates(FILE *outfile, bool includeLatParams){
  //JLM changed for convience. Nothing wrong with the function
  if (includeLatParams)
    fprintf(outfile,"fractional coordinates\n\n");

    for(int i = 1;i <= NMon;i++)
      for(int j = 0; j < Monomers[i].GetNumberOfAtoms(); j++)
	fprintf(outfile, "%-2s  %10.6f  %10.6f  %10.6f\n", 
		Monomers[i].GetSymbol(j).c_str(),
		Monomers[i].GetAtom(j).GetFractionalPosition()[0],
		Monomers[i].GetAtom(j).GetFractionalPosition()[1],
		Monomers[i].GetAtom(j).GetFractionalPosition()[2]);
    if (includeLatParams){
      fprintf(outfile,"\nlattice length  %10.6f %10.6f %10.6f\n",
	    UnitCellAxes[0],UnitCellAxes[1],UnitCellAxes[2]);
      fprintf(outfile,"lattice angles %10.6f %10.6f %10.6f\n",
	    UnitCellAngles[0],UnitCellAngles[1],UnitCellAngles[2]);
    }
}


// Shifts distance of each molecule relative to the cluster
// center of mass, controlled by "factor".  factor=1.0 means
// no change, factor < 1 shrinks the cluster, and factor > 1
// expands it.
void Cluster::AdjustIntermolecularSpacing(double factor) {
  
  printf("Expanding geometry by factor of %f\n",factor);
  
  // Create temp space for holding shifted centers of Mass.
  Vector old_mon_com(3), new_mon_com(3);
  
  for (int i=1;i<=NMon;i++) {// loop over monomers
    
    for (int dim=0;dim<3;dim++) {// loop over xyz
      old_mon_com[dim] = Monomers[i].GetCenterOfMass(dim);
      double diff = old_mon_com[dim] - CenterOfMass[dim];
      diff *= factor;
      new_mon_com[dim] = CenterOfMass[dim] + diff;
    }
    
    double old_dist = sqrt( pow(old_mon_com[0] - CenterOfMass[0],2)
			    + pow(old_mon_com[1] - CenterOfMass[1],2)
			    + pow(old_mon_com[2] - CenterOfMass[2],2) );
    
    double new_dist = sqrt( pow(new_mon_com[0] - CenterOfMass[0],2)
			    + pow(new_mon_com[1] - CenterOfMass[1],2)
			    + pow(new_mon_com[2] - CenterOfMass[2],2) );
    
    /*
     printf("\nMonomer %d\n",i);
     printf("old_com: (%f,%f,%f)\n",old_mon_com[0],old_mon_com[1],
     old_mon_com[2]);
     printf("new_com: (%f,%f,%f)\n",new_mon_com[0],new_mon_com[1],
     new_mon_com[2]);
     
     printf("Distance from Cluster COM to Mon COM.  Old: %f   New: %f\n",
     old_dist,new_dist);
    */
    Monomers[i].Translate(new_mon_com);

  }

  FILE *geom;
  string filename = "new_geom.xyz";
  if ((geom = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::AdjustIntermolecularSpacing : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }
  PrintXYZCartesian(geom);
  fclose(geom);
  
  printf("New geometry written to '%s'\n",filename.c_str());
  
   
  FILE *input;
  string input_file = "new_geom.in";
  if ((input = fopen(input_file.c_str(),"w"))==NULL) {
    printf("Cluster::AdjustIntermolecularSpacing : Cannot open file '%s'\n",
	   input_file.c_str());
    exit(1);
  }
  
  PrintInputFile(input);
  printf("\nNew input file written to '%s'\n",input_file.c_str());
  
  fclose(input);
  
}

void Cluster::AdjustCoord(int mode, Matrix RModeVector) { //Watit

	for (int i=1;i<=NMon;i++) {
		int Natoms = Monomers[i].GetNumberOfAtoms(); 
		Vector MonomerRModeVector(Natoms*3);
		for (int j=0;j<Natoms;j++) {
			for (int dim=0;dim<3;dim++) {
				MonomerRModeVector[dim+j*3] = RModeVector(j+(i-1)*Natoms,dim);
			}
		}
		Monomers[i].ShiftCoords(MonomerRModeVector);
	}

  	FILE *input;
  	char nametag[20];
  	sprintf(nametag,"fix_img_geom_%d.in",mode);
  	string input_file = nametag;
  	if ((input = fopen(input_file.c_str(),"w"))==NULL) {
    		printf("Cluster::AdjustCoord : Cannot open file '%s'\n",input_file.c_str());
    		exit(1);
  	}
  
  	PrintInputFile(input);
  	printf("\nNew input file written to '%s'\n",input_file.c_str());
  	fclose(input);
        exit(0);
}

void Cluster::UpdateTrajectoryFile(int step, bool new_trajectory) {
  
  printf("Updating trajectory file\n");
  ofstream traj;
  string filename = "traj.xyz";
  int Ntot = GetTotalNumberOfAtoms();
  int periodic_factor = 1;
  if ( Params::Parameters().IsPeriodic()) 
    periodic_factor = 8;
  
  // If starting new file
  if (new_trajectory) {
    traj.open(filename.c_str());
    traj << Ntot*periodic_factor << endl << "Initial Geometry, E = " << Energy_HMBI << endl;
    //fprintf(traj,"%d\nInitial Geometry, E = %15.9f\n",
    //	    Ntot*periodic_factor,Energy_HMBI);
  }
  
  // If appending another geometry to an existing trajectory
  else {
    traj.open(filename.c_str(),ios::app);
    traj << Ntot*periodic_factor << endl << "Step " << step << ", E = " << Energy_HMBI << endl;
    //fprintf(traj,"%d\nStep %d, E = %15.9f\n",
    //	    Ntot*periodic_factor,step,Energy_HMBI);
  }
  
  traj.precision(6);

  // Write the actual geometry
  if ( Params::Parameters().IsPeriodic()) {
    for (int k1=0;k1<2;k1++)
      for (int k2=0;k2<2;k2++)
	for (int k3=0;k3<2;k3++) {
	  double shift_x, shift_y, shift_z;
	  shift_x = k1*unit_cell[0][0] + k2*unit_cell[1][0] + k3*unit_cell[2][0];
	  shift_y = k1*unit_cell[0][1] + k2*unit_cell[1][1] + k3*unit_cell[2][1];
	  shift_z = k1*unit_cell[0][2] + k2*unit_cell[1][2] + k3*unit_cell[2][2];
	  //printf("cell < %d, %d, %d >.  Shift = %f, %f, %f\n",k1,k2,k3,shift_x,shift_y,shift_z);
	  
	  fflush(stdout);
	  
	  for (int i=0;i<Ntot;i++) {
	    double x = AtomicCoordinates[3*i] + shift_x;
	    double y = AtomicCoordinates[3*i+1] + shift_y;
	     double z = AtomicCoordinates[3*i+2] + shift_z;
	     traj << AtomicSymbols[i] << fixed << "   " << x << "   " << y << "   " << z << endl;
	  }
	}
  }
  else {
    for (int i=0;i<Ntot;i++) {
      double x = AtomicCoordinates[3*i];
      double y = AtomicCoordinates[3*i+1];
      double z = AtomicCoordinates[3*i+2];
       traj << AtomicSymbols[i] << fixed << "   " << x << "   " << y << "   " << z << endl;
    }
  }
}

void Cluster::ComputeDistanceMatrix() {
  
  for (int i=1;i<=NMon;i++) {
    printf("Monomer %d intramolecule distances\n",i);
    Monomers[i].ComputeIntramolecularDistances();
  }
  
  for (int i=1;i<=NDim;i++) {
    int m1 = Dimers[i].GetIndexA();
    int m2 = Dimers[i].GetIndexB();
    printf("Dimer (%d, %d) intermolecular distances\n",m1,m2);
    Dimers[i].ComputeIntermolecularDistances();
   }
  
  
}


// Print out the geometry containing N1xN2xN3 cells
void Cluster::WriteCrystalXYZFile(int N1, int N2, int N3) {
  
  FILE *xyz;
  string filename = "crystal.xyz";
   if ((xyz = fopen(filename.c_str(),"w"))==NULL) {
     printf("Cluster::WriteCrystalXYZFile : Cannot open file '%s'\n",
	    filename.c_str());
     exit(1);
   }
   
   
   // Count number of atoms, and print XYZ file header line
   int Natoms_cell = GetTotalNumberOfAtoms();
   int Natoms = Natoms_cell * N1 * N2 * N3;

   fprintf(xyz,"%d\n%d x %d x %d cell\n",Natoms,N1,N2,N3);

   for (int x=0;x<N1;x++) {
     for (int y=0;y<N2;y++) {
       for (int z=0;z<N3;z++) {

	 // Create copy of Monomers that we can translate as needed;
	 Monomer* ImageMonomers = new Monomer[NMon+1];
	 for (int i=1;i<=NMon;i++) {
	   ImageMonomers[i] = Monomers[i];
	 }
	
	 // Determine the shift from the central cell to the image cell
	 Vector shift(3);
	 for (int i=0;i<3;i++) {
	   shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	     z*unit_cell[2][i];
	 }

	 // Loop over monomers, shift them as needed, & print their
	 // coordinates
	 for (int imon=1;imon<=NMon;imon++) {
	   Vector new_com(3);
	   new_com = ImageMonomers[imon].GetCenterOfMass();

	   // Add the shift and translate the monomer
	   new_com += shift;

	   ImageMonomers[imon].Translate(new_com);
	   ImageMonomers[imon].PrintMonomerCartesian(xyz);
	 } // end loop over monomers
       } // end loop over z
     } // end loop over y
   } // end loop over x

   fclose(xyz);

 }


void Cluster::CreatePeriodicImageDimerList() {

   //reset all periodic symmetry factors of all dimers.
   for(int i=1;i<=NDim;i++){
     Dimers[i].SetPeriodicSymmetryFactor(0);
   }

   // Identify how far we have to go along each unit cell direction to 
   // stay within the cutoff.  Add 1 extra image cell to each, for good measure.
   double r_cutoff = Params::Parameters().GetLocalCutoff(0); // c0 cutoff
   int Nv[3] = {1,1,1}; // start with 1 image cell in each direction
   for (int i=0;i<3;i++) {
     double dist = 0;
     while (dist  < r_cutoff) {
       dist +=  unit_cell[i].Norm();
       Nv[i] += 1;
     }
     //printf("Nv[%d] = %d, dist = %f\n",i,Nv[i],dist);
   }

   // Prepare output xyz file for visualizing
   // Open the input file for writing
   FILE *xyz;
   string filename = "pbc.xyz";
   if ((xyz = fopen(filename.c_str(),"w"))==NULL) {
     printf("Cluster::CreatePeriodicImageDimerList() : Cannot open file '%s'\n",
	    filename.c_str());
     exit(1);
   }

   fprintf(xyz,"XXX\n\n"); // we set the number of atoms later

   // print central unit cell coords
   for (int i=1;i<=NMon;i++) 
     Monomers[i].PrintMonomerCartesian(xyz);
   int active_atoms = GetTotalNumberOfAtoms();

   // Create dynamic 6*N array, accesse as ImageList[i][j].  
   // The j-index is 6 items long: the Lattice Vector (3 items), the central
   // unit-cell monomer, the reference monomer that has been replicated,
   // and the multiplicity due to symmetry of that vector.
   // The other dimension, N grows dynamically to accommodate as many entries
   // as are needed to store the list that later generates all Images dimers.
   vector < vector<int> > ImageList(1,vector<int>(6,0));

   /*
     Start generating image cells.  Use a 2-step algorithm:  

     Pass 1: identifies the dimers we need to include between the
     central unit cell and its replicates.  These results are stored in
     ImageList.

     Pass 2: uses ImageList to construct the necessary dimers.  

     This is probably a slow algorithm, especially since it dynamically
     grows the ImageList.  But, Pass 1 lets us work with fairly small
     amounts of memory and to check for symmetry, etc, and to find how many
     dimers we actually need.  Then, we do the real part in Pass 2 with the 
     pruned list of known size.

     Pass 3: Further Prune number of dimers using space symmetry

   */

   bool UseSymmetry = Params::Parameters().UseCrystalSymmetry();

   // Pass 1
   bool skip; // flaging some.
   int keepit = 0;
   // Loop over the image cells, in both positive & negative directions
   for (int x=-Nv[0];x<=Nv[0];x++) 
     for (int y=-Nv[1];y<=Nv[1];y++)
       for (int z=-Nv[2];z<=Nv[2];z++) {
	 skip = false;

	 if (x==0 && y==0 && z==0) {
	   //printf("Skipping central unit cell, x=%d, y=%d, z=%d\n",x,y,z);
	   skip = true;
	 }

	 if (!skip) {
	   // Create copy of Monomers that we can translate as needed;
	   Monomer* ImageMonomers = new Monomer[NMon+1];
	   for (int i=1;i<=NMon;i++) {
	     ImageMonomers[i] = Monomers[i];

	     //printf("ImageMonomer %d before shifting:\n",i);
	     //ImageMonomers[i].PrintAll();
	   }


	   // Determine the shift from the central cell to the image cell
	   Vector shift(3);
	   for (int i=0;i<3;i++) {
	     shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	       z*unit_cell[2][i];
	   }

	   // Now translate the monomers
	   for (int imon=1;imon<=NMon;imon++) {

	     bool IsThisMonomerLocal = false;

	     Vector new_com(3);
	     new_com = ImageMonomers[imon].GetCenterOfMass();

	     //printf("Before shift\n");
	     //printf("Monomer %d\n",imon);
	     //Monomers[imon].PrintQChemCartesian();
	     //printf("Image Monomer %d\n",imon);
	     //ImageMonomers[imon].PrintQChemCartesian();
	     //printf("indexing: M: %d    IM: %d\n",Monomers[imon].GetIndex(),
	     //   ImageMonomers[imon].GetIndex());

	     // Add the shift and translate the monomer
	     new_com += shift;

	     ImageMonomers[imon].Translate(new_com);


	     //ImageMonomers[imon].SetIndex(153);
	     //printf("Monomer %d\n",imon);
	     //Monomers[imon].PrintQChemCartesian();
	     //printf("Image Monomer %d\n",imon);
	     //ImageMonomers[imon].PrintQChemCartesian();
	     //printf("indexing: M: %d    IM: %d\n",Monomers[imon].GetIndex(),
	     //	   ImageMonomers[imon].GetIndex());


	     // Now pair this image monomer with each monomer in central unit 
	     // cell and test distance relative to the cutoff
	     for (int jmon=1;jmon<=NMon;jmon++) {

	       Dimer Tmp;
	       //printf("jmon = %d, imon = %d\n",jmon,imon);
	       Tmp.Initialize(Monomers[jmon],ImageMonomers[imon]);
	       //printf("dimer from %d and image %d\n",jmon,imon);
	       //Tmp.PrintQChemCartesian();
	       if ( Tmp.GetDimerSeparation() < r_cutoff ) {
		 if ( Params::Parameters().PrintLevel() > 3 ) {
		   printf("K_vec = (%d,%d,%d)\n",x,y,z);
		   printf("(%d*,%d) Separation = %f. Keeping it.\n",
			  imon,jmon,Tmp.GetDimerSeparation() );
		 }
		 IsThisMonomerLocal = true;

		 bool AlreadyExistsBySymmetry = false;
		 if (UseSymmetry) {
		   // Test for symmetrically equivalent pair in ImageList
		   // Currently only test for translational symmetry.
		   for (int q=0;q<ImageList.size();q++) {
		     int qx = ImageList[q][0];
		     int qy = ImageList[q][1];
		     int qz = ImageList[q][2];
		     int qj = ImageList[q][3];
		     int qi = ImageList[q][4];

		     // Check translational symmetry: K(x,y,z) = -K, and swap
		     // order of imon/jmon.
		     if ( qx==-x && qy==-y && qz==-z && qi==jmon && qj==imon) {
		       if ( Params::Parameters().PrintLevel() )
			 printf("(%d,%d,%d) dimer (%d,%d) is equivalent to existing (%d,%d,%d) dimer (%d,%d) by symmetry\n",x,y,z,jmon,imon,qx,qy,qz,qj,qi);

		       // if match found, increment symmetry factor by 1
		       ImageList[q][5] += 1;
		       AlreadyExistsBySymmetry = true;
		     }

		   }
		 }

		 if (!AlreadyExistsBySymmetry) {
		   // Store info to later regenerate this list
		   int tmp[] = {x,y,z,jmon,imon,1};
		   vector<int> tmpvec (tmp,tmp+6);
		   if (keepit==0)
		     ImageList[0] = tmpvec; // first time, modify existing row
		   else {
		     ImageList.push_back(tmpvec); // add row to list
		   }
		   keepit++;
		 }
	       }
	     }

	     // Print out the coordinates.  Substitute based on whether or
	     // or not it is close enough to be "local".  If local, use "N"
	     // atoms.  Else, "Li" atoms.
	     if (IsThisMonomerLocal) {
	       ImageMonomers[imon].PrintMonomerCartesian(xyz);
	       active_atoms += ImageMonomers[imon].GetNumberOfAtoms();
	     }
	     //else
	     //  ImageMonomers[imon].PrintMonomerCartesian(xyz,"Li");

	   }
	   delete [] ImageMonomers;
	 }	  
       }
   //printf("Total number of image dimers = %d\n",keepit);
   if ( Params::Parameters().PrintLevel() > 1) {
     printf("List of saved image dimers\n");
     for (int i=0;i<ImageList.size();i++) {
       printf("%2d: (%d,%d,%d)  Dimer(%d,%d*)  %d\n",i,ImageList[i][0],
	      ImageList[i][1],ImageList[i][2],ImageList[i][3],ImageList[i][4],
	      ImageList[i][5]);
     }
   }

   // Close up geometry file
   fclose(xyz);
   // update number of atoms;
   string cmd = "sed -i -e s/XXX/";
   char count[10];
   sprintf(count,"%d",active_atoms);
   cmd += count;
   cmd += "/g " + filename;
   //printf("Runing: %s\n",cmd.c_str());
   system(cmd.c_str());

   // Pass 2: using symmetry to further weed ImageList
   Dimer* temp_DimerImages;
   temp_DimerImages = new Dimer[keepit+1];

   //tracking # of Dimers if no symmetry was used
   TotalDim_nosym += keepit;

   int Mon_index = NMon+1;
   int index=1;
   for (int idim=1;idim<=keepit;idim++) {

     //Mon_index++;
     // Read data from ImageList
     int x = ImageList[idim-1][0];
     int y = ImageList[idim-1][1];
     int z = ImageList[idim-1][2];
     int real_mon = ImageList[idim-1][3];
     int image_mon = ImageList[idim-1][4];
     int symfac = ImageList[idim-1][5];

     // Create copy of Monomer that we can translate
     Monomer ImageMon;
     ImageMon = Monomers[image_mon];

     // Determine the shift from the central cell to the image cell
     Vector shift(3);
     for (int i=0;i<3;i++) {
       shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	 z*unit_cell[2][i];
     }

     // Translate the monomer    
     Vector new_com(3);
     new_com = ImageMon.GetCenterOfMass();
     new_com += shift;
     ImageMon.Translate(new_com);
     ImageMon.FractionalTranslate(x,y,z);

     // Create the new dimer between real_mon in central cell and image_mon
     // in periodic cell
     ImageMon.SetIndex(Mon_index); // renumber it  
     //ImageMon.SetIndex(idim); // JLM definitely this is the one to change!  
     if (Params::Parameters().PrintLevel() > 3) 
       printf("Monomer %d, K_vec = (%d,%d,%d)\n",Mon_index,x,y,z);
     //changed idim with index
     temp_DimerImages[index].Initialize(Monomers[real_mon],ImageMon);

     // Link image monomer to one in central cell
     temp_DimerImages[index].SetImageCell(x,y,z);
     temp_DimerImages[index].SetReferenceMonomerIndex(image_mon);
     temp_DimerImages[index].SetSymmetryFactor(symfac);

     //Add to dimers that are symmetrical by translation to symmetry list.
     for(int j=2;j<=symfac;j++){
       vector<int> CopySymList = temp_DimerImages[index].GetSymmetryList();
       //vector<bool> CopyMonBList = temp_DimerImages[index].GetMonBList();
       vector<int> CopyKList = temp_DimerImages[index].GetSymmetricalImageCell();
       temp_DimerImages[index].MergeSymmetryList(CopySymList,CopyKList,false);
       
       vector<Vector> CopyAtom_Equivalency = temp_DimerImages[index].GetAtomEquivalency();
       temp_DimerImages[index].MergeEquivalencyList(CopyAtom_Equivalency,CopyAtom_Equivalency[0],false);

       vector<Matrix> CopyRotation = temp_DimerImages[index].GetRotationList();
       Matrix Rot(3,true);
       temp_DimerImages[index].AddToRotationList(Rot,CopyRotation,false,false);
     }
     //Using symmetry to determine if identical to other dimers.
     //First there is a check if any symmetry will be preformed
     // next checking non-image dimers then image imiges.
     // If it is identical, it will be overwrite by next dimer.

     if(!Params::Parameters().UseSpaceSymmetry()){
       //printf("(%d*=%d,%d) Dimer created.\n",image_mon,Mon_index,real_mon);
       //printf("shift = (%d,%d,%d)\n",x,y,z);
       Mon_index++;
       index++;

     }else{
       bool Unique=1;//varible to tell you if dimer identical
		     //to another dimer through symmetry

       if(temp_DimerImages[index].SymmetryCheck(Dimers, NDim+1, true))
	 Unique=0;
       else if(temp_DimerImages[index].SymmetryCheck(temp_DimerImages, index))
	 Unique=0;
       //increasing indexes if not symmetrical to any previous dimers
       if(!Params::Parameters().OldSymmetryDimers()){
         Mon_index++; //JLM: moved this here to keep ordering the same as the no symmetry case
       }
       if(Unique){
	 //if (Params::Parameters().PrintLevel() > 0 && Unique) 
	 //printf("(%d*=%d,%d) Dimer created.\n",image_mon,Mon_index,real_mon);
	 //printf("shift = (%d,%d,%d)\n",x,y,z);
         if(Params::Parameters().OldSymmetryDimers()){
           Mon_index++; //JLM: moved this here to keep ordering the same as the no symmetry case
         }
	 index++;
       }
     }

   }
   
   //Pass 3: Generating Dimers
   
   //for (int i=1;i<=NDim;i++){
   //  printf("d(%i,%i) is symmetrical to %i image dimers\n",
   //	   Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),Dimers[i].GetPeriodicSymmetryFactor());
   //}
   NDim_images = index-1;
   DimerImages = new Dimer[NDim_images+1];

   for (int i=1;i<=NDim_images;i++) { 
     DimerImages[i] = temp_DimerImages[i];
     if(Params::Parameters().PrintLevel() > 0)
       printf("The symmetry factor of d(%i,%i) is %i\n", 
	      DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB(),DimerImages[i].GetSymmetryFactor());
   }
   delete [] temp_DimerImages;

}

void Cluster::UpdateJobStatus(int ijob) {
  //if using monomer symmetry, run only "Unique" monomers, otherwise run all monomers
  int QM_Mon = NMon;
  if(Params::Parameters().UseMonomerSymmetry())
    QM_Mon = UniqueMon;
  
  //if symmetry is exploited by the MM, run only the "Unique" monomers in the MM, otherwise run all monomers. AIFF only
  int MM_Mon = NMon;
  if(Params::Parameters().UseMMSymmetry())
    MM_Mon = UniqueMon;

  // define some helpful points
  int QM_mon_start = 0;
  int QM_dim_start = QM_Mon;
  int MM_start = QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc;  
    
  //MolPro CBS Frequency calculations are determined entirely by 
  //if(Params::Parameters().DoCBS() && Params::Parameters().DoFreq() ){
  if(Params::Parameters().GetQMType() == 2 && Params::Parameters().DoFreq() ){
   QM_dim_start = 0;
   MM_start = 0;
    
    if(Params::Parameters().SingleFileMonomerHess()){
   	 QM_dim_start += QM_Mon;
	 MM_start += QM_Mon;  
   }else{
     for(int iMon = 1; iMon <= NMon; iMon++){
       if((Monomers[iMon].GetSymmetryFactor() != 0) || !Params::Parameters().UseMonomerSymmetry()){
	 QM_dim_start += 6*Monomers[iMon].GetNumberOfAtoms();
	 MM_start += 6*Monomers[iMon].GetNumberOfAtoms();
       }
     }
   }  
   double c0 = Params::Parameters().GetLocalCutoff(0);
   double c1 = Params::Parameters().GetLocalCutoff(1);
   for(int i=1;i<=NDim;i++){
     double separation =  Dimers[i].GetDimerSeparation();
     if(Dimers[i].GetSymmetryFactor() != 0 &&
	(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
       MM_start += 6*Dimers[i].GetNumberOfAtoms();
     }
     
   }
   for(int i=1;i<=NDim_images;i++){
     double separation =  DimerImages[i].GetDimerSeparation();
     if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
       MM_start += 6*DimerImages[i].GetNumberOfAtoms();
     }
   }
  }
  
  if(Params::Parameters().DoCCSDTCorrection()){
    if(Params::Parameters().DoFreq()){
      
      for(int iMon = 1; iMon <= NMon; iMon++){
	if((Monomers[iMon].GetSymmetryFactor() != 0) || !Params::Parameters().UseMonomerSymmetry()){
	  QM_dim_start += 6*Monomers[iMon].GetNumberOfAtoms();
	  MM_start += 6*Monomers[iMon].GetNumberOfAtoms();
	}
      }
      
      double c0 = Params::Parameters().GetLocalCutoff(0);
      double c1 = Params::Parameters().GetLocalCutoff(1);
      for(int i=1;i<=NDim;i++){
	double separation =  Dimers[i].GetDimerSeparation();  
	if(Dimers[i].GetSymmetryFactor() != 0 &&
	   (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	  MM_start += 2*(3*Dimers[i].GetNumberOfAtoms())*(3*Dimers[i].GetNumberOfAtoms()+1);
	}
      }
      
      for(int i=1;i<=NDim_images;i++){
	double separation =  DimerImages[i].GetDimerSeparation();
	if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	  MM_start += 2*(3*DimerImages[i].GetNumberOfAtoms())*(3*DimerImages[i].GetNumberOfAtoms()+1);
	}
      }
    }
    else{
      QM_dim_start += QM_Mon;;
      MM_start += QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc;
    }  


    //CCSD(T) counterpoise dimer forces are handle in seperate jobs
    if(Params::Parameters().DoForces() && Params::Parameters().DoCounterpoise() ){
      double c0 = Params::Parameters().GetLocalCutoff(0);
      double c1 = Params::Parameters().GetLocalCutoff(1);
      for(int i=1;i<=NDim;i++){
	double separation =  Dimers[i].GetDimerSeparation();
	if(Dimers[i].GetSymmetryFactor() != 0 &&
	   (!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0) ){
	  MM_start += 6*Dimers[i].GetNumberOfAtoms();
	}
      }
      for(int i=1;i<=NDim_images;i++){
	double separation =  DimerImages[i].GetDimerSeparation();
	if(!Params::Parameters().DoLocal2BodyTruncation()||separation <= c0){
	  MM_start += 6*DimerImages[i].GetNumberOfAtoms();
	}
      }
    }
    
  }

  if ( Params::Parameters().DoQMBenchmark() )
    MM_start++;
  int N_MM_jobs = 1 + QM_Mon + UniqueDim + NDim_images - NDim_trunc - NDimImages_trunc; // Tinker
  if (Params::Parameters().GetMMType()==2)  
    N_MM_jobs = MM_Mon; // AIFF
  //N_MM_jobs = NMon;
  int MM_end = MM_start + N_MM_jobs;
  
  if (ijob == QM_mon_start && Params::Parameters().UseFullQMOnly())
    printf("\nRunning 1 QM Monomer job:\n");  
  else if (ijob == QM_mon_start)
    printf("\nRunning %d QM Monomer jobs:\n",QM_dim_start);
    //printf("\nRunning %d QM Monomer jobs:\n",QM_Mon);
  
  // insert some spaces in output for pretty printing
  // Break each category into groups of 5 and rows of 50
  
  // QM monomer jobs
  if (ijob < QM_Mon && ijob%5==0 && ijob != QM_dim_start)
    printf(" ");
  
  if (ijob < QM_Mon && ijob%50==0 && ijob > 0 && ijob != QM_dim_start-1 )
    printf("  (%d)\n ",ijob);
  
  // QM dimer jobs
  if (ijob >= QM_dim_start && ijob < MM_start && (ijob-QM_dim_start)%5==0 )
    printf(" ");

  if (ijob >= QM_dim_start && ijob < MM_start && (ijob-QM_dim_start)%50==0 
      && (ijob-QM_dim_start) > 0 && ijob != MM_start-1)
    printf("  (%d)\n ",ijob-QM_dim_start);
  
  /*  
   if (ijob == MM_full)
   printf(" . (1)\n");
  */
  
  
  // MM monomer & dimer jobs
  if (ijob >= MM_start && (ijob-MM_start)%5==0 )
    printf(" ");
  
  if (ijob >= MM_start && (ijob-MM_start)%50==0 && (ijob-MM_start) > 0 &&
      ijob != MM_end-1)
    printf("  (%d)\n ",ijob-MM_start);
  
  
  // print a "." to represent a single job
  printf(".");
  
  if (ijob==QM_dim_start-1) printf("  (%d)\n",ijob+1);
  if (ijob==MM_start-1 && !Params::Parameters().DoQMBenchmark()) 
    printf("  (%d)\n\n",ijob-QM_dim_start+1);
  if (ijob==MM_start-2 && Params::Parameters().DoQMBenchmark() ) 
    printf("  (%d)\n",ijob-QM_dim_start+1);
  if (ijob==MM_start-1 && Params::Parameters().DoQMBenchmark()) 
     printf("  (1)\n\n");
  if (ijob==MM_end-1) printf("  (%d)\n\n",ijob-MM_start+1);
  
  
  else if (ijob==QM_dim_start-1)
     printf("\nRunning %d QM Dimer jobs:\n",MM_start-QM_dim_start);
    //printf("\nRunning %d QM Dimer jobs:\n",UniqueDim+NDim_images - NDim_trunc - NDimImages_trunc);
  else if ( ijob==MM_start-2 && Params::Parameters().DoQMBenchmark() )
    printf("\nRunning benchmark full system QM job:\n ");

  //else if (ijob==MM_full-1)
  //  printf("\nRunning Full system MM job:\n");
  else if (ijob==MM_start-1 && Params::Parameters().GetMMType() !=2 && !Params::Parameters().NeglectManyBody())
    printf("Running %d MM jobs (monomers, dimers, & full system):\n",N_MM_jobs);
  else if (ijob==MM_start-1 && Params::Parameters().GetMMType()==2)
    printf("Computing ab initio Polarizable Force Field Parameters for %d Monomers\n",N_MM_jobs);
  
  
  fflush(stdout);

}

void Cluster::OldUpdateJobStatus(int ijob) {
  
  //if using monomer symmetry, run only "Unique" monomers, otherwise run all monomers
  int QM_Mon = NMon;
  if(Params::Parameters().UseMonomerSymmetry())
    QM_Mon = UniqueMon;

  // define some helpful points
  int QM_mon_start = 0;
  int QM_dim_start = QM_Mon;
  int MM_full = QM_Mon + UniqueDim + NDim_images;
  int MM_start = MM_full + 1;
  int MM_end = MM_start + QM_Mon + UniqueDim + NDim_images;
  
  if (ijob == QM_mon_start)
    printf("\nRunning %d QM Monomer jobs:\n", QM_Mon);
  
  // insert some spaces in output for pretty printing
  // Break each category into groups of 5 and rows of 50
  
  // QM monomer jobs
  if (ijob < QM_Mon && ijob%5==0 && ijob != QM_dim_start)
    printf(" ");
  
  if (ijob < QM_Mon && ijob%50==0 && ijob > 0 && ijob != QM_dim_start-1 )
     printf("  (%d)\n ",ijob);
  
  // QM dimer jobs
  if (ijob >= QM_dim_start && ijob < MM_full && (ijob-QM_dim_start)%5==0 )
    printf(" ");
  
  if (ijob >= QM_dim_start && ijob < MM_full && (ijob-QM_dim_start)%50==0 
      && (ijob-QM_dim_start) > 0 && ijob != MM_full-1)
    printf("  (%d)\n ",ijob-QM_dim_start);
  
  if (ijob == MM_full)
    printf(" . (1)\n");
  
  // MM monomer & dimer jobs
  if (ijob >= MM_start && (ijob-MM_start)%5==0 )
    printf(" ");
  
  if (ijob >= MM_start && (ijob-MM_start)%50==0 && (ijob-MM_start) > 0 &&
      ijob != MM_end-1)
    printf("  (%d)\n ",ijob-MM_start);
  
  
  // print a "." to represent a single job
  if (ijob != NMon+NDim+NDim_images)
    printf(".");
  
  if (ijob==QM_dim_start-1) printf("  (%d)\n",ijob+1);
  if (ijob==MM_full-1) printf("  (%d)\n",ijob-QM_dim_start+1);
  if (ijob==MM_start-1) printf("\n");
  if (ijob==MM_end-1) printf("  (%d)\n\n",ijob-MM_start+1);
  
  
  else if (ijob==QM_dim_start-1)
    printf("\nRunning %d QM Dimer jobs:\n",NDim+NDim_images);
  else if (ijob==MM_full-1)
    printf("\nRunning Full system MM job:\n");
  else if (ijob==MM_start-1)
    printf("Running %d MM Monomer and Dimer jobs:\n",NMon+NDim+NDim_images);
  
  
  fflush(stdout);
  
}

// Look up index of dimer containing monomers imon1 & imon2
int Cluster::DimerLookup(int imon1, int imon2) {
  int index = -1;
  bool match = false;
  int i = 1;
  
  while (match == false && i<=NDim) {
    int indA = Dimers[i].GetIndexA();
    int indB = Dimers[i].GetIndexB();
    if ( (imon1==indA && imon2==indB) || (imon1==indB && imon2==indA) ) {
      match = true;
      index = i;
      //printf("Dimer (%d,%d) has index %d\n",imon1,imon2,index);
    }
    i++;
  }
  
  
  if (match)  
    return index;
  else {
    printf("Cluster::DimerLookup() - Dimer (%d,%d) not found.\n",imon1,imon2);
    exit(1);
  }
}


// Compute the classical force-field energy of the entire system using
// ab initio-derived parameters
double Cluster::ComputeAIFFEnergy() {
  
  time_t start_AIFF, stop_AIFF;
  start_AIFF = time(NULL);
  
  bool do_es = Params::Parameters().DoAIFFElectrostatics();
  bool do_ind = Params::Parameters().DoAIFFInduction();
  bool do_2b_disp = Params::Parameters().DoAIFF2BodyDispersion();
  bool do_3b_disp = Params::Parameters().DoAIFF3BodyDispersion();
  
  
  double Eaiff, E_es_pol, E_2b_disp, E_3b_disp;
  
  printf("\n -----------------------------------------------------------------\n");
  printf("  ** Compute ab initio force-field energy for the full system **\n");
  printf(" -----------------------------------------------------------------\n");
  
  if (Params::Parameters().IsPeriodic()) {
    // Compute the electrostatic + induction/polarization energy
    if (do_es || do_ind) 
      E_es_pol = ComputePeriodicMultipoleInteractions();
    
    // Compute the 2-body intermolecular dispersion energy for the
    // periodic system.
    if (do_2b_disp) 
      E_2b_disp = ComputePeriodicTwoBodyDispersion();
    
    // Compute the 3-body intermolecular dispersion energy for the
    // periodic system.
    if (do_3b_disp) 
      E_3b_disp = ComputeManyBodyDispersion();
  }
  
  else {
    // Compute the electrostatic + induction/polarization energy 
    if (do_es || do_ind)
      E_es_pol = ComputeClusterMultipoleInteractions();  
    //E_es_pol = 0.0; printf("HACK: E_es_pol = 0.0\n");
    
    // Compute the 2-body intermolecular dispersion energy for the
    // full cluster.
    if (do_2b_disp) 
      E_2b_disp = ComputeTwoBodyDispersion();
    
    // Compute the 3-body intermolecular dispersion energy for the
    // full cluster.
    if (do_3b_disp)
      E_3b_disp = ComputeManyBodyDispersion(); 
  }
  
  if (!do_es || !do_ind || !do_2b_disp || !do_3b_disp) {
    string terms[2];
    terms[0] = "no"; terms[1] = "yes";
    printf("  AIFF terms included:\n");
    printf("       Electrostatics: %s\n",terms[do_es].c_str());
    printf("            Induction: %s\n",terms[do_ind].c_str());
    printf("    2-body dispersion: %s\n",terms[do_2b_disp].c_str());
    printf("    3-body dispersion: %s\n",terms[do_3b_disp].c_str());
  }
  
  // Sum up the individual contributions
  Eaiff = E_es_pol + E_2b_disp + E_3b_disp;
  printf("\n*** Total AIFF energy for full system = %12.4f kJ/mol\n\n",Eaiff*HartreesToKJpermole);
  
  
  
  stop_AIFF = time(NULL);
  double AIFF_time = difftime(stop_AIFF, start_AIFF);
  printf("  Time to evaluate AIFF on the full system = %0f sec\n",AIFF_time);
  
  return Eaiff;
}


// Compute the classical electrostatics/induction for the entire
// cluster
double Cluster::ComputeClusterMultipoleInteractions() {
  
  time_t start_AIFF, stop_AIFF;
  start_AIFF = time(NULL);
  
  printf("\nStep 1: Compute classical electrostatic interactions for the full cluster.\n");
  
  // Note, we assume all the geometric interaction matrices Tab and
  // DampedTab have already been computed and are stored in the
  // corresponding Dimer objects.
  double Ees = 0.0, Eind = 0.0;
  
  /* Step 1: Compute permanent multipole contributions.  */
  
  // This is a purely pairwise additive effect, so for aperiodic
  // systems, we just sum the 2-body permanent multipole interactions
  // we have already computed.  Since there is no many-body
  // contribution from the permanent multipoles, the only point in
  // computing this is for the case where we perform QM calculations
  // only among local pairs.  In those cases, we need classical
  // electrostatics at long ranges.
  for (int i=1;i<=NDim;i++) {
    Ees += Dimers[i].GetMMElectrostaticEnergy(); // in hartrees
  }
  
  // If not induction term, we're done.  Return to the calling
  // routine.
  bool do_ind = Params::Parameters().DoAIFFInduction();
  if (!do_ind) { 
    return Ees;
  }
  
  /* Step 2: Compute induced multipole contributions.  Note, only
     intermolecular effects are considered.  Intramolecular induction
     has presumably been accounted for in determining the distributed
     polarizabilities. */
  
  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction energy for the full cluster.\n");
  }
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 Params::Parameters().GetDampingFactor());
  
  // Create storage space for induced moments on all atoms
  int Natoms = GetTotalNumberOfAtoms(); 
  Multipole *dQ = new Multipole[Natoms]; 
  Multipole *old_dQ = new Multipole[Natoms];
  
  int count = 0;
  for (int imon=1;imon<=NMon;imon++) {
    for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
      // use larger of two ranks (see
      // Dimer::ComputeMultipoleInteraction() for explanation)
      int Rmom = Monomers[imon].GetAtom(iA).GetMultipoleMoments().GetRank();
      int Rpol = Monomers[imon].GetAtom(iA).GetPolarizability().GetRank();
      int max_rank = max(Rmom,Rpol);

      dQ[count].Initialize(max_rank);
      old_dQ[count].Initialize(max_rank);
      count++;
    }
  }
  
  // Create a list of the starting index for each monomer in the
  // induced moment list.  Index the same way as in Monomer list:
  // start from 1.
  int *multipole_key = new int[NMon+1]; 
  multipole_key[1] = 0; // set the first one by hand
  for (int i=2;i<=NMon;i++) {
    multipole_key[i] = multipole_key[i-1] + Monomers[i-1].GetNumberOfAtoms();
  }
  
  // Get ready to begin iterations - initialize a few variables
  bool iterate = true;
  int cycle = 0;
  double ind_conv = Params::Parameters().GetInductionConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  if (Params::Parameters().PrintLevel() >= 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);
  
  double Eind_old = 0.0;
  Eind = 1000000.0;  // start with huge nonsense energy

  if (Params::Parameters().PrintLevel() >= 0) {
    printf("--------------------------------------------------\n");
    printf(" Cycle       E(induction)       Change\n");
    printf("--------------------------------------------------\n");
  }
  
  // Begin the iterations to find self-consistent induction
  while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {
    
    cycle++;
    
    // Save data from previous cycle and reset the variables for this cycle
    Eind_old = Eind;
    Eind = 0.0;
    for (int iA=0;iA<Natoms;iA++) {
      old_dQ[iA] = dQ[iA];
      dQ[iA].Set();
    }
    
    // Induce multipoles
    
    // loop over monomers to be induced - "inducee"
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
      
      // find where the induced multipole moments for this monomer start
      int offsetA = multipole_key[imonA]; 
      
      // loop over all monomers that pair with this dimer - "inducers"
      for (int imonB = 1; imonB <= NMon; imonB++) {
	
	if (imonA != imonB) {
	  int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	  // find where the induced multipole moments for this monomer start
	  int offsetB = multipole_key[imonB]; 
	  
	  // Set logical flag for if (imonA < imonB).  Important because
	  // Dimer objects are always stored with A<B, and there is
	  // some directionality implicit in the Tab interaction
	  // matrices used to induce the multipoles.  In other words,
	  // Tab for BA is transpose of Tab for AB, and we need to
	  // be sure to grab the proper one.
	  bool AB_order = true;
	  if (imonA > imonB) AB_order = false; 
	  
	  // Grab the index of the appropriate dimer -- need this to
	  // obtain the Tab matrices
	  int idimer = DimerLookup(imonA,imonB); 
	  
	  // Loop over atoms on each monomer
	  for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	    
	    Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments());
	    Polarizability PolA(Monomers[imonA].GetAtom(iA).GetPolarizability(),true);
	    int NpolA = PolA.GetLength();
	    int NmomA = QA.GetLength();
	    int dimA = min(NpolA,NmomA);
	    
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      
	      
	      Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	      int NmomB = QB.GetLength();
	      int NpolB = Monomers[imonB].GetAtom(iB).GetPolarizability().GetLength();
	      
	      Matrix Tab;
	      if (AB_order) {	
		Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	      }
	      else {
		Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		Tab.Transpose();
	      }
	      int dimT1 = Tab.GetRows();
	      int dimT2 = Tab.GetCols();
	      
	      // Induce multipoles on monomer A due to monomer B
	      // dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))  (B != A)
	      for (int a=0; a < NpolA; a++) // loop over elements of dQA
		for (int t=0; t < min(NpolA,dimT1); t++) 
		  for (int u=0; u < min(dimT2,NmomB); u++) { 
		    dQ[offsetA+iA](a) 
		      -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQ[offsetB+iB](u));
		  }
	      
	    } // end loop over inducer atoms
	  } // end loop over inducee atoms
	} // end if (imonA != imonB)
      } // end loop over inducer monomers
    } // end loop over inducee monomers
    
    
    // Now compute the energy contribution from the induced multipoles
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
      int offsetA = multipole_key[imonA]; 
      // loop over all monomers that pair with this dimer - "inducers"
      for (int imonB = 1; imonB <= NMon; imonB++) {
	if (imonA != imonB) {
	  int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	  int offsetB = multipole_key[imonB]; 
	  bool AB_order = true;
	  if (imonA > imonB) AB_order = false; 
	  
	  // Grab the index of the appropriate dimer -- need this to
	  // obtain the Tab matrices
	  int idimer = DimerLookup(imonA,imonB); 
	  
	  // Loop over atoms on each monomer
	  for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      
	      
	      Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	      
	      Matrix Tab;
	      if (AB_order) {
		Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		int dimT1 = Tab.GetRows();
		int dimT2 = Tab.GetCols();
		
		for (int t=0;t<min(dQ[offsetA+iA].GetLength(),dimT1);t++) 
		  for (int u=0;u<min(QB.GetLength(),dimT2);u++) {
		    Eind += 0.5*dQ[offsetA + iA](t)*Tab(t,u)*QB(u);
		  }
	      }
	      
	      else {
		Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		int dimT1 = Tab.GetRows();
		int dimT2 = Tab.GetCols();
		for (int t=0;t<min(dQ[offsetA + iA].GetLength(),dimT2);t++) 
		  for (int u=0;u<min(QB.GetLength(),dimT1);u++) {
		    Eind += 0.5*dQ[offsetA + iA](t)*Tab(u,t)*QB(u);
		  }	
	      } // end if (AB_order) ... else
	      
	    } // end loop over inducer atoms
	  } // end loop over inducee atoms
	} // end if (imonA != imonB)
      } // end loop over inducer monomers
    } // end loop over inducee monomers
    
    
    // Print out results for this cycle
    if (cycle==1 && Params::Parameters().PrintLevel() >= 0) 
      printf(" %3d       %12.6f     *********** kJ/mol\n",cycle,Eind*HartreesToKJpermole);
    else if (Params::Parameters().PrintLevel() >= 0) 
      printf(" %3d       %12.6f     %11.6f kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole);
    
    // Check convergence based on the energy change.
    if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv) {
      if (Params::Parameters().PrintLevel() >= 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies converged after %d iterations\n\n",cycle);
      iterate = false;
    }    
    
    if ( cycle == Params::Parameters().GetMaxPolarizationCycles()  && iterate == true ) {
      if (Params::Parameters().PrintLevel() >= 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies failed to converge after %d iterations\n\n",cycle);
      Params::Parameters().Warning();
     }
    
    fflush(stdout);
  } // end while loop

  
  // Hack: looking at what happens if we don't use damping for final induction energy
  if (Params::Parameters().GetMaxOptCycles()==1) {
    printf("Recomputing cluster induction energy without damping\n");
    Eind = 0.0;
    // Now compute the energy contribution from the induced multipoles
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
      int offsetA = multipole_key[imonA]; 
      // loop over all monomers that pair with this dimer - "inducers"
      for (int imonB = 1; imonB <= NMon; imonB++) {
	if (imonA != imonB) {
	  int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	  int offsetB = multipole_key[imonB]; 
	  bool AB_order = true;
	  if (imonA > imonB) AB_order = false; 
	  
	  // Grab the index of the appropriate dimer -- need this to
	  // obtain the Tab matrices
	  int idimer = DimerLookup(imonA,imonB); 
	  
	  // Loop over atoms on each monomer
	  for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

	      Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	      
	      Matrix Tab;
	      if (AB_order) {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
		int dimT1 = Tab.GetRows();
		int dimT2 = Tab.GetCols();
		
		for (int t=0;t<min(dQ[offsetA+iA].GetLength(),dimT1);t++) 
		  for (int u=0;u<min(QB.GetLength(),dimT2);u++) {
		    Eind += 0.5*dQ[offsetA + iA](t)*Tab(t,u)*QB(u);
		  }
	      }
	      
	      else {
		Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		int dimT1 = Tab.GetRows();
		int dimT2 = Tab.GetCols();
		for (int t=0;t<min(dQ[offsetA + iA].GetLength(),dimT2);t++) 
		  for (int u=0;u<min(QB.GetLength(),dimT1);u++) {
		    Eind += 0.5*dQ[offsetA + iA](t)*Tab(u,t)*QB(u);
		  }	
	      } // end if (AB_order) ... else
	      
	    } // end loop over inducer atoms
	  } // end loop over inducee atoms
	} // end if (imonA != imonB)
      } // end loop over inducer monomers
    } // end loop over inducee monomers
    
  }
  
  
  

  
  
  
  
  
  
  
  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles ***\n");
    count = 0;
    for (int imon=1;imon<=NMon;imon++) {
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	char label[50];
	sprintf(label,"Monomer %d, %s%d",imon,
		Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		Monomers[imon].GetAtom(iA).GetAtomIndex());
	string str = label;
	dQ[count].Print(str);
	count++;
      }
    }
  }

  delete [] dQ;
  delete [] old_dQ;
  delete [] multipole_key;
  
  
  
  // Store the results - in hartrees.
  E_Electrostatic_MM = Ees;
  E_Induction_MM = Eind;
  double Etot = Ees + Eind;
  
  // Print out a summary of the classical energy just calculated
  if (Params::Parameters().PrintLevel() >= 0) {
    printf("--------------------------------------------------------\n");
    printf("  Total Classical ES/Ind Energy for the Full Cluster\n");
    printf("--------------------------------------------------------\n");
    printf("        Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
    printf("        Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
    printf("      ---------------------------------------\n");
    printf("        Sum           = %12.6f kJ/mol\n",Etot*HartreesToKJpermole);
    printf("--------------------------------------------------------\n");
  }
  
  
  stop_AIFF = time(NULL);
  double AIFF_time = difftime(stop_AIFF, start_AIFF);
  printf("     Time to evaluate electrostatics & induction in the full system = %0f sec\n",AIFF_time);
  
  return Etot;
}


 Vector Cluster::ComputeClusterMultipoleGradient() {

   printf("\nComputing nuclear gradient of classical electrostatic interactions for\n    the full cluster.\n");

   // Initialize storage space for the gradient
   int Natoms = GetTotalNumberOfAtoms();
   Vector Grad(3*Natoms);

   // Create a list of the starting index for each monomer in the full
   // gradient.  Like the Monomers, the key indexes count from 1->NMon.
   // On the other hand, in the gradient, the indexing starts at 0.
   int *key = new int[NMon+1];
   key[1] = 0; // set the first one by hand
   for (int i=2;i<=NMon;i++) {
     key[i] = key[i-1] + 3*Monomers[i-1].GetNumberOfAtoms();
   }

   // Note, we assume all the geometric interaction matrices Tab and
   // DampedTab have already been computed and are stored in the
   // corresponding Dimer objects.
   if ( Params::Parameters().IsPeriodic() ) {
     printf("ERROR: Cluster::ComputeClusterMultipoleInteractions() - Periodic boundary conditions\n are not yet implemented\n");
     exit(1);
   }

   /* Step 1: Compute the nuclear gradient of the permanent multipole
      contributions.  */

   // This is a purely pairwise additive effect, so for aperiodic
   // systems, we just sum the 2-body permanent multipole interactions
   // we have already computed.  Since there is no many-body
   // contribution from the permanent multipoles, the only point in
   // computing this is for the case where we perform QM calculations
   // only among local pairs.  In those cases, we need classical
   // electrostatics at long ranges.

   double Ees = 0.0;
   for (int i=1;i<=NDim;i++) {

     // Energy:
     Ees += Dimers[i].GetMMElectrostaticEnergy(); // in hartrees

     // Gradient:
     Vector dimerGrad = Dimers[i].GetMM2BodyElectrostaticGradient(); 
     int indA = Dimers[i].GetIndexA();
     int indB = Dimers[i].GetIndexB();
     int NatomsA =  GetNumberOfAtoms(indA);
     int NatomsB =  GetNumberOfAtoms(indB);
     int startA = Grad_key[indA];
     int startB = Grad_key[indB];

     for (int j=0;j<3*NatomsA;j++) { // contrib for first monomer
       Grad[ startA + j ] += dimerGrad[j];
     }
     for (int j=0;j<3*NatomsB;j++) { // contrib for second monomer
       Grad[ startB + j ] += dimerGrad[ 3*NatomsA + j ];
     }
     printf("Dimer(%i,%i)\n",indA,indB);
     dimerGrad.PrintGradient("");

   }

   //Grad.Scale(HartreesToKJpermole);
   Grad.PrintGradient("Electrostatic Gradient, full cluster");
   //Grad.Scale(1.0/HartreesToKJpermole);

   /* Step 2: Compute the nuclear gradient of the induction energy */
   if (Params::Parameters().PrintLevel() > 0) {
     printf("\n");
     printf("Computing self-consistent induction energy for the full cluster.\n");
   }
   //printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
   // Params::Parameters().GetDampingFactor());

   // Initialize a few variables
   bool iterate = true;
   int cycle = 0;

   double ind_conv = Params::Parameters().GetInductionConvergence();
   double ind_gradconv = Params::Parameters().GetInductionGradConvergence();	
   //double Econv = 1.0e-7; // induction convergence threshold
   double Econv = 1.0/pow(10,ind_conv);
   double Gconv = 1.0/pow(10,ind_gradconv);
   if (Params::Parameters().PrintLevel() > 0) 
     printf("  -----ind_conv = %12.6f, Econv = %12.6f, ind_gradconv = %12.6f, Gconv = %12.6f\n", ind_conv,Econv,ind_gradconv,Gconv );

   double Eind_old = 0.0;
   double Eind;

   Vector Eind_grad(3*Natoms);
   Vector old_Eind_grad(3*Natoms);

   // Create storage space for induced moments on all atoms...
   Multipole *dQ = new Multipole[Natoms]; 
   Multipole *old_dQ = new Multipole[Natoms];

   // and their nuclear gradients
   Multipole *dQ_grad = new Multipole[3*Natoms*Natoms]; 
   Multipole *old_dQ_grad = new Multipole[3*Natoms*Natoms];

   int count = 0;
   for (int imon=1;imon<=NMon;imon++) {
     for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
       // use larger of two ranks (see
       // Dimer::ComputeMultipoleInteraction() for explanation)
       int Rmom = Monomers[imon].GetAtom(iA).GetMultipoleMoments().GetRank();
       int Rpol = Monomers[imon].GetAtom(iA).GetPolarizability().GetRank();
       int max_rank = max(Rmom,Rpol);

       dQ[count].Initialize(max_rank);
       old_dQ[count].Initialize(max_rank);

       for (int X=0;X<3*Natoms;X++) {
	 dQ_grad[3*Natoms*count + X].Initialize(max_rank);
	 old_dQ_grad[3*Natoms*count + X].Initialize(max_rank);
       }
       count++;
     }
   }

   // Create a list of the starting index for each monomer in the
   // induced moment & gradient lists.  Index multipole_key the same
   // way as in Monomer list: start from 1.
   int *multipole_key = new int[NMon+1]; 
   int *multipole_grad_key = new int[NMon+1]; 

   multipole_key[1] = 0; // set the first one by hand
   multipole_grad_key[1] = 0; // set the first one by hand
   for (int i=2;i<=NMon;i++) {
     multipole_key[i] = multipole_key[i-1] + Monomers[i-1].GetNumberOfAtoms();
     multipole_grad_key[i] = multipole_grad_key[i-1] 
       + 3*Natoms*Monomers[i-1].GetNumberOfAtoms();
   }

   // Get ready to begin iterations - initialize a few variables
   Eind = 1000000.0;  // start with huge nonsense energy

   if (Params::Parameters().PrintLevel() > 0) {
     printf("------------------------------------------------------------------------\n");
     printf(" Cycle     E(induction)       Change          |Grad|          Change   \n");
     printf("------------------------------------------------------------------------\n");
   }

   // Begin the iterations to find self-consistent induction
   while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {

     cycle++;

     // Save data from previous cycle and reset the variables for this cycle
     Eind_old = Eind;
     Eind = 0.0;
     old_Eind_grad = Eind_grad;
     Eind_grad.Set();

     for (int iA=0;iA<Natoms;iA++) {
       old_dQ[iA] = dQ[iA];
       dQ[iA].Set();

       for (int X=0;X<3*Natoms;X++) {
	 old_dQ_grad[3*Natoms*iA+X] = dQ_grad[3*Natoms*iA+X];
	 dQ_grad[3*Natoms*iA+X].Set();
       }
     }

     // Induce multipoles

     // loop over monomers to be induced - "inducee"
     for (int imonA = 1; imonA <= NMon; imonA++) {
       int NatomsA = Monomers[imonA].GetNumberOfAtoms();

       // find where the induced multipole moments for this monomer start
       int offsetA = multipole_key[imonA]; 
       int goffsetA = multipole_grad_key[imonA];

       // loop over all monomers that pair with this dimer - "inducers"
       for (int imonB = 1; imonB <= NMon; imonB++) {

	 if (imonA != imonB) {
	   int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	   // find where the induced multipole moments for this monomer start
	   int offsetB = multipole_key[imonB]; 
	   int goffsetB = multipole_grad_key[imonB];

	   // Set logical flag for if (imonA < imonB).  Important
	   // because Dimer objects are always stored with A<B, and
	   // directionality is implicit in the Tab interaction
	   // matrices used to induce the multipoles.  In other words,
	   // Tab for BA is transpose of Tab for AB, and we need to be
	   // sure to grab the proper one.
	   bool AB_order = true;
	   if (imonA > imonB) AB_order = false; 

	   // Grab the index of the appropriate dimer -- need this to
	   // obtain the Tab matrices
	   int idimer = DimerLookup(imonA,imonB); 
	   // Loop over atoms on each monomer
	   for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 

	     Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments());
	     Polarizability PolA(Monomers[imonA].GetAtom(iA).GetPolarizability(),true);
	     int NpolA = PolA.GetLength();
	     int NmomA = QA.GetLength();
	     int dimA = min(NpolA,NmomA);

	     for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

	       Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	       int NmomB = QB.GetLength();
	       int NpolB = Monomers[imonB].GetAtom(iB).GetPolarizability().GetLength();

	       // Read Tab, dTdX.  If iB < iA, we have to read T(iB,iA), and transpose.
	       // Like wise for dTdX.
	       Matrix Tab, dTdX;
	       if (AB_order) {
		 // for Tab
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));

		 // for dTdX:
		 Matrix small_dTdX(Dimers[idimer].GetDampedTabGradient(iA,iB));
		 int rows = small_dTdX.GetRows();
		 dTdX.Initialize(rows,3*Natoms);

		 for (int irows=0;irows<rows;irows++) 
		   for (int X=0;X<3;X++) {
		     dTdX(irows,offsetA*3 + iA*3 + X) = small_dTdX(irows,X);
		     dTdX(irows,offsetB*3 + iB*3 + X) = small_dTdX(irows,X+3); 
		   }

	       }
	       else {
		 // for Tab
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		 Tab.Transpose();

		 // for dTdX
		 Matrix small_dTdX(Dimers[idimer].GetDampedTabGradient(iB,iA));
		 int rows = small_dTdX.GetRows();
		 dTdX.Initialize(rows,3*Natoms);

		 // Loop over t/u separately to transpose those two indices in gradient
		 // We also have to "transpose" iA and iB in the X dimension.
		 for (int t=0;t<NmomB; t++)
		   for (int u=0;u<NmomA; u++) {
		     int tu = u*NmomB + t; // compound (t,u) index
		     int ut = t*NmomA + u; // compound (u,t) index

		     for (int X=0;X<3;X++) {
		       dTdX(ut,offsetB*3 + iB*3 + X) = small_dTdX(tu,X);
		       dTdX(ut,offsetA*3 + iA*3 + X) = small_dTdX(tu,X+3); 
		     }
		   }
	       }
	       int dimT1 = Tab.GetRows();
	       int dimT2 = Tab.GetCols();

	       // Induce multipoles on monomer A due to monomer B
	       // dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))  (B != A)
	       for (int a=0; a < NpolA; a++) // loop over elements of dQA
		 for (int t=0; t < min(NpolA,dimT1); t++) 
		   for (int u=0; u < min(dimT2,NmomB); u++) { 
		     // Induced moments:
		     dQ[offsetA+iA](a) 
		       -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQ[offsetB+iB](u));

		     // Gradient of induced moments:
		     int tu = u*Tab.GetRows() + t; // define a compound index for (t,u)
		     for (int X=0;X<3*Natoms;X++) {
		       dQ_grad[goffsetA + iA*3*Natoms + X](a) -= PolA(a,t)* 
		       (dTdX(tu,X)*(QB(u) + old_dQ[offsetB+iB](u)) 
			+ Tab(t,u)*old_dQ_grad[goffsetB + iB*3*Natoms + X](u));
		     }
		   }

	     } // end loop over inducer atoms
	   } // end loop over inducee atoms
	 } // end if (imonA != imonB)
       } // end loop over inducer monomers
     } // end loop over inducee monomers

     // Now compute the energy contribution from the induced multipoles
     for (int imonA = 1; imonA <= NMon; imonA++) {
       int NatomsA = Monomers[imonA].GetNumberOfAtoms();
       int offsetA = multipole_key[imonA]; 
       int goffsetA = multipole_grad_key[imonA];
       // loop over all monomers that pair with this dimer - "inducers"
       for (int imonB = 1; imonB <= NMon; imonB++) {
	 if (imonA != imonB) {
	   int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	   int offsetB = multipole_key[imonB]; 
	   int goffsetB = multipole_grad_key[imonB];
	   bool AB_order = true;
	   if (imonA > imonB) AB_order = false; 

	   // Grab the index of the appropriate dimer -- need this to
	   // obtain the Tab matrices
	   int idimer = DimerLookup(imonA,imonB); 

	   // Loop over atoms on each monomer
	   for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	     Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments());
	     int NmomA = QA.GetLength();

	     for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      
	       Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	       int NmomB = QB.GetLength();

	       // Read Tab, dTdX.  If iB < iA, we have to read T(iB,iA), and transpose.
	       // Like wise for dTdX.
	       Matrix Tab, dTdX;
	       if (AB_order) {
		 // for Tab
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));

		 // for dTdX:
		 Matrix small_dTdX(Dimers[idimer].GetDampedTabGradient(iA,iB));
		 int rows = small_dTdX.GetRows();
		 dTdX.Initialize(rows,3*Natoms);

		 for (int irows=0;irows<rows;irows++) 
		   for (int X=0;X<3;X++) {
		     dTdX(irows,offsetA*3 + iA*3 + X) = small_dTdX(irows,X);
		     dTdX(irows,offsetB*3 + iB*3 + X) = small_dTdX(irows,X+3); 
		   }

	       }
	       else {
		 // for Tab
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		 Tab.Transpose();

		 // for dTdX
		 Matrix small_dTdX(Dimers[idimer].GetDampedTabGradient(iB,iA));
		 int rows = small_dTdX.GetRows();
		 dTdX.Initialize(rows,3*Natoms);

		 // Loop over t/u separately to transpose those two indices in gradient
		 // We also have to "transpose" iA and iB in the X dimension.
		 for (int t=0;t<NmomB; t++)
		   for (int u=0;u<NmomA; u++) {
		     int tu = u*NmomB + t; // compound (t,u) index
		     int ut = t*NmomA + u; // compound (u,t) index

		     for (int X=0;X<3;X++) {
		       dTdX(ut,offsetB*3 + iB*3 + X) = small_dTdX(tu,X);
		       dTdX(ut,offsetA*3 + iA*3 + X) = small_dTdX(tu,X+3); 
		     }
		   }
	       }
	       int dimT1 = Tab.GetRows();
	       int dimT2 = Tab.GetCols();

	       for (int t=0; t<min(dQ[offsetA+iA].GetLength(),dimT1); t++) 
		 for (int u=0; u<min(QB.GetLength(),dimT2); u++) {
		   // Energy:
		   Eind += 0.5*dQ[offsetA + iA](t)*Tab(t,u)*QB(u);

		   // Gradient:
		   int tu = u*Tab.GetRows() + t; // define a compound index for (t,u)
		   for (int X=0;X<3*Natoms;X++) {
		     Eind_grad[X] +=
		       0.5*dQ_grad[goffsetA + iA*3*Natoms + X](t)*Tab(t,u)*QB(u)
		       + 0.5*dQ[offsetA + iA](t)*dTdX(tu,X)*QB(u);
		   }
		 }

	     } // end loop over inducer atoms
	   } // end loop over inducee atoms
	 } // end if (imonA != imonB)
       } // end loop over inducer monomers
     } // end loop over inducee monomers

     double Gradnorm = Eind_grad.Norm();
     Vector tmp(Eind_grad);
     tmp -= old_Eind_grad;
     double dGrad = tmp.Norm();

     // Print out results for this cycle
     if (cycle==1 && Params::Parameters().PrintLevel() > 0) 
       printf(" %3d     %12.6f       *********    %12.6f      ********* kJ/mol\n",cycle,
	      Eind*HartreesToKJpermole, Gradnorm*HartreesToKJpermole);
     else if (Params::Parameters().PrintLevel() > 0) 
       printf(" %3d     %12.6f     %11.6f    %12.6f    %11.6f kJ/mol\n",cycle,
	      Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole, 
	      Gradnorm*HartreesToKJpermole, dGrad*HartreesToKJpermole);

     // Check convergence based on the energy change (in kJ/mol).
     if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv 
	 && dGrad*HartreesToKJpermole < Gconv) {
       if (Params::Parameters().PrintLevel() > 0) 
	 printf("------------------------------------------------------------------------\n");
       printf("  Induction energies & gradients converged after %d iterations\n",cycle);
       if (Params::Parameters().PrintLevel() > 0) 
	 printf("\n");
       iterate = false;
     }    

     if ( cycle == Params::Parameters().GetMaxPolarizationCycles() && iterate == true ) {
       if (Params::Parameters().PrintLevel() > 0) 
	 printf("------------------------------------------------------------------------\n");
       printf("  Induction energies & gradients failed to converge after %d iterations\n",
	      cycle);
       Params::Parameters().Warning();
       if (Params::Parameters().PrintLevel() > 0) 
	 printf("\n");
     }
   } // end while loop

   //Eind_grad.Scale(HartreesToKJpermole);
   //Eind_grad.PrintGradient("Eind gradient");
   //Eind_grad.Scale(1.0/HartreesToKJpermole);

   // Add induction contribution to the full MM gradient
   Grad += Eind_grad;


   if (Params::Parameters().PrintLevel() > 1) {
     // Print out the final multipole moments
     printf(" *** Final induced multipoles ***\n");
     count = 0;
     for (int imon=1;imon<=NMon;imon++) {
       for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	 char label[50];
	 sprintf(label,"Monomer %d, %s%d",imon,
		 Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		 Monomers[imon].GetAtom(iA).GetAtomIndex());
	 string str = label;
	 dQ[count].Print(str);
	 count++;
       }
     }
   }


   // Store the results - in hartrees.
   E_Electrostatic_MM = Ees;
   E_Induction_MM = Eind;
   Energy_MM = Ees + Eind;

   // Print out a summary of the classical energy just calculated
   if (Params::Parameters().PrintLevel() > 0) {
     printf("------------------------------------------------------------\n");
     printf(" Total Classical Interaction Energy for the Full Cluster\n");
     printf("------------------------------------------------------------\n");
     printf("        Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
     printf("        Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
     printf("      ---------------------------------------\n");
     printf("        Total         = %12.6f kJ/mol\n",(Ees+Eind)*HartreesToKJpermole);
     printf("--------------------------------------------------------\n");
   }

   delete [] key;
   delete [] dQ;
   delete [] old_dQ;
   delete [] dQ_grad;
   delete [] old_dQ_grad;
   delete [] multipole_key;
   delete [] multipole_grad_key;

   Grad.PrintGradient("MM Gradiant");
     //Grad.Scale(HartreesToKJpermole);
     //PrintGradient("Full cluster MM gradient",Grad);
     //Grad.Scale(1.0/HartreesToKJpermole);

   //reducing elements of gradient using symmetry.

   //this key will indicate which atoms are in the gradient
   int gradkey[UniqueAtoms];
   int gradindex=0;
   for(int i=1; i<NMon;i++){
     if(Monomers[i].GetSymmetryFactor()!=0){
       for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
	 gradkey[gradindex] = Monomers[i].GetAtom(j).GetGlobalIndex();
	 gradindex++;
       }
     }
   }
   //Reducing the size of the gradient restricted  under symmetry.
   Vector ReducedGrad(3*UniqueAtoms);
   for(int iMon=1;iMon<=NMon;iMon++){
     Matrix RotMat = Monomers[iMon].GetRotationMatrix();
     // RotMat.Print("Rotation Matrix");
     for(int iatom=0; iatom<Monomers[iMon].GetNumberOfAtoms();iatom++){
       int SymGlobalIndex = Monomers[iMon].GetAtom(iatom).GetSymmetricalAtom();//Global index of the atom that the iAtom is symmetrical too
       for(int k=0; k<UniqueAtoms;k++){
	 if(SymGlobalIndex==gradkey[k]){//matching Global index of the atom it is symmetric to the index in the key
	   int GlobIndex = Monomers[iMon].GetAtom(iatom).GetGlobalIndex() - 1;//Global index used to get entree from Grad

	   printf("GlobIndex = %i SymGlobalIndex = %i\n",GlobIndex+1,SymGlobalIndex);

	   //rotating elements of the Gradient and placing it into the reduced gradient 
	   ReducedGrad[3*k]   += RotMat(0,0)*Grad[3*GlobIndex] + RotMat(0,1)*Grad[3*GlobIndex+1] + RotMat(0,2)*Grad[3*GlobIndex+2];
	   ReducedGrad[3*k+1] += RotMat(1,0)*Grad[3*GlobIndex] + RotMat(1,1)*Grad[3*GlobIndex+1] + RotMat(1,2)*Grad[3*GlobIndex+2];
	   ReducedGrad[3*k+2] += RotMat(2,0)*Grad[3*GlobIndex] + RotMat(2,1)*Grad[3*GlobIndex+1] + RotMat(2,2)*Grad[3*GlobIndex+2];
	 }
       }
     }
   }

   ReducedGrad.PrintGradient("Reduced Gradient");
   //exit(0);
   return ReducedGrad;
 }



 void Cluster::PreconvergeInducedMultipoleMomentsInUnitCell(Multipole dQ[], Multipole old_dQ[], int multipole_key[]) {

   time_t start_AIFF, stop_AIFF;
   start_AIFF = time(NULL);

   int Natoms = GetTotalNumberOfAtoms(); 

   if (Params::Parameters().PrintLevel() > 0) {
     printf("\n");
     printf("Computing self-consistent induction energy within just the unit cell.\n");
   }
   printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	  Params::Parameters().GetDampingFactor());

   // Get ready to begin iterations - initialize a few variables
   bool iterate = true;
   int cycle = 0;
   double ind_conv = Params::Parameters().GetInductionConvergence();
   double Econv = 1.0/pow(10,ind_conv);
   if (Params::Parameters().PrintLevel() >= 0) 
     printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);

   double Eind_old = 0.0;
   double Eind = 1000000.0;  // start with huge nonsense energy

   if (Params::Parameters().PrintLevel() >= 0) {
     printf("--------------------------------------------------\n");
     printf(" Cycle       E(induction)       Change\n");
     printf("--------------------------------------------------\n");
   }

   // Begin the iterations to find self-consistent induction
   while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {

     cycle++;

     // Save data from previous cycle and reset the variables for this cycle
     Eind_old = Eind;
     Eind = 0.0;
     for (int iA=0;iA<Natoms;iA++) {
       old_dQ[iA] = dQ[iA];
       dQ[iA].Set();
     }

     // Induce multipoles

     // loop over monomers to be induced - "inducee"
     for (int imonA = 1; imonA <= NMon; imonA++) {
       int NatomsA = Monomers[imonA].GetNumberOfAtoms();

       // find where the induced multipole moments for this monomer start
       int offsetA = multipole_key[imonA]; 

       // loop over all monomers that pair with this dimer - "inducers"
       for (int imonB = 1; imonB <= NMon; imonB++) {

	 if (imonA != imonB) {
	   int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	   // find where the induced multipole moments for this monomer start
	   int offsetB = multipole_key[imonB]; 

	   // Set logical flag for if (imonA < imonB).  Important because
	   // Dimer objects are always stored with A<B, and there is
	   // some directionality implicit in the Tab interaction
	   // matrices used to induce the multipoles.  In other words,
	   // Tab for BA is transpose of Tab for AB, and we need to
	   // be sure to grab the proper one.
	   bool AB_order = true;
	   if (imonA > imonB) AB_order = false; 

	   // Grab the index of the appropriate dimer -- need this to
	   // obtain the Tab matrices
	   int idimer = DimerLookup(imonA,imonB); 

	   // Loop over atoms on each monomer
	   for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 

	     Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments());
	     Polarizability PolA(Monomers[imonA].GetAtom(iA).GetPolarizability(),true);
	     int NpolA = PolA.GetLength();
	     int NmomA = QA.GetLength();
	     int dimA = min(NpolA,NmomA);

	     for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

	       Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	       int NmomB = QB.GetLength();
	       int NpolB = Monomers[imonB].GetAtom(iB).GetPolarizability().GetLength();

	       Matrix Tab;
	       if (AB_order) {
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	       }
	       else {
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		 Tab.Transpose();
	       }
	       int dimT1 = Tab.GetRows();
	       int dimT2 = Tab.GetCols();

	       // Induce multipoles on monomer A due to monomer B
	       // dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))  (B != A)
	       for (int a=0; a < NpolA; a++) // loop over elements of dQA
		 for (int t=0; t < min(NpolA,dimT1); t++) 
		   for (int u=0; u < min(dimT2,NmomB); u++) { 
		     dQ[offsetA+iA](a) 
		       -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQ[offsetB+iB](u));
		   }

	     } // end loop over inducer atoms
	   } // end loop over inducee atoms
	 } // end if (imonA != imonB)
       } // end loop over inducer monomers
     } // end loop over inducee monomers


     // Now compute the energy contribution from the induced multipoles
     for (int imonA = 1; imonA <= NMon; imonA++) {
       int NatomsA = Monomers[imonA].GetNumberOfAtoms();
       int offsetA = multipole_key[imonA]; 
       // loop over all monomers that pair with this dimer - "inducers"
       for (int imonB = 1; imonB <= NMon; imonB++) {
	 if (imonA != imonB) {
	   int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	   int offsetB = multipole_key[imonB]; 
	   bool AB_order = true;
	   if (imonA > imonB) AB_order = false; 

	   // Grab the index of the appropriate dimer -- need this to
	   // obtain the Tab matrices
	   int idimer = DimerLookup(imonA,imonB); 

	   // Loop over atoms on each monomer
	   for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	     for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

	       Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());

	       Matrix Tab;
	       if (AB_order) {
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		 int dimT1 = Tab.GetRows();
		 int dimT2 = Tab.GetCols();

		 for (int t=0;t<min(dQ[offsetA+iA].GetLength(),dimT1);t++) 
		   for (int u=0;u<min(QB.GetLength(),dimT2);u++) {
		     Eind += 0.5*dQ[offsetA + iA](t)*Tab(t,u)*QB(u);
		   }
	       }

	       else {
		 Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		 int dimT1 = Tab.GetRows();
		 int dimT2 = Tab.GetCols();
		 for (int t=0;t<min(dQ[offsetA + iA].GetLength(),dimT2);t++) 
		   for (int u=0;u<min(QB.GetLength(),dimT1);u++) {
		     Eind += 0.5*dQ[offsetA + iA](t)*Tab(u,t)*QB(u);
		   }	
	       } // end if (AB_order) ... else

	     } // end loop over inducer atoms
	   } // end loop over inducee atoms
	 } // end if (imonA != imonB)
       } // end loop over inducer monomers
     } // end loop over inducee monomers


     // Print out results for this cycle
     if (cycle==1 && Params::Parameters().PrintLevel() >= 0) 
       printf(" %3d       %12.6f     *********** kJ/mol\n",cycle,Eind*HartreesToKJpermole);
     else if (Params::Parameters().PrintLevel() >= 0) 
       printf(" %3d       %12.6f     %11.6f kJ/mol\n",cycle,
	      Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole);

     // Check convergence based on the energy change.
     if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv) {
       if (Params::Parameters().PrintLevel() >= 0) 
	 printf("--------------------------------------------------\n");
       printf("  Unit cell induction energies converged after %d iterations\n\n",cycle);
       iterate = false;
     }    

     if ( cycle == Params::Parameters().GetMaxPolarizationCycles()  && iterate == true ) {
       if (Params::Parameters().PrintLevel() >= 0) 
	 printf("--------------------------------------------------\n");
       printf("  Unit cell induction energies failed to converge after %d iterations\n\n",cycle);
       Params::Parameters().Warning();
     }

     fflush(stdout);
   } // end while loop

   if (Params::Parameters().PrintLevel() > 1) {
     // Print out the final multipole moments
     printf(" *** Converged induced multipoles in the non-periodic unit cell ***\n");
     int count = 0;
     for (int imon=1;imon<=NMon;imon++) {
       for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	 char label[50];
	 sprintf(label,"Monomer %d, %s%d",imon,
		 Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		 Monomers[imon].GetAtom(iA).GetAtomIndex());
	 string str = label;
	 dQ[count].Print(str);
	 count++;
       }
     }
   }

   // Copy the multipole moments over to the old ones, and zero the dQ
   // array
   for (int iA=0;iA<Natoms;iA++) {
     old_dQ[iA] = dQ[iA];
     dQ[iA].Set();
   }

   stop_AIFF = time(NULL);
   double AIFF_time = difftime(stop_AIFF, start_AIFF);
   printf("     Time to evaluate unit cell induction = %0f sec\n",AIFF_time);

 }


 // Compute the classical electrostatics/induction for the entire
 // cluster with periodic boundary conditions
 double Cluster::ComputePeriodicMultipoleInteractions() {

   printf("\nStep 1: Computing classical electrostatic interactions for the periodic system.\n");

   /* Step 1: Compute converged induced multipole moments for central
      unit cell atoms.  Note, only intermolecular effects are
      considered.  Intramolecular induction has presumably been
      accounted for in determining the distributed polarizabilities.

      For the periodic boundary condition case, we fully
      self-consistently induce within the central unit cell, just as in
      the non-periodic case.  But we also we create a sphere of image
      monomers around the central unit cell.  These image monomers
      induce moments on the central unit cell monomers.  After each
      cycle, the induced moments on the image monomers are constrained
      to be identical to those on the central unit cell monomers.  So
      we actually never explicitly induce moments on the image
      monomers.  We just grab the induced moments from the central unit
      cell molecules.
   */  
   double Ees = 0.0, Eind = 0.0;  

   bool do_ind = Params::Parameters().DoAIFFInduction();

   double deltaa = 0.0;

   if (do_ind) { 

     time_t start_ind_mom_time, stop_ind_mom_time;
     start_ind_mom_time = time(NULL);
     //printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
     // Params::Parameters().GetDampingFactor());


     // Create storage space for induced moments on all atoms in central
     // unit cell The induced moments on image monomers are identical to
     // these, so no need to store those separately
     int Natoms = GetTotalNumberOfAtoms(); 
     Multipole *dQ = new Multipole[Natoms];  // final storage of the latest multipole moments
     Multipole *old_dQ = new Multipole[Natoms]; // stores multipole moments from previous iteration
     Multipole *tmp_dQ = new Multipole[Natoms]; // temporarily stores moments from the current iteration

     int count = 0;
     for (int imon=1;imon<=NMon;imon++) {
       for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	 // use larger of two ranks (see
	 // Dimer::ComputeMultipoleInteraction() for explanation)
	 int Rmom = Monomers[imon].GetAtom(iA).GetMultipoleMoments().GetRank();
	 int Rpol = Monomers[imon].GetAtom(iA).GetPolarizability().GetRank();
	 int max_rank = max(Rmom,Rpol);

	 dQ[count].Initialize(max_rank);
	 old_dQ[count].Initialize(max_rank);
	 tmp_dQ[count].Initialize(max_rank);
	 count++;
       }
     }

     // Create a list of the starting index for each central unit cell
     // monomer in the induced moment list.  Index the same way as in
     // Monomer list: start from 1.
     int *multipole_key = new int[NMon+1]; 
     multipole_key[1] = 0; // set the first one by hand
     for (int i=2;i<=NMon;i++) {
       multipole_key[i] = multipole_key[i-1] + Monomers[i-1].GetNumberOfAtoms();
     }

     // Optionally find the converged moments in the (non-periodic)
     // unit cell first before adding the effect of image monomers.  In
     // practice, it doesn't seem to do much to speed convergence.
     if ( Params::Parameters().PreconvergeInducedMomentsInUnitCell() ) {
       PreconvergeInducedMultipoleMomentsInUnitCell(dQ, old_dQ, multipole_key);
     }     

     double r_cutoff = Params::Parameters().GetMaxPolarizationRadius();
     printf("\n");
     printf(" - Step 1a: Compute self-consistent induced multipole moments in large finite cluster.\n");
     printf("        Cutoff = %.1f Angstroms\n",r_cutoff);
     //  }
     // Start wall clock timer
     time_t start_time, stop_time;
     start_time = time(NULL);
     CreatePeriodicImageMonomerList(r_cutoff);
     // Stop the timer and print out the time
     stop_time = time(NULL);
     double elapsed_time = difftime(stop_time,start_time);
     printf("  Time to generate image monomer list = %.2f seconds\n",elapsed_time);

     if (Params::Parameters().PrintLevel() > 0)
       printf("Considering %d image monomers (cutoff = %.1f) for inducing multipoles\n",
	      NMon_images,Params::Parameters().GetMaxPolarizationRadius());

     // Get ready to begin iterations - initialize a few variables
     bool iterate = true;
     int cycle = 0;

     double ind_conv = Params::Parameters().GetInductionConvergence();
     double Econv = 1.0/pow(10,ind_conv);
     if (Params::Parameters().PrintLevel() > 0) 
       printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);

     double Eind_old = 0.0;
     Eind = 1000000.0;  // start with huge nonsense energy

     int Nskip = 0;

     printf("--------------------------------------------------\n");
     printf(" Cycle       E(induction)       Change\n");
     printf("--------------------------------------------------\n");
     double iter_scaling = Params::Parameters().GetInductionIterScaling();
     if (iter_scaling != 1.0) 
       printf("  Note: Scaling induced multipole moment iterations by %.2f\n",iter_scaling);

     // Begin the iterations to find self-consistent induction
     while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {

       cycle++;

       // Save data from previous cycle and reset the variables for this cycle
       Eind_old = Eind;
       Eind = 0.0;
       for (int iA=0;iA<Natoms;iA++) {
	 if (!(Params::Parameters().PreconvergeInducedMomentsInUnitCell() && cycle==1)){
	   old_dQ[iA] = dQ[iA];
	 }
	 else {
	   if (iA==0) printf("Using earlier guess for induced multipole moments\n");
	 }
	 dQ[iA].Set();
	 tmp_dQ[iA].Set();
       }

       // Induce multipoles

       // loop over central unit cell monomers to be induced - "inducee"
       for (int imonA = 1; imonA <= NMon; imonA++) {
	 int NatomsA = Monomers[imonA].GetNumberOfAtoms();

	 // find where the induced multipole moments for this monomer start
	 int offsetA = multipole_key[imonA]; 

	 // loop over all monomers that pair with this dimer - "inducers"
	 //    includes both central unit cell and image monomers
	 for (int imonB = 1; imonB <= NMon+NMon_images; imonB++) {
	   //printf("imonA = %d, imonB = %d\n",imonA,imonB);
	   // Grab the reference monomer, which we will use a lot.
	   bool IsImageMonomer;
	   int ref_monB;
	   if (imonB > NMon) {
	     IsImageMonomer = true;
	     ref_monB = MonomerImages[imonB-NMon].GetReferenceMonomerIndex();
	   }
	   else {
	     IsImageMonomer = false;
	     ref_monB = Monomers[imonB].GetReferenceMonomerIndex();
	   }      

	   if (imonA != imonB) {
	     int NatomsB = Monomers[ref_monB].GetNumberOfAtoms();
	     // find where the induced multipole moments for this monomer start
	     int offsetB = multipole_key[ref_monB]; 

	     // Set logical flag for if (imonA < imonB).  Important
	     // because Dimer objects are always stored with A<B, and
	     // there is some directionality implicit in the Tab
	     // interaction matrices used to induce the multipoles.  In
	     // other words, Tab for BA is transpose of Tab for AB, and
	     // we need to be sure to grab the proper one.  This only
	     // matters for dimers that are fully contained in the
	     // central unit cell.  If one of the monomers is an image
	     // monomer, it will be monomer B and must have a higher
	     // index than monomer A.
	     bool AB_order = true;
	     if (imonA > imonB) AB_order = false; 

	     // Grab the index of the appropriate dimer -- need this to
	     // obtain the Tab matrices.  If it involves an image
	     // monomer, we create the dimer and evaluate the Tab
	     // matrices right now.
	     //
	     // Future work: have to decide on whether it is better to
	     // store the Tab matrices in memory, on disk, or compute
	     // repeatedly on the fly.  Right now, just recompute them
	     // as needed.
	     int idimer; 	  
	     Dimer tmp_dimer;
	     bool TreatThisDimer = false;
	     if ( !IsImageMonomer) {
	       idimer = DimerLookup(imonA,imonB); 
	       TreatThisDimer = true;

              //if seperation of dimers is greater than c0, its Damped interaction Matrices do not exist
	      double seperation = Dimers[idimer].GetDimerSeparation();
              double c0 = Params::Parameters().GetLocalCutoff(0);
	      if(!Params::Parameters().DoLocal2BodyTruncation() || seperation >= c0)
	       Dimers[idimer].BuildDampedTabInteractionMatrices();
	     }
	     else {
	       // Create the dimer
	       tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB-NMon]);
	       //if (tmp_dimer.GetDimerSeparation() <= Params::Parameters().GetMaxPolarizationRadius() ) {
	       tmp_dimer.BuildDampedTabInteractionMatrices();
	       TreatThisDimer = true;
	       //}	    

	     }
	     if (!TreatThisDimer) Nskip++;

	     if (TreatThisDimer) {

	       // Loop over atoms on each monomer
	       for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 

		 Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments());
		 Polarizability PolA(Monomers[imonA].GetAtom(iA).GetPolarizability(),true);
		 int NpolA = PolA.GetLength();
		 int NmomA = QA.GetLength();
		 int dimA = min(NpolA,NmomA);

		 for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      
		   // Note these are always grabbed from the reference central unit cell monomer.
		   // This constrains image monomers to have same moments
		   Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());
		   int NmomB = QB.GetLength();
		   int NpolB = Monomers[ref_monB].GetAtom(iB).GetPolarizability().GetLength();

		   Matrix Tab;
		   if (AB_order) {
		     // If it's a real monomer, we just grab precomputed Tab
		     if (!IsImageMonomer) {
	             Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		     // Otherwise, we grab the Tab we just computed above
		   
		     }else {
		       Tab.Initialize(tmp_dimer.GetDampedTabInteractionMatrix(iA,iB));
		     }
		   }

		   else {  // this case only occurs for both monomers in
		     // central unit cell
		     Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		     Tab.Transpose();
		   }
		   int dimT1 = Tab.GetRows();
		   int dimT2 = Tab.GetCols();


		   // Induce multipoles on monomer A due to monomer B
		   // dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))  (B != A)
		   for (int a=0; a < NpolA; a++) {// loop over elements of dQA
		     for (int t=0; t < min(NpolA,dimT1); t++) {
		       for (int u=0; u < min(dimT2,NmomB); u++) { 

			 tmp_dQ[offsetA+iA](a) 
			   -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQ[offsetB+iB](u));

			 //original code, no iteration scaling
			 /*
			   dQ[offsetA+iA](a) 
			   -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQ[offsetB+iB](u));
			 */
		       } // end loop over u
		     } // end loop over t
		   } // end loop over a
		 } // end loop over inducer atoms iB
	       } // end loop over inducee atoms iA
	     } // end if (TreatThisDimer)
	   } // end if (imonA != imonB)
	 } // end loop over inducer monomers
       } // end loop over inducee monomers


       // Now update the induced multipoles, with scaling if requested.
       // The scaling mixes the new induced multipoles with those from
       // the previous iteration.  This can help cases with poor
       // convergence.
       count = 0;
       iter_scaling = Params::Parameters().GetInductionIterScaling();
       if (cycle==1) iter_scaling = 1.0; // No scaling on the first cycle
       for (int imon=1;imon<=NMon;imon++) {
	 for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	   // If no scaling, just copy over the multipole moments from tmp_dQ to dQ
	   if (iter_scaling ==1.0) {
	     dQ[count] = tmp_dQ[count];
	   }
	   else { // if scaling

	     // Grab the multipole moments
	     Multipole New = tmp_dQ[count];
	     Multipole Old = old_dQ[count];

	     // Mix the old and the new, as appropriate
	     New.Scale(iter_scaling);
	     Old.Scale(1.0-iter_scaling);
	     dQ[count] = New + Old;
	   }
	   count++;
	 }
       }



       // Now compute the energy contribution from the induced multipoles
       // Loop over all monomers in the central unit cell
       for (int imonA = 1; imonA <= NMon; imonA++) {
	 int NatomsA = Monomers[imonA].GetNumberOfAtoms();
	 int offsetA = multipole_key[imonA]; 
	 // loop over all monomers that pair with this dimer - "inducers"
	 //    These are central unit cell and image monomers
	 for (int imonB = 1; imonB <= NMon+NMon_images; imonB++) {

	   // Classify the monomer as image or central unit cell
	   bool IsImageMonomer;
	   int ref_monB;
	   if (imonB > NMon) {
	     IsImageMonomer = true;
	     ref_monB = MonomerImages[imonB-NMon].GetReferenceMonomerIndex();
	   }
	   else {
	     IsImageMonomer = false;
	     ref_monB = Monomers[imonB].GetReferenceMonomerIndex();
	   }      

	   if (imonA != imonB) {
	     int NatomsB = Monomers[ref_monB].GetNumberOfAtoms();
	     int offsetB = multipole_key[ref_monB]; 
	     bool AB_order = true;
	     if (imonA > imonB) AB_order = false; 

	     // Grab the index of the appropriate dimer -- need this to
	     // obtain the Tab matrices
	     int idimer; 	  
	     Dimer tmp_dimer;
	     if ( !IsImageMonomer) 
	       idimer = DimerLookup(imonA,imonB);
	     else {
	       // Create the dimer
	       tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB-NMon]);
	       tmp_dimer.BuildDampedTabInteractionMatrices();
	     }

	     // Loop over atoms on each monomer
	     for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	       for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

		 Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());

		 Matrix Tab;
		 if (AB_order) {
		   // If it's a real monomer, we just grab precomputed Tab
		   if (!IsImageMonomer) 
		     Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		   // Otherwise, we grab the Tab we just computed above
		   else {
		     Tab.Initialize(tmp_dimer.GetDampedTabInteractionMatrix(iA,iB));
		   }
		   int dimT1 = Tab.GetRows();
		   int dimT2 = Tab.GetCols();

		   for (int t=0;t<min(dQ[offsetA+iA].GetLength(),dimT1);t++) 
		     for (int u=0;u<min(QB.GetLength(),dimT2);u++) {
		       Eind += 0.5*dQ[offsetA + iA](t)*Tab(t,u)*QB(u);
		     }
		 }

		 else {
		   Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		   int dimT1 = Tab.GetRows();
		   int dimT2 = Tab.GetCols();
		   for (int t=0;t<min(dQ[offsetA + iA].GetLength(),dimT2);t++) 
		     for (int u=0;u<min(QB.GetLength(),dimT1);u++) {
		       Eind += 0.5*dQ[offsetA + iA](t)*Tab(u,t)*QB(u);
		     }	
		 } // end if (AB_order) ... else

	       } // end loop over inducer atoms
	     } // end loop over inducee atoms
	   } // end if (imonA != imonB)
	 } // end loop over inducer monomers imonB
       } // end loop over inducee monomers imonA


       // Print out results for this cycle
       if (cycle==1 && Params::Parameters().PrintLevel() >= 0) 
	 printf(" %3d       %12.6f     *********** kJ/mol\n",cycle,Eind*HartreesToKJpermole);
       else if (Params::Parameters().PrintLevel() >= 0) 
	 printf(" %3d       %12.6f     %11.6f kJ/mol\n",cycle,
		Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole);

       // Check convergence based on the energy change.
       if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv) {
	 if (Params::Parameters().PrintLevel() >= 0) 
	   printf("--------------------------------------------------\n");
	 printf("  Induced multipole moments converged after %d iterations\n\n",cycle);
	 iterate = false;
       }    

       if ( cycle == Params::Parameters().GetMaxPolarizationCycles()  && iterate == true ) {
	 if (Params::Parameters().PrintLevel() > 0) 
	   printf("--------------------------------------------------\n");
	 printf("  Induced multipole moments failed to converge after %d iterations\n\n",cycle);
	 Params::Parameters().Warning();
       }

       fflush(stdout);
     } // end while loop



     //if (Params::Parameters().PrintLevel() > 1) {
     // Print out the final multipole moments
     printf(" *** Final induced multipoles ***\n");
     count = 0;
     for (int imon=1;imon<=NMon;imon++) {
       for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	 char label[50];
	 sprintf(label,"Monomer %d, %s%d",imon,
		 Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		 Monomers[imon].GetAtom(iA).GetAtomIndex());
	 string str = label;
	 dQ[count].Print(str);
	 //by shuhao
	 Monomers[imon].GetAtom(iA).SetInduceMultipoleMoments(dQ[count]);
	 //printf("-----SetInduceMultipoleMoments for this atom-----\n");
	 //Monomers[imon].GetAtom(iA).GetInduceMultipoleMoments().Print(str);   
	 count++;
       }
     }
     //}

     delete [] dQ;
     delete [] old_dQ;
     delete [] multipole_key;

     // Print out the induction energy just calculated based on the cluster originated by the polarization cutoff
     if (Params::Parameters().PrintLevel() > 0) {
       printf("------------------------------------------------------------\n");
       printf(" Induction Energy for the Cluster in the polarization cutoff\n");
       printf("------------------------------------------------------------\n");
       printf("        Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
       printf("------------------------------------------------------------\n");
     }
     stop_ind_mom_time = time(NULL);
     double ind_mom_time = difftime(stop_ind_mom_time, start_ind_mom_time);
     printf("Time to obtain induced moments = %.2f seconds\n",ind_mom_time);


     /* Step 1b: Compute the undamped energy in the short-range distance (~10 bohr or 6 Angstrom) 

     For this case, we also calculate the fully
     self-consistently induce energy in the polarization cut-off just as in
     the non-periodic case using the undamped tab matrix.  We create a sphere of image
     monomers around the central unit cell.  These image monomers
     induce moments on the central unit cell monomers.  After each
     cycle, the induced moments on the image monomers are constrained
     to be identical to those on the central unit cell monomers.  So we
     actually never explicitly induce moments on the image monomers.  We just
     grab the induced moments from the central unit cell molecules.
     */
     double Ees1 = 0.0, Eind1 = 0.0;  


     printf("\n");
     printf(" - Step 1b: Computing undamped induction energy with damped induced multipole moments.\n");
     // Since we are re-using the same induced multipole moments, we just
     // recompute the energy, this time using the undamped Tab matrices.

     // Now compute the energy contribution from the induced multipoles
     // Loop over all monomers in the central unit cell
     for (int imonA = 1; imonA <= NMon; imonA++) {
       int NatomsA = Monomers[imonA].GetNumberOfAtoms();
       //int offsetA = multipole_key1[imonA]; 
       // loop over all monomers that pair with this dimer - "inducers"
       //    These are central unit cell and image monomers
       for (int imonB = 1; imonB <= NMon+NMon_images; imonB++) {

	 // Classify the monomer as image or central unit cell
	 bool IsImageMonomer;
	 int ref_monB;
	 if (imonB > NMon) {
	   IsImageMonomer = true;
	   ref_monB = MonomerImages[imonB-NMon].GetReferenceMonomerIndex();
	 }
	 else {
	   IsImageMonomer = false;
	   ref_monB = Monomers[imonB].GetReferenceMonomerIndex();
	 }      

	 if (imonA != imonB) {
	   int NatomsB = Monomers[ref_monB].GetNumberOfAtoms();
	   //int offsetB = multipole_key1[ref_monB]; 
	   bool AB_order = true;
	   if (imonA > imonB) AB_order = false; 

	   // Grab the index of the appropriate dimer -- need this to
	   // obtain the Tab matrices
	   int idimer; 	  
	   Dimer tmp_dimer;
	   if ( !IsImageMonomer) {
	     idimer = DimerLookup(imonA,imonB); 

             //if seperation of dimers is greater than c0, its interaction Matrices do not exist
	     double seperation = Dimers[idimer].GetDimerSeparation();
             double c0 = Params::Parameters().GetLocalCutoff(0);
	     if(!Params::Parameters().DoLocal2BodyTruncation() || seperation >= c0)
	      Dimers[idimer].BuildTabInteractionMatrices();
	   }
	   else {
	     // Create the dimer
	     tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB-NMon]);
	     //tmp_dimer.BuildDampedTabInteractionMatrices();
	     tmp_dimer.BuildTabInteractionMatrices();
	   }

	   // Loop over atoms on each monomer
	   for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on inducee 
	     Multipole dQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments()); // by shuhao Mar. 2011

	     for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on inducers      

	       Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());

	       Matrix Tab;
	       if (AB_order) {
		 // If it's a real monomer, we just grab precomputed Tab
		 if (!IsImageMonomer) 

		   //Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		   Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
		 // Otherwise, we grab the Tab we just computed above
		 else {
		   //Tab.Initialize(tmp_dimer.GetDampedTabInteractionMatrix(iA,iB));
		   Tab.Initialize(tmp_dimer.GetTabInteractionMatrix(iA,iB));
		 }
		 int dimT1 = Tab.GetRows();
		 int dimT2 = Tab.GetCols();

		 for (int t=0;t<min(dQA.GetLength(),dimT1);t++) 
		   for (int u=0;u<min(QB.GetLength(),dimT2);u++) {
		     Eind1 += 0.5*dQA(t)*Tab(t,u)*QB(u);
		   }
	       }

	       else {
		 //Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		 Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
		 int dimT1 = Tab.GetRows();
		 int dimT2 = Tab.GetCols();
		 for (int t=0;t<min(dQA.GetLength(),dimT2);t++) 
		   for (int u=0;u<min(QB.GetLength(),dimT1);u++) {
		     Eind1 += 0.5*dQA(t)*Tab(u,t)*QB(u);
		   }	
	       } // end if (AB_order) ... else

	     } // end loop over inducer atoms
	   } // end loop over inducee atoms
	 } // end if (imonA != imonB)
       } // end loop over inducer monomers
     } // end loop over inducee monomers

     // Print out the induction energy just calculated based on the cluster originated by the polarization cutoff
     if (Params::Parameters().PrintLevel() >= 0) {
       printf("------------------------------------------------------------------------------\n");
       printf(" Undamped Induction Energy for the Cluster defined by the polarization cutoff\n");
       printf("------------------------------------------------------------------------------\n");
       printf("        Undamped Induction     = %12.6f kJ/mol\n",Eind1*HartreesToKJpermole);
       printf("------------------------------------------------------------------------------\n");
     }

     deltaa = Eind - Eind1;
     printf("\nDifference between damped and undamped short-range induction energies = %.3f kJ/mol\n\n",
	    deltaa*HartreesToKJpermole);

     fflush(stdout);
   }

   // if we skipped induction, we need to set the moments to zero to
   // make the Ewald code work.  This is a memory-wasting hack, but
   // given that we virtually always want to use induction, that's ok.
   else{
     for (int imon=1;imon<=NMon;imon++) {
       for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	 int Rmom = Monomers[imon].GetAtom(iA).GetMultipoleMoments().GetRank();
	 int Rpol = Monomers[imon].GetAtom(iA).GetPolarizability().GetRank();
	 int max_rank = max(Rmom,Rpol);

	 Multipole dQ; 
	 dQ.Initialize(max_rank);
	 Monomers[imon].GetAtom(iA).SetInduceMultipoleMoments(dQ);
       }
     }
   }


   time_t start_ewald_recip_time, stop_ewald_recip_time;
   time_t start_ewald_direc_time, stop_ewald_direc_time;
   time_t start_self_intra_time, stop_self_intra_time;


   start_ewald_recip_time = time(NULL);
   /* Step 1c: Compute the total permanent and induce multipole contributions of Lattice_MM.  */

   printf(" - Step 1c: Compute Ewald sum over permanent and induced multipole moments\n");
   double StaticERecip = 0.0;
   double InduceERecip1 = 0.0;
   double undampedInduceERecip1 = 0.0;

   double StaticERecip_kn0 = 0.0;
   double InduceERecip1_kn0 = 0.0;
   double undampedInduceERecip1_kn0 = 0.0;

   // Get damping factor for the induction ewald sumation 
   double beta_damp = Params::Parameters().GetDampingFactor(); // damping is also needed for the ewald sumation of induction energy


   // Decide whether we use manual Ewald summation parameters or
   // auto-determine them.  By default, nX, nY, nZ, kX, kY, kZ, &
   // kappa_param = -1.  If that's true, then we auto-determine them.
   // If not, then we use the manual values, in which case all 6
   // parameters & kappa need to be provided by the user.  Note: if
   // only kappa_param > 0, then we can still autodetermine the
   // cutoffs.
   int nX,nY,nZ,kX,kY,kZ;
   int sml_nX1, sml_nY1, sml_nZ1, sml_kX1, sml_kY1, sml_kZ1;
   int sml_nX2, sml_nY2, sml_nZ2, sml_kX2, sml_kY2, sml_kZ2;
   int sml_nX3, sml_nY3, sml_nZ3, sml_kX3, sml_kY3, sml_kZ3;
   double shift1,shift2,shift3;
   double kappa_param;

   bool auto_determine_ewald_parameters = true;

   // Grab the Ewald summation cutoffs: Direct space
   nX = Params::Parameters().GetDirecLoopX();
   nY = Params::Parameters().GetDirecLoopY();
   nZ = Params::Parameters().GetDirecLoopZ();

   // Grab the Ewald summation cutoffs: Reciprocal space
   kX = Params::Parameters().GetRecipLoopX();
   kY = Params::Parameters().GetRecipLoopY();
   kZ = Params::Parameters().GetRecipLoopZ();   


   if ( nX>0 || nY>0 || nZ>0 || kX>0 || kY>0 || kZ>0) {
     auto_determine_ewald_parameters = false;
   }

   if (auto_determine_ewald_parameters) {
     // GJB: Algorithm for determining summation cutoffs in Ewald sum.
     // Based on Frenkel & Smit, Ch. 12.1.5

     // Grab lengths of each real space unit cell vector
     double a = GetUnitCellParameter("a");
     double b = GetUnitCellParameter("b");
     double c = GetUnitCellParameter("c");

     double accuracy_fac = Params::Parameters().GetEwaldAccuracy();
     kappa_param = Params::Parameters().GetEwaldKappa();

     printf("  Ewald summation cutoffs determined automatically.  Accuracy factor = %f\n",accuracy_fac);

     if (kappa_param < 0.0) {
       // Auto-determine kappa.  Optimal efficiency seems to come from
       // summing over the same number of cells in real (nX,nY,nZ) and
       // reciprocal space (kX,kY,kZ).  One can show that if nX=kX, 
       // kappa = sqrt(pi)/a, where a is the lattice parameter.  Similar
       // equations exist for the Y and Z components.

       // Here, we take an average value of the optimal kappa for each of
       // X, Y, and Z.
       kappa_param = sqrt(pi)*(1.0/a + 1.0/b + 1.0/c)/3.0;
       // Store this value.  Right now, we have two unique copies due
       // to code modifications.  Will clean that up soon.  The EwaldKappa
       // parameter is the one we want to use.
       Params::Parameters().SetEwaldKappa(kappa_param);

       printf("    Using optimal Ewald kappa parameter = %f\n",kappa_param);
     }
     else {
       printf("    Using user-defined Ewald kappa parameter = %f\n",kappa_param);
     }

     // Determine real space cutoffs    
     double rc = accuracy_fac/kappa_param;
     nX = (int) ceil(rc/a);
     nY = (int) ceil(rc/b);
     nZ = (int) ceil(rc/c);

     // Determine reciprocal space cutoffs:
     kX = (int) ceil(accuracy_fac*a*kappa_param/pi);
     kY = (int) ceil(accuracy_fac*b*kappa_param/pi);
     kZ = (int) ceil(accuracy_fac*c*kappa_param/pi);

     printf("    (nX, nY, nZ) = (%d, %d, %d)    (kX, kY, kZ) = (%d, %d, %d)\n",nX,nY,nZ,kX,kY,kZ);

     // also determine cutoffs if we were to use a slightly smaller
     // accuracy factors.  We do this for 3 different smaller values
     // values of accuracy_fac: shift1, shift2, and shift3.  This data
     // gets used to gauge the convergence error in our Ewald sum.
     shift1 = accuracy_fac - 1.0;
     shift2 = accuracy_fac - 2.0;
     shift3 = accuracy_fac - 3.0;

     // put some lower bounds so we don't get negative 
     if (shift1 < 0.0) shift1 = 0.0;
     if (shift2 < 0.0) shift2 = 0.0;
     if (shift3 < 0.0) shift3 = 0.0;

     double sml_rc1 = shift1/kappa_param;
     sml_nX1 = (int) ceil(sml_rc1/a);
     sml_nY1 = (int) ceil(sml_rc1/b);
     sml_nZ1 = (int) ceil(sml_rc1/c);
     sml_kX1 = (int) ceil(shift1*a*kappa_param/pi);
     sml_kY1 = (int) ceil(shift1*b*kappa_param/pi);
     sml_kZ1 = (int) ceil(shift1*c*kappa_param/pi);

     double sml_rc2 = shift2/kappa_param;
     sml_nX2 = (int) ceil(sml_rc2/a);
     sml_nY2 = (int) ceil(sml_rc2/b);
     sml_nZ2 = (int) ceil(sml_rc2/c);
     sml_kX2 = (int) ceil(shift2*a*kappa_param/pi);
     sml_kY2 = (int) ceil(shift2*b*kappa_param/pi);
     sml_kZ2 = (int) ceil(shift2*c*kappa_param/pi);

     double sml_rc3 = shift3/kappa_param;
     sml_nX3 = (int) ceil(sml_rc3/a);
     sml_nY3 = (int) ceil(sml_rc3/b);
     sml_nZ3 = (int) ceil(sml_rc3/c);
     sml_kX3 = (int) ceil(shift3*a*kappa_param/pi);
     sml_kY3 = (int) ceil(shift3*b*kappa_param/pi);
     sml_kZ3 = (int) ceil(shift3*c*kappa_param/pi);

   }
   else {
     // Grab kappa and verify that all parameters have been specified
       bool error_flag = false;

       kappa_param = Params::Parameters().GetEwaldKappa();
       if (kappa_param < 0.0) error_flag = true;

       // Have the cutoffs been defined?
       if (nX<0 || nY<0 || nZ<0 || kX<0 || kY<0 || kZ<0) error_flag = true;

       // Exit if one or more parameters is undefined.
       if (error_flag) {
	 printf("\nEWALD ERROR: When using manual Ewald parameter specification, *all* parameters\n");
	 printf("             must be set explicitly in the input file.\n");
	 printf("   Current values (value of -1 => indicates that the parameter needs to be specified):\n");
	 printf("             (nX, nY, nZ) = (%d, %d, %d)\n",nX,nY,nZ);
	 printf("             (kX, kY, kZ) = (%d, %d, %d)\n",kX,kY,kZ);
	 printf("             EWALD_KAPPA = %f\n",kappa_param);
	 printf("   If you do not understand this message, you might want to use automated parameter\n");
	 printf("   determination associated with the EWALD_ACCURACY parameter instead.\n");
	 printf("\nExiting.\n\n");
	 exit(1);      
       }
       printf("  Ewald summation cutoffs set manually, kappa parameter = %f\n",kappa_param);
   }



   printf("  Begin the Ewald summation in reciprocal space\n");
   //printf("  --kX=%d,kY=%d,kZ=%d\n",kX,kY,kZ);

   fflush(stdout);

   // GJB: for testing convergence
   double Ees_recip_test1 = 0.0, Eind_recip_test1 = 0.0, Ees_real_test1 = 0.0, Eind_real_test1 = 0.0;
   double Ees_recip_test2 = 0.0, Eind_recip_test2 = 0.0, Ees_real_test2 = 0.0, Eind_real_test2 = 0.0;
   double Ees_recip_test3 = 0.0, Eind_recip_test3 = 0.0, Ees_real_test3 = 0.0, Eind_real_test3 = 0.0;

   // grab information for building the RecipTab; and calculate the
   // convergence factor alphaa used in the self-energy correction,
   // which is the same as that in RecipTab and DirecTab
   double CellV = cell_volume;
   double V = CellV*AngToBohr*AngToBohr*AngToBohr;

   Vector RecipCellx = reciprocal_cell[0];
   Vector RecipCelly = reciprocal_cell[1];
   Vector RecipCellz = reciprocal_cell[2];

   // Loop over all monomers in unit cell for A
   for (int imonA = 1; imonA <= NMon; imonA++) {
     int NatomsA = Monomers[imonA].GetNumberOfAtoms();
     Vector RotVecA(Monomers[imonA].GetRotationVector());
     double RotAngA = Monomers[imonA].GetRotationAngle();

     Matrix RotMatA = Monomers[imonA].GetRotationMatrix();

     // Loop over atoms on each monomer A
     for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
       Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
       Multipole InduceQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments()); // convergent induced multipole for the induction energy 
       int NmomA = QA.GetLength();
       int InduceNmomA = InduceQA.GetLength();

       // Loop over all monomers in unit cell for B
       for (int imonB = 1; imonB <= NMon; imonB++) {
	 int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	 Vector RotVecB(Monomers[imonB].GetRotationVector());
	 double RotAngB = Monomers[imonB].GetRotationAngle();

	 Matrix RotMatB = Monomers[imonB].GetRotationMatrix();

	 // Loop over atoms on each monomer B
	 for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	   Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	   Multipole InduceQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	   int NmomB = QB.GetLength();
	   int InduceNmomB = InduceQB.GetLength();

	   // Loop over reciprocal lattice vectors
	   for (int kx = -kX; kx<=kX; kx++){
	     for (int ky = -kY; ky<=kY; ky++){
	       for (int kz = -kZ; kz<=kZ; kz++){                  

		 if (kx*kx+ky*ky+kz*kz==0){ //if |kn|=0
		   // the case of Ltot = 2 has the finite limit for the
		   // case of |kn|=0, and there are two cases of Ltot
		   // =2: one is dipole-dipole; the other is
		   // charge-quadrupole;

		   // in our code we build the RecipTab_kn0 matrix with
		   // |kn|=0: for Ltot < 2, the |kn| = 0 terms will
		   // diverge, so we set the F0, F1 = 0;

		   // for Ltot > 2, the |Kn| = 0 term is zero in the
		   // process of building the RecipTab matrix;

		   // for Ltot = 2, please see the atom.C for the
		   // construction of RecipTab_kn0

		   // this |kn|=0 case for any molecular crystal
		   // introduces only very small energy (~ 10^-4
		   // kJ/mol), but in physics of Ewald summation, this
		   // term should be included

		   Matrix RecipTabAB_kn0;
		   RecipTabAB_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotMatA,Monomers[imonB].GetAtom(iB),
												RotMatB, kx, ky, kz, CellV,
												RecipCellx, RecipCelly,  RecipCellz,-999.0);

		   //RecipTabAB_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
		   //									       Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
		   //								       kx, ky, kz, CellV,
		   //								       RecipCellx, RecipCelly, RecipCellz, -999.0);
		   int dimTAB1_kn0 = RecipTabAB_kn0.GetRows();
		   int dimTAB2_kn0 = RecipTabAB_kn0.GetCols();
		   // Loop over multipole order
		   for (int t0=0;t0<min(NmomA,dimTAB1_kn0);t0++)  
		     for (int u0=0;u0<min(NmomB,dimTAB2_kn0);u0++) {
		       StaticERecip_kn0 += 0.50*QA(t0)*RecipTabAB_kn0(t0,u0)*QB(u0); // lattice electrostatic energy for permanent multipoles
		       undampedInduceERecip1_kn0 += 0.50*InduceQA(t0)*RecipTabAB_kn0(t0,u0)*QB(u0); // induction contribution, undamped 
		     } 
		 } // endif |kn|=0

		 else { // |kn|!=0
		   Matrix RecipTabAB;
		   //RecipTabAB = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
		   //								       Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
		   //							       kx, ky, kz, CellV,
		   //							       RecipCellx, RecipCelly, RecipCellz,-999.0);
		   RecipTabAB = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotMatA,Monomers[imonB].GetAtom(iB),
											RotMatB, kx, ky, kz, CellV,
											RecipCellx, RecipCelly,  RecipCellz,-999.0);

		   int dimTAB1 = RecipTabAB.GetRows();
		   int dimTAB2 = RecipTabAB.GetCols();
		   for (int t=0;t<min(NmomA,dimTAB1);t++)  
		     for (int u=0;u<min(NmomB,dimTAB2);u++) {

		       double dE_es = 0.50*QA(t)*RecipTabAB(t,u)*QB(u);
		       double dE_ind = 0.50*InduceQA(t)*RecipTabAB(t,u)*QB(u);

		       StaticERecip += dE_es; // lattice electrostatic energy for permanent multipoles
		       undampedInduceERecip1 += dE_ind; // induction contribution

		       // Also find contribution from outermost sets of k vectors to obtain
		       // a rough estimate of the convergence error in our Ewald sum.  GJB
		       if ( abs(kx) > sml_kX1 || abs(ky) > sml_kY1 || abs(kz) > sml_kZ1 
			    && auto_determine_ewald_parameters) {
			 Ees_recip_test1 += dE_es;
			 Eind_recip_test1 += dE_ind;
		       }
		       if ( abs(kx) > sml_kX2 || abs(ky) > sml_kY2 || abs(kz) > sml_kZ2 
			    && auto_determine_ewald_parameters) {
			 Ees_recip_test2 += dE_es;
			 Eind_recip_test2 += dE_ind;
		       }
		       if ( abs(kx) > sml_kX3 || abs(ky) > sml_kY3 || abs(kz) > sml_kZ3 
			    && auto_determine_ewald_parameters) {
			 Ees_recip_test3 += dE_es;
			 Eind_recip_test3 += dE_ind;
		       }
		     }   

		 } // end else (|kn|!=0)

	       } // end loop over kx
	     } // end loop over ky
	   } // end loop over kz
	 } // end loop over atoms on monomerA
       } // end loop over atoms on monomerB
     } // end loop over A monomers
   } // end loop over B monomers

   stop_ewald_recip_time = time(NULL);
   double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
   printf("  Time to ewald summation in reciprocal space = %.2f seconds\n",ewald_recip_time);

   // Check convergence in reciprocal space
   double Recip_conv1 = Ees_recip_test1 + Eind_recip_test1;
   double Recip_conv2 = Ees_recip_test2 + Eind_recip_test2;
   double Recip_conv3 = Ees_recip_test3 + Eind_recip_test3;

   start_ewald_direc_time = time(NULL);


   // calculate the Ewald energy in the direct space 
   double StaticEDirec = 0.0;
   double InduceEDirec1 = 0.0;
   double undampedInduceEDirec1 = 0.0;

   printf("  Begin the Ewald summation in direct space\n");  
   //printf("  --nX=%d,nY=%d,nZ=%d\n",nX,nY,nZ);

   Vector UnitCellx = unit_cell[0];
   Vector UnitCelly = unit_cell[1];
   Vector UnitCellz = unit_cell[2];


   for (int imonA = 1; imonA <= NMon; imonA++) {
     int NatomsA = Monomers[imonA].GetNumberOfAtoms();
     //Vector RotVecA(Monomers[imonA].GetRotationVector());
     //double RotAngA = Monomers[imonA].GetRotationAngle();

     Matrix RotMatA = Monomers[imonA].GetRotationMatrix();

       for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
	 Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // should be the total multipole moments of induce and permanant multipoe
	 Multipole InduceQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments());
	 int NmomA = QA.GetLength();
	 int InduceNmomA = InduceQA.GetLength();

	 for (int imonB = 1; imonB <= NMon; imonB++) {
	   int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	   //Vector RotVecB(Monomers[imonB].GetRotationVector());
	   //double RotAngB = Monomers[imonB].GetRotationAngle();

	   Matrix RotMatB = Monomers[imonB].GetRotationMatrix();

	   for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	     Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	     Multipole InduceQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	     int NmomB = QB.GetLength();
	     int InduceNmomB = InduceQB.GetLength();

	     for (int nx = -nX; nx<=nX; nx++){
	       for (int ny = -nY; ny<=nY; ny++){
		 for (int nz = -nZ; nz<=nZ; nz++){

		   if (nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB);
		   //{               
		   //printf("Skipping self-interaction real-space term\n");
		   //} 

		   else { // |rAB-rn|!=0
		     Matrix DirecTabAB;
		     //yoni:using rotation matrix instant of rotation angle/vector
		     //DirecTabAB = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
		     // 								 Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
		     //									 nx, ny, nz, CellV,
		     //								 UnitCellx, UnitCelly, UnitCellz, -999.0);
		     DirecTabAB = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotMatA, Monomers[imonB].GetAtom(iB),RotMatB,
											  nx,ny,nz,CellV,UnitCellx, UnitCelly, UnitCellz, -999.0);
		     int DdimTAB1 = DirecTabAB.GetRows();
		     int DdimTAB2 = DirecTabAB.GetCols();

		     for (int t5=0;t5<min(NmomA,DdimTAB1);t5++) 
		       for (int u5=0;u5<min(NmomB,DdimTAB2);u5++) {

			 double dE_es = 0.50*QA(t5)*DirecTabAB(t5,u5)*QB(u5);
			 double dE_ind = 0.50*InduceQA(t5)*DirecTabAB(t5,u5)*QB(u5);

			 StaticEDirec += dE_es;
			 undampedInduceEDirec1 += dE_ind;

			 // Also find contribution from outermost sets of k vectors to obtain
			 // a rough estimate of the convergence error in our Ewald sum.  GJB
			 if ( abs(nx) > sml_nX1 || abs(ny) > sml_nY1 || abs(nz) > sml_nZ1 
			      && auto_determine_ewald_parameters) {
			   Ees_real_test1 += dE_es;
			   Eind_real_test1 += dE_ind;
			 }
			 if ( abs(nx) > sml_nX2 || abs(ny) > sml_nY2 || abs(nz) > sml_nZ2
			      && auto_determine_ewald_parameters) {
			   Ees_real_test2 += dE_es;
			   Eind_real_test2 += dE_ind;
			 }
			 if ( abs(nx) > sml_nX3 || abs(ny) > sml_nY3 || abs(nz) > sml_nZ3 
			      && auto_determine_ewald_parameters) {
			   Ees_real_test3 += dE_es;
			   Eind_real_test3 += dE_ind;
			 }
		       }

		   } // end else |rAB-rn|!=0
		 } // end loop nx
	       } // end loop ny
	     } // end loop nz
	   } // end loop over atoms monomerA
	 } // end loop over atoms monomerB
       } // end loop over first monomers
     } // end loop over second monomers  

   stop_ewald_direc_time = time(NULL);
   double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
   printf("  Time to ewald summation in direct space = %.2f seconds\n",ewald_direc_time);


   // Check convergence in real space space
   double Real_conv1 = (Ees_real_test1 + Eind_real_test1);
   double Real_conv2 = (Ees_real_test2 + Eind_real_test2);
   double Real_conv3 = (Ees_real_test3 + Eind_real_test3);

   // self-interaction energy corection and intramolecular correction

   start_self_intra_time = time(NULL);

   printf("  Begin the self-interaction and intramolecular energy correction\n");
   double Uself = 0.0; 
   double Uself_induce1 = 0.0;
   double Uself_Electrostatic = 0.0;
   double UIntramolecule = 0.0;
   double UIntramolecule_induce1 = 0.0;
   double undampedUIntramolecule_induce1 = 0.0;
   double UIntramolecule_Electrostatic = 0.0;


   // constant for the calculation Uself
   const double cccc = 1.0/3.0;
   double V_V = pow(V,cccc);
   double kappa = kappa_param/V_V;
   double alphaa = kappa*kappa;
   double ssff = kappa/sqrt(3.14159265359);
   printf("  -----alphaa = %12.6f, kappa = %12.6f, ssff = %12.6f\n", alphaa,kappa,ssff);

   double Perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
   //printf("--Perm = %12.6f\n", Perm);



   /*
  // just for test the E-self using another method, has commented out
  for (int Imon = 1; Imon <= NMon; Imon++) {
  int NAtoms = Monomers[Imon].GetNumberOfAtoms();
  for (int IA=0;IA<NAtoms;IA++)  {// loop over atoms on monomer
  Multipole QA1(Monomers[Imon].GetAtom(IA).GetMultipoleMoments()); // should be the total multipole moments of induce and permanant multipoe
  Multipole InduceQA1(Monomers[Imon].GetAtom(IA).GetInduceMultipoleMoments());
  Uselff = Uselff -ssff*(QA1(0)+InduceQA1(0))*(QA1(0)+InduceQA1(0))/Perm;
  ///printf("---Uselff = %12.6f kJ/mol\n",Uselff*HartreesToKJpermole);
  }
  }
  */
  
  
  for (int Imon = 1; Imon <= NMon; Imon++) {
    int NAtoms = Monomers[Imon].GetNumberOfAtoms();

    //Vector RotVec =(Monomers[Imon].GetRotationVector());
    //double RotAng = Monomers[Imon].GetRotationAngle();
    Matrix RotationMatrix = Monomers[Imon].GetRotationMatrix();

    for (int IA=0;IA<NAtoms;IA++)  {// loop over atoms on monomer 
      Multipole QA1(Monomers[Imon].GetAtom(IA).GetMultipoleMoments()); 
      Multipole InduceQA1(Monomers[Imon].GetAtom(IA).GetInduceMultipoleMoments());
      int NmomA1 = QA1.GetLength();
      int InduceNmomA1 = InduceQA1.GetLength();
      for (int IB=0;IB<NAtoms;IB++)  {// loop over atoms on monomer 
	Multipole QB1(Monomers[Imon].GetAtom(IB).GetMultipoleMoments()); 
	Multipole InduceQB1(Monomers[Imon].GetAtom(IB).GetInduceMultipoleMoments());
	int NmomB1 = QB1.GetLength();
	int InduceNmomB1 = InduceQB1.GetLength();
	
	// if IA==IB, do the self-energy correction
	if (IA==IB){
	  //printf("----calculate the self-energy----\n");
	  Uself_Electrostatic -= ssff*QA1(0)*QB1(0)/Perm;
	  Uself_induce1 -= ssff*InduceQA1(0)*QB1(0)/Perm;
	} 
	// else do the intramolecular-energy correction for the case that A and B are atoms in the same molecule in the central cell 
	else {
	  //printf("----calculate the intramolecular-energy----\n");
	  
	  
	  // Vector RotVecA = RotVec;
	  //double RotAngA = RotAng;
	  //Vector RotVecB = RotVec;
	  //double RotAngB = RotAng;
	  
	  Matrix RotationMatrixA = RotationMatrix;
	  Matrix RotationMatrixB = RotationMatrix;
	  Matrix DirecTabA1B1;
	  //Matrix DampedDirecTabA1B1; 
	  
          
	  //this interaction Tab matrix is for the case where A and B
	  //are atoms in the same molecule the explicit form of the
	  //interaction energy between A and B, so we use
	  //BuildDirecInteractionMatrix (not reciprocal or direct
	  //matrix)

	  DirecTabA1B1 = Monomers[Imon].GetAtom(IA).BuildInteractionMatrix(RotationMatrixA,
				 Monomers[Imon].GetAtom(IB),
				 RotationMatrixB, -999.0);

	  //DirecTabA1B1 = Monomers[Imon].GetAtom(IA).BuildInteractionMatrix(RotVecA, RotAngA,
	  //                                                                 Monomers[Imon].GetAtom(IB), RotVecB, RotAngB,-999.0);
	  
	  //DampedDirecTabA1B1 = Monomers[Imon].GetAtom(IA).BuildInteractionMatrix(RotVecA, RotAngA,
	  //                                                               Monomers[Imon].GetAtom(IB), RotVecB, RotAngB, beta_damp);
	  
	  
	  int DdimTA1B1 = DirecTabA1B1.GetRows();
	  int DdimTA1B2 = DirecTabA1B1.GetCols();
	  
	  for (int tt5=0;tt5<min(NmomA1,DdimTA1B1);tt5++) // for permanent multipoles
	    for (int uu5=0;uu5<min(NmomB1,DdimTA1B2);uu5++) {
	      UIntramolecule_Electrostatic -= 0.50*QA1(tt5)*DirecTabA1B1(tt5,uu5)*QB1(uu5);
	      //UIntramolecule_induce1 -= 0.50*InduceQA1(tt5)*DampedDirecTabA1B1(tt5,uu5)*QB1(uu5);
	      undampedUIntramolecule_induce1 -= 0.50*InduceQA1(tt5)*DirecTabA1B1(tt5,uu5)*QB1(uu5);
	    }
	} // end else
      } // end loop atoms A1
    } // end loop atoms  A2
  } // end loop monomers in the central cell
  
    //Uself = Uself_Electrostatic + Uself_induce1;
    //UIntramolecule = UIntramolecule_Electrostatic + UIntramolecule_induce1;
  
  stop_self_intra_time = time(NULL);
  double self_intra_time = difftime(stop_self_intra_time, start_self_intra_time);
  printf("  Time for self-interaction and intramolecular energy correction = %.0f seconds\n", self_intra_time);
  
  
  // print all values for debug
  printf("\n");
  printf("Summary of Electrostatic & Induction Energies\n");
  printf("--------------------------------------------------------\n");
  printf(" Permanent Multipole contribution:\n");
  printf("--------------------------------------------------------\n");
  printf("     Static_EDirect     = %12.4f kJ/mol\n",StaticEDirec*HartreesToKJpermole);
  printf("     Static_ERecip      = %12.4f kJ/mol\n",StaticERecip*HartreesToKJpermole);
  printf("     Static_ERecipkn0   = %12.4f kJ/mol\n",StaticERecip_kn0*HartreesToKJpermole);
  printf("     Static_UIntramol   = %12.4f kJ/mol\n",UIntramolecule_Electrostatic*HartreesToKJpermole);
  printf("     Static_Uself       = %12.4f kJ/mol\n",Uself_Electrostatic*HartreesToKJpermole);
  // printf("--------------------------------------------------------\n");
  // printf("     Induce_EDirec      = %12.4f kJ/mol\n",InduceEDirec1*HartreesToKJpermole);
  // printf("     Induce_ERecip      = %12.4f kJ/mol\n",InduceERecip1*HartreesToKJpermole);
  // printf("     Induce_ERecipkn0   = %12.4f kJ/mol\n",InduceERecip1_kn0*HartreesToKJpermole);
  // printf("     Induce_UIntramol   = %12.4f kJ/mol\n",UIntramolecule_induce1*HartreesToKJpermole);
  // printf("     Induce_Uself       = %12.4f kJ/mol\n",Uself_induce1*HartreesToKJpermole);
  
  printf("--------------------------------------------------------\n");
  printf(" Induction contribution (undamped):\n");
  printf("--------------------------------------------------------\n");
  printf("     Induce_EDirect     = %12.4f kJ/mol\n",undampedInduceEDirec1*HartreesToKJpermole);
  printf("     Induce_ERecip      = %12.4f kJ/mol\n",undampedInduceERecip1*HartreesToKJpermole);
  printf("     Induce_ERecipkn0   = %12.4f kJ/mol\n",undampedInduceERecip1_kn0*HartreesToKJpermole);
  printf("     Induce_UIntramol   = %12.4f kJ/mol\n",undampedUIntramolecule_induce1*HartreesToKJpermole);
  printf("     Induce_Uself       = %12.4f kJ/mol\n",Uself_induce1*HartreesToKJpermole);
  // printf("--------------------------------------------------------\n");
  // printf("     Total_ERecip       = %12.4f kJ/mol\n",Ewald_Recip*HartreesToKJpermole);
  // printf("     Total_EDirec       = %12.4f kJ/mol\n",Ewald_Direct*HartreesToKJpermole);
  // printf("     Total_UIntramol    = %12.4f kJ/mol\n",UIntramolecule*HartreesToKJpermole);
  // printf("     Total_Uself        = %12.4f kJ/mol\n",Uself*HartreesToKJpermole);
  
  
  
  // Store the results - in hartrees.
  double E_E_MM = StaticEDirec + StaticERecip + UIntramolecule_Electrostatic + Uself_Electrostatic;
  double E_E_MM_kn0 = StaticEDirec + StaticERecip + UIntramolecule_Electrostatic + Uself_Electrostatic + StaticERecip_kn0;
  //double E_I_MM = InduceEDirec1  + InduceERecip1    +  Uself_induce1  + UIntramolecule_induce1;
  //double E_I_MM_kn0 = InduceEDirec1  + InduceERecip1  +  InduceERecip1_kn0  +  Uself_induce1  + UIntramolecule_induce1;
  double undamped_ewald = undampedInduceEDirec1  + undampedInduceERecip1  +  undampedInduceERecip1_kn0  +  Uself_induce1  +  undampedUIntramolecule_induce1;
  
  Lattice_E_Electrostatic_MM = E_E_MM_kn0;
  //double Lattice_E_Induction_MM_old = E_I_MM_kn0;
  Lattice_E_Induction_MM = undamped_ewald + deltaa;
  
  //double Lattice_E_Electrostatic_MM_nokn0 = E_E_MM;
  //double Lattice_E_Induction_MM_nokn0 = E_I_MM;
  
  
  
    //double Etot1 = Lattice_E_Electrostatic_MM_nokn0 + Lattice_E_Induction_MM_nokn0;
  
  double Etot = Lattice_E_Electrostatic_MM + Lattice_E_Induction_MM;
  
  // Print out a summary of the classical energy of a crystal just calculated  
  if (Params::Parameters().PrintLevel() >= 0) {    
    printf("------------------------------------------------------------\n");    
    printf(" Total classical electrostatic + induction Energy\n");    
    printf("------------------------------------------------------------\n");    
    printf("     E_ewald(undamped)            = %12.4f kJ/mol\n",undamped_ewald*HartreesToKJpermole);
    printf("     E_SR(damped)-E_SR(undamped)  = %12.4f kJ/mol\n",deltaa*HartreesToKJpermole);
    printf("     E_ewald(undamped) + E_SR(damped)-E_SR(undamped) = %12.4f kJ/mol\n",(undamped_ewald+deltaa)*HartreesToKJpermole);
    printf("     induction energy in the polarization cut-off using damping i.e. E_SR(damped) =  %12.4f kJ/mol\n", Eind*HartreesToKJpermole);
    printf("     Lattice_Electrostatic        = %12.4f kJ/mol\n",Lattice_E_Electrostatic_MM*HartreesToKJpermole);
    printf("     Lattice_Induction            = %12.4f kJ/mol\n",Lattice_E_Induction_MM*HartreesToKJpermole);
    printf("   ---------------------------------------\n");    
    printf("     Total                        = %12.4f kJ/mol\n",Etot*HartreesToKJpermole);    
    // printf("     Lattice_Total_nokn0          = %12.4f kJ/mol\n",Etot1*HartreesToKJpermole);
    // printf("     estimated Ewald error        = %12.1e kJ/mol\n",(fabs(Recip_conv)+fabs(Real_conv))*HartreesToKJpermole);    
    printf("--------------------------------------------------------\n");  } 
  


  // Print out some error analysis for the Ewald sum convergence when
  // using the automatic Ewald parameter determination.  Note, the reason we print out
  // several different values is that the Ewald sum convergence seems to exhibit
  // artificial plateaus for small ranges of EWALD_ACCURACY.  Looking at multiple
  // points allows us to get a better sense of how it's behaving.  GJB
  if (Params::Parameters().PrintLevel() >= 0 && auto_determine_ewald_parameters ) {   
    double accuracy_fac = Params::Parameters().GetEwaldAccuracy();
    
    printf("\n Analysis of Ewald sum convergence (kJ/mol):\n");
    printf(" ----------------------------------------------\n");
    printf("   Ewald Accuracy         dE(Ewald) kJ/mol\n");
    printf(" ----------------------------------------------\n");
    printf("      %4.1f              %12.4f\n",shift3,(fabs(Recip_conv3)+fabs(Real_conv3))*HartreesToKJpermole);
    printf("      %4.1f              %12.4f\n",shift2,(fabs(Recip_conv2)+fabs(Real_conv2))*HartreesToKJpermole);
    printf("      %4.1f              %12.4f\n",shift1,(fabs(Recip_conv1)+fabs(Real_conv1))*HartreesToKJpermole);
    printf(" ----------------------------------------------\n");
    printf("      %4.1f (current)    %12.4f\n",accuracy_fac, Etot*HartreesToKJpermole);
    printf(" ----------------------------------------------\n");
    printf(" Note: In general, all 3 dE values should exhibit acceptable convergence.\n\n");
  }


  return Etot;
}

void Cluster::ComputeLatticeParamGrad(Vector& Grad) {

  /*
  if(!Params::Parameters().UseCartesianCoordinates()){

    Grad.PrintGradient("Grad Bohr");
  //converting grad into fractional coordinates
    Grad.Scale(1/BohrToAng);

    Grad.PrintGradient("Grad Ang");
    Grad = ConvertBetweenFractionAndCartesianCoordinates(Grad,1,1,1);

    Grad.PrintGradient("Grad Fractional");
  }
  */

  //cell params are frozen
  if(Params::Parameters().FreezeUnitCellParams()){

    for(int i = 0; i < 6; i++){
      Grad[3*UniqueAtoms + i] = 0;
    }
  }
  else{

    //first calculate the gradient due to change in lattice params for the full system (central cell) at MM
    Vector Grad_Lat_Full_MM(6);
    
    Vector full_lat_plus(6), full_lat_minus(6);
    bool neglect_mb = Params::Parameters().NeglectManyBody();
    double dE_dx1=0, dE_dx2=0, dE_dx3=0, dE_dy1=0, dE_dy2=0, dE_dy3=0, dE_dz1=0, dE_dz2=0, dE_dz3=0;

    //monomer symmetry untilized
    bool UseMonSym = Params::Parameters().UseMonomerSymmetry();
    if ( !neglect_mb ) {  
      
      /*
      if(Params::Parameters().GetMMType() == 5){//Read Crystal Lattice Gradient
       
        string path = Params::Parameters().GetMMPath();
   
        //Path of the quasiharmonic calculations
        if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
          path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

        // Set up the filename, with the full path.  File is 'full.out'
        string filename;
        filename = path + "/full.out";
    
        // Open the force file
        ifstream infile;
        infile.open( filename.c_str() );
        if ( !infile.is_open() ) {
          printf("Cluster::ComputeLatticeParamGrad : Cannot open file '%s'\n",
	        filename.c_str());
          exit(1);
         }
         while ( !infile.eof() ) {
           string line;
           getline(infile,line);
           if(line == "GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR"){
             getline(infile,line);//discard first line 
             getline(infile,line);//discard second line 
             getline(infile,line);//discard third line 
             
             //first vector
             getline(infile,line);
             istringstream iss(line);
             iss >> dE_dx1;
             iss >> dE_dy1;
             iss >> dE_dz1;

             //second vector
             getline(infile,line);
             istringstream iss2(line);
             iss2 >> dE_dx2;
             iss2 >> dE_dy2;
             iss2 >> dE_dz2;

             //third vector
             getline(infile,line);
             istringstream iss3(line);
             iss3 >> dE_dx3;
             iss3 >> dE_dy3;
             iss3 >> dE_dz3;
             break;
           }
         }
         infile.close();  // Close the file
      }else{*/
         
        for (int i=0;i<6;i++) {                 
	  full_lat_plus[i] = New_ComputeLatticeParamGradFullMM(i,-0.001);
        }
        for (int i=0;i<6;i++) {
	  full_lat_minus[i] = New_ComputeLatticeParamGradFullMM(i,0.001);
        }               
        Grad_Lat_Full_MM = full_lat_minus;
        Grad_Lat_Full_MM -= full_lat_plus;
      
        for (int i=0;i<3;i++) {
	  Grad_Lat_Full_MM[i] *= 1.0 / (2.0*0.001*AngToBohr);
        }
        for (int i=3;i<6;i++) {   
	  Grad_Lat_Full_MM[i] *= 1.0 / (2.0*0.01*DegreesToRadians);
        }
      //}
    }                 
    else {
      for (int i=0;i<6;i++) {
        Grad_Lat_Full_MM[i] = 0.0;             
       }
    }
        

    cout <<  " Grad_Lat_Full_MM[0] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[0] << "\n";
    cout <<  " Grad_Lat_Full_MM[1] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[1] << "\n";
    cout <<  " Grad_Lat_Full_MM[2] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[2] << "\n";
    cout <<  " Grad_Lat_Full_MM[3] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[3] << "\n";
    cout <<  " Grad_Lat_Full_MM[4] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[4] << "\n";
    cout <<  " Grad_Lat_Full_MM[5] = " << setprecision(11) << showpoint << Grad_Lat_Full_MM[5] << "\n";


    int Ntot = GetTotalNumberOfAtoms();
    //int *tmp_k_vec = NULL;
    double c0 = Params::Parameters().GetLocalCutoff(0);
    double c1 = Params::Parameters().GetLocalCutoff(1); 
    
    double a = UnitCellAxes[0];
    double b = UnitCellAxes[1];
    double c = UnitCellAxes[2];
    
    double alpha = UnitCellAngles[0];          
    double beta = UnitCellAngles[1];
    double gamma = UnitCellAngles[2]; 
    
    a *= AngToBohr;
    b *= AngToBohr;
    c *= AngToBohr;          
    alpha *= DegreesToRadians;
    beta *= DegreesToRadians;
    gamma *= DegreesToRadians;  
    
    double x1,x2,x3,y1,y2,y3,z1,z2,z3,beta_term,gamma_term, v;
    
    beta_term = ( cos(alpha) - cos(beta)*cos(gamma) ) / sin(gamma) ;
    gamma_term = sqrt( 1 - cos(beta)*cos(beta) - beta_term*beta_term ) ;
    v = gamma_term*sin(gamma);
    
    x1 = a;
    x2 = b*cos(gamma);
    x3 = c*cos(beta);  
    y1 = 0.0;
    y2 = b*sin(gamma);       
    y3 = c*beta_term;                        
    z1 = 0.0;
    z2 = 0.0;
    z3 = c*gamma_term;

    double nmb_scafac;
    if ( neglect_mb ) {
      nmb_scafac = 0.0;
    }
    else {
      nmb_scafac = 1.0;
    }
    

    //In some highly symmetrical crystals, atoms are shifted preserve symmetry
    //monomer contribution
    for(int i=1; i<=NMon;i++){

      for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){

	//The atom gradients have be rotated for the purposes of the nuclear gradient, we are reversing the rotation here
	Vector QM_Grad(3);
	Vector MM_Grad(3);
	for(int xyz=0;xyz<3;xyz++){
	  QM_Grad[xyz] = Monomers[i].GetQMGradient()[3*j+xyz];
	  MM_Grad[xyz] = Monomers[i].GetMMGradient()[3*j+xyz];
	}
	//Matrix Rot = Monomers[i].GetRotationMatrix();
	//Rot.Transpose();
	//QM_Grad = Rot.MatrixTimesVector(QM_Grad);
	//MM_Grad = Rot.MatrixTimesVector(MM_Grad);

	//Atom Fractional Coordinate
	//Vector Fract = Monomers[i].GetAtom(j).GetFractionalPosition();
	//Translation vector of the space group operator
	Vector Fract = Monomers[i].GetAtom(j).GetFractTranslationVector(true);
        for(int k = 0;k<3;k++){
	  Fract[k]  -= Monomers[i].GetFractTranslationVector(true)[k];
          Fract[k] += Monomers[i].GetAtom(j).GetShiftVector()[k];
	}	

	//x vector
	dE_dx1 += Fract[0] * (QM_Grad[0] - MM_Grad[0]);
	dE_dx2 += Fract[1] * (QM_Grad[0] - MM_Grad[0]);
	dE_dx3 += Fract[2] * (QM_Grad[0] - MM_Grad[0]);

	//y vector
	dE_dy1 += Fract[0] * (QM_Grad[1] - MM_Grad[1]);
	dE_dy2 += Fract[1] * (QM_Grad[1] - MM_Grad[1]);
	dE_dy3 += Fract[2] * (QM_Grad[1] - MM_Grad[1]);
	
	//z vector
	dE_dz1 += Fract[0] * (QM_Grad[2] - MM_Grad[2]);
	dE_dz2 += Fract[1] * (QM_Grad[2] - MM_Grad[2]);
	dE_dz3 += Fract[2] * (QM_Grad[2] - MM_Grad[2]);	
      }//end of loop over atoms
    }//end of lop over monomers



    
    //double dkbose;
    //for (int i=1;i<=3;i++) {

    //adding contribution of the non-image and image dimers symmetrical to the non-image dimers
    for (int i=1; i<=NDim;i++){

      if(Dimers[i].GetSymmetryFactor()!=0){
	int indexA = Dimers[i].GetIndexA(); 
	int indexB = Dimers[i].GetIndexB();
	int Na = Monomers[indexA].GetNumberOfAtoms(); 
	int Nb = Monomers[indexB].GetNumberOfAtoms();
	double damp = Dimers[i].GetDampingFactor(c0,c1);
	int ListSize = ( Dimers[i].GetSymmetryList().size() + Dimers[i].GetPeriodicSymmetryList().size())/2;


	for(int j=0; j<ListSize;j++){


	  //These are the indices of the monomer on the SymmetryList and PeriodicSymmetryList 
	  int SymA; 
	  int SymB;	  
	  
	  //determining which monomer in dimer is outside unit cell
	  //if both monomers are in unit cell (for a non-image dimer), WhichMon is -1
	  //int WhichMonB;
	  int k_vec[3]; 
	  //The symmetrical equivalent dimers that are periodic images have their contrubution to the gradient scaled by 1/2
	  float per_sca;

          //Rotation Matrix for the Dimers on the symmetry List
	  Matrix Rot;

	  //List of symmetrical atoms on the dimers on the Symmetry Dimer List
	  Vector SymAtomList;
		
	  if(j<Dimers[i].GetSymmetryList().size()/2){
	    SymA = Dimers[i].GetSymmetryList()[2*j];
	    SymB = Dimers[i].GetSymmetryList()[2*j+1];
	    k_vec[0] = 0;
	    k_vec[1] = 0;
	    k_vec[2] = 0;
	    //WhichMonB = -1;
	    per_sca = 1;
	    Rot = Dimers[i].GetRotationList()[j];
	    Rot.Transpose();
	    SymAtomList = Dimers[i].GetAtomEquivalency()[j];
	  }
	  else{
	    int k = j - Dimers[i].GetSymmetryList().size()/2;
	    SymA = Dimers[i].GetPeriodicSymmetryList()[2*k];
	    SymB = Dimers[i].GetPeriodicSymmetryList()[2*k+1];

	    k_vec[0] = Dimers[i].GetSymmetricalImageCell()[3*k];
	    k_vec[1] = Dimers[i].GetSymmetricalImageCell()[3*k+1];
	    k_vec[2] = Dimers[i].GetSymmetricalImageCell()[3*k+2];
	    per_sca = 0.5;
	    Rot = Dimers[i].GetPeriodicRotationList()[k];
	    Rot.Transpose();
	    SymAtomList = Dimers[i].GetPeriodicAtomEquivalency()[k];
	    //if(Dimers[i].GetMonBList()[k])
	    //  WhichMonB  =  1;
	    //else
	    //  WhichMonB = 0;
	  }

	  //first print statement, print elements before adding to it.
	  //printf("d(%i,%i)\n",SymA,SymB);
	  //printf("using d(%i,%i)\n",indexA,indexB);
	  //printf("before dE_dx1 = %f\n",dE_dx1);

	  //if(WhichMonB==0){//first atom in the symmetrical equvalient dimer is outside unit cell
	  //  int index = Dimers[i].GetPeriodicSymmetryList()[2*j];
	  //  mon = Monomers[index];
	    //  Startgrad = 0;
	    //  endloop = Na;
	  //}else if(WhichMonB==1){
	  //  int index = Dimers[i].GetPeriodicSymmetryList()[2*j+1];
	  //   mB = Monomers[index];
	  //  Startgrad = 3*Na;
	  //  endloop = Nb;
	  // }
	  //Matrix Rot = mB.GetRotationMatrix();
	  //Rot.Transpose();


	  //printf("d(%i,%i)\n",SymA,SymB);
	  //printf("using d(%i,%i)\n",indexA,indexB);
	  //printf("before dE_dz2 = %f\n",dE_dz2);

	  //SymAtomList.Print("SymAtomList");

	  //determines if SymA is symmetrical to MonA or MonB by determining
	  //Need to know if MonA and MonB have different number to atoms
	  bool FirstMonA = false;
	  for(int atom = 0; atom < Monomers[indexA].GetNumberOfAtoms();atom++){
	    if(SymAtomList[atom] == 1) FirstMonA = true;
	  }


	  for(int iatom=0; iatom<Na+Nb; iatom++){

	    int mon;
	    //int Startgrad;
	    //int LastAtom;
	    int tmp_k_vec[3];

	    //Fractional Translation
	    Vector Fract; 

	    //Vector Coord;
	    //when iatom is on monA
	    if(iatom<Na){
	      //mon = SymA;
	      //Startgrad = 0;
	      //LastAtom = Na;
	      //mon = SymA;
	      //Fract = Monomers[mon].GetAtom(iatom).GetFractionalPosition();
	      //Coord = Monomers[indexA].GetAtom(iatom).GetPosition();

	      /*
              if(WhichMonB==0){
		tmp_k_vec[0] = k_vec[0];
		tmp_k_vec[1] = k_vec[1];
		tmp_k_vec[2] = k_vec[2];
	      }else{
		tmp_k_vec[0] = 0;
		tmp_k_vec[1] = 0;
		tmp_k_vec[2] = 0;
	      }
              */
	    }
	    //when iatom is on monB
	    else{
	      //mon = SymB;
	      //Startgrad = 3*Na;
	      //LastAtom = Nb;
	      //mon = SymB;
	      //Fract = Monomers[mon].GetAtom(iatom-Na).GetFractionalPosition();
	      //Coord = Monomers[indexB].GetAtom(iatom-Na).GetPosition();

	      /*
	      if(WhichMonB==1){
		tmp_k_vec[0] = k_vec[0];
		tmp_k_vec[1] = k_vec[1];
		tmp_k_vec[2] = k_vec[2];
	      }else{
		tmp_k_vec[0] = 0;
		tmp_k_vec[1] = 0;
		tmp_k_vec[2] = 0;
	      }
              */
	    }
	    


	    int DimerAtom = (int) SymAtomList[iatom] - 1;
	    //The atom index on the monomer
	    int Symi;

	    //Determining atom index on monomer
	    if(FirstMonA){
	      if(DimerAtom < Na ){
		Symi = DimerAtom;  
		mon = SymA;
		tmp_k_vec[0] = 0;
		tmp_k_vec[1] = 0;
		tmp_k_vec[2] = 0;
	      }else{
		Symi = DimerAtom - Na;
		mon = SymB;
		tmp_k_vec[0] = k_vec[0];
		tmp_k_vec[1] = k_vec[1];
		tmp_k_vec[2] = k_vec[2];
	      }
	    }else{
	      if(DimerAtom < Nb ){
		Symi = DimerAtom;  
		mon = SymA;
		tmp_k_vec[0] = 0;
		tmp_k_vec[1] = 0;
		tmp_k_vec[2] = 0;
	      }
	      else{
		Symi = DimerAtom - Nb;
		mon = SymB;
		tmp_k_vec[0] = k_vec[0];
		tmp_k_vec[1] = k_vec[1];
		tmp_k_vec[2] = k_vec[2];
	      }
	    }

	    //Translation vector of the space group operator
	    Fract = Monomers[mon].GetAtom(Symi).GetFractTranslationVector(true);
            for(int j = 0;j<3;j++){
              Fract[j] += Monomers[mon].GetAtom(Symi).GetShiftVector()[j];
	    }	
             /*
	    //if monomer symmetry is employed, using gradient of the symmetrical unique monomer instead
	    if(UseMonSym)
	      mon = Monomers[mon].GetSymmetricalMonomer();

            
	    //find which atom in the dimers on the symmetry list is symmetrical to in the unique dimer
	    for(int k = 0; k<LastAtom;k++){
	      if(iatom<Na && Symi == -1){
		if(Monomers[mon].GetAtom(iatom).GetSymmetricalAtom() == Monomers[indexA].GetAtom(k).GetSymmetricalAtom()){
		  Symi = k;
		}
	      }
	      else if(Symi == -1){
		//atomB is the atom index of iatom on monB
		int atomB = iatom - Na;
		if(Monomers[mon].GetAtom(k).GetSymmetricalAtom() == Monomers[indexB].GetAtom(atomB).GetSymmetricalAtom()){
		  Symi = k;
		}
	      }  
	    }//end of loop over k
	    if(Symi==-1){
	      printf("Error::Cluster::ComputeLatticeParamGrad():: Cannot find a symmetrical match for atom %i on m%i to an atom on m%i",
		     iatom,mon,indexA);
	      exit(0);
	    }
            */

	    //loops over xyz coordinates for the rotation matrix used under symmetry
	    for(int xyz = 0;xyz<3;xyz++){

	      //x vector
	      dE_dx1 += per_sca * damp * Rot(0,xyz) * (tmp_k_vec[0] + Fract[0]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] -  nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );
	      
	      dE_dx2 += per_sca * damp * Rot(0,xyz) * (tmp_k_vec[1] + Fract[1]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );

	      dE_dx3 += per_sca * damp * Rot(0,xyz) * (tmp_k_vec[2] + Fract[2]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );
	      
	      //y vector
	      dE_dy1 += per_sca * damp * Rot(1,xyz) * (tmp_k_vec[0] + Fract[0]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );

	      dE_dy2 += per_sca * damp * Rot(1,xyz) * (tmp_k_vec[1] +  Fract[1]) * 
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );

	      dE_dy3 += per_sca * damp * Rot(1,xyz) * (tmp_k_vec[2] + Fract[2]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );

	      //z vector 
	      dE_dz1 +=  per_sca * damp * Rot(2,xyz) * (tmp_k_vec[0] + Fract[0]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );      
	    
	      dE_dz2 +=  per_sca * damp* Rot(2,xyz) * (tmp_k_vec[1] + Fract[1]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );    

	      dE_dz3 +=  per_sca * damp * Rot(2,xyz) * (tmp_k_vec[2] + Fract[2]) *
		( Dimers[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * Dimers[i].GetMMGradient()[3*iatom+xyz] );

	      /*
	        if(xyz==0){
		  printf("d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		  printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		  printf("per_sca = %f\n",per_sca);
		  printf("DimerAtom = %i\n",DimerAtom);
	          printf("SymA = %i SymB = %i\n",SymA, SymB);
		  printf("atom = %i\n",iatom);
		  printf("3*atom+xyz = %i\n",3*iatom+xyz);
		  printf("damp = %f\n",damp);
		  printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
		  printf("Fract[0] = %f\n",Fract[0]);
		  printf("Dimer QM = %f\n",Dimers[i].GetQMGradient()[3*iatom+xyz]);
		  printf("dE_dx1 = %f\n",dE_dx1);
		  printf("\n");
		}
		*/	      

	      if (Params::Parameters().FindSpatialDampingGradient() ) {

		//x vectors
		dE_dx1 += per_sca * Rot(0,xyz) * (tmp_k_vec[0] + Fract[0]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );

		dE_dx2 += per_sca * Rot(0,xyz) * (tmp_k_vec[1] + Fract[1]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
		
		dE_dx3 += per_sca * Rot(0,xyz) * (tmp_k_vec[2] + Fract[2]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );	
		
		//y vector
		dE_dy1 += per_sca * Rot(1,xyz) * (tmp_k_vec[0] + Fract[0]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
	  
		dE_dy2 += per_sca * Rot(1,xyz) * (tmp_k_vec[1] +  Fract[1]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
		
		dE_dy3 += per_sca * Rot(1,xyz) * (tmp_k_vec[2] + Fract[2]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );

		//z vector
		dE_dz1 += per_sca * Rot(2,xyz) * (tmp_k_vec[0] + Fract[0]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
		
		dE_dz2 += per_sca * Rot(2,xyz) * (tmp_k_vec[1] + Fract[1]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
		
		dE_dz3 += per_sca * Rot(2,xyz) * (tmp_k_vec[2] + Fract[2]) 
		  * Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		  * (Dimers[i].GetQMIntEnergy() - nmb_scafac * Dimers[i].GetMMIntEnergy() );
                /*
	        if(xyz==0){
		  printf("d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		  printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		  printf("per_sca = %f\n",per_sca);
		  //printf("WhichMonB = %i\n",WhichMonB);
	          printf("SymA = %i SymB = %i\n",SymA, SymB);
		  printf("atom = %i\n",DimerAtom);
		  printf("3*DimerAtoms+xyz = %i\n",3*DimerAtom+xyz);
		  printf("Damp Grad = %f\n",Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*DimerAtom+xyz]);
		  printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
		  printf("Fract[0] = %f\n",Fract[0]);
		  printf("Int QM = %f\n",Dimers[i].GetQMIntEnergy());
		  printf("dE_dx1 = %f\n",dE_dx1);
		  printf("\n");
		}
		*/

	      }//end of if statement to determine if the spatial damping is used
	      
	    }//end over loop over xyz

	    //monomer energy already subtracted for counterpoise corrected 
	    if(!Params::Parameters().DoCounterpoise()){

              //x vector
	      dE_dx1 -= per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi];
	      dE_dx2 -= per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi];
	      dE_dx3 -= per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi];
	      
	      //y vector
	      dE_dy1 -= per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi+1];
	      dE_dy2 -= per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi+1];
	      dE_dy3 -= per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi+1];

	      //z vector 
	      dE_dz1 -=  per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi+2];  
              dE_dz2 -=  per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi+2];    
	      dE_dz3 -=  per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi+2];

	    }
	    /*
		  printf("d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		  printf("per_sca = %f\n",per_sca);
		  //printf("WhichMonB = %i\n",WhichMonB);
	          printf("SymA = %i SymB = %i\n",SymA, SymB);
		  printf("atom = %i\n",Symi);
		  printf("3*Symi = %i\n",3*Symi);
		  printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
		  printf("Fract[0] = %f\n",Fract[0]);
		  printf("damp = %f\n",damp);
		  printf("Grad Monomer = %f\n",Monomers[mon].GetQMGradient()[3*Symi]);
		  printf("dE_dx1 = %f\n",dE_dx1);
		  printf("\n");
	    */
	    

	    //For the monomer MM of the two-body interaction, a negative is being subtrated therefore it is positive.
	    //If the MM uses qchem and counterpoise corrected, the two the monomer contribution is already subtrated
	    if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3){

              //x vector
              dE_dx1 += per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi]; 
	      dE_dx2 += per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi];
	      dE_dx3 += per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi];
	      
	      //y vector
	      dE_dy1 += per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi+1];
	      dE_dy2 += per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi+1];
	      dE_dy3 += per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi+1];

	      //z vector 
	      dE_dz1 +=  per_sca * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi+2];  
	      dE_dz2 +=  per_sca * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi+2];   
	      dE_dz3 +=  per_sca * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi+2];

	    }

	  }//end of loop over iatom

	  // printf("after dE_dz2 = %f\n",dE_dz2);

	}//end of loop over j
      }//end of if statement for symmetry factor not being zero
    }//end of loop over i


    //adding contribution of the image dimers symmetrical to the image dimers
    for (int i=1; i<=NDim_images;i++){

      int indexA = DimerImages[i].GetIndexA();
      int indexB = DimerImages[i].GetReferenceMonomerIndex();

      int imageB = DimerImages[i].GetIndexB();

      int Na = Monomers[indexA].GetNumberOfAtoms();
      int Nb = Monomers[indexB].GetNumberOfAtoms();
      double damp = DimerImages[i].GetDampingFactor(c0,c1);
      Vector Fract;
      
      for(int j=0; j<DimerImages[i].GetSymmetryList().size()/2;j++){

	//Rotation Dimer for dimers on symmetry list
	Matrix Rot = DimerImages[i].GetRotationList()[j];
	Rot.Transpose();
	
	//the symmetrical dimers on SymmetryList
	int SymA = DimerImages[i].GetSymmetryList()[2*j];
	int SymB = DimerImages[i].GetSymmetryList()[2*j+1];
	
	//printf("d(%i,%i)\n",SymA,SymB);
	//printf("using d(%i,%i)\n",indexA,imageB);
	//printf("before dE_dz2 = %f\n",dE_dz2);

	//determining which monomer in dimer is outiside unit cell
	//if both monomers are in unit cell (for a non-image dimer), WhichMon is -1
	//bool WhichMonB = DimerImages[i].GetMonBList()[j];
	int k_vec[3];
	k_vec[0] = DimerImages[i].GetSymmetricalImageCell()[3*j];
	k_vec[1] = DimerImages[i].GetSymmetricalImageCell()[3*j+1];
	k_vec[2] = DimerImages[i].GetSymmetricalImageCell()[3*j+2];
	


	//determines if SymA is symmetrical to indexA or indexB by determining
	bool FirstMonA = false;
	for(int atom = 0; atom < Monomers[indexA].GetNumberOfAtoms();atom++){
	  if( (int)DimerImages[i].GetAtomEquivalency()[j][atom]== 1)
	    FirstMonA = true;
	}


	for(int iatom=0; iatom<Na+Nb; iatom++){


	  int mon;
	  //int Startgrad;
	  //int LastAtom;
	  int tmp_k_vec[3];
	  //determining which monomer in dimer is outiside unit cell
	  //bool WhichMonB =  DimerImages[i].GetMonBList()[j];

	  //when iatom is on monA
	  if(iatom<Na){
	    //mon = SymA;
	    //Startgrad = 0;
	    //LastAtom = Na;	     

	    //Fract = Monomers[mon].GetAtom(iatom).GetFractionalPosition();
	    /*
	    if(WhichMonB==0){
	      tmp_k_vec[0] = k_vec[0];
	      tmp_k_vec[1] = k_vec[1];
	      tmp_k_vec[2] = k_vec[2];
	    }else{
	      tmp_k_vec[0] = 0;
	      tmp_k_vec[1] = 0;
	      tmp_k_vec[2] = 0;
	    }
	    */
	  }
	  //when iatom is on monB
	  else{
	    //mon = SymB;
	    //Startgrad = 3*Na;
	    //LastAtom = Nb;

	    //Fract = Monomers[mon].GetAtom(iatom-Na).GetFractionalPosition();
	    /*
	    if(WhichMonB==1){
	      tmp_k_vec[0] = k_vec[0];
	      tmp_k_vec[1] = k_vec[1];
	      tmp_k_vec[2] = k_vec[2];
	    }else{
	      tmp_k_vec[0] = 0;
	      tmp_k_vec[1] = 0;
	      tmp_k_vec[2] = 0;
	    }
	    */
	  }

	  //that iatom is symmetrical equivalent to
	  int DimerAtom = (int)DimerImages[i].GetAtomEquivalency()[j][iatom] - 1;
	  
	  //The atom index on the monomer
	  int Symi;

	  //Determining atom index on monomer
	  if(FirstMonA){
	    if(DimerAtom < Na ){
	      Symi = DimerAtom;  
	      mon = SymA;
	      tmp_k_vec[0] = 0;
	      tmp_k_vec[1] = 0;
	      tmp_k_vec[2] = 0;
	    }else{
	      Symi = DimerAtom - Na;
	      mon = SymB;
	      tmp_k_vec[0] = k_vec[0];
	      tmp_k_vec[1] = k_vec[1];
	      tmp_k_vec[2] = k_vec[2];
	    }
	  }else{
	    if(DimerAtom < Nb ){
	      Symi = DimerAtom;  
	      mon = SymA;
	      tmp_k_vec[0] = 0;
	      tmp_k_vec[1] = 0;
	      tmp_k_vec[2] = 0;
	    }
	    else{
	      Symi = DimerAtom - Nb;
	      mon = SymB;
	      tmp_k_vec[0] = k_vec[0];
	      tmp_k_vec[1] = k_vec[1];
	      tmp_k_vec[2] = k_vec[2];
	    }
	  }

	  //Translation vector of the space group operator
	  Fract = Monomers[mon].GetAtom(Symi).GetFractTranslationVector(true);
          for(int j = 0;j<3;j++){
            Fract[j] += Monomers[mon].GetAtom(Symi).GetShiftVector()[j];
	  }	

	  //if monomer symmetry is employed, using gradient of the symmetrical unique monomer instead
	  if(UseMonSym)
	    mon = Monomers[mon].GetSymmetricalMonomer();
	  /*
	  int Symi = -1;

	  //find which atom in the dimers on the symmetry list is symmetrical to in the unique dimer
	  for(int k = 0; k<LastAtom;k++){
	    if(iatom<Na && Symi==-1){
	      if(Monomers[indexA].GetAtom(iatom).GetSymmetricalAtom() == Monomers[mon].GetAtom(k).GetSymmetricalAtom()){
		Symi = k;
	      }
	    }
	    else if(Symi==-1){
	      //atomB is the atom index of iatom on monB
	      int atomB = iatom - Na;
	      if(Monomers[indexB].GetAtom(atomB).GetSymmetricalAtom() == Monomers[mon].GetAtom(k).GetSymmetricalAtom()){
		Symi = k;
	      }
	      
	    }  
	  }//end of loop over k
	  if(Symi==-1){
	    printf("Error::Cluster::ComputeLatticeParamGrad:: Cannot find a symmetrical match for atom %i on m%i to an atom on m%i",
		   iatom,mon,indexB);
	    exit(0);
	  }
	  */

	  //loops over xyz coordinates for the rotation matrix used under symmetry
	  for(int xyz = 0;xyz<3;xyz++){

	    //x vector
	    dE_dx1 += 0.5 * damp * Rot(0,xyz) * (tmp_k_vec[0] + Fract[0]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );

	    dE_dx2 += 0.5 * damp * Rot(0,xyz) * (tmp_k_vec[1] + Fract[1]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz]);
	    
	    dE_dx3 += 0.5 * damp * Rot(0,xyz) * (tmp_k_vec[2] + Fract[2]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );	

	    //y vector
	    dE_dy1 += 0.5 * damp * Rot(1,xyz) * (tmp_k_vec[0] + Fract[0]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );
	    
	    dE_dy2 += 0.5 * damp * Rot(1,xyz) * (tmp_k_vec[1] + Fract[1]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );

	    dE_dy3 += 0.5 * damp * Rot(1,xyz) * (tmp_k_vec[2] + Fract[2]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );

	    //z vector 
	    dE_dz1 += 0.5 * damp * Rot(2,xyz) * (tmp_k_vec[0] + Fract[0]) *	      
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] ); 
	    
	    dE_dz2 +=  0.5 * damp * Rot(2,xyz) * (tmp_k_vec[1] + Fract[1]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] );          

	    dE_dz3 +=  0.5 * damp  * Rot(2,xyz) * (tmp_k_vec[2] + Fract[2]) * 
	      ( DimerImages[i].GetQMGradient()[3*iatom+xyz] - nmb_scafac * DimerImages[i].GetMMGradient()[3*iatom+xyz] ); 
	     
            /*
	    if(xyz==0){
	      printf("d(%i,%i)\n",DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB());
	      printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
	      //printf("WhichMonB = %i\n",WhichMonB);
	      printf("SymA = %i SymB = %i\n",SymA, SymB);
	      printf("atom = %i\n",DimerAtom);
	      printf("3*DimerAtom+xyz = %i\n",3*iatom+xyz);
	      printf("damp = %f\n",damp);
	      printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
	      printf("Fract[0] = %f\n",Fract[0]);
	      printf("Dimer QM = %f\n",DimerImages[i].GetQMGradient()[3*iatom+xyz]);
	      printf("dE_dx1 = %f\n",dE_dx1);
	      printf("\n");
	    }
            */
	    
	    if (Params::Parameters().FindSpatialDampingGradient() ) {

	      //x vector
	      dE_dx1 += 0.5 *  Rot(0,xyz) * (tmp_k_vec[0] + Fract[0]) 
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );	      
	      
	      dE_dx2 += 0.5 *  Rot(0,xyz) * (tmp_k_vec[1] + Fract[1])
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
	      
	      
	      dE_dx3 += 0.5 *  Rot(0,xyz) * (tmp_k_vec[2] + Fract[2])
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );	      
	      
	      //y vector
	      dE_dy1 += 0.5  * Rot(1,xyz) * (tmp_k_vec[0] + Fract[0]) 
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );	
	      
	      dE_dy2 += 0.5 *  Rot(1,xyz) * (tmp_k_vec[1] + Fract[1]) 
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
	      
	      dE_dy3 += 0.5 * Rot(1,xyz) * (tmp_k_vec[2] + Fract[2])  
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
	      
	      
	      //z vector
	      dE_dz1 += 0.5 * Rot(2,xyz)  * (tmp_k_vec[0] + Fract[0]) 
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
	      
	      dE_dz2 += 0.5 * Rot(2,xyz) * (tmp_k_vec[1] + Fract[1])
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
		
	      dE_dz3 += 0.5 *  Rot(2,xyz) * (tmp_k_vec[2] + Fract[2])
		* DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz]
		* ( DimerImages[i].GetQMIntEnergy() - nmb_scafac * DimerImages[i].GetMMIntEnergy() );
	      /*
                if(xyz==0 ){
		  printf("d(%i,%i)\n",DimerImages[i].GetIndexA(),Dimers[i].GetIndexB());
		  printf("Rot(0,%i) = %f\n",xyz,Rot(0,xyz));
		  //printf("WhichMonB = %i\n",WhichMonB);
	          printf("SymA = %i SymB = %i\n",SymA, SymB);
		  printf("atom = %i\n",DimerAtom);
		  printf("3*DimerAtoms+xyz = %i\n",3*DimerAtom+xyz);
		  printf("Damp Grad = %f\n",DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*DimerAtom+xyz]);
		  printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
		  printf("Fract[0] = %f\n",Fract[0]);
		  printf("int QM = %f\n", DimerImages[i].GetQMIntEnergy());
		  printf("dE_dx1 = %f\n",dE_dx1);
		  printf("\n");
		}
	      */
	      
	    }// if statement to determine if damping gradient is to be used.
	  }//end over loop over xyz


	    //monomer energy already subtracted for counterpoise corrected 
	    if(!Params::Parameters().DoCounterpoise()){

	      //x vector
	      dE_dx1 -= 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi] ;
	      dE_dx2 -= 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi] ;
	      dE_dx3 -= 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi] ;	
	      
	      //y vector
	      dE_dy1 -= 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi+1] ;
	      dE_dy2 -= 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi+1] ;
	      dE_dy3 -= 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi+1] ;
	          
	      //z vector 
	      dE_dz1 -= 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetQMGradient()[3*Symi+2] ;
	      dE_dz2 -= 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetQMGradient()[3*Symi+2] ;
	      dE_dz3 -= 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetQMGradient()[3*Symi+2] ;
	    }
	   
            /*
	      printf("d(%i,%i)\n",DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB());
	      //printf("WhichMonB = %i\n",WhichMonB);
	      printf("SymA = %i SymB = %i\n",SymA, SymB);
	      printf("atom = %i\n",Symi);
	      printf("3*Symi = %i\n",3*Symi);
	      printf("damp = %f\n",damp);
	      printf("tmp_k_vec[0] = %i\n",tmp_k_vec[0]);
	      printf("Fract[0] = %f\n",Fract[0]);
	      printf("Monomer QM = %f\n",Monomers[mon].GetQMGradient()[3*Symi]);
	      printf("dE_dx1 = %f\n",dE_dx1);
	      printf("\n");
            */

	    //For the monomer MM of the two-body interaction, a negative is subtrated therefore it is positive.
	    //If the MM uses qchem and counterpoise corrected, the two monomer contributions are already subtracted
	    if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3){
	      //x vector
	      dE_dx1 += 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi] ;
	      dE_dx2 += 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi] ;
	      dE_dx3 += 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi] ;	
	      
	      //y vector
	      dE_dy1 += 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi+1] ;
	      dE_dy2 += 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi+1] ;
	      dE_dy3 += 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi+1] ;	      
	      
	      //z vector 
	      dE_dz1 += 0.5 * damp * (tmp_k_vec[0] + Fract[0]) * Monomers[mon].GetMMGradient()[3*Symi+2] ;
	      dE_dz2 += 0.5 * damp * (tmp_k_vec[1] + Fract[1]) * Monomers[mon].GetMMGradient()[3*Symi+2] ;
	      dE_dz3 += 0.5 * damp * (tmp_k_vec[2] + Fract[2]) * Monomers[mon].GetMMGradient()[3*Symi+2] ;
	    }

	  
	}//end of loop over iatom

	//printf("after dE_dz2 = %f\n",dE_dz2);
	//printf("\n");

      }//end of loop over j
    }//end of loop over i

    cout << " dE_dx1 = " << dE_dx1 << "\n";
    cout << " dE_dy1 = " << dE_dy1 << "\n";
    cout << " dE_dz1 = " << dE_dz1 << "\n";
    cout << " dE_dx2 = " << dE_dx2 << "\n";
    cout << " dE_dy2 = " << dE_dy2 << "\n";
    cout << " dE_dz2 = " << dE_dz2 << "\n";
    cout << " dE_dx3 = " << dE_dx3 << "\n";
    cout << " dE_dy3 = " << dE_dy3 << "\n";
    cout << " dE_dz3 = " << dE_dz3 << "\n";

    /*
 //test by yoni
    Vector TestVector(9);
    TestVector[0] = dE_dx1;
    TestVector[1] = dE_dy1;
    TestVector[2] = dE_dz1;
    TestVector[3] = dE_dx2;
    TestVector[4] = dE_dy2;
    TestVector[5] = dE_dz2;
    TestVector[6] = dE_dx3;
    TestVector[7] = dE_dy3;
    TestVector[8] = dE_dz3;

    LatticeParamsGradFromVectorsGrad(TestVector).Print("Test Function");

    */

    //if (Params::Parameters().GetConstraintStr()=="constantpressure") {
    //at constant pressure
    if(Params::Parameters().FindEnthalpy()){
      Matrix ExtPressureMat(3,false);
      ExtPressureMat = FormGradientExternalPressureTerm();
      
      ExtPressureMat.PrintHessian("ExtPressureMat");
      dE_dx1 += ExtPressureMat.Element(0,0);
      dE_dx2 += ExtPressureMat.Element(1,0);
      dE_dx3 += ExtPressureMat.Element(2,0);
      dE_dy2 += ExtPressureMat.Element(1,1);
      dE_dy3 += ExtPressureMat.Element(2,1);
      dE_dz3 += ExtPressureMat.Element(2,2);
      
      cout << " dE_dx1 after gept = " << dE_dx1 << "\n";
      cout << " dE_dx2 after gept = " << dE_dx2 << "\n";
      cout << " dE_dx3 after gept = " << dE_dx3 << "\n";
      cout << " dE_dy1 after gept = " << dE_dy1 << "\n";
      cout << " dE_dy2 after gept = " << dE_dy2 << "\n";
      cout << " dE_dy3 after gept = " << dE_dy3 << "\n";
      cout << " dE_dz1 after gept = " << dE_dz1 << "\n";
      cout << " dE_dz2 after gept = " << dE_dz2 << "\n";
      cout << " dE_dz3 after gept = " << dE_dz3 << "\n";
    }


    //printf("Grun init = %i\n",Quasiharmonic::quasiharmonic().IsGruneisenInitialized());
    //Vibrational contribution to the gradient using QuasiHarmonic Approximation
    if((Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())){

      //Isotropic Gruneisen
      if(Params::Parameters().UseVolumeGruneisen()){
	Matrix VibGradMat;
	VibGradMat = ComputeQuasiHarmonicGradient();
	VibGradMat.Print("VibGradMat");
	dE_dx1 += VibGradMat.Element(0,0);
	dE_dx2 += VibGradMat.Element(1,0);
	dE_dx3 += VibGradMat.Element(2,0);
	dE_dy2 += VibGradMat.Element(1,1);
	dE_dy3 += VibGradMat.Element(2,1);
	dE_dz3 += VibGradMat.Element(2,2);
      }//Anisotropic Gruneisen
      else{
	Vector VibGrad;
	VibGrad = Quasiharmonic::quasiharmonic().ComputeAnisotropicGradient();

	printf("Helmholtz Gradient\n");
	printf("dF/da = %f\n",VibGrad[0]);
	printf("dF/db = %f\n",VibGrad[1]);
	printf("dF/dc = %f\n",VibGrad[2]);
	printf("dF/dalpha = %f\n",VibGrad[3]);
	printf("dF/dbeta = %f\n",VibGrad[4]);
	printf("dF/dgamma = %f\n\n",VibGrad[5]);

	//Adding Helmholtz Gradient to Full MM gradient since both are already expressed in terms of params 
	//rather than vectors 
	Grad_Lat_Full_MM[0] += VibGrad[0];
	Grad_Lat_Full_MM[1] += VibGrad[1];
	Grad_Lat_Full_MM[2] += VibGrad[2];
	Grad_Lat_Full_MM[3] += VibGrad[3];
	Grad_Lat_Full_MM[4] += VibGrad[4];
	Grad_Lat_Full_MM[5] += VibGrad[5];
      }
    }
    
    Grad[3*UniqueAtoms] = (dE_dx1 + Grad_Lat_Full_MM[0]);
    
    Grad[3*UniqueAtoms+1] = (dE_dx2*cos(gamma) + dE_dy2*sin(gamma) + Grad_Lat_Full_MM[1]);
    
    Grad[3*UniqueAtoms+2] = (cos(beta)*( dE_dx3 - dE_dy3 * cos(gamma)/sin(gamma) )+ dE_dy3 * cos(alpha)/sin(gamma)+ dE_dz3 * sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)/sin(gamma)/sin(gamma)-cos(alpha)*cos(alpha)/sin(gamma)/sin(gamma)-cos(beta)*cos(beta)/sin(gamma)/sin(gamma) )+ Grad_Lat_Full_MM[2]);
    
    Grad[3*UniqueAtoms+3] = (c /sin(gamma)*sin(alpha) * ( -dE_dy3 + dE_dz3 * ( cos(alpha)-cos(beta)*cos(gamma) ) / sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) ) ) + Grad_Lat_Full_MM[3]);
    
    Grad[3*UniqueAtoms+4] = (c*sin(beta) * (-dE_dx3 + dE_dy3 * cos(gamma)/sin(gamma) + dE_dz3 * ( cos(beta)-cos(alpha)*cos(gamma) ) / sin(gamma) / sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) ) ) + Grad_Lat_Full_MM[4]);
    
    Grad[3*UniqueAtoms+5] = (b * dE_dy2 * cos(gamma) + c * dE_dy3 /sin(gamma)* ( -cos(alpha)*cos(gamma)/sin(gamma) + cos(beta)/sin(gamma) ) - b * dE_dx2 *sin(gamma) + (c * dE_dz3 * ( cos(alpha)*cos(alpha)*cos(gamma) + cos(beta)*cos(beta)*cos(gamma) - 0.5*cos(alpha)*cos(beta)*( 3 + cos(2*gamma) ) )  /sin(gamma)/sin(gamma) ) /  sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) ) + Grad_Lat_Full_MM[5]);


    /*cout << "Final Lattice Gradient Terms \n";     
    for (int i=0;i<6;i++) {
      cout << "Grad["<< 3*UniqueAtoms  + i <<"] = " << Grad[3*UniqueAtoms+i] << "\n";
    }*/          

    //locking lattice params under symmetry
    if(b_locked_by_a){
      Grad[3*UniqueAtoms] += Grad[3*UniqueAtoms+1];
      Grad[3*UniqueAtoms+1] = 0.0;
    }
    if(c_locked_by_a){
      Grad[3*UniqueAtoms] += Grad[3*UniqueAtoms+2];
      Grad[3*UniqueAtoms+2] = 0.0;
    }
    else if(c_locked_by_b){
      Grad[3*UniqueAtoms+1] += Grad[3*UniqueAtoms+2];
      Grad[3*UniqueAtoms+2] = 0.0;
    }
    if(beta_locked_by_alpha){
      Grad[3*UniqueAtoms+3] += Grad[3*UniqueAtoms+4];
      Grad[3*UniqueAtoms+4] = 0.0;
    }
    if(gamma_locked_by_alpha){
      Grad[3*UniqueAtoms+3] += Grad[3*UniqueAtoms+5];
      Grad[3*UniqueAtoms+5] = 0.0;
    }
    else if(gamma_locked_by_beta){
      Grad[3*UniqueAtoms+4] += Grad[3*UniqueAtoms+5];
      Grad[3*UniqueAtoms+5] = 0.0;
    }

    //for lattice angles at are 90 or 120, their gradient is set to zero
    if(lock_alpha)
      Grad[3*UniqueAtoms+3] = 0.0;
    if(lock_beta)
      Grad[3*UniqueAtoms+4] = 0.0;
    if(lock_gamma)
      Grad[3*UniqueAtoms+5] = 0.0;

    //PrintGradient("Grad after locking Lattice Params under symmetry",Grad);

    cout << "Final Lattice Gradient Terms \n";     
    for (int i=0;i<6;i++) {
      cout << "Grad["<< 3*UniqueAtoms + i <<"] = " << Grad[3*UniqueAtoms+i] << "\n";  
    } 

      /*
      for (int j=1;j<=NDim_images;j++) {
	tmp_k_vec = DimerImages[j].GetImageCell();
	double symfac = 1.0; //
	if ( Params::Parameters().UseCrystalSymmetry() )
	  symfac *= DimerImages[j].GetSymmetryFactor();
	
	
	int indexA = DimerImages[j].GetIndexA();
	int indexB = DimerImages[j].GetIndexB();
	int indexB_ref = DimerImages[j].GetReferenceMonomerIndex();
	int Na = DimerImages[j].GetMonomerA().GetNumberOfAtoms();
	int Nb = DimerImages[j].GetMonomerB().GetNumberOfAtoms();
	Monomer mB = Monomers[indexB_ref];
	Monomer mA = Monomers[indexA];
	double damp = DimerImages[j].GetDampingFactor(c0,c1);
	double sep = DimerImages[j].GetDimerSeparation();
	
	
	for (int k=Na;k<Na+Nb;k++) {
	  
	  if (i==1) {
	    dE_dx1 +=  0.5 * damp * symfac * tmp_k_vec[i-1]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-1]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    
	    dE_dx2 +=  0.5 * damp * symfac * tmp_k_vec[i]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    
	    dE_dx3 +=  0.5 * damp * symfac * tmp_k_vec[i+1]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i+1]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	  }           
	  
	  if (i==2) { 
	    dE_dy1 += 0.5 * damp * symfac * tmp_k_vec[i-2]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-2]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    //0.0;

	    dE_dy2 +=  0.5 * damp * symfac * tmp_k_vec[i-1]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-1]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    
	    dE_dy3 +=  0.5 * damp * symfac * tmp_k_vec[i]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	  }
	  
	  if (i==3) {       
	    dE_dz1 +=  0.5 * damp * symfac * tmp_k_vec[i-3]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-3]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    //0.0;              
	    
	    dE_dz2 +=  0.5 * damp * symfac * tmp_k_vec[i-2]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-2]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	    //0.0;     
	    
	    dE_dz3 +=  0.5 * damp * symfac * tmp_k_vec[i-1]  * ( DimerImages[j].GetQMGradient()[3*k+i-1]-mB.GetQMGradient()[3*(k-Na)+i-1] )
	      - 0.5 * damp * symfac * tmp_k_vec[i-1]  * nmb_scafac * ( DimerImages[j].GetMMGradient()[3*k+i-1]-mB.GetMMGradient()[3*(k-Na)+i-1] );
	  }
	  
	  
	  if (Params::Parameters().FindSpatialDampingGradient() ) {
	    if (i==1) {
	      dE_dx1 += 0.5* symfac * tmp_k_vec[i-1] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      
	      dE_dx2 += 0.5* symfac * tmp_k_vec[i] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      
	      dE_dx3 += 0.5* symfac * tmp_k_vec[i+1] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	    }
	    if (i==2) {
	      dE_dy1 += 0.5* symfac * tmp_k_vec[i-2] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      //0.0;
	      
	      dE_dy2 += 0.5* symfac * tmp_k_vec[i-1] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      
	      dE_dy3 += 0.5* symfac * tmp_k_vec[i] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	    }
	    if (i==3) {
	      dE_dz1 += 0.5* symfac * tmp_k_vec[i-3] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      //0.0;
	      
	      dE_dz2 += 0.5* symfac * tmp_k_vec[i-2] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	      //0.0;
	      
	      dE_dz3 += 0.5* symfac * tmp_k_vec[i-1] *DimerImages[j].GetSpatialDampingFunctionGradient(c0,c1)[3*k+i-1]
		* (DimerImages[j].GetQMIntEnergy()- nmb_scafac * DimerImages[j].GetMMIntEnergy() );
	    }
	  }
	}
      }
    }
    
    cout << " dE_dx1 = " << dE_dx1 << "\n";
    cout << " dE_dy1 = " << dE_dy1 << "\n";
    cout << " dE_dz1 = " << dE_dz1 << "\n";
    cout << " dE_dx2 = " << dE_dx2 << "\n";
    cout << " dE_dy2 = " << dE_dy2 << "\n";
    cout << " dE_dz2 = " << dE_dz2 << "\n";
    cout << " dE_dx3 = " << dE_dx3 << "\n";
    cout << " dE_dy3 = " << dE_dy3 << "\n";
    cout << " dE_dz3 = " << dE_dz3 << "\n";
    
    double dU11_da = 1.0;
    double dU21_db = cos(gamma);
    double dU21_dgamma = -b*sin(gamma);
    double dU22_db = sin(gamma);
    double dU22_dgamma = b*cos(gamma);
    double dU31_dc = cos(beta);
    double dU31_dbeta = -c*sin(beta);
    double dU32_dc = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
    double dU32_dalpha = -c*sin(alpha)/sin(gamma);
    double dU32_dbeta = c*sin(beta)*cos(gamma)/sin(gamma);
    double dU32_dgamma = c*(cos(beta) - beta_term*cos(gamma)/sin(gamma));
    double dU33_dc = gamma_term;
    double dU33_dalpha = - beta_term * dU32_dalpha / gamma_term;
    double dU33_dbeta =  ( cos(beta)*sin(beta) - beta_term*dU32_dbeta ) / gamma_term;
    double dU33_dgamma = -  beta_term*dU32_dgamma / gamma_term;
    
    cout << " dU11_da = " << dU11_da << "\n";
    cout << " dU21_db = "	<< dU21_db << "\n";
    cout << " dU21_dgamma = " << dU21_dgamma << "\n";
    cout << " dU22_db = "	<< dU22_db << "\n";
    cout << " dU22_dgamma = " << dU22_dgamma << "\n";
    cout << " dU31_dc = "	<< dU31_dc << "\n";
    cout << " dU31_dbeta = "	<< dU31_dbeta << "\n";
    cout << " dU32_dc = "	<< dU32_dc << "\n";
    cout << " dU32_dalpha = "	<< dU32_dalpha << "\n";
    cout << " dU32_dbeta = "	<< dU32_dbeta << "\n";
    cout << " dU32_dgamma = "	<< dU32_dgamma << "\n";
    cout << " dU33_dc = "	<< dU33_dc << "\n";
    cout << " dU33_dalpha = "	<< dU33_dalpha << "\n";
    cout << " dU33_dbeta = "	<< dU33_dbeta << "\n";
    cout << " dU33_dgamma = "	<< dU33_dgamma << "\n";
   
    
    if ( Params::Parameters().GetLatticeHessian() || Params::Parameters().GetElasticConstants() ) {
      Store_dE_dv_q(dE_dx1, dE_dy1, dE_dz1, dE_dx2, dE_dy2, dE_dz2, dE_dx3, dE_dy3, dE_dz3);
    }
    
    
      */
    
    
    /*

    if (Params::Parameters().ComputeStressTensor()==1 ) {
      
      Matrix Full_MM_Stress_tensor(3,3);// = Compute_dE_MM_PBC_stress_tensor();
      FormStressTensor(dE_dx1, dE_dy1, dE_dz1,
		       dE_dx2, dE_dy2, dE_dz2,
		       dE_dx3, dE_dy3, dE_dz3,
		       Full_MM_Stress_tensor);
      
      
	//FormStressTensor(dE_dx1, dE_dy1, dE_dz1,     
	//dE_dx2, dE_dy2, dE_dz2,
	//dE_dx3, dE_dy3, dE_dz3,
	//Grad_Lat_Full_MM[0],Grad_Lat_Full_MM[1],Grad_Lat_Full_MM[2],
	//Grad_Lat_Full_MM[3],Grad_Lat_Full_MM[4],Grad_Lat_Full_MM[5],
	//dU11_da, dU21_db, dU21_dgamma,dU22_db, dU22_dgamma,
	//dU31_dc, dU31_dbeta,
	//dU32_dc, dU32_dalpha, dU32_dbeta, dU32_dgamma,
	//dU33_dc, dU33_dalpha, dU33_dbeta, dU33_dgamma,
	//a, b, c, alpha, beta, gamma); 
      
      //stress_tensor.PrintHessian("stress_tensor");

      */  
    }    


  //return Grad;
}  

double Cluster::New_ComputeLatticeParamGradFullMM(int i, double delta) {

  if(!Params::Parameters().GetMMType()==1){
    printf("Cluster::New_ComputeLatticeParamGradFullMM. Has not been adapted for any other MM type other than Tinker. ");
    exit(1);
  }

  //return zero for params locked by symmetry 
  if(b_locked_by_a && i == 1) 
    return 0.0;
  else if((c_locked_by_a ||c_locked_by_b) && i == 2)
    return 0.0;
  else if(lock_alpha && i == 3)
    return 0.0;
  else if((beta_locked_by_alpha || lock_beta) && i == 4)
    return 0.0;
  else if((gamma_locked_by_alpha || gamma_locked_by_beta || lock_gamma) && i == 5)
    return 0.0;

  string cmd;
  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  cmd = "cd " + job_path;
  cmd += "; ";      

  double full_lat_energy = 0.0;

//if TINKER      
  if (Params::Parameters().GetMMType()==1)  { //tinker

    //First Command: Make a copy of full.key
    cmd += "cp full.key full_copy.key; ";
    cmd += "cp full.xyz full_copy.xyz; ";

    //Second Command: Delete the existing full.key
//    cmd += "rm -f full.key; ";

    // Run these 2 commands
    system(cmd.c_str());

    //Rewrite the full.key with one of the lattice params modified
    string keyfile = job_path + "/full.key";
    //string keyfile = Params::Parameters().GetMMPath() + "/full.key";

    /* Overwrite the keyfile */
    // Open the file for writing, write the Tinker rem section to it,
    // and close the file.
    FILE *key;

    if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
      printf("Cluster::New_ComputeLatticeParamGradFullMM : Cannot open file '%s'\n",
              keyfile.c_str());
      exit(1);
    }

    fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
    double lat[6];
    lat[0] = unit_cell[0].Norm();
    lat[1] = unit_cell[1].Norm();
    lat[2] = unit_cell[2].Norm();

    lat[3] = RadiansToDegrees*
        acos(unit_cell[1].DotProduct( unit_cell[2] ) / (lat[1]*lat[2]));
    lat[4]  = RadiansToDegrees*
        acos(unit_cell[0].DotProduct( unit_cell[2] ) / (lat[0]*lat[2]));
    lat[5] = RadiansToDegrees*
        acos(unit_cell[0].DotProduct( unit_cell[1] ) / (lat[0]*lat[1]));

    if (i >= 6) 
      cout << "Cluster::New_ComputeLatticeParamGradFullMM(i,delta) Error: i cannot exceed 6 \n";
    else {
      if (i<3) {
        lat[i] = lat[i] + delta;
      }
      else {
        lat[i] = lat[i] + delta*10;
      }
    }

  //preserving symmetry 
  if(b_locked_by_a) 
    lat[1] = lat[0];
  if(c_locked_by_a) 
    lat[2] = lat[0];
  if(c_locked_by_b)
    lat[2] = lat[1];
  if(beta_locked_by_alpha)
    lat[4] = lat[3];
  if(gamma_locked_by_alpha)
    lat[5] = lat[3];
  if(gamma_locked_by_beta)
    lat[5] = lat[4];

    /*
    printf("Tinker lattice parameters:\n");
    printf("a = %f, b = %f, c = %f\n",lat[0],lat[1],lat[2]);
    printf("alpha = %f, beta = %f, gamma = %f\n",lat[3],lat[4],lat[5]);
    */
    fprintf(key,"# Periodic boundary conditions\n");
    fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",lat[0],lat[1],lat[2]);
    fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",lat[3],lat[4],lat[5]);
    fprintf(key,"EWALD\t\tTRUE\n");

    // If we want vacuum boundary conditions for the ewald sum, set
    // EWALD-BOUNDARY flag.  For tinfoil (infinite dielectric) boundary
    // conditions, we omit this keyword.  For non-polar unit cells, it
    // doesn't matter.  But for polar ones, tinfoil boundary conditions
    // seem to work better.
    if (!Params::Parameters().TinFoilBoundaryConditions())
      fprintf(key,"EWALD-BOUNDARY\tTRUE\n");     
//    fprintf(key,"EWALD-CUTOFF\t\t10.0\n");
   
    fclose(key);

    //changes coordinates of the full system to preserve symmetry
    string xyzfile = job_path  + "/full.xyz";
    //string xyzfile = Params::Parameters().GetMMPath() + "/full.xyz";
    /* Overwrite the filename */
    // Open the file for writing, write the coodinates
    // and close the file.
    FILE *xyz;
    if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
      printf("Cluster::New_ComputeLatticeParamGradFullMM : Cannot open file '%s'\n",
              xyzfile.c_str());
      exit(1);
    }

  //Get the coordinates ofthe symmetry unique monomers and use that determine the coordinates for the rest of the system
  Vector Coords = GetSymmetryUniqueCoordinates(true);
  for(int j=0;j<6;j++)
    Coords[3*UniqueAtoms+j] = lat[j];
  
  MaintainCartesianSymmetry(Coords,true);
  //Coords.PrintGradient("Coords");
  //exit(0);
  //Coords = CartFromFractAndLatticeParams(Coords,true,lat[0],lat[1],lat[2],lat[3],lat[4],lat[5]);
  Coords = CartFromFractAndLatticeParams(Coords,true);
  Coords = GetSymmetryImposedCoordinates(Coords,true,true);
  //Coords.PrintGradient("Coords"); 
  //Coords = CartFromFractAndLatticeParams(Coords,false,lat[0],lat[1],lat[2],lat[3],lat[4],lat[5]);
  Coords = CartFromFractAndLatticeParams(Coords,false); 
  //Coords.PrintGradient("Coords"); 

  //Coords = CartFromFractAndLatticeParams(Coords,false);
  Vector TempCoord(Coords.GetLength()-6);
  for(int j=0;j<Coords.GetLength()-6;j++){
    TempCoord[j] = Coords[j];
  }
  //TempCoord.PrintGradient("TempCoord");
  Coords = TempCoord;
  //Coords.PrintGradient("Coords");
  //exit(0);

  //printf("i = %i delta = %f\n",i,delta); 
  //Coords.PrintGradient("Coords");

  //Creating new tinker xyz file
  PrintTinkerCartesian(Coords,xyz);
      
  fclose(xyz);
  //exit(0);

  string cmd2;
  // First command, change to local directory
  cmd2 = "cd " + job_path;
  cmd2 += "; ";

  //First command: run full.xyz with new full.key
  cmd2 += "analyze full.xyz e > full_lat.out;";  
  
  //replace full.xyz and full.key  with the original
  cmd2 += "cp full_copy.xyz full.xyz;";
  cmd2 += "cp full_copy.key full.key;";
  
  // Switch back to base directory
  cmd2 += "cd " + Params::Parameters().GetBasePath();
  
  //run these  commands
  system(cmd2.c_str());    
  
  string out_filename = job_path + "/full_lat.out";  

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::New_ComputeLatticeParamGradFullMM(i,delta) Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;   
  while ( !infile.eof() ) {
    getline(infile,line);  
    string match = line.substr(0,23);       
    if ( match==" Total Potential Energy" ) {
      istringstream iss(line);
      string tmp;          
      // throw away "text tags"       
      iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp;
      iss >> full_lat_energy; // read the energy
    }


  }

  full_lat_energy /= HartreesToKcalpermole; // convert to hartrees
  //cout << "full_lat_energy = " << full_lat_energy <<"\n";
  // Close the force file
  infile.close();
 

  }
  return  full_lat_energy;
}

//allows for convertion between Cartesian and Fractional coordinates with custom latticeparameter 
Vector Cluster::CartFromFractAndLatticeParams(Vector StartCoords,bool reverse){
//Vector Cluster::CartFromFractAndLatticeParams(Vector StartCoords,bool reverse,double a,double  b,double c,
//					      double alpha,double beta,double gamma){

  //alpha*= DegreesToRadians;
  //beta*= DegreesToRadians;
  //gamma*= DegreesToRadians;

  int Natom = StartCoords.GetLength()/3 -2;
  double a = StartCoords[3*Natom];
  double b = StartCoords[3*Natom+1];
  double c = StartCoords[3*Natom+2];
  double alpha = StartCoords[3*Natom+3]*DegreesToRadians;
  double beta = StartCoords[3*Natom+4]*DegreesToRadians;
  double gamma = StartCoords[3*Natom+5]*DegreesToRadians;					      

  double  beta_term = ( cos(alpha) - cos(beta)*cos(gamma) ) / sin(gamma) ;
  double  gamma_term = sqrt( 1 - cos(beta)*cos(beta) - beta_term*beta_term ) ;

  Matrix Convertion(3,3);
  Convertion(0,0) = a;
  Convertion(0,1) = b*cos(gamma);
  Convertion(0,2) = c*cos(beta);
  Convertion(1,0) = 0.0;
  Convertion(1,1) = b*sin(gamma);
  Convertion(1,2) = c*beta_term;
  Convertion(2,0) = 0.0;
  Convertion(2,1) = 0.0;
  Convertion(2,2) = c*gamma_term;
  
  if(reverse)
    Convertion.Inverse();
  //LatticeVectors.Print("LatticeVectors");

  Vector FinalCoords(StartCoords.GetLength());
  for(int i=0;i<StartCoords.GetLength()/3-2;i++){
    for(int j=0;j<3;j++){
      FinalCoords[3*i] += Convertion(0,j)*StartCoords[3*i+j];
      FinalCoords[3*i+1] += Convertion(1,j)*StartCoords[3*i+j];
      FinalCoords[3*i+2] += Convertion(2,j)*StartCoords[3*i+j];
    }
  }

  for(int j=0;j<6;j++)
   FinalCoords[3*Natom+j] = StartCoords[3*Natom+j];

  //FinalCoords.PrintGradient("FinalCoords");
  //exit(0);
  return FinalCoords;
}

void Cluster::ComputeLatticeParamGradFullMM() {
  string cmd;
  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath();  
  cmd = "cd " + job_path;
  cmd += "; "; 

//if TINKER 
  if (Params::Parameters().GetMMType()==1)  { //tinker

    // First command, run the job

    cmd += "xtalmin full.xyz 1000000 > full_xtalmin.out; ";

    // Second and third commands, remove force.dat and full.xyz_2
    cmd += "rm -f force.dat; ";
    cmd += "rm -f full.xyz_*; ";

    // Fourth command, switch back to base directory
    cmd += "cd " + Params::Parameters().GetBasePath();
 
    //run these 3 commands
    system(cmd.c_str());
  }
// Same set for AIFF not yet included as of July 2010
//  ReadLatticeParamGradFullMM( grad_lat_full_MM);
}

Vector Cluster::ReadLatticeParamGradFullMM(Vector grad_lat_full_MM) {

  //Set up the xtalmin filename, with the full path.  File is full_xtalmin.out
  string job_path = Params::Parameters().GetMMPath();

  string out_filename = job_path + "/full_xtalmin.out";
  // Open the xtalmin file
  ifstream infile;   
  infile.open( out_filename.c_str() );       
  if ( !infile.is_open() ) {
    printf("Cluster::ComputeLatticeParamGradFullMM : Cannot open file '%s'\n",
           out_filename.c_str());
    exit(1);
  }    

  //Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,36);

    // Search for the six Gradient components depending on lattice params
    if ( match==" Gradient due to the change in Abox " ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[0]; // read energy
    }
    if ( match==" Gradient due to the change in Bbox " ) {            
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[1]; // read energy
    }
    if ( match==" Gradient due to the change in Cbox " ) {            
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[2]; // read energy
    }
    if ( match==" Gradient due to the change in Alpha" ) {            
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[3]; // read energy
    }
    if ( match==" Gradient due to the change in Beta " ) {            
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[4]; // read energy
    }
    if ( match==" Gradient due to the change in Gamma" ) {            
      istringstream iss(line);
      string tmp;
      for (int i=0;i<7;i++)
        iss >> tmp; // throw away text
      iss >> grad_lat_full_MM[5]; // read energy
    }
  }
  // Close the file_xtalmin.out file
  infile.close();
  grad_lat_full_MM.Scale(0.529177249/627.5095);
  return grad_lat_full_MM;
}

Vector Cluster::GetCurrentCoordinates(bool GetLatticeParams) {

  if ( (((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly())
       && GetLatticeParams) {
    Vector AtomicCoordinatesPlusLattice(AtomicCoordinates.GetLength()+6);
    for (int i=0;i<AtomicCoordinates.GetLength();i++) {
      AtomicCoordinatesPlusLattice[i]=AtomicCoordinates[i];
    }
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()]=UnitCellAxes[0];
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+1]=UnitCellAxes[1];
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+2]=UnitCellAxes[2];
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+3]=UnitCellAngles[0];
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+4]=UnitCellAngles[1];
    AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+5]=UnitCellAngles[2];
/*
    if (Params::Parameters().FreezeUnitCellAngles() ) {
      AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+3] *= 0.0 ;
      AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+4] *= 0.0 ;
      AtomicCoordinatesPlusLattice[AtomicCoordinates.GetLength()+5] *= 0.0 ;
    }
*/
    return AtomicCoordinatesPlusLattice;
  }
  else
    return AtomicCoordinates;
}

Vector Cluster::GetLastAcceptedCoordinates(bool GetLatticeParams) {

  if ( (((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly())
       && GetLatticeParams) {
    int tmp = AcceptedAtomicCoordinates.GetLength();
    Vector AcceptedAtomicCoordinatesPlusLattice(tmp+6);

    for (int i=0;i<AcceptedAtomicCoordinates.GetLength();i++) {
      AcceptedAtomicCoordinatesPlusLattice[i]=AcceptedAtomicCoordinates[i];
    }

    AcceptedAtomicCoordinatesPlusLattice[tmp]=AcceptedUnitCellAxes[0];
    AcceptedAtomicCoordinatesPlusLattice[tmp+1]=AcceptedUnitCellAxes[1];
    AcceptedAtomicCoordinatesPlusLattice[tmp+2]=AcceptedUnitCellAxes[2];
    AcceptedAtomicCoordinatesPlusLattice[tmp+3]=AcceptedUnitCellAngles[0];
    AcceptedAtomicCoordinatesPlusLattice[tmp+4]=AcceptedUnitCellAngles[1];
    AcceptedAtomicCoordinatesPlusLattice[tmp+5]=AcceptedUnitCellAngles[2];

    return AcceptedAtomicCoordinatesPlusLattice;
  }
  else
    return AcceptedAtomicCoordinates;
}

double Cluster::UpdateLatticeParamsAndComputeEnergy(int i, Vector coords, double delta) {
  double E_Lat;
  coords.Print("coords inside UpdateLatticeParamsAndComputeEnergy");
  fflush(stdout);

/*
  Vector New_UnitCellParams(6);
  UnitCellAxes.Print(" UnitCellAxes at start");
  UnitCellAngles.Print(" UnitCellAngles at start");
  fflush(stdout);
  for (int j=0; j<3;j++) {
    New_UnitCellParams[j]= UnitCellAxes[j];
  }
  for (int j=0; j<3;j++) {
    New_UnitCellParams[j+3]= UnitCellAngles[j];
  }
  New_UnitCellParams.Print("New_UnitCellParams before initialization");
  New_UnitCellParams[i] += delta;
  New_UnitCellParams.Print("New_UnitCellParams after initialization");
  fflush(stdout);
*/
  // Storing the old unit cell vector list    
  Vector* new_unit_cell;
  new_unit_cell = new Vector[3];                
  for (int j=0;j<3;j++) {
    new_unit_cell[j] = unit_cell[j];             
    unit_cell[j].Print("unit_cell");
    new_unit_cell[j].Print("new_unit_cell");
    fflush(stdout);
  }
/* 
  //Re-calculating the unit cell vector with changed axes/angles
  double a = New_UnitCellParams[0];
  double b = New_UnitCellParams[1];
  double c = New_UnitCellParams[2];
  double alpha = New_UnitCellParams[3];
  double beta = New_UnitCellParams[4];
  double gamma = New_UnitCellParams[5];
*/


  int Natoms = (coords.GetLength() - 6) / 3;
  double a = coords[3*Natoms];
  double b = coords[3*Natoms+1];
  double c = coords[3*Natoms+2];
  double alpha = coords[3*Natoms+3];
  double beta = coords[3*Natoms+4];
  double gamma = coords[3*Natoms+5];

  double alpha_rad = alpha*DegreesToRadians;
  double beta_rad = beta*DegreesToRadians;
  double gamma_rad = gamma*DegreesToRadians;

  double beta_term =             
    (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
  double gamma_term =      
    sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);

  /*
  cout << "a = " << a << "\n"; 
  cout << "b = " << b << "\n";
  cout << "c = " << c << "\n";
  cout << "alpha = " << alpha << "\n";
  cout << "beta = " << beta << "\n";
  cout << "gamma = " << gamma << "\n";
  fflush(stdout);
  */

  // v1
  unit_cell[0][0] = a;
  unit_cell[0][1] = 0;       
  unit_cell[0][2] = 0;

  // v2
  unit_cell[1][0] = b*cos(gamma_rad);
  unit_cell[1][1] = b*sin(gamma_rad);
  unit_cell[1][2] = 0;

  // v3          
  unit_cell[2][0] = c*cos(beta_rad);
  unit_cell[2][1] = c*beta_term;
  unit_cell[2][2] = c*gamma_term;

  SetNewCoordinates(coords);
  RunJobsAndComputeEnergy();
  fflush(stdout);
  E_Lat = GetHMBIEnergy();
  cout << "E_Lat["<<i<<"] =" << setprecision(11) << showpoint << E_Lat <<"\n";
  //Revert back to original unit_cell matrix
  for (int j=0;j<3;j++) {
    unit_cell[j].Print("unit_cell before reverting");   
    unit_cell[j] = new_unit_cell[j];  
    unit_cell[j].Print("unit_cell after reverting");
    fflush(stdout);
  }
//  New_UnitCellParams[i] -= delta;
//  New_UnitCellParams.Print("New_UnitCellParams after reverting");
  fflush(stdout);
  return E_Lat;
}

void Cluster::UpdateLatticeParams(Vector coords) {

  //Lattice Parameters are stored in the last 6 entries of coords
  int Natoms = (coords.GetLength()-6)/3;

  //coords.Print("coord in UpdateLatticeParams\n");

  //int Natoms = GetTotalNumberOfAtoms();
  //int UniqueAtoms = GetNumberOfUniqueAtoms();
 
  /*
   cout << "before \n" <<"\n";
    cout << "a = " << UnitCellAxes[0] << "\n"; 
    cout << "b = " << UnitCellAxes[1] << "\n";
    cout << "c = " << UnitCellAxes[2] << "\n";
    cout << "alpha = " << UnitCellAngles[0] << "\n";
    cout << "beta = " << UnitCellAngles[1] << "\n";
    cout << "gamma = " << UnitCellAngles[2] << "\n";
    cout << "\n";
    fflush(stdout);
  */

  /*
  //saving the old lattice vectors
  Matrix old_unit_cell(3,3);
  // v1
  old_unit_cell(0,0) = unit_cell[0][0];
  old_unit_cell(0,1) = unit_cell[0][1];
  old_unit_cell(0,2) = unit_cell[0][2];
  
  // v2
  old_unit_cell(1,0) = unit_cell[1][0];
  old_unit_cell(1,1) = unit_cell[1][1];             
  old_unit_cell(1,2) = unit_cell[1][2];

  // v3
  old_unit_cell(2,0) = unit_cell[2][0];
  old_unit_cell(2,1) = unit_cell[2][1];
  old_unit_cell(2,2) = unit_cell[2][2];


  //saving old lattice parameters
  Vector OldUnitCellAxes = UnitCellAxes;  
  Vector OldUnitCellAngles = UnitCellAngles;

  */



  if (!Params::Parameters().FreezeUnitCellParams() && Params::Parameters().UseLatticeSymmetry()) { //JLM
    /*printf("Using lattice symmetry in update lattice params!\n\n");
    printf("b_locked_by_a = %i\n",b_locked_by_a);
    printf("c_locked_by_a = %i\n",c_locked_by_a);
    printf("c_locked_by_b = %i\n",c_locked_by_b);
    printf("lock_alpha = %i\n",lock_alpha);
    printf("lock_beta = %i\n",lock_beta);
    printf("lock_gamma = %i\n",lock_gamma);
    printf("beta_locked_by_alpha = %i\n",beta_locked_by_alpha);
    printf("gamma_locked_by_alpha = %i\n",gamma_locked_by_alpha);
    printf("gamma_locked_by_beta = %i\n",gamma_locked_by_beta);*/

    for (int i=0; i<3;i++) {    
      UnitCellAxes[i] = coords[3*Natoms+i];
      //UnitCellAngles[i] = coords[3*UniqueAtoms+3+i];
    }

    //locking lattice axis by symmetry
    if(b_locked_by_a){
      UnitCellAxes[1] = UnitCellAxes[0];
    }
    if(c_locked_by_a){
      UnitCellAxes[2] = UnitCellAxes[0];
    }
    else if(c_locked_by_b){
      UnitCellAxes[2] = UnitCellAxes[1];
    }


    //if lattice parameters are not 90 or 120, alter them
    if(!lock_alpha)
      UnitCellAngles[0] = coords[3*Natoms+3];
    if(!lock_beta) 
      UnitCellAngles[1] = coords[3*Natoms+4];
    if(!lock_gamma)
      UnitCellAngles[2] = coords[3*Natoms+5];
    
    //locking angles by symmetry
    if(beta_locked_by_alpha)
      UnitCellAngles[1] = UnitCellAngles[0];
    if(gamma_locked_by_alpha)
      UnitCellAngles[2] = UnitCellAngles[0];
    else if(gamma_locked_by_beta)
      UnitCellAngles[2] = UnitCellAngles[1];
     
  }
  else {
    for (int i=0; i<3;i++) {
      //if (FreezeLatParams[i] == 0) 
        UnitCellAxes[i] = coords[3*Natoms+i];
    }
    for (int i=0; i<3;i++) {  
      //if (FreezeLatParams[i+3] == 0)
        UnitCellAngles[i] = coords[3*Natoms+3+i];
    }
  }

  //cout << "Params::Parameters().GetNumericalStiffnessTensor_Effective() =" << Params::Parameters().GetNumericalStiffnessTensor_Effective() << "\n";
  //if (!Params::Parameters().GetNumericalStressTensor_Effective() && !Params::Parameters().GetNumericalStiffnessTensor_Effective() ) {
    //Re-calculating the unit cell vector with changed axes/angles
    double a = UnitCellAxes[0];
    double b = UnitCellAxes[1];
    double c = UnitCellAxes[2];
    double alpha = UnitCellAngles[0];
    double beta = UnitCellAngles[1];
    double gamma = UnitCellAngles[2];
    double alpha_rad = alpha*DegreesToRadians;
    double beta_rad = beta*DegreesToRadians;
    double gamma_rad = gamma*DegreesToRadians;
    double beta_term =
      (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
    double gamma_term =
      sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);
    
    
    /*cout << "after \n" <<"\n";
    cout << "a = " << a << "\n"; 
    cout << "b = " << b << "\n";
    cout << "c = " << c << "\n";
    cout << "alpha = " << alpha << "\n";
    cout << "beta = " << beta << "\n";
    cout << "gamma = " << gamma << "\n";
    fflush(stdout);*/
     
    


    // v1
    unit_cell[0][0] = a;
    unit_cell[0][1] = 0;
    unit_cell[0][2] = 0;

    // v2
    unit_cell[1][0] = b*cos(gamma_rad);
    unit_cell[1][1] = b*sin(gamma_rad);             
    unit_cell[1][2] = 0;

    // v3
    unit_cell[2][0] = c*cos(beta_rad);
    unit_cell[2][1] = c*beta_term;    
    unit_cell[2][2] = c*gamma_term;
    //}

    cell_volume = unit_cell[0][0]*unit_cell[1][1]*unit_cell[2][2];

    cout << "a = " << UnitCellAxes[0] << "\n"; 
    cout << "b = " << UnitCellAxes[1] << "\n";
    cout << "c = " << UnitCellAxes[2] << "\n";
    cout << "alpha = " << UnitCellAngles[0] << "\n";
    cout << "beta = " << UnitCellAngles[1] << "\n";
    cout << "gamma = " << UnitCellAngles[2] << "\n";
    cout << "cell_volume = " << cell_volume << " \n";
    cout << "\n";

}

//Get the Gradient for the lattice params by transforming the lattice vectors.
Vector Cluster::LatticeParamsGradFromVectorsGrad(Vector VectorsGrad){ //JLM

  //Check to make sure that VectorsGrad is the right size
  if(VectorsGrad.GetLength() != 9){
    printf("Error::LatticeParamsGradFromVectorsGrad() VectorsGrad is not length of 9. Length = %i\n",VectorsGrad.GetLength());
    exit(0);
  }
  
  //lattice params values
  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];
  double alpha = UnitCellAngles[0];
  double beta = UnitCellAngles[1];
  double gamma =  UnitCellAngles[2];


  /*cout << "a = " << a << "\n"; 
  cout << "b = " << b << "\n";
  cout << "c = " << c << "\n";
  cout << "alpha = " << alpha << "\n";
  cout << "beta = " << beta << "\n";
  cout << "gamma = " << gamma << "\n";
  cout << "\n";*/


  a *= AngToBohr;
  b *= AngToBohr;
  c *= AngToBohr;          
  alpha *= DegreesToRadians;
  beta *= DegreesToRadians;
  gamma *= DegreesToRadians;


  //lattice vector values
  double dE_dx1 = VectorsGrad[0];
  double dE_dy1 = VectorsGrad[1];
  double dE_dz1 = VectorsGrad[2];
  double dE_dx2 = VectorsGrad[3];
  double dE_dy2 = VectorsGrad[4];
  double dE_dz2 = VectorsGrad[5];
  double dE_dx3 = VectorsGrad[6];
  double dE_dy3 = VectorsGrad[7];
  double dE_dz3 = VectorsGrad[8];

  //preform the transformation
  Vector ParamsGrad(6);

  ParamsGrad[0] = dE_dx1;
  ParamsGrad[1] = dE_dx2*cos(gamma) + dE_dy2*sin(gamma);  
  ParamsGrad[2] = cos(beta)*( dE_dx3 - dE_dy3 * cos(gamma)/sin(gamma) )+ dE_dy3 * cos(alpha)/sin(gamma)
		   + dE_dz3 * sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)/sin(gamma)/sin(gamma)-cos(alpha)*cos(alpha)/sin(gamma)/sin(gamma)-cos(beta)*cos(beta)/sin(gamma)/sin(gamma) );

  ParamsGrad[3] = c /sin(gamma)*sin(alpha) * ( -dE_dy3 + dE_dz3 * ( cos(alpha)-cos(beta)*cos(gamma) ) 
						/ sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) ) );
    
  ParamsGrad[4] = c*sin(beta) * (-dE_dx3 + dE_dy3 * cos(gamma)/sin(gamma)+ dE_dz3 * ( cos(beta)-cos(alpha)*cos(gamma) ) / sin(gamma)
				  / sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) ) );
    
  ParamsGrad[5] = b * dE_dy2 * cos(gamma) + c * dE_dy3 /sin(gamma)* ( -cos(alpha)*cos(gamma)/sin(gamma) + cos(beta)/sin(gamma) ) - b * dE_dx2 *sin(gamma)
    + (c * dE_dz3 * ( cos(alpha)*cos(alpha)*cos(gamma) + cos(beta)*cos(beta)*cos(gamma) - 0.5*cos(alpha)*cos(beta)*( 3 + cos(2*gamma) ) )  /sin(gamma)/sin(gamma) )
    /  sqrt( -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)+2*cos(alpha)*cos(beta)*cos(gamma)+sin(gamma)*sin(gamma) );

  return ParamsGrad;
}

//if xyz coordinates of atom are frozen or locked by symmetry
//alter symmetry so it is maintained.
void Cluster::MaintainCartesianSymmetry(Vector& Coord,bool LatticeIncluded){

  double tolerance = Params::Parameters().GetSymmetryTolerance() ;

  //Subtract the lattice length from N if they are included in the coord vector
  int N = Coord.GetLength()/3;;
  if(Params::Parameters().UseLatticeSymmetry() || Params::Parameters().UseSpaceSymmetry()){

  if(LatticeIncluded){
    N -= 2;
    
    if (Params::Parameters().UseLatticeSymmetry()) {
    //printf("Using lattice symmetry in MaintainCartesianSymmetry!\n\n"); //JLM

      if (Params::Parameters().FreezeUnitCellParams() ) { 
       for (int i=0; i<3;i++) {    
          Coord[3*N+i] = UnitCellAxes[i];
          Coord[3*N+3+i] = UnitCellAngles[i] ;
        }
      }

      //locking lattice axis by symmetry
      if(b_locked_by_a){
        Coord[3*N+1] = Coord[3*N];
        //UnitCellAxes[1] = UnitCellAxes[0];
      }
      if(c_locked_by_a){
        Coord[3*N+2] = Coord[3*N];
        //UnitCellAxes[2] = UnitCellAxes[0];
      }
      else if(c_locked_by_b){
        Coord[3*N+2] = Coord[3*N+1];
        //UnitCellAxes[2] = UnitCellAxes[1];
      }


      //if lattice parameters are 90 or 120, do not alter them
      if(lock_alpha)
        Coord[3*N+3] = UnitCellAngles[0] ;
        //UnitCellAngles[0] = coords[3*Natoms+3];
      if(lock_beta) 
        Coord[3*N+4] = UnitCellAngles[1] ;
        //UnitCellAngles[1] = coords[3*Natoms+4];
      if(lock_gamma)
        Coord[3*N+5] = UnitCellAngles[2] ;
        //UnitCellAngles[2] = coords[3*Natoms+5];
    
      //locking angles by symmetry
      if(beta_locked_by_alpha)
        Coord[3*N+4] = Coord[3*N+3];
        //UnitCellAngles[1] = UnitCellAngles[0];
      if(gamma_locked_by_alpha)
        Coord[3*N+5] = Coord[3*N+3];
        //UnitCellAngles[2] = UnitCellAngles[0];
      else if(gamma_locked_by_beta)
        Coord[3*N+6] = Coord[3*N+4];
        //UnitCellAngles[2] = UnitCellAngles[1];
     
   }
   //Coord.PrintGradient("Coord");
   //exit(0);
  }


  //int Coordkey[UniqueAtoms];
  int Monomerkey[UniqueAtoms];
  int Atomkey[UniqueAtoms];
  int counter = 0;
  for(int i=1; i<=NMon;i++){
    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
      if(Monomers[i].GetAtom(j).InAsymmetricUnit()){
	//Coordkey[counter] = Monomers[i].GetAtom(j).GetGlobalIndex();
	Monomerkey[counter] = i;
	Atomkey[counter] = j;
	counter++;
      }
    }
  }
  
  /*
  printf("Monomerkey = ");
  for(int i = 0;i<UniqueAtoms;i++)
    printf("%i ",Monomerkey[i]);
  printf("\n\n");

  printf("Atomkey = ");
  for(int i = 0;i<UniqueAtoms;i++)
    printf("%i ",Atomkey[i]);
  printf("\n\n");
  */

  for(int i = 0; i<N; i++){
    //int Index = Coordkey[i];
    int iMon = Monomerkey[i];
    int iAtom = Atomkey[i];

    //printf("i = %i iMon = %i iAtom = %i\n",i,iMon, iAtom);
    //printf("Change SignY %i Change SignZ %i\n",
    //   Monomers[iMon].GetAtom(iAtom).ChangeYSign(),
    //   Monomers[iMon].GetAtom(iAtom).ChangeZSign()); 
    
    
    if(Monomers[iMon].GetAtom(iAtom).IsAtomFrozen()){
      //Coord[3*i] = 0.0;
      //Coord[3*i+1] = 0.0;
      //Coord[3*i+2] = 0.0;
      Coord[3*i] = SymmetryUniqueCoordinates[3*i];
      Coord[3*i+1] = SymmetryUniqueCoordinates[3*i+1];
      Coord[3*i+2] = SymmetryUniqueCoordinates[3*i+2];
      
    }
    
    //SPECIAL CASE: Hexagonal cell with an atom is a fixed fractional x but with a cartesian coordinate with a y degree of freedom.
    //The updated x cartesian coordinate (which is not a degree of freedom) is shifted so y can be optimized but the x fractional coordinate
    // is a constant.
    //Only tested for R-3c phase IV CO2
    //delta_x = 1/tan(gamma) * delta_y where gamma = 120 degrees
    if((Monomers[iMon].GetAtom(iAtom).FreezeX()  || Monomers[iMon].GetAtom(iAtom).IsXLocked() == 1) && 
       !(Monomers[iMon].GetAtom(iAtom).FreezeY() || Monomers[iMon].GetAtom(iAtom).IsYLocked() ) &&
       (UnitCellAngles[2] - 120) < tolerance){
      double delta =  Coord[3*i+1] - SymmetryUniqueCoordinates[3*i+1];      
      if(Monomers[iMon].GetAtom(iAtom).ChangeXSign())
	delta = -delta;
      if(LatticeIncluded)
	delta *= 1/tan(Coord[3*N+5]*DegreesToRadians);
      else
	delta *= 1/tan(UnitCellAngles[2]*DegreesToRadians);
      Coord[3*i] = SymmetryUniqueCoordinates[3*i] + delta;      
    }
    //alter or freeze x under symmetry
    else if(Monomers[iMon].GetAtom(iAtom).FreezeX()){
      Coord[3*i] = SymmetryUniqueCoordinates[3*i];
   }else if(Monomers[iMon].GetAtom(iAtom).IsXLocked() == 2){
      double delta =  Coord[3*i+2] - SymmetryUniqueCoordinates[3*i+2];
      if(Monomers[iMon].GetAtom(iAtom).ChangeXSign())
	delta = -delta;
      Coord[3*i] = SymmetryUniqueCoordinates[3*i] + delta;  
    }
    else if(Monomers[iMon].GetAtom(iAtom).IsXLocked() == 1){
      double delta =  Coord[3*i+1] - SymmetryUniqueCoordinates[3*i+1];
      if(Monomers[iMon].GetAtom(iAtom).ChangeXSign())
	delta = -delta;
      Coord[3*i] = SymmetryUniqueCoordinates[3*i] + delta;
    }


    //alter or freeze y under symmetry
    if(Monomers[iMon].GetAtom(iAtom).FreezeY()){
      Coord[3*i+1] = SymmetryUniqueCoordinates[3*i+1];
    }
    else if(Monomers[iMon].GetAtom(iAtom).IsYLocked()){
      double delta =  Coord[3*i+2] - SymmetryUniqueCoordinates[3*i+2];
      if(Monomers[iMon].GetAtom(iAtom).ChangeYSign())
	delta = -delta;
      Coord[3*i+1] = SymmetryUniqueCoordinates[3*i+1] + delta;  	
    }

    //freeze z under symmetry
    if(Monomers[iMon].GetAtom(iAtom).FreezeZ()){
      Coord[3*i+2] = SymmetryUniqueCoordinates[3*i+2];
    }

    /*
    
    //alter or freeze y coordinates to maintain symmetry 
    if(Monomers[iMon].GetAtom(iAtom).FreezeY()){
      Coord[3*i+1] = SymmetryUniqueCoordinates[3*i+1];
    }else if(Monomers[iMon].GetAtom(iAtom).IsYLocked()){
      double delta =  Coord[3*i]- SymmetryUniqueCoordinates[3*i];
      if(Monomers[iMon].GetAtom(iAtom).ChangeYSign())
	delta = -delta;
      Coord[3*i+1] = SymmetryUniqueCoordinates[3*i+1] + delta;
    }

    //alter or freeze z coordinates to maintain symmetry 
    if(Monomers[iMon].GetAtom(iAtom).FreezeZ()){
      Coord[3*i+2] = SymmetryUniqueCoordinates[3*i+2];
    }
    else if(Monomers[iMon].GetAtom(iAtom).IsZLocked() == 1){
      double delta =  Coord[3*i]- SymmetryUniqueCoordinates[3*i];
      if(Monomers[iMon].GetAtom(iAtom).ChangeZSign())
	delta = -delta; 
      Coord[3*i+2] = SymmetryUniqueCoordinates[3*i+2] + delta;
    }
    else if(Monomers[iMon].GetAtom(iAtom).IsZLocked() == 1){
      double delta =  Coord[3*i+1]- SymmetryUniqueCoordinates[3*i+1];
      if(Monomers[iMon].GetAtom(iAtom).ChangeZSign())
	delta = -delta;
      Coord[3*i+2] = SymmetryUniqueCoordinates[3*i+2] + delta;
    }

    */
 
  //if atoms need to be shifted with lattice parameters to maintain symmetry do so known
  //Note that only atoms that are not a degree of freedom will this happen to
    if(LatticeIncluded){
      double ShiftX = 0.0;
      double ShiftY = 0.0;
      double ShiftZ = 0.0;
      if(fabs(Monomers[iMon].GetAtom(iAtom).ShiftByA()) > tolerance){
	ShiftX = Coord[3*N] - UnitCellAxes[0]; 
	ShiftX *= Monomers[iMon].GetAtom(iAtom).ShiftByA();
	Coord[3*i] += ShiftX;
      }
      if(fabs(Monomers[iMon].GetAtom(iAtom).ShiftByB()) > tolerance){
	ShiftY = Coord[3*N+1] - UnitCellAxes[1];
	ShiftY *= Monomers[iMon].GetAtom(iAtom).ShiftByB();
	Coord[3*i+1] += ShiftY;
      }
      if(fabs(Monomers[iMon].GetAtom(iAtom).ShiftByC()) > tolerance){
	ShiftZ = Coord[3*N+2] - UnitCellAxes[2];
	ShiftZ *= Monomers[iMon].GetAtom(iAtom).ShiftByC();
	Coord[3*i+2] += ShiftZ;
      }
      
      /*
      printf("\nChange a = %f b = %f c = %f\n",
	     Coord[3*N] - UnitCellAxes[0],
	     Coord[3*N+1] - UnitCellAxes[1],
	     Coord[3*N+2] - UnitCellAxes[2]);
      printf("Mon%i Atom%i ShiftX = %f ShiftY = %f ShiftZ = %f\n",
	     iMon,iAtom,ShiftX,ShiftY,ShiftZ);
      */
    }


  }
 

  } //Ending if statement to check is Symmetry is activated

  //Coord.PrintGradient("Coord after shift");
  //exit(0);


}




//Converts between Fraction and Cartesion. InFractional varable states whether the coordinates are currently in fractional coordinates. If it is convert to Cartesian,
//Otherwise convert into fraction coordinates.
//Can only be use for periodic systems
Vector Cluster::ConvertBetweenFractionAndCartesianCoordinates(Vector Coord, bool InFractional, bool GradIncludeLattice, bool Transpose){

  int Natom;

  //Only Tinker has the lattice parameters included in the coordinates
  if(GradIncludeLattice && Params::Parameters().GetMMType() != 2){
    //The last few values in Coord are either the lattice parameters or the gradient for the lattice parameters
    Natom = (Coord.GetLength()-6)/3;
  }else{
    Natom = Coord.GetLength()/3;
  }

  /*f
  //Lattice Params
  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];
  double alpha = UnitCellAngles[0]*DegreesToRadians;
  double beta = UnitCellAngles[1]*DegreesToRadians;
  double gamma = UnitCellAngles[2]*DegreesToRadians;

  double v  = sqrt(1 - pow(cos(alpha),2) - pow(cos(beta),2) - pow(cos(gamma),2) + 2*cos(alpha)*cos(beta)*cos(gamma));
  */

  //Matrix to convertion between the coordinate systems
  Matrix Convertion(3,3);

  //building matrix to convert from fraction to cartesian.
  Convertion.SetColumnVector(unit_cell[0],0);
  Convertion.SetColumnVector(unit_cell[1],1);
  Convertion.SetColumnVector(unit_cell[2],2);
  //Convertion.Print("Convertion");
  //Convertion.Print("Current lattice vectors");

  //use inverse if coordinates are in cartesian coordinates
  if(!InFractional){
    Convertion.Inverse();
  }

  //Transpose the transformation matrix. Needed for the transformation of the gradients 
  if(Transpose && InFractional){
    Convertion.Inverse();
    Convertion.Transpose();
  }else if(Transpose)
    Convertion.Transpose();
  /*
    if(InFractional){
    
    Convertion(0,0) = a;
    Convertion(0,1) = b*cos(gamma);
    Convertion(0,2) = c*cos(beta);
    Convertion(1,0) = 0;
    Convertion(1,1) = b*sin(gamma);
    Convertion(1,2) = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
    Convertion(2,0) = 0;
    Convertion(2,1) = 0;
    Convertion(2,2) = c*v/(sin(gamma));

  }
  else{
    Convertion(0,0) = 1/a;
    Convertion(0,1) = -cos(gamma)/(a*sin(gamma));
    Convertion(0,2) = (cos(alpha)*cos(gamma)-cos(beta))/(a*v*sin(gamma));
    Convertion(1,0) = 0;
    Convertion(1,1) = 1/(b*sin(gamma));
    Convertion(1,2) = (cos(beta)*cos(gamma)-cos(alpha))/(b*v*sin(gamma));
    Convertion(2,0) = 0;
    Convertion(2,1) = 0;
    Convertion(2,2) = sin(gamma)/(c*v);
  }
  */

  //Convertion.Print("Convertion");

  //storing Coord as a matrix since it is easier to do convertion with a matrix
  Matrix TmpCoord(3,Natom);
  for(int i=0;i<Natom;i++)
    for(int xyz=0;xyz<3;xyz++)
      TmpCoord(xyz,i) = Coord[3*i+xyz];

  //TmpCoord.Transpose();
  //TmpCoord.Print("Before TmpCoord");
  //TmpCoord.Transpose();

  //Convertion
  TmpCoord = Convertion.Multiply(TmpCoord);

  //Putting coodinates back into vector form
   for(int i=0;i<Natom;i++)
     for(int xyz=0;xyz<3;xyz++)
       Coord[3*i+xyz] = TmpCoord(xyz,i);

   return Coord;
  
}

//Alters the lattice vectors while preserving symmetry 
//Does not update the lattice parameters.
void Cluster::AlterLatticeVectors(int VectorEntry,double ChangedBy){

  //alterning Lattice Vectors
  // v1
  if(VectorEntry == 0)
    unit_cell[0][0] += ChangedBy;
  else if(VectorEntry == 1)
    unit_cell[0][1] += ChangedBy;
  else if(VectorEntry == 2)
    unit_cell[0][2] += ChangedBy;
  // v2
  else if(VectorEntry == 3)
    unit_cell[1][0] += ChangedBy;
  else if(VectorEntry == 4)
    unit_cell[1][1] += ChangedBy;
  else if(VectorEntry == 5)
    unit_cell[1][2] += ChangedBy;
  // v3
  else if(VectorEntry == 6)
    unit_cell[2][0] += ChangedBy;
  else if(VectorEntry == 7)
    unit_cell[2][1] += ChangedBy;
  else if(VectorEntry == 8)
    unit_cell[2][2] += ChangedBy;
  else{
    printf("Error::Cluster::AlterLatticeVectors() Lattice Vector Entry %i not recognized.\n",
	   VectorEntry);
    exit(0);
  }

  //Matrix to convertion between the coordinate systems
  Matrix Convertion(3,3);

  //building matrix to convert from fraction to cartesian.
  Convertion.SetColumnVector(unit_cell[0],0);
  Convertion.SetColumnVector(unit_cell[1],1);
  Convertion.SetColumnVector(unit_cell[2],2);
  //Convertion.Inverse();


    printf("\nUnit cell vectors: in Angstroms\n");
    printf("A: (%f,%f,%f)\n",unit_cell[0][0],unit_cell[0][1],unit_cell[0][2]);
    printf("B: (%f,%f,%f)\n",unit_cell[1][0],unit_cell[1][1],unit_cell[1][2]);
    printf("C: (%f,%f,%f)\n",unit_cell[2][0],unit_cell[2][1],unit_cell[2][2]);
    //shuhao: reciprocal cell vector
    printf("\nUnit reciprocal  cell vectors: in Angstroms^-1\n");
    printf("A*: (%f,%f,%f)\n",reciprocal_cell[0][0],reciprocal_cell[0][1],reciprocal_cell[0][2]);
    printf("B*: (%f,%f,%f)\n",reciprocal_cell[1][0],reciprocal_cell[1][1],reciprocal_cell[1][2]);
    printf("C*: (%f,%f,%f)\n",reciprocal_cell[2][0],reciprocal_cell[2][1],reciprocal_cell[2][2]);


  //Convertion.Print("Convertion");

  int Natom = GetTotalNumberOfAtoms();

  //Getting the coordinates for the rest of the system base on symmetry
  Vector Coords = GetSymmetryUniqueCoordinates();
  Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
  Coords = GetSymmetryImposedCoordinates(Coords,true);
  Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
  Coords.PrintGradient("Coords");

  /*
  //storing Coord as a matrix since it is easier to do convertion with a matrix
  Matrix TmpCoord(3,Natom);
  for(int i=0;i<Natom;i++)
    for(int xyz=0;xyz<3;xyz++)
      TmpCoord(xyz,i) = Coord[3*i+xyz];

  TmpCoord.Transpose();
  TmpCoord.Print("Fract Coord");
  TmpCoord.Transpose();

  //Convertion
  TmpCoord = Convertion.Multiply(TmpCoord);  

  //Putting coodinates back into vector form
   for(int i=0;i<Natom;i++)
     for(int xyz=0;xyz<3;xyz++)
       Coord[3*i+xyz] = TmpCoord(xyz,i);

  TmpCoord.Transpose();
  TmpCoord.Print("Coord");
  TmpCoord.Transpose();
  */

   SetNewCoordinates(Coords);

}

void Cluster::ChangeCellVolume(){
  
  Vector Coords = GetSymmetryUniqueCoordinates(true);
  
  
  Coords.PrintGradient("Symmetry Imposed Coord before changing volume");
  
  
  //if the lattice parameters are fixed, unfix and refix after changing volume
  bool fixed_lattice = 0;
  if(Params::Parameters().FreezeUnitCellParams()){
    fixed_lattice = 1;
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","0");
  }
  
  if(!(Params::Parameters().GetMMType() == 2)){
    
    if(fabs(Params::Parameters().ChangeVolumeBy()) > 0.000001){
      double new_volume = cell_volume + Params::Parameters().ChangeVolumeBy();
      
      printf("\nChanging cell volume by %f\n",
	     Params::Parameters().ChangeVolumeBy()); 
      
      for(int i=0;i<3;i++)
	Coords[3*UniqueAtoms+i] = UnitCellAxes[i]*pow(new_volume/cell_volume,1.0/3.0);
    }
    
    
    if(fabs(Params::Parameters().ChangeParameterA()) > 0.000001){
      
      printf("\nChanging a  by %f\n",
	     Params::Parameters().ChangeParameterA()); 
      
      Coords[3*UniqueAtoms] += Params::Parameters().ChangeParameterA(); 
      
    }
    
    if(fabs(Params::Parameters().ChangeParameterB()) > 0.000001){
      printf("\nChanging b  by %f\n",
	     Params::Parameters().ChangeParameterA()); 
      Coords[3*UniqueAtoms+1] += Params::Parameters().ChangeParameterB(); 
      
    }
    
    if(fabs(Params::Parameters().ChangeParameterC()) > 0.000001){
      printf("\nChanging c  by %f\n",
	     Params::Parameters().ChangeParameterC()); 
      Coords[3*UniqueAtoms+2] += Params::Parameters().ChangeParameterC(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterAlpha()) > 0.000001){
      printf("\nChanging alpha  by %f\n",
	     Params::Parameters().ChangeParameterAlpha()); 
      Coords[3*UniqueAtoms+3] += Params::Parameters().ChangeParameterAlpha(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterBeta()) > 0.000001){
      printf("\nChanging beta  by %f\n",
	     Params::Parameters().ChangeParameterBeta()); 
      Coords[3*UniqueAtoms+4] += Params::Parameters().ChangeParameterBeta(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterGamma()) > 0.000001){
      printf("\nChanging gamma  by %f\n",
	     Params::Parameters().ChangeParameterBeta()); 
      Coords[3*UniqueAtoms+5] += Params::Parameters().ChangeParameterGamma(); 
    }
    
    
    
    //Setting the new lattice parameters and Altering volume coordinating to symmetry
    Cluster::MaintainCartesianSymmetry(Coords,true);
    Cluster::cluster().UpdateLatticeParams(Coords);
    Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
    Coords = GetSymmetryImposedCoordinates(Coords,true);
    Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
    
    
  }else{
    
    if(fabs(Params::Parameters().ChangeVolumeBy()) > 0.000001){
      double new_volume = cell_volume + Params::Parameters().ChangeVolumeBy();
      
      printf("\nChanging cell volume by %f\n",
	     Params::Parameters().ChangeVolumeBy()); 
      
      for(int i=0;i<3;i++)
	UnitCellAxes[i] = UnitCellAxes[i]*pow(new_volume/cell_volume,1.0/3.0);
    }
    
    
    if(fabs(Params::Parameters().ChangeParameterA()) > 0.000001){
      
      printf("\nChanging a  by %f\n",
	     Params::Parameters().ChangeParameterA()); 
      
      UnitCellAxes[0] += Params::Parameters().ChangeParameterA(); 
      
    }
    
    if(fabs(Params::Parameters().ChangeParameterB()) > 0.000001){
      printf("\nChanging b  by %f\n",
	     Params::Parameters().ChangeParameterA()); 
      UnitCellAxes[1] += Params::Parameters().ChangeParameterB(); 
      
    }
    
    if(fabs(Params::Parameters().ChangeParameterC()) > 0.000001){
      printf("\nChanging c  by %f\n",
	     Params::Parameters().ChangeParameterC()); 
      UnitCellAxes[2] += Params::Parameters().ChangeParameterC(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterAlpha()) > 0.000001){
      printf("\nChanging alpha  by %f\n",
	     Params::Parameters().ChangeParameterAlpha()); 
      UnitCellAngles[0] += Params::Parameters().ChangeParameterAlpha(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterBeta()) > 0.000001){
      printf("\nChanging beta  by %f\n",
	     Params::Parameters().ChangeParameterBeta()); 
      UnitCellAngles[1]  += Params::Parameters().ChangeParameterBeta(); 
    }
    
    if(fabs(Params::Parameters().ChangeParameterGamma()) > 0.000001){
      printf("\nChanging gamma  by %f\n",
	     Params::Parameters().ChangeParameterBeta()); 
      UnitCellAngles[2] += Params::Parameters().ChangeParameterGamma(); 
    }
    
    Cluster::MaintainCartesianSymmetry(Coords,false);
    
    //This part usually handled by UpdateLatticeParams but since Coord doesn't contain lattice parameters in AIFF, it is done here
    double a = UnitCellAxes[0];
    double b = UnitCellAxes[1];
    double c = UnitCellAxes[2];
    double alpha = UnitCellAngles[0];
    double beta = UnitCellAngles[1];
    double gamma = UnitCellAngles[2];
    double alpha_rad = alpha*DegreesToRadians;
    double beta_rad = beta*DegreesToRadians;
    double gamma_rad = gamma*DegreesToRadians;
    double beta_term =
      (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
    double gamma_term =
      sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);
    
    /*
      cout << "after \n" <<"\n";
      cout << "a = " << a << "\n"; 
      cout << "b = " << b << "\n";
      cout << "c = " << c << "\n";
      cout << "alpha = " << alpha << "\n";
      cout << "beta = " << beta << "\n";
      cout << "gamma = " << gamma << "\n";
      fflush(stdout);
    */    
    
    
    
    // v1
    unit_cell[0][0] = a;
    unit_cell[0][1] = 0;
    unit_cell[0][2] = 0;
    
    // v2
    unit_cell[1][0] = b*cos(gamma_rad);
    unit_cell[1][1] = b*sin(gamma_rad);             
    unit_cell[1][2] = 0;
    
    // v3
    unit_cell[2][0] = c*cos(beta_rad);
    unit_cell[2][1] = c*beta_term;    
    unit_cell[2][2] = c*gamma_term;
    //}
    
    cell_volume = unit_cell[0][0]*unit_cell[1][1]*unit_cell[2][2];

    cout << "a = " << UnitCellAxes[0] << "\n"; 
    cout << "b = " << UnitCellAxes[1] << "\n";
    cout << "c = " << UnitCellAxes[2] << "\n";
    cout << "alpha = " << UnitCellAngles[0] << "\n";
    cout << "beta = " << UnitCellAngles[1] << "\n";
    cout << "gamma = " << UnitCellAngles[2] << "\n";
    cout << "cell_volume = " << cell_volume << " \n";
    cout << "\n";

    Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,false,false);
    Coords = GetSymmetryImposedCoordinates(Coords,true,false);
    Coords = ConvertBetweenFractionAndCartesianCoordinates(Coords,true,false);

  }

  //Coords.PrintGradient("Symmetry Imposed Coord");

  
  //fix lattice parameters were fixed, refixing them.
  if(fixed_lattice)
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
  
  //Changing coordinates to preserve symmetry
  SetNewCoordinates(Coords);
  //Coords.PrintGradient("Coord under new volume");

}


//this function translates the center of mass of monomers when the lattice parameters changes while maintaining symmetry
void Cluster::TranslateUsingLatticeParams(Vector OldUnitCellAngles,Vector OldUnitCellAxes){
  
  //The Old Lattice Params 
  double a = OldUnitCellAxes[0];
  double b = OldUnitCellAxes[1];
  double c = OldUnitCellAxes[2];
  double alpha = OldUnitCellAngles[0]*DegreesToRadians;
  double beta = OldUnitCellAngles[1]*DegreesToRadians;
  double gamma = OldUnitCellAngles[2]*DegreesToRadians;
  double v  = sqrt(1 - pow(cos(alpha),2) - pow(cos(beta),2) - pow(cos(gamma),2) + 2*cos(alpha)*cos(beta)*cos(gamma));

  //Conversion from cartesian to fractional coordinates
  Matrix CartToFrac(3,3);
  CartToFrac(0,0) = 1/a;
  CartToFrac(0,1) = -cos(gamma)/(a*sin(gamma));
  CartToFrac(0,2) = (cos(alpha)*cos(gamma)-cos(beta))/(a*v*sin(gamma));
  CartToFrac(1,0) = 0;
  CartToFrac(1,1) = 1/(b*sin(gamma));
  CartToFrac(1,2) = (cos(beta)*cos(gamma)-cos(alpha))/(b*v*sin(gamma));
  CartToFrac(2,0) = 0;
  CartToFrac(2,1) = 0;
  CartToFrac(2,2) = sin(gamma)/(c*v);
  //Center of monomer coordinates for each monomer 
  Matrix MonomerCOM(3,NMon);
  for(int i=0;i<NMon;i++){
    MonomerCOM(0,i) = Monomers[i+1].GetCenterOfMass(0);
    MonomerCOM(1,i) = Monomers[i+1].GetCenterOfMass(1);
    MonomerCOM(2,i) = Monomers[i+1].GetCenterOfMass(2);
  }

  //store for later
  Matrix OldCOM = MonomerCOM;

  printf("Center of Mass of Monomers before shifting\n");
  for(int i=0;i<NMon;i++){
    printf("m%i %f %f %f\n",
	   i+1,MonomerCOM(0,i),MonomerCOM(1,i),MonomerCOM(2,i));
  } 

  MonomerCOM = CartToFrac.Multiply(MonomerCOM);


  //Changing The Lattice Params to the New Ones 
  a = UnitCellAxes[0];
  b = UnitCellAxes[1];
  c = UnitCellAxes[2];
  alpha = UnitCellAngles[0]*DegreesToRadians;
  beta = UnitCellAngles[1]*DegreesToRadians;
  gamma = UnitCellAngles[2]*DegreesToRadians;
  v  = sqrt(1 - pow(cos(alpha),2) - pow(cos(beta),2) - pow(cos(gamma),2) + 2*cos(alpha)*cos(beta)*cos(gamma));
  
  //Conversion from fractional to cartesian coordinates
  Matrix FracToCart(3,3);
  FracToCart(0,0) = a;
  FracToCart(0,1) = b*cos(gamma);
  FracToCart(0,2) = c*cos(beta);
  FracToCart(1,0) = 0;
  FracToCart(1,1) = b*sin(gamma);
  FracToCart(1,2) = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
  FracToCart(2,0) = 0;
  FracToCart(2,1) = 0;
  FracToCart(2,2) = c*v/(sin(gamma));

  MonomerCOM = FracToCart.Multiply(MonomerCOM);
  Matrix NewCOM = MonomerCOM;

  printf("Center of Mass of Monomers after shifting\n");
  for(int i=0;i<NMon;i++){
    printf("m%i %f %f %f\n",
	   i+1,MonomerCOM(0,i),MonomerCOM(1,i),MonomerCOM(2,i));
  } 

  //Shifting monomers so that all the symmetrical unique monomers coordinates change while symmetry is preserved
  for(int i=0;i<NMon;i++){
    int SymMon = Monomers[i+1].GetSymmetricalMonomer() - 1;
    double ShiftBack;
    for(int xyz=0;xyz<3;xyz++){
      ShiftBack = NewCOM(xyz,SymMon) - OldCOM(xyz,SymMon);
      // MonomerCOM(xyz,i) = NewCOM(xyz,i) - ShiftBack;
    }
    
    
  }

  //shifting the center of mass of all the monomers
  for(int i=0;i<NMon;i++) 
    Monomers[i+1].Translate(MonomerCOM.GetColumnVector(i));

  for(int i=1;i<=NMon;i++){
    printf("\nm%i\n",i);
    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++)
      printf("%s %f %f %f\n",Monomers[i].GetSymbol(j).c_str(), Monomers[i].GetAtom(j).GetPosition(0),
	     Monomers[i].GetAtom(j).GetPosition(1), Monomers[i].GetAtom(j).GetPosition(2));
  }

  ReadCurrentCoordinates();
}

// scalar triple product of a,b,c gives the volume.
// Basically the determinant of vectors a,b,c
double Cluster::UnitCellVolume() {

  /*
  double volume=0.0;
  Vector abox(3), bbox(3), cbox(3); 
  abox[0] =  unit_cell[0][0];
  abox[1] =  unit_cell[0][1]; 
  abox[2] =  unit_cell[0][2]; 
  bbox[0] =  unit_cell[1][0]; 
  bbox[1] =  unit_cell[1][1]; 
  bbox[2] =  unit_cell[1][2]; 
  cbox[0] =  unit_cell[2][0]; 
  cbox[1] =  unit_cell[2][1]; 
  cbox[2] =  unit_cell[2][2]; 
//  abox.PrintGradient("abox");
//  bbox.PrintGradient("bbox");
//  cbox.PrintGradient("cbox");

//  cout << "abox while calculating volume = " << abox[0] << "\n";
//  cout << "Initially, volume = " << volume <<"\n";
  volume = (abox.CrossProduct(bbox)).DotProduct(cbox);
//  abox.CrossProduct(bbox).PrintGradient("a x b");
  cout << "volume = " << volume << "\n";



  return volume;*/
  return cell_volume;

}


// Similar to CreatePeriodicImageDimerList, but it creates a list of monomers
// for inducing multipoles on central unit cell molecules.  
void Cluster::CreatePeriodicImageMonomerList(double r_cutoff) {
  // Identify how far we have to go along each unit cell direction to 
  // stay within the cutoff.  Add 1 extra image cell to each, for good measure.
  //double r_cutoff = Params::Parameters().GetMaxPolarizationRadius(); 
  int Nv[3] = {1,1,1}; // start with 1 image cell in each direction
  for (int i=0;i<3;i++) {
      double dist = 0;
      while (dist  < r_cutoff) {
	dist +=  unit_cell[i].Norm();
	Nv[i] += 1;
      }
      //printf("Nv[%d] = %d, dist = %f\n",i,Nv[i],dist);
    }

  // Prepare output xyz file for visualizing
  // Open the input file for writing
  FILE *xyz;
  string filename = "image_monomers.xyz";
  if ((xyz = fopen(filename.c_str(),"w"))==NULL) {
    printf("Cluster::CreatePeriodicImageMonomerList() : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  fprintf(xyz,"XXX\n\n"); // we set the number of atoms later

  // print central unit cell coords
  for (int i=1;i<=NMon;i++) 
    Monomers[i].PrintMonomerCartesian(xyz);
  int active_atoms = GetTotalNumberOfAtoms();



  // Create dynamic 6*N array, access as ImageList[i][j].  The j-index
  // is 4 items long: the Lattice Vector (3 items) and the reference
  // monomer that has been replicated.  The other dimension, N grows
  // dynamically to accomodate as many entries as are needed to store
  // the list that later generates all Images dimers.
  vector < vector<int> > ImageList(1,vector<int>(4,0));
  
  /*
    Start generating image cells.  Use a 2-step algorithm:  

    Pass 1: identifies the monomers we need to include.  These results
    are stored in ImageList.

    Pass 2: uses ImageList to construct the necessary dimers.  

    This is probably a slow algorithm, especially since it dynamically
    grows the ImageList.  But, Pass 1 lets us work with fairly small
    amounts of memory and to check for symmetry, etc, and to find how
    many dimers we actually need.  Then, we do the real part in Pas	double damp = DimerImages[j].GetDampingFactor(c0,c1);s 2
    with the pruned list of known size.

  */

  // Pass 1
  bool skip; // flag for ignoring some.
  int kept_it = 0;
  // Loop over the image cells, in both positive & negative directions
  for (int x=-Nv[0];x<=Nv[0];x++) 
    for (int y=-Nv[1];y<=Nv[1];y++)
      for (int z=-Nv[2];z<=Nv[2];z++) {
	skip = false;

	if (x==0 && y==0 && z==0) {
	  //printf("Skipping central unit cell, x=%d, y=%d, z=%d\n",x,y,z);
	  skip = true;
	}

	if (!skip) {
	  // Create copy of Monomers that we can translate as needed;
	  Monomer* ImageMonomers = new Monomer[NMon+1];
	  for (int i=1;i<=NMon;i++) {
	    ImageMonomers[i] = Monomers[i];
	    
	    //printf("ImageMonomer %d before shifting:\n",i);
	    //ImageMonomers[i].PrintAll();
	  }


	  // Determine the shift (translation vector) from the central
	  // cell to the image cell
	  Vector shift(3);
	  for (int i=0;i<3;i++) {
	    shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	      z*unit_cell[2][i];
	  }

	  // Translate the monomers from the central unit cell to the
	  // image cell
	  for (int imon=1;imon<=NMon;imon++) {

	    bool KeepThisImageMonomer = false;

	    Vector new_com(3);
	    new_com = ImageMonomers[imon].GetCenterOfMass();
	    
	    // Add the shift and translate the monomer
	    new_com += shift;
	    ImageMonomers[imon].Translate(new_com);

	    // Now pair this image monomer with each monomer in central unit 
	    // cell and test distance relative to the cutoff
	    for (int jmon=1;jmon<=NMon;jmon++) {

	      Dimer Tmp;
	      //printf("jmon = %d, imon = %d\n",jmon,imon);
	      Tmp.Initialize(Monomers[jmon],ImageMonomers[imon]);
	      //printf("dimer from %d and image %d\n",jmon,imon);
	      //Tmp.PrintQChemCartesian();
	      if ( Tmp.GetDimerSeparation() < r_cutoff ) {
		if ( Params::Parameters().PrintLevel() > 3 ) {
		  printf("K_vec = (%d,%d,%d)\n",x,y,z);
		  printf("(%d*,%d) Separation = %f. Keeping it.\n",
			 imon,jmon,Tmp.GetDimerSeparation() );
		}
		KeepThisImageMonomer = true;

	      }
	    }
	    if ( KeepThisImageMonomer )  {
	      // Store info to later regenerate this list
	      int tmp[] = {x,y,z,imon};
	      vector<int> tmpvec (tmp,tmp+4);
	      if (kept_it==0)
		ImageList[0] = tmpvec; // first time, modify existing row
	      else {
		ImageList.push_back(tmpvec); // add row to list
	      }
	      kept_it++;  //counts how many image monomers we have kept
	      
	    }

	    // Print out the coordinates to the file for visualization.
	    if (KeepThisImageMonomer) {
	      ImageMonomers[imon].PrintMonomerCartesian(xyz);
	      active_atoms += ImageMonomers[imon].GetNumberOfAtoms();
	    }
	  }
	  delete [] ImageMonomers;
	}	  
      }
  //printf("Total number of image monomers = %d\n",kept_it);

  if ( Params::Parameters().PrintLevel() ) {
    printf("List of saved image monomers\n");
    for (int i=0;i<ImageList.size();i++) {
      printf("%2d: (%d,%d,%d) of Mon %d\n",i,ImageList[i][0],
	     ImageList[i][1],ImageList[i][2],ImageList[i][3]);
    }
  }

  // Close up geometry file
  fclose(xyz);
  // update number of atoms;
  string cmd = "sed -i -e s/XXX/";
  char count[10];
  sprintf(count,"%d",active_atoms);
  cmd += count;
  cmd += "/g " + filename;
  //printf("Running: %s\n",cmd.c_str());
  system(cmd.c_str());
  

  // Pass 2: Now actually create image monomers using the ImageList
  // from Pass 1.

  MonomerImages = new Monomer[kept_it+1];

  int Mon_index = NMon;  // Monomer counter, image counters start from NMon+1
  int imon;
  for (imon=1;imon<=kept_it;imon++) {
    Mon_index++;
    // Read data from ImageList
    int x = ImageList[imon-1][0];
    int y = ImageList[imon-1][1];
    int z = ImageList[imon-1][2];
    int image_mon = ImageList[imon-1][3];

    // Create copy of Monomer that we can translate
    Monomer ImageMon;
    ImageMon = Monomers[image_mon];

    // Determine the shift from the central cell to the image cell
    Vector shift(3);
    for (int i=0;i<3;i++) {
      shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	z*unit_cell[2][i];
    }


    // Translate the monomer    
    Vector new_com(3);
    new_com = ImageMon.GetCenterOfMass();
    new_com += shift;
    ImageMon.Translate(new_com);
    if (Params::Parameters().PrintLevel() > 3) 
      printf("K_vec = (%d,%d,%d)\n",x,y,z);
    
    // Renumber it & set reference monomer index
    ImageMon.SetReferenceMonomerIndex(ImageMon.GetIndex());
    ImageMon.SetIndex(Mon_index); 
    
    MonomerImages[imon] = ImageMon;

    // JDH CHECK THIS:
    MonomerImages[imon].SetUnitCellIndex(1,x);
    MonomerImages[imon].SetUnitCellIndex(2,y);
    MonomerImages[imon].SetUnitCellIndex(3,z);
    if ( Params::Parameters().UseElectrostaticEmbedding()  ) { //JDH
      for (int iatom = 0; iatom < MonomerImages[imon].GetNumberOfAtoms() ; iatom++ ) {
	MonomerImages[imon].GetAtom(iatom).SetMultipoleMoments( Monomers[ MonomerImages[imon].GetReferenceMonomerIndex() ].GetAtom(iatom).GetMultipoleMoments() );
      }
    }
    
 
  }
  
  printf("  %d image monomers created\n",imon-1);
  NMon_images = imon-1;

} 

void Cluster::DebugLatticeGradient() {

  // Grab the unit cell parameters
  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];

  double alpha = UnitCellAngles[0];
  double beta = UnitCellAngles[1];
  double gamma = UnitCellAngles[2];

  printf("\nUnit cell parameters:\n");
  printf("\tA = %.4f     B = %.4f    C = %.4f\n",a,b,c);
  printf("\talpha = %.2f beta = %.2f gamma = %.2f\n\n",alpha,beta,gamma);

  // Convert distances to Bohr, angles to radians
  a *= AngToBohr;
  b *= AngToBohr;
  c *= AngToBohr;
  alpha *= DegreesToRadians;
  beta *= DegreesToRadians;
  gamma *= DegreesToRadians;
  
  printf("Unit cell parameters after unit conversion (now in bohr & radians):\n");
  printf("\tA = %.4f     B = %.4f    C = %.4f\n",a,b,c);
  printf("\talpha = %.4f beta = %.4f gamma = %.4f\n\n",alpha,beta,gamma);


  // Compute beta_term and gamma_terms
  double beta_term = (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma);
  double gamma_term = sqrt(1.0 - cos(beta)*cos(beta) - beta_term*beta_term);
  printf("beta_term = %f, gamma_term = %f\n",beta_term,gamma_term);

  // Compute gradients of beta_term w.r.t. angles
  double dB_dalpha = - sin(alpha) / sin(gamma);
  double dB_dbeta = sin(beta)*cos(gamma)/sin(gamma);
  double dB_dgamma = cos(beta) - beta_term*cos(gamma)/sin(gamma);
  printf("dB/dalpha = %f, dB/dbeta = %f, dB/dgamma = %f\n",dB_dalpha, dB_dbeta, dB_dgamma);


  // Compute gradients of gamma_Term w.r.t. angles
  double dG_dalpha = - beta_term * dB_dalpha / gamma_term;
  double dG_dbeta = ( cos(beta)*sin(beta) - beta_term*dB_dbeta ) / gamma_term;
  double dG_dgamma = - beta_term*dB_dgamma / gamma_term;
  printf("dG/dalpha = %f, dG/dbeta = %f, dG/dgamma = %f\n",dG_dalpha, dG_dbeta, dG_dgamma);


  // Now compute the coefficients for each term in the overall gradients.
  // dE_da:
  double dE_da_x1 = 1.0;

  // dE_db:
  double dE_db_x2 = cos(gamma);
  double dE_db_y2 = sin(gamma);

  // dE_dc:
  double dE_dc_x3 = cos(beta);
  double dE_dc_y3 = beta_term;
  double dE_dc_z3 = gamma_term;

  // dE_dalpha:
  double dE_dalpha_y3 = c*dB_dalpha;
  double dE_dalpha_z3 = c*dG_dalpha;

  // dE_dbeta:
  double dE_dbeta_x3 = -c*sin(beta);
  double dE_dbeta_y3 = c*dB_dbeta;
  double dE_dbeta_z3 = c*dG_dbeta;

  // dE_dgamma:
  double dE_dgamma_x2 = -b*sin(gamma);
  double dE_dgamma_y2 = b*cos(gamma);
  double dE_dgamma_y3 = c*dB_dgamma;
  double dE_dgamma_z3 = c*dG_dgamma;

  // Print out the results:

  printf("\nGradient      Coeff\n");
  printf("              Ex      Ey      Ez\n");
  printf("dE/da     = %10.6f (dE/dx1)\n",dE_da_x1);
  printf("dE/db     = %10.6f (dE/dx2) + %10.6f (dE/dy2)\n",dE_db_x2, dE_db_y2);
  printf("dE/dc     = %10.6f (dE/dx3) + %10.6f (dE/dy3) + %10.6f (dE/dz3)\n",
	 dE_dc_x3, dE_dc_y3, dE_dc_z3);
  printf("\n");
  printf("dE/dalpha = %10.6f (dE/dy3) + %10.6f (dE/dz3)\n", dE_dalpha_y3, dE_dalpha_z3);
  printf("dE/dbeta  = %10.6f (dE/dx3) + %10.6f (dE/dy3) + %10.6f (dE/dz3)\n",
	 dE_dbeta_x3, dE_dbeta_y3, dE_dbeta_z3);
  printf("dE/dgamma = %10.6f (dE/dx2) + %10.6f (dE/dy2) + %10.6f (dE/dy3) + %10.6f (dE/dz3)\n",
	 dE_dgamma_x2, dE_dgamma_y2, dE_dgamma_y3, dE_dgamma_z3);
  printf("\n");


  // Compute the nine gradient components (dE/dq), where q is any of x1, y1, z1, x2, y2, z2, 
  // x3, y3, or z3.

  Vector Grad(6); // Initialize the vector for the final gradient;
  Grad.Set();

  // Matrix for storing dE/dx1, dE/dy1, etc. Array column 1 -> x1,y1,z1, column 2 -> x2,y2,z2, etc.
  Matrix dEdq(3,3); 
  dEdq.Set();

  // grab the cutoffs for the QM -> MM damping
  double c0 = Params::Parameters().GetLocalCutoff(0);
  double c1 = Params::Parameters().GetLocalCutoff(1);

  double scafac; 
  
  if (Params::Parameters().NeglectManyBody() || Params::Parameters().TinkerDebug() ) {
    scafac = 0.0; // neglect MM many-body terms
  }
  else {
    scafac = 1.0; // treat them as usual
  }

  for (int i=1;i<=NDim_images;i++) {
    int indexA = DimerImages[i].GetIndexA();
    int indexB = DimerImages[i].GetIndexB();
    int indexB_ref = DimerImages[i].GetReferenceMonomerIndex();
    int Na = DimerImages[i].GetMonomerA().GetNumberOfAtoms();
    int Nb = DimerImages[i].GetMonomerB().GetNumberOfAtoms();

    int minA = DimerImages[i].GetMinimumDistanceAtomA();
    int minB = DimerImages[i].GetMinimumDistanceAtomB();

    double damp = DimerImages[i].GetDampingFactor(c0,c1);
    double R = DimerImages[i].GetDimerSeparation();
    //printf("DimerImage %d, R = %f, damping factor = %f\n",i, R, damp); 
    double symfac = 1.0; 
    if ( Params::Parameters().UseCrystalSymmetry() )
      symfac *= DimerImages[i].GetSymmetryFactor();

    // 1) Compute the derivatives dE/dQ w.r.t. global X,Y,Z.  dE/dX, dE/dY, dE/dZ
    // dE/dQ = sum_i sum_k sum_atoms (dE(ik)/dQ){QM} - (dE(k)/dQ){QM} - (dE(ik)/dQ){MM} + (dE(k)/dQ){MM}
    
    double dEdX = 0.0;
    double dEdY = 0.0;
    double dEdZ = 0.0;
   
    /******************************************************************/
    // For extra debugging, do finite difference gradient of the
    // damping function with respect to x,y,z.  This can be compared
    // against the analytical version below, which we actually use.
    
    // Grab the Monomers in this dimer
    Monomer MonA = DimerImages[i].GetMonomerA();
    Monomer MonB = DimerImages[i].GetMonomerB();
    Monomer refB = Monomers[DimerImages[i].GetReferenceMonomerIndex()];
    
    // Initialize the gradient of the damping function
    Matrix dDamp_dQ(Nb,3);
    dDamp_dQ.Set();

    // loop over atoms in Monomer B
    for (int j=0;j<Nb;j++) {
	double delta = 0.001;

      // loop over x,y,z
      for (int q=0;q<3;q++) {
	// Shift atom j by +delta
	Vector xyz = MonB.GetAtom(j).GetPosition();
	xyz[q] += delta;
	MonB.GetAtom(j).SetPosition(xyz);
	DimerImages[i].UpdateObjects(MonA,MonB,refB);
	double damp_plus = DimerImages[i].GetDampingFactor(c0,c1);
	xyz[q] -= delta;
	MonB.GetAtom(j).SetPosition(xyz);
	
	// Shift atom j by -delta
	xyz[q] -= delta;
	MonB.GetAtom(j).SetPosition(xyz);
	DimerImages[i].UpdateObjects(MonA,MonB,refB);
	double damp_minus = DimerImages[i].GetDampingFactor(c0,c1);
	xyz[q] += delta;
	MonB.GetAtom(j).SetPosition(xyz);
	DimerImages[i].UpdateObjects(MonA,MonB,refB);

	dDamp_dQ(j,q) = (damp_plus - damp_minus)/(2*delta*AngToBohr);
      }
    }
    //dDamp_dQ.Print("dDamp_dQ");

    // End finite difference gradient of the damping function.
    /******************************************************************/

    
    // debugging: set dQM = 0.0 to turn off QM term, or 1.0 to do normal gradient
    double dQM = 1.0;

    // Compute gradient contribution arising from image dimer
    // interactions
    for (int j=0;j<Nb;j++) {
      int jx = 3*j;
      int jy = 3*j + 1;
      int jz = 3*j + 2;
      // create appropriate indices for dimer gradient, where offset to second monomer is needed.
      int djx = jx + 3*Na; 
      int djy = jy + 3*Na;
      int djz = jz + 3*Na;

      // Grab the positions of the atoms that determine the intermolecular separation
      double x_i = DimerImages[i].GetMonomerA().GetAtom(minA).GetPosition(0);
      double y_i = DimerImages[i].GetMonomerA().GetAtom(minA).GetPosition(1);
      double z_i = DimerImages[i].GetMonomerA().GetAtom(minA).GetPosition(2);

      double x_j = DimerImages[i].GetMonomerB().GetAtom(minB).GetPosition(0);
      double y_j = DimerImages[i].GetMonomerB().GetAtom(minB).GetPosition(1);
      double z_j = DimerImages[i].GetMonomerB().GetAtom(minB).GetPosition(2);


      // For the derivative of the damping factor term, we only get a
      // contribution if atom j corresponds to the atom involved
      // defining the minimum distance.  Use a Kronecker delta(j,minB)
      // for this:
      double delta = 0.0;
      if (j == minB) delta = 1.0;

      // Compute dDamp(R)/dQ
      double dDamp_dXj, dDamp_dYj, dDamp_dZj;
      if (R <= c1 || R >= c0) { // outside the damping region
	dDamp_dXj = 0.0;
	dDamp_dYj = 0.0;
	dDamp_dZj = 0.0;
      }
      // near the c1 cutoffs, can need to avoid overflow from huge exp() term
      else if (2.0*fabs(c1-c0)/(c1-R) - fabs(c1-c0)/(R-c0) > 50){ 
	double value = 2.0*fabs(c1-c0)/(c1-R) - fabs(c1-c0)/(R-c0);
	printf("Avoiding blow-up term in dDamp_dQ.  c1, c0, R = %f, %f, %f  and exponential = e^%f\n",c1,c0,R,value);
	// in this case, damp*damp*exp() is very nearly = damp, since damp = 1/(1+exp()).
	dDamp_dXj = delta*damp
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (x_i-x_j)/(R*AngToBohr);
	dDamp_dYj = delta*damp
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (y_i-y_j)/(R*AngToBohr);
	dDamp_dZj = delta*damp
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (z_i-z_j)/(R*AngToBohr);
      }
      else {
	// Note: when differentiating w.r.t. R, we need distances in
	// a.u., so we convert R in the final division from Ang to Bohr.
	dDamp_dXj = delta*damp*damp*exp(2.0*fabs(c1-c0)/(c1-R) - fabs(c1-c0)/(R-c0))
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (x_i-x_j)/(R*AngToBohr);
	dDamp_dYj = delta*damp*damp*exp(2.0*fabs(c1-c0)/(c1-R) - fabs(c1-c0)/(R-c0))
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (y_i-y_j)/(R*AngToBohr);
	dDamp_dZj = delta*damp*damp*exp(2.0*fabs(c1-c0)/(c1-R) - fabs(c1-c0)/(R-c0))
	  *( 2.0*fabs(c1-c0)/((c1-R)*(c1-R)) + fabs(c1-c0)/((R-c0)*(R-c0)) ) * (z_i-z_j)/(R*AngToBohr);
      }


      Monomer mB = Monomers[indexB_ref];
      // add contributions from each atom in this dimer
      dEdX += 0.5*symfac*damp*( dQM*( DimerImages[i].GetQMGradient()[djx] - mB.GetQMGradient()[jx] )
				- scafac*( DimerImages[i].GetMMGradient()[djx] - mB.GetMMGradient()[jx] ) )
	+ 0.5*symfac*dDamp_dXj*(dQM*DimerImages[i].GetQMIntEnergy() - scafac*DimerImages[i].GetMMIntEnergy());

      dEdY += 0.5*symfac*damp*( dQM*( DimerImages[i].GetQMGradient()[djy] - mB.GetQMGradient()[jy] )
				- scafac*( DimerImages[i].GetMMGradient()[djy] - mB.GetMMGradient()[jy] ) )
	+ 0.5*symfac*dDamp_dYj*(dQM*DimerImages[i].GetQMIntEnergy() - scafac*DimerImages[i].GetMMIntEnergy());    

      dEdZ += 0.5*symfac*damp*( dQM*( DimerImages[i].GetQMGradient()[djz] - mB.GetQMGradient()[jz] )
				- scafac*( DimerImages[i].GetMMGradient()[djz] - mB.GetMMGradient()[jz] ) )
	+ 0.5*symfac*dDamp_dZj*(dQM*DimerImages[i].GetQMIntEnergy() - scafac*DimerImages[i].GetMMIntEnergy());
    }
    //printf("End dEdX = %15.9f, dEdY = %15.9f, dEdZ = %15.9f\n",dEdX,dEdY,dEdZ);

    // 2) Relate dE/dQ to dE/dq via (dE/dx_i) = n_i (dE/dX), & similar for dE/dy_i, dE/dz_i
    // Components for derivatives along each component of the 3 unit cell vectors
    for (int vec=0;vec<3;vec++) {
      // add contributions for each image dimer
      int n = DimerImages[i].GetImageCell()[vec];
      //printf("i = %d, n(%d) = %d\n",i,vec,n);
      dEdq(0,vec) += (double) n*dEdX;
      dEdq(1,vec) += (double) n*dEdY;
      dEdq(2,vec) += (double) n*dEdZ;
    }
  }
  
  //dEdq.Print("dEdq =");

  // Now compute the gradient in terms of the 6 lattice parameters
  // dE/da
  Grad[0] = dE_da_x1*dEdq(0,0);
  // dE/db
  Grad[1] = dE_db_x2*dEdq(0,1) + dE_db_y2*dEdq(1,1);
  // dE/dc
  Grad[2] = dE_dc_x3*dEdq(0,2) + dE_dc_y3*dEdq(1,2) + dE_dc_z3*dEdq(2,2);
  // dE/dalpha
  Grad[3] = dE_dalpha_y3*dEdq(1,2) + dE_dalpha_z3*dEdq(2,2);
  // dE/dbeta
  Grad[4] = dE_dbeta_x3*dEdq(0,2) + dE_dbeta_y3*dEdq(1,2) + dE_dbeta_z3*dEdq(2,2);
  // dE/dgamma
  Grad[5] = dE_dgamma_x2*dEdq(0,1) + dE_dgamma_y2*dEdq(1,1) + dE_dgamma_y3*dEdq(1,2) 
    + dE_dgamma_z3*dEdq(2,2);

  //Grad.Print("Lattice parameter gradient without MM-PBC term:\n");

  if (! Params::Parameters().NeglectManyBody()) {
    // Compute the lattice parameter gradient of the full crystal with MM.
    // We do this by finite difference.
    Vector Grad_MM_PBC = ComputeMMCrystalGradient();
    Grad += Grad_MM_PBC;
  }
  Grad.Print("Final Lattice parameter gradient:\n");


  fflush(stdout);

  // Stress tensor
  printf("Computing the stress tensor\n");
  Matrix Stress = ComputeStressTensor(dEdq);


}

// Compute MM lattice parameter gradient via finite difference of the
// periodic crystal MM energy
Vector Cluster::ComputeMMCrystalGradient() {
  
  // Initialize storage for the lattice parameters & gradient
  Vector LatticeParams(6);
  Vector LatticeGradient(6);
  LatticeGradient.Set();

  // Define the finite difference steps for abc and alpha/beta/gamma
  double delta = 0.001; // step size, in Angstroms
  double deltaTheta = 0.01; // step size in degrees

  string *params;
  params = new string[6];
  params[0] = "a"; params[1] = "b"; params[2] = "c";
  params[3] = "alpha", params[4] = "beta"; params[5] = "gamma";

  FILE *xyz;
  // xyzfile & keyfile have full paths, infile and outfile do not.
  string xyzfile, keyfile, infile, outfile;

  // Note, the Tinker coordinate file never changes during the finite
  // difference. Only the key file does.

  // Grab the current lattice parameters.  Store them in a vector as
  // a,b,c,alpha,beta,gamma.
  LatticeParams[0] = unit_cell[0].Norm(); // a
  LatticeParams[1] = unit_cell[1].Norm(); // b
  LatticeParams[2] = unit_cell[2].Norm(); // c
  
  LatticeParams[3] = RadiansToDegrees*
    acos(unit_cell[1].DotProduct( unit_cell[2] ) / (LatticeParams[1]*LatticeParams[2])); // alpha
  LatticeParams[4] = RadiansToDegrees*
    acos(unit_cell[0].DotProduct( unit_cell[2] ) / (LatticeParams[0]*LatticeParams[2])); // beta
  LatticeParams[5] = RadiansToDegrees*
    acos(unit_cell[0].DotProduct( unit_cell[1] ) / (LatticeParams[0]*LatticeParams[1])); // gamma

  // Create temp space for new lattice params we generate via finite
  // difference
  Vector NewLatticeParams(6);

  // Loop over 6 lattice parameters
  for (int i=0;i<6;i++) {

    /* Step one lattice parameter in the negative direction */

    // Set up filenames
    infile = params[i] + "_neg.xyz";
    outfile = params[i] + "_neg.out";
    xyzfile = Params::Parameters().GetMMPath() + "/" + params[i] + "_neg.xyz"; 
    keyfile = Params::Parameters().GetMMPath() + "/" + params[i] + "_neg.key"; 

    // Print the Tinker XYZ file
    if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
      printf("Monomer::ComputeMMCrystalGradient() : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    PrintTinkerCartesian(xyz);
    fclose(xyz);

    // Copy the lattice parameters
    Vector NewLatticeParams(6);
    NewLatticeParams = LatticeParams;

    if (i < 3) {
      NewLatticeParams[i] -= delta;
    }
    else {
      NewLatticeParams[i] -= deltaTheta;
    }

    // Print the key file
    FILE *key;
    if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateTinkerJob : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
    fprintf(key,"# Periodic boundary conditions\n");
    fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",NewLatticeParams[0],
	    NewLatticeParams[1],NewLatticeParams[2]);
    fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",NewLatticeParams[3],
	    NewLatticeParams[4],NewLatticeParams[5]);
    fprintf(key,"EWALD\t\tTRUE\n");
    if (!Params::Parameters().TinFoilBoundaryConditions())
      fprintf(key,"EWALD-BOUNDARY\tTRUE\n");
    fclose(key);

    // Run Tinker to compute the energy
    string job_path = Params::Parameters().GetMMPath(); 
    string cmd = "cd " + job_path;
    cmd += "; ";

    // Second command, run the job
    cmd += "analyze " + infile;
    cmd += " e > " + outfile;
    
    // Third command, switch back to base directory
    cmd += "; cd " + Params::Parameters().GetBasePath();
    system(cmd.c_str());

    // Read Energy from the output file
    double E1 = ReadTinkerEnergy(outfile); // in au

    /* Step the same lattice parameter in the positive direction */

    // Set up filenames
    infile = params[i] + "_pos.xyz";
    outfile = params[i] + "_pos.out";
    xyzfile = Params::Parameters().GetMMPath() + "/" + params[i] + "_pos.xyz"; 
    keyfile = Params::Parameters().GetMMPath() + "/" + params[i] + "_pos.key"; 

    // Print the Tinker XYZ file
    if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
      printf("Monomer::ComputeMMCrystalGradient() : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    PrintTinkerCartesian(xyz);
    fclose(xyz);

    // Copy the lattice parameters
    NewLatticeParams.Set();
    NewLatticeParams = LatticeParams;

    if (i < 3) {
      NewLatticeParams[i] += delta;
    }
    else {
      NewLatticeParams[i] += deltaTheta;
    }

    // Print the key file
    if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
      printf("Cluster::CreateTinkerJob : Cannot open file '%s'\n",
	     keyfile.c_str());
      exit(1);
    }
    fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
    fprintf(key,"# Periodic boundary conditions\n");
    fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",NewLatticeParams[0],
	    NewLatticeParams[1],NewLatticeParams[2]);
    fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",NewLatticeParams[3],
	    NewLatticeParams[4],NewLatticeParams[5]);
    fprintf(key,"EWALD\t\tTRUE\n");
    if (!Params::Parameters().TinFoilBoundaryConditions())
      fprintf(key,"EWALD-BOUNDARY\tTRUE\n");
    fclose(key);

    // Run Tinker to compute the energy
    job_path = Params::Parameters().GetMMPath(); 
    cmd = "cd " + job_path;
    cmd += "; ";

    // Second command, run the job
    cmd += "analyze " + infile;
    cmd += " e > " + outfile;
    
    // Third command, switch back to base directory
    cmd += "; cd " + Params::Parameters().GetBasePath();
    system(cmd.c_str());

    // Read Energy from the output file
    double E2 = ReadTinkerEnergy(outfile); // in au

    // Now finite difference the gradient for this parameter.
    if (i < 3) {
      LatticeGradient[i] = (E2 - E1) / (2*delta*AngToBohr);
    }
    else{
      LatticeGradient[i] = (E2 - E1) / (2*deltaTheta*DegreesToRadians);
    }
  }
  // Print the final periodic MM lattice gradient
  //LatticeGradient.Print("FDIFF PBC-MM lattice gradient\n");
  return LatticeGradient;
}


Matrix Cluster::ComputeStressTensor(Matrix LatticeGradients) {

  Matrix Stress(3,3);

  printf("WARNING:  This stress tensor routine neglects the MM pbc term\n");

  // Grab components of the unit cell vectors.  Convert to Bohr
  double a1x = unit_cell[0][0]*AngToBohr;
  double a1y = unit_cell[0][1]*AngToBohr;
  double a1z = unit_cell[0][2]*AngToBohr;

  double a2x = unit_cell[1][0]*AngToBohr;
  double a2y = unit_cell[1][1]*AngToBohr;
  double a2z = unit_cell[1][2]*AngToBohr;

  double a3x = unit_cell[2][0]*AngToBohr;
  double a3y = unit_cell[2][1]*AngToBohr;
  double a3z = unit_cell[2][2]*AngToBohr;

  printf("a1x = %f, a1y = %f, a1z = %f (bohr)\n",a1x,a1y,a1z);
  printf("a2x = %f, a2y = %f, a2z = %f (bohr)\n",a2x,a2y,a2z);
  printf("a3x = %f, a3y = %f, a3z = %f (bohr)\n",a3x,a3y,a3z);

  // Grab components of the lattice gradients - first index = xyz,
  // second = a1,a2,a3.  Units = hartress/bohr.
  double dEda1x = LatticeGradients(0,0);
  double dEda1y = LatticeGradients(1,0);
  double dEda1z = LatticeGradients(2,0);

  double dEda2x = LatticeGradients(0,1);
  double dEda2y = LatticeGradients(1,1);
  double dEda2z = LatticeGradients(2,1);

  double dEda3x = LatticeGradients(0,2);
  double dEda3y = LatticeGradients(1,2);
  double dEda3z = LatticeGradients(2,2);

  printf("dEda1x = %f, dEda1y = %f, dEda1z = %f\n",dEda1x,dEda1y,dEda1z);
  printf("dEda2x = %f, dEda2y = %f, dEda2z = %f\n",dEda2x,dEda2y,dEda2z);
  printf("dEda3x = %f, dEda3y = %f, dEda3z = %f\n",dEda3x,dEda3y,dEda3z);

  // Get Volume, convert from Ang^3 -> Bohr^3.
  double V = cell_volume * AngToBohr * AngToBohr * AngToBohr;

  // Compute components of the stress tensor
  Stress(0,0) = (dEda1x*a1x + dEda2x*a2x + dEda3x*a3x) / V;
  Stress(0,1) = (dEda1x*a1y + dEda2x*a2y + dEda3x*a3y) / V;
  Stress(0,2) = (dEda1x*a1z + dEda2x*a2z + dEda3x*a3z) / V;

  Stress(1,0) = (dEda1y*a1x + dEda2y*a2x + dEda3y*a3x) / V;
  Stress(1,1) = (dEda1y*a1y + dEda2y*a2y + dEda3y*a3y) / V;
  Stress(1,2) = (dEda1y*a1z + dEda2y*a2z + dEda3y*a3z) / V;

  Stress(2,0) = (dEda1z*a1x + dEda2z*a2x + dEda3z*a3x) / V;
  Stress(2,1) = (dEda1z*a1y + dEda2z*a2y + dEda3z*a3y) / V;
  Stress(2,2) = (dEda1z*a1z + dEda2z*a2z + dEda3z*a3z) / V;

  printf("Stress Tensor:\n");
  for (int i=0;i<3;i++) {
    printf("%15.9f  %15.9f  %15.9f\n",Stress(i,0), Stress(i,1), Stress(i,2));
  }

  return Stress;
}

Matrix Cluster::ComputeHMBIHessian() { // added on august 30th, 2010


  //For debugging purposes ONLY. Show single element of a hessian as fragments are added. Set to -1 -1 when not debugging hessian
  int ShowRow = -1;
  int ShowColumn = -1;





  //for (int i=1;i<=NDim;i++) {
  //  int indexA = Dimers[i].GetIndexA();
  //  int indexB = Dimers[i].GetIndexB();
  //  cout << "i, indexA, indexB = " << i << ", " << indexA << ", " << indexB << "\n";
  //}  
  fflush(stdout);
  
  // write the formalism here
  
  int Ntot = GetTotalNumberOfAtoms();                 
  Matrix hess(3*Ntot, 3*Ntot);
  //Matrix hess(3*UniqueAtoms, 3*UniqueAtoms);
  double scafac2=1.0;
  double scafac=1.0;
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();

/*  if (Params::Parameters().TinkerDebug() )
    scafac2 = 0.0;
  else
    scafac2 = 1.0;
*/
/*
  if (Params::Parameters().NeglectManyBody() ) {
    scafac = 0.0; // neglect many-body terms
  }
  else {
    scafac = 1.0; // treat them as usual
  }
*/
 // Contribution from QM hessian
 if (Params::Parameters().UseFullQMOnly()) {
   hess = Hess_QM;
 }
 else{
  // Contribution from full cluster MM hessian
  if (!Params::Parameters().NeglectManyBody()) {
    hess = Hess_MM;
    //hess.Set(); // partial nmb
  }

  if (Params::Parameters().NeglectManyBody()) {
    // No MM case.  Set all MM Hessians to zero.
    // This wastes memory, but it was simple.
    // Zero gradients because gradients are used for hessian calculations.
    hess.Set();
    scafac = 0.0;
    for (int i=1;i<=NMon;i++) {
      Matrix tmp_mat(3*Monomers[i].GetNumberOfAtoms(), 3*Monomers[i].GetNumberOfAtoms());
      tmp_mat.Set();
      Monomers[i].SetMMHessian(tmp_mat);

      Vector tmp_grad(3*Monomers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Monomers[i].SetMMGradient(tmp_grad);
    }
    for (int i=1;i<=NDim;i++) {
      Matrix tmp_mat(3*Dimers[i].GetNumberOfAtoms(), 3*Dimers[i].GetNumberOfAtoms());
      tmp_mat.Set();
      Dimers[i].SetMMHessian(tmp_mat);

      Vector tmp_grad(3*Dimers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Dimers[i].SetMMGradient(tmp_grad);

    }
    for (int i=1;i<=NDim_images;i++) {
      Matrix tmp_mat(3*DimerImages[i].GetNumberOfAtoms(), 3*DimerImages[i].GetNumberOfAtoms());
      tmp_mat.Set();
      DimerImages[i].SetMMHessian(tmp_mat);

      Vector tmp_grad(3*DimerImages[i].GetNumberOfAtoms());
      tmp_grad.Set();
      DimerImages[i].SetMMGradient(tmp_grad);
    }

  }

  if(Params::Parameters().BuildForceFieldOnly() ){
    // No QM case.  Set all QM Hessians to zero.
    // This wastes memory, but it was simple.
    // Zero gradients because gradients are used for hessian calculations.
    // hess.Set();
    scafac = 0.0;
    for (int i=1;i<=NMon;i++) {
      Matrix tmp_mat(3*Monomers[i].GetNumberOfAtoms(), 3*Monomers[i].GetNumberOfAtoms());
      tmp_mat.Set();
      Monomers[i].SetQMHessian(tmp_mat);

      Vector tmp_grad(3*Monomers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Monomers[i].SetQMGradient(tmp_grad);
    }
    for (int i=1;i<=NDim;i++) {
      Matrix tmp_mat(3*Dimers[i].GetNumberOfAtoms(), 3*Dimers[i].GetNumberOfAtoms());
      tmp_mat.Set();
      Dimers[i].SetQMHessian(tmp_mat);

      Vector tmp_grad(3*Dimers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Dimers[i].SetQMGradient(tmp_grad);

    }
    for (int i=1;i<=NDim_images;i++) {
      Matrix tmp_mat(3*DimerImages[i].GetNumberOfAtoms(), 3*DimerImages[i].GetNumberOfAtoms());
      tmp_mat.Set();
      DimerImages[i].SetQMHessian(tmp_mat);

      Vector tmp_grad(3*DimerImages[i].GetNumberOfAtoms());
      tmp_grad.Set();
      DimerImages[i].SetQMGradient(tmp_grad);
    }
  }

  //Case where only the full MM contributes to the Hessian
  // Zero gradients because gradients are used for hessian calculations.
  if(Params::Parameters().UseFullMMOnly()){
    for (int i=1;i<=NMon;i++) {
      Matrix tmp_mat(3*Monomers[i].GetNumberOfAtoms(), 3*Monomers[i].GetNumberOfAtoms());
      tmp_mat.Set();
      Monomers[i].SetQMHessian(tmp_mat);
      Monomers[i].SetMMHessian(tmp_mat);

      Vector tmp_grad(3*Monomers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Monomers[i].SetQMGradient(tmp_grad);
      Monomers[i].SetMMGradient(tmp_grad);
    }
    for (int i=1;i<=NDim;i++) {
      Matrix tmp(3*Dimers[i].GetNumberOfAtoms(), 3*Dimers[i].GetNumberOfAtoms());
      tmp.Set();
      Dimers[i].SetQMHessian(tmp);
      Dimers[i].SetMMHessian(tmp);

      Vector tmp_grad(3*Dimers[i].GetNumberOfAtoms());
      tmp_grad.Set();
      Dimers[i].SetQMGradient(tmp_grad);
      Dimers[i].SetMMGradient(tmp_grad);
    }
    for (int i=1;i<=NDim_images;i++) {
      Matrix tmp(3*DimerImages[i].GetNumberOfAtoms(), 3*DimerImages[i].GetNumberOfAtoms());
      tmp.Set();
      DimerImages[i].SetQMHessian(tmp);
      DimerImages[i].SetMMHessian(tmp);
      
      Vector tmp_grad(3*DimerImages[i].GetNumberOfAtoms());
      tmp_grad.Set();
      DimerImages[i].SetQMGradient(tmp_grad);
      DimerImages[i].SetMMGradient(tmp_grad);
    }
  }

  /*
  if (!(Params::Parameters().TinkerDebug()) )                                      
    hess = Hess_MM;
  else {                                             
    // No MM case.  Set all MM gradients to zero.
    // This wastes memory, but it was simple.
    hess.Set();
  }

  if (Params::Parameters().Do_fdTinkerHessian() ) {
    Vector origvec = GetCurrentCoordinates();
    hess = GetFiniteDifferenceTinkerHessian( origvec );
    PrintHessian("FD Tinker Hessian",hess);
  }
  */
  
  //PrintHessian("Hessian initially",hess);


  hess.Scale(scafac);

  // Create a list of the starting index for each monomer in the
  // full hessian.  Like the Monomers, The key indexes count from 1->NMon.
  // On the other hand, in the hessian, the indexing starts at (0,0).
  //int *row_key = new int[NMon+1];
  int row_key[NMon+1];
  int key[NMon+1];
  //int *key = new int[NMon+1];
  row_key[1] = 0;
  key[1] = 0; // set the first one by hand
  for (int i=2;i<=NMon;i++) {
    row_key[i] = row_key[i-1] + 3*Monomers[i-1].GetNumberOfAtoms();
    key[i] = row_key[i]*3*Ntot + row_key[i];
  }


  //list of atoms in the hessian
  //int atom_key[UniqueAtoms];
  //element to start 
  //int start_key[UniqueMon];
  //start_key[0] = 0;
  //int key_index = 0;
  //for(int i=1;i<=NMon;i++)
  //  if(Monomers[i].GetSymmetryFactor() != 0)
  //    for(int j=0;j<Monomers[i].GetNumberOfAtoms();j++){
  //	atom_key[key_index] = Monomers[i].GetAtom(j).GetGlobalIndex();
	//start_key[key_index] = start_key[key_index-1] + 3;
  //	key_index++;
  //    }

  //print atom_key
  //  printf("atom_key = ");
  //  for(int i=0;i<UniqueAtoms;i++)
  //    printf("%i ",atom_key[i]);
  //  printf("]\n");

  // 1-body contributions
  for (int i=1;i<=NMon;i++) {
    int Na = Monomers[i].GetNumberOfAtoms();
    int start = key[i];
    for (int j=0;j<3*Na;j++) {
      for (int k=0; k<3*Na;k++) {
	if ( !Params::Parameters().TinkerDebug() )
	  hess.Element(row_key[i]+j,row_key[i]+k) += scafac2*Monomers[i].GetQMHessian().Element(j,k);
	//if (Params::Parameters().GetMMType() != 2) // if not AIFF
	hess.Element(row_key[i]+j,row_key[i]+k) += -scafac*Monomers[i].GetMMHessian().Element(j,k);
	

	if( (row_key[i]+j == ShowRow) && (row_key[i]+k == ShowColumn) ){
	  printf("m%i\n",Monomers[i].GetIndex());
	  printf("QM(%i,%i) =  %f\n",row_key[i]+j,row_key[i]+k,
	   Monomers[i].GetQMHessian().Element(j,k));
	printf("MM(%i,%i) =  %f\n",row_key[i]+j,row_key[i]+k,
	   Monomers[i].GetMMHessian().Element(j,k));
	printf(" hess = %15.9f\n", hess.Element(row_key[i]+j,row_key[i]+k) );
	}
	
      }
    }
  }

  //PrintHessian("HMBI Hessian after 1-body contributation",hess);

  // 2-body interaction contributions
  
  /*
    The number of dimers calulated were reduced under symmetry.
    Each dimer had a list of symmetrically equivalent dimers.
    The full list of dimers are  being recreated and elements of the caluated dimer
    hessian are rotated to find the hessian of the equivalent dimer
  */
  
  //printf("Hess(3,3) = %f\n",hess(3,3));

  //dimer 12 hessian create from eqvalent dimer 34.
  //H12=R13*H34*R43' where R13 and R43 are rotation matrices.
  for (int i=1;i<=NDim;i++) {           
    int Na = Dimers[i].GetMonomerA().GetNumberOfAtoms();
    int Nb = Dimers[i].GetMonomerB().GetNumberOfAtoms();

    //Finding rotation matrix to map monomer A and B to there symmetrical unique monomer
    //Not the full rotation but used to make it.
    int index3 = Dimers[i].GetIndexA();        
    int index4 = Dimers[i].GetIndexB();        
    //Matrix Rot3= Monomers[index3].GetRotationMatrix();
    //Matrix Rot4= Monomers[index4].GetRotationMatrix();

    //int startA = key[indexA];
    //int startB = key[indexB];
    double damp;
    double c1=0, c0=0;
    if ( Params::Parameters().DoLocal2BodyTruncation() || Params::Parameters().IsPeriodic()  ) {
      c0 = Params::Parameters().GetLocalCutoff(0);  
      c1 = Params::Parameters().GetLocalCutoff(1);
      damp = Dimers[i].GetDampingFactor(c0,c1);
      //cout << "daya damp\n";
    }
    else {
      damp = 1.0;
    }
 
    fflush(stdout);
    
    //cout << "i = " << i<< "\n";
    //cout << "separation = " << Dimers[i].GetMonomerA().FindDistance( Dimers[i].GetMonomerB() ) << "\n";
    //cout << "damp = " << damp << "\n";
    //cout<< " index3, index4 = " << index3 << ", " << index4 << "\n";
    //Dimers[i].GetSpatialDampingFunctionHessian(c0,c1).PrintHessian("damping hessian");
    fflush(stdout);
    if(Dimers[i].GetSymmetryFactor() && damp >= pow(10.0,-6) ){
      
      //Getting list of all dimers that are symmetrical to this dimer. Rotating the
      //elements and adding them to the right entry.
      //Each entry in the list repressent a monomer so two entries are needed for a dimer.
      vector<int> SymmetryList = Dimers[i].GetSymmetryList();
      vector<int> PeriodicSymmetryList = Dimers[i].GetPeriodicSymmetryList();
      int NumberOfSymmetricalDimers = (SymmetryList.size() + PeriodicSymmetryList.size())/2;
      for(int j=0; j<NumberOfSymmetricalDimers;j++){
	//index of the monomers in the symmetrical equivalent dimer
	int index1;
	int index2;
        Matrix Rot;
	
        //Index of symmetrical atom on symmetry List
        Vector SymAtomList;
	
	//The symmetrical equivalent dimers that are periodic images have their contribution to the hessian scaled by 1/2
	float per_sca;
	
	//Symmetrical equivalent dimer is a not periodic image
	if(j<SymmetryList.size()/2){
	  //2*j is the first molecule in the dimer. 2*j+1 is the second monomer in the dimer
	  //printf("j=%i\n",j);
	  index1 = SymmetryList[2*j];
	  index2 = SymmetryList[2*j+1];
	  Rot = Dimers[i].GetRotationList()[j];
	  SymAtomList = Dimers[i].GetAtomEquivalency()[j];
	  per_sca = 1.0;
	}
	//Symmetrical equivalent dimer is a periodic image
	else{
	  //2*k is the first molecule in the dimer. 2*k+1 is the second monomer in the dimer
	  int k = j - SymmetryList.size()/2;
	  //printf("k=%i j=%i\n",k,j);
	  index1 = PeriodicSymmetryList[2*k];
	  index2 = PeriodicSymmetryList[2*k+1];
	  Rot = Dimers[i].GetPeriodicRotationList()[k];
	  SymAtomList = Dimers[i].GetPeriodicAtomEquivalency()[k];
	  per_sca = 0.5;
	}
	
	Rot.Transpose();
	//Rot.Print("Rot");

	//SymAtomList.Print("SymAtomList");
	//Finding rotation matrix to map monomer A and B to there symmetrical unique monomer
	//Not the full rotation but used to make it.
	//int startA = key[index1];
	//int startB = key[index2];

	//printf("d(%i,%i) mapping from d(%i,%i)\n",
	//       index1,index2,index3,index4);
	//fflush(stdout);
	//Matrix Rot1 = Monomers[index1].GetRotationMatrix();
	//Matrix Rot2 = Monomers[index2].GetRotationMatrix();
		
	//Piece arising from first monomer. 
	//Rotation under symmetry H12 = R13*H34*R24' where R13=R24 and R13=R1'*R3
	//Matrix Rot13 = Rot1.Multiply(Rot3,2);
	//Matrix Rot42 = Rot3.Multiply(Rot1,2);
	//printf("piece arising from first monomer\n");
	//Rot13.Print("Rot13");
	//Rot42.Print("Rot42");


	//Determines if SymA is symmetrical to MonA or MonB by determining which monomer is symmetrical atom 1 in dimer.
	//Need to know if MonA and MonB have different number to atoms
	bool FirstAInDim = false;
	for(int atom = 0; atom < Monomers[index3].GetNumberOfAtoms();atom++){
	  if(SymAtomList[atom] == 1) FirstAInDim = true;
	  //printf("atom = %i maxatom = %i index = %i SymAtomList = %i \n",atom,Monomers[index3].GetNumberOfAtoms(),
	  //   index3,(int)SymAtomList[atom]);
	}
          
        //SymAtomList.Print("SymAtomList");
	//printf("FirstAInDim = %i\n",FirstAInDim);

	for(int iatom=0;iatom<Na+Nb;iatom++){
	  
          int MonA;
	  int Symi = (int) SymAtomList[iatom] - 1;
          
          //Determines if monA and monB are the same monomer. If they are, the monomer hessian contributes to two body hessian.
          bool same_monomer;
          bool FirstA = 0;
  
	  //Determining atom index on monomer
	  if(FirstAInDim){
	    if(Symi < Na ){
	      MonA = index1;
	      FirstA = 1;
	    }else{
	      MonA = index2;
	      Symi -= Na;
	    }
	  }else{
	    if(Symi < Nb ){
	      MonA = index1;
	      FirstA = 1;
	    }
	    else{
	      MonA = index2;
	      Symi -= Nb;
	    }
	  }
	  
	  for(int jatom=0;jatom<Na+Nb;jatom++){
	    
	    
	    int MonB;
	    int Symj = (int) SymAtomList[jatom] - 1;
	    if(FirstAInDim){
	      if(Symj < Na ){
		MonB = index1;
		if(FirstA) same_monomer = 1;
		else same_monomer = 0;
	      }else{
		MonB = index2;
		Symj -= Na;
		if(!FirstA) same_monomer = 1;
		else same_monomer = 0;
	      }
	    }else{
	      if(Symj < Nb ){
		MonB = index1;
		if(FirstA) same_monomer = 1;
		else same_monomer = 0;
	      }
	      else{
		MonB = index2;
		Symj -= Nb;
		if(!FirstA) same_monomer = 1;
		else same_monomer = 0;
	      }
	    }
	    
	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //including dimer contributation
		    //printf("element (%i %i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //printf("add to element (%i,%i)\n",row_key[index1]+3*iatom+xyz1,row_key[index1]+3*jatom+xyz2);
		    //fflush(stdout);
		    
		    hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += per_sca*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*Dimers[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4) - scafac*Dimers[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		    
		    if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		      printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		      printf("dimer hess");
                      printf("j= %i\n",j);
                      printf("element = %i, %i\n",row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2);
                      Rot.Print("Rot");
		      printf(" damp = %f\n",damp);
		      printf(" per_sca = %f\n",per_sca);
		      printf(" iatom = %i jatom = %i\n",iatom,jatom);
                      printf(" xyz3 = %i xyz4 = %i\n",xyz3,xyz4);
                      printf(" Symi = %i Symj = %i\n",Symi,Symj);
		      printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		      printf(" Unrotated QM = %15.9f\n", Dimers[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" Unrotated MM = %15.9f\n", Dimers[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
                      printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		      printf(" rotated QM= %15.9f\n", Rot(xyz1,xyz3)*Rot(xyz2,xyz4)*Dimers[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" rotated MM= %15.9f\n", Rot(xyz1,xyz3)*Rot(xyz2,xyz4)*Dimers[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) );
		    }
		    
 
		    //damping factor contributation to the Hessian.
		    //Note that gradients are already partially rotated so the Rot1 matrix is used instend of Rot13 or Rot42
		    if(Params::Parameters().DoLocal2BodyTruncation()){
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += per_sca
			*Rot(xyz2,xyz4)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
			*Rot(xyz1,xyz3)*(scafac2*Dimers[i].GetQMGradient()[3*iatom+xyz3]-scafac*Dimers[i].GetMMGradient()[3*iatom+xyz3]);

		      
		    if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		        printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
                        printf("grad 1\n");
                        Rot.Print("Rot");
		        printf(" rot1(xyz4,xyz2) = %15.9f\n", Rot(xyz2,xyz4));
		        printf(" rot1(xyz3,xyz1) = %15.9f\n", Rot(xyz1,xyz3));
		        printf(" Damp Grad 1 = %15.9f\n", Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
			printf(" Dimer Grad QM 1 = %15.9f\n", Dimers[i].GetQMGradient()[3*iatom+xyz3]);
		        printf(" Grad QM 1 product = %15.9f\n",  Rot(xyz2,xyz4)*Rot(xyz1,xyz3)*
			       Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]*
		               Dimers[i].GetQMGradient()[3*iatom+xyz3]);
		        printf(" Dimer Grad MM 1 = %15.9f\n", Dimers[i].GetMMGradient()[3*iatom+xyz3]);
		        printf(" Grad MM 1 product = %15.9f\n", Rot(xyz2,xyz4)*Rot(xyz1,xyz3)*
		               Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]*
		               Dimers[i].GetMMGradient()[3*iatom+xyz3]);
		        printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[index1]+3*iatom+xyz1, row_key[index1]+3*jatom+xyz2));
		      }     
		      

		       hess.Element(row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2) += per_sca
			 *Rot(xyz1,xyz3)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
			 *Rot(xyz2,xyz4)*(scafac2*Dimers[i].GetQMGradient()[3*jatom+xyz4]-scafac*Dimers[i].GetMMGradient()[3*jatom+xyz4]);
		      
		       
		       if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		         printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
                         printf("grad 2\n");  
                         Rot.Print("Rot");                   
			 printf(" Damp Grad 2 =%15.9f\n", Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
			 printf(" Dimer Grad QM 2 = %15.9f\n", Dimers[i].GetQMGradient()[3*jatom+xyz4]);
			 printf(" Grad QM 2 product test= %15.9f\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)*
				Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
				*Dimers[i].GetQMGradient()[3*jatom+xyz4]);
			 printf(" Dimer MM 2 = %15.9f\n", Dimers[i].GetMMGradient()[3*jatom+xyz4]);
			 printf(" Grad MM 2 product = %15.9f\n", Rot(xyz2,xyz4)*Rot(xyz1,xyz3)       
				*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
				*Dimers[i].GetMMGradient()[3*jatom+xyz4]);
			 printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[index1]+3*iatom+xyz1, row_key[index1]+3*jatom+xyz2));
		      }
		       

		      hess.Element(row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2) += per_sca*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			*Dimers[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
			*(scafac2*Dimers[i].GetQMIntEnergy()-scafac*Dimers[i].GetMMIntEnergy());

		      
		      if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		        printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
                        printf("Damp Hess\n");
                        Rot.Print("Rot");
                        printf(" Damp Hess = %15.9f\n",Dimers[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4));
		        printf(" QM interaction = %15.9f MM interaction = %15.9f \n",Dimers[i].GetQMIntEnergy(),Dimers[i].GetMMIntEnergy());
		        printf(" Damp Hess product= %15.9f\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		               *Dimers[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
		               *(Dimers[i].GetQMIntEnergy()-Dimers[i].GetMMIntEnergy()));
		        printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[index1]+3*iatom+xyz1, row_key[index1]+3*jatom+xyz2));
		      }
		      
		      
		    }
		      
		    
		  }//loop over xyz4

	          //Rotating the damping grad but not the monomer grad
		  if(Params::Parameters().DoLocal2BodyTruncation()){

		    //QM Monomers:already subtracted if counterpoise corrected
		    if(!Params::Parameters().DoCounterpoise())
		      hess.Element(row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2) -= per_sca*scafac2*
			(Rot(xyz2,xyz3)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*Monomers[MonA].GetQMGradient()[3*Symi+xyz1]
			 +Rot(xyz1,xyz3)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*Monomers[MonB].GetQMGradient()[3*Symj+xyz2]);  
		    
		    if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3) > 0.0001 ))){
		      printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		      printf("Monomer Gradient Damp Gradient Term\n");
		      printf(" xyz1 = %i xyz2 = %i xyz3 = %i\n",xyz1,xyz2,xyz3);
		      printf(" rot(xyz2,xyz3) = %15.9f\n", Rot(xyz2,xyz3));
		      printf(" Monomer QM 1 = %15.9f\n", Monomers[MonA].GetQMGradient()[3*Symi+xyz1]);
		      printf(" Monomer QM 1 product = %15.9f\n", Rot(xyz2,xyz3)*Monomers[MonA].GetQMGradient()[3*Symi+xyz1]*
		             Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		      printf(" rot(xyz1,xyz3) = %15.9f\n", Rot(xyz1,xyz3));
		      printf(" Monomer QM 2 = %15.9f\n", Monomers[MonB].GetQMGradient()[3*Symj+xyz2]);
		      printf(" Monomer QM 2 product = %15.9f\n", Rot(xyz1,xyz3)*Monomers[MonB].GetQMGradient()[3*Symj+xyz2]*
		             Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2));
		    }
		    
		    
		    //MM Monomers:already subtrated if qchem is used and counterpoise corrected
		    if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		      hess.Element(row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2) += per_sca*scafac*
			(Rot(xyz2,xyz3)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*Monomers[MonA].GetMMGradient()[3*Symi+xyz1]
			 +Rot(xyz1,xyz3)*Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*Monomers[MonB].GetMMGradient()[3*Symj+xyz2]);	

		    
		    if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3) > 0.0001 ))){
 		      printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		      printf("Monomer Gradient\n");
		      printf(" xyz1 = %i xyz2 = %i xyz3 = %i\n",xyz1,xyz2,xyz3);
		      printf(" rot(xyz2,xyz3) = %15.9f\n", Rot(xyz2,xyz3));
		      printf(" Monomer MM 1 = %15.9f\n", Monomers[MonA].GetMMGradient()[3*Symi+xyz1]);
		      printf(" Monomer MM 1 product = %15.9f\n", Rot(xyz2,xyz3)*Monomers[MonA].GetMMGradient()[3*Symi+xyz1]*
			     Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*Symj+xyz3]);
		      printf(" rot(xyz1,xyz3) = %15.9f\n", Rot(xyz1,xyz3));
		      printf(" Momomer MM 2 = %15.9f\n", Monomers[MonB].GetMMGradient()[3*Symj+xyz2]);
		      printf(" Monomer MM 2 product = %15.9f\n", Rot(xyz3,xyz1)*Monomers[MonB].GetMMGradient()[3*Symj+xyz2]*
			     Dimers[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1,row_key[MonB]+3*Symj+xyz2));
		    }
		    
		    
		  }

		}//loop over xyz3
		
		 // Monomer hessian contribution is zero if iatom and jatom are not on the same monomer
		
		if(same_monomer){
		 
		  if(!Params::Parameters().DoCounterpoise())
		    //QM Monomers
		    hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) -= per_sca*damp
		      *scafac2*Monomers[MonA].GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		  
                  if( (row_key[MonA]+3*Symi+xyz1 == ShowRow) && (row_key[MonB]+3*Symj+xyz2 == ShowColumn)){
 		      printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		       printf("Monomer QM Hess = %15.9f\n",Monomers[MonA].GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		  }
		  
		  
		  //MM Monomers
		  if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		    hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += per_sca*damp*scafac*
		      Monomers[MonA].GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		  
		  if( (row_key[MonA]+3*Symi+xyz1 == ShowRow) && (row_key[MonB]+3*Symj+xyz2 == ShowColumn) ){
		    printf("using d(%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		   printf(" Monomer MM Hess = %15.9f\n",Monomers[index1].GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		   printf(" Current Hess Sum = %15.9f\n\n\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) );
                 }
		  

		}//end of if statement for MonA == MonB
	      }//loop over xyz2
	    }//loop over xyz1
	    
	  }//loop over jatom 
	}//loop over iatom
	
	
      }//loop over entries in Symmetry List
    }//if statement for the symmetry factor
    
  }//loop over dimers
  
  //PrintHessian("HMBI Hessian before PBC contribution",hess);  
  
  // 2-body contributions due to periodic boundary conditions
  if (Params::Parameters().IsPeriodic() ) { 

    //dimer 12 hessian create from eqvalent dimer 34.
    //H12=R13*H34*R43' where R13 and R43 are rotation matrices.
    for (int i=1;i<=NDim_images;i++) {    
      
      int Na = DimerImages[i].GetMonomerA().GetNumberOfAtoms();
      int Nb = DimerImages[i].GetMonomerB().GetNumberOfAtoms();       
      
      //Finding rotation matrix to map monomer A and B to there symmetrical unique monomer
      //Not the full rotation but used to make it.   
      int index3 = DimerImages[i].GetIndexA();        
      int index4 = DimerImages[i].GetReferenceMonomerIndex();   
      //Matrix Rot3= Monomers[index3].GetRotationMatrix();
      //Matrix Rot4= Monomers[index4].GetRotationMatrix();
      
      //damping functation
      double c1=0, c0=0;
      c0 = Params::Parameters().GetLocalCutoff(0);
      c1 = Params::Parameters().GetLocalCutoff(1);
      double  damp = DimerImages[i].GetDampingFactor(c0,c1);
      
      //Getting list of all dimers that are symmetrical to this dimer. Rotating the
      //elements and adding them to the right entry.
      //Each entry in the list repressent a monomer so two entries are needed for a dimer.
      vector<int> SymmetryList = DimerImages[i].GetSymmetryList();
      if(damp >= pow(10.0,-6)){
	for(int j=0; j<SymmetryList.size()/2;j++){
	  
	  //2*j is the first molecule in the dimer. 2*j+1 is the second monomer in the dimer
	  int index1 = SymmetryList[2*j];
	  int index2 = SymmetryList[2*j+1];
	  
	  Matrix Rot = DimerImages[i].GetRotationList()[j];
	  Rot.Transpose();
	  Vector SymAtomList = DimerImages[i].GetAtomEquivalency()[j];
	  
	  //printf("j= %i index1 = %i index2 = %i\n",j,index1,index2);
	  //fflush(stdout);
	  //int startA = key[index1];
	  //int startB = key[index2];
	  
	  //Matrix Rot1 = Monomers[index1].GetRotationMatrix();
	  //Matrix Rot2 = Monomers[index2].GetRotationMatrix();
	  
	  //printf("d(%i,%i) mapping from d(%i,%i)\n",
	  //       index1,index2,index3,DimerImages[i].GetIndexB());
	  //printf("row_key[index1] =%i row_key[index2]= %i\n",row_key[index1],row_key[index2]);
	  
	  //Rotation under symmetry H12 = R13*H34*R24' where R13=R24 and R13=R1'*R3
	  
	  
	  //Matrix Rot13 = Rot1.Multiply(Rot3,2);
	  //Matrix Rot42 = Rot3.Multiply(Rot1,2);
	  //printf("piece arising from first monomer\n");
	  //Rot13.Print("Rot13");
	  //Rot42.Print("Rot42");

	//Determines if SymA is symmetrical to MonA or MonB by determining which monomer is symmetrical atom 1 in dimer.
	//Need to know if MonA and MonB have different number to atoms
	  bool FirstAInDim = false;
	  for(int atom = 0; atom < Monomers[index3].GetNumberOfAtoms();atom++){
	    if(SymAtomList[atom] == 1) FirstAInDim = true;
	  }

          
          //SymAtomList.Print("SymAtomList");
	  //printf("FirstAInDim = %i\n",FirstAInDim);

	  for(int iatom=0;iatom<Na+Nb;iatom++){
	    
	    
            //Determines if monA and monB are the same monomer. If they are, the monomer hessian contributes to two body hessian.
            bool same_monomer;
            bool FirstA = 0;

	    //Determining atom index on monomer
	    int MonA;
	    int Symi = (int) SymAtomList[iatom] - 1;
	    if(FirstAInDim){
	      if(Symi < Na ){
		MonA = index1;
		FirstA = 1;
	      }else{
		MonA = index2;
		Symi -= Na;
	      }
	    }else{
	      if(Symi < Nb ){
		MonA = index1;
		FirstA = 1;
	      }
	      else{
		MonA = index2;
		Symi -= Nb;
	      }
	    }
	    
	    for(int jatom=0;jatom<Na+Nb;jatom++){
        
	      int MonB;
	      int Symj = (int) SymAtomList[jatom] - 1;
	      if(FirstAInDim){
		if(Symj < Na ){
		  MonB = index1;
		  if(FirstA) same_monomer = 1;
		  else same_monomer = 0;
		}else{
		  MonB = index2;
		  Symj -= Na;
		  if(!FirstA) same_monomer = 1;
		  else same_monomer = 0;
		}
	      }else{
		if(Symj < Nb ){
		  MonB = index1;
		  if(FirstA) same_monomer = 1;
		  else same_monomer = 0;
		}
		else{
		  MonB = index2;
		  Symj -= Nb;
		  if(!FirstA) same_monomer = 1;
		  else same_monomer = 0;
		}
	      }      	      
	      
	      for(int xyz1=0;xyz1<3;xyz1++){
		for(int xyz2=0;xyz2<3;xyz2++){
		  //xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		  for(int xyz3=0;xyz3<3;xyz3++){
		    for(int xyz4=0;xyz4<3;xyz4++){
		      //including dimer contributation
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			*(scafac2*DimerImages[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4) 
			  -scafac*DimerImages[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));

		      if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
			printf("using d(%i,%i)\n",DimerImages[i].GetIndexA(),DimerImages[i].GetIndexB());
			printf("FirstA = %i\n",FirstA);
			printf("same_monomer = %i\n",same_monomer);
			printf("MonA = %i MonB = %i\n",MonA,MonB);
                        Rot.Print("Rot");
			printf(" damp = %15.9f\n",damp);
			printf("iatom = %i jatom = %i\n",iatom,jatom);
			printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
			printf(" Unrotated QM = %15.9f\n", DimerImages[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
			printf(" Unratated MM = %15.9f\n", DimerImages[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
			printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
			printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
			printf(" rotated QM = %f\n", damp*DimerImages[i].GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4)*Rot(xyz1,xyz3)*Rot(xyz2,xyz4));
			printf(" rotated MM = %f\n", damp*DimerImages[i].GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4)*Rot(xyz1,xyz3)*Rot(xyz2,xyz4));
			printf(" Symi = %i Symj = %i\n",Symi,Symj);
			printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		      }
		      

		      //damping factor contributation to the Hessian.
		      //Note that gradients are already partially rotated so the Rot1 matrix is used instend of Rot13 or Rot42
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5
			*Rot(xyz1,xyz3)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
			*Rot(xyz2,xyz4)*(scafac2*DimerImages[i].GetQMGradient()[3*jatom+xyz4]-scafac*DimerImages[i].GetMMGradient()[3*jatom+xyz4]);
		     
		      
		      if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
			printf(" Rot(xyz1,xyz3) = %15.9f\n",Rot(xyz1,xyz3));
			printf(" Rot(xyz2,xyz4) = %15.9f\n",Rot(xyz2,xyz4));
			printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
			printf(" Damp Grad 1 =  %15.9f\n",DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
			printf(" Dimer QM Grad 1 = %15.9f\n",DimerImages[i].GetQMGradient()[3*jatom+xyz4]);
	                printf(" Dimer QM Grad 1 product = %15.9f\n",  Rot(xyz2,xyz4)*Rot(xyz1,xyz3)*
			       DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
			       *DimerImages[i].GetQMGradient()[3*jatom+xyz4]);
			printf(" Dimer MM Grad 1 = %15.9f\n",DimerImages[i].GetMMGradient()[3*jatom+xyz4]);
			printf(" Dimer QM Grad 1 product = %15.9f\n",  Rot(xyz2,xyz4)*Rot(xyz1,xyz3)*
			       DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
			       *DimerImages[i].GetMMGradient()[3*jatom+xyz4]);
			printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		      }
		      
		      
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5
			*Rot(xyz2,xyz4)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
			*Rot(xyz1,xyz3)*(scafac2*DimerImages[i].GetQMGradient()[3*iatom+xyz3]-scafac*DimerImages[i].GetMMGradient()[3*iatom+xyz3]);
		      
		      
		      if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		        printf(" Damp Grad 2 = %15.9f\n", DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
		        printf(" Dimer Grad 2 = %15.9f\n",DimerImages[i].GetQMGradient()[3*iatom+xyz3]);
			printf(" Dimer QM Grad 2 Product = %15.9f\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			       *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
			       *DimerImages[i].GetQMGradient()[3*iatom+xyz3]);
			printf(" Dimer MM Grad 2 = %15.9f\n",DimerImages[i].GetMMGradient()[3*jatom+xyz4]);
			printf(" Dimer QM Grad 2 Product = %15.9f\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			       *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
			       *DimerImages[i].GetMMGradient()[3*iatom+xyz3]);
			printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) );
		      }
		      
		      
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			*DimerImages[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
			*(scafac2*DimerImages[i].GetQMIntEnergy()-scafac*DimerImages[i].GetMMIntEnergy());
		
		      
		      if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2) == ShowColumn && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		        printf(" 3*iatom+xyz3 = %i\n", 3*iatom+xyz3);
		        printf(" 3*jatom+xyz4 = %i\n", 3*jatom+xyz4);
		        printf(" QM interaction = %15.9f MM interaction = %15.9f \n",DimerImages[i].GetQMIntEnergy(),DimerImages[i].GetMMIntEnergy());
		        printf(" Hess = %15.9f \n",DimerImages[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4));
		        printf(" QM Hess product = %15.9f product\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			       *DimerImages[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
			       *DimerImages[i].GetQMIntEnergy());
			printf(" MM Hess product = %15.9f product\n",Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
			       *DimerImages[i].GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
			       *DimerImages[i].GetMMIntEnergy());
			printf("Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		      }
		      
		      
		    }//loop over xyz4
		    
		    //Rotating gradient monomers contribation
		    
		    //QM:already subtracted if  counterpoise corrected
		    if(!Params::Parameters().DoCounterpoise())
	              hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) -= 0.5*scafac2*
			(Rot(xyz2,xyz3)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*Monomers[MonA].GetQMGradient()[3*Symi+xyz1]
			 +Rot(xyz1,xyz3)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*Monomers[MonB].GetQMGradient()[3*Symj+xyz2]);
		    
		    
		    if((row_key[MonA]+3*Symi+xyz1) == ShowRow && (row_key[MonB]+3*Symj+xyz2 == ShowColumn) && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		      printf(" Monomer QM 1 = %15.9f\n", Monomers[MonA].GetQMGradient()[3*Symi+xyz1]);
		      printf(" Monomer QM 1 product = %15.9f\n", Rot(xyz2,xyz3)
			     *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]
			     *Monomers[MonA].GetQMGradient()[3*Symi+xyz1]);
		      printf(" Monomer QM 2 = %15.9f\n", Monomers[MonB].GetQMGradient()[3*Symj+xyz2]);
		      printf(" Monomer QM 2 product = %15.9f\n", Rot(xyz1,xyz3)
			     *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
			     *Monomers[MonB].GetQMGradient()[3*Symj+xyz2]);
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		    }
		    
		    
		      //MM:already subtracted if qchem is used and counterpoise corrected
		    if(!Params::Parameters().DoCounterpoise()|| !Params::Parameters().GetMMType()!=3)
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5*scafac*
			(Rot(xyz2,xyz3)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*Monomers[MonA].GetMMGradient()[3*Symi+xyz1]
			 +Rot(xyz1,xyz3)*DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*Monomers[MonB].GetMMGradient()[3*Symj+xyz2]);
		    
		    
		    if((row_key[MonA]+3*Symi+xyz1 == ShowRow) && (row_key[MonB]+3*Symj+xyz2 == ShowColumn) && (fabs(Rot(xyz1,xyz3)) > 0 || fabs(Rot(xyz2,xyz3))) ){	
		      printf(" Monomer MM 1 = %15.9f\n", Monomers[MonA].GetMMGradient()[3*Symi+xyz1]);
		      printf(" Monomer MM 1 product = %15.9f\n", Rot(xyz3,xyz1)*Monomers[MonA].GetMMGradient()[3*Symi+xyz1]
		   	     *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		      printf(" Monomer MM 2 = %15.9f\n", Monomers[MonB].GetMMGradient()[3*Symj+xyz2]);
		      printf(" Monomer MM 2 product = %15.9f\n", Rot(xyz1,xyz3)*Monomers[MonB].GetMMGradient()[3*Symj+xyz2]
		  	     *DimerImages[i].GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		    }
		    

		  }//loop over xyz3
		  
		  
		  // Monomer hessian contribution is zero if iatom and jatom are not on the same monomer
		  
		  if(same_monomer){
		    //QM Monomers
		    if(!Params::Parameters().DoCounterpoise())
		      hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) -= 0.5*damp*scafac2
		        *Monomers[MonA].GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		    
		    //MM Monomers   
		    if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3)
		      hess.Element(  row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2) += 0.5*damp*scafac
		        *Monomers[MonA].GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		    
		    if( (row_key[MonA]+3*Symi+xyz1 == ShowRow) && (row_key[MonB]+3*Symj+xyz2 == ShowRow) ){
		      printf(" Monomer Hess QM = %15.9f\n",Monomers[MonA].GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		      printf(" Monomer Hess MM = %15.9f\n",Monomers[MonA].GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		      printf("Current Hess Sum = %15.9f\n\n", hess.Element( row_key[MonA]+3*Symi+xyz1, row_key[MonB]+3*Symj+xyz2));
		    }
		    
		  }
		  
		}//loop over xyz2
	      }//loop over xyz1
	    }//loop over jatom
	  }//loop over iatom

	}//loop over entries of symmetry list
      }//if statement checking if damp greater or equal to zero
    }//loop over image dimers
  }//if statement for periodicity
 }//if statement for full QM only

  //PrintHessian("HMBIHessian",hess);
  //delete [] key;

  
  return hess;
}

void Cluster::PrintHessian(string title, Matrix hess) {

  int dim = hess.GetRows();
  cout << "dim = " << dim <<"\n";
  printf("%s\n",title.c_str());
  for (int i=0;i<(dim-dim%9)/9;i++) {
    for (int j1=0;j1<1;j1++) {
      cout << "  \t";
      for (int k1=9*i;k1<9*i+9;k1++) {
        cout <<  setprecision(6) << showpoint << fixed << right << setw(15)<< k1 << "\t";
      } 
    cout << "\n";
    }
    for (int j=0;j<dim;j++) {
      cout << j << "\t";
      for (int k=9*i;k<9*i+9;k++) {
        cout <<  setprecision(6) << showpoint << fixed << right << setw(15)<< hess.Element(j,k) << "\t";
      }
    cout << "\n";
    }
  }
  if (dim%9 != 0) {
    for (int j1=0;j1<1;j1++) {
      cout << "  \t";
      for (int k1=(dim-dim%9);k1<dim;k1++) {
        cout <<  setprecision(6) << showpoint << fixed << right << setw(15)<< k1 << "\t";
      }             
    cout << "\n";
    }  
    for (int j=0;j<dim;j++) {
      cout << j << "\t";
      for (int k=(dim-dim%9);k<dim;k++) {
        cout <<  setprecision(6) << showpoint << fixed << right << setw(15)<< hess.Element(j,k) << "\t";                     
      }
    cout << "\n";
    }
  }

}

void Cluster::SetMMHessian() {
  
  if (Params::Parameters().GetMMType() == 1){  // Tinker


    //If Tinker does not include Dr. Berans correction to the
    //polarization hessian. Find hessian through finite difference
    if(Params::Parameters().Do_fdTinkerHessian()){
    //if(Params::Parameters().IsPeriodic())
      if(Params::Parameters().RunJobs()){
	printf("Finding PBC Hessian through finite difference\n");
	Hess_MM = CreateAndRunFiniteDifferenceHessian();
      }else
	Hess_MM = ReadFiniteDifferenceHessian();
    }else
      Hess_MM = ReadHessian(2);
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    cout << "AIFF Hessian not yet implemented \n";
    //Hess_MM = ComputeClusterMultipoleHessian();        
  }     
  else if (Params::Parameters().GetMMType() == 3) { // QChem
    Hess_MM = ReadHessian(1);             
  }
  else if (Params::Parameters().GetMMType() == 5) { // 
    Hess_MM = ReadCrystalFrequencies();
  }

  MM_Hess_Init = 1;
}

void Cluster::SetQMHessian() {
  
  if (Params::Parameters().GetQMType() == 5) { // Quantum Espresso
    //Hess_QM = ReadQuantumEspressoHessian();
    ReadQuantumEspressoHessian();
    //PrintHessian("New Full QM Hessian",Hess_QM);
    QM_Hess_Init = 1;
  }
  else if (Params::Parameters().GetQMType() == 8) { // DFTB
    //Hess_QM = ReadQuantumEspressoHessian();
    ReadDFTBHessian();
    //PrintHessian("New Full QM Hessian",Hess_QM);
    QM_Hess_Init = 1;
  }
  else{
    printf("Unknown QM Type = %i. Exiting...\n",Params::Parameters().GetQMType());
    exit(0);
    fflush(stdout);
  }

}

void Cluster::ReadQuantumEspressoHessian() {

  string path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  int Natoms = GetTotalNumberOfAtoms();
  Matrix hess(3*Natoms, 3*Natoms);
  //hess.Set();
  //hess.Print("testing hessian set");

  string filename;

  if(!Params::Parameters().UsePhonopyFreq()){
     // Set up the filename, with the full path.
    filename = path + "/dyn000.G";

    // Open the file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadQuantumEspressoHessian : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    string line; 
    while ( !infile.eof() ) {
      getline(infile,line);
      //printf("%s\n",line.c_str()); 
      // Search for the dynamical matrix  
      //Our final Hessian should be in units of hartrees/bohr^2
      if ( line.find("Dynamical  Matrix in cartesian axes") != -1 ) {
        getline(infile,line); // throw away header line
        getline(infile,line); // skip blank line
        getline(infile,line); // skip qpoints line
        for (int i=0;i<Natoms;i++) {
          for (int j=0;j<Natoms;j++) {
            getline(infile,line); //Always skip the atom-atom interaction line
            for (int m=0;m<3;m++) {
              getline(infile,line); //Always skip the atom-atom interaction line
              istringstream iss(line);
              for (int k=0;k<3;k++) {
                string imag;
                iss >> hess.Element(3*i+m,3*j+k); // Store the hessian elements
                iss >> imag;
                //printf("hess.Element(%d,%d): %f\nimag = %s\n",3*i+m,3*j+k,hess.Element(3*i+m,3*j+k),imag.c_str());
              }
            }

          }
        }
        break;
      }
    }

    infile.close();  // Close the file

    //PrintHessian("Before scaling",hess);

    //Convert units to hartree/bohr^2
    //Currently in Ry/bohr^2 so just divide by 2
    //hess.Scale(BohrToAng*BohrToAng*EVToHartrees);
    hess.Scale(0.5);

    //un-Mass weighted Hessian (i.e. the HMBI Hessian)
    Hess_HMBI = hess;
    //PrintHessian("Before mass scaling",hess);
    //printf("\n");
    for (int i=0;i<Natoms;i++) {
      double m1 = AtomicMasses[i];
      for (int j=0; j<Natoms;j++) {
        double m2 = AtomicMasses[j];
        double Gij = sqrt(m1*m2);
        for (int p=3*i;p<3*(i+1);p++) {
          for (int q=3*j;q<3*(j+1);q++) {
            hess(p,q) /= Gij;
          }       
        }
      }
    }

    //Mass-weighted Hessian
    Hess_QM = hess;
    //PrintHessian("New Full QM Hessian",hess);
    //fflush(stdout);
    //printf("\n");

  }
  else{
    // Set up the filename, with the full path.
    filename = path + "/qpoints_gamma.yaml";
  
    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadQuantumEspressoHessian : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    string line; 
    while ( !infile.eof() ) {
      getline(infile,line);
    //printf("%s\n",line.c_str()); 
      // Search for the dynamical matrix  
      //Our final Hessian should be in units of hartrees/bohr^2
      if ( line.substr(0,19)=="  dynamical_matrix:" ) {
        //getline(infile,line); // throw away header line
        for (int i=0;i<3*Natoms;i++) {
          getline(infile,line);
          istringstream iss(line);
          string tmp1,tmp2;
          iss >> tmp1; // throw away the matrix row index
          iss >> tmp2; // throw away the matrix column index
          for (int j=0;j<3*Natoms;j++) {
            string imag;
            iss >> hess.Element(i,j); // Store the hessian elements
            iss >> imag;
            //printf("hess.Element(%d,%d): %f\nimag = %s\n",i,j,hess.Element(i,j),imag.c_str());
          }
        }
        break;
      }
    }

    infile.close();  // Close the file

    //Currently in Ry/bohr^2 so just divide by 2
    hess.Scale(0.5);

    //Mass-weighted Hessian
    Hess_QM = hess;
    //PrintHessian("Before mass scaling",hess);
    //printf("\n");
    for (int i=0;i<Natoms;i++) {
      double m1 = AtomicMasses[i];
      for (int j=0; j<Natoms;j++) {
        double m2 = AtomicMasses[j];
        double Gij = sqrt(m1*m2);
        for (int p=3*i;p<3*(i+1);p++) {
          for (int q=3*j;q<3*(j+1);q++) {
            hess(p,q) *= Gij;
          }       
        }
      }
    }

    //un-Mass weighted Hessian (i.e. the HMBI Hessian)
    Hess_HMBI = hess;
    //PrintHessian("New Full QM Hessian",hess);
    //fflush(stdout);
    //printf("\n");
  }

}

void Cluster::ReadDFTBHessian() {

  string path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  int Natoms = GetTotalNumberOfAtoms();
  Matrix hess(3*Natoms, 3*Natoms);
  //hess.Set();
  //hess.Print("testing hessian set");

  printf("Correctly found DFTB Hessian read-in!\n");
  fflush(stdout);

  string filename;
  // Set up the filename, with the full path.
  filename = path + "/hessian_gamma.out";

  // Open the file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadDFTBHessian : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  int count=0;
  string line; 
  while ( !infile.eof() ) {
    Vector tmp(4);
    getline(infile,line);
    istringstream iss(line);
    iss >> tmp[0];
    iss >> tmp[1];
    iss >> tmp[2];
    iss >> tmp[3];

      for (int i=0;i<3*Natoms;i++) {
        for (int j=0;j<3*Natoms;j++) {
          //getline(infile,line); //Always skip the atom-atom interaction line
          if(count==4){
            //printf("Correctly trying to progress the line!\n");
            //printf("i = %i\tj = %i\n",i,j);
            //fflush(stdout);
            getline(infile,line); //Always skip the atom-atom interaction line
            //printf("Got new line!: %s\n",line.c_str());
            //fflush(stdout);
            istringstream iss(line);
            iss >> tmp[0];
            iss >> tmp[1];
            iss >> tmp[2];
            iss >> tmp[3];
            count=0;
          }
          //printf("line %s\n",line.c_str());
          //fflush(stdout);

          //iss >> hess.Element(i,j); // Store the hessian elements
          hess.Element(i,j) = tmp[count]; // Store the hessian elements
          count++;

          //printf("hess.Element(%d,%d): %f\n",i,j,hess.Element(i,j));
        }
      }
      break;

  }

  infile.close();  // Close the file

  //PrintHessian("Before scaling",hess);
  //fflush(stdout);
  //exit(0);
  //Convert units to hartree/bohr^2
  //Currently in Ry/bohr^2 so just divide by 2
  //hess.Scale(BohrToAng*BohrToAng*EVToHartrees);
  //hess.Scale(0.5);

  //un-Mass weighted Hessian (i.e. the HMBI Hessian)
  Hess_HMBI = hess;
  //PrintHessian("Before mass scaling",hess);
  //printf("\n");
  for (int i=0;i<Natoms;i++) {
    double m1 = AtomicMasses[i];
    for (int j=0; j<Natoms;j++) {
      double m2 = AtomicMasses[j];
      double Gij = sqrt(m1*m2);
      for (int p=3*i;p<3*(i+1);p++) {
        for (int q=3*j;q<3*(j+1);q++) {
          hess(p,q) /= Gij;
        }       
      }
    }
  }

  //Mass-weighted Hessian
  Hess_QM = hess;
  //PrintHessian("New Full QM Hessian",hess);
  //fflush(stdout);
  //printf("\n");

}

// Main routine for reading in hessian
// type 1 = qchem style, type 2 = tinker style
Matrix Cluster::ReadHessian(int type) {

  //string path = Params::Parameters().GetMMPath();
  string path = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  int Natoms = GetTotalNumberOfAtoms();
  Matrix hess(3*Natoms, 3*Natoms);

  if ( Params::Parameters().NeglectManyBody() )
    hess.Set();
  else {

    string filename;
    // Set up the filename, with the full path.  File is 'full.freq'
    if (type == 2) // Tinker MM job
      filename = path + "/full.freq";
    else
      filename = path + "/full.out";

    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadHessian : Cannot open file '%s'\n",
             filename.c_str());
      exit(1);
    }

    if (type == 2) { // look in the tinker .freq file
      string line;
      while ( !infile.eof() ) {
        getline(infile,line);
        // Search for the SCF hessian
        if ( line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):" ) {
          //getline(infile,line); // throw away header line
          for (int i=0;i<3*Natoms;i++) {
            for (int j=0;j<3*Natoms;j++) {
              getline(infile,line);
              istringstream iss(line);
              string tmp1,tmp2;
              iss >> tmp1; // throw away the matrix row index
              iss >> tmp2; // throw away the matrix column index
              iss >> hess.Element(i,j); // Store the hessian elements
            }
          }
          break;
        }
      }
    }

    else { //look in the qchem output
      string line;
      while ( !infile.eof() ) {
        getline(infile,line);
        // Search for the SCF hessian
/*
        if ( line.substr(0,30)==" Mass-Weighted Hessian Matrix:" ) {
          for (int k=0;k<(3*Natoms-3*Natoms%6)/6;k++) {
            getline(infile,line); // throw away the 1st blank line
            getline(infile,line); // throw away the 2nd blank line
            //while (line.substr(0,15)!=" Gradient time:") {

            for (int i=0;i<3*Natoms;i++) {
              getline(infile,line);
              istringstream iss(line);
              for (int j=6*k;j<6*k+6;j++) {                 
                iss >> hess.Element(i,j); // Store the hessian elements
              }
            }  
          }  
          getline(infile,line); // throw away the 1st blank line
          getline(infile,line); // check if this line signals end of hessian print or not
          if (line.substr(0,52) != " Translations and Rotations Projected Out of Hessian") {
            for (int i=0;i<3*Natoms;i++) {
              getline(infile,line);
              istringstream iss(line);
              string tmp;
              iss >> tmp; // throw away the atom index
              for (int j=(3*Natoms-3*Natoms%6);j<3*Natoms;j++) {
                iss >> hess.Element(i,j); // Store the hessian elements
              }
            }
            break;
          }
*/
      if ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) {
        for (int k=0;k<(3*Natoms-3*Natoms%6)/6;k++) {         
          getline(infile,line); // throw away header line                  
          //while (line.substr(0,15)!=" Gradient time:") {

          for (int i=0;i<3*Natoms;i++) {
            getline(infile,line);
            istringstream iss(line);
            string tmp;
            iss >> tmp; // throw away the atom index
            for (int j=6*k;j<6*k+6;j++) {
              iss >> hess.Element(i,j); // Store the hessian elements
            }                
          }                
        }
  	getline(infile,line); // check if this line signals end of hessian print or not
        if (line.substr(0,15) != " Gradient time:") {             
          for (int i=0;i<3*Natoms;i++) {             
            getline(infile,line);
            istringstream iss(line);
            string tmp;
            iss >> tmp; // throw away the atom index
            for (int j=(3*Natoms-3*Natoms%6);j<3*Natoms;j++) {
              iss >> hess.Element(i,j); // Store the hessian elements
            }
          }
          break;
        }    

          else {
            break;
          }
        }
      }
    }


/*
    else { //look in the qchem output
      string line;
      while ( !infile.eof() ) {
        getline(infile,line);
        // Search for the SCF hessian
        if ( line==" Hessian of the SCF Energy" ) {
          getline(infile,line); // throw away header line
          while (line.substr(0,15)!=" Gradient time:") {               

            for (int i=0;i<3*Natoms;i++) {     
              getline(infile,line);
              istringstream iss(line);
              string tmp;
              iss >> tmp; // throw away the atom index
              for (int j=0;j<6;j++) {
                iss >> hess.Element(i,j); // Store the hessian elements
              }
            }                 
          }
          break;               
        }      
      }
    }
*/
    infile.close();  // Close the file

  }

  //PrintHessian("New Full MM Hessian",hess);
  return hess;
}

Matrix Cluster::ReadCrystalFrequencies(){

  int Natoms = GetTotalNumberOfAtoms();

  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //  path = Params::Parameters().GetMMPath();
  //  path = Params::Parameters().GetHessianMMPath();
  
  //eigenvector but stored in hess matrix
  Matrix hess(3*Natoms,3*Natoms);

    string filename;
    // Set up the filename, with the full path.  File is 'full.out'
    filename = path + "/full.out";

    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Cluster::ReadCrystalFrequencies : Cannot open file '%s'\n",
             filename.c_str());
      exit(1);
    }

    string line;
    int iFreq = 0;
    while ( !infile.eof() ) {

      getline(infile,line);
      if(line.substr(0,13) == " FREQ(CM**-1)"){	
	int j;
	if(iFreq - 3*Natoms <= 6)
	  j = 6;
	else
	  j = 3*Natoms - iFreq;

	istringstream iss(line);
	string tmp; 
	iss >> tmp; // throw away first string
	for(int i=0;i<j;i++){
	  iss >> Freqs[iFreq+i];
	  //printf("Freq[%i] = %f\n",iFreq,Freqs[iFreq]);
 
	}
	getline(infile,line);

	//Read Eigenvectors;
	for(int i=0;i<Natoms;i++){

	  //X component of eigenvector
	  getline(infile,line);
	  //printf("%s\n",line.c_str());
	  istringstream iss(line);
	  string tmp;  
	  iss >> tmp; // throw away "AT." string
	  iss >> tmp; // throw away atom index
	  iss >> tmp; // throw away element symbol
	  iss >> tmp; // throw away "X" string
	  for(int k=0;k<j;k++){
	    iss >> hess(3*i,iFreq+k);
	    //printf("hess = %f\n",hess(3*i,iFreq+k));
	  }

	  //Y component of eigenvector
	  getline(infile,line);
	  iss.clear();
	  iss.str(line);
	  //printf("%s\n",line.c_str());
	  iss >> tmp; // throw away "Y" string
	  for(int k=0;k<j;k++){
	    iss >> hess(3*i+1,iFreq+k);
	    //printf("hess = %f\n",hess(3*i+1,iFreq+k));
	  }
		  
	  //Z component of eigenvector
	  getline(infile,line);
	  iss.clear();
	  iss.str(line);
	  //printf("%s\n",line.c_str());
	  iss >> tmp; // throw away "Z" string
	  for(int k=0;k<j;k++){
	    iss >> hess(3*i+2,iFreq+k);
	    //printf("hess = %f\n",hess(3*i+2,iFreq+k));
	  } 

	}

	iFreq += 6;
      }
    }
    //Freqs.Print("Freqs");
    //hess.PrintHessian("hess");

    //exit(0);
    infile.close();  // Close the file
    return hess;
}

//Gets finite difference Hessian for tinker
Matrix Cluster::ReadFiniteDifferenceHessian(){

  int Natoms = GetTotalNumberOfAtoms();

  string path = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //  path = Params::Parameters().GetMMPath();
  //  path = Params::Parameters().GetHessianMMPath();
  path += "/full";
  
  Matrix hess(3*Natoms,3*Natoms);
  
  
  if ( Params::Parameters().NeglectManyBody() )
    hess.Set();
  else {
    for(int i=0;i<3*Natoms;i++){
      Vector GradPlus(3*Natoms);
      Vector GradMinus(3*Natoms);
      
      //creating string from i for the filenames
      char num[10];
      sprintf(num,"%d",i);

      //Set path for the plus file
      string  filename = path + "/full" + num + "+.force";
      // Open the force file
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Cluster::ReadFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      //Set values for the plus grad
      // Read in the data - search for the "Nuclear forces:" string
      string line;
      while ( !infile.eof() ) {
	getline(infile,line);
	// Search for final q-chem energy
	if ( line==" Nuclear forces:" ) {
	  getline(infile,line); // throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	    iss >> GradPlus[3*j+k]; // Store the gradient elements
	    }
	  }
	  break;
	}
      }
      
      infile.close();  // Close the file

      //Set values for the minus grad
      //Set path for the plus file
      filename = path + "/full" + num + "-.force";
      // Open the force file
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Cluster::ReadFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      //Set values for the plus grad
      // Read in the data - search for the "Nuclear forces:" string
      while ( !infile.eof() ) {
	getline(infile,line);
	// Search for final q-chem energy
	if ( line==" Nuclear forces:" ) {
	  getline(infile,line); // throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> GradMinus[3*j+k]; // Store the gradient elements
	    }
	  }
	  break;
	}
      }
      
      infile.close();  // Close the file
      
      GradMinus -= GradPlus;
      GradMinus.Scale(1.0/(2*0.001*AngToBohr));

      //for(int j=0;j<3*Natoms;j++)
      //hess(j,i) = GradMinus[j];
      hess.SetColumnVector(GradMinus,i);
      
    }//end of loop over i
  }		  

  return hess;
}

//Create, Run and Read the Finite Difference Hessian with tinker 
//This was used before Dr. Beran corrected the polarization in tinker.
Matrix Cluster::CreateAndRunFiniteDifferenceHessian(){


  //Previous created
  //CreateFiniteDifferenceTinkerJobs(3);

  int Natoms = GetTotalNumberOfAtoms();

  //  path = Params::Parameters().GetMMPath();
  //  path = Params::Parameters().GetHessianMMPath();
  string path = Params::Parameters().GetHessianMMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  path += "/full";

  Matrix hess(3*Natoms,3*Natoms);
  
  
  if ( Params::Parameters().NeglectManyBody() )
    hess.Set();
  else {


    //Loop to run and read jobs
    for(int i=0;i<3*Natoms;i++){
      Vector GradPlus(3*Natoms);
      Vector GradMinus(3*Natoms);
      
      //creating string from i for the filenames
      char num[10];
      sprintf(num,"%d",i);


      //Run Plus jobs
      string cmd = RunFiniteDifferenceTinkerJob(i,true);
      system(cmd.c_str()); 

      //Set path for the plus file
      string  filename = path + "/full" + num + "+.force";
      // Open the force file
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Cluster::ReadFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      //Set values for the plus grad
      // Read in the data - search for the "Nuclear forces:" string
      string line;
      while ( !infile.eof() ) {
	getline(infile,line);
	// Search for final q-chem energy
	if ( line==" Nuclear forces:" ) {
	  getline(infile,line); // throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	    iss >> GradPlus[3*j+k]; // Store the gradient elements
	    }
	  }
	  break;
	}
      }
      
      infile.close();  // Close the file


      //run minus job
      cmd = RunFiniteDifferenceTinkerJob(i,false);
      system(cmd.c_str()); 

      //Set values for the minus grad
      //Set path for the plus file
      filename = path + "/full" + num + "-.force";

      // Open the force file
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Cluster::ReadFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      //Set values for the plus grad
      // Read in the data - search for the "Nuclear forces:" string
      while ( !infile.eof() ) {
	getline(infile,line);
	// Search for final q-chem energy
	if ( line==" Nuclear forces:" ) {
	  getline(infile,line); // throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> GradMinus[3*j+k]; // Store the gradient elements
	    }
	  }
	  break;
	}
      }
      
      infile.close();  // Close the file
      
      GradMinus -= GradPlus;
      GradMinus.Scale(1.0/(2*0.001*AngToBohr));

      //for(int j=0;j<3*Natoms;j++)
      //hess(j,i) = GradMinus[j];
      hess.SetColumnVector(GradMinus,i);
      
    }//end of loop over i
  }		  

  return hess;

}

void Cluster::ComputeHarmonicFrequencies(Matrix& Hessian) {
  

  Matrix Hess;
  int Nimag = 0;
  int Nfreq;
  //Crystal Frequency read from crystal output
  if(Params::Parameters().GetMMType() == 5){
    Hess = Hessian;

    Nfreq = Freqs.GetLength(); 
    //Vector freqs = Freqs;
    // Count number of negative eigenvalues (imaginary freqs)
    for (int i=0; i<Nfreq;i++) {
      if (Freqs[i] < 0.0) {
        Nimag++;
      }
    }


  }
  else{
    printf("\nComputing harmonic vibrational frequencies.\n");
    
    //Hessian already weighted
    //Matrix mwHess = ComputeMassWeightedHessian(Hessian);
    
    // Project out the translational & rotational modes
    printf("Projecting out translational and rotational modes\n");
    Hess = ProjectOutTranslationAndRotation(Hessian);

    // Diagonalize the mass-weighted Hessian
    Freqs = Hess.Diagonalize();
    
    // Print out the matrix of eigenvectors
    //PrintHessian("MW Eigenvectors\n",Hess);
    Hess = UnMassWeightedEigenVectors(Hess);
    
    //PrintHessian("UNMW Eigenvectors\n",Hess);
    SetRMode(Hess); // Watit - use UnmassHessian as Normal Modes to compute intensity 
       
    for ( int i=0; i<Hess.GetCols();i++) {
      for ( int j=0; j<Hess.GetRows();j++) {
	Hess.Element(j,i) /= Hess.GetColumnVector(i).Norm();
      }
    }
    //Zero out any really tiny frequencies & Convert the frequencies to wavenumbers // Watit
    Nzero = 0;
    Nfreq = Freqs.GetLength();
    double c1 = 5.89141e-7;
    double c2 = 3.94562;//convert wavenumber in atomic units into cm^-1
    for (int i=0; i<Nfreq;i++) {
      //Zero out any really tiny frequencies
      if (fabs(Freqs[i]) < 1.0e-6) {
        Freqs[i] = 0.0;
        Nzero += 1; 
      }
      //Count number of negative eigenvalues (imaginary freqs)
      if (Freqs[i] < 0.0) {
        Nimag++;
        double imag = sqrt(fabs(Freqs[i])/c1)*c2;
        Freqs[i] = -imag;
      }
      if (Freqs[i] > 0.0 ) {
        Freqs[i] = sqrt(fabs(Freqs[i])/c1)*c2;
      }
    }
    //Vector freqs = Freqs;
  }

	// Watit - Compute IR and Raman Intensities
	//ReadHessianData();
	if (!Params::Parameters().DoQuasiHarmonic()&&Params::Parameters().GetQMType()==3) {

		ReadIntensityData();

        	Matrix RMode = GetRMode();
        	RMode.Print("Normal Modes");

		int NVib = Nfreq - Nzero;
		int NIRActive=0;
		int NRamanActive=0;
		Matrix Vib(NVib,3);
		Matrix RModeVector(Nfreq/3,3);
		int k = 0;
		for (int i=0; i<Nfreq;i++) {
			if (Freqs[i] != 0.0 ) {
				Vib(k,0) = Freqs[i];
				Vib(k,1) = IntIR[i];
				Vib(k,2) = IntRaman[i];
				k++;
   			}
			if (Freqs[i] < 0.0) {
				for (int j=0; j<Nfreq/3; j++) {
					for ( int l=0; l<3; l++) {
						RModeVector(j,l) = GetRMode_value(l+3*j,i);
					}
				}

                                if (!Params::Parameters().NoFreqDisplacement()) {
                                  printf("IMAGINARY FREQUENCY FOUND!!! Exiting HMBI\n");
    				  AdjustCoord(i,RModeVector);
                                }

			}
  		}
		RModeVector.Print("img Normal Mode Vector");

  		SetVib(Vib);
		PlotIntensityData(NVib);

  		printf("\n***************************************************************\n");
		printf("*                   Vibrational Analysis                      *\n");
		printf("***************************************************************\n\n");
		if (Params::Parameters().DoRaman()) {
			printf("  Frequencies (cm-1)        IR Intensity     Raman Intensity\n");
			for(int i=0; i<NVib; i++) {
				for(int j=0; j<3; j++) {
					printf("%20.3f",Vib(i,j));
 				}
				printf("\n");
			}
		} else {
			printf("  Frequencies (cm-1)        IR Intensity\n");
			for(int i=0; i<NVib; i++) {
				for(int j=0; j<2; j++) {
					printf("%20.3f",Vib(i,j));
				}
				printf("\n");
			}
		}
		printf("\n");
		printf("***************************************************************\n");
	
	}


  Vector orig_freqs( Freqs );
  Freqs.SortLargestFirst();
  Vector Orig_Index(Nfreq);
  Orig_Index = CompareSortedOrigVecs(orig_freqs,Freqs);
  Hess = SortColumns(Hess,Orig_Index);
  
  //Unlike previous versions of HMBI, the columns of Hess are now reorginized in order of corresponding freqs.
  // Still need an order list for PrintNormalModes
  Orig_Index.IncrementalEntries(true);

  //Freqs.Print("Freqs");
  //orig_freqs.Print("orig_freqs");
  //Orig_Index.Print("Orig_Index");
  //exit(0);
  bool molecule = Params::Parameters().IsMolecule();

  // Count number of negative eigenvalues (imaginary freqs)
  for (int i=0; i<Nfreq;i++) {
    if (Freqs[i] < 0.0) {
      printf("Imaginary frequency: %.2f\n",Freqs[i]);
    }
  }
  printf("%d imaginary frequencies found\n\n",Nimag);
  int Nonzero = Freqs.Nonzero();

	if (!Params::Parameters().DoRaman()&&!(Params::Parameters().GetQMType()==3)) {
		//Count number of negative eigenvalues (imaginary freqs)
		Matrix RModeVector(Nfreq/3,3);
		for (int i=0;i<Nfreq;i++) {
			if (Freqs[i] < 0.0) {
				printf("Imaginary frequency: %.2f\n",Freqs[i]);
				for (int j=0; j<Nfreq/3; j++) {
					for ( int l=0; l<3; l++) {
						RModeVector(j,l) = GetRMode_value(l+3*j,i);
					}
				}
                                if (!molecule && !Params::Parameters().NoFreqDisplacement()){
				  AdjustCoord(i,RModeVector);
                                }
			}
		}
	  	printf("Harmonic Vibrational frequencies (cm-1): %d real frequencies\n",Nonzero);
	  	for (int i=0;i<Nonzero/6+1;i++) {
    			for (int j=0;j<6;j++) {
      				if ((6*i+j) < Nonzero && Freqs[6*i+j]>0.0) {
        				printf("  %8.2f",Freqs[6*i+j]);
	      			} else break;              
    			}
    			printf("\n");
  		}
	}


  //printf("\n Ref cell = %i\n",Quasiharmonic::quasiharmonic().IsReferenceCellInitialized());



  

  //Saving reference modes for the quasiharmonic approximating
  if(Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().IsSupercellJob()){
    if(!Quasiharmonic::quasiharmonic().IsReferenceCellInitialized()){
      Quasiharmonic::quasiharmonic().SetReferenceMode(Hess);
      PrintNormalModes(Freqs, Nfreq, Hess, Nonzero, Nimag, Orig_Index);
    }
    else { //Reordering freq so they are the same order as in the reference cell
      printf("\n\nReference was set\n");
      orig_freqs = Freqs;
      Freqs = Quasiharmonic::quasiharmonic().OrderFrequencies(Freqs,Hess);
      Orig_Index = CompareSortedOrigVecs(orig_freqs,Freqs);
      string File_Name = Quasiharmonic::quasiharmonic().GetHessianType();
      PrintNormalModes(Freqs, Nfreq, Hess, Nonzero, Nimag, Orig_Index,File_Name);
      //exit(0);
   }
  }
  else
    PrintNormalModes(Freqs, Nfreq, Hess, Nonzero, Nimag, Orig_Index);

  ComputeVibrationalContribution(Freqs);

  //Zero point energy is part of the vibrational contribution
  //double ZPE = ComputeZeroPointEnergy(Freqs);     


  //if doing quasiharmonic approximation and have imaginary frequencies, stop calculation unless its a translational mode.
  int nonpos = Nfreq - Nonzero;
  //printf("Nfreq = %i Nonzero = %i nonpos = %i\n",Nfreq,Nonzero,nonpos);
  //Freqs.Print("Freq");
  if(Params::Parameters().DoQuasiHarmonic()){
    if(Nimag > 0 && nonpos > 3){
      printf("ERROR::Cluster::ComputeHarmonicFrequencies: Imaginary frequencies produced during Quasi-harmonic calculations. Cell is not relaxed\n. See molden file.\n");
      printf("Nimag = %i\n",Nimag);
      exit(0);
    }
  }


}

Matrix Cluster::ComputeMassWeightedHessian(Matrix& Hessian) {

  int Natoms = GetTotalNumberOfAtoms();
  //int Natoms = UniqueAtoms;
  Matrix mwHess = Hessian;    

  /*
  int AtomCounter = 0;
  double* UniqueMasses = new double[UniqueAtoms];
  for(int imon=1;imon<=NMon;imon++)
    if(Monomers[imon].GetSymmetryFactor()!=0)
      for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms();iatom++){
	printf("setting UniqueMasses %i\n",AtomCounter);
	printf("m%i atom %i\n",imon,iatom);
	fflush(stdout);
	UniqueMasses[AtomCounter] = Monomers[imon].GetAtomicMass(iatom);
	AtomCounter++;
      }
   
   printf("UniqueMasses = [");
  for(int i=0;i<UniqueAtoms;i++)
    printf("%f ",UniqueMasses[i]);
  printf("]\n");

  for(int i=0;i<UniqueAtoms;i++){
    double m1 = UniqueMasses[i];
    for (int j=0; j<UniqueAtoms; j++){
      double m2 = UniqueMasses[j];
      double Gij = 1.0/sqrt(m1*m2);
      for (int p=3*i; p<3*(i+1); p++)
	for(int q=3*j;q<3*(j+1);q++)
	  mwHess(p,q) *= Gij;
    }
  }
  */
      

  for (int i=0;i<Natoms;i++) {
    double m1 = AtomicMasses[i];
    for (int j=0; j<Natoms;j++) {
      double m2 = AtomicMasses[j];
      double Gij = 1.0/sqrt(m1*m2);
      for (int p=3*i;p<3*(i+1);p++) {
        for (int q=3*j;q<3*(j+1);q++) {
          mwHess(p,q) *= Gij;
        }       
      }
    }
  }

  //delete [] UniqueMasses;
  
  return mwHess;
}

Matrix Cluster::UnMassWeightedEigenVectors(Matrix& Hessian) {
  int Natoms = GetTotalNumberOfAtoms();    
  for (int i=0;i<Natoms;i++) {
    double m1 = AtomicMasses[i];
    for (int j=0; j<3*Natoms;j++) {
      for (int p=3*i;p<3*(i+1);p++) {
        Hessian.Element(p,j) *= 1/sqrt(m1);
      }
    }
  }      
  return Hessian;
}  

// Projects out translation & rotation from the mass-weighted hessian.
Matrix Cluster::ProjectOutTranslationAndRotation(Matrix& Hessian) {

  int Natoms = GetTotalNumberOfAtoms();
  bool molecule = Params::Parameters().IsMolecule();
  Matrix Pvib(3*Natoms,true);
  /* Construct the Translational Projectors | Tx >, | Ty >, | Tz > */
  Matrix Tx(3*Natoms,1), Ty(3*Natoms,1), Tz(3*Natoms,1);
  Tx.Set(); Ty.Set(); Tz.Set();           

  for (int i=0;i<Natoms;i++) {
    Tx(3*i,0) = sqrt(AtomicMasses[i])/sqrt(TotalMass);
    Ty(3*i+1,0) = sqrt(AtomicMasses[i])/sqrt(TotalMass);
    Tz(3*i+2,0) = sqrt(AtomicMasses[i])/sqrt(TotalMass);
  }         

  // Form outer products |Tx><Tx|, |Ty><Ty|, & |Tz><Tz|
  Matrix TxTx = Tx.Multiply(Tx,3);
  Matrix TyTy = Ty.Multiply(Ty,3);
  Matrix TzTz = Tz.Multiply(Tz,3); 

  // Ptrans = |Tx><Tx| + |Ty><Ty| + |Tz><Tz|
  Matrix Ptrans = TxTx;
  Ptrans += TyTy;
  Ptrans += TzTz;             
  Pvib -= Ptrans;   

  if ( (!Params::Parameters().IsPeriodic() && (!Params::Parameters().UseFullQMOnly()) || molecule) ) { //periodic systems are not rotation invariant //JLM

    /* Construct the Rotational Projector */
    Matrix Ra(3*Natoms,1), Rb(3*Natoms,1), Rc(3*Natoms,1);     
    Ra.Set(); Rb.Set(); Rc.Set();

    // Compute and diagonalize the inertia tensor
    Matrix Inertia = ComputeInertiaTensor();
    Vector Imom(3);             
    Imom = Inertia.Diagonalize();  

    // Print out the inertia tensor results
    printf("\nMoments and Principal Axes of Inertia (atomic units):\n");
    printf("               1          2         3\n");
    printf("Moments:  %9.5f  %9.5f  %9.5f\n",Imom[0], Imom[1], Imom[2]);
    printf("   X      %9.5f  %9.5f  %9.5f\n", Inertia(0,0), Inertia(0,1), Inertia(0,2));
    printf("   Y      %9.5f  %9.5f  %9.5f\n", Inertia(1,0), Inertia(1,1), Inertia(1,2));
    printf("   Z      %9.5f  %9.5f  %9.5f\n", Inertia(2,0), Inertia(2,1), Inertia(2,2));
    printf("\n");

    // Construct |Ra>, |Rb>, and |Rc>
    // |Ra>[i] = sqrt(m_i)*r_i x (Vc x Vb) = sqrt(m_i)*r_i x (-Va)
    // |Rb>[i] = sqrt(m_i)*r_i x (Va x Vc) = sqrt(m_i)*r_i x (-Vb)
    // |Rc>[i] = sqrt(m_i)*r_i x (Vb x Va) = sqrt(m_i)*r_i x (-Vc)      
    // r_i is the xyz coordinates of the i-th atom.  Ra/Rb/Rc are 3*N long.

    // Grab the Principle axes vectors, and take their negative.
    // Need to add in a check to handle linear molecules correctly
    Vector Va;
    //Checking up to 3 decimal comparison match
    if( trunc(1000. * Imom[0]) == trunc(1000. * 0.000) ) {
      Va.Initialize(3);
    }
    else {
      Va = Inertia.GetColumnVector(0);
    }
    Vector Vb = Inertia.GetColumnVector(1);
    Vector Vc = Inertia.GetColumnVector(2);               
    Va.Scale(-1.0);                         
    Vb.Scale(-1.0);
    Vc.Scale(-1.0);               

    // form Ra, Rb, Rc.
    for (int iatom=0;iatom<Natoms;iatom++) {              
      Vector tmp(3), R(3);
      double mass = AtomicMasses[iatom];
      for (int i=0;i<3;i++) {
        R[i] = (AtomicCoordinates[3*iatom+i] - CenterOfMass[i])*AngToBohr;
      }
      Vector tmpA = R.CrossProduct(Va);
      Vector tmpB = R.CrossProduct(Vb);
      Vector tmpC = R.CrossProduct(Vc);       

      tmpA *= sqrt(mass);
      tmpB *= sqrt(mass);       
      tmpC *= sqrt(mass);

      for (int i=0;i<3;i++) {               
        Ra(3*iatom+i,0) = tmpA[i];
        Rb(3*iatom+i,0) = tmpB[i];
        Rc(3*iatom+i,0) = tmpC[i];
      }
    }  
    // Eliminate small values in Ra/Rb/Rc before we normalize.  Can
    // occasionally have numerical issues if we don't.
    for (int i=0;i<3*Natoms;i++) {         
      if (fabs(Ra(i,0)) < 1.0e-8) Ra(i,0) = 0.0;
      if (fabs(Rb(i,0)) < 1.0e-8) Rb(i,0) = 0.0;      
      if (fabs(Rc(i,0)) < 1.0e-8) Rc(i,0) = 0.0;
    } 

    // Normalize the vectors
    Vector tmp = Ra.GetColumnVector(0);   
    if(tmp.Norm() != 0.0)
     tmp.Normalize();
    Ra.SetColumnVector(tmp,0);  
     tmp = Rb.GetColumnVector(0);
    if(tmp.Norm() != 0.0)
     tmp.Normalize();
    Rb.SetColumnVector(tmp,0);
    tmp = Rc.GetColumnVector(0);
    if(tmp.Norm() != 0.0)
      tmp.Normalize();
    Rc.SetColumnVector(tmp,0);          

    // Form outer products |Ra><Ra|, |Rb><Rb|, & |Rc><Rc|
    Matrix RaRa = Ra.Multiply(Ra,3);
    Matrix RbRb = Rb.Multiply(Rb,3);
    Matrix RcRc = Rc.Multiply(Rc,3);    

    // Prot = |Ra><Ra| + |Rb><Rb| + |Rc><Rc|              
    Matrix Prot = RaRa;   
    Prot += RbRb;
    Prot += RcRc;            
    //Prot.Print("Final Rot projector");
    Pvib -= Prot;               //Periodic systems are not rotation-invariant
  }  //End periodic system if statement

  /* Create projector that leaves only pure vibration:
     Pvib = I - Ptrans - Prot;       
  */
  // Apply the projector to project out translation & rotation:  Hvib = P'HP

  /*
  //print dimensions 
  printf("Dimensions of Hessian %i %i\n",
  	 Hessian.GetRows(),Hessian.GetCols());
  //printf("Dimensions of Pvib %i %i\n",
	 Pvib.GetRows(),Pvib.GetCols());
  */

  //Pvib.PrintHessian("Pvib");

  Matrix HP = Hessian.Multiply(Pvib,1);  // H*P

  Matrix Hvib = Pvib.Multiply(HP,2); // P'*(H*P)

  return Hvib;                  
}

Matrix Cluster::ComputeInertiaTensor() {

  //  int Natoms = GetTotalNumberOfAtoms();
  int Natoms = GetNumberOfAtoms(Params::Parameters().GetFrequencyForMonomer());
  Matrix Inertia(3,3); // for the inertia tensor
  Inertia.Set();

  for (int iatom=0;iatom<Natoms;iatom++) {
    // grab coordinates for this atom & shift to center-of-mass origin
    double x = (AtomicCoordinates[3*iatom] - CenterOfMass[0])*AngToBohr;
    double y = (AtomicCoordinates[3*iatom+1] - CenterOfMass[1])*AngToBohr;
    double z = (AtomicCoordinates[3*iatom+2] - CenterOfMass[2])*AngToBohr;
    double m = AtomicMasses[iatom];
    
    // Construct the inertia tensor in the center-of-mass coordinate system
    Inertia(0,0) += m*(y*y + z*z);
    Inertia(1,1) += m*(x*x + z*z);
    Inertia(2,2) += m*(x*x + y*y);
    
    Inertia(0,1) -= m*x*y;
    Inertia(0,2) -= m*x*z;
    Inertia(1,2) -= m*y*z;
  }
  // Symmetrize the matrix
  Inertia(1,0) = Inertia(0,1);
  Inertia(2,0) = Inertia(0,2);
  Inertia(2,1) = Inertia(1,2);

  return Inertia;
}

void Cluster::ComputeQuasiHarmonicFrequencies(){

  if(Params::Parameters().UseVolumeGruneisen()){
    double reference_volume = Quasiharmonic::quasiharmonic().GetReferenceVolume();
    Vector reference_frequencies = Quasiharmonic::quasiharmonic().GetReferenceFrequencies();
    Vector Gruneisen_parameters = Quasiharmonic::quasiharmonic().GetGruneisenParameters();
    
   //change frequency length to include frequency past the gamma point
    if(Params::Parameters().IsSupercellJob())
      Freqs.Initialize(reference_frequencies.GetLength());
    
    //for(int i = 0; i < 3*GetTotalNumberOfAtoms();i++){
    for(int i = 0; i < reference_frequencies.GetLength();i++){
      if(fabs(Gruneisen_parameters[i]) < 1.0e-12 )
        Freqs[i] = 0.0;
      else
        Freqs[i] = reference_frequencies[i]*(pow(cell_volume/reference_volume,-Gruneisen_parameters[i]) );
    }
    
    //Freqs.Print("QuasiHarmonic Freq");
    
  }else 
    Freqs = Quasiharmonic::quasiharmonic().ComputeAnisotropicFrequencies();
  

  //if doing a optimization only print frequencies at the end of optimization
  if(!Params::Parameters().DoForces() && (Params::Parameters().UseDLFind()
                                        ||Params::Parameters().UseKNitro()
                                        ||Params::Parameters().UseConjugateGradient()) )
    PrintFrequencies();
  /*
  int Nonzero = Freqs.Nonzero();
  printf("QuasiHarmonic Vibrational frequencies (cm-1): %d real frequencies\n",Nonzero);
  for (int i=0;i<Nonzero/6+1;i++) {
    for (int j=0;j<6;j++) {
      if ((6*i+j) < Nonzero && Freqs[6*i+j]>0.0) {
        printf("  %8.2f",Freqs[6*i+j]);
      }
      else
        break;              
    }
    printf("\n");
  }
  */

  //Freqs.Print("Freqs");
  ComputeVibrationalContribution(Freqs);


  //Zero point energy is part of the vibrational contribution
  //double ZPE = ComputeZeroPointEnergy(Freqs);

}
// Returns vibrational energy contribution in kJ/mol
void Cluster::ComputeVibrationalContribution(Vector freqs){
  
  double temperature = Params::Parameters().GetTemperature();
  double inverse_temp = 1/(k_boltz*temperature);
  //printf("T = %4.0f K  1/(kT) = %8.3f J^-1\n",temperature,inverse_temp);
  double vib_energy = 0.0;
  double vib_entropy = 0.0;
  
  //constant volume heat capacity
  double heat_capacity_v = 0.0;


  //macroscopic gruneisen parameter, not the same gruneisen parameter used by quasiharmonic approximation
  double macro_gruneisen = 0.0;

  for (int i=0;i<freqs.GetLength();i++)
    if(freqs[i] > 0 && temperature != 0){


      vib_energy += Na/inverse_temp/1000*log(1-exp(-freqs[i]*c_light*100*h_planck*inverse_temp));
      //printf("Single vib_energy = %8.5f in kJ/mol\n",Na/inverse_temp/1000*log(1-exp(-freqs[i]*c_light*100*h_planck*inverse_temp)));
 

            //Second calculate the entropy "S"
            //1st term
            vib_entropy += -1.0*R_gas_constant*log(1-exp(-1 * h_planck * c_light*100 * freqs[i] * Na / R_gas_constant / temperature) );
            //2nd term (check sign)
       	    vib_entropy += 1.0 / temperature * h_planck * c_light*100.0 * freqs[i] *Na
                       *exp(-1.0* h_planck * c_light*100.0 * freqs[i] * Na / R_gas_constant / temperature)
                       /( -1*exp(-1* h_planck * c_light*100.0 * freqs[i] * Na / R_gas_constant / temperature) + 1.0);

      double cv_mode = pow((c_light*100.0 * freqs[i]*h_planck/temperature/ (exp(freqs[i]*c_light*100*h_planck*inverse_temp ) - 1 )),2)
	              * Na/k_boltz * exp( freqs[i]*c_light*100*h_planck*inverse_temp);

      heat_capacity_v += cv_mode;

      if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable()){
	double micro_gruneisen = Quasiharmonic::quasiharmonic().GetGruneisenParameters()[i];
        macro_gruneisen += micro_gruneisen*cv_mode;
       }


      //vib_entropy += h_planck*freqs[i]*c_light*100.0*Na/temperature/(exp(freqs[i]*c_light*100*h_planck*inverse_temp)-1);
      //vib_entropy -= R_gas_constant*log(1-exp(-freqs[i]*c_light*100*h_planck*inverse_temp));      
       //vib_entropy += h_planck*freqs[i]*c_light*100.0*Na/temperature/(exp(freqs[i]*c_light*100*h_planck*inverse_temp)-1);
      //printf("vib_entropy = %f\n",h_planck*freqs[i]*c_light*100.0*Na/temperature/(1-exp(freqs[i]*c_light*100*h_planck*inverse_temp)));
    }

  //vib_entropy -= vib_energy*1000/temperature;

  if(!Params::Parameters().GetSuperCellTotalKPoint())
    printf("\nHelmholtz vibrational energy (at %.0f K) excluding ZPE = %8.3f kJ/mol\n",
	   temperature,vib_energy);
  else
    printf("\nHelmholtz vibrational energy (at %.0f K) excluding ZPE = %8.3f kJ/mol\n",
	   temperature,vib_energy/double(Params::Parameters().GetSuperCellTotalKPoint()));

  vib_energy += ComputeZeroPointEnergy(Freqs);

  //macro grunseisen = sum_i(Cv_mode_i * micro grunseisen_i / Cv)
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
    macro_gruneisen  /= heat_capacity_v;



  //normalizing the vibrational energy and entropy if frequencies beyond the gamma point are used.
  if(Params::Parameters().GetSuperCellTotalKPoint()){
    //printf("\n\n Total k point = %i\n\n",Params::Parameters().GetSuperCellTotalKPoint());
    vib_energy /= double(Params::Parameters().GetSuperCellTotalKPoint());
    vib_entropy /= double(Params::Parameters().GetSuperCellTotalKPoint());
    heat_capacity_v /= double(Params::Parameters().GetSuperCellTotalKPoint());
  }


  


  printf("total Helmholtz vibrational energy (at %.0f K)  = %8.3f kJ/mol\n",
	 temperature,vib_energy);
  printf("total entropy (at %.0f K) = %8.3f J/(mol*K)\n",
	 temperature,vib_entropy);

  Energy_Vib = vib_energy/2625.5;
  Entropy = vib_entropy/1000.0/2625.5;
  Cv =  heat_capacity_v;
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
    gruneisen_parameter = macro_gruneisen;

}

void Cluster::PrintFrequencies(){

  int Nonzero = Freqs.Nonzero();

  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
    printf("QuasiHarmonic Vibrational frequencies (cm-1): %d real frequencies\n",Nonzero);
  else
   printf("Harmonic Vibrational frequencies (cm-1): %d real frequencies\n",Nonzero);   
  for (int i=0;i<Nonzero/6+1;i++) {
    for (int j=0;j<6;j++) {
      if ((6*i+j) < Nonzero && Freqs[6*i+j]>0.0) {
        printf("  %8.2f",Freqs[6*i+j]);
      }
      else
        break;              
    }
    printf("\n");
  }

}

void Cluster::ReadHessianData() {

cout << "ReadHessianData" << endl;

	int pointer[NMon+1];
	int natoms = 0;
        for(int imon=1;imon<=NMon;imon++) {
                pointer[imon] = 3*natoms;
		natoms += GetNumberOfAtoms(imon);
	}

        int DOFtotal = 3*natoms;
        Matrix Sum_HessQM(DOFtotal,DOFtotal);
	Matrix Sum_HessMM(DOFtotal,DOFtotal);
        int k = 0;
	int l = 0;
        for(int imon=1;imon<=NMon;imon++) {
		Monomers[imon].ReadHess();
                int DOF = 3*Monomers[imon].GetNumberOfAtoms();
                for(int i=0;i<DOF;i++) {
                        for(int j=0;j<DOF;j++) {
				l = j+(imon-1)*DOF;
//cout << "k, l = " << k << " " << l << endl;
                                Sum_HessQM(k,l) = Sum_HessQM(k,l) + Monomers[imon].GetHessQM_value(i,j);
				Sum_HessMM(k,l) = Sum_HessMM(k,l) + Monomers[imon].GetHessMM_value(i,j);
                        }
                        k += 1;
                }
        }
      	Sum_HessQM.PrintHessian("Sum_HessQM of 1-body");
	Sum_HessMM.PrintHessian("Sum_HessMM of 1-body");

	if(Params::Parameters().UseSpaceSymmetry()) {
		cout << "Using Symmetry to generate Hessian -- Watit" << endl;
                int A = 1;
                int B = 1;
                for(int idim=1;idim<=NDim;idim++) {
                        if(Dimers[idim].GetSymmetryFactor()!=0) {
                                if(Params::Parameters().DoLocal2BodyTruncation()&&(Dimers[idim].GetDimerSeparation()>Params::Parameters().GetLocalCutoff(0))) {
                                        break;
                                } else {
                                        Dimers[idim].ReadHess();
                                }
                                int ListSize = (Dimers[idim].GetSymmetryList().size() + Dimers[idim].GetPeriodicSymmetryList().size())/2;
                                for(int L=0; L<ListSize;L++){
                                        if(L<Dimers[idim].GetSymmetryList().size()/2){
                                                int dimA = Dimers[idim].GetSymmetryList()[2*L];
                                                int dimB = Dimers[idim].GetSymmetryList()[2*L+1];
                                                int NmonA = Monomers[dimA].GetNumberOfAtoms();
                                                int NmonB = Monomers[dimB].GetNumberOfAtoms();
                                                int DOFA = 3*NmonA;
                                                int DOFB = 3*NmonB;
                                                int c=0;
                                                Matrix RotT = Dimers[idim].GetRotationList()[L];
cout << "Dimer " << dimA << "," << dimB << endl;
RotT.Print("RotT");
                                                Matrix Rot = RotT;
                                                Rot.Transpose();
						Matrix RotHess(DOFA+DOFB,DOFA+DOFB);
						Matrix TwoDRotHess(DOFA+DOFB,DOFA+DOFB);
						for(int i=0;i<NmonA;i++) {
							for(int j=0;j<NmonB;j++) {
								
							}
							Matrix TempHessQMA(DOFA+DOFB,DOFA+DOFB);
							Matrix TempHessMMA(DOFA+DOFB,DOFA+DOFB);
							for(int j=0;j<DOFA;j++) {
								for(int m=0;m<DOFA;m++) {
									TempHessQMA(j,m) = Dimers[idim].GetHessQM_value(j,m);
									TempHessMMA(j,m) = Dimers[idim].GetHessMM_value(j,m);
								}
							}
							TempHessQMA.PrintHessian("Dimer HessQMA");
							TempHessMMA.PrintHessian("Dimer HessMMA");
						}
					} else {
                                                int K = L - Dimers[idim].GetSymmetryList().size()/2;
                                                int dimA = Dimers[idim].GetPeriodicSymmetryList()[2*K];
                                                int dimB = Dimers[idim].GetPeriodicSymmetryList()[2*K+1];
						int dimref = Dimers[idim].GetReferenceMonomerIndex();
                                                int NmonA = Monomers[dimA].GetNumberOfAtoms();
                                                int NmonB = Monomers[dimB].GetNumberOfAtoms();
                                                int DOFA = 3*NmonA;
                                                int DOFB = 3*NmonB;
                                                int c=0;
                                                Matrix RotT = Dimers[idim].GetPeriodicRotationList()[K];
cout << "Dimer " << dimA << "," << dimB << endl;
RotT.Print("RotT");
                                                Matrix Rot = RotT;
                                                Rot.Transpose();
                                                Matrix RotHess(DOFA+DOFB,DOFA+DOFB);
                                                Matrix TwoDRotHess(DOFA+DOFB,DOFA+DOFB);
                                                for(int i=0;i<NmonA;i++) {
                                                        Matrix TempHessQMA(DOFA+DOFB,DOFA+DOFB);
							Matrix TempHessMMA(DOFA+DOFB,DOFA+DOFB);
                                                        for(int j=0;j<DOFA;j++) {
                                                                for(int m=0;m<DOFA;m++) {
                                                                        TempHessQMA(j,m) = Dimers[idim].GetHessQM_value(j,m);
									TempHessMMA(j,m) = Dimers[idim].GetHessMM_value(j,m);
                                                                }
                                                        }
                                                TempHessQMA.PrintHessian("Dimer HessQMA");
						TempHessMMA.PrintHessian("Dimer HessMMA");
                                                }
					}
				}
			}
		}
	}

}

//***************************************************************************************************************************
// Watit - IR & Raman Intensities
//
// General Scheme Raman Intensity
// 	You need Unmassweight Hessian and Polarizability Derivative
// 	Unmassweight Hessian (RMode) is pulled from Cluster::ComputeHarmonicFrequencies
// 	Polarizability Derivative (PolD) is built in Cluster::ReadIntensityData using body interaction many (BIM) scheme
// 	Once you have both matrix, you can get Raman intensity for each normal mode by multiply each mode with PolD
// 	Then you follow QCHEM scheme, get symmetry and asymmetric tensor, and convert them to be Raman Intensity in A^4/AMU
// PolD Structure
//	Normally, PolD is a projection of dipole moment in electric field which is 3*3 matrix at each degree of freedom
//	Atom1	DOF1 = 	XX XY XZ
//	       		YX YY YZ
//	       		ZX ZY ZZ
//		DOF2 =  XX XY XZ
//			YX YY YZ
//			ZX ZY ZZ
//		DOF3 =  XX XY XZ
//			YX YY YZ
//			ZX ZY ZZ
//	Atom2	DOF1 =  XX XY XZ
//			...	
//	Computationally, we assume that the system is isotropic. the XY = YX, XZ = ZX and YZ = ZY
//	Then the PolD is print out as total degree of freedom by 6 (total_DOF * 6 element) in QChem
// 	where the 6 columns are XX, XY, YY, XZ, YZ, and ZZ respectively
// 	Atom1   DOF1 =  XX XY YY XZ YZ ZZ
// 		DOF2 =  XX XY YY XZ YZ ZZ
// 		DOF3 =  XX XY YY XZ YZ ZZ
// 	Atom2   DOF1 =  XX XY YY XZ YZ ZZ
// 			...
// BIM Scheme
// 	total PolD = Sum PolD_monomer + Sum (PolD_dimer - PolD_monomer) + 0.5*Sum (PolD_image_dimer - PolD_monomer)
//	Look at this row of matrix: mA = dimer A, dAB = dimer AB, idAB = image dimer AB
//	| total_PolD | | Sum PolD_monomer | |  Sum (PolD_dimer - PolD_monomer) | |   Sum (PolD_image_dimer - PolD_monomer) |
//	|   total    | |      	 mA 	  | | + dAB - mA + dAC - mA 	       | | + idAB - mA + dAC - mA 	     + ... |
//	|   total =  | |       	 mB       | | + dAB - mB	    + dBC - mB | | + idAB - mA            + dBC - mB + ... |
//	|   total    | |   	 mC       | |            + dAC - mC + dBC - mC | |             + dAC - mC + dBC - mC + ... |
//	PolD_dimer are collected as PolDA and PolDB due to incompretent of the coder, maybe need some clean up eventually
//	To match out dimer and monomer, you need a pointer to match them correctly
//	Note that, in PBC, we reduce redundant image dimers by half since it will give repeated result
//	It means that the real numbers of image_dimer must be recover by multily with SymFac, which is 2 
// Symmetry Scheme
// ---- Basic symmerty of PolD ---- 
// 	The operation undergoes 3D rotation
// 	For dAB to dAC as following
// 	PolD_AC = Rot.|Rot.PolDAB.Rot^t|
// 	1D rotation rotates PolDAB to local
// 	2D rotation rotates local to PolDAC
// 	3D rotation rotates DOF of a single atom of AB to AC
//	As mention in "PolD structure", PolDAB, first, need to be rebuilt as 3*3 at each DOF
//	Then operate with 3*3 rotation until getting 2D rotation matrix
//	Second, reduce 3*3 to 1*6 for each DOF which is 3*6 for each atom
//	Lastly, do final rotation and now PolD_AC is generated
// ---- Overview Scheme ----
//	Symmetry operation here follows Yoni scheme with some shortcut
//	Sum PolD_monomer is collect withouted symmetry because all monomer jobs are available to use right away
//	Sum PolD_dimer has two types including non-PBC and periodic dimers that equivalent to dimer in unit cell
//	Sum PolD_image_dimer is from the rest of dimers in PBC
//	Note that without symmetry periodic dimers is included in PolD_image_dimer
//
//	All dimers and image dimers are labeled by their corresponing symmetry factor (SymFac)
//		Only non-zero SymFac jobs will be created
//	 	The zero SymFac means that that dimer is symmetrical equivalent to one of unique dimers
//	For each unique dimer, there are list of dimers and rotational matrixs that symmetrically equivalent from this unique dimer
//	PolD of uncreated job will be created from this list by using rotation matrix according to "Basic symmetry of PolD"
//	Because rotating molecule sometimes results in change atomic position especially for the molecule contain inversion
//		Like C1 O2 O3 becomes C1 O3 O2
//	To correct this cause, we use GetAtomEquivalency function to switch order of 2D PolD before final rotation 
//
//***************************************************************************************************************************
void Cluster::ReadIntensityData() {

	// First, we make Polarizability Derivative Matrix
	// Sum Polarizability Derivative from Monomer
	int pointer[NMon+1];
	int natoms = 0;
//	cout << "pointer" << endl;
	for(int imon=1;imon<=NMon;imon++) {
		Monomers[imon].ReadIntensity();
		pointer[imon] = 3*natoms;
//		cout << pointer[imon] << endl;
		natoms += GetNumberOfAtoms(imon);
	}

	int DOFtotal = 3*natoms;
	Matrix Sum_DipD(DOFtotal,3);
	Matrix Sum_PolD(DOFtotal,6);
	int k = 0;
	for(int imon=1;imon<=NMon;imon++) {
		int DOF = 3*Monomers[imon].GetNumberOfAtoms();
		for(int i=0;i<DOF;i++) {
			for(int j=0;j<6;j++) {
				if(j<3) {
					Sum_DipD(k,j) = Sum_DipD(k,j) + Monomers[imon].GetDipD_value(i,j);
				}
				Sum_PolD(k,j) = Sum_PolD(k,j) + Monomers[imon].GetPolD_value(i,j);
			}
			k += 1;
		}
	}
//	cout << "check1" << endl;	
//	Sum_DipD.Print("Sum_DipD of 1-body");
//	Sum_PolD.Print("Sum_PolD of 1-body");
//	cout << "check2" << endl;
//                cout << "NDim " << NDim << endl;
//		cout << "NDim_images " << NDim_images << endl;
        if(Params::Parameters().UseSpaceSymmetry()) {
        	cout << "Using Symmetry to generate Vibrational intensities" << endl;

        	int A = 1;
        	int B = 1;
//		cout << "non-PBC" << endl;
        	for(int idim=1;idim<=NDim;idim++) {
                	//Dimers[idim].ReadIntensity();
                        if(Dimers[idim].GetSymmetryFactor()!=0) {
				if(Params::Parameters().DoLocal2BodyTruncation()&&(Dimers[idim].GetDimerSeparation()>Params::Parameters().GetLocalCutoff(0))) {
					break;
				} else {
					Dimers[idim].ReadIntensity();
				}
//				cout << "SymFac " << Dimers[idim].GetSymmetryFactor() << endl;
//				cout << "PeriodicSymFac " << Dimers[idim].GetPeriodicSymmetryFactor() << endl; 
				int ListSize = (Dimers[idim].GetSymmetryList().size() + Dimers[idim].GetPeriodicSymmetryList().size())/2;
                                for(int L=0; L<ListSize;L++){
                                        if(L<Dimers[idim].GetSymmetryList().size()/2){
						int dimA = Dimers[idim].GetSymmetryList()[2*L];
						int dimB = Dimers[idim].GetSymmetryList()[2*L+1];
						int NmonA = Monomers[dimA].GetNumberOfAtoms();
						int NmonB = Monomers[dimB].GetNumberOfAtoms();
			                        int DOFA = 3*NmonA;
               				        int DOFB = 3*NmonB;
						int c=0;
                                                Matrix RotT = Dimers[idim].GetRotationList()[L]; // Because Rot is prior transpose for hessian??
                                                Matrix Rot = RotT;
                                                Rot.Transpose();
						Matrix RotDipD(DOFA+DOFB,3);
						Matrix RotPolD(DOFA+DOFB,6);
						Matrix TwoDRotPolD(DOFA+DOFB,6);
						Matrix ThreeDRotPolD(3,6);
						for(int i=0;i<NmonA;i++) {
							Matrix DipDA_AtEachDOF(3,3);
							Matrix TempDipDA(3,3);
							Matrix PolDA_AtEachDOF(3,3);
							Matrix TempPolDA(3,3);
							for(int j=0;j<3;j++) {
								int k=j+i*3;
								for(int m=0;m<3;m++) {
									DipDA_AtEachDOF(j,m) = Dimers[idim].GetDipDA_value(k,m);
								}
							}
							TempDipDA = Rot.Multiply(DipDA_AtEachDOF).Multiply(RotT);
							for(int j=0;j<3;j++) {
								int k=j+i*3;
								int l=j+3*(Dimers[idim].GetAtomEquivalency()[L][c]-1);
								for(int m=0;m<3;m++) {
									RotDipD(l,m) = TempDipDA(j,m);
								}
								PolDA_AtEachDOF(0,0) = Dimers[idim].GetPolDA_value(k,0);
        	                                        	PolDA_AtEachDOF(0,1) = Dimers[idim].GetPolDA_value(k,1);
	        	                                        PolDA_AtEachDOF(0,2) = Dimers[idim].GetPolDA_value(k,3);
        	        	                               	PolDA_AtEachDOF(1,0) = Dimers[idim].GetPolDA_value(k,1);
                	        	                       	PolDA_AtEachDOF(1,1) = Dimers[idim].GetPolDA_value(k,2);
                        	        	               	PolDA_AtEachDOF(1,2) = Dimers[idim].GetPolDA_value(k,4);
                                	        	      	PolDA_AtEachDOF(2,0) = Dimers[idim].GetPolDA_value(k,3);
                         		               	       	PolDA_AtEachDOF(2,1) = Dimers[idim].GetPolDA_value(k,4);
                                        	        	PolDA_AtEachDOF(2,2) = Dimers[idim].GetPolDA_value(k,5);
								TempPolDA = Rot.Multiply(PolDA_AtEachDOF).Multiply(RotT);
								TwoDRotPolD(l,0) = TempPolDA(0,0);
								TwoDRotPolD(l,1) = (TempPolDA(0,1)+TempPolDA(1,0))/2;
								TwoDRotPolD(l,2) = TempPolDA(1,1);
								TwoDRotPolD(l,3) = (TempPolDA(0,2)+TempPolDA(2,0))/2;
								TwoDRotPolD(l,4) = (TempPolDA(1,2)+TempPolDA(2,1))/2;
								TwoDRotPolD(l,5) = TempPolDA(2,2);
							}
							c++;
						}
						
                                                for(int i=0;i<NmonB;i++) {
							Matrix DipDB_AtEachDOF(3,3);
							Matrix TempDipDB(3,3);
							Matrix PolDB_AtEachDOF(3,3);
                                                        Matrix TempPolDB(3,3);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                for(int m=0;m<3;m++) {
                                                                        DipDB_AtEachDOF(j,m) = Dimers[idim].GetDipDB_value(k,m);
                                                                }
                                                        }
                                                        TempDipDB = Rot.Multiply(DipDB_AtEachDOF).Multiply(RotT);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
								int l=j+3*(Dimers[idim].GetAtomEquivalency()[L][c]-1);
                                                                for(int m=0;m<3;m++) {
                                                                        RotDipD(l,m) = TempDipDB(j,m);
                                                                }
                                                                PolDB_AtEachDOF(0,0) = Dimers[idim].GetPolDB_value(k,0);
                                                                PolDB_AtEachDOF(0,1) = Dimers[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(0,2) = Dimers[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(1,0) = Dimers[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(1,1) = Dimers[idim].GetPolDB_value(k,2);
                                                                PolDB_AtEachDOF(1,2) = Dimers[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,0) = Dimers[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(2,1) = Dimers[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,2) = Dimers[idim].GetPolDB_value(k,5);
                                                                TempPolDB = Rot.Multiply(PolDB_AtEachDOF).Multiply(RotT);
								TwoDRotPolD(l,0) = TempPolDB(0,0);
                                                                TwoDRotPolD(l,1) = (TempPolDB(0,1)+TempPolDB(1,0))/2;
                                                                TwoDRotPolD(l,2) = TempPolDB(1,1);
                                                                TwoDRotPolD(l,3) = (TempPolDB(0,2)+TempPolDB(2,0))/2;
                                                                TwoDRotPolD(l,4) = (TempPolDB(1,2)+TempPolDB(2,1))/2;
                                                                TwoDRotPolD(l,5) = TempPolDB(2,2);
                                                        }
							c++;
						}
                                                for(int i=0;i<NmonA+NmonB;i++) {
                                                        Matrix TempPolD(3,6);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        TempPolD(j,k) = TwoDRotPolD(j+i*3,k);
                                                                }
                                                        }
                                                        ThreeDRotPolD = Rot.Multiply(TempPolD);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        RotPolD(j+i*3,k) = ThreeDRotPolD(j,k);
                                                                }
                                                        }
                                                }
						c=0;
						for(int i=0;i<DOFA;i++) {
							for(int j=0;j<6;j++) {
								A = pointer[dimA]+i;
                                                		if(j<3) {
                                                        		Sum_DipD(A,j) = Sum_DipD(A,j) + RotDipD(c,j)-Monomers[dimA].GetDipD_value(i,j);
                                                		}
								Sum_PolD(A,j) = Sum_PolD(A,j) + RotPolD(c,j)-Monomers[dimA].GetPolD_value(i,j);
							}
							c++;
						}
                                                for(int i=0;i<DOFB;i++) {
                                                        for(int j=0;j<6;j++) {
                                                                B = pointer[dimB]+i;
                                                                if(j<3) {
                                                                        Sum_DipD(B,j) = Sum_DipD(B,j) + RotDipD(c,j)-Monomers[dimB].GetDipD_value(i,j);
                                                                }
                                                                Sum_PolD(B,j) = Sum_PolD(B,j) + RotPolD(c,j)-Monomers[dimB].GetPolD_value(i,j);
                                                        }
							c++;
                                                }
                                        } else {
						//cout << "PeriodicSym" << endl;
						//cout << "PeriodicNSym" << Dimers[idim].GetPeriodicSymmetryList().size()/2<< endl;
						int K = L - Dimers[idim].GetSymmetryList().size()/2;
                                                int dimA = Dimers[idim].GetPeriodicSymmetryList()[2*K];
                                                int dimB = Dimers[idim].GetPeriodicSymmetryList()[2*K+1];
                                                int NmonA = Monomers[dimA].GetNumberOfAtoms();
                                                int NmonB = Monomers[dimB].GetNumberOfAtoms();
                                                int DOFA = 3*NmonA;
                                                int DOFB = 3*NmonB;
                                                int c=0;
                                                Matrix RotT = Dimers[idim].GetPeriodicRotationList()[K];
						Matrix Rot = RotT;
	       					Rot.Transpose();
						Matrix RotDipD(DOFA+DOFB,3);
                                                Matrix RotPolD(DOFA+DOFB,6);
                                                Matrix TwoDRotPolD(DOFA+DOFB,6);
                                                Matrix ThreeDRotPolD(3,6);
                                                for(int i=0;i<NmonA;i++) {
                                                        Matrix DipDA_AtEachDOF(3,3);
                                                        Matrix TempDipDA(3,3);
							Matrix PolDA_AtEachDOF(3,3);
                                                        Matrix TempPolDA(3,3);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                for(int m=0;m<3;m++) {
                                                                        DipDA_AtEachDOF(j,m) = Dimers[idim].GetDipDA_value(k,m);
                                                                }
                                                        }
                                                        TempDipDA = Rot.Multiply(DipDA_AtEachDOF).Multiply(RotT);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                int l=j+3*(Dimers[idim].GetPeriodicAtomEquivalency()[K][c]-1);
                                                                for(int m=0;m<3;m++) {
                                                                        RotDipD(l,m) = TempDipDA(j,m);
                                                                }
                                                                PolDA_AtEachDOF(0,0) = Dimers[idim].GetPolDA_value(k,0);
                                                                PolDA_AtEachDOF(0,1) = Dimers[idim].GetPolDA_value(k,1);
                                                                PolDA_AtEachDOF(0,2) = Dimers[idim].GetPolDA_value(k,3);
                                                                PolDA_AtEachDOF(1,0) = Dimers[idim].GetPolDA_value(k,1);
                                                                PolDA_AtEachDOF(1,1) = Dimers[idim].GetPolDA_value(k,2);
                                                                PolDA_AtEachDOF(1,2) = Dimers[idim].GetPolDA_value(k,4);
                                                                PolDA_AtEachDOF(2,0) = Dimers[idim].GetPolDA_value(k,3);
                                                                PolDA_AtEachDOF(2,1) = Dimers[idim].GetPolDA_value(k,4);
                                                                PolDA_AtEachDOF(2,2) = Dimers[idim].GetPolDA_value(k,5);
                                                                TempPolDA = Rot.Multiply(PolDA_AtEachDOF).Multiply(RotT);
                                                                TwoDRotPolD(l,0) = TempPolDA(0,0);
                                                                TwoDRotPolD(l,1) = (TempPolDA(0,1)+TempPolDA(1,0))/2;
                                                                TwoDRotPolD(l,2) = TempPolDA(1,1);
                                                                TwoDRotPolD(l,3) = (TempPolDA(0,2)+TempPolDA(2,0))/2;
                                                                TwoDRotPolD(l,4) = (TempPolDA(1,2)+TempPolDA(2,1))/2;
                                                                TwoDRotPolD(l,5) = TempPolDA(2,2);
                                                        }
                                                        c++;
                                                }

                                                for(int i=0;i<NmonB;i++) {
							Matrix DipDB_AtEachDOF(3,3);
							Matrix TempDipDB(3,3);
							Matrix PolDB_AtEachDOF(3,3);
                                                        Matrix TempPolDB(3,3);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                for(int m=0;m<3;m++) {
                                                                        DipDB_AtEachDOF(j,m) = Dimers[idim].GetDipDB_value(k,m);
                                                                }
                                                        }
                                                        TempDipDB = Rot.Multiply(DipDB_AtEachDOF).Multiply(RotT);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                int l=j+3*(Dimers[idim].GetPeriodicAtomEquivalency()[K][c]-1);
                                                                for(int m=0;m<3;m++) {
                                                                        RotDipD(l,m) = TempDipDB(j,m);
                                                                }
                                                                PolDB_AtEachDOF(0,0) = Dimers[idim].GetPolDB_value(k,0);
                                                                PolDB_AtEachDOF(0,1) = Dimers[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(0,2) = Dimers[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(1,0) = Dimers[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(1,1) = Dimers[idim].GetPolDB_value(k,2);
                                                                PolDB_AtEachDOF(1,2) = Dimers[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,0) = Dimers[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(2,1) = Dimers[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,2) = Dimers[idim].GetPolDB_value(k,5);
                                                                TempPolDB = Rot.Multiply(PolDB_AtEachDOF).Multiply(RotT);
                                                                TwoDRotPolD(l,0) = TempPolDB(0,0);
                                                                TwoDRotPolD(l,1) = (TempPolDB(0,1)+TempPolDB(1,0))/2;
                                                                TwoDRotPolD(l,2) = TempPolDB(1,1);
                                                                TwoDRotPolD(l,3) = (TempPolDB(0,2)+TempPolDB(2,0))/2;
                                                                TwoDRotPolD(l,4) = (TempPolDB(1,2)+TempPolDB(2,1))/2;
                                                                TwoDRotPolD(l,5) = TempPolDB(2,2);
                                                        }
                                                        c++;
                                                }
                                                for(int i=0;i<NmonA+NmonB;i++) {
                                                        Matrix TempPolD(3,6);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        TempPolD(j,k) = TwoDRotPolD(j+i*3,k);
                                                                }
                                                        }
                                                        ThreeDRotPolD = Rot.Multiply(TempPolD);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        RotPolD(j+i*3,k) = ThreeDRotPolD(j,k);
                                                                }
                                                        }
                                                }
                                                //RotPolD.Print("PeriodicRotPolD");
                                                c=0;
                                                for(int i=0;i<DOFA;i++) {
                                                        for(int j=0;j<6;j++) {
                                                                A = pointer[dimA]+i;
                                                                if(j<3) {
                                                                        Sum_DipD(A,j) = Sum_DipD(A,j) + 0.5*(RotDipD(c,j)-Monomers[dimA].GetDipD_value(i,j));
                                                                }
                                                                Sum_PolD(A,j) = Sum_PolD(A,j) + 0.5*(RotPolD(c,j)-Monomers[dimA].GetPolD_value(i,j));
                                                        }
                                                        c++;
                                                }
                                                for(int i=0;i<DOFB;i++) {
                                                        for(int j=0;j<6;j++) {
                                                                B = pointer[dimB]+i;
                                                                if(j<3) {
                                                                        Sum_DipD(B,j) = Sum_DipD(B,j) + 0.5*(RotDipD(c,j)-Monomers[dimB].GetDipD_value(i,j));
                                                                }
                                                                Sum_PolD(B,j) = Sum_PolD(B,j) + 0.5*(RotPolD(c,j)-Monomers[dimB].GetPolD_value(i,j));
                                                        }
                                                        c++;
                                                }
//                                                Sum_PolD.Print("PeriodicSum_PolD");
					}
                                }
                        }
                }
                if(Params::Parameters().IsPeriodic()) {
                        //cout << "ImagesPBC" <<endl;
                        for(int idim=1;idim<=NDim_images;idim++) {
//                                DimerImages[idim].ReadIntensity();
                                if(DimerImages[idim].GetSymmetryFactor()!=0){
					DimerImages[idim].ReadIntensity();
//					cout << "ImagesSymFac " << DimerImages[idim].GetSymmetryFactor() << endl;
                                        int ListSize = DimerImages[idim].GetSymmetryList().size()/2;
/*
					for(int i=0;i<DimerImages[idim].GetSymmetryList().size();i++){
						cout << DimerImages[idim].GetSymmetryList()[i] << endl;
					}
*/
//					cout << "ListSize " << ListSize << endl;
					for(int L=0; L<ListSize;L++){
                                                int dimA = DimerImages[idim].GetSymmetryList()[2*L];
                                                int dimB = DimerImages[idim].GetSymmetryList()[2*L+1];
                                                int NmonA = Monomers[dimA].GetNumberOfAtoms();
                                                int NmonB = Monomers[dimB].GetNumberOfAtoms();
                                                int DOFA = 3*NmonA;
                                                int DOFB = 3*NmonB;
                                                int c=0;
                                                Matrix RotT = DimerImages[idim].GetRotationList()[L];
                                                Matrix Rot = RotT;
                                                Rot.Transpose();
						Matrix RotDipD(DOFA+DOFB,3);
                                                Matrix RotPolD(DOFA+DOFB,6);
                                                Matrix TwoDRotPolD(DOFA+DOFB,6);
                                                Matrix ThreeDRotPolD(3,6);
                                                for(int i=0;i<NmonA;i++) {
                                                        Matrix DipDA_AtEachDOF(3,3);
                                                        Matrix TempDipDA(3,3);
							Matrix PolDA_AtEachDOF(3,3);		
                                                        Matrix TempPolDA(3,3);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                for(int m=0;m<3;m++) {
                                                                        DipDA_AtEachDOF(j,m) = DimerImages[idim].GetDipDA_value(k,m);
                                                                }
                                                        }
                                                        TempDipDA = Rot.Multiply(DipDA_AtEachDOF).Multiply(RotT);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                int l=j+3*(DimerImages[idim].GetAtomEquivalency()[L][c]-1);
                                                                for(int m=0;m<3;m++) {
                                                                        RotDipD(l,m) = TempDipDA(j,m);
                                                                }
                                                                PolDA_AtEachDOF(0,0) = DimerImages[idim].GetPolDA_value(k,0);
                                                                PolDA_AtEachDOF(0,1) = DimerImages[idim].GetPolDA_value(k,1);
                                                                PolDA_AtEachDOF(0,2) = DimerImages[idim].GetPolDA_value(k,3);
                                                                PolDA_AtEachDOF(1,0) = DimerImages[idim].GetPolDA_value(k,1);
                                                                PolDA_AtEachDOF(1,1) = DimerImages[idim].GetPolDA_value(k,2);
                                                                PolDA_AtEachDOF(1,2) = DimerImages[idim].GetPolDA_value(k,4);
                                                                PolDA_AtEachDOF(2,0) = DimerImages[idim].GetPolDA_value(k,3);
                                                                PolDA_AtEachDOF(2,1) = DimerImages[idim].GetPolDA_value(k,4);
                                                                PolDA_AtEachDOF(2,2) = DimerImages[idim].GetPolDA_value(k,5);
                                                                TempPolDA = Rot.Multiply(PolDA_AtEachDOF).Multiply(RotT);
                                                                TwoDRotPolD(l,0) = TempPolDA(0,0);
                                                                TwoDRotPolD(l,1) = (TempPolDA(0,1)+TempPolDA(1,0))/2;
                                                                TwoDRotPolD(l,2) = TempPolDA(1,1);
                                                                TwoDRotPolD(l,3) = (TempPolDA(0,2)+TempPolDA(2,0))/2;
                                                                TwoDRotPolD(l,4) = (TempPolDA(1,2)+TempPolDA(2,1))/2;
                                                                TwoDRotPolD(l,5) = TempPolDA(2,2);
                                                        }
                                                        c++;
                                                }
                                                for(int i=0;i<NmonB;i++) {
							Matrix DipDB_AtEachDOF(3,3);
							Matrix TempDipDB(3,3);
							Matrix PolDB_AtEachDOF(3,3);
                                                        Matrix TempPolDB(3,3);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                for(int m=0;m<3;m++) {
                                                                        DipDB_AtEachDOF(j,m) = DimerImages[idim].GetDipDB_value(k,m);
                                                                }
                                                        }
                                                        TempDipDB = Rot.Multiply(DipDB_AtEachDOF).Multiply(RotT);
                                                        for(int j=0;j<3;j++) {
                                                                int k=j+i*3;
                                                                int l=j+3*(DimerImages[idim].GetAtomEquivalency()[L][c]-1);
                                                                for(int m=0;m<3;m++) {
                                                                        RotDipD(l,m) = TempDipDB(j,m);
                                                                }
                                                                PolDB_AtEachDOF(0,0) = DimerImages[idim].GetPolDB_value(k,0);
                                                                PolDB_AtEachDOF(0,1) = DimerImages[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(0,2) = DimerImages[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(1,0) = DimerImages[idim].GetPolDB_value(k,1);
                                                                PolDB_AtEachDOF(1,1) = DimerImages[idim].GetPolDB_value(k,2);
                                                                PolDB_AtEachDOF(1,2) = DimerImages[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,0) = DimerImages[idim].GetPolDB_value(k,3);
                                                                PolDB_AtEachDOF(2,1) = DimerImages[idim].GetPolDB_value(k,4);
                                                                PolDB_AtEachDOF(2,2) = DimerImages[idim].GetPolDB_value(k,5);
                                                                TempPolDB = Rot.Multiply(PolDB_AtEachDOF).Multiply(RotT);
                                                                TwoDRotPolD(l,0) = TempPolDB(0,0);
                                                                TwoDRotPolD(l,1) = (TempPolDB(0,1)+TempPolDB(1,0))/2;
                                                                TwoDRotPolD(l,2) = TempPolDB(1,1);
                                                                TwoDRotPolD(l,3) = (TempPolDB(0,2)+TempPolDB(2,0))/2;
                                                                TwoDRotPolD(l,4) = (TempPolDB(1,2)+TempPolDB(2,1))/2;
                                                                TwoDRotPolD(l,5) = TempPolDB(2,2);
                                                        }
                                                        c++;
                                                }
                                                for(int i=0;i<NmonA+NmonB;i++) {
                                                        Matrix TempPolD(3,6);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        TempPolD(j,k) = TwoDRotPolD(j+i*3,k);
                                                                }
                                                        }
                                                        ThreeDRotPolD = Rot.Multiply(TempPolD);
                                                        for(int j=0;j<3;j++) {
                                                                for(int k=0;k<6;k++){
                                                                        RotPolD(j+i*3,k) = ThreeDRotPolD(j,k);
                                                                }
                                                        }
                                                }
                                                //RotPolD.Print("RotPolD");
                                                c=0;
                                                for(int i=0;i<DOFA;i++) {
                                                        for(int j=0;j<6;j++) {
                                                                A = pointer[dimA]+i;
                                                                if(j<3) {
                                                                        Sum_DipD(A,j) = Sum_DipD(A,j) + 0.5*(RotDipD(c,j)-Monomers[dimA].GetDipD_value(i,j));
                                                                }
                                                                Sum_PolD(A,j) = Sum_PolD(A,j) + 0.5*(RotPolD(c,j)-Monomers[dimA].GetPolD_value(i,j));
                                                        }
                                                        c++;
                                                }
                                                for(int i=0;i<DOFB;i++) {
                                                        for(int j=0;j<6;j++) {
                                                                B = pointer[dimB]+i;
                                                                if(j<3) {
                                                                        Sum_DipD(B,j) = Sum_DipD(B,j) + 0.5*(RotDipD(c,j)-Monomers[dimB].GetDipD_value(i,j));
                                                                }
                                                                Sum_PolD(B,j) = Sum_PolD(B,j) + 0.5*(RotPolD(c,j)-Monomers[dimB].GetPolD_value(i,j));
                                                        }
                                                        c++;
                                                }
					}
                                }
                        }
                }
	} else {//No symmetry
	
		int A = 1;
		int B = 1;
//		cout << "non-PBC" <<endl;
		for(int idim=1;idim<=NDim;idim++) {
			if(Dimers[idim].GetSymmetryFactor()!=0){
                                if(Params::Parameters().DoLocal2BodyTruncation()&&(Dimers[idim].GetDimerSeparation()>Params::Parameters().GetLocalCutoff(0))) {
                                        break;
                                } else {
                                        Dimers[idim].ReadIntensity();
                                }

				int dimA = Dimers[idim].GetIndexA();
				int dimB = Dimers[idim].GetIndexB();
				//int dimref = Dimers[idim].GetIndexRef();
				int DOFA = 3*Monomers[dimA].GetNumberOfAtoms();
				int DOFB = 3*Monomers[dimB].GetNumberOfAtoms();
				//int DOFref = 3*Monomers[dimref].GetNumberOfAtoms();
				//cout << "dimA = " << dimA << " dimB = " <<dimB << endl; //" dimref = " << dimref << endl;
				for(int i=0;i<DOFA;i++) {
					for(int j=0;j<6;j++) {
						A = pointer[dimA]+i;
						if(j<3) {
							Sum_DipD(A,j) = Sum_DipD(A,j) + Dimers[idim].GetDipDA_value(i,j) - Monomers[dimA].GetDipD_value(i,j);
						}
						Sum_PolD(A,j) = Sum_PolD(A,j) + Dimers[idim].GetPolDA_value(i,j) - Monomers[dimA].GetPolD_value(i,j);
					}
				}
				for(int i=0;i<DOFB;i++) {
					for(int j=0;j<6;j++) {
						B = pointer[dimB]+i;
						if(j<3) {
							Sum_DipD(B,j) = Sum_DipD(B,j) + Dimers[idim].GetDipDB_value(i,j) - Monomers[dimB].GetDipD_value(i,j);
						}
						Sum_PolD(B,j) = Sum_PolD(B,j) + Dimers[idim].GetPolDB_value(i,j) - Monomers[dimB].GetPolD_value(i,j);
					}
				}
			}
		}
		if(Params::Parameters().IsPeriodic()) {
//			cout << "PBC" <<endl;
			for(int idim=1;idim<=NDim_images;idim++) {
				if(DimerImages[idim].GetSymmetryFactor()!=0){
//					cout << "ImagesSymFac " << DimerImages[idim].GetSymmetryFactor() << endl;
					DimerImages[idim].ReadIntensity();
					int dimA = DimerImages[idim].GetIndexA();
					int dimB = DimerImages[idim].GetReferenceMonomerIndex();
					int DOFA = 3*Monomers[dimA].GetNumberOfAtoms();
					int DOFB = 3*Monomers[dimB].GetNumberOfAtoms();
					int SymFac = DimerImages[idim].GetSymmetryFactor();
//                           		cout << "dimA = " << dimA << " dimB = " <<DimerImages[idim].GetIndexB() <<endl;

					for(int i=0;i<DOFA;i++) {
        	                       		for(int j=0;j<6;j++) {
                                        		A = pointer[dimA]+i;
							if(j<3) {
								Sum_DipD(A,j) = Sum_DipD(A,j) + 0.5*SymFac*(DimerImages[idim].GetDipDA_value(i,j) - Monomers[dimA].GetDipD_value(i,j));
							}
                       			                Sum_PolD(A,j) = Sum_PolD(A,j) + 0.5*SymFac*(DimerImages[idim].GetPolDA_value(i,j) - Monomers[dimA].GetPolD_value(i,j));
                               			}
                       			}
	                       		for(int i=0;i<DOFB;i++) {
        	                       		for(int j=0;j<6;j++) {
                                        		B = pointer[dimB]+i;
							if(j<3) {
								Sum_DipD(B,j) = Sum_DipD(B,j) + 0.5*SymFac*(DimerImages[idim].GetDipDB_value(i,j) - Monomers[dimB].GetDipD_value(i,j));
							}
                       	                		Sum_PolD(B,j) = Sum_PolD(B,j) + 0.5*SymFac*(DimerImages[idim].GetPolDB_value(i,j) - Monomers[dimB].GetPolD_value(i,j));
						}
        	                       	}
                       		}
			}	
		}
	}
//        Sum_PolD_test.Print("Sum_PolD test unitcell only");
//        cout << "check3" << endl;
//	Sum_DipD.Print("HMBI Dipole Derivative");
//	Sum_PolD.Print("HMBI Polarizability Derivative");

	// Second, we pull UnmassWeighted Normal Modes
	//Matrix RMode = GetRMode();
	//RMode.Print("Normal Modes");
	// Third, we compute VPolD
	Matrix VDipD(DOFtotal,3);
	Matrix VPolD(DOFtotal,6);
	Matrix IRIntensity(DOFtotal,1);
	Matrix RamanIntensity(DOFtotal,1);
	double Vpp;
	double Vpq;
	double RInt;
	for (int i=0;i<DOFtotal;i++) {
		for (int j=0;j<6;j++) {
			if(j<3) {
				for (int k=0;k<DOFtotal;k++) {
					VDipD(i,j) += Sum_DipD(k,j)*GetRMode_value(k,i);
				}
			}
			for (int k=0;k<DOFtotal;k++) {
				VPolD(i,j) += Sum_PolD(k,j)*GetRMode_value(k,i);
			}
		}
		IntIR[i] = (VDipD(i,0)*VDipD(i,0)+VDipD(i,1)*VDipD(i,1)+VDipD(i,2)*VDipD(i,2))*975.12; // KM/MOL ??
		Vpp = (VPolD(i,0) + VPolD(i,2) + VPolD(i,5))/3;
		Vpq = 0.5*((VPolD(i,0)-VPolD(i,2))*(VPolD(i,0)-VPolD(i,2))+(VPolD(i,0)-VPolD(i,5))*(VPolD(i,0)-VPolD(i,5))+(VPolD(i,2)-VPolD(i,5))*(VPolD(i,2)-VPolD(i,5))+6*((VPolD(i,1)*VPolD(i,1))+VPolD(i,3)*VPolD(i,3)+VPolD(i,4)*VPolD(i,4)));
		RInt = 45*Vpp*Vpp + 7*Vpq;
		IntRaman[i] = RInt*0.078424; // A^4/AMU ??
	}


	//SetVPolD(VPolD);
	// GetVPolD_Matrix().Print("VPolD");

	// Fourth, we zero out small intensity
	int NIR;
	int NRaman;
	NIR = IntIR.GetLength();
	for(int i=0; i<NIR;i++) {
		if( fabs(IntIR[i]) < 1.0e-4) {
			IntIR[i] = 0.0;
		}
	}
	NRaman =  IntRaman.GetLength();
	for(int i=0; i<NRaman;i++) {
		if(fabs(IntRaman[i]) < 1.0e-4) {
			IntRaman[i] = 0.0;
		}
	}

	SetIntIR(IntIR);
	SetIntRaman(IntRaman);
	//IntIR.Print("IR Intensity");
}

// Watit
void Cluster::PlotIntensityData(int NVib) {

        double cvtf = 2.35482004503095; // 2sqrt(2LN2)
        double fwhm = 10; // Full Width at Half Maximum 
	double MaxFreq = Vib(NVib-1,0);
        double NormFac = 0.398942280401433; // 1/sqrt(2pi)
        double NSample = 2*fwhm*NVib;
        double PeakCutoff = 10;
        double sd=fwhm/cvtf;
	double IRMax = -9999999.0;
	double RamanMax = -9999999.0;

        string ir_output_file = "ir-plot.sh";
	FILE *ir_outfile;
        ir_outfile = fopen(ir_output_file.c_str(), "w");
        for(int i=0; i<NVib; i++) {
                if (IRMax < Vib(i,1)) { IRMax = Vib(i,1);}
        }
        fprintf(ir_outfile, "#!/usr/bin/gnuplot\n");
        fprintf(ir_outfile, "set termopt enhanced\n");
        fprintf(ir_outfile, "set style line 1 lw 1.5\n");
        fprintf(ir_outfile, "set font \"Helvetica,14\"\n");
        fprintf(ir_outfile, "set title \"Raman Spectroscopy\" font \"Bold-Helvetica,18\"\n");
        fprintf(ir_outfile, "set xlabel \"Frequency (cm^{-1})\"\n");
        fprintf(ir_outfile, "set ylabel \"Intensity (a.u.)\"\n");
        fprintf(ir_outfile, "set autoscale\n");
        fprintf(ir_outfile, "NormFac = %.6f\n",NormFac);
        fprintf(ir_outfile, "RamanMax = %.6f\n",IRMax);
        fprintf(ir_outfile, "sd = %.6f\n",sd);
        for(int i=0; i<NVib; i++) {
                fprintf(ir_outfile, "f%d(x)=100*%.6f*NormFac/RamanMax/sd*exp(-((x-%.6f)*(x-%.6f))/2/sd/sd)\n",i,Vib(i,1),Vib(i,0),Vib(i,0));
        }
        fprintf(ir_outfile, "f(x) = ");
        for(int i=0; i<NVib; i++) {
                if (i<(NVib-1)) {
                        fprintf(ir_outfile, "f%d(x)+",i);
                } else {
                        fprintf(ir_outfile, "f%d(x)\n",i);
                }
	}
        fprintf(ir_outfile, "set samples %.0f\n",NSample);
        fprintf(ir_outfile, "set term pngcairo enhanced size 800,600 font \"Helvetica,14\"\n");
        fprintf(ir_outfile, "set output \"ir-plot.png\"\n");
        fprintf(ir_outfile, "plot [x=0:%.0f]f(x) title \"Add title here\"\n",MaxFreq);
        fprintf(ir_outfile, "set table \"ir-plot.out\"\n");
        fprintf(ir_outfile, "replot\n");
        fprintf(ir_outfile, "unset table\n");
        fclose(ir_outfile);
        system("chmod +x ir-plot.sh");
//        system("./ir-plot.sh");

	if (Params::Parameters().DoRaman()) {
		string rm_output_file = "raman-plot.sh";
		FILE *rm_outfile;
		rm_outfile = fopen(rm_output_file.c_str(), "w");
		double MaxFreq = Vib(NVib-1,0);
		for(int i=0; i<NVib; i++) {
			if (RamanMax < Vib(i,2)) { RamanMax = Vib(i,2);}
		}
		fprintf(rm_outfile, "#!/usr/bin/gnuplot\n");
		fprintf(rm_outfile, "set termopt enhanced\n");
		fprintf(rm_outfile, "set style line 1 lw 1.5\n");
		fprintf(rm_outfile, "set font \"Helvetica,14\"\n");
		fprintf(rm_outfile, "set title \"Raman Spectroscopy\" font \"Bold-Helvetica,18\"\n");
		fprintf(rm_outfile, "set xlabel \"Frequency (cm^{-1})\"\n");
		fprintf(rm_outfile, "set ylabel \"Intensity (a.u.)\"\n");
		fprintf(rm_outfile, "set autoscale\n");
		fprintf(rm_outfile, "NormFac = %.6f\n",NormFac);
		fprintf(rm_outfile, "RamanMax = %.6f\n",RamanMax);
		fprintf(rm_outfile, "sd = %.6f\n",sd);
		for(int i=0; i<NVib; i++) {
			fprintf(rm_outfile, "f%d(x)=100*%.6f*NormFac/RamanMax/sd*exp(-((x-%.6f)*(x-%.6f))/2/sd/sd)\n",i,Vib(i,2),Vib(i,0),Vib(i,0));
		}
		fprintf(rm_outfile, "f(x) = ");
		for(int i=0; i<NVib; i++) {
			if (i<(NVib-1)) {
                        	fprintf(rm_outfile, "f%d(x)+",i);
			} else {
				fprintf(rm_outfile, "f%d(x)\n",i);
			}
		}
		fprintf(rm_outfile, "set samples %.0f\n",NSample);
		fprintf(rm_outfile, "set term pngcairo enhanced size 800,600 font \"Helvetica,14\"\n");
		fprintf(rm_outfile, "set output \"raman-plot.png\"\n");
		fprintf(rm_outfile, "plot [x=0:%.0f]f(x) title \"Add title here\"\n",MaxFreq);
		fprintf(rm_outfile, "set table \"raman-plot.out\"\n");
		fprintf(rm_outfile, "replot\n");
		fprintf(rm_outfile, "unset table\n");
		fclose(rm_outfile);
		system("chmod +x raman-plot.sh");
    	 	//system("./raman-plot.sh");

	}
}

// Returns the vibrational ZPE contribution, in kJ/mol
double Cluster::ComputeZeroPointEnergy(Vector freqs) {

  double zpe = 0.0;
  // ZPE = sum (h*c*nu_i), along with unit conversions to kJ/mol
  for (int i=0;i<freqs.GetLength();i++) {
    if (freqs[i] > 0)
      zpe += 0.5*h_planck*freqs[i]*c_light*100.0*Na/1000.0;
  }

  if(!Params::Parameters().GetSuperCellTotalKPoint())
    printf("Zero-point vibrational energy = %8.3f kJ/mol\n",zpe);
  else
    printf("Zero-point vibrational energy = %8.3f kJ/mol\n",zpe/double(Params::Parameters().GetSuperCellTotalKPoint()));

  return zpe;
}


//Print Vibrational Frequencies
void Cluster::PrintVibrationalFrequencies(Vector freq,string header,bool RemoveOldFile){

  //Create the frequency file if it does not exist
  FILE *freq_file;
  string infile = GetInputFilename();
  string freq_name = infile.substr(0,infile.size()-3);
  freq_name += ".freq";

  if(RemoveOldFile){
    if ((freq_file = fopen(freq_name.c_str(),"w"))==NULL){
      printf("Cluster::PrintVibrationalFrequencies : Cannot open file '%s'\n",freq_name.c_str());
      exit(1);
    }
  }
  else{
    if ((freq_file = fopen(freq_name.c_str(),"a"))==NULL){
      printf("Cluster::PrintVibrationalFrequencies : Cannot open file '%s'\n",freq_name.c_str());
      exit(1);
    }
  }
  //print frequencies
  fprintf(freq_file,"%s\n",header.c_str());
  for(int i=0;i<freq.GetLength();i++){
    //freqstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
    fprintf(freq_file,"    %10.6f\n",freq[i]);
  }
  fprintf(freq_file,"\n");

  fclose(freq_file);
}

/*

//read frequencies from file
void Cluster::ReadQuasiHarmonicFrequencies(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus){

  //Check to make sure all the sets of freqiencies vectors are the same length.
  //This does not check if the frequency in the frequency files of the right number of entries
  if(FreqPlus.GetLength() != reference_frequencies.GetLength() ){
    printf("Cluster::ReadQuasiHarmonicFrequencies : FreqPlus vector incorrect length \n  Length = %i instend of %i\n",
	   FreqPlus.GetLength(),reference_frequencies.GetLength());
    exit(1);
  }
  if(FreqMinus.GetLength() != reference_frequencies.GetLength()){
    printf("Cluster::ReadQuasiHarmonicFrequencies : FreqMinus vector incorrect length \n  Length = %i instend of %i \n",
	   FreqMinus.GetLength(),reference_frequencies.GetLength());
    exit(1);
  }

  //Get name frequencie file
  FILE *freq_file;
  string filename = GetInputFilename();
  string freq_name = filename.substr(0,filename.size()-3);
  freq_name += ".freq";
  
  //Open the frequency file 
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadQuasiHarmonicFrequencies : Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  }  

  //counts the number of geometry. Should be three
  int countgeom = 0;
  //read reference frequencies
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    //reference frequencies
    if(line.substr(0,20) == "ReferenceFrequencies"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discarding the "ReferenceFrequencies" string
      iss >> tmp;
      iss >> reference_volume;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<reference_frequencies.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> reference_frequencies[i];
      }
    }

    //plus frequencies
    if(line.substr(0,10) == "PlusVolume"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> VolumePlus;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<FreqPlus.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> FreqPlus[i];
      }
    }
  
    //plus frequencies
    if(line.substr(0,11) == "MinusVolume"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> VolumeMinus;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<FreqMinus.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> FreqMinus[i];
      }
    }
  }
  

  infile.close();

  if(countgeom != 3){
    printf("ERROR::Cluster::ReadQuasiHarmonicFrequencies() Incorrect number of volumes countgeom = %i\n",
	   countgeom);
    exit(1);
  }


  //for super cells make sure that all volumes have the same number of real freqencies
  int imag_num=0;
  int imag_plus=0;
  int imag_minus=0;
  for(int i = 0;i < reference_frequencies.GetLength(); i++)
    if(reference_frequencies[i] > 0.0)
       imag_num++;
  
  for(int i = 0; i < FreqPlus.GetLength(); i++)
    if(FreqPlus[i] > 0.0)
      imag_plus++;
   
  if(imag_num != imag_plus){
   printf("Cluster::ReadQuasiHarmonicFrequencies()::Error::Number of real frequencies differ between reference geometry and geometry of increased volume.\n");
   exit(0);
  }

  for(int i = 0; i < FreqMinus.GetLength(); i++)
    if(FreqMinus[i] > 0.0)
      imag_minus++;

  if(imag_num != imag_minus){
   printf("Cluster::ReadQuasiHarmonicFrequencies()::Error::Number of real frequencies differ between reference geometry and geometry of decreased volume.\n");
   exit(0);
  }

  //reference_frequencies.Print("reference frequencies");
  //FreqPlus.Print("FreqPlus");
  //FreqMinus.Print("FreqMinus");
  //exit(0);
}

*/

// GrunseisenParameter = -(dLn(freq)/dLn(V))
//void Cluster::SetGruneisenParameters(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus){
  

  /*
  //Check to make sure all the sets of freqiencies vectors are the same length.
  //This does not check if the frequency in the frequency files of the right number of entries
  if(FreqPlus.GetLength() != reference_frequencies.GetLength() ){
    printf("Cluster::SetGruneisenParameters() : FreqPlus vector incorrect length \n  Length = %i instend of %i\n",
	   FreqPlus.GetLength(),reference_frequencies.GetLength());
    exit(1);
  }
  if(FreqMinus.GetLength() != reference_frequencies.GetLength()){
    printf("Cluster::SetGruneisenParameters() : FreqMinus vector incorrect length \n  Length = %i instend of %i \n",
	   FreqMinus.GetLength(),reference_frequencies.GetLength());
    exit(1);
  }

  //check that volumes are non-zeros
  if(VolumePlus == 0 ||VolumeMinus == 0){
   printf("Cluster::SetGruneisenParameters() : Volumes have zero value VolumePlus = %f VolumeMinus = %f  \n",
	  VolumePlus,VolumeMinus);
    exit(1);
  }
  //check that there is a change in volume
  if(fabs(VolumePlus - VolumeMinus) < 0.000001){
   printf("Cluster::SetGruneisenParameters() : No change in volume Volume = %f  \n",
	  VolumePlus);
    exit(1);
  }
    

  //for(int i=0; i<3*GetTotalNumberOfAtoms(); i++)
  for(int i=0; i<Gruneisen_parameters.GetLength();i++)
    if(FreqPlus[i] < 0.00001 && FreqMinus[i] < 0.00001)
      Gruneisen_parameters[i] = 0.0; 
    else      Gruneisen_parameters[i] = log(FreqPlus[i]/FreqMinus[i]);

  
  
  Gruneisen_parameters.Scale(1/log(VolumeMinus/VolumePlus));
  Gruneisen_parameters.Print("Gruneisen_parameters");
  */

  //Quasiharmonic::quasiharmonic().SetIsotropicGruneisenParameters(FreqPlus,VolumePlus,FreqMinus,VolumeMinus);

  //gruneisen parameters have been initialized. Use them to find frequency
  //Gruneisen_init = 1;
//}


//void Cluster::InitializeQuasiHarmonicForSuperCell(){

  //int total_k_point = Params::Parameters().GetSuperCellTotalKPoint();

  //Gruneisen_parameters.Initialize(3*GetTotalNumberOfAtoms()*total_k_point);
  //reference_frequencies.Initialize(3*GetTotalNumberOfAtoms()*total_k_point);

  //printf("\nInitializing QuasiHarmonics for Supercell\n");
  //printf("total_k_point = %i\n",total_k_point);
  //printf("Natoms = %i\n",GetTotalNumberOfAtoms());
  //printf("Gruneisen_parameters length = %i\n",Gruneisen_parameters.GetLength());
  //printf("reference frequencies lenght = %i\n\n",reference_frequencies.GetLength());

//}
//Print Vibrational Analysis and normal modes in Molden Format into a new inpt_filename-molden.out file
void Cluster::PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index) {

  /* Create the xyz file */
  FILE *xyz;
  string infile = GetInputFilename();
  string molden_file = infile.substr(0,infile.size()-3);
  molden_file += "-molden.out";
  if ((xyz = fopen(molden_file.c_str(),"w"))==NULL) {
    printf("Cluster::PrintNormalModes : Cannot open file '%s'\n",
           molden_file.c_str());
    exit(1);
  }           

  ofstream moldenstream;
  moldenstream.open(molden_file.c_str());
  int Natoms = GetTotalNumberOfAtoms();
  //  moldenstream << "======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======\n";
  moldenstream << "[Molden Format]\n";
  moldenstream << "[Atoms] (Angs)\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << GetAtomicSymbol(i).c_str() << "  " <<  right << setw(5) << i+1 << "  "
         <<  right << setw(5) << GetAtomicNumber(i) << "  " 
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+1] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+2] << "\n";
  }

  moldenstream << "[FREQ]\n";
  for (int i=0; i < Nonzerofreqs; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  //zero frequency vibrations
  //if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++) {
      moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
    }
    //}

  moldenstream << "[FR-COORD]\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << GetAtomicSymbol(i).c_str() << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+1]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+2]*AngToBohr << "\n";
  }

  moldenstream << "[FR-NORM-COORD]\n";
  for (int i=0; i < Nonzerofreqs; i++)	{
    moldenstream << "vibration " << i+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  " 
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "            
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";            
    }
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++)  {
    moldenstream << "vibration " << i+Nimagfreqs-Nfreq+Nonzerofreqs+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  //if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++)  {          
      moldenstream << "vibration " << i+Nimagfreqs+1 << "\n";
      for (int k=0; k < Natoms; k++) {
        moldenstream <<  "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
	// }
    }
  }
//   moldenstream << "\n======= END OF MOLDEN-FORMATTED INPUT FILE =======\n";
  moldenstream.close();
}
//Print Vibrational Analysis and normal modes in Molden Format into a custom named file
void Cluster::PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index, string Name) {

  /* Create the xyz file */
  FILE *xyz;
  //string infile = GetInputFilename();
  //string molden_file = infile.substr(0,infile.size()-3);
  string molden_file = Name;
  molden_file += "-molden.out";
  if ((xyz = fopen(molden_file.c_str(),"w"))==NULL) {
    printf("Cluster::PrintNormalModes : Cannot open file '%s'\n",
           molden_file.c_str());
    exit(1);
  }           

  ofstream moldenstream;
  moldenstream.open(molden_file.c_str());
  int Natoms = GetTotalNumberOfAtoms();
  //  moldenstream << "======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======\n";
  moldenstream << "[Molden Format]\n";
  moldenstream << "[Atoms] (Angs)\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << GetAtomicSymbol(i).c_str() << "  " <<  right << setw(5) << i+1 << "  "
         <<  right << setw(5) << GetAtomicNumber(i) << "  " 
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+1] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+2] << "\n";
  }

  moldenstream << "[FREQ]\n";
  for (int i=0; i < Nonzerofreqs; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  //zero frequency vibrations
  //if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++) {
      moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
    }
    //}

  moldenstream << "[FR-COORD]\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << GetAtomicSymbol(i).c_str() << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+1]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << AtomicCoordinates[3*i+2]*AngToBohr << "\n";
  }

  moldenstream << "[FR-NORM-COORD]\n";
  for (int i=0; i < Nonzerofreqs; i++)	{
    moldenstream << "vibration " << i+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  " 
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "            
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";            
    }
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++)  {
    moldenstream << "vibration " << i+Nimagfreqs-Nfreq+Nonzerofreqs+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  //if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++)  {          
      moldenstream << "vibration " << i+Nimagfreqs+1 << "\n";
      for (int k=0; k < Natoms; k++) {
        moldenstream <<  "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
      }
    }
    //}
//   moldenstream << "\n======= END OF MOLDEN-FORMATTED INPUT FILE =======\n";
  moldenstream.close();
}

/*
void Cluster::PrintFrequencies(Vector freq;string filename){

  filename += .txt;
  File file;
  if (( = fopen(molden_file.c_str(),"w"))==NULL) {

  }


}
*/

// Compare original vector with sorted vector
Vector Cluster::CompareSortedOrigVecs(Vector orig, Vector sorted) {
  Vector original_index( orig.GetLength() );
  for (int m=0;m< original_index.GetLength(); m++) {
    original_index.Element(m) = -1;
  }

  for (int i=0;i<sorted.GetLength();i++) {
    for (int j=0;j<orig.GetLength();j++) {
      if  (sorted[i] == orig[j]) {
        //check if j-value not already stored
        int kcount =0; 
        for (int k=0;k<i;k++) {
          if (original_index[k] == j) {
            kcount++;
          }
        }
        if (kcount==0) {
          original_index[i] = j;
        }
      }            
    }               
  }                  
  return original_index;

}


Matrix Cluster::SortColumns(Matrix orig,Vector sorted) {

  if(orig.GetCols() != sorted.GetLength()){
    printf("Error::Matrix::SortColumns() orig column and sorted length mismatch orig column = %i, sort length = %i\n"
	   ,orig.GetCols(),sorted.GetLength());
  }


  Matrix SortedMatrix(orig.GetRows(),orig.GetCols());

  for(int i = 0;i < sorted.GetLength(); i++){
    int orig_col = int(sorted[i]);
    //printf("i = %i orig_col = %i\n",i,orig_col);
    SortedMatrix.SetColumnVector(orig.GetColumnVector(orig_col),i);
  }

  //sorted.Print("sorted");
  //orig.PrintHessian("orig");
  //SortedMatrix.PrintHessian("SortedMatrix");
  //exit(0);
  return SortedMatrix;

}

Matrix Cluster::ComputeInverseUnitCellMatrix() {
  
  Matrix Inv_unit_cell(3,3);
  /* this function computes the inverse of the 3x3 unit_cell matrix in units of bohr^-1*/
  /* this function is used for the formulation of the gradient in the
     presence of external pressure - KDN - aug 2011*/
  
  /* Requisites: unit_cell matrix */
  
  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];
  
  double alpha = UnitCellAngles[0];
  double beta = UnitCellAngles[1];       
  double gamma = UnitCellAngles[2];

  a *= AngToBohr;
  b *= AngToBohr;
  c *= AngToBohr;
  alpha *= DegreesToRadians;                 
  beta *= DegreesToRadians;                 
  gamma *= DegreesToRadians;

  //comment: inv_unit_cell is already in atomic units

  Inv_unit_cell.Initialize(3,3);
  Inv_unit_cell.Element(0,0) = 1/a;
  Inv_unit_cell.Element(0,1) = 0.0;
  Inv_unit_cell.Element(0,2) = 0.0;

  Inv_unit_cell.Element(1,0) = - cos(gamma) / sin(gamma) / a;
  Inv_unit_cell.Element(1,1) = 1 / sin(gamma) / b;            
  Inv_unit_cell.Element(1,2) = 0.0;    

  Inv_unit_cell.Element(2,0) = ( cos(alpha)*cos(gamma)/sin(gamma)-cos(beta)/sin(gamma) )
                        / a / sin(gamma)
                        / sqrt( 1-cos(beta)*cos(beta)
                               -1/sin(gamma)/sin(gamma)
                                *(cos(alpha)-cos(beta)*cos(gamma))
                               	*(cos(alpha)-cos(beta)*cos(gamma)) );

  Inv_unit_cell.Element(2,1) = ( cos(beta)*cos(gamma)/sin(gamma)-cos(alpha)/sin(gamma) )
                        / b / sin(gamma)
                        / sqrt( 1-cos(beta)*cos(beta)
                               -1/sin(gamma)/sin(gamma)               
                                *(cos(alpha)-cos(beta)*cos(gamma))   
                                *(cos(alpha)-cos(beta)*cos(gamma)) );

  Inv_unit_cell.Element(2,2) = 1 / c
                        / sqrt( 1-cos(beta)*cos(beta)
                               -1/sin(gamma)/sin(gamma)
                               	*(cos(alpha)-cos(beta)*cos(gamma))
                               	*(cos(alpha)-cos(beta)*cos(gamma)) );
  
  return Inv_unit_cell;

}

Matrix Cluster::FormExternalStressTensor() {
  /*             
    This routine will form the diagonal  external stress tensor using the external pressure
    provided by the user in atomic units (hartrees per bohr^3)
    This subroutine should be called only when the external pressure is provided and
    not the external stress tensor.
    In this case the stress tensor element = - pressure*kroneckardelta
  */

  Matrix ExternalStressTensor(3,3);

  double kronecker_delta = 0.0; 
  Matrix iden3(3,true);            
  ExternalStressTensor.Initialize(iden3);
  ExternalStressTensor.Scale(1*Params::Parameters().GetExternalPressure()*GigaPascalsToHartreesperBohrcube );

  //comment: ext. stress tensor is now in atomic units       
  return ExternalStressTensor;

}

Matrix Cluster::FormGradientExternalPressureTerm() {
  //GEPT stands for gradient's external pressure term
  Matrix GEPT(3,false);

  //now call FormExternalStressTensor() and ComputeInverseUnitCellMatrix().
  //this way the two matrices are formed at each opt step even though
  //externalstresstensor is constant throughout the opt.
  //check this during debug             

  //comment: GEPT is in atomic units

  //Matrix Inv_unit_cell = ComputeInverseUnitCellMatrix();
  GEPT = FormExternalStressTensor();
  //Inv_unit_cell.PrintHessian("inv_uc last");         
  GEPT.PrintHessian("EST last");   

  //Now multiplying the Pressure by the derivites of the volume with respect of the volume
  //Since the lattice vectors form a triangular matrix the determinant(which is also the volume) is the product of the diagonal entries.
  GEPT(0,0) *= unit_cell[1][1]*unit_cell[2][2]*AngToBohr*AngToBohr;
  GEPT(1,1) *= unit_cell[0][0]*unit_cell[2][2]*AngToBohr*AngToBohr;
  GEPT(2,2) *= unit_cell[0][0]*unit_cell[1][1]*AngToBohr*AngToBohr;

  //GEPT = ExternalStressTensor.Multiply(Inv_unit_cell,1);             
  GEPT.PrintHessian("gept");

  //GEPT.Transpose(); //look at the equation 3 in documentation for constant pressure opt
  //GEPT.PrintHessian("gept 2nd");

  //GEPT.Scale(UnitCellVolume()*AngToBohr*AngToBohr*AngToBohr ); //multiply the elements by the unit cell volume

  //GEPT.PrintHessian("gept last");
  return GEPT;

}

//Vibrational contribution to lattice vector gradients using quasiharmonic approximation
Matrix Cluster::ComputeQuasiHarmonicGradient(){

  /*
    gradient for lattice is found from using a chain rule derivative.

    dF/dv = dV/dv * sum_f( dF/df * df/dV )
  */
 

  double Temperature = Params::Parameters().GetTemperature();
  double h_kT = h_planck/(k_boltz*Temperature);
  double Vol = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  double Ref_Vol = Quasiharmonic::quasiharmonic().GetReferenceVolume()*AngToBohr*AngToBohr*AngToBohr;
  //double Ref_Vol = reference_volume*AngToBohr*AngToBohr*AngToBohr;

  //printf("h_kT = %E\n",h_kT);
  //printf("Vol = %f Bohr\n",Vol);
  //printf("Ref_Vol = %f Bohr\n",Ref_Vol);

  //derivative of vibrational contribuation with respect to volume
  double dF_dVol = 0;
  for(int i = 0; i<Freqs.GetLength() ; i++){
    if(Freqs[i] > 0.0 ){
    
      double dE = 0.0;//gradient contribution from frequency mode
      double freq_Hz = Freqs[i]*WaveNumberToFreq;
      double freq_ref_Hz =  Quasiharmonic::quasiharmonic().GetReferenceFrequencies()[i]*WaveNumberToFreq;
      //double freq_ref_Hz = reference_frequencies[i]*WaveNumberToFreq;
      double Grun = Quasiharmonic::quasiharmonic().GetGruneisenParameters()[i];
      //double Grun = Gruneisen_parameters[i];
      
      /*printf("Freqs[%i] = %f cm^-1\n",i,Freqs[i]);
      printf("Freqs[%i] = %E Hz\n",i,freq_Hz);
      printf("Ref_Freqs[%i] = %E cm^-1\n",i,Quasiharmonic::quasiharmonic().GetReferenceFrequencies()[i]);
      printf("Ref_Freqs[%i] = %E Hz\n",i,freq_ref_Hz);
      printf("Grun[%i] = %f\n",i,Grun);*/
      
      //derivative of Vibrational Energy with respect to frequency
      dE = 0.5;
      if(Temperature != 0)
        dE += 1/(exp(h_kT*freq_Hz) - 1);
      dE *= h_planck/HartreesToJ;

      //printf(" dF/df = %E\n",dE);
      
      //derivative of frequency with respect to volume using Quasiharmonic approximation
      dE *= -Grun*freq_ref_Hz*pow(Vol,-Grun-1)*pow(Ref_Vol,Grun);
      
      //printf(" df/dV = %E\n",-Grun*freq_ref_Hz*pow(Vol,-Grun-1)*pow(Ref_Vol,Grun));
      //printf(" dF[%i]/dV = %E\n",i,dE);
      
      dF_dVol += dE;
      //printf(" sum = %E\n\n",dF_dVol);
    }
  }


  printf("\nUnit cell vectors: in Angstroms\n");
  printf("A: (%f,%f,%f)\n",unit_cell[0][0],unit_cell[0][1],unit_cell[0][2]);
  printf("B: (%f,%f,%f)\n",unit_cell[1][0],unit_cell[1][1],unit_cell[1][2]);
  printf("C: (%f,%f,%f)\n",unit_cell[2][0],unit_cell[2][1],unit_cell[2][2]);

  //Including derivatives of the volume with respect to the lattice vectors
  Matrix dF_dv(3,3);
  dF_dv.Set();
  dF_dv(0,0) = unit_cell[1][1]*unit_cell[2][2]*AngToBohr*AngToBohr*dF_dVol;
  dF_dv(1,1) = unit_cell[0][0]*unit_cell[2][2]*AngToBohr*AngToBohr*dF_dVol;
  dF_dv(2,2) = unit_cell[0][0]*unit_cell[1][1]*AngToBohr*AngToBohr*dF_dVol;


   //normalizing gradient if frequencies beyond the gamma point were used.
   if(Params::Parameters().IsSupercellJob())
     dF_dv.Scale(1/double(Params::Parameters().GetSuperCellTotalKPoint()));


  /*printf("dV_dx1 = %f\n",unit_cell[1][1]*unit_cell[2][2]*AngToBohr*AngToBohr);
  printf("dV_dy2 = %f\n",unit_cell[0][0]*unit_cell[2][2]*AngToBohr*AngToBohr);
  printf("dV_dz3 = %f\n",unit_cell[0][0]*unit_cell[1][1]*AngToBohr*AngToBohr);*/

  //dF_dv.Print("dF_dv");
  
  return dF_dv;

}

// Compute Axilrod-Teller-Muto atomic 3-body dispersion - between all triplets
// of monomers.  This is the function used by the AIFF.
double Cluster::ComputeThreeBodyDispersion(string type) {
  time_t start_time, stop_time;
  start_time = time(NULL);

    double E3b_disp = 0.0;

  // Loop over triplets of monomers (trimers)
  for (int imon=1; imon<=NMon; imon++) {
    int NatomsI = GetNumberOfAtoms(imon);
    for (int jmon=imon+1; jmon<=NMon; jmon++) {
      int NatomsJ = GetNumberOfAtoms(jmon);
      for (int kmon=jmon+1; kmon<=NMon; kmon++) {
	int NatomsK = GetNumberOfAtoms(kmon);
	
	// Loop over individual atoms in the trimers
	for (int i=0;i<NatomsI;i++) {
	  Atom AtomI = Monomers[imon].GetAtom(i);
	  for (int j=0;j<NatomsJ;j++) {
	    Atom AtomJ = Monomers[jmon].GetAtom(j);
	    for (int k=0;k<NatomsK;k++) {
	      Atom AtomK = Monomers[kmon].GetAtom(k);
	      // Because we sum over only unique triplets, get rid of
	      // 1/6 factor.
	      E3b_disp += ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	    }
	  }
	}
      }
    }
  }

  printf("     Non-periodic contribution to 3-body dispersion = %.4f kJ/mol\n",
	 E3b_disp*HartreesToKJpermole);

  if ( Params::Parameters().IsPeriodic() ) {

    // Two cases:  1) 2 atoms in central cell, 1 in periodic image cell
    //             2) 1 atom in central cell, 2 in periodic image cell

    // Are my factors of 1/3 correct?
    // Logic: default: 1/6 sum(ijk) U_ijk (if we sum over *all* ijk, not
    // just unique triplets).  Now consider the periodic case:
    // If molecules 1 & 2 in unit cell, 3 is image, would get
    //       E <- (1/6) U_123 + (1/6) U_213 = (1/3) U_123
    // If molecule 1 in unit cell, 2 & 3 are images, would get
    //       E <- (1/6) U_123 + (1/6) U_132 = (1/3) U_123

    double case1 = 0.0, case2 = 0.0;

    double r_cutoff = Params::Parameters().GetMaxThreeBodyDispersionRadius();
    printf("     Creating list of periodic image monomers within radius %f\n",r_cutoff);
    CreatePeriodicImageMonomerList(r_cutoff); // create list of image monomers


    // Case 1: 2 atoms in central cell, 1 in periodic image cell
    printf("     Case 1: 2 monomers in unit cell, 1 image monomer\n"); fflush(stdout);
    for (int imon=1; imon<=NMon; imon++) {
      printf("       Monomers: %d, *, *\n",imon); fflush(stdout);
      int NatomsI = GetNumberOfAtoms(imon);
      for (int jmon=imon+1; jmon<=NMon; jmon++) {
	int NatomsJ = GetNumberOfAtoms(jmon);
	for (int kmon=0; kmon<=NMon_images; kmon++) {
	  //printf("Monomers: %d, %d, %d\n",imon,jmon,kmon);
	  int NatomsK = MonomerImages[kmon].GetNumberOfAtoms();

	  // Loop over individual atoms in the trimers
	  for (int i=0;i<NatomsI;i++) {
	    Atom AtomI = Monomers[imon].GetAtom(i);
	    for (int j=0;j<NatomsJ;j++) {
	      Atom AtomJ = Monomers[jmon].GetAtom(j);
	      for (int k=0;k<NatomsK;k++) {
		Atom AtomK = MonomerImages[kmon].GetAtom(k);
		
		double contrib = (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
		E3b_disp += contrib;
		case1 += contrib;
		//E3b_disp += (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	      }
	    }
	  }
	}
      }
    }

    printf("       Case 1 contribution = %.4f kJ/mol\n",case1*HartreesToKJpermole);

    // Case 2: 1 atom in central cell, 2 in periodic image cells
    printf("     Case 2: 1 monomer in unit cell, 2 image monomer\n"); fflush(stdout);
    for (int imon=1; imon<=NMon; imon++) {
      printf("       Monomers: %d, *, *\n",imon); fflush(stdout);
      int NatomsI = GetNumberOfAtoms(imon);
      for (int jmon=0; jmon<=NMon_images; jmon++) {
	int NatomsJ = MonomerImages[jmon].GetNumberOfAtoms();
	for (int kmon=jmon+1; kmon<=NMon_images; kmon++) {
	  int NatomsK = MonomerImages[kmon].GetNumberOfAtoms();

	  // Loop over individual atoms in the trimers
	  for (int i=0;i<NatomsI;i++) {
	    Atom AtomI = Monomers[imon].GetAtom(i);
	    for (int j=0;j<NatomsJ;j++) {
	      Atom AtomJ = MonomerImages[jmon].GetAtom(j);
	      for (int k=0;k<NatomsK;k++) {
		Atom AtomK = MonomerImages[kmon].GetAtom(k);

		double contrib = (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
		E3b_disp += contrib;
		case2 += contrib;
		//E3b_disp += (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	      }
	    }
	  }
	}
      }
    }
    printf("       Case 2 contribution = %.4f kJ/mol\n",case2*HartreesToKJpermole);
  }


  stop_time = time(NULL);
  double disp_time = difftime(stop_time, start_time);
  printf("     Time spent computing AIFF 3-body dispersion = %.0f seconds\n",
	 disp_time);
  return E3b_disp;
}


double Cluster::ComputeAxilrodTellerMutoDispersionTerm(Atom AtomI, Atom AtomJ, Atom AtomK, string type) {

  // Get geometrical parameters... in bohr and radians 
  double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
  double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
  double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;
  
  // use law of cosines to get cos(phi)
  double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
  double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
  double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;
  
  /*
    printf("Rij = %f, Rik = %f, Rjk = %f\n",Rij,Rik,Rjk);
    printf("PhiI = %f, PhiJ = %f, PhiK = %f,  Sum = %f\n",
	 acos(cosPhiI)*RadiansToDegrees,
	 acos(cosPhiJ)*RadiansToDegrees,
	 acos(cosPhiK)*RadiansToDegrees,
	 (acos(cosPhiI)+acos(cosPhiJ)+acos(cosPhiK))*RadiansToDegrees); 
  */
  
  // Get the C9 coefficient
  double C9ijk;
  if (type == "Tkatchenko") 
    C9ijk = AtomI.EstimateC9Coefficient(AtomJ, AtomK, "Tkatchenko"); // tabulated estimate
  else 
    C9ijk = AtomI.CasimirC9Coefficient(AtomJ, AtomK); // Casimir-Polder integration (preferred)

  // Get the van der Waals diameters
  double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
  double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
  double Dk = AtomK.LookupAtomicDispersionParameter("Rvdw");
  
  // Get Tang-Toennies-type 3-body damping factor
  // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
  double betaIJ = -0.31*(Di+Dj) + 3.43;
  double betaIK = -0.31*(Di+Dk) + 3.43;
  double betaJK = -0.31*(Dj+Dk) + 3.43;
  
  double F6ij = AtomI.TangToenniesDampingFactor(6,betaIJ,Rij);
  double F6ik = AtomI.TangToenniesDampingFactor(6,betaIK,Rik);
  double F6jk = AtomJ.TangToenniesDampingFactor(6,betaJK,Rjk);
  double damping = F6ij * F6ik * F6jk;
  
  // ATM 3-body dispersion contribution
  double Eatm = damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
    (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
  
  //printf("Rij = %f, Rik = %f, Rjk = %f   E3 = %17.12f kcal/mol\n",Rij,Rik,Rjk,Eatm*HartreesToKcalpermole);

  return Eatm;
}

void Cluster::ReadDispersionAtomTypes(ifstream& infile) {

  int Natoms = GetTotalNumberOfAtoms();
  DispersionAtomTypes = new string [Natoms];
  
  string line;
  // Rewind the file, just in case
  Rewind(infile);

  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,11)=="$dispersion") {
      int j=0;
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) == "$end") {
	  break;
	}
	else if ( !line.empty() ) {
	  istringstream iss(line);
	  string atom_num, type;
	  iss >> atom_num;
	  iss >> type;

	  // Set the atom_Type
	  DispersionAtomTypes[j] = type;
	  j++;
	}
      }
    }
  }
  infile.clear();

  //for (int i=0; i<Natoms;i++) {
  //  printf("Atom %d, type = |%s|\n",i,DispersionAtomTypes[i].c_str());    
  //}

  int iatom = 0;
  for (int imon=1;imon<=NMon;imon++) {
    for (int i=0;i<Monomers[imon].GetNumberOfAtoms();i++) {
      Monomers[imon].GetAtom(i).SetDispersionAtomType( DispersionAtomTypes[iatom] );
      iatom++;
    }
  }
  
}


void Cluster::ReadDispersionCoefficients(ifstream& infile) {

  int Natoms = GetTotalNumberOfAtoms();
  double* C6_UCHF = new double [Natoms];
  double* C6 = new double [Natoms];
  

  string line;
  // Rewind the file, just in case
  Rewind(infile);

  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,11)=="$dispersion") {
      int j=0;
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) == "$end") {
	  break;
	}
	else if ( !line.empty() ) {
	  istringstream iss(line);
	  string atom_num;
	  iss >> atom_num;
	  // Read the important data
	  iss >> C6_UCHF[j];
	  iss >> C6[j];
	  j++;
	}
      }
    }
  }
  infile.clear();

  printf("\nReading C6 dispersion coefficients\n");
  for (int i=0; i<Natoms;i++) {
    printf("Atom %d, C6(UCHF) = %f, new C6 = %f\n",i,C6_UCHF[i],C6[i]);    
  }
  printf("\n");


  int iatom = 0;
  for (int imon=1;imon<=NMon;imon++) {
    for (int i=0;i<Monomers[imon].GetNumberOfAtoms();i++) {
      Monomers[imon].GetAtom(i).SetUCHFDispersionCoefficients( C6_UCHF[iatom] );
      Monomers[imon].GetAtom(i).SetDispersionCoefficients( C6[iatom] );
      // Determine the dispersion atom type
      string atom_type;
      if (Monomers[imon].GetAtom(i).GetSymbol() == "H") atom_type = "Hs";
      else if (Monomers[imon].GetAtom(i).GetSymbol() == "C") atom_type = "Csp2";
      else if (Monomers[imon].GetAtom(i).GetSymbol() == "N") atom_type = "Nsp2sp3";
      else if (Monomers[imon].GetAtom(i).GetSymbol() == "O") atom_type = "Osp3";
      else {
	printf("Cluster::ReadDispersionCoefficients: Unknown Atom Type %s\n",(Monomers[imon].GetAtom(i).GetSymbol().c_str()));
	exit(1);
      }
      Monomers[imon].GetAtom(i).SetDispersionAtomType( atom_type );
      iatom++;
    }
  }

  

  delete [] C6;
  delete [] C6_UCHF;

}


double Cluster::ComputeEmpiricalTwoBodyDispersion() {

  double E_2body_disp = 0.0;

  for (int i=1;i<=NDim;i++) {
    E_2body_disp += Dimers[i].EstimateTwoBodyDispersion();
  }
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      E_2body_disp += 0.5*DimerImages[i].EstimateTwoBodyDispersion();
    }
  }

  printf("2-Body Dispersion (dimer routines): %f\n",
	 E_2body_disp*HartreesToKcalpermole);
  return E_2body_disp;
}


// Computes 2-body dispersion in a finite system (or purely within the
// unit cell)
double Cluster::ComputeTwoBodyDispersion() {

  time_t start_dispersion_time, stop_dispersion_time;
  start_dispersion_time = time(NULL);
  if ( ! Params::Parameters().IsPeriodic() ) {
    // if periodic, printing handled elsewhere
    printf("\nStep 2: Compute 2-body dispersion contribution\n"); 
  }

  double E_2body_disp = 0.0;
  //double E_2body_disp_C6 = 0.0;
  //double E_2body_disp_C8 = 0.0;
  //double E_2body_disp_C10 = 0.0;

  for (int i=1;i<=NDim;i++) {
    E_2body_disp += Dimers[i].ComputeTwoBodyDispersion();
    //E_2body_disp_C6 += Dimers[i].ComputeTwoBodyDispersionC6();
    //E_2body_disp_C8 += Dimers[i].ComputeTwoBodyDispersionC8();
    //E_2body_disp_C10 += Dimers[i].ComputeTwoBodyDispersionC10();
  }

  //printf("  2-Body Intermolecular Dispersion: %f kJ/mol, E_2body_disp_C6 = %f, E_2body_disp_C8 = %f,E_2body_disp_C10 = %f\n",
  //	 E_2body_disp*HartreesToKJpermole, E_2body_disp_C6*HartreesToKJpermole,E_2body_disp_C8*HartreesToKJpermole,E_2body_disp_C10*HartreesToKJpermole);


  stop_dispersion_time = time(NULL);
  if ( ! Params::Parameters().IsPeriodic() ) {
    // if periodic, we print out elsewhere
    printf("     Final Total dispersion  = %.4f kJ/mol\n",
	     E_2body_disp*HartreesToKJpermole);
    double dispersion_time = difftime(stop_dispersion_time, start_dispersion_time);
    printf("     Time spent computing AIFF 2-body dispersion = %.0f seconds\n\n", dispersion_time);
  }

  return E_2body_disp;
}

// Computes the two body dispersion interaction in an infinite lattice
double Cluster::ComputePeriodicTwoBodyDispersion() {

  time_t start_dispersion_time, stop_dispersion_time;
  start_dispersion_time = time(NULL);

  double E_2body_disp = 0.0;
  double E_2body_disp_C6 = 0.0;
  double E_2body_disp_C8 = 0.0;
  double E_2body_disp_C10 = 0.0;



  printf("\nStep 2: Compute 2-body dispersion contribution\n"); 
  // Step 1: Evaluate contribution within central unit cell
  printf(" - Step 2a: Evaluate contribution within central unit cell\n");
  E_2body_disp = ComputeTwoBodyDispersion();

  //for (int i=1;i<=NDim;i++) {
  //  E_2body_disp += Dimers[i].ComputeTwoBodyDispersion();
  //  E_2body_disp_C6 += Dimers[i].ComputeTwoBodyDispersionC6();
  //  E_2body_disp_C8 += Dimers[i].ComputeTwoBodyDispersionC8();
  //  E_2body_disp_C10 += Dimers[i].ComputeTwoBodyDispersionC10();
  //}

  // Step 2: Evaluate contribution between central unit cell monomers
  // and their periodic images.  Sum interactions up to a large cutoff
  // distance.
  

  // Identify how far we have to go along each unit cell direction to
  // stay within the cutoff.  Add 1 extra image cell to each, for good
  // measure.
  double E_2body_disp_pbc = 0.0;
  double E_2body_disp_pbc_C6 = 0.0;
  double E_2body_disp_pbc_C8 = 0.0;
  double E_2body_disp_pbc_C10 = 0.0;

  // Distance cutoff beyond which 2-body dispersion is neglected in AIFF
  double r_cutoff = Params::Parameters().GetMaxTwoBodyDispersionRadius();  // in Angstroms

  printf(" - Step 2b: Evaluate periodic 2-body dispersion interactions within %.2f Angstroms\n",r_cutoff);
  int Nv[3] = {1,1,1}; // start with 1 image cell in each direction
  for (int i=0;i<3;i++) {
      double dist = 0;
      while (dist  < r_cutoff) {
	dist +=  unit_cell[i].Norm();
	Nv[i] += 1;
      }
    }

  // Pass 1
  bool skip; // flag for ignoring some.
  int keepit = 0;
  // Loop over the image cells, in both positive & negative directions
  for (int x=-Nv[0];x<=Nv[0];x++) 
    for (int y=-Nv[1];y<=Nv[1];y++)
      for (int z=-Nv[2];z<=Nv[2];z++) {
	skip = false;

	if (x==0 && y==0 && z==0) {
	  skip = true; // Skip the central unit cell here.  It has
		       // been accounted for already.
	}

	if (!skip) {
	  // Create copy of Monomers that we can translate as needed;
	  Monomer* ImageMonomers = new Monomer[NMon+1];
	  for (int i=1;i<=NMon;i++) {
	    ImageMonomers[i] = Monomers[i];
	  }

	  // Determine the shift from the central cell to the image cell
	  Vector shift(3);
	  for (int i=0;i<3;i++) {
	    shift[i] = x*unit_cell[0][i] + y*unit_cell[1][i] + 
	      z*unit_cell[2][i];
	  }

	  // Now translate the monomers
	  for (int imon=1;imon<=NMon;imon++) {

	    bool IsThisMonomerLocal = false;

	    Vector new_com(3);
	    new_com = ImageMonomers[imon].GetCenterOfMass();
	    
	    // Add the shift and translate the monomer
	    new_com += shift;
	    ImageMonomers[imon].Translate(new_com);

	    // Now pair this image monomer with each monomer in central unit 
	    // cell and test distance relative to the cutoff
	    for (int jmon=1;jmon<=NMon;jmon++) {

	      Dimer Tmp;
	      //printf("jmon = %d, imon = %d\n",jmon,imon);
	      Tmp.Initialize(Monomers[jmon],ImageMonomers[imon]);
	      //printf("\ndimer from %d and image %d, distance = %f\n",jmon,imon,Tmp.GetDimerSeparation() );
	      //Tmp.PrintQChemCartesian();
	      if ( Tmp.GetDimerSeparation() < r_cutoff ) {
		if ( Params::Parameters().PrintLevel() > 3) {
		  printf("K_vec = (%d,%d,%d)\n",x,y,z);
		  printf("(%d*,%d) Separation = %f. Keeping it.\n",
		       imon,jmon,Tmp.GetDimerSeparation() );
		}
		E_2body_disp_pbc += Tmp.ComputeTwoBodyDispersion();
		//E_2body_disp_pbc_C6 += Tmp.ComputeTwoBodyDispersionC6();
		//E_2body_disp_pbc_C8 += Tmp.ComputeTwoBodyDispersionC8();
		//E_2body_disp_pbc_C10 += Tmp.ComputeTwoBodyDispersionC10();
		keepit++;
	      }
	    }
	  }
	  delete [] ImageMonomers;
	}
	
      }	  

  printf("     Total # of PBC dimers used in dispersion lattice sum = %d\n",keepit);
  printf("     Total central unit cell dispersion contribution = %.4f\n",
	 E_2body_disp*HartreesToKJpermole);
	//printf("     Total central unit cell dispersion contribution C6 = %.4f\n",
	 //E_2body_disp_C6*HartreesToKJpermole);
	 //printf("     Total central unit cell dispersion contribution c8 = %.4f\n",
	// E_2body_disp_C8*HartreesToKJpermole);
	 //printf("     Total central unit cell dispersion contribution c10 = %.4f\n",
	 //E_2body_disp_C10*HartreesToKJpermole);
	 
  printf("     Total dispersion contribution from periodic images = %.4f\n",
	 E_2body_disp_pbc*HartreesToKJpermole);
	//printf("     Total dispersion contribution from periodic images c6 = %.4f\n",
	// E_2body_disp_pbc_C6*HartreesToKJpermole);
	// printf("     Total dispersion contribution from periodic images c8 = %.4f\n",
	// E_2body_disp_pbc_C8*HartreesToKJpermole);
	 //printf("     Total dispersion contribution from periodic images c10 = %.4f\n",
	// E_2body_disp_pbc_C10*HartreesToKJpermole);
	  
	 

  E_2body_disp += 0.5*E_2body_disp_pbc;
  printf("     Final Total dispersion  = %.4f kJ/mol\n",
	 E_2body_disp*HartreesToKJpermole);
	
	//E_2body_disp_C6 += 0.5*E_2body_disp_pbc_C6;
 // printf("     dispersion from C6 contribution  = %.4f kJ/mol\n",
	// E_2body_disp_C6*HartreesToKJpermole);
	 
	//E_2body_disp_C8 += 0.5*E_2body_disp_pbc_C8;
 // printf("     dispersion from C8 contribution  = %.4f kJ/mol\n",
	// E_2body_disp_C8*HartreesToKJpermole);
	 
	//E_2body_disp_C10 += 0.5*E_2body_disp_pbc_C10;
 // printf("     dispersion from C10 contribution   = %.4f kJ/mol\n",
	// E_2body_disp_C10*HartreesToKJpermole);
	 
	//printf(" percentage analyses: C6/tot = %f, C8/tot = %f, C10/tot =%f\n",
	// E_2body_disp_C6/E_2body_disp, E_2body_disp_C8/E_2body_disp, E_2body_disp_C10/E_2body_disp); 
	
 
  stop_dispersion_time = time(NULL);
  double dispersion_time = difftime(stop_dispersion_time, start_dispersion_time);
  printf("  Time spent computing AIFF 2-body dispersion = %.0f seconds\n\n", dispersion_time);

  return E_2body_disp;
}


// Computes 3-body dispersion in either a finite or periodic system
double Cluster::ComputeManyBodyDispersion() {

  printf("Step 3: Evaluate 3-body dispersion interactions\n");

  double E_3body_disp = 0.0;
  double monomer_3b = 0.0;
  double dimer_3b = 0.0;
  double trimer_3b = 0.0;
  
  /*
   When we talk about 3-body interactions in the HMBI model, we are
   really interested in the dispersion between three monomers.  A
   broader definition of 3-body terms would include all 3-body
   interatomic (rather than intermolecular) interactions.  This
   broader definition includes:

   ** 3 atoms on a single monomer (monomer_3b) 
   ** 3 atoms spread over two monomers (dimer_3b) 
   ** 3 atoms spread over three monomers (trimer_3b) 

   In HMBI, the monomer_3b terms are treated quantum mechanically, so
   we don't need them.  The dimer_3b terms are also largely done with
   QM.  If we do local truncation and use MM for the long-range
   two-body interactions, there will be some MM dimer_3b terms.
   However, 3-body interactions die off quickly, so those terms are
   probably very small for any reasonable cutoff and can be neglected.
   So the HMBI model really only needs the trimer_3b terms.  However,
   for the sake of curiosity, this code here computes all of them, but
   we use only the trimer_3b terms.

  */

  // Monomer terms
  for (int i=1;i<=NMon;i++) {
    monomer_3b += Monomers[i].ComputeThreeBodyDispersion();
  }
  // Dimer terms
  for (int i=1;i<=NDim;i++) {
    dimer_3b += Dimers[i].ComputeThreeBodyDispersion();
  }
  // If periodic, add image dimer terms.
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      dimer_3b += DimerImages[i].ComputeThreeBodyDispersion();
    }
  }
  // Trimer terms - the ones we really care about
  trimer_3b = ComputeThreeBodyDispersion();


  double sum = monomer_3b + dimer_3b + trimer_3b;
  printf("\n--------------------------------\n");
  printf("  3-Body Dispersion Results\n");
  printf("--------------------------------\n");
  printf("  Monomers: %10.4f kJ/mol\n",monomer_3b*HartreesToKJpermole);
  printf("    Dimers:%10.4f kJ/mol\n",dimer_3b*HartreesToKJpermole);
  printf("   Trimers: %10.4f kJ/mol\n",trimer_3b*HartreesToKJpermole);
  printf("--------------------------------\n");
  printf("     Total: %10.4f kJ/mol\n",sum*HartreesToKJpermole);
  printf("--------------------------------\n");
  
  printf("Note: Only trimer contribution is added to HMBI energy\n\n");
  E_3body_disp = trimer_3b;

  return E_3body_disp;  
}

bool Cluster::IsFrozenUnitCellParam(int i) {
  return FreezeLatParams[i];
}


// Custom Basis Sections:
string Cluster::ReadHBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$H" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "")  rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadCBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$C" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadNBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$N" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadOBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$O" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadSBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$S" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadClBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Cl" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadIBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$I" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadSnBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Sn" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadPBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$P" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadKBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$K" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string Cluster::ReadNaBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Na" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadBrBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Br" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadFBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$F" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadWBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$W" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  
  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string Cluster::ReadVBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$V" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  
  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}
