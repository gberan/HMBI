// for getpid:
#include <sys/types.h>
#include <unistd.h>
// other:
#include "cluster.h"
#include "atom.h"  // by Ali
#include <sys/stat.h>
#include "opt.h"

#include "constants.h"
using namespace hmbi_constants;

#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */


// Define some tags for parallel job types
#define QCHEM_TAG 1
#define TINKER_TAG 2
#define CAMCASP_TAG 3
#define DIETAG 10
#define BUFFSIZE 1000

Cluster::Cluster() : spins(NULL), charges(NULL), Natoms_per_monomer(NULL), 
		     Monomers(NULL), Dimers(NULL), DimerImages(NULL),
		     AtomicSymbols(NULL), AtomicNumbers(NULL), AtomicMasses(NULL),
		     reciprocal_cell(NULL), unit_cell(NULL), 
		     Grad_key(NULL) {


  Energy_HMBI = 0.0;

}

void Cluster::Initialize(ifstream &infile, int Nproc) {

  // Print Header
  PrintHeader(); 

  // Search for job title string
  //Title = GetJobTitle(infile);
  //printf("Job: %s\n",Title.c_str() );

  // Print the input file back out to keep a record of what was run.
  PrintInputFile(infile);

  // Read the HMBI job parameters
  ReadHMBIParameters(infile);  // Note, QC_rem set later on.
  
  // Read the AIFF parameters  --- by Ali
  ReadAIFFParameters(infile); 

  // Create QM & MM directories, if needed
  CreateWorkingDirectories();

  // Print out information about external software packages in use
  PrintExternalPrograms();

  // Count how many monomers are present
  NMon = FindNumberOfMonomers(infile);
  NDim = NMon*(NMon-1)/2;
  printf("The system contains %d monomers, %d dimers\n\n",NMon,NDim);

  // Define the conversion factor for hartrees -> kJ/mol/monomer
  HartreeToKJpMM = HartreesToKJpermole/(double) NMon;

  // Find number of atoms in each monomer and in the total cluster
  CountNumberOfAtoms(infile);

  // Read in Q-Chem $rem section
  string qc_rem = ReadQChemRemSection(infile);
  Params::Parameters().SetParameter("QC_REM",qc_rem); 
  if (Params::Parameters().PrintLevel() > 0) 
    printf("Q-Chem $rem section:\n%s",qc_rem.c_str() );

  // Read in Q-Chem $basis section
  string qc_basis = ReadQChemBasisSection(infile);
  Params::Parameters().SetParameter("QC_BASIS",qc_basis);
  if (Params::Parameters().PrintLevel() > 0)
    printf("Q-Chem $basis section:\n%s",qc_basis.c_str() );

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

  // Read current coordinates into the vector AtomicCoordinates
  // This gets used if we do geometry optimizations, for example
  // Initialize the vector first
  AtomicCoordinates.Initialize( 3 * GetTotalNumberOfAtoms() ); 
  ReadCurrentCoordinates();

  // Generate the dimers.  Note, counting starts at 1;
  Dimers = new Dimer[NDim+1];
  int index = 1;
  for (int a=1;a<NMon;a++)
    for (int b=a+1;b<=NMon;b++) {
      //printf("index = %d, d(%d,%d)\n",index,a,b);
      Dimers[index].Initialize(Monomers[a], Monomers[b]);
      index++;
    }

  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the monomers and the dimers
    printf("*** Print out the monomers\n");
    for (int i=1;i<=NMon;i++) {
      Monomers[i].PrintAll();
    }
    printf("*** Print out the dimers\n");
    for (int i=1;i<=NDim;i++) {
      Dimers[i].PrintAll();
    }
  }

  // Create a list of the Atomic symbols
  SetAtomicSymbols();

  // GJB: exploratory
  if ( Params::Parameters().DoAtomicInductionDamping() ) {
    ReadAtomicInductionDampingFactors(infile);
  }

  // Read in the dispersion atom types if doing 3-body ATM dispersion
  if ( Params::Parameters().EstimateThreeBodyDispersion() ) {
    ReadDispersionAtomTypes(infile);
  }
  if ( Params::Parameters().DoMP2DispersionCorrection() ) {
    ReadDispersionCoefficients(infile);
  }

  // Handle periodic boundary conditions (PBC)
  if ( Params::Parameters().IsPeriodic() ) {
    // Read in the unit cell parameters
    if ( Params::Parameters().ReadLatticeVectors() )
      OldReadUnitCellDimensions(infile);
    else    
      ReadUnitCellDimensions(infile);

    // If using fractional coordinates, convert to Cartesian (Angstroms)
    if ( Params::Parameters().ReadFractionalCoordinates() ) {
      ConvertFractionalToCartesianCoordinates();
	exit(1);
    }

    // Turn on local truncation
    Params::Parameters().SetParameter("LOCAL_2_BODY","TRUE");


    
    // Identify real space periodic images
    CreatePeriodicImageDimerList();

  }
  else
    NDim_images = 0;

  fflush(stdout); // flush the output stream

  // Initialize Gradient memory, if needed
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    Grad_MM.Initialize( 3*GetTotalNumberOfAtoms() );
    
    // Create a list of the starting index for each monomer in the full
    // gradient.  Like the Monomers, the key indexes count from 1->NMon.
    // On the other hand, in the gradient, the indexing starts at 0.
    Grad_key = new int[NMon+1];
    Grad_key[1] = 0; // set the first one by hand
    for (int i=2;i<=NMon;i++) {
      Grad_key[i] = Grad_key[i-1] + 3*Monomers[i-1].GetNumberOfAtoms();
    }
  }



  if ( Params::Parameters().IsPeriodic() )
    WriteCrystalXYZFile(2,2,2);

  // Done with all the initialization.  
  // Run RunJobsAndComputeEnergy() to start doing the hard work.

}

// Destructor
Cluster::~Cluster() {

  delete [] Monomers;
  delete [] Dimers;

  delete [] spins;
  delete [] charges;
  delete [] Natoms_per_monomer;

  delete [] AtomicSymbols;
  delete [] AtomicNumbers;
  delete [] AtomicMasses;
  
  if ( Params::Parameters().IsPeriodic() ) {
    delete [] MonomerImages;
    delete [] DimerImages;
    delete [] unit_cell;
    delete [] reciprocal_cell;
    if ( Params::Parameters().DoAIFFInduction() ) {
      delete [] InducedMultipoles;
    }
  }

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    delete [] Grad_key;
  }

  if ( Params::Parameters().EstimateThreeBodyDispersion() ) {
    delete [] DispersionAtomTypes;
  }


}

// Creates & runs the jobs, computes the HMBI energy from the results
void Cluster::RunJobsAndComputeEnergy() {

  if ( Params::Parameters().DoCounterpoise() ) 
    printf("Including counterpoise correction for basis set superposition error\n     in dimer jobs.\n");


  printf("\n");

  // Create and Run QM and MM calculations, if requested 
  if ( Params::Parameters().RunJobs() ) {

    // Create the jobs
    if ( !Params::Parameters().BuildForceFieldOnly())
	CreateQMJobs();
    if ( !Params::Parameters().NeglectManyBody() )
      CreateMMJobs();

	printf("Number of dimer jobs = %d\n",NDim+NDim_images); //HACK
	//exit(1);

	if (Params::Parameters().CreateJobsOnly()) {
	  // create the shell scripts for camcasp, if needed, then exit.
	  if (Params::Parameters().GetMMType()==2) {
	    for (int i=1;i<=NMon;i++) {
	      Monomers[i].RunCamCaspJob(); // create the shell script
	    }
	  }
	  printf("QM and MM jobs created.  Exiting.\n");
	    exit(0);
	}

    // Now run the jobs
    RunJobs();
    
    // If we only wanted to run the QM & MM jobs, exit now.
    if ( Params::Parameters().RunJobsOnly()) {
      printf("QM and MM jobs complete.  Exiting.\n");
      exit(1);
    }


    // If we only wanted the force field, exit now
    if ( Params::Parameters().BuildForceFieldOnly()) {
      printf("Force field parameters have been calculated and saved.  Exiting.\n");
      exit(0);
    }
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
 
  if (! Params::Parameters().BuildForceFieldOnly()) {
    ReadQMResults();
  }
  if (! Params::Parameters().NeglectManyBody() ) {
    ReadMMResults();
  }
  
  /* Process QM and MM results */

  // Perform HMBI calculations
  Energy_HMBI = ComputeHMBIEnergy();
  if ( Params::Parameters().GetMMType() == 4)
    printf("EE-PA Energy = %15.9f hartrees\n",Energy_HMBI);
  else
    printf("HMBI Energy = %15.9f hartrees\n",Energy_HMBI);
  
  printf("\n");
  

  // A simple routine for playing with different possible local
  // truncation parameters,
  printf("Params: ScanLocalTruncationParams = %d\n",Params::Parameters().ScanLocalTruncationParameters() );
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
      
      printf("\n");
    }
  } // End ScanLocalTruncationParameters code



  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().OrientDebug() ) {
    printf("--------------------------------------------------\n");
    printf("Debug: Orient output:\n");
    double Edebug = ReadOrientEnergy();
    printf("--------------------------------------------------\n\n");
  }


  // Optionally, Compute the HMBI gradient
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    Grad_HMBI = ComputeHMBIGradient();
    Grad_HMBI.PrintGradient("\nFull HMBI Gradient");
    printf("\n");

    if (Params::Parameters().PrintLevel() > 0 )
      PrintGradient("\nFull HMBI Gradient",Grad_HMBI);
  }
}


void Cluster::PrintHeader() {

    printf("\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" *                                           *\n");
    printf(" *  Hybrid Many-Body Interaction QM/MM Code  *\n");
    printf(" *                                           *\n");
    printf(" *       by Greg Beran, Kaushik Nanda        *\n");
    printf(" *             and Ali Sebetci               *\n");
    printf(" *                                           *\n");
    printf(" *   University of California, Riverside     *\n");
    printf(" *                                           *\n");
    printf(" *                Nov. 2009                  *\n");
    printf(" *                                           *\n"); 
    printf(" * * * * * * * * * * * * * * * * * * * * * * *\n\n");
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
    printf("\n");
 
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



  // shuhao: creat the reciprocal cell vectors;

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

  // shuhao: creat the reciprocal cell vectors; 
  
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
    //printf("GetUnitCellParameter(): a = %f, b = %f, c = %f\n",a,b,c);
    //printf("GetUnitCellParameter(): alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);

    if (type== "a" || type == "A") {return a;}
    else if (type == "b" || type == "B") {return b;}
    else if (type == "c" || type == "C") {return c;}
    else if (type == "alpha" || type == "ALPHA") {return alpha;}
    else if (type == "beta" || type == "BETA") {return beta;}
    else if (type == "gamma" || type == "GAMMA") {return gamma;}
    else {
      printf("ERROR: Cluster::SetUnitCellParameter(): Unknown parameter %s\n",type.c_str());
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
      printf("ERROR: Cluster::SetUnitCellParameter(): Unknown parameter %s\n",type.c_str());
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


}

void Cluster::SetUnitCellVectors(Vector v1, Vector v2, Vector v3) {
  unit_cell[0] = v1;
  unit_cell[1] = v2;
  unit_cell[2] = v3;
}

// Converts unit cell fractional coordinates to Cartesian coordinates
// in Angstroms.  Assumes the lattice parameters have already been
// read and processed.
void Cluster::ConvertFractionalToCartesianCoordinates() {

  /*
  Vector v1 = GetUnitCellVector(0);
  Vector v2 = GetUnitCellVector(1);
  Vector v3 = GetUnitCellVector(2);

  // create array of the proper size
  int Natoms = GetTotalNumberOfAtoms();
  Vector CartesianCoords(3*Natoms);
  CartesianCoords.Set();
  
  for (int i=0;i<Natoms;i++) {
    CartesianCoords[3*i] = AtomicCoordinates[3*i]*v1[0] 

  }
  */

  printf("Converting from Fractional to Cartesian coordinates");

  double a = UnitCellAxes[0];
  double b = UnitCellAxes[1];
  double c = UnitCellAxes[2];

  double alpha = UnitCellAngles[0]*DegreesToRadians;
  double beta = UnitCellAngles[1]*DegreesToRadians;
  double gamma = UnitCellAngles[2]*DegreesToRadians;

  double beta_term = 
    (cos(alpha) - cos(beta)*cos(gamma) ) / sin(gamma);
  double gamma_term = 
    sqrt(1-cos(beta)*cos(beta) - beta_term*beta_term);

  // Create transformation matrix
  // http://en.wikipedia.org/wiki/Fractional_coordinates
  Matrix FracToCart(3,3);
  FracToCart.Set();
  FracToCart(0,0) = a;
  FracToCart(0,1) = b*cos(gamma);
  FracToCart(0,2) = c*cos(beta);
  FracToCart(1,0) = 0.0;
  FracToCart(1,1) = b*sin(gamma);
  FracToCart(1,2) = c*beta_term;
  FracToCart(2,0) = 0.0;
  FracToCart(2,1) = 0.0;
  FracToCart(2,2) = c*gamma_term;

  FracToCart.Print("Frac->Cart Transformation Matrix");

  // create array to store transformed coordinates
  int Natoms = GetTotalNumberOfAtoms();
  Vector CartesianCoords(3*Natoms);
  CartesianCoords.Set();

  for (int iatom=0; iatom<Natoms;iatom++) {
    Vector frac(3), cart(3);
    for (int j=0;j<3;j++) {
      frac[j] = AtomicCoordinates[3*iatom+j];
    }
    cart = FracToCart.MatrixTimesVector(frac);
        for (int j=0;j<3;j++) {
      frac[j] = CartesianCoords[3*iatom+j] = cart[j];
    }
  }

  printf("Fractional:\n");
  PrintXYZCartesian();
  SetNewCoordinates(CartesianCoords);
  printf("Cartesian:\n");
  PrintXYZCartesian();

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

void Cluster::ReadGeometry(ifstream& infile) {


  // Right now we have 2 options: Tinker or Q-Chem-style XYZ geometries
  if (Params::Parameters().GetMMType()==1) {
    ReadTinkerGeometry(infile);
  }
  else if (Params::Parameters().GetMMType()==2)  {
    ReadSimpleXYZGeometry(infile);
  }
  else
    ReadSimpleXYZGeometry(infile);

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
      // GJB: if need IP
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
      Monomers[i].Initialize( i, charge, spin, atoms, mon_xyz, natoms, atypes,
			      nconn, conn, Monomer_charges);
    }

    else { // no embedding charges case
      Monomers[i].Initialize( i, charge, spin, atoms, mon_xyz, natoms, atypes,
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
	infile >> IPs[index];
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
    Monomers[i].Initialize(i,charge,spin,atoms,mon_xyz,natoms,
			   Monomer_charges);
    }
    else { // no embedding charges case
    Monomers[i].Initialize(i,charge,spin,atoms,mon_xyz,natoms);
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
      if ( Params::Parameters().GetJobType() > 1 && Params::Parameters().GetJobType() < 100 )
	    rem += "jobtype = force\n";
	
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {

	  rem += line;
	  rem += "\n";
	}
	else {
	  // Insert keywords if doing gradients/hessians
	  if ( Params::Parameters().GetJobType() > 1 ) {
	    rem += "no_reorient = true\nsym_ignore = true\n";
	    // And for Raman intensities, if needed
	    //if ( Params::Parameters().DoRamanIntensities() ) {
	    //  rem += "doraman = true\n";
	    //}
	  }
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

// Print out the paths to any external programs we are using
void Cluster::PrintExternalPrograms() {

  printf("The following external software packages will be used:\n");

  // First QM
  printf("   Q-Chem: %s\n",getenv("QC"));

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

// Updates all Monomers and Dimers with new coordinates Also erases
// any energies and gradients unless ResetEnergies flag is false.
// Careful: this flag should be true unless you have a good reason for
// it not to be.
void Cluster::SetNewCoordinates(Vector NewCoords, bool ResetEnergies) {

  /*  Deal with Monomers first: */
  int i = 0;
  double xyz[3];
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

    if (ResetEnergies) 
      Monomers[imon].ResetEnergiesAndGradients();
    
  }   

  /* Update Dimers */

  // Dimers
  for (int i=1;i<=NDim;i++) {
    int m1 = Dimers[i].GetIndexA();
    int m2 = Dimers[i].GetIndexB();
    // Update dimers with modified monomer data
    //printf("m1 = %d, m2 = %d\n",m1,m2);
    Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);
    if (ResetEnergies) 
      Dimers[i].ResetEnergiesAndGradients();

  }

  if ( Params::Parameters().IsPeriodic() ) {
    delete [] DimerImages;
    // Identify real space periodic images - local fragments might have changed during optimization
    CreatePeriodicImageDimerList();
  }

  /* Full Cluster */
  // geometry is read directly from monomers.  So just reset the energies.
  if (ResetEnergies) 
    ResetEnergiesAndGradients();

  // Finally, reset the current coordinates in the cluster object
  ReadCurrentCoordinates();
}

void Cluster::CreateQMJobs() {

  // Monomers
  printf("Creating QM Monomer jobs...  ");
  fflush(stdout); // flush the output stream

  for (int i=1;i<=NMon;i++) {
    if (Params::Parameters().TinkerDebug()) {
      Monomers[i].CreateMMJob(Monomers, NMon);
    }
    else
      Monomers[i].CreateQChemJob(Monomers, NMon);
  }
  
  // Dimers
  printf("Creating QM Dimer jobs.\n");
  fflush(stdout); // flush the output stream
  for (int i=1;i<=NDim;i++) {
    if (Params::Parameters().TinkerDebug()) 
      Dimers[i].CreateMMJob(Monomers, NMon);
    else 
      Dimers[i].CreateQChemJob(Monomers, NMon);
  }
  
  // If doing periodic boundary conditions, create additional dimer jobs
  if ( Params::Parameters().IsPeriodic() ) {
    
    for (int i=1;i<=NDim_images;i++) {
      if (Params::Parameters().UseEmbeddingCharges() ) {
	printf("Cluster::CreateQMJobs() ERROR: Embedding charges not yet implemented for periodic systems\n");
	exit(1);
      }
      if (Params::Parameters().TinkerDebug()) 
	DimerImages[i].CreateMMJob(Monomers, NMon);
      else
	DimerImages[i].CreateQChemJob(Monomers,NMon);	
    }
  }
}

void Cluster::CreateMMJobs() {
    // Monomers
    printf("Creating MM Monomer jobs...  ");
    fflush(stdout); // flush the output stream
    for (int i=1;i<=NMon;i++) {
      Monomers[i].CreateMMJob(Monomers, NMon);
    }
    
    // Dimers
    printf("Creating MM Dimer jobs.\n\n");
    fflush(stdout); // flush the output stream
    for (int i=1;i<=NDim;i++) {
      Dimers[i].CreateMMJob(Monomers, NMon);
    }
    
    // If doing periodic boundary conditions, create additional dimer jobs
    if ( Params::Parameters().IsPeriodic() ) {
      
      for (int i=1;i<=NDim_images;i++) {
	if (Params::Parameters().UseEmbeddingCharges() ) {
	  printf("Cluster::CreateQMJobs() ERROR: Embedding charges not yet implemented for periodic systems\n");
	  exit(1);
	}
	
	DimerImages[i].CreateMMJob(Monomers,NMon);
      }
    }
    
    // Full Cluster
    if ( Params::Parameters().DoQMBenchmark() ) { // full cluster benchmark 
      CreateFullClusterQChemJob(false);
    }

    if ( Params::Parameters().GetMMType()==1 ) { // Tinker
      CreateFullClusterTinkerJob();
    }
    else if ( Params::Parameters().GetMMType()==2 ) { // AIFF
      // The AIFF energy for the full cluster is evaluated later, when
      // we read the MM energies.
    }
    else if ( Params::Parameters().GetMMType()==3 ) { // Q-Chem
      CreateFullClusterQChemJob(true);
    }
    else {
      printf("Cluster::CreateMMJobs: Unknown MM_type: %d\n",
	     Params::Parameters().GetMMType() );
      exit(1);
    }

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

void Cluster::CreateFullClusterTinkerJob() {
  // Full cluster
  string xyzfile = Params::Parameters().GetMMPath() + "/full.xyz"; 
  string keyfile = Params::Parameters().GetMMPath() + "/full.key"; 

  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Monomer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }
  PrintTinkerCartesian(xyz);
  fclose(xyz);


  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
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


void Cluster::CreateFullClusterQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.force'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  
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

  /* Step 1: create list of jobs to run */
  // Count the jobs
  int N_QM_jobs = NMon + NDim + NDim_images;
  int N_MM_jobs;
  if ( Params::Parameters().GetMMType() == 2 )
    N_MM_jobs = NMon; // for AIFF, only do monomer QM jobs here
  else
    N_MM_jobs = N_QM_jobs + 1; // otherwise,  do all of them.

  if ( Params::Parameters().DoQMBenchmark() )
    N_QM_jobs++;


  int Njobs = N_QM_jobs + N_MM_jobs;
  //printf("QM jobs = %d, MM = %d, total = %d\n",N_QM_jobs,N_MM_jobs,Njobs);

  // Create and Combine the job lists
  string *JobList = new string[Njobs];
  string *QMList = NULL;// = new string[N_QM_jobs];
  string *MMList = NULL;// = new string[N_MM_jobs];

  if (Params::Parameters().BuildForceFieldOnly()) {
    printf("Skipping QM jobs\n");
    Njobs -= N_QM_jobs;
    N_QM_jobs = 0;
  }
  else {
    QMList = RunQMJobs();
    for (int i=0;i<N_QM_jobs;i++) {
      JobList[i] = QMList[i];    
    }
  }
    
  if (Params::Parameters().NeglectManyBody() ) {
    printf("Many-body terms being neglected - Skipping MM jobs.\n");
    Njobs -= N_MM_jobs;
  }
  else {
    MMList = RunMMJobs();
  
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
  if (Nproc > 1 ) {

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
    if ( Params::Parameters().GetMMType() == 2) 
      ExpensiveJobs += NMon; // AIFF monomer jobs
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
      if (ijob < N_QM_jobs) 
	tag = QCHEM_TAG;
      else {
	if (Params::Parameters().GetMMType() == 1 )
	  tag = TINKER_TAG;
	else if (Params::Parameters().GetMMType() == 2 )
	  tag = CAMCASP_TAG;
	else if (Params::Parameters().GetMMType() == 3 )
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
      if (ijob < N_QM_jobs) 
	tag = QCHEM_TAG;
      else {
	if (Params::Parameters().GetMMType() == 1 )
	  tag = TINKER_TAG;
	else if (Params::Parameters().GetMMType() == 2 )
	  tag = CAMCASP_TAG;
	else if (Params::Parameters().GetMMType() == 3 )
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

    printf("Time spent on QM jobs: %0.f sec.    MM jobs: %0.f sec\n",QM_time,MM_time);


#ifdef PARALLEL
  }
#endif /* PARALLEL */
  
  if ( Params::Parameters().GetMMType() == 2  && Params::Parameters().OrientDebug() ) { // Ali
    system(RunOrientJob().c_str());
  }
  
  delete [] JobList;

}

// Returns a list of commands for running each job
string* Cluster::RunQMJobs() {

  int Njobs = NMon + NDim + NDim_images;
  if ( Params::Parameters().DoQMBenchmark() ) 
    Njobs++;

  string *JobList = new string[Njobs];
  
  int ijob = 0;

  // Monomers
  for (int i=1;i<=NMon;i++) {
    if (Params::Parameters().TinkerDebug()) 
      JobList[ijob] = Monomers[i].RunMMJob();
    else    
      JobList[ijob] = Monomers[i].RunQChemJob();
    ijob++;
  }

  // Dimers
  for (int i=1;i<=NDim;i++) {
    if (Params::Parameters().TinkerDebug()) 
      JobList[ijob] = Dimers[i].RunMMJob();
    else
      JobList[ijob] = Dimers[i].RunQChemJob();
    ijob++;
  }
  
  // If doing periodic boundary conditions, have additional dimer jobs
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1;i<=NDim_images;i++) {
      if (Params::Parameters().TinkerDebug()) 
	JobList[ijob] = DimerImages[i].RunMMJob();
      else    
	JobList[ijob] = DimerImages[i].RunQChemJob();
      ijob++;
    }
  }
  
  if (Params::Parameters().IsPeriodic() ) 
    printf("Need to run %d QM monomer jobs, %d dimer jobs, %d dimer image jobs.\n",NMon,NDim,NDim_images);
  else
    printf("Need to run %d QM monomer jobs and %d dimer jobs.\n",NMon,NDim);

  // If benchmarking with full QM job
  if ( Params::Parameters().DoQMBenchmark() ) {
    printf("Also running full cluster QM job as a benchmark\n");
    JobList[ijob] = RunFullClusterQChemJob(false);
    ijob++;
  }

  return JobList;
}

// Returns a list of commands for running each job
string* Cluster::RunMMJobs() {

  int Njobs = NMon + NDim + NDim_images + 1;
  if ( Params::Parameters().GetMMType() == 2 ) {// AIFF - monomers only for now
    Njobs = NMon;
  } 

  string *JobList = new string[Njobs];
  
  int ijob = 0;

  // Monomers
  for (int i=1;i<=NMon;i++) {
    JobList[ijob] = Monomers[i].RunMMJob();
    ijob++;
  }
  
  // If we are not using the AIFF, we also evaluate dimers and full
  // cluster here.  If we are using the AIFF, this happens later, when
  // we read the energies.  For the AIFF, only the monomer properties
  // are determined here.
  if ( Params::Parameters().GetMMType() != 2 ) {// not the AIFF
    // Dimers
    for (int i=1;i<=NDim;i++) {
      JobList[ijob] = Dimers[i].RunMMJob();
      ijob++;
    }
    
    // If doing periodic boundary conditions, have additional dimer jobs
    if ( Params::Parameters().IsPeriodic() ) {
      printf("Including MM Dimer Image jobs\n");
      for (int i=1;i<=NDim_images;i++) {
	JobList[ijob] = DimerImages[i].RunMMJob();
	ijob++;
      }
    }
  }

  // Full cluster
  if ( Params::Parameters().GetMMType() == 1 ) {// Tinker
    JobList[ijob] = RunFullClusterTinkerJob();
    ijob++;
  }
  else if ( Params::Parameters().GetMMType() == 3 ) {// QChem
    JobList[ijob] = RunFullClusterQChemJob(true);
    ijob++;
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

// Returns command to run the full cluster tinker job
string Cluster::RunFullClusterTinkerJob() {

  // The job:
  string infile = "full.xyz";
  string outfile = "full.out";

  // To Execute Tinker:
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
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
  
  return cmd;

}


void Cluster::ReadQMResults() {
  //printf("Cluster::ReadQMResults()\n");

  // Monomers
  for (int i=1;i<=NMon;i++) {
    if ( Params::Parameters().TinkerDebug() ) {
      Monomers[i].ReadMMResults();
      Monomers[i].SetQMEnergy( Monomers[i].GetMMEnergy() );
      Monomers[i].SetQMGradient( Monomers[i].GetMMGradient() );
    }
    else
      Monomers[i].ReadQMResults();
  }

  // Dimers
  for (int i=1;i<=NDim;i++) {
    int m1 = Dimers[i].GetIndexA();
    int m2 = Dimers[i].GetIndexB();
    // Update dimers with modified monomer data
    //printf("m1 = %d, m2 = %d\n",m1,m2);
    Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);
    if ( Params::Parameters().TinkerDebug() ) {
      Dimers[i].ReadMMResults();
    }
    else 
      Dimers[i].ReadQMResults();

    if ( Params::Parameters().TinkerDebug() ) {
      Dimers[i].SetQMGradient( Dimers[i].GetMMGradient() );
    }
  }

  // Additional dimers for periodic system
  if ( Params::Parameters().IsPeriodic() ) {
    //printf("Reading Dimer_image QM results\n");
    for (int i=1;i<=NDim_images;i++) {
      int m1 = DimerImages[i].GetIndexA();
      int m2 = DimerImages[i].GetIndexB();
      int ref_mon = DimerImages[i].GetReferenceMonomerIndex();
      
      // Since we don't store a separate list of all the image monomers, we
      // grab it directly from the dimer object.
      Monomer MonB = DimerImages[i].GetMonomerB();
      
      //printf("ReadQMResults: Dimer (%d,%d), where B refers to monomer %d\n",m1,m2,ref_mon);
      //fflush(stdout);
      DimerImages[i].UpdateObjects(Monomers[m1],MonB,Monomers[ref_mon]);
      //DimerImages[i].UpdateObjects(Monomers[m1],MonB);
      //if ( Params::Parameters().TinkerDebug() ) {
      //Dimers[i].ReadMMResults();
      //DimerImages[i].SetQMGradient( DimerImages[i].GetMMGradient() );
      //}
      //else
      //printf("ReadQMResults for the DimerImages\n"); fflush(stdout);
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

  // If we did a benchmark full cluster QM calculation
  if ( Params::Parameters().DoQMBenchmark() ) {
    Energy_QM = ReadQChemEnergy(false);
  }


}

void Cluster::ReadMMResults() {

  // Monomers
  for (int i=1;i<=NMon;i++) {
    Monomers[i].ReadMMResults();
  }

  // Dimers
  for (int i=1;i<=NDim;i++) {
    int m1 = Dimers[i].GetIndexA();
    int m2 = Dimers[i].GetIndexB();
    // Update dimers with modified monomer data
    Dimers[i].UpdateObjects(Monomers[m1],Monomers[m2]);
    Dimers[i].ReadMMResults();
  }

  // Additional dimers for periodic system
  if ( Params::Parameters().IsPeriodic() ) {

    /*
    if ( Params::Parameters().GetMMType()==2) {
	printf("Warning: Trying to do PBC AIFF electrostatics.  Not implemented\n");
	//exit(1);
    }
    */

    printf("Reading Dimer_image MM results\n");
    for (int i=1;i<=NDim_images;i++) {
      int m1 = DimerImages[i].GetIndexA();
      int m2 = DimerImages[i].GetIndexB();
      int ref_mon = DimerImages[i].GetReferenceMonomerIndex();
      
      Monomer MonB = DimerImages[i].GetMonomerB();
      DimerImages[i].UpdateObjects(Monomers[m1],MonB,Monomers[ref_mon]);
      DimerImages[i].ReadMMResults( Monomers[ref_mon] );
    }
    printf("Done updating MM in PBC objects\n"); fflush(stdout);
  }


  // Full Cluster
  // Read in the energy of the full cluster at the MM level
  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) { 
    SetMMGradient(); // get energies & gradient simultaneously in this case
    //PrintGradient("Full MM Gradient",Grad_MM);
  }
  else {
    SetMMEnergy();
    if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
      SetMMGradient();
     //PrintGradient("Full MM Gradient",Grad_MM);
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
 
  MM_Grad_Init = 1;
}

// Main routine for reading in gradient 
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
// type 1 = qchem style, type 2 = tinker style
Vector Cluster::ReadGradient(int type) {
  
  string path = Params::Parameters().GetMMPath();
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
/*
Vector Cluster::ReadQMCPGradient(int type) {

  string path = Params::Parameters().GetMMPath();
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
    if (type=="QM")
      dE = Dimers[i].GetQMIntEnergy();
    else {
      dE = Dimers[i].GetMMIntEnergy();
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

    // at present, we don't have any symmetry for non-PBC cases
    dE *= Dimers[i].GetSymmetryFactor();
    if (Params::Parameters().GetMMType() == 2) {
      dEes *= Dimers[i].GetSymmetryFactor();
      dEind *= Dimers[i].GetSymmetryFactor();
      dEdisp *= Dimers[i].GetSymmetryFactor();
    
      //dEdispc6 *= Dimers[i].GetSymmetryFactor();
      //dEdispc8 *= Dimers[i].GetSymmetryFactor();
      //dEdispc10 *= Dimers[i].GetSymmetryFactor();
    }
    if ( Params::Parameters().DoLocal2BodyTruncation() ) {
      double damp = Dimers[i].GetDampingFactor(c0,c1);
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
	if (type=="QM")
	  dE = DimerImages[i].GetQMIntEnergy();
	else {
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
	  /*
	  int iA = DimerImages[i].GetIndexA();
	  int iB = DimerImages[i].GetIndexB();
	  int *K_vector = NULL;
	  K_vector = DimerImages[i].GetImageCell();
	  
	  printf("E(%3d %3d ) = %9.3f  E(%3d) = %9.3f E(%3d) = %9.3f  Eint = %9.3f  Ref = %d, k=(%d,%d,%d)  sym_fac = %d\n",
		 iA,iB,DimerImages[i].GetMMEnergy()*2625.5, 
		 iA, DimerImages[i].GetMonomerA().GetMMEnergy()*2625.5,
		 iB, DimerImages[i].GetMonomerB().GetMMEnergy()*2625.5,
		 DimerImages[i].GetMMIntEnergy()*2625.5,
		 DimerImages[i].GetReferenceMonomerIndex(),K_vector[0], K_vector[1], K_vector[2], DimerImages[i].GetSymmetryFactor());
	  */
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
  if ( Params::Parameters().NeglectManyBody() || Params::Parameters().TinkerDebug() ) {
    printf("Neglecting many-body terms\n");
    scafac = 0.0; 
    Efull_MM = 0.0;
  }

  else {
    Efull_MM = GetMMEnergy();
    scafac = 1.0;
  }

  //HACK
  //Efull_MM = 0.0;
  //printf("HACK!!! Efull_MM = 0!\n");
  // End Hack

  double QM_1b, QM_2b;

  if (! Params::Parameters().BuildForceFieldOnly()) {
    QM_1b = GetTotalOneBodyEnergy("QM");
    QM_2b = GetTotalTwoBodyEnergy("QM");
  }
  else {
    QM_1b = 0.0; QM_2b = 0.0;
  }

  double MM_1b = GetTotalOneBodyEnergy("MM");
  double MM_2b = GetTotalTwoBodyEnergy("MM");


  //printf("HACK: zeroing out QM & MM 1 & 2-body terms\n");
  //QM_1b = 0.0;
  //QM_2b = 0.0;
  //MM_1b = 0.0;
  //MM_2b = 0.0;
  //Efull_MM = 0.0;

  double MM_mb = Efull_MM - scafac*MM_1b - scafac*MM_2b;
  
  printf("QM 1-body    = %15.9f hartrees\n",QM_1b);
  printf("QM 2-body    = %15.9f hartrees\n",QM_2b);

  printf("MM full      = %9.3f kJ/mol\n",Efull_MM*2625.5);
  printf("MM 1-body    = %9.3f kJ/mol\n",MM_1b*2625.5);
  printf("MM 2-body    = %9.3f kJ/mol\n",MM_2b*2625.5);
  printf("MM many-body = %9.3f kJ/mol\n",MM_mb*2625.5);

  // Combine QM 1- & 2-body with MM many-body to get HMBI energy
  energy = QM_1b + QM_2b + MM_mb;

  //double E12_QM = GetTotalOneBodyEnergy("QM") + GetTotalTwoBodyEnergy("QM");
  double E12_QM = QM_1b + QM_2b;
  printf("\nQM 1+2-body energy = %15.9f hartrees\n",E12_QM);

  //double E_mb_MM = Efull_MM - scafac*GetTotalOneBodyEnergy("MM")
  //  - scafac*GetTotalTwoBodyEnergy("MM");
  if ( Params::Parameters().GetMMType() != 4)
    printf("MM many-body energy = %.3f kJ/mol (%.3f kJ/mol per monomer)\n",
	   MM_mb*HartreesToKJpermole, MM_mb*HartreeToKJpMM);
  
  return energy;
}

void Cluster::PrintGradient(string title, const Vector& grad) {
  
  int Natoms = GetTotalNumberOfAtoms();

  printf("%s\n",title.c_str());
  for (int i=0;i<Natoms;i++) {
    printf("%2s %15.10f %15.10f %15.10f\n",GetAtomicSymbol(i).c_str(),
    	   grad[3*i],grad[3*i+1],grad[3*i+2]);
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
  */

  int Ntot = GetTotalNumberOfAtoms();
  Vector grad(3*Ntot);
  double scafac; 
  
  if (Params::Parameters().NeglectManyBody() || Params::Parameters().TinkerDebug() ) {
    scafac = 0.0; // neglect many-body terms
  }
  else {
    scafac = 1.0; // treat them as usual
  }
  
  // Contribution from full cluster MM gradient
  if (! Params::Parameters().NeglectManyBody() )
    grad = Grad_MM;
  else {
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
  }

  grad.Scale(scafac);

  // Create a list of the starting index for each monomer in the 
  // full gradient.  Like the Monomers, The key indexes count from 1->NMon.
  // On the other hand, in the gradient, the indexing starts at 0.
  int *key = new int[NMon+1];
  key[1] = 0; // set the first one by hand
  for (int i=2;i<=NMon;i++) {
    key[i] = key[i-1] + 3*Monomers[i-1].GetNumberOfAtoms();
  }


  // 1-body contributions
  for (int i=1;i<=NMon;i++) {
    int Na = Monomers[i].GetNumberOfAtoms();
    int start = key[i];
    for (int j=0;j<3*Na;j++) {
      grad[start+j] += Monomers[i].GetQMGradient()[j];
      //if (Params::Parameters().GetMMType() != 2) // if not AIFF
      grad[start+j] -= scafac*Monomers[i].GetMMGradient()[j] ;
    }
  }
  
  // grab the cutoffs for the QM -> MM damping
  double c0 = Params::Parameters().GetLocalCutoff(0);
  double c1 = Params::Parameters().GetLocalCutoff(1);

  // if ( !(Params::Parameters().DoCounterpoise() ) ) {
  // 2-body interaction contributions
  for (int i=1;i<=NDim;i++) {
    int Na = Dimers[i].GetMonomerA().GetNumberOfAtoms();
    int Nb = Dimers[i].GetMonomerB().GetNumberOfAtoms();
    int indexA = Dimers[i].GetIndexA();
    int indexB = Dimers[i].GetIndexB();
    int startA = key[indexA];
    int startB = key[indexB];

    double damp = Dimers[i].GetDampingFactor(c0,c1);
    //printf("Dimer %d, damping factor = %f\n",i,damp);

    // piece arising from first monomer
    /*    if (Params::Parameters().GetMMType() == 2) { // AIFF
      for (int j=0;j<3*Na;j++) {
	grad[startA+j] += 
	  Dimers[i].GetQMGradient()[j] - scafac*Dimers[i].GetMMGradient()[j];
      }
      for (int j=0;j<3*Nb;j++) {
	int k = j + 3*Na; // in dimer gradient, need offset to second monomer
	grad[startB+j] += 
	  Dimers[i].GetQMGradient()[k]- scafac*Dimers[i].GetMMGradient()[k];
      }

    }
    else {*/
      for (int j=0;j<3*Na;j++) {
	grad[startA+j] += 
	  damp*(( Dimers[i].GetQMGradient()[j] - Monomers[indexA].GetQMGradient()[j] )
		- scafac*(Dimers[i].GetMMGradient()[j] - Monomers[indexA].GetMMGradient()[j]) );
      }
      // piece arising from second monomer
      for (int j=0;j<3*Nb;j++) {
	int k = j + 3*Na; // in dimer gradient, need offset to second monomer
	grad[startB+j] += 
	  damp*(( Dimers[i].GetQMGradient()[k] - Monomers[indexB].GetQMGradient()[j])
		 - scafac*(Dimers[i].GetMMGradient()[k] - Monomers[indexB].GetMMGradient()[j]) );
      }
      //}
  }
  // }

/*
 if ( Params::Parameters().DoCounterpoise() ) {
  // 2-body interaction contributions
  for (int i=1;i<=NDim;i++) {            
    int Na = Dimers[i].GetMonomerA().GetNumberOfAtoms();
    int Nb = Dimers[i].GetMonomerB().GetNumberOfAtoms();
    int indexA = Dimers[i].GetIndexA();
    int indexB = Dimers[i].GetIndexB();           
    int startA = key[indexA];
    int startB = key[indexB];

      for (int j=0;j<3*Na;j++) {
        grad[startA+j] +=
          ( Dimers[i].GetQMCPGradient()[j] - Monomers[indexA].GetQMGradient()[j])
          - scafac*(Dimers[i].GetMMGradient()[j] - Monomers[indexA].GetMMGradient()[j]);
      }
      // piece arising from second monomer              
      for (int j=0;j<3*Nb;j++) {       
        int k = j + 3*Na; // in dimer gradient, need offset to second monomer
        grad[startB+j] +=    
          ( Dimers[i].GetQMCPGradient()[k] - Monomers[indexB].GetQMGradient()[j])
          - scafac*(Dimers[i].GetMMGradient()[k] - Monomers[indexB].GetMMGradient()[j]);
      }

  }
 }
*/
  //PrintGradient("HMBI Gradient before PBC contrib",grad);

  // 2-body contributions due to periodic boundary conditions
  if (Params::Parameters().IsPeriodic() ) {
      printf("Looking at PBC contributions to the gradient\n");
    for (int i=1;i<=NDim_images;i++) {
      int indexA = DimerImages[i].GetIndexA();
      int indexB = DimerImages[i].GetIndexB();
      int indexB_ref = DimerImages[i].GetReferenceMonomerIndex();
      int Na = DimerImages[i].GetMonomerA().GetNumberOfAtoms();
      int Nb = DimerImages[i].GetMonomerB().GetNumberOfAtoms();
      int startA = key[indexA];
      int startB = key[indexB_ref];

      double damp = DimerImages[i].GetDampingFactor(c0,c1);
      //printf("DimerImage %d, damping factor = %f\n",i,damp);

      double symfac = 1.0;
      if ( Params::Parameters().UseCrystalSymmetry() )
	symfac *= DimerImages[i].GetSymmetryFactor();
      //printf("Symmetry factor = %f\n",symfac);

      // piece arising from first monomer
      for (int j=0;j<3*Na;j++) {
	Monomer mA = Monomers[indexA];
	
	grad[startA+j] += 
	  0.5*symfac*damp*( ( DimerImages[i].GetQMGradient()[j] - mA.GetQMGradient()[j])
	   - scafac*( DimerImages[i].GetMMGradient()[j] - mA.GetMMGradient()[j] ) );	
      }

      // Piece arising from second monomer 
      for (int j=0;j<3*Nb;j++) {
	int k = j + 3*Na; // in dimer gradient, need offset to second monomer
	Monomer mB = Monomers[indexB_ref];
	
	grad[startB+j] += 
	  0.5*symfac*damp*( ( DimerImages[i].GetQMGradient()[k] - mB.GetQMGradient()[j] )
	   - scafac*( DimerImages[i].GetMMGradient()[k] - mB.GetMMGradient()[j] ) );
	
      }
    }

  
    printf("Calling debug lattice gradient routine\n");
    DebugLatticeGradient();

  }
  
  //PrintGradient("HMBI Gradient",grad);

  delete [] key;

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

void Cluster::PrintTinkerCartesian(FILE *outfile) {

  // Print tile line
  fprintf(outfile,"%d  Full Cluster\n",GetTotalNumberOfAtoms()); 
  int shift = 0;

  // apply shifts to ensure proper indexing of atom numbers/connectivity
  for (int i=1;i<=NMon;i++) {
    Monomers[i].PrintTinkerCartesian(shift,false,outfile);
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
  else if (Params::Parameters().GetMMType() == 2) {// Qchem2
    string tmp = Params::Parameters().GetQChemRem2();
    tmp.erase(0,5); //erase "$rem" & replace it with "$qchem2"
    fprintf(outfile,"$qchem2%s", tmp.c_str() );
  }
  

  if ( Params::Parameters().IsPeriodic() ) {
    fprintf(outfile,"$unit_cell\n%f %f %f\n%f %f %f\n$end\n",UnitCellAxes[0],
	    UnitCellAxes[1], UnitCellAxes[2], UnitCellAngles[0],
	    UnitCellAngles[1], UnitCellAngles[2]);
  }


  fprintf(outfile,"\n");

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


void Cluster::UpdateTrajectoryFile(int step, bool new_trajectory) {

  printf("Updating trajectory file\n");
  ofstream traj;
  string filename = "traj.xyz";
  int Ntot = GetTotalNumberOfAtoms();
  int periodic_factor = 1;
  if ( Params::Parameters().IsPeriodic() ) 
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
  if ( Params::Parameters().IsPeriodic() ) {
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
    printf("Cluster::CreatePeriodicImageDimerList() : Cannot open file '%s'\n",
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
	delete [] ImageMonomers;
      } // end loop over z
    } // end loop over y
  } // end loop over x


  
  fclose(xyz);

}


void Cluster::CreatePeriodicImageDimerList() {

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
  // The other dimension, N grows dynamically to accomodate as many entries
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

  */

  bool UseSymmetry = Params::Parameters().UseCrystalSymmetry();

  // Pass 1
  bool skip; // flag for ignoring some.
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
		if ( Params::Parameters().PrintLevel() ) {
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

  if ( Params::Parameters().PrintLevel() ) {
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
  //printf("Running: %s\n",cmd.c_str());
  system(cmd.c_str());
  



  // Pass 2: Now actually create the dimer list, using ImageList from
  // Pass 1 to generate the dimers

  DimerImages = new Dimer[keepit+1];

  int Mon_index = NMon;
  int idim;
  for (idim=1;idim<=keepit;idim++) {
    Mon_index++;
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

    // Create the new dimer between real_mon in central cell and image_mon
    // in periodic cell
    ImageMon.SetIndex(Mon_index); // renumber it  
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Monomer %d, K_vec = (%d,%d,%d)\n",Mon_index,x,y,z);
    DimerImages[idim].Initialize(Monomers[real_mon],ImageMon);
      
    // Link image monomer to one in central cell
    DimerImages[idim].SetReferenceMonomerIndex(image_mon);
    DimerImages[idim].SetImageCell(x,y,z);
    DimerImages[idim].SetSymmetryFactor(symfac);

    if (Params::Parameters().PrintLevel() > 0) 
      printf("(%d*=%d,%d) Dimer created.\n",image_mon,Mon_index,real_mon);
  }

  printf("%d image dimers created\n\n",idim-1);
  NDim_images = idim-1;

}

void Cluster::UpdateJobStatus(int ijob) {

  // define some helpful points

  
  int QM_mon_start = 0;
  int QM_dim_start = NMon;
  int MM_start = NMon + NDim + NDim_images;
  if ( Params::Parameters().DoQMBenchmark() )
    MM_start++;
  int N_MM_jobs = 1 + NMon + NDim + NDim_images; // Tinker
  if (Params::Parameters().GetMMType()==2)  
    N_MM_jobs = NMon; // AIFF
  int MM_end = MM_start + N_MM_jobs;

  if (ijob == QM_mon_start)
    printf("\nRunning %d QM Monomer jobs:\n",NMon);
 
  // insert some spaces in output for pretty printing
  // Break each category into groups of 5 and rows of 50

  // QM monomer jobs
  if (ijob < NMon && ijob%5==0 && ijob != QM_dim_start)
    printf(" ");

  if (ijob < NMon && ijob%50==0 && ijob > 0 && ijob != QM_dim_start-1 )
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
    printf("\nRunning %d QM Dimer jobs:\n",NDim+NDim_images);
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

  // define some helpful points
  int QM_mon_start = 0;
  int QM_dim_start = NMon;
  int MM_full = NMon + NDim + NDim_images;
  int MM_start = MM_full + 1;
  int MM_end = MM_start + NMon + NDim + NDim_images;

  if (ijob == QM_mon_start)
    printf("\nRunning %d QM Monomer jobs:\n",NMon);
 
  // insert some spaces in output for pretty printing
  // Break each category into groups of 5 and rows of 50

  // QM monomer jobs
  if (ijob < NMon && ijob%5==0 && ijob != QM_dim_start)
    printf(" ");

  if (ijob < NMon && ijob%50==0 && ijob > 0 && ijob != QM_dim_start-1 )
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


  double Eaiff=0.0, E_es=0.0, E_ind=0.0, E_es_pol=0.0, E_2b_disp=0.0, E_3b_disp=0.0;

  printf("\n -----------------------------------------------------------------\n");
  printf("  ** Compute ab initio force-field energy for the full system **\n");
  printf(" -----------------------------------------------------------------\n");

  fflush(stdout);

  if (Params::Parameters().IsPeriodic()) {
    // Compute periodic crystal induction/electrostatics


    if (do_es && do_ind && Params::Parameters().GetEwaldInductionType()==3) {
      // Old algorithm
      E_es_pol = ComputePeriodicMultipoleInteractions();
      E_es=E_es_pol; // hack to make sure the energy gets saved
    }

    else if (do_es && do_ind) {
      // Newer, more efficient matrix algorithms
      // First get the induced multipoles
      if (Params::Parameters().GetEwaldInductionType()==2 || Params::Parameters().GetEwaldInductionType()==4) {
	ComputePeriodicAIFFInducedMultipolesMadelung();
      }
      else if ( Params::Parameters().GetEwaldInductionType()==1) {
	ComputePeriodicAIFFInducedMultipoles(); // Find the induced moments
      }
      else{
	printf("Cluster::ComputeAIFFEnergy() ERROR: Unknown Ewald Induction type.\n"); exit(1);
      }

      // Now perform the Ewald sum
      // Compute how induction damping affects Ewald sum
      //double Eind_corr = ComputeEwaldDampingCorrection();
      ComputePeriodicAIFFElectrostaticsAndInduction();
      
      // Free up some memory
      if (Params::Parameters().GetEwaldInductionType()==4) {
	delete [] ReciprocalTabs;
	delete [] DirectTabs;
	delete [] DeltaDampTabs;
      }



      E_es = Lattice_E_Electrostatic_MM;
      E_ind = Lattice_E_Induction_MM;
      printf("New routines: E_es = %f, E_ind = %f\n",E_es*HartreesToKJpermole,E_ind*HartreesToKJpermole);
      E_es_pol = E_es + E_ind;
    }
    else if (do_es && !do_ind) {    
      E_es = ComputePeriodicAIFFElectrostatics(); // No induction
    }

    /*
    // Old, working, but slow routines
    if (do_es || do_ind) 
      E_es_pol = ComputePeriodicMultipoleInteractions();
    */

    // Compute the 2-body intermolecular dispersion energy for the
    // periodic system.
    if (do_2b_disp) E_2b_disp = ComputePeriodicTwoBodyDispersion();

    // Compute the 3-body intermolecular dispersion energy for the
    // periodic system.
    if (do_3b_disp) E_3b_disp = ComputeManyBodyDispersion();
  }

  else {
    // Compute the electrostatic + induction/polarization energy 
    printf("New Cluster ES/Ind algorithms:\n");
    if (do_es)
      E_es = ComputeClusterAIFFElectrostatics();

    if (do_ind)
      E_ind = ComputeClusterAIFFInduction();

    E_es_pol = E_es + E_ind;
    printf("New: E_es = %f, E_ind = %f, Sum = %f\n",E_es*HartreesToKJpermole,E_ind*HartreesToKJpermole,(E_es+E_ind)*HartreesToKJpermole);
    

    /* Old working, but slow routines
    printf("Old Cluster ES/Ind algorithms:\n");
    if (do_es || do_ind)
      E_es_pol = ComputeClusterMultipoleInteractions();  
    //E_es_pol = 0.0; printf("HACK: E_es_pol = 0.0\n");
    printf("Old: E_es_pol = %f\n",E_es_pol*HartreesToKJpermole);
    */


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
  Eaiff = E_es + E_ind + E_2b_disp + E_3b_disp;
  printf("E_es_pol = %f\n",E_es_pol*HartreesToKJpermole);
  printf("E_2b_disp = %f\n",E_2b_disp*HartreesToKJpermole);
  printf("E_3b_disp = %f\n",E_3b_disp*HartreesToKJpermole);
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

// Compute AIFF permanent multipolar electrostatics for full cluster.  Matrix algorithm.
double Cluster::ComputeClusterAIFFElectrostatics() {

  // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);

  // Note, we assume all the geometric interaction matrices Tab have
  // already been computed and are stored in the corresponding Dimer
  // objects.
  
  // This is a purely pairwise additive effect, so for aperiodic
  // systems, we just sum the 2-body permanent multipole interactions
  // we have already computed.  Since there is no many-body
  // contribution from the permanent multipoles, the only point in
  // computing this is for the case where we perform QM calculations
  // only among local pairs.  In those cases, we need classical
  // electrostatics at long ranges.

  double Ees = 0.0;
  printf("Compute classical electrostatic interactions for the full cluster.\n");
  for (int i=1;i<=NDim;i++) {
    Ees += Dimers[i].GetMMElectrostaticEnergy(); // in hartrees
  }
  
  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > -1)
      printf("Full cluster electrostatics wall time = %.5f seconds\n",elapsed_time);


  return Ees;
}

// Compute self-consistnet induction for the full cluster.
// Matrix-based algorithm.  See dimer version
// Dimer::ComputeAIFFInduction() for more thorough comments.  This
// algorithm stores the matrices fully in core.  May need more
// efficient version for dealing with really large clusters!
double Cluster::ComputeClusterAIFFInduction() {
  // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);

  double Eind = 0.0;
  printf("Computing self-consistent induction energy for the entire cluster.\n");

  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);

  // Create storage space for induced moments on all atoms.
  // In the process, also create some lists for array offsets used later on.  
  int Natoms = GetTotalNumberOfAtoms(); 
  int sizes[Natoms]; // gives dimensionality of the induced multipoles, etc for each atom
  int offset[Natoms]; // gives offset for given atom in matrices
  int Atom_count[NMon+1]; // gives atom offset for first atom of each new monomer
  int dim=0; // total size of the induced multipole array, etc.

  offset[0] = 0;
  //Atom_count[0] = 0; Atom_count[1] = 0;
  int iatom = 0;
  while (iatom<Natoms) {
    for (int imon=1;imon<=NMon;imon++) {
      Atom_count[imon] = iatom;
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	int AtomicNumber = Monomers[imon].GetAtom(iA).GetAtomicNumber();
	if (AtomicNumber == 1) {
	  sizes[iatom] = 3;
	}
	else {
	  sizes[iatom] = 8;
	}
	dim += sizes[iatom];
	if (iatom > 0) 
	  offset[iatom] = offset[iatom-1] + sizes[iatom-1];

	iatom++;
	//printf("atom %d, size = %d, offset = %d, current dim = %d\n",iatom-1,sizes[iatom-1],offset[iatom-1],dim);
      }
    }
  }

  int array_sizes = dim*dim + 3*dim; // Z(dim,dim), V0, V0_copy, dQ
  double req_memory = (double) array_sizes*8.0/(1024.0*1024.0);
  
  printf("------------------------------------------------------\n");
  printf("  Full Cluster Induction memory requirements:\n");
  printf("                    Number of sites = %d\n",Natoms);
  printf("     Number of inducible multipoles = %d\n",dim);
  printf("              Total memory required = %.2f MB\n",req_memory);
  printf("------------------------------------------------------\n");
  fflush(stdout);

  // Initialize storage for Z, dQ, and V0
  Vector V0(dim);
  Vector dQ(dim);
  Matrix Z(dim,dim);


  // Build the Z matrix.

  // Diagonal blocks first: inverse of polarizability.  Loop over
  // atoms on each monomer, invert polarizability, and put it in Z.

  // Loop over each atom on each monomer
  iatom=0;
  for (int imon=1;imon<=NMon;imon++) {
    for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {

      // Grab the non-zero block of the Polarizability matrix.  That means skipping the first row/col (charge-charge pol)
      // which is always zero.  Hence, the 1,1 offset in the GetBlock call.
      int Npols = sizes[iatom]; // dimensionality of the non-zero polarizabilities block.  
      Matrix PolMat = Monomers[imon].GetAtom(iA).GetPolarizability().GetPolarizabilities().GetBlock(Npols,Npols,1,1);
      PolMat.Inverse(); // Invert the polarizability.
      Z.SetBlock(PolMat,offset[iatom],offset[iatom]); // Set the block.

      iatom++;
    }
  }

  //Z.Print("Z with diagonal blocks");

  // Now set off-diagonal blocks between pairs of atoms in different monomers
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int imonB=imonA+1;imonB<=NMon;imonB++) {
      // Grab the index of the appropriate dimer -- need this to
      // obtain the Tab matrices
      int idimer = DimerLookup(imonA,imonB); 

      ////////////////// NEED TO GET NpolsA, NpolsB, offsetA, offsetB right.


      for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
	int indexA = Atom_count[imonA] + iA;
	int NpolsA = sizes[indexA];
	int offsetA = offset[indexA];
	//printf("imonA = %d, iA=%d, indexA=%d, NpolsA=%d, offsetA=%d\n",imonA,iA,indexA,NpolsA,offsetA);

	for (int iB=0; iB<Monomers[imonB].GetNumberOfAtoms(); iB++) {
	  int indexB = Atom_count[imonB] + iB;
	  int NpolsB = sizes[indexB];
	  int offsetB = offset[indexB];
	  //printf("imonB = %d, iB=%d, indexB=%d, NpolsB=%d, offsetB=%d\n",imonB,iB,indexB,NpolsB,offsetB);

	  Matrix Tab_block = Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB).GetBlock(NpolsA,NpolsB,1,1);
	  Z.SetBlock(Tab_block,offsetA,offsetB);
	  Tab_block.Transpose(); // Tba = Transpose(Tab)
	  Z.SetBlock(Tab_block,offsetB,offsetA);
	}
      }
    }
  }
  //Z.Print("Final Z matrix");

  // Build potential V0 due to permanent multipoles
  // For an atom "a" on MonA:
  // V0 = sum(atoms "b" on all other monomers MonB) Tab_tu * Qb_u
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {
      int indexA = Atom_count[imonA] + iA;
      int nQA = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetLength();
      Vector tmp(nQA);
      for (int imonB=1;imonB<=NMon;imonB++) {
	if (imonA != imonB) { // no intramolecular interactions

	  // Set logical flag for if (imonA < imonB).  Important
	  // because Dimer objects are always stored with A<B, and
	  // there is directionality implicit in the Tab interaction
	  // matrices: Tab for BA is the transpose of Tab for AB, and
	  // we need to be sure to grab the proper one.
	  bool AB_order = true;
	  if (imonA > imonB) AB_order = false; 

	  // Grab the index of the appropriate dimer -- need this to
	  // obtain the Tab matrices
	  int idimer = DimerLookup(imonA,imonB); 

	  for (int iB=0;iB<Monomers[imonB].GetNumberOfAtoms();iB++) {
	    Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());

	    // Grab Tab for the current pair of atoms
	    Matrix Tab;
	    if (AB_order) {
	      Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	    }
	    else {
	      Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
	      Tab.Transpose();
	    }

	    // Contract Tab with multipole moments on atom B: Tab*QB
	    tmp += Tab.MatrixTimesVector(QB.GetMoments());
	  }
	}
      }
      // Now grab the relevant elements of the tmp vector and store them in V0
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	V0[t+offsetA] = tmp[t+1];
      }
    }
  }

  // Now solve Z*dQ = -V0 for the induced moments dQ --> dQ =
  // Z^(-1)*V0.  

  // Can optionally switch to the expert equation solver, which is
  // numerically more stable and can provide error estimates, albeit
  // with roughly double the memory requirement.  I haven't yet found
  // cases where it matters.  Even in cases where the polarization catastrophe
  // occurs, the resulting induction energies don't seem to change.

  Vector V0_copy = V0; // create a copy of V0 which gets destroyed
		       // when we solve for dQ
  V0_copy.Scale(-1.0);

  Z.SolveLinearEquations(V0_copy); // Solve Z*dQ=-V0 (via LU decomposition)
  //Z.ExpertSolveLinearEquations(V0_copy); // Solve Z*dQ=-V0 (via LU decomposition)
  dQ = V0_copy;

  // Compute the induction energy from these induced moments: 
  // Eind = 0.5*dQA*V0
  Eind = 0.5*dQ.DotProduct(V0);
  //printf("Eind = %f\n",Eind*HartreesToKJpermole);



  // Transfer the sparse list of induced multipole moments to
  // non-sparse one.
  Multipole *dQall = new Multipole[Natoms];
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
      int rank = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetRank();
      int indexA = Atom_count[imonA] + iA;
      dQall[indexA].Initialize(rank);
      for (int t=0;t<sizes[indexA];t++) {
	dQall[indexA](t+1) = dQ[offset[indexA]+t];  // t+1 on LHS because always skip charge
      }
    }
  }
  // Optionally print out the induced multipole moments
  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles ***  New algorithm\n");
    int count = 0;
    for (int imon=1;imon<=NMon;imon++) {
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	char label[50];
	sprintf(label,"Monomer %d, %s%d",imon,
		Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		Monomers[imon].GetAtom(iA).GetAtomIndex());
	string str = label;
	dQall[count].Print(str);
	count++;
      }
    }
  }

  // Free up memory 
  delete [] dQall;

  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > -1)
      printf("Full cluster induction energy (matrix) wall time = %.5f seconds\n",elapsed_time);


  return Eind;
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

  }

  //Grad.Scale(HartreesToKJpermole);
  //Grad.PrintGradient("Electrostatic Gradient, full cluster");
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


    //Grad.Scale(HartreesToKJpermole);
    //PrintGradient("Full cluster MM gradient",Grad);
    //Grad.Scale(1.0/HartreesToKJpermole);
  return Grad;
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
		    if (!IsImageMonomer) 
		      Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
		    // Otherwise, we grab the Tab we just computed above
		    else {
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
	    if ( !IsImageMonomer) {
	      idimer = DimerLookup(imonA,imonB); 
	    }
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
  
  
    /* Step 1b: Compute the undamped induction energy in the short-range distance (~10 bohr or 6 Angstrom) 
       This allows us to assess the effect of damping on the induction energy in a finite cluster using
       the same set of self-consistently induced multipoles.  
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
  
  // GJB
  double Total_ERecip = 0.0;
  double Total_ERecip_kn0 = 0.0;

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
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      Multipole InduceQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments()); // convergent induced multipole for the induction energy 
      int NmomA = QA.GetLength();
      int InduceNmomA = InduceQA.GetLength();
      

      // GJB sum the two
      Multipole QA_tot(QA);
      QA_tot = QA + InduceQA;
      /*
      if (iA==0) {
	QA.Print("QA");
	InduceQA.Print("Induced QA");
	QA_tot.Print("Sum");
      }
      */

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  Multipole InduceQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	  int NmomB = QB.GetLength();
	  int InduceNmomB = InduceQB.GetLength();
          
	  // GJB sum the two
	  Multipole QB_tot(QB);
	  QB_tot = QB + InduceQB;

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
		  RecipTabAB_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
											       Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											       CellV, RecipCellx, RecipCelly, RecipCellz);

		  
		  
		  //printf("Mon %d, atom %d with Mon %d, atom %d\n",imonA,iA,imonB,iB);
		  //RecipTabAB_kn0.Print("RecipTabAB_kn0"); //GJB
		  
		  int dimTAB1_kn0 = RecipTabAB_kn0.GetRows();
		  int dimTAB2_kn0 = RecipTabAB_kn0.GetCols();
		  // Loop over multipole order
		  for (int t0=0;t0<min(NmomA,dimTAB1_kn0);t0++)  
		    for (int u0=0;u0<min(NmomB,dimTAB2_kn0);u0++) {
		      StaticERecip_kn0 += 0.50*QA(t0)*RecipTabAB_kn0(t0,u0)*QB(u0); // lattice electrostatic energy for permanent multipoles
		      undampedInduceERecip1_kn0 += 0.50*InduceQA(t0)*RecipTabAB_kn0(t0,u0)*QB(u0); // induction contribution, undamped 
		      Total_ERecip_kn0 += 0.50*QA_tot(t0)*RecipTabAB_kn0(t0,u0)*QB_tot(u0); // GJB: lattice electrostatic energy for total multipoles
		    } 
	
		} // endif |kn|=0
		
		else { // |kn|!=0
		  Matrix RecipTabAB;
		  RecipTabAB = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
										       Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
										       kx, ky, kz, CellV,
										       RecipCellx, RecipCelly, RecipCellz,-999.0);

		  int dimTAB1 = RecipTabAB.GetRows();
		  int dimTAB2 = RecipTabAB.GetCols();
		  for (int t=0;t<min(NmomA,dimTAB1);t++)  
		    for (int u=0;u<min(NmomB,dimTAB2);u++) {

		      double dE_es = 0.50*QA(t)*RecipTabAB(t,u)*QB(u);
		      double dE_ind = 0.50*InduceQA(t)*RecipTabAB(t,u)*QB(u);

		      StaticERecip += dE_es; // lattice electrostatic energy for permanent multipoles
		      undampedInduceERecip1 += dE_ind; // induction contribution
		      Total_ERecip += 0.50*QA_tot(t)*RecipTabAB(t,u)*QB_tot(u); // GJB: Total contribution

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
  double Total_EDirec = 0.0;

  printf("  Begin the Ewald summation in direct space\n");  
  //printf("  --nX=%d,nY=%d,nZ=%d\n",nX,nY,nZ);
  
  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];
  
  
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
      for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
	Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
	Multipole InduceQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments());
	int NmomA = QA.GetLength();
	int InduceNmomA = InduceQA.GetLength();
	

	// GJB sum the permanent and induced multipoles
	Multipole QA_tot(QA);
	QA_tot = QA + InduceQA;
	/*
	if (iA==0) {
	  QA.Print("QA");
	  InduceQA.Print("Induced QA");
	  QA_tot.Print("Sum");
	}
	*/
	for (int imonB = 1; imonB <= NMon; imonB++) {
	  int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	  Vector RotVecB(Monomers[imonB].GetRotationVector());
	  double RotAngB = Monomers[imonB].GetRotationAngle();
	  
	  for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	    Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	    Multipole InduceQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    int NmomB = QB.GetLength();
	    int InduceNmomB = InduceQB.GetLength();
            

	    // GJB sum the two
	    Multipole QB_tot(QB);
	    QB_tot = QB + InduceQB;
	    /*
	    if (iB==0) {
	      QB.Print("QB");
	      InduceQB.Print("Induced QB");
	      QB_tot.Print("Sum");
	    }
	    */
	    for (int nx = -nX; nx<=nX; nx++){
	      for (int ny = -nY; ny<=nY; ny++){
                for (int nz = -nZ; nz<=nZ; nz++){
		  
                  if (nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB);
		  //{               
		  //printf("Skipping self-interaction real-space term\n");
		  //} 
		  
                  else { // |rAB-rn|!=0
		    Matrix DirecTabAB;
		    DirecTabAB = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
											 Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											 nx, ny, nz, CellV,
											 UnitCellx, UnitCelly, UnitCellz, -999.0);
		    
		    int DdimTAB1 = DirecTabAB.GetRows();
		    int DdimTAB2 = DirecTabAB.GetCols();
		    
		    for (int t5=0;t5<min(NmomA,DdimTAB1);t5++) 
		      for (int u5=0;u5<min(NmomB,DdimTAB2);u5++) {
			
			double dE_es = 0.50*QA(t5)*DirecTabAB(t5,u5)*QB(u5);
			double dE_ind = 0.50*InduceQA(t5)*DirecTabAB(t5,u5)*QB(u5);

			StaticEDirec += dE_es;
			undampedInduceEDirec1 += dE_ind;
			
			Total_EDirec += 0.50*QA_tot(t5)*DirecTabAB(t5,u5)*QB_tot(u5);;

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
  
  double Total_Uself = 0.0;
  double Total_UIntramol = 0.0;

  
  // constant for the calculation Uself
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  double kappa = kappa_param/V_V;
  double alphaa = kappa*kappa;
  double ssff = kappa/sqrt(3.14159265359);
  printf("  -----alphaa = %12.6f, kappa = %12.6f, ssff = %12.6f\n", alphaa,kappa,ssff);
  
  double Perm=4.0*pi*epsilon*1000.0/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  printf("pi = %f, epsilon = %g, MtoAng=%g, AngToBohr=%g, ec=%g, Na=%g, HartreesToKJpermole=%g\n",pi,
	 epsilon,MtoAng,AngToBohr,ec,Na,HartreesToKJpermole);
  printf("--Perm = %12.6f\n", Perm);
  
  
    
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
    Vector RotVec(Monomers[Imon].GetRotationVector());
    double RotAng = Monomers[Imon].GetRotationAngle();
    for (int IA=0;IA<NAtoms;IA++)  {// loop over atoms on monomer 
      Multipole QA1(Monomers[Imon].GetAtom(IA).GetMultipoleMoments()); 
      Multipole InduceQA1(Monomers[Imon].GetAtom(IA).GetInduceMultipoleMoments());
      int NmomA1 = QA1.GetLength();
      int InduceNmomA1 = InduceQA1.GetLength();

      Multipole QA1_tot(QA1);
      QA1_tot = QA1 + InduceQA1;

      for (int IB=0;IB<NAtoms;IB++)  {// loop over atoms on monomer 
	Multipole QB1(Monomers[Imon].GetAtom(IB).GetMultipoleMoments()); 
	Multipole InduceQB1(Monomers[Imon].GetAtom(IB).GetInduceMultipoleMoments());
	int NmomB1 = QB1.GetLength();
	int InduceNmomB1 = InduceQB1.GetLength();
	
	// GJB sum the two
	Multipole QB1_tot(QB1);
	QB1_tot = QB1 + InduceQB1;

	// if IA==IB, do the self-energy correction
	if (IA==IB){
	  //printf("----calculate the self-energy----\n");
	  Uself_Electrostatic -= ssff*QA1(0)*QB1(0)/Perm;
	  Uself_induce1 -= ssff*InduceQA1(0)*QB1(0)/Perm;
	  Total_Uself -= ssff*QA1_tot(0)*QB1_tot(0)/Perm;
	} 
	// else do the intramolecular-energy correction for the case that A and B are atoms in the same molecule in the central cell 
	else {
	  //printf("----calculate the intramolecular-energy----\n");
	  
	  
	  Vector RotVecA = RotVec;
	  double RotAngA = RotAng;
	  Vector RotVecB = RotVec;
	  double RotAngB = RotAng;
	  
	  
	  Matrix DirecTabA1B1;
	  //Matrix DampedDirecTabA1B1; 
	  
          
	  //this interaction Tab matrix is for the case where A and B
	  //are atoms in the same molecule the explicit form of the
	  //interaction energy between A and B, so we use
	  //BuildDirecInteractionMatrix (not reciprocal or direct
	  //matrix)
	  DirecTabA1B1 = Monomers[Imon].GetAtom(IA).BuildInteractionMatrix(RotVecA, RotAngA,
									   Monomers[Imon].GetAtom(IB), RotVecB, RotAngB,-999.0);
	  
	  //DampedDirecTabA1B1 = Monomers[Imon].GetAtom(IA).BuildInteractionMatrix(RotVecA, RotAngA,
	  //                                                               Monomers[Imon].GetAtom(IB), RotVecB, RotAngB, beta_damp);
	  
	  
	  int DdimTA1B1 = DirecTabA1B1.GetRows();
	  int DdimTA1B2 = DirecTabA1B1.GetCols();
	  
	  for (int tt5=0;tt5<min(NmomA1,DdimTA1B1);tt5++) // for permanent multipoles
	    for (int uu5=0;uu5<min(NmomB1,DdimTA1B2);uu5++) {
	      UIntramolecule_Electrostatic -= 0.50*QA1(tt5)*DirecTabA1B1(tt5,uu5)*QB1(uu5);
	      //UIntramolecule_induce1 -= 0.50*InduceQA1(tt5)*DampedDirecTabA1B1(tt5,uu5)*QB1(uu5);
	      undampedUIntramolecule_induce1 -= 0.50*InduceQA1(tt5)*DirecTabA1B1(tt5,uu5)*QB1(uu5);
	      Total_UIntramol -= 0.50*QA1_tot(tt5)*DirecTabA1B1(tt5,uu5)*QB1_tot(uu5);

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

  
  double Esurf = 0.0; // Surface term = 0 if tin-foil boundary conditions
  if (! Params::Parameters().TinFoilBoundaryConditions()) { // Vacuum Boundary conditions
    Vector Qsurf(3);
    Vector Qsurf_ind(3);

    for (int Imon = 1; Imon <= NMon; Imon++) {
      int NAtoms = Monomers[Imon].GetNumberOfAtoms();
      for (int IA=0;IA<NAtoms;IA++)  {// loop over atoms on monomer 
	Multipole QA(Monomers[Imon].GetAtom(IA).GetMultipoleMoments()); 
	Multipole InduceQA(Monomers[Imon].GetAtom(IA).GetInduceMultipoleMoments());

	// Permanent contribution
	Vector pos = Monomers[Imon].GetAtom(IA).GetPosition();
	pos.Scale(QA(0)*AngToBohr); // 00 component times position, q*r
	Qsurf += pos; // 00 component times position, q*r
	Qsurf[0] += QA(1); // 11c component
	Qsurf[1] += QA(2); // 11s component
	Qsurf[2] += QA(3); // 10 component

	// Induced contribution - No induced charge... only dipole.
	Qsurf_ind[0] += InduceQA(1); // 11c component
	Qsurf_ind[1] += InduceQA(2); // 11s component
	Qsurf_ind[2] += InduceQA(3); // 10 component
      }
    }
    Qsurf.Print("Qsurf");
    Qsurf_ind.Print("Qsurf_ind");

    // Esurf = 2*pi/((2*eps + 1)V ) * | sum_A Q_00*rA + Q_1? | ^2
    // where Q_1? means dipole vector, rA is position vector for atom.

    double Esurf_perm = 2*pi/(3*V*Perm)*Qsurf.DotProduct(Qsurf);

    Qsurf += Qsurf_ind; // Combine the permanent and induced portions.
    Esurf = 2*pi/(3*V*Perm)*Qsurf.DotProduct(Qsurf);
    printf("Surface dipole contribution under vacuum boundary conditions = %f kJ/mol\n",Esurf*HartreesToKJpermole);
    printf("      Without induction = %f kJ/mol\n",Esurf_perm*HartreesToKJpermole);
  }
  
  
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
  
  
  printf("\nGJB: Totals\n");
  printf("     EDirect     = %12.4f kJ/mol\n",Total_EDirec*HartreesToKJpermole);
  printf("     ERecip      = %12.4f kJ/mol\n",Total_ERecip*HartreesToKJpermole);
  printf("     ERecipkn0   = %12.4f kJ/mol\n",Total_ERecip_kn0*HartreesToKJpermole);
  printf("     UIntramol   = %12.4f kJ/mol\n",Total_UIntramol*HartreesToKJpermole);
  printf("     Uself       = %12.4f kJ/mol\n",Total_Uself*HartreesToKJpermole);

  
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
    many dimers we actually need.  Then, we do the real part in Pass 2
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
		if ( Params::Parameters().PrintLevel() ) {
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
    if (Params::Parameters().PrintLevel() > 0) 
      printf("K_vec = (%d,%d,%d)\n",x,y,z);
    
    // Renumber it & set reference monomer index
    ImageMon.SetReferenceMonomerIndex(ImageMon.GetIndex());
    ImageMon.SetIndex(Mon_index); 
    
    MonomerImages[imon] = ImageMon;
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


void Cluster::ComputeHarmonicFrequencies(Matrix& Hessian) {

  printf("\nComputing harmonic vibrational frequencies.\n");

  // Mass weight the Hessian
  Matrix mwHess = ComputeMassWeightedHessian(Hessian);

  // Project out the translational & rotational modes
  printf("Projecting out translational and rotational modes\n");
  Matrix Hess;
  if (! Params::Parameters().IsPeriodic() ) {
    Hess = ProjectOutTranslationAndRotation(mwHess);
  }
  else {
    Hess = mwHess;
  }

  // Diagonalize the mass-weighted Hessian
  Vector freqs = Hess.Diagonalize();
  int Nfreq = freqs.GetLength();

  // Zero out any really tiny frequencies
  for (int i=0; i<Nfreq;i++) {
    if ( fabs(freqs[i]) < 1.0e-6) {
      freqs[i] = 0.0;
    }
  }

  // Define some constants used in frequency unit conversions
  double c1 = 5.89141e-7;
  double c2 = 3.94562;

  // Count number of negative eigenvalues (imaginary freqs)
  int Nimag = 0;
  for (int i=0; i<Nfreq;i++) {
    if (freqs[i] < 0.0) {
      Nimag++;
      double imag = sqrt(fabs(freqs[i])/c1)*c2;
      printf("Imaginary frequency: %.2f\n",imag);
      freqs[i] = 0.0;
    }
  }
  printf("%d imaginary frequencies found\n\n",Nimag);

  // Convert the frequencies to wavenumbers
  for (int i=0; i<Nfreq;i++) {
    freqs[i] = sqrt(fabs(freqs[i])/c1)*c2;
  }
  freqs.SortLargestFirst();
  int Nonzero = freqs.Nonzero();
  printf("Harmonic Vibrational frequencies (cm-1): %d real frequencies\n",Nonzero);
  for (int i=0;i<Nonzero/6+1;i++) {
    for (int j=0;j<6;j++) {
      if ((6*i+j) < Nonzero && freqs[6*i+j]>0.0) {
	printf("  %8.2f",freqs[6*i+j]);
      }
      else
	break;
    }
    printf("\n");
  }

  // Get the zero-point energy
  double ZPE = ComputeZeroPointEnergy(freqs);
}

Matrix Cluster::ComputeMassWeightedHessian(Matrix& Hessian) {

  //int Natoms = GetTotalNumberOfAtoms();
  int Natoms = Hessian.GetCols()/3;
  Matrix mwHess = Hessian;
 
  for (int i=0;i<Natoms;i++) {
    double m1 = AtomicMasses[i];
    for (int j=0; j<Natoms;j++) {
      double m2 = AtomicMasses[j];
      double Gij = 1.0/sqrt(m1*m2);

      for (int p=3*i;p<3*(i+1);p++)
	for (int q=3*j;q<3*(j+1);q++) {
	  mwHess(p,q) *= Gij;
	}
    }
  }
  return mwHess;
}

// Projects out translation & rotation from the mass-weighted hessian.
Matrix Cluster::ProjectOutTranslationAndRotation(Matrix& Hessian) {

  //int Natoms = GetTotalNumberOfAtoms();
  int Natoms = Hessian.GetCols()/3;

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
  Vector Va = Inertia.GetColumnVector(0);
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
  tmp.Normalize();
  Ra.SetColumnVector(tmp,0);

  tmp = Rb.GetColumnVector(0);
  tmp.Normalize();
  Rb.SetColumnVector(tmp,0);

  tmp = Rc.GetColumnVector(0);
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

  /* Create projector that leaves only pure vibration: 
     Pvib = I - Ptrans - Prot;
  */
  Matrix Pvib(3*Natoms,true);
  Pvib -= Ptrans;
  Pvib -= Prot;

  // Apply the projector to project out translation & rotation:  Hvib = P'HP
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

// Returns the vibrational ZPE contribution, in kJ/mol
double Cluster::ComputeZeroPointEnergy(Vector freqs) {

  double zpe = 0.0;
  // ZPE = sum (h*c*nu_i), along with unit conversions to kJ/mol
  for (int i=0;i<freqs.GetLength();i++) {
    if (freqs[i] > 0)
      zpe += 0.5*h_planck*freqs[i]*c_light*100.0*Na/1000.0;
  }
  printf("Zero-point vibrational energy = %8.3f kJ/mol\n",zpe);
  return zpe;
}

// Compute Axilrod-Teller-Muto atomic 3-body dispersion - between all triplets
// of monomers.  This is the function used by the AIFF.
double Cluster::ComputeThreeBodyDispersion(string type) {
  time_t start_time, stop_time;
  start_time = time(NULL);

  double E3b_disp = 0.0;


  // A few lines for saving geometries of dominant trimers to disk, printing contributions
  // Threshhold for printing out significant contributions
  double print_thresh = Params::Parameters().GetThreeBodyDispersionPrintThreshold();
  int itrimer = 0;

  // Open file for printing
  FILE *trimer_file;
  if (print_thresh < 100000.0) {
    string filename = "big_trimers.xyz";
    if ((trimer_file = fopen(filename.c_str(),"w"))==NULL) {
      printf("Cluster::ComputeThreeBodyDispersion() : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
  }
  
  // Now get to the real work

  // Loop over triplets of monomers (trimers)
  for (int imon=1; imon<=NMon; imon++) {
    int NatomsI = GetNumberOfAtoms(imon);
    for (int jmon=imon+1; jmon<=NMon; jmon++) {
      int NatomsJ = GetNumberOfAtoms(jmon);
      for (int kmon=jmon+1; kmon<=NMon; kmon++) {
	int NatomsK = GetNumberOfAtoms(kmon);

	double trimer_contrib = 0.0;
	
	// Loop over individual atoms in the trimers
	for (int i=0;i<NatomsI;i++) {
	  Atom AtomI = Monomers[imon].GetAtom(i);
	  for (int j=0;j<NatomsJ;j++) {
	    Atom AtomJ = Monomers[jmon].GetAtom(j);
	    for (int k=0;k<NatomsK;k++) {
	      Atom AtomK = Monomers[kmon].GetAtom(k);
	      // Because we sum over only unique triplets, get rid of
	      // 1/6 factor.
	      double contrib = ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	      E3b_disp += contrib;

	      // Also accumulate contribution from this trimer
	      trimer_contrib += contrib;
	    }
	  }
	}
	// Check trimer contribution
	if (fabs(trimer_contrib) > print_thresh) {
	  printf("Trimer (%d, %d, %d) contrib = %.4f kJ/mol\n",
		 imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	  // Write geometry to XYZ file for later viewing
	  itrimer++;
	  fprintf(trimer_file,"%d\n%d Trimer (%d, %d, %d) contrib = %.4f kJ/mol\n",
		  NatomsI+NatomsJ+NatomsK,itrimer,imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	  Monomers[imon].PrintMonomerCartesian(trimer_file);
	  Monomers[jmon].PrintMonomerCartesian(trimer_file);
	  Monomers[kmon].PrintMonomerCartesian(trimer_file);

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

	  double trimer_contrib = 0.0;

	  // Loop over individual atoms in the trimers
	  for (int i=0;i<NatomsI;i++) {
	    Atom AtomI = Monomers[imon].GetAtom(i);
	    for (int j=0;j<NatomsJ;j++) {
	      Atom AtomJ = Monomers[jmon].GetAtom(j);
	      for (int k=0;k<NatomsK;k++) {
		Atom AtomK = MonomerImages[kmon].GetAtom(k);
		
		double contrib = (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
		E3b_disp += contrib;
		trimer_contrib += contrib;
		case1 += contrib;
		//E3b_disp += (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	      }
	    }
	  }
	  // Check trimer contribution
	  if (fabs(trimer_contrib) > print_thresh) {
	    printf("Trimer (%d, %d, %d) contrib = %.4f\n",
		   imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	    // Write geometry to XYZ file for later viewing
	    itrimer++;
	    fprintf(trimer_file,"%d\n%d Trimer (%d, %d, %d) contrib = %.4f kJ/mol\n",
		    NatomsI+NatomsJ+NatomsK,itrimer,imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	    Monomers[imon].PrintMonomerCartesian(trimer_file);
	    Monomers[jmon].PrintMonomerCartesian(trimer_file);
	    MonomerImages[kmon].PrintMonomerCartesian(trimer_file);
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

	  double trimer_contrib = 0.0;

	  // Loop over individual atoms in the trimers
	  for (int i=0;i<NatomsI;i++) {
	    Atom AtomI = Monomers[imon].GetAtom(i);
	    for (int j=0;j<NatomsJ;j++) {
	      Atom AtomJ = MonomerImages[jmon].GetAtom(j);
	      for (int k=0;k<NatomsK;k++) {
		Atom AtomK = MonomerImages[kmon].GetAtom(k);

		double contrib = (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
		E3b_disp += contrib;
		trimer_contrib += contrib;
		case2 += contrib;
		//E3b_disp += (1.0/3.0)*ComputeAxilrodTellerMutoDispersionTerm(AtomI,AtomJ,AtomK, type);
	      }
	    }
	  }
	  // Check trimer contribution
	  if (fabs(trimer_contrib) > print_thresh) {
	    printf("Trimer (%d, %d, %d) contrib = %.4f\n",
		   imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	    // Write geometry to XYZ file for later viewing
	    itrimer++;
	    fprintf(trimer_file,"%d\n%d Trimer (%d, %d, %d) contrib = %.4f kJ/mol\n",
		    NatomsI+NatomsJ+NatomsK,itrimer,imon,jmon,kmon,trimer_contrib*HartreesToKJpermole);
	    Monomers[imon].PrintMonomerCartesian(trimer_file);
	    MonomerImages[jmon].PrintMonomerCartesian(trimer_file);
	    MonomerImages[kmon].PrintMonomerCartesian(trimer_file);
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

  if (print_thresh < 100000.0) fclose(trimer_file);

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

  // Compute the short-range damping function
  double damping = 1.0; // no damping
  if (Params::Parameters().DoDispersionDamping()) {
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
    damping = F6ij * F6ik * F6jk;
  }

  // ATM 3-body dispersion contribution
  double Eatm = damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
    (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
  
  //printf("Rij = %f, Rik = %f, Rjk = %f   E3 = %17.12f kcal/mol\n",Rij,Rik,Rjk,Eatm*HartreesToKcalpermole);

  return Eatm;
}

//GJB exploratory
void Cluster::ReadAtomicInductionDampingFactors(ifstream& infile) {

  string line;
  // Rewind the file, just in case
  Rewind(infile);

  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,18)=="$induction_damping") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) == "$end") {
	  break;
	}
	else if ( !line.empty() ) {
	  istringstream iss(line);
	  int atomic_number;
	  double damping_factor;
	  iss >> atomic_number;
	  iss >> damping_factor;
	  // Set the atom_Type
	  Params::Parameters().SetAtomicInductionDampingFactor(atomic_number,damping_factor);
	}
      }
    }
  }
  infile.clear();

  printf("Atomic Induction Damping factors:\n");
  for (int i=1;i<93;i++) {
    printf("%2d   %.2f\n",i,Params::Parameters().GetAtomicInductionDampingFactor(i));
  }


  
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
		if ( Params::Parameters().PrintLevel() ) {
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

  if ( !Params::Parameters().DoDispersionDamping() ) {
    printf("  Monomers: %10.4f kJ/mol (undamped)\n",monomer_3b*HartreesToKJpermole);
    printf("    Dimers: %10.4f kJ/mol (undamped)\n",dimer_3b*HartreesToKJpermole);
    printf("   Trimers: %10.4f kJ/mol (undamped)\n",trimer_3b*HartreesToKJpermole);
    printf("--------------------------------\n");
    printf("     Total: %10.4f kJ/mol (undamped)\n",sum*HartreesToKJpermole);


  }
  else {
    printf("  Monomers: %10.4f kJ/mol\n",monomer_3b*HartreesToKJpermole);
    printf("    Dimers: %10.4f kJ/mol\n",dimer_3b*HartreesToKJpermole);
    printf("   Trimers: %10.4f kJ/mol\n",trimer_3b*HartreesToKJpermole);
    printf("--------------------------------\n");
    printf("     Total: %10.4f kJ/mol\n",sum*HartreesToKJpermole);
  }
  printf("--------------------------------\n");

  printf("Note: Only trimer contribution is added to HMBI energy\n\n");
  E_3body_disp = trimer_3b;

  return E_3body_disp;  
}


