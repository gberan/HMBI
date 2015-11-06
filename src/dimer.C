#include "dimer.h"
#include "constants.h"
using namespace hmbi_constants;

Dimer::Dimer() : Tabs(NULL), DampedTabs(NULL), TabsGrad(NULL), 
		 DampedTabsGrad(NULL) {

 

  MM_Grad_Init = 0;
  QM_Grad_Init = 0;

}

Dimer::~Dimer() {

  // Special clean up if using the AIFF
  if ( Params::Parameters().GetMMType()==2 ) { // AIFF
    delete [] Tabs;
    delete [] DampedTabs;

    if ( MM_Grad_Init ) {
      delete [] TabsGrad;
      delete [] DampedTabsGrad;
    }
  }
  
}

void Dimer::Initialize(Monomer& m1, Monomer& m2) {

  // Assign monomers such that index(A) < index(B)
  if ( m1.GetIndex() < m2.GetIndex() ) {
    MonA = m1; MonB = m2; 
  }
  else {
    MonA = m2; MonB = m1;
  }

  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex();

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();

  // Set spin state.  Crude algorithm for determining dimer spin state
  // from monomers:
  // 2 singlets -> singlet
  // singlet + other -> other
  // non-singlet + non-singlet -> singlet if even number of electrons
  //      e.g. 2 triplets: 3 + 3 = 6 -> singlet
  // non-singlet + non-singlet -> doublet if odd number of electrons
  //      e.g. doublet + triplet = 2 + 3 = 5 -> doublet
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else if (spinA>1 && spinB>1 && (spinA+spinB)%2==0)  
    spin = 1;
  else if (spinA>1 && spinB>1 && (spinA+spinB)%2==1)
    spin = 2; 
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  QM_Grad_Init = 0;
  MM_Grad_Init = 0;

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }

  // PBC features - 0 by default
  is_image = false;
  reference_MonB = 0;
  K_vector[0] = 0;
  K_vector[1] = 0;
  K_vector[2] = 0;

  // No extra symmetry by default
  sym_fac = 1;

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];
  //printf("Separation = %f (based on atoms %d and %d)\n",Separation,minA,minB);
}

// This function updates the Monomers, along with all members derived
// from those, in case those objects have changed.  Ideally would do
// this more dynamically using pointers...
void Dimer::UpdateObjects(Monomer& m1, Monomer& m2) {

  MonA = m1; 
  MonB = m2; 

  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex();  

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();

  // Set spin state.  Right now, it can handle 2 singlets, or singlet + anything.
  // Need something smarter to handle e.g. for 2 triplets
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];
}

void Dimer::UpdateObjects(Monomer& m1, Monomer& m2, Monomer& ref_m2) {

  MonA = m1; 
  MonB = m2; 

  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex();
  int ref_indexB = ref_m2.GetIndex();

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();

  // Set spin state.  Right now, it can handle 2 singlets, or singlet + anything.
  // Need something smarter to handle e.g. for 2 triplets
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];

  // Get energies/gradients from non-image monomer
  MonB.SetQMEnergy( ref_m2.GetQMEnergy() );
  MonB.SetMMEnergy( ref_m2.GetMMEnergy() );

  if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    //printf("Dimer::UpdateObjects(): setting gradients from reference monomer %d\n",ref_m2.GetIndex()) ;
    MonB.SetQMGradient( ref_m2 );
    MonB.SetMMGradient( ref_m2 );
  }
    

  // If using the AIFF, need to get FF parameters
  if (Params::Parameters().GetMMType()==2) { 

    fflush(stdout);
    // Copy the local coordinates
    MonB.SetRotationAngle( ref_m2.GetRotationAngle() );
    double vec[3];
    for (int i=0;i<3;i++) 
      MonB.SetRotationVector(i,ref_m2.GetRotationVector()[i]);

    // Copy the multipole moments & polarizabilities
    int NatomsB = MonB.GetNumberOfAtoms();
    int Natoms_refB = ref_m2.GetNumberOfAtoms();

    if (NatomsB != Natoms_refB) {
      printf("Dimer::UpdateObjects: Image monomer and its reference have different numbers of atoms: %d and %d\n",NatomsB,Natoms_refB);
    }

    for (int iatom=0;iatom<NatomsB; iatom++) {    
      if (ref_m2.GetAtom(iatom).MultipolesAreInitialized() && ref_m2.GetAtom(iatom).PolarizabilityIsInitialized() ) {
	// Copy multipole moments, polarizabilities, and dispersion coefficients
	MonB.GetAtom(iatom).SetMultipoleMoments( ref_m2.GetAtom(iatom).GetMultipoleMoments() );
	MonB.GetAtom(iatom).SetPolarizability( ref_m2.GetAtom(iatom).GetPolarizability() );

	double iso_pol = ref_m2.GetAtom(iatom).GetIsotropicDipolePolarizability();
	MonB.GetAtom(iatom).SetIsotropicDipolePolarizability(iso_pol);

	double C6 = ref_m2.GetAtom(iatom).GetC6Coefficient();
	double C8 = ref_m2.GetAtom(iatom).GetC8Coefficient();
	double C10 = ref_m2.GetAtom(iatom).GetC10Coefficient();
	MonB.GetAtom(iatom).SetDispersionCoefficients( C6, C8, C10);
	MonB.GetAtom(iatom).SetUCHFDispersionCoefficients( C6, C8, C10);
	
	MonB.GetAtom(iatom).SetDispersionAtomType(ref_m2.GetAtom(iatom).GetDispersionAtomType());

	// Copy Frequency-dependent polarizabilities
	MonB.GetAtom(iatom).SetFreq_Pol_Dipole(ref_m2.GetAtom(iatom).GetFreq_Pol_Dipole());
	MonB.GetAtom(iatom).SetFreq_Pol_Quad(ref_m2.GetAtom(iatom).GetFreq_Pol_Quad());


	//printf("iatom = %d\n",iatom);
	//MonB.GetAtom(iatom).GetMultipoleMoments().Print("Multipole moments\n");
	//MonB.GetAtom(iatom).GetPolarizability().Print("Polarizability\n");
	
      }
    }
  }
}


void Dimer::CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job) {
  
  // Set up the filename, with the full path.  File is e.g. 'm1.force'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  filename += ".in";

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateQChemJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  // Print comment section
  fprintf(job,"$comment\nDimer (%d,%d)\n$end\n\n",indexA,indexB);

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

  fprintf(job,"%s\n",rem.c_str());  
  fprintf(job,"%s\n",basis.c_str());  

/*
  // for CP correction, insert jobtype = BSSE at end of rem section.  
  // Note: This will override any previous setting for jobtype.
  if ( Params::Parameters().DoCounterpoise() ) {
    int spot = rem.find("$end");
    string bsse = "jobtype = bsse\n";
    rem.insert(spot,bsse);
    //printf("New rem:\n%s\n",rem.c_str());
  }
*/

  if ( Params::Parameters().DoCounterpoise() ) {  //added on 9th april, 2010 -> to perform bsse correction manually 

//monomer job at monomer basis is not required to prevent doing the same job twice: this job is same as m*.out
/*    fprintf(job, "@@@\n\n");
    fprintf(job,"$comment\nMonomer %d at Monomer Basis\n$end\n\n",indexA);
    fprintf(job,"$molecule\n%d %d\n",MonA.GetChargeState(), MonA.GetSpinState());
    MonA.PrintMonomerCartesian(job);
    fprintf(job,"$end\n");           
    // Print $rem section
    if (MM_job)
      rem = Params::Parameters().GetQChemRem2();
    else
      rem = Params::Parameters().GetQChemRem();
    fprintf(job,"%s\n",rem.c_str());  
*/

    fprintf(job, "@@@\n\n");
    fprintf(job,"$comment\nMonomer %d at Dimer Basis\n$end\n\n",indexA);
    fprintf(job,"$molecule\n%d %d\n",MonA.GetChargeState(), MonA.GetSpinState());    
    MonA.PrintMonomerCartesian(job);
    MonB.PrintGhostMonomerCartesian(job);
    fprintf(job,"$end\n");
    // Print $rem section
    if (MM_job)
      rem = Params::Parameters().GetQChemRem2();
    else 
      rem = Params::Parameters().GetQChemRem();
    fprintf(job,"%s\n",rem.c_str());  

    fprintf(job,"%s\n",basis.c_str());  

//monomer job at monomer basis is not required to prevent doing the same job twice: this job is same as m*.out
/*
    fprintf(job, "@@@\n\n");
    fprintf(job,"$comment\nMonomer %d at Monomer Basis\n$end\n\n",indexB);
    fprintf(job,"$molecule\n%d %d\n",MonB.GetChargeState(), MonB.GetSpinState());    
    MonB.PrintMonomerCartesian(job);
    fprintf(job,"$end\n");
    // Print $rem section
    if (MM_job)
      rem = Params::Parameters().GetQChemRem2();
    else 
      rem = Params::Parameters().GetQChemRem();
    fprintf(job,"%s\n",rem.c_str());  
*/

    fprintf(job, "@@@\n\n");
    fprintf(job,"$comment\nMonomer %d at Dimer Basis\n$end\n\n",indexB);
    fprintf(job,"$molecule\n%d %d\n",MonB.GetChargeState(), MonB.GetSpinState());    
    MonA.PrintGhostMonomerCartesian(job);
    MonB.PrintMonomerCartesian(job);
    fprintf(job,"$end\n");
    // Print $rem section
    if (MM_job)
      rem = Params::Parameters().GetQChemRem2();
    else 
      rem = Params::Parameters().GetQChemRem();
    fprintf(job,"%s\n",rem.c_str());  

    fprintf(job,"%s\n",basis.c_str());  

  }

//  fprintf(job,"%s\n",rem.c_str());

  // Optionally print $external_charges section  -- doesn't work with CP correction
  if ( Params::Parameters().UseEmbeddingCharges() ) {
    if ( Params::Parameters().DoCounterpoise() ) {
      printf("ERROR: QChem external charges don't currently work with Counterpoise correction\n");
      exit(1);
    }
    fprintf(job,"$external_charges\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB) {
	Monomers[i].PrintQChemEmbeddingCharges(job);
      }
    }
    fprintf(job,"$end\n");
  }

  fclose(job);
}

// Wrapper for creating MM jobs
void Dimer::CreateMMJob(Monomer Monomers[], int NMon) {
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
    CreateTinkerJob(Monomers, NMon);
  }
  else if ( Params::Parameters().GetMMType()==2 ) { // AIFF
    // Do nothing.  Dimers get handled internally with
    // Dimer::ComputeMultipoleInteractions();
  }
  else if ( Params::Parameters().GetMMType()==3 ) { //QChem
    CreateQChemJob(Monomers, NMon, true);
  }

  else {
    printf("Monomer::CreateMMJob: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  }

}

// Creates Tinker MM job
void Dimer::CreateTinkerJob(Monomer Monomers[], int NMon) {

  // Set up the filenames, with the full path.  
  // Files are e.g. 'd1_2.xyz' and 'd1_2.key' because tinker hates
  // decimal points in filenames
  string filename = Params::Parameters().GetMMPath() + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d_%d",indexA,indexB);
  filename += label;
  string xyzfile = filename + ".xyz";
  string keyfile = filename + ".key";

  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Dimer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }

    if ( Params::Parameters().UseEmbeddingCharges() ) {
    // Count total atoms: monomer atoms + charges
    int Ntot = 0;
    for (int i=1;i<=NMon;i++) {
      Ntot += Monomers[i].GetNumberOfAtoms();
    }
    fprintf(xyz,"%d  Dimer (%d,%d) with embedding charges\n",
	    Ntot,indexA,indexB);

    // Print the geometry
    int shift = 0;
    for (int i=1;i<=NMon;i++) {
      if (i==indexA || i==indexB)
	Monomers[i].PrintTinkerCartesian(shift,false,xyz);
      else
	Monomers[i].PrintTinkerEmbeddingCharges(shift,xyz);

      // Shift for the next monomer
      shift += Monomers[i].GetNumberOfAtoms();
    }
  }
  else
    PrintTinkerCartesian(xyz);


  fclose(xyz);

  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
  if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
    printf("Dimer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }
  fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
  fclose(key);
  
}




// Returns command string for running QChem job
string Dimer::RunQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string infile = filename + ".in";
  string outfile = filename + ".out";

  // Execute Q-Chem
  // Q-Chem wants to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_prefix = "d";
  local_prefix += label;
  string local_outfile = local_prefix + ".out";
  string local_infile = local_prefix +  ".in";
 
  cmd += "qchem " + local_infile;
  cmd += " ";
  cmd += local_outfile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// Wrapper that controls which type of MM job is executed
// Returns command string for running the job
string Dimer::RunMMJob() {
  string job;

  if (Params::Parameters().GetMMType()==1) // Tinker
    job = RunTinkerJob();
  else if (Params::Parameters().GetMMType()==2) // AIFF
    job = "echo no dimer job for CamCasp";
  else if (Params::Parameters().GetMMType()==3) // QChem
    job = RunQChemJob(true);
  else {
    printf("ERROR: Dimer::RunMMJob(): Unknown MM type.\n");
    exit(1);
  }
  return job;
}

// Returns command string for running the job
string Dimer::RunTinkerJob() {
  // Set up the filenames, with no path info
  // Files are e.g. 'm1.xyz' and 'm1.key'
  string filename;
  char label[10]; // index ought to have fewer than 10 digits
  sprintf(label,"d%d_%d",indexA,indexB);
  filename = label;
  string infile = filename + ".xyz";
  string outfile = filename + ".out";

  // Execute Tinker
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "analyze " + infile;
  cmd += " e > ";
  cmd += outfile + "; ";

  // Third command, Rename output from d1_2.out -> d1.2.out
  string newfile = outfile;
  int index = newfile.find("_");
  newfile.replace(index,1,".");
  //printf("new file name = %s\n",newfile.c_str());
  cmd += "mv " + outfile;
  cmd += " " + newfile;
  cmd += "; ";

  /* Actual job running & checking of status have been moved to main.C 
     to simplify parallel implementation.

  printf("Running Tinker energy calculation on %s\n",infile.c_str());
  fflush(stdout); // flush the output stream
  //printf("Executing: %s\n",cmd.c_str());
  system(cmd.c_str());
  */

  // If doing force job, compute the gradient
  if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    // Run the job
    cmd += "minimize " + infile;
    cmd += " 100000000 > ";
    
    // Rename the force file
    string force_file = filename + ".force";
    index = force_file.find("_");
    force_file.replace(index,1,"."); // replace "_" with "." in name
    cmd += force_file + ";";

    // Remove tmp file and extra geom file created by minimize job
    cmd += "rm -f " + infile;
    cmd += "_2;";
  }

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  //  Actual job execution has been moved to main.C

  return cmd;
}

// Read in the Q-Chem energy
double Dimer::ReadQChemEnergy(bool MM_job) {
  double energy, cp_en[4] ={0};
  int cp_en_index = 0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();
  
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }

  double nuc_chg, chg_chg; // for if using embedding charges

  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);
    string match2 = line.substr(0,11);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      //printf("Std: MATCH FOUND!: %s\n",line.c_str());
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++) 
	iss >> tmp; // throw away text 
     iss >> cp_en[cp_en_index]; // read energy
     //printf("cp_en[%d] = %f\n",cp_en_index,cp_en[cp_en_index]);
    }
    if ( match2=="User input:" ) {
      cp_en_index++;
    }

    // Search for extra energies that arise when using embedding charges
    if ( match==" Nucleus-charge energy" ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> nuc_chg; // read energy 
    }

    // Search for extra energies that arise when using embedding charges
    if ( match==" Charge-charge energy " ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> chg_chg; // read energy 
    }

  }

  energy = cp_en[1];
  if ( (energy==0.0) || !(energy/energy == 1) ) {
    printf("Dimer::ReadQChemEnergy : Dimer Single Point Energy Error '%s'\n",
          out_filename.c_str());
    printf("energy = %f\n",energy);
    exit(1);
  }


  if ( Params::Parameters().UseEmbeddingCharges() ) {
    double old = energy;
    energy = energy - (nuc_chg + chg_chg);

    if ( Params::Parameters().PrintLevel() > 0 ) {
      printf("Embedding charges in use.  Correcting Final Energy\n");
      printf("Original:  %15.9f\n",old);
      printf("Nuc-chg =  %15.9f, Chg-chg = %15.9f\n",nuc_chg,chg_chg);
      printf("Corrected: %15.9f\n",energy);
    }
  }

  if ( Params::Parameters().PrintLevel() > 0) printf("QM Energy = %15.9f\n",energy);
 
  // Close the force file
  infile.close();
  return energy;
}





double Dimer::ReadQChemCounterpoiseInteractionEnergy() {

  double energy=0, cp_e[3]={0};
  int cp_energy_index=0; //goes from 0 to 2

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  //cout << " out_filename = " << out_filename << "\n";
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : Cannot open file '%s'\n",
           out_filename.c_str());
    exit(1);
  }

  //Remember that now our .out files has 3 energies and the corrected
  //dimer SPE is now E0-E1-E2 which we shall use to calculate
  //interaction energies.  Also, interaction energy will come out to
  //be equal to E0-E1-E3 (so no need to calculate monomer energies at
  //monomer basis.  cp_e0 = dimer energy, cp_e1 = mon1 energy at dimer
  //basis, cp_e2 = mon2 energy at dimer basis
  
  // Read in the data
  string line;           
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);
    string match2 = line.substr(0,3);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      //printf("CP: MATCH FOUND!!!  line=%s\n",line.c_str());
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++)
        iss >> tmp; // throw away text
      iss >> cp_e[cp_energy_index]; // read energy 
      //printf("cp_en[%d] = %f\n",cp_energy_index,cp_e[cp_energy_index]); 
      cp_energy_index++;
   }
  
  }
  
  if ( !(cp_energy_index == 3) || !(cp_e[0]/cp_e[0] == 1) || !(cp_e[1]/cp_e[1] == 1) || !(cp_e[2]/cp_e[2] == 1)  ) {
     printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : At Least One Single Point Energy Not Found or Incorrect '%s'\n",
           out_filename.c_str());
     exit(1);            
  }

  energy = cp_e[0] -cp_e[1]-cp_e[2];
  //cout << "energy = " << energy << "\n"; 
  if ( Params::Parameters().PrintLevel() > 0)
    printf("d(%d,%d) QM CP-corrected interaction Energy = %15.9f\n",
	   indexA,indexB,energy);
 
  // Convert kJ/mol -> hartrees
//  energy /= HartreesToKJpermole;
//  cout << "energy = " << energy << "\n";

  // Close the force file
  infile.close();

  return energy;

}



/* Following method of calculating CP corrected dimer energy is commented out
// Read in the Q-Chem counterpoise-corrected interaction energy
// Not sure how this works with embedding charges.
double Dimer::ReadQChemCounterpoiseInteractionEnergy() {

  double energy;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }


  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,12);

    // Search for final q-chem energy
    if ( match=="  DE, kJ/mol" ) {
            istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++) 
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }

  }

 
  if ( Params::Parameters().PrintLevel() > 0) 
    printf("QM CP-corrected interaction Energy = %15.9f\n",energy);
 
  // Convert kJ/mol -> hartrees
  energy /= HartreesToKJpermole;

  // Close the force file
  infile.close();

  return energy;
}
*/

void Dimer::ReadQMResults() {

  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().TinkerDebug() ) {
    Energy_QM = ReadTinkerEnergy();
  }
  else 
    Energy_QM = ReadQChemEnergy();

  if (Params::Parameters().DoCounterpoise()) {
    dEint_QM = ReadQChemCounterpoiseInteractionEnergy();
    //printf("Counterpoise-corrected interaction energy = %15.9f\n",dEint_QM);
  }
  else
    dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();

  if ( Params::Parameters().PrintLevel() > 0) {
    printf("dEint_QM = %15.9f  ",dEint_QM);
    printf("Energy_QM = %f, monA = %f, monB = %f\n",Energy_QM,
	   MonA.GetQMEnergy(),MonB.GetQMEnergy());
  }

  // Calculate MP2 energy correction, if needed
  double MP2_dispersion_correction;
  if ( Params::Parameters().DoMP2DispersionCorrection() ) {
    MP2_dispersion_correction = ComputeInteratomicMP2TwoBodyDispersionCorrection();

    printf("Standard MP2 interaction energy = %8.3f\n", 
	   dEint_QM*HartreesToKJpermole);
    printf("MP2 correction = %f kJ/mol\n",
	   MP2_dispersion_correction*HartreesToKJpermole);
    printf("Dispersion-corrected interaction energy = %8.3f\n",
	    (dEint_QM + MP2_dispersion_correction)*HartreesToKJpermole);

  }




  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    if ( Params::Parameters().TinkerDebug() ) 
      SetMMGradient();
    else 
      SetQMGradient();
  }
}

// Use "other" monomer as a reference for Mon B;
void Dimer::ReadQMResults(Monomer& otherB) {

  //printf("Using Monomer %d as a reference\n",otherB.GetIndex() ); fflush(stdout);

  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().TinkerDebug() ) {
    Energy_QM = ReadTinkerEnergy();
  }
  else 
    Energy_QM = ReadQChemEnergy();

  if (Params::Parameters().DoCounterpoise()) {
    dEint_QM = ReadQChemCounterpoiseInteractionEnergy();
    //printf("Counterpoise-corrected interaction energy = %15.9f\n",dEint_QM);
  }
  else {
    dEint_QM = Energy_QM - MonA.GetQMEnergy() - otherB.GetQMEnergy();
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("dEint_QM = %15.9f  ",dEint_QM);
    printf("Energy_QM = %f, monA = %f, monB = %f\n",Energy_QM,
	   MonA.GetQMEnergy(),otherB.GetQMEnergy());
  }

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    if ( Params::Parameters().TinkerDebug() ) {
      SetMMGradient();
    }
    else {
      SetQMGradient();
    }
    MonB.SetQMGradient(otherB);
  }
}



void Dimer::ReadMMResults() {

  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient(); // get energies & gradient simultaneously in this case
  }
  else {
    SetMMEnergy();
    if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
      SetMMGradient();
    }
  }

  dEint_MM = Energy_MM - MonA.GetMMEnergy() - MonB.GetMMEnergy();
   
}

// Version for use with PBC: Monomer B is a periodic image of
// reference monomer otherB
void Dimer::ReadMMResults(Monomer& otherB) {

  // Read in the Energies, compute 2-body interaction energies
  SetMMEnergy();
  dEint_MM = Energy_MM - MonA.GetMMEnergy() - otherB.GetMMEnergy();
   
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient();
    MonB.SetMMGradient(otherB);
  }
}

void Dimer::SetMMEnergy() {
  double energy;
  if ( Params::Parameters().GetMMType() == 1 ) // Tinker
    energy = ReadTinkerEnergy();
  else if ( Params::Parameters().GetMMType() == 2) { // AIFF
    energy = ComputeAIFFEnergy();
  }
  else if ( Params::Parameters().GetMMType() == 3) { // QChem
    energy = ReadQChemEnergy(true);
  }
  else {
    printf("Dimer::SetMMEnergy: Unknown MM_type: %d\n",Params::Parameters().GetMMType());
    exit(1);
  }

  Energy_MM = energy;
  //printf("MM Energy = %15.9f\n",Energy_MM);


}

double Dimer::ReadTinkerEnergy() {

  string path = Params::Parameters().GetMMPath();

  // Search for the "Polarization" string
  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"/d%d.%d",MonA.GetIndex(), MonB.GetIndex() );
  string filename = path + label;
  string out_filename = filename + ".out";
  

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadTinkerEnergy : Cannot open file '%s'\n",out_filename.c_str());
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

// Get the QM Gradient, wrapper routine
void Dimer::SetQMGradient() {
  string path = Params::Parameters().GetQMPath();
  Grad_QM = ReadGradient(path,1);
  QM_Grad_Init = 1;
}

// Get the MM Gradient, wrapper routine
void Dimer::SetMMGradient() {

  if (Params::Parameters().GetMMType() == 1) { // Tinker 
    string path = Params::Parameters().GetMMPath();
    Grad_MM = ReadGradient(path,2);
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    Grad_MM = ComputeMultipoleGradient();
  }
  else if (Params::Parameters().GetMMType() == 3) { //QChem
    string path = Params::Parameters().GetMMPath();
    Grad_MM = ReadGradient(path,1);
  }
  else {
    printf("Dimer::SetMMGradient() - MM type = %d is unknown\n",
	   Params::Parameters().GetMMType());
    exit(1);
  }

  MM_Grad_Init = 1;
}

// Main routine for reading in gradient files
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
Vector Dimer::ReadGradient(string path, int type) {

  Vector grad(3*Natoms), gradA( 3*Natoms ),  gradB( 3*Natoms );
  int cp_grad_index =0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.force'
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d",indexA,indexB);
  string filename = path + label;
  if (type == 2) // Tinker MM job
    filename += ".force";
  else // other
    filename += ".out";


  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadGradient : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }
  if ( !(Params::Parameters().DoCounterpoise()) || type==2 ) {
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
 }

 if ( Params::Parameters().DoCounterpoise() && type==1 ) {
   // Read in the data - search for the "Nuclear forces:" string
   string line;               
   while ( !infile.eof() ) {
     getline(infile,line);
     // Search for final q-chem energy
     if ( line==" Nuclear forces:" && (cp_grad_index == 0) ) {
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
       //      break;
       cp_grad_index++;
     }
     
     if ( line==" Nuclear forces:" && (cp_grad_index == 1) ) {
       getline(infile,line); // throw away header line
       for (int i=0;i<(MonA.GetNumberOfAtoms()+MonB.GetNumberOfAtoms());i++) {
	 getline(infile,line);
	 istringstream iss(line);
	 string tmp;
	 iss >> tmp; // throw away the atom index
	 for (int j=0;j<3;j++) {
	   iss >> gradA[3*i+j]; // Store the gradient elements
	 }
       }
       //      break;
       cp_grad_index++;
     }
     
     if ( line==" Nuclear forces:" && (cp_grad_index == 2) ) {
       getline(infile,line); // throw away header line
       for (int i=0;i<(MonA.GetNumberOfAtoms()+MonB.GetNumberOfAtoms());i++) {
	 getline(infile,line);
	 istringstream iss(line);
	 string tmp;
	 iss >> tmp; // throw away the atom index
	 for (int j=0;j<3;j++) {
	   iss >> gradB[3*i+j]; // Store the gradient elements
	 }
       }
       //      break;
       cp_grad_index++;
     }
     
     
   }
   
   if ( !(type == 2) && !(cp_grad_index == 3) ) {
     printf("Dimer::ReadQMGradient : At Least One Counterpoise Gradient Not Found or Incorrect '%s'\n",
	    filename.c_str());
     exit(1);
   }
   
   // Find the CP gradient.  We also add in the contributions from isolated
   // MonA and MonB, because these get subtracted off later.  
   for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
     grad[j] = grad[j] - gradA[j] - gradB[j] + MonA.GetQMGradient()[j];
   }
   for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
     int k = j + 3*MonA.GetNumberOfAtoms(); // in dimer gradient, need offset to second monomer
     grad[k] = grad[k] - gradA[k] - gradB[k] + MonB.GetQMGradient()[j];
   }
 } // end Counterpoise Gradient section
 
 infile.close();  // Close the file
 
 //PrintGradient("Dimer::ReadGradient(): Dimer Gradient",grad);
 return grad; 
}

/*Vector Dimer::ReadQMCPGradient(string path) {

  Vector grad(3*Natoms), gradA( MonA.GetNumberOfAtoms() ),  gradB( MonB.GetNumberOfAtoms() );
  int cp_grad_index =0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.force'
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d",indexA,indexB);         
  string filename = path + label;       
  filename += ".out";  


  // Open the force file   
  ifstream infile;         
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadGradient : Cannot open file '%s'\n",  
           filename.c_str());
    exit(1);
  }

  // Read in the data - search for the "Nuclear forces:" string
  string line;               
  while ( !infile.eof() ) {
    getline(infile,line);
    // Search for final q-chem energy
    if ( line==" Nuclear forces:" && (cp_grad_index == 0) ) {
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
//      break;
      cp_grad_index++;
    }

    if ( line==" Nuclear forces:" && (cp_grad_index == 1) ) {
      getline(infile,line); // throw away header line

      for (int i=0;i<MonA.GetNumberOfAtoms();i++) {                
        getline(infile,line);           
        istringstream iss(line);
        string tmp;
        iss >> tmp; // throw away the atom index
        for (int j=0;j<3;j++) {
          iss >> gradA[3*i+j]; // Store the gradient elements
        }
      }
//      break;
      cp_grad_index++;
    }

    if ( line==" Nuclear forces:" && (cp_grad_index == 2) ) {
      getline(infile,line); // throw away header line

      for (int i=0;i<MonB.GetNumberOfAtoms();i++) {                
        getline(infile,line);           
        istringstream iss(line);
        string tmp;
        iss >> tmp; // throw away the atom index
        for (int j=0;j<3;j++) {
          iss >> gradB[3*i+j]; // Store the gradient elements
        }
      }
//      break;
      cp_grad_index++;
    }


  }

  if ( !(cp_grad_index == 3) ) {
   printf("Dimer::ReadQMCPGradient : At Least One Gradient Not Found or Incorrect '%s'\n",
           filename.c_str());
     exit(1);
  }

  for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
    grad[j] = grad[j] - gradA[j];
  }

  for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
    int k = j + 3*MonA.GetNumberOfAtoms(); // in dimer gradient, need offset to second monomer
    grad[k] = grad[k] - gradB[j];
  }

  infile.close();  // Close the file

  //PrintGradient("Dimer Gradient",grad);
  return grad;
}

*/

// Print out the Cartesian coordinates in xyz format
void Dimer::PrintDimerCartesian() {
  MonA.PrintMonomerCartesian();
  MonB.PrintMonomerCartesian();
}

// Print out $molecule section for Q-Chem
void Dimer::PrintQChemCartesian(FILE *outfile) {
  fprintf(outfile,"$molecule\n%d %d\n",charge,spin);

/*  if ( Params::Parameters().DoCounterpoise() ) {
    fprintf(outfile,"--\n%d %d\n",MonA.GetChargeState(), MonA.GetSpinState());
  }
*/
   MonA.PrintMonomerCartesian(outfile);
/*
  if ( Params::Parameters().DoCounterpoise() ) {
    fprintf(outfile,"--\n%d %d\n",MonB.GetChargeState(), MonB.GetSpinState());
  }
*/
   MonB.PrintMonomerCartesian(outfile);
   fprintf(outfile,"$end\n\n");

 //  if ( Params::Parameters().DoCounterpoise() ) { // make special qchem input file to grab bsse correction manually (added on 9th april)
      


}

// Print out the Cartesian coordinates in Tinker xyz format
void Dimer::PrintTinkerCartesian(FILE *outfile) {
  // Print tile line
  fprintf(outfile,"%d  Dimer (%d,%d)\n",Natoms,indexA,indexB); 

  // Print first and second monomers
  MonA.PrintTinkerCartesian(0,false,outfile);
  // apply shift to second monomer to get proper indexing
  int shift = MonA.GetNumberOfAtoms();
  MonB.PrintTinkerCartesian(shift,false,outfile);

}



// Wrapper for printing QM gradient
void Dimer::PrintQMGradient(string title) {
  printf("QM_Grad_Init = %d\n",QM_Grad_Init);  
  if (QM_Grad_Init) 
    PrintGradient(title,Grad_QM);
  else 
    printf("QM Gradient not initialized\n");
}

// Wrapper for printing <M gradient
void Dimer::PrintMMGradient(string title) {
  printf("MM_Grad_Init = %d\n",MM_Grad_Init);
  if (MM_Grad_Init) 
    PrintGradient(title,Grad_MM);
  else 
    printf("MM Gradient not initialized\n");
}

// Print out Gradient as an Nx3 matrix
void Dimer::PrintGradient(string title, Vector grad) {
  //printf("Dimer::PrintGradient incomplete\n");

  printf("%s\n",title.c_str());
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  for (int i=0;i<Na;i++)
    printf("%2s %15.10f %15.10f %15.10f\n",MonA.GetSymbol(i).c_str(),
	   grad[3*i],grad[3*i+1],grad[3*i+2]);
  for (int i=Na;i<Na+Nb;i++)
    printf("%2s %15.10f %15.10f %15.10f\n",MonB.GetSymbol(i-Na).c_str(),
	   grad[3*i],grad[3*i+1],grad[3*i+2]);

  
}

void Dimer::PrintAll() {

  printf("-----------------------------------------------\n");
  printf("Dimer (%d,%d)      %d atoms\n\n",indexA,indexB,Natoms);

  PrintQChemCartesian();
  printf("Intermolecular separation = %.3f Angstroms\n\n",Separation);
  
  printf("Energies:\nQM: %15.9f\nMM: %15.9f\n",Energy_QM,Energy_MM);
  printf("Interaction Energies:\nQM: %15.9f\nMM: %15.9f\n\n",dEint_QM,dEint_MM);

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    printf("Forces:\n");
    PrintQMGradient("QM:\n");
    PrintMMGradient("MM:\n");
  }



  printf("-----------------------------------------------\n");
}

// Get the damping factor for smooth pairwise local truncations
// Note: This is NOT the same damping factor used to damp short range
// induction energies.  That one is stored in the Param class and
// computed in the Atom class.
double Dimer::GetDampingFactor(double c0, double c1) {
  /* 
     Applies a smooth bump function between distances c1 and c0
     to smoothly damp from 1 -> 0, according to:

     Subotnik et al, J. Chem. Phys. 128, 034103 (2008)

     c0 = distance after which damp = 0  (complete damping)
     c1 = distance before which damp = 1 (no damping)
     Between c1 and c0, we have smooth decay from 1 -> 0.
   */
  double damp;
  double R =  Separation; //distance between the two fragments

  if (R <= c1) {
    damp = 1.0; 
  }
  else if (R >= c0) {
    damp = 0.0; 
  }
  else {
    double tmp = 2*fabs(c1 - c0)/(c1 - R) - fabs(c1 - c0)/(R - c0);
    damp = 1.0/(1.0 + exp(tmp));
  }

  if ( Params::Parameters().PrintLevel() > 0 )
    printf("GetDampingFactor: d(%d,%d)  R = %.3f, damp = %.4f\n",
	   indexA,indexB,R,damp);
	
  return damp;

}

int Dimer::GetNzero(double c0, double c1) {
  /*
   this function accounts for an increment in Nzero whenever the separation between monomers is greater than c0

   */
  int Nzero;
  double R =  Separation; //distance between the two fragments

  if ( (R >= c0) && (R != c1) ) {
    Nzero = 1;     
  }
  else {
    Nzero = 0;         
  }

  return Nzero;

}

int Dimer::GetNfull(double c0, double c1) {        
  /*
   this function accounts for an increment in Nfull whenever the separation between monomers is less than c1

   */
  int Nfull;
  double R =  Separation; //distance between the two fragments

  if (R <= c1) {
    Nfull = 1;
  }
  else {            
    Nfull = 0;        
  }                

  return Nfull;

}

int Dimer::GetNdamp(double c0, double c1) {        
  /*
   this function accounts for an increment in Ndamp whenever the separation between monomers is between c1 and c0

   */
  int Ndamp;
  double R =  Separation; //distance between the two fragments

  if ( (R < c0) && (R > c1) ) {
    Ndamp = 1;
  }
  else {            
    Ndamp = 0;        
  }                

  return Ndamp;

}


// Compute distances between atoms in the 2 fragments
void Dimer::ComputeIntermolecularDistances() {

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  printf(" Mon 1    Mon 2       Distance (Ang)\n");
  for (int i=0;i<Na;i++) {
    string atom1 = MonA.GetAtom(i).GetSymbol();
    for (int j=0;j<Nb;j++) {		
      string atom2 = MonB.GetAtom(j).GetSymbol();
      
      double dist = MonA.GetAtom(i).GetInterAtomicDistance( MonB.GetAtom(j) );
      printf("%2s(%3d)  %2s(%3d)       %7.4f\n",atom1.c_str(),i,
	     atom2.c_str(),j+Na,dist);
    }
  }
  printf("\n");
}

// Evaluate the force field energy using ab initio-derived parameters
double Dimer::ComputeAIFFEnergy() {
  double Eaiff = 0.0, E_es = 0.0, E_pol = 0.0, E_2b_disp = 0.0;

  printf("\nEvaluating ab initio force field contribution for dimer (%d,%d)\n",
	   indexA,indexB);

  bool do_es = Params::Parameters().DoAIFFElectrostatics();
  bool do_ind = Params::Parameters().DoAIFFInduction();
  bool do_2b_disp = Params::Parameters().DoAIFF2BodyDispersion();

  
  // Compute the electrostatic + induction/polarization energy -- old, slower routine
  /*
  double E_es_pol = 0.0;
  if (do_es || do_ind) 
    E_es_pol = ComputeMultipoleInteractions(); 

  // For printing, grab these energies
  E_es = GetMMElectrostaticEnergy(); 
  E_pol = GetMMInductionEnergy(); 
  */

  // Compute the permanent electrostatics
  if (do_es)
    E_es = ComputeAIFFElectrostatics();

  // Compute the 2-body induction energy
  if (do_ind)
    E_pol = ComputeAIFFInduction();


  // Compute the 2-body intermolecular dispersion energy
  if (do_2b_disp)
    E_2b_disp = ComputeTwoBodyDispersion();


  printf("  Electrostatic energy for dimer (%d,%d) = %8.4f kJ/mol\n",indexA,indexB,E_es*HartreesToKJpermole);
  printf("  Induction energy for dimer (%d,%d) =     %8.4f kJ/mol\n",indexA,indexB,E_pol*HartreesToKJpermole);
  if (!Params::Parameters().DoDispersionDamping() ) 
    printf("  Dispersion energy for dimer (%d,%d) =    %8.4f kJ/mol (undamped)\n",indexA,indexB,E_2b_disp*HartreesToKJpermole);
  else
    printf("  Dispersion energy for dimer (%d,%d) =    %8.4f kJ/mol\n",indexA,indexB,E_2b_disp*HartreesToKJpermole);
  
  // Add up the individual contributions
  Eaiff = E_es + E_pol + E_2b_disp;

  printf("  Total AIFF energy for dimer (%d,%d) =    %8.4f kJ/mol\n",indexA,indexB,Eaiff*HartreesToKJpermole);

  return Eaiff;
}


// Build Tab matrix for interaction between 2 monomer multipoles and
// compute their electrostatic & induction energies
double Dimer::ComputeMultipoleInteractions() {
  double energy = 0.0;

  if (Params::Parameters().PrintLevel() > 0) {
    printf("\nComputing permanent multipole interactions for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  // Variables to store total permanent electrostatic, induction & total energies
  double Ees = 0.0, Eind = 0.0, Etot = 0.0;


  // Step 1: Build the geometric interaction matrix Tab.  This matrix
  // contains only geometry information and gets used for evaluating
  // both the permanant and induced multipole electrostatics.  We
  // build two versions.  The first, Tabs, stores the standard
  // matrix. The second, DampedTabs, stores a version that is damped.
  // Damping is necessary for computing induction energies.

  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor for this dimer

  
  /*
  // hack: different damping factors for water-water and glycine-water.
  // Assumes glycine is monomer 1.  Let that be set by damping factor.  Otherwise,
  // we assume beta_damp = 1.45 because it's water-water.
  if (indexA!=1 && indexB!=1) {
    beta_damp = 1.45;
    printf("*** HACK: Setting beta_damp = %.3f for dimer (%d,%d)\n",beta_damp,indexA,indexB);
  }
  */


  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  Tabs = new Matrix[NatomsA*NatomsB];
  DampedTabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the undamped version of Tab.  "-999.0" is a dummy
      // parameter to turn damping off.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      Tabs[iB*NatomsA+iA].Initialize(tmpTab);  

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, beta_damp);
      DampedTabs[iB*NatomsA+iA].Initialize(tmpTab);  
    }
  }


  // Step 2: Evaluate permanent and induced electrostatic interactions

  // Step 2a: Permanent multipole contributions
  //          E_es = sum(tu) QA(t)*Tab(t,u)*QB(u) 
  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
      
      // Evaluate the energy
      Vector tmp(QA.GetLength());
      // Tab*QB
      tmp = Tabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
      // QA*(Tab*QB)
      double Energy = QA.GetMoments().DotProduct(tmp);
      Ees += Energy;

      /* Debug printing: break energy down into separate contributions */
      if (Params::Parameters().PrintLevel() > 0) {
      // Create text label for each atom, and print the energy contribution
      char labelA[10], labelB[10];
      sprintf(labelA,"%s%d",MonA.GetAtom(iA).GetSymbol().c_str(),
	      MonA.GetAtom(iA).GetAtomIndex());
      sprintf(labelB,"%s%d",MonB.GetAtom(iB).GetSymbol().c_str(),
	      MonB.GetAtom(iB).GetAtomIndex());
      printf("Electrostatic contribution for %s(A)--%s(B) interaction = %f\n",
	     labelA, labelB,Energy*HartreesToKJpermole);
      }
      /* Debug printing: break energy down into separate contributions */
      if (Params::Parameters().PrintLevel() > 1) {
        Matrix Edecomp(5,5);
	Matrix Tab(Tabs[iB*NatomsA+iA]);
	
	for (int t=0;t<QA.GetLength();t++)
	  for (int u=0;u<QB.GetLength();u++) {
	    // chg-chg
	    if (t==0 && u==0) 
	      Edecomp(0,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-dip
	    else if (t==0 && u>=1 && u<=3) {
	      Edecomp(0,1) += QA(t)*Tab(t,u)*QB(u);
	    }
	    else if (u==0 && t>=1 && t<=3)
	      Edecomp(1,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-quad
	    else if (t==0 && u>=4 && u<=8)
	      Edecomp(0,2) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=4 && t<=8)
	      Edecomp(2,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-oct
	    else if (t==0 && u>=9 && u<=15)
	      Edecomp(0,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=9 && t<=15)
	      Edecomp(3,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-hexadec
	    else if (t==0 && u>=16 && u<=24)
	      Edecomp(0,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=16 && t<=24)
	      Edecomp(4,0) += QA(t)*Tab(t,u)*QB(u);
	    
	    // dip-dip
	    else if (t>=1 && t<=3 && u>=1 && u<=3)
	      Edecomp(1,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-quad
	    else if (t>=1 && t<=3 && u>=4 && u<=8)
	      Edecomp(1,2) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=4 && t<=8)
	      Edecomp(2,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-oct
	    else if (t>=1 && t<=3 && u>=9 && u<=15)
	      Edecomp(1,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=9 && t<=15)
	      Edecomp(3,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-hexadec
	    else if (t>=1 && t<=3 && u>=16 && u<=24)
	      Edecomp(1,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=16 && t<=24)
	      Edecomp(4,1) += QA(t)*Tab(t,u)*QB(u);
	    
	    // quad-quad
	    else if (t>=4 && t<=8 && u>=4 && u<=8)
	      Edecomp(2,2) += QA(t)*Tab(t,u)*QB(u);
	    // quad-oct
	    else if (t>=4 && t<=8 && u>=9 && u<=15)
	      Edecomp(2,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=4 && u<=8 && t>=9 && t<=15)
	      Edecomp(3,2) += QA(t)*Tab(t,u)*QB(u);
	    // quad-hexadec
	    else if (t>=4 && t<=8 && u>=16 && u<=24)
	      Edecomp(2,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=4 && u<=8 && t>=16 && t<=24)
	      Edecomp(4,2) += QA(t)*Tab(t,u)*QB(u);
	    
	    // oct-oct
	    else if (t>=9 && t<=15 && u>=9 && u<=15)
	      Edecomp(3,3) += QA(t)*Tab(t,u)*QB(u);
	    // oct-hexa
	    else if (t>=9 && t<=15 && u>=16 && u<=24)
	      Edecomp(3,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=9 && u<=15 && t>=16 && t<=24)
	      Edecomp(4,3) += QA(t)*Tab(t,u)*QB(u);
	    
	    // hexadec-hexadec
	    else if (t>=16 && t<=24 && u>=16 && u<=24)
	      Edecomp(4,4) += QA(t)*Tab(t,u)*QB(u);
	  }
	Edecomp.Scale(HartreesToKJpermole);
	Edecomp.Print("Electrostatic energy decomposition"); printf("\n");
      }
    }
  }   

  // If not doing induction, exit here
  bool do_ind = Params::Parameters().DoAIFFInduction();
  if (! do_ind) {
    // Store the electrostaic energy contibution, in hartrees.
    E_Electrostatic_MM = Ees;   
    return Ees;
  }

    
  // Step 2b: Induced 2-body multipole contributions
  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction energy for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);

  // Initialize a few variables
  bool iterate = true;
  int cycle = 0;
  double ind_conv = Params::Parameters().GetInductionConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  if (Params::Parameters().PrintLevel() > 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);

  double Eind_old = 0.0;

  // Create lists that will store Induced Moments for each atom on A
  // and B and initialize them to the proper sizes based on the higher
  // rank object: polarizabilities or multipoles.  This is important:
  // CamCasp sometimes prints polarizabilities with high ranks & lots
  // of zeros, even if a lower rank would do.  We need two copies of
  // each list.
  Multipole *dQA = new Multipole[NatomsA];
  Multipole *old_dQA = new Multipole[NatomsA];
  for (int iA=0;iA<NatomsA;iA++) {
    int Rmom = MonA.GetAtom(iA).GetMultipoleMoments().GetRank();
    int Rpol = MonA.GetAtom(iA).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    printf("MonA atom %d: Rmom = %d, Rpol = %d\n",iA,Rmom,Rpol);

    dQA[iA].Initialize(max_rank);
    old_dQA[iA].Initialize(max_rank);
  }
  Multipole *dQB = new Multipole[NatomsB];
  Multipole *old_dQB = new Multipole[NatomsB];
  for (int iB=0;iB<NatomsB;iB++) {
    int Rmom = MonB.GetAtom(iB).GetMultipoleMoments().GetRank();
    int Rpol = MonB.GetAtom(iB).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    printf("MonB atom %d: Rmom = %d, Rpol = %d\n",iB,Rmom,Rpol);

    dQB[iB].Initialize(max_rank);
    old_dQB[iB].Initialize(max_rank);
  }

  Eind = 1000000.0;  // start with huge nonsense energy

   // Start wall clock timer
    time_t start_time, stop_time;
    start_time = time(NULL);  

  if (Params::Parameters().PrintLevel() > 0) {
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
    double EindA = 0.0,  EindB = 0.0;
    for (int iA=0;iA<NatomsA;iA++) {
      old_dQA[iA] = dQA[iA];
      dQA[iA].Set();
    }
    for (int iB=0;iB<NatomsB;iB++) {
      old_dQB[iB] = dQB[iB];
      dQB[iB].Set();
    }


    // Induce multipoles
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      // Grab some data we need
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Polarizability PolA(MonA.GetAtom(iA).GetPolarizability(),true);
      int NpolA = PolA.GetLength();
      int NmomA = QA.GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	// Grab some data we need
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	Polarizability PolB(MonB.GetAtom(iB).GetPolarizability(),true);
	int NpolB = PolB.GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	int dimT1 = Tab.GetRows();
	int dimT2 = Tab.GetCols();

	// Induce multipoles on monomer A due to monomer B
	// dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))
	for (int a=0;a<NpolA;a++) {// loop over elements of dQA
	  for (int t=0;t<min(NpolA,dimT1);t++) {
	    for (int u=0;u<min(dimT2,NmomB);u++) {
	      dQA[iA](a) -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQB[iB](u));
	    }
	  }  
	}

	// Induce multipoles on monomer B due to monomer A
	// dQB(a) = dQB(a) - polB(a,t)*Tab(u,t)*(QA(u) + old_dQA(u))
	for (int a=0;a<NpolB;a++) {// loop over elements of dQA
	  for (int t=0;t<min(NpolB,dimT2);t++) {
	    for (int u=0;u<min(dimT1,NmomA);u++) {
	      dQB[iB](a) -= PolB(a,t)*Tab(u,t)*(QA(u) + old_dQA[iA](u));
	    }
	  }
	}
      }
    }
    
    // Compute the induction energy using the induced multipoles determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);

	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);
	  }
      }
    }
    Eind = EindA + EindB;
  
    // Print out results for this cycle
    if (cycle==1 && Params::Parameters().PrintLevel() > 0) 
      printf(" %3d       %12.6f     *********** kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole);
    else if (Params::Parameters().PrintLevel() > 0) 
      printf(" %3d       %12.6f     %11.6f kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole);

    // Check convergence based on the energy change (in kJ/mol).
    if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies converged after %d iterations\n",cycle);
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
      iterate = false;
    }    

    if ( cycle == Params::Parameters().GetMaxPolarizationCycles()  && iterate == true ) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies failed to converge after %d iterations\n",cycle);
      Params::Parameters().Warning();
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
    }
  }


    // Stop the timer and print out the time
    stop_time = time(NULL);
    double elapsed_time = difftime(stop_time,start_time);
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Iterative Induction energy wall time = %.5f seconds\n",elapsed_time);

  // Hack: looking at what happens if we don't use damping for final induction energy
  if (Params::Parameters().GetMaxOptCycles()==1) {
    printf("Recomputing dimer induction energy without damping\n");
    double EindA = 0.0;
    double EindB = 0.0;
    Eind = 0.0;
    // Compute the induction energy using the induced multipoles determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A
      
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();
      
      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      
	
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	
	
	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	
	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);
	  }
      }
    }
    Eind = EindA + EindB;
  }

  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles ***\n");
    for (int iA=0;iA<NatomsA;iA++) {
      char label[50];
      sprintf(label,"Monomer 1, %s%d",MonA.GetAtom(iA).GetSymbol().c_str(),
	      MonA.GetAtom(iA).GetAtomIndex());
      string str = label;
      
      
      dQA[iA].Print(str);
    }
    for (int iB=0;iB<NatomsB;iB++) {
      char label[50];
      sprintf(label,"Monomer 2, %s%d",MonB.GetAtom(iB).GetSymbol().c_str(),
	      MonB.GetAtom(iB).GetAtomIndex());
      string str = label;
      dQB[iB].Print(str);
    }
  }

  // Store the energy contibutions, in hartrees.
  E_Electrostatic_MM = Ees;
  E_Induction_MM = Eind;
  Etot = Ees + Eind;

  if (Params::Parameters().PrintLevel() > 0) {
    // Print out a summary of the classical energy just calculated
    printf("--------------------------------------------------------\n");
    printf(" Total Classical Interaction Energy for Dimer (%d,%d)\n",
	   indexA,indexB);
    printf("--------------------------------------------------------\n");
    printf("      Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
    printf("      Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
    printf("    ---------------------------------------\n");
    printf("      Total         = %12.6f kJ/mol\n",Etot*HartreesToKJpermole);
    printf("--------------------------------------------------------\n");
  }

  delete [] dQA;
  delete [] old_dQA;
  delete [] dQB;
  delete [] old_dQB;

  return Etot;
}

// Matrix algorithm for computing the AIFF electrostatics
double Dimer::ComputeAIFFElectrostatics() {

  if (Params::Parameters().PrintLevel() > 0) {
    printf("\nComputing permanent multipole interactions for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  double Ees = 0.0;

  // Step 1: Build the geometric interaction matrix Tab.  This matrix
  // contains only geometry information and gets used for evaluating
  // both the permanant multipole electrostatics.  No damping is
  // applied in Tab.
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  Tabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the undamped version of Tab.  "-999.0" is a dummy
      // parameter to turn damping off.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      Tabs[iB*NatomsA+iA].Initialize(tmpTab);  
    }
  }


  // Step 2: Evaluate permanent multipole electrostatic interactions
  //          E_es = sum(tu) QA(t)*Tab(t,u)*QB(u) 
  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
      
      // Evaluate the energy
      Vector tmp(QA.GetLength());
      // Tab*QB
      tmp = Tabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
      // QA*(Tab*QB)
      double Energy = QA.GetMoments().DotProduct(tmp);
      Ees += Energy;

    }
  }

  E_Electrostatic_MM = Ees;   
  return Ees;
}

// Matrix algorithm for computing the AIFF induction energy.  
double Dimer::ComputeAIFFInduction() {
  double Eind = 0.0;

  // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);

  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction energy for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor

  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);

  // Step 1: Build the damped geometric interaction matrix Tab for
  // each intermolecular pair of atoms
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  // They get stored as a vector of matrices.
  DampedTabs = new Matrix[NatomsA*NatomsB]; 
  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);
      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, beta_damp);

      DampedTabs[iB*NatomsA+iA].Initialize(tmpTab);  
    }
  }

  // Step 2: Evaluate induced multipole contributions for this dimer
  
  // Matrix algorithm: Need to solve: Z*dQ = -V0 
  // where V0 is a long vector containing the potential due to the
  // permanent multipole moments from the atoms on monA and monB, and
  // dQ is the vector of induced multipole moments we want to solve.

  // Z is a matrix, which has block structure that can be thought of
  // as:
  //      ---------               --
  //     | AA | AB |             | (pol^A)^{-1} if A==B
  // Z =  ---------    where Z = |
  //     | BA | BB |             | T^{AB}      if A!=B
  //      ---------               --
  //
  // In other words, the diagonal blocks AA and BB are themselves
  // block diagonal matrices consisting of the inverse of the
  // polarizability tensor on each atom.  The off diagonal elements
  // just look like blocks of the Tab matrix and its transpose.
  //
  // For efficiency, size the arrays appropriately based on the types
  // (rank) of multipoles we allow to be induced on each atom type.
  // In practice, that means dipoles on hydrogens and dipoles +
  // quadrupoles on non-hydrogen atoms.



  // Compute array sizes by counting up atom types on each monomer.
  // Induce dipoles (3) on H and dipoles & quadrupoles (3+5=8) 
  // on heavy atoms.  Size: 8*(heavy atoms) + 3*(H atoms)
  int sizes[NatomsA + NatomsB];
  int offset[NatomsA + NatomsB];

  offset[0] = 0;
  int dim=0;
  for (int iA=0;iA<NatomsA;iA++) {
    if (MonA.GetAtom(iA).GetAtomicNumber()==1) {
      sizes[iA] = 3;
      dim += 3;
    }
    else {
      sizes[iA] = 8; // other atoms - dipole & quadrupole blocks
      dim += 8;
    }
    if (iA > 0) 
      offset[iA] = offset[iA-1] + sizes[iA-1];
  }
  int dimA = dim; // save size of MonA atom blocks

  for (int iB=0;iB<NatomsB;iB++) {
    if (MonB.GetAtom(iB).GetAtomicNumber()==1) {
      sizes[iB+NatomsA] = 3;
      dim += 3;
    }
    else {
      sizes[iB+NatomsA] = 8; // other atoms - dipole & quadrupole blocks
      dim += 8;
    }
    offset[iB+NatomsA] = offset[iB+NatomsA-1] + sizes[iB+NatomsA-1];
  }

  // Initialize storage for Z, dQ, and V0
  Vector V0(dim);
  Vector dQ(dim);
  Matrix Z(dim,dim);
  

  // Build Z matrix.

  // Diagonal blocks first: inverse of polarizability.  Loop over
  // atoms on each monomer, invert polarizability, and put it in Z.
  for (int iA=0;iA<NatomsA;iA++) {

    // The code generally stores pols as 9x9 matrix, but some of those
    // rows/cols are zeroes, which make the matrix singular and
    // prevent inversion and would unnecessarily increase the
    // dimensionality of the matrices here.  Therefore, we grab just
    // the non-zero block.  This will be either a 3x3 dipole-dipole
    // block (H atoms) or the 8x8 block consisting of dipole and
    // quadrupole polarizabilities.
    int Npols = sizes[iA]; // dimensionality of the non-zero polarizabilities block.
 
    // Grab the right block of the Polarizability matrix.  always skip first row/col, since that
    // would be charge-charge, which is zero.  So use (1,1) as offset here.
    Matrix PolMat = MonA.GetAtom(iA).GetPolarizability().GetPolarizabilities().GetBlock(Npols,Npols,1,1);
    PolMat.Inverse();
    // Set the block.  
    Z.SetBlock(PolMat,offset[iA],offset[iA]);
  }
  for (int iB=0;iB<NatomsB;iB++) {

    int Npols = sizes[iB + NatomsA];

    // Grab the right block of the Polarizability matrix.  always skip first row/col, since that
    // would be charge-charge, which is zero.  So use (1,1) as offset here.
    Matrix PolMat = MonB.GetAtom(iB).GetPolarizability().GetPolarizabilities().GetBlock(Npols,Npols,1,1);
    PolMat.Inverse();
    Z.SetBlock(PolMat,offset[iB+NatomsA],offset[iB+NatomsA]);
  }

  // Now set off-diagonal blocks between pairs of atoms in different monomers
  //int offsetA, offsetB; 
  int offsetA_all = 0;
  for (int iA=0;iA<NatomsA;iA++) {
    // Find sizes, offset;
    int NpolsA = sizes[iA];
    int offsetA = offset[iA];

    for (int iB=0;iB<NatomsB;iB++) {
      // Find sizes, offset
      int NpolsB = sizes[iB + NatomsA];
      int offsetB = offset[NatomsA+iB];

      Matrix Tab_block = DampedTabs[iB*NatomsA+iA].GetBlock(NpolsA,NpolsB,1,1);
      Z.SetBlock(Tab_block,offsetA,offsetB);
      Tab_block.Transpose(); // Tba = Transpose(Tab)
      Z.SetBlock(Tab_block,offsetB,offsetA);
    }
  }

  // Build potential V0 due to permanent multipoles
  // For an atom "a" on MonA:
  // V0 = sum(atoms "b" on MonB) Tab_tu * Qb_u
  for (int iA=0;iA<NatomsA;iA++) {
    int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
    Vector tmp(nQA);
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
      // Tab*QB
      tmp += DampedTabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
    }
    
    // Now grab the relevant elements of the tmp vector and store them in V0
    int nelem = sizes[iA];
    int offsetA = offset[iA];
    for (int t=0;t<nelem;t++) {
      V0[t+offsetA] = tmp[t+1];
    }
  }
  // Repeat for potential felt at atom "b" on MonB due to MonA atoms
  for (int iB=0;iB<NatomsB;iB++) {
    int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
    Vector tmp(nQB);
    for (int iA=0;iA<NatomsA;iA++) {

      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      // Tba*QA = Tab'*QA
      Matrix Tab = DampedTabs[iB*NatomsA+iA];
      Tab.Transpose();
      tmp += Tab.MatrixTimesVector(QA.GetMoments());
     
    }
    
    // Now grab the relevant elements of the tmp vector and store them in V0
    int nelem = sizes[NatomsA+iB];
    int offsetB = offset[NatomsA+iB];
    for (int t=0;t<nelem;t++) {
      V0[t+offsetB] = tmp[t+1];
    }
  }


  // Now solve Z*dQ = -V0 for the induced moments dQ --> dQ =
  // Z^(-1)*V0.  
  Vector V0_copy = V0; // create a copy of V0 which gets destroyed
		       // when we solve for dQ
  V0_copy.Scale(-1.0);
  Z.SolveLinearEquations(V0_copy); // Solve Z*dQ=-V0 (via LU decomposition)
  dQ = V0_copy;

  // Compute the induction energy from these induced moments.
  // Eind = 0.5*dQ*V0
  Eind = 0.5*dQ.DotProduct(V0);
  //printf("Eind = %f\n",Eind*HartreesToKJpermole);
  E_Induction_MM = Eind;



  // Transfer the sparse list of induced multipole moments to
  // non-sparse one.
  Multipole *dQA = new Multipole[NatomsA];
  Multipole *dQB = new Multipole[NatomsB];
  for (int iA=0;iA<NatomsA;iA++) {
    int rank = MonA.GetAtom(iA).GetMultipoleMoments().GetRank();
    int offsetA = offset[iA];
    dQA[iA].Initialize(rank);
    for (int t=0;t<sizes[iA];t++) {
      dQA[iA](t+1) = dQ[offsetA+t];  // t+1 on LHS because always skip charge
    }
  }
  for (int iB=0;iB<NatomsB;iB++) {
    int rank = MonB.GetAtom(iB).GetMultipoleMoments().GetRank();
    int offsetB = offset[NatomsA+iB];
    dQB[iB].Initialize(rank);
    for (int t=0;t<sizes[iB+NatomsA];t++) {
      dQB[iB](t+1) = dQ[offsetB+t];  // t+1 on LHS because always skip charge
    }
  }


  // Optionally print out the induced multipole moments
  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles ***\n");
    for (int iA=0;iA<NatomsA;iA++) {
      char label[50];
      sprintf(label,"Monomer 1, %s%d",MonA.GetAtom(iA).GetSymbol().c_str(),
	      MonA.GetAtom(iA).GetAtomIndex());
      string str = label;
      
      
      dQA[iA].Print(str);
    }
    for (int iB=0;iB<NatomsB;iB++) {
      char label[50];
      sprintf(label,"Monomer 2, %s%d",MonB.GetAtom(iB).GetSymbol().c_str(),
	      MonB.GetAtom(iB).GetAtomIndex());
      string str = label;
      dQB[iB].Print(str);
    }
  }

  //printf("EindA = %f, EindB=%f, Eind = %f kJ/mol\n",EindA*HartreesToKJpermole,EindB*HartreesToKJpermole,Eind*HartreesToKJpermole);

  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > 0)
      printf("Dimer induction energy wall time = %.5f seconds\n",elapsed_time);
  
  
  delete [] dQA;
  delete [] dQB;

  return Eind;
}

// Build the nuclear gradient of the Tab matrix for interaction
// between 2 monomer multipoles.  Compute the forces due to
// electrostatics and induction.
Vector Dimer::ComputeMultipoleGradient() {

  Vector Grad(3*Natoms);

  // Variables to store total permanent electrostatic, induction & total energies
  double Ees = 0.0, Eind = 0.0, Etot = 0.0;

  printf("\nEvaluating ab initio force field energy & gradient for dimer (%d,%d)\n",
	   indexA,indexB);

  // Step 1: Build the gradient of the geometric interaction matrix Tab.  

  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor for this dimer

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  // In case these have already been allocated, delete them.
  // If not, we are just deleting NULL, which is harmless.
  delete [] Tabs;
  delete [] DampedTabs;
  delete [] TabsGrad;
  delete [] DampedTabsGrad;


  Tabs = new Matrix[NatomsA*NatomsB];
  DampedTabs = new Matrix[NatomsA*NatomsB];

  TabsGrad = new Matrix[NatomsA*NatomsB];
  DampedTabsGrad = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      /***Construct the Tab geometric interaction matrices and their
	  nuclear gradients ***/
      Matrix Tab(nQA,nQB);

      // Compute the undamped version of Tab.  "-999.0" is a dummy
      // parameter to turn damping off.
      Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      Tabs[iB*NatomsA+iA].Initialize(Tab);  

      // Compute the damped version of Tab and store it.
      Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, beta_damp);
      DampedTabs[iB*NatomsA+iA].Initialize(Tab);  


      /*** Now do the gradients of Tab and store them ***/
      Matrix Tab_grad(nQA*nQB,6);
      Tab = Tabs[iB*NatomsA+iA];

      // Compute the undamped version of dTab/dX.  "-999.0" is a dummy
      // parameter to turn damping off.
      Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotVecA,
                    RotAngA, MonB.GetAtom(iB), RotVecB, RotAngB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      TabsGrad[iB*NatomsA+iA].Initialize(Tab_grad);  


      // Compute the damped version of dTab/dX and store it.
      Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotVecA,
                     RotAngA, MonB.GetAtom(iB), RotVecB, RotAngB, beta_damp);
      DampedTabsGrad[iB*NatomsA+iA].Initialize(Tab_grad);  
    }

  // Step 2a: Compute the electrostatic interaction forces
  // Step 2a: Permanent multipole contributions
  //          dE/dX = sum(tu) QA(t)*dT(t,u)/dX*QB(u) 

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());

      int nQA = QA.GetLength();
      int nQB = QB.GetLength();

     
      // Evaluate the energy
      Vector tmp(nQA);
      // Tab*QB
      tmp = Tabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
      // QA*(Tab*QB)
      double Energy = QA.GetMoments().DotProduct(tmp);
      Ees += Energy;


      // And the gradient contribution:
      Matrix dTdX(TabsGrad[iB*NatomsA+iA]);
      /*
      if (iA==0 && iB==0) {
	//dTdX.Print("\ndTdX for O-O interaction");
	printf("\ndT/dX for O-O interaction\n");
	printf(" t   u      dAx        dAy        dAz        dBx        dBy        dBz\n");

	string types[25];
	types[0]="00";
	types[1]="1x"; types[2]="1y"; types[3]="1z";
	types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; 
	types[8]="22s";
	types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; 
	types[13]="32s"; 
	types[14]="33c"; types[15]="33s";
	types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; 
	types[20]="42s"; 
	types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";
	
	for (int u=0;u<nQB;u++) {
	  for (int t=0;t<nQA;t++) {
	    int tu = u*nQA + t;
	    printf("%3s %3s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		   types[t].c_str(),types[u].c_str(),dTdX(tu,0),dTdX(tu,1),
		   dTdX(tu,2),dTdX(tu,3), dTdX(tu,4),dTdX(tu,5));	    
	  }
	}
	printf("\n");
      }
      */

      // Create temporary storage for the gradient contribution
      // arising from a single pair of interacting atoms.
      Vector tmpGrad(6); // x1,y1,z1 and x2,y2,z2
      tmpGrad.Set();

      for (int X=0;X<6;X++) 
	for (int t=0;t<nQA;t++)
	  for (int u=0;u<nQB;u++) {
	    int tu = u*nQA + t; // define a compound index
	    tmpGrad[X] += QA(t)*dTdX(tu,X)*QB(u);
	  }

      // Incorporate this contribution into the full dimer gradient.
      NatomsA = MonA.GetNumberOfAtoms();
      for (int i=0;i<3;i++) {
	Grad[iA*3+i]  += tmpGrad[i]; // atom 1
	Grad[3*NatomsA+iB*3+i]  += tmpGrad[i+3]; // atom 2
	}
    }

  // We need to save this gradient contribution to build the full
  // cluster gradient later.
  Grad_Electrostatic = Grad;


  // Step 2b: Induced 2-body multipole gradient contributions
  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction nuclear gradient for dimer (%d,%d)\n",indexA,indexB);
  }
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction gradient.\n",beta_damp);

  // Initialize a few variables
  bool iterate = true;
  int cycle = 0;

  double ind_conv = Params::Parameters().GetInductionConvergence();
  double ind_gradconv = Params::Parameters().GetInductionGradConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  double Gconv = 1.0/pow(10,ind_gradconv);
  if (Params::Parameters().PrintLevel() > 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f, ind_gradconv = %12.6f, Gconv = %12.6f\n", ind_conv,Econv,ind_gradconv,Gconv );

  double Eind_old = 0.0;

  Vector Eind_grad(3*Natoms);
  Vector old_Eind_grad(3*Natoms);

  // Create lists that will store Induced Moments for each atom on A
  // and B and initialize them to the proper sizes based on the higher
  // rank object: polarizabilities or multipoles.  This is important:
  // CamCasp sometimes prints polarizabilities with high ranks & lots
  // of zeros, even if a lower rank would do.  We need two copies of
  // each list.
  Multipole *dQA = new Multipole[NatomsA];
  Multipole *old_dQA = new Multipole[NatomsA];

  // For the gradient, we have one multipole vector for each atom on a 
  // monomer.  Each of those vectors has 3*Natoms gradient components.
  // Overall, we have (3*Natoms)*NatomsA multipole vectors.
  Multipole *dQA_grad = new Multipole[3*Natoms*NatomsA];
  Multipole *old_dQA_grad = new Multipole[3*Natoms*NatomsA];

  for (int iA=0;iA<NatomsA;iA++) {
    int Rmom = MonA.GetAtom(iA).GetMultipoleMoments().GetRank();
    int Rpol = MonA.GetAtom(iA).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQA[iA].Initialize(max_rank);
    old_dQA[iA].Initialize(max_rank);

    for (int X=0;X<3*Natoms;X++) {
      dQA_grad[3*Natoms*iA + X].Initialize(max_rank);
      old_dQA_grad[3*Natoms*iA + X].Initialize(max_rank);
    }
  }

  Multipole *dQB = new Multipole[NatomsB];
  Multipole *old_dQB = new Multipole[NatomsB];

  Multipole *dQB_grad = new Multipole[3*Natoms*NatomsB];
  Multipole *old_dQB_grad = new Multipole[3*Natoms*NatomsB];


  for (int iB=0;iB<NatomsB;iB++) {
    int Rmom = MonB.GetAtom(iB).GetMultipoleMoments().GetRank();
    int Rpol = MonB.GetAtom(iB).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQB[iB].Initialize(max_rank);
    old_dQB[iB].Initialize(max_rank);
    for (int X=0;X<3*Natoms;X++) {
      dQB_grad[3*Natoms*iB + X].Initialize(max_rank);
      old_dQB_grad[3*Natoms*iB + X].Initialize(max_rank);
    }
  }

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
    double EindA = 0.0,  EindB = 0.0;
    old_Eind_grad = Eind_grad;
    Eind_grad.Set();

    for (int iA=0;iA<NatomsA;iA++) {
      old_dQA[iA] = dQA[iA];
      dQA[iA].Set();
      
      for (int X=0;X<3*Natoms;X++) {
	old_dQA_grad[3*Natoms*iA+X] = dQA_grad[3*Natoms*iA+X];
	dQA_grad[3*Natoms*iA+X].Set();
      }
    }
    for (int iB=0;iB<NatomsB;iB++) {
      old_dQB[iB] = dQB[iB];
      dQB[iB].Set();

      for (int X=0;X<3*Natoms;X++) {
	old_dQB_grad[3*Natoms*iB+X] = dQB_grad[3*Natoms*iB+X];
	dQB_grad[3*Natoms*iB+X].Set();
      }
    }

    // Induce multipoles
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      // Grab some data we need
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Polarizability PolA(MonA.GetAtom(iA).GetPolarizability(),true);
      int nQA = QA.GetLength();
      int NpolA = PolA.GetLength();
      int NmomA = QA.GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	// Grab some data we need
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	Polarizability PolB(MonB.GetAtom(iB).GetPolarizability(),true);
	int nQB = QB.GetLength();
	int NpolB = PolB.GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	int dimT1 = Tab.GetRows();
	int dimT2 = Tab.GetCols();

	/*
	if (iA==0 && iB==0 &&) {
	  Matrix dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	  dTdX.Scale(HartreesToKJpermole);
	  printf("\nDamped dT/dX for O-O interaction\n");
	  printf(" t   u      dAx        dAy        dAz        dBx        dBy        dBz\n");

	  string types[25];
	  types[0]="00";
	  types[1]="1x"; types[2]="1y"; types[3]="1z";
	  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; 
	  types[8]="22s";
	  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; 
	  types[13]="32s"; 
	  types[14]="33c"; types[15]="33s";
	  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; 
	  types[20]="42s"; 
	  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";
	  
	  for (int u=0;u<nQB;u++) {
	    for (int t=0;t<nQA;t++) {
	      int tu = u*nQA + t;
	      printf("%3s %3s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		     types[t].c_str(),types[u].c_str(),dTdX(tu,0),dTdX(tu,1),
		     dTdX(tu,2),dTdX(tu,3), dTdX(tu,4),dTdX(tu,5));	    
	    }
	  }
	  printf("\n");
	}
	*/

	// Grab dTdX, and map it from (Qa*Qb,6) -> (Qa*Qb,3*Natoms);
	Matrix small_dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	int rows = small_dTdX.GetRows();
	Matrix dTdX(rows,3*Natoms);
	for (int irows=0;irows<rows;irows++) 
	  for (int X=0;X<3;X++) {
	    dTdX(irows,iA*3+X) = small_dTdX(irows,X);
	    dTdX(irows, 3*NatomsA+iB*3+X) = small_dTdX(irows,X+3); 
	  }

	// Induce multipoles on monomer A due to monomer B
	// dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))
	for (int a=0;a<NpolA;a++) // loop over elements of dQA
	  for (int t=0;t<min(NpolA,dimT1);t++) 
	    for (int u=0;u<min(dimT2,NmomB);u++) {
	      // Induced moments:
	      dQA[iA](a) -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQB[iB](u));

	      // Gradient of induced moments:
	      int tu = u*nQA + t; // define a compound index for (t,u)
	      for (int X=0;X<3*Natoms;X++) {
		dQA_grad[3*Natoms*iA + X](a) -= PolA(a,t)* 
		  (dTdX(tu,X)*(QB(u) + old_dQB[iB](u)) 
		   + Tab(t,u)*old_dQB_grad[3*Natoms*iB + X](u));
	      }
	    }

	// Induce multipoles on monomer B due to monomer A
	// dQB(a) = dQB(a) - polB(a,t)*Tab(u,t)*(QA(u) + old_dQA(u))
	for (int a=0;a<NpolB;a++) // loop over elements of dQA
	  for (int u=0;u<min(NpolB,dimT2);u++) 
	    for (int t=0;t<min(dimT1,NmomA);t++) {
	      // Induced moments:
	      dQB[iB](a) -= PolB(a,u)*Tab(t,u)*(QA(t) + old_dQA[iA](t));

	      // Gradient of induced moments:
	      int tu = u*nQA + t; // define a compound index for (t,u)
	      for (int X=0;X<3*Natoms;X++) {
		dQB_grad[3*Natoms*iB + X](a) -= PolB(a,u)* 
		  (dTdX(tu,X)*(QA(t) + old_dQA[iA](t)) 
		   + Tab(t,u)*old_dQA_grad[3*Natoms*iA + X](t));
	      }
	    }
      }
    }

    // Compute the induction energy/gradient using the induced multipoles
    // determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int nQA = QA.GetLength();
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int nQB = QB.GetLength();
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);

	// Grab dTdX, and map it from (Qa*Qb,6) -> (Qa*Qb,3*Natoms);
	Matrix small_dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	int rows = small_dTdX.GetRows();
	Matrix dTdX(rows,3*Natoms);
	for (int X=0;X<3;X++) 
	  for (int irows=0;irows<rows;irows++) {
	    dTdX(irows,iA*3+X) = small_dTdX(irows,X);
	    dTdX(irows, 3*NatomsA+iB*3+X) = small_dTdX(irows,X+3); 
	  }
	    
	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    // Energy: 
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);

	    // Gradient:
	    int tu = u*nQA + t; // define a compound index for (t,u)
	    for (int X=0;X<3*Natoms;X++) {
	      // Contribution due to dQA
	      Eind_grad[X] += 0.5*dQA_grad[3*Natoms*iA + X](t)*Tab(t,u)*QB(u)
		+ 0.5*dQA[iA](t)*dTdX(tu,X)*QB(u);

	      // Contribution due to dQB
	      Eind_grad[X] += 0.5*QA(t)*Tab(t,u)*dQB_grad[3*Natoms*iB + X](u)
		+ 0.5*QA(t)*dTdX(tu,X)*dQB[iB](u);
	    }
	  }
      }
    }
    Eind = EindA + EindB;
  
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
  }

  // Store the energy contibutions, in hartrees.
  E_Electrostatic_MM = Ees;
  E_Induction_MM = Eind;
  Etot = Ees + Eind;
  Energy_MM = Etot;

  if (Params::Parameters().PrintLevel() > 0) {
    // Print out a summary of the classical energy just calculated
    printf("--------------------------------------------------------\n");
    printf(" Total Classical Interaction Energy for Dimer (%d,%d)\n",
	   indexA,indexB);
    printf("--------------------------------------------------------\n");
    printf("      Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
    printf("      Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
    printf("    ---------------------------------------\n");
    printf("      Total         = %12.6f kJ/mol\n",Etot*HartreesToKJpermole);
    printf("--------------------------------------------------------\n");
  }

  // Return the gradient
  //Eind_grad.Scale(HartreesToKJpermole);
  //Eind_grad.PrintGradient("Induction gradient");
  //Eind_grad.Scale(1.0/HartreesToKJpermole); 

  //Grad.Scale(HartreesToKJpermole); 
  //Grad.PrintGradient("\nElectrostatic Gradient");
  //Grad.Scale(1.0/HartreesToKJpermole); 

  // Combine the permanent and induced gradient contributions
  Grad += Eind_grad;

  //Grad.PrintGradient("Dimer gradient\n");

  delete [] dQA;
  delete [] old_dQA;
  delete [] dQB;
  delete [] old_dQB;

  delete [] dQA_grad;
  delete [] old_dQA_grad;
  delete [] dQB_grad;
  delete [] old_dQB_grad;


  return Grad;
}

void Dimer::BuildDampedTabInteractionMatrices() {

  //printf("\nBuilding Damped Tab Interaction matrices for dimer (%d,%d)\n",
  //	   indexA,indexB);

  // Get the damping factor
  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor for this dimer

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  DampedTabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, beta_damp);
      DampedTabs[iB*NatomsA+iA].Initialize(tmpTab);  
      /* if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole);
	tmpTab.Print("Tab - damped");
      } */
    }
}

void Dimer::BuildTabInteractionMatrices() {

  //printf("\nBuilding Damped Tab Interaction matrices for dimer (%d,%d)\n",
  //	   indexA,indexB);

  // Get the damping factor
  double beta_damp = -999.0; // damping factor for this dimer

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  Tabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Vector RotVecA(MonA.GetRotationVector());
      double RotAngA = MonA.GetRotationAngle();
      Vector RotVecB(MonB.GetRotationVector());
      double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
						       MonB.GetAtom(iB),
						       RotVecB,RotAngB, beta_damp);
      Tabs[iB*NatomsA+iA].Initialize(tmpTab);  
      /* if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole);
	tmpTab.Print("Tab - damped");
      } */
    }
}


double Dimer::ComputeInteratomicMP2TwoBodyDispersionCorrection() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp = 0.0;

  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    // Get the pure C6 coefficients for atom I
    double C6i_uchf = AtomI.GetC6Coefficient("UCHF");
    double C6i_new = AtomI.GetC6Coefficient();

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j); 

      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get the pure C6 coefficients for atom J
      double C6j_uchf = AtomJ.GetC6Coefficient("UCHF");
      double C6j_new = AtomJ.GetC6Coefficient();

      // Compute mixed C6 coefficients for this atom pair
      // C6ij = sqrt(C6i*C6j)
      double C6_uchf = sqrt(C6i_uchf*C6j_uchf);
      double C6_new = sqrt(C6i_new*C6j_new);

      // Get the interatomic distance
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Compute the short-range damping function
      // Get the van der Waals radii (in bohr)
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");

      bool DoTangToenniesDamping = false;// if false, use Aziz
      double F6ij;
      if (DoTangToenniesDamping) {
	// Get Tang-Toennies-type 3-body damping factor
	// Use empirical expression for beta: J Chem Phys 131, 094106 (2009).
	double betaIJ = -0.334*(Di+Dj) + 4.386;
	//printf("betaIJ = %f\n",betaIJ);
	F6ij = AtomI.TangToenniesDampingFactor(6,betaIJ,Rij);
      }
      else {// Aziz damping
	double S = 1.4;
	F6ij = exp(-(S*(Di+Dj)/Rij-1)*(S*(Di+Dj)/Rij-1));
      }

      printf("Dimer (%d,%d)... Atom pair (%d,%d) contrib = %f\n",indexA, indexB, 
	     i, j, -F6ij*(C6_new - C6_uchf)/pow(Rij,6.0));

      printf("(%d,%d): C6(uchf) = %f, C6(cks) = %f, Rij = %f, Dij = %f,  F6ij = %f\n",i,j,C6_uchf,C6_new,Rij,Di+Dj,F6ij);
      // Compute the dispersion contribution
      E2b_disp += -F6ij*(C6_new - C6_uchf)/pow(Rij,6.0);
      //printf("  contrib = %f\n",-F6ij*(C6_new - C6_uchf)/pow(Rij,6.0));
    }
  }
  printf("Dimer (%d,%d) dispersion correction: %f\n",indexA,indexB,E2b_disp);
  return E2b_disp;
}

// This function uses a list of predetermined c6 coeffs.  You should
// generally use the ab initio dispersion model found in
// Dimer::ComputeTwoBodyDispersion() instead.
double Dimer::EstimateTwoBodyDispersion() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      //double Rij = MonA.GetAtom(i).GetInterAtomicDistance(MonB.GetAtom(j))*AngToBohr;
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
      //printf("Rij = %f\n",Rij);
      
      // Get the C6 coefficient for this pair of atoms
      double C6ij = AtomI.EstimateC6Coefficient(AtomJ);


      // Compute the short-range damping function
      double damping = 1.0; // no damping
      if (Params::Parameters().DoDispersionDamping()) {      
	
	// Get the van der Waals diameters for the two atoms
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	
	// Get Tang-Toennies-type 2-body damping factor
	// Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	double beta = -0.33*(Di+Dj) + 4.39;
	damping = AtomI.TangToenniesDampingFactor(6,beta,Rij);
      }

      // E(2b) += - damping * C6ij/Rij^6 
      E2b_disp -= damping * C6ij / (pow(Rij,6));
    }
  }
  //printf("Approximate 2-body dispersion for dimer (%d,%d): %f in kcal/mol\n",indexA, indexB,
  // E2b_disp*HartreesToKcalpermole);
  return E2b_disp;
}

// Does C6 and C8 2-body dispersion.  Can do C10 as well, but that
// term is small and hard to get right, so we disable it.  The
// dispersion coefficients are calculated via ab initio calculation of
// the frequency-dependent polarizabilities and Casimir-Polder
// integration.
double Dimer::ComputeTwoBodyDispersion() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Compute the C6 & C8 dispersion terms via Casimir-Polder integration
      double C6ij = AtomI.CasimirC6Coefficient(AtomJ);
      double C8ij = AtomI.CasimirC8Coefficient(AtomJ);
      // We used to get C10 as well, but that contributes little
      // and is hard to get right.
      //double C10ij = AtomI.CasimirC10Coefficient(AtomJ);

      // printf("Dispersion coefficients for atom pair (%d,%d): C6 = %f, C8=%f\n",i,j,C6ij,C8ij);


      /*
      // Get the C6 coefficient for this pair of atoms

      double C6i = AtomI.GetC6Coefficient();
      double C6j = AtomJ.GetC6Coefficient();
      double Poli = AtomI.GetIsotropicDipolePolarizability();
      double Polj = AtomJ.GetIsotropicDipolePolarizability();

      // Use Tang (Phys Rev 177, 108-114, 1969) combination rule for
      // C6 coefficients.  Or see Stone's Theory of IMF, eqn 4.3.14
      // (pg 61 in the 2002 edition).
      //double C6ij = 2*Poli*Polj*C6i*C6j / (Poli*Poli*C6j + Polj*Polj*C6i);
      
      // Get the C8 coefficient for this pair of atoms
      double C8i = AtomI.GetC8Coefficient();
      double C8j = AtomJ.GetC8Coefficient();
      //printf("C8i=%f, C8j=%f\n",C8i,C8j);
      // Just use crude combination rule, Cab = sqrt(Caa*Cbb)
      // handle signs, in case coefficient is negative
      double sign = 1.0;
      if (C8i < 0.0) sign *= -1.0;
      if (C8j < 0.0) sign *= -1.0;
      double C8ij = sqrt( fabs(C8i) * fabs(C8j) );

      // Get the C10 coefficient for this pair of atoms
      double C10i = AtomI.GetC10Coefficient();
      double C10j = AtomJ.GetC10Coefficient();
      //printf("C10i=%f, C10j=%f\n",C10i,C10j);
      // Just use crude combination rule, Cab = sqrt(Caa*Cbb)
      // handle signs, in case coefficient is negative
      sign = 1.0;
      if (C10i < 0.0) sign *= -1.0;
      if (C10j < 0.0) sign *= -1.0;
      double C10ij = sqrt( fabs(C10i) * fabs(C10j) );
      //double C10ijnew = AtomI.CasimirC10Coefficient(AtomJ);

      // printf("C6ij  = %f (based on %f and %f)\n",C6ij,C6i,C6j);
      // printf("C8ij  = %f (based on %f and %f)\n",C8ij,C8i,C8j);
      // printf("C10ij = %f (based on %f and %f)\n",C10ij,C10i,C10j);
      */
      
      // Compute the short-range damping function
      double damping6 = 1.0; // no damping
      double damping8 = 1.0; // no damping
      if (Params::Parameters().DoDispersionDamping()) {   
	// Get Tang-Toennies-type 2-body damping factor
	// Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	// Frankly, it doesn't matter.  With typical values of the following,
	// the damping is already less than 0.1% by 4-5 Angstroms, which is
	// well-beyond our quantum cutoff.  
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	double beta = -0.33*(Di+Dj) + 4.39;
	damping6 = AtomI.TangToenniesDampingFactor(6,beta,Rij);
	damping8 = AtomI.TangToenniesDampingFactor(8,beta,Rij);
	//double damping10 = AtomI.TangToenniesDampingFactor(10,beta,Rij);
	//printf("beta = %f, Rij = %f, damping 6/8/10 = %f, %f, %f\n",beta,Rij,damping6,damping8,damping10);
      }

      // E(2b) += - damping6 * C6ij/Rij^6 - damping8 * C8ij/Rij^8 - damping10 * C10ij/Rij^10
      double d6 = -1.0*damping6 * C6ij / (pow(Rij,6));
      double d8 = -1.0*damping8 * C8ij / (pow(Rij,8));
      //double d10 = -1.0*damping10 * C10ij / (pow(Rij,10));

      //E2b_disp += d6 + d8 + d10;
      E2b_disp += d6 + d8;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp*HartreesToKJpermole);
  }

  // Store the energy, in hartrees
  E_2b_disp_MM = E2b_disp;

  return E2b_disp;

}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C6.
double Dimer::ComputeTwoBodyDispersionC6() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C6 = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Get the C6 coefficient for this pair of atoms
      double C6ij = AtomI.CasimirC6Coefficient(AtomJ);    

      // Compute the short-range damping function
      double damping6 = 1.0; // no damping
      if (Params::Parameters().DoDispersionDamping()) {   
	// Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	// Frankly, it doesn't matter.  With typical values of the following,
	// the damping is already less than 0.1% by 4-5 Angstroms, which is
	// well-beyond our quantum cutoff.  
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");  
	double beta = -0.33*(Di+Dj) + 4.39;
	damping6 = AtomI.TangToenniesDampingFactor(6,beta,Rij);
      }
	
      // E(2b) += - damping6 * C6ij/Rij^6 - damping8 * C8ij/Rij^8 - damping10 * C10ij/Rij^10
      double d6 = -1.0*damping6 * C6ij / (pow(Rij,6));
      E2b_disp_C6 += d6;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C6 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C6*HartreesToKJpermole);
  }
  return E2b_disp_C6;
}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C8.
double Dimer::ComputeTwoBodyDispersionC8() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C8 = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
      
      // Get the C8 coefficient for this pair of atoms
      double C8ij = AtomI.CasimirC8Coefficient(AtomJ);
      
      // Compute the short-range damping function
      double damping8 = 1.0; // no damping
      if (Params::Parameters().DoDispersionDamping()) {   
	// Get Tang-Toennies-type 2-body damping factor
	// Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	// Frankly, it doesn't matter.  With typical values of the following,
	// the damping is already less than 0.1% by 4-5 Angstroms, which is
	// well-beyond our quantum cutoff.  
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	double beta = -0.33*(Di+Dj) + 4.39;
	double damping8 = AtomI.TangToenniesDampingFactor(8,beta,Rij);
      }

      // E(2b) += - damping8 * C8ij/Rij^8 
      double d8 = -1.0*damping8 * C8ij / (pow(Rij,8));
      E2b_disp_C8 += d8 ;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C8 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C8*HartreesToKJpermole);
  }
  return E2b_disp_C8;
}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C10.
double Dimer::ComputeTwoBodyDispersionC10() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C10 = 0.0;
  
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
  
      // Get the C10 coefficient for this pair of atoms
      double C10ij = AtomI.CasimirC10Coefficient(AtomJ);

      // Compute the short-range damping function
      double damping10 = 1.0; // no damping
      if (Params::Parameters().DoDispersionDamping()) {   
	// Get Tang-Toennies-type 2-body damping factor
	// Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	// Frankly, it doesn't matter.  With typical values of the following,
	// the damping is already less than 0.1% by 4-5 Angstroms, which is
	// well-beyond our quantum cutoff.  
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	double beta = -0.33*(Di+Dj) + 4.39;
	damping10 = AtomI.TangToenniesDampingFactor(10,beta,Rij);
      }

      // E(2b) += - damping10 * C10ij/Rij^10
      double d10 = -1.0*damping10 * C10ij / (pow(Rij,10));
      E2b_disp_C10 += d10;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C10 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C10*HartreesToKJpermole);
  }
  return E2b_disp_C10;
}

// Does Axilrod-Teller-Muto 3-body dispersion, albeit with approximate
// tabulated C9 coefficients.
double Dimer::EstimateThreeBodyDispersion() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E3b_disp = 0.0;

  // Do as two cases: 

  /* Case 1: 2 atoms from MonA, 1 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=i+1;j<NatomsA;j++) {
      Atom AtomJ = MonA.GetAtom(j);
      for (int k=0;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

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
	double C9ijk = AtomI.EstimateC9Coefficient(AtomJ, AtomK, "Tkatchenko");
	
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
	//printf("Daab(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }

  /* Case 2: 1 atom from MonA, 2 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=0;j<NatomsB;j++) {
      Atom AtomJ = MonB.GetAtom(j);
      for (int k=j+1;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

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
	double C9ijk = AtomI.EstimateC9Coefficient(AtomJ, AtomK, "Tkatchenko");
	
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
	//printf("Dabb(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }
  return E3b_disp;

}


// Does Axilrod-Teller-Muto 3-body dispersion with ab initio C9
// coefficients computed via Casimir-Polder integration of the
// frequency-dependent polarizabilities.  The contributions here
// involve 2 atoms on one monomer and 1 on the other.
double Dimer::ComputeThreeBodyDispersion() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E3b_disp = 0.0;

  // Do as two cases: 

  /* Case 1: 2 atoms from MonA, 1 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=i+1;j<NatomsA;j++) {
      Atom AtomJ = MonA.GetAtom(j);
      for (int k=0;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;

	// Get the C9 coefficient
	double C9ijk = AtomI.CasimirC9Coefficient(AtomJ, AtomK);
	
	// Compute the short-range damping function
	double damping = 1.0; // no damping
	if (Params::Parameters().DoDispersionDamping()) {   
	  // Get the van der Waals diameters
	  // Note: these are empirical and should be replaced with
	  // something more rigorous.  In this crude version, the
	  // DispersionAtomType parameters should have already been set
	  // when computing the 2-body dispersion.
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
	//printf("Daab(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }

  /* Case 2: 1 atom from MonA, 2 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=0;j<NatomsB;j++) {
      Atom AtomJ = MonB.GetAtom(j);
      for (int k=j+1;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

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
	double C9ijk = AtomI.CasimirC9Coefficient(AtomJ, AtomK);


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
	//printf("Dabb(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }
  return E3b_disp;

}
