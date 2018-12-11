// for getpid:
#include <sys/types.h>
#include <unistd.h>
// other:
#include "quasiharmonic.h"
#include "cluster.h"
#include "supercell.h"
#include "dlf_interface.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
using namespace hmbi_constants;
//using namespace std;

Quasiharmonic::Quasiharmonic() {

  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();

  volume_ref = 0.0;
  a_ref = 0.0;
  b_ref = 0.0;
  c_ref = 0.0;
  alpha_ref = 0.0;
  beta_ref = 0.0;
  gamma_ref = 0.0;
  //Gruneisen_init = false;
  Coords_ref.Initialize(Cluster::cluster().GetSymmetryUniqueCoordinates().GetLength());
  Modes_ref.Initialize(3*Natoms,3*Natoms);
  ref_init = false;

  if(!Params::Parameters().IsSupercellJob()){
    total_k_point = 1;
    freq_ref.Initialize(3*Natoms);
    if(Params::Parameters().UseVolumeGruneisen())
      Gruneisen_parameters.Initialize(3*Natoms);
    else{
      Grun_a.Initialize(3*Natoms);
      Grun_b.Initialize(3*Natoms);
      Grun_c.Initialize(3*Natoms);
    }
  }else{

    total_k_point = Supercell::supercell().GetTotalKPoint();
    freq_ref.Initialize(3*Natoms*total_k_point);
    freq_ref_k.Initialize(3*Natoms);
    if(Params::Parameters().UseVolumeGruneisen())
      Gruneisen_parameters.Initialize(3*Natoms*total_k_point);
    else{
      Grun_a.Initialize(3*Natoms*total_k_point);
      Grun_b.Initialize(3*Natoms*total_k_point);
      Grun_c.Initialize(3*Natoms*total_k_point);
    }
  }

  hess_type = "";
}

//destruction operator
Quasiharmonic::~Quasiharmonic(){

}

void Quasiharmonic::Initialize(){

  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();

  //reference value
  Coords_ref.Initialize(Cluster::cluster().GetSymmetryUniqueCoordinates().GetLength());
  a_ref = 0.0;
  b_ref = 0.0;
  c_ref = 0.0;
  alpha_ref = 0.0;
  beta_ref = 0.0;
  gamma_ref = 0.0;
  volume_ref = 0.0;
  //Gruneisen_init = false;
  Modes_ref.Initialize(3*Natoms,3*Natoms);
  ref_init = false;

  if(!Params::Parameters().IsSupercellJob()){
    total_k_point = 1;
    freq_ref.Initialize(3*Natoms);
    if(Params::Parameters().UseVolumeGruneisen())
      Gruneisen_parameters.Initialize(3*Natoms);
    else{
      Grun_a.Initialize(3*Natoms);
      Grun_b.Initialize(3*Natoms);
      Grun_c.Initialize(3*Natoms);
    }
  }else{
    total_k_point = Supercell::supercell().GetTotalKPoint();
    freq_ref.Initialize(3*Natoms*total_k_point);
    freq_ref_k.Initialize(3*Natoms);
    if(Params::Parameters().UseVolumeGruneisen())
      Gruneisen_parameters.Initialize(3*Natoms*total_k_point);
    else{
      Grun_a.Initialize(3*Natoms*total_k_point);
      Grun_b.Initialize(3*Natoms*total_k_point);
      Grun_c.Initialize(3*Natoms*total_k_point);
    }
  }

  hess_type = "";

  //creating directories for each quasiharmonic step if saving them
  if(Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().Params::Parameters().AreQHFAvailable()){
    Cluster::cluster().CreateWorkingDirectories("reference",0);
    Cluster::cluster().CreateWorkingDirectories("reference",1);      
    
    if(Params::Parameters().UseVolumeGruneisen()){
      Cluster::cluster().CreateWorkingDirectories("PlusFreq",0);
      Cluster::cluster().CreateWorkingDirectories("PlusFreq",1);
      Cluster::cluster().CreateWorkingDirectories("MinusFreq",0);
      Cluster::cluster().CreateWorkingDirectories("MinusFreq",1);
    }else{
      Cluster::cluster().CreateWorkingDirectories("PlusFreqA",0);
      Cluster::cluster().CreateWorkingDirectories("PlusFreqA",1);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqA",0);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqA",1);
      Cluster::cluster().CreateWorkingDirectories("PlusFreqB",0);
      Cluster::cluster().CreateWorkingDirectories("PlusFreqB",1);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqB",0);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqB",1);
      Cluster::cluster().CreateWorkingDirectories("PlusFreqC",0);
      Cluster::cluster().CreateWorkingDirectories("PlusFreqC",1);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqC",0);
      Cluster::cluster().CreateWorkingDirectories("MinusFreqC",1);
    }
  }
}

void Quasiharmonic::DetermineReferenceFrequencies(ifstream &infile,int &imag_num)  {

  //name of file reference input file is stored and name of the directory (within the
  //respective QM/MM directory) where the terms are stored 
  hess_type = "reference";

  //change job type to Optimization and reread the qc_rem
  Params::Parameters().SetParameter("JOBTYPE","Opt");
  string qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
  Params::Parameters().SetParameter("QC_REM",qc_rem);
  int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  int nvarin, nvarin2, nspec;
  if(Params::Parameters().UseFullQMOnly()){
   nvarin =  3*UniqueAtoms+9;
   nvarin2 = 5*UniqueAtoms+15;
   nspec = 3*UniqueAtoms+9; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }
  else{
   nvarin =  3*UniqueAtoms+6;
   nvarin2 = 5*UniqueAtoms+10;
   nspec = 3*UniqueAtoms+6; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }

  /*int nvarin =  3*UniqueAtoms+6;
  int nvarin2 = 5*UniqueAtoms+10;// nframe*nat*3 + nweight + nmass + n_po_scaling
  int nspec = 2*UniqueAtoms+4;// nspec= nat + nz + 5*ncons + 2*nconn*/ //JLM original params
  int master = 1; // work is done on the master node
  printf("\nRelaxing reference Cell\n\n");


  //stop if already found reference frequencies
  if(Params::Parameters().NumberOfCompletedQuasiSteps() > 0){
    //reference value
    Coords_ref = Cluster::cluster().GetSymmetryUniqueCoordinates();
    a_ref = Coords_ref[3*UniqueAtoms];
    b_ref = Coords_ref[3*UniqueAtoms+1];
    c_ref = Coords_ref[3*UniqueAtoms+2];
    alpha_ref = Coords_ref[3*UniqueAtoms+3];
    beta_ref = Coords_ref[3*UniqueAtoms+4];
    gamma_ref = Coords_ref[3*UniqueAtoms+5];
    volume_ref = Cluster::cluster().UnitCellVolume();
    //set the volume to a constant
    Params::Parameters().SetConstantVolume(volume_ref);
    
    printf("\nSkipping reference. Already availible.\n");
    
    //reference frequencies already avalible
    if(!Params::Parameters().IsSupercellJob()) 
      ReadMoldenFrequencies();
  
  }else{

    if(Params::Parameters().RunJobs())
      dl_find_(&nvarin, &nvarin2, &nspec, &master);
    else if(Params::Parameters().GetMMType() != 5){//Forces not needed for crystal09 frequency calculation
     Params::Parameters().SetParameter("JOBTYPE","FORCE");        
     Cluster::cluster().RunJobsAndComputeEnergy();    
    }

    

    //OptimizeGeometry("SteepestDescent");
  
    //print optimized reference structure
    FILE *reference;
    string reference_file = "reference.in";
    if ((reference = fopen(reference_file.c_str(),"w"))==NULL) {
      printf("Error: Cannot open file '%s'\n",reference_file.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintInputFile(reference);
    printf("\nReference file written to '%s'\n",reference_file.c_str());
    fclose(reference);
    
    //reference value
    Coords_ref = Cluster::cluster().GetSymmetryUniqueCoordinates();
    a_ref = Coords_ref[3*UniqueAtoms];
    b_ref = Coords_ref[3*UniqueAtoms+1];
    c_ref = Coords_ref[3*UniqueAtoms+2];
    alpha_ref = Coords_ref[3*UniqueAtoms+3];
    beta_ref = Coords_ref[3*UniqueAtoms+4];
    gamma_ref = Coords_ref[3*UniqueAtoms+5];
    volume_ref = Cluster::cluster().UnitCellVolume();
    //set the volume to a constant
    Params::Parameters().SetConstantVolume(volume_ref);
    
    //Get the frequencies of the relaxed reference cell
    
    //get frequencies of the relaxed reference cell.
    printf("\nFind frequencies\n\n");
    
    Params::Parameters().SetParameter("JOBTYPE","HESSIAN");
    
    //changing the QC rem to perform frequency calculations
    qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    
    Cluster::cluster().RunJobsAndComputeEnergy();
    //Cluster::cluster().GetFrequency().Print("freqs");
    
    char vol_str[75];
    
    //only gamma point frequencies
    if(!Params::Parameters().IsSupercellJob()){
      
      //Setting the reference frequencies
      //Cluster::cluster().SetReferenceFrequencies(Cluster::cluster().GetFrequency());
      freq_ref = Cluster::cluster().GetFrequency();
      
      //Cluster::cluster().SetReferenceVolume(volume_ref);
      
      //print frequencies
      printf("volume_ref %f\n",volume_ref);
      sprintf(vol_str,"ReferenceFrequencies %.3f ",volume_ref);
      
      Cluster::cluster().PrintVibrationalFrequencies(freq_ref,vol_str,1);
      //Cluster::cluster().PrintVibrationalFrequencies(Cluster::cluster().GetFrequency(),volume_ref,1);
      
      //only print volume
      if(Params::Parameters().UseVolumeGruneisen()){
	sprintf(vol_str,"ReferenceFrequencies %15.9f ",volume_ref);   
      }//print lattice parameters
      else{
	sprintf(vol_str,"ReferenceFrequencies %15.9f %15.9f %15.9f",  a_ref,b_ref,c_ref);
      }
      
      Cluster::cluster().PrintVibrationalFrequencies(freq_ref,vol_str,1);
      
    }
    
    //use frequencies beyond gamma point
    else{
      
      //Updating supercell for new unit cell coordinates and parameters
      Supercell::supercell().UpdateSupercellCoordinates();
      
      if (!Params::Parameters().NeglectManyBody()) {
	printf("doing Supercell tinker job\n");
	//Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	Supercell::supercell().RunHMBISupercellTinkerHessianJob();
      }
      //Supercell::supercell().SetSupercellHessian(); //JLM Redundant calls
      Supercell::supercell().ComputeHMBISupercellHessian();
      Supercell::supercell().ComputePhonons();
      
      //set reference volumes
      //Cluster::cluster().SetReferenceVolume(volume_ref);
      //print frequencies
      printf("volume_ref %f\n",volume_ref);
      
      //only print volume
      if(Params::Parameters().UseVolumeGruneisen())
	sprintf(vol_str,"ReferenceFrequencies %15.9f ",volume_ref);
      //print lattice parameters
      else
	sprintf(vol_str,"ReferenceFrequencies %15.9f %15.9f %15.9f",  a_ref,b_ref,c_ref);
      
      //supercell phonons
      Matrix phonons_mat = Supercell::supercell().GetPhononMatrix();
      //Vector phonons_vec = phonons_mat.StackColumnsOfMatrix();
      freq_ref = phonons_mat.StackColumnsOfMatrix();
      //phonons_vec.Print("phonons_vec");
      
      //Cluster::cluster().SetReferenceFrequencies(phonons_vec);
      //Cluster::cluster().PrintVibrationalFrequencies(phonons_vec,vol_str,1);
      Cluster::cluster().PrintVibrationalFrequencies(freq_ref,vol_str,1);
      
      //counting number of imaginary frequencies
      for(int i = 0 ; i < 3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point; i++)
	if(freq_ref[i] < 0.0)
	  imag_num++;
      
    }
  }

}

//Quasi-harmonic approximation expanding and contracting isotropically
void Quasiharmonic::DetermineIsotropicFrequencies(ifstream &infile,double &VolumePlus,double &VolumeMinus,
						  Vector &FreqPlus,Vector &FreqMinus,int &imag_num){
  printf("\ndoing Quasi-harmonic approximation\n");
  
  //set reference values
  DetermineReferenceFrequencies(infile,imag_num);

  //Frozen lattice params for all opimizations while approximating the dependance of volume on the vibration.
  //Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
  

    
  //values used for optimizer
  int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  int nvarin, nvarin2, nspec;
  if(Params::Parameters().UseFullQMOnly()){
   nvarin =  3*UniqueAtoms+9;
   nvarin2 = 5*UniqueAtoms+15;
   nspec = 3*UniqueAtoms+9; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }
  else{
   nvarin =  3*UniqueAtoms+6;
   nvarin2 = 5*UniqueAtoms+10;
   nspec = 3*UniqueAtoms+6; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }

  /*int nvarin =  3*UniqueAtoms+6;
  int nvarin2 = 5*UniqueAtoms+10;// nframe*nat*3 + nweight + nmass + n_po_scaling
  int nspec = 2*UniqueAtoms+4;// nspec= nat + nz + 5*ncons + 2*nconn*/ //JLM original params
  int master = 1; // work is done on the master node
  
  //varying the volume isotropically
  double delta = 10.0;
  //increase volume by delta


  printf("Next geometry\n");
  
  //name of file Plus input file is stored and name of the directory (within the
  //respective QM/MM directory) where the terms are stored 
  hess_type = "PlusFreq";
  
  //change job type to Optimization and reread the qc_rem
  Params::Parameters().SetParameter("JOBTYPE","Opt");
  string qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
  Params::Parameters().SetParameter("QC_REM",qc_rem);
  
  Vector Coords = Coords_ref;
  
  //double volume = volume_ref + delta;
  //Coords.Print("Coords before");
  //Coords[3*UniqueAtoms] = a_ref*pow(volume/volume_ref,1.0/3.0); 
  //Coords[3*UniqueAtoms+1] = b_ref*pow(volume/volume_ref,1.0/3.0);
  //Coords[3*UniqueAtoms+2] = c_ref*pow(volume/volume_ref,1.0/3.0);
  double volume;

  //read coordinates from file for Plus geometry instead of increasing volume isotropically
  if(Params::Parameters().QuasiReadQuasiGeometries()){
    Coords = ReadGeometryFromInput(volume);
  }
  else{  //expand the cell isotropically
  
    volume = volume_ref + delta;
    Coords.Print("Coords before");
    Coords[3*UniqueAtoms] = a_ref*pow(volume/volume_ref,1.0/3.0); 
    Coords[3*UniqueAtoms+1] = b_ref*pow(volume/volume_ref,1.0/3.0);
    Coords[3*UniqueAtoms+2] = c_ref*pow(volume/volume_ref,1.0/3.0);


  }

  cout << "Plus structure's volume is " << volume << "\n";

  //Coords.Print("Coords after");
  //unfreezing lattice parameters so they can be altered
  Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","0");
  Cluster::cluster().MaintainCartesianSymmetry(Coords,true);
  Cluster::cluster().UpdateLatticeParams(Coords);
  Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
  Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
  Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
 
  //set the volume to a constant
  //Params::Parameters().SetConstantVolume(volume);
  
  //freezing lattice params
  Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
  
  //Coords.PrintGradient("full Coords");
  //change job type to Optimization and reread the qc_rem
  Params::Parameters().SetParameter("JOBTYPE","Opt");
  qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
  Params::Parameters().SetParameter("QC_REM",qc_rem);
  Cluster::cluster().SetNewCoordinates(Coords);
  printf("\nRelaxing Plus geometry\n\n");

  if(Params::Parameters().NumberOfCompletedQuasiSteps() > 1){
    
    printf("\nSkipping Plus geometry\n");
  }
  else{
    if(Params::Parameters().RunJobs()){
      dl_find_(&nvarin, &nvarin2, &nspec, &master);

      //print optimized structure
      FILE *file;
      string file_name = hess_type + ".in";
      if ((file = fopen(file_name.c_str(),"w"))==NULL) {
	printf("Error: Cannot open file '%s'\n",file_name.c_str());
	exit(1);
      } 
      Cluster::cluster().PrintInputFile(file);
      printf("\nfile written to '%s'\n",file_name.c_str());
      fclose(file);
      
    }else if(Params::Parameters().GetMMType() != 5){//Forces not needed for crystal09 frequency calculation
     Params::Parameters().SetParameter("JOBTYPE","FORCE");        
     Cluster::cluster().RunJobsAndComputeEnergy();   
    }
      
    //OptimizeGeometry("SteepestDescent");
  
    //get frequencies of the relaxed reference cell.
    printf("\nFind frequencies\n\n");
    Params::Parameters().SetParameter("JOBTYPE","HESSIAN");
 
    //changing the QC rem to perform frequency calculations
    qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().RunJobsAndComputeEnergy();
    //Cluster::cluster().GetFrequency().Print("freqs");
    
    //volume string for printing to frequency file
    char vol_str[30];
    
    //only gamma point frequencies
    if(!Params::Parameters().IsSupercellJob()){
      //Setting the Plus values
      FreqPlus = Cluster::cluster().GetFrequency();
      VolumePlus = volume;
      
      //print frequencies
      sprintf(vol_str,"PlusVolume %.3f ",volume);
      printf("%s\n",vol_str);
      Cluster::cluster().PrintVibrationalFrequencies(Cluster::cluster().GetFrequency(),vol_str,0);
    }
    //use frequencies beyond gamma point
    else{
      
      //Updating supercell for new unit cell coordinates and parameters
      Supercell::supercell().UpdateSupercellCoordinates();
      
      //use frequencies beyond gamma point
      if (!Params::Parameters().NeglectManyBody()) {
	printf("doing Supercell tinker job\n");
	//Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	Supercell::supercell().RunHMBISupercellTinkerHessianJob();
      }
      //Supercell::supercell().SetSupercellHessian(); //JLM Redundant calls
      Supercell::supercell().ComputeHMBISupercellHessian();
      Supercell::supercell().ComputePhonons();
      
      //Setting the Plus values
      Matrix phonons_mat = Supercell::supercell().GetPhononMatrix();
      FreqPlus = phonons_mat.StackColumnsOfMatrix();
      VolumePlus = volume;
     
      sprintf(vol_str,"PlusVolume %.3f ",volume);
      printf("%s\n",vol_str);
      Cluster::cluster().PrintVibrationalFrequencies(FreqPlus,vol_str,0);
      
      //counting number of imaginary frequencies
      int imag_plus=0;
      for(int i = 0 ; i < 3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point; i++){
	if(FreqPlus[i] < 0.0)
	  imag_plus++;
      }
      if(imag_plus != imag_num && imag_plus > 3){
	printf("Error::Quasiharmonic::DetermineIsotropicFrequencies Number of imaginary frequency changes when cell volumes is increased\n");
	printf("imag_plus = %i\n",imag_plus);
	//exit(0); //JLM
      }
      
    }
  }

  //name of file Minus input file is stored and name of the directory (within the
  //respective QM/MM directory) where the terms are stored 
  hess_type = "MinusFreq";
  
  //read coordinates from file for this geometry instead of increasing volume it isotropically
  if(Params::Parameters().QuasiReadQuasiGeometries()){
    Coords = ReadGeometryFromInput(volume);
  }else{  //decrease volume by delta
  
    Coords = Coords_ref;
    volume = volume_ref - delta;
    //Coords.Print("Coords before");
    Coords[3*UniqueAtoms] = a_ref*pow(volume/volume_ref,1.0/3.0); 
    Coords[3*UniqueAtoms+1] = b_ref*pow(volume/volume_ref,1.0/3.0);
    Coords[3*UniqueAtoms+2] = c_ref*pow(volume/volume_ref,1.0/3.0);
  }

  cout << "Minus Structure's volume is " << volume << "\n";
 
  //Coords.Print("Coords after");
  //unfreezing lattice parameters so they can be altered
  Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","0");
  Cluster::cluster().MaintainCartesianSymmetry(Coords,true);
  Cluster::cluster().UpdateLatticeParams(Coords);
  Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
  Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
  Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
  //Coords.PrintGradient("full Coords");   
  Cluster::cluster().SetNewCoordinates(Coords);
  
  //freezing lattice params
  Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
  
  printf("\nRelaxing Minus geometry\n\n");
  
  //change job type to Optimization and reread the qc_rem
  Params::Parameters().SetParameter("JOBTYPE","Opt");
  qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
  Params::Parameters().SetParameter("QC_REM",qc_rem);
  Cluster::cluster().SetNewCoordinates(Coords);


  if(Params::Parameters().NumberOfCompletedQuasiSteps() > 3){
    printf("\nSkipping Minus geometry\n");
  }
  else{
    if(Params::Parameters().RunJobs()){
      dl_find_(&nvarin, &nvarin2, &nspec, &master);

      //print optimized structure
      FILE *file;
    string file_name = hess_type + ".in";
    if ((file = fopen(file_name.c_str(),"w"))==NULL) {
      printf("Error: Cannot open file '%s'\n",file_name.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintInputFile(file);
    printf("\nfile written to '%s'\n",file_name.c_str());
    fclose(file);

    }
    else if(Params::Parameters().GetMMType() != 5){//Forces not needed for crystal09 frequency calculation
     Params::Parameters().SetParameter("JOBTYPE","FORCE");        
     Cluster::cluster().RunJobsAndComputeEnergy();    
    }
    //OptimizeGeometry("SteepestDescent");
    
    //Coords.PrintGradient("full Coords");
  
    //get frequencies of the relaxed reference cell.
    printf("\nFind frequencies\n\n");
    Params::Parameters().SetParameter("JOBTYPE","HESSIAN");
    
    //changing the QC rem to perform frequency calculations
    qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().RunJobsAndComputeEnergy();
    //Cluster::cluster().GetFrequency().Print("freqs");
    
    //print frequencies
    char vol_str[30];

    //only gamma point frequencies
    if(!Params::Parameters().IsSupercellJob()){  
      
      //print frequencies
      sprintf(vol_str,"MinusVolume %.3f ",volume);
      printf("%s\n",vol_str);
      Cluster::cluster().PrintVibrationalFrequencies(Cluster::cluster().GetFrequency(),vol_str,0);
      
      //Setting the Minus values
      FreqMinus = Cluster::cluster().GetFrequency();
      VolumeMinus = volume; 
      
    }
    //use frequencies beyond gamma point
    else{
    
      //Updating supercell for new unit cell coordinates and parameters
      Supercell::supercell().UpdateSupercellCoordinates();
      
      if (!Params::Parameters().NeglectManyBody()) {
	printf("doing Supercell tinker job\n");
	//Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	Supercell::supercell().RunHMBISupercellTinkerHessianJob();
      }
      //Supercell::supercell().SetSupercellHessian();  //JLM Redundant calls 
      Supercell::supercell().ComputeHMBISupercellHessian();
      Supercell::supercell().ComputePhonons();
      
      //Setting the Plus values
      Matrix phonons_mat = Supercell::supercell().GetPhononMatrix();
      FreqMinus = phonons_mat.StackColumnsOfMatrix();
      VolumePlus = volume;
      
      //print frequencies
      sprintf(vol_str,"MinusVolume %.3f ",volume);
      printf("%s\n",vol_str);
      Cluster::cluster().PrintVibrationalFrequencies(FreqMinus,vol_str,0);
      
      //counting number of imaginary frequencies
      int imag_minus=0;
      for(int i = 0 ; i < 3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point; i++){
	if(FreqMinus[i] < 0.0)
	  imag_minus++;
      }
      if(imag_minus != imag_num && imag_minus > 2){
	printf("Error::Number of imaginary frequency changes when cell volumes is decreased\n");
	printf("imag_minus = %i\n",imag_minus);
	//exit(0); //JLM
      }
    }
  }

}

//Used to find grunesian parameter for single lattice params
//Is currently broken. Do not use
void Quasiharmonic::DetermineLatticeFrequencies(ifstream &infile,double &LatticePlus,double &LatticeMinus,
						Vector &FreqPlus,Vector &FreqMinus,int &imag_num,int param){

  printf("\ndoing Quasi-harmonic approximation\n");

    //set reference values
  if(param == 0)  
    Quasiharmonic::quasiharmonic().DetermineReferenceFrequencies(infile,imag_num);


  //values used for optimizer
  int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  int nvarin, nvarin2, nspec;
  if(Params::Parameters().UseFullQMOnly()){
   nvarin =  3*UniqueAtoms+9;
   nvarin2 = 5*UniqueAtoms+15;
   nspec = 3*UniqueAtoms+9; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }
  else{
   nvarin =  3*UniqueAtoms+6;
   nvarin2 = 5*UniqueAtoms+10;
   nspec = 3*UniqueAtoms+6; //nspec should be: nat + nz + 5*ncons + 2*nconn + nat
  }

  /*int nvarin =  3*UniqueAtoms+6;
  int nvarin2 = 5*UniqueAtoms+10;// nframe*nat*3 + nweight + nmass + n_po_scaling
  int nspec = 2*UniqueAtoms+4;// nspec= nat + nz + 5*ncons + 2*nconn*/ //JLM original params
  int master = 1; // work is done on the master node

  double delta;
  //double delta_volume = 15.0;
  //double volume = volume_ref+delta_volume;
  Vector Coords = Coords_ref;


 //header for printing to frequency file
  char param_str[30];

  //frequency already found at this geometry if true
  bool SkipPlus = false;

  //get plus lattice length a
  if(param==0){
    delta = 0.1;
    Coords[3*UniqueAtoms] += delta; 
    //Coords[3*UniqueAtoms] = a_ref*(volume/volume_ref);
    LatticePlus = Coords[3*UniqueAtoms];

    //header for frequency file
    sprintf(param_str,"PlusA %15.9f ",LatticePlus);

    //name of file where frequencies and modes are printed 
    hess_type = "PlusFreqA";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 1)
      SkipPlus = true;
    else
      printf("\na altered to %f\n",LatticePlus);

  }
  //get plus for for lattice length b
  else if(param==1){
    delta = 0.1;
    Coords[3*UniqueAtoms+1] += delta; 
    //Coords[3*UniqueAtoms+1] = b_ref*(volume/volume_ref);

    LatticePlus = Coords[3*UniqueAtoms+1];


    //header for frequency file
    sprintf(param_str,"PlusB %15.9f ",LatticePlus);

    //name of file where frequencies and modes are printed 
    hess_type = "PlusFreqB";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 3)
      SkipPlus = true;
    else
      printf("\nb altered to %f\n",LatticePlus);   


  }
  //get plus for lattice length c
  else if(param==2){
    delta = 0.1;
    Coords[3*UniqueAtoms+2] += delta; 
    //Coords[3*UniqueAtoms+2] = c_ref*(volume/volume_ref);
    LatticePlus = Coords[3*UniqueAtoms+2];


    //header for frequency file
    sprintf(param_str,"PlusC %15.9f ",LatticePlus);

    //name of file where frequencies and modes are printed 
    hess_type = "PlusFreqC";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 5)
      SkipPlus = true;
    else
      printf("\nc altered to %f\n",LatticePlus);


  }//get plus for lattice angles
  else if (param>2 && param<6){
    printf("ERROR::Quasiharmonic::DetermineLatticeFrequencies()::Currently unable to determine grunesian parameters. for lattice angles\n");
    printf("param = %i\n",param);
    exit(0);
  }
  else{
    printf("ERROR::Quasiharmonic::DetermineLatticeFrequencies()::Unable to recognize lattice parameter.\n param = %i\n",param);
    exit(0);
  }

  //skip plus geometry because it is already read.
  if(SkipPlus){
    printf("\nSkip geometry %s\n",hess_type.c_str());
  }
  else{
    //unfreezing lattice parameters so they can be altered
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","0");
    Cluster::cluster().MaintainCartesianSymmetry(Coords,true);
    Cluster::cluster().UpdateLatticeParams(Coords);
    Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
    Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
    Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
    
    //freezing lattice params
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
    
    //Coords.PrintGradient("full Coords");
    //change job type to Optimization and reread the qc_rem
    Params::Parameters().SetParameter("JOBTYPE","Opt");
    string qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().SetNewCoordinates(Coords);
    printf("\nRelaxing Plus geometry\n\n");
    
    if(Params::Parameters().RunJobs())
      dl_find_(&nvarin, &nvarin2, &nspec, &master);
    else if(Params::Parameters().GetMMType() != 5){//Forces not needed for crystal09 frequency calculation
     Params::Parameters().SetParameter("JOBTYPE","FORCE");        
     Cluster::cluster().RunJobsAndComputeEnergy();    
    }
    //OptimizeGeometry("SteepestDescent");
    
    //print optimized structure
    FILE *file;
    string file_name = hess_type + ".in";
    if ((file = fopen(file_name.c_str(),"w"))==NULL) {
      printf("Error: Cannot open file '%s'\n",file_name.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintInputFile(file);
    printf("\nfile written to '%s'\n",file_name.c_str());
    fclose(file);
    
    //get frequencies of the relaxed reference cell.
    printf("\nFind frequencies\n\n");
    Params::Parameters().SetParameter("JOBTYPE","HESSIAN");
    
    //changing the QC rem to perform frequency calculations
    qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().RunJobsAndComputeEnergy();
    //Cluster::cluster().GetFrequency().Print("freqs");
    
    
    //only gamma point frequencies
    if(!Params::Parameters().IsSupercellJob()){
      
      //print frequencies
      FreqPlus = Cluster::cluster().GetFrequency();
      printf("%s\n",param_str);
      Cluster::cluster().PrintVibrationalFrequencies(Cluster::cluster().GetFrequency(),param_str,0);
    }
    //use frequencies beyond gamma point
    else{
      
      //Updating supercell for new unit cell coordinates and parameters
      Supercell::supercell().UpdateSupercellCoordinates();
      
      //use frequencies beyond gamma point
      if (!Params::Parameters().NeglectManyBody()) {
	printf("doing Supercell tinker job\n");
	//Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	Supercell::supercell().RunHMBISupercellTinkerHessianJob();
      }
      //Supercell::supercell().SetSupercellHessian();  //JLM Redundant calls
      Supercell::supercell().ComputeHMBISupercellHessian();
      Supercell::supercell().ComputePhonons();
      
      
      //Setting the Plus values
      Matrix phonons_mat = Supercell::supercell().GetPhononMatrix();
      FreqPlus = phonons_mat.StackColumnsOfMatrix();
      
      
      printf("%s\n",param_str);
      Cluster::cluster().PrintVibrationalFrequencies(FreqPlus,param_str,0);
      
      //counting number of imaginary frequencies
      int imag_plus=0;
      for(int i = 0 ; i < 3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point; i++){
	if(FreqPlus[i] > 0.0)
	  imag_plus++;
      }
      if(imag_plus != imag_num){
	printf("Error::Number of imaginary frequency changes when cell volumes is increased\n");
	exit(0);
      }
      
    }
  }

  //returning coords to reference values
  Coords = Coords_ref;
  //volume = volume_ref-delta_volume;
  
  //frequency already found at this geometry if true
  bool SkipMinus = false; 

  //get minus lattice length a
  if(param==0){
    
    Coords[3*UniqueAtoms] -= delta; 
    //Coords[3*UniqueAtoms] = a_ref*(volume/volume_ref);
    LatticeMinus = Coords[3*UniqueAtoms];



    //header for frequency file
    sprintf(param_str,"MinusA %15.9f ",LatticeMinus);

    //name of file where frequencies and modes are printed 
    hess_type = "MinusFreqA";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 2)
      SkipMinus = true;
    else
      printf("\na altered to %f\n",LatticeMinus);
  }
  //get minus for for lattice length b
  else if(param==1){

    Coords[3*UniqueAtoms+1] -= delta; 
    //Coords[3*UniqueAtoms+1] = a_ref*(volume/volume_ref);
    LatticeMinus = Coords[3*UniqueAtoms+1];

    //header for frequency file
    sprintf(param_str,"MinusB %15.9f ",LatticeMinus);

    //name of file where frequencies and modes are printed 
    hess_type = "MinusFreqB";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 4)
      SkipMinus = true;
    else
      printf("\nb altered to %f\n",LatticeMinus); 
  }
  //get minus for lattice length c
  else if(param==2){

    Coords[3*UniqueAtoms+2] -= delta; 
    //Coords[3*UniqueAtoms] = c_ref*(volume/volume_ref);
    LatticeMinus = Coords[3*UniqueAtoms+2];

    //header for frequency file
    sprintf(param_str,"MinusC %15.9f ",LatticeMinus);

    //name of file where frequencies and modes are printed 
    hess_type = "MinusFreqC";

    if(Params::Parameters().NumberOfCompletedQuasiSteps() > 6){
      SkipMinus = true;
    }
    else
      printf("\nc altered to %f\n",LatticeMinus);

  }//get minus for lattice angles
  else if (param>2 && param<6){
    printf("ERROR::Quasiharmonic::DetermineLatticeFrequencies()::Currently unable to determine grunesian parameters. for lattice angles\n");
    printf("param = %i\n",param);
    exit(0);
  }  
  else{
    printf("ERROR::Quasiharmonic::DetermineLatticeFrequencies()::Unable to reconize lattice parameter.\n param = %i\n",param);
    exit(0);
  }

  //skip minus geometry because it is already read.
  if(SkipMinus){
    printf("\nSkip geometry %s\n",hess_type.c_str());
  }
  else{
    //unfreezing lattice parameters so they can be altered
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","0");
    Cluster::cluster().MaintainCartesianSymmetry(Coords,true);
    Cluster::cluster().UpdateLatticeParams(Coords);
    Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
    Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
    Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
    
    //freezing lattice params
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","1");
    
    //Coords.PrintGradient("full Coords");
    //change job type to Optimization and reread teh qc_rem
    Params::Parameters().SetParameter("JOBTYPE","Opt");
    string qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().SetNewCoordinates(Coords);
    printf("\nRelaxing Minus geometry\n\n");
    
    if(Params::Parameters().RunJobs())
      dl_find_(&nvarin, &nvarin2, &nspec, &master);
    else if(Params::Parameters().GetMMType() != 5){//Forces not needed for crystal09 frequency calculation
     Params::Parameters().SetParameter("JOBTYPE","FORCE");        
     Cluster::cluster().RunJobsAndComputeEnergy();    
    }
    //OptimizeGeometry("SteepestDescent");
    
    //print optimized structure
    FILE *file;
    string file_name = hess_type + ".in";
    if ((file = fopen(file_name.c_str(),"w"))==NULL) {
      printf("Error: Cannot open file '%s'\n",file_name.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintInputFile(file);
    printf("\nfile written to '%s'\n",file_name.c_str());
    fclose(file);
    
    //get frequencies of the relaxed reference cell.
    printf("\nFind frequencies\n\n");
    Params::Parameters().SetParameter("JOBTYPE","HESSIAN");
    
    //changing the QC rem to perform frequency calculations
    qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
    Params::Parameters().SetParameter("QC_REM",qc_rem);
    Cluster::cluster().RunJobsAndComputeEnergy();
    //Cluster::cluster().GetFrequency().Print("freqs");
    
    
    //only gamma point frequencies
    if(!Params::Parameters().IsSupercellJob()){
      
      //print frequencies
      FreqMinus = Cluster::cluster().GetFrequency();
      printf("%s\n",param_str);
      Cluster::cluster().PrintVibrationalFrequencies(Cluster::cluster().GetFrequency(),param_str,0);
    }
    //use frequencies beyond gamma point
    else{
      
      //Updating supercell for new unit cell coordinates and parameters
      Supercell::supercell().UpdateSupercellCoordinates();
      
      //use frequencies beyond gamma point
      if (!Params::Parameters().NeglectManyBody()) {
	printf("doing Supercell tinker job\n");
	//Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	Supercell::supercell().RunHMBISupercellTinkerHessianJob();
      }
      //Supercell::supercell().SetSupercellHessian();  //JLM Redundant calls
      Supercell::supercell().ComputeHMBISupercellHessian();
      Supercell::supercell().ComputePhonons();
    
      
      //Setting the Plus values
      Matrix phonons_mat = Supercell::supercell().GetPhononMatrix();
      FreqMinus = phonons_mat.StackColumnsOfMatrix();
      
      
      printf("%s\n",param_str);
      Cluster::cluster().PrintVibrationalFrequencies(FreqMinus ,param_str,0);
      
      //counting number of imaginary frequencies
      int imag_minus=0;
      for(int i = 0 ; i < 3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point; i++){
	if(FreqPlus[i] > 0.0)
	  imag_minus++;
      }
      if(imag_minus != imag_num){
	printf("Error::Number of imaginary frequency changes when cell volumes is increased\n");
	exit(0);
      }
      
    }
  }
}

/*
void Quasiharmonic::ReadReferenceFrequencies(infile infile){

  int total_k_points = Params::Parameters().GetSuperCellTotalKPoint();
  int Natoms = Clusters::cluster().GetTotalNumberOfAtoms();
  if(!Params::Parameters().IsSupercellJob() &&
     ref_freq.GetLength() != 3*Natoms){
    printf("Error : Quasiharmonic::ReadReferenceFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   ref_freq.GetLength(),3*Natoms);
    exit(0);
  }else if(Params::Parameters().IsSupercellJob() &&
	   ref_freq.GetLength() != 3*Natoms*total_k_points){
    printf("Error : Quasiharmonic::ReadReferenceFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   ref_freq.GetLength(),3*Natoms);
    exit(0);
  }

  //Get name frequencie file
  FILE *freq_file;
  string filename = Cluster::cluster().GetInputFilename();
  string freq_name = filename.substr(0,filename.size()-3);
  freq_name += ".freq";
  
  //Open the frequency file 
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("Quasiharmonic::ReadIsotropicFrequencies: Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  }  
  

  for(int i=0;i<freq_ref.GetLength();i++){

  //read reference frequencies
  string line;
  getline(infile,line);

  istringstream iss(line);
  string tmp;//discarding the "ReferenceFrequencies" string
  iss >> tmp;
  iss >> volume_ref;
  //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){

    getline(infile,line);
    istringstream iss(line);
    iss >> freq_ref[i];
  }


}
*/


//read frequencies from file
void Quasiharmonic::ReadIsotropicFrequencies(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus){

  //Making sure that ref_freq is the right length
  int total_k_points = Params::Parameters().GetSuperCellTotalKPoint();
  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
  if(!Params::Parameters().IsSupercellJob() &&
     freq_ref.GetLength() != 3*Natoms){
    printf("Error::Quasiharmonic::ReadIsotropicFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   freq_ref.GetLength(),3*Natoms);
    exit(0);
  }else if(Params::Parameters().IsSupercellJob() &&
	   freq_ref.GetLength() != 3*Natoms*total_k_points){
    printf("Error::Quasiharmonic::ReadIsotropicFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   freq_ref.GetLength(),3*Natoms*total_k_points);
    exit(0);
  }

  //Check to make sure all the sets of freqiencies vectors are the same length.
  //This does not check if the frequency in the frequency files of the right number of entries
  if(FreqPlus.GetLength() != freq_ref.GetLength() ){
    printf("ERROR::Quasiharmonic::ReadIsotropicFrequencies : FreqPlus vector incorrect length \n  Length = %i instend of %i\n",
	   FreqPlus.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(FreqMinus.GetLength() != freq_ref.GetLength()){
    printf("ERROR::Quasiharmonic::ReadIsotropicFrequencies : FreqMinus vector incorrect length \n  Length = %i instend of %i \n",
	   FreqMinus.GetLength(),freq_ref.GetLength());
    exit(1);
  }


  //Get name of the frequency file
  FILE *freq_file;
  string filename = Cluster::cluster().GetInputFilename();
  string freq_name = filename.substr(0,filename.size()-3);
  freq_name += ".freq";
  
  //counts the number of geometry. Should be two (excluding reference frequencies)
  int countgeom = 0;

  //Open the frequency file 
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("ERROR:: Quasiharmonic::ReadIsotropicFrequencies: Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  }  

  string line;
  while ( !infile.eof() ) {
    //read reference frequencies
    getline(infile,line);

    //reference frequencies
    if(line.substr(0,20) == "ReferenceFrequencies"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discarding the "ReferenceFrequencies" string
      iss >> tmp;
      iss >> volume_ref;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<freq_ref.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> freq_ref[i];
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
    printf("ERROR::Quasihamonic::ReadIsotropicFrequencies() Incorrect number of volumes countgeom = %i\n",
	   countgeom);
    exit(1);
  }


  //for super cells make sure that all volumes have the same number of real freqencies
  // or for gamma point frequencies, no more than 3 non-positive frequencies (translational frequencies)
  int imag_num=0;
  int imag_plus=0;
  int imag_minus=0;
  for(int i = 0;i < freq_ref.GetLength(); i++)
    if(freq_ref[i] < 0.0)
       imag_num++;
  
  for(int i = 0; i < FreqPlus.GetLength(); i++)
    if(FreqPlus[i] < 0.0)
      imag_plus++;
   
  if(imag_num != imag_plus && (freq_ref.Nonpositive() > 3 && FreqPlus.Nonpositive() > 3) ){
   printf("ERROR::Quasihamonic::ReadIsotropicFrequencies() Number of real frequencies differ between reference geometry and geometry of increased volume.\n");
   exit(0);
  }

  for(int i = 0; i < FreqMinus.GetLength(); i++)
    if(FreqMinus[i] < 0.0)
      imag_minus++;

  if(imag_num != imag_minus && (freq_ref.GetLength()-freq_ref.Nonpositive() > 3 &&  FreqMinus.Nonpositive() > 3)){
   printf("Quasihamonic::ReadIsotropicFrequencies::Error::Number of real frequencies differ between reference geometry and geometry of decreased volume.\n");
   exit(0);
  }

}

//currently not accounting for lattice angles
//This function is broken. Do not use
void Quasiharmonic::ReadAnisotropicFrequencies(Vector& plus_a_freq,double& plus_a,
					       Vector& minus_a_freq,double& minus_a,
					       Vector& plus_b_freq,double& plus_b,
					       Vector& minus_b_freq,double& minus_b,
					       Vector& plus_c_freq,double& plus_c,
					       Vector& minus_c_freq,double& minus_c){
 
  //Making sure that ref_freq is the right length
  int total_k_points = Params::Parameters().GetSuperCellTotalKPoint();
  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
  if(!Params::Parameters().IsSupercellJob() &&
     freq_ref.GetLength() != 3*Natoms){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   freq_ref.GetLength(),3*Natoms);
    exit(0);
  }else if(Params::Parameters().IsSupercellJob() &&
	   freq_ref.GetLength() != 3*Natoms*total_k_points){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies():: ref_freq has %i entries instead of %i entries\n",
	   freq_ref.GetLength(),3*Natoms);
    exit(0);
  }
  
  //Check to make sure all the sets of freqiencies vectors are the same length.
  //This does not check if the frequency in the frequency files of the right number of entries
  if(plus_a_freq.GetLength() != freq_ref.GetLength() ){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : plus_a_freq vector incorrect length \n  Length = %i instend of %i\n",
	   plus_a_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(minus_a_freq.GetLength() != freq_ref.GetLength()){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : minus_a_freq vector incorrect length \n  Length = %i instend of %i \n",
	   minus_a_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(plus_b_freq.GetLength() != freq_ref.GetLength() ){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : plus_b_freq vector incorrect length \n  Length = %i instend of %i\n",
	   plus_b_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(minus_b_freq.GetLength() != freq_ref.GetLength()){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : minus_b_freq vector incorrect length \n  Length = %i instend of %i \n",
	   minus_b_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(plus_c_freq.GetLength() != freq_ref.GetLength() ){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : plus_c_freq vector incorrect length \n  Length = %i instend of %i\n",
	   plus_c_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(minus_c_freq.GetLength() != freq_ref.GetLength()){
    printf("Error::Quasiharmonic::ReadAnisotropicFrequencies : minus_c_freq vector incorrect length \n  Length = %i instend of %i \n",
	   minus_c_freq.GetLength(),freq_ref.GetLength());
    exit(1);
  }

  //Get name of frequency file
  FILE *freq_file;
  string filename = Cluster::cluster().GetInputFilename();
  string freq_name = filename.substr(0,filename.size()-3);
  freq_name += ".freq";
  
  //counts the number of geometries. Should be two (excluding reference frequencies)
  int countgeom = 0;
  
  //Open the frequency file 
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("Quasiharmonic::ReadAnisotropicFrequencies: Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  }  
  
  string line;
  while ( !infile.eof() ) {
    //read reference frequencies
    getline(infile,line);
    
    //reference frequencies
    if(line.substr(0,20) == "ReferenceFrequencies"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discarding the "ReferenceFrequencies" string
      iss >> tmp;
      iss >> a_ref;
      iss >> b_ref;
      iss >> c_ref;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<freq_ref.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> freq_ref[i];
      }
    }
    
    //plus a frequencies
    if(line.substr(0,5) == "PlusA"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> plus_a;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<plus_a_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> plus_a_freq[i];
      }
    }
    //minus a frequencies
    if(line.substr(0,6) == "MinusA"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> minus_a;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<minus_a_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> minus_a_freq[i];
      }
    }

    //plus b frequencies
    if(line.substr(0,5) == "PlusB"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> plus_b;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<plus_b_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> plus_b_freq[i];
      }
    }
    //minus b frequencies
    if(line.substr(0,6) == "MinusB"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> minus_b;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<plus_b_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> minus_b_freq[i];
      }
    }

    //plus c frequencies
    if(line.substr(0,5) == "PlusC"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> plus_c;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<plus_c_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> plus_c_freq[i];
      }
    }
    //minus b frequencies
    if(line.substr(0,6) == "MinusC"){
      countgeom++;
      istringstream iss(line);
      string tmp;//discard the "PlusVolume" string
      iss >> tmp;
      iss >> minus_c;
      //for(int i=0;i<3*GetTotalNumberOfAtoms();i++){
      for(int i=0;i<minus_c_freq.GetLength();i++){
	getline(infile,line);
	istringstream iss(line);
	iss >> minus_c_freq[i];
      }
    }


  }

  infile.close();

  if(countgeom != 7){
    printf("ERROR::Quasihamonic::ReadAnisotropicFrequencies() Incorrect number of volumes countgeom = %i\n",
	   countgeom);
    exit(1);
  }
  

  //for super cells make sure that all volumes have the same number of real freqencies
  // or for gamma point frequencies, no more than 3 non-positive frequencies (translational frequencies)
  //start with a
  int imag_num=0;
  int imag_plus=0;
  int imag_minus=0;

  for(int i = 0;i < freq_ref.GetLength(); i++)
    if(freq_ref[i] < 0.0)
       imag_num++;
  
  for(int i = 0; i < plus_a_freq.GetLength(); i++)
    if(plus_a_freq[i] < 0.0)
      imag_plus++;
   
  if(imag_num != imag_plus && (freq_ref.Nonpositive() > 3 && plus_a_freq.Nonpositive() > 3) ){
   printf("ERROR::Quasihamonic::ReadAnisotropicFrequencies() Number of real frequencies differ between reference geometry and geometry of increased lattice a.\n");
   exit(0);
  }

  for(int i = 0; i < minus_a_freq.GetLength(); i++)
    if(minus_a_freq[i] < 0.0)
      imag_minus++;

  if(imag_num != imag_minus && (freq_ref.Nonpositive() > 3 && minus_a_freq.Nonpositive() > 3) ){
   printf("Quasihamonic::ReadAnisotropicFrequencies::Error::Number of real frequencies differ between reference geometry and geometry of decrease lattice a.\n");
   exit(0);
  }  

  //reset for b
  imag_plus=0;
  imag_minus=0;

  for(int i = 0; i < plus_b_freq.GetLength(); i++)
    if(plus_b_freq[i] < 0.0)
      imag_plus++;

  if(imag_num != imag_plus && (freq_ref.Nonpositive() > 3 && plus_b_freq.Nonpositive() > 3) ){
   printf("ERROR::Quasihamonic::ReadAnisotropicFrequencies() Number of real frequencies differ between reference geometry and geometry of increased lattice b.\n");
   exit(0);
  }

  for(int i = 0; i < minus_b_freq.GetLength(); i++)
    if(minus_b_freq[i] < 0.0)
      imag_minus++;

  if(imag_num != imag_minus && (freq_ref.Nonpositive() > 3 && minus_b_freq.Nonpositive() > 3) ){
   printf("ERROR::Quasihamonic::ReadAnisotropicFrequencies() Number of real frequencies differ between reference geometry and geometry of decreased lattice b.\n");
   exit(0);
  }  

  //reset for c
  imag_plus=0;
  imag_minus=0;

  for(int i = 0; i < plus_c_freq.GetLength(); i++)
    if(plus_c_freq[i] < 0.0)
      imag_plus++;

  if(imag_num != imag_plus  && (freq_ref.Nonpositive() > 3 && plus_c_freq.Nonpositive() > 3) ){
   printf("ERROR::Quasihamonic::ReadAnisotropicFrequencies() Number of real frequencies differ between reference geometry and geometry of increased lattice c.\n");
   exit(0);
  }

  for(int i = 0; i < minus_c_freq.GetLength(); i++)
    if(minus_b_freq[i] < 0.0)
      imag_minus++;

  if(imag_num != imag_minus && (freq_ref.Nonpositive() > 3 && plus_c_freq.Nonpositive() > 3) ){
   printf("Quasihamonic::ReadAnisotropicFrequencies::Error::Number of real frequencies differ between reference geometry and geometry of decreased lattice c.\n");
   exit(0);
  }  


  /*
  printf("ref a = %f b = %f c = %f\n",a_ref,b_ref,c_ref);
  freq_ref.Print("reference frequencies");
  printf("PlusA = %f\n",plus_a);
  plus_a_freq.Print("plus_a_freq");
  printf("MinusA = %f\n",minus_a);
  minus_a_freq.Print("minus_a_freq");
 printf("PlusB = %f\n",plus_b);
  plus_b_freq.Print("plus_b_freq");
  printf("MinusB = %f\n",minus_b);
  minus_b_freq.Print("minus_b_freq");
 printf("PlusC = %f\n",plus_c);
  plus_c_freq.Print("plus_c_freq");
  printf("MinusC = %f\n",minus_c);
  minus_c_freq.Print("minus_c_freq");
  exit(0);
  */
  //Gruneisen_init = true;

}


//GrunseisenParameter = -(dLn(freq)/dLn(V)) or -(dLn(freq)/dLn(u)) for lattice lengths
void Quasiharmonic::SetGruneisenParameters(Vector& FreqPlus,double& VolumePlus,Vector& FreqMinus,double& VolumeMinus, int Param){

  //Check to make sure all the sets of freqiencies vectors are the same length.
  //This does not check if the frequency in the frequency files have the right number of entries
  if(FreqPlus.GetLength() != freq_ref.GetLength() ){
    printf("Quasihamonic::SetGruneisenParameters : FreqPlus vector incorrect length \n  Length = %i instead of %i\n",
	   FreqPlus.GetLength(),freq_ref.GetLength());
    exit(1);
  }
  if(FreqMinus.GetLength() != freq_ref.GetLength()){
    printf("Quasihamonic::SetGruneisenParameters : FreqMinus vector incorrect length \n  Length = %i instead of %i \n",
	   FreqMinus.GetLength(),freq_ref.GetLength());
    exit(1);
  }

  //check that volumes are non-zeros
  if(VolumePlus == 0 ||VolumeMinus == 0){
   printf("Quasihamonic::SetGruneisenParameters : Volumes have zero value VolumePlus = %f VolumeMinus = %f  \n",
	  VolumePlus,VolumeMinus);
    exit(1);
  }
  //check that there is a change in volume
  if(fabs(VolumePlus - VolumeMinus) < 0.000001){
   printf("Quasihamonic::SetGruneisenParameters : No change in volume Volume = %f  \n",
	  VolumePlus);
    exit(1);
  }
    
  //Set volume base, isotropic gruneisen parameters
  if(Param == 0){
    //for(int i=0; i<3*GetTotalNumberOfAtoms(); i++)
    for(int i=0; i<Gruneisen_parameters.GetLength();i++)
      if(freq_ref[i] < 0.00001 ||FreqPlus[i] < 0.00001 || FreqMinus[i] < 0.00001)
	Gruneisen_parameters[i] = 0.0; 
      else 
	Gruneisen_parameters[i] = log(FreqPlus[i]/FreqMinus[i]);
    
    Gruneisen_parameters.Scale(1/log(VolumeMinus/VolumePlus));

   //print gruneisen parameters
    FILE *grun_file;
    string filename = Cluster::cluster().GetInputFilename();
    string grun_name = filename.substr(0,filename.size()-3);
    grun_name += ".grun";
    if ((grun_file = fopen(grun_name.c_str(),"w"))==NULL){
      printf("Cluster::PrintVibrationalFrequencies : Cannot open file %s\n",grun_name.c_str());
      exit(1);
    }
    fprintf(grun_file,"%s\n","Gruneisen Parameters");
    for(int i=0;i<Gruneisen_parameters.GetLength();i++){
     fprintf(grun_file,"    %10.6f\n",Gruneisen_parameters[i]);
    }
    fprintf(grun_file,"\n");
    fclose(grun_file);
    Gruneisen_parameters.Print("Gruneisen_parameters");
  }
  //set Gruneisen parameters for a
  else if(Param == 1){
    for(int i=0; i<Grun_a.GetLength();i++)
      if(freq_ref[i] < 0.00001 ||FreqPlus[i] < 0.00001 || FreqMinus[i] < 0.00001)
	Grun_a[i] = 0.0; 
      else 
	Grun_a[i] = log(FreqPlus[i]/FreqMinus[i]);
    
    Grun_a.Scale(1/log(VolumeMinus/VolumePlus));
    Grun_a.Print("Grun_a");
    printf("\n");
  }
  //set Gruneisen parameters for b
  else if(Param == 2){
    for(int i=0; i<Grun_b.GetLength();i++)
      if(freq_ref[i] < 0.00001 ||FreqPlus[i] < 0.00001 || FreqMinus[i] < 0.00001)
	Grun_b[i] = 0.0; 
      else 
	Grun_b[i] = log(FreqPlus[i]/FreqMinus[i]);
    Grun_b.Scale(1/log(VolumeMinus/VolumePlus));
    Grun_b.Print("Grun_b");
    printf("\n");
  }
  //set Gruneisen parameters for c
  else if(Param == 3){
    for(int i=0; i<Grun_b.GetLength();i++)
      if(freq_ref[i] < 0.00001 ||FreqPlus[i] < 0.00001 || FreqMinus[i] < 0.00001)
	Grun_c[i] = 0.0; 
      else 
	Grun_c[i] = log(FreqPlus[i]/FreqMinus[i]);
    
    Grun_c.Scale(1/log(VolumeMinus/VolumePlus));
    Grun_c.Print("Grun_c");           
    printf("\n");
  }
  //set Gruneisen parameters for lattice angles
  else if(Param >= 3 && Param <= 6){
    printf("Error::Quasiharmonic::SetGruneisenParameters():: Gruneisen parameters not yet available for lattice length\n Param = %i\n",
	   Param);
    exit(0);
  }else{
    printf("Error::Quasiharmonic::SetGruneisenParameters():: Param not reconized Param = %i\n",Param);
    exit(0);
  }
  fflush(stdout);
  
}

//Isotropic QuasiHarmonic frequencies found in the cluster class
Vector Quasiharmonic::ComputeAnisotropicFrequencies(){

  //current lengths
  double a = Cluster::cluster().GetUnitCellAxes()[0];
  double b = Cluster::cluster().GetUnitCellAxes()[1];
  double c = Cluster::cluster().GetUnitCellAxes()[2];

  Vector Freqs(freq_ref.GetLength());
  
  for(int i=0;i<freq_ref.GetLength();i++)
   if(fabs(Grun_a[i]) < 1.0e-12 || fabs(Grun_b[i]) < 1.0e-12 ||fabs(Grun_c[i]) < 1.0e-12) 
    Freqs[i] = 0.0;
   else
    Freqs[i] = freq_ref[i]*pow(a/a_ref,-Grun_a[i])*pow(b/b_ref,-Grun_b[i])*pow(c/c_ref,-Grun_c[i]);

  //Freqs.Print("Freqs");

  return Freqs;

}


//reading reference frequencies and modes from molden file
void Quasiharmonic::ReadMoldenFrequencies(){

  Modes_ref.Set();
  freq_ref.Set();

  //Open the frequency file
  string freq_name = Cluster::cluster().GetInputFilename();
  freq_name = freq_name.substr(0,freq_name.size()-3);
  freq_name += "-molden.out";
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("ERROR:: Quasiharmonic::ReadMoldenFrequencies: Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  } 

  int iMode = 1;
  string line;
  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
  while ( !infile.eof() ) {

    getline(infile,line);

    //read frequency from molden file
    if(line == "[FREQ]"){
      for(int ifreq = 0;ifreq < 3*Natoms; ifreq++){
	getline(infile,line);
      	istringstream iss(line);
	iss >> freq_ref[ifreq];
      }
    }

    //header for the mode is "vibration %i" were %i is a number
    char ModeHeader[13];
    sprintf(ModeHeader,"vibration %d",iMode);

    //read vibrational modes from molden file
    if(ModeHeader == line){
      //printf("in if %s\n\n",ModeHeader);
      for(int iatom = 0;iatom < Natoms; iatom++){
	getline(infile,line);
	istringstream iss(line);	
	iss >> Modes_ref(3*iatom,iMode-1);
	iss >> Modes_ref(3*iatom+1,iMode-1);
	iss >> Modes_ref(3*iatom+2,iMode-1);
      }
      iMode++;
    }

  }


  infile.close();
  freq_ref.Print("freq_ref");
  //Modes_ref.PrintHessian("Modes_ref");
  //exit(0);
}

//reading frequencies from molden files for reference supercell frequencies and modes at a given k point
void Quasiharmonic::ReadMoldenFrequencies(int k_point){

  Modes_ref.Set();
  freq_ref_k.Set();
  
  //path for the reference molden files if saved
  string path = Params::Parameters().Params::Parameters().GetHessianPath();
  
  //Path of the reference molden files if saved
  //if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
  //    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  char k_string[3];
  sprintf(k_string,"%d",k_point);

  string freq_name = path + "/molden-" + k_string + ".out";
  //Open the frequency file 
  ifstream infile;
  infile.open( freq_name.c_str() );
  if ( !infile.is_open() ) {
    printf("ERROR:: Quasiharmonic::ReadMoldenFrequencies: Cannot open file '%s'\n",
	   freq_name.c_str());
    exit(1);
  } 

  int iMode = 1;
  string line;
  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
  while ( !infile.eof() ) {

    getline(infile,line);

    //read frequency from molden file
    if(line == "[FREQ]"){
      for(int ifreq = 0;ifreq < 3*Natoms; ifreq++){
	getline(infile,line);
      	istringstream iss(line);
	iss >> freq_ref_k[ifreq];
      }
    }

    //header for the mode is "vibration %i" were %i is a number
    char ModeHeader[13];
    sprintf(ModeHeader,"vibration %d",iMode);

    //read vibrational modes from molden file
    if(ModeHeader == line){
      //printf("in if %s\n\n",ModeHeader);
      for(int iatom = 0;iatom < Natoms; iatom++){
	getline(infile,line);
	istringstream iss(line);	
	iss >> Modes_ref(3*iatom,iMode-1);
	iss >> Modes_ref(3*iatom+1,iMode-1);
	iss >> Modes_ref(3*iatom+2,iMode-1);
      }
      iMode++;
    }

  }

  //freq_ref_k.Print("freq_ref_k");
  //Modes_ref.PrintHessian("Modes_ref");
}


//Isotropic lattice gradient handled in cluster.C
Vector Quasiharmonic::ComputeAnisotropicGradient(){

  double Temperature = Params::Parameters().GetTemperature();
  double h_kT = h_planck/(k_boltz*Temperature);


  //lattice length in Bohr
  double a = Cluster::cluster().GetUnitCellAxes()[0];
  double b = Cluster::cluster().GetUnitCellAxes()[1];
  double c = Cluster::cluster().GetUnitCellAxes()[2];

  //Freq in cm^-1
  Vector Freqs = Cluster::cluster().GetFrequency();

  //Helmholtz gradient
  Vector dF_dp(6);
  

  //dF/dP = sum(dF/df * df/dp)
  for(int i=0;i<Freqs.GetLength(); i++){
    if(Freqs[i] > 0.0 ||fabs(Grun_a[i]) < 1.0e-12 || fabs(Grun_b[i]) < 1.0e-12 ||fabs(Grun_c[i]) < 1.0e-12){
      double freq_Hz = Freqs[i]*WaveNumberToFreq;
      double freq_ref_Hz =  freq_ref[i]*WaveNumberToFreq;
      
      //derivative of Vibrational Energy with respect to frequency
      //  dF/df
      double dF_df = 0.5 + 1/(exp(h_kT*freq_Hz) - 1);
      dF_df *= h_planck/HartreesToJ;
      
      //printf("\ni = %i\n",i);
      // printf("f in cm^-1 = %f\n",Freqs[i]);
      //printf("f in Hz = %f\n",Freqs_Hz);
      //printf("dF_df = %f/n\n",dF_df);


      //derivative of frequency with respect to lattice params
      //  df/dp
      Vector df_dp(6);
      df_dp[0] = -freq_ref_Hz*Grun_a[i]*pow(a,-Grun_a[i]-1)*pow(a_ref,Grun_a[i])*pow(b/b_ref,-Grun_b[i])*pow(c/c_ref,-Grun_c[i])*BohrToAng;
      df_dp[1] = -freq_ref_Hz*Grun_b[i]*pow(b,-Grun_b[i]-1)*pow(b_ref,Grun_b[i])*pow(a/a_ref,-Grun_a[i])*pow(c/c_ref,-Grun_c[i])*BohrToAng;
      df_dp[2] = -freq_ref_Hz*Grun_c[i]*pow(c,-Grun_c[i]-1)*pow(c_ref,Grun_c[i])*pow(a/a_ref,-Grun_a[i])*pow(b/b_ref,-Grun_b[i])*BohrToAng;
      //currently not including gradient for lattice angles
      df_dp[3] = 0.0;
      df_dp[4] = 0.0;
      df_dp[5] = 0.0;
      
      //printf("freq_ref cm^-1 = %f\n",freq_ref[i]);
      //printf("freq_ref Hz = %f\n",freq_ref_Hz);
      //printf("Grun_a = %f\n",Grun_a[i]);
      //printf("

      //combining chain rule product
      dF_dp[0] += dF_df * df_dp[0];
      dF_dp[1] += dF_df * df_dp[1];
      dF_dp[2] += dF_df * df_dp[2];
      dF_dp[3] += dF_df * df_dp[3];
      dF_dp[4] += dF_df * df_dp[4];
      dF_dp[5] += dF_df * df_dp[5];

    }
  }

  //normalizing if doing lattice dynamics 
  if(Params::Parameters().IsSupercellJob())  
    dF_dp.Scale(1/double(Params::Parameters().GetSuperCellTotalKPoint()));

  return dF_dp;
}

//Orders frequencies based on mode overlap with reference modes
Vector Quasiharmonic::OrderFrequencies(Vector OriginalFreq,Matrix Modes){


  //Modes.PrintHessian("Modes");

  //if using supercell frequencies, use the reference frequencies for a particular k point
  Vector reference_freq;
  if(Params::Parameters().IsSupercellJob())
    reference_freq = freq_ref_k;
  else
    reference_freq = freq_ref;


  /*reference_freq.Print("reference_freq");
  printf("\n");
  OriginalFreq.Print("OriginalFreq");
  printf("\n");*/

  /*printf("size reference_freq = %i\tsize gamma reference_freq = %i\n",reference_freq.GetLength(),freq_ref.GetLength());
  fflush(stdout);
  exit(1);*/

  //Make sure that Freq_ref and Freq are the same length
  if(reference_freq.GetLength() != OriginalFreq.GetLength()){
    printf("Error::Quasiharmonic::OrderFrequecies():: Wrong number of frequencies Length= %i\n",
            OriginalFreq.GetLength());
    exit(0);
  }

  //Make sure Ref_Modes and Modes are the same size. 
  if(Modes.GetRows() != Modes_ref.GetRows() || Modes.GetCols() != Modes_ref.GetCols()){
    printf("Error::Quasiharmonic::OrderFrequencies():: Modes wrong size Row = %i Col = %i\n",
           Modes.GetRows(),Modes.GetCols());
    exit(0);
  }

  if(Modes.GetRows() != OriginalFreq.GetLength() || Modes.GetCols() != OriginalFreq.GetLength()){
    printf("Error::Quasiharmonic::OrderFrequencies():: Mismatch Modes with Freq size: Modes Rows = %i Col = %i,Freq Length = %i\n",
           Modes.GetRows(),Modes.GetCols(),OriginalFreq.GetLength());
  }

  /*Modes_ref.PrintHessian("Modes_ref");
  printf("\n");
  Modes.PrintHessian("Modes");
  fflush(stdout);*/

  //printf("size reference_freq = %i\tsize OrigFreq = %i\n",reference_freq.GetLength(),OriginalFreq.GetLength());
  //fflush(stdout);

  Vector Freq(OriginalFreq.GetLength());
  Freq.Set(-1000000.0); //default value
  int count_imaginary = 0;
  Matrix overlap(OriginalFreq.GetLength(),OriginalFreq.GetLength());
  //OriginalFreq.Print("Original Freq");
  //Location of original frequencies on Freq vector 
  /*int entries[OriginalFreq.GetLength()];
  for(int i = 0;i < OriginalFreq.GetLength();i++){
    entries[i] = -1;
    //printf("entries[%i] = %i\n",i,entries[i]);
  }*/

  //For each frequency mode-match against all available reference freqs.
  //Perhaps restrict this to go only over a certain freq range (100 cm^-1?)
  for(int i=0;i<OriginalFreq.GetLength();i++){
  //if(OriginalFreq[i] > 0.0001){
    //printf("i = %i Freq Length = %i\n",i,Freq.GetLength());
    //printf("i = %i Freq = %f\n",i,OriginalFreq[i]);
    if(OriginalFreq[i] < 0.0){
      count_imaginary += 1;
    }
    

    Vector ModeVector = Modes.GetColumnVector(i);
      
    //Testing not doing this JLM
    //Make first nonzero entry positive so modes are consistent
    /*int sign_spot = -1 ;
    for(int k=0;k<ModeVector.GetLength();k++){
      if(fabs(ModeVector[k] > 0.0001) && sign_spot == -1){
	sign_spot = k;
	//printf("i = %i sign_spot = %i\n",i,sign_spot);
      }
    }
    //in case of zero vector
    if(sign_spot == -1) sign_spot = 0;
      
    //Make first entry positive so modes are consistent
    if(ModeVector[sign_spot] < 0.0)
      ModeVector.Scale(-1.0);*/

    //Now perform the overlap with all the reference frequencies
    for(int j=0;j<reference_freq.GetLength();j++){
      //printf("i = %i j = %i Freq Length = %i\n",i,j,Freq.GetLength());
      Vector ModeVector_ref = Modes_ref.GetColumnVector(j);

      //JLM Testing not doing this
      //if(ModeVector_ref[sign_spot] < 0.0)
      //  ModeVector_ref.Scale(-1.0);      
	
      //dot product of the modes
      double dot_product;

      //Limit range over which the frequency match can occur??
      if(fabs(fabs(OriginalFreq[i])-fabs(reference_freq[j])) >= 25.0)
        dot_product = 0.0;
      else
        dot_product = ModeVector.DotProduct(ModeVector_ref);

      overlap(i,j) = fabs(dot_product);
      
    }

  }

  //printf("count_imaginary = %i\n",count_imaginary);
  //fflush(stdout);

  //Take care of imaginary freqs
  if(count_imaginary > 0) {
    for(int i=OriginalFreq.GetLength()-1; i>OriginalFreq.GetLength()-1-count_imaginary;i--){
      for(int j=OriginalFreq.GetLength()-1;j>OriginalFreq.GetLength()-1-count_imaginary;j--){
        if(i==j){
         //printf("Setting overlap(%i,%i) = %f\n",i,j,1.0);
         overlap(i,j) = 1.0;
        }
        else{
         //printf("Setting overlap(%i,%i) = %f\n",i,j,0.0);
         overlap(i,j) = 0.0;
        }
      }
    }
  }

  //fflush(stdout);
  //exit(1);
  /*overlap.PrintHessian("Overlap matrix");
  printf("\n");
  fflush(stdout);*/
  //exit(1);

  //Now to solve the assignment
  //Make a lapack call that maximizes to overlap for each frequency mode?
  //Will need to assign along each row (to assign original mode)

  //Trying Stable marriage problem
  //Both the reference and original frequencies will be un assigned
  //Essentially iterate over the reference frequency until each frequency is  
  //optimally assigned to some original frequency
  //The dot products will be pre-tabulated to ensure no additional math need be done
  
  int num_changes_made = OriginalFreq.GetLength();
  int count_iterations = 0;

  //holds the ref freq location the Orig freq is assigned to (int)
  Vector mapping_orig(OriginalFreq.GetLength());
  mapping_orig.Set(-1000);
  //mapping_orig.Print("Mapping orig");
  //printf("\n");
  //holds the orig freq location the ref freq is assigned to (int)
  Vector mapping_ref(OriginalFreq.GetLength());
  mapping_ref.Set(-1000);
  //holds the currently assigned overlaps of the reference freq (double)
  Vector current_ref_overlap(OriginalFreq.GetLength());
  Vector tmp_overlap;
  Vector sorted_tmp_overlap;
  
  //Orig Freqs on the rows
  //Reference Freqs on the columns
  while(num_changes_made > 0){
    num_changes_made = 0;
    //printf("Iteration: %i\n",count_iterations);
    //Iterate over each unassigned original frequency
    for(int i=0;i<OriginalFreq.GetLength();i++){
    //for(int i=0;i<24;i++){
      //int i = 0;

      //Only need to go through this if orig freq unassigned?
      if(mapping_orig[i] == -1000) {

        //Get overlap row cooresponding to the current original vector
        tmp_overlap = overlap.GetRowVector(i);
        //tmp_overlap.Print("Current Row");
        //fflush(stdout);
      
        //Now find a way to sort this vector...
        sorted_tmp_overlap = tmp_overlap;
        Vector map_ref_from_sorting = sorted_tmp_overlap.SortLargestFirstAndReturn();
        //sorted_tmp_overlap.Print("Sorted Row");

        /*if(count_iterations==3){
            printf("Attempting to assign freq %i!\n",mapping_orig[i]);
           // printf("potential_overlap  = %f\tcurrent overlap = %f\n",potential_overlap,current_overlap);
          sorted_tmp_overlap.Print("Freq 147");
           
        }*/

        double potential_overlap;
        double current_overlap;

        //Check the reference frequency for optimal overlap
        for(int j=0;j<sorted_tmp_overlap.GetLength();j++){

          potential_overlap = sorted_tmp_overlap[j];
          //since we sorted the list we need to 'unsort' it to correctly
          //map the reference frequency overlap to the proposed frequency overlap
          //int actual_location = tmp_overlap.Find(potential_overlap);
          //Possibly get rid of this if the call works and use only map_ref_from_sorting
          //int actual_location = int(map_ref_from_sorting[j]);
          //printf("actual_location  = %i\n",map_ref_from_sorting[j]);
          if(mapping_ref[map_ref_from_sorting[j]] != -1000)
            current_overlap = current_ref_overlap[map_ref_from_sorting[j]];
          else
            current_overlap = 0.0;
          
          //printf("potential_overlap  = %f\tcurrent overlap = %f\n",potential_overlap,current_overlap);

          //fflush(stdout);
          //exit(1);
          /*if(count_iterations==3){
            //printf("Attempting to assign freq %i!\n",mapping_orig[i]);
            printf("Ref freq = %i\n",map_ref_from_sorting[j]);
            printf("potential_overlap  = %f\tcurrent overlap = %f\n",potential_overlap,current_overlap);
            printf("mapping_ref[map_ref_from_sorting[j]] = %f\n",mapping_ref[map_ref_from_sorting[j]]);
            
          }*/

          if(potential_overlap >= current_overlap){
            /*printf("A better potential overlap was found!\n");
            printf("potential_overlap  = %f\tcurrent overlap = %f\n",potential_overlap,current_overlap);
            printf("Assigning ref freq %i to orig freq %i\n",map_ref_from_sorting[j],i);
            printf("Orig freq %i used to be assigned to ref freq %i, will now be assigned to ref freq %f\n",i,mapping_orig[i],map_ref_from_sorting[j]);
            printf("Ref freq %i used to be assigned to orig freq %f, will now be assigned to orig freq %i\n",map_ref_from_sorting[j],mapping_ref[map_ref_from_sorting[j]],i);
            printf("Previous Ref freq %i assigned to orig freq %f, will now be assigned to %i\n",map_ref_from_sorting[j], mapping_ref[map_ref_from_sorting[j]],-1000);
            printf("Logic check actual_location: %i mapping_ref[map_ref_from_sorting[j]]: %f mapping_orig[mapping_ref[map_ref_from_sorting[j]]]: %f\n",map_ref_from_sorting[j],mapping_ref[map_ref_from_sorting[j]],mapping_orig[mapping_ref[map_ref_from_sorting[j]]]);
            printf("Previous Ref overlap %f will now be assigned to %f\n",current_ref_overlap[map_ref_from_sorting[j]],potential_overlap);*/

            /*mapping_ref.Print("mapping_ref before reassignment");
            printf("\n");
            mapping_orig.Print("mapping_orig before reassignment");
            printf("\n");
            current_ref_overlap.Print("current_ref_overlap before reassignment");
            printf("\n");*/
            bool assign_freq = false;
            //If overlap with this ref freq is greater than the current overlap
            //then reassign the orig freq to this current ref.
            if(potential_overlap > current_overlap){
              assign_freq = true;
            }
            //If overlap with this ref freq is equal to the current overlap
            //and ref freq is unassigned then assign to the current orig freq
            else if(potential_overlap == current_overlap && int(mapping_ref[map_ref_from_sorting[j]]) == -1000){
              assign_freq = true;
            }
            //If overlap with this ref freq is equal to the current overlap
            //then assign the closer frequency value to this assignment
            else{
              double current_freq_diff;
              if(mapping_ref[map_ref_from_sorting[j]] != 1000)
                current_freq_diff = fabs(fabs(OriginalFreq[mapping_ref[map_ref_from_sorting[j]]])-fabs(reference_freq[map_ref_from_sorting[j]]));
              else
                current_freq_diff = 100000;

              double potential_freq_diff = fabs(fabs(OriginalFreq[i])-fabs(reference_freq[map_ref_from_sorting[j]]));
              if(potential_freq_diff < current_freq_diff){
                assign_freq = true;
              }
            }

            if(assign_freq){
              //Unassign old spot
              if(int(mapping_ref[map_ref_from_sorting[j]])!=-1000)
                mapping_orig[int(mapping_ref[map_ref_from_sorting[j]])] = -1000;
              //Assign new spots
              mapping_orig[i] = int(map_ref_from_sorting[j]);
              mapping_ref[map_ref_from_sorting[j]] = i;
              //store the new reference overlap
              current_ref_overlap[map_ref_from_sorting[j]] = potential_overlap;
              num_changes_made++;
              break;
            }

            /*printf("\nLocation 9: ref_mapping = %f orig_freq_mapping = %f current_ref_overlap = %f\n",mapping_ref[9],mapping_orig[9],current_ref_overlap[9]);
            printf("Location 10: ref_mapping = %f orig_freq_mapping = %f current_ref_overlap = %f\n\n\n",mapping_ref[10],mapping_orig[10],current_ref_overlap[10]);*/


            /*mapping_ref.Print("mapping_ref after reassignment");
            printf("\n");
            mapping_orig.Print("mapping_orig after reassignment");
            printf("\n");
            current_ref_overlap.Print("current_ref_overlap after reassignment");
            printf("\n");
            fflush(stdout);*/
            //Update changes made count
          }

        } //End loop over the potential reference overlaps
      } // End if statement checking if the reference freq is unassigned
    }
    
    //Check to see if iteration count has been exceeded
    if(count_iterations == 50){
      printf("Unable to assign all frequencies in 50 iterations. Exiting...\n");
      exit(1);
    }
    else {
      count_iterations++;
    }

    //printf("num_changes_made: %i\n",num_changes_made);
    
    //Ensure all freqs are assigned
    //if(num_changes_made==0 && )
  }
  
  //printf("Frequencies were assigned in %d iterations.\n\n",count_iterations);
  /*current_ref_overlap.Print("Current ref overlap");
  printf("\n");
  mapping_orig.Print("Orig freq mapping to ref freq location");
  printf("\n");
  mapping_ref.Print("Reference freq mapping to orig freq location");
  printf("\n");
  OriginalFreq.Print("Original Freqs");*/
  //printf("\nExiting...\n");
  //exit(1);

  //Once this is complete then assign original frequencies to the final Freq vector;

  //Now assign the new frequencies
  for(int i = 0; i<Freq.GetLength();i++){
    Freq[i] = OriginalFreq[mapping_ref[i]];
  }
  
  Vector check_freq = Freq;
  check_freq -= reference_freq;
  /*current_ref_overlap.Print("Current ref overlap");
  printf("\n");
  Freq.Print("New Freq assignments");
  printf("\n");
  reference_freq.Print("Reference Freq");
  printf("\n");
  check_freq.Print("Check Freq");
  printf("\n");
  fflush(stdout);
  exit(1);*/

  //making sure every spot in vector Freq is filled.
  for(int i = 0; i<Freq.GetLength();i++){
    if(fabs(Freq[i] + 1000000.0) < 0.000001){
      printf("Error::Quasiharmonic::OrderFrequencies() No frequency placed in entry %i in frequency vector\n",
	     i);
      Freq.Print("Freq");
      exit(0);
    }
  }
  
  //Vector error = reference_freq;
  //error -= Freq;
  /*for(int i = 0; i<Freq.GetLength();i++){
    printf("i = %i\treference_freq[i] = %f\tFreq[i] = %f\tError = %f\n",i,reference_freq[i],Freq[i],fabs(reference_freq[i])-fabs(Freq[i]));
  }
  fflush(stdout);*/

  return Freq;
  
}


//if two or more freqs have the lowest RMS modes with the same ref mode.
//this algorithm should find the placement of those modes.
void Quasiharmonic::MangageModeVectors(Matrix& Modes,Vector& OriginalFreq,Vector& Freq,vector<int> entered_entries){

  //Freq.Print("Freq at beginning of ManageModeVector()");


  //if using supercell frequencies, use the reference frequencies for a particular k point
  Vector reference_freq;
  if(Params::Parameters().IsSupercellJob())
    reference_freq = freq_ref_k;
  else
    reference_freq = freq_ref;

  //reference_freq.Print("Reference freq");


  //disbuted entry is the first entry in "entered_entries"
  Vector disbuted_mode = Modes_ref.GetColumnVector(entered_entries[0]);

  //Make first nonzero entry positive so modes are consistant
  int sign_spot = -1;
  for(int k=0;k<disbuted_mode.GetLength();k++)
    if(fabs(disbuted_mode[k] > 0.0001) && sign_spot == -1){
      sign_spot = k;
      //printf("entered_entries[0] = %i sign_spot = %i\n",entered_entries[0],sign_spot);
    }
  //in case of zero vector
  if(sign_spot == -1) sign_spot = 0;

  //printf("disbuted_mode %i",entered_entries[0]);
  //disbuted_mode.Print("");
  
  //entry with the lowest RMS
  //double lowest_RMS = 1000.0;

  //entry with the highest Dot Product
  double dot_max = 0.0;
  
  //entry which belongs in the disputed entry
  int right_entry = -1;

  //printf("Quasiharmonic::MangageModeVectors 3\n");
  
  int num_entries = entered_entries.size() - 1;
  int entry[num_entries];
  
  //determining which has the lowest RMS to the reference modes
  for(int i = 0; i < num_entries; i++){
    entry[i] = entered_entries[i+1];
    
    Vector thisMode = Modes.GetColumnVector(entry[i]);
    if(thisMode[sign_spot] <0)
      thisMode.Scale(-1.0);

    double dot = fabs(thisMode.DotProduct(disbuted_mode));
    if(dot > dot_max){
      dot_max = dot;
      right_entry = entry[i];
    }

    //mode of entry i
    //Vector diff_vector = Modes.GetColumnVector(entry[i]);
    //if(diff_vector[sign_spot] <0)
    // diff_vector.Scale(-1.0);
    
    //for(int j = 0;j<diff_vector.GetLength();j++){
    //  diff_vector[j] -= disbuted_mode[j];
    //}
    
    //double RMS = diff_vector.RMS();

    
    //printf("i = %i entry = %i,Freq = %f dot = %f dot_max = %f\n",
    //	   i,entry[i],Freq[i], dot, dot_max);
    
    //if(RMS < lowest_RMS){
    //  lowest_RMS = RMS;
    //  right_entry = entry[i];
    //}
  }
  
  if(right_entry == -1){
    printf("ERROR::Quasiharmonic::ManageModeVector() Unable to find frequency for entry  %i\n",
	   entered_entries[0]);
    
    exit(0);
  }

  Freq[entered_entries[0]] = OriginalFreq[right_entry];
  //Freq.Print("Freq after finding right_entry");

  
  // find placement for rest. Fill only unfilled spot (marked by being filled by 1000000.0
  for(int i = 0; i < num_entries;i++){
    //printf("entry= %i right_entry = %i\n",entry[i],right_entry);
    if(right_entry != entry[i]){

      //since the modes should have simpler RMS to the reference Modes, use difference in freq instead.
      //May not be the best way to handle this
      int min_freq_location = -1;
      //double min_freq_diff = 100000.0;
      double min_freq_diff = 100.0;
      
      double thisFreq = OriginalFreq[entry[i]];

      for(int j = 0; j < reference_freq.GetLength();j++){
        double otherFreq = reference_freq[j];

        //if(entry[i] == 64){
	//printf("entry[%i] = %i,j = %i, thisFreq = %f otherFreq = %f\n",i,entry[i],j,thisFreq,otherFreq);
	//printf("Freq[%i] == %f\n",j,Freq[j]);
	//printf("diff = %f min_freq_diff = %f min_freq_location = %i\n\n",fabs(thisFreq - otherFreq),min_freq_diff,min_freq_location);
	//printf("Freq[j] = %f\n\n",Freq[j]);
        //}
        if(fabs(Freq[j] + 1000000.0) < 0.000001 && (fabs(thisFreq - otherFreq) < min_freq_diff) ){
          min_freq_diff = fabs(thisFreq - otherFreq);
          min_freq_location = j;
	  //printf("Freq %f in %i\n",Freq[j],min_freq_location);
        }
      }
      
      if(min_freq_location == -1){
        printf("ERROR::Quasiharmonic::ManageModeVector() : Unable to find placement for Freq %f\n", thisFreq);
        exit(0);
      }
      Freq[min_freq_location] = OriginalFreq[entry[i]];
      //printf("Freq[%i] now contains %f\n",min_freq_location,Freq[min_freq_location]);
    }

  }

  //Freq.Print("Freq at end of ManageModeVector()");
  //exit(0);
}

Vector Quasiharmonic::ReadGeometryFromInput(double& volume){

  int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  Vector coords(3*UniqueAtoms+6);

  string input_name = hess_type + ".in";
  ifstream infile;
  infile.open( input_name.c_str() );
  if ( !infile.is_open() ) {
    printf("ERROR:: Quasiharmonic::ReadGeometryFromInput: Cannot open file '%s'\n",
	   input_name.c_str());
    exit(1);
  } 

  //printf("input_name = %s\n",input_name.c_str());

  int index = 0;
  int atom = 0;
  bool molec_sxn = false;

  //read input file
  int Ntot = Cluster::cluster().GetTotalNumberOfAtoms();

  string line;

  //read coordinates
  while ( !infile.eof() ) {
    getline(infile,line);
    // Identify molecule section, and read up to first monomer
    if (line.substr(0,9)=="$molecule") {
      molec_sxn = true;
    }
    if (line.substr(0,2)=="--" && molec_sxn==true) {
      index++; // increment monomer counter
      //discard charge of crystal
      getline(infile,line);
      for(int i=0;i<Cluster::cluster().GetNumberOfAtoms(index);i++){
	getline(infile,line);
	istringstream iss(line);
        double x,y,z;
	//if(Cluster::cluster().GetMonomerSymmetryFactor(index) != 0){
	//only count atoms in assymetric unit cell
	if(Cluster::cluster().GetMonomer(index).GetAtom(i).InAsymmetricUnit()){
	  string tmp;
	  //discard atom number if using tinker for the MM
	  if (Params::Parameters().GetMMType()==1)
	    iss >> tmp;
	  //discard atom symbol
	  iss >> tmp;
	  //read atom coodinates
	  iss >> coords[3*atom];
	  iss >> coords[3*atom+1];
	  iss >> coords[3*atom+2];
	  atom++;
	}
	
      }
      // If we have reached the end of the section, break
      if (line.substr(0,4)=="$end" && molec_sxn==true) {
	break;
      }
    }
  }

  //read lattice parameters
  Cluster::cluster().Rewind(infile);
  // Start reading the file
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,10)=="$unit_cell") {
      
      {
        // read in axis lengths
        getline(infile,line);
        istringstream iss(line);
        iss >> coords[3*atom];
        iss >> coords[3*atom+1];
        iss >> coords[3*atom+2];
      }

      {
        // read in angle length
        getline(infile,line);
        istringstream iss(line);
        iss >> coords[3*atom+3];
        iss >> coords[3*atom+4];
        iss >> coords[3*atom+5];
        break;
      }
    }
  }


  infile.close();


  double alpha_rad = coords[3*atom+3]*DegreesToRadians;
  double beta_rad = coords[3*atom+4]*DegreesToRadians;
  double gamma_rad = coords[3*atom+5]*DegreesToRadians;
  double beta_term =
    (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
  double gamma_term =
    sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);

  volume = 
    (coords[3*atom])
    *(coords[3*atom+1]*sin(gamma_rad))
    *(coords[3*atom+2]*gamma_term);

  //coords.PrintGradient("coords");
  return coords;

}
