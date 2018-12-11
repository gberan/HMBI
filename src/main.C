#include "main.h"
using namespace std;

/*
HMBI code.  Invoke as:
      hmbi <input file> [# of processors]

The number of processors is optional, and the default is 1.  To run
in code in parallel, the source code must be compiled with -DPARALLEL,
and MPI is required.
*/

main(int argc, char ** argv) {

  if (argc < 2) {
    printf("HMBI syntax:: hmbi <input file> [ # of processors (optional)]\n");
    exit(1);
  }
  
  int nproc = 1;
  if (argc > 2) {
    string tmp = argv[2];
    nproc = atoi(tmp.c_str());
  }

#ifdef PARALLEL
  int mynode=0,totalnodes;

  if (nproc > 1) {
    MPI_Init(&argc, &argv);
    MPI_Comm io_comm;
    
    MPI_File *fh, *fg;      
    MPI_Info info;       

    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode
    if (mynode==0) 
      printf("Running MPI parallel version with %d processors\n",nproc);

  }

  if (mynode == 0) {  
#endif /* PARALLEL */

    // Set the number of processors
    Params::Parameters().SetNumberOfProcessors(nproc);

    // Store the primary process ID of the job.  Used to create unique
    // scratch files
    int pid = getpid();
    Params::Parameters().SetPID(pid);

    //JLM
    /*stringstream ss;
    ss << pid;
    string str = ss.str();
    string tmp_str = "/scratch/$USER/" + str; 
    printf("scratch path set to: %s\n",tmp_str.c_str());
    Params::Parameters().SetScratchPath(tmp_str.c_str());*/

    // Start wall clock timer
    time_t start_time, stop_time;
    start_time = time(NULL);

    // Open input file
    const int maxlength = 50;
    char filename[maxlength];
    strcpy(filename, argv[1]);

    ifstream infile;
    infile.open(filename);
    assert(infile.is_open()); 

    // Initialize Cluster object
    Cluster::cluster().Initialize(infile,nproc,filename);
    printf("Cluster initialization complete\n");

    //Entering Start time into the cluster
    Cluster::cluster().SetStartTime(start_time);
    bool supercell_Init = false;
    Supercell::supercell().SetSupercellInit(false);
    //SetInfileName(infile);
    //SetFilename(filename);
    //yoni: test Create supercell input
    if(Params::Parameters().CreateSupercellInput()){
      Supercell::supercell().Initialize(infile,nproc,filename);
      Supercell::supercell().CreateSupercellInputFile();
      exit(0);
    }

    //If analyze_only = true, do not optimize
    if(!Params::Parameters().RunJobs())
      Params::Parameters().SetParameter("MAX_OPT_CYCLES","1");

    if(Params::Parameters().IsSupercellJob() && !supercell_Init && Params::Parameters().GetQMType()==8){ //DFTB
      Supercell::supercell().Initialize(infile,nproc,filename);
      supercell_Init = true;
    }

    if(Params::Parameters().DoQuasiHarmonic()){
      
      printf("Doing quasiharmonic\n");

      if(Params::Parameters().IsSupercellJob() && !supercell_Init){
	Supercell::supercell().Initialize(infile,nproc,filename);
        supercell_Init = true;
      }
	
      Quasiharmonic::quasiharmonic().Initialize();

      if(!Params::Parameters().IsPeriodic() && !Params::Parameters().UseFullQMOnly()){
	printf("ERROR::Cannot perform Quasi-harmonic calulations on non-periodic systems.\n");
	exit(0);
	
      }
      if(Params::Parameters().GetMMType() != 1 && Params::Parameters().GetMMType() != 5 &&
	 !Params::Parameters().AreQHFAvailable() && !Params::Parameters().UseFullQMOnly() ){
	printf("ERROR::Quasi-harmonic has not been adapted to any other many-body type other than Tinker or crystal.\n");
	exit(0);
      }
      //supercell frequencies may have imaginary frequencies. making sure that number of them doesn't charge since only real frequencies are used.
      int imag_num=0;

      //gruneisen parameters only use values
      if(Params::Parameters().UseVolumeGruneisen()){
	
	//Quasi-harmonics approximates the dependance of the vibrations on volume. 
	//We need the vibrations at three volumes of the (relaxed) crystal for this approximation.
	
	double VolumePlus;
	double VolumeMinus;
	Vector FreqPlus;
	Vector FreqMinus;
	int total_k_point;//if only gamma point frequencies are used. total_k_point = 1

	
	//Quasi-harmonic using only the gamma point frequencies
	if( !Params::Parameters().IsSupercellJob()){
	  total_k_point = 1;
	  FreqPlus.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());
	  FreqMinus.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());
	}
	//Quasi-harmonic includes frequencies beyond those in the gamma point
	else{
	  total_k_point = Supercell::supercell().GetTotalKPoint();
	  FreqPlus.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);
	  FreqMinus.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);

          //Cluster::cluster().InitializeQuasiHarmonicForSuperCell();

	  //Quasiharmonic::quasiharmonic().ReadMoldenFrequencies(0);
	}
	
	
	if(!Params::Parameters().AreQHFAvailable()){
          //Create the frequency file if it does not exist
	  Quasiharmonic::quasiharmonic().DetermineIsotropicFrequencies(infile,VolumePlus,VolumeMinus,FreqPlus,FreqMinus,imag_num);
          string infile = Cluster::cluster().GetInputFilename();
          string freq_name = infile.substr(0,infile.size()-3);
          freq_name += ".freq";
          printf("\n\nQuasiharmonic frequencies are in %s\n Set ARE_QHF_AVAILABLE = true to use.\n ",freq_name.c_str());
	  exit(0);
	  Params::Parameters().SetParameter("ARE_QHF_AVAILABLE","TRUE");
	}else Quasiharmonic::quasiharmonic().ReadIsotropicFrequencies(FreqPlus,VolumePlus,FreqMinus,VolumeMinus);
		
	Quasiharmonic::quasiharmonic().SetGruneisenParameters(FreqPlus,VolumePlus,FreqMinus,VolumeMinus,0);
	//Quasiharmonic::quasiharmonic().InitializingGruneisen(true);



      }
      
      //seperate grunessian parameters for each lattice parameter
      else{	


	double plus_a,minus_a,plus_b,minus_b,plus_c,minus_c;

	Vector plus_a_freq;
	Vector minus_a_freq;

	Vector plus_b_freq;
	Vector minus_b_freq;

	Vector plus_c_freq;
	Vector minus_c_freq;	

	int total_k_point;//if only gamma point frequencies are used. total_k_point = 1
	int imag_num=0;//supercell frequencies may have imaginary frequencies. making sure that number of them doesn't charge since
	//only real frequencies are used.

	//Quasi-harmonic using only the gamma point frequencies
	if( !Params::Parameters().IsSupercellJob()){
	  total_k_point = 1;
	  plus_a_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());
	  minus_a_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());

	  plus_b_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());
	  minus_b_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());

	  plus_c_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());
	  minus_c_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms());

	}
	//Quasi-harmonic includes frequencies beyond those in the gamma point
	else{
          if(!supercell_Init){
	    Supercell::supercell().Initialize(infile,nproc,filename);
            supercell_Init = true;
          }
	  total_k_point = Supercell::supercell().GetTotalKPoint();
	  plus_a_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);
	  minus_a_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);

	  plus_b_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);
	  minus_b_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);

	  plus_c_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);
	  minus_c_freq.Initialize(3*Cluster::cluster().GetTotalNumberOfAtoms()*total_k_point);

	  //Cluster::cluster().InitializeQuasiHarmonicForSuperCell();
	}
	if(!Params::Parameters().AreQHFAvailable()){
	  //Get frequency dependency of lattice lengths
	  Quasiharmonic::quasiharmonic().DetermineLatticeFrequencies(infile,plus_a,minus_a,plus_a_freq,minus_b_freq,imag_num,0);
	  Quasiharmonic::quasiharmonic().DetermineLatticeFrequencies(infile,plus_b,minus_b,plus_b_freq,minus_b_freq,imag_num,1);
	  Quasiharmonic::quasiharmonic().DetermineLatticeFrequencies(infile,plus_c,minus_c,plus_c_freq,minus_c_freq,imag_num,2);
	  Params::Parameters().SetParameter("ARE_QHF_AVAILABLE","TRUE");
          exit(0);
	}
	else Quasiharmonic::quasiharmonic().ReadAnisotropicFrequencies(plus_a_freq,plus_a,minus_a_freq,minus_a,
								       plus_b_freq,plus_b, minus_b_freq,minus_b,
								       plus_c_freq,plus_c,minus_c_freq,minus_c);
	
	Quasiharmonic::quasiharmonic().SetGruneisenParameters(plus_a_freq,plus_a,minus_a_freq,minus_a,1);
	Quasiharmonic::quasiharmonic().SetGruneisenParameters(plus_b_freq,plus_b, minus_b_freq,minus_b,2);
	Quasiharmonic::quasiharmonic().SetGruneisenParameters(plus_c_freq,plus_c,minus_c_freq,minus_c,3);
	//Quasiharmonic::quasiharmonic().InitializingGruneisen(true);
      }
    }
  
    if (!Params::Parameters().DoForces())
      Params::Parameters().SetUseDLFind(false); 

    // Compute the Energy (and forces, if requested)
    if (!Params::Parameters().UseDLFind()) {
      if ( Params::Parameters().GetJobTypeStr() == "nmr" ) { //JDH
	Cluster::cluster().PredictMagneticProperties(); //JDH 
      } 
      else if (!Params::Parameters().GetUseLocalOptimizer()){
	Cluster::cluster().RunJobsAndComputeEnergy();
      } 
    }
 
    if ( Params::Parameters().PrintLevel() > 0)
      Cluster::cluster().ComputeDistanceMatrix();

    // For debugging: Compute finite difference gradients
    bool do_finite_difference_grad = Params::Parameters().UseFiniteDifferenceGradients();
    if (do_finite_difference_grad && Params::Parameters().DoForces()) {
      
      Vector Grad = Cluster::cluster().GetHMBIGradient();
      Grad.Print("Analytical HMBI Gradient");


      int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
      int NMon = Cluster::cluster().GetNumberOfMonomers();
      Vector Eplus(3*UniqueAtoms), Eminus(3*UniqueAtoms);
      
      int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
      //Vector Eplus(3*Natoms), Eminus(3*Natoms);
      Vector Original_coords = Cluster::cluster().GetCurrentCoordinates();
      Vector Original_unique_coords = Cluster::cluster().GetSymmetryUniqueCoordinates();
      
      
      Vector LatticeGradient(6);
      if (Params::Parameters().DoLatticeVectorFiniteDifference()){
	LatticeGradient.Initialize(9);
      }
      
      Params::Parameters().SetParameter("JOBTYPE","ENERGY");
      
      //Original_Coord.PrintGradient("Original_Coords");
      
      // Do the finite difference
      double delta = 0.001; // step size, in Angstroms
      if(true){


	//Original_coords.PrintGradient("Original_coord");
	Eminus = Cluster::cluster().EnergyOfFiniteDifference(delta,true);
	Eplus =  Cluster::cluster().EnergyOfFiniteDifference(delta,false);
      }
      //Eminus.Print("Eminus");
      //Eplus.Print("Eplus");   


      // Now do finite difference over the lattice parameters, if appropriate
      double deltaTheta = 0.005;

      if(true){
	if (Params::Parameters().IsPeriodic() && !Params::Parameters().DoLatticeVectorFiniteDifference()) {

	  string *params;
	  params = new string[6];
	  params[0] = "a"; params[1] = "b"; params[2] = "c";
	  params[3] = "alpha", params[4] = "beta"; params[5] = "gamma";
	  

	  //Vector Coords = Original_coords; // used to reset coordinates
	  Vector Coords = Original_unique_coords;
          //Coords.PrintGradient("Coords");
          //exit(0);
	  for (int i=0;i<6;i++) {
	    double current = Cluster::cluster().GetUnitCellParameter(params[i]);
	    double new1,new2;

	    //skip if lattice param is fixed or locked under symmetry
	    bool skip = 0;
	    if(i == 1 && Cluster::cluster().IsBLocked())
	      skip = 1;
	    else if(i == 2 && Cluster::cluster().IsCLocked())
	      skip = 1;
	    else if(i == 3 && Cluster::cluster().IsAlphaLocked())
	      skip = 1;
	    else if(i == 4 && Cluster::cluster().IsBetaLocked())
	      skip = 1;
	    else if(i == 5 && Cluster::cluster().IsGammaLocked())
	      skip = 1;

	    if(!skip) {
	      if (i < 3) {
		printf("Shifting coord %s by -%f\n",params[i].c_str(),delta);
		new1 = current - delta;
	      }
	      else {
		printf("Shifting coord %s by -%f\n",params[i].c_str(),deltaTheta);
		new1 = current - deltaTheta;
	      }
	    
	      //Updating Coordinates and running jobs
	      Coords[3*UniqueAtoms+i] = new1;

              //Shifting atoms if necessary to preserve symmetry
              //none of these  atoms shifted should be a degree of freedom.
              Cluster::cluster().MaintainCartesianSymmetry(Coords,true);

	      //Set the new lattice parameters;
	      Cluster::cluster().UpdateLatticeParams(Coords);
	      //Getting the coordinates for the rest of the system base on symmetry
	      Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
	      Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
	      Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
	      
	      /*
		Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,0,1);
		Cluster::cluster().UpdateLatticeParams(Coords);
		Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,1,1); 
		
	      */
	      
	      Coords.PrintGradient("Coords");
	      Cluster::cluster().SetNewCoordinates(Coords);
	      Cluster::cluster().RunJobsAndComputeEnergy();
	      double E1 = Cluster::cluster().GetHMBIEnergy();
	      
	      //Reset Lattice Params
	      printf("Reseting Coords\n");
	      Coords = Original_unique_coords;
	      Cluster::cluster().UpdateLatticeParams(Original_coords);
	      Cluster::cluster().SetNewCoordinates(Original_coords);
	      
	      if (i < 3) {
		printf("Shifting coord %s by +%f\n",params[i].c_str(),delta);
		new2 = current + delta;
	      }
	      else {
		printf("Shifting coord %s by +%f\n",params[i].c_str(),deltaTheta);
		new2 = current + deltaTheta;
	      }
	      Coords[3*UniqueAtoms+i] = new2;
	      

              //Shifting atoms if necessary to preserve symmetry
              //none of these  atoms shifted should be a degree of freedom.
              Cluster::cluster().MaintainCartesianSymmetry(Coords,true);
	      
	      Cluster::cluster().UpdateLatticeParams(Coords);
	      //Getting the coordinates for the rest of the system base on symmetry
	      Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,false,true);
	      Coords = Cluster::cluster().GetSymmetryImposedCoordinates(Coords,true);
	      Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,true,true);
	      
	      /*
		Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,0,1);
		Cluster::cluster().UpdateLatticeParams(Coords);
		Coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(Coords,1,1);
	      */ 
	      Coords.PrintGradient("Coords");
	      Cluster::cluster().SetNewCoordinates(Coords);
	      Cluster::cluster().RunJobsAndComputeEnergy();
	      double E2 = Cluster::cluster().GetHMBIEnergy();
	      
	      //Reset Lattice Params
	      printf("Reseting Coords\n");
	      Coords = Original_unique_coords;;
	      Cluster::cluster().UpdateLatticeParams(Original_coords);
	      Cluster::cluster().SetNewCoordinates(Original_coords);
	    
	      
	      
	      if (i < 3) {
		LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*AngToBohr);
	      }
	      else{
		LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*DegreesToRadians);
	      }
	    }
	    else
	      LatticeGradient[i] = 0.0;
	  }	  
	}
      

	//double current = Cluster::cluster().GetUnitCellParameter(params[i]);
	//double new1,new2;
	//if (i < 3) {
	//  printf("Shifting coord %s by -%f\n",params[i].c_str(),delta);
	//  new1 = current - delta;
	//}
	//else {
	//  printf("Shifting coord %s by -%f\n",params[i].c_str(),deltaTheta);
	//  new1 = current - deltaTheta;
	//}
	//Cluster::cluster().SetUnitCellParameter(params[i],new1);
	//Cluster::cluster().SetNewCoordinates(Coords);
	//Cluster::cluster().RunJobsAndComputeEnergy();
	//double E1 = Cluster::cluster().GetHMBIEnergy();
      
	// if (i < 3) {
	//  printf("Shifting coord %s by +%f\n",params[i].c_str(),delta);
	//  new2 = current + delta;
	//}
	//else {
	//  printf("Shifting coord %s by +%f\n",params[i].c_str(),deltaTheta);
	//  new2 = current + deltaTheta;
	//}
	// Update the coordinates & get the energy
	//Cluster::cluster().SetUnitCellParameter(params[i],new2);
	//Cluster::cluster().SetNewCoordinates(Coords);
	//Cluster::cluster().RunJobsAndComputeEnergy();
	//double E2 = Cluster::cluster().GetHMBIEnergy();
	
	//if (i < 3) {
	//  LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*AngToBohr);
	//}
	//else{
	//  LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*DegreesToRadians);
	//}
	
	// Reset the lattice parameter to its original value
	//Cluster::cluster().SetUnitCellParameter(params[i],current);
	
	//Finite Difference for the lattice vectors
	if (Params::Parameters().IsPeriodic() && Params::Parameters().DoLatticeVectorFiniteDifference()) {
  
	  string *VectorEntry;
	  VectorEntry = new string[9];
	  VectorEntry[0] = "x1"; VectorEntry[1] = "y1"; VectorEntry[2] = "z1";
	  VectorEntry[3] = "x2"; VectorEntry[4] = "y2"; VectorEntry[5] = "z2";
	  VectorEntry[6] = "x3"; VectorEntry[7] = "y3"; VectorEntry[8] = "z3";
	  
	  Vector Coords = Original_coords; // used to reset coordinates
	  
	  
	  //counter over the vector
	  int i = 0;
	  
	  for(int v=0;v<3;v++){
	    for(int e=0;e<3;e++){
	      
	      //altering the lattice vectors
	      printf("Shifting coord %s by -%f\n",VectorEntry[i].c_str(),delta);
	      Cluster::cluster().AlterLatticeVectors(i,-delta);
	      Cluster::cluster().RunJobsAndComputeEnergy();
	      double E1 = Cluster::cluster().GetHMBIEnergy();
	      
	      //reset lattice vectors and coordinates using the lattice parameters
	      Cluster::cluster().UpdateLatticeParams(Original_coords);
	      Cluster::cluster().SetNewCoordinates(Original_coords);
	      
	      //altering the lattice vectors
	      printf("Shifting coord %s by +%f\n",VectorEntry[i].c_str(),delta);
	      Cluster::cluster().AlterLatticeVectors(i,delta);
	      Cluster::cluster().RunJobsAndComputeEnergy();
	      double E2 = Cluster::cluster().GetHMBIEnergy(); 
	      
	      //reset lattice vectors and coordinates using the lattice parameters
	      Cluster::cluster().UpdateLatticeParams(Original_coords);
	      Cluster::cluster().SetNewCoordinates(Original_coords);
	      
	      //finding Lattice Gradient
	      LatticeGradient[i] = (E2 - E1) / ((2*delta)*AngToBohr);
	      i++;
	    }
	  }
	}
      }
      // Form the actual gradient, G = (Eminus - Eplus)/(2*delta)
      Grad.Set();
      Grad = Eminus;
      Grad -= Eplus;

      Grad.Scale(1.0/(2*delta*AngToBohr));
      Grad.Print("Finite Difference HMBI Gradient");
      
      LatticeGradient.Print("Finite Difference Lattice Gradient");
      
      //Transform Lattice Vector Gradient into Lattice Parameters
      if (Params::Parameters().IsPeriodic() && Params::Parameters().DoLatticeVectorFiniteDifference()) {
      Vector ParamGradient = Cluster::cluster().LatticeParamsGradFromVectorsGrad(LatticeGradient);
      ParamGradient.Print("Finite Difference Lattice Parameter Gradient");
      }
      
      printf("*** Finished with finite difference.");
      exit(1);
      // End finite difference debug code
    }

    // Optimize geometry, if requested
    if ( Params::Parameters().DoForces() ) {
      // DL-Find is our default optimization library
	
	if (Params::Parameters().UseDLFind()) {
	  // Use DL-FIND for geometry optimization
          int nvarin, nvarin2, nspec;
	  int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();
	  nvarin = 3*UniqueAtoms;
	  nvarin2 = 5*UniqueAtoms;// nframe*nat*3 + nweight + nmass + n_po_scaling
	  nspec = 3*UniqueAtoms; // nspec= nat + nz + 5*ncons + 2*nconn + nat
	  
          if(Params::Parameters().UseFullQMOnly()){
            nvarin = 3*UniqueAtoms+9;
	    nvarin2 = 5*UniqueAtoms+15;
	    nspec = 3*UniqueAtoms+9;
          }
	  else if(Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2){
	    nvarin = 3*UniqueAtoms+6;
	    nvarin2 = 5*UniqueAtoms+10;
	    nspec = 3*UniqueAtoms+6;
          }
	  int master = 1; // work is done on the master node
	  printf("Using DL-FIND for optimization\n");
	  dl_find_(&nvarin, &nvarin2, &nspec, &master);

          if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
            Cluster::cluster().PrintFrequencies();

	  // Final energy printing doesn't work, because final
	  //dlf_put_coords erases energies...  
	  double energy =	Cluster::cluster().GetHMBIEnergy();
	  printf("HMBI Final Energy = %15.9f\n",energy);

	  printf("\n");
          if ( Params::Parameters().PrintLevel() > 0)
	    Cluster::cluster().ComputeDistanceMatrix();

	  // Save a copy of the new geometry in a new input file.
	  FILE *input;
	  string input_file = "new_geom.in";
	  if ((input = fopen(input_file.c_str(),"w"))==NULL) {
	    printf("OptimizeGeometry() : Cannot open file '%s'\n",input_file.c_str());
	    exit(1);
          }
	
	  Cluster::cluster().PrintInputFile(input);
	  printf("\nNew input file written to '%s'\n",input_file.c_str());
	  fclose(input);

        }

#ifdef KNITRO
      // Can also use KNITRO.  Experimental at this point.
      else if (Params::Parameters().UseKNitro()) {

	Cluster::cluster().UpdateTrajectoryFile(0,true);

	printf("Using K-Nitro optimizer\n");
	// Define all the parameters Knitro needs
	KTR_context *kc;
	int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
	int n = 3*Natoms; // degrees of freedom
	int objGoal = 0; // minimize
	int objType = 0; // non-linear function
	// X variables are unbounded
	double *xLoBnds, *xUpBnds;
	xLoBnds = new double[n];
	xUpBnds = new double[n];
	for (int i=0;i<n;i++) {
	  xLoBnds[i] = -KTR_INFBOUND;
	  xUpBnds[i] = KTR_INFBOUND;
	}
	
	// Constraints
	int m = 0;
	int *cType = NULL;
	double *cLoBnds = NULL;
	double *cUpBnds = NULL;
	// Jacobian details
	int nnzJ = 0;
	int *jacIndexVars = NULL;
	int *jacIndexCons = NULL;
	
	// Hessian details
	int nnzH = 0;
	int *hessRows = NULL;
	int *hessCols = NULL;
	
	double* lambdaInitial = NULL;

	// Grab initial geometry and energy
	double energy =	Cluster::cluster().GetHMBIEnergy();
	Vector tmp_coords = Cluster::cluster().GetCurrentCoordinates();
	// convert positions to Bohr
	tmp_coords.Scale(AngToBohr);
	double *x;
	x = new double[n];
	for (int i=0;i<3*Natoms;i++) {
	  x[i] = tmp_coords[i];
	}
	
	// Initialize the solver
	kc = KTR_new();
	if (kc == NULL)
	  {
	    printf ("Failed to find a Ziena license.\n");
	    return( -1 );
	  }
	
	// Attach the callback functions to this solver
	//    this tells the solver how to get energies/gradients
	if (KTR_set_func_callback (kc, &knitro_callback_energy_and_gradient) != 0)
	  exit( -1 );
	if (KTR_set_grad_callback (kc, &knitro_callback_energy_and_gradient) != 0)
	  exit( -1 );
	if (KTR_set_hess_callback (kc, &knitro_callback_get_hessian) != 0)
	  exit( -1 );
	//    if (KTR_set_newpoint_callback (kc, &callbackNewPoint) != 0)
	//    exit( -1 );
	
	// Initialize everything
	int nStatus;

	// Set custom options
	nStatus = KTR_set_int_param(kc, KTR_PARAM_ALG, 2); // 0=auto-choose algorithm
	nStatus = KTR_set_int_param(kc, KTR_PARAM_BLASOPTION, 0); // MKL library
	nStatus = KTR_set_int_param(kc, KTR_PARAM_DEBUG, 0); // Debug printing
	nStatus = KTR_set_int_param(kc, KTR_PARAM_HESSOPT, 6); // 6 = L-BFGS hessian
	// Play with the hessian size to see what works well:
	nStatus = KTR_set_int_param(kc, KTR_PARAM_LMSIZE, 10); // L-BFGS hessian size
	nStatus = KTR_set_int_param(kc, KTR_PARAM_MAXIT,Params::Parameters().GetMaxOptCycles());// Max Iterations

	// Convergence criterion:
	nStatus = KTR_set_double_param(kc, KTR_PARAM_OPTTOL, 1.0e-6);
	nStatus = KTR_set_double_param(kc, KTR_PARAM_FEASTOL, 1.0e-6);
	nStatus = KTR_set_double_param(kc, KTR_PARAM_OPTTOLABS, 0.0e0); // not sure
	nStatus = KTR_set_double_param(kc, KTR_PARAM_FEASTOLABS, 0.0e0); // not sure

	// Output
	nStatus = KTR_set_int_param(kc, KTR_PARAM_OUTLEV,6);// Print level 0-6
	nStatus = KTR_set_int_param(kc, KTR_PARAM_OUTMODE,2); // print to both output & logfile

	// For ALG=1:
	nStatus = KTR_set_int_param(kc, KTR_PARAM_BAR_MAXBACKTRACK, 10); // default = 3
	nStatus = KTR_set_int_param(kc, KTR_PARAM_BAR_MAXREFACTOR, 5); // default = 3
	nStatus = KTR_set_int_param(kc, KTR_PARAM_BAR_MURULE, 0); // default = 0 = auto -> 4
	nStatus = KTR_set_int_param(kc, KTR_PARAM_BAR_DIRECTINTERVAL, 0); // 0
	

	nStatus = KTR_init_problem( kc, n, objGoal, objType, 
					xLoBnds, xUpBnds,
					m, cType, cLoBnds, cUpBnds,
					nnzJ, jacIndexVars, jacIndexCons,
					nnzH, hessRows, hessCols, x, lambdaInitial);
	
	
	// Knitro keeps problem specs, so we can free up memory.
	delete [] xLoBnds;
	delete [] xUpBnds;
	
	double* lambda;
	lambda = new double[n+m];
	
	// Call the solver
	nStatus = KTR_solve(kc, x, lambda, 0, &energy, NULL, NULL, NULL, 
			     NULL, NULL, NULL);
	
	
	// Clean up
	KTR_free(&kc);

	delete [] x;
	delete [] lambda;

      }
#endif
      else if (Params::Parameters().GetUseLocalOptimizer()){
        string type;
        if (Params::Parameters().UseConjugateGradient()) {
	  printf("\nUsing Conjugate Gradient for optimization\n");
	  type = "ConjugateGradients";
        }
        else if (Params::Parameters().UseSteepestDescent()) {
	  printf("\nUsing Steepest Descent for optimization\n");
	  type = "SteepestDescent";
        }
        else {
          printf("\nUsing L-BFGS for optimization\n");
	  type = "LBFGS";
        }
       Opt().OptimizeGeometry(type);
      }
    }
    //Performing Hessian Calculation after optimization or if force calculations were not available
    Params::Parameters().SetParameter("IS_OPT_COMPLETE", "1");
    if (Params::Parameters().DoFreqAfterOpt() ) {
      Params::Parameters().SetParameter("JOBTYPE","HESSIAN"); //this turns Do_Forces off and Do_Freq on
    }

    if (  Params::Parameters().DoFreqAfterOpt() &&
         Params::Parameters().IsOptOver()  && Params::Parameters().DoFreq()
	 && !(Params::Parameters().DoFiniteDifferenceFreqs() )
	 && !(Params::Parameters().DoQuasiHarmonic())) {
      string qc_rem = Cluster::cluster().ReadQChemRemSection(infile);
      Params::Parameters().SetParameter("QC_REM",qc_rem);
      
      if (Params::Parameters().AreHessiansAvailable() ) {
	Params::Parameters().SetParameter("ANALYZE_ONLY","TRUE");
      }
      else {
	Params::Parameters().SetParameter("ANALYZE_ONLY","FALSE");
      }
      Cluster::cluster().RunJobsAndComputeEnergy();

      // maybe here we initialize the supercell object
      // remember that supercell currently only features for running a tinker supercell job
      if (Params::Parameters().IsSupercellJob() ) {
        if(!supercell_Init){
	  Supercell::supercell().Initialize(infile,nproc,filename);
          supercell_Init = true;
	  printf("Supercell initialization complete\n");
        }
	
	fflush(stdout);
	
	if (!Params::Parameters().Supercell_Analyze_Only() &&
	    !Params::Parameters().NeglectManyBody()) {
	  printf("doing Supercell tinker job\n");
	  //Supercell::supercell().RunHMBISupercellTinkerFDHessianJob();
	  Supercell::supercell().RunHMBISupercellTinkerHessianJob();
	}
        //printf("\ntrying to set supercell hessian (will probably break here)\n");
        //fflush(stdout);
	//Supercell::supercell().SetSupercellHessian();//really don't need two calls to this
        printf("\ntrying to compute hmbi supercell hessian\n");
        fflush(stdout);
	Supercell::supercell().ComputeHMBISupercellHessian();
        printf("\ntrying to compute phonons\n");
        fflush(stdout);
	Supercell::supercell().ComputePhonons();
	if (Params::Parameters().GetThermalPropertiesUsingPhonons() ) {
	  if(Params::Parameters().GetThermalPropertiesAtMultipleTemperatures()){
	    Supercell::supercell().ComputeThermalPropertiesUsingPhonons(infile);
	  }else{
	    Supercell::supercell().ComputeVibrationalEnergy();
	    double Energy_HMBI = Cluster::cluster().ComputeHMBIEnergy();
	    printf("HMBI Energy = %15.9f hartrees\n",Energy_HMBI);


	    //Vibrational Contribution
	    if( (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()) && 
		(Params::Parameters().DoFreq() ||
		 (Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) ){
	      double Temperature = Params::Parameters().GetTemperature();
	      double Entropy = Cluster::cluster().GetEntropy();
	      double gruneisen_parameter = Cluster::cluster().GetGruneisenParameter();
	      double Cv = Cluster::cluster().GetCv();
	      double Energy_Vib = Cluster::cluster().GetHelmholtzEnergy();
	      printf("\n\n");
	      printf("------Thermodynamic Properaties------------------------\n");
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
	      printf("-----------------------------------------------------");
	    }
	    

	  }
	}
      } 
    }
    
    //Finite difference Hessian
    if ( Params::Parameters().GetJobTypeStr() == "hessian" &&
	 Params::Parameters().DoFiniteDifferenceFreqs()){
      printf("Computing finite difference Hessian\n");  fflush(stdout);


      Params::Parameters().SetParameter("JOBTYPE","FORCE");

      Vector Original_coords = Cluster::cluster().GetCurrentCoordinates();
      int imon = Params::Parameters().GetFrequencyForMonomer();
      printf("Computing the vibrational frequencies for Monomer %d\n",imon);
      Matrix Hessian = GetFiniteDifferenceHessian(Original_coords,imon);
      Hessian = Cluster::cluster().ComputeMassWeightedHessian(Hessian);
      Cluster::cluster().ComputeHarmonicFrequencies(Hessian);

    }
    
    infile.close();

    
    if (Params::Parameters().CheckWarnings() > 0) {
      printf("\nHMBI job completed, but with %d warning(s).\n",Params::Parameters().CheckWarnings());
    }
    else {
      printf("\nHMBI job successful!!!\n");
    } 

#ifdef PARALLEL    
    // Tell all nodes to shut down
    if (nproc > 1) {
      MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
      for (int rank = 1; rank < totalnodes; rank++) {
	MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
      }
    }
#endif /* PARALLEL */

    // Stop the timer and print out the time
    stop_time = time(NULL);
    double elapsed_time = difftime(stop_time,start_time);
    printf("Total job wall time = %.0f seconds\n",elapsed_time);

#ifdef PARALLEL
  }
  
  else {
    /* Slave node controller: just accept & run the jobs */

    int MPI_PrintLevel = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode      

    while (1) {
      char job[BUFFSIZE];
      int success;
      MPI_Status status;     
      while (1) {
	success = 0;
	// run the job
	MPI_Recv(&job,BUFFSIZE,MPI_CHAR,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	
	/* Check the tag of the received message. */
	if (status.MPI_TAG == DIETAG) {
	  printf("%d: Received DIETAG\n",mynode);
	  break;
	}
	else if (status.MPI_TAG == QCHEM_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: Q-chem job: %s\n",mynode,job);
	}
        else if (status.MPI_TAG == QCHEMCP_TAG) {
          if (MPI_PrintLevel > 0 )
            printf("%d: Q-chem job: %s\n",mynode,job);
        }
	else if (status.MPI_TAG == MOLPRO_TAG) {
          if (MPI_PrintLevel > 0 )
            printf("%d: MolPro job: %s\n",mynode,job);	  
	}
	else if (status.MPI_TAG == TINKER_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: Tinker job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == CAMCASP_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: CamCasp job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == G09_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: G09 job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == QUANTUM_ESPRESSO_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: Quantum Espresso job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == PSI4_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: PSI4 job: %s\n",mynode,job);
        }
	else if (status.MPI_TAG == DFTB_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: DFTB+ job: %s\n",mynode,job);
	}
	else {
	  printf("Unknown MPI tag\n");
	}

	// Go ahead and run the job 
	system(job);
	string jobstr = job;
	bool IsCPDimer = CheckIfCPDimerJob(jobstr,status.MPI_TAG);
	success = CheckIfJobSuccessful(jobstr,status.MPI_TAG,IsCPDimer);
	//yoni: I believe this part is causing my hmbi jobs to stop while doing a MM job and without ending process
        if (success==0) {
          //rerun the job.one more try.note that we dont use "while" here.if job fails twice, it is upto the user
          // note that at the present time, it is assumed that the jobs need not be created again if they fail the first time
	  cout << "job " << job << " failed, so rerunning it once\n";
	  system(job);       
	  success = CheckIfJobSuccessful(jobstr,status.MPI_TAG,IsCPDimer);
	  cout << "job " << job << " success = " << success << "\n";
	  }  
	//MPI_Send(&success,1,MPI_INT,0,0,MPI_COMM_WORLD);
	MPI_Send(&success,1,MPI_INT,0,status.MPI_TAG,MPI_COMM_WORLD);
      }
      if (MPI_PrintLevel > 0 )
	printf("MPI: Node %d exiting\n",mynode);
      break;
    } 
  }
  // Finalize MPI and exit
  fflush(stdout);

  //Added extra call to delete all nodes (helps clean up better). JLM
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
  for (int rank = 1; rank < totalnodes; rank++) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }
  MPI_Finalize();
#endif /* PARALLEL */
}

Vector GetFiniteDifferenceGradient(Vector Original_coords) {

  int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
  Vector Eplus(3*Natoms), Eminus(3*Natoms), Grad(3*Natoms);
  
  // Do the finite difference
  double delta = 0.001; // step size
  for (int i=0;i<3*Natoms;i++) {
    
    printf("Shifting coord %d by -%f\n",i,delta);
    Vector Coords = Original_coords;
    Coords[i] -= delta;
    // Update the coordinates & get the energy
    Cluster::cluster().SetNewCoordinates(Coords);
    Cluster::cluster().RunJobsAndComputeEnergy();
    Eplus[i] = Cluster::cluster().GetHMBIEnergy();
    
    printf("Shifting coord %d by +%f\n",i,delta);
    Coords[i] += 2.0*delta;
    // Update the coordinates & get the energy
    Cluster::cluster().SetNewCoordinates(Coords);
    Cluster::cluster().RunJobsAndComputeEnergy();
    Eminus[i] = Cluster::cluster().GetHMBIEnergy();    
  }

  // Form the actual gradient in hartrees/bohr
  Grad.Set();
  Grad = Eminus;
  Grad -= Eplus;
  Grad.Scale(1.0/(2*delta*AngToBohr)); 
  //Grad.Print("Finite Difference HMBI Gradient");
  
  return Grad;
}

// Use analytical nuclear gradient to compute the hessian via finite
// difference
Matrix GetFiniteDifferenceHessian(Vector Original_coords, int imon) {

  //.int UniqueAtoms = Cluster::cluster().GetNumberOfUniqueAtoms();

  //Key to determine which symmetrical unique atoms are in the hessian.
  //int HessKey[UniqueAtoms];
  //int NMon = Cluster::cluster().GetNumberOfMonomers();
  //for(int iMon=1;iMon<NMon;iMon++){
  //  if(Cluster::cluster().GetMonomer(iMon)
  //}



  int Natoms, istart;
  if (imon == 0) {
    Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
    istart = 0;
  }
  else {
    Natoms = Cluster::cluster().GetNumberOfAtoms(imon);
    istart = 0;
    for (int q=1;q<imon;q++) {
      istart += 3*Cluster::cluster().GetNumberOfAtoms(q);
    }
  }
  printf("Natoms = %d\n",Natoms);
  //printf("istart = %d\n",istart);

  Matrix Hess(3*Natoms,3*Natoms);
  Hess.Set();
  
  // Do the finite difference
  double delta = 0.001; // step size
  for (int i=istart;i<istart+3*Natoms;i++) {
    
    printf("Shifting coord %d by -%f\n",i,delta);
    Vector Coords = Original_coords;
    Coords[i] -= delta;
    // Update the coordinates & get the gradient
    Cluster::cluster().SetNewCoordinates(Coords);
    Cluster::cluster().RunJobsAndComputeEnergy();
    Vector Gplus = Cluster::cluster().GetHMBIGradient();
    
    Gplus.PrintGradient("Gplus:");

    printf("Shifting coord %d by +%f\n",i,delta);
    Coords[i] += 2.0*delta;
    // Update the coordinates & get the gradient
    Cluster::cluster().SetNewCoordinates(Coords);
    Cluster::cluster().RunJobsAndComputeEnergy();
    Vector Gminus = Cluster::cluster().GetHMBIGradient();    

    Gminus.PrintGradient("Gminus:");

    // Use the two gradients to determine a column of the Hessian
    // H(i,*) = (Gminus - Gplus)/(2*delta)
    Gminus -= Gplus;
    Gminus.Scale(1.0/(2*delta*AngToBohr)); 

    //Gminus.Print("Grad");

    if (imon == 0)
      Hess.SetColumnVector(Gminus,i);
    else {
      for (int j=0;j<3*Natoms;j++) {
	Hess(j,i-istart) = Gminus[j+istart];
      }
    }
  }
  
  // Form the actual hessian in atomic units, hartree/bohr^2

  Hess.PrintHessian("Finite Difference HMBI Hessian");
  
  return Hess;

}

bool CheckIfCPDimerJob(string job, int job_type) {
 
  bool IsCPDimerJob = 0;
  if (job_type == QCHEMCP_TAG) {
    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.       
    int tmp1 = job.find("; qchem");
    int tmp2 = job.find(".out;",tmp1);
    string trunc1 = job.substr(tmp1,tmp2-tmp1+4);
    int tmp3 = trunc1.rfind(" ");   
    string filename = trunc1.substr(tmp3+1);
    string tmp_d = filename.substr(0,1);
    //  fflush(stdout);

    if (tmp_d == "d") {
      IsCPDimerJob = 1;
    }
  }
  else {
    IsCPDimerJob = 0;
  }
  fflush(stdout);
  return IsCPDimerJob;

}

int CheckIfJobSuccessful(string job, int job_type, bool IsCPDimerJob) {

  int success = 0;

  if (job_type == QCHEM_TAG || job_type == QCHEMCP_TAG) {
    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.
    int tmp1 = job.find("; qchem");
    int tmp2 = job.find(".out;",tmp1);
    string trunc1 = job.substr(tmp1,tmp2-tmp1+4);
    int tmp3 = trunc1.rfind(" ");
    string filename = trunc1.substr(tmp3+1);
    // Now make filename version with full path 
    string outfile = job.substr(3,tmp1-3);
    outfile += "/" + filename;
    int tmp4 = outfile.find(".out");
    string inputfile = outfile.substr(0,tmp4) + ".in";

    // Open Q-Chem output file, search for "JobOver = TRUE"
    string line;
    ifstream jobout;
    int tmp_count;

    fflush(stdout);
    if ( (IsCPDimerJob==1) && (job_type == QCHEMCP_TAG) ) {
      tmp_count = 3;   
     }
    else {
      tmp_count = 1;
    }

    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      if ((line=="GJB about to print energy:") ||(line.find("Thank you") != string::npos)) {
        tmp_count += -1;
      } 
    }

    if (tmp_count==0) {
      success = 1;
    }
 
    jobout.close();
    if (!success) {
      printf("Warning: Q-Chem job %s failed\n",filename.c_str());
    }
  }

  else if (job_type == MOLPRO_TAG) {
    //printf("Currently assume MolPro jobs are successful.\n  Need to implement proper test.\n");

    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.

    //printf("job = %s\n",job.c_str());

    int tmp1 = job.find(";molpro");
    int tmp2 = job.find(".inp");
    //printf("tmp1 = %i,tmp2 = %i\n",tmp1,tmp2);

    string filename = job.substr(tmp1+8,tmp2-tmp1-8);
    string outfile = filename + ".out";
      

    //printf("filename = %s\n",filename.c_str());
    //printf("outfile = %s\n",outfile.c_str());

    // Test for " Variable memory released" in molecule output file
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      //printf("line = %s\n",line.c_str());
      //JLM old way of checking
      //string match = line.substr(0,25);
      //if ( match==" Variable memory released" ) {

      if ( line.find("Variable memory released") != -1 ) {
	success = 1;
	//printf("Job is sucessful\n");
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: Molpro job %s failed\n",filename.c_str());
    }



    //success = 1;

  }
   else if (job_type == TINKER_TAG) {
    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.
    fflush(stdout);
    int tmp1 = job.find("analyze");
    int tmp2 = job.find(".out;",tmp1);
    string trunc1 = job.substr(tmp1,tmp2-tmp1+4);
    int tmp3 = trunc1.rfind(" ");
    string filename = trunc1.substr(tmp3+1);
    // May need to replace "_" with ".".  
    int index = filename.find("_");
    if (index != filename.npos) 
      filename.replace(index,1,".");

    // Now make filename version with full path 
    string outfile = job.substr(3,tmp1-5);
    outfile += "/" + filename;

    //printf("filename = |%s|\n",filename.c_str() );
    //printf("path+filename = |%s|\n",outfile.c_str() );

    // Test for "Total Potential Energy" in tinker output file
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      string match = line.substr(0,23);
      if ( match==" Total Potential Energy" ) {
	success = 1;
      }
    }
    jobout.close();

    if (!success) {
      printf("Warning: Tinker job %s failed\n",filename.c_str());
    }

  }
  else if (job_type == CAMCASP_TAG) {
    //printf("Currently assume CamCasp jobs are successful.\n  Need to implement proper test.\n");
    success = 1;

    if (!success) {
      printf("Warning: CamCasp job failed\n");
    }
  }
  else if(job_type == CRYSTAL_TAG){
    //assume crystal job is sucessful
    success = 1;
    
    if(!success){
      printf("Warning: Crystal job failed\n");
    }
  }
  else if(job_type == G09_TAG){

    
    int tmp1 = job.find("; g09");
    int tmp2 = job.find(".com;");
    string filename = job.substr(tmp1+6,tmp2-tmp1-6);
    string outfile = job.substr(3,tmp1-3);
    outfile += "/" + filename + ".log";

    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());

    int CheckRaman = 0;
    while ( !jobout.eof() ) {
      getline(jobout,line);
      if (line.find("raman") != string::npos) {
	CheckRaman++;
      }
      if (line.find("PolarDeriv") != string::npos) {
	CheckRaman++;
      }
      if (line.find("Normal termination") != string::npos) {
        success = 1;
      } else if (line.find("Error termination") != string::npos && CheckRaman > 1) {
	printf("Warning: G09 job %s failed but PolarDeriv is calculated\n",filename.c_str());
	success = 1;
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: G09 job %s failed\n",filename.c_str());
    }
  }
  else if(job_type == QUANTUM_ESPRESSO_TAG){
    //assume crystal job is successful
    //printf("Currently assumes Quantum Espresso jobs are successful.\n  Need to implement proper test.\n");
    printf("Checking if Quantum Espresso job is successful.\n");
    //int tmp1 = job.find("-i ");
    //int tmp2 = job.find(".in");
    //printf("tmp1 = %i,tmp2 = %i\n",tmp1,tmp2);

    //string filename = job.substr(tmp1+3,tmp2-tmp1-3);
    //string outfile = filename + ".out";
    string outfile = "fullQE.out"; 

    //printf("filename = %s\n",filename.c_str());
    printf("outfile = %s\n",outfile.c_str());
    fflush(stdout);

    // Test for " Variable memory released" in molecule output file
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      //printf("line = %s\n",line.c_str());
      //JLM old way of checking
      //string match = line.substr(0,25);
      //if ( match==" Variable memory released" ) {

      if ( line.find("JOB DONE") != -1 ) {
	success = 1;
	printf("Job is sucessful\n");
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: Quantum Espresso job %s failed\n",outfile.c_str());
    }

   // success = 1;
  }
  else if (job_type == PSI4_TAG) {

    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.

    int tmp1 = job.find("; psi4");
    int tmp2 = job.find(".in;");
   // string filename = job.substr(tmp1+7,tmp2-tmp1-7);
    string filename = job.substr(3,tmp1-3);
  //  filename += "mp2d/";
    filename += "/";
    filename += job.substr(tmp1+7,tmp2-tmp1-7);
    string outfile = filename + ".out";

    // Test for "*** Psi4 exiting successfully. Buy a developer a beer! " in molecule output file
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      //printf("line = %s\n",line.c_str());
      string match = line.substr(0,55);
      if ( match=="*** Psi4 exiting successfully. Buy a developer a beer!" ) {
	success = 1;
	//printf("Job is sucessful\n");
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: PSI4 job %s failed\n",filename.c_str());
    }

    //success = 1;

  }
  else if(job_type == DFTB_TAG){
    //assume crystal job is successful
    //printf("Currently assumes DFTB+ jobs are successful.\n  Need to implement proper test.\n");
    printf("Checking if DFTB+ job is successful.\n");
    //int tmp1 = job.find("-i ");
    //int tmp2 = job.find(".in");
    //printf("tmp1 = %i,tmp2 = %i\n",tmp1,tmp2);

    //string filename = job.substr(tmp1+3,tmp2-tmp1-3);
    //string filename=""
    //string outfile = filename + ".out";
    string outfile = "dftb_out.hsd";

    //printf("filename = %s\n",filename.c_str());
    printf("outfile = %s\n",outfile.c_str());
    fflush(stdout);

    // Test for " Variable memory released" in molecule output file
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      //printf("line = %s\n",line.c_str());
      if ( line.find("Total Energy") != -1 ) {
	success = 1;
	//printf("Job is sucessful\n");
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: DFTB+ job %s failed\n",outfile.c_str());
    }
    //success = 1;
  }
  else {
    printf("Error: CheckIfJobSuccessful: Unknown job type: Tag = %d\n",job_type);
    exit(1);
  }

  fflush(stdout);

  return success;
}
