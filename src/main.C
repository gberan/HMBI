#include "main.h"

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
    Cluster::cluster().Initialize(infile,nproc);
    printf("Cluster initialization complete\n");

    if (! Params::Parameters().DoForces())
      Params::Parameters().SetUseDLFind(false);

    // Compute the Energy (and forces, if requested)
    if (!Params::Parameters().UseDLFind()) {
      Cluster::cluster().RunJobsAndComputeEnergy();
    }

    if ( Params::Parameters().PrintLevel() > 0)
      Cluster::cluster().ComputeDistanceMatrix();

    // For debugging: Compute finite difference stress tensor
    bool do_finite_difference_stress_tensor = false;
    if (do_finite_difference_stress_tensor) {
      double delta = 0.001; // step size, in Angstroms

      Vector a1(3), a2(3), a3(3); // original unit cell vectors
      Vector new_a1(3),new_a2(3), new_a3(3); // deformed unit cell vectors
      Matrix Strain(3,3), Stress(3,3); 
      Stress.Set();
      
      // Grab original unit cell vectors
      a1 = Cluster::cluster().GetUnitCellVector(0);
      a2 = Cluster::cluster().GetUnitCellVector(1);
      a3 = Cluster::cluster().GetUnitCellVector(2);
      printf("Original lattice vectors (Ang):\n");
      printf("a1 = (%f, %f, %f)\n",a1[0],a1[1],a1[2]);
      printf("a2 = (%f, %f, %f)\n",a2[0],a2[1],a2[2]);
      printf("a3 = (%f, %f, %f)\n",a3[0],a3[1],a3[2]);

      // Get cell volume, convert to Bohr^3
      double V = Cluster::cluster().GetCellvolume() * AngToBohr * AngToBohr * AngToBohr;

      // Get the atomic coords.  We use these when triggering the
      // reset of the dimer images, even though we don't change the atomic
      // positions at all.
      Vector Coords = Cluster::cluster().GetCurrentCoordinates();

	// Loop over elements of the strain tensor e_ij
	for (int i=0;i<3;i++) {
	  for (int j=0;j<3;j++) {

	    // Set one element of the strain tensor to +delta
	    printf("Setting Strain(%d,%d) to +%f\n",i,j,delta);
	    Strain.Set(); new_a1.Set(); new_a2.Set(); new_a3.Set();
	    Strain(i,j) = delta;
	    Strain.Print("Strain tensor");
	    
	    // Deform the lattice vectors: new_ai = ai + Strain * ai
	    new_a1 = Strain.MatrixTimesVector(a1);
	    new_a1 += a1;
	    new_a2 = Strain.MatrixTimesVector(a2);
	    new_a2 += a2;
	    new_a3 = Strain.MatrixTimesVector(a3);
	    new_a3 += a3;
	    
	    printf("Updated lattice vectors (Ang):\n");
	    printf("a1 = (%f, %f, %f)\n",new_a1[0],new_a1[1],new_a1[2]);
	    printf("a2 = (%f, %f, %f)\n",new_a2[0],new_a2[1],new_a2[2]);
	    printf("a3 = (%f, %f, %f)\n",new_a3[0],new_a3[1],new_a3[2]);

	    // Set the new lattice vectors and get the energy
	    Cluster::cluster().SetUnitCellVectors(new_a1, new_a2, new_a3);
	    Cluster::cluster().SetNewCoordinates(Coords);
	    Cluster::cluster().RunJobsAndComputeEnergy();
	    double E1 = Cluster::cluster().GetHMBIEnergy();

	    // Repeat for strain tensor element -delta
	    printf("Setting Strain(%d,%d) to -%f\n",i,j,delta);
	    Strain.Set(); new_a1.Set(); new_a2.Set(); new_a3.Set();
	    Strain(i,j) = -delta;
	    Strain.Print("Strain tensor");
	    
	    // Deform the lattice vectors: new_ai = ai + Strain * ai
	    new_a1 = Strain.MatrixTimesVector(a1);
	    new_a1 += a1;
	    new_a2 = Strain.MatrixTimesVector(a2);
	    new_a2 += a2;
	    new_a3 = Strain.MatrixTimesVector(a3);
	    new_a3 += a3;
	    
	    printf("Updated lattice vectors (Ang):\n");
	    printf("a1 = (%f, %f, %f)\n",new_a1[0],new_a1[1],new_a1[2]);
	    printf("a2 = (%f, %f, %f)\n",new_a2[0],new_a2[1],new_a2[2]);
	    printf("a3 = (%f, %f, %f)\n",new_a3[0],new_a3[1],new_a3[2]);

	    // Set the new lattice vectors and get the energy
	    Cluster::cluster().SetUnitCellVectors(new_a1, new_a2, new_a3);
	    Cluster::cluster().SetNewCoordinates(Coords);
	    Cluster::cluster().RunJobsAndComputeEnergy();
	    double E2 = Cluster::cluster().GetHMBIEnergy();

	    // Final Stress tensor element stress(i,j) = 1/V dE/dstrain(i,j)
	    Stress(i,j) = (1.0/V)*(E1 - E2)/(2*delta);
	  }
	}

	printf("Finite Difference Stress Tensor:\n");
	for (int i=0;i<3;i++) {
	  printf("%15.9f  %15.9f  %15.9f\n",Stress(i,0), Stress(i,1), Stress(i,2));
	}
    }



  
    // For debugging: Compute finite difference gradients
    bool do_finite_difference_grad = Params::Parameters().UseFiniteDifferenceGradients();
    if (do_finite_difference_grad && Params::Parameters().DoForces()) {

      Vector Grad = Cluster::cluster().GetHMBIGradient();
      Grad.Print("Analytical HMBI Gradient");
      int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
      Vector Eplus(3*Natoms), Eminus(3*Natoms);
      Vector Original_coords = Cluster::cluster().GetCurrentCoordinates();

      Vector LatticeGradient(6);
      LatticeGradient.Set();

      // Do the finite difference
      double delta = 0.001; // step size, in Angstroms
      /*
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
      */
      // Now do finite difference over the lattice parameters, if appropriate
      double deltaTheta = 0.01;
      if (Params::Parameters().IsPeriodic() ) {
	string *params;
	params = new string[6];
	params[0] = "a"; params[1] = "b"; params[2] = "c";
	params[3] = "alpha", params[4] = "beta"; params[5] = "gamma";

	Vector Coords = Original_coords; // we need these to trigger the reset of the
	// dimer images, even if we don't change the atomic positions at all.

	for (int i=0;i<6;i++) {
	  double current = Cluster::cluster().GetUnitCellParameter(params[i]);
	  double new1,new2;
	  if (i < 3) {
	    printf("Shifting coord %s by -%f\n",params[i].c_str(),delta);
	    new1 = current - delta;
	  }
	  else {
	    printf("Shifting coord %s by -%f\n",params[i].c_str(),deltaTheta);
	    new1 = current - deltaTheta;
	  }
	  Cluster::cluster().SetUnitCellParameter(params[i],new1);
	  Cluster::cluster().SetNewCoordinates(Coords);
	  Cluster::cluster().RunJobsAndComputeEnergy();
	  double E1 = Cluster::cluster().GetHMBIEnergy();

	  if (i < 3) {
	    printf("Shifting coord %s by +%f\n",params[i].c_str(),delta);
	    new2 = current + delta;
	  }
	  else {
	    printf("Shifting coord %s by +%f\n",params[i].c_str(),deltaTheta);
	    new2 = current + deltaTheta;
	  }
	  // Update the coordinates & get the energy
	  Cluster::cluster().SetUnitCellParameter(params[i],new2);
	  Cluster::cluster().SetNewCoordinates(Coords);
	  Cluster::cluster().RunJobsAndComputeEnergy();
	  double E2 = Cluster::cluster().GetHMBIEnergy();

	  if (i < 3) {
	    LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*AngToBohr);
	  }
	  else{
	    LatticeGradient[i] = (E2 - E1) / ((new2 - new1)*DegreesToRadians);
	  }

	  // Reset the lattice parameter to its original value
	  Cluster::cluster().SetUnitCellParameter(params[i],current);

	}

      }

      // Form the actual gradient, G = (Eminus - Eplus)/(2*delta)
      Grad.Set();
      Grad = Eminus;
      Grad -= Eplus;
      Grad.Scale(1.0/(2*delta*AngToBohr));

      Grad.Print("Finite Difference HMBI Gradient");

      LatticeGradient.Print("Finite Difference Lattice Parameter Gradient");

      printf("*** Finished with finite difference.");
      exit(1);
    }
    // End finite difference debug code

    
    // Optimize geometry, if requested
    if ( Params::Parameters().DoForces() ) {
      if (Params::Parameters().UseDLFind()) {
	// Use DL-FIND for geometry optimization
	int Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
	int nvarin = 3*Natoms;
	int nvarin2 = 5*Natoms;// nframe*nat*3 + nweight + nmass + n_po_scaling
	int nspec = 2*Natoms; // nspec= nat + nz + 5*ncons + 2*nconn
	int master = 1; // work is done on the master node
	printf("Using DL-FIND for optimization\n");
	dl_find_(&nvarin, &nvarin2, &nspec, &master);

	// Final energy printing doesn't work, because final
	//dlf_put_coords erases energies...  
	double energy =	Cluster::cluster().GetHMBIEnergy();
	printf("HMBI Final Energy = %15.9f\n",energy);

	printf("\n");
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
      else{
	
	string type = "SteepestDescent";
	//string type = "ConjugateGradients";
	OptimizeGeometry(type );
      }
    }



    if ( Params::Parameters().GetJobTypeStr() == "hessian" ) {
      printf("Computing finite difference Hessian\n");  fflush(stdout);
      Vector Original_coords = Cluster::cluster().GetCurrentCoordinates();
      int imon = Params::Parameters().GetFrequencyForMonomer();
      printf("Computing the vibrational frequencies for Monomer %d\n",imon);
      Matrix Hessian = GetFiniteDifferenceHessian(Original_coords,imon);

      Cluster::cluster().ComputeHarmonicFrequencies(Hessian);

	

    }
    
    infile.close();

    
    if (Params::Parameters().CheckWarnings() > 0) {
      printf("\nHMBI job completed, but with %d warning(s).\n",Params::Parameters().CheckWarnings());
    }
    else {
      printf("\nHMBI job succesful!!!\n");
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
	  //printf("%d: Received DIETAG\n",mynode);
	  break;
	}
	else if (status.MPI_TAG == QCHEM_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: Q-chem job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == TINKER_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: Tinker job: %s\n",mynode,job);
	}
	else if (status.MPI_TAG == CAMCASP_TAG) {
	  if (MPI_PrintLevel > 0 )
	    printf("%d: CamCasp job: %s\n",mynode,job);
	}
	else {
	  printf("Unknown MPI tag\n");
	}

	// Go ahead and run the job 
	system(job);
	string jobstr = job;
	success = CheckIfJobSuccessful(jobstr,status.MPI_TAG);
	
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
  MPI_Finalize();
#endif /* PARALLEL */
}

void OptimizeGeometry(string type) {

  // Initialize geometry opt parameters
  int opt_cycle = 1;
  double dE = 1000, Gnorm = 1000; //nonsense initial values
  double econv = 1.0e-6;
  double gconv = 3.0e-4;
  double stepsize = 0.25;
  bool reset; // for resetting CG algorithm

  // Initialize some empty vectors with same dimensions as the Gradient
  Vector StepDir(Cluster::cluster().GetHMBIGradient(),false);
  Vector Grad_current(Cluster::cluster().GetHMBIGradient(),false);
  Vector StepDir_old(Cluster::cluster().GetHMBIGradient(),false);
  Vector Grad_old(Cluster::cluster().GetHMBIGradient(),false);


  if (type == "SteepestDescent") {
    while (opt_cycle <= Params::Parameters().GetMaxOptCycles() && 
	   (fabs(dE) > econv  || fabs(Gnorm) > gconv ) ) {
      
      
      //printf("XXX\nXXX Starting next opt cycle\n");
      
      // before the step
      double E_current = Cluster::cluster().GetHMBIEnergy();
      Grad_current = Cluster::cluster().GetHMBIGradient();
      Vector Coords_current = Cluster::cluster().GetCurrentCoordinates();
      
      //printf("XXX Current Gradient: %12.6f\n",Grad_current.Norm());
      
      // If this is the first cycle, do some special printing
      if (opt_cycle == 1) {
	printf("Cycle %d: Energy = %15.9f   |Grad| = %12.6f\n",
	       0, E_current, Grad_current.Max(true));
	printf("Cycle %d:     dE = %15.9f  |dGrad| = %12.6f\n",
	       0, 0.0, 0.0);
	Cluster::cluster().UpdateTrajectoryFile(0,true);
      }
      
      // use Gradient printer to show coords
      if (Params::Parameters().PrintLevel() > 0 )
	Cluster::cluster().PrintGradient("Original coordinates",Coords_current);
      
      // Create empty coords array with proper size
      Vector Coords_new(Coords_current, false); 
      
      // Take optimization step
      //Grad_current.Scale(-1.0);
	SteepestDescent(Coords_new, Coords_current, Grad_current, stepsize);
      
      // use Gradient printer to show coords
      if (Params::Parameters().PrintLevel() > 0 )
	Cluster::cluster().PrintGradient("New coordinates",Coords_new);
      
      // Update the coordinates in the cluster object
      Cluster::cluster().SetNewCoordinates(Coords_new);
      
      // Create the new jobs, run them, and get the HMBI energy
      Cluster::cluster().RunJobsAndComputeEnergy();
      
      // Print output
      double E = Cluster::cluster().GetHMBIEnergy();
      Vector Grad( Cluster::cluster().GetHMBIGradient() );
      
      //printf("XXX New Gradient: %12.6f\n",Grad.Max(true));
      
      Gnorm = Grad.Max(true);
      printf("Cycle %d: Energy = %15.9f   |Grad| = %10.6f\n",opt_cycle,
	     E,Gnorm);
      dE = E - E_current;
      Vector dGrad(Grad);
      dGrad -= Grad_current;
      double dG = dGrad.Max(true);
      
      printf("Cycle %d:     dE = %15.9f  |dGrad| = %10.6f  step = %8.3f\n",
	     opt_cycle, dE, dG, stepsize);
      
      Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
      
      
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
      
      // Adjust step size
      if (opt_cycle > 1) {
	// if stepped too far, backup and shrink stepsize
	if (dE > 0.0) {
	  printf("Cycle     Back-up.  Decreasing step size for next cycle.\n");
	  stepsize /= 1.5;
	  // Back up
	  Cluster::cluster().SetNewCoordinates(Coords_current);
	  Cluster::cluster().RunJobsAndComputeEnergy();
	  Grad_current = Grad_old;
	}
	else {
	  printf("Cycle    Increasing step size for next cycle.\n");
	  stepsize *= 1.2;
	}
	
	if (stepsize > 2.0)
	  stepsize = 2.0;
      }
      


      //stepsize = 0.5;
      opt_cycle++;
    }
  }
  else if (type == "ConjugateGradients") {
    // Start optimization cycles
    while (opt_cycle <= Params::Parameters().GetMaxOptCycles() && 
	   (fabs(dE) > econv  || fabs(Gnorm) > gconv ) ) {
      
      reset = false;
      
      //printf("XXX\nXXX Starting next opt cycle\n");
      
      // For CG optimizer, need to save previous gradient
      if (opt_cycle > 1) {
	Grad_old = Grad_current;
	StepDir_old = StepDir;
	
	//printf("Backing up Gradient and StepDir\n");
	//printf("XXX |Grad_old| = %12.6f, |StepDir_old| = %12.6f\n",
	//Grad_old.Norm(),StepDir_old.Norm());
      }
      
      // Grab Energy, Gradient, and XYZ coordinates for geometry
      // before the step
      double E_current = Cluster::cluster().GetHMBIEnergy();
      Grad_current = Cluster::cluster().GetHMBIGradient();
      Vector Coords_current = Cluster::cluster().GetCurrentCoordinates();
      
      //printf("XXX Current Gradient: %12.6f\n",Grad_current.Norm());
      
      // If this is the first cycle, do some special printing
      if (opt_cycle == 1) {
	printf("Cycle %d: Energy = %15.9f   |Grad| = %12.6f\n",
	       0, E_current, Grad_current.Max(true));
	printf("Cycle %d:     dE = %15.9f  |dGrad| = %12.6f\n",
	       0, 0.0, 0.0);
	Cluster::cluster().UpdateTrajectoryFile(0,true);
      }
      
      // use Gradient printer to show coords
      if (Params::Parameters().PrintLevel() > 0 )
	Cluster::cluster().PrintGradient("Original coordinates",Coords_current);
      
      // Create empty coords array with proper size
      Vector Coords_new(Coords_current, false); 
      
      if (opt_cycle % 20 == 0) {
	reset = true;
	//if (stepsize > 0.5)
	stepsize = 0.25;
      }
      
      // Take optimization step
      if (opt_cycle==1) {
	// Use steepest descent for first step, sets StepDir
	Grad_current.Scale(-1.0);
	SteepestDescent(Coords_new, Coords_current, Grad_current, stepsize);
	StepDir = Grad_current;
	StepDir.Scale(-1.0);
      }
      else
	ConjugateGradients(Coords_new, StepDir, Coords_current, Grad_current, 
			   Grad_old, StepDir_old,stepsize, reset);
      
      // use Gradient printer to show coords
      if (Params::Parameters().PrintLevel() > 0 )
	Cluster::cluster().PrintGradient("New coordinates",Coords_new);
      
      // Update the coordinates in the cluster object
      Cluster::cluster().SetNewCoordinates(Coords_new);
      
      // Create the new jobs, run them, and get the HMBI energy
      Cluster::cluster().RunJobsAndComputeEnergy();
      
      // Print output
      double E = Cluster::cluster().GetHMBIEnergy();
      Vector Grad( Cluster::cluster().GetHMBIGradient() );
      
      //printf("XXX New Gradient: %12.6f\n",Grad.Max(true));
      
      Gnorm = Grad.Max(true);
      printf("Cycle %d: Energy = %15.9f   |Grad| = %10.6f\n",opt_cycle,
	     E,Gnorm);
      dE = E - E_current;
      Vector dGrad(Grad);
      dGrad -= Grad_current;
      double dG = dGrad.Max(true);
      
      printf("Cycle %d:     dE = %15.9f  |dGrad| = %10.6f  step = %8.3f\n",
	     opt_cycle, dE, dG, stepsize);
      
      Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
      
      
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
      
      // Adjust step size
      if (opt_cycle > 1) {
	// if stepped too far, backup and shrink stepsize
	if (dE > 0.0) {
	  printf("Cycle     Back-up.  Decreasing step size for next cycle.\n");
	  stepsize /= 1.5;
	  // Back up
	  Cluster::cluster().SetNewCoordinates(Coords_current);
	  Cluster::cluster().RunJobsAndComputeEnergy();
	  Grad_current = Grad_old;
	  StepDir = StepDir_old;
	  
	  
	}
	else {
	  printf("Cycle    Increasing step size for next cycle.\n");
	  stepsize *= 1.2;
	}
	
	if (stepsize > 2.0)
	  stepsize = 2.0;
      }
      


      //stepsize = 0.5;
      opt_cycle++;
    }
  }
    
  printf("Cycle %d: Opt completed.  dE = %15.9f, |Grad| = %10.6f\n",opt_cycle-1,
	 dE,Gnorm);
  
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

    if (imon == 0)
      Hess.SetColumnVector(Gminus,i);
    else {
      for (int j=0;j<3*Natoms;j++) {
	Hess(j,i-istart) = Gminus[j+istart];
      }
    }
  }
  
  // Form the actual hessian in atomic units, hartree/bohr^2

  Hess.Print("Finite Difference HMBI Hessian");
  
  return Hess;

}


int CheckIfJobSuccessful(string job, int job_type) {

  int success = 0;

  if (job_type == QCHEM_TAG) {
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

    // Open Q-Chem output file, search for "JobOver = TRUE"
    string line;
    ifstream jobout;
    jobout.open(outfile.c_str());
    assert(jobout.is_open());
    while ( !jobout.eof() ) {
      getline(jobout,line);
      if (line==" JobOver = TRUE") {
	success = 1;
      }
    }
    jobout.close();
    if (!success) {
      printf("Warning: Q-Chem job %s failed\n",filename.c_str());
    }
  }

  else if (job_type == TINKER_TAG) {
    // Determine the output file name from the job string
    // need to do a bunch of pruning of the strings.
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
  else {
    printf("Error: CheckIfJobSuccessful: Unknown job type: Tag = %d\n",job_type);
    exit(1);
  }

  fflush(stdout);

  return success;
}
