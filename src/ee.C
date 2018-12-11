#include "ee.h"




/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////         Charge Embedding             ////////////////////
////////////////////////////////////////////////////////////////////////////////////


// JDH Creates the list of commands for running the jobs and submitts them...
void RunElectrostaticEmbeddingJobs() {
  int Njobs = Cluster::cluster().GetNumberOfMonomers();
  int NMon  = Cluster::cluster().GetNumberOfMonomers();
  // Step 1: Create list of G09 jobs to run
  string *JobList = new string[NMon];

  int ijob = 0;
  for (int i=1; i<=NMon; i++ ) {

    if ( Params::Parameters().GetMMType()==99 ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = Cluster::cluster().GetMonomer(i).RunG09GDMAJob();
      } else {
	printf("ERROR GDMA only compatible with G09\n");
	exit(1);
      }
    } else if ( Params::Parameters().GetMMType()==98 ) { 
      
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = Cluster::cluster().GetMonomer(i).RunG09ChelpGJob(); 
      } else if (  Params::Parameters().GetQMPackage() == "DALTON"  )  {
	JobList[ijob] = Cluster::cluster().GetMonomer(i).RunOrcaChelpGJob(); 
      } else {
	printf("ERROR ChelpG only supported for Dalton (using orca) and G09 - \n");
	exit(1);
      }


    } else if ( Params::Parameters().GetMMType()==97 ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = Cluster::cluster().GetMonomer(i).RunG09HirshfeldJob(); 
      } else {
	printf("ERROR Hirshfeld only compatible with G09\n");
	exit(1);
      }
    }

    ijob++;
  }
  
  // Step 2: Run the G09 Jobs

#ifdef PARALLEL
  int Nproc = Params::Parameters().GetNumberOfProcessors();
  if (Nproc > (NMon+1) ) { // DEBUG
    printf("NOTE: %d processors will be idle during charge embedding calculation\n", Nproc - NMon);
    Nproc = NMon+1;
  } // END DEBUG
  if (Nproc > 1 ) {
    // Run the Parallel version...
    
    printf("The Charge Embedding jobs will be run in parallel on %d slave processors \n", Nproc-1);
    
    // load up MPI data locally
    MPI_Status status;
    int mynode,totalnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode

    ijob = 0;
    int success = 0; // outcome of jobs;
    int tag, rank; // tag = type of job, rank = processor number
    int result; 


    //DEBUG
    printf("DEBUG: totalnodes = %d,  mynode = %d\n", totalnodes, mynode);

    // Seed a job to each processor initially
    //for (rank=1;rank<totalnodes;rank++) {
    for (rank=1;rank<Nproc;rank++) {

      // Set the tag (JDH: FINISH THIS!!!)
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	tag = G09_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	tag = ORCA_TAG;
      } else {
	printf("ERROR in RunElectrostaticEmbeddingJobs() - not prepared to run QM jobs using this package\n");
	exit(1);
      }


      char job[BUFFSIZE]; // 
      sprintf(job,"%s",JobList[ijob].c_str());
      MPI_Send(&job,BUFFSIZE,MPI_CHAR,rank,tag,MPI_COMM_WORLD);
      ijob++;
    }

    int jobs_run = 0;

    // Now continue running jobs until all are done
    while (ijob < NMon) {
      MPI_Recv(&result,1,MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD,&status);
      
      success += result;

      // Set the tag (JDH: FINISH THIS!!!)
      Cluster::cluster().UpdateJobStatus(jobs_run);
      
      jobs_run++;
    
      char job[BUFFSIZE];
      sprintf(job,"%s",JobList[ijob].c_str());
	
      MPI_Send(&job,BUFFSIZE,MPI_CHAR, status.MPI_SOURCE,tag,MPI_COMM_WORLD);
      ijob++;
    }
    
    // Collect all remaining results
    //for (int rank = 1; rank < totalnodeuse; rank++) {
    for (int rank = 1; rank < Nproc; rank++) {
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      Cluster::cluster().UpdateJobStatus(jobs_run);
      jobs_run++;

      success += result;
    }

    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      printf("Finished G09 Jobs for electrostatic embedding\n");
      if (success == NMon)
	printf("All %d G09 jobs completed successfully\n", success);
      else {
	printf("Error: Only %d of %d G09 jobs completed successfully\n",success,Njobs);
	printf("Correct failed jobs by hand...Then re-run \n");
	exit(1);
      }
    } else if ( Params::Parameters().GetQMPackage() == "ORCA" ||Params::Parameters().GetQMPackage() == "DALTON"  ) {
      printf("Finished ORCA Jobs for electrostatic embedding\n");
      if (success == NMon)
	printf("All %d ORCA jobs completed successfully\n", success);
      else {
	printf("Error: Only %d of %d ORCA jobs completed successfully\n",success,Njobs);
	printf("Correct failed jobs by hand...Then re-run \n");
	exit(1);
      }
    }

    


  } else {
#endif // PARALLEL 
    
    // Serial version:

    ijob = 0;

    while (ijob < NMon ) {
      system(JobList[ijob].c_str() );
      ijob++;
    }

#ifdef PARALLEL 
  }
#endif // PARALLEL

  if ( Params::Parameters().GetMMType()==99 ) {
    // Step 3: run the formchk utility on each of the *.chk file
    for ( int i=1;i<=NMon;i++ ){ 
      system(Cluster::cluster().GetMonomer(i).RunFormChk().c_str() );
    }
    
    // Step 4: Create list of GDMA jobs to run
    ijob = 0;
    for (int i=1; i<=NMon; i++ ) {
      JobList[ijob] = Cluster::cluster().GetMonomer(i).RunGDMAJob();
      ijob++;
    }
  
  // Step 5: Run the gdma jobs

#ifdef PARALLEL
    Nproc = Params::Parameters().GetNumberOfProcessors();
    if (Nproc > 1 ) {
      // Run the Parallel version...
      printf("The G09 jobs for GDMA will be run in parallel on %d slave processors \n", Nproc-1);
      
      // load up MPI data locally
      MPI_Status status;
      int mynode,totalnodes;
      MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
      MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode
      
      ijob = 0;
      int success = 0; // outcome of jobs;
      int tag, rank; // tag = type of job, rank = processor number
      int result; 
      
      // Seed a job to each processor initially
      for (rank=1;rank<totalnodes;rank++) {
	
	// Set the tag (JDH: FINISH THIS!!!)
	tag = GDMA_TAG;
	
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
	
	// Set the tag (JDH: FINISH THIS!!!)
	Cluster::cluster().UpdateJobStatus(jobs_run);
	
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
	
	Cluster::cluster().UpdateJobStatus(jobs_run);
	jobs_run++;
	
	success += result;
      }
      
      printf("Done running parallel G09 jobs for charge embedding\n");
      
      if (success == Njobs)
	printf("All %d charge embedding jobs completed successfully\n", success);
      else {
	printf("Error: Only %d of %d charge emb. jobs completed successfully\n",success,Njobs);
	printf("Correct failed jobs by hand...Then run (SOMETHING... JDH: FINISH THIS!!!)\n");
	exit(1);
      }
      
      
    } else {
#endif // PARALLEL 
      
      // Serial version:
      
      ijob = 0;
      
      while (ijob < NMon ) {
	system(JobList[ijob].c_str() );
	ijob++;
      }
      
#ifdef PARALLEL 
    }
#endif // PARALLEL
  } // end GDMA


  // Step 6: Format the .mom file s.t. HMBI can read them...

  // Open up a file in the MM path and then essentially just read over
  // the portions of the .mom file in the QM directory that we just made
  // formating them such that HMBI can read them in.
  if ( Params::Parameters().GetMMType()==99 ) {
    for (int i=1;i<=NMon;i++) {
      Cluster::cluster().GetMonomer(i).ModifyGDMAMomFile();
    }
  }
  

  delete [] JobList;
}
  

// Non-periodic
void AssignMonomerListsForElectrostaticEmbedding(int NMon, Monomer Monomers[] ) {

  double LARGE_NUMBER = 9999999;
  
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  
  // UNIT CELL MONOMERS
  for (int imon=1;imon<=NMon;imon++ ) {
    double separation = LARGE_NUMBER;
    for ( int a=1; a<=NMon_jobs; a++) {
      separation = min( separation, Monomers[a].FindDistance( Monomers[imon]).Element(0) );
    }

    Monomers[imon].SetUseInEmbedding(false);
    if ( separation < Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(true);
    }
  }
    /*
    Monomers[imon].SetUseInEmbedding(true);
    if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(false);
    }
    */
  
  // IMAGE MONOMERS
  //  for (int imon=1;imon<=NMon_images;imon++ ) // {
  //   double separation = LARGE_NUMBER;
  //   for (int a=1; a<=NMon_jobs;a++) {
  //     separation = min( separation, Monomers[a].FindDistance( MonomerImages[imon]).Element(0) );
  //   }
  //   MonomerImages[imon].SetUseInEmbedding(true);
  //   if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
  //     MonomerImages[imon].SetUseInEmbedding(false);
  //   }

  //   if ( Params::Parameters().UseEwald() ) {
  //     MonomerImages[imon].SetUseInEwald(true);
  //     if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) ) ) {
  // 	MonomerImages[imon].SetUseInEwald(false);
  //     }
  //     // MonomerImages[imon].SetUseInEwaldBuffer(true);
  //     // if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) + Params::Parameters().GetEwaldProbeShellBuffer() ) ) {
  //     // 	MonomerImages[imon].SetUseInEwaldBuffer(false);
  //     // }

       
  //   }
    
  // }
  
}



// Periodic 
void AssignMonomerListsForElectrostaticEmbedding(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]) {



  
  double LARGE_NUMBER = 9999999;
  
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  
  // UNIT CELL MONOMERS
  for (int imon=1;imon<=NMon;imon++ ) {
    double separation = LARGE_NUMBER;
    for ( int a=1; a<=NMon_jobs; a++) {
      separation = min( separation, Monomers[a].FindDistance( Monomers[imon]).Element(0) );
    }

    Monomers[imon].SetUseInEmbedding(false);
    if ( separation < Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(true);
    }
    
    /*
    Monomers[imon].SetUseInEmbedding(true);
    if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(false);
    }
    */

    if ( Params::Parameters().UseEwald() ) {
      Monomers[imon].SetUseInEwald(true);
      if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) ) ) {
	Monomers[imon].SetUseInEwald(false);
      }
      //Monomers[imon].SetUseInEwaldBuffer(true);
      // if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) + Params::Parameters().GetEwaldProbeShellBuffer() ) ) {
      // 	Monomers[imon].SetUseInEwaldBuffer(false);
      // }

    }
  }
  
  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++ ) {
    double separation = LARGE_NUMBER;
    for (int a=1; a<=NMon_jobs;a++) {
      separation = min( separation, Monomers[a].FindDistance( MonomerImages[imon]).Element(0) );
    }
    MonomerImages[imon].SetUseInEmbedding(true);
    if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      MonomerImages[imon].SetUseInEmbedding(false);
    }

    if ( Params::Parameters().UseEwald() ) {
      MonomerImages[imon].SetUseInEwald(true);
      if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) ) ) {
	MonomerImages[imon].SetUseInEwald(false);
      }
      // MonomerImages[imon].SetUseInEwaldBuffer(true);
      // if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) + Params::Parameters().GetEwaldProbeShellBuffer() ) ) {
      // 	MonomerImages[imon].SetUseInEwaldBuffer(false);
      // }

       
    }
    
  }
  
}



// Periodic 
void AssignMonomerListsForTwoBodyElectrostaticEmbedding(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]) {



  
  double LARGE_NUMBER = 9999999;
  
  int NMon_jobs;
  //if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
  NMon_jobs = NMon;
  //} else {
  //  NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  //}

  
  // UNIT CELL MONOMERS
  for (int imon=1;imon<=NMon;imon++ ) {
    double separation = LARGE_NUMBER;
    for ( int a=1; a<=NMon_jobs; a++) {
      separation = min( separation, Monomers[a].FindDistance( Monomers[imon]).Element(0) );
    }

    Monomers[imon].SetUseInEmbedding(false);
    if ( separation < Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(true);
    }
    
    /*
    Monomers[imon].SetUseInEmbedding(true);
    if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      Monomers[imon].SetUseInEmbedding(false);
    }
    */

    if ( Params::Parameters().UseEwald() ) {
      Monomers[imon].SetUseInEwald(true);
      if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) ) ) {
	Monomers[imon].SetUseInEwald(false);
      }
      //Monomers[imon].SetUseInEwaldBuffer(true);
      // if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) + Params::Parameters().GetEwaldProbeShellBuffer() ) ) {
      // 	Monomers[imon].SetUseInEwaldBuffer(false);
      // }

    }
  }
  
  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++ ) {
    double separation = LARGE_NUMBER;
    for (int a=1; a<=NMon_jobs;a++) {
      separation = min( separation, Monomers[a].FindDistance( MonomerImages[imon]).Element(0) );
    }
    MonomerImages[imon].SetUseInEmbedding(true);
    if ( separation > Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
      MonomerImages[imon].SetUseInEmbedding(false);
    }

    if ( Params::Parameters().UseEwald() ) {
      MonomerImages[imon].SetUseInEwald(true);
      if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) ) ) {
	MonomerImages[imon].SetUseInEwald(false);
      }
      // MonomerImages[imon].SetUseInEwaldBuffer(true);
      // if ( separation > (max(Params::Parameters().GetClusterCutoff(),Params::Parameters().GetTwoBodyCutoff()) + Params::Parameters().GetEwaldProbeShellBuffer() ) ) {
      // 	MonomerImages[imon].SetUseInEwaldBuffer(false);
      // }

       
    }
    
  }
  
}




void print(const vector<int> & v, const char * msg)
{
    int size = v.size();

    for (int i = 0; i < size; ++i)
        cout << v[i] << " ";

    cout << msg << endl;
}

void SortAaccordingtoB( Vector& A, Vector& B) {

  srand(time(0));
  
  vector<double> a(A.GetLength()), b(B.GetLength());

  if ( A.GetLength() != B.GetLength() ) {
    printf("ERROR trying to sort A according to B, but length(A) != length(B)\n");
    exit(1);
  }
  
  for (int i = 0; i < A.GetLength(); ++i)
    {
      a[i] = A.Element(i); //%rand() % 10;
      b[i] = B.Element(i); //%i;
    }
  
  //  print(a, "<- A");
  //print(b, "<- B");

  //printf("DEBUG: start sort\n"); fflush(stdout);
  sort(a.begin(), a.end(), MyComparator(b));
  //printf("DEBUG: end sort\n"); fflush(stdout);
  
  for (int i = 0; i < A.GetLength(); ++i) {
    B.Element(i) = b[i];
    A.Element(i) = a[i];
  }

  
  //print(b, "<- B (sorted)");

}



//////////////////////////////
//////////////////////////////
//////////////////////////////
//////////////////////////////
////EWALD EMBEDDING///////////
//////////////////////////////
//////////////////////////////
//////////////////////////////
//////////////////////////////


Matrix ComputeEwaldCharges(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Vector *unit_cell) { //, Vector *ewald_itr ) {


  //int mkl_cntr = 0;
  //int max_mkl_cntr = 10;
  
  //ewald_itr[0].Element(0) = ewald_itr[0].Element(0) + 1;
  
  double AngToBohr = 1.889725989;
  double pi =  3.14159265359;
  double eps = pow(8.85418781762,-12); // C^2/m
  double ec = pow(1.60217733,-19); // C per au
  double Na = pow(6.0221367,23); // Avogadro's number
  double MtoAng = pow(1.0,10.00);
  //double perm = 4*pi*eps*1000/(MtoAng*pow(ec,2)*Na)*4.184; // Gives kkcal/mol for electrostatic potential energy with r in Angstrom
  double perm = 0.003011467140768;


  /////////////// REDO ////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // STEP 1: Compute ewald potential 
  ComputeEwaldPotentialAtAtomCenters( NMon,  Monomers, NMon_images, MonomerImages, unit_cell);


  // STEP 2: Determine the number of Probe points.
  int NProbe = 0;
  // MONOMERS
  for (int imon=1;imon<=NMon;imon++ ) {
    if ( Monomers[imon].GetUseInEwald() ) {
      NProbe += Monomers[imon].GetNumberOfAtoms(); 
    }
    
  }
  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++ ) {
    if ( MonomerImages[imon].GetUseInEwald() ) {
      NProbe += MonomerImages[imon].GetNumberOfAtoms(); // Atom-centered test pt.
    }
  }
  // At this point, NProbe = total number of atoms in the QM region (max of 2bd,cluster cutoff)

  // STEP 3: Build Ewald Probe Point Array
  Matrix TmpEwaldSuperCellTestPoints(NProbe,3);
  Vector EwaldPotential(NProbe);
  int itr = 0;

  // UNIT CELL MONOMERS
  for ( int imon=1; imon<=NMon; imon++) {
    //if ( Monomers[imon].GetUseInEwaldBuffer() ) {
    if ( Monomers[imon].GetUseInEwald() ) {
      for ( int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	for ( int x =0;x<3;x++ ) {
	  TmpEwaldSuperCellTestPoints.Element(itr,x) = Monomers[imon].GetAtom(iatom).GetPosition(x);
	}

	EwaldPotential.Element(itr) = Monomers[imon].GetAtom(iatom).GetEwaldPotential();
	
	itr++;	
      }
    }
  }

  // IMAGE MONOMERS
  for ( int imon=1; imon<NMon_images; imon++ ) {
    if ( MonomerImages[imon].GetUseInEwald() ) {
      for ( int iatom=0; iatom < MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	for ( int x =0;x<3;x++ ) {
	  TmpEwaldSuperCellTestPoints.Element(itr,x) = MonomerImages[imon].GetAtom(iatom).GetPosition(x);
	}
	
	EwaldPotential.Element(itr) = Monomers[ MonomerImages[imon].GetReferenceMonomerIndex()  ].GetAtom(iatom).GetEwaldPotential();

	itr++;

	

      }
    }
  }

  Matrix TestPoints = TmpEwaldSuperCellTestPoints; // replace TmpEwaldSuperCellTestPoints with TestPoints throughout code (to do)


  //OPTIONALLY PRINT A COLLECTION OF F ATOMS TO VISUALIZE TEST POINTS:
  if ( Params::Parameters().PrintLevel() > 0 ) {
    FILE *xyz;
    string filename = "probe_points.xyz";
    if ((xyz = fopen(filename.c_str(),"w"))==NULL) {
      printf("EWALD TEST POINTS: Cannot open file '%s'\n",
  	     filename.c_str());
      exit(1);
    }

    fprintf(xyz,"%d\n F Test point illustration\n", NProbe);

    for ( int i=0; i<NProbe; i++) {
      fprintf(xyz,"F \t %f \t %f \t %f\n", TestPoints.Element(i,0), TestPoints.Element(i,1), TestPoints.Element(i,2) );
      

    }
    fclose(xyz);
  }


  // STEP 4: Build Grid of Fit Points
  Matrix SphereGrid = GetDistanceOrderedAtomCenteredGrid(NMon, Monomers, NMon_images, MonomerImages, NProbe, unit_cell);
  

  // STEP 5: Determine the fixed potential
  Vector FixedPotential = CalculateFixedPotentialAtAtomCenter(NMon, Monomers, NMon_images, MonomerImages, TestPoints);


  // STEP 6: Calculate difference between V(ewald) and V(fixed) = Fit Potential
  Vector FitPotential(NProbe);
  for (int i=0;i<NProbe; i++ ) {
    FitPotential.Element(i) = EwaldPotential.Element(i) - FixedPotential.Element(i);
  }

  // STEP 7: Perform least squares fit:
  Vector OptCharges = CalculateOptimiumCharges(TestPoints, SphereGrid, FitPotential);
  

  // STEP 8: Create EwaldCharges Matrix:
  Matrix EwaldCharges(SphereGrid.GetRows(),4);
  for ( int i=0;i<SphereGrid.GetRows();i++) {
    for (int j=0;j<3;j++) {
      EwaldCharges.Element(i,j) = SphereGrid.Element(i,j);
    }
    //EwaldCharges.Element(i,3) = FitPotential.Element(i);
    EwaldCharges.Element(i,3) = OptCharges.Element(i);
  }


  // OPTIONAL CODE START HERE

  // DEBUG: Test the Least Squares Solution
  Vector CheckFitPotential(TestPoints.GetRows());
  Matrix A( TestPoints.GetRows(), SphereGrid.GetRows() );
  for ( int j=0; j<TestPoints.GetRows();j++ ) {
    for ( int i=0;i<SphereGrid.GetRows();i++ ) {
      Vector R(3);
      for (int r=0;r<3;r++) {
	R.Element(r) = TestPoints.Element(j,r) - SphereGrid.Element(i,r);
      }
      double dist = sqrt( R.DotProduct(R) );
      A.Element(j,i) = 1 / ( dist * perm );
    }
  }


  Matrix coeff = A;

  CheckFitPotential = coeff.MatrixTimesVector(OptCharges);

  
  // STEP 9: Assign the charges to monomer objects for 
  // itr = 0;
  // for ( int imon=1;imon<=NMon;imon++) {
  //   if ( Monomers[imon].GetUseInEwald() ) {
  //     for ( int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++ ) {
  // 	Monomers[imon].GetAtom(iatom).SetEwaldFitPotential(EwaldPotential.Element(itr));
  // 	Monomers[imon].GetAtom(iatom).SetFixedPotential(FixedPotential.Element(itr));
  // 	itr++;
	
  //     }
  //   }
  // }
  // for ( int imon=1;imon<=NMon_images;imon++) {
  //   if ( MonomerImages[imon].GetUseInEwald() ) {
  //     for ( int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++ ) {

  // 	MonomerImages[imon].GetAtom(iatom).SetEwaldFitPotential(EwaldPotential.Element(itr));
  // 	MonomerImages[imon].GetAtom(iatom).SetFixedPotential(FixedPotential.Element(itr));
  // 	itr++;

  //     }
  //   }
  // }

  Vector PotFitDiff(NProbe);
  Vector PotFitPercentDiff(NProbe);


  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("Print Detailed Ewald Performance summary\n");
    printf("Vewald \t Vfix \t Vfit \t VOpt  || CHARGE  \n");
  }
  
  for ( int i=0;i<NProbe;i++) {
    PotFitDiff.Element(i) = EwaldPotential.Element(i) - ( CheckFitPotential.Element(i) + FixedPotential.Element(i)  );
    PotFitPercentDiff.Element(i) = abs(PotFitDiff.Element(i) / EwaldPotential.Element(i)) * 100;

    if ( Params::Parameters().PrintLevel() > 1 ) {
      printf("%f \t %f \t %f \t %f \t %f\n", EwaldPotential.Element(i), FixedPotential.Element(i), FitPotential.Element(i), CheckFitPotential.Element(i), OptCharges.Element(i) );
    }

  }

  printf("Ewald Performance:\n");
  printf("\t RMSD SCRMP vs. EWALD Potential: %f \n", PotFitDiff. RMS() );
  printf("\t RMSD SCRMP vs. EWALD %% diff. : %f \n", PotFitPercentDiff. RMS() );
  

  
  return EwaldCharges;
  
}

Matrix GetDistanceOrderedAtomCenteredGrid(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int Ntest, Vector *unit_cell) {//, double percent ) {


  // Get test points such that they lie on atom centers on monomers surrounding the monomers included in the MonomerImages[] file:
  
  // NOTE: we need to go out far enough such that we have TestPoints.GetRows() atoms included, 5A should be more than enough, but lets put in a loop just to make sure:
  double shell_size = 5;
  //  bool fin = false;
  //  while( fin == false ) {

  // Get test points such that they lie on atom centers on monomers surrounding the monomers included in the MonomerImages[] file:

  double INNER_CUTOFF = max(
			  Params::Parameters().GetTwoBodyCutoff(),
			  Params::Parameters().GetElectrostaticEmbeddingCutoff()
			  );

  double OUTER_CUTOFF = INNER_CUTOFF + shell_size;

 

  // Loop through imagedimer get get all dimers outside emb. cutoff and
  // inside emb.cutoff + NShellFit

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  
  // This is copied from Cluster::CreatePeriodicImageMonomerList() (almost)
  ////////////////////////////////////////////////////////////////////////
  
  // Identify how far we have to go along each unit cell direction to 
  // stay within the cutoff.  Add 1 extra image cell to each, for good measure.
  //double r_cutoff = Params::Parameters().GetMaxPolarizationRadius(); 
  int Nv[3] = {1,1,1}; // start with 1 image cell in each direction
  for (int i=0;i<3;i++) {
    double dist = 0;
    while (dist  < OUTER_CUTOFF) {
      dist +=  unit_cell[i].Norm();
      Nv[i] += 1;
    }
    //printf("Nv[%d] = %d, dist = %f\n",i,Nv[i],dist);
  }
  
  int active_atoms=0;
  int active_non_H_atoms=0;
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
	    for (int jmon=1;jmon<=NMon_jobs;jmon++) {

	      Dimer Tmp;
	      //printf("jmon = %d, imon = %d\n",jmon,imon);
	      Tmp.Initialize(Monomers[jmon],ImageMonomers[imon]);
	      //printf("dimer from %d and image %d\n",jmon,imon);
	      //Tmp.PrintQChemCartesian();
	      if (  Tmp.GetDimerSeparation() < OUTER_CUTOFF ) {
		KeepThisImageMonomer = true;
	      }

	      if ( Tmp.GetDimerSeparation() < INNER_CUTOFF  ) {
		KeepThisImageMonomer = false;
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
	      active_atoms += ImageMonomers[imon].GetNumberOfAtoms();

	      for ( int iatom=0; iatom < ImageMonomers[imon].GetNumberOfAtoms(); iatom++) {
		if ( ImageMonomers[imon].GetAtom(iatom).GetSymbol() != "H" ) {
		  active_non_H_atoms++;
		}
	      }
	      
	    }
	  }
	  delete [] ImageMonomers;
	}	  
      }
  //printf("Total number of image monomers = %d\n",kept_it);

  // if ( active_atoms < TestPoints.GetRows() ) {
  //   printf("WARNING ee.C::increse EWALD_FIT_N_SHELLS not enough fit points included in current cutoff\n");
  // }


  bool use_heavy_atom_only = true;
  int NFitPoints;
  if ( use_heavy_atom_only ) {
    NFitPoints = active_non_H_atoms;
  } else {
    NFitPoints = active_atoms;
  }
  
  
  
  Monomer *EmbeddingMonomers;
  EmbeddingMonomers = new Monomer[kept_it+1];
  //Matrix ShellGrid(NFitPoints,3);
  int ichrg = 0;
  
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

    
    // Renumber it & set reference monomer index
    ImageMon.SetReferenceMonomerIndex(ImageMon.GetIndex());
    ImageMon.SetIndex(Mon_index); 
    
    EmbeddingMonomers[imon] = ImageMon;

    EmbeddingMonomers[imon].SetUnitCellIndex(1,x);
    EmbeddingMonomers[imon].SetUnitCellIndex(2,y);
    EmbeddingMonomers[imon].SetUnitCellIndex(3,z);
    //MonomerImages[imon].SetUnitCellIndexX(x);
    //MonomerImages[imon].SetUnitCellIndexY(y);
    //MonomerImages[imon].SetUnitCellIndexZ(z);

    // for ( int iatom=0; iatom< EmbeddingMonomers[imon].GetNumberOfAtoms(); iatom++) {
    //   if ( ichrg > NFitPoints ) {
    // 	printf("ERROR something is wrong\n"); fflush(stdout);
    //   }
    //   if ( use_heavy_atom_only ) {
    // 	if ( EmbeddingMonomers[imon].GetAtom(iatom).GetSymbol() != "H" ) {
    // 	  for ( int i=0; i<3;i++) {
    // 	    ShellGrid.Element(ichrg,i) = EmbeddingMonomers[imon].GetAtom(iatom).GetPosition(i);
    // 	  }
    // 	  ichrg++;
    // 	}
    //   } else {
    // 	for ( int i=0; i<3;i++) {
    // 	  ShellGrid.Element(ichrg,i) = EmbeddingMonomers[imon].GetAtom(iatom).GetPosition(i);
    // 	}
    // 	ichrg++;
    //   }
    // }

    
  }


  // Now we get a couple vectors of atom lists
  int Natoms_tot = 0;
  for ( imon=1; imon<=kept_it;imon++) {
    for ( int iatom=0; iatom < EmbeddingMonomers[imon].GetNumberOfAtoms(); iatom++ ) {
      if ( use_heavy_atom_only ) {
	if ( EmbeddingMonomers[imon].GetAtom(iatom).GetSymbol() != "H" ) {
	  Natoms_tot++;
	}
      } else {
	Natoms_tot++;
      }
    }
  }

  Matrix AtomPositionList(Natoms_tot,3);
  Vector AtomIndexList(Natoms_tot);
  Vector AtomDistanceList(Natoms_tot);
  double LARGE_NUMBER = 9999999;

  int itr = 0;
  for ( imon=1; imon<=kept_it;imon++) {
    for ( int iatom=0; iatom < EmbeddingMonomers[imon].GetNumberOfAtoms(); iatom++ ) {
      if ( use_heavy_atom_only ) {
	if ( EmbeddingMonomers[imon].GetAtom(iatom).GetSymbol() != "H" ) {
	  AtomPositionList.Element(itr,0) = EmbeddingMonomers[imon].GetAtom(iatom).GetPosition(0);
	  AtomPositionList.Element(itr,1) = EmbeddingMonomers[imon].GetAtom(iatom).GetPosition(1);
	  AtomPositionList.Element(itr,2) = EmbeddingMonomers[imon].GetAtom(iatom).GetPosition(2);
	  AtomIndexList.Element(itr) = itr;

	  // Get distance to atom:
	  double separation = LARGE_NUMBER;
	  for ( int a=1; a<=NMon_jobs; a++) {
	    for ( int b=0; b<Monomers[a].GetNumberOfAtoms(); b++) {
	      separation = min(separation, Monomers[a].GetAtom(b).GetInterAtomicDistance( EmbeddingMonomers[imon].GetAtom(iatom)) );
	    }
	    //separation = min( separation, Monomers[a].FindDistance( Monomers[imon]).Element(0) );
	  }
	  
	  AtomDistanceList.Element(itr) = separation;
	  itr++;
	}
      } else {
	// Don't worry about this, the exact position doesn't really matter and only placing charges on heavy atoms surrounding the emb. cutoff gives a more uniform distribution
      }
    }
  }
  
  SortAaccordingtoB(AtomIndexList, AtomDistanceList);

  // Find the number of atoms for the final grid
  //int Ngrid_pts = 0;
  //for ( int i=0; i<Natoms_tot; i++) {
  //  if ( AtomDistanceList.Element(i) < 
  //}

  //int test = ceil( TestPoints.GetRows() + percent*TestPoints.GetRows());
  int test = Ntest; 
  int Ngrid_pts = min( Natoms_tot, test );

  Matrix ShellGrid(Ngrid_pts,3);
  for ( int i=0; i<Ngrid_pts; i++) {
    ShellGrid.Element(i,0) = AtomPositionList.Element(AtomIndexList.Element(i) ,0 );
    ShellGrid.Element(i,1) = AtomPositionList.Element(AtomIndexList.Element(i) ,1 );
    ShellGrid.Element(i,2) = AtomPositionList.Element(AtomIndexList.Element(i) ,2 );
  }

  
  // Now lets print out all of the charges so that we can see them
  if ( Params::Parameters().PrintLevel() > 0 ) {
    FILE *xyz;
    string filename = "ewald_shell.xyz";
    if ((xyz = fopen(filename.c_str(),"w"))==NULL) {
      printf("EE::()  : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }

    fprintf(xyz,"%d\n", ShellGrid.GetRows() );
    fprintf(xyz,"H atoms placed on every node in the shell grid used to fit Ewald potnential\n");
    for (int i=0;i<ShellGrid.GetRows(); i++) {
      fprintf(xyz,"H \t %f \t %f \t %f\n", ShellGrid.Element(i,0), ShellGrid.Element(i,1), ShellGrid.Element(i,2)  );
    }
    fclose(xyz);

  }

  //  fin = true;

  if ( ShellGrid.GetRows() < Ntest ) {
    printf("ERROR: ee.C -> GetDistanceOrderedAtomCenteredGrid(): \n \t The shell size of 5A didn't provide enough opt points,\\ increase to 6A and rerun code. (this should never happer)\n");
    exit(1);
  }
  
  
  //} // end while loop fin = flase


         
  return ShellGrid;
  delete [] EmbeddingMonomers;


}

Vector CalculateFixedPotentialAtAtomCenter(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Matrix TestPoints) {

  double perm = 0.003011467140768; //f to give kkcal/mol/q (q in au, R in Angstrom)\n", perm);
  Vector FixedPotential(TestPoints.GetRows());

  
  for ( int itest=0;itest<TestPoints.GetRows();itest++) {
    // ADD Contribution from INSIDE unit cell:
    for ( int imon=1;imon<=NMon;imon++ ) {
      if ( Monomers[imon].GetUseInEmbedding() ) {
	for (int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	  Vector R(3);
	  for ( int r=0;r<3;r++) {
	    R.Element(r) = TestPoints.Element(itest,r) - Monomers[imon].GetAtom(iatom).GetPosition(r);
	  }
	  double dist = sqrt( R.DotProduct(R) );
	  if ( dist > 0.00 ) {
	    FixedPotential.Element(itest) += Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0) / (dist * perm);
	  } else {
	  }
	}
      }
    }

    // ADD Contribution from OUTSIDE unit unit cell:
    for ( int imon=1;imon<=NMon_images;imon++ ) {
      if ( MonomerImages[imon].GetUseInEmbedding() ) {
	for (int iatom=0; iatom < MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	  Vector R(3);
	  for ( int r=0;r<3;r++) {
	    R.Element(r) = TestPoints.Element(itest,r) - MonomerImages[imon].GetAtom(iatom).GetPosition(r);
	  }
	  double dist = sqrt( R.DotProduct(R) );
	  if ( dist > 0.00 ) {
	    FixedPotential.Element(itest) += MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0) / (dist * perm);
	  } else {
	    //printf("DEBUG: FIXED pot. evaluated at charge %d\n", itest );
	  }
	}
      }
    }

  }

  return FixedPotential;
}


Vector CalculateOptimiumCharges(Matrix TestPoints, Matrix SphereGrid, Vector FitPotential ){
  double perm = 0.003011467140768;//printf("DEBUG: value of perm = %f to give kkcal/mol/q (q in au, R in Angstrom)\n", perm);

  
  
  Matrix A( TestPoints.GetRows(), SphereGrid.GetRows() );
  for ( int j=0; j<TestPoints.GetRows();j++ ) {
    for ( int i=0;i<SphereGrid.GetRows();i++ ) {
      Vector R(3);
      for (int r=0;r<3;r++) {
	R.Element(r) = TestPoints.Element(j,r) - SphereGrid.Element(i,r);
      }
      double dist = sqrt( R.DotProduct(R) );
      A.Element(j,i) = 1 / ( dist * perm );
    }
  }


  Matrix coeff = A;

  Vector OptCharges(max(TestPoints.GetRows(),SphereGrid.GetRows()));
  for ( int i=0;i< SphereGrid.GetRows(); i++) {
    OptCharges.Element(i) = FitPotential.Element(i);
  }


  bool pass = false;
  A.SolveLeastSquares( OptCharges, &pass);


  if (!pass) {
    printf("\t Expert solver failed, no longer modifying NProbe and Nfit...\n With current algorithim this shouldn't happen. Which means this is bad, very bad...\n" );
    exit(1);
  }


  return OptCharges;

  

}


Matrix BuildEwaldTestPointGrid(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]) {

  // DETERMINE THE NUMBER OF TEST POINTS IN UNIT CELL
  int NTest = 0;
  Matrix EwaldUnitCellTestPoints; // matrix contining the coordinates for all test points within the unit cell:

  for ( int imon=1; imon<=NMon; imon++ ) {
    NTest += Monomers[imon].GetNumberOfAtoms();
  }

  EwaldUnitCellTestPoints.Initialize(NTest,3);

  // ADD ATOM CENTERED TEST POINTS
  int itr = 0;
  for ( int imon=1; imon<=NMon; imon++) {
    for ( int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {

      for ( int x =0;x<3;x++ ) {
	EwaldUnitCellTestPoints.Element(itr,x) = Monomers[imon].GetAtom(iatom).GetPosition(x);
      }
      itr++;
     
    }
  }

  return EwaldUnitCellTestPoints; 
  
}




void ComputeEwaldPotentialAtAtomCenters( int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Vector *unit_cell ) {

  double AngToBohr = 1.889725989;
  double pi =  3.14159265359;
  double eps = pow(8.85418781762,-12); // C^2/m
  double ec = pow(1.60217733,-19); // C per au
  double Na = pow(6.0221367,23); // Avogadro's number
  double MtoAng = pow(1.0,10.00);
  //double perm = 4*pi*eps*1000/(MtoAng*pow(ec,2)*Na)*4.184; // Gives kkcal/mol for electrostatic potential energy with r in Angstrom
  double perm = 0.003011467140768; // Not sure why this isn't working... DEBUG
  //printf("DEBUG: value of perm = %f to give kkcal/mol/q (q in au, R in Angstrom)\n", perm);


  Vector *reciprocal_cell;
  reciprocal_cell = new Vector[3];
  for (int i=0;i<3;i++) {
    //unit_cell[i].Initialize(3);
    reciprocal_cell[i].Initialize(3);
  }
  
  // CREATE reciprocal cell vectors;
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

  double CellVol = unit_cell[0][0]*unit_cell[1][1]*unit_cell[2][2];

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
  

  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("\n\n\tComputing Ewald potential at unit cell atom centers\n" );
  }
  //printf("\tCell contains: %d monomers\n", NMon);

 
  Matrix TestPoints = BuildEwaldTestPointGrid(NMon, Monomers, NMon_images, MonomerImages);
  

  //Vector UnitCellx, Vector UnitCelly, Vector UnitCellz, Vector RecipCellx, Vector RecipCelly, Vector RecipCellz ) {

  int NTest = TestPoints.GetRows(); 

    // double alpha = acos(unit_cell[1].DotProduct(unit_cell[2])/(b*c))*RadiansToDegrees;
    // double beta =  acos(unit_cell[0].DotProduct(unit_cell[2])/(a*c))*RadiansToDegrees;
    // double gamma = acos(unit_cell[0].DotProduct(unit_cell[1])/(a*b))*RadiansToDegrees;
    // printf("GetUnitCellParameter(): a = %f, b = %f, c = %f\n",a,b,c);
    // printf("GetUnitCellParameter(): alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);

  //Volume of a column vectors that make triangular matrix is product of the main diagonal -- DEBUG wait, does this always hold?
  double cell_volume = unit_cell[0][0]*unit_cell[1][1]*unit_cell[2][2];
  // e0 in a.u. = 0.07957747
  // Define some useful constants:


  // GREG's:
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  //double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr

  double accuracy_fac = Params::Parameters().GetEwaldAccuracy();
  double kappa_param = Params::Parameters().GetEwaldKappa();

  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("  Ewald summation cutoffs determined automatically.  Accuracy factor = %f\n",accuracy_fac);
  }
  ///////////////

  // Grab the Ewald summation cutoffs: Direct space
  int nX = Params::Parameters().GetDirecLoopX();
  int nY = Params::Parameters().GetDirecLoopY();
  int nZ = Params::Parameters().GetDirecLoopZ();
  
  // Grab the Ewald summation cutoffs: Reciprocal space
  int kX = Params::Parameters().GetRecipLoopX();
  int kY = Params::Parameters().GetRecipLoopY();
  int kZ = Params::Parameters().GetRecipLoopZ(); 

  bool auto_determine_ewald_parameters = true;
  if ( nX>0 || nY>0 || nZ>0 || kX>0 || kY>0 || kZ>0) {
    auto_determine_ewald_parameters = false;
  }

  double a = unit_cell[0].Norm();
  double b = unit_cell[1].Norm();
  double c = unit_cell[2].Norm();
  
  if (auto_determine_ewald_parameters) {
    // Adapted GJB: Algorithm for determining summation cutoffs in Ewald sum.
    // Based on Frenkel & Smit, Ch. 12.1.5

    double accuracy_fac = Params::Parameters().GetEwaldAccuracy();
    kappa_param = Params::Parameters().GetEwaldKappa();

    //printf("  Ewald summation cutoffs determined automatically.  Accuracy factor = %f\n",accuracy_fac);

    if (kappa_param < 0.0) {
    
      // Auto-determine kappa.  Optimal efficiency seems to come from
      // summing over the same number of cells in real (nX,nY,nZ) and
      // reciprocal space (kX,kY,kZ).  One can show that if nX=kX, 
      // kappa = sqrt(pi)/a, where a is the lattice parameter.  Similar
      // equations exist for the Y and Z components.
      
      // Here, we take an average value of the optimal kappa for each of
      // X, Y, and Z.
      kappa_param = sqrt(pi)*(1.0/a + 1.0/b + 1.0/c)/3.0;
      // Store this value.  
      Params::Parameters().SetEwaldKappa(kappa_param);
      if ( Params::Parameters().PrintLevel() >= 1 ) {
	printf("    Using optimal Ewald kappa parameter = %f\n",kappa_param);
      }
    }
    else {
      if ( Params::Parameters().PrintLevel() > 1 ) {
	printf("    Using user-defined Ewald kappa parameter = %f\n",kappa_param);
      }
    }
    
    //printf("     Ewald accuracy = %f\n", accuracy_fac);

    double rc = accuracy_fac/kappa_param;
    nX = (int) ceil(rc/a);
    nY = (int) ceil(rc/b);
    nZ = (int) ceil(rc/c);

    // Determine reciprocal space cutoffs:
    kX = (int) ceil(accuracy_fac*a*kappa_param/pi);
    kY = (int) ceil(accuracy_fac*b*kappa_param/pi);
    kZ = (int) ceil(accuracy_fac*c*kappa_param/pi);

    if ( Params::Parameters().PrintLevel() >= 1 ) {
      printf("    (nX, nY, nZ) = (%d, %d, %d)    (kX, kY, kZ) = (%d, %d, %d)\n",nX,nY,nZ,kX,kY,kZ);
    }

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

  /////////////



  double CellV = cell_volume;
  if ( Params::Parameters().PrintLevel() > 2 ) {
    printf("Volume of cell: %f\n", cell_volume);
  }
  
  Vector Vewald;
  Vewald.Initialize(NTest);
  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("Performing Real Space sum\n"); fflush(stdout);
  }
  // Loop over test points:
  for (int itest = 0; itest < NTest; itest++ ) {
    if ( Params::Parameters().PrintLevel() > 2 ) {
      printf("Computing Real Space contribtuion to Ewald potential at test point %d of %d\n",itest+1, NTest); fflush(stdout);
    }
    // Loop over monomers in unit cell:
    for (int imonA = 1; imonA <= NMon; imonA++) {
      //if ( imonA != GetIndex() ) { 
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
      
      // Loop over atoms on each monomer
      for (int iA=0; iA<NatomsA; iA++) {
	
	// loop over vectors:
	for (int nx = -nX; nx<=nX; nx++){
	  for (int ny = -nY; ny<=nY; ny++){
	    for (int nz = -nZ; nz<=nZ; nz++){

	      Vector R;
	      R.Initialize(3);
	      for (int r=0;r<3;r++) {
		R.Element(r) = TestPoints.Element(itest,r) - (Monomers[imonA].GetAtom(iA).GetPosition(r) + (nx * unit_cell[0][r] ) + (ny * unit_cell[1][r] ) + (nz * unit_cell[2][r]  ) );
		
	      }
	      double dist = sqrt( R.DotProduct(R) );
	      
	      if (nx*nx+ny*ny+nz*nz == 0  && dist == 0) {
		Vewald.Element(itest) -= 2*kappa_param/sqrt(pi) * Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetMoments().Element(0)/perm;


		
	      } else {
		Vewald.Element(itest) +=  Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetMoments().Element(0) * erfc( kappa_param * dist) / (dist * perm);

	      }
	    }
	  }
	}
	
      } // loop over atoms iA loop
    } // loop over monomer A-loop
  } // loop over test points: itest


  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("Performing Reciprocal Space\n");
  }
  Vector Vrecip;
  Vrecip.Initialize(NTest);
  for (int itest=0; itest < NTest; itest++ ) {
    if ( Params::Parameters().PrintLevel() > 2) {
      printf("Computing Reciprocal Space contribtuion to Ewald potential at test point %d of %d\n",itest+1, NTest); fflush(stdout);
    }
    // Loop over monomers in unit cell:
    for (int imonA = 1; imonA <= NMon; imonA++) {
      // Loop over the atoms on each monomer:
      for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
	// loop over reciprocal lattice vectors
	for (int kx = -kX; kx<=kX; kx++){
	  for (int ky = -kY; ky<=kY; ky++){
	    for (int kz = -kZ; kz<=kZ; kz++){ 
	      // Let's try including all of the terms:
	      if ( !(kx==0 && ky==0 && kz==0) ) {
		Vector R,K;
		K.Initialize(3);
		R.Initialize(3);
		for ( int r=0; r<3; r++ ) {
		  K.Element(r) = (kx * reciprocal_cell[0][r] + ky * reciprocal_cell[1][r] + kz * reciprocal_cell[2][r] );
		  R.Element(r) = (TestPoints.Element(itest,r) - Monomers[imonA].GetAtom(iA).GetPosition(r));
		}
		double k2 = K.DotProduct(K);
		// 1/4 in exp
		Vrecip.Element(itest) +=  (4*pi)/(cell_volume * perm) * Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetMoments().Element(0) * cos(K.DotProduct(R)) * exp( -k2/(4 * pow(kappa_param,2) ))  * 1/(k2);
	      }
	    }
	  }
	}
      }
    }
  }


  // Add terms together:
  for ( int i=0; i<NTest; i++ ) {
    Vewald.Element(i) +=  Vrecip.Element(i);
  }

  // Save the Ewald Potential at atom sites to objects in the atom class inside Monomers
  int itr = 0;
  for ( int imon=1; imon<=NMon; imon++) {
    for ( int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
      Monomers[imon].GetAtom(iatom).SetEwaldPotential( Vewald.Element(itr) );
      itr++;
      
    }
  }
  
 

   // Correction for boundary:
  Vector z(3);
  for (int imonA = 1; imonA <= NMon; imonA++) {
    for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
      for ( int x=0;x<3;x++) {
	z.Element(x) = z.Element(x) + Monomers[imonA].GetAtom(iA).GetPosition(x)*Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetMoments().Element(0);
      }
    }
  }

  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("CORRECTION FOR VACUUM BOUNDARY: %f\n", 2*pi/(3 * cell_volume * perm) * z.DotProduct(z) ); fflush(stdout);
  }
 
  //return Vewald;

}


//////////////////////////////
//////////////////////////////
//////////////////////////////
//////////////////////////////
////TWO-BODY EMBEDDING////////
//////////////////////////////
//////////////////////////////
//////////////////////////////
//////////////////////////////





//void DetermineTwoBodyElectrostaticEmbeddingEnvironment(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) // {
//   // Let's construct a code that allows us to look at charge transfer effects at the two-body level

  
//   AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers );

  
  
//   // CREATE CHELPG DIMER JOBS:
//   for ( int idim=1; idim<=NDim; idim++) {
//     //printf("\n\nDIMER %d \n", idim); fflush(stdout);
//     if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
//       if ( Params::Parameters().GetMMType()==97 ) {
// 	Dimers[idim].CreateG09HirshfeldJob(); 
//       } else if  ( Params::Parameters().GetMMType()==98  ) {
//        	Dimers[idim].CreateG09ChelpGJob(); 
//       } 
//     } 
//   }




//   // RUN CHELPG DIMER JOBS:
//   if ( Params::Parameters().RunJobs() ) {
//     for ( int idim=1; idim<=NDim; idim++) {
//       if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
// 	printf("\t Running D(%d,%d)...\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
// 	if ( Params::Parameters().GetMMType()==97 ) {
// 	  system( Dimers[idim].RunG09HirshfeldJob().c_str()  );
// 	} else if  ( Params::Parameters().GetMMType()==98  ) {
// 	  system( Dimers[idim].RunG09ChelpGJob().c_str()  );
// 	}
//       }
//     }
//   }




//   // Read Charges from Dimer Jobs:
//   // DIMER(S)
//   for (int idim=1;idim<=NDim;idim++) {
//     // if ( USE_DIMER ) {
//     if ( Dimers[idim].GetUseInTwoBodyChargeCalculation()  ) {
//       if ( Params::Parameters().GetMMType()==97 ) {
// 	Dimers[idim].ReadHirshfeldCharges();
//       } else if ( Params::Parameters().GetMMType()==98  ) {
// 	Dimers[idim].ReadChelpGCharges();
//       }
//     }
//   }

//   ComputeTwoBodyCharges(NMon, Monomers, NDim, Dimers);
// }


// This is a beta version to see if its even important to include charges:
//void DetermineTwoBodyElectrostaticEmbeddingEnvironment(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NMon_images, Monomer MonomerImages[], int NDim_images, Dimer DimerImages[] ) // {
//   // Let's construct a code that allows us to look at charge transfer effects at the two-body level

//   printf("DEBUG: two-body charge code for periodic system initializing...\n"); 

//   // CREATE CHELPG DIMER JOBS:
//   for ( int idim=1; idim<=NDim; idim++) {
//     //printf("DIMER %d \n", idim);
//     //if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
//     if ( Dimers[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
//       printf("DEBUG XXX: Creating d(%d,%d) \n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() ); fflush(stdout);
//       if ( Params::Parameters().GetMMType()==97 ) {
// 	Dimers[idim].CreateG09HirshfeldJob(NMon, Monomers, NMon_images, MonomerImages );
//       } else if  ( Params::Parameters().GetMMType()==98  ) {
//        	Dimers[idim].CreateG09ChelpGJob(NMon, Monomers, NMon_images, MonomerImages);
//       } 
//     } 
//   }
//   if ( Params::Parameters().IsPeriodic() ) {
//     for ( int idim=1; idim<=NDim_images; idim++ ) {
//       //if ( DimerImages[idim].GetUseInTwoBodyChargeCalculation() ) {
//       if ( DimerImages[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
// 	printf("DEBUG XXX: Creating image d(%d,%d) \n", DimerImages[idim].GetIndexA(), DimerImages[idim].GetIndexB() ); fflush(stdout);
// 	if ( Params::Parameters().GetMMType()==97 ) {
// 	  DimerImages[idim].CreateG09HirshfeldJob(NMon, Monomers, NMon_images, MonomerImages );
// 	} else if  ( Params::Parameters().GetMMType()==98  ) {
// 	  DimerImages[idim].CreateG09ChelpGJob(NMon, Monomers, NMon_images, MonomerImages);
// 	}
	
//       } 
//     }
//   }

//   printf("DEBUG XXX: NDim_images = %d\n", NDim_images );
  

//   // RUN CHELPG DIMER JOBS:
//   if ( Params::Parameters().RunJobs() ) {
//     for ( int idim=1; idim<=NDim; idim++) {
//       //if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
//       if ( Dimers[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
// 	printf("\t Running D(%d,%d)...\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
// 	if ( Params::Parameters().GetMMType()==97 ) {
// 	  system( Dimers[idim].RunG09HirshfeldJob().c_str()  );
// 	} else if  ( Params::Parameters().GetMMType()==98  ) {
// 	  system( Dimers[idim].RunG09ChelpGJob().c_str()  );
// 	}
//       }
//     }
//     if ( Params::Parameters().IsPeriodic() ) {
//       for ( int idim=1; idim<=NDim_images; idim++ ) {
// 	//if ( DimerImages[idim].GetUseInTwoBodyChargeCalculation() ) {
// 	if ( DimerImages[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
// 	  printf("\t Running D(%d,%d)...\n", DimerImages[idim].GetIndexA(), DimerImages[idim].GetIndexB() );
// 	  if ( Params::Parameters().GetMMType()==97 ) {
// 	    system( Dimers[idim].RunG09HirshfeldJob().c_str()  );
// 	  } else if  ( Params::Parameters().GetMMType()==98  ) {
// 	    system( DimerImages[idim].RunG09ChelpGJob().c_str()  );
// 	  }
// 	}
//       }
//     }
//   }



//   // Read Charges from Dimer Jobs:
//   // DIMER(S)
//   for (int idim=1;idim<=NDim;idim++) {
//     // if ( USE_DIMER ) {
//     //if ( Dimers[idim].GetUseInTwoBodyChargeCalculation()  ) {
//     if ( Dimers[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
//       if ( Params::Parameters().GetMMType()==97 ) {
// 	Dimers[idim].ReadHirshfeldCharges();
//       } else if ( Params::Parameters().GetMMType()==98  ) {
// 	Dimers[idim].ReadChelpGCharges();
//       }
      
//     }
//   }
  

//   if ( Params::Parameters().IsPeriodic() ) {
//     for ( int idim=1; idim<=NDim_images; idim++ ) {
//       //if ( DimerImages[idim].GetUseInTwoBodyChargeCalculation()  ) {
//       if ( DimerImages[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
// 	printf("DEBUG: reading in dimer(%d)\n", idim);
// 	if ( Params::Parameters().GetMMType()==97 ) {
// 	  DimerImages[idim].ReadHirshfeldCharges();
// 	} else if ( Params::Parameters().GetMMType()==98  ) {
// 	  DimerImages[idim].ReadChelpGCharges();
// 	}
	
//       }
//     }
//   }
  
  
//   ComputeTwoBodyCharges(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);

//   printf("DEBUG: two-body charge code finalizing...\n"); fflush(stdout);

// }



//BETA ComputeTwoBodyCharges
//void ComputeTwoBodyCharges(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) // {

//   printf("Computing Two-Body Charges\n"); fflush(stdout);

//   int NMon_jobs;
//   // if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
//   //   NMon_jobs = NMon;
//   // } else {
//   //   NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
//   // }

//   // printf("Number of monomer jobs: %d\n", NMon_jobs ); fflush(stdout);


//   // We need to compute two-body charges for EVERY monoemr in the unit cell: NMon_jobs = NMon
//   NMon_jobs = NMon;

  

//   double tmp_charge; // = new double[4]; // apparently this need to be longer
//   // TWO-BODY CONTRIBUTIONS (within unit cell)
//   for ( int idim=1;idim<=NDim;idim++) {
//     if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
//       for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {

// 	printf("Dimer %d:(%d,%d), atom %d\n",idim, Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB(), iatom); fflush(stdout);
	
// 	tmp_charge = Dimers[idim].GetMonomerA().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	tmp_charge -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	tmp_charge += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	
// 	Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);


//       }

//       if ( Dimers[idim].GetMonomerB().GetIndex() <= Params::Parameters().GetNumAsymmetricMonomers() ) {
// 	for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
// 	  tmp_charge = Dimers[idim].GetMonomerB().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	  
// 	  Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);

// 	}
//       } // end monB
      
//     }
//   }

//   //printf("DEBUG XXX: finished adding the charges to monomer objects\n");fflush(stdout);
  
//   // Compute Charges through Two-Body
//   for (int imon=1; imon<=NMon_jobs; imon++) {
//     for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
//       //printf("DEBUG atom %d\n", iatom); fflush(stdout);
 
//       //printf("DEBUG 6\n"); fflush(stdout);
//       tmp_charge = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
//       tmp_charge += Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody();

//       //printf("DEBUG monomer %d, atom %d\n", imon, iatom); fflush(stdout);
      
//       Monomers[imon].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);
//       //printf("DEBUG 7\n"); fflush(stdout);
//       //Multipole Charge(1, chelpg_charges);
//       //Monomers[imon].GetAtom(iatom).SetMultipoleMoments(Charge);
      
//     }
//   }


//   double *chelpg_charges = new double[4]; // apparently this need so
  
//   printf("DEBUG: Warning, setting two-body charges to one-body object to quickly make NMR jobs for the Ewald paper\n");
//   for (int imon=1; imon<=NMon_jobs; imon++) {
//     for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
//       Monomers[imon].GetAtom(iatom).SetEmbeddedCharge( Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody() );

//       chelpg_charges[0] = Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody();
//       Multipole Charge(1, chelpg_charges);
//       Monomers[imon].GetAtom(iatom).SetMultipoleMoments(Charge);
	  
//     }
//   }
  

  

//   // Print the charges:
//   printf("\n Printing the One and Two-Body charges for each atom in the asymmetric unit\n"); fflush(stdout);

//   printf("Monomer \t Atom \t Type \t One-Body \t Two-Body \n"); fflush(stdout);
//   for (int imon=1;imon<= NMon_jobs; imon++ ) {
//     for (int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++ ) {
//       printf(" %d \t\t %d \t %s \t %f \t %f  \n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0),  Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody() ); fflush(stdout);
//       //printf("\n\nMonomer: %d \t Atom: %d \t Type: %s \t %f \t %f\n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );

//     }
//   }

// }




//BETA 
//void ComputeTwoBodyCharges(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]) // {
  
//   printf("Computing Two-Body Charges\n"); fflush(stdout);

//   int NMon_jobs;
//   // if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
//   //   NMon_jobs = NMon;
//   // } else {
//   //   NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
//   // }

//   // printf("Number of monomer jobs: %d\n", NMon_jobs );

//   // We need to compute two-body charges for EVERY monoemr in the unit cell: NMon_jobs = NMon
//   NMon_jobs = NMon;


//   //AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );

//   double tmp_charge; // = new double[4]; // apparently this need to be longer
//   // TWO-BODY CONTRIBUTIONS (within unit cell)
//   for ( int idim=1;idim<=NDim;idim++) {
    
//     //if ( Dimers[idim].GetUseInTwoBodyChargeCalculation() ) {
//     if ( Dimers[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
      
//       for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	
	
// 	tmp_charge = Dimers[idim].GetMonomerA().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	tmp_charge -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	tmp_charge += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	
// 	Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);


//       }

//       if ( Dimers[idim].GetMonomerB().GetIndex() <= NMon_jobs ) {
// 	for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
// 	  tmp_charge = Dimers[idim].GetMonomerB().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	  
// 	  Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);

// 	}
//       } // end monB
      
//     }
//   }

  
//   for (int idim=1; idim<=NDim_images;idim++) {
//     //if ( DimerImages[idim].GetUseInTwoBodyChargeCalculation()  ) {
//     if ( DimerImages[idim].GetDimerSeparation() <= Params::Parameters().GetTwoBodyCutoff() ) {
//       for (int iatom=0; iatom< DimerImages[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
// 	tmp_charge = DimerImages[idim].GetMonomerA().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	tmp_charge -= Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);

// 	//DEBUG: print the 2-body contridution:
// 	// if ( DimerImages[idim].GetIndexA() == 1 ) {
// 	//   if (  iatom == 0 ) {
// 	//     printf("DEBUG: Two-body contribution for D(%d,%d) atom 1 = %f\n", DimerImages[idim].GetIndexA(), DimerImages[idim].GetIndexB(), tmp_charge );
// 	//     printf("\t (two-body charge running total = %f\n",Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody() );
// 	//   }
// 	// }
	
// 	tmp_charge += Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	
// 	Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);

	
	
// 	//Matrix tmp = DimerImages[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
// 	//tmp -= Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
// 	//tmp += Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
// 	//Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
//       }
//       // Check to see if we need to read in the data for monomer B
//       if ( DimerImages[idim].GetMonomerB().GetIndex() <= NMon_jobs ) {
// 	for (int iatom=0; iatom< DimerImages[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
// 	  tmp_charge = DimerImages[idim].GetMonomerB().GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge -= Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
// 	  tmp_charge += Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetEmbeddedChargeTwoBody();

	  
// 	  Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);



// 	  //Matrix tmp = DimerImages[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
// 	  //tmp -= Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
// 	  //tmp += Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
// 	  //Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	  
// 	} 
//       }// end monB
//     }
//   }
  

//   // Compute Charges through Two-Body
//   for (int imon=1; imon<=NMon_jobs; imon++) {
//     for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
//       //printf("DEBUG atom %d\n", iatom); fflush(stdout);
 
//       //printf("DEBUG 6\n"); fflush(stdout);
//       tmp_charge = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0);
//       tmp_charge += Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody();

      
//       Monomers[imon].GetAtom(iatom).SetEmbeddedChargeTwoBody(tmp_charge);

//       //printf("DEBUG 7\n"); fflush(stdout);
//       //Multipole Charge(1, chelpg_charges);
//       //Monomers[imon].GetAtom(iatom).SetMultipoleMoments(Charge);
      
//     }
//   }




  
//   // Print the charges:
//   printf("\n Printing the One and Two-Body charges for each atom in the asymmetric unit\n"); fflush(stdout);

//   printf("Monomer \t Atom \t Type \t One-Body \t Two-Body \n");
//   for (int imon=1;imon<= NMon_jobs; imon++ ) {
//     for (int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++ ) {
//       printf(" %d \t\t %d \t %s \t %f \t %f  \n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0),  Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody() );
//       //printf("\n\nMonomer: %d \t Atom: %d \t Type: %s \t %f \t %f\n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );

//     }
//   }


//   double *chelpg_charges = new double[4]; // apparently this is needed...

//   printf("DEBUG: Warning, setting two-body charges to one-body object to quickly make NMR jobs for the Ewald paper\n");
//   for (int imon=1; imon<=NMon_jobs; imon++) {
//     for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
//       Monomers[imon].GetAtom(iatom).SetEmbeddedCharge( Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody() );


//       chelpg_charges[0] = Monomers[imon].GetAtom(iatom).GetEmbeddedChargeTwoBody();
//       Multipole Charge(1, chelpg_charges);
//       Monomers[imon].GetAtom(iatom).SetMultipoleMoments(Charge);
//     }
//   } 


  

// }
