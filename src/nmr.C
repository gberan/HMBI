#include "nmr.h"


#include <algorithm> // JDH - need this for 'sort' function
#include <vector>

using namespace hmbi_constants;
#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */





//void PredictMagneticProperties() // { UNFNISHED - still in cluster.C

//   printf("\n\n\n");
//   printf("***************************************************\n");
//   printf("*   Fragment-Based Magnetic Propery Prediction    *\n");
//   printf("***************************************************\n\n");



//   int NMon = Cluster::cluster().GetNumberOfMonomers();


//   // Electrostatic Embedding:
//   if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    
//     // BUILD JOBS
//     if ( Params::Parameters().GetMMType()==99 ) {
//       for (int i=1;i<=NMon;i++) {
// 	if (Params::Parameters().RunJobs() || Params::Parameters().KeepNMRFiles() ) {
// 	  Cluster::cluster().GetMonomer(i).CreateGDMAJob();
// 	}
//       }
//     } else if ( Params::Parameters().GetMMType()==98 ) {
//       for (int i=1; i<=NMon;i++) {
// 	if (Params::Parameters().RunJobs() || Params::Parameters().KeepNMRFiles() ) {
// 	  Cluster::cluster().GetMonomer(i).CreateG09ChelpGJob(); 
// 	}
//       }
//     } else {
//       printf("The NMR code only uses GDMA or ChelpG to save a lot of time! \n");
//       exit(1);
//     }

//     // RUN JOBS
//     if ( Params::Parameters().RunJobs() ) {
//       RunElectrostaticEmbeddingJobs();
//     }

//     // READ EE DATA
//     if ( Params::Parameters().GetMMType()==99 ) {
//       // Now read in the Multipole Moment Data:
//       //printf("\nReading multipole moments from file...\n\n"); fflush(stdout);
//       for (int imon=1;imon<=NMon;imon++) {
// 	Cluster::cluster().GetMonomer(imon).ReadMultipoleMoments(); fflush(stdout);
//       }
//     } else if (Params::Parameters().GetMMType()==98 ) {
//       //printf("\nReading ChelpG charges from file...\n\n"); fflush(stdout);
//       for (int imon=1;imon<=NMon;imon++) {
// 	Cluster::cluster().GetMonomer(imon).ReadCHelpGCharges();
//       }
//     }

    
//     // Self-Consistent EE (optional)
//     if ( Params::Parameters().UseSelfConsistentEmbedding() ) {
//       if ( Params::Parameters().RunJobs() ) {
// 	DetermineSelfConsistentElectrostaticEmbeddingEnvironment();
//       }
//     }

 
//     exit(1);

//   } // END Electrostatic Embedding:

// }


// Misc. Functions

void AssignBasisSetsToAtoms(int NMon, Monomer Monomers[]  ) {

  double LARGE_NUMBER = 999999;


  printf("\n\nUsing Locally Dense Basis set:  \n");
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel1().c_str(), Params::Parameters().GetMixedBasisCutOff() );
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel2().c_str(), Params::Parameters().GetMixedBasisCutOff2() );
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel3().c_str(), Params::Parameters().GetMixedBasisCutOff3() );
  
  
  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("\n\nAll atoms in asymmetric unit assigned basis: %s  \n", Params::Parameters().GetMixedBasisLevel1().c_str() );
  }
  for (int imon=1;imon<=Params::Parameters().GetNumAsymmetricMonomers();imon++) {
    for (int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(1);
    }
  }
  
  
  // UNIT CELL MONOMERS
  for (int imon=Params::Parameters().GetNumAsymmetricMonomers();imon<=NMon;imon++) {
    for (int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      double separation = LARGE_NUMBER;
      for (int i=1;i<=Params::Parameters().GetNumAsymmetricMonomers();i++) {
	for (int a=0;a<Monomers[i].GetNumberOfAtoms();a++) {
	  separation = min(separation, Monomers[i].GetAtom(a).GetInterAtomicDistance( Monomers[imon].GetAtom(iatom) ) );
	}
      }
      // Assign basis set based on separation distance
      if ( separation <= Params::Parameters().GetMixedBasisCutOff() ) {
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(1);
      } else if ( separation <= Params::Parameters().GetMixedBasisCutOff2() ) {
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(2);
      } else { 
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(3);
      }
      
    }
  }
  
  
  
}




void AssignBasisSetsToAtoms(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ) {

  double LARGE_NUMBER = 999999;


  printf("\n\nUsing Locally Dense Basis set:  \n");
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel1().c_str(), Params::Parameters().GetMixedBasisCutOff() );
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel2().c_str(), Params::Parameters().GetMixedBasisCutOff2() );
  printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel3().c_str(), Params::Parameters().GetMixedBasisCutOff3() );
  
  
  if ( Params::Parameters().PrintLevel() > 1 ) {
    printf("\n\nAll atoms in asymmetric unit assigned basis: %s  \n", Params::Parameters().GetMixedBasisLevel1().c_str() );
  }
  for (int imon=1;imon<=Params::Parameters().GetNumAsymmetricMonomers();imon++) {
    for (int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(1);
    }
  }
  
  
  // UNIT CELL MONOMERS
  for (int imon=Params::Parameters().GetNumAsymmetricMonomers();imon<=NMon;imon++) {
    for (int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      double separation = LARGE_NUMBER;
      for (int i=1;i<=Params::Parameters().GetNumAsymmetricMonomers();i++) {
	for (int a=0;a<Monomers[i].GetNumberOfAtoms();a++) {
	  separation = min(separation, Monomers[i].GetAtom(a).GetInterAtomicDistance( Monomers[imon].GetAtom(iatom) ) );
	}
      }
      // Assign basis set based on separation distance
      if ( separation <= Params::Parameters().GetMixedBasisCutOff() ) {
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(1);
      } else if ( separation <= Params::Parameters().GetMixedBasisCutOff2() ) {
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(2);
      } else { 
	Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(3);
      }
      
    }
  }
  
  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++) {
    for (int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
      double separation = LARGE_NUMBER;
      for (int i=1;i<=Params::Parameters().GetNumAsymmetricMonomers();i++) {
	for (int a=0;a<Monomers[i].GetNumberOfAtoms();a++) {
	  separation = min(separation, Monomers[i].GetAtom(a).GetInterAtomicDistance( MonomerImages[imon].GetAtom(iatom) ) );
	}
      }
      
      // now, based on the distance, set the atom to a MixedBasisRegion...
      if ( separation <= Params::Parameters().GetMixedBasisCutOff() ) {
	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(1);
      } else if ( separation <= Params::Parameters().GetMixedBasisCutOff2() ) {
	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(2);
      } else { 
	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(3);
      }
    }
  }
  
  
}


//void AssignBasisSetsToAtomsTwoBodyCharge(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ) // {

//   double LARGE_NUMBER = 999999;


//   printf("\n\nUsing Locally Dense Basis set:  \n");
//   printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel1().c_str(), Params::Parameters().GetMixedBasisCutOff() );
//   printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel2().c_str(), Params::Parameters().GetMixedBasisCutOff2() );
//   printf("Basis %s out to %f Anstroms from Aymmetric unit\n", Params::Parameters().GetMixedBasisLevel3().c_str(), Params::Parameters().GetMixedBasisCutOff3() );
  
  
//   if ( Params::Parameters().PrintLevel() > 1 ) {
//     printf("\n\nAll atoms in asymmetric unit assigned basis: %s  \n", Params::Parameters().GetMixedBasisLevel1().c_str() );
//   }
//   for (int imon=1;imon<=NMon;imon++) {
//     for (int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
//       Monomers[imon].GetAtom(iatom).SetMixedBasisRegion(1);
//     }
//   }
  
 
  
//   // IMAGE MONOMERS
//   for (int imon=1;imon<=NMon_images;imon++) {
//     for (int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
//       double separation = LARGE_NUMBER;
//       for (int i=1;i<=NMon;i++) {
// 	for (int a=0;a<Monomers[i].GetNumberOfAtoms();a++) {
// 	  separation = min(separation, Monomers[i].GetAtom(a).GetInterAtomicDistance( MonomerImages[imon].GetAtom(iatom) ) );
// 	}
//       }
      
//       // now, based on the distance, set the atom to a MixedBasisRegion...

      
      
//       if ( separation <= Params::Parameters().GetMixedBasisCutOff() ) {
// 	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(1);
//       } else if ( separation <= Params::Parameters().GetMixedBasisCutOff2() ) {
// 	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(2);
	
//       } else { 
// 	MonomerImages[imon].GetAtom(iatom).SetMixedBasisRegion(3);
	
//       }
//     }
//   }
  
  
// }


void CreateClusterJobs(int NMon, Monomer Monomers[]  )  {

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  printf("\n\n Creating %d Cluster Job(s) with cutoff: %f \n",NMon_jobs, Params::Parameters().GetClusterCutoff());

  if ( Params::Parameters().GetNumAsymmetricMonomers() > NMon ) {
    printf("ERROR in CreateClusterJob() inconsistent number of asymmetric monomers\n");
    exit(1);
  }

  
  // CHECK FOR SPIN COMPATIBILITY (will need to update this ASAP for EPR stuff!!)
  for (int i=1; i<=NMon_jobs;i++) {
    if ( Monomers[i].GetSpinState() != 1 ) {
      printf("ERROR: Cluster::CreateQMNMRJobs() cluster/fragment method \n");
      printf("not ready for monomers with non-singlet spin states\n");
      exit(1);
    }
    if (Params::Parameters().PrintLevel() > 0 ) {
      printf("Including Monomer %d in cluster\n", Monomers[i].GetIndex() ); 
    }
  }
    

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

    // MAKE CLUSTER JOB
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      CreateClusterG09Job(NMon, Monomers, total_charge, iclust);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
      CreateClusterDALTONJob(NMon, Monomers, total_charge, iclust);
    } else {
      printf("ERROR: CreateClusterJob() QM Package %s not yet supported\n", Params::Parameters().GetQMPackage().c_str() );
    }
    
  }
  
}



void CreateClusterJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] )  {

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  printf("\n\n Creating %d Cluster Job(s) with cutoff: %f \n",NMon_jobs, Params::Parameters().GetClusterCutoff());

  if ( Params::Parameters().GetNumAsymmetricMonomers() > NMon ) {
    printf("ERROR in CreateClusterJob() inconsistent number of asymmetric monomers\n");
    exit(1);
  }

  
  // CHECK FOR SPIN COMPATIBILITY (will need to update this ASAP for EPR stuff!!)
  for (int i=1; i<=NMon_jobs;i++) {
    if ( Monomers[i].GetSpinState() != 1 ) {
      printf("ERROR: Cluster::CreateQMNMRJobs() cluster/fragment method \n");
      printf("not ready for monomers with non-singlet spin states\n");
      exit(1);
    }
    if (Params::Parameters().PrintLevel() > 0 ) {
      printf("Including Monomer %d in cluster\n", Monomers[i].GetIndex() ); 
    }
  }
    

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
    for (int imon=1; imon<=NMon_images;imon++) {
      MonomerImages[imon].SetUseInClusterCalculation(false);
 
      if ( Monomers[iclust].FindDistance( MonomerImages[imon] ).Element(0) <= Params::Parameters().GetClusterCutoff() ) {
	MonomerImages[imon].SetUseInClusterCalculation(true);
	total_charge += MonomerImages[imon].GetChargeState();
      }
    }

    // MAKE CLUSTER JOB
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      CreateClusterG09Job(NMon, Monomers, NMon_images, MonomerImages, total_charge, iclust);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
      CreateClusterDALTONJob(NMon, Monomers, NMon_images, MonomerImages, total_charge, iclust);
    } else {
      printf("ERROR: CreateClusterJob() QM Package %s not yet supported\n", Params::Parameters().GetQMPackage().c_str() );
    }
    
  }
  
}

void CreateClusterJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Matrix EwaldCharges ) {

    int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  printf("\n\n Creating %d Cluster Job(s) with cutoff: %f \n",NMon_jobs, Params::Parameters().GetClusterCutoff());

  if ( Params::Parameters().GetNumAsymmetricMonomers() > NMon ) {
    printf("ERROR in CreateClusterJob() inconsistent number of asymmetric monomers\n");
    exit(1);
  }

  
  // CHECK FOR SPIN COMPATIBILITY (will need to update this ASAP for EPR stuff!!)
  for (int i=1; i<=NMon_jobs;i++) {
    if ( Monomers[i].GetSpinState() != 1 ) {
      printf("ERROR: Cluster::CreateQMNMRJobs() cluster/fragment method \n");
      printf("not ready for monomers with non-singlet spin states\n");
      exit(1);
    }
    if (Params::Parameters().PrintLevel() > 0 ) {
      printf("Including Monomer %d in cluster\n", Monomers[i].GetIndex() ); 
    }
  }
    

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
    for (int imon=1; imon<=NMon_images;imon++) {
      MonomerImages[imon].SetUseInClusterCalculation(false);
 
      if ( Monomers[iclust].FindDistance( MonomerImages[imon] ).Element(0) <= Params::Parameters().GetClusterCutoff() ) {
	MonomerImages[imon].SetUseInClusterCalculation(true);
	total_charge += MonomerImages[imon].GetChargeState();
      }
    }

    // MAKE CLUSTER JOB
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      CreateClusterG09Job(NMon, Monomers, NMon_images, MonomerImages, total_charge, EwaldCharges, iclust);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
      CreateClusterDALTONJob(NMon, Monomers, NMon_images, MonomerImages, total_charge, EwaldCharges, iclust);
    } else {
      printf("ERROR: CreateClusterJob() QM Package %s not yet supported\n", Params::Parameters().GetQMPackage().c_str() );
    }
    
  }


}


void CreateClusterG09Job(int NMon, Monomer Monomers[], int total_charge, int iclust ) {
  
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  std::ostringstream ss;
  ss << iclust;
  if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
    //sprintf(filename,"%s/cluster.%d.com",path.c_str(),iclust);
    filename = path + "/" + "cluster." + ss.str() + ".com"; 
  } else {
    filename = path + "/" + "tmp" + ".com"; // HACK
  }
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("CreateClusterG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Cluster fragment composed of all monomomers within %12.2f angstrom of the first %d monomers in the input file\n", Params::Parameters().GetClusterCutoff(), Params::Parameters().GetNumAsymmetricMonomers()  );
  fprintf(job,"\n");
  
  fprintf(job,"%d %d\n", total_charge, 1 ); // NOTE: spin hard-coded to singlet for now...


  // UNIT CELL
  Monomers[iclust].PrintMonomerCartesian(job);
  for (int imon=1;imon<=NMon;imon++ ) {
    if ( Monomers[imon].GetUseInClusterCalculation() && imon != iclust ) {
      Monomers[imon].PrintMonomerCartesian(job);
    }
  }


  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    fprintf(job,"\n"); // use a blank line to separate

    for ( int i=1;i<=NMon; i++) {
      if ( !Monomers[i].GetUseInClusterCalculation() && Monomers[i].GetUseInEmbedding() ) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }

    
    // for (int i=Params::Parameters().GetNumAsymmetricMonomers();i<=NMon;i++) {
    //   if ( !Monomers[i].GetUseInClusterCalculation() && Monomers[iclust].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
    // 	Monomers[i].PrintEmbeddingCharges(job);
    //   }
    // }
    
     
  } // END Params::Parameters().UseElectrostaticEmbedding()



  // FULL CUSTOM BASIS (not compatible with LDBS)
  if ( Params::Parameters().CustomBasis() == 2 ){
    // USE Custom Basis for all atoms in the cluster:
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    // Inside unit cell:
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    printf("WARNING: Custom Basis not specified for atom type: %s, using mixed basis level 1 as default. \n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );
	  }
	  
	  fprintf(job,"****\n");
	  
	}
      }
    } // End loop over unit cell atoms
    
     
  }
  
  
  //  MIXED BASIS 
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Params::Parameters().CustomBasis() == 0 ) {
	    // Now print the basis based on the region...
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	  } else if ( Params::Parameters().CustomBasis() == 1) {
	    //printf("DEBUG: USING CUSTOM BASIS!!!\n");
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      printf("DEBUG: USING CUSTOM BASIS!!!\n");
	      if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	      } else {
		fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	      }
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	    
	  }
	  fprintf(job,"****\n");
	}
      }
    }
  }
  
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job); 
  
}



void CreateClusterG09Job(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, int iclust ) {
  
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  std::ostringstream ss;
  ss << iclust;
  if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
    //sprintf(filename,"%s/cluster.%d.com",path.c_str(),iclust);
    filename = path + "/" + "cluster." + ss.str() + ".com"; 
  } else {
    filename = path + "/" + "tmp" + ".com"; // HACK
  }
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("CreateClusterG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Cluster fragment composed of all monomomers within %12.2f angstrom of the first %d monomers in the input file\n", Params::Parameters().GetClusterCutoff(), Params::Parameters().GetNumAsymmetricMonomers()  );
  fprintf(job,"\n");
  
  fprintf(job,"%d %d\n", total_charge, 1 ); // NOTE: spin hard-coded to singlet for now...


  // UNIT CELL
  Monomers[iclust].PrintMonomerCartesian(job);
  for (int imon=1;imon<=NMon;imon++ ) {
    if ( Monomers[imon].GetUseInClusterCalculation() && imon != iclust ) {
      Monomers[imon].PrintMonomerCartesian(job);
    }
  }

  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++ ) {
    if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
      MonomerImages[imon].PrintMonomerCartesian(job);
    }
  }

  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    fprintf(job,"\n"); // use a blank line to separate


    for ( int i=1;i<=NMon; i++) {
      if ( !Monomers[i].GetUseInClusterCalculation() && Monomers[i].GetUseInEmbedding() ) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }

    
    // for (int i=Params::Parameters().GetNumAsymmetricMonomers();i<=NMon;i++) {
    //   if ( !Monomers[i].GetUseInClusterCalculation() && Monomers[iclust].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
    // 	Monomers[i].PrintEmbeddingCharges(job);
    //   }
    // }

    for (int i=1; i<=NMon_images;i++ ) {
      if ( !MonomerImages[i].GetUseInClusterCalculation() && MonomerImages[i].GetUseInEmbedding() ) {
	MonomerImages[i].PrintEmbeddingCharges(job);
      }
    }
    
    
    // for (int i=1; i<=NMon_images;i++ ) {
    //   if ( MonomerImages[i].GetUseInEmbedding() ) {
    // 	if ( !MonomerImages[i].GetUseInClusterCalculation() && Monomers[iclust].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
    // 	  MonomerImages[i].PrintEmbeddingCharges(job);
    // 	}
    //   }
    // }
      
    
  } // END Params::Parameters().UseElectrostaticEmbedding()



  // FULL CUSTOM BASIS (not compatible with LDBS)
  if ( Params::Parameters().CustomBasis() == 2 ){
    // USE Custom Basis for all atoms in the cluster:
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    // Inside unit cell:
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    printf("WARNING: Custom Basis not specified for atom type: %s, using mixed basis level 1 as default. \n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );
	  }
	  
	  fprintf(job,"****\n");
	  
	}
      }
    } // End loop over unit cell atoms
    
    // LATTICE ATOMS:
    for (int imon=1;imon<=NMon_images;imon++) {
      if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  
	  if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	  
	  fprintf(job,"****\n");
	}
      }
    }
    
  }
  
  
  //  MIXED BASIS 
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Params::Parameters().CustomBasis() == 0 ) {
	    // Now print the basis based on the region...
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	  } else if ( Params::Parameters().CustomBasis() == 1) {
	    //printf("DEBUG: USING CUSTOM BASIS!!!\n");
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      printf("DEBUG: USING CUSTOM BASIS!!!\n");
	      if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	      } else {
		fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	      }
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	    
	  }
	  fprintf(job,"****\n");
	}
      }
    }
    
    
    for (int imon=1;imon<=NMon_images;imon++) {
      if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  // Now print the basis based on the region...
	  if ( Params::Parameters().CustomBasis() == 0 ) {
	    if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	  } else if ( Params::Parameters().CustomBasis() == 1) {
	    if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	      } else {
		fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	      }
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	    
	    
	  }
	  fprintf(job,"****\n");
	}
      }
    }
  }
  
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job); 
  
}

void CreateClusterG09Job(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, Matrix EwaldCharges, int iclust ) {

    string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  std::ostringstream ss;
  ss << iclust;
  if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
    //sprintf(filename,"%s/cluster.%d.com",path.c_str(),iclust);
    filename = path + "/" + "cluster." + ss.str() + ".com"; 
  } else {
    filename = path + "/" + "tmp" + ".com"; // HACK
  }
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("CreateClusterG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Cluster fragment composed of all monomomers within %12.2f angstrom of the first %d monomers in the input file\n", Params::Parameters().GetClusterCutoff(), Params::Parameters().GetNumAsymmetricMonomers()  );
  fprintf(job,"\n");
  
  fprintf(job,"%d %d\n", total_charge, 1 ); // NOTE: spin hard-coded to singlet for now...


  // UNIT CELL
  Monomers[iclust].PrintMonomerCartesian(job);
  for (int imon=1;imon<=NMon;imon++ ) {
    if ( Monomers[imon].GetUseInClusterCalculation() && imon != iclust ) {
      Monomers[imon].PrintMonomerCartesian(job);
    }
  }

  // IMAGE MONOMERS
  for (int imon=1;imon<=NMon_images;imon++ ) {
    if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
      MonomerImages[imon].PrintMonomerCartesian(job);
    }
  }

  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    fprintf(job,"\n"); // use a blank line to separate
    
    for (int i=Params::Parameters().GetNumAsymmetricMonomers();i<=NMon;i++) {
      if ( !Monomers[i].GetUseInClusterCalculation() && Monomers[iclust].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }
    
    for (int i=1; i<=NMon_images;i++ ) {
      if ( MonomerImages[i].GetUseInEmbedding() ) {
	if ( !MonomerImages[i].GetUseInClusterCalculation() && Monomers[iclust].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	  MonomerImages[i].PrintEmbeddingCharges(job);
	}
      }
    }
    
    // Now Print the Ewald Charges:
    // Now print the charges on the Ewald Sphere:
    for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
      if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	fprintf(job, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      } else {
	fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      }
    }


  } // END Params::Parameters().UseElectrostaticEmbedding()



  // FULL CUSTOM BASIS (not compatible with LDBS)
  if ( Params::Parameters().CustomBasis() == 2 ){
    // USE Custom Basis for all atoms in the cluster:
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    // Inside unit cell:
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    printf("WARNING: Custom Basis not specified for atom type: %s, using mixed basis level 1 as default. \n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );
	  }
	  
	  fprintf(job,"****\n");
	  
	}
      }
    } // End loop over unit cell atoms
    
    // LATTICE ATOMS:
    for (int imon=1;imon<=NMon_images;imon++) {
      if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  
	  if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	  
	  fprintf(job,"****\n");
	}
      }
    }
    
  }
  
  
  //  MIXED BASIS 
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    // loop over all the atoms and print the appropriate basis...
    for (int imon=1;imon<=NMon;imon++) {
      if ( Monomers[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  if ( Params::Parameters().CustomBasis() == 0 ) {
	    // Now print the basis based on the region...
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	  } else if ( Params::Parameters().CustomBasis() == 1) {
	    //printf("DEBUG: USING CUSTOM BASIS!!!\n");
	    if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      printf("DEBUG: USING CUSTOM BASIS!!!\n");
	      if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	      } else {
		fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	      }
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	    
	  }
	  fprintf(job,"****\n");
	}
      }
    }
    
    
    for (int imon=1;imon<=NMon_images;imon++) {
      if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
	for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	  atom_counter++;
	  fprintf(job,"%d 0\n", atom_counter);
	  
	  // Now print the basis based on the region...
	  if ( Params::Parameters().CustomBasis() == 0 ) {
	    if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	  } else if ( Params::Parameters().CustomBasis() == 1) {
	    if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	      if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	      } else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
		fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	      } else {
		fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	      }
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	    } else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    } else {
	      fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	    }
	    
	    
	  }
	  fprintf(job,"****\n");
	}
      }
    }
  }
  
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job); 


}


void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int total_charge, int iclust) { 

  // make the dalton cluster job...
   string path;
   path = Params::Parameters().GetQMPath() + "/";
   string filename;
   std::ostringstream ss;
   ss << iclust;
   if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
     filename = path + "/" + "cluster." + ss.str() + ".dal";
   } 


   // Open the input file for writing
   FILE *job;
   if (( job = fopen(filename.c_str(),"w"))==NULL) {
     printf("CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
     exit(1);
   }
   
   fprintf(job,"ATOMBASIS\n");
   fprintf(job,"\tCluster Job\n");
   fprintf(job,"\tProbably don't need this line...\n");
   fflush(stdout);

   // Now we need to find the total number of atoms in the cluster...
   int tot_nAtoms = 0;
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       //printf("DEBUG: using monomer: %d in Dalton cluster \n", Monomers[imon].GetIndex());
       tot_nAtoms += Monomers[imon].GetNumberOfAtoms();
     }
   }

   fprintf(job,"Atomtypes=%d Spherical Angstrom Nosymmetry\n",tot_nAtoms);

   // PRINT ATOMS IN UNIT CELL
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	 
	 // There is probably a better way to do this, but I'm in a hurry
	 double charge;
	 if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" ) {
	   charge = 1.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" ) {
	   charge = 6.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" ) {
	   charge = 7.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	   charge = 8.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "F" ) {
	   charge = 9.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" ) {
	   charge = 15.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" ) {
	   charge = 16.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	   charge = 17.0;
	 } else {
	   printf("ERROR: building dalton cluster job, unknown atom type..\n");
	   printf("Don't panic, just grep this line and add the atom you want in\n");
	   exit(1);
	 }
	 
	 if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
	 }
	 Monomers[imon].GetAtom(iatom).PrintQChemCartesian(job);
       }
     }
   }

    
  fprintf(job,"\n");

  
  // now print out the $dalton section:
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job);
  
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    
    path = Params::Parameters().GetQMPath() + "/";
    if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
      filename = path + "/" + "cluster." + ss.str() + ".pot";
    }

    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }


    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if ( Monomers[imon].GetUseInEmbedding() ) {
	if ( !Monomers[imon].GetUseInClusterCalculation() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }
    
      
    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");
    
    int counter = 1;
    
    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    
    
    // Print the Mutlipoles Section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    // Rank 1 (dipole)
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if ( !Monomers[imon].GetUseInClusterCalculation()  )  {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
      
    } // end rank 1
    
      // Rank 2 
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  if ( !Monomers[imon].GetUseInClusterCalculation()   ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank>=2) {
		Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	      }
	      counter++;
	    }
	  }
	}
      }
      
    } // end rank 2
    
      // Rank 3: (Octupoles)
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if ( Monomers[imon].GetUseInEmbedding()  ) {
	  if ( !Monomers[imon].GetUseInClusterCalculation()   ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++)  {
	      int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank >=3) {
		Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	      }
	      counter++;
	    }
	  }
	}
      }
      
    } //end rank 3
    
    if ( Params::Parameters().GetChargeEmbeddingRank() > 3 ) {
      printf("ERROR: multipole embedding for NMR calculations only goes up to rank 3: octupole\n");
      exit(1);
    }
    
    fclose(job);
    

  } // end electsotatic embedding
}


void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, int iclust ) { 

  // make the dalton cluster job...
   string path;
   path = Params::Parameters().GetQMPath() + "/";
   string filename;
   std::ostringstream ss;
   ss << iclust;
   if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
     filename = path + "/" + "cluster." + ss.str() + ".dal";
   } 


   // Open the input file for writing
   FILE *job;
   if (( job = fopen(filename.c_str(),"w"))==NULL) {
     printf("CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
     exit(1);
   }
   
   fprintf(job,"ATOMBASIS\n");
   fprintf(job,"\tCluster Job\n");
   fprintf(job,"\tProbably don't need this line...\n");
   fflush(stdout);

   // Now we need to find the total number of atoms in the cluster...
   int tot_nAtoms = 0;
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       //printf("DEBUG: using monomer: %d in Dalton cluster \n", Monomers[imon].GetIndex());
       tot_nAtoms += Monomers[imon].GetNumberOfAtoms();
     }
   }
   for (int imon=1;imon<=NMon_images;imon++ ) {
     if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
       //printf("DEBUG: using monomer: %d in Dalton cluster \n", MonomerImages[imon].GetIndex());
       tot_nAtoms += MonomerImages[imon].GetNumberOfAtoms();
     }
   }

   fprintf(job,"Atomtypes=%d Spherical Angstrom Nosymmetry\n",tot_nAtoms);

   // PRINT ATOMS IN UNIT CELL
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	 
	 // There is probably a better way to do this, but I'm in a hurry
	 double charge;
	 if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" ) {
	   charge = 1.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" ) {
	   charge = 6.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" ) {
	   charge = 7.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	   charge = 8.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "F" ) {
	   charge = 9.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" ) {
	   charge = 15.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" ) {
	   charge = 16.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	   charge = 17.0;
	 } else {
	   printf("ERROR: building dalton cluster job, unknown atom type..\n");
	   printf("Don't panic, just grep this line and add the atom you want in\n");
	   exit(1);
	 }
	 
	 if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
	 }
	 Monomers[imon].GetAtom(iatom).PrintQChemCartesian(job);
       }
     }
   }

   
  for (int imon=1;imon<=NMon_images;imon++) {
    //printf("\t\tDEBUG: MONOMER: %d\n", imon); fflush(stdout);
    if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
      for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	
	// There is probably a better way to do this, but I'm in a hurry
	double charge;
	if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" ) {
	  charge = 1.0;
	}	else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" ) {
	  charge = 6.0;
	}	else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" ) {
	  charge = 7.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	  charge = 8.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "F" ) {
	  charge = 9.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" ) {
	  charge = 15.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" ) {
	  charge = 16.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	  charge = 17.0;
	} else {
	  printf("ERROR: building dalton cluster job, unknown atom type..\n");
	  printf("Don't panic, just grep this line and add the atom you want in\n");
	  exit(1);
	}
	
	//printf("DEBUG: in image monomer loop mon %d, atom %d\n", imon, iatom);fflush(stdout);
	
	if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
	MonomerImages[imon].GetAtom(iatom).PrintQChemCartesian(job);
	
	//printf("DEBUG: in image monomer loop mon %d, atom %d - end loop\n", imon, iatom);fflush(stdout);
      }
    }
  }
  
  fprintf(job,"\n");

  
  // now print out the $dalton section:
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job);
  
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    
    path = Params::Parameters().GetQMPath() + "/";
    if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
      filename = path + "/" + "cluster." + ss.str() + ".pot";
    }

    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }


    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if ( Monomers[imon].GetUseInEmbedding() ) {
	if ( !Monomers[imon].GetUseInClusterCalculation() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }
    
    for (int imon=1; imon<= NMon_images ; imon++) {
      if (MonomerImages[imon].GetUseInEmbedding() ) {
	if ( !MonomerImages[imon].GetUseInClusterCalculation() ) {
	  total_atoms += MonomerImages[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++){
	    fflush(stdout);
	    int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) {
	      high++;
	    }
	  }
	}
      }
    }
      
    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");
    
    int counter = 1;
    
    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    
    // Print the mult. emb. crds. for the lattice
    for (int imon=1; imon<=NMon_images ; imon++) {
      if ( !MonomerImages[imon].GetUseInClusterCalculation()  ){ 
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  //printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	  for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    
    // Print the Mutlipoles Section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    for (int imon=1; imon<=NMon_images; imon++) {
      if ( !MonomerImages[imon].GetUseInClusterCalculation()  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	    Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    // Rank 1 (dipole)
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if ( !Monomers[imon].GetUseInClusterCalculation()  )  {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
      
      for (int imon=1; imon<= NMon_images ; imon++) {
	if ( !MonomerImages[imon].GetUseInClusterCalculation()  ) {
	  if (MonomerImages[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
    } // end rank 1
    
      // Rank 2 
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  if ( !Monomers[imon].GetUseInClusterCalculation()   ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank>=2) {
		Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	      }
	      counter++;
	    }
	  }
	}
      }
      
      for (int imon=1; imon<=NMon_images ; imon++) {
	if ( !MonomerImages[imon].GetUseInClusterCalculation()  ) {
	  if ( MonomerImages[imon].GetUseInEmbedding() ){
	    for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	      int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank>=2) {
		Vector moments = Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	      }
	      counter++;
	    }
	  }
	}
      }
    } // end rank 2
    
      // Rank 3: (Octupoles)
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if ( Monomers[imon].GetUseInEmbedding()  ) {
	  if ( !Monomers[imon].GetUseInClusterCalculation()   ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++)  {
	      int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank >=3) {
		Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	      }
	      counter++;
	    }
	  }
	}
      }
      
      for (int imon=1; imon<=NMon_images; imon++)  {
	if ( !MonomerImages[imon].GetUseInClusterCalculation()  ) {
	  if ( MonomerImages[imon].GetUseInEmbedding() ){
	    for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	      int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	      if(rank >=3)  {
		Vector moments = Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
		//Vector moments = Multipole::Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
		fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	      }
	      counter++;
	    }
	  }
	}
      }
    } //end rank 3
    
    if ( Params::Parameters().GetChargeEmbeddingRank() > 3 ) {
      printf("ERROR: multipole embedding for NMR calculations only goes up to rank 3: octupole\n");
      exit(1);
    }
    
    fclose(job);
    

  } // end electsotatic embedding
}


void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, Matrix EwaldCharges, int iclust ) { 

  // make the dalton cluster job...
   string path;
   path = Params::Parameters().GetQMPath() + "/";
   string filename;
   std::ostringstream ss;
   ss << iclust;
   if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
     //sprintf(filename,"%scluster.%d.mol", path.c_str(), iclust);
     filename = path + "/" + "cluster." + ss.str() + ".dal";
   } 

   // Open the input file for writing
   FILE *job;
   if (( job = fopen(filename.c_str(),"w"))==NULL) {
     printf("CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
     exit(1);
   }
   
   fprintf(job,"ATOMBASIS\n");
   fprintf(job,"\tCluster Job\n");
   fprintf(job,"\tProbably don't need this line...\n");
   fflush(stdout);

   // Now we need to find the total number of atoms in the cluster...
   int tot_nAtoms = 0;
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       //printf("DEBUG: using monomer: %d in Dalton cluster \n", Monomers[imon].GetIndex());
       tot_nAtoms += Monomers[imon].GetNumberOfAtoms();
     }
   }
   for (int imon=1;imon<=NMon_images;imon++ ) {
     if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
       //printf("DEBUG: using monomer: %d in Dalton cluster \n", MonomerImages[imon].GetIndex());
       tot_nAtoms += MonomerImages[imon].GetNumberOfAtoms();
     }
   }

   fprintf(job,"Atomtypes=%d Spherical Angstrom Nosymmetry\n",tot_nAtoms);

   // PRINT ATOMS IN UNIT CELL
   for (int imon=1;imon<=NMon;imon++ ) {
     if ( Monomers[imon].GetUseInClusterCalculation() ) {
       for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
	 
	 // There is probably a better way to do this, but I'm in a hurry
	 double charge;
	 if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "H" ) {
	   charge = 1.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "C" ) {
	   charge = 6.0;
	 }	else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "N" ) {
	   charge = 7.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	   charge = 8.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "F" ) {
	   charge = 9.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "P" ) {
	   charge = 15.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "S" ) {
	   charge = 16.0;
	 } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	   charge = 17.0;
	 } else {
	   printf("ERROR: building dalton cluster job, unknown atom type..\n");
	   printf("Don't panic, just grep this line and add the atom you want\n");
	   exit(1);
	 }
	 
	 if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
	 } else if ( Monomers[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	   fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
	 }
	 Monomers[imon].GetAtom(iatom).PrintQChemCartesian(job);
       }
     }
   }

   
  for (int imon=1;imon<=NMon_images;imon++) {
    //printf("\t\tDEBUG: MONOMER: %d\n", imon); fflush(stdout);
    if ( MonomerImages[imon].GetUseInClusterCalculation() ) {
      for ( int iatom=0;iatom<MonomerImages[imon].GetNumberOfAtoms();iatom++) {
	
	// There is probably a better way to do this, but I'm in a hurry
	double charge;
	if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "H" ) {
	  charge = 1.0;
	}	else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "C" ) {
	  charge = 6.0;
	}	else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "N" ) {
	  charge = 7.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	  charge = 8.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "F" ) {
	  charge = 9.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "P" ) {
	  charge = 15.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "S" ) {
	  charge = 16.0;
	} else if ( MonomerImages[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	  charge = 17.0;
	} else {
	  printf("ERROR: building dalton cluster job, unknown atom type..\n");
	  printf("Don't panic, just grep this line and add the atom you want in\n");
	  exit(1);
	}
	
	//printf("DEBUG: in image monomer loop mon %d, atom %d\n", imon, iatom);fflush(stdout);
	
	if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonomerImages[imon].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
	MonomerImages[imon].GetAtom(iatom).PrintQChemCartesian(job);
	
	//printf("DEBUG: in image monomer loop mon %d, atom %d - end loop\n", imon, iatom);fflush(stdout);
      }
    }
  }
  
  fprintf(job,"\n");
  //fclose(job);
  
  // if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
  //   filename = path + "/" + "cluster" + ".dal";
  // } else {
  //   filename = path + "/" + "tmp" + ".dal"; //print for debugging
  // }
  
  // if (( job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }
  
  // now print out the $dalton section:
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  
  fprintf(job,"\n"); // blank line to end the file...
  
  fclose(job);
  
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    
    path = Params::Parameters().GetQMPath() + "/";
    if (Params::Parameters().RunJobs() || Params::Parameters().CreateJobsOnly() ) {
      filename = path + "/" + "cluster." + ss.str() + ".pot";
    }

    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateClusterDALTONJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
      
    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if ( Monomers[imon].GetUseInEmbedding() ) {
	if ( !Monomers[imon].GetUseInClusterCalculation() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }
    
    for (int imon=1; imon<= NMon_images ; imon++) {
      if (MonomerImages[imon].GetUseInEmbedding() ) {
	if ( !MonomerImages[imon].GetUseInClusterCalculation() ) {
	  total_atoms += MonomerImages[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++){
	    fflush(stdout);
	    int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) {
	      high++;
	    }
	  }
	}
      }
    }

    total_atoms += EwaldCharges.GetRows();
    
    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");
    
    int counter = 1;
      
    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	      counter++;
	  }
	}
      }
    }
    
    // Print the mult. emb. crds. for the lattice
    for (int imon=1; imon<=NMon_images ; imon++) {
      if ( !MonomerImages[imon].GetUseInClusterCalculation()  ){ 
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  //printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	  for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    

    // Print the mult. emb. crds. for the Ewald lattice
    for ( int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++ ) {
      fprintf(job, "X  %4.8f  %4.8f  %4.8f\n", EwaldCharges.Element(ichrg,0), EwaldCharges.Element(ichrg,1), EwaldCharges.Element(ichrg,2) );
    }
    

    // Print the Mutlipoles Section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if ( !Monomers[imon].GetUseInClusterCalculation()  ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    for (int imon=1; imon<=NMon_images; imon++) {
      if ( !MonomerImages[imon].GetUseInClusterCalculation()  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	    Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    for ( int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++) {
      fprintf(job, "%d  %4.5f\n",counter,EwaldCharges.Element(ichrg,3));
      counter++;
    }

    
    if ( Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      printf("Ewald Embedding is only compatible with rank=0 (charge) embedding\n");
      exit(1);
    }

    fclose(job);
    
    
  } // end electsotatic embedding
}



void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {


  printf("\n\nBuilding Monomer and Dimer Jobs\n");
  printf("\tTwo-Body cut off: %f (angstrom)\n", Params::Parameters().GetTwoBodyCutoff());

  
  
  double LARGE_NUMBER = 999999.9;

  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }

  // MONOMER JOBS
  for (int imon=1; imon<=NTargetMonomers; imon++) {
    Monomers[imon].SetUseInEmbedding(true);
    if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
      Monomers[imon].CreateQChemJob(Monomers, NMon );
    } else if ( Params::Parameters().GetQMPackage() == "G09") {
      Monomers[imon].CreateG09Job(Monomers, NMon );
    } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
      Monomers[imon].CreateOrcaJob(Monomers, NMon );
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
      printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version molpro\n");
      exit(1);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
      Monomers[imon].CreateDaltonJob(Monomers, NMon );
    } else {
      printf("ERROR: CreateMonomerAndDimerJobs() unknown QM package: %s\n", Params::Parameters().GetQMPackage().c_str() );
      exit(1);
    }
  }

  int NDimJobs = 0;
  AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers );

  
  // DIMER JOBS
  for ( int idim=1; idim<=NDim; idim++) {
    //printf("DIMER %d \n", idim);
    if ( Dimers[idim].GetUseInTwoBodyCalculation() ) {
      NDimJobs++;
      if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	Dimers[idim].CreateQChemJob(Monomers, NMon );
      } else if ( Params::Parameters().GetQMPackage() == "G09") {
	Dimers[idim].CreateG09Job(Monomers, NMon );
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	Dimers[idim].CreateOrcaJob(Monomers, NMon );
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
	printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
	exit(1);
	//  Dimers[i].CreateMolProJob( Monomers, NMon);
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	Dimers[idim].CreateDaltonJob( Monomers, NMon );
      }
    } else {
      //printf("\t NOT USING DIMER\n");
    }
  }


  printf("\nCreated %d Monomer and %d Dimer Jobs\n", NTargetMonomers, NDimJobs);

}








void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) {


  printf("\n\nBuilding Monomer and Dimer Jobs For a Periodic System\n");
  printf("\tTwo-Body cut off: %f (angstrom)\n", Params::Parameters().GetTwoBodyCutoff());

  
  
  double LARGE_NUMBER = 999999.9;

  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }

  // MONOMER JOBS
  for (int imon=1; imon<=NTargetMonomers; imon++) {
    Monomers[imon].SetUseInEmbedding(true);
    if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
      Monomers[imon].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images);
    } else if ( Params::Parameters().GetQMPackage() == "G09") {
      Monomers[imon].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images);
    } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
      Monomers[imon].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images);
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
      printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
      exit(1);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
      Monomers[imon].CreateDaltonJob(Monomers, NMon, MonomerImages, NMon_images);
    } else {
      printf("ERROR: CreateMonomerAndDimerJobs() unknown QM package: %s\n", Params::Parameters().GetQMPackage().c_str() );
      exit(1);
    }
  }

  int NDimJobs = 0;
  AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );

  // DIMER JOBS
  for ( int idim=1; idim<=NDim; idim++) {
    //printf("DEBUG XXX TEST dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
    //printf("DIMER %d \n", idim);

    bool USE_DIMER = true;
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
	//printf("DEBUG XXX skip dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
      USE_DIMER = false;
    }

    if ( USE_DIMER ) {
      NDimJobs++;
      Dimers[idim].SetUseInTwoBodyCalculation(true);
      if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	//Dimers[idim].CreateQChemJob(Monomers, NMon );
	Dimers[idim].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "G09") {
	//Dimers[idim].CreateG09Job(Monomers, NMon );
	Dimers[idim].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	Dimers[idim].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
	printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
	exit(1);
	//  Dimers[i].CreateMolProJob( Monomers, NMon);
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	//Dimers[idim].CreateDaltonJob( Monomers, NMon );
	Dimers[idim].CreateDaltonJob(Monomers, NMon, MonomerImages, NMon_images);
      }
    } 
  }
  
  // IMAGE DIMER LOOP
  for ( int idim=1; idim<=NDim_images; idim++) {

    // bool USE_DIMER = true;
    // if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
    //   if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
    // 	USE_DIMER = false;
    //   }
    // }

    // double separation = LARGE_NUMBER;
    // for (int imon=1; imon<=NTargetMonomers; imon++ ) {
    //   separation = min( separation, Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0) );
    // }

    
    // if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
    //   USE_DIMER = false;
    // } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
    //   USE_DIMER = false;
    // }
    
    
    // if (USE_DIMER ) {
    if ( DimerImages[idim].GetUseInTwoBodyCalculation() ) {
      NDimJobs++;
      DimerImages[idim].SetUseInTwoBodyCalculation(true);
      if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	DimerImages[idim].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "G09") {
	DimerImages[idim].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	DimerImages[idim].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images);
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
	printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
	exit(1);
	//  Dimers[i].CreateMolProJob( Monomers, NMon);
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	DimerImages[idim].CreateDaltonJob( Monomers, NMon, MonomerImages, NMon_images);
      }

    }
  }


  


  printf("\nCreated %d Monomer and %d Dimer Jobs\n", NTargetMonomers, NDimJobs);

}

//void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer ImageMonomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[], Matrix EwaldCharges, int NMon_scell, Monomer SuperCellMonomers[] ) {
void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[], Matrix EwaldCharges ) {
  printf("\n\nBuilding Monomer and Dimer Jobs using Electrostatic Emebedding based on the Ewald potential\n");
  printf("\tTwo-Body cut off: %f (angstrom)\n", Params::Parameters().GetTwoBodyCutoff());

  // We largely copy 'CreateMonomerAndDimerJobs' but we add in the charges on the ewald shere:
  
  double LARGE_NUMBER = 999999.9;

  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }

  // MONOMER JOBS
  for (int imon=1; imon<=NTargetMonomers; imon++) {
    Monomers[imon].SetUseInEmbedding(true);
    if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
      //Monomers[imon].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges); //UNFINISHED
    } else if ( Params::Parameters().GetQMPackage() == "G09") {
      Monomers[imon].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
    } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
      Monomers[imon].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
      printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
      exit(1);
    } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
      Monomers[imon].CreateDaltonJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges); 
    } else {
      printf("ERROR: CreateMonomerAndDimerJobs() unknown QM package: %s\n", Params::Parameters().GetQMPackage().c_str() );
      exit(1);
    }
  }

  int NDimJobs = 0;

  // DIMER JOBS
  for ( int idim=1; idim<=NDim; idim++) {
    //printf("DEBUG XXX TEST dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
    
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
	//printf("DEBUG XXX skip dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
      USE_DIMER = false;
    }

    if ( USE_DIMER) {
      NDimJobs++;
      Dimers[idim].SetUseInTwoBodyCalculation(true);
      if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	//Dimers[idim].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      } else if ( Params::Parameters().GetQMPackage() == "G09") {
	Dimers[idim].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	Dimers[idim].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
	printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
	exit(1);
	//  Dimers[i].CreateMolProJob( Monomers, NMon);
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	Dimers[idim].CreateDaltonJob( Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      }

    }

  } 

  // IMAGE DIMER LOOP
  for ( int idim=1; idim<=NDim_images; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( DimerImages[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 ) { // && separation <= Params::Parameters().GetClusterCutoff() ) {
      
      // Actually need to see if it's within the cluster cut off for all of the cluster
      // jobs, if its not then we still need to make the dimer job:
      bool inside = true;
      for (int imon=1; imon<=NTargetMonomers; imon++ ) {
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
      NDimJobs++;
      DimerImages[idim].SetUseInTwoBodyCalculation(true);
      if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	//DimerImages[idim].CreateQChemJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      } else if ( Params::Parameters().GetQMPackage() == "G09") {
	DimerImages[idim].CreateG09Job(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	DimerImages[idim].CreateOrcaJob(Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);	     } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
	printf("ERROR: CreateMonomerAndDimerJobs() MOLPRO not implimented in current version\n");
	exit(1);
	//  Dimers[i].CreateMolProJob( Monomers, NMon);
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	DimerImages[idim].CreateDaltonJob( Monomers, NMon, MonomerImages, NMon_images, EwaldCharges);
      }

    }

  } 

  printf("\nCreated %d Monomer and %d Dimer Jobs uisng Ewald Embedding\n", NTargetMonomers, NDimJobs);
  
  

}

void RunFragmentJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {

  // CREATE LIST OF JOBS TO RUN
  int Njobs=0;
  int NMon_jobs;
  int NDim_jobs;

  // MONOMER
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    Njobs += NMon;
  } else {
    Njobs += Params::Parameters().GetNumAsymmetricMonomers();
  }
  NMon_jobs = Njobs;
  
  // DIMER
  for (int i=1;i<=NDim;i++ ) {
    if ( Dimers[i].GetUseInTwoBodyCalculation() ) {
      Njobs++;
    }
  }
  
  
  NDim_jobs = Njobs - NMon_jobs;
  
  // CLUSTER
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
      Njobs += NMon;
    } else {
      Njobs += Params::Parameters().GetNumAsymmetricMonomers();
    }
  }
  
  

  printf("Total QM Jobs:  %d \n", Njobs);
  printf("\tMonomer Jobs: %d \n", NMon_jobs);
  printf("\tDimer jobs:   %d \n", NDim_jobs);
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    printf("\tCluster Jobs: %d \n", NMon_jobs);
  }


  string *JobList = new string[Njobs];

  //JobList = GetFragmentJobList(Njobs, NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);
  JobList = GetFragmentJobList(Njobs, NMon, Monomers, NDim, Dimers );
 
  
  printf("\n\nRunning fragment jobs:\n");

#ifndef PARALLEL
  // Serial  Version 
  int ijob = 0;
  
  time_t QM_start_time, QM_stop_time;
  time_t MM_start_time, MM_stop_time;
    
  QM_start_time = time(NULL); // start the QM timer
  
  while (ijob < Njobs) {
    QM_stop_time = time(NULL);

    printf("execute job %d: %s \n",ijob,JobList[ijob].c_str());
    system(JobList[ijob].c_str());
    UpdateJobStatus(ijob,NMon_jobs,NDim_jobs);
    ijob++;
  }
  
  QM_stop_time = time(NULL);
  
  double QM_time = difftime(QM_stop_time,QM_start_time);
  

  printf("\nTime spent on QM Embedded Dipole NMR jobs: %0.f sec.\n",QM_time);
  
#else
  // Parallel Version
  int Nproc = Params::Parameters().GetNumberOfProcessors();
  if (Nproc == 1) {
    //run all the jobs in serial
    int ijob = 0;
    time_t QM_start_time, QM_stop_time;
    time_t MM_start_time, MM_stop_time;
    QM_start_time = time(NULL); // start the QM timer
    while (ijob < Njobs) {
      QM_stop_time = time(NULL);
      //printf("execute job %d: %s \n",ijob,JobList[ijob].c_str());
      system(JobList[ijob].c_str());
      UpdateJobStatus(ijob,NMon_jobs,NDim_jobs);
      ijob++;
    }
    QM_stop_time = time(NULL);
    double QM_time = difftime(QM_stop_time,QM_start_time);
    printf("\nTime spent on QM Embedded Dipole NMR jobs: %0.f sec.\n",QM_time);
  } else if (Nproc > 1) {
    // Actually running the jobs in parallel:
    printf("JDH: The jobs will be run in parallel on %d slave processors\n",Nproc-1);
    
    // load up MPI data locally
    MPI_Status status;
    int mynode,totalnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode

    // Create a few variables
    int ijob = 0; // tracks job number
    int success = 0; // outcome of jobs;
    int tag, rank; // tag = type of job, rank = processor number
    int result; 

    // Seed a job to each processor initially
    for (rank=1;rank<totalnodes;rank++) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	tag = G09_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	tag = MOLPRO_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	tag = ORCA_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	tag = DALTON_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" )  {
	tag = QCHEM_TAG;
      } else {
	printf("ERROR unknown QM Package\n");
	exit(1);
      }

      char job[BUFFSIZE];
      sprintf(job,"%s",JobList[ijob].c_str());
      MPI_Send(&job,BUFFSIZE,MPI_CHAR,rank,tag,MPI_COMM_WORLD);
      ijob++;
    }
    
    int jobs_run = 0;

    // Now continue running the jobs until all are done:
    while (ijob < Njobs) {
      MPI_Recv(&result,1,MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD,&status);
      
      success += result;

      if (Params::Parameters().GetQMPackage() == "G09" ) {
	tag = G09_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	tag = MOLPRO_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	tag = ORCA_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	tag = DALTON_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" )  {
	tag = QCHEM_TAG;
      } else {
	printf("ERROR unknown QM Package\n");
	exit(1);
      }
      


      UpdateJobStatus(jobs_run, NMon_jobs, NDim_jobs);
      jobs_run++;
      char job[BUFFSIZE];
      sprintf(job,"%s", JobList[ijob].c_str() );

      MPI_Send(&job,BUFFSIZE,MPI_CHAR, status.MPI_SOURCE,tag,MPI_COMM_WORLD);
      ijob++;
    }

    // Collect all remaining results
    // Collect all remaining results
    for (int rank = 1; rank < totalnodes; rank++) {
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      UpdateJobStatus(jobs_run, NMon_jobs, NDim_jobs);
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


  } // end of parallel code
#endif
    


}

void RunFragmentJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) {

  // CREATE LIST OF JOBS TO RUN
  int Njobs=0;
  int NMon_jobs;
  int NDim_jobs;

  // MONOMER
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    Njobs += NMon;
  } else {
    Njobs += Params::Parameters().GetNumAsymmetricMonomers();
  }
  NMon_jobs = Njobs;
  
  // DIMER
  for (int i=1;i<=NDim;i++ ) {
    if ( Dimers[i].GetUseInTwoBodyCalculation() ) {
      Njobs++;
    }
  }
  
  for (int i=1;i<=NDim_images;i++ ) {
    if ( DimerImages[i].GetUseInTwoBodyCalculation() ) {
      Njobs++;
    }
  }
  
  NDim_jobs = Njobs - NMon_jobs;
  
  // CLUSTER
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
      Njobs += NMon;
    } else {
      Njobs += Params::Parameters().GetNumAsymmetricMonomers();
    }
  }
  
  
  // end periodic
  // else {
  //   printf("ERROR non-periodic code is unfinished\n");
  //   exit(1);
  //   // UNFINISHED!!!!!!!!
  //   if ( Params::Parameters().GetClusterCutoff() > 0 ) {
  
  //   }  else {
  //     Njobs = NMon + NDim;
  //   }
  // }

  printf("Total QM Jobs:  %d \n", Njobs);
  printf("\tMonomer Jobs: %d \n", NMon_jobs);
  printf("\tDimer jobs:   %d \n", NDim_jobs);
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    printf("\tCluster Jobs: %d \n", NMon_jobs);
  }


  string *JobList = new string[Njobs];

  JobList = GetFragmentJobList(Njobs, NMon, Monomers, NDim, Dimers, NDim_images, DimerImages);

  
  printf("\n\nRunning fragment jobs:\n");

#ifndef PARALLEL
  // Serial  Version 
  int ijob = 0;
  
  time_t QM_start_time, QM_stop_time;
  time_t MM_start_time, MM_stop_time;
    
  QM_start_time = time(NULL); // start the QM timer
  
  while (ijob < Njobs) {
    QM_stop_time = time(NULL);

    printf("execute job %d: %s \n",ijob,JobList[ijob].c_str());
    system(JobList[ijob].c_str());
    UpdateJobStatus(ijob,NMon_jobs,NDim_jobs);
    ijob++;
  }
  
  QM_stop_time = time(NULL);
  
  double QM_time = difftime(QM_stop_time,QM_start_time);
  

  printf("\nTime spent on QM Embedded Dipole NMR jobs: %0.f sec.\n",QM_time);
  
#else
  // Parallel Version
  int Nproc = Params::Parameters().GetNumberOfProcessors();
  if (Nproc == 1) {
    //run all the jobs in serial
    int ijob = 0;
    time_t QM_start_time, QM_stop_time;
    time_t MM_start_time, MM_stop_time;
    QM_start_time = time(NULL); // start the QM timer
    while (ijob < Njobs) {
      QM_stop_time = time(NULL);
      //printf("execute job %d: %s \n",ijob,JobList[ijob].c_str());
      system(JobList[ijob].c_str());
      UpdateJobStatus(ijob,NMon_jobs,NDim_jobs);
      ijob++;
    }
    QM_stop_time = time(NULL);
    double QM_time = difftime(QM_stop_time,QM_start_time);
    printf("\nTime spent on QM Embedded Dipole NMR jobs: %0.f sec.\n",QM_time);
  } else if (Nproc > 1) {
    // Actually running the jobs in parallel:
    printf("JDH: The jobs will be run in parallel on %d slave processors\n",Nproc-1);
    
    // load up MPI data locally
    MPI_Status status;
    int mynode,totalnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);     // get mynode

    // Create a few variables
    int ijob = 0; // tracks job number
    int success = 0; // outcome of jobs;
    int tag, rank; // tag = type of job, rank = processor number
    int result; 

    // Seed a job to each processor initially
    for (rank=1;rank<totalnodes;rank++) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	tag = G09_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	tag = MOLPRO_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	tag = ORCA_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	tag = DALTON_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" )  {
	tag = QCHEM_TAG;
      } else {
	printf("ERROR unknown QM Package\n");
	exit(1);
      }


      char job[BUFFSIZE];
      sprintf(job,"%s",JobList[ijob].c_str());
      MPI_Send(&job,BUFFSIZE,MPI_CHAR,rank,tag,MPI_COMM_WORLD);
      ijob++;
    }
    
    int jobs_run = 0;

    // Now continue running the jobs until all are done:
    while (ijob < Njobs) {
      MPI_Recv(&result,1,MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD,&status);
      
      success += result;

      if (Params::Parameters().GetQMPackage() == "G09" ) {
	tag = G09_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	tag = MOLPRO_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	tag = ORCA_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	tag = DALTON_TAG;
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" )  {
	tag = QCHEM_TAG;
      } else {
	printf("ERROR unknown QM Package\n");
	exit(1);
      }



      UpdateJobStatus(jobs_run, NMon_jobs, NDim_jobs);
      jobs_run++;
      char job[BUFFSIZE];
      sprintf(job,"%s", JobList[ijob].c_str() );

      MPI_Send(&job,BUFFSIZE,MPI_CHAR, status.MPI_SOURCE,tag,MPI_COMM_WORLD);
      ijob++;
    }

    // Collect all remaining results
    // Collect all remaining results
    for (int rank = 1; rank < totalnodes; rank++) {
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      UpdateJobStatus(jobs_run, NMon_jobs, NDim_jobs);
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


  } // end of parallel code
#endif
    


}


string* GetFragmentJobList(int Njobs, int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {

  string *JobList = new string[Njobs];

  int ijob=0;

  // MONOMER(S)
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  for (int i=1; i<= NMon_jobs; i++) {
    if (Params::Parameters().GetQMPackage() == "G09" ) {
      JobList[ijob] = Monomers[i].RunG09Job();
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
      JobList[ijob] = Monomers[i].RunMolProJob();
    } else if ( Params::Parameters().GetQMPackage() == "QCHEM") {
      JobList[ijob] = Monomers[i].RunQChemJob();
    } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
      JobList[ijob] = Monomers[i].RunOrcaJob();
    } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
      JobList[ijob] = Monomers[i].RunDaltonJob();
    } else {
      printf("ERROR GetFragmentJobList: not prepared to automatically run %s jobs...\n", Params::Parameters().GetQMPackage().c_str() );
      exit(1);
    }
    ijob++;
  }
  
  
  // DIMER(S)
  for (int i=1;i<=NDim;i++ ) {
    if ( Dimers[i].GetUseInTwoBodyCalculation() ) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = Dimers[i].RunG09Job();
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	JobList[ijob] = Dimers[i].RunMolProJob();
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	JobList[ijob] = Dimers[i].RunQChemJob();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	JobList[ijob] = Dimers[i].RunOrcaJob();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	JobList[ijob] = Dimers[i].RunDaltonJob();
      }

      ijob++;
    }
  }
  
  
  
  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for (int i=1;i<= NMon_jobs; i++) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = RunG09ClusterJob(i);
      } else if (Params::Parameters().GetQMPackage() == "DALTON" ) {
	JobList[ijob] = RunDaltonClusterJob(i);
      } else if (Params::Parameters().GetQMPackage() == "ORCA" ) {
	//JobList[ijob] = RunOrcaClusterJob(i);
	printf("Cluster ORCA not yet ready\n");
	exit(1);
      } else {
	printf("Only prepared to run G09 or DALTON: ORCA(almost)???? cluster jobs automatically\n");
	exit(1);
      }
      ijob++;
    }
  }
  
  return JobList;

}



// This function only makes 
string* GetFragmentJobList(int Njobs, int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) {

  string *JobList = new string[Njobs];

  int ijob=0;

  // MONOMER(S)
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  for (int i=1; i<= NMon_jobs; i++) {
    if (Params::Parameters().GetQMPackage() == "G09" ) {
      JobList[ijob] = Monomers[i].RunG09Job();
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO") {
      JobList[ijob] = Monomers[i].RunMolProJob();
    } else if ( Params::Parameters().GetQMPackage() == "QCHEM") {
      JobList[ijob] = Monomers[i].RunQChemJob();
    } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
      JobList[ijob] = Monomers[i].RunOrcaJob();
    } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
      JobList[ijob] = Monomers[i].RunDaltonJob();
    } else {
      printf("ERROR GetFragmentJobList: not prepared to automatically run %s jobs...\n", Params::Parameters().GetQMPackage().c_str() );
      exit(1);
    }
    ijob++;
  }
  
  
  // DIMER(S)
  for (int i=1;i<=NDim;i++ ) {
    if ( Dimers[i].GetUseInTwoBodyCalculation() ) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = Dimers[i].RunG09Job();
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	JobList[ijob] = Dimers[i].RunMolProJob();
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	JobList[ijob] = Dimers[i].RunQChemJob();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	JobList[ijob] = Dimers[i].RunOrcaJob();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	JobList[ijob] = Dimers[i].RunDaltonJob();
      }

      ijob++;
    }
  }
  
  for (int i=1;i<=NDim_images;i++ ) {
    if ( DimerImages[i].GetUseInTwoBodyCalculation() ) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = DimerImages[i].RunG09Job();
      } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	JobList[ijob] = DimerImages[i].RunMolProJob();
      } else if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
	JobList[ijob] = DimerImages[i].RunQChemJob();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	JobList[ijob] = DimerImages[i].RunOrcaJob();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	JobList[ijob] = DimerImages[i].RunDaltonJob();
      }

      ijob++;
    }
  }
  
  
  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for (int i=1;i<= NMon_jobs; i++) {
      if (Params::Parameters().GetQMPackage() == "G09" ) {
	JobList[ijob] = RunG09ClusterJob(i);
      } else if (Params::Parameters().GetQMPackage() == "DALTON" ) {
	JobList[ijob] = RunDaltonClusterJob(i);
      } else if (Params::Parameters().GetQMPackage() == "ORCA" ) {
	//JobList[ijob] = RunOrcaClusterJob(i);
	printf("Cluster ORCA not yet ready\n");
	exit(1);
      } else {
	printf("Only prepared to run G09 or ORCA(almost)???? cluster jobs automatically\n");
	exit(1);
      }
      ijob++;
    }
  }
  

  return JobList;
  
}

string RunDaltonClusterJob(int cluster_number) {
  string path;
  path = Params::Parameters().GetQMPath() + "/";

  std::ostringstream ss;
  ss << cluster_number;

  string infile = path + "/" + "cluster." + ss.str();
  
  string local_infile = "cluster." + ss.str();
  string local_outfile = "cluster." + ss.str() + ".out";


  string cmd = "cd " + path;
  cmd += "; ";
  cmd += "dalton -mb 2000 -noarch -nobackup -o " + local_outfile + " " + local_infile; 
  cmd += "; ";
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

string RunG09ClusterJob(int cluster_number) {
  string path;
  path = Params::Parameters().GetQMPath() + "/";

  std::ostringstream ss;
  ss << cluster_number;
  string infile = path + "/" + "cluster." + ss.str() + ".com";
  //  sprintf(infile,"%scluster.%d.com",path.c_str(), cluster_number);

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";
  
  // Second command, run the job
  string local_infile;
  std::ostringstream ssb;
  ssb << cluster_number;
  local_infile = "cluster." + ssb.str() + ".com";
  
  cmd += "g09 " + local_infile;
  cmd += "; ";

  // Third command, swith back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;

}

void ReadEFGData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {

  // main routine for reading EFG data
  // grabs the Electric Filed Gradient tensors from G09 output
  // stores it as matrix objects (using the NMR stuff for now)

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  
  // MONOMER(S)
  for ( int imon=1; imon <= NMon_jobs; imon++) {
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      Monomers[imon].ReadG09EFGData();
    } else {
      printf("EFG calcs required the use of G09\n");
    }
  }

  // DIMER(S)
  for (int i=1;i<=NDim;i++) { 
    if ( Dimers[i].GetUseInTwoBodyCalculation()  ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Dimers[i].ReadG09EFGData(); 
      } else {
	printf("EFG calcs only implemented with G09\n");
	exit(1);
      }
    }
  }


  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for ( int i=1; i <= NMon_jobs; i++) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Monomers[i].ReadG09ClusterEFGData();
      } else {
	printf("EFG calcs required the use of G09\n");
      }
    }
  }

}

void ReadEFGData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]){

  // main routine for reading EFG data
  // grabs the Electric Filed Gradient tensors from G09 output
  // stores it as matrix objects (using the NMR stuff for now)

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  
  // MONOMER(S)
  for ( int imon=1; imon <= NMon_jobs; imon++) {
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      Monomers[imon].ReadG09EFGData();
    } else {
      printf("EFG calcs required the use of G09\n");
    }
  }

  // DIMER(S)
  for (int i=1;i<=NDim;i++) { 
    if ( Dimers[i].GetUseInTwoBodyCalculation()  ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Dimers[i].ReadG09EFGData(); 
      } else {
	printf("EFG calcs only implemented with G09\n");
	exit(1);
      }
    }
  }

  for (int i=1;i<=NDim_images;i++) { 
    if ( DimerImages[i].GetUseInTwoBodyCalculation()  ) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	DimerImages[i].ReadG09EFGData(); 
      } else {
	printf("EFG calcs only implemented with G09\n");
	exit(1);
      }
    }
  }
  


  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for ( int i=1; i <= NMon_jobs; i++) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Monomers[i].ReadG09ClusterEFGData();
      } else {
	printf("EFG calcs required the use of G09\n");
      }
    }
  }

}



void ComputeEFGTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {


  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }


  // Something strange with the initialization: zero out for now (fix later)
  Matrix tmp(3,3);
  for ( int imon=1;imon<=NMon_jobs;imon++) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }

  // TWO-BODY CONTRIBUTIONS (within unit cell)
  for ( int idim=1;idim<=NDim;idim++) {
    if ( Dimers[idim].GetUseInTwoBodyCalculation() ) {
      for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	Matrix tmp = Dimers[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	
	tmp -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	tmp += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );

      }

      if ( Dimers[idim].GetMonomerB().GetIndex() <= Params::Parameters().GetNumAsymmetricMonomers() ) {
	for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	  Matrix tmp = Dimers[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
	  tmp -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	  tmp += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	  Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	}
      } // end monB
      
    }
  }



  // COMPUTE FRAGMENT EFG TENSORS 
  for ( int imon=1; imon<=NMon_jobs; imon++ ) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Matrix tmp;
      
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	tmp = Monomers[imon].GetAtom(iatom).GetCluster3x3Tensor(); 
      } else {
	tmp = Monomers[imon].GetAtom(iatom).GetMonomer3x3Tensor();
      }
      tmp += Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }


  // PRINT EFG TENSORS
  printf("\n\nPrinting Final EFG Tensor Calculations \n");

  for (int imon=1;imon<= NMon_jobs; imon++ ) {
    for (int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++ ) {
      printf("\n\nMonomer: %d \t Atom: %d \t Type: %s\n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );

      Matrix tmp = Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Matrix tmpT = Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      tmpT.Transpose();
      Matrix tmpS = tmp + tmpT;
      tmpS *= 0.5;

      Vector eigs = tmpS.FindEigenvalues();


      printf("EFG Tensor @ 2-body level (MHz?):\n");

      double v33 = 0.00;
      if ( Params::Parameters().UseScaledEFGTensors() ) {

	// NEW METHOD: using the raw off diagonal and the diagonal values with the iso contribution subtracted off: Use this method if there is an issue reading in the diagonal components:
	printf("Isotropic: (not computed) \n" );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Tensor (with iso subtracted from diagonal terms):\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
	
	//double v33 = 0.00;
	for (int i = 0;i<=2; i++ ) {
	  if ( abs(eigs.Element(i) ) > v33  ) {
	    v33 = abs(eigs.Element(i) );
	  }
	}

      } else {
      // OLD METHOD: using the entire tensor in standard form:
	double iso = tmp.Trace() / 3;

	printf("Isotropic: %f\n", iso );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Tensor:\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
	
	
	//double v33 = 0.00;
	for (int i = 0;i<=2; i++ ) {
	  if ( abs(eigs.Element(i) - iso) > v33  ) {
	    v33 = abs(eigs.Element(i) - iso);
	  }
	}

      }

      if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	printf("17-O Cq: %f\n",v33 * 6.057392911 );
      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	printf("35-Cl Cq: %f\n",v33 * 19.1867568 );
      }
      

    }
  }

  printf("\n\n");

}



void ComputeEFGTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]) {


  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }


  // Something strange with the initialization: zero out for now (fix later)
  Matrix tmp(3,3);
  for ( int imon=1;imon<=NMon_jobs;imon++) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }

  // TWO-BODY CONTRIBUTIONS (within unit cell)
  for ( int idim=1;idim<=NDim;idim++) {
    if ( Dimers[idim].GetUseInTwoBodyCalculation() ) {
      for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	Matrix tmp = Dimers[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	
	tmp -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	tmp += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );

      }

      if ( Dimers[idim].GetMonomerB().GetIndex() <= Params::Parameters().GetNumAsymmetricMonomers() ) {
	for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	  Matrix tmp = Dimers[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
	  tmp -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	  tmp += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	  Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	}
      } // end monB
      
    }
  }

  // TWO-BODY CONTRIBUTIONS (with image dimers)
  for (int idim=1;idim<=NDim_images;idim++) {
    if ( DimerImages[idim].GetUseInTwoBodyCalculation()  ) {
      for (int iatom=0; iatom< DimerImages[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	Matrix tmp = DimerImages[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	tmp -= Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	tmp += Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
      }
      // Check to see if we need to read in the data for monomer B
      if ( DimerImages[idim].GetMonomerB().GetIndex() <= Params::Parameters().GetNumAsymmetricMonomers() ) {
	for (int iatom=0; iatom< DimerImages[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	  Matrix tmp = DimerImages[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
	  tmp -= Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor();
	  tmp += Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	  Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	  
	} 
      }// end monB
    }
  }
  

  // COMPUTE FRAGMENT EFG TENSORS 
  for ( int imon=1; imon<=NMon_jobs; imon++ ) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Matrix tmp;
      
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	tmp = Monomers[imon].GetAtom(iatom).GetCluster3x3Tensor(); 
      } else {
	tmp = Monomers[imon].GetAtom(iatom).GetMonomer3x3Tensor();
      }
      tmp += Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }


  // PRINT EFG TENSORS
  printf("\n\nPrinting Final EFG Tensor Calculations \n");

  for (int imon=1;imon<= NMon_jobs; imon++ ) {
    for (int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++ ) {
      printf("\n\nMonomer: %d \t Atom: %d \t Type: %s\n", imon, iatom+1, Monomers[imon].GetAtom(iatom).GetSymbol().c_str() );

      Matrix tmp = Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Matrix tmpT = Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      tmpT.Transpose();
      Matrix tmpS = tmp + tmpT;
      tmpS *= 0.5;

      Vector eigs = tmpS.FindEigenvalues();


      printf("EFG Tensor @ 2-body level (MHz?):\n");

      double v33 = 0.00;
      if ( Params::Parameters().UseScaledEFGTensors() ) {

	// NEW METHOD: using the raw off diagonal and the diagonal values with the iso contribution subtracted off: Use this method if there is an issue reading in the diagonal components:
	printf("Isotropic: (not computed) \n" );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Tensor (with iso subtracted from diagonal terms):\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
	
	//double v33 = 0.00;
	for (int i = 0;i<=2; i++ ) {
	  if ( abs(eigs.Element(i) ) > v33  ) {
	    v33 = abs(eigs.Element(i) );
	  }
	}

      } else {
      // OLD METHOD: using the entire tensor in standard form:
	double iso = tmp.Trace() / 3;

	printf("Isotropic: %f\n", iso );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Tensor:\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
	
	
	//double v33 = 0.00;
	for (int i = 0;i<=2; i++ ) {
	  if ( abs(eigs.Element(i) - iso) > v33  ) {
	    v33 = abs(eigs.Element(i) - iso);
	  }
	}

      }

      if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "O" ) {
	printf("17-O Cq: %f\n",v33 * 6.057392911 );
      } else if ( Monomers[imon].GetAtom(iatom).GetSymbol() == "Cl" ) {
	printf("35-Cl Cq: %f\n",v33 * 19.1867568 );
      }
      

    }
  }

  printf("\n\n");

}


void ReadNMRData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {

  double LARGE_NUMBER = 999999.9;
  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }
  // main routine for reading NMR data
  // grabs the Electric Filed Gradient tensors from G09 output
  // stores it as matrix objects (using the NMR stuff for now)

  AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers );
  
  
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  
  // MONOMER(S)
  for ( int imon=1; imon <= NMon_jobs; imon++) {
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      Monomers[imon].ReadG09NMRdata();
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
      Monomers[imon].ReadMolProNMRdata();
    } else if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
      Monomers[imon].ReadQChemNMRdata(); 
    } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
      Monomers[imon].ReadDaltonNMRdata(); 
    } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
      Monomers[imon].ReadOrcaNMRdata(); 
    } else {
      printf("QM Pacakge not supported\n");
      exit(1);
    }
  }

  printf("DEBUG reading dimer data\n"); 
  
  // DIMER(S)
  for (int idim=1;idim<=NDim;idim++) {
    bool USE_DIMER = true;
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
	//printf("DEBUG XXX skip dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
      USE_DIMER = false;
    }



    
    if ( USE_DIMER  ) {
      printf("DEBUG: reading in dimer(%d)\n", idim);
      if ( Params::Parameters().GetQMPackage() == "G09") {
	Dimers[idim].ReadG09NMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	Dimers[idim].ReadDaltonNMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	Dimers[idim].ReadOrcaNMRdata();
      } else {
	printf("QM Pacakge not supported\n");
      } 
    }
  }


  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for ( int i=1; i <= NMon_jobs; i++) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Monomers[i].ReadG09ClusterNMRData();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	Monomers[i].ReadDaltonClusterNMRData();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	//Monomers[i].ReadOrcaClusterNMRData();
      } else {
	printf("QM Package not supported\n");
      }
    }
  }

}





void ReadNMRData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]){

  double LARGE_NUMBER = 999999.9;
  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }

  // main routine for reading NMR data
  // grabs the Electric Filed Gradient tensors from G09 output
  // stores it as matrix objects (using the NMR stuff for now)

  AssignDimerListsForFragmentCalculation(NMon, Monomers, NDim, Dimers, NDim_images, DimerImages );
  
  
  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }

  
  // MONOMER(S)
  for ( int imon=1; imon <= NMon_jobs; imon++) {
    if ( Params::Parameters().GetQMPackage() == "G09" ) {
      Monomers[imon].ReadG09NMRdata();
    } else if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
      Monomers[imon].ReadMolProNMRdata();
    } else if ( Params::Parameters().GetQMPackage() == "QCHEM" ) {
      Monomers[imon].ReadQChemNMRdata(); 
    } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
      Monomers[imon].ReadDaltonNMRdata(); 
    } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
      Monomers[imon].ReadOrcaNMRdata(); 
    } else {
      printf("QM Pacakge not supported\n");
      exit(1);
    }
  }

  printf("DEBUG reading dimer data\n"); 
  
  // DIMER(S)
  for (int idim=1;idim<=NDim;idim++) {
    bool USE_DIMER = true;
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
	//printf("DEBUG XXX skip dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
      USE_DIMER = false;
    }

    if ( USE_DIMER ) {
    
      //if ( Dimers[idim].GetUseInTwoBodyCalculation()  ) {
      if ( Params::Parameters().GetQMPackage() == "G09") {
	Dimers[idim].ReadG09NMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	Dimers[idim].ReadDaltonNMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	Dimers[idim].ReadOrcaNMRdata();
      } else {
	printf("QM Pacakge not supported\n");
      } 
    }
  }

  for (int idim=1;idim<=NDim_images;idim++) {
    if ( DimerImages[idim].GetUseInTwoBodyCalculation()  ) {
      if ( Params::Parameters().GetQMPackage() == "G09") {
	DimerImages[idim].ReadG09NMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON") {
	DimerImages[idim].ReadDaltonNMRdata();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA") {
	DimerImages[idim].ReadOrcaNMRdata();
      } else {
	printf("QM Pacakge not supported\n");
      } 
    }
  }
 


  // CLUSTER(S)
  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    for ( int i=1; i <= NMon_jobs; i++) {
      if ( Params::Parameters().GetQMPackage() == "G09" ) {
	Monomers[i].ReadG09ClusterNMRData();
      } else if ( Params::Parameters().GetQMPackage() == "DALTON" ) {
	Monomers[i].ReadDaltonClusterNMRData();
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	//Monomers[i].ReadOrcaClusterNMRData();
	printf("OCRA cluster not finished\n");
	exit(1);
      } else {
	printf("QM Package not supported\n");
      }
    }
  }


}


void ComputeNMRShieldingTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {


  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  
  // Something strange with the initialization: zero out for now (fix later)
  Matrix tmp(3,3);
  for ( int imon=1;imon<=NMon_jobs;imon++) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }


  // COMPUTE FRAGMENT NMR TENSORS 
  for ( int imon=1; imon<=NMon_jobs; imon++ ) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Matrix tmp;
      
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	tmp = Monomers[imon].GetAtom(iatom).GetCluster3x3Tensor(); 
      } else {
	tmp = Monomers[imon].GetAtom(iatom).GetMonomer3x3Tensor();
      }
      tmp += Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }

  
  // LOOP OVER NMon_jobs
  for ( int jmon=1; jmon<= NMon_jobs; jmon++) {
    printf("------------------------------------------------------------------------\n");
    printf("\n\nPrinting Final Shielding Tensors for Monomer %d\n", jmon);



    // OPTIONALLY PRINT 1-BODY OR CLUSTER CONTRIBUTIONS
    if ( Params::Parameters().PrintLevel() >= 1) {
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	printf("\nCluster Contribution:\n");
      } else {
	printf("\n1-Body Contribution:\n");
      }

      for (int iatom=0; iatom < Monomers[jmon].GetNumberOfAtoms(); iatom++ ) {
	printf("\nMonomer: %d \t Atom: %d \t Type: %s\n", jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str() );
	
	Matrix tmp = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
	Matrix tmpT = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
	tmpT.Transpose();
	Matrix tmpS = tmp + tmpT;
	tmpS *= 0.5;
	
	Vector eigs = tmpS.FindEigenvalues();
	
	double iso = tmp.Trace() / 3;
	
	printf("Shielding Tensor @ 2-body level (ppm):\n");
	printf("Isotropic: %f\n", iso );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Shielding Tensor:\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
      }


    }
    

    // TWO-BODY CONTRIBUTIONS (within unit cell)
    for ( int idim=1;idim<=NDim;idim++) {
      if ( Dimers[idim].GetUseInTwoBodyCalculation()  ) {
	//printf("DIMER JOBS: \n");
	if ( Dimers[idim].GetMonomerA().GetIndex() == jmon ) {

	  for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	    Matrix tmp = Dimers[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	    
	    tmp -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor(); // tmp now contains only the 2-Body contribution

	    // OPTIONALLY PRINT INDIVIDUAL TWO-BODY CONTRIBUTIONS
	    if ( Params::Parameters().PrintLevel() > 1 ) {
	      Matrix tmpT = tmp;
	      tmpT.Transpose();
	      Matrix tmpS = tmp + tmpT;
	      tmpS *= 0.5;
	      
	      Vector eigs = tmpS.FindEigenvalues();
	      
	      double iso = tmp.Trace() / 3;
	      
	      printf("\n\nDimer(%d,%d) Separation: %f (angstrom)\n",Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex(), Dimers[idim].GetDimerSeparation() );
	      printf("Monomer: %d \t Atom: %d \t Type: %s Dimer(%d,%d) Contribution to the Shielding Tensor (ppm):\n",jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str(), Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex() );
	      printf("Isotropic: %f\n", iso );
	      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	      printf("Total Shielding Tensor:\n");
	      for ( int i=0; i<3;i++ ) {
		printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	      }
	    }


	    tmp += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	    Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );

	    
	  }
	}

	if ( Dimers[idim].GetMonomerB().GetIndex() == jmon ) {
	  
	  for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	    Matrix tmp = Dimers[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
	    tmp -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor(); // tmp now contains only the 2-Body contribution
	    
	    // OPTIONALLY PRINT INDIVIDUAL TWO-BODY CONTRIBUTIONS
	    if ( Params::Parameters().PrintLevel() > 1 ) {
	      Matrix tmpT = tmp;
	      tmpT.Transpose();
	      Matrix tmpS = tmp + tmpT;
	      tmpS *= 0.5;
	      
	      Vector eigs = tmpS.FindEigenvalues();
	      
	      double iso = tmp.Trace() / 3;
	      
	      printf("\n\nDimer(%d,%d) Separation: %f (angstrom)\n",Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex(), Dimers[idim].GetDimerSeparation() );
	      printf("Monomer: %d \t Atom: %d \t Type: %s Dimer(%d,%d) Contribution to the Shielding Tensor (ppm):\n",jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str(), Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex() );
	      printf("Isotropic: %f\n", iso );
	      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	      printf("Total Shielding Tensor:\n");
	      for ( int i=0; i<3;i++ ) {
		printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	      }
	    }
	    
	    
	    tmp += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	    Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	  }
	  
	} // END USE DIMER B
	
	
      }
    } // END LOOP OVER DIMERS
    
    
   
    // PRINT NMR SHIELDING TENSORS 

    printf("\n\nFinal Shielding Tensors:\n");
    for (int iatom=0; iatom < Monomers[jmon].GetNumberOfAtoms(); iatom++ ) {
      printf("\nMonomer: %d \t Atom: %d \t Type: %s\n", jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str() );
      
      Matrix tmp = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Matrix tmpT = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
      tmpT.Transpose();
      Matrix tmpS = tmp + tmpT;
      tmpS *= 0.5;
      
      Vector eigs = tmpS.FindEigenvalues();
      
      double iso = tmp.Trace() / 3;
      
      printf("Shielding Tensor @ 2-body level (ppm):\n");
      printf("Isotropic: %f\n", iso );
      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
      printf("Total Shielding Tensor:\n");
      for ( int i=0; i<3;i++ ) {
	printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
      }
    }
    

  } // END LOOP OVER NMon_jobs (jmon)


  printf("\n\n");


}


void ComputeNMRShieldingTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]) {

  double LARGE_NUMBER = 999999.9;

  int NMon_jobs;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0) {
    NMon_jobs = NMon;
  } else {
    NMon_jobs = Params::Parameters().GetNumAsymmetricMonomers();
  }
  
  // Something strange with the initialization: zero out for now (fix later)
  Matrix tmp(3,3);
  for ( int imon=1;imon<=NMon_jobs;imon++) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }


  // COMPUTE FRAGMENT NMR TENSORS 
  for ( int imon=1; imon<=NMon_jobs; imon++ ) {
    for ( int iatom=0;iatom<Monomers[imon].GetNumberOfAtoms();iatom++) {
      Matrix tmp;
      
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	tmp = Monomers[imon].GetAtom(iatom).GetCluster3x3Tensor(); 
      } else {
	tmp = Monomers[imon].GetAtom(iatom).GetMonomer3x3Tensor();
      }
      tmp += Monomers[imon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Monomers[imon].GetAtom(iatom).SetTwoBody3x3Tensor( tmp);
    }
  }

  
  // LOOP OVER NMon_jobs
  for ( int jmon=1; jmon<= NMon_jobs; jmon++) {
    printf("------------------------------------------------------------------------\n");
    printf("\n\nPrinting Final Shielding Tensors for Monomer %d\n", jmon);



    // OPTIONALLY PRINT 1-BODY OR CLUSTER CONTRIBUTIONS
    if ( Params::Parameters().PrintLevel() >= 1) {
      if ( Params::Parameters().GetClusterCutoff() > 0 ) {
	printf("\nCluster Contribution:\n");
      } else {
	printf("\1-Body Contribution:\n");
      }

      for (int iatom=0; iatom < Monomers[jmon].GetNumberOfAtoms(); iatom++ ) {
	printf("\nMonomer: %d \t Atom: %d \t Type: %s\n", jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str() );
	
	Matrix tmp = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
	Matrix tmpT = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
	tmpT.Transpose();
	Matrix tmpS = tmp + tmpT;
	tmpS *= 0.5;
	
	Vector eigs = tmpS.FindEigenvalues();
	
	double iso = tmp.Trace() / 3;
	
	printf("Shielding Tensor @ 2-body level (ppm):\n");
	printf("Isotropic: %f\n", iso );
	printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	printf("Total Shielding Tensor:\n");
	for ( int i=0; i<3;i++ ) {
	  printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	}
      }


    }
    

    // TWO-BODY CONTRIBUTIONS (within unit cell)
    for ( int idim=1;idim<=NDim;idim++) {
      
      bool USE_DIMER = true;
      if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
	if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	  USE_DIMER = false;
	  //printf("DEBUG XXX skip dimers: %d,%d\n", Dimers[idim].GetIndexA(), Dimers[idim].GetIndexB() );
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
      
      if ( USE_DIMER ) {
	
	
      //if ( Dimers[idim].GetUseInTwoBodyCalculation()  ) {
	//printf("DIMER JOBS: \n");
      if ( Dimers[idim].GetMonomerA().GetIndex() == jmon ) {

	  for (int iatom=0; iatom< Dimers[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	    Matrix tmp = Dimers[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	    
	    tmp -= Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor(); // tmp now contains only the 2-Body contribution

	    // OPTIONALLY PRINT INDIVIDUAL TWO-BODY CONTRIBUTIONS
	    if ( Params::Parameters().PrintLevel() > 1 ) {
	      Matrix tmpT = tmp;
	      tmpT.Transpose();
	      Matrix tmpS = tmp + tmpT;
	      tmpS *= 0.5;
	      
	      Vector eigs = tmpS.FindEigenvalues();
	      
	      double iso = tmp.Trace() / 3;
	      
	      printf("\n\nDimer(%d,%d) Separation: %f (angstrom)\n",Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex(), Dimers[idim].GetDimerSeparation() );
 	      printf("Monomer: %d \t Atom: %d \t Type: %s Dimer(%d,%d) Contribution to the Shielding Tensor (ppm):\n",jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str(), Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex() );
	      printf("Isotropic: %f\n", iso );
	      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	      printf("Total Shielding Tensor:\n");
	      for ( int i=0; i<3;i++ ) {
		printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	      }
	    }


	    tmp += Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	    Monomers[Dimers[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );

	    
	  }
	}

	if ( Dimers[idim].GetMonomerB().GetIndex() == jmon ) {
	  
	  for (int iatom=0; iatom< Dimers[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	    Matrix tmp = Dimers[idim].GetMonomerB().GetAtom(iatom).GetTwoBody3x3Tensor();
	    tmp -= Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor(); // tmp now contains only the 2-Body contribution
	    
	    // OPTIONALLY PRINT INDIVIDUAL TWO-BODY CONTRIBUTIONS
	    if ( Params::Parameters().PrintLevel() > 1 ) {
	      Matrix tmpT = tmp;
	      tmpT.Transpose();
	      Matrix tmpS = tmp + tmpT;
	      tmpS *= 0.5;
	      
	      Vector eigs = tmpS.FindEigenvalues();
	      
	      double iso = tmp.Trace() / 3;
	      
	      printf("\n\nDimer(%d,%d) Separation: %f (angstrom)\n",Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex(), Dimers[idim].GetDimerSeparation() );
	      printf("Monomer: %d \t Atom: %d \t Type: %s Dimer(%d,%d) Contribution to the Shielding Tensor (ppm):\n",jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str(), Dimers[idim].GetMonomerA().GetIndex(), Dimers[idim].GetMonomerB().GetIndex() );
	      printf("Isotropic: %f\n", iso );
	      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	      printf("Total Shielding Tensor:\n");
	      for ( int i=0; i<3;i++ ) {
		printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	      }
	    }
	    
	    
	    tmp += Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	    Monomers[Dimers[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	  }
	  
	} // END USE DIMER B
	
	
      }
    } // END LOOP OVER DIMERS
    
    // TWO-BODY CONTRIBUTIONS (with image dimers)
    for (int idim=1;idim<=NDim_images;idim++) {
      if ( DimerImages[idim].GetUseInTwoBodyCalculation()  ) {

	if ( DimerImages[idim].GetMonomerA().GetIndex() == jmon  ) {

	  for (int iatom=0; iatom< DimerImages[idim].GetMonomerA().GetNumberOfAtoms(); iatom++ ) {
	    Matrix tmp = DimerImages[idim].GetMonomerA().GetAtom(iatom).GetTwoBody3x3Tensor();
	    tmp -= Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetMonomer3x3Tensor(); // tmp now contains only the 2-Body contribution

	    
	    // OPTIONALLY PRINT INDIVIDUAL TWO-BODY CONTRIBUTIONS
	    if ( Params::Parameters().PrintLevel() > 1 ) {
	      Matrix tmpT = tmp;
	      tmpT.Transpose();
	      Matrix tmpS = tmp + tmpT;
	      tmpS *= 0.5;
	      
	      Vector eigs = tmpS.FindEigenvalues();
	      
	      double iso = tmp.Trace() / 3;
	      
	      printf("\n\nImage Dimer(%d,%d) Separation: %f (angstrom)\n",DimerImages[idim].GetMonomerA().GetIndex(), DimerImages[idim].GetMonomerB().GetIndex(), DimerImages[idim].GetDimerSeparation() );
	      printf("Monomer: %d \t Atom: %d \t Type: %s Image Dimer(%d,%d) Contribution to the Shielding Tensor (ppm):\n",jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str(), DimerImages[idim].GetMonomerA().GetIndex(), DimerImages[idim].GetMonomerB().GetIndex() );
	      printf("Isotropic: %f\n", iso );
	      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
	      printf("Total Shielding Tensor:\n");
	      for ( int i=0; i<3;i++ ) {
		printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
	      }
	    }

	    tmp += Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).GetTwoBody3x3Tensor();
	    Monomers[DimerImages[idim].GetMonomerA().GetIndex()].GetAtom(iatom).SetTwoBody3x3Tensor( tmp );
	  }
	  
	  // if ( DimerImages[idim].GetMonomerB().GetIndex() == jmon ) {

	  // // Check to see if we need to read in the data for monomer B

	  //   for (int iatom=0; iatom< DimerImages[idim].GetMonomerB().GetNumberOfAtoms(); iatom++ ) {
	  //     Matrix tmp = DimerImages[idim].GetMonomerB().GetAtom(iatom).GetNMRTwoBodyShieldingTensor();
	  //     tmp -= Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetNMRMonomerShieldingTensor(); 
	  //     tmp += Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).GetNMRTwoBodyShieldingTensor();
	  //     Monomers[DimerImages[idim].GetMonomerB().GetIndex()].GetAtom(iatom).SetNMRTwoBodyShieldingTensor( tmp );
	      
	  //   } 
	  // }// end monB
	  
	}
      }
    }
    
    
   
    // PRINT NMR SHIELDING TENSORS 

    printf("\n\nFinal Shielding Tensors:\n");
    for (int iatom=0; iatom < Monomers[jmon].GetNumberOfAtoms(); iatom++ ) {
      printf("\nMonomer: %d \t Atom: %d \t Type: %s\n", jmon, iatom+1, Monomers[jmon].GetAtom(iatom).GetSymbol().c_str() );
      
      Matrix tmp = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
      Matrix tmpT = Monomers[jmon].GetAtom(iatom).GetTwoBody3x3Tensor();
      tmpT.Transpose();
      Matrix tmpS = tmp + tmpT;
      tmpS *= 0.5;
      
      Vector eigs = tmpS.FindEigenvalues();
      
      double iso = tmp.Trace() / 3;
      
      printf("Shielding Tensor @ 2-body level (ppm):\n");
      printf("Isotropic: %f\n", iso );
      printf("E-vals: %f \t %f \t %f\n", eigs.Element(0), eigs.Element(1), eigs.Element(2) );
      printf("Total Shielding Tensor:\n");
      for ( int i=0; i<3;i++ ) {
	printf("%f \t %f \t %f\n",tmp.Element(i,0), tmp.Element(i,1), tmp.Element(i,2)  );
      }
    }
    

  } // END LOOP OVER NMon_jobs (jmon)


  printf("\n\n");


}
 

 void UpdateJobStatus(int ijob, int NMon_jobs, int NDim_jobs) {
  
  if ( ijob == 1 ) 
    printf("\nRunning %d QM Monomer jobs:\n", NMon_jobs);  

  // insert some spaces in output for pretty printing
  // Break each category into groups of 5 and rows of 50
  
  // QM monomer jobs
  if (ijob < NMon_jobs && ijob%5==0 && ijob )
    printf(" ");
  
  if (ijob < NMon_jobs && ijob%50==0 && ijob > 0 && ijob )
    printf("  (%d)\n ",ijob);
  
  // QM dimer jobs
  if (ijob >= NMon_jobs+1 && ijob  )
    printf(" ");

  if (ijob >= NMon_jobs+1 &&  (ijob - NMon_jobs + 1)%50==0 
      && (ijob-NMon_jobs+1) > 0 )
    printf("  (%d)\n ",ijob-NMon_jobs+1);
  
  // print a "." to represent a single job
  printf(".");

  
  if (ijob==NMon_jobs) printf("  (%d)\n",ijob+1);

  if (ijob==NMon_jobs) 
    printf("\nRunning %d Dimer jobs:\n", NDim_jobs);
  

  if ( ijob > NMon_jobs + NDim_jobs )
    printf("\nRunning Cluster job\n");


  if ( Params::Parameters().GetClusterCutoff() > 0 ) {
    if ( ijob == 2*NMon_jobs + NDim_jobs) {
      printf("\n");
    } 
  } else {
    if ( ijob == NMon_jobs + NDim_jobs) {
      printf("\n");
    } 
  }

  fflush(stdout);
}


void AssignDimerListsForFragmentCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ) {
  
  double LARGE_NUMBER = 999999.9;

  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }
  
  printf("Assign Dimer and Monomer Lists using %d monomers\n", NTargetMonomers);
  
  // DIMER JOBS
  for ( int idim=1; idim<=NDim; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    
    
    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    }
    // else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
    //   USE_DIMER = false;
    // }

    // if ( USE_DIMER) {
    //   Dimers[idim].SetUseInTwoBodyCalculation(true);
    //   //printf("DEBUG JDH: Setting dimer job (%d,%d) in two-body (before cluster considerations\n", Dimers[idim].GetMonomerA(),Dimers[idim].GetMonomerB()  );
    // }
  }

  // We want to exclude dimers that are within the cluster cutoff of all aysmmetric monomers
  if ( Params::Parameters().GetClusterCutoff() != 0 )  {
    for ( int idim=1; idim<=NDim; idim++) {
      if ( Dimers[idim].GetUseInTwoBodyCalculation() ) {
  	//printf("DEBUG JDH: using dimer job (%d,%d) in two-body (before cluster considerations\n", Dimers[idim].GetMonomerA(),Dimers[idim].GetMonomerB()  );
  	bool keepit = false;
  	for (int imon=1; imon<=NTargetMonomers; imon++ ) {
  	  if ( Dimers[idim].GetIndexA() == imon ) {
  	    double separation = Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0);
  	    if ( separation > Params::Parameters().GetClusterCutoff() ) {
  	      keepit = true;
  	    }
  	  }
  	}
  	// if ( !keepit ) {
  	//   //printf("DEBUG JDH: excluding dimer job (%d,%d) \n", Dimers[idim].GetMonomerA(),Dimers[idim].GetMonomerB()  );
  	//   Dimers[idim].SetUseInTwoBodyCalculation(false);
  	// }
      }
    }
  }


  // Identify dimers for use in TwoBody Charge Calculation:
  // if ( Params::Parameters().UseTwoBodyCharge() ) {
  //   for ( int idim=1; idim<=NDim; idim++) {
  //     Dimers[idim].SetUseInTwoBodyChargeCalculation(true); // we actually want to use all dimers in the unit cell
  //     // if ( Dimers[idim].GetDimerSeparation() > Params::Parameters().GetTwoBodyCutoff()) {
  //     // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(false);
  //     // } else {
  //     // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(true);
  //     // }
  //   }
  // }
  
  
}



//void AssignDimerListsForFragmentCalculation(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) {
void AssignDimerListsForFragmentCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) {

  double LARGE_NUMBER = 999999.9;

  int NTargetMonomers = 0;
  if ( Params::Parameters().GetNumAsymmetricMonomers() == 0 ) {
    NTargetMonomers = NMon;
  } else if ( Params::Parameters().GetNumAsymmetricMonomers() <= NMon ) {
    NTargetMonomers = Params::Parameters().GetNumAsymmetricMonomers();
  } else {
    printf("ERROR: CreateMonomerAndDimerJobs() the number of asymmetric monomers specified exceeds NMon!!!\n");
    exit(1);
  }

  // DIMER JOBS
  for ( int idim=1; idim<=NDim; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( Dimers[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    }

    // else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
    //   USE_DIMER = false;
    // }

    if ( USE_DIMER) {
      Dimers[idim].SetUseInTwoBodyCalculation(true);
    }
  }


  // Identify dimers for use in TwoBody Charge Calculation:
  // if ( Params::Parameters().UseTwoBodyCharge() ) {
  //   for ( int idim=1; idim<=NDim; idim++) {
  //     if ( Dimers[idim].GetDimerSeparation() > Params::Parameters().GetTwoBodyCutoff()) {
  // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(false);
  //     } else {
  // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(true);
  //     }
  //   }
  // }

  // IMAGE DIMER LOOP
  for ( int idim=1; idim<=NDim_images; idim++) {
    bool USE_DIMER = true;
    
    if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
      if ( DimerImages[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
	USE_DIMER = false;
      }
    }

    double separation = LARGE_NUMBER;
    for (int imon=1; imon<=NTargetMonomers; imon++ ) {
      separation = min( separation, Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0) );
    }

    if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
      USE_DIMER = false;
    } else if ( Params::Parameters().GetClusterCutoff() != 0 ) { // && separation <= Params::Parameters().GetClusterCutoff() ) {
      
      // Actually need to see if it's within the cluster cut off for all of the cluster
      // jobs, if its not then we still need to make the dimer job:
      bool inside = true;
      for (int imon=1; imon<=NTargetMonomers; imon++ ) {
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


  
  

  // We want to exclude dimers that are within the cluster cutoff of all aysmmetric monomers
  if ( Params::Parameters().GetClusterCutoff() != 0 )  {

    // UNIT CELL DIMERS
    for ( int idim=1; idim<=NDim; idim++) {
      if ( Dimers[idim].GetUseInTwoBodyCalculation() ) {
	bool keepit = false;
	for (int imon=1; imon<=NTargetMonomers; imon++ ) {
	  if ( Dimers[idim].GetIndexA() == imon ) {
	    double separation = Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0);
	    if ( separation > Params::Parameters().GetClusterCutoff() ) {
	      keepit = true;
	    }
	  }
	}
	if ( !keepit ) {
	  Dimers[idim].SetUseInTwoBodyCalculation(false);
	}
      }
    }

    // IMAGE DIMERS
    for ( int idim=1; idim<=NDim_images; idim++) {
      if ( DimerImages[idim].GetUseInTwoBodyCalculation() ) {
	bool keepit = false;
	for (int imon=1; imon<=NTargetMonomers; imon++ ) {
	  if ( DimerImages[idim].GetIndexA() == imon ) {
	    double separation = Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0);
	    if ( separation > Params::Parameters().GetClusterCutoff() ) {
	      keepit = true;
	    }
	  }
	}
	if ( !keepit ) {
	  DimerImages[idim].SetUseInTwoBodyCalculation(false);
	}
      }
    }
    
  }



  

}

  

//void AssignDimerListsForFragmentTwoBodyChargeCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ) // {
//   double LARGE_NUMBER = 999999.9;

//   int NTargetMonomers = 0;

//   NTargetMonomers = NMon;


//   // DIMER JOBS
//   for ( int idim=1; idim<=NDim; idim++) {
//     Dimers[idim].SetUseInTwoBodyCalculation(true);
//     // bool USE_DIMER = true;
    
//     // if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
//     //   if ( Dimers[idim].GetMonomerA().GetIndex() > NMon ) {
//     // 	USE_DIMER = false;
//     //   }
//     // }

//     // double separation = LARGE_NUMBER;
//     // for (int imon=1; imon<=NTargetMonomers; imon++ ) {
//     //   separation = min( separation, Monomers[imon].FindDistance( Dimers[idim].GetMonomerB()  ).Element(0) );
//     // }

//     // if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
//     //   USE_DIMER = false;
//     // }

//     // // else if ( Params::Parameters().GetClusterCutoff() != 0 && separation <= Params::Parameters().GetClusterCutoff() ) {
//     // //   USE_DIMER = false;
//     // // }

//     // if ( USE_DIMER) {
//     //   Dimers[idim].SetUseInTwoBodyCalculation(true);
//     // }
//   }


//   // Identify dimers for use in TwoBody Charge Calculation:
//   // if ( Params::Parameters().UseTwoBodyCharge() ) {
//   //   for ( int idim=1; idim<=NDim; idim++) {
//   //     if ( Dimers[idim].GetDimerSeparation() > Params::Parameters().GetTwoBodyCutoff()) {
//   // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(false);
//   //     } else {
//   // 	Dimers[idim].SetUseInTwoBodyChargeCalculation(true);
//   //     }
//   //   }
//   // }

//   // IMAGE DIMER LOOP
//   // for ( int idim=1; idim<=NDim_images; idim++) {
//   //   bool USE_DIMER = true;
    
//   //   // if ( Params::Parameters().GetNumAsymmetricMonomers() != 0 ) {
//   //   //   if ( DimerImages[idim].GetMonomerA().GetIndex() > Params::Parameters().GetNumAsymmetricMonomers() ) {
//   //   // 	USE_DIMER = false;
//   //   //   }
//   //   // }

//   //   double separation = LARGE_NUMBER;
//   //   for (int imon=1; imon<=NTargetMonomers; imon++ ) {
//   //     separation = min( separation, Monomers[imon].FindDistance( DimerImages[idim].GetMonomerB()  ).Element(0) );
//   //   }

//   //   if ( separation > Params::Parameters().GetTwoBodyCutoff()) {
//   //     USE_DIMER = false;
//   //   } 

//   //   if ( USE_DIMER) {
//   //     DimerImages[idim].SetUseInTwoBodyCalculation(true);
//   //   }

//   // } 


//   // Identify dimers for use in TwoBody Charge Calculation:
//   if ( Params::Parameters().UseTwoBodyCharge() ) {
//     for ( int idim=1; idim<=NDim_images; idim++) {
//       if ( DimerImages[idim].GetDimerSeparation() > Params::Parameters().GetTwoBodyCutoff()) {
// 	DimerImages[idim].SetUseInTwoBodyChargeCalculation(false);
//       } else {
// 	DimerImages[idim].SetUseInTwoBodyChargeCalculation(true);
//       }
//     }
//   }
  

  

// }
