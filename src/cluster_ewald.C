#include "cluster.h"
#include "constants.h"
using namespace hmbi_constants;

// Figure out the direct and reciprocal space cutoffs for the Ewald
// sum, and return the appropriate value of the kappa parameter that
// determines the partitioning between direct and reciprocal space.
double Cluster::DetermineDirectAndRecipSpaceCutoffs(int &nX, int &nY, int &nZ, int &kX, int &kY, int &kZ) {

  // Decide whether we use manual Ewald summation parameters or
  // auto-determine them.  By default, nX, nY, nZ, kX, kY, kZ, &
  // kappa_param = -1.  If that's true, then we auto-determine them.
  // If not, then we use the manual values, in which case all 6
  // parameters & kappa need to be provided by the user.  Note: if
  // only kappa_param > 0, then we can still autodetermine the
  // cutoffs.
 
  double kappa_param;
 
  // Grab the Ewald summation cutoffs: Direct space
  nX = Params::Parameters().GetDirecLoopX();
  nY = Params::Parameters().GetDirecLoopY();
  nZ = Params::Parameters().GetDirecLoopZ();
  
  // Grab the Ewald summation cutoffs: Reciprocal space
  kX = Params::Parameters().GetRecipLoopX();
  kY = Params::Parameters().GetRecipLoopY();
  kZ = Params::Parameters().GetRecipLoopZ();   
  
  bool auto_determine_ewald_parameters = true;
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
      Params::Parameters().SetEwaldKappa(kappa_param); // Save this value.  
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
 
  return kappa_param;
}

// Computes the AIFF induced multipole moments for the periodic system due to a
// large, *finite* set of periodic image multipoles inducing on the
// monomers in the unit cell.  Alternatively, one can embed in the
// Ewald potential, but that's not done here.
void Cluster::ComputePeriodicAIFFInducedMultipoles() {
  // This routine is similar to the ComputeClusterAIFFInduction.
  // The Z matrix is identical. The key difference occurs in the
  // potential on the right-hand side.  Here, this potential includes
  // contributions from periodic image monomers within the
  // polarization cutoff, and the RHS must be updated iteratively with the
  // new induced multipole moments that contribue to the crystal potential.

  // It might be worthwhile separating the function that builds V0,
  // since the rest doesn't change in a non-periodic system.  Then we
  // could just have separate routines for evaluating finite V0 or
  // Ewald V0, for instance.

 // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);

  double Eind = 0.0;
  double r_cutoff = Params::Parameters().GetMaxPolarizationRadius();
  printf("\nComputing self-consistent induction energy for the entire cluster\n   with cutoff of %.1f Angstroms.\n",r_cutoff);
  CreatePeriodicImageMonomerList(r_cutoff);

  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);


  // Create storage space for induced moments on all atoms.
  // In the process, also create some lists for array offsets used later on.  
  int Natoms = GetTotalNumberOfAtoms(); 
  int NatomsPBC = 0;
  for (int imon=1;imon<=NMon_images;imon++)
    NatomsPBC += MonomerImages[imon].GetNumberOfAtoms(); 

  printf("We have %d central cell monomers (%d atoms) and %d image monomers (%d atoms)\n",NMon,Natoms,NMon_images,NatomsPBC);

  int sizes[Natoms]; // gives dimensionality of the induced multipoles, etc for each atom
  int offset[Natoms]; // gives offset for given atom in matrices
  int Atom_count[NMon+1]; // gives atom offset for first atom of each new monomer
  int dim=0; // total size of the induced multipole array, etc.

  offset[0] = 0;
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

  // Memory estimate... doesn't include memory for Tab matrices.
  int array_sizes = 2*dim*dim + 3*dim + Natoms*25 + 25*25*Natoms*NatomsPBC; // Z, Z_copy, V0, fullV, dQ, dQall 
  double req_memory = (double) array_sizes*8.0/(1024.0*1024.0);
  
  printf("------------------------------------------------------\n");
  printf("  Periodic Crystal Induction memory requirements:\n");
  printf("                    Number of sites = %d\n",Natoms);
  printf("     Number of inducible multipoles = %d\n",dim);
  printf("              Total memory required = %.2f MB\n",req_memory);
  printf("------------------------------------------------------\n");
  fflush(stdout);

  // Initialize storage for Z, dQ, and V0
  Vector V0(dim);
  Vector dQ(dim);
  Matrix Z(dim,dim), Z_copy(dim,dim);
  Multipole *dQall = new Multipole[Natoms];
  Matrix *DampedTabs = new Matrix[Natoms*NatomsPBC]; 
  // Could shrink down memory & comp effor by shrinking DampedTabs to only blocks that are needed.



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

  // Now set Z's off-diagonal blocks (between pairs of atoms in
  // different monomers)
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int imonB=imonA+1;imonB<=NMon;imonB++) {
      // Grab the index of the appropriate dimer -- need this to
      // obtain the Tab matrices
      int idimer = DimerLookup(imonA,imonB); 

      for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
	int indexA = Atom_count[imonA] + iA;
	int NpolsA = sizes[indexA];
	int offsetA = offset[indexA];

	for (int iB=0; iB<Monomers[imonB].GetNumberOfAtoms(); iB++) {
	  int indexB = Atom_count[imonB] + iB;
	  int NpolsB = sizes[indexB];
	  int offsetB = offset[indexB];

	  Matrix Tab_block = Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB).GetBlock(NpolsA,NpolsB,1,1);
	  Z.SetBlock(Tab_block,offsetA,offsetB);
	  Tab_block.Transpose(); // Tba = Transpose(Tab)
	  Z.SetBlock(Tab_block,offsetB,offsetA);
	}
      }
    }
  }


  // Build potential V0 due to permanent multipoles on monomers lying
  // in the central unit cell and on periodic images within the
  // polarization cutoff.

  // For an atom "a" on MonA:
  // V0 = sum(atoms "b" on all other monomers MonB) Tab_tu * Qb_u

  // First get contributions due to other monomers in the central unit
  // cell
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
      // Grab the relevant elements of the tmp vector and store them in V0
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	V0[t+offsetA] = tmp[t+1];
      }

    }
  }

  // Now loop over contributions to central unit cell monomer A from
  // periodic image monomers B.
  int idimer=0;
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {
      int indexA = Atom_count[imonA] + iA;
      int nQA = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetLength();
      Vector tmp(nQA);
      for (int imonB=1;imonB<=NMon_images;imonB++) {
	  
	  int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	  Dimer tmp_dimer;
	  tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB]);
	  tmp_dimer.BuildDampedTabInteractionMatrices();
	  
	  for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {
	    Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());
	    Matrix Tab = tmp_dimer.GetDampedTabInteractionMatrix(iA,iB);
	    
	    DampedTabs[idimer] = Tab;
	    idimer++;

	    // Contract Tab with multipole moments on atom B: Tab*QB
	    tmp += Tab.MatrixTimesVector(QB.GetMoments());
	}
      }
      // Store the relevant elements of the tmp vector in V0
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	V0[t+offsetA] += tmp[t+1];
      }
    }
  }
  
  // If we were solving the full, infinitely large set of equations
  // for all molecules in the crystal interacting with the central
  // unit cell, we would just solve Z*dQ = V0 once.  However, we don't
  // do that.  So we need to capture the fact that the potential due
  // to the crystal isn't known until we know the induced moments.  So
  // we solve iteratively for dV.  Put another way, we're effectively
  // inverting a small sub-block of the Z matrix instead of the whole
  // thing, so we need to iteratively include the effects that are
  // missed by that (full self-consistency with the periodic image
  // multipoles).

  // Get ready to iterate induction contrib from crystal potential.
  // We have to keep updating the potential from periodic image monomers to
  // reflect the moments being induced on them.
  double ind_conv = Params::Parameters().GetInductionConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  if (Params::Parameters().PrintLevel() >= 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);
  int iter = 0;
  double deltaE = 100.0;
  double Eind_old; 


  while (fabs(deltaE) > Econv && iter <= Params::Parameters().GetMaxPolarizationCycles()) {
    iter++;
    Eind_old = Eind;
    Vector fullV = V0;

    if (iter > 1) { // in first iteration, induced multipoles are zero, so skip this.

      //fullV.Set(); printf("Debug: zeroing out V0 before adding induced contrib\n");

      // Loop over contributions to central unit cell monomer A from
      // periodic image monomers B.
      idimer=0;
      for (int imonA=1;imonA<=NMon;imonA++) {
	for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {
	  int indexA = Atom_count[imonA] + iA;
	  int nQA = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetLength();
	  Vector dV(nQA);
	  for (int imonB=1;imonB<=NMon_images;imonB++) {
	    
	    int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	    
	    for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {
	      int rankB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments().GetRank());
	      int indexB = Atom_count[ref_monB] + iB;
	      Multipole dQB;
	      dQB.Initialize(rankB);
	      for (int t=0;t<sizes[indexB];t++) {
		dQB(t+1) = dQall[indexB](t+1);
	      }
	      Matrix Tab = DampedTabs[idimer];
	      
	      // Contract Tab with induced multipole moments on atom B: Tab*dQB
	      dV += Tab.MatrixTimesVector(dQB.GetMoments());
	      idimer++;
	    }
	  }
	  // Store the relevant elements of the tmp vector in V0
	  int nelem = sizes[indexA];
	  int offsetA = offset[indexA];
	  for (int t=0;t<nelem;t++) {
	  fullV[t+offsetA] += dV[t+1];
	  }
	}
      }
    }
    
    // Now solve Z*dQ = -V for the induced moments dQ --> dQ =
    // Z^(-1)*V.  After iteration #1, Could use tricks to solve for
    // for subsequent RHS vectors that exploit LU decomposition of Z
    // created on the first pass, but it seems like the bottleneck
    // occurs in evaluating V0, not in solving Z*dQ=-V.  Doing so would also
    // eliminate the need for Z_copy, but again, the RAM bottleneck right now
    // is in DampedTabs right now.

    //fullV.Print("Cluster potential");
    fullV.Scale(-1.0);
    Z_copy = Z;
    Z_copy.SolveLinearEquations(fullV); // Solve Z*dQ=-V (via LU decomposition)
    dQ = fullV;
    //dQ.Print("Corresponding induced moments");

    // Compute the induction energy from these induced moments: 
    // Eind = 0.5*dQA*V0
    Eind = 0.5*dQ.DotProduct(V0);
    if (iter==1) {
      printf("iter = %d, Eind = %f kJ/mol\n",iter,Eind*HartreesToKJpermole);
    }
    else{
      deltaE = (Eind - Eind_old)*HartreesToKJpermole;
      printf("iter = %d, Eind = %f kJ/mol, DeltaE = %f kJ/mol\n",iter,Eind*HartreesToKJpermole,deltaE);
    }
    fflush(stdout);
    
    
    // Transfer the sparse list of induced multipole moments to
    // non-sparse one.
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
    
  }
  // Optionally print out the induced multipole moments
  if (Params::Parameters().PrintLevel() > -1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles: Finite Cluster ***\n");
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
  

  // Store the final induced multipole moment list in InducedMultipoles;
  int count = 0;
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
      Monomers[imonA].GetAtom(iA).SetInduceMultipoleMoments(dQall[count]);
      count++;
    }
  }
    

  // Free up memory 
  delete [] dQall;  // could rewrite to get rid of dQall... see Madelung version for example.
  delete [] DampedTabs;
  
  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > -1)
      printf("Full cluster induction energy (matrix) wall time = %.5f seconds\n",elapsed_time);


  //return Eind;

}


// Compute the classical electrostatics for the entire periodic
// crystal with periodic boundary conditions using multipolar Ewald
// summation.
double Cluster::ComputePeriodicAIFFElectrostatics() {

  printf("\nComputing permanent crystal multipole electrostatics via Ewald summation.\n");
 
  double Ees = 0.0;

  time_t start_ewald_recip_time, stop_ewald_recip_time;
  time_t start_ewald_direc_time, stop_ewald_direc_time;
  time_t start_self_intra_time, stop_self_intra_time;
  
  // Determine the Ewald summation cutoffs in direct (nX/nY/nZ) and
  // reciprocal space (kX,kY,kZ) and the corresponding kappa value
  int nX,nY,nZ,kX,kY,kZ;
  double kappa_param = DetermineDirectAndRecipSpaceCutoffs(nX,nY,nZ,kX,kY,kZ);



  // (1) Reciprocal space contribution
  start_ewald_recip_time = time(NULL);
  printf("  Begin the Ewald summation in reciprocal space\n");
  fflush(stdout);
  double Erecip = 0.0, Erecip_kn0 = 0.0;

  // grab information for building the RecipTab; and calculate the
  // convergence factor alphaa used in the self-energy correction,
  // which is the same as that in RecipTab and DirecTab
  double CellV = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  
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
      int NmomA = QA.GetLength();
      
      Vector V0(NmomA); // stores the potential for atom iA due to all other atoms

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();

	  Matrix RecipTab(NmomA,NmomB);

	  // Loop over reciprocal lattice vectors
	  for (int kx = -kX; kx<=kX; kx++){
	    for (int ky = -kY; ky<=kY; ky++){
	      for (int kz = -kZ; kz<=kZ; kz++){                  
		if (kx*kx+ky*ky+kz*kz!=0){ //if |kn|1=0
		  RecipTab += Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
										      Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
										      kx, ky, kz, cell_volume,
										      RecipCellx, RecipCelly, RecipCellz,-999.0);
		} 
	      } // end loop over kz
	    } // end loop over ky
	  } // end loop over kx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += RecipTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B
      Erecip += 0.50*QA.GetMoments().DotProduct(V0);
    } // end loop over atoms iA
  } // end loop over monomers A



 
  // (2) Reciprocal space kn=0 term.  
  // The case of L=2 has a finite limit when |kn|=0 that needs to be
  // computed.  There are two cases of L=2: dipole-dipole &
  // charge-quadrupole.
  Erecip_kn0=0.0;
  // Loop over all monomers in unit cell for A
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      int NmomA = QA.GetLength();
      Vector V0(NmomA);

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();

	  // Build RecipTab_kn0 matrix terms for L=2.	  
	  Matrix RecipTab_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
											    Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											    cell_volume, RecipCellx, RecipCelly, RecipCellz);
	  // Increment potential for atom iA due to all periodic images of atom QB
	  // use AddTo function instead of += since we have vector size mismatch
	  int max_elem = min(NmomB-1,8);
	  V0.AddTo(RecipTab_kn0.MatrixTimesVector(QB.GetMoments().GetRange(0,max_elem)),true); 
		  
	} // end loop over atoms iB
      } // end loop over monomers B

      Erecip_kn0 += 0.50*QA.GetMoments().DotProduct(V0);

    } // end loop over atoms iA
  } // end loop over monomers A

  //printf("New Erecip_kn0 = %f\n",Erecip_kn0*HartreesToKJpermole);

  stop_ewald_recip_time = time(NULL);
  double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
  printf("  Time to evaluate reciprocal space Ewald contributions = %.0f seconds\n",ewald_recip_time);
  
  
  // Alternate k=0 term based on Manby book chapter, pg 179.
  if (Params::Parameters().UseGlobalCoordinates() ) {// eqns valid only if global coords
    double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
    double V = cell_volume*AngToBohr*AngToBohr*AngToBohr;
    double constant = 2.0*pi/(3.0*CellV*perm);

    Vector Q1(3); 
    double Q00=0.0, Q20=0.0; 
    // Sum up Q00, Q1*, Q20
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
        for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
	  Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
	  Q00 += QA.GetMoments()[0];
	  Q1[0] += QA.GetMoments()[1];
	  Q1[1] += QA.GetMoments()[2];
	  Q1[2] += QA.GetMoments()[3];
	  Q20 += QA.GetMoments()[4];
	}
    }
    double Ekn0 = constant*(Q1.DotProduct(Q1) + 4*Q00*Q20);
    printf("My Erecip_kn0 = %f, Manby Ekn0 = %f kJ/mol\n",Erecip_kn0*HartreesToKJpermole,Ekn0*HartreesToKJpermole);
  }

  // (3) Direct space contribution.  
  start_ewald_direc_time = time(NULL);
  printf("  Compute the direct space Ewald contribution.\n");  
  fflush(stdout);

  double Edirect=0.0; 
  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];

  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();    

    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      Vector V0(NmomA);
      
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
	
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
          
	  Matrix DirectTab(NmomA,NmomB);

	  // loop over periodic image cells
	  for (int nx = -nX; nx<=nX; nx++){
	    for (int ny = -nY; ny<=nY; ny++){
	      for (int nz = -nZ; nz<=nZ; nz++){
		
		if (!(nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB)) { // |rAB-rn|!=0
		  Matrix tmpTab = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
											  Monomers[imonB].GetAtom(iB), 
											  RotVecB, RotAngB,  nx, ny, nz, cell_volume,
											  UnitCellx, UnitCelly, UnitCellz, -999.0);
		  DirectTab.AddBlockToMatrix(tmpTab,0,0);
		}
	      } // end loop over nz
	    } // end loop over ny
	  } // end loop over nx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += DirectTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B
      
      // We now have complete V0 for atom iA, so we can compute the electrostatic energy contrib.
      double Ees = 0.5*QA.GetMoments().DotProduct(V0);
      Edirect += Ees;

    } // end loop over atoms iA
  } // end loop over monomers A

  stop_ewald_direc_time = time(NULL);
  double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
  printf("  Time to evaluate direct space Ewald contribution = %.1f seconds\n",ewald_direc_time);


  // Compute (4) the self-interaction energy and (5) intramolecular
  // energy correction.  These are necessary because the Ewald
  // summation includes some interactions we don't want.  The
  // self-interaction term is for the charges interacting with
  // themselves, while the intramolecular term is used to subtract out
  // interactions beteween atoms in the same molecule in the central
  // unit cell.
  
  // Eself = -\sqrt(kappa^2/pi) \sum_atoms |Q_00|^2  
  // Eintra = -0.5*QA * Tab * QB
  //   (negative signs for both because they cancel out other terms
  //    in the total Ewald energy)

  start_self_intra_time = time(NULL);
  
  printf("  Begin the self-interaction and intramolecular energy correction\n");
  // constant for the calculation Uself
  double V_V = pow(CellV,1.0/3.0);
  double kappa = kappa_param/V_V;
  double alphaa = kappa*kappa; // not used...
  double Perm=4.0*pi*epsilon*1000.0/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  double ssff = kappa/sqrt(pi)/Perm;
  printf("  -----alphaa = %12.6f, kappa = %12.6f, Perm = %12.6f, ssff = %12.6f\n", alphaa,kappa,Perm,ssff);
  

  double Eself=0.0, Eintra=0.0;
  
  for (int imon = 1; imon <= NMon; imon++) { // loop over each monomer in unit cell
    int NAtoms = Monomers[imon].GetNumberOfAtoms();
    Vector RotVec(Monomers[imon].GetRotationVector());
    double RotAng = Monomers[imon].GetRotationAngle();
    for (int iA=0;iA<NAtoms;iA++)  {// loop over the atoms on this monomer 
      Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();

      Vector V0(NmomA); // stores the potential

      for (int iB=0;iB<NAtoms;iB++)  {// second loop over the atoms on this monomer 
	Multipole QB(Monomers[imon].GetAtom(iB).GetMultipoleMoments()); 
	int NmomB = QB.GetLength();
	
	// if iA==iB: self-energy correction
	if (iA==iB){
	  Eself -= ssff*QA(0)*QA(0);
	} 
	// if iA!=iB: intramolecular correction
	else {
	  // Want the standard Tab matrix for intramolecular correction.
	  Matrix Tab = Monomers[imon].GetAtom(iA).BuildInteractionMatrix(RotVec, RotAng,
									 Monomers[imon].GetAtom(iB), 
									 RotVec, RotAng,-999.0);

	  V0 += Tab.MatrixTimesVector(QB.GetMoments()); // increment the potential
	} 
      } // end loop over atoms iB

      Eintra -= 0.5*QA.GetMoments().DotProduct(V0);

    } // end loop over atoms iA
  } // end loop over monomers in the central cell
    
  stop_self_intra_time = time(NULL);
  double self_intra_time = difftime(stop_self_intra_time, start_self_intra_time);
  printf("  Time for self-interaction and intramolecular energy correction = %.0f seconds\n", self_intra_time);
  
  
  // (6) Surface term: This term equals zero under tinfoil boundary
  // conditions (our default), but it's nonzero under vacuum boundary
  // conditions.

  // Esurf = 2*pi/((2*eps + 1)V ) * | sum_A Q_00*rA + Q_1x | ^2
  // where Q_1x means dipole vector, rA is position vector for atom, and eps is the
  // dielectric (=1 for vacuum, infinity for tin-foil)
  double Esurf = 0.0; // Surface term = 0 if tin-foil boundary conditions
  string boundary_type = "Tinfoil";
  if (! Params::Parameters().TinFoilBoundaryConditions()) { // Vacuum Boundary conditions
    boundary_type = "Vacuum";
    Vector Qsurf(3);

    if (! Params::Parameters().UseGlobalCoordinates() ) {
      printf("Ewald error: Must set USE_GLOBAL_COORDINATES = TRUE to use vacuum boundary conditions\n");
      exit(1);
    }

    for (int imon = 1; imon <= NMon; imon++) {
      int NAtoms = Monomers[imon].GetNumberOfAtoms();
      for (int iA=0;iA<NAtoms;iA++)  {// loop over atoms on monomer 
	Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 

	// Permanent contribution
	Vector pos = Monomers[imon].GetAtom(iA).GetPosition();
	pos.Scale(QA(0)*AngToBohr); // 00 component times position, q*r
	Qsurf += pos; // 00 component times position, q*r
	Qsurf[0] += QA(1); // 11c component
	Qsurf[1] += QA(2); // 11s component
	Qsurf[2] += QA(3); // 10 component
      }
    }
    //Qsurf.Print("Qsurf");

    double Esurf_perm = 2*pi/(3*CellV*Perm)*Qsurf.DotProduct(Qsurf);

    Esurf = 2*pi/(3*CellV*Perm)*Qsurf.DotProduct(Qsurf);
    //printf("Surface dipole contribution under vacuum boundary conditions = %f kJ/mol\n",Esurf*HartreesToKJpermole);
  }

  // Combine all the Ewald terms, in Hartrees.
  double Etot = Edirect + Erecip + Erecip_kn0 + Eintra + Eself + Esurf;
  Lattice_E_Electrostatic_MM = Etot;

  // print all values for debug
  printf("\n");
  printf("--------------------------------------------------------\n");
  printf(" Permanent Multipole Electrostatic contributions:\n");
  printf("--------------------------------------------------------\n");
  printf("     Direct space       = %12.4f kJ/mol\n",Edirect*HartreesToKJpermole);
  printf("     Reciprocal space   = %12.4f kJ/mol\n",Erecip*HartreesToKJpermole);
  printf("     Reciprocal (k=0)   = %12.4f kJ/mol\n",Erecip_kn0*HartreesToKJpermole);
  printf("     Intramolecular     = %12.4f kJ/mol\n",Eintra*HartreesToKJpermole);
  printf("     Self-interaction   = %12.4f kJ/mol\n",Eself*HartreesToKJpermole);
  printf("     Surface term       = %12.4f kJ/mol (%s)\n",Esurf*HartreesToKJpermole,boundary_type.c_str());
  printf("--------------------------------------------------------\n");
  printf("     Total              = %12.4f kJ/mol\n",Etot*HartreesToKJpermole);
  printf("--------------------------------------------------------\n");

  
  return Etot;
}



// Compute the classical electrostatics for the entire periodic
// crystal with periodic boundary conditions using multipolar Ewald
// summation.  This routine requires that the induced multipoles are
// already computed and stored in the Cluster class member object
// InducedMultipoles.
void Cluster::ComputePeriodicAIFFElectrostaticsAndInduction() {

  printf("\nComputing permanent and induced crystal multipole electrostatics via Ewald summation.\n");

  double Ees = 0.0;

  time_t start_ewald_recip_time, stop_ewald_recip_time;
  time_t start_ewald_direc_time, stop_ewald_direc_time;
  time_t start_self_intra_time, stop_self_intra_time;
  
  // Determine the Ewald summation cutoffs in direct (nX/nY/nZ) and
  // reciprocal space (kX,kY,kZ) and the corresponding kappa value
  int nX,nY,nZ,kX,kY,kZ;
  double kappa_param = DetermineDirectAndRecipSpaceCutoffs(nX,nY,nZ,kX,kY,kZ);

  // (1) Reciprocal space contribution
  start_ewald_recip_time = time(NULL);
  printf("  Begin the Ewald summation in reciprocal space\n");
  fflush(stdout);
  double Erecip = 0.0, Erecip_kn0 = 0.0;
  double Erecip_ind = 0.0, Erecip_kn0_ind = 0.0;

  // grab information for building the RecipTab; and calculate the
  // convergence factor alphaa used in the self-energy correction,
  // which is the same as that in RecipTab and DirecTab
  double CellV = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  
  Vector RecipCellx = reciprocal_cell[0];
  Vector RecipCelly = reciprocal_cell[1];
  Vector RecipCellz = reciprocal_cell[2];

  // Loop over all monomers in unit cell for A
  int icount = 0;
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      int NmomA = QA.GetLength();
      //Multipole dQA(InducedMultipoles[icount]);
      Multipole dQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments());
      icount++; // increment this now for the next atom.
      Vector V0(NmomA); // stores the potential on atom iA due to all permanent multipoles


      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  
	  Matrix RecipTab(NmomA,NmomB);

	  // Loop over reciprocal lattice vectors
	  for (int kx = -kX; kx<=kX; kx++){
	    for (int ky = -kY; ky<=kY; ky++){
	      for (int kz = -kZ; kz<=kZ; kz++){                  
		if (kx*kx+ky*ky+kz*kz!=0){ //if |kn|1=0
		  RecipTab += Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
										      Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
										      kx, ky, kz, cell_volume,
										      RecipCellx, RecipCelly, RecipCellz,-999.0);
		} 
	      } // end loop over kz
	    } // end loop over ky
	  } // end loop over kx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += RecipTab.MatrixTimesVector(QB.GetMoments());
	} // end loop over atoms iB
      } // end loop over monomers B
      Erecip += 0.50*QA.GetMoments().DotProduct(V0);
      Erecip_ind += 0.50*dQA.GetMoments().DotProduct(V0);
    } // end loop over atoms iA
  } // end loop over monomers A

  printf("Erecip: Permanent = %f, Induced = %f\n",Erecip*HartreesToKJpermole,Erecip_ind*HartreesToKJpermole);

 
  // (2) Reciprocal space kn=0 term.  
  // The case of L=2 has a finite limit when |kn|=0 that needs to be
  // computed.  There are two cases of L=2: dipole-dipole &
  // charge-quadrupole.
  Erecip_kn0=0.0; Erecip_kn0_ind = 0.0;
  // Loop over all monomers in unit cell for A
  icount = 0;
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      int NmomA = QA.GetLength();
      //Multipole dQA(InducedMultipoles[icount]);
      Multipole dQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments());
      icount++; // increment this now for the next atom.
      Vector V0(NmomA);

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();

	  // Build RecipTab_kn0 matrix terms for L=2.	  
	  Matrix RecipTab_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
											    Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											    cell_volume, RecipCellx, RecipCelly, RecipCellz);
	  // Increment potential for atom iA due to all periodic images of atom QB
	  // use AddTo function instead of += since we have vector size mismatch
	  int max_elem = min(NmomB-1,8);
	  V0.AddTo(RecipTab_kn0.MatrixTimesVector(QB.GetMoments().GetRange(0,max_elem)),true); 
		  
	} // end loop over atoms iB
      } // end loop over monomers B
      Erecip_kn0 += 0.50*QA.GetMoments().DotProduct(V0);
      Erecip_kn0_ind += 0.50*dQA.GetMoments().DotProduct(V0);

    } // end loop over atoms iA
  } // end loop over monomers A

  //printf("New Erecip_kn0 = %f\n",Erecip_kn0*HartreesToKJpermole);

  stop_ewald_recip_time = time(NULL);
  double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
  printf("  Time to evaluate reciprocal space Ewald contributions = %.0f seconds\n",ewald_recip_time);
  
  
  // Alternate k=0 term based on Manby book chapter, pg 179.
  if (Params::Parameters().UseGlobalCoordinates() ) {// eqns valid only if global coords
    double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
    double V = cell_volume*AngToBohr*AngToBohr*AngToBohr;
    double constant = 2.0*pi/(3.0*CellV*perm);

    Vector Q1(3); 
    double Q00=0.0, Q20=0.0; 
    // Sum up Q00, Q1*, Q20
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
        for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
	  Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
	  Q00 += QA.GetMoments()[0];
	  Q1[0] += QA.GetMoments()[1];
	  Q1[1] += QA.GetMoments()[2];
	  Q1[2] += QA.GetMoments()[3];
	  Q20 += QA.GetMoments()[4];
	}
    }
    double Ekn0 = constant*(Q1.DotProduct(Q1) + 4*Q00*Q20);
    printf("My Erecip_kn0 = %f, Manby Ekn0 = %f kJ/mol\n",Erecip_kn0*HartreesToKJpermole,Ekn0*HartreesToKJpermole);
  }

  // (3) Direct space contribution.  
  start_ewald_direc_time = time(NULL);
  printf("  Compute the direct space Ewald contribution.\n");  
  fflush(stdout);

  double Edirect=0.0, Edirect_ind=0.0;
  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];

  icount = 0;
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();    

    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      //Multipole dQA(InducedMultipoles[icount]);
      Multipole dQA(Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments());
      icount++; // increment this now for the next atom.

      Vector V0(NmomA);
      
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
	
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
          
	  Matrix DirectTab(NmomA,NmomB);

	  // loop over periodic image cells
	  for (int nx = -nX; nx<=nX; nx++){
	    for (int ny = -nY; ny<=nY; ny++){
	      for (int nz = -nZ; nz<=nZ; nz++){
		
		if (!(nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB)) { // |rAB-rn|!=0
		  Matrix tmpTab = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
											  Monomers[imonB].GetAtom(iB), 
											  RotVecB, RotAngB,  nx, ny, nz, cell_volume,
											  UnitCellx, UnitCelly, UnitCellz, -999.0);
		  DirectTab.AddBlockToMatrix(tmpTab,0,0);
		}
	      } // end loop over nz
	    } // end loop over ny
	  } // end loop over nx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += DirectTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B
      
      // We now have complete V0 for atom iA, so we can compute the electrostatic energy contrib.
      Edirect += 0.5*QA.GetMoments().DotProduct(V0);
      Edirect_ind += 0.5*dQA.GetMoments().DotProduct(V0);
      

    } // end loop over atoms iA
  } // end loop over monomers A

  stop_ewald_direc_time = time(NULL);
  double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
  printf("  Time to evaluate direct space Ewald contribution = %.1f seconds\n",ewald_direc_time);


  // Compute (4) the self-interaction energy and (5) intramolecular
  // energy correction.  These are necessary because the Ewald
  // summation includes some interactions we don't want.  The
  // self-interaction term is for the charges interacting with
  // themselves, while the intramolecular term is used to subtract out
  // interactions beteween atoms in the same molecule in the central
  // unit cell.
  
  // Eself = -\sqrt(kappa^2/pi) \sum_atoms |Q_00|^2  
  // Eintra = -0.5*QA * Tab * QB
  //   (negative signs for both because they cancel out other terms
  //    in the total Ewald energy)

  start_self_intra_time = time(NULL);
  
  printf("  Begin the self-interaction and intramolecular energy correction\n");
  // constant for the calculation Uself
  double V_V = pow(CellV,1.0/3.0);
  double kappa = kappa_param/V_V;
  double alphaa = kappa*kappa; // not used...
  double Perm=4.0*pi*epsilon*1000.0/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  double ssff = kappa/sqrt(pi)/Perm;
  printf("  -----alphaa = %12.6f, kappa = %12.6f, Perm = %12.6f, ssff = %12.6f\n", alphaa,kappa,Perm,ssff);
  

  double Eself=0.0, Eintra=0.0, Eintra_ind=0.0;
  icount=0;
  for (int imon = 1; imon <= NMon; imon++) { // loop over each monomer in unit cell
    int NAtoms = Monomers[imon].GetNumberOfAtoms();
    Vector RotVec(Monomers[imon].GetRotationVector());
    double RotAng = Monomers[imon].GetRotationAngle();
    for (int iA=0;iA<NAtoms;iA++)  {// loop over the atoms on this monomer 
      Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      //Multipole dQA(InducedMultipoles[icount]);
      Multipole dQA(Monomers[imon].GetAtom(iA).GetInduceMultipoleMoments());
     icount++; // increment this now for the next atom.
      Vector V0(NmomA); // stores the potential

      for (int iB=0;iB<NAtoms;iB++)  {// second loop over the atoms on this monomer 
	Multipole QB(Monomers[imon].GetAtom(iB).GetMultipoleMoments()); 
	int NmomB = QB.GetLength();
	
	// if iA==iB: self-energy correction
	// Note: no induced charges, so no self-energy term for induced multipoles
	if (iA==iB){
	  Eself -= ssff*QA(0)*QA(0);
	} 
	// if iA!=iB: intramolecular correction
	else {
	  // Want the standard Tab matrix for intramolecular correction.
	  Matrix Tab = Monomers[imon].GetAtom(iA).BuildInteractionMatrix(RotVec, RotAng,
									 Monomers[imon].GetAtom(iB), 
									 RotVec, RotAng,-999.0);

	  V0 += Tab.MatrixTimesVector(QB.GetMoments()); // increment the potential
	} 
      } // end loop over atoms iB
      Eintra -= 0.5*QA.GetMoments().DotProduct(V0);
      Eintra_ind -= 0.5*dQA.GetMoments().DotProduct(V0);

    } // end loop over atoms iA
  } // end loop over monomers in the central cell
    
  stop_self_intra_time = time(NULL);
  double self_intra_time = difftime(stop_self_intra_time, start_self_intra_time);
  printf("  Time for self-interaction and intramolecular energy correction = %.0f seconds\n", self_intra_time);
  
  
  // (6) Surface term: This term equals zero under tinfoil boundary
  // conditions (our default), but it's nonzero under vacuum boundary
  // conditions.  This derivation only holds if we're doing
  // global coordinates.  If using local coords, we need to rotate the
  // dipole contributions to the global axis frame.

  // Esurf = 2*pi/((2*eps + 1)V ) * | sum_A Q_00*rA + Q_1x | ^2
  // where Q_1x means dipole vector, rA is position vector for atom, and eps is the
  // dielectric (=1 for vacuum, infinity for tin-foil)

  // Is induction term correct?  Surface term depends on net dipole

  double Esurf = 0.0; // Surface term = 0 if tin-foil boundary conditions
  string boundary_type = "Tinfoil";
  if (! Params::Parameters().TinFoilBoundaryConditions()) { // Vacuum Boundary conditions
    boundary_type = "Vacuum";
    Vector Qsurf(3);

    if (! Params::Parameters().UseGlobalCoordinates() ) {
      printf("Ewald error: Must set USE_GLOBAL_COORDINATES = TRUE to use vacuum boundary conditions\n");
      exit(1);
    }

    icount = 0;
    for (int imon = 1; imon <= NMon; imon++) {
      int NAtoms = Monomers[imon].GetNumberOfAtoms();
      for (int iA=0;iA<NAtoms;iA++)  {// loop over atoms on monomer 
	Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 
	//Multipole dQA(InducedMultipoles[icount]);
	Multipole dQA(Monomers[imon].GetAtom(iA).GetInduceMultipoleMoments());
	// Permanent contribution
	Vector pos = Monomers[imon].GetAtom(iA).GetPosition();
	pos.Scale(QA(0)*AngToBohr); // 00 component times position, q*r
	Qsurf += pos; // 00 component times position, q*r
	Qsurf[0] += QA(1); // 11c component
	Qsurf[1] += QA(2); // 11s component
	Qsurf[2] += QA(3); // 10 component

	// Induced contribution - No induced charge... only dipole.
	Qsurf[0] += dQA(1); // 11c component
	Qsurf[2] += dQA(2); // 11s component
	Qsurf[3] += dQA(3); // 10 component
      }
    }
    //Qsurf.Print("Qsurf");
    Esurf = 2*pi/(3*CellV*Perm)*Qsurf.DotProduct(Qsurf);
  }

  // Compute correction for what would have happened if we had damped induction.
  double Eind_correction = ComputeEwaldDampingCorrection();

  // Combine all the Ewald terms, in Hartrees.
  double Etot = Edirect + Erecip + Erecip_kn0 + Eintra + Eself + Esurf; // includes net surface term
  double Etot_ind = Edirect_ind + Erecip_ind + Erecip_kn0_ind + Eintra_ind + Eind_correction; // no Eself term

  Lattice_E_Electrostatic_MM = Etot;
  Lattice_E_Induction_MM = Etot_ind;

  // print all values for debug
  printf("\n");
  printf("--------------------------------------------------------\n");
  printf(" Permanent Multipole Electrostatic contributions:\n");
  printf("--------------------------------------------------------\n");
  printf("     Direct space       = %12.4f kJ/mol\n",Edirect*HartreesToKJpermole);
  printf("     Reciprocal space   = %12.4f kJ/mol\n",Erecip*HartreesToKJpermole);
  printf("     Reciprocal (k=0)   = %12.4f kJ/mol\n",Erecip_kn0*HartreesToKJpermole);
  printf("     Intramolecular     = %12.4f kJ/mol\n",Eintra*HartreesToKJpermole);
  printf("     Self-interaction   = %12.4f kJ/mol\n",Eself*HartreesToKJpermole);
  printf("     Net Surface term   = %12.4f kJ/mol (%s)\n",Esurf*HartreesToKJpermole,boundary_type.c_str());
  printf("--------------------------------------------------------\n");
  printf("     Permanent Total    = %12.4f kJ/mol\n",Etot*HartreesToKJpermole);
  printf("--------------------------------------------------------\n");

  printf("\n");
  printf("--------------------------------------------------------\n");
  printf(" Induced Multipole Electrostatic contributions:\n");
  printf("--------------------------------------------------------\n");
  printf("     Direct space       = %12.4f kJ/mol\n",Edirect_ind*HartreesToKJpermole);
  printf("     Reciprocal space   = %12.4f kJ/mol\n",Erecip_ind*HartreesToKJpermole);
  printf("     Reciprocal (k=0)   = %12.4f kJ/mol\n",Erecip_kn0_ind*HartreesToKJpermole);
  printf("     Intramolecular     = %12.4f kJ/mol\n",Eintra_ind*HartreesToKJpermole);
  //printf("     Surface term       =     not implemented\n",Esurf_ind*HartreesToKJpermole,boundary_type.c_str());
  printf("     Damping correction = %12.4f kJ/mol\n",Eind_correction*HartreesToKJpermole);
  printf("--------------------------------------------------------\n");
  printf("     Induced Total      = %12.4f kJ/mol\n",Etot_ind*HartreesToKJpermole);
  printf("--------------------------------------------------------\n");

 
  //return Etot;
}


// Computes the AIFF induced multipole moments for the periodic system 
// using a full Ewald/Madelung potential, rather than a finite cluster.
void Cluster::ComputePeriodicAIFFInducedMultipolesMadelung() {
  // This routine is similar to the ComputeClusterAIFFInduction.
  // The Z matrix is identical. The key difference occurs in the
  // potential on the right-hand side.  Here, this potential includes
  // contributions from periodic image monomers within the
  // polarization cutoff, and the RHS must be updated iteratively with the
  // new induced multipole moments that contribue to the crystal potential.

  // It might be worthwhile separating the function that builds V0,
  // since the rest doesn't change in a non-periodic system.  Then we
  // could just have separate routines for evaluating finite V0 or
  // Ewald V0, for instance.

 // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);

  double Eind = 0.0;
  printf("\n\nComputing induced multipoles using the Madelung potential\n");

  /*
  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);
  */

  // Create storage space for induced moments on all atoms.
  // In the process, also create some lists for array offsets used later on.  
  int Natoms = GetTotalNumberOfAtoms(); 
  printf("We have %d central cell monomers (%d atoms)\n",NMon,Natoms);

  double r_cutoff = Params::Parameters().GetMaxPolarizationRadius();
  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor
  CreatePeriodicImageMonomerList(r_cutoff);
  printf("\nIncluding induction damping correction in Madelung potential.\n");
  printf("     Tang-Toennies damping factor = %3f bohr.\n",beta_damp);
  printf("     Cluster cutoff = %.1f Angstroms.\n",r_cutoff);
  printf("     Number of image monomers = %d\n",NMon_images);
  

  int sizes[Natoms]; // gives dimensionality of the induced multipoles, etc for each atom
  int offset[Natoms]; // gives offset for given atom in matrices
  int Atom_count[NMon+1]; // gives atom offset for first atom of each new monomer
  int dim=0; // total size of the induced multipole array, etc.

  offset[0] = 0;
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

  // Memory estimate... doesn't include memory for Tab matrices.
  int array_sizes = 2*dim*dim + 3*dim + Natoms*25; // Z, Z_copy, V0, fullV, dQ, dQall 
  double req_memory = (double) array_sizes*8.0/(1024.0*1024.0);
  
  printf("------------------------------------------------------\n");
  printf("INACCURATE CHART... FIX OR DELETE THIS\n");
  printf("  Periodic Crystal Induction memory requirements:\n");
  printf("                    Number of sites = %d\n",Natoms);
  printf("     Number of inducible multipoles = %d\n",dim);
  printf("              Total memory required = %.2f MB\n",req_memory);
  printf("------------------------------------------------------\n");
  fflush(stdout);

  // Initialize storage for Z, dQ, and V0
  Vector V0(dim);
  Vector dQ(dim);
  Matrix Z(dim,dim), Z_copy(dim,dim);


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

  // Now set Z's off-diagonal blocks (between pairs of atoms in
  // different monomers)
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int imonB=imonA+1;imonB<=NMon;imonB++) {
      // Grab the index of the appropriate dimer -- need this to
      // obtain the Tab matrices
      int idimer = DimerLookup(imonA,imonB); 

      for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
	int indexA = Atom_count[imonA] + iA;
	int NpolsA = sizes[indexA];
	int offsetA = offset[indexA];

	for (int iB=0; iB<Monomers[imonB].GetNumberOfAtoms(); iB++) {
	  int indexB = Atom_count[imonB] + iB;
	  int NpolsB = sizes[indexB];
	  int offsetB = offset[indexB];

	  Matrix Tab_block = Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB).GetBlock(NpolsA,NpolsB,1,1);
	  Z.SetBlock(Tab_block,offsetA,offsetB);
	  Tab_block.Transpose(); // Tba = Transpose(Tab)
	  Z.SetBlock(Tab_block,offsetB,offsetA);
	}
      }
    }
  }


  // Build potential V0 due to permanent multipoles on monomers lying
  // in the central unit cell and on periodic images within the
  // polarization cutoff.
    time_t start_1, stop_1;
    start_1 = time(NULL);
    BuildEwaldInteractionMatrices();
    stop_1 = time(NULL);
    double time_1 = difftime(stop_1, start_1);
    printf("  Time to evaluate Ewald interaction matrices = %.0f seconds\n",time_1);
    V0 = ComputeMadelungPotentialForInduction(false);

  
  // If we were solving the full, infinitely large set of equations
  // for all molecules in the crystal interacting with the central
  // unit cell, we would just solve Z*dQ = V0 once.  However, we don't
  // do that.  So we need to capture the fact that the potential due
  // to the crystal isn't known until we know the induced moments.  So
  // we solve iteratively for dV.  Put another way, we're effectively
  // inverting a small sub-block of the Z matrix instead of the whole
  // thing, so we need to iteratively include the effects that are
  // missed by that (full self-consistency with the periodic image
  // multipoles).

  // Get ready to iterate induction contrib from crystal potential.
  // We have to keep updating the potential from periodic image monomers to
  // reflect the moments being induced on them.
  double ind_conv = Params::Parameters().GetInductionConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  if (Params::Parameters().PrintLevel() >= 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);
  int iter = 0;
  double deltaE = 100.0;
  double Eind_old; 
  while (fabs(deltaE) > Econv && iter <= Params::Parameters().GetMaxPolarizationCycles()) {
    iter++;
    Eind_old = Eind;

    // Build the Madelung potential, with or without induced contribution
    Vector fullV;
    bool induced_moments_only;
    if (iter==1) { // first iteration, induced moments are zero.
      fullV = V0;
    }
    else {  // fullV = V0 + V(induced,PBC)
      bool induced_moments_only = true; // now they are non-zero, so include them.
      fullV = ComputeMadelungPotentialForInduction(induced_moments_only);
      fullV += V0;
    }

    // Now solve Z*dQ = -V for the induced moments dQ --> dQ =
    // Z^(-1)*V.  After iteration #1, Could use tricks to solve for
    // for subsequent RHS vectors that exploit LU decomposition of Z
    // created on the first pass, but it seems like the bottleneck
    // occurs in evaluating V0, not in solving Z*dQ=-V.  Doing so would also
    // eliminate the need for Z_copy, but again, the RAM bottleneck right now
    // is in DampedTabs right now.
    fullV.Scale(-1.0);
    Z_copy = Z;
    Z_copy.SolveLinearEquations(fullV); // Solve Z*dQ=-V (via LU decomposition)
    dQ = fullV;
    
    // Compute the induction energy from these induced moments: 
    // Eind = 0.5*dQA*V0
    Eind = 0.5*dQ.DotProduct(V0);
    if (iter==1) {
      printf("Madelung: iter = %d, Eind = %f kJ/mol\n",iter,Eind*HartreesToKJpermole);
    }
    else{
      deltaE = (Eind - Eind_old)*HartreesToKJpermole;
      printf("Madelung: iter = %d, Eind = %f kJ/mol, DeltaE = %f kJ/mol\n",iter,Eind*HartreesToKJpermole,deltaE);
    }
    fflush(stdout);
    
    
    // Store the current induced multipoles in non-sparse format
    // non-sparse one.
    for (int imonA=1;imonA<=NMon;imonA++) {
      for (int iA=0; iA<Monomers[imonA].GetNumberOfAtoms(); iA++) {
	int rank = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetRank();
	int indexA = Atom_count[imonA] + iA;
	Multipole NonSparsedQ;
	NonSparsedQ.Initialize(rank);
	for (int t=0;t<sizes[indexA];t++) {
	  NonSparsedQ(t+1) = dQ[offset[indexA]+t];  // t+1 on LHS because always skip charge
	}
	Monomers[imonA].GetAtom(iA).SetInduceMultipoleMoments(NonSparsedQ);
      }
    }
    
  }


  // Optionally print out the induced multipole moments
  if (Params::Parameters().PrintLevel() > -1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles: Madelung ***\n");
    for (int imon=1;imon<=NMon;imon++) {
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	char label[50];
	sprintf(label,"Monomer %d, %s%d",imon,
		Monomers[imon].GetAtom(iA).GetSymbol().c_str(),
		Monomers[imon].GetAtom(iA).GetAtomIndex());
	string str = label;
	Monomers[imon].GetAtom(iA).GetInduceMultipoleMoments().Print(str);
      }
    }
  }
  
  
  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > -1)
      printf("Full cluster induction energy (matrix) wall time = %.5f seconds\n",elapsed_time);


  //return Eind;

}



// Compute the Madelung potential at each atom A in the unit cell for
// use in periodic multipolar induction.  The resulting Vmad is
// sparse, with only elements corresponding to dipole moments and
// sometimes quadrupoles, since those are the terms we induce.
// Specifically, Hydrogen = dipole only (3 terms), while heavy atoms =
// dipole and quadrupole (8 terms).  If induced moments only, we only
// calculate the contribution arising from induced moments in the
// self-consistent induction (i.e. the induced moments on atoms
// outside the central unit cell).
Vector Cluster::SlowComputeMadelungPotentialForInduction(bool induced_moments_only) {

  Vector Vmad;

  time_t start_ewald_recip_time, stop_ewald_recip_time;
  time_t start_ewald_direc_time, stop_ewald_direc_time;
  time_t start_self_intra_time, stop_self_intra_time;
  

  // Do some accounting to determine how much storage we need for Vmad
  // based on the atom types and how many terms we store for each.
  int Natoms = GetTotalNumberOfAtoms(); 
  int sizes[Natoms]; // gives dimensionality of the induced multipoles, etc for each atom
  int offset[Natoms]; // gives offset for given atom in matrices
  int Atom_count[NMon+1]; // gives atom offset for first atom of each new monomer
  int dim=0; // total size of the induced multipole array, etc.

  offset[0] = 0;
  int iatom = 0;
  while (iatom<Natoms) {
    for (int imon=1;imon<=NMon;imon++) {
      Atom_count[imon] = iatom;
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	int AtomicNumber = Monomers[imon].GetAtom(iA).GetAtomicNumber();
	if (AtomicNumber == 1) {
	  sizes[iatom] = 3; // dipole moments only
	}
	else {
	  sizes[iatom] = 8; // dipoles + quadrupoles
	}
	dim += sizes[iatom];
	if (iatom > 0) 
	  offset[iatom] = offset[iatom-1] + sizes[iatom-1];

	iatom++;
	//printf("atom %d, size = %d, offset = %d, current dim = %d\n",iatom-1,sizes[iatom-1],offset[iatom-1],dim);
      }
    }
  }

  // Allocate the memory for the Madelung potential
  Vmad.Initialize(dim);

  // Determine the Ewald summation cutoffs in direct (nX/nY/nZ) and
  // reciprocal space (kX,kY,kZ) and the corresponding kappa value
  int nX,nY,nZ,kX,kY,kZ;
  double kappa_param = DetermineDirectAndRecipSpaceCutoffs(nX,nY,nZ,kX,kY,kZ);

  // *** Now begin evaluating the potential *** //

  // (1) Reciprocal space contribution
  start_ewald_recip_time = time(NULL);
  //printf("  Evaluate the reciprocal space contribution to the Madelung potential\n");
  fflush(stdout);

  // grab information for building the RecipTab; and calculate the
  // convergence factor alphaa used in the self-energy correction,
  // which is the same as that in RecipTab and DirecTab
  double CellV = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  
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
      int NmomA = QA.GetLength();
       
      Vector V0(NmomA); // stores the potential for atom iA due to all other atoms

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }

	  int NmomB = QB.GetLength();

	  Matrix RecipTab(NmomA,NmomB);

	  // Loop over reciprocal lattice vectors
	  for (int kx = -kX; kx<=kX; kx++){
	    for (int ky = -kY; ky<=kY; ky++){
	      for (int kz = -kZ; kz<=kZ; kz++){                  
		if (kx*kx+ky*ky+kz*kz!=0){ //if |kn|1=0
		  RecipTab += Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
										      Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
										      kx, ky, kz, cell_volume,
										      RecipCellx, RecipCelly, RecipCellz,-999.0);
		} 
	      } // end loop over kz
	    } // end loop over ky
	  } // end loop over kx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += RecipTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B

      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A


  // (2) Reciprocal space kn=0 term.  
  // The case of L=2 has a finite limit when |kn|=0 that needs to be
  // computed.  There are two cases of L=2: dipole-dipole &
  // charge-quadrupole.
  // Loop over all monomers in unit cell for A
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      int NmomA = QA.GetLength();
      Vector V0(NmomA);

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }

	  // Build RecipTab_kn0 matrix terms for L=2.	  
	  Matrix RecipTab_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
											    Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											    cell_volume, RecipCellx, RecipCelly, RecipCellz);
	  // Increment potential for atom iA due to all periodic images of atom QB
	  // use AddTo function instead of += since we have vector size mismatch
	  int max_elem = min(NmomB-1,8);
	  V0.AddTo(RecipTab_kn0.MatrixTimesVector(QB.GetMoments().GetRange(0,max_elem)),true); 

	} // end loop over atoms iB
      } // end loop over monomers B

      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A
  

  stop_ewald_recip_time = time(NULL);
  double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
  //printf("  Time to evaluate reciprocal space Madelung contributions = %.0f seconds\n",ewald_recip_time);

  // (3) Direct space contribution.  
  start_ewald_direc_time = time(NULL);
  //printf("  Evaluate the direct space contribution to the Madelung potential\n");
  fflush(stdout);

  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];

  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();    

    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      Vector V0(NmomA);
      
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
	
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }
          
	  Matrix DirectTab(NmomA,NmomB);

	  // loop over periodic image cells
	  for (int nx = -nX; nx<=nX; nx++){
	    for (int ny = -nY; ny<=nY; ny++){
	      for (int nz = -nZ; nz<=nZ; nz++){
		
		if (!(nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB)) { // |rAB-rn|!=0
		  Matrix tmpTab = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
											  Monomers[imonB].GetAtom(iB), 
											  RotVecB, RotAngB,  nx, ny, nz, cell_volume,
											  UnitCellx, UnitCelly, UnitCellz, -999.0);
		  DirectTab.AddBlockToMatrix(tmpTab,0,0);
		}
	      } // end loop over nz
	    } // end loop over ny
	  } // end loop over nx

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += DirectTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B
      
      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A

  stop_ewald_direc_time = time(NULL);
  double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
  //printf("  Time to evaluate direct space Madelung contribution = %.1f seconds\n",ewald_direc_time);

  // Compute (4) the self-interaction energy and (5) intramolecular
  // energy correction.  These are necessary because the Ewald
  // summation includes some interactions we don't want.  The
  // self-interaction term is for the charges interacting with
  // themselves, while the intramolecular term is used to subtract out
  // interactions beteween atoms in the same molecule in the central
  // unit cell.
  
  start_self_intra_time = time(NULL);

  //printf("  Begin the self-interaction and intramolecular energy correction\n");
  // constant for the calculation Uself
  double V_V = pow(CellV,1.0/3.0);
  double kappa = kappa_param/V_V;
  double alphaa = kappa*kappa; // not used...
  double Perm=4.0*pi*epsilon*1000.0/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  double ssff = kappa/sqrt(pi)/Perm;
  //printf("  -----alphaa = %12.6f, kappa = %12.6f, Perm = %12.6f, ssff = %12.6f\n", alphaa,kappa,Perm,ssff);
  
  for (int imon = 1; imon <= NMon; imon++) { // loop over each monomer in unit cell
    int NAtoms = Monomers[imon].GetNumberOfAtoms();
    Vector RotVec(Monomers[imon].GetRotationVector());
    double RotAng = Monomers[imon].GetRotationAngle();
    for (int iA=0;iA<NAtoms;iA++)  {// loop over the atoms on this monomer 
      Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();

      Vector V0(NmomA); // stores the potential

      for (int iB=0;iB<NAtoms;iB++)  {// second loop over the atoms on this monomer 
	Multipole QB(Monomers[imon].GetAtom(iB).GetMultipoleMoments()); 
	int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imon].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }
	
	// if iA==iB: self-energy correction.  
	// Actually this charge-only self-energy term doesn't
	// contribute to the Madelung potential for induction.  Could skip
	// it, but it's trivial and leave it in case this routine ever
	// gets adapted to a full Madelung potential code.
	if (iA==iB){
	  V0[0] -= 2.0*ssff*QB(0);
	} 
	// if iA!=iB: intramolecular correction
	else {
	  // Want the standard Tab matrix for intramolecular correction.
	  Matrix Tab = Monomers[imon].GetAtom(iA).BuildInteractionMatrix(RotVec, RotAng,
									 Monomers[imon].GetAtom(iB), 
									 RotVec, RotAng,-999.0);

	  // increment the potential.  It's subtracted because in the end we are canceling
	  // the inclusion of intramolecular contributions from the earlier Ewald terms.
	  V0 -= Tab.MatrixTimesVector(QB.GetMoments()); 
	} 
      } // end loop over atoms iB

      // Now store contribution to Vmad
      int indexA = Atom_count[imon] + iA;
      int offsetA = offset[indexA];
      int nelem = sizes[indexA];
      //printf("Storing %d elements of Vmad(self/intra) for mon %d, atom %d, at %d\n",nelem,imon,iA,offsetA);
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     
    } // end loop over atoms iA
  } // end loop over monomers in the central cell

    
  stop_self_intra_time = time(NULL);
  double self_intra_time = difftime(stop_self_intra_time, start_self_intra_time);
  //printf("  Time for self-interaction and intramolecular energy correction = %.0f seconds\n", self_intra_time);


  // (6) Compute correction due to induction damping.  We use undamped
  // Tab terms in the Ewald sum, but for induction, we want to correct
  // for short-range damping.  Since this effect is purely a
  // short-range one, we can evaluate the correction using
  // conventional Tab and DampedTab matrices in a finite cluster.
  // This term only involves induced moments, not permanent ones.  The
  // correction looks like 
  // V = V(damped) - V(undamped) = (DampedTab - Tab)*QB
    for (int imonA = 1; imonA <= NMon; imonA++) {
      int NatomsA = Monomers[imonA].GetNumberOfAtoms();
      
      for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
	int NmomA = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetLength(); 
	Vector damped_V0(NmomA);

	if (! induced_moments_only) { // Case where we are treating permanent multipoles
	
	// Loop over monomers in the central unit cell
	for (int imonB = 1; imonB <= NMon; imonB++) {
	  if (imonA != imonB) {
	    int NatomsB = Monomers[imonB].GetNumberOfAtoms();

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
	    
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	      
	      Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());       
	      // Grab Tab for the current pair of atoms
	      Matrix Tab, damped_Tab;
	      if (AB_order) {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
		damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	      }
	      else {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
		Tab.Transpose();
		
		damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		damped_Tab.Transpose();
	      }
	      
	      // Compute (damped_Tab - Tab)*QB
	      damped_Tab -= Tab;
	      damped_V0 += damped_Tab.MatrixTimesVector(QB.GetMoments()); 
	      
	    
	    }  // end loop over atom iB
	  } // end if imonA != imonB
	} // end loop over monomer imonB

	}
	else {
	  // Handle the induced moments part.  Since we only want
	  // Madelung potential terms from the induced moments outside
	  // the central unit cell, we subtract out the undamped
	  // pieces which get included in steps (1)-(3), but we
	  // *don't* add in the corresponding damped ones.  In other
	  // words, we remove these contributions entirely.

	  // Loop over monomers in the central unit cell
	  for (int imonB = 1; imonB <= NMon; imonB++) {
	    if (imonA != imonB) {
	      int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	      
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
	      
	      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
		Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
		
		// Grab Tab for the current pair of atoms
		Matrix Tab;
		if (AB_order) {
		  Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
		}
		else {
		  Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
		  Tab.Transpose();
		}
		// Compute -Tab*dQB
		Tab.Scale(-1.0);
		damped_V0 += Tab.MatrixTimesVector(dQB.GetMoments()); 
		
		
	      }  // end loop over atom iB
	    } // end if imonA != imonB
	  } // end loop over monomer imonB
	}
	    
	// Now loop over periodic image monomers

	for (int imonB=1;imonB<=NMon_images;imonB++) {
	  int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	  Dimer tmp_dimer;
	  tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB]);
	  tmp_dimer.BuildTabInteractionMatrices();
	  tmp_dimer.BuildDampedTabInteractionMatrices();
	  for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {

	    Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());
	    if (induced_moments_only) {
	      QB.Set();
	      Multipole dQB(Monomers[ref_monB].GetAtom(iB).GetInduceMultipoleMoments());
	      QB += dQB;
	    }

	    Matrix Tab = tmp_dimer.GetTabInteractionMatrix(iA,iB);
	    Matrix damped_Tab = tmp_dimer.GetDampedTabInteractionMatrix(iA,iB);
	    
	    // Compute (damped_Tab - Tab)*QB
	    damped_Tab -= Tab;
	    damped_V0 += damped_Tab.MatrixTimesVector(QB.GetMoments());
	  }
	} // end loop over monomer imonB


	// Now store contribution to Vmad
	int indexA = Atom_count[imonA] + iA;
	int offsetA = offset[indexA];
	int nelem = sizes[indexA];
	for (int t=0;t<nelem;t++) {
	  Vmad[t+offsetA] += damped_V0[t+1];
      }   

      } // end loop over atom iA
    } // end loop over monomer imonA
  

  /*
  // Debug: Compute Ewald energy = 0.5*QA*Vmad
  double Eewald = 0.0;
  for (int imon = 1; imon <= NMon; imon++) { // loop over each monomer in unit cell
    int NAtoms = Monomers[imon].GetNumberOfAtoms();
    for (int iA=0;iA<NAtoms;iA++)  {// loop over the atoms on this monomer 
      Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 

      int NmomA = QA.GetLength();

      int indexA = Atom_count[imon] + iA;
      int offsetA = offset[indexA];
      int nelem = sizes[indexA];

      // Grab relevant sub-vectors of Vmad and QA
      Vector V0 = Vmad.GetRange(offsetA,offsetA+nelem-1);
      Vector qA = QA.GetMoments().GetRange(1,nelem);
      Eewald += 0.5*qA.DotProduct(V0);
    }
  }
  printf("Ewald sum from Madelung potential = %f\n",Eewald*HartreesToKJpermole);
  */

  return Vmad;
}

// Compute the difference in energy for a large cluster with and without damping.
double Cluster::ComputeEwaldDampingCorrection() {

  // Start wall clock timer
  time_t start_time, stop_time;
  start_time = time(NULL);
  
  double Eind = 0.0;
  double r_cutoff = Params::Parameters().GetMaxPolarizationRadius();
  double beta_damp = Params::Parameters().GetDampingFactor(); // damping factor
  CreatePeriodicImageMonomerList(r_cutoff);
  printf("\nComputing correction for induction damping.\n");
  printf("     Tang-Toennies damping factor = %3f bohr.\n",beta_damp);
  printf("     Cluster cutoff = %.1f Angstroms.\n",r_cutoff);

  // In the process, also create some lists for array offsets used later on.  
  int Natoms = GetTotalNumberOfAtoms(); 
  int NatomsPBC = 0;
  for (int imon=1;imon<=NMon_images;imon++) {
    NatomsPBC += MonomerImages[imon].GetNumberOfAtoms(); 
  }
 
  double Edamped = 0.0;

  // Must loop over both monomers in the central unit cell and periodic images
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {
      Multipole dQA = Monomers[imonA].GetAtom(iA).GetInduceMultipoleMoments();
      int nQA = dQA.GetLength();

      Vector damped_V0(nQA);

      // First loop over partner molecules B in the central unit cell
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
	    Matrix Tab, damped_Tab;
	    if (AB_order) {
	      Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
	      damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	    }
	    else {
	      Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
	      Tab.Transpose();

	      damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
	      damped_Tab.Transpose();
	    }

	    // Compute (damped_Tab - Tab)*QB
	    damped_Tab -= Tab;
	    damped_V0 += damped_Tab.MatrixTimesVector(QB.GetMoments());

	  } // end loop over atom iB
	} // end if (imonA != imonB)
      } // end loop over monB

      // Now loop over periodic image monomers B
      for (int imonB=1;imonB<=NMon_images;imonB++) {
	  
	  int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	  Dimer tmp_dimer;
	  tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB]);
	  tmp_dimer.BuildTabInteractionMatrices();
	  tmp_dimer.BuildDampedTabInteractionMatrices();
	  
	  for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {
	    Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());
	    Matrix Tab = tmp_dimer.GetTabInteractionMatrix(iA,iB);
	    Matrix damped_Tab = tmp_dimer.GetDampedTabInteractionMatrix(iA,iB);

	    // Compute (damped_Tab - Tab)*QB
	    damped_Tab -= Tab;
	    damped_V0 += damped_Tab.MatrixTimesVector(QB.GetMoments());
	}
      }

      Edamped += 0.5*dQA.GetMoments().DotProduct(damped_V0);
    } // end loop over atom iA
  } // end loop over monA

  printf("Damping correction = %.4f\n",Edamped*HartreesToKJpermole);

  // Stop the timer and print out the time
  stop_time = time(NULL);
  double elapsed_time = difftime(stop_time,start_time);
  if (Params::Parameters().PrintLevel() > -1)
      printf("Induction damping correction wall time = %.5f seconds\n",elapsed_time);

  return Edamped;

}




// Compute the Madelung potential at each atom A in the unit cell for
// use in periodic multipolar induction.  The resulting Vmad is
// sparse, with only elements corresponding to dipole moments and
// sometimes quadrupoles, since those are the terms we induce.
// Specifically, Hydrogen = dipole only (3 terms), while heavy atoms =
// dipole and quadrupole (8 terms).  If induced moments only, we only
// calculate the contribution arising from induced moments in the
// self-consistent induction (i.e. the induced moments on atoms
// outside the central unit cell).  
// This routine requires you to call Cluster::BuildEwaldInteractionMatrices() 
// first (once only) to create the various Tab matrices needed here.
// Note: We precompute and reuse RecipTab, DirecTab, and the damping correction,
// but we recompute the Tab terms for terms (2), (4) and (5) on the fly since
// they are really cheap.
Vector Cluster::ComputeMadelungPotentialForInduction(bool induced_moments_only) {

  Vector Vmad;

  time_t start_total, stop_total;

  time_t start_ewald_recip_time, stop_ewald_recip_time;
  time_t start_ewald_direc_time, stop_ewald_direc_time;
  time_t start_self_intra_time, stop_self_intra_time;
  time_t start_damping_time, stop_damping_time;
  
  start_total = time(NULL);

  // Do some accounting to determine how much storage we need for Vmad
  // based on the atom types and how many terms we store for each.
  int Natoms = GetTotalNumberOfAtoms(); 
  int sizes[Natoms]; // gives dimensionality of the induced multipoles, etc for each atom
  int offset[Natoms]; // gives offset for given atom in matrices
  int Atom_count[NMon+1]; // gives atom offset for first atom of each new monomer
  int dim=0; // total size of the induced multipole array, etc.

  offset[0] = 0;
  int iatom = 0;
  while (iatom<Natoms) {
    for (int imon=1;imon<=NMon;imon++) {
      Atom_count[imon] = iatom;
      for (int iA=0; iA<Monomers[imon].GetNumberOfAtoms(); iA++) {
	int AtomicNumber = Monomers[imon].GetAtom(iA).GetAtomicNumber();
	if (AtomicNumber == 1) {
	  sizes[iatom] = 3; // dipole moments only
	}
	else {
	  sizes[iatom] = 8; // dipoles + quadrupoles
	}
	dim += sizes[iatom];
	if (iatom > 0) 
	  offset[iatom] = offset[iatom-1] + sizes[iatom-1];

	iatom++;
	//printf("atom %d, size = %d, offset = %d, current dim = %d\n",iatom-1,sizes[iatom-1],offset[iatom-1],dim);
      }
    }
  }

  // Allocate the memory for the Madelung potential
  Vmad.Initialize(dim);

   // *** Now begin evaluating the potential *** //

  // (1) Reciprocal space contribution
  start_ewald_recip_time = time(NULL);
  //printf("  Evaluate the reciprocal space contribution to the Madelung potential\n");
  fflush(stdout);

  // grab information for building the RecipTab; and calculate the
  // convergence factor alphaa used in the self-energy correction,
  // which is the same as that in RecipTab and DirecTab
  double CellV = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  
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
      int NmomA = QA.GetLength();
       
      Vector V0(NmomA); // stores the potential for atom iA due to all other atoms

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    //printf("mon %d, atom %d Number of elements: QB = %d, dQB = %d\n",imonB,iB,QB.GetLength(),dQB.GetLength());
	    QB += dQB;
	  }

	  int NmomB = QB.GetLength();
	  
	  int indexAB = EwaldLookupAtomPair(imonA,iA,imonB,iB);
	  Matrix RecipTab = ReciprocalTabs[indexAB];

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += RecipTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B

      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A


  // (2) Reciprocal space kn=0 term.  
  // The case of L=2 has a finite limit when |kn|=0 that needs to be
  // computed.  There are two cases of L=2: dipole-dipole &
  // charge-quadrupole.
  // Loop over all monomers in unit cell for A
  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();
    
    // Loop over atoms on each monomer A
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); // permanant multipole for the eletrostatic energy
      int NmomA = QA.GetLength();
      Vector V0(NmomA);

      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }

	  // Build RecipTab_kn0 matrix terms for L=2.	  
	  Matrix RecipTab_kn0 = Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix_kn0(RotVecA, RotAngA,
											    Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
											    cell_volume, RecipCellx, RecipCelly, RecipCellz);
	  // Increment potential for atom iA due to all periodic images of atom QB
	  // use AddTo function instead of += since we have vector size mismatch
	  int max_elem = min(NmomB-1,8);
	  V0.AddTo(RecipTab_kn0.MatrixTimesVector(QB.GetMoments().GetRange(0,max_elem)),true); 
		  
	} // end loop over atoms iB
      } // end loop over monomers B

      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      //printf("Storing %d elements of Vmad(recip_kn0) for mon %d, atom %d, at %d\n",nelem,imonA,iA,offsetA);
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A
  

  stop_ewald_recip_time = time(NULL);
  double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
  //printf("  Time to evaluate reciprocal space Madelung contributions = %.0f seconds\n",ewald_recip_time);

  // (3) Direct space contribution.  
  start_ewald_direc_time = time(NULL);
  //printf("  Evaluate the direct space contribution to the Madelung potential\n");
  fflush(stdout);

  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];

  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();    

    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      Vector V0(NmomA);
      
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
	
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }
          
	  int indexAB = EwaldLookupAtomPair(imonA,iA,imonB,iB);
	  Matrix DirectTab = DirectTabs[indexAB];

	  // Increment potential for atom iA due to all periodic images of atom QB
	  V0 += DirectTab.MatrixTimesVector(QB.GetMoments());

	} // end loop over atoms iB
      } // end loop over monomers B
      
      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int nelem = sizes[indexA];
      int offsetA = offset[indexA];
      //printf("Storing %d elements of Vmad(direct) for mon %d, atom %d, at %d\n",nelem,imonA,iA,offsetA);
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     

    } // end loop over atoms iA
  } // end loop over monomers A

  stop_ewald_direc_time = time(NULL);
  double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
  //printf("  Time to evaluate direct space Madelung contribution = %.1f seconds\n",ewald_direc_time);

  // Compute (4) the self-interaction energy and (5) intramolecular
  // energy correction.  These are necessary because the Ewald
  // summation includes some interactions we don't want.  The
  // self-interaction term is for the charges interacting with
  // themselves, while the intramolecular term is used to subtract out
  // interactions beteween atoms in the same molecule in the central
  // unit cell.
  
  start_self_intra_time = time(NULL);

  //printf("  Begin the self-interaction and intramolecular energy correction\n");
  // constant for the calculation Uself
  double V_V = pow(CellV,1.0/3.0);
  double kappa_param = Params::Parameters().GetEwaldKappa();
  double kappa = kappa_param/V_V;
  double alphaa = kappa*kappa; // not used...
  double Perm=4.0*pi*epsilon*1000.0/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  double ssff = kappa/sqrt(pi)/Perm;
  //printf("  -----alphaa = %12.6f, kappa = %12.6f, Perm = %12.6f, ssff = %12.6f\n", alphaa,kappa,Perm,ssff);
  
  for (int imon = 1; imon <= NMon; imon++) { // loop over each monomer in unit cell
    int NAtoms = Monomers[imon].GetNumberOfAtoms();
    Vector RotVec(Monomers[imon].GetRotationVector());
    double RotAng = Monomers[imon].GetRotationAngle();
    for (int iA=0;iA<NAtoms;iA++)  {// loop over the atoms on this monomer 
      Multipole QA(Monomers[imon].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();

      Vector V0(NmomA); // stores the potential

      for (int iB=0;iB<NAtoms;iB++)  {// second loop over the atoms on this monomer 
	Multipole QB(Monomers[imon].GetAtom(iB).GetMultipoleMoments()); 
	int NmomB = QB.GetLength();
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imon].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }
	
	// if iA==iB: self-energy correction.  
	// Actually this charge-only self-energy term doesn't
	// contribute to the Madelung potential for induction.  Could skip
	// it, but it's trivial and leave it in case this routine ever
	// gets adapted to a full Madelung potential code.
	if (iA==iB){
	  V0[0] -= 2.0*ssff*QB(0);
	} 
	// if iA!=iB: intramolecular correction
	else {
	  // Want the standard Tab matrix for intramolecular correction.
	  Matrix Tab = Monomers[imon].GetAtom(iA).BuildInteractionMatrix(RotVec, RotAng,
									 Monomers[imon].GetAtom(iB), 
									 RotVec, RotAng,-999.0);

	  // increment the potential.  It's subtracted because in the end we are canceling
	  // the inclusion of intramolecular contributions from the earlier Ewald terms.
	  V0 -= Tab.MatrixTimesVector(QB.GetMoments()); 
	} 
      } // end loop over atoms iB

      // Now store contribution to Vmad
      int indexA = Atom_count[imon] + iA;
      int offsetA = offset[indexA];
      int nelem = sizes[indexA];
      //printf("Storing %d elements of Vmad(self/intra) for mon %d, atom %d, at %d\n",nelem,imon,iA,offsetA);
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += V0[t+1];
      }     
    } // end loop over atoms iA
  } // end loop over monomers in the central cell

    
  stop_self_intra_time = time(NULL);
  double self_intra_time = difftime(stop_self_intra_time, start_self_intra_time);
  //printf("  Time for self-interaction and intramolecular energy correction = %.0f seconds\n", self_intra_time);


  // (6) Compute correction due to induction damping.  We use undamped
  // Tab terms in the Ewald sum, but for induction, we want to correct
  // for short-range damping.  Since this effect is purely a
  // short-range one, we can evaluate the correction using
  // conventional Tab and DampedTab matrices in a finite cluster.
  // This term only involves induced moments, not permanent ones.  The
  // correction looks like 
  // V = V(damped) - V(undamped) = (DampedTab - Tab)*QB

  start_damping_time = time(NULL);

  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      int NmomA = Monomers[imonA].GetAtom(iA).GetMultipoleMoments().GetLength(); 
      Vector damped_V0(NmomA);
      
      if (! induced_moments_only) { // Case where we are treating permanent multipoles
	
	// Loop over monomers in the central unit cell
	for (int imonB = 1; imonB <= NMon; imonB++) {
	  if (imonA != imonB) {
	    int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	    
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
	    
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	      
	      Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());       
	      // Grab Tab for the current pair of atoms
	      Matrix Tab, damped_Tab;
	      if (AB_order) {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
		damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iA,iB));
	      }
	      else {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
		Tab.Transpose();
		
		damped_Tab.Initialize(Dimers[idimer].GetDampedTabInteractionMatrix(iB,iA));
		damped_Tab.Transpose();
	      }
	      
	      // Compute (damped_Tab - Tab)*QB
	      damped_Tab -= Tab;
	      damped_V0 += damped_Tab.MatrixTimesVector(QB.GetMoments()); 
	      
	    }  // end loop over atom iB
	  } // end if imonA != imonB
	} // end loop over monomer imonB
	
      }
      else {
	// Handle the induced moments part.  Since we only want
	// Madelung potential terms from the induced moments outside
	// the central unit cell, we subtract out the undamped
	// pieces which get included in steps (1)-(3), but we
	// *don't* add in the corresponding damped ones.  In other
	// words, we remove these contributions entirely.
	
	// Loop over monomers in the central unit cell
	for (int imonB = 1; imonB <= NMon; imonB++) {
	  if (imonA != imonB) {
	    int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	    
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
	    
	    for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	      Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	      
	      // Grab Tab for the current pair of atoms
	      Matrix Tab;
	      if (AB_order) {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iA,iB));
	      }
	      else {
		Tab.Initialize(Dimers[idimer].GetTabInteractionMatrix(iB,iA));
		Tab.Transpose();
	      }
	      // Compute -Tab*dQB
	      Tab.Scale(-1.0);
	      damped_V0 += Tab.MatrixTimesVector(dQB.GetMoments()); 
	      
	      
	    }  // end loop over atom iB
	  } // end if imonA != imonB
	} // end loop over monomer imonB
      }
      
      // Now loop over periodic image monomer contributions.  We've
      // already summed the Tabs exploiting the periodicity of the QB,
      // so we can just contract with QB from central cell monomers.
      for (int imonB=1;imonB<=NMon;imonB++) {
	for (int iB=0;iB<Monomers[imonB].GetNumberOfAtoms();iB++) {
	  
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  if (induced_moments_only) {
	    QB.Set();
	    Multipole dQB(Monomers[imonB].GetAtom(iB).GetInduceMultipoleMoments());
	    QB += dQB;
	  }
	  
	  int indexAB = EwaldLookupAtomPair(imonA,iA,imonB,iB);
	  Matrix DeltaDampTab = DeltaDampTabs[indexAB];
	  
	  damped_V0 += DeltaDampTab.MatrixTimesVector(QB.GetMoments()); 
	}
      } // end loop over monomer imonB
      
      //damped_V0.Print("damped V0 after PBC image monomers");
      
      // Now store contribution to Vmad
      int indexA = Atom_count[imonA] + iA;
      int offsetA = offset[indexA];
      int nelem = sizes[indexA];
      for (int t=0;t<nelem;t++) {
	Vmad[t+offsetA] += damped_V0[t+1];
      }   
      
    } // end loop over atom iA
  } // end loop over monomer imonA
  
  stop_damping_time = time(NULL);
  double damping_time = difftime(stop_damping_time, start_damping_time);
  //printf("  Time for damping correction = %.0f seconds\n", damping_time);
  
  stop_total = time(NULL);
  double time_total = difftime(stop_total, start_total);
  //printf("  Time to evaluate Madelung potential once Tabs are available = %.0f seconds\n",time_total);

  return Vmad;
}



void Cluster::BuildEwaldInteractionMatrices() {

  printf("Computing Ewald Interaction matrices\n");

  // Allocate storage.  Need Tab matrix for all possible pairs of atoms.
  int Natoms = GetTotalNumberOfAtoms(); 
  ReciprocalTabs = new Matrix[Natoms*Natoms];
  DirectTabs = new Matrix[Natoms*Natoms];
  DeltaDampTabs = new Matrix[Natoms*Natoms];  

  // Also create key for locating correct


  time_t start_ewald_recip_time, stop_ewald_recip_time;
  time_t start_ewald_direc_time, stop_ewald_direc_time;
  time_t start_damping_time, stop_damping_time;
  



  // Determine the Ewald summation cutoffs in direct (nX/nY/nZ) and
  // reciprocal space (kX,kY,kZ) and the corresponding kappa value
  int nX,nY,nZ,kX,kY,kZ;
  double kappa_param = DetermineDirectAndRecipSpaceCutoffs(nX,nY,nZ,kX,kY,kZ);

  // *** Now begin evaluating the potential *** //

  // (1) Reciprocal space interaction matrix
  start_ewald_recip_time = time(NULL);
  printf("  Evaluate the reciprocal space interaction matrices\n");

  // grab information for building the RecipTab; and calculate the
  // convergence factor alphaa used in the self-energy correction,
  // which is the same as that in RecipTab and DirecTab
  double CellV = cell_volume*AngToBohr*AngToBohr*AngToBohr;
  
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
      int NmomA = QA.GetLength();
       
      // Loop over all monomers in unit cell for B
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
        
	// Loop over atoms on each monomer B
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B            
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  
	  Matrix tmpRecipTab(NmomA,NmomB);

	  // Loop over reciprocal lattice vectors
	  for (int kx = -kX; kx<=kX; kx++){
	    for (int ky = -kY; ky<=kY; ky++){
	      for (int kz = -kZ; kz<=kZ; kz++){                  
		if (kx*kx+ky*ky+kz*kz!=0){ //if |kn|1=0
		  tmpRecipTab += Monomers[imonA].GetAtom(iA).BuildRecipInteractionMatrix(RotVecA, RotAngA,
										      Monomers[imonB].GetAtom(iB), RotVecB, RotAngB,
										      kx, ky, kz, cell_volume,
										      RecipCellx, RecipCelly, RecipCellz,-999.0);
		} 
	      } // end loop over kz
	    } // end loop over ky
	  } // end loop over kx

	  // Now store this contribution
	  int indexAB = EwaldLookupAtomPair(imonA,iA,imonB,iB);
	  ReciprocalTabs[indexAB].Initialize(tmpRecipTab);

	} // end loop over atoms iB
      } // end loop over monomers B

    } // end loop over atoms iA
  } // end loop over monomers A


  stop_ewald_recip_time = time(NULL);
  double ewald_recip_time = difftime(stop_ewald_recip_time, start_ewald_recip_time);
  printf("  Time to evaluate reciprocal space interaction matrices = %.0f seconds\n",ewald_recip_time);


  // (3) Direct space contribution.  
  start_ewald_direc_time = time(NULL);
  //printf("  Evaluate the direct space contribution to the Madelung potential\n");
  fflush(stdout);

  Vector UnitCellx = unit_cell[0];
  Vector UnitCelly = unit_cell[1];
  Vector UnitCellz = unit_cell[2];

  for (int imonA = 1; imonA <= NMon; imonA++) {
    int NatomsA = Monomers[imonA].GetNumberOfAtoms();
    Vector RotVecA(Monomers[imonA].GetRotationVector());
    double RotAngA = Monomers[imonA].GetRotationAngle();    

    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on monomer A
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();
      Vector V0(NmomA);
      
      for (int imonB = 1; imonB <= NMon; imonB++) {
	int NatomsB = Monomers[imonB].GetNumberOfAtoms();
	Vector RotVecB(Monomers[imonB].GetRotationVector());
	double RotAngB = Monomers[imonB].GetRotationAngle();
	
	for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on monomer B
	  Multipole QB(Monomers[imonB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
          
	  Matrix tmpDirectTab(NmomA,NmomB);

	  // loop over periodic image cells
	  for (int nx = -nX; nx<=nX; nx++){
	    for (int ny = -nY; ny<=nY; ny++){
	      for (int nz = -nZ; nz<=nZ; nz++){
		
		if (!(nx*nx+ny*ny+nz*nz==0 && imonA==imonB && iA==iB)) { // |rAB-rn|!=0
		  Matrix tmpTab = Monomers[imonA].GetAtom(iA).BuildDirecInteractionMatrix(RotVecA, RotAngA,
											  Monomers[imonB].GetAtom(iB), 
											  RotVecB, RotAngB,  nx, ny, nz, cell_volume,
											  UnitCellx, UnitCelly, UnitCellz, -999.0);
		  tmpDirectTab.AddBlockToMatrix(tmpTab,0,0);
		}
	      } // end loop over nz
	    } // end loop over ny
	  } // end loop over nx

	  // Now store this contribution
	  int indexAB = EwaldLookupAtomPair(imonA,iA,imonB,iB);
	  DirectTabs[indexAB].Initialize(tmpDirectTab);


	} // end loop over atoms iB
      } // end loop over monomers B
    } // end loop over atoms iA
  } // end loop over monomers A
  
  stop_ewald_direc_time = time(NULL);
  double ewald_direc_time = difftime(stop_ewald_direc_time, start_ewald_direc_time);
  printf("  Time to evaluate direct space interaction matrices = %.1f seconds\n",ewald_direc_time);
  



  // Tab's for the damping correction.
  start_damping_time = time(NULL);
  //     First initialize the memory
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {
      Multipole QA(Monomers[imonA].GetAtom(iA).GetMultipoleMoments()); 
      int NmomA = QA.GetLength();

      for (int imonB=1;imonB<=NMon_images;imonB++) {
	int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {
	  Multipole QB(Monomers[ref_monB].GetAtom(iB).GetMultipoleMoments());
	  int NmomB = QB.GetLength();
	  
	  int indexAB = EwaldLookupAtomPair(imonA,iA,ref_monB,iB);
	  Matrix tmp(NmomA,NmomB);
	  DeltaDampTabs[indexAB] = tmp;
	}
      }
    }
  }


  //     Now sum up the matrices
  for (int imonA=1;imonA<=NMon;imonA++) {
    for (int iA=0;iA<Monomers[imonA].GetNumberOfAtoms();iA++) {

      for (int imonB=1;imonB<=NMon_images;imonB++) {
	
	int ref_monB = MonomerImages[imonB].GetReferenceMonomerIndex();
	Dimer tmp_dimer;
	tmp_dimer.Initialize(Monomers[imonA],MonomerImages[imonB]);
	tmp_dimer.BuildTabInteractionMatrices();
	tmp_dimer.BuildDampedTabInteractionMatrices();
	
	for (int iB=0;iB<Monomers[ref_monB].GetNumberOfAtoms();iB++) {
	  Matrix Tab = tmp_dimer.GetTabInteractionMatrix(iA,iB);
	  Matrix damped_Tab = tmp_dimer.GetDampedTabInteractionMatrix(iA,iB);
	  damped_Tab -= Tab;
	  
	  int indexAB = EwaldLookupAtomPair(imonA,iA,ref_monB,iB);
	  DeltaDampTabs[indexAB] += damped_Tab;
	}
      }
    }
  }
  
  stop_damping_time = time(NULL);
  double damping_time = difftime(stop_damping_time, start_damping_time);
  printf("  Time to evaluate Ewald damping correction interaction matrices = %.0f seconds\n",damping_time);
 


}

// Maps monomer number/atom number pair onto a 1-D index.
int Cluster::EwaldLookupAtomPair( int imonA, int iatomA, int imonB, int iatomB ) {

  // First create list of atom number where each monomer starts.
  int monomer_list[NMon+1];
  monomer_list[1] = 0; // first monomer starts at zero.
  for (int i=2;i<=NMon;i++) {
    monomer_list[i] = monomer_list[i-1] + Natoms_per_monomer[i-1];
    //printf("monomer_list[%d] = %d\n",i,monomer_list[i]);
  }

  // Now determine the index for each atom in the A/B pair....
  int indexA = monomer_list[imonA] + iatomA;
  int indexB = monomer_list[imonB] + iatomB;
  // and the joint index.  
  int Natoms = GetTotalNumberOfAtoms(); 
  int indexAB = indexA*Natoms + indexB;
  //printf("Atom %d.%d indexA = %d, Atom %d.%d indexB = %d, final index = %d\n",imonA,iatomA,indexA,imonB,iatomB,indexB,index);

  return indexAB;
}
