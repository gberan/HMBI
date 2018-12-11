#include "dlf_interface.h"

/*
Interface for running DL-FIND optimizers

KN and GJB, 11/09

*/

void dlf_get_params_(int *nvar, int *nvar2, int *nspec, double coords[], 
		     double coords2[], int spec[], int *ierr, 
		     double *tolerance, int *printl, int *maxcycle,
		     int *maxene, int *tatoms, int *icoord, int *iopt,
		     int *iline, double *maxstep, double *scalestep,
		     int *lbfgs_mem, int *nimage, double *nebk, int *dump,
		     int *restart, int *nz_i, int *ncons_i, int *nconn_i,
		     int *update, int *maxupd, double *delta, double *soft,
		     int *inithessian, int *carthessian, int *tsrel, 
		     int *maxrot, double *tolrot, int *nframe, int *nmass,
		     int *nweight, double *timestep, double *fric0, 
		     double *fricfac, double *fricp, int *imultistate,
		     int *state_i, int *state_j, double *pf_c1,
		     double *pf_c2, double *gp_c3, double *gp_c4,
		     double *ln_t1, double *ln_t2, int *printfile,
		     double *tolerance_e, double *tolerance_max_g,
		     double *tolerance_rms_g, double *distort, int *massweight,
		     double *minstep, int *maxdump, int *task,
		     double *temperature, int *po_pop_size,
		     double *po_radius, double *po_contraction,
		     double *po_tolerance_r, double *po_tolerance_g,
		     int *po_distribution, int *po_maxcycle,
		     int *po_init_pop_size, int *po_reset,
		     double *po_mutation_rate, double *po_death_rate,
		     double *po_scalefac, int *po_nsave, int *ntasks,
		     int *tdlf_farm, int *n_po_scaling,
     		     double *neb_climb_test, double *neb_freeze_test, int *nzero,
    		     int *coupled_states, int *qtsflag,
     		     int *imicroiter, int *maxmicrocycle, 
                     double *init_tr, int *micro_esp_fit){
  //JLM should check the default settings for the bottom 3 lines of the call (not yet set)

  //int Natoms =Cluster::cluster().GetTotalNumberOfAtoms();
  int Natoms = Cluster::cluster().GetNumberOfUniqueAtoms();

  //Vector vec_coords( Cluster::cluster().GetCurrentCoordinates() ,true );
  //Vector vec_coords( Cluster::cluster().GetSymmetryUniqueCoordinates(), true);
  Vector vec_coords = Cluster::cluster().GetSymmetryUniqueCoordinates() ;
  //vec_coords.Print(" vec_coords inside get_params \n");  
  //*nvar = 3*Cluster::cluster().GetTotalNumberOfAtoms();

  if( Params::Parameters().UseFullQMOnly() )
    *nvar = 3*Cluster::cluster().GetNumberOfUniqueAtoms() + 9;
  else if( Params::Parameters().IsPeriodic() && (Params::Parameters().GetMMType() != 2))
    *nvar = 3*Cluster::cluster().GetNumberOfUniqueAtoms() + 6;
  else
    *nvar = 3*Cluster::cluster().GetNumberOfUniqueAtoms();
  *nvar2 = 0;

  // Create 'spec' array.  The first Natoms elements use integers to
  // label each atom based on its monomer---the "residues" in
  // DL-FIND.  The second Natoms elements give the atomic number of
  // each atom.  Additional elements defining constraints or
  // connectivity are possible, but I haven't looked into that yet.
  // GJB
  int ispec = 0;
  // Label the "residues"
  // If using hybrid delocalized (i.e. separate deloc for each monomer, and
  // cartesian in between) assign based on monomer number:
  if (Params::Parameters().GetCoordType() == 1 || Params::Parameters().GetCoordType() == 2) {
    for (int imon=1;imon<=Cluster::cluster().GetNumberOfMonomers();imon++) {
    //  int natoms = Cluster::cluster().GetNumberOfAtoms(imon);
    //  for (int j=0;j<natoms;j++) {
    //	spec[ispec] = imon;
    //  ispec++;
      int natoms = Cluster::cluster().GetMonomer(imon).GetNumberOfAtoms();
      for (int j=0;j<natoms;j++){
	if(Cluster::cluster().GetMonomer(imon).GetAtom(j).InAsymmetricUnit()){
	  spec[ispec] = imon;
	  ispec++;
	}
      }
    }
    //lattice parameters not on a fragment
    if ( Params::Parameters().UseFullQMOnly()) {
      spec[ispec] = 0;
      //      spec[ispec] = Cluster::cluster().GetNumberOfMonomers() + 1;
      ispec++;
      spec[ispec] = 0;
      //      spec[ispec] = Cluster::cluster().GetNumberOfMonomers() + 1;
      ispec++;
      spec[ispec] = 0;
      //      spec[ispec] = Cluster::cluster().GetNumberOfMonomers() + 1;
      ispec++;
    }
    else if ( Params::Parameters().IsPeriodic()  && Params::Parameters().GetMMType() != 2) {
      spec[ispec] = 0;
      //      spec[ispec] = Cluster::cluster().GetNumberOfMonomers() + 1;
      ispec++;
      spec[ispec] = 0;
      //      spec[ispec] = Cluster::cluster().GetNumberOfMonomers() + 1;
      ispec++;
    }
  }
  else {
    printf("No fragments\n");
    
    // Otherwise, all atoms are on the same residue:
    for (int i=0;i<Natoms;i++) {
      spec[ispec] = 0;
      ispec++;
    }
    //lattice parameters not on a fragment
    if (Params::Parameters().UseFullQMOnly()) {
      spec[ispec] = 1;
      ispec++;
      spec[ispec] = 1;
      ispec++;
      spec[ispec] = 1;
      ispec++;
    }
    else if ( Params::Parameters().IsPeriodic() && Params::Parameters().GetMMType() != 2) {
      spec[ispec] = 1;
      ispec++;
      spec[ispec] = 1;
      ispec++;
    }
  }

  // Set the atomic numbers.
  for (int imon=1;imon<=Cluster::cluster().GetNumberOfMonomers();imon++) {
    for(int iatom=0;iatom<Cluster::cluster().GetNumberOfAtoms(imon);iatom++){
      if(Cluster::cluster().GetMonomer(imon).GetAtom(iatom).InAsymmetricUnit()){//Only use atoms of the symmetrical unique monomers
	//spec[ispec] = Cluster::cluster().GetAtomicNumber(i);
	spec[ispec] = Cluster::cluster().GetMonomer(imon).GetAtomicNumber(iatom);
	ispec++;
      }
    }
  }

  if ( Params::Parameters().UseFullQMOnly() ) {  
    spec[ispec] = 0;
    ispec++;
    spec[ispec] = 0;
    ispec++;
    spec[ispec] = 0;
    ispec++;
  }
  else if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) {  
    spec[ispec] = 0;
    ispec++;
    spec[ispec] = 0;
    ispec++;
  }

  //  printf("ispec = %d   nspec = %d\n",ispec,*nspec);
  // for (int q=0;q<ispec;q++) 
  //  printf("spec[%d] = %d\n",q,spec[q]);

  // Get the coordinates.  In bohrs (which DL-FIND wants) 
  Vector tmp_coords = Cluster::cluster().GetSymmetryUniqueCoordinates();
  if ( Params::Parameters().UseFullQMOnly()) {
    //converting Angs to Bohr
    Vector tmp_vec;
    Vector new_coords(3*Natoms+9);
    for(int i=0;i<3*Natoms;i++) {
      new_coords[i] = tmp_coords[i];
    }
    for(int i=0;i<3;i++) {
      tmp_vec = Cluster::cluster().GetUnitCellVector(i);
      for(int j=0;j<3;j++) {
        new_coords[3*Natoms + i*3 + j] = tmp_vec[j];
      }
    }
    new_coords.Scale(AngToBohr);
    for(int i=0;i<3*Natoms+9;i++) {
      coords[i] = new_coords[i];
    }
  }
  else if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) {
      /*printf("coords before saving\n");
      for (int q=0;q<Natoms;q++){
        printf("atom %i: %f %f %f\n",q,tmp_coords[3*q],tmp_coords[3*q+1],tmp_coords[3*q+2]);
        fflush(stdout);
      }
      printf("a b c = %f %f %f \n", tmp_coords[3*Natoms],tmp_coords[3*Natoms+1],tmp_coords[3*Natoms+2]);
      printf("alpha beta gamma = %f %f %f \n\n", tmp_coords[3*Natoms+3],tmp_coords[3*Natoms+4],tmp_coords[3*Natoms+5]);
      fflush(stdout);*/

    //converting Angs to  Bohr and Lattice Angles into Radians
    for(int i=0;i<3*Natoms+3;i++)
      tmp_coords[i] *= AngToBohr;
    for(int i=0;i<3;i++)
      tmp_coords[3*Natoms+3+i] *= DegreesToRadians;
    for (int i=0;i<3*Natoms+6;i++) {
      coords[i] = tmp_coords[i];
    }

  }
  else{
    tmp_coords.Scale(AngToBohr);
    for (int i=0;i<3*Natoms;i++) {
      coords[i] = tmp_coords[i];
    }
  }

  /*
  printf("coords after saving\n");
  for (int q=0;q<Natoms;q++) 
     printf("atom %i: %f %f %f\n",q,coords[3*q],coords[3*q+1],coords[3*q+2]);
  printf("a b c = %f %f %f \n", coords[3*Natoms],coords[3*Natoms+1],coords[3*Natoms+2]);
  printf("alpha beta gamma = %f %f %f \n", coords[3*Natoms+3],coords[3*Natoms+4],coords[3*Natoms+5]);
  fflush(stdout);
  //exit(0);
  */

  *ierr = 0;
  *printl = 3;
  *maxcycle = Params::Parameters().GetMaxOptCycles();
  *maxene = 2*Params::Parameters().GetMaxOptCycles();
  *tatoms = 1;   // it is ivar and not tatoms

  //icoord: coordinate system type. 
  //0 = cartesian; 
  //1 = hdlc (hybrid delocalized internal coordinates)
  *icoord = Params::Parameters().GetCoordType();
  *iopt = Params::Parameters().GetOptType(); 
  *iline = Params::Parameters().GetLineSearchType(); 
  *maxstep = 1.0;
  *scalestep = 0.5;
  *lbfgs_mem = -1;// how many previous steps to store
  *nimage  = -1;
  *nebk = -1.0;
  *dump = -1;
  *restart = -1;
  *ncons_i = -1;
  *nconn_i = -1;
  *update = 3;
  *maxupd = 1;
  *delta = -1;
  *soft  = 1e20;
  *inithessian = -1;
  *carthessian = -1;
  *tsrel = false;
  *maxrot = -1;
  *tolrot = -1e20;
  *nframe = 1;  // 1 for standard opt jobs.  Higher for NEB, for example.

  if (Params::Parameters().UseFullQMOnly()) {
    *nmass = Cluster::cluster().GetNumberOfUniqueAtoms() + 3;
    *nweight = Cluster::cluster().GetNumberOfUniqueAtoms() + 3;
    *nz_i = Cluster::cluster().GetNumberOfUniqueAtoms() + 3;
    //JLM
  }
  else if (Params::Parameters().IsPeriodic()) {
    *nmass = Cluster::cluster().GetNumberOfUniqueAtoms() + 2;
    *nweight = Cluster::cluster().GetNumberOfUniqueAtoms() + 2;
    *nz_i = Cluster::cluster().GetNumberOfUniqueAtoms() + 2;
  }
  else{
    *nmass = Cluster::cluster().GetNumberOfUniqueAtoms();
    *nweight = Cluster::cluster().GetNumberOfUniqueAtoms();
    *nz_i = Cluster::cluster().GetNumberOfUniqueAtoms();
  }

  *timestep = -1;
  *fric0 = -1;
  *fricfac =-1;
  *fricp = -1;
  *imultistate =0;
  *state_i = 1;
  *state_j = 2;
  *pf_c1 = 5.0;
  *pf_c2 = 5.0;
  *gp_c3 = 1.0;
  *gp_c4 = 0.9;
  *ln_t1 = 1.0e-4;
  *ln_t2 = 1.0;
  *printfile = 3;
  *tolerance_e = 100.0e-8; // E_tol: 100x10^-8  q-chem value
  *tolerance = 300.0e-6; // Grad tol: qchem value = 300x10^-6
  *distort = 0.;
  *massweight = 0;
  *minstep = -1;
  *maxdump = -1;
  *task = -1;
  *temperature = -1;
  *po_pop_size = -1;
  *po_radius = -1;
  *po_contraction = -1;
  *po_tolerance_r = -1;
  *po_tolerance_g = -1; 
  *po_distribution = -1;
  *po_maxcycle = -1;
  *po_init_pop_size = -1;
  *po_reset = -1;
  *po_mutation_rate = -1;
  *po_death_rate = -1;
  *po_scalefac = -1;
  *po_nsave = -1;
  *ntasks = 1;
  *tdlf_farm = 1;
  *n_po_scaling = 0;
  *init_tr = Params::Parameters().GetStupidVariable();


  //yoni::Coverages for the max and rms gradient component
  double ave_sym =0.0;
  int mon_count = 0; //counts number of monomers with non-zero symmetry factors
  
  for(int imon=1;imon<=Cluster::cluster().GetNumberOfMonomers();imon++)
    if(Cluster::cluster().GetMonomerSymmetryFactor(imon)!=0){
      ave_sym += Cluster::cluster().GetMonomerSymmetryFactor(imon);
      mon_count++;
    }
  ave_sym /= mon_count;

  //JLM
  //if(!Params::Parameters().UseSpaceSymmetry()){
    ave_sym = 1;
    mon_count = 1;
  //}
  //*tolerance_max_g = ave_sym * 4.5000e-4;
  //*tolerance_rms_g = ave_sym * 3.000e-4;

  /*Setting Convergence Watit
  if(Params::Parameters().GetConvergence()==0) {
    *tolerance_max_g = 1.4e-3;
    *tolerance_rms_g = 9.0e-4;
  }
  else if(Params::Parameters().GetConvergence()==2) {
    *tolerance_max_g = 1.5e-5;
    *tolerance_rms_g = 1.0e-5;
  }
  else if(Params::Parameters().GetConvergence()==3) {
    *tolerance_max_g = 2e-6;
    *tolerance_rms_g = 1e-6;
  }
  else {
    *tolerance_max_g = 4.5000e-4;
    *tolerance_rms_g = 3.000e-4;
  }*/
  //Set Convergence Cameron
  *tolerance_e = Params::Parameters().GetETol();
  *tolerance_max_g = Params::Parameters().GetMaxg();
  *tolerance_rms_g = Params::Parameters().GetMaxRms();
  printf("tolerance_max_g = %1.4e\n", *tolerance_max_g);
  printf("tolerance_rms_g = %1.4e\n", *tolerance_rms_g);
  printf("tolerance_e = %1.4e\n", *tolerance_e);

  //resetting the cycle counter
  Params::Parameters().ResetOptCycle();

  printf("Finished initializing DLFind\n");
  fflush(stdout);
}

void dlf_get_gradient_(int *nvarin, double *coords,double *energy, 
		       double *gradient,int *iimage, int *dlf_status) {

  *dlf_status=0;
  Params::Parameters().IncrementOptCycle();
  int Natoms = Cluster::cluster().GetNumberOfUniqueAtoms();

  if(Params::Parameters().UnFreezeParamsAfter() > 0 &&
     Params::Parameters().UnFreezeParamsAfter() <= Params::Parameters().GetOptCycle()){
    Params::Parameters().SetParameter("FREEZE_UNITCELLPARAMS","FALSE");
    printf("Params No longer frozen\n");
    printf("Cycle %i\n",Params::Parameters().GetOptCycle());
    printf("Unfreeze after %i\n",Params::Parameters().UnFreezeParamsAfter());
  }

  //printf("dlf_get_gradient\n");
  Vector vec_coords(coords,*nvarin);
  Vector all_coordinates;
  Vector tmp_vec;

  if (Params::Parameters().UseFullQMOnly() && Params::Parameters().GetQMType() == 5) { //Quantum Espresso
    //printf("Inside FullQMOnly \n");
    //fflush(stdout);

    vec_coords.Scale(BohrToAng);
    Vector old_coords = Cluster::cluster().GetLastAcceptedCoordinates();


    /*old_coords.PrintGradient("Last Accepted coordinates");
    printf("\n");
    fflush(stdout);

    vec_coords.PrintGradient("New coordinates before coupling");
    printf("\n");
    fflush(stdout);*/

    for(int i=0;i<old_coords.GetLength()-6;i++){
      old_coords[i] = vec_coords[i];
    }

    //printf("got old coords \n");
    //fflush(stdout);

    Matrix tmp(3,3);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++)
        tmp(i,j) = vec_coords[3*Natoms+3*i+j];

    //Need to convert from the 3 vector form to the 6 parameter form
    Vector lattice_param_new = Cluster::cluster().ComputeCellParametersFromLatticeVectors(tmp);

    //Enforcing symmetry
    Cluster::cluster().MaintainCartesianSymmetry(old_coords,true);

    //printf("Turning old_coords into fractional coords\n\n");
    //fflush(stdout);
    
    //Turn to fractional to get lattice gradient coupling
    old_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(old_coords,false,true);

    //replace old lattice with new one
    for(int i=0;i<6;i++){
      old_coords[3*Natoms+i] = lattice_param_new[i];
    }

    Cluster::cluster().UpdateLatticeParams(old_coords);

    //When including lattice symmetry, need to ensure symmetry hasn't locked certain parameters together.
    //If so we need to update the lattice parameters again...
    //First update lattice params again
    Vector axes = Cluster::cluster().GetUnitCellAxes();
    Vector angles = Cluster::cluster().GetUnitCellAngles();

    //replace old lattice with new one
    for(int i=0;i<3;i++){
      old_coords[3*Natoms+i] = axes[i];
      old_coords[3*Natoms+3+i] = angles[i];
    }

    //Also need to update the lattice vectors
    Vector tmp2(3);
    for(int i=0;i<3;i++){
      tmp2 = Cluster::cluster().GetUnitCellVector(i);
      for (int j=0; j<3;j++){
        vec_coords[3*Natoms+3*i+j] = tmp2[j];
      }
    }

  //Gives all the coordinates for the central unit cell
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(old_coords,true);

    //Transform back to cartesian coords
    all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,true,true);
    old_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(old_coords,true,true);

    //Need to give DLFIND the coupled coordinates
    Vector dlfind_coords= vec_coords;
    for (int i=0;i<3*Natoms;i++)
      dlfind_coords[i] = old_coords[i];

    dlfind_coords.Scale(AngToBohr);
    for (int i=0;i<3*Natoms+9;i++) {
      coords[i] = dlfind_coords[i];
    }

    //Gives all the coordinates for the central unit cell
    /*all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(old_coords,true);

    //Transform back to cartesian coords

    //This gives all the atoms in the central unit cell
    all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,true,true);*/

    /*all_coordinates.PrintGradient("New coordinates after coupling");
    printf("\n");
    fflush(stdout);*/

    //Need a new list that includes only the atoms in the assymetric unit cell
    //get the Coord_key for the Coordinates
    /*int coord_key[Natoms];
    int j = 0;
    int NMon = Cluster::cluster().GetNumberOfMonomers();
    for(int iMon=1;iMon<=NMon;iMon++){
      //if(Monomers[iMon].GetSymmetryFactor()!=0)
      for(int iatom=0; iatom<Cluster::cluster().GetNumberOfAtoms(iMon);iatom++){
	if(Cluster::cluster().GetMonomer(iMon).GetAtom(iatom).InAsymmetricUnit()){
	  coord_key[j] = Cluster::cluster().GetMonomer(iMon).GetAtom(iatom).GetGlobalIndex();
	  j++;
	}
      }
    }*/
    

    /*printf("coord_key =[");
    for(int i=0;i<Natoms;i++)
      printf(" %i",coord_key[i]);
    printf("]\n");
    fflush(stdout);*/
    //exit(0);
    //old_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(old_coords,true,true);
    //Update only the symmetrically unique atoms
    /*for (int i = 0; i < Natoms; i++){
      for( int k = 0; k<3; k++){
        //grad[3*i + k] = Grad_QM[3*(gradkey[i]-1)+k];

        vec_coords[3*i + k] = all_coordinates[3*(coord_key[i]-1)+k];
        //printf("vec_coords[%i] = %f\n",3*i + k,vec_coords[3*i + k]);
        //fflush(stdout);
      }
    }

    //Need to give DLFIND the coupled coordinates
    Vector dlfind_coords= vec_coords;
    dlfind_coords.Scale(AngToBohr);
    for (int i=0;i<3*Natoms+9;i++) {
      coords[i] = dlfind_coords[i];
    }*/

  }
  else if (Params::Parameters().UseFullQMOnly() && Params::Parameters().GetQMType() == 8) { //DFTB+
    //printf("Inside FullQMOnly \n");
    //fflush(stdout);

    vec_coords.Scale(BohrToAng);
    Vector old_coords = Cluster::cluster().GetLastAcceptedCoordinates();


    for(int i=0;i<old_coords.GetLength()-6;i++){
      old_coords[i] = vec_coords[i];
    }

    //printf("got old coords \n");
    //fflush(stdout);

    Matrix tmp(3,3);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++)
        tmp(i,j) = vec_coords[3*Natoms+3*i+j];

    //Need to convert from the 3 vector form to the 6 parameter form
    Vector lattice_param_new = Cluster::cluster().ComputeCellParametersFromLatticeVectors(tmp);

    for(int i=0; i<6; i++){
      old_coords[old_coords.GetLength()-6+i] = lattice_param_new[i];
    }

    //Correcting coordinate locked under symmetry 
    Cluster::cluster().MaintainCartesianSymmetry(old_coords,true); 
  
    //Set the new lattice parameters
    Cluster::cluster().UpdateLatticeParams(old_coords);

    //Getting the coordinates for the rest of the system based on symmetry

    old_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(old_coords,false,true);
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(old_coords,true);
    all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,true,true);

  }
  else if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)){
    for(int i=0;i<vec_coords.GetLength()-3;i++)
      vec_coords[i] *= BohrToAng;
    for(int i=vec_coords.GetLength()-3; i < vec_coords.GetLength(); i++)
      vec_coords[i] *= RadiansToDegrees;

    //Correcting coordinate locked under symmetry 
    Cluster::cluster().MaintainCartesianSymmetry(vec_coords,true); 
  
    //Set the new lattice parameters
    Cluster::cluster().UpdateLatticeParams(vec_coords);

    //Getting the coordinates for the rest of the system based on symmetry

    vec_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(vec_coords,false,true);
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,true);
    all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,true,true);
  }
  else{
    vec_coords.Scale(BohrToAng); 
    Cluster::cluster().MaintainCartesianSymmetry(vec_coords,false);
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,false);
  } 
  //vec_coords.PrintGradient("vec_coords in dlf_get_gradient_() after converting units");
  //uses the coordinates of the symmetrical unique atoms to get the coordinates of all the atoms in the unit cell.

  Cluster::cluster().SetNewCoordinates(all_coordinates);

  //Cluster::cluster().PrintGradient("New coordinates",all_coordinates);  
  Cluster::cluster().RunJobsAndComputeEnergy();
  *energy = Cluster::cluster().GetHMBIEnergy();
  Vector tmp_grad(*nvarin);
  if ( Params::Parameters().UseFiniteDifferenceGradients() ) {
    tmp_grad = GetFiniteDifferenceGradient(vec_coords);
    tmp_grad.PrintGradient("Finite Difference Gradient");
  }
  else {
    tmp_grad = Cluster::cluster().GetHMBIGradient();
    //tmp_grad.PrintGradient("Analytical Gradient");
  }

  for (int i=0;i<*nvarin;i++) {
    gradient[i] = tmp_grad[i];
  }

  Vector vec_gradient(gradient,*nvarin);
  //vec_gradient.PrintGradient("vec_gradient in get_gradient");

  // Print out some information
  double Gnorm = tmp_grad.RMS();
  printf("\nCycle %d: Energy = %15.9f   |Grad| = %10.6f\n\n",
	 Params::Parameters().GetOptCycle(),*energy,Gnorm);
  //Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
  fflush(stdout);

  //int Natoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  // printf("coords in get_gradient\n");
  // for (int q=0;q<Natoms;q++) 
  //    printf("atom %i: %f %f %f\n",q,coords[3*q],coords[3*q+1],coords[3*q+2]);
  //printf("a b c = %f %f %f \n", coords[3*Natoms],coords[3*Natoms+1],coords[3*Natoms+2]);
  //printf("alpha beta gamma = %f %f %f \n", coords[3*Natoms+3],coords[3*Natoms+4],coords[3*Natoms+5]);

  //print optimization step
  Cluster::cluster().MakeOptimizationStepFile();
  // create restart file
  FILE *restart;
  string restart_file = "restart.in";
  if ((restart = fopen(restart_file.c_str(),"w"))==NULL) {
    printf("dlf_put_coords_ : Cannot open file '%s'\n",restart_file.c_str());
    exit(1);
  } 
  Cluster::cluster().PrintInputFile(restart);
  printf("\nRestart file written to '%s'\n",restart_file.c_str());

  fclose(restart);

  //print fractional coordinates
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    FILE *fract;
    string fract_file = "fract_coord.frac";
    if ((fract = fopen(fract_file.c_str(),"w"))==NULL) {
      printf("dlf_put_coords_ : Cannot open file '%s'\n",fract_file.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintFractionalCoordinates(fract);
    fclose(fract);
  }

}

void dlf_get_hessian_(int *nvar, double *coords, double *hessian, 
		      int *dlf_status) {

  printf("dlf_get_hessian():: Hessian is not yet implemented\n");
  exit(1);

}

void dlf_get_multistate_gradients_(int *nvarin, double *coords, 
				   double *energy, double *gradient, 
				   double *coup, int *needcoupling, 
				   int *iimage, int *dlf_status) {

  printf("dlf_get_multistate_gradient(): This feature is not implemented\n");
  exit(1);
}

//Everything done here was moved to dlf_get_gradient_
void dlf_put_coords_(int *nvar, int* mode, double* energy, 
		     double *coords, int *iam) {


  /*  mode = 1: optimized coordinates
      mode = 2: Transition mode: coordinates displaced relative to 
                those in mode 1
      <0 : NEB image of number -mode

      iam is dummy argument for now
  */

  //printf("dlf_put_coords: Storing coordinates\n");
  /*fflush(stdout);

  Vector vec_coords(coords,*nvar);
  Vector all_coordinates;*/
  //Cluster::cluster().PrintGradient("New coordinates before conversion to angstrom",vec_coords);


  // Convert back to Angstroms and degrees
  //if ( ((Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly() ) {
  //if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2) ) {
    /*
    for(int i=0;i<vec_coords.GetLength()-3;i++)
      vec_coords[i] *= BohrToAng;
    for(int i=vec_coords.GetLength()-3; i < vec_coords.GetLength(); i++)
      vec_coords[i] *= RadiansToDegrees;
    */

    //vec_coords.PrintGradient("coords in dlf_put_coords_() before converting units");
    //printf("\n");
    //Convert to  Angstroms and Degrees
    /*for(int i=0;i<vec_coords.GetLength()-3;i++)
      vec_coords[i] *= BohrToAng;
    for(int i=0;i<3;i++)
      vec_coords[vec_coords.GetLength()-3+i] *=  RadiansToDegrees;

    //Correcting coordinate locked under symmetry 
    Cluster::cluster().MaintainCartesianSymmetry(vec_coords,true);
    //Vector UniqueMonomerCoord = Cluster::cluster().GetFullMonomerCoordinates(vec_coords,true);
    //exit(0);


    //vec_coords.PrintGradient("Coord in dlf_put_coords_");

    //Set the new lattice parameters
    Cluster::cluster().UpdateLatticeParams(vec_coords);

    //Getting the coordinates for the rest of the system base on symmetry
    vec_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(vec_coords,false,true);
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,true);
    all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,true,true);
    
    //uses the coordinates of the symmetrically unique atoms to get the fractional coordinates of all the atoms in the unit cell.
    //all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,1);
    //getting the cartesian coordinates
    //all_coordinates = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coordinates,1,1);

  }
  else{

    vec_coords.Scale(BohrToAng);
    Cluster::cluster().MaintainCartesianSymmetry(vec_coords,false);
    all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,0); 
  }


  //vec_coords.PrintGradient("coords in dlf_get_coords_() after converting units");


  Vector vec_coords(coords,*nvar);
  vec_coords.Scale(BohrToAng);

  //vec_coords.PrintGradient("Print 13: vec_coords scaled incorrectly");
  //Again correct for the incorrect scaling of the unitcell angles
  if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2) ) {
    for (int i=vec_coords.GetLength()-3; i < vec_coords.GetLength(); i++) {
      vec_coords[i] = vec_coords[i]/BohrToAng/DegreesToRadians; // correct for the previous error and also change radians to degrees
    }
  }
  */
/*
  Cluster::cluster().SetNewCoordinates(all_coordinates,false);

  //print optimization step
  Cluster::cluster().MakeOptimizationStepFile();
  // create restart file
  FILE *restart;
  string restart_file = "restart.in";
  if ((restart = fopen(restart_file.c_str(),"w"))==NULL) {
    printf("dlf_put_coords_ : Cannot open file '%s'\n",restart_file.c_str());
    exit(1);
  } 
  Cluster::cluster().PrintInputFile(restart);
  printf("\nRestart file written to '%s'\n",restart_file.c_str());

  fclose(restart);

  //print fractional coordinates
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    FILE *fract;
    string fract_file = "fract_coord.frac";
    if ((fract = fopen(fract_file.c_str(),"w"))==NULL) {
      printf("dlf_put_coords_ : Cannot open file '%s'\n",fract_file.c_str());
      exit(1);
    } 
    Cluster::cluster().PrintFractionalCoordinates(fract);
    fclose(fract);
  }

*/
  //int Natoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  //  printf("coords in put_coords\n");
  // for (int q=0;q<Natoms;q++) 
  //    printf("atom %i: %f %f %f\n",q,coords[3*q],coords[3*q+1],coords[3*q+2]);
      //printf("a b c = %f %f %f \n", coords[3*Natoms],coords[3*Natoms+1],coords[3*Natoms+2]);
  //printf("alpha beta gamma = %f %f %f \n", coords[3*Natoms+3],coords[3*Natoms+4],coords[3*Natoms+5]);

}
				    
void dlf_error_() {
  printf("DL-FIND has returned a fatal error.  Exiting...\n"); 
  exit(1);
}

void dlf_update_() {
  printf("dlf_update(): This function is unnecessary for our purposes.  It does nothing.\n");
}


//created by yoni to deallocate all parameters
void dlf_deallocate_(int *nvar, int *nvar2, int *nspec, double coords[], 
		     double coords2[], int spec[], int *ierr, 
		     double *tolerance, int *printl, int *maxcycle,
		     int *maxene, int *tatoms, int *icoord, int *iopt,
		     int *iline, double *maxstep, double *scalestep,
		     int *lbfgs_mem, int *nimage, double *nebk, int *dump,
		     int *restart, int *nz_i, int *ncons_i, int *nconn_i,
		     int *update, int *maxupd, double *delta, double *soft,
		     int *inithessian, int *carthessian, int *tsrel, 
		     int *maxrot, double *tolrot, int *nframe, int *nmass,
		     int *nweight, double *timestep, double *fric0, 
		     double *fricfac, double *fricp, int *imultistate,
		     int *state_i, int *state_j, double *pf_c1,
		     double *pf_c2, double *gp_c3, double *gp_c4,
		     double *ln_t1, double *ln_t2, int *printfile,
		     double *tolerance_e, double *tolerance_max_g,
		     double *tolerance_rms_g, double *distort, int *massweight,
		     double *minstep, int *maxdump, int *task,
		     double *temperature, int *po_pop_size,
		     double *po_radius, double *po_contraction,
		     double *po_tolerance_r, double *po_tolerance_g,
		     int *po_distribution, int *po_maxcycle,
		     int *po_init_pop_size, int *po_reset,
		     double *po_mutation_rate, double *po_death_rate,
		     double *po_scalefac, int *po_nsave, int *ntasks,
		     int *tdlf_farm, int *n_po_scaling){

  delete nvar;
  delete nvar2;
  delete nspec; 
  delete[] coords; 
  delete[] coords2;
  delete[] spec;
  delete ierr;
  delete tolerance; 
  delete printl;
  delete maxcycle;
  delete maxene;
  delete tatoms;
  delete icoord; 
  delete iopt;
  delete iline;
  delete maxstep;
  delete scalestep;
  delete lbfgs_mem; 
  delete nimage;
  delete nebk;
  delete dump;
  delete restart;
  delete nz_i;
  delete ncons_i;
  delete nconn_i;
  delete update;
  delete maxupd; 
  delete delta;
  delete soft;
  delete inithessian;
  delete carthessian;
  delete tsrel; 
  delete maxrot;
  delete tolrot;
  delete nframe;
  delete nmass;
  delete nweight;
  delete timestep;
  delete fric0; 
  delete fricfac;
  delete fricp; 
  delete imultistate;
  delete state_i;
  delete state_j;
  delete pf_c1;
  delete pf_c2;
  delete gp_c3;
  delete gp_c4;
  delete ln_t1; 
  delete ln_t2;
  delete printfile;
  delete tolerance_e;
  delete tolerance_max_g;
  delete tolerance_rms_g;
  delete distort;
  delete massweight;
  delete minstep; 
  delete maxdump;
  delete task;
  delete temperature;
  delete po_pop_size;
  delete po_radius;
  delete po_contraction;
  delete po_tolerance_r;
  delete po_tolerance_g;
  delete po_distribution;
  delete po_maxcycle;
  delete po_init_pop_size;
  delete po_reset;
  delete po_mutation_rate;
  delete po_death_rate;
  delete po_scalefac;
  delete po_nsave; 
  delete ntasks;
  delete tdlf_farm;
  delete n_po_scaling;
}

