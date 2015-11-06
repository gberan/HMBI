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
		     double *tolerance_e, double *distort, int *massweight,
		     double *minstep, int *maxdump, int *task,
		     double *temperature, int *po_pop_size,
		     double *po_radius, double *po_contraction,
		     double *po_tolerance_r, double *po_tolerance_g,
		     int *po_distribution, int *po_maxcycle,
		     int *po_init_pop_size, int *po_reset,
		     double *po_mutation_rate, double *po_death_rate,
		     double *po_scalefac, int *po_nsave, int *ntasks,
		     int *tdlf_farm, int *n_po_scaling) {
  
  int Natoms =Cluster::cluster().GetTotalNumberOfAtoms();

  Vector vec_coords( Cluster::cluster().GetCurrentCoordinates() ,true );
  //vec_coords.Print(" vec_coords inside get_params \n");  
  *nvar = 3*Cluster::cluster().GetTotalNumberOfAtoms();
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
      int natoms = Cluster::cluster().GetNumberOfAtoms(imon);
      for (int j=0;j<natoms;j++) {
	spec[ispec] = imon;
      ispec++;
      }
    }
  }
  else {
    // Otherwise, all atoms are on the same residue:
    for (int i=0;i<Natoms;i++) {
      spec[ispec] = 1;
      ispec++;
    }
  }
  // Set the atomic numbers.
  for (int i=0;i<Natoms;i++) {
    spec[ispec] = Cluster::cluster().GetAtomicNumber(i);
    ispec++;
  }
  /*
  printf("ispec = %d   nspec = %d\n",ispec,*nspec);
  for (int q=0;q<ispec-1;q++) 
    printf("spec[%d] = %d\n",q,spec[q]);
  */

  // Get the coordinates.  Convert them to bohr (which DL-FIND wants).
  Vector tmp_coords = Cluster::cluster().GetCurrentCoordinates();
  tmp_coords.Scale(AngToBohr);
  for (int i=0;i<3*Natoms;i++) {
    coords[i] = tmp_coords[i];
  }

  *ierr = 0;
  //*tolerance = 300.0e-6; // Grad tol: qchem value = 300x10^-6
  *tolerance = -1.0; // default DL-find value
  *printl = 3;
  *maxcycle = Params::Parameters().GetMaxOptCycles();
  *maxene = 2*Params::Parameters().GetMaxOptCycles();
  *tatoms = 1;   // it is ivar and not tatoms
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
  *nz_i = Cluster::cluster().GetTotalNumberOfAtoms();
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
  *nmass = Cluster::cluster().GetTotalNumberOfAtoms();
  *nweight = Cluster::cluster().GetTotalNumberOfAtoms();
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
  //*tolerance_e = 100.0e-8; // E_tol: 100x10^-8  q-chem value
  *tolerance_e = -1.0; // default dl-find value
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
}

void dlf_get_gradient_(int *nvarin, double *coords, double *energy, 
		       double *gradient, int *iimage, int *dlf_status) {

  *dlf_status=0;
  Params::Parameters().IncrementOptCycle();

  //printf("dlf_get_gradient\n");
  Vector vec_gradient;
  Vector vec_coords(coords,*nvarin);	
  // Convert back to Angstroms
  vec_coords.Scale(BohrToAng);
  //Cluster::cluster().PrintGradient("New coordinates",vec_coords);    

  Cluster::cluster().SetNewCoordinates(vec_coords);
  Cluster::cluster().RunJobsAndComputeEnergy();
  *energy = Cluster::cluster().GetHMBIEnergy();
  fflush(stdout);

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




  // Print out some information
  double Gnorm = tmp_grad.RMS();
  printf("\nCycle %d: Energy = %15.9f   |Grad| = %10.6f\n\n",
	 Params::Parameters().GetOptCycle(),*energy,Gnorm);
  //Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
  fflush(stdout);
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


void dlf_put_coords_(int *nvar, int* mode, double* energy, 
		     double *coords, int *iam) {

  /*  mode = 1: optimized coordinates
      mode = 2: Transition mode: coordinates displaced relative to 
                those in mode 1
      <0 : NEB image of number -mode

      iam is dummy argument for now
  */

  //printf("dlf_put_coords: Storing coordinates\n");

  Vector vec_coords(coords,*nvar);
  vec_coords.Scale(BohrToAng);
  Cluster::cluster().SetNewCoordinates(vec_coords,false);

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

}
				    
void dlf_error_() {
  printf("DL-FIND has returned a fatal error.  Exiting...\n"); 
  exit(1);
}

void dlf_update_() {
  printf("dlf_update(): This function is unnecessary for our purposes.  It does nothing.\n");
}

