#ifndef _dlf_interface_h
#define _dlf_interface_h

#include <iostream>
#include <stdio.h>
#include <fstream.h>
#include <cassert>
#include <iomanip.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "monomer.h"
#include "dimer.h"
#include "cluster.h"
#include "vector.h"
#include "opt.h"
#include "main.h"

#include "constants.h"
using namespace hmbi_constants;

// Library routines for DL-FIND
//
// nvarin:  3*Natoms
// nvarin2: repeat for second array (coords2)
// nspec:   size of int array spec
// master:  1=master node in parallel calc, 0 otherwise
//
extern "C" void dl_find_(int *nvarin, int *nvarin2, int *nspec, int *master);

// These next functions are C++ functions that get called by DL-FIND.
// To prevent C++ from mangling the function names at compile time, we
// have to treat them as 'extern "C"' as well.
#ifdef __cplusplus 
extern "C" {
#endif

  void dlf_get_params_(int *n, int *nvar2, int *nsp, double coords[], 
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
		       int *tdlf_farm, int *n_po_scaling);

  void dlf_get_gradient_(int *nvarin, double *coords, double *energy, 
			 double *gradient, int *iimage, int *status);

  void dlf_get_hessian_(int *nvar, double *coords, double *hessian, 
			int *status);

  void dlf_get_multistate_gradients_(int *nvarin, double *coords, 
				     double *energy, double *gradient, 
				     double *coup, int *needcoupling, 
				     int *iimage, int *status);

  void dlf_put_coords_(int *nvar, int* mode, double* energy, double *coords, 
  		      int *iam);

  void dlf_error_();

  void dlf_update_();
#ifdef __cplusplus
}
#endif
// end DL-FIND library routines


#endif
