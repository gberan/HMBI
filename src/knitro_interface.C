#ifdef KNITRO
#include "knitro_interface.h"

/*
Interface for running K-NITRO optimizers

GJB 7/11

*/

int knitro_callback_energy_and_gradient(const int evalRequestCode, 
					 const int n,
					 const int m, 
					 const int nnzJ, 
					 const int nnzH, 
					 const double * const x, 
					 const double * lambda,
					 double * energy,
					 double * c,
					 double * gradient,
					 double * jac,
					 double * hessian,
					 double * hessVector,
					 void * userParams) {

  if (evalRequestCode == KTR_RC_EVALFC) {
    
    Params::Parameters().IncrementOptCycle();
    
    //printf("Knitro: Getting energy & gradient\n");
    Vector vec_gradient;
    Vector vec_coords(n);	
    
    // copy coordinates to a vector for use in Cluster
    for (int i=0;i<n;i++) {
      vec_coords[i] = x[i];
    }

    vec_coords.Print("Nuclear coords returned to code:");

    // Convert back to Angstroms
    vec_coords.Scale(BohrToAng);
    Cluster::cluster().PrintGradient("New coordinates",vec_coords);    
    
    Cluster::cluster().SetNewCoordinates(vec_coords);
    Cluster::cluster().RunJobsAndComputeEnergy();
    *energy = Cluster::cluster().GetHMBIEnergy();
    fflush(stdout);
    
    Vector tmp_grad(n);
    if ( Params::Parameters().UseFiniteDifferenceGradients() ) {
      tmp_grad = GetFiniteDifferenceGradient(vec_coords);
      tmp_grad.PrintGradient("Finite Difference Gradient");
    }
    else {
      tmp_grad = Cluster::cluster().GetHMBIGradient();
      //tmp_grad.PrintGradient("Analytical Gradient");
    }
    
    for (int i=0;i<n;i++) {
      gradient[i] = tmp_grad[i];
    }
    
    // Print out some information
    double Gnorm = tmp_grad.RMS();
    printf("\nCycle %d: Energy = %15.9f   |Grad| = %10.6f\n\n",
	   Params::Parameters().GetOptCycle(),*energy,Gnorm);
    Cluster::cluster().UpdateTrajectoryFile(Params::Parameters().GetOptCycle());
    fflush(stdout);
    return(0);
    
  }
  else if (evalRequestCode == KTR_RC_EVALGA) {
    //printf("Knitro requested gradient\n");
    return(0);
  }
  else{
    printf ("*** callbackEvalFCorGA incorrectly called with eval code %d\n",
	    evalRequestCode);
    return( -1 );
  }
}

int knitro_callback_get_hessian(const int evalRequestCode, 
					 const int n,
					 const int m, 
					 const int nnzJ, 
					 const int nnzH, 
					 const double * const x, 
					 const double * lambda,
					 double * energy,
					 double * c,
					 double * gradient,
					 double * jac,
					 double * hessian,
					 double * hessVector,
					 void * userParams){

  printf("knitro_callback_get_hessian():: Analytical hessian unavailable\n");

  if (evalRequestCode == KTR_RC_EVALH || evalRequestCode == KTR_RC_EVALHV) {
    return(0);
  }
  else{
    printf ("*** callbackEvalHess incorrectly called with eval code %d\n",
	    evalRequestCode);
    return( -1 );
  }

}

#endif
