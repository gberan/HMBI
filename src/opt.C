#include "opt.h"


bool Optimize() {

  double conv_tol = 1.0e-5;
  printf("conv_tol = %f\n",conv_tol);
  
  // dummy return until we write something smarter
  return 1;
}


void SteepestDescent(Vector& Rout, Vector& Rin, Vector& Grad) {

  /*
    Input: 
    Rin - initial coords
    Grad - Gradient at Rin
    N - size of Rin


    Output:
    Rout - new coordinates

   */

  double stepsize = 0.5; // Fixed stepsize, for now
  
  // Step: R(n+1) = R(n) + stepsize*(-Grad(n))
  Rout = Grad;
  Rout.Scale(-1.0*stepsize);
  Rout += Rin;
    
}



void SteepestDescent(Vector& Rout, Vector& Rin, Vector& Grad, 
		     double stepsize) {

  /*
    Input: 
    Rin - initial coords
    Grad - Gradient at Rin
    N - size of Rin


    Output:
    Rout - new coordinates

   */
  
  // Step: R(n+1) = R(n) + stepsize*(-Grad(n))
  Rout = Grad;
  Rout.Scale(-1.0*stepsize);
  Rout += Rin;
    
}

// Returns Rout, and StepDir_out
void ConjugateGradients(Vector& Rout, Vector& StepDir_out, Vector& Rin, 
			Vector& Grad, Vector& Grad_old, Vector& StepDir_old,
			double stepsize, bool reset) {

  printf("Cycle  : Conjugate Gradients\n");

  // Very wasteful with memory.  Clean up.

  //Grad.Print("CG: Gradient");
  //Grad_old.Print("CG: Old Gradient");
  //StepDir_old.Print("CG: Old Step Direction");

  // Compute Polak-Ribiere (PR) beta: 
  //          beta_i = G_i'(G_i - G_{i-1}) / G_{i-1}' G_{i-1}

  double beta;

  if (!reset) {
    Vector dG(Grad);
    dG -= Grad_old;
    beta = Grad.DotProduct(dG);
    beta /= Grad_old.DotProduct(Grad_old);
    
    //printf("beta = %f\n",beta);
  }
  else {
    printf("Cycle : Reseting beta\n");
    beta = 0.0;
  }

  // Step direction: d_i = -g_i + beta_i * d_{i-1};
  Vector StepDir(StepDir_old);
  StepDir.Scale(beta);
  StepDir -= Grad;

  //StepDir.Print("CG: New Step Direction");
  StepDir_out = StepDir;

  StepDir.Scale(-1.0);

  // Take actual step.  Use Steepest Descent routine to do it.
  SteepestDescent(Rout, Rin, StepDir,stepsize);
  /*
  double stepsize = 0.5;
  Rout = StepDir;
  Rout.Scale(stepsize);
  Rout += Rin;
  */



			  
}
