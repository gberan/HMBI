#include "opt.h"

// Constructor sets default parameters
Opt::Opt() {
  // Initialize geometry optimization parameters
  opt_cycle = 1;
  Natoms = Cluster::cluster().GetNumberOfUniqueAtoms();
  dE = 1000, Gnorm = 1000; //nonsense initial values
  Gmax = 1000, StepNorm = 1000, StepMax = 1000; // Newer optimizer uses 5 convergence criteria for tighter (reliable??) convergence (KDN)
  reduced_Gnorm = 1000, reduced_Gmax = 1000, constraint_norm = 1000;
  e_current = 0.0, e_old = 0.0; 
  pred = 0.0, tr_max = 1.25, tr_min = 0.00001, tr_init = 0.75;
  numStepAccepted=numStepsRejected=0;

  //econv = 1.0e-6;
  //gconv = 3.0e-4;
  econv = 100.0e-8;
  gconv = 300.0e-6; //rms gradient convergence criterion
  Gmax_conv = 4.5e-4;
  StepMax_conv = 1.8e-3;
  StepNorm_conv = 1.2e-3;
  stepsize = 0.25;
  stepsize_new = stepsize_old = Params::Parameters().GetInitialStepSize();

  reduced_Gconv = 8.4e-5;
  reduced_Gmax_conv = 8.4e-5;
  constraint_conv = 6.4e-5;
  reset=false; // for resetting CG algorithm
  armijoConditionMet = false;
  wolfePassed = false;
  lbfgs_sd_step = false;
  type_opt = "";
  
  // Initialize some empty vectors with same dimensions as the Gradient
  stepDir = Cluster::cluster().GetHMBIGradient();
  grad_current = Cluster::cluster().GetHMBIGradient();
  stepDir_old = Cluster::cluster().GetHMBIGradient();
  grad_old = Cluster::cluster().GetHMBIGradient();
  if ( (Params::Parameters().IsPeriodic()  && (Params::Parameters().GetMMType() != 2)) || Params::Parameters().UseFullQMOnly()){
    coords_new = Cluster::cluster().GetSymmetryUniqueCoordinates(true);
    coords_old = Cluster::cluster().GetSymmetryUniqueCoordinates(true);
    for(int i=coords_new.GetLength()-3; i < coords_new.GetLength(); i++){
      coords_new[i] *= DegreesToRadians;
      coords_old[i] *= DegreesToRadians;
    }
  }
  else{
    coords_new = Cluster::cluster().GetSymmetryUniqueCoordinates(false);
    coords_old = Cluster::cluster().GetSymmetryUniqueCoordinates(false);
  }

}

Opt::~Opt() {

}

/*bool Optimize::Optimize() {

  double conv_tol = 1.0e-5;
  printf("conv_tol = %f\n",conv_tol);
  
  // dummy return until we write something smarter
  return 1;
}*/

void Opt::OptimizeGeometry(string type) {

  //printf("passed string type: %s\n",type.c_str());
  //exit(0);
  string infile = Cluster::cluster().GetInputFilename();

  double ave_sym =0.0;
  int mon_count = 0; //counts number of monomers with non-zero symmetry factors
  
  if (Params::Parameters().UseSpaceSymmetry()) {
    for(int imon=1;imon<=Cluster::cluster().GetNumberOfMonomers();imon++)
      if(Cluster::cluster().GetMonomerSymmetryFactor(imon)!=0){
        ave_sym += Cluster::cluster().GetMonomerSymmetryFactor(imon);
        mon_count++;
      }

    ave_sym /= mon_count;
  
    printf("ave_sym = %f mon_count = %i\n",ave_sym,mon_count);
    Gmax_conv *= ave_sym;
    printf("Gmax_conv = %f\n",Gmax_conv);
    gconv *= ave_sym;
    printf("tolerance_rms_g = %f\n",gconv);
  }
  
  if (type == "SteepestDescent") {
    type_opt = "SD";
    SteepestDescent();
  }
  else if (type == "ConjugateGradients") {
    type_opt = "CG";
    ConjugateGradients();
  }
  else if (type == "LBFGS") {
    type_opt = "LBFGS";
    LBFGS();
  }
  else {
    printf("string type not recognized. type = %s. Exiting...\n",type.c_str());
    exit(1);
  }

}

void Opt::SteepestDescent() {

  //double stepsize_old = 0.65;
  //double stepsize_new = 0.65;
  double g_min, tmp;
  Vector all_coords;
  printf("Doing SteepestDescent\n");
  //printf("opt_cycle = %i dE = %f Gnorm = %f\n",opt_cycle,dE,Gnorm);
  //printf("maxcycle = %i, econv =%f, gconv = %f\n",Params::Parameters().GetMaxOptCycles(),econv,gconv);

  while (opt_cycle <= Params::Parameters().GetMaxOptCycles() && 
	   (fabs(dE) > econv  || fabs(Gnorm) > gconv  || fabs(Gmax) > Gmax_conv) ){


    //grad_current.PrintGradient("grad_current");
    //printf("\n");
    //grad_old.PrintGradient("grad_old");
    //printf("\n");

    coords_new.PrintGradient("coords_new");
    printf("\n");
    coords_old.PrintGradient("coords_old");
    printf("\n");

    //make next step
    // Step: R(n+1) = R(n) + stepsize*(-Grad(n))
    if (opt_cycle > 1) {

      if (opt_cycle > 2 ) {
        stepsize_new = CheckStep(stepsize_old, stepDir_old);
        printf("stepsize_old = %f, stepsize_new = %f\n",stepsize_old, stepsize_new);
        printf("\n");
        coords_new.PrintGradient("coords_new after check");
        printf("\n");
      }

      //coords_old.PrintGradient("coords_old");
      printf("\n");
      stepDir = grad_current;
      stepDirNorm = stepDir.Norm();
      if(opt_cycle == 2)
        stepsize_new *= stepDirNorm;
      stepDir.Scale(-1.0*stepsize_new/stepDirNorm);
      coords_new.PrintGradient("coords_new grad scaled");
      coords_new += stepDir;
      stepDir.Scale(1.0/stepsize_new);
      printf("\n");
      coords_new.PrintGradient("coords_new");
      printf("\n");

      tmp_coords = coords_new;

      /*if (Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
        CoupleGradients();
      }
      else if (Params::Parameters().IsPeriodic()) {
        Cluster::cluster().UpdateLatticeParams(coords_new);
      }*/

      if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
        for(int i=tmp_coords.GetLength()-3; i < tmp_coords.GetLength(); i++)
          tmp_coords[i] *= RadiansToDegrees;
      }
      
      if (Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
        CoupleGradients();
        all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
      }
      else if (Params::Parameters().IsPeriodic()) {
        Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,true);
        Cluster::cluster().UpdateLatticeParams(tmp_coords);
        tmp_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(tmp_coords,false,true);
        all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
        all_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coords,true,true);
      }
      else {
        Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,false);
        all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,false);
      }

      /*if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
        for(int i=all_coords.GetLength()-3; i < all_coords.GetLength(); i++)
          all_coords[i] *= DegreesToRadians;
      }*/

      //coords_new = Cluster::cluster().GetSymmetryUniqueCoordinates(true);
      
      //printf("\n");
      //coords_new.PrintGradient("coords_new");
      //printf("\n");

      //Cluster::cluster().MaintainCartesianSymmetry(vec_coords,false);
      //all_coordinates = Cluster::cluster().GetSymmetryImposedCoordinates(vec_coords,false);

      // Update the coordinates in the cluster object
      //Cluster::cluster().SetNewCoordinates(coords_new);
      Cluster::cluster().SetNewCoordinates(all_coords);

    }
    //grad_current.PrintGradient("Actual gradient");


    e_old = e_current;
    grad_old = grad_current;
    stepDir_old = stepDir;
    stepDirNorm_old = stepDirNorm;
    //coords_old = coords_new; 
    if(opt_cycle > 2){
      stepsize_old = stepsize_new;
    }
    //else if(opt_cycle <= 2)
    //  coords_old = coords_new;

    // Create the new jobs, run them, and get the HMBI energy
    Cluster::cluster().RunJobsAndComputeEnergy();
    e_current = Cluster::cluster().GetHMBIEnergy();
    grad_current = Cluster::cluster().GetHMBIGradient();

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      for(int i=grad_current.GetLength()-3; i < grad_current.GetLength(); i++)
        grad_current[i] *= DegreesToRadians;
    }

    //update optimization exit params
    dE = e_current - e_old;
    Gmax = grad_current.Max(true);
    //g_min = grad_current.Min(true);
    Gnorm=grad_current.RMS();
    //dGrad -= grad_current;
    //double dG = dGrad.Max(true);  
    //don't use dE for convergence for first step 
    if(opt_cycle == 1){
      dE = 0.0;
    }

    printf("e_old = %15.12f, e_current = %15.12f, stepsize_old = %15.12f\n", e_old, e_current, stepsize_old); 
    //need to think of a smart way to vary stepsize
    /*if (dE > 0.0) { 
      printf("dE = %f, rejecting step", dE);
      e_current = e_old;
      grad_current = grad_old;
      Gmax = grad_current.Max(true);
      Gnorm=grad_current.RMS();
    }*/
      
    //printf("Cycle %d:     dE = %15.9f  |dGrad| = %10.6f  step = %8.3f\n",opt_cycle, dE, dG, stepsize);
    printf("Cycle %d:     dE = %15.12f  Max(Grad) = %10.9f  RMS(Grad) = %10.9f\n",opt_cycle, dE, fabs(Gmax), fabs(Gnorm));
      
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
    opt_cycle++;

    //print fractional coordinates
    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      FILE *fract;
      string fract_file = "fract_coord.frac";
      if ((fract = fopen(fract_file.c_str(),"w"))==NULL) {
        printf("Opt:LBFGS : Cannot open file '%s'\n",fract_file.c_str());
        exit(1);
      } 
      Cluster::cluster().PrintFractionalCoordinates(fract);
      fclose(fract);
    }
  }  

  if(fabs(dE) <= econv  && fabs(Gnorm) <= gconv  && fabs(Gmax) <= Gmax_conv )
    printf("Converged!\n");

}


void Opt::ConjugateGradients() {

  /*Conjugate Gradient formula:
  1st iteration do Steepest Descent
  From there on out do the following:
  1. Calculate steepest direction (/delta(x_n) = -gradient(f(x_n)))
  2. Compute b_n ***see below
  3. Update the conjugate direction s_n = /delta(x_n) + b_n*s_n-1
  4. Perform a line search (optimize /alpha_n = min(f(x_n + /alpha*s_n))
  5. Update the position: x_n+1 = x_n + /alpha_n*s_n

  Fletcher-Reeves:
            /delta(x^T_n)/delta(x_n)
  b_n =  --------------------------------
          /delta(x^T_n-1)/delta(x_n-1) 

  Polak-Ribiere: (used in Quantum Espresso)

            /delta(x^T_n)(/delta(x_n) - /delta(x_n-1))
  b_n =  ------------------------------------------------
                   /delta(x^T_n-1)/delta(x_n-1) 
  */

  //double stepsize_new = 0.75, stepsize_old = 0.75;
  double g_min;
  //printf("Doing Conjugate Gradient\n");
  printf("maxcycle = %i, econv =%f, gconv = %f\n",Params::Parameters().GetMaxOptCycles(),econv,gconv);
  double beta_num, beta_denom;
  reset = false;
  Vector beta;
  Vector all_coords;

  while (opt_cycle <= Params::Parameters().GetMaxOptCycles() && 
	   (fabs(dE) > econv  || fabs(Gnorm) > gconv  || fabs(Gmax) > Gmax_conv) ){

    stepDir_old.PrintGradient("stepDir_old");
    printf("\n");
    coords_old.PrintGradient("coords_old");
    printf("\n");
    grad_old.PrintGradient("grad_old");
    printf("\n");
    
    if (opt_cycle-1 == 1) {
      printf("\nrunning a Steepest Descent step\n");
      //Making our first stepsize = stepDir.Norm()
      stepDir = grad_current;
      stepsize_new = stepsize_old = stepDir.Norm()*stepsize_new;
      stepDirNorm = stepDir.Norm();
      stepDir.Scale(-1.0*stepsize_new/stepDirNorm);

      coords_new = coords_old;
      coords_new += stepDir;

      //pred = 0.5*grad_old.DotProduct(stepDir_old);
      //pred += e_old;
      stepDir.Scale(1.0/stepsize_new);
      //reset = false;
    }
    else if (opt_cycle > 2) {

      beta_denom = grad_old.DotProduct(grad_old);

      //Fletcher-Reeves
      beta_num = grad_current.DotProduct(grad_current);

      //Polak-Ribiere
      //beta = grad_current;
      //beta -= grad_old;
      //beta_num = grad_current.DotProduct(beta);

      beta = stepDir_old;
      beta_num /= beta_denom;
      if(beta_num < 0.0)
        beta_num = 0.0;

      //beta.Scale(beta_num); 
      beta.Scale(beta_num*stepDirNorm_old); 

      printf("beta numerator = %f, beta denom = %f\n",beta_num,beta_denom);

      stepDir = grad_current;
      stepDir.Scale(-1.0);
      stepDir += beta;
      stepDirNorm = stepDir.Norm();
      stepDir.Scale(1.0/stepDirNorm);

      stepDir.PrintGradient("stepDir before check");
      printf("\n");
      coords_new.PrintGradient("coords_new before check");
      printf("\n");
     
      stepsize_new = CheckStep(stepsize_old, stepDir_old);
      printf("stepsize_old = %f, stepsize_new = %f\n",stepsize_old, stepsize_new);

      stepDir.PrintGradient("stepDir after check");
      printf("\n");
      coords_new.PrintGradient("coords_new after check");
      printf("\n");

      stepDir.Scale(stepsize_new);
      coords_new += stepDir;
      stepDir.Scale(1.0/stepsize_new);

    } 

    /*if ((Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())) && opt_cycle > 1){
      printf("found this\n");
      CoupleGradients();
    }
    else if (Params::Parameters().IsPeriodic() && opt_cycle > 1) {
      Cluster::cluster().UpdateLatticeParams(coords_new);
    }*/

 
    tmp_coords = coords_new;

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      for(int i=tmp_coords.GetLength()-3; i < tmp_coords.GetLength(); i++)
        tmp_coords[i] *= RadiansToDegrees;
    }
 
    if (Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
      CoupleGradients();
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
    }
    else if (Params::Parameters().IsPeriodic()) {
      Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,true);
      Cluster::cluster().UpdateLatticeParams(tmp_coords);
      tmp_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(tmp_coords,false,true);
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
      all_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coords,true,true);
    }
    else {
      Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,false);
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,false);
    }

    Cluster::cluster().SetNewCoordinates(all_coords);


    stepDir.PrintGradient("stepDir current");
    printf("\n");
    coords_new.PrintGradient("coords new");
    printf("\n");

    // Update the coordinates in the cluster object
    //Cluster::cluster().SetNewCoordinates(coords_new);

    e_old = e_current;
    grad_old = grad_current; 
    stepDir_old = stepDir; 
    stepDirNorm_old = stepDirNorm;
    stepsize_old = stepsize_new;

    // Create the new jobs, run them, and get the HMBI energy
    Cluster::cluster().RunJobsAndComputeEnergy();
    e_current = Cluster::cluster().GetHMBIEnergy();
    grad_current = Cluster::cluster().GetHMBIGradient();

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      for(int i=grad_current.GetLength()-3; i < grad_current.GetLength(); i++)
        grad_current[i] *= DegreesToRadians;
    }

    printf("e_current = %15.12f, e_old = %15.12f\n", e_current, e_old);
    grad_current.PrintGradient("grad_current");
    printf("\n");


    //update optimization exit params
    dE = e_current - e_old;
    Gmax = grad_current.Max(true);
    //g_min = grad_current.Min(true);
    Gnorm=grad_current.RMS();
    //dGrad -= grad_current;
    //double dG = dGrad.Max(true); 
 
    //grad_current.Scale(-1.0);

    //don't use dE for convergence for first step 
    if(opt_cycle == 1){
      dE = 0.0;
    }    
    /*else if(opt_cycle%30 == 0 ){
      //printf("running a steepest descent step to clear memory issues\n");
      reset = true;
    }*/
    /*else if(opt_cycle == 2 )
      coords_old = coords_new;*/

    printf("Cycle %d:     dE = %15.12f  Max(Grad) = %10.6f  RMS(Grad) = %f\n",opt_cycle, dE, fabs(Gmax), fabs(Gnorm));
      
    Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
     
    if (dE<=0.0000) { 
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

      //print fractional coordinates
      if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
        FILE *fract;
        string fract_file = "fract_coord.frac";
        if ((fract = fopen(fract_file.c_str(),"w"))==NULL) {
          printf("Opt:LBFGS : Cannot open file '%s'\n",fract_file.c_str());
          exit(1);
        } 
        Cluster::cluster().PrintFractionalCoordinates(fract);
        fclose(fract);
      }
    }

    opt_cycle++;
  }

  if(fabs(dE) <= econv  && fabs(Gnorm) <= gconv  && fabs(Gmax) <= Gmax_conv )
    printf("Converged!\n");

}

void Opt::LBFGS() {
  /*
  Quasi-Newton method of optimization

  search direction:
  s_k = -H_k * g_k
  where g_k is the gradient and H_k is the inverse of the hessian
  Note my gradients are negative so our search direction will be 
  s_k = H_k * g_k
  */
  
  //Number of steps we will save
  int m = 4;
  double /*stepsize_new = 0.9, stepsize_old = 0.9,*/ tmp, tmp2;
  Vector all_coords;
  coordinates.Initialize(grad_current.GetLength(),m+1);  
  gradients.Initialize(grad_current.GetLength(),m+1); 
  coord_diff.Initialize(grad_current.GetLength(),m); 
  grad_diff.Initialize(grad_current.GetLength(),m); 
  roe.Initialize(m);

  
  /*coordinates.Print("coords");  
  gradients.Print("gradients"); 
  coord_diff.Print("coord diff"); 
  grad_diff.Print("grad diff"); 
  roe.Print("roe"); 
  exit(0);*/

  while (opt_cycle <= Params::Parameters().GetMaxOptCycles() && 
	   (fabs(dE) > econv  || fabs(Gnorm) > gconv  || fabs(Gmax) > Gmax_conv ) ){//JLM added 1 more check || dE >= 0.0000000

    //Loop to check the previous step and form a new one
    if (opt_cycle > 2){
      stepsize_new = CheckStep(stepsize_new, stepDir_old);

      coords_new.PrintGradient("coords new after check");

      //tmp = stepDir.Norm();
      //tmp2 = min(stepDir.Norm(),stepsize_new);
      //stepDir.Scale(stepsize_new/tmp);
      stepDir.Scale(stepsize_new);
      coords_new += stepDir;
      stepDir.Scale(1.0/stepsize_new);
      //stepDir.Scale(tmp/stepsize_new);
      
    }

    if(opt_cycle == 2){
      //Run a steepest descent step first and then begin LBFGS algorithm
      //SteepestDescent(coords_new, coords_old, grad_current, stepsize_new);
      lbfgs_sd_step = true;
      coordinates.SetColumnVector(coords_old,0);
      gradients.SetColumnVector(grad_current,0);

      //Making our first stepsize = stepDir.Norm()
      stepDir = grad_current;
      stepsize_new = stepsize_old = stepDir.Norm()*stepsize_new;
      stepDirNorm = stepDir.Norm();
      stepDir.Scale(1.0/stepDirNorm);
      stepDir.Scale(-1.0*stepsize_new);

      coords_new = coords_old;
      coords_new += stepDir;

      //pred = 0.5*grad_old.DotProduct(stepDir_old);
      //pred += e_old;
      stepDir.Scale(1.0/stepsize_new);

    }

    /*if ((Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())) && opt_cycle > 1){
      CoupleGradients();
    }
    else if (Params::Parameters().IsPeriodic() && opt_cycle > 1) {
      Cluster::cluster().UpdateLatticeParams(coords_new);
    }*/


    tmp_coords = coords_new;

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      for(int i=tmp_coords.GetLength()-3; i < tmp_coords.GetLength(); i++)
        tmp_coords[i] *= RadiansToDegrees;
    }

    if (Params::Parameters().AreGradientsCoupled() && (Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly())){
      CoupleGradients();
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
    }
    else if (Params::Parameters().IsPeriodic()) {
      Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,true);
      Cluster::cluster().UpdateLatticeParams(tmp_coords);
      tmp_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(tmp_coords,false,true);
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,true);
      all_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(all_coords,true,true);
    }
    else {
      Cluster::cluster().MaintainCartesianSymmetry(tmp_coords,false);
      all_coords= Cluster::cluster().GetSymmetryImposedCoordinates(tmp_coords,false);
    }

    Cluster::cluster().SetNewCoordinates(all_coords);

    stepDir.PrintGradient("stepDir current");
    printf("\n");
    coords_new.PrintGradient("coords new");
    printf("\n");

    // Update the coordinates in the cluster object
    //Cluster::cluster().SetNewCoordinates(coords_new);

    //Update old positions and old gradients
    e_old = e_current;
    grad_old = grad_current;
    //coords_old = coords_new; 
    stepDir_old = stepDir; 
    stepDirNorm_old = stepDirNorm;
    stepsize_old = stepsize_new;

    grad_old.PrintGradient("grad_old");
    printf("\n");

    // Create the new jobs, run them, and get the HMBI energy
    Cluster::cluster().RunJobsAndComputeEnergy();
    e_current = Cluster::cluster().GetHMBIEnergy();
    grad_current = Cluster::cluster().GetHMBIGradient();

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      for(int i=grad_current.GetLength()-3; i < grad_current.GetLength(); i++)
        grad_current[i] *= DegreesToRadians;
    }

    printf("e_current = %f, e_old = %f\n", e_current, e_old);
    grad_current.PrintGradient("grad_current");
    printf("\n");


    //update optimization exit params
    dE = e_current - e_old;
    Gmax = grad_current.Max(true);
    //g_min = grad_current.Min(true);
    Gnorm=grad_current.RMS();
    //dGrad -= grad_current;
    //double dG = dGrad.Max(true); 
 
    //grad_current.Scale(-1.0);

    //don't use dE for convergence for first step 
    if(opt_cycle == 1){
      dE = 0.0;
    }


    printf("Cycle %d:     dE = %15.12f  Max(Grad) = %10.6f  RMS(Grad) = %f\n",opt_cycle, dE, fabs(Gmax), fabs(Gnorm));
      
    Cluster::cluster().UpdateTrajectoryFile(opt_cycle);
     
    if (dE<=0.0000) { 
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

      //print fractional coordinates
      if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
        FILE *fract;
        string fract_file = "fract_coord.frac";
        if ((fract = fopen(fract_file.c_str(),"w"))==NULL) {
          printf("Opt:LBFGS : Cannot open file '%s'\n",fract_file.c_str());
          exit(1);
        } 
        Cluster::cluster().PrintFractionalCoordinates(fract);
        fclose(fract);
      }

    }

    opt_cycle++;
  }

  if(fabs(dE) <= econv  && fabs(Gnorm) <= gconv  && fabs(Gmax) <= Gmax_conv) // && dE <= 0.0
    printf("Converged!\n");
  else
    printf("Not converged\n");
}

void Opt::CoupleGradients() {

  //Correcting coordiate locked under symmetry 
  //Cluster::cluster().MaintainCartesianSymmetry(coords_new,true);
  //coords_new.PrintGradient("Cart symm");
  //Need to convert to fractional coordinates using the old lattice parameters first
  //Vector tmp_coords = coords_new;
  printf("\n");
  tmp_coords.PrintGradient("before coupling");
  printf("\n");
  for (int i = 3*Natoms; i<3*Natoms+6; i++){
    tmp_coords[i] = coords_old[i];
  }
  //convert to fractional coordinates
  tmp_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(tmp_coords,false,true);
  tmp_coords.PrintGradient("old params, fract");
  printf("\n");

  //update new lattice parameters
  Cluster::cluster().UpdateLatticeParams(coords_new);
  printf("\n");
  for (int i = 3*Natoms; i<3*Natoms+6; i++){
    tmp_coords[i] = coords_new[i];
  }
  tmp_coords.PrintGradient("new params, fract");
  printf("\n");
  //convert back to cartesian
  tmp_coords = Cluster::cluster().ConvertBetweenFractionAndCartesianCoordinates(tmp_coords,true,true);
  tmp_coords.PrintGradient("after coupling");
  printf("\n");
  //coords_new.PrintGradient("Cart coupled grad");
  //getting full coordinates from the symmetry imposed coordinates
  //coords_new = Cluster::cluster().GetSymmetryImposedCoordinates(coords_new,true);
  //coords_new.PrintGradient("Symm imposed");

}

double Opt::CheckStep( double alpha_in, Vector& step_dir ) {

//General check of the previous step taken

//Ensures the step length (/alpha) never becomes too large
//Done adjusting alpha if the following is satisfied:
//f(x + /alpha_i*step_dir) - f(x)  <= /alpha_i*t
//t = c*step_dir*old_grad and 0</tau,c<1
//Until satisfied set /alpha_i = /tau*/alpha_i-1

  //armijoConditionMet = false;
  double alpha_out;
  double c1 = 0.0001, c2 = 0.10, tau = 0.5, alpha_keep = 0.75;
  double t=0.0, proj, oldproj, svar, alpha, tmp;
  bool stepAccepted=false;
  wolfePassed = false;
  armijoConditionMet=false;

  step_dir.Print("Step dir 1");

  grad_old.Print("grad_old");

  if (dE > 0.00000){
    printf("Rejecting step, new energy %6.12f larger\n", dE);
  }
  else{
    //Now check the Armijo-Golstein condition
    t = c1*step_dir.DotProduct(grad_old); 
    if(abs(t) < pow(10,-9))
      t = 0.0;
    if(dE <= (t*alpha_in)){
     printf("Armijo condition met!\n");
     armijoConditionMet=true;
     stepAccepted = true;
     alpha_out = alpha_in;

     //Next check the Wolfe condition
     t = step_dir.DotProduct(grad_current);
     //printf("step_dir.DotProduct(grad_current) = %6.12f\n",t);
     t = c2*step_dir.DotProduct(grad_old);
     //printf("c2*step_dir.DotProduct(grad_old) = %6.12f\n",t);
     
     if (type_opt == "LBFGS") {
       c2 = 0.90;
       //Strong Wolfe condition
       // |step_dir*grad_current| <= c2*|step_dir*grad_old|
       if( abs(step_dir.DotProduct(grad_current)) <= c2*abs(step_dir.DotProduct(grad_old)) ){
         printf("Strong Wolfe condition met! Accepting step\n");
         wolfePassed = true;
       }
       else{
         printf("Strong Wolfe condition failed. Increasing step size for the next iteration\n");
         //if(alpha_keep == alpha_in)
          // alpha_keep *= 1.1;

         //alpha_out = alpha_in+alpha_in*0.1;
         //stepAccepted = false;
       }
     }
     else{
       //Wolfe condition
       // step_dir*grad_current <= c2*step_dir*grad_old
       if( (step_dir.DotProduct(grad_current)) <= (c2*step_dir.DotProduct(grad_old)) ){
         printf("Wolfe condition met! Accepting step\n");
         wolfePassed = true;
       }
       else{
         printf("Wolfe condition failed. Increasing step size for the next iteration\n");
         //alpha_out = alpha_in+alpha_in*0.1;
       }
     }
    }
    else{   
     printf("Armijo condition not met, Rejecting Step\n");
     //reset = true;
     //alpha_out = alpha_in*tau;
    }
  }

  if(stepAccepted){
    numStepAccepted++;
    coords_old = coords_new;
    if(lbfgs_sd_step)
      lbfgs_sd_step = false;
    if(type_opt == "LBFGS"){
      //alpha_out = alpha_keep;
      UpdateHessian();
      /*if(stepDirNorm == stepDirNorm_old)
        alpha_out = alpha_in;
      else
        alpha_out = TrustRadius(alpha_in);*/
    }
    //else{
      alpha_out = TrustRadius(alpha_in);
    //}
    //else {
      //alpha_out = alpha_in;

      //if (!wolfePassed)
      //  alpha_out *= 1.1;

      //if(numStepAccepted%5 == 0){
      //  alpha_out*=1.25;
      //}
    //}

  }
  else {

    grad_current = grad_old;
    //if(type_opt == "LBFGS"){
      //alpha_out = stepDir_old.Norm();
      //alpha_out /= 4.0;
      //alpha_out = TrustRadius(alpha_in);
    //}
    //else{
      //alpha_out = alpha_in*tau;
      
    //}

    coords_new = coords_old;
    e_current = e_old;
    stepDir = stepDir_old;
    stepDirNorm = stepDirNorm_old;
    numStepAccepted = 0;
    alpha_out = TrustRadius(alpha_in);
  }

  printf("alpha in = %f, alpha out = %f, stepDir.Norm = %f\n", alpha_in, alpha_out, stepDirNorm);
  printf("t = %15.12f, e_old = %f, e_current = %f\n", t, e_old, e_current);
  
  return alpha_out;
  
}

double Opt::TrustRadius(double tr) {
  /*Goal is to create radius that the program may step across next
  Built based off paper from Ayers (2014) */
  double test, r1 = 0.1, r2 = 0.75, c1 = 0.5, c2 = 1.5, tmp;
  double dEOs, den, scalar = 1.0;
  //printf("Began trust radius calculation!\n");
  
  //tmp = stepDir_old.Norm();
  if(dE > 0.0){
    numStepsRejected++;
    if(numStepsRejected == 3 && type_opt != "SD" && !lbfgs_sd_step){
      printf("Running a Steepest Descent Step instead\n");
      lbfgs_sd_step = true;
      stepDir = grad_current;
      stepDirNorm = stepDir.Norm();
      stepDir.Scale(-1.0/stepDirNorm);
      printf("Reseting trust radius\n");
      tr = 0.65*stepDirNorm;
      numStepsRejected = 0;
    }
    else{
      //tr = min(stepDir_old.Norm(),tr)*0.5;
      //stepDir_old.Scale(0.5);
      dEOs = tr*grad_old.DotProduct(stepDir_old);//(tr/tmp)
      den = dE - dEOs;
      tr = (-0.5*dEOs*tr)/den;

      if(tr < tr_min){
        //tr = tr_min;
        //tr = max(tr_min,stepDir_old.Norm());
        tr = max(tr_min,stepDirNorm_old);
        printf("Warning! trust radius became lower than the tr_min\n");
      }
    }
  }
  else{
    numStepsRejected = 0;
    //note: wouldn't be here is armijo condition wasn't met 
    if( armijoConditionMet && (stepDir_old.Norm() > tr)) scalar = 1.5;
    else scalar = 1.1;    

    if(wolfePassed) scalar*=1.75;

    tr = min(tr_max,scalar*tr);
    tr = min(tr,stepDirNorm);
    //tr = min(tr,stepDir_old.Norm());

    if(tr < tr_min){
      //tr = tr_min;
      //tr = max(tr_min,stepDir.Norm());
      tr = max(tr_min,stepDirNorm);
      printf("Warning! trust radius became lower than the tr_min\n");
    }
    printf("tr = %f scalar = %f stepDir.Norm() = %f\n", tr, scalar, stepDirNorm);

  }


  
  /*if(dE > 0.0){
    tr = min(stepDir_old.Norm(),tr)*0.5;
    stepDir_old.Scale(0.5);
  }
  else{
    //stepDir_old.Scale(tr/tmp);
    //pred = 0.5*grad_old.DotProduct(stepDir_old);
    //pred = grad_old.DotProduct(stepDir_old)+0.5*stepDir_old.DotProduct(stepDir_old);
    pred = (tr-0.5*pow(tr,2))*grad_old.DotProduct(stepDir_old);
    //pred += e_old;
    //stepDir_old.Scale(tmp/tr);
    //test = -1.0*dE/(e_old - pred);
    test = -1.0*dE/(pred);
    printf("pred = %f e_old = %f dE = %f pred - e_old =%f test = %f\n", pred, e_old, dE,e_old - pred, test);

    if(pred > 0.9)
      tr = min(c2*tr, tr_max);
    else if(pred <0.1)
      tr = c1*tmp;
      


    //printf("test = %f \n",test);
    //Energy based trust radius **not as strict as the gradient based trust radius
    //if( test < r1){
    //  tr = c1*tmp;
    //}
    //else if(test > r2){
    //  tr = min(c2*tr, tr_max);
   // }
    //else{
   //   tr = stepDirNorm;
    //}
  }*/



  //printf("new tr = %f stepDir.Norm() = %f tr/tmp = %f\n", tr, stepDir.Norm(), tr/stepDir.Norm());
  return tr;


}

void Opt::UpdateHessian() {
  
  Vector tmp, q, z, alpha, beta, tmp2;
  double tmpNum, hess, tmpNum2;
  alpha.Initialize(coord_diff.GetCols());
  beta.Initialize(coord_diff.GetCols());

  /*coordinates.Print("coords");  
  gradients.Print("gradients"); 
  coord_diff.Print("coord diff"); 
  grad_diff.Print("grad diff"); 
  roe.Print("roe");   */

  //Calculate new values
  /*tmp = coords_old;
  tmp -= coordinates.GetColumnVector(1);

  tmp2 = grad_current;
  tmp2 -= gradients.GetColumnVector(1);    

  tmpNum = tmp.DotProduct(tmp2);

 if (tmpNum < 0.0){
   printf("Running a SteepestDescent Step instead\n");
   stepDir = grad_current;
   stepDir.Scale(-1.0);
 }
 else{*/

  //First move everything back one column
  for (int i=coordinates.GetCols()-1; i>0; i-- ) {
    coordinates.SetColumnVector(coordinates.GetColumnVector(i-1), i);
    gradients.SetColumnVector(gradients.GetColumnVector(i-1), i);
  }
  coordinates.SetColumnVector(coords_old,0);
  gradients.SetColumnVector(grad_current,0);

  for (int i=coord_diff.GetCols()-1; i>0; i-- ){
    coord_diff.SetColumnVector(coord_diff.GetColumnVector(i-1), i);
    grad_diff.SetColumnVector(grad_diff.GetColumnVector(i-1), i);
    roe[i] = roe.Element(i-1);
  }
    
  //Calculate new values
  tmp = coords_old;
  //tmp.Print("coord_diff col vec 1");
  //printf("tmp length = %d, coords_old = %d,\n", tmp.GetLength(), coords_old.GetLength()); 
  tmp -= coordinates.GetColumnVector(1);
  coord_diff.SetColumnVector(tmp,0);
printf("\n");
  tmp.PrintGradient("coord_diff");
printf("\n");
  tmp = grad_current;
  tmp -= gradients.GetColumnVector(1);    
  grad_diff.SetColumnVector(tmp,0);

  tmp.PrintGradient("grad_diff");
printf("\n");
  tmpNum = coord_diff.GetColumnVector(0).DotProduct(grad_diff.GetColumnVector(0));
  //roe[0] = fabs(1.0 / tmpNum);
  roe[0] = 1.0 / tmpNum;

  /*coordinates.Print("coords");  
  gradients.Print("gradients"); 
  coord_diff.Print("coord diff"); 
  grad_diff.Print("grad diff"); */
  printf("\n");
  roe.Print("roe"); 
  printf("\n");

  //Now begin calculating everything
  q = grad_current;

  //for(int i = coord_diff.GetCols()-1; i>0; i--) {
  for(int i = coord_diff.GetCols()-1; i>=0; i--) {
    alpha[i] = (coord_diff.GetColumnVector(i)).DotProduct(q);
    alpha[i] *= roe[i];
    //printf("i = %d, alpha = %f\n", i,alpha[i]);    
    tmp = grad_diff.GetColumnVector(i);
    tmp.Scale(alpha[i]);
    q -= tmp;
    //q.Print("q");
    //printf("\n");
  }
  
  /*tmpNum = grad_diff.GetColumnVector(1).DotProduct(grad_diff.GetColumnVector(1));
  //If this is the first time building the "hessian" then make it an identity
  if (tmpNum == 0.0000)
    hess = 1.0;
  else
    hess = (coord_diff.GetColumnVector(1)).DotProduct(grad_diff.GetColumnVector(1))/tmpNum;*/

  tmpNum = grad_diff.GetColumnVector(0).DotProduct(grad_diff.GetColumnVector(0));
  hess = (coord_diff.GetColumnVector(0)).DotProduct(grad_diff.GetColumnVector(0))/tmpNum;

  //printf("\nhess = %15.12f\n",hess);

  z = q;
  z.Scale(hess);
  //z.Print("z");

  //for(int i = 1; i<=coord_diff.GetCols()-1; i++) {
  for(int i = 0; i<=coord_diff.GetCols()-1; i++) {
    beta[i] = (grad_diff.GetColumnVector(i)).DotProduct(z);
    beta[i] *= roe[i];

    tmp = coord_diff.GetColumnVector(i);
    tmp.Scale(alpha[i] - beta[i]);
    z += tmp;

    //printf("\n");
    //printf("i = %d, alpha = %f, beta = %f\n", i,alpha[i],beta[i]);
    //z.Print("z");
  }

  z.Scale(-1.0);
  stepDir = z;
  //if(stepDir.DotProduct(grad_old) > 0.0);
  //  stepDir.Scale(-1.0); 
 //}
  stepDirNorm = stepDir.Norm();
  stepDir.Scale(1.0/stepDirNorm);

}

void Opt::SteepestDescent(Vector& Rout, Vector& Rin, Vector& Grad, 
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



void Opt::WriteHessianToFile(string& filename, int& Natoms, Vector& current_coords, Matrix& Hessian) {

  string datafile = Params::Parameters().GetTmpFilesPath() + "/" + filename;
  /* Create an empty file */
  FILE *data;

  if ((data = fopen(datafile.c_str(),"w"))==NULL) {
    printf("HMBI : Cannot open file '%s'\n",
        datafile.c_str());
    exit(1);
  }

  fprintf(data,"Current HMBI Coordiantes\n");
  for (int i=0;i<Natoms;i++) {
    fprintf(data,"%2s %15.10f %15.10f %15.10f\n",Cluster::cluster().GetAtomicSymbol(i).c_str(),
           current_coords[3*i],current_coords[3*i+1],current_coords[3*i+2]);
  }
  if ( (Params::Parameters().IsPeriodic() ) && (Params::Parameters().GetMMType() != 2) ) {
    fprintf(data,"UnitCellLength %15.10f %15.10f %15.10f\n",
           current_coords[3*Natoms],current_coords[3*Natoms+1],current_coords[3*Natoms+2]);
    fprintf(data,"UnitCellAngles %15.10f %15.10f %15.10f\n",
           current_coords[3*Natoms+3],current_coords[3*Natoms+4],current_coords[3*Natoms+5]);
  }

  fprintf(data,"HMBI analytical hessian at lower level of theory for geometry optimization:\n");
  int cols = Hessian.GetCols();
  int rows = Hessian.GetRows();

  for (int i=0;i<(cols-cols%9)/9;i++) {    
    for (int j1=0;j1<1;j1++) {
      fprintf(data,"  \t");
      for (int k1=9*i;k1<9*i+9;k1++) {
        fprintf(data,"%15.6d \t", k1);
      }
      fprintf(data,"\n");
    }
    for (int j=0;j<rows;j++) {
      fprintf(data,"%d \t", j);
      for (int k=9*i;k<9*i+9;k++) {
        fprintf(data,"%15.6f \t", Hessian.Element(j,k) );
      }
      fprintf(data,"\n");
    }
  }
  if (cols%9 != 0) {
    for (int j1=0;j1<1;j1++) {
      fprintf(data,"  \t");
      for (int k1=(cols-cols%9);k1<cols;k1++) {
        fprintf(data,"%15.6d \t", k1);             
      }
      fprintf(data,"\n"); 
    }       
    for (int j=0;j<rows;j++) {
      fprintf(data,"%d \t", j);
      for (int k=(cols-cols%9);k<cols;k++) { 
        fprintf(data,"%15.6f \t", Hessian.Element(j,k) );
      }
      fprintf(data,"\n");
    }
  }

  fprintf(data,"End of File\n");
  fclose(data);

}



void Opt::ReadHessianFromFile(string& filename, int& Natoms, Matrix& Hessian) { //reads any matrix from a file

  string datafile = Params::Parameters().GetTmpFilesPath() + "/" + filename;

  // Open the data file            
  ifstream infile;
  infile.open( datafile.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadHessian : Cannot open file '%s'\n",        
           datafile.c_str());                  
    exit(1);
  }

  string line;
  while ( !infile.eof() ) {
    getline(infile,line);              
    if ( line.substr(0,77)=="HMBI analytical hessian at lower level of theory for geometry optimization:\n" ) {
      for (int k=0;k<(3*Natoms-3*Natoms%9)/9;k++) {
        getline(infile,line); // throw away header line

        for (int i=0;i<3*Natoms;i++) {
          getline(infile,line);
          istringstream iss(line);
          string tmp;
          iss >> tmp; // throw away the atom index
          for (int j=9*k;j<9*k+9;j++) {
            iss >> Hessian.Element(i,j); // Store the hessian elements
          }
        }
      }
      getline(infile,line); // check if this line signals end of hessian print or not
      if (line.substr(0,11) != "End of File") {
        for (int i=0;i<3*Natoms;i++) {
          getline(infile,line);  
          istringstream iss(line);
          string tmp;
          iss >> tmp; // throw away the atom index
          for (int j=(3*Natoms-3*Natoms%9);j<3*Natoms;j++) {
            iss >> Hessian.Element(i,j); // Store the hessian elements
          }
        }                  
        break;
      }

      else {
        break;
      }
    }
  }
}



