#ifndef _opt_h
#define _opt_h
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "params.h"
#include "vector.h"
#include "matrix.h"
#include "cluster.h"
#include "constants.h"
using namespace hmbi_constants;


class Opt {

 private:
  int opt_cycle, Natoms, numStepAccepted, numStepsRejected;
  double dE, Gnorm, Gmax, StepNorm, StepMax, stepsize_new, stepsize_old; 
  double reduced_Gnorm, reduced_Gmax, constraint_norm;
  double econv, gconv, Gmax_conv, StepMax_conv, StepNorm_conv;
  double stepsize, reduced_Gconv, reduced_Gmax_conv, constraint_conv;
  double e_current, e_old, stepDirNorm_old, stepDirNorm;
  double pred, tr_max, tr_min, tr_init;
  bool reset, armijoConditionMet, wolfePassed, lbfgs_sd_step;
  string type_opt;
  Vector stepDir;
  Vector grad_current;
  Vector stepDir_old;
  Vector grad_old;
  Vector coords_new;
  Vector coords_old;
  Vector tmp_coords;

  Matrix coordinates;  
  Matrix gradients; 
  Matrix coord_diff; 
  Matrix grad_diff; 
  Vector roe;
  //bool Optimize();
 public:

  Opt();  // The real work is done by ReadInputFile
  ~Opt(); 
  void OptimizeGeometry(string type);

  //void SteepestDescent(const Vector& Rout, const Vector& Rin, const Vector& Gradient);
  void SteepestDescent();
  void SteepestDescent(Vector& Rout, Vector& Rin, Vector& Gradient, double stepsize);
  void ConjugateGradients();
  void LBFGS();

  void CoupleGradients();
  double TrustRadius(double tr);
  double CheckStep(double alpha_in, Vector& step_dir);
  void UpdateHessian();

  /*void ConjugateGradients(Vector& Rout, Vector& StepDir_out, Vector& Rin, 
			Vector& Grad, Vector& Grad_old, Vector& StepDir_old,
			double stepsize, bool reset);

  void ConstrainedLBFGSOpt(Vector& Rout, Vector& Rin,
                        Vector& Grad, Vector& Grad_old,
                        Vector& delta_q, bool reset,
                        Matrix& ykmat, Matrix& skmat, int M, int opt_cycle,
                        Vector& delta_x, Vector& diff_grad, Vector& rho, bool SD_reset,
                        Vector& delta_y, Matrix& T_mat_b, Matrix& T_mat_ti,
                        Vector& reduced_grad_ti, Vector& reduced_grad_ti_last,
                        double* trust_radius, double* svar);

  Vector TwoLoopRecursion(Vector& reduced_grad_ti_last, Matrix& Init_Inv_Hess,
                        Matrix& skmat, Matrix& ykmat, Vector& rho,
                        int M, int opt_cycle );

  Matrix ConstructTMatrix();

  int StepAcceptance(double E_new, double E_old, Vector& stepsize, Vector& grad_new, Vector& grad_old);

  int StepAcceptance_new(double E_new, double E_old, Vector& stepsize, Vector& grad_new, Vector& grad_old,
                       double merit_fn_new, double merit_fn_old, Vector& merit_fn_grad, Vector& merit_fn_grad_old,
                       Matrix& T_mat_b, Vector& delta_y, double* Mu);

  void lbfgs_store(Matrix& skmat, Matrix& ykmat, Vector& dR, Vector& Diff_Grad, Vector& rho, int opt_cycle, int M);

  Vector GetOldGradVector(Matrix& ykmat, Matrix& skmat, int opt_cycle, int M);

  void lbfgs_matrix_reshuffle(bool reshuffle, Matrix& skmat, Matrix& ykmat, Vector& rho, int M, int N);

  void DoNothing(int opt_cycle);

  Vector Get_dR_Given_delta_x(Vector& delta_x, Vector& delta_y, Matrix& T_mat_b, Matrix& T_mat_ti);

  Vector Get_delta_x_Given_dR(Vector& dR, Vector& delta_y, Matrix& T_mat_b, Matrix& T_mat_ti);
*/
  void WriteHessianToFile(string& filename, int& Natoms, Vector& current_coords, Matrix& Hessian);

  void ReadHessianFromFile(string& filename, int& Natoms, Matrix& Hessian);

  /*Matrix GetTmpMMNucHessian();

  Matrix GetTmpMMLatHessian();

  double LatBlockApproxHessian_FDGrad(int k, int i, double deltak, double deltai, string pm);*/

};

#endif
