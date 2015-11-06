#ifndef _opt_h
#define _opt_h
#include <stdio.h>
#include <fstream.h>
#include <cassert>
#include <iomanip.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "params.h"
#include "vector.h"

bool Optimize();
//void SteepestDescent(const Vector& Rout, const Vector& Rin, const Vector& Gradient);
void SteepestDescent(Vector& Rout, Vector& Rin, Vector& Gradient);
void SteepestDescent(Vector& Rout, Vector& Rin, Vector& Gradient, double stepsize);


void ConjugateGradients(Vector& Rout, Vector& StepDir_out, Vector& Rin, 
			Vector& Grad, Vector& Grad_old, Vector& StepDir_old,
			double stepsize, bool reset);

#endif
