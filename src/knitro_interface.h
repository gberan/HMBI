#ifdef KNITRO
#ifndef _knitro_interface_h
#define _knitro_interface_h

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
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
#include "knitro.h"

#include "constants.h"
using namespace hmbi_constants;

// Function prototypes for K-NITRO to operate in callback mode.

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
					 void * userParams);

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
				 void * userParams);

				

#endif
#endif
