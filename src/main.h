#ifndef _main_h
#define _main_h

#include <iostream>
#include <stdio.h>
#include <fstream.h>
#include <cassert>
#include <iomanip.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <unistd.h> //for getcwd
// following, plus iostream, were all I had before
#include <string>
using std::string;
#include "monomer.h"
#include "dimer.h"
#include "atom.h"
#include "cluster.h"
#include "vector.h"
#include "opt.h"
#include "dlf_interface.h"
#include <time.h>

#include "constants.h"
using namespace hmbi_constants;

// Some extra includes/defines/functions for parallel version
#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */

// Define some tags for job types - mostly for MPI
#define QCHEM_TAG 1
#define TINKER_TAG 2
#define CAMCASP_TAG 3
#define AIFF_DISP_TAG 4
#define DIETAG 10
#define BUFFSIZE 10000

// Other function prototypes
void OptimizeGeometry(string type);
void DebugVectorClass();
int CheckIfJobSuccessful(string job, int job_type);

Vector GetFiniteDifferenceGradient(Vector coords);
Matrix GetFiniteDifferenceHessian(Vector coords, int imon=0);

#endif
