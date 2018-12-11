#ifndef _main_h
#define _main_h

#include <cstring>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
using std::string;
#include "atom.h"
#include "monomer.h"
#include "dimer.h"
#include "params.h"
#include "cluster.h"
#include "vector.h"
#include "opt.h"
#include "nmr.h"
#include "supercell.h"
#include <algorithm>
#include "quasiharmonic.h"
#include "dlf_interface.h"
#include "knitro.h"
#include "knitro_interface.h"
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
#define QCHEMCP_TAG 5
#define MOLPRO_TAG 6
#define CRYSTAL_TAG 7
#define G09_TAG 8
#define GDMA_TAG 9
#define DIETAG 10
#define DALTON_TAG 11
#define QUANTUM_ESPRESSO_TAG 12
#define ORCA_TAG 13
#define PSI4_TAG 14
#define DFTB_TAG 15
#define BUFFSIZE 10000

//Opt *GeomOpt;

// Other function prototypes
void DebugVectorClass();
bool CheckIfCPDimerJob(string job, int job_type);
int CheckIfJobSuccessful(string job, int job_type,bool IsCPDimer);

Vector GetFiniteDifferenceGradient(Vector coords);
Matrix GetFiniteDifferenceHessian(Vector coords, int imon=0);

#endif
