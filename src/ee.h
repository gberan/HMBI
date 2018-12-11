#ifndef _ee_h
#define _ee_h
#include <stdio.h>
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
#include "nmr.h"
using namespace hmbi_constants;

#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */
// Define some tags for parallel job types
#define QCHEM_TAG 1
#define TINKER_TAG 2
#define CAMCASP_TAG 3
#define QCHEMCP_TAG 5
#define MOLPRO_TAG 6
#define CRYSTAL_TAG 7
#define G09_TAG 8
#define GDMA_TAG 9
#define DALTON_TAG 11
#define DIETAG 10
#define ORCA_TAG 13


// Electrostatic Embedding Functions:
void DetermineTwoBodyElectrostaticEmbeddingEnvironment(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void DetermineTwoBodyElectrostaticEmbeddingEnvironment(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NMon_images, Monomer MonomerImages[], int NDim_images, Dimer DimerImages[]  );

//void ComputeTwoBodyCharges(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
//void ComputeTwoBodyCharges(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]);



void CreateDaltonPotFile(int k, int type);

void RunElectrostaticEmbeddingJobs();

//void DetermineSelfConsistentElectrostaticEmbeddingEnvironment(); (still in cluster.C)

void ComputeSelfConsistentEmbeddingCharges();

void AssignMonomerListsForElectrostaticEmbedding(int NMon, Monomer Monomers[]);
void AssignMonomerListsForElectrostaticEmbedding(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]);
void AssignMonomerListsForTwoBodyElectrostaticEmbedding(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]);

// EWALD
Matrix ComputeEwaldCharges(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Vector *unit_cell); // DONE

Matrix BuildEwaldTestPointGrid(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[]); // DONE

void ComputeEwaldPotentialAtAtomCenters( int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Vector *unit_cell); // DONE


Matrix GetDistanceOrderedAtomCenteredGrid(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int Ntest, Vector *unit_cell); // DONE

Vector CalculateFixedPotentialAtAtomCenter(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Matrix TestPoints);// DONE

Vector CalculateOptimiumCharges(Matrix TestPoints, Matrix SphereGrid, Vector FitPotential );// DONE


struct MyComparator
{
    const vector<double> & value_vector;

    MyComparator(const vector<double> & val_vec):
        value_vector(val_vec) {}

    bool operator()(double i1, double i2)
    {
        return value_vector[i1] < value_vector[i2];
    }
};

void print(const vector<int> & v, const char * msg);

void SortAaccordingtoB( Vector& A, Vector& B);



/* // Multipole embedding code: */
/* void CreateDaltonPotFile(int k, int type); */


/* } */


#endif
