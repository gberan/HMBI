#ifndef _nmr_h
#define _nmr_h
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <iostream>
using std::string;
#include "params.h"
#include "ee.h"
#include "vector.h"
#include "matrix.h"
#include "cluster.h"
#include "constants.h"
#include "main.h"
using namespace hmbi_constants;

// Define some tags for parallel job types
#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */

//void PredictMagneticProperties();  //UNFNISHED - still in cluster.C

void AssignBasisSetsToAtoms(int NMon, Monomer Monomers[]);
void AssignBasisSetsToAtoms(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] );
//void AssignBasisSetsToAtomsTwoBodyCharge(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] );

void CreateClusterJobs(int NMon, Monomer Monomers[] );
void CreateClusterJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] );
void CreateClusterJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], Matrix EwaldCharges );

void CreateClusterG09Job(int NMon, Monomer Monomers[], int total_charge, int iclust );
void CreateClusterG09Job(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, int iclust );
void CreateClusterG09Job(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, Matrix EwaldCharges, int iclust );

void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int total_charge, int iclust ); 
void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, int iclust ); 
void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int total_charge, Matrix EwaldCharges, int iclust );

//void CreateClusterDALTONJob(int NMon, Monomer Monomers[], int NMon_images, Monomer ImageMonomers[], int total_charge, Matrix EwaldCharges, int iclust ); 

void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ); 
void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ); 
void CreateMonomerAndDimerJobs(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[], Matrix EwaldCharges );

void AssignDimerListsForFragmentCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void AssignDimerListsForFragmentCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] );
//void AssignDimerListsForFragmentTwoBodyChargeCalculation(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] );

void RunFragmentJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[]);
void RunFragmentJobs(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] );

string* GetFragmentJobList(int Njobs, int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] ); 
string* GetFragmentJobList(int Njobs, int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[] ); 

string RunG09ClusterJob(int cluster_number);

string RunDaltonClusterJob(int cluster_number); 

void UpdateJobStatus(int ijob, int NMon_jobs, int NDim_jobs);

void ReadEFGData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void ReadEFGData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]);

void ComputeEFGTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void ComputeEFGTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]);

void ReadNMRData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void ReadNMRData(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]);

void ComputeNMRShieldingTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[] );
void ComputeNMRShieldingTensors(int NMon, Monomer Monomers[], int NDim, Dimer Dimers[], int NDim_images, Dimer DimerImages[]);

/* // Electrostatic Embedding Functions: */

/* void RunElectrostaticEmbeddingJobs(); */


/* // Ewlad - NMR Functions: */

/* //Vector ComputeChargesFromEwaldPotential(Matrix TestPoints, int NTest, Monomer Monomers[], int NMon, double cell_volume, double a, double b, double c, Vector UnitCellx, Vector UnitCelly, Vector UnitCellz, Vector RecipCellx, Vector RecipCelly, Vector RecipCellz ); */
/* Vector ComputeEwaldPotentialAtTestPoints( Matrix TestPoints, int NTest, Monomer Monomers[], int NMon, double cell_volume, double a, double b, double c, Vector *unit_cell, Vector *reciprocal_cell ); */

/* 					  //Vector UnitCellx, Vector UnitCelly, Vector UnitCellz, Vector RecipCellx, Vector RecipCelly, Vector RecipCellz ); */


/* Matrix ComputeEwaldCharges( Matrix TestPoints, int NTest, Vector EwaldPotential, Vector *unit_cell, Monomer Monomers[], int NMon ); */

#endif 
