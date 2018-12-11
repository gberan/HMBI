#ifndef _input_h
#define _input_h
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

// Shouldn't need to define this here again!?
/* void Rewind(ifstream& infile) { */
/*   infile.clear(); */
/*   infile.seekg(0,ios::beg); */
/* } */

// Read sections for different QM packages:
string ReadG09Section(ifstream& infile);
string ReadDaltonSection(ifstream& infile);
// adding string for psi4
string ReadPSI4Section(ifstream& infile);

// Read sections for custom basis sets:
string ReadHBasis( ifstream& infile);
string ReadCBasis( ifstream& infile);
string ReadNBasis( ifstream& infile);
string ReadOBasis( ifstream& infile);
string ReadSBasis( ifstream& infile);
string ReadClBasis( ifstream& infile);
string ReadIBasis( ifstream& infile);
string ReadSnBasis( ifstream& infile);
string ReadPBasis( ifstream& infile);
string ReadKBasis( ifstream& infile);
string ReadNaBasis( ifstream& infile);
string ReadBrBasis( ifstream& infile);
string ReadFBasis( ifstream& infile);



// Read sections for different QM packages:
/* string Cluster::ReadG09Section(ifstream& infile); */
/* string Cluster::ReadDaltonSection(ifstream& infile); */


/* // Read sections for custom basis sets: */
/* string Cluster::ReadHBasis( ifstream& infile); */
/* string Cluster::ReadCBasis( ifstream& infile); */
/* string Cluster::ReadNBasis( ifstream& infile); */
/* string Cluster::ReadOBasis( ifstream& infile); */
/* string Cluster::ReadSBasis( ifstream& infile); */
/* string Cluster::ReadClBasis( ifstream& infile); */
/* string Cluster::ReadIBasis( ifstream& infile); */
/* string Cluster::ReadSnBasis( ifstream& infile); */
/* string Cluster::ReadPBasis( ifstream& infile); */
/* string Cluster::ReadKBasis( ifstream& infile); */
/* string Cluster::ReadNaBasis( ifstream& infile); */
/* string Cluster::ReadBrBasis( ifstream& infile); */
/* string Cluster::ReadFBasis( ifstream& infile); */










#endif
