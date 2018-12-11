#ifndef _solvate_h
#define _solvate_h
#include <iostream>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
using std::string;
#include "vector.h" // my vector class, "Vector"
#include "monomer.h"

void AddMolecules(int Nsolvent, int Natoms, string *AtSym, Vector *XYZ );
void AddMolecules(Monomer &Unsolvated, int Nsolvent);
double CheckDistances(Vector *New_Solvent, int Natoms_solvent, Vector *XYZ, int Natoms);
Monomer MergeTwoMonomers(Monomer MonA, Monomer MonB);

double frand();

// File I/O
Monomer ReadCoordinates(ifstream &infile, int Natoms);
void WriteXYZfile(string filename, Monomer m, string title);
void WriteTinkerXYZfile(string filename, Monomer m, string title);




// Create reference solvent molecules
Monomer CreateAmoebaWaterMolecule();
Monomer CreateAmoebaBenzeneMolecule();



#endif

