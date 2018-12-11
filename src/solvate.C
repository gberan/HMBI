#include "solvate.h"
#include "constants.h"
using namespace hmbi_constants;
#include <time.h>




main(int argc, char ** argv) {

  if (argc < 3) {
    printf("Solvate syntax:: solvate <input geom> <# of solvent molecules to add>\n");
    exit(1);
  }

  int N = atoi(argv[2]); // convert char to integer

  printf("Adding %d solvent molecules\n",N);

  // Initialize random seed for random number generator
  srand( time(NULL) ); 
  //srand( 200 ); 

  // Initialize Parameters
  //Params Parameters;
  Params::Parameters().SetParameter("MM_CODE","TINKER");

  // Read XYZ file containing coordinates of current system
  // Open input file
  const int maxlength = 100;
  char filename[maxlength];
  strcpy(filename, argv[1]);
  ifstream infile;
  infile.open(filename);
  assert(infile.is_open()); 
  
  // Read Natoms and comment line
  string line;
  int Natoms;
  getline(infile,line);
  {
    istringstream iss(line);
    iss >> Natoms;
  }
  string title;
  //getline(infile,title); // grab comment line
  printf("Natoms = %d\n",Natoms);

  // Read the input geometry
  Monomer Unsolvated = ReadCoordinates(infile,Natoms);
  Unsolvated.PrintTinkerCartesian();


  // Add N solvent molecules
  AddMolecules(Unsolvated, N);
  Unsolvated.PrintTinkerCartesian();

  // Write new geometry to disk
  string newgeom = "newgeom.xyz";
  char label[10];
  sprintf(label,"%d",N);
  title = "Solvated system";
  //title += " + "; title += label;
  //title += " solvent molecules";
  WriteTinkerXYZfile(newgeom,Unsolvated,title);



  infile.close();
  //delete [] AtSym;
  //delete [] XYZ;

}

void AddMolecules(Monomer &Unsolvated, int Nsolvent) {

  // Create reference solvent molecule
  Monomer Solvent = CreateAmoebaWaterMolecule();
  //Monomer Solvent = CreateAmoebaBenzeneMolecule();
  
  int n_added = 0;
  while (n_added < Nsolvent) {

    // Choose atom from existing geometry;
    int Natoms = Unsolvated.GetNumberOfAtoms();
    int starting_atom = (int) floor(Natoms*frand());
    //printf("Starting atom = %d\n",starting_atom);
    Vector Origin(3);
    Origin = Unsolvated.GetAtom(starting_atom).GetPosition();
    //Origin.Print("Origin");
    
    // Choose distance from starting atom.
    // Base distance = 2A, and use random noise to pick +/- 10%
    double distance = frand()*0.10; // +/- 10% noise
    if (frand() < 0.5) distance *= -1.0; // choose sign of noise
    distance = 2*(distance+1.0);  // now scale factor is noisy
    //printf("distance = %f\n",distance);
    
    // Choose Position of new molecule
    // first, create vector with length 'distance'
    Vector Position(3);
    for (int i=0;i<3;i++) {
      int sign = 1;
      if (frand() < 0.5) 
	sign = -1;
      Position[i] = sign*frand();
    }
    Position.Normalize();
    Position.Scale(distance);
    // shift this vector relative to the position of our starting atom
    Position += Origin;
    
    // Choose Rotation axis and angle for new molecule
    double Theta = 360.0*frand();
    Vector Rot_axis(3);
    for (int i=0;i<3;i++) {
      int sign = 1;
      if (frand() < 0.5) 
	sign = -1;
      Rot_axis[i] = sign*frand();
    }
    Rot_axis.Normalize();
    
    // Rotate the solvent molecule and translate to its new position
    Monomer New_Solvent(Solvent);
    for (int i=0;i<New_Solvent.GetNumberOfAtoms();i++) {
      Vector tmp = Solvent.GetAtom(i).GetPosition().RotateAboutAxis3D(Theta,Rot_axis);
      tmp += Position;
      New_Solvent.GetAtom(i).SetPosition(tmp);
    }
    
    //printf("Final solvent position:\n");
    //New_Solvent.PrintTinkerCartesian();
    
    // Check distances with all others to make sure we aren't
    // catastrophically close to any other atoms
    int a,b; //scratch variables, unused
    Vector tmp = Unsolvated.FindDistance(New_Solvent);
    double min_dist = tmp[0];
    double cutoff = 1.7;
    if (min_dist < cutoff) {
      printf("New solvent molecule is too close to other atoms (minimum distance = %.2f A.  Trying again.\n",min_dist);
    }
    else {
      printf("New solvent molecule %d accepted.  It is %.2f A from the nearest neighbor\n",n_added+1,min_dist);
      

      // Update the geometry by merging the two monomers
      Monomer newMonomer;
      newMonomer = MergeTwoMonomers(Unsolvated,New_Solvent);
      Natoms = newMonomer.GetNumberOfAtoms(); 
      Unsolvated = newMonomer; // replace original monomer with new one

      n_added++;
    }   
  }

}

double CheckDistances(Vector *New_Solvent, int Natoms_solvent, Vector *XYZ, int Natoms) {
  double min_dist = 100000.0;;

    for (int i=0;i<Natoms_solvent;i++)
      for (int j=0;j<Natoms;j++)
      {
	Vector Diff(New_Solvent[i],true);
	Diff -= XYZ[j];
	double dist = Diff.Norm();
	if (dist < min_dist)
	  min_dist = dist;
      }
  return min_dist;

}

// Floating point pseudorandom number generator: returns number
// between 0 and 1
double frand() {
  double number = rand()/(1.0*RAND_MAX);
  return number;
}

// Write an XYZ file
void WriteXYZfile(string filename, Monomer m, string title) {
  
  FILE *outfile;
  if ((outfile = fopen(filename.c_str(),"w"))==NULL) {
    printf("Solvate: Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  int Natoms = m.GetNumberOfAtoms();
  fprintf(outfile,"%d\n%s\n",Natoms,title.c_str());
  m.PrintMonomerCartesian(outfile);

  fclose(outfile);
  printf("New geometry written to disk\n");
}

// Write a Tinker-style XYZ file
void WriteTinkerXYZfile(string filename, Monomer m, string title) {
  
  FILE *outfile;
  if ((outfile = fopen(filename.c_str(),"w"))==NULL) {
    printf("Solvate: Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  int Natoms = m.GetNumberOfAtoms();
  fprintf(outfile,"%d  %s\n",Natoms,title.c_str());
  m.PrintTinkerCartesian(0,false,outfile);

  fclose(outfile);
  printf("New geometry written to disk\n");
}


Monomer ReadCoordinates(ifstream& infile, int Natoms) {
  string line;

  string *symbols = new string [Natoms];
  double *xyz = new double [3*Natoms];
  int *atom_types = new int [Natoms];
  int *Nconnected = new int [Natoms];
  int *Connectivity = new int [6*Natoms]; // up to 6 connections per atom

  int atom=0;
  while ( atom<Natoms && !infile.eof() ) {
    int ind;
    getline(infile,line);
    istringstream iss(line);
    iss >> ind; // we don't really use the atom index here
    iss >> symbols[atom];
    iss >> xyz[3*atom];
    iss >> xyz[3*atom+1];
    iss >> xyz[3*atom+2];
    iss >> atom_types[atom];

    // Read Connectivity
    int item = 0;
    while (iss >> Connectivity[atom*6+item] ) {
      //printf(" (i=%d) %d ",atom*6+item,Connectivity[atom*6+item]);
      item++;
    }
    //printf("  Total of %d connections\n",item);
    Nconnected[atom] = item; // set number of connections
    atom++; // i

  }

  Monomer Geom;
  int charge = 0; int spin = 1;
string type = "type1";//yoni : added from Kaushik's code: added type to function 
 Geom.Initialize(1,0,charge,spin,type,symbols,xyz,Natoms,atom_types,Nconnected,Connectivity);

  delete [] symbols;
  delete [] xyz;
  delete [] atom_types;
  delete [] Nconnected;
  delete [] Connectivity;

  return Geom;

}

Monomer MergeTwoMonomers(Monomer MonA, Monomer MonB) {
  
  Monomer AB; // to store combined monomer

  // Count number of atoms
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  int Ntot = Na+Nb;

  int charge = MonA.GetChargeState() + MonB.GetChargeState();
  int spin = 1;  // assume singlet.  For now anyways, it doesn't matter.
  // we don't use spin for anything.
  string type = "type1";//yoni : added from Kaushik's code  

  // Arrays to store all key info for combined system
  string *symbols = new string [Ntot];
  double *xyz = new double [3*Ntot];
  int *atom_types = new int [Ntot];
  int *Nconnected = new int [Ntot];
  int *Connectivity = new int [6*Ntot];

  // Grab info from first monomer
  for (int iatom=0;iatom<Na; iatom++) {
    symbols[iatom] = MonA.GetSymbol(iatom);
    for (int q=0;q<3;q++) 
      xyz[3*iatom+q] = MonA.GetAtom(iatom).GetPosition(q);
    atom_types[iatom] = MonA.GetAtom(iatom).GetMMAtomType();
    Nconnected[iatom] = MonA.GetAtom(iatom).GetNumberOfConnectedAtoms();
    for (int q=0;q<Nconnected[iatom];q++) 
      Connectivity[6*iatom+q] = MonA.GetAtom(iatom).GetConnectivity(q);   
  }

  // Grab info from second monomer
  for (int iatom=0;iatom<Nb; iatom++) {
    symbols[iatom+Na] = MonB.GetSymbol(iatom);
    for (int q=0;q<3;q++) 
      xyz[3*(Na+iatom)+q] = MonB.GetAtom(iatom).GetPosition(q);
    atom_types[iatom+Na] = MonB.GetAtom(iatom).GetMMAtomType();
    Nconnected[iatom+Na] = MonB.GetAtom(iatom).GetNumberOfConnectedAtoms();
    for (int q=0;q<Nconnected[iatom+Na];q++) 
      Connectivity[6*(Na+iatom)+q] = MonB.GetAtom(iatom).GetConnectivity(q) + Na;   
  }
 
  //yoni: global index not used because it is not need
  AB.Initialize(1,0,charge,spin,type,symbols,xyz,Ntot,atom_types,
		Nconnected,Connectivity);//yoni : added type from Kaushik's code

  return AB;
}


// Create Water solvent molecule for the amoeba force field
Monomer CreateAmoebaWaterMolecule() {
  // B3LYP/6-311++G(3df,2p) geom, Standard Nuclear Orientation
  int Natoms = 3;
  int charge = 0;
  int spin = 1;
  string type = "type1";//yoni : added from Kaushik's code
  string symbols[3] = {"O","H","H"};
  double xyz[9] = {0.000000,  0.000000, -0.111140, 
		   0.763309,  0.000000,  0.472345,
		   -0.763309,  0.000000,  0.472345};
  int atom_types[3] = {22,23,23};
  int Nconnected[3] = {2,1,1};
  int Connectivity[18];
  Connectivity[0] = 2; Connectivity[1] = 3; // O atom
  Connectivity[6] = 1; // H1 atom
  Connectivity[12] = 1; // H2 atom

  Monomer Water;
  //yoni:Global indexing not used because it is not needed
  Water.Initialize(1,0,charge,spin,type,symbols,xyz,Natoms,atom_types,Nconnected,Connectivity);//yoni : added type from Kaushik's code
  printf("Generic Water:\n");
  Water.PrintTinkerCartesian(); 

  return Water;
}


// Create Benzene solvent molecule for the amoeba force field
Monomer CreateAmoebaBenzeneMolecule() {
  // B3LYP/6-311++G(3df,2p) geom, Standard Nuclear Orientation
  int Natoms = 12;
  int charge = 0;
  int spin = 1;
  string type = "type1";//yoni : added from Kaushik's code
  string symbols[12] = {"C","C","C","H","H","H","C","C","C","H","H","H"};
  double xyz[3*12] = {-0.576018,    1.288201,   -0.025687,
		   -1.060122,    0.317183,    0.878326,
		    0.480171,    0.969821,   -0.908527,
		   -1.023274,    2.285684,   -0.049863,
		   -1.884052,    0.541215,    1.564214,
		    0.845568,    1.714997,   -1.620928,
		    0.577982,   -1.287212,    0.018176,
		    1.061381,   -0.315465,   -0.886083,
		   -0.477586,   -0.967240,    0.900462,
		    1.026328,   -2.283479,    0.040809,
		    1.884678,   -0.538765,   -1.573088,
		   -0.841837,   -1.712579,    1.612116};

  int atom_types[12] = {168,168,168,169,169,169,
			168,168,168,169,169,169};
  int Nconnected[12] = {3,3,3,1,1,1,3,3,3,1,1,1};
  int Connectivity[6*12];
  Connectivity[0] = 2; Connectivity[1] = 3; Connectivity[2] = 4; // C1
  Connectivity[6] = 1; Connectivity[7] = 9; Connectivity[8] = 5; // C2
  Connectivity[12] = 1; Connectivity[13] = 8; Connectivity[14] = 6; // C3
  Connectivity[18] = 1; // H1
  Connectivity[24] = 2; // H2
  Connectivity[30] = 3; // H3
  Connectivity[36] = 8; Connectivity[37] = 9; Connectivity[38] = 10; // C4
  Connectivity[42] = 3; Connectivity[43] = 7; Connectivity[44] = 11; // C5
  Connectivity[48] = 2; Connectivity[49] = 7; Connectivity[50] = 12; // C6
  Connectivity[54] = 7; // H4
  Connectivity[60] = 8; // H5
  Connectivity[66] = 9; // H6


  Monomer Benzene;
//yoni: added type from Kaushik's code
  //yoni:Global indexing not used because it is not needed
  Benzene.Initialize(1,0,charge,spin,type,symbols,xyz,Natoms,atom_types,Nconnected,Connectivity);
  printf("Generic Benzene:\n");
  Benzene.PrintTinkerCartesian();  

  return Benzene;
}
