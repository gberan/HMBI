#include <string>
using std::string;
#include "monomer.h"
#include <sstream>  // by Ali
#include "constants.h"
using namespace hmbi_constants;

Monomer::Monomer() : Monomer_index(0), ref_mon(0), mN_("m0"), spin(1), charge(0), 
		     Natoms(0), Atom_List(NULL), IonizationPotential(0.0),
		     MonomerMass(0.0), QM_Job_Complete(false),
		     MM_Job_Complete(false), Energy_QM(0.0), Energy_MM(0.0),
		     QM_Grad_Init(0), MM_Grad_Init(0), RotationAngle(0.0) 
{
  // The constructor contents got moved to the Initialize member function.
  // empty constructor, just initializes pointers to zero.

  CenterOfMass.Initialize(3);
  RotationVector.Initialize(3);
  

}

Monomer::Monomer(const Monomer &other) : Monomer_index(other.Monomer_index),
					 ref_mon(other.ref_mon),
					 mN_(other.mN_), 
					 spin(other.spin),
					 charge(other.charge), 
					 Natoms(other.Natoms), 
					 Atom_List(NULL), 
					 IonizationPotential(other.IonizationPotential),
					 MonomerMass(other.MonomerMass),
					 QM_Job_Complete(other.QM_Job_Complete),
					 MM_Job_Complete(other.MM_Job_Complete),
					 Energy_QM(other.Energy_QM),
					 Energy_MM(other.Energy_MM),
					 QM_Grad_Init(other.QM_Grad_Init),
					 MM_Grad_Init(other.MM_Grad_Init),
					 RotationAngle(other.RotationAngle) {

  CenterOfMass = other.CenterOfMass;
  //FindCenterOfMass();

  // Copy over Atom_List
  if (Natoms > 0) {
    Atom_List = new Atom[Natoms];
    for (int i=0;i<Natoms;i++)
      Atom_List[i] = other.Atom_List[i];
  }

  // Copy over Gradients
  QM_Grad_Init = other.QM_Grad_Init;
  if (QM_Grad_Init) {
    Grad_QM = other.Grad_QM;
  }
  MM_Grad_Init = other.MM_Grad_Init;
  if (MM_Grad_Init) {
    Grad_MM = other.Grad_MM;
  }
  
  RotationVector = other.RotationVector;
}

Monomer::~Monomer() {
  delete [] Atom_List;
}

// Initialization, w/ no feature for MM details
void Monomer::Initialize(int ind, int charge_in, int spin_in, string atoms[], 
double *xyz, int natoms) {

  Monomer_index = ind;
  ref_mon = ind; // by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel(); 
  charge = charge_in;
  spin = spin_in;
  Natoms = natoms;

  MonomerMass = 0.0;

  // Read in the geometry
  delete [] Atom_List;
  Atom_List = new Atom[Natoms];
  //Atom *TempAtoms = new Atom[Natoms];
  //Atom_List = TempAtoms;
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].Initialize(i+1,atoms[i],xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
    MonomerMass += Atom_List[i].GetAtomicMass();
  }

  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali
}

// Initialization, w/ no feature for MM details
void Monomer::Initialize(int ind, int charge_in, int spin_in, string atoms[], 
			 double *xyz, int natoms, Vector charges) {

  Monomer_index = ind;
  ref_mon = ind; // by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel(); 
  charge = charge_in;
  spin = spin_in;
  Natoms = natoms;

  MonomerMass = 0.0;

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].Initialize(i+1,atoms[i],xyz[3*i], xyz[3*i+1], xyz[3*i+2], charges[i]);
    MonomerMass += Atom_List[i].GetAtomicMass();
  }

  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali
}


// Initialization for QM and Tinker atoms
void Monomer::Initialize(int ind, int charge_in, int spin_in, string atoms[], 
double *xyz, int natoms, int *atom_types, int *Nconnected, 
			 int* connectivity) {

  Monomer_index = ind;
  ref_mon = ind;// by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel();
  charge = charge_in;
  spin = spin_in;
  Natoms = natoms;

  MonomerMass = 0.0;

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  double coord[3];
  for (int i=0;i<Natoms;i++) {
    coord[0] = xyz[3*i];
    coord[1] = xyz[3*i+1];
    coord[2] = xyz[3*i+2];
    Atom_List[i].Initialize(i+1, atoms[i], coord, atom_types[i], 
			    Nconnected[i], &connectivity[6*i]);
    
    MonomerMass += Atom_List[i].GetAtomicMass();
  }

  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali
    
}

// Initialization for QM and Tinker atoms, with embedding point charges
void Monomer::Initialize(int ind, int charge_in, int spin_in, string atoms[], 
double *xyz, int natoms, int *atom_types, int *Nconnected, 
			 int* connectivity, Vector charges) {

  Monomer_index = ind;
  ref_mon = ind;// by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel();
  charge = charge_in;
  spin = spin_in;
  Natoms = natoms;

  MonomerMass = 0.0;

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  double coord[3];
  for (int i=0;i<Natoms;i++) {
    coord[0] = xyz[3*i];
    coord[1] = xyz[3*i+1];
    coord[2] = xyz[3*i+2];
    Atom_List[i].Initialize(i+1, atoms[i], coord, atom_types[i], 
			    Nconnected[i], &connectivity[6*i], charges[i]);
    
    MonomerMass += Atom_List[i].GetAtomicMass();
  }

  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali

}

// set the monomer string label
void Monomer::SetLabel() {
  mN_ = "m";
  char label[10];
  sprintf(label,"%d",Monomer_index);
  mN_ += label;
}

// Find the center of mass
void Monomer::FindCenterOfMass() {
  /*     COM = Sum(m_i*r_i)/Sum(m_i)   */

  CenterOfMass.Set();

  for (int dim=0;dim<3;dim++) { // loop over xyz
    for (int i=0;i<Natoms;i++) { // loop over atoms
      double pos = GetAtom(i).GetCoordinate(dim);
      double mass = GetAtom(i).GetAtomicMass();
      CenterOfMass[dim] += mass*pos/MonomerMass;
    }      
  }

  if ( Params::Parameters().PrintLevel() > 2 ) 
    printf("Monomer %d: Center of Mass = (%f, %f, %f)\n",GetIndex(),
	 CenterOfMass[0],CenterOfMass[1],CenterOfMass[2]);

}

// -------------------------------------------------------
// by Ali --- begin
//
// Find the local coordinates
void Monomer::FindLocalCoord(bool use_old_axes) {

  bool Use_Global_Coords = Params::Parameters().UseGlobalCoordinates();
  bool Use_Ali_Scheme = false; // mostly outdated
  // Ali's scheme does guarantee consistent coordinate systems for
  // 3-atom molecules, though, so we use it for mapping out, e.g. the
  // water force field parameters as a function of geometry.
  if (Params::Parameters().BuildForceFieldOnly())
    Use_Ali_Scheme = true; 

  Vector x_axis(3), y_axis(3), z_axis(3), pos_vec(3);
  Vector RotVec(3);
  double RotAngle;

//if Natoms < 2, use glob coords, since we cannot define local coords with gram schmidt for a monoatomic species
  if (Use_Global_Coords || (Natoms < 2) ) {
    printf("if (Monomer_index==1) Using global coordinates for each monomer.\n");
    // use global cartesian axes (1,0,0), (0,1,0), and (0,0,1)
    x_axis[0] = 1.0;
    y_axis[1] = 1.0;
    z_axis[2] = 1.0;
    
    RotAngle = 0.0;
    // RotVec is an arbitrary vector - (1,0,0)
    RotVec[0] = 1.0; 

  }
  else if (Use_Ali_Scheme) {
    printf("Using Ali's local coordinate scheme\n");
    // First atom is the origin of the local coord system
    // Positive local z axis is from the first to the second atom
    z_axis[0] = GetAtom(1).GetCoordinate(0) - GetAtom(0).GetCoordinate(0);
    z_axis[1] = GetAtom(1).GetCoordinate(1) - GetAtom(0).GetCoordinate(1);
    z_axis[2] = GetAtom(1).GetCoordinate(2) - GetAtom(0).GetCoordinate(2);
    z_axis.Normalize();    
        
    // zx plane is the plane of the first three atoms
    // define a tentative x axis from the first to the third atom
    x_axis[0] = GetAtom(2).GetCoordinate(0) - GetAtom(0).GetCoordinate(0);
    x_axis[1] = GetAtom(2).GetCoordinate(1) - GetAtom(0).GetCoordinate(1);
    x_axis[2] = GetAtom(2).GetCoordinate(2) - GetAtom(0).GetCoordinate(2);
    
    // y axis is the z cross tentative x
    y_axis[0] = z_axis[1]*x_axis[2] - z_axis[2]*x_axis[1];
    y_axis[1] = z_axis[2]*x_axis[0] - z_axis[0]*x_axis[2];
    y_axis[2] = z_axis[0]*x_axis[1] - z_axis[1]*x_axis[0];
    y_axis.Normalize();    
    
    // true x axis is the cross product of y and z
    x_axis[0] = y_axis[1]*z_axis[2] - y_axis[2]*z_axis[1];
    x_axis[1] = y_axis[2]*z_axis[0] - y_axis[0]*z_axis[2];
    x_axis[2] = y_axis[0]*z_axis[1] - y_axis[1]*z_axis[0];
    x_axis.Normalize(); 

    /*
    printf("Local coordinate axes:\n");
    printf("X: %10.6f %10.6f %10.6f\n",x_axis[0],x_axis[1],x_axis[2]);
    printf("Y: %10.6f %10.6f %10.6f\n",y_axis[0],y_axis[1],y_axis[2]);
    printf("Z: %10.6f %10.6f %10.6f\n",z_axis[0],z_axis[1],z_axis[2]);
    */

    // To rotate the global coord sys to the local one (for Orient program), 
    // we need a rotation angle,  
    // and a vector about which the rotation takes place.    
    RotAngle = acos((x_axis[0]+y_axis[1]+z_axis[2]-1)/2);
    //if (Params::Parameters().PrintLevel() > 0) 
    //printf("Rotation Angle (in radian)     %10.6f          \n", RotAngle);  	  
    
    RotAngle = RotAngle*RadiansToDegrees;
    //if (Params::Parameters().PrintLevel() > 0) 
    printf("Rotation Angle (in degrees)     %10.6f          \n", RotAngle);  	  
    
    double De;
    De = sqrt((y_axis[2]-z_axis[1])*(y_axis[2]-z_axis[1])+(z_axis[0]-x_axis[2])*(z_axis[0]-x_axis[2])+(x_axis[1]-y_axis[0])*(x_axis[1]-y_axis[0]));
    if (De==0.0)
      De=0.000000001;
    
    
    RotVec[0] = (y_axis[2]-z_axis[1])/De; 
    RotVec[1] = (z_axis[0]-x_axis[2])/De; 
    RotVec[2] = (x_axis[1]-y_axis[0])/De;

  }
  else {
    if (Monomer_index==1) printf("Finding local coordinates for each monomer.\n");

    //Gram-Schmidt's orthogonalization process - added on 4th Feb, 2010 (kdn)
    
    /* explanation */
    
    // z-axis will be the vector joining the two atoms.  x-axis will
    // be the found by projecting the position vector of the 1st atom
    // onto the z-axis cross product of the two axes will give the
    // y-axis.  remember the right-hand rule for cross product.
    
    //Start

    //Forming normalized z-axis
    z_axis[0] = GetAtom(1).GetCoordinate(0) - GetAtom(0).GetCoordinate(0);
    z_axis[1] = GetAtom(1).GetCoordinate(1) - GetAtom(0).GetCoordinate(1);
    z_axis[2] = GetAtom(1).GetCoordinate(2) - GetAtom(0).GetCoordinate(2);
    z_axis.Normalize();


    //Forming normalized x-axis using Gram-Schmidt

    // pos_vec = position vector of the 0th atom --> not essential anymore
    pos_vec[0] = GetAtom(0).GetCoordinate(0);
    pos_vec[1] = GetAtom(0).GetCoordinate(1);
    pos_vec[2] = GetAtom(0).GetCoordinate(2);

    /* GJB I commented this out based on Shuhao's version.  Hope that was correct!
    //Forming normalized x-axis using Gram-Schmidt if the 1st atom is
    //origin or if the origin and the two atoms fall on the same line
    //(crossproduct =0), generate and use a random vector

    if ( (fabs(pos_vec.CrossProduct(z_axis)[0])<0.0001) && (fabs(pos_vec.CrossProduct(z_axis)[1])<0.0001) && (fabs(pos_vec.CrossProduct(z_axis)[2])<0.0001) ) {
      pos_vec[0] = rand();
      pos_vec[1] = rand(); 
      pos_vec[2] = rand(); 
      pos_vec.Normalize();
      //pos_vec.Print("pos_vec_random");

    //random vec should be printed out for reference??
    //  pos_vec.Print("random position vector used");
     } 
    */
    //Forming normalized x-axis using Gram-Schmidt if the 1st atom is
    //origin or if the origin and the two atoms fall on the same line
    //(dotproduct =0), generate and use a random vector
   if ( pos_vec.DotProduct(z_axis) == 0.00 ) {
      pos_vec[0] = rand();
      pos_vec[1] = rand();
      pos_vec[2] = rand();
      pos_vec.Normalize();
      //pos_vec.Print("pos_vec_random");

      //random vec should be printed out for reference??
      //  pos_vec.Print("random position vector used");
   }

    // using the projection of pos_vec onto the z-axis to find x-axis
    x_axis = z_axis;
    x_axis *= -1*pos_vec.DotProduct(z_axis)/ z_axis.DotProduct(z_axis);
    x_axis += pos_vec;

    // y axis is the z cross tentative x
    y_axis[0] = z_axis[1]*x_axis[2] - z_axis[2]*x_axis[1];
    y_axis[1] = z_axis[2]*x_axis[0] - z_axis[0]*x_axis[2];
    y_axis[2] = z_axis[0]*x_axis[1] - z_axis[1]*x_axis[0];
    y_axis.Normalize();

    // true x axis is the cross product of y and z (it actually gives back the same x-axis formed above)
    x_axis[0] = y_axis[1]*z_axis[2] - y_axis[2]*z_axis[1];
    x_axis[1] = y_axis[2]*z_axis[0] - y_axis[0]*z_axis[2];
    x_axis[2] = y_axis[0]*z_axis[1] - y_axis[1]*z_axis[0];
    x_axis.Normalize();

    //printf("Local coordinate axes:\n");
    //printf("X: %10.6f %10.6f %10.6f\n",x_axis[0],x_axis[1],x_axis[2]);
    //printf("Y: %10.6f %10.6f %10.6f\n",y_axis[0],y_axis[1],y_axis[2]);
    //printf("Z: %10.6f %10.6f %10.6f\n",z_axis[0],z_axis[1],z_axis[2]);

    //End : Gram-Schmidt

    // To rotate the global coord sys to the local one (for Orient program), 
    // we need a rotation angle,  
    // and a vector about which the rotation takes place.    
    RotAngle = acos((x_axis[0]+y_axis[1]+z_axis[2]-1)/2);
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Rotation Angle (in radian)     %10.6f          \n", RotAngle);  	  
    
    RotAngle = RotAngle*RadiansToDegrees;
    if (Params::Parameters().PrintLevel() > 0) 
      printf("Rotation Angle (in degrees)     %10.6f          \n", RotAngle);  	  
    
    double De;
    De = sqrt((y_axis[2]-z_axis[1])*(y_axis[2]-z_axis[1])+(z_axis[0]-x_axis[2])*(z_axis[0]-x_axis[2])+(x_axis[1]-y_axis[0])*(x_axis[1]-y_axis[0]));
    if (De==0.0)
      De=0.000000001;
    
    
    RotVec[0] = (y_axis[2]-z_axis[1])/De; 
    RotVec[1] = (z_axis[0]-x_axis[2])/De; 
    RotVec[2] = (x_axis[1]-y_axis[0])/De;
  }
  
  // If this is our first time here, store the coordinates
  //  if ( !use_old_axes) {
  SetRotationAngle(RotAngle); 
  RotationVector = RotVec;
  /*  }
    else { // after first energy calculation, use ones previously stored
    printf("Using previously stored local coordinates for monomer %d\n",Monomer_index);
    RotAngle = RotationAngle;
    RotVec = RotationVector;
    }
  */

  if (Params::Parameters().PrintLevel() > 0) 
    printf("Rotation vector      %10.6f          %10.6f           %10.6f          \n", 
	   RotVec[0], RotVec[1], RotVec[2]);  	  
  
  // local coordinates in local coord system
  double loc_xyz[3*Natoms];
  
  // local coordinates in global coord system
  double loc_xyz_g[3*Natoms]; 

  if (Use_Global_Coords) {
    for (int i=0;i<Natoms;i++) {
      loc_xyz[3*i] = GetAtom(i).GetCoordinate(0);
      loc_xyz[3*i+1] = GetAtom(i).GetCoordinate(1);
      loc_xyz[3*i+2] = GetAtom(i).GetCoordinate(2);
    }
  }
  else {

    loc_xyz[0] = 0.000000;  // local coords of the first atom
    loc_xyz[1] = 0.000000;
    loc_xyz[2] = 0.000000;
    
    for (int i=1;i<Natoms;i++) {
      
      // local x in global coord sys
      loc_xyz_g[3*i] = GetAtom(i).GetCoordinate(0) - GetAtom(0).GetCoordinate(0);
      
      // local y in global coord sys
      loc_xyz_g[3*i+1] = GetAtom(i).GetCoordinate(1) - GetAtom(0).GetCoordinate(1);
      
      // local z in global coord sys
      loc_xyz_g[3*i+2] = GetAtom(i).GetCoordinate(2) - GetAtom(0).GetCoordinate(2);
      
      //local x in local coord sys (dot product of x_axis and local x in global coord)
      loc_xyz[3*i] = x_axis[0]*loc_xyz_g[3*i] + x_axis[1]*loc_xyz_g[3*i+1] + x_axis[2]*loc_xyz_g[3*i+2];
      
      //local y in local coord sys (dot product of y_axis and local y in global coord)
      loc_xyz[3*i+1] = y_axis[0]*loc_xyz_g[3*i] + y_axis[1]*loc_xyz_g[3*i+1] + y_axis[2]*loc_xyz_g[3*i+2];
      
      //local z in local coord sys (dot product of z_axis and local z in global coord)
      loc_xyz[3*i+2] = z_axis[0]*loc_xyz_g[3*i] + z_axis[1]*loc_xyz_g[3*i+1] + z_axis[2]*loc_xyz_g[3*i+2];
      
    }    
  }
  //const double AngsToBohr = 1.889725989;

  for (int i=0;i<Natoms;i++) {
    double loc[3];
    loc[0] = loc_xyz[3*i]; 
    loc[1] = loc_xyz[3*i+1]; 
    loc[2] = loc_xyz[3*i+2]; 
    GetAtom(i).SetLocalPosition(loc);

    if (Params::Parameters().PrintLevel() > 0) 
      printf("%s      %10.6f   %10.6f   %10.6f  Angstroms \n", 
	     GetAtom(i).GetSymbol().c_str(), 
	     GetAtom(i).GetLocalPosition(0),
	     GetAtom(i).GetLocalPosition(1),
	     GetAtom(i).GetLocalPosition(2));
  }

}


// by Ali --- end
// ------------------------------------------



double Monomer::ReadQChemEnergy(bool MM_job) {

  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  
  string out_filename = path + "/" + mN_ + ".out"; 
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::SetQChemEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }
  
  double nuc_chg, chg_chg; // for if using embedding charges

  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++)
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }

    // Search for extra energies that arise when using embedding charges
    if ( match==" Nucleus-charge energy" ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> nuc_chg; // read energy 
    }

    // Search for extra energies that arise when using embedding charges
    if ( match==" Charge-charge energy " ) {
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> chg_chg; // read energy 
    }
  }

  if ( Params::Parameters().UseEmbeddingCharges() ) {
    double old = energy;
    energy = energy - (nuc_chg + chg_chg);

    if ( Params::Parameters().PrintLevel() > 0 ) {
      printf("Embedding charges in use.  Adjusting the final energy\n");
      printf("Original:  %15.9f\n",old);
      printf("Nuc-chg =  %15.9f, Chg-chg = %15.9f\n",nuc_chg,chg_chg);
      printf("Corrected: %15.9f\n",energy);
    }
  }

  if ( Params::Parameters().PrintLevel() > 0) printf("QM Energy = %15.9f\n",energy);
 
  // Close the force file
  infile.close();
  return energy;
}

void Monomer::SetMMEnergy() {
  double energy;
  if ( Params::Parameters().GetMMType()==1 ) // Tinker
    energy = ReadTinkerEnergy();
  else if ( Params::Parameters().GetMMType()==2 ) // AIFF
    energy = 0.0; // AIFF computes only 2-body & higher terms.  
  else if ( Params::Parameters().GetMMType()==3) // QChem
    energy = ReadQChemEnergy(true);
  else {
    printf("Monomer::SetMMEnergy: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  }

  Energy_MM = energy;
  //printf("MM Energy = %15.9f\n",Energy_MM);

}

double Monomer::ReadTinkerEnergy() {
  double energy ;

  // if embedding charges, need to read polarization energy line
  //if ( Params::Parameters().UseEmbeddingCharges() ) {
  if ( 1 ) {
   
    // Set up the filename, with the full path.  File is e.g. 'm3.out'
    string out_filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".out"; 

    // Open the energy file
    ifstream infile;
    infile.open( out_filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadTinkerEnergy : Cannot open file '%s'\n",out_filename.c_str());
      exit(1);
    }
    
    // Read in the data
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);

      string match = line.substr(0,23);
      if ( match==" Total Potential Energy" ) {
	istringstream iss(line);
	string tmp; 
	// throw away "text tags"
	iss >> tmp; iss >> tmp; iss >> tmp; iss>> tmp; 
	iss >> energy; // read the energy
      }
      

    }
    
    if ( Params::Parameters().PrintLevel() > 0 ) 
      printf("Tinker Energy = %f kcal/mol\n",energy);
    
    energy /= HartreesToKcalpermole; // convert to hartrees
    
    // Close the force file
    infile.close();
    }

  else {
    // Polarization energy = 0 for a monomer without embedding
    double energy = 0.0;
  }

  return energy;
}

// Read the atom-centered distributed multipole moments, as computed
// by CamCasp
void Monomer::ReadMultipoleMoments() {
  int Natoms_found = 0;

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".mom"; 

  // Open the moments file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadMultipoleMoments(): Cannot open monomer %d multipole file '%s'",
	   Monomer_index,filename.c_str());
    exit(1);
  }
  
  // Create a list of atom names which is used to track our progress as
  // we read the moments file
  string AtmInt[Natoms];
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
  
  // Read in moments file
  string line;
  int i=0;
  string entryAtm;
  int nmom, rank;
    
  while ( !infile.eof() ) {
    int iatom;
    getline(infile,line);
    i++;
    if (i>2){
      
      istringstream iss(line);
      string entry;
      iss >> entry;
      bool AtomFound = false;
      for(int n=0;n<Natoms;n++) {
	if (StringToUpper(entry)==StringToUpper(AtmInt[n])) {
	  AtomFound = true;
	  entryAtm = entry;
	  iatom = n;
	}
      }
      // Starting a new atom.
      if (AtomFound) {
	Natoms_found++;

	// Parse the header line to grab the rank.
	iss >> entry; // garbage
	iss >> entry; // garbage
	iss >> entry; // garbage
	iss >> entry; // garbage
	iss >> rank;
	

	// figure out how many moments/lines we need to read
	if (rank == 0) { nmom = 1;}
	else if (rank == 1) {nmom = 4;}
	else if (rank == 2) {nmom = 9;}
	else if (rank == 3) {nmom = 16;}
	else if (rank == 4) {nmom = 25;}
	else {
	  printf("Error: Monomer::ReadMultipoleMoments(): Invalid rank (%d) for Monomer %d\n",rank,Monomer_index);
	  exit(1);
	}

	double *moments =  new double[nmom];
	for (int imom=0;imom<nmom;imom++) {
	  infile >> moments[imom];
	}

	// Initialize the Multipole object & attach it to an atom
	Multipole Moments(rank, moments);
	GetAtom(iatom).SetMultipoleMoments(Moments);

	if (Params::Parameters().PrintLevel() > 0) {
	  // Print out the multipole moments we just read
	  char mon_ind[10];
	  sprintf(mon_ind,"%d",Monomer_index);
	  string title = "Moments for Monomer ";
	  title += mon_ind;
	  title += ", Atom " + AtmInt[iatom];
	  Moments.Print(title);
	}

	delete [] moments;

      } // end if (AtomFound)
    }  // end if (i>2) 
  } // end while loop

  infile.close();


  if (Natoms_found != Natoms) {
    printf("ERROR: Monomer::ReadMultipoleMoments(): Only found %d of %d expected atoms for monomer %d.\n",Natoms_found,Natoms,Monomer_index);
    exit(1);
  }
}

// Read the distributed polarizabilities for the AIFF, as computed by
// CamCasp.
void Monomer::ReadPolarizabilities() {

  int Natoms_found = 0;

  int rank = 2; // for now, default to rank 2.  Need a mechanism to
		// handle higher ranks.
  int npol;
  if (rank == 0) { npol = 1;}
  else if (rank == 1) {npol = 4;}
  else if (rank == 2) {npol = 9;}
  else if (rank == 3) {npol = 16;}
  else if (rank == 4) {npol = 25;}

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".pol"; 

  // Open the polarizabilities file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadPolarizabilities(): For Monomer %d, cannot open file '%s'\n",
	   Monomer_index,filename.c_str());
    exit(1);
  }
 
  // Create a list of atom names which is used to track our progress as
  // we read the polarizability file
  string AtmInt[Natoms];
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
    
  // Read in the polarizabilities
  string line;
  int iline = 0;
  string entryAtm;
  while ( !infile.eof() ) {
    iline++;
    int iatom;
    getline(infile,line);
    
    istringstream iss(line);
    string entry;
    iss >> entry;
    bool AtomFound = false;
    bool entryBool = false;
    for (int n=0;n<Natoms;n++) {
      if (StringToUpper(entry)==StringToUpper(AtmInt[n])) {
	AtomFound = true;
	entryAtm = entry;
	iatom = n;
      }
    }
    if (AtomFound){
      Natoms_found++;

      // Now start reading the polarizabilities.
      // they are an (npol x npol) array, initially get stored in 
      // a long vector in row order.  
      double *pols =  new double[npol*npol]; 
      for (int ipol=0;ipol<npol*npol;ipol++) {
	infile >> pols[ipol];
      }

      // Now transpose to put it in column order.. col1, col2, col3, etc.
      double *polsT = new double[npol*npol];
      for (int i=0;i<npol;i++) 
	for (int j=0;j<npol;j++) 
	  polsT[i+j*npol] = pols[j + i*npol];

      // Initialize the atomic polarizabilities
      Polarizability Pols(rank,polsT);
      GetAtom(iatom).SetPolarizability(Pols);

      delete [] pols;
      delete [] polsT;            
    } 
  }

  infile.close();  

  if (Natoms_found != Natoms) {
    printf("ERROR: Monomer::ReadPolarizabilities(): Only found %d of %d expected atoms for monomer %d.\n",Natoms_found,Natoms,Monomer_index);
    exit(1);
  }

}

// Read the frequency-dependent uncoupled polarizabilities as computed
// by CamCasp.  These get used to compute dispersion coefficients.
void Monomer::ReadFreqPolarizabilities() {
  // We typically use freq-dependent polarizabilities that have been
  // evaluated at 10 frequencies.
  int nfreq = 10;
  Vector FreqPoldipole(nfreq);
  Vector FreqPolquad(nfreq);
  int Natoms_found = 0;

  printf("  Reading isotropic frequency-dependent polarizabilities.\n");

  
  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".freq_pol"; 
  
  // Open the frequency-dependent polarizabilities file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadFreqPolarizabilities(): Cannot open monomer %d frequency-dependent polarizability file '%s'",Monomer_index, filename.c_str());
    exit(1);
  }
  
  // Create a list of atom names which is used to track our progress
  // as we read through the polarizability file
  string AtmInt[Natoms]; // define which atoms we will read
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
  
  // Read in frequency-dependent polarizability file
  string line; // define a string variable used by "getline" function
  int i=0;
  string entryAtm; // define for printing 
  
  while ( !infile.eof() ) {
    int iatom;
    string text = "Site";    
    
    getline(infile,line);// read by line
    i++;
    if (i>4){     
      istringstream iss(line);
      string site;
      string atom;
      
      iss >> site; // site
      iss >> atom; // atom 
      bool AtomFound = false;
      for(int n=0;n<Natoms;n++) {
	if (site == text && StringToUpper(atom)==StringToUpper(AtmInt[n]) ) {
	  AtomFound = true;
	  entryAtm = atom;
	  iatom = n;
	}
      }
      // We have found an atom in the file, now read the polarizabilities
      if (AtomFound) {
	Natoms_found++;
	FreqPoldipole.Set();
	FreqPolquad.Set();
	// First, read in the frequency-dependent isotropic
	// dipole-dipole polarizability
	getline(infile,line);
	istringstream iss(line);
	string rank1;
	string rank2;
	string text1 = "10"; // dipole component, l=1, m=0
	iss >> rank1;
	iss >> rank2;
	
	if (rank1 == text1 && rank2 == text1) {
	  double *freqpoldipole = new double[nfreq];
	  for (int ifreq=0;ifreq<nfreq;ifreq++) {
	    infile >> freqpoldipole[ifreq];
	  }
	  // Attach the dipole polarizability data to the atom
	  FreqPoldipole.Initialize(freqpoldipole,nfreq);
	  GetAtom(iatom).SetFreq_Pol_Dipole(FreqPoldipole); 
	  delete [] freqpoldipole; 
 
	  // Now read a bunch of lines to jump ahead to the quadrupole
	  // components
	  getline(infile,line); 
	  getline(infile,line);    
	  getline(infile,line);    
	  getline(infile,line);     
	  getline(infile,line);    
	  getline(infile,line);   
	  getline(infile,line);    
	  getline(infile,line);
	  
	  // Read in the frequency-dependent isotropic
	  // quadrupole-quadrupole polarizability.  Note, if we didn't
	  // go to high enough rank, these will be zero and are not
	  // present in the polarizability file.  In that case, the if
	  // statement is never true.

	  istringstream iss(line); // istringstream function for reading word by word
	  string text2 = "20"; // quadrupole component, l=2, m=0

	  string rank3;
	  string rank4;
	  
	  iss >> rank3;
	  iss >> rank4;
	  
	  if (rank3 == text2 && rank4 == text2) {
	    double *freqpolquad = new double[nfreq];
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpolquad[ifreq];
	    }

	    // Attach the quadrupole polarizability data to the atom
	    FreqPolquad.Initialize(freqpolquad,nfreq);
	    GetAtom(iatom).SetFreq_Pol_Quad(FreqPolquad);
	    delete [] freqpolquad;
	  } // end if rank3 == text2 && rank4 == text2
	}// end if rank1 == text1 && rank2 == text1  

	if (Params::Parameters().PrintLevel() > 0) {
	  // Print out the frequency-dependent polarizability we just read
	  char mon_ind[10];
	  sprintf(mon_ind,"%d",Monomer_index);
	  string title1 = "Frequency dependent dipole polarizability for Monomer";
	  title1 += mon_ind;
	  title1 += ", Atom " + AtmInt[iatom];
	  FreqPoldipole.Print(title1);

	  string title2 = "Frequency dependent quadrupole polarizability for Monomer ";
	  title2 += mon_ind;
	  title2 += ", Atom " + AtmInt[iatom];	 
	  FreqPolquad.Print(title2);
	}

      } // end if (AtomFound) 
    } // end if (i>4) 
  } // end while loop 
  
  infile.close(); 

  if (Natoms_found != Natoms) {
    printf("ERROR: Monomer::ReadFreqPolarizabilities(): Only found %d of %d expected atoms for monomer %d.\n",Natoms_found,Natoms,Monomer_index);
    exit(1);
  }

}

// Read the diagonal anisotropic frequency-dependent uncoupled
// polarizabilities as computed by CamCasp, then average them to get
// isotropic ones.  These get used to compute dispersion coefficients.
void Monomer::ReadDiagonalAnisotropicFreqPolarizabilities() {

  printf("  Reading diagonal anisotropic frequency-dependent polarizabilities and averaging them.\n");

  // We typically use freq-dependent polarizabilities that have been
  // evaluated at 10 frequencies.
  int nfreq = 10;

  // These next arrays store the final isotropic polarizabilities
  Vector FreqDepDipolePol(nfreq);
  Vector FreqDepQuadrupolePol(nfreq);

  // counter used to make sure we find data for each atom in the monomer
  int Natoms_found = 0;

  // Temporary arrays to store the 3 diagonal components of the
  // anisotropic dipole-dipole polarizability at each frequency:
  Vector freqpol_10(nfreq);
  Vector freqpol_11s(nfreq);
  Vector freqpol_11c(nfreq);
  // and the quadrupole-quadrupole components:
  Vector freqpol_20(nfreq);
  Vector freqpol_21c(nfreq);
  Vector freqpol_21s(nfreq);
  Vector freqpol_22c(nfreq);
  Vector freqpol_22s(nfreq);

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".freq_pol"; 
  
  // Open the frequency-dependent polarizabilities file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadFreqPolarizabilities(): Cannot open monomer %d frequency-dependent polarizability file '%s'",Monomer_index, filename.c_str());
    exit(1);
  }
  
  // Create a list of atom names which is used to track our progress
  // as we read through the polarizability file
  string AtmInt[Natoms]; // define which atoms we will read
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
  
  // Read in frequency-dependent polarizability file
  string line; // define a string variable used by "getline" function
  int i=0;
  string entryAtm; // define for printing 
  
  while ( !infile.eof() ) {
    int iatom;
    string text = "Site";    
    
    getline(infile,line);// read by line
    i++;
    if (i>4){     
      istringstream iss(line);
      string site;
      string atom;
      
      iss >> site; // site
      iss >> atom; // atom 
      bool AtomFound = false;
      for(int n=0;n<Natoms;n++) {
	if (site == text && StringToUpper(atom)==StringToUpper(AtmInt[n]) ) {
	  AtomFound = true;
	  entryAtm = atom;
	  iatom = n;
	}
      }
      // We have found an atom in the file, now read the polarizabilities
      if (AtomFound) {
	Natoms_found++;

	// Clean up the arrays before we start
	FreqDepDipolePol.Set();
	FreqDepQuadrupolePol.Set();

	freqpol_10.Set();
	freqpol_11c.Set();
	freqpol_11s.Set();
	freqpol_20.Set();
	freqpol_21c.Set();
	freqpol_21s.Set();
	freqpol_22c.Set();
	freqpol_22s.Set();

	string r10 = "10", r11c = "11c", r11s = "11s";
	string r20 = "20", r21c = "21c", r21s = "21s", r22c = "22c", r22s = "22s";

	bool read = true;
	while (read == true) {
	  getline(infile,line);

	  // Check if we are done with this atom
	  if (line=="End") {
	    read = false;
	    break;
	  }
	  
	  // Otherwise, check if we have a diagonal polarizability element
	  istringstream iss(line);
	  string rank1, rank2;
	  iss >> rank1;
	  iss >> rank2;
	  
	  //printf("rank1 = %s, rank2 = %s\n",rank1.c_str(), rank2.c_str());

	  // Dipole polarizabilities:

	  // Case: 10 10 component
	  if (rank1==r10 && rank2==r10) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_10[ifreq];
	    }
	    //freqpol_10.Print("10 10 components:");
	  }

	  // Case: 11c 11c component
	  if (rank1==r11c && rank2==r11c) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_11c[ifreq];
	    }
	    //freqpol_11c.Print("11c 11c components:");
	  }

	  // Case: 11s 11s component
	  if (rank1==r11s && rank2==r11s) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_11s[ifreq];
	    }
	    //freqpol_11s.Print("11s 11s components:");
	  }

	  // And quadrupole polarizabilities:

	  // Case: 20 20 component
	  if (rank1==r20 && rank2==r20) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_20[ifreq];
	    }
	    //freqpol_20.Print("20 20 components:");
	  }

	  // Case: 21c 21c component
	  if (rank1==r21c && rank2==r21c) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_21c[ifreq];
	    }
	    //freqpol_21c.Print("21c 21c components:");
	  }

	  // Case: 21s 21s component
	  if (rank1==r21s && rank2==r21s) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_21s[ifreq];
	    }
	    //freqpol_21s.Print("21c 21c components:");
	  }

	  // Case: 22c 22c component
	  if (rank1==r22c && rank2==r22c) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_22c[ifreq];
	    }
	    //freqpol_22c.Print("22c 22c components:");
	  }

	  // Case: 22s 22s component
	  if (rank1==r22s && rank2==r22s) {
	    for (int ifreq=0;ifreq<nfreq;ifreq++) {
	      infile >> freqpol_22s[ifreq];
	    }
	    //freqpol_22s.Print("22c 22c components:");
	  }


	} // end while loop (reading relevant pols for this atom)

	// Now average them to get isotropic quantities
	// Dipole-dipole:
	FreqDepDipolePol = freqpol_10;
	FreqDepDipolePol += freqpol_11c;
	FreqDepDipolePol += freqpol_11s;
	FreqDepDipolePol.Scale(1.0/3.0);

	
	// Quadrupole-quadrupole:
	FreqDepQuadrupolePol = freqpol_20;
	FreqDepQuadrupolePol += freqpol_21c;
	FreqDepQuadrupolePol += freqpol_21s;
	FreqDepQuadrupolePol += freqpol_22c;
	FreqDepQuadrupolePol += freqpol_22s;
	FreqDepQuadrupolePol.Scale(1.0/5.0);

	// Attach the isotropic values to the atom
	GetAtom(iatom).SetFreq_Pol_Dipole(FreqDepDipolePol); 
	GetAtom(iatom).SetFreq_Pol_Quad(FreqDepQuadrupolePol);

	// Optional printing
	if (Params::Parameters().PrintLevel() > 0) {
	  printf("Monomer %d, atom %s\n",Monomer_index, atom.c_str() );
	  FreqDepDipolePol.Print("Isotropic frequency-dependent dipole polarizability");
	  FreqDepQuadrupolePol.Print("Isotropic frequency-dependent quadrupole polarizability");
	}

	  
      } // end if (AtomFound) 
    } // end if (i>4) 
  } // end while loop 
  
  infile.close(); 

  if (Natoms_found != Natoms) {
    printf("ERROR: Monomer::ReadDiagonalAnisotropicFreqPolarizabilities(): Only found %d of %d expected atoms for monomer %d.\n",Natoms_found,Natoms,Monomer_index);
    exit(1);
  }

}
  

// Read the isotropic dipole polarizability for the AIFF, as computed
// by CamCasp.  It gets used to compute the C9 triple dipole 3-body
// dispersion term.  
//
// This function is now obsolete, since we calculate the dispersion
// coefficients via Casimir-Polder integration.
void Monomer::ReadIsotropicDipolePolarizability() {

  int Natoms_found = 0;

  double isotropic_dipole_polarizability = 0.0;

  int rank = 2; // for now, default to rank 2.  Need a mechanism to
		// handle higher ranks.
  int npol;
  if (rank == 0) { npol = 1;}
  else if (rank == 1) {npol = 4;}
  else if (rank == 2) {npol = 9;}
  else if (rank == 3) {npol = 16;}
  else if (rank == 4) {npol = 25;}

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".iso_pol"; 

  // Open the polarizabilities file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadIsotropicDipolePolarizability(): For Monomer %d, cannot open file '%s'\n",
	   Monomer_index,filename.c_str());
    exit(1);
  }
 
  // Create a list of atom names which is used to track our progress as
  // we read the polarizability file
  string AtmInt[Natoms];
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
    
  // Read in the polarizabilities
  string line;
  int iline = 0;
  string entryAtm;
  while ( !infile.eof() ) {
    iline++;
    int iatom;
    getline(infile,line);
    
    istringstream iss(line);
    string entry;
    iss >> entry;
    bool AtomFound = false;
    bool entryBool = false;
    for (int n=0;n<Natoms;n++) {
      if (StringToUpper(entry)==StringToUpper(AtmInt[n])) {
	AtomFound = true;
	entryAtm = entry;
	iatom = n;
      }
    }
    if (AtomFound){
      Natoms_found++;
      // Now start reading the polarizabilities.  they are in an (npol
      // x npol) array.  We want the element (2,2).
      double scratch;
      for (int ipol=0;ipol<npol*npol;ipol++) {
	if (ipol == npol + 1)
	  infile >> isotropic_dipole_polarizability;
	else
	  infile >> scratch;
      }

      printf("Isotropic Dipole Polarizability = %f\n",isotropic_dipole_polarizability);

      // Initialize the atomic isotropic dipole polarizability
      GetAtom(iatom).SetIsotropicDipolePolarizability(isotropic_dipole_polarizability);
    } 
  }
  infile.close();  


  if (Natoms_found != Natoms) {
    printf("ERROR: Monomer::ReadFreqPolarizabilities(): Only found %d of %d expected atoms for monomer %d.\n",Natoms_found,Natoms,Monomer_index);
    exit(1);
  }

}

// This function was used in an older version of the AIFF to read in
// dispersion coefficients form CamCasp.  We now calculate these coefficients
// directly from the frequency-dependent polarizabilities.  So this
// function is no longer needed.
void Monomer::ReadDispersionCoefficients() {

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".disp"; 

  // Open the moments file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadDispersionCoefficients(): Cannot open dispersion coefficient file '%s'",
	   filename.c_str());
    exit(1);
  }
  
  // Create a list of atom names which is used to track our progress as
  // we read the moments file
  string AtmInt[Natoms];
  for (int i=0;i<Natoms;i++){
    char label[10];
    AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
  }
  
  // Read in dispersion coefficients file
  string line;
  int i=0;
  string entryAtm;
  
  printf("\nAb initio two-body dispersion coefficients for monomer %d:\n",
	 Monomer_index);
  while ( !infile.eof() ) {
    int iatom;
    double C6=0.0; double C8=0.0; double C10=0.0;

    // Parse atom name line
    getline(infile,line);
    if (strlen(line.c_str())==0) break;
    istringstream iss(line);
    string atom1, atom2;
    iss >> atom1; // atom 1
    iss >> atom2; // atom 2
    //printf("atom1 = %s, atom2 = %s\n", atom1.c_str(), atom2.c_str());
    if (atom1 == atom2 && StringToUpper(AtmInt[i]) == StringToUpper(atom1)) {
      //printf("MATCH! Found atom %s: %s %s\n",AtmInt[i].c_str(),atom1.c_str(),atom2.c_str());
      getline(infile,line);

      istringstream disp(line);
      string junk;
      disp >> junk; disp >> junk; disp >> junk; // ignore some initial fields
      disp >> C6;
      disp >> junk;  // if this fails (i.e. at end of line), disp.eof() set to TRUE
      if (! disp.eof()) { // if there are more entries on the line still
	disp >> C8;
	disp >> junk;
	disp >> C10;
      }
      GetAtom(i).SetDispersionCoefficients(C6,C8,C10);
      GetAtom(i).PrintDispersionCoefficients();

      // For now, we use empirical damping factors that depend on the vdW radii.
      // To specify these, we use some empirical atom types.  These are set here.
      // Future work should try to do better to model these terms.

      // Determine the dispersion atom type
      string atom_type;
      if (GetAtom(i).GetSymbol() == "H") atom_type = "Hs";
      else if (GetAtom(i).GetSymbol() == "C") atom_type = "Csp2";
      else if (GetAtom(i).GetSymbol() == "N") atom_type = "Nsp2sp3";
      else if (GetAtom(i).GetSymbol() == "O") atom_type = "Osp3";
      else if (GetAtom(i).GetSymbol() == "He") atom_type = "He";
      else if (GetAtom(i).GetSymbol() == "Ne") atom_type = "Ne";
      else if (GetAtom(i).GetSymbol() == "Ar") atom_type = "Ar";
      else if (GetAtom(i).GetSymbol() == "Kr") atom_type = "Kr";
      else {
	printf("Monomer::ReadDispersionCoefficients: Unknown Atom Type %s\n",GetAtom(i).GetSymbol().c_str());
	exit(1);
      }
      GetAtom(i).SetDispersionAtomType( atom_type );
      i++;
    } // end while loop 
  }
  infile.close();

}

// Give each atom a type
void Monomer::SetEmpiricalAtomDispersionType() {

  for (int i=0;i<Natoms;i++) {
    string atom_type;
      if (GetAtom(i).GetSymbol() == "H") atom_type = "Hs";
      else if (GetAtom(i).GetSymbol() == "C") atom_type = "Csp2";
      else if (GetAtom(i).GetSymbol() == "N") atom_type = "Nsp2sp3";
      else if (GetAtom(i).GetSymbol() == "O") atom_type = "Osp3";
      else if (GetAtom(i).GetSymbol() == "He") atom_type = "He";
      else if (GetAtom(i).GetSymbol() == "Ne") atom_type = "Ne";
      else if (GetAtom(i).GetSymbol() == "S") atom_type = "Ssp3";
      else if (GetAtom(i).GetSymbol() == "Ar") atom_type = "Ar";
      else if (GetAtom(i).GetSymbol() == "Kr") atom_type = "Kr";
      else {
	printf("Monomer::SetEmpiricalAtomDispersionType: Unknown Atom Type %s\n",GetAtom(i).GetSymbol().c_str());
	exit(1);
      }
      GetAtom(i).SetDispersionAtomType( atom_type );
  }
}


// Create a monomer qchem input file
void Monomer::CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job) {
  // We pass in Monomers/NMon so it can access list of other monomers.
  // This is necessary for embedding charges, for example.


  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();  
  string filename = path + "/" + mN_ + ".in"; 

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateQChemJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  // Print comment section
  fprintf(job,"$comment\nMonomer %d\n$end\n\n",Monomer_index);

  // Print $molecule section
  PrintQChemCartesian(job);

  // Print $rem section
  if (MM_job)
    fprintf(job,"%s\n",Params::Parameters().GetQChemRem2().c_str());
  else
    fprintf(job,"%s\n",Params::Parameters().GetQChemRem().c_str());
  fprintf(job,"%s\n",Params::Parameters().GetQChemBasis().c_str());

  // Optionally print $external_charges section
  if ( Params::Parameters().UseEmbeddingCharges() ) {
    fprintf(job,"$external_charges\n");
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index) {
	Monomers[i].PrintQChemEmbeddingCharges(job);
      }
    }
    fprintf(job,"$end\n");
  }
  
  
  fclose(job);
}

// Wrapper for creating MM jobs
void Monomer::CreateMMJob(Monomer Monomers[], int NMon) {
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
    CreateTinkerJob(Monomers, NMon);
  }
  else if ( Params::Parameters().GetMMType()==2 ) { //Orient
    CreateCamCaspJob();  // by Ali
  }
  else if ( Params::Parameters().GetMMType()==3 ) {// QChem job
    CreateQChemJob(Monomers, NMon, true);
  }
  else {
    printf("Monomer::CreateMMJob: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  }


}

//-----------------------------------------------------------------
// begin  by Ali
// Creates CamCasp job
void Monomer::CreateCamCaspJob() {
  
//  printf("Create CamCasp Job \n");
  
  string AIFFBasisSet = Params::Parameters().GetAIFFBasisSet();
  string CamCaspHome = Params::Parameters().GetCamCaspHome();

  // Set up the filename, with the full path.  
  string cltfile = Params::Parameters().GetMMPath() + "/" + mN_ + ".clt"; 

  /* Create the clt file */
  FILE *clt;
  if ((clt = fopen(cltfile.c_str(),"w"))==NULL) {
    printf("Monomer::CreateCamCaspJob : Cannot open file '%s'\n",
	   cltfile.c_str());
    exit(1);
  }

   fprintf(clt,"! monomer properties calculation\n"); 
   fprintf(clt,"! ==============================\n\n"); 
   fprintf(clt,"Global\n"); 
   fprintf(clt,"  Units Angstrom Degrees\n"); 
   fprintf(clt,"  Overwrite Yes\n"); 
   fprintf(clt,"  CamCASP %s \n", CamCaspHome.c_str()); 
   fprintf(clt,"End\n\n"); 
   fprintf(clt,"Molecule monomer %s\n",mN_.c_str()); 
   fprintf(clt,"  I.P. %f\n",IonizationPotential); // for Camcasp 5.6
   fprintf(clt,"CHARGE %d\n",charge);

   for (int i=0;i<Natoms;i++) {

       fprintf(clt,"%s%d      %5d     %10.6f          %10.6f           %10.6f         Type   %s%d   \n", GetAtom(i).GetSymbol().c_str(), GetAtom(i).GetAtomIndex(), GetAtom(i).GetAtomicNumber(), GetAtom(i).GetLocalPosition(0), GetAtom(i).GetLocalPosition(1), GetAtom(i).GetLocalPosition(2), GetAtom(i).GetSymbol().c_str(), GetAtom(i).GetAtomIndex());  	  
   }

   fprintf(clt,"End\n\n"); 
   fprintf(clt,"Files\n"); 
   fprintf(clt,"  Molecule monomer\n"); 
   fprintf(clt,"  Basis %s \n", AIFFBasisSet.c_str()); 
   //fprintf(clt,"  Options Tests\n"); // TURN THIS OFF... it hurts polarizabilitiess & dispersion coeffs
   //fprintf(clt,"  File-prefix monomer_%s \n",AIFFBasisSet.c_str()); 
   fprintf(clt,"  File-prefix %s\n",mN_.c_str()); 
   fprintf(clt,"  Orient file\n"); 
   fprintf(clt,"  Process file\n"); 
   fprintf(clt,"  Sites file\n"); 
   fprintf(clt,"  Interface file\n"); 
   fprintf(clt,"End\n"); 
   fprintf(clt,"Finish\n"); 
  
   fclose(clt);

}
// end
// 
// -----------------------------------------------------------------

// Creates Tinker MM job
void Monomer::CreateTinkerJob(Monomer Monomers[], int NMon) {

  // Set up the filenames, with the full path.  
  // Files are e.g. 'm1.xyz' and 'm1.key'
  string xyzfile = Params::Parameters().GetMMPath() + "/" + mN_ + ".xyz";
  string keyfile = Params::Parameters().GetMMPath() + "/" + mN_ + ".key";


  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Monomer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }

  if ( Params::Parameters().UseEmbeddingCharges() ) {
    // Count total atoms: monomer atoms + charges
    int Ntot = 0;
    for (int i=1;i<=NMon;i++) {
      Ntot += Monomers[i].GetNumberOfAtoms();
    }
    fprintf(xyz,"%d  Monomer %d with embedding charges\n",Ntot,Monomer_index);

    // Print the geometry
    int shift = 0;
    for (int i=1;i<=NMon;i++) {
      if (i==Monomer_index)
	PrintTinkerCartesian(shift,false,xyz);
      else
	Monomers[i].PrintTinkerEmbeddingCharges(shift,xyz);

      // Shift for the next monomer
      shift += Monomers[i].GetNumberOfAtoms();
    }
  }
  else
    PrintTinkerCartesian(0,true,xyz);

  fclose(xyz);
  


  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
  if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
    printf("Monomer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }
  fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
  fclose(key);
  
}

// Returns command string for running the job
string Monomer::RunQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();
  
  string infile = path + "/" + mN_ + ".in"; 
  string outfile = path + "/" + mN_ + ".out"; 


  // Execute Q-Chem
  // Q-Chem wants to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_outfile = mN_ + ".out";  
  string local_infile = mN_ + ".in";
 
  cmd += "qchem " + local_infile;
  cmd += " ";
  cmd += local_outfile;
  cmd += "; ";

  /*
  // Rename force file, if appropriate
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    string force_file = mN_ + ".force";
    cmd += "mv -f force.dat " + force_file;
    cmd += ";";
  }
  */

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// Wrapper that controls which type of MM job is executed
// Returns command string for running the job
string Monomer::RunMMJob() {
  string job;

  if (Params::Parameters().GetMMType()==1) // Tinker
    job = RunTinkerJob();
  else if (Params::Parameters().GetMMType()==2) // Orient by Ali
    job = RunCamCaspJob();  // by Ali
  else if (Params::Parameters().GetMMType()==3) // QChem
    job = RunQChemJob(true);  
  else {
    printf("ERROR: Monomer::RunMMJob(): Tinker is only MM program currently.\n");
    exit(1);
  }
  return job;
}


// GJB: Create a unix shell script to run the CamCasp job.
// We used to pass a long command string, but this ended up
// getting so long that it was easy to run into buffer overflows.
// Now we just store the name of the shell script to execute.
string Monomer::RunCamCaspJob() {

  // convert double to string
  string IP;
  std::stringstream IP_out;
  //IP_out << IonizationPotential;
  IP_out << GetIonizationPotential();
  IP = IP_out.str();

  string AIFFBasisSet = Params::Parameters().GetAIFFBasisSet();

  if (Params::Parameters().PrintLevel() > 0) {
    printf("Monomer::RunCamCaspJob(): m%d IP is %s\n", Monomer_index,IP.c_str());
    printf("Monomer::RunCamCaspJob(): m%d Basis Set is %s \n", Monomer_index, AIFFBasisSet.c_str());
  }

  // Read the path where MM files are stored
  string job_path = Params::Parameters().GetMMPath(); 

  // Set logfile name
  string logfile = mN_ + ".log";

  string scriptfilename = job_path + "/" + mN_ + ".sh";
  ofstream script;
  script.open(scriptfilename.c_str());
  script << "#!/bin/sh\n\n";

  // Set working directory
  // process id is used to give unique name to the temporary files
  char pid[10];
  sprintf(pid,"%d",Params::Parameters().GetPID());
  
  // Read scratch directory.  Requires that $QCSCRATCH is set.
  // Might want to generalize it in case we are not using Q-Chem
  string scr = getenv("QCSCRATCH");
  //printf("scr = %s\n",scr.c_str());

  string wdir = scr + "/" + mN_ + "_";
  wdir += pid;
  //printf("Working dir = %s\n",wdir.c_str());

  // First command, change to local directory
  //string cmd = "cd " + job_path;

  string cmd;

  // Script: Create a directory for the monomer and move monomer clt file
  // to this directory.  Then change to this working directory.

  script << "# Define directories" << endl;
  script << "SCRATCH_DIR=" + wdir + "\n";
  script << "JOB_DIR=" + job_path + "\n" << endl;
  
  script << "# Create the working directory and change to it" << endl;
  cmd = "mkdir $SCRATCH_DIR\n";
  cmd += "cp " + mN_ + ".clt $SCRATCH_DIR\n";
  cmd += "cd $SCRATCH_DIR\n";

  script << cmd << endl;
  // Script: Now start running the job
  script << "# Compute distributed multipole moments and polarizabilities" << endl;

  cmd = "cluster < " + mN_ + ".clt > " + "$JOB_DIR/" + logfile + "\n";
  // Eliminate IP flag for for CamCasp 5.6
  //cmd += "runprops " + mN_ + " -clt " + mN_ + ".clt -ip " + IP + " -d " 
  //  + mN_ + "_" + pid + " -q bg >> " + "$JOB_DIR/" + logfile + "\n"; 
  cmd += "runprops " + mN_ + " -clt " + mN_ + ".clt -d " 
    + mN_ + "_" + pid + " -q bg >> " + "$JOB_DIR/" + logfile + "\n"; 
  cmd += "cd " + mN_ + "_" + pid + "\n";
  cmd += "localize " + mN_ + " --pol " + mN_ + "_NL4_ >> $JOB_DIR/" 
    + logfile + " 2>&1 \n";
  // note above: the "2>&1" bit redirects stderr as well as stdout.  In bash, >>& doesn't
  // work.  So this equivalent has to be used instead.
  script << cmd << endl;

  script << "# Save the multipole moments and polarizabilities" << endl;
  // Store the multipole moments and polarizabilities in a simple file
  cmd = "cp OUT/" + mN_ + "_DMA2_L4.mom $JOB_DIR/" + mN_ + ".mom\n";
  cmd += "cp " + mN_ + "_ref_wt4_L2_static.pol $JOB_DIR/" + mN_ + ".pol\n";
  script << cmd;

  if ( Params::Parameters().UseDiagonalAnisotropicPolarizabilities() ) {
    // Also grab  the (uncoupled) frequency-dependent polarizabilities, used to
    // compute the dispersion coefficients.
    cmd = "cp " + mN_ + "_ref_wt4_L2.casimir $JOB_DIR/" + mN_ + ".freq_pol\n";
    script << cmd;
  }
  else {
    script << endl;
    script << "# Compute isotropic polarizabilities and dispersion coefficients" << endl;
    // Now compute the isotropic dispersion coefficients.  Need to re-localize
    // the polarizabilities with flag --isotropic.
    cmd = "localize " + mN_ + " --pol " + mN_ + "_NL4_ --isotropic >> $JOB_DIR/" 
      + logfile + " 2>&1\n";
    // The next two files are actually not important any more, since we now 
    // read the frequency-dependent polarizabilities and compute the dispersion
    // coefficients ourselves.
    cmd += "cp " + mN_ + "_wt4_L2iso_C10iso.pot $JOB_DIR/" + mN_ + ".disp\n";
    cmd += "cp " + mN_ + "_ref_wt4_L2_static.pol $JOB_DIR/" + mN_ + ".iso_pol\n";
    // Also grab  the (uncoupled) frequency-dependent polarizabilities, used to
    // compute the dispersion coefficients.
    cmd += "cp " + mN_ + "_ref_wt4_L2.casimir $JOB_DIR/" + mN_ + ".freq_pol\n";
    script << cmd;
  }

  script << endl << "# Clean up" << endl;
  // Switch back to base directory
  cmd = "cd " + Params::Parameters().GetBasePath();

  // remove the temporary orient files
  if (! Params::Parameters().OrientDebug() ) {
    cmd += "\nrm -rf $SCRATCH_DIR\n";
  }
  script << cmd << endl;

  // Close the script file
  script.close();

  // Create the cmd to be executed by our code
  cmd = "cd " + job_path + "; sh " + scriptfilename + "; cd " 
    + Params::Parameters().GetBasePath();

  return cmd;
}

// Returns command string for running the job
string Monomer::RunTinkerJob() {
  // Set up the filenames, with no path info
  // Files are e.g. 'm1.xyz' and 'm1.out'
  string infile = mN_ + ".xyz"; 
  string outfile = mN_ + ".out"; 


  // Execute Tinker
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath(); 
  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "analyze " + infile;
  cmd += " e > ";
  cmd += outfile + "; ";

  // If doing force job, compute the gradient
  if ( Params::Parameters().DoForces()  || Params::Parameters().DoFiniteDifferenceFreqs()) {
    // Run the job
    cmd += "minimize " + infile;
    cmd += " 100000000 > ";
    string force_file = mN_ + ".force";
    cmd += force_file + ";";

    // Remove tmp file and extra geom file created by minimize job
    cmd += "rm -f " + infile;
    cmd += "_2;";
  }

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  // Actual job running & checking of status have been moved to main.C 
  // to simplify parallel implementation.

  return cmd;
}

void Monomer::ReadQMResults() {

  // Read in the Energy
  Energy_QM = ReadQChemEnergy();

  double MP2_dispersion_correction;
  if ( Params::Parameters().DoMP2DispersionCorrection() )
    MP2_dispersion_correction = ComputeInteratomicMP2TwoBodyDispersionCorrection();

  // Optionally Read in the Forces
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetQMGradient();
  }
}

void Monomer::ReadMMResults() {

  // If we're using the AIFF, we read in the monomer force field
  // parameters here.  SetMMEnergy either sets the monomer energy zero
  // (AIFF) or reads from an output file (other types).

  if ( Params::Parameters().GetMMType()==2 ) {// AIFF
    // Read in the ab initio force field parameters
    ReadMultipoleMoments();
    ReadPolarizabilities();
    // Next two have been replaced with improved AIFF routines
    //ReadIsotropicDipolePolarizability();
    //ReadDispersionCoefficients();

    if ( Params::Parameters().UseDiagonalAnisotropicPolarizabilities() )
      ReadDiagonalAnisotropicFreqPolarizabilities();
    else 
      ReadFreqPolarizabilities();// by shuhao
  
    SetEmpiricalAtomDispersionType(); // used to figure out dispersion damping
  }

  // Read in the Energy
  SetMMEnergy();

  // Optionally Read in the Forces
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient();
  }
}

// Get the QM Gradient, wrapper routine
void Monomer::SetQMGradient() {
  string path = Params::Parameters().GetQMPath();
  Grad_QM = ReadGradient(path,1);
  QM_Grad_Init = 1;
}

// Get the MM Gradient, wrapper routine
void Monomer::SetMMGradient() {
  if (Params::Parameters().GetMMType() == 1) { // Tinker 
    string path = Params::Parameters().GetMMPath();
    Grad_MM = ReadGradient(path,2);
    MM_Grad_Init = 1;
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    // AIFF gradient handled elsewhere
    if (!MM_Grad_Init) Grad_MM.Initialize(3*Natoms); // create empty array
    MM_Grad_Init = 1;
  }
  else if (Params::Parameters().GetMMType() == 3) { // QChem
    string path = Params::Parameters().GetMMPath();
    Grad_MM = ReadGradient(path,1);
    MM_Grad_Init = 1;
  }
  else {
    printf("Monomer::SetMMGradient() - MM type = %d is unknown\n",
	   Params::Parameters().GetMMType());
    exit(1);
  }

}

// Get the QM Gradient, copied from another Monomer
void Monomer::SetQMGradient(Monomer& other) {
  if (other.QM_Grad_Init) {
    Grad_QM = other.GetQMGradient();
    QM_Grad_Init = 1;
  }
}

// Get the MM Gradient, copied from another Monomer
void Monomer::SetMMGradient(Monomer& other) {
  if (other.MM_Grad_Init) {
    Grad_MM = other.GetMMGradient();
    MM_Grad_Init = 1;
  }
}

// Main routine for reading in gradient files
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
Vector Monomer::ReadGradient(string path, int type) {

  // Initialize the gradient
  Vector grad(3*Natoms);

  // Set up the filename, with the full path.  File is e.g. 'm1.force'
  string filename = path + "/" + mN_;
  if (type == 2) // Tinker MM job
    filename += ".force";
  else // other
    filename += ".out";
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadGradient : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  // Read in the data - search for the "Nuclear forces:" string
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    // Search for final q-chem energy
    if ( line==" Nuclear forces:" ) {
      getline(infile,line); // throw away header line

      for (int i=0;i<Natoms;i++) {
	getline(infile,line);
	istringstream iss(line);
	string tmp;
	iss >> tmp; // throw away the atom index
	for (int j=0;j<3;j++) {
	  iss >> grad[3*i+j]; // Store the gradient elements
	}
      }
      break;
    }
  }

  infile.close();  // Close the file

  //PrintGradient(" MonomerGradient",grad);
  return grad;
}

// Find the minimum distance between this and another monomer Mon
Vector Monomer::FindDistance(const Monomer &MonB) {

  //printf("Finding distance\n"); fflush(stdout);

  // Choose 1 of the following options for defining the distance 
  bool use_min_dist = true; // separation between 2 closest atoms
  bool use_com_dist = false; // separation between centers of mass

  double min_dist = 9999999999.0;

  int atomA=0, atomB=0;

  Vector Result(3); Result.Set();

  // Minimum Distance criterion
  if (use_min_dist) {
    double dist;
    int NatomsB = MonB.Natoms;//GetNumberOfAtoms();      
    for (int i=0;i<Natoms;i++) {
      for (int j=0;j<NatomsB;j++) {             
        dist = Atom_List[i].GetInterAtomicDistance(MonB.Atom_List[j]);
        //dist = Atom_List[i].GetInterAtomicDistance(MonB.GetAtom(j));
        //printf("dist = %f   min_dist = %f\n",dist,min_dist);
        if (dist < min_dist) {
          min_dist = dist;
	  atomA = i;
	  atomB = j;
          //printf("min_dist updated to %f for atom pair %d,%d\n",min_dist,i,j);
	}
      }
    }
    //printf("Final min_dist = %f, between %d and %d\n",min_dist,atomA,atomB);
    Result[0] = min_dist;
    Result[1] = (double) atomA;
    Result[2] = (double) atomB;

  }

  // Center-of-Mass Distance criterion
  else if (use_com_dist) {
    Vector Acom = GetCenterOfMass();
    Vector Bcom = MonB.CenterOfMass; //GetCenterOfMass();
    Acom -= Bcom;
    min_dist = Acom.Norm();
    //printf("Dist = %f\n",min_dist);
    Result[0] = min_dist;
  }

  return Result;
}

// Print out the Cartesian coordinates in xyz format
void Monomer::PrintMonomerCartesian(FILE *outfile) {
  for (int i=0;i<Natoms;i++) {
    GetAtom(i).PrintQChemCartesian(outfile);
  }
}

void Monomer::PrintGhostMonomerCartesian(FILE *outfile) {
  for (int i=0;i<Natoms;i++) {
    fprintf(outfile, "@%s  %10.6f  %10.6f  %10.6f \n", 
      GetAtom(i).GetSymbol().c_str(), GetAtom(i).GetCoordinate(0),
      GetAtom(i).GetCoordinate(1), GetAtom(i).GetCoordinate(2));
  }
}

// Print out the Cartesian coordinates in xyz format
void Monomer::PrintMonomerCartesian(FILE *outfile, string symbol) {
  for (int i=0;i<Natoms;i++) {
    GetAtom(i).PrintQChemCartesian(outfile,symbol);
  }
}

// Print out $molecule section for Q-Chem
void Monomer::PrintQChemCartesian(FILE *outfile) {
  fprintf(outfile,"$molecule\n%d %d\n",charge,spin);
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintQChemCartesian(outfile);
  }
  fprintf(outfile,"$end\n");
  
}

// Print out the Cartesian coordinates in Tinker xyz format
void Monomer::PrintTinkerCartesian(int shift, bool monomer, FILE *outfile) {
  
  if (monomer)
    fprintf(outfile,"%d  (Monomer %d)\n",Natoms, Monomer_index); 

  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintTinkerCartesian(shift,outfile);
  }
}


// Print out the Tinker-style coordinates for embedding charges
void Monomer::PrintTinkerEmbeddingCharges(int shift, FILE *outfile) {
  
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintTinkerEmbeddingCharge(shift,outfile);
  }
}


// Wrapper for printing QM gradient
void Monomer::PrintQMGradient(string title) {
  //printf("QM_Grad_Init = %d\n",QM_Grad_Init);  
  if (QM_Grad_Init) 
    PrintGradient(title,Grad_QM);
  else 
    printf("QM Gradient not initialized\n");
}

// Wrapper for printing <M gradient
void Monomer::PrintMMGradient(string title) {
  printf("MM_Grad_Init = %d\n",MM_Grad_Init);
  if (MM_Grad_Init) 
    PrintGradient(title,Grad_MM);
  else 
    printf("MM Gradient not initialized\n");
}

// Print out Gradient as an Nx3 matrix
void Monomer::PrintGradient(string title, Vector grad) {

  printf("%s\n",title.c_str());
  for (int i=0;i<Natoms;i++) {
    printf("%2s %15.10f %15.10f %15.10f\n",Atom_List[i].GetSymbol().c_str(),
    	   grad[3*i],grad[3*i+1],grad[3*i+2]);
  }
}

void Monomer::PrintAll() {

  printf("-----------------------------------------------\n");
  printf("Monomer %d           %d atoms\n\n",Monomer_index,Natoms);
  PrintQChemCartesian();

  printf("\nEmbedding charges:\n");
  PrintQChemEmbeddingCharges();

  printf("Energies:\nQM: %15.9f\nMM: %15.9f\n\n",Energy_QM,Energy_MM);

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    printf("Forces:\n");
    PrintQMGradient("QM:\n");
    PrintMMGradient("MM:\n");
  }
  printf("-----------------------------------------------\n");


}

void Monomer::PrintQChemEmbeddingCharges(FILE* outfile) {

  //fprintf(outfile,"$external_charges\n");
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintQChemEmbeddingCharge(outfile);
  }
  //fprintf(outfile,"$end\n\n");
  

}

//Translate molecule such that it has new_com center of mass.
void Monomer::Translate(Vector new_com) {

  // find coordinates of each atom relative to the center of mass (com),
  // shift the com, and adjust the atom coordinates to the new com.
  double old_coord[3]; double new_coord[3];
  
  for (int i=0;i<Natoms;i++) {
    for (int dim=0;dim<3;dim++) {
      old_coord[dim] = GetAtom(i).GetCoordinate(dim);
      new_coord[dim] = old_coord[dim] - GetCenterOfMass(dim) + new_com[dim];
    }    
    Atom_List[i].SetPosition(new_coord);
  }
  
  // Update the center of mass
  FindCenterOfMass();
   
}


void Monomer::ComputeIntramolecularDistances() {

  printf(" Atom 1   Atom 2       Distance (Ang)\n");
  for (int i=0;i<Natoms;i++) {
    string atom1 = Atom_List[i].GetSymbol();
    for (int j=i+1;j<Natoms;j++) {		
      string atom2 = Atom_List[j].GetSymbol();
      
      double dist = Atom_List[i].GetInterAtomicDistance( Atom_List[j] );
      printf("%2s(%3d)  %2s(%3d)       %7.4f\n",atom1.c_str(),i,
	     atom2.c_str(),j,dist);
    }
  }
  printf("\n");
}


// Overload "=" operator
// Assumes none of the arrays have been initialized in the copy
Monomer& Monomer::operator=(const Monomer& other) {

  if (this!=&other) {
    Monomer_index = other.Monomer_index;
    ref_mon = other.ref_mon;
    mN_ = other.mN_;
    spin = other.spin;
    charge = other.charge;
    Natoms = other.Natoms;

    // Copy over Atom_List
    // First delete old list, just in case. It should either be
    // allocated or NULL, so no harm in using delete
    delete [] Atom_List; 
    Atom_List = new Atom[Natoms];
    for (int i=0;i<Natoms;i++)
      Atom_List[i] = other.Atom_List[i];

    IonizationPotential = other.IonizationPotential;
    
    MonomerMass = other.MonomerMass;
    CenterOfMass = other.CenterOfMass;

    QM_Job_Complete = other.QM_Job_Complete;
    MM_Job_Complete = other.MM_Job_Complete;
    
    Energy_QM = other.Energy_QM;
    Energy_MM = other.Energy_MM;
    
    QM_Grad_Init = other.QM_Grad_Init;
    if (QM_Grad_Init) {
      Grad_QM = other.Grad_QM;
    }

    MM_Grad_Init = other.MM_Grad_Init;
    if (MM_Grad_Init) {
      Grad_MM = other.Grad_MM;
    }

    RotationAngle = other.RotationAngle;
    RotationVector = other.RotationVector;
  }
  return *this;
}


double Monomer::ComputeInteratomicMP2TwoBodyDispersionCorrection() {

  double E2b_disp = 0.0;

  for (int i=0;i<Natoms;i++) {
    Atom AtomI = Atom_List[i];

    // Get the pure C6 coefficients for atom I
    double C6i_uchf = AtomI.GetC6Coefficient("UCHF");
    double C6i_new = AtomI.GetC6Coefficient();

    for (int j=i+1;j<Natoms;j++) {
      Atom AtomJ = Atom_List[j];

      // Get the pure C6 coefficients for atom J
      double C6j_uchf = AtomJ.GetC6Coefficient("UCHF");
      double C6j_new = AtomJ.GetC6Coefficient();

      // Compute mixed C6 coefficients for this atom pair
      // C6ij = sqrt(C6i*C6j)
      double C6_uchf = sqrt(C6i_uchf*C6j_uchf);
      double C6_new = sqrt(C6i_new*C6j_new);

      // Get the interatomic distance
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Compute the damping function - requires vdW diameters
      double F6ij = 1.0;
      if (Params::Parameters().DoDispersionDamping()) {

	// Get the van der Waals diameters
	double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");

	bool DoTangToenniesDamping = true;
	if (DoTangToenniesDamping) {
	  // Get Tang-Toennies-type 3-body damping factor
	  // Use empirical expression for beta: J Chem Phys 131, 094106 (2009).
	  double betaIJ = -0.334*(Di+Dj) + 4.386;
	  F6ij = AtomI.TangToenniesDampingFactor(6,betaIJ,Rij);
	}
	else {// Aziz damping
	  double S = 1.4;
	  F6ij = exp(-(S*(Di+Dj)/Rij-1)*(S*(Di+Dj)/Rij-1));
	}
      }

      //printf("(%d,%d): C6(uchf) = %f, C6(cks) = %f, Rij = %f.  F6ij = %f\n",i,j,C6_uchf,C6_new,Rij, F6ij);
      // Compute the dispersion contribution
      E2b_disp += -F6ij*(C6_new - C6_uchf)/pow(Rij,6);

    }
  }
  printf("Monomer %d dispersion correction: %f\n",Monomer_index,E2b_disp);
  return E2b_disp;
}

// Obsolete: An empirical estimate for interatomic 3-body dispersion.
// Use Monomer::ComputeThreeBodyDispersion() instead.
double Monomer::ComputeInteratomicThreeBodyDispersion() {
  double E3b_disp = 0.0;

  for (int i=0;i<Natoms;i++) {
    Atom AtomI = Atom_List[i];
    for (int j=i+1;j<Natoms;j++) {
      Atom AtomJ = Atom_List[j];
      for (int k=j+1;k<Natoms;k++) {
	Atom AtomK = Atom_List[k];
	
	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;
	
	// use law of cosines to get cos(phi)
	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;

	// Get the C9 coefficient
	double C9ijk = AtomI.EstimateC9Coefficient(AtomJ, AtomK, "Tkatchenko");
	
	// Compute the short-range damping function
	double damping = 1.0; // no damping
	if (Params::Parameters().DoDispersionDamping()) {

	  // Get the van der Waals diameters
	  double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	  double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	  double Dk = AtomK.LookupAtomicDispersionParameter("Rvdw");
	  
	  // Get Tang-Toennies-type 3-body damping factor
	  // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	  double betaIJ = -0.31*(Di+Dj) + 3.43;
	  double betaIK = -0.31*(Di+Dk) + 3.43;
	  double betaJK = -0.31*(Dj+Dk) + 3.43;
	  
	  double F6ij = AtomI.TangToenniesDampingFactor(6,betaIJ,Rij);
	  double F6ik = AtomI.TangToenniesDampingFactor(6,betaIK,Rik);
	  double F6jk = AtomJ.TangToenniesDampingFactor(6,betaJK,Rjk);
	  double damping = F6ij * F6ik * F6jk;
	  //printf("M: F6ij = %f, F6ik = %f, F6jk = %f, damping =%f\n",F6ij, F6ik, F6jk,damping);
	}

	// ATM 3-body dispersion contribution
	//printf("M(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
      }
    }
  }
  
  //printf("Monomer %d: Interatomic 3-body dispersion: %f in kcal/mol\n",
  // Monomer_index,E3b_disp*HartreesToKcalpermole);
  return E3b_disp;
}

// AIFF function for computing 3-body Axilrod-Teller dispersion.  The
// C9 coefficient comes from Casimir-Polder integration.
double Monomer::ComputeThreeBodyDispersion() {
  double E3b_disp = 0.0;

  for (int i=0;i<Natoms;i++) {
    Atom AtomI = Atom_List[i];
    for (int j=i+1;j<Natoms;j++) {
      Atom AtomJ = Atom_List[j];
      for (int k=j+1;k<Natoms;k++) {
	Atom AtomK = Atom_List[k];
	
	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;
	
	// use law of cosines to get cos(phi)
	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;
	
	// Get the C9 coefficient
	double C9ijk = AtomI.CasimirC9Coefficient(AtomJ, AtomK);

	// Compute the short-range damping function
	double damping = 1.0; // no damping
	if (Params::Parameters().DoDispersionDamping()) {

	  // Get the van der Waals diameters
	  double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
	  double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
	  double Dk = AtomK.LookupAtomicDispersionParameter("Rvdw");
	  
	  // Get Tang-Toennies-type 3-body damping factor
	  // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
	  double betaIJ = -0.31*(Di+Dj) + 3.43;
	  double betaIK = -0.31*(Di+Dk) + 3.43;
	  double betaJK = -0.31*(Dj+Dk) + 3.43;
	  
	  double F6ij = AtomI.TangToenniesDampingFactor(6,betaIJ,Rij);
	  double F6ik = AtomI.TangToenniesDampingFactor(6,betaIK,Rik);
	  double F6jk = AtomJ.TangToenniesDampingFactor(6,betaJK,Rjk);
	  double damping = F6ij * F6ik * F6jk;
	}

	// ATM 3-body dispersion contribution
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
      }
    }
  }
  
  return E3b_disp;
}
