#include <string>
using std::string;
#include "monomer.h"
#include "quasiharmonic.h"
#include <sstream>  // by Ali
#include "constants.h"
#include <sys/stat.h>
using namespace hmbi_constants;

Monomer::Monomer() : Monomer_index(0), ref_mon(0), mN_("m0"), spin(1),sym_fac(1), charge(0), 
		     Natoms(0), Atom_List(NULL), Sym_Atom(NULL), IonizationPotential(0.0),
		     MonomerMass(0.0), QM_Job_Complete(false),
		     MM_Job_Complete(false), Energy_QM(0.0), Energy_MM(0.0),
		     QM_Grad_Init(0), MM_Grad_Init(0), QM_Hess_Init(0), MM_Hess_Init(0),
		     RotationAngle(0.0){
  // The constructor contents got moved to the Initialize member function.
  // empty constructor, just initializes pointers to zero.

  CenterOfMass.Initialize(3);
  RotationVector.Initialize(3);
  Symmetry_Rotation.Initialize(3,3);
  Symmetry_Rotation.Set_Iden();
  
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    Fractional_Rotation.Initialize(3,3);
    Fractional_Rotation.Set_Iden();
    Fractional_Translate.Initialize(3);
  }
  
  //Local_Rotation.Initialize(3,3);
  //Local_Rotation.Set_Iden();
  //Opt_Center.Initialize(3);

}

Monomer::Monomer(const Monomer &other) : Monomer_index(other.Monomer_index),
					 ref_mon(other.ref_mon),
					 mN_(other.mN_), 
					 spin(other.spin),
					 charge(other.charge),
 					 Monomer_type(other.Monomer_type), 
					 Natoms(other.Natoms), 
                                         unique_atoms(other.unique_atoms),
					 Atom_List(NULL), 
                                         Sym_Atom(NULL),
					 sym_fac(other.sym_fac),
					 SymMon(other.SymMon),
					 IonizationPotential(other.IonizationPotential),
					 MonomerMass(other.MonomerMass),
					 QM_Job_Complete(other.QM_Job_Complete),
					 MM_Job_Complete(other.MM_Job_Complete),
					 Energy_QM(other.Energy_QM),
					 Energy_MM(other.Energy_MM),
					 QM_Grad_Init(other.QM_Grad_Init),
					 MM_Grad_Init(other.MM_Grad_Init),
                                         QM_Hess_Init(other.QM_Hess_Init),
                                         MM_Hess_Init(other.MM_Hess_Init),
					 RotationAngle(other.RotationAngle)
					{
  CenterOfMass = other.CenterOfMass;
  //FindCenterOfMass();
  // Copy over Atom_List and symmetry list
  if (Natoms > 0) {
    Atom_List = new Atom[Natoms];
    Sym_Atom = new int[Natoms];
    for (int i=0;i<Natoms;i++){
      Atom_List[i] = other.Atom_List[i];
      Sym_Atom[i] = other.Sym_Atom[i];
    }
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
  
  // Copy over the Hessians
  QM_Hess_Init = other.QM_Hess_Init;
  if (QM_Hess_Init) {
    Hess_QM = other.Hess_QM;
  }
  MM_Hess_Init = other.MM_Hess_Init;
  if (MM_Hess_Init) {
    Hess_MM = other.Hess_MM;
  }

  RotationVector = other.RotationVector;
  Symmetry_Rotation = other.Symmetry_Rotation;
  //Local_Rotation = other.Local_Rotation;
  //Opt_Center = other.Opt_Center;
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    Fractional_Rotation = other.Fractional_Rotation;
    Fractional_Translate = other.Fractional_Translate;
  }
}



Monomer::~Monomer() {
  delete [] Atom_List;
  delete [] Sym_Atom;
}

// Initialization, w/ no feature for MM details
void Monomer::Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, string atoms[], 
			 double *xyz, int natoms) {
  Monomer_index = ind;
  ref_mon = ind; // by default, it is not a PBC image or a supercell monomer, so ref_mon is just itself. 
  SetLabel(); 
  charge = charge_in;
  spin = spin_in;
  Monomer_type = type;
  Natoms = natoms;
  unique_atoms = natoms;
  MonomerMass = 0.0;
  SymMon = ind;
  delete [] Sym_Atom;
  Sym_Atom = new int[Natoms];

  // Read in the geometry
  delete [] Atom_List;
  Atom_List = new Atom[Natoms];

  //Atom *TempAtoms = new Atom[Natoms];
  //Atom_List = TempAtoms;
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].Initialize(i+1,natoms_previous,atoms[i],xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
    MonomerMass += Atom_List[i].GetAtomicMass();
    Sym_Atom[i] = i;
  }

  FindCenterOfMass();

  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali

  //Initialize Gradient
  QM_Grad_Init = 0;
  MM_Grad_Init = 0;
  if(Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().DoFreq()){
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }

  //Initialize Hessian
  MM_Hess_Init = 0;               
  QM_Hess_Init = 0;
  if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
    Hess_QM.Initialize(3*Natoms, 3*Natoms);
    Hess_MM.Initialize(3*Natoms, 3*Natoms);
  }

  na_ = 0;
  nb_ = 0;
  nc_ = 0;
  
}

// Initialization, w/ no feature for MM details
void Monomer::Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, string atoms[], 
			 double *xyz, int natoms, Vector charges) {
  Monomer_index = ind;
  ref_mon = ind; // by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel(); 
  charge = charge_in;
  spin = spin_in;
  Monomer_type = type;
  Natoms = natoms;
  unique_atoms = natoms;
  MonomerMass = 0.0;

  //symmetrical monomer itself by default
  SymMon = ind;
  Sym_Atom = new int[Natoms];

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].Initialize(i+1,natoms_previous,atoms[i],xyz[3*i], xyz[3*i+1], xyz[3*i+2], charges[i]);
    MonomerMass += Atom_List[i].GetAtomicMass();
    Sym_Atom[i] = i;
  }
  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali

  //Initialize Gradient
  if(Params::Parameters().DoForces() || Params::Parameters().DoFreq()){
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }

  na_ = 0;
  nb_ = 0;
  nc_ = 0;

}


// Initialization for QM and Tinker atoms
void Monomer::Initialize(int ind,int natoms_previous, int charge_in, int spin_in, string type, 
			 string atoms[], double *xyz, int natoms, int *atom_types, int *Nconnected, 
			 int* connectivity)  {
  Monomer_index = ind;
  ref_mon = ind;// by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel();
  charge = charge_in;
  spin = spin_in;
  Monomer_type = type;
  Natoms = natoms;
  unique_atoms = natoms;
  MonomerMass = 0.0;

  SymMon = ind;
  Sym_Atom = new int[Natoms];

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  double coord[3];
  for (int i=0;i<Natoms;i++) {
    coord[0] = xyz[3*i];
    coord[1] = xyz[3*i+1];
    coord[2] = xyz[3*i+2];
    Atom_List[i].Initialize(i+1,natoms_previous, atoms[i], coord, atom_types[i], 
			    Nconnected[i], &connectivity[6*i]);
    MonomerMass += Atom_List[i].GetAtomicMass();
    Sym_Atom[i] = i;
  }

  FindCenterOfMass();
  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali

  //Initialize Gradient
  if(Params::Parameters().DoForces() || Params::Parameters().DoFreq()){
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }
  
}

// Initialization for QM and Tinker atoms, with embedding point charges
void Monomer::Initialize(int ind,int natoms_previous,int charge_in, int spin_in,
			 string type, string atoms[], double *xyz, int natoms,
			 int *atom_types, int *Nconnected, int* connectivity, Vector charges) {
  Monomer_index = ind;
  ref_mon = ind;// by default, it is not a PBC image, so ref_mon is just itself. 
  SetLabel();
  charge = charge_in;
  spin = spin_in;
  Monomer_type = type;
  Natoms = natoms;
  unique_atoms = natoms;
  MonomerMass = 0.0;
  SymMon = ind;
  Sym_Atom = new int[Natoms];

  // Read in the geometry
  Atom_List = new Atom[Natoms];
  double coord[3];
  for (int i=0;i<Natoms;i++) {
    coord[0] = xyz[3*i];
    coord[1] = xyz[3*i+1];
    coord[2] = xyz[3*i+2];



    Atom_List[i].Initialize(i+1, natoms_previous, atoms[i], coord, atom_types[i], 
			    Nconnected[i], &connectivity[6*i], charges[i]);

    
    MonomerMass += Atom_List[i].GetAtomicMass();
    Sym_Atom[i] = i;
  }

  FindCenterOfMass();

  if (Params::Parameters().GetMMType()==2)
    FindLocalCoord(); // by Ali

  //Initialize Gradient
  if(Params::Parameters().DoForces() || Params::Parameters().DoFreq()){
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }

  na_ = 0;
  nb_ = 0;
  nc_ = 0;

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
  CenterOfMass.Set(); //in case value is already set
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
  //if (Params::Parameters().BuildForceFieldOnly())
  //  Use_Ali_Scheme = true; 

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
    //z_axis[0] = GetAtom(1).GetCoordinate(0) - GetAtom(0).GetCoordinate(0);
    //z_axis[1] = GetAtom(1).GetCoordinate(1) - GetAtom(0).GetCoordinate(1);
    //z_axis[2] = GetAtom(1).GetCoordinate(2) - GetAtom(0).GetCoordinate(2);
    //z_axis.Normalize();
    // Forming normalized x-axis using Gram-Schmidt
      // pos_vec = position vector of the 0th atom --> not essential anymore
      //pos_vec[0] = GetAtom(0).GetCoordinate(0);
    //pos_vec[1] = GetAtom(0).GetCoordinate(1);
    //pos_vec[2] = GetAtom(0).GetCoordinate(2);
   
    //Forming normalized x-axis using Gram-Schmidt if the 1st atom is
    //origin or if the origin and the two atoms fall on the same line
    //(crossproduct =0), generate and use a random vector

    //if ( (fabs(pos_vec.CrossProduct(z_axis)[0])<0.0001) && (fabs(pos_vec.CrossProduct(z_axis)[1])<0.0001) && (fabs(pos_vec.CrossProduct(z_axis)[2])<0.0001) ) {
    //pos_vec[0] = rand();
    //pos_vec[1] = rand(); 
    //pos_vec[2] = rand(); 
    //pos_vec.Normalize();
      //pos_vec.Print("pos_vec_random");

    //random vec should be printed out for reference??
    //  pos_vec.Print("random position vector used");
    // } 
    //Commented this out because using cross products are a better way to determine in the first 
    //atoms are on the same line instand of dot products
    //Forming normalized x-axis using Gram-Schmidt if the 1st atom is
    //origin or if the origin and the two atoms fall on the same line
    //(dotproduct =0), generate and use a random vector
  /* if ( pos_vec.DotProduct(z_axis) == 0.00 ) {
      pos_vec[0] = rand();
      pos_vec[1] = rand();
      pos_vec[2] = rand();
      pos_vec.Normalize();


      //pos_vec.Print("pos_vec_random");

      //random vec should be printed out for reference??
      //  pos_vec.Print("random position vector used");
   }*/

    // using the projection of pos_vec onto the z-axis to find x-axis
    /*
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
    */
    //printf("Local coordinate axes:\n");
    //printf("X: %10.6f %10.6f %10.6f\n",x_axis[0],x_axis[1],x_axis[2]);
    //printf("Y: %10.6f %10.6f %10.6f\n",y_axis[0],y_axis[1],y_axis[2]);
    //printf("Z: %10.6f %10.6f %10.6f\n",z_axis[0],z_axis[1],z_axis[2]);
    

    
    //End : Gram-Schmidt


    // To rotate the global coord sys to the local one (for Orient program), 
    // we need a rotation angle,  
    // and a vector about which the rotation takes place.
    
    /* Matrix Inertia(3,3);
    Inertia.Set();
    printf("COM coordinates\n");
    //    printf("Center of Charge is:%f\n",CenterOfCharge[0]);
    //    printf("Center of Charge is:%f\n",CenterOfCharge[1]);
    //    printf("Center of Charge is:%f\n",CenterOfCharge[2]);
    //shift origin to the center of mass
      for (int iatom=0;iatom<Natoms;iatom++) {
      double x = GetAtom(iatom).GetCoordinate(0) - CenterOfMass[0]; 
      double y = GetAtom(iatom).GetCoordinate(1) - CenterOfMass[1];
      double z = GetAtom(iatom).GetCoordinate(2) - CenterOfMass[2];
      double m = GetAtom(iatom).GetAtomicMass();
      //double m = GetAtom(iatom).GetAtomicNumber();
      //    printf("Atomic number is: %f\n",m);
      printf("%f %f %f\n",x,y,z);
      //printf("x is %f\n",x);
      //printf("y is %f\n",y);
      //printf("z is %f\n",z);
      //printf("GetAtom(iatom).GetGoordinate is %f\n",GetAtom(iatom).GetCoordinate())
      //construct the inertia tensor matrix
      Inertia(0,0) += m*(y*y + z*z);
      Inertia(1,1) += m*(x*x + z*z);
      Inertia(2,2) += m*(x*x + y*y);
      
      Inertia(0,1) -= m*x*y;
      Inertia(0,2) -= m*x*z;
      Inertia(1,2) -= m*y*z;
    }
    Inertia(1,0) = Inertia(0,1);
    Inertia(2,0) = Inertia(0,2);
    Inertia(2,1) = Inertia(1,2);

    Inertia.Print("Inertia Tensor");

    //Diagonalize the matrix and obtain the eigenvectors, take the eigenvectors as the x,y,z local coordinates
    Vector Imom(3);
    Imom = Inertia.Diagonalize();
    Inertia.Print("Inertia vectors");
    Imom.Print("Imom is:");
    x_axis = Inertia.GetColumnVector(0);
    y_axis = Inertia.GetColumnVector(1);
    z_axis = Inertia.GetColumnVector(2);

    //For the zero vector case, such as linear molecule
    if (x_axis[0]<0.000001 && x_axis[1]<0.000001 && x_axis[2]<0.000001) {
      x_axis[0] = y_axis[1]*z_axis[2] - y_axis[2]*z_axis[1];
      x_axis[1] = y_axis[2]*z_axis[0] - y_axis[0]*z_axis[2];
      x_axis[2] = y_axis[0]*z_axis[1] - y_axis[1]*z_axis[0];
    }

    if (y_axis[0]<0.000001 && y_axis[1]<0.000001 && y_axis[2]<0.000001) {
      y_axis[0] = z_axis[1]*x_axis[2] - z_axis[2]*x_axis[1];
      y_axis[1] = z_axis[2]*x_axis[0] - z_axis[0]*x_axis[2];
      y_axis[2] = z_axis[0]*x_axis[1] - z_axis[1]*x_axis[0];
    }

    if (z_axis[0]<0.000001 && z_axis[1]<0.000001 && z_axis[2]<0.000001) {
      z_axis[0] = x_axis[1]*y_axis[2] - x_axis[2]*y_axis[1];
      z_axis[1] = x_axis[2]*y_axis[0] - x_axis[0]*y_axis[2];
      z_axis[2] = x_axis[0]*y_axis[1] - x_axis[1]*y_axis[0];
    }


    x_axis.Normalize();
    y_axis.Normalize();
    z_axis.Normalize();

    Vector temp(3);
    //For the case De and Rotation Vector equal 0

    //Y,X,Z
        if ((y_axis[2]-z_axis[1]<0.000001)&&(z_axis[0]-x_axis[2]<0.000001)&&(x_axis[1]-y_axis[0]<0.000001)) {
      temp = x_axis;
      x_axis = y_axis;
      y_axis = temp;
    }
    //X,Z,Y
    if ((y_axis[2]-z_axis[1]<0.000001)&&(z_axis[0]-x_axis[2]<0.000001)&&(x_axis[1]-y_axis[0]<0.000001)) {
      temp = x_axis;
      x_axis = y_axis;
      y_axis = z_axis;
      z_axis = temp;
    }
    //Z,Y,X
    if ((y_axis[2]-z_axis[1]<0.000001)&&(z_axis[0]-x_axis[2]<0.000001)&&(x_axis[1]-y_axis[0]<0.000001)) {
      temp = x_axis;
      x_axis = y_axis;
      y_axis = z_axis;
      z_axis = temp;
    }
    //Y,Z,X
    if ((y_axis[2]-z_axis[1]<0.000001)&&(z_axis[0]-x_axis[2]<0.000001)&&(x_axis[1]-y_axis[0]<0.000001)) {
      temp = x_axis;
      x_axis = y_axis;
      y_axis = temp;
    }
    //Z,X,Y
    if ((y_axis[2]-z_axis[1]<0.000001)&&(z_axis[0]-x_axis[2]<0.000001)&&(x_axis[1]-y_axis[0]<0.000001)) {
      temp = x_axis;
      x_axis = y_axis;
      y_axis = z_axis;
      z_axis = temp;
      }
    
    x_axis.Print("x axis");
    y_axis.Print("y axis");
    z_axis.Print("z axis");
    Symmetry_Rotation.SetRowVector(x_axis,0);
    Symmetry_Rotation.SetRowVector(y_axis,1);
    Symmetry_Rotation.SetRowVector(z_axis,2);
    Symmetry_Rotation.Print("Rotation Matrix");


    RotAngle = acos((x_axis[0]+y_axis[1]+z_axis[2]-1)/2);
    //if (Params::Parameters().PrintLevel() > 0) 
      printf("Rotation Angle (in radian)     %10.6f          \n", RotAngle);  	  
    
    RotAngle = RotAngle*RadiansToDegrees;
    //if (Params::Parameters().PrintLevel() > 0) 
      printf("Rotation Angle (in degrees)     %10.6f          \n", RotAngle);  	  
    
    double De;
    De = sqrt((y_axis[2]-z_axis[1])*(y_axis[2]-z_axis[1])+(z_axis[0]-x_axis[2])*(z_axis[0]-x_axis[2])+(x_axis[1]-y_axis[0])*(x_axis[1]-y_axis[0]));
    //De = 2*sin(RotAngle*DegreesToRadians);
    if (De==0.0)
      De=0.000000001;
    
    
    RotVec[0] = (y_axis[2]-z_axis[1])/De; 
    RotVec[1] = (z_axis[0]-x_axis[2])/De; 
    RotVec[2] = (x_axis[1]-y_axis[0])/De;*/
  }


  // If this is our first time here, store the coordinates
  //if ( !use_old_axes) {
  //RotVec[0] = 1;
  //RotAngle = 0;
  //SetRotationAngle(RotAngle); 
  //RotationVector = RotVec;
  //}
  //else { // after first energy calculation, use ones previously stored
  //  printf("Using previously stored local coordinates for monomer %d\n",Monomer_index);
  //  RotAngle = RotationAngle;
  //  RotVec = RotationVector;
  //}

  //if (Params::Parameters().PrintLevel() > 0) 
  //  printf("Rotation vector      %10.6f          %10.6f           %10.6f          \n", 
  //	   RotVec[0], RotVec[1], RotVec[2]);  	  
  
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

    // loc_xyz[0] = 0.000000;  // local coords of the first atom
    //loc_xyz[1] = 0.000000;
    //loc_xyz[2] = 0.000000;
    for (int i=0;i<Natoms;i++) {
      
      // local x in global coord sys
      loc_xyz[3*i] = GetAtom(i).GetCoordinate(0) - CenterOfMass[0];
      // local y in global coord sys
      loc_xyz[3*i+1] = GetAtom(i).GetCoordinate(1) - CenterOfMass[1];
      // local z in global coord sys
      loc_xyz[3*i+2] = GetAtom(i).GetCoordinate(2) - CenterOfMass[2];

      //no longer need, now using center of mass coordinates for local coordinates.
      /*loc_xyz_g[3*i] = GetAtom(i).GetCoordinate(0) - CenterOfMass[0];
      
      // local y in global coord sys
      loc_xyz_g[3*i+1] = GetAtom(i).GetCoordinate(1) - CenterOfMass[1];
      
      // local z in global coord sys
      loc_xyz_g[3*i+2] = GetAtom(i).GetCoordinate(2) - CenterOfMass[2];
      
      //local x in local coord sys (dot product of x_axis and local x in global coord)
      loc_xyz[3*i] = x_axis[0]*loc_xyz_g[3*i] + x_axis[1]*loc_xyz_g[3*i+1] + x_axis[2]*loc_xyz_g[3*i+2];
      
      //local y in local coord sys (dot product of y_axis and local y in global coord)
      loc_xyz[3*i+1] = y_axis[0]*loc_xyz_g[3*i] + y_axis[1]*loc_xyz_g[3*i+1] + y_axis[2]*loc_xyz_g[3*i+2];
      
      //local z in local coord sys (dot product of z_axis and local z in global coord)
      loc_xyz[3*i+2] = z_axis[0]*loc_xyz_g[3*i] + z_axis[1]*loc_xyz_g[3*i+1] + z_axis[2]*loc_xyz_g[3*i+2];
      */
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

double Monomer::ReadMolProEnergy() {

  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  path = Params::Parameters().GetQMPath();  

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string out_filename = path + "/" + mN_ + ".out"; // Change to suffix of the molpro job - which I think we will just keep as .out...

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadMolProEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }
  
  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,17);

    // Search for final Molpro energy
    if ( match==" SETTING EMONOMER" ) { 
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }
  }
    
  if ( Params::Parameters().PrintLevel() > 0) printf("MolPro Obtained QM Monomer Energy = %15.9f\n",energy);
  
  // Close the output file
  infile.close();

  //Get CCSD(T) correction
  if(Params::Parameters().DoCCSDTCorrection()){
    out_filename = path + "/" + mN_ + ".CCSDT.out"; // Change to suffix of the molpro job - which I think we will just keep as .out...
  
    // Open the energy file		
    infile.open( out_filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadMolProEnergy : Cannot open file '%s'\n",
	     out_filename.c_str());
      exit(1);
    }
      

    // Read in the data
    double value = 0.0;
    while ( !infile.eof() ) {
      getline(infile,line);
      string match = line.substr(0,17);
      
      // Search for final Molpro energy
      if ( match==" SETTING EMONOMER" ) { 
	istringstream iss(line);
	string tmp;
	for (int i=0;i<3;i++)
	  iss >> tmp; // throw away text 
	iss >> value; // read CCSD(T) correction
	
      }

    }
    energy += value; //add the CCSD(T) correction to the energy
    infile.close(); 
  }
  
  //printf("%s energy = %f\n",mN_.c_str(),energy);
  return energy;

}

double Monomer::ReadQChemEnergy(bool MM_job) {

  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath(); 

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  string out_filename = path + "/" + mN_ + ".out"; 


  //using symmetry for monomers
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();

    //sym_fac=0 if the monomer is not the "Unique" monomer
    if(sym_fac == 0 && UseMonSym){
      char label[10];
      sprintf(label,"%d", SymMon);
      out_filename = path + "/m" + label + ".out";
    }


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
   
  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


    // Set up the filename, with the full path.  File is e.g. 'm3.out'
    string out_filename = path + "/" + mN_ + ".out";
    //string out_filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".out";
  
    //using symmetry for monomers
    bool UseMonSym = Params::Parameters().UseMonomerSymmetry();

    //sym_fac=0 if the monomer is not the "Unique" monomer
    if(sym_fac == 0 && UseMonSym){
      char label[10];
      sprintf(label,"%d", SymMon);
      //out_filename = Params::Parameters().GetMMPath() + "/m" + label + ".out";
     out_filename = path + "/m" + label + ".out";
    }

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


// JDH
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

// Read the atom-centered distributed multipole moments, as computed
// by CamCasp
void Monomer::ReadMultipoleMoments(Monomer& SymMonomer) {
  int Natoms_found = 0;

  // Set up the filename, with the full path.  
  string filename = Params::Parameters().GetMMPath() + "/" + mN_ + ".mom";

  //use symmetry for monomers
  bool UseMMSym = Params::Parameters().UseMMSymmetry();

  //use the input file for the monomer this monomer is symmetrical to
  if(!sym_fac && UseMMSym){
    char label[10];
    sprintf(label,"%d", SymMon);
    filename = Params::Parameters().GetMMPath() + "/m" + label + ".mom";;
  }

  // Open the moments file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Cluster::ReadMultipoleMoments(): Cannot open monomer %d multipole file '%s'",
	   Monomer_index,filename.c_str());
    exit(0);
  }

  // Create a list of atom names which is used to track our progress as
  // we read the moments file.

  // Due to using symmetry atoms are not nessarally in the correct order.
  // Using global index number to match atoms in file to atoms in monomer. 
  string AtmInt[Natoms];

  //printf("Atoms in AtomInt is [");
  for (int i=0;i<Natoms;i++){
    char label[10];
    if(UseMMSym)
      AtmInt[i] = SymMonomer.GetAtom(i).GetSymbol().c_str();    
    else
      AtmInt[i] = GetAtom(i).GetSymbol().c_str(); 
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;
    //printf("%s ",AtmInt[i].c_str());
  }
  //printf("]\n");

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
	  //AtomFound = true;
	  entryAtm = entry;
	  //iatom = n;

	  
	  if(UseMMSym){//if the MM is exploiting symmetry

	      iatom = GetSymmetricalAtom(n);
	      AtomFound = true;
	    //matching atoms in the symmetrical monomer to this monomer
	    //int j=0;
	    //while(!AtomFound&&j<Natoms){

	      //if symmetry is being exploited, that the order in input file may not be the same on list,
	      //correcting for that
	      //if(GetAtom(j).GetSymmetricalAtom() == SymMonomer.GetAtom(n).GetGlobalIndex()){
	      //iatom = j;
	    //}
	      //j++;
	      //}//end while j < Natoms
	  } else{//symmetry not exploited
	    AtomFound = true;
	    iatom = n; 
	  }
	}

	//if(!AtomFound){
	//printf("Error::Monomer::ReadMultipoleMoments:: Cannot find a symmetrical atom for %s in m%d\n",
	//       AtmInt[n].c_str(),Monomer_index);
	//exit(0);
	//}

      }// end of for n<Natoms


      
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
void Monomer::ReadPolarizabilities(Monomer& SymMonomer) {

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

  //use symmetry for monomers
  bool UseMMSym = Params::Parameters().UseMMSymmetry();

  //use the input file for the monomer this monomer is symmetrical to
  if(!sym_fac && UseMMSym){
    char label[10];
    sprintf(label,"%d", SymMon);
    filename = Params::Parameters().GetMMPath() + "/" +"m"+label+ ".pol";
    }
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
  //printf("Atoms in AtomInt is [");
  string AtmInt[Natoms];
  for (int i=0;i<Natoms;i++){
    char label[10];
    if(UseMMSym)
      AtmInt[i] = SymMonomer.GetAtom(i).GetSymbol().c_str();
    else
      AtmInt[i] = GetAtom(i).GetSymbol().c_str();
    sprintf(label,"%d",GetAtom(i).GetAtomIndex());
    AtmInt[i] += label;   
    //printf("%s ",AtmInt[i].c_str());

  }  
  //printf("]\n");

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
	//AtomFound = true;
	entryAtm = entry;
	//iatom = n;

	if(UseMMSym){
	  
	  iatom = GetSymmetricalAtom(n);
	  AtomFound = true;
	  //matching atoms in the symmetrical monomer to this monomer
	  //int i=0;
	  //while(!AtomFound && i<Natoms){
	  //  if(GetAtom(i).GetSymmetricalAtom() == SymMonomer.GetAtom(n).GetGlobalIndex()){
	  //    AtomFound = true;
	  //    iatom = i;
	  //    printf("iatom = %i\n",iatom);
	  //  }
	  //  i++;
	  //}
	}else{
	  AtomFound = true;
	  iatom = n;
	}
	//if(!AtomFound){
	//  printf("Error::Monomer::ReadMultipoleMoments:: Cannot find a symmetrical atom for %s in m%d\n",
	//	 AtmInt[n].c_str(),Monomer_index);
	// exit(0);
	//}
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

  //use the input file for the monomer this monomer is symmetrical to
  if(!sym_fac){  
    char label[10];
    sprintf(label,"%d", SymMon);
    filename =  Params::Parameters().GetMMPath() + "/m" + label + ".freq_pol";
    }
  
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
void Monomer::ReadDiagonalAnisotropicFreqPolarizabilities(Monomer& SymMonomer) {

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

  //flag for exploting symmetry for the MM terms
  bool UseMMSym = Params::Params::Parameters().UseMMSymmetry();
  
  //read the input monomer of the file it is symmetrical too if the MM is exploiting symmetry.
  if(!sym_fac && UseMMSym){
    char label[10];
    sprintf(label, "%d",SymMon);
    filename = Params::Parameters().GetMMPath() + "/m" + label +".freq_pol";
    }
  
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
    if(UseMMSym)
      AtmInt[i] = SymMonomer.GetAtom(i).GetSymbol().c_str();
    else
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
	  entryAtm = atom;
	  if(UseMMSym){

	    iatom = GetSymmetricalAtom(n);
	    AtomFound = true;
	    //matching atoms in the symmetrical monomer to this monomer
	    //int i=0;
	    //while(!AtomFound||i<Natoms){
	      //if(GetAtom(i).GetSymmetricalAtom() == SymMonomer.GetAtom(n).GetGlobalIndex()){
	      //iatom = i;
	      //AtomFound = true;
	      //}
	      //i++;
	      //}
	  }else{
	    iatom = n;
	    AtomFound = true;
	  }
	  //if(!AtomFound){
	  //  printf("Error::Monomer::ReadFreqPolarizabilities:: Cannot find a symmetrical atom for %s in m%d\n",
	  //	   AtmInt[n].c_str(),Monomer_index);
	  // exit(0);
	  //}
	  
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
      else if (GetAtom(i).GetSymbol() == "Ar") atom_type = "Ar";
      else if (GetAtom(i).GetSymbol() == "Kr") atom_type = "Kr";
      else {
	printf("Monomer::SetEmpiricalAtomDispersionType: Unknown Atom Type %s\n",GetAtom(i).GetSymbol().c_str());
	exit(1);
      }
      GetAtom(i).SetDispersionAtomType( atom_type );
  }
}

//Wrapper for create QM jobs
void Monomer::CreateQMJob(Monomer Monomers[], int NMon){

  //Flag for using symmetry for the Monomers
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();
  
  if ( Params::Parameters().GetQMType()==1 ) { // Qchem
    CreateQChemJob(Monomers, NMon);
  } 
  else if(Params::Parameters().GetQMType()==2) { // MolPro
    if(Params::Parameters().DoFreq()) {
      if(Params::Parameters().SingleFileMonomerHess()) {
        CreateMolProJob(Monomers, NMon);	 
      }
      else {
        CreateFiniteDifferenceMolProJob(Monomers, NMon, false);
	if(Params::Parameters().DoCCSDTCorrection()) {
	  CreateFiniteDifferenceMolProJob(Monomers, NMon, true);
	}
      }
    } 
    else {
      CreateMolProJob(Monomers, NMon);
      if(Params::Parameters().DoCCSDTCorrection()) {
        CreateMolProCCSDTJob(Monomers, NMon);
      }
    }
  }
  else if( Params::Parameters().GetQMType()==3){ // G09
    CreateG09Job(Monomers, NMon);
  }
  else if( Params::Parameters().GetQMType()==7){ // PSI4
    CreatePSI4Job(Monomers, NMon);
  }
  else {
    printf("Monomer::CreateQMJob: Unknown QM_type: %d\n",Params::Parameters().GetQMType() );
    exit(1);
  }
}


// Create a monomer qchem input file
void Monomer::CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job) {
  // We pass in Monomers/NMon so it can access list of other monomers.
  // This is necessary for embedding charges, for example.

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path;
  if (MM_job && !Params::Parameters().DoFreq()) 
    path = Params::Parameters().GetMMPath();
  else if(!MM_job && !Params::Parameters().DoFreq())
    path = Params::Parameters().GetQMPath(); 
  else if(MM_job && Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianMMPath();
  else if(!MM_job && Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

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

// Overload the CreateQChemJob for periodic systems with charge embedding  JDH NEW 
void Monomer::CreateQChemJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath();  
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

  fprintf(job,"%s\n",Params::Parameters().GetQChemRem().c_str());

  fprintf(job,"%s\n",Params::Parameters().GetQChemBasis().c_str());

  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    // includes a $basis section that list the basis to be used for each atom if requested.
    // Check to see if monomer is within the cut-off distance from the first monomer
    double dist = Monomers[1].FindDistance( Monomers[Monomer_index] ).Element(0);
    //for (int imon=1;imon<=NMon;imon++) {
    //  dist = min(dist , Monomers[1].FindDistance( Monomers[] ) );
    //}

    fprintf(job,"$basis\n");
    // Now loop over all the atoms and include the appropriate basis
    for (int iatom=1;iatom<=Monomers[Monomer_index].GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%s %d\n", GetAtom(iatom-1).GetSymbol().c_str(), iatom);
      if ( dist <= Params::Parameters().GetMixedBasisCutOff() ) {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
      }
      fprintf(job,"****\n");
    }
    fprintf(job,"$end\n");

  }

  // Optionally print $external_charges section
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"$external_charges\n");
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }

    // Now print the image charges .
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( GetIndex() != MonomerImages[i].GetIndex() ) {
	  if ( MonomerImages[i].GetUseInEmbedding() ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    } else {
      printf("DEBUG Monomer::CreateQChemJob called the periodic code for creating a monomer, but you are not using a perioid system \n");
      exit(1);
    }

    fprintf(job,"$end\n");
  }
  fclose(job);
}

void Monomer::CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon) {
  string path;
  path = Params::Parameters().GetQMPath();
  string charge_file;
  charge_file = path + "/" + mN_ + ".ch";
  string charge_file_local = mN_ + ".ch";

  
  FILE *chrg;
  if (( chrg = fopen(charge_file.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolProChargeEmbeddingJob : cannot open file '%s'\n",charge_file.c_str());
    exit(1);
  }
  fprintf(chrg,"Embedding charges for monomer %d\n", Monomer_index ); 
  //count the charges
  int number_charges = 0;
  for (int i=1;i<=NMon;i++) {
    if ( i != Monomer_index) {
	number_charges += Monomers[i].GetNumberOfAtoms();
    }
  }

  fprintf(chrg,"%d\n", number_charges);
  
  // Now Print the charges to the molpro charge file:
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(chrg);	
    }
  }
  
  fclose(chrg);
}

void Monomer::CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images) {
  string path;
  path = Params::Parameters().GetQMPath();
  string charge_file;
  charge_file = path + "/" + mN_ + ".ch";
  string charge_file_local = mN_ + ".ch";

  FILE *chrg;
  if (( chrg = fopen(charge_file.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolProChargeEmbeddingJob : cannot open file '%s'\n",charge_file.c_str());
    exit(1);
  }
  fprintf(chrg,"Embedding charges for monomer %d\n", Monomer_index ); 
  //count the charges
  int number_charges = 0;
  for (int i=1;i<=NMon;i++) {
    if ( i != Monomer_index) {
      number_charges += Monomers[i].GetNumberOfAtoms();
    }
  }

  // Add in the number of charges from periodic images if requested...
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1; i<=NMon_images; i++) {
      if ( MonomerImages[i].GetUseInEmbedding() ) {
	number_charges += MonomerImages[i].GetNumberOfAtoms();
      }
    }
    
  }
  
  fprintf(chrg,"%d\n", number_charges);
  
  // Now Print the charges to the molpro charge file:
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(chrg);	
    }
  }
  
  // Print the charges from periodic images if needed...
  if ( Params::Parameters().IsPeriodic() ) {
    for (int i=1; i<=NMon_images; i++) {
      if ( GetIndex() != MonomerImages[i].GetIndex() ) {
	if ( MonomerImages[i].GetUseInEmbedding() ) {
	  MonomerImages[i].PrintEmbeddingCharges(chrg);
	}
      }
    }
  }
  
  fclose(chrg);
}

void Monomer::CreateMolProJob(Monomer Monomers[],int NMon){

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath(); 
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/" + mN_ + ".inp"; 

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolproJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  //input header
  fprintf(job,"*** Monomer %d\n",Monomer_index);
  //fprintf(job,"memory,500 m\n")
  fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
  fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
  fprintf(job,"nogprint;\n");
  fprintf(job,"gdirect;\n");
  fprintf(job,"symmetry,nosym;\n");
  fprintf(job,"orient,noorient;\n");
  fprintf(job,"angstrom\n");

  //Print Cartesian coords
  fprintf(job,"GEOMETRY={\n");
  PrintMolProMonomerCastesian(job,1) ;
  fprintf(job,"}\n\n");

  //charge and spin (spin = multiplicity - 1)
  fprintf(job,"SET,CHARGE=%i\n",charge);
  fprintf(job,"SET,SPIN=%i\n\n",spin-1);


  //CBS limit calculations
  if(Params::Parameters().DoCBS()){
    //first basis set
    //dimer section
    fprintf(job,"!First basis\n");
    fprintf(job,"!-----\n");
    fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
    fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    fprintf(job,"E_HF1 = energy\n");
    fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    fprintf(job,"E_MP2_1 = energy\n\n");



    //second basis set
    fprintf(job,"!Second basis\n");
    fprintf(job,"!-----\n");
    fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
    fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    fprintf(job,"E_HF2 = energy\n");
    fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    fprintf(job,"E_MP2_2 = energy\n\n");

    //CBS Extrapole
    int basis1 = Params::Parameters().CBSBasis1();
    int basis2 = Params::Parameters().CBSBasis2();
    if(basis1 >= basis2){
      printf("Monomer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
      exit(0);
    }
    fprintf(job,"!CBS extrapolate\n");
    fprintf(job,"!-----\n");
    fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	   basis2,basis1);
    fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
	   basis2,basis1,basis2,basis1);
    fprintf(job,"Emonomer = E_HF_inf + E_Coor_inf\n");
  }
  else{
    //print rem section
    fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
    fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
    
    //print energy in output
    
    fprintf(job,"Emonomer=energy\n");
    fprintf(job,"\n");
    //fprintf(job,"!summary\n");
    //fprintf(job,"!--------\n");
    //fprintf(job,"show,Emonomer\n");
    
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    
  }
  fclose(job);
}

void Monomer::CreateMolProCCSDTJob(Monomer Monomers[],int NMon){

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath(); 
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/" + mN_ + ".CCSDT.inp";
  //printf("%s\n",filename.c_str()); 

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolCCSDTProJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  //input header
  fprintf(job,"*** Monomer %d CCSD(T) correction\n",Monomer_index);
  //fprintf(job,"memory,500 m\n")
  fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
  fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
  fprintf(job,"symmetry,nosym;\n");
  fprintf(job,"orient,noorient;\n");
  fprintf(job,"angstrom\n");  

  //Print Cartesian coords
  fprintf(job,"GEOMETRY={\n");
  PrintMolProMonomerCastesian(job,1) ;
  fprintf(job,"}\n\n");


  //charge and spin (spin = multiplicity - 1)
  fprintf(job,"SET,CHARGE=%i\n",charge);
  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

  //basis
  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
  //MP2 
  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
  fprintf(job,"E_MP2 = energy\n");
  if(Params::Parameters().DoForces()){
    fprintf(job,"FORCE\n");
  }
  if(Params::Parameters().DoFreq() ){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
  }
  //CCSDT
  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
  fprintf(job,"E_CCSDT = energy\n");
  if(Params::Parameters().DoForces()){
    fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
  }
  if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
  }
  fprintf(job,"\n");
  fprintf(job,"Emonomer = E_CCSDT - E_MP2\n");
  fclose(job);
}

void Monomer::CreateFiniteDifferenceMolProJob(Monomer Monomers[],int NMon,bool CCSDT){

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath(); 
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  //directory for finite difference steps
  path += "/" + mN_;

  //CCSD(T) jobs
  if(CCSDT){
    //directory for finite difference steps
    path += ".CCSDT";
  }

  //create monomer directory if it doesn't already exist
  struct stat st;
  if( stat(path.c_str(),&st) != 0){
    string cmd;
    //printf("Directory %s not found.  Creating it.\n",
    //   dimerpath.c_str());
    cmd = "mkdir " + path;	 
    //printf("%s\n",cmd.c_str());
    system(cmd.c_str());
  }

  FILE *job;
  string filename;
  for(int i=0;i<3*Natoms;i++){

    char num[10];
    sprintf(num,"%d",i);

    //Plus Geometry
    string filename = path + "/" + mN_ + "+" + num + ".inp";
    // Open the input file for writing
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateFiniteDifferenceMolProJob() : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
    //input header
    fprintf(job,"*** Monomer %d\n",Monomer_index);
    //fprintf(job,"memory,500 m\n")
    fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
    fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
    fprintf(job,"symmetry,nosym;\n");
    fprintf(job,"orient,noorient;\n");
    fprintf(job,"angstrom\n");
    
    //Print Cartesian coords
    fprintf(job,"GEOMETRY={\n");
    Vector Coords = GetCoordinates();
    Coords[i] += 0.001;
    PrintMolProMonomerCastesian(job,Coords,1);
    fprintf(job,"}\n\n");
    
    //charge and spin (spin = multiplicity - 1)
    fprintf(job,"SET,CHARGE=%i\n",charge);
    fprintf(job,"SET,SPIN=%i\n\n",spin-1);


    //CCSD(T) correction
    if(CCSDT){
      //basis
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
      //MP2 contribution
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
      fprintf(job,"E_MP2 = energy\n");
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");	
      }
      //CCSDT
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
      fprintf(job,"E_CCSDT = energy\n");
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
      }
      fprintf(job,"\n");
      fprintf(job,"Emonomer = E_CCSDT - E_MP2\n");
      
    //CBS limit calculations
    }else if(Params::Parameters().DoCBS()){
      //first basis set
      fprintf(job,"!First basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_HF1 = energy\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_MP2_1 = energy\n\n");
      
      
      
      //second basis set
      fprintf(job,"!Second basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_HF2 = energy\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_MP2_2 = energy\n\n");
      
      //CBS Extrapole
      int basis1 = Params::Parameters().CBSBasis1();
      int basis2 = Params::Parameters().CBSBasis2();
      if(basis1 >= basis2){
	printf("Monomer::CreateFiniteDifferenceMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	exit(0);
      }
      fprintf(job,"!CBS extrapolate\n");
      fprintf(job,"!-----\n");
      fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	      basis2,basis1);
      fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
	      basis2,basis1,basis2,basis1);
      fprintf(job,"Emonomer = E_HF_inf + E_Coor_inf\n");
    }

    else{
      //print rem section
      fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
      fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
      
      //print energy in output
      fprintf(job,"Emonomer=energy\n");
      fprintf(job,"\n");
      //fprintf(job,"!summary\n");
      //fprintf(job,"!--------\n");
      //fprintf(job,"show,Emonomer\n");
      
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      
    }
    fclose(job);

    //Minus Geometry
    filename = path + "/" + mN_ + "-" + num + ".inp";
    // Open the input file for writing
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateFiniteDifferenceMolProJob() : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
    //input header
    fprintf(job,"*** Monomer %d\n",Monomer_index);
    //fprintf(job,"memory,500 m\n")
    fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
    fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
    fprintf(job,"symmetry,nosym;\n");
    fprintf(job,"orient,noorient;\n");
    fprintf(job,"angstrom\n");
    
    //Print Cartesian coords
    fprintf(job,"GEOMETRY={\n");
    Coords[i] -= 0.002;
    PrintMolProMonomerCastesian(job,Coords,1);
    fprintf(job,"}\n\n");

    //charge and spin (spin = multiplicity - 1)
    fprintf(job,"SET,CHARGE=%i\n",charge);
    fprintf(job,"SET,SPIN=%i\n",spin-1);

        //CCSD(T) correction
    if(CCSDT){
      //basis
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
      //MP2 contribution
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
      fprintf(job,"E_MP2 = energy\n");
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");	
      }
      //CCSDT
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
      fprintf(job,"E_CCSDT = energy\n");
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
      }
      fprintf(job,"\n");
      fprintf(job,"Emonomer = E_CCSDT - E_MP2\n");
      
    //CBS limit calculations
    }else if(Params::Parameters().DoCBS()){
      //first basis set
      //dimer section
      fprintf(job,"!First basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_HF1 = energy\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_MP2_1 = energy\n\n");
      
      
      
      //second basis set
      fprintf(job,"!Second basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_HF2 = energy\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      fprintf(job,"E_MP2_2 = energy\n\n");
      
      //CBS Extrapole
      int basis1 = Params::Parameters().CBSBasis1();
      int basis2 = Params::Parameters().CBSBasis2();
      if(basis1 >= basis2){
	printf("Monomer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	exit(0);
      }
      fprintf(job,"!CBS extrapolate\n");
      fprintf(job,"!-----\n");
      fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	      basis2,basis1);
      fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
	      basis2,basis1,basis2,basis1);
      fprintf(job,"Emonomer = E_HF_inf + E_Coor_inf\n");
    }
    else{
      //print rem section
      fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
      fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
      
      //print energy in output
      fprintf(job,"Emonomer=energy\n");
      fprintf(job,"\n");
      //fprintf(job,"!summary\n");
      //fprintf(job,"!--------\n");
      //fprintf(job,"show,Emonomer\n");
      
      if(Params::Parameters().DoFreq()){
	if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	  fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
	else
	  fprintf(job,"FORCE\n");
      }
      
    }
    fclose(job);
  }
  /*
  //creating csh script to run all jobs
  filename = path + "/" + mN_ + ".csh";

  // Open the input file for writing
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateMolProJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  fprintf(job,"#!/bin/csh\n\n");
  fprintf(job,"rm -f *out *xml\n");
  fprintf(job,"foreach s(m*inp)\n");
  fprintf(job,"  molpro $s\n");
  fprintf(job,"end\n");
  fclose(job);
  */
}

// Wrapper for creating MM jobs
void Monomer::CreateMMJob(Monomer Monomers[], int NMon) {

  //Flag for using symmetry for the Monomers: used by Tinker
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();

  //Flag for using symmetry for the AIFF
  bool UseAIFFSym = Params::Parameters().UseMMSymmetry();
  
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
    if(sym_fac!=0||!UseMonSym)//symmetry factor != 0 or monomer symmetry 
      CreateTinkerJob(Monomers, NMon);
  }
  else if ( Params::Parameters().GetMMType()==2 ) { //Orient
    if(sym_fac!=0 || !UseAIFFSym)
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


  string path;
  if(Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianMMPath();
  else
   path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  // Set up the filenames, with the full path.  
  // Files are e.g. 'm1.xyz' and 'm1.key'
  string xyzfile = path + "/" + mN_ + ".xyz";
  string keyfile = path + "/" + mN_ + ".key";;

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



// Wrapper that controls which type of QM job is executed
// Returns command string for running the job
string Monomer::RunQMJob(bool CCSDT) {
  string job;

  if (Params::Parameters().GetQMType()==1) // QChem
    job = RunQChemJob();
  else if (Params::Parameters().GetQMType()==2){ // MolPro
    if(Params::Parameters().DoFreq() && !Params::Parameters().SingleFileMonomerHess())
      job = RunFiniteDifferenceMolProJob(); 
    else
      job = RunMolProJob(CCSDT);
  }
  else if(Params::Parameters().GetQMType()==3) { // G09
    job = RunG09Job();
  }  
  else if(Params::Parameters().GetQMType()==7) { //PSI4
    job = RunPSI4Job();
  }
  else {
    printf("ERROR: Monomer::RunQMJob(): Unknown QM program. QM_Type = %i\n",Params::Parameters().GetQMType());
    exit(1);
  }
  return job;
}

// Returns command string for running the job
string Monomer::RunQChemJob(bool MM_job) {

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path;
  if (MM_job && !Params::Parameters().DoFreq() ) {
      path = Params::Parameters().GetMMPath();
  }
  if (!MM_job && !Params::Parameters().DoFreq() )  {
      path = Params::Parameters().GetQMPath();
  }
  if (MM_job && Params::Parameters().DoFreq() ) {
      path = Params::Parameters().GetHessianMMPath();
  }
  else if (!MM_job && Params::Parameters().DoFreq() )  {
      path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  
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

string Monomer::RunMolProJob(bool CCSDT) {

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath();
  if (Params::Parameters().DoFreq() )  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string infile = path + "/" + mN_ + ".inp"; 
  string outfile = path + "/" + mN_ + ".out";
  string xmlfile =  path + "/" + mN_ + ".xml";
  if(CCSDT){
    infile =  path + "/" + mN_ + ".CCSDT.inp"; 
    outfile = path + "/" + mN_ + ".CCSDT.out";
    xmlfile = path + "/" + mN_ + ".CCSDT.xml";
  }
  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  //remove previous .out and .xml
  cmd += "rm -f " + outfile + " " + xmlfile + ";";
  
  //run molpro
  cmd += "molpro " + infile + " --no-xml-output;";

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  
  //printf("%s\n",cmd.c_str());
  return cmd;

}

// Returns command string for running Molpro jobs
string Monomer::RunFiniteDifferenceMolProJob() {
  
  // Set up the filename, with the full path. 
  string path = Params::Parameters().GetQMPath(); 
  if(Params::Parameters().DoFreq() )
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasihamN_ rmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  path += "/" + mN_ ;

  string cshfile = path + "/" + mN_ + ".csh"; 

  //printf("%s\n",path.c_str());
  //printf("%s\n",cshfile.c_str());

  // First command, change to local directory 
  string cmd = "cd " + path + ";";
   
  //run shell script
  cmd += "csh " + cshfile + ";";

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("%s\n",cmd.c_str());

  return cmd;
}

// Returns string for running individual 
string Monomer::RunFiniteDifferenceMolProJobs(int i,bool plus,bool CCSDT){

  // Set up the filename, with the full path. 
  string path = Params::Parameters().GetQMPath(); 
  if(Params::Parameters().DoFreq() )
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasihamN_ rmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  path += "/" + mN_ ;

  //CCSD(T) jobs
  if(CCSDT){
    //directory for finite difference steps
    path += ".CCSDT";
  }

  char num[10];
  sprintf(num,"%d",i);

  string infile,outfile,xmlfile,logfile;
  if(plus){
    infile = path + "/" + mN_ + "+" +  num + ".inp";
    outfile = path + "/" + mN_ + "+" +  num + ".out";  
    xmlfile = path + "/" + mN_ + "+" +  num + ".xml"; 
    logfile = path + "/" + mN_ + "+" +  num + ".log"; 
  }else{
    infile = path + "/" + mN_ + "-" +  num + ".inp";
    outfile = path + "/" + mN_ + "-" +  num + ".out";  
    xmlfile = path + "/" + mN_ + "-" +  num + ".xml";  
    logfile = path + "/" + mN_ + "+" +  num + ".log";
  }
  
  // First command, change to local directory 
  string cmd = "cd " + path + ";";


  //deleting prevous outputs
  ifstream check;
  check.open( outfile.c_str() );
  if ( check.is_open() ) {
  cmd += "rm -f " + outfile + ";";
  }
  check.close();

  check.open( xmlfile.c_str() );
  if ( check.is_open() ) {
    cmd += "rm -f " + xmlfile + ";";
  }
  check.close();

  check.open( logfile.c_str() );
  if ( check.is_open() ) {
    cmd += "rm -f " + logfile + ";";
  }
  check.close();

  //run shell script
  cmd += "molpro " + infile + ";";

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("%s\n\n",cmd.c_str());

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
    printf("ERROR: Monomer::RunMMJob(): Unknown MM program. MM_Type = %i\n",
	   Params::Parameters().GetMMType());
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
  if(Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs())
    job_path = Params::Parameters().GetHessianMMPath();


  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "analyze " + infile;
  cmd += " e > ";
  cmd += outfile + "; ";

  // If doing force job, compute the gradient
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    // Run the job
    cmd += "minimize " + infile;
    cmd += " 100000000 > ";
    string force_file = mN_ + ".force";
    cmd += force_file + ";";

    // Remove tmp file and extra geom file created by minimize job
    cmd += "rm -f " + infile;
    cmd += "_2;";
  }
    // If doing freq job, compute the freq
    if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
    // Run the job
    cmd += "vibrate " + infile;    
    cmd += " 1 > ";
    string freq_file = mN_ + ".freq";
    cmd += freq_file + ";";

    // Remove tmp file and extra geom file created by vibrate job
    cmd += "rm -f " + mN_ + ".001;"; 
    //cmd += "_2;";
  }

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  // Actual job running & checking of status have been moved to main.C 
  // to simplify parallel implementation.

  return cmd;
}

void Monomer::ReadQMResults(Monomer& SymMonomer) {
  // Read in the Energy
  if(Params::Parameters().GetQMType() == 1)//Qchem
    Energy_QM = ReadQChemEnergy();
  else if(Params::Parameters().GetQMType() == 2) // MolPro
    Energy_QM = ReadMolProEnergy();
  else if(Params::Parameters().GetQMType() == 3) // G09
    Energy_QM = ReadG09Energy();
  else if(Params::Parameters().GetQMType() == 7) // PSI4
    Energy_QM = ReadPSI4Energy();
  else {
    printf("Monomer::ReadQMResults() - QM type = %d is unknown\n",Params::Parameters().GetQMType());
    exit(1);
  }

  double MP2_dispersion_correction;
  if ( Params::Parameters().DoMP2DispersionCorrection() )
    MP2_dispersion_correction = ComputeInteratomicMP2TwoBodyDispersionCorrection();

  // Optionally, read in the Forces
  if ( Params::Parameters().DoForces() ||Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetQMGradient();
  }

  // Optionally, read in the Hessian
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() ){
    if(sym_fac !=0 || !Params::Parameters().UseMonomerSymmetry()){
      SetQMHessian();
    }else{
      //seperate function for setting hessian for monomers under symmetry
      Hess_QM = SetSymmetricalHessian(1,SymMonomer);
    }
  }

}

void Monomer::ReadMMResults(Monomer& SymMonomer) {


  // If we're using the AIFF, we read in the monomer force field
  // parameters here.  SetMMEnergy either sets the monomer energy zero
  // (AIFF) or reads from an output file (other types).

  if ( Params::Parameters().GetMMType()==2 ) {// AIFF
    // Read in the ab initio force field parameters
    ReadMultipoleMoments(SymMonomer);


    ReadPolarizabilities(SymMonomer);

    // Next two have been replaced with improved AIFF routines
    //ReadIsotropicDipolePolarizability();
    //ReadDispersionCoefficients();

    if ( Params::Parameters().UseDiagonalAnisotropicPolarizabilities() )
      ReadDiagonalAnisotropicFreqPolarizabilities(SymMonomer);
    else 
      ReadFreqPolarizabilities();// by shuhao
  
    SetEmpiricalAtomDispersionType(); // used to figure out dispersion damping
  }

  // Read in the Energy
  SetMMEnergy();

  // Optionally read in the Forces
  if ( Params::Parameters().DoForces() /*||Params::Parameters().DoFreq()*/ || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient();
  }

  // Optionally read in the hessian
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() &&
       ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs())) {
    if(sym_fac != 0 || !Params::Parameters().UseMonomerSymmetry())  
      SetMMHessian(); 
    else
      //seperate function for setting hessian for monomers under symmetry
      Hess_MM = SetSymmetricalHessian(2,SymMonomer);
  }  
}

// Get the QM Gradient, wrapper routine
void Monomer::SetQMGradient() {
  string path = Params::Parameters().GetQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  Grad_QM = ReadGradient(path,1);
  QM_Grad_Init = 1;
}

// Get the MM Gradient, wrapper routine
void Monomer::SetMMGradient() {
  if (Params::Parameters().GetMMType() == 1) { // Tinker 
    string path = Params::Parameters().GetMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
       path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

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
  Vector TempGrad(3*Natoms);

  //Getting the path, filename may be changed if monomer symmetry is being used
  string filename = path + "/" + mN_;

  //if symmetry is being used
  bool UseMonSym = Params::Parameters().UseMonomerSymmetry();

  //if monomer symmetry is being used and this monomer is not the "Unique" monomer
  if(!sym_fac && UseMonSym){
    char label[10];
    sprintf(label, "%d",SymMon);
    filename = path + "/m" + label;
  }
  //Getting the extention
  if (type == 2) // Tinker MM job
    filename += ".force";
  else if(Params::Parameters().GetQMType()==3) // G09
    filename += ".log";
  else // other
    filename += ".out";
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadGradient : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  // Read in the data - search for the "Nuclear forces:" string
  string line;

  //for CBS, number of Grad read
  int NumberGrad = 0;
  while ( !infile.eof() ) {
    getline(infile,line);
      // Search for final q-chem or tinker energy
    if(Params::Parameters().GetQMType() == 1 || type == 2) {
      if ( line==" Nuclear forces:" ) {
	getline(infile,line); // throw away header line
	
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp; //for CBS, number of Grad read

	  iss >> tmp; // throw away the atom index
	  for (int j=0;j<3;j++) {
	    iss >> TempGrad[3*i+j];//gradients will be rotated if need due to symmetry 
	    // Store the gradient elements
	  }
	}
	break;
      }
    }//End Qchem and Tinker.
    else if(Params::Parameters().GetQMType() == 2){ //Molpro
      //Extrapolate gradients to the CBS limit
      if(Params::Parameters().DoCBS()){
	double basis1 = double(Params::Parameters().CBSBasis1());
	double basis2 = double(Params::Parameters().CBSBasis2());
	double ExpTerm = exp(1.54*(basis2-basis1));
	

	if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	  getline(infile,line);// throw away header line
	  NumberGrad++;
	  for (int i=0;i<Natoms;i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=0;j<3;j++) {
	      double entry;
	      iss >> entry;
	      //HF in the first basis  
	      if(NumberGrad == 1)
		entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      //MP2 in the first basis
	      else if(NumberGrad==2)		  
		entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      //HF in the second basis
	      else if(NumberGrad==3)
		entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      //MP2 in the second basis
	      else if(NumberGrad==4)
		entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      TempGrad[3*i+j] +=entry;//gradients will be rotated if needed due to symmetry 

	      // Store the gradient elements
	    }

	  }	
	  if(NumberGrad == 4){
	    break;
	  }
	}

      }


      else{
	if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	  getline(infile,line);// throw away header line
	  
	  for (int i=0;i<Natoms;i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=0;j<3;j++) {
	      iss >> TempGrad[3*i+j];//gradients will be rotated if need due to symmetry 
	      // Store the gradient elements
	    }
	  }
	  break;
	}
      }
    }
    else if(Params::Parameters().GetQMType()==3) { //G09 Watit
      if(line.find("Forces (Hartrees/Bohr)") != string::npos) {
	getline(infile,line);
	getline(infile,line);
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  istringstream iss(line);
          string tmp;
	  iss >> tmp;
	  iss >> tmp;
	  for (int j=0;j<3;j++) {
	    iss >> TempGrad[3*i+j];
	    TempGrad[3*i+j]=-TempGrad[3*i+j];
	  }
	}
      }
    }
    else if(Params::Parameters().GetQMType()==7) { //PSI4 CSG
      if(line.find("-Total Gradient:") != string::npos) {
        getline(infile,line);
	getline(infile,line);
	//getline(infile,line);
	//cout << line << endl; // make sure I'm printing out the gradients
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  //cout << line << endl; // make sure I'm printing out the gradients
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp;
//	  cout << tmp << endl;
	  for (int j=0;j<3;j++) {
            // psi4 prints the gradients in kcal/mol*A but hmbi wants hartree/bohr
            // Old way
	    iss >> TempGrad[3*i+j];
            // PSI4 prints the gradients in kcal/mol*A
            // Need to convert gradient to hartree/bohr
            double temp = TempGrad[3*i+j]*(1/(627.509*1.8897259886));
            TempGrad[3*i+j] = temp;
	    //TempGrad[3*i+j]=-TempGrad[3*i+j];
	  //  cout << TempGrad[3*i+j] << endl;
          }
        }
      }
    }
  }
  infile.close();

  //R.Print("Rotation Matrix");

  //Include the CCSD(T) Correction
  //CCSD(T) has not been adapted to q-chem 
  if(Params::Parameters().DoCCSDTCorrection() && type == 1){
    //Getting the path, filename may be changed if monomer symmetry is being used
    filename = path + "/" + mN_;
    
    //if monomer symmetry is being used and this monomer is not the "Unique" monomer
    if(!sym_fac && UseMonSym){
      char label[10];
      sprintf(label, "%d",SymMon);
      filename = path + "/m" + label;
    }
    
    filename += ".CCSDT.out";

    // Open the force file
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadGradient : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
    
    int NumberGrad = 0;
    while ( !infile.eof() ) {
      getline(infile,line);
      if(line.substr(0,60) ==" Atom          dE/dx               dE/dy               dE/dz"){
	NumberGrad++;
	if(NumberGrad==1)
	  getline(infile,line);// throw away header line
	for (int i=0;i<Natoms;i++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  for (int j=0;j<3;j++) {
	    double entry;
	    iss >> entry;
	    //change the sign of the MP2 part
	    if(NumberGrad == 1)
	      entry *= -1;
	    TempGrad[3*i+j] +=entry;//gradients will be rotated if needed due to symmetry 
	  }
	}
      }
    }
    //Wrong number of Gradients
    if(!(NumberGrad==2)){
      printf("Monomer::ReadGradient : At Least One Gradient Not Found or Incorrect '%s'\n",
	     filename.c_str());
      printf("NumberGrad = %i\n",NumberGrad);
      exit(1);
    }
    infile.close();

  }

  //rotating monomers due to symmetry
  Matrix R = Symmetry_Rotation;
  R.Transpose();
  for (int i=0;i<Natoms;i++){
      for(int j=0;j<3;j++){
	if(!UseMonSym){//if only the "Unique" monomer output file was made, rotate that monomer into the orientation of the this monomer
	  grad[3*i+j] = TempGrad[3*i+j];
	}else{
	  grad[3*i+j] = R(j,0)*TempGrad[3*Sym_Atom[i]]+R(j,1)*TempGrad[3*Sym_Atom[i]+1]+R(j,2)*TempGrad[3*Sym_Atom[i]+2];
        }
      }
   }

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

// Print out the Embedding Charges for molpro 
void Monomer::PrintMolProEmbeddingCharges(FILE* outfile) {
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintMolProEmbeddingCharges(outfile);
  }
}

// Print out the Cartesian coordinate for molpro
void Monomer::PrintMolProMonomerCastesian(FILE *outfile,int start){
  for (int i=0;i<Natoms;i++) {
    fprintf(outfile, "%s%i  %10.6f  %10.6f  %10.6f \n", 
      GetAtom(i).GetSymbol().c_str(),i+start,
      GetAtom(i).GetCoordinate(0),GetAtom(i).GetCoordinate(1), 
      GetAtom(i).GetCoordinate(2));
  }
}

// Print out the Cartesian coordinate for molpro using coordinates definded by CCoords
void Monomer::PrintMolProMonomerCastesian(FILE *outfile,Vector Coords,int start){

  if(Coords.GetLength() != 3*Natoms){
    printf("ERROR:Monomer:PrintMolProMonomerCastesian() : Length of Coords not same as 3*Natom. Coords =%i\n",
	   Coords.GetLength());
    exit(0);
  }

  for (int i=0;i<Natoms;i++) {
    fprintf(outfile, "%s%i  %10.6f  %10.6f  %10.6f \n", 
      GetAtom(i).GetSymbol().c_str(),i+start,
	    Coords[3*i],Coords[3*i+1],Coords[3*i+2]);
  }
}

// Print out the Cartesian coordinate for PSI4 //CSG
void Monomer::PrintPSI4MonomerCartesian(FILE *outfile,int start){
  for (int i=0;i<Natoms;i++) {
    fprintf(outfile, "%s%i  %10.6f  %10.6f  %10.6f \n", 
      GetAtom(i).GetSymbol().c_str(),i+start,
      GetAtom(i).GetCoordinate(0),GetAtom(i).GetCoordinate(1), 
      GetAtom(i).GetCoordinate(2));
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

Vector Monomer::GetCoordinates(){

  Vector Coords(3*Natoms);
  for(int i=0;i<Natoms;i++)
    for(int xyz=0;xyz<3;xyz++)
      Coords[3*i+xyz] = Atom_List[i].GetPosition(xyz);

  return Coords;

}

//Find if momoners are symmetrical to other monomers
//Symmetrical Rotation Matrix found here may be corrected in 
//Cluster::CorrectingMonomerRotations
bool Monomer::SymmetryCheck(Monomer Monomers[]){
  

  //center of mass coordinates and inertia tenser
  Matrix COM(Natoms,3);
  COM.Set();
  //inertia tensor which will be used to determine how to rotate the molecule
  Matrix Inertia(3,3);
  Inertia.Set();
  // printf("center of for m%i are: %f %f %f\n",Monomer_index,GetCenterOfMass(0),GetCenterOfMass(1),GetCenterOfMass(2));
  //printf("center of mass coordinates are\n");
  for(int i=0;i<Natoms;i++) {
    COM(i,0) = Atom_List[i].GetPosition(0)-GetCenterOfMass(0);// x center of mass
    COM(i,1) = Atom_List[i].GetPosition(1)-GetCenterOfMass(1);// y center of mass
    COM(i,2) = Atom_List[i].GetPosition(2)-GetCenterOfMass(2);// z center of mass
    // printf( "%-2s %10.6f  %10.6f  %10.6f\n",
    //	    Atom_List[i].GetSymbol().c_str(),COM(i,0),COM(i,1),COM(i,2));
    
    //inertia tensor
    double mass = Atom_List[i].GetAtomicMass();
    Inertia(0,0) += mass*(COM(i,1)*COM(i,1)+COM(i,2)*COM(i,2));//m(y*y + z*z)
    Inertia(1,1) += mass*(COM(i,0)*COM(i,0)+COM(i,2)*COM(i,2));//m(x*x + z*z);
    Inertia(2,2) += mass*(COM(i,0)*COM(i,0)+COM(i,1)*COM(i,1));//m(x*x + y*y);
    Inertia(0,1) -= mass*COM(i,0)*COM(i,1);//m*x*y;
    Inertia(0,2) -= mass*COM(i,0)*COM(i,2);//m*x*z;
    Inertia(1,2) -= mass*COM(i,1)*COM(i,2);//m*y*z;
  }
  Inertia(1,0) = Inertia(0,1);
  Inertia(2,0) = Inertia(0,2);
  Inertia(2,1) = Inertia(1,2);
  //diagonalize inertia matrix in amu*A^2
  Vector MomentsOfInertia = Inertia.Diagonalize();
  //printf("the moments of inertia for m%i are: %10.6f  %10.6f  %10.6f\n",Monomer_index, MomentsOfInertia[0],MomentsOfInertia[1],MomentsOfInertia[2]);
  //The principal axes of rotation (for an asymmetric top)	
  Vector x_axis = Inertia.GetColumnVector(0);
  Vector y_axis = Inertia.GetColumnVector(1);
  Vector z_axis = Inertia.GetColumnVector(2);

  //Inertia.Print("Inertia matrix");
  //Inertia.Inverse();
  //Inertia.Print("Inverse");
 
  //Set the axes to the principal axes of rotation to make local axes.
  Matrix LocalCoor(Natoms,3);
  //printf("Using dot products,Internal coordinates are for monomer %i:\n", Monomer_index);
  for(int i=0;i<Natoms;i++) {
    LocalCoor(i,0) = COM(i,0)*x_axis[0]+COM(i,1)*x_axis[1]+COM(i,2)*x_axis[2];
    LocalCoor(i,1) = COM(i,0)*y_axis[0]+COM(i,1)*y_axis[1]+COM(i,2)*y_axis[2];
    LocalCoor(i,2) = COM(i,0)*z_axis[0]+COM(i,1)*z_axis[1]+COM(i,2)*z_axis[2];
    //printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	    Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  }


  //Comparing local coordinates to other monomers to determine symmetry
  //Difference schemes for Spherical, Symmetric, Asymetric and Linear Tops

  bool EquivalencyFound = 0; //bool to tell us if an equivalent monomer has been found
  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance 
  
  //Find equivalency for single atom monomers (NOT IMPLEMENTED)
  if(Natoms<2){
    printf("Error::Monomer::SymmetryCheck:: Symmetry has not been incorporated for single atom monomers.\n Turn off symmetry with \"SPACE_SYMMETRY = false\" and restart job\n");
    exit(0);
  }
  //Check symmetry for Spherical Tops monomer
  else if(fabs(MomentsOfInertia[2]-MomentsOfInertia[0]) < tolerance){

    Inertia.Transpose();
    LocalCoor = OrientateSphericalTop(LocalCoor,Inertia);

    if(SymmetryForSphericalTop(Monomers,LocalCoor,Inertia))
      EquivalencyFound = 1;
  }
  //Check symmetry for Symmetrical or Linear Tops monomer
  else if(fabs(MomentsOfInertia[0]-MomentsOfInertia[1])< 0.000001 ||
     fabs(MomentsOfInertia[1]- MomentsOfInertia[2])< 0.000001){
   
    //Set inertia axis assocated with non degenerate moment to the z-axis if not already.
    if(fabs(MomentsOfInertia[1]- MomentsOfInertia[2])< 0.000001){
      Vector x_new = z_axis;
      Vector z_new = x_axis;
      //Inertia.Print("Inertia matrix");
      Inertia.ReorderColumns(0,2);
      //Inertia.Print("Reordered matrix");
      //printf("Switch x and z axis\n");
      for(int i=0;i<Natoms;i++) {
	LocalCoor(i,0) = COM(i,0)*x_new[0]+COM(i,1)*x_new[1]+COM(i,2)*x_new[2];
	LocalCoor(i,2) = COM(i,0)*z_new[0]+COM(i,1)*z_new[1]+COM(i,2)*z_new[2];
	//printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
	//	Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
      }
    }

    //Checking special case of a linear top, rotated it the same say as a spherical top but use COM coordinates
    if(fabs(MomentsOfInertia[0]) < 0.000001){

      //Discarding inertia matrix since COM are used instead of inertia axis
      Inertia.Set_Iden();
      LocalCoor = OrientateLinearTop(COM,Inertia);
      //LocalCoor.Print("LocalCoor");
      //Inertia.Print("Rotation");

     if(SymmetryForLinearTop(Monomers,LocalCoor,Inertia)){
       //printf("match for monomer %i was found\n",Monomer_index);
       EquivalencyFound = 1;
      }
    }  
    else{
     Inertia.Transpose();
     //Inertia.Print("Inertia");
     LocalCoor = OrientateSymmetricalTop(LocalCoor, Inertia);
     //Inertia.Print("Inertia after rotated");
     Inertia.Transpose();

      //Checking Symmetry for symmetrical top
      if(SymmetryForSymmetricTop(Monomers,LocalCoor,Inertia)){
        //printf("match for monomer %i was found\n",Monomer_index);
        EquivalencyFound = 1;
      }
    }
  }
  //Check symmetry for asymmetrical tops
  else if(SymmetryForAsymmetricTop(Monomers,LocalCoor,Inertia)){
    EquivalencyFound = 1;
    //printf("match for monomer %i was found\n",Monomer_index);
  }
  return EquivalencyFound;
}


//determine which asymmetrical top or linear molecules are symmetrical equivalent
bool Monomer::SymmetryForAsymmetricTop(Monomer Monomers[],Matrix LocalCoor,Matrix Inertia){


  //second monomer
  Matrix Inertia2(3,3); 
  Inertia2.Set();
  for(int i=1;i<Monomer_index;i++){
    if(Natoms == Monomers[i].GetNumberOfAtoms() && Monomers[i].GetSymmetryFactor()){//check to see that the two monomers have the same number of atoms
      //printf("compare to monomer %i with a symfac of %i\n", Monomers[i].GetIndex(),Monomers[i].GetSymmetryFactor());
      Matrix COM2(Monomers[i].GetNumberOfAtoms(),3);
      for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	COM2(j,0) = Monomers[i].Atom_List[j].GetPosition(0)-Monomers[i].GetCenterOfMass(0);// x center of mass
	COM2(j,1) = Monomers[i].Atom_List[j].GetPosition(1)-Monomers[i].GetCenterOfMass(1);// y center of mass
	COM2(j,2) = Monomers[i].Atom_List[j].GetPosition(2)-Monomers[i].GetCenterOfMass(2);// y center of mass	
	//inertia tensor
	double mass = Monomers[i].Atom_List[j].GetAtomicMass();
	Inertia2(0,0) += mass*(COM2(j,1)*COM2(j,1)+COM2(j,2)*COM2(j,2));//m(y*y + z*z)
	Inertia2(1,1) += mass*(COM2(j,0)*COM2(j,0)+COM2(j,2)*COM2(j,2));//m(x*x + z*z);
	Inertia2(2,2) += mass*(COM2(j,0)*COM2(j,0)+COM2(j,1)*COM2(j,1));//m(x*x + y*y);
	Inertia2(0,1) -= mass*COM2(j,0)*COM2(j,1);//m*x*y;
	Inertia2(0,2) -= mass*COM2(j,0)*COM2(j,2);//m*x*z;
	Inertia2(1,2) -= mass*COM2(j,1)*COM2(j,2);//m*y*z;
      }
      Inertia2(1,0) = Inertia2(0,1);
      Inertia2(2,0) = Inertia2(0,2);
      Inertia2(2,1) = Inertia2(1,2);
      //diagonalize inertia matrix
      Vector MomentOfInertia = Inertia2.Diagonalize();

      //determines how close the atoms have to be
      double tolerance = Params::Parameters().GetSymmetryTolerance();
      //Special case of a linear top, switching x- and z-axis
      if(fabs(MomentOfInertia[0]) < tolerance)
        Inertia2.ReorderColumns(0,2);


      //use for the local axes
      Vector x_axis = Inertia2.GetColumnVector(0);
      Vector y_axis = Inertia2.GetColumnVector(1);
      Vector z_axis = Inertia2.GetColumnVector(2);
      
      Matrix LocalCoor2(Monomers[i].GetNumberOfAtoms(),3);
     //printf("LocalCoor2 for m%i\n",i);
      for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	LocalCoor2(j,0) = COM2(j,0)*x_axis[0]+COM2(j,1)*x_axis[1]+COM2(j,2)*x_axis[2];
	LocalCoor2(j,1) = COM2(j,0)*y_axis[0]+COM2(j,1)*y_axis[1]+COM2(j,2)*y_axis[2];
	LocalCoor2(j,2) = COM2(j,0)*z_axis[0]+COM2(j,1)*z_axis[1]+COM2(j,2)*z_axis[2];
	//printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
	//	Monomers[i].Atom_List[j].GetSymbol().c_str(), LocalCoor2(j,0),LocalCoor2(j,1),LocalCoor2(j,2));
      }


      //determine if identical
      //these are used to perform reflection
      for(int XReflex=1;XReflex >=-1;XReflex-=2){//reflection about YZ plane
	for(int YReflex=1;YReflex >=-1;YReflex-=2){//reflection about XZ plane
	  for(int ZReflex=1;ZReflex >=-1;ZReflex-=2){//reflection about XY plane
	    //printf("%i, %i, %i\n",XReflex,YReflex,ZReflex);
	    Matrix Reflex(3,3);
	    Reflex.Set();
	    Reflex(0,0) = XReflex;
	    Reflex(1,1) = YReflex;
	    Reflex(2,2) = ZReflex;

	    //Reflex.Print("Reflex");
	    Matrix OperatedLocal2 = LocalCoor2.Multiply(Reflex);
	    //printf("operated matrix is:\n");
	    //for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	    //  printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
	    //	      Monomers[i].Atom_List[j].GetSymbol().c_str(),OperatedLocal2(j,0),OperatedLocal2(j,1),OperatedLocal2(j,2));	
	    //}
	    int j=0;//counter for atoms of molecule one
	    bool nosym = 0;//one if molecule are not identical
	    while(j<Natoms && nosym == 0){//this is a loop over the atoms of monomer 1
	      bool match = 0;//this indicates whether there is a matching atom
	      int k = 0;
	      while(k<Natoms && match ==0){//looping over atoms of second monomer, the two monomers have the same number of atoms
		if(Atom_List[j].GetSymbol() == Monomers[i].Atom_List[k].GetSymbol()){
		  double x_diff = fabs(LocalCoor(j,0) - OperatedLocal2(k,0));
		  double y_diff = fabs(LocalCoor(j,1) - OperatedLocal2(k,1));
		  double z_diff = fabs(LocalCoor(j,2) - OperatedLocal2(k,2));
		  if(x_diff < tolerance && y_diff < tolerance && z_diff < tolerance)
		    match = 1;//identical atom found
		}
		k++;
	      }		
	      if(match == 0)
		nosym = 1;
	      j++;
	    }				
	    if(nosym == 0){
	      sym_fac = 0;
	      Monomers[i].IncreaseSymmetryFactor();

	      //Setting which monomer this monomer is symmetrical to
	      SymMon = Monomers[i].GetIndex();
	      
	      // char label[10];
	      //sprintf(label, "%d", Monomers[i].GetIndex());
	      //SymMon = "m";
	      //SymMon += label;
	      

	      //Inertia.Print("Inertia");
	      //Inertia2.Print("Inertia2");
	      //Reflex.Print("Reflex");


	      //Matrix to rotate this monomer into its symmetrical equivalent monomer
	      Matrix Rotation = Reflex.Multiply(Inertia,4);
	      Rotation = Inertia2.Multiply(Rotation);
	      //Rotation.Print("Rotation using Inertia");
	      
	      //setting rotation matrix
	      SetRotationMatrix(Rotation);

	      // if(Params::Parameters().DoForces()){
	      FindSymmetricalEquivalentAtoms(Monomers[i]);

	      //FindMapping(Monomers[i]);
	      //}
	      return 1;
	    }//end of if statement for nosym
	  }//end of reflection about XY plane
	}//end of reflection about XZ plane
      }//end of reflection about XY

    }
  }
 
  //Find symmetry fo the symmetrically unique monomers
  DeterminePointGroup(LocalCoor,Inertia,0);
 
  return 0;
}

//determine which spherical tops monomers are equivalent
bool Monomer::SymmetryForSphericalTop(Monomer Monomers[],Matrix LocalCoor,Matrix Rot){

  //All previous monomers will be rotated in an identical matter as this monomer
  //in the function Monomer::SymmetryCheck()

  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance

  Matrix LocalCoor2(Natoms,3);//LocalCoordinates for monomer checking against
  for(int i=1;i<Monomer_index;i++){//looping over all previous monomers
    if(Natoms == Monomers[i].GetNumberOfAtoms() && Monomers[i].GetSymmetryFactor()){//check to see that the two monomers have the same number of atoms

      //printf("Center of mass coordinates for monomer m%i:\n",i);

      //Finding center of mass coordinates
      //Not using Inertia Tenser because spherical tops do not have unique Principle Inertia Vectors due to the degeneracy of the Moment
      for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	LocalCoor2(j,0) = Monomers[i].Atom_List[j].GetPosition(0)-Monomers[i].GetCenterOfMass(0);// x center of mass
	LocalCoor2(j,1) = Monomers[i].Atom_List[j].GetPosition(1)-Monomers[i].GetCenterOfMass(1);// y center of mass
	LocalCoor2(j,2) = Monomers[i].Atom_List[j].GetPosition(2)-Monomers[i].GetCenterOfMass(2);// y center of mass
	//printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
	//	Monomers[i].Atom_List[j].GetSymbol().c_str(),LocalCoor2(j,0),LocalCoor2(j,1),LocalCoor2(j,2));

      }

      Matrix Rot2(3,3);
      Rot2.Set_Iden();
      //printf("Compare symmetry with m%i\n",i);
      LocalCoor2 = OrientateSphericalTop(LocalCoor2,Rot2);

      //determine if identical
      //these are used to perform reflection
      for(int XReflex=1;XReflex >=-1;XReflex-=2){//reflection about YZ plane
	for(int YReflex=1;YReflex >=-1;YReflex-=2){//reflection about XZ plane
	  for(int ZReflex=1;ZReflex >=-1;ZReflex-=2){//reflection about XY plane
	    //printf("%i, %i, %i\n",XReflex,YReflex,ZReflex);
	    Matrix Reflex(3,3);
	    Reflex.Set();
	    Reflex(0,0) = XReflex;
	    Reflex(1,1) = YReflex;
	    Reflex(2,2) = ZReflex;
	    Matrix OperatedLocal2= LocalCoor2.Multiply(Reflex);
	      /*printf("operated matrix is:\n");
		for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
		printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
		Monomers[i].Atom_List[j].GetSymbol().c_str(),OperatedLocal2(j,0),OperatedLocal2(j,1),OperatedLocal2(j,2));	
		}*/
	    int j=0;//counter for atoms of molecule one
	    bool nosym = 0;//one if molecule are not identical
	    while(j<Natoms && nosym == 0){//this is a loop over the atoms of monomer 1
	      bool match = 0;//this indicates whether there is a matching atom
	      int k = 0;
	      while(k<Natoms && match ==0){//looping over atoms of second monomer, the two monomers have the same number of atoms
		if(Atom_List[j].GetSymbol() == Monomers[i].Atom_List[k].GetSymbol()){
		  double x_diff = fabs(LocalCoor(j,0) - OperatedLocal2(k,0));
		  double y_diff = fabs(LocalCoor(j,1) - OperatedLocal2(k,1));
		  double z_diff = fabs(LocalCoor(j,2) - OperatedLocal2(k,2));
		  if(x_diff < tolerance && y_diff < tolerance && z_diff < tolerance)
		    match = 1;//identical atom found
		}
		k++;
	      }		
	      if(match == 0)
		nosym = 1;
	      j++;
	    }				
	    if(nosym == 0){
	      //printf("monomer %i is the identical to monomer %i\n", Monomer_index,Monomers[i].GetIndex());
	      sym_fac = 0;

	      //Rotation Matrix
              Matrix Rotation = Reflex.Multiply(Rot,2);
	      Rotation = Rot2.Multiply(Rotation,2);

	      //Rotation.Print("Rotation");
	      Monomers[i].IncreaseSymmetryFactor();

	      //setting rotation matrix
	      SetRotationMatrix(Rotation);
	      
	      //Setting which monomer this monomer is symmetrical to
	      SymMon = Monomers[i].GetIndex();
	       
	      // char label[10];
	      //sprintf(label, "%d", Monomers[i].GetIndex());
	      //SymMon = "m";
	      //SymMon += label;
	      
	      // if(Params::Parameters().DoForces()){
	      FindSymmetricalEquivalentAtoms(Monomers[i]);

	      //FindMapping(Monomers[i]);
	      
	      //}
	      return 1;
	    }//end of if statement for nosym
	  }//end of reflection about XY plane
	}//end of reflection about XZ plane
      }//end of reflection about XY 
    }//end if statement for check for symmetry factor
  }//end of loop over monomers

  //Find symmetry for the symmetrically unique monomers
  DeterminePointGroup(LocalCoor,Rot,3);
  
  return 0;

}

//determine which Linear tops monomers are equivalent
bool Monomer::SymmetryForLinearTop(Monomer Monomers[],Matrix LocalCoor,Matrix Rot){

  //All previous monomers will be rotated in an identical matter as this monomer
  //in the function Monomer::SymmetryCheck()

  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance

  Matrix LocalCoor2(Natoms,3);//LocalCoordinates for monomer checking against
  for(int i=1;i<Monomer_index;i++){//looping over all previous monomers
    if(Natoms == Monomers[i].GetNumberOfAtoms() && Monomers[i].GetSymmetryFactor()){//check to see that the two monomers have the same number of atoms

      //printf("Center of mass coordinates for monomer m%i:\n",i);

      //Finding center of mass coordinates
      //Not using Inertia Tenser because spherical tops do not have unique Principle Inertia Vectors due to the degeneracy of the Moment
      for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	LocalCoor2(j,0) = Monomers[i].Atom_List[j].GetPosition(0)-Monomers[i].GetCenterOfMass(0);// x center of mass
	LocalCoor2(j,1) = Monomers[i].Atom_List[j].GetPosition(1)-Monomers[i].GetCenterOfMass(1);// y center of mass
	LocalCoor2(j,2) = Monomers[i].Atom_List[j].GetPosition(2)-Monomers[i].GetCenterOfMass(2);// y center of mass
	//printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
	//	Monomers[i].Atom_List[j].GetSymbol().c_str(),LocalCoor2(j,0),LocalCoor2(j,1),LocalCoor2(j,2));

      }

      Matrix Rot2(3,3);
      Rot2.Set_Iden();
      //printf("Compare symmetry with m%i\n",i);
      LocalCoor2 = OrientateLinearTop(LocalCoor2,Rot2);

      //determine if identical
      //these are used to perform reflection
      for(int XReflex=1;XReflex >=-1;XReflex-=2){//reflection about YZ plane
	for(int YReflex=1;YReflex >=-1;YReflex-=2){//reflection about XZ plane
	  for(int ZReflex=1;ZReflex >=-1;ZReflex-=2){//reflection about XY plane
	    //printf("%i, %i, %i\n",XReflex,YReflex,ZReflex);
	    Matrix Reflex(3,3);
	    Reflex.Set();
	    Reflex(0,0) = XReflex;
	    Reflex(1,1) = YReflex;
	    Reflex(2,2) = ZReflex;
	    Matrix OperatedLocal2= LocalCoor2.Multiply(Reflex);
	      /*printf("operated matrix is:\n");
		for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
		printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
		Monomers[i].Atom_List[j].GetSymbol().c_str(),OperatedLocal2(j,0),OperatedLocal2(j,1),OperatedLocal2(j,2));	
		}*/
	    int j=0;//counter for atoms of molecule one
	    bool nosym = 0;//one if molecule are not identical
	    while(j<Natoms && nosym == 0){//this is a loop over the atoms of monomer 1
	      bool match = 0;//this indicates whether there is a matching atom
	      int k = 0;
	      while(k<Natoms && match ==0){//looping over atoms of second monomer, the two monomers have the same number of atoms
		if(Atom_List[j].GetSymbol() == Monomers[i].Atom_List[k].GetSymbol()){
		  double x_diff = fabs(LocalCoor(j,0) - OperatedLocal2(k,0));
		  double y_diff = fabs(LocalCoor(j,1) - OperatedLocal2(k,1));
		  double z_diff = fabs(LocalCoor(j,2) - OperatedLocal2(k,2));
		  if(x_diff < tolerance && y_diff < tolerance && z_diff < tolerance)
		    match = 1;//identical atom found
		}
		k++;
	      }		
	      if(match == 0)
		nosym = 1;
	      j++;
	    }				
	    if(nosym == 0){
	      //printf("monomer %i is the identical to monomer %i\n", Monomer_index,Monomers[i].GetIndex());
	      //Rot.Print("Rot");
	      //Reflex.Print("Reflex");
	      //Rot2.Print("Rot2");
	      
	      sym_fac = 0;

	      //Rotation Matrix
              Matrix Rotation = Reflex.Multiply(Rot,2);
	      Rotation = Rot2.Multiply(Rotation,2);

	      //Rotation.Print("Rotation");
	      Monomers[i].IncreaseSymmetryFactor();

	      //setting rotation matrix
	      SetRotationMatrix(Rotation);
	      
	      //Setting which monomer this monomer is symmetrical to
	      SymMon = Monomers[i].GetIndex();
	       
	      // char label[10];
	      //sprintf(label, "%d", Monomers[i].GetIndex());
	      //SymMon = "m";
	      //SymMon += label;
	      
	      // if(Params::Parameters().DoForces()){
	      FindSymmetricalEquivalentAtoms(Monomers[i]);

	      //FindMapping(Monomers[i]);
	      
	      //}
	      return 1;
	    }//end of if statement for nosym
	  }//end of reflection about XY plane
	}//end of reflection about XZ plane
      }//end of reflection about XY 
    }//end if statement for check for symmetry factor
  }//end of loop over monomers

  //Find symmetry for the symmetrically unique monomers
  DeterminePointGroup(LocalCoor,Rot,1);
  
  return 0;

}


//determine which symmetrical tops monomers are equivalent
bool Monomer::SymmetryForSymmetricTop(Monomer Monomers[],Matrix LocalCoor,Matrix Inertia){

  //All previous monomers will be rotated in an identical matter as this monomer
  //in the function Monomer::SymmetryCheck()

  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance

  Matrix Inertia2(3,3);//inertia tensor for second monomer
  Inertia2.Set();
  for(int i=1;i<Monomer_index;i++){//looping over all previous monomers
    if(Natoms == Monomers[i].GetNumberOfAtoms() && Monomers[i].GetSymmetryFactor()){//check to see that the two monomers have the same number of atoms
      //printf("compare to monomer %i with a symfac of %i\n", Monomers[i].GetIndex(),Monomers[i].GetSymmetryFactor());
      Matrix COM2(Monomers[i].GetNumberOfAtoms(),3);
      for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	COM2(j,0) = Monomers[i].Atom_List[j].GetPosition(0)-Monomers[i].GetCenterOfMass(0);// x center of mass
	COM2(j,1) = Monomers[i].Atom_List[j].GetPosition(1)-Monomers[i].GetCenterOfMass(1);// y center of mass
	COM2(j,2) = Monomers[i].Atom_List[j].GetPosition(2)-Monomers[i].GetCenterOfMass(2);// z center of mass	
	//inertia tensor
	double mass = Monomers[i].Atom_List[j].GetAtomicMass();
	Inertia2(0,0) += mass*(COM2(j,1)*COM2(j,1)+COM2(j,2)*COM2(j,2));//m(y*y + z*z)
	Inertia2(1,1) += mass*(COM2(j,0)*COM2(j,0)+COM2(j,2)*COM2(j,2));//m(x*x + z*z);
	Inertia2(2,2) += mass*(COM2(j,0)*COM2(j,0)+COM2(j,1)*COM2(j,1));//m(x*x + y*y);
	Inertia2(0,1) -= mass*COM2(j,0)*COM2(j,1);//m*x*y;
	Inertia2(0,2) -= mass*COM2(j,0)*COM2(j,2);//m*x*z;
	Inertia2(1,2) -= mass*COM2(j,1)*COM2(j,2);//m*y*z;
      }
      Inertia2(1,0) = Inertia2(0,1);
      Inertia2(2,0) = Inertia2(0,2);
      Inertia2(2,1) = Inertia2(1,2);
      //diagonalize inertia matrix
      Vector MomentsOfInertia2 = Inertia2.Diagonalize();
      //Inertia2.Print("Inertia2");
      
      bool skip_monomer = 1;//critiary to know if we should skip monomer if not a symmetrical top

      //Making sure that the monomer is a symmetrical top by seeing if two moments are degenerate.
      //The degeneracy of the first two moments are degenerate. 
      //It is possible that the second and third are the degenerate ones. This will be check next
      if(fabs(MomentsOfInertia2[0]-MomentsOfInertia2[1]) < tolerance)
	skip_monomer=0;

      //Checking to see if the second and third monomers are degenerate
      //Also making sure that the inertia axis assocated with the non-degenerate moment is along the z-axis
      if(fabs(MomentsOfInertia2[1]-MomentsOfInertia2[2]) < tolerance){
	skip_monomer=0; 
	//switching x and z
	Inertia2.ReorderColumns(0,2);
	
      }

      //check to make sure monomer is not spherical top. Only two of the moments should be degenerate
      if(fabs(MomentsOfInertia2[0]-MomentsOfInertia2[1]) < tolerance &&
	 fabs(MomentsOfInertia2[1]-MomentsOfInertia2[2]) < tolerance)
	skip_monomer=1;
      
      //only continue with this monomer is a symmetrical top
      if(!skip_monomer){
	Vector x_axis = Inertia2.GetColumnVector(0);
	Vector y_axis = Inertia2.GetColumnVector(1);
	Vector z_axis = Inertia2.GetColumnVector(2);
	Matrix LocalCoor2(Monomers[i].GetNumberOfAtoms(),3);
	for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
	  LocalCoor2(j,0) = COM2(j,0)*x_axis[0]+COM2(j,1)*x_axis[1]+COM2(j,2)*x_axis[2];
	  LocalCoor2(j,1) = COM2(j,0)*y_axis[0]+COM2(j,1)*y_axis[1]+COM2(j,2)*y_axis[2];
	  LocalCoor2(j,2) = COM2(j,0)*z_axis[0]+COM2(j,1)*z_axis[1]+COM2(j,2)*z_axis[2];
	}//end of loop over atoms

	//printf("Camparing to m%i\n",i);
	Inertia2.Transpose();
	LocalCoor2 = OrientateSymmetricalTop(LocalCoor2,Inertia2);
	Inertia2.Transpose();
	//determine if identical
	//these are used to perform reflection
	for(int XReflex=1;XReflex >=-1;XReflex-=2){//reflection about YZ plane
	  for(int YReflex=1;YReflex >=-1;YReflex-=2){//reflection about XZ plane
	    for(int ZReflex=1;ZReflex >=-1;ZReflex-=2){//reflection about XY plane
	      //printf("%i, %i, %i\n",XReflex,YReflex,ZReflex);
	      Matrix Reflex(3,3);
	      Reflex.Set();
	      Reflex(0,0) = XReflex;
	      Reflex(1,1) = YReflex;
	      Reflex(2,2) = ZReflex;
	      Matrix OperatedLocal2= LocalCoor2.Multiply(Reflex);
	      /*printf("operated matrix is:\n");
		for(int j=0; j<Monomers[i].GetNumberOfAtoms(); j++){
		printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
		Monomers[i].Atom_List[j].GetSymbol().c_str(),OperatedLocal2(j,0),OperatedLocal2(j,1),OperatedLocal2(j,2));	
		}*/
	      int j=0;//counter for atoms of molecule one
	      bool nosym = 0;//one if molecules are not identical
	      while(j<Natoms && nosym == 0){//this is a loop over the atoms of monomer 1
		bool match = 0;//this indicates whether there is a matching atom
		int k = 0;
		while(k<Natoms && match ==0){//looping over atoms of second monomer, the two monomers have the same number of atoms
		  if(Atom_List[j].GetSymbol() == Monomers[i].Atom_List[k].GetSymbol()){
		    double x_diff = fabs(LocalCoor(j,0) - OperatedLocal2(k,0));
		    double y_diff = fabs(LocalCoor(j,1) - OperatedLocal2(k,1));
		    double z_diff = fabs(LocalCoor(j,2) - OperatedLocal2(k,2));
		    if(x_diff < tolerance && y_diff < tolerance && z_diff < tolerance)
		      match = 1;//identical atom found
		  }
		  k++;
		}		
		if(match == 0)
		  nosym = 1;
		j++;
	      }				
	      if(nosym == 0){
		//printf("monomer %i is the identical to monomer %i\n", Monomer_index,Monomers[i].GetIndex());
		sym_fac = 0;
		Monomers[i].IncreaseSymmetryFactor();
		
		//Setting which monomer this monomer is symmetrical to
		SymMon = Monomers[i].GetIndex();
		
		// char label[10];
		//sprintf(label, "%d", Monomers[i].GetIndex());
		//SymMon = "m";
		//SymMon += label;

		//Inertia.Print("Inertia");
		//Inertia2.Print("Inertia2");
		//Reflex.Print("Reflex");

		//Matrix to rotate this monomer into its symmetrical equivalent monomer
		Matrix Rotation = Reflex.Multiply(Inertia,4);
		Rotation = Inertia2.Multiply(Rotation);
		//Rotation.Print("Rotation using Inertia");

		//setting rotation matrix
		SetRotationMatrix(Rotation);
		FindSymmetricalEquivalentAtoms(Monomers[i]);

		// if(Params::Parameters().DoForces()){
		//FindMapping(Monomers[i]);
		
		//}
		return 1;
	      }//end of if statement for nosym
	    }//end of reflection about XY plane
	  }//end of reflection about XZ plane
	}//end of reflection about XY
	
      }//end of if statement checking if monomers should be skipped
    }//end of if statement check monomer symmetry factor
  }//end of loop over monomers


  //Find symmetry for the symmetrically unique monomers
  DeterminePointGroup(LocalCoor,Inertia,2);
  return 0;
}



//Here Spherical Tops are orientated in such a way so that determining symmetry will be easy.
//Only the Td and Oh point groups are programmed. The other point groups which are Spherical Tops are rare.
Matrix Monomer::OrientateSphericalTop(Matrix LocalCoor, Matrix &Rotation){

  //Tolerance
  double tolerance = Params::Parameters().GetSymmetryTolerance();
  //Atom rotated into the Z axis
  int Zatom;
  //Atom with the shortest distance will be placed on the z-axis
  double low_dist = 9999;
  
  //Finding atom with the shortest distance to the origin but not at the origin. 
  //This atom is rotated into the z-axis
  for(int i=0;i<Natoms;i++){
    double dist = pow(LocalCoor(i,0),2)+pow(LocalCoor(i,1),2)+pow(LocalCoor(i,2),2);
    dist = sqrt(dist);
    //printf("i = %i dist = %f\n",i,dist); 
    if(dist < low_dist && dist > tolerance){
      Zatom = i;
      low_dist = dist;
    }
  }

  //printf("Zatom = %i",Zatom);
  
  //Rotate so the atom is in the XZ plane 
  Matrix Rot(3,3);
  double angle;
  if(fabs(LocalCoor(Zatom,1))<tolerance)
    angle= pi/2;
  else
    angle = atan(LocalCoor(Zatom,1)/LocalCoor(Zatom,2));
  Rot.Set_Iden();
  Rot(1,1)=cos(angle);
  Rot(1,2)=-sin(angle);
  Rot(2,1)=sin(angle);
  Rot(2,2)=cos(angle);
  //printf("Zatom = %i\n",Zatom);
  //printf("angle = %f\n",angle);
  //Rot.Print("Rot");
  LocalCoor.Transpose();
  LocalCoor = Rot.Multiply(LocalCoor);
  LocalCoor.Transpose();


  //Matrix that is the product of rotations that get the monomer from the COM coordinates into its local coordinates.
  Rotation = Rot.Multiply(Rotation);
  //printf("Coordinates after rotating atom into XZ plane\n");
  //for(int i=0;i<Natoms;i++) {
  // printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
  //	    Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  //}
  
  //Rotate so the atom is along the Z-axis.
  if(fabs(LocalCoor(Zatom,0)) < tolerance)
    angle = pi/2;
  else
    angle = atan(LocalCoor(Zatom,0)/LocalCoor(Zatom,2));
  Rot.Set_Iden();
  Rot(0,0)=cos(angle);
  Rot(0,2)=-sin(angle);
  Rot(2,0)=sin(angle);
  Rot(2,2)=cos(angle);   
  //printf("Zatom = %i\n",Zatom);
  //printf("angle = %f\n",angle);
  //Rot.Print("Rot");
  LocalCoor.Transpose();
  LocalCoor = Rot.Multiply(LocalCoor);
  LocalCoor.Transpose();
  

  //Matrix that is the product of rotations that get the monomer from the COM coordinates into its local coordinates.
  Rotation = Rot.Multiply(Rotation);  
  //printf("Coordinates after rotating atom into Z-axis\n");
  //for(int i=0;i<Natoms;i++) {
  //  printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
  //	    Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  //}
  

  //Determine symmetrically equivalent atom to the "Zatom" 
  //This atom wil be rotated into the XZ plane.
  int XZatom;
    //first pick atom
    for(int i=0;i<Natoms;i++){
      bool skip_atom = 0;//bool to skip atom
      if(i==Zatom){// skip the Zatom
	skip_atom=1;
	//printf("%i==%i\n",i,Zatom);
      }
      if(GetAtom(i).GetSymbol()!=GetAtom(Zatom).GetSymbol()){// the atom must be the same as Zatom
	skip_atom=1;
	//printf("%i %s!=%s\n",
	//       i,GetAtom(i).GetSymbol().c_str(),GetAtom(Zatom).GetSymbol().c_str());
      }
      double ZCoord = pow(LocalCoor(i,0),2) + pow(LocalCoor(i,1),2);//skip any atom along the z axis
      ZCoord = sqrt(ZCoord);
 
      if(fabs(ZCoord) < tolerance){//skip atoms along the z axis
	skip_atom=1;
	//printf("%i ZCoord = %f\n",i,ZCoord);
      }

      double dist = pow(LocalCoor(i,0),2) + pow(LocalCoor(i,1),2) + pow(LocalCoor(i,2),2);//distance must be the same as the Zatom
      dist = sqrt(dist);
      if(fabs(dist-low_dist) > tolerance){
	skip_atom=1;
	//printf("%i dist = %f low_dist = %f\n",
	//       i,dist,low_dist);
      }
      if(!skip_atom){
	XZatom = i;
	break;
      }
    }
      
    //printf("XZatom = %i\n",XZatom);
    
    //Rotate "XZatom" into XZ plane
    if(fabs(LocalCoor(XZatom,1)) < tolerance)
      angle = pi/2;
    else
      angle = atan(LocalCoor(XZatom,1)/LocalCoor(XZatom,0));
    Rot.Set_Iden();
    Rot(0,0)=cos(angle);
    Rot(0,1)=sin(angle);
    Rot(1,0)=-sin(angle);
    Rot(1,1)=cos(angle);  
    //printf("XZatom = %i\n",XZatom);
    //printf("angle = %f\n",angle);
    //Rot.Print("Rot");
    LocalCoor.Transpose();
    LocalCoor = Rot.Multiply(LocalCoor);
    LocalCoor.Transpose();

    //Matrix that is the product of rotations that get the monomer from the COM coordinates into its local coordinates.
    Rotation = Rot.Multiply(Rotation);
    //printf("Coordinates after rotating atom into XZ plane\n");
    //for(int i=0;i<Natoms;i++) {
    //  printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	      Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
    //
    
    //Rotation.Print("Full Rotation");
    return LocalCoor;
    
}


//Linear tops will be orientated in such a way so that determining symmetry will be easy.
Matrix Monomer::OrientateLinearTop(Matrix LocalCoor, Matrix &Rotation){

  //Tolerance
  double tolerance = Params::Parameters().GetSymmetryTolerance();
  //Atom to be rotated into the Z axis
  int atom;
  

  //Find the atom with the shortest distance from the origin but not at the origin. This atom will be rotated into the z-axis
  for(int i=0;i<Natoms;i++){
    double dist = pow(LocalCoor(i,0),2)+pow(LocalCoor(i,1),2)+pow(LocalCoor(i,2),2);
    dist = sqrt(dist);
    //printf("i = %i dist = %f\n",i,dist); 
    if(dist > tolerance){
    atom = i;
    break;
    }
  }

  //printf("atom = %i",atom);
  
  //Rotate so atom "atom" is in the XZ plane 
  Matrix Rot(3,3);
  double angle;
  //if(fabs(LocalCoor(atom,1))<tolerance)// Watit Change here
  if(fabs(LocalCoor(atom,2))<tolerance)
    angle= pi/2;
  else
    angle = atan(LocalCoor(atom,1)/LocalCoor(atom,2));
  Rot.Set_Iden();
  Rot(1,1)=cos(angle);
  Rot(1,2)=-sin(angle);
  Rot(2,1)=sin(angle);
  Rot(2,2)=cos(angle);
  //printf("Zatom = %i\n",Zatom);
  //printf("angle = %f\n",angle);
  //Rot.Print("Rot");
  LocalCoor.Transpose();
  LocalCoor = Rot.Multiply(LocalCoor);
  LocalCoor.Transpose();


  //Matrix that is the product of rotations that get the monomer from the COM coordinates into its local coordinates.
  Rotation = Rot.Multiply(Rotation);
  //printf("Coordinates after rotating atom into XZ plane\n");
  //for(int i=0;i<Natoms;i++) {
  // printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
  //	    Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  //}
  
  //Rotate so the atom "atom" is in the Z-axis.
  //if(fabs(LocalCoor(atom,1))<tolerance)// Watit Change here
  if(fabs(LocalCoor(atom,2))<tolerance)
    angle = pi/2;
  else
    angle = atan(LocalCoor(atom,0)/LocalCoor(atom,2));
  Rot.Set_Iden();
  Rot(0,0)=cos(angle);
  Rot(0,2)=-sin(angle);
  Rot(2,0)=sin(angle);
  Rot(2,2)=cos(angle);   
  //printf("atom = %i\n",atom);
  //printf("angle = %f\n",angle);
  //Rot.Print("Rot");
  LocalCoor.Transpose();
  LocalCoor = Rot.Multiply(LocalCoor);
  LocalCoor.Transpose();
  

  //Matrix that is the product of rotations that get the monomer from the COM coordinates into its local coordinates.
  Rotation = Rot.Multiply(Rotation);
  
  return LocalCoor;
    
}


//Orient symmetrical top so that atom with the smallest XY projection is in the XZ plane.
//This will orientate the symmetrical top in a way that determining symmetry will be easy. 
Matrix Monomer::OrientateSymmetricalTop(Matrix LocalCoor, Matrix &PrevRot){

  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance
 
  //LocalCoor.Print("LocalCoor in Orientate function");


  //Finding atom with the lowest projection into the XY plane. 
  //If the right inertia axis is aligned along the z-axis, it will always be the same atom
  //or an symmetrically equivalent atom
  int low_atom;//atom with the lowest projection into the XY plane
  double low_proj = 1000;//lowest projection into the XY plane
  for(int i=0;i<Natoms;i++){
    bool skip_atom = 0;//critera for skipping atom
    double XY_proj = pow(LocalCoor(i,0),2)+pow(LocalCoor(i,1),2);
    XY_proj = sqrt(XY_proj);
    //printf("XY_proj = %f\n",XY_proj);
    if(XY_proj<tolerance){
      skip_atom = 1;
      //printf("skip atom %i\n",i);
    }
    if(XY_proj<low_proj && !skip_atom){
      //printf("counting %i\n",i);
      low_atom = i;
      low_proj = XY_proj;
    }
  }	
  //printf("low_atom = %i;low_proj = %f\n",
  //     low_atom, low_proj);
  
  //rotating molecule with shortest projection into the XZ plane
  Matrix Rot(3,3);
  double angle;
  if(fabs(LocalCoor(low_atom,0))<tolerance)
    angle= pi/2;
  else
    angle = atan(LocalCoor(low_atom,1)/LocalCoor(low_atom,0));
  Rot.Set_Iden();
  //printf("Look_coor = %f\n",LocalCoor(low_atom,1));
  //printf("angle = %f\n", angle*180/pi);
  Rot(0,0)=cos(angle);
  Rot(0,1)=sin(angle);
  Rot(1,0)=-sin(angle);
  Rot(1,1)=cos(angle);
  LocalCoor.Transpose();
  LocalCoor = Rot.Multiply(LocalCoor);
  LocalCoor.Transpose();
  
  //printf("Coordinates after rotating atom into XZ plane\n");
  //for(int i=0;i<Natoms;i++) {
  // printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
  //	    Atom_List[i].GetSymbol().c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  //}  
  
  //This is required for determining the Rotation Matrix
  // PrevRot.Print("PrevRot");
  //Rot.Print("Rot");
  PrevRot = Rot.Multiply(PrevRot);
  //PrevRot.Print("PrevRot * Rot");

  return LocalCoor;

}

//Determines which atoms are the symmetrical equivalent atom on the symmetrically unique mononers
void Monomer::FindSymmetricalEquivalentAtoms(Monomer& UniqueMon){

  //Get Center of Mass Coordinates and Apply rotation matrix
  Matrix COM(3,Natoms);
  COM.Set();

  for(int i=0;i<Natoms;i++) {
    COM(0,i) = Atom_List[i].GetPosition(0)-GetCenterOfMass(0);// x center of mass
    COM(1,i) = Atom_List[i].GetPosition(1)-GetCenterOfMass(1);// y center of mass
    COM(2,i) = Atom_List[i].GetPosition(2)-GetCenterOfMass(2);// y center of mass
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //		       GetAtom(i).GetSymbol().c_str(), COM(0,i),COM(1,i),COM(2,i));
  }

  //Get Center of Mass Coordinates for symmetrically unique monomer
  Matrix UniqueCOM(3,Natoms);
  UniqueCOM.Set();
  for(int i=0;i<Natoms;i++) {
    UniqueCOM(0,i) = UniqueMon.Atom_List[i].GetPosition(0)-UniqueMon.GetCenterOfMass(0);// x center of mass
    UniqueCOM(1,i) = UniqueMon.Atom_List[i].GetPosition(1)-UniqueMon.GetCenterOfMass(1);// y center of mass
    UniqueCOM(2,i) = UniqueMon.Atom_List[i].GetPosition(2)-UniqueMon.GetCenterOfMass(2);// y center of mass
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //		       UniqueMon.GetAtom(i).GetSymbol().c_str(), UniqueCOM(0,i),UniqueCOM(1,i),UniqueCOM(2,i));
  }
  
  //Apply rotation operator to non symmetrical unique monomer is it is in the same orienation
  //as the symmetrical unique monomer
  Matrix Rot = GetRotationMatrix();
  Matrix AfterT = Rot.Multiply(COM);
  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be

  //finding the equivalent atoms
  int i = 0;//index for atoms on this monomer
  while(i < Natoms){
    int MatchFound = 0;//tells you if the symmetrical equivalent atom was found
    int j=0;//index for atoms on the UniqueMon
    while(j<Natoms && MatchFound==0){
      	double x_diff = fabs(AfterT(0,i)-UniqueCOM(0,j));
	double y_diff = fabs(AfterT(1,i)-UniqueCOM(1,j));
	double z_diff = fabs(AfterT(2,i)-UniqueCOM(2,j));
	if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance){
	  MatchFound = 1;

	  //setting iatom to be symmerically equivalent to jatom or that atom jatom is equivalent to in the asymmetrical unit 
	  GetAtom(i).SetSymmetricalAtom(UniqueMon.GetAtom(j).GetSymmetricalAtom());
	  //GetAtom(i).SetSymmetricalAtom(UniqueMon.GetAtom(j).GetGlobalIndex());

	  /*
	    Setting atom rotation matrix as the rotatation operator to the atom in the asymmetrical unit cell from
	    the monomer rotation matrix and the Rotation matrix for jatom to atom in the asymmetrical unit cell.
	    Note that is rotation matrix may be corrected later by Cluster::CorrectingMonomerRotations().
	  */
	  Matrix Rot = UniqueMon.GetAtom(j).GetRotationMatrix().Multiply(Symmetry_Rotation);
	  GetAtom(i).SetRotationMatrix(Rot);

	  //Setting the equivalency between atoms between the monomers if rotated by Symmetry_Rotation
	  Sym_Atom[i] = j;
	}//end of if statement
	j++;
    }//end of loop over j
    i++;

    if(MatchFound ==0){
      printf("Error::Monomer::FindSymmetricalEquivalentAtoms:: Atom %i on m%i has no equivalent atom on m%i\n",
	     GetAtom(i).GetAtomIndex(),GetIndex(),UniqueMon.GetIndex());
      exit(1);
    }
  }//end of loop over i

}

//This function determines the point group symmetry of symmetrical unique monomers.
//The point group operations which atoms in the monomer are symmetrical equivalent and  
//in the crystal's asymmetric unit.
//Not all point groups are programmed in. 
void Monomer::DeterminePointGroup(Matrix LocalCoor,Matrix RotToLocal,int molecule_type){

  //if molecule_type = 0, molecule is an asymetrcial top
  //                 = 1, molecule is linear top
  //                 = 2, molecule is a symmetrical top             
  //                 = 3, molecule is a spherical top  
  

  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be

  string PointGroup = "";
  
  //operators

  //inversion
  Matrix inverse(3,true);
  inverse(0,0) = -1.0;
  inverse(1,1) = -1.0;
  inverse(2,2) = -1.0;
  //inverse.Print("inverse");
  
  //reflection in the yz plane
  Matrix reflect_yz(3,true);
  reflect_yz(0,0) = -1.0;
  //reflect_yz.Print("reflect_yz"); 
  
  //reflection in the xz plane
  Matrix reflect_xz(3,true);
  reflect_xz(1,1) = -1.0;
  //reflect_xz.Print("reflect_xz");
  
  //reflection in the yz plane
  Matrix reflect_xy(3,true);
  reflect_xy(2,2) = -1.0;
  //reflect_xy.Print("reflect_xy");
  
  //C2 rotation about the x-axis
  Matrix C2_x(3,true);
  C2_x(1,1) = -1.0;
  C2_x(2,2) = -1.0;
  //C2_x.Print("C2_x");
  
  //C2 rotation about the y-axis
  Matrix C2_y(3,true);
  C2_y(0,0) = -1.0;
  C2_y(2,2) = -1.0;
  //C2_y.Print("C2_y");
  
  //C2 rotation about the z-axis
  Matrix C2_z(3,true);
  C2_z(0,0) = -1.0;
  C2_z(1,1) = -1.0;
  //C2_z.Print("C2_z");

  //the rest of the rotations. Only goes up to C8. The molecule will be oriented so this rotation is about the z-axis
  
  //C3 rotation
  double angle = 2*pi/3;
  Matrix C3(3,true);
  C3(0,0) = cos(angle);
  C3(1,1) = cos(angle);
  C3(0,1) = -sin(angle);
  C3(1,0) = sin(angle);
  //C3.Print("C3");
  
  
  //C4 rotation
  angle = 2*pi/4;
  Matrix C4(3,true);
  C4(0,0) = cos(angle);
  C4(1,1) = cos(angle);
  C4(0,1) = -sin(angle);
  C4(1,0) = sin(angle);
  //C4.Print("C4");  
  
  //C5 rotation
  angle = 2*pi/5;
  Matrix C5(3,true);
  C5(0,0) = cos(angle);
  C5(1,1) = cos(angle);
  C5(0,1) = -sin(angle);
  C5(1,0) = sin(angle);
  //C5.Print("C5");
  
  //C6 rotation
  angle = 2*pi/6;
  Matrix C6(3,true);
  C6(0,0) = cos(angle);
  C6(1,1) = cos(angle);
  C6(0,1) = -sin(angle);
  C6(1,0) = sin(angle);
  //C6.Print("C6");

  //C7 rotation
  angle = 2*pi/7;
  Matrix C7(3,true);
  C7(0,0) = cos(angle);
  C7(1,1) = cos(angle);
  C7(0,1) = -sin(angle);
  C7(1,0) = sin(angle);
  //C7.Print("C6");

  //C8 rotation
  angle = 2*pi/8;
  Matrix C8(3,true);
  C8(0,0) = cos(angle);
  C8(1,1) = cos(angle);
  C8(0,1) = -sin(angle);
  C8(1,0) = sin(angle);
  //C8.Print("C8");
  
  //improper rotations, no S1 or S2 because those are the same as an reflection and inversion respectively
  Matrix S3 = reflect_xy.Multiply(C3);
  Matrix S4 = reflect_xy.Multiply(C4);
  Matrix S5 = reflect_xy.Multiply(C5);
  Matrix S6 = reflect_xy.Multiply(C6);
  Matrix S7 = reflect_xy.Multiply(C7);
  Matrix S8 = reflect_xy.Multiply(C8);
  
  LocalCoor.Transpose();
  RotToLocal.Transpose();

  //Assymetrical top point group
  //possible point groups: C1 Ci C2 Cs D2 C2v C2h D2h 
  if(molecule_type == 0){
    if ( Params::Parameters().PrintLevel() > 0)
      printf("\n molecule is an asymmetrical top\n");  

    //test inversion
    bool CanInverse = 0;
    Matrix OperatedCoor = inverse.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor))
      CanInverse = 1;
    
    //test reflection
    int NumberOfReflects = 0;
    OperatedCoor = reflect_xy.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor))
      NumberOfReflects += 1;
    OperatedCoor = reflect_xz.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor))  
      NumberOfReflects +=  1;
    OperatedCoor = reflect_yz.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor))
      NumberOfReflects += 1;
    
    //test C2 rotation
    int NumberOfC2 = 0;
    OperatedCoor = C2_x.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      NumberOfC2 += 1; 
    OperatedCoor = C2_y.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      NumberOfC2 += 1;
    OperatedCoor = C2_z.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      NumberOfC2 += 1;
    if ( Params::Parameters().PrintLevel() > 0)
      printf("Inverse = %i #reflects = %i #C2 = %i\n",CanInverse,NumberOfReflects,NumberOfC2);

    //determining point group
    if(CanInverse && NumberOfReflects == 3 && NumberOfC2 == 3)
      PointGroup = "D2h";
    else if(!CanInverse && NumberOfReflects == 0 && NumberOfC2 == 3)
      PointGroup = "D2";
    else if(!CanInverse && NumberOfReflects == 2 && NumberOfC2 == 1)
      PointGroup = "C2v";
    else if(CanInverse && NumberOfReflects == 1 && NumberOfC2 == 1)
      PointGroup = "C2h";
    else if(!CanInverse && NumberOfReflects == 0 && NumberOfC2 == 1)
      PointGroup = "C2";
    else if(CanInverse && NumberOfReflects == 0 && NumberOfC2 == 0)
      PointGroup = "Ci";
    else if(!CanInverse && NumberOfReflects == 1 && NumberOfC2 == 0)
      PointGroup = "Cs";
    else if(!CanInverse && NumberOfReflects == 0 && NumberOfC2 == 0)
      PointGroup = "C1";
    else {
      printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for asymmetrical top with %i inversion centers, %i reflection planes, and %i C2 axis\n",CanInverse,NumberOfReflects,NumberOfC2);
      exit(0);
    }

  }//Linear top point group
  else if(molecule_type == 1){
  //possible point group: C(inf)v or D(inf)h
    if ( Params::Parameters().PrintLevel() > 0)
      printf("\n molecule is an linear top\n");

    //Is there a reflection in the xy plane
    bool ReflexExist;
   
    Matrix OperatedCoor = reflect_xy.Multiply(LocalCoor);    
    //LocalCoor.Print("LocalCoor");
    //OperatedCoor.Print("OperatedCoor");

    ReflexExist = MatchCoordinates(LocalCoor,OperatedCoor);
  
    //Defining point group
    if(ReflexExist){
      PointGroup = "D(inf)h";
    }
    else{
      PointGroup = "C(inf)v";
    }
  
  }//Symmetrical top point groups
  else if (molecule_type == 2){
  //possible point groups: Cnv Cnh D2d Dnd Cn Dn S(2*n) where n > 2
    if ( Params::Parameters().PrintLevel() > 0)
      printf("\n molecule is a symmetrical top\n");
    
    //test inversion
    //bool CanInverse = 0;
    //Matrix OperatedCoor = inverse.Multiply(LocalCoor); 
    //if(MatchCoordinates(LocalCoor,OperatedCoor))
    //  CanInverse = 1;
 

    //test C2 rotation, only checking the ones perpendicular to the z-axis
    int NumberOfC2 = 0;
    Matrix OperatedCoor = C2_x.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      NumberOfC2 += 1; 
    OperatedCoor = C2_y.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      NumberOfC2 += 1;

    //test for C2_z, axis should be aligned along the z-axis
    bool CanC2_z = 0;
    OperatedCoor = C2_z.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      CanC2_z = 1;

    /*
    LocalCoor.Transpose();
    OperatedCoor.Transpose();
    LocalCoor.Print("LocalCoor");
    OperatedCoor.Print("OperatedCoor");
    C2_z.Print("C2_z");
    LocalCoor.Transpose();
    OperatedCoor.Transpose();
    */


    //test C3 rotation, axis should be aligned along the z-axis
    bool CanC3 = 0;
    OperatedCoor = C3.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      CanC3 = 1; 

    //test C4 rotation, axis should be aligned along the z-axis
    bool CanC4 = 0;
    OperatedCoor = C4.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      CanC4 = 1; 

    //test C5 rotation, axis should be aligned along the z-axis
    bool CanC5 = 0;
    OperatedCoor = C5.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      CanC5 = 1; 

    //test C6 rotation, axis should be aligned along the z-axis
    bool CanC6 = 0;
    OperatedCoor = C6.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor)) 
      CanC6 = 1; 

   
    //test reflection
    bool CanReflectXY = 0; 
    OperatedCoor = reflect_xy.Multiply(LocalCoor);
    if(MatchCoordinates(LocalCoor,OperatedCoor))
      CanReflectXY = 1;   

 
    int OtherReflects = 0;
    OperatedCoor = reflect_xz.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor))  
      OtherReflects +=  1;
    OperatedCoor = reflect_yz.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor))
      OtherReflects += 1;  

 
    
    //improper rotations
    bool improper_axis = 0;
    bool CanS4 = 0;
    OperatedCoor = S4.Multiply(LocalCoor);  
    if(MatchCoordinates(LocalCoor,OperatedCoor)) {
      CanS4 = 1;
      improper_axis = 1;
    }
    bool CanS5 = 0;
    OperatedCoor = S5.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor)) {
      CanS5 = 1;
      improper_axis = 1;
    }

    bool CanS6 = 0;
    OperatedCoor = S6.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor)) {
      CanS6 = 1;
      improper_axis = 1;
    }
      
    bool CanS7 = 0;
    OperatedCoor = S7.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor)) {
      CanS7 = 1;
      improper_axis = 1;
    }  
    bool CanS8 = 0;
    OperatedCoor = S8.Multiply(LocalCoor); 
    if(MatchCoordinates(LocalCoor,OperatedCoor)) {
      CanS8 = 1;
      improper_axis = 1;
    }  
    if ( Params::Parameters().PrintLevel() > 0){
      printf("CanReflectXY = %i NumberOfC2 = %i  OtherReflects = %i\n",CanReflectXY,NumberOfC2, OtherReflects);
      printf("CanS4 = %i CanS6 = %i \n",CanS4, CanS6);
    }

    if(CanReflectXY){//Dnh or Cnh groups
      if(NumberOfC2 > 0){//Dnh groups
	if(CanC6) PointGroup = "D6h";
	else if(CanC5) PointGroup = "D5h";
	else if(CanC4) PointGroup = "D4h";
	else if(CanC3) PointGroup = "D3h";
	else{
	  printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for symmetrical top with %i hortizontal C2 axis, a horizontal reflection, and no C(3-6) axis\n",
                 NumberOfC2);
	  exit(0);
	}
      }else{//Cnh groups
       	if(CanC6) PointGroup = "C6h";
	else if(CanC5) PointGroup = "C5h";
	else if(CanC4) PointGroup = "C4h";
	else if(CanC3) PointGroup = "C3h";
	else{
	  printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for symmetrical top with no C2 axis, a horizontal reflection, and no C(3-6) axis\n");
	  exit(0);
	}
      }	
    }
    else if(OtherReflects){//Dnd or CnV
      if(improper_axis){//Dnd
	if(CanC4) PointGroup = "D4d";
	else if(CanC3) PointGroup = "D3d";
	else if(CanC2_z) PointGroup = "D2d";
	else{//Dnd
	  printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for symmetrical top with a vertical C2 axis, an improper axis, and no C(2-4) axis\n");
	  exit(0);
	}
      }
      else {//Cnv
	if(CanC6) PointGroup = "C6v";
	else if(CanC5) PointGroup = "C5v";
	else if(CanC4) PointGroup = "C4v";
	else if(CanC3) PointGroup = "C3v";
	else{
	  printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for symmetrical top with a vertical C2 axis and no C(3-6) axis\n");
	  exit(0);
	}
      }
    }
    else{//Cn and S2n point groups
      if(CanS8) PointGroup = "S8";
      else if(CanS6) PointGroup = "S6";
      else if(CanS4) PointGroup = "S4";
      else if(CanC6) PointGroup = "C6";
      else if(CanC5) PointGroup = "C5";
      else if(CanC4) PointGroup = "C4";
      else if(CanC3) PointGroup = "C3";
      else{
	printf("Error::Monomer::DeterminePointGroup() : Unable to determine Point Group for symmetrical top with  C(3-6) or S(4-8) axis\n");
	exit(0);
      }	  
      if ( Params::Parameters().PrintLevel() > 0)
        printf("Point group = %s\n",PointGroup.c_str());      
    }
  }//spherical top
  else if(molecule_type == 3){
  //possible point group Td or Oh
    printf("\n molecule is an spherical top\n");
    printf("Currently unable to determine symmetry for a spherical top\n");
    exit(0);

  }
  else{
    printf("ERROR::Monomer::DeterminePointGroup() molecule type not recognized molecule_type = %i\n",molecule_type);
    exit(0);
 }
  if ( Params::Parameters().PrintLevel() > 0)
    printf("Point group is %s\n",PointGroup.c_str());
  
  vector<Matrix> Ops;

  //Get operator for point group
  //Not all point groups are programmed in. Add as needed.
  if(PointGroup == "C1"){
    //No useful operations in point group 
    return;
  }else if(PointGroup == "D(inf)h"){
    //not including any of the C(inf) or S(inf) operator. 
    Ops.push_back(inverse);
    Ops.push_back(C2_y);
    Ops.push_back(C2_x);
    Ops.push_back(reflect_xy);
  }else if(PointGroup == "Ci"){
    Ops.push_back(inverse);
  }else if(PointGroup == "Cs"){
    //don't know the plane of the reflection so including all of them
    Ops.push_back(reflect_yz);
    Ops.push_back(reflect_xz);
    Ops.push_back(reflect_xy);

  }else if(PointGroup == "C2v"){
    //don't know the rotation axis so including all of them
    Ops.push_back(C2_x);
    Ops.push_back(C2_y);
    Ops.push_back(C2_z);

    //don't know the plane of the reflection so including all of them
    Ops.push_back(reflect_yz);
    Ops.push_back(reflect_xz);
    Ops.push_back(reflect_xy);

    //yoni
    //C3v tested with phase I NH3, symmetry code is incompatible with this crystal due to inability to shift lattice parameters and keep same 
    //cartesian coordinates of atoms in an assymetrical cell without breaking symmetry.
    //C3v therefore has not be reliably tested


     }else if(PointGroup == "C3v"){
     
     //C3 and C3_2 rotation
    Ops.push_back(C3);
    Ops.push_back(C3.Multiply(C3));

    //all the reflection
    Ops.push_back(reflect_xz);
    Ops.push_back(reflect_xz.Multiply(C3));
    Ops.push_back(reflect_xz.Multiply(C3).Multiply(C3));


  }else{
    printf("ERROR::Monomer::DeterminePointGroup() Point group \"%s\" not programmed in yet. Include \"SPACE_SYMMETRY = false\" in input or programmed in point group\n",
	   PointGroup.c_str());
    exit(0);
  }

  //use operations and match atoms
  for(int i=0;i<Ops.size();i++){
    Matrix OperatedCoor = Ops[i].Multiply(LocalCoor);

    for(int j = 0;j<Natoms;j++){
      int k=0;
      bool AtomMatched = false;
      for(int k = 0; k<Natoms;k++){
	//double x_diff = fabs(LocalCoor(0,j)-OperatedCoor(0,k));
	//double y_diff = fabs(LocalCoor(1,j)-OperatedCoor(1,k));
	//double z_diff = fabs(LocalCoor(2,j)-OperatedCoor(2,k)); 
	double x_diff = fabs(LocalCoor(0,k)-OperatedCoor(0,j));
	double y_diff = fabs(LocalCoor(1,k)-OperatedCoor(1,j));
	double z_diff = fabs(LocalCoor(2,k)-OperatedCoor(2,j)); 

	bool SameAtomType = false;
	if(GetSymbol(j) == GetSymbol(k)) {
	  SameAtomType = true;
	  
	}
	
	if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance && SameAtomType){

	  //printf("Atom %i is sym to Atom %i\n",Atom_List[j].GetGlobalIndex(),Atom_List[k].GetGlobalIndex());

	  if(Atom_List[j].GetSymmetricalAtom() >  Atom_List[k].GetGlobalIndex()){
	    Atom_List[j].SetSymmetricalAtom(Atom_List[k].GetGlobalIndex());
	    Matrix Rot = Ops[i].Multiply(RotToLocal);
	    RotToLocal.Transpose();
	    Rot = RotToLocal.Multiply(Rot);
	    RotToLocal.Transpose();
	    Atom_List[j].SetRotationMatrix(Rot);
	  } 
	}      
      }//end of loop over k
    }//end of loop over j
  }//end of loop over i

  //RotToLocal.Print("RotToLocal");

  //Setting the # of symmetrical unique atoms in the monomer and important information related to the point group to these atoms.
  unique_atoms = 0;
  for(int i=0; i<Natoms;i++){
    if ( Params::Parameters().PrintLevel() > 0)
      printf("%i sym to %i\n",Atom_List[i].GetGlobalIndex(),Atom_List[i].GetSymmetricalAtom());
    //Atom_List[i].GetRotationMatrix().Print("Rot");

    if(Atom_List[i].InAsymmetricUnit()){
      unique_atoms++;
      

      //if(PointGroup == "D(inf)h" || PointGroup == "Ci"){

      Vector COM(3);
      COM[0] = Atom_List[i].GetPosition(0) - GetCenterOfMass(0);
      COM[1] = Atom_List[i].GetPosition(1) - GetCenterOfMass(1);    
      COM[2] = Atom_List[i].GetPosition(2) - GetCenterOfMass(2);   
      
      //printf("Atom %i COM = %f %f %f\n",Atom_List[i].GetGlobalIndex(),COM[0],COM[1],COM[2]);
      
      //Freezing atom if it is at the COM(center of mass) since that may break Point Group symmetry
      //Here we are determining if certain atoms or coordinates need to be fixed or linked in order to preserve
      //point (and space) group symmetry
      //Further freezing and linking of atoms occurs in Monomer::SetFractTranslationVector().
      double R = sqrt(LocalCoor(0,i)*LocalCoor(0,i)+LocalCoor(1,i)*LocalCoor(1,i)+LocalCoor(2,i)*LocalCoor(2,i));
      //Freezing atom at COM coordinate origin
      if(R < tolerance){
	Atom_List[i].SetAtomToFrozen(1);
	
      }else{



        //Freezing x coordinate if it is zero
        if(fabs(COM[0]) < tolerance) {
          //printf("Atom %i x COM coordinate is zero, freezing\n",Atom_List[i].GetGlobalIndex());
          Atom_List[i].SetFreezeX(1);
        }
	//If the absolute value of the x and z COM coordiantes for an atom match, they are linked by symmetry
	//Then checking if the sign is the same.	
	else if(fabs((fabs(COM[2]) - fabs(COM[0]))) < tolerance){
	  //Atom_List[i].LockZ(1);
	  Atom_List[i].LockX(2);
	  if(fabs(COM[2] - COM[0]) > tolerance) {
	    //Atom_List[i].SetZGradSign(1);
	    Atom_List[i].SetXGradSign(1);
	  }
	}
	//Assuming that if the absolute value of the x and y COM coordiantes for an atom match, that they are linked by symmetry
	//Then checking if the sign is the same 
	else if(fabs((fabs(COM[1]) - fabs(COM[0]))) < tolerance){
	  //Atom_List[i].LockY(1);
	  Atom_List[i].LockX(1);
	  if(fabs(COM[1] - COM[0]) > tolerance){
	    //Atom_List[i].SetYGradSign(1);
	    Atom_List[i].SetXGradSign(1);

	  }
	}


        //Freezing y coordinate if it is zero
        if(fabs(COM[1]) < tolerance){
          Atom_List[i].SetFreezeY(1);
	//this particular part of the if statement (testing if Y is locked by Z) is untested  
	}else if(fabs((fabs(COM[2]) - fabs(COM[1]))) < tolerance){
	  //Atom_List[i].LockZ(2);
	  Atom_List[i].LockY(1);
	  if(fabs(COM[2] - COM[1]) > tolerance){
	    //Atom_List[i].LockZ(1);
	    Atom_List[i].SetYGradSign(1);	   
	  }
	}

        //freezing z coordinate if it is zero
        if(fabs(COM[2]) < tolerance){
          //printf("Atom %i z COM coordinate is zero, freezing\n",Atom_List[i].GetGlobalIndex());
          Atom_List[i].SetFreezeZ(1);
        }
      }


      /*
      printf("Mon %i atom %i\n",Monomer_index,i);
      printf("Atom freezing = % x freeze = %i y freeze = %i z freeze = %i\n",
         Atom_List[i].IsAtomFrozen(),Atom_List[i].FreezeX(),Atom_List[i].FreezeY(),
         Atom_List[i].FreezeZ());
      printf("LockX = %i LockY = %i\n",Atom_List[i].IsXLocked(),Atom_List[i].IsYLocked());
      printf("ChangeX = %i ChangeY = %i\n",Atom_List[i].ChangeXSign(),Atom_List[i].ChangeYSign());
      */
    }//end of if statement checking if atom is in asymmetrical unit cell  
  }//end of loop over i
  if ( Params::Parameters().PrintLevel() > 0)
    printf("Number of unique_atoms = %i\n",unique_atoms);
  //exit(0);
}

//For turning off symmetry.
//Since atoms within a monomer may have already be found to be
//symmetrically equivalent, there symmetry info is reset.
void Monomer::ResetSymmetry(){
  sym_fac = 1;
  SymMon = Monomer_index;
  unique_atoms = Natoms;
  Symmetry_Rotation.Set_Iden();
  Fractional_Rotation.Set_Iden();
  Fractional_Translate.Set();
  for(int i=0;i<Natoms;i++){
    Sym_Atom[i] = i;
    GetAtom(i).ResetSymmetry();
  }
}

//Check to see if the coordinates of this momomer are unchanged after symmetry operation
bool Monomer::MatchCoordinates(Matrix Coord1,Matrix Coord2){
  
  //Coord1.Print("Coord1");
  //Coord2.Print("Coord2");


  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be

  bool CoordMatched = true;
  int i = 0;//index for atoms on this monomer
  while(i < Natoms && CoordMatched){
    int j=0;
    bool AtomMatched = false;
    while(j < Natoms && !AtomMatched){
     double x_diff = fabs(Coord1(0,i)-Coord2(0,j));
     double y_diff = fabs(Coord1(1,i)-Coord2(1,j));
     double z_diff = fabs(Coord1(2,i)-Coord2(2,j)); 
     bool SameAtomType = false;
     if(GetSymbol(i) == GetSymbol(j)) {
       SameAtomType = true;
   
    }
     if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance && SameAtomType){
       AtomMatched = true;
       //printf("atom %i = %s atom %i = %s\n",i,GetSymbol(i).c_str(),j,GetSymbol(j).c_str());
       //printf(" i %f %f %f \n",Coord1(0,i),Coord1(1,i),Coord1(2,i));
       //printf(" j %f %f %f \n\n",Coord2(0,j),Coord2(0,j),Coord2(2,j));   
     }
     
     if(j == Natoms-1 && !AtomMatched)
       CoordMatched = false;
     
    j++;      
   }//end of loop over j
   i++;

  }//end of loop over i

  return CoordMatched;

}

//Finding the translational operator in fractional coordinates to the unique monomer
void Monomer::SetFractTranslationVector(Monomer& UniqueMon){

  //Use the first atom and its symmetrical equivalent on UniqueMon
  Vector FractCoord = Atom_List[0].GetFractionalPosition();
  

  Vector SymFractCoord = UniqueMon.GetAtom(Sym_Atom[0]).GetFractionalPosition();
  //
  
  //Find translation
  Vector OperatedCoord = Fractional_Rotation.MatrixTimesVector(FractCoord);
  for(int i=0;i<3;i++)
    Fractional_Translate[i] = SymFractCoord[i] - OperatedCoord[i];

  //printf("Sym_Atom = %i\n",Sym_Atom);
  //FractCoord.Print("FractCoord");
  //SymFractCoord.Print("SymFractCoord");
  //OperatedCoord.Print("OperatedCoord");
  //Fractional_Translate.Print("Fractional_Translate");
  
  Matrix Inverse = Fractional_Rotation;
  Inverse.Inverse();
  Vector Back_Translate = Inverse.MatrixTimesVector(Fractional_Translate);
  Back_Translate.Scale(-1.0);
  ///Back_Translate.Print("Back_Translate");

  double tolerance = Params::Parameters().GetSymmetryTolerance();//tolerance
  
  //checking translation for all atoms
  for(int i=0;i<Natoms;i++){
    //does the translation work for all atoms
    bool translation_found = false;
    //Get the coordinates of the atoms
    Vector FractCoord = Atom_List[i].GetFractionalPosition();
    
    SymFractCoord = UniqueMon.GetAtom(Sym_Atom[i]).GetFractionalPosition();
    //apply symmetry rotation 
    OperatedCoord = Fractional_Rotation.MatrixTimesVector(FractCoord);
    
    //printf("i=%i Sym_Atom = %i\n",i,Sym_Atom);
    //FractCoord.Print("FractCoord");
    //SymFractCoord.Print("SymFractCoord");
    //OperatedCoord.Print("OperatedCoord");
    
    for(int j=0;j<3;j++){
      OperatedCoord[j] += Fractional_Translate[j];
      double diff =  OperatedCoord[j] - SymFractCoord[j];
      //printf("%i %i j=%i OperatedCoord=%f SymFractCoord=%f, diff = %f\n",
      //	     i,Sym_Atom,j,OperatedCoord[j],SymFractCoord[j],diff);
      if(fabs(diff) > tolerance){
	printf("ERROR::Monomer::SetFractTranslationVector:: Cannot find the translation operator between m%i and m%i\n",
	       Monomer_index,SymMon);
	exit(0);
      }
    }
    
    
    //Determining which coordinates are shifted by lattice parameters
    if(GetAtom(i).InAsymmetricUnit() 
       && Params::Parameters().UseSpaceSymmetry()){

      //If X or Y fractional coordinates are fixed to a 1/6,1/5,1/4,1/3,1/2 or these numbers multipled by an integer,
      //do none of the other check because they can get a false positive for their check and overid that check.
      //bool XFractToValue = 0;
      //bool YFractToValue = 0;
      bool XFixToValue = 0;
      bool YFixToValue = 0;

      //Determining if the x,y,z coordinates are a degree of freedom.
      //If they are, that coordiante can not be shifted by the lattice parameters
      bool XDOF = 1;
      bool YDOF = 1;
      bool ZDOF = 1;
      if(GetAtom(i).IsAtomFrozen()){
        XDOF = 0;
        YDOF = 0;
        ZDOF = 0;  
      }
      if(GetAtom(i).IsXLocked() || GetAtom(i).FreezeX()) {
        XDOF = 0;
      }
      if(GetAtom(i).IsYLocked() || GetAtom(i).FreezeY()){
        YDOF = 0;
      }
      if( GetAtom(i).FreezeZ() /*|| GetAtom(i).IsZLocked() ||*/){
        ZDOF = 0;
      }

      //if a fractional coordinate is zero, freeze it
      if((fabs(FractCoord[0]) < tolerance && !XDOF) && (fabs(FractCoord[1]) < tolerance) && !YDOF && (fabs(FractCoord[2]) < tolerance && !ZDOF) ){
	GetAtom(i).SetAtomToFrozen(1);
	XDOF = 0;
	YDOF = 0;
	ZDOF = 0;
      }
      else if(fabs(FractCoord[0]) < tolerance && XDOF){
	GetAtom(i).SetFreezeX(1);
	XDOF = 0;
      }
      else if(fabs(FractCoord[1]) < tolerance && YDOF){
	GetAtom(i).SetFreezeY(1);
	YDOF = 0;
      }
      else if(fabs(FractCoord[2]) < tolerance && ZDOF){
	GetAtom(i).SetFreezeZ(1);
	ZDOF = 0;
      }

      //printf("Mon%i atom%i XDOF = %i YDOF = %i ZDOF = %i\n",
      //        Monomer_index,i,XDOF,YDOF,ZDOF); 
     

      //A fractional coordinate is 1/6,1/5,1/4,1/3,1/2 or this numbers multipled by an integer, assume that it must be shifting by lattice parameters,
      //This coordinate can no longer be degree of freedom for optimzation
      for(int j = 1; j<=6;j++){
	
	//shift vector for frozen x
	if(fabs(FractCoord[0] - 0.25*double(j))  < tolerance && !XDOF){
	  GetAtom(i).SetShiftByA(0.25*double(j));
	  XFixToValue = 1;
	}else if(fabs(FractCoord[0] + 0.25*double(j)) < tolerance && !XDOF){
	  GetAtom(i).SetShiftByA(-0.25*double(j));
	  XFixToValue = 1;
	}else if(fabs(FractCoord[0] - 1/6*double(j)) < tolerance && !XDOF){
	  GetAtom(i).SetShiftByA(1/6*double(j));
	  XFixToValue = 1;
	}else if(fabs(FractCoord[0] + 1/6*double(j)) < tolerance && !XDOF){
	  GetAtom(i).SetShiftByA(-1/6*double(j));
	  XFixToValue = 1;	  
	}
	  
	//shift vector for frozen y
	if(fabs(FractCoord[1] - 0.25*double(j))  < tolerance && !YDOF){
	  GetAtom(i).SetShiftByB(0.25*double(j));
	  YFixToValue = 1;
	  //YFractToValue = 1;
	}else if(fabs(FractCoord[1] + 0.25*double(j)) < tolerance && !YDOF){
	  GetAtom(i).SetShiftByB(-0.25*double(j));
	  YFixToValue = 1;
	  //YFractToValue = 1;
	}else if(fabs(FractCoord[1] - 1/6*double(j)) < tolerance && !YDOF){
	  GetAtom(i).SetShiftByB(1/6*double(j));
	  YFixToValue = 1;
	  //YFractToValue = 1;
	}else if(fabs(FractCoord[1] + 1/6*double(j)) < tolerance && !YDOF){
	  GetAtom(i).SetShiftByB(-1/6*double(j));
	  YFixToValue = 1;
	  //YFractToValue = 1;
	}
       

	//shift vector for frozen z
	if(fabs(FractCoord[2] - 0.25*double(j))  < tolerance && !ZDOF){
	  GetAtom(i).SetShiftByC(0.25*double(j));
	  //ZFractToValue = 1;
	}else if(fabs(FractCoord[2] + 0.25*double(j)) < tolerance && !ZDOF){
	  GetAtom(i).SetShiftByC(-0.25*double(j));
	  //ZFractToValue = 1;
	}else if(fabs(FractCoord[2] - 1/6*double(j)) < tolerance && !ZDOF){
	  GetAtom(i).SetShiftByC(1/6*double(j));
	  //ZFractToValue = 1;
  	}else if(fabs(FractCoord[2] + 1/6*double(j)) < tolerance && !ZDOF){
	  GetAtom(i).SetShiftByC(-1/6*double(j));
	  //ZFractToValue = 1;
	}
      }
      //double xy_diff = FractCoord[1] - FractCoord[0];
      //double xz_diff = FractCoord[2] - FractCoord[0];
      //double yz_diff = FractCoord[2] - FractCoord[1];
      
      //double xy_sum = FractCoord[1] + FractCoord[0];
      //double xz_sum = FractCoord[2] + FractCoord[0];
      //double yz_sum = FractCoord[2] + FractCoord[1];

      double xy_diff = FractCoord[0] - FractCoord[1];  
      double xz_diff = FractCoord[0] - FractCoord[2];
      double yz_diff = FractCoord[1] - FractCoord[2];
      
      double xy_sum = FractCoord[0] + FractCoord[1];
      double xz_sum = FractCoord[0] + FractCoord[2];
      double yz_sum = FractCoord[1] + FractCoord[2];

      
      //Z shifting by C to match the x coordinates should overid shifting by Z to match y
      //bool ZLinkByX = 0;

      //X shifting by A to match the Z coordinates should overid shifting by x to match y
      bool XLinkByZ = 0;
      
      
      //If two fractional coordinates differ by only 0.5*j, assume that it must be shifting by lattice parameters,
      //This coordinate can no longer be degree of freedom for optimzation.
      for(int j=-3; j<=3;j++){

	//if criteria is met when j = 0, symmetry preservation already taken care off by lock_y and lock_z varibles.
	if(j != 0){
	  //printf("0.5*j = %f\n", double(j)*0.5);
	  

	  if(fabs(xz_diff - double(j)*0.5) < tolerance  && !XFixToValue && !XDOF){
	    //printf("ShiftByC xz_sum\n");
	    GetAtom(i).SetShiftByA(double(j)*0.5);
	    XLinkByZ = 1;
	  }
	  else if(fabs(xy_diff - double(j)*0.5) < tolerance && !XFixToValue && !XDOF && !XLinkByZ){
	    GetAtom(i).SetShiftByA(double(j)*0.5);
	  }

	  if(fabs(yz_diff - double(j)*0.5) < tolerance && !YFixToValue  && !YDOF ){
	    GetAtom(i).SetShiftByB(double(j)*0.5);
	  }

	  if(fabs(xz_sum - double(j)*0.5) < tolerance  && !XFixToValue && !XDOF){
	    //GetAtom(i).SetShiftByC(0.5*double(j));
	    //printf("ShiftByC xz_sum\n");
	    GetAtom(i).SetShiftByA(0.5*double(j));
	    XLinkByZ = 1;
	  }
	  else if(fabs(xy_sum - double(j)*0.5) < tolerance && !XFixToValue && !XDOF && !XLinkByZ){
	    GetAtom(i).SetShiftByA(double(j)*0.5);
	  }

	  if(fabs(yz_sum - double(j)*0.5) < tolerance && !YFixToValue  && !YDOF ){
	    GetAtom(i).SetShiftByB(double(j)*0.5);
	  }
	  

	  /*
	  if(fabs(xy_diff - double(j)*0.5) < tolerance && !YFractToValue && !YDOF){
	    GetAtom(i).SetShiftByB(0.5*double(j));
	    //printf("ShiftByB xy_diff\n");
	  }
	  if(fabs(xz_diff - double(j)*0.5) < tolerance && !ZFractToValue && !ZDOF){
	    GetAtom(i).SetShiftByC(0.5*double(j));
	    ZLinkByX = 1;
	    //printf("ShiftByC xz_diff\n");
	    
	  }else if(fabs(yz_diff - double(j)*0.5) < tolerance && !ZLinkByX && !ZFractToValue && !ZDOF){
	    GetAtom(i).SetShiftByC(0.5*double(j));
	    //printf("ShiftByC yz_diff\n");
	  }
	  if(fabs(xy_sum - double(j)*0.5) < tolerance  && !YFractToValue && !YDOF) {
	    GetAtom(i).SetShiftByB(0.5*double(j));
	    // printf("ShiftByB xy_sum\n");
	  }
	  if(fabs(xz_sum - double(j)*0.5) < tolerance && !ZFractToValue && !ZDOF){
	    GetAtom(i).SetShiftByC(0.5*double(j));
	    //printf("ShiftByC xz_sum\n");
	    ZLinkByX = 1;
	  }else if(fabs(yz_sum - double(j)*0.5) < tolerance && !ZLinkByX  && !ZFractToValue && !ZDOF){
	    GetAtom(i).SetShiftByC(0.5*double(j));
	    //printf("ShiftByC yz_sum\n");
	  }
	  */
	  
	} 
      }
      // printf("Mon%i Atom%i shiftbya = %f shiftbyb = %f shiftbyc = %f\n",
      //       Monomer_index,i,GetAtom(i).ShiftByA(),
      //     GetAtom(i).ShiftByB(),GetAtom(i).ShiftByC());
    }   
    
    //Determining the translation vector from the symmetrical atom in the asymetric unit (non-monomer)
    else{

      //Global index for symmetrical atom in the unit cell
      int SymAtomGlob = GetAtom(i).GetSymmetricalAtom();
      
      //index of SymAtomGlob in the monomer
      int SymAtom = -1;
      


      for(int j=0;j<Natoms;j++){
      
	if(UniqueMon.GetAtom(j).GetGlobalIndex() == SymAtomGlob){
	  SymAtom = j;
	}
      }
	

      if(SymAtom == -1){
	printf("Monomer::SetFractTranslationVector : Cannot determine which atom atom%i is symmetrical to in asymetric cell\n",i);
	exit(0);
      }
      //printf("i = %i SymAtom = %i SymAtomGlob = %i\n",i,SymAtom,SymAtomGlob);
      

      
      Vector FractCoord = GetAtom(i).GetFractionalPosition();
      Vector SymFractCoord = UniqueMon.GetAtom(SymAtom).GetFractionalPosition();
      
      
      //Find translation
      Vector Atom_Translate(3);
      Matrix Rot = GetAtom(i).GetFractionalRotationMatrix();
      //Rot.Inverse();
      Vector OperatedCoord = Rot.MatrixTimesVector(FractCoord);


      //Rot.Print("Rot");
      //printf("\nmon%i atom%i\n",Monomer_index,i);  
      //FractCoord.Print("FractCoord");
      //SymFractCoord.Print("SymFractCoord");
      //OperatedCoord.Print("OperatedCoord");

      for(int j=0;j<3;j++)
	Atom_Translate[j] =  SymFractCoord[j] - OperatedCoord[j];
      GetAtom(i).SetFractTranslationVector(Atom_Translate);


      //Since some atoms in the asymetrical unit shift fractional coordinates in order
      //to preserve symmetry while preserving degrees of freedoms Cartesian coordinates,
      //corresponding shifts are  found here
      Vector ShiftVector = UniqueMon.GetAtom(SymAtom).GetShiftVector();
      Rot.Inverse();
      ShiftVector = Rot.MatrixTimesVector(ShiftVector);
      GetAtom(i).SetShiftVector(ShiftVector);


    }
    if (Params::Parameters().PrintSymmetryInfo()) {
      printf("\nmon%i atom%i\n",Monomer_index,i);  
      //GetRotationMatrix().Print("Cart Mon Rot");
      GetAtom(i).GetFractionalRotationMatrix().Print("Rot");
      //GetAtom(i).GetRotationMatrix().Print("Cart Rot");
      GetAtom(i).GetFractTranslationVector(false).Print("atom Translation");
      //GetAtom(i).GetFractTranslationVector(true).Print("atom Back Translation");
      GetAtom(i).GetShiftVector().Print("Shift Vector");;
      //GetFractTranslationVector(true).Print("mon Back Translation");
    }
  }


}

//Function returns the translation part of the symmetry operator
Vector Monomer::GetFractTranslationVector(bool FromUniqueToThis){


  Vector TransOp;

  //The translational operator is part of the symmetry operator going from the symmetry unique monomer to this one.
  if(FromUniqueToThis){
    Matrix InverseRot = Fractional_Rotation;
    InverseRot.Inverse();
    //InverseRot.Print("InverseRot");
    TransOp = InverseRot.MatrixTimesVector(Fractional_Translate);
    TransOp.Scale(-1.0);
    //TransOp.Print("TransOp");
  }else 
    TransOp = Fractional_Translate;


  //TransOp.Print("TransOp");

  /*
  //test out ReverseTranslation
  double tolerance = Params::Parameters().GetSymmetryTolerance();//tolerance
  printf("m%i\n",Monomer_index);
  for(int i=0;i<Natoms;i++){
    Vector FractCoord = Atom_List[i].GetFractionalPosition();

    
    int Sym_Atom = -1;
    for(int j=0; j<Natoms;j++)
      if(Atom_List[i].GetSymmetricalAtom() == UniqueMon.GetAtom(j).GetGlobalIndex())
	Sym_Atom = j;
    if(Sym_Atom == -1){
      printf("ERROR::Monmomer::SetFractTranslationVector:: Unable to determine Symmetrical atom to atom %i on m%i \n",
	     Atom_List[i].GetSymmetricalAtom(),Monomer_index);
    }

    Vector SymFractCoord = UniqueMon.GetAtom(Sym_Atom).GetFractionalPosition();
    
    Vector OperatedCoord = InverseRot.MatrixTimesVector(SymFractCoord);
    for(int j=0;j<3;j++)
      OperatedCoord[j] += ReverseTranslation[j];
    printf("i=%i Sym_Atom = %i\n",i,Sym_Atom);
    FractCoord.Print("FractCoord");
    SymFractCoord.Print("SymFractCoord");
    OperatedCoord.Print("OperatedCoord");
    

    
    for(int j=0;j<3;j++){
      double diff =  OperatedCoord[j] - FractCoord[j];
      //printf("%i %i j=%i OperatedCoord=%f FractCoord=%f, diff = %f\n",
      //	     i,Sym_Atom,j,OperatedCoord[j],FractCoord[j],diff);
      if(fabs(diff) > tolerance){
	printf("ERROR::Monomer::GetFractReverseTranslationVector:: Cannot find the translation operator between m%i and m%i\n",
	       Monomer_index,SymMon);
	exit(0);
      }
    }
  }
    */
  
  return TransOp;


}

//OUTDATED, NO LONGER USED.
//determines what the transformation matrix is to the unique symmetrical monomer
//Ax = y
//where A is the transforming matrix from the monomer to the unique one.
void Monomer::FindMapping(Monomer& UniqueMon){
  
  printf("Finding mapping to monomer %i to monomer %i\n",
	 GetIndex(), UniqueMon.GetIndex());
  bool RotationFound=0;
  int MonomerType=1;//used to determine the method to find the transformation matrix
                    //0 for non-planar monomers
                    //1 for planar monomers
                    //2 for linear monomers
                    //3 for single atom monomers
                    //set to planar by default



  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be
  //Center of mass coordinates 
  Matrix COM(3,Natoms);
  COM.Set();
  for(int i=0;i<Natoms;i++) {
    COM(0,i) = Atom_List[i].GetPosition(0)-GetCenterOfMass(0);// x center of mass
    COM(1,i) = Atom_List[i].GetPosition(1)-GetCenterOfMass(1);// y center of mass
    COM(2,i) = Atom_List[i].GetPosition(2)-GetCenterOfMass(2);// y center of mass
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //		       GetAtom(i).GetSymbol().c_str(), COM(0,i),COM(1,i),COM(2,i));
  }


  
  if(Natoms==1)
    MonomerType=3;//single atom monomer
  else if(Natoms==2)
    MonomerType=2;//two atom monomers are always linear
  else if(Natoms==3)
    MonomerType=1;//three atom monomers may be linear or planar but will have to check
  else{

    //matrix of atoms, multiple inverse by the matrix made by the unique monomer to find the transformation matrix 
    Matrix AtomMatrix(3,3);

    //find the atom symbol to be used to match ones on unique matrix
    string Symbol[3];

    //this is to determine which atoms form the matrix 
    int first = 0;
    int second;
    int third;

    //the determinant is used to see if the molecule is planar
    double det = 0;
    
    //Looping over all the atoms trying if a matrix made
    // from their coordinates will make a linear independent system
    while((first < Natoms) && (fabs(det) < tolerance)){
      second = first+1;
      while((second < Natoms) && (fabs(det) < tolerance)){
	third = second+1;
	while((third < Natoms) && (fabs(det) < tolerance)){
	  //Making matrix from the coordinates from 3 atoms
	  AtomMatrix.SetColumnVector(COM.GetColumnVector(first),0);
	  AtomMatrix.SetColumnVector(COM.GetColumnVector(second),1);
	  AtomMatrix.SetColumnVector(COM.GetColumnVector(third),2);
	  Symbol[0] = GetAtom(first).GetSymbol();
	  Symbol[1] = GetAtom(second).GetSymbol();
	  Symbol[2] = GetAtom(third).GetSymbol();
	  //printf("Atoms of the matrix is \n");
	  // for(int i=0;i<3;i++){
	  //  printf("%10.6f  %10.6f  %10.6f\n",
	  //	   AtomMatrix(i,0),AtomMatrix(i,1),AtomMatrix(i,2));
	  //}
	  //the determinant function overwrites the matrix
	  //creating a temp matrix to prevent that
	  Matrix tempMat = AtomMatrix;
	  det = tempMat.Determinant();
	  //printf("Determinant is %f \n", det);
	  if(fabs(det) > tolerance){
	    MonomerType=0;//non-planar monomers
	  }
	  third++;
	}
	second++;
      }
      first++;
    }
    //Find the transformation matrix for non-planar molecules
    if(MonomerType==0)
      RotationFound =  MappingMatrixForNonPlanar(UniqueMon,COM,AtomMatrix,Symbol);
  }
  
  //Next if monomer is not non-planar, determine if it is planar or linear
  if(MonomerType==1){
    bool IsLinear=1;//varible to determine if molecule is linear, set to true by default
    //if molecule is not linear, the cross product between any of its atoms will be nonzero
    //only need two atoms for planar molecules
    int first=0;
    int second;
    Matrix AtomMatrix(3,2);//this matrix will made from the atoms of the first matrix
    //will be rotated to match the coordinates of the second 
    
    //find the atom symbol to be used to match ones on unique matrix
    string Symbol[3];
    
    while(first<Natoms && IsLinear==1){
      second=first+1;
      while(second<Natoms && IsLinear==1){
	//if the cross product between the any two atoms is nonzero, the molecule is not linear
	Vector FirstAtom = COM.GetColumnVector(first);
	Vector SecondAtom = COM.GetColumnVector(second);
	double CrossMag = FirstAtom.CrossProduct(SecondAtom).Norm();
	if(CrossMag > 0.001){
	  AtomMatrix.SetColumnVector(FirstAtom,0);
	  AtomMatrix.SetColumnVector(SecondAtom,1);
	  Symbol[0] = GetAtom(first).GetSymbol();
	  Symbol[1] = GetAtom(second).GetSymbol();
	  IsLinear=0;
	}
	second++;
      } 
      first++;
    }
    if(IsLinear==0){
      //this part is to find the transformation matrix for planar molecules
      RotationFound = MappingMatrixForPlanar(UniqueMon,COM,AtomMatrix,Symbol);
    }else
      MonomerType=2;
  }
  
  if(MonomerType==2){
    RotationFound = MappingMatrixForLinear(UniqueMon,COM);
  }
  
  if(RotationFound==1){
    printf("Rotation found\n");
    /*
    //Cleaning up any numerical noice in rotation matirix
    Matrix Rotation = GetRotationMatrix();
    Matrix AdjustedRotation(3,3);
    AdjustedRotation.Set();
    printf("tolerance = %f\n", tolerance);
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	if(fabs(Rotation(i,j)-1) < tolerance)//R(i,j)=1
	  AdjustedRotation(i,j)=1;
	else if(fabs(Rotation(i,j)) < tolerance)//R(i,j)=0
	  AdjustedRotation(i,j)=0;
	else if(fabs(Rotation(i,j)+1) < tolerance)//R(i,j)=-1
	  AdjustedRotation(i,j)=-1;
	else if(fabs(Rotation(i,j)-0.5) < tolerance)//R(i,j)=0.5 if rotated by 30 or 60 degrees
	  AdjustedRotation(i,j) = 0.5;
	else if(fabs(Rotation(i,j)+0.5) < tolerance)//R(i,j)=-0.5 if rotated by 30 or 60 degrees
	  AdjustedRotation(i,j) = -0.5;
	else if(fabs(Rotation(i,j)-0.866025) < tolerance)//R(i,j) = 0.866025 if rotated by 30 or 60 degrees
	  AdjustedRotation(i,j) = 0.866025;
	else if(fabs(Rotation(i,j)+0.866025) < tolerance)//R(i,j) = -0.886025 if rotated by 30 or 60 degrees
	  AdjustedRotation(i,j) = -0.866025;
	else
	  AdjustedRotation(i,j) = Rotation(i,j);
    AdjustedRotation.Print("Rotation");
    SetRotationMatrix(AdjustedRotation);
    */
  }
  else{
    printf("No Rotation found");
    exit(0);
  }
}


//OUTDATED, NO LONGER USED.
bool Monomer::MappingMatrixForNonPlanar(Monomer& UniqueMon,Matrix COM, Matrix AtomMatrix, string Symbol[]){

  Matrix RotationMatrix(3,3);

  //Center of mass unique monomer
  Matrix UniqueCOM(3,Natoms);
  UniqueCOM.Set();
  //printf("Unique atom\n");
  for(int i=0;i<Natoms;i++){
    UniqueCOM(0,i) = UniqueMon.GetAtom(i).GetPosition(0)- UniqueMon.GetCenterOfMass(0);// x center of mass
    UniqueCOM(1,i) = UniqueMon.GetAtom(i).GetPosition(1)-UniqueMon.GetCenterOfMass(1);// y center of mass
    UniqueCOM(2,i) = UniqueMon.GetAtom(i).GetPosition(2)-UniqueMon.GetCenterOfMass(2);// z center of mass
    // printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
	   //	   UniqueMon.GetAtom(i).GetSymbol().c_str(), UniqueCOM(0,i),UniqueCOM(1,i),UniqueCOM(2,i));
  }
  
  //magnitude
  double mag1 = AtomMatrix.GetColumnVector(0).Norm();
  double mag2 = AtomMatrix.GetColumnVector(1).Norm();
  double mag3  = AtomMatrix.GetColumnVector(2).Norm();
  //printf("Distance for COM of the Atoms used to make matrix are %f %f %f\n",mag1,mag2,mag3);
  //fflush(stdout);
  Matrix UniqueAtomMatrix(3,3);

  //Three loops over the atoms of the unique monomer.
  //Trying to make a 3x3 matrix using equivalent atoms to
  //matrix maded by the first monomer. Using these to matrices
  //the rotation matrix will be found.
  double  det = 0;//determinant for the unique monomer
  int first = 0;//first atom in the 3X3 matrix
  bool DoesNotMatch = 0;//this will tell us if the transformation is correct
  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be

  do{//looping over the "first" atom
    int second = 0;//first atom in the 3X3 matrix
    do{//looping over the "second" atom
      int third = 0;//third atom in the 3X3 matrix
      do{//looping over the "third" atom
	//printf("atoms looking at are %i,%i,%i\n",
	//       first,second,third);
	
	//looking for conditions where it is known that the unique matrix will not be mapped to.
	bool skiploop = 0;
	double UniqueMag1 = UniqueCOM.GetColumnVector(first).Norm(); 
	double UniqueMag2 = UniqueCOM.GetColumnVector(second).Norm(); 
	double UniqueMag3 = UniqueCOM.GetColumnVector(third).Norm();
	//if any of the atoms in the matrices are the same, the mapping matrix will be incorrect
	if(first==second || second==third || first==third)
	  skiploop = 1;
	
	//if the elements of the atoms in the matrices are not the same, the mapping will be incorrect 
	if(UniqueMon.GetAtom(first).GetSymbol()!=Symbol[0] || UniqueMon.GetAtom(second).GetSymbol()!=Symbol[1]||UniqueMon.GetAtom(third).GetSymbol()!=Symbol[2])
	  skiploop = 1;
	
	//if the length of the atoms between the two molecules are not the same, the mapping matrix will be incorrect.
	if(fabs(UniqueMag1-mag1)>tolerance || fabs(UniqueMag2-mag2)>tolerance ||fabs(UniqueMag3-mag3)>tolerance)
	  skiploop = 1;
	
	if(skiploop)
	  DoesNotMatch=1;
	else{
	  //printf("magnitude of the unique matrix is %f %f %f\n",mag1,mag2,mag3);
	  UniqueAtomMatrix.SetColumnVector(UniqueCOM.GetColumnVector(first),0);
	  UniqueAtomMatrix.SetColumnVector(UniqueCOM.GetColumnVector(second),1);
	  UniqueAtomMatrix.SetColumnVector(UniqueCOM.GetColumnVector(third),2);
	  // printf("Atoms of the unique matrix is \n");
	  // for(int i=0;i<3;i++){
	  //   printf("%10.6f  %10.6f  %10.6f\n",
	  //	   UniqueAtomMatrix(i,0),UniqueAtomMatrix(i,1),UniqueAtomMatrix(i,2));
	  //}
	  Matrix tempUnique = UniqueAtomMatrix;
	  det = tempUnique.Determinant();
	  //delete [] tempUnique; how do you delete
	  //printf("Determinant of unique matrix %f\n",det);
	  fflush(stdout);
	  Matrix tempinverse = AtomMatrix;
	  tempinverse.Inverse();
	  RotationMatrix = UniqueAtomMatrix.Multiply(tempinverse);
	  //printf("Rotation Matrix is\n");
	  //for(int i=0;i<3;i++){
	  //  printf("%10.6f  %10.6f  %10.6f\n",
	  //   RotationMatrix(i,0),RotationMatrix(i,1),RotationMatrix(i,2));
	  //}
	  Matrix AfterT = RotationMatrix.Multiply(COM);
	  // printf("Using Rotation Matrix is\n");
	  //for(int i=0;i<Natoms;i++){
	  //  printf("%-2s %10.6f  %10.6f  %10.6f\n",
	  //	   GetAtom(i).GetSymbol().c_str(), AfterT(0,i),AfterT(1,i),AfterT(2,i));
	  //}
	  int i=0;//atoms on first molecule
	  DoesNotMatch=0;//used to identify whether the rotation was correct
	  while(i<Natoms && DoesNotMatch==0){
	    int j = 0;//atoms on unique molecule
	    bool match=0;//match for each atoms
	    while(j<Natoms && match==0){
	      if(GetAtom(i).GetSymbol() == UniqueMon.GetAtom(j).GetSymbol()){
		double x_diff = fabs(AfterT(0,i)-UniqueCOM(0,j));
		double y_diff = fabs(AfterT(1,i)-UniqueCOM(1,j));
		double z_diff = fabs(AfterT(2,i)-UniqueCOM(2,j));
		//printf("i=%i and j=%i\n",i,j);
		// printf("differance of atom %i and atom %i is %f %f %f\n",
		//	   GetAtom(i).GetAtomIndex(),UniqueMon.GetAtom(j).GetAtomIndex(), x_diff,y_diff,z_diff);
		//printf("Atom %f %f %f,Unique %f %f %f\n",
		//	   AfterT(0,i),AfterT(1,i),AfterT(2,i),UniqueCOM(0,j),UniqueCOM(1,j),UniqueCOM(2,j));
		if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance){
		  match=1;//identical atom found
		  GetAtom(i).SetSymmetricalAtom(UniqueMon.GetAtom(j).GetGlobalIndex());
		  //printf("atom %-2s is idential to unique atom %-2s\n",
			 //GetAtom(i).GetSymbol().c_str(), UniqueMon.GetAtom(j).GetSymbol().c_str());
		}
	      }
	      j++;
	    }
	    if(match==0){
	      DoesNotMatch=1;
	    }
	    i++;	       
	  }
	}
	third++;    
      }while(third<Natoms && DoesNotMatch==1);
      second++;
    }while(second<Natoms && DoesNotMatch==1);
    first++;
  }while(first<Natoms && DoesNotMatch==1);

  if(DoesNotMatch==0){
    SetRotationMatrix(RotationMatrix);
    RotationMatrix.Print("Rotation Matrix is");
    return 1;
  }else
    return 0;
}

//OUTDATED, NO LONGER USED.
bool Monomer::MappingMatrixForPlanar(Monomer& UniqueMon, Matrix COM, Matrix AtomMatrix, string Symbol[]){


  Matrix RotationMatrix(3,3);
  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be

  //AtomMatrix.Print("");
  //This function will find a series of rotations to find out to tranform the
  //"first" monomer to the "unique" one if the monomers are planar.
  // First both monomers will be rotated into xy plane. Then find a 2D 
  // transformation matrix (equivalent to a rotation around the Z-axis).
  // The operation will have the following form
  //
  // T2D*Rx2*Rz1*Rx1*m1= Rx4*Rz2*Rx3*m2
  //
  //A total transformation matrix is the product of all the operations
  //
  //T(total) = Rx3'*Rz2'*Rx4'*T2D*Rx2*Rz1*Rx1
  //
  //T(total)*m1=m2


  //AtomMatrix.Print("AtomMatrix");

  //Rotating the first monomer to the XY plane
  //Step One: rotating the molecule along the X-axiz so the first atom is in the xy plane
  Matrix Rx1(3,3);
  double angle;
  if(fabs(AtomMatrix(1,0))<tolerance)
    angle= pi/2;
  else
    angle = atan(AtomMatrix(2,0)/AtomMatrix(1,0));

  Rx1.Set_Iden();
  Rx1(1,1)=cos(angle);
  Rx1(1,2)=sin(angle);
  Rx1(2,1)=-sin(angle);
  Rx1(2,2)=cos(angle);
  AtomMatrix = Rx1.Multiply(AtomMatrix);

  //AtomMatrix.Print("AtomMatrix after Rx1");

  //Step Two: rotating along the Z-axis so the first atom is parallel to x-axis
  Matrix Rz1(3,3);
  if(fabs(AtomMatrix(0,0))<tolerance)
    angle= pi/2;
  else
    angle = atan(AtomMatrix(1,0)/AtomMatrix(0,0));
  
  Rz1.Set_Iden();
  Rz1(0,0)=cos(angle);
  Rz1(0,1)=sin(angle);
  Rz1(1,0)=-sin(angle);
  Rz1(1,1)=cos(angle);
  AtomMatrix = Rz1.Multiply(AtomMatrix);

  //AtomMatrix.Print("AtomMatrix after Rz1");

  // Step Three: rotate along x-axis so entire molecule is in xy plane
  Matrix Rx2(3,3);
  if(fabs(AtomMatrix(1,1)) < tolerance){
    angle = pi/2;
  }else{
    angle = atan(AtomMatrix(2,1)/AtomMatrix(1,1));
  }
  Rx2.Set_Iden();
  Rx2(1,1)=cos(angle);
  Rx2(1,2)=sin(angle);
  Rx2(2,1)=-sin(angle);
  Rx2(2,2)=cos(angle);
  AtomMatrix = Rx2.Multiply(AtomMatrix);

  // AtomMatrix.Print("AtomMatrix after Rx2");


  //Multiply all the transpose to determine the matrix that rotates the atoms from the plane
  //into it current coordinates.
  Matrix FirstRotation = Rx2.Multiply(Rz1);
  FirstRotation = FirstRotation.Multiply(Rx1);
  //FirstRotation.Transpose();
  //FirstRotation.Print("First rotation");
  
  //this is to be used later for finding the rotation when both monomers
  //are rotated into the xy-axis
  Matrix squareAtomMatrix(3,3);
  squareAtomMatrix.Set_Iden();
  squareAtomMatrix.SetColumnVector(AtomMatrix.GetColumnVector(0),0);
  squareAtomMatrix.SetColumnVector(AtomMatrix.GetColumnVector(1),1); 

  //Center of mass unique monomer
  Matrix UniqueCOM(3,Natoms);
  UniqueCOM.Set();
  //printf("Unique atom\n");
  for(int i=0;i<Natoms;i++){
    UniqueCOM(0,i) = UniqueMon.GetAtom(i).GetPosition(0)- UniqueMon.GetCenterOfMass(0);// x center of mass
    UniqueCOM(1,i) = UniqueMon.GetAtom(i).GetPosition(1)-UniqueMon.GetCenterOfMass(1);// y center of mass
    UniqueCOM(2,i) = UniqueMon.GetAtom(i).GetPosition(2)-UniqueMon.GetCenterOfMass(2);// z center of mass
    // printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   UniqueMon.GetAtom(i).GetSymbol().c_str(), UniqueCOM(0,i),UniqueCOM(1,i),UniqueCOM(2,i));
  }

  //magnitude
  double mag1 = AtomMatrix.GetColumnVector(0).Norm();
  double mag2 = AtomMatrix.GetColumnVector(1).Norm();

  int first = 0;//first atom in the 3X2 matrix
  bool DoesNotMatch = 0;//this will tell us if the transformation done is correct
  do{//first atom loop
    int second = 0;
    do{//second atom loop

      // printf("atoms looking at are %i,%i\n",
      //      first,second);
	
	Vector Atom1 =  UniqueCOM.GetColumnVector(first);
	Vector Atom2 =  UniqueCOM.GetColumnVector(second);
	
	//looking for conditions were it is known that the unique matrix will not be mapped to.
	bool skiploop = 0;
	double UniqueMag1 = Atom1.Norm(); 
	double UniqueMag2 = Atom2.Norm();
	double CrossMag = Atom1.CrossProduct(Atom2).Norm();
	
	//if any of the atoms in the matrices are the same, the mapping matrix will be incorrect
	if(first==second)
	  skiploop = 1;
	
	//if the elements of the atoms in the matrices are not the same, the mapping will be incorrect 
	if(UniqueMon.GetAtom(first).GetSymbol()!=Symbol[0] || UniqueMon.GetAtom(second).GetSymbol()!=Symbol[1])
	  skiploop = 1;

	//if the cross Product of the two vectors are the same, they are parallel
	if(CrossMag < 0.001)
	  skiploop = 1;
	
	//if the length of the atoms between the two molecules are not the same, the mapping matrix will be incorrect.
	if(fabs(UniqueMag1-mag1)>tolerance||fabs(UniqueMag2-mag2)>tolerance)
	  skiploop = 1;

	if(skiploop)
	  DoesNotMatch=1;
	else{

	  Matrix UniqueAtomMatrix(3,2);
	  UniqueAtomMatrix.SetColumnVector(Atom1,0);
	  UniqueAtomMatrix.SetColumnVector(Atom2,1);
	  //UniqueAtomMatrix.Print("");

	  //performing the same opertations the atoms of the unique monomer as the other ones
	  //we are determining how to rotate the unique monomer so it is in the xy plane

	  //UniqueAtomMatrix.Print("UniqueAtomMatrix");

	  //Step One: rotating the along the X-axiz so the first atom is in the xy plan
	  Matrix Rx1(3,3);
	  double angle;
	  if(fabs(UniqueAtomMatrix(1,0))<tolerance)
	    angle= pi/2;
	  else
	    angle = atan(UniqueAtomMatrix(2,0)/UniqueAtomMatrix(1,0));

	  Rx1.Set_Iden();
	  Rx1(1,1)=cos(angle);
	  Rx1(1,2)=sin(angle);
	  Rx1(2,1)=-sin(angle);
	  Rx1(2,2)=cos(angle);
	  UniqueAtomMatrix = Rx1.Multiply(UniqueAtomMatrix);

	  // UniqueAtomMatrix.Print("UniqueAtomMatrix after Rx1");

	  //Step Two: rotating along the Z-axis so the first atom is parallel to x-axis
	  Matrix Rz1(3,3);
	  if(fabs(UniqueAtomMatrix(0,0))<tolerance)
	    angle= pi/2;
	  else
	    angle = atan(UniqueAtomMatrix(1,0)/UniqueAtomMatrix(0,0));

	  Rz1.Set_Iden();
	  Rz1(0,0)=cos(angle);
	  Rz1(0,1)=sin(angle);
	  Rz1(1,0)=-sin(angle);
	  Rz1(1,1)=cos(angle);
	  UniqueAtomMatrix = Rz1.Multiply(UniqueAtomMatrix);

	  //UniqueAtomMatrix.Print("UniqueAtomMatrix after Rz1");
	  
	  // Step Three: rotate along x-axis so entire molecule is in the xy plane
	  Matrix Rx2(3,3);
	  if(fabs(UniqueAtomMatrix(1,1)) < tolerance){
	    angle = pi/2;

	  }else{
	    angle = atan(UniqueAtomMatrix(2,1)/UniqueAtomMatrix(1,1));
	  }
	  Rx2.Set_Iden();
	  Rx2(1,1)=cos(angle);
	  Rx2(1,2)=sin(angle);
	  Rx2(2,1)=-sin(angle);
	  Rx2(2,2)=cos(angle);
	  UniqueAtomMatrix = Rx2.Multiply(UniqueAtomMatrix);

	  //UniqueAtomMatrix.Print("UniqueAtomMatrix after Rx2");

	  //Multiply all the transpose to determine the matrix that rotates the atoms from the plane
	  //into it current coordinates.
	  Matrix SecondRotation = Rx2.Multiply(Rz1);
	  SecondRotation = SecondRotation.Multiply(Rx1);

	  Matrix temp_COM = SecondRotation.Multiply(UniqueCOM);// used for testing only
	  printf("UniqueCom rotated into XY plane of m%i\n",UniqueMon.GetIndex());
	  for(int i=0;i<Natoms;i++){
	    printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
		   UniqueMon.GetAtom(i).GetSymbol().c_str(), temp_COM(0,i),temp_COM(1,i),temp_COM(2,i));
	  } 

	  

	  //Now that both molecules are in the xy-plane
	  //A 2D transformation matrix will be found.
	  //This matrix will have the form of a rotation about
	  //the z-axis when the 3rd dimension is used.


	  //T*X = Y
	  //T= Y*X^-1
	  Matrix inverseUnique(3,3);
	  inverseUnique.Set_Iden();
	  inverseUnique.SetColumnVector(UniqueAtomMatrix.GetColumnVector(0),0);
	  inverseUnique.SetColumnVector(UniqueAtomMatrix.GetColumnVector(1),1);
	  inverseUnique.Inverse();		
	  Matrix Transform = squareAtomMatrix.Multiply(inverseUnique);

	  //putting all the operators together
	  //(T*R1)*m2 = R2*m1
	  //(R2^-1*T*R1)*m2 = m1
	  SecondRotation.Transpose();
	  RotationMatrix = SecondRotation.Multiply(Transform);
	  RotationMatrix = RotationMatrix.Multiply(FirstRotation);

	  //squareAtomMatrix.Print("squareAtomMatrix");
	  //UniqueAtomMatrix.Print("UnqiueAtomMatrix");
	  //inverseUnique.Print("InverseUnique");
	  
	  // FirstRotation.Print("First Rotation");
	  //SecondRotation.Print("Second Rotation");
	  //Transform.Print("Rotation in xy plane");
	  RotationMatrix.Print("Total Rotation");

	  //checking if RotationMatrix
	  Matrix AfterT = RotationMatrix.Multiply(COM);
	  //for(int i=0;i<Natoms;i++){
	  //  printf("%-2s %10.6f  %10.6f  %10.6f\n",
	  //	   GetAtom(i).GetSymbol().c_str(), AfterT(0,i),AfterT(1,i),AfterT(2,i));
	  //}

	  int i =0;//atoms on current monomer
	  DoesNotMatch=0;//used to identify whether the rotation was correct
	  while(i<Natoms && DoesNotMatch==0){
	    int j = 0;//atoms on unique molecule
	    bool match=0;//match for each atoms
	    while(j<Natoms && match==0){
	      if(GetAtom(i).GetSymbol() == UniqueMon.GetAtom(j).GetSymbol()){
		double x_diff = fabs(AfterT(0,i)-UniqueCOM(0,j));
		double y_diff = fabs(AfterT(1,i)-UniqueCOM(1,j));
		double z_diff = fabs(AfterT(2,i)-UniqueCOM(2,j));
		//printf("i=%i and j=%i\n",i,j);
		//printf("differance of atom %i and atom %i is %f %f %f\n",
		//	   GetAtom(i).GetAtomIndex(),UniqueMon.GetAtom(j).GetAtomIndex(), x_diff,y_diff,z_diff);
		//printf("Atom %f %f %f,Unique %f %f %f\n",
	       	//   AfterT(0,i),AfterT(1,i),AfterT(2,i),COM(0,j),COM(1,j),COM(2,j));
		if(x_diff<tolerance && y_diff<tolerance && z_diff<tolerance){
		  match=1;//identical atom found
		  GetAtom(i).SetSymmetricalAtom(UniqueMon.GetAtom(j).GetGlobalIndex());
		  // printf("atom %i is idential to unique atom %i\n",
		  //	   GetAtom(i).GetGlobalIndex(), UniqueMon.GetAtom(j).GetSymmetricalAtom());
		}
	      }
	      j++;
	    }
	    if(match==0){
	      DoesNotMatch=1;
	    }
	    i++;	       
	  }
	}
	second++;
    }while(second<Natoms && DoesNotMatch==1);
    first++;
  }while(first<Natoms && DoesNotMatch==1);

  if(DoesNotMatch==0){
    SetRotationMatrix(RotationMatrix);
    return 1;
  }else
    
    return 0;
}

//OUTDATED, NO LONGER USED.
bool Monomer::MappingMatrixForLinear(Monomer& UniqueMon, Matrix COM){
  Matrix RotationMatrix(3,3);
  RotationMatrix.Set();
  
  double tolerance = Params::Parameters().GetSymmetryTolerance();//determines how close the atoms have to be
  //Center of mass unique monomer
  Matrix UniqueCOM(3,Natoms);
  UniqueCOM.Set();
  //printf("Unique atom\n");
  for(int i=0;i<Natoms;i++){
    UniqueCOM(0,i) = UniqueMon.GetAtom(i).GetPosition(0)- UniqueMon.GetCenterOfMass(0);// x center of mass
    UniqueCOM(1,i) = UniqueMon.GetAtom(i).GetPosition(1)-UniqueMon.GetCenterOfMass(1);// y center of mass
    UniqueCOM(2,i) = UniqueMon.GetAtom(i).GetPosition(2)-UniqueMon.GetCenterOfMass(2);// z center of mass
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //   UniqueMon.GetAtom(i).GetSymbol().c_str(), UniqueCOM(0,i),UniqueCOM(1,i),UniqueCOM(2,i));
  }

  bool DoesNotMatch = 0;
  int iAtom = 0;//atoms of the molecule symmetrical to the unique molecule 
  do{//loop over iAtom
    int jAtom = 0;//atoms of the unique molecule
    do{//loop over jAtom

      //Reseting DoesNotMatch to zero
      DoesNotMatch = 0;
      bool skip = 0;
      // printf("Atom 1 = %s Atom 2 = %s\n", GetSymbol(iAtom).c_str(), UniqueMon.GetSymbol(jAtom).c_str());
      if(GetSymbol(iAtom) != UniqueMon.GetSymbol(jAtom))
	skip=1;
      //If the distance of the center of mass is zero, skip
      double length1 = COM.GetColumnVector(iAtom).Norm();
      double length2 = UniqueCOM.GetColumnVector(jAtom).Norm();
      // printf("length1=%f length2=%f\n", length1,length2);
      if(fabs(length1)<tolerance){
	 skip=1;
	 //printf("length1=0\n");
      }
      if(fabs(length2)<tolerance){
	skip=1;
	//printf("length2=0\n");
      }
      if(fabs(length1-length2)>tolerance){
	skip=1;
	//printf("length1 does not equal length2\n");
      }
      if(skip==1){
	DoesNotMatch=1;
	//printf("DoesNotMatch =1\n");
      }else{
	//Use to find the diagonal of the rotational matrix using the division of the coordinates of the two atoms
	double dX=1;
	double dY=1;
	double dZ=1;
	if(COM(0,iAtom)!=0 && COM(0,jAtom)!=0)
	  dX = UniqueCOM(0,jAtom)/COM(0,iAtom);
	if(COM(1,iAtom)!=0 && COM(1,jAtom)!=0)
	  dY = UniqueCOM(1,jAtom)/COM(1,iAtom);
	if(COM(2,iAtom)!=0 && COM(2,jAtom)!=0)
	   dZ = UniqueCOM(2,jAtom)/COM(2,iAtom);

	//COM.Print("COM");
	//printf("x1=%f y1=%f z1=%f\n",COM(0,iAtom),COM(1,iAtom),COM(2,iAtom));
	//printf("x2=%f y2=%f z2=%f\n",UniqueCOM(0,jAtom),UniqueCOM(1,jAtom),UniqueCOM(2,jAtom));
	//printf("dx =%f dy=%f dz=%f\n", dX,dY,dZ);

	RotationMatrix(0,0) = dX;
	RotationMatrix(1,1) = dY;
	RotationMatrix(2,2) = dZ;
	RotationMatrix.Print("Linear Rotation");
	Matrix AfterT = RotationMatrix.Multiply(COM);
	  
	int i=0;//atoms on first molecule
	DoesNotMatch=0;//used to identify whether the rotation was correct
	while(i<Natoms && DoesNotMatch==0){
	  int j = 0;//atoms on unique molecule
	  bool match=0;//match for each atoms
	  while(j<Natoms && match==0){
	    if(GetAtom(i).GetSymbol() == UniqueMon.GetAtom(j).GetSymbol()){
	      double x_diff = fabs(AfterT(0,i)-UniqueCOM(0,j));
	      double y_diff = fabs(AfterT(1,i)-UniqueCOM(1,j));
	      double z_diff = fabs(AfterT(2,i)-UniqueCOM(2,j));
	      //printf("i=%i and j=%i\n",i,j);
	      // printf("differance of atom %i and atom %i is %f %f %f\n",
	      //	   GetAtom(i).GetAtomIndex(),UniqueMon.GetAtom(j).GetAtomIndex(), x_diff,y_diff,z_diff);
	      //printf("Atom %f %f %f,Unique %f %f %f\n",
	      //	   AfterT(0,i),AfterT(1,i),AfterT(2,i),UniqueCOM(0,j),UniqueCOM(1,j),UniqueCOM(2,j));
	      if(x_diff<0.0001 && y_diff<0.0001 && z_diff<0.0001){
		match=1;//identical atom found
		GetAtom(i).SetSymmetricalAtom(UniqueMon.GetAtom(j).GetGlobalIndex());
		//printf("atom %-2s is identical to unique atom %-2s\n",
		//GetAtom(i).GetSymbol().c_str(), UniqueMon.GetAtom(j).GetSymbol().c_str());
	      }
	    }
	    j++;
	  }
	  if(match==0){
	    DoesNotMatch=1;
	  }
	  i++;	       
	}
      }
      jAtom++;
    }while(DoesNotMatch==1 && jAtom<Natoms);
    iAtom++;
  }while(DoesNotMatch==1 && iAtom<Natoms);

  if(DoesNotMatch==0){
    SetRotationMatrix(RotationMatrix);
    return 1;
  }else
    return 0;
}

// Print out the Cartesian coordinates in Tinker xyz format
void Monomer::PrintTinkerCartesian(int shift, bool monomer, FILE *outfile) {
  
  if (monomer)
    fprintf(outfile,"%d  (Monomer %d)\n",Natoms, Monomer_index); 

  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintTinkerCartesian(shift,outfile);
  }
}

// Print out Cartesian coordinates different than the assign coordinates Tinker xyz format
void Monomer::PrintTinkerCartesian(Vector coords,int shift,bool monomer, FILE *outfile) {

  if (monomer)
    fprintf(outfile,"%d  (Monomer %d)\n",Natoms, Monomer_index); 

  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintTinkerCartesian(coords[3*i],coords[3*i+1],coords[3*i+2],shift,outfile);
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

// Wrapper for printing MM gradient
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

void Monomer::ShiftCoords(Vector MonomerRModeVector) { // Watit

	Vector RModeVector = MonomerRModeVector;
	Vector old_coord(3); Vector new_coord(3);
	for (int i=0;i<Natoms;i++) {
		for (int dim=0;dim<3;dim++) {
			old_coord[dim] = GetAtom(i).GetCoordinate(dim);
			new_coord[dim] = old_coord[dim]+ RModeVector[dim+i*3];
		}
		Atom_List[i].SetPosition(new_coord);
	}
}

void Monomer::FractionalTranslate(int x,int y,int z){

  //printf("m%i\n",Monomer_index);

  //translate the fractional coordiantes by the number of unit cells over.
  Vector new_coord(3);
  for(int i=0;i<Natoms;i++){
    new_coord[0] = GetAtom(i).GetFractionalPosition(0) + double (x);
    new_coord[1] = GetAtom(i).GetFractionalPosition(1) + double (y);
    new_coord[2] = GetAtom(i).GetFractionalPosition(2) + double (z);
    Atom_List[i].SetFractionalPosition(new_coord);
    //printf("%i %f %f %f\n",i,new_coord[0],new_coord[1],new_coord[2]);
  }

  //alter the symmetry translational vector
  //untested to see if it produces the correct translational vector
  Fractional_Translate[0] += double(x);
  Fractional_Translate[1] += double(y);
  Fractional_Translate[2] += double(z);
  //Fractional_Translate.Print("Fractional_Translate");

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
    Monomer_type = other.Monomer_type;
    Natoms = other.Natoms;
    unique_atoms = other.unique_atoms;

    sym_fac = other.sym_fac;
    SymMon = other.SymMon;
    Symmetry_Rotation = other.Symmetry_Rotation;

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      Fractional_Rotation = other.Fractional_Rotation;
      Fractional_Translate = other.Fractional_Translate;
    }

    // Copy over Atom_List and Sym_Atom
    // First delete old list, just in case. It should either be
    // allocated or NULL, so no harm in using delete
    delete [] Atom_List; 
    Atom_List = new Atom[Natoms];
    delete [] Sym_Atom;
    Sym_Atom = new int[Natoms];
    for (int i=0;i<Natoms;i++){
      Atom_List[i] = other.Atom_List[i];
      Sym_Atom[i] = other.Sym_Atom[i];
    }
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

    QM_Hess_Init = other.QM_Hess_Init;         
    if (QM_Hess_Init) {    
      Hess_QM = other.Hess_QM;      
    }

    MM_Hess_Init = other.MM_Hess_Init;
    if (MM_Hess_Init) {
      Hess_MM = other.Hess_MM;
    }
    RotationAngle = other.RotationAngle;
    RotationVector = other.RotationVector;
  }
  return *this;
}

Matrix Monomer::ReadHessian(string path, int type) {

  Matrix monomer_hess( 3*Natoms, 3*Natoms );

  // Set up the filename, with the full path.  File is e.g. 'm1.freq'
  string filename = path + "/" + mN_;
  if (type == 2) // Tinker MM job                  
    filename += ".freq";
  else if(Params::Parameters().GetQMType()==3) // G09
    filename += ".log";      
  else // Qchem or Molpro
    filename += ".out";

  // Open the freq file               
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadHessian : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }                            

  if (type == 2) { // look in the tinker .freq file
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      // Search for the SCF hessian
      if ( line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):" ) {
        //getline(infile,line); // throw away header line
        for (int i=0;i<3*Natoms;i++) {
          for (int j=0;j<3*Natoms;j++) {
            getline(infile,line);
            istringstream iss(line);
            string tmp1,tmp2;
            iss >> tmp1; // throw away the matrix row index
            iss >> tmp2; // throw away the matrix column index
            iss >> monomer_hess.Element(i,j); // Store the hessian elements
          }
        }
        break;
      }
    }
  }

  else if(Params::Parameters().GetQMType() == 1){ //look in the qchem output
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      // Search for the SCF hessian

      if ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) {
        for (int k=0;k<(3*Natoms-3*Natoms%6)/6;k++) {
          getline(infile,line); // throw away header line           
          //while (line.substr(0,15)!=" Gradient time:") {

          for (int i=0;i<3*Natoms;i++) {
            getline(infile,line);
            istringstream iss(line);
            string tmp;
            iss >> tmp; // throw away the atom index
            for (int j=6*k;j<6*k+6;j++) {            
              iss >> monomer_hess.Element(i,j); // Store the hessian elements
            }
          }
        }
        getline(infile,line); // check if this line signals end of hessian print or not
        if (line.substr(0,15) != " Gradient time:") {
          for (int i=0;i<3*Natoms;i++) {
            getline(infile,line);
            istringstream iss(line);
            string tmp;
            iss >> tmp; // throw away the atom index
            for (int j=(3*Natoms-3*Natoms%6);j<3*Natoms;j++) {
              iss >> monomer_hess.Element(i,j); // Store the hessian elements
            }
          }
          break;
        }
        else {
          break;
        }
      }
    }
  }
  else if(Params::Parameters().GetQMType() == 2){
    if(Params::Parameters().DoCBS()){
      double basis1 = double(Params::Parameters().CBSBasis1());
      double basis2 = double(Params::Parameters().CBSBasis2());
      double ExpTerm = exp(1.54*(basis2-basis1));
      string line;
      int grad_index=0;
     while ( !infile.eof() ) {
       getline(infile,line);
       if ( line.substr(0,16)==" Force Constants"){
	 grad_index++;
	 for (int k=0;k<=(3*Natoms-3*Natoms%5)/5;k++) {  
	   getline(infile,line); // throw away header line
	   for (int i=5*k;i<3*Natoms;i++){
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=5*k;j<5*k+5;j++) {   
		double entry = 0.0;   
		if(i >= j){
		  //HF in the first basis
		  if(grad_index==1)
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  //MP2 in the first basis
		  else if(grad_index==2)
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  else if(grad_index==3)
		   entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  else if(grad_index==4)
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  monomer_hess.Element(i,j) = entry;// Store the hessian elements
		}
	      }
	   }
	 }
       }
     }
    }
    else{
      string line;
      while ( !infile.eof() ) {
	getline(infile,line);
	if ( line.substr(0,16)==" Force Constants"){
	  for (int k=0;k<=(3*Natoms-3*Natoms%5)/5;k++) {  
	    getline(infile,line); // throw away header line
	    for (int i=5*k;i<3*Natoms;i++){
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=5*k;j<5*k+5;j++) {      
		if(i >= j)
		  iss >> monomer_hess.Element(i,j); // Store the hessian elements
	      }
	    }
	  }
	}
      }
    }
    //copying lower triangle
    for(int i = 0;i<3*Natoms;i++)
      for(int j = i; j<3*Natoms;j++)
	monomer_hess.Element(i,j) = monomer_hess.Element(j,i);
  }
  else if(Params::Parameters().GetQMType()==3){ // G09 Watit
    string line;
    while (!infile.eof()) {
      getline(infile,line);
      if(line.substr(0,17)==" Hessian entering") {
	int maxk=(3*Natoms-3*Natoms%5)/5;
	if(3*Natoms%5!=0) {
	  maxk+=1;
	}
	for(int k=0;k<maxk; k++) {
	  getline(infile,line);
	  for(int i=5*k;i<3*Natoms;i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp;
	    for(int l=0;l<5;l++) {
	      int j=5*k+l;
	      stringstream ss;
	      string tmp2;
	      iss >> tmp2;
	      if(tmp2.length()!=0&&tmp2.length()<14) {
		replace(tmp2.begin(),tmp2.end(),'D','E');
		ss << tmp2;
		ss >> monomer_hess.Element(i,j);
		//monomer_hess.Element(i,j)=-monomer_hess.Element(i,j);
		monomer_hess.Element(j,i)=monomer_hess.Element(i,j);
	      }
	    }
	  }
        }
      }
    }				
  } 
  else{
    printf("Monomer::ReadHessian() - QM type = %d is unknown\n",Params::Parameters().GetQMType());
    exit(1);
  }


  //Including CCSD(T) correction, molpro only
  if(Params::Parameters().DoCCSDTCorrection() && type == 1 ) {
    //monomer_hess.PrintHessian("Before CCSD(T)");
    Matrix CCSDT_hess = ReadFiniteDifferenceCCSDTHessian(path);
    for(int i=0;i<3*Natoms;i++){
      for(int j=0;j<3*Natoms;j++){
	monomer_hess(i,j) += CCSDT_hess(i,j);
      }
    }
    //monomer_hess.PrintHessian("monomer_hess");
  }

  infile.close();  // Close the file
/*
  for (int i=0;i<3*Natoms;i++) {  
    for (int j=0;j<3*Natoms;j++) {
      cout << monomer_hess.Element(i,j) << "\t";
    }
    cout << "\n";
  }
*/

  return monomer_hess;
}

Matrix Monomer::ReadFiniteDifferenceMolProHessian(string path){


  Matrix monomer_hess( 3*Natoms, 3*Natoms );

  //directory for finite difference steps
  path += "/" + mN_;

  for(int i=0;i<3*Natoms;i++){

    Vector GradPlus(3*Natoms);
    Vector GradMinus(3*Natoms);

    // The plus ouput file; 
    char num[10];
    sprintf(num,"%d",i);
    string filename = path + "/" + mN_ + "+" + num + ".out";

    // Open the .out file
    ifstream infile;                 
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadFiniteDifferenceMolProHessian : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }

    int grad_index = 0;
    string line;

    //Using the CBS limit
    if(Params::Parameters().DoCBS()){
      double basis1 = double(Params::Parameters().CBSBasis1());
      double basis2 = double(Params::Parameters().CBSBasis2());
      double ExpTerm = exp(1.54*(basis2-basis1));

      while ( !infile.eof() ) {
	getline(infile,line);
	if(line.substr(0,60) ==" Atom          dE/dx               dE/dy               dE/dz"){ 
	    grad_index++;
	  if(!Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	    getline(infile,line);// throw away header line
	  for (int j=0;j<Natoms;j++) {
	    double entry;
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> entry;
	      //HF in the first basis
	      if(grad_index == 1){
		entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      }
	      //MP2 in the first basis
	      else if(grad_index == 2){
		entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      }
	      //HF in the second basis
	      else if(grad_index == 3){
		entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      } 
	      //MP2 in the second basis
	      else if(grad_index == 4){
		entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      }

	      GradPlus[3*j+k] += entry;
	      // Store the gradient elements
	    }
	  }
	  //printf("grad_index = %i\n",grad_index);
	  //break;
	}	
      }

      if (!(grad_index == 4)) {
	printf("Monomer::ReadFiniteDifferenceMolProHessian At Least One Gradient for CBS Not Found or Incorrect '%s'\n",
	       filename.c_str());
	printf("grad_index = %i\n",grad_index);
	exit(1);
      }
    }else{
      while ( !infile.eof() ) {
	getline(infile,line);
	if(line.substr(0,60)==" Atom          dE/dx               dE/dy               dE/dz"){
	  if(!Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	    getline(infile,line);// throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> GradPlus[3*j+k];
	      // Store the gradient elements
	    }
	   }
	  grad_index++;
	  //break;
	}
      }
    }

     infile.close();
     //GradPlus.Print("GradPlus");

    // The plus output file;  
    filename = path + "/" + mN_ + "-" + num + ".out";

    // Open the .out file                 
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadFiniteDifferenceMolProHessian : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }

    grad_index = 0;

    //Using the CBS limit
    if(Params::Parameters().DoCBS()){
      double basis1 = double(Params::Parameters().CBSBasis1());
      double basis2 = double(Params::Parameters().CBSBasis2());
      double ExpTerm = exp(1.54*(basis2-basis1));

      while ( !infile.eof() ) {
	getline(infile,line);
	if(line.substr(0,60)==" Atom          dE/dx               dE/dy               dE/dz"){ 
	  grad_index++;
	  if(!Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	    getline(infile,line);// throw away header line
	  for (int j=0;j<Natoms;j++) {
	    double entry;
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> entry;
	      //HF in the first basis
	      if(grad_index == 1){
		entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      }
	      //MP2 in the first basis
	      else if(grad_index == 2){
		entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
	      }
	      //HF in the second basis
	      else if(grad_index == 3){
		entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      } 
	      //MP2 in the second basis
	      else if(grad_index == 4){
		entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
	      }

	      GradMinus[3*j+k] += entry;
	      // Store the gradient elements
	    }
	  }

	  //break;
	}
	//printf("grad_index = %i\n",grad_index);
	
      }
      
      if (!(grad_index == 4)) {
	printf("Monomer::ReadFiniteDifferenceMolProHessian At Least One Gradient for CBS Not Found or Incorrect '%s'\n",
	       filename.c_str());
	printf("grad_index = %i\n",grad_index);
	exit(1);
      }
      
    }else{
      while ( !infile.eof() ) {
	getline(infile,line);
	if(line.substr(0,60)==" Atom          dE/dx               dE/dy               dE/dz"){
	  if(!Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	    getline(infile,line);// throw away header line
	  
	  for (int j=0;j<Natoms;j++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int k=0;k<3;k++) {
	      iss >> GradMinus[3*j+k];
	      // Store the gradient elements
	    }
	   }
	  grad_index++;
	  //break;
	}
      }
    }
    infile.close();

    //GradMinus.Print("GradMinus");
    GradPlus -= GradMinus;
    GradPlus.Scale(1/(0.002*AngToBohr));
    //GradPlus.Print("Hess Column");
    monomer_hess.SetColumnVector(GradPlus,i);
  }

  //Including CCSD(T) correction, molpro only
  if(Params::Parameters().DoCCSDTCorrection() ) {
    //monomer_hess.PrintHessian("Before CCSD(T)");
    Matrix CCSDT_hess = ReadFiniteDifferenceCCSDTHessian(path);
    for(int i=0;i<3*Natoms;i++)
      for(int j=0;j<3*Natoms;j++)
	monomer_hess(i,j) += CCSDT_hess(i,j);
    //monomer_hess.PrintHessian("monomer_hess");
  }

  //printf("m%i\n",GetIndex());
  //monomer_hess.PrintHessian("monomer_hess");



  return monomer_hess;
}

//CCSD(T) correction has only been implemented for molpro
Matrix Monomer::ReadFiniteDifferenceCCSDTHessian(string path){

  Matrix monomer_hess(3*Natoms,3*Natoms);

  path = path + ".CCSDT";

  for(int i=0;i<3*Natoms;i++){
    
    int grad_index = 0;
    Vector GradPlus(3*Natoms);
    Vector GradMinus(3*Natoms);
    
    // The plus ouput file; 
    char num[10];
    sprintf(num,"%d",i);
    string filename = path + "/" + mN_ + "+" + num + ".out";
    //printf("filename = %s\n",filename.c_str());
    
    // Open the .out file
    ifstream infile;                 
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("ReadFiniteDifferenceCCSDTHessian : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
    
    while ( !infile.eof() ) {
      string line;
      getline(infile,line);
      if(line.substr(0,60) ==" Atom          dE/dx               dE/dy               dE/dz"){	
	grad_index++;
	if(grad_index==1){
	  getline(infile,line);// throw away header line
	}
	for (int j=0;j<Natoms;j++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  for (int k=0;k<3;k++) {
	    double entry;
	    iss >> entry;
	    //change the sign of the MP2 part
	    if(grad_index == 1)
	      entry *= -1;
	    GradPlus[3*j+k] += entry;//gradients will be rotated if needed due to symmetry
	  }	  
	}
      } 
    }   

    infile.close();
    //GradPlus.PrintGradient("GradPlus");

      if (!(grad_index == 2)) {
	printf("Monomer::ReadFiniteDifferenceCCSDTHessian() At Least One Gradient for CBS Not Found or Incorrect '%s'\n",
	       filename.c_str());
	printf("grad_index = %i\n",grad_index);
	exit(1);
      }


    filename = path + "/" + mN_ + "-" + num + ".out";
    //printf("filename = %s\n",filename.c_str());
    
    // Open the .out file                 
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Monomer::ReadFiniteDifferenceCCSDTHessian() : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
    
    grad_index=0;

    while ( !infile.eof() ) {
      string line;
      getline(infile,line);
      if(line.substr(0,60) ==" Atom          dE/dx               dE/dy               dE/dz"){	
	grad_index++;
	if(grad_index==1)
	  getline(infile,line);// throw away header line
	
	for (int j=0;j<Natoms;j++) {
	  getline(infile,line);
	  istringstream iss(line);
	  string tmp;
	  iss >> tmp; // throw away the atom index
	  for (int k=0;k<3;k++) {
	    double entry;
	    iss >> entry;
	    //change the sign of the MP2 part
	    if(grad_index == 1)
	      entry *= -1;
	    GradMinus[3*j+k] += entry;//gradients will be rotated if needed due to symmetry 
	  }	  
	}
      }

    }
    if (!(grad_index == 2)) {
      printf("Monomer::ReadFiniteDifferenceCCSDTHessian() At Least One Gradient for CBS Not Found or Incorrect '%s'\n",
	     filename.c_str());
      printf("grad_index = %i\n",grad_index);
      exit(1);
    }
    
    infile.close();
    //GradMinus.PrintGradient("GradMinus");
    GradPlus -= GradMinus;
    GradPlus.Scale(1/(0.002*AngToBohr));
    //GradPlus.Print("Hess Column");
    
    monomer_hess.SetColumnVector(GradPlus,i);
  }
  //monomer_hess.PrintHessian("CCSD(T)_hess");
  return monomer_hess;

}

//Copying and rotating hessian for symmetrically equivalent monomers
Matrix Monomer::SetSymmetricalHessian(int type,Monomer& Sym_Mon){

  Matrix monomer_hess( 3*Natoms, 3*Natoms);
  Matrix Sym_hess( 3*Natoms, 3*Natoms );
  Matrix Rot = Symmetry_Rotation;
  Rot.Transpose();
  
  if (type == 1) {//qchem or molpro
    Sym_hess = Sym_Mon.GetQMHessian();
    QM_Hess_Init = 1;    
  }
  if (type == 2) {// tinker
    Sym_hess = Sym_Mon.GetMMHessian();
    MM_Hess_Init = 1;    
  }
  
  //Element of the Hessian under Symmetry H12 = R13*H34*R42; R42 is the transpose of R24
  for(int i=0; i<Natoms;i++){
    for(int j=0; j<Natoms;j++){  
      for(int xyz1=0;xyz1<3;xyz1++)
	for(int xyz2=0;xyz2<3;xyz2++)
	  for(int xyz3=0;xyz3<3;xyz3++)
	    for(int xyz4=0;xyz4<3;xyz4++){
	      //printf(" matrix entry (%i,%i) by (%i,%i) \n",3*i+xyz1,3*j+xyz2,3*i+xyz3,3*j+xyz4);
	      //printf(" 3*i = %i 3*j = %i xyz1 = %i xyz2 = %i xyz3 = %i xyz4 =%i\n",3*i,3*j,xyz1,xyz2,xyz3,xyz4); 
	      //monomer_hess(3*Symi+xyz1,3*Symj+xyz2) += Rot(xyz1,xyz3)*Sym_hess(3*i+xyz3,3*j+xyz4)*Rot(xyz2,xyz4);  
              monomer_hess(3*i+xyz1,3*j+xyz2) += Rot(xyz1,xyz3)*Sym_hess(3*Sym_Atom[i]+xyz3,3*Sym_Atom[j]+xyz4)*Rot(xyz2,xyz4); 
	    }
      
    }
  }

  //printf("m%i\n",Monomer_index);
  //monomer_hess.PrintHessian("Rotated Hess");
  return monomer_hess;
}

void Monomer::SetMMHessian() {
  if (Params::Parameters().GetMMType() == 1) { // Tinker
    string path = Params::Parameters().GetHessianMMPath();

   //Path of the quasiharmonic calculations
   if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Hess_MM = ReadHessian(path,2); 
    MM_Hess_Init = 1;  
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    // AIFF hessian handled elsewhere ??
    if (!MM_Hess_Init) Hess_MM.Initialize(3*Natoms); // create empty array
    MM_Hess_Init = 1;         
  }
  else if (Params::Parameters().GetMMType() == 3) { // QChem
    string path = Params::Parameters().GetHessianMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Hess_MM = ReadHessian(path,1);
    MM_Hess_Init = 1;
  }
  else {
    printf("Monomer::SetMMHessian() - MM type = %d is unknown\n",
           Params::Parameters().GetMMType());
    exit(1);
  }

}


// Get the QM Hessian, wrapper routine
void Monomer::SetQMHessian() {           
  string path = Params::Parameters().GetHessianQMPath();     
/*  if (Params::Parameters().TinkerDebug() ) {
    Hess_QM = ReadHessian(path,2);
    QM_Hess_Init = 1;
  }
  else {
*/

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //MolPro hessian found by finite difference
  if(Params::Parameters().GetQMType() == 2 && !Params::Parameters().SingleFileMonomerHess()){
      Hess_QM = ReadFiniteDifferenceMolProHessian(path);
  }else
    Hess_QM = ReadHessian(path,1);
  QM_Hess_Init = 1;      
//  }
}

// Get the QM Hessian, copied from another Monomer 
void Monomer::SetQMHessian(Monomer& other) {
  if (other.QM_Hess_Init) {
    Hess_QM = other.GetQMHessian();
    QM_Hess_Init = 1;
  }
}

// Get the MM Hessian, copied from another Monomer
void Monomer::SetMMHessian(Monomer& other) {
  if (other.MM_Hess_Init) {
    Hess_MM = other.GetMMHessian();
    MM_Hess_Init = 1;
  }
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
      // Get the van der Waals diameters
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");

      bool DoTangToenniesDamping = true;
      double F6ij;
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

	// ATM 3-body dispersion contribution
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
      }
    }
  }
  
  return E3b_disp;
}

// JDH use GDMA inplace of CAMCASP for speed
void Monomer::CreateGDMAJob() {
  /*
    Running GDMA and creating the multipole moment files is a multiple step process:
    1) Create and run g09 files using g09 header, but saving a .chk file
    2) using 'formchk' gaussian utility to create a formated check point file
    3) create and fun a gdma input file that reads in data from the formated checkpoint file
    4) Run GDMA to create the .mom file
    5) reformat the .mom file to match the format necessary for HMBI
   */

  // Create G09 input file:
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".gdma.com";
  string check_file = mN_ + ".gdma.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateGDMAJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for GDMA run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }
  
  fclose(job);

  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  cmd = "sed -i s/charge// ";
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

  // Create the gdma input file:
  filename = path + "/" + mN_ + ".gdma.in";
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateGDMAJob : cannot open gdma input file '%s'\n",filename.c_str());
    exit(1);
  }
  
  fprintf(job,"Title \"GDMA input file for monomer %d \" \n", Monomer_index);
  check_file = mN_ + ".gdma.fchk";
  fprintf(job,"File %s\n", check_file.c_str() );
  fprintf(job,"\n");
  fprintf(job,"Angstrom\n");
  fprintf(job,"Multipoles\n");
  fprintf(job,"  switch 4\n");
  fprintf(job,"  Limit 4\n");
  fprintf(job,"  Limit 1 H\n");
  fprintf(job,"  Radius H 0.325\n");
  fprintf(job,"  Punch m%d.mom\n", Monomer_index);
  fprintf(job,"Start\n");
  fprintf(job,"\n");
  fprintf(job,"Finish\n");

  fclose(job);

}


// JDH: takes the .mom file in the QM path created with GDMA
// and creates a modified .mom file in the MM path 
// that HMBI can read:
void Monomer::ModifyGDMAMomFile(){
  string qm_path = Params::Parameters().GetQMPath();
  string mm_path = Params::Parameters().GetMMPath();
  
  string in_file = qm_path + "/" + mN_ + ".mom";
  string out_file = mm_path + "/" + mN_ + ".mom";

  // Open the outfile for writing...
  FILE *outfile;
  if (( outfile = fopen(out_file.c_str(),"w"))==NULL) {
    printf("Monomer::ModifyGDMAMomFile : cannot open mom outfile '%s'\n",out_file.c_str());
    exit(1);
  }

  // Open the infile for reading...
  ifstream infile;
  infile.open( in_file.c_str() );
  if ( !infile.is_open() ) {
    printf("Monomer::ModifyGDMAMomFile : Cannot open mom input file '%s'\n",
	   out_file.c_str());
    exit(1);
  }

  // Print the first few lines to the new .mom files
  fprintf(outfile,"! monomer\n" );
  fprintf(outfile,"!  Basis: (need to read out of g09.gdma file)\n");
  fprintf(outfile,"\n");

  // Now start reading through the gdma .mom file and pull out the stuff we need...
  // Remember: I suck at C++ file I/O
  int iatom = 0;
  string line, newline;
  size_t start_pos;
  int numH = 0;
  int numC = 0;
  int numS = 0;
  int numN = 0;
  int numO = 0;
  int numCl = 0;
  int numF = 0;
  getline(infile,line);
  getline(infile,line);
  getline(infile,line);
  getline(infile,line); // Throw away the first 4 lines...
  // Turns out we don't number the atoms individually - *sigh*
  int num = 0;
  while ( !infile.eof() ) {
    getline(infile,line);
    //printf("%s\n", line.c_str() );
    if ( line.length() > 0 ) {
      string match = line.substr(0,2);
      if ( match  == "C " ) {
	//numC++;
	num++;
	string oldlabel = "C ";
	ostringstream tmp;
	tmp << "C";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}

      } else if ( match == "H " ) {
	num++;
	string oldlabel = "H ";
	ostringstream tmp;
	tmp << "H";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=3;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
	
      } else if ( match == "S " ) {
	num++;
	string oldlabel = "S ";
	ostringstream tmp;
	tmp << "S";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
      } else if ( match == "N " ){ 
	num++;
	string oldlabel = "N ";
	ostringstream tmp;
	tmp << "N";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
      } else if (match == "O ") {
	num++;
	string oldlabel = "O ";
	ostringstream tmp;
	tmp << "O";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
      } else if (match == "Cl ") {
	num++;
	string oldlabel = "Cl ";
	ostringstream tmp;
	tmp << "Cl";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
      } else if (match == "F ") {
	num++;
	string oldlabel = "F ";
	ostringstream tmp;
	tmp << "F";
	tmp << num;
	string newlabel = tmp.str();

	start_pos = line.find(oldlabel);
	line.replace(start_pos, oldlabel.length(), newlabel); 

	fprintf(outfile,"%s\n", line.c_str());

	for (int i=1;i<=8;i++) {
	  getline(infile,line);
	  fprintf(outfile,"%s\n",line.c_str() );
	}
      } else {
	printf("ERROR Monomer::ModifyGDMAMomFile Sorry, I didn't include this atom type, but don't panic just go to the function and add it in..\n");
	exit(1);
      }
      
    }
  }
  
  infile.close();
  fclose(outfile);
}

void Monomer::PrintEmbeddingCharges(FILE* outfile) {
  for (int i=0;i<Natoms;i++) {
    Atom_List[i].PrintEmbeddingCharges(outfile);
  }
}

// JDH + Watit
void Monomer::CreateG09Job( Monomer Monomers[], int NMon, bool MM_job ) {

  //printf("ERROR: Monomer::CreateG09Job( Monomer Monomers[], int NMon )  funtion not finished: complete basis section, and make sure ee is treated properly\n");
  //exit(1);
  
  string path;
  if (MM_job && !Params::Parameters().DoFreq() ) {
    path = Params::Parameters().GetMMPath();
  }
  if (!MM_job && !Params::Parameters().DoFreq() )  {
    path = Params::Parameters().GetQMPath();
  }
  if (MM_job && Params::Parameters().DoFreq() ) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if (!MM_job && Params::Parameters().DoFreq() )  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/" + mN_ + ".com"; 


  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  string gaussianHeader; // Watit
  gaussianHeader = Params::Parameters().GetGaussianHeader();
  if(Params::Parameters().DoForces()) {
    gaussianHeader += " force\n";
  }
  if(Params::Parameters().DoFreq()&&!Params::Parameters().DoRaman()) {
    gaussianHeader += " freq=noraman IOp(7/33=3)\n";
  }
  if(Params::Parameters().DoRaman()&&!Params::Parameters().DoForces()) {
    gaussianHeader += " freq=raman IOp(7/33=3)\n";
  }

  fprintf(job,"%s\n", gaussianHeader.c_str() );
  fprintf(job,"Monomer %d\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  // Optionally print $external_charges section
  if ( Params::Parameters().UseElectrostaticEmbedding() ) {

    fprintf(job,"\n"); // use a blank line to separate geom from charge list
    
    // Print Charges from the other monomers in the unit cell (or hmbi input file)
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }
  }

  // Optionally print the mixed basis section...
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);
      if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");

    }

  } // end full custom basis section:

  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      //printf("DEBUG: Basis Region: %d\n", GetAtom(iatom).GetMixedBasisRegion() ); 
      //printf("DEBUG: Basis 1: %s\n", Params::Parameters().GetNMRMixedBasisLevel1() ); fflush(stdout);
      fprintf(job,"%d 0\n", iatom+1);

      if ( Params::Parameters().CustomBasis() == 0 ) {
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	// Only use custom basis sets on asymmetric monomers:
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  // Awkward... but loop through all of the custom basis sets
	  // change this to a numerical reference later:
	  if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } 


      fprintf(job,"****\n");
    }

  }


  fprintf(job,"\n");
  
  fclose(job);
}



void Monomer::CreateG09Job( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".com"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d\n", Monomer_index);
  fprintf(job,"\n");

 
  

  //printf("DEBUG: %s\n", Params::Parameters().GetGaussianHeader().c_str() );


  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);


  // Optionally print $external_charges section
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n"); // use a blank line to separate geom from charge list
    
    // Print Charges from the other monomers in the unit cell (or hmbi input file)
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index  && Monomers[Monomer_index].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }
    
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images;i++) {
	if ( MonomerImages[i].GetUseInEmbedding() && Monomers[Monomer_index].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	  //printf("Printing Charges for MonomerImages[%d] = %d\n", i, MonomerImages[i].GetIndex() ); fflush(stdout);
	  //for (int iatom=0; iatom<MonomerImages[i].GetNumberOfAtoms(); iatom++ ){
	  //printf("Atom = %d  charge = %f\n", iatom,MonomerImages[i].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);
	  //}
	  //MonomerImages[i].PrintQChemEmbeddedDipoles(job);
	  
	  MonomerImages[i].PrintEmbeddingCharges(job);
	}
      }
    } else {
      printf("DEBUG: you called the periodic version of CreateG09NMRJob for a periodic system, but you're not using PBCs!\n");
      exit(1);
    }
    
  }  // end use in NMR charge embeding...
  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);
      if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");

    }

  } // end full custom basis section:


  // Optionally print the mixed basis section...
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      //printf("DEBUG: Basis Region: %d\n", GetAtom(iatom).GetMixedBasisRegion() ); 
      //printf("DEBUG: Basis 1: %s\n", Params::Parameters().GetNMRMixedBasisLevel1() ); fflush(stdout);
      fprintf(job,"%d 0\n", iatom+1);

      if ( Params::Parameters().CustomBasis() == 0 ) {
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	// Only use custom basis sets on asymmetric monomers:
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  // Awkward... but loop through all of the custom basis sets
	  // change this to a numerical reference later:
	  if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } 


      fprintf(job,"****\n");
    }

  }

  //printf("DEBUG: after print the mixed basis to the G09 file...\n");
  //fflush(stdout);

  fprintf(job,"\n");
  
  fclose(job);



}


void Monomer::CreateG09Job( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ) {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".com"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d\n", Monomer_index);
  fprintf(job,"\n");

 
  

  //printf("DEBUG: %s\n", Params::Parameters().GetGaussianHeader().c_str() );


  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);


  // Optionally print $external_charges section
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n"); // use a blank line to separate geom from charge list
    
    // Print Charges from the other monomers in the unit cell (or hmbi input file)
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index  && Monomers[Monomer_index].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }
    
    for (int i=1;i<=NMon_images;i++) {
      if ( MonomerImages[i].GetUseInEmbedding() && Monomers[Monomer_index].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
	//printf("Printing Charges for MonomerImages[%d] = %d\n", i, MonomerImages[i].GetIndex() ); fflush(stdout);
	//for (int iatom=0; iatom<MonomerImages[i].GetNumberOfAtoms(); iatom++ ){
	//printf("Atom = %d  charge = %f\n", iatom,MonomerImages[i].GetAtom(iatom).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);
	//}
	//MonomerImages[i].PrintQChemEmbeddedDipoles(job);
	
	MonomerImages[i].PrintEmbeddingCharges(job);
      }
    }
    
    // Now print the charges on the Ewald Sphere:
    for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
      if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	fprintf(job, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      } else {
	fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      }
    }

    
  }  // end use in NMR charge embeding...
  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);
      if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");

    }

  } // end full custom basis section:


  // Optionally print the mixed basis section...
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      //printf("DEBUG: Basis Region: %d\n", GetAtom(iatom).GetMixedBasisRegion() ); 
      //printf("DEBUG: Basis 1: %s\n", Params::Parameters().GetNMRMixedBasisLevel1() ); fflush(stdout);
      fprintf(job,"%d 0\n", iatom+1);

      if ( Params::Parameters().CustomBasis() == 0 ) {
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	// Only use custom basis sets on asymmetric monomers:
	if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  // Awkward... but loop through all of the custom basis sets
	  // change this to a numerical reference later:
	  if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } 


      fprintf(job,"****\n");
    }

  }

  //printf("DEBUG: after print the mixed basis to the G09 file...\n");
  fflush(stdout);

  fprintf(job,"\n");
  
  fclose(job);



}

void Monomer::CreateG09HirshfeldJob(){
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".Hirshfeld.com";
  //  string check_file = mN_ + ".chelpG.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09HirshfeldJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for ChelpG run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      fprintf(job,"V %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);

  //printf("DEBUG: Monomer %d\n", GetIndex() ); fflush(stdout);

  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("Shouldn't need custom vdw \n");
    exit(1);
  } else if (custom_vdw_radii == false ) {
    cmd = "sed -i s/charge/'pop=Hirshfeld'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}

void Monomer::CreateG09HirshfeldJob(Monomer Monomers[], int NMon ) {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".Hirshfeld.com";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09HirshfeldJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for CHELPG SCE run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // Print Charge Section:
  // Print Charges from the other monomers in the unit cell (or hmbi input file)
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(job);
    }
  }

  fprintf(job,"\n");
  

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }

  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  
  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'charge pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("Shouldn't need custom vdw \n");
    exit(1);
  } else {
    cmd = "sed -i s/charge/'charge pop=Hirshfeld'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}


void Monomer::CreateG09HirshfeldJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".Hirshfeld.com";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09HirshfeldJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for CHELPG SCE run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // Print Charge Section:
  // Print Charges from the other monomers in the unit cell (or hmbi input file)
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(job);
    }
  }

  for (int i=1;i<=NMon_images;i++) {
    if ( MonomerImages[i].GetUseInEmbedding() ) {
      MonomerImages[i].PrintEmbeddingCharges(job);
    }
  }
  fprintf(job,"\n");
  

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }

  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  
  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'charge pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("Shouldn't need custom vdw \n");
    exit(1);
  } else {
    cmd = "sed -i s/charge/'charge pop=Hirshfeld'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}



void Monomer::CreateG09ChelpGJob() {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".chelpG.com";
  //  string check_file = mN_ + ".chelpG.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09ChelpGJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for ChelpG run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      fprintf(job,"V %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }

  }

  fprintf(job,"\n");
  fclose(job);

  //printf("DEBUG: Monomer %d\n", GetIndex() ); fflush(stdout);

  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else if (custom_vdw_radii == false ) {
    cmd = "sed -i s/charge/'pop=ChelpG'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}


void Monomer::CreateG09ChelpGJob(Monomer Monomers[], int NMon) {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".chelpG.com";
  //  string check_file = mN_ + ".chelpG.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09ChelpGJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for CHELPG SCE run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // Print Charge Section:
  // Print Charges from the other monomers in the unit cell (or hmbi input file)
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(job);
    }
  }

  fprintf(job,"\n");
  

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }

  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      fprintf(job,"V %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
    

  }

  fprintf(job,"\n");
  fclose(job);


  
  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'charge pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else {
    cmd = "sed -i s/charge/'charge pop=ChelpG'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}

void Monomer::CreateG09ChelpGJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".chelpG.com";
  //  string check_file = mN_ + ".chelpG.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09ChelpGJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for CHELPG SCE run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // Print Charge Section:
  // Print Charges from the other monomers in the unit cell (or hmbi input file)
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(job);
    }
  }

  for (int i=1;i<=NMon_images;i++) {
    if ( MonomerImages[i].GetUseInEmbedding() ) {
      MonomerImages[i].PrintEmbeddingCharges(job);
    }
  }
  fprintf(job,"\n");
  

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }

  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      fprintf(job,"V %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }


  }

  fprintf(job,"\n");
  fclose(job);


  
  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'charge pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else {
    cmd = "sed -i s/charge/'charge pop=ChelpG'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}

void Monomer::CreateG09ChelpGJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ) {
  // Create G09 input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".chelpG.com";
  //  string check_file = mN_ + ".chelpG.chk";

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09ChelpGJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }

  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Monomer %d Gaussian input for CHELPG SCE run\n", Monomer_index);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"%d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);

  fprintf(job,"\n");

  // Print Charge Section:
  // Print Charges from the other monomers in the unit cell (or hmbi input file)
  for (int i=1;i<=NMon;i++) {
    if (i != Monomer_index) {
      Monomers[i].PrintEmbeddingCharges(job);
    }
  }

  for (int i=1;i<=NMon_images;i++) {
    if ( MonomerImages[i].GetUseInEmbedding() ) {
      MonomerImages[i].PrintEmbeddingCharges(job);
    }
  }

  // Now print the Ewald Charges:
  for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
    fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
  }

  fprintf(job,"\n");
  

  // If using a mixed basis add in the level 1 basis set portion
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%d 0\n", iatom+1);


      if ( Params::Parameters().CustomBasis() > 0) {
	
	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	} else if ( GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	  fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }


      fprintf(job,"****\n");
    }
    fprintf(job,"\n");
  }

  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      fprintf(job,"V %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }

  }

  fprintf(job,"\n");
  fclose(job);


  
  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";
  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'charge pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else {
    cmd = "sed -i s/charge/'charge pop=ChelpG'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());

}


// Create PSI4 Input Files // CSG
void Monomer::CreatePSI4Job(Monomer Monomers[],int NMon) {//GJB for debug

  //Set up the filename, with the full path. File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath();
  string filename = path + "/" + mN_ + ".in";

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreatePSI4Job : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  //input header
  fprintf(job,"# Monomer %d\n",Monomer_index);
  fprintf(job,"memory %i mb\n\n",Params::Parameters().MemoryUse());

  //Print Cartesian Coordinates
  fprintf(job,"molecule monomer {\n");
  fprintf(job,"%i %i\n\n",charge, spin);
  PrintPSI4MonomerCartesian(job,1);
  fprintf(job,"}\n\n");

  //Print REM section 
  fprintf(job,"%s\n", Params::Parameters().GetPSI4Rem().c_str() );

  //fprintf(job,"Emonomer=energy\n\n");
  fclose(job);
}

void Monomer::ReadDaltonNMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();
  string out_filename = path + "/" + mN_ + ".out";
  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadDaltonNMRdata : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }
  
  int iatom = 0;
  Matrix tmpTensor(3,3);
  string line;
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 20 ) {
      string match = line.substr(1,23);
      
      //printf("DEBUG: match = %s \n", match.c_str());
      
      if ( match == " Total shielding tensor") {
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	
	string trash;
	for (int j=0;j<3;j++ ){
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> tmpTensor(j,0);
	  iss >> tmpTensor(j,1);
	  iss >> tmpTensor(j,2);
	  getline(infile,line);
	}

	if ( abs(tmpTensor.Trace()) < 0.00001 ) {
	  printf("ERROR: Monomer::ReadDaltonNMRdata: zero tensor detected for monomer %d, atom %d \n", Monomer_index,iatom++ );
	  exit(1);
	}
	
	GetAtom(iatom).SetMonomer3x3Tensor(tmpTensor);
	GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	iatom++;
	fflush(stdout);

	//tmpTensor.Print("DEBUG: Monomer::ReadDalton");

      }
      
    }
  }
  
}

// JDH
void Monomer::ReadQChemNMRdata() {
  printf("Haven't made Monomer::ReadQChemNMRdata() yet\n");
  exit(1);
}

void Monomer::ReadG09ClusterEFGData() {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename = path + "/";
  char label[20];
  sprintf(label,"cluster.%d.log", GetIndex() );
  filename += label;
  
  ifstream infile;
  infile.open( filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("ReadG09ClusterEFGData(): Cannot open file '%s'\n",
	   filename.c_str() );
    exit(1);
  }
  
  Matrix diagonal(GetNumberOfAtoms(),3);
  Matrix offdiagonal(GetNumberOfAtoms(),3);
  string line;

  // Read in the tensor data
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 30 ) {
      if ( Params::Parameters().UseScaledEFGTensors() ) {
	string match = line.substr(19,52);
	if ( match == "3XX-RR        3YY-RR        3ZZ-RR" ) {
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<GetNumberOfAtoms();i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> diagonal(i,0);
	    iss >> diagonal(i,1);
	    iss >> diagonal(i,2);
	    getline(infile,line);	    
	  }
	}
	
	match = line.substr(21,50);
	if ( match == "XY            XZ            YZ") {
	  getline(infile,line); 
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<GetNumberOfAtoms();i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> offdiagonal(i,0);
	    iss >> offdiagonal(i,1);
	    iss >> offdiagonal(i,2);
	    getline(infile,line);
	  }
	  
	}
	
      } else {
	string match = line.substr(21,50);
	if ( match == "XX            YY            ZZ") {
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<GetNumberOfAtoms();i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> diagonal(i,0);
	    iss >> diagonal(i,1);
	    iss >> diagonal(i,2);
	    getline(infile,line);
	  }
	}
	
	if ( match == "XY            XZ            YZ") {
	  getline(infile,line); 
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<GetNumberOfAtoms();i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> offdiagonal(i,0);
	    iss >> offdiagonal(i,1);
	    iss >> offdiagonal(i,2);
	    getline(infile,line);
	  }
	  
	}
      }
      
    }
  }
  
  infile.close();
  
  for ( int iatom=0;iatom < GetNumberOfAtoms(); iatom++) {
    Matrix tmpTensor(3,3);
    tmpTensor(0,0) = diagonal(iatom,0);
    tmpTensor(1,1) = diagonal(iatom,1);
    tmpTensor(2,2) = diagonal(iatom,2);
    
    tmpTensor(0,1) = offdiagonal(iatom,0);
    tmpTensor(0,2) = offdiagonal(iatom,1);
    tmpTensor(1,2) = offdiagonal(iatom,2); 
    
    tmpTensor(1,0) = tmpTensor(0,1);
    tmpTensor(2,0) = tmpTensor(0,2);
    tmpTensor(2,1) = tmpTensor(1,2);
    
    GetAtom(iatom).SetCluster3x3Tensor(tmpTensor);
    fflush(stdout);
  }

}

void Monomer::ReadG09ClusterNMRData() {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename = path + "/";
  char label[20];
  sprintf(label,"cluster.%d.log", GetIndex() );
  filename += label;
  
  ifstream infile;
  infile.open( filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("ReadG09ClusterNMRData(): Cannot open file '%s'\n",
	   filename.c_str() );
    exit(1);
  }

  // CHECK IF MP2 (should really just parse the $g09 input :/)
  int iatom = 0;
  string line;
  bool MP2_JOB = false;
  string search_match = "SCF GIAO Magnetic shielding tensor";
  while ( !infile.eof() ) {
    getline(infile,line);
    
    if ( line.length() > 30 ) {
      //printf("%s\n", line.substr(2,4).c_str() );
      
      // Search for MP2 NMR section 
      string match = line.substr(1,3);
      //printf("%s\n", match.c_str() );
      
      if ( match  == "MP2" ) {
	//printf("%s\n", match.c_str() );
	MP2_JOB = true;
	search_match = "MP2 GIAO Magnetic shielding tensor";
      }
      
    }
  }

  infile.close(); // To rewind
  infile.open( filename.c_str() );

  Matrix tmpTensor(3,3);
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 30 ) {
      string match = line.substr(1,34);
      
      if ( match  == search_match ) {
	// Loop over all the atoms and read in the tensor:
	for (int iatom=0;iatom<GetNumberOfAtoms();iatom++) {
	  getline(infile,line); // Throw away frist line...
	  getline(infile,line);
	  
	  string trash;
	  for (int j=0;j<3;j++ ){
	    istringstream iss(line);
	    iss >> trash;
	    iss >> tmpTensor(j,0);
	    //printf("tmpTensor(%d,0) = %f\n", j, tmpTensor(j,0) );
	    iss >> trash;
	    iss >> tmpTensor(j,1);
	    iss >> trash;
	    iss >> tmpTensor(j,2);
	    getline(infile,line);
	  }
	  
	  GetAtom(iatom).SetCluster3x3Tensor(tmpTensor);

	}
      }
    }
  }
  
  infile.close();

}


void Monomer::ReadDaltonClusterNMRData() {

  string path;
  path = Params::Parameters().GetQMPath();
  string filename = path + "/";
  char label[20];
  sprintf(label,"cluster.%d.out", GetIndex() );
  filename += label;
  
  ifstream infile;
  infile.open( filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Cluster::NMRCluster : Cannot open file '%s'\n",
	   filename.c_str() );
    exit(1);
  }

  Matrix tmpTensor(3,3);
  string line;
  int iatom=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 20 ) {
      string match = line.substr(1,23);
      
      if ( match == " Total shielding tensor" && iatom < GetNumberOfAtoms()) {
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);
	
	string trash;
	for (int j=0;j<3;j++ ){
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> tmpTensor(j,0);
	  iss >> tmpTensor(j,1);
	  iss >> tmpTensor(j,2);
	  getline(infile,line);
	}
	GetAtom(iatom).SetCluster3x3Tensor(tmpTensor);
	iatom++;
      }
    }
    
    if ( iatom == GetNumberOfAtoms() )
      break;

  }

  infile.close();

}

void Monomer::ReadG09EFGData() {
  string path;
  path = Params::Parameters().GetQMPath();  
  
  string out_filename = path + "/" + mN_ + ".log";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadG09EFGData : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  Matrix diagonal(Natoms,3);
  Matrix offdiagonal(Natoms,3);
  string line;

  // Read in the tensor data
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 30 ) {
      if ( Params::Parameters().UseScaledEFGTensors() ) {
	// Get the Diagaonal Terms:
	string match = line.substr(19,52);
	if ( match == "3XX-RR        3YY-RR        3ZZ-RR" ) {
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<Natoms;i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> diagonal(i,0);
	    iss >> diagonal(i,1);
	    iss >> diagonal(i,2);
	    getline(infile,line);
	  }
	  
	  // getline(infile,line); 
	  // getline(infile,line);
	  // getline(infile,line);
	  // getline(infile,line);
	  // getline(infile,line);
	  
	  // for ( int i=0;i<Natoms;i++) {
	  //   istringstream iss(line);
	  //   iss >> trash;
	  //   iss >> trash;
	  //   iss >> offdiagonal(i,0);
	  //   iss >> offdiagonal(i,1);
	  //   iss >> offdiagonal(i,2);
	  //   getline(infile,line);
	  // }
	}
      } else {
	// Get all terms:
	string match = line.substr(21,50);   
	if ( match == "XX            YY            ZZ") {
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<Natoms;i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> diagonal(i,0);
	    iss >> diagonal(i,1);
	    iss >> diagonal(i,2);
	    getline(infile,line);
	  }
	  
	  getline(infile,line); 
	  getline(infile,line);
	  getline(infile,line);
	  getline(infile,line);
	  getline(infile,line);
	  
	  for ( int i=0;i<Natoms;i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> offdiagonal(i,0);
	    iss >> offdiagonal(i,1);
	    iss >> offdiagonal(i,2);
	    getline(infile,line);
	  }
	}
      }

    } // end if > 30
  } //end while 

  infile.close();

  // Now get the Off-diagonal elements in the event we pulled only the Scaled diagonal components of the EFG tensors
  if ( Params::Parameters().UseScaledEFGTensors() ) {
    infile.open( out_filename.c_str()  );
    while ( !infile.eof() ) {
      getline(infile,line);
      if ( line.length() > 30 ) {
	string match = line.substr(21,50);  
	if ( match == "XX            YY            ZZ") {
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<Natoms;i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> trash;
	    iss >> trash;
	    iss >> trash;
	    getline(infile,line);
	  }
	  
	  getline(infile,line); 
	  getline(infile,line);
	  getline(infile,line);
	  getline(infile,line);
	  getline(infile,line);
	  
	  for ( int i=0;i<Natoms;i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> offdiagonal(i,0);
	    iss >> offdiagonal(i,1);
	    iss >> offdiagonal(i,2);
	    getline(infile,line);
	  }
	}
      }
    } // end while loop
    infile.close();
  }


  // diagonal.Print("HERE is the EFG"); Looks good
  //printf("\n");
  //offdiagonal.Print("Off diagonal");

  // Now we reconstrct the tensor, let's put it in the NMR stuff for now.
  for ( int i=0;i<Natoms;i++) {
    Matrix tmpTensor(3,3);
    tmpTensor(0,0) = diagonal(i,0);
    tmpTensor(1,1) = diagonal(i,1);
    tmpTensor(2,2) = diagonal(i,2);
    
    tmpTensor(0,1) = offdiagonal(i,0);
    tmpTensor(0,2) = offdiagonal(i,1);
    tmpTensor(1,2) = offdiagonal(i,2); 
    
    tmpTensor(1,0) = tmpTensor(0,1);
    tmpTensor(2,0) = tmpTensor(0,2);
    tmpTensor(2,1) = tmpTensor(1,2);
    
    GetAtom(i).SetMonomer3x3Tensor(tmpTensor);
    GetAtom(i).SetTwoBody3x3Tensor(tmpTensor);
    fflush(stdout);
  }
}

void Monomer::ReadG09NMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();  
  
  string out_filename = path + "/" + mN_ + ".log";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadG09NMRdata : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  int iatom = 0;
  string line;
  bool MP2_JOB = false;
  string search_match = "SCF GIAO Magnetic shielding tensor";
  bool normal_term = false;
  bool convergence = true;
  while ( !infile.eof() ) {
      getline(infile,line);
      
       if ( line.length() > 40 ) {
	string match = line.substr(1,33);
	if ( match == "Normal termination of Gaussian 09") {
	  normal_term = true;
	}
	match = line.substr(1,40);
	if (match == ">>>>>>>>>> Convergence criterion not met" ) {
	  convergence = false;
	}
      }


      if ( line.length() > 30 ) {
  	//printf("%s\n", line.substr(2,4).c_str() );
	
  	// Search for MP2 NMR section 
  	string match = line.substr(1,3);
  	//printf("%s\n", match.c_str() );

  	if ( match  == "MP2" ) {
  	  //printf("%s\n", match.c_str() );
  	  MP2_JOB = true;
	  search_match = "MP2 GIAO Magnetic shielding tensor";
  	}

      }
  }

  if ( normal_term == false) {
    printf("Monomer %d failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }

  if (  convergence == false) {
    printf("Monomer %d is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  infile.close();

  // Now actually read in the shielding tensor data...
  infile.open( out_filename.c_str() );
  Matrix tmpTensor(3,3);
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 30 ) {
      //printf("DEBUG: print the sub string: %s\n", line.substr(1,34).c_str() ); fflush(stdout);
	
      // Search for MP2 NMR section 
      string match = line.substr(1,34);
      
      if ( match  == search_match ) {
	//printf("DEBUG: here is the string from the file: %s\n", match.c_str() ); fflush(stdout);
	
	// Loop over all the atoms and read in the tensor:
	int iatom = 0;
	for (int i=1;i<=Natoms;i++) {
	  getline(infile,line); // Throw away frist line...
	  getline(infile,line);
	  //printf("DEBUG: first line of tensor: %s \n", line.c_str() ); fflush(stdout);
	  string trash;
	  for (int j=0;j<3;j++ ){
	    istringstream iss(line);
	    iss >> trash;
	    iss >> tmpTensor(j,0);
	    //printf("tmpTensor(%d,0) = %f\n", j, tmpTensor(j,0) );
	    iss >> trash;
	    iss >> tmpTensor(j,1);
	    iss >> trash;
	    iss >> tmpTensor(j,2);
	    getline(infile,line);
	  }
	  

	  if ( abs(tmpTensor.Trace()) < .00001 ) {
	    printf("ERROR: Monomer::ReadG09NMRdata: zero tensor detected for monomer %d, atom %d \n", Monomer_index,iatom++ );
	    exit(1);
	  }

	  GetAtom(iatom).SetMonomer3x3Tensor(tmpTensor);
	  GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	  iatom++;
	  fflush(stdout);
	  //printf("DEBUG: Monomer::ReadG09NMRdata() Shielding Tensor for Monomer: %s, atom %d \n", mN_.c_str(), iatom + 1 );
	  //tmpTensor.Print("DEBUG: Monomer::ReadG09NMRdata()"); fflush(stdout); 
	}

      }
      
    }

  }
  	
  infile.close();
  

}

double Monomer::ReadG09Energy() {
	double energy = 0.0;

  	// Set up the filename, with the full path.  File is e.g. 'm1.out'
  	string path = Params::Parameters().GetQMPath();  
  	string out_filename = path + "/" + mN_ + ".log";
	string qm_method;
	string Header = Params::Parameters().GetGaussianHeader().c_str();
	transform(Header.begin(), Header.end(), Header.begin(), ::toupper);

	if(Header.find("MP2")!=string::npos) {
		qm_method = "MP2";
	} else if(Header.find("HF")!=string::npos) {
		qm_method = "HF";
	} else {
		printf("Monomer::ReadG09Energy : Unknown QM method '%s'\n",out_filename.c_str());
	}

  	// Open the energy file
  	ifstream infile;
  	infile.open(out_filename.c_str());
  	if(!infile.is_open()) {
    		printf("Monomer::ReadG09Energy : Cannot open file '%s'\n",
	   	out_filename.c_str());
    		exit(1);
  	}

    	// Read in the data 
  	string line;
  	while(!infile.eof()) {
    		getline(infile,line);
    		if(qm_method=="MP2") { // Watit
			if(line.find("EUMP2") != string::npos) {
				istringstream iss(line);
				stringstream ss;
				string tmp;
				for(int i=0;i<5;i++) {
					iss >> tmp;
				}
				iss >> tmp;
				int D = tmp.find("D");
				replace(tmp.begin(),tmp.end(),'D','E');
				ss << tmp;
				ss >> energy;
			}
		} else {
			string match = line.substr(0,9);

    			// Search for final g09 energy for HF
    			if(match==" SCF Done") { 
	      			istringstream iss(line);
      				string tmp;
      				for(int i=0;i<4;i++) {
					iss >> tmp; // throw away text 
				}
				iss >> energy; // read energy  
			}
    		}
  	}
 
  	// Close the force file
  	infile.close();
  	return energy;

}

string Monomer::RunGDMAJob() {
  string path;
  path = Params::Parameters().GetQMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";
  
  string infile = path + "/" + mN_ + ".gdma.in";
  string local_infile = mN_ + ".gdma.in";
  string local_outfile = mN_ + ".gdma.out";
  
  cmd += "gdma < " + local_infile + " > " + local_outfile;
  cmd += "; ";
  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
  
}

string Monomer::RunG09Job(bool MM_job) {
  string path;
  if(MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetMMPath();
  }
  if(!MM_job && !Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetQMPath();
  }
  if (MM_job && Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if(!MM_job && Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string infile = path + "/" + mN_ + ".com";

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile;
  local_infile = mN_ + ".com";

  cmd += "g09 " + local_infile;
  cmd += "; ";
  
  // Third command, swith back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// JDH: Returns a command string for running a GMDA G09 Job
string Monomer::RunG09GDMAJob() {
  string path;
  path = Params::Parameters().GetQMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".gdma.com";
  string local_infile = mN_ + ".gdma.com";
  
  cmd += "g09 " + local_infile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;

}

string Monomer::RunFormChk() {
  string path;
  path = Params::Parameters().GetQMPath();
  
  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".gdma.chk";
  string local_infile = mN_ + ".gdma.chk";
  string local_outfile = mN_ + ".gdma.fchk";

  cmd += "formchk " + local_infile + " " + local_outfile;
  cmd += "; ";
  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

string Monomer::RunG09HirshfeldJob() {
  string path;
  path = Params::Parameters().GetMMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".Hirshfeld.com";
  string local_infile = mN_ + ".Hirshfeld.com";
  
  cmd += "g09 " + local_infile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

string Monomer::RunG09ChelpGJob() {
  string path;
  path = Params::Parameters().GetMMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".chelpG.com";
  string local_infile = mN_ + ".chelpG.com";
  
  cmd += "g09 " + local_infile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// Run a PSI4 job
string Monomer::RunPSI4Job(bool MM_job) {
  string path;
  if(MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetMMPath();
  }
  if(!MM_job && !Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetQMPath();
  }
  if (MM_job && Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if(!MM_job && Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string infile = path + "/" + mN_ + ".in";

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile;
  local_infile = mN_ + ".in";

  cmd += "psi4 " + local_infile;
  cmd += "; ";
  
  // Third command, swith back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  printf("%s\n",cmd.c_str());
  return cmd;
}

//PSI4 read energy //CSG
double Monomer::ReadPSI4Energy() {

  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'm1.out'
  string path;
  path = Params::Parameters().GetQMPath();  

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string out_filename = path + "/" + mN_ + ".out"; // Change to suffix of the molpro job - which I think we will just keep as .out...

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) { 
    printf("Monomer::ReadPSI4Energy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }
  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,18);

    // Search for final Molpro energy
    if ( match=="  \"CURRENT ENERGY\"" ) { 
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
        iss >> energy; // read energy
    }
  }
    
  if ( Params::Parameters().PrintLevel() > 0) printf("PSI4 Obtained QM Monomer Energy = %15.9f\n",energy);
  
  // Close the output file
  infile.close();
  
  
  return energy;

}

//JDH
void Monomer::ReadMolProNMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".molpro.out";

  string cmd;

  
  // First we check to see if it is an MP2 job so we know which tensors to read:
  cmd = "cd " + path;
  cmd += "; ";
  char label[120];

  bool UseLMP2 = true;
  if ( UseLMP2 ) {
    //printf("DEBUG XXX Monomer::Reading molpro data: MP2\n");
    // MP2 Job
    // Run Script to create tmp.txt file
    cmd = "cd " + path;
    cmd += "; ";
    sprintf(label,"awk -v n=2 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";


    sprintf(label," | awk '{print $2, $3, $4}' > file1.out");
    cmd += label;
    cmd += "; ";
    sprintf(label,"awk -v n=3 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";
    sprintf(label," | awk '{print $2, $3, $4}' > file2.out");
    cmd += label;
    cmd += "; ";

    sprintf(label,"awk -v n=4 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";
    sprintf(label," | awk '{print $2, $3, $4}' > file3.out");
    cmd += label;
    cmd += "; ";
  } else {
    //printf("DEBUG XXX Monomer::Reading molpro data: HF\n");
    // HF job
    cmd = "cd " + path;
    cmd += "; ";
    sprintf(label,"awk -v n=2 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";
    sprintf(label," | awk '{print $2, $3, $4}' > file1.out");
    cmd += label;
    cmd += "; ";

    
    
    sprintf(label,"awk -v n=3 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";
    sprintf(label," | awk '{print $2, $3, $4}' > file2.out");
    cmd += label;
    cmd += "; ";

    sprintf(label,"awk -v n=4 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./" + mN_ + ".molpro.out";
    sprintf(label," | awk '{print $2, $3, $4}' > file3.out");
    cmd += label;
    cmd += "; ";
  }
  
   system(cmd.c_str() ); // run the jobs to make file1,2,3 which contain the 1st,2nd,3rd lines of the tensors for each atom
  
   
  ifstream infile1, infile2, infile3;
  string outfile1 = path + "/file1.out";
  string outfile2 = path + "/file2.out";
  string outfile3 = path + "/file3.out";
  
  infile1.open( outfile1.c_str() );
  infile2.open( outfile2.c_str() );
  infile3.open( outfile3.c_str() );
  
  string line1, line2, line3;
  Matrix tmpTensor(3,3);
  for (int iatom=0;iatom<Natoms; iatom++) {
    getline(infile1,line1);
    getline(infile2,line2);
    getline(infile3,line3);
    
    istringstream iss1(line1);
    istringstream iss2(line2);
    istringstream iss3(line3);
    
    iss1 >> tmpTensor(0,0);
    iss1 >> tmpTensor(0,1);
    iss1 >> tmpTensor(0,2);
    
    iss2 >> tmpTensor(1,0);
    iss2 >> tmpTensor(1,1);
    iss2 >> tmpTensor(1,2);
    
    iss3 >> tmpTensor(2,0);
    iss3 >> tmpTensor(2,1);
    iss3 >> tmpTensor(2,2);
    


    GetAtom(iatom).SetMonomer3x3Tensor(tmpTensor);
  }

  infile1.close();
  infile2.close();
  infile3.close();

  // now we cleanup the temp files
  cmd = "cd " + path;
  cmd += "; ";
  cmd += "rm file1.out file2.out file3.out";
  cmd += "; ";
  cmd += "cd " + Params::Parameters().GetBasePath();
  system( cmd.c_str() );

}


int Monomer::GetUnitCellIndex( int i){

  if ( i == 1) {
    return na_;
  } else if ( i == 2 ) {
    return nb_; 
  } else if ( i == 3 ) {
    return nc_;
  } else  {
    printf("ERROR Monomer::GetUnitCellIndex(i) i=1:a, i=2:b and i=3:c \n");
    exit(1);
  }
}

void Monomer::SetUnitCellIndex(int i, int index) {
  if ( i==1 ) {
    na_ = index;
  } else if ( i == 2 ) {
    nb_ = index; 
  } else if ( i == 3 ) {
    nc_ = index;
  } else  {
    printf("ERROR Monomer::SetUnitCellIndex(i, index) i=1:a, i=2:b and i=3:c \n");
    exit(1);
  }

}

void Monomer::ReadHirshfeldCharges(){
  string path;
  path = Params::Parameters().GetMMPath();  
  
  string out_filename = path + "/" + mN_ + ".Hirshfeld.log";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadHirshfeldCharges : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  int iatom = 0;
  string line;
  //string search_match = "Charges from ESP fit, RMS="; // 26 characters
  string search_match = "Hirshfeld charges, spin densities"; // 33 characters
  bool normal_term = false;
  bool convergence = true;
  while ( !infile.eof() ) {
    getline(infile,line);
    
    if ( line.length() > 40 ) {
      string match = line.substr(1,33);
      if ( match == "Normal termination of Gaussian 09") {
	normal_term = true;
      }
      match = line.substr(1,40);
      if (match == ">>>>>>>>>> Convergence criterion not met" ) {
	convergence = false;
      }
    }
  }
  
  if ( normal_term == false) {
    printf("Hirshfeld calculation for Monomer %d failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  if (  convergence == false) {
    printf("Hirshfeld calculation Monomer %d is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  infile.close();
  
  infile.open( out_filename.c_str() );
  double *chelpg_charges = new double[4]; // apparently this need so
  chelpg_charges[1]=0;
  chelpg_charges[2]=0;
  chelpg_charges[3]=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 25 ) {
      // Search for Hirshfeld Charges Section
      string match = line.substr(1,33);
      if ( match  == search_match ) {
	getline(infile,line);// Throw away first line
	getline(infile,line);
	//Extract CHelpG Charges and add them to the array
	string trash;
	fflush(stdout);
	for (int i=0;i<Natoms;i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);
	  //printf("CHelpG charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] );

	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  GetAtom(i).SetMultipoleMoments(Charge);

	  //printf("DEBUG XXX Hirshfeld charge: %f\n", chelpg_charges[0] ); fflush(stdout);

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}
      }
    }
    
    
  } // end while loop

  delete [] chelpg_charges;


}


void Monomer::ReadG09ChelpGCharges() {
  string path;
  path = Params::Parameters().GetMMPath();  
  
  string out_filename = path + "/" + mN_ + ".chelpG.log";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadChelpGCharges : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  int iatom = 0;
  string line;
  string search_match = "Charges from ESP fit, RMS=";
  bool normal_term = false;
  bool convergence = true;
  while ( !infile.eof() ) {
    getline(infile,line);
    
    if ( line.length() > 40 ) {
      string match = line.substr(1,33);
      if ( match == "Normal termination of Gaussian 09") {
	normal_term = true;
      }
      match = line.substr(1,40);
      if (match == ">>>>>>>>>> Convergence criterion not met" ) {
      	convergence = false;
      }

      
    }
  }


  if ( convergence == false ) {
    infile.close();
    infile.open( out_filename.c_str() );
    while ( !infile.eof() ) {
      getline(infile,line);
      
      if ( line.length() > 40 ) {
	string match = line.substr(1,33);
	if ( match == "Normal termination of Gaussian 09") {
	  normal_term = true;
	}
	
	
	match = line.substr(1,29);
	if ( match == "Linear equations converged to" ) {
	  getline(infile,line);
	  match = line.substr(1,9);
	  if ( match == "SCF Done:"  ) {
	    convergence = true;
	  }
	}
	
      }
    }
  }
  


  
  
  if ( normal_term == false) {
    printf("ChelpG calculation for Monomer %d failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  if (  convergence == false) {
    printf("ChelpG calculation Monomer %d is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  infile.close();
  
  infile.open( out_filename.c_str() );
  double *chelpg_charges = new double[4]; // apparently this need so
  chelpg_charges[1]=0;
  chelpg_charges[2]=0;
  chelpg_charges[3]=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 25 ) {
      // Search for CHelpG Charges Section
      string match = line.substr(1,26);
      if ( match  == search_match ) {
	// Throw away first two lines
	getline(infile,line);
	getline(infile,line);
	getline(infile,line);

	//Extract CHelpG Charges and add them to the array
	string trash;
	fflush(stdout);
	for (int i=0;i<Natoms;i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);
	  //printf("CHelpG charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] );

	  if ( abs( chelpg_charges[0] ) < 0.000001 ) {
	    printf("Possible Error: encountered a ChelpG charge of zero!\n");
	    exit(1);
	  }


	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  GetAtom(i).SetMultipoleMoments(Charge);
	  

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}
      }
    }
    
    
  } // end while loop

  delete [] chelpg_charges;


}


void Monomer::ReadOrcaChelpGCharges(){
  printf("Reading in Chelp G chages...\n");

  string path;
  path = Params::Parameters().GetMMPath();  
  
  string out_filename = path + "/" + mN_ + ".chelpG.out";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadOrcaChelpGCharges : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  int iatom = 0;
  string line;
  string search_match = "CHELPG Charges";
  bool normal_term = false;
  bool convergence = true;
  while ( !infile.eof() ) {
    getline(infile,line);
    
    if ( line.length() > 40 ) {
      string match = line.substr(0,61);
      
      if ( match == "                             ****ORCA TERMINATED NORMALLY****") {
	normal_term = true;
      }
      // match = line.substr(1,40);
      // if (match == ">>>>>>>>>> Convergence criterion not met" ) {
      // 	convergence = false;
      // }
    }
  }
  
  if ( normal_term == false) {
    printf("ChelpG calculation for Monomer %d failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }
  
  // if (  convergence == false) {
  //   printf("ChelpG calculation Monomer %d is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", GetIndex() );
  //   exit(1);
  // }
  
  infile.close();
  
  infile.open( out_filename.c_str() );
  double *chelpg_charges = new double[4]; // apparently this need so
  chelpg_charges[1]=0;
  chelpg_charges[2]=0;
  chelpg_charges[3]=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 10 ) {
      

      // Search for CHelpG Charges Section
      string match = line.substr(0,14);
      //printf("%s\n", match.c_str() );
      if ( match  == search_match ) {
	
	
	// Throw away first  lines
	getline(infile,line);
	getline(infile,line);

	//Extract CHelpG Charges and add them to the array
	string trash;
	fflush(stdout);
	for (int i=0;i<Natoms;i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);

	  //printf("ChelpG charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] );

	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  GetAtom(i).SetMultipoleMoments(Charge);
	  

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}
      }
    }
    
    
  } // end while loop

  delete [] chelpg_charges;


}

void Monomer::ReadOrcaNMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();  
  
  string out_filename = path + "/" + mN_ + ".out";

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Monomer::ReadOrcaNMRdata : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  int iatom = 0;
  string line;
  string search_match = "CHEMICAL SHIFTS";
  bool normal_term = false;
  bool convergence = true;
  while ( !infile.eof() ) {
      getline(infile,line);
      

      if ( line.length() > 40 ) {
	string match = line.substr(0,61);
	
	if ( match == "                             ****ORCA TERMINATED NORMALLY****") {
	  normal_term = true;
	}
      }

  }

  if ( normal_term == false) {
    printf("Monomer %d failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndex() );
    exit(1);
  }

  
  infile.close();

  // Now actually read in the shielding tensor data...
  infile.open( out_filename.c_str() );
  Matrix tmpTensor(3,3);
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 10 ) {
      //printf("%s\n", line.substr(0,15).c_str() ); fflush(stdout);
	
      // Search for MP2 NMR section 
      string match = line.substr(0,15);
      
      
      if ( match  == search_match ) {
	//printf("DEBUG: here is the string from the file: %s\n", match.c_str() ); fflush(stdout);
	
	for (int skip = 1; skip <= 5; skip++){
	  getline(infile,line); // Throw away frist 19 lines...
	}

	

	// Loop over all the atoms and read in the tensor:
	int iatom = 0;
	for (int i=1;i<=Natoms;i++) {
	  for (int skip = 1; skip <= 14; skip++){
	    getline(infile,line); // Throw away frist 14 lines...
	  }
	  



	  getline(infile,line);
	  //printf("DEBUG: first line of tensor: %s \n", line.c_str() ); fflush(stdout);
	  string trash;
	  for (int j=0;j<3;j++ ){
	    istringstream iss(line);
	    iss >> tmpTensor(j,0);
	    iss >> tmpTensor(j,1);
	    iss >> tmpTensor(j,2);
	    getline(infile,line);
	  }
	  

	  if ( abs(tmpTensor.Trace()) < .00001 ) {
	    printf("ERROR: Monomer::ReadOrcaNMRdata: zero tensor detected for monomer %d, atom %d \n", Monomer_index,iatom++ );
	    exit(1);
	  }

	  GetAtom(iatom).SetMonomer3x3Tensor(tmpTensor);
	  GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	  iatom++;
	  fflush(stdout);


	  // Thow away lines at the end
	  for (int skip = 1; skip <= 9; skip++){
	    getline(infile,line); // Throw away frist 14 lines...
	  }
	  
	  //printf("DEBUG: Monomer::ReadOrcaNMRdata() Shielding Tensor for Monomer: %s, atom %d \n", mN_.c_str(), iatom + 1 );
	  //tmpTensor.Print("DEBUG: Monomer::ReadG09NMRdata()"); fflush(stdout); 
	}

      }
      
    }

  }
  	
  infile.close();
  

}

// Watit - Hess, IR, Raman Intensities

void Monomer::ReadHess(){

        Matrix HessQM(3*Natoms, 3*Natoms);
        Matrix HessMM(3*Natoms, 3*Natoms);
	string path;
	string out_filename;

	// QM Hessian
        if(Params::Parameters().GetQMType()==3) { // G09
		path = Params::Parameters().GetHessianQMPath();
		out_filename = path + "/" + mN_ + ".log";
        }           

        string line;
        ifstream infile;
        int DOF = 3*Natoms; // Define degrees of freedom

        infile.open(out_filename.c_str());
        if(!infile.is_open()) {
                printf("Monomer::ReadHess QM : Cannot open file '%s'\n",out_filename.c_str());
                exit(1);
        }

	if(Params::Parameters().GetQMType()==3) { //G09 Watit
		string line;
		while (!infile.eof()) {
			getline(infile,line);
			if(line.substr(0,17)==" Hessian entering") {
				int maxk=(3*Natoms-3*Natoms%5)/5;
				if(3*Natoms%5!=0) {
					maxk+=1;
				}
				for(int k=0;k<maxk; k++) {
					getline(infile,line);
					for(int i=5*k;i<3*Natoms;i++) {
						getline(infile,line);
						istringstream iss(line);
						string tmp;
						iss >> tmp;
						for(int l=0;l<5;l++) {
							int j=5*k+l;
							stringstream ss;
							string tmp2;
							iss >> tmp2;
							if(tmp2.length()!=0&&tmp2.length()<14) {
								replace(tmp2.begin(),tmp2.end(),'D','E');
								ss << tmp2;
								ss >> HessQM(i,j);
								HessQM(j,i)=HessQM(i,j);
							}
						}
					}
				}
			}
		}
	}
	infile.close();
	SetHessQM(HessQM);
	//HessQM.PrintHessian("1-body QM Hessian");

	//MM Hessian
        if(Params::Parameters().GetMMType()==1) { // Tinker
                path = Params::Parameters().GetHessianMMPath();
                out_filename = path + "/" + mN_ + ".freq";
        }

	infile.open(out_filename.c_str());
        if(!infile.is_open()) {
                printf("Monomer::ReadHess MM : Cannot open file '%s'\n",out_filename.c_str());
                exit(1);
        }
	
	if(Params::Parameters().GetMMType()==1) { // Tinker
		string line;
		while (!infile.eof()) {
			getline(infile,line);
			if(line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):") {
				for(int i=0;i<3*Natoms;i++) {
					for(int j=0;j<3*Natoms;j++) {
						getline(infile,line);
						istringstream iss(line);
						string tmp1,tmp2;
						iss >> tmp1;
						iss >> tmp2;
						iss >> HessMM(i,j);
					}
				}
				break;
			}
		}
	}
        infile.close();
        SetHessMM(HessMM);
        //HessMM.PrintHessian("1-body MM Hessian");

}


void Monomer::ReadIntensity() {
	// Set path to hessian_files/qm
	string path;
	path = Params::Parameters().GetHessianQMPath();
	string out_filename = path + "/" + mN_;

	if(Params::Parameters().GetQMType()==1) { // QChem
		out_filename += ".out";
	} else if(Params::Parameters().GetQMType()==3) { // G09
		out_filename += ".log";
	}

	// Read PolD from output file
	string line;
	ifstream infile;
	int DOF = 3*Natoms; // Define degrees of freedom

	infile.open(out_filename.c_str());
	if(!infile.is_open()) {
		printf("Monomer::ReadIntensity : Cannot open file '%s'\n",out_filename.c_str());
		exit(1);
	}
	Matrix DipD(DOF,3);
	Matrix PolD(DOF,6);
	int check_DipD=0;
	int check_PolD=0;
	if(Params::Parameters().GetQMType()==1) { // QChem
		while (getline (infile,line)) {
			if(line.find("DipDeriv") != string::npos) {
				check_DipD++;
				getline(infile,line);
				for(int i=0;i<DOF;i++) {
					getline(infile,line);
					istringstream iss(line);
					int trash;
					iss >> trash;
					for (int j=0;j<3;j++) {
						iss >> DipD(i,j);
					}
				}
			}
			if(Params::Parameters().DoRaman()) {
  				if(line.find("PolDeriv") != string::npos) {
					check_PolD++;
					getline(infile,line);
					for(int i=0;i<DOF;i++) {
						getline(infile,line);
						istringstream iss(line);
						int trash;
						iss >> trash;
						for (int j=0;j<6;j++) {
							iss >> PolD(i,j);
						}
   					}
  				}
		 	}
		}
	} else if(Params::Parameters().GetQMType()==3) { // G09
		while(getline(infile,line)) {
			if(line.substr(0,17)==" Hessian entering") {
				while(getline(infile,line)) {
                        		if(line.substr(0,12)==" DipoleDeriv") {
                                		check_DipD++;
                                		for(int i=0;i<DOF;i++) {
							string tmp;
							for (int j=0;j<3;j++) {
								stringstream ss;
								tmp = line.substr(j*15+16,15);
								replace(tmp.begin(),tmp.end(),'D','E');
								ss << tmp;
								ss >> DipD(i,j);
								//DipD(i,j)=-DipD(i,j);
							}
							getline(infile,line);
                                        	}
                        			if(Params::Parameters().DoRaman()) {
								if(line.substr(0,11)==" PolarDeriv") {
                                        			check_PolD++;
        	                                		for(int i=0;i<DOF;i++) {
									string tmp;
									for (int j=0;j<3;j++) {
										stringstream ss;
										tmp = line.substr(j*15+16,15);
										replace(tmp.begin(),tmp.end(),'D','E');
										ss << tmp;
										ss >> PolD(i,j);
										//PolD(i,j)=-PolD(i,j);
									}
									getline(infile,line);
                                               			        for (int j=3;j<6;j++) {
                                                       		        	stringstream ss;
                                                                                tmp = line.substr((j-3)*15+16,15);
                                                                		replace(tmp.begin(),tmp.end(),'D','E');
	                                                                	ss << tmp;
        	                                                        	ss >> PolD(i,j);
										//PolD(i,j)=-PolD(i,j);
									}
									getline(infile,line);
                        	                                }
                                	                }
						}
                                        }
                                }
                        }
                }
	}
        //cout << "check " << check_DipD << endl;
        if(check_DipD < 1) {
               	printf("Monomer::ReadIntensity : Cannot find Dipole Derivative for '%s'\n",out_filename.c_str());
               	exit(1);
	}
        if(Params::Parameters().DoRaman()&&check_PolD < 1) {
                printf("Monomer::ReadIntensity : Cannot find PolDerivative for '%s'\n",out_filename.c_str());
	        exit(1);
        }
	infile.close();
	SetDipD(DipD);
	SetPolD(PolD);

	//printf("m%i\n",GetIndex());		
	//DipD.Print("DipoleD");
}



void Monomer::CreateOrcaChelpGJob() {
  // Create ORCA input file:
  string path;
  path = Params::Parameters().GetMMPath();
  string filename;
  filename = path + "/" + mN_ + ".chelpG.inp";


  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateOrcaChelpGJob : cannot open g09 input file '%s'\n",filename.c_str());
    exit(1);
  }


  //  fprintf(job,"%%chk=%s\n", check_file.c_str() );
  //fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  //fprintf(job,"Monomer %d Gaussian input for ChelpG run\n", Monomer_index);
  //fprintf(job,"\n");

  
  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);
  fprintf(job,"*\n");
  fprintf(job,"\n");

  // If using a mixed basis add in the level 1 basis set portion
  // if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
  //   for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
  //     fprintf(job,"%d 0\n", iatom+1);


  //     if ( Params::Parameters().CustomBasis() > 0) {
	
  // 	if ( GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
  // 	} else if ( GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
  // 	} else if ( GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
  // 	} else if ( GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );  
  // 	} else {
  // 	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
  // 	}
	
  //     } else {
  // 	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
  //     }


  //     fprintf(job,"****\n");
  //   }
  //   fprintf(job,"\n");
  // }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  // COMMENT THIS OUT FOR ORCA
  /*
  for ( int iatom=0;iatom<GetNumberOfAtoms(); iatom++) {
    if ( GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  */

  fclose(job);

  //printf("DEBUG: Monomer %d\n", GetIndex() ); fflush(stdout);

  // a bit hackish, but now we use some system calls to remove the 
  // 'charge' and 'nmr' keywords from the G09 input files
  //string job_path = Params::Parameters().GetMMPath(); 
  //string cmd = "cd " + job_path;
  //cmd += "; ";

  string cmd;
  //if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
  // } else if (custom_vdw_radii == false ) {
  //   cmd = "sed -i s/charge/'pop=ChelpG'/ ";
  // }
  //cmd += filename.c_str();
  //system(cmd.c_str());
  cmd = "sed -i s/NMR/CHELPG/ ";
  cmd += filename.c_str();
  system(cmd.c_str());;

}



void Monomer::CreateDaltonJob( Monomer Monomers[], int NMon) {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".dal"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tMonomer: %d\n", Monomer_index);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Spherical Angstrom Nosymmetry\n",GetNumberOfAtoms() );

  for ( int iatom=0;iatom<GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    } else if ( GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    } else if ( GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if ( GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line Monomer::CreateDaltonJob() and add the atom you want in\n");
      exit(1);
    }
    
    if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    }
    GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");

  // fclose(job);

  // filename = path + "/" + mN_ + ".dal";
  // // Open the input file for writing
  // if (( job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }

  //fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);


  // Finally, make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    filename = path + "/" + mN_ + ".pot";
    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {//prevents inclusion of  QM monomer in pot file
	total_atoms += Monomers[imon].GetNumberOfAtoms();
	for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	  int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	  if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	    high++;
	  }
	}
      }
    }


    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");


    int counter=1;

    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != Monomer_index ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }

    // Print the Mutlipole section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }

    // Rank 1
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if (  imon != Monomer_index )  {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	    counter++;
	  }
	}
      }
      
    } // end rank 1
    

    // Rank 2
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if (  imon != Monomer_index ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank>=2) {
	      Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	    }
	    counter++;
	  }
	}
      }
      
    } // end rank 2


    // RANK 3 Terms: Octupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if (  imon != Monomer_index ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++)  {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >=3) {
	      Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	    }
	    counter++;
	  }
	}
      }
      
    } //end rank 3
    

    if ( Params::Parameters().GetChargeEmbeddingRank() > 3 ) {
      printf("ERROR: multipole embedding for NMR calculations only goes up to rank 3: octupole\n");
      exit(1);
    }

    fclose(job);

  } // END CHARGE EMBEDDING SECTION
}



void Monomer::CreateDaltonJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {


  //printf("CHECK TO MAKE SURE MULTIPOLE ARE PRINTING PROPERLY\n");
  //exit(1);

  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".dal"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tMonomer: %d\n", Monomer_index);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Charge=%d Spherical Angstrom Nosymmetry\n",GetNumberOfAtoms(), GetChargeState() );

  for ( int iatom=0;iatom<GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if ( GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;  
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton monomer job, unknown atom type..\n");
      printf("Don't panic, just grep this line in Monomer::CreateDaltonJob() and add the atom you want\n");
      exit(1);
    }
    
    if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    }
    GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");

  // fclose(job);

  // filename = path + "/" + mN_ + ".dal";
  // // Open the input file for writing
  // if (( job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }
  
  // fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);

  // Finally, make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    filename = path + "/" + mN_ + ".pot";
    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {//prevents inclusion of  QM monomer in pot file
	total_atoms += Monomers[imon].GetNumberOfAtoms();
	for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	  int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	  if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	    high++;
	  }
	}
      }
    }

    for (int imon=1; imon<= NMon_images ; imon++) {
      if (MonomerImages[imon].GetUseInEmbedding() ) {
	total_atoms += MonomerImages[imon].GetNumberOfAtoms();
	for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++){
	  fflush(stdout);
	  int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	  if(rank >= 2) {
	    high++;
	  }
	}
      }
    }
    
    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");
    
    int counter=1;
    
    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != Monomer_index ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    // Print the mult. emb. crds. for the lattice
    for (int imon=1; imon<=NMon_images ; imon++) {
      if ( MonomerImages[imon].GetUseInEmbedding() ){
	//printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	  fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
	  counter++;
	}
      }
    }
    
    // Print the Mutlipole section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    for (int imon=1; imon<=NMon_images; imon++) {
      if ( MonomerImages[imon].GetUseInEmbedding() ){
	for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	  Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	  fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	  counter++;
	}
      }
    }

    // Rank 1
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if (  imon != Monomer_index )  {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	    counter++;
	  }
	}
      }
      
      for (int imon=1; imon<= NMon_images ; imon++) {
	if (MonomerImages[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	    counter++;
	  }
	}
      }
    } // end rank 1

    // Rank 2
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if (  imon != Monomer_index ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank>=2) {
	      Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	    }
	    counter++;
	  }
	}
      }
      
      for (int imon=1; imon<=NMon_images ; imon++) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	    int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank>=2) {
	      Vector moments = Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[4],moments[5],moments[6],moments[7],moments[8],moments[9]);
	    }
	    counter++;
	  }
	}
      }
    } // end rank 2
    
    // RANK 3 Terms: Octupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if (  imon != Monomer_index ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++)  {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >=3) {
	      Vector moments = Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(Monomers[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	    }
	    counter++;
	  }
	}
      }
      
      for (int imon=1; imon<=NMon_images; imon++)  {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >=3)  {
	      Vector moments = Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
	      //Vector moments = Multipole::Multipole().Spherical_to_Cartesian(MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments());
	      fprintf(job, "%d  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f  %4.5f\n",counter,moments[10],moments[11],moments[12],moments[13],moments[14],moments[15],moments[16],moments[17],moments[18],moments[19]);
	    }
	    counter++;
	  }
	}
      }
    } //end rank 3
    

    if ( Params::Parameters().GetChargeEmbeddingRank() > 3 ) {
      printf("ERROR: multipole embedding for NMR calculations only goes up to rank 3: octupole\n");
      exit(1);
    }

    fclose(job);

  }  // END PRINT MULTIPOLE EMBEDDING SECTION

}

void Monomer::CreateDaltonJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ) {


  //printf("CHECK TO MAKE SURE MULTIPOLE ARE PRINTING PROPERLY\n");
  //exit(1);

  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".dal"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tMonomer: %d\n", Monomer_index);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Charge=%d Spherical Angstrom Nosymmetry\n",GetNumberOfAtoms(), GetChargeState() );

  for ( int iatom=0;iatom<GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if ( GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton monomer job, unknown atom type..\n");
      printf("Don't panic, just grep this line in Monomer::CreateDaltonJob() and add the atom you want\n");
      exit(1);
    }
    
    if ( GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    }
    GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");

  // fclose(job);

  // filename = path + "/" + mN_ + ".dal";
  // // Open the input file for writing
  // if (( job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }
  
  // fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);

  // Finally, make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    filename = path + "/" + mN_ + ".pot";
    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateDaltonJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    
    
    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {//prevents inclusion of  QM monomer in pot file
	total_atoms += Monomers[imon].GetNumberOfAtoms();
	for(int iatom=0; iatom < Monomers[imon].GetNumberOfAtoms(); iatom++) {
	  int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	  if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	    high++;
	  }
	}
      }
    }

    for (int imon=1; imon<= NMon_images ; imon++) {
      if (MonomerImages[imon].GetUseInEmbedding() ) {
	total_atoms += MonomerImages[imon].GetNumberOfAtoms();
      }
    }

    total_atoms += EwaldCharges.GetRows(); //Add in Ewald charges
    
    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");
    
    int counter=1;
    
    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != Monomer_index ){ 
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", Monomers[imon].GetAtom(iatom).GetSymbol().c_str(), Monomers[imon].GetAtom(iatom).GetPosition(0),Monomers[imon].GetAtom(iatom).GetPosition(1),Monomers[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    // Print the mult. emb. crds. for the lattice
    for (int imon=1; imon<=NMon_images ; imon++) {
      if ( MonomerImages[imon].GetUseInEmbedding() ){
	//printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	  fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
	  counter++;
	}
      }
    }

    for (int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++) {
      fprintf(job, "X  %4.8f  %4.8f  %4.8f\n", EwaldCharges.Element(ichrg,0), EwaldCharges.Element(ichrg,1), EwaldCharges.Element(ichrg,2) );
    }
    
    // Print the Mutlipole section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != Monomer_index ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    
    for (int imon=1; imon<=NMon_images; imon++) {
      if ( MonomerImages[imon].GetUseInEmbedding() ){
	for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	  Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	  fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	  counter++;
	}
      }
    }

    for (int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++) {
      fprintf(job, "%d  %4.5f\n",counter,EwaldCharges.Element(ichrg,3));
      counter++;
    }
    

    if ( Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      printf("ERROR: multipole embedding for NMR calculations only goes up to rank 0 when using Ewald code\n");
      exit(1);
    }

    fclose(job);

  }  // END PRINT MULTIPOLE EMBEDDING SECTION

}


string Monomer::RunDaltonJob() {
  string path;
  path = Params::Parameters().GetQMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_;
  string local_infile = mN_;
  string local_outfile = mN_ + ".out";
  

  // adapt this for parallel runs: dalton -N $NP -o mN_.out mN_.dal
  cmd += "dalton -mb 2000 -noarch -nobackup -o " + local_outfile + " " + local_infile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();


  return cmd;
}


string Monomer::RunOrcaChelpGJob() {
  string path;
  path = Params::Parameters().GetMMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".chelpG.inp";
  string local_infile = mN_ + ".chelpG.inp";
  string local_outfile = mN_ + ".chelpG.out";
  
  cmd += "orca " + local_infile + " > " + local_outfile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

string Monomer::RunOrcaJob() {
  string path;
  path = Params::Parameters().GetQMPath();

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  string infile = path + "/" + mN_ + ".inp";
  string local_infile = mN_ + ".inp";
  string local_outfile = mN_ + ".out";
  
  cmd += "orca " + local_infile + " > " + local_outfile; 
  cmd += "; ";

  // switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}



void Monomer::CreateOrcaJob( Monomer Monomers[], int NMon ) {

  //printf("ERROR: Monomer::CreateG09Job( Monomer Monomers[], int NMon )  funtion not finished: complete basis section, and make sure ee is treated properly\n");
  //exit(1);
  
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".inp"; 


  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateOrcaJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }


  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    string charge_file;
    charge_file = mN_ + ".chrg";
    fprintf(job,"%%pointcharges \"%s\"\n",charge_file.c_str());
  }

  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);
  


  

  fprintf(job,"*\n");
  fprintf(job,"\n");

 
  fclose(job);


  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    printf("ERROR: ORCA not functional with Electrostatic Embedding\n");
    exit(1);
    path = Params::Parameters().GetQMPath();
    filename = path + "/" + mN_ + ".chrg"; 

    // Open the input file for writing
    FILE *job;
    if (( job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Monomer::CreateOrcaJob : cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    int num_charges = 0;
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index) {
	num_charges += Monomers[i].GetNumberOfAtoms();
      }
    }

    fprintf(job,"%d\n", num_charges);
    // Optionally print $external_charges section
    for (int i=1;i<=NMon;i++) {
      if (i != Monomer_index) {
	Monomers[i].PrintEmbeddingCharges(job);
      }
    }
    fclose(job);
  }



}


void Monomer::CreateOrcaJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".inp"; 

  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateOrcaJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);
  

  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {

    // Optionally print $external_charges section
    printf("ERROR ORCA not compatible with electrostatic embedding\n");
    exit(1);

  }
  
  fprintf(job,"*\n");
  fprintf(job,"\n");
  
  fclose(job);

}



void Monomer::CreateOrcaJob( Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges )  {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  filename = path + "/" + mN_ + ".inp"; 



  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateOrcaJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", GetChargeState(), GetSpinState() );
  PrintMonomerCartesian(job);
  

  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {

    // Optionally print $external_charges section
    printf("ERROR ORCA not compatible with electrostatic embedding\n");
    exit(1);

  }
  
  // // Optionally print $external_charges section
  // for (int i=1;i<=NMon;i++) {
  //   if (i != Monomer_index  && Monomers[Monomer_index].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
  //     Monomers[i].PrintEmbeddingCharges(job);
  //   }
  // }

  // for (int i=1;i<=NMon_images;i++) {
  //   if ( MonomerImages[i].GetUseInEmbedding() && Monomers[Monomer_index].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff() ) {
  //     MonomerImages[i].PrintEmbeddingCharges(job);
  //   }
  // }

  // // Charges on the Ewald Sphere:
  // for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
  //   if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
  //     fprintf(job, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
  //   } else if ( Params::Parameters().GetQMPackage() == "G09" ){
  //     fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
  //   } else if ( Params::Parameters().GetQMPackage() == "ORCA" ){
  //     fprintf(job, "Q  %10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
  //   } 


  // }

  
  fprintf(job,"*\n");
  fprintf(job,"\n");



 
  
  fclose(job);



}
