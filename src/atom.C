#include <string>
using std::string;
#include "atom.h"
#include "constants.h"
using namespace hmbi_constants;

Atom::Atom() : AtSym("XX"), AtNum(-1), connectivity(NULL) {
  // empty constructor

  init_pols = false;
  init_multipoles = false;
  init_InduceMultipoles = false;// by shuhao

  Rotation.Initialize(3,3);//Rotation matrix is the identity by default; by yoni
  Rotation.Set_Iden();

  Fractional_Translation.Initialize(3);
  Shift_Vector.Initialize(3);

  Monomer3x3Tensor.Initialize(3,3); 
  TwoBody3x3Tensor.Initialize(3,3);
  Cluster3x3Tensor.Initialize(3,3);
}


// Copy constructor
Atom::Atom(const Atom &other) {
  //printf("Using Atom copy constructor\n"); fflush(stdout);
  if (this!=&other) {  
    AtSym = other.AtSym;
    AtNum = other.AtNum;
    mass = other.mass;
    DispersionAtomType_ = other.DispersionAtomType_;
    C6 = other.C6;
    C8 = other.C8;
    C10 = other.C10;
    C6_UCHF = other.C6_UCHF;
    C8_UCHF = other.C8_UCHF;
    C10_UCHF = other.C10_UCHF;
    isotropic_dipole_polarizability = other.isotropic_dipole_polarizability;

    xyz.Initialize(3);
    xyz = other.xyz;

    local_xyz.Initialize(3);
    local_xyz = other.local_xyz;
    
    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      fractional_xyz.Initialize(3);//by yoni
      fractional_xyz = other.fractional_xyz;//by yoni
    }
    
    freq_pol_dipole.Initialize(10); // by shuhao
    freq_pol_dipole = other.freq_pol_dipole;// by shuhao
    
    freq_pol_quad.Initialize(10); // by shuhao
    freq_pol_quad = other.freq_pol_quad;// by shuhao

    point_charge = other.point_charge;

    atom_index = other.atom_index;
    global_index = other.global_index; //by yoni

    //Symmetry information
    freeze_atom = other.freeze_atom;//by yoni
    change_x_sign = other.change_x_sign;//by yoni
    change_y_sign = other.change_y_sign;//by yoni
    
    freeze_x = other.freeze_x;//by yoni
    freeze_y = other.freeze_y;//by yoni
    freeze_z = other.freeze_z;//by yoni
    lock_x = other.lock_y;//by yoni
    lock_y = other.lock_y;//by yoni


    lock_z = other.lock_z;//by yoni    
    change_z_sign = other.change_z_sign;//by yoni


    Shift_Vector = other.Shift_Vector;//by yoni
    Rotation = other.Rotation; //by yoni
    Fractional_Translation = other.Fractional_Translation;//by yoni
    Sym_Atom = other.Sym_Atom; //by yoni

    



    MM_atom_type = other.MM_atom_type;
    Nconnected = other.Nconnected;

    

    // connectivity
    if ( Nconnected > 0 ) {
      connectivity = new int[Nconnected];
      for (int i=0;i<Nconnected;i++) {
	connectivity[i] = other.connectivity[i];
      }
    }
    else {
      connectivity = NULL;
    }

    MultipoleMoments = other.MultipoleMoments;
    Pols = other.Pols;

    // NMR stuff...
    Monomer3x3Tensor = other.Monomer3x3Tensor;
    TwoBody3x3Tensor = other.TwoBody3x3Tensor;
    Cluster3x3Tensor = other.Cluster3x3Tensor;
    MixedBasisRegion_ = other.MixedBasisRegion_;
    EwaldPotential = other.EwaldPotential;
    
    EmbeddedCharge = other.EmbeddedCharge;
    EmbeddedChargeTwoBody = other.EmbeddedChargeTwoBody;

  init_pols = other.init_pols;
  init_multipoles = other.init_multipoles;
  init_InduceMultipoles = other.init_InduceMultipoles;

  }
  else {
    printf("Atom:: Using copy constructor to initialize self\n");
  }

}

Atom::~Atom() {
  delete [] connectivity;
}

void Atom::Initialize(int index,int previous_atoms, string symbol, double x, double y, double z, double charge)  {



  atom_index = index;
  global_index = previous_atoms+index;//by yoni

  //symmetry
  Sym_Atom = global_index;//symmetrical to itself unless otherwise set; by yoni
  freeze_atom = 0;//by yoni
  lock_x = 0;//
  lock_y = 0;//by yoni
  change_y_sign = 0;//by yoni

  freeze_x = 0;//by yoni
  freeze_y = 0;//by yoni
  freeze_z = 0;//by yoni


  //get rid of later
  change_z_sign = 0;//by yoni
  lock_z = 0;//by yoni

  // printf("%s has a global index of %i\n",
  //	 symbol.c_str(),global_index); 
  if (symbol.length() > 2) {
    printf("Atom::Atom Error: Atomic symbol too long\n");
    exit(1);
  }

  // Set AtSym
  // Set string and proper capitalization
  AtSym = symbol;
  AtSym[0] = toupper(AtSym[0]);
  AtSym[1] = tolower(AtSym[1]);

  // Set point charge value
  point_charge = charge;

  Nconnected = 0; // not necessary unless using tinker
  //connectivity = new int[2]; // allocate trivial amount of space,
			      // just as a place holder so we can
			      // safely delete it later.
  connectivity = NULL;
 
  // Set coordinates
  xyz.Initialize(3);
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;

  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      // Fractional Coordinates will be set later
      fractional_xyz.Initialize(3);
    }
  local_xyz.Initialize(3); // empty initialization.  Gets set later

  freq_pol_dipole.Initialize(10); // by shuhao
  freq_pol_quad.Initialize(10); // by shuhao
  
  SetAtomicMass();

  C6 = 0.0; C8=0.0; C10=0.0;
  C6_UCHF = 0.0; C8_UCHF = 0.0; C10_UCHF = 0.0; 
  isotropic_dipole_polarizability = 0.0;

  //NMR Stuff...
  EmbeddedCharge = 9999;
  Monomer3x3Tensor.Initialize(3,3);
  TwoBody3x3Tensor.Initialize(3,3);
  Cluster3x3Tensor.Initialize(3,3);
  EwaldPotential = 0;
  EmbeddedCharge = 0;
  EmbeddedChargeTwoBody =0;
  
}

void Atom::Initialize(int index,int previous_atoms, string symbol, double coord[3], int type, 
		      int Nconnect, int connect[], double charge) {

  atom_index = index;
  global_index = previous_atoms+index;//by yoni
  //printf("%s has a global index of %i\n",
  //	 symbol.c_str(),global_index);

  //symmetry
  Sym_Atom = global_index;//symmetrical to itself unless otherwise set; by yoni
  freeze_atom = 0;//by yoni
  lock_x = 0;//by yoni
  lock_y = 0;//by yoni
  change_x_sign = 0;//by yoni
  change_y_sign = 0;//by yoni

  freeze_x = 0;//by yoni
  freeze_y = 0;//by yoni
  freeze_z = 0;//by yoni


  //get rid of later
  change_z_sign = 0;//by yoni
  lock_z = 0;//by yoni

  // Set AtSym
  if (symbol.length() > 2) {
    printf("Atom::Atom Error: Atomic symbol too long: %s\n",symbol.c_str());
    exit(1);
  }

  // Set string and proper capitalization
  AtSym = symbol;
  AtSym[0] = toupper(AtSym[0]);
  AtSym[1] = tolower(AtSym[1]);

  // Set point charge value
  point_charge = charge;

  // Set coordinates
  xyz.Initialize(coord,3);
  local_xyz.Initialize(3); // empty initialization.  Gets set later
  if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
    fractional_xyz.Initialize(3);
  }

  freq_pol_dipole.Initialize(10); // by shuhao
  freq_pol_quad.Initialize(10); // by shuhao

  // Set MM type and connectivity;
  MM_atom_type = type;

  Nconnected = Nconnect;

  connectivity = new int[Nconnected];
  for (int i=0;i<Nconnected;i++) {
    connectivity[i] = connect[i];
  }

  SetAtomicMass();

  C6 = 0.0; C8 = 0.0; C10 = 0.0;
  C6_UCHF = 0.0; C8_UCHF = 0.0; C10_UCHF = 0.0; 
  isotropic_dipole_polarizability = 0.0;

  // JDH Initialize the magnetic properties stuff...
  EmbeddedCharge = 9999;
  Monomer3x3Tensor.Initialize(3,3);
  TwoBody3x3Tensor.Initialize(3,3);
  Cluster3x3Tensor.Initialize(3,3);
  EwaldPotential = 0;
  EmbeddedCharge = 0;
  EmbeddedChargeTwoBody = 0;
  
  //NMRMonomerShieldingTensor.Initialize(3,3);
  //NMRTwoBodyShieldingTensor.Initialize(3,3);
  //NMRClusterShieldingTensor.Initialize(3,3);

}

// Set new coordinates for an atom
void Atom::SetPosition(double new_xyz[3]) {
  for (int i=0;i<3;i++)
    xyz[i] = new_xyz[i];
}

// Set new coordinates for an atom
void Atom::SetPosition(Vector new_xyz) {
  if (new_xyz.GetLength() != 3) {
    printf("Atom::SetPosition() ERROR: Cartesian coordinate vector has length %d\n",new_xyz.GetLength());
    exit(1);
  }

  xyz = new_xyz;
}


// Set local coordinates for an atom  by Ali
void Atom::SetLocalPosition(double new_loc_xyz[3]) {
  for (int i=0;i<3;i++)
    local_xyz[i] = new_loc_xyz[i];
}

//Set fractional coordinates for an atom by Yoni
void Atom::SetFractionalPosition(Vector new_fract){
    if (new_fract.GetLength() != 3) {
    printf("Atom::SetFractionalPosition() ERROR: Cartesian coordinate vector has length %d\n",new_fract.GetLength());
    exit(1);
  }

  fractional_xyz = new_fract;
}

//For turning off symmetry if the crystal belongs to a space group with
//no symmetrically equivalent monomers but the monomers point
//group allows allows for equivalent atoms.
void Atom::ResetSymmetry(){
  Sym_Atom = global_index;
  freeze_atom = 0;
  lock_x = 0;
  lock_y = 0;
  change_x_sign = 0;
  change_y_sign = 0;

  freeze_x = 0;
  freeze_y = 0;
  freeze_z = 0;
  Rotation.Initialize(3,3);
  Rotation.Set_Iden();
  Fractional_Translation.Initialize(3);
  Shift_Vector.Initialize(3);

  //get rid of
  lock_z = 0;
  change_z_sign = 0;

}

void Atom::SetFreq_Pol_Dipole(Vector new_freq_pol_dipole) {
  if (new_freq_pol_dipole.GetLength() != 10) {
    printf("Atom::SetFreq_Pol_Dipole ERROR: SetFreq_Pol_Dipole vector has length %d\n",new_freq_pol_dipole.GetLength());
    exit(1);
  }
  freq_pol_dipole = new_freq_pol_dipole;
  
}

void Atom::SetFreq_Pol_Quad(Vector new_freq_pol_quad) {
  if (new_freq_pol_quad.GetLength() != 10) {
    printf("Atom::SetFreq_Pol_Quad ERROR: SetFreq_Pol_Quad vector has length %d\n",new_freq_pol_quad.GetLength());
    exit(1);
  }
  freq_pol_quad = new_freq_pol_quad;
  
}

int Atom::GetAtomicNumber() {
  // first row is a dummy to let us count from 1.
  string AtomicSymbols[101] = {
    "XX",
    "H","He",
    "Li","Be","B","C","N","O","F","Ne", // period 2
    "Na","Mg","Al","Si","P","S","Cl","Ar", // period 3
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn", // period 4
      "Ga","Ge","As","Se","Br","Kr", // period 4
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
      "In","Sn","Sb","Te","I","Xe", //period 5
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
      "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb",
      "Bi","Po","At","Rn",// period 6
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
      "Es","Fm"};// period 7

  int i, at_num = 0;
  for (i=1;i<101;i++) {
    if ( AtSym==AtomicSymbols[i]) {
      //printf("Match: %s and %s\n",AtSym.c_str(),AtomicSymbols[i].c_str());
      at_num = i;
    }
  }
  return at_num;
}

void Atom::SetAtomicMass() {
  // First mass is a dummy one
  double AtomicMasses[101] = {0.0, // dummy placeholder
    1.00783 ,   4.00260,   7.01600,   9.01218,  11.00931,  // B
    12.00000,  14.00307,  15.99491,  18.99840,  19.99244,  // Ne
    22.9898 ,  23.98504,  26.98153,  27.97693,  30.97376,  // P
    31.97207,  34.96885,  39.948  ,  38.96371,  39.96259,  // Ca
    44.95592,  47.90   ,  50.9440 ,  51.9405 ,  54.9380 ,  // Mn
    55.9349 ,  58.9332 ,  57.9353 ,  62.9296 ,  63.9291 ,  // Zn
    68.9257 ,  73.9219 ,  74.9216 ,  79.9165 ,  78.9183 ,  // Br
    83.80   ,  84.9117 ,  87.9056 ,  88.9059 ,  89.9043 ,  // Zr
    92.9060 ,  97.9055 ,  98.9062 , 101.9037 , 102.9048 ,  // Rh
    105.9032, 106.90509, 113.9036 , 114.9041 , 118.69   ,  // Sn
    120.9038, 129.9067 , 126.9044 , 131.9042 , 132.9051 ,  // Cs
    137.9050, 138.9061 , 139.9053 , 140.9074 , 141.9075 ,  // Nd
    144.913 , 151.9195 , 152.9209 , 157.9241 , 159.9250 ,  // Tb
    163.9288, 164.9303 , 165.9304 , 168.9344 , 173.9390 ,  // Yb
    174.9409, 179.9468 , 180.9480 , 183.9510 , 186.9560 ,  // Re
    192.    , 192.9633 , 194.9648 , 196.9666 , 201.9706 ,  // Hg
    204.9745, 207.9766 , 208.9804 , 208.9825 , 209.987  ,  // At
    222.0175, 223.0198 , 226.0254 , 227.0278 , 232.0382 ,  // Th
    231.0359, 238.0508 , 237.0480 , 244.064  , 243.0614 ,  // Am
    247.070 , 247.0702 , 251.080  , 254.0881 , 257.095  }; // Fm 

  //mass of deuterium
  if(Params::Parameters().IsDeuterated() && GetAtomicNumber() == 1){
    mass = 2.01410178;
  }
  else mass = AtomicMasses[GetAtomicNumber()];
  //printf("%i mass = %f\n",GetAtomicNumber(),mass);

}


//Function returns the translation part of the symmetry operator
Vector Atom::GetFractTranslationVector(bool FromUniqueToThis){


  Vector TransOp;
  
  //The translational operator is part of the symmetry operator going from the symmetry unique monomer to this one.
  if(FromUniqueToThis){
    Matrix InverseRot = Fractional_Rotation;
    InverseRot.Inverse();
    //InverseRot.Print("InverseRot");
    TransOp = InverseRot.MatrixTimesVector(Fractional_Translation);
    TransOp.Scale(-1.0);
    //TransOp.Print("TransOp");
  }else 
    TransOp = Fractional_Translation;
  
  return TransOp;

}


// Compute the distance between this atom and another
double Atom::GetInterAtomicDistance(const Atom& atomB) {
  double dist, dx, dy, dz;
  dx = xyz[0] - atomB.xyz[0];
  dy = xyz[1] - atomB.xyz[1];
  dz = xyz[2] - atomB.xyz[2];
  dist = sqrt( dx*dx + dy*dy + dz*dz );

  return dist;
}
// Computes the Tab matrix for electrostatic interactions. See Stone's
// "Theory of Intermolecular Forces" book, which is where all the
// equations in this subroutine come from.

//this function has been replaced with one that repressents rotation in matrix
//representation instend of axis/angle notation
Matrix Atom::BuildInteractionMatrix(Vector thisRotVec, double thisRotAng, 
				    Atom& other, Vector otherRotVec,
				    double otherRotAng, double beta_damp) {


 double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 

  //Matrix thisTranspose = thisRotation;
  //Matrix otherTranspose = otherRotation;
  //thisTranspose.Transpose();
  //otherTranspose.Transpose();



  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix Tab(nQA,nQB);

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);

  // Find RAB = RA - RB and the distance;
  double Rnorm = GetInterAtomicDistance(other)*AngToBohr;
  // Predefine Rnorm^x here for simplicity in later equations
  double Rnorm2 = Rnorm*Rnorm;
  double Rnorm3 = Rnorm2*Rnorm;
  double Rnorm4 = Rnorm3*Rnorm;
  double Rnorm5 = Rnorm4*Rnorm;

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-5.
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  for (int n=1;n<=5;n++) {
    if (if_damp)
      damp[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm);
    else
      damp[n] = 1.0;
  }

  //if (atom_index==1 && other.atom_index==1) {
  //  printf("Rnorm = %.8f\n",Rnorm);
  //  damp.Print("Tab damping factors\n");
  // }

  
  // Apply the damping factor (which =1 if no damping)
  Rnorm /= damp[1];
  Rnorm2 /= damp[2];
  Rnorm3 /= damp[3];
  Rnorm4 /= damp[4];
  Rnorm5 /= damp[5];

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  Vector eAB(RB);
  eAB -= RA;
  eAB.Normalize();


  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    //eA[i] = thisTranspose.MatrixTimesVector(unit_mat.GetColumnVector(i));
    eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    //eB[i] = otherTranspose.MatrixTimesVector(unit_mat.GetColumnVector(i));
    eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  }


  // Define rA, rB, cAB
  // rA = eA dot eAB... component of eA lying along A->B axis
  // rB = eA dot eBA = eA dot (-eAB) ... component of eB lying along A->B axis
  // cAB(i,j) = rAi dot rBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(eAB);
    rB[i] = -1.0*eB[i].DotProduct(eAB);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1.0/Rnorm;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]/Rnorm2;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]/Rnorm2;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3.0*rA[t-1]*rB[u-1] + cAB(t-1,u-1))/Rnorm3;
	  if (t < dim2) 
	    Tab(u,t) = (3.0*rB[t-1]*rA[u-1] + cAB(u-1,t-1))/Rnorm3;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3.0*raz*raz - 1.0)/Rnorm3;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)/Rnorm3;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz/Rnorm3;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz/Rnorm3;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)/Rnorm3;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby/Rnorm3;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5.0*pow(raz,3.0) - 3.0*raz)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5.0*pow(rbz,3.0) - 3.0*rbz)/Rnorm4;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5.0*raz*raz - 1.0)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5.0*rbz*rbz - 1.0)/Rnorm4;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5.0*raz*raz - 1.0)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5.0*rbz*rbz - 1.0)/Rnorm4;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)/Rnorm4;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz/Rnorm4;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3.0*ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3.0*rby*rby)/Rnorm4;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3.0*rax*rax - ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3.0*rbx*rbx - rby*rby)/Rnorm4;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15.0*raz*raz*rB[u-1] + 6.0*raz*cAB(2,u-1) 
			- 3.0*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15.0*rbz*rbz*rA[u-1] + 6.0*rbz*cAB(u-1,2) 
			  - 3.0*rA[u-1])/Rnorm4;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			      + 5.0*rax*raz*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
				+ 5.0*rbx*rbz*rA[u-1])/Rnorm4;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			      + 5.0*ray*raz*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
				+ 5.0*rby*rbz*rA[u-1])/Rnorm4;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5.0*(rax*rax-ray*ray)*rB[u-1] 
				  + 2.0*rax*cAB(0,u-1) 
				  - 2.0*ray*cAB(1,u-1))/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				    + 2.0*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))/Rnorm4;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			      + ray*cAB(0,u-1))/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
				+ rby*cAB(u-1,0))/Rnorm4;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)/Rnorm5;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)/Rnorm5;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)/Rnorm5;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)/Rnorm5;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby/Rnorm5;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)/Rnorm5;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)/Rnorm5;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))/Rnorm5;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)/Rnorm5;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))/Rnorm5;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))/Rnorm5;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))/Rnorm5;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))/Rnorm5;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))/Rnorm5;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))/Rnorm5;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))/Rnorm5;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)/Rnorm5;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))/Rnorm5;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))/Rnorm5;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))/Rnorm5;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))/Rnorm5;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))/Rnorm5;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))/Rnorm5;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))/Rnorm5;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))/Rnorm5;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))/Rnorm5;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))/Rnorm5;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))/Rnorm5;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))/Rnorm5;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))/Rnorm5;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))/Rnorm5;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  return Tab;

}

// Computes the Tab matrix for electrostatic interactions. See Stone's
// "Theory of Intermolecular Forces" book, which is where all the
// equations in this subroutine come from.
Matrix Atom::BuildInteractionMatrix(Matrix thisRotation,Atom& other, Matrix otherRotation, double beta_damp){

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  

  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 


  //if symmetry is not being exploited by the MM, then elements of the interaction matrix do not have to be rotated
  if(!Params::Parameters().UseMMSymmetry()){
    thisRotation.Set_Iden();
    otherRotation.Set_Iden();
  }

  thisRotation.Transpose();
  otherRotation.Transpose();

  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix Tab(nQA,nQB);

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);

  // Find RAB = RA - RB and the distance;
  double Rnorm = GetInterAtomicDistance(other)*AngToBohr;
  // Predefine Rnorm^x here for simplicity in later equations
  double Rnorm2 = Rnorm*Rnorm;
  double Rnorm3 = Rnorm2*Rnorm;
  double Rnorm4 = Rnorm3*Rnorm;
  double Rnorm5 = Rnorm4*Rnorm;

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-5.
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  for (int n=1;n<=5;n++) {
    if (if_damp)
      damp[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm);
    else
      damp[n] = 1.0;
  }

  //if (atom_index==1 && other.atom_index==1) {
  //  printf("Rnorm = %.8f\n",Rnorm);
  //  damp.Print("Tab damping factors\n");
  // }

  
  // Apply the damping factor (which =1 if no damping)
  Rnorm /= damp[1];
  Rnorm2 /= damp[2];
  Rnorm3 /= damp[3];
  Rnorm4 /= damp[4];
  Rnorm5 /= damp[5];

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  Vector eAB(RB);
  eAB -= RA;
  eAB.Normalize();


  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];
  
  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = thisRotation.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = otherRotation.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  }


  // Define rA, rB, cAB
  // rA = eA dot eAB... component of eA lying along A->B axis
  // rB = eA dot eBA = eA dot (-eAB) ... component of eB lying along A->B axis
  // cAB(i,j) = rAi dot rBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(eAB);
    rB[i] = -1.0*eB[i].DotProduct(eAB);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1.0/Rnorm;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]/Rnorm2;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]/Rnorm2;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3.0*rA[t-1]*rB[u-1] + cAB(t-1,u-1))/Rnorm3;
	  if (t < dim2) 
	    Tab(u,t) = (3.0*rB[t-1]*rA[u-1] + cAB(u-1,t-1))/Rnorm3;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3.0*raz*raz - 1.0)/Rnorm3;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)/Rnorm3;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz/Rnorm3;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz/Rnorm3;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)/Rnorm3;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray/Rnorm3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby/Rnorm3;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5.0*pow(raz,3.0) - 3.0*raz)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5.0*pow(rbz,3.0) - 3.0*rbz)/Rnorm4;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5.0*raz*raz - 1.0)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5.0*rbz*rbz - 1.0)/Rnorm4;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5.0*raz*raz - 1.0)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5.0*rbz*rbz - 1.0)/Rnorm4;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)/Rnorm4;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz/Rnorm4;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3.0*ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3.0*rby*rby)/Rnorm4;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3.0*rax*rax - ray*ray)/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3.0*rbx*rbx - rby*rby)/Rnorm4;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15.0*raz*raz*rB[u-1] + 6.0*raz*cAB(2,u-1) 
			- 3.0*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15.0*rbz*rbz*rA[u-1] + 6.0*rbz*cAB(u-1,2) 
			  - 3.0*rA[u-1])/Rnorm4;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			      + 5.0*rax*raz*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
				+ 5.0*rbx*rbz*rA[u-1])/Rnorm4;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			      + 5.0*ray*raz*rB[u-1])/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
				+ 5.0*rby*rbz*rA[u-1])/Rnorm4;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5.0*(rax*rax-ray*ray)*rB[u-1] 
				  + 2.0*rax*cAB(0,u-1) 
				  - 2.0*ray*cAB(1,u-1))/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				    + 2.0*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))/Rnorm4;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			      + ray*cAB(0,u-1))/Rnorm4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
				+ rby*cAB(u-1,0))/Rnorm4;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)/Rnorm5;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)/Rnorm5;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)/Rnorm5;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)/Rnorm5;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby/Rnorm5;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)/Rnorm5;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)/Rnorm5;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))/Rnorm5;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)/Rnorm5;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))/Rnorm5;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))/Rnorm5;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))/Rnorm5;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))/Rnorm5;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))/Rnorm5;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))/Rnorm5;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))/Rnorm5;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))/Rnorm5;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)/Rnorm5;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))/Rnorm5;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))/Rnorm5;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))/Rnorm5;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))/Rnorm5;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))/Rnorm5;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))/Rnorm5;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))/Rnorm5;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))/Rnorm5;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))/Rnorm5;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))/Rnorm5;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))/Rnorm5;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))/Rnorm5;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))/Rnorm5;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))/Rnorm5;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))/Rnorm5;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))/Rnorm5;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))/Rnorm5;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  return Tab;

}


// shuhao 2010_Aug
// Computes the reciprocal space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"

//yoni Rotations now use matrix representation instead of axes/angle representation
Matrix Atom::BuildRecipInteractionMatrix(Matrix thisRotMat, Atom& other, Matrix otherRotMat,
					 int kx, int ky, int kz, double CellV, Vector RecipCellx, 
					 Vector RecipCelly, Vector RecipCellz, double beta_damp){

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
 
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

 // printf("nQA = %d, nQB = %d\n",nQA,nQB);

  Matrix Tab(nQA,nQB);

  // Grab global position of lattice vectors in reciprocal space, kn, switch to a.u.
  // kn = kxA*+kyB*+kzC*
  Vector kn(3);
  kn[0]=kx*RecipCellx[0]+ky*RecipCelly[0]+kz*RecipCellz[0];
  kn[1]=kx*RecipCellx[1]+ky*RecipCelly[1]+kz*RecipCellz[1];
  kn[2]=kx*RecipCellx[2]+ky*RecipCelly[2]+kz*RecipCellz[2];
  kn.Scale(1/AngToBohr); // the unit of reciprocal vector is length^-1
 // printf("-----kn[0]  = %12.6f,kn[1]  = %12.6f,kn[2]  = %12.6f Bohr\n", kn[0],kn[1],kn[2]);

  //Find the module of kn, |kn|
  double Rkn = sqrt(kn[0]*kn[0]+kn[1]*kn[1]+kn[2]*kn[2]);
 // double inversRkn = 1.0/Rkn;
 // printf("-----Rkn = %12.6f Bohr^-1\n", Rkn);

    // the vectorof  AB
     Vector RA(xyz);
     Vector RB(other.xyz);

     RA.Scale(AngToBohr);
     RB.Scale(AngToBohr);

     Vector rAB(RB);
     rAB -= RA;

 //  printf("-----rAB[0]  = %12.6f,rAB[1]  = %12.6f,rAB[2]  = %12.6f Bohr\n", rAB[0],rAB[1],rAB[2]);
 // calculate the (kn)dot(rAB)
  double kndotrAB = kn[0]*rAB[0]+kn[1]*rAB[1]+kn[2]*rAB[2];
 //  printf("-----kndotrAB  = %12.6f\n", kndotrAB);
  
  //predefine |kn|^x here for simplicity in later equations
  double Rkn2 = Rkn*Rkn;
  double Rkn3 = Rkn2*Rkn;
  double Rkn4 = Rkn3*Rkn;
  double Rkn5 = Rkn4*Rkn;
 // printf("-----Rkn2 = %12.6f Bohr^-2, Rkn3 = %12.6f Bohr^-3,Rkn4 = %12.6f Bohr^-4, Rkn5 = %12.6f Bohr^-5\n", Rkn2,Rkn3,Rkn4,Rkn5);


  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

    // calculate the lattice induction and electrostatic separately using different convergent factors; 
    // for the lattice induction, using a small convergent factor to push the reciprocal energy into direct space for the damping;
    double damping_kappa;
    if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;

    kappa_param /= damping_kappa;

  // calculation of G(kn) which is the total factor of convergence for reciprocal space 
  
  // the unit of Lattice.getvolume() is bohr^3
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  // take the advantage of the symmetry of reciprocal space and the sum can be only over a hemisphere of reciprocal space omitting |kn|=0 (saving the computation cost), 
  // and note the coeffcient need to be changed to VV = 8.0*pi/V;  
  double VV = 4.0*pi/V;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V = %12.6f Bohr^3,VV = %12.6f Bohr^-3,V_V = %12.6f Bohr,VVV = %12.6f Bohr^2\n", V,VV,V_V,VVV);
  
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V;
  double alpha = kappa*kappa;
  double alpha4 = 4.0*alpha;
 // printf("-----alpha = %12.6f, alpha4 = %12.6f\n", alpha, alpha4);
  
  double Gkn = VV*exp(-Rkn2/alpha4);
 // printf("Gkn = %12.6f \n",Gkn); 
  
  
  // calculate the factor related with Ltot (-1)^L/(2L-1)!! * cos^L(kndotrAb) * |kn|^L;
  // note that cos^L(kndotrAb) is the Lth differential cos(kndotrAB).
  // multipling by Gkn is the total factor
  
  // Ltot = 0: (-1)!! = 1, (-1)^0 = 1, 0th differential cos(kndotrAB) = cos(kndotrAB), |kn|^0 = 1. 
  double F0 = cos(kndotrAB)*Gkn/Rkn2;
 // printf("F0 = %12.6f \n",F0);
  
  // Ltot = 1: (1)!! = 1, (-1)^1 = -1, 1th differential cos(kndotrAB) = -sin(kndotrAB), |kn|^1 = Rkn.
  double F1 = sin(kndotrAB)*Gkn/Rkn;
 //  printf("F1 = %12.6f \n",F1);
  
  // Ltot = 2: (3)!! = 3, (-1)^2 = 1, 2th differential cos(kndotrAB) = -cos(kndotrAB), |kn|^2 = Rkn2.
  double F2 = -cos(kndotrAB)*Gkn/3.0;
 //  printf("F2 = %12.6f \n",F2);
  
  // Ltot = 3: (5)!! = 15, (-1)^3 = -1, 3th differential cos(kndotrAB) = sin(kndotrAB), |kn|^3 = Rkn3.
  double F3 = -sin(kndotrAB)*Rkn*Gkn/15.0;
 //  printf("F3 = %12.6f \n",F3);
   
  // Ltot = 4: (7)!! = 105, (-1)^4 = 1, 4th differential cos(kndotrAB) = cos(kndotrAB), |kn|^4 = Rkn4.
  double F4 = cos(kndotrAB)*Rkn2*Gkn/105.0;
 // printf("F4 = %12.6f \n",F4); 

  //Define ekn = kn/norm(kn),the unit vector from global origin -> kn
  Vector ekn(kn);
  ekn.Normalize();   
 //  printf("-----ekn[0]=%12.6f,ekn[1]=%12.6f,ekn[2]=%12.6f Bohr\n", ekn[0],ekn[1],ekn[2]);

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  //Vector eAB(RB);
  //eAB -= RA;
  //eAB.Normalize();

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  thisRotMat.Transpose();
  otherRotMat.Transpose();
  //if symmetry is not being exploited by the MM, then elements of the interaction matrix do not have to be rotated
  if(!Params::Parameters().UseMMSymmetry()){
    thisRotMat.Set_Iden();
    otherRotMat.Set_Iden();
  }

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    //eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i] = thisRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = otherRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  
  //  eA[i].Print("local axis system for atom A");
  //  eB[i].Print("local axis system for atom B");

  } 

  
    // Define rA, rB, cAB
  // knA = rA = eA dot ekn... component of eA lying along kn
  // knB = rB = eB dot (-ekn) ... component of eB lying along kn
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(ekn);
    rB[i] = -1.0*eB[i].DotProduct(ekn);

    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB

  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out 
  // by shuhao, in kn case, no information of 
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }  


  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];


 // rA.Print("rA");
 // rB.Print("rB");
 //  cAB.Print("cAB");

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

 // printf("dim1 = %d, dim2 = %d\n", dim1, dim2);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction 
  // by shuhao:  the factor related with Ltot: [(-1)^Ltot/(2Ltot-1)!! * |Kn|^Ltot * (the Ltotth differential cos(kn dot rAB))]  
  // by shuhao:  the factor not related with Ltot G(kn) added at the end    
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1*F0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*F1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*F1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*F2;
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*F2;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*F2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*F2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*F2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*F2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*F2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*F2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*F3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*F3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*F3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*F3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*F3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*F3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*F3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*F3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*F3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*F3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*F4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)*F4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)*F4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*F4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby*F4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*F4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*F4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*F4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*F4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*F4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*F4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*F4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*F4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*F4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*F4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*F4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*F4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*F4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*F4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*F4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*F4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*F4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*F4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*F4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*F4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*F4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*F4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*F4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*F4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*F4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon // by shuhao multiply by Gkn
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  
 //Tab.Print("TabRecip inside atom.C");
  
  return Tab;

}// shuhao 2010_Aug
// Computes the reciprocal space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"
Matrix Atom::BuildRecipInteractionMatrix(Vector thisRotVec, double thisRotAng, 
				    Atom& other, Vector otherRotVec,
				    double otherRotAng, int kx, int ky, int kz, double CellV,
                                    Vector RecipCellx, Vector RecipCelly, Vector RecipCellz, double beta_damp) {

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
 
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

 // printf("nQA = %d, nQB = %d\n",nQA,nQB);

  Matrix Tab(nQA,nQB);

  // Grab global position of lattice vectors in reciprocal space, kn, switch to a.u.
  // kn = kxA*+kyB*+kzC*
  Vector kn(3);
  kn[0]=kx*RecipCellx[0]+ky*RecipCelly[0]+kz*RecipCellz[0];
  kn[1]=kx*RecipCellx[1]+ky*RecipCelly[1]+kz*RecipCellz[1];
  kn[2]=kx*RecipCellx[2]+ky*RecipCelly[2]+kz*RecipCellz[2];
  kn.Scale(1/AngToBohr); // the unit of reciprocal vector is length^-1
 // printf("-----kn[0]  = %12.6f,kn[1]  = %12.6f,kn[2]  = %12.6f Bohr\n", kn[0],kn[1],kn[2]);

  //Find the module of kn, |kn|
  double Rkn = sqrt(kn[0]*kn[0]+kn[1]*kn[1]+kn[2]*kn[2]);
 // double inversRkn = 1.0/Rkn;
 // printf("-----Rkn = %12.6f Bohr^-1\n", Rkn);

    // the vectorof  AB
     Vector RA(xyz);
     Vector RB(other.xyz);

     RA.Scale(AngToBohr);
     RB.Scale(AngToBohr);

     Vector rAB(RB);
     rAB -= RA;

 //  printf("-----rAB[0]  = %12.6f,rAB[1]  = %12.6f,rAB[2]  = %12.6f Bohr\n", rAB[0],rAB[1],rAB[2]);
 // calculate the (kn)dot(rAB)
  double kndotrAB = kn[0]*rAB[0]+kn[1]*rAB[1]+kn[2]*rAB[2];
 //  printf("-----kndotrAB  = %12.6f\n", kndotrAB);
  
  //predefine |kn|^x here for simplicity in later equations
  double Rkn2 = Rkn*Rkn;
  double Rkn3 = Rkn2*Rkn;
  double Rkn4 = Rkn3*Rkn;
  double Rkn5 = Rkn4*Rkn;
 // printf("-----Rkn2 = %12.6f Bohr^-2, Rkn3 = %12.6f Bohr^-3,Rkn4 = %12.6f Bohr^-4, Rkn5 = %12.6f Bohr^-5\n", Rkn2,Rkn3,Rkn4,Rkn5);


  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

    // calculate the lattice induction and electrostatic separately using different convergent factors; 
    // for the lattice induction, using a small convergent factor to push the reciprocal energy into direct space for the damping;
    double damping_kappa;
    if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;

    kappa_param /= damping_kappa;

  // calculation of G(kn) which is the total factor of convergence for reciprocal space 
  
  // the unit of Lattice.getvolume() is bohr^3
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  // take the advantage of the symmetry of reciprocal space and the sum can be only over a hemisphere of reciprocal space omitting |kn|=0 (saving the computation cost), 
  // and note the coefficient need to be changed to VV = 8.0*pi/V;  
  double VV = 4.0*pi/V;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V = %12.6f Bohr^3,VV = %12.6f Bohr^-3,V_V = %12.6f Bohr,VVV = %12.6f Bohr^2\n", V,VV,V_V,VVV);
  
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V;
  double alpha = kappa*kappa;
  double alpha4 = 4.0*alpha;
 // printf("-----alpha = %12.6f, alpha4 = %12.6f\n", alpha, alpha4);
  
  double Gkn = VV*exp(-Rkn2/alpha4);
 // printf("Gkn = %12.6f \n",Gkn); 
  
  
  // calculate the factor related with Ltot (-1)^L/(2L-1)!! * cos^L(kndotrAb) * |kn|^L;
  // note that cos^L(kndotrAb) is the Lth differential cos(kndotrAB).
  // multipling by Gkn is the total factor
  
  // Ltot = 0: (-1)!! = 1, (-1)^0 = 1, 0th differential cos(kndotrAB) = cos(kndotrAB), |kn|^0 = 1. 
  double F0 = cos(kndotrAB)*Gkn/Rkn2;
 // printf("F0 = %12.6f \n",F0);
  
  // Ltot = 1: (1)!! = 1, (-1)^1 = -1, 1th differential cos(kndotrAB) = -sin(kndotrAB), |kn|^1 = Rkn.
  double F1 = sin(kndotrAB)*Gkn/Rkn;
 //  printf("F1 = %12.6f \n",F1);
  
  // Ltot = 2: (3)!! = 3, (-1)^2 = 1, 2th differential cos(kndotrAB) = -cos(kndotrAB), |kn|^2 = Rkn2.
  double F2 = -cos(kndotrAB)*Gkn/3.0;
 //  printf("F2 = %12.6f \n",F2);
  
  // Ltot = 3: (5)!! = 15, (-1)^3 = -1, 3th differential cos(kndotrAB) = sin(kndotrAB), |kn|^3 = Rkn3.
  double F3 = -sin(kndotrAB)*Rkn*Gkn/15.0;
 //  printf("F3 = %12.6f \n",F3);
   
  // Ltot = 4: (7)!! = 105, (-1)^4 = 1, 4th differential cos(kndotrAB) = cos(kndotrAB), |kn|^4 = Rkn4.
  double F4 = cos(kndotrAB)*Rkn2*Gkn/105.0;
 // printf("F4 = %12.6f \n",F4); 

  //Define ekn = kn/norm(kn),the unit vector from global origin -> kn
  Vector ekn(kn);
  ekn.Normalize();   
 //  printf("-----ekn[0]=%12.6f,ekn[1]=%12.6f,ekn[2]=%12.6f Bohr\n", ekn[0],ekn[1],ekn[2]);

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  //Vector eAB(RB);
  //eAB -= RA;
  //eAB.Normalize();

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  
  //  eA[i].Print("local axis system for atom A");
  //  eB[i].Print("local axis system for atom B");

  } 

  
    // Define rA, rB, cAB
  // knA = rA = eA dot ekn... component of eA lying along kn
  // knB = rB = eB dot (-ekn) ... component of eB lying along kn
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(ekn);
    rB[i] = -1.0*eB[i].DotProduct(ekn);

    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB

  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out 
  // by shuhao, in kn case, no information of 
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }  


  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];


 // rA.Print("rA");
 // rB.Print("rB");
 //  cAB.Print("cAB");

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

 // printf("dim1 = %d, dim2 = %d\n", dim1, dim2);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction 
  // by shuhao:  the factor related with Ltot: [(-1)^Ltot/(2Ltot-1)!! * |Kn|^Ltot * (the Ltotth differential cos(kn dot rAB))]  
  // by shuhao:  the factor not related with Ltot G(kn) added at the end    
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1*F0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*F1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*F1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*F2;
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*F2;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*F2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*F2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*F2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*F2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*F2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*F2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*F3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*F3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*F3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*F3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*F3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*F3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*F3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*F3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*F3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*F3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*F4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)*F4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)*F4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*F4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby*F4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*F4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*F4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*F4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*F4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*F4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*F4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*F4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*F4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*F4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*F4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*F4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*F4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*F4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*F4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*F4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*F4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*F4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*F4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*F4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*F4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*F4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*F4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*F4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*F4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*F4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon // by shuhao multiply by Gkn
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  
 //Tab.Print("TabRecip inside atom.C");
  
  return Tab;

}


// shuhao 2010_Aug
// Computes the reciprocal space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"
// This interaction matrix is used for the case kn = 0 in the reciprocal space

//yoni::updated so the rotation is repressented using matrix notation instend of axis/angle notation
Matrix Atom::BuildRecipInteractionMatrix_kn0(Matrix thisRotMat,Atom& other, Matrix otherRotMat,
					     int kx, int ky, int kz, double CellV,
					     Vector RecipCellx, Vector RecipCelly, Vector RecipCellz, double beta_damp){

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
  
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix Tab(nQA,nQB);

  // Grab global position of lattice vectors in reciprocal space, kn, switch to a.u.
  // kn = kxA*+kyB*+kzC*
  Vector kn(3);
  kn[0]=kx*RecipCellx[0]+ky*RecipCelly[0]+kz*RecipCellz[0];
  kn[1]=kx*RecipCellx[1]+ky*RecipCelly[1]+kz*RecipCellz[1];
  kn[2]=kx*RecipCellx[2]+ky*RecipCelly[2]+kz*RecipCellz[2];
  kn.Scale(1/AngToBohr); // the unit of reciprocal vector is length^-1
  //printf("-----kn[0]  = %12.6f,kn[1]  = %12.6f,kn[2]  = %12.6f Bohr\n", kn[0],kn[1],kn[2]);

  //Find the module of kn, |kn|
  double Rkn = sqrt(kn[0]*kn[0]+kn[1]*kn[1]+kn[2]*kn[2]);
  //printf("-----Rkn = %12.6f Bohr^-1\n", Rkn);
  
  //predefine |kn|^x here for simplicity in later equations
  double Rkn2 = Rkn*Rkn;
  double Rkn3 = Rkn2*Rkn;
  double Rkn4 = Rkn3*Rkn;
  double Rkn5 = Rkn4*Rkn;
  //printf("-----Rkn2 = %12.6f Bohr^-2, Rkn3 = %12.6f Bohr^-3,Rkn4 = %12.6f Bohr^-4, Rkn5 = %12.6f Bohr^-5\n", Rkn2,Rkn3,Rkn4,Rkn5);


  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  double damping_kappa;

  if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;

  kappa_param /= damping_kappa;


  // calculation of G(kn) which is the total factor of convergence for reciprocal space 
  
  // the unit of Lattice.getvolume() is bohr^3
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  //double VV = 8.0*pi/V;
  // for |kn| = 0 the extra factor of 2 introduced in Equation (7) (Mol. Phy. Vol 106, 2008, Page 1567-1589) must be left out.
  // and note for the case of |kn| != 0, the summation  can take advantage of the symmetry of reciprocal space and implies that the sum is over a hemisphere of  reciprocal space omitting |kn| = 0
  // this means for |kn| = 0, VV = 4.0*pi/V; for |kn| != 0, VV = 8.0*pi/V, and sum is from 0 < kn <= kmax not -kmax <= kn <= kmax; 
  double VV = 4.0*pi/V;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V = %12.6f Bohr^3,VV = %12.6f Bohr^-3,V_V = %12.6f Bohr,VVV = %12.6f Bohr^2\n", V,VV,V_V,VVV);
  
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V;
  double alpha = kappa*kappa;
  double alpha4 = 4.0*alpha;
  //printf("-----alpha = %12.6f, alpha4 = %12.6f\n", alpha, alpha4);
  
  //double Gkn = VV*exp(-Rkn2/alpha4)/Rkn2;
  //printf("Gkn = %12.6f \n",Gkn);
  
  
  // the vectorof  AB
     Vector RA(xyz);
     Vector RB(other.xyz);
  
     RA.Scale(AngToBohr);
     RB.Scale(AngToBohr);

     Vector rAB(RB);
     rAB -= RA;  
    
   //  printf("-----rAB[0]  = %12.6f,rAB[1]  = %12.6f,rAB[2]  = %12.6f Bohr\n", rAB[0],rAB[1],rAB[2]);
  
  // calculate the (kn)dot(rAB)
  double kndotrAB = kn[0]*rAB[0]+kn[1]*rAB[1]+kn[2]*rAB[2];
  //  printf("-----kndotrAB  = %12.6f\n", kndotrAB);
  
  // calculate the factor related with Ltot (-1)^L/(2L-1)!! * cos^L(kndotrAb) * |kn|^L;
  // note that cos^L(kndotrAb) is the Lth differential cos(kndotrAB).
  // multiple by Gkn is the total factor
  
  // Ltot = 0: (-1)!! = 1, (-1)^0 = 1, 0th differential cos(kndotrAB) = cos(kndotrAB), |kn|^0 = 1. 
  double F0 = 0.0;
  // printf("F0 = %12.6f \n",F0);
  
  // Ltot = 1: (1)!! = 1, (-1)^1 = -1, 1th differential cos(kndotrAB) = -sin(kndotrAB), |kn|^1 = Rkn.
  double F1 = 0.0;
  // printf("F1 = %12.6f \n",F1);
  
  // Ltot = 2: (3)!! = 3, (-1)^2 = 1, 2th differential cos(kndotrAB) = -cos(kndotrAB), |kn|^2 = Rkn2.
  double F2 = VV/3.0;
  // printf("F2 = %12.6f \n",F2);
  
  // Ltot = 3: (5)!! = 15, (-1)^3 = -1, 3th differential cos(kndotrAB) = sin(kndotrAB), |kn|^3 = Rkn3.
  //double F3 = -sin(kndotrAB)*Rkn*VV*exp(-Rkn2/alpha4)/15.0;
  double F3 = 0.0;
  // printf("F3 = %12.6f \n",F3);
   
  // Ltot = 4: (7)!! = 105, (-1)^4 = 1, 4th differential cos(kndotrAB) = cos(kndotrAB), |kn|^4 = Rkn4.
  //double F4 = cos(kndotrAB)*Rkn2*VV*exp(-Rkn2/alpha4)/105.0;
  double F4 = 0.0;
  // printf("F4 = %12.6f \n",F4);
  

  //Define ekn = kn/norm(kn),the unit vector from global origin -> kn
  Vector ekn(kn);
  //ekn.Normalize(); 
  // because in this case, the kn are (0 0 0) , so no way ekn.Normalize() ; 
  // printf("-----ekn[0]=%12.6f,ekn[1]=%12.6f,ekn[2]=%12.6f Bohr\n", ekn[0],ekn[1],ekn[2]);

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];  

  
  thisRotMat.Transpose();
  otherRotMat.Transpose();
  //if symmetry is not being exploited by the MM, then elements of the interaction matrix do not have to be rotated
  if(!Params::Parameters().UseMMSymmetry()){
    thisRotMat.Set_Iden();
    otherRotMat.Set_Iden();
  }


  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = thisRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = otherRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  } 

  
    // Define rA, rB, cAB
  // knA = rA = eA dot ekn... component of eA lying along kn
  // knB = rB = eB dot (-ekn) ... component of eB lying along kn
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(ekn);
    rB[i] = -1.0*eB[i].DotProduct(ekn);

    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  //cAB.Print("cAB");
 
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out 
  // by shuhao, in kn case, no information of 
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }  


  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction 
  // by shuhao:  the factor related with Ltot: [(-1)^Ltot/(2Ltot-1)!! * |Kn|^Ltot * (the Ltotth differential cos(kn dot rAB))]  
  // by shuhao:  the factor not related with Ltot G(kn) added at the end    
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1*F0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*F1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*F1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*F2;
        //  printf("---Dipole-Dipole interaction Tab terms---\n");
        //  printf("---Tab(%d,%d)=(3*rA[%d-1]*rB[%d-1]+cAB(%d-1,%d-1))*F2=(3*%12.6f*%12.6f+%12.6f)*%12.6f=%12.6f---\n",t,u,t,u,t,u,rA[t-1],rB[u-1],cAB(t-1,u-1),F2,Tab(t,u));
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*F2;
        //  printf("---Tab(%d,%d)=(3*rB[%d-1]*rA[%d-1]+cAB(%d-1,%d-1))*F2=(3*%12.6f*%12.6f+%12.6f)*%12.6f=%12.6f---\n",u,t,t,u,u,t,rB[t-1],rA[u-1],cAB(u-1,t-1),F2,Tab(u,t));
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*F2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*F2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*F2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*F2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*F2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*F2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*F3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*F3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*F3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*F3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*F3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*F3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*F3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*F3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*F3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*F3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*F4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)*F4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)*F4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*F4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby*F4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*F4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*F4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*F4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*F4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*F4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*F4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*F4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*F4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*F4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*F4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*F4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20.0*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*F4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*F4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*F4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20.0*raz*rbx*cAB(2,0) 
				- 20.0*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*F4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*F4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*F4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*F4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*F4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*F4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*F4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*F4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*F4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*F4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*F4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon // by shuhao multiply by Gkn
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  return Tab;

}

// shuhao 2010_Aug
// Computes the reciprocal space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"
// This interaction matrix is used for the case kn = 0 in the reciprocal space
Matrix Atom::BuildRecipInteractionMatrix_kn0(Vector thisRotVec, double thisRotAng, 
				    Atom& other, Vector otherRotVec,
				    double otherRotAng, int kx, int ky, int kz, double CellV,
                                    Vector RecipCellx, Vector RecipCelly, Vector RecipCellz, double beta_damp) {

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
  
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix Tab(nQA,nQB);

  // Grab global position of lattice vectors in reciprocal space, kn, switch to a.u.
  // kn = kxA*+kyB*+kzC*
  Vector kn(3);
  kn[0]=kx*RecipCellx[0]+ky*RecipCelly[0]+kz*RecipCellz[0];
  kn[1]=kx*RecipCellx[1]+ky*RecipCelly[1]+kz*RecipCellz[1];
  kn[2]=kx*RecipCellx[2]+ky*RecipCelly[2]+kz*RecipCellz[2];
  kn.Scale(1/AngToBohr); // the unit of reciprocal vector is length^-1
  //printf("-----kn[0]  = %12.6f,kn[1]  = %12.6f,kn[2]  = %12.6f Bohr\n", kn[0],kn[1],kn[2]);

  //Find the module of kn, |kn|
  double Rkn = sqrt(kn[0]*kn[0]+kn[1]*kn[1]+kn[2]*kn[2]);
  //printf("-----Rkn = %12.6f Bohr^-1\n", Rkn);
  
  //predefine |kn|^x here for simplicity in later equations
  double Rkn2 = Rkn*Rkn;
  double Rkn3 = Rkn2*Rkn;
  double Rkn4 = Rkn3*Rkn;
  double Rkn5 = Rkn4*Rkn;
  //printf("-----Rkn2 = %12.6f Bohr^-2, Rkn3 = %12.6f Bohr^-3,Rkn4 = %12.6f Bohr^-4, Rkn5 = %12.6f Bohr^-5\n", Rkn2,Rkn3,Rkn4,Rkn5);


  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  double damping_kappa;

  if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;

  kappa_param /= damping_kappa;


  // calculation of G(kn) which is the total factor of convergence for reciprocal space 
  
  // the unit of Lattice.getvolume() is bohr^3
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  //double VV = 8.0*pi/V;
  // for |kn| = 0 the extra factor of 2 introduced in Equation (7) (Mol. Phy. Vol 106, 2008, Page 1567-1589) must be left out.
  // and note for the case of |kn| != 0, the summation  can take advantage of the symmetry of reciprocal space and implies that the sum is over a hemisphere of  reciprocal space omitting |kn| = 0
  // this means for |kn| = 0, VV = 4.0*pi/V; for |kn| != 0, VV = 8.0*pi/V, and sum is from 0 < kn <= kmax not -kmax <= kn <= kmax; 
  double VV = 4.0*pi/V;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V = %12.6f Bohr^3,VV = %12.6f Bohr^-3,V_V = %12.6f Bohr,VVV = %12.6f Bohr^2\n", V,VV,V_V,VVV);
  
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V;
  double alpha = kappa*kappa;
  double alpha4 = 4.0*alpha;
  //printf("-----alpha = %12.6f, alpha4 = %12.6f\n", alpha, alpha4);
  
  //double Gkn = VV*exp(-Rkn2/alpha4)/Rkn2;
  //printf("Gkn = %12.6f \n",Gkn);
  
  
  // the vectorof  AB
     Vector RA(xyz);
     Vector RB(other.xyz);
  
     RA.Scale(AngToBohr);
     RB.Scale(AngToBohr);

     Vector rAB(RB);
     rAB -= RA;  
    
   //  printf("-----rAB[0]  = %12.6f,rAB[1]  = %12.6f,rAB[2]  = %12.6f Bohr\n", rAB[0],rAB[1],rAB[2]);
  
  // calculate the (kn)dot(rAB)
  double kndotrAB = kn[0]*rAB[0]+kn[1]*rAB[1]+kn[2]*rAB[2];
  //  printf("-----kndotrAB  = %12.6f\n", kndotrAB);
  
  // calculate the factor related with Ltot (-1)^L/(2L-1)!! * cos^L(kndotrAb) * |kn|^L;
  // note that cos^L(kndotrAb) is the Lth differential cos(kndotrAB).
  // multiple by Gkn is the total factor
  
  // Ltot = 0: (-1)!! = 1, (-1)^0 = 1, 0th differential cos(kndotrAB) = cos(kndotrAB), |kn|^0 = 1. 
  double F0 = 0.0;
  // printf("F0 = %12.6f \n",F0);
  
  // Ltot = 1: (1)!! = 1, (-1)^1 = -1, 1th differential cos(kndotrAB) = -sin(kndotrAB), |kn|^1 = Rkn.
  double F1 = 0.0;
  // printf("F1 = %12.6f \n",F1);
  
  // Ltot = 2: (3)!! = 3, (-1)^2 = 1, 2th differential cos(kndotrAB) = -cos(kndotrAB), |kn|^2 = Rkn2.
  double F2 = VV/3.0;
  // printf("F2 = %12.6f \n",F2);
  
  // Ltot = 3: (5)!! = 15, (-1)^3 = -1, 3th differential cos(kndotrAB) = sin(kndotrAB), |kn|^3 = Rkn3.
  //double F3 = -sin(kndotrAB)*Rkn*VV*exp(-Rkn2/alpha4)/15.0;
  double F3 = 0.0;
  // printf("F3 = %12.6f \n",F3);
   
  // Ltot = 4: (7)!! = 105, (-1)^4 = 1, 4th differential cos(kndotrAB) = cos(kndotrAB), |kn|^4 = Rkn4.
  //double F4 = cos(kndotrAB)*Rkn2*VV*exp(-Rkn2/alpha4)/105.0;
  double F4 = 0.0;
  // printf("F4 = %12.6f \n",F4);
  

  //Define ekn = kn/norm(kn),the unit vector from global origin -> kn
  Vector ekn(kn);
  //ekn.Normalize(); 
  // because in this case, the kn are (0 0 0) , so no way ekn.Normalize() ; 
  // printf("-----ekn[0]=%12.6f,ekn[1]=%12.6f,ekn[2]=%12.6f Bohr\n", ekn[0],ekn[1],ekn[2]);

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  } 

  
    // Define rA, rB, cAB
  // knA = rA = eA dot ekn... component of eA lying along kn
  // knB = rB = eB dot (-ekn) ... component of eB lying along kn
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(ekn);
    rB[i] = -1.0*eB[i].DotProduct(ekn);

    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  //cAB.Print("cAB");
 
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out 
  // by shuhao, in kn case, no information of 
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }  


  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */

  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction 
  // by shuhao:  the factor related with Ltot: [(-1)^Ltot/(2Ltot-1)!! * |Kn|^Ltot * (the Ltotth differential cos(kn dot rAB))]  
  // by shuhao:  the factor not related with Ltot G(kn) added at the end    
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1*F0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*F1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*F1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*F2;
        //  printf("---Dipole-Dipole interaction Tab terms---\n");
        //  printf("---Tab(%d,%d)=(3*rA[%d-1]*rB[%d-1]+cAB(%d-1,%d-1))*F2=(3*%12.6f*%12.6f+%12.6f)*%12.6f=%12.6f---\n",t,u,t,u,t,u,rA[t-1],rB[u-1],cAB(t-1,u-1),F2,Tab(t,u));
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*F2;
        //  printf("---Tab(%d,%d)=(3*rB[%d-1]*rA[%d-1]+cAB(%d-1,%d-1))*F2=(3*%12.6f*%12.6f+%12.6f)*%12.6f=%12.6f---\n",u,t,t,u,u,t,rB[t-1],rA[u-1],cAB(u-1,t-1),F2,Tab(u,t));
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*F2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*F2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*F2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*F2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*F2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*F2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*F2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*F3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*F3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*F3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*F3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*F3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*F3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*F3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*F3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*F3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*F3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*F3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*F3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*F4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*rax*pow(raz,3.0) - 3*rax*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rbx*pow(rbz,3.0) - 3*rbx*rbz)*F4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7*ray*pow(raz,3.0) - 3*ray*raz)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*rby*pow(rbz,3.0) - 3*rby*rbz)*F4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7*raz*raz - 1.0)*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*F4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7*raz*raz - 1.0)*rax*ray*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7*rbz*rbz - 1.0)*rbx*rby*F4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*F4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*F4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*F4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*F4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*F4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*F4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*F4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*F4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*F4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*F4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*F4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*F4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20.0*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*F4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*F4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*F4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20.0*raz*rbx*cAB(2,0) 
				- 20.0*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*F4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*F4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*F4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*F4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*F4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*F4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*F4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*F4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*F4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*F4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*F4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*F4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*F4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*F4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*F4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*F4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon // by shuhao multiply by Gkn
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
  return Tab;

}

// shuhao 2010_Aug
// Computes the direct space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"


//updated by Yoni to use Rotation Matrices instead of Rotation Vectors and Angles

Matrix Atom::BuildDirecInteractionMatrix(Matrix thisRotation, Atom& other, Matrix otherRotation, int nx, int ny, int nz, 
                                         double CellV, Vector UnitCellx, Vector UnitCelly, Vector UnitCellz, double beta_damp) {

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole; 
  //printf("--perm = %12.6f\n", perm);
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
  
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  // printf("nQA = %d, nQB = %d\n",nQA,nQB);

  Matrix Tab(nQA,nQB);

  
  // Grab global position of lattice vectors in direct space, rn, switch to a.u.
  // rn = nxA+nyB+nzC
  
  Vector rn(3);

  rn[0]=nx*UnitCellx[0]+ny*UnitCelly[0]+nz*UnitCellz[0];
  rn[1]=nx*UnitCellx[1]+ny*UnitCelly[1]+nz*UnitCellz[1];
  rn[2]=nx*UnitCellx[2]+ny*UnitCelly[2]+nz*UnitCellz[2];

  rn.Scale(AngToBohr);
  // rn.Print("rn");

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);
  
  Vector rAB(RB);
  rAB -= RA;
  //rAB.Print("rAB");

  // rAB - rn
  Vector rABn(rAB);
  rABn -= rn;
  //rABn.Print("rABn");

  // Find RAB = RA - RB and the distance;
  //double RABn = sqrt(rABn.DotProduct(rABn));
   double RABn =sqrt(rABn[0]*rABn[0]+rABn[1]*rABn[1]+rABn[2]*rABn[2]);
   //printf("-----RABn  = %12.6f Bohr\n", RABn);
  // Predefine Rnorm^x here for simplicity in later equations
  double RABn2 = RABn*RABn;
  double RABn3 = RABn2*RABn;
  double RABn4 = RABn3*RABn;
  double RABn5 = RABn4*RABn;
  double RABn6 = RABn5*RABn;
  //printf("-----RABn = %12.6f Bohr, RABn2 = %12.6f Bohr^2, RABn3 = %12.6f Bohr^3,RABn4 = %12.6f Bohr^4, RABn5 = %12.6f Bohr^5\n", RABn, RABn2,RABn3,RABn4,RABn5);

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-5.
  double damping_kappa;
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  for (int n=1;n<=5;n++) {
    if (if_damp)
      damp[n] = TangToenniesDampingFactor(n, beta_damp, RABn);
    else
      damp[n] = 1.0;
  } 

  if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;
 

  //if (atom_index==1 && other.atom_index==1) {
  //  printf("Rnorm = %.8f\n",Rnorm);
  //  damp.Print("Tab damping factors\n");
  // }


  // Apply the damping factor (which =1 if no damping)
  double DRABn = RABn/damp[1];
  double DRABn2 = RABn2/damp[2];
  double DRABn3 = RABn3/damp[3];
  double DRABn4 = RABn4/damp[4];
  double DRABn5 = RABn5/damp[5];
  kappa_param /= damping_kappa;
  
  //printf("-----kappa_param = %12.6f, DRABn = %12.6f Bohr, DRABn2 = %12.6f Bohr^2, DRABn3 = %12.6f Bohr^3, DRABn4 = %12.6f Bohr^4, DRABn5 = %12.6f Bohr^5\n", kappa_param,DRABn,DRABn2,DRABn3,DRABn4,DRABn5);


  // Define erABn = rABn/norm, the unit vector from AB -> rn
  Vector erABn(rABn);
  erABn.Normalize();
  //erABn.Print("erABn");

  // It is convenient to redefine the error function in terms of the function FL[p]
  double V1 = CellV;
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V1 = %12.6f Angstrom^3,V = %12.6f Bohr^3,VVV = %12.6f Bohr^2\n", V1,V,VVV);
 
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V; 
  double alpha = kappa*kappa;
  double p = alpha*RABn2;
  double pp = kappa*RABn;
  double cc1 = 2.0*kappa/sqrt(pi);
  double cc2 = (4.0/3.0)*pow(kappa,3.0)/sqrt(pi);
  double cc3 = (8.0/15.0)*pow(kappa,5.0)/sqrt(pi);
  double cc4 = (16.0/105.0)*pow(kappa,7.0)/sqrt(pi);
  //printf("-----alpha = %12.6f, p = %12.6f, pp = %12.6f\n", alpha, p, pp);
   const double cc = 105.0/32.0;
  //printf("-----erf(pp) = %12.6f, pow(pp,3.0) = %12.6f, pow(pp,5.0) = %12.6f, pow(p,3.0)= %12.6f, pow(p,4.0)= %12.6f, exp(-p) = %12.6f\n", erf(pp), pow(pp,3.0), pow(pp,5.0),pow(p,3.0),pow(p,4.0),exp(-p));  


  double f0 = 1.0/DRABn - erf(pp)/DRABn;
  double f1 = 1.0/DRABn2 - erf(pp)/DRABn2 + cc1*exp(-p)/DRABn;
  double f2 = 1.0/DRABn3 - erf(pp)/DRABn3 + cc1*exp(-p)/DRABn2 + cc2*exp(-p);
  double f3 = 1.0/DRABn4 - erf(pp)/DRABn4 + cc1*exp(-p)/DRABn3 + cc2*exp(-p)/DRABn + cc3*exp(-p)*RABn;
  double f4 = 1.0/DRABn5 - erf(pp)/DRABn5 + cc1*exp(-p)/DRABn4 + cc2*exp(-p)/DRABn2 + cc3*exp(-p) + cc4*exp(-p)*RABn2;

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  //if symmetry is not being exploited by the MM, then elements of the interaction matrix do not have to be rotated
  if(!Params::Parameters().UseMMSymmetry()){
    thisRotation.Set_Iden();
    otherRotation.Set_Iden();
  }

  thisRotation.Transpose();
  otherRotation.Transpose();

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    //eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i] = thisRotation.MatrixTimesVector(unit_mat.GetColumnVector(i));
    eA[i].Normalize();

    eB[i].Initialize(3);
    // eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);    
    eB[i] = otherRotation.MatrixTimesVector(unit_mat.GetColumnVector(i));
    eB[i].Normalize();

    //eA[i].Print("eA");
    //eB[i].Print("eB");
  }


  // Define rA, rB, cAB
  // rA = eA dot erABn... component of eA lying along rAB-rn axis
  // rB = eA dot (-erABn) ... component of eB lying along rn-rAB axis
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(erABn);
    rB[i] = -1.0*eB[i].DotProduct(erABn);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */
  


  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  //printf("dim1 = %d, dim2 = %d",dim1,dim2);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1.0*f0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*f1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*f1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*f2;
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*f2;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*f2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*f2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*f2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*f2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*f2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*f2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*f3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*f3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*f3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*f3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*f3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*f3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*f3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*f3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*f3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*f3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*f3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*f3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*f3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*f4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*f4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*rax*pow(raz,3.0) - 3*rax*raz)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*rbx*pow(rbz,3.0) - 3*rbx*rbz)*f4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*ray*pow(raz,3.0) - 3*ray*raz)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*rby*pow(rbz,3.0) - 3*rby*rbz)*f4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7.0*raz*raz - 1.0)*(rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7.0*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*f4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7.0*raz*raz - 1.0)*rax*ray*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7.0*rbz*rbz - 1.0)*rbx*rby*f4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*f4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*f4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*f4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*f4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*f4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*f4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*f4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7.0*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7.0*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*f4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7.0*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7.0*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*f4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*f4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7.0*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7.0*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*f4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*f4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*f4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*f4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*f4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*f4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*f4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*f4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*f4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*f4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*f4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*f4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*f4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*f4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*f4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*f4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*f4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*f4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*f4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*f4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*f4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*f4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*f4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }
  */
   
  //Tab.Print("TabDirect inside atom.C");

  return Tab;

}


// shuhao 2010_Aug
// Computes the direct space Tab matrix for electrostatic interactions for periodic cell. See Stone's
// "Theory of Intermolecular Forces" book and M.Leslie's "Molecular Physics Vol. 106, 1567-1578, 2008"
Matrix Atom::BuildDirecInteractionMatrix(Vector thisRotVec, double thisRotAng, 
				    Atom& other, Vector otherRotVec,
				    double otherRotAng, int nx, int ny, int nz, double CellV, 
                                    Vector UnitCellx, Vector UnitCelly, Vector UnitCellz, double beta_damp) {

  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm=4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole; 
  //printf("--perm = %12.6f\n", perm);
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 

  // get the convergence factor
  double kappa_param = Params::Parameters().GetEwaldKappa();
  
  // Allocate storage for Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  // printf("nQA = %d, nQB = %d\n",nQA,nQB);

  Matrix Tab(nQA,nQB);

  
  // Grab global position of lattice vectors in direct space, rn, switch to a.u.
  // rn = nxA+nyB+nzC
  
  Vector rn(3);

  rn[0]=nx*UnitCellx[0]+ny*UnitCelly[0]+nz*UnitCellz[0];
  rn[1]=nx*UnitCellx[1]+ny*UnitCelly[1]+nz*UnitCellz[1];
  rn[2]=nx*UnitCellx[2]+ny*UnitCelly[2]+nz*UnitCellz[2];

  rn.Scale(AngToBohr);
  // rn.Print("rn");

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);
  
  Vector rAB(RB);
  rAB -= RA;
  //rAB.Print("rAB");

  // rAB - rn
  Vector rABn(rAB);
  rABn -= rn;
  //rABn.Print("rABn");

  // Find RAB = RA - RB and the distance;
  //double RABn = sqrt(rABn.DotProduct(rABn));
   double RABn =sqrt(rABn[0]*rABn[0]+rABn[1]*rABn[1]+rABn[2]*rABn[2]);
   //printf("-----RABn  = %12.6f Bohr\n", RABn);
  // Predefine Rnorm^x here for simplicity in later equations
  double RABn2 = RABn*RABn;
  double RABn3 = RABn2*RABn;
  double RABn4 = RABn3*RABn;
  double RABn5 = RABn4*RABn;
  double RABn6 = RABn5*RABn;
  //printf("-----RABn = %12.6f Bohr, RABn2 = %12.6f Bohr^2, RABn3 = %12.6f Bohr^3,RABn4 = %12.6f Bohr^4, RABn5 = %12.6f Bohr^5\n", RABn, RABn2,RABn3,RABn4,RABn5);

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-5.
  double damping_kappa;
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;

  for (int n=1;n<=5;n++) {
    if (if_damp)
      damp[n] = TangToenniesDampingFactor(n, beta_damp, RABn);
    else
      damp[n] = 1.0;
  } 

  if (if_damp)
     damping_kappa = 200.0;
    else
     damping_kappa = 1.0;
 

  //if (atom_index==1 && other.atom_index==1) {
  //  printf("Rnorm = %.8f\n",Rnorm);
  //  damp.Print("Tab damping factors\n");
  // }


  // Apply the damping factor (which =1 if no damping)
  double DRABn = RABn/damp[1];
  double DRABn2 = RABn2/damp[2];
  double DRABn3 = RABn3/damp[3];
  double DRABn4 = RABn4/damp[4];
  double DRABn5 = RABn5/damp[5];
  kappa_param /= damping_kappa;
  
  //printf("-----kappa_param = %12.6f, DRABn = %12.6f Bohr, DRABn2 = %12.6f Bohr^2, DRABn3 = %12.6f Bohr^3, DRABn4 = %12.6f Bohr^4, DRABn5 = %12.6f Bohr^5\n", kappa_param,DRABn,DRABn2,DRABn3,DRABn4,DRABn5);


  // Define erABn = rABn/norm, the unit vector from AB -> rn
  Vector erABn(rABn);
  erABn.Normalize();
  //erABn.Print("erABn");

  // It is convenient to redefine the error function in terms of the function FL[p]
  double V1 = CellV;
  double V = CellV*AngToBohr*AngToBohr*AngToBohr;
  const double cccc = 1.0/3.0;
  double V_V = pow(V,cccc);
  //printf("-----V1 = %12.6f Angstrom^3,V = %12.6f Bohr^3,VVV = %12.6f Bohr^2\n", V1,V,VVV);
 
  // alpha is a positive constant which determines the convergence of the direct and reciprocal space sums.
  // The term -pi/(alpha*V) cancels out for electroneutral cells if the same screening parameter alpha is used for all terms.
  double kappa = kappa_param/V_V; 
  double alpha = kappa*kappa;
  double p = alpha*RABn2;
  double pp = kappa*RABn;
  double cc1 = 2.0*kappa/sqrt(pi);
  double cc2 = (4.0/3.0)*pow(kappa,3.0)/sqrt(pi);
  double cc3 = (8.0/15.0)*pow(kappa,5.0)/sqrt(pi);
  double cc4 = (16.0/105.0)*pow(kappa,7.0)/sqrt(pi);
  //printf("-----alpha = %12.6f, p = %12.6f, pp = %12.6f\n", alpha, p, pp);
   const double cc = 105.0/32.0;
  //printf("-----erf(pp) = %12.6f, pow(pp,3.0) = %12.6f, pow(pp,5.0) = %12.6f, pow(p,3.0)= %12.6f, pow(p,4.0)= %12.6f, exp(-p) = %12.6f\n", erf(pp), pow(pp,3.0), pow(pp,5.0),pow(p,3.0),pow(p,4.0),exp(-p));  


  double f0 = 1.0/DRABn - erf(pp)/DRABn;
  double f1 = 1.0/DRABn2 - erf(pp)/DRABn2 + cc1*exp(-p)/DRABn;
  double f2 = 1.0/DRABn3 - erf(pp)/DRABn3 + cc1*exp(-p)/DRABn2 + cc2*exp(-p);
  double f3 = 1.0/DRABn4 - erf(pp)/DRABn4 + cc1*exp(-p)/DRABn3 + cc2*exp(-p)/DRABn + cc3*exp(-p)*RABn;
  double f4 = 1.0/DRABn5 - erf(pp)/DRABn5 + cc1*exp(-p)/DRABn4 + cc2*exp(-p)/DRABn2 + cc3*exp(-p) + cc4*exp(-p)*RABn2;

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();

    //eA[i].Print("eA");
    //eB[i].Print("eB");
  }


  // Define rA, rB, cAB
  // rA = eA dot erABn... component of eA lying along rAB-rn axis
  // rB = eA dot (-erABn) ... component of eB lying along rn-rAB axis
  // cAB(i,j) = eAi dot eBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(erABn);
    rB[i] = -1.0*eB[i].DotProduct(erABn);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB
  
  // if nQA > nQB, need to swap/transpose arrays to make indexing work
  // out
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  /*
  if (atom_index==1 && other.atom_index==1) {
    printf("Rnorm2 = %f\n",Rnorm2);
    rA.Print("rA");
    rB.Print("rB");
  }
  */
  


  // make sure larger dimension runs first in loop
  int dim1 = max(nQA,nQB);
  int dim2 = min(nQA,nQB);

  //printf("dim1 = %d, dim2 = %d",dim1,dim2);

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Now begin the actual matrix construction
  for (int t=0;t<dim1;t++){
    for (int u=0;u<dim2;u++) {

      // note, whenever computing mixed-rank terms, e.g. 20_00, we
      // check that dim2 is large enough to handle the Tab(u,t) case
      // as well as the Tab(t,u) case.

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	Tab(t,u) = 1.0*f0;
      }
      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3) { // 1*_00 (chg-dip)
	if (u==0) {
	  Tab(t,u) = rA[t-1]*f1;
	  if (t < dim2) 
	    Tab(u,t) = rB[t-1]*f1;
	}
	// Dipole-dipole terms - Ltot = 2
	else if (u>=1 && u<=3) {// 1*_1* (dip-dip) 
	  Tab(t,u) = (3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))*f2;
	  if (t < dim2) 
	    Tab(u,t) = (3*rB[t-1]*rA[u-1] + cAB(u-1,t-1))*f2;
	}
      }
      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00 
	Tab(t,u) = 0.5*(3*raz*raz - 1.0)*f2;
	if (t < dim2) Tab(u,t) = 0.5*(3*rbz*rbz - 1.0)*f2;
      }
      else if (t==5 && u==0) { // 21c_00
	Tab(t,u) = sqrt(3.0)*rax*raz*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rbz*f2;
      }
      else if (t==6 && u==0) { // 21s_00
	Tab(t,u) = sqrt(3.0)*ray*raz*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rby*rbz*f2;
      }
      else if (t==7 && u==0) { // 22c_00
	Tab(t,u) = sqrt(3.0)/2.0*(rax*rax-ray*ray)*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(rbx*rbx-rby*rby)*f2;
      }
      else if (t==8 && u==0) { // 22s_00
	Tab(t,u) = sqrt(3.0)*rax*ray*f2;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*rbx*rby*f2;
      }
      
      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	Tab(t,u) = 0.5*(5*pow(raz,3.0) - 3*raz)*f3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(5*pow(rbz,3.0) - 3*rbz)*f3;
      }
      else if (t==10 && u==0) { // 31c_00
	Tab(t,u) = sqrt(6.0)/4.0*rax*(5*raz*raz - 1.0)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rbx*(5*rbz*rbz - 1.0)*f3;
      }
      else if (t==11 && u==0) { // 31s_00
	Tab(t,u) = sqrt(6.0)/4.0*ray*(5*raz*raz - 1.0)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*rby*(5*rbz*rbz - 1.0)*f3;
      }
      else if (t==12 && u==0) { // 32c_00
	Tab(t,u) = sqrt(15.0)/2.0*raz*(rax*rax - ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*rbz*(rbx*rbx - rby*rby)*f3;
      }
      else if (t==13 && u==0) { // 32s_00
	Tab(t,u) = sqrt(15.0)*rax*ray*raz*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*rbx*rby*rbz*f3;
      }
      else if (t==14 && u==0) { // 33c_00
	Tab(t,u) = sqrt(10.0)/4.0*rax*(rax*rax - 3*ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rbx*(rbx*rbx - 3*rby*rby)*f3;
      }
      else if (t==15 && u==0) { // 33s_00
	Tab(t,u) = sqrt(10.0)/4.0*ray*(3*rax*rax - ray*ray)*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*rby*(3*rbx*rbx - rby*rby)*f3;
      }

      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	Tab(t,u) = 0.5*(15*raz*raz*rB[u-1] + 6*raz*cAB(2,u-1) 
			- 3*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(15*rbz*rbz*rA[u-1] + 6*rbz*cAB(u-1,2) 
			- 3*rA[u-1])*f3;
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	Tab(t,u) = sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz 
			    + 5*rax*raz*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rbx*cAB(u-1,2) + cAB(u-1,0)*rbz 
			    + 5*rbx*rbz*rA[u-1])*f3;
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	Tab(t,u) = sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz 
			    + 5*ray*raz*rB[u-1])*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(rby*cAB(u-1,2) + cAB(u-1,1)*rbz 
			    + 5*rby*rbz*rA[u-1])*f3;
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	Tab(t,u) = sqrt(3.0)/2.0*(5*(rax*rax-ray*ray)*rB[u-1] 
				+ 2*rax*cAB(0,u-1) - 2*ray*cAB(1,u-1))*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)/2.0*(5*(rbx*rbx-rby*rby)*rA[u-1] 
				+ 2*rbx*cAB(u-1,0) - 2*rby*cAB(u-1,1))*f3;
      }
      else if (t==8 && u>=1 && u<=3) { // 22s_1*
	Tab(t,u) = sqrt(3.0)*(5*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			    + ray*cAB(0,u-1))*f3;
	if (t < dim2) 
	  Tab(u,t) = sqrt(3.0)*(5*rbx*rby*rA[u-1] + rbx*cAB(u-1,1) 
			    + rby*cAB(u-1,0))*f3;
      }

      // Charge-hexadecapole terms - Ltot = 4   (untested)
      else if (t==16 && u==0) { // 40_00
	Tab(t,u) = 0.125*(35*pow(raz,4.0) - 30*raz*raz + 3)*f4;
	if (t < dim2) 
	  Tab(u,t) = 0.125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)*f4;
      }
      else if (t==17 && u==0) { // 41c_00
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*rax*pow(raz,3.0) - 3*rax*raz)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*rbx*pow(rbz,3.0) - 3*rbx*rbz)*f4;
      }
      else if (t==18 && u==0) { // 41s_00
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*ray*pow(raz,3.0) - 3*ray*raz)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*rby*pow(rbz,3.0) - 3*rby*rbz)*f4;
      }
      else if (t==19 && u==0) { // 42c_00
	Tab(t,u) = sqrt(5.0)/4.0*(7.0*raz*raz - 1.0)*(rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/4.0*(7.0*rbz*rbz - 1.0)*(rbx*rbx-rby*rby)*f4;
      }
      else if (t==20 && u==0) { // 42s_00
	Tab(t,u) = sqrt(5.0)/2.0*(7.0*raz*raz - 1.0)*rax*ray*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(5.0)/2.0*(7.0*rbz*rbz - 1.0)*rbx*rby*f4;
      }
      else if (t==21 && u==0) { // 43c_00
	Tab(t,u) = sqrt(70.0)/4.0*rax*raz*(rax*rax-3*ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rbx*rbz*(rbx*rbx-3*rby*rby)*f4;
      }
      else if (t==22 && u==0) { // 43s_00
	Tab(t,u) = sqrt(70.0)/4.0*ray*raz*(3*rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(70.0)/4.0*rby*rbz*(3*rbx*rbx-rby*rby)*f4;
      }
      else if (t==23 && u==0) { // 44c_00
	Tab(t,u) = sqrt(35.0)/8.0*(pow(rax,4.0) - 6*rax*rax*ray*ray
				 + pow(ray,4.0))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/8.0*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				   + pow(rby,4.0))*f4;
      }
      else if (t==24 && u==0) { // 44s_00
	Tab(t,u) = sqrt(35.0)/2.0*rax*ray*(rax*rax-ray*ray)*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(35.0)/2.0*rbx*rby*(rbx*rbx-rby*rby)*f4;
      }

      // Dipole-Octopole - Ltot = 4
      else if (t==9 && u>=1 && u<=3) {// 30_1*
	Tab(t,u) = 0.5*(35*pow(raz,3.0)*rB[u-1] + 15*raz*raz*cAB(2,u-1) 
			- 15*raz*rB[u-1] - 3*cAB(2,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = 0.5*(35*pow(rbz,3.0)*rA[u-1] + 15*rbz*rbz*cAB(u-1,2) 
			  - 15*rbz*rA[u-1] - 3*cAB(u-1,2))*f4;
      }
      else if (t==10 && u>=1 && u<=3) {// 31c_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*rax*raz*raz*rB[u-1] + 5*raz*raz*cAB(0,u-1) 
				+ 10*rax*raz*cAB(2,u-1) - 5*rax*rB[u-1] 
				- cAB(0,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rbx*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,0)
				  + 10*rbx*rbz*cAB(u-1,2) - 5*rbx*rA[u-1] 
				  - cAB(u-1,0))*f4;
      }
      else if (t==11 && u>=1 && u<=3) {// 31s_1*
	Tab(t,u) = sqrt(6.0)/4.0*(35*ray*raz*raz*rB[u-1] + 5*raz*raz*cAB(1,u-1) 
				+ 10*ray*raz*cAB(2,u-1)	- 5*ray*rB[u-1] 
				- cAB(1,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(6.0)/4.0*(35*rby*rbz*rbz*rA[u-1] + 5*rbz*rbz*cAB(u-1,1)
				  + 10*rby*rbz*cAB(u-1,2) - 5*rby*rA[u-1] 
				  - cAB(u-1,1))*f4;
      }
      else if (t==12 && u>=1 && u<=3) {// 32c_1*
	Tab(t,u) = sqrt(15.0)/2.0*((rax*rax-ray*ray)*(7.0*raz*rB[u-1] + cAB(2,u-1))
				 + 2*raz*(rax*cAB(0,u-1) 
					  - ray*cAB(1,u-1)))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)/2.0*((rbx*rbx-rby*rby)*(7.0*rbz*rA[u-1] 
						      + cAB(u-1,2)) 
				   + 2*rbz*(rbx*cAB(u-1,0) 
					    - rby*cAB(u-1,1)))*f4;
      }
      else if (t==13 && u>=1 && u<=3) {// 32s_1*
	Tab(t,u) = sqrt(15.0)*(rax*ray*(7.0*raz*rB[u-1] + cAB(2,u-1))
			     + raz*(rax*cAB(1,u-1) + ray*cAB(0,u-1)))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(15.0)*(rbx*rby*(7.0*rbz*rA[u-1] + cAB(u-1,2))
			       + rbz*(rbx*cAB(u-1,1) + rby*cAB(u-1,0)))*f4;
	    }
      else if (t==14 && u>=1 && u<=3) {// 33c_1*
	Tab(t,u) = sqrt(10.0)/4.0*(7.0*pow(rax,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				 - 21*rax*ray*ray*rB[u-1]
				 - 6*rax*ray*cAB(1,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(7.0*pow(rbx,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,0)
				 - 21*rbx*rby*rby*rA[u-1]
				   - 6*rbx*rby*cAB(u-1,1))*f4;
      }
      else if (t==15 && u>=1 && u<=3) {// 33s_1*
	Tab(t,u) = sqrt(10.0)/4.0*(-7.0*pow(ray,3.0)*rB[u-1] 
				 + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				 + 21*rax*rax*ray*rB[u-1]
				 + 6*rax*ray*cAB(0,u-1))*f4;
	if (t < dim2) 
	  Tab(u,t) = sqrt(10.0)/4.0*(-7.0*pow(rby,3.0)*rA[u-1] 
				   + 3*(rbx*rbx-rby*rby)*cAB(u-1,1)
				   + 21*rbx*rbx*rby*rA[u-1]
				   + 6*rbx*rby*cAB(u-1,0))*f4;
      }
    
      // Quadrupole-quadrupole terms - Ltot = 4
      // Note: Stone's book arranged these with u >= t, but I use
      // the opposite convention.  So my Tab(u,t) = his Tab(t,u),
      // and vice-versa.
      else if (t==4 && u==4) { // 20_20
	Tab(t,u) = 0.75*(35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz 
			 + 20*raz*rbz*cAB(2,2) + 2*cAB(2,2)*cAB(2,2)
			 + 1)*f4;
      }
      else if (t==5 && u==4) { // 20_21c
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				+ 10*raz*rbx*cAB(2,2) + 10*raz*rbz*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,2))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*raz - 5*rax*raz 
				+ 10*rbz*rax*cAB(2,2) + 10*rbz*raz*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(2,2))*f4;
      }
      else if (t==6 && u==4) { // 20_21s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rby*rbz - 5*rby*rbz 
				+ 10*raz*rby*cAB(2,2) + 10*raz*rbz*cAB(2,1) 
				+ 2*cAB(2,1)*cAB(2,2))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*ray*raz - 5*ray*raz 
				+ 10*rbz*ray*cAB(2,2) + 10*rbz*raz*cAB(1,2) 
				+ 2*cAB(1,2)*cAB(2,2))*f4;
      }
      else if (t==7 && u==4) { // 20_22c 
	Tab(u,t) = sqrt(3.0)/4.0*(35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby 
				- 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*cAB(2,0) 
				- 20*raz*rby*cAB(2,1) + 2*cAB(2,0)*cAB(2,0)
				- 2*cAB(2,1)*cAB(2,1))*f4;
	Tab(t,u) = sqrt(3.0)/4.0*(35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray 
				- 5*rax*rax + 5*ray*ray + 20*rbz*rax*cAB(0,2) 
				- 20*rbz*ray*cAB(1,2) + 2*cAB(0,2)*cAB(0,2)
				- 2*cAB(1,2)*cAB(1,2))*f4;
      }
      else if (t==8 && u==4) { // 20_22s
	Tab(u,t) = sqrt(3.0)/2.0*(35*raz*raz*rbx*rby - 5*rbx*rby 
				+ 10*raz*rbx*cAB(2,1) + 10*raz*rby*cAB(2,0) 
				+ 2*cAB(2,0)*cAB(2,1))*f4;
	Tab(t,u) = sqrt(3.0)/2.0*(35*rbz*rbz*rax*ray - 5*rax*ray
				+ 10*rbz*rax*cAB(1,2) + 10*rbz*ray*cAB(0,2) 
				+ 2*cAB(0,2)*cAB(1,2))*f4;
      }
      else if (t==5 && u==5) { // 21c_21c
	Tab(t,u) = (35*rax*raz*rbx*rbz + 5*rax*rbx*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,0) + 5*raz*rbx*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,0) + cAB(0,0)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,0))*f4;
      }
      else if (t==6 && u==5) { // 21c_21s
	Tab(u,t) = (35*rax*raz*rby*rbz + 5*rax*rby*cAB(2,2) 
		    + 5*rax*rbz*cAB(2,1) + 5*raz*rby*cAB(0,2) 
		    + 5*raz*rbz*cAB(0,1) + cAB(0,1)*cAB(2,2) 
		    + cAB(0,2)*cAB(2,1))*f4;
	Tab(t,u) = (35*rbx*rbz*ray*raz + 5*rbx*ray*cAB(2,2) 
		    + 5*rbx*raz*cAB(1,2) + 5*rbz*ray*cAB(2,0) 
		    + 5*rbz*raz*cAB(1,0) + cAB(1,0)*cAB(2,2) 
		    + cAB(2,0)*cAB(1,2))*f4;
      }
      else if (t==7 && u==5) { // 21c_22c 
	Tab(u,t) = 0.5*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby 
			+ 10*rax*rbx*cAB(2,0) - 10*rax*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(0,0) - 10*raz*rby*cAB(0,1)
			+ 2*cAB(0,0)*cAB(2,0) - 2*cAB(0,1)*cAB(2,1))*f4;
	Tab(t,u) = 0.5*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray 
			+ 10*rbx*rax*cAB(0,2) - 10*rbx*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,0) - 10*rbz*ray*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,2) - 2*cAB(1,0)*cAB(1,2))*f4;
      }
      else if (t==8 && u==5) { // 21c_22s
	Tab(u,t) = (35*rax*raz*rbx*rby + 5*rax*rbx*cAB(2,1) 
		    + 5*rax*rby*cAB(2,0) + 5*raz*rbx*cAB(0,1) 
		    + 5*raz*rby*cAB(0,0) + cAB(0,0)*cAB(2,1) 
		    + cAB(0,1)*cAB(2,0))*f4;
	Tab(t,u) = (35*rbx*rbz*rax*ray + 5*rbx*rax*cAB(1,2)
		    + 5*rbx*ray*cAB(0,2) + 5*rbz*rax*cAB(1,0) 
		    + 5*rbz*ray*cAB(0,0) + cAB(0,0)*cAB(1,2) 
		    + cAB(1,0)*cAB(0,2))*f4;
      }
      else if (t==6 && u==6) { // 21s_21s
	Tab(t,u) = (35*ray*raz*rby*rbz + 5*ray*rby*cAB(2,2) 
		    + 5*ray*rbz*cAB(2,1) + 5*raz*rby*cAB(1,2) 
		    + 5*raz*rbz*cAB(1,1) + cAB(1,1)*cAB(2,2) 
		    + cAB(1,2)*cAB(2,1))*f4;
      }
      else if (t==7 && u==6) { // 21s_22c
	Tab(u,t) = 0.5*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby 
			+ 10*ray*rbx*cAB(2,0) - 10*ray*rby*cAB(2,1) 
			+ 10*raz*rbx*cAB(1,0) - 10*raz*rby*cAB(1,1) 
			+ 2*cAB(1,0)*cAB(2,0) - 2*cAB(1,1)*cAB(2,1))*f4;
	Tab(t,u) = 0.5*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray 
			+ 10*rby*rax*cAB(0,2) - 10*rby*ray*cAB(1,2) 
			+ 10*rbz*rax*cAB(0,1) - 10*rbz*ray*cAB(1,1)
			+ 2*cAB(0,1)*cAB(0,2) - 2*cAB(1,1)*cAB(1,2))*f4;
      }
      else if (t==8 && u==6) { // 21s_22s       
	Tab(u,t) = (35*ray*raz*rbx*rby + 5*ray*rbx*cAB(2,1) 
		    + 5*ray*rby*cAB(2,0) + 5*raz*rbx*cAB(1,1) 
		    + 5*raz*rby*cAB(1,0) + cAB(1,0)*cAB(2,1) 
		    + cAB(1,1)*cAB(2,0))*f4;
	Tab(t,u) = (35*rby*rbz*rax*ray + 5*rby*rax*cAB(1,2) 
		    + 5*rby*ray*cAB(0,2) + 5*rbz*rax*cAB(1,1) 
		    + 5*rbz*ray*cAB(0,1) + cAB(0,1)*cAB(1,2) 
		    + cAB(1,1)*cAB(0,2))*f4;
      }
      else if (t==7 && u==7) { // 22c_22c
	Tab(t,u) = 0.25*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby 
			 - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby 
			 + 20*rax*rbx*cAB(0,0) - 20*rax*rby*cAB(0,1)
			 - 20*ray*rbx*cAB(1,0) + 20*ray*rby*cAB(1,1) 
			 + 2*cAB(0,0)*cAB(0,0) - 2*cAB(0,1)*cAB(0,1)
			 - 2*cAB(1,0)*cAB(1,0) + 2*cAB(1,1)*cAB(1,1))*f4;
      }
      else if (t==8 && u==7) { // 22c_22s
	Tab(u,t) = 0.5*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby 
			+ 10*rax*rbx*cAB(0,1) + 10*rax*rby*cAB(0,0) 
			- 10*ray*rbx*cAB(1,1) - 10*ray*rby*cAB(1,0) 
			+ 2*cAB(0,0)*cAB(0,1) - 2*cAB(1,0)*cAB(1,1))*f4;
	Tab(t,u) = 0.5*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			+ 10*rbx*rax*cAB(1,0) + 10*rbx*ray*cAB(0,0) 
			- 10*rby*rax*cAB(1,1) - 10*rby*ray*cAB(0,1) 
			+ 2*cAB(0,0)*cAB(1,0) - 2*cAB(0,1)*cAB(1,1))*f4;
      }
      else if (t==8 && u==8) { // 22s_22s
	Tab(t,u) = (35*rax*ray*rbx*rby + 5*rax*rbx*cAB(1,1) 
		    + 5*rax*rby*cAB(1,0) + 5*ray*rbx*cAB(0,1) 
		    + 5*ray*rby*cAB(0,0) + cAB(0,0)*cAB(1,1) 
		    + cAB(0,1)*cAB(1,0))*f4;
      }
      else  {
	if (t > u) {
	  //printf("WARNING:: Multipole interactions of type %s...%s not implemented.\n",
	  // types[t].c_str(),types[u].c_str());
	  Tab(t,u) = 0.0;
	  if (t < dim2) 
	    Tab(u,t) = 0.0;
	}
	//exit(1);
      }
    }
  }

  // Divide by 4*pi*epsilon
  Tab.Scale(1.0/perm);
  
  // Undo transposes, if necessary
  if (nQB > nQA) {
    Tab.Transpose();
    Vector tmp(rA);
    rA = rB;
    rB = tmp;
    cAB.Transpose();
  }

  delete [] eA;
  delete [] eB;

  /*
  if (atom_index==1 && other.atom_index==1 && if_damp==0) {
    Tab.Print("Tab inside atom.C");
  }d
  */
   
  //Tab.Print("TabDirect inside atom.C");

  return Tab;

}

// Computes the gradient of the Tab matrix for electrostatic
// interactions. See Stone's "Theory of Intermolecular Forces" book
// and Mol. Phys. 82, 411 (1994) for details.  Note, however, that I
// have remapped those expressions onto the variables used here for Tab.
Matrix Atom::BuildInteractionMatrixGradient(Matrix Tab, Matrix thisRotMat, 
					    Atom& other, Matrix otherRotMat,
					    double beta_damp) {
  /*
    
  dT/dX, where X is any Cartesian displacement of one of the two atoms
  (xA, yA, zA, xB, yB, zB) is given by:

  dT/dX = (dT/dq)*(dq/dX)

  where q is a list of coordinates in which Tab is expressed.  In
  particular, we focus on R*R, rax, ray, raz, rbx, rby, rbz.  The
  cAB(*,*) derivatives all end up with zero contribution, because
  d(cAB(r,s))/dX = 0 for all r and s.

  We build dT/dq and dq/dX separately, and then contract them to form dT/dX.

  The variables q are ordered as:
  0 - R*R, 1 - rax, 2 - ray, 3 - raz, 4 - rbx, 5 - rby, 6 - rbz

  In order to apply damping (as we do for induction energies), we end up needing
  the undifferentiated & *undamped* Tab matrix.  If no damping is applied, the
  Tab matrix is still used, but gives zero contribution.

  */


  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm =4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 

  // Allocate storage for dTdX, dTdq, and dq/dX
  int size_q = 7; // number of variables q with which we differentiate Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix dTdX(nQA*nQB,6); // 6 cartesians describe the two atoms.
  Matrix dTdq(nQA*nQB,size_q);
  Matrix dqdX(size_q,6);

  dTdX.Set(); dTdq.Set(); dqdX.Set();

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);

  // Find RAB = RA - RB and the distance;
  double Rnorm = GetInterAtomicDistance(other)*AngToBohr;
  // Predefine Rnorm^x here for simplicity in later equations
  double Rnorm2 = Rnorm*Rnorm;
  double Rnorm3 = Rnorm2*Rnorm;
  double Rnorm4 = Rnorm3*Rnorm;
  double Rnorm5 = Rnorm4*Rnorm;
  double Rnorm6 = Rnorm5*Rnorm;
  double Rnorm7 = Rnorm6*Rnorm;

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  Vector eAB(RB);
  eAB -= RA;
  Vector R = eAB;
  eAB.Normalize();

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = thisRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = otherRotMat.MatrixTimesVector(unit_mat.GetColumnVector(i));
    //eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  }

  // Define rA, rB, cAB
  // rA = eA dot eAB... component of eA lying along A->B axis
  // rB = eA dot eBA = eA dot (-eAB) ... component of eB lying along A->B axis
  // cAB(i,j) = rAi dot rBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(eAB);
    rB[i] = -1.0*eB[i].DotProduct(eAB);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-7.
  Vector damp_grad(6); // d(damp)/dR
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;
  for (int n=1;n<=5;n++) {
    if (if_damp) {
      // damping factor
      damp[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm, 0); 
      // d(damp)/dR
      damp_grad[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm, 1); 
    }
    else {
      damp[n] = 1.0; 
      damp_grad[n] = 0.0;
    }
  }
  //printf("Rnorm = %.8f\n",Rnorm);
  //damp.Print("damping factors\n");
  //damp_grad.Print("d(damp)/dR");

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  double cxx = cAB(0,0);
  double cxy = cAB(0,1);
  double cxz = cAB(0,2);
  double cyx = cAB(1,0);
  double cyy = cAB(1,1);
  double cyz = cAB(1,2);
  double czx = cAB(2,0);
  double czy = cAB(2,1);
  double czz = cAB(2,2);

  // Build dq/dX
  for (int X=0;X<3;X++) { // loop over derivs in xyz directions

    // d(R*R)/dX - row 0
    dqdX(0,X) = -2*R[X]; // w.r.t. Xa
    dqdX(0,X+3) = 2*R[X]; // w.r.t. Xb
  
    for (int q=0;q<3;q++) { // loop over components of rA and rB
      // d(ra*)/dX - rows 1-3
      dqdX(q+1,X) = rA[q]*R[X]/Rnorm2 - eA[q][X]/Rnorm; // w.r.t. Xa
      dqdX(q+1,X+3) = -rA[q]*R[X]/Rnorm2 + eA[q][X]/Rnorm; // w.r.t. Xb
      
      // d(rbx)/dX - rows 4-6
      dqdX(q+4,X) = rB[q]*R[X]/Rnorm2 + eB[q][X]/Rnorm; // w.r.t. Xa
      dqdX(q+4,X+3) = -rB[q]*R[X]/Rnorm2 - eB[q][X]/Rnorm; // w.r.t. Xb
    }
  }


  //if (atom_index==1 && other.atom_index==1 && if_damp) 
  //  dqdX.Print("\ndamped dq/dX.  q = rows, xyz of each atom = cols.");
  //else if (atom_index==1 && other.atom_index==1 && !if_damp)
  //  dqdX.Print("\n dq/dX.  q = rows, xyz of each atom = cols.");
	   

  // Build dT/dq

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Note, derivatives w.r.t. cAB are neglected here, even when
  // nonzero, since the corresponding dq/dX terms are zero.  In other
  // words, they don't contribute to the overall gradient.

  // Now begin the actual matrix construction
  for (int t=0;t<nQA;t++){
    for (int u=0;u<nQB;u++) {
      int tu = u*nQA + t; // define a compound index for (t,u)
      int ut = t*nQA + u; // (u,t) index, which is used in a few places

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	// only 1 non-zero contribution
	dTdq(tu,0) = -1.0/(2.0*Rnorm3)*damp[1] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[1];  // q=R*R
      }

      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3 && u==0) { // 1*_00 (dip-chg)
	dTdq(tu,0) = -rA[t-1]/Rnorm4*damp[2] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[2]; // q=R*R
	dTdq(tu,t) = 1.0/Rnorm2*damp[2] ;  // q = rA*
      }
      else if (t==0 && u>=1 && u<=3) { // 00_1* (chg-dip)
	dTdq(tu,0) = -rB[u-1]/Rnorm4*damp[2] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[2]; // q=R*R
	dTdq(tu,u+3) = 1.0/Rnorm2*damp[2] ;  // q = rB*
      }

      // Dipole-dipole terms - Ltot = 2
      else if (t>=1 && t<=3 && u>=1 && u<=3) {// 1*_1* (dip-dip) 
	dTdq(tu,0) = -1.5*(3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,t) = 3.0*rB[u-1]/Rnorm3*damp[3]; // q = rA*
	dTdq(tu,u+3) = 3.0*rA[t-1]/Rnorm3*damp[3]; // q = rB*
      }

      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00
	dTdq(tu,0) = -0.75*(3*raz*raz-1)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,3) = 3.0*raz/Rnorm3*damp[3]; // q=raz
      }
      else if (t==5 && u==0) { // 21c_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*rax*raz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*raz/Rnorm3*damp[3]; // q=rax
	dTdq(tu,3) = sqrt(3.0)*rax/Rnorm3*damp[3]; // q=raz
      }
      else if (t==6 && u==0) { // 21s_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*ray*raz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,2) = sqrt(3.0)*raz/Rnorm3*damp[3]; // q=ray
	dTdq(tu,3) = sqrt(3.0)*ray/Rnorm3*damp[3]; // q=raz
      }
      else if (t==7 && u==0) { // 22c_00
	dTdq(tu,0) = -0.75*sqrt(3.0)*(rax*rax-ray*ray)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*rax/Rnorm3*damp[3]; // q=rax
	dTdq(tu,2) = -sqrt(3.0)*ray/Rnorm3*damp[3];// q=ray
      }
      else if (t==8 && u==0) { // 22s_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*rax*ray/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*ray/Rnorm3*damp[3]; // q=rax
	dTdq(tu,2) = sqrt(3.0)*rax/Rnorm3*damp[3];// q=ray
      }
      
      else if (t==0 && u==4) { // 00_20
	dTdq(tu,0) = -0.75*(3*rbz*rbz-1)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,6) = 3.0*rbz/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==5) { // 00_21c
	dTdq(tu,0) = -1.5*sqrt(3.0)*rbx*rbz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rbz/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,6) = sqrt(3.0)*rbx/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==6) { // 00_21s
	dTdq(tu,0) = -1.5*sqrt(3.0)*rby*rbz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,5) = sqrt(3.0)*rbz/Rnorm3*damp[3]; // q=rby
	dTdq(tu,6) = sqrt(3.0)*rby/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==7) { // 00_22c
	dTdq(tu,0) = -0.75*sqrt(3.0)*(rbx*rbx-rby*rby)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rbx/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,5) = -sqrt(3.0)*rby/Rnorm3*damp[3];// q=rby
      }
      else if (t==0 && u==8) { // 00_22s
	dTdq(tu,0) = -1.5*sqrt(3.0)*rbx*rby/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rby/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,5) = sqrt(3.0)*rbx/Rnorm3*damp[3];// q=rby
      }

      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	dTdq(tu,0) = -(5.0*raz*raz*raz - 3.0*raz)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,3) = 0.5*(15*raz*raz - 3.0)/Rnorm4*damp[4]; // q=raz
      }
      else if (t==10 && u==0) { // 31c_00
	dTdq(tu,0) = -0.5*sqrt(6.0)*rax*(5.0*raz*raz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(6.0)*(5.0*raz*raz - 1.0)/Rnorm4*damp[4]; // q=rax
	dTdq(tu,3) = 2.5*sqrt(6.0)*rax*raz/Rnorm4*damp[4]; // q=raz
      }
      else if (t==11 && u==0) { // 31s_00
	dTdq(tu,0) = -0.5*sqrt(6.0)*ray*(5.0*raz*raz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(6.0)*(5.0*raz*raz - 1.0)/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = 2.5*sqrt(6.0)*ray*raz/Rnorm4*damp[4]; // q=raz
      }
      else if (t==12 && u==0) { // 32c_00
	dTdq(tu,0) = -sqrt(15.0)*raz*(rax*rax-ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*rax*raz/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = -sqrt(15.0)*ray*raz/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(15.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=raz
      }
      else if (t==13 && u==0) { // 32s_00
	dTdq(tu,0) = -2*sqrt(15.0)*rax*ray*raz/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*ray*raz/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = sqrt(15.0)*rax*raz/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = sqrt(15.0)*rax*ray/Rnorm4*damp[4]; // q=raz
      }
      else if (t==14 && u==0) { // 33c_00
	dTdq(tu,0) = -0.5*sqrt(10.0)*(rax*rax*rax 
				    - 3*rax*ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(10.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(10.0)*rax*ray/Rnorm4*damp[4]; // q=ray	
      }
      else if (t==15 && u==0) { // 33s_00
	dTdq(tu,0) = -0.5*sqrt(10.0)*ray*(3*rax*rax - ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(10.0)*rax*ray/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(10.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=ray
      }

      else if (t==0 && u==9) { // 00_30
	dTdq(tu,0) = -(5.0*rbz*rbz*rbz - 3*rbz)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,6) = 0.5*(15.0*rbz*rbz - 3)/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==10) { // 00_31c
	dTdq(tu,0) = -0.5*sqrt(6.0)*rbx*(5*rbz*rbz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(6.0)*(5*rbz*rbz - 1.0)/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,6) = 2.5*sqrt(6.0)*rbx*rbz/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==11) { // 00_31s
	dTdq(tu,0) = -0.5*sqrt(6.0)*rby*(5*rbz*rbz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(6.0)*(5*rbz*rbz - 1.0)/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = 2.5*sqrt(6.0)*rby*rbz/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==12) { // 00_32c
	dTdq(tu,0) = -sqrt(15.0)*rbz*(rbx*rbx-rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*rbx*rbz/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = -sqrt(15.0)*rby*rbz/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(15.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==13) { // 00_32s
	dTdq(tu,0) = -2*sqrt(15.0)*rbx*rby*rbz/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*rby*rbz/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = sqrt(15.0)*rbx*rbz/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = sqrt(15.0)*rbx*rby/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==14) { // 00_33c
	dTdq(tu,0) = -0.5*sqrt(10.0)*(rbx*rbx*rbx 
				    - 3*rbx*rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(10.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(10.0)*rbx*rby/Rnorm4*damp[4]; // q=rby	
      }
      else if (t==0 && u==15) { // 00_33s
	dTdq(tu,0) = -0.5*sqrt(10.0)*rby*(3*rbx*rbx - rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(10.0)*rbx*rby/Rnorm4*damp[4]; // q=rby
	dTdq(tu,5) = 0.75*sqrt(10.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbx
      }


      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	dTdq(tu,0) = -(15.0*raz*raz*rB[u-1] + 6.0*raz*cAB(2,u-1) - 
		       3.0*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,3) = 0.5*(30.0*raz*rB[u-1] + 6.0*cAB(2,u-1))/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = 0.5*(15.0*raz*raz - 3.0)/Rnorm4*damp[4]; // q = rB*
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz +
				 5.0*rax*raz*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*(cAB(2,u-1) + 5.0*raz*rB[u-1])/Rnorm4*damp[4]; // q=rax
	dTdq(tu,3) = sqrt(3.0)*(cAB(0,u-1) + 5.0*rax*rB[u-1])/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = sqrt(3.0)*5.0*rax*raz/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz +
				 5.0*ray*raz*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,2) = sqrt(3.0)*(cAB(2,u-1) + 5.0*raz*rB[u-1])/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = sqrt(3.0)*(cAB(1,u-1) + 5.0*ray*rB[u-1])/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = sqrt(3.0)*5.0*ray*raz/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	dTdq(tu,0) = -sqrt(3.0)*(5.0*(rax*rax-ray*ray)*rB[u-1] + 
			       2.0*rax*cAB(0,u-1) 
			       - 2.0*ray*cAB(1,u-1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(3.0)*(10*rax*rB[u-1] 
				  + 2.0*cAB(0,u-1))/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(3.0)*(-10*ray*rB[u-1] 
				  - 2.0*cAB(1,u-1))/Rnorm4*damp[4]; // q=ray
	dTdq(tu,u+3) = 0.5*5.0*sqrt(3.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==8 && u>=1 && u<=3) { // 21s_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(5.0*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
				 + ray*cAB(0,u-1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*(5.0*ray*rB[u-1] + cAB(1,u-1))/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = sqrt(3.0)*(5.0*rax*rB[u-1] + cAB(0,u-1))/Rnorm4*damp[4]; // q=ray
	dTdq(tu,u+3) = sqrt(3.0)*5.0*rax*ray/Rnorm4*damp[4]; // q=rB*
      }

      else if (u==4 && t>=1 && t<=3) { // 1*_20
	dTdq(tu,0) = -(15*rbz*rbz*rA[t-1] + 6*rbz*cAB(t-1,2) - 
		       3*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,6) = 0.5*(30*rbz*rA[t-1] + 6*cAB(t-1,2))/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = 0.5*(15*rbz*rbz - 3)/Rnorm4*damp[4]; // q = rA*
      }
      else if (u==5 && t>=1 && t<=3) { // 1*_21c
	dTdq(tu,0) = -2*sqrt(3.0)*(rbx*cAB(t-1,2) + cAB(t-1,0)*rbz +
				 5*rbx*rbz*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*(cAB(t-1,2) + 5*rbz*rA[t-1])/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,6) = sqrt(3.0)*(cAB(t-1,0) + 5*rbx*rA[t-1])/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = sqrt(3.0)*5*rbx*rbz/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==6 && t>=1 && t<=3) { // 1*_21s
	dTdq(tu,0) = -2*sqrt(3.0)*(rby*cAB(t-1,2) + cAB(t-1,1)*rbz +
				 5*rby*rbz*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,5) = sqrt(3.0)*(cAB(t-1,2) + 5*rbz*rA[t-1])/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = sqrt(3.0)*(cAB(t-1,1) + 5*rby*rA[t-1])/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = sqrt(3.0)*5*rby*rbz/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==7 && t>=1 && t<=3) { // 1*_22c
	dTdq(tu,0) = -sqrt(3.0)*(5*(rbx*rbx-rby*rby)*rA[t-1] + 
			       2*rbx*cAB(t-1,0) 
			       - 2*rby*cAB(t-1,1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(3.0)*(10*rbx*rA[t-1] 
				  + 2*cAB(t-1,0))/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(3.0)*(-10*rby*rA[t-1] 
				  - 2*cAB(t-1,1))/Rnorm4*damp[4]; // q=rb
	dTdq(tu,t) = 0.5*5*sqrt(3.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==8 && t>=1 && t<=3) { // 1*_21s
	dTdq(tu,0) = -2*sqrt(3.0)*(5*rbx*rby*rA[t-1] + rbx*cAB(t-1,1) 
				 + rby*cAB(t-1,0))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*(5*rby*rA[t-1] + cAB(t-1,1))/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = sqrt(3.0)*(5*rbx*rA[t-1] + cAB(t-1,0))/Rnorm4*damp[4]; // q=rby
	dTdq(tu,t) = sqrt(3.0)*5*rbx*rby/Rnorm4*damp[4]; // q=rA*
      }
      
      // Charge-hexadecapole terms - Ltot = 4
      else if (t==16 && u==0) { // 40_00
	dTdq(tu,0) = -0.3125*(35*pow(raz,4.0) - 30*raz*raz + 3)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 2.5*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==17 && u==0) { // 41c_00
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rax*pow(raz,3.0) 
				      - 3*rax*raz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(10.0)*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.75*sqrt(10.0)*(7.0*rax*raz*raz - rax)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==18 && u==0) { // 41s_00
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*ray*pow(raz,3.0) 
				      - 3*ray*raz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(10.0)*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.75*sqrt(10.0)*(7.0*ray*raz*raz - ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==19 && u==0) { // 42c_00
	dTdq(tu,0) = -5/8.0*sqrt(5.0)*(7.0*raz*raz - 1.0)
	  *(rax*rax - ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*rax/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -0.5*sqrt(5.0)*(7.0*raz*raz-1)*ray/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 3.5*sqrt(5.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==20 && u==0) { // 42s_00
	dTdq(tu,0) = -1.25*sqrt(5.0)*(7.0*raz*raz-1)*rax*ray/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*ray/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*rax/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 7.0*sqrt(5.0)*rax*ray*raz/Rnorm5*damp[5]; // q=ray
      }
      else if (t==21 && u==0) { // 43c_00
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rax*raz
	  *(rax*rax-3*ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(70.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(70.0)*rax*ray*raz/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.25*sqrt(70.0)*rax*(rax*rax-3*ray*ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==22 && u==0) { // 43s_00
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*ray*raz
	  *(3*rax*rax-ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(70.0)*rax*ray*raz/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(70.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.25*sqrt(70.0)*ray*(3*rax*rax-ray*ray)/Rnorm5*damp[5]; // q=raz
	
      }
      else if (t==23 && u==0) { // 44c_00
	dTdq(tu,0) = -5/16.0*sqrt(35.0)*(pow(rax,4.0) - 6*rax*rax*ray*ray 
				       + pow(ray,4.0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(35.0)*(pow(rax,3.0) - 3*rax*ray*ray)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(35.0)*(pow(ray,3.0) - 3*rax*rax*ray)/Rnorm5*damp[5]; // q=ray
	
      }
      else if (t==24 && u==0) { // 44s_00
	dTdq(tu,0) = -1.25*sqrt(35.0)*rax*ray*(rax*rax-ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(35.0)*(3*rax*rax*ray - pow(ray,3.0))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(35.0)*(pow(rax,3.0) - 3*ray*ray*rax)/Rnorm5*damp[5]; // q=ray
      }


      else if (t==0 && u==16) { // 00_40
	dTdq(tu,0) = -0.3125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,6) = 2.5*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==17) { // 00_41c
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rbx*pow(rbz,3.0) 
				      - 3*rbx*rbz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(10.0)*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.75*sqrt(10.0)*(7.0*rbx*rbz*rbz - rbx)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==18) { // 00_41s
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rby*pow(rbz,3.0) 
				      - 3*rby*rbz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(10.0)*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.75*sqrt(10.0)*(7*rby*rbz*rbz - rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==19) { // 00_42c
	dTdq(tu,0) = -5/8.0*sqrt(5.0)*(7*rbz*rbz - 1.0)
	  *(rbx*rbx - rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rbx/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -0.5*sqrt(5.0)*(7*rbz*rbz-1)*rby/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 3.5*sqrt(5.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==20) { // 00_42s
	dTdq(tu,0) = -1.25*sqrt(5.0)*(7*rbz*rbz-1)*rbx*rby/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rby/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rbx/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 7*sqrt(5.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rby
      }
      else if (t==0 && u==21) { // 00_43c
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rbx*rbz
	  *(rbx*rbx-3*rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(70.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(70.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.25*sqrt(70.0)*rbx*(rbx*rbx-3*rby*rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==22) { // 00_43s
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rby*rbz
	  *(3*rbx*rbx-rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(70.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.75*sqrt(70.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.25*sqrt(70.0)*rby*(3*rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbz
	
      }
      else if (t==0 && u==23) { // 00_44c
	dTdq(tu,0) = -5/16.0*sqrt(35.0)*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				       + pow(rby,4.0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(35.0)*(pow(rbx,3.0) - 3*rbx*rby*rby)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(35.0)*(pow(rby,3.0) - 3*rbx*rbx*rby)/Rnorm5*damp[5]; // q=rby
	
      }
      else if (t==0 && u==24) { // 00_44s
	dTdq(tu,0) = -1.25*sqrt(35.0)*rbx*rby*(rbx*rbx-rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(35.0)*(3*rbx*rbx*rby - pow(rby,3.0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(35.0)*(pow(rbx,3.0) - 3*rby*rby*rbx)/Rnorm5*damp[5]; // q=rby
      }

      // Dipole-octupole terms - Ltot = 4
      else if (t==9 && u>=1 && u<=3) { // 30_1*
	dTdq(tu,0) = -1.25*(35*pow(raz,3.0)*rB[u-1] 
			    + 15*raz*raz*cAB(2,u-1)
			    - 15*raz*rB[u-1] - 3*cAB(2,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 7.5*(7*raz*raz*rB[u-1] + 2*raz*cAB(2,u-1) 
			  - rB[u-1])/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 2.5*(7*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=rB*
      } 
      else if (t==10 && u>=1 && u<=3) { // 31c_1*
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rax*raz*raz*rB[u-1]
				     + 5*raz*raz*cAB(0,u-1)
				     + 10*rax*raz*cAB(2,u-1)
				     - 5*rax*rB[u-1] 
				     - cAB(0,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(6.0)*(35*raz*raz*rB[u-1] + 10*raz*cAB(2,u-1) 
				   - 5*rB[u-1])/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.5*sqrt(6.0)*(35*rax*raz*rB[u-1] + 5*raz*cAB(0,u-1)
				  + 5*rax*cAB(2,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 0.25*sqrt(6.0)*(35*rax*raz*raz - 5*rax)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==11 && u>=1 && u<=3) { // 31s_1*
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*ray*raz*raz*rB[u-1]
				     + 5*raz*raz*cAB(1,u-1)
				     + 10*ray*raz*cAB(2,u-1)
				     - 5*ray*rB[u-1] 
				     - cAB(1,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(6.0)*(35*raz*raz*rB[u-1] + 10*raz*cAB(2,u-1) 
				   - 5*rB[u-1])/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(6.0)*(35*ray*raz*rB[u-1] + 5*raz*cAB(1,u-1)
				  + 5*ray*cAB(2,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 0.25*sqrt(6.0)*(35*ray*raz*raz - 5*ray)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==12 && u>=1 && u<=3) { // 32c_1*
	dTdq(tu,0) = -1.25*sqrt(15.0)*((rax*rax-ray*ray)
				     *(7*raz*rB[u-1] + cAB(2,u-1)) 
				     + 2*raz*(rax*cAB(0,u-1) 
					      - ray*cAB(1,u-1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*(rax*(7*raz*rB[u-1] + cAB(2,u-1))
			       + raz*cAB(0,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -sqrt(15.0)*(ray*(7*raz*rB[u-1] + cAB(2,u-1))
				+ raz*cAB(1,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(15.0)*(7*(rax*rax-ray*ray)*rB[u-1] 
				   + 2*(rax*cAB(0,u-1) 
					- ray*cAB(1,u-1)))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 3.5*sqrt(15.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; //q=rb*
	
      }
      else if (t==13 && u>=1 && u<=3) { // 32s_1*
	dTdq(tu,0) = -2.5*sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
				    + raz*(rax*cAB(1,u-1) 
					   + ray*cAB(0,u-1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*(ray*(7*raz*rB[u-1] + cAB(2,u-1)) 
			       + raz*cAB(1,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = sqrt(15.0)*(rax*(7*raz*rB[u-1] + cAB(2,u-1)) 
			       + raz*cAB(0,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = sqrt(15.0)*(7*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			       + ray*cAB(0,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 7*sqrt(15.0)*rax*ray*raz/Rnorm5*damp[5]; //q=rb*
      }
      else if (t==14 && u>=1 && u<=3) { // 33c_1*
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7*pow(rax,3.0)*rB[u-1] 
				      + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				      - 21*rax*ray*ray*rB[u-1]
				      - 6*rax*ray*cAB(1,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(10.0)*(7*rax*rax*rB[u-1] 
				    + 2*rax*cAB(0,u-1) 
				    - 7*ray*ray*rB[u-1] 
				    - 2*ray*cAB(1,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(10.0)*(ray*cAB(0,u-1) + 7*rax*ray*rB[u-1] 
				    + rax*cAB(1,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,u+3) = 1.75*sqrt(10.0)*(pow(rax,3.0) 
				      - 3*rax*ray*ray)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==15 && u>=1 && u<=3) { // 33s_1*
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(-7*pow(ray,3.0)*rB[u-1] 
				      + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				      + 21*rax*rax*ray*rB[u-1]
				      + 6*rax*ray*cAB(0,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(10.0)*(rax*cAB(1,u-1) 
				   + 7*rax*ray*rB[u-1] 
				   + ray*cAB(0,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(10.0)*(-7*ray*ray*rB[u-1] 
				    - 2*ray*cAB(1,u-1)
				    + 7*rax*rax*rB[u-1]
				    + 2*rax*cAB(0,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,u+3) = 1.75*sqrt(10.0)*(3*rax*rax*ray 
				      - pow(rax,3.0))/Rnorm5*damp[5]; // q=rb*
      }

      else if (u==9 && t>=1 && t<=3) { // 1*_30
	dTdq(tu,0) = -1.25*(35*pow(rbz,3.0)*rA[t-1] 
			    + 15*rbz*rbz*cAB(t-1,2)
			    - 15*rbz*rA[t-1] - 3*cAB(t-1,2))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,6) = 7.5*(7*rbz*rbz*rA[t-1] + 2*rbz*cAB(t-1,2) 
			  - rA[t-1])/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 2.5*(7*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=ra*
      } 
      else if (u==10 && t>=1 && t<=3) { // 1*_31c
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rbx*rbz*rbz*rA[t-1]
				     + 5*rbz*rbz*cAB(t-1,0)
				     + 10*rbx*rbz*cAB(t-1,2)
				     - 5*rbx*rA[t-1] 
				     - cAB(t-1,0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(6.0)*(35*rbz*rbz*rA[t-1] + 10*rbz*cAB(t-1,2) 
				   - 5*rA[t-1])/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.5*sqrt(6.0)*(35*rbx*rbz*rA[t-1] + 5*rbz*cAB(t-1,0)
				  + 5*rbx*cAB(t-1,2))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 0.25*sqrt(6.0)*(35*rbx*rbz*rbz - 5*rbx)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==11 && t>=1 && t<=3) { // 1*_31s
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rby*rbz*rbz*rA[t-1]
				     + 5*rbz*rbz*cAB(t-1,1)
				     + 10*rby*rbz*cAB(t-1,2)
				     - 5*rby*rA[t-1] 
				     - cAB(t-1,1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(6.0)*(35*rbz*rbz*rA[t-1] + 10*rbz*cAB(t-1,2) 
				   - 5*rA[t-1])/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(6.0)*(35*rby*rbz*rA[t-1] + 5*rbz*cAB(t-1,1)
				  + 5*rby*cAB(t-1,2))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 0.25*sqrt(6.0)*(35*rby*rbz*rbz - 5*rby)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==12 && t>=1 && t<=3) { // 1*_32c
	dTdq(tu,0) = -1.25*sqrt(15.0)*((rbx*rbx-rby*rby)
				     *(7*rbz*rA[t-1] + cAB(t-1,2)) 
				     + 2*rbz*(rbx*cAB(t-1,0) 
					      - rby*cAB(t-1,1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*(rbx*(7*rbz*rA[t-1] + cAB(t-1,2))
			       + rbz*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -sqrt(15.0)*(rby*(7*rbz*rA[t-1] + cAB(t-1,2))
				+ rbz*cAB(t-1,1))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(15.0)*(7*(rbx*rbx-rby*rby)*rA[t-1] 
				   + 2*(rbx*cAB(t-1,0) 
					- rby*cAB(t-1,1)))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 3.5*sqrt(15.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; //q=ra*
	
      }
      else if (u==13 && t>=1 && t<=3) { // 1*_32s
	dTdq(tu,0) = -2.5*sqrt(15.0)*(rbx*rby*(7*rbz*rA[t-1] + cAB(t-1,2))
				    + rbz*(rbx*cAB(t-1,1) 
					   + rby*cAB(t-1,0)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*(rby*(7*rbz*rA[t-1] + cAB(t-1,2)) 
			       + rbz*cAB(t-1,1))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = sqrt(15.0)*(rbx*(7*rbz*rA[t-1] + cAB(t-1,2)) 
			       + rbz*cAB(t-1,0))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = sqrt(15.0)*(7*rbx*rby*rA[t-1] + rbx*cAB(t-1,1) 
			       + rby*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 7*sqrt(15.0)*rbx*rby*rbz/Rnorm5*damp[5]; //q=ra*
      }
      else if (u==14 && t>=1 && t<=3) { // 1*_33c
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7*pow(rbx,3.0)*rA[t-1] 
				      + 3*(rbx*rbx-rby*rby)*cAB(t-1,0)
				      - 21*rbx*rby*rby*rA[t-1]
				      - 6*rbx*rby*cAB(t-1,1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(10.0)*(7*rbx*rbx*rA[t-1] 
				    + 2*rbx*cAB(t-1,0) 
				    - 7*rby*rby*rA[t-1] 
				    - 2*rby*cAB(t-1,1))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(10.0)*(rby*cAB(t-1,0) + 7*rbx*rby*rA[t-1] 
				    + rbx*cAB(1,t-1))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,t) = 1.75*sqrt(10.0)*(pow(rbx,3.0) 
				    - 3*rbx*rby*rby)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==15 && t>=1 && t<=3) { // 1*_33s
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(-7*pow(rby,3.0)*rA[t-1] 
				      + 3*(rbx*rbx-rby*rby)*cAB(t-1,1)
				      + 21*rbx*rbx*rby*rA[t-1]
				      + 6*rbx*rby*cAB(t-1,0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(10.0)*(rbx*cAB(t-1,1) 
				   + 7*rbx*rby*rA[t-1] 
				   + rby*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.75*sqrt(10.0)*(-7*rby*rby*rA[t-1] 
				    - 2*rby*cAB(t-1,1)
				    + 7*rbx*rbx*rA[t-1]
				    + 2*rbx*cAB(t-1,0))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,t) = 1.75*sqrt(10.0)*(3*rbx*rbx*rby 
				    - pow(rby,3.0))/Rnorm5*damp[5]; // q=ra*
      }



      // Quadrupole-quadrupole terms - Ltot = 4
      else if (t==4 && u==4) { // 20_20
	dTdq(tu,0) = -15/8.0*(35*raz*raz*rbz*rbz - 5*raz*raz 
			      - 5*rbz*rbz + 20*raz*rbz*czz 
			      + 2*czz*czz + 1)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.75*(70*raz*rbz*rbz - 10*raz 
			   + 20*rbz*czz)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,6) = 0.75*(70*raz*raz*rbz - 10*rbz 
			   + 20*raz*czz)/Rnorm5*damp[5]; // q=rbz
      }

      else if (t==4 && u==5) { // 20_21c
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				    + 10*raz*rbx*czz + 10*raz*rbz*czx
				    + 2*czx*czz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rbx*rbz + 10*rbx*czz 
				  + 10*rbz*czx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*sqrt(3.0)*(35*raz*raz*rbz - 5*rbz 
				  + 10*raz*czz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.5*sqrt(3.0)*(35*raz*raz*rbx - 5*rbx 
				  + 10*raz*czx)/Rnorm5*damp[5]; // q=rbz

	// 21c_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*rax*raz - 5*rax*raz 
				    + 10*rbz*rax*czz + 10*rbz*raz*cxz
				    + 2*cxz*czz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*rax*raz + 10*rax*czz 
				  + 10*raz*cxz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*sqrt(3.0)*(35*rbz*rbz*raz - 5*raz 
				  + 10*rbz*czz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,3) = 0.5*sqrt(3.0)*(35*rbz*rbz*rax - 5*rax 
				  + 10*rbz*cxz)/Rnorm5*damp[5]; // q=raz
      }

      else if (t==4 && u==6) { // 20_21s
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rby*rbz - 5*rby*rbz 
				    + 10*raz*rby*czz + 10*raz*rbz*czy
				    + 2*czy*czz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rby*rbz + 10*rby*czz 
				  + 10*rbz*czy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = 0.5*sqrt(3.0)*(35*raz*raz*rbz - 5*rbz 
				  + 10*raz*czz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(3.0)*(35*raz*raz*rby - 5*rby
				  + 10*raz*czy)/Rnorm5*damp[5]; // q=rbz

	// 21s_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*ray*raz - 5*ray*raz 
				    + 10*rbz*ray*czz + 10*rbz*raz*cyz
				    + 2*cyz*czz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*ray*raz + 10*ray*czz 
				  + 10*raz*cyz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,2) = 0.5*sqrt(3.0)*(35*rbz*rbz*raz - 5*raz 
				  + 10*rbz*czz)/Rnorm5*damp[5]; // q=ray
	dTdq(ut,3) = 0.5*sqrt(3.0)*(35*rbz*rbz*ray - 5*ray
				  + 10*rbz*cyz)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==4 && u==7) { // 20_22c
	dTdq(tu,0) = -0.625*sqrt(3.0)*(35*raz*raz*rbx*rbx 
				     - 35*raz*raz*rby*rby 
				     - 5*rbx*rbx + 5*rby*rby 
				     + 20*raz*rbx*czx
				     - 20*raz*rby*czy + 2*czx*czx
				     - 2*czy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.25*sqrt(3.0)*(70*raz*rbx*rbx 
				   - 70*raz*rby*rby
				   + 20*rbx*czx 
				   - 20*rby*czy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.25*sqrt(3.0)*(70*raz*raz*rbx - 10*rbx 
				   + 20*raz*czx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.25*sqrt(3.0)*(-70*raz*raz*rby + 10*rby
				   - 20*raz*czy)/Rnorm5*damp[5]; // q=rby
	
	// 22c_20
	dTdq(ut,0) = -0.625*sqrt(3.0)*(35*rbz*rbz*rax*rax 
				     - 35*rbz*rbz*ray*ray 
				     - 5*rax*rax + 5*ray*ray 
				     + 20*rbz*rax*cxz
				     - 20*rbz*ray*cyz + 2*cxz*cxz
				     - 2*cyz*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.25*sqrt(3.0)*(70*rbz*rax*rax 
				   - 70*rbz*ray*ray
				   + 20*rax*cxz 
				   - 20*ray*cyz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.25*sqrt(3.0)*(70*rbz*rbz*rax - 10*rax 
				   + 20*rbz*cxz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.25*sqrt(3.0)*(-70*rbz*rbz*ray + 10*ray
				   - 20*rbz*cyz)/Rnorm5*damp[5]; // q=ray

      }
      else if (t==4 && u==8) { // 20_22s
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rbx*rby 
				    - 5*rbx*rby
				    + 10*raz*rbx*czy
				    + 10*raz*rby*czx 
				    + 2*czx*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rbx*rby + 10*rbx*czy 
				  + 10*rby*czx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*sqrt(3.0)*(35*raz*raz*rby - 5*rby 
				  + 10*raz*czy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(3.0)*(35*raz*raz*rbx - 5*rbx 
				  + 10*raz*czx)/Rnorm5*damp[5]; // q=rby

	// 22s_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*rax*ray 
				    - 5*rax*ray
				    + 10*rbz*rax*cyz
				    + 10*rbz*ray*cxz 
				    + 2*cxz*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*rax*ray + 10*rax*cyz 
				  + 10*ray*cxz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*sqrt(3.0)*(35*rbz*rbz*ray - 5*ray 
				  + 10*rbz*cyz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*sqrt(3.0)*(35*rbz*rbz*rax - 5*rax 
				  + 10*rbz*cxz)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==5 && u==5) { // 21c_21c
	dTdq(tu,0) = -2.5*(35*rax*raz*rbx*rbz + 5*rax*rbx*czz 
			   + 5*rax*rbz*czx + 5*raz*rbx*cxz 
			   + 5*raz*rbz*cxx + cxx*czz 
			   + cxz*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rbx*rbz + 5*rbx*czz + 5*rbz*czx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rbx*rbz + 5*rbx*cxz + 5*rbz*cxx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*rax*raz*rbz + 5*rax*czz + 5*raz*cxz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = (35*rax*raz*rbx + 5*rax*czx + 5*raz*cxx)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==5 && u==6) { // 21c_21s
	dTdq(tu,0) = -2.5*(35*rax*raz*rby*rbz + 5*rax*rby*czz 
			   + 5*rax*rbz*czy + 5*raz*rby*cxz 
			   + 5*raz*rbz*cxy + cxy*czz 
			   + cxz*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rby*rbz + 5*rby*czz + 5*rbz*czy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rby*rbz + 5*rby*cxz + 5*rbz*cxy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = (35*rax*raz*rbz + 5*rax*czz + 5*raz*cxz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = (35*rax*raz*rby + 5*rax*czy + 5*raz*cxy)/Rnorm5*damp[5]; // q=rbz

	// 21s_21c
	dTdq(ut,0) = -2.5*(35*rbx*rbz*ray*raz + 5*rbx*ray*czz 
			   + 5*rbx*raz*cyz + 5*rbz*ray*czx 
			   + 5*rbz*raz*cyx + cyx*czz 
			   + czx*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = (35*rbz*ray*raz + 5*ray*czz + 5*raz*cyz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = (35*rbx*ray*raz + 5*ray*czx + 5*raz*cyx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,2) = (35*rbx*rbz*raz + 5*rbx*czz + 5*rbz*czx)/Rnorm5*damp[5]; // q=ray
	dTdq(ut,3) = (35*rbx*rbz*ray + 5*rbx*cyz + 5*rbz*cyx)/Rnorm5*damp[5]; // q=raz
      }
      
      else if (t==5 && u==7) { // 21c_22c
	dTdq(tu,0) = -1.25*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby
			    + 10*rax*rbx*czx - 10*rax*rby*czy 
			    + 10*raz*rbx*cxx - 10*raz*rby*cxy
			    + 2*cxx*czx - 2*cxy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*(35*raz*rbx*rbx - 35*raz*rby*rby 
			  + 10*rbx*czx - 10*rby*czy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.5*(35*rax*rbx*rbx - 35*rax*rby*rby
			  + 10*rbx*cxx - 10*rby*cxy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*(70*rax*raz*rbx + 10*rax*czx 
			  + 10*raz*cxx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*(-70*rax*raz*rby - 10*rax*czy 
			  - 10*raz*cxy)/Rnorm5*damp[5]; // q=rby
	
	// 22c_21c     
	dTdq(ut,0) = -1.25*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray
			    + 10*rbx*rax*cxz - 10*rbx*ray*cyz 
			    + 10*rbz*rax*cxx - 10*rbz*ray*cyx
			    + 2*cxx*cxz - 2*cyx*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = 0.5*(35*rbz*rax*rax - 35*rbz*ray*ray 
			  + 10*rax*cxz - 10*ray*cyz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = 0.5*(35*rbx*rax*rax - 35*rbx*ray*ray
			  + 10*rax*cxx - 10*ray*cyx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*(70*rbx*rbz*rax + 10*rbx*cxz 
			  + 10*rbz*cxx)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*(-70*rbx*rbz*ray - 10*rbx*cyz 
			  - 10*rbz*cyx)/Rnorm5*damp[5]; // q=ray
	
      }
      else if (t==5 && u==8) { // 21c_22s
	dTdq(tu,0) = -2.5*(35*rax*raz*rbx*rby + 5*rax*rbx*czy
			   + 5*rax*rby*czx + 5*raz*rbx*cxy
			   + 5*raz*rby*cxx + cxx*czy 
			   + cxy*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rbx*rby + 5*rbx*czy + 5*rby*czx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rbx*rby + 5*rbx*cxy + 5*rby*cxx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*rax*raz*rby + 5*rax*czy + 5*raz*cxy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*rax*raz*rbx + 5*rax*czx + 5*raz*cxx)/Rnorm5*damp[5]; // q=rby
	
	// 22s_21c     
	dTdq(ut,0) = -2.5*(35*rbx*rbz*rax*ray + 5*rbx*rax*cyz
			   + 5*rbx*ray*cxz + 5*rbz*rax*cyx
			   + 5*rbz*ray*cxx + cxx*cyz 
			   + cyx*cxz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = (35*rbz*rax*ray + 5*rax*cyz + 5*ray*cxz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = (35*rbx*rax*ray + 5*rax*cyx + 5*ray*cxx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = (35*rbx*rbz*ray + 5*rbx*cyz + 5*rbz*cyx)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = (35*rbx*rbz*rax + 5*rbx*cxz + 5*rbz*cxx)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==6 && u==6) { // 21s_21s
	dTdq(tu,0) = -2.5*(35*ray*raz*rby*rbz + 5*ray*rby*czz
			   + 5*ray*rbz*czy + 5*raz*rby*cyz
			   + 5*raz*rbz*cyy + cyy*czz 
			   + cyz*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = (35*raz*rby*rbz + 5*rby*czz + 5*rbz*czy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = (35*ray*rby*rbz + 5*rby*cyz + 5*rbz*cyy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = (35*ray*raz*rbz + 5*ray*czz + 5*raz*cyz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = (35*ray*raz*rby + 5*ray*czy + 5*raz*cyy)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==6 && u==7) { // 21s_22c
	dTdq(tu,0) = -1.25*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby
			    + 10*ray*rbx*cxz - 10*ray*rby*czy
			    + 10*raz*rbx*cyx - 10*raz*rby*cyy +
			    2*cyx*czx - 2*cyy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.5*(35*raz*rbx*rbx - 35*raz*rby*rby
			  + 10*rbx*czx - 10*rby*czy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*(35*ray*rbx*rbx - 35*ray*rby*rby
			  + 10*rbx*cyx - 10*rby*cyy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*(70*ray*raz*rbx + 10*ray*czx 
			  + 10*raz*cyx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -0.5*(70*ray*raz*rby + 10*ray*czy
			   + 10*raz*cyy)/Rnorm5*damp[5]; // q=rby

	// 22c_21s
	dTdq(ut,0) = -1.25*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray
			    + 10*rby*rax*czx - 10*rby*ray*cyz
			    + 10*rbz*rax*cxy - 10*rbz*ray*cyy +
			    2*cxy*cxz - 2*cyy*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,5) = 0.5*(35*rbz*rax*rax - 35*rbz*ray*ray
			  + 10*rax*cxz - 10*ray*cyz)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,6) = 0.5*(35*rby*rax*rax - 35*rby*ray*ray
			  + 10*rax*cxy - 10*ray*cyy)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*(70*rby*rbz*rax + 10*rby*cxz 
			  + 10*rbz*cxy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = -0.5*(70*rby*rbz*ray + 10*rby*cyz
			   + 10*rbz*cyy)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==6 && u==8) { // 21s_22s
	dTdq(tu,0) = -2.5*(35*ray*raz*rbx*rby + 5*ray*rbx*czy
			   + 5*ray*rby*czx + 5*raz*rbx*cyy
			   + 5*raz*rby*cyx + cyx*czy 
			   + cyy*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = (35*raz*rbx*rby + 5*rbx*czy + 5*rby*czx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = (35*ray*rbx*rby + 5*rbx*cyy + 5*rby*cyx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*ray*raz*rby + 5*ray*czy + 5*raz*cyy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*ray*raz*rbx + 5*ray*czx + 5*raz*cyx)/Rnorm5*damp[5]; // q=rby

	// 22s_21s
	dTdq(ut,0) = -2.5*(35*rby*rbz*rax*ray + 5*rby*rax*cyz
			   + 5*rby*ray*cxz + 5*rbz*rax*cyy
			   + 5*rbz*ray*cxy + cxy*cyz 
			   + cyy*cxz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,5) = (35*rbz*rax*ray + 5*rax*cyz + 5*ray*cxz)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,6) = (35*rby*rax*ray + 5*rax*cyy + 5*ray*cxy)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = (35*rby*rbz*ray + 5*rby*cyz + 5*rbz*cyy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = (35*rby*rbz*rax + 5*rby*cxz + 5*rbz*cxy)/Rnorm5*damp[5]; // q=ray

      }
      else if (t==7 && u==7) { // 22c_22c
	dTdq(tu,0) = -0.625*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby
			     - 35*ray*ray*rbx*rbx
			     + 35*ray*ray*rby*rby
			     + 20*rax*rbx*cxx - 20*rax*rby*cxy
			     - 20*ray*rbx*cyx + 20*ray*rby*cyy
			     + 2*cxx*cxx - 2*cxy*cxy 
			     - 2*cyx*cyx + 2*cyy*cyy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*(70*rax*rbx*rbx - 70*rax*rby*rby
			   + 20*rbx*cxx - 20*rby*cxy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.25*(-70*ray*rbx*rbx + 70*ray*rby*rby
			   - 20*rbx*cyx + 20*rby*cyy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = 0.25*(70*rax*rax*rbx - 70*ray*ray*rbx
			   + 20*rax*cxx - 20*ray*cyx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.25*(-70*rax*rax*rby + 70*ray*ray*rby
			   - 20*rax*cxy + 20*ray*cyy)/Rnorm5*damp[5]; // q=rby
      }
      else if (t==7 && u==8) { // 22c_22s
	dTdq(tu,0) = -1.25*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby
			    + 10*rax*rbx*cxy + 10*rax*rby*cxx
			    - 10*ray*rbx*cyy - 10*ray*rby*cyx
			    + 2*cxx*cxy - 2*cyx*cyy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*(70*rax*rbx*rby + 10*rbx*cxy 
			  + 10*rby*cxx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*(-70*ray*rbx*rby - 10*rbx*cyy 
			  - 10*rby*cyx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = 0.5*(35*rax*rax*rby - 35*ray*ray*rby
			  + 10*rax*cxy - 10*ray*cyy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*(35*rax*rax*rbx - 35*ray*ray*rbx
			  + 10*rax*cxx - 10*ray*cyx)/Rnorm5*damp[5]; // q=rby

	// 22s_22c
	dTdq(ut,0) = -1.25*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			    + 10*rbx*rax*cyx + 10*rbx*ray*cxx
			    - 10*rby*rax*cyy - 10*rby*ray*cxy
			    + 2*cxx*cyx - 2*cxy*cyy)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = 0.5*(70*rbx*rax*ray + 10*rax*cyx 
			  + 10*ray*cxx)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,5) = 0.5*(-70*rby*rax*ray - 10*rax*cyy 
			  - 10*ray*cxy)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,1) = 0.5*(35*rbx*rbx*ray - 35*rby*rby*ray
			  + 10*rbx*cyx -10*rby*cyy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*(35*rbx*rbx*rax - 35*rby*rby*rax
			  + 10*rbx*cxx - 10*rby*cxy)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==8 && u==8) { // 22s_22s
	dTdq(tu,0) = -2.5*(35*rax*ray*rbx*rby + 5*rax*rbx*cyy
			   + 5*rax*rby*cyx + 5*ray*rbx*cxy
			   + 5*ray*rby*cxx + cxx*cyy 
			   + cxy*cyx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*ray*rbx*rby + 5*rbx*cyy + 5*rby*cyx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = (35*rax*rbx*rby + 5*rbx*cxy + 5*rby*cxx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = (35*rax*ray*rby + 5*rax*cyy + 5*ray*cxy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*rax*ray*rbx + 5*rax*cyx + 5*ray*cxx)/Rnorm5*damp[5]; // q=rby
      }



    }
  }

  // Build dT/dX = dT/dq * dq/dX;
  dTdX = dTdq.Multiply(dqdX);
  dTdX.Scale(1.0/perm);// Divide by 4*pi*epsilon
  
  //dTdX.Print("\ndT/dX");

  delete [] eA;
  delete [] eB;

  return dTdX;
}


//Outdated, now using rotation matrix instend of angle/axis notation

// Computes the gradient of the Tab matrix for electrostatic
// interactions. See Stone's "Theory of Intermolecular Forces" book
// and Mol. Phys. 82, 411 (1994) for details.  Note, however, that I
// have remapped those expressions onto the variables used here for Tab.
Matrix Atom::BuildInteractionMatrixGradient(Matrix Tab, Vector thisRotVec, 
					    double thisRotAng, 
					    Atom& other, Vector otherRotVec,
					    double otherRotAng, 
					    double beta_damp){
  /*
    
  dT/dX, where X is any Cartesian displacement of one of the two atoms
  (xA, yA, zA, xB, yB, zB) is given by:

  dT/dX = (dT/dq)*(dq/dX)

  where q is a list of coordinates in which Tab is expressed.  In
  particular, we focus on R*R, rax, ray, raz, rbx, rby, rbz.  The
  cAB(*,*) derivatives all end up with zero contribution, because
  d(cAB(r,s))/dX = 0 for all r and s.

  We build dT/dq and dq/dX separately, and then contract them to form dT/dX.

  The variables q are ordered as:
  0 - R*R, 1 - rax, 2 - ray, 3 - raz, 4 - rbx, 5 - rby, 6 - rbz

  In order to apply damping (as we do for induction energies), we end up needing
  the undifferentiated & *undamped* Tab matrix.  If no damping is applied, the
  Tab matrix is still used, but gives zero contribution.

  */


  // Start out with some preparation
  
  // Define the permativity constant 4*pi*epsilon. 
  // in units of hartrees per bohr  
  double perm =4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na)*HartreesToKJpermole;
  // in units of kJ/mol per bohr
  //double perm = 4*pi*epsilon*1000/(MtoAng*AngToBohr*ec*ec*Na); 

  // Allocate storage for dTdX, dTdq, and dq/dX
  int size_q = 7; // number of variables q with which we differentiate Tab
  int nQA = GetMultipoleMoments().GetLength();
  int nQB = other.GetMultipoleMoments().GetLength();

  Matrix dTdX(nQA*nQB,6); // 6 cartesians describe the two atoms.
  Matrix dTdq(nQA*nQB,size_q);
  Matrix dqdX(size_q,6);

  dTdX.Set(); dTdq.Set(); dqdX.Set();

  // Grab global position of each multipole expansion site, switch to a.u.
  Vector RA(xyz);
  Vector RB(other.xyz);

  RA.Scale(AngToBohr);
  RB.Scale(AngToBohr);

  // Find RAB = RA - RB and the distance;
  double Rnorm = GetInterAtomicDistance(other)*AngToBohr;
  // Predefine Rnorm^x here for simplicity in later equations
  double Rnorm2 = Rnorm*Rnorm;
  double Rnorm3 = Rnorm2*Rnorm;
  double Rnorm4 = Rnorm3*Rnorm;
  double Rnorm5 = Rnorm4*Rnorm;
  double Rnorm6 = Rnorm5*Rnorm;
  double Rnorm7 = Rnorm6*Rnorm;

  // Define eAB = (RB-RA)/norm(RB-RA), the unit vector from A -> B
  Vector eAB(RB);
  eAB -= RA;
  Vector R = eAB;
  eAB.Normalize();

  // Define some helpful geometric vectors.
  // eA, eB are the unit vectors defining the local coordinate systems
  // of A and B in terms of the global coordinates
  Matrix unit_mat(3,true); // 3x3 Identity matrix
  Vector unitX = unit_mat.GetColumnVector(1);
  Vector *eA, *eB;
  eA = new Vector[3];
  eB = new Vector[3];

  for (int i=0;i<3;i++) {
    eA[i].Initialize(3);
    eA[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*thisRotAng,thisRotVec);
    eA[i].Normalize();

    eB[i].Initialize(3);
    eB[i] = unit_mat.GetColumnVector(i).RotateAboutAxis3D(-1.0*otherRotAng,otherRotVec);
    eB[i].Normalize();
  }

  // Define rA, rB, cAB
  // rA = eA dot eAB... component of eA lying along A->B axis
  // rB = eA dot eBA = eA dot (-eAB) ... component of eB lying along A->B axis
  // cAB(i,j) = rAi dot rBj
  Vector rA(3), rB(3);
  Matrix cAB(3,3), tmpA(3,3), tmpB(3,3);

  for (int i=0;i<3;i++) {
    rA[i] = eA[i].DotProduct(eAB);
    rB[i] = -1.0*eB[i].DotProduct(eAB);  
    
    tmpA.SetColumnVector(eA[i],i);
    tmpB.SetColumnVector(eB[i],i);
  }

  cAB = tmpA.Multiply(tmpB,2); // cAB = tmpA'*tmpB

  // Damping factor - using Tang-Toennies damping factor.  Requires
  // parameter beta_damp that must be specified earlier.
  Vector damp(6); // for convenience, we use indexing 1-7.
  Vector damp_grad(6); // d(damp)/dR
  bool if_damp = true;
  if (beta_damp == -999.0)
    if_damp = false;
  for (int n=1;n<=5;n++) {
    if (if_damp) {
      // damping factor
      damp[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm, 0); 
      // d(damp)/dR
      damp_grad[n] = TangToenniesDampingFactor(n, beta_damp, Rnorm, 1); 
    }
    else {
      damp[n] = 1.0; 
      damp_grad[n] = 0.0;
    }
  }
  //printf("Rnorm = %.8f\n",Rnorm);
  //damp.Print("damping factors\n");
  //damp_grad.Print("d(damp)/dR");

  // make some handy aliases
  double rax = rA[0];
  double ray = rA[1];
  double raz = rA[2];
  double rbx = rB[0];
  double rby = rB[1];
  double rbz = rB[2];

  double cxx = cAB(0,0);
  double cxy = cAB(0,1);
  double cxz = cAB(0,2);
  double cyx = cAB(1,0);
  double cyy = cAB(1,1);
  double cyz = cAB(1,2);
  double czx = cAB(2,0);
  double czy = cAB(2,1);
  double czz = cAB(2,2);

  // Build dq/dX
  for (int X=0;X<3;X++) { // loop over derivs in xyz directions

    // d(R*R)/dX - row 0
    dqdX(0,X) = -2*R[X]; // w.r.t. Xa
    dqdX(0,X+3) = 2*R[X]; // w.r.t. Xb
  
    for (int q=0;q<3;q++) { // loop over components of rA and rB
      // d(ra*)/dX - rows 1-3
      dqdX(q+1,X) = rA[q]*R[X]/Rnorm2 - eA[q][X]/Rnorm; // w.r.t. Xa
      dqdX(q+1,X+3) = -rA[q]*R[X]/Rnorm2 + eA[q][X]/Rnorm; // w.r.t. Xb
      
      // d(rbx)/dX - rows 4-6
      dqdX(q+4,X) = rB[q]*R[X]/Rnorm2 + eB[q][X]/Rnorm; // w.r.t. Xa
      dqdX(q+4,X+3) = -rB[q]*R[X]/Rnorm2 - eB[q][X]/Rnorm; // w.r.t. Xb
    }
  }


  //if (atom_index==1 && other.atom_index==1 && if_damp) 
  //  dqdX.Print("\ndamped dq/dX.  q = rows, xyz of each atom = cols.");
  //else if (atom_index==1 && other.atom_index==1 && !if_damp)
  //  dqdX.Print("\n dq/dX.  q = rows, xyz of each atom = cols.");
	   

  // Build dT/dq

  /* A helpful key: indexing for t/u used here:
     0 - 00   
     1 - 1x    2 - 1y    3 - 1z   
     4 - 20    5 - 21c   6 - 21s   7 - 22c   8 - 22s
     9 - 30   10 - 31c  11 - 31s  12 - 32c  13 - 32s  14 - 33c  15 - 33s
     16 - 40  17 - 41c  18 - 41s  19 - 42c  20 - 42s  21 - 43c  22 - 43s  23 - 44c  24 - 44s
  */
  string types[25];
  types[0]="00";
  types[1]="1x"; types[2]="1y"; types[3]="1z";
  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; types[8]="22s";
  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; types[13]="32s"; 
  types[14]="33c"; types[15]="33s";
  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; types[20]="42s"; 
  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";

  // Note, derivatives w.r.t. cAB are neglected here, even when
  // nonzero, since the corresponding dq/dX terms are zero.  In other
  // words, they don't contribute to the overall gradient.

  // Now begin the actual matrix construction
  for (int t=0;t<nQA;t++){
    for (int u=0;u<nQB;u++) {
      int tu = u*nQA + t; // define a compound index for (t,u)
      int ut = t*nQA + u; // (u,t) index, which is used in a few places

      // Charge-charge term - Ltot = 0
      if (t==0 && u==0)  {// 00_00 (chg-chg)
	// only 1 non-zero contribution
	dTdq(tu,0) = -1.0/(2.0*Rnorm3)*damp[1] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[1];  // q=R*R
      }

      // Charge-dipole terms - Ltot = 1
      else if (t>=1 && t<=3 && u==0) { // 1*_00 (dip-chg)
	dTdq(tu,0) = -rA[t-1]/Rnorm4*damp[2] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[2]; // q=R*R
	dTdq(tu,t) = 1.0/Rnorm2*damp[2] ;  // q = rA*
      }
      else if (t==0 && u>=1 && u<=3) { // 00_1* (chg-dip)
	dTdq(tu,0) = -rB[u-1]/Rnorm4*damp[2] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[2]; // q=R*R
	dTdq(tu,u+3) = 1.0/Rnorm2*damp[2] ;  // q = rB*
      }

      // Dipole-dipole terms - Ltot = 2
      else if (t>=1 && t<=3 && u>=1 && u<=3) {// 1*_1* (dip-dip) 
	dTdq(tu,0) = -1.5*(3*rA[t-1]*rB[u-1] + cAB(t-1,u-1))/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,t) = 3.0*rB[u-1]/Rnorm3*damp[3]; // q = rA*
	dTdq(tu,u+3) = 3.0*rA[t-1]/Rnorm3*damp[3]; // q = rB*
      }

      // Charge-quadrupole terms - Ltot = 2
      else if (t==4 && u==0) { // 20_00
	dTdq(tu,0) = -0.75*(3*raz*raz-1)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,3) = 3.0*raz/Rnorm3*damp[3]; // q=raz
      }
      else if (t==5 && u==0) { // 21c_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*rax*raz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*raz/Rnorm3*damp[3]; // q=rax
	dTdq(tu,3) = sqrt(3.0)*rax/Rnorm3*damp[3]; // q=raz
      }
      else if (t==6 && u==0) { // 21s_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*ray*raz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,2) = sqrt(3.0)*raz/Rnorm3*damp[3]; // q=ray
	dTdq(tu,3) = sqrt(3.0)*ray/Rnorm3*damp[3]; // q=raz
      }
      else if (t==7 && u==0) { // 22c_00
	dTdq(tu,0) = -0.75*sqrt(3.0)*(rax*rax-ray*ray)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*rax/Rnorm3*damp[3]; // q=rax
	dTdq(tu,2) = -sqrt(3.0)*ray/Rnorm3*damp[3];// q=ray
      }
      else if (t==8 && u==0) { // 22s_00
	dTdq(tu,0) = -1.5*sqrt(3.0)*rax*ray/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*ray/Rnorm3*damp[3]; // q=rax
	dTdq(tu,2) = sqrt(3.0)*rax/Rnorm3*damp[3];// q=ray
      }
      
      else if (t==0 && u==4) { // 00_20
	dTdq(tu,0) = -0.75*(3*rbz*rbz-1)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,6) = 3.0*rbz/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==5) { // 00_21c
	dTdq(tu,0) = -1.5*sqrt(3.0)*rbx*rbz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rbz/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,6) = sqrt(3.0)*rbx/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==6) { // 00_21s
	dTdq(tu,0) = -1.5*sqrt(3.0)*rby*rbz/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,5) = sqrt(3.0)*rbz/Rnorm3*damp[3]; // q=rby
	dTdq(tu,6) = sqrt(3.0)*rby/Rnorm3*damp[3]; // q=rbz
      }
      else if (t==0 && u==7) { // 00_22c
	dTdq(tu,0) = -0.75*sqrt(3.0)*(rbx*rbx-rby*rby)/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rbx/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,5) = -sqrt(3.0)*rby/Rnorm3*damp[3];// q=rby
      }
      else if (t==0 && u==8) { // 00_22s
	dTdq(tu,0) = -1.5*sqrt(3.0)*rbx*rby/Rnorm5*damp[3] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[3]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*rby/Rnorm3*damp[3]; // q=rbx
	dTdq(tu,5) = sqrt(3.0)*rbx/Rnorm3*damp[3];// q=rby
      }

      // Charge-octupole terms - Ltot = 3
      else if (t==9 && u==0) { // 30_00
	dTdq(tu,0) = -(5.0*raz*raz*raz - 3.0*raz)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,3) = 0.5*(15*raz*raz - 3.0)/Rnorm4*damp[4]; // q=raz
      }
      else if (t==10 && u==0) { // 31c_00
	dTdq(tu,0) = -0.5*sqrt(6.0)*rax*(5.0*raz*raz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(6.0)*(5.0*raz*raz - 1.0)/Rnorm4*damp[4]; // q=rax
	dTdq(tu,3) = 2.5*sqrt(6.0)*rax*raz/Rnorm4*damp[4]; // q=raz
      }
      else if (t==11 && u==0) { // 31s_00
	dTdq(tu,0) = -0.5*sqrt(6.0)*ray*(5.0*raz*raz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(6.0)*(5.0*raz*raz - 1.0)/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = 2.5*sqrt(6.0)*ray*raz/Rnorm4*damp[4]; // q=raz
      }
      else if (t==12 && u==0) { // 32c_00
	dTdq(tu,0) = -sqrt(15.0)*raz*(rax*rax-ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*rax*raz/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = -sqrt(15.0)*ray*raz/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(15.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=raz
      }
      else if (t==13 && u==0) { // 32s_00
	dTdq(tu,0) = -2*sqrt(15.0)*rax*ray*raz/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*ray*raz/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = sqrt(15.0)*rax*raz/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = sqrt(15.0)*rax*ray/Rnorm4*damp[4]; // q=raz
      }
      else if (t==14 && u==0) { // 33c_00
	dTdq(tu,0) = -0.5*sqrt(10.0)*(rax*rax*rax 
				    - 3*rax*ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(10.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(10.0)*rax*ray/Rnorm4*damp[4]; // q=ray	
      }
      else if (t==15 && u==0) { // 33s_00
	dTdq(tu,0) = -0.5*sqrt(10.0)*ray*(3*rax*rax - ray*ray)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(10.0)*rax*ray/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(10.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=ray
      }

      else if (t==0 && u==9) { // 00_30
	dTdq(tu,0) = -(5.0*rbz*rbz*rbz - 3*rbz)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,6) = 0.5*(15.0*rbz*rbz - 3)/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==10) { // 00_31c
	dTdq(tu,0) = -0.5*sqrt(6.0)*rbx*(5*rbz*rbz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(6.0)*(5*rbz*rbz - 1.0)/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,6) = 2.5*sqrt(6.0)*rbx*rbz/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==11) { // 00_31s
	dTdq(tu,0) = -0.5*sqrt(6.0)*rby*(5*rbz*rbz - 1.0)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(6.0)*(5*rbz*rbz - 1.0)/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = 2.5*sqrt(6.0)*rby*rbz/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==12) { // 00_32c
	dTdq(tu,0) = -sqrt(15.0)*rbz*(rbx*rbx-rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*rbx*rbz/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = -sqrt(15.0)*rby*rbz/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(15.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==13) { // 00_32s
	dTdq(tu,0) = -2*sqrt(15.0)*rbx*rby*rbz/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*rby*rbz/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = sqrt(15.0)*rbx*rbz/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = sqrt(15.0)*rbx*rby/Rnorm4*damp[4]; // q=rbz
      }
      else if (t==0 && u==14) { // 00_33c
	dTdq(tu,0) = -0.5*sqrt(10.0)*(rbx*rbx*rbx 
				    - 3*rbx*rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(10.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(10.0)*rbx*rby/Rnorm4*damp[4]; // q=rby	
      }
      else if (t==0 && u==15) { // 00_33s
	dTdq(tu,0) = -0.5*sqrt(10.0)*rby*(3*rbx*rbx - rby*rby)/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(10.0)*rbx*rby/Rnorm4*damp[4]; // q=rby
	dTdq(tu,5) = 0.75*sqrt(10.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rbx
      }


      // Dipole-quadrupole terms - Ltot = 3
      else if (t==4 && u>=1 && u<=3) { // 20_1*
	dTdq(tu,0) = -(15.0*raz*raz*rB[u-1] + 6.0*raz*cAB(2,u-1) - 
		       3.0*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,3) = 0.5*(30.0*raz*rB[u-1] + 6.0*cAB(2,u-1))/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = 0.5*(15.0*raz*raz - 3.0)/Rnorm4*damp[4]; // q = rB*
      }
      else if (t==5 && u>=1 && u<=3) { // 21c_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(rax*cAB(2,u-1) + cAB(0,u-1)*raz +
				 5.0*rax*raz*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*(cAB(2,u-1) + 5.0*raz*rB[u-1])/Rnorm4*damp[4]; // q=rax
	dTdq(tu,3) = sqrt(3.0)*(cAB(0,u-1) + 5.0*rax*rB[u-1])/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = sqrt(3.0)*5.0*rax*raz/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==6 && u>=1 && u<=3) { // 21s_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(ray*cAB(2,u-1) + cAB(1,u-1)*raz +
				 5.0*ray*raz*rB[u-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,2) = sqrt(3.0)*(cAB(2,u-1) + 5.0*raz*rB[u-1])/Rnorm4*damp[4]; // q=ray
	dTdq(tu,3) = sqrt(3.0)*(cAB(1,u-1) + 5.0*ray*rB[u-1])/Rnorm4*damp[4]; // q=raz
	dTdq(tu,u+3) = sqrt(3.0)*5.0*ray*raz/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==7 && u>=1 && u<=3) { // 22c_1*
	dTdq(tu,0) = -sqrt(3.0)*(5.0*(rax*rax-ray*ray)*rB[u-1] + 
			       2.0*rax*cAB(0,u-1) 
			       - 2.0*ray*cAB(1,u-1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(3.0)*(10*rax*rB[u-1] 
				  + 2.0*cAB(0,u-1))/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(3.0)*(-10*ray*rB[u-1] 
				  - 2.0*cAB(1,u-1))/Rnorm4*damp[4]; // q=ray
	dTdq(tu,u+3) = 0.5*5.0*sqrt(3.0)*(rax*rax-ray*ray)/Rnorm4*damp[4]; // q=rB*
      }
      else if (t==8 && u>=1 && u<=3) { // 21s_1*
	dTdq(tu,0) = -2.0*sqrt(3.0)*(5.0*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
				 + ray*cAB(0,u-1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,1) = sqrt(3.0)*(5.0*ray*rB[u-1] + cAB(1,u-1))/Rnorm4*damp[4]; // q=rax
	dTdq(tu,2) = sqrt(3.0)*(5.0*rax*rB[u-1] + cAB(0,u-1))/Rnorm4*damp[4]; // q=ray
	dTdq(tu,u+3) = sqrt(3.0)*5.0*rax*ray/Rnorm4*damp[4]; // q=rB*
      }

      else if (u==4 && t>=1 && t<=3) { // 1*_20
	dTdq(tu,0) = -(15*rbz*rbz*rA[t-1] + 6*rbz*cAB(t-1,2) - 
		       3*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,6) = 0.5*(30*rbz*rA[t-1] + 6*cAB(t-1,2))/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = 0.5*(15*rbz*rbz - 3)/Rnorm4*damp[4]; // q = rA*
      }
      else if (u==5 && t>=1 && t<=3) { // 1*_21c
	dTdq(tu,0) = -2*sqrt(3.0)*(rbx*cAB(t-1,2) + cAB(t-1,0)*rbz +
				 5*rbx*rbz*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*(cAB(t-1,2) + 5*rbz*rA[t-1])/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,6) = sqrt(3.0)*(cAB(t-1,0) + 5*rbx*rA[t-1])/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = sqrt(3.0)*5*rbx*rbz/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==6 && t>=1 && t<=3) { // 1*_21s
	dTdq(tu,0) = -2*sqrt(3.0)*(rby*cAB(t-1,2) + cAB(t-1,1)*rbz +
				 5*rby*rbz*rA[t-1])/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,5) = sqrt(3.0)*(cAB(t-1,2) + 5*rbz*rA[t-1])/Rnorm4*damp[4]; // q=rby
	dTdq(tu,6) = sqrt(3.0)*(cAB(t-1,1) + 5*rby*rA[t-1])/Rnorm4*damp[4]; // q=rbz
	dTdq(tu,t) = sqrt(3.0)*5*rby*rbz/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==7 && t>=1 && t<=3) { // 1*_22c
	dTdq(tu,0) = -sqrt(3.0)*(5*(rbx*rbx-rby*rby)*rA[t-1] + 
			       2*rbx*cAB(t-1,0) 
			       - 2*rby*cAB(t-1,1))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(3.0)*(10*rbx*rA[t-1] 
				  + 2*cAB(t-1,0))/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(3.0)*(-10*rby*rA[t-1] 
				  - 2*cAB(t-1,1))/Rnorm4*damp[4]; // q=rb
	dTdq(tu,t) = 0.5*5*sqrt(3.0)*(rbx*rbx-rby*rby)/Rnorm4*damp[4]; // q=rA*
      }
      else if (u==8 && t>=1 && t<=3) { // 1*_21s
	dTdq(tu,0) = -2*sqrt(3.0)*(5*rbx*rby*rA[t-1] + rbx*cAB(t-1,1) 
				 + rby*cAB(t-1,0))/Rnorm6*damp[4] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[4]; // q=R*R
	dTdq(tu,4) = sqrt(3.0)*(5*rby*rA[t-1] + cAB(t-1,1))/Rnorm4*damp[4]; // q=rbx
	dTdq(tu,5) = sqrt(3.0)*(5*rbx*rA[t-1] + cAB(t-1,0))/Rnorm4*damp[4]; // q=rby
	dTdq(tu,t) = sqrt(3.0)*5*rbx*rby/Rnorm4*damp[4]; // q=rA*
      }
      
      // Charge-hexadecapole terms - Ltot = 4
      else if (t==16 && u==0) { // 40_00
	dTdq(tu,0) = -0.3125*(35*pow(raz,4.0) - 30*raz*raz + 3)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 2.5*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==17 && u==0) { // 41c_00
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rax*pow(raz,3.0) 
				      - 3*rax*raz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(10.0)*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.75*sqrt(10.0)*(7.0*rax*raz*raz - rax)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==18 && u==0) { // 41s_00
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*ray*pow(raz,3.0) 
				      - 3*ray*raz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(10.0)*(7.0*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.75*sqrt(10.0)*(7.0*ray*raz*raz - ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==19 && u==0) { // 42c_00
	dTdq(tu,0) = -5/8.0*sqrt(5.0)*(7.0*raz*raz - 1.0)
	  *(rax*rax - ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*rax/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -0.5*sqrt(5.0)*(7.0*raz*raz-1)*ray/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 3.5*sqrt(5.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==20 && u==0) { // 42s_00
	dTdq(tu,0) = -1.25*sqrt(5.0)*(7.0*raz*raz-1)*rax*ray/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*ray/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(5.0)*(7.0*raz*raz-1)*rax/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 7.0*sqrt(5.0)*rax*ray*raz/Rnorm5*damp[5]; // q=ray
      }
      else if (t==21 && u==0) { // 43c_00
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rax*raz
	  *(rax*rax-3*ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(70.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(70.0)*rax*ray*raz/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.25*sqrt(70.0)*rax*(rax*rax-3*ray*ray)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==22 && u==0) { // 43s_00
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*ray*raz
	  *(3*rax*rax-ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(70.0)*rax*ray*raz/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(70.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.25*sqrt(70.0)*ray*(3*rax*rax-ray*ray)/Rnorm5*damp[5]; // q=raz
	
      }
      else if (t==23 && u==0) { // 44c_00
	dTdq(tu,0) = -5/16.0*sqrt(35.0)*(pow(rax,4.0) - 6*rax*rax*ray*ray 
				       + pow(ray,4.0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(35.0)*(pow(rax,3.0) - 3*rax*ray*ray)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(35.0)*(pow(ray,3.0) - 3*rax*rax*ray)/Rnorm5*damp[5]; // q=ray
	
      }
      else if (t==24 && u==0) { // 44s_00
	dTdq(tu,0) = -1.25*sqrt(35.0)*rax*ray*(rax*rax-ray*ray)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*sqrt(35.0)*(3*rax*rax*ray - pow(ray,3.0))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*sqrt(35.0)*(pow(rax,3.0) - 3*ray*ray*rax)/Rnorm5*damp[5]; // q=ray
      }


      else if (t==0 && u==16) { // 00_40
	dTdq(tu,0) = -0.3125*(35*pow(rbz,4.0) - 30*rbz*rbz + 3)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,6) = 2.5*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==17) { // 00_41c
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rbx*pow(rbz,3.0) 
				      - 3*rbx*rbz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(10.0)*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.75*sqrt(10.0)*(7.0*rbx*rbz*rbz - rbx)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==18) { // 00_41s
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7.0*rby*pow(rbz,3.0) 
				      - 3*rby*rbz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(10.0)*(7.0*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.75*sqrt(10.0)*(7*rby*rbz*rbz - rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==19) { // 00_42c
	dTdq(tu,0) = -5/8.0*sqrt(5.0)*(7*rbz*rbz - 1.0)
	  *(rbx*rbx - rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rbx/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -0.5*sqrt(5.0)*(7*rbz*rbz-1)*rby/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 3.5*sqrt(5.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==20) { // 00_42s
	dTdq(tu,0) = -1.25*sqrt(5.0)*(7*rbz*rbz-1)*rbx*rby/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rby/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(5.0)*(7*rbz*rbz-1)*rbx/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 7*sqrt(5.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rby
      }
      else if (t==0 && u==21) { // 00_43c
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rbx*rbz
	  *(rbx*rbx-3*rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(70.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(70.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.25*sqrt(70.0)*rbx*(rbx*rbx-3*rby*rby)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==0 && u==22) { // 00_43s
	dTdq(tu,0) = -5/8.0*sqrt(70.0)*rby*rbz
	  *(3*rbx*rbx-rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(70.0)*rbx*rby*rbz/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.75*sqrt(70.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.25*sqrt(70.0)*rby*(3*rbx*rbx-rby*rby)/Rnorm5*damp[5]; // q=rbz
	
      }
      else if (t==0 && u==23) { // 00_44c
	dTdq(tu,0) = -5/16.0*sqrt(35.0)*(pow(rbx,4.0) - 6*rbx*rbx*rby*rby 
				       + pow(rby,4.0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(35.0)*(pow(rbx,3.0) - 3*rbx*rby*rby)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(35.0)*(pow(rby,3.0) - 3*rbx*rbx*rby)/Rnorm5*damp[5]; // q=rby
	
      }
      else if (t==0 && u==24) { // 00_44s
	dTdq(tu,0) = -1.25*sqrt(35.0)*rbx*rby*(rbx*rbx-rby*rby)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.5*sqrt(35.0)*(3*rbx*rbx*rby - pow(rby,3.0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(35.0)*(pow(rbx,3.0) - 3*rby*rby*rbx)/Rnorm5*damp[5]; // q=rby
      }

      // Dipole-octupole terms - Ltot = 4
      else if (t==9 && u>=1 && u<=3) { // 30_1*
	dTdq(tu,0) = -1.25*(35*pow(raz,3.0)*rB[u-1] 
			    + 15*raz*raz*cAB(2,u-1)
			    - 15*raz*rB[u-1] - 3*cAB(2,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 7.5*(7*raz*raz*rB[u-1] + 2*raz*cAB(2,u-1) 
			  - rB[u-1])/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 2.5*(7*pow(raz,3.0) - 3*raz)/Rnorm5*damp[5]; // q=rB*
      } 
      else if (t==10 && u>=1 && u<=3) { // 31c_1*
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rax*raz*raz*rB[u-1]
				     + 5*raz*raz*cAB(0,u-1)
				     + 10*rax*raz*cAB(2,u-1)
				     - 5*rax*rB[u-1] 
				     - cAB(0,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*sqrt(6.0)*(35*raz*raz*rB[u-1] + 10*raz*cAB(2,u-1) 
				   - 5*rB[u-1])/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.5*sqrt(6.0)*(35*rax*raz*rB[u-1] + 5*raz*cAB(0,u-1)
				  + 5*rax*cAB(2,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 0.25*sqrt(6.0)*(35*rax*raz*raz - 5*rax)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==11 && u>=1 && u<=3) { // 31s_1*
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*ray*raz*raz*rB[u-1]
				     + 5*raz*raz*cAB(1,u-1)
				     + 10*ray*raz*cAB(2,u-1)
				     - 5*ray*rB[u-1] 
				     - cAB(1,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.25*sqrt(6.0)*(35*raz*raz*rB[u-1] + 10*raz*cAB(2,u-1) 
				   - 5*rB[u-1])/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(6.0)*(35*ray*raz*rB[u-1] + 5*raz*cAB(1,u-1)
				  + 5*ray*cAB(2,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 0.25*sqrt(6.0)*(35*ray*raz*raz - 5*ray)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==12 && u>=1 && u<=3) { // 32c_1*
	dTdq(tu,0) = -1.25*sqrt(15.0)*((rax*rax-ray*ray)
				     *(7*raz*rB[u-1] + cAB(2,u-1)) 
				     + 2*raz*(rax*cAB(0,u-1) 
					      - ray*cAB(1,u-1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*(rax*(7*raz*rB[u-1] + cAB(2,u-1))
			       + raz*cAB(0,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -sqrt(15.0)*(ray*(7*raz*rB[u-1] + cAB(2,u-1))
				+ raz*cAB(1,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*sqrt(15.0)*(7*(rax*rax-ray*ray)*rB[u-1] 
				   + 2*(rax*cAB(0,u-1) 
					- ray*cAB(1,u-1)))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 3.5*sqrt(15.0)*raz*(rax*rax-ray*ray)/Rnorm5*damp[5]; //q=rb*
	
      }
      else if (t==13 && u>=1 && u<=3) { // 32s_1*
	dTdq(tu,0) = -2.5*sqrt(15.0)*(rax*ray*(7*raz*rB[u-1] + cAB(2,u-1))
				    + raz*(rax*cAB(1,u-1) 
					   + ray*cAB(0,u-1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = sqrt(15.0)*(ray*(7*raz*rB[u-1] + cAB(2,u-1)) 
			       + raz*cAB(1,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = sqrt(15.0)*(rax*(7*raz*rB[u-1] + cAB(2,u-1)) 
			       + raz*cAB(0,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = sqrt(15.0)*(7*rax*ray*rB[u-1] + rax*cAB(1,u-1) 
			       + ray*cAB(0,u-1))/Rnorm5*damp[5]; // q=raz
	dTdq(tu,u+3) = 7*sqrt(15.0)*rax*ray*raz/Rnorm5*damp[5]; //q=rb*
      }
      else if (t==14 && u>=1 && u<=3) { // 33c_1*
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7*pow(rax,3.0)*rB[u-1] 
				      + 3*(rax*rax-ray*ray)*cAB(0,u-1)
				      - 21*rax*ray*ray*rB[u-1]
				      - 6*rax*ray*cAB(1,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.75*sqrt(10.0)*(7*rax*rax*rB[u-1] 
				    + 2*rax*cAB(0,u-1) 
				    - 7*ray*ray*rB[u-1] 
				    - 2*ray*cAB(1,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = -1.5*sqrt(10.0)*(ray*cAB(0,u-1) + 7*rax*ray*rB[u-1] 
				    + rax*cAB(1,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,u+3) = 1.75*sqrt(10.0)*(pow(rax,3.0) 
				      - 3*rax*ray*ray)/Rnorm5*damp[5]; // q=rB*
      }
      else if (t==15 && u>=1 && u<=3) { // 33s_1*
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(-7*pow(ray,3.0)*rB[u-1] 
				      + 3*(rax*rax-ray*ray)*cAB(1,u-1)
				      + 21*rax*rax*ray*rB[u-1]
				      + 6*rax*ray*cAB(0,u-1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 1.5*sqrt(10.0)*(rax*cAB(1,u-1) 
				   + 7*rax*ray*rB[u-1] 
				   + ray*cAB(0,u-1))/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.75*sqrt(10.0)*(-7*ray*ray*rB[u-1] 
				    - 2*ray*cAB(1,u-1)
				    + 7*rax*rax*rB[u-1]
				    + 2*rax*cAB(0,u-1))/Rnorm5*damp[5]; // q=ray
	dTdq(tu,u+3) = 1.75*sqrt(10.0)*(3*rax*rax*ray 
				      - pow(rax,3.0))/Rnorm5*damp[5]; // q=rb*
      }

      else if (u==9 && t>=1 && t<=3) { // 1*_30
	dTdq(tu,0) = -1.25*(35*pow(rbz,3.0)*rA[t-1] 
			    + 15*rbz*rbz*cAB(t-1,2)
			    - 15*rbz*rA[t-1] - 3*cAB(t-1,2))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,6) = 7.5*(7*rbz*rbz*rA[t-1] + 2*rbz*cAB(t-1,2) 
			  - rA[t-1])/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 2.5*(7*pow(rbz,3.0) - 3*rbz)/Rnorm5*damp[5]; // q=ra*
      } 
      else if (u==10 && t>=1 && t<=3) { // 1*_31c
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rbx*rbz*rbz*rA[t-1]
				     + 5*rbz*rbz*cAB(t-1,0)
				     + 10*rbx*rbz*cAB(t-1,2)
				     - 5*rbx*rA[t-1] 
				     - cAB(t-1,0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.25*sqrt(6.0)*(35*rbz*rbz*rA[t-1] + 10*rbz*cAB(t-1,2) 
				   - 5*rA[t-1])/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.5*sqrt(6.0)*(35*rbx*rbz*rA[t-1] + 5*rbz*cAB(t-1,0)
				  + 5*rbx*cAB(t-1,2))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 0.25*sqrt(6.0)*(35*rbx*rbz*rbz - 5*rbx)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==11 && t>=1 && t<=3) { // 1*_31s
	dTdq(tu,0) = -5/8.0*sqrt(6.0)*(35*rby*rbz*rbz*rA[t-1]
				     + 5*rbz*rbz*cAB(t-1,1)
				     + 10*rby*rbz*cAB(t-1,2)
				     - 5*rby*rA[t-1] 
				     - cAB(t-1,1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,5) = 0.25*sqrt(6.0)*(35*rbz*rbz*rA[t-1] + 10*rbz*cAB(t-1,2) 
				   - 5*rA[t-1])/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(6.0)*(35*rby*rbz*rA[t-1] + 5*rbz*cAB(t-1,1)
				  + 5*rby*cAB(t-1,2))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 0.25*sqrt(6.0)*(35*rby*rbz*rbz - 5*rby)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==12 && t>=1 && t<=3) { // 1*_32c
	dTdq(tu,0) = -1.25*sqrt(15.0)*((rbx*rbx-rby*rby)
				     *(7*rbz*rA[t-1] + cAB(t-1,2)) 
				     + 2*rbz*(rbx*cAB(t-1,0) 
					      - rby*cAB(t-1,1)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*(rbx*(7*rbz*rA[t-1] + cAB(t-1,2))
			       + rbz*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -sqrt(15.0)*(rby*(7*rbz*rA[t-1] + cAB(t-1,2))
				+ rbz*cAB(t-1,1))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(15.0)*(7*(rbx*rbx-rby*rby)*rA[t-1] 
				   + 2*(rbx*cAB(t-1,0) 
					- rby*cAB(t-1,1)))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 3.5*sqrt(15.0)*rbz*(rbx*rbx-rby*rby)/Rnorm5*damp[5]; //q=ra*
	
      }
      else if (u==13 && t>=1 && t<=3) { // 1*_32s
	dTdq(tu,0) = -2.5*sqrt(15.0)*(rbx*rby*(7*rbz*rA[t-1] + cAB(t-1,2))
				    + rbz*(rbx*cAB(t-1,1) 
					   + rby*cAB(t-1,0)))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = sqrt(15.0)*(rby*(7*rbz*rA[t-1] + cAB(t-1,2)) 
			       + rbz*cAB(t-1,1))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = sqrt(15.0)*(rbx*(7*rbz*rA[t-1] + cAB(t-1,2)) 
			       + rbz*cAB(t-1,0))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = sqrt(15.0)*(7*rbx*rby*rA[t-1] + rbx*cAB(t-1,1) 
			       + rby*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbz
	dTdq(tu,t) = 7*sqrt(15.0)*rbx*rby*rbz/Rnorm5*damp[5]; //q=ra*
      }
      else if (u==14 && t>=1 && t<=3) { // 1*_33c
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(7*pow(rbx,3.0)*rA[t-1] 
				      + 3*(rbx*rbx-rby*rby)*cAB(t-1,0)
				      - 21*rbx*rby*rby*rA[t-1]
				      - 6*rbx*rby*cAB(t-1,1))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 0.75*sqrt(10.0)*(7*rbx*rbx*rA[t-1] 
				    + 2*rbx*cAB(t-1,0) 
				    - 7*rby*rby*rA[t-1] 
				    - 2*rby*cAB(t-1,1))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -1.5*sqrt(10.0)*(rby*cAB(t-1,0) + 7*rbx*rby*rA[t-1] 
				    + rbx*cAB(1,t-1))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,t) = 1.75*sqrt(10.0)*(pow(rbx,3.0) 
				    - 3*rbx*rby*rby)/Rnorm5*damp[5]; // q=ra*
      }
      else if (u==15 && t>=1 && t<=3) { // 1*_33s
	dTdq(tu,0) = -5/8.0*sqrt(10.0)*(-7*pow(rby,3.0)*rA[t-1] 
				      + 3*(rbx*rbx-rby*rby)*cAB(t-1,1)
				      + 21*rbx*rbx*rby*rA[t-1]
				      + 6*rbx*rby*cAB(t-1,0))/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,4) = 1.5*sqrt(10.0)*(rbx*cAB(t-1,1) 
				   + 7*rbx*rby*rA[t-1] 
				   + rby*cAB(t-1,0))/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.75*sqrt(10.0)*(-7*rby*rby*rA[t-1] 
				    - 2*rby*cAB(t-1,1)
				    + 7*rbx*rbx*rA[t-1]
				    + 2*rbx*cAB(t-1,0))/Rnorm5*damp[5]; // q=rby
	dTdq(tu,t) = 1.75*sqrt(10.0)*(3*rbx*rbx*rby 
				    - pow(rby,3.0))/Rnorm5*damp[5]; // q=ra*
      }



      // Quadrupole-quadrupole terms - Ltot = 4
      else if (t==4 && u==4) { // 20_20
	dTdq(tu,0) = -15/8.0*(35*raz*raz*rbz*rbz - 5*raz*raz 
			      - 5*rbz*rbz + 20*raz*rbz*czz 
			      + 2*czz*czz + 1)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.75*(70*raz*rbz*rbz - 10*raz 
			   + 20*rbz*czz)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,6) = 0.75*(70*raz*raz*rbz - 10*rbz 
			   + 20*raz*czz)/Rnorm5*damp[5]; // q=rbz
      }

      else if (t==4 && u==5) { // 20_21c
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rbx*rbz - 5*rbx*rbz 
				    + 10*raz*rbx*czz + 10*raz*rbz*czx
				    + 2*czx*czz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rbx*rbz + 10*rbx*czz 
				  + 10*rbz*czx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*sqrt(3.0)*(35*raz*raz*rbz - 5*rbz 
				  + 10*raz*czz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = 0.5*sqrt(3.0)*(35*raz*raz*rbx - 5*rbx 
				  + 10*raz*czx)/Rnorm5*damp[5]; // q=rbz

	// 21c_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*rax*raz - 5*rax*raz 
				    + 10*rbz*rax*czz + 10*rbz*raz*cxz
				    + 2*cxz*czz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*rax*raz + 10*rax*czz 
				  + 10*raz*cxz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*sqrt(3.0)*(35*rbz*rbz*raz - 5*raz 
				  + 10*rbz*czz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,3) = 0.5*sqrt(3.0)*(35*rbz*rbz*rax - 5*rax 
				  + 10*rbz*cxz)/Rnorm5*damp[5]; // q=raz
      }

      else if (t==4 && u==6) { // 20_21s
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rby*rbz - 5*rby*rbz 
				    + 10*raz*rby*czz + 10*raz*rbz*czy
				    + 2*czy*czz)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rby*rbz + 10*rby*czz 
				  + 10*rbz*czy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = 0.5*sqrt(3.0)*(35*raz*raz*rbz - 5*rbz 
				  + 10*raz*czz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = 0.5*sqrt(3.0)*(35*raz*raz*rby - 5*rby
				  + 10*raz*czy)/Rnorm5*damp[5]; // q=rbz

	// 21s_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*ray*raz - 5*ray*raz 
				    + 10*rbz*ray*czz + 10*rbz*raz*cyz
				    + 2*cyz*czz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*ray*raz + 10*ray*czz 
				  + 10*raz*cyz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,2) = 0.5*sqrt(3.0)*(35*rbz*rbz*raz - 5*raz 
				  + 10*rbz*czz)/Rnorm5*damp[5]; // q=ray
	dTdq(ut,3) = 0.5*sqrt(3.0)*(35*rbz*rbz*ray - 5*ray
				  + 10*rbz*cyz)/Rnorm5*damp[5]; // q=raz
      }
      else if (t==4 && u==7) { // 20_22c
	dTdq(tu,0) = -0.625*sqrt(3.0)*(35*raz*raz*rbx*rbx 
				     - 35*raz*raz*rby*rby 
				     - 5*rbx*rbx + 5*rby*rby 
				     + 20*raz*rbx*czx
				     - 20*raz*rby*czy + 2*czx*czx
				     - 2*czy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.25*sqrt(3.0)*(70*raz*rbx*rbx 
				   - 70*raz*rby*rby
				   + 20*rbx*czx 
				   - 20*rby*czy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.25*sqrt(3.0)*(70*raz*raz*rbx - 10*rbx 
				   + 20*raz*czx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.25*sqrt(3.0)*(-70*raz*raz*rby + 10*rby
				   - 20*raz*czy)/Rnorm5*damp[5]; // q=rby
	
	// 22c_20
	dTdq(ut,0) = -0.625*sqrt(3.0)*(35*rbz*rbz*rax*rax 
				     - 35*rbz*rbz*ray*ray 
				     - 5*rax*rax + 5*ray*ray 
				     + 20*rbz*rax*cxz
				     - 20*rbz*ray*cyz + 2*cxz*cxz
				     - 2*cyz*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.25*sqrt(3.0)*(70*rbz*rax*rax 
				   - 70*rbz*ray*ray
				   + 20*rax*cxz 
				   - 20*ray*cyz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.25*sqrt(3.0)*(70*rbz*rbz*rax - 10*rax 
				   + 20*rbz*cxz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.25*sqrt(3.0)*(-70*rbz*rbz*ray + 10*ray
				   - 20*rbz*cyz)/Rnorm5*damp[5]; // q=ray

      }
      else if (t==4 && u==8) { // 20_22s
	dTdq(tu,0) = -1.25*sqrt(3.0)*(35*raz*raz*rbx*rby 
				    - 5*rbx*rby
				    + 10*raz*rbx*czy
				    + 10*raz*rby*czx 
				    + 2*czx*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,3) = 0.5*sqrt(3.0)*(70*raz*rbx*rby + 10*rbx*czy 
				  + 10*rby*czx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*sqrt(3.0)*(35*raz*raz*rby - 5*rby 
				  + 10*raz*czy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*sqrt(3.0)*(35*raz*raz*rbx - 5*rbx 
				  + 10*raz*czx)/Rnorm5*damp[5]; // q=rby

	// 22s_20
	dTdq(ut,0) = -1.25*sqrt(3.0)*(35*rbz*rbz*rax*ray 
				    - 5*rax*ray
				    + 10*rbz*rax*cyz
				    + 10*rbz*ray*cxz 
				    + 2*cxz*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,6) = 0.5*sqrt(3.0)*(70*rbz*rax*ray + 10*rax*cyz 
				  + 10*ray*cxz)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*sqrt(3.0)*(35*rbz*rbz*ray - 5*ray 
				  + 10*rbz*cyz)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*sqrt(3.0)*(35*rbz*rbz*rax - 5*rax 
				  + 10*rbz*cxz)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==5 && u==5) { // 21c_21c
	dTdq(tu,0) = -2.5*(35*rax*raz*rbx*rbz + 5*rax*rbx*czz 
			   + 5*rax*rbz*czx + 5*raz*rbx*cxz 
			   + 5*raz*rbz*cxx + cxx*czz 
			   + cxz*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rbx*rbz + 5*rbx*czz + 5*rbz*czx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rbx*rbz + 5*rbx*cxz + 5*rbz*cxx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*rax*raz*rbz + 5*rax*czz + 5*raz*cxz)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,6) = (35*rax*raz*rbx + 5*rax*czx + 5*raz*cxx)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==5 && u==6) { // 21c_21s
	dTdq(tu,0) = -2.5*(35*rax*raz*rby*rbz + 5*rax*rby*czz 
			   + 5*rax*rbz*czy + 5*raz*rby*cxz 
			   + 5*raz*rbz*cxy + cxy*czz 
			   + cxz*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rby*rbz + 5*rby*czz + 5*rbz*czy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rby*rbz + 5*rby*cxz + 5*rbz*cxy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = (35*rax*raz*rbz + 5*rax*czz + 5*raz*cxz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = (35*rax*raz*rby + 5*rax*czy + 5*raz*cxy)/Rnorm5*damp[5]; // q=rbz

	// 21s_21c
	dTdq(ut,0) = -2.5*(35*rbx*rbz*ray*raz + 5*rbx*ray*czz 
			   + 5*rbx*raz*cyz + 5*rbz*ray*czx 
			   + 5*rbz*raz*cyx + cyx*czz 
			   + czx*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = (35*rbz*ray*raz + 5*ray*czz + 5*raz*cyz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = (35*rbx*ray*raz + 5*ray*czx + 5*raz*cyx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,2) = (35*rbx*rbz*raz + 5*rbx*czz + 5*rbz*czx)/Rnorm5*damp[5]; // q=ray
	dTdq(ut,3) = (35*rbx*rbz*ray + 5*rbx*cyz + 5*rbz*cyx)/Rnorm5*damp[5]; // q=raz
      }
      
      else if (t==5 && u==7) { // 21c_22c
	dTdq(tu,0) = -1.25*(35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby
			    + 10*rax*rbx*czx - 10*rax*rby*czy 
			    + 10*raz*rbx*cxx - 10*raz*rby*cxy
			    + 2*cxx*czx - 2*cxy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*(35*raz*rbx*rbx - 35*raz*rby*rby 
			  + 10*rbx*czx - 10*rby*czy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = 0.5*(35*rax*rbx*rbx - 35*rax*rby*rby
			  + 10*rbx*cxx - 10*rby*cxy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*(70*rax*raz*rbx + 10*rax*czx 
			  + 10*raz*cxx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*(-70*rax*raz*rby - 10*rax*czy 
			  - 10*raz*cxy)/Rnorm5*damp[5]; // q=rby
	
	// 22c_21c     
	dTdq(ut,0) = -1.25*(35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray
			    + 10*rbx*rax*cxz - 10*rbx*ray*cyz 
			    + 10*rbz*rax*cxx - 10*rbz*ray*cyx
			    + 2*cxx*cxz - 2*cyx*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = 0.5*(35*rbz*rax*rax - 35*rbz*ray*ray 
			  + 10*rax*cxz - 10*ray*cyz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = 0.5*(35*rbx*rax*rax - 35*rbx*ray*ray
			  + 10*rax*cxx - 10*ray*cyx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*(70*rbx*rbz*rax + 10*rbx*cxz 
			  + 10*rbz*cxx)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*(-70*rbx*rbz*ray - 10*rbx*cyz 
			  - 10*rbz*cyx)/Rnorm5*damp[5]; // q=ray
	
      }
      else if (t==5 && u==8) { // 21c_22s
	dTdq(tu,0) = -2.5*(35*rax*raz*rbx*rby + 5*rax*rbx*czy
			   + 5*rax*rby*czx + 5*raz*rbx*cxy
			   + 5*raz*rby*cxx + cxx*czy 
			   + cxy*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*raz*rbx*rby + 5*rbx*czy + 5*rby*czx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,3) = (35*rax*rbx*rby + 5*rbx*cxy + 5*rby*cxx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*rax*raz*rby + 5*rax*czy + 5*raz*cxy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*rax*raz*rbx + 5*rax*czx + 5*raz*cxx)/Rnorm5*damp[5]; // q=rby
	
	// 22s_21c     
	dTdq(ut,0) = -2.5*(35*rbx*rbz*rax*ray + 5*rbx*rax*cyz
			   + 5*rbx*ray*cxz + 5*rbz*rax*cyx
			   + 5*rbz*ray*cxx + cxx*cyz 
			   + cyx*cxz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = (35*rbz*rax*ray + 5*rax*cyz + 5*ray*cxz)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,6) = (35*rbx*rax*ray + 5*rax*cyx + 5*ray*cxx)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = (35*rbx*rbz*ray + 5*rbx*cyz + 5*rbz*cyx)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = (35*rbx*rbz*rax + 5*rbx*cxz + 5*rbz*cxx)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==6 && u==6) { // 21s_21s
	dTdq(tu,0) = -2.5*(35*ray*raz*rby*rbz + 5*ray*rby*czz
			   + 5*ray*rbz*czy + 5*raz*rby*cyz
			   + 5*raz*rbz*cyy + cyy*czz 
			   + cyz*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = (35*raz*rby*rbz + 5*rby*czz + 5*rbz*czy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = (35*ray*rby*rbz + 5*rby*cyz + 5*rbz*cyy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,5) = (35*ray*raz*rbz + 5*ray*czz + 5*raz*cyz)/Rnorm5*damp[5]; // q=rby
	dTdq(tu,6) = (35*ray*raz*rby + 5*ray*czy + 5*raz*cyy)/Rnorm5*damp[5]; // q=rbz
      }
      else if (t==6 && u==7) { // 21s_22c
	dTdq(tu,0) = -1.25*(35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby
			    + 10*ray*rbx*cxz - 10*ray*rby*czy
			    + 10*raz*rbx*cyx - 10*raz*rby*cyy +
			    2*cyx*czx - 2*cyy*czy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = 0.5*(35*raz*rbx*rbx - 35*raz*rby*rby
			  + 10*rbx*czx - 10*rby*czy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = 0.5*(35*ray*rbx*rbx - 35*ray*rby*rby
			  + 10*rbx*cyx - 10*rby*cyy)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = 0.5*(70*ray*raz*rbx + 10*ray*czx 
			  + 10*raz*cyx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = -0.5*(70*ray*raz*rby + 10*ray*czy
			   + 10*raz*cyy)/Rnorm5*damp[5]; // q=rby

	// 22c_21s
	dTdq(ut,0) = -1.25*(35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray
			    + 10*rby*rax*czx - 10*rby*ray*cyz
			    + 10*rbz*rax*cxy - 10*rbz*ray*cyy +
			    2*cxy*cxz - 2*cyy*cyz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,5) = 0.5*(35*rbz*rax*rax - 35*rbz*ray*ray
			  + 10*rax*cxz - 10*ray*cyz)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,6) = 0.5*(35*rby*rax*rax - 35*rby*ray*ray
			  + 10*rax*cxy - 10*ray*cyy)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = 0.5*(70*rby*rbz*rax + 10*rby*cxz 
			  + 10*rbz*cxy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = -0.5*(70*rby*rbz*ray + 10*rby*cyz
			   + 10*rbz*cyy)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==6 && u==8) { // 21s_22s
	dTdq(tu,0) = -2.5*(35*ray*raz*rbx*rby + 5*ray*rbx*czy
			   + 5*ray*rby*czx + 5*raz*rbx*cyy
			   + 5*raz*rby*cyx + cyx*czy 
			   + cyy*czx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,2) = (35*raz*rbx*rby + 5*rbx*czy + 5*rby*czx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,3) = (35*ray*rbx*rby + 5*rbx*cyy + 5*rby*cyx)/Rnorm5*damp[5]; // q=raz
	dTdq(tu,4) = (35*ray*raz*rby + 5*ray*czy + 5*raz*cyy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*ray*raz*rbx + 5*ray*czx + 5*raz*cyx)/Rnorm5*damp[5]; // q=rby

	// 22s_21s
	dTdq(ut,0) = -2.5*(35*rby*rbz*rax*ray + 5*rby*rax*cyz
			   + 5*rby*ray*cxz + 5*rbz*rax*cyy
			   + 5*rbz*ray*cxy + cxy*cyz 
			   + cyy*cxz)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,5) = (35*rbz*rax*ray + 5*rax*cyz + 5*ray*cxz)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,6) = (35*rby*rax*ray + 5*rax*cyy + 5*ray*cxy)/Rnorm5*damp[5]; // q=rbz
	dTdq(ut,1) = (35*rby*rbz*ray + 5*rby*cyz + 5*rbz*cyy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = (35*rby*rbz*rax + 5*rby*cxz + 5*rbz*cxy)/Rnorm5*damp[5]; // q=ray

      }
      else if (t==7 && u==7) { // 22c_22c
	dTdq(tu,0) = -0.625*(35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby
			     - 35*ray*ray*rbx*rbx
			     + 35*ray*ray*rby*rby
			     + 20*rax*rbx*cxx - 20*rax*rby*cxy
			     - 20*ray*rbx*cyx + 20*ray*rby*cyy
			     + 2*cxx*cxx - 2*cxy*cxy 
			     - 2*cyx*cyx + 2*cyy*cyy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.25*(70*rax*rbx*rbx - 70*rax*rby*rby
			   + 20*rbx*cxx - 20*rby*cxy)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.25*(-70*ray*rbx*rbx + 70*ray*rby*rby
			   - 20*rbx*cyx + 20*rby*cyy)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = 0.25*(70*rax*rax*rbx - 70*ray*ray*rbx
			   + 20*rax*cxx - 20*ray*cyx)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.25*(-70*rax*rax*rby + 70*ray*ray*rby
			   - 20*rax*cxy + 20*ray*cyy)/Rnorm5*damp[5]; // q=rby
      }
      else if (t==7 && u==8) { // 22c_22s
	dTdq(tu,0) = -1.25*(35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby
			    + 10*rax*rbx*cxy + 10*rax*rby*cxx
			    - 10*ray*rbx*cyy - 10*ray*rby*cyx
			    + 2*cxx*cxy - 2*cyx*cyy)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = 0.5*(70*rax*rbx*rby + 10*rbx*cxy 
			  + 10*rby*cxx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = 0.5*(-70*ray*rbx*rby - 10*rbx*cyy 
			  - 10*rby*cyx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = 0.5*(35*rax*rax*rby - 35*ray*ray*rby
			  + 10*rax*cxy - 10*ray*cyy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = 0.5*(35*rax*rax*rbx - 35*ray*ray*rbx
			  + 10*rax*cxx - 10*ray*cyx)/Rnorm5*damp[5]; // q=rby

	// 22s_22c
	dTdq(ut,0) = -1.25*(35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray
			    + 10*rbx*rax*cyx + 10*rbx*ray*cxx
			    - 10*rby*rax*cyy - 10*rby*ray*cxy
			    + 2*cxx*cyx - 2*cxy*cyy)/Rnorm7*damp[5] 
	  + Tab(u,t)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(ut,4) = 0.5*(70*rbx*rax*ray + 10*rax*cyx 
			  + 10*ray*cxx)/Rnorm5*damp[5]; // q=rbx
	dTdq(ut,5) = 0.5*(-70*rby*rax*ray - 10*rax*cyy 
			  - 10*ray*cxy)/Rnorm5*damp[5]; // q=rby
	dTdq(ut,1) = 0.5*(35*rbx*rbx*ray - 35*rby*rby*ray
			  + 10*rbx*cyx -10*rby*cyy)/Rnorm5*damp[5]; // q=rax
	dTdq(ut,2) = 0.5*(35*rbx*rbx*rax - 35*rby*rby*rax
			  + 10*rbx*cxx - 10*rby*cxy)/Rnorm5*damp[5]; // q=ray
      }
      else if (t==8 && u==8) { // 22s_22s
	dTdq(tu,0) = -2.5*(35*rax*ray*rbx*rby + 5*rax*rbx*cyy
			   + 5*rax*rby*cyx + 5*ray*rbx*cxy
			   + 5*ray*rby*cxx + cxx*cyy 
			   + cxy*cyx)/Rnorm7*damp[5] 
	  + Tab(t,u)/(2.0*Rnorm)*damp_grad[5]; // q=R*R
	dTdq(tu,1) = (35*ray*rbx*rby + 5*rbx*cyy + 5*rby*cyx)/Rnorm5*damp[5]; // q=rax
	dTdq(tu,2) = (35*rax*rbx*rby + 5*rbx*cxy + 5*rby*cxx)/Rnorm5*damp[5]; // q=ray
	dTdq(tu,4) = (35*rax*ray*rby + 5*rax*cyy + 5*ray*cxy)/Rnorm5*damp[5]; // q=rbx
	dTdq(tu,5) = (35*rax*ray*rbx + 5*rax*cyx + 5*ray*cxx)/Rnorm5*damp[5]; // q=rby
      }



    }
  }

  // Build dT/dX = dT/dq * dq/dX;
  dTdX = dTdq.Multiply(dqdX);
  dTdX.Scale(1.0/perm);// Divide by 4*pi*epsilon
  
  //dTdX.Print("\ndT/dX");

  delete [] eA;
  delete [] eB;

  return dTdX;
}

// Compute the Tang-Toennies damping factor and its nuclear derivatives
double Atom::TangToenniesDampingFactor(int n, double beta, double R, 
				       int ideriv) {
  // The Tang-Toennies damping factor for dispersion/induction
  // F_n = 1 - exp(-beta*R) sum [(beta*R)^k / k!]
  //
  // n is the order of the term being damped, R is the distance
  // between the two interacting sites, and beta is an empirical
  // parameter.  
  //
  // By default, ideriv = 0, meaning we just want the damping factor.
  // If ideriv = 1, we compute dF_n/dR.  Higher derivatives are not
  // implemented at present.

  double damp;
  double kfactorial;
  
  if (ideriv == 0) {
    damp = 1.0;
    for (int k=0; k<=n; k++) {
      // Evaluate k!
      if (k==0) 
	kfactorial = 1.0;
      else
	kfactorial *= k;
      
      // Evaluate the contribution to the damping factor
      damp -= exp(-beta*R)*pow(beta*R,k)/kfactorial;
    }
  }
  else if (ideriv == 1) {
    for (int k=0; k<=n; k++) {
      // Evaluate k!
      if (k==0) 
	kfactorial = 1.0;
      else
	kfactorial *= k;
    }
    damp = beta*exp(-beta*R)*pow(beta*R,n)/kfactorial;
  }
  else {
    printf("Atom::TangToenniesDampingFactor() - Unknown derivative order: %d\n",ideriv);
    exit(1);
  }
  
  return damp;
}


void Atom::PrintTinkerCartesian(int shift, FILE *outfile) {
  // shift if it is not the first fragment in a list


  // Print everything up to the connectivity
  fprintf(outfile,"%2d  %2s  %10.6f  %10.6f  %10.6f  %d",atom_index+shift, 
	  AtSym.c_str(), xyz[0], xyz[1], xyz[2], MM_atom_type);


  // Print the connectivity
  for (int i=0;i<Nconnected;i++) {
    fprintf(outfile,"  %d",connectivity[i]+shift);
  }

  // End the line
  fprintf(outfile,"\n");

}

//Print coordinates for atom that differ for its current coordates in Tinker format
void Atom::PrintTinkerCartesian(double x, double y, double z, int shift, FILE *outfile) {

  // Print everything up to the connectivity
    // Print everything up to the connectivity
  fprintf(outfile,"%2d  %2s  %10.6f  %10.6f  %10.6f  %d",atom_index+shift, 
	  AtSym.c_str(), x, y, z, MM_atom_type);

  // Print the connectivity
  for (int i=0;i<Nconnected;i++) {
    fprintf(outfile,"  %d",connectivity[i]+shift);
  }

  // End the line
  fprintf(outfile,"\n");

}

void Atom::PrintTinkerEmbeddingCharge(int shift, FILE *outfile) {
  // shift if it is not the first fragment in a list
  
  // Modify atomic symbol to indicate charge. Solely for sake of clarity.
  string AtomName = "Q" + AtSym; 

  // Print everything up to the connectivity
  // note: shift MM_atom_type by 300 for embedding charges
  fprintf(outfile,"%2d  %2s  %10.6f  %10.6f  %10.6f  %d",atom_index+shift, 
	  AtomName.c_str(), xyz[0], xyz[1], xyz[2], MM_atom_type+300);

  // Print the connectivity
  for (int i=0;i<Nconnected;i++) {
    fprintf(outfile,"  %d",connectivity[i]+shift);
  }

  // End the line
  fprintf(outfile,"\n");


}

// Overload '=' operator
Atom& Atom::operator=(const Atom& other) {

  if (this!=&other) {
    AtSym = other.AtSym;
    AtNum = other.AtNum;
    mass = other.mass;

    // AIFF parameters
    DispersionAtomType_ = other.DispersionAtomType_;
    C6 = other.C6;
    C8 = other.C8;
    C10 = other.C10;
    C6_UCHF = other.C6_UCHF;
    C8_UCHF = other.C8_UCHF;
    C10_UCHF = other.C10_UCHF;
    isotropic_dipole_polarizability = other.isotropic_dipole_polarizability;

    // Atomic coordinates
    xyz.Initialize(3);
    xyz = other.xyz;

    local_xyz.Initialize(3);
    local_xyz = other.local_xyz;

    if(Params::Parameters().IsPeriodic() || Params::Parameters().UseFullQMOnly()){
      fractional_xyz.Initialize(3);
      fractional_xyz = other.fractional_xyz;
    }

    freq_pol_dipole.Initialize(10); // by shuhao
    freq_pol_dipole = other.freq_pol_dipole;// by shuhao
    
    freq_pol_quad.Initialize(10); // by shuhao
    freq_pol_quad = other.freq_pol_quad;// by shuhao

    point_charge = other.point_charge;


    atom_index = other.atom_index;
    MM_atom_type = other.MM_atom_type;
    Nconnected = other.Nconnected;
    
    //symmetry by yoni
    global_index = other.global_index;
    Sym_Atom = other.Sym_Atom;
    Rotation = other.Rotation;

    // connectivity
    if ( Nconnected > 0 ) {
      connectivity = new int[Nconnected];
      for (int i=0;i<Nconnected;i++) {
	connectivity[i] = other.connectivity[i];
      }
    }
    else {
      connectivity = NULL;
      //connectivity = new int[2];
      //connectivity[0] = 0;
    }

    MultipoleMoments = other.MultipoleMoments;
    Pols = other.Pols;
    init_multipoles = other.init_multipoles;
    init_pols = other.init_pols;
    init_InduceMultipoles = other.init_InduceMultipoles;

    // JDH magnetic properties stuff...
    Monomer3x3Tensor = other.Monomer3x3Tensor;
    TwoBody3x3Tensor = other.TwoBody3x3Tensor;
    Cluster3x3Tensor = other.Cluster3x3Tensor;
    EwaldPotential = other.EwaldPotential;

    EmbeddedCharge = other.EmbeddedCharge;
    EmbeddedChargeTwoBody = other.EmbeddedChargeTwoBody;
    
    //NMRMonomerShieldingTensor = other.NMRMonomerShieldingTensor;
    //NMRTwoBodyShieldingTensor = other.NMRTwoBodyShieldingTensor;
    //NMRClusterShieldingTensor = other.NMRClusterShieldingTensor;

  }
  return *this;
}

double Atom::LookupAtomicDispersionParameter(string property) {

  string *names = new string [26];
  Vector Pol(26), tabC6(26), tabC9(26), Rvdw(26);
  // note: called these "tabCX" for tabulated, to distinguish from AIFF CX coeffs

  // Parameters are from von Lilienfeld and Tkatchenko,
  // J. Chem. Phys. 132, 234109 (2010).
  names[0] = "Hfree";   Pol[0]=4.50;  tabC6[0]=6.5;   tabC9[0]=21.6;  Rvdw[0]=3.10; 
  names[1] = "Hs";      Pol[1]=2.75;  tabC6[1]=2.42;  tabC9[1]=4.91;  Rvdw[1]=2.63; 
  names[2] = "He";      Pol[2]=1.38;  tabC6[2]=1.46;  tabC9[2]=1.47;  Rvdw[2]=2.65; 
  names[3] = "Cfree";   Pol[3]=12.0;  tabC6[3]=46.6;  tabC9[3]=373.;  Rvdw[3]=3.59; 
  names[4] = "Csp";     Pol[4]=9.73;  tabC6[4]=30.6;  tabC9[4]=199.;  Rvdw[4]=3.35; 
  names[5] = "Csp2";    Pol[5]=9.67;  tabC6[5]=30.3;  tabC9[5]=195.;  Rvdw[5]=3.34; 
  names[6] = "Csp3";    Pol[6]=8.64;  tabC6[6]=24.1;  tabC9[6]=139.;  Rvdw[6]=3.22; 
  names[7] = "Nfree";   Pol[7]=7.40;  tabC6[7]=24.2;  tabC9[7]=117.;  Rvdw[7]=3.34; 
  names[8] = "Nsp2sp3"; Pol[8]=6.36;  tabC6[8]=17.9;  tabC9[8]=74.4;  Rvdw[8]=3.18; 
  names[9] = "Ofree";   Pol[9]=5.40;  tabC6[9]=15.6;  tabC9[9]=52.6;  Rvdw[9]=3.19; 
  names[10] = "Osp2";   Pol[10]=4.92; tabC6[10]=13.0; tabC9[10]=39.8; Rvdw[10]=3.09; 
  names[11] = "Osp3";   Pol[11]=4.81; tabC6[11]=12.4; tabC9[11]=37.1; Rvdw[11]=3.07; 
  names[12] = "Ffree";  Pol[12]=3.80; tabC6[12]=9.52; tabC9[12]=24.2; Rvdw[12]=3.04; 
  names[13] = "Fsp3";   Pol[13]=3.46; tabC6[13]=7.89; tabC9[13]=18.3; Rvdw[13]=2.95; 
  names[14] = "Ne";     Pol[14]=2.67; tabC6[14]=6.38; tabC9[14]=12.0; Rvdw[14]=2.91; 
  names[15] = "Sifree"; Pol[15]=37.0; tabC6[15]=305.; tabC9[15]=8550; Rvdw[15]=4.20; 
  names[16] = "Sisp3";  Pol[16]=25.6; tabC6[16]=146.; tabC9[16]=2846; Rvdw[16]=3.72;
  names[17] = "Pfree";  Pol[17]=25.0; tabC6[17]=185.; tabC9[17]=3561; Rvdw[17]=4.01; 
  names[18] = "Sfree";  Pol[18]=19.6; tabC6[18]=134.; tabC9[18]=1925; Rvdw[18]=3.86; 
  names[19] = "Ssp3";   Pol[19]=18.2; tabC6[19]=115.; tabC9[19]=1532; Rvdw[19]=3.76; 
  names[20] = "Clfree"; Pol[20]=15.0; tabC6[20]=94.6; tabC9[20]=1014; Rvdw[20]=3.71; 
  names[21] = "Clsp3";  Pol[21]=14.6; tabC6[21]=89.4; tabC9[21]=932.; Rvdw[21]=3.68; 
  names[22] = "Ar";     Pol[22]=11.1; tabC6[22]=64.3; tabC9[22]=518.; Rvdw[22]=3.55; 
  names[23] = "Brfree"; Pol[23]=20.0; tabC6[23]=162.; tabC9[23]=2511; Rvdw[23]=3.93; 
  names[24] = "Brsp3";  Pol[24]=19.5; tabC6[24]=155.; tabC9[24]=2340; Rvdw[24]=3.90; 
  names[25] = "Kr";     Pol[25]=16.8; tabC6[25]=130.; tabC9[25]=1572; Rvdw[25]=3.82; 


  int match_index = -1;
  bool match = false;
  for (int i=0;i<26;i++) {
    if (names[i] == DispersionAtomType_) {
      match_index = i;
      match = true;
      break;
    }
  }
  if (! match) {
      printf("Atom::LookupDispersionParameters(): Atom Type |%s| not found\n",DispersionAtomType_.c_str());
      exit(1);
    }


  delete [] names;
  if (property == "Pol") return Pol[match_index];
  else if (property == "C6") return tabC6[match_index];
  else if (property == "C9") return tabC9[match_index];
  else if (property == "Rvdw") return Rvdw[match_index];
  else {
    printf("Atom::LookupDispersionParameters(): Unknown property %s\n",property.c_str());
    exit(1);
  }

}


double Atom::EstimateC6Coefficient(Atom other) {

  // Grab polarizabilities
  double polA = LookupAtomicDispersionParameter("Pol");
  double polB = other.LookupAtomicDispersionParameter("Pol");

  // Grab pure atomic C6 coefficients
  double C6aa = LookupAtomicDispersionParameter("C6");
  double C6bb = other.LookupAtomicDispersionParameter("C6");

  double Pa = C6aa*polB/polA;
  double Pb = C6bb*polA/polB;

  double C6ab = 2*C6aa*C6bb / ( Pa + Pb );
  //printf("C6ab = %f\n",C6ab);
  return C6ab;
}

double Atom::EstimateC9Coefficient(Atom other1, Atom other2, string type) {

  double C9ijk = 0.0;

  if (type == "Tkatchenko") {
    // Grab polarizabilities
    double polI = LookupAtomicDispersionParameter("Pol");
    double polJ = other1.LookupAtomicDispersionParameter("Pol"); 
    double polK = other2.LookupAtomicDispersionParameter("Pol"); 
    
    // Grab pure atomic C9 coefficients
    double C9iii = LookupAtomicDispersionParameter("C9");
    double C9jjj = other1.LookupAtomicDispersionParameter("C9");
    double C9kkk = other2.LookupAtomicDispersionParameter("C9");
    
    // Compute the C9 coefficient for this triplet ijk
    double Pi = C9iii * polJ * polK / (polI*polI);
    double Pj = C9jjj * polI * polK / (polJ*polJ);
    double Pk = C9kkk * polI * polJ / (polK*polK);
    
    C9ijk = (8.0/3.0)*Pi*Pj*Pk*(Pi+Pj+Pk)/((Pi+Pj)*(Pj+Pk)*(Pk+Pi));
    //printf("new C9ijk = %f\n",C9ijk);
  }
  else { // default
    // Use the Tang (Phys Rev 177, 108-114, 1969) rule for estimating
    // C9 from the three individual C6 coefficients and the isotropic
    // dipole polarizabilities.  This is also described in Stone's
    // Theory of IMF book, Chapter 9.

    // Grab the isotropic dipole polarizabilities
    double polI = GetIsotropicDipolePolarizability();
    double polJ = other1.GetIsotropicDipolePolarizability();
    double polK = other2.GetIsotropicDipolePolarizability();
    
    // Grab the individual atom C6 coefficients.  These are typically
    // computed via Casimir-Polder integration of the
    // frequency-dependent polarizabilities.
    double C6ii = GetC6Coefficient();
    double C6jj = other1.GetC6Coefficient();
    double C6kk = other2.GetC6Coefficient();    

    // Compute a set of useful intermediates.
    double Si = C6ii*polJ*polK/polI;
    double Sj = C6jj*polI*polK/polJ;
    double Sk = C6kk*polI*polJ/polK;
    
    // Evaluate C9
    C9ijk = 2*Si*Sj*Sk*(Si+Sj+Sk)/((Si+Sj)*(Sj+Sk)*(Sk+Si));

  }


  return C9ijk;
}

Vector Atom::GaussLegendreWeights(double x1, double x2, int n) {
  // Adapted from Numerical Recipes routine gauleg, Chapter 4.5
  // It computes the n abcissas x values and their corresponding weights w 
  // on which the quadrature should be performed.  The input values x1 & x2
  // are the limits of integration.  Once you have w and x computed, the
  // integral is approximated as:
  //
  //   int_x1^x2 = sum_i^n f(x_i)*w(i)
  //
  // This routine also factors in the extra terms to the weights due to the
  // change of variables that maps the frequency omega to t according to:
  //
  //                omega = omega0 * (1+x)/(1-x)
  //
  // In particular, the integration variable maps d(omega) = 2*omega0/(1-x)^2 dx
  //
  // This change of variables converts the semi-infinite
  // Casimir-Polder integral (i.e. from 0 to infinity) to one between
  // -1 and 1.  See Misquitta and Stone, Mol. Phys. 106, 1631-1643
  // (2008).

  double omega0 = 0.5; // base frequency in the change of variables.
		       // See the Misquitta and Stone paper mentioned
		       // above.

  Vector Weights(n);

  if (n==10) {
    // Use precomputed weights if n=10.  In practice, this is probably the only
    // part of this function we use.  But the rest is here, just in case.
    Weights[0] = 0.017111419760829;
    Weights[1] = 0.042964786325738;
    Weights[2] = 0.077678726763695;
    Weights[3] = 0.131054117334672;
    Weights[4] = 0.223896873021545;
    Weights[5] = 0.407948854240782;
    Weights[6] = 0.838730580639330;
    Weights[7] = 2.131641821331806;
    Weights[8] = 8.208052005690481;
    Weights[9] = 97.920920814877775;
  }
  else {
    // Determine the weights on the fly
    // Notation: 
    // * x here are the values of the frequencies at which the
    //   polarizabilities are evaluated, albeit in a transformed 
    //   coordinate.
    // * w are the raw Gauss-Legendre weights for each value of x.
    // * tm1sq = (1-x)^2  
    Vector x(n+1); 
    Vector w(n+1); // the raw Gauss-Legendre weights
    Vector tm1sq(n+1); // (1-x)^2.  Used for the change of variables.
    
    Vector tmpWeights(n+1);

    double EPS = 3.0e-11;
    double p1, p2, p3, pp;
    
    int m = (n+1)/2;
    double xm = 0.5*(x2+x1);
    double xl = 0.5*(x2-x1);
    
    for (int i=1;i<=m;i++) {
      double z = cos(pi*(i-0.25)/(n+0.5));
      double z1 = 100.0; // nonsense initial value to ensure we enter the while loop.
      
      while ( fabs(z-z1) > EPS ) {
	p1 = 1.0;
	p2 = 0.0;
	
	for (int j=1;j<=n;j++) {
	  p3 = p2;
	  p2 = p1;
	  p1 = ((2.0*j-1.0)*z*p2 - (j-1.0)*p3)/j;
        }
	
	pp = n*(z*p1-p2)/(z*z-1.0);
	z1 = z;
	z = z1-p1/pp; 
      }//end while
      
      x[i] = xm-xl*z;
      x[n+1-i] = xm+xl*z;
      
      tm1sq[i] = (1.0 - x[i])*(1.0 - x[i]);
      tm1sq[n+1-i] = (1.0 - x[n+1-i])*(1.0 - x[n+1-i]);
      
      w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
      w[n+1-i] = w[i];
      
      // Scale weights with extra terms due to the change of variables from
      // omega to x.
      tmpWeights[i]=w[i]*2.0*omega0/tm1sq[i];
      tmpWeights[n+1-i]=w[n+1-i]*2.0*omega0/tm1sq[n+1-i];
    } 
    
    // The Numerical Recipes routine uses fortran-style indexing
    // (i.e. arrays start at 1).  So now we shift to C-style (start at
    // 0).
    for (int q=0;q<n;q++) {
      Weights[q] = tmpWeights[q+1];
    }
  }
  //Weights.Print("Gauss-Legendre Weights\n");
  return Weights; 
}

// Compute the C6 dispersion coefficient by Casimir-Polder integration
double Atom::CasimirC6Coefficient(Atom other) {	
  double C6 = 0.0;
  double CPint = 0.0;
  int Nquad = 10;
  Vector weights = GaussLegendreWeights(-1.0, 1.0, Nquad);

  // These are computed via Casimir-Polder integration of the
  // (uncoupled) isotropic frequency-dependent polarizabilities using
  // Gauss-Legendre quadrature.

  // Grab the frequency-dependent isotropic polarizabilities
  Vector alphaI = GetFreq_Pol_Dipole();
  Vector alphaJ = other.GetFreq_Pol_Dipole();
  //alphaI.Print("alphaI\n");
  //alphaJ.Print("alphaJ\n");
  

  for (int i=0;i<Nquad;i++) {
    CPint = CPint + weights[i]*alphaI[i]*alphaJ[i];
  }    
  // Evaluate C6, the values are from the Casimir-Polder integration   
  C6 = 6.0*CPint/(2.0*pi);
  //printf("C6ij=%f\n",C6);
  

  return C6;	
}

// Compute the C8 dispersion coefficient by Casimir-Polder integration
double Atom::CasimirC8Coefficient(Atom other) {	
  double C8 = 0.0;
  double CPint = 0.0;
  int Nquad = 10;
  Vector weights = GaussLegendreWeights(-1.0, 1.0, Nquad);
  
  // These are computed via Casimir-Polder integration of the
  // (uncoupled) isotropic frequency-dependent polarizabilities using
  // Gauss-Legendre quadrature.
  
  // Grab the frequency-dependent isotropic polarizabilities
  Vector alphaI = GetFreq_Pol_Dipole();
  Vector alphaIQ = GetFreq_Pol_Quad(); 
  
  Vector alphaJ = other.GetFreq_Pol_Dipole();
  Vector alphaJQ = other.GetFreq_Pol_Quad();
  
  for (int i=0;i<Nquad;i++) {
    CPint = CPint + weights[i]*alphaIQ[i]*alphaJ[i]+weights[i]*alphaI[i]*alphaJQ[i];
  }    
  C8 = 15.0*CPint/(2.0*pi);  
  //printf("C8ij=%f\n",C8);

  return C8;	
}

// Compute the C10 dispersion coefficient by Casimir-Polder
// integration This routine only captures the piece that comes from
// quadrupole-quadrupole.  The dipole-octopole terms are probably
// similar in magnitude, but they are neglected here since we usually
// don't have dipole-octopole polarizabilities available.  Overall,
// it's probably best to skip the C10 terms.
double Atom::CasimirC10Coefficient(Atom other) {	
  double C10 = 0.0;
  double CPint = 0.0;
  int Nquad = 10;
  Vector weights = GaussLegendreWeights(-1.0, 1.0, Nquad);
  
  // These are computed via Casimir-Polder integration of the
  // (uncoupled) isotropic frequency-dependent polarizabilities using
  // Gauss-Legendre quadrature.
  
  // Grab the frequency-dependent isotropic polarizabilities
  Vector alphaIQ = GetFreq_Pol_Quad();
  Vector alphaJQ = other.GetFreq_Pol_Quad();
  
  for (int i=0;i<Nquad;i++) {
    CPint = CPint + weights[i]*alphaIQ[i]*alphaJQ[i];
  }    
  // Evaluate C9, the values are from the Casimir-Polder integration   
  C10 = 70.0*CPint/(2.0*pi);
  //printf("C10ij=%f\n",C10);

  return C10;	
}

// Compute the 3-body Axilrod-Teller C9 dispersion coefficient
double Atom::CasimirC9Coefficient(Atom other1, Atom other2) {	
  double C9 = 0.0;
  double CPint = 0.0;
  int Nquad = 10;
  Vector weights = GaussLegendreWeights(-1.0, 1.0, Nquad);
  
  // These are computed via Casimir-Polder integration of the
  // (uncoupled) isotropic frequency-dependent polarizabilities using
  // Gauss-Legendre quadrature.
  
  // Grab the frequency-dependent isotropic polarizabilities
  Vector alphaI = GetFreq_Pol_Dipole();
  Vector alphaJ = other1.GetFreq_Pol_Dipole();
  Vector alphaK = other2.GetFreq_Pol_Dipole();
  
  for (int i=0;i<Nquad;i++) {
    CPint = CPint + weights[i]*alphaI[i]*alphaJ[i]*alphaK[i];
  }    
  // Evaluate C9, the values are from the Casimir-Polder integration   
  C9 = 6.0*CPint/(2.0*pi);
     
  return C9;	
}


void Atom::SetMonomer3x3Tensor(Matrix tensor3x3) {
  if (tensor3x3.GetRows() != 3 || tensor3x3.GetCols() != 3 ) {
    printf("ERROR: Atom::SetMonomer3x3Tensor requires a 3x3 tensor\n");
    exit(1);
  }

  Monomer3x3Tensor = tensor3x3;

}

void Atom::SetTwoBody3x3Tensor(Matrix tensor3x3) {
  if (tensor3x3.GetRows() != 3 || tensor3x3.GetCols() != 3 ) {
    printf("ERROR: Atom::SetTwoBody3x3Tensor requires a 3x3 tensor\n");
    exit(1);
  }
  TwoBody3x3Tensor = tensor3x3;
}

void Atom::SetCluster3x3Tensor(Matrix tensor3x3) {
  if (tensor3x3.GetRows() != 3 || tensor3x3.GetCols() != 3 ) {
    printf("ERROR: Atom::SetCluster3x3Tensor requires a 3x3 tensor\n");
    exit(1);
  }
  Cluster3x3Tensor = tensor3x3;
}

void Atom::PrintEmbeddingCharges(FILE *outfile) {
  // CHARGE:
  if ( Params::Parameters().UseElectrostaticEmbedding() ) {
    if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
      fprintf(outfile, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", xyz[0],xyz[1],
	      xyz[2], GetMultipoleMoments().GetMoments().Element(0)  );
    } else {
      fprintf(outfile, "%10.6f  %10.6f  %10.6f  %10.6f\n", xyz[0],xyz[1],
	    xyz[2], GetMultipoleMoments().GetMoments().Element(0) ); 
    }
  }
}

double Atom::GetVanDerWaalsRadius() {
  
  if ( GetSymbol() == "H" ) {
    van_der_waals_radius = 1.2;
  } else if ( GetSymbol() == "C" ) {
    van_der_waals_radius = 1.7;
  } else if ( GetSymbol() == "N" ) {
    van_der_waals_radius = 1.55;
  } else if ( GetSymbol() == "O" ) {
    van_der_waals_radius = 1.52;
  } else if ( GetSymbol() == "F" ) {
    van_der_waals_radius = 1.47; 
  } else if ( GetSymbol() == "Cl" ) {
    van_der_waals_radius = 1.75;
  } else if ( GetSymbol() == "Br" ) {
    van_der_waals_radius = 1.85;
  } else if ( GetSymbol() == "P" ) {
    van_der_waals_radius = 1.8;
  } else if ( GetSymbol() == "S" ) {
    van_der_waals_radius = 1.8;
  } else if ( GetSymbol() == "Sn" ) {
    van_der_waals_radius = 2.17;
  } else if ( GetSymbol() == "Na" ) {
    van_der_waals_radius = 2.27;
  } else if ( GetSymbol() == "K" ) {
    van_der_waals_radius = 2.75;
  } else if ( GetSymbol() == "He" ) {
    van_der_waals_radius = 1.4;
  } else {
    printf("Atom::GetVanDerWaalsRadius()  Radius not tabulated\n");
    exit(1);
  }


  return van_der_waals_radius;

}
