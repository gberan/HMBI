#include <string>
using std::string;
#include "multipole.h"

Multipole::Multipole(): Rank(-1), Nmom(0) {

}

// primary constructor
Multipole::Multipole(int rank, double* moments) : Rank(rank), Nmom(0) {

  // Figure out how many spherical tensor moments this rank corresponds to.
  Nmom = 0;
  for (int i=0;i<=rank;i++) 
    Nmom += 2*i + 1;

  // Note, dipole moments are read as 10 (z), 11c (x), 11s (y).  We will store them
  // as x,y,z, or 11c, 11s, 10, so we need to reorder them.
  double z = moments[1];
  double x = moments[2];
  double y = moments[3];
  moments[1] = x;
  moments[2] = y;
  moments[3] = z;

  // Initialize the moments vector
  Vector tmp(moments, Nmom);
  

  Moments.Initialize(Nmom);
  Moments = tmp;

}

// Copy constructor
Multipole::Multipole(const Multipole& other, bool copy_data) : Rank(other.Rank), Nmom(other.Nmom) {

  // Read the ranks
  //Rank = other.GetRank();
  //Nmom = other.GetLength();

  // Initialize the vectors
  Moments.Initialize(Nmom);

  if (copy_data) 
    Moments = other.GetMoments();
  else 
    Moments.Set();

}

Multipole::~Multipole() {

}

void Multipole::Initialize(int rank) {

  Rank = rank;
  // Figure out how many spherical tensor moments this rank corresponds to.
  Nmom = 0;
  for (int i=0;i<=rank;i++) 
    Nmom += 2*i + 1;
  
  Moments.Initialize(Nmom);

}

void Multipole::Initialize(const Multipole& other, bool copy_data) {

   // Read the ranks
  Rank = other.GetRank();
  Nmom = other.GetLength();

  // Initialize the vectors
  Moments.Initialize(Nmom);

  if (copy_data) 
    Moments = other.Moments;
  else 
    Moments.Set();

}


// Multiply the multipole moments by a double
Multipole& Multipole::Scale(double factor) {
  Moments.Scale(factor);
  return *this;
}

// Add to multipole vectors of equivalent rank
Multipole Multipole::operator+(const Multipole& other) {

  // First test that the ranks of the two multipole vectors are
  // identical
  if (other.Rank != Rank) {
    printf("ERROR: Multipole addition requires both multipoles have the same rank\n");
    exit(1);
  }

  // Create an empty multipole object with the proper rank
  Multipole result;
  result.Initialize(Rank);

  result.Moments = Moments;
  result.Moments += other.Moments;

  return result;

}



// Overload the '=' operator
Multipole& Multipole::operator=(const Multipole& other) {
  if (this!=&other) {
    // Read the ranks
    Rank = other.Rank;
    Nmom = other.Nmom;
    
    // Initialize the moments
    if (Nmom > 0) {
      Moments = other.Moments;
    }
  }
  return *this;
}

// Print out the multipole moments with pretty formatting
void Multipole::Print(string title) {  

  if (Rank < 0) {
    printf("Error: Multipole::Print(): Unknown Rank = %d\n",Rank); 
    exit(1);
  }

  printf("%s:  Rank %d\n",title.c_str(),Rank);
  printf("  Charge:       (00)  %10.6f\n",Moments[0]);
  if (Rank >=1) {
    printf("  Dipole:       (10)  %10.6f  (11c) %10.6f  (11s) %10.6f\n",
	   Moments[3],Moments[1],Moments[2]);
  }
  if (Rank >=2) {
    printf("  Quadrupole:   (20)  %10.6f  (21c) %10.6f  (21s) %10.6f  (22c) %10.6f  (22s) %10.6f\n",
	   Moments[4],Moments[5],Moments[6],Moments[7],Moments[8]);
  }
  if (Rank >=3) {
    printf("  Octopole:     (30)  %10.6f  (31c) %10.6f  (31s) %10.6f  (32c) %10.6f  (32s) %10.6f\n",
	   Moments[9],Moments[10],Moments[11],Moments[12],Moments[13]);
    printf("                (33c) %10.6f  (33s) %10.6f\n",
	   Moments[14],Moments[15]);
  }
  if (Rank >=4) {
    printf("  Hexadecapole: (40)  %10.6f  (41c) %10.6f  (41s) %10.6f  (42c) %10.6f  (42s) %10.6f\n",
	   Moments[16],Moments[17],Moments[18],Moments[19],Moments[20]);
    printf("                (43c) %10.6f  (43s) %10.6f  (44c) %10.6f  (44s) %10.6f\n",
	   Moments[21],Moments[22],Moments[23],Moments[24]);
  }

  printf("\n");
}

// JDH use in Dalton multipole emb. 
Vector Multipole::Spherical_to_Cartesian(Multipole Spherical){
//converts from spherical moments to cartesian moments
  Vector S_mom = Spherical.GetMomentsVector();
  int dim,rank = Spherical.GetRank();
  if(rank==0){//dim of cartesian moments vector
    dim=1;
  }
  else if(rank==1){
    dim=4;
  }
  else if(rank==2){
    dim=10;
  }
  else if(rank==3){
    dim=20;
  }
  else if(rank==4){
    dim=35;
  }
    double rt2 = 1.4142135623730950488, rt3 = 1.7320508075688772935, rt8=2.8284271247461900976, rt5 = 2.2360679774997896964, rt7 = 2.6457513110645905905;
    double rt10=rt2*rt5, rt12=rt2*rt3*rt2, rt16=rt8*rt2,rt24=rt12*rt2, rt35=rt5*rt7, rt42=rt2*rt3*rt7, rt70=rt2*rt5*rt7;
//precompute some roots because they are handy
//note the moments are stored as xx xy xz yy yz zz or xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz for use in dalton
//Not the weird order Stone uses!
  Vector C_mom(dim);
  //charge and dipole are identical between spherical and cartesian

  if(rank>=0)
    C_mom[0]=S_mom[0];
  if(rank>=1){
    C_mom[1]=S_mom[1];
    C_mom[2]=S_mom[2];
    C_mom[3]=S_mom[3];
  }
  if(rank>=2){
    C_mom[4]=-(1.0/2.0)*S_mom[4]+(1.0/2.0)*rt3*S_mom[7];//xx
    C_mom[5]=(1.0/2.0)*rt3*S_mom[8];//xy
    C_mom[6]=(1.0/2.0)*rt3*S_mom[5];//xz
    C_mom[7]=-(1.0/2.0)*S_mom[4]-(1.0/2.0)*rt3*S_mom[7];//yy
    C_mom[8]=(1.0/2.0)*rt3*S_mom[6];//yz
    C_mom[9]=S_mom[4];//zz
  }
  if(rank>=3){
    C_mom[10]=(rt5/rt8)*S_mom[14]-(rt3/rt8)*S_mom[10];//xxx
    C_mom[11]=(rt5/rt8)*S_mom[15]-(1.0/rt24)*S_mom[11];//xxy
    C_mom[12]=(rt5/rt12)*S_mom[12]-(1.0/2.0)*S_mom[9];//xxz
    C_mom[13]=(rt5/rt8)*S_mom[14]-(1.0/rt24)*S_mom[10];//xyy
    C_mom[14]=(rt5/rt12)*S_mom[13];//xyz
    C_mom[15]=(rt2/rt3)*S_mom[10];//xzz
    C_mom[16]=(rt5/rt8)*S_mom[15]-(rt3/rt8)*S_mom[11];//yyy
    C_mom[17]=-(rt5/rt12)*S_mom[12]-(1.0/2.0)*S_mom[9];//yyz
    C_mom[18]=(rt2/rt3)*S_mom[11];//yzz
    C_mom[19]=S_mom[9];//zzz
  }
  if(rank>=4){
    C_mom[20]=(3.0/8.0)*S_mom[16]-(1.0/4.0)*rt5*S_mom[19]+(1.0/32.0)*rt35*S_mom[23];//xxxx
    C_mom[21]=(1.0/8.0)*(-rt5*S_mom[20]+rt35*S_mom[24]);//xxxy
    C_mom[22]=(1.0/16.0)*(-rt10*S_mom[17]+rt70*S_mom[21]);//xxxz
    C_mom[23]=(1.0/8.0)*S_mom[16]-(1.0/32.0)*rt35*S_mom[23];//xxyy
    C_mom[24]=(1.0/16.0)*(-rt10*S_mom[18]+rt70*S_mom[22]);//xxyz
    C_mom[25]=-(1.0/2.0)*S_mom[16]+(1.0/4.0)*rt5*S_mom[19];//xxzz
    C_mom[26]=-(1.0/8.0)*(rt5*S_mom[20]+rt35*S_mom[24]);//xyyy
    C_mom[27]=-(1.0/16.0)*(rt10*S_mom[17]+rt70*S_mom[21]);//xyyz
    C_mom[28]=(1.0/4.0)*rt5*S_mom[20];//xyzz
    C_mom[29]=(5.0/8.0)*S_mom[17];//xzzz
    C_mom[30]=(3.0/8.0)*S_mom[16]+(1.0/4.0)*rt5*S_mom[19]+(1.0/32.0)*rt35*S_mom[23];//yyyy
    C_mom[31]=-(1.0/16.0)*(3.0*rt10*S_mom[18]+rt70*S_mom[22]);//yyyz
    C_mom[32]=-(1.0/2.0)*S_mom[16]-(1.0/4.0)*rt5*S_mom[19];//yyzz
    C_mom[33]=(rt5/rt8)*S_mom[18];//yzzz
    C_mom[34]=S_mom[16];//zzzz
  }
 
  return C_mom;
}
