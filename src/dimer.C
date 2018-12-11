#include "dimer.h"
#include "quasiharmonic.h"
#include "constants.h"
#include <sys/stat.h>
using namespace hmbi_constants;

Dimer::Dimer() : Tabs(NULL), DampedTabs(NULL), TabsGrad(NULL), 
		 DampedTabsGrad(NULL) {
 
  MM_Grad_Init = 0;
  QM_Grad_Init = 0;
  MM_Hess_Init = 0;
  QM_Hess_Init = 0;
  Grad_Electrostatic_Init = 0;
  Tabs_Init = 0;
  DampedTabs_Init = 0;
  TabsGrad_Init = 0;
  DampedTabsGrad_Init = 0;

}

// Overload "=" operator
// Assumes none of the arrays have been initialized in the copy
Dimer& Dimer::operator=(const Dimer& other) {
  if (this!=&other) {
    // copy the monomers forming the dimer
    MonA = other.MonA;
    MonB = other.MonB;

    //copy the indices of these two monomers
    indexA = other.indexA;
    indexB = other.indexB;
    
    // if dimer image copy the ref_mon index
    reference_MonB = other.reference_MonB;

    // copy the types of monomers
    typeA = other.typeA;
    typeB = other.typeB;

    //copy the boolean for images
    is_image = other.is_image;

    // copy the K_vector[3]
    K_vector[0] = other.K_vector[0];
    K_vector[1] = other.K_vector[1];
    K_vector[2] = other.K_vector[2];

    //copy the sym info 
    sym_fac = other.sym_fac;
    sym_fac_period = other.sym_fac_period;
    List_Sym_Dim = other.List_Sym_Dim;
    List_Sym_Periodic = other.List_Sym_Periodic;
    //MonB_List = other.MonB_List;
    k_list = other.k_list;
    Rotation = other.Rotation;
    Rotation_Periodic = other.Rotation_Periodic;
    Atom_Equivalency = other.Atom_Equivalency;
    Atom_Equivalency_Periodic = other.Atom_Equivalency_Periodic;

    Energy_QM = other.Energy_QM;
    Energy_MM = other.Energy_MM;
    dEint_QM = other.dEint_QM;
    dEint_MM = other.dEint_MM;
    E_Electrostatic_MM = other.E_Electrostatic_MM;
    E_Induction_MM = other.E_Induction_MM;
    E_2b_disp_MM = other.E_2b_disp_MM;
    OrientEnergy = other.OrientEnergy;

    Separation = other.Separation;
    minA = other.minA;
    minB = other.minB;

    spin = other.spin;
    charge = other.charge;
    Natoms = other.Natoms;

    //none of the following have been tested when initialized
    QM_Grad_Init = other.QM_Grad_Init;
//    if(Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().IsSupercellJob())
//       (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()))
    //if(Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
    //   (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()) ||
    //   (Params::Parameters().DoFreq() && Params::Parameters().IsSupercellJob()) )
    if(QM_Grad_Init) 
       Grad_QM = other.Grad_QM;;

    MM_Grad_Init = other.MM_Grad_Init;
//    if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().IsSupercellJob())
//	 (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()))
    //if(Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
    //   (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()) ||
    //   (Params::Parameters().DoFreq() && Params::Parameters().IsSupercellJob()) )
    if(MM_Grad_Init)
      Grad_MM = other.Grad_MM;

    Grad_Electrostatic_Init = other.Grad_Electrostatic_Init;
    if(Grad_Electrostatic_Init)
      Grad_Electrostatic = other.Grad_Electrostatic;

    Tabs_Init = other.Tabs_Init;
    if(other.Tabs_Init){
      int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
      Tabs = new Matrix[MultiplyAtoms];
      for(int i=0;i<MultiplyAtoms;i++){
	Tabs[i] = other.Tabs[i];
      }
    }

    DampedTabs_Init = other.DampedTabs_Init;
    if(other.DampedTabs_Init){
      int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();		
      DampedTabs = new Matrix[MultiplyAtoms];
      for(int i=0;i<MultiplyAtoms;i++){
	DampedTabs[i] = other.DampedTabs[i];
      }
    }

    TabsGrad_Init = other.TabsGrad_Init;
    if(other.TabsGrad_Init){
      int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
      for(int i=0;i<MultiplyAtoms;i++){
	TabsGrad[i] = other.TabsGrad[i];
      }
    }

    DampedTabsGrad_Init = other.DampedTabsGrad_Init;
    if(other.DampedTabsGrad_Init){
      int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
      for(int i=0;i<MultiplyAtoms;i++){
	DampedTabsGrad[i] = other.DampedTabsGrad[i];
      }
    }

    QM_Hess_Init = other.QM_Hess_Init;
    //if( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() )){
    if(QM_Hess_Init){
      int N = 3*(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  Hess_QM(i,j) = other.Hess_QM(i,j);
    }

    MM_Hess_Init = other.MM_Hess_Init;
    //if( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() )){
    if(MM_Hess_Init){
      int N = 3*(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  Hess_MM(i,j) = other.Hess_MM(i,j);
    }

  }
  return *this;
}

Dimer::Dimer(const Dimer &other) :     spin(other.spin),
                                         charge(other.charge),
                                         Natoms(other.Natoms),
                                         Energy_QM(other.Energy_QM),
                                         Energy_MM(other.Energy_MM),
                                         QM_Grad_Init(other.QM_Grad_Init),
                                         MM_Grad_Init(other.MM_Grad_Init),
                                         QM_Hess_Init(other.QM_Hess_Init),
                                         MM_Hess_Init(other.MM_Hess_Init)  
{
  // copy the monomers forming the dimer
  MonA = other.MonA;
  MonB = other.MonB;
  
  //copy the indices of these two monomers
  indexA = other.indexA;
  indexB = other.indexB;

  // if dimer image copy the ref_mon index
  reference_MonB = other.reference_MonB; 

  // copy the types of monomers
  typeA = other.typeA;
  typeB = other.typeB;

  //copy the boolean for images
  is_image = other.is_image;

  // copy the K_vector[3];
  K_vector[0] = other.K_vector[0];
  K_vector[1] = other.K_vector[1];
  K_vector[2] = other.K_vector[2];

  //copy the symmetry info
  sym_fac = other.sym_fac;
  sym_fac_period = other.sym_fac_period;
  List_Sym_Dim = other.List_Sym_Dim;
  List_Sym_Periodic = other.List_Sym_Periodic;
  //MonB_List = other.MonB_List;
  k_list = other.k_list;
  Rotation = other.Rotation;
  Rotation_Periodic = other.Rotation_Periodic;
  Atom_Equivalency = other.Atom_Equivalency;
  Atom_Equivalency_Periodic = other.Atom_Equivalency_Periodic;


  Energy_QM = other.Energy_QM;
  Energy_MM = other.Energy_MM;
  dEint_QM = other.dEint_QM;
  dEint_MM = other.dEint_MM;
  E_Electrostatic_MM = other.E_Electrostatic_MM;
  E_Induction_MM = other.E_Induction_MM;
  E_2b_disp_MM = other.E_2b_disp_MM;
  OrientEnergy = other.OrientEnergy; 

  minA = other.minA;
  minB = other.minB;
  Separation = other.Separation;

  // Copy over Gradients
  QM_Grad_Init = other.QM_Grad_Init;   
  if(QM_Grad_Init){
// if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().DoFreq()) { //JLM
  //if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
   //    (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable())){
 // if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
   //    (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()) ||
     //  (Params::Parameters().DoFreq() && Params::Parameters().IsSupercellJob()) ){
    Grad_QM = other.Grad_QM;
  }

  MM_Grad_Init = other.MM_Grad_Init;
  if(MM_Grad_Init){
//if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().DoFreq()) { //JLM
 // if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
  //     (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable())){  
 // if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || 
   //    (Params::Parameters().DoFreq() && !Params::Parameters().DoQuasiHarmonic() && !Params::Parameters().AreQHFAvailable()) ||
     //  (Params::Parameters().DoFreq() && Params::Parameters().IsSupercellJob()) ){     
    Grad_MM = other.Grad_MM;   
  }

  Grad_Electrostatic_Init = other.Grad_Electrostatic_Init;
  if(Grad_Electrostatic_Init)
    Grad_Electrostatic = other.Grad_Electrostatic;
  
  Tabs_Init = other.Tabs_Init;
  if(other.Tabs_Init){
    int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
    Tabs = new Matrix[MultiplyAtoms];
    for(int i=0;i<MultiplyAtoms;i++){
      Tabs[i] = other.Tabs[i];
    }
  }
  
  DampedTabs_Init = other.DampedTabs_Init;
  if(other.DampedTabs_Init){
    int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();	   
    DampedTabs = new Matrix[MultiplyAtoms];
    for(int i=0;i<MultiplyAtoms;i++){
      DampedTabs[i] = other.DampedTabs[i];
    }
  }
  TabsGrad_Init = other.TabsGrad_Init;
  if(other.TabsGrad_Init){
    int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
    for(int i=0;i<MultiplyAtoms;i++){
      TabsGrad[i] = other.TabsGrad[i];
    }
  }
  DampedTabsGrad_Init = other.DampedTabsGrad_Init;
  if(other.DampedTabsGrad_Init){
    int MultiplyAtoms = MonA.GetNumberOfAtoms()*MonB.GetNumberOfAtoms();
    for(int i=0;i<MultiplyAtoms;i++){
      DampedTabsGrad[i] = other.DampedTabsGrad[i];
    }
  }
  

  // Copy over the Hessians
  QM_Hess_Init = other.QM_Hess_Init;
  if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
    Hess_QM = other.Hess_QM;
  }
  MM_Hess_Init = other.MM_Hess_Init;
  if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() )) {
    Hess_MM = other.Hess_MM;
  }


}

Dimer::~Dimer() {

  // Special clean up if using the AIFF
  if ( Params::Parameters().GetMMType()==2 ) { // AIFF
    delete [] Tabs;
    delete [] DampedTabs;

    if ( MM_Grad_Init ) {
      delete [] TabsGrad;
      delete [] DampedTabsGrad;
    }
  }
}

void Dimer::Initialize(Monomer& m1, Monomer& m2) {
  // Assign monomers such that index(A) < index(B)
  if ( m1.GetIndex() < m2.GetIndex() ) {
    MonA = m1; MonB = m2; 
  }
  else {
    MonA = m2; MonB = m1;
  }

  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex();

  typeA = MonA.GetType();
  typeB = MonB.GetType();

  //energy
  Energy_QM = 0;
  Energy_MM = 0;
  dEint_QM = 0;
  dEint_MM = 0;
  E_Electrostatic_MM = 0;
  E_Induction_MM = 0;
  E_2b_disp_MM = 0;

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();


  // Set spin state.  Crude algorithm for determining dimer spin state
  // from monomers:
  // 2 singlets -> singlet
  // singlet + other -> other
  // non-singlet + non-singlet -> singlet if even number of electrons
  //      e.g. 2 triplets: 3 + 3 = 6 -> singlet
  // non-singlet + non-singlet -> doublet if odd number of electrons
  //      e.g. doublet + triplet = 2 + 3 = 5 -> doublet
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else if (spinA>1 && spinB>1 && (spinA+spinB)%2==0)  
    spin = 1;
  else if (spinA>1 && spinB>1 && (spinA+spinB)%2==1)
    spin = 2; 
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  QM_Grad_Init = 0;
  MM_Grad_Init = 0;

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() || Params::Parameters().DoFreq()) {
    Grad_QM.Initialize(3*Natoms);
    Grad_MM.Initialize(3*Natoms);
  }

  MM_Hess_Init = 0;               
  QM_Hess_Init = 0;
  if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {
    Hess_QM.Initialize(3*Natoms, 3*Natoms);
    Hess_MM.Initialize(3*Natoms, 3*Natoms);
  }

  // PBC features - 0 by default
  is_image = false;
  reference_MonB = 0;
  K_vector[0] = 0;
  K_vector[1] = 0;
  K_vector[2] = 0;

  // No extra symmetry by default
  sym_fac = 1;
  sym_fac_period = 0;
  //Symmetry list should be empty but clearing it just in case
  List_Sym_Dim.clear();
  List_Sym_Periodic.clear();
  // setting Dimer Symmetrical to itself
  List_Sym_Dim.push_back(indexA);
  List_Sym_Dim.push_back(indexB);
  //printf("items in List_Sym_Dim ");
  //for(int i=0;i<List_Sym_Dim.size();i++)
  //  printf("%i ",List_Sym_Dim[i]);
  //printf("\n");

  //setting up the rotation matrix of all symmetrical dimers
  Rotation.clear();
  Rotation_Periodic.clear();
  Matrix Rot(3,3);
  Rot.Set_Iden();
  Rotation.push_back(Rot);

  //Atoms are equivalence to itself
  Atom_Equivalency.clear();
  Atom_Equivalency_Periodic.clear();
  Vector equivalence(Natoms);
  equivalence.Incremental();
  Atom_Equivalency.push_back(equivalence);

  //MonB_List.clear();
  k_list.clear();

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];
  //printf("Separation = %f (based on atoms %d and %d)\n",Separation,minA,minB);

  UseInTwoBody_ = false; // JDH default initialization behavior
  //  UseInTwoBodyCharge_ = false; // JDH default initialization behavior
  
}

// This function updates the Monomers, along with all members derived
// from those, in case those objects have changed.  Ideally would do
// this more dynamically using pointers...
void Dimer::UpdateObjects(Monomer& m1, Monomer& m2) {
  
  if ( Params::Parameters().GetJobTypeStr() == "nmr" )  {
    // Hack... preserve dimer NMR data here.  More reason to do the key
    // code bits with pointers, so we don't worry about overwriting data
    // we want to keep!  The problem here is that dimer NMR data is
    // stored on each individual monomer atom.  So it gets overwritten when
    // we update the monomers inside each dimer.
    printf("Updating NMR data for dimer (%d,%d)\n",indexA,indexB);
    /*
    Matrix *tmpNMR;
    tmpNMR = new Matrix[Natoms];

    int NatomsA = MonA.GetNumberOfAtoms();
    int NatomsB = MonB.GetNumberOfAtoms();
    for (int iatom = 0;iatom < NatomsA; iatom++) { //  MonA
      if ( MonA.GetAtom(0).IfNMRShieldingTensor()) {
	tmpNMR[iatom].Initialize(3,3);
	tmpNMR[iatom] = MonA.GetAtom(iatom).GetNMRShieldingTensor();
	//tmpNMR[iatom].Print("PLAIN NRMShieldingTensor");
      }
    }
    for (int iatom = 0;iatom < NatomsB; iatom++) { //  MonB
      if ( MonB.GetAtom(0).IfNMRShieldingTensor() ) {
	tmpNMR[iatom+NatomsA].Initialize(3,3);
	tmpNMR[iatom+NatomsA] = MonB.GetAtom(iatom).GetNMRShieldingTensor();
      }
    } 
    */
    MonA = m1; 
    MonB = m2; 

  } else {
    MonA = m1; 
    MonB = m2; 
  }
  
  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex(); 
  typeA = MonA.GetType(); 
  typeB = MonB.GetType();

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();

  // Set spin state.  Right now, it can handle 2 singlets, or singlet + anything.
  // Need something smarter to handle e.g. for 2 triplets
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];
  
}

void Dimer::UpdateObjects(Monomer& m1, Monomer& m2, Monomer& ref_m2) {

  MonA = m1; 
  MonB = m2; 

  indexA = MonA.GetIndex();
  indexB = MonB.GetIndex();
  int ref_indexB = ref_m2.GetIndex();
  typeA = MonA.GetType(); 
  typeB = MonB.GetType();

  // Set basic variables
  Natoms = MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms();
  charge = MonA.GetChargeState() + MonB.GetChargeState();

  // Set spin state.  Right now, it can handle 2 singlets, or singlet + anything.
  // Need something smarter to handle e.g. for 2 triplets
  int spinA = MonA.GetSpinState();
  int spinB = MonB.GetSpinState(); 
  if ( spinA==1 && spinB==1 ) spin = spinA;
  else if (spinA==1 && spinB!=1) spin = spinB;
  else if (spinB==1 && spinA!=1) spin = spinA;
  else {
    printf("Error: Dimer::Initialize() : Don't know how to handle combination of monomer spin states %d and %d\n",spinA,spinB);
    exit(1);
  }

  // Find minimum separation between two monomers
  Vector tmp = MonA.FindDistance(MonB);
  Separation = tmp[0];
  minA = (int) tmp[1];
  minB = (int) tmp[2];

  // Get energies/gradients from non-image monomer
  MonB.SetQMEnergy( ref_m2.GetQMEnergy() );
  MonB.SetMMEnergy( ref_m2.GetMMEnergy() );

  if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    //printf("Dimer::UpdateObjects(): setting gradients from reference monomer %d\n",ref_m2.GetIndex()) ;
    MonB.SetQMGradient( ref_m2 );
    MonB.SetMMGradient( ref_m2 );
  }
    
  if (Params::Parameters().GetJobTypeStr() == "hessian" && !Params::Parameters().DoFiniteDifferenceFreqs()) {
    //printf("Dimer::UpdateObjects(): setting hessians from reference monomer %d\n",ref_$
    MonB.SetQMHessian( ref_m2 );
    MonB.SetMMHessian( ref_m2 );
  }

  //Copy Rotational Matrix
  MonB.SetRotationMatrix(ref_m2.GetRotationMatrix());

  // If using the AIFF, need to get FF parameters
  if (Params::Parameters().GetMMType()==2) { 

    fflush(stdout);
    // Copy the local coordinates
    MonB.SetRotationAngle( ref_m2.GetRotationAngle() );
    double vec[3];
    for (int i=0;i<3;i++) 
      MonB.SetRotationVector(i,ref_m2.GetRotationVector()[i]);

    // Copy the multipole moments & polarizabilities
    int NatomsB = MonB.GetNumberOfAtoms();
    int Natoms_refB = ref_m2.GetNumberOfAtoms();

    if (NatomsB != Natoms_refB) {
      printf("Dimer::UpdateObjects: Image monomer and its reference have different numbers of atoms: %d and %d\n",NatomsB,Natoms_refB);
    }

    for (int iatom=0;iatom<NatomsB; iatom++) {    
      if (ref_m2.GetAtom(iatom).MultipolesAreInitialized() && ref_m2.GetAtom(iatom).PolarizabilityIsInitialized() ) {
	// Copy multipole moments, polarizabilities, and dispersion coefficients
	MonB.GetAtom(iatom).SetMultipoleMoments( ref_m2.GetAtom(iatom).GetMultipoleMoments() );
	MonB.GetAtom(iatom).SetPolarizability( ref_m2.GetAtom(iatom).GetPolarizability() );

	double iso_pol = ref_m2.GetAtom(iatom).GetIsotropicDipolePolarizability();
	MonB.GetAtom(iatom).SetIsotropicDipolePolarizability(iso_pol);

	double C6 = ref_m2.GetAtom(iatom).GetC6Coefficient();
	double C8 = ref_m2.GetAtom(iatom).GetC8Coefficient();
	double C10 = ref_m2.GetAtom(iatom).GetC10Coefficient();
	MonB.GetAtom(iatom).SetDispersionCoefficients( C6, C8, C10);
	MonB.GetAtom(iatom).SetDispersionAtomType(ref_m2.GetAtom(iatom).GetDispersionAtomType());

	// Copy Frequency-dependent polarizabilities
	MonB.GetAtom(iatom).SetFreq_Pol_Dipole(ref_m2.GetAtom(iatom).GetFreq_Pol_Dipole());
	MonB.GetAtom(iatom).SetFreq_Pol_Quad(ref_m2.GetAtom(iatom).GetFreq_Pol_Quad());


	//printf("iatom = %d\n",iatom);
	//MonB.GetAtom(iatom).GetMultipoleMoments().Print("Multipole moments\n");
	//MonB.GetAtom(iatom).GetPolarizability().Print("Polarizability\n");
	
      }
    }
  }
}

Vector Dimer::GetCurrentCoordinates(){
  
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  Vector Coords(3*(Na+Nb));

  //Get coordinates for Monomer A
  for(int i=0;i<Na;i++){
    for(int xyz=0;xyz<3;xyz++)
      Coords[3*i+xyz] = MonA.GetAtom(i).GetPosition(xyz);
  }

  //Get coordinates for Monomer B
  for(int i=0;i<Nb;i++){
    for(int xyz=0;xyz<3;xyz++)
      Coords[3*(Na+i)+xyz] = MonB.GetAtom(i).GetPosition(xyz);
  }
  
  return Coords;
}

//Symmetry matrix found in symmetry functions may be altered in
//Cluster::CorrectingDimerRotations()
bool Dimer::SymmetryCheck(Dimer Dimers[],int index, bool ImageToNonImage){
  MonA.FindCenterOfMass();
  MonB.FindCenterOfMass();
  //find the center of mass of dimer
  // printf("finding symmetry for Dimer (%i,%i)\n",indexA, indexB);
  Vector COM(3);//center of mass vector of dimer
  COM.Set();
  COM[0]= (MonA.GetCenterOfMass(0)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(0)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[1]= (MonA.GetCenterOfMass(1)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(1)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[2]= (MonA.GetCenterOfMass(2)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(2)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  //printf("center of mass is %10.6f  %10.6f  %10.6f\n",COM[0],COM[1],COM[2]);;
  //printf("number of dimers is %i\n",DimerNumberOfAtoms);
  
  //printf("number of atoms is dimer is %i\n", Natoms);
  
  //Finding center of mass coordinates and Inertia tensor
  Matrix COMCoord(Natoms,3);//center of mass coordinates for atoms
  Matrix Inertia(3,3);//Inertia tensor
  Inertia.Set();

  //printf("COM for D(%i,%i)\n",
  //MonA.GetIndex(),MonB.GetIndex());
  //MonA 	
  for(int i=0; i< MonA.GetNumberOfAtoms();i++){
    COMCoord(i,0) = MonA.GetAtom(i).GetPosition(0)-COM[0];
    COMCoord(i,1) = MonA.GetAtom(i).GetPosition(1)-COM[1];
    COMCoord(i,2) = MonA.GetAtom(i).GetPosition(2)-COM[2];
    double mass = MonA.GetAtom(i).GetAtomicMass();
    Inertia(0,0) += mass*(COMCoord(i,1)*COMCoord(i,1)+COMCoord(i,2)*COMCoord(i,2));//m*(y^2+z^2)
    Inertia(1,1) += mass*(COMCoord(i,0)*COMCoord(i,0)+COMCoord(i,2)*COMCoord(i,2));//m*(x^2+z^2)
    Inertia(2,2) += mass*(COMCoord(i,0)*COMCoord(i,0)+COMCoord(i,1)*COMCoord(i,1));//m*(x^2+z^2)
    Inertia(0,1) += -mass*(COMCoord(i,0)*COMCoord(i,1));//-m*x*y
    Inertia(0,2) += -mass*(COMCoord(i,0)*COMCoord(i,2));//-m*x*z
    Inertia(1,2) += -mass*(COMCoord(i,1)*COMCoord(i,2));//-m*x*z
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   MonA.GetAtom(i).GetSymbol().c_str(), COMCoord(i,0),COMCoord(i,1),COMCoord(i,2));	
  }
  //MonB 
  for(int i=0;i<MonB.GetNumberOfAtoms();i++){
    int j = MonA.GetNumberOfAtoms()+i;
    COMCoord(j,0)= MonB.GetAtom(i).GetPosition(0)-COM[0];
    COMCoord(j,1)= MonB.GetAtom(i).GetPosition(1)-COM[1];
    COMCoord(j,2)= MonB.GetAtom(i).GetPosition(2)-COM[2];
    double mass = MonB.GetAtom(i).GetAtomicMass();	
    Inertia(0,0) += mass*(COMCoord(j,1)*COMCoord(j,1)+COMCoord(j,2)*COMCoord(j,2));//m*(y^2+z^2)
    Inertia(1,1) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,2)*COMCoord(j,2));//m*(x^2+z^2)
    Inertia(2,2) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,1)*COMCoord(j,1));//m*(x^2+z^2)
    Inertia(0,1) += -mass*(COMCoord(j,0)*COMCoord(j,1));//-m*x*y
    Inertia(0,2) += -mass*(COMCoord(j,0)*COMCoord(j,2));//-m*x*z
    Inertia(1,2) += -mass*(COMCoord(j,1)*COMCoord(j,2));//-m*x*z
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   MonB.GetAtom(i).GetSymbol().c_str(), COMCoord(j,0),COMCoord(j,1),COMCoord(j,2));
  }
  Inertia(1,0) = Inertia(0,1);
  Inertia(2,0) = Inertia(0,2);
  Inertia(2,1) = Inertia(1,2);
  //diagonalizing
  Vector MomentsOfInertia = Inertia.Diagonalize();
  //printf("the moments of inertia are: %10.6f  %10.6f  %10.6f\n",MomentsOfInertia[0],MomentsOfInertia[1],MomentsOfInertia[2]);

  //new local axes	
  Vector x_axis = Inertia.GetColumnVector(0);
  Vector y_axis = Inertia.GetColumnVector(1);
  Vector z_axis = Inertia.GetColumnVector(2);
  //printf("Vectors of inertia are\n");
  //for(int i=0;i<3;i++){
  //	printf("%10.6f  %10.6f  %10.6f\n",x_axis[i],y_axis[i],z_axis[i]);
  //}
  
  //Dot product solutions
  Matrix LocalCoor(Natoms,3);
  string AtomList[Natoms]; 
  // printf("local coordinates\n");
  //printf("The new Internal coordinates are for dimer(%i,%i):\n",indexA,indexB);
  for(int i=0;i<MonA.GetNumberOfAtoms();i++){
    AtomList[i]= MonA.GetSymbol(i);
    LocalCoor(i,0) = COMCoord(i,0)*x_axis[0]+COMCoord(i,1)*x_axis[1]+COMCoord(i,2)*x_axis[2];
    LocalCoor(i,1) = COMCoord(i,0)*y_axis[0]+COMCoord(i,1)*y_axis[1]+COMCoord(i,2)*y_axis[2];
    LocalCoor(i,2) = COMCoord(i,0)*z_axis[0]+COMCoord(i,1)*z_axis[1]+COMCoord(i,2)*z_axis[2];
    //printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	AtomList[i].c_str(), LocalCoor(i,0),LocalCoor(i,1),LocalCoor(i,2));
  }
  for(int i=0;i<MonB.GetNumberOfAtoms();i++){
    int j = MonA.GetNumberOfAtoms()+i;
    AtomList[j]= MonB.GetSymbol(i);
    LocalCoor(j,0) = COMCoord(j,0)*x_axis[0]+COMCoord(j,1)*x_axis[1]+COMCoord(j,2)*x_axis[2];
    LocalCoor(j,1) = COMCoord(j,0)*y_axis[0]+COMCoord(j,1)*y_axis[1]+COMCoord(j,2)*y_axis[2];
    LocalCoor(j,2) = COMCoord(j,0)*z_axis[0]+COMCoord(j,1)*z_axis[1]+COMCoord(j,2)*z_axis[2];
    //printf( "%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	AtomList[j].c_str(), LocalCoor(j,0),LocalCoor(j,1),LocalCoor(j,2));
  }
  if(DetermineSymmetry(Dimers, LocalCoor, Inertia, AtomList, index,ImageToNonImage)){
    //printf("match for dimer(%i,%i) was found\n",indexA,indexB);
    return 1;
  }else{
    //printf("no match found\n");
    return 0;
  }
}
bool Dimer::DetermineSymmetry(Dimer Dimers[],Matrix LocalCoor,Matrix ThisInertia,
			      string AtomList[], int index, 
			      bool ImageToNonImage) {

  //second dimer
  double tolerence = Params::Parameters().GetSymmetryTolerance();
  //printf("tolerence is  %f\n",tolerence);
  
  for(int i=1;i<index;i++) {
    if((Natoms == Dimers[i].GetNumberOfAtoms()) && (Dimers[i].GetSymmetryFactor()!= 0)){
      Monomer SecMonA = Dimers[i].MonA;
      Monomer SecMonB = Dimers[i].MonB;
      //SecMonA.FindCenterOfMass();
      //SecMonB.FindCenterOfMass();
      //printf("compare to dimer (%i,%i)\n",Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
      Vector COM(3);
      COM[0] = (SecMonA.GetCenterOfMass(0)*SecMonA.GetMonomerMass()+SecMonB.GetCenterOfMass(0)*SecMonB.GetMonomerMass())/(SecMonA.GetMonomerMass()+SecMonB.GetMonomerMass());
      COM[1]= (SecMonA.GetCenterOfMass(1)*SecMonA.GetMonomerMass()+SecMonB.GetCenterOfMass(1)*SecMonB.GetMonomerMass())/(SecMonA.GetMonomerMass()+SecMonB.GetMonomerMass());
      COM[2]= (SecMonA.GetCenterOfMass(2)*SecMonA.GetMonomerMass()+SecMonB.GetCenterOfMass(2)*SecMonB.GetMonomerMass())/(SecMonA.GetMonomerMass()+SecMonB.GetMonomerMass());
      
      
      //finding moment of center of mass coordinates and moment of Inertia
      Matrix COMCoord(Natoms,3);
      Matrix InertiaSym(3,3);
      InertiaSym.Set();
      for(int j=0; j<SecMonA.GetNumberOfAtoms();j++){
	COMCoord(j,0) = SecMonA.GetAtom(j).GetPosition(0)-COM[0];
	COMCoord(j,1) = SecMonA.GetAtom(j).GetPosition(1)-COM[1];
	COMCoord(j,2) = SecMonA.GetAtom(j).GetPosition(2)-COM[2];
	double mass = SecMonA.GetAtom(j).GetAtomicMass();	
	InertiaSym(0,0) += mass*(COMCoord(j,1)*COMCoord(j,1)+COMCoord(j,2)*COMCoord(j,2));//m*(y^2+z^2)
	InertiaSym(1,1) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,2)*COMCoord(j,2));//m*(x^2+z^2)
	InertiaSym(2,2) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,1)*COMCoord(j,1));//m*(x^2+z^2)
	InertiaSym(0,1) += -mass*(COMCoord(j,0)*COMCoord(j,1));//-m*x*y
	InertiaSym(0,2) += -mass*(COMCoord(j,0)*COMCoord(j,2));//-m*x*z
	InertiaSym(1,2) += -mass*(COMCoord(j,1)*COMCoord(j,2));//-m*x*z
	//printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
	//      SecMonA.GetSymbol(j).c_str(), COMCoord(j,0),COMCoord(j,1),COMCoord(j,2));
      }
      //next for MonB adding the values for MonB  
      for(int j=0;j<SecMonB.GetNumberOfAtoms();j++){
	int k = SecMonA.GetNumberOfAtoms()+j;
	COMCoord(k,0)= SecMonB.GetAtom(j).GetPosition(0)-COM[0];
	COMCoord(k,1)= SecMonB.GetAtom(j).GetPosition(1)-COM[1];
	COMCoord(k,2)= SecMonB.GetAtom(j).GetPosition(2)-COM[2];
	double mass = SecMonB.GetAtom(j).GetAtomicMass();	
	InertiaSym(0,0) += mass*(COMCoord(k,1)*COMCoord(k,1)+COMCoord(k,2)*COMCoord(k,2));//m*(y^2+z^2)
	InertiaSym(1,1) += mass*(COMCoord(k,0)*COMCoord(k,0)+COMCoord(k,2)*COMCoord(k,2));//m*(x^2+z^2)
	InertiaSym(2,2) += mass*(COMCoord(k,0)*COMCoord(k,0)+COMCoord(k,1)*COMCoord(k,1));//m*(x^2+z^2)
	InertiaSym(0,1) += -mass*(COMCoord(k,0)*COMCoord(k,1));//-m*x*y
	InertiaSym(0,2) += -mass*(COMCoord(k,0)*COMCoord(k,2));//-m*x*z
	InertiaSym(1,2) += -mass*(COMCoord(k,1)*COMCoord(k,2));//-m*x*z
	//printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
	//      SecMonB.GetSymbol(j).c_str(), COMCoord(k,0),COMCoord(k,1),COMCoord(k,2));		
      }
      InertiaSym(1,0) = InertiaSym(0,1);
      InertiaSym(2,0) = InertiaSym(0,2);
      InertiaSym(2,1) = InertiaSym(1,2);

      //diagonalizing
      Vector MomentsOfInertia = InertiaSym.Diagonalize();
      Vector x_axis = InertiaSym.GetColumnVector(0);
      Vector y_axis = InertiaSym.GetColumnVector(1);
      Vector z_axis = InertiaSym.GetColumnVector(2);
      //Dot product solutions
      Matrix LocalCoor2(Natoms,3);
      string AtomList2[Natoms];
      for(int j=0;j<SecMonA.GetNumberOfAtoms();j++){
	AtomList2[j] = SecMonA.GetSymbol(j);
	LocalCoor2(j,0) = COMCoord(j,0)*x_axis[0]+COMCoord(j,1)*x_axis[1]+COMCoord(j,2)*x_axis[2];
	LocalCoor2(j,1) = COMCoord(j,0)*y_axis[0]+COMCoord(j,1)*y_axis[1]+COMCoord(j,2)*y_axis[2];
	LocalCoor2(j,2) = COMCoord(j,0)*z_axis[0]+COMCoord(j,1)*z_axis[1]+COMCoord(j,2)*z_axis[2];
	//printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
	//       SecMonA.GetSymbol(j).c_str(), LocalCoor2(j,0),LocalCoor2(j,1),LocalCoor2(j,2));
      }
      for(int j=0;j<SecMonB.GetNumberOfAtoms();j++){
	int k = SecMonA.GetNumberOfAtoms()+j;
	AtomList2[k] = SecMonB.GetSymbol(j);
	LocalCoor2(k,0) = COMCoord(k,0)*x_axis[0]+COMCoord(k,1)*x_axis[1]+COMCoord(k,2)*x_axis[2];
	LocalCoor2(k,1) = COMCoord(k,0)*y_axis[0]+COMCoord(k,1)*y_axis[1]+COMCoord(k,2)*y_axis[2];
	LocalCoor2(k,2) = COMCoord(k,0)*z_axis[0]+COMCoord(k,1)*z_axis[1]+COMCoord(k,2)*z_axis[2];
	//printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
	//        SecMonB.GetSymbol(j).c_str(), LocalCoor2(k,0),LocalCoor2(k,1),LocalCoor2(k,2));
      }


      //Detmining Symmetry
      //these are used to preform reflection
      for(int XReflex=1;XReflex >=-1;XReflex-=2){
	for(int YReflex=1;YReflex >=-1;YReflex-=2){
	  for(int ZReflex=1;ZReflex >=-1;ZReflex-=2){

	    //Equivalent atoms on the dimer;
	    Vector equivalency(Natoms);

	    //printf("%i, %i, %i\n",XReflex,YReflex,ZReflex);
	    Matrix Reflex(3,3);
	    Reflex.Set();
	    Reflex(0,0) = XReflex;
	    Reflex(1,1) = YReflex;
	    Reflex(2,2) = ZReflex;
	    Matrix OperatedLocal2= LocalCoor2.Multiply(Reflex);
	    //LocalCoor2.Print("LocalCoor2");
	    //OperatedLocal2.Print("OperatedLocal2");
	    //LocalCoor.Print("LocalCoor");


	    int j=0;//counter for atoms of molecule one
	    bool nosym = 0;//one if molecule are not identical
	    while(j<Natoms && nosym == 0){//this is a loop over the atoms of dimer 1
	      bool match = 0;//this indicates whether there is a matching atom
	      int k = 0;
	      while(k<Natoms && match ==0){//looping over atoms of second monomer
		if(AtomList[j] == AtomList2[k]){
		  double x_diff = fabs(LocalCoor(j,0) - OperatedLocal2(k,0));
		  double y_diff = fabs(LocalCoor(j,1) - OperatedLocal2(k,1));
		  double z_diff = fabs(LocalCoor(j,2) - OperatedLocal2(k,2));
		  if(x_diff < tolerence && y_diff < tolerence && z_diff < tolerence){
		    match = 1;//identical atom found
                    equivalency[k] = j+1;

		    //printf("Match for on dimer(%i,%i) %f %f %f\n"
		    //	   ,Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),OperatedLocal2(k,0),OperatedLocal2(k,1),OperatedLocal2(k,2));
		    //printf(" dimer(%i,%i) %f %f %f\n"
		    //	   ,indexA,indexB,LocalCoor(j,0),LocalCoor(j,1),LocalCoor(j,2));
		    //printf("x_diff = %f,y_diff = %f,z_diff = %f\n",x_diff,y_diff,z_diff);
		    //printf("LocalCoor(j,2) - OperatedLocal2(k,2)) = %f\n", LocalCoor(j,2) - OperatedLocal2(k,2));

		}

		}
		k++;
	      }		
	      if(match == 0)
		nosym = 1;
	      j++;
	    }				
	    if(nosym == 0){
	      if(!ImageToNonImage){


		Dimers[i].MergeEquivalencyList(GetAtomEquivalency(),equivalency,false);
		
		//Create Rotation Matrix
		Matrix Rotation = Reflex.Multiply(ThisInertia,4);
		Rotation = InertiaSym.Multiply(Rotation);


		//if the first dimer is an image dimer and the second is not
		//setting symmetry factor for symmetrical dimers where both 
		//the dimers are either both image or non-image dimers
		int sym = sym_fac + Dimers[i].GetSymmetryFactor();
		Dimers[i].SetSymmetryFactor(sym);
		AddToSymmetryList(Dimers[i],Rotation,ImageToNonImage);

		//printf("D(%i,%i) symmetrical to D(%i,%i)\n",
		//       indexA,indexB,Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
		  
	      }else{



		Dimers[i].MergeEquivalencyList(GetAtomEquivalency(),equivalency,true);

		//Create Rotation Matrix
		Matrix Rotation = Reflex.Multiply(ThisInertia,4);
		Rotation = InertiaSym.Multiply(Rotation);

		//setting the symmetry factor for symmetrical dimers 
	        //where one of the dimers is an image not the other is not 
		int sym = sym_fac + Dimers[i].GetPeriodicSymmetryFactor();
		Dimers[i].SetPeriodicSymmetryFactor(sym);
		AddToSymmetryList(Dimers[i],Rotation,ImageToNonImage);

		//printf("D(%i,%i) symmetrical to D(%i,%i)\n",
		//      indexA,indexB,Dimers[i].GetIndexA(),Dimers[i].GetIndexB());
	      }
	      sym_fac=0;
	      return 1;
	    }
	  }
	}
      }
    }
  }
  return 0;
}

//This function is mostly out of date but does contain two subfunctions that are still used.
void Dimer::AddToSymmetryList(Dimer& SymDimer,Matrix Rot,bool ImageToNonImage){

  //printf("Determining which monomer MonA in the dimer is symmetrical to.\n"); 


  /*
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  Vector COM(3);
  Matrix COMA(3,Na);
  //finding center of mass for dimer
  COM[0] = (MonA.GetCenterOfMass(0)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(0)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[1]= (MonA.GetCenterOfMass(1)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(1)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[2]= (MonA.GetCenterOfMass(2)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(2)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());

  //printf("Center of mass of MonA\n");
  for(int i=0;i<Na;i++){
    COMA(0,i) = MonA.GetAtom(i).GetPosition(0)-COM[0];
    COMA(1,i) = MonA.GetAtom(i).GetPosition(1)-COM[1];
    COMA(2,i) = MonA.GetAtom(i).GetPosition(2)-COM[2];
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   MonA.GetSymbol(i).c_str(), COMA(0,i),COMA(1,i),COMA(2,i));
  }

  Monomer SymA = SymDimer.GetMonomerA();
  Monomer SymB = SymDimer.GetMonomerB();
  int SymNa = SymA.GetNumberOfAtoms();
  int SymNb = SymB.GetNumberOfAtoms();

  //finding center of mass for SymDimer
  Vector SymCOM(3);
  Matrix SymCOMA(3,SymNa);
  Matrix SymCOMB(3,SymNb);
  //Matrix SymRotA = SymA.GetRotationMatrix();
  //Matrix SymRotB = SymB.GetRotationMatrix();

  //finding center of mass for SymDimer
  SymCOM[0] = (SymA.GetCenterOfMass(0)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(0)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());
  SymCOM[1] = (SymA.GetCenterOfMass(1)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(1)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());
  SymCOM[2] = (SymA.GetCenterOfMass(2)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(2)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());

  //printf("Center of mass of SymMonA\n");
  for(int i=0;i<SymNa;i++){
    SymCOMA(0,i) = SymA.GetAtom(i).GetPosition(0)-SymCOM[0];
    SymCOMA(1,i) = SymA.GetAtom(i).GetPosition(1)-SymCOM[1];
    SymCOMA(2,i) = SymA.GetAtom(i).GetPosition(2)-SymCOM[2];
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   SymA.GetSymbol(i).c_str(), SymCOMA(0,i),SymCOMA(1,i),SymCOMA(2,i));
  }
  //printf("Center of mass of SymMonB\n");
  for(int i=0;i<SymNb;i++){
    SymCOMB(0,i) = SymB.GetAtom(i).GetPosition(0)-SymCOM[0];
    SymCOMB(1,i) = SymB.GetAtom(i).GetPosition(1)-SymCOM[1];
    SymCOMB(2,i) = SymB.GetAtom(i).GetPosition(2)-SymCOM[2];
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   SymB.GetSymbol(i).c_str(), SymCOMB(0,i),SymCOMB(1,i),SymCOMB(2,i));
  }

  //Making rotation matrix from the product of the matrices and checking to see if it is the correct rotation
  //Matrix RotAA = SymRotA.Multiply(RotA,2);
  //Matrix RotBA = SymRotB.Multiply(RotA,2);
  bool notsymA = false;//determine if MonA is symmetrical equivalent in the dimer to SymA
  bool notsymB = false;//determine if MonA is symmetrical equivalent in the dimer to SymB


  //Applying rotation to determine correct way to rotate MonA
  Matrix Rot_COM = Rot.Multiply(COMA);

 
  double tolerence = Params::Parameters().GetSymmetryTolerance();



  //checking to see if monA in the dimer is symmetrical to SymA in the dimer
  int i=0;//counter for atoms of molecule one
  if(Na!=SymNa)
    notsymA = true;
  while(i<Na && notsymA == false){//this is a loop over the atoms of dimer 1
    bool match = 0;//this indicates whether there is a matching atom
    int j = 0;
    while(j<Na && match ==0){//looping over atoms of second monomer
      if(MonA.GetSymbol(i)==SymA.GetSymbol(j)){
	//double x_diff = fabs(Rot_COMA(0,i) - SymCOMA(0,j));
	//double y_diff = fabs(Rot_COMA(1,i) - SymCOMA(1,j));
	//double z_diff = fabs(Rot_COMA(2,i) - SymCOMA(2,j));
	double x_diff = fabs(Rot_COM(0,i) - SymCOMA(0,j));
	double y_diff = fabs(Rot_COM(1,i) - SymCOMA(1,j));
	double z_diff = fabs(Rot_COM(2,i) - SymCOMA(2,j));
	if(x_diff < tolerence && y_diff < tolerence && z_diff < tolerence){
	  match = true;//identical atom found
	  //printf("SymA found be between %i and %i\n",i,j);

	}
      }
      j++;
    }
    if(match == 0)
      notsymA = true;
    i++;
  }

  
  //checking to see if monA in the dimer is symmetrical equivalent to SymB in the dimer if not symmetrical to SymA
  if(notsymA){
    int i=0;//counter for atoms of molecule one
    if(Na!=SymNb)
      notsymB = true;
    while(i<Na && notsymB == false){
      //printf("Looking at atom %i Symbol = %s\n",i,MonA.GetSymbol(i).c_str());
      bool match = 0;//this indicates whether there is a matching atom
      int j = 0;
      while(j<Na && match ==0){//looping over atoms of second monomer
	//printf("Does atom %i Symbol = %s match\n",j,SymB.GetSymbol(j).c_str());
	if(MonA.GetSymbol(i)==SymB.GetSymbol(j)){
	  //printf("Rot_COMB %f  SymCOMB %f\n",Rot_COMB(0,i), SymCOMB(0,j));
	  //double x_diff = fabs(Rot_COMB(0,i) - SymCOMB(0,j));
	  //double y_diff = fabs(Rot_COMB(1,i) - SymCOMB(1,j));
	  //double z_diff = fabs(Rot_COMB(2,i) - SymCOMB(2,j));
	  double x_diff = fabs(Rot_COM(0,i) - SymCOMB(0,j));
	  double y_diff = fabs(Rot_COM(1,i) - SymCOMB(1,j));
	  double z_diff = fabs(Rot_COM(2,i) - SymCOMB(2,j));
	  if(x_diff < tolerence && y_diff < tolerence && z_diff < tolerence){
	  match = true;//identical atom found
	  //printf("SymB found be between %i and %i\n",i,j);
	  }
	}
	j++;
      }
      if(match == 0){
	printf("No match for for atom %i\n",i);
	notsymB = true;
      }
      i++; 
    }
  }
  */

   //Merge Symmetrical list of the two dimers together.
  //If MonA in the dimer in symmetrically equivalent to SymA in SymDimer, invert this dimer symmetry list 
  //if(!notsymA ||!notsymB){
  SymDimer.AddToRotationList(Rot,Rotation,ImageToNonImage,false);
  SymDimer.MergeSymmetryList(List_Sym_Dim,/*MonB_List,*/k_list,ImageToNonImage);

    /*
  }else{
    printf("Error::Dimer::AddToSymmetryList(): Cannot determine Symmetry between d(%i,%i) and d(%i,%i)\n",
	   indexA,indexB,SymDimer.GetIndexA(),SymDimer.GetIndexB());
    Rot.Print("Rot");
    COMA.Transpose();
    SymCOMA.Transpose();
    SymCOMB.Transpose();

    COMA.Print("COMA");
    SymCOMA.Print("SymCOMA");
    SymCOMB.Print("SymCOMB");

    exit(0);
  }
    */
}

//add to the RotationList by copying it. The symfac should be increased before using this function
void Dimer::AddToRotationList(Matrix Rot, vector<Matrix> OtherMatrixList, bool ImageToNonImage,bool InverseList){

  
  if(!ImageToNonImage){
    if(!InverseList)
      for(int i=0;i<OtherMatrixList.size();i++){
	Matrix NewRot = Rot.Multiply(OtherMatrixList[i]);
	Rotation.push_back(NewRot);
      }
    else
      for(int i=OtherMatrixList.size()-1;i>0;i--){
	Matrix NewRot = Rot.Multiply(OtherMatrixList[i]);
	Rotation.push_back(NewRot);
      }
  }
  else{
    if(!InverseList){
      for(int i=0;i<OtherMatrixList.size();i++){
	Matrix NewRot = Rot.Multiply(OtherMatrixList[i]);
	Rotation_Periodic.push_back(NewRot);
      }
    }
    else{
      for(int i=OtherMatrixList.size()-1;i>0;i--){
	Matrix NewRot = Rot.Multiply(OtherMatrixList[i]);
	Rotation_Periodic.push_back(NewRot);
      }
    }
  }


}

//merges atom equivalancy lists
//Symmetry List may be altered in Cluster::CorrectingMonomerRotations()
void Dimer::MergeEquivalencyList(vector<Vector> OtherList,Vector ThisList, bool ImageToNonImage){

  //if(indexA==1 && indexB==2 && ImageToNonImage){
   // printf("D(%i,%i)\n",indexA,indexB);
   // ThisList.Print("ThisList");
    
    //for(int i=0;i<OtherList.size();i++)
     // OtherList[i].Print("OtherList");
  //}

  //rearrange OtherList so it the atoms on it are equivalent to atoms on this dimer
  for(int i=0;i<OtherList.size();i++){
    Vector Equivalency(OtherList[i].GetLength());
    for(int j=0;j<OtherList[i].GetLength();j++)
      for(int k=0;k<ThisList.GetLength();k++)
	if(OtherList[i][j] == ThisList[k]){
	  Equivalency[k] = j+1;
	}
    //if(indexA==1 && indexB==2 && ImageToNonImage){
      //Equivalency.Print("Equivalency");
    //}

    //making sure every value in Equivalency is filled (non-zero)
    int nonzero = Equivalency.Nonzero();
    if(nonzero < Equivalency.GetLength()){
      printf("Dimer::MergeEquivalencyList:: Unable to merge atom equivalency lists");
      exit(0);
    }

    if(!ImageToNonImage)
      Atom_Equivalency.push_back(Equivalency);
    else
      Atom_Equivalency_Periodic.push_back(Equivalency);
  }	

}

//this function merges both the dimer symmetry list and the MonB_List which list which dimers on the list are outside the unit cell
void Dimer::MergeSymmetryList(vector<int> OtherSymList/*,vector<bool> OtherMonBList*/,vector<int> OtherKList, bool ImageToNonImage){

 
  //if both dimers are image or non image dimers.
  if(!ImageToNonImage){
      for(int i=0;i<OtherSymList.size();i++)
	List_Sym_Dim.push_back(OtherSymList[i]);
    //If one of the dimers which symmetry list is being an image dimer and the other is not
  }else{
      for(int i=0;i<OtherSymList.size();i++)
	List_Sym_Periodic.push_back(OtherSymList[i]);

  }

  //printf("D(%i,%i) updated List_Sym_Dim ",indexA,indexB);
  //for(int i=0;i<List_Sym_Dim.size();i++)
  //  printf("%i ",List_Sym_Dim[i]);
  //printf("\n");

  //merging K_List
  //This list will tell us which monomers on the symmetry lists are outside the unit cell
  //if(!InverseList)
    for(int i=0; i<OtherKList.size();i++){
      k_list.push_back(OtherKList[i]);
      //MonB_List.push_back(OtherMonBList[i]);
      //k_list.push_back(OtherKList[3*i]);
      //k_list.push_back(OtherKList[3*i+1]);
      //k_list.push_back(OtherKList[3*i+2]);
    }
    /*
  else
      for(int i=OtherMonBList.size()-1; i>=0 ;i--){
	if(OtherMonBList[i]==1)
	  MonB_List.push_back(0);
	else
	  MonB_List.push_back(1);
	k_list.push_back(OtherKList[3*i]);
	k_list.push_back(OtherKList[3*i+1]);
	k_list.push_back(OtherKList[3*i+2]);
      }
    */
 
}

void Dimer::ResetSymmetryList(){
  List_Sym_Dim.clear();
  List_Sym_Dim.push_back(indexA);
  List_Sym_Dim.push_back(indexB);

  List_Sym_Periodic.clear();
  //MonB_List.clear();
  k_list.clear();
  //printf("D(%i,%i): after reset size of List_Sym_Dim %i\n",
  //	 indexA,indexB,List_Sym_Dim.size());
  //printf("D(%i,%i) reset List_Sym_Dim to ",indexA,indexB);
  //for(int i=0;i<List_Sym_Dim.size();i++)
  //  printf("%i ",List_Sym_Dim[i]);
  //printf("\n");
}


//After correcting the rotation matrix, determine atom equivalency base on that 
void Dimer::AlterSymmetryList(Dimer& SymDimer,int ListNumb,Matrix Rot,bool ImageToNonImage){

  


  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  Vector COM(3);
  Matrix CoordCOM(3,Na+Nb);
  string AtomList[Na+Nb];


  //finding center of mass for dimer
  COM[0] = (MonA.GetCenterOfMass(0)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(0)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[1]= (MonA.GetCenterOfMass(1)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(1)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[2]= (MonA.GetCenterOfMass(2)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(2)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());

  //printf("Center of mass of MonA\n");
  for(int i=0;i<Na;i++){
    AtomList[i] = MonA.GetAtom(i).GetSymbol();
    CoordCOM(0,i) = MonA.GetAtom(i).GetPosition(0)-COM[0];
    CoordCOM(1,i) = MonA.GetAtom(i).GetPosition(1)-COM[1];
    CoordCOM(2,i) = MonA.GetAtom(i).GetPosition(2)-COM[2];



    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   MonA.GetSymbol(i).c_str(), COMA(0,i),COMA(1,i),COMA(2,i));
  }

  for(int i=0;i<Nb;i++){
    int j = Na +i;
    AtomList[j] = MonB.GetAtom(i).GetSymbol();
    CoordCOM(0,j) = MonB.GetAtom(i).GetPosition(0)-COM[0];
    CoordCOM(1,j) = MonB.GetAtom(i).GetPosition(1)-COM[1];
    CoordCOM(2,j) = MonB.GetAtom(i).GetPosition(2)-COM[2];

  }

  Monomer SymA = SymDimer.GetMonomerA();
  Monomer SymB = SymDimer.GetMonomerB();
  int SymNa = SymA.GetNumberOfAtoms();
  int SymNb = SymB.GetNumberOfAtoms();

  //finding center of mass for SymDimer
  string SymAtomList[SymNa+SymNb];
  Vector SymCOM(3);
  Matrix SymCoordCOM(3,SymNa+SymNb);

  //finding center of mass for SymDimer
  SymCOM[0] = (SymA.GetCenterOfMass(0)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(0)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());
  SymCOM[1] = (SymA.GetCenterOfMass(1)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(1)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());
  SymCOM[2] = (SymA.GetCenterOfMass(2)*SymA.GetMonomerMass()+SymB.GetCenterOfMass(2)*SymB.GetMonomerMass())/(SymA.GetMonomerMass()+SymB.GetMonomerMass());

  //printf("Center of mass of SymMonA\n");
  for(int i=0;i<SymNa;i++){
    SymAtomList[i] = SymA.GetAtom(i).GetSymbol();
    SymCoordCOM(0,i) = SymA.GetAtom(i).GetPosition(0)-SymCOM[0];
    SymCoordCOM(1,i) = SymA.GetAtom(i).GetPosition(1)-SymCOM[1];
    SymCoordCOM(2,i) = SymA.GetAtom(i).GetPosition(2)-SymCOM[2];

  }
  //printf("Center of mass of SymMonB\n");
  for(int i=0;i<SymNb;i++){
    int j = Na + i;
    SymAtomList[j] = SymB.GetAtom(i).GetSymbol();
    SymCoordCOM(0,j) = SymB.GetAtom(i).GetPosition(0)-SymCOM[0];
    SymCoordCOM(1,j) = SymB.GetAtom(i).GetPosition(1)-SymCOM[1];
    SymCoordCOM(2,j) = SymB.GetAtom(i).GetPosition(2)-SymCOM[2];
  }

  //Making rotation matrix from the product of the matrices and checking to see if it is the correct rotation
  //Matrix RotAA = SymRotA.Multiply(RotA,2);
  //Matrix RotBA = SymRotB.Multiply(RotA,2);

  
  Matrix RotCoord = Rot.Multiply(CoordCOM);

  /*
  CoordCOM.Transpose();
  RotCoord.Transpose();
  SymCoordCOM.Transpose();
  Rot.Print("\nRot");
  CoordCOM.Print("CoordCOM");
  RotCoord.Print("RotCoord");
  SymCoordCOM.Print("SymCoordCOM");
  CoordCOM.Transpose();
  RotCoord.Transpose();
  SymCoordCOM.Transpose();
  */

  //bool SymAFound = false;//equivalency monomer for MonA found
  // bool SymBFound = false;//equivalency monomer for MonB found


  Vector AtomEquivalency(Na+Nb);

  double tolerence = Params::Parameters().GetSymmetryTolerance();

  int j=0;//counter for atoms of molecule one
  bool nosym = 0;//Symmetry cannot be found
  while(j<Na+Nb && nosym == 0){//this is a loop over the atoms of dimer 1
    bool match = 0;//this indicates whether there is a matching atom
    int k = 0;
    while(k<Na+Nb && match ==0){//looping over atoms of second monomer
      if(AtomList[j] == SymAtomList[k]){
	double x_diff = fabs(SymCoordCOM(0,j) - RotCoord(0,k));
	double y_diff = fabs(SymCoordCOM(1,j) - RotCoord(1,k));
	double z_diff = fabs(SymCoordCOM(2,j) - RotCoord(2,k));
	if(x_diff < tolerence && y_diff < tolerence && z_diff < tolerence){
	  match = 1;//identical atom found
	  AtomEquivalency[k] = j+1;
	  
	  //printf("Match for on dimer(%i,%i) %f %f %f\n"
	  //	   ,Dimers[i].GetIndexA(),Dimers[i].GetIndexB(),OperatedLocal2(k,0),OperatedLocal2(k,1),OperatedLocal2(k,2));
	  // printf(" dimer(%i,%i) %f %f %f\n"
	  //   ,indexA,indexB,LocalCoor(j,0),LocalCoor(j,1),LocalCoor(j,2));
	  //printf("x_diff = %f,y_diff = %f,z_diff = %f\n",x_diff,y_diff,z_diff);
	  //printf("LocalCoor(j,2) - OperatedLocal2(k,2)) = %f\n", LocalCoor(j,2) - OperatedLocal2(k,2));
	  
	}
	
      }
      k++;
    }		
    if(match == 0)
      nosym = 1;
    j++;
  }

  if(nosym){
    printf("Dimer::AlterSymmetryList() :: Cannot correct dimer (%i,%i) symmetry list\n",indexA,indexB);
    exit(0);
  }


  //Replacing rotation and atom equivalency
  if(ImageToNonImage){
    //printf("Not ImageToNonImage");
    //Atom_Equivalency_Periodic.Print("AtomEquivalency");
    Rotation_Periodic[ListNumb] = Rot;
    Atom_Equivalency_Periodic[ListNumb] = AtomEquivalency;
    //Atom_Equivalency_Periodic()[ListNumb].Print("AtomEquivalency");
    //exit(0);

  }
  else{
    //AtomEquivalency.Print("AtomEquivalency");
    //printf("Not ImageToNonImage");
    //printf("ListNumb = %i\n",ListNumb);
    //Rotation[ListNumb].Print("Rotation before");
    //Rot.Print("Rot");
    Rotation[ListNumb] = Rot;
    //Rotation[ListNumb].Print("Rotation after");
    Atom_Equivalency[ListNumb] = AtomEquivalency;
    //GetAtomEquivalency()[ListNumb].Print("AtomEquivalency");

  }





}

bool Dimer::IsDimerPlanarOrLinear(){
  double tolerance = Params::Parameters().GetSymmetryTolerance();//Tolerance 

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  MonA.FindCenterOfMass();
  MonB.FindCenterOfMass();
  //find the center of mass of dimer
  Vector COM(3);//center of mass vector of dimer
  COM.Set();
  COM[0]= (MonA.GetCenterOfMass(0)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(0)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[1]= (MonA.GetCenterOfMass(1)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(1)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  COM[2]= (MonA.GetCenterOfMass(2)*MonA.GetMonomerMass()+MonB.GetCenterOfMass(2)*MonB.GetMonomerMass())/(MonA.GetMonomerMass()+MonB.GetMonomerMass());
  //printf("center of mass is %10.6f  %10.6f  %10.6f\n",COM[0],COM[1],COM[2]);;
  //printf("number of dimers is %i\n",DimerNumberOfAtoms);
  
  //printf("number of atoms is dimer is %i\n", Natoms);
  
  //Finding center of mass coordinates and Inertia tensor
  Matrix COMCoord(3,Natoms);//center of mass coordinates for atoms
  

  //MonA 	
  for(int i=0; i< Na;i++){
    COMCoord(0,i) = MonA.GetAtom(i).GetPosition(0)-COM[0];
    COMCoord(1,i) = MonA.GetAtom(i).GetPosition(1)-COM[1];
    COMCoord(2,i) = MonA.GetAtom(i).GetPosition(2)-COM[2];
    //double mass = MonA.GetAtom(i).GetAtomicMass();
    //Inertia(0,0) += mass*(COMCoord(i,1)*COMCoord(i,1)+COMCoord(i,2)*COMCoord(i,2));//m*(y^2+z^2)
    //Inertia(1,1) += mass*(COMCoord(i,0)*COMCoord(i,0)+COMCoord(i,2)*COMCoord(i,2));//m*(x^2+z^2)
    //Inertia(2,2) += mass*(COMCoord(i,0)*COMCoord(i,0)+COMCoord(i,1)*COMCoord(i,1));//m*(x^2+z^2)
    //Inertia(0,1) += -mass*(COMCoord(i,0)*COMCoord(i,1));//-m*x*y
    //Inertia(0,2) += -mass*(COMCoord(i,0)*COMCoord(i,2));//-m*x*z
    //Inertia(1,2) += -mass*(COMCoord(i,1)*COMCoord(i,2));//-m*y*z
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    //	   MonA.GetAtom(i).GetSymbol().c_str(), COMCoord(i,0),COMCoord(i,1),COMCoord(i,2));	
  }
  //MonB 
  for(int i=0;i<Nb;i++){
    int j = Na+i;
    COMCoord(0,j)= MonB.GetAtom(i).GetPosition(0)-COM[0];
    COMCoord(1,j)= MonB.GetAtom(i).GetPosition(1)-COM[1];
    COMCoord(2,j)= MonB.GetAtom(i).GetPosition(2)-COM[2];
    //double mass = MonB.GetAtom(i).GetAtomicMass();	
    //Inertia(0,0) += mass*(COMCoord(j,1)*COMCoord(j,1)+COMCoord(j,2)*COMCoord(j,2));//m*(y^2+z^2)
    //Inertia(1,1) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,2)*COMCoord(j,2));//m*(x^2+z^2)
    //Inertia(2,2) += mass*(COMCoord(j,0)*COMCoord(j,0)+COMCoord(j,1)*COMCoord(j,1));//m*(x^2+z^2)
    //Inertia(0,1) += -mass*(COMCoord(j,0)*COMCoord(j,1));//-m*x*y
    //Inertia(0,2) += -mass*(COMCoord(j,0)*COMCoord(j,2));//-m*x*z
    //Inertia(1,2) += -mass*(COMCoord(j,1)*COMCoord(j,2));//-m*y*z
    //printf("%-2s  %10.6f  %10.6f  %10.6f\n", 
    // MonB.GetAtom(i).GetSymbol().c_str(), COMCoord(j,0),COMCoord(j,1),COMCoord(j,2));
  }
  //Inertia(1,0) = Inertia(0,1);
  //Inertia(2,0) = Inertia(0,2);
  //Inertia(2,1) = Inertia(1,2);
  //diagonalizing
  //Inertia.Print("Inertia");
  //Vector MomentsOfInertia = Inertia.Diagonalize();
  

  //printf("d(%i,%i)\n",MonA.GetIndex(),MonB.GetIndex());
  //printf("MomentsOfInertia = %f %f %f\n",MomentsOfInertia[0],MomentsOfInertia[1],MomentsOfInertia[2]);
  //exit(0);
  //bool AsymmetricalTop;
  //if(fabs(MomentsOfInertia[2]-MomentsOfInertia[0]) < tolerance ||
  //   fabs(MomentsOfInertia[1]-MomentsOfInertia[0]) < tolerance  ||
  //   fabs(MomentsOfInertia[2]-MomentsOfInertia[1]) < tolerance){
  //  AsymmetricalTop =  false;
  //}
//else 
//  AsymmetricalTop = true;
//printf("AsymmetricalTop = %i\n",AsymmetricalTop);

  bool PlanarOrLinear = 1;

  //Check if planar or linear using determent of atoms treated as vectors

  Matrix ThreeAtomCoord(3,3);//Coordinates of three atoms
  ThreeAtomCoord.Set();
  for(int iatom=0;iatom < Na+Nb;iatom++){
    for(int jatom=iatom+1;jatom <Na+Nb;jatom++){
      for(int katom=jatom+1;katom<Na+Nb;katom++){
	Matrix ThreeAtomCoord(3,3);//Coordinates of three atoms
	ThreeAtomCoord.Set();
	ThreeAtomCoord.SetColumnVector(COMCoord.GetColumnVector(iatom),0);
	ThreeAtomCoord.SetColumnVector(COMCoord.GetColumnVector(jatom),1);
	ThreeAtomCoord.SetColumnVector(COMCoord.GetColumnVector(katom),2);

	//printf("iatom = %i jatom = %i katom =%i\n",iatom,jatom,katom);

	//if(GetIndexA() == 3 && GetIndexB() == 6)
	//ThreeAtomCoord.Print("ThreeAtomCoord");
	if(fabs(ThreeAtomCoord.Determinant()) > tolerance){
	  //if(GetIndexA() == 3 && GetIndexB() == 6)
	    //printf("Dimer is not linear or planar\n");
	  PlanarOrLinear = 0;
	  break;
	}
       
      }
      if(!PlanarOrLinear)
	break;
    }
    if(!PlanarOrLinear)
      break;
  }


  return PlanarOrLinear;
}


void Dimer::SetReferenceMonomerIndex(int index){

  reference_MonB = index;


  //commented out for debug purposes
  //replacing index B with reference_MonB in the symmetry list
  List_Sym_Dim.erase(List_Sym_Dim.begin()+1);
  List_Sym_Dim.push_back(index);

  //Setting MonB as the monomer outside the unit cell on the MonB_list
  //MonB_List.push_back(1);

  //Initialize k_list
  k_list.push_back(K_vector[0]);
  k_list.push_back(K_vector[1]);
  k_list.push_back(K_vector[2]);
  
  //printf("d(%i,%i)\n",indexA,indexB);
  //printf("K_vector = %i, %i, %i \n",K_vector[0],K_vector[1],K_vector[2]);
  //printf("k_list = %i, %i, %i \n",k_list[0],k_list[1],k_list[2]);
    
}

void Dimer::CreateQMJob(Monomer Monomers[], int NMon) {

  if ( Params::Parameters().GetQMType()==1 ) // QChem
    CreateQChemJob(Monomers, NMon);
  else if (Params::Parameters().GetQMType()==2 ){ //Molpro
    //molpro hessian found from finite difference
    if(Params::Parameters().DoFreq() ){
      if(Params::Parameters().DoEnergyFiniteDifferenceFreqs())
	CreateEnergyFiniteDifferenceHessianJobs(Monomers,NMon,false);
      else
	CreateFiniteDifferenceMolProJob(Monomers,NMon,false);
      if(Params::Parameters().DoCCSDTCorrection())
	CreateEnergyFiniteDifferenceHessianJobs(Monomers,NMon,true);
    }else{
      CreateMolProJob(Monomers,NMon);
      if(Params::Parameters().DoCCSDTCorrection()){
	CreateMolProCCSDTJob(Monomers,NMon);
        //Finite Difference Forces handled for CCSD(T)
        if(Params::Parameters().DoForces())
	  CreateFiniteDifferenceMolProJob(Monomers,NMon,true);
      }
    }
  }
  else if(Params::Parameters().GetQMType()==3) { // G09
    CreateG09Job(Monomers, NMon);
  }
  else if(Params::Parameters().GetQMType()==7) { // PSI4
    CreatePSI4Job(Monomers, NMon);
  }
  else{
    printf("ERROR:Dimer::CreateQMJob: Unknown QM_type: %d\n",
	   Params::Parameters().GetQMType() );
    exit(1);
  }
  



}
			
void Dimer::CreateQChemJob(Monomer Monomers[], int NMon, bool MM_job) {
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
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

  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  filename += ".in";

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateQChemJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  // Print comment section
  fprintf(job,"$comment\nDimer (%d,%d)\n$end\n\n",indexA,indexB);

  // Print $molecule section
  PrintQChemCartesian(job);
  
  // Print $rem section
  string rem;
  string basis;
  if (MM_job)
    rem = Params::Parameters().GetQChemRem2();
  else
    rem = Params::Parameters().GetQChemRem();

  basis = Params::Parameters().GetQChemBasis();

  fprintf(job,"%s\n",rem.c_str());  
  fprintf(job,"%s\n",basis.c_str());  

  if ( Params::Parameters().DoCounterpoise() ) { //JLM 05/13/16
    //adding ghost monomer A
    fprintf(job, "\n@@@\n");
    PrintQChemGhostCartesian(job,true,false);
    fprintf(job,"%s\n",rem.c_str());  
    fprintf(job,"%s\n",basis.c_str());

    //adding ghost monomer B
    fprintf(job, "\n@@@\n");
    PrintQChemGhostCartesian(job,false,true);
    fprintf(job,"%s\n",rem.c_str());  
    fprintf(job,"%s\n",basis.c_str());
    
  }
  
  // Optionally print $external_charges section  -- doesn't work with CP correction
  if ( Params::Parameters().UseEmbeddingCharges() && !Params::Parameters().DoCounterpoise()) {
    fprintf(job,"$external_charges\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB) {
	Monomers[i].PrintQChemEmbeddingCharges(job);
      }
    }
    fprintf(job,"$end\n");
  }

  fclose(job);
}

void Dimer::CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images) {

  string path;
  path = Params::Parameters().GetQMPath();
  string charge_file;
  char chrg_label[30];
  sprintf(chrg_label,"d%d.%d.ch", indexA, indexB);
  charge_file = path + "/" + chrg_label;
  string charge_file_local = chrg_label;
  
  FILE *chrg;
  if (( chrg = fopen(charge_file.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolProEmbeddedDipoleJob : cannot open file '%s'\n",charge_file.c_str());
    exit(1);
  }
  fprintf(chrg,"Embedding charges for dimer\n");

  //count the charges
  int number_charges = 0;
  for (int i=1;i<=NMon;i++) {
    if ( i != indexA && i != indexB) {
      if ( Monomers[i].GetUseInEmbedding() ) {
	number_charges += Monomers[i].GetNumberOfAtoms();
      }
    }
  }

  // Now count up the image charges:
  for (int i=1;i<=NMon_images;i++) {
    if ( i != indexA && i != indexB  ) {
      // Check to see if we are using a charge cut-off
      if ( MonomerImages[i].GetUseInEmbedding() ) {
	number_charges += MonomerImages[i].GetNumberOfAtoms();
      }
    }
  }

  fprintf(chrg,"%d\n", number_charges);

  for (int i=1;i<=NMon;i++) {
    if (i != indexA && i != indexB ) {
      if ( Monomers[i].GetUseInEmbedding() ) {
	Monomers[i].PrintMolProEmbeddingCharges(chrg);
      }
    }
  }
  
  // Now Add in the image charges:
  for (int i=1;i<=NMon_images;i++) {
    if ( i != indexA && i != indexB  ) {
      if ( MonomerImages[i].GetUseInEmbedding() ) {
	MonomerImages[i].PrintMolProEmbeddingCharges(chrg);
      } 
    } 
  }

  
  fclose(chrg);
}

void Dimer::CreateMolProChargeEmbeddingJob(Monomer Monomers[], int NMon) {
  string path;
  path = Params::Parameters().GetQMPath();
  string charge_file;
  char chrg_label[30];
  sprintf(chrg_label,"d%d.%d.ch", indexA, indexB);
  charge_file = path + "/" + chrg_label;
  string charge_file_local = chrg_label;

  FILE *chrg;
  if (( chrg = fopen(charge_file.c_str(),"w"))==NULL) {
    printf("Monomer::CreateMolProEmbeddedDipoleJob : cannot open file '%s'\n",charge_file.c_str());
    exit(1);
  }
  fprintf(chrg,"Embedding charges for dimer\n");

  //count the charges
  int number_charges = 0;
  for (int i=1;i<=NMon;i++) {
    if ( i != indexA && i != indexB) {
      number_charges += Monomers[i].GetNumberOfAtoms();
    }
  }
  
  fprintf(chrg,"%d\n", number_charges);

  for (int i=1;i<=NMon;i++) {
    if (i != indexA && i != indexB ) {
      Monomers[i].PrintMolProEmbeddingCharges(chrg);
    }
  }
  fclose(chrg);

}

void Dimer::CreateMolProJob(Monomer Monomers[], int NMon) {
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();
  if (Params::Parameters().DoFreq() )  {
      path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  filename += ".inp";

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateMolProJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  //input header
  
  fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
  //fprintf(job,"memory,500 m\n");
  fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
  fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
  fprintf(job,"nogprint;\n");
  fprintf(job,"gdirect;\n");
  fprintf(job,"symmetry,nosym;\n");
  fprintf(job,"orient,noorient;\n");
  fprintf(job,"angstrom\n");

  //Print Cartesian coords
  fprintf(job,"GEOMETRY={\n");
  MonA.PrintMolProMonomerCastesian(job,1);
  MonB.PrintMolProMonomerCastesian(job,MonA.GetNumberOfAtoms()+1);
  fprintf(job,"}\n\n");

  //print rem section
  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
 
  if(Params::Parameters().DoCounterpoise()){
    if(Params::Parameters().DoCBS()){
      
      //dimer section using the first basis set
      fprintf(job,"\n!dimer using first basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"SET,CHARGE=%i\n",charge);
      fprintf(job,"SET,SPIN=%i\n",spin-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      fprintf(job,"E_cc_HF1 = energy;\n");
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_cc_MP2_1 = energy;\n");
      
      //ghost monomer A using the first basis set
      fprintf(job,"\n!Ghost Monomer A using first basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccA_HF1 = energy;\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccA_MP2_1 = energy;\n");
      
      //monomer B using the first basis set
      fprintf(job,"\n!Ghost Monomer B using first basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonB.GetNumberOfAtoms();i++)
        fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccB_HF1 = energy;\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );  
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccB_MP2_1 = energy;\n");
      
      //counterpoise correction for first basis set
      fprintf(job,"\n!counterpoise correction for first basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"E_HF1 = (E_cc_HF1-E_ccA_HF1-E_ccB_HF1)\n");
      fprintf(job,"E_MP2_1 = (E_cc_MP2_1-E_ccA_MP2_1-E_ccB_MP2_1)\n\n");
      
      //Second basis set
      fprintf(job,"%s", Params::Parameters().GetMolProCBSRem().c_str() );
      
      //dimer section using the second basis set
      fprintf(job,"\n!dimer using second basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy\n");
      fprintf(job,"SET,CHARGE=%i\n",charge);
      fprintf(job,"SET,SPIN=%i\n",spin-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_cc_HF2 = energy;\n");
  
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_cc_MP2_2 = energy;\n");
      
      //monomer A using the second basis set
      fprintf(job,"\n!Ghost Monomer A using second basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccA_HF2 = energy;\n");
      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccA_MP2_2 = energy;\n");
      
      //monomer B using the first basis set
      fprintf(job,"\n!Ghost Monomer B using second basis\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonB.GetNumberOfAtoms();i++)
        fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
      fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccB_HF2 = energy;\n");

      fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      } 
      fprintf(job,"E_ccB_MP2_2 = energy;\n");
      
      //counterpoise correction for second basis set
      fprintf(job,"\n!counterpoise correction for second basis\n");
      fprintf(job,"!-------\n");
      fprintf(job,"E_HF2 = (E_cc_HF2-E_ccA_HF2-E_ccB_HF2)\n");
      fprintf(job,"E_MP2_2 = (E_cc_MP2_2-E_ccA_MP2_2-E_ccB_MP2_2)\n\n");
      
      //CBS extrapolation
      int basis1 = Params::Parameters().CBSBasis1();
      int basis2 = Params::Parameters().CBSBasis2();
      if(basis1 >= basis2){
	printf("ERROR:Dimer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	exit(0);
      }
      fprintf(job,"!CBS extrapolate\n");
      fprintf(job,"!-----\n");
      fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	      basis2,basis1);
      fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
	      basis2,basis1,basis2,basis1);
      fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
    }
    else{
      
      //dimer section
      fprintf(job,"!dimer\n");
      fprintf(job,"!-----\n");
      fprintf(job,"SET,CHARGE=%i\n",charge);
      fprintf(job,"SET,SPIN=%i\n",spin-1);
      fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_cc = energy\n\n");
      
      
      //monomer A
      fprintf(job,"!Ghost Monomer A\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
      fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
      if(Params::Parameters().DoForces()){
        fprintf(job,"FORCE\n");
      }
      
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccA = energy\n\n");
      
      //monomer B
      fprintf(job,"!Ghost Monomer B\n");
      fprintf(job,"!-----\n");
      fprintf(job,"dummy");
      for(int i=0;i<MonB.GetNumberOfAtoms();i++)
        fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
      fprintf(job,"\n");
      fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
      fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
      fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
      if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      fprintf(job,"E_ccB = energy\n\n");

      fprintf(job,"!summary\n");
      fprintf(job,"!-------\n");
      fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
      fprintf(job,"\n");
      
    }


  } 
  else{
    
    //Perform the CBS extrapolation
     if(Params::Parameters().DoCBS()){
       
       //first basis set, basis info already printed
       fprintf(job,"SET,CHARGE=%i\n",charge);
       fprintf(job,"SET,SPIN=%i\n",spin-1);
       fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
       fprintf(job,"SET,CHARGE=%i\n",charge);
       fprintf(job,"SET,SPIN=%i\n",spin-1);
       fprintf(job,"E_HF1 = energy\n");
       if(Params::Parameters().DoForces()){
	 fprintf(job,"FORCE\n");
       }
       if(Params::Parameters().DoFreq()){
	 fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
       }
       fprintf(job,"SET,CHARGE=%i\n",charge);
       fprintf(job,"SET,SPIN=%i\n",spin-1);
       fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
       fprintf(job,"E_MP2_1 = energy\n");
       if(Params::Parameters().DoForces()){
	 fprintf(job,"FORCE\n");
       }   
       if(Params::Parameters().DoFreq()){
	 fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
       }    
       fprintf(job,"\n");

       //second basis set
       fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
       fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
       fprintf(job,"E_HF2 = energy\n");
       if(Params::Parameters().DoForces()){
	 fprintf(job,"FORCE\n");
       }
       if(Params::Parameters().DoFreq()){
	 fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
       }
       fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
       fprintf(job,"E_MP2_2 = energy\n");
       if(Params::Parameters().DoForces()){
	 fprintf(job,"FORCE\n");
       }
       if(Params::Parameters().DoFreq()){
	 fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
       }
       fprintf(job,"\n");

       //Extrapole energy
       int basis1 = Params::Parameters().CBSBasis1();
       int basis2 = Params::Parameters().CBSBasis2();
       if(basis1 >= basis2){
	 printf("Dimer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	 exit(0);
       }
       fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	       basis2,basis1);
       fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
	       basis2,basis1,basis2,basis1);
       fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
     }
     else{
       fprintf(job,"SET,CHARGE=%i\n",charge);
       fprintf(job,"SET,SPIN=%i\n",spin-1);
       fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
       
       //print energy in output
       fprintf(job,"Edimer=energy\n");
       fprintf(job,"\n");
       
       if(Params::Parameters().DoForces()){
	 fprintf(job,"FORCE\n");
       }
       
       if(Params::Parameters().DoFreq()){
	 fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
       }
       
     }
   }
  fclose(job);
}

void Dimer::CreateMolProCCSDTJob(Monomer Monomers[0], int NMon){

  // Set up the filename, with the full path.  File is e.g. 'm1.in'
  string path = Params::Parameters().GetQMPath(); 
  if (Params::Parameters().DoFreq())
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
     path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  filename += ".CCSDT.inp";

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateMolCCSDTProJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  //input header
  fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
  //fprintf(job,"memory,500 m\n");
  fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
  fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
  fprintf(job,"nogprint\n");
  //fprintf(job,"gdirect\n");
  fprintf(job,"symmetry,nosym;\n");
  fprintf(job,"orient,noorient;\n");
  fprintf(job,"angstrom\n");

  //Print Cartesian coords
  fprintf(job,"GEOMETRY={\n");
  MonA.PrintMolProMonomerCastesian(job,1);
  MonB.PrintMolProMonomerCastesian(job,MonA.GetNumberOfAtoms()+1);
  fprintf(job,"}\n\n");


  //Forces calculations Created in CreateFiniteDifferenceMolProJob()
  if(Params::Parameters().DoCounterpoise()){

    //dimer section
    fprintf(job,"!dimer\n");
    fprintf(job,"!-----\n");
    //basis
    fprintf(job,"SET,CHARGE=%i\n",charge);
    fprintf(job,"SET,SPIN=%i\n",spin-1);
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
    //MP2 
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
    fprintf(job,"E_cc_MP2 = energy\n");
    /*
    if(Params::Parameters().DoForces()){
        fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
        fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    */
    //CCSDT
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
    fprintf(job,"E_cc_CCSDT = energy\n");
    /*
    if(Params::Parameters().DoForces()){
        fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
        fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    */

    //monomer A
    fprintf(job,"\n!Ghost Monomer A\n");
    fprintf(job,"!-----\n");
    fprintf(job,"dummy");
    for(int i=0;i<MonA.GetNumberOfAtoms();i++)
      fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
    fprintf(job,"\n");
    fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
    fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
    //MP2 
    fprintf(job,"%s", Params::Parameters().GetMolProCCSDTMP2().c_str() );
    fprintf(job,"E_ccA_MP2 = energy\n");

    /*
    if(Params::Parameters().DoForces()){
      fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
      fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
     }
    */

    //CCSDT
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
    fprintf(job,"E_ccA_CCSDT = energy\n");
    /*
    if(Params::Parameters().DoForces()){
        fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
        fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    */

    //monomer B
    fprintf(job,"\n!Ghost Monomer B\n");
    fprintf(job,"!-----\n");
    fprintf(job,"dummy");
    for(int i=0;i<MonB.GetNumberOfAtoms();i++)
      fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
    fprintf(job,"\n");
    fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
    fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
    fprintf(job,"E_ccB_MP2 = energy\n"); 
    /*
    if(Params::Parameters().DoForces()){
	fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
	fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
    */
      //CCSDT
      fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
      fprintf(job,"E_ccB_CCSDT = energy\n");
      /*
      if(Params::Parameters().DoForces()){
          fprintf(job,"FORCE\n");
      }
      if(Params::Parameters().DoFreq()){
          fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
      }
      */
      
      fprintf(job,"\nE_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
      fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n"); 
  }
  else{
    //Charge and spin
    fprintf(job,"SET,CHARGE=%i\n",charge);
    fprintf(job,"SET,SPIN=%i\n",spin-1);
    //basis
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
    //MP2 
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
    fprintf(job,"E_MP2 = energy\n");
    if(Params::Parameters().DoForces()){
        fprintf(job,"FORCE\n");
    }
    if(Params::Parameters().DoFreq()){
        fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
    fprintf(job,"\n");
    //CCSDT
    fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
    fprintf(job,"E_CCSDT = energy\n");
    if(Params::Parameters().DoForces()){
      fprintf(job,"{FORCE,NUMERICAL,RSTEP = 0.001889726}\n");
    }
    if(Params::Parameters().DoFreq()){
        fprintf(job,"{frequencies, step=0.001889726; noproject; print,hessian; coord,3n;}\n");
    }
  }
  fprintf(job,"\n");
  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
  fclose(job);
}

void Dimer::CreateFiniteDifferenceMolProJob(Monomer Monomers[], int NMon,Vector Coord, bool CCSDT){

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  string cmd;

  
  string pathstring = Params::Parameters().GetQMPath();;
  if(Params::Parameters().DoFreq())
    pathstring = Params::Parameters().GetHessianQMPath();
  

   //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  /*
  if(!LatticeParams){
    char num[10];
     sprintf(num,"%d",Coord);
     pathstring += "+" + num;
  }else{ 
     if(Coord == 0) pathstring += "+a";
     else if(Coord == 1) pathstring += "+b";
     else if(Coord == 2) pathstring += "+c";
     else if(Coord == 3) pathstring += "+alpha;
     else if(Coord == 4) pathstring += "+beta";
     else pathstring += "+gamma";
  }
  */

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = pathstring + "/" + dimername;
  //printf("dimername = %s\n dimerpath = %s\n",dimername.c_str(),dimerpath.c_str());

  if(CCSDT)
    dimerpath += ".CCSDT";

  printf("path = %s\n",pathstring.c_str());
  exit(0);
      
 
  //create dimer director if doesn't already exist
  struct stat st;
  if( stat(pathstring.c_str(),&st) != 0){
    //printf("Directory %s not found.  Creating it.\n",
    //   dimerpath.c_str());
    cmd = "mkdir " + pathstring;	 
    //printf("%s\n",cmd.c_str());
    system(cmd.c_str());
  }  
  
  //List of Symmetrical equivalent atoms in the asymetrical cell previous shifted.
  //Used to make sure that each symmetrical shift only happens once.
  vector<int> SymList;
  for(int i=0;i<Na+Nb;i++){
    
    int SymAtom;
    if(i<Na){
      SymAtom = MonA.GetAtom(i).GetSymmetricalAtom();
    }else{
      int j = i - Na;
      SymAtom = MonA.GetAtom(j).GetSymmetricalAtom();
    }
    
   
    
  }

  
}


//Use if gradients or hessians are found with finite difference.
//Must be used for counterpoise corrected numerical gradients and hessians (i.e. CCSD(T) gradients and all counterpoise correction hessians)
void Dimer::CreateFiniteDifferenceMolProJob(Monomer Monomers[], int NMon,bool CCSDT) {

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  string cmd;


  string pathstring = Params::Parameters().GetQMPath();;
  if(Params::Parameters().DoFreq())
    pathstring = Params::Parameters().GetHessianQMPath();

   //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();



  
  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = pathstring + "/" + dimername;

  if(CCSDT)
    dimerpath += ".CCSDT";



  //create dimer director if doesn't already exist
  struct stat st;
  if( stat(dimerpath.c_str(),&st) != 0){
    //printf("Directory %s not found.  Creating it.\n",
    //   dimerpath.c_str());
    cmd = "mkdir " + dimerpath;	 
    //printf("%s\n",cmd.c_str());
    system(cmd.c_str());
  }

  //else {
  //  printf("Directory %s exists.\n",
  //   dimerpath.c_str());
  //}
  
  //cmd = "cd " + dimerpath;
  //printf("%s\n",cmd.c_str());
  //system(cmd.c_str());
  
  FILE *job;
  string filename;
  for(int i=0;i<3*(Na+Nb);i++){
    // The plus input file; 
    char num[10];
    sprintf(num,"%d",i);
    filename = dimerpath + "/" + dimername + "+" + num + ".inp";
    // Open the input file for writing
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateFiniteDifferenceMolProJob() : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    //input header
    fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
    fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
    fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
    fprintf(job,"symmetry,nosym;\n");
    fprintf(job,"orient,noorient;\n");
    fprintf(job,"angstrom\n");


  
    //Set Coordinates of file
    fprintf(job,"GEOMETRY={\n"); 
    //MonA
    Vector Coords = MonA.GetCoordinates();
    if(i < 3*Na)
      Coords[i] -= 0.001;
    MonA.PrintMolProMonomerCastesian(job,Coords,1);
    //MonB
    Coords = MonB.GetCoordinates();
    if(i >= 3*Na)
      Coords[i-3*Na] -= 0.001;
    MonB.PrintMolProMonomerCastesian(job,Coords,Na+1);
    fprintf(job,"}\n");
    
    if(Params::Parameters().DoCounterpoise()){

      //Can only handle finite difference gradients
      if(CCSDT){

	//dimer section
	fprintf(job,"!dimer\n");
	fprintf(job,"!-----\n");
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",charge);
	fprintf(job,"SET,SPIN=%i\n",spin-1);
	//basis
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	//MP2 
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	fprintf(job,"E_cc_MP2 = energy\n");
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_cc_CCSDT = energy\n");
	//monomer A
	fprintf(job,"\n!Ghost Monomer A\n");
	fprintf(job,"!-----\n");
	fprintf(job,"dummy");
	for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	  fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	fprintf(job,"\n");
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	//MP2 
	fprintf(job,"%s", Params::Parameters().GetMolProCCSDTMP2().c_str() );
	fprintf(job,"E_ccA_MP2 = energy\n");
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_ccA_CCSDT = energy\n");
	//monomer B
	fprintf(job,"\n!Ghost Monomer B\n");
	fprintf(job,"!-----\n");
	fprintf(job,"dummy");
	for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	  fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	fprintf(job,"\n");
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	fprintf(job,"E_ccB_MP2 = energy\n"); 
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_ccB_CCSDT = energy\n");
	
	fprintf(job,"\nE_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	fprintf(job,"Edimer = E_CCSDT - E_MP2\n");

      }
      else{

	//CBS limit calculations
	if(Params::Parameters().DoCBS()){
	  
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  
	  //dimer section using the first basis set
	  fprintf(job,"\n!dimer using first basis\n");
	  fprintf(job,"!-----\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_cc_HF1 = energy;\n");
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_MP2_1 = energy;\n");
	  
	  //monomer A using the first basis set
	  fprintf(job,"\n!Ghost Monomer A using first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_HF1 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_MP2_1 = energy;\n");
	  
	  //monomer B using the first basis set
	  fprintf(job,"\n!Ghost Monomer B using first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_HF1 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );  
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_MP2_1 = energy;\n");
	  
	  //counterpoise correction for first basis set
	  fprintf(job,"\n!counterpoise correction for first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF1 = (E_cc_HF1-E_ccA_HF1-E_ccB_HF1)\n");
	  fprintf(job,"E_MP2_1 = (E_cc_MP2_1-E_ccA_MP2_1-E_ccB_MP2_1)\n\n");
	  
	  //Second basis set
	  fprintf(job,"%s", Params::Parameters().GetMolProCBSRem().c_str() );
	  
	  //dimer section using the second basis set
	  fprintf(job,"\n!dimer using second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_HF2 = energy;\n");
	  
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_MP2_2 = energy;\n");
	  
	  //monomer A using the second basis set
	  fprintf(job,"\n!Ghost Monomer A using second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_HF2 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_MP2_2 = energy;\n");
	  
	  //monomer B using the first basis set
	  fprintf(job,"\n!Ghost Monomer B using second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_HF2 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  } 
	  fprintf(job,"E_ccB_MP2_2 = energy;\n");
	  
	  //counterpoise correction for second basis set
	  fprintf(job,"\n!counterpoise correction for second basis\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"E_HF2 = (E_cc_HF2-E_ccA_HF2-E_ccB_HF2)\n");
	  fprintf(job,"E_MP2_2 = (E_cc_MP2_2-E_ccA_MP2_2-E_ccB_MP2_2)\n\n");
	  
	  //CBS extrapolation
	  int basis1 = Params::Parameters().CBSBasis1();
	  int basis2 = Params::Parameters().CBSBasis2();
	  if(basis1 >= basis2){
	    printf("ERROR:Dimer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	    exit(0);
	  }
	  fprintf(job,"!CBS extrapolate\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		  basis2,basis1);
	  fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		  basis2,basis1,basis2,basis1);
	  fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
	  
	  
	}
	else{
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  
	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc = energy\n\n");
	  
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA = energy\n\n");
	  
	  //monomer B
	  fprintf(job,"!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	}
      }
    } 
    else{
      //CBS limit calculations
      if(Params::Parameters().DoCBS()){
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",charge);
	fprintf(job,"SET,SPIN=%i\n",spin-1);
	
	//first basis set
	fprintf(job,"!First basis\n");
	fprintf(job,"!-----\n");
	fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_HF1 = energy\n");
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_MP2_1 = energy\n\n");
      
	//second basis set
	fprintf(job,"!Second basis\n");
	fprintf(job,"!-----\n");
	fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
	fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_HF2 = energy\n");
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_MP2_2 = energy\n\n");
	//CBS Extrapole
	int basis1 = Params::Parameters().CBSBasis1();
	int basis2 = Params::Parameters().CBSBasis2();
	if(basis1 >= basis2){
	  printf("Dimer::CreateFiniteDifferenceMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	  exit(0);
	}
	fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		basis2,basis1);
	fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		basis2,basis1,basis2,basis1);
	fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
      }
      else{
	//print rem section
	fprintf(job,"%s", Params::Parameters().GetMolProRem().c_str() );
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );    
	
	//print energy in output
	fprintf(job,"Edimer=energy\n");
	fprintf(job,"\n");
	
	if(Params::Parameters().DoFreq())
	  fprintf(job,"FORCE\n");
      }
    }
    fclose(job);
    
    // The minus input file;
    filename = dimerpath + "/" + dimername + "-"  + num + ".inp";
    
    // Open the input file for writing
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateMolProJob : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }
    //input header
    fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
    fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
    fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
    fprintf(job,"symmetry,nosym;\n");
    fprintf(job,"orient,noorient;\n");
    fprintf(job,"angstrom\n");

    //Set Coordinates of file
    fprintf(job,"GEOMETRY={\n"); 
    //MonA
    Coords = MonA.GetCoordinates();
    if(i < 3*Na)
      Coords[i] += 0.001;
    MonA.PrintMolProMonomerCastesian(job,Coords,1);
    //MonB
    Coords = MonB.GetCoordinates();
    if(i >= 3*Na)
      Coords[i-3*Na] += 0.001;
    MonB.PrintMolProMonomerCastesian(job,Coords,Na+1);
    fprintf(job,"}\n");
    
    if(Params::Parameters().DoCounterpoise()){


      //Can only handle finite difference gradients
      if(CCSDT){

	//dimer section
	fprintf(job,"!dimer\n");
	fprintf(job,"!-----\n");

	//basis
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	//MP2 	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",charge);
	fprintf(job,"SET,SPIN=%i\n",spin-1);
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	fprintf(job,"E_cc_MP2 = energy\n");
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_cc_CCSDT = energy\n");
	//monomer A
	fprintf(job,"\n!Ghost Monomer A\n");
	fprintf(job,"!-----\n");
	fprintf(job,"dummy");
	for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	  fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	fprintf(job,"\n");
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	//MP2 
	fprintf(job,"%s", Params::Parameters().GetMolProCCSDTMP2().c_str() );
	fprintf(job,"E_ccA_MP2 = energy\n");
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_ccA_CCSDT = energy\n");
	//monomer B
	fprintf(job,"\n!Ghost Monomer B\n");
	fprintf(job,"!-----\n");
	fprintf(job,"dummy");
	for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	  fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	fprintf(job,"\n");
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	fprintf(job,"E_ccB_MP2 = energy\n"); 
	//CCSDT
	fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	fprintf(job,"E_ccB_CCSDT = energy\n");
	
	fprintf(job,"\nE_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	fprintf(job,"Edimer = E_CCSDT - E_MP2\n");

      }
      else{
	//CBS limit calculations
	if(Params::Parameters().DoCBS()){
	  
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  
	  //dimer section using the first basis set
	  fprintf(job,"\n!dimer using first basis\n");
	  fprintf(job,"!-----\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_cc_HF1 = energy;\n");
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_MP2_1 = energy;\n");
	  
	  //monomer A using the first basis set
	  fprintf(job,"\n!Ghost Monomer A using first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_HF1 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_MP2_1 = energy;\n");
	  
	  //monomer B using the first basis set
	  fprintf(job,"\n!Ghost Monomer B using first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_HF1 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );  
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_MP2_1 = energy;\n");
	  
	  //counterpoise correction for first basis set
	  fprintf(job,"\n!counterpoise correction for first basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF1 = (E_cc_HF1-E_ccA_HF1-E_ccB_HF1)\n");
	  fprintf(job,"E_MP2_1 = (E_cc_MP2_1-E_ccA_MP2_1-E_ccB_MP2_1)\n\n");
	  
	  //Second basis set
	  fprintf(job,"%s", Params::Parameters().GetMolProCBSRem().c_str() );
	  
	  //dimer section using the second basis set
	  fprintf(job,"\n!dimer using second basis\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy\n");
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_HF2 = energy;\n");
	  
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc_MP2_2 = energy;\n");
	
	  //monomer A using the second basis set
	  fprintf(job,"\n!Ghost Monomer A using second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_HF2 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA_MP2_2 = energy;\n");
	  
	  //monomer B using the first basis set
	  fprintf(job,"\n!Ghost Monomer B using second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB_HF2 = energy;\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  } 
	  fprintf(job,"E_ccB_MP2_2 = energy;\n");
	  
	  //counterpoise correction for second basis set
	  fprintf(job,"\n!counterpoise correction for second basis\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"E_HF2 = (E_cc_HF2-E_ccA_HF2-E_ccB_HF2)\n");
	  fprintf(job,"E_MP2_2 = (E_cc_MP2_2-E_ccA_MP2_2-E_ccB_MP2_2)\n\n");
	  
	  //CBS extrapolation
	  int basis1 = Params::Parameters().CBSBasis1();
	  int basis2 = Params::Parameters().CBSBasis2();
	  if(basis1 >= basis2){
	    printf("ERROR:Dimer::CreateMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	    exit(0);
	  }
	  fprintf(job,"!CBS extrapolate\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		  basis2,basis1);
	  fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		  basis2,basis1,basis2,basis1);
	  fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
	}
	else{
	  
	  //print rem section
	  fprintf(job,"%s", Params::Parameters().GetMolProRem().c_str() );
	  
	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n",spin-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  
	  if(Params::Parameters().DoFreq()){
	fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_cc = energy\n\n");
      
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccA = energy\n\n");
	  
	  //monomer B
	  fprintf(job,"!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //change and spin (spin = multiplicity -1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  if(Params::Parameters().DoFreq()){
	    fprintf(job,"FORCE\n");
	  }
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	}
      }
    } 
    else{

      //CBS limit calculations
      if(Params::Parameters().DoCBS()){
	//change and spin (spin = multiplicity -1)
	fprintf(job,"SET,CHARGE=%i\n",charge);
	fprintf(job,"SET,SPIN=%i\n",spin-1);
	//first basis set
	fprintf(job,"!First basis\n");
	fprintf(job,"!-----\n");
	fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_HF1 = energy\n");
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_MP2_1 = energy\n\n");
      
	//second basis set
	fprintf(job,"!Second basis\n");
	fprintf(job,"!-----\n");
	fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
	fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_HF2 = energy\n");
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	}
	fprintf(job,"E_MP2_2 = energy\n\n");
	//CBS Extrapole
	int basis1 = Params::Parameters().CBSBasis1();
	int basis2 = Params::Parameters().CBSBasis2();
	if(basis1 >= basis2){
	  printf("Dimer::CreateFiniteDifferenceMolProJob(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	  exit(0);
	}
	fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		basis2,basis1);
	fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		basis2,basis1,basis2,basis1);
	fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
      }
      else{
	//print rem section
	fprintf(job,"%s", Params::Parameters().GetMolProRem().c_str() );
	fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() );  
	
	//print energy in output
	fprintf(job,"Edimer=energy\n");
	fprintf(job,"\n");

	if(Params::Parameters().DoFreq()){
	  fprintf(job,"FORCE\n");
	} 
      }
    }
    fclose(job); 
  }

  /*
  //creating csh script to run all jobs
  filename = dimerpath + "/" + dimername + ".csh";

  // Open the input file for writing
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateMolProJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  fprintf(job,"#!/bin/csh\n\n");
  fprintf(job,"rm -f *out *xml\n");
  fprintf(job,"foreach s(d*inp)\n");
  fprintf(job,"  molpro $s\n");
  fprintf(job,"end\n");
  fclose(job);
  */

}

//Create CCSD(T) correction finite difference calculations in molpro
void Dimer::CreateEnergyFiniteDifferenceHessianJobs(Monomer Monomers[], int NMon,bool CCSDT) {
   
 int Na = MonA.GetNumberOfAtoms();
 int Nb = MonB.GetNumberOfAtoms();
 Vector CoordA = MonA.GetCoordinates();
 Vector CoordB = MonB.GetCoordinates();
 string cmd;
 

  string pathstring = Params::Parameters().GetQMPath();;
  if(Params::Parameters().DoFreq())
    pathstring = Params::Parameters().GetHessianQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    pathstring += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = pathstring + "/" + dimername;
  //printf("dimername = %s\n dimerpath = %s\n",dimername.c_str(),dimerpath.c_str());
  
  if(CCSDT)
    dimerpath += ".CCSDT";

  //create dimer director if doesn't already exist
  struct stat st;
  if( stat(dimerpath.c_str(),&st) != 0){
    //printf("Directory %s not found.  Creating it.\n",
    //   dimerpath.c_str());
    cmd = "mkdir " + dimerpath;	 
    //printf("%s\n",cmd.c_str());
    system(cmd.c_str());
  }



  for(int i=0;i<3*(Na+Nb);i++){
    for(int j=0;j<=i;j++){

      char num[10];
      sprintf(num,"%d",i); 
      char num2[10];
      sprintf(num2,"%d",j);    
      
      //Second derivitive plus plus step
      string filename = dimerpath + "/" + dimername + "+" + num + "+" + num2 + ".inp";
      //printf("filename = %s\n",filename.c_str());
      
      // Open the input file for writing
      FILE *job;
      if ((job = fopen(filename.c_str(),"w"))==NULL) {
	printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs() : Cannot open file '%s'\n",filename.c_str());
	exit(1);
      }
      
      //input header
      fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
      fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
      fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
      fprintf(job,"symmetry,nosym;\n");
      fprintf(job,"orient,noorient;\n");
      fprintf(job,"angstrom\n");
      
      //Set Coordinates of file
      fprintf(job,"GEOMETRY={\n"); 
      //MonA
      Vector Coord = CoordA;
      if(i<3*Na)
	Coord[i] += 0.001;
      if(j<3*Na)
	Coord[j] += 0.001;
      MonA.PrintMolProMonomerCastesian(job,Coord,1);
      
      //MonB
      Coord = CoordB;
      if(i >= 3*Na)
	Coord[i-3*Na] += 0.001;
      if(j >= 3*Na)
	Coord[j-3*Na] += 0.001;
      MonB.PrintMolProMonomerCastesian(job,Coord,Na+1);
      fprintf(job,"}\n\n");      
      
      if(Params::Parameters().DoCounterpoise()){
	if(CCSDT){
	  
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());

	  //Dimer job
	  fprintf(job,"\n!dimer\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_cc_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_cc_CCSDT = energy\n");
	  fprintf(job,"\n");  

	  //monomer A
	  fprintf(job,"\n!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);


	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccA_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccA_CCSDT = energy\n");
	  fprintf(job,"\n");

	  //monomer B
	  fprintf(job,"\n!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);

	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccB_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccB_CCSDT = energy\n");
	  fprintf(job,"\n\n");

	  fprintf(job,"E_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	  fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");

	}
	//CBS limit calculations :: Not tested
	else if(Params::Parameters().DoCBS()){
	  printf("ERROR::Dimer::CreateEnergyFiniteDifferenceHessianJobs() This function has not been implemented for CBS limit extrapolation\n");
	  exit(0);
	  
	}else{
	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );

	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_cc = energy\n\n");
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);

	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_ccA = energy\n\n");
	  //monomer B
	  fprintf(job,"!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	  
	}
      }
      else{
	//CCSD(T) correction
	if(CCSDT){
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_CCSDT = energy\n");
	  fprintf(job,"\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	  
	  //CBS limit calculations :: Not tested
	}else if(Params::Parameters().DoCBS()){
	  printf("ERROR::Dimer::CreateEnergyFiniteDifferenceHessianJobs() This function has not been implemented for CBS limit extrapolation\n");
	  exit(0);	  	  
     	}
	//No CBS CCSD(T) job:: Not tested
	else{

	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  
	  //print energy in output
	  fprintf(job,"Edimer=energy\n");
	  fprintf(job,"\n");
	  
	}
      }
      fclose(job);
      
      //Second derivitive plus minus step
      filename = dimerpath + "/" + dimername + "+" + num + "-" + num2 + ".inp";
      //printf("filename = %s\n",filename.c_str());
      // Open the input file for writing
      if ((job = fopen(filename.c_str(),"w"))==NULL) {
	printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs() : Cannot open file '%s'\n",filename.c_str());
	exit(1);
      }      
      
      //input header
      fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
      fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
      fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
      fprintf(job,"symmetry,nosym;\n");
      fprintf(job,"orient,noorient;\n");
      fprintf(job,"angstrom\n");
      
      //Set Coordinates of file
      fprintf(job,"GEOMETRY={\n"); 
      //MonA
      Coord = CoordA;
      if(i<3*Na)
	Coord[i] += 0.001;
      if(j<3*Na)
	Coord[j] -= 0.001;
      MonA.PrintMolProMonomerCastesian(job,Coord,1);
      
      //MonB
      Coord = CoordB;
      if(i >= 3*Na)
	Coord[i-3*Na] += 0.001;
      if(j >= 3*Na)
	Coord[j-3*Na] -= 0.001;
      MonB.PrintMolProMonomerCastesian(job,Coord,Na+1);
      fprintf(job,"}\n\n"); 
      
      if(Params::Parameters().DoCounterpoise()){
	if(CCSDT){
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());

	  //Dimer job
	  fprintf(job,"\n!dimer\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_cc_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_cc_CCSDT = energy\n");
	  fprintf(job,"\n");  

	  //monomer A
	  fprintf(job,"\n!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccA_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccA_CCSDT = energy\n");
	  fprintf(job,"\n");

	  //monomer B
	  fprintf(job,"\n!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccB_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccB_CCSDT = energy\n");
	  fprintf(job,"\n\n");

	  fprintf(job,"E_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	  fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");

	}
	//CBS limit calculations :: Not tested
	else if(Params::Parameters().DoCBS()){
	  printf("ERROR::Dimer::CreateEnergyFiniteDifferenceHessianJobs() This function has not been implemented for CBS limit extrapolation\n");
	  exit(0);
	  
	}else{
	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );

	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_cc = energy\n\n");
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_ccA = energy\n\n");
	  //monomer B
	  fprintf(job,"!Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	  
	}

      }else{

        //CCSD(T) correction
        if(CCSDT){

	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_CCSDT = energy\n");
	  fprintf(job,"\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	
   	  //CBS limit calculations :: Not tested
        }else if(Params::Parameters().DoCBS()){
	
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //first basis set
	  fprintf(job,"!First basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF1 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_1 = energy\n\n");
	
	  //second basis set
	  fprintf(job,"!Second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF2 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_2 = energy\n\n");
	
	  //CBS Extrapole
	  int basis1 = Params::Parameters().CBSBasis1();
	  int basis2 = Params::Parameters().CBSBasis2();
	  if(basis1 >= basis2){
	    printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	    exit(0);
	  }
	  fprintf(job,"!CBS extrapolate\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
	 	  basis2,basis1);
	  fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		  basis2,basis1,basis2,basis1);
	  fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
	
	  //Not CBS CCSD(T) job:: Not tested
        }else{

	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);

	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	
	  //print energy in output
	  fprintf(job,"Edimer=energy\n");
	  fprintf(job,"\n");
	  //fprintf(job,"!summary\n");
	  //fprintf(job,"!--------\n");
	  //fprintf(job,"show,Emonomer\n");
        }
      }
      fclose(job);
      
      
      //Second derivitive minus plus step
      filename = dimerpath + "/" + dimername + "-" + num + "+" + num2 + ".inp";
      //printf("filename = %s\n",filename.c_str());
      // Open the input file for writing
      if ((job = fopen(filename.c_str(),"w"))==NULL) {
	printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs() : Cannot open file '%s'\n",filename.c_str());
	exit(1);
      }  
      
      //input header
      fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
      fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
      fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
      fprintf(job,"symmetry,nosym;\n");
      fprintf(job,"orient,noorient;\n");
      fprintf(job,"angstrom\n");
      
      //Set Coordinates of file
      fprintf(job,"GEOMETRY={\n"); 
      //MonA
      Coord = CoordA;
      if(i<3*Na)
	Coord[i] -= 0.001;
      if(j<3*Na)
	Coord[j] += 0.001;
      MonA.PrintMolProMonomerCastesian(job,Coord,1);
      
      //MonB
      Coord = CoordB;
      if(i >= 3*Na)
	Coord[i-3*Na] -= 0.001;
      if(j >= 3*Na)
	Coord[j-3*Na] += 0.001;
      MonB.PrintMolProMonomerCastesian(job,Coord,Na+1);
      fprintf(job,"}\n\n"); 
      
      if(Params::Parameters().DoCounterpoise()){
	if(CCSDT){
	  
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());

	  //Dimer job
	  fprintf(job,"\n!dimer\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_cc_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_cc_CCSDT = energy\n");
	  fprintf(job,"\n");  

	  //monomer A
	  fprintf(job,"\n!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  //charge and spin (spin = multiplicity - 1)
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccA_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccA_CCSDT = energy\n");
	  fprintf(job,"\n");

	  //monomer B
	  fprintf(job,"\n!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccB_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccB_CCSDT = energy\n");
	  fprintf(job,"\n\n");

	  fprintf(job,"E_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	  fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	}
	//CBS limit calculations :: Not tested
	else if(Params::Parameters().DoCBS()){
	  printf("ERROR::Dimer::CreateEnergyFiniteDifferenceHessianJobs() This function has not been implemented for CBS limit extrapolation\n");
	  exit(0);
	  
	}else{
      	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_cc = energy\n\n");
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_ccA = energy\n\n");
	  //monomer B
	  fprintf(job,"!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	  
	}
      }else{
        //CCSD(T) correction
        if(CCSDT){
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_CCSDT = energy\n");
	  fprintf(job,"\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	
	  //CBS limit calculations :: Not tested
         }else if(Params::Parameters().DoCBS()){
	
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //first basis set
	  fprintf(job,"!First basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF1 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_1 = energy\n\n");
	
	  //second basis set
	  fprintf(job,"!Second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF2 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_2 = energy\n\n");
	
	  //CBS Extrapole
	  int basis1 = Params::Parameters().CBSBasis1();
	  int basis2 = Params::Parameters().CBSBasis2();
	  if(basis1 >= basis2){
	    printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	    exit(0);
	  }
	  fprintf(job,"!CBS extrapolate\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		  basis2,basis1);
	  fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		  basis2,basis1,basis2,basis1);
	  fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
	
	  //Not CBS CCSD(T) job:: Not tested
        }else{

	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	
	  //print energy in output
	  fprintf(job,"Edimer=energy\n");
	  fprintf(job,"\n");
        }
      }
      fclose(job);
      
      
      //Second derivitive minus plus step
      filename = dimerpath + "/" + dimername + "-" + num + "-" + num2 + ".inp";
      //printf("filename = %s\n",filename.c_str());
      // Open the input file for writing
      if ((job = fopen(filename.c_str(),"w"))==NULL) {
	printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs() : Cannot open file '%s'\n",filename.c_str());
	exit(1);
      }  
      
      //input header
      fprintf(job,"*** Dimer (%d,%d) \n",indexA,indexB);
      fprintf(job,"memory,%i m\n",Params::Parameters().MemoryUse());
      fprintf(job,"gthresh,orbital=1d-8,energy=1d-8,grid=1d-8\n");
      fprintf(job,"symmetry,nosym;\n");
      fprintf(job,"orient,noorient;\n");
      fprintf(job,"angstrom\n");
      
      
      //Set Coordinates of file
      fprintf(job,"GEOMETRY={\n"); 
      //MonA
      Coord = CoordA;
      if(i<3*Na)
	Coord[i] -= 0.001;
      if(j<3*Na)
	Coord[j] -= 0.001;
      MonA.PrintMolProMonomerCastesian(job,Coord,1);
      
      //MonB
      Coord = CoordB;
      if(i >= 3*Na)
	Coord[i-3*Na] -= 0.001;
      if(j >= 3*Na)
	Coord[j-3*Na] -= 0.001;
      MonB.PrintMolProMonomerCastesian(job,Coord,Na+1);
      fprintf(job,"}\n\n"); 
      
      if(Params::Parameters().DoCounterpoise()){
	if(CCSDT){
	  
	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());

	  //Dimer job
	  fprintf(job,"\n!dimer\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_cc_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_cc_CCSDT = energy\n");
	  fprintf(job,"\n");  

	  //monomer A
	  fprintf(job,"\n!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccA_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccA_CCSDT = energy\n");
	  fprintf(job,"\n");

	  //monomer B
	  fprintf(job,"\n!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_ccB_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_ccB_CCSDT = energy\n");
	  fprintf(job,"\n\n");

	  fprintf(job,"E_MP2 = E_cc_MP2 - E_ccA_MP2 - E_ccB_MP2\n");
	  fprintf(job,"E_CCSDT = E_cc_CCSDT - E_ccA_CCSDT - E_ccB_CCSDT\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	//CBS limit calculations :: Not tested
	}else if(Params::Parameters().DoCBS()){
	  printf("ERROR::Dimer::CreateEnergyFiniteDifferenceHessianJobs() This function has not been implemented for CBS limit extrapolation\n");
	  exit(0);
	  
	}else{
	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );      

	  //dimer section
	  fprintf(job,"!dimer\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_cc = energy\n\n");
	  //monomer A
	  fprintf(job,"!Ghost Monomer A\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonA.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonA.GetAtom(i).GetSymbol().c_str(),i+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonB.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonB.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	  fprintf(job,"E_ccA = energy\n\n");
	  //monomer B
	  fprintf(job,"!Ghost Monomer B\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"dummy");
	  for(int i=0;i<MonB.GetNumberOfAtoms();i++)
	    fprintf(job,",%s%i",MonB.GetAtom(i).GetSymbol().c_str(),i+MonA.GetNumberOfAtoms()+1);
	  fprintf(job,"\n");
	  fprintf(job,"SET,CHARGE=%i\n",MonA.GetChargeState());
	  fprintf(job,"SET,SPIN=%i\n\n",MonA.GetSpinState()-1);
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );  
	  fprintf(job,"E_ccB = energy\n\n");
	  
	  fprintf(job,"!summary\n");
	  fprintf(job,"!-------\n");
	  fprintf(job,"Edimer= (E_cc-E_ccA-E_ccB)\n");
	  fprintf(job,"\n");
	  
	}
      }else{
        //CCSD(T) correction
        if(CCSDT){

	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //basis
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTBasis().c_str());
	  //MP2 contribution
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTMP2().c_str());
	  fprintf(job,"E_MP2 = energy\n");
	  //CCSDT
	  fprintf(job,"%s",Params::Parameters().GetMolProCCSDTInst().c_str());
	  fprintf(job,"E_CCSDT = energy\n");
	  fprintf(job,"\n");
	  fprintf(job,"Edimer = E_CCSDT - E_MP2\n");
	
	  //CBS limit calculations :: Not tested
        }else if(Params::Parameters().DoCBS()){
	 
	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);	  

	  //first basis set
	  fprintf(job,"!First basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF1 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_1 = energy\n\n");
	 
	  //second basis set
	  fprintf(job,"!Second basis\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"%s\n", Params::Parameters().GetMolProCBSRem().c_str() );
	  fprintf(job,"%s",Params::Parameters().GetMolProHFRem().c_str() );
	  fprintf(job,"E_HF2 = energy\n");
	  fprintf(job,"%s", Params::Parameters().GetMolProInst().c_str() ); 
	  fprintf(job,"E_MP2_2 = energy\n\n");
	
	  //CBS Extrapole
	  int basis1 = Params::Parameters().CBSBasis1();
	  int basis2 = Params::Parameters().CBSBasis2();
	  if(basis1 >= basis2){
	    printf("Dimer::CreateEnergyFiniteDifferenceHessianJobs(): basis1 cannot be a larger basis2 during CBS extrapolation\n");
	    exit(0);
	  }
	  fprintf(job,"!CBS extrapolate\n");
	  fprintf(job,"!-----\n");
	  fprintf(job,"E_HF_inf = E_HF2 + (E_HF2 - E_HF1)/(2.71828^(1.54*(%i-%i)) - 1)\n",
		  basis2,basis1);
	  fprintf(job,"E_Coor_inf = (%i^3*(E_MP2_2-E_HF2)-%i^3*(E_MP2_1-E_HF1))/(%i^3-%i^3)\n",
		  basis2,basis1,basis2,basis1);
	  fprintf(job,"Edimer = E_HF_inf + E_Coor_inf\n");
	
	  //No CBS CCSD(T) job:: Not tested
        }else{

	  fprintf(job,"SET,CHARGE=%i\n",charge);
	  fprintf(job,"SET,SPIN=%i\n\n",spin-1);

	  //print rem section
	  fprintf(job,"%s\n", Params::Parameters().GetMolProRem().c_str() );
	  fprintf(job,"%s\n", Params::Parameters().GetMolProInst().c_str() );
	
	  //print energy in output
	  fprintf(job,"Edimer=energy\n");
	  fprintf(job,"\n");
	  //fprintf(job,"!summary\n");
	  //fprintf(job,"!--------\n");
	  //fprintf(job,"show,Emonomer\n");
        }
      }
      fclose(job);
    }
  }
}



// Wrapper for creating MM jobs
void Dimer::CreateMMJob(Monomer Monomers[], int NMon) {
  if ( Params::Parameters().GetMMType()==1 ) { // Tinker
    CreateTinkerJob(Monomers, NMon);
  }
  else if ( Params::Parameters().GetMMType()==2 ) { // AIFF
    // Do nothing.  Dimers get handled internally with
    // Dimer::ComputeMultipoleInteractions();
  }
  else if ( Params::Parameters().GetMMType()==3 ) { //QChem
    CreateQChemJob(Monomers, NMon, true);
  }

  else {
    printf("ERROR:Dimer::CreateMMJob: Unknown MM_type: %d\n",
	   Params::Parameters().GetMMType() );
    exit(1);
  }
}

// Creates Tinker MM job
void Dimer::CreateTinkerJob(Monomer Monomers[], int NMon) {

  // Set up the filenames, with the full path.  
  // Files are e.g. 'd1_2.xyz' and 'd1_2.key' because tinker hates
  // decimal points in filenames
  string filename;
  if (!Params::Parameters().DoFreq() ) {
    filename = Params::Parameters().GetMMPath(); 
  }
  else {
    filename = Params::Parameters().GetHessianMMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    filename += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
 
  filename += "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d_%d",indexA,indexB);
  filename += label;
  string xyzfile = filename + ".xyz";
  string keyfile = filename + ".key";

  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Dimer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }

    if ( Params::Parameters().UseEmbeddingCharges() ) {
    // Count total atoms: monomer atoms + charges
    int Ntot = 0;
    for (int i=1;i<=NMon;i++) {
      Ntot += Monomers[i].GetNumberOfAtoms();
    }
    fprintf(xyz,"%d  Dimer (%d,%d) with embedding charges\n",
	    Ntot,indexA,indexB);

    // Print the geometry
    int shift = 0;
    for (int i=1;i<=NMon;i++) {
      if (i==indexA || i==indexB)
	Monomers[i].PrintTinkerCartesian(shift,false,xyz);
      else
	Monomers[i].PrintTinkerEmbeddingCharges(shift,xyz);

      // Shift for the next monomer
      shift += Monomers[i].GetNumberOfAtoms();
    }
  }
  else
    PrintTinkerCartesian(xyz);


  fclose(xyz);

  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
  if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
    printf("Dimer::CreateTinkerJob : Cannot open file '%s'\n",
	   keyfile.c_str());
    exit(1);
  }
  fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
  fclose(key);
  
}

// Wrapper that controls which type of QM job is executed
// Returns command string for running the job
string Dimer::RunQMJob(bool CCSDT) {
  string job;

  if (Params::Parameters().GetQMType()==1) // QChem
    job = RunQChemJob();
  else if (Params::Parameters().GetQMType()==2){ // MolPro
    if(Params::Parameters().DoFreq())
      job = RunFiniteDifferenceMolProJob();
    else
      job = RunMolProJob(CCSDT);
  } 
  else if(Params::Parameters().GetQMType()==3) { // G09
      job = RunG09Job();
  }
  else if(Params::Parameters().GetQMType()==7) { // PSI4
      job = RunPSI4Job();
  }
  else {
    printf("ERROR: Dimer::RunQMJob(): Unknown QM type.\n");
    exit(1);
  }
  return job;
}



// Returns command string for running QChem job
string Dimer::RunQChemJob(bool MM_job) {
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
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

  string filename = path + "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string infile = filename + ".in";
  string outfile = filename + ".out";

  // Execute Q-Chem
  // Q-Chem wants to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_prefix = "d";
  local_prefix += label;
  string local_outfile = local_prefix + ".out";
  string local_infile = local_prefix +  ".in";
 
  cmd += "qchem " + local_infile;
  cmd += " ";
  cmd += local_outfile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;
}

// Returns command string for running QChem job
string Dimer::RunMolProJob(bool CCSDT) {
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
  string path = Params::Parameters().GetQMPath(); ;
  if(Params::Parameters().DoFreq() )  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string filename =  "d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string infile = path + "/" + filename + ".inp";
  string outfile = path + "/" + filename + ".out";
  string xmlfile = path + "/" + filename + ".xml";
  if(CCSDT){
    infile = path + "/" + filename + ".CCSDT.inp";
    outfile = path + "/" + filename + ".CCSDT.out";
    xmlfile = path + "/" + filename + ".CCSDT.xml";
  }

  // First command, change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  //remove previous .out and .xml
  cmd += "rm -f " + outfile + " " + xmlfile + ";";
  
  //run molpro
  cmd += "molpro " + infile + " --no-xml-output;" ;

  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("%s\n",cmd.c_str());
  return cmd;

}


// Returns command string for running Molpro Finite Difference job
string Dimer::RunFiniteDifferenceMolProJob() {
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
  string path = Params::Parameters().GetQMPath(); ;
  if(Params::Parameters().DoFreq() )
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  path += "/" + dimername;
  string cshfile = dimername + ".csh";

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

// Returns command string for running MolPro job
string Dimer::RunFiniteDifferenceMolProJob(int i,bool Plus,bool CCSDT) {

  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
  string path = Params::Parameters().GetQMPath(); ;
  if(Params::Parameters().DoFreq() )
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  path += "/" + dimername;


  if(CCSDT)
    path += ".CCSDT";

  char num[10];
  sprintf(num,"%d",i);


  string dimerfile;
  if(Plus)
    dimerfile = path + "/" + dimername + "+" + num ;
  else
    dimerfile = path + "/" + dimername + "-" + num ;    

  string infile = dimerfile + ".inp";
  string outfile = dimerfile + ".out";
  string xmlfile = dimerfile + ".xml";


  // First command, change to local directory 
  string cmd = "cd " + path + ";";
  
  //remove previous .out and .xml
  cmd += "rm -f " + outfile + " " + xmlfile + ";";

  //run molpro job
  cmd += "molpro " + infile + ";";
  
  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("%s\n",cmd.c_str());
  //exit(0);
  return cmd;
   
}

// Returns command string for running MolPro job
string Dimer::RunFiniteDifferenceMolProJob(int i,bool Plus1,int j,bool Plus2, bool CCSDT) {

  // Set up the filename, with the full path.  File is e.g. 'd1.2.in'
  string path = Params::Parameters().GetQMPath(); ;
  if(Params::Parameters().DoFreq() )
    path = Params::Parameters().GetHessianQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  path += "/" + dimername;


  if(CCSDT)
    path += ".CCSDT";

  string Sign1 = "-";
  string Sign2 = "-";
  if(Plus1) Sign1 = "+";
  if(Plus2) Sign2 = "+";

  char num1[10];
  sprintf(num1,"%d",i);
  char num2[10];
  sprintf(num2,"%d",j);

  string dimerfile = dimername + Sign1 + num1 + Sign2 + num2;
   

  string infile = path + "/" + dimerfile + ".inp";
  string outfile = path + "/" + dimerfile + ".out";
  string xmlfile = path + "/" + dimerfile + ".xml";
  string logfile = path + "/" + dimerfile + ".log";

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

  //run molpro job
  cmd += "molpro " + infile + ";";
  
  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  //printf("%s\n\n",cmd.c_str());

  return cmd;
}

//Creating PSI4 dimer jobs // CSG
void Dimer::CreatePSI4Job(Monomer Monomers[], int NMon) {

  //Set up the filename, with the full path. File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20];
  sprintf(label, "%d.%d.in",indexA,indexB);
  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreatePSI4Job : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  //input header

  fprintf(job, "# Dimer (%d,%d) \n",indexA,indexB);
  fprintf(job,"memory %i mb\n",Params::Parameters().MemoryUse());
  
  //Print Cartesian coords
  fprintf(job,"molecule dimer {\n");
  fprintf(job,"%i %i\n\n",charge,spin);
  MonA.PrintPSI4MonomerCartesian(job,1);
  fprintf(job,"--\n");
  fprintf(job,"%i %i\n\n",charge, spin);
  MonB.PrintPSI4MonomerCartesian(job,MonA.GetNumberOfAtoms()+1);
  fprintf(job,"}\n\n");

  // Print REM section
  fprintf(job,"%s\n", Params::Parameters().GetPSI4Rem().c_str() );

  fclose(job);
}

// Wrapper that controls which type of MM job is executed
// Returns command string for running the job
string Dimer::RunMMJob() {
  string job;

  if (Params::Parameters().GetMMType()==1) // Tinker
    job = RunTinkerJob();
  else if (Params::Parameters().GetMMType()==2) // AIFF
    job = "echo no dimer job for CamCasp";
  else if (Params::Parameters().GetMMType()==3) // QChem
    job = RunQChemJob(true);
  else {
    printf("ERROR: Dimer::RunMMJob(): Unknown MM type.\n");
    exit(1);
  }
  return job;
}

// Returns command string for running the job
string Dimer::RunTinkerJob() {
  // Set up the filenames, with no path info
  // Files are e.g. 'm1.xyz' and 'm1.key'
  string filename;
  char label[10]; // index ought to have fewer than 10 digits
  sprintf(label,"d%d_%d",indexA,indexB);
  filename = label;
  string infile = filename + ".xyz";
  string outfile = filename + ".out";

  // Execute Tinker
  // Need to be in the local directory, so we have to use local
  // paths, not global ones, and change to the proper directory.

  // First command, change to local directory
  string job_path = Params::Parameters().GetMMPath();
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() ) {
    job_path = Params::Parameters().GetHessianMMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
 
  string cmd = "cd " + job_path;
  cmd += "; ";

  // Second command, run the job
  cmd += "analyze " + infile;
  cmd += " e > ";
  cmd += outfile + "; ";

  // Third command, Rename output from d1_2.out -> d1.2.out
  string newfile = outfile;
  int index = newfile.find("_");
  newfile.replace(index,1,".");
  //printf("new file name = %s\n",newfile.c_str());
  cmd += "mv " + outfile;
  cmd += " " + newfile;
  cmd += "; ";

  /* Actual job running & checking of status have been moved to main.C 
     to simplify parallel implementation.

  printf("Running Tinker energy calculation on %s\n",infile.c_str());
  fflush(stdout); // flush the output stream
  //printf("Executing: %s\n",cmd.c_str());
  system(cmd.c_str());
  */

  // If doing force job, compute the gradient
  if (Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    // Run the job
    cmd += "minimize " + infile;
    cmd += " 100000000 > ";
    
    // Rename the force file
    string force_file = filename + ".force";
    index = force_file.find("_");
    force_file.replace(index,1,"."); // replace "_" with "." in name
    cmd += force_file + ";";

    // Remove tmp file and extra geom file created by minimize job
    cmd += "rm -f " + infile;
    cmd += "_2;";
  }
  if ( Params::Parameters().DoFreq() && !( Params::Parameters().DoFiniteDifferenceFreqs() ) ) {   
    // Run the job
    cmd += "vibrate " + infile;     
    cmd += " 1 > ";        

    // Rename the freq file                
    string freq_file = filename + ".freq";               
    index = freq_file.find("_");
    freq_file.replace(index,1,"."); // replace "_" with "." in name
    cmd += freq_file + ";";            


    // Remove tmp file and extra geom file created by vibrate job
    cmd += "rm -f *.001;"; 
//    cmd += "_2;";
  }


  // Final command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  //  Actual job execution has been moved to main.C

  return cmd;
}

// Read in the Q-Chem energy
double Dimer::ReadQChemEnergy(bool MM_job) {
  double energy, cp_en[4] ={0};
  int cp_en_index = 0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path;
  if (MM_job) 
    path = Params::Parameters().GetMMPath();
  else 
    path = Params::Parameters().GetQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }

  double nuc_chg, chg_chg; // for if using embedding charges

  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);
    string match2 = line.substr(0,11);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      //printf("Std: MATCH FOUND!: %s\n",line.c_str());
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++) 
	iss >> tmp; // throw away text 
     iss >> cp_en[cp_en_index]; // read energy
     //printf("cp_en[%d] = %f\n",cp_en_index,cp_en[cp_en_index]);
    }
    if ( match2=="User input:" ) {
      cp_en_index++;
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

  energy = cp_en[1];
  if ( (energy==0.0 && sym_fac != 0) || (!(energy/energy == 1)&& sym_fac != 0)) {
    printf("Dimer::ReadQChemEnergy : Dimer Single Point Energy Error '%s'\n",
          out_filename.c_str());
    printf("energy = %f\n",energy);
    exit(1);
  }


  if ( Params::Parameters().UseEmbeddingCharges() ) {
    double old = energy;
    energy = energy - (nuc_chg + chg_chg);

    if ( Params::Parameters().PrintLevel() > 0 ) {
      printf("Embedding charges in use.  Correcting Final Energy\n");
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

//Read MolPro Energy
double Dimer::ReadMolProEnergy(){

  double energy = 0.0;
  //Path of QM
  string path = Params::Parameters().GetQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::SetMolProEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;

  int ReadEnergy = 0;

  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,15);
    
    // Search for final MolPro energy
    if ( match==" SETTING EDIMER" ) { 
      ReadEnergy++;
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }
  }
    

  if(ReadEnergy == 0){
    printf("Dimer::SetMolProEnergy : Cannot find energy in '%s'\n",
	   out_filename.c_str());
    exit(1);
  }



  if ( Params::Parameters().PrintLevel() > 0) printf("MolPro Obtained QM Dimer Energy = %15.9f\n",energy);

  // Close the output file
  infile.close();

  //Get CCSD(T) correction
  if(Params::Parameters().DoCCSDTCorrection()){
   out_filename = filename + ".CCSDT.out";

   ReadEnergy = 0;


    // Open the energy file
    ifstream infile;
    infile.open( out_filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Dimer::SetMolProEnergy : Cannot open file '%s'\n",
	     out_filename.c_str());
      exit(1);
    }
    double value = 0;
    while ( !infile.eof() ) {
      getline(infile,line);
      string match = line.substr(0,15);
     // Search for final MolPro energy
      if ( match==" SETTING EDIMER" ) { //JDH change this to what Greg made the string
	ReadEnergy++;
        istringstream iss(line);
        string tmp;
        for (int i=0;i<3;i++)
	  iss >> tmp; // throw away text 
        iss >> value; // read CCSD(T) correction
      }
    }


    if(ReadEnergy == 0){
      printf("Dimer::SetMolProEnergy : Cannot find energy in '%s'\n",
	     out_filename.c_str());
      exit(1);
    }

    //printf("d(%d.%d)",indexA,indexB);
    //printf("pre CCSDT energy = %f\n",energy);
    //printf("CCSDT energy = %f\n",value);
    energy += value; //add the CCSD(T) correction to the energy
    
    //printf("CCSDT energy = %f\n\n",energy);
    infile.close();

  }

  return energy;

}

//Read PSI4 Energy // CSG
double Dimer::ReadPSI4Energy(){

  double energy = 0.0;
  //Path of QM
  string path = Params::Parameters().GetQMPath();
  
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  //File name
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::SetPSI4Energy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;

  int ReadEnergy = 0;

  /*while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,39);
    // Search for final PSI4 energy
    if ( match=="  \"COUNTERPOISE CORRECTED TOTAL ENERGY\"" ) { 
      ReadEnergy++;
      istringstream iss(line);
      string tmp;
      for (int i=0;i<5;i++)
	iss >> tmp; // throw away text 
        iss >> energy; // read energy
        cout << energy << endl; // Just to make sure the right energy is printed
    }
  }*/

  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,18);
    // Search for final PSI4 energy
    if ( match=="  \"CURRENT ENERGY\"" ) { 
      ReadEnergy++;
      istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++)
	iss >> tmp; // throw away text 
        iss >> energy; // read energy
        cout << energy << endl; // Just to make sure the right energy is printed
    }
  }
    

  if(ReadEnergy == 0){
    printf("Dimer::SetPSI4Energy : Cannot find energy in '%s'\n",
	   out_filename.c_str());
    exit(1);
  }



  if ( Params::Parameters().PrintLevel() > 0) printf("PSI4 Obtained QM Dimer Energy = %15.9f\n",energy);

  // Close the output file
  infile.close();
  
  //cout << energy << endl;
  return energy;

}

//Watit
double Dimer::ReadG09Energy() {
  double energy = 0.0;
  string path = Params::Parameters().GetQMPath();
  int dirsize = Params::Parameters().GetBasePath().size();
  string qm_method = path.substr(dirsize);
  transform(qm_method.begin(), qm_method.end(), qm_method.begin(), ::toupper);	
  string out_filename;
  char label[20];
  sprintf(label,"%d.%d.log", indexA, indexB );
  out_filename = path + "/d" + label;

  ifstream infile;
  infile.open( out_filename.c_str() );
  if(!infile.is_open()) {
    printf("Dimer::ReadG09Energy : Cannot open file '%s'\n",
    out_filename.c_str());
    exit(1);
  }

  string line;
  if(Params::Parameters().DoCounterpoise()) {
    while(!infile.eof()) {
      getline(infile,line);
      if(line.find("Counterpoise corrected") != string::npos) {
        istringstream iss(line);
        string tmp;
        for(int i=0;i<4;i++) {
          iss >> tmp; // throw away text 
        }
        iss >> energy; // read energy
      }
    }	
  } else {
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
  }

  //printf("d%i.%i\n",indexA,indexB);
  //cout << "energy " << energy << endl;

  infile.close();
  return energy;
}

double Dimer::ReadQChemCounterpoiseInteractionEnergy() {

  double energy=0, cp_e[3]={0};
  int cp_energy_index=0; //goes from 0 to 2

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  string filename = path + "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  //cout << " out_filename = " << out_filename << "\n";
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : Cannot open file '%s'\n",
           out_filename.c_str());
    exit(1);
  }

  //Remember that now our .out files has 3 energies and the corrected
  //dimer SPE is now E0-E1-E2 which we shall use to calculate
  //interaction energies.  Also, interaction energy will come out to
  //be equal to E0-E1-E3 (so no need to calculate monomer energies at
  //monomer basis.  cp_e0 = dimer energy, cp_e1 = mon1 energy at dimer
  //basis, cp_e2 = mon2 energy at dimer basis
  
  // Read in the data
  string line;           
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,22);
    string match2 = line.substr(0,3);

    // Search for final q-chem energy
    if ( match==" Q-Chem Final Energy =" ) {
      //printf("CP: MATCH FOUND!!!  line=%s\n",line.c_str());
      istringstream iss(line);
      string tmp;
      for (int i=0;i<4;i++)
        iss >> tmp; // throw away text
      iss >> cp_e[cp_energy_index]; // read energy 
      //printf("cp_en[%d] = %f\n",cp_energy_index,cp_e[cp_energy_index]); 
      cp_energy_index++;
   }
  
  }
  
  if ( !(cp_energy_index == 3) || !(cp_e[0]/cp_e[0] == 1) || !(cp_e[1]/cp_e[1] == 1) || !(cp_e[2]/cp_e[2] == 1)  ) {
     printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : At Least One Single Point Energy Not Found or Incorrect '%s'\n",
           out_filename.c_str());
     exit(1);            
  }

  energy = cp_e[0] -cp_e[1]-cp_e[2];
  //cout << "energy = " << energy << "\n"; 
  if ( Params::Parameters().PrintLevel() > 0)
    printf("d(%d,%d) QM CP-corrected interaction Energy = %15.9f\n",
	   indexA,indexB,energy);

  // Convert kJ/mol -> hartrees
//  energy /= HartreesToKJpermole;
//  cout << "energy = " << energy << "\n";

  // Close the force file
  infile.close();

  return energy;

}



/* Following method of calculating CP corrected dimer energy is commented out
// Read in the Q-Chem counterpoise-corrected interaction energy
// Not sure how this works with embedding charges.
double Dimer::ReadQChemCounterpoiseInteractionEnergy() {

  double energy;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  string path = Params::Parameters().GetQMPath();
  string filename = path + "/d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  string out_filename = filename + ".out";
  
  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemCounterpoiseInteractionEnergy : Cannot open file '%s'\n",
	   out_filename.c_str());
    exit(1);
  }


  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    string match = line.substr(0,12);

    // Search for final q-chem energy
    if ( match=="  DE, kJ/mol" ) {
            istringstream iss(line);
      string tmp;
      for (int i=0;i<3;i++) 
	iss >> tmp; // throw away text 
      iss >> energy; // read energy
    }

  }

 
  if ( Params::Parameters().PrintLevel() > 0) 
    printf("QM CP-corrected interaction Energy = %15.9f\n",energy);
 
  // Convert kJ/mol -> hartrees
  energy /= HartreesToKJpermole;

  // Close the force file
  infile.close();

  return energy;
}
*/

void Dimer::ReadQMResults() {

 if ( Params::Parameters().PrintLevel() > 0) {
      //shows which dimers and energy of dimer 
    printf("Dimer(%d,%d) energies \n", MonA.GetIndex(), MonB.GetIndex());
  }
  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().TinkerDebug() ) {
    Energy_QM = ReadTinkerEnergy();
  }
  else if (sym_fac != 0){//only find energy if symmetry factor is not zero
    
    if(Params::Parameters().GetQMType() == 1){//Qchem
    	Energy_QM = ReadQChemEnergy();
     if (Params::Parameters().DoCounterpoise()) {
        dEint_QM = ReadQChemCounterpoiseInteractionEnergy();
     //printf("Counterpoise-corrected interaction energy = %15.9f\n",dEint_QM);
     }
     else
        dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
    }
    else if(Params::Parameters().GetQMType() == 2){//MolPro
      	Energy_QM = ReadMolProEnergy(); 
        if(Params::Parameters().DoCounterpoise())
          dEint_QM = Energy_QM;
        else
	  dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
    }
    else if(Params::Parameters().GetQMType()==3) { // G09
      Energy_QM = ReadG09Energy();
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
    }
    else if(Params::Parameters().GetQMType()==7) { // PSI4
      Energy_QM = ReadPSI4Energy();
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
    }
    else {
      printf("ERROR::Dimer::SetQMResults() : Unknown QM_type: %d\n",Params::Parameters().GetQMType());
      exit(1);
    }
  }else{
    dEint_QM = 0;
    Energy_QM = 0;
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("dEint_QM = %15.9f  ",dEint_QM);
    printf("Energy_QM = %f, monA = %f, monB = %f\n",Energy_QM, MonA.GetQMEnergy(),MonB.GetQMEnergy());
  }

  //Set Gradient
  if ( Params::Parameters().DoForces() ||Params::Parameters().DoFreq() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    if ( Params::Parameters().TinkerDebug() ) 
      SetMMGradient();
    else 
      SetQMGradient();
  }
  //Set Hessian
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() &&
       !(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) {
    if ( Params::Parameters().TinkerDebug() )
      SetMMHessian();
    else{
      SetQMHessian();
    }
  }

}

// Use "other" monomer as a reference for Mon B;
void Dimer::ReadQMResults(Monomer& otherB) {

  //printf("Using Monomer %d as a reference\n",otherB.GetIndex() ); fflush(stdout);

  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().TinkerDebug() ) {
    Energy_QM = ReadTinkerEnergy();
  }
  else if( Params::Parameters().GetQMType() == 1){//QChem
    Energy_QM = ReadQChemEnergy();

    if (Params::Parameters().DoCounterpoise()) {
      dEint_QM = ReadQChemCounterpoiseInteractionEnergy();
      //printf("Counterpoise-corrected interaction energy = %15.9f\n",dEint_QM);
    }
    else {
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - otherB.GetQMEnergy();
    }
    
  }
  else if( Params::Parameters().GetQMType() == 2){ //MolPro
    Energy_QM = ReadMolProEnergy(); 
    if(Params::Parameters().DoCounterpoise())
      dEint_QM = Energy_QM;
    else
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
  }
  else if(Params::Parameters().GetQMType()==3) { // G09
      Energy_QM = ReadG09Energy();
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
  }
  else if(Params::Parameters().GetQMType()==7) {//PSI4
      Energy_QM = ReadPSI4Energy();
      dEint_QM = Energy_QM - MonA.GetQMEnergy() - MonB.GetQMEnergy();
  }
  else {
    printf("ERROR::Dimer::SetQMResults() : Unknown QM_type: %d\n",Params::Parameters().GetQMType());
    exit(1);
  }


  if ( Params::Parameters().PrintLevel() > 0) {
    printf("dEint_QM = %15.9f  ",dEint_QM);
    printf("Energy_QM = %f, monA = %f, monB = %f\n",Energy_QM,
	   MonA.GetQMEnergy(),otherB.GetQMEnergy());
    }
  
  //Set Gradient
  if ( Params::Parameters().DoForces() ||Params::Parameters().DoFreq() || Params::Parameters().DoFiniteDifferenceFreqs()) {
    if ( Params::Parameters().TinkerDebug() ) {
      SetMMGradient();
    }
    else {
      SetQMGradient();
    }
    MonB.SetQMGradient(otherB);
  }

  //Set Hessian
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() &&
       !(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable()) ) {
    if ( Params::Parameters().TinkerDebug() ) {
      SetMMHessian();
    }
    else {
      SetQMHessian();
    }
    MonB.SetQMHessian(otherB);
  }   
}



void Dimer::ReadMMResults() {

  // Read in the Energies, compute 2-body interaction energies
  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient(); // get energies & gradient simultaneously in this case
  }
  else {
    SetMMEnergy();
    if (Params::Parameters().DoForces()  || Params::Parameters().DoFreq() ||Params::Parameters().DoFiniteDifferenceFreqs()) {
      SetMMGradient();
    }
  }

  //Read in the Energies and Hessian
  if ( Params::Parameters().GetMMType()==2 && Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() &&
    !(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) {
    SetMMHessian(); // get energies & gradient simultaneously in this case
  }
  else {
    SetMMEnergy();
    if (Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs() && 
      !(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())) {
      SetMMHessian();
    }
  }

  dEint_MM = Energy_MM - MonA.GetMMEnergy() - MonB.GetMMEnergy();
   
}

// Version for use with PBC: Monomer B is a periodic image of
// reference monomer otherB
void Dimer::ReadMMResults(Monomer& otherB) {

  // Read in the Energies, compute 2-body interaction energies
  SetMMEnergy();

  dEint_MM = Energy_MM - MonA.GetMMEnergy() - otherB.GetMMEnergy();


  //Find Gradient if requested
  if ( Params::Parameters().DoForces() || Params::Parameters().DoFreq() ||Params::Parameters().DoFiniteDifferenceFreqs() ) {
    SetMMGradient();
    MonB.SetMMGradient(otherB);
  }

  //Find Hessian if requested
  if ( Params::Parameters().DoFreq() && !Params::Parameters().DoFiniteDifferenceFreqs()) {
    SetMMHessian();
    MonB.SetMMHessian(otherB);
  }

}

void Dimer::SetMMEnergy() {
  double energy;
  if ( Params::Parameters().GetMMType() == 1 ) // Tinker
    energy = ReadTinkerEnergy();
  else if ( Params::Parameters().GetMMType() == 2) { // AIFF
    //energy = ComputeMultipoleInteractions();
    energy = ComputeAIFFEnergy();
  }
  else if ( Params::Parameters().GetMMType() == 3) { // QChem
    energy = ReadQChemEnergy(true);
  }
  else {
    printf("Dimer::SetMMEnergy: Unknown MM_type: %d\n",Params::Parameters().GetMMType());
    exit(1);
  }

  Energy_MM = energy;
  //printf("MM Energy = %15.9f\n",Energy_MM);


}

double Dimer::ReadTinkerEnergy() {

  string path = Params::Parameters().GetMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  // Search for the "Polarization" string
  double energy = 0.0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"/d%d.%d",MonA.GetIndex(), MonB.GetIndex() );
  string filename = path + label;
  string out_filename = filename + ".out";

  // Open the energy file
  ifstream infile;
  infile.open( out_filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadTinkerEnergy : Cannot open file '%s'\n",out_filename.c_str());
    exit(1);
  }

  // Read in the data
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);

    /*
    string match = line.substr(0,13);

    if ( match==" Polarization" ) {
      istringstream iss(line);
      string tmp; 
      iss >> tmp; // throw away "Polarization" tag
      iss >> energy; // read the energy
    }
    */
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
  return energy;
}

// Get the QM Gradient, wrapper routine
void Dimer::SetQMGradient() {
  string path = Params::Parameters().GetQMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  Grad_QM = ReadGradient(path,1);
  QM_Grad_Init = 1;
}

// Get the MM Gradient, wrapper routine
void Dimer::SetMMGradient() {

  if (Params::Parameters().GetMMType() == 1) { // Tinker 
    string path = Params::Parameters().GetMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Grad_MM = ReadGradient(path,2);
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    Grad_MM = ComputeMultipoleGradient();
  }
  else if (Params::Parameters().GetMMType() == 3) { //QChem
    string path = Params::Parameters().GetMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Grad_MM = ReadGradient(path,1);
  }
  else {
    printf("Dimer::SetMMGradient() - MM type = %d is unknown\n",
	   Params::Parameters().GetMMType());
    exit(1);
  }

  MM_Grad_Init = 1;
}

// Main routine for reading in gradient files
// Assumes Nx4 structure, where each row has an index, then GX, GY, GZ.
Vector Dimer::ReadGradient(string path, int type) {
  
  Vector grad(3*Natoms), gradA( 3*Natoms ),  gradB( 3*Natoms );//TempGrad(3*Natoms),TempGradA(3*Natoms),TempGradB(3*Natoms);
  int cp_grad_index =0;
  
  // Set up the filename, with the full path.  File is e.g. 'd1.2.force'
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d",indexA,indexB);
  fflush(stdout);
  string filename = path + label;
  if (type == 2) // Tinker MM job
    filename += ".force";
  else if(Params::Parameters().GetQMType()==3) // G09
    filename += ".log";
  else {// QChem or Molpro or PSI4
    filename += ".out";
  }
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadGradient : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  if((!(Params::Parameters().DoCounterpoise())&&!(Params::Parameters().GetQMType()==3))&&!(Params::Parameters().GetQMType()==7)||type==2) {
  //if ( !(Params::Parameters().DoCounterpoise()) || type==2 ) {
    // Read in the data - search for the "Nuclear forces:" string
    
    //for CBS, number of Grad read
    int NumberGrad = 0;
    
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      
      
      // Search for final q-chem or tinker energy
      if(Params::Parameters().GetQMType() == 1 ||type==2 ){ //Qchem and Tinker.
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
		grad[3*i+j] +=entry;//gradients will be rotated if need due to symmetry 
		
		// Store the gradient elements
	      }
	      
	    }	
	    if(NumberGrad == 4){
	      break;
	    }
	  }  
	} //End CBS
        else{
	  if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	    getline(infile,line);// throw away header line
	    
	    for (int i=0;i<Natoms;i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=0;j<3;j++) {
		iss >> grad[3*i+j];//gradients will be rotated if need due to symmetry 
		// Store the gradient elements
	      }
	    }
	    break;
	  }
	}
      }
      //printf("d(%i,%i)\n",GetIndexA(),GetIndexB());
      //RotationA.Print("RotationA");
      //RotationB.Print("RotationB");
      //PrintGradient(" DimerGradient",grad);
      fflush(stdout);
    }
  } //End Molpro
  
  if ( Params::Parameters().DoCounterpoise() && type==1 && !(Params::Parameters().GetQMType()==3) && !(Params::Parameters().GetQMType()==7) ) {
    if(Params::Parameters().GetQMType() == 1){
      // Read in the data - search for the "Nuclear forces:" string
      string line;               
      while ( !infile.eof() ) {
        getline(infile,line);
        // Search for final q-chem energy
        if ( line==" Nuclear forces:" && (cp_grad_index == 0)) {
	  getline(infile,line); // throw away header line
	  for (int i=0;i<Natoms;i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=0;j<3;j++) {
              iss >> grad[3*i+j];  // Store the gradient elements
	    }
	  }
	  cp_grad_index++;
        }
	
        if ( line==" Nuclear forces:" && (cp_grad_index == 1) ) {
	  getline(infile,line); // throw away header line
	  for (int i=0;i<(MonA.GetNumberOfAtoms()+MonB.GetNumberOfAtoms());i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=0;j<3;j++) {
              iss >> gradA[3*i+j]; // Store the gradient elements
	    }
	  }
	  cp_grad_index++;
        }
	
        if ( line==" Nuclear forces:" && (cp_grad_index == 2) ) {
	  getline(infile,line); // throw away header line
	  for (int i=0;i<(MonA.GetNumberOfAtoms()+MonB.GetNumberOfAtoms());i++) {
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=0;j<3;j++){
              iss >> gradB[3*i+j]; // Store the gradient elements
	    }
	  }
	  cp_grad_index++;
        }
      }
      
    } 
    else if(Params::Parameters().GetQMType() == 2){
      if(Params::Parameters().DoCBS()){
	double basis1 = double(Params::Parameters().CBSBasis1());
	double basis2 = double(Params::Parameters().CBSBasis2());
	double ExpTerm = exp(1.54*(basis2-basis1));
	string line;
	while ( !infile.eof() ) {
	  getline(infile,line);	  
	  if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	    cp_grad_index++;	    
	    getline(infile,line);//discard first line
	    for (int i=0;i<Natoms;i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=0;j<3;j++) {
		double entry;
		iss >> entry;
		//HF in the first basis for dimer 
		if(cp_grad_index == 1){
		  entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  grad[3*i+j] +=entry;
		}
		//MP2 in the first basis for dimer
		else if(cp_grad_index == 2){
		  entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  grad[3*i+j] += entry;
		}
		//HF in the first basis with Ghost MonA
		else if(cp_grad_index == 3){
		  entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  gradA[3*i+j] +=entry;
		}
		//MP2 in the first basis with Ghost MonA
		else if(cp_grad_index == 4){
		  entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  gradA[3*i+j] += entry;
		}
		//HF in the first basis with Ghost MonB
		else if(cp_grad_index == 5){
		  entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  gradB[3*i+j] +=entry;
		}
		//MP2 in the first basis with Ghost MonB
		else if(cp_grad_index == 6){
		  entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  gradB[3*i+j] += entry;
		}
		//HF in the second basis for dimer
		else if(cp_grad_index == 7){
		  entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  grad[3*i+j] +=entry;
		}
		//MP2 in the second basis for dimer 
		else if(cp_grad_index == 8){
		  entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  grad[3*i+j] +=entry;
		}
		//HF in the second basis with Ghost MonA
		else if(cp_grad_index == 9){
		  entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  gradA[3*i+j] +=entry;
		}
		//MP2 in the second basis for Ghost MonA
		else if(cp_grad_index == 10){
		  entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  gradA[3*i+j] +=entry;
		}
		//HF in the second basis with Ghost MonB
		else if(cp_grad_index ==11){
		  entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  gradB[3*i+j] +=entry;
		}
		//MP2 in the second basis for Ghost MonB
		else if(cp_grad_index == 12){
		  entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  gradB[3*i+j] +=entry;
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
	  if(line==" Atom          dE/dx               dE/dy               dE/dz" && (cp_grad_index == 0)){
	    getline(infile,line);// throw away header line 
	    for (int i=0;i<Natoms;i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=0;j<3;j++) {
		iss >> grad[3*i+j]; 
		// Store the gradient elements
	      }
	    }  
	    cp_grad_index++;
	  }
	  
	  if(line==" Atom          dE/dx               dE/dy               dE/dz" && (cp_grad_index == 1)){
	    getline(infile,line);// throw away header line 
	    for (int i=0;i<Natoms;i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=0;j<3;j++) {
		iss >> gradA[3*i+j];
		// Store the gradient elements
	      }
	    } 
	    cp_grad_index++;
	  }
	  
	  if(line==" Atom          dE/dx               dE/dy               dE/dz" && (cp_grad_index == 2)){
	    getline(infile,line);// throw away header line 
	    for (int i=0;i<Natoms;i++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		iss >> tmp; // throw away the atom index
		for (int j=0;j<3;j++) {
		  iss >> gradB[3*i+j];
		  // Store the gradient elements
		}
	    }  
	    cp_grad_index++;
	  }
	  
	    
	}
      }
    }
    if((!(cp_grad_index==12)&&Params::Parameters().DoCBS())||(!(cp_grad_index==3)&&!Params::Parameters().DoCBS())&&!(Params::Parameters().GetQMType()==3)) {
      printf("Dimer::ReadGradient : At Least One Counterpoise Gradient Not Found or Incorrect '%s'\n",filename.c_str());
      printf("cp_grad_index = %i\n",cp_grad_index);
      exit(1);
    }
    
    //PrintGradient("Dimer::ReadGradient(): full  Dimer Gradient",grad);    
    //PrintGradient("Dimer::ReadGradient(): Dimer A Gradient",gradA);  
    //PrintGradient("Dimer::ReadGradient(): Dimer B Gradient",gradB); 

    // Find the CP gradient.
    for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
      grad[j] = grad[j] - gradA[j] - gradB[j];
    }
    for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
      int k = j + 3*MonA.GetNumberOfAtoms();
      grad[k] = grad[k] - gradA[k] - gradB[k];
    }
  }// end Counterpoise Gradient section

  if(type==1&&Params::Parameters().GetQMType()==3) { //G09 Watit
    string line;
    while(!infile.eof()) {
      getline(infile,line);
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
            iss >> grad[3*i+j];
	    grad[3*i+j]=-grad[3*i+j];
          }
        }
      }
    }

    if(Params::Parameters().DoCounterpoise()) {
      for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
        grad[j] = grad[j] - MonA.GetQMGradient()[j];
      }
      for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
        int k = j + 3*MonA.GetNumberOfAtoms();
        grad[k] = grad[k] - MonB.GetQMGradient()[j];
      }
    }
  }

  if(type==1&&Params::Parameters().GetQMType()==7) { //PSI4 CSG
    string line;
    while(!infile.eof()) {
      getline(infile,line);
      cout << line << endl;
      if(line.find("-Total Gradient:") !=string::npos) {
        getline(infile,line);
        getline(infile,line);
        //getline(infile,line);
        for (int i=0;i<Natoms;i++) {
          getline(infile,line);
          istringstream iss(line);
          string tmp;
	  iss >> tmp;
          for (int j=0;j<3;j++) {
            // Psi4 prints gradients in kcal/mol*A, but HMBI want hartree/bohr
            // how it was originally
            iss >> grad[3*i+j];
            // PSI4 gradients are printed in kcal/mol*A
            // Need to convert gradients to hartree/bohr
            double temp = grad[3*i+j]*(1/(627.509*1.8897259886));
            grad[3*i+j] = temp;
	    //grad[3*i+j]=-grad[3*i+j];
	//    cout << grad[3*i+j] << endl;
          }
        }
      }
    }

    if(Params::Parameters().DoCounterpoise()) {
      for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
        grad[j] = grad[j] - MonA.GetQMGradient()[j];
      }
      for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
        int k = j + 3*MonA.GetNumberOfAtoms();
        grad[k] = grad[k] - MonB.GetQMGradient()[j];
      }
    }
  }

  //printf("Gradient for dimer(%i,%i)\n", indexA, indexB);
  //PrintGradient("Test 2 Dimer Gradient",grad);


  infile.close();  // Close the file


  //Include the CCSD(T) Correction
  //CCSD(T) has not been adapted to q-chem 
  if(type == 1 && Params::Parameters().DoCCSDTCorrection()){


    if(Params::Parameters().DoCounterpoise()){
      //PrintGradient("Dimer::ReadGradient(): before CCSD(T)",grad);
      Vector CCSDTGradient = ReadFiniteDifferenceMolProGradient(path,true);
      for(int i=0;i<3*Natoms;i++)
	grad[i] += CCSDTGradient[i];
      
    }
    else{
      filename = path + label;
      filename += ".CCSDT.out";
      // Open the force file  
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Dimer::ReadGradient : Cannot open file '%s'\n",filename.c_str());
	exit(1);
      }
      int NumberGrad = 0;
      //counterpoise corrected CCSD(T) gradient handled by  ReadFiniteDifferenceMolProGradient()
      while ( !infile.eof() ) {
	string line;
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
	      grad[3*i+j] +=entry;//gradients will be rotated if needed due to symmetry 
	      
            }
          }
        }
      }

      //Wrong number of Gradients
      if(!(NumberGrad==2)){
	 printf("Dimer::ReadGradient : At Least One Gradient Not Found or Incorrect '%s'\n",
		filename.c_str());
	 printf("NumberGrad = %i\n",NumberGrad);
	 exit(1);
      }
      infile.close();  // Close the file
      //PrintGradient("Dimer::ReadGradient(): Dimer Gradient",grad);
    }
  }

  
  //PrintGradient("Dimer::ReadGradient(): Dimer Gradient",grad);	

  return grad; 
}

//Currently used to find CCSD(T) correction gradient using counterpoise correction
Vector Dimer::ReadFiniteDifferenceMolProGradient(string path,bool CCSDT){

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  Vector Grad(3*(Na+Nb));

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = path + "/" + dimername;

  if(CCSDT)
    dimerpath += ".CCSDT";



  for(int i=0;i<3*(Na+Nb);i++){
    double PlusEnergy;
    double MinusEnergy;

    // The plus input file; 
    char num[10];
    sprintf(num,"%d",i);
    string filename = dimerpath + "/" + dimername + "+" + num + ".out";

    // Open the .out file
    ifstream infile;                 
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Dimer::ReadFiniteDifferenceMolProGradient : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      string match = line.substr(0,15);
      
      // Search for final MolPro energy
      if ( match==" SETTING EDIMER" ) {
	istringstream iss(line);
	string tmp;
	for (int i=0;i<3;i++)
	  iss >> tmp; // throw away text 
	iss >> PlusEnergy; // read energy
      }
    }
    infile.close();

    // The minus input file; 
    filename = dimerpath + "/" + dimername + "-" + num + ".out";
    // Open the .out file 
    ifstream CCSDTfile;                
    CCSDTfile.open( filename.c_str() );
    if ( !CCSDTfile.is_open() ) {
      printf("Dimer::ReadFiniteDifferenceMolProHessian : Cannot open file '%s'\n",
	     filename.c_str());
      exit(1);
    }

    while ( !CCSDTfile.eof() ) {
      getline(CCSDTfile,line);
      string match = line.substr(0,15);
      
      // Search for final MolPro energy
      if ( match==" SETTING EDIMER" ) { //JDH change this to what Greg made the string
	istringstream iss(line);
	string tmp;
	for (int i=0;i<3;i++)
	  iss >> tmp; // throw away text 
	iss >> MinusEnergy; // read energy
      }
    }
    CCSDTfile.close();
    
    Grad[i] = MinusEnergy - PlusEnergy;

  } 
  Grad.Scale(1.0/(2*0.001*AngToBohr));
  //Grad.PrintGradient("CCSD(T)");

  return Grad;

}

/*Vector Dimer::ReadQMPGradient(string path) {

  Vector grad(3*Natoms), gradA( MonA.GetNumberOfAtoms() ),  gradB( MonB.GetNumberOfAtoms() );
  int cp_grad_index =0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.force'
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d",indexA,indexB);         
  string filename = path + label;       
  filename += ".out";  


  // Open the force file   
  ifstream infile;         
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadGradient : Cannot open file '%s'\n",  
           filename.c_str());
    exit(1);
  }

  // Read in the data - search for the "Nuclear forces:" string
  string line;               
  while ( !infile.eof() ) {
    getline(infile,line);
    // Search for final q-chem energy
    if ( line==" Nuclear forces:" && (cp_grad_index == 0) ) {
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
//      break;
      cp_grad_index++;
    }

    if ( line==" Nuclear forces:" && (cp_grad_index == 1) ) {
      getline(infile,line); // throw away header line

      for (int i=0;i<MonA.GetNumberOfAtoms();i++) {                
        getline(infile,line);           
        istringstream iss(line);
        string tmp;
        iss >> tmp; // throw away the atom index
        for (int j=0;j<3;j++) {
          iss >> gradA[3*i+j]; // Store the gradient elements
        }
      }
//      break;
      cp_grad_index++;
    }

    if ( line==" Nuclear forces:" && (cp_grad_index == 2) ) {
      getline(infile,line); // throw away header line

      for (int i=0;i<MonB.GetNumberOfAtoms();i++) {                
        getline(infile,line);           
        istringstream iss(line);
        string tmp;
        iss >> tmp; // throw away the atom index
        for (int j=0;j<3;j++) {
          iss >> gradB[3*i+j]; // Store the gradient elements
        }
      }
//      break;
      cp_grad_index++;
    }


  }

  if ( !(cp_grad_index == 3) ) {
   printf("Dimer::ReadQMCPGradient : At Least One Gradient Not Found or Incorrect '%s'\n",
           filename.c_str());
     exit(1);
  }

  for (int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
    grad[j] = grad[j] - gradA[j];
  }

  for (int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
    int k = j + 3*MonA.GetNumberOfAtoms(); // in dimer gradient, need offset to second monomer
    grad[k] = grad[k] - gradB[j];
  }

  infile.close();  // Close the file

  //PrintGradient("Dimer Gradient",grad);
  return grad;
}

*/

// Print out the Cartesian coordinates in xyz format
void Dimer::PrintDimerCartesian() {
  MonA.PrintMonomerCartesian();
  MonB.PrintMonomerCartesian();
}

// Print out the Cartesian coordinates in xyz format
void Dimer::PrintDimerCartesian(FILE *outfile) {
  MonA.PrintMonomerCartesian(outfile);
  MonB.PrintMonomerCartesian(outfile);
}

// Print out $molecule section for Q-Chem
void Dimer::PrintQChemCartesian(FILE *outfile) {
  fprintf(outfile,"$molecule\n%d %d\n",charge,spin);

/*  if ( Params::Parameters().DoCounterpoise() ) {
    fprintf(outfile,"--\n%d %d\n",MonA.GetChargeState(), MonA.GetSpinState());
  }
*/
   MonA.PrintMonomerCartesian(outfile);
/*
  if ( Params::Parameters().DoCounterpoise() ) {
    fprintf(outfile,"--\n%d %d\n",MonB.GetChargeState(), MonB.GetSpinState());
  }
*/
   MonB.PrintMonomerCartesian(outfile);
   fprintf(outfile,"$end\n\n");

 //  if ( Params::Parameters().DoCounterpoise() ) { // make special qchem input file to grab bsse correction manually (added on 9th april)
      


}

// Print out $molecule section for Q-Chem
void Dimer::PrintQChemGhostCartesian(FILE *outfile,bool printGhostMonA,bool printGhostMonB) {
  fprintf(outfile,"$molecule\n%d %d\n",charge,spin);
   
   if ( printGhostMonA )
     MonA.PrintGhostMonomerCartesian(outfile);
   else
     MonA.PrintMonomerCartesian(outfile);

   if ( printGhostMonB ) 
     MonB.PrintGhostMonomerCartesian(outfile);
   else
     MonB.PrintMonomerCartesian(outfile);

   fprintf(outfile,"$end\n\n");


}

// Print out the Cartesian coordinates in Tinker xyz format
void Dimer::PrintTinkerCartesian(FILE *outfile) {
  // Print tile line
  fprintf(outfile,"%d  Dimer (%d,%d)\n",Natoms,indexA,indexB); 

  // Print first and second monomers
  MonA.PrintTinkerCartesian(0,false,outfile);
  // apply shift to second monomer to get proper indexing
  int shift = MonA.GetNumberOfAtoms();
  MonB.PrintTinkerCartesian(shift,false,outfile);

}



// Wrapper for printing QM gradient
void Dimer::PrintQMGradient(string title) {
  printf("QM_Grad_Init = %d\n",QM_Grad_Init);  
  if (QM_Grad_Init) 
    PrintGradient(title,Grad_QM);
  else 
    printf("QM Gradient not initialized\n");
}

// Wrapper for printing MM gradient
void Dimer::PrintMMGradient(string title) {
  printf("MM_Grad_Init = %d\n",MM_Grad_Init);
  if (MM_Grad_Init) 
    PrintGradient(title,Grad_MM);
  else 
    printf("MM Gradient not initialized\n");
}

// Print out Gradient as an Nx3 matrix
void Dimer::PrintGradient(string title, Vector grad) {
  //printf("Dimer::PrintGradient incomplete\n");

  printf("%s\n",title.c_str());
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  for (int i=0;i<Na;i++)
    printf("%2s %15.10f %15.10f %15.10f\n",MonA.GetSymbol(i).c_str(),grad[3*i],grad[3*i+1],grad[3*i+2]);
  for (int i=Na;i<Na+Nb;i++)
    printf("%2s %15.10f %15.10f %15.10f\n",MonB.GetSymbol(i-Na).c_str(),grad[3*i],grad[3*i+1],grad[3*i+2]);

  
}

void Dimer::PrintAll() {

  printf("-----------------------------------------------\n");
  printf("Dimer (%d,%d)      %d atoms\n\n",indexA,indexB,Natoms);

  PrintQChemCartesian();
  printf("Intermolecular separation = %.3f Angstroms\n\n",Separation);
  
  printf("Energies:\nQM: %15.9f\nMM: %15.9f\n",Energy_QM,Energy_MM);
  printf("Interaction Energies:\nQM: %15.9f\nMM: %15.9f\n\n",dEint_QM,dEint_MM);

  if ( Params::Parameters().DoForces() || Params::Parameters().DoFiniteDifferenceFreqs() ) {
    printf("Forces:\n");
    PrintQMGradient("QM:\n");
    PrintMMGradient("MM:\n");
  }



  printf("-----------------------------------------------\n");
}

// Get the damping factor for smooth pairwise local truncations
// Note: This is NOT the same damping factor used to damp short range
// induction energies.  That one is stored in the Param class and
// computed in the Atom class.
double Dimer::GetDampingFactor(double c0, double c1) {
  /* 
     Applies a smooth bump function between distances c1 and c0
     to smoothly damp from 1 -> 0, according to:

     Subotnik et al, J. Chem. Phys. 128, 034103 (2008)

     c0 = distance after which damp = 0  (complete damping)
     c1 = distance before which damp = 1 (no damping)
     Between c1 and c0, we have smooth decay from 1 -> 0.
   */

  //check that cutoff0 is not less than cutoff1
  if(c0 < c1){
    printf("ERROR::Dimer::GetDampingFactor(): cutoff0 is less than cutoff1. Alter cutoffs and rerun.\n");
    exit(0);
  }

  double damp;
  double R =  Separation; //distance between the two fragments

  if (R <= c1) {
    damp = 1.0; 
  }
  else if (R >= c0) {
    damp = 0.0; 
  }
  else {
    double tmp = 2*fabs(c1 - c0)/(c1 - R) - fabs(c1 - c0)/(R - c0);
    damp = 1.0/(1.0 + exp(tmp));
  }

  if ( Params::Parameters().PrintLevel() > 0 )
    printf("GetDampingFactor: d(%d,%d)  R = %.3f, damp = %.4f\n",
	   indexA,indexB,R,damp);
	
  return damp;

}


Vector Dimer::GetSpatialDampingFunctionGradient(double c0, double c1) {


  /*
    This subroutine returns the gradient of the spatial damping function
    that should be included in the HMBI gradient formalism
  */

  //check that cutoff0 is not less than cutoff1
  if(c0 < c1){
    printf("ERROR:Dimer:SpatialDampingFunctionGradient(): cutoff0 is less than cutoff1. Alter cutoffs and rerun.\n");
    exit(0);
  }

  //initialize this gradient
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  //Vector grad_temp(3*(Na+Nb));
  Vector grad_damping(3*(Na+Nb));
  double R =  Separation; //distance between the two fragments
  double x1,y1,z1,x2,y2,z2; // xyz coords of the two atoms defining the minimum distance between the two monomers

  // Now first determine the two atoms which define the shortest distance between the two monomers
  // So, we build subroutines called MonA.MinDistAtomIndex_A(MonB) and MonA.MinDistAtomIndex_B(MonB)
  // in monomer.C and then we already used these to define mindist_atomA and mindist_atomB,
  // which are the indices of these two atoms in the respective monomers.
  // We now use these indices to build the grad_damping
  // Note that only 6 terms in this gradient will be non-zero.

  //printf("Damp Grad for D(%i,%i)\n",indexA,indexB);
  //printf("Damp = %f\n",GetDampingFactor(c0,c1));
  // printf("QM Interaction = %f\n",GetQMIntEnergy());
  //printf("MM Interaction = %f\n",GetMMIntEnergy());
  // printf("minA = %i, minB = %i\n",minA,minB);
  //fflush(stdout);

  x1 =  MonA.GetAtom(minA).GetPosition(0);
  y1 =  MonA.GetAtom(minA).GetPosition(1);
  z1 =  MonA.GetAtom(minA).GetPosition(2);
  x2 =  MonB.GetAtom(minB).GetPosition(0);
  y2 =  MonB.GetAtom(minB).GetPosition(1);
  z2 =  MonB.GetAtom(minB).GetPosition(2);
/*
  cout << " mindist_atomA = " << mindist_atomA << "\n";
  cout << " mindist_atomB = " << mindist_atomB << "\n";

  cout << "x1 = " <<x1 << "\n";
  cout << "y1 = " <<y1 << "\n"; 
  cout << "z1 = " <<z1 << "\n";
  cout << "x2 = " <<x2 << "\n";
  cout << "y2 = " <<y2 << "\n";               
  cout << "z2 = " <<z2 << "\n";
*/


  // Getting rotation matrix
  //Matrix RotA = MonA.GetRotationMatrix();
  //Matrix RotB = MonB.GetRotationMatrix();

  // RotA.Print("RotA");
  //RotB.Print("RotB");

  /*
     The following terms are non-zero
     3*mindist_atomA, 3*mindist_atomA+1, 3*mindist_atomA+2,
     3*Na+3*mindist_atomB, 3*Na+3*mindist_atomB+1, 3*Na+3*mindist_atomB+2
  */
  // Now write the expressions for these 6 terms using brute force method

  if ( !(R<=c1) && !(R>=c0) ) { //these 6 terms are zero if  R>=c0 and R<=c1
    double exp1 = (2/(c1-R) - 1/(R-c0))*fabs(c0-c1);
/*
    cout << "exp1 = " << exp1 << "\n";

    grad_damping[3*mindist_atomA] = -1*exp(exp1)*(x1-x2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1);
    cout << "1st grad_damping[" << 3*mindist_atomA << "] = " << grad_damping[3*mindist_atomA] << "\n";
    grad_damping[3*mindist_atomA] /= (1+exp(exp1));
    cout << "2nd grad_damping[" << 3*mindist_atomA << "] = " << grad_damping[3*mindist_atomA] << "\n";
    grad_damping[3*mindist_atomA] /= (1+exp(exp1));
    cout << "3rd grad_damping[" << 3*mindist_atomA << "] = " << grad_damping[3*mindist_atomA] << "\n";
    grad_damping[3*mindist_atomA] /= R;
    cout << "4th grad_damping[" << 3*mindist_atomA << "] = " << grad_damping[3*mindist_atomA] << "\n";
*/    
    if ( fabs( exp(exp1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1) ) == std::numeric_limits<double>::infinity() ) { 
    // this logic has been added to avoid infinity/infinity ratio.. 
    // the idea is that if the numerator is infinity, the denominator which is always greater than the numerator must make the ratio zero (kdn)
      grad_damping[3*minA] = 0.0;
      grad_damping[3*minA+1] = 0.0;
      grad_damping[3*minA+2] = 0.0;
      grad_damping[3*Na+3*minB] = 0.0;
      grad_damping[3*Na+3*minB+1] = 0.0;
      grad_damping[3*Na+3*minB+2] = 0.0;
    }
    else {
      /*
      grad_temp[3*minA] = -1*exp(exp1)*(x1-x2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_temp[3*minA+1] = -1*exp(exp1)*(y1-y2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_temp[3*minA+2] = -1*exp(exp1)*(z1-z2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_temp[3*Na+3*minB] = -1*exp(exp1)*(x2-x1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_temp[3*Na+3*minB+1] = -1*exp(exp1)*(y2-y1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_temp[3*Na+3*minB+2] = -1*exp(exp1)*(z2-z1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);

      //applying rotation operator
      grad_damping[3*minA] = RotA(0,0)*grad_temp[3*minA]+RotA(0,1)*grad_temp[3*minA+1]+RotA(0,2)*grad_temp[3*minA+2];
      grad_damping[3*minA+1] = RotA(1,0)*grad_temp[3*minA]+RotA(1,1)*grad_temp[3*minA+1]+RotA(1,2)*grad_temp[3*minA+2];
      grad_damping[3*minA+2] = RotA(2,0)*grad_temp[3*minA]+RotA(2,1)*grad_temp[3*minA+1]+RotA(2,2)*grad_temp[3*minA+2];

      grad_damping[3*Na+3*minB] = RotB(0,0)*grad_temp[3*Na+3*minB]+RotB(0,1)*grad_temp[3*Na+3*minB+1]+RotB(0,2)*grad_temp[3*Na+3*minB+2];
      grad_damping[3*Na+3*minB+1] = RotB(1,0)*grad_temp[3*Na+3*minB]+RotB(1,1)*grad_temp[3*Na+3*minB+1]+RotB(1,2)*grad_temp[3*Na+3*minB+2];
      grad_damping[3*Na+3*minB+2] = RotB(2,0)*grad_temp[3*Na+3*minB]+RotB(2,1)*grad_temp[3*Na+3*minB+1]+RotB(2,2)*grad_temp[3*Na+3*minB+2];
      */
      grad_damping[3*minA] = -1*exp(exp1)*(x1-x2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_damping[3*minA+1] = -1*exp(exp1)*(y1-y2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_damping[3*minA+2] = -1*exp(exp1)*(z1-z2)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
      grad_damping[3*Na+3*minB] = -1*exp(exp1)*(x2-x1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
     grad_damping[3*Na+3*minB+1] = -1*exp(exp1)*(y2-y1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);
     grad_damping[3*Na+3*minB+2] = -1*exp(exp1)*(z2-z1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1)/((1+exp(exp1))*(1+exp(exp1))*R*AngToBohr);


    }
  }
  if ((c1==0) && (c0==0) ) {
    grad_damping.Set(); //no local truncation
  }
  //cout << "R = " << R <<"\n";
  //printf("d(%i,%i)\n",GetIndexA(),GetIndexB());
  //grad_temp.PrintGradient("grad_damping non rotated");
  //grad_damping.PrintGradient("grad_damping rotated");
  return grad_damping;
}

Matrix Dimer::GetSpatialDampingFunctionHessian(double c0, double c1) {

  //check that cutoff0 is not less than cutoff1
  if(c0 < c1){
    printf("ERROR::Dimer::GetSpatialDampingFunctionHessian(): cutoff0 is less than cutoff1. Alter cutoffs and rerun.\n");
    exit(0);
  }

  //initialize this gradient
  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  Matrix hess_damping(3*(Na+Nb),false);
  double R =  Separation; //distance between the two fragments
  double x1,y1,z1,x2,y2,z2; // xyz coords of the two atoms defining the minimum distance between the two monomers

  // Now first determine the two atoms which define the shortest distance between the two monomers
  // So, we build subroutines called MonA.MinDistAtomIndex_A(MonB) and MonA.MinDistAtomIndex_B(MonB)
  // in monomer.C and then we already used these to define mindist_atomA and mindist_atomB,
  // which are the indices of these two atoms in the respective monomers.
  // We now use these indices to build the grad_damping
  // Note that only 6 terms in this gradient will be non-zero.

  //cout << " minA = " << minA << "\n";
  //cout << " minB = " << minB << "\n";

  x1 =  MonA.GetAtom(minA).GetPosition(0);
  y1 =  MonA.GetAtom(minA).GetPosition(1);
  z1 =  MonA.GetAtom(minA).GetPosition(2);
  x2 =  MonB.GetAtom(minB).GetPosition(0);
  y2 =  MonB.GetAtom(minB).GetPosition(1);
  z2 =  MonB.GetAtom(minB).GetPosition(2);

  if ( !(R<=c1) && !(R>=c0) ) { //these 6 terms are zero if  R>=c0 and R<=c1
    double exp1 = (2/(c1-R) - 1/(R-c0))*fabs(c0-c1);

    if ( fabs( exp(exp1)*(1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1) ) == std::numeric_limits<double>::infinity() ) {
    // this logic has been added to avoid infinity/infinity ratio..
    // the idea is that if the numerator is infinity, the denominator which is always greater than the numerator must make the ratio zero (kdn)
      hess_damping.Set();

    }
    else {
      double t1 = exp(exp1);
      double tx2 = (x1-x2);
      double ty2 = (y1-y2);
      double tz2 = (z1-z2);
      double t3 = (1/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R))*fabs(c0-c1);
      double t4 = 1/((1+exp(exp1))*(1+exp(exp1))*R);

      double dt1_dx1 = exp(exp1)*tx2*t3/R ;
      double dt1_dy1 = exp(exp1)*ty2*t3/R ;  
      double dt1_dz1 = exp(exp1)*tz2*t3/R ;  
      double dt1_dx2 = exp(exp1)*(-tx2)*t3/R ;  
      double dt1_dy2 = exp(exp1)*(-ty2)*t3/R ;  
      double dt1_dz2 = exp(exp1)*(-tz2)*t3/R ;  

      double dtx2_dx1 = 1;
      double dtx2_dy1 = 0;
      double dtx2_dz1 = 0;
      double dtx2_dx2 = -1;
      double dtx2_dy2 = 0;
      double dtx2_dz2 = 0;

      double dty2_dx1 = 0;
      double dty2_dy1 = 1;
      double dty2_dz1 = 0;
      double dty2_dx2 = 0;
      double dty2_dy2 = -1;
      double dty2_dz2 = 0;

      double dtz2_dx1 = 0;
      double dtz2_dy1 = 0;
      double dtz2_dz1 = 1;
      double dtz2_dx2 = 0; 
      double dtz2_dy2 = 0;
      double dtz2_dz2 = -1;

      double dt3_dx1 = 2*tx2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;
      double dt3_dy1 = 2*ty2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;
      double dt3_dz1 = 2*tz2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;
      double dt3_dx2 = -2*tx2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;
      double dt3_dy2 = -2*ty2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;
      double dt3_dz2 = -2*tz2*(1/(c0-R)/(c0-R)/(c0-R) + 2/(c1-R)/(c1-R)/(c1-R))*fabs(c0-c1)/R ;

      double dt4_dx1 = tx2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;
      double dt4_dy1 = ty2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;
      double dt4_dz1 = tz2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;
      double dt4_dx2 = -tx2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;
      double dt4_dy2 = -ty2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;
      double dt4_dz2 = -tz2*(-1-t1-2*t1*t3*R) / (1+exp(exp1)) / (1+exp(exp1)) / (1+exp(exp1)) /R/R/R;

      // -1 comes from the -1 in grad_damping which none of t1,t2,t3,t4 accommodates
      //1st row
      hess_damping.Element(3*minA,3*minA) = - dt1_dx1*tx2*t3*t4 - dtx2_dx1*t1*t3*t4 - dt3_dx1*t1*tx2*t4 - dt4_dx1*t1*tx2*t3 ;
      hess_damping.Element(3*minA,3*minA+1) = - dt1_dy1*tx2*t3*t4 - dtx2_dy1*t1*t3*t4 - dt3_dy1*t1*tx2*t4 - dt4_dy1*t1*tx2*t3 ;
      hess_damping.Element(3*minA,3*minA+2) = - dt1_dz1*tx2*t3*t4 - dtx2_dz1*t1*t3*t4 - dt3_dz1*t1*tx2*t4 - dt4_dz1*t1*tx2*t3 ;
      hess_damping.Element(3*minA,3*Na+3*minB) = - dt1_dx2*tx2*t3*t4 - dtx2_dx2*t1*t3*t4 - dt3_dx2*t1*tx2*t4 - dt4_dx2*t1*tx2*t3 ;
      hess_damping.Element(3*minA,3*Na+3*minB+1) = - dt1_dy2*tx2*t3*t4 - dtx2_dy2*t1*t3*t4 - dt3_dy2*t1*tx2*t4 - dt4_dy2*t1*tx2*t3 ;
      hess_damping.Element(3*minA,3*Na+3*minB+2) = - dt1_dz2*tx2*t3*t4 - dtx2_dz2*t1*t3*t4 - dt3_dz2*t1*tx2*t4 - dt4_dz2*t1*tx2*t3 ;
       
      //2nd row
      hess_damping.Element(3*minA+1,3*minA) = - dt1_dx1*ty2*t3*t4 - dty2_dx1*t1*t3*t4 - dt3_dx1*t1*ty2*t4 - dt4_dx1*t1*ty2*t3 ;
      hess_damping.Element(3*minA+1,3*minA+1) = - dt1_dy1*ty2*t3*t4 - dty2_dy1*t1*t3*t4 - dt3_dy1*t1*ty2*t4 - dt4_dy1*t1*ty2*t3 ;
      hess_damping.Element(3*minA+1,3*minA+2) = - dt1_dz1*ty2*t3*t4 - dty2_dz1*t1*t3*t4 - dt3_dz1*t1*ty2*t4 - dt4_dz1*t1*ty2*t3 ;
      hess_damping.Element(3*minA+1,3*Na+3*minB) = - dt1_dx2*ty2*t3*t4 - dty2_dx2*t1*t3*t4 - dt3_dx2*t1*ty2*t4 - dt4_dx2*t1*ty2*t3 ;
      hess_damping.Element(3*minA+1,3*Na+3*minB+1) = - dt1_dy2*ty2*t3*t4 - dty2_dy2*t1*t3*t4 - dt3_dy2*t1*ty2*t4 - dt4_dy2*t1*ty2*t3 ;
      hess_damping.Element(3*minA+1,3*Na+3*minB+2) = - dt1_dz2*ty2*t3*t4 - dty2_dz2*t1*t3*t4 - dt3_dz2*t1*ty2*t4 - dt4_dz2*t1*ty2*t3 ;

      //3rd row
      hess_damping.Element(3*minA+2,3*minA) = - dt1_dx1*tz2*t3*t4 - dtz2_dx1*t1*t3*t4 - dt3_dx1*t1*tz2*t4 - dt4_dx1*t1*tz2*t3 ;
      hess_damping.Element(3*minA+2,3*minA+1) = - dt1_dy1*tz2*t3*t4 - dtz2_dy1*t1*t3*t4 - dt3_dy1*t1*tz2*t4 - dt4_dy1*t1*tz2*t3 ;
      hess_damping.Element(3*minA+2,3*minA+2) = - dt1_dz1*tz2*t3*t4 - dtz2_dz1*t1*t3*t4 - dt3_dz1*t1*tz2*t4 - dt4_dz1*t1*tz2*t3 ;
      hess_damping.Element(3*minA+2,3*Na+3*minB) = - dt1_dx2*tz2*t3*t4 - dtz2_dx2*t1*t3*t4 - dt3_dx2*t1*tz2*t4 - dt4_dx2*t1*tz2*t3 ;
      hess_damping.Element(3*minA+2,3*Na+3*minB+1) = - dt1_dy2*tz2*t3*t4 - dtz2_dy2*t1*t3*t4 - dt3_dy2*t1*tz2*t4 - dt4_dy2*t1*tz2*t3 ;
      hess_damping.Element(3*minA+2,3*Na+3*minB+2) = - dt1_dz2*tz2*t3*t4 - dtz2_dz2*t1*t3*t4 - dt3_dz2*t1*tz2*t4 - dt4_dz2*t1*tz2*t3 ;

      //4th row
      hess_damping.Element(3*Na+3*minB,3*minA) =  dt1_dx1*tx2*t3*t4 + dtx2_dx1*t1*t3*t4 + dt3_dx1*t1*tx2*t4 + dt4_dx1*t1*tx2*t3 ;
      hess_damping.Element(3*Na+3*minB,3*minA+1) =  dt1_dy1*tx2*t3*t4 + dtx2_dy1*t1*t3*t4 + dt3_dy1*t1*tx2*t4 + dt4_dy1*t1*tx2*t3 ;
      hess_damping.Element(3*Na+3*minB,3*minA+2) =  dt1_dz1*tx2*t3*t4 + dtx2_dz1*t1*t3*t4 + dt3_dz1*t1*tx2*t4 + dt4_dz1*t1*tx2*t3 ;
      hess_damping.Element(3*Na+3*minB,3*Na+3*minB) =  dt1_dx2*tx2*t3*t4 + dtx2_dx2*t1*t3*t4 + dt3_dx2*t1*tx2*t4 + dt4_dx2*t1*tx2*t3 ;
      hess_damping.Element(3*Na+3*minB,3*Na+3*minB+1) =  dt1_dy2*tx2*t3*t4 + dtx2_dy2*t1*t3*t4 + dt3_dy2*t1*tx2*t4 + dt4_dy2*t1*tx2*t3 ;
      hess_damping.Element(3*Na+3*minB,3*Na+3*minB+2) =  dt1_dz2*tx2*t3*t4 + dtx2_dz2*t1*t3*t4 + dt3_dz2*t1*tx2*t4 + dt4_dz2*t1*tx2*t3 ;

      //5th row
      hess_damping.Element(3*Na+3*minB+1,3*minA) =  dt1_dx1*ty2*t3*t4 + dty2_dx1*t1*t3*t4 + dt3_dx1*t1*ty2*t4 + dt4_dx1*t1*ty2*t3 ;
      hess_damping.Element(3*Na+3*minB+1,3*minA+1) =  dt1_dy1*ty2*t3*t4 + dty2_dy1*t1*t3*t4 + dt3_dy1*t1*ty2*t4 + dt4_dy1*t1*ty2*t3 ;
      hess_damping.Element(3*Na+3*minB+1,3*minA+2) =  dt1_dz1*ty2*t3*t4 + dty2_dz1*t1*t3*t4 + dt3_dz1*t1*ty2*t4 + dt4_dz1*t1*ty2*t3 ;
      hess_damping.Element(3*Na+3*minB+1,3*Na+3*minB) =  dt1_dx2*ty2*t3*t4 + dty2_dx2*t1*t3*t4 + dt3_dx2*t1*ty2*t4 + dt4_dx2*t1*ty2*t3 ;
      hess_damping.Element(3*Na+3*minB+1,3*Na+3*minB+1) =  dt1_dy2*ty2*t3*t4 + dty2_dy2*t1*t3*t4 + dt3_dy2*t1*ty2*t4 + dt4_dy2*t1*ty2*t3 ;
      hess_damping.Element(3*Na+3*minB+1,3*Na+3*minB+2) =  dt1_dz2*ty2*t3*t4 + dty2_dz2*t1*t3*t4 + dt3_dz2*t1*ty2*t4 + dt4_dz2*t1*ty2*t3 ;

      //6th row
      hess_damping.Element(3*Na+3*minB+2,3*minA) =  dt1_dx1*tz2*t3*t4 + dtz2_dx1*t1*t3*t4 + dt3_dx1*t1*tz2*t4 + dt4_dx1*t1*tz2*t3 ;
      hess_damping.Element(3*Na+3*minB+2,3*minA+1) =  dt1_dy1*tz2*t3*t4 + dtz2_dy1*t1*t3*t4 + dt3_dy1*t1*tz2*t4 + dt4_dy1*t1*tz2*t3 ;
      hess_damping.Element(3*Na+3*minB+2,3*minA+2) =  dt1_dz1*tz2*t3*t4 + dtz2_dz1*t1*t3*t4 + dt3_dz1*t1*tz2*t4 + dt4_dz1*t1*tz2*t3 ;
      hess_damping.Element(3*Na+3*minB+2,3*Na+3*minB) =  dt1_dx2*tz2*t3*t4 + dtz2_dx2*t1*t3*t4 + dt3_dx2*t1*tz2*t4 + dt4_dx2*t1*tz2*t3 ;
      hess_damping.Element(3*Na+3*minB+2,3*Na+3*minB+1) =  dt1_dy2*tz2*t3*t4 + dtz2_dy2*t1*t3*t4 + dt3_dy2*t1*tz2*t4 + dt4_dy2*t1*tz2*t3 ;
      hess_damping.Element(3*Na+3*minB+2,3*Na+3*minB+2) =  dt1_dz2*tz2*t3*t4 + dtz2_dz2*t1*t3*t4 + dt3_dz2*t1*tz2*t4 + dt4_dz2*t1*tz2*t3 ;

    }       

  }
  if ( (c1==0) && (c0==0) ) {
    hess_damping.Set(); //no local truncation
  }

  hess_damping.Scale(1/AngToBohr/AngToBohr); //correct unit is (1/bohr)^2
  //hess_damping.Print("final hess_damping");
  return hess_damping;         

}

int Dimer::GetNzero(double c0, double c1) {
  /*
   this function accounts for an increment in Nzero whenever the separation between monomers is greater than c0

   */
  int Nzero;
  double R =  Separation; //distance between the two fragments

  if ( (R >= c0) && (R != c1) ) {
    Nzero = 1;     
  }
  else {
    Nzero = 0;         
  }

  return Nzero;

}

int Dimer::GetNfull(double c0, double c1) {        
  /*
   this function accounts for an increment in Nfull whenever the separation between monomers is less than c1

   */
  int Nfull;
  double R =  Separation; //distance between the two fragments

  if (R <= c1) {
    Nfull = 1;
  }
  else {            
    Nfull = 0;        
  }                

  return Nfull;

}

int Dimer::GetNdamp(double c0, double c1) {        
  /*
   this function accounts for an increment in Ndamp whenever the separation between monomers is between c1 and c0

   */
  int Ndamp;
  double R =  Separation; //distance between the two fragments

  if ( (R < c0) && (R > c1) ) {
    Ndamp = 1;
  }
  else {            
    Ndamp = 0;        
  }                

  return Ndamp;

}


// Compute distances between atoms in the 2 fragments
void Dimer::ComputeIntermolecularDistances() {

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  printf(" Mon 1    Mon 2       Distance (Ang)\n");
  for (int i=0;i<Na;i++) {
    string atom1 = MonA.GetAtom(i).GetSymbol();
    for (int j=0;j<Nb;j++) {		
      string atom2 = MonB.GetAtom(j).GetSymbol();
      
      double dist = MonA.GetAtom(i).GetInterAtomicDistance( MonB.GetAtom(j) );
      printf("%2s(%3d)  %2s(%3d)       %7.4f\n",atom1.c_str(),i,
	     atom2.c_str(),j+Na,dist);
    }
  }
  printf("\n");
}

// Evaluate the force field energy using ab initio-derived parameters
double Dimer::ComputeAIFFEnergy() {
  double Eaiff = 0.0, E_2b_disp = 0.0, E_es_pol = 0.0, E_es = 0.0, E_pol = 0.0;

  printf("\nEvaluating ab initio force field contribution for dimer (%d,%d)\n",
	   indexA,indexB);

  bool do_es = Params::Parameters().DoAIFFElectrostatics();
  bool do_ind = Params::Parameters().DoAIFFInduction();
  bool do_2b_disp = Params::Parameters().DoAIFF2BodyDispersion();

  // Compute the electrostatic + induction/polarization energy
  if (do_es || do_ind) 
    E_es_pol = ComputeMultipoleInteractions(); // the real call


  // For printing, grab these energies
  E_es = GetMMElectrostaticEnergy(); 


  E_pol = GetMMInductionEnergy(); 

  // Compute the 2-body intermolecular dispersion energy
  if (do_2b_disp)
    E_2b_disp = ComputeTwoBodyDispersion();

  printf("  Electrostatic energy for dimer (%d,%d) = %8.4f kJ/mol\n",indexA,indexB,E_es*HartreesToKJpermole);
  printf("  Induction energy for dimer (%d,%d) =     %8.4f kJ/mol\n",indexA,indexB,E_pol*HartreesToKJpermole);
  printf("  Dispersion energy for dimer (%d,%d) =    %8.4f kJ/mol\n",indexA,indexB,E_2b_disp*HartreesToKJpermole);
  
  // Add up the individual contributions
  Eaiff = E_es_pol + E_2b_disp; 

  printf("  Total AIFF energy for dimer (%d,%d) =    %8.4f kJ/mol\n",indexA,indexB,Eaiff*HartreesToKJpermole);

  return Eaiff;
}


// Build Tab matrix for interaction between 2 monomer multipoles and
// compute their electrostatic & induction energies
double Dimer::ComputeMultipoleInteractions() {
  double energy = 0.0;

  if (Params::Parameters().PrintLevel() > 0) {
    printf("\nComputing permanent multipole interactions for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  // Variables to store total permanent electrostatic, induction & total energies
  double Ees = 0.0, Eind = 0.0, Etot = 0.0;

  // Step 1: Build the geometric interaction matrix Tab.  This matrix
  // contains only geometry information and gets used for evaluating
  // both the permanant and induced multipole electrostatics.  We
  // build two versions.  The first, Tabs, stores the standard
  // matrix. The second, DampedTabs, stores a version that is damped.
  // Damping is necessary for computing induction energies.

  double beta_damp = GetAIFFDampingFactor(); // damping factor for this dimer
  //double beta_damp = Params::Parameters().GetDampingFactor(); // no longer using

  // hack: different damping factors for water-water and glycine-water.
  // Assumes glycine is monomer 1.  Let that be set by damping factor.  Otherwise,
  // we assume beta_damp = 1.45 because it's water-water.
  /*
  if (indexA!=1 && indexB!=1) {
    beta_damp = 1.45;
    printf("*** HACK: Setting beta_damp = %.3f for dimer (%d,%d)\n",beta_damp,indexA,indexB);
  }
  */

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  Tabs = new Matrix[NatomsA*NatomsB];
  DampedTabs = new Matrix[NatomsA*NatomsB];
  Tabs_Init = 1;// Tabs has been initialized
  DampedTabs_Init = 1;// DampedTab has been initialized

  Matrix RotationMatrixA = MonA.GetRotationMatrix();
  Matrix RotationMatrixB = MonB.GetRotationMatrix();

  //RotationMatrixA.Print("Interaction MatrixA");
  //RotationMatrixB.Print("Interaction MatrixB");
  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {
  
      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      // Vector RotVecA(MonA.GetRotationVector());
      //double RotAngA = MonA.GetRotationAngle();
      //Vector RotVecB(MonB.GetRotationVector());
      //double RotAngB = MonB.GetRotationAngle();



      Matrix tmpTab(nQA,nQB);

      // Compute the undamped version of Tab.  "-999.0" is a dummy
      // parameter to turn damping off.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
      						       MonB.GetAtom(iB),RotationMatrixB, -999.0);
      //tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, -999.0);

      // Store this Tab matrix as one element in our list Tabs.
      Tabs[iB*NatomsA+iA].Initialize(tmpTab);  
      /*if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole); 
	tmpTab.Print("Tab - undamped");
	}*/

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
      					       MonB.GetAtom(iB),RotationMatrixB, beta_damp);
      //tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, beta_damp);

      DampedTabs[iB*NatomsA+iA].Initialize(tmpTab);  

      /* if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole);
	tmpTab.Print("Tab - damped");

      } */
    }
  }
  // Step 2: Evaluate permanent and induced electrostatic interactions

  // Step 2a: Permanent multipole contributions
  //          E_es = sum(tu) QA(t)*Tab(t,u)*QB(u) 
  for (int iA=0;iA<NatomsA;iA++) {
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
      
      // Evaluate the energy
      Vector tmp(QA.GetLength());
      // Tab*QB
      tmp = Tabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
      // QA*(Tab*QB)
      double Energy = QA.GetMoments().DotProduct(tmp);
      Ees += Energy;

      /* Debug printing: break energy down into separate contributions */
      if (Params::Parameters().PrintLevel() > 0) {
      // Create text label for each atom, and print the energy contribution
      char labelA[10], labelB[10];
      sprintf(labelA,"%s%d",MonA.GetAtom(iA).GetSymbol().c_str(),
	      MonA.GetAtom(iA).GetAtomIndex());
      sprintf(labelB,"%s%d",MonB.GetAtom(iB).GetSymbol().c_str(),
	      MonB.GetAtom(iB).GetAtomIndex());
      printf("Electrostatic contribution for %s(A)--%s(B) interaction = %f\n",
	     labelA, labelB,Energy*HartreesToKJpermole);
      }
      /* Debug printing: break energy down into separate contributions */
      if (Params::Parameters().PrintLevel() > 1) {
        Matrix Edecomp(5,5);
	Matrix Tab(Tabs[iB*NatomsA+iA]);
	
	for (int t=0;t<QA.GetLength();t++)
	  for (int u=0;u<QB.GetLength();u++) {
	    // chg-chg
	    if (t==0 && u==0) 
	      Edecomp(0,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-dip
	    else if (t==0 && u>=1 && u<=3) {
	      Edecomp(0,1) += QA(t)*Tab(t,u)*QB(u);
	    }
	    else if (u==0 && t>=1 && t<=3)
	      Edecomp(1,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-quad
	    else if (t==0 && u>=4 && u<=8)
	      Edecomp(0,2) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=4 && t<=8)
	      Edecomp(2,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-oct
	    else if (t==0 && u>=9 && u<=15)
	      Edecomp(0,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=9 && t<=15)
	      Edecomp(3,0) += QA(t)*Tab(t,u)*QB(u);
	    // chg-hexadec
	    else if (t==0 && u>=16 && u<=24)
	      Edecomp(0,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u==0 && t>=16 && t<=24)
	      Edecomp(4,0) += QA(t)*Tab(t,u)*QB(u);
	    
	    // dip-dip
	    else if (t>=1 && t<=3 && u>=1 && u<=3)
	      Edecomp(1,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-quad
	    else if (t>=1 && t<=3 && u>=4 && u<=8)
	      Edecomp(1,2) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=4 && t<=8)
	      Edecomp(2,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-oct
	    else if (t>=1 && t<=3 && u>=9 && u<=15)
	      Edecomp(1,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=9 && t<=15)
	      Edecomp(3,1) += QA(t)*Tab(t,u)*QB(u);
	    // dip-hexadec
	    else if (t>=1 && t<=3 && u>=16 && u<=24)
	      Edecomp(1,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=1 && u<=3 && t>=16 && t<=24)
	      Edecomp(4,1) += QA(t)*Tab(t,u)*QB(u);
	    
	    // quad-quad
	    else if (t>=4 && t<=8 && u>=4 && u<=8)
	      Edecomp(2,2) += QA(t)*Tab(t,u)*QB(u);
	    // quad-oct
	    else if (t>=4 && t<=8 && u>=9 && u<=15)
	      Edecomp(2,3) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=4 && u<=8 && t>=9 && t<=15)
	      Edecomp(3,2) += QA(t)*Tab(t,u)*QB(u);
	    // quad-hexadec
	    else if (t>=4 && t<=8 && u>=16 && u<=24)
	      Edecomp(2,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=4 && u<=8 && t>=16 && t<=24)
	      Edecomp(4,2) += QA(t)*Tab(t,u)*QB(u);
	    
	    // oct-oct
	    else if (t>=9 && t<=15 && u>=9 && u<=15)
	      Edecomp(3,3) += QA(t)*Tab(t,u)*QB(u);
	    // oct-hexa
	    else if (t>=9 && t<=15 && u>=16 && u<=24)
	      Edecomp(3,4) += QA(t)*Tab(t,u)*QB(u);
	    else if (u>=9 && u<=15 && t>=16 && t<=24)
	      Edecomp(4,3) += QA(t)*Tab(t,u)*QB(u);
	    
	    // hexadec-hexadec
	    else if (t>=16 && t<=24 && u>=16 && u<=24)
	      Edecomp(4,4) += QA(t)*Tab(t,u)*QB(u);
	  }
	Edecomp.Scale(HartreesToKJpermole);
	Edecomp.Print("Electrostatic energy decomposition"); printf("\n");
      }
    }
  }  

  // If not doing induction, exit here
  bool do_ind = Params::Parameters().DoAIFFInduction();
  if (! do_ind) {
    // Store the electrostaic energy contibution, in hartrees.
    E_Electrostatic_MM = Ees;   
    return Ees;
  }
    
  // Step 2b: Induced 2-body multipole contributions
  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction energy for dimer (%d,%d)\n",
	   indexA,indexB);
  }
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction.\n",
	 beta_damp);

  // Initialize a few variables
  bool iterate = true;
  int cycle = 0;
  double ind_conv = Params::Parameters().GetInductionConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  if (Params::Parameters().PrintLevel() > 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f\n", ind_conv,Econv);

  double Eind_old = 0.0;

  // Create lists that will store Induced Moments for each atom on A
  // and B and initialize them to the proper sizes based on the higher
  // rank object: polarizabilities or multipoles.  This is important:
  // CamCasp sometimes prints polarizabilities with high ranks & lots
  // of zeros, even if a lower rank would do.  We need two copies of
  // each list.
  Multipole *dQA = new Multipole[NatomsA];
  Multipole *old_dQA = new Multipole[NatomsA];
  for (int iA=0;iA<NatomsA;iA++) {
    int Rmom = MonA.GetAtom(iA).GetMultipoleMoments().GetRank();
    int Rpol = MonA.GetAtom(iA).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQA[iA].Initialize(max_rank);
    old_dQA[iA].Initialize(max_rank);
  }
  Multipole *dQB = new Multipole[NatomsB];
  Multipole *old_dQB = new Multipole[NatomsB];
  for (int iB=0;iB<NatomsB;iB++) {
    int Rmom = MonB.GetAtom(iB).GetMultipoleMoments().GetRank();
    int Rpol = MonB.GetAtom(iB).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQB[iB].Initialize(max_rank);
    old_dQB[iB].Initialize(max_rank);
  }

  Eind = 1000000.0;  // start with huge nonsense energy

  

  if (Params::Parameters().PrintLevel() > 0) {
    printf("--------------------------------------------------\n");
    printf(" Cycle       E(induction)       Change\n");
    printf("--------------------------------------------------\n");
  }
  // Begin the iterations to find self-consistent induction
  while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {
  
    cycle++;

    // Save data from previous cycle and reset the variables for this cycle
    Eind_old = Eind;
    Eind = 0.0;
    double EindA = 0.0,  EindB = 0.0;
    for (int iA=0;iA<NatomsA;iA++) {
      old_dQA[iA] = dQA[iA];
      dQA[iA].Set();
    }
    for (int iB=0;iB<NatomsB;iB++) {
      old_dQB[iB] = dQB[iB];
      dQB[iB].Set();
    }


    // Induce multipoles
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      // Grab some data we need
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Polarizability PolA(MonA.GetAtom(iA).GetPolarizability(),true);
      int NpolA = PolA.GetLength();
      int NmomA = QA.GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	// Grab some data we need
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	Polarizability PolB(MonB.GetAtom(iB).GetPolarizability(),true);
	int NpolB = PolB.GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	int dimT1 = Tab.GetRows();
	int dimT2 = Tab.GetCols();

	// Induce multipoles on monomer A due to monomer B
	// dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))
	for (int a=0;a<NpolA;a++) // loop over elements of dQA
	  for (int t=0;t<min(NpolA,dimT1);t++) 
	    for (int u=0;u<min(dimT2,NmomB);u++) {
	      dQA[iA](a) -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQB[iB](u));
	    }

	// Induce multipoles on monomer B due to monomer A
	// dQB(a) = dQB(a) - polB(a,t)*Tab(u,t)*(QA(u) + old_dQA(u))
	for (int a=0;a<NpolB;a++) // loop over elements of dQA
	  for (int t=0;t<min(NpolB,dimT2);t++) 
	    for (int u=0;u<min(dimT1,NmomA);u++) {
	      dQB[iB](a) -= PolB(a,t)*Tab(u,t)*(QA(u) + old_dQA[iA](u));
	    }
      }
    }
    
    // debug
    //dQA[0].GetMoments().Print("\ndQA");
    //dQB[0].GetMoments().Print("\n\ndQB");
    


    // Compute the induction energy using the induced multipoles determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);

	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);
	  }
      }
    }
    Eind = EindA + EindB;
  
    // Print out results for this cycle
    if (cycle==1 && Params::Parameters().PrintLevel() > 0) 
      printf(" %3d       %12.6f     *********** kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole);
    else if (Params::Parameters().PrintLevel() > 0) 
      printf(" %3d       %12.6f     %11.6f kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole);

    // Check convergence based on the energy change (in kJ/mol).
    if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies converged after %d iterations\n",cycle);
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
      iterate = false;
    }    

    if ( cycle == Params::Parameters().GetMaxPolarizationCycles()  && iterate == true ) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("--------------------------------------------------\n");
      printf("  Induction energies failed to converge after %d iterations\n",cycle);
      Params::Parameters().Warning();
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
    }
  }



  // Hack: looking at what happens if we don't use damping for final induction energy
  if (Params::Parameters().GetMaxOptCycles()==1) {
    printf("Recomputing dimer induction energy without damping\n");
    double EindA = 0.0;
    double EindB = 0.0;
    Eind = 0.0;
    // Compute the induction energy using the induced multipoles determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A
      
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();
      
      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      
	
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	
	
	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	
	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);
	  }
      }
    }
    Eind = EindA + EindB;
  }





  if (Params::Parameters().PrintLevel() > 1) {
    // Print out the final multipole moments
    printf(" *** Final induced multipoles ***\n");
    for (int iA=0;iA<NatomsA;iA++) {
      char label[50];
      sprintf(label,"Monomer 1, %s%d",MonA.GetAtom(iA).GetSymbol().c_str(),
	      MonA.GetAtom(iA).GetAtomIndex());
      string str = label;
      
      
      dQA[iA].Print(str);
    }
    for (int iB=0;iB<NatomsB;iB++) {
      char label[50];
      sprintf(label,"Monomer 2, %s%d",MonB.GetAtom(iB).GetSymbol().c_str(),
	      MonB.GetAtom(iB).GetAtomIndex());
      string str = label;
      dQB[iB].Print(str);
    }
  }

  // Store the energy contibutions, in hartrees.
  E_Electrostatic_MM = Ees;
  E_Induction_MM = Eind;
  Etot = Ees + Eind;

  if (Params::Parameters().PrintLevel() > 0) {
    // Print out a summary of the classical energy just calculated
    printf("--------------------------------------------------------\n");
    printf(" Total Classical Interaction Energy for Dimer (%d,%d)\n",
	   indexA,indexB);
    printf("--------------------------------------------------------\n");
    printf("      Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
    printf("      Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
    printf("    ---------------------------------------\n");
    printf("      Total         = %12.6f kJ/mol\n",Etot*HartreesToKJpermole);
    printf("--------------------------------------------------------\n");
  }

  delete [] dQA;
  delete [] old_dQA;
  delete [] dQB;
  delete [] old_dQB;

  return Etot;

}

// Build the nuclear gradient of the Tab matrix for interaction
// between 2 monomer multipoles.  Compute the forces due to
// electrostatics and induction.
Vector Dimer::ComputeMultipoleGradient() {

  Vector Grad(3*Natoms), GradRotated(3*Natoms);

  // Variables to store total permanent electrostatic, induction & total energies
  double Ees = 0.0, Eind = 0.0, Etot = 0.0;

  printf("\nEvaluating ab initio force field energy & gradient for dimer (%d,%d)\n",
	   indexA,indexB);

// Step 1: Build the gradient of the geometric interaction matrix Tab.  

  //yoni :added damp factor changed for multiple damps from Kaushik's code
  double beta_damp = GetAIFFDampingFactor(); // damping factor for this dimer
  //double beta_damp = Params::Parameters().GetDampingFactor(); //

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();


  // In case these have already been allocated, delete them.
  // If not, we are just deleting NULL, which is harmless.
  delete [] Tabs;
  delete [] DampedTabs;
  delete [] TabsGrad;
  delete [] DampedTabsGrad;


  Tabs = new Matrix[NatomsA*NatomsB];
  DampedTabs = new Matrix[NatomsA*NatomsB];

  TabsGrad = new Matrix[NatomsA*NatomsB];
  DampedTabsGrad = new Matrix[NatomsA*NatomsB];

  //TabsGrad and DampedTabsGrad have been initialized
  TabsGrad_Init = 1;
  DampedTabsGrad_Init = 1;

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Matrix RotationMatrixA = MonA.GetRotationMatrix();
      Matrix RotationMatrixB = MonB.GetRotationMatrix();
      //Vector RotVecA(MonA.GetRotationVector());
      //double RotAngA = MonA.GetRotationAngle();
      //Vector RotVecB(MonB.GetRotationVector());
      //double RotAngB = MonB.GetRotationAngle();
  
      /***Construct the Tab geometric interaction matrices and their
	  nuclear gradients ***/
      Matrix Tab(nQA,nQB);

      // Compute the undamped version of Tab.  "-999.0" is a dummy
      // parameter to turn damping off.
      Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
						       MonB.GetAtom(iB),
						       RotationMatrixB, -999.0);
      //Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
      //					       MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      Tabs[iB*NatomsA+iA].Initialize(Tab);  

      // Compute the damped version of Tab and store it.      
      Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
      					       MonB.GetAtom(iB),
      					       RotationMatrixB, beta_damp);

      // Tab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
      //					       MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, beta_damp);
      DampedTabs[iB*NatomsA+iA].Initialize(Tab);  


      /*** Now do the gradients of Tab and store them ***/
      Matrix Tab_grad(nQA*nQB,6);
      Tab = Tabs[iB*NatomsA+iA];

      // Compute the undamped version of dTab/dX.  "-999.0" is a dummy
      // parameter to turn damping off.
      //Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotVecA,
      //              RotAngA, MonB.GetAtom(iB), RotVecB, RotAngB, -999.0);
      Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotationMatrixA,
								 MonB.GetAtom(iB),
								 RotationMatrixB, -999.0);
      // Store this Tab matrix as one element in our list Tabs.
      TabsGrad[iB*NatomsA+iA].Initialize(Tab_grad); 

      // Compute the damped version of dTab/dX and store it.
      //Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotVecA,
      //               RotAngA, MonB.GetAtom(iB), RotVecB, RotAngB, beta_damp);
      Tab_grad = MonA.GetAtom(iA).BuildInteractionMatrixGradient(Tab,RotationMatrixA,
								 MonB.GetAtom(iB),
								 RotationMatrixB, beta_damp);
      DampedTabsGrad[iB*NatomsA+iA].Initialize(Tab_grad);  

    }

  // Step 2a: Compute the electrostatic interaction forces
  // Step 2a: Permanent multipole contributions
  //          dE/dX = sum(tu) QA(t)*dT(t,u)/dX*QB(u) 

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());

      int nQA = QA.GetLength();
      int nQB = QB.GetLength();

     
      // Evaluate the energy
      Vector tmp(nQA);
      // Tab*QB
      tmp = Tabs[iB*NatomsA+iA].MatrixTimesVector(QB.GetMoments());
      // QA*(Tab*QB)
      double Energy = QA.GetMoments().DotProduct(tmp);
      Ees += Energy;


      // And the gradient contribution:
      Matrix dTdX(TabsGrad[iB*NatomsA+iA]);
      /*
      if (iA==0 && iB==0) {
	//dTdX.Print("\ndTdX for O-O interaction");
	printf("\ndT/dX for O-O interaction\n");
	printf(" t   u      dAx        dAy        dAz        dBx        dBy        dBz\n");

	string types[25];
	types[0]="00";
	types[1]="1x"; types[2]="1y"; types[3]="1z";
	types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; 
	types[8]="22s";
	types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; 
	types[13]="32s"; 
	types[14]="33c"; types[15]="33s";
	types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; 
	types[20]="42s"; 
	types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";
	
	for (int u=0;u<nQB;u++) {
	  for (int t=0;t<nQA;t++) {
	    int tu = u*nQA + t;
	    printf("%3s %3s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		   types[t].c_str(),types[u].c_str(),dTdX(tu,0),dTdX(tu,1),
		   dTdX(tu,2),dTdX(tu,3), dTdX(tu,4),dTdX(tu,5));	    
	  }
	}
	printf("\n");
      }
      */

      // Create temporary storage for the gradient contribution
      // arising from a single pair of interacting atoms.
      Vector tmpGrad(6); // x1,y1,z1 and x2,y2,z2
      tmpGrad.Set();

      for (int X=0;X<6;X++) 
	for (int t=0;t<nQA;t++)
	  for (int u=0;u<nQB;u++) {
	    int tu = u*nQA + t; // define a compound index
	    tmpGrad[X] += QA(t)*dTdX(tu,X)*QB(u);
	  }

      // Incorporate this contribution into the full dimer gradient.
      NatomsA = MonA.GetNumberOfAtoms();
      for (int i=0;i<3;i++) {
	Grad[iA*3+i]  += tmpGrad[i]; // atom 1
	Grad[3*NatomsA+iB*3+i]  += tmpGrad[i+3]; // atom 2
	}
    }

  // We need to save this gradient contribution to build the full
  // cluster gradient later.
  Grad_Electrostatic = Grad;

  //printf("Dimer(%i,%i)\n",indexA,indexB);
  //Grad.PrintGradient("Electrostatics");

  Grad_Electrostatic_Init = 1;


  // Step 2b: Induced 2-body multipole gradient contributions
  if (Params::Parameters().PrintLevel() > 0) {
    printf("\n");
    printf("Computing self-consistent induction nuclear gradient for dimer (%d,%d)\n",indexA,indexB);
  }
  printf("  Tang-Toennies damping factor of %.3f applied to self-consisted induction gradient.\n",beta_damp);

  // Initialize a few variables
  bool iterate = true;
  int cycle = 0;

  double ind_conv = Params::Parameters().GetInductionConvergence();
  double ind_gradconv = Params::Parameters().GetInductionGradConvergence();
  double Econv = 1.0/pow(10,ind_conv);
  double Gconv = 1.0/pow(10,ind_gradconv);
  if (Params::Parameters().PrintLevel() > 0) 
    printf("  -----ind_conv = %12.6f, Econv = %12.6f, ind_gradconv = %12.6f, Gconv = %12.6f\n", ind_conv,Econv,ind_gradconv,Gconv );

  double Eind_old = 0.0;

  Vector Eind_grad(3*Natoms);
  Vector old_Eind_grad(3*Natoms);

  // Create lists that will store Induced Moments for each atom on A
  // and B and initialize them to the proper sizes based on the higher
  // rank object: polarizabilities or multipoles.  This is important:
  // CamCasp sometimes prints polarizabilities with high ranks & lots
  // of zeros, even if a lower rank would do.  We need two copies of
  // each list.
  Multipole *dQA = new Multipole[NatomsA];
  Multipole *old_dQA = new Multipole[NatomsA];

  // For the gradient, we have one multipole vector for each atom on a 
  // monomer.  Each of those vectors has 3*Natoms gradient components.
  // Overall, we have (3*Natoms)*NatomsA multipole vectors.
  Multipole *dQA_grad = new Multipole[3*Natoms*NatomsA];
  Multipole *old_dQA_grad = new Multipole[3*Natoms*NatomsA];

  for (int iA=0;iA<NatomsA;iA++) {
    int Rmom = MonA.GetAtom(iA).GetMultipoleMoments().GetRank();
    int Rpol = MonA.GetAtom(iA).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQA[iA].Initialize(max_rank);
    old_dQA[iA].Initialize(max_rank);

    for (int X=0;X<3*Natoms;X++) {
      dQA_grad[3*Natoms*iA + X].Initialize(max_rank);
      old_dQA_grad[3*Natoms*iA + X].Initialize(max_rank);
    }
  }

  Multipole *dQB = new Multipole[NatomsB];
  Multipole *old_dQB = new Multipole[NatomsB];

  Multipole *dQB_grad = new Multipole[3*Natoms*NatomsB];
  Multipole *old_dQB_grad = new Multipole[3*Natoms*NatomsB];


  for (int iB=0;iB<NatomsB;iB++) {
    int Rmom = MonB.GetAtom(iB).GetMultipoleMoments().GetRank();
    int Rpol = MonB.GetAtom(iB).GetPolarizability().GetRank();
    int max_rank = max(Rmom,Rpol);

    dQB[iB].Initialize(max_rank);
    old_dQB[iB].Initialize(max_rank);
    for (int X=0;X<3*Natoms;X++) {
      dQB_grad[3*Natoms*iB + X].Initialize(max_rank);
      old_dQB_grad[3*Natoms*iB + X].Initialize(max_rank);
    }
  }

  Eind = 1000000.0;  // start with huge nonsense energy

  if (Params::Parameters().PrintLevel() > 0) {
    printf("------------------------------------------------------------------------\n");
    printf(" Cycle     E(induction)       Change          |Grad|          Change   \n");
    printf("------------------------------------------------------------------------\n");
  }
  // Begin the iterations to find self-consistent induction
  while (iterate && (cycle < Params::Parameters().GetMaxPolarizationCycles()) ) {
  
    cycle++;

    // Save data from previous cycle and reset the variables for this cycle
    Eind_old = Eind;
    Eind = 0.0;
    double EindA = 0.0,  EindB = 0.0;
    old_Eind_grad = Eind_grad;
    Eind_grad.Set();

    for (int iA=0;iA<NatomsA;iA++) {
      old_dQA[iA] = dQA[iA];
      dQA[iA].Set();
      
      for (int X=0;X<3*Natoms;X++) {
	old_dQA_grad[3*Natoms*iA+X] = dQA_grad[3*Natoms*iA+X];
	dQA_grad[3*Natoms*iA+X].Set();
      }
    }
    for (int iB=0;iB<NatomsB;iB++) {
      old_dQB[iB] = dQB[iB];
      dQB[iB].Set();

      for (int X=0;X<3*Natoms;X++) {
	old_dQB_grad[3*Natoms*iB+X] = dQB_grad[3*Natoms*iB+X];
	dQB_grad[3*Natoms*iB+X].Set();
      }
    }

    // Induce multipoles
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      // Grab some data we need
      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      Polarizability PolA(MonA.GetAtom(iA).GetPolarizability(),true);
      int nQA = QA.GetLength();
      int NpolA = PolA.GetLength();
      int NmomA = QA.GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	// Grab some data we need
	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	Polarizability PolB(MonB.GetAtom(iB).GetPolarizability(),true);
	int nQB = QB.GetLength();
	int NpolB = PolB.GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);
	int dimT1 = Tab.GetRows();
	int dimT2 = Tab.GetCols();

	/*
	if (iA==0 && iB==0 &&) {
	  Matrix dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	  dTdX.Scale(HartreesToKJpermole);
	  printf("\nDamped dT/dX for O-O interaction\n");
	  printf(" t   u      dAx        dAy        dAz        dBx        dBy        dBz\n");

	  string types[25];
	  types[0]="00";
	  types[1]="1x"; types[2]="1y"; types[3]="1z";
	  types[4]="20"; types[5]="21c"; types[6]="21s"; types[7]="22c"; 
	  types[8]="22s";
	  types[9]="30"; types[10]="31c"; types[11]="31s"; types[12]="32c"; 
	  types[13]="32s"; 
	  types[14]="33c"; types[15]="33s";
	  types[16]="40"; types[17]="41c"; types[18]="41s"; types[19]="42c"; 
	  types[20]="42s"; 
	  types[21]="43c"; types[22]="43s"; types[23]="44c"; types[24]="44s";
	  
	  for (int u=0;u<nQB;u++) {
	    for (int t=0;t<nQA;t++) {
	      int tu = u*nQA + t;
	      printf("%3s %3s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		     types[t].c_str(),types[u].c_str(),dTdX(tu,0),dTdX(tu,1),
		     dTdX(tu,2),dTdX(tu,3), dTdX(tu,4),dTdX(tu,5));	    
	    }
	  }
	  printf("\n");
	}
	*/

	// Grab dTdX, and map it from (Qa*Qb,6) -> (Qa*Qb,3*Natoms);
	Matrix small_dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	int rows = small_dTdX.GetRows();
	Matrix dTdX(rows,3*Natoms);
	for (int irows=0;irows<rows;irows++) 
	  for (int X=0;X<3;X++) {
	    dTdX(irows,iA*3+X) = small_dTdX(irows,X);
	    dTdX(irows, 3*NatomsA+iB*3+X) = small_dTdX(irows,X+3); 
	  }

	// Induce multipoles on monomer A due to monomer B
	// dQA(a) = dQA(a) - polA(a,t)*Tab(t,u)*(QB(u)+ old_dQB(u))
	for (int a=0;a<NpolA;a++) // loop over elements of dQA
	  for (int t=0;t<min(NpolA,dimT1);t++) 
	    for (int u=0;u<min(dimT2,NmomB);u++) {
	      // Induced moments:
	      dQA[iA](a) -= PolA(a,t)*Tab(t,u)*(QB(u) + old_dQB[iB](u));

	      // Gradient of induced moments:
	      int tu = u*nQA + t; // define a compound index for (t,u)
	      for (int X=0;X<3*Natoms;X++) {
		dQA_grad[3*Natoms*iA + X](a) -= PolA(a,t)* 
		  (dTdX(tu,X)*(QB(u) + old_dQB[iB](u)) 
		   + Tab(t,u)*old_dQB_grad[3*Natoms*iB + X](u));
	      }
	    }

	// Induce multipoles on monomer B due to monomer A
	// dQB(a) = dQB(a) - polB(a,t)*Tab(u,t)*(QA(u) + old_dQA(u))
	for (int a=0;a<NpolB;a++) // loop over elements of dQA
	  for (int u=0;u<min(NpolB,dimT2);u++) 
	    for (int t=0;t<min(dimT1,NmomA);t++) {
	      // Induced moments:
	      dQB[iB](a) -= PolB(a,u)*Tab(t,u)*(QA(t) + old_dQA[iA](t));

	      // Gradient of induced moments:
	      int tu = u*nQA + t; // define a compound index for (t,u)
	      for (int X=0;X<3*Natoms;X++) {
		dQB_grad[3*Natoms*iB + X](a) -= PolB(a,u)* 
		  (dTdX(tu,X)*(QA(t) + old_dQA[iA](t)) 
		   + Tab(t,u)*old_dQA_grad[3*Natoms*iA + X](t));
	      }
	    }
      }
    }

    // Compute the induction energy/gradient using the induced multipoles
    // determined for this cycle
    for (int iA=0;iA<NatomsA;iA++)  {// loop over atoms on A

      Multipole QA(MonA.GetAtom(iA).GetMultipoleMoments());
      int nQA = QA.GetLength();
      int NmomA = QA.GetLength();
      int NpolA = MonA.GetAtom(iA).GetPolarizability().GetLength();

      for (int iB=0;iB<NatomsB;iB++)  {// loop over atoms on B      

	Multipole QB(MonB.GetAtom(iB).GetMultipoleMoments());
	int nQB = QB.GetLength();
	int NpolB = MonB.GetAtom(iB).GetPolarizability().GetLength();
	int NmomB = QB.GetLength();	

	Matrix Tab(DampedTabs[iB*NatomsA+iA]);

	// Grab dTdX, and map it from (Qa*Qb,6) -> (Qa*Qb,3*Natoms);
	Matrix small_dTdX(DampedTabsGrad[iB*NatomsA+iA]);
	int rows = small_dTdX.GetRows();
	Matrix dTdX(rows,3*Natoms);
	for (int X=0;X<3;X++) 
	  for (int irows=0;irows<rows;irows++) {
	    dTdX(irows,iA*3+X) = small_dTdX(irows,X);
	    dTdX(irows, 3*NatomsA+iB*3+X) = small_dTdX(irows,X+3); 
	  }
	    
	for (int t=0;t<NmomA;t++) 
	  for (int u=0;u<NmomB;u++) {
	    // Energy: 
	    EindA += 0.5*dQA[iA](t)*Tab(t,u)*QB(u);
	    EindB += 0.5*QA(t)*Tab(t,u)*dQB[iB](u);

	    // Gradient:
	    int tu = u*nQA + t; // define a compound index for (t,u)
	    for (int X=0;X<3*Natoms;X++) {
	      // Contribution due to dQA
	      Eind_grad[X] += 0.5*dQA_grad[3*Natoms*iA + X](t)*Tab(t,u)*QB(u)
		+ 0.5*dQA[iA](t)*dTdX(tu,X)*QB(u);

	      // Contribution due to dQB
	      Eind_grad[X] += 0.5*QA(t)*Tab(t,u)*dQB_grad[3*Natoms*iB + X](u)
		+ 0.5*QA(t)*dTdX(tu,X)*dQB[iB](u);
	    }
	  }
      }
    }
    Eind = EindA + EindB;
  
    double Gradnorm = Eind_grad.Norm();
    Vector tmp(Eind_grad);
    tmp -= old_Eind_grad;
    double dGrad = tmp.Norm();

    // Print out results for this cycle
    if (cycle==1 && Params::Parameters().PrintLevel() > 0) 
      printf(" %3d     %12.6f       *********    %12.6f      ********* kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole, Gradnorm*HartreesToKJpermole);
    else if (Params::Parameters().PrintLevel() > 0) 
      printf(" %3d     %12.6f     %11.6f    %12.6f    %11.6f kJ/mol\n",cycle,
	     Eind*HartreesToKJpermole,(Eind-Eind_old)*HartreesToKJpermole, 
	     Gradnorm*HartreesToKJpermole, dGrad*HartreesToKJpermole);

    // Check convergence based on the energy change (in kJ/mol).
    if (fabs(Eind - Eind_old)*HartreesToKJpermole < Econv 
	&& dGrad*HartreesToKJpermole < Gconv) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("------------------------------------------------------------------------\n");
      printf("  Induction energies & gradients converged after %d iterations\n",cycle);
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
      iterate = false;
    }    

    if ( cycle == Params::Parameters().GetMaxPolarizationCycles() && iterate == true ) {
      if (Params::Parameters().PrintLevel() > 0) 
	printf("------------------------------------------------------------------------\n");
      printf("  Induction energies & gradients failed to converge after %d iterations\n",
	     cycle);
      Params::Parameters().Warning();
      if (Params::Parameters().PrintLevel() > 0) 
	printf("\n");
    }
  }

  // Store the energy contibutions, in hartrees.
  E_Electrostatic_MM = Ees;
  E_Induction_MM = Eind;
  Etot = Ees + Eind;
  Energy_MM = Etot;

  if (Params::Parameters().PrintLevel() > 0) {
    // Print out a summary of the classical energy just calculated
    printf("--------------------------------------------------------\n");
    printf(" Total Classical Interaction Energy for Dimer (%d,%d)\n",
	   indexA,indexB);
    printf("--------------------------------------------------------\n");
    printf("      Electrostatic = %12.6f kJ/mol\n",Ees*HartreesToKJpermole);
    printf("      Induction     = %12.6f kJ/mol\n",Eind*HartreesToKJpermole);
    printf("    ---------------------------------------\n");
    printf("      Total         = %12.6f kJ/mol\n",Etot*HartreesToKJpermole);
    printf("--------------------------------------------------------\n");
  }

  // Return the gradient
  //Eind_grad.Scale(HartreesToKJpermole);
  //Eind_grad.PrintGradient("Induction gradient");
  //Eind_grad.Scale(1.0/HartreesToKJpermole); 

  //Grad.Scale(HartreesToKJpermole); 
  //Grad.PrintGradient("\nElectrostatic Gradient");
  //Grad.Scale(1.0/HartreesToKJpermole); 

  // Combine the permanent and induced gradient contributions
  Grad += Eind_grad;

  //Grad.PrintGradient("Dimer gradient\n");


  printf("Nonrotated grad of d(%i,%i)",indexA,indexB);
  Grad.PrintGradient(" ");
  
  
  //Rotating Gradient due to symmetry
  Matrix RotationA = MonA.GetRotationMatrix();
  Matrix RotationB = MonB.GetRotationMatrix();
  int SymA = MonA.GetSymmetricalMonomer();
  int SymB = MonB.GetSymmetricalMonomer();
  int NatomMonA = MonA.GetNumberOfAtoms();
  RotationA.Print("RotationA");
  RotationB.Print("RotationB");
  fflush(stdout);
  for(int i=0;i<NatomMonA;i++){
    for(int j=0;j<3;j++){
      GradRotated[3*i+j]= RotationA(j,0)*Grad[3*i]+RotationA(j,1)*Grad[3*i+1]+RotationA(j,2)*Grad[3*i+2];
    }
  }
  for(int i=NatomMonA;i<Natoms;i++){
    for(int j=0;j<3;j++){
      GradRotated[3*i+j]= RotationB(j,0)*Grad[3*i]+RotationB(j,1)*Grad[3*i+1]+RotationB(j,2)*Grad[3*i+2];
    }
  }
  
  // GradRotated.PrintGradient(" ");




  delete [] dQA;
  delete [] old_dQA;
  delete [] dQB;
  delete [] old_dQB;

  delete [] dQA_grad;
  delete [] old_dQA_grad;
  delete [] dQB_grad;
  delete [] old_dQB_grad;

  //return now only for testing purposes
  return GradRotated;
  //return Grad;
}

void Dimer::BuildDampedTabInteractionMatrices() {

  //printf("\nBuilding Damped Tab Interaction matrices for dimer (%d,%d)\n",
  //	   indexA,indexB);

  // Get the damping factor
  double beta_damp = GetAIFFDampingFactor();
  //double beta_damp = Params::Parameters().GetDampingFactor(); // no longer using

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  DampedTabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer
      Matrix RotationMatrixA = MonA.GetRotationMatrix();
      Matrix RotationMatrixB = MonB.GetRotationMatrix();
      //Vector RotVecA(MonA.GetRotationVector());
      //double RotAngA = MonA.GetRotationAngle();
      //Vector RotVecB(MonB.GetRotationVector());
      //double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
						       MonB.GetAtom(iB),
						       RotationMatrixB, beta_damp);
      //tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
      //						       MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, beta_damp);
      DampedTabs[iB*NatomsA+iA].Initialize(tmpTab);  

      // printf("D(%i,%i) Damped Tabs\n",indexA,indexB);
      //DampedTabs[iB*NatomsA+iA].PrintHessian("");
      //fflush(stdout);

      /* if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole);
	tmpTab.Print("Tab - damped");
      } */
    }
}

void Dimer::BuildTabInteractionMatrices() {

  //printf("\nBuilding Damped Tab Interaction matrices for dimer (%d,%d)\n",
  //	   indexA,indexB);

  // Get the damping factor
  double beta_damp = -999.0; // damping factor for this dimer

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  Tabs = new Matrix[NatomsA*NatomsB];

  for (int iA=0;iA<NatomsA;iA++)
    for (int iB=0;iB<NatomsB;iB++) {

      // create Tab[iA*NatomsA+iB] matrix
      int nQA = MonA.GetAtom(iA).GetMultipoleMoments().GetLength();
      int nQB = MonB.GetAtom(iB).GetMultipoleMoments().GetLength();
      
      // Grab info about local coordinate systems for each monomer      
      Matrix RotationMatrixA = MonA.GetRotationMatrix();
      Matrix RotationMatrixB = MonB.GetRotationMatrix();
      //Vector RotVecA(MonA.GetRotationVector());
      //double RotAngA = MonA.GetRotationAngle();
      //Vector RotVecB(MonB.GetRotationVector());
      //double RotAngB = MonB.GetRotationAngle();

      Matrix tmpTab(nQA,nQB);

      // Compute the damped version of Tab and store it.
      tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotationMatrixA,
						       MonB.GetAtom(iB),
						       RotationMatrixB, beta_damp);
      //tmpTab = MonA.GetAtom(iA).BuildInteractionMatrix(RotVecA,RotAngA,
      //					       MonB.GetAtom(iB),
      //					       RotVecB,RotAngB, beta_damp);
      Tabs[iB*NatomsA+iA].Initialize(tmpTab);  
      /* if (iA==0 && iB==0) {
	tmpTab.Scale(HartreesToKJpermole);
	tmpTab.Print("Tab - damped");
      } */
    }
}

Matrix Dimer::ReadHessian(string path, int type) {

  int NA = MonA.GetNumberOfAtoms();
  int NB = MonB.GetNumberOfAtoms();

  Matrix dimer_hess( 3*(NA+NB), 3*(NA+NB) ), hessA( 3*(NA+NB), 3*(NA+NB) ), hessB( 3*(NA+NB), 3*(NA+NB) );
  int cp_hess_index = 0;

  // Set up the filename, with the full path.  File is e.g. 'd1.2.out'
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d",indexA,indexB);               
  string filename = path + label;                  

  if (type == 2) // Tinker MM job  
    filename += ".freq";
  else if(Params::Parameters().GetQMType()==3) // G09
    filename += ".log";
  else // other
    filename += ".out";                           

  // Open the .freq/.out file
  ifstream infile;                 
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadHessian : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }
  
  if (type == 2) { // look in the tinker .freq file
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      // Search for the SCF hessian
      if ( line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):" ) {  
        //getline(infile,line); // throw away header line
        for (int i=0;i<3*(NA+NB);i++) {
          for (int j=0;j<3*(NA+NB);j++) {
            getline(infile,line);
            istringstream iss(line);
            string tmp1,tmp2;
            iss >> tmp1; // throw away the matrix row index
            iss >> tmp2; // throw away the matrix column index
            iss >> dimer_hess.Element(i,j); // Store the hessian elements
          }
        }
        break;
      }
    }
  }

  else if (Params::Parameters().GetQMType() == 1){//look in the qchem output
    
    if ( !(Params::Parameters().DoCounterpoise() ) ) {
      string line;
      while ( !infile.eof() ) {
	getline(infile,line);
	// Search for the SCF hessian
	if ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) {
	  for (int k=0;k<(3*(NA+NB)-3*(NA+NB)%6)/6;k++) {         
	    getline(infile,line); // throw away header line                  
	    //while (line.substr(0,15)!=" Gradient time:") {
	    
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=6*k;j<6*k+6;j++) {
		iss >> dimer_hess.Element(i,j); // Store the hessian elements
	      }                
	    }                
	  }
	  getline(infile,line); // check if this line signals end of hessian print or not
	  if (line.substr(0,15) != " Gradient time:") {             
	    for (int i=0;i<3*(NA+NB);i++) {             
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=(3*(NA+NB)-3*(NA+NB)%6);j<3*(NA+NB);j++) {
		iss >> dimer_hess.Element(i,j); // Store the hessian elements
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
    else if ( (Params::Parameters().DoCounterpoise() ) && (type==1) ) {
      string line;
      while ( !infile.eof() ) {      
	getline(infile,line);        
	// Search for the SCF hessian
	if ( ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) && (cp_hess_index == 0) ) {
	  for (int k=0;k<(3*(NA+NB)-3*(NA+NB)%6)/6;k++) {
	    getline(infile,line); // throw away header line
	    //while (line.substr(0,15)!=" Gradient time:") {
	    
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=6*k;j<6*k+6;j++) {
		iss >> dimer_hess.Element(i,j); // Store the hessian elements
	      }
	    }
	  }
	  getline(infile,line); // check if this line signals end of hessian print or not
	  if ( (line.substr(0,15) != " Gradient time:") || (line.substr(0,18) != " Writing REM_CC_EA") ) {           
	    for (int i=0;i<3*(NA+NB);i++) {                 
	      getline(infile,line);
	      istringstream iss(line);     
	      string tmp;           
	      iss >> tmp; // throw away the atom index
	      for (int j=(3*(NA+NB)-3*(NA+NB)%6);j<3*(NA+NB);j++) {
              iss >> dimer_hess.Element(i,j); // Store the hessian elements
	      }
	    }
	    cp_hess_index++;
	  }
	}
	if ( ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) && (cp_hess_index == 1) ) {
	  for (int k=0;k<(3*(NA+NB)-3*(NA+NB)%6)/6;k++) {
	    getline(infile,line); // throw away header line
	    //while (line.substr(0,15)!=" Gradient time:") {
	    
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=6*k;j<6*k+6;j++) {
		iss >> hessA.Element(i,j); // Store the hessian elements
	      }
	    }
	  }
	  getline(infile,line); // check if this line signals end of hessian print or not
	  if ( (line.substr(0,15) != " Gradient time:") || (line.substr(0,18) != " Writing REM_CC_EA") ) {
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=(3*(NA+NB)-3*(NA+NB)%6);j<3*(NA+NB);j++) {
		iss >> hessA.Element(i,j); // Store the hessian elements
	      }
	    }
	    cp_hess_index++;
	  }
	}
	if ( ( (line.substr(0,26)==" Hessian of the SCF Energy") || (line.substr(0,15) == " Final Hessian.") ) && (cp_hess_index == 2) ) {
	  for (int k=0;k<(3*(NA+NB)-3*(NA+NB)%6)/6;k++) {
	    getline(infile,line); // throw away header line
	    //while (line.substr(0,15)!=" Gradient time:") {
	    
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int j=6*k;j<6*k+6;j++) {
		iss >> hessB.Element(i,j); // Store the hessian elements
	      }
	    }
	  }
	  getline(infile,line); // check if this line signals end of hessian print or not
	  if ( (line.substr(0,15) != " Gradient time:") || (line.substr(0,18) != " Writing REM_CC_EA") ) {
	    for (int i=0;i<3*(NA+NB);i++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;    
	      iss >> tmp; // throw away the atom index
	      for (int j=(3*(NA+NB)-3*(NA+NB)%6);j<3*(NA+NB);j++) {
		iss >> hessB.Element(i,j); // Store the hessian elements
	      }
	    }
	    cp_hess_index++;
	  }
	}
	
	
      }                  
      if ( !(type == 2) && !(cp_hess_index == 3) ) {
	cout << cp_hess_index <<"\n";
	printf("Dimer::ReadHessian : At Least One Counterpoise Gradient Not Found or Incorrect '%s'\n",filename.c_str());
	exit(1);
      }
      
      //hessA.PrintHessian("hessA");
      //hessB.PrintHessian("hessB");
      //dimer_hess.PrintHessian("dimer_hess");
      //MonA.GetQMHessian().PrintHessian("MonA");
      //MonB.GetQMHessian().PrintHessian("MonB");                
      
      //1st block (top-left)
      for (int l=0;l<3*MonA.GetNumberOfAtoms();l++) {
	for (int m=0;m<3*MonA.GetNumberOfAtoms();m++) {
	  //dimer_hess.Element(l,m) = dimer_hess.Element(l,m) - hessA.Element(l,m) - hessB.Element(l,m) + MonA.GetQMHessian().Element(l,m);
	  dimer_hess.Element(l,m) = dimer_hess.Element(l,m) - hessA.Element(l,m) - hessB.Element(l,m);
	}
      }
      //2nd block (bottom-left)
      for (int l=0;l< 3*MonA.GetNumberOfAtoms();l++) {
	for (int m=3*MonA.GetNumberOfAtoms();m< (3*MonA.GetNumberOfAtoms()+3*MonB.GetNumberOfAtoms() );m++) {
	  dimer_hess.Element(l,m) = dimer_hess.Element(l,m) - hessA.Element(l,m) - hessB.Element(l,m);
	}
      }
      //3rd block (top-right)
      for (int l=3*MonA.GetNumberOfAtoms();l< (3*MonA.GetNumberOfAtoms()+3*MonB.GetNumberOfAtoms() );l++) {
	for (int m=0;m<3*MonA.GetNumberOfAtoms();m++) {
	  dimer_hess.Element(l,m) = dimer_hess.Element(l,m) - hessA.Element(l,m) - hessB.Element(l,m);
	}
      }
      //4th block (bottom-right)
      for (int l=0;l< 3*MonB.GetNumberOfAtoms() ;l++) {
	for (int m=0;m< 3*MonB.GetNumberOfAtoms();m++) {
	  int l1 = l + 3*MonA.GetNumberOfAtoms();
	  int m1 = m + 3*MonA.GetNumberOfAtoms();
	  //dimer_hess.Element(l1,m1) = dimer_hess.Element(l1,m1) - hessA.Element(l1,m1) - hessB.Element(l1,m1) + MonB.GetQMHessian().Element(l,m);
	  dimer_hess.Element(l1,m1) = dimer_hess.Element(l1,m1) - hessA.Element(l1,m1) - hessB.Element(l1,m1);
	}
      }
    
    }
  }
  else if(Params::Parameters().GetQMType() == 2){//look in the MolPro output
   
    string line;
    while ( !infile.eof() ) {
      getline(infile,line);
      if ( line.substr(0,16)==" Force Constants"){
	for (int k=0;k<=(3*(NA+NB)-3*(NA+NB)%5)/5;k++) {  
	  getline(infile,line); // throw away header line
	  for (int i=5*k;i<3*(NA+NB);i++){
	    getline(infile,line);
	    istringstream iss(line);
	    string tmp;
	    iss >> tmp; // throw away the atom index
	    for (int j=5*k;j<5*k+5;j++) {      
	      if(i >= j)
		iss >> dimer_hess.Element(i,j); // Store the hessian elements
	      
	    }//end of loop over j
	  }//end of loop over i
	}//end of loop over k
      }//end of if statement for string " Force Constants"
    }//end of file
    
      //coping lower triangle
    for(int i = 0;i<3*(NA+NB);i++)
      for(int j = i; j<3*(NA+NB);j++)
	  dimer_hess.Element(i,j) = dimer_hess.Element(j,i);
  }
  else if(Params::Parameters().GetQMType()==3) { // G09 Watit
    string line;
    while (!infile.eof()) {
      getline(infile,line);
      if(Params::Parameters().DoCounterpoise()) {
        if(line.find("Counterpoise corrected") != string::npos) {
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
                      ss >> dimer_hess.Element(i,j);
                      dimer_hess.Element(j,i)=dimer_hess.Element(i,j);
		      //HessQM(i,j)=dimer_hess.Element(i,j);
		      //HessQM(j,i)=HessQM(i,j);
                    }
		  }
                }
              }
            }
          }
	}
      } 
      else {
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
                  ss >> dimer_hess.Element(i,j);
                  dimer_hess.Element(j,i)=dimer_hess.Element(i,j);
                  //HessQM(i,j)=dimer_hess.Element(i,j);
                  //HessQM(j,i)=HessQM(i,j);
                }
              }
            }
          }
        }
      }
    }
    if(Params::Parameters().DoCounterpoise()) {
      for(int i=0;i<3*MonA.GetNumberOfAtoms();i++) {
        for(int j=0;j<3*MonA.GetNumberOfAtoms();j++) {
          dimer_hess.Element(i,j) = dimer_hess.Element(i,j) - MonA.GetQMHessian().Element(i,j);
        }
      }
      for(int i=0;i<3*MonB.GetNumberOfAtoms();i++) {
        int k = i + 3*MonA.GetNumberOfAtoms();
        for(int j=0;j<3*MonB.GetNumberOfAtoms();j++) {
          int l = j + 3*MonA.GetNumberOfAtoms();
          dimer_hess.Element(k,l) = dimer_hess.Element(k,l) - MonB.GetQMHessian().Element(i,j);
        }
      }
    }
  }
  else{
    printf("Dimer::ReadHessian() - QM type = %d is unknown\n",Params::Parameters().GetQMType());
    exit(1);
    
  }
  

  infile.close();  // Close the file
  /*
    for (int i=0;i<3*(NA+NB);i++) {
    for (int j=0;j<3*(NA+NB);j++) {
    cout << dimer_hess.Element(i,j) << "\t";
    }
    cout << "\n";
    }
  */

  //Including CCSD(T) correction, molpro only
  if(Params::Parameters().DoCCSDTCorrection() && type == 1){
    //printf("d%i.%i\n",indexA,indexB);
    //dimer_hess.PrintHessian("dimer_hess before CCSD(T)");
    Matrix CCSDT_hess = ReadEnergyFiniteDifferenceHessian(path,true);
    
    for(int i=0; i<3*(NA+NB);i++){
      for(int j=0; j<3*(NA+NB);j++){
	dimer_hess(i,j) += CCSDT_hess(i,j);
      }
    }
    //dimer_hess.PrintHessian("dimer hess");
    //exit(0);

  }


  //string title = "d" + char(indexA) + "." + char(indexB) + " hessian";
  //printf("d%i.%i\n",indexA,indexB);
  //dimer_hess.PrintHessian("hessian");
  return dimer_hess;               
}

Matrix Dimer::ReadFiniteDifferenceMolProHessian(string path){ 

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();

  Matrix dimer_hess(3*(Na+Nb),3*(Na+Nb));

  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = path + "/" + dimername; 
  

  if(Params::Parameters().DoEnergyFiniteDifferenceFreqs()){
   dimer_hess = ReadEnergyFiniteDifferenceHessian(path,false);
  }
  else{

    for(int i=0;i<3*(Na+Nb);i++){
      Vector GradPlus(3*(Na+Nb));
      Vector GradMinus(3*(Na+Nb));
      
      // The plus input file; 
      char num[10];
      sprintf(num,"%d",i);
      string filename = dimerpath + "/" + dimername + "+" + num + ".out";
      // Open the .out file
      ifstream infile;                 
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Dimer::ReadFiniteDifferenceMolProHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }   
      int cp_hess_index = 0;
      string line;
      
      if(Params::Parameters().DoCBS()){
	double basis1 = double(Params::Parameters().CBSBasis1());
	double basis2 = double(Params::Parameters().CBSBasis2());
	double ExpTerm = exp(1.54*(basis2-basis1));
	if(Params::Parameters().DoCounterpoise()){
	  Vector grad(3*(Na+Nb));
	  Vector gradA(3*(Na+Nb));
	  Vector gradB(3*(Na+Nb));
	  while ( !infile.eof() ) {
	    getline(infile,line);
	    if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	      cp_hess_index++;
	      getline(infile,line);//discard first line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;  
		
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double entry;
		  iss >> entry;
		  //HF for the dimer in the first basis 
		  if(cp_hess_index == 1){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
		  }
		  //MP2 for the dimer in the first basis 
		  else if(cp_hess_index == 2){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry; 
		  }
		  //HF with ghost MonA in the first basis
		  else if(cp_hess_index == 3){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry;
		  }
		  //Mp2 with ghost MonA in the first basis
		  else if(cp_hess_index == 4){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry; 
		  }
		  //HF with ghost MonB in the first basis
		  else if(cp_hess_index == 5){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		  //Mp2 with ghost MonB in the first basis
		  else if(cp_hess_index == 6){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry; 
		  }
		  //HF for dimer in the second basis
		  else if(cp_hess_index == 7){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
                }
		  //MP2 for dimer in the second basis
		  else if(cp_hess_index == 8){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
		  }
		  //HF with ghost MonA in the second basis
		  else if(cp_hess_index == 9){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry;
                }
		  //MP2 with ghost MonA in the second basis
		  else if(cp_hess_index == 10){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry;
		  }
		  //HF with ghost MonA in the second basis
		  else if(cp_hess_index == 11){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		  //MP2 with ghost MonA in the second basis
		  else if(cp_hess_index == 12){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		}
	      }
	    }
	  }
	  
	  
	  
	  if(cp_hess_index != 12){
	    printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		   filename.c_str());
	    printf("cp_hess_index = %i\n",cp_hess_index);
	    exit(1);
	  }
	  //grad.PrintGradient("grad");
	  //gradA.PrintGradient("gradA");
	  //gradB.PrintGradient("gradB");
	  for(int j=0;j<3*(Na+Nb);j++)
	    GradPlus[j] = grad[j] - gradA[j] - gradB[j];
	  
	  //GradPlus.PrintGradient("GradPlus");
	  
	}else{
	  while ( !infile.eof() ) {
	    getline(infile,line);
	    if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	      cp_hess_index++;
	      getline(infile,line);//discard first line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double entry;
		  iss >> entry;
		  //HF in the first basis
		  if(cp_hess_index == 1){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  //MP2 in the first basis  
		  else if(cp_hess_index == 2){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  //HF in the second basis
		  else if(cp_hess_index == 3){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  //MP2 in the second basis
		  else if(cp_hess_index == 4){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  GradPlus[3*j+k] += entry;
		}
		
	      }
	      
	    }	
	  }
	  if(cp_hess_index != 4){
	    printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		   filename.c_str());
	    printf("cp_hess_index = %i\n",cp_hess_index);
	    exit(1);
	  }
	}
	//GradPlus.Print("GradPlus"); 
      }
      else{
	
	while ( !infile.eof() ) {
	  getline(infile,line);
	  if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 0){
	    getline(infile,line);// throw away header line
	    
	    for (int j=0;j<Na+Nb;j++) {
	      getline(infile,line);
	      istringstream iss(line);
	      string tmp;
	      iss >> tmp; // throw away the atom index
	      for (int k=0;k<3;k++) {
		iss >> GradPlus[3*j+k];
		// Store the gradient elements
	      }
	    }
	    cp_hess_index++;
	    //break;
	  }
	  if(Params::Parameters().DoCounterpoise()){
	    if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 1){
	      getline(infile,line);// throw away header line
	      for (int j=0;j<Natoms;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double value;
		  iss >> value;
		  GradPlus[3*j+k] -= value;
		  // Store the gradient elements
		}
	      }
	      cp_hess_index++;
	      //break;
	    }
	    if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 2){
	      getline(infile,line);// throw away header line
	      for (int j=0;j<Natoms;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double value;
		  iss >> value;
		  GradPlus[3*j+k] -= value;
		  // Store the gradient elements
		}
	      }
	      cp_hess_index++;
	      //break;
	    }	
	  }
	}
	if((cp_hess_index != 3 && Params::Parameters().DoCounterpoise()) ||
	   (cp_hess_index != 1 && !Params::Parameters().DoCounterpoise()) ){
	  printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		 filename.c_str());
	  printf("cp_hess_index = %i\n",cp_hess_index);
	  exit(1);
	}
      }
      //GradPlus.PrintGradient("GradPlus");
      //printf("cp_hess_index = %i\n",cp_hess_index);
      infile.close();  // Close the file
      
      // The minus input file; 
      num[10];
      sprintf(num,"%d",i);
      filename = dimerpath + "/" + dimername + "-" + num + ".out";
      
      // Open the .out file        
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("Dimer::ReadFiniteDifferenceMolProHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }
      
      cp_hess_index = 0;
      
      
      if(Params::Parameters().DoCBS()){
	double basis1 = double(Params::Parameters().CBSBasis1());
	double basis2 = double(Params::Parameters().CBSBasis2());
	double ExpTerm = exp(1.54*(basis2-basis1));
	if(Params::Parameters().DoCounterpoise()){
	  Vector grad(3*(Na+Nb));
	  Vector gradA(3*(Na+Nb));
	  Vector gradB(3*(Na+Nb));
	  
	  while ( !infile.eof() ) {
	    getline(infile,line);
	    if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	      cp_hess_index++;
	      getline(infile,line);//discard first line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double entry;
		  iss >> entry;
		  //HF for the dimer in the first basis 
		  if(cp_hess_index == 1){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
		  }
		  //MP2 for the dimer in the first basis 
		  else if(cp_hess_index == 2){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry; 
		  }
		  //HF with ghost MonA in the first basis
		  else if(cp_hess_index == 3){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry;
		  }
		  //Mp2 with ghost MonA in the first basis
		  else if(cp_hess_index == 4){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry; 
		  }
		  //HF with ghost MonB in the first basis
		  else if(cp_hess_index == 5){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		  //Mp2 with ghost MonB in the first basis
		  else if(cp_hess_index == 6){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry; 
		  }
		  //HF for dimer in the second basis
		  else if(cp_hess_index == 7){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
		  }
		  //MP2 for dimer in the second basis
		  else if(cp_hess_index == 8){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    grad[3*j+k] += entry;
		  }
		  //HF with ghost MonA in the second basis
		  else if(cp_hess_index == 9){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradA[3*j+k] += entry;
		  }
		  //MP2 with ghost MonA in the second basis
		  else if(cp_hess_index == 10){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
                  gradA[3*j+k] += entry;
		  }
		  //HF with ghost MonA in the second basis
		  else if(cp_hess_index == 11){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		  //MP2 with ghost MonA in the second basis
		  else if(cp_hess_index == 12){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    gradB[3*j+k] += entry;
		  }
		}
	      }
	    }	  
	  }
	  
	  if(cp_hess_index != 12){
	    printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		   filename.c_str());
	    printf("cp_hess_index = %i\n",cp_hess_index);
	    exit(1);
	  }  
	  
	  
	  //grad.PrintGradient("grad");
	  //gradA.PrintGradient("gradA");
	  //gradB.PrintGradient("gradB");
	  for(int j=0;j<3*(Na+Nb);j++)
	    GradMinus[j] = grad[j] - gradA[j] - gradB[j];
	  
	  //GradMinus.PrintGradient("GradMinus");
	  
	}else{
	  while ( !infile.eof() ) {
	    getline(infile,line);
	    if(line==" Atom          dE/dx               dE/dy               dE/dz"){
	      cp_hess_index++;
	      getline(infile,line);//discard first line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double entry;
		  iss >> entry;
		  //HF in the first basis
		  if(cp_hess_index == 1){
		    entry *= -1/(ExpTerm-1)+pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  //MP2 in the first basis  
		  else if(cp_hess_index == 2){
		    entry *= -pow(basis1,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  else if(cp_hess_index == 3){
		    entry *= (ExpTerm)/(ExpTerm-1)-pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		  }
		  else if(cp_hess_index == 4){
		    entry *= pow(basis2,3)/(pow(basis2,3)-pow(basis1,3));
		    
		  }
		  GradMinus[3*j+k] += entry;
		}
		
	      }
	      
	    }	
	  }
	  if(cp_hess_index != 4){
	    printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		   filename.c_str());
	    printf("cp_hess_index = %i\n",cp_hess_index);
	    exit(1);
	  }  
	}
	
      }else{
	while ( !infile.eof() ) {
	  getline(infile,line);
	  if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 0){
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
	    cp_hess_index++;
	    //break;
	  }
	  
	  if(Params::Parameters().DoCounterpoise()){
	    if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 1){
	      getline(infile,line);// throw away header line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double value;
		  iss >> value;
		  GradMinus[3*j+k] -= value;
		// Store the gradient elements
		}
	      }
	      cp_hess_index++;
	      //break;
	    }
	    if(line==" Atom          dE/dx               dE/dy               dE/dz" && cp_hess_index == 2){
	      getline(infile,line);// throw away header line
	      for (int j=0;j<Na+Nb;j++) {
		getline(infile,line);
		istringstream iss(line);
		string tmp;
		iss >> tmp; // throw away the atom index
		for (int k=0;k<3;k++) {
		  double value;
		  iss >> value;
		  GradMinus[3*j+k] -= value;
		  // Store the gradient elements
		}
	      }
	      cp_hess_index++;
	      //break;
	    }
	    
	  }
	}
	
	if((cp_hess_index != 3 && Params::Parameters().DoCounterpoise()) ||
	   (cp_hess_index != 1 && !Params::Parameters().DoCounterpoise())){
	  printf("Dimer::ReadFiniteDifferenceHessian : At Least One Gradient Not Found or Incorrect '%s'\n",
		 filename.c_str());
	  printf("cp_hess_index = %i\n",cp_hess_index);
	  exit(1);
	}  
      }
      
      //GradMinus.Print("GradMinus");
      //printf("cp_hess_index = %i\n",cp_hess_index);
      infile.close();  // Close the file
      
      //take difference and divide by coordinate difference
      GradMinus -= GradPlus;
      GradMinus.Scale(1.0/(2*0.001*AngToBohr));
      //GradMinus.PrintGradient("dGrad");
      dimer_hess.SetColumnVector(GradMinus,i);
      
    }
    
  }
  
  //Including CCSD(T) correction, molpro only
  if(Params::Parameters().DoCCSDTCorrection()){
    //dimer_hess.PrintHessian("dimer_hess before CCSD(T)");
    Matrix CCSDT_hess = ReadEnergyFiniteDifferenceHessian(path,true);
    
    for(int i=0; i<3*(Na+Nb);i++)
      for(int j=0; j<3*(Na+Nb);j++)
	dimer_hess(i,j) += CCSDT_hess(i,j);
    //dimer_hess.PrintHessian("dimer hess");


  }
  //if(GetIndexA() == 2 && GetIndexB() == 44){

  //printf("d(%i,%i)\n",GetIndexA(),GetIndexB());
  //dimer_hess.PrintHessian("dimer hess");
    //}
  return dimer_hess;
}

//Implemented for molpro only
Matrix Dimer::ReadEnergyFiniteDifferenceHessian(string path,bool CCSDT) {

  int Na = MonA.GetNumberOfAtoms();
  int Nb = MonB.GetNumberOfAtoms();
  Matrix dimer_hess(3*(Na+Nb),3*(Na+Nb));


  string dimername =  "d"; 
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  dimername += label;
  string dimerpath = path + "/" + dimername;
  if(CCSDT)
    dimerpath += ".CCSDT";


  //printf("path = %s\n",path.c_str());
  for(int i=0;i<3*(Na+Nb);i++){
    for(int j=0;j<=i;j++){
      
      double energy1,energy2,energy3,energy4;

      char num[10];
      sprintf(num,"%d",i); 
      char num2[10];
      sprintf(num2,"%d",j);    
      
      //Second derivitive plus plus step
      string filename = dimerpath + "/" + dimername + "+" + num + "+" + num2 + ".out";
      //printf("filename = %s\n",filename.c_str());

      // Open the energy file
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      // Search for MolPro energy for finite difference step
      string line;
      int number_energies = 0;
      while ( !infile.eof() ) {
	getline(infile,line);
	string match = line.substr(0,15);
	

	// Search for final MolPro energy
	if ( match==" SETTING EDIMER" ) { 
	  number_energies++;
	  istringstream iss(line);
	  string tmp;
	  for (int k=0;k<3;k++)
	    iss >> tmp; // throw away text 
	  iss >> energy1; // read energy
	}
      }
      infile.close();
      if(number_energies != 1){
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Incorrect number of energies found in %s\n",filename.c_str());
	printf("number of energies = %i",number_energies);
        exit(1);
      }

      //Second derivitive plus minus step
      filename = dimerpath + "/" + dimername + "+" + num + "-" + num2 + ".out";
      //printf("filename = %s\n",filename.c_str());

      // Open the energy file
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      // Search for MolPro energy for finite difference step
      number_energies = 0;
      while ( !infile.eof() ) {
	getline(infile,line);
	string match = line.substr(0,15);
	

	// Search for final MolPro energy
	if ( match==" SETTING EDIMER" ) { 
	  number_energies++;
	  istringstream iss(line);
	  string tmp;
	  for (int k=0;k<3;k++)
	    iss >> tmp; // throw away text 
	  iss >> energy2; // read energy
	}
      }
      infile.close();
      if(number_energies != 1){
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Incorrect number of energies found in %s\n",filename.c_str());
	printf("number of energies = %i",number_energies);
        exit(1);
      }

     //Second derivitive minus plus step
      filename = dimerpath + "/" + dimername + "-" + num + "+" + num2 + ".out";
      //printf("filename = %s\n",filename.c_str());

      // Open the energy file
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      // Search for MolPro energy for finite difference step
      number_energies = 0;
      while ( !infile.eof() ) {
	getline(infile,line);
	string match = line.substr(0,15);
	

	// Search for final MolPro energy
	if ( match==" SETTING EDIMER" ) { 
	  number_energies++;
	  istringstream iss(line);
	  string tmp;
	  for (int k=0;k<3;k++)
	    iss >> tmp; // throw away text 
	  iss >> energy3; // read energy
	}
      }
      infile.close();
      if(number_energies != 1){
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Incorrect number of energies found in %s\n",filename.c_str());
	printf("number of energies = %i",number_energies);
        exit(1);
      }

     //Second derivitive minus  step
      filename = dimerpath + "/" + dimername + "-" + num + "-" + num2 + ".out";
      //printf("filename = %s\n",filename.c_str());

      // Open the energy file
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Cannot open file '%s'\n",
	       filename.c_str());
	exit(1);
      }

      // Search for MolPro energy for finite difference step
      number_energies = 0;
      while ( !infile.eof() ) {
	getline(infile,line);
	string match = line.substr(0,15);
	

	// Search for final MolPro energy
	if ( match==" SETTING EDIMER" ) { 
	  number_energies++;
	  istringstream iss(line);
	  string tmp;
	  for (int k=0;k<3;k++)
	    iss >> tmp; // throw away text 
	  iss >> energy4; // read energy
	}
      }
      infile.close();
      if(number_energies != 1){
	printf("ERROR:Dimer::ReadEnergyFiniteDifferenceHessian : Incorrect number of energies found in %s\n",filename.c_str());
	printf("number of energies = %i",number_energies);
        exit(1);
      }


      //Get hessian from energies
      double FirstDev1 = (energy1 - energy2)/(0.002*AngToBohr);
      double FirstDev2 = (energy3 - energy4)/(0.002*AngToBohr);
      dimer_hess(i,j) = (FirstDev1 - FirstDev2)/(.002*AngToBohr);
      //printf("FirstDev1 = %15.9f FirstDev2 = %15.9f\n",FirstDev1,FirstDev2);
      //printf("dimer(%i,%i) = %15.9f\n",i,j,dimer_hess(i,j));
    }
  }

  //Get the upper triangle of the Hessian
  for(int i=0;i<3*(Na+Nb);i++){
    for(int j=i+1;j<3*(Na+Nb);j++){
      //printf("i = %i j = %i\n",i,j);
      dimer_hess(i,j) = dimer_hess(j,i);
    }
  }

  //dimer_hess.PrintHessian("dimer_hess CCSD(T)");
  return dimer_hess;

}


void Dimer::SetMMHessian() {  

  if (Params::Parameters().GetMMType() == 1) { // Tinker
    string path = Params::Parameters().GetHessianMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Hess_MM = ReadHessian(path,2);
  }
  else if (Params::Parameters().GetMMType() == 2) { // AIFF
    cout << " AIFF hessian not yet implemented \n";
    exit(1);
//    Hess_MM = ComputeMultipoleHessian(); // ??               
  }
  else if (Params::Parameters().GetMMType() == 3) { //QChem 
    string path = Params::Parameters().GetHessianMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    Hess_MM = ReadHessian(path,1);
  }
  else {
    printf("Dimer::SetMMHessian() - MM type = %d is unknown\n",
           Params::Parameters().GetMMType());
    exit(1);
  }

  MM_Hess_Init = 1;  
}

// Get the QM Hessian, wrapper routine
void Dimer::SetQMHessian() { 
  string path = Params::Parameters().GetHessianQMPath();
  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


  //MolPro hessian found by finite difference
  if(Params::Parameters().GetQMType() == 2)
    Hess_QM = ReadFiniteDifferenceMolProHessian(path);
  else
    Hess_QM = ReadHessian(path,1);
  QM_Hess_Init = 1;
}

// This function uses a list of predetermined c6 coeffs.  You should
// generally use the ab initio dispersion model found in
// Dimer::ComputeTwoBodyDispersion() instead.
double Dimer::EstimateTwoBodyDispersion() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      //double Rij = MonA.GetAtom(i).GetInterAtomicDistance(MonB.GetAtom(j))*AngToBohr;
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
      //printf("Rij = %f\n",Rij);
      
      // Get the C6 coefficient for this pair of atoms
      double C6ij = AtomI.EstimateC6Coefficient(AtomJ);

      // Get the van der Waals diameters for the two atoms
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
      
      // Get Tang-Toennies-type 2-body damping factor
      // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
      double beta = -0.33*(Di+Dj) + 4.39;
      double damping = AtomI.TangToenniesDampingFactor(6,beta,Rij);

      // E(2b) += - damping * C6ij/Rij^6 
      E2b_disp -= damping * C6ij / (pow(Rij,6));
    }
  }
  //printf("Approximate 2-body dispersion for dimer (%d,%d): %f in kcal/mol\n",indexA, indexB,
  // E2b_disp*HartreesToKcalpermole);
  return E2b_disp;
}

// Does C6 and C8 2-body dispersion.  Can do C10 as well, but that
// term is small and hard to get right, so we disable it.  The
// dispersion coefficients are calculated via ab initio calculation of
// the frequency-dependent polarizabilities and Casimir-Polder
// integration.
double Dimer::ComputeTwoBodyDispersion() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp = 0.0;
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Compute the C6 & C8 dispersion terms via Casimir-Polder integration
      double C6ij = AtomI.CasimirC6Coefficient(AtomJ);
      double C8ij = AtomI.CasimirC8Coefficient(AtomJ);

     
      // We used to get C10 as well, but that contributes little
      // and is hard to get right.
      //double C10ij = AtomI.CasimirC10Coefficient(AtomJ);

      //printf("Dispersion coefficients for atom pair (%d,%d): C6 = %f, C8=%f\n",i,j,C6ij,C8ij);


      /*
      // Get the C6 coefficient for this pair of atoms

      double C6i = AtomI.GetC6Coefficient();
      double C6j = AtomJ.GetC6Coefficient();
      double Poli = AtomI.GetIsotropicDipolePolarizability();
      double Polj = AtomJ.GetIsotropicDipolePolarizability();

      // Use Tang (Phys Rev 177, 108-114, 1969) combination rule for
      // C6 coefficients.  Or see Stone's Theory of IMF, eqn 4.3.14
      // (pg 61 in the 2002 edition).
      //double C6ij = 2*Poli*Polj*C6i*C6j / (Poli*Poli*C6j + Polj*Polj*C6i);
      
      // Get the C8 coefficient for this pair of atoms
      double C8i = AtomI.GetC8Coefficient();
      double C8j = AtomJ.GetC8Coefficient();
      //printf("C8i=%f, C8j=%f\n",C8i,C8j);
      // Just use crude combination rule, Cab = sqrt(Caa*Cbb)
      // handle signs, in case coefficient is negative
      double sign = 1.0;
      if (C8i < 0.0) sign *= -1.0;
      if (C8j < 0.0) sign *= -1.0;
      double C8ij = sqrt( fabs(C8i) * fabs(C8j) );

a      // Get the C10 coefficient for this pair of atoms
      double C10i = AtomI.GetC10Coefficient();
      double C10j = AtomJ.GetC10Coefficient();
      //printf("C10i=%f, C10j=%f\n",C10i,C10j);
      // Just use crude combination rule, Cab = sqrt(Caa*Cbb)
      // handle signs, in case coefficient is negative
      sign = 1.0;
      if (C10i < 0.0) sign *= -1.0;
      if (C10j < 0.0) sign *= -1.0;
      double C10ij = sqrt( fabs(C10i) * fabs(C10j) );
      //double C10ijnew = AtomI.CasimirC10Coefficient(AtomJ);

      // printf("C6ij  = %f (based on %f and %f)\n",C6ij,C6i,C6j);
      // printf("C8ij  = %f (based on %f and %f)\n",C8ij,C8i,C8j);
      // printf("C10ij = %f (based on %f and %f)\n",C10ij,C10i,C10j);
      */
      
      // Get Tang-Toennies-type 2-body damping factor
      // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
      // Frankly, it doesn't matter.  With typical values of the following,
      // the damping is already less than 0.1% by 4-5 Angstroms, which is
      // well-beyond our quantum cutoff.  
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
      double beta = -0.33*(Di+Dj) + 4.39;
      double damping6 = AtomI.TangToenniesDampingFactor(6,beta,Rij);
      double damping8 = AtomI.TangToenniesDampingFactor(8,beta,Rij);
      //double damping10 = AtomI.TangToenniesDampingFactor(10,beta,Rij);
      //printf("beta = %f, Rij = %f, damping 6/8/10 = %f, %f, %f\n",beta,Rij,damping6,damping8,damping10);

      // E(2b) += - damping6 * C6ij/Rij^6 - damping8 * C8ij/Rij^8 - damping10 * C10ij/Rij^10
      double d6 = -1.0*damping6 * C6ij / (pow(Rij,6));
      double d8 = -1.0*damping8 * C8ij / (pow(Rij,8));
      //double d10 = -1.0*damping10 * C10ij / (pow(Rij,10));

      //E2b_disp += d6 + d8 + d10;
      E2b_disp += d6 + d8;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp*HartreesToKJpermole);
  }

  // Store the energy, in hartrees
  E_2b_disp_MM = E2b_disp;

  return E2b_disp;

}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C6.
double Dimer::ComputeTwoBodyDispersionC6() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C6 = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;

      // Get the C6 coefficient for this pair of atoms
      double C6ij = AtomI.CasimirC6Coefficient(AtomJ);    

      // Get Tang-Toennies-type 2-body damping factor
      // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
      // Frankly, it doesn't matter.  With typical values of the following,
      // the damping is already less than 0.1% by 4-5 Angstroms, which is
      // well-beyond our quantum cutoff.  
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");  
      double beta = -0.33*(Di+Dj) + 4.39;
      double damping6 = AtomI.TangToenniesDampingFactor(6,beta,Rij);

      // E(2b) += - damping6 * C6ij/Rij^6 - damping8 * C8ij/Rij^8 - damping10 * C10ij/Rij^10
      double d6 = -1.0*damping6 * C6ij / (pow(Rij,6));
      E2b_disp_C6 += d6;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C6 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C6*HartreesToKJpermole);
  }
  return E2b_disp_C6;
}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C8.
double Dimer::ComputeTwoBodyDispersionC8() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C8 = 0.0;
 
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
      
      // Get the C8 coefficient for this pair of atoms
      double C8ij = AtomI.CasimirC8Coefficient(AtomJ);
      
      // Get Tang-Toennies-type 2-body damping factor
      // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
      // Frankly, it doesn't matter.  With typical values of the following,
      // the damping is already less than 0.1% by 4-5 Angstroms, which is
      // well-beyond our quantum cutoff.  
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
      double beta = -0.33*(Di+Dj) + 4.39;
      double damping8 = AtomI.TangToenniesDampingFactor(8,beta,Rij);

      // E(2b) += - damping8 * C8ij/Rij^8 
      double d8 = -1.0*damping8 * C8ij / (pow(Rij,8));
      E2b_disp_C8 += d8 ;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C8 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C8*HartreesToKJpermole);
  }
  return E2b_disp_C8;
}

// Version of Dimer::ComputeTwoBodyDispersion() that does only C10.
double Dimer::ComputeTwoBodyDispersionC10() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E2b_disp_C10 = 0.0;
  
  for (int i=0;i<NatomsA; i++) {
    Atom AtomI = MonA.GetAtom(i);

    for (int j=0; j<NatomsB; j++) {
      Atom AtomJ = MonB.GetAtom(j);
      //printf("Dimer (%d,%d)... Atom pair (%d,%d)\n",indexA, indexB, i, j);

      // Get interatomic distance, in bohr
      double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
  
      // Get the C10 coefficient for this pair of atoms
      double C10ij = AtomI.CasimirC10Coefficient(AtomJ);

      // Get Tang-Toennies-type 2-body damping factor
      // Use empirical expression for beta: J Chem Phys 132, 234109 (2010).
      // Frankly, it doesn't matter.  With typical values of the following,
      // the damping is already less than 0.1% by 4-5 Angstroms, which is
      // well-beyond our quantum cutoff.  
      double Di = AtomI.LookupAtomicDispersionParameter("Rvdw");
      double Dj = AtomJ.LookupAtomicDispersionParameter("Rvdw");
      double beta = -0.33*(Di+Dj) + 4.39;
      double damping10 = AtomI.TangToenniesDampingFactor(10,beta,Rij);

      // E(2b) += - damping10 * C10ij/Rij^10
      double d10 = -1.0*damping10 * C10ij / (pow(Rij,10));
      E2b_disp_C10 += d10;
    }
  }
  if ( Params::Parameters().PrintLevel() > 0) {
    printf("  2-body dispersion from C10 for dimer (%d,%d): %f kJ/mol\n",indexA, indexB,
	   E2b_disp_C10*HartreesToKJpermole);
  }
  return E2b_disp_C10;
}

// Does Axilrod-Teller-Muto 3-body dispersion, albeit with approximate
// tabulated C9 coefficients.
double Dimer::EstimateThreeBodyDispersion() {

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E3b_disp = 0.0;

  // Do as two cases: 

  /* Case 1: 2 atoms from MonA, 1 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=i+1;j<NatomsA;j++) {
      Atom AtomJ = MonA.GetAtom(j);
      for (int k=0;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;
	/*
	printf("Rij = %f, Rik = %f, Rjk = %f\n",Rij,Rik,Rjk);
	printf("PhiI = %f, PhiJ = %f, PhiK = %f,  Sum = %f\n",
	       acos(cosPhiI)*RadiansToDegrees,
	       acos(cosPhiJ)*RadiansToDegrees,
	       acos(cosPhiK)*RadiansToDegrees,
	       (acos(cosPhiI)+acos(cosPhiJ)+acos(cosPhiK))*RadiansToDegrees); 
	*/


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
	
	// ATM 3-body dispersion contribution
	//printf("Daab(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }

  /* Case 2: 1 atom from MonA, 2 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=0;j<NatomsB;j++) {
      Atom AtomJ = MonB.GetAtom(j);
      for (int k=j+1;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;
	/*
	printf("Rij = %f, Rik = %f, Rjk = %f\n",Rij,Rik,Rjk);
	printf("PhiI = %f, PhiJ = %f, PhiK = %f,  Sum = %f\n",
	       acos(cosPhiI)*RadiansToDegrees,
	       acos(cosPhiJ)*RadiansToDegrees,
	       acos(cosPhiK)*RadiansToDegrees,
	       (acos(cosPhiI)+acos(cosPhiJ)+acos(cosPhiK))*RadiansToDegrees); 
	*/

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
	
	// ATM 3-body dispersion contribution
	//printf("Dabb(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }
  return E3b_disp;

}


// Does Axilrod-Teller-Muto 3-body dispersion with ab initio C9
// coefficients computed via Casimir-Polder integration of the
// frequency-dependent polarizabilities.  The contributions here
// involve 2 atoms on one monomer and 1 on the other.
double Dimer::ComputeThreeBodyDispersion() {
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  double E3b_disp = 0.0;

  // Do as two cases: 

  /* Case 1: 2 atoms from MonA, 1 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=i+1;j<NatomsA;j++) {
      Atom AtomJ = MonA.GetAtom(j);
      for (int k=0;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;

	// Get the C9 coefficient
	double C9ijk = AtomI.CasimirC9Coefficient(AtomJ, AtomK);
	
	// Get the van der Waals diameters
	// Note: these are empirical and should be replaced with
	// something more rigorous.  In this crude version, the
	// DispersionAtomType parameters should have already been set
	// when computing the 2-body dispersion.
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
	//printf("Daab(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }

  /* Case 2: 1 atom from MonA, 2 from MonB */

  // Loop over individual atoms in the trimers
  for (int i=0;i<NatomsA;i++) {
    Atom AtomI = MonA.GetAtom(i);
    for (int j=0;j<NatomsB;j++) {
      Atom AtomJ = MonB.GetAtom(j);
      for (int k=j+1;k<NatomsB;k++) {
	Atom AtomK = MonB.GetAtom(k);

	// Get geometrical parameters... in bohr and radians 
	double Rij = AtomI.GetInterAtomicDistance(AtomJ)*AngToBohr;
	double Rik = AtomI.GetInterAtomicDistance(AtomK)*AngToBohr;
	double Rjk = AtomJ.GetInterAtomicDistance(AtomK)*AngToBohr;

	double cosPhiI =  (Rij*Rij + Rik*Rik - Rjk*Rjk)/(2.0*Rij*Rik) ;
	double cosPhiJ =  (Rij*Rij + Rjk*Rjk - Rik*Rik)/(2.0*Rij*Rjk) ;
	double cosPhiK =  (Rjk*Rjk + Rik*Rik - Rij*Rij)/(2.0*Rik*Rjk) ;
	/*
	printf("Rij = %f, Rik = %f, Rjk = %f\n",Rij,Rik,Rjk);
	printf("PhiI = %f, PhiJ = %f, PhiK = %f,  Sum = %f\n",
	       acos(cosPhiI)*RadiansToDegrees,
	       acos(cosPhiJ)*RadiansToDegrees,
	       acos(cosPhiK)*RadiansToDegrees,
	       (acos(cosPhiI)+acos(cosPhiJ)+acos(cosPhiK))*RadiansToDegrees); 
	*/

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
	//printf("Dabb(%d,%d,%d) contrib -> %f kcal/mol\n",
	//       i,j,k,damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	//       (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3))*HartreesToKcalpermole);
	E3b_disp += damping * C9ijk * (3*cosPhiI*cosPhiJ*cosPhiK + 1) /
	  (pow(Rij,3)*pow(Rik,3)*pow(Rjk,3));
	
      }
    }
  }
  return E3b_disp;

}

double Dimer::GetAIFFDampingFactor() {

  string montype[Params::Parameters().GetNumberMonomerTypes()];
  for (int i=0;i< Params::Parameters().GetNumberMonomerTypes();i++) {
    montype[i] = Params::Parameters().GetMonTypeArray()[i];
  }  

  int NoMonTypes =  Params::Parameters().GetNumberMonomerTypes();  
  double beta_damp; 

  //if (NoMonTypes == 1) {//yoni :commented out if statement so that is only one mention to state dampening factor
  //  beta_damp = Params::Parameters().GetDampingFactor(); // damping factor for this dimer
  //}
  //else {
    //cout << "MonAType  MonBType = " << MonA.GetType()<< "  " << MonB.GetType()<< "\n"; 
    string MonAType = MonA.GetType();        
    string MonBType = MonB.GetType();

    //cout << "MonAType  MonBType = " <<  MonAType << "  " <<  MonBType << "\n";
/*
    string montype[NoMonTypes];
    for (int i=0;i<NoMonTypes;i++) {
      montype[i] = Params::Parameters().GetMonTypeArray().Element(0,i);
    }
*/
//    string montype[NoMonTypes] = Params::Parameters().GetMonTypeArray();
//    StringMatrix montype(Params::Parameters().GetMonTypeArray());
//    montype = Params::Parameters().GetMonTypeArray();
//    Params::Parameters().GetMonTypeArray().Print(" Params::Parameters().GetMonTypeArray() 666");
//    montype.Print("montype 666");

    /*for (int i=0;i<NoMonTypes;i++) {
      cout << montype[i] << "\t";
      //cout << montype.Element(0,i) << "\t";
      fflush(stdout);
    }*/

    //fflush(stdout);
//    cout << montype.Element(1,i);
    int MonAType_index, MonBType_index;
    for (int i=0;i<NoMonTypes;i++) {  
      if ( MonAType.compare( montype[i] ) == 0 ) {
        MonAType_index = i;            
      }            
    }

    for (int j=0;j<NoMonTypes;j++) {
      if ( MonBType.compare( montype[j] ) == 0 ) {
        MonBType_index = j;                     
      }
    }


    beta_damp = Params::Parameters().GetDimerDampingFactor().Element( MonAType_index,  MonBType_index); // damping factor for this dimer
//  }       
    
  return beta_damp;
}

// JDH: overload for electrostatic embedding in periodic systems
void Dimer::CreateQChemJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {
  // First make the standard Q-Chem Dimer Job:
  string path = Params::Parameters().GetQMPath();
  
  string filename = path + "/d";
  char label[20]; // index ought to have fewer than 10 digits
  sprintf(label,"%d.%d",indexA,indexB);
  filename += label;
  filename += ".in";
  
  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateQChemJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  // Print comment section
  fprintf(job,"$comment\nDimer (%d,%d)\n$end\n\n",indexA,indexB);

  // Print $molecule section
  PrintQChemCartesian(job);

  // Print $rem section
  string rem = Params::Parameters().GetQChemRem();
  string basis = Params::Parameters().GetQChemBasis();

  fprintf(job,"%s\n",rem.c_str());  
  fprintf(job,"%s\n",basis.c_str());

  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    // We split things up by monomers...
    
    // First check and see what basis we want to use for monA:
    // remember
    fprintf(job,"$basis\n");
    double dist = Monomers[1].FindDistance( MonA ).Element(0);
    for (int iatom=1;iatom<=MonA.GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%s %d\n", MonA.GetAtom(iatom-1).GetSymbol().c_str(), iatom);
      if ( dist <= Params::Parameters().GetMixedBasisCutOff() ) {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
      }
      fprintf(job,"****\n");
    }
    // Now check monB
    dist = Monomers[1].FindDistance( MonB ).Element(0);
    for (int iatom=1;iatom<=MonB.GetNumberOfAtoms(); iatom++) {
      fprintf(job,"%s %d\n", MonB.GetAtom(iatom-1).GetSymbol().c_str(), iatom + MonA.GetNumberOfAtoms() );
      if ( dist <= Params::Parameters().GetMixedBasisCutOff() ) {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
      }
      fprintf(job,"****\n");
    }
    fprintf(job,"$end\n"); // end the qchem $basis section
  }
  
  // Optionally print $external_charges section
  if ( Params::Parameters().UseElectrostaticEmbedding() )  {
    fprintf(job,"$external_charges\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }

    // Now print the image charges if we are using PBCs...
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( i != indexA && i != indexB  ) {
	  if ( MonomerImages[i].GetUseInEmbedding() ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    }
    fprintf(job,"$end\n");
  }
  fclose(job);

}


void Dimer::CreateG09Job(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {

  string path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.com",indexA,indexB);
  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateG09Job : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  string gaussianHeader;
  gaussianHeader = Params::Parameters().GetGaussianHeader();
  fprintf(job,"%s\n", gaussianHeader.c_str() );
  fprintf(job,"Dimer %d,%d\n", indexA, indexB);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  if(Params::Parameters().DoCounterpoise()) { // Watit
    fprintf(job,"%d,%d %d,%d %d,%d\n", MonA.GetChargeState()+MonB.GetChargeState(), spin_state, MonA.GetChargeState(), MonA.GetSpinState(), MonB.GetChargeState(), MonB.GetSpinState());
    for(int i=0; i<MonA.GetNumberOfAtoms(); i++) {
      fprintf(job, "%s  %10.6f  %10.6f  %10.6f 1\n",MonA.GetAtom(i).GetSymbol().c_str(), MonA.GetAtom(i).GetCoordinate(0), MonA.GetAtom(i).GetCoordinate(1), MonA.GetAtom(i).GetCoordinate(2));
     }
     for(int i=0; i<MonB.GetNumberOfAtoms(); i++) {
       fprintf(job, "%s  %10.6f  %10.6f  %10.6f 2\n",MonB.GetAtom(i).GetSymbol().c_str(), MonB.GetAtom(i).GetCoordinate(0), MonB.GetAtom(i).GetCoordinate(1), MonB.GetAtom(i).GetCoordinate(2));
     }
  } 
  else {
    fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state);
    MonA.PrintMonomerCartesian(job);
    MonB.PrintMonomerCartesian(job);
  }
  
  // Print charges if requested
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if ( i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }

    // Now print the image charges if we are using PBCs...
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( MonomerImages[i].GetIndex() != indexA && MonomerImages[i].GetIndex() != indexB  ) {
	  if ( MonomerImages[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    }

  }

 if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.

    int atom_counter = 0;

    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:


  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {

	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  } // end mixed basis section
  
  
  fprintf(job,"\n");

  fclose(job);
  
}

void Dimer::CreateG09Job(Monomer Monomers[], int NMon, bool MM_job ) {

  string path;
  if(MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetMMPath();
  }
  if(!MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetQMPath();
  }
  if (MM_job && Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if(!MM_job && Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetHessianQMPath();
  }

  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.com",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateG09Job : Cannot open file '%s'\n",filename.c_str());
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
  if(Params::Parameters().DoCounterpoise()) {
    gaussianHeader += " counterpoise=2\n";
  }

  fprintf(job,"%s\n", gaussianHeader.c_str() );
  fprintf(job,"Dimer %d,%d\n", indexA, indexB);
  fprintf(job,"\n");


  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  if(Params::Parameters().DoCounterpoise()) { // Watit
    fprintf(job,"%d,%d %d,%d %d,%d\n", MonA.GetChargeState()+MonB.GetChargeState(), spin_state, MonA.GetChargeState(), MonA.GetSpinState(), MonB.GetChargeState(), MonB.GetSpinState());
    for(int i=0; i<MonA.GetNumberOfAtoms(); i++) {
      fprintf(job, "%s  %10.6f  %10.6f  %10.6f 1\n",MonA.GetAtom(i).GetSymbol().c_str(), MonA.GetAtom(i).GetCoordinate(0), MonA.GetAtom(i).GetCoordinate(1), MonA.GetAtom(i).GetCoordinate(2));
    }
    for(int i=0; i<MonB.GetNumberOfAtoms(); i++) {
      fprintf(job, "%s  %10.6f  %10.6f  %10.6f 2\n",MonB.GetAtom(i).GetSymbol().c_str(), MonB.GetAtom(i).GetCoordinate(0), MonB.GetAtom(i).GetCoordinate(1), MonB.GetAtom(i).GetCoordinate(2));
    }
  } 
  else {
    // Print charge/spin and cartesian crds
    fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state);
    MonA.PrintMonomerCartesian(job);
    MonB.PrintMonomerCartesian(job);
  }

  // Print charges if needed
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }
  }
  
  // Optionally apply mixed basis definition

  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.

    int atom_counter = 0;

    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:

  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {

	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	printf("ERROR Dimer::CreateG09Job this shouldn't happen..\n");
	exit(1);
	//region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  } // end mixed basis section
  

  fprintf(job,"\n");

  fclose(job);
  
}

void Dimer::CreateG09Job(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ) {

  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.com",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateG09Job : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  string gaussianHeader;
  gaussianHeader = Params::Parameters().GetGaussianHeader();
  fprintf(job,"%s\n", gaussianHeader.c_str() );
  fprintf(job,"Dimer %d,%d\n", indexA, indexB);
  fprintf(job,"\n");

  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  
  // Print charges if requested
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if ( i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }

    // Now print the image charges 
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( MonomerImages[i].GetIndex() != indexA && MonomerImages[i].GetIndex() != indexB  ) {
	  if ( MonomerImages[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    }

    // Now print the charges to mimic the ewald potential
    for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
      if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	fprintf(job, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      } else {
	fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      }
    }


  }

 if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }

    fprintf(job,"\n"); // Blank line to separate from crds.

    int atom_counter = 0;

    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:


  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) {
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {

	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  } // end mixed basis section
  
  
  fprintf(job,"\n");

  fclose(job);
  
}

string Dimer::RunG09Job(bool MM_job) {
  string path;

  if(MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetMMPath();
  }
  if(!MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetQMPath();
  }
  if(MM_job && Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if(!MM_job && Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20];
  sprintf(label,"%d.%d.com", indexA, indexB);
  filename += label;

  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile = "d";
  local_infile += label;
  cmd += "g09 " + local_infile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

string Dimer::RunG09ChelpGJob() {
  string path;
  path = Params::Parameters().GetMMPath();
  string filename = path + "/d";
  char label[20];
  sprintf(label,"%d.%d.chelpG.com", indexA, indexB);
  filename += label;

  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile = "d";
  local_infile += label;
  cmd += "g09 " + local_infile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

string Dimer::RunG09HirshfeldJob() {
  string path;
  path = Params::Parameters().GetMMPath();
  string filename = path + "/d";
  char label[20];
  sprintf(label,"%d.%d.Hirshfeld.com", indexA, indexB);
  filename += label;

  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile = "d";
  local_infile += label;
  cmd += "g09 " + local_infile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

string Dimer::RunDaltonJob() {
  string path;
  path = Params::Parameters().GetQMPath();
  string filename = path + "/d";
  char label[20];
  sprintf(label,"%d.%d", indexA, indexB);
  filename += label;

  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile = "d";
  local_infile += label;
  string local_outfile = "d";
  local_outfile += label;
  local_outfile += ".out";


  cmd += "dalton -mb 2000 -noarch -nobackup -o " + local_outfile + " " + local_infile;
  cmd += "; ";

  // adapt this for parallel runs: dalton -N $NP -o mN_.out mN_.dal
  
  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();

  return cmd;

}


string Dimer::RunOrcaJob() {
  string path;
  path = Params::Parameters().GetQMPath();
  //string filename = path + "/d";
  
  char label[20];
  sprintf(label,"d%d.%d.inp", indexA, indexB);
  string inputfile = label;

  string outfilename = path + "/d";
  sprintf(label,"d%d.%d.out", indexA, indexB);
  string outputfile = label;


  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  cmd += "orca " + inputfile + " > " + outputfile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

string Dimer::RunPSI4Job(bool MM_job) {
  string path;

  if(MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetMMPath();
  }
  if(!MM_job && !Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetQMPath();
  }
  if(MM_job && Params::Parameters().DoFreq()) {
    path = Params::Parameters().GetHessianMMPath();
  }
  else if(!MM_job && Params::Parameters().DoFreq())  {
    path = Params::Parameters().GetHessianQMPath();
  }

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) {
    path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  }

  string filename = path + "/d";
  char label[20];
  sprintf(label,"%d.%d.in", indexA, indexB);
  filename += label;

  // First change to local directory
  string cmd = "cd " + path;
  cmd += "; ";

  // Second command, run the job
  string local_infile = "d";
  local_infile += label;
  cmd += "psi4 " + local_infile;
  cmd += "; ";

  // Third command, switch back to base directory
  cmd += "cd " + Params::Parameters().GetBasePath();
  return cmd;
}

void Dimer::ReadMolProNMRdata() {
  printf("ERROR Dimer::ReadMolProNMRdata() has not been debugged\n");
  exit(1);




  Matrix tmpTensor(3,3);
  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();

  string path;
  path = Params::Parameters().GetQMPath();
  string filename;
  char label[120];
  char local_filename[20];
  sprintf(local_filename,"d%d.%d.molpro.out", indexA, indexB);
  filename = path + local_filename;

  string cmd;
  
  // Frist we check to see if it is an MP2 job so we know which tensors to read:
  cmd = "cd " + path;
  cmd += "; ";

  if ( 4 == 5 ) {
 
    // MP2 Job
    // Run Script to create tmp.txt file
    cmd = "cd " + path;
    cmd += "; ";
    sprintf(label,"awk -v n=2 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
    sprintf(label," | awk '{print $2, $3, $4}' > file1.out");
    cmd += label;
    cmd += "; ";
    sprintf(label,"awk -v n=3 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
    sprintf(label," | awk '{print $2, $3, $4}' > file2.out");
    cmd += label;
    cmd += "; ";

    sprintf(label,"awk -v n=4 '/LMP2 shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
    sprintf(label," | awk '{print $2, $3, $4}' > file3.out");
    cmd += label;
    cmd += "; ";
  } else {

    // HF job
    cmd = "cd " + path;
    cmd += "; ";
    sprintf(label,"awk -v n=2 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
    sprintf(label," | awk '{print $2, $3, $4}' > file1.out");
    cmd += label;
    cmd += "; ";
    sprintf(label,"awk -v n=3 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
    sprintf(label," | awk '{print $2, $3, $4}' > file2.out");
    cmd += label;
    cmd += "; ";

    sprintf(label,"awk -v n=4 '/HF shielding tensor \\(ppm\\)/{queue[NR+n]} NR in queue' ");
    cmd += label;
    cmd += "./";
    cmd += local_filename;
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
   for (int iatom = 0; iatom < (NatomsA + NatomsB); iatom++ ) {
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
     


     if (iatom < NatomsA) {
      MonA.GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
    } else {
      MonB.GetAtom(iatom-NatomsA).SetTwoBody3x3Tensor(tmpTensor);
    }

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

void Dimer::ReadDaltonNMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();
  string out_filename;
  char label[20];
  sprintf(label,"%d.%d.out", indexA, indexB );

  out_filename = path + "/d" + label;

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadDaltonNMRdata : Cannot open file '%s'\n",
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

	printf("DEBUG: tmpTensor trace = %f\n", abs(tmpTensor.Trace()) );
	if ( abs(tmpTensor.Trace()) < .001 ) {
	  printf("ERROR: Dimer::ReadDaltonNMRdata: zero tensor detected for Dimer (%d,%d), atom %d \n", indexA, indexB ,iatom++ );
	  exit(1);
	}
	
	if ( iatom <  MonA.GetNumberOfAtoms() ) {
	  MonA.GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	} else {
	  MonB.GetAtom(iatom - MonA.GetNumberOfAtoms() ).SetTwoBody3x3Tensor(tmpTensor);
	}

	iatom++;
	fflush(stdout);

	//tmpTensor.Print("DEBUG: Monomer::ReadDalton");

      }
      
    }
  }

  if ( abs(tmpTensor.Trace()) < .001 ) {
    printf("ERROR: Dimer::ReadDaltonNMRdata: zero tensor detected for Dimer (%d,%d), atom %d \n", indexA, indexB ,iatom++ );
    exit(1);
  }

  infile.close();
}

void Dimer::ReadG09NMRdata() {
  string path;
  path = Params::Parameters().GetQMPath();
  string out_filename;
  char label[20];
  sprintf(label,"%d.%d.log", indexA, indexB );

  out_filename = path + "/d" + label;

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadG09NMRdata : Cannot open file '%s'\n",
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
    printf("Dimer (%d,%d) failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndexA(), GetIndexB() );
    exit(1);
  }

  if (  convergence == false) {
    printf("Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", GetIndexA(), GetIndexB() );
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
	for (int i=1;i<=(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms() );i++) {
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
	    printf("WARNING: Dimer::ReadG09NMRdata: zero tensor detected for dimer (%d,%d), atom %d \n", indexA, indexB ,iatom++ );
	    exit(1);
	  }
	  
	  //tmpTensor.Print("Printing the two-body tensor"); OKAY

	  if ( iatom < MonA.GetNumberOfAtoms() ) {
	    MonA.GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	  } else {
	    MonB.GetAtom(iatom - MonA.GetNumberOfAtoms() ).SetTwoBody3x3Tensor(tmpTensor);
	  }

	  

	  iatom++;
	  fflush(stdout);
	  //printf("DEBUG: Dimer::ReadG09NMRdata() Shielding Tensor for Monomer: %d, atom %d \n", indexA, iatom + 1 );
	  //tmpTensor.Print("DEBUG: Dimer::ReadG09NMRdata()"); fflush(stdout); 
	  //printf("DEBUG: iso %f\n", tmpTensor.Trace() / 3);
	}

      }
      
    }

  }

  if ( abs(tmpTensor.Trace()) < .00001 ) {
    printf("ERROR: Dimer::ReadG09NMRdata: zero tensor detected for dime (%d,%d), atom %d \n", indexA, indexB ,iatom++ );
    exit(1);
  }
  	
  infile.close();

}

void Dimer::ReadG09EFGData() {
  string path;
  path = Params::Parameters().GetQMPath();
  string out_filename;
  char label[20];
  sprintf(label,"%d.%d.log", indexA, indexB );

  out_filename = path + "/d" + label;

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadG09NMRdata : Cannot open file '%s'\n",
	   out_filename.c_str() );
    exit(1);
  }

  Matrix diagonal(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms(),3);
  Matrix offdiagonal(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms(),3);
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
	  for ( int i=0;i<(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());i++) {
	    istringstream iss(line);
	    iss >> trash;
	    iss >> trash;
	    iss >> diagonal(i,0);
	    iss >> diagonal(i,1);
	    iss >> diagonal(i,2);
	    getline(infile,line);
	  }
	}

      } else {
	string match = line.substr(21,50);
      
	if ( match == "XX            YY            ZZ") {
	  
	  getline(infile,line); // Throw away the first line
	  getline(infile,line);
	  string trash;
	  for ( int i=0;i<(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());i++) {
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
	  
	  for ( int i=0;i<(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());i++) {
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
	  // All ready read in the diagonal terms:
	  for ( int i=0;i<(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());i++) {
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
	  
	  for ( int i=0;i<(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms());i++) {
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
    } // end while 
    infile.close();
  }

  // Now we reconstruct the tensors
  for (int i=0; i < (MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms()); i++) {
    Matrix tmp(3,3);
    tmp(0,0) = diagonal(i,0);
    tmp(1,1) = diagonal(i,1);
    tmp(2,2) = diagonal(i,2);
    
    tmp(0,1) = offdiagonal(i,0);
    tmp(0,2) = offdiagonal(i,1);
    tmp(1,2) = offdiagonal(i,2); 
    
    tmp(1,0) = tmp(0,1);
    tmp(2,0) = tmp(0,2);
    tmp(2,1) = tmp(1,2);
    if ( i< MonA.GetNumberOfAtoms() ) {
      MonA.GetAtom(i).SetTwoBody3x3Tensor(tmp);
    } else {
      MonB.GetAtom(i - MonA.GetNumberOfAtoms() ).SetTwoBody3x3Tensor(tmp);
    }

  }

}

void Dimer::ReadQChemNMRdata() {
  // main routine for reading in NMR data, Grabs the 'total shielding tensor' from qchem
  // largely based on prof. Beran's code for reading from file

  // Temporary storage for 3x3 NMR Shielding tensor
  Matrix tmpTensor(3,3);

  int NatomsA = MonA.GetNumberOfAtoms();
  int NatomsB = MonB.GetNumberOfAtoms();
  
  // Set up the filename, with the full path.  File is e.g. 'm1.force'
  string path = Params::Parameters().GetQMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.qchem",indexA,indexB);

  string filename = path + label + ".out";
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadQChemElectrostaticEmbeddedNMRdata : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  // Read in the data - search for the "NMR-SHIELDING TENSORS" string
  // to know where it starts.
  // Note: this is a really sloppy routine.  Assumes no format changes,
  // rigid fields/spacing in some cases, etc.  Ought to be smarter!
  int iatom=0;
  string line;
  while ( !infile.eof() ) {
    getline(infile,line);
    // Search for NMR section
    if (line.length() > 40) {
      string match = line.substr(20,21);
      if ( match=="NMR-SHIELDING TENSORS") {
	for (int i=0;i<31;i++) {// feed forward to first atom
	  getline(infile,line);
	}
	iatom++;
	tmpTensor.Set();
	for (int i=0;i<3;i++) {// read the actual tensor
	  getline(infile,line);
	  istringstream iss(line);
	  iss >> tmpTensor(i,0);
	  iss >> tmpTensor(i,1);
	  iss >> tmpTensor(i,2);
	}
	
       
	fflush(stdout);
	MonA.GetAtom(iatom-1).SetTwoBody3x3Tensor(tmpTensor);
	
	while (iatom < Natoms) {
	  for (int i=0;i<29;i++) // feed forward to next atom's tensor
	    getline(infile,line);
	  iatom++;
	  tmpTensor.Set();
	  for (int i=0;i<3;i++) {// read the actual tensor
	    getline(infile,line);
	    istringstream iss(line);
	    iss >> tmpTensor(i,0);
	    iss >> tmpTensor(i,1);
	    iss >> tmpTensor(i,2);
	  }
	  //printf("atom %d\n",iatom); //JDH: debuging
	  //tmpTensor.Print("Shielding Tensor:"); //JDH: debuging, this looks good...
	  fflush(stdout); // flush the output stream
	  if (iatom <= NatomsA) {
	    MonA.GetAtom(iatom-1).SetTwoBody3x3Tensor(tmpTensor);
	  }
	  else {
	    MonB.GetAtom(iatom-NatomsA-1).SetTwoBody3x3Tensor(tmpTensor);
	  }
	}
      }
    }
    if (line.length() > 34) {
      string match = line.substr(2,33);
      if ( match=="Summary of detailed contributions") {
	//printf("Found end!  Read NMR tensors for %d atoms\n",iatom);
	break;
      }
    }
  }

  infile.close(); 
}


void Dimer::CreateG09HirshfeldJob() {

  //printf("DEBUG: creating hirshfeld job for dimer(%d,%d)...\n",indexA,indexB); fflush(stdout);
  
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.Hirshfeld",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09HirshfeldJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for Hirshfeld run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

  // Print External Charges? Maybe later...

  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<MonA.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      region = MonB.GetAtom(iatom).GetMixedBasisRegion();
      

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("custom vdw radii and hershfeld?\n");
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




void Dimer::CreateG09HirshfeldJob(int NMon, Monomer Monomers[] ){
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.Hirshfeld",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09HirshfeldJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for Hirshfeld run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

  // Print External Charges?
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }
  }
  

  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	printf("You shouldn't be here...\n");
	exit(1);
	// Uncomment this when we are dealing with Periodic systems
	//region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("custom vdw radii and hershfeld?\n");
    exit(1);
  } else if (custom_vdw_radii == false ) {
    cmd = "sed -i s/charge/'charge pop=Hirshfeld'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());


}

void Dimer::CreateG09HirshfeldJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ){
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.Hirshfeld",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for Hirshfeld run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

  // Print External Charges? Maybe later...

  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	// Uncomment this when we are dealing with Periodic systems
	region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  string cmd;
  if ( custom_vdw_radii ) {
    //cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
    printf("vdw with hershfeld\n");
    exit(1);
  } else if (custom_vdw_radii == false ) {
    cmd = "sed -i s/charge/'charge pop=Hirshfeld'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());


}


void Dimer::CreateG09ChelpGJob() {

  printf("DEBUG: Mixed basis set won't work properly with  Dimer::CreateG09ChelpGJob() \n");
  exit(1);

  
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.chelpG",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateG09ChelpGJob : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for ChelpG run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);



  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<MonA.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( MonA.GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}



	// else {
	//   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() ); // Always assigne level 1 basis set to monomer in unit cell for two-body charge calculations
	// }


	

      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      region = MonB.GetAtom(iatom).GetMixedBasisRegion();
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  }



  
  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);

  


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



void Dimer::CreateG09ChelpGJob(int NMon, Monomer Monomers[] ) {
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.chelpG",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for ChelpG run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

  // Print External Charges
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }
  }
  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() ); 
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }

	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}


	// else {
	//   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );// Always assigne level 1 basis set to monomer in unit cell for two-body charge calculations
	// }





      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	printf("You shouldn't be here...\n");
	exit(1);
	// Uncomment this when we are dealing with Periodic systems
	//region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }

      fprintf(job,"****\n");
    }
    
    
    
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else if (custom_vdw_radii == false ) {
    cmd = "sed -i s/charge/'charge pop=ChelpG'/ ";
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());


}


void Dimer::CreateG09ChelpGJob(int NMon, Monomer Monomers[], int NMon_images, Monomer MonomerImages[] ){

  //printf("DEBUG XXX Building Two-Body Charge Dimer(%d,%d) \n", indexA, indexB);

  
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.chelpG",indexA,indexB);

  string filename = path + label + ".com";
  
  // Open the input file for writing
  FILE *job;
  if (( job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Monomer::CreateG09Job : cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"%s\n", Params::Parameters().GetGaussianHeader().c_str() );
  fprintf(job,"Dimer (%d,%d) Gaussian input for ChelpG run\n", indexA, indexB);
  fprintf(job,"\n");

 
  // Print charge/spin and cartesian crds
  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateG09Job, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"%d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

  // Print External Charges? Maybe later...

  
  if ( Params::Parameters().CustomBasis() == 2 ) {
    // Use Custom basis for all atoms, not compatible with mixxed basis
    if ( Params::Parameters().GetMixedBasisCutOff() > 0 ) { 
      printf("ERROR, custom basis on all atoms requested, not compatible with mixed basis\n");
      exit(1);
    }
    
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms(); iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
    }

    for (int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetNBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetOBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetSBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetClBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetIBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetSnBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetPBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetKBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetNaBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
      } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );	
      } else {
	fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
      }
      fprintf(job,"****\n");
      
    }

  } // end full custom basis section:
   
   
  // Optionally apply mixed basis definition
  if ( Params::Parameters().GetMixedBasisCutOff() > 0 ){
    fprintf(job,"\n"); // Blank line to separate from crds.
    
    int atom_counter = 0;
    
    // loop over all the atoms and print the appropriate basis...
    for ( int iatom=0;iatom<Monomers[indexA].GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	
	if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
	  
	  if ( MonA.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	  
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}

	// else {
	//   fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() ); // Always assigne level 1 basis set to monomer in unit cell for two-body charge calculations
	// }

      }

      fprintf(job,"****\n");
    }
    
    for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
      atom_counter++;
      fprintf(job,"%d 0\n", atom_counter);

      int region;
      // Get region...
      if ( MonB.GetIndex() > NMon) {
	// Uncomment this when we are dealing with Periodic systems
	region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
      } else {
	region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
      }

      //printf("DEBUG XXX: D(%d,%d) region = %d \n", indexA, indexB, region);
      
      // Now print the basis based on the region...
      if ( Params::Parameters().CustomBasis() == 0 ) { 
	if ( region == 1 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      } else if ( Params::Parameters().CustomBasis() == 1) {
	if ( region == 1 ) {
	  if ( MonB.GetAtom(iatom).GetSymbol() == "H" && Params::Parameters().GetHBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetHBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "C" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetCBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "N" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetOBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetClBasis().c_str() );  
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetIBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetSnBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetPBasis().c_str() ); 
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetKBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" && Params::Parameters().GetCBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetNaBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "Br" && Params::Parameters().GetBrBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetBrBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" && Params::Parameters().GetFBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetFBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" && Params::Parameters().GetWBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetWBasis().c_str() );
	  } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" && Params::Parameters().GetVBasis() != "" ) {
	    fprintf(job,"%s\n", Params::Parameters().GetVBasis().c_str() );  
	  } else {
	    fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel1().c_str() );
	  }
	} else if ( region == 2 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel2().c_str() );
	} else if ( region == 3 ) {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	} else {
	  fprintf(job,"%s\n", Params::Parameters().GetMixedBasisLevel3().c_str() );
	}
      }
 
      fprintf(job,"****\n");
    }
    
    
    
  }


  bool custom_vdw_radii = false;
  // The ChelpG implementation in G09 doesn't contain VDW radii
  // for all the nuclei we care about, so we have to loop through
  // and determin where we need a special VDW section for the job:
  // List of Nuclei we care about yet G09 does not have tabulated
  // VDW radii for:  Br, I, Sn, K  
  // Loop through all of the atoms and if one of these atoms appears, 
  // add a VDW radii definition to the G09 input file:

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms(); iatom++) {
    if ( MonA.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonA.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }
  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms(); iatom++) {
    if ( MonB.GetAtom(iatom).GetSymbol() == "Br"  ) {
      fprintf(job,"Br %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I"  ) {
      fprintf(job,"I %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Sn" ) {
      fprintf(job,"Sn %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      fprintf(job,"K %f\n", MonB.GetAtom(iatom).GetVanDerWaalsRadius() );
      custom_vdw_radii = true;
    }
  }

  fprintf(job,"\n");
  fclose(job);


  string cmd;
  if ( custom_vdw_radii ) {
    cmd = "sed -i s/charge/'pop=(ChelpG,ReadRadii,Mk)'/ ";
  } else if (custom_vdw_radii == false ) {
    //cmd = "sed -i s/charge/'charge pop=ChelpG'/ ";
    cmd = "sed -i s/charge/'pop=ChelpG'/ "; // Will need to add additional function to include charge emebedding for two-body SCE calculations
  }
  cmd += filename.c_str();
  system(cmd.c_str());
  cmd = "sed -i s/nmr// ";
  cmd += filename.c_str();
  system(cmd.c_str());


}

void Dimer::ReadHirshfeldCharges() {
  
  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.Hirshfeld",indexA,indexB);

  string filename = path + label + ".log";
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadHirshfledCharges : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  
  int iatom = 0;
  string line;
  string search_match = "Hirshfeld charges, spin densities"; // 33 characters
  bool normal_term = false;
  bool convergence = true;
  int convergence_issues = 0;
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
	convergence_issues++;
      }

      
      
    }
  }
  
  if ( normal_term == false) {
    printf("Hirshfeld calculation for Dimer (%d,%d) failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );
    exit(1);
  }

  if ( convergence_issues > 1  ) {
    printf("Hirshfeld calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );
    exit(1);
  } else if ( convergence_issues == 1) {
    printf("Warning Hirshfeld calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );

  }

  
  // if (  convergence == false) {
  //   printf("Hirshfeld calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );
  //   exit(1);
  // }
  
  infile.close();
  

  infile.open( filename.c_str() );
  double *chelpg_charges = new double[4]; // apparently this need so
  chelpg_charges[0]=0.0;
  chelpg_charges[1]=0;
  chelpg_charges[2]=0;
  chelpg_charges[3]=0;
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 25 ) {
      // Search for Hirshfelf Charges Section
      string match = line.substr(1,33);
      if ( match  == search_match ) {
	getline(infile,line);// Throw away first line
	getline(infile,line);

	//Extract Hirshfeld Charges and add them to the array
	string trash;
	fflush(stdout);
	for (int i=0;i<MonA.GetNumberOfAtoms();i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);

	  printf("Hirshfeld charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] ); fflush(stdout);
	  
	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  MonA.GetAtom(i).SetMultipoleMoments(Charge);
	  
	  //printf("Length of multipole: %d\n",  MonA.GetAtom(i).GetMultipoleMoments().GetLength() );

	  printf("Hirshfeld charge extracted: Atom %d: Charge = %f\n\n", i+1, MonA.GetAtom(i).GetMultipoleMoments().GetMoments().Element(0)   );

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}

	// Now read in the MonB charges just in case:
	for (int i=0;i<MonB.GetNumberOfAtoms();i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);

	  printf("Hirshfeld charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] ); fflush(stdout);
	  
	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  MonB.GetAtom(i).SetMultipoleMoments(Charge);
	  
	  //printf("Length of multipole: %d\n",  MonA.GetAtom(i).GetMultipoleMoments().GetLength() );

	  printf("Hirshfeld charge extracted: Atom %d: Charge = %f\n\n", i+1, MonB.GetAtom(i).GetMultipoleMoments().GetMoments().Element(0)   );

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}

	
      }
    }
    
    
  } // end while loop

  delete [] chelpg_charges;

  infile.close();

}

void Dimer::ReadChelpGCharges() {

  string path = Params::Parameters().GetMMPath();
  char label[20]; // index ought to have fewer than 20 digits
  sprintf(label,"/d%d.%d.chelpG",indexA,indexB);

  string filename = path + label + ".log";
  
  // Open the force file
  ifstream infile;
  infile.open( filename.c_str() );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadChelpGCharges : Cannot open file '%s'\n",
	   filename.c_str());
    exit(1);
  }

  
  int iatom = 0;
  string line;
  string search_match = "Charges from ESP fit, RMS=";
  bool normal_term = false;
  bool convergence = true;
  int convergence_issues = 0;
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
	convergence_issues++;
	//printf("DEBUG: here is the string! %s\n", line.c_str() ); 
      }
    }
  }

  //  printf("DEBUG nuber of convergence issues: %d\n", convergence_issues );
  
  if ( normal_term == false) {
    printf("ChelpG calculation for Dimer (%d,%d) failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );
    exit(1);
  }

  if ( convergence_issues > 1  ) {
    //    printf("convergence_issues = %d\n", convergence_issues);
    printf("\nChelpG calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n\n", MonA.GetIndex(), MonB.GetIndex() );
    exit(1);
  } else if ( convergence_issues == 1) {
    printf("\nWarning ChelpG calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n\n", MonA.GetIndex(), MonB.GetIndex() );

  }  


  
  // if (  convergence == false) {
  //   printf("ChelpG calculation Dimer (%d,%d) is experiencing SCF convergence issues:\n this can often be addressed by reading in a converged\n density using a smaller basis set.\n Once fixed, re-run hmbi using 'analyze_only = true'\n", MonA.GetIndex(), MonB.GetIndex() );
  //   exit(1);
  // }
  
  infile.close();
  

  infile.open( filename.c_str() );
  double *chelpg_charges = new double[4]; // apparently this need so
  chelpg_charges[0]=0.0;
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
	for (int i=0;i<MonA.GetNumberOfAtoms();i++) { 
	  //printf("DEBUG: Reading in charges for atom %d\n", i); fflush(stdout);
	  
	  //chelpg_chaclrges[0] = MonA.GetAtom(i).GetMultipoleMoments().GetMoments().Element(0); 

	  //printf("DEBUG: read charges for saved atom %d\n", i); fflush(stdout);
	  

	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);

	  //printf("ChelpG charge extracted: Atom %d: Charge = %f\n", i+1, chelpg_charges[0] ); fflush(stdout);

	  if ( !(abs(chelpg_charges[0]) > 0.0 && abs(chelpg_charges[0]) < 100) ) {
	    printf("ERROR reading in ChelpG Charges for Dimer (%d,%d)\n",indexA,indexB);
	    exit(1);
	  }
	  
	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  MonA.GetAtom(i).SetMultipoleMoments(Charge);
	  
	  //printf("Length of multipole: %d\n",  MonA.GetAtom(i).GetMultipoleMoments().GetLength() );

	  //printf("ChelpG charge extracted: Atom %d: Charge = %f\n\n", i+1, MonA.GetAtom(i).GetMultipoleMoments().GetMoments().Element(0)   );

	  //printf("DEBUG: monomer %s: ensure CHelpG charge was assigned to multipole momoent properly: charge = %f\n", mN_.c_str() , GetAtom(i).GetMultipoleMoments().GetMoments().Element(0) ); fflush(stdout);

	}

	// Read in MonB charges...
	for (int i=0;i<MonB.GetNumberOfAtoms();i++) { 
	  istringstream iss(line);
	  iss >> trash;
	  iss >> trash;
	  iss >> chelpg_charges[0];
	  getline(infile,line);

	  // Initialize the Multipole object with rank=0 and moments=CHelpG charge
	  Multipole Charge(1, chelpg_charges);
	  MonB.GetAtom(i).SetMultipoleMoments(Charge);
	  
	}
	
      }
    }
    
    
  } // end while loop

  delete [] chelpg_charges;

  infile.close();

}

// Watit - Raman intensities
void Dimer::ReadHess() {

        int NA = MonA.GetNumberOfAtoms();
        int NB = MonB.GetNumberOfAtoms();
        Matrix HessQM(3*(NA+NB),3*(NA+NB));
        Matrix HessMM(3*(NA+NB),3*(NA+NB));
        string path;
        string out_filename;
	char label[20];
        if(Params::Parameters().GetQMType()==3) { // G09
                path = Params::Parameters().GetHessianQMPath();
		sprintf(label,"%d.%d.log",indexA,indexB);
                out_filename = path + "/d" + label;
        }	

	string line;
        ifstream infile;
        int DOFA = 3*NA; // Define degrees of freedom
        int DOFB = 3*NB;

	infile.open(out_filename.c_str());
        if(!infile.is_open()) {
                printf("Dimer::ReadHess QM : Cannot open file '%s'\n",out_filename.c_str());
                exit(1);
        }

	if(Params::Parameters().GetQMType()==3) {
		string line;
		while (!infile.eof()) {
                        getline(infile,line);
                        if(Params::Parameters().DoCounterpoise()) {
                                if(line.find("Counterpoise corrected") != string::npos) {
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
                        } else {
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
	}
	infile.close();
	SetHessQM(HessQM);
	//HessQM.PrintHessian("2-body QM Hessian");

	//MM Hessian
	if(Params::Parameters().GetMMType()==1) { // Tinker
                path = Params::Parameters().GetHessianMMPath();
                sprintf(label,"%d.%d.freq",indexA,indexB);
                out_filename = path + "/d" + label;
        }

        infile.open(out_filename.c_str());
        if(!infile.is_open()) {
                printf("Dimer::ReadHess MM : Cannot open file '%s'\n",out_filename.c_str());
                exit(1);
        }
	
	if(Params::Parameters().GetMMType()==1) {
                string line;
                while(!infile.eof()) {
                        getline(infile,line);
			if(line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):") {
                                for(int i=0;i<3*(NA+NB);i++) {
                                        for(int j=0;j<3*(NA+NB);j++) {
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
        //HessMM.PrintHessian("2-body MM Hessian");

}

void Dimer::ReadIntensity() {
	// Set path to hessian_files/qm
	string path;
	path = Params::Parameters().GetHessianQMPath();
	string out_filename;
	char label[20];
	if(Params::Parameters().GetQMType()==1) { // QChem
		sprintf(label,"%d.%d.out",indexA,indexB);
	} else if(Params::Parameters().GetQMType()==3) { // G09
		sprintf(label,"%d.%d.log",indexA,indexB);
	}
	out_filename = path + "/d" + label;

	// Read PolD from output file
	string line;
	ifstream infile;
	int DOFA = 3*MonA.GetNumberOfAtoms(); // Define degrees of freedom
	int DOFB = 3*MonB.GetNumberOfAtoms();

	infile.open(out_filename.c_str());
	if(!infile.is_open()) {
		printf("Dimer::ReadIntensity : Cannot open file '%s'\n",out_filename.c_str());
		exit(1);
	}
        Matrix DipDA(DOFA,3);
        Matrix DipDB(DOFB,3);
	Matrix PolDA(DOFA,6);
	Matrix PolDB(DOFB,6);
	int check_DipD=0;
	int check_PolD=0;
	if(Params::Parameters().GetQMType()==1) { // QChem
		while(getline (infile,line)) {
			if(line.find("DipDeriv") != string::npos) {
				check_DipD++;
                        	getline(infile,line);
                        	for(int i=0;i<DOFA;i++) {
                                	getline(infile,line);
                                	istringstream iss(line);
                                	int trash;
                                	iss >> trash;
                                	for (int j=0;j<3;j++) {
                                        	iss >> DipDA(i,j);
                                	}
                        	}
                        	for(int i=0;i<DOFB;i++) {
                                	getline(infile,line);
                                	istringstream iss(line);
                                	int trash;
                                	iss >> trash;
                                	for (int j=0;j<3;j++) {
                                     		iss >> DipDB(i,j);
                                	}
                        	}
                	}
			if(Params::Parameters().DoRaman()) {
				if(line.find("PolDeriv") != string::npos) {
					check_PolD++;
					getline(infile,line);
					for(int i=0;i<DOFA;i++) {
						getline(infile,line);
						istringstream iss(line);
						int trash;
						iss >> trash;
						for (int j=0;j<6;j++) {
							iss >> PolDA(i,j);
						}
					}
					for(int i=0;i<DOFB;i++) {
						getline(infile,line);
						istringstream iss(line);
						int trash;
						iss >> trash;
						for (int j=0;j<6;j++) {
							iss >> PolDB(i,j);
						}
					}
				}
			}
		}
	} else if(Params::Parameters().GetQMType()==3) { // G09
                while(getline(infile,line)) {
                        if(line.substr(0,17)==" Hessian entering") {
				while(getline(infile,line)) {
                                	if(line.substr(0,15)==" DipoleDeriv   ") {
                                        	check_DipD++;
        	                                for(int i=0;i<DOFA;i++) {
                                	                string tmp;
        	                                        for (int j=0;j<3;j++) {
                	                                        stringstream ss;
								tmp = line.substr(j*15+16,j*15+31);
                                	                        replace(tmp.begin(),tmp.end(),'D','E');
                                        	                ss << tmp;
                                                	        ss >> DipDA(i,j);
								//DipDA(i,j)=-DipDA(i,j);
	                                                }
							getline(infile,line);
        	                                }
                	                        for(int i=0;i<DOFB;i++) {
                                        	        string tmp;
                                                	for (int j=0;j<3;j++) {
	                                                        stringstream ss;
								tmp = line.substr(j*15+16,j*15+31);
                	                                        replace(tmp.begin(),tmp.end(),'D','E');
                        	                                ss << tmp;
                                	                        ss >> DipDB(i,j);
								//DipDB(i,j)=-DipDB(i,j);
                                        	        }
							getline(infile,line);
	                                        }
        	                                if(Params::Parameters().DoRaman()) {
                		                                if(line.substr(0,15)==" PolarDeriv    ") {
                        		                        check_PolD++;
                                        		        for(int i=0;i<DOFA;i++) {
	                                    		                string tmp;
                                        	          	      	for (int j=0;j<3;j++) {
                                                	                	stringstream ss;
										tmp = line.substr(j*15+16,j*15+31);
                                                                		replace(tmp.begin(),tmp.end(),'D','E');
                                                                		ss << tmp;
	                                                                	ss >> PolDA(i,j);
										//PolDA(i,j)=-PolDA(i,j);
        	                                                	}
                	                                        	getline(infile,line);
                        	                                	for(int j=3;j<6;j++) {
                                	                                	stringstream ss;
										tmp = line.substr((j-3)*15+16,(j-3)*15+31);
                                                	                	replace(tmp.begin(),tmp.end(),'D','E');
                                                        	        	ss << tmp;
                                                                		ss >> PolDA(i,j);
										//PolDA(i,j)=-PolDA(i,j);
	                                                        	}
									getline(infile,line);	
        	                                        	}
                	                                	for(int i=0;i<DOFB;i++) {
                                        	                	string tmp;
                                                	        	for (int j=0;j<3;j++) {
                                                        	        	stringstream ss;
										tmp = line.substr(j*15+16,j*15+31);
                                                                		replace(tmp.begin(),tmp.end(),'D','E');
	                                                                	ss << tmp;
        	                                                        	ss >> PolDB(i,j);
										//PolDB(i,j)=-PolDB(i,j);
                	                                        	}
                        	                                	getline(infile,line);
                                	                        	for (int j=3;j<6;j++) {
                                        	                        	stringstream ss;
										tmp = line.substr((j-3)*15+16,(j-3)*15+31);
                                                        	        	replace(tmp.begin(),tmp.end(),'D','E');
                                                                		ss << tmp;
                                                                		ss >> PolDB(i,j);
										//PolDB(i,j)=-PolDB(i,j);
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
	//DipDA.Print("DipDA");
	//DipDB.Print("DipDB");
	//cout << "check " << check_DipD << endl;
        if(check_DipD < 1) {
                printf("Dimer::ReadIntensity : Cannot find Dipole Derivative for '%s'\n",out_filename.c_str());
                exit(1);
	}
	if(Params::Parameters().DoRaman()&&check_PolD < 1) {
		printf("Dimer::ReadIntensity : Cannot find Polarizability Derivative for '%s'\n",out_filename.c_str());
		exit(1);
	}
	infile.close();
        SetDipDA(DipDA);
        SetDipDB(DipDB);
	SetPolDA(PolDA);
	SetPolDB(PolDB);

	//printf("d%i.%i\n",indexA,indexB);
	//DipDA.Print("DipDA");
	//DipDB.Print("DipDB");
}


void Dimer::CreateDaltonJob(Monomer Monomers[], int NMon){


  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.dal",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tDimer: %d,%d\n", indexA, indexB);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Spherical Angstrom Nosymmetry\n", MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms() );

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonA.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;  
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonA.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    
    if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    } else {
      exit(1);
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    }
    MonA.GetAtom(iatom).PrintQChemCartesian(job);
  }

  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonB.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;  
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonB.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    
    if ( Monomers[indexB].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( Monomers[indexB].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( Monomers[indexB].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    }
    MonB.GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);


  // Make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    filename = path + "/d";
    sprintf(label,"%d.%d.pot",indexA,indexB);
    filename += label;
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }



    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != indexA && imon != indexB ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom<  Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }

    

    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");

    int counter = 1;

    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != indexA && imon != indexB ){ 
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
      if (  imon != indexA && imon != indexB  ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }
    

    // RANK 1 Terms: Dipole
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if (  imon != indexA && imon != indexB )  {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
      
    } // end RANK 1

    // RANK 2 Terms: Quadrupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if (  imon != indexA && imon != indexB  ) {
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
      
    } // end RANK 2


    // RANK 3 Terms: Octupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if (  imon != indexA && imon != indexB   ) {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
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
      }
      

    } //end rank 3

    fclose(job);

  } // end of charge embedding

}



void Dimer::CreateDaltonJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ){
  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.dal",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateDaltonJob, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tDimer: %d,%d\n", indexA, indexB);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Charge=%d Spin=%d Spherical Angstrom Nosymmetry\n", MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms(), MonA.GetChargeState() + MonB.GetChargeState(), spin_state );

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonA.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if (MonA.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if (MonA.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonA.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    
    if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    } else {
      printf("This should never happen, atom not mapped to basis set region... \n");
      exit(1);
    }
    
    MonA.GetAtom(iatom).PrintQChemCartesian(job);
  }

  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {

    // Get region...
    int region = 9999999;
    if ( MonB.GetIndex() > NMon ) {
      region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
    } else {
      region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
    }


    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonB.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if (MonB.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if (MonB.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonB.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    

    //string basis = Params::Parameters().GetNMRMixedBasisLevel1().c_str();
    //std::transform(basis.begin(), basis.end(), basis.begin(), ::tolower);
    //cout << "DEBUG XXX: string: " << basis  << endl; // this works...

    if ( region == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( region == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( region == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    } else {
      printf("ERROR: Dimer::CreateDaltonJob (overloaded for PBC) undefined basis set region for mixed basis calcualtion: \n");
      printf("error found in dimer (%d,%d) atom %d\n", indexA, indexB, iatom );
      printf("requested a region = %d\n", region);
      exit(1);
    }
    MonB.GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");
  // fclose(job);

  
  // filename = path + "/d";
  // sprintf(label,"%d.%d.dal",indexA,indexB);
  // filename += label;

  // if ((job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }

  // fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);


  // Make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    filename = path + "/d";
    sprintf(label,"%d.%d.pot",indexA,indexB);
    filename += label;
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }



    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != indexA && imon != indexB ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom<  Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }

    for (int imon=1; imon<= NMon_images ; imon++) {
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
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
    }
    

    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");

    int counter = 1;

    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != indexA && imon != indexB ){ 
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
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  //printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	  for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
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
      if (  imon != indexA && imon != indexB  ) {
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
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	    Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }

    // RANK 1 Terms: Dipole
    if (Params::Parameters().GetChargeEmbeddingRank() > 0 ) {
      fprintf(job, "ORDER 1\n");
      fprintf(job, "%d\n",total_atoms);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)   {
	if (  imon != indexA && imon != indexB )  {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<Monomers[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
      
      for (int imon=1; imon<= NMon_images ; imon++) {
	if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
	  if (MonomerImages[imon].GetUseInEmbedding() ) {
	    for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	      Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	      fprintf(job, "%d   %4.5f  %4.5f  %4.5f\n",counter,moments[1],moments[2],moments[3]);
	      counter++;
	    }
	  }
	}
      }
    } // end RANK 1

    // RANK 2 Terms: Quadrupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 1 ) {
      fprintf(job, "ORDER 2\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++) {
	if (  imon != indexA && imon != indexB  ) {
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
	if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
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
      }
    } // end RANK 2


    // RANK 3 Terms: Octupole
    if (Params::Parameters().GetChargeEmbeddingRank() > 2 ) {
      fprintf(job, "ORDER 3\n");
      fprintf(job, "%d\n",high);
      counter=1;
      for (int imon=1; imon<=NMon; imon++)  {
	if (  imon != indexA && imon != indexB   ) {
	  if ( Monomers[imon].GetUseInEmbedding() ) {
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
      }
      
      for (int imon=1; imon<=NMon_images; imon++)  {
	if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
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
      }
    } //end rank 3

    fclose(job);

  } // end of charge embedding

}

void Dimer::CreateDaltonJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ){
  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.dal",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }

  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateDaltonJob, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  fprintf(job,"ATOMBASIS\n");
  fprintf(job,"\tDimer: %d,%d\n", indexA, indexB);
  fprintf(job,"\tProbably don't need this line...\n");

  fprintf(job,"Atomtypes=%d Charge=%d Spin=%d Spherical Angstrom Nosymmetry\n", MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms(), MonA.GetChargeState() + MonB.GetChargeState(), spin_state );

  for ( int iatom=0;iatom<MonA.GetNumberOfAtoms();iatom++) {
    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonA.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonA.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if (MonA.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if (MonA.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonA.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonA.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    
    if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( Monomers[indexA].GetAtom(iatom).GetMixedBasisRegion() == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    } else {
      printf("This should never happen, atom not mapped to basis set region... \n");
      exit(1);
    }
    
    MonA.GetAtom(iatom).PrintQChemCartesian(job);
  }

  for ( int iatom=0;iatom<MonB.GetNumberOfAtoms();iatom++) {

    // Get region...
    int region = 9999999;
    if ( MonB.GetIndex() > NMon ) {
      region = MonomerImages[ MonB.GetIndex() - NMon].GetAtom(iatom).GetMixedBasisRegion();
    } else {
      region = Monomers[ MonB.GetIndex() ].GetAtom(iatom).GetMixedBasisRegion();
    }


    // There is probably a better way to do this, but I'm in a hurry
    double charge;
    if ( MonB.GetAtom(iatom).GetSymbol() == "H" ) {
      charge = 1.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "C" ) {
      charge = 6.0;
    }	else if ( MonB.GetAtom(iatom).GetSymbol() == "N" ) {
      charge = 7.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "O" ) {
      charge = 8.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "F" ) {
      charge = 9.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "P" ) {
      charge = 15.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "S" ) {
      charge = 16.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cl" ) {
      charge = 17.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else if (MonB.GetAtom(iatom).GetSymbol() == "Na" ) {
      charge = 11.0;
    } else if (MonB.GetAtom(iatom).GetSymbol() == "I" ) {
      charge = 53.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "V" ) {
      charge = 23.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "W" ) {
      charge = 74.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "Cs" ) {
      charge = 55.0;
    } else if ( MonB.GetAtom(iatom).GetSymbol() == "K" ) {
      charge = 19.0;
    } else {
      printf("ERROR: building dalton dimer job, unknown atom type: %s\n", MonB.GetAtom(iatom).GetSymbol().c_str()  );
      printf("Don't panic, just grep this line and add the atom you want in\n");
      exit(1);
    }
    

    //string basis = Params::Parameters().GetNMRMixedBasisLevel1().c_str();
    //std::transform(basis.begin(), basis.end(), basis.begin(), ::tolower);
    //cout << "DEBUG XXX: string: " << basis  << endl; // this works...

    if ( region == 1 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel1().c_str() );
    } else if ( region == 2 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel2().c_str() );
    } else if ( region == 3 ) {
      fprintf(job,"Charge=%f Atoms=1 Basis=%s\n", charge, Params::Parameters().GetMixedBasisLevel3().c_str() );
    } else {
      printf("ERROR: Dimer::CreateDaltonJob (overloaded for PBC) undefined basis set region for mixed basis calcualtion: \n");
      printf("error found in dimer (%d,%d) atom %d\n", indexA, indexB, iatom );
      printf("requested a region = %d\n", region);
      exit(1);
    }
    MonB.GetAtom(iatom).PrintQChemCartesian(job);
  }

  fprintf(job,"\n");
  // fclose(job);

  
  // filename = path + "/d";
  // sprintf(label,"%d.%d.dal",indexA,indexB);
  // filename += label;

  // if ((job = fopen(filename.c_str(),"w"))==NULL) {
  //   printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
  //   exit(1);
  // }

  // fprintf(job,"\n");
  fprintf(job,"%s\n", Params::Parameters().GetDaltonSection().c_str() );
  fprintf(job,"\n"); // blank line to end the file...
  fclose(job);


  // Make the potential file for electrostatic embedding
  if ( Params::Parameters().UseElectrostaticEmbedding() ) { 
    filename = path + "/d";
    sprintf(label,"%d.%d.pot",indexA,indexB);
    filename += label;
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateDaltonJob : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }



    // Count the number of atoms:
    int total_atoms=0;
    int high=0;
    //for (int imon=1; imon<=min(NMon,Params::Parameters().GetNumAsymmetricMonomers() ); imon++) {
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != indexA && imon != indexB ) {
	if ( Monomers[imon].GetUseInEmbedding() ) {
	  total_atoms += Monomers[imon].GetNumberOfAtoms();
	  for(int iatom=0; iatom<  Monomers[imon].GetNumberOfAtoms(); iatom++) {
	    int rank = Monomers[imon].GetAtom(iatom).GetMultipoleMoments().GetRank();
	    if(rank >= 2) { // if rank is > 2 dealing with higher order moments
	      high++;
	    }
	  }
	}
      }
    }

    for (int imon=1; imon<= NMon_images ; imon++) {
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
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
    }

    total_atoms += EwaldCharges.GetRows();

    fprintf(job, "! Charge embedding file Generated in HMBI\n");
    fprintf(job, "@COORDINATES\n");
    fprintf(job, "%d\n",total_atoms);
    fprintf(job, "AA\n");

    int counter = 1;

    // Print the multipole emb. crds for inside the unit cell
    for (int imon=1; imon<=NMon; imon++) {
      if ( imon != indexA && imon != indexB ){ 
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
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  //printf("DEUBG: \t \t Printing crds for image monomer %d\n", imon);
	  for(int iatom=0; iatom< MonomerImages[imon].GetNumberOfAtoms(); iatom++) {
	    fprintf(job, "%s  %4.8f  %4.8f  %4.8f\n", MonomerImages[imon].GetAtom(iatom).GetSymbol().c_str(), MonomerImages[imon].GetAtom(iatom).GetPosition(0),MonomerImages[imon].GetAtom(iatom).GetPosition(1),MonomerImages[imon].GetAtom(iatom).GetPosition(2));
	    counter++;
	  }
	}
      }
    }
    
    
    // Print the mult. emb. crds. for the Ewald lattice
    for ( int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++ ) {
      fprintf(job, "X  %4.8f  %4.8f  %4.8f\n", EwaldCharges.Element(ichrg,0), EwaldCharges.Element(ichrg,1), EwaldCharges.Element(ichrg,2) );
    }
    
    

    // Print the Mutlipole section of the .pot file:
    fprintf(job, "@MULTIPOLES\n");
    fprintf(job, "ORDER 0\n");
    fprintf(job, "%d\n",total_atoms);
    counter=1;
    // Always print charge terms:
    for (int imon=1; imon<=NMon; imon++) {
      if (  imon != indexA && imon != indexB  ) {
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
      if ( MonomerImages[imon].GetIndex() != indexA && MonomerImages[imon].GetIndex() != indexB  ) {
	if ( MonomerImages[imon].GetUseInEmbedding() ){
	  for(int iatom=0; iatom<MonomerImages[imon].GetNumberOfAtoms() ; iatom++) {
	    Vector moments = MonomerImages[imon].GetAtom(iatom).GetMultipoleMoments().GetMoments();
	    fprintf(job, "%d  %4.5f\n",counter,moments[0]);
	    counter++;
	  }
	}
      }
    }

    for ( int ichrg=0; ichrg<EwaldCharges.GetRows(); ichrg++) {
      fprintf(job, "%d  %4.5f\n",counter,EwaldCharges.Element(ichrg,3));
      counter++;
    }


    fclose(job);

  } // end of charge embedding

}


void Dimer::ReadOrcaNMRdata(){
  string path;
  path = Params::Parameters().GetQMPath();
  string out_filename;
  char label[20];
  sprintf(label,"%d.%d.out", indexA, indexB );

  out_filename = path + "/d" + label;

  ifstream infile;
  infile.open( out_filename.c_str()  );
  if ( !infile.is_open() ) {
    printf("Dimer::ReadG09NMRdata : Cannot open file '%s'\n",
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
    printf("Dimer (%d,%d) failed, please run the monomer job\n and then re-run hmbi using 'analyze_only = true'\n", GetIndexA(), GetIndexB() );
    exit(1);
  }


  infile.close();

  // Now actually read in the shielding tensor data...
  infile.open( out_filename.c_str() );
  Matrix tmpTensor(3,3);
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.length() > 10 ) {
      //printf("DEBUG: print the sub string: %s\n", line.substr(1,34).c_str() ); fflush(stdout);
	
      // Search for MP2 NMR section 
      string match = line.substr(0,15);
      
      if ( match  == search_match ) {
	for (int skip = 1; skip <= 5; skip++){
	  getline(infile,line); // Throw away frist 19 lines...
	}

	// Loop over all the atoms and read in the tensor:
	int iatom = 0;
	for (int i=1;i<=(MonA.GetNumberOfAtoms() + MonB.GetNumberOfAtoms() );i++) {
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
	    printf("WARNING: Dimer::ReadOrcaNMRdata: zero tensor detected for dimer (%d,%d), atom %d \n", indexA, indexB ,iatom++ );
	    exit(1);
	  }
	  
	  //tmpTensor.Print("Printing the two-body tensor"); OKAY

	  if ( iatom < MonA.GetNumberOfAtoms() ) {
	    MonA.GetAtom(iatom).SetTwoBody3x3Tensor(tmpTensor);
	  } else {
	    MonB.GetAtom(iatom - MonA.GetNumberOfAtoms() ).SetTwoBody3x3Tensor(tmpTensor);
	  }

	  // Thow away lines at the end
	  for (int skip = 1; skip <= 9; skip++){
	    getline(infile,line); // Throw away frist 14 lines...
	  }

	  iatom++;
	  //fflush(stdout);
	  //printf("\nDEUBG: d(%d,%d):\n", indexA, indexB);
	  //printf("DEBUG: Dimer::ReadOrcaNMRdata() Shielding Tensor for Monomer: %d, atom %d \n", indexA, iatom  );
	  //tmpTensor.Print("DEBUG: Dimer::ReadOrcaNMRdata()"); fflush(stdout); 
	  //printf("DEBUG: iso %f\n", tmpTensor.Trace() / 3);
	}

      }
      
    }

  }

  	
  infile.close();

}



void Dimer::CreateOrcaJob(Monomer Monomers[], int NMon ) {

  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.inp",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateOrcaJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }


  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    string charge_file;
    sprintf(label,"d%d.%d.chrg",indexA,indexB);
    charge_file = label;
    
    fprintf(job,"%%pointcharges \"%s\"\n",charge_file.c_str());

  }
  fprintf(job,"\n");


  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateOrcaJob, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);
  

 
  fprintf(job,"*\n");
  fprintf(job,"\n");
  


  fclose(job);


  //Print charges if needed ( these are placed in a separate file)
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    printf("Orca with Electrostatic Embedding not functional, but code is ready to go\n");
    exit(1);

    filename = path + "/d";
    sprintf(label,"%d.%d.chrg",indexA,indexB);
    filename += label;

    FILE *job;
    if ((job = fopen(filename.c_str(),"w"))==NULL) {
      printf("Dimer::CreateOrcaJob : Cannot open file '%s'\n",filename.c_str());
      exit(1);
    }

    int num_charges = 0;
    for ( int i=1;i<=NMon;i++ ) {
      if (i != indexA && i != indexB ) {
  	if ( Monomers[i].GetUseInEmbedding() ) {
	  num_charges += Monomers[i].GetNumberOfAtoms();
	}
      }
    }

    fprintf(job,"%d\n", num_charges);
    for (int i=1;i<=NMon;i++) {
      if (i != indexA && i != indexB ) {
  	if ( Monomers[i].GetUseInEmbedding() ) {
  	  Monomers[i].PrintEmbeddingCharges(job);
  	}
      }
    }

    fclose(job);
  }
  
  //  Optionally apply mixed basis definition
  
}



void Dimer::CreateOrcaJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images ) {



  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.inp",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateOrcaJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }



  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  fprintf(job,"\n");


  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateOrcaJob, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);





  
  // Print charges if requested
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    printf("Orca with Electrostatic Embedding not functional, but code is ready to go\n");
    exit(1);
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if ( i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }

    // Now print the image charges if we are using PBCs...
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( MonomerImages[i].GetIndex() != indexA && MonomerImages[i].GetIndex() != indexB  ) {
	  if ( MonomerImages[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    }

  }



  // Optionally apply mixed basis definition


  fprintf(job,"*\n");
  fprintf(job,"\n");

  fclose(job);
  
}



void Dimer::CreateOrcaJob(Monomer Monomers[], int NMon, Monomer MonomerImages[], int NMon_images, Matrix EwaldCharges ) {

  printf("Orca with Electrostatic Embedding not functional, but code is ready to go\n");
  exit(1);

  string path;
  path = Params::Parameters().GetQMPath();

  string filename = path + "/d";
  char label[20]; 
  sprintf(label,"%d.%d.inp",indexA,indexB);

  filename += label;

  // Open the input file for writing
  FILE *job;
  if ((job = fopen(filename.c_str(),"w"))==NULL) {
    printf("Dimer::CreateOrcaJob : Cannot open file '%s'\n",filename.c_str());
    exit(1);
  }


  

  fprintf(job,"%s\n", Params::Parameters().GetOrcaHeader().c_str() );
  fprintf(job,"\n");


  int spin_state;
  if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 1 ) {
    spin_state = 1;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 1 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 1 && MonB.GetSpinState() == 2 ) {
    spin_state = 2;
  } else if ( MonA.GetSpinState() == 2 && MonB.GetSpinState() == 2 ) {
    spin_state = 1;
  } else {
    printf("ERROR: Dimer::CreateOrcaJob, didn't plan on that sort of spin state...\n");
    exit(1);
  }

  // Print charge/spin and cartesian crds
  fprintf(job,"* xyz %d %d\n", MonA.GetChargeState() + MonB.GetChargeState() , spin_state );
  MonA.PrintMonomerCartesian(job);
  MonB.PrintMonomerCartesian(job);



  
  // Print charges if requested
  if ( Params::Parameters().UseElectrostaticEmbedding()  ) {
    fprintf(job,"\n");
    for (int i=1;i<=NMon;i++) {
      if ( i != indexA && i != indexB ) {
	if ( Monomers[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(Monomers[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	  Monomers[i].PrintEmbeddingCharges(job);
	}
      }
    }

    // Now print the image charges 
    if ( Params::Parameters().IsPeriodic() ) {
      for (int i=1;i<=NMon_images; i++ ) {
	if ( MonomerImages[i].GetIndex() != indexA && MonomerImages[i].GetIndex() != indexB  ) {
	  if ( MonomerImages[i].GetUseInEmbedding() && Monomers[indexA].FindDistance(MonomerImages[i]).Element(0) <= Params::Parameters().GetElectrostaticEmbeddingCutoff()  ) {
	    MonomerImages[i].PrintEmbeddingCharges(job);
	  }
	}
      }
    }

    // Now print the charges to mimic the ewald potential
    for ( int i=0; i< EwaldCharges.GetRows(); i++ ) {
      if ( Params::Parameters().GetQMPackage() == "MOLPRO" ) {
	fprintf(job, "%10.6f,%10.6f,%10.6f,%10.6f,0\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      } else if ( Params::Parameters().GetQMPackage() == "G09" ) {
	fprintf(job, "%10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      } else if ( Params::Parameters().GetQMPackage() == "ORCA" ) {
	fprintf(job, "Q  %10.6f  %10.6f  %10.6f  %10.6f\n", EwaldCharges.Element(i,0), EwaldCharges.Element(i,1), EwaldCharges.Element(i,2), EwaldCharges.Element(i,3) );
      }
    }


  }


  fprintf(job,"*\n");
  fprintf(job,"\n");

  fclose(job);
  
}



