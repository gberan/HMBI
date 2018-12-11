#include <sys/types.h>
#include <unistd.h>
// other:
#include "supercell.h"
#include "quasiharmonic.h"
#include "atom.h"  // by Ali
#include <sys/stat.h>
#include "opt.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
using namespace hmbi_constants;

#ifdef PARALLEL
#include <mpi.h>
#endif /* PARALLEL */

// Define some tags for parallel job types
#define QCHEM_TAG 1
#define TINKER_TAG 2
#define CAMCASP_TAG 3
#define QCHEMCP_TAG 5
#define DIETAG 10

Supercell::Supercell() : Supercell_Size(NULL), Supercell_Monomers(NULL) {}

//
void Supercell::Initialize(ifstream &infile, int Nproc, string input_filename) {

  //printf("Supercell is being initialized!!!\n");
  //fflush(stdout);
  //initialize the size of the supercell first
  Supercell_Size = Params::Parameters().GetSupercellSize();
  //Supercell_Size.Print("Supercell_Size in Supercell");
  printf("supercell size 0: %f   1: %f   2: %f\n",Supercell_Size[0],Supercell_Size[1],Supercell_Size[2]);
  fflush(stdout);
  Ncells = int(Supercell_Size[0]) * int(Supercell_Size[1]) * int(Supercell_Size[2]);
  //find #atoms in unitcell and then in supercell
  Natoms = ( Cluster::cluster().GetCurrentCoordinates().GetLength()-6 )/3;
  //printf("Natoms = %d  Ncells = %d\n",Natoms,Ncells);
  //fflush(stdout);
  Supercell_Natoms = Ncells*Natoms;
  //Supercell_Natoms = int(Supercell_Size[0] * Supercell_Size[1] * Supercell_Size[2] * Natoms);

  //FIND #NMons in unitcell and then in supercell
  NMon = Cluster::cluster().GetNumberOfMonomers();
  Supercell_NMon = Ncells*NMon;
  //read unit cell dimensions ... do we need this?????
  unit_cell = new Vector[3];
  unit_cell[0] = Cluster::cluster().GetUnitCellVector(0);
  unit_cell[1] = Cluster::cluster().GetUnitCellVector(1);                       
  unit_cell[2] = Cluster::cluster().GetUnitCellVector(2);

  //now create coordinates of the supercell
  // 1st 3*Natoms coordinates are for the central unit cell ... for (not a whole lot of) convenience
  CreateSupercellCoordinates(); // needed for DFTB
  AtomicMasses = new double[Natoms];
  SetAtomicSymbols();
  CreateSupercellFractionalCoordinates();
  CreateSupercellUnitCell();
  if ( !Params::Parameters().GetThermalPropertiesUsingPhonons() ) {

    ReadReciprocalSpaceBoundaryPointsAndCreateGrid(infile);
//    CreateReciprocalSpaceSamplingGrid(infile);
  }
  else {
    CreateReciprocalSpaceSamplingGrid(infile);
  }
  // create supercell monomers
  CreateSupercellMonomerList();
  //now create full MM PBC geometry input file
  //If doing quasiharmonic approximation, create tinker job later.
  if(!Params::Parameters().NeglectManyBody() && !Params::Parameters().DoQuasiHarmonic())
    CreateSupercellTinkerJob();

  phonopy_mapping = CreatePhonopyMapping();
  supercell_init = true;
  //printf("Finished initializing Supercell\n");
  //fflush(stdout);
}


//Update Supercell if the lattice parameters or unit cell coordinates have been altered.
void Supercell::UpdateSupercellCoordinates(){

  //deleting and recreating attributates that may have altered.
  delete [] Supercell_Monomers;
  delete [] Supercell_AtomicCoordinates;
  delete [] Supercell_FractionalCoordinates;
   Monomer* Supercell_Monomers;
   Vector* Supercell_AtomicCoordinates;
   Vector* Supercell_FractionalCoordinates;

  printf("\nUpdate Supercell\n");

  //read unit cell dimensions ... do we need this?????
  unit_cell[0] = Cluster::cluster().GetUnitCellVector(0);
  unit_cell[1] = Cluster::cluster().GetUnitCellVector(1);                      
  unit_cell[2] = Cluster::cluster().GetUnitCellVector(2);    

  //now create coordinates of the supercell
  // 1st 3*Natoms coordinates are for the central unitcell ... for (not a whole lot of) convenience
  CreateSupercellCoordinates(); // not needed at this time (august 9th, 2012)
  CreateSupercellFractionalCoordinates();

  // create supercell monomers
  CreateSupercellMonomerList();
  
  //now create full MM PBC geometry input file 
  if(!Params::Parameters().NeglectManyBody())
    CreateSupercellTinkerJob();

}


Supercell::~Supercell() {                   

//  delete *Supercell_Size;
  delete [] Supercell_Monomers;
  delete [] Supercell_AtomicCoordinates;
  delete [] Supercell_FractionalCoordinates;
  delete [] AtomicMasses;
  delete [] AtomicMasses_lp;
  delete [] unit_cell;
  delete [] recip_boundary_pts;
  delete [] recip_grid_pts;
  delete [] plot_break;
  
}



void Supercell::CreateSupercellCoordinates() {

  Vector atomic_coords = Cluster::cluster().GetCurrentCoordinates();  
  //atomic_coords.Print("atomic_coords");
  int Natoms = (atomic_coords.GetLength()-6 )/3;

//  Supercell_AtomicCoordinates.Initialize(3*Supercell_Natoms);
//  Vector* Supercell_AtomicCoordinates; 


  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);


  int na_l, nb_l, nc_l, na_u, nb_u, nc_u; 
  // these are the lower/upper limits for unitcell indices ... eg: -1 to 1 for a supercell of size 3x3x3, 0 to 1 for a supercell of size 2x2x2
  //now depending on whether na, nb and nc are odd/even, we have different values for these limits
  // if even
  if (na%2==0) {na_l = na/2-1; na_u=na/2;}
  if (nb%2==0) {nb_l = nb/2-1; nb_u=nb/2;}
  if (nc%2==0) {nc_l = nc/2-1; nc_u=nc/2;}
  // if odd
  if (na%2==1) {na_l = (na-1)/2; na_u = (na-1)/2;}
  if (nb%2==1) {nb_l = (nb-1)/2; nb_u = (nb-1)/2;}
  if (nc%2==1) {nc_l = (nc-1)/2; nc_u = (nc-1)/2;}

  Supercell_AtomicCoordinates =	new Vector[na*nb*nc]; // this way of defining is awesome

  int count_sc=0;
  for (int i=-na_l;i<=na_u;i++) {
    for (int j=-nb_l;j<=nb_u;j++) {
      for (int k=-nc_l;k<=nc_u;k++) {

        Supercell_AtomicCoordinates[count_sc].Initialize(3*Natoms);
        //Imp Note for future: subtract na_l, nb_l, nc_l to from these indices to get actual indices

        //now let's find the shift and add (or translate)
      	Vector shift(3);
       	for (int n=0;n<3;n++)	{
       	  shift[n] = i*unit_cell[0][n] + j*unit_cell[1][n] + k*unit_cell[2][n];
        }

        for (int m=0;m<Natoms;m++) {
          Supercell_AtomicCoordinates[count_sc][3*m] = atomic_coords[3*m];  
          Supercell_AtomicCoordinates[count_sc][3*m+1] = atomic_coords[3*m+1];
          Supercell_AtomicCoordinates[count_sc][3*m+2] = atomic_coords[3*m+2];

          Supercell_AtomicCoordinates[count_sc][3*m] += shift[0];
          Supercell_AtomicCoordinates[count_sc][3*m+1] += shift[1];
          Supercell_AtomicCoordinates[count_sc][3*m+2] += shift[2];

        }
        count_sc++;
      }
    }
  }

}

string Supercell::GetDFTBGeomInput() {
  //printf("Attemping to create DFTB geometry input\n");
  int tmp; 
  int atomicNum;
  int lengthArray;
  bool match = false;
  string typeAtom;
  string geomInput = " TypeNames = {";
  
  //printf("geomInput = %s\n",geomInput.c_str());
  
  vector<int> check;
  vector<string> atomTypes;
  int numUniqueAtoms;

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  //printf("Correctly found the GetDFTBGeomInput in Supercell!!!\n");
  //fflush(stdout);
  Vector tmp_vec;

  //printf("Supercell_NMon = %i\n",Supercell_NMon);
  //fflush(stdout);

  for(int i=1;i<=Supercell_NMon;i++){
     for(int j=0;j<Supercell_Monomers[i].GetNumberOfAtoms();j++){
      atomicNum = Supercell_Monomers[i].GetAtom(j).GetAtomicNumber();
      typeAtom = Supercell_Monomers[i].GetAtom(j).GetSymbol();
      //printf("Monomer %d, atom %d = %d\n",i,j,atomicNum);
      //fflush(stdout);
      lengthArray=check.size();

      if(i==1 && j==0) {
        match=false;
      }
      else{
        for (int k=0; k<lengthArray; k++){
          //printf("check[k] = %d\natomicNum = %d\n",check[k],atomicNum);
          //fflush(stdout);
          if(check[k]==atomicNum){
            match=true;
          }
        }
      }
      if(!match){
        numUniqueAtoms++;
        check.push_back(atomicNum);
        atomTypes.push_back(typeAtom);
      }
      match=false;
    }
  }
  
  //Make the beginning of the atom type list
  for (int k=0; k<atomTypes.size(); k++){
    //printf("k = %i\tatomTypes[k] = %s\n",k,atomTypes[k].c_str());
    geomInput = geomInput + " \"" + atomTypes[k] + "\"";
  }

  //Now construct the actual geometry input...
  geomInput = geomInput + " }\n  TypesAndCoordinates [Angstrom] = {\n";
  //printf("geomInput = %s\n",geomInput.c_str());
  //fflush(stdout);
  //exit(0);

  for (int i=1;i<=Supercell_NMon;i++){
    for (int j=0;j<Supercell_Monomers[i].GetNumberOfAtoms();j++){
      typeAtom = Supercell_Monomers[i].GetAtom(j).GetSymbol();

      //Locate atom Type in the list of atomTypes...
      for (int k=0; k<atomTypes.size(); k++){
        //printf("check[k] = %d\natomicNum = %d\n",check[k],atomicNum);
        //fflush(stdout);
        if(atomTypes[k]==typeAtom){
          stringstream ss;
          ss << k+1;
          string str = ss.str();
          geomInput = geomInput + "    " + str + "\t";
        }
      }
      Vector tmp_xyz = Supercell_Monomers[i].GetAtom(j).GetPosition();
      //construct string with current atom Type num, and coordinates
      for( int l=0; l< 3; l++){
        stringstream ss;
        ss << tmp_xyz[l];
        string str = ss.str();
        geomInput = geomInput + str + " ";
      }
      geomInput = geomInput + "\n";
    }
  }
  geomInput = geomInput + "  }\n";
  //printf("geomInput = \n%s\n",geomInput.c_str());
  //fflush(stdout);
  //exit(0);

  return geomInput;

}

//JLM Need to check carefull that this unit cell is centered on odd numbered cells...
//Not sure how to do this yet.
void Supercell::CreateSupercellUnitCell() {

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);
  
  // Initialize the unit cell vector list
  supercell_unit_cell = new Vector[3];
  for (int i=0;i<3;i++) {
    supercell_unit_cell[i].Initialize(3);
  }

  double a = Cluster::cluster().GetUnitCellParameter("a");
  double b = Cluster::cluster().GetUnitCellParameter("b");
  double c = Cluster::cluster().GetUnitCellParameter("c");
  double alpha = Cluster::cluster().GetUnitCellParameter("alpha");
  double beta = Cluster::cluster().GetUnitCellParameter("beta");
  double gamma = Cluster::cluster().GetUnitCellParameter("gamma");

  a *= na;
  b *= nb;
  c *= nc;
  double alpha_rad = alpha*DegreesToRadians;
  double beta_rad = beta*DegreesToRadians;
  double gamma_rad = gamma*DegreesToRadians;

  printf("Supercell params: a= %d\tb = %d\tc = %d\nalpha = %d\tbeta = %d\tgamma = %d\n",a,b,c,alpha,beta,gamma);
  fflush(stdout);

  double beta_term =             
    (cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad) ) / sin(gamma_rad);
  double gamma_term =      
    sqrt(1-cos(beta_rad)*cos(beta_rad) - beta_term*beta_term);

  // v1
  supercell_unit_cell[0][0] = a;
  supercell_unit_cell[0][1] = 0;       
  supercell_unit_cell[0][2] = 0;

  // v2
  supercell_unit_cell[1][0] = b*cos(gamma_rad);
  supercell_unit_cell[1][1] = b*sin(gamma_rad);
  supercell_unit_cell[1][2] = 0;

  // v3          
  supercell_unit_cell[2][0] = c*cos(beta_rad);
  supercell_unit_cell[2][1] = c*beta_term;
  supercell_unit_cell[2][2] = c*gamma_term;

}

void Supercell::CreateSupercellMonomerList() {

/* 
   This subroutine will create the monomer list in different unit cells of the supercell.
   Once we have the monomers created for the supercell, we can deal with anything!!
   This subroutine derives heavily from Cluster::CreatePeriodicImageList().
*/

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);
  //int NMon = Cluster::cluster().GetNumberOfMonomers();
  int NDim = Cluster::cluster().GetNumberOfDimers();
  int NDim_images = Cluster::cluster().GetNumberOfDimerImages();


  // maybe we pass this in supercell.h so that we can call these monomers whenever we want
  //Supercell_Monomers = new Monomer[na*nb*nc*(NMon+1)]; // awesome but 3 years late?
  Supercell_Monomers = new Monomer[na*nb*nc*NMon+1];
  // being wasteful by copying the primitive unitcell monomers here but that's okay for the time being


  int na_l, nb_l, nc_l, na_u, nb_u, nc_u;
  // these are the lower/upper limits for unitcell indices ... eg: -1 to 1 for a supercell of size 3x3x3, 0 to 1 for a supercell of size 2x2x2
  //now depending on whether na, nb and nc are odd/even, we have different values for these limits
  // if even
  if (na%2==0) {na_l = na/2-1; na_u=na/2;}
  if (nb%2==0) {nb_l = nb/2-1; nb_u=nb/2;}
  if (nc%2==0) {nc_l = nc/2-1; nc_u=nc/2;}
  // if odd                  
  if (na%2==1) {na_l = (na-1)/2; na_u = (na-1)/2;}
  if (nb%2==1) {nb_l = (nb-1)/2; nb_u = (nb-1)/2;}
  if (nc%2==1) {nc_l = (nc-1)/2; nc_u = (nc-1)/2;}    

  int count_sc=0;
  int mon_count = 1; //added by yoni
  for (int i=-na_l;i<=na_u;i++) {
    for (int j=-nb_l;j<=nb_u;j++) {
      for (int k=-nc_l;k<=nc_u;k++) {
        // store the count_sc which corresponds to the primitive cell (0,0,0) for future
        if ((i==0) && (j==0) & (k==0)) {
          primitive_cell_index = count_sc;
        }

	Vector shift(3);
	for (int n=0;n<3;n++) {
	  shift[n] = i*unit_cell[0][n] + j*unit_cell[1][n] + k*unit_cell[2][n];
	}

        for (int m=1;m<=NMon;m++) {
	  // Determine the shift from the central cell to the image cell       
	  //Supercell_Monomers[count_sc*(NMon+1)+m] = Cluster::cluster().GetMonomer(m);            
	  //Supercell_Monomers[count_sc*(NMon+1)+m].SetReferenceMonomerIndex(m);

	  //Make monomer
	  Supercell_Monomers[mon_count] = Cluster::cluster().GetMonomer(m);

	  //Set the reference monomer
	  Supercell_Monomers[mon_count].SetReferenceMonomerIndex(m);

          // bool IsThisMonomerLocal = false; // do we need this??

	  // Now translate the monomers
          Vector new_com(3);
          new_com = Supercell_Monomers[mon_count].GetCenterOfMass();

          // Add the shift and translate the monomer
          new_com += shift;
          Supercell_Monomers[mon_count].Translate(new_com);
	  
	  //Supercell_Monomers[mon_count].PrintMonomerCartesian();

	  mon_count++;
      
        } 
	/*
        Vector shift(3);
        for (int n=0;n<3;n++) {
          shift[n] = i*unit_cell[0][n] + j*unit_cell[1][n] + k*unit_cell[2][n];
        }

        // Now translate the monomers
        for (int imon=1;imon<=NMon;imon++) {

          // bool IsThisMonomerLocal = false; // do we need this??

          Vector new_com(3);
          new_com = Supercell_Monomers[count_sc*(NMon+1)+imon].GetCenterOfMass();

          // Add the shift and translate the monomer
          new_com += shift;

          Supercell_Monomers[count_sc*(NMon+1)+imon].Translate(new_com);
        }
	  */
        count_sc++;
	
      }
    }
  }
  
}

//return the supercell mapping vector
Vector Supercell::CreatePhonopyMapping(){
  //JLM
  //Need to create a list that can map Phonopy-generated
  //atoms to hmbi-generated atoms. Ultimately
  //will be used to correctly place frequencies
  Vector map(Supercell_Natoms);
  Vector hmbi(Supercell_Natoms);
  Vector phonopy(Supercell_Natoms);
  
  //Get size of supercell
  Vector size = Params::Parameters().GetSupercellSize();
  //Separate into x,y,z
  int na,nb,nc;
  na=size[0];
  nb=size[1];
  nc=size[2];
  Vector h_na(na);
  Vector h_nb(nb);
  Vector h_nc(nc);
  Vector p_na(na);
  Vector p_nb(nb);
  Vector p_nc(nc);

  //First need to take care of the re-ordering that HMBI undergoes that Phonopy does not
  // if Even, then they are the same
  if (na%2==0) {for (int i=0; i<na; i++) {h_na[i] = i; p_na[i] = i;}}
  if (nb%2==0) {for (int i=0; i<nb; i++) {h_nb[i] = i; p_nb[i] = i;}}
  if (nc%2==0) {for (int i=0; i<nc; i++) {h_nc[i] = i; p_nc[i] = i;}}

  // if Odd, HMBI is reordered but phonopy is not
  int reorder;
  if (na%2==1) {
    reorder=na/2;
    for (int i=0; i<na; i++) {
      h_na[reorder] = i;
      p_na[i] = i;

      if(i==(na/2)) reorder = (na/2)-1;
      else if(i>=(na/2)) reorder--;
      else reorder++;
    }
  }
  if (nb%2==1) {
    reorder=nb/2;
    for (int i=0; i<nb; i++) {
      h_nb[reorder] = i;
      p_nb[i] = i;

      if(i==(nb/2)) reorder = (nb/2)-1;
      else if(i>=(nb/2)) reorder--;
      else reorder++;
    }
  }
  if (nc%2==1) {
    reorder=nc/2;
    for (int i=0; i<nc; i++) {
      h_nc[reorder] = i;
      p_nc[i] = i;

      if(i==(nc/2)) reorder = (nc/2)-1;
      else if(i>(nc/2)) reorder--;
      else reorder++;
    }
  }

  //useful variables
  string hold1,hold2,hold3;
  int hold;
  int count=0;
  int count2=0;
  int h_i, h_j, h_k;
  int p_i, p_j, p_k;

  //First create the HMBI list
  //For each supercell translate the atoms in the central unit cell
  for (int i=0; i<na; i++) {
   h_i = h_na.Element(i);
   for (int j=0; j<nb; j++) {
    h_j = h_nb.Element(j);
    for (int k=0; k<nc; k++) {
     h_k = h_nc.Element(k);
     for (int a=0; a<Natoms; a++) {
       stringstream ss;//create a stringstream
       ss << a << h_i << h_j << h_k;
       hold3 = ss.str();
       hold = atoi(hold3.c_str());
       hmbi[count] = hold;
       count++;
     }
    }
   }
  }

  //Now create the Phonopy list
  //for each atom in the central unit cell, translate into proper Supercell
  for (int a=0; a<Natoms; a++) {
   for (int i=0; i<nc; i++) {
    p_i = p_nc.Element(i);
    for (int j=0; j<nb; j++) {
     p_j = p_nb.Element(j);
     for (int k=0; k<na; k++) {
       p_k = p_na.Element(k);
       stringstream ss;//create a stringstream
       ss << a << p_k << p_j << p_i;
       hold3 = ss.str();
       hold = atoi(hold3.c_str());
       phonopy[count2] = hold;
       count2++;
     }
    }
   }
  }

  int tmp,tmp2;
  //Now iterate over the Phonopy list and find what element corresponds
  //each element in the HMBI list
  for (int s=0; s<Supercell_Natoms; s++) {  
    //The site 
    tmp = phonopy[s];
    tmp2 = hmbi.Find(tmp);
    map[s]=tmp2;
  }

  return map;
}


void Supercell::CreateSupercellTinkerJob() {

  //Path for tinker calculations
  string path = Params::Parameters().GetHessianMMPath();
  

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()){ 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
      printf("adding to path  %s\n",Quasiharmonic::quasiharmonic().GetHessianType().c_str());
   }

  // First command, change to local directory
  string xyzfile = path + "/supercell.xyz";
  string keyfile = path + "/supercell.key";
  //xyzfile = Params::Parameters().GetHessianMMPath() + "/supercell.xyz";
  //keyfile = Params::Parameters().GetHessianMMPath() + "/supercell.key";
  //printf("xyzfile = %s\n",xyzfile.c_str());
  //printf("keyfile = %s\n",keyfile.c_str());

  /* Create the xyz file */
  FILE *xyz;
  if ((xyz = fopen(xyzfile.c_str(),"w"))==NULL) {
    printf("Supercell::CreateSupercellTinkerJob : Cannot open file '%s'\n",
           keyfile.c_str());
    exit(1);  
  }                  
  PrintTinkerCartesian(xyz);
  fclose(xyz);                          

  string tinkerjobstr =  Supercell_Natoms + "\n";

  /* Create the keyfile */
  // Open the file for writing, write the Tinker rem section to it,
  // and close the file.
  FILE *key;
  if ((key = fopen(keyfile.c_str(),"w"))==NULL) {
    printf("Cluster::CreateTinkerJob : Cannot open file '%s'\n",    
           keyfile.c_str());
    exit(1);
  }
  fprintf(key,"%s\n", Params::Parameters().GetTinkerRem().c_str() );
  
  
  double a,b,c,alpha,beta,gamma;
  GetSuperCellLatticeParameters(a,b,c,alpha,beta,gamma);

  printf("Tinker lattice parameters:\n");
  printf("a = %f, b = %f, c = %f\n",a,b,c);
  printf("alpha = %f, beta = %f, gamma = %f\n",alpha,beta,gamma);

  fprintf(key,"# Periodic boundary conditions\n");         

  fprintf(key,"A-AXIS\t\t%f\nB-AXIS\t\t%f\nC-AXIS\t\t%f\n",a,b,c);
  fprintf(key,"ALPHA\t\t%f\nBETA\t\t%f\nGAMMA\t\t%f\n",alpha,beta,gamma);
  fprintf(key,"EWALD\t\tTRUE\n");
  // If we want vacuum boundary conditions for the ewald sum, set
  // EWALD-BOUNDARY flag.  For tinfoil (infinite dielectric) boundary
  // conditions, we omit this keyword.  For non-polar unit cells, it
  // doesn't matter.  But for polar ones, tinfoil boundary conditions
  // seem to work better.
  if (!Params::Parameters().TinFoilBoundaryConditions())
      fprintf(key,"EWALD-BOUNDARY\tTRUE\n");

  fclose(key);            

}


//return the supercell parameters
void Supercell::GetSuperCellLatticeParameters(double& a, double& b, double& c,double& alpha, double& beta, double& gamma){

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  a = na*unit_cell[0].Norm();
  b = nb*unit_cell[1].Norm();
  c = nc*unit_cell[2].Norm();

  alpha = RadiansToDegrees*
      acos(unit_cell[1].DotProduct( unit_cell[2] ) * nb * nc / (b*c));
  beta  = RadiansToDegrees*
      acos(unit_cell[0].DotProduct( unit_cell[2] ) * na * nc / (a*c));             
  gamma = RadiansToDegrees*
      acos(unit_cell[0].DotProduct( unit_cell[1] ) * na * nb / (a*b));  


}


void Supercell::PrintTinkerCartesian(FILE *outfile) {

  // Print tile line
  fprintf(outfile,"%d  Supercell\n", Supercell_Natoms);
  int shift = 0;

  //int NMon = Cluster::cluster().GetNumberOfMonomers();

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  int na_l, nb_l, nc_l, na_u, nb_u, nc_u;
  // these are the lower/upper limits for unitcell indices ... eg: -1 to 1 for a supercell of size 3x3x3, 0 to 1 for a supercell of size 2x$
  //now depending on whether na, nb and nc are odd/even, we have different values for these limits
  // if even
  if (na%2==0) {na_l = na/2-1; na_u=na/2;}
  if (nb%2==0) {nb_l = nb/2-1; nb_u=nb/2;}
  if (nc%2==0) {nc_l = nc/2-1; nc_u=nc/2;}
  // if odd
  if (na%2==1) {na_l = (na-1)/2; na_u = (na-1)/2;}
  if (nb%2==1) {nb_l = (nb-1)/2; nb_u = (nb-1)/2;}
  if (nc%2==1) {nc_l = (nc-1)/2; nc_u = (nc-1)/2;}

  // apply shifts to ensure proper indexing of atom numbers/connectivity
  int count_sc=0;
  int mon_count = 1;
  for (int i=-na_l;i<=na_u;i++) {
    for (int j=-nb_l;j<=nb_u;j++) {
      for (int k=-nc_l;k<=nc_u;k++) {


        for (int m=1;m<=NMon;m++) {
          //Supercell_Monomers[count_sc*(NMon+1)+m].PrintAll();
          //fflush(stdout);

          //Supercell_Monomers[count_sc*(NMon+1)+m].PrintTinkerCartesian(shift,false,outfile);
          //shift += Supercell_Monomers[count_sc*(NMon+1)+m].GetNumberOfAtoms();

	  Supercell_Monomers[mon_count].PrintTinkerCartesian(shift,false,outfile);
	  shift += Supercell_Monomers[mon_count].GetNumberOfAtoms();
	  
	  mon_count++;
        }

        count_sc++;
      }
    }
  }

}



void Supercell::RunHMBISupercellTinkerHessianJob() {

  string job_path;
//  if (Params::Parameters().GetJobTypeStr() == "phonon") {  // run hessian (vibrate) job
    job_path = Params::Parameters().GetHessianMMPath();
 

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();


    string cmd2 = "cd " + job_path;
    cmd2 += "; ";

    // Run the job
    cmd2 += "vibrate supercell.xyz 1 > supercell.freq;";
    // Remove tmp file and extra geom file created by vibrate job             
    cmd2 += "rm -f supercell.001; ";

    // Switch back to base directory      
    cmd2 += "cd " + Params::Parameters().GetBasePath();

    //printf("Executing: %s\n",cmd2.c_str());    
    system(cmd2.c_str());

//  }

}

void Supercell::SetSupercellHessian() {

  if (Params::Parameters().GetQMType() == 5){  // Quantum Espresso JLM   
    SupercellHess_QM = ReadSupercellHessian(3);
    QM_SupercellHess_Init = 1;
  }
  else if (Params::Parameters().GetQMType() == 8){  // DFTB JLM   
    SupercellHess_QM = ReadSupercellHessian(3);
    QM_SupercellHess_Init = 1;
  }
  else if (Params::Parameters().GetMMType() == 1){  // Tinker   
    //SupercellHess_MM = ReadSupercellFDHessian(2);
    SupercellHess_MM = ReadSupercellHessian(2);
    MM_SupercellHess_Init = 1;
  }
  else{
    printf("Error:Supercell::SetSupercellHessian(): \nSupercell frequency calculations have only been adapted for\n tinker (MM) or for QE (QM)\n");
    fflush(stdout);
    exit(0);
  }


}


// Main routine for reading in supercell hessian (similar to ReadHessian in cluster.C)
// type 1 = qchem style, type 2 = tinker style, type 3 = QE style
Matrix Supercell::ReadSupercellHessian(int type) {

  string path;  
  //printf("Supercell_Natoms = %d\n", Supercell_Natoms);
  Matrix supercell_hess(3*Supercell_Natoms, 3*Supercell_Natoms);

  if ( Params::Parameters().UseFullQMOnly() ){ //JLM

    path = Params::Parameters().GetHessianQMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    //printf("Now trying to read the full Quantum Espresso hessian\n");
    Natoms = Cluster::cluster().GetTotalNumberOfAtoms();
    //printf("Natoms in ReadSupercellHessian= %d\n", Natoms);

    if(Params::Parameters().UsePhononFreq()) {
      //****Note: Phonon will be an un-mass weghted hessian might need to mass-weight
      int sizeA = Supercell_Size[0];
      int sizeB = Supercell_Size[1];
      int sizeC = Supercell_Size[2];
      int primitive_cell_index;
      if((sizeA*sizeB*sizeC)%2 == 1)
        primitive_cell_index = ((sizeA*sizeB*sizeC)/2);
      else
        primitive_cell_index = 0;

      // Set up the filename, with the full path.
      for (int a=1;a<=sizeA;a++) {
        for (int b=1;b<=sizeB;b++) {
          for (int c=1;c<=sizeC;c++) {

            stringstream ss;
            ss << a-1 << b-1 << c-1;
            string str = ss.str();
            string filename = path + "/dyn" + str+".G";
            //printf("filename: %s\n",filename.c_str());
            
            // Open the force file
            ifstream infile;
            infile.open( filename.c_str() );
            if ( !infile.is_open() ) {
              printf("Supercell::ReadHessian : Cannot open file '%s'\n",filename.c_str());
              exit(1);
            }
            
            string line; 
            while ( !infile.eof() ) {
              getline(infile,line);
              //printf("%s\n",line.c_str()); 
              // Search for the dynamical matrix  
              //Our final Hessian should be in units of hartrees/bohr^2
              if ( line.find("Dynamical  Matrix in cartesian axes") != -1 ) {
                getline(infile,line); // throw away header line
                getline(infile,line); // skip blank line
                getline(infile,line); // skip qpoints line
                for (int i=0;i<Natoms;i++) {
                  for (int j=0;j<Natoms;j++) {
                    getline(infile,line); //Always skip the atom-atom interaction line
                    for (int m=0;m<3;m++) {
                      getline(infile,line); //Always skip the atom-atom interaction line
                      istringstream iss(line);
                      for (int k=0;k<3;k++) {
                        string imag;
                        //iss >> hess.Element(3*i+m,3*j+k); // Store the hessian elements
                        iss >> supercell_hess.Element(3*Natoms*primitive_cell_index+3*i+m,3*Natoms*primitive_cell_index+3*j+k);
                        iss >> imag;
                        //printf("hess.Element(%d,%d): %f\nimag = %s\n",3*i+m,3*j+k,hess.Element(3*i+m,3*j+k),imag.c_str());
                      }
                    }
                  }
                }
                break;
              }
            }
            
            infile.close();  // Close the file
            if((primitive_cell_index == (sizeA*sizeB*sizeC)-1) && ((sizeA*sizeB*sizeC) != (a*b*c)) ){
              primitive_cell_index=0;
            }
            else{
              primitive_cell_index++;
            }
          }
        }
      }

      //Cluster::cluster().PrintHessian("Before scaling",supercell_hess);
      //fflush(stdout);

      //Our final Hessian should be in units of hartrees/bohr^2
      //Currently in Ry/bohr^2 so just divide by 2
      supercell_hess.Scale(0.5);

    }
    else if(Params::Parameters().UseDFTBFreq()) {

      printf("Correctly found the DFTB Freq read-in in Supercell\n");
      fflush(stdout);

      //Matrix hess(3*Supercell_Natoms, 3*Supercell_Natoms);
      //hess.Set();
      //hess.Print("testing hessian set");

      string filename;
      // Set up the filename, with the full path.
      filename = path + "/hessian_phn.out";

      // Open the file
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
        printf("Cluster::ReadDFTBHessian : Cannot open file '%s'\n",filename.c_str());
        exit(1);
      }

      int count=0;
      string line; 
      while ( !infile.eof() ) {
        Vector tmp(4);
        getline(infile,line);
        istringstream iss(line);
        iss >> tmp[0];
        iss >> tmp[1];
        iss >> tmp[2];
        iss >> tmp[3];

        for (int i=0;i<3*Supercell_Natoms;i++) {
          for (int j=0;j<3*Supercell_Natoms;j++) {
            //getline(infile,line); //Always skip the atom-atom interaction line
            if(count==4){
              //printf("Correctly trying to progress the line!\n");
              //printf("i = %i\tj = %i\n",i,j);
              //fflush(stdout);
              getline(infile,line); //Always skip the atom-atom interaction line
              //printf("Got new line!: %s\n",line.c_str());
              //fflush(stdout);
              istringstream iss(line);
              iss >> tmp[0];
              iss >> tmp[1];
              iss >> tmp[2];
              iss >> tmp[3];
              count=0;
            }
            //printf("line %s\n",line.c_str());
            //fflush(stdout);

            //iss >> hess.Element(i,j); // Store the hessian elements
            supercell_hess.Element(i,j) = tmp[count]; // Store the hessian elements
            count++;

            //printf("hess.Element(%d,%d): %f\n",i,j,hess.Element(i,j));
          }
        }
        break;

      }

      infile.close();  // Close the file

      //supercell_hess.Print("Supercell set");

    //PrintHessian("Before scaling",hess);
    //fflush(stdout);
    //exit(0);
    //Convert units to hartree/bohr^2
    //Currently in Ry/bohr^2 so just divide by 2
    //hess.Scale(BohrToAng*BohrToAng*EVToHartrees);
    //hess.Scale(0.5);

  //un-Mass weighted Hessian (i.e. the HMBI Hessian)
  /*Hess_HMBI = hess;
  //PrintHessian("Before mass scaling",hess);
  //printf("\n");
  for (int i=0;i<Natoms;i++) {
    double m1 = AtomicMasses[i];
    for (int j=0; j<Natoms;j++) {
      double m2 = AtomicMasses[j];
      double Gij = sqrt(m1*m2);
      for (int p=3*i;p<3*(i+1);p++) {
        for (int q=3*j;q<3*(j+1);q++) {
          hess(p,q) /= Gij;
        }       
      }
    }
  }

  //Mass-weighted Hessian
  Hess_QM = hess;*/

      //exit(0);
    }
    else if(Params::Parameters().UsePhonopyFreq()) {

      // Set up the filename, with the full path.
      string filename = path + "/FORCE_CONSTANTS_SPG";

      // Open the force file
      ifstream infile;
      infile.open( filename.c_str() );
      if ( !infile.is_open() ) {
        printf("Supercell::ReadHessian : Cannot open file '%s'\n",filename.c_str());
        exit(1);
      }
    
      string line; 
      int dim1, dim2;
      string tmp1,tmp2,tmp3;

      //Get dimensions of the force_constants
      getline(infile,line);
      istringstream iss(line);
      iss >> tmp1;
      iss >> tmp2;
      dim1=atoi(tmp1.c_str());
      dim2=atoi(tmp2.c_str());

      //Need to address ordering differences between phonopy and HMBI
      //Phonopy orders their atom displacements in x,y,z
      //HMBI orders our atom displacements from z,y,x
      //This always needs to be addressed
      Vector mapping=phonopy_mapping;

      //i goes over dim1 in FORCE_CONSTANTS_SPG
      //j goes over dim2 in FORCE_CONSTANTS_SPG
      //k handles the three elements of each line
      for (int i=0;i<dim1;i++) {
        for (int j=0;j<dim2;j++) {
        //for (int j=0;j<2;j+=1) {
          getline(infile,line);

          int tmpi=mapping.Element(i);
          int tmpj=mapping.Element(j);
          //printf("i: %i\t j: %i\n",i,j);
          //printf("mapping i: %i\t j: %i\n",tmpi,tmpj);
          //fflush(stdout);  
 
          for (int k=0;k<3;k++) { //Leave part unchanged; Just need to change mapping
            getline(infile,line);
            istringstream iss(line);
            //This is the line in FORCE_CONSTANTS_SPG
            //Now need to correctly map this position to
            //Where it should be in HMBI

            iss >> supercell_hess.Element(3*tmpi+k,3*tmpj);
            iss >> supercell_hess.Element(3*tmpi+k,3*tmpj+1);
            iss >> supercell_hess.Element(3*tmpi+k,3*tmpj+2);

            //printf("hess.Element(%d,%d):\t%f\n",3*tmpi+k,3*tmpj,supercell_hess.Element(3*tmpi+k,3*tmpj));
            //printf("hess.Element(%d,%d):\t%f\n",3*tmpi+k,3*tmpj+1,supercell_hess.Element(3*tmpi+k,3*tmpj+1));
            //printf("hess.Element(%d,%d):\t%f\n",3*tmpi+k,3*tmpj+2,supercell_hess.Element(3*tmpi+k,3*tmpj+2));
            //fflush(stdout);

            //printf("hess.Element(%d,%d):\t%f\n",i+k,j,supercell_hess.Element(i+k,j));
            //printf("hess.Element(%d,%d):\t%f\n",i+k,j+1,supercell_hess.Element(i+k,j+1));
            //printf("hess.Element(%d,%d):\t%f\n",i+k,j+2,supercell_hess.Element(i+k,j+2));
          }
        }
      }

      infile.close();  // Close the file

      //Cluster::cluster().PrintHessian("Before scaling",supercell_hess);
      //fflush(stdout);

      //Our final Hessian should be in units of hartrees/bohr^2
      //Currently in Ry/bohr^2 so just divide by 2
      supercell_hess.Scale(0.5);
    }
    else {
      printf("Unrecognized QM Freq type in. Exiting...");
      fflush(stdout);
      exit(0);
    }

    // FD hessians could be non-symmetric due to numerical noise, so symmetrize them
    for (int i=0;i<3*Supercell_Natoms;i++) {
      for (int j=i;j<3*Supercell_Natoms;j++) {
        supercell_hess.Element(i,j) += supercell_hess.Element(j,i); 
        supercell_hess.Element(i,j) *= 0.5;
      }
    }

    for (int i=1;i<3*Supercell_Natoms;i++) {
      for (int j=0;j<i;j++) {
        supercell_hess.Element(i,j) *= 0.0; 
        supercell_hess.Element(i,j) += supercell_hess.Element(j,i);
      }
    }

   
    //Cluster::cluster().PrintHessian("Total supercell",supercell_hess);
    //fflush(stdout);
  }
  else if ( Params::Parameters().NeglectManyBody() ){
    supercell_hess.Set();
  }
  else {

    path = Params::Parameters().GetHessianMMPath();

    //Path of the quasiharmonic calculations
    if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();
  
    string filename;
    // Set up the filename, with the full path.  File is 'supercell.freq'
    if (type == 2) // Tinker MM job
      filename = path + "/supercell.freq";  
    else{
      printf("Error:Supercell::ReadSupercellHessian(): Supercell MM frequency calculations has only been adapted for tinker\n");
      exit(0);
    }

    // Open the force file                
    ifstream infile;               
    infile.open( filename.c_str() );                 
    if ( !infile.is_open() ) {
      printf("Supercell::ReadHessian : Cannot open file '%s'\n",
             filename.c_str());            
      exit(1);                
    }

    if (type == 2) { // look in the tinker .freq file      
      string line;
      while ( !infile.eof() ) {                  
        getline(infile,line);
        // Search for the SCF hessian
        if ( line.substr(0,40)==" Hessian Matrix (in hartrees/Bohr/Bohr):" ) {
   	  //getline(infile,line); // throw away header line
          for (int i=0;i<3*Supercell_Natoms;i++) {
            for (int j=0;j<3*Supercell_Natoms;j++) {
              getline(infile,line);
              istringstream iss(line);
              string tmp1,tmp2;       
              iss >> tmp1; // throw away the matrix row index
              iss >> tmp2; // throw away the matrix column index
              iss >> supercell_hess.Element(i,j); // Store the hessian elements
            }   
          }
          break;
        }
      }
    }
    infile.close();  // Close the file

  }

  //Cluster::cluster().PrintHessian("Supercell Hessian", supercell_hess);
  return supercell_hess;
}



void Supercell::ComputeHMBISupercellHessian() {
 

 double full_pbc_term=1.0; //for debugging ... turning off/on of the full pbc mm term
 double scafac2=1.0; //for debugging ... turning off/on of the qm term
 double scafac=1.0; //for debugging ... turning off/on of the 1-/2-body mm term

 SetSupercellHessian();
 if (Params::Parameters().UseFullQMOnly()) {  // Contribution from QM hessian
   hess = SupercellHess_QM;
 }
 else { // Contribution from full cluster MM gradient
  hess = SupercellHess_MM;
  hess.Scale(full_pbc_term);  


  // now correcting only the terms involving the central unit cell ... always remember this in the future that this is not exactly the supercell hmbi hessian
  // as only some of the terms related to the primitive cell are corrected
  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  int na_l, nb_l, nc_l, na_u, nb_u, nc_u;
  // these are the lower/upper limits for unitcell indices ... eg: -1 to 1 for a supercell of size 3x3x3, 0 to 1 for a supercell of size 2x2x2
  //now depending on whether na, nb and nc are odd/even, we have different values for these limits
  // if even
  if (na%2==0) {na_l = na/2-1; na_u=na/2;}                  
  if (nb%2==0) {nb_l = nb/2-1; nb_u=nb/2;}
  if (nc%2==0) {nc_l = nc/2-1; nc_u=nc/2;}                  
  // if odd
  if (na%2==1) {na_l = (na-1)/2; na_u = (na-1)/2;}
  if (nb%2==1) {nb_l = (nb-1)/2; nb_u = (nb-1)/2;}     
  if (nc%2==1) {nc_l = (nc-1)/2; nc_u = (nc-1)/2;}               


  //int NMon = Cluster::cluster().GetNumberOfMonomers();
  int NDim = Cluster::cluster().GetNumberOfDimers();
  int NDim_images = Cluster::cluster().GetNumberOfDimerImages();

  // Create a list of the starting index for each monomer in the
  // full hessian.  Like the Monomers, The key indexes count from 1->NMon.
  // On the other hand, in the hessian, the indexing starts at (0,0).
  //int *row_key = new int[na*nb*nc*(NMon+1)];
  //for (int i=0;i<na*nb*nc*(NMon+1);i++) {
  //  row_key[i] *= 0;
  //}
  int row_key[Supercell_NMon];
  //int *key = new int[Supercell_NMon];


  //int *key = new int[na*nb*nc*(NMon+1)];
  //for (int j=0;j<na*nb*nc;j++) {
    //row_key[j*(NMon+1)+1] = j*3*Natoms;
  //key[0] = 0; // set the first one by hand
  row_key[0] = 0;
  row_key[1] = 3*Supercell_Monomers[1].GetNumberOfAtoms();


  printf("Supercell_Mon = %i\n",Supercell_NMon);
  //printf("row_key = [ 0 %i",row_key[1]);
  for (int i=2;i<=Supercell_NMon;i++) {
    row_key[i] = row_key[i-1] + 3*Supercell_Monomers[i-1].GetNumberOfAtoms();
    //printf(" %i",row_key[i]);
    //fflush(stdout);
    }
  //printf("]\n");

  //printf("Hessian (351,351) before mon contr =%f\n",hess.Element(351,351));

  // 1-body contributions
  // debug status: success
  //printf("Starting one body contributions\n");
  for (int i=1;i<=NMon;i++) {
    int Na = Cluster::cluster().GetMonomer(i).GetNumberOfAtoms();
    //int start = key[i];
    //Cluster::cluster().GetMonomer(i).PrintAll();

    int start = row_key[primitive_cell_index*NMon+i-1];
    //printf("primitive_cell_index = %i\n",primitive_cell_index);
    //printf("row_key = %i\n",start);
  

    for (int j=0;j<3*Na;j++) {
      for (int k=0; k<3*Na;k++) {

	hess.Element(start+j,start+k) += scafac2*Cluster::cluster().GetMonomer(i).GetQMHessian().Element(j,k);
	
	hess.Element(start+j,start+k) -= scafac*Cluster::cluster().GetMonomer(i).GetMMHessian().Element(j,k);


	//if(row_key[primitive_cell_index*(NMon+1)+i]+j == 351 &&
	//  row_key[primitive_cell_index*(NMon+1)+i]+k == 351){
	//  printf("Initial Hessian = %f\n",hess.Element(row_key[primitive_cell_index*(NMon+1)+i]+j,row_key[primitive_cell_index*(NMon+1)+i]+k));
	//}

	  //if(row_key[primitive_cell_index*(NMon+1)+i]+j == 351 &&
	  //   row_key[primitive_cell_index*(NMon+1)+i]+k == 351){
	  //  printf("Super Hess Monomer %i\n",i);
	  //  printf("element (%i,%i)\n",j,k);
	  //  printf("pimitive_cell_index = %i\n",primitive_cell_index);
	  //  printf("QM Monomer Hess = %f\n",Cluster::cluster().GetMonomer(i).GetQMHessian().Element(j,k));
	  //  printf("Hess = %f\n",hess.Element(row_key[primitive_cell_index*(NMon+1)+i]+j,row_key[primitive_cell_index*(NMon+1)+i]+k));

	  //}

	  //if(row_key[primitive_cell_index*(NMon+1)+i]+j == 351 &&
	  //   row_key[primitive_cell_index*(NMon+1)+i]+k == 351){
	  //  printf("MM Monomer Hess = %f\n",Cluster::cluster().GetMonomer(i).GetMMHessian().Element(j,k));
	  //  printf("Hess = %f\n",hess.Element(row_key[primitive_cell_index*(NMon+1)+i]+j,row_key[primitive_cell_index*(NMon+1)+i]+k));
	 // }
      } 
    }
  }
  //printf("Finished one body contributions. Beginning 2. NDim = %i\n",NDim);
  //printf("Before dimers element 957,1647 %f\n",hess.Element(957,1647));

  // 2-body interaction contributions    
  for (int i=1;i<=NDim;i++) {
    Dimer Dimers( Cluster::cluster().GetDimer(i) );
    int Na = Dimers.GetMonomerA().GetNumberOfAtoms();
    int Nb = Dimers.GetMonomerB().GetNumberOfAtoms();
 
    int index3 = Dimers.GetIndexA();        
    int index4 = Dimers.GetIndexB();
    Monomer Mon3 = Cluster::cluster().GetMonomer(index3);
    Monomer Mon4 = Cluster::cluster().GetMonomer(index4);

    //Matrix Rot3 = Mon3.GetRotationMatrix();
    //Matrix Rot4 = Mon4.GetRotationMatrix();


    //int indexA = Dimers.GetIndexA();
    //int indexB = Dimers.GetIndexB();
    
    //Monomer MonA = Dimers.GetMonomerA();//( Cluster::cluster().GetMonomer(indexA) );
    //Monomer MonB = Dimers.GetMonomerB();//( Cluster::cluster().GetMonomer(indexB) );
    //Monomer MonA = Cluster::cluster().GetMonomer(indexA);
    //Monomer MonB = Cluster::cluster().GetMonomer(indexB);


    //MonA.PrintAll();
    //MonB.PrintAll();
    
    //Dimers.PrintAll();
    
    //int startA = row_key[primitive_cell_index*NMon+indexA-1];
    //int startB = row_key[primitive_cell_index*NMon+indexB-1];

    double damp;
    double c1=0, c0=0;
    if ( Params::Parameters().DoLocal2BodyTruncation() || Params::Parameters().IsPeriodic()  ) {
      c0 = Params::Parameters().GetLocalCutoff(0);  
      c1 = Params::Parameters().GetLocalCutoff(1);
      damp = Dimers.GetDampingFactor(c0,c1);
      //cout << "daya damp\n";
    }
    else {
      damp = 1.0;
    }
   
    fflush(stdout);

    //if ( Params::Parameters().DoLocal2BodyTruncation() || Params::Parameters().IsPeriodic()  ) {
    if(Dimers.GetSymmetryFactor() && damp >= pow(10.0,-6) ){
      
      //if(damp < fabs(0.000001)){
      //	Matrix tmp(3*Na+3*Nb,3*Na+3*Nb);
      //	Dimers.SetQMHessian(tmp);
      //	Dimers.SetMMHessian(tmp);
      //}
      
      //}
      //else {
      //damp = 1.0;
      //}
      
      //printf("\nd(%i,%i)\n",indexA,indexB);
      //printf("startA = %i startB = %i\n",startA,startB);
      //Dimers.GetQMHessian().PrintHessian("Dimer Hessian");
      //fflush(stdout);
      
      //printf("MonA\n");
      
      //Getting list of all dimers that are symmetrical to this dimer. Rotating the
      //elements and adding them to the right entry.
      //Each entry in the list repressent a monomer so two entries are needed for a dimer.
      vector<int> SymmetryList = Dimers.GetSymmetryList();
      vector<int> PeriodicSymmetryList = Dimers.GetPeriodicSymmetryList();
      int NumberOfSymmetricalDimers = (SymmetryList.size() + PeriodicSymmetryList.size())/2;
      
      for(int j=0; j<NumberOfSymmetricalDimers;j++){
	//index of the monomers in the symmetrical equivalent dimer
	int index1;
	int index2;
	Monomer MonA;
	Monomer MonB;
	Matrix Rot;
	Vector SymAtomList;

        //Is MonA symmetrical to the first monomer in the dimer(Mon3)
	bool FirstAInDim = false; 
	
	//The symmetrical equivalent dimers that are periodic images have their contribution to the hessian scaled by 1/2
	float per_sca;

	int startA;
	int startB;
	int startA_image;//equal to startA if dimer on symetry dimer list is a non-periodic dimer
	int startB_image;// equal to startB if dimer on symetry dimer list is a non-periodic dimer

	//Piece arising from first monomer. 
	//Symmetrical equivalent dimer is a not periodic image
	if(j<SymmetryList.size()/2){
	  //2*j is the first molecule in the dimer. 2*j+1 is the second monomer in the dimer
	  //printf("j=%i\n",j);
	  index1 = SymmetryList[2*j];
	  index2 = SymmetryList[2*j+1];
	  //Mon1 = Cluster::cluster().GetMonomer(index1);
	  //Mon2 = Cluster::cluster().GetMonomer(index2);
	  per_sca = 1.0;

	//Rotation under symmetry H12 = R13*H34*R24' where R13=R24 and R13=R1'*R3
	  Rot = Dimers.GetRotationList()[j];
          Rot.Transpose();
	  SymAtomList = Dimers.GetAtomEquivalency()[j];


          //printf("d(%i,%i)\n",index3,index4);
          //SymAtomList.Print("SymAtomList");
          //fflush(stdout);
	  
	  //Determines if Mon3 is symmetrical to Mon1 or Mon2 by determining
	  //Need to know if Mon3 and Mon4 have different number to atoms
	  for(int atom = 0; atom < Mon3.GetNumberOfAtoms();atom++){
	    if(SymAtomList[atom] == 1) FirstAInDim = true;
	  }
	  
	  //Determining which monomer is symetrical to Mon3
	  if(FirstAInDim){
	    MonA = Cluster::cluster().GetMonomer(index1);
	    MonB = Cluster::cluster().GetMonomer(index2);
	    startA = row_key[primitive_cell_index*NMon+index1-1];
	    startB = row_key[primitive_cell_index*NMon+index2-1];
	    
	  }else{
	    MonA = Cluster::cluster().GetMonomer(index2);
	    MonB = Cluster::cluster().GetMonomer(index1);	   
	    startA = row_key[primitive_cell_index*NMon+index2-1];
	    startB = row_key[primitive_cell_index*NMon+index1-1];
	  }
	  startA_image = startA;
	  startB_image = startB;

	}
	//Symmetrical equivalent dimer is a periodic image
	else{

	  //2*k is the first molecule in the dimer. 2*k+1 is the second monomer in the dimer
	  int k = j - SymmetryList.size()/2;
	  //printf("k=%i j=%i\n",k,j);
	  index1 = PeriodicSymmetryList[2*k];
	  index2 = PeriodicSymmetryList[2*k+1];
	  //Mon1 = Cluster::cluster().GetMonomer(index1);
	  //Mon2 = Cluster::cluster().GetMonomer(index2);
	  per_sca = 0.5;
	 

	  int cell_index[3];
	  cell_index[0] = Dimers.GetSymmetricalImageCell()[3*k];
	  cell_index[1] = Dimers.GetSymmetricalImageCell()[3*k+1];
	  cell_index[2] = Dimers.GetSymmetricalImageCell()[3*k+2];
	  //uci = unitcell index for unitcell from inversion
	  // uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
	  int uci = GetCellIndex(na,nb,nc,cell_index[0],cell_index[1],cell_index[2]); //uci = unitcell index
	  int uci_mirror = GetCellIndex(na,nb,nc,-1*cell_index[0],-1*cell_index[1],-1*cell_index[2]);	

	//Rotation under symmetry H12 = R13*H34*R24' where R13=R24 and R13=R1'*R3
 	  Rot = Dimers.GetPeriodicRotationList()[k];
	  Rot.Transpose();
	  SymAtomList = Dimers.GetPeriodicAtomEquivalency()[k];

	  //Determines if Mon3 is symmetrical to Mon1 or Mon2 by determining
	  //Need to know if Mon3 and Mon4 have different number to atoms

	  for(int atom = 0; atom < Mon3.GetNumberOfAtoms();atom++){
	    if(SymAtomList[atom] == 1) FirstAInDim = true;
	    //printf("atom = %i maxatom = %i index = %i SymAtomList = %i \n",atom,Monomers[index3].GetNumberOfAtoms(),
	    //   index3,(int)SymAtomList[atom]);
	  }
	  if(FirstAInDim){
	    MonA = Cluster::cluster().GetMonomer(index1);
	    MonB = Cluster::cluster().GetMonomer(index2);
	    startA = row_key[primitive_cell_index*NMon+index1-1];
	    startB = row_key[primitive_cell_index*NMon+index2-1];
	    startA_image = row_key[uci_mirror*NMon+index1-1];	  
	    startB_image = row_key[uci*NMon+index2-1];
	  }else{
	    MonA = Cluster::cluster().GetMonomer(index2);
	    MonB = Cluster::cluster().GetMonomer(index1);
	    startA = row_key[primitive_cell_index*NMon+index2-1];   
	    startB = row_key[primitive_cell_index*NMon+index1-1]; 
	    startA_image = row_key[uci*NMon+index2-1];
	    startB_image = row_key[uci_mirror*NMon+index1-1];
	  }
	  
	  /*
	  //Mon1 is outside the unit cell
	  if(Dimers.GetMonBList()[k] == 0){
	  //uci = unitcell index for unitcell from inversion
	  // uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
	  startA = row_key[primitive_cell_index*NMon+index1-1];
	    startB = row_key[primitive_cell_index*NMon+index2-1];
	    startA_image = row_key[uci*NMon+index1-1];	  
	    startB_image = row_key[uci_mirror*NMon+index2-1];
	    }
	    //Mon2 is outside the unit cell
	    else{
	    //uci = unitcell index for unitcell from inversion
	    // uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
	    startA = row_key[primitive_cell_index*NMon+index1-1];
	    startB = row_key[primitive_cell_index*NMon+index2-1];	  
	    startA_image = row_key[uci_mirror*NMon+index1-1];	  
	    startB_image = row_key[uci*NMon+index2-1];
	    }
	  */
	  
	  
	}  

	//Piece arising from first monomer. 
	for(int iatom=0;iatom<Na;iatom++){
	  
	  //determine which atom iatom is symmetrical equivalent on MonA	  
	  int Symi = (int) SymAtomList[iatom] - 1;
	  if(!FirstAInDim){
	    Symi-= Nb;
	  }

	  for(int jatom=0;jatom<Na;jatom++){
	    
	    //determine which atom jatom is symmetrical equivalent on MonA
	    int Symj = (int) SymAtomList[jatom] - 1;
	    if(!FirstAInDim)
	      Symj -= Nb;
	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*Dimers.GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4) - scafac*Dimers.GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 & fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		   /*   printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		      printf("monA\n");
		      printf(" damp = %f\n",damp);
		      printf(" per_sca = %f\n",per_sca);
		      printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		      printf(" Unrotated QM = %15.9f\n", Dimers.GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" Unrotated MM = %15.9f\n", Dimers.GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		      printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		      printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		      printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2)  );
		    }	*/	   


		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca
		      *Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*Dimers.GetQMGradient()[3*jatom+xyz4]-scafac*Dimers.GetMMGradient()[3*jatom+xyz4]);
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //   printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot(xyz3,xyz1) = %f\n",Rot(xyz1,xyz3));
		    // printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //   printf(" Grad Damp iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		    //  printf(" Grad QM Dimer jatom = %f\n",Dimers.GetQMGradient()[3*jatom+xyz4]);
		    //  printf(" Grad MM Dimer jatom = %f\n",Dimers.GetMMGradient()[3*jatom+xyz4]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}	   

		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca
		      *Rot(xyz2,xyz4)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*Dimers.GetQMGradient()[3*iatom+xyz3]-scafac*Dimers.GetMMGradient()[3*iatom+xyz3]);
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*jatom+xyz4,3*iatom+xyz3);
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Grad Damp Symi = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
		    //  printf(" Grad QM Dimer Symj = %f\n",Dimers.GetQMGradient()[3*iatom+xyz3]);
		    //  printf(" Grad MM Dimer Symj = %f\n",Dimers.GetMMGradient()[3*iatom+xyz3]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}	   

		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
		      *(scafac2*Dimers.GetQMIntEnergy()-scafac*Dimers.GetMMIntEnergy());
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4));
		    // printf(" Int QM Dimer = %f\n",Dimers.GetQMIntEnergy());
		    //  printf(" Int MM Dimer = %f\n",Dimers.GetMMIntEnergy());
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}	   	  
		    
		    //}
		    //if (!Params::Parameters().Do_fdTinkerHessian() ) { 
		    
		    
		  }//end of loop over xyz4

		  //Rotating gradient of QM monomers contributation
		  
		  //QM Monomers:already subtracted if counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise())
		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) -= per_sca*scafac2*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonA.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonA.GetQMGradient()[3*Symi+xyz1]);
		  
		  //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001)  ){
		  //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //  printf("monA\n");
		  //  printf(" per_sca = %f\n",per_sca);
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element - (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer Symj = %f\n",MonA.GetQMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" Grad QM Monomer Symi  = %f\n",MonA.GetQMGradient()[3*Symi+xyz1]);
		  // printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		  //}	 
		  
		  //MM Monomers:already subtrated if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		    hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca*scafac*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonA.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonA.GetMMGradient()[3*Symi+xyz1]);
		  
		  // if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //   printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //  printf("monA\n");
		  //  printf(" per_sca = %f\n",per_sca);
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element - (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp Symj = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" Grad MM Monomer i = %f\n",MonA.GetMMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp iSym = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" Grad MM Monomer jatom = %f\n",MonA.GetMMGradient()[3*Symi+xyz2]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		  //}	
		    
		}//end of loop over xyz3
		
		//Including Monomer contribation
		//Don't have to be rotated so outside xyz3 and xyz4 loop
		  
		//QM Monomers
		if(!Params::Parameters().DoCounterpoise())
		  hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) -= per_sca*damp*scafac2*
		    MonA.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);

		//if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 ){
		//  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		//  printf("monA\n");
		//  printf(" per_sca = %f\n",per_sca);
		//  printf(" Element = (%i,%i)\n",3*iatom+xyz1,3*jatom+xyz2);
		//  printf(" damp = %f\n",damp);
		//  printf(" Mon QM Hess = %f\n",MonA.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		//  printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		//}
	    
		  //MM Monomers
		if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		  hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) += per_sca*damp*scafac*
		    MonA.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		
		//if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 ){
		//  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		//  printf("monA\n");
		//  printf(" per_sca = %f\n",per_sca);
		//  printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		//  printf(" damp = %f\n",damp);
		//  printf(" Mon MM Hess = %f\n",MonA.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		//  printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		//}
		
		//}
		
	      }//end of loop over xyz2
	    }//end of loop over xyz1
	  }//end of loop over jatom
	}//end of loop over iatom

	//Piece arising from second monomer. 
	for(int iatom=0;iatom<Nb;iatom++){
	  
	  //determine which atom iatom is symmetrical equivalent on MonB
	  int Symi = (int) SymAtomList[Na+iatom] - 1;
	  if(FirstAInDim){
	    Symi-= Na;
	  }
	  
	  //determine which atom jatom is symmetrical equivalent on MonB
	  for(int jatom=0;jatom<Nb;jatom++){
	    
	    int Symj = (int) SymAtomList[Na+jatom] - 1;
	    if(FirstAInDim){
	      Symj-= Na;
	    }


	    // piece arising from second monomer
	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //if ( !Params::Parameters().TinkerDebug() ) {
  
		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) += per_sca*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*Dimers.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4)
                        -scafac*Dimers.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //   printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Unrotated QM = %15.9f\n", Dimers.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", Dimers.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2)  );
		    //}	
	    
		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) +=  per_sca
		      *Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*Dimers.GetQMGradient()[3*Na+3*jatom+xyz4]-scafac*Dimers.GetMMGradient()[3*Na+3*jatom+xyz4]);
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Grad Damp iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		    //  printf(" Grad QM Dimer jatom = %f\n",Dimers.GetQMGradient()[3*Na+3*jatom+xyz4]);
		    //  printf(" Grad MM Dimer jatom = %f\n",Dimers.GetMMGradient()[3*Na+3*jatom+xyz4]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    // }
	 
		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) += per_sca
		      *Rot(xyz2,xyz4)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*Dimers.GetQMGradient()[3*Na+3*iatom+xyz3]-scafac*Dimers.GetMMGradient()[3*Na+3*iatom+xyz3]);
		    

		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz4,3*Na+3*iatom+xyz3);
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Grad Damp jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		    //  printf(" Grad QM Dimer iatom = %f\n",Dimers.GetQMGradient()[3*Na+3*iatom+xyz4]);
		    //  printf(" Grad MM Dimer iatom = %f\n",Dimers.GetMMGradient()[3*Na+3*iatom+xyz4]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    //}	 	   

		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) +=  per_sca*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4)
		      *(scafac2*Dimers.GetQMIntEnergy()-scafac*Dimers.GetMMIntEnergy());

		    //  if(startA+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //    printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //    printf("monA\n");
		    //    printf(" per_sca = %f\n",per_sca);
		    //	printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //    printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //    printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //    printf(" Damp Hess = %f\n",Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //    printf(" Int QM  = %f\n",Dimers.GetQMIntEnergy());
		    //    printf(" Int MM = %f\n",Dimers.GetMMIntEnergy());
		    //    printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    //  }
	      //}
		    
		    
		    //if (!Params::Parameters().Do_fdTinkerHessian() ) {   
		    

		  }//end of xyz4

		  //Rotating gradient of monomers contributation

		  //QM Monomers:already subtracted if counterpoise corrected
                  if(!Params::Parameters().DoCounterpoise())
		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) -= per_sca*scafac2*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonB.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonB.GetQMGradient()[3*Symi+xyz1]);

		  //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //  printf("monB\n");
		  //  printf(" per_sca = %f\n",per_sca);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element - (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer Symj = %f\n",MonB.GetQMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //  printf(" Grad QM Monomer Symi = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		  //}	 

		  //MM Monomers:already subtracted if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		    hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) += per_sca*scafac*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonB.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonB.GetMMGradient()[3*Symi+xyz1]);
		  
		  //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //  printf("monB\n");
		  //  printf(" per_sca = %f\n",per_sca);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		   // printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Grad Damp iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" Grad MM Monomer Symj = %f\n",MonB.GetMMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer jatom = %f\n",MonB.GetMMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		  //}
		  
		}//end of xyz3	
		

		 //Including Monomer contributation
		 //Don't have to be rotated so outside xyz3 and xyz4 loop

		 //QM Monomers
		 if(!Params::Parameters().DoCounterpoise())
		   hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2) -= per_sca*damp*scafac2
		     *MonB.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);

		 //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 ){
		 //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		 //  printf("monB\n");
		 //  printf(" per_sca = %f\n",per_sca);
		 //  printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symi+xyz2);
		 //  printf(" damp = %f\n",damp);
		 //  printf(" Mon QM Hess = %f\n",MonB.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		 // printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		 //}

		 //MM Monomers
		 if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		   hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2) += per_sca*damp*scafac*
		     MonB.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
	    
		 //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558){
		 //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		 //  printf("monB\n");
		 //  printf(" per_sca = %f\n",per_sca);
		 //  printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		 //  printf(" damp = %f\n",damp);
		 //  printf(" Mon MM Hess = %f\n",MonB.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		 //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		 //}

	    
	    //}


	      }//end of xyz2
	    }//end of xyz1
	  }//end of loop over jatom
	}//end of loop over iatom

	//1st interaction block contribution
	for(int iatom=0;iatom<Na;iatom++){
	  	
	  //determine which atom iatom is symmetrical equivalent on MonA
          int Symi = (int) SymAtomList[iatom] - 1;
	  if(!FirstAInDim){
	    Symi-= Nb;
	  }

	  for(int jatom=0;jatom<Nb;jatom++){
	    
	    //determine which atom jatom is symmetrical equivalent on MonB
            int Symj = (int) SymAtomList[Na+jatom] - 1;
	    if(FirstAInDim){
	      Symj-= Na;
	    }
	    
	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //if ( !Params::Parameters().TinkerDebug() ) {
		    hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) += per_sca*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*Dimers.GetQMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4)
                        -scafac*Dimers.GetMMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    

		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("1st\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Unrotated QM = %15.9f\n", Dimers.GetQMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", Dimers.GetMMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) );
		    //}	
		    
		    hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) += per_sca
		      *Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*Dimers.GetQMGradient()[3*Na+3*jatom+xyz4]-scafac*Dimers.GetMMGradient()[3*Na+3*jatom+xyz4]);
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //   printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("1st\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" damp Grad iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		    //  printf(" QM Dimer Grad jatom = %15.9f\n",Dimers.GetQMGradient()[3*Na+3*jatom+xyz4] );
		    //  printf(" MM Dimer Grad jatom = %15.9f\n",Dimers.GetMMGradient()[3*Na+3*jatom+xyz4] );
		    //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		    //}			    
		    
		    hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) += per_sca
		      *Rot(xyz2,xyz4)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*Dimers.GetQMGradient()[3*iatom+xyz3]-scafac*Dimers.GetMMGradient()[3*iatom+xyz3]);
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("1st\n");
		    //  printf(" Element = (%i,%i)\n",3*Symi+xyz3,3*Na+3*Symj+xyz4);
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" damp Grad jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]);
		    //  printf(" QM Dimer Grad Symi = %15.9f\n",Dimers.GetQMGradient()[3*iatom+xyz3]);
		    //  printf(" MM Dimer Grad Symi = %15.9f\n",Dimers.GetMMGradient()[3*iatom+xyz3]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		    //}	

		    hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) += per_sca*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*Na+3*jatom+xyz4)
			 *(scafac2*Dimers.GetQMIntEnergy()-scafac*Dimers.GetMMIntEnergy());
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    // printf("1st\n");
		    //   printf(" per_sca = %f\n",per_sca);
		    //   printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" QM Int = %f\n",Dimers.GetQMIntEnergy());
		    //  printf(" MM Int = %f\n",Dimers.GetMMIntEnergy());
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) );
		    //}

		    //if (!Params::Parameters().Do_fdTinkerHessian() ) {  
		    
		  
		
		  }//loop over xyz4
                
                  //QM Monomers:already subtracted if counterpoise corrected
	          if(!Params::Parameters().DoCounterpoise())
	            hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) -= per_sca*scafac2*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonB.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*scafac2*MonA.GetQMGradient()[3*Symi+xyz1]);
		  
		  //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //    printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //    printf("1st\n");
		   //   printf(" per_sca = %f\n",per_sca);
		   //   printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz2);
		   //   printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		   //   printf(" Rot = %f\n",Rot(xyz1,xyz3));
		   //   printf(" Rot = %f\n",Rot(xyz2,xyz3));
		   //   printf(" damp Grad iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		   //   printf(" QM Mon Grad Symj = %f\n",MonB.GetQMGradient()[3*Symj+xyz2] );
		   //   printf(" damp Grad jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		   //   printf(" QM Mon Grad Symi = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		   //   printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		   // }			    
		    
		    //MM Monomers:already subtracted if qchem is used and counterpoise corrected
		    if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3 )
		      hess.Element(startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) += per_sca*scafac2*
			(Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonB.GetMMGradient()[3*Symj+xyz2]
			 +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*scafac2*MonA.GetMMGradient()[3*Symi+xyz1]);
			
		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("1st\n");
		    //  printf(" per_sca = %f\n",per_sca);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Symj+xyz2);
		    //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		    //  printf(" damp Grad iatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		    //  printf(" QM Mon Grad Symj = %f\n",MonB.GetMMGradient()[3*Symj+xyz2] );
		    //  printf(" damp Grad jatom = %f\n",Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		    //  printf(" QM Mon Grad Symi = %f\n",MonB.GetMMGradient()[3*Symi+xyz1]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		    //}		

		}//loop over xyz3
	      }//loop over xyz2
	    }//loop over xyz1
	  }//loop over jatom
	}//loop over iatom 

	//2nd interaction block contribution ;
	for(int iatom=0;iatom<Nb;iatom++){
	  
	  //determine which atom iatom is symmetrical equivalent on MonB
          int Symi = (int) SymAtomList[Na+iatom] - 1;
	  if(FirstAInDim){
	    Symi-= Na;
	  }

  
	  for(int jatom=0;jatom<Na;jatom++){
	    
	    //determine which atom jatom is symmetrical equivalent on MonA
            int Symj = (int) SymAtomList[jatom] - 1;
	    if(!FirstAInDim){
	      Symj-= Nb;
	    }

	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //if ( !Params::Parameters().TinkerDebug() ) {
	      
		    //including dimer contribution
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += per_sca*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*Dimers.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4)
			-scafac*Dimers.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4));
	    
		    //if(startB+3*Symi+xyz1 == 1874 && startA_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" QM Dimer Hess = %f\n",Dimers.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" MM Dimer Hess = %f\n",Dimers.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}
		    
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += per_sca
		      *Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*Dimers.GetQMGradient()[3*jatom+xyz4]-scafac*Dimers.GetMMGradient()[3*jatom+xyz4]);
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startA_image+3*Symj+xyz2 == 558  && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element(%i,%i)\n",3*Na+3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Grad Damp iatom = %f\n", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		    //  printf(" QM Dimer Grad jatom = %f\n",Dimers.GetQMGradient()[3*jatom+xyz4]);
		    //  printf(" MM Dimer Grad jatom = %f\n",Dimers.GetMMGradient()[3*jatom+xyz4]);
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}
		    
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += per_sca
		      *Rot(xyz2,xyz4)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*Dimers.GetQMGradient()[3*Na+3*iatom+xyz3]-scafac*Dimers.GetMMGradient()[3*Na+3*iatom+xyz3]);
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startA_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element(%i,%i)\n",3*jatom+xyz4,3*Na+3*iatom+xyz3);
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		     // printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		     // printf(" Grad Damp Symj = %f\n", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
		    //  printf(" QM Dimer Grad Symi = %f\n",Dimers.GetQMGradient()[3*Na+3*iatom+xyz3]);
		    //  printf(" MM Dimer Grad Symi = %f\n",Dimers.GetMMGradient()[3*Na+3*iatom+xyz3]);
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}
		    
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += per_sca*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*jatom+xyz4)
		      *(scafac2*Dimers.GetQMIntEnergy()-scafac*Dimers.GetMMIntEnergy());
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startA_image+3*jatom+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element(%i,%i)\n",3*Na+3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",Dimers.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*jatom+xyz4) );
		    //  printf(" QM Int = %f\n",Dimers.GetQMIntEnergy());
		    //  printf(" MM Int = %f\n",Dimers.GetMMIntEnergy());
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*jatom+xyz2));
		    //}
		    
		    //}
		    //if (!Params::Parameters().Do_fdTinkerHessian() ) {
		    
		  
		  }//loop over xyz4

                  //Rotating gradient of monomers contribution
		    
		  //QM Monomers:already subtracted if  counterpoise corrected
                  if(!Params::Parameters().DoCounterpoise())
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) -= per_sca*scafac2*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonA.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonB.GetQMGradient()[3*Symi+xyz1]); 
		  
		  //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*jatom+xyz2 == 558  && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //  printf("2nd\n");
		  //  printf(" element(%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" element(%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		   // printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		   // printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		   // printf(" Grad Damp iatom = %f\n", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		   // printf(" QM Monomer Grad Symj = %f\n",MonA.GetQMGradient()[3*Symj+xyz2]);
		   // printf(" Grad Damp jatom = %f\n", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		   // printf(" QM Monomer Grad Symi = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		   // printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		  //}	
		    
	          //MM Monomers:already subtracted if qchem is used and counterpoise corrected
	          if(!Params::Parameters().DoCounterpoise()  || !Params::Parameters().GetMMType()!=3)
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += per_sca*scafac*
		      (Rot(xyz1,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonA.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonB.GetMMGradient()[3*Symi+xyz1]); 
		  
		  //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*jatom+xyz2 == 558  && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //    printf("using d(%i,%i)\n",Dimers.GetIndexA(),Dimers.GetIndexB());
		  //     printf("2nd\n");
		  //    printf(" element(%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //    printf(" element(%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //    printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //    printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //    printf(" Grad Damp iatom = %f", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //    printf(" MM Monomer Grad Symj = %f\n",MonA.GetMMGradient()[3*Symj+xyz2]);
		  //    printf(" Grad Damp jatom = %f", Dimers.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //    printf(" MM Monomer Grad Symi = %f\n",MonB.GetMMGradient()[3*Symi+xyz1]);
		  //    printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		  // }


		}//loop over xyz3
	      }//loop over xyz2
	    }//loop over xyz1
	  }//loop over jatom
	}//loop over iatom
      }//loop over entries in Symmetry List
    }//if statement for the symmetry factor

  }//loop over dimers
  //printf("Finished two body contributions\n");
  //fflush(stdout);

  //Cluster::cluster().PrintHessian("HMBI supercell Hessian before periodic",hess);

  for (int i=1;i<=NDim_images;i++) {
    //    for (int i=1;i<=1;i++) {
    Dimer DimerImages( Cluster::cluster().GetDimerImage(i) );
    //int indexA = DimerImages.GetIndexA();
    //int indexB = DimerImages.GetIndexB();
    //int indexB_ref = DimerImages.GetReferenceMonomerIndex();

    int index3 = DimerImages.GetIndexA();        
    int index4 = DimerImages.GetReferenceMonomerIndex();   
    int indexB = DimerImages.GetIndexB();

    //Matrix Rot3 = Cluster::cluster().GetMonomer(index3).GetRotationMatrix();
    //Matrix Rot4 = Cluster::cluster().GetMonomer(index4).GetRotationMatrix();

    //Monomer MonA = DimerImages.GetMonomerA();//Cluster::cluster().GetMonomer(indexA);
    //Monomer MonB = DimerImages.GetMonomerB();//Cluster::cluster().GetMonomer(indexB_ref);
    Monomer Mon3 = Cluster::cluster().GetMonomer(index3);
    Monomer Mon4 = Cluster::cluster().GetMonomer(index4);
    int Na = DimerImages.GetMonomerA().GetNumberOfAtoms();
    int Nb = DimerImages.GetMonomerB().GetNumberOfAtoms();
    //int startA = key[indexA];     
    //int startB = key[indexB_ref];
    
    //MonA.PrintAll();
    //MonB.PrintAll();
    //DimerImages.PrintAll();
    
    //int *cell_index = NULL;
    //cell_index = DimerImages.GetImageCell();
    
    //damping functation
    double c1=0, c0=0;
    c0 = Params::Parameters().GetLocalCutoff(0);
    c1 = Params::Parameters().GetLocalCutoff(1);
    double  damp = DimerImages.GetDampingFactor(c0,c1);

    //int uci = GetCellIndex(na,nb,nc,cell_index[0],cell_index[1],cell_index[2]); //uci = unitcell index
    //int uci_mirror = GetCellIndex(na,nb,nc,-1*cell_index[0],-1*cell_index[1],-1*cell_index[2]); 
    vector<int> SymmetryList = DimerImages.GetSymmetryList();
    if(damp >= pow(10.0,-6)){
      for(int j=0; j<SymmetryList.size()/2;j++){
    
	//2*j is the first molecule in the dimer. 2*j+1 is the second monomer in the dimer
	int index1 = SymmetryList[2*j];
	int index2 = SymmetryList[2*j+1];
	Monomer MonA;
	Monomer MonB;	
  
	int cell_index[3];
	cell_index[0] = DimerImages.GetSymmetricalImageCell()[3*j];
	cell_index[1] = DimerImages.GetSymmetricalImageCell()[3*j+1];
	cell_index[2] = DimerImages.GetSymmetricalImageCell()[3*j+2];	
	int uci = GetCellIndex(na,nb,nc,cell_index[0],cell_index[1],cell_index[2]); //uci = unitcell index
	int uci_mirror = GetCellIndex(na,nb,nc,-1*cell_index[0],-1*cell_index[1],-1*cell_index[2]);
	//uci = unitcell index for unitcell from inversion
	// uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
      	
        //Rotation under symmetry H12 = R13*H34*R24' where R13=R24 and R13=R1'*R3
	Matrix Rot = DimerImages.GetRotationList()[j];
        Rot.Transpose();
        Vector SymAtomList = DimerImages.GetAtomEquivalency()[j];

        //Is MonA symmetrical to the first monomer in the dimer(Mon3)
        bool FirstAInDim = false;

	//Determines if Mon3 is symmetrical to Mon1 or Mon2 by determining
	//Need to know if Mon3 and Mon4 have different number to atoms
        for(int atom = 0; atom < Mon3.GetNumberOfAtoms();atom++){
	  if(SymAtomList[atom] == 1) FirstAInDim = true;
	}

	int startA;
	int startB;
	int startA_image;
	int startB_image;

	//MonA is symmetrical to Mon3
	if(FirstAInDim){

          MonA = Cluster::cluster().GetMonomer(index1);
	  MonB = Cluster::cluster().GetMonomer(index2);	

	  //uci = unitcell index for unitcell from inversion
	  // uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
	  startA = row_key[primitive_cell_index*NMon+index1-1];
          startB = row_key[primitive_cell_index*NMon+index2-1];
          startA_image = row_key[uci_mirror*NMon+index1-1];	  
          startB_image = row_key[uci*NMon+index2-1];
          


	  }
	  //MonB is symmetrical to Mon3
	  else{

	    MonA = Cluster::cluster().GetMonomer(index2);
	    MonB = Cluster::cluster().GetMonomer(index1);		
	    
	    //uci = unitcell index for unitcell from inversion
	    // uci_mirror would not probably work when we have supercells like 2x2x3 or 2x2x2 but the current logic works for 3x3x3, 1x3x5, etc.
	    startA = row_key[primitive_cell_index*NMon+index2-1];
	    startB = row_key[primitive_cell_index*NMon+index1-1];	  
	    startA_image = row_key[uci*NMon+index2-1];	  
	    startB_image = row_key[uci_mirror*NMon+index1-1];
	  }

	//printf("d(%i,%i) k = %i %i %i\n",index1,index2,cell_index[0],cell_index[1],cell_index[2]);

	//Piece arising from first monomer. 
	for(int iatom=0;iatom<Na;iatom++){
	  
	  //determine which atom iatom is symmetrical equivalent on MonA
	  int Symi = (int) SymAtomList[iatom] - 1;
	  if(!FirstAInDim){
	    Symi-= Nb;
	  }

	  for(int jatom=0;jatom<Na;jatom++){
	      
	    //determine which atom jatom is symmetrical equivalent on MonA
	    int Symj = (int) SymAtomList[jatom] - 1;
	    if(!FirstAInDim){
	      Symj-= Nb;
	    }
  
	    double c1=0, c0=0;
	    c0 = Params::Parameters().GetLocalCutoff(0);
	    c1 = Params::Parameters().GetLocalCutoff(1);
	    double  damp = DimerImages.GetDampingFactor(c0,c1);

	    // piece arising from first monomer
	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //if ( !Params::Parameters().TinkerDebug() ) {	    
		    
		    hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) += 0.5*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*DimerImages.GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4) 
			-scafac*DimerImages.GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));	
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		     // printf(" Unrotated QM = %15.9f\n", DimerImages.GetQMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", DimerImages.GetMMHessian().Element(3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot42 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 )  );
		    //}	
	

		    //}
		    //if (!Params::Parameters().Do_fdTinkerHessian() ) {
		    
		    hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 )  += 0.5
		      *Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*DimerImages.GetQMGradient()[3*jatom+xyz4]-scafac*DimerImages.GetMMGradient()[3*jatom+xyz4]);
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		    //  printf(" Grad QM Dimer jatom = %f\n",DimerImages.GetQMGradient()[3*jatom+xyz4]);
		    //  printf(" Grad MM Dimer jatom = %f\n",DimerImages.GetMMGradient()[3*jatom+xyz4]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}
		     

		    hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) +=  0.5
		      *Rot(xyz2,xyz4)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*DimerImages.GetQMGradient()[3*iatom+xyz3]-scafac*DimerImages.GetMMGradient()[3*iatom+xyz3]);

		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz1,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monA\n");
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Grad Damp jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
		    //  printf(" Grad QM Dimer iatom = %f\n",DimerImages.GetQMGradient()[3*iatom+xyz3]);
		    //  printf(" Grad MM Dimer iatom = %f\n",DimerImages.GetMMGradient()[3*iatom+xyz3]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}	 

		    hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) += 0.5*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4)
		      *(scafac2*DimerImages.GetQMIntEnergy()-scafac*DimerImages.GetMMIntEnergy());

		    //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monA\n");
		    //   printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot13(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot24(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" Int QM Dimer = %f\n",DimerImages.GetQMIntEnergy());
		    //  printf(" Int MM Dimer = %f\n",DimerImages.GetMMIntEnergy());
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		    //}	   	 

		    

		    //}
		  }//loop over xyz4

                  //Rotating gradient monomers contribation

		  //QM:already subtracted if  counterpoise corrected
                  if(!Params::Parameters().DoCounterpoise()) 
	            hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) -= 0.5*scafac2*
	              (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonA.GetQMGradient()[3*Symj+xyz2]
		      +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonA.GetQMGradient()[3*Symi+xyz1]);

		  //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("monA\n");
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer Symj = %f\n",MonA.GetQMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" Grad QM Monomer jatom = %f\n",MonA.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		  //}
	
		  //MM:already subtracted if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise()|| !Params::Parameters().GetMMType()!=3)
		    hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) += 0.5*scafac*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonA.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonA.GetMMGradient()[3*Symi+xyz1]);
		  
		  //if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 && (Rot(xyz1,xyz3) > 0.0001 || Rot(xyz2,xyz3) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("monA\n");
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" Grad MM Monomer Symj = %f\n",MonA.GetMMGradient()[3*Symj+xyz1]);
		  //  printf(" Grad Damp jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" Grad MM Monomer Symi = %f\n",MonA.GetMMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		  //}	 

		}//loop over xyz3

		//Including Monomer contribution
		//Don't have to be rotated so outside xyz3 and xyz4 loop
		
		//QM Monomers
		if(!Params::Parameters().DoCounterpoise())
		  hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) -= 0.5*damp*scafac2
		    *MonA.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);

		//if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 ){
		//  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		//  printf("monA\n");
		//  printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		//  printf(" damp = %f\n",damp);
		//  printf(" Mon QM Hess = %f\n",MonA.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		//  printf(" Current Hess Sum = %15.9f\n\n",hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2)  );
		//}

		  //MM Monomers   
		if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3)
		  hess.Element( startA+3*Symi+xyz1, startA+3*Symj+xyz2 ) +=  0.5*damp*scafac
		    *MonA.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);

		//if(startA+3*Symi+xyz1 == 1874 && startA+3*Symj+xyz2 == 558 ){
		//    printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		//    printf("monA\n");
		//    printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		//    printf(" damp = %f\n",damp);
		//    printf(" Mon MM Hess = %f\n",MonA.GetMMHessian().Element(3*iatom+xyz1,3*jatom+xyz2));
		//    printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startA+3*Symi+xyz1,startA+3*Symj+xyz2) );
		//  }
		    
	      }//loop over xyz2
	    }//loop over xyz1
	  }//loop over jatom
	}//loop over iatom 

	// piece arising from second monomer
	for(int iatom=0;iatom<Nb;iatom++){
	  	  
	  //determine which atom iatom is symmetrical equivalent on MonB
	  int Symi = (int) SymAtomList[Na+iatom] - 1;
	  if(FirstAInDim){
	    Symi-= Na;
	  }

	  for(int jatom=0;jatom<Nb;jatom++){

	    //determine which atom jatom is symmetrical equivalent on MonB
	    int Symj = (int) SymAtomList[Na+jatom] - 1;
	    if(FirstAInDim){
	      Symj-= Na;
	    }

	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //if ( !Params::Parameters().TinkerDebug() ) {
		    
		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2) += 0.5*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*DimerImages.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4)
                        -scafac*DimerImages.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
	            //  printf("monB\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Unrotated QM = %15.9f\n", DimerImages.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", DimerImages.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		    // printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startB+3*Symi+xyz1, startB+3*Symj+xyz2)  );
		    //}	

		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2 ) += 0.5
		      *Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*DimerImages.GetQMGradient()[3*Na+3*jatom+xyz4]-scafac*DimerImages.GetMMGradient()[3*Na+3*jatom+xyz4]);
		    
		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		    //  printf(" Grad QM Dimer jatom = %f\n",DimerImages.GetQMGradient()[3*Na+3*jatom+xyz4]);
		    //  printf(" Grad MM Dimer jatom = %f\n",DimerImages.GetMMGradient()[3*Na+3*jatom+xyz4]);
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    //}

		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2 ) += 0.5
		      *Rot(xyz2,xyz4)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*DimerImages.GetQMGradient()[3*Na+3*iatom+xyz3]-scafac*DimerImages.GetMMGradient()[3*Na+3*iatom+xyz3]);

		    //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Grad Damp jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]);
		    //  printf(" Grad QM Dimer iatom = %f\n",DimerImages.GetQMGradient()[3*Na+3*iatom+xyz3]);
		    //  printf(" Grad MM Dimer iatom = %f\n",DimerImages.GetMMGradient()[3*Na+3*iatom+xyz3]);
		    // printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    //}	 
		    
		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2 ) += 0.5*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4)
		      *(DimerImages.GetQMIntEnergy()-DimerImages.GetMMIntEnergy()); 
		    
		    //if(startA+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("monB\n");
		    //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Rot13(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot24(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Int QM  = %f\n",DimerImages.GetQMIntEnergy());
		    //  printf(" Int MM = %f\n",DimerImages.GetMMIntEnergy());
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		    //}
		
	
		    //}
		  }//loop over xyz4
		  
		  //Rotating gradient of monomers contributation

		  //QM:already subtracted if counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise())
		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2 ) -= 0.5*scafac2*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonB.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonB.GetQMGradient()[3*Symi+xyz1]);
		  
		  //if(startB+3*Symi+xyz1 == 1968 && startB+3*Symj+xyz2 == 1098 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("monB\n");
		  //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer Symj = %f\n",MonB.GetQMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //  printf(" Grad QM Monomer Symi = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		  //}

		    //MM:already subtracted if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise() || !Params::Parameters().GetMMType()!=3)
		    hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2 ) += 0.5*scafac*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonB.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonB.GetMMGradient()[3*Symi+xyz1]);
		  
		  //if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //   printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("monB\n");
		  //  printf(" Element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot2(xyz2,xyz4) = %f\n",Rot(xyz2,xyz3));		      
		  //  printf(" Grad Damp iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" Grad QM Monomer Symj = %f\n",MonB.GetMMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //  printf(" Grad MM Monomer iatom = %f\n",MonB.GetMMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startB+3*Symi+xyz1,startB+3*Symj+xyz2) );
		  //}

		}//loop over xyz3
		
		//Including Monomer contributation
	        //Don't have to be rotated so outside xyz3 and xyz4 loop

	        //QM Monomers
		if(!Params::Parameters().DoCounterpoise())
		  hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2) -= 0.5*damp*scafac2
		    *MonB.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		
		//if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 ){
		//   printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		//   printf("monB\n");
		//   printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		//   printf(" damp = %f\n",damp);
		//   printf(" Mon QM Hess = %f\n",MonB.GetQMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		//   printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startB+3*Symi+xyz1, startB+3*Symj+xyz2) );
		//}
     
		 //MM Monomers
		if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		  hess.Element( startB+3*Symi+xyz1, startB+3*Symj+xyz2) += 0.5*damp*scafac
		    *MonB.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2);
		
		//if(startB+3*Symi+xyz1 == 1874 && startB+3*Symj+xyz2 == 558 ){
		//  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		//  printf("monB\n");
		//  printf(" Element = (%i,%i)\n",3*Symi+xyz1,3*Symj+xyz2);
		//  printf(" damp = %f\n",damp);
		//  printf(" Mon MM Hess = %f\n",MonB.GetMMHessian().Element(3*Symi+xyz1,3*Symj+xyz2));
		//  printf(" Current Hess Sum = %15.9f\n\n", hess.Element(startB+3*Symi+xyz1, startB+3*Symj+xyz2) );
		//}
		

	      }//loop over xyz2
	    }//loop over xyz1
	  }//loop over jatom
	}//loop over iatom
	
	// there are also 2 additional pieces arising from the interaction between the 2 monomers (see diagram below)
	
	// [.....monA contribution........|...1st interaction block....]
	//[.....2nd interaction block....|....monB contribution.......]

	//1st interaction block contribution
	for(int iatom=0;iatom<Na;iatom++){

	  //determine which atom iatom is symmetrical equivalent on MonA
	  int Symi = (int) SymAtomList[iatom] - 1;
	  if(!FirstAInDim){
	    Symi-= Nb;
	  }	  

	  for(int jatom=0;jatom<Nb;jatom++){

	    //determine which atom jatom is symmetrical equivalent on MonB
	    int Symj = (int) SymAtomList[Na+jatom] - 1;
	    if(FirstAInDim){
	      Symj-= Na;
	    }

	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		//xyz3 and xyz4 indices are need to rotate elements of dimers hessian under symmetry
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    //including dimer contribution
		    //if ( !Params::Parameters().TinkerDebug() ) {
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) += 0.5*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*DimerImages.GetQMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4)-scafac*DimerImages.GetMMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    

		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //   printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf(" 1st\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Unrotated QM = %15.9f\n", DimerImages.GetQMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", DimerImages.GetMMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) );
		    //}

		    //damping factor contributation to the Hessian
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) +=  0.5
		      *Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*DimerImages.GetQMGradient()[3*Na+3*jatom+xyz4]-scafac*DimerImages.GetMMGradient()[3*Na+3*jatom+xyz4]);

		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf(" 1st\n");
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Unrotated QM = %15.9f\n", DimerImages.GetQMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Unrotated MM = %15.9f\n", DimerImages.GetMMHessian().Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Element of Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2) );
		    //}
		    
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) += 0.5
		      *Rot(xyz2,xyz4)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*DimerImages.GetQMGradient()[3*iatom+xyz3]-scafac*DimerImages.GetMMGradient()[3*iatom+xyz3]);

		    //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf(" 1st\n");
		    //  printf(" Element = (%i,%i)\n",3*jatom+xyz4,3*Na+3*iatom+xyz3);
		    //  printf(" damp Grad jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz4]);
		    //  printf(" QM Dimer Grad iatom = %15.9f\n",DimerImages.GetQMGradient()[3*Na+3*iatom+xyz3] );
		    //  printf(" MM Dimer Grad iatom = %15.9f\n",DimerImages.GetMMGradient()[3*Na+3*iatom+xyz3] );
		    //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    // printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		    //}	   
		    
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) += 0.5*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*Na+3*jatom+xyz4)
		      *(DimerImages.GetQMIntEnergy()-DimerImages.GetMMIntEnergy());
		    
		    // if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001 ){
		    //   printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf(" 1st\n");
		    //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*jatom+xyz4);
		    //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess   = %f\n",DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*iatom+xyz3,3*Na+3*jatom+xyz4));
		    //  printf(" Int QM = %15.9f\n",DimerImages.GetQMIntEnergy());
		    // printf(" Int MM = %15.9f\n",DimerImages.GetMMIntEnergy());
		    //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		    //}	

		    
		    
		  }//end of loop over xyz4

		  //Rotating gradient of monomers contribution
		    
		  //QM:already subtracted if counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise())
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) -= 0.5*scafac2*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonB.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonA.GetQMGradient()[3*Symi+xyz1]);
		    
		  //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf(" 1st\n");
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" damp Grad iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" QM Mon Grad Symj = %f\n",MonB.GetQMGradient()[3*Symj+xyz2] );
		  //  printf(" damp Grad jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //  printf(" QM Mon Grad jatom = %f\n",MonA.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		  //}		

		  //MM:already subtracted if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		    hess.Element( startA+3*Symi+xyz1, startB_image+3*Symj+xyz2) += 0.5*scafac*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]*MonB.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]*MonA.GetMMGradient()[3*Symi+xyz1]);
		    //}

		  //if(startA+3*Symi+xyz1 == 1874 && startB_image+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001)){
		  // printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		   // printf("1st\n");
		  //  printf(" Element = (%i,%i)\n",3*iatom+xyz3,3*Na+3*Symj+xyz2);
		  //  printf(" Element = (%i,%i)\n",3*Na+3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Element of Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Element of Rot(xyz2,xyz4) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" damp Grad iatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*iatom+xyz3]);
		  //  printf(" QM Mon Grad Symj = %f\n",MonB.GetQMGradient()[3*Symj+xyz3] );
		  //  printf(" damp Grad jatom = %f\n",DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*jatom+xyz3]);
		  //  printf(" QM Mon Grad Symi = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Current Hess Sum = %15.9f\n\n", hess.Element( startA+3*Symi+xyz1,startB_image+3*Symj+xyz2)  );
		  //}
		  
		}//end of loop over xyz3
	      }//end of loop over xyz2
	    }//end of loop over xyz1
	  }//end of loop over jatom
	}//end of loop over iatom

	// 2th interaction block contribution
	for(int iatom=0;iatom<Nb;iatom++){
	  
	  //determine which atom iatom is symmetrical equivalent on MonB
	  int Symi = (int) SymAtomList[Na+iatom] - 1;
	  if(FirstAInDim){
	    Symi-= Na;
	  }	
	  
	  for(int jatom=0;jatom<Na;jatom++){
	   
	    //determine which atom jatom is symmetrical equivalent on MonA
	    int Symj = (int) SymAtomList[jatom] - 1;
	    if(!FirstAInDim){
	      Symj-= Nb;
	    }

	    for(int xyz1=0;xyz1<3;xyz1++){
	      for(int xyz2=0;xyz2<3;xyz2++){
		for(int xyz3=0;xyz3<3;xyz3++){
		  for(int xyz4=0;xyz4<3;xyz4++){
		    
		    //if ( !Params::Parameters().TinkerDebug() ) {
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) += 0.5*damp*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *(scafac2*DimerImages.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4)
			- scafac*DimerImages.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4)); 
		    
		    //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*Symj+xyz2 == 558  && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element ( %i,%i)\n",3*Na+3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" damp = %f\n",damp);
		    //  printf(" Rot = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot = %f\n",Rot(xyz2,xyz4));
		    //  printf(" QM Dimer Hess = %f\n",DimerImages.GetQMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" MM Dimer Hess = %f\n",DimerImages.GetMMHessian().Element(3*Na+3*iatom+xyz3,3*jatom+xyz4));
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}

		    //damping factor contribution to the Hessian
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) +=  0.5
		      *Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]
		      *Rot(xyz2,xyz4)*(scafac2*DimerImages.GetQMGradient()[3*jatom+xyz4]-scafac*DimerImages.GetMMGradient()[3*jatom+xyz4]);

		    //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*jatom+xyz2 == 558  && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element =( %i,%i)\n",3*Na+3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Grad Damp iatom = %f\n", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		    //  printf(" QM Dimer Grad jatom = %f\n",DimerImages.GetQMGradient()[3*jatom+xyz4]);
		    //  printf(" MM Dimer Grad jatom = %f\n",DimerImages.GetMMGradient()[3*jatom+xyz4]);
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}

		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) +=  0.5
		      *Rot(xyz2,xyz4)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]
		      *Rot(xyz1,xyz3)*(scafac2*DimerImages.GetQMGradient()[3*Na+3*iatom+xyz3]-scafac*DimerImages.GetMMGradient()[3*Na+3*iatom+xyz3]);

		    //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*Symj+xyz2 == 558 && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element =( %i,%i)\n",3*jatom+xyz4,3*iatom+xyz3);
		    //  printf(" Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Grad Damp jatom = %f\n", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz4]);
		    //  printf(" QM Dimer Grad iatom = %f\n",DimerImages.GetQMGradient()[3*Na+3*iatom+xyz3]);
		    //  printf(" MM Dimer Grad iatom = %f\n",DimerImages.GetMMGradient()[3*Na+3*iatom+xyz3]);
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}

		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) +=  0.5*Rot(xyz1,xyz3)*Rot(xyz2,xyz4)
		      *DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*jatom+xyz4)
		      *(DimerImages.GetQMIntEnergy()-DimerImages.GetMMIntEnergy());
		    
		    //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*Symj+xyz2 == 558  && fabs(Rot(xyz1,xyz3)) > 0.0001 && fabs(Rot(xyz2,xyz4)) > 0.0001){
		    //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		    //  printf("2nd\n");
		    //  printf(" element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*jatom+xyz4);
		    //  printf(" Rot13 = %f\n",Rot(xyz1,xyz3));
		    //  printf(" Rot24 = %f\n",Rot(xyz2,xyz4));
		    //  printf(" Damp Hess = %f\n",DimerImages.GetSpatialDampingFunctionHessian(c0,c1).Element(3*Na+3*iatom+xyz3,3*jatom+xyz4) );
		    //  printf(" QM Int = %f\n",DimerImages.GetQMIntEnergy());
		    //  printf(" MM Int = %f\n",DimerImages.GetMMIntEnergy());
		    //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		    //}

		    
		  }//end of loop over xyz4
		  //Rotating gradient of monomers contribution
		  
		  //QM:already subtracted if counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise())
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) -= 0.5*scafac2*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonA.GetQMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonB.GetQMGradient()[3*Symi+xyz1]);
		  
		  //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*Symj+xyz2 == 558 && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001) ){
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("2nd\n");
		  //  printf(" element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" element = (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f\n", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" QM Monomer Grad Symj = %f\n",MonA.GetQMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" QM Monomer Grad iatom = %f\n",MonB.GetQMGradient()[3*Symi+xyz1]);
		  //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		  //}
		  
		  //MM:already subtracted if qchem is used and counterpoise corrected
		  if(!Params::Parameters().DoCounterpoise() || Params::Parameters().GetMMType()!=3)
		    hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2) +=  0.5*scafac*
		      (Rot(xyz1,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]*MonA.GetMMGradient()[3*Symj+xyz2]
		       +Rot(xyz2,xyz3)*DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]*MonB.GetMMGradient()[3*Symi+xyz1]);
		  //}
		  
		  //if(startB+3*Symi+xyz1 == 1874  && startA_image+3*Symj+xyz2 == 558  && (fabs(Rot(xyz1,xyz3)) > 0.0001 || fabs(Rot(xyz2,xyz3)) > 0.0001)) {
		  //  printf("using d(%i,%i)\n",DimerImages.GetIndexA(),DimerImages.GetIndexB());
		  //  printf("2nd\n");
		  //   printf(" element = (%i,%i)\n",3*Na+3*iatom+xyz3,3*Symj+xyz2);
		  //  printf(" element = (%i,%i)\n",3*jatom+xyz3,3*Symi+xyz1);
		  //  printf(" Rot(xyz1,xyz3) = %f\n",Rot(xyz1,xyz3));
		  //  printf(" Rot(xyz2,xyz3) = %f\n",Rot(xyz2,xyz3));
		  //  printf(" Grad Damp iatom = %f", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*Na+3*iatom+xyz3]);
		  //  printf(" MM Monomer Grad Symj = %f\n",MonA.GetMMGradient()[3*Symj+xyz2]);
		  //  printf(" Grad Damp jatom = %f\n", DimerImages.GetSpatialDampingFunctionGradient(c0,c1)[3*jatom+xyz3]);
		  //  printf(" MM Monomer Grad iatom = %f\n",MonB.GetMMGradient()[3*Symi+xyz1]);
		  //  printf(" Hess = %f\n\n",hess.Element(startB+3*Symi+xyz1, startA_image+3*Symj+xyz2));
		  //}
		  
		}//end of loop over xyz3
              }//end of loop over xyz2
	    }//end of loop over xyz1
	  }//end of loop over jatom
	}//end of loop over iatom
      }//loop over entries of symmetry list
    }//end of if statement for damp > 10^-6
  }//end of loop over imagine dimers
 }//end of if statement to see if this is a QE supercell hessian or not
  //Cluster::cluster().PrintHessian("HMBI supercell Hessian",hess);
  
  //delete [] key;

}          


int Supercell::GetCellIndex(int na, int nb, int nc, int ind_a, int ind_b, int ind_c) {

  int count_sc=0;

  int na_l, nb_l, nc_l, na_u, nb_u, nc_u;
  // these are the lower/upper limits for unitcell indices ... eg: -1 to 1 for a supercell of size 3x3x3, 0 to 1 for a supercell of size 2x2x2
  //now depending on whether na, nb and nc are odd/even, we have different values for these limits
  // if even
  if (na%2==0) {na_l = na/2-1; na_u=na/2;}
  if (nb%2==0) {nb_l = nb/2-1; nb_u=nb/2;}
  if (nc%2==0) {nc_l = nc/2-1; nc_u=nc/2;}
  // if odd
  if (na%2==1) {na_l = (na-1)/2; na_u = (na-1)/2;}
  if (nb%2==1) {nb_l = (nb-1)/2; nb_u = (nb-1)/2;}
  if (nc%2==1) {nc_l = (nc-1)/2; nc_u = (nc-1)/2;}

  // now if the image monomer of the dimer image is beyond the supercell, then we need to map such a monomer onto the appropriate monomer in the supercell
  // example: if image monomer is in (-2,1,1) cell and supercell size is 3x3x3, then it is mapped onto the monomer in (-2-(-2)/|-2|*3,1,1) = (1,1,1) cell 
  // so now spatial truncation scheme and supercell size are independent as should be the case (If not, you are wrong)

  int index_a = ind_a;
  int index_b =	ind_b;
  int index_c =	ind_c;

  while ( (index_a > na_u) || (index_a < -1*na_l) ) {
    index_a += -1*(ind_a/abs(ind_a)*na) ;
  }

  while ( (index_b > nb_u) || (index_b < -1*nb_l) ) {
    index_b += -1*(ind_b/abs(ind_b)*nb) ;
  }
  while ( (index_c > nc_u) || (index_c < -1*nc_l) ) {
    index_c += -1*(ind_c/abs(ind_c)*nc) ;
  }


  for (int i=-na_l;i<=na_u;i++) {
    for (int j=-nb_l;j<=nb_u;j++) {         
      for (int k=-nc_l;k<=nc_u;k++) {

        if ((i==index_a) && (j==index_b) && (k==index_c)) {

          goto out_of_loops;
        }
        count_sc++;  

      }
    }
  }

  out_of_loops:

  if (count_sc == na*nb*nc) { //something is wrong
    // this case should not arise anymore since we neglect such dimer images but lookout for it during debugging
    // this case arised when there was a dimerimage with a cellindex=(-2,-1,-1) in a 3x3x3 supercell job which allows indices from -1 to +1 only
    // this case was possible because the truncation cutoffs were large but the unitcell lengths were compartively smaller

    cout << "Inconsistent supercell size, the truncation cutoffs and the unitcell lengths. Bbye!!!\n";
    exit(1);
  }

  return count_sc;
}




void Supercell::ReadReciprocalSpaceBoundaryPointsAndCreateGrid(ifstream& infile) {

    string line;
  // Rewind the file, just in case
  Cluster::cluster().Rewind(infile);           

  // First find the number of such boundary points the user specified in the input file (1st line of this "$" section
  while ( !infile.eof() ) { 
    getline(infile,line);
    if (line.substr(0,24)=="$reciprocal_space_points") {
      getline(infile,line);
      istringstream iss(line);        
      iss >> tot_recip_boundary_pts;
      iss >> recip_stepsize;
      break;
    }
  }

  // now initialize a set of vectors based on this number
  recip_boundary_pts = new Vector[tot_recip_boundary_pts];                  
  plot_break = new bool[tot_recip_boundary_pts];
  for (int i=0;i<tot_recip_boundary_pts;i++)
    recip_boundary_pts[i].Initialize(3);

  // Rewind the file, just in case
  Cluster::cluster().Rewind(infile);

  // now read the reciprocal space boundary points in the user-specified order
  while ( !infile.eof() ) {       
    getline(infile,line);
    if (line.substr(0,24)=="$reciprocal_space_points") { 
      getline(infile,line); // this line is the total number of such points and stepsize... not desired anymore

      for (int i=0;i<tot_recip_boundary_pts;i++) {
        getline(infile,line);
        istringstream iss2(line);  
        iss2 >> recip_boundary_pts[i][0];
       	iss2 >>	recip_boundary_pts[i][1];
       	iss2 >>	recip_boundary_pts[i][2];
        iss2 >> plot_break[i];

      }

      break;
    }
  }

  // now create grid of k-space points for phonon dispersion curve (if desired)
  CreateReciprocalSpaceGrid();

}


void Supercell::CreateReciprocalSpaceGrid() {

  total_k_points = 1;
  for (int i=0;i<(tot_recip_boundary_pts-1);i++) {

    //first find the direction along which we move in the k-space using the current k-vec and the next k-vec
    Vector direction = recip_boundary_pts[i+1];
    direction -= recip_boundary_pts[i];

    //now find the length of this vector
    double length = direction.Norm();
    // now find the total number of points along this direction based on the stepsize and whether the band-plot requires a break 
    if (!plot_break[i]) {
      total_k_points += int(length / recip_stepsize);
    }
    else {
      total_k_points += 1;  // if we break we add just one point which is the next boundary point
    }

  }


  //store total_k_point to be used for other parts of the program
  Params::Parameters().SetSuperCellTotalKPoint(total_k_points);

  // now we need to initialize a total of "total_k_points" vectors in the k-space
  recip_grid_pts = new Vector[total_k_points];
  for (int i=0;i<total_k_points;i++) {
    recip_grid_pts[i].Initialize(3);
  }

  // 1st recip_grid_pts is the 1st zone boundary point specified
  recip_grid_pts[0] = recip_boundary_pts[0];
  int total_k_points_local=1;
  for (int i=0;i<(tot_recip_boundary_pts-1);i++) {
    int total_k_points_local_tmp=0;
    //first find the direction along which we move in the k-space using the current k-vec and the next k-vec
    Vector direction = recip_boundary_pts[i+1];
    direction -= recip_boundary_pts[i];
    direction.Print("direction");

    //now find the length of this vector
    double length = direction.Norm();

    double stepsize;
    if (!plot_break[i]) {
      total_k_points_local_tmp += int(length / recip_stepsize);
      stepsize = recip_stepsize;
    }
    else {
      total_k_points_local_tmp += 1;  // if we break we add just one point which is the next boundary point    
      stepsize = length;
    }
    // now normalize the direction vector
    direction.Normalize();
    direction.Print("direction after normalize");

    for (int j=total_k_points_local;j<(total_k_points_local+total_k_points_local_tmp);j++) {
      recip_grid_pts[j] = direction;
      recip_grid_pts[j].Scale(stepsize*(j-total_k_points_local+1));
      recip_grid_pts[j] += recip_boundary_pts[i];

    }
    total_k_points_local += total_k_points_local_tmp;
    
  }

}


void Supercell::CreateSupercellFractionalCoordinates() {
  double a = unit_cell[0].Norm();
  double b = unit_cell[1].Norm();
  double c = unit_cell[2].Norm();
  double alpha = acos( unit_cell[1].DotProduct(unit_cell[2]) / b / c);
  double beta = acos( unit_cell[0].DotProduct(unit_cell[2]) / a / c);   
  double gamma = acos( unit_cell[1].DotProduct(unit_cell[0]) / a / b);
  double sqrt_term = sqrt(1-cos(alpha)*cos(alpha)
                     -cos(beta)*cos(beta)
                     -cos(gamma)*cos(gamma)
                     + 2 * cos(alpha) * cos(beta) * cos(gamma));
  //form matrix to transform cartesian to fractional coordinates
  double mat[9];           
  mat[0] = 1/a;
  mat[1] = -cos(gamma)/a/sin(gamma);
  mat[2] = ( cos(alpha)*cos(gamma) - cos(beta) ) / a / sin(gamma) / sqrt_term;
  mat[3] = 0;
  mat[4] = 1/ b/ sin(gamma);
  mat[5] = (cos(beta)*cos(gamma)-cos(alpha)) / b / sin(gamma) / sqrt_term;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = sin(gamma)/ c/ sqrt_term;

  //cout << "Printing the Transformation Matrix that converts Cartesian to Fractional Coordinates. \n";
  //for (int i=0;i<3;i++) {
  //  for (int j=0;j<3;j++) {
  //    cout << setw(15) << mat[3*i+j] << "\t";
  //  }
  //  cout << "\n";
  //}

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  Supercell_FractionalCoordinates = new Vector[na*nb*nc]; // this way of defining is awesome

  //cout << "Natoms = " << Natoms << "\n";
  //fflush(stdout);

  for (int count_sc=0;count_sc<na*nb*nc;count_sc++) {
    Supercell_FractionalCoordinates[count_sc].Initialize(3*Natoms);
    for (int m=0;m<Natoms;m++) {


      Supercell_FractionalCoordinates[count_sc][3*m] += Supercell_AtomicCoordinates[count_sc][3*m]*mat[0];
      Supercell_FractionalCoordinates[count_sc][3*m] += Supercell_AtomicCoordinates[count_sc][3*m+1]*mat[1];
      Supercell_FractionalCoordinates[count_sc][3*m] += Supercell_AtomicCoordinates[count_sc][3*m+2]*mat[2];

      Supercell_FractionalCoordinates[count_sc][3*m+1] += Supercell_AtomicCoordinates[count_sc][3*m]*mat[3]; 
      Supercell_FractionalCoordinates[count_sc][3*m+1] += Supercell_AtomicCoordinates[count_sc][3*m+1]*mat[4];
      Supercell_FractionalCoordinates[count_sc][3*m+1] += Supercell_AtomicCoordinates[count_sc][3*m+2]*mat[5];

      Supercell_FractionalCoordinates[count_sc][3*m+2] += Supercell_AtomicCoordinates[count_sc][3*m]*mat[6];
      Supercell_FractionalCoordinates[count_sc][3*m+2] += Supercell_AtomicCoordinates[count_sc][3*m+1]*mat[7];
      Supercell_FractionalCoordinates[count_sc][3*m+2] += Supercell_AtomicCoordinates[count_sc][3*m+2]*mat[8];


    }
    //Supercell_FractionalCoordinates[count_sc].Print("Supercell_FractionalCoordinates");

  }

}


void Supercell::ComputePhonons() {

  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  double c1 = 5.89141e-7;
  double c2 = 3.94562;

  int control1 = 0, control2 = 0 ;

  Matrix prim_hess_mat = Cluster::cluster().GetHMBIHessian();
  //prim_hess_mat.Print("prim_hess_mat");
  //printf("\n");
  Matrix mw_prim_hess_mat = Cluster::cluster().ComputeMassWeightedHessian(prim_hess_mat);
  //mw_prim_hess_mat.Print("mw_prim_hess_mat"); 

  //fflush(stdout);
  Matrix Error_mat1 = mw_prim_hess_mat.BuildMatrixFromBlocks(mw_prim_hess_mat);
  //Error_mat1.Print("Error_mat1 init");
  //fflush(stdout);

  Matrix Error_mat2(6*Natoms,6*Natoms);

  phonon_mat.Initialize(3*Natoms,total_k_points);
  //for (int i=0;i<total_k_points;i++) {  
  //  recip_grid_pts[i].Print("recip_grid_pts[i]");    
  //}
  for (int i=0;i<total_k_points;i++) {// i is the k_point index
    
    Matrix dyna_mat(6*Natoms,6*Natoms);  //square matrix .. 6*Natoms to store the complex numbers
    Vector eigvals_tmp;
    Vector eigvals(3*Natoms);
    
    for (int f=0;f<Natoms;f++) {
      for (int g=0;g<Natoms;g++) {
	
        for (int h=0;h<na*nb*nc;h++) { //index of cells
	  
          Vector disp_vec(3);
          disp_vec[0] = Supercell_FractionalCoordinates[h][0]-Supercell_FractionalCoordinates[primitive_cell_index][0];  
          disp_vec[1] = Supercell_FractionalCoordinates[h][1]-Supercell_FractionalCoordinates[primitive_cell_index][1];
          disp_vec[2] = Supercell_FractionalCoordinates[h][2]-Supercell_FractionalCoordinates[primitive_cell_index][2];

	  
          double the_dot_product = disp_vec.DotProduct(recip_grid_pts[i]) ;
	  
          double mass_fac_f = 1.0, mass_fac_g = 1.0;
          double m1 = AtomicMasses[f];
          double m2 = AtomicMasses[g];
	  
          if (m1 == 0.0) {
            mass_fac_f = 0;
            m1 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting
          }
	  
          if (m2 == 0.0) {
            mass_fac_g = 0;
       	    m2 += 1.0 ;	// adding 1.0 to avoid 0/0 while mass-weighting
          }
	  
          for (int j=0;j<3;j++) { //index of coordinates
       	    for (int k=0;k<3;k++) { //index of coordinates
	      
              dyna_mat.Element(3*f+j,3*g+k) += mass_fac_g*mass_fac_f * 1/( sqrt(m1*m2) )
	                   	* 0.5*( hess.Element( 3*Natoms*primitive_cell_index + 3*f+j, h*3*Natoms + 3*g+k)
                                      +hess.Element( 3*Natoms*primitive_cell_index + 3*g+k, h*3*Natoms + 3*f+j) )
                               * cos(-2*pi*the_dot_product);

              dyna_mat.Element(3*f+3*Natoms+j,3*g+k) += mass_fac_g*mass_fac_f * 1/( sqrt(m1*m2) )
                               * 0.5*( hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
       	       	       	       	      +hess.Element( 3*Natoms*primitive_cell_index + 3*g+k,h*3*Natoms + 3*f + j) )
                               * sin(-2*pi*the_dot_product);

              dyna_mat.Element(3*f+j,3*g+3*Natoms+k) -= mass_fac_g*mass_fac_f * 1/( sqrt(m1*m2) )
                               * 0.5*( hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
       	       	       	       	      +hess.Element( 3*Natoms*primitive_cell_index + 3*g+k,h*3*Natoms + 3*f + j) )
                               * sin(-2*pi*the_dot_product);

              dyna_mat.Element(3*f+3*Natoms+j,3*g+3*Natoms+k) += mass_fac_g*mass_fac_f * 1/( sqrt(m1*m2) )
                               * 0.5*( hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
       	       	       	       	      +hess.Element( 3*Natoms*primitive_cell_index + 3*g+k,h*3*Natoms + 3*f + j) )
                               * cos(-2*pi*the_dot_product);

            }
          }
	  
        }
      }
    }
    //dyna_mat.Print("dyna_mat before non-analytical term");

    //if gamma point, store the dyna_mat into gamma_dyna_mat
    if (recip_grid_pts[i].DotProduct(recip_grid_pts[i]) == 0.0 ) {

      //cout << "this gamma-point calc. should appear only once\n";
      //fflush(stdout);
      Error_mat1 -= dyna_mat ;
      //Error_mat1.Print("Error_mat1 before scaling");
      
      //      Error_mat1.Scale(0.0);
      Error_mat1.Scale(1.0/(na*nb*nc));
      //Error_mat1.Print("Error_mat1");
      
      Matrix reduced_Error_mat2(mw_prim_hess_mat);// = Cluster::cluster().ProjectOutTranslationAndRotation(mw_prim_hess_mat);
      reduced_Error_mat2.Set();
      //reduced_Error_mat2.Print("reduced_Error_mat2 init");
      reduced_Error_mat2 -= mw_prim_hess_mat;
      //reduced_Error_mat2.Print("reduced_Error_mat2");
      Error_mat2 = reduced_Error_mat2.BuildMatrixFromBlocks(reduced_Error_mat2);

      //Error_mat2.Print("Error_mat2 before scaling");

      Error_mat2.Scale(1.0/(na*nb*nc));
      //Error_mat2.Print("Error_mat2");

    }

   //if(i==primitive_cell_index){
   //  printf( "priv cell = %i\n",primitive_cell_index);
   //  printf ("k_points = %i\n",i);
   //  dyna_mat.PrintHessian("Final Dynamical Matrix");
   // }

   /*if(i==14){
     printf ("k_points = %i\n",i);
     recip_grid_pts[i].Print("recip points");
     printf("\n");
     dyna_mat.PrintHessian("Final Dynamical Matrix");
     printf("\n");
     fflush(stdout);
     exit(0);
    }*/

    eigvals_tmp = dyna_mat.Diagonalize();
    
    //eigvals_tmp.Print("Eigvals");
    //printf("\n");

    UnMassWeightedEigenVectors(dyna_mat);

    int Nfreq = eigvals_tmp.GetLength()/2;

    for ( int k=0; k<dyna_mat.GetCols();k++) {
      double norm_col =  dyna_mat.GetColumnVector(k).Norm();
      for ( int j=0; j<dyna_mat.GetRows();j++) {
        dyna_mat.Element(j,k) /= norm_col;
      }
    }
    
    
    // now using the filter properly for ice crystal where eight odd numbered monomers are important and not all 16 monomers ...
    
    //cout << "Printing the intensities\n";

    Vector sum_sq_vec(6*Natoms);
    for ( int k=0; k<dyna_mat.GetCols();k++) {
      //      if (k%2 == 0 ) { //only worry about the even numbered columns due to artificial degeneracy in freqs
      Vector dynacol = dyna_mat.GetColumnVector(k);
      double sum1 = 0.0;
      for (int mon_i=1; mon_i<=Cluster::cluster().GetNumberOfMonomers(); mon_i++) {
        if (mon_i%2 == 1) {
          double mon_at = Cluster::cluster().GetMonomer(mon_i).GetNumberOfAtoms() ;
          for (int j=0; j<3*mon_at; j++) {
	    
            sum1 += dynacol[(int)((mon_i-1)*mon_at+j)] * dynacol[(int)((mon_i-1)*mon_at+j)] 
                 + dynacol[(int)(dyna_mat.GetRows()/2 + (mon_i-1)*mon_at+j)] * dynacol[(int)(dyna_mat.GetRows()/2 + (mon_i-1)*mon_at+j)] ;
	    
          }
        }
      }
      if (sum1 < 0.01) {
        sum_sq_vec[k] += 1.0;
      }
      else {
        sum_sq_vec[k] += 1.0;
	
      }
      //cout << "eigenvector " << k << " has intensity = " << sum1 << "\n";
      
    }

    
    //sum_sq_vec.Print("sum_sq_vec");
    
    // now eigvals_tmp has degenerate values ... since we diagonalize the complex matrix of 3N size with a 6N real matrix
    for (int l=0;l<6*Natoms;l++) {
      if (l%2 == 0) {
        if (eigvals_tmp[l] < 0.0) { 
          eigvals[l/2] = - sqrt( fabs(eigvals_tmp[l] * sum_sq_vec[l])/c1)*c2 ;
        }
        else {
          eigvals[l/2] = sqrt( fabs(eigvals_tmp[l] * sum_sq_vec[l])/c1)*c2 ;
        }
      }
    }
    
    //if(i==0){
    //  eigvals_tmp.Print("eigvals_tmp");
    //	eigvals.Print("eig");
    // dyna_mat.PrintHessian("Normal Mode Eigenvectors");
    //}
    
    //dyna_mat goes from a 6Nx6N to a 3Nx3N
    dyna_mat = GetOnlyRealMatrixComponent(dyna_mat,eigvals);

    //eigvals.Print("eigvals");
    //phonon_mat.SetColumnVector(eigvals,i);
    //dyna_mat.Print("Trying to print Dynamat");
    //printf("\n");
  fflush(stdout);
/*
time for printing out the molden-format normal modes (april 17th 2013 - KDN)

one thing that is different from gamma-point vibrational analysis is that the normal modes could have imaginary components,
even though we have real frequencies.
 so the question is "How do we visualize these complex normal modes"

for starters, we just worry about the real part of these complex normal modes (coming out of the Cosine)

since dyna_mat is already in the unmassweighted normal mode form, we simply pick the 1st block of this matrix.
*/

    Vector orig_eigvals( eigvals );
    eigvals.SortLargestFirst();

    //int Nfreq = eigvals.GetLength();
    //printf("Nfreq = %i\n",Nfreq);
    Vector Orig_Index( Nfreq );
    Orig_Index = Cluster::cluster().CompareSortedOrigVecs(orig_eigvals,eigvals);
    //Vector orig_freqs( Freqs );
    dyna_mat = Cluster::cluster().SortColumns(dyna_mat,Orig_Index);
    
    //Unlike previous versions of HMBI, the columns of dyna_mat are now re-orginized in order of corresponding freqs.
    // Still need an order list for PrintNormalModes
    Orig_Index.IncrementalEntries(true);
    // Count number of negative eigenvalues (imaginary freqs)
    int Nimag = 0;
    for (int bb=0; bb<Nfreq;bb++) {
      if (eigvals[bb] < 0.0) {
        Nimag++;
      }
    }

    //printf("%d imaginary frequencies found\n\n",Nimag);
    //fflush(stdout);
    int Nonzero = eigvals.Nonzero();

    //Saving reference modes for the quasiharmonic approximating
    if(Params::Parameters().DoQuasiHarmonic()){
      if(!Quasiharmonic::quasiharmonic().IsReferenceCellInitialized()){
	//Quasiharmonic::quasiharmonic().SetReferenceMode();
	PrintNormalModes(eigvals, Nfreq, dyna_mat, Nonzero, Nimag, Orig_Index, i); 
      }
      //Reordering freq so they are the same order as in the reference cell
      else {
	//printf("\n\nReference was set\n");
        Quasiharmonic::quasiharmonic().ReadMoldenFrequencies(i);
	Vector orig_freqs = eigvals;
	eigvals = Quasiharmonic::quasiharmonic().OrderFrequencies(eigvals,dyna_mat);
	Orig_Index = Cluster::cluster().CompareSortedOrigVecs(orig_freqs,eigvals);
	if(Params::Parameters().SaveQuasiHarmonicCalculations()){
	  string file_name = Quasiharmonic::quasiharmonic().GetHessianType();
	  PrintNormalModes(eigvals, Nfreq, dyna_mat, Nonzero, Nimag, Orig_Index,i,file_name);
	}
	//exit(0);
      }
    }

    phonon_mat.SetColumnVector(eigvals,i);


    //if(i==0){
    //  dyna_mat.PrintHessian("Real Eigenvectors");
    //}       
    
    //PrintNormalModes(eigvals, Nfreq, real_dynamat, Nonzero, Nimag, Orig_Index, i);
    
    if(!Params::Parameters().DoQuasiHarmonic())
      PrintNormalModes(eigvals, Nfreq, dyna_mat, Nonzero, Nimag, Orig_Index, i);    

    eigvals.Set();
    eigvals_tmp.Set();
    dyna_mat.Set();
  }
  

  //phonon_mat.Transpose();
  //cout << "Printing Phonons...\n";
  //for (int kk=0;kk<(3*Natoms+2);kk++) {
  //  cout << "0  ";
  //}
  //cout << "\n";

//  Cluster::cluster().PrintHessian("phonon_mat",phonon_mat);
//   phonon_mat.PrintHessian("phonon_mat");
//  phonon_mat.Print("phonon_mat");

  //for (int jjj=0;jjj<total_k_points;jjj++) {             
  // cout << jjj << "  ";
  // for (int kkk=0;kkk<3*Natoms;kkk++) {
  //   cout <<  phonon_mat.Element(jjj,kkk) << "  ";
  // }            
  // cout << "  0\n";
  //}

//  cout << total_k_points - 1 << "  ";
//  for (int kk=0;kk<(3*Natoms+1);kk++) {
//    cout << "0  ";
//  }
//  cout << "\n";
  //phonon_mat.Transpose();

  //phonon_mat.PrintHessian("phonon_mat");

  if(Params::Parameters().DoQuasiHarmonic() && !Quasiharmonic::quasiharmonic().IsReferenceCellInitialized())
    Quasiharmonic::quasiharmonic().SetReferenceMode();

  //Calculate Density of States Calculations
  ComputePDOS(phonon_mat);

}

//similar to cluster but in supercell
void Supercell::SetAtomicSymbols() {

  for (int i=0;i<Natoms;i++) {
    AtomicMasses[i] = Cluster::cluster().GetAtomicMasses(i);
  }

}


void Supercell::CreateReciprocalSpaceSamplingGrid(ifstream& infile) {

  int shrink1, shrink2, shrink3; //the number of uniformly spaced grid points along each reciprocal space unit vector
  
  //JLM If not specified then it should correspond to the Supercell grid size;
  shrink1 = int(Supercell_Size[0]);
  shrink2 = int(Supercell_Size[1]);
  shrink3 = int(Supercell_Size[2]);

  string line;
  // Rewind the file, just in case                
  Cluster::cluster().Rewind(infile);

  // First find the number of such boundary points the user specified in the input file (1st line of this "$" section
  while ( !infile.eof() ) {
    getline(infile,line);
    if (line.substr(0,24)=="$reciprocal_space_points" ) {
      getline(infile,line);
      istringstream iss(line);
      iss >> shrink1;
      iss >> shrink2;
      iss >> shrink3;
      break;
    }
  }

  //cout << "shrink = " << shrink1 << " " << shrink2 << " " << shrink3 << "\n";

  // now initialize a set of vectors based on this number
  total_k_points = shrink1*shrink2*shrink3;
  // cout << "total_k_points = " << total_k_points << "\n";

  //store total_k_point to be used for other parts of the program
  Params::Parameters().SetSuperCellTotalKPoint(total_k_points);
  

  // now we need to initialize a total of "total_k_points" vectors in the k-space
  recip_grid_pts = new Vector[total_k_points];
  for (int i=0;i<total_k_points;i++) {                      
    recip_grid_pts[i].Initialize(3);
  }

  int count_grid=0;

  printf("Recip grid points\n");

  //Now make the uniformly spaced grid according to Monkhorst-Pack scheme (Phys. Rev. B, Volume 13, No. 12, Page 5188, 1976)
  //As of 27th August 2012, we do not consider any symmetry while creating this simple grid
  for (int i=1;i<=shrink1;i++) {
    for (int j=1;j<=shrink2;j++) {
      for (int k=1;k<=shrink3;k++) {


        recip_grid_pts[count_grid][0] = (2*i-shrink1-1.0)/2/shrink1 ;
        recip_grid_pts[count_grid][1] = (2*j-shrink2-1.0)/2/shrink2 ;
        recip_grid_pts[count_grid][2] = (2*k-shrink3-1.0)/2/shrink3 ;

	printf("%f\t%f\t%f\n",recip_grid_pts[count_grid][0],recip_grid_pts[count_grid][1],recip_grid_pts[count_grid][2]);        
        count_grid++;
      }
    }
  }



}


void Supercell::ComputeThermalPropertiesUsingPhonons(ifstream& infile) {

  double Temp_min,Temp_max,Temp_step;
  string line;
  // Rewind the file, just in case
  Cluster::cluster().Rewind(infile);


  //default values
  Temp_min = 0.0;
  Temp_max = 298;
  Temp_step = 1.0;

  // First find the number of such boundary points the user specified in the input file (1st line of this "$" section
  while ( !infile.eof() ) {     
    getline(infile,line);               
    if (line.substr(0,19)=="$thermal_parameters" ) {
      getline(infile,line);
      istringstream iss(line);              
      iss >> Temp_min;
      iss >> Temp_max;
      iss >> Temp_step;
      break;
    }  
  }  

  double Temp_current = Temp_min;
  printf("\n------Thermodynamic Properaties-----------------------\n");
  cout << setw(16) << "T [K] " << setw(16) << "F [kJ/mol] " << setw(16) << "S [J/mol/K] " << setw(16) << "E [kJ/mol] " << setw(16) << "C_v [J/mol/K] \n";
  while (Temp_current <= Temp_max) {

    double Helmholtz_E = 0.0;
    double Vib_Entropy = 0.0;
    double Vib_Internal_E = 0.0;
    double Sp_Heat_Cap_V = 0.0;


    for (int i=0;i<(3*Natoms);i++) {
      for (int j=0;j<total_k_points;j++) {


        // for thermal properties, we consider only real frequencies and neglect imaginary frequencies (need to cite this ... )
        if ( (phonon_mat.Element(i,j) > 0.0) && !((i<3) && (j==171)) ) {

          //First calculate Helmholtz free energy "F" 
          //1st term
          Helmholtz_E += 0.5 * h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na/1000.0;
          //2nd term
     	  Helmholtz_E += R_gas_constant*Temp_current*log(1-exp(-1.0 * h_planck *c_light*100.0* phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) ) / 1000.0;

          //Second calculate the entropy "S"
          //1st term
          Vib_Entropy += -1.0*R_gas_constant*log(1-exp(-1 * h_planck * c_light*100 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) );
          //2nd term (check sign)
          if (Temp_current != 0.0) {
       	  Vib_Entropy += 1.0 / Temp_current * h_planck * c_light*100.0 * phonon_mat.Element(i,j) *Na
                     *exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current)
                     /( -1*exp(-1* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) + 1.0);
          }
          else {
            Vib_Entropy += 0.0;
          }

       	  //Third calculate the energy "E"          
       	  //1st term
       	  Vib_Internal_E += 0.5 * h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na/1000.0;
          //2nd term
          Vib_Internal_E += h_planck * c_light*100.0 * phonon_mat.Element(i,j) *Na / 1000.0
                     * exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current)
                     /( -1*exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) + 1.0);


          //Fourth calculate the specific heat capacity at constant volume "C_v"
          //Only term
       	  if (Temp_current != 0.0) { 
            Sp_Heat_Cap_V += R_gas_constant*Temp_current*Temp_current/Temp_current/Temp_current
                          * ( h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current)
                          * ( h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current)
                          * exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current)
                          / (-1* exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) + 1.0)
       	       	       	  / (-1* exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / Temp_current) + 1.0);
       	  }
       	  else {
       	    Sp_Heat_Cap_V += 0.0;
       	  }



        }

        else {
          Sp_Heat_Cap_V += 0.0;
          Helmholtz_E += 0.0;
          Vib_Entropy += 0.0;
          Vib_Internal_E += 0.0;
        }
      }
    }


    Helmholtz_E *= 1.0/total_k_points;
    Vib_Entropy *= 1.0/total_k_points;
    Vib_Internal_E *= 1.0/total_k_points;
    Sp_Heat_Cap_V *= 1.0/total_k_points;

    cout << setw(16) << Temp_current << setw(16) << Helmholtz_E << setw(16) << Vib_Entropy << setw(16) << Vib_Internal_E << setw(16) << Sp_Heat_Cap_V << "  final\n";
    fflush(stdout);
    Temp_current += Temp_step;
  }
  printf("-----------------------------------------------------\n");
}

//Get vibrational energy at a specific temperature
void Supercell::ComputeVibrationalEnergy(){

  double temp = Params::Parameters().GetTemperature();
  double inverse_temp = 1/(k_boltz*temp);
  double Helmholtz_E = 0.0;
  double Vib_Entropy = 0.0;
  double Cv = 0.0;

  //macroscopic gruneisen parameter, not the same gruneisen parameter used by quasiharmonic approximation
  //double macro_gruneisen = 0.0;


  for (int i=0;i<(3*Natoms);i++) {
    for (int j=0;j<total_k_points;j++) {
     
      // for thermal properties, we consider only real frequencies and neglect imaginary frequencies (need to cite this ... )
      if ( (phonon_mat.Element(i,j) > 0.0) && !((i<3) && (j==171)) ) {
	
	//First calculate Helmholtz free energy "F" 
	//1st term
	Helmholtz_E += 0.5 * h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na/1000.0;
	
	if (temp != 0.0) {
	  //2nd term
	  Helmholtz_E += R_gas_constant*temp*log(1-exp(-1.0 * h_planck *c_light*100.0* phonon_mat.Element(i,j) * Na / R_gas_constant / temp) ) / 1000.0;

          //Second calculate the entropy "S
	  //1st term
	  Vib_Entropy += -1.0*R_gas_constant*log(1-exp(-1 * h_planck * c_light*100 * phonon_mat.Element(i,j) * Na / R_gas_constant / temp) );
	  //2nd term (check sign)
       	  Vib_Entropy += 1.0 / temp * h_planck * c_light*100.0 * phonon_mat.Element(i,j) *Na
	    *exp(-1.0* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / temp)
	    /( -1*exp(-1* h_planck * c_light*100.0 * phonon_mat.Element(i,j) * Na / R_gas_constant / temp) + 1.0);
	}

	//Constant volume heat capacity
	double Cv_mode = pow((c_light*100.0 * phonon_mat.Element(i,j) *h_planck/temp/ (exp(phonon_mat.Element(i,j)*c_light*100*h_planck*inverse_temp ) - 1 )),2)
	  * Na/k_boltz * exp( phonon_mat.Element(i,j) *c_light*100*h_planck*inverse_temp);
	Cv += Cv_mode;

	//if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable()){
	//  double micro_gruneisen = Quasiharmonic::quasiharmonic().GetGruneisenParameters()[i];
	//  macro_gruneisen *= micro_gruneisen*Cv_mode;
	//}


       //Next calculate entropy
       //1st term
       //if(temp == 0.0)
	//entropy += h_planck*phonon_mat.Element(i,j)*c_light*100.0*Na/temp/(exp(phonon_mat.Element(i,j)*c_light*100*h_planck*inverse_temp)-1);

       //printf("term 1 freq = %f entropy = %f\n",phonon_mat.Element(i,j),entropy);

       //2nd term
       //entropy -= R_gas_constant*log(1-exp(-phonon_mat.Element(i,j)*c_light*100*h_planck*inverse_temp));
       
      //printf("term 2 freq = %f entropy = %f\n",phonon_mat.Element(i,j),entropy);
      //exit(0);
      }
    }
  }
  //if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
  // macro_gruneisen /= Cv;


  Helmholtz_E *= 1.0/double(total_k_points);
  Vib_Entropy *= 1.0/double(total_k_points);
  Cv *= 1.0/double(total_k_points);

  printf("Helmholtz Energy  (at %.0f K)  = %9.3f kJ/mol\n",
	   temp,Helmholtz_E);
  printf("Entropy (at %.0f K) = %9.3f kJ/mol\n",
           temp,Vib_Entropy); 

  Cluster::cluster().SetHelmholtzEnergy(Helmholtz_E/2625.5);
  Cluster::cluster().SetEntropy(Vib_Entropy/2625.5/1000);
  Cluster::cluster().SetCv(Cv);
  
  //if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().AreQHFAvailable())
  //  Cluster::cluster().SetGruneisenParameter(macro_gruneisen);
}

void Supercell::UnMassWeightedEigenVectors(Matrix& Hessian) {

  for (int i=0;i<Natoms;i++) { 
    double m1 = AtomicMasses[i];
    for (int j=0; j<6*Natoms;j++) {     
      for (int p=3*i;p<3*(i+1);p++) {
        Hessian.Element(p,j) *= 1/sqrt(m1);
        Hessian.Element(3*Natoms+p,j) *= 1/sqrt(m1);
      }
    }
  }
//  return Hessian;
}


void Supercell::Compute_Elastic_Square_Brackets() {

//  constant int const_natoms = Natoms;
  double C_one[3][3][3][16][16] = {0.0};   
  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);  

  for (int f=0;f<Natoms;f++) {          
    for (int g=0;g<Natoms;g++) {
      for (int h=0;h<na*nb*nc;h++) {      

        Vector disp_vec(3);
        disp_vec[0] = Supercell_FractionalCoordinates[h][3*g+0]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+0]; 
        disp_vec[1] = Supercell_FractionalCoordinates[h][3*g+1]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+1];
        disp_vec[2] = Supercell_FractionalCoordinates[h][3*g+2]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+2];

        double mass_fac_f = 1.0, mass_fac_g = 1.0;

        double m1 = AtomicMasses[f];
        double m2 = AtomicMasses[g];          

        if (m1 == 0.0) {
          mass_fac_f = 0;
          m1 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting
        }     
        if (m2 == 0.0) {
          mass_fac_g = 0;                
          m2 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting
        }


        for (int j=0;j<3;j++) {
          for (int k=0;k<3;k++) {
            for (int l=0;l<3;l++) {

              C_one[j][k][l][f][g] -= mass_fac_f*mass_fac_g*2*pi/( sqrt(m1*m2) )
                               * hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
                               * disp_vec[l];

            }
          }
	}
      }
    }
  }

  double C_two[3][3][3][3][16][16] = {0.0};
//  int na = Supercell_Size[0];
//  int nb = Supercell_Size[1]; 
//  int nc = Supercell_Size[2];   

  for (int f=0;f<Natoms;f++) {
    for (int g=0;g<Natoms;g++) {
      for (int h=0;h<na*nb*nc;h++) {

        Vector disp_vec(3);
        disp_vec[0] = Supercell_FractionalCoordinates[h][3*g+0]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+0];
        disp_vec[1] = Supercell_FractionalCoordinates[h][3*g+1]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+1];
        disp_vec[2] = Supercell_FractionalCoordinates[h][3*g+2]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+2];

        double mass_fac_f = 1.0, mass_fac_g = 1.0;  

        double m1 = AtomicMasses[f];  
        double m2 = AtomicMasses[g];               

        if (m1 == 0.0) {               
          mass_fac_f = 0;               
          m1 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }                    
        if (m2 == 0.0) {  
          mass_fac_g = 0;     
          m2 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }  

        for (int j=0;j<3;j++) {
          for (int k=0;k<3;k++) {
            for (int l=0;l<3;l++) {
              for (int m=0;m<3;m++) {


                C_two[j][k][l][m][f][g] -= mass_fac_f*mass_fac_g*4*pi*pi/( sqrt(m1*m2) )
                               * hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
                               * disp_vec[l] * disp_vec[m] ;

              }
            }
          }
	}
      }
      if ((f==0) && (g==0)) {                
      cout << "\n C_two \n";
        cout << C_two[0][0][0][0][f][g] << "  " << C_two[0][0][0][1][f][g] << "  " << C_two[0][0][0][2][f][g] << "  "
             << C_two[0][0][1][0][f][g] << "  " << C_two[0][0][1][1][f][g] << "  " << C_two[0][0][1][2][f][g] << "  "
	     << C_two[0][0][2][0][f][g] << "  " << C_two[0][0][2][1][f][g] << "  " << C_two[0][0][2][2][f][g] << "\n";

        cout << C_two[0][1][0][0][f][g] << "  " << C_two[0][1][0][1][f][g] << "  " << C_two[0][1][0][2][f][g] << "  "
             << C_two[0][1][1][0][f][g] << "  " << C_two[0][1][1][1][f][g] << "  " << C_two[0][1][1][2][f][g] << "  "
             << C_two[0][1][2][0][f][g] << "  " << C_two[0][1][2][1][f][g] << "  " << C_two[0][1][2][2][f][g] << "\n";

        cout << C_two[0][2][0][0][f][g] << "  " << C_two[0][2][0][1][f][g] << "  " << C_two[0][2][0][2][f][g] << "  "
             << C_two[0][2][1][0][f][g] << "  " << C_two[0][2][1][1][f][g] << "  " << C_two[0][2][1][2][f][g] << "  "
             << C_two[0][2][2][0][f][g] << "  " << C_two[0][2][2][1][f][g] << "  " << C_two[0][2][2][2][f][g] << "\n";

        cout << C_two[1][0][0][0][f][g] << "  " << C_two[1][0][0][1][f][g] << "  " << C_two[1][0][0][2][f][g] << "  "
             << C_two[1][0][1][0][f][g] << "  " << C_two[1][0][1][1][f][g] << "  " << C_two[1][0][1][2][f][g] << "  "
             << C_two[1][0][2][0][f][g] << "  " << C_two[1][0][2][1][f][g] << "  " << C_two[1][0][2][2][f][g] << "\n";

        cout << C_two[1][1][0][0][f][g] << "  " << C_two[1][1][0][1][f][g] << "  " << C_two[1][1][0][2][f][g] << "  "
             << C_two[1][1][1][0][f][g] << "  " << C_two[1][1][1][1][f][g] << "  " << C_two[1][1][1][2][f][g] << "  "        
             << C_two[1][1][2][0][f][g] << "  " << C_two[1][1][2][1][f][g] << "  " << C_two[1][1][2][2][f][g] << "\n";       

        cout << C_two[1][2][0][0][f][g] << "  " << C_two[1][2][0][1][f][g] << "  " << C_two[1][2][0][2][f][g] << "  "
             << C_two[1][2][1][0][f][g] << "  " << C_two[1][2][1][1][f][g] << "  " << C_two[1][2][1][2][f][g] << "  "
             << C_two[1][2][2][0][f][g] << "  " << C_two[1][2][2][1][f][g] << "  " << C_two[1][2][2][2][f][g] << "\n";

        cout << C_two[2][0][0][0][f][g] << "  " << C_two[2][0][0][1][f][g] << "  " << C_two[2][0][0][2][f][g] << "  "
             << C_two[2][0][1][0][f][g] << "  " << C_two[2][0][1][1][f][g] << "  " << C_two[2][0][1][2][f][g] << "  "
             << C_two[2][0][2][0][f][g] << "  " << C_two[2][0][2][1][f][g] << "  " << C_two[2][0][2][2][f][g] << "\n";

        cout << C_two[2][1][0][0][f][g] << "  " << C_two[2][1][0][1][f][g] << "  " << C_two[2][1][0][2][f][g] << "  "
             << C_two[2][1][1][0][f][g] << "  " << C_two[2][1][1][1][f][g] << "  " << C_two[2][1][1][2][f][g] << "  "
             << C_two[2][1][2][0][f][g] << "  " << C_two[2][1][2][1][f][g] << "  " << C_two[2][1][2][2][f][g] << "\n";

        cout << C_two[2][2][0][0][f][g] << "  " << C_two[2][2][0][1][f][g] << "  " << C_two[2][2][0][2][f][g] << "  "
             << C_two[2][2][1][0][f][g] << "  " << C_two[2][2][1][1][f][g] << "  " << C_two[2][2][1][2][f][g] << "  "
             << C_two[2][2][2][0][f][g] << "  " << C_two[2][2][2][1][f][g] << "  " << C_two[2][2][2][2][f][g] << "\n";

        fflush(stdout);
      }  
    }
  }




  double square_brackets[3][3][3][3] = {0.0};

  double vol = Cluster::cluster().UnitCellVolume();
  for (int j=0;j<3;j++) {     
    for (int k=0;k<3;k++) {
      for (int l=0;l<3;l++) {
        for (int m=0;m<3;m++) {
          for (int f=0;f<Natoms;f++) {      
            for (int g=0;g<Natoms;g++) {

              double m1 = AtomicMasses[f];
              double m2 = AtomicMasses[g];

              square_brackets[j][k][l][m] += 1/(8*pi*pi*vol) * sqrt(m1*m2) * C_two[j][k][l][m][f][g];

            }
          }
          cout << square_brackets[j][k][l][m] << "\t";
        }
      }
      cout << "\n";
    }
  }
}


void Supercell::Compute_Elastic_Round_Brackets() {

  Matrix C_zero_reduced(3*(Natoms-1),3*(Natoms-1));
  int na = int(Supercell_Size[0]);
  int nb = int(Supercell_Size[1]);
  int nc = int(Supercell_Size[2]);

  Matrix gamma_dynamat = Cluster::cluster().GetHMBIHessian();
  Matrix mw_gamma_dynamat = Cluster::cluster().ComputeMassWeightedHessian(gamma_dynamat);
  Matrix POTR_gamma_dynamat = Cluster::cluster().ProjectOutTranslationAndRotation(mw_gamma_dynamat);



  for (int f=1;f<Natoms;f++) {
    for (int g=1;g<Natoms;g++) {
      for (int h=0;h<na*nb*nc;h++) {

        double mass_fac_f = 1.0, mass_fac_g = 1.0;  

        double m1 = AtomicMasses[f];  
        double m2 = AtomicMasses[g];               

        if (m1 == 0.0) {               
          mass_fac_f = 0;               
          m1 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }                    
        if (m2 == 0.0) {  
          mass_fac_g = 0;     
          m2 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }  

        for (int j=0;j<3;j++) {
          for (int k=0;k<3;k++) {

            C_zero_reduced.Element(3*(f-1)+j,3*(g-1)+k) += POTR_gamma_dynamat.Element(3*(f)+j,3*(g)+k);
//            C_zero_reduced.Element(3*(f-1)+j,3*(g-1)+k) += 1/( sqrt(m1*m2) ) * hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k);

          }
	}
      }
    }
  }

  C_zero_reduced.Inverse(); //this is the gamma matrix in Born and Huang (page 234)

  Matrix Gamma_matrix(3*Natoms,3*Natoms);
  for (int f=1;f<Natoms;f++) {      
    for (int g=1;g<Natoms;g++) {
      for (int j=0;j<3;j++) {       
        for (int k=0;k<3;k++) {     

          Gamma_matrix.Element(3*f+j,3*g+k) = C_zero_reduced.Element(3*(f-1)+j,3*(g-1)+k);

        }
      }
    }
  }



  double C_one[3][3][3][16][16] = {0.0};

  for (int f=0;f<Natoms;f++) {
    for (int g=0;g<Natoms;g++) {
      for (int h=0;h<na*nb*nc;h++) {

        Vector disp_vec(3);
        disp_vec[0] = Supercell_FractionalCoordinates[h][3*g+0]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+0];
        disp_vec[1] = Supercell_FractionalCoordinates[h][3*g+1]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+1];
        disp_vec[2] = Supercell_FractionalCoordinates[h][3*g+2]-Supercell_FractionalCoordinates[primitive_cell_index][3*f+2];

        double mass_fac_f = 1.0, mass_fac_g = 1.0;  

        double m1 = AtomicMasses[f];  
        double m2 = AtomicMasses[g];               

        if (m1 == 0.0) {               
          mass_fac_f = 0;               
          m1 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }                    
        if (m2 == 0.0) {  
          mass_fac_g = 0;     
          m2 += 1.0 ; // adding 1.0 to avoid 0/0 while mass-weighting  
        }  

 	for (int j=0;j<3;j++) {
          for (int k=0;k<3;k++) {
            for (int l=0;l<3;l++) {

              C_one[j][k][l][f][g] -= mass_fac_f*mass_fac_g*2*pi/( sqrt(m1*m2) )
                               * hess.Element( 3*Natoms*primitive_cell_index + 3*f+j,h*3*Natoms + 3*g + k)
                               * disp_vec[l];

            }                
          }
        }
      }
    }
  }




  double round_brackets[3][3][3][3] = {0.0};

  cout << "round brackets\n";
  double vol = Cluster::cluster().UnitCellVolume();
  for (int j=0;j<3;j++) {            
    for (int k=0;k<3;k++) {
      for (int l=0;l<3;l++) {
        for (int m=0;m<3;m++) {
          for (int f=0;f<Natoms;f++) {
            for (int g=0;g<Natoms;g++) {

              double m1 = AtomicMasses[f];
              double m2 = AtomicMasses[g];

              for (int p=0;p<3;p++) {  
                for (int q=0;q<3;q++) {

                  double sum1_parenthesis = 0.0, sum2_parenthesis = 0.0;
                  for (int r=0;r<Natoms;r++) {
                    double m3 = AtomicMasses[r];
                    sum1_parenthesis += sqrt(m3) * C_one[p][j][k][f][r];
                  }
                  for (int s=0;s<Natoms;s++) {
                    double m4 = AtomicMasses[s];  
                    sum2_parenthesis += sqrt(m4) * C_one[q][l][m][g][s];
       	       	  }

                  round_brackets[j][k][l][m] -= 1/(4*pi*pi*vol) * Gamma_matrix.Element(3*f+p,3*g+q)* sum1_parenthesis * sum2_parenthesis;


                }
              }

            }
          }
          cout << round_brackets[j][k][l][m]  << "\t";
        }
      }
      cout << "\n";
      fflush(stdout);
    }
  }
}


void Supercell::RunHMBISupercellTinkerFDHessianJob() {

  string job_path;
//  if (Params::Parameters().GetJobTypeStr() == "phonon") {  // run hessian (vibrate) job
    job_path = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      job_path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

    string cmd2 = "cd " + job_path;
    cmd2 += "; ";

    // Run the job
    cmd2 += "testhess supercell.xyz Y Y G 0.00001 Y > supercell.freq;";

    // Remove tmp file and extra geom file created by vibrate job
    cmd2 += "rm -f supercell.001; ";

    // Switch back to base directory
    cmd2 += "cd " + Params::Parameters().GetBasePath();

    //printf("Executing: %s\n",cmd2.c_str());
    system(cmd2.c_str());

//  }


}

// Main routine for reading in supercell hessian (similar to ReadHessian in cluster.C
// type 1 = qchem style, type 2 = tinker style
Matrix Supercell::ReadSupercellFDHessian(int type) {

  string path;
  path = Params::Parameters().GetHessianMMPath();

  //Path of the quasiharmonic calculations
  if(Params::Parameters().DoQuasiHarmonic() && Params::Parameters().SaveQuasiHarmonicCalculations() && !Params::Parameters().AreQHFAvailable()) 
      path += "/" + Quasiharmonic::quasiharmonic().GetHessianType();

  Matrix supercell_hess(3*Supercell_Natoms, 3*Supercell_Natoms);

  if ( Params::Parameters().NeglectManyBody() )         
    supercell_hess.Set();
  else {

    string filename;
    // Set up the filename, with the full path.  File is 'supercell.freq'

    if (type == 2) // Tinker MM job
      filename = path + "/supercell.freq";

    // Open the force file
    ifstream infile;
    infile.open( filename.c_str() );
    if ( !infile.is_open() ) {
      printf("Supercell::ReadHessian : Cannot open file '%s'\n",
             filename.c_str());
      exit(1);
    }

    if (type == 2) { // look in the tinker .freq file
      string line,line2;
      while ( !infile.eof() ) {

        getline(infile,line);
        // Search for the SCF hessian
        if ( line.substr(0,30)==" 3x3 Hessian Block for Atoms :" ) {
          istringstream iss(line);
          string tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
          int at1,at2;
          iss >> tmp1;
          iss >> tmp2;
       	  iss >> tmp3;
       	  iss >> tmp4;
       	  iss >> tmp5;
       	  iss >> tmp6;
          iss >> at1;
          iss >> at2;

          getline(infile,line2);  //blank line

          for (int i=0;i<3;i++) {
            string line1;
            getline(infile,line1);
            istringstream iss(line1);
            string tmp7;
            iss >> tmp7;
            iss >> supercell_hess.Element(3*(at1-1)+i,3*(at2-1)); // Store the hessian elements
       	    iss >> supercell_hess.Element(3*(at1-1)+i,3*(at2-1)+1); // Store the hessian elements  
       	    iss >> supercell_hess.Element(3*(at1-1)+i,3*(at2-1)+2); // Store the hessian elements  
          }

        }
      }
    }

    infile.close();  // Close the file

  }
  //PrintHessian("New Full MM Hessian", supercell_hess);

  // now FD hessians could be non-symmetric due to numerical noise, so symmetrize them

  for (int i=0;i<3*Supercell_Natoms;i++) {
    for (int j=i;j<3*Supercell_Natoms;j++) {
      supercell_hess.Element(i,j) += supercell_hess.Element(j,i); 
      supercell_hess.Element(i,j) *= 0.5;
    }
  }

  for (int i=1;i<3*Supercell_Natoms;i++) {
    for (int j=0;j<i;j++) {
      supercell_hess.Element(i,j) *= 0.0; 
      supercell_hess.Element(i,j) += supercell_hess.Element(j,i);
    }
  }

  supercell_hess.Scale(BohrToAng*BohrToAng/HartreesToKcalpermole);

  return supercell_hess; 
}       


void Supercell::ComputePDOS_old(Matrix phonon_matrix) {

  /*
    1) First find the max./min. freq in the phonon_matrix
    2) Then for each step of 1 cm^-1 in the range "max-min", find the number of freqs. (SIMPLE"S" LOL)
    3) On second thoughts, min. freq is not required. We always start with zero.
  */

  double maxfreq = 0.0;
  for (int i=0;i<total_k_points;i++) {
    for (int j=0;j<3*Natoms;j++) {
      if (phonon_matrix.Element(j,i) > maxfreq) {

        maxfreq = phonon_matrix.Element(j,i);
      }
    }
  }
  //cout << "maxfreq = " << maxfreq << "\n";

  cout << "Listing Density of States as function of frequency...\n"; 
  //following code to count the DOS for each freq. step could be more efficient??
  for (int k=0;k<(maxfreq+6.0);k=(k+5)) {

    int freq_count = 0;
    for (int i=0;i<total_k_points;i++) {
      for (int j=0;j<3*Natoms;j++) {      

        if ( (phonon_matrix.Element(j,i) >= k) && (phonon_matrix.Element(j,i) < (k+5)) ) {
          freq_count++;
        }
      }
    }
  }

}

void Supercell::ComputePDOS(Matrix phonon_matrix) {

  /*
    1) We are going to use the Normal Distribution curve to plot the PDOS in the following manner.
    2) Normal Distribution Function is a Gaussian: f(w) = (1/(sigma*sqrt(2*Pi))) exp( -(w-w0)**2 / (2*sigma**2) )
    3) So the idea really is to get an array of "w0" around which these gaussians are centered.
       This array needs to span all the frequencies "w" calculated using the Monkhorst-Pack scheme.
    4) Then for each "w0" we calculate and then add the contribution of each discrete calculated phonon "w" using the Normal Distribution Function.
    5) For starters, "w0" are separated by a constant stepsize, let's say " (w_max - w_min) / 30000" which is about 0.1 (cm)^-1  if "w_min=0" and "w_max=3000".
       Maybe we hard-code w0 = 0.1 (cm)^-1 .
  */


  //phonon_matrix.PrintHessian("phonon matrix");

  double maxfreq = 0.0;
  for (int i=0;i<total_k_points;i++) {
    for (int j=0;j<3*Natoms;j++) {
      if (phonon_matrix.Element(j,i) > maxfreq) {

        maxfreq = phonon_matrix.Element(j,i);
      }
    }
  }    
  //cout << "maxfreq = " << maxfreq << "\n";

  double minfreq = 1000.0; //arbitrarily
  for (int i=0;i<total_k_points;i++) {
    for (int j=0;j<3*Natoms;j++) {  
      if (phonon_matrix.Element(j,i) < minfreq) {

        minfreq = phonon_matrix.Element(j,i);
      }
    }
  }
  //cout << "minfreq = " << minfreq << "\n";

  double sigma = 5;       
  int num_w0 = int ((maxfreq - minfreq) / sigma) ; ;
  Vector w0(num_w0 + 20);

  //Now, start initializing the freq_centers "w0"
  // "w0" spans "minfreq - sigma" to "maxfreq + sigma"

  for (int i=0;i<(num_w0 + 20);i++) {

    if (i==0) {
      w0[i] = minfreq - 10*sigma ;
    }
    else {
      w0[i] = w0[i-1] + sigma;
    }

  } 

  Vector PDOS(num_w0 + 20);

  for (int i=0;i<(num_w0 + 20);i++) { 

    for (int j=0;j<total_k_points;j++) {         
      for (int k=0;k<3*Natoms;k++) {

        PDOS[i] += ( 1/sigma/sqrt(2.0*pi) ) * exp( - (phonon_matrix.Element(k,j) - w0[i]) * (phonon_matrix.Element(k,j) - w0[i]) / 2.0 / sigma / sigma ) ;

      }
    }
  }

  //Normalize the PDOS w.r.t. the total number of phonons
  PDOS.Scale(1.0/(total_k_points) );

  //cout << "Printing PDOS:\n";
  //for (int i=0;i<(num_w0 + 20);i++) {        
  //  cout << i << "  " << w0[i] << "  " << PDOS[i] << "\n";
  //}  

  //printf("\n\n\n");
  //fflush(stdout);

  string PDOSfile = Params::Parameters().GetHessianPath() + "/PDOS.dat";

  /*Create PDOS.dat file*/
  FILE *PDOSdata;
  if ((PDOSdata = fopen(PDOSfile.c_str(),"w"))==NULL) {
      printf("Supercell : Cannot open file '%s'\n",
        PDOSfile.c_str());
        exit(1);
  }

  fprintf(PDOSdata,"# Sigma =  %f \n", sigma);

  for (int i=0;i<(num_w0 + 20);i++) {
    fprintf(PDOSdata," %f  %f \n", w0[i], PDOS[i]);
  }  


}    

//Due to the way the dynamical matrix is created (6N*6N matrix separating the real and imaginary components),  
//the eigenvectors maps to double degenerate eigenvalues. Therefore the matrix and halved without losing any information.
//The imaginary component is than disgarded.
Matrix Supercell::GetOnlyRealMatrixComponent(Matrix FullMatrix,Vector eigvals){

  if(FullMatrix.GetRows() != FullMatrix.GetCols()){
    printf("Error::Supercell::GetOnlyRealMatrixComponent() FullMatrix not a square matrix. Row = %i Col = %i\n",
	   FullMatrix.GetRows(), FullMatrix.GetCols());
    exit(0);
  }

  if(FullMatrix.GetRows()%2 != 0){
    printf("Error::Supercell::GetOnlyRealMatrixComponent() FullMatrix has an odd number rows. Row = %i\n",
	   FullMatrix.GetRows());
    exit(0);
  }

  int FullRow = FullMatrix.GetRows();
  int RealRow = FullRow/2;
  Matrix RealMatrix(RealRow);

  
  //Of the two degenerate eignvectors the one greatest number of non-zero real elements is included in RealMatrix
  for(int i=0;i<RealRow;i++){
    int EvenSum = 0;
    int OddSum = 0;
   
    for(int j=0;j<RealRow;j++){
      if(fabs(FullMatrix(j,2*i)) > 0.00001){
	EvenSum++;
      }
      if(fabs(FullMatrix(j,2*i+1)) > 0.00001){
	OddSum++;
      }
    }
    

    //if there are degenerate modes, than imaginary and real component will not nessarily be next to each other.
    if(EvenSum == OddSum && EvenSum == 0){
      //FullMatrix.PrintHessian("FullMatrix");

      double thisFreq = eigvals[i]; 
      
      for(int j=0;j<RealRow;j++){
	double otherFreq = eigvals[j]; 
	//printf("thisFreq[%i] = %f otherFreq[%i] = %f\n",
	//     i,thisFreq,j,otherFreq);
	if((fabs(thisFreq - otherFreq) < 0.000001) && (thisFreq != otherFreq)){
	  //printf("Freq[%i] = %f and Freq[%i] %f match\n",i,thisFreq,j,otherFreq);
	  Vector OldOdd = FullMatrix.GetColumnVector(2*i+1);
	  Vector NewOdd = FullMatrix.GetColumnVector(2*j+1);
	  //FullMatrix.PrintHessian("FullMatrix before");
	  //OldOdd.Print("OldOdd");
	  //NewOdd.Print("NewOdd");
	  FullMatrix.SetColumnVector(NewOdd,2*i+1);
	  FullMatrix.SetColumnVector(OldOdd,2*j+1);
	  //FullMatrix.PrintHessian("FullMatrix");
	  //exit(0);

	  
	  break;
	}
      }
      for(int j=0;j<RealRow;j++)
	  RealMatrix(j,i) = FullMatrix(j,2*i+1);  

    }
    else{
      
      if (EvenSum >= OddSum)
	for(int j=0;j<RealRow;j++)
	  RealMatrix(j,i) = FullMatrix(j,2*i);
      else
	for(int j=0;j<RealRow;j++)
	  RealMatrix(j,i) = FullMatrix(j,2*i+1);
    }
  }
  //RealMatrix.PrintHessian("RealMatrix");  
  return RealMatrix;
}

//void Supercell::ComputeElasticConstants() {

//Print Vibrational Analysis and normal modes in Molden Format into a new inpt_filename-molden.out file
void Supercell::PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index, int k_index) {

  Vector atomic_coords = Cluster::cluster().GetCurrentCoordinates();

  stringstream ss;
  ss << k_index;
  string k_index_str = ss.str();

  /* Create the xyz file */
  FILE *xyz;
//  string infile = GetInputFilename();
//  string molden_file = infile.substr(0,infile.size()-3);
  string molden_file = Params::Parameters().GetHessianPath();
  molden_file += "/molden-";
  molden_file += k_index_str;
  molden_file += ".out";
  if ((xyz = fopen(molden_file.c_str(),"w"))==NULL) {
    printf("Cluster::PrintNormalModes : Cannot open file '%s'\n",
           molden_file.c_str());                        
    exit(1);
  }

  ofstream moldenstream;
  moldenstream.open(molden_file.c_str());
//  int Natoms = GetTotalNumberOfAtoms();
  //  moldenstream << "======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======\n";
  moldenstream << "[Molden Format]\n";           
  moldenstream << "[Atoms] (Angs)\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << Cluster::cluster().GetAtomicSymbol(i).c_str() << "  " <<  right << setw(5) << i+1 << "  "
         <<  right << setw(5) << Cluster::cluster().GetAtomicNumber(i) << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+1] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+2] << "\n";
  }

  moldenstream << "[FREQ]\n";
  for (int i=0; i < Nonzerofreqs; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  //zero frequency vibrations
  if (Params::Parameters().PrintLevel() > 0) {                 
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++) {
      moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
    }
  }

  moldenstream << "[FR-COORD]\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << Cluster::cluster().GetAtomicSymbol(i).c_str() << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+1]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+2]*AngToBohr << "\n";
  }

  moldenstream << "[FR-NORM-COORD]\n";               
  for (int i=0; i < Nonzerofreqs; i++)  {
    moldenstream << "vibration " << i+1 << "\n";        
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++)  {
    moldenstream << "vibration " << i+Nimagfreqs-Nfreq+Nonzerofreqs+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++)  {
      moldenstream << "vibration " << i+Nimagfreqs+1 << "\n";
      for (int k=0; k < Natoms; k++) {
        moldenstream <<  "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
      }
    }
  }
//   moldenstream << "\n======= END OF MOLDEN-FORMATTED INPUT FILE =======\n";              
  moldenstream.close();
} 

//Print Vibrational Analysis and normal modes in Molden Format a custom  directory
void Supercell::PrintNormalModes(Vector freq, int Nfreq, Matrix UnMWHess, int Nonzerofreqs, int Nimagfreqs, Vector orig_index, int k_index, string file_name) {

  Vector atomic_coords = Cluster::cluster().GetCurrentCoordinates();

  stringstream ss;
  ss << k_index;
  string k_index_str = ss.str();

  /* Create the xyz file */
  FILE *xyz;
//  string infile = GetInputFilename();
//  string molden_file = infile.substr(0,infile.size()-3);
  string molden_file = Params::Parameters().GetHessianPath();
  molden_file += "/" + file_name;
  //molden_file += "/molden-";
  molden_file += k_index_str;
  molden_file += ".out";
  if ((xyz = fopen(molden_file.c_str(),"w"))==NULL) {
    printf("Cluster::PrintNormalModes : Cannot open file '%s'\n",
           molden_file.c_str());                        
    exit(1);
  }

  ofstream moldenstream;
  moldenstream.open(molden_file.c_str());
//  int Natoms = GetTotalNumberOfAtoms();
  //  moldenstream << "======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======\n";
  moldenstream << "[Molden Format]\n";           
  moldenstream << "[Atoms] (Angs)\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << Cluster::cluster().GetAtomicSymbol(i).c_str() << "  " <<  right << setw(5) << i+1 << "  "
         <<  right << setw(5) << Cluster::cluster().GetAtomicNumber(i) << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+1] << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+2] << "\n";
  }

  moldenstream << "[FREQ]\n";
  for (int i=0; i < Nonzerofreqs; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++) {
    moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
  }
  //zero frequency vibrations
  if (Params::Parameters().PrintLevel() > 0) {                 
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++) {
      moldenstream << setprecision(8) << showpoint <<  right << setw(16) << freq[i] << "\n";
    }
  }

  moldenstream << "[FR-COORD]\n";
  for (int i = 0; i < Natoms; i++) {
    moldenstream << right << setw(5) << Cluster::cluster().GetAtomicSymbol(i).c_str() << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+1]*AngToBohr << "  "
         << setprecision(8) << showpoint <<  right << setw(15) << atomic_coords[3*i+2]*AngToBohr << "\n";
  }

  moldenstream << "[FR-NORM-COORD]\n";               
  for (int i=0; i < Nonzerofreqs; i++)  {
    moldenstream << "vibration " << i+1 << "\n";        
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  for (int i=Nfreq-Nimagfreqs; i < Nfreq; i++)  {
    moldenstream << "vibration " << i+Nimagfreqs-Nfreq+Nonzerofreqs+1 << "\n";
    for (int k=0; k < Natoms; k++) {
      moldenstream <<  "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
           << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
    }
  }
  if (Params::Parameters().PrintLevel() > 0) {
    for (int i=Nonzerofreqs; i<Nfreq-Nimagfreqs; i++)  {
      moldenstream << "vibration " << i+Nimagfreqs+1 << "\n";
      for (int k=0; k < Natoms; k++) {
        moldenstream <<  "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+1, int(orig_index[i]) ) << "  "
             << setprecision(8) << showpoint <<  right << setw(15) << UnMWHess.Element(k*3+2, int(orig_index[i]) ) << "\n";
      }
    }
  }
//   moldenstream << "\n======= END OF MOLDEN-FORMATTED INPUT FILE =======\n";              
  moldenstream.close();
}  

//make an input file from the super cell
void Supercell::CreateSupercellInputFile(){

  //file name
  FILE *input;
  string input_file = "supercell.in";
  if ((input = fopen(input_file.c_str(),"w"))==NULL) {
    printf("CreateSupercellInputFile : Cannot open file '%s'\n",input_file.c_str());
    exit(1);
  } 
  
  fprintf(input,"%s\n\n",Params::Parameters().GetHMBIRem().c_str());
  fprintf(input,"$molecule\n");
  //assume that the spin of the super cell is the same as the cluster
  fprintf(input,"%i %i\n",Ncells*Cluster::cluster().GetClusterCharge(),
	  Cluster::cluster().GetClusterSpin());

  if(Params::Parameters().GetMMType() == 1){
    
    //print monomer
    //fprintf(input,"found tinker input.\n Supercell_NMon = %i\n",Supercell_NMon);
    for(int iMon=1;iMon<=Supercell_NMon;iMon++){
      
      int ref_mon = Supercell_Monomers[iMon].GetReferenceMonomerIndex();
      //fprintf(input,"ref mon = %i\n",ref_mon);
      
      fprintf(input,"--\n%i %i\n", Cluster::cluster().GetMonomerChange(ref_mon),
	      Cluster::cluster().GetMonomerSpin(ref_mon));
      for(int iatom=0;iatom<Supercell_Monomers[iMon].GetNumberOfAtoms();iatom++){
	Supercell_Monomers[iMon].GetAtom(iatom).PrintTinkerCartesian(0,input);
      }
    }
      
    fprintf(input,"$end\n\n");

    // Print $qchem section - erase "$rem" and replace it with "$qchem"
    string tmp = Params::Parameters().GetQChemRem();
    tmp.erase(0,5);
    fprintf(input,"$qchem%s", tmp.c_str() );


    fprintf(input,"$tinker%s$end\n",Params::Parameters().GetTinkerRem().c_str());
    double a,b,c,alpha,beta,gamma;
    GetSuperCellLatticeParameters(a,b,c,alpha,beta,gamma);
    fprintf(input,"\n$unit_cell\n%f %f %f\n%f %f %f\n$end",
	    a,b,c,alpha,beta,gamma);
 
  }else{
    printf("ERROR::Supercell::CreateSupercellInputFile() : Not adapted for any other type other than tinker\n");
    exit(0);
  }
  
  fclose(input);

}
