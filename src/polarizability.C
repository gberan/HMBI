#include <string>
using std::string;
#include "polarizability.h"

Polarizability::Polarizability() : Rank(-1), Npol(0) {
  // empty constructor  
}

// primary constructor
Polarizability::Polarizability(int rank, double* pols) : Rank(rank), Npol(0) {

  // Figure out how many spherical tensor moments this rank corresponds to.
  Npol = 0;
  for (int i=0;i<=rank;i++) 
    Npol += 2*i + 1;
  
  // Initialize the moments vector
  Matrix tmp(pols,Npol,Npol);

  // Just as for polarizabilities, we use a slightly different order
  // to store them than CamCasp does.  We use xyz order, instead of
  // 10, 11c, 11s order.  To get this, reorder rows and columns 1,2,3
  // -> 2,3,1.
  if (Rank > 0) {
    // Start with [0,1,2,3,...]
    tmp.ReorderRows(1,2); // [0,2,1,3,...]
    tmp.ReorderRows(2,3); // [0,2,3,1,...]
    
    // Repeat for columns
    // Start with [0,1,2,3,...]
    tmp.ReorderColumns(1,2); // [0,2,1,3,...]
    tmp.ReorderColumns(2,3); // [0,2,3,1,...]

  }

  Pols.Initialize(tmp);
}

// Copy constructor
Polarizability::Polarizability(const Polarizability& other, bool copy_data) :
  Rank(other.Rank), Npol(other.Npol) {

  // Initialize the vectors
  if (Npol > 0) {
    Pols.Initialize(Npol);
    
    if (copy_data) 
      Pols = other.GetPolarizabilities();
    else 
      Pols.Set();
  }
}

Polarizability::~Polarizability() {
 
}

// Overload the '=' operator
Polarizability& Polarizability::operator=(const Polarizability& other) {

  if (this!=&other) {
    // Read the ranks
    Rank = other.Rank;
    Npol = other.Npol;
    
    // Initialize the moments
    if (Npol > 0) {
      //Pols.Initialize(Npol,Npol); // initialization should be unneccessary now,
                                  // but I haven't tested that
      Pols = other.Pols;
    }
  }

  return *this;
}

// Print out the polarizability 
void Polarizability::Print(string title) {  


  if (Rank < 0) {
    printf("Error: Polarizability::Print(): Unknown Rank = %d\n",Rank); 
    exit(1);
  }

  printf("%s:  Rank %d\n",title.c_str(),Rank);
  Pols.Print("");

}
