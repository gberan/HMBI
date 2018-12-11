#ifndef _vector_h
#define _vector_h

#include <iostream>
#include <ostream>
#include <math.h>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::string;

/* 
   A simple vector class

   GJB 2/09
*/

class Vector {

  int dim;
  double *vec;
  bool allocated_;

 public:
  // Constructors & Destructor
  Vector();// Empty constructor, no allocation
  Vector(double *data, int N); // Standard constructor
  Vector(int N); // Initialize to zero vector of length N
  Vector(const Vector& other, bool if_data=true); // Copy constructor
  ~Vector();

  void Initialize(int N);
  void Initialize(double *data, int N);

  // Initialize vector elements to 0.0
  void Set() {for (int i=0;i<dim;i++) {vec[i] = 0.0;}};

  // Initialize vector so every element in vector has the same value
  void Set(double value) {for(int i=0;i<dim;i++) {vec[i] = value;}};

  // Return the number of vector elements
  int GetLength() const {return dim;};
  // Return pointer to the full vector, use with care
  double *GetVectorPtr() const {return vec;};
  // Return single vector element, counting from 0 -> dim-1
  double& Element(int i) {return vec[i];};
  // Nicer routine for returning element
  double& operator[](int i) const
  {
    if (i > dim) {
      printf("ERROR: Vector element %d out of range (dim = %d)\n",i,dim);
      exit(1);
    }
    return vec[i];
  }

  // Sort in descending order
  void SortLargestFirst();
  //returns the mapping of the original locations in the vector
  Vector SortLargestFirstAndReturn();

  //Finds the FIRST location an element in the Vector
  //returns the element number
  int Find(double find);

  // Print out vector
  void Print(string title);
  void PrintGradient(string title);

  // Compute Dot product of vector with another
  double DotProduct(const Vector& other);
  // Compute Outer product of this vector with another: this * other
  //Matrix OuterProduct(const Vector& other);
  // Compute Cross product of this vector with another: this x other
  Vector CrossProduct(const Vector& other);

  // Return the 2-norm of the vector
  double Norm();
  // Return Norm/sqrt(dim)
  double RMS();

  double Max(bool AbsVal=false);
  double Min(bool AbsVal=false);
  
  // Multiply Vector by a scalar
  void Scale(double coeff) const;
  // Normalize vector
  void Normalize();
  // Count number of positive elements
  int Nonzero();
  //Count number of non positive elements
  int Nonpositive(){return GetLength() - Nonzero();};

  // Rotate vector by 'angle' degrees about an arbitrary 'axis'
  Vector RotateAboutAxis3D( double angle, Vector axis);

  // Set vector with incrementally increasing entries from a starting value
  void Incremental(int Start = 1);

  //Vector with incremental integer values, increasing or decreasing
  void IncrementalEntries(bool increasing);

  // Overload operators
  Vector& operator=(const Vector& other);
  Vector& operator+=(const Vector& other);
  Vector& operator-=(const Vector& other);

  Vector& operator*=(double coeff); // untested, may not work
  Vector& operator*(double coeff); //untested, may not work

  //friend Vector operator+( Vector V1, Vector V2)

};

#endif
