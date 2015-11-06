#include "vector.h"
#include "constants.h"
using namespace hmbi_constants;

// Empty constructor
Vector::Vector() : vec(NULL), allocated_(false) {

}

// Constructor: Initialize with data
Vector::Vector(double *data, int N) : dim(N), vec(NULL), allocated_(true) {
  vec = new double[dim];
  memcpy(vec,data,sizeof(double)*dim);

}

// Constructor: Initialize blank vector
Vector::Vector(int N) : dim(N), vec(NULL), allocated_(true) {
  vec = new double[dim];
  Set(); // initialize to zero
}

// Copy constructor
Vector::Vector(const Vector &other, bool if_data) : dim(other.dim), 
  vec(NULL), allocated_(true) {

  vec = new double[dim];
  if (if_data)
    memcpy(vec,other.vec,sizeof(double)*dim);
  else
    Set();
}

// Destructor
Vector::~Vector() {
  if (allocated_) 
      delete [] vec;
}

// Create empty vector, alternative to the constructor
void Vector::Initialize(int N) {
  dim = N;
  vec = new double[dim];
  allocated_ = true;
  Set(); // initialize to zero
}

// Create vector from standard C array
void Vector::Initialize(double *data, int N) {
  dim = N;
  vec = new double[dim];
  memcpy(vec,data,sizeof(double)*dim);
  allocated_ = true;
}

// Return vector containing a subset of sequential elements from first
// to last, inclusive.
Vector Vector::GetRange(int first, int last) {

  //printf("Vector::GetRange(): dim=%d, first=%d, last=%d\n",dim,first,last);  

  if (last <= first) {
    printf("ERROR: Vector::GetRange(): Last element (%d) must be greater than first one (%d).\n",first,last);
    exit(1);
  }

  if (last < 0 || first < 0) {
    printf("ERROR: Vector::GetRange(): First (%d) and last (%d) element must be positive integers.\n",first,last);
    exit(1);
  }

  if (first >= dim || last >= dim) {
    printf("ERROR: Vector::GetRange(): First (%d) and last (%d) element must be within the range 0->%d.\n",first,last,dim-1);
    exit(1);
  }

  int size = last - first + 1;
  Vector subset(size);
  int j=0;
  for (int i=first; i<=last; i++) {
    subset[j] = vec[i];
    j++;
  }
  return subset;
}

// Sort the vector from largest to smallest.
void Vector::SortLargestFirst() {

  // Create a copy;
  double *tmp = new double[dim];
  memcpy(tmp,vec,sizeof(double)*dim);

  double max;
  int jmax;

  // Begin sorting
  for (int i=0;i<dim;i++) {
    max = -9999999.0;
    jmax = -1;
    for (int j=0;j<dim;j++) {
      if (tmp[j] > max) {
	max = tmp[j];
	jmax = j;
      }
    }
    vec[i] = max;
    tmp[jmax] = -9999999.0;
  }

}

// Print out Vector
void Vector::Print(string title) {
  printf("%s    Length = %d\n",title.c_str(),dim);
  for (int i=0;i<dim;i++) 
    printf("%.9f\n",vec[i]);
}

// Print out Vector horizontally, N_per_line numbers per line
void Vector::PrintTranspose(string title, int N_per_line) {
  printf("%s    Length = %d\n",title.c_str(),dim);

  int current_line = 0;
  for (int i=0;i<dim;i++) {
    current_line++;
    printf("%10.6f ",vec[i]);
    if (current_line % N_per_line == 0 || (i+1)==dim)
      printf("\n");
  }
  printf("\n");
}


// Print Gradient Vector with xyz columns
void Vector::PrintGradient(string title) {

  if (dim%3) {
    printf("ERROR: Vector::PrintGradient(): dim = %d, dim_mod(3) = %d\n",dim,dim%3);
    exit(1);
  }

  printf("%s  Length = %d\n",title.c_str(),dim);
  for (int i=0;i<dim/3;i++) 
    printf("%.6f %.6f %.6f\n",vec[3*i],vec[3*i+1],vec[3*i+2]);

}

double Vector::DotProduct(const Vector& other) {
  double dot = 0.0;
  for (int i=0;i<dim;i++)
    dot += vec[i]*other.vec[i];
  return dot;
}

Vector Vector::CrossProduct(const Vector& other) {
  // First check that both vectors are 3-D
  if ( dim!=3  || other.dim!=3 ) {
    printf("ERROR: Vector::CrossProduct(): Both vector lengths should be 3\n");
    exit(1);
  }

  Vector prod(3);
  prod[0] = vec[1]*other.vec[2] - vec[2]*other.vec[1];
  prod[1] = vec[2]*other.vec[0] - vec[0]*other.vec[2];
  prod[2] = vec[0]*other.vec[1] - vec[1]*other.vec[0];
  return prod;
}


// Compute 2-norm of a vector
double Vector::Norm() {

  double norm = sqrt(DotProduct(*this));
  return norm;
}

double Vector::RMS() {
  double rms = Norm()/sqrt(1.0*dim);
  return rms;
}

// Find maximum value in vector.  If AbsVal=true, use absolute values
double Vector::Max(bool AbsVal) {

  double max = 0.0;

  for (int i=0;i<dim;i++) {
    double value;
    if (AbsVal) 
      value = fabs(vec[i]);
    else
      value = vec[i];
    if (value > max && !AbsVal)
      max = value;
    if (value > fabs(max) && AbsVal)
      max = vec[i];
  }

  //printf("Max value in vector = %f\n",max);
  return max;
}


void Vector::Scale(double coeff) const {
  // could do faster version, BLAS DSCAL(?)
  for (int i=0;i<dim;i++)
    vec[i] *= coeff;
}

// Normalize vector
void Vector::Normalize() {
  double norm = Norm();
  if (norm > 1.0e-8) 
    Scale(1.0/norm);
  else {
    printf("Vector::Normalize(): cannot normalize, norm = %f\n",norm);
    exit(1);
  }
}

int Vector::Nonzero() {

  int nonzero = 0;
  double thresh = 1.0e-12;
  for (int i=0;i<dim;i++) {
    if (vec[i] > thresh) nonzero++;
  }
  return nonzero;
}

Vector Vector::RotateAboutAxis3D( double angle, Vector axis) {

  // Convert angle to radians & evaluate some trigonometry functions
  double theta = angle*DegreesToRadians;
  double c = cos(theta);
  double s = sin(theta);
  double t = 1.0 - cos(theta);

  // Normalize the axis of rotation
  axis.Normalize();
  double X = axis[0];
  double Y = axis[1];
  double Z = axis[2];

  // Create the rotation matrix
  double RotMat[3][3]; // 3x3 rotation matrix
  RotMat[0][0] = t*X*X + c;
  RotMat[0][1] = t*X*Y + s*Z;
  RotMat[0][2] = t*X*Z - s*Y;

  RotMat[1][0] = t*X*Y - s*Z;
  RotMat[1][1] = t*Y*Y + c;
  RotMat[1][2] = t*Y*Z + s*X;

  RotMat[2][0] = t*X*Z + s*Y;
  RotMat[2][1] = t*Y*Z - s*X;
  RotMat[2][2] = t*Z*Z + c;

  // Find the rotated vector: new = RotMat*vec
  Vector newvec(3);  // initialized to zero
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) 
      newvec[i] += RotMat[i][j]*vec[j];
  
  //newvec.Print("Rotated Vector");

  return newvec;
}

// Vector addition.  Can handle vectors of different sizes by effectively
// padding the shorter vector with zeroes at the end.
Vector Vector::Add(const Vector& other, bool size_mismatch_ok) {

  if (!size_mismatch_ok && dim != other.dim) {
    printf("Vector::Add() Error: vector size mismatch!  Either make sizes match or set size_mismatch_ok=true\n");
    exit(1);
  }
  else {
    // Find dimension of new vector as larger of two vectors being added
    int new_dim;
    if (dim >= other.dim) new_dim = dim;
    else new_dim = other.dim;
    Vector Result(new_dim);

    // Add the two vectors
    for (int i=0;i<dim;i++)
      Result[i] = vec[i];
    for (int i=0;i<other.dim;i++)
      Result[i] += other[i];
    return Result;
  }

}


// Vector addition.  Adds vector "other" to the current vector.  Can
// handle vectors of different sizes, as long as the "other" vector is
// the shorter one, in which case it just adds the shorter one to the
// first N elements of the longer one.
Vector Vector::AddTo(const Vector& other, bool size_mismatch_ok) {

  if (!size_mismatch_ok && dim != other.dim) {
    printf("Vector::Add() Error: vector size mismatch!  Either make sizes match or set size_mismatch_ok=true\n");
    exit(1);
  }

  if (dim < other.dim) {
      printf("Vector::AddTo() Error: other vector cannot be longer than result vector.\n");
      exit(1);
  }

  for (int i=0;i<other.dim;i++) {
    vec[i] += other.vec[i]; // BLAS call could be faster
  }

  return *this;
}



// Overload "=" operator
Vector& Vector::operator=(const Vector& other) {


 // Copy over the data
  if (this!=&other) {
    // Initialize the vector, if needed

    //if (!allocated_ ) 
    //  Initialize(other.dim);
    if (other.dim > 0) {
      delete [] vec;
      Initialize(other.dim);
      // Copy over the data
      memcpy(vec,other.vec,sizeof(double)*dim);
    }
    else {
      dim = 0;
      vec = NULL;
      allocated_ = false;
    }
  }
  return *this;
} 

// Overload "+=" operator
Vector& Vector::operator+=(const Vector& other) {
  if (dim == other.dim) {
    for (int i=0;i<dim;i++)
      vec[i] += other.vec[i]; // BLAS call could be faster
  }
  else {
    printf("Vector::+= Error: vector size mismatch!  Use Vector::Add to combine vectors of different lengths\n");
    exit(1);
  }
  return *this;
}

// Overload "-=" operator
Vector& Vector::operator-=(const Vector& other) {

  if (dim == other.dim) {
    for (int i=0;i<dim;i++)
      vec[i] -= other.vec[i]; // BLAS call could be faster
  }
  else {
    printf("Vector::-= Error: vector size mismatch!\n");
    exit(1);
  }
  return *this;
}

// Overload "*=" operator to scalar multiplication
Vector& Vector::operator*=(double coeff) {
  // could do faster version, BLAS DSCAL(?)
  if (allocated_) {
    for (int i=0;i<dim;i++)
      vec[i] *= coeff;
  }
  else{
    printf("Vector::*= Error: vector not allocated!\n");
    exit(1);
  }
  return *this;
}

