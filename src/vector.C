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

// Overloaded sort function that returns the original
// location of each element
// Sort the vector from largest to smallest.
Vector Vector::SortLargestFirstAndReturn() {

  // Create a copy;
  double *tmp = new double[dim];
  Vector orig_location(dim);
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
    orig_location[i] = jmax;
    tmp[jmax] = -9999999.0;
  }

  return orig_location;

}

// Find the first matching element in Vector
// By default returns the first element if there is no match
int Vector::Find(double find) {
  bool found=false;
  int elem=0;
  for (int i=0;i<dim;i++) {
    if(vec[i]==find){
      //printf("Found element %i\t at position %i\n",vec[i],i);
      elem=i;
      found=true;
      break;
    }
  }

  if(!found){
    printf("Warning match not found!!! Recheck the vector\nExiting HMBI");
    exit(0);
  }

  return elem;
}

// Print out Vector
void Vector::Print(string title) {
  printf("%s    Length = %d\n",title.c_str(),dim);
  for (int i=0;i<dim;i++) 
    printf("%.6f\n",vec[i]);
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


// Find minimum value in vector.  If AbsVal=true, use absolute values
double Vector::Min(bool AbsVal) {

  double min = 1000000.0;

  for (int i=0;i<dim;i++) {
    double value;
    if (AbsVal) 
      value = fabs(vec[i]);
    else
      value = vec[i];
    if (value < min && !AbsVal)
      min = value;
    if (value < fabs(min) && AbsVal)
      min = vec[i];
  }

  //printf("Max value in vector = %f\n",max);
  return min;
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

  // printf("RotMat\n");
  //for(int i=0;i<3;i++){
  //  printf("%f %f %f\n",
  //	   RotMat[i][0],RotMat[i][1],RotMat[i][2]);
  //}

  // Find the rotated vector: new = RotMat*vec
  Vector newvec(3);  // initialized to zero
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) 
      newvec[i] += RotMat[i][j]*vec[j];
  
  //newvec.Print("Rotated Vector");
  fflush(stdout);

  return newvec;
}

// Set vector with incrementally increasing entries from a starting value
void Vector::Incremental(int Start){

  for(int i=0;i<dim;i++){
    vec[i] = double (i+Start);
  }

}

// Set vector with incrementally increasing or decreasing values.
void Vector::IncrementalEntries(bool increasing){

  for(int i = 0; i < dim;i++)
    if(increasing)
      vec[i] = double(i);
    else
      vec[dim-i] = double(i); 
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
    printf("Vector::+= Error: vector size mismatch!\n");
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

//Overload "*" operator for scalar multiplication
Vector& Vector::operator*( double coeff) {
  //Vector new_vec(dim);

  for (int x=0;x<dim;x++) {
    vec[x] *= coeff;
  }
  return *this; 
}
