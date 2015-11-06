#ifndef _matrix_h
#define _matrix_h

#include <iostream>
#include <math.h>
#include <string>
//#include <vector>
#include "vector.h"
using std::string;

// LAPACK routines:
extern "C" void dgemm_(char *TransA, char *TransB, int *M, int *N, int *K,
		       double *alpha, double *A, int *LDA, double *B,
		       int *LDB, double *beta, double *C, int *LDC);

extern "C" void dsyev_(char *job,  char *type, int *N, double *A, int *LDAs, 
		       double* lambda, double *work, int *lwork, int *info); 

extern "C" void dgetrf_(int *M, int *N, double *A, int *LDA, double *IPIV, 
			int *INFO);

extern "C" void dgetri_(int *N, double *A, int *LDA, double *IPIV, 
			double *work, int *lwork, int *info);

extern "C" void dgesvd_(char *jobu, char *jobvt, int *M, int *N, double *A, 
			int *LDA, double *S, double *U, int *LDU, double *VT,
			int *LDVT, double *work, int *lwork, int *info);

extern "C" void dgesv_(int *N, int *NHRS, double *A, int* LDA, int *ipiv, 
		       double* B, int* LDB, int *info);

extern "C" void dgesvx_(char *FACT, char *TRANS, int *N, int *NHRS, double *A, 
			int* LDA, double *AF, int *LDAF, int *IPIV, char *EQUED, 
			double *R, double *C, double *B, int *LDB, double *X, 
			int* LDX, double *RCOND, double *FERR, double *BERR, 
			double *WORK, int* IWORK, int *INFO);



/*
Matrix class
*/

class Matrix {
  
   int rows;
   int cols;
   double *matrix;
   bool allocated_;

   // private copy constructor to prevent compiler from creating its own
   //Matrix(const Matrix& other);

 public:
  // Constructors & Destructor
  Matrix(); //Empty constructor, no allocation
  Matrix(double *data, int M, int N); //Standard
  Matrix(int M, int N); //initialize all values to zero with given dimensions

  Matrix(const Matrix& other, bool if_data=true); //copy constructor
  //Matrix(int N); // square matrix of dim N
  Matrix(int N, bool if_iden=false); //square (or identity) matrix
//  Matrix(const Vector& other, bool if_iden=true);
  ~Matrix();

  void Initialize(int M, int N);
  void Initialize(const Matrix& other);


  bool IsAllocated() {return allocated_;};

  //initialize matrix elements to 0.0

  void Set() {
    for (int x=0; x < rows*cols; x++) 
      matrix[x] = 0.0;
  };

  void Set_Iden() {
    if (rows != cols) {
      printf("ERROR: Matrix::Set_Iden(): matrix is not square\n");
      exit(1);
    }
    Set();
    for (int x=0; x < rows; x++) 
      matrix[x + rows*x] = 1.0;
  };


  //Return the dimensions of the vector elements
  int GetDimensions() {return rows, cols;};

  int GetRows() const {return rows;};
  int GetCols() const {return cols;};
  
  //return pointer to the full matrix
  double *GetMatrixPtr() const {return matrix;};

  // return single vector element
  double& Element(int i, int j) {return matrix[rows*j + i];};

  double& operator()(int i, int j) const
  {
    if (i >= rows || j >= cols) {
      printf("ERROR: Matrix element (%d,%d) is out of bounds.\n",i,j);
      exit(1);
    }
    return matrix[rows*j + i];
  } 
  //nicer routine for getting element?? did not quite understand  

  // print out matrix
  void Print(string title);

  /* Linear algebra routines */

  // matrix multiplication with matrix
  Matrix Multiply(Matrix& B, int type=1);
  // Diagonalization.  If computed, eigenvectors overwrite original matrix.
  // The eigenvalues are stored in the output Vector.
  Vector FindEigenvalues(); // gets eigenvalues only
  Vector Diagonalize(); // gets eigenvalues & eigenvectors (as column vectors).
  Vector SVD(Matrix &U, Matrix &VT, bool find_condition=false); // computes Singular Value Decomposition
  void Inverse(); // compute matrix inverse for square matrix
  void SolveLinearEquations(Vector &b); // solve A*x = b, solution vector overwrites b.  Destroys A.
  void ExpertSolveLinearEquations(Vector &b); // solve A*x = b, numerically fancier version.  


  double TwoNorm(); // compute matrix 2 norm

  //scale with a scalar
  void Scale(double coeff) const;
  // Overwrite matrix with its transpose
  void Transpose();

  

  // dot product of vectors can be a special case of 1-D matrices
 
  //overload operators
  Matrix& operator=(const Matrix& other);
  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  Matrix& operator*=(double coeff);

  Matrix operator+(const Matrix& other); //test this
  Matrix operator-(const Matrix& other); //test this
  Matrix operator*(Matrix& other); //test this
//  Matrix& operator=(const Vector& other, int i);
//  Matrix& operator=(Vector other);
//  Matrix& RowEqualsVec(Matrix& matrix, Vector& vec, int i);

  Vector MatrixTimesVector(const Vector& b);
  // Extract a vector from the matrix
  Vector GetColumnVector(int icol);
  Vector GetRowVector(int irow);

  // Input vector from the matrix
  void SetColumnVector(Vector vec, int icol);
  void SetRowVector(Vector vec, int irow);

  // Set a block of a matrix, with top left corner at (irow,icol)
  Matrix GetBlock(int nrows, int ncols, int row_offset, int col_offset);
  void SetBlock(const Matrix& block, int row_offset, int col_offset);
  void AddBlockToMatrix(const Matrix& block, int row_offset, int col_offset);

  void ReorderRows(int i, int j);
  void ReorderColumns(int i, int j);


};

#endif
