#include <iostream.h>
#include "matrix.h"
#include <vector>
#include "vector.h"
using namespace std;


//Empty Constructor (Why??)
Matrix::Matrix() : rows(0), cols(0), matrix(NULL), allocated_(false) {

}

//Constructor: Initialize with data
Matrix::Matrix(double *data, int M, int N) : rows(M), cols(N), matrix(NULL), allocated_(true) {

  matrix =  new double[rows*cols];
  memcpy(matrix, data, sizeof(double)*rows*cols);
}


//Constructor: Initialize blank matrix
Matrix::Matrix(int M, int N) :  rows(M), cols(N), matrix(NULL), allocated_(true) {

  matrix = new double[rows*cols];
  Set();  // defined in matrix.h
}

//Copy constructor
Matrix::Matrix(const Matrix &other, bool if_data) : rows(other.rows), cols(other.cols), matrix(NULL), allocated_(other.allocated_)
{ 
  if (rows*cols > 0) {  
    matrix =  new double[rows*cols];
    if (if_data)
      memcpy(matrix, other.matrix, sizeof(double)*rows*cols);
    else
      Set();
  }
}


//a square matrix initialized to zero (default) or identity.
Matrix::Matrix(int N, bool if_iden) : rows(N), cols(N), matrix(NULL), allocated_(true) {
  matrix =  new double[rows*cols];

  if (if_iden) {
    Set_Iden(); // function defined in matrix.h for creating identity matrix
  }
  else
    Set();
}


//Destructor
Matrix::~Matrix() {
  if (allocated_) 
    delete [] matrix;
}


//Create empty Matrix, alternative to the constructor
void Matrix::Initialize(int M, int N) {
  rows = M;
  cols = N;

  matrix =  new double[rows*cols];
  allocated_ = true;

  Set(); //initialize to zero
}

//Create empty Matrix, alternative to the constructor
void Matrix::Initialize(const Matrix& other) {
  rows = other.rows;
  cols = other.cols;

  matrix = new double[rows*cols];
  allocated_ = true;
  memcpy(matrix, other.matrix, sizeof(double)*rows*cols);

}

// Print out Matrix
void Matrix::Print(string title) {
 
  printf("%s    ROWS = %d     COLUMNS = %d\n", title.c_str(), rows, cols); //what is title.c_str()??
  for (int x=0;x<rows;x++) {
    for (int y=0;y<cols;y++) {
      printf("%8.5f ",Element(x,y));
      //printf("%5.2f ",Element(x,y));
    }
    printf("\n"); 
  }
}  

Matrix Matrix::Multiply(Matrix &B, int type) {
  /*
    type: 1: C=A*B (default)
          2: C=A'*B
	  3: C=A*B'
	  4: C=A'*B'


     SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
           BETA, C, LDC)
     
     DGEMM performs: C = alpha*A*B + beta*C
     This routine assumes alpha=1 and beta=0;
     
     // Callling parameters:

     M = rows of product C
     N = cols of product C
     K = contraction dimension
     ALPHA = 1.0 here
     A = matrix A
     LDA = rows of A on entry
     B = matrix B
     LDB = rows of B on entry
     BETA = 0.0 here
     C = stores output matrix.
     LDC = rows of C on entry

  */
  
  // Check dimensions
  bool ok_dim = true;
  if (type==1 && cols != B.rows ) // C=A*B
    ok_dim = false;
  else if (type==2 && rows != B.rows ) // C=A'*B
    ok_dim = false;
  else if (type==3 && cols != B.cols ) // C=A*B'
    ok_dim = false;
  else if (type==4 && rows != B.cols ) // C=A'*B'
    ok_dim = false;
  
  if (!ok_dim) {
    printf("ERROR: Matrix::Multiply(): dimensions do not match for type %d\n",type);
    exit(1);
  }

  // Assign the dimensions.  dim1/dim2 are rows/cols for product matrix
  // and we contract over dimN.
  int dim1,dim2,dimN;
  //standard:
  dim1 = rows;
  dimN = cols; // = B.rows
  dim2 = B.cols;
  char TransA = 'N';
  char TransB = 'N';

  if (type==2) {
    dim1 = cols;
    dimN = rows;
    TransA = 'T';
  }
  else if (type==3) {
    dim2 = B.rows;
    TransB = 'T';
  }   
  else if (type==4) {
    dim1 = cols;
    dimN = rows;
    dim2 = B.rows;
    TransA = 'T';
    TransB = 'T';
  }

  double alpha = 1.0;
  double beta = 0.0;
  Matrix C(dim1,dim2);
  C.Set();

  // Call LAPACK matrix-matrix multiply
  dgemm_(&TransA, &TransB, &dim1, &dim2, &dimN, &alpha, matrix, &rows,
	 B.matrix, &B.rows, &beta, C.matrix, &dim1);

  return C;
}

// Compute eigenvalues for a real, symmetric matrix. The input matrix
// is preserved in this routine.
Vector Matrix::FindEigenvalues() {

  char job = 'N';
  char type = 'U';
  Vector evals(rows);
  int lwork = 3*rows-1;
  int info;
  double *work = new double[lwork];
  
  /*
    Lapack call information:
  
    DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    
    Parameters:
    JOBZ - 'N' = eigenvalues, 'V' = eigenvalues & eigenvectors
    UPLO - 'U' or 'L': Upper or lower triangular matrix supplied
    N    - size of the matrix
    A    - the matrix
    LDA  - number of rows of A
    W    - stores output eigenvalues, size N vector
    WORK - scratch space, size LWORK
    LWORK - 3*rows-1 in size;
    INFO - output code, 0 = success.  
  */

  dsyev_(&job, &type, &rows, matrix, &rows, evals.GetVectorPtr(), work, &lwork, &info); 

  if (info != 0) {
    printf("ERROR: Matrix::FindEigenvalues(): Diagonalization failed! DSYEV Info = %d\n",info);
    exit(1);
  }
   
  delete [] work;
  return evals;
}


// Find eigenvalues & eigenvectors for a real, symmetric matrix
// eigenvalues are returned in a vector, and each column of the
// original matrix is replaced with the corresponding eigenvectors.
// THE INPUT MATRIX IS DESTROYED HERE.
Vector Matrix::Diagonalize() {
  char job = 'V';
  char type = 'U';
  Vector evals(rows);
  int lwork = 3*rows-1;
  int info;
  double *work = new double[lwork];

  // Call LAPACK to do the real work
  dsyev_(&job, &type, &rows, matrix, &rows, evals.GetVectorPtr(), work, &lwork, &info); 

  if (info != 0) {
    printf("ERROR: Matrix::FindEigenvalues(): Diagonalization failed! DSYEV Info = %d\n",info);
    exit(1);
  }
  delete [] work;
  return evals;
}


// Perform singular value decomposition, A = U*SIGMA*V'.  
// Returns vector of singular values, and also U and V(transpose).
Vector Matrix::SVD(Matrix &U, Matrix &VT, bool find_condition) {

  char job = 'A';
  int dim = min(rows,cols);
  Vector sigma(dim);
  int info;
  // allocate temporary scratch space
  int tmp1 = 3*min(rows,cols) + max(rows,cols);
  int tmp2 = 5*min(rows,cols);
  int lwork = max(tmp1,tmp2);
  double *work = new double[lwork];

  // Verify size of U/VT, or allocate memory, as needed
  if (!U.allocated_) {
    U.Initialize(rows,rows);
  }
  else {
    if (U.rows != rows || U.cols != rows) {
      printf("Matrix::SVD - Input U matrix needs to be %dx%d or unallocated\n",
	     rows,rows);
      exit(1);
    }
  }
  if (!VT.allocated_) {
    VT.Initialize(cols,cols);
  }
  else {
    if (VT.rows != cols || VT.cols != cols) {
      printf("Matrix::SVD - Input VT matrix needs to be %dx%d or unallocated\n",
	     cols,cols);
      exit(1);
    }
  }
    
  // Call LAPACK to do the real work
  dgesvd_(&job, &job, &rows, &cols, matrix, &rows, sigma.GetVectorPtr(),
	  U.matrix, &rows, VT.matrix, &cols, work, &lwork, &info);

  if (info != 0) {
    printf("ERROR: Matrix::SVD(): Singular Value Decomposition failed! DGESVD Info = %d\n",info);
    exit(1);
  }
  

  if ( find_condition ) {
    double cond = sigma[0]/sigma[dim-1];
    printf("Condition number = %f\n",cond);
  }

  delete [] work;
  return sigma;
}

// Compute inverse of a square matrix, overwrites input matrix
void Matrix::Inverse() {

  /*
   Step 1: Perform LU factorization:
      DGETRF( M, N, A, LDA, IPIV, INFO )
   
   Step 2: Find Inverse
      DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  */

  if (rows != cols) {
    printf("ERROR: Matrix::Inverse() - matrix must be diagonalizable.  rows = %d, cols = %d\n",rows,cols);
    exit(1);
  }

  int N = rows;
  int info;
  double *ipiv = new double[N]; // pivot points from DGETRF, used in DGETRI
  double *work = new double[N]; // scratch space

  // LU decomposition
  dgetrf_(&N, &N, matrix, &N, ipiv, &info);

  if (info != 0) {
    printf("ERROR: Matrix::Inverse() - LU decomposition failed.  DGETRF Info = %d\n",info);
    exit(1);
  }
    
  // Compute inverse
  dgetri_(&N, matrix, &N, ipiv, work, &N, &info);
  
  if (info != 0) {
    printf("ERROR: Matrix::Inverse() - Inversion failed.  DGETRI Info = %d\n",info);
    exit(1);
  }
  
  delete [] ipiv;
  delete [] work;

}

// Solve A*x=b using LU decomposition.  Assume a single
// right-hand-side vector b.  Returns the vector x.
void Matrix::SolveLinearEquations(Vector &b) {
  int info;
  int NRHS = 1; // only one right-hand side (single b vector)
  int N = rows; // Number of linear equations = # of rows of the matrix
  int *ipiv = new int[N]; // The pivot indices that define the
			  // permutation matrix P;

  // Check dimensions: This routine requires a square matrix A, and cols(A) = rows(b).
  if ( rows != cols) {
    printf("Matrix::SolveLinearEquations() requires a square matrix: %d x %d\n",rows,cols);
    exit(1);
  }
  if ( cols != b.GetLength() ) {
    printf("Matrix::SolveLinearEquations(): Dimensions of b vector (%d) incompatible with A (%d x %d)\n",b.GetLength(),rows,cols);
    exit(1);
  }
  
  // Call LAPACK to do the real work
  dgesv_(&N, &NRHS, matrix, &rows, ipiv, b.GetVectorPtr(), &rows, &info);

  if (info != 0) {
    printf("ERROR: Matrix::SolveLinearEquations(): Equation solver failed! DGESV Info = %d\n",info);
    exit(1);
  }

  delete [] ipiv;
}

// Solve A*x=b using LU decomposition, expert version.  Includes
// iterative refinement and error estimates.  Useful if the matrix is
// potentially nearly singular or has some other pathology.  Assume a
// single right-hand-side vector b.  Returns the vector x.  This
// version requires roughly double the RAM of the non-expert version.
void Matrix::ExpertSolveLinearEquations(Vector &b) {
  int info;
  char FACT= 'N'; // no equilibration, and A has not been factorized.
		  // Maybe want to play with FACT='E'?
  char TRANS = 'N'; // normal (no transpose) for A.
  char EQUED; 
  int N = rows; // Number of linear equations = # of rows of the matrix

  int NRHS = 1; // only one right-hand side (single b vector).  Note,
		// if we ever want to increase NRHS, be sure to look
		// at ferr/berr allocation below.

  int *IPIV = new int[N]; // The pivot indices that define the
			  // permutation matrix P;

  double *AF = new double[rows*cols];
  double *R = new double[N];
  double *C = new double[N];
  double *X = new double[N*NRHS];
  double *WORK = new double[4*N];
  double rcond,ferr,berr; // if NHRS > 1, ferr and berr need to be arrays
  int *IWORK = new int[N];


  // Check dimensions: This routine requires a square matrix A, and cols(A) = rows(b).
  if ( rows != cols) {
    printf("Matrix::SolveLinearEquations() requires a square matrix: %d x %d\n",rows,cols);
    exit(1);
  }
  if ( cols != b.GetLength() ) {
    printf("Matrix::SolveLinearEquations(): Dimensions of b vector (%d) incompatible with A (%d x %d)\n",b.GetLength(),rows,cols);
    exit(1);
  }

  // Call LAPACK to do the real work
  dgesvx_(&FACT, &TRANS, &N, &NRHS, matrix, &N, AF, &N, IPIV, &EQUED, R, C,
	 b.GetVectorPtr(), &N, X, &N, &rcond, &ferr, &berr, WORK, IWORK, &info);
  
  if (info != 0) {
    printf("ERROR: Matrix::SolveLinearEquations(): Equation solver failed! DGESV Info = %d\n",info);
    exit(1);
  }

  printf("Expert Linear Equation Solver:\n");
  printf("   Condition Number = %f\n",rcond);
  printf("     Backward Error = %f\n",ferr);
  printf("      Forward Error = %f\n",berr);


  // Copy data from X to b
  for (int i=0;i<N*NRHS;i++) {
    b[i] = X[i];
  }

  // Free up memory
  delete [] IPIV;
  delete [] AF;
  delete [] R;
  delete [] C;
  delete [] X;
  delete [] WORK;
  delete [] IWORK;
}


double Matrix::TwoNorm() {
  int dim = min(rows,cols);
  Matrix U, VT;
  Vector sigma(dim);
  sigma = SVD(U,VT);
  return sigma[0];
}

void Matrix::Scale(double coeff) const {
 
  for (int x=0; x < rows*cols; x++) 
      matrix[x] *= coeff;
}

void Matrix::Transpose() {

  int dim1 = rows;
  int dim2 = cols;

  // backup the data
  double *tmp = new double[rows*cols];
  memcpy(tmp,matrix,sizeof(double)*rows*cols);

  for (int i=0;i<dim1;i++)
    for (int j=0;j<dim2;j++) {
      matrix[i*dim2 + j] = tmp[j*dim1 + i];
    }

  rows = dim2;
  cols = dim1;
  
  delete [] tmp;

}

//Starting operator-overloading here-on

// '=' operator will overwrite old matrix entirely
Matrix& Matrix::operator=(const Matrix& other) {

  if (this!=&other) {
    rows = other.rows;
    cols = other.cols;
    if (allocated_) 
      delete [] matrix;
    
    matrix = new double[rows*cols];
    memcpy(matrix,other.matrix,sizeof(double)*rows*cols);
    allocated_ = true;
  }

  return *this;
}

//Overloading "+=" operator

Matrix& Matrix::operator+=(const Matrix& other) {
  if ((rows == other.rows) && (cols == other.cols)) {
    for (int x=0; x < rows*cols; x++)  
        matrix[x] += other.matrix[x];
  }
  else {
    printf("Matrix::+= Error: matrix size mismatch!\n");
    exit(1);
  }
  return *this;
}


//  overloading "-=" operator
Matrix&	Matrix::operator-=(const Matrix& other)	{
  if ((rows == other.rows) && (cols == other.cols)) {
    for (int x=0; x < rows*cols; x++)
        matrix[x] -= other.matrix[x];
  }
  else {
    printf("Matrix::-= Error: matrix size mismatch!\n");
    exit(1);
  }
  return *this;
}


//overloading "*=" for scalar multiplication
Matrix& Matrix::operator*=(double coeff) {
  
  for (int x=0; x < rows*cols; x++)
      matrix[x] *= coeff;
  return *this;
} 


// overloading "+" for matrix addition
Matrix Matrix::operator+(const Matrix& other) {
  if ((rows == other.rows) && (cols == other.cols)) {
   Matrix C(rows,cols);

    for (int x=0; x < rows; x++)
      for (int y=0; y < cols; y++)
        C(x,y) = Element(x,y) + other(x,y); 

    return C;
  }
  else {
    printf("Matrix::+ Error: matrix size mismatch!\n");
    exit(1);
  }
}

// overloading "-" for matrix addition
Matrix Matrix::operator-(const Matrix& other) {
  if ((rows == other.rows) && (cols == other.cols)) {
   Matrix C(rows,cols);

    for (int x=0; x < rows; x++)
      for (int y=0; y < cols; y++)
        C(x,y) = Element(x,y) - other(x,y); 

    return C;
  }
  else {
    printf("Matrix::+ Error: matrix size mismatch!\n");
    exit(1);
  }
}


//Overload the "*" operator for simple matrix-matrix multiply C = A*other
Matrix Matrix::operator*(Matrix& other) {

  if (cols == other.rows) {     
   Matrix C(rows,other.cols);

   /*
   for (int x=0; x < rows; x++) {
     for (int y=0; y < other.cols; y++) {
       for (int z=0; z < cols; z++) {                       
	 C(x,y) += Element(x,z)*other(z,y);
       }
     }
     }*/
   C = Multiply(other,1);

   return C;
  }
  
  else {
    printf("The two matrices have incorrect dimensionalities for matrix multiplication\n");
    exit(1);
  }
}

// vec = A*b
Vector Matrix::MatrixTimesVector(const Vector& b) {
  
  int dimN = b.GetLength();
  if (cols != dimN) {
    printf("ERROR: Matrix::MatrixTimesVector(): Incompatible dimensions\n");
    exit(1);
  }

  Vector prod(rows);
  prod.Set();

  for (int i=0;i<rows;i++)
    for (int j=0;j<dimN;j++) {
      prod[i] += Element(i,j)*b[j];
    }

  return prod;
}

Vector Matrix::GetColumnVector(int icol) {

  double *tmp = new double[rows];
  for (int i=0;i<rows;i++) {
    tmp[i] = Element(i,icol);
  }

  Vector vec(tmp,rows);

  delete [] tmp;

  return vec;
}

Vector Matrix::GetRowVector(int irow) {

  double *tmp = new double[cols];
  for (int i=0;i<cols;i++) {
    tmp[i] = Element(irow,i);
  }

  Vector vec(tmp,cols);

  delete [] tmp;

  return vec;
}

void Matrix::SetColumnVector(Vector vec, int icol) {

  for (int i=0;i<rows;i++) {
    Element(i,icol) = vec[i];
  }
}

void Matrix::SetRowVector(Vector vec, int irow) {

  for (int i=0;i<cols;i++) {
    Element(irow,i) = vec[i];
  }
}

// row_offset and col_offset specify the matrix element for the top
// left corner of the block matrix.  nrows and ncols specify the size
// of the block.
Matrix Matrix::GetBlock(int nrows, int ncols, int row_offset, int col_offset) {

  // First check sizing... is the block too large?
  if ( (nrows + row_offset > rows ) || (ncols + col_offset) > cols) {
    printf("ERROR: Matrix::GetBlock(): Matrix block of size %d x %d from (%d,%d) is too large\n",nrows,ncols,row_offset,col_offset);
    printf("       Your original matrix is only %d x %d.\n",rows,cols);
    exit(1);
  }

  Matrix block(nrows,ncols);
  for (int j=0; j<ncols;j++) {
    for (int i=0; i<nrows;i++) {
      block.Element(i,j) = Element(i+row_offset,j+row_offset);
    }
  }

  return block;
}


// row_offset and col_offset specify the matrix element for the top
// left corner of the block matrix
void Matrix::SetBlock(const Matrix& block, int row_offset, int col_offset) {

  // First check sizing... will the block fit?
  if ( (block.rows + row_offset > rows ) || (block.cols + col_offset) > cols) {
    printf("ERROR: Matrix::SetBlock(): Matrix block of size %d x %d will not fit at (%d,%d)\n",block.rows,block.cols,row_offset,col_offset);
    printf("       Your target matrix is only %d x %d.\n",rows,cols);
    exit(1);
  }


  // Now set elements in the target matrix using those in the block
  for (int j=0; j<block.cols;j++) {
    for (int i=0; i<block.rows;i++) {
      Element(i+row_offset,j+col_offset) = block(i,j);
    }
  }

}

// row_offset and col_offset specify the matrix element for the top
// left corner of the block matrix
void Matrix::AddBlockToMatrix(const Matrix& block, int row_offset, int col_offset) {

  // First check sizing... will the block fit?
  if ( (block.rows + row_offset > rows ) || (block.cols + col_offset) > cols) {
    printf("ERROR: Matrix::AddBlocktoMatrix(): Matrix block of size %d x %d will not fit at (%d,%d)\n",block.rows,block.cols,row_offset,col_offset);
    printf("       Your target matrix is only %d x %d.\n",rows,cols);
    exit(1);
  }


  // Now set elements in the target matrix using those in the block
  for (int j=0; j<block.cols;j++) {
    for (int i=0; i<block.rows;i++) {
      Element(i+row_offset,j+col_offset) += block(i,j);
    }
  }

}


void Matrix::ReorderRows(int i, int j) {

  Vector tmp(GetRowVector(i)); // back up row i
  SetRowVector(GetRowVector(j),i); // move row j to row i
  SetRowVector(tmp,j); // put old row i in row j

}

void Matrix::ReorderColumns(int i, int j) {

  Vector tmp(GetColumnVector(i)); // back up column i
  SetColumnVector(GetColumnVector(j),i); // move column j to column i
  SetColumnVector(tmp,j); // put old column i in column j

}
