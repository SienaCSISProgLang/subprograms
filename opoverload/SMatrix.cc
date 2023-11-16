#include "UtilException.h"
#include "SMatrix.h"
#include "SVector.h"
#include "SVectorDouble.h"
#include "SVector3.h"

#include <assert.h>

#ifdef _SCOREC_NewCompiler
#include <iostream>
#endif

#ifdef _SCOREC_NewCompiler
namespace SCOREC_Util {
#endif

void invmat(double **a, double **y, int n);

// SMatrixValue ops
SMatrixValue::SMatrixValue()
{v=0; Rows = Cols = 0;}

SMatrixValue::SMatrixValue(const SMatrixValue& mat)
{
  init(mat.Rows,mat.Cols,mat.v);
}

SMatrixValue::SMatrixValue(int r, int c, double **val)
{
  init(r,c,val);
}

void SMatrixValue::init(int r, int c, double **val)
{
  int i,j;
  v = new double *[Rows=r];
  Cols = c;
  double *help = new double[Cols*Rows];
  double *help2 = &help[0];
  double *help3 = val[0];
  for (i=0; i!= Rows; i++){
    v[i] = &help[i*Cols];
    for (j=0; j != Cols; j++) *(help2++) = *(help3++);
  }
}


SMatrixValue::~SMatrixValue()
{
  if(v)
  {
    if(v[0])
      delete [] v[0];
    delete [] v;
  }
}

// SMatrix ops
SMatrix::SMatrix()
  : M(new SMatrixValue)
{ }

SMatrix::SMatrix(int si, int sj)
  : M(new SMatrixValue)
{ setSize(si,sj); }

SMatrix::SMatrix(int si, int sj, double value)
  : M(new SMatrixValue)
{
  setSize(si,sj);
  double **v = M->v;
  int k = 0;
  for(int i=0; i < si; i++){
    for(int j=0; j < sj; j++, k++)
      v[i][j] = value;
  }
}

SMatrix::SMatrix(int si, int sj, const double *values)
  : M(new SMatrixValue)
{
  setSize(si,sj);
  double **v = M->v;
  int k = 0;
  for(int i=0; i < si; i++){
    for(int j=0; j < sj; j++, k++)
      v[i][j] = values[k];
  }
}

SMatrix::SMatrix(int r, int c, double ** m)
  : M(new SMatrixValue)
{
  M->v=m;
  M->Rows=r;
  M->Cols=c;
}

SMatrix::SMatrix(int si, int sj, SMatrixType type)
  : M(new SMatrixValue)
{
  setSize(si,sj);
  double **v = M->v;
  switch(type){
  case Ident:
    for(int i=0; i < si; i++){
      for(int j=0; j < sj; j++){
	v[i][j] = (i==j ? 1.0 : 0.0);
      }
    }
    break;
  }
}


void SMatrix::setRow(int row, double *values)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }
  double *help = &(M->v[row][0]);
  for(int i=0; i < M->Cols; i++)
    *(help++) = values[i];

}
void SMatrix::setColumn(int column, double *values, int minRow, int maxRow)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }
  if (maxRow == 0) maxRow = M->Rows-1;
  int cols = M->Cols;
  double *help = &(M->v[minRow][column]);
  int stop = maxRow - minRow + 1;
  for(int i=0; i < stop; i++)
  {
    *help = values[i];
    help += cols;
  }

}


int SMatrix::numCols() const
{ return M->Cols; }

int SMatrix::numRows() const
{ return M->Rows; }


double * & SMatrix::operator [] (int n)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }
  return( M->v[n]);
}

void SMatrix::setSize(int r, int c)
{

  if(M->isShared()){
    M = new SMatrixValue();
  }
  double **v = M->v = new double *[M->Rows=r];
  M->Cols = c;
  double *help = new double[c*r];
  for (int i=0; i!= r; i++)
    v[i] = &help[i*c];
}


SMatrix SMatrix::operator * (const SMatrix &mat) const
{
  const int Rows = M->Rows;
  const int Cols = M->Cols;
  double **v = M->v;
  double **mv = mat.M->v;

  int i,j,k,r,c;
  assert(Cols == mat.numRows());
  double ** newv = new double *[r=Rows];
  c = mat.numCols();
  double *help = new double[r*c];
  double *help2;
  double *help3;
  double *help4 = &help[0];

  for (i=0; i != r; i++)
  {
    newv[i] = &help[i*c];
    for (j=0; j != c; j++)
    {
      *help4 = 0.;
      help2 = v[i];
      help3 = &mv[0][j];
      for (k=0; k != Cols; k++)
      {
        *help4 += *(help2++) * *help3;
        help3+=c;
      }
      help4++;
    }
  }
  return SMatrix(r,c,newv);

}  


// Multiply M1 * M2 where the sizes are such that M2 can be replaced
// by the result, i.e. M1 has to be quadratic.
void SMatrix::multInPlace(SMatrix &mat)
{
  const int Cols = M->Cols;
  double **v = M->v;
  double **mv = mat.M->v;

  int i,j,k,c;
  const int Rows = M->Rows;
  assert(Rows == Cols);
  assert(Cols == mat.numRows());
  c = mat.numCols();
  double *help = new double[Cols];
  double *help2;
  double *help3, *help3a;

  for (i=0; i != c; i++)
  {
    // first save a column of mat since it will be overwritten
    help3a = help3 = &mv[0][i];
    for (int ii=0; ii != Cols; ii++)
    {
      help[ii] = *help3a;
      help3a += c;
    }
    for (j=0; j != Cols; j++)
    {
      help2 = &v[j][0];
      *help3 = 0.0;
      for (k=0; k != Cols; k++)
        *help3 += *(help2++) * help[k];
      help3 += c;
    }
  }
  delete help;
}

/*
This function multiplies B^T * C * B, where C is a 9x9 matrix and 
B is the full Gradient B Matrix as it returned by
DofMatrix Gradient(const Interpolation3d<DofVector> &interp, const SPoint3 &pt)
(see documentation for the function Gradient())

We use a special multiplication routine since the matrix B is stored sparse
without any zero values. This multiplication routine accounts for that.
*/

// this routine really shouldn't be in this library, mwb.

/*
SMatrix SMatrix::multFullGradBCB3D(SMatrix &B)
{
  double **v = M->v;
  //const int Rows = M->Rows;
  //const int Cols = M->Cols;
  int i,k;
  int r3 = B.numCols();
  int r = r3*3;
  double **newv = new double *[r]; // This will store pointer to help
  double *help = new double[r*r];  // This will be the final matrix
  double *help4 = help;
  double help_local[9];

  double *v1,*v2,*v3,b1,b2,b3,*B1,*B2,*B3;
  // Compute B^T*C line by line and multiply each line with B to get
  // the final matrix -> avoid to allocate memory for intermediate matrix
  for (i=0; i<r3; i++)
  {
    v1 = &v[0][0];
    v2 = &v[3][0];
    v3 = &v[4][0];
    b1 = B[0][i];
    b2 = B[3][i];
    b3 = B[4][i];
    for (k=0; k<9; k++)
      help_local[k] = b1 * *(v1++) + b2 * *(v2++) + b3 * *(v3++);
    // Line of intermediate matrix is done, now multiply with B
    B1 = &B[0][0];
    B2 = &B[3][0];
    B3 = &B[4][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[0]
                                    + *(B2++) * help_local[3]
                                    + *(B3++) * help_local[4];
    B1 = &B[1][0];
    B2 = &B[5][0];
    B3 = &B[6][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[1]
                                    + *(B2++) * help_local[5]
                                    + *(B3++) * help_local[6];
    B1 = &B[2][0];
    B2 = &B[7][0];
    B3 = &B[8][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[2]
                                    + *(B2++) * help_local[7]
                                    + *(B3++) * help_local[8];
  }
  for (i=0; i<r3; i++)
  {
    v1 = &v[1][0];
    v2 = &v[5][0];
    v3 = &v[6][0];
    b1 = B[1][i];
    b2 = B[5][i];
    b3 = B[6][i];
    for (k=0; k<9; k++)
      help_local[k] = b1 * *(v1++) + b2 * *(v2++) + b3 * *(v3++);
    // Line of intermediate matrix is done, now multiply with B
    B1 = &B[0][0];
    B2 = &B[3][0];
    B3 = &B[4][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[0]
                                    + *(B2++) * help_local[3]
                                    + *(B3++) * help_local[4];
    B1 = &B[1][0];
    B2 = &B[5][0];
    B3 = &B[6][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[1]
                                    + *(B2++) * help_local[5]
                                    + *(B3++) * help_local[6];
    B1 = &B[2][0];
    B2 = &B[7][0];
    B3 = &B[8][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[2]
                                    + *(B2++) * help_local[7]
                                    + *(B3++) * help_local[8];
  }
  for (i=0; i<r3; i++)
  {
    v1 = &v[2][0];
    v2 = &v[7][0];
    v3 = &v[8][0];
    b1 = B[2][i];
    b2 = B[7][i];
    b3 = B[8][i];
    for (k=0; k<9; k++)
      help_local[k] = b1 * *(v1++) + b2 * *(v2++) + b3 * *(v3++);
    // Line of intermediate matrix is done, now multiply with B
    B1 = &B[0][0];
    B2 = &B[3][0];
    B3 = &B[4][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[0]
                                    + *(B2++) * help_local[3]
                                    + *(B3++) * help_local[4];
    B1 = &B[1][0];
    B2 = &B[5][0];
    B3 = &B[6][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[1]
                                    + *(B2++) * help_local[5]
                                    + *(B3++) * help_local[6];
    B1 = &B[2][0];
    B2 = &B[7][0];
    B3 = &B[8][0];
    for (k=0; k<r3; k++) *(help4++) = *(B1++) * help_local[2]
                                    + *(B2++) * help_local[7]
                                    + *(B3++) * help_local[8];
  }


  for (i=0; i<r; i++) newv[i] = &help[r*i];
  return SMatrix(r,r,newv);
}
*/

SMatrix SMatrix::multBCBgeometricTerm(const SMatrix &B) const
{
  int r = B.numCols();
  int i,j;

  double **newv = new double *[r];
  double *help = new double[r*r];
  double *help2 = help;
  for (int k=0; k<r; k++)
  {
    newv[k] = &help[k*r];
    for (int l=0; l<r; l++)
      *(help2++) = 0.;
  }

  double *v1 = &(M->v[0][0]);
  double V = *v1;
  double V1 = *(v1+1),V2 = *(v1+2),V3 = *(v1+3),V4 = *(v1+4);
  double V5 = *(v1+5),V6 = *(v1+6),V7 = *(v1+7),V8 = *(v1+8);
  double *Bj = B.M->v[0];
  double *Bi0, *Bi1, *Bi2, *Bi3, *Bi4, *Bi5,*Bi6,*Bi7,*Bi8;
  int r1 = r+1;
  for (i=0; i<r; i+=3)
  {
    Bi0 = Bj + i; Bi1 = Bi0 + r; Bi2 = Bi1 + r;
    Bi3 = Bi2 + r1; Bi4 = Bi3 + r; Bi5 = Bi4 + r;
    Bi6 = Bi5 + r1; Bi7 = Bi6 + r; Bi8 = Bi7 + r;
    for (j=i; j<r; j+=3)
    {
      help[i*r+j]
        =(*Bi0 * V  + *Bi1 * V3+ *Bi2 * V6)* *(Bj+j)
        +(*Bi0 * V1 + *Bi1 * V4+ *Bi2 * V7)* *(Bj+r+j)
        +(*Bi0 * V2 + *Bi1 * V5+ *Bi2 * V8)* *(Bj+r+r+j);
       help[(i+1)*r+j+1]
         =(*Bi3 * V  + *Bi4 * V3+ *Bi5 * V6)* *(Bj+3*r+j+1)
         +(*Bi3 * V1 + *Bi4 * V4+ *Bi5 * V7)* *(Bj+4*r+j+1)
         +(*Bi3 * V2 + *Bi4 * V5+ *Bi5 * V8)* *(Bj+5*r+j+1);
       help[(i+2)*r+j+2]
         =(*Bi6 * V  + *Bi7 * V3+ *Bi8 * V6)* *(Bj+6*r+j+2)
         +(*Bi6 * V1 + *Bi7 * V4+ *Bi8 * V7)* *(Bj+7*r+j+2)
         +(*Bi6 * V2 + *Bi7 * V5+ *Bi8 * V8)* *(Bj+8*r+j+2);

    }
  }
  // account for symmetry
  for (i=0; i<r; i+=3)
  {
    for (j=i; j<r; j+=3)
    {
      help[j*r+i] = help[i*r+j];
      help[(j+1)*r+i+1] = help[(i+1)*r+j+1];
      help[(j+2)*r+i+2] = help[(i+2)*r+j+2];

    }
  }
  return SMatrix(r,r,newv);
}

/*
This works for all B^T*C*B
*/
SMatrix SMatrix::multBCB(const SMatrix &B) const
{
  double **v = M->v;
  //const int Rows = M->Rows;
  //const int Cols = M->Cols;

  int i,j,k,l,m;
  double t;
  int r = B.numCols();
  int c = numCols();
  double **newv = new double *[r];
  double *help = new double[r*r];
  double *local_row = new double[c];
  double *help_local;
  double *help4 = help;
  double *help_v;
  for (k=0; k<r; k++)
  {
    help_local = local_row;
    newv[k] = &help[r*k];
    t = B.getValue(0,k);
    //t = B[0][k];
    help_v = v[0];
    for (i=0; i != c; i++) *(help_local++)=t* *(help_v++);
    for (j=1; j != c; j++)
    {
      help_local = local_row;
      //t = B[j][k];
      t = B.getValue(j,k);
      help_v = v[j];
      for (i=0; i != c; i++)
        *(help_local++)+=t* *(help_v++);
    }

    help_local = local_row;
    //double *help_B = B[0];
    const double *help_B = B.getLine(0);
    help4 = newv[k];
    t = *local_row;
    for (l=0; l != r; l++) *(help4++)=t* *(help_B++);
    for (m=1; m != c; m++)
    {
      t = *(++help_local);
      help4 = newv[k];
      for (l=0; l != r; l++)
        *(help4++)+=t* *(help_B++);
    }
  }
  delete local_row;
  return SMatrix(r,r,newv);
}

/*
This works for cymmetric C
B^T*C*B
*/
SMatrix SMatrix::multBsymCB(const SMatrix &B) const
{
  double **v = M->v;
  //const int Rows = M->Rows;
  //const int Cols = M->Cols;

  int i,j,k,l,m;
  double t;
  int r = B.numCols();
  int c = numCols();
  double **newv = new double *[r];
  double *help = new double[r*r];
  double *local_row = new double[c];
  double *help_local;
  double *help4 = help;
  double *help_v;
  for (k=0; k<r; k++)
  {
    help_local = local_row;
    newv[k] = &help[r*k];
    t = B.getValue(0,k);
    help_v = v[0];
    for (i=0; i != c; i++) *(help_local++)=t* *(help_v++);
    for (j=1; j != c; j++)
    {
      help_local = local_row;
      t = B.getValue(j,k);
      help_v = v[j];
      for (i=0; i != c; i++)
        *(help_local++)+=t* *(help_v++);
    }

    help_local = local_row;
    double *help_B = B.getRow(0)+k;
    help4 = newv[k]+k;
    t = *local_row;
    for (l=0; l != r-k; l++) *(help4++)=t* *(help_B++);
    for (m=1; m != c; m++)
    {
      t = *(++help_local);
      help4 = newv[k]+k;
      help_B+=k;
      for (l=0; l != r-k; l++) *(help4++)+= t* *(help_B++);
    }
// Now symmetrize
    double *help5;
    help4 = newv[k]+k+1;
    help5 = newv[k]+r+k;
    for (l=1; l != r-k; l++)
    {
      *help5 = (*help4++);
      help5+=r;
    }
  }
  delete local_row;
  return SMatrix(r,r,newv);
}

void SMatrix::invert()
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  SMatrixValue temp(*(this->M));
  invmat(temp.v,M->v,M->Rows);
}

void SMatrix::transpose()
{
  int i,j;

  if(M->isShared())
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
 
  int nr = M->Rows, nc = M->Cols;

  // We need a temporary to hold the complete matrix. I don't know
  // of an algorithm to transpose an arbitrary rectangular matrix in
  // place using only a double as a temporary.
  double *temp = new double[nr*nc];

  // copy the matrix data into temporary storage space
  double *help1 = temp;
  double *help = M->v[0];
  for (i=0; i<nr*nc; i++)
    *(help1++) = *(help++);

  // reassign the number of columns/rows
  M->Rows = nc; M->Cols = nr;

  // copy the matrix data transposed into the original matrix
  help = M->v[0];
  help1 = temp;
  for(i = 0; i < nr; i++)
    for(j = 0; j < nc; j++)
      // at this point we don't have a valid operator [][] anymore, i.e.
      // we have to do some pointer arithmetic to get it right.
      *(help+i+j*nr) = *(help1++);

  // reassign the pointer to the rows.
  // delete only v which is the pointer to the rows, don't delete v[0]
  // since that points to the contents of the matrix
  double **v = M->v;
  delete v;
  // need now pointer array for the rows (which are the columns of the
  // old matrix)
  v = new double *[nc];
  for (i=0; i<nc; i++)
    v[i] = &help[i*nr];

  // we deleted M->v, so assign the new value
  M->v = v;

  // clean up
  delete temp;
}

SMatrix inverse(const SMatrix &m)
{
  SMatrix ret(m.numRows(),m.numCols());
  SMatrixValue temp(*(m.M));
  invmat(temp.v,ret.M->v,m.numRows());
  return ret;
}

SMatrix tranpose(const SMatrix &m)
{
  SMatrix ret(m.numCols(),m.numRows());

  int nr = m.numRows(), nc = m.numCols();
  for(int i = 0; i < nr; i++){
    for(int j = 0; j < nc; j++)
      ret.M->v[j][i] = m.M->v[i][j];
  }
  return ret;
}

SMatrix inverseTranspose(const SMatrix &m)
  // this isn't really any more efficent than tranpose(inverse(m))
  // should be able to make more efficent at some point
{
  SMatrix inv(m.numRows(),m.numCols());
  invmat(m.M->v,inv.M->v,m.numRows());
  
  SMatrix ret(m.numCols(),m.numRows());
  int nr = m.numRows(), nc = m.numCols();
  for(int i = 0; i < nr; i++){
    for(int j = 0; j < nc; j++)
      ret.M->v[j][i] = inv.M->v[i][j];
  }
  return ret;

}

SMatrix operator + (const SMatrix &m1, const SMatrix &m2)
{
  return SMatrix(m1) += m2;
}

SMatrix operator - (const SMatrix &m1, const SMatrix &m2)
{
  return SMatrix(m1) -= m2;
}

SMatrix & SMatrix::operator = (double val)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  
  for(i=0; i != Rows; i++){
    for(j=0; j != Cols; j++)
      v[i][j] = val;
  }
  return *this;
}

// This function computes M = fact1 * M + fact2 * mat
// The result is stored in M, i.e. no memory is allocated
// unless the matrix was shared and a new copy has to be split off.
void SMatrix::multAddInPlace(double fact1, const SMatrix &mat, double fact2)
{
  //const int Rows = M->Rows;
  //const int Cols = M->Cols;

  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  //double **v = M->v;

  assert(numCols() == mat.numCols());
  assert(numRows() == mat.numRows());

  double *help = M->v[0];
  double *help2 = mat.M->v[0];
  int stop = M->Cols*M->Rows;
  for (int i=0; i != stop; i++) 
  {
    *help = *help * fact1 + *(help2++) * fact2;
    help++;
  }
}  

SMatrix & SMatrix::operator += (const SMatrix &mat)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i;

#ifdef _SCOREC_NewCompiler
  if(numCols() != mat.numCols() || numRows() != mat.numRows())
    throw InvalidOperation(__FILE__,__LINE__);
#endif
    //"error in SMatrix::operator +=\n";  //error
  double *help = v[0];
  double *help2 = mat.M->v[0];
  int stop = Cols*Rows;
  for (i=0; i != stop; i++) 
    *(help++) += *(help2++);
  return *this;
}  

SMatrix & SMatrix::operator -= (const SMatrix &mat)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i;

#ifdef _SCOREC_NewCompiler
  if(numCols() != mat.numCols() || numRows() != mat.numRows())
    throw InvalidOperation(__FILE__,__LINE__);
#endif
    //"error in SMatrix::operator +=\n";  //error
  double *help = v[0];
  double *help2 = mat.M->v[0];
  int stop = Cols*Rows;
  for (i=0; i != stop; i++) 
    *(help++) -= *(help2++);
  return *this;
}  

SMatrix & SMatrix::operator *=(double mult)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i;

  double *help = v[0];
  int stop = Cols*Rows;
  for (i=0; i<stop; i++) 
    *(help++) *= mult;
  return *this;
}  

SMatrix & SMatrix::operator /=(double mult)
{
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }

  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i;

  double *help = v[0];
  int stop = Cols*Rows;
  for (i=0; i<stop; i++) 
    *(help++) /= mult;
  return *this;
}  

SMatrix SMatrix::tmul (const SMatrix &mat) const
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j,k,c;
  assert(Rows == mat.numRows());
  double ** newv = new double *[Cols];
  c = mat.numCols();
  double *help = new double[c*Cols];
  double *help4 = help;
  double *help3;
  double *helpvStart1 = &v[0][0];
  double *help3Start0 = &mat.M->v[0][0];
  for (i=0; i!= Cols; i++){
    newv[i] = &help[i*c];;
    double *help3Start = help3Start0;
    for (j=0; j != c; j++){
      help3 = help3Start++;
      *help4 = 0.0;
      double *helpv = helpvStart1;
      for(k=0; k != Rows; k++) 
      {
	*help4 += *helpv * *help3;
        helpv+=Cols; 
        help3+=c;
      }
      help4++;
    }
    helpvStart1++;
  }
  return SMatrix(Cols,c,newv);
  
}  

SMatrix SMatrix::mulTranspose(const SMatrix &mat) const
     // a*bT
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j,k,r,c;
  assert(Cols == mat.numCols());
  double ** newv = new double *[r=Rows];
  c = mat.numRows();
  double *help = new double[r*c];
  for (i=0; i!= r; i++){
    newv[i] = &help[i*c];;
    for (j=0; j != c; j++){
      newv[i][j]=0.0;
      for(k=0; k != Cols; k++) 
	newv[i][j] += v[i][k] * mat.M->v[j][k];
    }
  }
  return SMatrix(r,c,newv);
}  


SMatrix operator * (const SMatrix &m, double num)
{
  int i,j,r,c;

  c = m.numCols(); r = m.numRows();
  double ** newv = new double *[r];
  double *help = new double[r*c];
  double *help2 = help;
  double *help3 = m.M->v[0];
  for (i=0; i!= r; i++){
    newv[i] = &help[i*c];
    for (j=0; j != c; j++) *(help2++)= *(help3++)*num;
  }
  return SMatrix(r,c,newv);
}  

// The sgi wants this one even though we already have the function
// taking in a const SMatrix
SMatrix operator * (SMatrix &m, double num)
{
  int i,j,r,c;

  c = m.numCols(); r = m.numRows();
  double ** newv = new double *[r];
  double *help = new double[r*c];
  double *help2 = help;
  double *help3 = m.M->v[0];
  for (i=0; i!= r; i++){
    newv[i] = &help[i*c];
    for (j=0; j != c; j++) *(help2++)= *(help3++)*num;
  }
  return SMatrix(r,c,newv);
}  

SMatrix operator / (const SMatrix &m, double num)
{
  int i,j,r,c;

  c = m.numCols(); r = m.numRows();
  double ** newv = new double *[r];
  double *help = new double[r*c];
  double *help2 = help;
  double *help3 = m.M->v[0];
  for (i=0; i!= r; i++){
    newv[i] = &help[i*c];
    for (j=0; j != c; j++) 
      *(help2++)= *(help3++)/num;
  }
  return SMatrix(r,c,newv);
}  


SMatrix operator * (double num, const SMatrix &m)
{
  int i,j,r,c;

  c = m.numCols(); r = m.numRows();
  double ** newv = new double *[r];
  double *help = new double[r*c];
  double *help2 = help;
  double *help3 = m.M->v[0];
  for (i=0; i!= r; i++){
    newv[i] = &help[i*c];
    for (j=0; j != c; j++) *(help2++)= *(help3++)*num;
  }
  return SMatrix(r,c,newv);
}  

SVectorDouble SMatrix::operator * (const SVectorDouble &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  if(Cols != vec.size())
    ;  //error
  double * newv = new double[Rows];
  for(i=0;i != Rows; i++){
    newv[i]=0.0;
    for(j=0; j != Cols; j++)
      newv[i] += v[i][j]*vec.v[j];
  }
  return SVectorDouble(Rows,newv,1);  // SVector now owns newv
}


SVector<double> SMatrix::operator * (const SVector<double> &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Cols == vec.size());
  double * newv = new double[Rows];
  for(i=0;i != Rows; i++){
    newv[i]=0.0;
    for(j=0; j != Cols; j++)
      newv[i] += v[i][j]*vec.v[j];
  }
  return SVector<double>(Rows,newv,1);  // SVector now owns newv
}

SVector<double> SMatrix::tmul(const SVector<double> &vec)
{
  double **v = M->v;
  double *Vecv = vec.v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Rows == vec.size());
  double * newv = new double[Cols];
  for(i=0;i != Cols; i++){
    newv[i]=0.0;
    double *VecPointer = Vecv;
    double *helpv = &v[0][i];
    for(j=0; j != Rows; j++)
    {
      newv[i] += *helpv * *(VecPointer++);
      helpv += Cols;
    }
  }
  return SVector<double>(Cols,newv,1);  // SVector now owns newv

}


SVectorDouble SMatrix::tmul(const SVectorDouble &vec)
{
  double **v = M->v;
  double *Vecv = vec.v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Rows == vec.size());
  double * newv = new double[Cols];
  for(i=0;i != Cols; i++){
    newv[i]=0.0;
    double *VecPointer = Vecv;
    double *helpv = &v[0][i];
    for(j=0; j != Rows; j++)
    {
      newv[i] += *helpv * *(VecPointer++);
      helpv += Cols;
    }
  }
  return SVectorDouble(Cols,newv,1);  // SVector now owns newv

}

/*
SVector<double> SMatrix::tmul(const FVector<6> &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Rows == 6);
  double * newv = new double[Cols];
  for(i=0;i != Cols; i++){
    newv[i]=0.0;
    for(j=0; j != Rows; j++)
      newv[i] += v[j][i]*vec(j);
  }
  return SVector<double>(Cols,newv,1);  // SVector now owns newv

}
*/


SVector<SVector3> SMatrix::operator * (const SVector<SVector3> &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  if(Cols != vec.size())
    ;  //error
  SVector3 * newv = new SVector3[Rows];
  for(i=0;i != Rows; i++){
    newv[i][0]=0.0; newv[i][1]=0.0; newv[i][2]=0.0;
    for(j=0; j != Cols; j++)
      newv[i] += v[i][j]*vec.v[j];
  }
  return SVector<SVector3>(Cols,newv,1);  // SVector now owns newv

}

#ifndef IBM

SVector<SVector<double> > SMatrix::operator * (const SVector<SVector<double> > &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Cols == vec.size());
  int vsize = vec(0).size();
  SVector<double> * newv = new SVector<double> [Rows];
  for(i=0;i != Rows; i++){
    newv[i].setSize(vsize);
    newv[i] = 0;
    for(j=0; j != Cols; j++)
      newv[i] += v[i][j]*vec.v[j];
  }
  return SVector<SVector<double> >(Cols,newv,1);  // SVector now owns newv
}

SVector<SVectorDouble> SMatrix::operator * (const SVector<SVectorDouble> &vec)
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;
  assert(Cols == vec.size());
  int vsize = vec(0).size();
  SVectorDouble * newv = new SVectorDouble [Rows];
  for(i=0;i != Rows; i++){
    newv[i].setSize(vsize);
    newv[i] = 0;
    for(j=0; j != Cols; j++)
      newv[i] += v[i][j]*vec.v[j];
  }
  return SVector<SVectorDouble>(Cols,newv,1);  // SVector now owns newv
}
#endif

ostream &
operator <<(ostream &out, const SMatrix &m)
{
  int i,j;

  for(i=0; i != m.numRows(); i++){
    if(i)
      out << "\n[ ";
    else
      out << "[[ ";
    for(j=0; j != m.numCols(); j++)
      out << m(i,j) << " ";
    out << "]";
  }
  out << "]";
  return out;
}

void SMatrix::print() const
{
  double **v = M->v;
  const int Rows = M->Rows;
  const int Cols = M->Cols;

  int i,j;

  for(i=0; i != Rows; i++){
    for(j=0; j != Cols; j++)
      cout << v[i][j] << " ";
    cout << "\n";
  }
  cout << "\n";
}

#ifdef _SCOREC_NewCompiler
} // end of namespace SCOREC_Util
#endif

