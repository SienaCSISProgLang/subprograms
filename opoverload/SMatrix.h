/* Copyright (c) 1997 by the Scientific Computation Research Center,
 *                           Rensselaer Polytechnic Institute,
 *                           Troy, NY 12180
 */

#ifndef H_SMatrix
#define H_SMatrix

#ifdef _SCOREC_NewCompiler
#include <iosfwd>
using namespace std;
#else
#include <iostream.h>
#endif

#include "RCObject.h"
#include "RCPtr.h"

#ifdef _SCOREC_NewCompiler
namespace SCOREC_Util {
#endif

class SVector3;
template<class T> class SVector;
class SVectorDouble;

enum SMatrixType {
  Ident
  };

class SMatrixValue : public RCObject {
//friend class SVector<double>;
public:
  SMatrixValue();
  SMatrixValue(const SMatrixValue&);
  SMatrixValue(int r, int c, double **val);
  ~SMatrixValue();

  friend class RCPtr<SMatrixValue>;
  friend class SMatrix;

  friend SMatrix inverse(const SMatrix &);
  friend SMatrix tranpose(const SMatrix &);
  friend SMatrix inverseTranspose(const SMatrix &);
  friend SMatrix operator * (const SMatrix &, double);
  // Need to explicitly define this one for the sgi
  friend SMatrix operator * (SMatrix &, double);
  friend SMatrix operator * (double, const SMatrix &);
  friend SMatrix operator / (const SMatrix &, double);

protected:
  void init(int r, int c, double **val);
  double **v;
  int Rows;
  int Cols;
};

/** A matrix of doubles */
class SMatrix{
//friend class SVector<double>;
//friend SMatrix multBsymCB(FMatrix<6,6> &, SMatrix &);
//friend SMatrix multBsymCB(FMatrix<3,3> &, SMatrix &);
//friend SMatrix multBsymCB(FMatrix<1,1> &, SMatrix &);
//friend SMatrix multBCBgeometricTerm(FMatrix<3,3> &C, SMatrix &B);
public:
  SMatrix();
  ///
  SMatrix(int si, int sj); 
  ///
  SMatrix(int si, int sj, double value); 
  ///
  SMatrix(int si, int sj, const double *values);
  ///
  SMatrix(int si, int sj, SMatrixType);

  // queries
  ///
  // This is the same as operator (), but since it has only read access
  // we don't have to split off another instance of the matrix
  double getValue(int, int) const; 
  ///
  int numCols() const;
  ///
  int numRows() const;

  /// assignment
  SMatrix & operator = (double );

  /// element assignment
  double & operator()(int i, int j);
  /// element access
  double operator()(int i,int j) const;

  /// invert in place
  void invert();     
  /// transpose in place
  void transpose();  

  ///
  void setColumn(int, double *, int minRow = 0, int maxRow = 0);
  void setRow(int, double *);
  ///
  friend SMatrix inverse(const SMatrix &);
  ///
  friend SMatrix tranpose(const SMatrix &);
  ///
  friend SMatrix inverseTranspose(const SMatrix &);

  ///
  friend SMatrix operator + (const SMatrix &, const SMatrix &);
  ///
  friend SMatrix operator - (const SMatrix &, const SMatrix &);

  ///
  SMatrix operator * (const SMatrix &) const;
  ///
  SMatrix & operator += (const SMatrix &);
  ///
  SMatrix & operator -= (const SMatrix &);

  ///
  // compute M*fact1 + mat*fact2 and store result in M
  void multAddInPlace(double fact1, const SMatrix &mat, double fact2);
  ///
  SMatrix tmul(const SMatrix &) const;
  void multInPlace(SMatrix &);
  ///
  SMatrix mulTranspose(const SMatrix &) const;
  SMatrix multBCB(const SMatrix &) const;
  SMatrix multBCBgeometricTerm(const SMatrix &) const;
  SMatrix multBsymCB(const SMatrix &) const;

  // SMatrix-SVector operators
  ///
  SVector<double> operator * (const SVector<double> &);
  SVectorDouble operator * (const SVectorDouble &);
  ///
  SVector<double> tmul(const SVector<double> &);
 // SVector<double> tmul(const FVector<6> &);
  SVectorDouble tmul(const SVectorDouble &);
  ///
  SVector<SVector3> operator * (const SVector<SVector3> &vec);
  SVector<SVectorDouble > operator * (const SVector<SVectorDouble > &vec);
#ifndef IBM
  SVector<SVector<double> > operator * (const SVector<SVector<double> > &vec);
#endif

  // SMatrix-scalar operators
  ///
  friend SMatrix operator * (const SMatrix &, double);
  // The sgi needs this one to be happy
  friend SMatrix operator * (SMatrix &, double);
  ///
  friend SMatrix operator * (double, const SMatrix &);
  ///
  friend SMatrix operator / (const SMatrix &, double);

  ///
  SMatrix & operator *= (double );
  ///
  SMatrix & operator /= (double );

  double * & operator [] (int n);
  // Same as before, but gives only read access, i.e. we don't have to
  // check for reference counting
  const double * getLine(int i) const;
  ///
  void setSize(int,int);

  void print() const; // debugging  
friend ostream & operator <<(ostream &, const SMatrix &);
  
private:
  double *getRow(int) const; // returns the address of a row
  RCPtr<SMatrixValue> M;
  SMatrix(int r, int c, double ** m);
};

inline double & SMatrix::operator()(int i, int j)
{ 
  if(M->isShared()){
    M = new SMatrixValue(M->Rows,M->Cols,M->v);
  }
  return M->v[i][j]; 
}

// This is almost the same as above, but it gives only read access
// so no need to split of another instance
inline double SMatrix::getValue(int i, int j) const
{ return M->v[i][j]; }

// Again, a possibility to get read access to a line that saves the check
// needed for the reference counting
inline const double * SMatrix::getLine(int i) const
{
  return (M->v[i]);
}

inline double SMatrix::operator()(int i,int j) const
{ return M->v[i][j]; }

// This does the same as operator [], but since it has only
// read access we don't have to split off a new instance if the
// matrix is shared.
inline double * SMatrix::getRow(int i) const
{ return M->v[i]; }

#ifdef _SCOREC_NewCompiler
} // end of namespace SCOREC_Util
#endif


#endif

