/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralMatrix.h,v 1.3 2004/11/24 20:41:59 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include "Vector.h"
#include <vector>
using namespace std;
vector<int>nullVec(0);

class GeneralMatrix;

class ConstGeneralMatrix {
  friend class GeneralMatrix;
  protected:
    ConstVector data;
    int rows;
    int cols;
    int ld;
  public:
    ConstGeneralMatrix(const double* d, int m, int n)
      : data(d, m*n), rows(m), cols(n), ld(m) {}
    ConstGeneralMatrix(const GeneralMatrix& m);
    ConstGeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols);
    ConstGeneralMatrix(const ConstGeneralMatrix& m, int i, int j, int nrows, int ncols);
    virtual ~ConstGeneralMatrix() {}
    
    const double& get(int i, int j) const
      {return data[j*ld+i];}
    int numRows() const {return rows;}
    int numCols() const {return cols;}
    int getLD() const {return ld;}
    const double* base() const {return data.base();}
    const ConstVector& getData() const {return data;}
    
    double getNormInf() const;
    double getNorm1() const;
    /* x = scalar(a)*x + scalar(b)*this*d */
    void multVec(double a, Vector& x, double b, const ConstVector& d) const;
    /* x = scalar(a)*x + scalar(b)*this'*d */
    void multVecTrans(double a, Vector& x, double b, const ConstVector& d) const;
    /* x = x + this*d */
    void multaVec(Vector& x, const ConstVector& d) const
      {multVec(1.0, x, 1.0, d);}
    /* x = x + this'*d */
    void multaVecTrans(Vector& x, const ConstVector& d) const
      {multVecTrans(1.0, x, 1.0, d);}
    /* x = x - this*d */
    void multsVec(Vector& x, const ConstVector& d) const
      {multVec(1.0, x, -1.0, d);}
    /* x = x - this'*d */
    void multsVecTrans(Vector& x, const ConstVector& d) const
      {multVecTrans(1.0, x, -1.0, d);}
    /* m = inv(this)*m */
    void multInvLeft(GeneralMatrix& m) const;
    /* m = inv(this')*m */
    void multInvLeftTrans(GeneralMatrix& m) const;
    /* d = inv(this)*d */
    void multInvLeft(Vector& d) const;
    /* d = inv(this')*d */
    void multInvLeftTrans(Vector& d) const;
    
    bool isFinite() const;
    /** Returns true of the matrix is exactly zero. */
    bool isZero() const;
    virtual void print() const;
  protected:
    void multInvLeft(const char* trans, int mrows, int mcols, int mld, double* d) const;
  };


class GeneralMatrix {
  friend class ConstGeneralMatrix;
  protected:
    Vector data;
    int rows;
    int cols;
    int ld;
    GeneralMatrix * tmpGMp;
public:
    GeneralMatrix(int m, int n)
      : data(m*n), rows(m), cols(n), ld(m) {}
    GeneralMatrix(const double* d, int m, int n)
      : data(d, m*n), rows(m), cols(n), ld(m) {}
    GeneralMatrix(double* d, int m, int n)
      : data(d, m*n), rows(m), cols(n), ld(m) {}
    GeneralMatrix(const GeneralMatrix& m);
    GeneralMatrix(const ConstGeneralMatrix& m);
    GeneralMatrix(const GeneralMatrix&m, const char* dummy); // transpose
    GeneralMatrix(const ConstGeneralMatrix&m, const char* dummy); // transpose
    GeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols);
    GeneralMatrix(GeneralMatrix& m, int i, int j, int nrows, int ncols);
    /* this = a*b */
    GeneralMatrix(const GeneralMatrix& a, const GeneralMatrix& b);
    /* this = a*b' */
    GeneralMatrix(const GeneralMatrix& a, const GeneralMatrix& b, const char* dum);
    /* this = a'*b */
    GeneralMatrix(const GeneralMatrix& a, const char* dum, const GeneralMatrix& b);
    /* this = a'*b */
    GeneralMatrix(const GeneralMatrix& a, const char* dum1,
				  const GeneralMatrix& b, const char* dum2);

    /* generate new matrix b as subset or whole of matrix a but 
    reordered by vrows and vcols as Matlab b=a(vrows,vcols) 
    ignores non-positive elements and passing zero length vector is equivalent to 
    Matlab operator ":" = all elements of that dimension in its order */

    GeneralMatrix(const GeneralMatrix& a, const vector<int>& vrows, const vector<int>& vcols);
    GeneralMatrix(const ConstGeneralMatrix& a, const vector<int> &vrows, const vector<int> &vcols)
      {GeneralMatrix(GeneralMatrix(a), vrows, vcols);};

    virtual ~GeneralMatrix();
    const GeneralMatrix& operator=(const GeneralMatrix& m)
     { 
      if (this!=&m)
        {
        if (rows==m.rows && cols==m.cols && ld==m.ld && data.base()!=m.data.base() )
          copy(m);
        else
          {
          data=m.data; rows=m.rows; cols=m.cols; ld=m.ld;
          }
        }
      return *this;
      }
   
    /* emulate Matlab repmat: new matrix = multv*multh*this */
    GeneralMatrix& repmat(int multv, int multh);

    const double& get(int i, int j) const
      {return data[j*ld+i];}
    double& get(int i, int j)
      {return data[j*ld+i];}
    int numRows() const {return rows;}
    int numCols() const {return cols;}
    int getLD() const {return ld;}
    double* base() {return data.base();}
    const double* base() const {return data.base();}
    Vector& getData() {return data;}
    const Vector& getData() const {return data;}
    
    double getNormInf() const
      {return ConstGeneralMatrix(*this).getNormInf();}
    double getNorm1() const
      {return ConstGeneralMatrix(*this).getNorm1();}
    
    /* place matrix m to the position (i,j) */
    void place(const ConstGeneralMatrix& m, int i, int j);
    void place(const GeneralMatrix& m, int i, int j)
      {place(ConstGeneralMatrix(m), i, j);}
    
    /* this = a*b */
    void mult(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b);
    void mult(const GeneralMatrix& a, const GeneralMatrix& b)
      {mult(ConstGeneralMatrix(a), ConstGeneralMatrix(b));}
    
    /* this = this + scalar*a*b */
    void multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b,
      double mult=1.0);
    void multAndAdd(const GeneralMatrix& a, const GeneralMatrix& b,
      double mult=1.0)
      {multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), mult);}
    
    /* this = this + scalar*a*b' */
    void multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b,
      const char* dum, double mult=1.0);
    void multAndAdd(const GeneralMatrix& a, const GeneralMatrix& b,
      const char* dum, double mult=1.0)
      {multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), dum, mult);}
    
    /* this = this + scalar*a'*b */
    void multAndAdd(const ConstGeneralMatrix& a, const char* dum, const ConstGeneralMatrix& b,
      double mult=1.0);
    void multAndAdd(const GeneralMatrix& a, const char* dum, const GeneralMatrix& b,
      double mult=1.0)
      {multAndAdd(ConstGeneralMatrix(a), dum, ConstGeneralMatrix(b), mult);}
    
    /* this = this + scalar*a'*b' */
    void multAndAdd(const ConstGeneralMatrix& a, const char* dum1,
      const ConstGeneralMatrix& b, const char* dum2, double mult=1.0);
    void multAndAdd(const GeneralMatrix& a, const char* dum1,
      const GeneralMatrix& b, const char* dum2, double mult=1.0)
      {multAndAdd(ConstGeneralMatrix(a), dum1, ConstGeneralMatrix(b),dum2, mult);}
    
    /* this = this + scalar*a*a' */
    void addOuter(const ConstVector& a, double mult=1.0);
    void addOuter(const Vector& a, double mult=1.0)
      {addOuter(ConstVector(a), mult);}
    
    /* this = this * m */
    void multRight(const ConstGeneralMatrix& m);
    void multRight(const GeneralMatrix& m)
      {multRight(ConstGeneralMatrix(m));}
    
    /* this = m * this */
    void multLeft(const ConstGeneralMatrix& m);
    void multLeft(const GeneralMatrix& m)
      {multLeft(ConstGeneralMatrix(m));}
    
    /* this = this * m' */
    void multRightTrans(const ConstGeneralMatrix& m);
    void multRightTrans(const GeneralMatrix& m)
      {multRightTrans(ConstGeneralMatrix(m));}
    
    /* this = m' * this */
    void multLeftTrans(const ConstGeneralMatrix& m);
    void multLeftTrans(const GeneralMatrix& m)
      {multLeftTrans(ConstGeneralMatrix(m));}

    /* Matlab element product: this = this .*m */
    void multElements(const GeneralMatrix& m);
    void multElements(const ConstGeneralMatrix& m)
      {multElements(GeneralMatrix(m));}


    /* this = this * m^(-1)  */
    void multInvRight(GeneralMatrix&m);
    
    /* x = scalar(a)*x + scalar(b)*this*d */
    void multVec(double a, Vector& x, double b, const ConstVector& d) const;
//      {ConstGeneralMatrix(*this).multVec(a, x, b, d);}
    
    /* x = scalar(a)*x + scalar(b)*this'*d */
    void multVecTrans(double a, Vector& x, double b, const ConstVector& d) const
      {ConstGeneralMatrix(*this).multVecTrans(a, x, b, d);}
    
    /* x = x + this*d */
    void multaVec(Vector& x, const ConstVector& d) const
      {ConstGeneralMatrix(*this).multaVec(x, d);}
    
    /* x = x + this'*d */
    void multaVecTrans(Vector& x, const ConstVector& d) const
      {ConstGeneralMatrix(*this).multaVecTrans(x, d);}
    
    /* x = x - this*d */
    void multsVec(Vector& x, const ConstVector& d) const
      {ConstGeneralMatrix(*this).multsVec(x, d);}
    
    /* x = x - this'*d */
    void multsVecTrans(Vector& x, const ConstVector& d) const
      {ConstGeneralMatrix(*this).multsVecTrans(x, d);}
    
    /* this = zero */
    void zeros();
    
    /** this = unit (on main diagonal) */
    void unit();
    
    /* this = NaN */
    void nans();
    
    /* this = Inf */
    void infs();
    
    /* this = scalar*this */
    void mult(double a);
    
    /* this = this + scalar*m */
    void add(double a, const ConstGeneralMatrix& m);
    void add(double a, const GeneralMatrix& m)
      {add(a, ConstGeneralMatrix(m));}
    
    /* this = this + scalar*m' */
    void add(double a, const ConstGeneralMatrix& m, const char* dum);
    void add(double a, const GeneralMatrix& m, const char* dum)
      {add(a, ConstGeneralMatrix(m), dum);}
    
    /* Returns true if this and m matrices are different for more than tolerance tol */
    bool isDiff(const GeneralMatrix& m, const double tol)const;
    bool isDiffSym(const GeneralMatrix& m, const double tol)const;
    bool isDiffUpprTriang(const GeneralMatrix& m, const double tol=0.0)const
      {return isDiffSym(m, tol);}

    bool isFinite() const
      {return (ConstGeneralMatrix(*this)).isFinite();}
    
    bool isZero() const
      {return (ConstGeneralMatrix(*this)).isZero();}
    
    virtual void print() const
      {ConstGeneralMatrix(*this).print();}

    void copyColumns(const GeneralMatrix& m, int istart, int iend, int ito);
    void copyColumns(const ConstGeneralMatrix& m, int istart, int iend, int ito)
      {copyColumns(GeneralMatrix( m), istart, iend, ito);};


    /* emulates Matlab command A(a,b)=B(c,d) where a,b,c,d are Matlab index vectors starting from 1 or ":") */
    static void AssignByVectors(GeneralMatrix& a, const vector<int>& vToRows, const vector<int>& vToCols
      , const GeneralMatrix& b, const vector<int>& vrows, const vector<int>& vcols);
    static void AssignByVectors(GeneralMatrix& a, const vector<int>& vToRows, const vector<int>& vToCols
      , const ConstGeneralMatrix& b, const vector<int> &vrows, const vector<int> &vcols)
      {AssignByVectors(a, vToRows, vToCols, GeneralMatrix(b), vrows, vcols);};
    void AssignByVectors(const vector<int>& vToRows, const vector<int>& vToCols
      , const GeneralMatrix& b, const vector<int>& vrows, const vector<int>& vcols)
      {AssignByVectors( *this, vToRows, vToCols, b, vrows, vcols);};
    void AssignByVectors( const vector<int>& vToRows, const vector<int>& vToCols
      , const ConstGeneralMatrix& b, const vector<int> &vrows, const vector<int> &vcols)
      {AssignByVectors(*this, vToRows, vToCols, GeneralMatrix(b), vrows, vcols);};
    void AssignByVectors( const GeneralMatrix& b, const vector<int>& vrows, const vector<int>& vcols)
      {AssignByVectors( *this, nullVec, nullVec, b, vrows, vcols);};
    void AssignByVectors( const ConstGeneralMatrix& b, const vector<int> &vrows, const vector<int> &vcols)
      {AssignByVectors(*this, nullVec, nullVec, GeneralMatrix(b), vrows, vcols);};
    void AssignByVectors(const vector<int>& vToRows, const vector<int>& vToCols, const GeneralMatrix& b)
      {AssignByVectors( *this, vToRows, vToCols, b, nullVec, nullVec);};
    void AssignByVectors( const vector<int>& vToRows, const vector<int>& vToCols, const ConstGeneralMatrix& b)
      {AssignByVectors(*this, vToRows, vToCols, GeneralMatrix(b), nullVec, nullVec);};


private:
  void copy(const ConstGeneralMatrix& m, int ioff , int joff );
  void copy(const GeneralMatrix& m, int ioff , int joff )
    {copy(ConstGeneralMatrix(m), ioff, joff);}
  void copy(const ConstGeneralMatrix& m);
  void copy(const GeneralMatrix& m);
  void gemm(const char* transa, const ConstGeneralMatrix& a,
    const char* transb, const ConstGeneralMatrix& b,
    double alpha, double beta);
  void gemm(const char* transa, const GeneralMatrix& a,
    const char* transb, const GeneralMatrix& b,
    double alpha, double beta)
    {gemm(transa, ConstGeneralMatrix(a), transb, ConstGeneralMatrix(b),
    alpha, beta);}
  
  /* this = this * op(m) (without whole copy of this) */
  void gemm_partial_right(const char* trans, const ConstGeneralMatrix& m,
    double alpha, double beta);
  void gemm_partial_right(const char* trans, const GeneralMatrix& m,
    double alpha, double beta)
    {gemm_partial_right(trans, ConstGeneralMatrix(m), alpha, beta);}
  
  /* this = op(m) *this (without whole copy of this) */
  void gemm_partial_left(const char* trans, const ConstGeneralMatrix& m,
    double alpha, double beta);
  void gemm_partial_left(const char* trans, const GeneralMatrix& m,
    double alpha, double beta)
    {gemm_partial_left(trans, ConstGeneralMatrix(m), alpha, beta);}
  
  /* number of rows/columns for copy used in gemm_partial_* */
  static int md_length;
};





#endif /* GENERAL_MATRIX_H */


// Local Variables:
// mode:C++
// End:
