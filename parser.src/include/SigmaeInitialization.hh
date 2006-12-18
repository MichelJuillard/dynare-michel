#ifndef _SIGMAEINITIALIZATION_HH
#define _SIGMAEINITIALIZATION_HH

using namespace std;

#include <string>
#include <vector>

#include "Statement.hh"

//! Stores a Sigma_e statement
class SigmaeStatement : public Statement
{
public:
  //! Matrix form (lower or upper triangular) enum
  enum matrix_form_type
    {
      eLower = 0,              //!< Lower triangular matrix
      eUpper = 1               //!< Upper triangular matrix
    };
  //! Type of a matrix row
  typedef vector<string> row_type;
  //! Type of a complete matrix
  typedef vector<row_type> matrix_type;

  //! An exception indicating that a matrix is neither upper triangular nor lower triangular
  class MatrixFormException
  {
  };
private:
  //! The matrix
  const matrix_type matrix;
  //! Matrix form (lower or upper)
  const matrix_form_type matrix_form;

  //! Returns the type (upper or lower triangular) of a given matrix
  /*! Throws an exception if it is neither upper triangular nor lower triangular */
  static matrix_form_type determineMatrixForm(const matrix_type &matrix) throw (MatrixFormException);

public :
  SigmaeStatement(const matrix_type &matrix_arg) throw (MatrixFormException);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

#endif
