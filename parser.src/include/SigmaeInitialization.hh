#ifndef _SIGMAEINITIALIZATION_HH
#define _SIGMAEINITIALIZATION_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/13/2004
  \par This file defines the SigmaeInitialization class.
*/
//------------------------------------------------------------------------------
#include <string>
#include <sstream>
#include <vector>
//------------------------------------------------------------------------------
/*!
  \class  SigmaeInitialization
  \brief  Handles Sigma_e command
*/
class SigmaeInitialization
{

  /*! Matrix type enum */
  enum MatrixType
    {
      eLower = 0,                  //!< Lower matrix
      eUpper = 1                   //!< Upper matrix
    };
private :
  /*!  A row of Sigma_e */
  std::vector<std::string>        row;
  //!  The hole matrix Sigma_e */
  std::vector<std::vector<std::string> >    matrix;
  /*! Output of this class */
  std::ostringstream        *output;
  /*! Matrix type (eLower(lower) or eUpper) */
  MatrixType          type;

  /*! Check that the matrix is triangular or square */
  void  CheckMatrix(void);
  /*! Print matrix to output */
  void  SetMatrix(void);

public :
  /*! Constructor */
  SigmaeInitialization();
  /*! Destructor */
  ~SigmaeInitialization();
  /*! Pointer to error function of parser class */
  void (* error) (const char* m);
  /*!
    Set output reference
    \param iOutput : reference to an ostringstream
  */
  void  setOutput(std::ostringstream* iOutput);
  /*!
    Add an expression to current row
    \param expression : a string expression
  */
  void  AddExpression(std::string expression);
  /*! Add current row to matrix and clears the row */
  void  EndOfRow();
  /*! Check matrix and print it to output */
  void  set(void);
};
//------------------------------------------------------------------------------
#endif
