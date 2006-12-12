#include "SigmaeInitialization.hh"

SigmaeStatement::SigmaeStatement(const matrix_type &matrix_arg) throw (MatrixFormException) :
  matrix(matrix_arg),
  matrix_form(determineMatrixForm(matrix))
{
}

SigmaeStatement::matrix_form_type
SigmaeStatement::determineMatrixForm(const matrix_type &matrix) throw (MatrixFormException)
{
  unsigned int nbe;
  int inc;
  matrix_form_type type;
  // Checking if first or last row has one element.
  if (matrix.front().size() == 1)
    {
      inc = 1;
      nbe = 2;
      type = eLower;
    }
  else if (matrix.back().size() == 1)
    {
      inc = -1;
      nbe = matrix.front().size()-1;
      type = eUpper;
    }
  else
    throw MatrixFormException();

  // Checking if matrix is triangular (upper or lower):
  // each row has one element more or less than the previous one
  // and first or last one has one element.
  matrix_type::const_iterator ir;
  for (ir = matrix.begin(), ir++; ir != matrix.end(); ir++, nbe += inc )
    if (ir->size() != nbe)
      throw MatrixFormException();

  return type;
}

void
SigmaeStatement::writeOutput(ostream &output) const
{
  unsigned int ic, ic1;
  unsigned int ir, ir1;

  output << "M_.Sigma_e = [...\n";
  for (ir = 0; ir < matrix.size(); ir++)
    {
      output << "\t";
      for (ic = 0; ic < matrix.size(); ic++)
        {
          if (ic >= ir && matrix_form == eUpper)
            {
              ic1 = ic-ir;
              ir1 = ir;
            }
          else if (ic < ir && matrix_form == eUpper)
            {
              ic1 = ir-ic;
              ir1 = ic;
            }
          else if (ic > ir && matrix_form == eLower)
            {
              ic1 = ir;
              ir1 = ic;
            }
          else if (ic <= ir && matrix_form == eLower)
            {
              ic1 = ic;
              ir1 = ir;
            }

          output << matrix[ir1][ic1] << " ";
        }
      output << ";...\n";
    }
  output << "\t];...\n";
}
