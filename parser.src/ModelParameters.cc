/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the ModelParemeters class methodes.
*/

#include "ModelParameters.hh"

ModelParameters::ModelParameters() : eq_nbr(0),
                                     exo_nbr(0), var_exo_nbr(0),
                                     exo_det_nbr(0), var_exo_det_nbr(0),
                                     endo_nbr(0), var_endo_nbr(0),
                                     parameter_nbr(0), local_parameter_nbr(0),
                                     recur_nbr(0),
                                     max_lag(0), max_lead(0),
                                     max_endo_lag(0), max_endo_lead(0),
                                     max_exo_lag(0), max_exo_lead(0),
                                     max_exo_det_lag(0), max_exo_det_lead(0),
                                     max_recur_lag(0), max_recur_lead(0)
{
  // Empty
}
