/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the ModelParemeters class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include "ModelParameters.h"
#include "limits.h"
//------------------------------------------------------------------------------
int ModelParameters::model_nbr = 0;
int ModelParameters::eq_nbr = 0;
int ModelParameters::exo_nbr = 0;
int ModelParameters::var_exo_nbr = 0;
int ModelParameters::exo_det_nbr = 0;
int ModelParameters::var_exo_det_nbr = 0;
int ModelParameters::endo_nbr = 0;
int ModelParameters::var_endo_nbr = 0;
int ModelParameters::parameter_nbr = 0;
int ModelParameters::local_parameter_nbr = 0;
int ModelParameters::lagged_nbr = 0;
int ModelParameters::static_nbr = 0;
int ModelParameters::forward_nbr = 0;
int ModelParameters::both_nbr = 0;
int ModelParameters::recur_nbr = 0;
int ModelParameters::max_lag = 0;
int ModelParameters::max_lead = 0;
int ModelParameters::max_endo_lag = 0;
int ModelParameters::max_endo_lead = 0;
int ModelParameters::max_exo_lag = 0;
int ModelParameters::max_exo_lead = 0;
int ModelParameters::max_exo_det_lag = 0;
int ModelParameters::max_exo_det_lead = 0;
int ModelParameters::max_recur_lag = 0;
int ModelParameters::max_recur_lead = 0;
using namespace std;
//------------------------------------------------------------------------------
ModelParameters::ModelParameters()
{
	// Empty
}
//------------------------------------------------------------------------------
ModelParameters::~ModelParameters()
{
	// Empty
}
//------------------------------------------------------------------------------
