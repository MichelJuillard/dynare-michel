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
int ModelParameters::endo_nbr = 0;
int ModelParameters::var_endo_nbr = 0;
int ModelParameters::parameter_nbr = 0;
int ModelParameters::lagged_nbr = 0;
int ModelParameters::static_nbr = 0;
int ModelParameters::forward_nbr = 0;
int ModelParameters::both_nbr = 0;
int ModelParameters::recur_nbr = 0;
int ModelParameters::max_lag=INT_MIN;
int ModelParameters::max_lead=INT_MIN;
using namespace std;
//int ModelParameters::endo_min_lag=INT_MAX;
//int ModelParameters::endo_max_lag=INT_MIN;
//int ModelParameters::exo_min_lag=INT_MAX;
//int ModelParameters::exo_max_lag=INT_MIN;
//int ModelParameters::recur_min_lag=INT_MAX;
//int ModelParameters::recur_max_lag=INT_MIN;
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
