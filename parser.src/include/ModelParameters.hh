#ifndef _MODELPARAMETERS_HH
#define _MODELPARAMETERS_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/13/2004
  \par This file defines the ModelParameters class.
*/
//------------------------------------------------------------------------------
/*!
  \class  ModelParameters
  \brief  Stors model parameters
*/
class ModelParameters
{
public :
  //! Constructor
  ModelParameters();
  /*!  Number of equations*/
  int eq_nbr;
  /*!  Number of declared Exogenous variables */
  int exo_nbr;
  /*!  Number of Exogenous variables that apear in model equations*/
  int var_exo_nbr;
  /*!  Number of deterministic Exogenous variables */
  int exo_det_nbr;
  int var_exo_det_nbr;
  /*!  Number of declared Endogenous variables */
  int endo_nbr;
  /*!  Number of Endogenous variables that apear in model equations*/
  int var_endo_nbr;
  /*!  Number of parameters */
  int parameter_nbr;
  /*!  Number of local parameters */
  int local_parameter_nbr;
  /*!  Number of recursive variables */
  int recur_nbr;

  int max_lag;
  int max_lead;
  int max_endo_lag;
  int max_endo_lead;
  int max_exo_lag;
  int max_exo_lead;
  int max_exo_det_lag;
  int max_exo_det_lead;
  int max_recur_lag;
  int max_recur_lead;
};
//------------------------------------------------------------------------------
#endif
