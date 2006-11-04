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
		/*! Constructor */
		ModelParameters();
		/*! Destructor */
		~ModelParameters();
		/*!  Number of equations*/
		static int eq_nbr;
		/*!  Number of declared Exogenous variables */
		static int exo_nbr;
		/*!  Number of Exogenous variables that apear in model equations*/
		static int var_exo_nbr;
		/*!  Number of deterministic Exogenous variables */
		static int exo_det_nbr;
		static int var_exo_det_nbr;
		/*!  Number of declared Endogenous variables */
		static int endo_nbr;
		/*!  Number of Endogenous variables that apear in model equations*/
		static int var_endo_nbr;
		/*!  Number of parameters */
		static int parameter_nbr;
		/*!  Number of local parameters */
		static int local_parameter_nbr;
		/*!  Number of recursive variables */
		static int recur_nbr;		
		
		static int max_lag;
		static int max_lead;
		static int max_endo_lag;
		static int max_endo_lead;
		static int max_exo_lag;
		static int max_exo_lead;
		static int max_exo_det_lag;
		static int max_exo_det_lead;
		static int max_recur_lag;
		static int max_recur_lead;
};
//------------------------------------------------------------------------------
#endif
