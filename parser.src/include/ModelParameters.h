#ifndef MODELPARAMETERS_H
#define MODELPARAMETERS_H
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
		/*! Number of sub-models */
		static int model_nbr;
		/*!  Number of equations*/
		static int eq_nbr;
		/*!  Number of declared Exogenous variables */
		static int exo_nbr;
		/*!  Number of Exogenous variables that apear in model equations*/
		static int var_exo_nbr;
		/*!  Number of deterministic Exogenous variables */
		static int exo_det_nbr;
		/*!  Number of declared Endogenous variables */
		static int endo_nbr;
		/*!  Number of Endogenous variables that apear in model equations*/
		static int var_endo_nbr;
		/*!  Number of parameters */
		static int parameter_nbr;
		/*!  Number of lagged variables */
		static int lagged_nbr;
		/*!  Number of static variables */
		static int static_nbr;
		/*!  Number of forward-looking variables */
		static int forward_nbr;
		/*!  Number of variables that are both lagged and forward-looking */
		static int both_nbr;
		/*!  Number of equations */
		static int equation_nbr;
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

		/*!  Minimum lag for endogenous variables */
		//static int endo_min_lag;	
		/*! Maximum lag  for endogenous variables*/
		//static int endo_max_lag;
		/*! Miniimum lag  for exogenous variables*/
		//static int exo_min_lag;
		/*! Maximum lag  for exogenous variables*/
		//static int exo_max_lag;
		/*! Minimum lag  for recusive variables*/
		//static int recur_min_lag;
		/*! Maximum lag  for recursive variables*/
		//static int recur_max_lag;
};
//------------------------------------------------------------------------------
#endif
