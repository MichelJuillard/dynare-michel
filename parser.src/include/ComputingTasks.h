#ifndef COMPUTINGTASKS_H
#define COMPUTINGTASKS_H
//------------------------------------------------------------------------------
/** \file 
 * \version 1.0
 * \date 12/16/2003
 * \par This file defines the ComputingTasks class .
 */
//------------------------------------------------------------------------------
#include <sstream>
//------------------------------------------------------------------------------
#include "TmpSymbolTable.h"
#include "SymbolTable.h"
//------------------------------------------------------------------------------
/*! 
 \class EstimationParams
 \brief EstimationParams
 */
struct EstimationParams{
	int type;
	std::string name;
	std::string init_val;
	std::string prior;
	std::string low_bound;
	std::string up_bound;
	std::string mean;
	std::string std;
	std::string p3;
	std::string p4;
	std::string jscale;
	
	EstimationParams()
	{
		clear();
	}
	void clear()
	{
		type = 0;
		name = "";
		init_val = "NaN";
		prior = "NaN";
		low_bound = "-Inf";
		up_bound = "Inf";
		mean = "NaN";
		std = "NaN";
		p3 = "NaN";
		p4 = "NaN";
		jscale = "NaN";
	}
};

/*!
 \class  ComputingTasks
 \brief  This class concerns following Dynare commands :
 steady, check, simul, stoch_simul, dynare_estimation,  calib, osr and olr 
*/
class ComputingTasks
{
	private :
		//! Output of this class
		std::ostringstream	*output;

	public :		
		//! Constructor
		ComputingTasks();
		//! Destructor
		~ComputingTasks();
		/*! Pointer to error function of parser class */                               
		void (* error) (const char* m);
		/*! Estimation parameters */
		EstimationParams	*EstimParams;
		/*!
			Sets output reference
			\param iOutput : reference to an ostringstream
		*/
		void 	setOutput(std::ostringstream* iOutput); 
		// Prints "steady;" to output
		void 	setSteady(void);						
		// Prints "check;" to output
		void 	setCheck(void);							
		// Prints "simul(oo_.dr);" to output
		void 	setSimul(void);	
		//! Prints tmp1 and "stoch_simul(var_list_);"  to output
		void 	setStochSimul(std::string tmp1);				
		/*! Sets an option by printing 
			option_.<iName> = <iValue>
		*/
		void 	setOption(std::string iName, std::string iValue);
		/*! Prints "dynare_estimation;" */
		void 	runEstimation(std::string);							
		/*! Prints some estimation initialisation */
		void 	setEstimationInit(void);
		//! Prints optimization options */
		void 	setOptimOptions(std::string str1, std::string str2, int task);
		/*! Prints estimated elements */	
		void 	setEstimatedElements(void);
		void 	setEstimationStandardError(void);
		void    set_trend_element(std::string, std::string);

		void 	setCalibInit(void);
		
		void 	setCalibVariance(void);
		
		void 	setCalibCovariance(void);
		
		void 	setCalibAutoCorrelation(void);
		
		void 	setCalib(void);
		// write "osr(var_list_,osr_params_,optim_weights_);" 
		void 	setOsr(std::string tmp1);				
		// writes "olr(var_list_,olr_inst_,obj_var_,optim_weights_);"
		void 	setOlr(std::string tmp1, std::string tmp2);	 													
		
		void 	setOptimWeightsInit(void);
		
		void 	setOptimWeights1(void);
		
		void 	setOptimWeights2(void);
		
		void	set(void);
		
		static 	std::string	get(void);
};
//------------------------------------------------------------------------------
#endif
