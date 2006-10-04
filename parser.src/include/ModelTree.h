#ifndef MODELTREE_H
#define MODELTREE_H
//------------------------------------------------------------------------------
/*! \file
 \version 1.0
 \date 04/13/2003
 \par This file defines the ModelTree class.
*/
//------------------------------------------------------------------------------
#include <string>
#include <vector>
#include <list>
#include <stack>
#include <sstream>
#include <fstream>
//#include <ext/hash_map>
//------------------------------------------------------------------------------
#include "SymbolTable.h"
#include "OperatorTable.h"
#include "NumericalConstants.h"
#include "ModelTypes.h"
#include "DataTree.h"
/*New*/
#include "BlockTriangular.h"
/*EndNew*/
//------------------------------------------------------------------------------
/*!
 \class  ModelTree
 \brief  Stores a model equations and derivatives.
*/
class ModelTree : public DataTree
{
	private :
		/*! Stores ID of equations and their derivatives */
		std::vector<std::vector<DerivativeIndex> > mDerivativeIndex;

		/*! Stores ID of equations and their 1st ordre derivatives */
		//vector<DerivativeIndex> mDerivative1Index;
		/*! Stores ID of equations and their 2nd ordre derivatives */
		//vector<DerivativeIndex> mDerivative2Index;
		/*! Only used for C++ standalone version */
		std::vector<double> 			mResults;
		/*! Output for static model */
		std::ostringstream			StaticOutput;
		/*! Output for dynamic model */
		std::ostringstream			DynamicOutput;
		/*! Output for main file */
		static std::ostringstream	output;
		/*! Output file stream for static model */
		std::ofstream				mStaticModelFile;
		/*! Output file stream for dynamic model */
		std::ofstream				mDynamicModelFile;
		/*!
		   A token is writen as a temporary expression
		   if its cost is greater than min_cost
		*/
		int 					min_cost;
		/*! left and right parentheses can be (,[ or ),] */
		char 					lpar, rpar;
		/*! Name of parameter variables ("params" for C output, and M_.params for Matlab) */
		std::string 			param_name;

	private :
		/*! Computes argument derivative */
		inline NodeID 		DeriveArgument(NodeID iArg, Type iType, int iVarID);
		/*! Sets derivative adresse of token in modeltree */
		inline void 		setDerivativeAdress(NodeID iTokenID, NodeID iDerivative,int iVarID);
		/*! Gets output argument of terminal token */
		inline std::string 		getArgument(NodeID id, Type type, EquationType  iEquationType);
		/*! Gets expression of part of model tree */
		inline std::string		getExpression(NodeID StartID, EquationType  iEquationType, int iEquationID = -1);
		/*!
		  Return true if node is to be writen as temp result
		  depending on its time cost and number of times it is referenced in model tree
		*/
		inline bool		writeAsTemp(NodeID id);
		/*! Remove unreferenced tokens after a derivation */
		//void 		RemoveUnref(int iBeginID, int iEndID, int iOrder);
		/*! Decrementing reference couter of unreferenced tokens in last derivative */
		void 		DecrementUnref(int iBeginID, int iEndID, int iOrder);
		inline int optimize(NodeID id);
	public :
		/*! When Jacobian (vs endogenous) is writen this flag is set to true */
		bool		computeJacobian;
		/*! When Jacobian (vs endogenous and exogenous) is writen this flag is set to true */
		bool		computeJacobianExo;
		/*! When Hessian  is writen this flag is set to true */
		bool		computeHessian;
		/*! Constructor */
		ModelTree();
		/*! Destructor */
		~ModelTree();
		/*! Opens output M files (1st and 2nd derivatives) */
		void OpenMFiles(std::string iModelFileName1, std::string iModelFileName2 = "");
		/*! Opens output C files (1st and 2nd derivatives) */
		void OpenCFiles(std::string iModelFileName1, std::string iModelFileName2 = "");
		/*! Saves output string into output M files */
		void SaveMFiles();

		/*! Saves output string into output C files */
		/*New*/
        void SaveCFiles(Model_Block* ModelBlock);
        /*EndNew*/
		/*! Computes derivatives of ModelTree */
		void 		derive(int iOrder);
		/*!
		  Writes output file for static model :
			- equations
			- 1st order derivatives with respect to endogenous variables (without lags)
		*/
		std::string 		setStaticModel(void);
		/*!
		  Writes output file for dynamic stochastic model :
		  - equations
		  - 1st order and 2nd order derivatives with respect to endogenous, exogenous, exogenous_det (in specific order)
		*/
		std::string 		setDynamicModel(void);

        /*New*/
		std::string  setDynamicModel_New(Model_Block* ModelBlock, simple* Index_Var_IM);
		/*EndNew*/


		/*! Writes initialization of various Matlab variables */
		void 		ModelInitialization(void);
		/*! Returns string output for main file */
		static 		std::string get();
};
//------------------------------------------------------------------------------
#endif
