#ifndef _MODELTYPES_HH
#define _MODELTYPES_HH
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/13/2004
 \par This file defines types for ModelTree class.
*///------------------------------------------------------------------------------
#include <string>
#include <list>
//#include <ext/hash_map>
#include <map>

typedef class MetaToken* NodeID;

/*!
  \struct  Derivative
  \brief  Derivative structure 
*/
//------------------------------------------------------------------------------
struct Derivative
{
	/*! ID of variable in VariableTable with respect to which token is derived */
	int variable_id; 			
	/*! "Forward" address of starting token of derivative */
	int derivative_address; 	
	/*! Constructor */
	inline Derivative(int iVarID = -1, int iDerivAdr = -1, Type iDerivType = eUNDEF)
	{
		variable_id = iVarID;
		derivative_address = iDerivAdr;
	}
};
/*!
  \struct  DerivativeIndex
  \brief  Derivative index structure 
*/
//------------------------------------------------------------------------------
struct DerivativeIndex{
	/*!/ id of "equal token" in mModelTree */
	NodeID token_id;
	/*! id of equation being derived */
	int equation_id;
	/*! column number of Jacobian or Hessian matrix */
	int derivators; 
	/*! Constructor */
	inline DerivativeIndex(NodeID iToken_id, int iEquationID, int iDerivator)
	{
		token_id = iToken_id;
		equation_id = iEquationID;
		derivators = iDerivator;
	}
};
/*!
  \struct  MToken
  \brief  Basic token structure, essentially for computing hash keys
*/
//------------------------------------------------------------------------------
struct MToken
{
  public :
	/*! ID of first operand */
	NodeID 	id1;
	/*! Type of first operand */
	Type 	type1;
	/*! ID of second operand */
	NodeID 	id2;
	/*! Operator code */	
	int 	op_code;
	/*! Constructor */
	inline MToken(NodeID iId1 = NULL, Type iType1 = eUNDEF, NodeID iId2 = NULL, int iOpCode = -1)
	{ 
		id1 = iId1;
        type1 = iType1;
		id2 = iId2;
		op_code = iOpCode;

	}
	/*! Copy constructor */
	inline MToken(const MToken& t)
	{
		id1 = t.id1;
		id2 = t.id2;
		type1 = t.type1;
		op_code = t.op_code;
	}
	/*! Destructor */
	~MToken() {	}

    /*! Computes the hash key */
    inline std::string Key() const
    {
        char key[100];
        snprintf(key, 100, "%p %p %d %d", id1, id2, type1, op_code);
        // std::cout << key << std::endl;
        return std::string(key);
    }
};
//------------------------------------------------------------------------------
/*!
  \struct  MetaToken
  \brief  Meta token structure 
*/
struct MetaToken : public MToken
{
 public :
 	/*! Index of token starting from zero */
 	int 						idx;
	/*!
	 Number of times a given token is referenced?
	 Each element of vector refer to derivative 0, 1, 2, etc
	*/
	std::vector<int> 				reference_count;
	/*! Commulative computation time cost */
	int 							cost;
	/*!
		If the token can be a temporary result, this flag is set to 1.
		If it can not be a temporary result (operation with constants etc)
		it is set to 0.
		After writing corresponding expression for the first time it is set to -1.
	*/
	int							tmp_status;
	/*! Name of unkown functions */
	std::string 						func_name;
	/*! Address of first order partial derivative */
	std::map<int, NodeID, std::less<int> > 	d1;
	/*! Expression of token, exp[0] for static model
		exp[1] for dynamic model */
	std::string						exp[2];
	/*! Flag : operand 1 is followed or not */
	bool						followed1;
	/*! Flag : operand 2 is followed or not */
	bool						followed2;
	/*! Operator name */
	std::string  					op_name;
    /*! set to one when left node has been treated */
    int left_done;
    /*! set to one when right node has been treated */
    int right_done;
    /*! set to one if closing parenthesis after token */
    int close_parenthesis;

	/*! Cast MetaToken -> MToken */
	inline operator MToken() const					
	{
		MToken lToken(id1,type1,id2, op_code);
		return lToken;
	}	
	/*! Equal operator */
	inline MetaToken operator = (const MetaToken& t1)
	{
		id1 = t1.id1;
		type1 = t1.type1;
		id2 = t1.id2;
		op_code = t1.op_code;
		op_name = t1.op_name;
		reference_count = t1.reference_count;
		cost = t1.cost;
		func_name = t1.func_name;		
		d1 = t1.d1;
		followed1 = t1.followed1; 
		followed2 = t1.followed2;
		exp[0] = exp[0];
		exp[1] = exp[1];
		idx = t1.idx;
		tmp_status = t1.tmp_status;
		left_done = t1.left_done;
		right_done = t1.right_done;
		close_parenthesis = t1.close_parenthesis;
		//cout << "Creation by equality\n";
		return *this;
	}
	/*! Constructor */
	inline MetaToken(NodeID iId1 = NULL, Type iType1 = eUNDEF,NodeID iId2= NULL, int iOpCode = -1)
		: MToken(iId1, iType1,iId2,iOpCode)
	{
		cost = 0;
		followed1 = followed2 = false; 
		tmp_status = 0;
		left_done = 0;
		right_done = 0;
		close_parenthesis = 0;
	}
	/*! Copy constructor */
	inline MetaToken(const MetaToken& mt)
		: MToken(mt)
	{
		followed1 = mt.followed1;
		followed2 = mt.followed2;
		cost = mt.cost;
		reference_count = mt.reference_count;
		func_name = mt.func_name;
		d1 = mt.d1;
		op_name = mt.op_name;
		exp[0] = mt.exp[0];
		exp[1] = mt.exp[1];
		idx = mt.idx;
		tmp_status = mt.tmp_status;
		left_done = mt.left_done;
		right_done = mt.right_done;
	}
	/*! Destructor */
	~MetaToken()
	{
	}

    /*! Sets derivative with respect to iVarID equal to iDerivative */
    inline void setDerivativeAddress(NodeID iDerivative, int iVarID)
    {
        d1[iVarID] = iDerivative;
    }
};
//------------------------------------------------------------------------------
/*!
  \brief  Less operator define comparaison function for map class
*/

inline	bool operator < (const MToken &t1, const MToken &t2)
{
		if 	(t1.op_code != t2.op_code)
		{
			return (t1.op_code < t2.op_code);
		}
		if (t1.id1 != t2.id1)
		{
			return (t1.id1 < t2.id1);
		}
		if (t1.id2 != t2.id2)
		{
			return (t1.id2 < t2.id2);
		}
		if (t1.type1 != t2.type1)
			return (t1.type1 < t2.type1);		
		else return false;

}

//------------------------------------------------------------------------------
/*! Equation type enum */
enum EquationType
{
	eStaticEquations = 0, 	//!< Equation of static model
	eDynamicEquations =	1,	//!<  Equation of dynamic model
	eStaticDerivatives = 2,	//!< Derivative of static model
	eDynamicDerivatives = 3	//!< Derivative of dynamic model
};
//------------------------------------------------------------------------------
#endif
