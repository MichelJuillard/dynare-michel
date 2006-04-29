/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the NumericalInitialization class methodes.
*/
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
#include "NumericalInitialization.h"
//------------------------------------------------------------------------------
//ostringstream	NumericalInitialization::output;	
//------------------------------------------------------------------------------
NumericalInitialization::NumericalInitialization()
{
	//Empty
}
//------------------------------------------------------------------------------
NumericalInitialization::~NumericalInitialization()
{
	//Empty
}
//------------------------------------------------------------------------------
void NumericalInitialization::setOutput(ostringstream* iOutput)
{
	output = iOutput;	
}
//------------------------------------------------------------------------------
void  NumericalInitialization::SetConstant (string name, string expression)
{	  	

	//Testing if symbol exists 
	if (!SymbolTable::Exist(name))
    {
      		string msg = "Unknown parameter: " + name;
      		(* error) (msg.c_str());
    }
    // Testing symbol type
  	if (SymbolTable::getType(name) != eParameter)
    {
      		string msg = "Non-parameter used as a parameter: " + name;
      		(* error) (msg.c_str());
    }
    // Writing expression
    *output << "M_.params( " << SymbolTable::getID(name)+1 << " ) = " << expression << ";\n";
    *output << name << " = M_.params( " << SymbolTable::getID(name)+1 << " );\n";
    // Deleting expression
    //TODO
		
} 
//------------------------------------------------------------------------------
void  NumericalInitialization::SetLocalConstant (string name, string expression)
{	  	
  
  //Testing if symbol exists 
  if (SymbolTable::Exist(name))
    {
      if (SymbolTable::getType(name) != eLocalParameter)
	{
	  string msg = "Using existing symbol " + name + 
	    " as local parameter name";
	  (* error) (msg.c_str());
	}
    }
  else
    {
      SymbolTable::AddSymbolDeclar(name,eLocalParameter,name);
    }

  // Writing expression
  *output << name << " = " << expression << ";\n";
} 
//------------------------------------------------------------------------------
void  NumericalInitialization::BeginInitval (void) 
{
	
	// Writing a Matlab comment
	*output << "%\n% INITVAL instructions \n%\n";
	// Writing initval block to set initial values for variables 
	*output << "options_.initval_file = 0;\nendval_=0;\n";
	if(ModelParameters::recur_nbr > 0)
	{
		*output << "recurs_ = zeros(" << ModelParameters::recur_nbr << ", 1);\n";
	}
}
//------------------------------------------------------------------------------
void  NumericalInitialization::SetInit (string name, string expression)
{
	
	
	//Testing if symbol exists 
	if (!SymbolTable::Exist(name))
    {
    	string msg = "Unknown parameter: " + name;
   		(* error) (msg.c_str());
    }
    Type 	type = SymbolTable::getType(name);
    int 	id = SymbolTable::getID(name);
    // Writing instrcuction that set initial or terminal value 
    // for a variable 
    if (type == eEndogenous)
    {
    	*output << "oo_.steady_state( " << id+1 << " ) = " << expression << ";\n";
    }
    else if (type == eExogenous)
    {
    	*output << "oo_.exo_steady_state( " << id+1 << " ) = " << expression << ";\n";
    }
    else if (type == eExogenousDet)
    {
    	*output << "oo_.exo_det_steady_state( " << id+1 << " ) = " << expression << ";\n";
    }
    // Testing if symbol is a variable (Exogenousous deterministic or recursive) 
    else if ( type != eRecursiveVariable ) 
    {
    	cout << "Error : Non-variable symbol used in INITVAL: " << name << endl;
	}
}
//------------------------------------------------------------------------------
void  NumericalInitialization::EndInitval(void) 
{
	*output << "oo_.y_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];\n";
  	*output << "if M_.exo_nbr > 0;\n";
  	*output << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];\n";
  	*output <<"end;\n";
  	*output << "if M_.exo_det_nbr > 0;\n";
  	*output << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];\n";
  	*output <<"end;\n";
}
//------------------------------------------------------------------------------
void  NumericalInitialization::BeginEndval (void) 
{
	// Writing a Matlab comment
	*output << "%\n% ENDVAL instructions\n%\n";
	// Writing endval block to set terminal values for variables 
	*output << "ys0_= oo_.steady_state;\nex0_ = oo_.exo_steady_state;\nrecurs0_ = recurs_;\nendval_ = 1;\n";

}
//------------------------------------------------------------------------------
void  NumericalInitialization::EndEndval (void) 
{
	
	*output << "oo_.y_simul = [oo_.y_simul oo_.steady_state*ones(1,M_.maximum_lead+M_.maximum_lead)];\n";
  	*output << "if M_.exo_nbr > 0;\n";
  	*output << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];\n";
  	*output <<"end;\n";
  	*output << "if M_.exo_det_nbr > 0;\n";
  	*output << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];\n";
  	*output <<"end;\n";
}
//------------------------------------------------------------------------------
void  NumericalInitialization::BeginHistval (void) 
{
	// Writing a Matlab comment
	*output << "%\n% HISTVAL instructions\n%\n";
	
}
//------------------------------------------------------------------------------
void  NumericalInitialization::SetHist (string name, int lag, string expression) 
{

	
	//Testing if symbol exists 
	if (!SymbolTable::Exist(name))
    {
    	string msg = "Unknown parameter: " + name;
      	(* error) (msg.c_str());
    }
    Type 	type = SymbolTable::getType(name);
    int 	id = SymbolTable::getID(name);
    // Testing symbol type
	if (type == eEndogenous)
    {
    	*output << "oo_.y_simul( " << id+1 << ", MP.max_lag + " << lag + 1 << ") = " << expression << ";\n";
    }
    else if (type == eExogenous)
    {
    	*output << "oo_.exo_simul( MP.max_lag + " << lag + 1 << ", " << id+1 << " ) = " << expression << ";\n";
    }
    // Tetsting if symbol is a variable (Exogenousous deterministic or recursive) 
    else if ((type != eExogenousDet) || 
    	 (type != eRecursiveVariable))   		
    {
    	string msg = "Non-variable symbol : " + name;
   		(* error) (msg.c_str());
	}
    // Deleting expression
    // TODO
    
/////////////////////////////////
/*
  char buffer[200];
  int offset;

  offset = v->var_ptr-var_list;
  if (v->endo_exo == 1)
    {

      sprintf(buffer,"oo_.y_simul(%d,M_.maximum_lag+(%s))=",v->nbr+1,lag);
    }
  else if (v->endo_exo == 0)
    {
      initval_check1[offset] = 1;
      sprintf(buffer,"oo_.exo_simul(M_.maximum_lag+(%s),%d)=",lag,v->nbr+1);
    }

  str_output(buffer);

  p_expression(q);
  str_output(";\n");
*/
}
//------------------------------------------------------------------------------
/*
string NumericalInitialization::get(void)
{
	return output.str();
}
*/
//------------------------------------------------------------------------------
