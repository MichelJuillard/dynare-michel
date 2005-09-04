/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the SymbolTable class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
//------------------------------------------------------------------------------
#include "SymbolTable.h"

using namespace std;
//------------------------------------------------------------------------------
map<string, Symbol, less<string> > SymbolTable::symboltable;// = *(new symbolmap);
vector< vector<string> > SymbolTable::name_table = *(new vector< vector<string> >);	
vector< vector<string> > SymbolTable::tex_name_table = *(new vector< vector<string> >);	
ostringstream	SymbolTable::output;	
void (* SymbolTable::error) (const char* ) = NULL;
//------------------------------------------------------------------------------
SymbolTable::SymbolTable()
{
	  name_table.resize(10);
	  tex_name_table.resize(10);
}
//------------------------------------------------------------------------------
SymbolTable::~SymbolTable()
{
	  // Empty
}
//------------------------------------------------------------------------------
int SymbolTable::AddSymbol(string name,Type type, string tex_name)
{
  symboltable[name].type = type;
  symboltable[name].referenced = eNotReferenced; 
  name_table[(int) type].push_back(name);  
  tex_name_table[(int) type].push_back(tex_name);  
  
  switch (type)
    {
    case eExogenous:
	symboltable[name].id = ModelParameters::exo_nbr;
	ModelParameters::exo_nbr++;   
	return  ModelParameters::exo_nbr-1;
    case eExogenousDet:
	symboltable[name].id = ModelParameters::exo_det_nbr;
	ModelParameters::exo_det_nbr++;
	return  ModelParameters::exo_det_nbr-1;
    case eEndogenous:
	symboltable[name].id = ModelParameters::endo_nbr;
	ModelParameters::endo_nbr++;
	return  ModelParameters::endo_nbr-1;
    case eParameter:
	symboltable[name].id = ModelParameters::parameter_nbr;
	ModelParameters::parameter_nbr++;
	return  ModelParameters::parameter_nbr-1;
    case eRecursiveVariable:
	symboltable[name].id = ModelParameters::recur_nbr;
	ModelParameters::recur_nbr++;
	return  ModelParameters::recur_nbr-1;
    default:
  // should never happen
      return -1;
    }
  
}
//------------------------------------------------------------------------------
int SymbolTable::AddSymbolDeclar(string name,Type type, string tex_name)
{
	//Testing if the symbol exist in the map 
	if ( !Exist(name) ){
		//The symbol dosn't exist, adding it
		return AddSymbol(name,type, tex_name);
	}
	else{   
		//The symbol exists, testing its type
		if (symboltable[name].type == type){
			cout << "Warning : symbol " << name << " declared more than once.\n";
			return getID(name);		}
		else{
			string msg = "symbol " + name + " declared more than once with different types.";
			(* error) (msg.c_str());
			return -1;
		}
	}

}
//------------------------------------------------------------------------------
void SymbolTable::AddSymbolRange(string name,int nbr,Type type, string tex_name)
{
}
//------------------------------------------------------------------------------
void SymbolTable::AddLag(string name,int lag)
{
	//Testing if the symbol exist in the map 
	if ( !Exist(name) ){
		return;
	}
	Type type = symboltable[name].type;
	if ((type == eEndogenous) || (type == eExogenous))
	{
		symboltable[name].lags.push_back(lag);
	}
}
//------------------------------------------------------------------------------
void  SymbolTable::ResetType(string name,Type new_type)
{
	symboltable[name].type = new_type;
}
//------------------------------------------------------------------------------
void  SymbolTable::SetReferenced(string name)
{
	symboltable[name].referenced = eReferenced;
}
//------------------------------------------------------------------------------
void SymbolTable::SetOutOfScoop(string name)
{
	symboltable[name].referenced = eOutOfScoop;
}
//------------------------------------------------------------------------------
Reference 	SymbolTable::isReferenced(std::string name)
{
	return symboltable[name].referenced;
}
//------------------------------------------------------------------------------
void SymbolTable::clean()
{
	string unused;
	bool   warning = false;
	vector<Type>	types(3);
	vector<int>		nb_type(3);
	
	types[0] = eEndogenous;
	types[1] = eExogenous;
	types[2] = eExogenousDet;
	
	nb_type[0] = ModelParameters::endo_nbr;
	nb_type[1] = ModelParameters::exo_nbr;
	nb_type[2] = ModelParameters::exo_det_nbr;
	
	// Removing unused variables
	for (int t = 0; t < 3; t++)
	{
		// Checking if all variables are used
		for (int s1 = 0; s1 < nb_type[t]; s1++)
		{
			string name = getNameByID(types[t],s1);
			string tex_name = getTexNameByID(types[t],s1);
			if (isReferenced(name) == eNotReferenced)
			{
				symboltable.erase(name);
				vector<string>::iterator it;
				it = find(name_table[types[t]].begin(), name_table[types[t]].end(), name);
				name_table[types[t]].erase(it);
				it = find(tex_name_table[types[t]].begin(), tex_name_table[types[t]].end(), tex_name);
				tex_name_table[types[t]].erase(it);
				nb_type[t]--;
				unused += "fprintf(1,'%-30s";
				switch(types[t])
				  {
				  case eEndogenous 	: unused += "Endogenous variable\\n','";break;
				  case eExogenous 	: unused += "Exogenous variable\\n','";break;
				  case eExogenousDet 	: unused += "Exogenous deterministic variable\\n','";break;
				  default : ;
				  }
				unused += name;
				unused += "');\n";
				warning = true;
				for (int s2 = s1; s2 < nb_type[t]; s2++)
				{
					name = getNameByID(types[t],s2);
					// Decrementing symbol table ids in ST
					symboltable[name].id--;
				}
				s1--;
			}
		}
	}
	ModelParameters::endo_nbr 		= nb_type[0];
	ModelParameters::exo_nbr 		= nb_type[1];
	ModelParameters::exo_det_nbr	= nb_type[2];
	/*
	// Checking if unused parameters
	for (int s1 = 0; s1 < ModelParameters::parameter_nbr; s1++)
	{
		string name = getNameByID(eParameter,s1);
		if (isReferenced(name) == eNotReferenced)
		{
				unused += "fprintf(1,'%-30sParameter\\n','";
				unused += name;
				unused += "');\n";

			warning = true;
		}
	}
*/
	if (warning)
	{
		output << "fprintf(1,'Warning : symbol(s) :\\n');\n";
		output << unused;
		output << "fprintf(1,'are declared but not used in the model equations. ');\n";
		output << "reply = input('Continue? [y]\\\\n ','s');\n";
		output << "if isempty(reply), reply='y'; end;\n";
		output << "if strcmpi(reply(1),'n'),\n  return;\nend\n";
	}
	//PrintSymbolTable();
}
string SymbolTable::get()
{
	if (ModelParameters::exo_nbr > 0)
	{
		output << "M_.exo_names = '" << getNameByID(eExogenous, 0) << "';\n";
		output << "M_.exo_names_tex = '" << getTexNameByID(eExogenous, 0) << "';\n";
		for (int id = 1; id < ModelParameters::exo_nbr; id++)
		{
			output << "M_.exo_names = str2mat(M_.exo_names,'" << getNameByID(eExogenous, id) << "');\n";
			output << "M_.exo_names_tex = str2mat(M_.exo_names_tex,'" << getTexNameByID(eExogenous, id) << "');\n";
		}
	}
	if (ModelParameters::exo_det_nbr > 0)
	{
		output << "lgxdet_ = '" << getNameByID(eExogenousDet, 0) << "';\n";
		output << "lgxdet_tex_ = '" << getTexNameByID(eExogenousDet, 0) << "';\n";
		for (int id = 1; id < ModelParameters::exo_det_nbr; id++)
		{
			output << "lgxdet_ = str2mat(lgxdet_,'" << getNameByID(eExogenousDet, id) << "');\n";
			output << "lgxdet_tex_ = str2mat(lgxdet_tex_,'" << getTexNameByID(eExogenousDet, id) << "');\n";
		}
	}
	if (ModelParameters::endo_nbr > 0)
	{
		output << "M_.endo_names = '" << getNameByID(eEndogenous, 0) << "';\n";
		output << "M_.endo_names_tex = '" << getTexNameByID(eEndogenous, 0) << "';\n";
		for (int id = 1; id < ModelParameters::endo_nbr; id++)
		{
			output << "M_.endo_names = str2mat(M_.endo_names,'" << getNameByID(eEndogenous, id) << "');\n";
			output << "M_.endo_names_tex = str2mat(M_.endo_names_tex,'" << getTexNameByID(eEndogenous, id) << "');\n";
		}
	}
	if (ModelParameters::recur_nbr > 0)
	{
		output << "M_.recur_names = '" << getNameByID(eRecursiveVariable, 0) << "';\n";
		output << "M_.recur_names_tex = '" << getTexNameByID(eRecursiveVariable, 0) << "';\n";
		for (int id = 1; id < ModelParameters::recur_nbr; id++)
		{
			output << "M_.recur_names = str2mat(M_.recur_names,'" << getNameByID(eRecursiveVariable, id) << "');\n";
			output << "M_.recur_names_tex = str2mat(M_.recur_names_tex,'" << getTexNameByID(eRecursiveVariable, id) << "');\n";
		}
	}
	if (ModelParameters::parameter_nbr > 0)
	{
		output << "M_.param_names = '" << getNameByID(eParameter, 0) << "';\n";
		output << "M_.param_names_tex = '" << getTexNameByID(eParameter, 0) << "';\n";
		for (int id = 1; id < ModelParameters::parameter_nbr; id++)
		{
			output << "M_.param_names = str2mat(M_.param_names,'" << getNameByID(eParameter, id) << "');\n";
			output << "M_.param_names_tex = str2mat(M_.param_names_tex,'" << getTexNameByID(eParameter, id) << "');\n";
		}
	}

	output << "M_.exo_det_nbr = " << ModelParameters::exo_det_nbr << ";\n";	
	output << "M_.exo_nbr = " << ModelParameters::exo_nbr << ";\n";
	output << "M_.Sigma_e = zeros(" << ModelParameters::exo_nbr
	       << ", " << ModelParameters::exo_nbr << ");\n";
	output << "M_.endo_nbr = " << ModelParameters::endo_nbr << ";\n";
	output << "M_.recur_nbr = " << ModelParameters::recur_nbr << ";\n";
	output << "M_.param_nbr = " << ModelParameters::parameter_nbr << ";\n";
	return output.str();
}


//------------------------------------------------------------------------------
/*
ostream & operator << (ostream & output, Symbol s)
{
	switch (s.type) {
    		case eExogenous    : output << "Exogenous    , "; break;
    		case eExogenousDet : output << "ExogenousDet ,"; break;
    		case eEndogenous   : output << "Endogenous   , "; break;
    		case eRecursiveVariable : output << "RecursiveVariable , "; break;
    		case eParameter : output << "Parameter , "; break;
     		case eLoopIndex : output << "LoopIndex , "; break; 
    	}
	output << s.id;
	switch (s.referenced) {
		case eReferenced    : cout << ", Referenced, "; break; 
		case eNotReferenced : cout << ", NotReferenced, "; break; 
		case eOutOfScoop    : cout << ", OutOfScoop, "; break; 

	}
	return output;
}
//------------------------------------------------------------------------------
void SymbolTable::PrintSymbolTable(void)
{
	map<string, Symbol, less<string> >::iterator iter;
	cout << "\nCurrent Map :\n";
	for (iter = symboltable.begin();iter != symboltable.end(); iter++)
		cout << (*iter).first << " , " << (*iter).second << endl;;
}
//------------------------------------------------------------------------------
*/
  
