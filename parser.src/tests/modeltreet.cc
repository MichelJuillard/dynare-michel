#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
using namespace std;

#include "d_tab.h"
#include "SymbolTable.h"
#include "VariableTable.h"
#include "ModelTree.h"
#include "NumericalConstants.h"


int main(void)
{
	SymbolTable 		st; 
	VariableTable 		vt;
	NumericalConstants 	nc;
	ModelTree			model;
	vector<int>			t(20);
	
	//Adding 2 different symbols with AddSymbolDeclar	
	SymbolTable::AddSymbolDeclar("c",eExogenous);
	//SymbolTable::PrintSymbolTable();
	SymbolTable::AddSymbolDeclar("k",eEndogenous);
	//SymbolTable::PrintSymbolTable();
	SymbolTable::AddSymbolDeclar("aa",eParameter);
	//SymbolTable::PrintSymbolTable();
	SymbolTable::AddSymbolDeclar("x",eExogenous);
	SymbolTable::AddSymbolDeclar("alph",eParameter);
	SymbolTable::AddSymbolDeclar("delt",eParameter);

	VariableTable::AddVariable("k",0);	
	VariableTable::AddVariable("x",0);
	VariableTable::AddVariable("c",0);

	//VariableTable::AddVariable("k",1);
	//VariableTable::AddVariable("y",0);


	t[0] = model.AddToken("aa"); 
	t[1] = model.AddToken("x");
	t[2] = model.AddToken("k");
	t[3] = model.AddToken(Argument(t[0], eTempResult), 
						  Argument(t[1], eTempResult), TIMES);
	t[4] = model.AddToken("alph");
	t[5] = model.AddToken(Argument(t[2], eTempResult), 
						  Argument(t[4], eTempResult), POWER);
	t[6] = model.AddToken(Argument(t[3], eTempResult), 
						  Argument(t[5], eTempResult), TIMES);					  
						  
	t[7] = model.AddToken("delt");
	//t[8] = model.AddToken("1");
	t[9] = model.AddToken(Argument(1, eTempResult), 
						  Argument(t[7], eTempResult), MINUS);
	t[10] = model.AddToken(Argument(t[9], eTempResult), 
						  Argument(t[2], eTempResult), TIMES);
						  
	t[11] = model.AddToken(Argument(t[2], eTempResult), 
						  UMINUS);

	t[12] = model.AddToken(Argument(t[11], eTempResult), 
						  Argument(t[6], eTempResult), PLUS);
	t[13] = model.AddToken(Argument(t[12], eTempResult), 
						  Argument(t[10], eTempResult), PLUS);	
						  
	t[14] = model.AddToken("c");				  					  
	t[15] = model.AddToken(Argument(t[14], eTempResult), 
						  Argument(t[13], eTempResult), EQUAL);						  

	model.derive(2);
	cout << model.getDynamicModel() << endl;	
}
