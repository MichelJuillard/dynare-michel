#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

#include "SymbolTable.h"



int main(void)
{
	SymbolTable st;
	ModelParameters mparams;
	
	//Adding 3 different symbols with AddSymbolDeclar
	cout << "\nTrying  AddSymbolDeclar Symbol 1, RecursiveVariable";
	SymbolTable::AddSymbolDeclar("Symbol 1",eRecursiveVariable);
	SymbolTable::PrintSymbolTable();
	cout << "\nTrying  AddSymbolDeclar Symbol 2, ExogenousDet";
	SymbolTable::AddSymbolDeclar("Symbol 2",eExogenousDet);
	SymbolTable::PrintSymbolTable();
	cout << "\nTrying  AddSymbolDeclar Symbol 3, Endogenous";
	SymbolTable::AddSymbolDeclar("Symbol 3",eEndogenous);
	SymbolTable::PrintSymbolTable();
	
	// Adding an existing symbol,same type, with AddSymbolDeclar
	cout << "\nTrying  AddSymbolDeclar Symbol 3, Endogenous";
	SymbolTable::AddSymbolDeclar("Symbol 3", eEndogenous);
	SymbolTable::PrintSymbolTable();
	// Adding an existing symbol, type is different, with AddSymbolDeclar
	cout << "\nTrying  AddSymbolDeclar Symbol 3, Exogenous";
	SymbolTable::AddSymbolDeclar("Symbol 3", eExogenous);
	SymbolTable::PrintSymbolTable();
	
	// Adding a new symbol with AddSymbolInline
	cout << "\nTrying  AddSymbolInline Symbol 4, RecursiveVariable";
	SymbolTable::AddSymbolInline("Symbol 4",eRecursiveVariable);
	SymbolTable::PrintSymbolTable();
		
	// Adding an existing non referenced symbol,same type, with AddSymbolInline
	cout << "\nTrying  AddSymbolInline Symbol 2, ExogenousDet";
	SymbolTable::AddSymbolInline("Symbol 2",eExogenousDet);
	SymbolTable::PrintSymbolTable();
	// Adding an existing non referenced symbol, type is different, with AddSymbolInline
	cout << "\nTrying  AddSymbolInline Symbol 1, Endogenous";
	SymbolTable::AddSymbolInline("Symbol 1",eEndogenous);
	SymbolTable::PrintSymbolTable();
	
	// Adding an existing referenced symbol, same type, with AddSymbolInline
	cout << "\nTrying  AddSymbolInline Symbol 4, Exogenous";
	SymbolTable::AddSymbolInline("Symbol 4",eExogenous);
	SymbolTable::PrintSymbolTable();
}
