#ifndef DATATREEDEF_H
#define DATATREEDEF_H
#include <iostream>
#include <string>
#include <vector>
using namespace std;
//------------------------------------------------------------------------------
#include "VariableTable.h"
#include "d_tab.h"
#include "NumericalConstants.h"
#include "DataTree.h"
//------------------------------------------------------------------------------
inline int DataTree::PushToken(int iArg1,int iOpCode, int iArg2, Type iType1)
{
	//cout << "ModeltreeSize = " << mModelTree.size() << endl;
	
	//cout << "PushToken 3\n"; 
	MetaToken 	lToken(iArg1, iType1, iArg2, iOpCode);

	if (iOpCode != NoOpCode)
	{
		lToken.op_name = operator_table.str(iOpCode);
	}
	else
	{
		lToken.op_name = "";	 
	}
	lToken.reference_count.resize(current_order+1,0);
	//lToken.reference_count[current_order] = 1;
	//for (int i=0; i < current_order; i++)
	//	lToken.reference_count[current_order] += lToken.reference_count[i];
	mModelTree.push_back(lToken);
	//Updating refernce counters
	if(iType1 == eTempResult && iOpCode != EQUAL)
	{
		//cout << mModelTree.size()-1 << ":" << iArg1.id << ":" << mModelTree[iArg1.id].reference_count[current_order] << endl;
		int s = mModelTree[iArg1].reference_count.size();
		for (int i = s; i <= current_order; i++)
		{
			int rc = mModelTree[iArg1].reference_count[i-1];
			mModelTree[iArg1].reference_count.push_back(rc);
		}
		mModelTree[iArg1].reference_count[current_order]++;
		//cout << mModelTree.size()-1 << ":" << iArg1.id << ":" << mModelTree[iArg1.id].reference_count[current_order] << endl;
		
	}
	if(iArg2 != -1 && iOpCode != EQUAL)
	{
		//cout << mModelTree.size()-1 << ":" << iArg2.id << ":" << mModelTree[iArg2.id].reference_count[current_order] << endl;
		int s = mModelTree[iArg2].reference_count.size();
		for (int i = s; i <= current_order; i++)
		{
			int rc = mModelTree[iArg2].reference_count[i-1];
			mModelTree[iArg2].reference_count.push_back(rc);
		}
		mModelTree[iArg2].reference_count[current_order]++;
		//cout << mModelTree.size()-1 << ":" << iArg2.id << ":" << mModelTree[iArg2.id].reference_count[current_order] << endl;
	}
	int ID = mModelTree.size()-1;
	mIndexOfTokens.insert(make_pair(lToken,ID));
	//cout << "PushToken 3 : end\n";
	return ID;
	
}
//------------------------------------------------------------------------------
inline bool DataTree::Exist(MToken iToken)
{
	map<MToken,int>::iterator iter;
	iter =  mIndexOfTokens.find(iToken);
	//Testing if token exists 
	if (iter == mIndexOfTokens.end()) return false;
	else return true;	
}
inline int DataTree::AddTerminal(int iArg1, Type iType1)
{
	if (iType1 == eTempResult)
	{
		return iArg1;
	}
	else
	{
		MetaToken 	lToken(iArg1, iType1, -1, NoOpCode);
		if (Exist(lToken))
		{			
			int ID = mIndexOfTokens[lToken];
			return ID;
		}
		return PushToken(iArg1, NoOpCode, -1, iType1);
	}
}
inline int DataTree::AddTerminal(string iArgName, int iLag)
{
	int		id;
	Type	type;
	
	if (!SymbolTable::Exist(iArgName))
	{
		cout << "ModelTree::AddToken : Error : Unknown symbol: " << iArgName << endl;
		exit(-1);
	}
	type = SymbolTable::getType(iArgName);
	if (type != eUNDEF)
	{
		if (type == eEndogenous ||
		  type == eExogenousDet  ||
		  type == eExogenous  ||
		  type == eRecursiveVariable)
		  {
			id = VariableTable::getID(iArgName,iLag);
			if (id == -1)
			{
				cout << "ModelTree::AddToken : Error : Unknown variable " << iArgName << " from VariableTable\n";
				exit(-1);
			}
		  }
		else
			id = SymbolTable::getID(iArgName);

	}
	else
	{
		cout << "ModelTree::AddToken : Error : Unknown parameter: " << iArgName << endl;
		exit(-1);
	}
	return AddTerminal(id, type);
}
inline int DataTree::AddPlus(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, PLUS);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 && iArg2 != 0)
	{
		return PushToken(iArg1,PLUS,iArg2);
	}
	else if (iArg1 != 0)
	{
	  return iArg1;
	}
	else if (iArg2 != 0)
	{
	  return iArg2;
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddMinus(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, MINUS);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 && iArg2 != 0)
	{
		return PushToken(iArg1,MINUS, iArg2);
	}
	else if (iArg1 != 0)
	{
		return iArg1;
	}
	else if (iArg2 != 0)
	{
		return PushToken(iArg2, UMINUS);
	}
	else
	{
		return 0;
	}
}
inline int DataTree::AddUMinus(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, UMINUS);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 )
	{
		return PushToken(iArg1, UMINUS);
	}
	else
	{
	  return 0;
	}

}
inline int DataTree::AddTimes(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, TIMES);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 && iArg1 != 1 && iArg2 != 0 && iArg2 != 1)
	{
		return PushToken(iArg1, TIMES, iArg2);
	}
	else if (iArg1 != 0 && iArg1 != 1 && iArg2 == 1)
	{
		return iArg1;
	}
	else if (iArg2 != 0 && iArg2 != 1 && iArg1 == 1)
	{
		return iArg2;
	}
	else if (iArg2 == 1 && iArg1 == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}	
}
inline int DataTree::AddDivide(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, DIVIDE);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 && iArg2 != 0 && iArg2 != 1)
	{
	  return PushToken(iArg1,DIVIDE,iArg2);
	}
	else if (iArg2 == 1)
	{
	  return iArg1;
	}
	else if (iArg1 == 0 && iArg2 != 0)
	{
	  return 0;
	}
	else
	{
	  cout << "DIVIDE 0/0 non available\n";
	  exit(-1);
	}
}
inline int DataTree::AddPower(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, POWER);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0 && iArg2 != 0 && iArg2 != 1)
	{
	  return PushToken(iArg1,POWER,iArg2);
	}
	else if (iArg2 == 1)
	{
	  return iArg1;
	}
	else if (iArg2 == 0)
	{
	  return 1;
	}
	else
	{
	  return 0;
	}

}
inline int DataTree::AddExp(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, EXP);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}

	if (iArg1 != 0)
	{
	  return PushToken(iArg1, EXP);
	}
	else
	{
	  return 1;
	}
}
inline int DataTree::AddLog(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, LOG);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0 && iArg1 != 1)
	{
	  return PushToken(iArg1, LOG);
	}
	else if (iArg1 == 1)
	{
		return 0;
	}
	else
	{
	  cout << "log(0) isn't available\n";
	  exit(-1);
	}
}
inline int DataTree::AddLog10(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, LOG10);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0 && iArg1 != 1)
	{
	  return PushToken(iArg1, LOG);
	}
	else if (iArg1 == 1)
	{
		return 0;
	}			
	else
	{
	  cout << "log10(0) isn't available\n";
	  exit(-1);
	}
}
inline int DataTree::AddCos(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, COS);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, COS);
	}
	else
	{
	  return 1;
	}
}
inline int DataTree::AddSin(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, SIN);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, SIN);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddTan(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, TAN);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, TAN);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddACos(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ACOS);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	return PushToken(iArg1, COS);
}
inline int DataTree::AddASin(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ASIN);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, SIN);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddATan(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ATAN);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, ATAN);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddCosH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, COSH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, COSH);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddSinH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, SINH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, SINH);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddTanH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, TANH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, TANH);
	}
	  else
	{
	  return 0;
	}
}
inline int DataTree::AddACosH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ACOSH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	return PushToken(iArg1, ACOSH);
}
inline int DataTree::AddASinH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ASINH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, ASINH);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddATanH(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, ATANH);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, ATANH);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddSqRt(int iArg1)
{
	MetaToken 	lToken(iArg1, eTempResult, -1, SQRT);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	if (iArg1 != 0)
	{
	  return PushToken(iArg1, SQRT);
	}
	else
	{
	  return 0;
	}
}
inline int DataTree::AddEqual(int iArg1, int iArg2)
{
	MetaToken 	lToken(iArg1, eTempResult, iArg2, EQUAL);

	if (Exist(lToken))
	{			
		int ID = mIndexOfTokens[lToken];
		return ID;
	}
	return PushToken(iArg1,EQUAL,iArg2);
}
#endif
