#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>

using namespace std;
#include "SigmaeInitialization.h"

int main(void)
{
	SigmaeInitialization siginit;
	
	siginit.AddExpression("00");
	siginit.EndOfRow();
	siginit.AddExpression("10");
	siginit.AddExpression("11");
	siginit.EndOfRow();
	siginit.AddExpression("20");
	siginit.AddExpression("21");
	siginit.AddExpression("22");
	siginit.EndOfRow();
	siginit.AddExpression("30");
	siginit.AddExpression("31");
	siginit.AddExpression("32");
	siginit.AddExpression("33");
	siginit.EndOfRow();
	
	cout << siginit.get();
}
