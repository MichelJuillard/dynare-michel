#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
using namespace std;

int mexPrintf(/*const string*/const char* str, ...);
void mexErrMsgTxt(const string str);
void* mxMalloc(int amount);
void* mxRealloc(void* to_extend, int amount);
void mxFree(void* to_release);
void mexEvalString(const string str);
