#include "mex_interface.hh"
#include <cstring>
#include <sstream>

using namespace std;

int
mexPrintf(string str, ...)
{
  va_list vl;
  size_t found, found_p=0;
  found=str.find_first_of("%");
  va_start(vl,str);
  while (found!=string::npos)
  {
    ostringstream tmp_out("");
    //tmp_out.clear();
    tmp_out << "%";
    char c = str[found+1];
    while((c>='0' and c<='9') or c=='.')
      {
        tmp_out << c;
        found++;
        c = str[found+1];
      }
    tmp_out << c;
    switch(c)
      {
        case 'd':
          printf(str.substr(found_p, found-found_p).c_str());
          printf(tmp_out.str().c_str(),va_arg(vl,int));
          break;
        case 'e':
        case 'f':
        case 'g':
          printf(str.substr(found_p, found-found_p).c_str());
          printf(tmp_out.str().c_str(),va_arg(vl,double));
          break;
        case 's':
          printf(str.substr(found_p, found-found_p).c_str());
          printf(tmp_out.str().c_str(),va_arg(vl,char*));
          break;
        case 'x':
          printf(str.substr(found_p, found-found_p).c_str());
          printf(tmp_out.str().c_str(),va_arg(vl,int));
          break;
      }
    found_p = found+2;
    found=str.find_first_of("%",found_p);
   }
  printf(str.substr(found_p, str.size()-found_p+1).c_str());
  return 0;
}

void
mexErrMsgTxt(const string str)
{
  perror(str.c_str());
}

void
mxFree(void* to_release)
{
  free(to_release);
}

void*
mxMalloc(int amount)
{
  return malloc(amount);
}

void*
mxRealloc(void* to_extend, int amount)
{
  return realloc(to_extend, amount);
}


void
mexEvalString(const string str)
{
}
