#ifndef _INTERFACE_HH
#define _INTERFACE_HH

using namespace std;

namespace interfaces
{
  string comment();
  string delete_file(string s);
  string file_exist(string s);
  string compile(string s);
  string function_close();
  string function_file_extension();
  string strvcat(string s1, string s2);
  string load_model_function_files(string filename);
}
#endif
