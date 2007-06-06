#include <string>

#include "Interface.hh"

namespace interfaces
{
  string comment()
  {
    return "% ";
  }

  string delete_file(string s)
  {
    return "delete " + s;
  }

  string file_exist(string s)
  {
    return "exist('" + s + "')";
  }

  string compile(string s)
  {
    return "mex -O " + s + "\n";
  }

  string function_close()
  {
    return "";
  }

  string function_file_extension()
  {
    return ".m";
  }

  string strvcat(string s1, string s2)
  {
    return "strvcat(" + s1 + "," + s2 + ")";
  }

  string load_model_function_files(string filename)
  {
    return "";
  }
}
