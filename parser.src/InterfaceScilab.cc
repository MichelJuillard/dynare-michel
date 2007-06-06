#include <string>

#include "Interface.hh"

namespace interfaces
{
  string comment()
  {
    return "// ";
  }

  string delete_file(string s)
  {
    return "mdelete " + s;
  }

  string file_exist(string s)
  {
    return "file_exist(" + s + ")";
  }

  string compile(string s)
  {
    return "";
  }

  string function_close()
  {
    return "endfunction\n";
  }

  string function_file_extension()
  {
    return ".sci";
  }

  string strvcat(string s1, string s2)
  {
    return "[" + s1 + ";" + s2 + "]";
  }

  string load_model_function_files(string filename)
  {
    return "getf('" + filename + "_static.sci');\ngetf('" + filename + "_dynamic.sci');\n";
  }
}
