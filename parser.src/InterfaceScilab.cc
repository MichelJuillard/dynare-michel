#include <string>

#include "Interface.hh"

namespace interfaces
{
  std::string comment()
  {
    return "// ";
  }

  std::string delete_file(std::string s)
  {
    return "mdelete " + s;
  }

  std::string file_exist(std::string s)
  {
    return "file_exist(" + s + ")";
  }

  std::string compile(std::string s)
  {
    return "";
  }

  std::string function_close()
  {
    return "endfunction\n";
  }

  std::string function_file_extension()
  {
    return ".sci";
  }

  std::string strvcat(std::string s1, std::string s2)
  {
    return "[" + s1 + ";" + s2 + "]";
  }

  std::string load_model_function_files(std::string filename)
  {
    return "getf('" + filename + "_static.sci');\ngetf('" + filename + "_dynamic.sci');\n";
  }
}
