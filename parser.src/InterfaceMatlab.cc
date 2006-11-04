#include <string>

#include "Interface.hh"

namespace interfaces
{
  std::string comment()
  {
    return "% ";
  }

  std::string delete_file(std::string s)
  {
    return "delete " + s;
  }

  std::string file_exist(std::string s)
  {
    return "exist(" + s + ")";
  }

  std::string compile(std::string s)
  {
    return "mex -O" + s + "\n";
  }

  std::string function_close()
  {
    return "";
  }

  std::string function_file_extension()
  {
    return ".m";
  }

  std::string strvcat(std::string s1, std::string s2)
  {
    return "strvcat(" + s1 + "," + s2 + ")";
  }

  std::string load_model_function_files(std::string filename)
  {
    return "";
  }
}
