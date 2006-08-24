#include <fstream>
enum eInterface{
  eMatlab = 0,
  eScilab = 1
};

class interfaces {
 public:
  static eInterface type;
  interfaces(eInterface t){type = t;}
  static inline void set_interface(eInterface t){type = t;}
  static inline std::string comment(void){
    if (type == eMatlab) return "% ";
    else if (type == eScilab) return "// ";
    else return "";
  }
  static inline std::string delete_file(std::string s)
    {
    if (type == eMatlab) 
      return "delete " + s;
    else if (type == eScilab) 
      return "mdelete " + s;
    else return "";
    }      
  static inline std::string file_exist(std::string s)
    {
    if (type == eMatlab) 
      return "exist(" + s + ")";
    else if (type == eScilab) 
      return "file_exist(" + s + ")";
    else return "";
    }
  static inline std::string compile(std::string s)      
    {
    if (type == eMatlab) 
      return "mex -O" + s + "\n";
    else return "";
    }
  static inline std::string function_close(void){
    if (type == eScilab) return "endfunction\n";
    else return "";
  }
  static inline std::string function_file_extension(void){
    if (type == eMatlab) return ".m";
    else if (type == eScilab) return ".sci";
    else return "";
  }
  static inline std::string strvcat(std::string s1,std::string s2)
    {
      if (type == eMatlab) return "strvcat("+s1+","+s2+")";
      else if (type == eScilab) return "[" + s1 + ";" + s2 + "]";
      else return "";
    }
  static inline std::string load_model_function_files(std::string filename)
    {
      if (type == eScilab)
	return "getf('" + filename + "_static.sci');\ngetf('" + filename + "_dynamic.sci');\n";
      else return "";
    }
};

