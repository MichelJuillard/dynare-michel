#ifndef _INTERFACE_HH
#define _INTERFACE_HH

namespace interfaces {
    std::string comment();
    std::string delete_file(std::string s);
    std::string file_exist(std::string s);
    std::string compile(std::string s);
    std::string function_close();
    std::string function_file_extension();
    std::string strvcat(std::string s1, std::string s2);
    std::string load_model_function_files(std::string filename);
}

#endif
