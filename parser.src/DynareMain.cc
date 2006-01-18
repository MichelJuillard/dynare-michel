/*! \file DynareMain.cc
 \version 1.0
 \date 04/26/2004
 \par Main file of Dynare Parser.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

using namespace std;

#include "DynareParser.h"  
#include "DynareScanner.h"
#include "OutputFile.h"
//------------------------------------------------------------------------------
/*!  main function
    \brief Main function of Dynare.
    \param argc Number of command line argumetns from runtime system 
    \param argv Command line arguments from runtime system     
*/
int main(int argc, char** argv) 
{
  OutputFile output_file;  
  ostringstream output;
  int retval = 0;
  try {
    if (argc <2)
    {
    	cout << "Missing model file\n";
    	cout << "Dynare usage : dynare model_file [debug]\n";
    	exit(-1);
    }	

    // create a new scan buffer if a file name was passed to the program 
    ylmm::basic_buffer* input = 0;
    input = new ylmm::basic_buffer(argv[1],false);

    // Create the scanner and parser
    dynare::scanner s(input);
    dynare::parser  p(s);

    // Create the messenger and set it in the parser and scanner 
    ylmm::basic_messenger<ylmm::basic_lock> out;
    s.messenger(out);
    p.messenger(out);
    // Sets the file name in parser instance
    p.set_file_name(argv[1]);
    // Sets string output of scanner
    s.setoutput(&output);
    // Sets string output of parser
    p.setoutput(&output);
    for (int arg=2; arg < argc ; arg++)
    {
    	if (string(argv[arg]) == string("debug"))
    		p.tracing(true); 
    	else if (string(argv[arg]) == string("noclearall"))
    		output_file.clear_all = false;
    }
    // Parse the input
    cout << "Starting Dynare ...\n";
    cout << "Parsing your model file ...\n";
    retval = p.parse(0);
    if (retval != 0)	exit(-1);
    // Execute final instructions
    p.finish();
    string name = argv[1];
    name.erase(name.size()-4,4);
    // Opening and init main Output file (M file)
    output_file.Open(name+".m");
    // Writing remaining string output to output file 
    output_file.Save(output);
	  	
  }
  // Handeling parser and scanner exeptions
  catch (std::exception& e) {
    cout << e.what() << std::endl;
    return 0;
  }
  cout << "Parsing done\n";
  cout << "Starting Matlab computing ...\n";
  return retval;
}
//------------------------------------------------------------------------------
