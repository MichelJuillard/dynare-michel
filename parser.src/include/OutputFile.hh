#ifndef _OUTPUTFILE_HH
#define _OUTPUTFILE_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/26/2004
  \par This file defines the OutputFile class.
*/
//------------------------------------------------------------------------------
#include <fstream>
#include <string>
//------------------------------------------------------------------------------
/*!
  \class  OutputFile
  \brief  Handles opening, writing and closing of output file
*/
class OutputFile
{
private :
  /*! Output file stream */
  ofstream  mOutputFile;
public :
  /*! Flag if set to true, execute "clear all"*/
  bool    clear_all;

  /*! Constructor */
  OutputFile();
  /*! Destructor */
  ~OutputFile();
  /*! Opens a given file and writes some initialization */
  void Open(string iFileName);
  /*! Writes output data from SymbolTable and passed strings to output file */
  void Save(ostringstream& iOutput);
};
//------------------------------------------------------------------------------
#endif
