#ifndef SHOCKS_H
#define SHOCKS_H
//------------------------------------------------------------------------------
/** \file 
 * \version 1.0
 * \date 04/26/2004
 * \par This file defines the Shocks class.
 */
//------------------------------------------------------------------------------
#include <string>
#include <sstream>
#include <list>
#include <vector>
//------------------------------------------------------------------------------
/*! Shock type enum */
enum shocktype
{
	eDeterministic = 0, //!< Deterministic shock
	eSTDerror = 1,		//!< STD error shock
	eVariance = 2,		//!< Variance shock
	eCovariance = 3,	//!< Covariance shock
	eCorrelation = 4	//!< Corelated shock
}; 
//------------------------------------------------------------------------------
/*!
  \class  Shocks
  \brief  Handles Shocks command  
*/
class Shocks
{
	private :		
		/*! 
		 \class ShockElement
		 \brief Shock element strcuture 
		 */
		struct ShockElement{
		  std::string period1;
		  std::string period2;
		  std::string value;
		};
		/*! 
		 \class ShockElement
		 \brief Shock Structure 
		 */
		struct Shock{
		  int 		id1;
		  int 		id2;
		  shocktype type;
		  std::list<ShockElement> shock_elems;
		  std::string value;
		};		
		/*! Output string of this class */
		std::ostringstream		*output;
		/*! Vector of begin period range */
		std::vector<std::string>		mPeriod1;
		/*! vector of end period range */
		std::vector<std::string>		mPeriod2;
		/*! vector of shock values */
		std::vector<std::string>		mValues;
	public :
		/*! Constructor */
		Shocks();
		/*! Destructor */
		~Shocks();
		/*! Pointer to error function of parser class */                               
		void (* error) (const char* m);
		/*!
			Set output reference
			\param iOutput : reference to an ostringstream
		*/
		void 	setOutput(std::ostringstream* iOutput); 
		/*! Initialize shocks (or mshocks, for multiplicative shocks) block */
		void 	BeginShocks(void);								
		/*! Sets a shock */
		void 	AddShock(shocktype type, int id1, int id2 = 0, std::string value = "");
		/*! Adds a period rage */
		void 	AddPeriod(std::string p1, std::string p2);
		/*! Adds a period */
		void	AddPeriod(std::string p1);
		/*! Adds a value */
		void	AddValue(std::string value);
};
//------------------------------------------------------------------------------
#endif
