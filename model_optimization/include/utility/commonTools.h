#pragma once
/**
Common utility functions.
*/
#include "tstring.h"
#include <strstream>
#include <cstdarg>

namespace SU{
	
/** 
* Format a string object
*/
_tString formatString(const char* format, ...);
/**
* Converts a type to a _tString.
* Convert a type to a string by streaming it. Requires that there's an ostream
* inserter available for type T.
*/
template <class T>
_tString ToString(T num)
{
	_tString strNum = _T("");

   {
      std::strstream buf;

	   buf << num << std::ends;

#ifdef _UNICODE
      std::string temp = buf.str();

      USES_CONVERSION;

      strNum = A2W(temp.c_str());
#else
	   strNum = buf.str();
#endif
	   buf.freeze(false);
   }

   return strNum;
}
}