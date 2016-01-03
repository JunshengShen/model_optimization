#include <memory>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>

namespace SU{

	_tString formatString(const char* format, ...)
	{
		char buf[1024];
		va_list arglist;
		va_start(arglist, format);
		_vsnprintf(buf, 1024, format, arglist);
		va_end(arglist);
		return _tString(buf);
	}
}