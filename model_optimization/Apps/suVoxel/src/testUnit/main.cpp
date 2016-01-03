#include "../config.h"
#include <suOctree.h>

//#define NODEBUG 1


//TestUnit Example
UTFUNC(debug2Console)
{
	debugstream.sink (suDebugSinkConsole::sOnly);

	cdebug << "add";
	
	return;
}
int _tmain(int argc, _TCHAR* argv[])
{
	bool state = suUnitTest::gOnly().run ();
	suUnitTest::gOnly().dumpResults (std::cout);

	return state;
}
