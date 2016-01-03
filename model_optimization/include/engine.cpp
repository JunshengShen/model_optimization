#include <engine.h>
/////////////////////////////
//      Core               //
/////////////////////////////
#include <Engines\Core\Stdhdr.cpp>
#include <Engines\Core\UnicodeString.cpp>
#include <Engines\Core\String.cpp>

#include <Engines\Math\RGBColor.cpp>
#include <Engines\Math\SpaceVector.cpp>

/////////////////////////////
//      Utilities          //
/////////////////////////////
#include <utility/timing.cpp>
#include <utility/unitTest.cpp>
#include <utility/commonTools.cpp>
#include <Engines/Tools/Console.cpp>

/////////////////////////////
//          Math           //
/////////////////////////////
#include <Engines\Math\Matrix4.cpp>
#include <Engines\Math\PCA.cpp>

/////////////////////////////
//          Depth          //
/////////////////////////////
#ifdef USE_OPENCV
#include <Engines\Depth\cvUtils.cpp>
#endif

/****************************/
/*     Graphics Objects     */
/****************************/
#ifdef USE_ANN
#include <Engines\Graphics Objects\KDTreeN.cpp>      //General purpose nearest-neighbor structure
#include <Engines\Graphics Objects\KDTree2.cpp>      //General purpose nearest-neighbor structure
//#include <Engines\Graphics Objects\KDTree3.cpp>     //General purpose nearest-neighbor structure
#endif
