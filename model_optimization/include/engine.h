#pragma once
/////////////////////////////
//      Core               //
/////////////////////////////

//Forward declarations for several classes/structures
#include <Engines/Core/ClassList.h>

//All #includes that are generic or written by other people
#include <Engines\Core\ExternalFiles.h>

//A nice std::vector like vector structure written by Matthew Fisher
#include <Engines\Core\Vector.h>

//Collection of useful constants, macros, and functions
#pragma warning (disable : 4996)  //to disable warning from Microsoft, becuase we use a lot common "unsafe" function like "fopen"...
#include <Engines\Core\Stdhdr.h>

//Generic 32-bit RGBA color structure.  This early include is needed by SpaceVector.
#include <Engines\Math\RGBColor.h>

//SpaceVector, which defines 2D, 3D and 4D vectors and rectangles.  This early include is needed by Grid.
#include <Engines\Math\SpaceVector.h>

//A std::string like string class written by Matthew Fisher, used here to make "num-string convertion" more convinient.
#include <Engines\Core\String.h>
#include <Engines\Core\UnicodeString.h>


/////////////////////////////
//      Utilities          //
/////////////////////////////
#include <utility/unitTest.h>
#include <utility/commonTools.h>
#include <Engines/Tools/Console.h>


/////////////////////////////
//          Math           //
/////////////////////////////
#include <Engines\Math\Matrix4.h>
#include <Engines\Math\DenseMatrix.h>
#include <Engines\Math\PCA.h>

/////////////////////////////
//          Depth          //
/////////////////////////////
#ifdef USE_OPENCV
#include <Engines\Depth\cvUtils.h>
#endif
/****************************/
/*     Graphics Objects     */
/****************************/
#ifdef  USE_ANN
#include <Engines\Graphics Objects\KDTreeN.h>      //General purpose nearest-neighbor structure
#include <Engines\Graphics Objects\KDTree2.h>      //General purpose nearest-neighbor structure
//#include <Engines\Graphics Objects\KDTree3.h>      //General purpose nearest-neighbor structure
#endif