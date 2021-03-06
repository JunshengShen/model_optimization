#pragma once
/*
ExternalFiles.h

All #includes that are part of C++ core, STL, or other libraries
*/

#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0500
#endif

// Windows
#include <windows.h>
#include <winuser.h>

#ifdef USE_OPENGL
// OpenGL
#include <gl\gl.h>
#include <gl\glu.h>
#endif

// C/C++
#include <math.h>
#include <memory.h>
#include <time.h>
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <exception>

// STL
#include <map>
#include <queue>
#include <set>
#include <hash_map>
#include <hash_set>

// TAUCS (sparse matrix package)
#ifdef USE_TAUCS
extern "C" {
#include "Taucs.h"
};
#endif

#ifdef USE_ZLIB
// ZLib (compression)
#include "zlib.h"

// PNG saving and loading (requires zlib)
#ifdef USE_PNG
#include "png.h"
#endif

#endif

#ifdef USE_SDL
#include "SDL\SDL_image.h"
#endif

//using namespace std avoids typing std::
using namespace std;
using namespace stdext;
