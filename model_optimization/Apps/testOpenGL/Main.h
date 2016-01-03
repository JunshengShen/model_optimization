/*
Main.h
Written by Alec Yao

Main.h is included by all source files and includes every header file in the correct order.
*/
#pragma once

//
// Config.h includes a series of #defines used to control compiling options
//
#include "Config.h"

//
// Engine.h includes everything that rarely changes between applications, such as vector/Matrix4 libraries,
// OpenGL/DirectX graphics devices, software rasterizers, etc.
#include <Engines.h>