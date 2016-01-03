#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <queue>
#include <vector>

// 关于Mesh方面的头文件
#undef min
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

using namespace OpenMesh::Smoother;

// 相关定义等
struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3f Point;
	typedef OpenMesh::Vec3f Normal;  //point normal
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  suMesh;

