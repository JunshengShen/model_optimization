#pragma once
/*!
* \file suMesh.h
* \date 2015/01/20 
*
* \author Yuan Yao
* Contact: fly2mars@gmail.com
*
* \brief This file define a mesh based on OpenMesh
*
* TODO: Add addictive properties
*
* \note
**/
#define   _USE_MATH_DEFINES
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

/*
 * \class suMesh
 * \date 2015 / 01 / 20 *
 * \author YuanYao fly2mars@gmail.com
 *
 * \brief This file define a mesh based on OpenMesh
 *
 * TODO: Add addictive properties
 *
 * \note
 **/
class suMesh{
public:
	suMesh();

	void load(const char* filename);
	void save(const char * filename);
	void rescale();
	void append(const Mesh & m);
	Mesh & operator= (const Mesh& m);	

private:

};