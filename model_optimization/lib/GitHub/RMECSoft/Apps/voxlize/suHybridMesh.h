#pragma once

/*!
* \file suHybridMesh.h
* \date 2015/01/20
*
* \author Wenyu Jin
* \recoded Yuan Yao 
* Contact: fly2mars@gmail.com
*
* \brief This file define a hybrid mesh to generate voxel
*        I rewrite this class and add a new strategy to check if a point is inside of a mesh.
* TODO: 
* \note
**/
#include "config.h"
#include "suMesh.h"
#include "suOctree.h"
#include "../MeshQuery/mesh_query.h"

class suHybridMesh
{
public:
	typedef enum {
		MESHSTATE_NODEFINE = -1,      //No meaning
		MESHSTATE_READ,               //Loaded
		MESHSTATE_WRITE,              //Saved
		MESHSTATE_EXECUTE             //Include two mesh(in Boolean ...)
	} MESHSTATE;
public:
	suHybridMesh(){ Clear(); }
	~suHybridMesh();

	int LoadMeshFile(const char *pFileName) throw (suException);
	int SaveMeshFile(const char *pFileName) throw (suException);

	//partition functions
	int PartitionSpace(unsigned int level) throw (suException);
	int FillCell() throw(suException);
	int GenLeafNodeVector() throw (suException);
	size_t LabelBoundaryNode() throw (suException);
	int labelBoundaryNode(suCell* pNode);
	int LabelInAndOutNode();
	bool insideMeshByQeury(suCell* pNode);
	int ExtendLabelNode();
	void extendLabelNode(suCell* pNode);

	//I/O functions

	/*Save voxels label into a VTK file
	 *\note If @dx, @dy, and @dz set to zero, the original leaf nodes of oct tree are exported.
	 *      Otherwise, exported node size with (@dx, @dy, @dz).
	*/
	int SaveAsVTK(const char *pVTKFileName, const char *pVTKHead = "UnKnown Name", float dx=0, float dy=0, float dz=0);


	void Clear();

	Mesh&     GetMesh(){ return mesh_; }

public:
	Mesh        mesh_;        //Current mesh
	suOctree    octree_;      // 网格造型对应的八叉树 
		
	Mesh::Point bbMin_;       //Bounding box point
	Mesh::Point bbMax_;       //Bounding box point

	Rectangle3f bbox_;        //Bounding box 3d

	vector<suCell*>  LeafCellVector;      // 八叉树的全部叶子节点
	vector<suCell*>  BoundaryNodeVector;  // store all boundary nodes.
	vector<suCell*>  ExteriorNodeVector;  // store all exterior nodes.
	vector<suCell*>  InteriorNodeVector;  // store all nodes inside of mesh.

	MeshObject *pMeshObject_;

	int         meshState_;

	float       length_;
	float       width_;
	float       height_;
	float       dx_;
	float       dy_;
	float       dz_;
	int         level_;

};
