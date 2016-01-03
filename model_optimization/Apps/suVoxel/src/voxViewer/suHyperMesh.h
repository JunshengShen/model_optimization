#pragma once
/*
Copyright (c) 2015, Yuan Yao.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Example:

*/
#include <suBaseObj.h>
#include <suMesh.h>
#include "suObjOctree.h"
#include "nodeData.h"
#include "../MeshQuery/mesh_query.h"


/*
* 这里定义的是八叉树中数据区存储数据的数据结构
* 数据区存储的数据包括：
* 顶点句柄
* 顶点在其他网格造型中的相对位置：
*
*  ---------+------------------+------------------------------------------
*      意义 |  code            | 说明
*  ---------+------------------+------------------------------------------
* 	  未定义 |  'U'             |
* 		内部 |  'I'             |
* 		表面 |  'S'             |
* 		外部 |  'O'             |
*  ---------+------------------+--------------+---------------------------
*/

class suHyperMesh : public suBaseObject
{
public:
	suHyperMesh();
	~suHyperMesh();

	bool LoadMeshFromFile(const char *pFileName);
	void createMeshObject();

	int  PartitionSpace(unsigned int level);
	/* FillNodeData() 
	 * Add mesh vertex into leaf node
	 */
	int  FillNodeData();
	/* LabelBoundaryNode() 
	 * Iteratly label all the node which intersecting with triangles
	 * and including vertex.
	 **/
	int  LabelBoundaryNode();
	int  LabelBoundaryNeighbors();
	bool isBoundaryNode(suObejctOctree<NodeData>::OctreeNode *pNode);
	int  isInteriorNode(suObejctOctree<NodeData>::OctreeNode *pNode);
	bool floodfill();
	bool recursionLabel(suObejctOctree<NodeData>::OctreeNode *pNode, NODE_LABEL label);

	bool saveVTK(const char *pVTKFileName, const char *pVTKHead = "UnKnown Name",
		float dx = 0, float dy = 0, float dz = 0);
	bool test(float dx = 0,float dy = 0,float dz=0);
	bool test1(float dx = 0, float dy = 0, float dz = 0);
	virtual void draw(int draw_type = 0);

	void clear();
	bool isLoad(){ return isLoad_; }

private:
	/* @return 0 if the cube center is on the positive half space (outside) of triangle.
	 * @return 1 of the cube center is on the negative side, but the projection is outside
	 * the triangle.
	 * \note: not finished
	 */	
	int trigCubeSign(suMesh::FaceHandle &fh, suObejctOctree<NodeData>::OctreeNode *pNode);
	// Find all faces around current nodes in the neighbors.
	void addNeighborFaces(std::vector<suObejctOctree<NodeData>::OctreeNode*>  &neigbors,
		std::vector<suMesh::FaceHandle> &FaceHandleQueue);
	bool insideMesh(suObejctOctree<NodeData>::OctreeNode *pNode);	
public:
	suObejctOctree<NodeData> octree;
	suMesh             mesh;	
	suMesh::Point      bbMin, bbMax;

private:
	bool  isLoad_;
	std::vector<suObejctOctree<NodeData>::OctreeNode *>  ExteriorNodeVector;
	std::vector<suObejctOctree<NodeData>::OctreeNode *>  InteriorNodeVector;

	MeshObject *pMeshObject_;
};