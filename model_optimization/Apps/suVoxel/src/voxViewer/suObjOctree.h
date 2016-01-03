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
=======================================
    suObejctOctree<NodeData>     octreeForDraw_;
    octreeForDraw_.SetBBox(bbMin.data(), bbMax.data());
    /// add mesh to
    for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
   {
      octreeForDraw_.AddCell(mesh_.point(v_it), level);
    }
*/
#include <GL/glew.h>
#include <assert.h>
#include <suBaseObj.h>
#include <suOctree.h>
#include "nodeData.h"


///* -- Cube -- */
#define CUBE_NUM_VERT           8
#define CUBE_NUM_FACES          6
#define CUBE_NUM_EDGE_PER_FACE  4
#define CUBE_VERT_PER_OBJ       (CUBE_NUM_FACES*CUBE_NUM_EDGE_PER_FACE)
#define CUBE_VERT_ELEM_PER_OBJ  (CUBE_VERT_PER_OBJ*3)
#define CUBE_VERT_PER_OBJ_TRI   (CUBE_VERT_PER_OBJ+CUBE_NUM_FACES*2)    /* 2 extra edges per face when drawing quads as triangles */


//============================================================
//class CallbackDraw is a traverse callback function
//used to draw a octave space segmentation with OpenGL
//============================================================
class CallbackDraw : public suOctree<NodeData>::Callback
{
public:
	int count;
	void drawBox(GLfloat *pVert, GLubyte *pIndice)
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glVertexPointer(3, GL_FLOAT, 0, pVert);
		glNormalPointer(GL_FLOAT, 0, cube_n);

		for (int i = 0; i < CUBE_NUM_FACES; i++)
		{
			glDrawElements(GL_LINE_LOOP, CUBE_NUM_EDGE_PER_FACE, GL_UNSIGNED_BYTE, cube_vi + i*CUBE_NUM_EDGE_PER_FACE);
		}
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	//used to generate openGL data for rendering octree.
	virtual bool operator()(const float min[3], const float max[3], suOctree<NodeData>::OctreeNode *pCurNode)
	{
		if (pCurNode->label_ == NON_LEAF_CELL  
			|| pCurNode->label_ == LEAF_CELL                   //for test: only draw bounary cells
			|| pCurNode->label_ == UNDEFINE_CELL
			//|| pCurNode->label_ == INTERIOR_CELL  //test
		    || pCurNode->label_ == EXTERIOR_CELL  //test
			//|| pCurNode->label_ == BOUNDARY_CELL  //test
			//|| pCurNode->label_ == BOUNDARY_CELL_SPECIAL  //test
			)
		{
			return true; 
		}
		OpenMesh::Vec3f pmin(min[0], min[1], min[2]);
		OpenMesh::Vec3f pmax(max[0], max[1], max[2]);

		OpenMesh::Vec3f szNode = pmax - pmin;

		std::vector<OpenMesh::Vec3f> vertArr_;	vertArr_.resize(8);
		vertArr_[0] = pmin;
		vertArr_[1] = pmin + OpenMesh::Vec3f(0.0f, szNode[1], 0.0f);
		vertArr_[2] = pmin + OpenMesh::Vec3f(szNode[0], szNode[1], 0.0f);
		vertArr_[3] = pmin + OpenMesh::Vec3f(szNode[0], 0.0f, 0.0f);

		vertArr_[4] = pmin + OpenMesh::Vec3f(0.0f, 0.0f, szNode[2]);
		vertArr_[5] = pmin + OpenMesh::Vec3f(0.0f, szNode[1], szNode[2]);
		vertArr_[6] = pmax;
		vertArr_[7] = pmin + OpenMesh::Vec3f(szNode[0], 0.0f, szNode[2]);

		std::vector<float> vvert;
		for (int i = 0; i < 8; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				vvert.push_back(vertArr_[i][k]);
			}
		}
		GLfloat *pVert = &vvert[0];

		drawBox(pVert, cube_vi);
		// Subdivide cell if needed.
		return true;
	}
public:
	static GLubyte cube_vi[24/*CUBE_VERT_PER_OBJ*/];
	static GLfloat cube_n[18/*CUBE_NUM_FACES * 3*/];
};



//============================================================
//class CallbackMeshSpaceLabel is a traverse callback function
//used to label octree space by mesh.
//============================================================
class CallbackMeshSpacePartition : public suOctree<NodeData>::Callback
{
public:

	//used to generate openGL data for rendering octree.
	virtual bool operator()(const float min[3], const float max[3], suOctree<NodeData>::OctreeNode *pCurNode)
	{
		//pCurNode->label_ = NON_LEAF_CELL;  //may not need
		//pCurNode = new suOctree<NodeData>::OctreeNode();
		//pCurNode->level_ = pCurNode->level_ + 1;
		//pCurNode->parent_ = pCurNode;

		return true;
	}

};


template<class V>
class suObejctOctree : public suOctree<V>, public suBaseObject
{
public:
	std::vector<OctreeNode*> leafNodesArr;     //leaf nodes array	
public:
	virtual void draw(int draw_type = 0)
	{		
		CallbackDraw ct;
		traverse(&ct);
	}
	//A new octree generation method
	void traversePartitionSpace(OctreeNode* currNode, Callback* callback = NULL)
	{
		if (currNode->level_ == level_)
		{
			currNode->label_ = LEAF_CELL;
			leafNodesArr.push_back(currNode);        
			return;
		}

		if (currNode->children_[0] == NULL)
		{
			currNode->label_ = NON_LEAF_CELL;
			for (int i = 0; i < 8; i++)
			{
				currNode->children_[i] = new OctreeNode();
				currNode->children_[i]->level_ = currNode->level_ + 1;
				currNode->children_[i]->parent_ = currNode;

				//generate code
				currNode->children_[i]->xLocCode_ = ((i & 0x1) << currNode->level_) | currNode->xLocCode_;
				currNode->children_[i]->yLocCode_ = (((i & 0x2) >> 1) << currNode->level_) | currNode->yLocCode_;
				currNode->children_[i]->zLocCode_ = (((i & 0x4) >> 2) << currNode->level_) | currNode->zLocCode_;

				if (callback)
				{					
					float min[3];
					float max[3];

					//call back
					callback->operator()(min, max, currNode->children_[i]);
				}
				traversePartitionSpace(currNode->children_[i], callback);
			}
		}
	}



	int Generate(unsigned int level)
	{
		leafNodesArr.clear();

		clear();		
		setLevel(level);

		CallbackMeshSpacePartition ct;
				
		root_ = new OctreeNode();
		traversePartitionSpace(root_ , &ct);

		return 1;
	}


};
