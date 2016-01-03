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
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

//Voxel type
typedef enum
{
	EXTERIOR_CELL = 0,           // 外部节点
	BOUNDARY_CELL = 1,           // 表面节点，此表面节点中有三角面
	BOUNDARY_CELL_SPECIAL = 2,   // 表面节点，此表面节点中无三角面，只是三角面的边穿过了该节点
	INTERIOR_CELL = 3,           // 内部节点
	UNDEFINE_CELL = 4,           // 未定义节点
	NON_LEAF_CELL = 5,           // 非叶子节点
	LEAF_CELL = 6                // 重复标定，有更好方法吗？
} NODE_LABEL;



template <class V>
class suOctree
{
public:
	struct Point
	{
		float x;
		float y;
		float z;
		Point():x(0),y(0),z(0){}
		Point(const Point& p2) : x(p2.x), y(p2.y), z(p2.z) {}
		Point(float in_x, float in_y, float in_z) : x(in_x), y(in_y), z(in_z) {}
		Point(const float p2[3]) : x(p2[0]), y(p2[1]), z(p2[2]) {}
		Point& operator=(const Point& p2) { x = p2.x; y = p2.y; z = p2.z; return *this; }
		operator float*() { return &x; }
		operator const float*() const { return &x; }
		Point operator+(const Point& p2) const { return Point(x + p2.x, y + p2.y, z + p2.z); }
		Point operator-(const Point& p2) const { return Point(x - p2.x, y - p2.y, z - p2.z); }
		Point operator*(float f) const { return Point(x*f, y*f, z*f); }
		Point operator/(float f) const { return Point(x/f, y/f, z/f); }
		bool operator< (const Point& p2) const { return x < p2.x && y < p2.y && z < p2.z; }
		bool operator> (const Point& p2) const { return x > p2.x && y > p2.y && z > p2.z; }
		bool operator>=(const Point& p2) const { return x >= p2.x && y >= p2.y && z >= p2.z; }
		bool operator<=(const Point& p2) const { return x <= p2.x && y <= p2.y && z <= p2.z; }
	};

	struct OctreeNode
	{
		V nodeData_;
		OctreeNode* children_[8];
		OctreeNode* parent_;

		unsigned int   xLocCode_;     // X code
		unsigned int   yLocCode_;     // Y code
		unsigned int   zLocCode_;     // Z code
		unsigned int   level_;        // current level, in 32 system unsigned bit is
		                              // represented by 32 bit, so MAX(level) = 32
		NODE_LABEL     label_;

		OctreeNode() :xLocCode_(0), yLocCode_(0), zLocCode_(0), level_(0), label_(UNDEFINE_CELL), parent_(0)
		{
			for (int i = 0; i < 8; i++)
				children_[i] = 0;
		}
		virtual ~OctreeNode()
		{
			for (int i = 0; i < 8; i++)
			if (children_[i])
				delete children_[i];
			parent_ = 0;
		}
	};

	class Callback
	{
	public:
		// Return value: true = continue; false = abort.
		virtual bool operator()(const float min[3], const float max[3], OctreeNode* pCurNode) = 0;
	};
	
public:
	Point min_;
	Point max_;
	Point cellSize_;
	OctreeNode* root_;
	unsigned level_;

public:
	suOctree() :root_(0){}
	suOctree(float minPoint[3], float maxPoint[3], float cellSize[3]) : _min(minPoint), _max(maxPoint), _cellSize(cellSize), _root(0) {}
	virtual ~suOctree(){ clear(); }

	suOctree& operator=(suOctree& Octree);

	void setLevel(int level){ level_ = level; }
	
	void traverse(Callback* callback)
	{
		if (!callback) return;

		traverseRecursive(callback, min_, max_, root_);
	}

	void traverseRecursive(Callback* callback, const Point& currMin, const Point& currMax, OctreeNode* currNode)
	{
		if (!currNode)
			return;
		// leaf node
		if (!currNode->children_[0])
		{
			bool shouldContinue = callback->operator()(currMin, currMax, currNode);
			if (!shouldContinue)
				return;
		}

		Point delta = currMax - currMin;
		Point mid = (delta * 0.5f) + currMin;
		traverseRecursive(callback, currMin, mid, currNode->children_[0]);
		traverseRecursive(callback, Point(mid.x, currMin.y, currMin.z),
			Point(currMax.x, mid.y, mid.z), currNode->children_[1]);
		traverseRecursive(callback, Point(currMin.x, mid.y, currMin.z),
			Point(mid.x, currMax.y, mid.z), currNode->children_[2]);
		traverseRecursive(callback, Point(mid.x, mid.y, currMin.z),
			Point(currMax.x, currMax.y, mid.z), currNode->children_[3]);
		traverseRecursive(callback, Point(currMin.x, currMin.y, mid.z),
			Point(mid.x, mid.y, currMax.z), currNode->children_[4]);
		traverseRecursive(callback, Point(mid.x, currMin.y, mid.z),
			Point(currMax.x, mid.y, currMax.z), currNode->children_[5]);
		traverseRecursive(callback, Point(currMin.x, mid.y, mid.z),
			Point(mid.x, currMax.y, currMax.z), currNode->children_[6]);
		traverseRecursive(callback, mid, currMax, currNode->children_[7]);
	}

	void clear()
	{
		if(root_ != NULL)delete root_;
		root_ = NULL;
	}

	class Iterator
	{
	public:
		Iterator getChild(int i)
		{
			return Iterator(currNode_->children_[i]);
		}
		V* getData()
		{
			if (currNode_)
				return &currNode_->nodeData_;
			else return NULL;
		}
	protected:
		OctreeNode* currNode_;
		Iterator(OctreeNode* node) : currNode_(node) {}
		friend class Octree;
	};

	Iterator getIterator()
	{
		return Iterator(root_);
	}

	//Get & Set bounding box
	void  SetBBox(float minPoint[3], float maxPoint[3])
	{
		clear();
		min_ = Point(minPoint);
		max_ = Point(maxPoint);
		root_ = new OctreeNode();
	}
	void GetBBox(float &xMin, float &yMin, float &zMin, float &xMax, float &yMax, float &zMax)
	{
		xMin = min_.x;
		yMin = min_.y;
		zMin = min_.z;
		xMax = max_.x;
		yMax = max_.y;
		zMax = max_.z;
		return;
	}

	OctreeNode* GetRoot() const { return root_; }

	int Get6Neighbour(OctreeNode* pCurrentNode, vector<OctreeNode*> &nodesArr)
	{
		if (!pCurrentNode) 
			return -1;
		if (pCurrentNode != GetTreeNode(pCurrentNode->level_,
			pCurrentNode->xLocCode_, pCurrentNode->yLocCode_, pCurrentNode->zLocCode_))
			return -2;

		nodesArr.clear();

		OctreeNode*   pNeighbourNode = NULL;
				
		pNeighbourNode = getupneighbor_x(pCurrentNode);
		if (pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);		
		pNeighbourNode = getdownneighbor_x(pCurrentNode);
		if (NULL != pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);
				
		pNeighbourNode = getupneighbor_y(pCurrentNode);
		if (NULL != pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);		
		pNeighbourNode = getdownneighbor_y(pCurrentNode);
		if (NULL != pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);
				
		pNeighbourNode = getupneighbor_z(pCurrentNode);
		if (NULL != pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);		
		pNeighbourNode = getdownneighbor_z(pCurrentNode);
		if (NULL != pNeighbourNode)
			nodesArr.push_back(pNeighbourNode);

		return nodesArr.size();		
	}
	int Get18Neighbour(OctreeNode* pCurrentNode, vector<OctreeNode*> &CellVector);
	int Get26Neighbour(OctreeNode* pCurrentNode, vector<OctreeNode*> &CellVector);


	OctreeNode* GetTreeNode(unsigned int level, unsigned int   xLocCode, unsigned int   yLocCode, unsigned int  zLocCode){
		OctreeNode*  pCurrentNode = root_;
		if (!pCurrentNode)return NULL;

		for (unsigned int i = 0; i <= level; i++)
		{
			if (!pCurrentNode->children_[0])
			{
				// 检测是否是需要的节点
				if (level == pCurrentNode->level_    &&
					xLocCode == pCurrentNode->xLocCode_ &&
					yLocCode == pCurrentNode->yLocCode_ &&
					zLocCode == pCurrentNode->zLocCode_)
					return pCurrentNode;
				else
					return NULL;
			}
			unsigned int locCode = (unsigned int)1 << i;

			unsigned int index = ((xLocCode & locCode) >> i) |
			                    (((yLocCode & locCode) >> i) << 1) |
				                (((zLocCode & locCode) >> i) << 2);

			pCurrentNode = pCurrentNode->children_[index];
		}

		return NULL;
	}
	OctreeNode* GetTreeNode(unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode)
	{
		OctreeNode*  pCurrentNode = root_;
		if (!pCurrentNode) return NULL;

		for (unsigned int i = 0; i < 32; i++)
		{
			if (pCurrentNode->children_[0] == NULL)
				break;

			unsigned int locCode = (unsigned int)1 << i;
			// 跟三个方向上的编码获取孩子节点的游标，是suOctree::Generate()中生成方向编码的逆过程
			unsigned int index = ((xLocCode & locCode) >> i) |
				(((yLocCode & locCode) >> i) << 1) |
				(((zLocCode & locCode) >> i) << 2);

			pCurrentNode = pCurrentNode->children_[index];
		}

		return pCurrentNode;
	}
	/* GetTreeNode(x,y,z) return the node pointer
	 * where the point(x,y,z) belong to.
	 * to test
	 */
	OctreeNode* GetTreeNode(float x, float y, float z)
	{
		OctreeNode  *pCurrentNode = root_;
		Point p = Point(x, y, z);
		Point curNodeMax = max_;
		Point curNodeMin = min_;
		Point delta    = curNodeMax - curNodeMin;
		Point midPoint = curNodeMin + delta / 2;

		unsigned int xLocCode = 1;
		unsigned int yLocCode = 1;
		unsigned int zLocCode = 1;
				
		// If current point is not in the Boundary Box.
		if (p<curNodeMin || p>curNodeMax)
			return NULL;

		while (pCurrentNode)
		{			
			if (!pCurrentNode->children_[0])
				return pCurrentNode;

			// In current node, determine where the point belong to.		
			xLocCode = (p.x >= midPoint.x) ? 1 : 0;
			yLocCode = (p.y >= midPoint.y) ? 1 : 0;
			zLocCode = (p.z >= midPoint.z) ? 1 : 0;

			pCurrentNode = pCurrentNode->children_[ ((zLocCode << 2) | (yLocCode << 1) | xLocCode)];

			//update current node information
			delta = delta / 2.0;
			curNodeMin.x = curNodeMin.x + delta.x * xLocCode;
			curNodeMin.y = curNodeMin.y + delta.y * yLocCode;
			curNodeMin.z = curNodeMin.z + delta.z * zLocCode;
			midPoint = curNodeMin + delta / 2.0;
		}

		return pCurrentNode;
	}
	/*!
	* \function  GetNodeLocation
	*
	* \brief  return node index on x, y and z.	
	* \param  @pNode  a undefined node or a special-boundary node	
	*         @x, y, z  index on x, y and z.
	* \return @true: if a special-boundary nodes is found. false: if no such node is found.
	* \author Yuan Yao
	* \date 2015/10/26
	*/
	bool GetNodeLocation(OctreeNode* pNode, int &x, int &y, int &z)
	{
		Point p;
		if (GetLocation(pNode, p.x, p.y, p.z) < 0) return false;

		Point sizeNode;
		sizeNode = (max_ - min_) / pow(2, pNode->level_);

		p = p - min_;
		x = p.x / sizeNode.x;
		y = p.y / sizeNode.y;
		z = p.z / sizeNode.z;
		return true;
	}

	int GetLocation(OctreeNode* pCurrentNode, float &x, float &y, float &z)
	{
		if (!pCurrentNode)  return -1;
		if (pCurrentNode != GetTreeNode(pCurrentNode->level_, pCurrentNode->xLocCode_, pCurrentNode->yLocCode_, pCurrentNode->zLocCode_))
			return -2;

		float  XCellCurrentMin = min_.x, XCellCurrentMax = max_.x;
		float  YCellCurrentMin = min_.y, YCellCurrentMax = max_.y;
		float  ZCellCurrentMin = min_.z, ZCellCurrentMax = max_.z;

		float xBBox = max_.x - min_.x;
		float yBBox = max_.y - min_.y;
		float zBBox = max_.z - min_.z;

		for (unsigned int i = 0; i < pCurrentNode->level_; i++)
		{
			int nDevide = i + 1;
			float  Xdis = xBBox / (float)(0x1 << nDevide);  //(i+1)*2
			float  Ydis = yBBox / (float)(0x1 << nDevide);
			float  Zdis = zBBox / (float)(0x1 << nDevide);

			if (pCurrentNode->xLocCode_ & (0x1 << i))
				XCellCurrentMin += Xdis;
			else
				XCellCurrentMax -= Xdis;

			if (pCurrentNode->yLocCode_ & (0x1 << i))
				YCellCurrentMin += Ydis;
			else
				YCellCurrentMax -= Ydis;

			if (pCurrentNode->zLocCode_ & (0x1 << i))
				ZCellCurrentMin += Zdis;
			else
				ZCellCurrentMax -= Zdis;
		}

		x = (XCellCurrentMin + XCellCurrentMax) / 2.0f;
		y = (YCellCurrentMin + YCellCurrentMax) / 2.0f;
		z = (ZCellCurrentMin + ZCellCurrentMax) / 2.0f;

		return 0;
	}

private:
	// 对某个方向的编码进行加/减码操作
	unsigned int increaseposition(unsigned int LocCode, unsigned int level)
	{
		unsigned int CurrentLocCode = LocCode;
		High2Low32bit(CurrentLocCode);
		CurrentLocCode = CurrentLocCode + (unsigned int)(0x1 << (32 - level));
		High2Low32bit(CurrentLocCode);

		return CurrentLocCode;
	}
	unsigned int decreaseposition(unsigned int LocCode)
	{
		unsigned int CurrentLocCode = LocCode;
				
		High2Low32bit(CurrentLocCode);		
		CurrentLocCode -= 1;		
		High2Low32bit(CurrentLocCode);

		return CurrentLocCode;
	}

	// Determine if node is on the upper boundary of axis.(for each x,y and z)
	int isinupboundary(unsigned int LocCode, unsigned int level)
	{
		for (unsigned int i = 0; i < level; i++)
		{
			if ((LocCode & 0x1) == 0)	return 0;
			LocCode = LocCode >> 1;
		}
		return 1;
	}
	int  isindownboundary(unsigned int LocCode)
	{
		return (0 == (LocCode ? 1 : 0) );
	}

public:
	OctreeNode* getupneighbor_x(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isinupboundary(pCurrentNode->xLocCode_, pCurrentNode->level_))
		{
			unsigned int xLocCode = increaseposition(pCurrentNode->xLocCode_, pCurrentNode->level_);
			pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode_, pCurrentNode->zLocCode_);
		}
		return pNeighbourNode;
	}
	OctreeNode* getdownneighbor_x(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isindownboundary(pCurrentNode->xLocCode_))
		{
			unsigned int xLocCode = decreaseposition(pCurrentNode->xLocCode_);
			pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode_, pCurrentNode->zLocCode_);
		}
		return pNeighbourNode;
	}
	OctreeNode* getupneighbor_y(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isinupboundary(pCurrentNode->yLocCode_, pCurrentNode->level_))
		{
			unsigned int yLocCode = increaseposition(pCurrentNode->yLocCode_, pCurrentNode->level_);
			pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode_, pCurrentNode->yLocCode_, pCurrentNode->zLocCode_);
		}
		return pNeighbourNode;
	}
	OctreeNode* getdownneighbor_y(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isindownboundary(pCurrentNode->yLocCode_))
		{
			unsigned int yLocCode = decreaseposition(pCurrentNode->yLocCode_);
			pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode_, yLocCode, pCurrentNode->zLocCode_);
		}
		return pNeighbourNode;
	}

	OctreeNode* getupneighbor_z(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isinupboundary(pCurrentNode->zLocCode_, pCurrentNode->level_))
		{
			unsigned int zLocCode = increaseposition(pCurrentNode->zLocCode_, pCurrentNode->level_);
			pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode_, pCurrentNode->yLocCode_, zLocCode);
		}
		return pNeighbourNode;
	}
	OctreeNode* getdownneighbor_z(OctreeNode* pCurrentNode)
	{
		OctreeNode*   pNeighbourNode = NULL;
		if (!isindownboundary(pCurrentNode->zLocCode_))
		{
			unsigned int zLocCode = decreaseposition(pCurrentNode->zLocCode_);
			pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode_, pCurrentNode->yLocCode_, zLocCode);
		}
		return pNeighbourNode;
	}

	/*Generate(unsigned int level)
	 *\brief Generate octave tree with a tree level.
	 *\return: err code
	 *
	**/
	int Generate(unsigned int level)
	{
		if (root_ != NULL)
		{
			delete root_;
		}

		//Check bounding box
		Point delta = max_ - min_;
		if (delta.x < 0.00000001 && delta.y < 0.00000001 && delta.z < 0.00000001)
		{
			return 0;
		}

		OctreeNode*   pCurrentNode = NULL;
		queue<OctreeNode*> NodeQueue;

		root_ = new OctreeNode();

		NodeQueue.push(root_);
		while (!NodeQueue.empty())
		{
			pCurrentNode = NodeQueue.front();
			NodeQueue.pop();

			if (pCurrentNode->level_ == level)	continue;    //access a leaf node

			pCurrentNode->label_ = NON_LEAF_CELL;
			//pCurrentNode->nodeData_ = NULL;

			for (unsigned int k = 0; k < 8; k++)
			{
				pCurrentNode->children_[k] = new OctreeNode();

				//Location code:
			    //             z：1-0 <----> up-bottom  y: 1-0 <----> left-right  x: 1-0 <----> front-back
				//             order: leaf-root
				pCurrentNode->children_[k]->xLocCode_ = ((k & 0x1) << pCurrentNode->level_) | pCurrentNode->xLocCode_;
				pCurrentNode->children_[k]->yLocCode_ = (((k & 0x2) >> 1) << pCurrentNode->level_) | pCurrentNode->yLocCode_;
				pCurrentNode->children_[k]->zLocCode_ = (((k & 0x4) >> 2) << pCurrentNode->level_) | pCurrentNode->zLocCode_;

				pCurrentNode->children_[k]->level_ = pCurrentNode->level_ + 1;
				pCurrentNode->children_[k]->parent_ = pCurrentNode;

				NodeQueue.push(pCurrentNode->children_[k]);
			}
		}
		return 1;
	}

	private:
		/* Exchange high bit and low bit one by one
		* 111001 -> 100111
		**/
		void High2Low32bit(unsigned int &CurrentLocCode)
		{
			for (unsigned int i = 0; i < 16; i++)
			{
				unsigned int bit_high = (CurrentLocCode & ((unsigned)0x1 << (31 - i))) >> (31 - i);
				unsigned int bit_low = (CurrentLocCode & ((unsigned)0x1 << i)) >> i;

				if (bit_high == bit_low)
					continue;

				if (bit_high == 1)
				{
					CurrentLocCode -= (unsigned)0x1 << (31 - i);
					CurrentLocCode += (unsigned)0x1 << i;
				}
				else
				{
					CurrentLocCode += (unsigned)0x1 << (31 - i);
					CurrentLocCode -= (unsigned)0x1 << i;
				}
			}
		}
		

};

