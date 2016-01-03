#pragma once

#include <iostream>
#include <queue>
#include <vector>
using namespace std;


// 体素的类型
typedef enum
{
	EXTERIOR_CELL = 0,           // 外部节点 
	BOUNDARY_CELL = 1,           // 表面节点，此表面节点中有三角面
	BOUNDARY_CELL_SPECIAL = 2,   // 表面节点，此表面节点中无三角面，只是三角面的边穿过了该节点
	INTERIOR_CELL = 3,           // 内部节点  
	UNDEFINE_CELL = 4,           // 未定义节点	
	NON_LEAF_CELL = 5            // 非叶子节点	
} CELLLABEL;

// 体素的定义
class suCell {
public:
	unsigned int   xLocCode;     // X方向的编码
	unsigned int   yLocCode;     // Y方向的编码
	unsigned int   zLocCode;     // Z方向的编码
	unsigned int   level;        // 当前节点层级，由于每个方向上的编码只有32位，因此最高层级只能是32级
	suCell   *pParent;     // 当前节点的父节点
	suCell   *pChildren;   // 当前节点的孩子节点

	CELLLABEL      label;        // 当前节点的类型
	void           *data;        // 其他数据，目前还没有定，主要包括：
	//     当前节点类型：外部、内部、边缘
	//     体素内部的三角网格标记：是否可访问
	//     体素内部的三角网格引用
public:
	// 构造函数
	suCell()
	{
		xLocCode = 0;
		yLocCode = 0;
		zLocCode = 0;
		level = 0;
		pParent = NULL;
		pChildren = NULL;

		label = UNDEFINE_CELL;
		data = NULL;

	}
};

class suOctree
{

public:
	suOctree();
	suOctree(suOctree& Octree);
	suOctree& operator=(suOctree& Octree);
	~suOctree();

	// 设置/获取八叉树所在的空间
	int  SetBox(float xMin, float yMin, float zMin, float xMax, float yMax, float zMax, int signal = 0);
	void GetBox(float &xMin, float &yMin, float &zMin, float &xMax, float &yMax, float &zMax)
	{
		xMin = this->xMin;
		yMin = this->yMin;
		zMin = this->zMin;
		xMax = this->xMax;
		yMax = this->yMax;
		zMax = this->zMax;

		return;
	}

	// 将空间按指定大小(dx X dy X dz)划分空间
	int Generate(float dx, float dy, float dz);

	// 将空间按指定大小(dx X dy X dz)划分空间
	int Generate(unsigned int level);

	// 销毁八叉树
	void DeleteTree();

	// 复制八叉树
	suCell* Clone();

	// 获取/设置八叉树的根节点
	suCell* GetRoot() const { return pRoot; }
	int SetRoot(suCell* pCurrentNode)
	{
		if (NULL != pCurrentNode)
		{
			DeleteTree();
			pRoot = pCurrentNode;

			return 0;
		}

		return -1;
	}

	// 获取当前节点的6领域
	int Get6Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// 获取当前节点的18领域
	int Get18Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// 获取当前节点的26领域
	int Get26Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// 获取指定位置的节点
	suCell* GetTreeNode(unsigned int level, unsigned int   xLocCode, unsigned int   yLocCode, unsigned int  zLocCode);
	// 获取指定编码方式的叶子节点
	suCell* GetTreeNode(unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode);
	// 获取指定点所在节点的
	suCell* GetTreeNode(float x, float y, float z);

	// 获取指定节点的中心位置
	int GetLocation(suCell* pCurrentNode, float &x, float &y, float &z);

private:

	float xMin, yMin, zMin;
	float xMax, yMax, zMax;
	suCell* pRoot;

	// 销毁八叉树
	void freeTree(suCell* pCurrentRoot);

	// 对某个方向的编码进行加/减码操作
	unsigned int increaseposition(unsigned int LocCode, unsigned int level);
	unsigned int decreaseposition(unsigned int LocCode);

	// 是不是在最上界位置上
	int isinupboundary(unsigned int LocCode, unsigned int level);
public:

	// 获取上/下界方向的领域
	suCell* getupneighbor_x(suCell* pCurrentNode);
	suCell* getdownneighbor_x(suCell* pCurrentNode);

	suCell* getupneighbor_y(suCell* pCurrentNode);
	suCell* getdownneighbor_y(suCell* pCurrentNode);

	suCell* getupneighbor_z(suCell* pCurrentNode);
	suCell* getdownneighbor_z(suCell* pCurrentNode);
};

