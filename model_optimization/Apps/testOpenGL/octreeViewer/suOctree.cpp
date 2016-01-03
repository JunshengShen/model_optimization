#include "suOctree.h"
#include "../config.h"

// 对数据区数据结构的定义
//#include "slot.h"
//-------------------------------------------------------------------------------------------------------
/*
 * 构造函数
 */

// 参数
//
// 说明
//     默认构造函数
suOctree::suOctree()
{
	xMin = 0.0f;
	yMin = 0.0f;
	zMin = 0.0f;
	xMax = 0.0f;
	yMax = 0.0f;
	zMax = 0.0f;

	pRoot = NULL;
}

// 参数
//     Octree     拷贝对象
// 说明
//     拷贝构造函数
suOctree::suOctree(suOctree& Octree)
{
	if (NULL != (pRoot = Octree.Clone()))
	{
		xMin = Octree.xMin;
		yMin = Octree.yMin;
		zMin = Octree.zMin;
		xMax = Octree.xMax;
		yMax = Octree.yMax;
		zMax = Octree.zMax;
	}
	else
	{
		xMin = 0.0f;
		yMin = 0.0f;
		zMin = 0.0f;
		xMax = 0.0f;
		yMax = 0.0f;
		zMax = 0.0f;
	}
}

// 参数
//
// 说明
//     赋值构造函数
suOctree& suOctree::operator=(suOctree& Octree)
{
	if (NULL != (pRoot = Octree.Clone()))
	{
		xMin = Octree.xMin;
		yMin = Octree.yMin;
		zMin = Octree.zMin;
		xMax = Octree.xMax;
		yMax = Octree.yMax;
		zMax = Octree.zMax;
	}
	else
	{
		xMin = 0.0f;
		yMin = 0.0f;
		zMin = 0.0f;
		xMax = 0.0f;
		yMax = 0.0f;
		zMax = 0.0f;
	}

	return *this;
}
//-------------------------------------------------------------------------------------------------------

// 析构函数	
suOctree::~suOctree()
{
	// 释放八叉树所占空间
	DeleteTree();
}
//-------------------------------------------------------------------------------------------------------

// 设置八叉树所在的空间
// 参数
//     xMin, yMin, zMin       最低点坐标
//     xMax, yMax, zMax       最高点坐标
//     signal                 设置标志：1  只设置空间的大小，不为八叉树的根节点分配空间
//                                      0  设置空间的大小，同时为八叉树的根节点分配空间  
// 返回值
//     -1       参数不合法
//     -2       分配跟节点失败  
//     0        设置成功
int suOctree::SetBox(float xMin, float yMin, float zMin, float xMax, float yMax, float zMax, int signal)
{
	if (xMin >= xMax || yMin >= yMax || zMin >= zMax || (0 != signal && 1 != signal))
		return -1;

	this->xMin = xMin;
	this->yMin = yMin;
	this->zMin = zMin;
	this->xMax = xMax;
	this->yMax = yMax;
	this->zMax = zMax;

	if (1 == signal)
	{
		pRoot = NULL;
		return 0;
	}

	// 节点初始化已经在suCell的初始函数中实现了
	pRoot = new suCell;
	if (NULL == pRoot)
		return -2;

	return 0;
}
//-------------------------------------------------------------------------------------------------------

// 将空间按指定大小(dx X dy X dz)划分空间
// 参数
//     dx X dy X dz   停止分裂时候，suCell在X,Y,Z三个方向的尺度
// 返回值
//     0         生成八叉树成功
//     1         生成八叉树失败
//    -1         八叉树的根节点为空，不能进行空间划分
//    -2         指定大小小于等于0
//    -3         内存分配失败
// 说明
//     用队列的方式，采用层级访问的方式一层一层的分裂体素
//     默认分裂条件：当前体素的大小不能超过2*2*2
int suOctree::Generate(float dx, float dy, float dz)
{
	// 八叉树的根节点不能为空
	if (NULL == pRoot)
		return -1;

	// 指定的大小不能为空或者为0
	if (dx <= 0.0f || dy <= 0.0f || dz <= 0.0f)
		return -2;

	suCell*        pCurrentNode = NULL;
	suCell*        pChildNode = NULL;
	queue<suCell*> CellQueue;

	pRoot->label = UNDEFINE_CELL;
	pRoot->data = NULL;


	// 计算空间盒子的大小
	float xDis = xMax - xMin;
	float yDis = yMax - yMin;
	float zDis = zMax - zMin;

	CellQueue.push(pRoot);
	while (false == CellQueue.empty())
	{
		// 获取队列中的头节点，并将该节点从队列中剔除
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		// 检测当前体素是否需要分裂       只是为了测试
		const float levelDis = (float)(1 << pCurrentNode->level);
		if (xDis / levelDis <= dx && yDis / levelDis <= dy && zDis / levelDis <= dz)
			continue;

		// 分裂当前体素
		pCurrentNode->pChildren = new suCell[8];
		if (NULL == pCurrentNode->pChildren)
			return -3;


		pCurrentNode->label = NON_LEAF_CELL;
		pCurrentNode->data = NULL;


		for (unsigned int i = 0; i < 8; i++)
		{
			pChildNode = pCurrentNode->pChildren + i;

			// 方向编码：  z：1-0 <----> 上-下  y: 1-0 <----> 左-右  x: 1-0 <----> 前-后
			//         zyx       意义 
			//         000       下-左-前
			//         001       下-左-后
			//         010       下-右-前 
			//         011       下-右-后 
			//         100       上-左-前
			//         101       上-左-后
			//         110       上-右-前
			//         111		 上-右-后	
			pChildNode->xLocCode = ((i & 0x1) << pCurrentNode->level) | pCurrentNode->xLocCode;
			pChildNode->yLocCode = (((i & 0x2) >> 1) << pCurrentNode->level) | pCurrentNode->yLocCode;
			pChildNode->zLocCode = (((i & 0x4) >> 2) << pCurrentNode->level) | pCurrentNode->zLocCode;

			pChildNode->level = pCurrentNode->level + 1;
			pChildNode->pParent = pCurrentNode;
			pChildNode->pChildren = NULL;


			pChildNode->label = UNDEFINE_CELL;
			pChildNode->data = NULL;

			CellQueue.push(pChildNode);
		}
	}

	return 0;
}
//-------------------------------------------------------------------------------------------------------

// 将空间划分为指定层级的八叉树
// 参数
//     level     八叉树的层级
// 返回值
//     0         生成八叉树成功
//     1         生成八叉树失败
//    -1         八叉树的根节点为空，不能进行空间划分
//    -2         指定大小小于等于0
//    -3         内存分配失败
// 说明
//     用队列的方式，采用层级访问的方式一层一层的分裂体素
//     默认分裂条件：当前体素的大小不能超过2*2*2
int suOctree::Generate(unsigned int level)
{
	// 八叉树的根节点不能为空
	if (NULL == pRoot)
		return -1;

	suCell*        pCurrentNode = NULL;
	suCell*        pChildNode = NULL;
	queue<suCell*> CellQueue;


	pRoot->label = UNDEFINE_CELL;
	pRoot->data = NULL;


	CellQueue.push(pRoot);
	while (! CellQueue.empty())
	{
		// 获取队列中的头节点，并将该节点从队列中剔除
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		// 如果当前节点的层级已大于最高级别，则停止分裂
		if (pCurrentNode->level >= level)
			continue;

		// 分裂当前体素
		///分列时，可考首先虑体的性质，实现多分辨率的网格划分
		pCurrentNode->pChildren = new suCell[8];
		if (NULL == pCurrentNode->pChildren)
			return -3;


		pCurrentNode->label = NON_LEAF_CELL;
		pCurrentNode->data = NULL;


		for (unsigned int i = 0; i < 8; i++)
		{
			pChildNode = pCurrentNode->pChildren + i;

			// 方向编码：  z：1-0 <----> 上-下  y: 1-0 <----> 左-右  x: 1-0 <----> 前-后
			//         zyx       意义 
			//         000       下-左-前
			//         001       下-左-后
			//         010       下-右-前 
			//         011       下-右-后 
			//         100       上-左-前
			//         101       上-左-后
			//         110       上-右-前
			//         111		 上-右-后	
			//所生成编码与论文中的高低位相反
			pChildNode->xLocCode = ((i & 0x1) << pCurrentNode->level) | pCurrentNode->xLocCode;
			pChildNode->yLocCode = (((i & 0x2) >> 1) << pCurrentNode->level) | pCurrentNode->yLocCode;
			pChildNode->zLocCode = (((i & 0x4) >> 2) << pCurrentNode->level) | pCurrentNode->zLocCode;

			pChildNode->level = pCurrentNode->level + 1;
			pChildNode->pParent = pCurrentNode;
			pChildNode->pChildren = NULL;

			pChildNode->label = UNDEFINE_CELL;
			pChildNode->data = NULL;


			CellQueue.push(pChildNode);
		}	

	}

	return 0;
}
//-------------------------------------------------------------------------------------------------------

// 复制八叉树
// 参数
//
// 返回值
//     NULL      内存分配失败   或者原八叉树本身就是空树
//     非NULL    复制八叉树的根节点
// 说明
//     
suCell* suOctree::Clone()
{
	suCell*        pCloneRoot = NULL;           // 克隆树的根节点
	suCell*        pCloneNode = NULL;           // 克隆树的当前节点
	suCell*        pCurrentNode = NULL;         // 原始树的当前节点  
	queue<suCell*> CellQueue;                   // 原始树的访问队列
	queue<suCell*> CloneQueue;                  // 克隆树的访问队列

	if (NULL == pRoot)
		return NULL;

	// 根节点初始化在构造函数中实现的
	pCloneRoot = new suCell;
	if (NULL == pCloneRoot)
		return NULL;
	// 克隆树的孩子节点相关数据的复制
	pCloneRoot->xLocCode = pRoot->xLocCode;
	pCloneRoot->yLocCode = pRoot->yLocCode;
	pCloneRoot->zLocCode = pRoot->zLocCode;
	pCloneRoot->level = pRoot->level;

	// 注意以下两项不能直接赋值
	pCloneRoot->pParent = NULL;
	pCloneRoot->pChildren = NULL;

	pCloneRoot->label = UNDEFINE_CELL;
	pCloneRoot->data = NULL;


	CellQueue.push(pRoot);
	CloneQueue.push(pCloneRoot);
	while (false == CellQueue.empty())
	{
		// 获取队列中的头节点，并将该节点从队列中剔除
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		pCloneNode = CloneQueue.front();
		CloneQueue.pop();

		if (NULL != pCurrentNode->pChildren)
		{
			pCloneNode->pChildren = new suCell[8];
			// 如果分配失败的话，应该把已经分配的空间收回再返回
			if (NULL == pCloneNode->pChildren)
			{
				freeTree(pCloneRoot);
				delete pCloneRoot;
				return NULL;
			}

			for (int i = 0; i < 8; i++)
			{
				// 原始树的孩子节点
				suCell*      pChildNode = pCurrentNode->pChildren + i;

				// 克隆树的孩子节点
				suCell*      pCloneChild = pCloneNode->pChildren + i;

				// 克隆树的孩子节点相关数据的复制
				pCloneChild->xLocCode = pChildNode->xLocCode;
				pCloneChild->yLocCode = pChildNode->yLocCode;
				pCloneChild->zLocCode = pChildNode->zLocCode;
				pCloneChild->level = pChildNode->level;

				// 注意以下两项不能直接赋值
				pCloneChild->pParent = pCloneNode;
				pCloneChild->pChildren = NULL;


				pCloneChild->label = UNDEFINE_CELL;
				pCloneChild->data = NULL;


				CellQueue.push(pChildNode);
				CloneQueue.push(pCloneChild);
			}
		}
	}

	return pCloneRoot;
}
//-------------------------------------------------------------------------------------------------------

#define  Highe2Low32bit(CurrentLocCode)                                                        \
for (unsigned int i = 0; i < 16; i++)                                                      \
{                                                                                          \
	unsigned int bit_high = (CurrentLocCode & ((unsigned)0x1 << (31 - i))) >> (31 - i);    \
	unsigned int bit_low = (CurrentLocCode & ((unsigned)0x1 << i)) >> i;                  \
	\
if (bit_high == bit_low)                                                               \
	continue;                                                                          \
	\
if (bit_high == 1)                                                                     \
{                                                                                      \
	CurrentLocCode -= (unsigned)0x1 << (31 - i);                                       \
	CurrentLocCode += (unsigned)0x1 << i;                                              \
}                                                                                     \
 else                                                                                 \
{                                                                                      \
	CurrentLocCode += (unsigned)0x1 << (31 - i);                                       \
	CurrentLocCode -= (unsigned)0x1 << i;                                              \
}                                                                                      \
}

// 对某个方向的编码进行加/减码操作
// 说明
//    参见simple and efficient traversal methods for quadtrees and octrees
//    加/减码 不能直接按照论文中的方法，主要是这里所采用的编码方式不同：
//    论文中是高位优先，这里是低位优先
unsigned int suOctree::increaseposition(unsigned int LocCode, unsigned int level)
{
	unsigned int CurrentLocCode = LocCode;

	// 编码反转过去
	Highe2Low32bit(CurrentLocCode);

	// 加码操作
	CurrentLocCode = CurrentLocCode + (unsigned int)(0x1 << (32 - level));

	// 编码反转回来
	Highe2Low32bit(CurrentLocCode);

	return CurrentLocCode;
}

unsigned int suOctree::decreaseposition(unsigned int LocCode)
{
	unsigned int CurrentLocCode = LocCode;

	// 编码反转过去
	Highe2Low32bit(CurrentLocCode);

	// 减码操作
	CurrentLocCode -= 1;

	// 编码反转回来
	Highe2Low32bit(CurrentLocCode);

	return CurrentLocCode;
}
//-------------------------------------------------------------------------------------------------------

// 是不是在最上界位置上
int suOctree::isinupboundary(unsigned int LocCode, unsigned int level)
{
	for (unsigned int i = 0; i < level; i++)
	{
		if (0 == (LocCode & 0x1))
			return 0;

		LocCode = LocCode >> 1;
	}

	return 1;
}
//-------------------------------------------------------------------------------------------------------
#define isindownboundary(LocCode)            \
	(0 == LocCode ? 1 : 0)
//-------------------------------------------------------------------------------------------------------

// 获取上界方向的邻域：x方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  上界方向的指针  上界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getupneighbor_x(suCell* pCurrentNode)
{
	unsigned int xLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	// 是不是最上界
	if (0 == isinupboundary(pCurrentNode->xLocCode, pCurrentNode->level))
	{
		xLocCode = increaseposition(pCurrentNode->xLocCode, pCurrentNode->level);
		pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}

// 获取下界方向的邻域：x方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  下界方向的指针  下界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getdownneighbor_x(suCell* pCurrentNode)
{
	unsigned int xLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	// 是不是最下界
	if (0 == isindownboundary(pCurrentNode->xLocCode))
	{
		xLocCode = decreaseposition(pCurrentNode->xLocCode);
		pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}
//-------------------------------------------------------------------------------------------------------

// 获取上/下界方向的邻域：y方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  上界方向的指针  上界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getupneighbor_y(suCell* pCurrentNode)
{
	unsigned int yLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	if (0 == isinupboundary(pCurrentNode->yLocCode, pCurrentNode->level))
	{
		yLocCode = increaseposition(pCurrentNode->yLocCode, pCurrentNode->level);
		pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode, yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}

// 获取下界方向的邻域：y方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  下界方向的指针  下界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getdownneighbor_y(suCell* pCurrentNode)
{
	unsigned int yLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	if (0 == isindownboundary(pCurrentNode->yLocCode))
	{
		yLocCode = decreaseposition(pCurrentNode->yLocCode);
		pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode, yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}
//-------------------------------------------------------------------------------------------------------

// 获取上/下界方向的邻域：z方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  上界方向的指针  上界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getupneighbor_z(suCell* pCurrentNode)
{
	unsigned int zLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	if (0 == isinupboundary(pCurrentNode->zLocCode, pCurrentNode->level))
	{
		zLocCode = increaseposition(pCurrentNode->zLocCode, pCurrentNode->level);
		pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode, pCurrentNode->yLocCode, zLocCode);
	}
	return pNeighbourNode;
}

// 获取下界方向的邻域：z方向
// 参数
//     pCurrentNode       当前节点的指针 
// 返回值
//     非NULL  下界方向的指针  下界方向的邻域存在
//     NULL    参数错误、不是树中的节点、已经在上界位置
// 说明
// 
suCell* suOctree::getdownneighbor_z(suCell* pCurrentNode)
{
	unsigned int zLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	if (0 == isindownboundary(pCurrentNode->zLocCode))
	{
		zLocCode = decreaseposition(pCurrentNode->zLocCode);
		pNeighbourNode = GetTreeNode(pCurrentNode->xLocCode, pCurrentNode->yLocCode, zLocCode);
	}

	return pNeighbourNode;
}
//-------------------------------------------------------------------------------------------------------

// 获取当前节点的6领域
// 参数
//     pCurrentNode       当前节点的指针
//     CellVector         存储领域节点的向量
// 返回值
//     -1      当前节点指针为空
//     -2      当前节点指针不是树中的节点
//      >= 0   获取领域节点的数量
//  说明
//
int suOctree::Get6Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode 必须是树中的节点
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		throw suException("Illegal node.");

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;

	// x方向 加码
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// x方向 减码
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	// y方向 加码
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// y方向 减码
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	// z方向 加码
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z方向 减码
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// 获取当前节点的18领域
int suOctree::Get18Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode 必须是树中的节点
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		return -2;

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;
	suCell*          pNextNeighbourNode = NULL;

	// x方向 加码
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);
		// pNeighbourNode y 和 z两个方向的加/减码
		// y方向 加码
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// y方向 减码
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// x方向 减码
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// pNeighbourNode y 和 z两个方向的加/减码
		// y方向 加码
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// y方向 减码
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

	}

	// y方向 加码
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// y方向 减码
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// z方向 加码
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z方向 减码
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// 获取当前节点的26领域
int suOctree::Get26Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode 必须是树中的节点
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		return -2;

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;
	suCell*          pNextNeighbourNode = NULL;
	suCell*          pSubNextNeighbourNode = NULL;

	// x方向 加码
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);
		// pNeighbourNode y 和 z两个方向的加/减码
		// y方向 加码
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);
			// pNextNeighbourNode z方向的加/减码
			// z方向 加码
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z方向 减码
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}
		// y方向 减码
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z方向的加/减码
			// z方向 加码
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z方向 减码
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// x方向 减码
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// pNeighbourNode y 和 z两个方向的加/减码
		// y方向 加码
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z方向的加/减码
			// z方向 加码
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z方向 减码
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}
		// y方向 减码
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z方向的加/减码
			// z方向 加码
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z方向 减码
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// y方向 加码
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// y方向 减码
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z方向 加码
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z方向 减码
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// z方向 加码
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z方向 减码
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// 获取指定位置的节点
// 参数
//      level                           指定节点所在的层级
//      xLocCode, yLocCode, zLocCode    指定节点在三个方向上的编码
// 返回值
//      NULL                            查找不成功  
//      非NULL                          查找成功，返回指定节点的指针
// 说明
//
suCell* suOctree::GetTreeNode(unsigned int level, unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode)
{
	suCell*  pCurrentNode = pRoot;

	// 从根节点开始查找
	if (NULL == pCurrentNode)
		return NULL;

	for (unsigned int i = 0; i <= level; i++)
	{
		if (NULL == pCurrentNode->pChildren)
		{
			// 检测是否是需要的节点
			if (level == pCurrentNode->level    &&
				xLocCode == pCurrentNode->xLocCode &&
				yLocCode == pCurrentNode->yLocCode &&
				zLocCode == pCurrentNode->zLocCode)
				return pCurrentNode;
			else
				return NULL;
		}

		unsigned int locCode = (unsigned int)1 << i;
		// 跟三个方向上的编码获取孩子节点的游标，是suOctree::Generate()中生成方向编码的逆过程
		unsigned int index = ((xLocCode & locCode) >> i) |
			(((yLocCode & locCode) >> i) << 1) |
			(((zLocCode & locCode) >> i) << 2);

		pCurrentNode = pCurrentNode->pChildren + index;
	}

	return NULL;
}

// 获取指定编码方式的叶子节点
// 参数
//     xLocCode, yLocCode, zLocCode    指定节点在三个方向上的编码
// 返回值
//      NULL                            查找不成功  
//      非NULL                          查找成功，返回指定节点的指针
// 说明
//
suCell* suOctree::GetTreeNode(unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode)
{
	suCell*  pCurrentNode = pRoot;

	// 从根节点开始查找
	if (NULL == pCurrentNode)
		return NULL;

	// 按最大查找
	for (unsigned int i = 0; i < 32; i++)
	{
		if (NULL == pCurrentNode->pChildren)
			break;

		unsigned int locCode = (unsigned int)1 << i;
		// 跟三个方向上的编码获取孩子节点的游标，是suOctree::Generate()中生成方向编码的逆过程
		unsigned int index = ((xLocCode & locCode) >> i) |
			(((yLocCode & locCode) >> i) << 1) |
			(((zLocCode & locCode) >> i) << 2);

		pCurrentNode = pCurrentNode->pChildren + index;
	}

	return pCurrentNode;
}
//-------------------------------------------------------------------------------------------------------

// 获取指定点所在节点的
// 参数
//     x, y, z     指定点的坐标
// 返回值
//      NULL                            查找不成功  
//      非NULL                          查找成功，返回指定节点的指针
// 说明
//
suCell* suOctree::GetTreeNode(float x, float y, float z)
{
	suCell*  pCurrentNode = pRoot;
	float  XCellCurrentMin = xMin, XCellCurrentMax = xMax;
	float  YCellCurrentMin = yMin, YCellCurrentMax = yMax;
	float  ZCellCurrentMin = zMin, ZCellCurrentMax = zMax;

	// 不在给定的空间中，直接返回
	if (xMin > x || xMax < x || yMax < y || yMin > y || zMax < z || zMin > z)
		return NULL;

	while (NULL != pCurrentNode)
	{
		// 当前节点就是所求质之节点
		if (NULL == pCurrentNode->pChildren)
			return pCurrentNode;

		// 当前的节点的尺寸		
		float  Xdis = (xMax - xMin) / (float)(0x1 << (pCurrentNode->level + 0));
		float  Ydis = (yMax - yMin) / (float)(0x1 << (pCurrentNode->level + 0));
		float  Zdis = (zMax - zMin) / (float)(0x1 << (pCurrentNode->level + 0));

		if (0 != pCurrentNode->level)
		{
			if (pCurrentNode->xLocCode & (0x1 << (pCurrentNode->level - 1)))
				XCellCurrentMin += Xdis / 1.0f;
			else
				XCellCurrentMax -= Xdis / 1.0f;

			if (pCurrentNode->yLocCode & (0x1 << (pCurrentNode->level - 1)))
				YCellCurrentMin += Ydis / 1.0f;
			else
				YCellCurrentMax -= Ydis / 1.0f;

			if (pCurrentNode->zLocCode & (0x1 << (pCurrentNode->level - 1)))
				ZCellCurrentMin += Zdis / 1.0f;
			else
				ZCellCurrentMax -= Zdis / 1.0f;

		}
		// 求解当前节点的中心节点
		float mid_x = XCellCurrentMin + Xdis / 2.0f;
		float mid_y = YCellCurrentMin + Ydis / 2.0f;
		float mid_z = ZCellCurrentMin + Zdis / 2.0f;

		// 确定在当前节点中的位置编码
		unsigned int xLocCode = 1;
		if (mid_x >= x)
			xLocCode = 0;
		unsigned int yLocCode = 1;
		if (mid_y >= y)
			yLocCode = 0;
		unsigned int zLocCode = 1;
		if (mid_z >= z)
			zLocCode = 0;

		pCurrentNode = pCurrentNode->pChildren + ((zLocCode << 2) | (yLocCode << 1) | xLocCode);
	}

	return pCurrentNode;
}
//-------------------------------------------------------------------------------------------------------

// 将子树所占的空间释放出来
// 参数
//     pCurrentRoot    子树的根节点
// 返回值
//
//  说明
//     
void suOctree::freeTree(suCell* pCurrentRoot)
{
	if (NULL != pCurrentRoot)
	{
		if (NULL != pCurrentRoot->pChildren)
		{
			for (unsigned int i = 0; i < 8; i++)
				freeTree(pCurrentRoot->pChildren + i);

			delete[]pCurrentRoot->pChildren;
		}
	}
}
//-------------------------------------------------------------------------------------------------------

// 销毁八叉树
// 参数
//
// 返回值
//
// 说明
// 
void suOctree::DeleteTree()
{
	if (NULL == pRoot)
		return;


	// 释放数据区所占空间
	queue<suCell*> CellQueue;                   // 按照层级访问八叉树
	CellQueue.push(pRoot);
	while (false == CellQueue.empty())
	{
		// 获取队列中的头节点，并将该节点从队列中剔除
		suCell* pCurrentNode = CellQueue.front();
		CellQueue.pop();

		if (NULL == pCurrentNode)
			continue;

		// 非叶节点，将该节点标记为NON_LEAF_CELL
		if (NULL != pCurrentNode->pChildren)
		{
			for (unsigned int i = 0; i < 8; i++)
			{
				suCell* pChildNode = pCurrentNode->pChildren + i;
				CellQueue.push(pChildNode);
			}
			continue;
		}

		// 边界节点，同时数据区有效
		if (NULL != pCurrentNode->data)
		{
			// 根据不同情况获取数据，并释放空间
			if (BOUNDARY_CELL == pCurrentNode->label)
			{
				//PSLOTVECTOR pSlot = (PSLOTVECTOR)pCurrentNode->data;
				//delete pSlot;
			}
			else
			{
				//PSLOTFACE pSlot = (PSLOTFACE)pCurrentNode->data;
				//delete pSlot;
			}
		}
	}


	freeTree(pRoot);
	delete pRoot;
	pRoot = NULL;

	return;
}
//-------------------------------------------------------------------------------------------------------

// 获取指定节点的中心位置
// 参数
//		pCurrentNode    指定节点的指针
//      x, y, z         指定节点所求的空间坐标
// 返回值
//		-1              参数错误
//      -2              指定节点不在当前的树
//      0               成功
// 说明
// 
int suOctree::GetLocation(suCell* pCurrentNode, float &x, float &y, float &z)
{
	// 参数检查
	if (NULL == pCurrentNode)
		return -1;

	// 检查指定节点是否在当前的书中
	if (pCurrentNode != GetTreeNode(pCurrentNode->level, pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		return -2;

	float  XCellCurrentMin = xMin, XCellCurrentMax = xMax;
	float  YCellCurrentMin = yMin, YCellCurrentMax = yMax;
	float  ZCellCurrentMin = zMin, ZCellCurrentMax = zMax;

	for (unsigned int i = 0; i < pCurrentNode->level; i++)
	{
		float  Xdis = (xMax - xMin) / (float)(0x1 << (i + 1));
		float  Ydis = (yMax - yMin) / (float)(0x1 << (i + 1));
		float  Zdis = (zMax - zMin) / (float)(0x1 << (i + 1));

		if (pCurrentNode->xLocCode & (0x1 << i))
			XCellCurrentMin += Xdis;
		else
			XCellCurrentMax -= Xdis;

		if (pCurrentNode->yLocCode & (0x1 << i))
			YCellCurrentMin += Ydis;
		else
			YCellCurrentMax -= Ydis;

		if (pCurrentNode->zLocCode & (0x1 << i))
			ZCellCurrentMin += Zdis;
		else
			ZCellCurrentMax -= Zdis;
	}

	x = (XCellCurrentMin + XCellCurrentMax) / 2.0f;
	y = (YCellCurrentMin + YCellCurrentMax) / 2.0f;
	z = (ZCellCurrentMin + ZCellCurrentMax) / 2.0f;

	return 0;
}
//-------------------------------------------------------------------------------------------------------