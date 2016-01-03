#include "suOctree.h"
#include "../config.h"

// �����������ݽṹ�Ķ���
//#include "slot.h"
//-------------------------------------------------------------------------------------------------------
/*
 * ���캯��
 */

// ����
//
// ˵��
//     Ĭ�Ϲ��캯��
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

// ����
//     Octree     ��������
// ˵��
//     �������캯��
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

// ����
//
// ˵��
//     ��ֵ���캯��
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

// ��������	
suOctree::~suOctree()
{
	// �ͷŰ˲�����ռ�ռ�
	DeleteTree();
}
//-------------------------------------------------------------------------------------------------------

// ���ð˲������ڵĿռ�
// ����
//     xMin, yMin, zMin       ��͵�����
//     xMax, yMax, zMax       ��ߵ�����
//     signal                 ���ñ�־��1  ֻ���ÿռ�Ĵ�С����Ϊ�˲����ĸ��ڵ����ռ�
//                                      0  ���ÿռ�Ĵ�С��ͬʱΪ�˲����ĸ��ڵ����ռ�  
// ����ֵ
//     -1       �������Ϸ�
//     -2       ������ڵ�ʧ��  
//     0        ���óɹ�
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

	// �ڵ��ʼ���Ѿ���suCell�ĳ�ʼ������ʵ����
	pRoot = new suCell;
	if (NULL == pRoot)
		return -2;

	return 0;
}
//-------------------------------------------------------------------------------------------------------

// ���ռ䰴ָ����С(dx X dy X dz)���ֿռ�
// ����
//     dx X dy X dz   ֹͣ����ʱ��suCell��X,Y,Z��������ĳ߶�
// ����ֵ
//     0         ���ɰ˲����ɹ�
//     1         ���ɰ˲���ʧ��
//    -1         �˲����ĸ��ڵ�Ϊ�գ����ܽ��пռ仮��
//    -2         ָ����СС�ڵ���0
//    -3         �ڴ����ʧ��
// ˵��
//     �ö��еķ�ʽ�����ò㼶���ʵķ�ʽһ��һ��ķ�������
//     Ĭ�Ϸ�����������ǰ���صĴ�С���ܳ���2*2*2
int suOctree::Generate(float dx, float dy, float dz)
{
	// �˲����ĸ��ڵ㲻��Ϊ��
	if (NULL == pRoot)
		return -1;

	// ָ���Ĵ�С����Ϊ�ջ���Ϊ0
	if (dx <= 0.0f || dy <= 0.0f || dz <= 0.0f)
		return -2;

	suCell*        pCurrentNode = NULL;
	suCell*        pChildNode = NULL;
	queue<suCell*> CellQueue;

	pRoot->label = UNDEFINE_CELL;
	pRoot->data = NULL;


	// ����ռ���ӵĴ�С
	float xDis = xMax - xMin;
	float yDis = yMax - yMin;
	float zDis = zMax - zMin;

	CellQueue.push(pRoot);
	while (false == CellQueue.empty())
	{
		// ��ȡ�����е�ͷ�ڵ㣬�����ýڵ�Ӷ������޳�
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		// ��⵱ǰ�����Ƿ���Ҫ����       ֻ��Ϊ�˲���
		const float levelDis = (float)(1 << pCurrentNode->level);
		if (xDis / levelDis <= dx && yDis / levelDis <= dy && zDis / levelDis <= dz)
			continue;

		// ���ѵ�ǰ����
		pCurrentNode->pChildren = new suCell[8];
		if (NULL == pCurrentNode->pChildren)
			return -3;


		pCurrentNode->label = NON_LEAF_CELL;
		pCurrentNode->data = NULL;


		for (unsigned int i = 0; i < 8; i++)
		{
			pChildNode = pCurrentNode->pChildren + i;

			// ������룺  z��1-0 <----> ��-��  y: 1-0 <----> ��-��  x: 1-0 <----> ǰ-��
			//         zyx       ���� 
			//         000       ��-��-ǰ
			//         001       ��-��-��
			//         010       ��-��-ǰ 
			//         011       ��-��-�� 
			//         100       ��-��-ǰ
			//         101       ��-��-��
			//         110       ��-��-ǰ
			//         111		 ��-��-��	
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

// ���ռ仮��Ϊָ���㼶�İ˲���
// ����
//     level     �˲����Ĳ㼶
// ����ֵ
//     0         ���ɰ˲����ɹ�
//     1         ���ɰ˲���ʧ��
//    -1         �˲����ĸ��ڵ�Ϊ�գ����ܽ��пռ仮��
//    -2         ָ����СС�ڵ���0
//    -3         �ڴ����ʧ��
// ˵��
//     �ö��еķ�ʽ�����ò㼶���ʵķ�ʽһ��һ��ķ�������
//     Ĭ�Ϸ�����������ǰ���صĴ�С���ܳ���2*2*2
int suOctree::Generate(unsigned int level)
{
	// �˲����ĸ��ڵ㲻��Ϊ��
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
		// ��ȡ�����е�ͷ�ڵ㣬�����ýڵ�Ӷ������޳�
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		// �����ǰ�ڵ�Ĳ㼶�Ѵ�����߼�����ֹͣ����
		if (pCurrentNode->level >= level)
			continue;

		// ���ѵ�ǰ����
		///����ʱ���ɿ�������������ʣ�ʵ�ֶ�ֱ��ʵ����񻮷�
		pCurrentNode->pChildren = new suCell[8];
		if (NULL == pCurrentNode->pChildren)
			return -3;


		pCurrentNode->label = NON_LEAF_CELL;
		pCurrentNode->data = NULL;


		for (unsigned int i = 0; i < 8; i++)
		{
			pChildNode = pCurrentNode->pChildren + i;

			// ������룺  z��1-0 <----> ��-��  y: 1-0 <----> ��-��  x: 1-0 <----> ǰ-��
			//         zyx       ���� 
			//         000       ��-��-ǰ
			//         001       ��-��-��
			//         010       ��-��-ǰ 
			//         011       ��-��-�� 
			//         100       ��-��-ǰ
			//         101       ��-��-��
			//         110       ��-��-ǰ
			//         111		 ��-��-��	
			//�����ɱ����������еĸߵ�λ�෴
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

// ���ư˲���
// ����
//
// ����ֵ
//     NULL      �ڴ����ʧ��   ����ԭ�˲���������ǿ���
//     ��NULL    ���ư˲����ĸ��ڵ�
// ˵��
//     
suCell* suOctree::Clone()
{
	suCell*        pCloneRoot = NULL;           // ��¡���ĸ��ڵ�
	suCell*        pCloneNode = NULL;           // ��¡���ĵ�ǰ�ڵ�
	suCell*        pCurrentNode = NULL;         // ԭʼ���ĵ�ǰ�ڵ�  
	queue<suCell*> CellQueue;                   // ԭʼ���ķ��ʶ���
	queue<suCell*> CloneQueue;                  // ��¡���ķ��ʶ���

	if (NULL == pRoot)
		return NULL;

	// ���ڵ��ʼ���ڹ��캯����ʵ�ֵ�
	pCloneRoot = new suCell;
	if (NULL == pCloneRoot)
		return NULL;
	// ��¡���ĺ��ӽڵ�������ݵĸ���
	pCloneRoot->xLocCode = pRoot->xLocCode;
	pCloneRoot->yLocCode = pRoot->yLocCode;
	pCloneRoot->zLocCode = pRoot->zLocCode;
	pCloneRoot->level = pRoot->level;

	// ע�����������ֱ�Ӹ�ֵ
	pCloneRoot->pParent = NULL;
	pCloneRoot->pChildren = NULL;

	pCloneRoot->label = UNDEFINE_CELL;
	pCloneRoot->data = NULL;


	CellQueue.push(pRoot);
	CloneQueue.push(pCloneRoot);
	while (false == CellQueue.empty())
	{
		// ��ȡ�����е�ͷ�ڵ㣬�����ýڵ�Ӷ������޳�
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		pCloneNode = CloneQueue.front();
		CloneQueue.pop();

		if (NULL != pCurrentNode->pChildren)
		{
			pCloneNode->pChildren = new suCell[8];
			// �������ʧ�ܵĻ���Ӧ�ð��Ѿ�����Ŀռ��ջ��ٷ���
			if (NULL == pCloneNode->pChildren)
			{
				freeTree(pCloneRoot);
				delete pCloneRoot;
				return NULL;
			}

			for (int i = 0; i < 8; i++)
			{
				// ԭʼ���ĺ��ӽڵ�
				suCell*      pChildNode = pCurrentNode->pChildren + i;

				// ��¡���ĺ��ӽڵ�
				suCell*      pCloneChild = pCloneNode->pChildren + i;

				// ��¡���ĺ��ӽڵ�������ݵĸ���
				pCloneChild->xLocCode = pChildNode->xLocCode;
				pCloneChild->yLocCode = pChildNode->yLocCode;
				pCloneChild->zLocCode = pChildNode->zLocCode;
				pCloneChild->level = pChildNode->level;

				// ע�����������ֱ�Ӹ�ֵ
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

// ��ĳ������ı�����м�/�������
// ˵��
//    �μ�simple and efficient traversal methods for quadtrees and octrees
//    ��/���� ����ֱ�Ӱ��������еķ�������Ҫ�����������õı��뷽ʽ��ͬ��
//    �������Ǹ�λ���ȣ������ǵ�λ����
unsigned int suOctree::increaseposition(unsigned int LocCode, unsigned int level)
{
	unsigned int CurrentLocCode = LocCode;

	// ���뷴ת��ȥ
	Highe2Low32bit(CurrentLocCode);

	// �������
	CurrentLocCode = CurrentLocCode + (unsigned int)(0x1 << (32 - level));

	// ���뷴ת����
	Highe2Low32bit(CurrentLocCode);

	return CurrentLocCode;
}

unsigned int suOctree::decreaseposition(unsigned int LocCode)
{
	unsigned int CurrentLocCode = LocCode;

	// ���뷴ת��ȥ
	Highe2Low32bit(CurrentLocCode);

	// �������
	CurrentLocCode -= 1;

	// ���뷴ת����
	Highe2Low32bit(CurrentLocCode);

	return CurrentLocCode;
}
//-------------------------------------------------------------------------------------------------------

// �ǲ��������Ͻ�λ����
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

// ��ȡ�Ͻ緽�������x����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �Ͻ緽���ָ��  �Ͻ緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
// 
suCell* suOctree::getupneighbor_x(suCell* pCurrentNode)
{
	unsigned int xLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	// �ǲ������Ͻ�
	if (0 == isinupboundary(pCurrentNode->xLocCode, pCurrentNode->level))
	{
		xLocCode = increaseposition(pCurrentNode->xLocCode, pCurrentNode->level);
		pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}

// ��ȡ�½緽�������x����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �½緽���ָ��  �½緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
// 
suCell* suOctree::getdownneighbor_x(suCell* pCurrentNode)
{
	unsigned int xLocCode = 0;
	suCell*        pNeighbourNode = NULL;

	// �ǲ������½�
	if (0 == isindownboundary(pCurrentNode->xLocCode))
	{
		xLocCode = decreaseposition(pCurrentNode->xLocCode);
		pNeighbourNode = GetTreeNode(xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode);
	}

	return pNeighbourNode;
}
//-------------------------------------------------------------------------------------------------------

// ��ȡ��/�½緽�������y����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �Ͻ緽���ָ��  �Ͻ緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
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

// ��ȡ�½緽�������y����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �½緽���ָ��  �½緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
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

// ��ȡ��/�½緽�������z����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �Ͻ緽���ָ��  �Ͻ緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
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

// ��ȡ�½緽�������z����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ�� 
// ����ֵ
//     ��NULL  �½緽���ָ��  �½緽����������
//     NULL    �������󡢲������еĽڵ㡢�Ѿ����Ͻ�λ��
// ˵��
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

// ��ȡ��ǰ�ڵ��6����
// ����
//     pCurrentNode       ��ǰ�ڵ��ָ��
//     CellVector         �洢����ڵ������
// ����ֵ
//     -1      ��ǰ�ڵ�ָ��Ϊ��
//     -2      ��ǰ�ڵ�ָ�벻�����еĽڵ�
//      >= 0   ��ȡ����ڵ������
//  ˵��
//
int suOctree::Get6Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode ���������еĽڵ�
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		throw suException("Illegal node.");

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;

	// x���� ����
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// x���� ����
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	// y���� ����
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// y���� ����
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	// z���� ����
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z���� ����
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// ��ȡ��ǰ�ڵ��18����
int suOctree::Get18Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode ���������еĽڵ�
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		return -2;

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;
	suCell*          pNextNeighbourNode = NULL;

	// x���� ����
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);
		// pNeighbourNode y �� z��������ļ�/����
		// y���� ����
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// y���� ����
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// x���� ����
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// pNeighbourNode y �� z��������ļ�/����
		// y���� ����
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// y���� ����
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);

	}

	// y���� ����
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// y���� ����
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// z���� ����
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z���� ����
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// ��ȡ��ǰ�ڵ��26����
int suOctree::Get26Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector)
{
	// pCurrentNode ���������еĽڵ�
	if (NULL == pCurrentNode)
		return -1;

	if (pCurrentNode != GetTreeNode(pCurrentNode->level,
		pCurrentNode->xLocCode, pCurrentNode->yLocCode, pCurrentNode->zLocCode))
		return -2;

	CellVector.clear();

	suCell*          pNeighbourNode = NULL;
	suCell*          pNextNeighbourNode = NULL;
	suCell*          pSubNextNeighbourNode = NULL;

	// x���� ����
	pNeighbourNode = getupneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);
		// pNeighbourNode y �� z��������ļ�/����
		// y���� ����
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);
			// pNextNeighbourNode z����ļ�/����
			// z���� ����
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z���� ����
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}
		// y���� ����
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z����ļ�/����
			// z���� ����
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z���� ����
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// x���� ����
	pNeighbourNode = getdownneighbor_x(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// pNeighbourNode y �� z��������ļ�/����
		// y���� ����
		pNextNeighbourNode = getupneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z����ļ�/����
			// z���� ����
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z���� ����
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}
		// y���� ����
		pNextNeighbourNode = getdownneighbor_y(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
		{
			CellVector.push_back(pNextNeighbourNode);

			// pNextNeighbourNode z����ļ�/����
			// z���� ����
			pSubNextNeighbourNode = getupneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
			// z���� ����
			pSubNextNeighbourNode = getdownneighbor_z(pNextNeighbourNode);
			if (NULL != pSubNextNeighbourNode)
				CellVector.push_back(pSubNextNeighbourNode);
		}

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// y���� ����
	pNeighbourNode = getupneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}
	// y���� ����
	pNeighbourNode = getdownneighbor_y(pCurrentNode);
	if (NULL != pNeighbourNode)
	{
		CellVector.push_back(pNeighbourNode);

		// z���� ����
		pNextNeighbourNode = getupneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
		// z���� ����
		pNextNeighbourNode = getdownneighbor_z(pNeighbourNode);
		if (NULL != pNextNeighbourNode)
			CellVector.push_back(pNextNeighbourNode);
	}

	// z���� ����
	pNeighbourNode = getupneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);
	// z���� ����
	pNeighbourNode = getdownneighbor_z(pCurrentNode);
	if (NULL != pNeighbourNode)
		CellVector.push_back(pNeighbourNode);

	return CellVector.size();
}
//-------------------------------------------------------------------------------------------------------

// ��ȡָ��λ�õĽڵ�
// ����
//      level                           ָ���ڵ����ڵĲ㼶
//      xLocCode, yLocCode, zLocCode    ָ���ڵ������������ϵı���
// ����ֵ
//      NULL                            ���Ҳ��ɹ�  
//      ��NULL                          ���ҳɹ�������ָ���ڵ��ָ��
// ˵��
//
suCell* suOctree::GetTreeNode(unsigned int level, unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode)
{
	suCell*  pCurrentNode = pRoot;

	// �Ӹ��ڵ㿪ʼ����
	if (NULL == pCurrentNode)
		return NULL;

	for (unsigned int i = 0; i <= level; i++)
	{
		if (NULL == pCurrentNode->pChildren)
		{
			// ����Ƿ�����Ҫ�Ľڵ�
			if (level == pCurrentNode->level    &&
				xLocCode == pCurrentNode->xLocCode &&
				yLocCode == pCurrentNode->yLocCode &&
				zLocCode == pCurrentNode->zLocCode)
				return pCurrentNode;
			else
				return NULL;
		}

		unsigned int locCode = (unsigned int)1 << i;
		// �����������ϵı����ȡ���ӽڵ���α꣬��suOctree::Generate()�����ɷ������������
		unsigned int index = ((xLocCode & locCode) >> i) |
			(((yLocCode & locCode) >> i) << 1) |
			(((zLocCode & locCode) >> i) << 2);

		pCurrentNode = pCurrentNode->pChildren + index;
	}

	return NULL;
}

// ��ȡָ�����뷽ʽ��Ҷ�ӽڵ�
// ����
//     xLocCode, yLocCode, zLocCode    ָ���ڵ������������ϵı���
// ����ֵ
//      NULL                            ���Ҳ��ɹ�  
//      ��NULL                          ���ҳɹ�������ָ���ڵ��ָ��
// ˵��
//
suCell* suOctree::GetTreeNode(unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode)
{
	suCell*  pCurrentNode = pRoot;

	// �Ӹ��ڵ㿪ʼ����
	if (NULL == pCurrentNode)
		return NULL;

	// ��������
	for (unsigned int i = 0; i < 32; i++)
	{
		if (NULL == pCurrentNode->pChildren)
			break;

		unsigned int locCode = (unsigned int)1 << i;
		// �����������ϵı����ȡ���ӽڵ���α꣬��suOctree::Generate()�����ɷ������������
		unsigned int index = ((xLocCode & locCode) >> i) |
			(((yLocCode & locCode) >> i) << 1) |
			(((zLocCode & locCode) >> i) << 2);

		pCurrentNode = pCurrentNode->pChildren + index;
	}

	return pCurrentNode;
}
//-------------------------------------------------------------------------------------------------------

// ��ȡָ�������ڽڵ��
// ����
//     x, y, z     ָ���������
// ����ֵ
//      NULL                            ���Ҳ��ɹ�  
//      ��NULL                          ���ҳɹ�������ָ���ڵ��ָ��
// ˵��
//
suCell* suOctree::GetTreeNode(float x, float y, float z)
{
	suCell*  pCurrentNode = pRoot;
	float  XCellCurrentMin = xMin, XCellCurrentMax = xMax;
	float  YCellCurrentMin = yMin, YCellCurrentMax = yMax;
	float  ZCellCurrentMin = zMin, ZCellCurrentMax = zMax;

	// ���ڸ����Ŀռ��У�ֱ�ӷ���
	if (xMin > x || xMax < x || yMax < y || yMin > y || zMax < z || zMin > z)
		return NULL;

	while (NULL != pCurrentNode)
	{
		// ��ǰ�ڵ����������֮�ڵ�
		if (NULL == pCurrentNode->pChildren)
			return pCurrentNode;

		// ��ǰ�Ľڵ�ĳߴ�		
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
		// ��⵱ǰ�ڵ�����Ľڵ�
		float mid_x = XCellCurrentMin + Xdis / 2.0f;
		float mid_y = YCellCurrentMin + Ydis / 2.0f;
		float mid_z = ZCellCurrentMin + Zdis / 2.0f;

		// ȷ���ڵ�ǰ�ڵ��е�λ�ñ���
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

// ��������ռ�Ŀռ��ͷų���
// ����
//     pCurrentRoot    �����ĸ��ڵ�
// ����ֵ
//
//  ˵��
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

// ���ٰ˲���
// ����
//
// ����ֵ
//
// ˵��
// 
void suOctree::DeleteTree()
{
	if (NULL == pRoot)
		return;


	// �ͷ���������ռ�ռ�
	queue<suCell*> CellQueue;                   // ���ղ㼶���ʰ˲���
	CellQueue.push(pRoot);
	while (false == CellQueue.empty())
	{
		// ��ȡ�����е�ͷ�ڵ㣬�����ýڵ�Ӷ������޳�
		suCell* pCurrentNode = CellQueue.front();
		CellQueue.pop();

		if (NULL == pCurrentNode)
			continue;

		// ��Ҷ�ڵ㣬���ýڵ���ΪNON_LEAF_CELL
		if (NULL != pCurrentNode->pChildren)
		{
			for (unsigned int i = 0; i < 8; i++)
			{
				suCell* pChildNode = pCurrentNode->pChildren + i;
				CellQueue.push(pChildNode);
			}
			continue;
		}

		// �߽�ڵ㣬ͬʱ��������Ч
		if (NULL != pCurrentNode->data)
		{
			// ���ݲ�ͬ�����ȡ���ݣ����ͷſռ�
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

// ��ȡָ���ڵ������λ��
// ����
//		pCurrentNode    ָ���ڵ��ָ��
//      x, y, z         ָ���ڵ�����Ŀռ�����
// ����ֵ
//		-1              ��������
//      -2              ָ���ڵ㲻�ڵ�ǰ����
//      0               �ɹ�
// ˵��
// 
int suOctree::GetLocation(suCell* pCurrentNode, float &x, float &y, float &z)
{
	// �������
	if (NULL == pCurrentNode)
		return -1;

	// ���ָ���ڵ��Ƿ��ڵ�ǰ������
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