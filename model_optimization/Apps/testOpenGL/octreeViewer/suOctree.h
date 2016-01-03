#pragma once

#include <iostream>
#include <queue>
#include <vector>
using namespace std;


// ���ص�����
typedef enum
{
	EXTERIOR_CELL = 0,           // �ⲿ�ڵ� 
	BOUNDARY_CELL = 1,           // ����ڵ㣬�˱���ڵ�����������
	BOUNDARY_CELL_SPECIAL = 2,   // ����ڵ㣬�˱���ڵ����������棬ֻ��������ıߴ����˸ýڵ�
	INTERIOR_CELL = 3,           // �ڲ��ڵ�  
	UNDEFINE_CELL = 4,           // δ����ڵ�	
	NON_LEAF_CELL = 5            // ��Ҷ�ӽڵ�	
} CELLLABEL;

// ���صĶ���
class suCell {
public:
	unsigned int   xLocCode;     // X����ı���
	unsigned int   yLocCode;     // Y����ı���
	unsigned int   zLocCode;     // Z����ı���
	unsigned int   level;        // ��ǰ�ڵ�㼶������ÿ�������ϵı���ֻ��32λ�������߲㼶ֻ����32��
	suCell   *pParent;     // ��ǰ�ڵ�ĸ��ڵ�
	suCell   *pChildren;   // ��ǰ�ڵ�ĺ��ӽڵ�

	CELLLABEL      label;        // ��ǰ�ڵ������
	void           *data;        // �������ݣ�Ŀǰ��û�ж�����Ҫ������
	//     ��ǰ�ڵ����ͣ��ⲿ���ڲ�����Ե
	//     �����ڲ������������ǣ��Ƿ�ɷ���
	//     �����ڲ���������������
public:
	// ���캯��
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

	// ����/��ȡ�˲������ڵĿռ�
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

	// ���ռ䰴ָ����С(dx X dy X dz)���ֿռ�
	int Generate(float dx, float dy, float dz);

	// ���ռ䰴ָ����С(dx X dy X dz)���ֿռ�
	int Generate(unsigned int level);

	// ���ٰ˲���
	void DeleteTree();

	// ���ư˲���
	suCell* Clone();

	// ��ȡ/���ð˲����ĸ��ڵ�
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

	// ��ȡ��ǰ�ڵ��6����
	int Get6Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// ��ȡ��ǰ�ڵ��18����
	int Get18Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// ��ȡ��ǰ�ڵ��26����
	int Get26Neighbour(suCell* pCurrentNode, vector<suCell*> &CellVector);

	// ��ȡָ��λ�õĽڵ�
	suCell* GetTreeNode(unsigned int level, unsigned int   xLocCode, unsigned int   yLocCode, unsigned int  zLocCode);
	// ��ȡָ�����뷽ʽ��Ҷ�ӽڵ�
	suCell* GetTreeNode(unsigned int  xLocCode, unsigned int  yLocCode, unsigned int  zLocCode);
	// ��ȡָ�������ڽڵ��
	suCell* GetTreeNode(float x, float y, float z);

	// ��ȡָ���ڵ������λ��
	int GetLocation(suCell* pCurrentNode, float &x, float &y, float &z);

private:

	float xMin, yMin, zMin;
	float xMax, yMax, zMax;
	suCell* pRoot;

	// ���ٰ˲���
	void freeTree(suCell* pCurrentRoot);

	// ��ĳ������ı�����м�/�������
	unsigned int increaseposition(unsigned int LocCode, unsigned int level);
	unsigned int decreaseposition(unsigned int LocCode);

	// �ǲ��������Ͻ�λ����
	int isinupboundary(unsigned int LocCode, unsigned int level);
public:

	// ��ȡ��/�½緽�������
	suCell* getupneighbor_x(suCell* pCurrentNode);
	suCell* getdownneighbor_x(suCell* pCurrentNode);

	suCell* getupneighbor_y(suCell* pCurrentNode);
	suCell* getdownneighbor_y(suCell* pCurrentNode);

	suCell* getupneighbor_z(suCell* pCurrentNode);
	suCell* getdownneighbor_z(suCell* pCurrentNode);
};

