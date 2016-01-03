#ifndef SLOT_H_
#define SLOT_H_

#include "suMesh.h"

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
typedef struct SSlotVector
{
	vector<Mesh::VertexIter>	 VertexVector;
	vector<char>	             LabelVector;
	
	// 获取数据区存储数据的个数
	inline unsigned int GetSize() const {return VertexVector.size();}
	
	// 往数据区添加数据项
	// 参数
	//      CurrentVertexIter 指定三角面句柄 
	//      Label             该三角面的标记，默认不知道('U')
	// 返回值
	//   0
	// 参数
	//    三角面动态数组 与 标记动态数组 的数组必须对应
	//    操作是相互绑定：对一个操作，就必须对另外一个有同样的操作
	int AddElement(Mesh::VertexIter  CurrentVertexIter, char Label = 'U')
	{
		VertexVector.push_back(CurrentVertexIter);
		LabelVector.push_back(Label);
		
		return 0;
	}

	// 清空全部的数据
	void ClearAll()
	{
		VertexVector.clear();
		LabelVector.clear();
		return;
	}
	
	// 设置指定三角面的关系
	// 参数
	//     CurrentVertexIter  指定三角面的句柄
	//     Label              指定三角面在另外的网格造型中位置：内部(i)、相交(s)、外部(o)
	// 返回值
	//     -1                 Label 的值有错误
	//     -2                 三角面不再当前的SLOT中
	//      0                 操作成功
	// 说明
	//
	int SetLabel(Mesh::VertexIter  CurrentVertexIter, char Label)
	{
		if (Label != 'I' && Label != 'S' && Label != 'O')
			return -1;
			
		unsigned int index = 0;
		for (; index < VertexVector.size(); index++)
		{
			if (CurrentVertexIter == VertexVector[index])
				break;
		}
		if (index == VertexVector.size())
			return -2;
		
		LabelVector[index] = Label;
		
		return 0;
	}	
} SLOTVECTOR, *PSLOTVECTOR;
//--------------------------------------------------------------------------------------------

typedef struct SSlotFace
{
	vector<Mesh::FaceHandle>	 FacesVector;
	
	// 获取数据区存储数据的个数
	inline unsigned int GetSize() const {return FacesVector.size();}
	
	// 往数据区添加数据项
	// 参数
	//      CurrentFaceIter   指定三角面句柄 
	// 返回值
	//   0
	// 参数
	//
	//
	int AddElement(Mesh::FaceHandle  CurrentFaceHandle)
	{
		FacesVector.push_back(CurrentFaceHandle);
		
		return 0;
	}

	// 清空全部的数据
	void ClearAll()
	{
		FacesVector.clear();
		return;
	}	
} SLOTFACE, *PSLOTFACE;

#endif  //
