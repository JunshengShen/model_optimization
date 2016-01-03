#ifndef SLOT_H_
#define SLOT_H_

#include "suMesh.h"

/*
 * ���ﶨ����ǰ˲������������洢���ݵ����ݽṹ
 * �������洢�����ݰ�����
 * ������
 * �������������������е����λ�ã�
 *
 *  ---------+------------------+------------------------------------------   
 *      ���� |  code            | ˵��
 *  ---------+------------------+------------------------------------------   
 * 	  δ���� |  'U'             |
 * 		�ڲ� |  'I'             |
 * 		���� |  'S'             | 
 * 		�ⲿ |  'O'             |
 *  ---------+------------------+--------------+---------------------------   
 */
typedef struct SSlotVector
{
	vector<Mesh::VertexIter>	 VertexVector;
	vector<char>	             LabelVector;
	
	// ��ȡ�������洢���ݵĸ���
	inline unsigned int GetSize() const {return VertexVector.size();}
	
	// �����������������
	// ����
	//      CurrentVertexIter ָ���������� 
	//      Label             ��������ı�ǣ�Ĭ�ϲ�֪��('U')
	// ����ֵ
	//   0
	// ����
	//    �����涯̬���� �� ��Ƕ�̬���� ����������Ӧ
	//    �������໥�󶨣���һ���������ͱ��������һ����ͬ���Ĳ���
	int AddElement(Mesh::VertexIter  CurrentVertexIter, char Label = 'U')
	{
		VertexVector.push_back(CurrentVertexIter);
		LabelVector.push_back(Label);
		
		return 0;
	}

	// ���ȫ��������
	void ClearAll()
	{
		VertexVector.clear();
		LabelVector.clear();
		return;
	}
	
	// ����ָ��������Ĺ�ϵ
	// ����
	//     CurrentVertexIter  ָ��������ľ��
	//     Label              ָ�������������������������λ�ã��ڲ�(i)���ཻ(s)���ⲿ(o)
	// ����ֵ
	//     -1                 Label ��ֵ�д���
	//     -2                 �����治�ٵ�ǰ��SLOT��
	//      0                 �����ɹ�
	// ˵��
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
	
	// ��ȡ�������洢���ݵĸ���
	inline unsigned int GetSize() const {return FacesVector.size();}
	
	// �����������������
	// ����
	//      CurrentFaceIter   ָ���������� 
	// ����ֵ
	//   0
	// ����
	//
	//
	int AddElement(Mesh::FaceHandle  CurrentFaceHandle)
	{
		FacesVector.push_back(CurrentFaceHandle);
		
		return 0;
	}

	// ���ȫ��������
	void ClearAll()
	{
		FacesVector.clear();
		return;
	}	
} SLOTFACE, *PSLOTFACE;

#endif  //
