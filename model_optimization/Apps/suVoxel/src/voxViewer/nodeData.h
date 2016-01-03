#pragma once
#include <suMesh.h>

class NodeData{
public:
	int AddElement(suMesh::VertexIter  CurrentVertexIter, char Label = 'U')
	{
		VertexVector.push_back(CurrentVertexIter);
		LabelVector.push_back(Label);
		return 0;
	}
	int AddElement(suMesh::FaceHandle faceHandle)
	{
		FaceVector.push_back(faceHandle);
		return 0;
	}
	void ClearAll()
	{
		VertexVector.clear();
		FaceVector.clear();
		LabelVector.clear();
		return;
	}
	int SetLabel(suMesh::VertexIter  CurrentVertexIter, char Label)
	{
		if (Label != 'I' && Label != 'S' && Label != 'O')
			return -1;

		for (unsigned int index = 0; index < VertexVector.size(); index++)
		{
			if (VertexVector[index] == CurrentVertexIter)
			{
				LabelVector[index] = Label;
				return 0;
			}
		}
		return -2;
	}
public:
	//vertex in the node
	std::vector<suMesh::VertexIter>     VertexVector;
	//face handle in the node
	std::vector<suMesh::FaceHandle>     FaceVector;
	std::vector<char>                   LabelVector;
};