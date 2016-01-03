// ----------------------------------------------------------------------------
/*
 *   @brief      mesh���ػ�����
 *               1. ������ݺ�
 *               2. �趨�˲�����ȼ���ڵ�ߴ磬���ݽڵ�ߴ����˲������
 *               3. ���ְ˲���������Ҷ�ڵ��д洢��������ݣ���������Ƭ�涥�㣬
 *                  ���������������е�λ�ã�bool����ʱʹ�ã������ǿ��������
 *				   �ṹ�д洢��ɫ�����ϵ�������Ϣ��
 *               4. �궨�ڵ㣬ʹ��Jonathan Shewchuk�ṩ��MeshQuery���ж������ڲ��ṹ��
 *   @output     ���ػ��ļ���vtk��������: ������� / ��������
 *   @author     Yuan Yao - fly2mars@gmail.com
 *   @date       2015/01/18
 *   @todo       
 *               1. Add volume computing
 *               2. Recode octree
 *               3. Testing
 *   @license    This file is part of the 3D Printing Cloud Platform library, RMEC.
 ***/
#include "config.h"
#include "suOptionHelper.h"
#include "suHybridMesh.h"

using namespace std;

void voxelization(suHybridMesh& mesh, suOptionValue& option)
{
	//patition
	if (option.eDevideType == eByTreeDeep)
	{
		mesh.PartitionSpace(option.nTreeDeep);
	}
	//Get leaf nodes
	if (mesh.GenLeafNodeVector() < 0) return ;
	std::cout << "Leaf node size: " << mesh.LeafCellVector.size() << std::endl;

	//label boundary
	mesh.LabelBoundaryNode();

	mesh.LabelInAndOutNode();
	mesh.ExtendLabelNode();

	std::string sOutputFile = option.sOutput;
	mesh.SaveAsVTK(sOutputFile.c_str(), sOutputFile.c_str());
}

// ----------------------------------------------------------------------------
/*
*   @brief      �������ػ�����mesh���
*               1. ���ݷֱ������ػ�
*               2. ͳ���ڲ��ڵ�������
*               3. ͳ�Ʊ߽�ڵ�����
*               4. �������
*   @output     ���(��λmm^2)
*   @author     Yuan Yao - fly2mars@gmail.com
*   @date       2015/02/08
*   @todo
*   @license    This file is part of the 3D Printing Cloud Platform library, RMEC.
***/
double computeVolume(suHybridMesh &mesh)
{	
	//compute a reasonable octree level
	double minLen = 0.01; 
	double len = Math::Max(mesh.width_, mesh.length_, mesh.height_);
	unsigned int num = Math::Ceiling(len / minLen);

	int level = Math::Ceiling(log(num) / log(8));
	mesh.PartitionSpace(level);

	//Get leaf nodes
	if (mesh.GenLeafNodeVector() < 0) return 0;	

	//label boundary
	mesh.LabelBoundaryNode();

	mesh.LabelInAndOutNode();
	mesh.ExtendLabelNode();

	unsigned int nCells = 0;
	for (unsigned int i = 0; i < mesh.LeafCellVector.size(); i++)
	{
		if (mesh.LeafCellVector[i]->label == EXTERIOR_CELL)
		{
			nCells++;
		}
	}
	return mesh.length_ * mesh.height_ * mesh.width_ - nCells * mesh.dx_ * mesh.dy_ * mesh.dz_;
}

int _tmain(int argc, _TCHAR* argv[])
{
	////parse option
	suCMDParser option(argv, argc);
	suOptionValue optionValue(&option);

	if (!optionValue.getParams()) return -1;

	//model file name
	std::string sInputFile = optionValue.sInput;

	suHybridMesh hMesh;
	try
	{
		hMesh.LoadMeshFile(sInputFile.c_str());
		
		//output
		if (optionValue.eFunction == eCOMPUTE_VOLUME)
		{
			std::cout << computeVolume(hMesh);
		}
		else if(optionValue.eFunction == eVOXELIZE){
			voxelization(hMesh, optionValue);
		}

	}
	catch (suException e){
		cout << e.what() << std::endl;
	}
	catch (...){
		cout << "Unknown exception!";
	}
	return 1;
}
