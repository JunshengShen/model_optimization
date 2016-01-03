// ----------------------------------------------------------------------------
/*
 *   @brief      mesh体素化工具
 *               1. 计算包容盒
 *               2. 设定八叉树深度计算节点尺寸，根据节点尺寸计算八叉树深度
 *               3. 划分八叉树，各个叶节点中存储自身的数据（包含三角片面顶点，
 *                  顶点在其它网格中的位置（bool运算时使用），我们可以在这个
 *				   结构中存储颜色、材料等其它信息）
 *               4. 标定节点，使用Jonathan Shewchuk提供的MeshQuery库判定网格内部结构。
 *   @output     体素化文件，vtk表达，可用于: 计算体积 / 生成网格
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
*   @brief      利用体素化计算mesh体积
*               1. 根据分辨率体素化
*               2. 统计内部节点数据量
*               3. 统计边界节点数量
*               4. 计算体积
*   @output     体积(单位mm^2)
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
