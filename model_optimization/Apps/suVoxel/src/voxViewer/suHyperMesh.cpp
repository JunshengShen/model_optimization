#include "suHyperMesh.h"
#include <exception>
#include "overlap.h"
#include "../config.h"
#include "voxel_output.h"
#include "envolution.h"
#include "evo_lib.h"
#include "MeshViewer.h"
#include <Windows.h>

using namespace std;
bool suHyperMesh::LoadMeshFromFile(const char *pFileName)
{
	if (!pFileName)
		return false;

	if (isLoad_)
	{
		clear();
		isLoad_ = false;
	}

	if (OpenMesh::IO::read_mesh(mesh, pFileName))
	{
		suMesh::ConstVertexIter  v_it(mesh.vertices_begin()),
			v_end(mesh.vertices_end());

		bbMin = bbMax = mesh.point(v_it);
		for (; v_it != v_end; ++v_it)
		{
			bbMin.minimize(mesh.point(v_it));
			bbMax.maximize(mesh.point(v_it));
		}

		// compute face & vertex normals
		mesh.update_normals();

		//Create a mesh object (by mesh_query lib)
		createMeshObject();
		isLoad_ = true;
	}

	return isLoad_;
}

void suHyperMesh::clear()
{
	if (pMeshObject_)
	{
		delete pMeshObject_;
		pMeshObject_ = 0;
	}
	octree.clear();
	mesh.clear();
}

int suHyperMesh::PartitionSpace(unsigned int level)
{
	octree.clear();	
	float _bbMin[3];
	float _bbMax[3];
	for (int i = 0; i < 3; i++)
	{
		_bbMin[i] = bbMin[i];
		_bbMax[i] = bbMax[i];

	}
	octree.SetBBox(_bbMin, _bbMax);
	octree.Generate(level);

	FillNodeData();

	return 0;
}
/* FillNodeData add mesh vertex to corresponding node.
* The label of this node is also set to BOUNDARY_CELL.
* Note: In this step, boundary node seemed
*/
int suHyperMesh::FillNodeData()
{
	if (!octree.root_) return -1;

	suMesh::VertexIter   v_it(mesh.vertices_begin());
	suMesh::VertexIter   v_end(mesh.vertices_end());

	for (; v_it != v_end; ++v_it)
	{
		suMesh::Point   point = mesh.point(v_it);

		suOctree<NodeData>::OctreeNode  *pNode = octree.GetTreeNode(point[0], point[1], point[2]);
		if (pNode)
		{
			pNode->nodeData_.AddElement(v_it);
			pNode->label_ = BOUNDARY_CELL;
		}
	}

	return 0;
}


int suHyperMesh::LabelBoundaryNode()
{
	std::vector<suObejctOctree<NodeData>::OctreeNode*> &leafNodesArr = octree.leafNodesArr;
	bool isLabeled = false;
	do
	{
		isLabeled = false;
		for (size_t i = 0; i < leafNodesArr.size(); i++)
		{
			// 1. leaf node with mesh vertex(already labeled in FillNodeData)
			if (leafNodesArr[i]->label_ == BOUNDARY_CELL
				//|| leafNodesArr[i]->label_ == BOUNDARY_CELL_SPECIAL
				) continue;

			// 2. transverse all leaf nodes
			//      -label node which is interact with mesh faces	
			else if (isBoundaryNode(leafNodesArr[i]))
				isLabeled = true;
		}
	} while (isLabeled);

	return 0;
}

bool suHyperMesh::isBoundaryNode(suObejctOctree<NodeData>::OctreeNode *pNode)
{
	// labeled by neighbors
	std::vector<suObejctOctree<NodeData>::OctreeNode*>  neigbors;
	if (octree.Get6Neighbour(pNode, neigbors) == 0)
	{
		return false;
	}

	std::vector<suMesh::FaceHandle> FaceHandleQueue;

	if (pNode->label_ == BOUNDARY_CELL_SPECIAL)
	{
		FaceHandleQueue = pNode->nodeData_.FaceVector;
	}
	/// to avoid repeatly add face to node data.
	int nStartIndex = FaceHandleQueue.size();


	// Find all faces around current nodes in the neighbors.
	addNeighborFaces(neigbors, FaceHandleQueue);
	if (FaceHandleQueue.empty())
	{
		return 0;
	}

	// Use 3D Triangle-Box Overlap to test relationship between current box and triangles
	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf
	/// get box center
	float boxCenter[3];
	octree.GetLocation(pNode, boxCenter[0], boxCenter[1], boxCenter[2]);

	// get "radius" of current node
	float boxHalfSize[3];
	float minX, minY, minZ;
	float maxX, maxY, maxZ;
	octree.GetBBox(minX, minY, minZ, maxX, maxY, maxZ);
	boxHalfSize[0] = (maxX - minX) / (float)((unsigned int)1 << (pNode->level_ + 1));
	boxHalfSize[1] = (maxY - minY) / (float)((unsigned int)1 << (pNode->level_ + 1));
	boxHalfSize[2] = (maxZ - minZ) / (float)((unsigned int)1 << (pNode->level_ + 1));

	float triverts[3][3];

	// triangle is overlapping with current node?
	bool bIsAdd = false;
	for (size_t i = nStartIndex; i < FaceHandleQueue.size() - nStartIndex; i++)
	{
		// for each triangle
		suMesh::FaceHandle curFaceHandle = FaceHandleQueue[i];
		suMesh::FaceVertexIter fv_it = mesh.fv_iter(curFaceHandle);
		suMesh::Point    point; 

		for (int j = 0; j < 3; j++)
		{
			point = mesh.point(fv_it++);
			triverts[j][0] = point[0];
			triverts[j][1] = point[1];
			triverts[j][2] = point[2];
		}
		// if overlapping with current node
		if (TriBoxOverlap(boxCenter, boxHalfSize, triverts))
		{
			pNode->label_ = BOUNDARY_CELL_SPECIAL;

			pNode->nodeData_.AddElement(FaceHandleQueue[i]);
			bIsAdd = true;

			




		}
	}

	return bIsAdd;
}


int suHyperMesh::LabelBoundaryNeighbors()
{
	ExteriorNodeVector.clear();
	InteriorNodeVector.clear();

	std::vector<suObejctOctree<NodeData>::OctreeNode*> &leafNodesArr = octree.leafNodesArr;
	for (size_t i = 0; i < leafNodesArr.size(); i++)
	{
		if (leafNodesArr[i]->label_ == BOUNDARY_CELL
			|| leafNodesArr[i]->label_ == BOUNDARY_CELL_SPECIAL
			)
		{
			suObejctOctree<NodeData>::OctreeNode* pNode = leafNodesArr[i];

			std::vector<suObejctOctree<NodeData>::OctreeNode*>  neighbors;
			if (octree.Get6Neighbour(pNode, neighbors) < 1) continue;
					

			// Labeling neighbor nodes
			for (int j = 0; j < (int)neighbors.size(); j++)
			{
				if (neighbors[j]->label_ == BOUNDARY_CELL
					|| neighbors[j]->label_ == BOUNDARY_CELL_SPECIAL
					|| neighbors[j]->label_ == EXTERIOR_CELL)  continue;
							
				//int iRet = isInteriorNode(neighbors[j]);
				int iRet = insideMesh(neighbors[j]);							
				if ( iRet)
				{
					neighbors[j]->label_ = INTERIOR_CELL;
					InteriorNodeVector.push_back(neighbors[j]);
				}
				else 
				{					
					neighbors[j]->label_ = EXTERIOR_CELL;
					ExteriorNodeVector.push_back(neighbors[j]);
				}
			}
		}
	}
	std::cout << "Label " << ExteriorNodeVector.size() << " exterior nodes and " << 
		    InteriorNodeVector.size() << " interior node.\n";
	return ExteriorNodeVector.size() + InteriorNodeVector.size();
}


//todo: 
//using cube & triangle relation to determine inside or outside
//Function trigCubeSign is not finished.
int suHyperMesh::isInteriorNode(suObejctOctree<NodeData>::OctreeNode *pNode)
{
	int iRet = -3;
	// If it is not an undefined leaf node.
	if (pNode->label_ != LEAF_CELL) return iRet;

	
	// labeled by neighbors
	std::vector<suObejctOctree<NodeData>::OctreeNode*>  neigbors;
	if (octree.Get6Neighbour(pNode, neigbors) == 0)
	{
		return iRet;
	}

	std::vector<suMesh::FaceHandle> FaceHandleQueue;

	addNeighborFaces(neigbors, FaceHandleQueue);

	for (size_t i = 0; i < FaceHandleQueue.size(); i++)
	{
		iRet = trigCubeSign(FaceHandleQueue[i], pNode);
		if (iRet == 1)
			return iRet;
	}
	return iRet;

}

void suHyperMesh::draw(int draw_type /*= 0*/)
{
	octree.draw();
}

suHyperMesh::suHyperMesh()
{
	isLoad_ = 0;
	pMeshObject_ = 0;

	mesh.request_face_normals();
	mesh.request_vertex_normals();
}

void suHyperMesh::addNeighborFaces(std::vector<suObejctOctree<NodeData>::OctreeNode*> &neigbors, std::vector<suMesh::FaceHandle> &FaceHandleQueue)
{
	for (unsigned int i = 0; i < neigbors.size(); i++)
	{
		suObejctOctree<NodeData>::OctreeNode* pLeafNode = neigbors[i];

		/// find neighbors who are boundary nodes		
		if (pLeafNode->label_ == BOUNDARY_CELL)
		{
			NodeData &nodeData = pLeafNode->nodeData_;
			for (int j = 0; j < nodeData.VertexVector.size(); j++)
			{
				/// find neighbor faces by vertex;
				suMesh::VertexIter vIter = nodeData.VertexVector[j];;

				for (suMesh::VertexFaceIter vfIter = mesh.vf_iter(vIter.handle()); vfIter; vfIter++)
				{
					// add new faces into current node					
					bool bExist = false;
					for (size_t ii = 0; ii < FaceHandleQueue.size(); ii++)
					{
						if (vfIter.handle() == FaceHandleQueue[ii])
						{
							bExist = true;
						}
					}
					if (!bExist)
					{
						FaceHandleQueue.push_back(vfIter.handle());
					}
				}
			}
		}
		else if (pLeafNode->label_ == BOUNDARY_CELL_SPECIAL)
		{
			//add new face from special boundary node
			NodeData &nodeData = pLeafNode->nodeData_;
			for (size_t j = 0; j < nodeData.FaceVector.size(); j++)
			{
				// add new faces into current node					
				bool bExist = false;
				for (size_t ii = 0; ii < FaceHandleQueue.size(); ii++)
				{
					if (nodeData.FaceVector[j] == FaceHandleQueue[ii])
					{
						bExist = true;
					}
				}
				if (!bExist)
				{
					FaceHandleQueue.push_back(nodeData.FaceVector[j]);
				}
			}
		}
	}
}

int suHyperMesh::trigCubeSign(suMesh::FaceHandle &fh, suObejctOctree<NodeData>::OctreeNode *pNode)
{
	float boxCenter[3];
	octree.GetLocation(pNode, boxCenter[0], boxCenter[1], boxCenter[2]);
    
	Vec3f nodeCenter(boxCenter[0], boxCenter[1], boxCenter[2]);

	Vec3f triverts[3];
	suMesh::FaceVertexIter fv_it = mesh.fv_iter(fh);

	for (int i = 0; i < 3; i++ )
	{
		suMesh::Point    point = mesh.point(fv_it);
		triverts[i].x = point[0];
		triverts[i].y = point[1];
		triverts[i].z = point[2];
		fv_it++;
	}

	suMesh::Point nf = mesh.normal(*fv_it);
	Vec3f n(nf[0],nf[1],nf[2]);
	Vec3f d = nodeCenter - triverts[0];
	if (Vec3f::Dot(n,d) > 0)
	{
		return 0;
	}

	/*float normN = n.Length();
	d = nodeCenter - (Vec3f::Dot(n, d) / normN) * n;
	Vec3f n0 = Vec3f::Cross(triverts[1] - d, triverts[2] - d);*/
	
	return 1;
}

/**
*\brief  A popular way to check whether or not a ray has
*        crossed a mesh is to count the number of ray/triangle
*        intersections.
*\param
*         @pNode: An undefined node
*         @x @y @z   indicate the ray orientation, eg. (0,0,1).
*\return  node label
*\
*/
bool suHyperMesh::insideMesh(suObejctOctree<NodeData>::OctreeNode *pNode)
{
	float x, y, z;
	if (octree.GetLocation(pNode, x, y, z) != 0) return false;
	Vec3f pCenter = Vec3f(x, y, z);

	double center[3];
	center[0] = x;
	center[1] = y;
	center[2] = z;

	return point_inside_mesh(center, pMeshObject_);
	
}

bool suHyperMesh::floodfill()
{
	//floodfill interior node
	std::vector<suObejctOctree<NodeData>::OctreeNode *>  outsides = ExteriorNodeVector;
	std::vector<suObejctOctree<NodeData>::OctreeNode *>  insides = InteriorNodeVector;
	for (size_t i = 0; i < InteriorNodeVector.size(); i++)
	{	
		recursionLabel(InteriorNodeVector[i], INTERIOR_CELL);
	}

	for (size_t i = 0; i < ExteriorNodeVector.size(); i++)
	{		
		recursionLabel(ExteriorNodeVector[i], EXTERIOR_CELL);
	}
	return true;
}

bool suHyperMesh::recursionLabel(suObejctOctree<NodeData>::OctreeNode *pNode, NODE_LABEL label)
{
	std::vector<suObejctOctree<NodeData>::OctreeNode*>  neighbors;
	if (octree.Get6Neighbour(pNode, neighbors) <= 0)
	{
		if (pNode->label_ == LEAF_CELL)
		{
			pNode->label_ = label;
			if (label == EXTERIOR_CELL) ExteriorNodeVector.push_back(pNode);
			else if (label == INTERIOR_CELL) InteriorNodeVector.push_back(pNode);
		}
			
		return false;
	}

	for (int j = 0; j < neighbors.size(); j++)
	{
		//std::cout << j;
		if (neighbors[j]->label_ == LEAF_CELL)
		{
			//std::cout << "leaf: " << label << std::endl;
			neighbors[j]->label_ = label;

			if (label == EXTERIOR_CELL) ExteriorNodeVector.push_back(pNode);
			else if (label == INTERIOR_CELL) InteriorNodeVector.push_back(pNode);

			recursionLabel(neighbors[j], label);
		}
		//std::cout << neighbors.size() << std::endl;
	}
	return true;

}

void suHyperMesh::createMeshObject()
{
	if (pMeshObject_ != NULL)
	{
		delete pMeshObject_;
	}
	double *vertex = new double[3 * mesh.n_vertices()];
	int *trig = new int[3 * mesh.n_faces()];

	int idx = 0;
	for (suMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
		for (int jj = 0; jj < 3; jj++){
			vertex[idx++] = mesh.point(*v_it)[jj];
		}
	}
	idx = 0;
	for (suMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		suMesh::FaceVertexIter v_it = mesh.fv_begin(f_it.handle());
		for (; v_it != mesh.fv_end(f_it.handle()); ++v_it)
		{

			trig[idx++] = v_it.handle().idx();
		}
	}
	pMeshObject_ = construct_mesh_object(mesh.n_vertices(), vertex, mesh.n_faces(), trig);
}



suHyperMesh::~suHyperMesh()
{
	clear();
	
}
int voxel_min;

bool suHyperMesh::test1(float dx, float dy, float dz)
{
	int oct_level;
	cin >> oct_level;
	
	string file;
	cin >> file;
	//if (!file)
	//	return false;

	if (isLoad_)
	{
		clear();
		isLoad_ = false;
	}

	if (OpenMesh::IO::read_mesh(mesh, file))
	{
		suMesh::ConstVertexIter  v_it(mesh.vertices_begin()),
			v_end(mesh.vertices_end());

		bbMin = bbMax = mesh.point(v_it);
		for (; v_it != v_end; ++v_it)
		{
			bbMin.minimize(mesh.point(v_it));
			bbMax.maximize(mesh.point(v_it));
		}

		// compute face & vertex normals
		mesh.update_normals();

		//Create a mesh object (by mesh_query lib)
		createMeshObject();
		isLoad_ = true;
	}



	return isLoad_;	
	//输入中要加入材料属性与要求   暂时不做 12.29
	//初始化加入到循环中（体素化  元胞初始化）然后开始演化(基于规则的推理过程）
	//演化结束调用外部工具 等待结束 读入结果
	//外部评估 评估之后与输入参数比较 更新状态
	//判断 多个规则并行判断  把规则做成规则行书（不是循环条件）  现在有三个规则
	//反馈做成一个状态

	//while(do_)
	 ///if(not_init)
	  ////init()
	 ///evo()
	 ///oofem()
	 ///do_=if_do()
		
	//read_point_information("d://KITTEN.out.m0.1.vtu");

	//return 0;
}

bool suHyperMesh::test(float dx,float dy, float dz)
{
	
	float minx, miny, minz;
	float maxx, maxy, maxz;
	octree.GetBBox(minx, miny, minz, maxx, maxy, maxz);

	suObejctOctree<NodeData>::Point vMin(minx, miny, minz);
	suObejctOctree<NodeData>::Point vMax(maxx, maxy, maxz);
	suObejctOctree<NodeData>::Point vSize = vMax - vMin;

	if (dx == 0)
	{
		float nSize = pow(2, octree.level_);
		dx = vSize.x / nSize;
		dy = vSize.y / nSize;
		dz = vSize.z / nSize;
	}

	unsigned int Xdim = (unsigned int)(vSize.x / dx);
	unsigned int Ydim = (unsigned int)(vSize.y / dy);
	unsigned int Zdim = (unsigned int)(vSize.z / dz);

	/*vtk << "DIMENSIONS " << Xdim << " " << Ydim << " " << Zdim << endl;

	// 数据区的其他信息
	vtk << "ASPECT_RATIO 1 1 1" << endl;
	vtk << "ORIGIN 0 0 0" << endl;
	vtk << "POINT_DATA " << Xdim * Ydim * Zdim << endl;
	vtk << "SCALARS volume_scalars char 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl << endl;*/
	int count = 0;//在后面输出txt文件时输出体素的编号
	int voxel_number = 0;
	cout << octree.level_ << endl;
	string test_dir;
	test_dir = "D:\\oofem\\build2.3\\Debug\\test.txt";
	string test1_dir;
	test1_dir = "D:\\oofem\\build2.3\\Debug\\test1.txt";
	string test2_dir;
	test2_dir = "D:\\oofem\\build2.3\\Debug\\test2.txt";
	fstream outfile;
	outfile.open(test_dir, ios::out);
	outfile.close();








	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suObejctOctree<NodeData>::OctreeNode *pNode = octree.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);//体素中心的坐标，能够据此计算出编号
				if (NULL != pNode)
				{
					// 内部节点输出 1 
					if (INTERIOR_CELL == pNode->label_)
					{
						//vtk << 1 << " ";
						//fstream outfile;
						//outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY<< "  " << IndexZ;
						//outfile << "LSpace " << count++ << "	 nodes  8 ";
						//voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						//voxel_output asd(pChildNode, level);
						//outfile.close();
						voxel_number++;

						//asd->output_point();
						//delete asd;

					}
					// 边界节点输出 2
					else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
					{
						//vtk << 2 << " ";
						//fstream outfile;
						//outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
						//outfile << "LSpace " << count++ << "	 nodes  8 ";
						//voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						//outfile.close();
						voxel_number++;
						//voxel_output asd(pChildNode, level);
						//asd->output_point();
						//delete asd;
					}
					// 外部节点输出 0
					//else if (EXTERIOR_CELL == pNode->label_)
					//vtk << 0 << " ";

					// 其他情况输出 3    这种情况应该是不出现的
					//else
					//vtk << 3 << " ";
					//}
					//else
					//vtk << 4 << " ";
					//}
					//vtk << endl;
					//}
					//vtk << endl;
				}
			}
		}
	}
	//fstream outfile;

	outfile.open(test_dir, ios::app);
	outfile << "Majnun.out" << endl;
	outfile << "test of Brick elements with nlgeo 1(strain is the Green-Lagrangian strain) rotated as a rigid body" << endl;
	outfile << "#NonLinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#LinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#nsteps 5 rtolv 1.e-6 stiffMode 1 controlmode 1 maxiter 100" << endl;
	outfile << "#vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "#domain 3d" << endl;
	outfile << "#OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "LinearStatic nsteps 3 nmodules 1" << endl;
	outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "domain 3d" << endl;
	outfile << "OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "ndofman " << pow((pow(2, octree.level_) + 1), 3) << " nelem " << voxel_number << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 " << endl;
	outfile.close();
	cout << voxel_number;
	coor(vSize.x, vSize.y, vSize.z, octree.level_);





	vector<ca_envolution> auto_cell(pow(pow(2, octree.level_), 3));
	vector<int> ca(voxel_number);



	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suObejctOctree<NodeData>::OctreeNode *pNode = octree.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);//体素中心的坐标，能够据此计算出编号
				if (NULL != pNode)
				{
					// 内部节点输出 1 
					if (INTERIOR_CELL == pNode->label_)
					{
						//vtk << 1 << " ";
						fstream outfile;
						outfile.open(test_dir, ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY<< "  " << IndexZ;
						outfile << "LSpace " << ++count << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						int temp = merton(IndexX, IndexY, IndexZ, octree.level_);
						//auto_cell[temp].x = IndexX;
						//auto_cell[temp].y = IndexY;
						//auto_cell[temp].z = IndexZ;
						auto_cell[temp].level = octree.level_;
						auto_cell[temp].merton = temp;
						auto_cell[temp].lable = 0;
						auto_cell[temp].out = 1;
						ca[count - 1] = temp;
						//voxel_output asd(pChildNode, level);
						outfile.close();


						asd->output_point();
						delete asd;

					}
					// 边界节点输出 2
					else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
					{
						//vtk << 2 << " ";
						fstream outfile;
						outfile.open(test_dir, ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
						outfile << "LSpace " << ++count << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						int temp = merton(IndexX, IndexY, IndexZ, octree.level_);
						//auto_cell[temp].x = IndexX;
						//auto_cell[temp].y = IndexY;
						//auto_cell[temp].z = IndexZ;
						auto_cell[temp].level = octree.level_;
						auto_cell[temp].merton = temp;
						auto_cell[temp].lable = 1;
						auto_cell[temp].location = 1;
						auto_cell[temp].out = 1;
						ca[count - 1] = temp;
						outfile.close();

						//voxel_output asd(pChildNode, level);
						asd->output_point();
						delete asd;
					}
					// 外部节点输出 0
					//else if (EXTERIOR_CELL == pNode->label_)
						//vtk << 0 << " ";

					// 其他情况输出 3    这种情况应该是不出现的
					//else
						//vtk << 3 << " ";
				}
				//else
					//vtk << 4 << " ";
			}
			//vtk << endl;
		}
		//vtk << endl;
	}
	//fstream outfile;
	outfile.open(test_dir, ios::app);
	outfile << "SimpleCS 1" << endl << "IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0" << endl << "BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.0" << endl << "BoundaryCondition  2 loadTimeFunction 1 prescribedvalue 0.5" << endl << "PiecewiseLinFunction 1 npoints 2 t 2 0. 1000. f(t) 2 0. 1000." << endl;
	outfile.close();

	/*for (int i = 0; i < count; i++)
	{
	cout << ca[i] << " "<<i<<" ";
	}*/




	/*cout<<endl;
	cout << auto_cell[3882].location << endl;
	cout << auto_cell[3874].location << endl;
	cout << auto_cell[3889].location << endl;
	cout << auto_cell[3873].location << endl;
	cout << auto_cell[3879].location << endl;
	cout << auto_cell[3847].location << endl;
	*/
	int add = 2;
	for (;;)
	{
		int location_temp = 0;
		for (int i = 0; i < count; i++)//遍历所有体素
		{
			if (auto_cell[ca[i]].location == 0)//判断本身是否为0
			{

				int *n_merton = six_n_merton(ca[i], octree.level_);//输入莫顿序  输出6-领域莫顿序    六位数组
				/*int compare = 0;
				for (int j = 0; j < 6; j++)//计算出六领域中loaction最大的
				{
				if (auto_cell[n_merton[j]].location > compare);
				compare = auto_cell[n_merton[j]].location;
				}
				if (compare > 0)//location是否赋值  是则赋值   计数加1
				{
				auto_cell[ca[i]].location_temp = compare + 1;
				location_temp++;
				//cout << ca[i] << " " << auto_cell[ca[i]].location_temp << " " << auto_cell[ca[i]].location << "       ";
				}
				else
				cout << auto_cell[n_merton[0]].location << auto_cell[n_merton[1]].location << auto_cell[n_merton[2]].location << auto_cell[n_merton[3]].location << auto_cell[n_merton[4]].location << auto_cell[n_merton[5]].location;
				*/
				if ((auto_cell[n_merton[0]].location + auto_cell[n_merton[1]].location + auto_cell[n_merton[2]].location + auto_cell[n_merton[3]].location + auto_cell[n_merton[4]].location + auto_cell[n_merton[5]].location) != 0)
				{
					auto_cell[ca[i]].location_temp = add;

				}
			}
		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].location_temp == add)
			{
				auto_cell[ca[i]].location = add;
				//cout << auto_cell[ca[i]].location << " " ;
				auto_cell[ca[i]].location_temp = 0;
				location_temp++;
			}
		}
		//cout << "location_temp="<<location_temp << endl;
		if (location_temp == 0)
			break;
		add++;
	}

	int void_voxel = 0;
	
	cout << endl;
	cout << void_voxel << endl;
	//输入中要加入材料属性与要求
	//初始化加入到循环中（体素化  元胞初始化）然后开始演化(基于规则的推理过程）
	//演化结束调用外部工具 等待结束 读入结果
	//外部评估 评估之后与输入参数比较 更新状态
	//判断 多个规则并行判断  把规则做成规则行书（不是循环条件）  现在有三个规则
	//反馈做成一个状态
	for (int k = 1; k <= 6; k++)
	{
		
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].lable == 0 && auto_cell[ca[i]].out == 1)
			{
				int *n_merton = six_n_merton(ca[i], octree.level_);
				/*if ((auto_cell[ca[i]].location >= auto_cell[ca[0]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[1]].location) &&
				(auto_cell[ca[i]].location >= auto_cell[ca[2]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[3]].location) &&
				(auto_cell[ca[i]].location >= auto_cell[ca[4]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[5]].location))
				{

				auto_cell[ca[i]].out_temp = 0;
				}*/

				int max = 0;
				for (int j = 0; j < 26; j++)
				{
					if ((auto_cell[n_merton[j]].location>max) && (auto_cell[n_merton[j]].out == 1))
						max = auto_cell[n_merton[j]].location;
				}
				delete n_merton;
				if (auto_cell[ca[i]].location >= max)
				{
					auto_cell[ca[i]].out_temp = 0;//

				}//cout << max << "  " << auto_cell[ca[i]].location << "  " << auto_cell[ca[i]].out_temp << "||";
			}

		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out_temp == 0)
			{
				auto_cell[ca[i]].out = 0;
				auto_cell[ca[i]].out_temp = 1;

			}

		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 0)
			{
				void_voxel++;

			}

		}
		voxel_min = count;
		if (voxel_min > void_voxel)
		{
			voxel_min = void_voxel;
		}
		cout << k << ":" << void_voxel << ";" << count << endl;
		void_voxel = 0;


		cout << "min=" << count - voxel_min;




		fstream outstl;
		outstl.open("d://STL.stl", ios::out);
		outstl << "solid \"sample\"" << endl;
		outstl.close();
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 1)
			{
				int *n_merton = six_n_merton(ca[i], octree.level_);
				/*cout << "(" << ca[i] << ")" << endl;

				for (int j = 0; j < 6; j++)
				{
				cout << "merton:" << n_merton[j] << " ";
				}
				*/
				int *code_ = code(ca[i], octree.level_);
				/*for (int j = 0; j < 3; j++)
				{
				cout << "code:" << code_[j] << " ";
				}*/

				if (code_[0] == (pow(2, octree.level_) - 1))
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1)*dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] * dz) << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) *dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[0]].out == 0)//x+
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1)*dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] * dz) << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) *dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();



				}
				if (code_[0] == 0)
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "-1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] * dz) << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "-1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[1]].out == 0)//x-
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "-1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] * dz) << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "-1 0 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				if (code_[1] == (pow(2, octree.level_) - 1))
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[2]].out == 0)//y+
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				if (code_[1] == 0)
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 -1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 -1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[3]].out == 0)//y-
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 -1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1)*dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 -1 0" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				if (code_[2] == (pow(2, octree.level_) - 1))
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[4]].out == 0)//z+
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << code_[2] * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << code_[2] * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				if (code_[2] == 0)
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
				else if (auto_cell[n_merton[5]].out == 0)//z-
				{
					outstl.open("d://STL.stl", ios::app);
					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;


					outstl << "  facet normal " << "0 0 1" << endl;
					outstl << "    outer loop" << endl;
					outstl << "      vertex " << (code_[0] - 1) * dx << ' ' << (code_[1] - 1)*dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << (code_[1] - 1) * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "      vertex " << code_[0] * dx << ' ' << code_[1] * dy << ' ' << (code_[2] - 1) * dz << endl;
					outstl << "    endloop" << endl;
					outstl << "  endfacet" << endl;

					outstl.close();
				}
			}
		}
		outstl.open("d://STL.stl", ios::app);
		outstl << "endsolid \"sample\"";






		int voxel_number1 = 0;


		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 1)
			{
				voxel_number1++;
				//fstream outfile;
				//outfile.open("d://test1.txt", ios::app);
				//int *asddd = code(ca[i], octree.level_);
				//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
				//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
				//outfile << "LSpace " << ++fuck << "	 nodes  8 ";
				//voxel_output* asd = new voxel_output(asddd[0], asddd[1], asddd[2], octree.level_);
				//outfile.close();
				//voxel_output asd(pChildNode, level);
				//asd->output_point1();
				//delete asd;
			}
		}

		outfile.open(test1_dir, ios::out);
		outfile << "Majnun.out" << endl;
		outfile << "test of Brick elements with nlgeo 1(strain is the Green-Lagrangian strain) rotated as a rigid body" << endl;
		outfile << "#NonLinearStatic  nmsteps 1 nsteps 1 " << endl;
		outfile << "#LinearStatic  nmsteps 1 nsteps 1 " << endl;
		outfile << "#nsteps 5 rtolv 1.e-6 stiffMode 1 controlmode 1 maxiter 100" << endl;
		outfile << "#vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
		outfile << "#domain 3d" << endl;
		outfile << "#OutputManager tstep_all dofman_all element_all" << endl;
		outfile << "LinearStatic nsteps 3 nmodules 1" << endl;
		outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
		outfile << "domain 3d" << endl;
		outfile << "OutputManager tstep_all dofman_all element_all" << endl;
		outfile << "ndofman " << pow((pow(2, octree.level_) + 1), 3) << " nelem " << voxel_number1 << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 " << endl;
		outfile.close();
		coor1(vSize.x, vSize.y, vSize.z, octree.level_);




		int temp_count = 0;
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 1)
			{
				fstream outfile;
				outfile.open(test1_dir, ios::app);
				int *asddd = code(ca[i], octree.level_);
				//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
				//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
				outfile << "LSpace " << ++temp_count << "	 nodes  8 ";
				voxel_output* asd = new voxel_output(asddd[0], asddd[1], asddd[2], octree.level_);
				outfile.close();
				//voxel_output asd(pChildNode, level);
				asd->output_point1();
				delete asd;
			}
		}



		outfile.open(test1_dir, ios::app);
		outfile << "SimpleCS 1" << endl << "IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0" << endl << "BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.0" << endl << "BoundaryCondition  2 loadTimeFunction 1 prescribedvalue 0.5" << endl << "PiecewiseLinFunction 1 npoints 2 t 2 0. 1000. f(t) 2 0. 1000." << endl;
		outfile.close();

















		outfile.open(test2_dir, ios::out);
		outfile << "Majnun.out" << endl;
		outfile << "Simple bending of a cantilever beam, quadratic elements." << endl;
		outfile << "LinearStatic nsteps 2 controllmode 1 rtolv 1.e-3 nmodules 1" << endl;
		outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 6 1 2 4 5 27 28 stype 1" << endl;
		outfile << "domain 3d" << endl;
		outfile << "OutputManager tstep_all dofman_all element_all" << endl;
		outfile << "ndofman " << pow((pow(2, octree.level_) + 1), 3) << " nelem " << voxel_number1 * 5 << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1" << endl;
		outfile.close();
		coor2(vSize.x, vSize.y, vSize.z, octree.level_);




		int number = 0;//一个六面体体素生成四个四个四面体体素
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 1)
			{
				fstream outfile;
				//outfile.open("d://test1.txt", ios::app);
				int *asddd = code(ca[i], octree.level_);
				//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
				//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
				//outfile << "LSpace " << ++fuck << "	 nodes  8 ";
				voxel_output* asd = new voxel_output(asddd[0], asddd[1], asddd[2], octree.level_);
				//outfile.close();
				//voxel_output asd(pChildNode, level);
				asd->output_point2(number);
				delete asd;
				number++;
			}
		}



		outfile.open(test2_dir, ios::app);
		outfile << "SimpleCS 1 thick 1.0 width 1.0" << endl;
		outfile << "IsoLE 1 d 0.0102 E 200.0 n 0.394  tAlpha 0.0000780" << endl;
		outfile << "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0" << endl;
		outfile << "ConstantSurfaceLoad 2 ndofs 3 loadType 2 Components 3 0.0 -84.118 0.0 loadTimeFunction 1" << endl;
		outfile << "ConstantFunction 1 f(t) 1.0" << endl;
		outfile.close();
		system("cd /d D:\\oofem\\build2.3\\Debug &oofem.exe&oofem -f test1.txt");
		read_point_information("D:\\oofem\\build2.3\\Debug\\Majnun.out.m0.1.vtu");
	}
	//system("cd /d D:\\oofem\\build2.3\\Debug &oofem.exe&oofem -f test1.txt");
	return 0;
}




bool
suHyperMesh::saveVTK(const char *pVTKFileName, const char *pVTKHead, 
	float dx, float dy, float dz)
{
	if (NULL == pVTKFileName || NULL == pVTKHead)
		return false;

	std::ofstream vtk(pVTKFileName);

	if (!vtk)
		return false;

	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << pVTKHead << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET STRUCTURED_POINTS" << endl;

	// 
	float minx, miny, minz;
	float maxx, maxy, maxz;
	octree.GetBBox(minx, miny, minz, maxx, maxy, maxz);

	suObejctOctree<NodeData>::Point vMin(minx, miny, minz);
	suObejctOctree<NodeData>::Point vMax(maxx, maxy, maxz);
	suObejctOctree<NodeData>::Point vSize = vMax - vMin;
		
	if (dx == 0)
	{
		float nSize = pow(2, octree.level_);
		dx = vSize.x / nSize;
		dy = vSize.y / nSize;
		dz = vSize.z / nSize;		
	}

	unsigned int Xdim = (unsigned int)(vSize.x / dx);
	unsigned int Ydim = (unsigned int)(vSize.y / dy);
	unsigned int Zdim = (unsigned int)(vSize.z / dz);

	vtk << "DIMENSIONS " << Xdim << " " << Ydim << " " << Zdim << endl;

	// 数据区的其他信息
	vtk << "ASPECT_RATIO 1 1 1" << endl;
	vtk << "ORIGIN 0 0 0" << endl;
	vtk << "POINT_DATA " << Xdim * Ydim * Zdim << endl;
	vtk << "SCALARS volume_scalars char 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl << endl;
	int count = 0;//在后面输出txt文件时输出体素的编号
	int voxel_number = 0;
	cout << octree.level_<<endl;
	
	fstream outfile;
	outfile.open("d://test.txt", ios::out);
	outfile.close();








	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suObejctOctree<NodeData>::OctreeNode *pNode = octree.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);//体素中心的坐标，能够据此计算出编号
				if (NULL != pNode)
				{
					// 内部节点输出 1 
					if (INTERIOR_CELL == pNode->label_)
					{
						//vtk << 1 << " ";
						//fstream outfile;
						//outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY<< "  " << IndexZ;
						//outfile << "LSpace " << count++ << "	 nodes  8 ";
						//voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						//voxel_output asd(pChildNode, level);
						//outfile.close();
						voxel_number++;

						//asd->output_point();
						//delete asd;

					}
					// 边界节点输出 2
					else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
					{
						//vtk << 2 << " ";
						//fstream outfile;
						//outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
						//outfile << "LSpace " << count++ << "	 nodes  8 ";
						//voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						//outfile.close();
						voxel_number++;
						//voxel_output asd(pChildNode, level);
						//asd->output_point();
						//delete asd;
					}
					// 外部节点输出 0
					//else if (EXTERIOR_CELL == pNode->label_)
					//vtk << 0 << " ";

					// 其他情况输出 3    这种情况应该是不出现的
					//else
					//vtk << 3 << " ";
					//}
					//else
					//vtk << 4 << " ";
					//}
					//vtk << endl;
					//}
					//vtk << endl;
				}
			}
		}
	}
	//fstream outfile;

	outfile.open("d://test.txt", ios::app);
	outfile << "Majnun.out" << endl;
	outfile << "test of Brick elements with nlgeo 1(strain is the Green-Lagrangian strain) rotated as a rigid body" << endl;
	outfile << "#NonLinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#LinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#nsteps 5 rtolv 1.e-6 stiffMode 1 controlmode 1 maxiter 100" << endl;
	outfile << "#vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "#domain 3d" << endl;
	outfile << "#OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "LinearStatic nsteps 3 nmodules 1" << endl;
	outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "domain 3d" << endl;
	outfile << "OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "ndofman " << pow((pow(2,octree.level_)+1),3)<< " nelem " << voxel_number << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 " << endl;
	outfile.close();
	cout << voxel_number;
	coor(vSize.x, vSize.y, vSize.z, octree.level_);





 	vector<ca_envolution> auto_cell(pow(pow(2, octree.level_), 3));
	vector<int> ca(voxel_number);



	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suObejctOctree<NodeData>::OctreeNode *pNode = octree.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);//体素中心的坐标，能够据此计算出编号
				if (NULL != pNode)
				{
					// 内部节点输出 1 
					if (INTERIOR_CELL == pNode->label_)
					{
						vtk << 1 << " ";
						fstream outfile;
						outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY<< "  " << IndexZ;
						outfile << "LSpace " << ++count << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ,octree.level_);
						int temp = merton(IndexX, IndexY, IndexZ, octree.level_);
						//auto_cell[temp].x = IndexX;
						//auto_cell[temp].y = IndexY;
						//auto_cell[temp].z = IndexZ;
						auto_cell[temp].level = octree.level_;
						auto_cell[temp].merton = temp;
						auto_cell[temp].lable = 0;
						ca[count - 1] = temp;
						//voxel_output asd(pChildNode, level);
						outfile.close();
					
						
						asd->output_point();
						delete asd;		

					}
					// 边界节点输出 2
					else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
					{
						vtk << 2 << " ";
						fstream outfile;
						outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
						outfile << "LSpace " << ++count << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, octree.level_);
						int temp = merton(IndexX, IndexY, IndexZ, octree.level_);
						//auto_cell[temp].x = IndexX;
						//auto_cell[temp].y = IndexY;
						//auto_cell[temp].z = IndexZ;
						auto_cell[temp].level = octree.level_;
						auto_cell[temp].merton = temp;
						auto_cell[temp].lable = 1;
						auto_cell[temp].location = 1;
						ca[count - 1] = temp;
						outfile.close();
				
						//voxel_output asd(pChildNode, level);
						asd->output_point();
						delete asd;
					}
					// 外部节点输出 0
					else if (EXTERIOR_CELL == pNode->label_)
						vtk << 0 << " ";

					// 其他情况输出 3    这种情况应该是不出现的
					else
						vtk << 3 << " ";
				}
				else
					vtk << 4 << " ";
			}
			vtk << endl;
		}
		vtk << endl;
	}
	//fstream outfile;
	outfile.open("d://test.txt", ios::app);
	outfile << "SimpleCS 1" << endl << "IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0" << endl << "BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.0" << endl << "BoundaryCondition  2 loadTimeFunction 1 prescribedvalue 0.5" << endl << "PiecewiseLinFunction 1 npoints 2 t 2 0. 1000. f(t) 2 0. 1000." << endl;
	outfile.close();

	/*for (int i = 0; i < count; i++)
	{
		cout << ca[i] << " "<<i<<" ";
	}*/



	
	/*cout<<endl;
	cout << auto_cell[3882].location << endl;
	cout << auto_cell[3874].location << endl;
	cout << auto_cell[3889].location << endl;
	cout << auto_cell[3873].location << endl;
	cout << auto_cell[3879].location << endl;
	cout << auto_cell[3847].location << endl;
	*/
	int add = 2;
	for (;;)
	{
		int location_temp = 0;
		for (int i = 0; i < count; i++)//遍历所有体素
		{
			if (auto_cell[ca[i]].location == 0)//判断本身是否为0
			{

				int *n_merton = six_n_merton(ca[i], octree.level_);//输入莫顿序  输出6-领域莫顿序    六位数组
				/*int compare = 0;
				for (int j = 0; j < 6; j++)//计算出六领域中loaction最大的
				{
					if (auto_cell[n_merton[j]].location > compare);
					compare = auto_cell[n_merton[j]].location;
				}
				if (compare > 0)//location是否赋值  是则赋值   计数加1
				{
					auto_cell[ca[i]].location_temp = compare + 1;
					location_temp++;
					//cout << ca[i] << " " << auto_cell[ca[i]].location_temp << " " << auto_cell[ca[i]].location << "       ";
				}
				else
					cout << auto_cell[n_merton[0]].location << auto_cell[n_merton[1]].location << auto_cell[n_merton[2]].location << auto_cell[n_merton[3]].location << auto_cell[n_merton[4]].location << auto_cell[n_merton[5]].location;
				*/
				if ((auto_cell[n_merton[0]].location + auto_cell[n_merton[1]].location + auto_cell[n_merton[2]].location + auto_cell[n_merton[3]].location + auto_cell[n_merton[4]].location + auto_cell[n_merton[5]].location) != 0)
				{
					auto_cell[ca[i]].location_temp = add;
					
				}
			}
		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].location_temp == add)
			{
				auto_cell[ca[i]].location = add;
				//cout << auto_cell[ca[i]].location << " " ;
				auto_cell[ca[i]].location_temp = 0;
				location_temp++;
			}
		}
		//cout << "location_temp="<<location_temp << endl;
		if (location_temp == 0)
			break;
		add++;
	}
	
	int void_voxel = 0;
	float porosity[6];
	cout << endl;
	cout << void_voxel << endl;

	for (int k = 1; k <= 6; k++)
	{
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].lable == 0 && auto_cell[ca[i]].out == 1)
			{
				int *n_merton = six_n_merton(ca[i], octree.level_);
				/*if ((auto_cell[ca[i]].location >= auto_cell[ca[0]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[1]].location) &&
					(auto_cell[ca[i]].location >= auto_cell[ca[2]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[3]].location) &&
					(auto_cell[ca[i]].location >= auto_cell[ca[4]].location) && (auto_cell[ca[i]].location >= auto_cell[ca[5]].location))
					{
					auto_cell[ca[i]].out_temp = 0;
					}*/

				int max = 0;
				for (int j = 0; j < 6; j++)
				{
					if ((auto_cell[n_merton[j]].location>max) && (auto_cell[n_merton[j]].out == 1))
						max = auto_cell[n_merton[j]].location;
				}
				delete n_merton;
				if (auto_cell[ca[i]].location >= max)
				{
					auto_cell[ca[i]].out_temp = 0;
					
				}//cout << max << "  " << auto_cell[ca[i]].location << "  " << auto_cell[ca[i]].out_temp << "||";
			}
		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out_temp == 0)
			{
				auto_cell[ca[i]].out = 0;
				auto_cell[ca[i]].out_temp = 1;
				
			}

		}
		for (int i = 0; i < count; i++)
		{
			if (auto_cell[ca[i]].out == 0)
			{
				void_voxel++;

			}

		}
		porosity[k - 1] = void_voxel / count;
		cout << k << ":" << void_voxel << ";" <<count << endl;
		void_voxel = 0;
	}




	



	int voxel_number1 = 0;


	for (int i = 0; i < count; i++)
	{
		if (auto_cell[ca[i]].out == 1)
		{
			voxel_number1++;
			//fstream outfile;
			//outfile.open("d://test1.txt", ios::app);
			//int *asddd = code(ca[i], octree.level_);
			//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
			//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
			//outfile << "LSpace " << ++fuck << "	 nodes  8 ";
			//voxel_output* asd = new voxel_output(asddd[0], asddd[1], asddd[2], octree.level_);
			//outfile.close();
			//voxel_output asd(pChildNode, level);
			//asd->output_point1();
			//delete asd;
		}
	}

	outfile.open("d://test1.txt", ios::out);
	outfile << "Majnun.out" << endl;
	outfile << "test of Brick elements with nlgeo 1(strain is the Green-Lagrangian strain) rotated as a rigid body" << endl;
	outfile << "#NonLinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#LinearStatic  nmsteps 1 nsteps 1 " << endl;
	outfile << "#nsteps 5 rtolv 1.e-6 stiffMode 1 controlmode 1 maxiter 100" << endl;
	outfile << "#vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "#domain 3d" << endl;
	outfile << "#OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "LinearStatic nsteps 3 nmodules 1" << endl;
	outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << endl;
	outfile << "domain 3d" << endl;
	outfile << "OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "ndofman " << pow((pow(2, octree.level_) + 1), 3) << " nelem " << voxel_number1 << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 " << endl;
	outfile.close();
	coor1(vSize.x, vSize.y, vSize.z, octree.level_);




	int fuck = 0;
	for (int i = 0; i < count; i++)
	{
		if (auto_cell[ca[i]].out==1 )
		{
			fstream outfile;
			outfile.open("d://test1.txt", ios::app);
			int *asddd = code(ca[i], octree.level_);
			//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
			//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
			outfile << "LSpace " << ++fuck << "	 nodes  8 ";
			voxel_output* asd = new voxel_output(asddd[0],asddd[1], asddd[2], octree.level_);
			outfile.close();
			//voxel_output asd(pChildNode, level);
			asd->output_point1();
			delete asd;
		}
	}

	

	outfile.open("d://test1.txt", ios::app);
	outfile << "SimpleCS 1" << endl << "IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0" << endl << "BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.0" << endl << "BoundaryCondition  2 loadTimeFunction 1 prescribedvalue 0.5" << endl << "PiecewiseLinFunction 1 npoints 2 t 2 0. 1000. f(t) 2 0. 1000." << endl;
	outfile.close();

















	outfile.open("d://test2.txt", ios::out);
	outfile << "Majnun.out" << endl;
	outfile << "Simple bending of a cantilever beam, quadratic elements." << endl;
	outfile << "LinearStatic nsteps 2 controllmode 1 rtolv 1.e-3 nmodules 1" << endl;
	outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 6 1 2 4 5 27 28 stype 1" << endl;
	outfile << "domain 3d" << endl;
	outfile << "OutputManager tstep_all dofman_all element_all" << endl;
	outfile << "ndofman " << pow((pow(2, octree.level_) + 1), 3)<<" nelem "<<voxel_number1*5<<" ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1" << endl;
	outfile.close();
	coor2(vSize.x, vSize.y, vSize.z, octree.level_);




	int number = 0;//一个六面体体素生成四个四个四面体体素
	for (int i = 0; i < count; i++)
	{
		if (auto_cell[ca[i]].out == 1)
		{
			fstream outfile;
			//outfile.open("d://test1.txt", ios::app);
			int *asddd = code(ca[i], octree.level_);
			//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
			//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
			//outfile << "LSpace " << ++fuck << "	 nodes  8 ";
			voxel_output* asd = new voxel_output(asddd[0], asddd[1], asddd[2], octree.level_);
			//outfile.close();
			//voxel_output asd(pChildNode, level);
			asd->output_point2(number);
			delete asd;
			number++;
		}
	}



	outfile.open("d://test2.txt", ios::app);
	outfile << "SimpleCS 1 thick 1.0 width 1.0" << endl;
	outfile << "IsoLE 1 d 0.0102 E 200.0 n 0.394  tAlpha 0.0000780" << endl;
	outfile << "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0" << endl;
	outfile << "ConstantSurfaceLoad 2 ndofs 3 loadType 2 Components 3 0.0 -84.118 0.0 loadTimeFunction 1" << endl;
	outfile << "ConstantFunction 1 f(t) 1.0" << endl;
	outfile.close();



	/*
	vector<int> shit(voxel_number);
	for (int i = 0; i < voxel_number; i++)
	{
		shit[i] = i;
	}
	for (int i = 0; i < voxel_number; i++)
	{
		cout << shit[i] << " ";
	}*/












	/*


	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suObejctOctree<NodeData>::OctreeNode *pNode = octree.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);//体素中心的坐标，能够据此计算出编号
				if (NULL != pNode)
				{
					// 内部节点输出 
					if (INTERIOR_CELL == pNode->label_)
					{

						fstream outfile;
						outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY<< "  " << IndexZ;
						outfile << "LSpace " << count++ << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, 3);
						//voxel_output asd(pChildNode, level);


						asd->output_point();
						delete asd;

					}
					// 边界节点输出
					else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
					{

						fstream outfile;
						outfile.open("d://test.txt", ios::app);
						//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
						//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
						outfile << "LSpace " << count++ << "	 nodes  8 ";
						voxel_output* asd = new voxel_output(IndexX, IndexY, IndexZ, 3);

						//voxel_output asd(pChildNode, level);
						asd->output_point();
						delete asd;
					}


				}
			}
		}
	}


	



	*/

	/*
	cout << endl;

	for (int i = 0; i <= 5; i++)
		cout << porosity[i] << "  ";

		*/






	return 0;
}

