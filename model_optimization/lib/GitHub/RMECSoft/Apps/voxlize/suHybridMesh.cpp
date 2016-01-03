#include "config.h"
#include "suHybridMesh.h"
#include "cellData.h"
#include "overlap.h"
#include <map>

int suHybridMesh::LoadMeshFile(const char *pFileName) throw (suException)
{
	if (NULL == pFileName)
		throw suException("File is Null!");

	if (MESHSTATE_NODEFINE != meshState_)
		throw suException("Mesh is already loaded!");

	if (!Utility::FileExists(String(pFileName)))
		throw suException("Can't find " + std::string(pFileName));

	if (OpenMesh::IO::read_mesh(mesh_, pFileName))
	{
		meshState_ = MESHSTATE_READ;

		//Get bounding box
		Mesh::VertexIter  v_it(mesh_.vertices_begin());
		Mesh::VertexIter  v_end(mesh_.vertices_end());

		bbMin_ = bbMax_ = mesh_.point(v_it);
		for (; v_it != v_end; ++v_it)
		{
			bbMin_.minimize(mesh_.point(v_it));
			bbMax_.maximize(mesh_.point(v_it));
		}
		bbox_.Min = Vec3f(bbMin_[0], bbMin_[1], bbMin_[2]);
		bbox_.Max = Vec3f(bbMax_[0], bbMax_[1], bbMax_[2]);

		length_ = bbMax_[0] - bbMin_[0];
		width_ = bbMax_[1] - bbMin_[1];
		height_ = bbMax_[2] - bbMin_[2];

		std::cout << "length: " << length_ << std::endl;
		std::cout << "width : " << width_ << std::endl;
		std::cout << "height: " << height_ << std::endl;

		//Create a mesh object (throw mesh_query lib)
		double *vertex = new double[3 * mesh_.n_vertices()];
		int *trig = new int[3 * mesh_.n_faces()];

		int idx = 0;
		for (Mesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it){
			for (int jj = 0; jj < 3; jj++){
				vertex[idx++] = mesh_.point(*v_it)[jj];
			}
		}
		idx = 0;
		for (Mesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it)
		{
			Mesh::FaceVertexIter v_it = mesh_.fv_begin(f_it.handle());
			for (; v_it != mesh_.fv_end(f_it.handle()); ++v_it)
			{
				trig[idx++] = v_it.handle().idx();
			}
		}
		pMeshObject_ = construct_mesh_object(mesh_.n_vertices(), vertex, mesh_.n_faces(), trig);

		return 0;
	}
	else
		throw suException("Mesh read error!");
}

int suHybridMesh::SaveMeshFile(const char *pFileName)
{
	return 0;
}

int suHybridMesh::PartitionSpace(unsigned int level) throw (suException)
{
	if (meshState_ == MESHSTATE_NODEFINE)
		throw(suException("No mesh loaded!"));

	level_ = level;

	if (level_ < 0) throw suException("octave tree's level less than 0!");

	if (0 != octree_.SetBox(bbMin_[0], bbMin_[1], bbMin_[2], bbMax_[0], bbMax_[1], bbMax_[2]))
		throw suException("octree setbox() error!");

	if (0 != octree_.Generate(level))
		throw suException("octree Generate() error!");

	float _nCount = (float)pow(2, level_);
	dx_ = length_ / _nCount;
	dy_ = width_ / _nCount;
	dz_ = height_ / _nCount;

	std::cout << "dx: " << dx_ << std::endl;
	std::cout << "dy: " << dz_ << std::endl;
	std::cout << "dz: " << dy_ << std::endl;

	FillCell();
	return 0;
}

int suHybridMesh::FillCell() throw (suException)
{
	if (octree_.GetRoot() == NULL)
		throw suException("Can't get octree root node!");

	// take each vertex into the corresponding node.
	Mesh::VertexIter   v_it(mesh_.vertices_begin());
	Mesh::VertexIter   v_end(mesh_.vertices_end());

	for (; v_it != v_end; ++v_it)
	{		
		Mesh::Point   point = mesh_.point(v_it);

		suCell*  pCurrentCell = octree_.GetTreeNode(point[0], point[1], point[2]);
		if (pCurrentCell != NULL)
		{			
			if (pCurrentCell->data == NULL)
			{
				CellDataVertex *pCellVertex = new CellDataVertex;
				pCellVertex->AddElement(v_it);
				pCurrentCell->data = (void *)pCellVertex;
			}
			else
			{
				CellDataVertex *pSlotVector = (CellDataVertex*)(pCurrentCell->data);
				pSlotVector->AddElement(v_it);
			}
		}
	}

	return 0;
}

int suHybridMesh::GenLeafNodeVector() throw (suException)
{
	if (meshState_ == MESHSTATE_NODEFINE)
		throw suException("No mesh loaded!");

	// 情新获取
	LeafCellVector.clear();

	suCell *pCurrentNode = octree_.GetRoot();
	if (NULL == pCurrentNode)
		throw suException("No octree exist!");

	// 按照层级访问八叉树
	queue<suCell*> CellQueue;
	CellQueue.push(pCurrentNode);

	while (!CellQueue.empty())
	{
		// 获取队列中的头节点，并将该节点从队列中剔除
		pCurrentNode = CellQueue.front();
		CellQueue.pop();

		if (!pCurrentNode) continue;

		// 将叶子节点保存到动态数组中
		if (pCurrentNode->pChildren != NULL)
		{
			for (unsigned int i = 0; i < 8; i++)
			{
				suCell* pChildNode = pCurrentNode->pChildren + i;
				CellQueue.push(pChildNode);
			}
		}
		else
			LeafCellVector.push_back(pCurrentNode);
	}

	return LeafCellVector.size();
}

/*!
* \function  LabelBoundaryNode()
*
* \brief  Label all all boundary/special-boundary nodes 
*         and store them into a vector.
*
* \return @Number of boundary nodes.
* \author  Yuan Yao
* \date
*/
size_t suHybridMesh::LabelBoundaryNode() throw (suException)
{	
	BoundaryNodeVector.clear();

	//Label boundary nodes
	for (unsigned int i = 0; i < LeafCellVector.size(); i++){
		suCell *pLeafNode = LeafCellVector[i];
		if (pLeafNode->data != NULL){
			pLeafNode->label = BOUNDARY_CELL;
			BoundaryNodeVector.push_back(pLeafNode);
		}

	}

	// Label special-boundary node which have triangles inside.
	unsigned int isFound = 0;
	do {
		// Labeling until no new special-boundary node is founded.
		isFound = 0;

		// label node intersect with triangle
		for (unsigned int i = 0; i < LeafCellVector.size(); i++)
		{			
			suCell *pLeafNode = LeafCellVector[i];
			if (pLeafNode->label == BOUNDARY_CELL)
				continue;
						
			if (labelBoundaryNode(pLeafNode)){
				isFound = 1;
				BoundaryNodeVector.push_back(pLeafNode);			
			}
		}
	} while (isFound);

	return BoundaryNodeVector.size();
}

/*!
* \function  labelBoundaryNode(suCell* pBoundaryNode)
*
* \brief  Label a special-boundary nodes
*         and store corresponding triangles.
* \param  @pNode  a undefined node or a special-boundary node
* \return @true: if a special-boundary nodes is found.
* \note   This function is iteratively invoked until all corresponding triangles are found.
*         
*          
* \author  Wenyu Jin
* \change by Yuan Yao
* \date   
*/
int suHybridMesh::labelBoundaryNode(suCell* pNode)
{
	vector<suCell*>  NeighbourCell;
	if (octree_.Get6Neighbour(pNode, NeighbourCell) < 0)
		return 0;

	/// Store all triangles that (maybe) intersect with this node.
	//All related face to current node.	
	
	//保存当前节点已有的相交三角形，这里我使用了一个map[面索引,面句柄]的数据结构来存储。
	map<int, Mesh::FaceHandle> FaceHandleMap;
	//保存本次新加入的相交三角形
	vector<Mesh::FaceHandle> FaceHandleVector;
	
	
	CellDataFace* pAllSlotFace = NULL;
	if (pNode->label == BOUNDARY_CELL_SPECIAL)
	{
		
		pAllSlotFace = (CellDataFace*)(pNode->data);   //If a special-boundary's data is exist, 
	                                                   //try to find new triangles intersect with the node.
		assert(pAllSlotFace);
		for (unsigned int i = 0; i < pAllSlotFace->GetSize();i++)
		{
			FaceHandleMap[pAllSlotFace->FacesVector[i].idx()] = pAllSlotFace->FacesVector[i];
		}
	}	

	// Check the nearest domain
	// add faces in this domain into FaceHandleVector.
	for (unsigned int i = 0; i < NeighbourCell.size(); i++)
	{		
		suCell* pLeafNode = NeighbourCell[i];
		if (pLeafNode->label == UNDEFINE_CELL)
			continue;
				
		if (pLeafNode->label == BOUNDARY_CELL)
		{
			//Add faces binding of the nearest boundary node into FaceHandleVector.
			CellDataVertex *pSlotVector = (CellDataVertex*)(pLeafNode->data);
			for (unsigned int j = 0; j < pSlotVector->GetSize(); j++)
			{	
				Mesh::VertexIter CurrentVertexIter = pSlotVector->VertexVector[j];
				Mesh::VertexFaceIter vf_it = mesh_.vf_iter(CurrentVertexIter.handle());

				for (; vf_it; ++vf_it)
				{					
					if (FaceHandleMap.find(vf_it.handle().idx()) == FaceHandleMap.end())
					{
						FaceHandleVector.push_back(vf_it.handle() );
					}
				}
			}
		}
		else if (pLeafNode->label == BOUNDARY_CELL_SPECIAL){
			//Add faces binding the nearest special-boundary node.
			CellDataFace* pSlotFace = (CellDataFace*)pLeafNode->data;
			
			assert(pSlotFace);

			for (unsigned int j = 0; j < pSlotFace->GetSize(); j++)
			{
				if (FaceHandleMap.find(pSlotFace->FacesVector[j].idx() ) == FaceHandleMap.end())
				{
					FaceHandleVector.push_back(pSlotFace->FacesVector[j]);
				}
			}
		}
	}
	
	//no intersect triangles
	if (FaceHandleVector.empty())
		return 0;

	// 获取当前节点的中心
	float boxcenter[3];
	octree_.GetLocation(pNode, boxcenter[0], boxcenter[1], boxcenter[2]);

	// 获取当前节点长度的一半
	float boxhalfsize[3];
	float minX, minY, minZ;
	float maxX, maxY, maxZ;
	octree_.GetBox(minX, minY, minZ, maxX, maxY, maxZ);
	boxhalfsize[0] = (maxX - minX) / (float)((unsigned int)1 << (pNode->level + 1));
	boxhalfsize[1] = (maxY - minY) / (float)((unsigned int)1 << (pNode->level + 1));
	boxhalfsize[2] = (maxZ - minZ) / (float)((unsigned int)1 << (pNode->level + 1));

	// 获取三角面的三个顶点坐标
	float triverts[3][3];

	// 标记是否往当前节点中添加了新的三角面
	unsigned int isAdd = 0;
	for (unsigned int ii = 0; ii < FaceHandleVector.size();ii++)
	{
		// 获取当前三角面
		Mesh::FaceHandle CurrentFaceHandle = FaceHandleVector[ii];		
		Mesh::FaceVertexIter fv_it = mesh_.fv_iter(CurrentFaceHandle);
		// 顶点1
		Mesh::Point    point = mesh_.point(fv_it);
		triverts[0][0] = point[0];
		triverts[0][1] = point[1];
		triverts[0][2] = point[2];
		// 顶点2
		point = mesh_.point(++fv_it);
		triverts[1][0] = point[0];
		triverts[1][1] = point[1];
		triverts[1][2] = point[2];
		// 顶点3
		point = mesh_.point(++fv_it);
		triverts[2][0] = point[0];
		triverts[2][1] = point[1];
		triverts[2][2] = point[2];

		// 三角面与体素之间存在重叠
		if (1 == TriBoxOverlap(boxcenter, boxhalfsize, triverts))
		{
			// 当前节点被判定为特殊表面节点
			pNode->label = BOUNDARY_CELL_SPECIAL;

			// 当前节点中将会添加没有处理的三角面
			isAdd = 1;

			// 将三角面添加到节点数据区
			if (pNode->data == NULL)
			{
				CellDataFace*  pSlotFace = new CellDataFace;
				if (NULL == pSlotFace)
					return -1;

				pSlotFace->AddElement(CurrentFaceHandle);

				pNode->data = (void *)pSlotFace;
			}
			else
			{
				CellDataFace*  pSlotFace = (CellDataFace*)pNode->data;

				bool _bRepeat = false;
				for (unsigned int ii = 0; ii < pSlotFace->FacesVector.size(); ii++)
				{
					if (pSlotFace->FacesVector[ii].idx() == CurrentFaceHandle.idx())
					{
						_bRepeat = true;
					}
				}
				if (!_bRepeat)
				{
					pSlotFace->AddElement(CurrentFaceHandle);
				}
				
			}
		}
	}

	if (!isAdd)
		return 0;	
	else
		return 1;
}

int suHybridMesh::SaveAsVTK(const char *pVTKFileName, const char *pVTKHead /*= "UnKnown Name"*/, float dx/*=0*/, float dy/*=0*/, float dz/*=0*/)
{
	if (NULL == pVTKFileName || NULL == pVTKHead)
		return -1;

	float FLOAT_ZERO = 0.0000001f;

	dx = (dx < FLOAT_ZERO) ? dx_ : dx;
	dy = (dy < FLOAT_ZERO) ? dy_ : dy;
	dz = (dz < FLOAT_ZERO) ? dz_ : dz;
	ofstream vtk(pVTKFileName);
	//stringstream vtk;

	if (!vtk) throw suException(std::string("Can't open ") + std::string(pVTKFileName));

	// VTK File Header
	vtk << "# vtk DataFile Version 2.0" << "\n";
	vtk << pVTKHead << "\n";
	vtk << "ASCII" << "\n";
	vtk << "DATASET STRUCTURED_POINTS" << "\n";

	// VTK File Body
	float minx, miny, minz;
	float maxx, maxy, maxz;
	octree_.GetBox(minx, miny, minz, maxx, maxy, maxz);

	unsigned int Xdim = (unsigned int)(length_ / dx);
	unsigned int Ydim = (unsigned int)(width_ / dy);
	unsigned int Zdim = (unsigned int)(height_ / dz);

	vtk << "DIMENSIONS " << Xdim << " " << Ydim << " " << Zdim << "\n";

	// 数据区的其他信息
	vtk << "ASPECT_RATIO 1 1 1" << "\n";
	vtk << "ORIGIN 0 0 0" << "\n";
	vtk << "POINT_DATA " << Xdim * Ydim * Zdim << "\n";
	vtk << "SCALARS volume_scalars char 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n" << "\n";

	for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
	{
		float CellCenter_Z = minz + IndexZ * dz + dz / 2.0f;

		for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
		{
			float CellCenter_Y = miny + IndexY * dy + dy / 2.0f;

			for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
			{
				float CellCenter_X = minx + IndexX * dx + dx / 2.0f;

				suCell *pCell = octree_.GetTreeNode(CellCenter_X, CellCenter_Y, CellCenter_Z);

				if (NULL != pCell)
				{
					// 内部节点输出 
					if (INTERIOR_CELL == pCell->label)
						vtk << INTERIOR_CELL << " ";

					// 边界节点输出
					else if (BOUNDARY_CELL == pCell->label || BOUNDARY_CELL_SPECIAL == pCell->label)
						vtk << BOUNDARY_CELL << " ";

					// 外部节点输出
					else if (EXTERIOR_CELL == pCell->label)
						vtk << EXTERIOR_CELL << " ";

					// debug:其他情况输出,这种情况应该是不出现的
					else
						vtk << UNDEFINE_CELL << " ";
				}
			}
			vtk << "\n";
		}
		vtk << "\n";
	}

	return 0;
}

void suHybridMesh::Clear()
{
	meshState_ = MESHSTATE_NODEFINE;

	length_ = 0;
	width_ = 0;
	height_ = 0;
	dx_ = 0;;
	dy_ = 0;
	dz_ = 0;
	level_ = 0;

	pMeshObject_ = NULL;

	BoundaryNodeVector.clear();
	ExteriorNodeVector.clear();
	InteriorNodeVector.clear();
	LeafCellVector.clear();
}
/*!
 * \function
 *
 * \brief Label in/out nodes near the boundary nodes.
 *
 * \return @Number of exterior nodes.
 * \author  Yuan Yao
 * \date
 */
int suHybridMesh::LabelInAndOutNode()
{
	ExteriorNodeVector.clear();
	InteriorNodeVector.clear();

	for (size_t i = 0; i < BoundaryNodeVector.size(); i++)
	{
		suCell* pBoundaryNode = BoundaryNodeVector[i];

		//Get 6-neighbor nodes
		vector<suCell*>  nnNode;
		if (octree_.Get6Neighbour(pBoundaryNode, nnNode) == 0)
			return -3;

		//Set node type for 6-neighbor nodes.
		for (int j = 0; j < (int)nnNode.size(); j++)
		{
			if (nnNode[j]->label == UNDEFINE_CELL)
			{

				//if (insideMesh(nnNode[j], 0, 0, 1))
				if (insideMeshByQeury(nnNode[j]) )
				{
					nnNode[j]->label = INTERIOR_CELL;
					InteriorNodeVector.push_back(nnNode[j]);
				}
				else
				{
					nnNode[j]->label = EXTERIOR_CELL;
					ExteriorNodeVector.push_back(nnNode[j]);
				}
			}
		}

	}
	return ExteriorNodeVector.size() + InteriorNodeVector.size();
}


/**
*\brief  使用区域增长标记所有未定义节点.
*
*\return  1
*\
*/
int suHybridMesh::ExtendLabelNode()
{
	for (unsigned int i = 0; i < ExteriorNodeVector.size(); i++)
	{
		extendLabelNode(ExteriorNodeVector[i]);
	}
	for (unsigned int i = 0; i < InteriorNodeVector.size(); i++)
	{
		extendLabelNode(InteriorNodeVector[i]);
	}

	return 1;
}
/**
*\brief  通过递归增长，标记邻域所有相邻节点。
*
*\return  1
*\
*/
void suHybridMesh::extendLabelNode(suCell* pNode)
{
	vector<suCell*>  NeighbourCell;
	if (octree_.Get6Neighbour(pNode, NeighbourCell) == 0) return;

	for (int j = 0; j < (int)NeighbourCell.size(); j++)
	{
		if (NeighbourCell[j]->label == UNDEFINE_CELL)
		{
			NeighbourCell[j]->label = pNode->label;
			extendLabelNode(NeighbourCell[j]);
		}
	}
}
/**
*\brief  使用mesh_query判断节点中心是否位于mesh中
*
*\return  1
*\
*/
bool suHybridMesh::insideMeshByQeury(suCell* pNode)
{
	float x, y, z;
	if (octree_.GetLocation(pNode, x, y, z) != 0) return false;
		
	Vec3f pCenter = Vec3f(x, y, z);

	double center[3];
	center[0] = x;
	center[1] = y;
	center[2] = z;

	return point_inside_mesh(center, pMeshObject_);
}

suHybridMesh::~suHybridMesh()
{
	if (pMeshObject_ != NULL)
	{
		destroy_mesh_object(pMeshObject_);
	}
}
