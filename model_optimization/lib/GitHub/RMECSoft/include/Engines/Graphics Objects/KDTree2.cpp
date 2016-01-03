/*
KDTree2.cpp
Written by Matthew Fisher

A 3D KD-tree that supports fast points-within-rectangle queries
*/

#ifdef USE_ANN
#ifdef USE_KDTREE
#pragma comment(lib, "ANN.lib")

#ifdef USE_OPENCV

KDTree2::KDTree2()
{
    nnIdx = NULL;
    dists = NULL;
    queryPt = NULL;
    dataPts = NULL;
    kdTree = NULL;
}

KDTree2::~KDTree2()
{
    FreeMemory();
}

void KDTree2::FreeMemory()
{
    if(nnIdx)
    {
        delete[] nnIdx;
        nnIdx = NULL;
    }
    if(dists)
    {
        delete[] dists;
        dists = NULL;
    }
    if(kdTree)
    {
        delete kdTree;
        kdTree = NULL;
    }
    if(queryPt)
    {
        annDeallocPt(queryPt);
        queryPt = NULL;
    }
    if(dataPts)
    {
        annDeallocPts(dataPts);
        dataPts = NULL;
    }
}


void KDTree2::BuildTree( const suVector<cv::Point> &Points )
{
	FreeMemory();
	UINT PointCount = Points.Length();
	//Console::WriteString(String("Building KD tree, ") + String(PointCount) + String(" points..."));
	queryPt = annAllocPt(2); // allocate query point
	dataPts = annAllocPts(PointCount, 2); // allocate data points
	nnIdx = new ANNidx[KDTree2MaxK];  // allocate near neigh indices
	dists = new ANNdist[KDTree2MaxK]; // allocate near neighbor dists
	for(UINT i = 0; i < PointCount; i++)
	{
		dataPts[i][0] = Points[i].x;
		dataPts[i][1] = Points[i].y;
	}

	kdTree = new ANNkd_tree( // build search structure
		dataPts,    // the data points
		PointCount, // number of points
		2);         // dimension of space
	//Console::WriteString(String("done.\n"));
}

UINT KDTree2::Nearest(const cv::Point &Pos)
{
    suVector<UINT> Result(1);
    KNearest(Pos, 1, Result, 0.0f);
    return Result[0];
}

void KDTree2::KNearest(const cv::Point &Pos, UINT k, suVector<UINT> &Result, float Epsilon)
{
    Assert(k <= KDTree2MaxK, "k too large");
	queryPt[0] = Pos.x;
	queryPt[1] = Pos.y;
    kdTree->annkSearch(   // search
        queryPt,          // query point
        k,                // number of near neighbors
        nnIdx,            // nearest neighbors (returned)
        dists,            // distance (returned)
        Epsilon);        // error bound

    if(Result.Length() < k)
    {
        Result.ReSize(k);
    }
    for(UINT i = 0; i < k; i++)
    {
        Result[i] = nnIdx[i];
    }
}

void KDTree2::KNearest(const cv::Point &Pos, UINT k, UINT *Result, float Epsilon)
{
    Assert(k <= KDTree2MaxK, "k too large");
	queryPt[0] = Pos.x;
	queryPt[1] = Pos.y;
    kdTree->annkSearch( // search
        queryPt,        // query point
        k,                // number of near neighbors
        nnIdx,            // nearest neighbors (returned)
        dists,            // distance (returned)
        Epsilon);        // error bound
    for(UINT i = 0; i < k; i++)
    {
        Result[i] = nnIdx[i];	
    }
}

void KDTree2::WithinDistance(const cv::Point &Pos, float Radius, suVector<UINT> &Result) const
{
    queryPt[0] = Pos.x;
	queryPt[1] = Pos.y;

    int NeighborCount = kdTree->annkFRSearch(
        queryPt,
        Radius * Radius,
        KDTree2MaxK,
        nnIdx,
        dists,
        0.0f);

    Result.ReSize(Math::Min(UINT(NeighborCount), KDTree2MaxK));
    for(UINT i = 0; int(i) < Result.Length(); i++)
    {
        Result[i] = nnIdx[i];
    }
}

//void KDTree2::BuildTree(const suVector<Vec2i> &Points)
//{
//	FreeMemory();
//	UINT PointCount = Points.Length();
//	Console::WriteString(String("Building KD tree, ") + String(PointCount) + String(" points..."));
//	queryPt = annAllocPt(2); // allocate query point
//	dataPts = annAllocPts(PointCount, 2); // allocate data points
//	nnIdx = new ANNidx[KDTree2MaxK];  // allocate near neigh indices
//	dists = new ANNdist[KDTree2MaxK]; // allocate near neighbor dists
//	for(UINT i = 0; i < PointCount; i++)
//	{
//		dataPts[i][0] = Points[i].x;
//		dataPts[i][0] = Points[i].y;
//	}
//
//	kdTree = new ANNkd_tree( // build search structure
//		dataPts,    // the data points
//		PointCount, // number of points
//		2);         // dimension of space
//	Console::WriteString(String("done.\n"));
//}

//UINT KDTree2::Nearest(const Vec2i &Pos)
//{
//	suVector<UINT> Result(1);
//	KNearest(Pos, 1, Result, 0.0f);
//	return Result[0];
//}
//
//void KDTree2::KNearest(const Vec2i &Pos, UINT k, suVector<UINT> &Result, float Epsilon)
//{
//	Assert(k <= KDTree2MaxK, "k too large");
//	queryPt[0] = Pos.x;
//	queryPt[1] = Pos.y;
//	kdTree->annkSearch(   // search
//		queryPt,          // query point
//		k,                // number of near neighbors
//		nnIdx,            // nearest neighbors (returned)
//		dists,            // distance (returned)
//		Epsilon);        // error bound
//
//	if(Result.Length() < k)
//	{
//		Result.ReSize(k);
//	}
//	for(UINT i = 0; i < k; i++)
//	{
//		Result[i] = nnIdx[i];
//	}
//}
//
//void KDTree2::KNearest(const Vec2i &Pos, UINT k, UINT *Result, float Epsilon)
//{
//	Assert(k <= KDTree2MaxK, "k too large");
//	queryPt[0] = Pos.x;
//	queryPt[1] = Pos.y;
//	kdTree->annkSearch( // search
//		queryPt,        // query point
//		k,                // number of near neighbors
//		nnIdx,            // nearest neighbors (returned)
//		dists,            // distance (returned)
//		Epsilon);        // error bound
//	for(UINT i = 0; i < k; i++)
//	{
//		Result[i] = nnIdx[i];
//	}
//}
//
//void KDTree2::WithinDistance(const Vec2i &Pos, float Radius, suVector<UINT> &Result) const
//{
//	queryPt[0] = Pos.x;
//	queryPt[1] = Pos.y;
//
//	int NeighborCount = kdTree->annkFRSearch(
//		queryPt,
//		Radius * Radius,
//		KDTree2MaxK,
//		nnIdx,
//		dists,
//		0.0f);
//
//	Result.ReSize(Math::Min(UINT(NeighborCount), KDTree2MaxK));
//	for(UINT i = 0; int(i) < Result.Length(); i++)
//	{
//		Result[i] = nnIdx[i];
//	}
//}
#endif
#endif
#endif