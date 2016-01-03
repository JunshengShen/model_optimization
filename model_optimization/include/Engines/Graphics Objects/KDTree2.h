/*
KDTree2.h
Written by Matthew Fisher

A 3D KD-tree that supports fast points-within-rectangle queries
*/

#ifdef USE_ANN
#ifdef USE_KDTREE
#ifdef USE_OPENCV
#include <ANN/ANN.h>

const UINT KDTree2MaxK = 20000;

class KDTree2
{
public:
    KDTree2();
    ~KDTree2();
    void FreeMemory();

    
	void BuildTree(const suVector<cv::Point> &Points);	
	void KNearest(const cv::Point &Pos, UINT k, suVector<UINT> &Result, float Epsilon);
	void KNearest(const cv::Point &Pos, UINT k, UINT *Result, float Epsilon);
	void WithinDistance(const cv::Point &Pos, float Radius, suVector<UINT> &Result) const;
	UINT Nearest(const cv::Point &Pos);
	__forceinline cv::Point GetPoint(UINT Index)
	{
		cv::Point Result;
		Result.x = static_cast<UINT>(dataPts[Index][0]);
		Result.y = static_cast<UINT>(dataPts[Index][1]);
		return Result;
	}
	/*
	void BuildTree(const suVector<Vec2i> &Points);
	void KNearest(const Vec2i &Pos, UINT k, suVector<UINT> &Result, float Epsilon);
	void KNearest(const Vec2i &Pos, UINT k, UINT *Result, float Epsilon);
	void WithinDistance(const Vec2i &Pos, float Radius, suVector<UINT> &Result) const;
	UINT Nearest(const Vec2i &Pos);
	__forceinline Vec2i GetPoint(UINT Index)
	{
	Vec2i Result;
	Result.x = static_cast<int>(dataPts[Index][0]);
	Result.y = static_cast<int>(dataPts[Index][1]);
	return Result;
	}*/

private:
    /*__forceinline Vec2i GetDataPoint(UINT PointIndex)
    {
		Vec2i Result;
		Result.x = static_cast<int>(dataPts[PointIndex][0]);
		Result.y = static_cast<int>(dataPts[PointIndex][1]);
		return Result;
    }*/
	__forceinline cv::Point GetDataPoint(UINT PointIndex)
	{
		cv::Point Result;
		Result.x = static_cast<int>(dataPts[PointIndex][0]);
		Result.y = static_cast<int>(dataPts[PointIndex][1]);
		return Result;
	}
    ANNpointArray    dataPts; // data points
    ANNpoint         queryPt; // query point
    ANNidxArray      nnIdx;   // near neighbor indices
    ANNdistArray     dists;   // near neighbor distances
    ANNkd_tree*      kdTree;  // search structure
    //Mutex          _Lock;
};
#endif
#endif
#endif