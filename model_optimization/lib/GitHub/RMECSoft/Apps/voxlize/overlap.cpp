/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Moller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/*                                                      */
/* 2010-02-16   金文玉将此代码修改，然后应用到自己的任务*/
/*              修改包括：编码方式、变量命名、必要的注释*/
/*                                                      */
/*                                                      */
/********************************************************/
#include "overlap.h"
//------------------------------------------------------------------------------------------------------

#define OVERLAP_CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 

//------------------------------------------------------------------------------------------------------

#define OVERLAP_DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

//------------------------------------------------------------------------------------------------------

#define OVERLAP_SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 

//------------------------------------------------------------------------------------------------------

#define OVERLAP_FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;
//------------------------------------------------------------------------------------------------------

static int planeBoxOverlap(float normal[3],float d, float maxbox[3])
{
	int   q;
	float vmin[3],vmax[3];
  
	// 根据向量<normal[0], normal[1], normal[2]>确定BOX中最高顶点和最顶点
	for (q = OVERLAP_X; q <= OVERLAP_Z; q++)
	{
		if (normal[q] > 0.0f)
		{
			vmin[q] = -maxbox[q];
			vmax[q] = maxbox[q];
		} else
		{
			vmin[q] = maxbox[q];
			vmax[q] = -maxbox[q];
		}
	}
	
	// 相对最低点在平面外
	if (OVERLAP_DOT(normal,vmin) + d > 0.0f) 
		return 0;
	
	// 相对最高点在平面外或者上
	if (OVERLAP_DOT(normal,vmax) + d >= 0.0f) 
		return 1;
  
	return 0;
}
//------------------------------------------------------------------------------------------------------

// X-tests
#define OVERLAP_AXISTEST_X01(a, b, fa, fb)			                   \
	p0 = a*v0[OVERLAP_Y] - b*v0[OVERLAP_Z];			       	           \
	p2 = a*v2[OVERLAP_Y] - b*v2[OVERLAP_Z];			       	           \
    if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}      \
	rad = fa * boxhalfsize[OVERLAP_Y] + fb * boxhalfsize[OVERLAP_Z];   \
	if (min > rad || max < -rad) return 0;

#define OVERLAP_AXISTEST_X2(a, b, fa, fb)			                  \
	p0 = a * v0[OVERLAP_Y] - b * v0[OVERLAP_Z];			              \
	p1 = a * v1[OVERLAP_Y] - b * v1[OVERLAP_Z];			       	      \
    if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}     \
	rad = fa * boxhalfsize[OVERLAP_Y] + fb * boxhalfsize[OVERLAP_Z];  \
	if (min > rad || max < -rad) return 0;

//------------------------------------------------------------------------------------------------------

// Y-tests 
#define OVERLAP_AXISTEST_Y02(a, b, fa, fb)			                  \
	p0 = -a * v0[OVERLAP_X] + b * v0[OVERLAP_Z];		      	      \
	p2 = -a * v2[OVERLAP_X] + b * v2[OVERLAP_Z];	       	       	  \
    if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}     \
	rad = fa * boxhalfsize[OVERLAP_X] + fb * boxhalfsize[OVERLAP_Z];  \
	if (min > rad || max < -rad) return 0;

#define OVERLAP_AXISTEST_Y1(a, b, fa, fb)			                  \
	p0 = -a * v0[OVERLAP_X] + b * v0[OVERLAP_Z];		      	      \
	p1 = -a * v1[OVERLAP_X] + b * v1[OVERLAP_Z];       	       	      \
    if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}     \
	rad = fa * boxhalfsize[OVERLAP_X] + fb * boxhalfsize[OVERLAP_Z];  \
	if (min > rad || max < -rad) return 0;

//------------------------------------------------------------------------------------------------------

// Z-tests

#define OVERLAP_AXISTEST_Z12(a, b, fa, fb)			                  \
	p1 = a * v1[OVERLAP_X] - b * v1[OVERLAP_Y];			              \
	p2 = a * v2[OVERLAP_X] - b * v2[OVERLAP_Y];			       	      \
    if (p2 < p1) {min = p2; max = p1;} else {min = p1; max = p2;}     \
	rad = fa * boxhalfsize[OVERLAP_X] + fb * boxhalfsize[OVERLAP_Y];  \
	if (min > rad || max < -rad) return 0;

#define OVERLAP_AXISTEST_Z0(a, b, fa, fb)			                  \
	p0 = a * v0[OVERLAP_X] - b * v0[OVERLAP_Y];				          \
	p1 = a * v1[OVERLAP_X] - b * v1[OVERLAP_Y];			              \
    if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}     \
	rad = fa * boxhalfsize[OVERLAP_X] + fb * boxhalfsize[OVERLAP_Y];  \
	if (min > rad || max < -rad) return 0;

//------------------------------------------------------------------------------------------------------

// 判断三角面是否与BOX相交
// 参数
//     boxcenter     BOX的中心坐标
//     boxhalfsize   BOX在X\Y\Z三个方向长度的一半 
//     triverts      三角面的三个顶点坐标，三个顶点分别是triverts[0]\triverts[1]\triverts[2]
// 返回值
//     1             三角面与BOX存在重叠
//     0             三角面与BOX不存在重叠
// 说明
//    use separating axis theorem to test overlap between triangle and box 
//    need to test for overlap in these directions: 
//      1) the {x,y,z}-directions (actually, since we use the AABB of the triangle 
//         we do not even need to test these) 
//      2) normal of the triangle 
//      3) crossproduct(edge from tri, {x,y,z}-directin) 
//         this gives 3x3=9 more tests 
//    具体参见文章Fast 3D Triangle-Box Overlap Testing
int TriBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3])
{
	// 平移之后，三个顶点的坐标
	float v0[3], v1[3], v2[3];
	// 最大最小值
	float min, max;
	// 平面的方程的常量项
	float d;
	// 三个方向的投影值
	float p0, p1, p2;
	// 投影长度
	float rad;
	// 到中心的距离
	float fex, fey, fez;  
	// 三角面在新坐标之下的的法向量、三条边
	float normal[3], e0[3], e1[3], e2[3];

	// This is the fastest branch on Sun 
	// move everything so that the boxcenter is in (0,0,0)
	OVERLAP_SUB(v0, triverts[0], boxcenter);
	OVERLAP_SUB(v1, triverts[1], boxcenter);
	OVERLAP_SUB(v2, triverts[2], boxcenter);

	// compute triangle edges
	OVERLAP_SUB(e0,v1,v0);      /* tri edge 0 */
	OVERLAP_SUB(e1,v2,v1);      /* tri edge 1 */
	OVERLAP_SUB(e2,v0,v2);      /* tri edge 2 */

	// Bullet 3:  
	//  test the 9 tests first (this was faster)
	fex = fabs(e0[OVERLAP_X]);
	fey = fabs(e0[OVERLAP_Y]);
	fez = fabs(e0[OVERLAP_Z]);
	OVERLAP_AXISTEST_X01(e0[OVERLAP_Z], e0[OVERLAP_Y], fez, fey);
	OVERLAP_AXISTEST_Y02(e0[OVERLAP_Z], e0[OVERLAP_X], fez, fex);
	OVERLAP_AXISTEST_Z12(e0[OVERLAP_Y], e0[OVERLAP_X], fey, fex);

	fex = fabs(e1[OVERLAP_X]);
	fey = fabs(e1[OVERLAP_Y]);
	fez = fabs(e1[OVERLAP_Z]);
	OVERLAP_AXISTEST_X01(e1[OVERLAP_Z], e1[OVERLAP_Y], fez, fey);
	OVERLAP_AXISTEST_Y02(e1[OVERLAP_Z], e1[OVERLAP_X], fez, fex);
	OVERLAP_AXISTEST_Z0(e1[OVERLAP_Y], e1[OVERLAP_X], fey, fex);

	fex = fabs(e2[OVERLAP_X]);
	fey = fabs(e2[OVERLAP_Y]);
	fez = fabs(e2[OVERLAP_Z]);
	OVERLAP_AXISTEST_X2(e2[OVERLAP_Z], e2[OVERLAP_Y], fez, fey);
	OVERLAP_AXISTEST_Y1(e2[OVERLAP_Z], e2[OVERLAP_X], fez, fex);
	OVERLAP_AXISTEST_Z12(e2[OVERLAP_Y], e2[OVERLAP_X], fey, fex);

	// Bullet 1: 
	//  first test overlap in the {x,y,z}-directions
	//  find min, max of the triangle each direction, and test for overlap in
	//  that direction -- this is equivalent to testing a minimal AABB around 
	//  the triangle against the AABB 	

	// test in X-direction 
	OVERLAP_FINDMINMAX(v0[OVERLAP_X],v1[OVERLAP_X],v2[OVERLAP_X],min,max);
	if (min > boxhalfsize[OVERLAP_X] || max < -boxhalfsize[OVERLAP_X]) 
		return 0;

	// test in Y-direction
	OVERLAP_FINDMINMAX(v0[OVERLAP_Y],v1[OVERLAP_Y],v2[OVERLAP_Y],min,max);
	if (min > boxhalfsize[OVERLAP_Y] || max < -boxhalfsize[OVERLAP_Y]) 
		return 0;

	// test in Z-direction
	OVERLAP_FINDMINMAX(v0[OVERLAP_Z],v1[OVERLAP_Z],v2[OVERLAP_Z],min,max);
	if (min > boxhalfsize[OVERLAP_Z] || max < -boxhalfsize[OVERLAP_Z]) 
		return 0;

	// Bullet 2:
	//  test if the box intersects the plane of the triangle
	//  compute plane equation of triangle: normal*x+d=0
	OVERLAP_CROSS(normal,e0,e1);
	d = -OVERLAP_DOT(normal,v0);  /* plane eq: normal.x+d=0 */
	if (0 == planeBoxOverlap(normal, d, boxhalfsize)) 
		return 0;

	// box and triangle overlaps
	return 1;
}
//------------------------------------------------------------------------------------------------------