// UnitTest.cpp : 定义控制台应用程序的入口点。
// Yuan Yao,
// 2012/12/12
// 单元测试例程。
// 每个UTFUNC模块定义了一个suUnitTestFunction类的继承类和相应的test()函数。
// 并将自身加入到suUnitTestFunction的成员变量中维护。

#include "stdafx.h"
#include <engine.h>


//#include <fstream>
//#include <sstream>


//测试std::vector的例程
//参考： http://www.cplusplus.com/reference/stl/vector/
UTFUNC(test_std_vector)
{
	std::vector<int> arr;
	for (int i=0; i<10; i++)
	{
		arr.push_back(i);
	}

	for (int i=0; i<10; i++)
	{
		std::cout << "arr[" << i << "]=" << arr[i] << std::endl;
	}
	VERIFY(true);
}

//测试std::list的例程
//参考 http://www.cplusplus.com/reference/stl/list/
UTFUNC(test_std_list)
{
	std::list<int> mylist;
	std::list<int>::iterator it;

	// set some initial values:
	for (int i=1; i<=5; i++) mylist.push_back(i); // 1 2 3 4 5

	it = mylist.begin();
	++it;       // it points now to number 2           ^

	mylist.insert (it,10);                        // 1 10 2 3 4 5

	// "it" still points to number 2                      ^
	mylist.insert (it,2,20);                      // 1 10 20 20 2 3 4 5

	--it;       // it points now to the second 20            ^

	std::vector<int> myvector (2,30);
	mylist.insert (it,myvector.begin(),myvector.end());
	// 1 10 20 30 30 20 2 3 4 5
	//               ^
	std::cout << "mylist contains:";
	for (it=mylist.begin(); it!=mylist.end(); it++)
		std::cout << " " << *it;
	std::cout << std::endl;

	VERIFY(true);
}

//测试std::map的例程
//参考 http://www.cplusplus.com/reference/stl/map/
UTFUNC(test_std_map)
{
	std::map<std::string,int> mymap;
	std::map<std::string,int>::iterator it;

	mymap["first"]=50;
	mymap["second"]=100;
	mymap["third"]=150;
	mymap["forth"]=200;

	it=mymap.find("second");
	mymap.erase (it);
	mymap.erase (mymap.find("forth"));

	// print content:
	std::cout << "elements in mymap:" << std::endl;
	std::cout << "first => " << mymap.find("first")->second << std::endl;
	std::cout << "third => " << mymap.find("third")->second << std::endl;

	VERIFY(true);
}

/** 
*定义一个简单的三维点对象
*/

//class Vec3f{
//public:
//	Vec3f(){x_=0; y_=0; z_=0;}
//	Vec3f(float x, float y, float z){x_=x; y_=y; z_=z;}
//	bool operator ==(const Vec3f& rhs)    //重载==运算符
//	{
//		return (this->x_ == rhs.x_) && (this->y_ == rhs.y_) && (this->z_ == rhs.z_);
//	}
//	Vec3f& operator +=(const Vec3f& rhs)
//	{
//		this->x_ += rhs.x_;
//		this->y_ += rhs.y_;
//		this->z_ += rhs.z_;
//		return *this;
//	}
//public:
//	float x_;
//	float y_;
//	float z_;
//};

/** 
*定义个简单的六面体对象
*/
class Hexahedral{
public:
	Hexahedral(){vertices_.resize(8);}
	int& getVertice(int idx){return vertices_[idx];}
	Hexahedral& operator = (Hexahedral&  rhs)   //重载操作符"="
	{
		this->vertices_.clear();
		this->vertices_.resize(8);
				
		for (int i=0; i<8; i++)
		{
			this->vertices_[i] = rhs.getVertice(i);
		}
	}
private:
	std::vector<int> vertices_;   //顶点id
};


//测试生成六面体网格
std::vector<Vec3f> points;   //数据结构：六面体顶点点云，在这里定义是为了使下一个测试单元可以使用其中存储的数据。
unsigned int nSize = 0;      //顶点数目
UTFUNC(Grid_Generation)
{
	const int width = 10;
	const int height = 10;
	const int thick = 10;

	nSize = (width+1)*(thick+1)*(height+1);
	
	Hexahedral hexahedrals[width][thick][height];   //数据结构：六面体数组,用三维坐标作为六面体的存取id，
	                                                //          例如位于原点的第一个六面体为hexahedrals[0][0][0]

	//初始化坐标四面体顶点坐标
	//生成一个包含width*height*thick个六面体的网格对象
	//六面体的尺寸为1*1*1
	
	//生成所有六面体空间节点，用vector存储节点数组
	for (int z=0; z<=height; z++)
	{
		for (int y=0; y<=thick; y++)
		{
			for (int x=0; x<=width; x++)
			{
				points.push_back(Vec3f(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) ) );
			}
		}
	}

	//生成六面体数组,针对六面体的顶点索引顺序见doc/六面体结构.vsd。
	std::vector<Vec3f> trans;
	for (int z=0; z<=1; z++)
	{
		for (int y=0; y<=1; y++)
		{
			for (int x=0; x<=1; x++)
			{
				trans.push_back(Vec3f(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) ));
			}
		}
	}
	
	Hexahedral h_;
	for (int z=0; z<height; z++)
	{
		for (int y=0; y<thick; y++)
		{
			for (int x=0; x<width; x++)
			{				
				for (int i=0; i<8; i++)
				{
					//从点云中查找当前顶点的索引
					//带入太多的循环，需要优化
					int idx = 0;
					Vec3f v = Vec3f(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) );
					v += trans[i];
					for (;idx<(int)points.size(); idx++)
					{
						if (v == points[idx])
						{
							break;
						}
					}

					hexahedrals[x][y][z].getVertice(i) = idx;
				}

			}
		}
	}
}

//测试写入.fly格式的文件
//参考：http://en.wikipedia.org/wiki/PLY_(file_format)

UTFUNC(Write_FLY_Format)
{

	_tString strFileName = "r:/test.ply";
	std::ofstream File(strFileName.c_str());

	unsigned int vc = nSize;

	File << "ply" << std::endl;
	File << "format ascii 1.0" << std::endl;
	File << "element vertex " << vc << std::endl;
	File << "property float x" << std::endl;
	File << "property float y" << std::endl;
	File << "property float z" << std::endl;
	File << "element face " << 0 << std::endl;
	File << "property list uchar int vertex_indices" << std::endl;
	File << "end_header" << std::endl;

	for(unsigned i = 0; i < nSize; i++)
	{
		File << points[i].x << "\t" << points[i].y << "\t" << points[i].z << std::endl;
	}
}

UTFUNC(Test_PCA)
{
	PCA<float> pca;
	DenseMatrix<float> data;
	data.LoadFromFile("PCA_data.txt");
	data = data.Transpose();
	//data.SaveToMATLAB("r:/matlab_data.txt");
	pca.InitFromPointMatrix(data);
	pca.ReducedDimension(0.9);
	
}
int _tmain(int argc, _TCHAR* argv[])
{
	bool state = suUnitTest::gOnly().run ();
	suUnitTest::gOnly().dumpResults (std::cout);

	return state;
}

