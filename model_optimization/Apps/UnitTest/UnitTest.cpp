// UnitTest.cpp : �������̨Ӧ�ó������ڵ㡣
// Yuan Yao,
// 2012/12/12
// ��Ԫ�������̡�
// ÿ��UTFUNCģ�鶨����һ��suUnitTestFunction��ļ̳������Ӧ��test()������
// ����������뵽suUnitTestFunction�ĳ�Ա������ά����

#include "stdafx.h"
#include <engine.h>


//#include <fstream>
//#include <sstream>


//����std::vector������
//�ο��� http://www.cplusplus.com/reference/stl/vector/
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

//����std::list������
//�ο� http://www.cplusplus.com/reference/stl/list/
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

//����std::map������
//�ο� http://www.cplusplus.com/reference/stl/map/
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
*����һ���򵥵���ά�����
*/

//class Vec3f{
//public:
//	Vec3f(){x_=0; y_=0; z_=0;}
//	Vec3f(float x, float y, float z){x_=x; y_=y; z_=z;}
//	bool operator ==(const Vec3f& rhs)    //����==�����
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
*������򵥵����������
*/
class Hexahedral{
public:
	Hexahedral(){vertices_.resize(8);}
	int& getVertice(int idx){return vertices_[idx];}
	Hexahedral& operator = (Hexahedral&  rhs)   //���ز�����"="
	{
		this->vertices_.clear();
		this->vertices_.resize(8);
				
		for (int i=0; i<8; i++)
		{
			this->vertices_[i] = rhs.getVertice(i);
		}
	}
private:
	std::vector<int> vertices_;   //����id
};


//������������������
std::vector<Vec3f> points;   //���ݽṹ�������嶥����ƣ������ﶨ����Ϊ��ʹ��һ�����Ե�Ԫ����ʹ�����д洢�����ݡ�
unsigned int nSize = 0;      //������Ŀ
UTFUNC(Grid_Generation)
{
	const int width = 10;
	const int height = 10;
	const int thick = 10;

	nSize = (width+1)*(thick+1)*(height+1);
	
	Hexahedral hexahedrals[width][thick][height];   //���ݽṹ������������,����ά������Ϊ������Ĵ�ȡid��
	                                                //          ����λ��ԭ��ĵ�һ��������Ϊhexahedrals[0][0][0]

	//��ʼ�����������嶥������
	//����һ������width*height*thick����������������
	//������ĳߴ�Ϊ1*1*1
	
	//��������������ռ�ڵ㣬��vector�洢�ڵ�����
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

	//��������������,���������Ķ�������˳���doc/������ṹ.vsd��
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
					//�ӵ����в��ҵ�ǰ���������
					//����̫���ѭ������Ҫ�Ż�
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

//����д��.fly��ʽ���ļ�
//�ο���http://en.wikipedia.org/wiki/PLY_(file_format)

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

