//类voxel_output中能够生成一个体素的八个顶点的编码   构造函数中输入目标体素的指针       output_point()函数输出八个点的坐标在TXT文件中  
//函数coor输入x,y,z的最大值与分割次数  输出所有点的坐标
#include "suOctree.h"
#include <fstream>
#include <cmath>

using namespace std;

class voxel_output
{
private:
	int x_code;
	int y_code;
	int z_code;
	double point_code[8];
	int level;

public:
	voxel_output(suObejctOctree<NodeData>::OctreeNode* current) //将目标体素中的参数对应赋值
	{
		x_code = current->xLocCode_ + 1;
		y_code = current->yLocCode_ + 1;
		z_code = current->zLocCode_ + 1;
		level = current->level_;

	}


	void output_point()//计算与输出点的编码
	{

		double q = pow(2, level) + 1;
		point_code[0] = (x_code - 1)*q + y_code + z_code *q*q;
		point_code[1] = (x_code - 1)*q + y_code + 1 + z_code *q*q;
		point_code[2] = x_code *q + y_code + 1 + z_code *q*q;
		point_code[3] = x_code *q + y_code + z_code *q*q;
		
		point_code[4] = (x_code - 1)*q + y_code + (z_code - 1)*q*q;
		point_code[5] = (x_code - 1)*q + y_code + 1 + (z_code - 1)*q*q;
		point_code[6] = x_code *q + y_code + 1 + (z_code - 1)*q*q;
		point_code[7] = x_code *q + y_code + (z_code - 1)*q*q;
		
		fstream outfile;
		outfile.open("d://test.txt", ios::app);
		for (int i = 0; i <= 7; i++)
		{
			outfile << point_code[i] << "  ";
		}
		//outfile <<"    "<< x_code << "  " << y_code << "  " << z_code;
		outfile <<"    	mat 1 crossSect 1	nlgeo 1"<< endl;

		outfile.close();
	}

};





void coor(float max_x, float max_y, float max_z, int level)//三个包装盒的尺寸  分割次数
{
	double lines = pow(2, level);//计算每一行有多少体素   dxdydz分别是三个方向上每个体素的尺寸
	float dx = max_x / lines;
	float dy = max_y / lines;
	float dz = max_z / lines;

	fstream outfile;
	outfile.open("d://test2.txt", ios::out);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//输出每个节点的坐标
	{
		for (int j = 0; j < lines; j++)
		{
			for (int k = 0; k < lines; k++)
			{
				outfile << number++ << "  " << j *dx << "  " << k *dy << "  " << i *dz << endl;

			}
		}
	}
	outfile.close();
}

/*fstream outfile;
outfile.open("d://test.txt", ios::app);
outfile << "LSpace " << count++ << "	 nodes  ";
voxel_output* asd = new voxel_output(pNode);
//voxel_output asd(pChildNode, level);


asd->output_point();
delete asd;
*/