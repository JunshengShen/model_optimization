
//类voxel_output中能够生成一个体素的八个顶点的编码   构造函数中输入目标体素的指针       output_point()函数输出八个点的坐标在TXT文件中  
//函数coor输入x,y,z的最大值与分割次数  输出所有点的坐标
//coor函数覆盖原来的TXT文件
//voxel_output函数在使用时需要先out一个TXT文件
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
	voxel_output(int x_, int y_, int z_, int levell) 
	{
		x_code = x_ + 1;
		y_code = y_ + 1;
		z_code = z_ + 1;
		level = levell;
	};


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
	void output_point1()//计算与输出点的编码
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
		outfile.open("d://test1.txt", ios::app);
		for (int i = 0; i <= 7; i++)
		{
			outfile << point_code[i] << "  ";
		}
		//outfile <<"    "<< x_code << "  " << y_code << "  " << z_code;
		outfile << "    	mat 1 crossSect 1	nlgeo 1" << endl;

		outfile.close();
	}



	void output_point2(int number_)//计算与输出点的编码
	{
		double q = pow(2, level) + 1;
		point_code[0] = (x_code - 1)*q + y_code + z_code *q*q;
		point_code[1] = (x_code - 1)*q + y_code + 1 + z_code *q*q;
		point_code[3] = x_code *q + y_code + 1 + z_code *q*q;
		point_code[2] = x_code *q + y_code + z_code *q*q;

		point_code[4] = (x_code - 1)*q + y_code + (z_code - 1)*q*q;
		point_code[5] = (x_code - 1)*q + y_code + 1 + (z_code - 1)*q*q;
		point_code[7] = x_code *q + y_code + 1 + (z_code - 1)*q*q;
		point_code[6] = x_code *q + y_code + (z_code - 1)*q*q;

		fstream outfile;
		outfile.open("d://test2.txt", ios::app);
		/*for (int i = 0; i <= 7; i++)
		{
			outfile << point_code[i] << "  ";
		}
		//outfile <<"    "<< x_code << "  " << y_code << "  " << z_code;
		outfile << "    	mat 1 crossSect 1	nlgeo 1" << endl;*/


		outfile << "LTRSpace  " << 5 * number_ + 1 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[3] << "  "
			<< point_code[2] << "  " << point_code[6] << "  mat 1 crossSect 1	NIP  1" << endl;
		outfile << "LTRSpace  " << 5 * number_ + 2 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[3] << "  "
			<< point_code[6] << "  " << point_code[5] << "  mat 1 crossSect 1	NIP  1" << endl;
		outfile << "LTRSpace  " << 5 * number_ + 3 << "  " << "nodes  4	  " << point_code[3] << "  " << point_code[6] << "  "
			<< point_code[5] << "  " << point_code[7] << "  mat 1 crossSect 1	NIP  1" << endl;
		outfile << "LTRSpace  " << 5 * number_ + 4 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[1] << "  "
			<< point_code[3] << "  " << point_code[5] << "  mat 1 crossSect 1	NIP  1" << endl;
		outfile << "LTRSpace  " << 5 * number_ + 5 << "  " << "nodes  4	  " << point_code[4] << "  " << point_code[6] << "  "
			<< point_code[5] << "  " << point_code[0] << "  mat 1 crossSect 1	NIP  1" << endl;









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
	outfile.open("d://test.txt", ios::app);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//输出每个节点的坐标
	{
		for (int j = 0; j < lines; j++)
		{
			for (int k = 0; k < lines; k++)
			{
				
				if (k == 0)
				{
					outfile <<"node   "<< number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << "       bc 3 2 2 2" << endl;
				}
				else
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << endl;
			}
		}
	}
	outfile.close();
}











void coor1(float max_x, float max_y, float max_z, int level)//三个包装盒的尺寸  分割次数
{
	double lines = pow(2, level);//计算每一行有多少体素   dxdydz分别是三个方向上每个体素的尺寸
	float dx = max_x / lines;
	float dy = max_y / lines;
	float dz = max_z / lines;

	fstream outfile;
	outfile.open("d://test1.txt", ios::app);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//输出每个节点的坐标
	{
		for (int j = 0; j < lines; j++)
		{
			for (int k = 0; k < lines; k++)
			{

				if (k == 0)
				{
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << "       bc 3 2 2 2" << endl;
				}
				else
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << endl;
			}
		}
	}
	outfile.close();
}






void coor2(float max_x, float max_y, float max_z, int level)//三个包装盒的尺寸  分割次数
{
	double lines = pow(2, level);//计算每一行有多少体素   dxdydz分别是三个方向上每个体素的尺寸
	float dx = max_x / lines;
	float dy = max_y / lines;
	float dz = max_z / lines;

	fstream outfile;
	outfile.open("d://test2.txt", ios::app);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//输出每个节点的坐标
	{
		for (int j = 0; j < lines; j++)
		{
			for (int k = 0; k < lines; k++)
			{

				if (k == 0)
				{
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << "       bc 3 1 1 1" << endl;
				}
				else
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << endl;
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

/*
int lines = 2;
if (divided_times >= 1)
for (int i = 1; i <= divided_times; i++)
{
lines += pow(2, i - 1);
}
lines = lines*lines*lines;
for (long double i = 0; i < lines; i++)
{
int line_number = 2;//行，列，层的数量
if (divided_times >= 1)
for (int i = 1; i <= divided_times; i++)
{
line_number += pow(2, (i - 1));
}
int line_number_z;
line_number_z = i / (line_number*line_number) + 1;
int top_x_y = i - (line_number_z - 1)*(line_number*line_number);
int line_number_x;
line_number_x = top_x_y / line_number + 1;
int line_number_y;
if (line_number_x>1)
line_number_y = top_x_y - (line_number_x - 1)*line_number + 1;
else line_number_y = top_x_y + 1;
float dx = mx / (line_number - 1);
float dy = my / (line_number - 1);
float dz = mz / (line_number - 1);
float  coo_x = dx*(line_number_x - 1);
float  coo_y = dy*(line_number_y - 1);
float  coo_z = dy*(line_number_z - 1);
cout << i + 1 << "     " << coo_x << "  " << coo_y << "  " << coo_z << endl;
outfile << i + 1 << "     " << coo_x << "  " << coo_y << "  " << coo_z << endl;
}

*/