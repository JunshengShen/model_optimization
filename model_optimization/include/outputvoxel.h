//��voxel_output���ܹ�����һ�����صİ˸�����ı���   ���캯��������Ŀ�����ص�ָ��       output_point()��������˸����������TXT�ļ���  
//����coor����x,y,z�����ֵ��ָ����  ������е������
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
	voxel_output(suObejctOctree<NodeData>::OctreeNode* current) //��Ŀ�������еĲ�����Ӧ��ֵ
	{
		x_code = current->xLocCode_ + 1;
		y_code = current->yLocCode_ + 1;
		z_code = current->zLocCode_ + 1;
		level = current->level_;

	}


	void output_point()//�����������ı���
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





void coor(float max_x, float max_y, float max_z, int level)//������װ�еĳߴ�  �ָ����
{
	double lines = pow(2, level);//����ÿһ���ж�������   dxdydz�ֱ�������������ÿ�����صĳߴ�
	float dx = max_x / lines;
	float dy = max_y / lines;
	float dz = max_z / lines;

	fstream outfile;
	outfile.open("d://test2.txt", ios::out);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//���ÿ���ڵ������
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