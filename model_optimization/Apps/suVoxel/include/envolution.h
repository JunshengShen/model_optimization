//�����ﶨ���ݻ�����
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
int merton(int x, int y, int z, int level);//����Ī����

int *code(int merton, int level);


int *six_n_merton(int merton_code, int level);

struct ca_envolution
{
	//int x, y, z;  Ī�������x,y,z��Ϣ  ����Ҫ
	int location = 0;
	int merton;
	int level;
	bool lable=0;//�涨1Ϊ�߽�    0Ϊ�ڲ�
	int location_temp = 0;
	bool out = 0;//1Ϊ��� 0Ϊ�����  Ҳ���Ǵ������
	bool out_temp = 1;
	//ca_envolution(int x_, int y_, int z_, int level_,int merton_,int lable_) :x(x_), y(y_), z(z_), level(level_), merton(merton_),lable(lable_){}
};


int asc2(char a);

struct point_displacement_vector
{
	double x_displacement_vector;
	double y_displacement_vector;
	double z_displacement_vector;
};

double trans(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8, char a9, char a10, char a11, char a12, char a13);
void read_point_information(string address);