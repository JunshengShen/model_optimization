//�����ﶨ���ݻ�����

int merton(int x, int y, int z, int level)//����Ī����
{
	int merton = 0;
	for (int i = level - 1; i >= 0; i--)
	{
		merton += (z >> i) << (3 * i + 2);
		z = z - ((z >> i) << i);
		merton += (y >> i) << (3 * i + 1);
		y = y - ((y >> i) << i);
		merton += (x >> i) << (3 * i);
		x = x - ((x >> i) << i);
	}
	return merton;
}
int *code(int merton, int level)//����һ��Ī���� ������ӦX Y Z���
{
	int *code = new int[3]{0, 0, 0};
	for (int i = level - 1; i >= 0; i--)
	{

		code[2] += (merton >> (i * 3 + 2)) << i;
		merton -= (merton >> (i * 3 + 2)) << (i * 3 + 2);
		code[1] += (merton >> (i * 3 + 1)) << i;
		merton -= (merton >> (i * 3 + 1)) << (i * 3 + 1);
		code[0] += (merton >> (i * 3)) << i;
		merton -= (merton >> (i * 3)) << (i * 3);
	}
	return code;
}
int *six_n_merton(int merton_code,int level)//�����������Ī����
{
	int  *six = new int[26];  //0,1��2,3,4,5�ֱ�ΪX+ X- Y+ Y- Z+ Z-
	int *code_ = new int[3];
	code_ = code(merton_code, level);
	
	six[0] = merton(code_[0] + 1, code_[1], code_[2], level);
	six[1] = merton(code_[0] - 1, code_[1], code_[2], level);
	six[2] = merton(code_[0], code_[1] + 1, code_[2], level);
	six[3] = merton(code_[0], code_[1] - 1, code_[2], level);
	six[4] = merton(code_[0], code_[1], code_[2] + 1, level);
	six[5] = merton(code_[0], code_[1], code_[2] - 1, level);

	six[6] = merton(code_[0] + 1, code_[1]+1, code_[2]-1, level);
	six[7] = merton(code_[0]-1 , code_[1]-1, code_[2]-1, level);
	six[8] = merton(code_[0]-1, code_[1]+1, code_[2]-1, level);
	six[9] = merton(code_[0] +1, code_[1]-1, code_[2]-1, level);
	six[10] = merton(code_[0] , code_[1]+1, code_[2]-1, level);
	six[11] = merton(code_[0] , code_[1]-1, code_[2]-1, level);
	six[12] = merton(code_[0]+1 , code_[1], code_[2]-1, level);
	six[13] = merton(code_[0]-1 , code_[1], code_[2]-1, level);
	six[14] = merton(code_[0]+1 , code_[1]+1, code_[2]+1, level);
	six[15] = merton(code_[0]-1 , code_[1]-1, code_[2]+1, level);
	six[16] = merton(code_[0]-1 , code_[1]+1, code_[2]+1, level);
	six[17] = merton(code_[0]+1 , code_[1]-1, code_[2]+1, level);
	six[18] = merton(code_[0] , code_[1]+1, code_[2]+1, level);
	six[19] = merton(code_[0] , code_[1]-1, code_[2]+1, level);
	six[20] = merton(code_[0] +1, code_[1], code_[2]+1, level);
	six[21] = merton(code_[0] -1, code_[1], code_[2]+1, level);
	six[22] = merton(code_[0]+1 , code_[1]+1, code_[2], level);
	six[23] = merton(code_[0]-1 , code_[1]-1, code_[2], level);
	six[24] = merton(code_[0]+1 , code_[1]-1, code_[2], level);
	six[25] = merton(code_[0]-1 , code_[1]+1, code_[2], level);

	

	return six;
	delete six;
}

struct ca_envolution
{
	//int x, y, z;  Ī�������x,y,z��Ϣ  ����Ҫ
	int location = 0;
	int merton;
	int level;
	bool lable=0;//�涨1Ϊ�߽�    0Ϊ�ڲ�
	int location_temp = 0;
	bool out = 1;//1Ϊ��� 0Ϊ�����  Ҳ���Ǵ������
	bool out_temp = 1;
	//ca_envolution(int x_, int y_, int z_, int level_,int merton_,int lable_) :x(x_), y(y_), z(z_), level(level_), merton(merton_),lable(lable_){}
};