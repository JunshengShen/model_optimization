#include "envolution.h"
vector<point_displacement_vector> all_pont_displacement_vector;
int merton(int x, int y, int z, int level)
{
	int merton = 0;
	for (int i = level - 1; i >= 0; i--)
	{
		merton += (x >> i) << (3 * i + 2);
		x = x - ((x >> i) << i);
		merton += (y >> i) << (3 * i + 1);
		y = y - ((y >> i) << i);
		merton += (z >> i) << (3 * i);
		z = z - ((z >> i) << i);
	}
	return merton;
}
int *code(int merton, int level)//输入一个莫顿序 输出其对应X Y Z编号
{
	int *code = new int[3]{0, 0, 0};
	for (int i = level - 1; i >= 0; i--)
	{

		code[0] += (merton >> (i * 3 + 2)) << i;
		merton -= (merton >> (i * 3 + 2)) << (i * 3 + 2);
		code[1] += (merton >> (i * 3 + 1)) << i;
		merton -= (merton >> (i * 3 + 1)) << (i * 3 + 1);
		code[2] += (merton >> (i * 3)) << i;
		merton -= (merton >> (i * 3)) << (i * 3);
	}
	return code;
}

int *six_n_merton(int merton_code, int level)//计算六领域的莫顿序
{
	int  *six = new int[26];  //0,1，2,3,4,5分别为X+ X- Y+ Y- Z+ Z-
	int *code_ = new int[3];
	code_ = code(merton_code, level);

	six[0] = merton(code_[0] + 1, code_[1], code_[2], level);
	six[1] = merton(code_[0] - 1, code_[1], code_[2], level);
	six[2] = merton(code_[0], code_[1] + 1, code_[2], level);
	six[3] = merton(code_[0], code_[1] - 1, code_[2], level);
	six[4] = merton(code_[0], code_[1], code_[2] + 1, level);
	six[5] = merton(code_[0], code_[1], code_[2] - 1, level);

	six[6] = merton(code_[0] + 1, code_[1] + 1, code_[2] + 1, level);
	six[7] = merton(code_[0] - 1, code_[1] - 1, code_[2] + 1, level);

	six[8] = merton(code_[0] + 1, code_[1] - 1, code_[2] + 1, level);
	six[9] = merton(code_[0] - 1, code_[1] + 1, code_[2] + 1, level);



	six[10] = merton(code_[0] + 1, code_[1], code_[2] + 1, level);
	six[11] = merton(code_[0] - 1, code_[1], code_[2] + 1, level);
	six[12] = merton(code_[0], code_[1] + 1, code_[2] + 1, level);
	six[13] = merton(code_[0], code_[1] - 1, code_[2] + 1, level);

	six[14] = merton(code_[0] + 1, code_[1] + 1, code_[2], level);
	six[15] = merton(code_[0] + 1, code_[1] - 1, code_[2], level);
	six[16] = merton(code_[0] - 1, code_[1] + 1, code_[2], level);
	six[17] = merton(code_[0] - 1, code_[1] - 1, code_[2], level);

	six[18] = merton(code_[0] + 1, code_[1] + 1, code_[2] - 1, level);
	six[19] = merton(code_[0] + 1, code_[1] - 1, code_[2] - 1, level);
	six[20] = merton(code_[0], code_[1] + 1, code_[2] - 1, level);
	six[21] = merton(code_[0], code_[1] - 1, code_[2] - 1, level);
	six[22] = merton(code_[0] - 1, code_[1] + 1, code_[2] - 1, level);
	six[23] = merton(code_[0] - 1, code_[1] - 1, code_[2] - 1, level);
	six[24] = merton(code_[0], code_[1] + 1, code_[2] - 1, level);
	six[25] = merton(code_[0], code_[1] - 1, code_[2] - 1, level);
	return six;
	delete six;
}

int asc2(char a)
{
	switch (a)
	{
	case '0':
		return 0;
		break;
	case '1':
		return 1;
		break;
	case '2':
		return 2;
		break;
	case '3':
		return 3;
		break;
	case '4':
		return 4;
		break;
	case '5':
		return 5;
		break;
	case '6':
		return 6;
		break;
	case '7':
		return 7;
		break;
	case '8':
		return 8;
		break;
	case '9':
		return 9;
		break;
	}
}
double trans(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8, char a9, char a10, char a11, char a12, char a13)
{
	double a = 0;
	a += asc2(a1);
	a += 0.1*asc2(a3);
	a += 0.01*asc2(a4);
	a += 0.001*asc2(a5);
	a += 0.0001*asc2(a6);
	a += 0.00001*asc2(a7);
	a += 0.000001*asc2(a8);
	int b = 0;
	b += asc2(a11) * 100 + asc2(a12) * 10 + asc2(a13);
	if (a10 == '-')
	{
		b = b*(-1);
	}
	a = a*pow(10, b);
	return a;
}
void read_point_information(string address)
{
	fstream in;
	char read_temp;
	in.open(address, ios::in);
	while (!in.eof())
	{
		in >> read_temp;
		if (read_temp == 'N')
		{
			in >> read_temp;
			if (read_temp == 'a')
			{
				in >> read_temp;
				if (read_temp == 'm')
				{
					in >> read_temp;
					if (read_temp == 'e')
					{
						in >> read_temp;
						if (read_temp == '=')
						{
							in >> read_temp;
							if (read_temp == '"')
							{
								in >> read_temp;
								if (read_temp == 'D')
								{
									in >> read_temp;
									if (read_temp == 'i')
									{
										for (;;)
										{
											in >> read_temp;
											//if (read_temp == '<')
											//break;
											int break_ = 0;
											if (read_temp == '>')
											{
												point_displacement_vector inf_temp;
												double temp_vector[3];
												char a[13];
												int count = 0;

												for (;;)
												{
													//in >> read_temp;
													in >> read_temp;
													if (read_temp == '<')
													{
														break_++;
														break;
													}
													if (read_temp == '-')
													{
														for (int i = 0; i < 13; i++)
														{
															in >> a[i];
														}
														temp_vector[count] = -1 * trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12]);
													}
													else
													{
														a[0] = read_temp;
														for (int i = 1; i < 13; i++)
														{
															in >> a[i];
														}
														temp_vector[count] = trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12]);
													}
													count++;
													if (count == 3)
													{
														count = 0;
														inf_temp.x_displacement_vector = temp_vector[0];
														inf_temp.y_displacement_vector = temp_vector[1];
														inf_temp.z_displacement_vector = temp_vector[2];
														all_pont_displacement_vector.push_back(inf_temp);
														cout << inf_temp.x_displacement_vector << " " << inf_temp.y_displacement_vector << " " << inf_temp.z_displacement_vector << endl;
													}
												}
												break;
											}
											if (break_)
												break;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	in.close();
}
