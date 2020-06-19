#include"GPS_location.h"
#define GM 3.986005e14
#define We 7.292115e-5
#define row 36
#define PI 3.141592653589793
/*===========================
*���϶����˳������õ��ĳ���
*GM:�������������͵��������ĳ˻�
*We:������ת��ƽ�����ٶ�
*row���洢�������������
*PI��Բ����
*============================*/

/*���������๹�캯���ľ���ʵ�֣����벢�洢���ݵ���������ļ���·��*/
datareading::datareading(string str_,string stre_)  
{
	strin = str_;
	strine = stre_;
}

/*==============================================
*�������ĺ���rank����ð�����򷨽�ʵ��һ�������ڿ�״
*������������š�
*
*row���洢�����������
*startline����Ҫʵ�����������������ʼ�С�
*endline����Ҫʵ�����������������ֹ�С�
*p_data��ָ��洢���ݵ������ָ�롣
================================================*/
void datareading::rank(int startline, int endline)
{
	double k;
	for (int i = startline; i < endline; i++)
	{
		for (int j = i + 1; j < endline+1; j++)
		{
			double k1, k2;
			k1 = p_data[j][3];
			k2 = p_data[j][4] * 3600 + p_data[j][5] * 60 + p_data[j][6];

			double k1_, k2_;
			k1_ = p_data[i][3];
			k2_ = p_data[i][4] * 3600 + p_data[i][5] * 60 + p_data[i][6];

			if (k1 < k1_ || (k1 == k1_&&k2<k2_))
			{
				for (int l = 0; l < row; l++)
				{
					k = p_data[i][l];
					p_data[i][l] = p_data[j][l];
					p_data[j][l] = k;
				}
			}
		}
	}
}

/*==================================================
*��������datapretreatment������ʵ���������ݵĲ��ִ���
*1.��������������ļ����������
*2.����ǰ���е�ͷ�ļ�������������������
*  �Ŀ�ѧ������ת��Ϊc++��
*  �Ŀ�ѧ�����ı�﷽����
*3.ͳ�������������ĵ�һ������ĳ��ʱ�������ݿ����Ŀ��
===================================================*/
void  datareading::datapretreatment()
{
	int j=0;
	string line;
	ifstream infile(strin);
	ofstream outfile(strine);
	infile.seekg(307, ios_base::beg);
	while (!infile.eof())
	{
	
			getline(infile, line);
			for (int i = 0; i < line.length(); i++)
			{
				if (line[i] == 'D')
				{
					line[i] = 'e';
					j++;
				}
			}
			outfile << line << endl;
		
	}
	outfile.close();
	infile.close();
	piece = j / 29;
}

/*=======================================================
*��������datacopy������ʵ�ֻ����rank�����������������ݵ�����
*
*1��2������һ����̬���顣
*2��3�������������ݴ��ļ����뵽�����Ķ�̬�����С�
*3��4����PRN��ͬ��һ������������һ��
*4��5������rank����ʵ��PRN��ͬ�����ݰ�ʱ���Ⱥ�˳���������
=========================================================*/
void datareading::datacopy()
{
	p_data = new double *[piece];                       //1
	for (int i = 0; i < piece; i++)
	{
		p_data[i] = new double[36];
	}                                                   //2

	ifstream infile(strine);

	for (int j = 0; j < piece;j++)
	{
		for (int k = 0; k < 36; k++)
		{
			infile >> p_data[j][k];
		}
	}                                                   //3

	double  k;
	for (int i = 0; i < piece-1;)
	{
		for (int j = i + 1; j < piece; j++)
		{
			if (p_data[i][0] == p_data[j][0] && abs(i - j) == 1)
				i++;
			if (p_data[i][0] == p_data[j][0] && abs(i - j) > 1)
			{
				for (int l = 0; l < 36; l++)
				{
					k = p_data[i + 1][l];
					p_data[i + 1][l] = p_data[j][l];
					p_data[j][l] = k;
				}
				i++;
			}
		}
		i++;
	}                                                    //4

	int position_fro, positon_beh;
	for (int i = 0; i < piece-1; i++)
	{
		if (i == 0)
		 position_fro = 0;
		if (p_data[i][0] != p_data[i + 1][0])
		{
			positon_beh = i;
			rank(position_fro, positon_beh);
			position_fro = positon_beh+1;
		}
		if (i == piece - 2)
		{
			positon_beh = piece - 1;
			rank(position_fro, positon_beh);
		}

	}                                                    //5

/*	ofstream fout("C:\\Users\\tension fly\\Desktop\\920.txt");
	for (int j = 0; j < piece; j++)
	{
		for (int k = 0; k < 36; k++)
		{
			fout << p_data[j][k] << "  ";
		}
		fout << endl << endl;
	}
*/
}

/*���������๹�캯���ľ���ʵ�֣����벢�洢������������ļ���·�����洢���������ָ��*/
gpslocation::gpslocation(string str_,double **p_)
{
	strout = str_;
	p2data = p_;
}

/*===================================
*��������deal������ʵ������λ�õļ��㣺
*1��2�������Ķ��弰������Ĵ���
*2��3���õ�������⿪���շ��̡�
*3��4�������ǵ������жϼ�ƽ����ǵļ��㡣
*4��5���㶯������
*5��6��������ƽ������x,y��
*6��7��˲ʱ������ľ��ȡ�
*7��8��������Ŀռ�����ϵ��X,Y,Z��
*8��9���������������ָ���ļ���
=====================================*/
void gpslocation::deal(int PRN,double t)
{
	double n0, n, M, E, gd,cosf,sinf, f, u_;      //1
	n0 = sqrt(GM) / pow(sqrtA, 3);
	n = delta_n + n0;
	M = M0 + n*(t - t0);
	                                              //2
	E = M;
	for (;;)
	{
		gd = E;
		E = M + e*sin(E);
		if (fabs(gd - E) < 1e-5)
			break;
	}
	                                             //3
	cosf = (cos(E) - e) / (1 - e*cos(E));
	sinf = sqrt(1 - e*e)*sin(E) / (1 - e*cos(E));
	f = atan(sqrt(1 - e*e)*sin(E) / (cos(E) - e));

	if (cosf < 0)
		f += PI;

	u_ = w + f;
	                                             //4
	double delta_u, delta_r, delta_i, u, r, i;

	delta_u = Cuc*cos(2 * u_) + Cus*sin(2 * u_);
	delta_r = Crc*cos(2 * u_) + Crs*sin(2 * u_);
	delta_i = Cic*cos(2 * u_) + Cis*sin(2 * u_);

	u = u_ + delta_u;
	r = sqrtA*sqrtA*(1 - e*cos(E)) + delta_r;
	i = i0 + delta_i + idot*(t - t0);
	                                             //5
	double x, y;

	x = r*cos(u);
	y = r*sin(u);
	                                             //6
	double L;

	L = omega0 + (omegadot - We)*t - omegadot*t0;
	                                             //7
	double X, Y, Z;

	X = x*cos(L) - y*cos(i)*sin(L);
	Y = x*sin(L) + y*cos(i)*cos(L);
	Z = y*sin(i);
	                                            //8
	ofstream outfile(strout, ios::app);
	outfile<<PRN<<" "<<t<<" "<< X << "  " << Y << "  " << Z << endl;
	outfile.close();
	                                            //9
}

/*===========================================
*��������batching����������deal����ʵ��ÿ������
*ÿ���30���λ�ü��㼰�����
==============================================*/
void gpslocation::batching()
{
	int j,prn;
	double t;
	for (int i = 0; i < piece-1;i++)
	{
		if (p2data[i][0]==p2data[i + 1][0])
		{
			j =int (((p2data[i + 1][3] - p2data[i][3]) * 86400 + (p2data[i + 1][4] - p2data[i][4]) * 3600 + (p2data[i + 1][5] - p2data[i][5]) * 60 + (p2data[i + 1][6] - p2data[i][6])) / 30 +0.5);
			p2data[i + 1][18] = p2data[i][18] + (p2data[i + 1][3] - p2data[i][3]) * 86400 + (p2data[i + 1][4] - p2data[i][4]) * 3600 + (p2data[i + 1][5] - p2data[i][5]) * 60 + (p2data[i + 1][6] - p2data[i][6]);
			
			Crs = p2data[i][11];
			delta_n = p2data[i][12];
			M0 = p2data[i][13];
			Cuc = p2data[i][14];
			e = p2data[i][15];
			Cus = p2data[i][16];
			sqrtA = p2data[i][17];
			t0 = p2data[i][18];
			Cic = p2data[i][19];
			omega0 = p2data[i][20];
			Cis = p2data[i][21];
			i0 = p2data[i][22];
			Crc = p2data[i][23];
			w = p2data[i][24];
			omegadot = p2data[i][25];
			idot = p2data[i][26];

			prn = (int)p2data[i][0];
			for (int k = 0; k < j; k++)
			{
				t = t0 + k * 30;
				deal(prn,t);
			}
		}
	
	}
	cout << "�������" << endl;

}

/*��������������ʵ�����ʵ��������غ����ĵ��ã�ʵ����Ӧ�Ĵ����ܡ�*/
int main()
{
	datareading read("C:\\Users\\tension fly\\Desktop\\����.txt","C:\\Users\\tension fly\\Desktop\\����e.txt");
	read.datapretreatment();
	read.datacopy();

	gpslocation location("C:\\Users\\tension fly\\Desktop\\resualt_xyz.txt", p_data);
	location.batching();

	system("pause");
	return 0;
}