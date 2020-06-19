#ifndef GPS_LOCATION_H_
#define GPS_LOCATION_H_
#include"string"
#include"iostream"
#include"fstream"
#include"iomanip"
using namespace std;
double **p_data;   /*ָ�򵼺��������ݴ洢�����ָ��*/
int piece = 0;     /*�����������ĵ�һ������ĳ��ʱ�������ݿ�*/

/*=============================================
*��������һ�����ݶ�ȡ���������Ķ��壬���а���
*���Ա�Ķ��塣
*
*strin�����ڴ洢�������ĵ�·��
*strine:���ڴ洢�����ĵ����������ݵ����·��
*
*datareading����Ϊ��Ĺ��캯��
*datapretreatment��������ʵ�ֵ����������ݵ�Ԥ����
*rank�������ڵ����������ݵ��ڲ�����
*datacopy�����������ݵĴ��������ݵ���һ������洢
===============================================*/
class datareading
{
private:
	string strin,strine;
	
public:
	datareading(string str_,string stre_);
	void datapretreatment();
	void rank(int startline, int endline);
	void datacopy();
};

/*=============================================
*��������һ�����Ƕ�λ�������Ķ��壬���а������Ա
*�Ķ��塣
*
*strout�����ڼ�������λ�õ����·���Ĵ洢
*Crs,Cuc,Cus,Cic,Cis,Crc,delta_n,omegadot, idot
*��Ϊ�Ÿ��㶯����
*M0��e,sqrtA,omega0,i0,w���������ǹ������
*t0���ο�ʱ��
*
*gpslocationΪ�๹�캯��
*deal�������ڽ��е��������ṩ���ݵĵ��μ���
*batching�����������ݵ���������
*==============================================*/
class gpslocation
{
private:
	string strout;
	double Crs, delta_n, M0, Cuc, e, Cus, sqrtA, t0, Cic, omega0, Cis, i0, Crc, w, omegadot, idot;
	double **p2data;

public:
	gpslocation(string str_,double **p_);
	void deal(int PRN, double t);
	void batching();
};

#endif