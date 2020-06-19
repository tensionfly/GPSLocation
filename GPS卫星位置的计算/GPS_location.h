#ifndef GPS_LOCATION_H_
#define GPS_LOCATION_H_
#include"string"
#include"iostream"
#include"fstream"
#include"iomanip"
using namespace std;
double **p_data;   /*指向导航电文数据存储数组的指针*/
int piece = 0;     /*所给导航电文的一个卫星某个时间点的数据块*/

/*=============================================
*接下来是一个数据读取及处理的类的定义，其中包括
*类成员的定义。
*
*strin：用于存储导航电文的路径
*strine:用于存储处理后的导航电文数据的输出路径
*
*datareading函数为类的构造函数
*datapretreatment函数用于实现导航电文数据的预处理
*rank函数用于导航电文数据的内部排序
*datacopy函数用于数据的处理及将数据导入一个数组存储
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
*接下来是一个卫星定位计算的类的定义，其中包括类成员
*的定义。
*
*strout：用于计算卫星位置的输出路径的存储
*Crs,Cuc,Cus,Cic,Cis,Crc,delta_n,omegadot, idot
*：为九个摄动参数
*M0，e,sqrtA,omega0,i0,w：六个卫星轨道参数
*t0：参考时刻
*
*gpslocation为类构造函数
*deal函数用于进行导航电文提供数据的单次计算
*batching函数用于数据的批量计算
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