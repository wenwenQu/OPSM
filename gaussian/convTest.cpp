/*//这里调用的"convol.h"好像和algo.h冲突
#include <iostream>
#include <ctime>
#include "convol.h"
using namespace std;

//===================== 以下为DP_PMF、PMFCheck测试函数 =====================

extern int size_th;//这里要调试出一个最好的size_th

int main1()
{//调试size_th的值
	int n=10000;
	int sn=10;
	srand((unsigned)time(0));
	clock_t start, end;

	double** pr=new double*[sn];
	for(int j=0; j<sn; j++)
	{
		pr[j]=new double[n];
		for(int i=0; i<n; i++)
		{
			pr[j][i]=(double)rand()/RAND_MAX;
		}
	}

	cout<<"Naive DP_PMF"<<endl;
	start=clock();
	for(int j=0; j<sn; j++) double* r1=DP_PMF(pr[j], n);
	end=clock();
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	cout<<"DC_PMF0 (NoPrune) by FFT"<<endl;
	start=clock();
	for(int j=0; j<sn; j++) double* r2=DC_PMF0(pr[j], n);
	end=clock();
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	//for(int th=8192; th>=16; th/=2)//粗粒度实验: 结果size_th = 512 -> 256时最好
	for(int th=1000; th>=100; th-=50)//细粒度实验: 结果size_th = 600 -> 400时最好
	{
		size_th=th;
		cout<<"------ size_th = "<<size_th<<" ------"<<endl;

		cout<<"DC_PMF (NoPrune) by FFT"<<endl;
		start=clock();
		for(int j=0; j<sn; j++) double* r2=DC_PMF(pr[j], n);
		end=clock();
		cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
	}

	return 0;
}
*/

/*
int main()
{
	int n=2000;
	srand((unsigned)time(0));
	clock_t start, end;

	double* pr=new double[n];
	for(int i=0; i<n; i++)
	{
		pr[i]=(double)rand()/RAND_MAX;
	}

	cout<<"Naive DP_PMF"<<endl;
	start=clock();
	double* r1=DP_PMF(pr, n);
	end=clock();
	for(int i=0; i<n+1; i++) cout<<r1[i]<<"\t";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	cout<<"DC_PMF (NoPrune) by FFT"<<endl;
	start=clock();
	double* r2=DC_PMF(pr, n);
	end=clock();
	for(int i=0; i<n+1; i++) cout<<r2[i]<<"\t";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	cout<<"PMFCheck (Prune) by FFT"<<endl;
	start=clock();
	double* r3;
	PMFCheck(r3, pr, n, n/10, 2);
	end=clock();
	for(int i=0; i<n+1; i++) cout<<r3[i]<<"\t";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	cout<<"PMFCheck0 (Prune) by FFT"<<endl;
	start=clock();
	double* r4;
	PMFCheck0(r4, pr, n, n/10, 2);
	end=clock();
	for(int i=0; i<n+1; i++) cout<<r3[i]<<"\t";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	return 0;
}
//*/

//===================== 以下为convolution测试函数 =====================
/*
int main()
{
	int n1=1000, n2=1000;
	srand((unsigned)time(0));

	double* v1=new double[n1];
	for(int i=0; i<n1; i++)
	{
		v1[i]=(double)rand()/RAND_MAX;
	}

	double* v2=new double[n2];
	for(int i=0; i<n2; i++)
	{
		v2[i]=(double)rand()/RAND_MAX;
	}

	int n=n1+n2-1;
	clock_t start, end;

	cout<<"Convolution by FFT"<<endl;
	start=clock();
	double* r1=convFFT(v1, n1, v2, n2);
	end=clock();
	for(int i=0; i<n; i++) cout<<r1[i]<<" ";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	cout<<"Naive Convolution"<<endl;
	start=clock();
	double* r2=convNaive(v1, n1, v2, n2);
	end=clock();
	for(int i=0; i<n; i++) cout<<r2[i]<<" ";
	cout<<endl;
	cout<<"Running time: "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

	return 0;
}
*/

/*
int main()
{
	double v1[]={0.23, 0.42, 0.678, 0.45, 0.96, 0.678, 0.56};
	double v2[]={0.44, 0.67, 0.12, 0.67, 0.77, 0.82};
	int n1=7;
	int n2=6;

	int n=n1+n2-1;
	double* r=convFFT(v1, n1, v2, n2);
	for(int i=0; i<n; i++) cout<<r[i]<<" ";

	return 0;
}
//*/

/*
int main()
{
	double v1[]={0.23, 0.42, 0.678, 0.45, 0.96, 0.678, 0.56};
	double v2[]={0.44, 0.67, 0.12, 0.67, 0.77, 0.82};
	int n1=7;
	int n2=6;

	int n=n1+n2-1;
	double* r=convNaive(v1, n1, v2, n2);
	for(int i=0; i<n; i++) cout<<r[i]<<" ";

	return 0;
}
//*/

/*//产生0-1之间的随机数
int main()
{
	srand((unsigned)time(0));
	while(1)
	{
		cout <<(double)rand()/RAND_MAX;
		cout <<endl;
	}
	return 0;
}
*/

//===============================================================
