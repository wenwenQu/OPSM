#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stack>
#include <algorithm>
#include <assert.h>

#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <boost/math/distributions/poisson.hpp>
using namespace boost::math;
using namespace std;

#include "datStruct.h"
#include "io.h"
#include "convol.h"

#define EXP 2.71828182845904523536
#define PI 3.14159265397


double RationalApproximation(double t)
{
	// Abramowitz and Stegun formula 26.2.23.
	// The absolute value of the error should be less than 4.5 e-4.
	double c[] = { 2.515517, 0.802853, 0.010328 };
	double d[] = { 1.432788, 0.189269, 0.001308 };
	return t - ((c[2] * t + c[1])*t + c[0]) /
		(((d[2] * t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
	if (p <= 0.0 || p >= 1.0)
	{
		stringstream os;
		os << "Invalid input argument (" << p
			<< "); must be larger than 0 but less than 1.";
		throw invalid_argument(os.str());
	}

	// See article above for explanation of this section.
	if (p < 0.5)
	{
		// F^-1(p) = - G^-1(p)
		return -RationalApproximation(sqrt(-2.0*log(p)));
	}
	else
	{
		// F^-1(p) = G^-1(1-p)
		return RationalApproximation(sqrt(-2.0*log(1 - p)));
	}
}


void getGu(Range& rng, double& u, double& sigma){
	u = (rng.l + rng.r)/2;

}

//=============== 以下为DP求概率的算法 ===============

//给定t1<...<ti对应的DPMatrix (in), 求t1<...<ti<t(i+1)对应的DPMatrix (out)
void appendRng(DPMatrix &out, DPMatrix &in, Range& rng)//rng=t(i+1)
{

	out.seqLen = in.seqLen + 1;
	out.seq = new Range[out.seqLen];
	for (int i = 0; i < in.seqLen; i++) out.seq[i] = in.seq[i];
	out.seq[in.seqLen] = rng;
	//设置列信息，列数+splits vec
	vector<double> newSplits;
	newSplits.reserve(in.rngNum + total_std); //wenwen: why??
	int pos1 = 0;//in.splits的下标
	int pos2 = 0;//0表示rng.l, 1表示rng.r
	//merge "in.splits"和"{rng.l, rng.r}" (去重复)  	//update by wenwen

	while (pos1 <= in.rngNum && pos2 < total_std)
	{
		double v1 = in.splits[pos1];
		double v2 = rng.u+(pos2-std_num)*rng.sigma;//((pos2 == 0) ? rng.l : rng.r);
		//cout << v2 << endl;
		if (v1 < v2)
		{
			newSplits.push_back(v1);
			pos1++;
		}
		else//v1>=v2
		{
			newSplits.push_back(v2);
			pos2++;
			if (v1 == v2) pos1++;//去重复
		}
	}
	while (pos1 <= in.rngNum)
	{
		double v1 = in.splits[pos1];
		newSplits.push_back(v1);
		pos1++;
	}
	while (pos2 < total_std)
	{
		double v2 = rng.u+(pos2-std_num)*rng.sigma;//((pos2 == 0) ? rng.l : rng.r);
		newSplits.push_back(v2);
		pos2++;
	}


	//设置列数
	out.rngNum = newSplits.size() - 1;
	//以下newSplits -> out.splits
	out.splits = new double[newSplits.size()];
	vector<double>::iterator it;
	pos1 = 0;//归并完pos1没用了，这里重用该变量
	for (it = newSplits.begin(); it != newSplits.end(); it++)
	{
		out.splits[pos1] = *it;
		pos1++;
	}
	//以下设置DP数组
	out.DP = new double*[out.seqLen];
	for (int i = 0; i < out.seqLen; i++) out.DP[i] = new double[out.rngNum];
	//以下获取概率数组(对于每个ti一个entry)

	//update by wenwen
	out.pr = new Polynomial<double>**[out.rngNum];
	int pos = -1;
	for (int i = 0; i < out.rngNum; i++) {
		out.pr[i] = new Polynomial<double>* [out.seqLen];

		//reuse pr[pos][j] of "in"
		if(out.splits[i] == in.splits[pos + 1])
			pos++;
		for (int j = 0; j < out.seqLen-1; j++){
			if(pos >=0 && pos <in.rngNum)
				out.pr[i][j] = in.pr[pos][j]; // here pr[pos][j] of "in" may be NULL
			else
				out.pr[i][j] = NULL; // out of t(j)'s range
		}

		//for the last t
		int j = out.seqLen-1;
		if(!out.seq[j].contains(newSplits[i], newSplits[i + 1]))
			out.pr[i][j] = NULL; // out of t(j)'s range
		int left = 0, right = 2*std_num;
		int middle=0;

		while (left <= right){
			middle = left+(right-left)/2;
			if(newSplits[i] >= (out.seq[j].u+(middle - std_num )*out.seq[j].sigma))
				left = middle+1;
			else if(newSplits[i] < (out.seq[j].u+(middle - std_num )*out.seq[j].sigma))
				right = middle-1;
		}

		if (right < 0 || right >= 2 * std_num)
			out.pr[i][j] = NULL;
		else out.pr[i][j] = &(out.seq[j].pdf(right));
	}


	//以下设置DP元素值
	for (int i = 0; i < out.seqLen - 1; i++)//-1因为最后一行没有东西可以重用
	{
		int pos = 0;//在in.DP中列的位置, 此时对应于<=in.splits[pos+1]的范围内
		for (int j = 0; j < out.rngNum; j++)
		{
			//Range cur(out.splits[j], out.splits[j + 1]);
			//range结果是<=cur.r的范围内

			//by wenwen: pay as you go???
			if (pos < in.rngNum && out.splits[j + 1] == in.splits[pos + 1])
			{//可以重用
				out.DP[i][j] = in.DP[i][pos];
				pos++;
			}
			else
			{//计算该entry，此时前面需要用到的entry应该已经计算好了
				//double gap = cur.gap();
				if (j == 0)//DP对应于y=0的base case
				{
					/*
					double term = 1;
					int k = i;//从i往前扫
					while (out.seq[k].contains(out.splits[j], out.splits[j + 1]))//在cur区间上pdf不为0
					{

						term = term * pr[k] * gap / (i - k + 1);

						k--;
						if (k < 0) break;  out.DP[i][j] = term;
					}

					if (k < 0)

					else out.DP[i][j] = 0;
					*/

					out.DP[i][j] = multiIntegral(out.pr[j], i+1, out.splits[j], out.splits[j + 1]); //当pr[j]为0，递归函数会直接返回0，因此不需要判断
				}
				else
				{
					//以下计算概率部分
					double sum = out.DP[i][j - 1];
					double term = 1;//递推项：subrange上的概率
					int k = i;//从i往前扫
					while (out.seq[k].contains(out.splits[j], out.splits[j + 1]))
					{
						term = multiIntegral(out.pr[j]+k, i-k+1, out.splits[j], out.splits[j + 1]);
						sum += k > 0 ?  term * out.DP[k - 1][j - 1] : term;

						k--;
						if (k < 0) break;
					}
					out.DP[i][j] = sum;
				}
			}
		}
	}
	//计算最后一行
	for (int j = 0; j < out.rngNum; j++)
	{
		int i = out.seqLen - 1;
		//Range cur(out.splits[j], out.splits[j + 1]);
		//double gap = cur.gap();
		if (j == 0)//DP对应于y=0的base case
		{
			/*
			double term = 1;
			int k = i;//从i往前扫
			while (out.seq[k].contains(out.splits[j], out.splits[j + 1]) && pr[k] > 0)//在cur区间上pdf不为0
			{
				term = term * pr[k] * gap / (i - k + 1);
				k--;
				if (k < 0) break;
			}
			if (k < 0) out.DP[i][j] = term;
			else out.DP[i][j] = 0;
			*/
			out.DP[i][j] = multiIntegral(out.pr[j], i+1, out.splits[j], out.splits[j + 1]);
		}
		else
		{
			//以下计算概率部分
			double sum = out.DP[i][j - 1];
			double term = 1;//递推项：subrange上的概率
			int k = i;//从i往前扫
			while (out.seq[k].contains(out.splits[j], out.splits[j + 1]))
			{
				term = multiIntegral(out.pr[j]+k, i-k+1, out.splits[j], out.splits[j + 1]);
				//term = term * pr[k] * gap / (i - k + 1);
				sum += k > 0 ? term * out.DP[k - 1][j - 1] : term;
				k--;
				if (k < 0) break;
			}
			out.DP[i][j] = sum;
		}
	}
}

/*//以下函数废弃：使用按行单步计算的方法，动态增加subrange
//此外这些函数没有去除重复的split点
void setSubRanges(vector<double> &output, vector<Range> &input)
{
	//input为一组range
	//output为由range决定的子range的split points
	output.clear();
	output.reserve(2*input.size());
	vector<Range>::iterator it;
	for(it=input.begin(); it!=input.end(); it++)
	{
		output.push_back(it->l);
		output.push_back(it->r);
	}
	sort(output.begin(), output.end());
}

void setSubRanges(vector<double> &output, Range* row, int col)
{
	//input为矩阵的一行, col为列数
	//output为由range决定的子range的split points
	output.clear();
	output.reserve(2*col);
	for(int i=0; i<col; i++)
	{
		output.push_back(row[i].l);
		output.push_back(row[i].r);
	}
	sort(output.begin(), output.end());
}

//以下函数直接计算DP (naive)
void DP(DPMatrix& re, vector<Range>& ranges)
{
	vector<double> splits;
	setSubRanges(splits, ranges);
	//以下设置re.splits
	re.rngNum=splits.size()-1;
	re.splits=new double[re.rngNum+1];
	vector<double>::iterator it;
	int pos=0;
	for(it=splits.begin(); it!=splits.end(); it++)
	{
		re.splits[pos]=*it;
		pos++;
	}
	//以下设置re.seq
	re.seqLen=ranges.size();
	re.seq=new Range[re.seqLen];
	for(int i=0; i<re.seqLen; i++) re.seq[i]=ranges[i];
	//以下获取概率数组(对于每个ti一个entry)
	double pr[re.seqLen];
	for(int i=0; i<re.seqLen; i++) pr[i]=re.seq[i].pdf();
	//以下设置re.DP
	re.DP=new double*[re.seqLen];
	for(int i=0; i<re.seqLen; i++) re.DP[i]=new double[re.rngNum];
	//计算
	for(int i=0; i<re.seqLen; i++)
	{
		for(int j=0; j<re.rngNum; j++)
		{
			Range cur(re.splits[j], re.splits[j+1]);
			double gap=cur.gap();
			if(j==0)//base case
			{
				double term=1;
				int k=i;//从i往前扫
				while(re.seq[k].contains(cur) && pr[k]>0)//在cur区间上pdf不为0
				{
					term = term*pr[k]*gap/(i-k+1);
					k--;
					if(k<0) break;
				}
				if(k<0) re.DP[i][j]=term;
				else re.DP[i][j]=0;
			}
			else
			{
				//以下计算概率部分
				double sum=re.DP[i][j-1];
				double term=1;//递推项：subrange上的概率
				int k=i;//从i往前扫
				while(re.seq[k].contains(cur) && pr[k]>0)//在cur区间上pdf不为0
				{
					term = term*pr[k]*gap/(i-k+1);
					sum += k>0?term*re.DP[k-1][j-1]:term;
					k--;
					if(k<0) break;
				}
				re.DP[i][j]=sum;
			}
		}
	}
}
//*/

//=============== 以下为DFS Expected Support ===============
Range rngOP(Range& lastRng, Range& nextRng)
{
	double l1 = lastRng.l;
	double l2 = nextRng.l;
	return Range(l1 > l2 ? l1 : l2, nextRng.r);
}

void DFS_ExpSup(vector<int>& seq, vector<Row>& db, vector<int>& dict, Matrix& mat, int minrow, int mincol, ofstream& fout)
{//检查当前seq是否frequent, 然后递归
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*
	cout<<"Checking seq = ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<"..."<<endl;
	*/
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//mat为数据矩阵
	vector<Row> db1;
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1.push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();
		}
	}
	if (sum0 < minrow) return;//CntPrune
	double sum = 0;//用于refinement
	int cnt = db1.size();
	for (int i = 0; i < cnt; i++)
	{
		db1[i].mat = new DPMatrix();
		appendRng(*db1[i].mat, *db1_pt[i].mat, mat.M[db1[i].rowID][col]);//必须seq长度>1
		sum += db1[i].mat->getPr();//概率累加

		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		//by wenwen : here logic???
		if (sum0 + sum < minrow)//Partial CntPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete db1[j].mat;
			}
			return;
		}
	}
	if (sum >= minrow)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputExpSup(seq, db1, sum, fout);//输出pattern
		vector<int> dict1;
		//构造新字典: 去掉当前最后字符
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			int cur = *it;
			if (cur != col) dict1.push_back(cur);
		}
		for (vector<int>::iterator it = dict1.begin(); it != dict1.end(); it++)
		{
			//by wenwen: can we use seq continuely, first push_back, then pop to avoid copy?
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ExpSup(seq1, db1, dict1, mat, minrow, mincol, fout);//递归
		}
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1.begin(); it != db1.end(); it++)
	{
		delete (*it).mat;
	}
}

void DFS_ExpSup(Matrix &mat, int minrow, int mincol, ofstream& fout)
{//外部调用接口
	//对于matrix格式，没有missing data，因此长度为1的都frequent
	//这里假设 (1) minrow <= row; (2) mincol > 1
	int row = mat.row;
	int col = mat.col;

	for (int i = 0; i < col; i++)
	{
		//cout<<"Processing "<<i<<" ..."<<endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vector<int> seq; // by wenwen: why use seq here?
		seq.push_back(i);

		vector<Row> db;
		db.reserve(row);
		for (int j = 0; j < row; j++)
		{
			Row cur;
			cur.rowID = j;
			cur.lastRng = mat.M[j][i];
			cur.mat = new DPMatrix(cur.lastRng);
			db.push_back(cur);
		}

		vector<int> dict;
		//构造新字典: 去掉当前最后字符
		for (int j = 0; j < col; j++)
		{
			if (j != i) dict.push_back(j);
		}
		//递归
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ExpSup(seq1, db, dict, mat, minrow, mincol, fout);
		}
		//释放db
		for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
		{
			delete (*it).mat;
		}
	}
}

//=============== 以下为DFS Probablisticaly Frequent ===============

void DFS_ProbFreq(vector<int>& seq, vector<Row>& db, vector<int>& dict, Matrix& mat, int minrow, int mincol, double th_prob, ofstream& fout)
{//检查当前seq是否frequent, 然后递归
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*
	cout<<"Checking seq = ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<"..."<<endl;
	//*/
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double prod = minrow * th_prob;
	//mat为数据矩阵
	vector<Row> db1;
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1.push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();
		}
	}
	double sum = 0;//用于refinement
	int cnt = db1.size();
	if (cnt < minrow) return;//CntPrune
	if (sum0 < prod) return;//MarkovPrune
	double* vec = new double[cnt];//用于算convolution的向量
	for (int i = 0; i < cnt; i++)
	{
		db1[i].mat = new DPMatrix();
		appendRng(*db1[i].mat, *db1_pt[i].mat, mat.M[db1[i].rowID][col]);//必须seq长度>1
		vec[i] = db1[i].mat->getPr();
		sum += vec[i];//概率累加
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + sum < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放convolution计算用的向量
			delete[] vec;
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete db1[j].mat;
			}
			return;
		}
	}

	//BEGIN: ExpPrune
	double delta = (minrow - sum - 1) / sum;

	if (delta > 0)
	{
		bool expPrune = false;
		if (delta >= 2 * EXP - 1)
		{
			if (pow(2, -sum * delta) < th_prob) expPrune = true;
		}
		else
		{
			if (exp(-sum * delta*delta / 4) < th_prob) expPrune = true;
		}
		if (expPrune)
		{
			//cout<<"ExpPruned"<<endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//释放convolution计算用的向量
			delete[] vec;
			//手动释放db1的DP矩阵等信息
			for (vector<Row>::iterator it = db1.begin(); it != db1.end(); it++)
			{
				delete (*it).mat;
			}
			return;
		}
	}
	//END: ExpPrune

	//计算是否probFreq
	double *buf;
	bool freq = PMFCheck(buf, vec, cnt, minrow, th_prob);
	//delete[] buf;

	if (freq)
	{
		if (seq.size() >= mincol) outputFreq(seq, db1, th_prob, fout);//输出pattern
		vector<int> dict1;
		//构造新字典: 去掉当前最后字符
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			int cur = *it;
			if (cur != col) dict1.push_back(cur);
		}
		for (vector<int>::iterator it = dict1.begin(); it != dict1.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreq(seq1, db1, dict1, mat, minrow, mincol, th_prob, fout);//递归
		}
	}
	//释放convolution计算用的向量
	delete[] vec;
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1.begin(); it != db1.end(); it++)
	{
		delete (*it).mat;
	}
}

void DFS_ProbFreq(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{//外部调用接口
	//对于matrix格式，没有missing data，因此长度为1的都frequent
	//这里假设 (1) minrow <= row; (2) mincol > 1
	int row = mat.row;
	int col = mat.col;

	for (int i = 0; i < col; i++)
	{
		//cout<<"Processing "<<i<<" ..."<<endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vector<int> seq;
		seq.push_back(i);

		vector<Row> db;
		db.reserve(row);
		for (int j = 0; j < row; j++)
		{
			Row cur;
			cur.rowID = j;
			cur.lastRng = mat.M[j][i];
			cur.mat = new DPMatrix(cur.lastRng);
			db.push_back(cur);
		}

		vector<int> dict;
		//构造新字典: 去掉当前最后字符
		for (int j = 0; j < col; j++)
		{
			if (j != i) dict.push_back(j);
		}
		//递归
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreq(seq1, db, dict, mat, minrow, mincol, th_prob, fout);
		}
		//释放db
		for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
		{
			delete (*it).mat;
		}
	}
}

//=============== 以下为BFS Expected Support ===============

void bottomCheck_ExpSup(vector<int> seq, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, ofstream& fout)
{//处理当前pattern，输出frequent的pattern，并将该pattern插入nextRoot
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*
	cout<<"Checking seq = ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<"..."<<endl;
	//*/
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//这里是新的check，需要各种pruning
	if (subseqExists(curRoot, seq) == false) return;//首先进行subsequence pruning
	///////////////////////////////
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1->push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();//利用上次的pr作为上界
		}
	}
	if (sum0 < minrow)//CntPrune (利用上次的pr作为上界)
	{
		delete db1;
		return;
	}

	int cnt = db1->size();
	double sum = 0;
	for (int i = 0; i < cnt; i++)
	{
		(*db1)[i].mat = new DPMatrix();
		appendRng(*(*db1)[i].mat, *db1_pt[i].mat, mat.M[(*db1)[i].rowID][col]);//必须seq长度>1
		sum += (*db1)[i].mat->getPr();//概率累加
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + sum < minrow)//Partial CntPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete (*db1)[j].mat;
			}
			delete db1;
			return;
		}
	}
	if (sum >= minrow)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputExpSup(seq, *db1, sum, fout);//输出pattern
		addSeq(nextRoot, seq);
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void recursiveCheck_ExpSup(int k, int level, vector<int> seq, Node *curNode, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, vector<int>* dict, ofstream& fout)
{//处理当前节点并递归
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*
	cout<<"Checking seq = ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<"..."<<endl;
	//*/
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//recursiveCheck走的是以前的老路(根据curNode->ch)，不需要pruning(也pruning不掉)
	//同样的，不需要输出pattern
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			row1.mat = new DPMatrix();
			appendRng(*(row1.mat), *(it->mat), mat.M[row][col]);
			db1->push_back(row1);
		}
	}
	vector<Node*> chList = curNode->chList;
	//递归
	if (level < k - 2)
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ExpSup(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, fout);
			dbStack.pop();
		}
	}
	else if (level == k - 2)
	{
		//收集dict
		vector<int> *dict = new vector<int>;
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			dict->push_back(child->last);
		}
		//递归
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ExpSup(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, fout);
			dbStack.pop();
		}
		//k-2层释放dict，后面的不用管
		delete dict;
	}
	else if (level == k - 1)
	{
		//递归
		for (vector<int>::iterator it = dict->begin(); it != dict->end(); it++)
		{
			if ((*it) != seq.back())//去掉已经出现的列
			{
				vector<int> seq1 = seq;
				seq1.push_back(*it);
				dbStack.push(db1);
				bottomCheck_ExpSup(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, fout);
				dbStack.pop();
			}
		}
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void Apriori_ExpSup(Matrix &mat, int minrow, int mincol, ofstream& fout)
{
	int row = mat.row;
	int col = mat.col;

	//Length 1 tree: 由于假设没有missing value而且minrow<=row，不用检查是否frequent
	Node *curRoot = initTree(), *nextRoot;
	for (int i = 0; i < col; i++)
	{
		curRoot->appendChild(new Node(i));
	}
	//pattern growth
	for (int k = 2; !emptyTree(curRoot); k++)
	{
		//printTree(curRoot);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//cout<<"Generating Length "<<k<<" Patterns ..."<<endl;
		//初始化nextRoot
		nextRoot = initTree();
		//从curRoot生成nextRoot
		stack<vector<Row>*> dbStack;
		//处理第一层
		for (vector<Node*>::iterator it = curRoot->chList.begin(); it != curRoot->chList.end(); it++)
		{
			int i = (*it)->last;
			vector<int> seq;
			seq.push_back(i);

			vector<Row>* db = new vector<Row>;
			db->reserve(row);
			for (int j = 0; j < row; j++)
			{
				Row cur;
				cur.rowID = j;
				cur.lastRng = mat.M[j][i];
				cur.mat = new DPMatrix(cur.lastRng);
				db->push_back(cur);
			}
			vector<int> *dict = new vector<int>;
			if (k == 2)
			{
				for (int j = 0; j < col; j++)
				{
					if (j != i)
					{
						vector<int> seq1;
						seq1.push_back(i);
						seq1.push_back(j);
						dbStack.push(db);
						bottomCheck_ExpSup(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, fout);
						dbStack.pop();
					}
				}
			}
			else if (k == 3)
			{
				//构造字典
				vector<int> *dict = new vector<int>;
				for (vector<Node*>::iterator dit = (*it)->chList.begin(); dit != (*it)->chList.end(); dit++)
				{
					Node* child = *dit;
					dict->push_back(child->last);
				}
				//利用dict递归
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ExpSup(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, fout);
					dbStack.pop();
				}
				delete dict;
			}
			else
			{
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ExpSup(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, fout);
					dbStack.pop();
				}
			}
			//释放db
			for (vector<Row>::iterator it1 = (*db).begin(); it1 != (*db).end(); it1++)
			{
				delete (*it1).mat;
			}
			delete db;
		}
		//更新当前tree
		delete curRoot;
		curRoot = nextRoot;
	}
	delete curRoot;
}

//=============== 以下为BFS Probablistically Frequent ===============

void bottomCheck_ProbFreq(vector<int> seq, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, double th_prob, ofstream& fout)
{//处理当前pattern，输出frequent的pattern，并将该pattern插入nextRoot
	//这里是新的check，需要各种pruning
	if (subseqExists(curRoot, seq) == false) return;//首先进行subsequence pruning
	///////////////////////////////
	double prod = minrow * th_prob;
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1->push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();//利用上次的pr作为上界
		}
	}
	double sum = 0;//用于refinement
	int cnt = db1->size();
	if (cnt < minrow)//CntPrune
	{
		delete db1;
		return;
	}
	if (sum0 < prod)//MarkovPrune
	{
		delete db1;
		return;
	}
	double* vec = new double[cnt];//用于算convolution的向量
	for (int i = 0; i < cnt; i++)
	{
		(*db1)[i].mat = new DPMatrix();
		appendRng(*(*db1)[i].mat, *db1_pt[i].mat, mat.M[(*db1)[i].rowID][col]);//必须seq长度>1
		vec[i] = (*db1)[i].mat->getPr();
		sum += vec[i];//概率累加
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + sum < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放convolution计算用的向量
			delete[] vec;
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete (*db1)[j].mat;
			}
			delete db1;
			return;
		}
	}

	//BEGIN: ExpPrune
	double delta = (minrow - sum - 1) / sum;

	if (delta > 0)
	{
		bool expPrune = false;
		if (delta >= 2 * EXP - 1)
		{
			if (pow(2, -sum * delta) < th_prob) expPrune = true;
		}
		else
		{
			if (exp(-sum * delta*delta / 4) < th_prob) expPrune = true;
		}
		if (expPrune)
		{
			//释放convolution计算用的向量
			delete[] vec;
			//手动释放db1的DP矩阵等信息
			for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
			{
				delete (*it).mat;
			}
			delete db1;
			return;
		}
	}
	//END: ExpPrune

	//计算是否probFreq
	double *buf;
	bool freq = PMFCheck(buf, vec, cnt, minrow, th_prob);
	//	delete[] buf;

	if (freq)
	{
		if (seq.size() >= mincol) outputFreq(seq, *db1, th_prob, fout);//输出pattern
		addSeq(nextRoot, seq);
	}

	//释放convolution计算用的向量
	delete[] vec;
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void recursiveCheck_ProbFreq(int k, int level, vector<int> seq, Node *curNode, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, vector<int>* dict, double th_prob, ofstream& fout)
{//处理当前节点并递归
	//recursiveCheck走的是以前的老路(根据curNode->ch)，不需要pruning(也pruning不掉)
	//同样的，不需要输出pattern
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			row1.mat = new DPMatrix();
			appendRng(*(row1.mat), *(it->mat), mat.M[row][col]);
			db1->push_back(row1);
		}
	}
	vector<Node*> chList = curNode->chList;
	//递归
	if (level < k - 2)
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreq(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, fout);
			dbStack.pop();
		}
	}
	else if (level == k - 2)
	{
		//收集dict
		vector<int> *dict = new vector<int>;
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			dict->push_back(child->last);
		}
		//递归
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreq(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, fout);
			dbStack.pop();
		}
		//k-2层释放dict，后面的不用管
		delete dict;
	}
	else if (level == k - 1)
	{
		//递归
		for (vector<int>::iterator it = dict->begin(); it != dict->end(); it++)
		{
			if ((*it) != seq.back())//去掉已经出现的列
			{
				vector<int> seq1 = seq;
				seq1.push_back(*it);
				dbStack.push(db1);
				bottomCheck_ProbFreq(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, fout);
				dbStack.pop();
			}
		}
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void Apriori_ProbFreq(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{
	int row = mat.row;
	int col = mat.col;

	//Length 1 tree: 由于假设没有missing value而且minrow<=row，不用检查是否frequent
	Node *curRoot = initTree(), *nextRoot;
	for (int i = 0; i < col; i++)
	{
		curRoot->appendChild(new Node(i));
	}
	//pattern growth
	for (int k = 2; !emptyTree(curRoot); k++)
	{
		//cout<<"Generating Length "<<k<<" Patterns ..."<<endl;
		//初始化nextRoot
		nextRoot = initTree();
		//从curRoot生成nextRoot
		stack<vector<Row>*> dbStack;
		//处理第一层
		for (vector<Node*>::iterator it = curRoot->chList.begin(); it != curRoot->chList.end(); it++)
		{
			int i = (*it)->last;
			vector<int> seq;
			seq.push_back(i);

			vector<Row>* db = new vector<Row>;
			db->reserve(row);
			for (int j = 0; j < row; j++)
			{
				Row cur;
				cur.rowID = j;
				cur.lastRng = mat.M[j][i];
				cur.mat = new DPMatrix(cur.lastRng);
				db->push_back(cur);
			}
			if (k == 2)
			{
				for (int j = 0; j < col; j++)
				{
					if (j != i)
					{
						vector<int> seq1;
						seq1.push_back(i);
						seq1.push_back(j);
						dbStack.push(db);
						bottomCheck_ProbFreq(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, fout);
						dbStack.pop();
					}
				}
			}
			else if (k == 3)
			{
				//构造字典
				vector<int> *dict = new vector<int>;
				for (vector<Node*>::iterator dit = (*it)->chList.begin(); dit != (*it)->chList.end(); dit++)
				{
					Node* child = *dit;
					dict->push_back(child->last);
				}
				//利用dict递归
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreq(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, fout);
					dbStack.pop();
				}
				delete dict;
			}
			else
			{
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreq(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, fout);
					dbStack.pop();
				}
			}
			//释放db
			for (vector<Row>::iterator it1 = (*db).begin(); it1 != (*db).end(); it1++)
			{
				delete (*it1).mat;
			}
			delete db;
		}
		//更新当前tree
		delete curRoot;
		curRoot = nextRoot;
	}
	delete curRoot;
}

//comment: 只有bottomCheck函数是不同的
//db排序只对pruning+输出有用：目前只对dfs函数进行了排序


//------------------------------------------approximate PMF

void bottomCheck_ProbFreqApprox(vector<int> seq, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, double th_prob, double tm, ofstream& fout)
{//处理当前pattern，输出frequent的pattern，并将该pattern插入nextRoot
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //这里是新的check，需要各种pruning
	if (subseqExists(curRoot, seq) == false) return;//首先进行subsequence pruning
													///////////////////////////////
	double prod = minrow * th_prob;

	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1->push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();//利用上次的pr作为上界
		}
	}
	if (sum0 < minrow)//CntPrune (利用上次的pr作为上界)
	{
		delete db1;
		return;
	}
	if (sum0 < prod)
	{
		delete db1;
		return;
	}

	int cnt = db1->size();
	double mu = 0;
	double segmaSquare = 0;
	double tmpProb;
	for (int i = 0; i < cnt; i++)
	{
		(*db1)[i].mat = new DPMatrix();
		appendRng(*(*db1)[i].mat, *db1_pt[i].mat, mat.M[(*db1)[i].rowID][col]);//必须seq长度>1
		tmpProb = (*db1)[i].mat->getPr();
		mu += tmpProb;//概率累加
		segmaSquare = segmaSquare + tmpProb * (1 - tmpProb);
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + mu < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete (*db1)[j].mat;
			}
			delete db1;
			return;
		}
	}
	double t = (minrow - 0.5 - mu) / sqrt(segmaSquare);

	if (t <= tm)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputFreq(seq, *db1, th_prob, fout);//输出pattern
		addSeq(nextRoot, seq);
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void recursiveCheck_ProbFreqApprox(int k, int level, vector<int> seq, Node *curNode, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, vector<int>* dict, double th_prob, double tm, ofstream& fout)
{//处理当前节点并递归
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //recursiveCheck走的是以前的老路(根据curNode->ch)，不需要pruning(也pruning不掉)
 //同样的，不需要输出pattern
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			row1.mat = new DPMatrix();
			appendRng(*(row1.mat), *(it->mat), mat.M[row][col]);
			db1->push_back(row1);
		}
	}
	vector<Node*> chList = curNode->chList;
	//递归
	if (level < k - 2)
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreqApprox(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, tm, fout);
			dbStack.pop();
		}
	}
	else if (level == k - 2)
	{
		//收集dict
		vector<int> *dict = new vector<int>;
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			dict->push_back(child->last);
		}
		//递归
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreqApprox(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, tm, fout);
			dbStack.pop();
		}
		//k-2层释放dict，后面的不用管
		delete dict;
	}
	else if (level == k - 1)
	{
		//递归
		for (vector<int>::iterator it = dict->begin(); it != dict->end(); it++)
		{
			if ((*it) != seq.back())//去掉已经出现的列
			{
				vector<int> seq1 = seq;
				seq1.push_back(*it);
				dbStack.push(db1);
				bottomCheck_ProbFreqApprox(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, tm, fout);
				dbStack.pop();
			}
		}
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}



void Apriori_ProbFreqApprox(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{
	int row = mat.row;
	int col = mat.col;

	double tm = NormalCDFInverse((1 - th_prob) / sqrt(2 * PI));
	//cout << tm << endl;

	//Length 1 tree: 由于假设没有missing value而且minrow<=row，不用检查是否frequent
	Node *curRoot = initTree(), *nextRoot;
	for (int i = 0; i < col; i++)
	{
		curRoot->appendChild(new Node(i));
	}
	//pattern growth
	for (int k = 2; !emptyTree(curRoot); k++)
	{
		//printTree(curRoot);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//cout << "Generating Length " << k << " Patterns ..." << endl;
		//初始化nextRoot
		nextRoot = initTree();
		//从curRoot生成nextRoot
		stack<vector<Row>*> dbStack;
		//处理第一层
		for (vector<Node*>::iterator it = curRoot->chList.begin(); it != curRoot->chList.end(); it++)
		{
			int i = (*it)->last;
			vector<int> seq;
			seq.push_back(i);

			vector<Row>* db = new vector<Row>;
			db->reserve(row);
			for (int j = 0; j < row; j++)
			{
				Row cur;
				cur.rowID = j;
				cur.lastRng = mat.M[j][i];
				cur.mat = new DPMatrix(cur.lastRng);
				db->push_back(cur);
			}
			vector<int> *dict = new vector<int>;
			if (k == 2)
			{
				for (int j = 0; j < col; j++)
				{
					if (j != i)
					{
						vector<int> seq1;
						seq1.push_back(i);
						seq1.push_back(j);
						dbStack.push(db);
						bottomCheck_ProbFreqApprox(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, tm, fout);
						dbStack.pop();
					}
				}
			}
			else if (k == 3)
			{
				//构造字典
				vector<int> *dict = new vector<int>;
				for (vector<Node*>::iterator dit = (*it)->chList.begin(); dit != (*it)->chList.end(); dit++)
				{
					Node* child = *dit;
					dict->push_back(child->last);
				}
				//利用dict递归
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreqApprox(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, tm, fout);
					dbStack.pop();
				}
				delete dict;
			}
			else
			{
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreqApprox(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, tm, fout);
					dbStack.pop();
				}
			}
			//释放db
			for (vector<Row>::iterator it1 = (*db).begin(); it1 != (*db).end(); it1++)
			{
				delete (*it1).mat;
			}
			delete db;
		}
		//更新当前tree
		delete curRoot;
		curRoot = nextRoot;
	}
	delete curRoot;
}


//--------------------------------------------DFS_Prob_Freq_Approximation

void DFS_ProbFreqApprox(vector<int>& seq, vector<Row>& db, vector<int>& dict, Matrix& mat, int minrow, int mincol, double th_prob, double tm, ofstream& fout)
{//检查当前seq是否frequent, 然后递归
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double prod = minrow * th_prob;
	//mat为数据矩阵
	vector<Row> db1;
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1.push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();
		}
	}
	//double sum = 0;//用于refinement
	double mu = 0;
	double segmaSquare = 0;
	double tmpProb;
	int cnt = db1.size();
	if (cnt < minrow) return;//CntPrune
	if (sum0 < prod) return;//MarkovPrune
//	double* vec = new double[cnt];//用于算convolution的向量
	for (int i = 0; i < cnt; i++)
	{
		db1[i].mat = new DPMatrix();
		appendRng(*db1[i].mat, *db1_pt[i].mat, mat.M[db1[i].rowID][col]);//必须seq长度>1


		tmpProb = (db1)[i].mat->getPr();
		mu += tmpProb;//概率累加
		segmaSquare = segmaSquare + tmpProb * (1 - tmpProb);
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + mu < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放convolution计算用的向量
			//delete[] vec;
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete db1[j].mat;
			}
			return;
		}
	}
	double t = (minrow - 0.5 - mu) / sqrt(segmaSquare);
	if (t <= tm)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputFreq(seq, db1, th_prob, fout);//输出pattern
		vector<int> dict1;
		//构造新字典: 去掉当前最后字符
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			int cur = *it;
			if (cur != col) dict1.push_back(cur);
		}
		for (vector<int>::iterator it = dict1.begin(); it != dict1.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreqApprox(seq1, db1, dict1, mat, minrow, mincol, th_prob, tm, fout);//递归
		}
	}

	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1.begin(); it != db1.end(); it++)
	{
		delete (*it).mat;
	}
}

void DFS_ProbFreqApprox(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{//外部调用接口
 //对于matrix格式，没有missing data，因此长度为1的都frequent
 //这里假设 (1) minrow <= row; (2) mincol > 1
	int row = mat.row;
	int col = mat.col;
	double mu_m = NormalCDFInverse((1 - th_prob) / sqrt(2 * PI));
	for (int i = 0; i < col; i++)
	{
		//cout << "Processing " << i << " ..." << endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vector<int> seq;
		seq.push_back(i);

		vector<Row> db;
		db.reserve(row);
		for (int j = 0; j < row; j++)
		{
			Row cur;
			cur.rowID = j;
			cur.lastRng = mat.M[j][i];
			cur.mat = new DPMatrix(cur.lastRng);
			db.push_back(cur);
		}

		vector<int> dict;
		//构造新字典: 去掉当前最后字符
		for (int j = 0; j < col; j++)
		{
			if (j != i) dict.push_back(j);
		}
		//递归
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreqApprox(seq1, db, dict, mat, minrow, mincol, th_prob, mu_m, fout);
		}
		//释放db
		for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
		{
			delete (*it).mat;
		}
	}
}





//--------------------------------------------DFS_Prob_Freq_Approximation_poi

void DFS_ProbFreqApproxPoi(vector<int>& seq, vector<Row>& db, vector<int>& dict, Matrix& mat, int minrow, int mincol, double th_prob, double mu_m, ofstream& fout)
{//检查当前seq是否frequent, 然后递归
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	double prod = minrow * th_prob;
	//mat为数据矩阵
	vector<Row> db1;
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1.push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();
		}
	}
	//double sum = 0;//用于refinement
	double mu = 0;
	double segmaSquare = 0;
	double tmpProb;
	int cnt = db1.size();
	if (cnt < minrow) return;//CntPrune
	if (sum0 < prod) return;//MarkovPrune
	//double* vec = new double[cnt];//用于算convolution的向量
	for (int i = 0; i < cnt; i++)
	{
		db1[i].mat = new DPMatrix();
		appendRng(*db1[i].mat, *db1_pt[i].mat, mat.M[db1[i].rowID][col]);//必须seq长度>1


		tmpProb = (db1)[i].mat->getPr();
		mu += tmpProb;//概率累加
		segmaSquare = segmaSquare + tmpProb * (1 - tmpProb);
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + mu < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放convolution计算用的向量
			//delete[] vec;
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete db1[j].mat;
			}
			return;
		}
	}
	if (mu >= mu_m)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputFreq(seq, db1, th_prob, fout);//输出pattern
		vector<int> dict1;
		//构造新字典: 去掉当前最后字符
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			int cur = *it;
			if (cur != col) dict1.push_back(cur);
		}
		for (vector<int>::iterator it = dict1.begin(); it != dict1.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreqApproxPoi(seq1, db1, dict1, mat, minrow, mincol, th_prob, mu_m, fout);//递归
		}
	}

	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1.begin(); it != db1.end(); it++)
	{
		delete (*it).mat;
	}
}

void DFS_ProbFreqApproxPoi(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{//外部调用接口
 //对于matrix格式，没有missing data，因此长度为1的都frequent
 //这里假设 (1) minrow <= row; (2) mincol > 1
	int row = mat.row;
	int col = mat.col;
	double mu_m = 0;


	//*****************solve equation for mu_m 


	double eps = 10.0e-6;
	
	double x = 1;
	
	if (abs(1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) - th_prob) < eps*th_prob) {
		mu_m = x;
	}
	else {
		double ub = 0;
		double lb = 0;

		if (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) > th_prob) {
			while (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) > th_prob) {

				ub = x;
				x = x / 2;

			}
			lb = x;
		}
		else {
			while (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) < th_prob) {
				lb = x;
				x = x * 2;
			}
			ub = x;
		}


		double mid = mid = (lb + ub) / 2;
		int counter = 0;
		while (abs(1 - cdf(boost::math::poisson_distribution<>(mid), minrow-1) - th_prob) > eps*th_prob) {
			if (1 - cdf(boost::math::poisson_distribution<>(mid), minrow-1) > th_prob) {
				ub = mid;
			}
			else {
				lb = mid;
			}
			mid = (lb + ub) / 2;
		}

		mu_m = mid;
	}
	//*****************end of solve equation for mu_m





	for (int i = 0; i < col; i++)
	{
		//cout << "Processing " << i << " ..." << endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vector<int> seq;
		seq.push_back(i);

		vector<Row> db;
		db.reserve(row);
		for (int j = 0; j < row; j++)
		{
			Row cur;
			cur.rowID = j;
			cur.lastRng = mat.M[j][i];
			cur.mat = new DPMatrix(cur.lastRng);
			db.push_back(cur);
		}

		vector<int> dict;
		//构造新字典: 去掉当前最后字符
		for (int j = 0; j < col; j++)
		{
			if (j != i) dict.push_back(j);
		}
		//递归
		for (vector<int>::iterator it = dict.begin(); it != dict.end(); it++)
		{
			vector<int> seq1 = seq;
			seq1.push_back(*it);
			DFS_ProbFreqApproxPoi(seq1, db, dict, mat, minrow, mincol, th_prob, mu_m, fout);
		}
		//释放db
		for (vector<Row>::iterator it = db.begin(); it != db.end(); it++)
		{
			delete (*it).mat;
		}
	}
}



//------------------------------------------approximate PMF POI

void bottomCheck_ProbFreqApproxPoi(vector<int> seq, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, double th_prob, double mu_m, ofstream& fout)
{//处理当前pattern，输出frequent的pattern，并将该pattern插入nextRoot
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //这里是新的check，需要各种pruning
	if (subseqExists(curRoot, seq) == false) return;//首先进行subsequence pruning
													///////////////////////////////
	double prod = minrow * th_prob;

	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	vector<Row> db1_pt;//记录db1中记录对应的db中记录
	double sum0 = 0;
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			db1->push_back(row1);
			db1_pt.push_back(*it);
			sum0 += it->mat->getPr();//利用上次的pr作为上界
		}
	}
	if (sum0 < minrow)//CntPrune (利用上次的pr作为上界)
	{
		delete db1;
		return;
	}
	if (sum0 < prod)
	{
		delete db1;
		return;
	}

	int cnt = db1->size();
	double mu = 0;
	double segmaSquare = 0;
	double tmpProb;
	for (int i = 0; i < cnt; i++)
	{
		(*db1)[i].mat = new DPMatrix();
		appendRng(*(*db1)[i].mat, *db1_pt[i].mat, mat.M[(*db1)[i].rowID][col]);//必须seq长度>1
		tmpProb = (*db1)[i].mat->getPr();
		mu += tmpProb;//概率累加
		segmaSquare = segmaSquare + tmpProb * (1 - tmpProb);
		sum0 -= db1_pt[i].mat->getPr();
		//虽然可以检查是否sum>=minrow，但是为了递归时db1正确设定，不进行这步pruning
		if (sum0 + mu < prod)//MarkovPrune, 这个是可以的，因为后面db1用不到了
		{
			//释放已经分配的mat
			for (int j = 0; j <= i; j++)
			{
				delete (*db1)[j].mat;
			}
			delete db1;
			return;
		}
	}
	
	if (mu > mu_m)//expSup>=minrow?
	{
		if (seq.size() >= mincol) outputFreq(seq, *db1, th_prob, fout);//输出pattern
		addSeq(nextRoot, seq);
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}

void recursiveCheck_ProbFreqApproxPoi(int k, int level, vector<int> seq, Node *curNode, Node *curRoot, Node *nextRoot, stack<vector<Row>*>& dbStack, Matrix& mat, int minrow, int mincol, vector<int>* dict, double th_prob, double mu_m, ofstream& fout)
{//处理当前节点并递归
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 /*
 cout<<"Checking seq = ";
 for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
 {
 cout<<*it<<" ";
 }
 cout<<"..."<<endl;
 //*/
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //recursiveCheck走的是以前的老路(根据curNode->ch)，不需要pruning(也pruning不掉)
 //同样的，不需要输出pattern
	vector<Row>* db = dbStack.top();
	vector<Row>* db1 = new vector<Row>();
	int col = seq.back();
	for (vector<Row>::iterator it = (*db).begin(); it != (*db).end(); it++)
	{
		int row = it->rowID;
		Range lastRng = mat.M[row][col];
		//用lastRng进行pruning
		Range curRng = rngOP(it->lastRng, lastRng);
		if (curRng.r > curRng.l)
		{
			//range不为空(没有等号, 即使等于, 概率也为0)
			Row row1;
			row1.rowID = row;
			row1.lastRng = curRng;
			row1.mat = new DPMatrix();
			appendRng(*(row1.mat), *(it->mat), mat.M[row][col]);
			db1->push_back(row1);
		}
	}
	vector<Node*> chList = curNode->chList;
	//递归
	if (level < k - 2)
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreqApproxPoi(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, mu_m, fout);
			dbStack.pop();
		}
	}
	else if (level == k - 2)
	{
		//收集dict
		vector<int> *dict = new vector<int>;
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			dict->push_back(child->last);
		}
		//递归
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			Node* child = *it;
			vector<int> seq1 = seq;
			seq1.push_back(child->last);
			dbStack.push(db1);
			recursiveCheck_ProbFreqApproxPoi(k, level + 1, seq1, child, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, mu_m, fout);
			dbStack.pop();
		}
		//k-2层释放dict，后面的不用管
		delete dict;
	}
	else if (level == k - 1)
	{
		//递归
		for (vector<int>::iterator it = dict->begin(); it != dict->end(); it++)
		{
			if ((*it) != seq.back())//去掉已经出现的列
			{
				vector<int> seq1 = seq;
				seq1.push_back(*it);
				dbStack.push(db1);
				bottomCheck_ProbFreqApproxPoi(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, mu_m, fout);
				dbStack.pop();
			}
		}
	}
	//手动释放db1的DP矩阵等信息
	for (vector<Row>::iterator it = db1->begin(); it != db1->end(); it++)
	{
		delete (*it).mat;
	}
	delete db1;
}



void Apriori_ProbFreqApproxPoi(Matrix &mat, int minrow, int mincol, double th_prob, ofstream& fout)
{

	// interface for main function 
	int row = mat.row;
	int col = mat.col;


	double mu_m = 0;


	//*****************solve equation for mu_m 


	double eps = 10.0e-6;

	double x = 1;

	if (abs(1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) - th_prob) < eps*th_prob) {
		mu_m = x;
	}
	else {
		double ub = 0;
		double lb = 0;

		if (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) > th_prob) {
			while (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) > th_prob) {

				ub = x;
				x = x / 2;

			}
			lb = x;
		}
		else {
			while (1 - cdf(boost::math::poisson_distribution<>(x), minrow-1) < th_prob) {
				lb = x;
				x = x * 2;
			}
			ub = x;
		}


		double mid = mid = (lb + ub) / 2;
		int counter = 0;
		while (abs(1 - cdf(boost::math::poisson_distribution<>(mid), minrow-1) - th_prob) > eps*th_prob) {
			if (1 - cdf(boost::math::poisson_distribution<>(mid), minrow-1) > th_prob) {
				ub = mid;
			}
			else {
				lb = mid;
			}
			mid = (lb + ub) / 2;
		}

		mu_m = mid;
	}
	//*****************end of solve equation for mu_m


	//Length 1 tree: 由于假设没有missing value而且minrow<=row，不用检查是否frequent
	Node *curRoot = initTree(), *nextRoot;
	for (int i = 0; i < col; i++)
	{
		curRoot->appendChild(new Node(i));
	}
	//pattern growth
	for (int k = 2; !emptyTree(curRoot); k++)
	{
		//printTree(curRoot);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//cout << "Generating Length " << k << " Patterns ..." << endl;
		//初始化nextRoot
		nextRoot = initTree();
		//从curRoot生成nextRoot
		stack<vector<Row>*> dbStack;
		//处理第一层
		for (vector<Node*>::iterator it = curRoot->chList.begin(); it != curRoot->chList.end(); it++)
		{
			int i = (*it)->last;
			vector<int> seq;
			seq.push_back(i);

			vector<Row>* db = new vector<Row>;
			db->reserve(row);
			for (int j = 0; j < row; j++)
			{
				Row cur;
				cur.rowID = j;
				cur.lastRng = mat.M[j][i];
				cur.mat = new DPMatrix(cur.lastRng);
				db->push_back(cur);
			}
			vector<int> *dict = new vector<int>;
			if (k == 2)
			{
				for (int j = 0; j < col; j++)
				{
					if (j != i)
					{
						vector<int> seq1;
						seq1.push_back(i);
						seq1.push_back(j);
						dbStack.push(db);
						bottomCheck_ProbFreqApproxPoi(seq1, curRoot, nextRoot, dbStack, mat, minrow, mincol, th_prob, mu_m, fout);
						dbStack.pop();
					}
				}
			}
			else if (k == 3)
			{
				//构造字典
				vector<int> *dict = new vector<int>;
				for (vector<Node*>::iterator dit = (*it)->chList.begin(); dit != (*it)->chList.end(); dit++)
				{
					Node* child = *dit;
					dict->push_back(child->last);
				}
				//利用dict递归
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreqApproxPoi(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, dict, th_prob, mu_m, fout);
					dbStack.pop();
				}
				delete dict;
			}
			else
			{
				for (vector<Node*>::iterator it1 = (*it)->chList.begin(); it1 != (*it)->chList.end(); it1++)
				{
					int j = (*it1)->last;
					vector<int> seq1;
					seq1.push_back(i);
					seq1.push_back(j);
					dbStack.push(db);
					recursiveCheck_ProbFreqApproxPoi(k, 2, seq1, *it1, curRoot, nextRoot, dbStack, mat, minrow, mincol, NULL, th_prob, mu_m, fout);
					dbStack.pop();
				}
			}                                                                                            
			//释放db
			for (vector<Row>::iterator it1 = (*db).begin(); it1 != (*db).end(); it1++)
			{
				delete (*it1).mat;
			}
			delete db;
		}
		//更新当前tree
		delete curRoot;
		curRoot = nextRoot;
	}
	delete curRoot;
}
