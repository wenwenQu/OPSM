#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
using namespace std;

#include "datStruct.h"

void meanStdDev(const vector<double>& vec, double& mean, double& stddev)
{
	int s = vec.size();
	double sum{ 0. }, sqsum{ 0. };

	for (int i = 0; i < s; i++) {
		sum += vec[i];
		sqsum += vec[i] * vec[i];
	}

	mean = sum / s;
	double sqdev = sqsum / (s - 1) - mean * mean*s / (s - 1);
    assert(sqdev>0);
	stddev = sqrt(sqdev);
}


//输入文件格式：
//第一行rowNum, colNum, repeat count
//第二行开始，每行一个gene: x11 x12 x13 x14 x21 x22 x23 x24, ...

Matrix loadMatrix(string inFile){
	int row, col, cnt;
	ifstream fin(inFile.c_str());
	
	fin>>row>>col>>cnt;
	cout<<"num_rows = "<<row<<endl;
	cout<<"num_cols = "<<col<<endl;

	Range** mat = new Range*[row];
	for(int i=0; i<row; i++){
		mat[i]=new Range[col/cnt];
		for(int j=0; j<col; j=j+cnt){
			double u,sigma;
			vector<double> x(cnt);
			vector<double> x2(7);
			vector<double> y(7);
			for(int k = 0; k<cnt; k++){
				fin>> x[k];
			}
			bool same = true;
			for (int k = 0; k < cnt-1; k++) {
				if (x[k + 1] != x[k]) {
					same = false;
					break;
				}
			}
			if (same) {
				u = x[0];
				sigma = 0.01;
			}
			else {
				meanStdDev(x, u, sigma);
			}
				for (int k = 0; k < total_std; k++) {
					x2[k] = u + (k - std_num) * sigma;
					y[k] = (1 / (sigma*sqrt(2 * M_PI)))*(exp(-((x2[k] - u)*(x2[k] - u)) / (2 * sigma*sigma)));
				}
				mat[i][j / cnt].u = u;
				mat[i][j / cnt].sigma = sigma;
				mat[i][j / cnt].cs = spline(x2, y);
				mat[i][j / cnt].l = u - std_num * sigma;
				mat[i][j / cnt].r = u + std_num * sigma;
		}
	}
	fin.close();
	cout << "Finish reading input matrix"<< endl;
	return Matrix(mat, row, col/cnt);
}

//第一行seq,expSup
//第二行"rowID,
void outputExpSup(vector<int>& seq, vector<Row>& db, double expSup, ofstream& fout)
{
	//以下输出部分可以灵活改动
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		fout<<*it<<" ";
	}
	fout<<","<<expSup;
	fout<<endl;

	//sort(db.begin(), db.end(), rowComp);//############ 排序后输出 #############
	for(vector<Row>::iterator it=db.begin(); it!=db.end(); it++)
	{
		fout<<it->rowID<<","<<it->mat->getPr()<<" ";
	}
	fout<<endl;
}

bool rowComp(Row& row1, Row& row2)
{
	return row1.mat->getPr() > row2.mat->getPr();//概率从大到小
}

void outputFreq(vector<int>& seq, vector<Row>& db, double expSup, ofstream& fout)
{
	//以下输出部分可以灵活改动
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		fout<<*it<<" ";
	}
	fout<<","<<expSup;
	fout<<endl;

	//sort(db.begin(), db.end(), rowComp);//############ 排序后输出 #############
	for(vector<Row>::iterator it=db.begin(); it!=db.end(); it++)
	{
		fout<<it->rowID<<","<<it->mat->getPr()<<" ";
	}
	fout<<endl;
}

/*
void outputExpSup(vector<int>& seq, vector<Row>& db, double expSup, ofstream& fout)
{
	//以下输出部分可以灵活改动
	fout<<"seq: ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		fout<<*it<<" ";
	}
	fout<<" | expSup = "<<expSup;
	fout<<endl;

	for(vector<Row>::iterator it=db.begin(); it!=db.end(); it++)
	{
		fout<<"("<<it->rowID<<", "<<it->mat->getPr()<<") ";
	}
	fout<<endl;
}

bool rowComp(Row& row1, Row& row2)
{
	return row1.mat->getPr() > row2.mat->getPr();//概率从大到小
}

void outputFreq(vector<int>& seq, vector<Row>& db, double th_prob, ofstream& fout)
{
	//以下输出部分可以灵活改动
	fout<<"seq: ";
	for(vector<int>::iterator it=seq.begin(); it!=seq.end(); it++)
	{
		fout<<*it<<" ";
	}
	fout<<endl;

	sort(db.begin(), db.end(), rowComp);//############ 排序后输出 #############
	for(vector<Row>::iterator it=db.begin(); it!=db.end(); it++)
	{
		fout<<"("<<it->rowID<<", "<<it->mat->getPr()<<") ";
	}
	fout<<endl;
}
*/
