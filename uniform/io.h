#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include "datStruct.h"

//输入文件格式：
//第一行rowNum, colNum
//第二行开始，每行一个gene: l1 r1 l2 r2 ...
//表示ranges: [l1, r1], [l2 ,r2], ...
Matrix loadMatrix(string inFile){
	int row, col;
	ifstream fin(inFile.c_str());
	
	fin>>row>>col;
	cout<<"num_rows = "<<row<<endl;
	cout<<"num_cols = "<<col<<endl;

	Range** mat = new Range*[row];
	for(int i=0; i<row; i++){
		mat[i]=new Range[col];
		for(int j=0; j<col; j++){
			double l, r;
			fin>>l>>r;
			/*
			if (l == r) {
				cerr << "The interval cannot be a point, the endpoints of the interval cannot have the same value" << endl;
				cerr << "Please check inverval [" << l << ", " << r << "]" << "in row = " << i << "\tcol = " << j << " in the input file" << endl;
				//cerr << ""
			}				
			assert(l != r);
			*/
			mat[i][j].set(l, r);
		}
	}
	fin.close();
	cout << "Finish reading input matrix"<< endl;
	return Matrix(mat, row, col);
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