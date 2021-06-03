#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <assert.h>
using namespace std;

double MISSING_VAL = -999999;

struct Range {
	double l, r;//如果是missing value, 则规定l=r=MISSING_VAL

	Range() {}

	Range(double l, double r)
	{
		this->l = l;
		this->r = r;
	}

	void set(double l, double r)
	{
		this->l = l;
		this->r = r;
	}

	bool isMissing()
	{
		return (l == MISSING_VAL);
	}

	string toString()
	{
		if (isMissing()) return "?";
		stringstream ss;
		ss << "[" << l << ", " << r << "]";
		return ss.str();
	}

	bool contains(Range sub)
	{
		return l <= sub.l && r >= sub.r;
	}

	double gap()
	{
		return (r - l);
	}

	double pdf()
	{
		if (l == r) return 0;
		return 1 / (r - l);
	}
};

struct Matrix
{
	Range** M;
	int row;
	int col;

	Matrix(Range** M, int row, int col)
	{
		this->M = M;
		this->row = row;
		this->col = col;
	}

	~Matrix()
	{
		for (int i = 0; i < row; i++) delete[] M[i];
		delete[] M;
	}

	void report()
	{
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				cout << M[i][j].toString() << "\t";
			}
			cout << endl;
		}
	}
};

struct DPMatrix
{
	int seqLen;//order的长度
	Range* seq;//order sequence
	int rngNum;//subrange数量
	double* splits;//subrange分隔点
	double** DP;//DP矩阵

	DPMatrix()
	{
		//不能用默认的，不然析构函数会报错
		seqLen = 0;
		seq = NULL;
		rngNum = 0;
		splits = NULL;
		DP = NULL;
	}

	DPMatrix(Range& rng)//DP对应于x=0的base case
	{
		//计算第一行时调用
		seqLen = 1;
		rngNum = 1;

		seq = new Range[1];
		seq[0] = rng;

		splits = new double[2];
		splits[0] = rng.l;
		splits[1] = rng.r;

		DP = new double*[1];
		DP[0] = new double[1];
		DP[0][0] = 1;
	}

	double getPr()
	{
		return DP[seqLen - 1][rngNum - 1];
	}

	void report()
	{
		for (int i = 0; i < seqLen; i++)
		{
			for (int j = 0; j < rngNum; j++)
			{
				cout << DP[i][j] << "\t";
			}
			cout << endl;
		}
	}

	~DPMatrix()
	{
		for (int i = 0; i < seqLen; i++) delete[] DP[i];
		if (DP != NULL) delete[] DP;
		if (seq != NULL) delete[] seq;
		if (splits != NULL) delete[] splits;
		//cout<<seqLen<<" x "<<rngNum<<" freed !"<<endl;//###################
	}
};

struct Row
{//矩阵的一行，维护算法需要的辅助信息
	int rowID;//行数
	Range lastRng;//用于词典pruning
	DPMatrix* mat; //概率计算用的DP矩阵, 为了防止出了block后就释放，需要动态分配+手动释放
};

//================= 以下结构用于prefixTree =================
struct Node
{
	int last;//序列的最后一个元素，root为-1 (level 0)
	vector<Node*> chList;//子节点

	Node(int last)
	{
		this->last = last;
	}

	~Node()
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			delete *it;
		}
	}

	void appendChild(Node* child)
	{
		chList.push_back(child);
	}

	Node* goToChild(int last)
	{
		for (vector<Node*>::iterator it = chList.begin(); it != chList.end(); it++)
		{
			if ((*it)->last == last) return *it;
		}
		return NULL;
	}
};

Node* initTree()
{
	//return root
	return new Node(-1);
}

bool emptyTree(Node* root)
{
	return (root->chList.size()) == 0;
}

bool exists(Node* root, vector<int> seq)
{
	Node* cur = root;
	for (vector<int>::iterator it = seq.begin(); it != seq.end(); it++)
	{
		int last = *it;
		Node* next = cur->goToChild(last);
		if (next == NULL) return false;
		cur = next;
	}
	return true;
}

bool subseqExists(Node* root, vector<int> seq)
{//anti-monoticity pruning
	int size = seq.size();
	for (int i = 0; i < size; i++)
	{
		vector<int> subseq;
		for (int j = 0; j < size; j++)
		{
			if (j != i) subseq.push_back(seq[j]);
		}
		if (!exists(root, subseq)) return false;
	}
	return true;
}

//前面的函数用于curRoot
//下面的函数用于nextRoot
void addSeq(Node* root, vector<int> seq)
{
	Node* cur = root;
	for (vector<int>::iterator it = seq.begin(); it != seq.end(); it++)
	{
		int last = *it;
		Node* next = cur->goToChild(last);
		if (next == NULL)
		{
			next = new Node(last);
			cur->appendChild(next);
		}
		cur = next;
	}
}

//以下用于调试
void printTree(int level, Node* cur)
{
	for (int i = 0; i < level; i++) cout << "=== ";
	cout << cur->last << endl;
	vector<Node*> list = cur->chList;
	for (vector<Node*>::iterator it = list.begin(); it != list.end(); it++)
	{
		printTree(level + 1, *it);
	}
}

void printTree(Node* root)
{
	vector<Node*> list = root->chList;
	for (vector<Node*>::iterator it = list.begin(); it != list.end(); it++)
	{
		printTree(1, *it);
	}
}
