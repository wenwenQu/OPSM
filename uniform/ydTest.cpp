#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib> 
#include <string>
using namespace std;

#include "io.h"
#include "algo.h"
#include "systemI.h"
extern std::atomic<bool> global_end_label;

extern double tt1, tt2, tt3; //@@@@@@@@@@@@@@
int main(int argc, char** argv)
{


	//argv[1] input file name
	//argv[2] method name
	//argv[3] min row
	//argv[4] min col
	//argv[5] optional th prob

	
	cout << "*************************************************************************" << endl;
	cout << "Please refer to the paper/read_me.txt for method names and more info." << endl;
	cout << "Method: Apri_ES  DFS_ES" << endl;
	cout << "opsm.exe input_file_name method_name min_row min_col " << endl << endl;
	cout << "Method: Apri_PF  DFS_PF  Apri_PFA  DFS_PFA   Apri_PFP  DFS_PFP" << endl;
	cout << "opsm.exe input_file_name method_name min_row min_col th_prob" << endl;
	cout << "*************************************************************************" << endl;
	

	cout << argv[0] << endl << argv[1] << endl;
	clock_t s1, e1;
	s1 = clock();
	Matrix mat = loadMatrix(argv[1]);
	e1 = clock();
	cout << "loading time: " << (double)(e1 - s1) / CLOCKS_PER_SEC << " s" << endl;


	string str = argv[2];
	int minrow = atoi(argv[3]);
	int mincol = atoi(argv[4]);
	double th_prob;
	int th_prob_fileName = th_prob * 100;
	cout << "minrow= " << minrow << "   mincol=" << mincol << "   thprob=" << endl;

	char outputfileTime[1000];
	char outputfile[1000];
	char outputfilePeakMem[1000];

	if (str == "Apri_ES" || str == "DFS_ES") {
		sprintf(outputfile, "opsm_output_%s_r%d_c%d.txt", argv[2], minrow, mincol);
		sprintf(outputfileTime, "runtime_opsm_output_%s_r%d_c%d.txt", argv[2], minrow, mincol);
		sprintf(outputfilePeakMem, "maxmem_opsm_output_%s_r%d_c%d.txt", argv[2], minrow, mincol);
	}
	else {
		th_prob = atof(argv[5]);
		sprintf(outputfilePeakMem, "maxmem_opsm_output_%s_r%d_c%d.txt", argv[2], minrow, mincol);
		sprintf(outputfile, "opsm_output_%s_r%d_c%d_p%d.txt", argv[2], minrow, mincol, th_prob_fileName);
		sprintf(outputfileTime, "runtime_opsm_output_%s_r%d_c%d_p%d.txt", argv[2], minrow, mincol, th_prob_fileName);
	}

	ofstream fout(outputfile);
	ofstream foutPeakMem(outputfilePeakMem);
	ofstream foutTime(outputfileTime);


	//GetCurrentPid();
	//thread t = thread(info, GetCurrentPid(), ref(foutPeakMem)); //lauch a thread to record memory

	clock_t start, end;
	start = clock();
	cout << "start Running" << endl;

	if (str == "Apri_ES") 
	{
		Apriori_ExpSup(mat, minrow, mincol, fout);
	}

	else if (str == "DFS_ES")
	{
		DFS_ExpSup(mat, minrow, mincol, fout);
	}
	else if (str == "Apri_PF")
	{
		Apriori_ProbFreq(mat, minrow, mincol, th_prob, fout);
	}

	else if (str == "DFS_PF")
	{
		DFS_ProbFreq(mat, minrow, mincol, th_prob, fout);
	}

	else if (str == "Apri_PFA")
	{
		Apriori_ProbFreqApprox(mat, minrow, mincol, th_prob, fout);
	}

	else if (str == "DFS_PFA")
	{
		DFS_ProbFreqApprox (mat, minrow, mincol, th_prob, fout);
	}

	//PFA Poisson 
	else if (str == "DFS_PFP")
	{
		DFS_ProbFreqApproxPoi(mat, minrow, mincol, th_prob, fout);
	}
	else if (str == "Apri_PFP")
	{
		Apriori_ProbFreqApproxPoi(mat, minrow, mincol, th_prob, fout);
	}
	
	end = clock();
	cout << "Running time: " << (double)(end - start) / CLOCKS_PER_SEC << " s" << endl;
	


	foutTime << (double)(end - start) / CLOCKS_PER_SEC << " s" << endl;//@@@@@@@@@@@@@@@@@
	foutTime.close();
	
	//global_end_label = false;
	//t.join();
	foutPeakMem.close();
	
	fout.close();

	return 0;
}