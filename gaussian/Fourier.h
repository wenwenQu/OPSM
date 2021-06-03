//使用 http://www.codeproject.com/Articles/9388/How-to-implement-the-FFT-algorithm
#pragma once

class CFourier
{
public:
	double pi;
	double *vector;
	CFourier(void);
	~CFourier(void);
	// FFT 1D
	//number_of_samples 必须是2^n
	void FFT(double data[], unsigned long number_of_samples);//data是re,re,re...
	void IFFT(double vec[], unsigned long number_of_samples);//vec是re,im,re,im...
	//结果都在vector中
};
