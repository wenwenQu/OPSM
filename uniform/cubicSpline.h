/*
 * cubicSpline.h
 *
 *  Created on: 2021Äê4ÔÂ16ÈÕ
 *      Author: vinny
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_


#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include "Polynomial.h"
using namespace std;

typedef vector<double> vec;

typedef Polynomial<double> SplineSet;

vector<SplineSet> spline(vec &x, vec &y)
{
    int n = x.size()-1;
	vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

    vec alpha;
    alpha.push_back(0);
    for(int i = 1; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

    vec c(n+1);
    vec l(n+1);
    vec mu(n+1);
    vec z(n+1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for(int j = n-1; j >= 0; --j)
    {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
    }


	vector<SplineSet> output_set(n);
    for(int i = 0; i < n; ++i)
    {
    	double x2 = x[i]*x[i];
    	double x3 = x2*x[i];
    	vector<double>& coef = output_set[i].data;
		coef.resize(4);
    	coef[0] = a[i]-b[i]*x[i]+c[i]*x2-d[i]*x3;
    	coef[1] = b[i]-2*c[i]*x[i]+3*d[i]*x2;
    	coef[2] = c[i]-3*d[i]*x[i];
    	coef[3] = d[i];
    	output_set[i].normalize(coef);
    }
    return output_set;
}

void getcuber(vec& x, vec& y, vector<SplineSet>&/*out*/ cs)
{
    for(int i = 0; i < x.size(); ++i)
    {
        x[i] = i;
        y[i] = sin(i);
    }

    cs = spline(x, y);
    //for(int i = 0; i < cs.size(); ++i)
    //    cout << cs[i].d << "\t" << cs[i].c << "\t" << cs[i].b << "\t" << cs[i].a << endl;
}



#endif /* CUBICSPLINE_H_ */
