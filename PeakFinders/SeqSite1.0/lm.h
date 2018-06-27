#ifndef LM_H
#define LM_H

#include <vector>
#include <cmath>
#include "seqsite.h"
#include "getMax.h"
#include "global.h"
#include "parameter.h"
extern parameter param;

using namespace std;

class lm
{
public:
	lm(void);  // linear model without intercept: y = slope * x + epsilon

	void solve(vector<double> y);
	bool slipSolve(vector<double> y, char strand);
	void generateX();
	bool slipSolve2BS(vector<double> iny, char strand);

	int lm_start;	// lm start of x
	int lm_end;		// lm end of x

	double result[2];
	vector<double> slopeVec;
	vector<double> R2Vec;
	vector<vector<double> > slopeMat;
	vector<vector<double> > R2Mat;
	vector<int> distVec;

protected:
	vector<double> min_y;
	vector<double> x;
	vector<double> xx;
  double min_slope;
	double x_max;
	int x_max_pos;
	double R2; // R squared
	double slope;


public:
	~lm(void);
};

#endif
