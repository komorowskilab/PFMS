#include "lm.h"

#define DEBUG_LM

lm::lm(void)
{
}

lm::~lm(void)
{
}

void lm::solve(vector<double> y)
{
	double Sxy = 0;
	double Sxx = 0;
	double Syy = 0;
	int i;
	int n = x.size();
	int m = y.size();
	if (n != m)
	{
		printf("Please make sure that the two vectors are of the same length!\n");
		exit(1);
	}

	for (i=0; i<n; i++)
	{
		Sxx += x[i] * x[i];
		Syy += y[i] * y[i];
		Sxy += x[i] * y[i];
	}

	slope = Sxy / Sxx;
	R2 = slope * Sxy / Syy;

	result[0] = R2;
	result[1] = slope;
}

bool lm::slipSolve(vector<double> iny, char strand)  // strand: +(f) / -(r)
{
	//int n = x.size();
	double Sxy = 0;
	double Sxx = 0;
	double Syy = 0;
	double tmp1, tmp2;
	int i, j, m;
	vector<double> y;
	int lm_arm = param.lmArm;
	for (i=0; i<lm_arm; i++)
		y.push_back(0);
	for (i=0; i<(int)iny.size(); i++)
		y.push_back(iny[i]);
	for (i=0; i<lm_arm; i++)
		y.push_back(0);
	m = y.size();

	lm_start = 0;						// start of x
	lm_end = lm_start + param.lmLength;	// end of x
	if ( lm_end > m)
	{
		printf("length(x)(%d) > length(y)(%d)\n", lm_end, m);
		return false;
	}
	for (i=lm_start; i<lm_end; i++)
	{
		Sxx += x[i] * x[i];
	}
	
	slopeVec.clear();
	R2Vec.clear();

	if (strand == '+' || strand == 'f')
	{
		for (j=(m-1); j>=(lm_end-1); j--)
		{
			Syy = 0;
			Sxy = 0;
			for (i=lm_start; i<lm_end; i++)
			{
				Syy += y[j-i] * y[j-i];
				Sxy += x[i] * y[j-i];
			}
			tmp1 = Sxy / Sxx;
			tmp2 = tmp1 * Sxy / Syy;
			slopeVec.push_back(tmp1);
			R2Vec.push_back(tmp2);
		}
	}
	else if (strand == '-' || strand == 'r')
	{
		for (j=0; j<=(m-lm_end); j++)
		{
			Syy = 0;
			Sxy = 0;
			for (i=lm_start; i<lm_end; i++)
			{
				Syy += y[j+i] * y[j+i];
				Sxy += x[i] * y[j+i];
			}
			tmp1 = Sxy / Sxx;
			tmp2 = tmp1 * Sxy / Syy;
			slopeVec.push_back(tmp1);
			R2Vec.push_back(tmp2);
		}
	}
	else
		return false;

#ifdef DEBUG_LM1
	for (i=0; i<=(m - lm_end); i++)
	{
		printf("Pos:%d\tslope: %f\tR2: %f\n", i, slopeVec[i], R2Vec[i]);
	}
#endif

  for(i=0; i<(int)slopeVec.size(); i++)
    slopeVec[i] = slopeVec[i] / min_slope;
	return true;
}

void lm::generateX()
{
	double lambda = param.fragLengthMean / param.fragLengthVar;
	double alpha = param.fragLengthMean * lambda;
	double half_alpha = alpha / 2.0;
	int i;
	int vectorLength = (int) ( 2 * param.fragLengthMean + 100);
	double tmp;
	x.push_back(0);
	for (i=1; i<vectorLength; i++)
	{
		tmp = pow(i, (half_alpha - 1)) / exp(i * lambda) / 100000;
		x.push_back(tmp);
	}
#ifdef DEBUG_LM1
	for (i=0; i<vectorLength; i++)
	{
		printf("%f\t", x[i]);
	}
	printf("\n");
#endif
	getMax(x, x_max, x_max_pos);
  min_y.resize(x.size());
  for (i=0; i<(int)min_y.size(); i++) 
    min_y[i] = 0;
  min_y[x_max_pos] = 1.0;

  solve(min_y);
  min_slope = slope;

  //printVec(x);
  //printVec(min_y);
  //printf("pos: %d, max_x: %f, min_R2: %f, min_slope: %f\n", x_max_pos, x_max,R2, min_slope);

}

bool lm::slipSolve2BS(vector<double> iny, char strand)  // strand: +(f) / -(r)
{
	//int n = x.size();
	int dist, maxDist = 100;
	double Sxy = 0;
	double Sxx = 0;
	double Syy = 0;
	double tmp1, tmp2;
	double sumTmp;
	int i, j, m;
	vector<double> y;
	int lm_arm = param.lmArm;
	for (i=0; i<lm_arm; i++)
		y.push_back(0);
	for (i=0; i<(int)iny.size(); i++)
		y.push_back(iny[i]);
	for (i=0; i<lm_arm; i++)
		y.push_back(0);
	m = y.size();

	lm_start = 0;						// start of x
	lm_end = lm_start + param.lmLength;	// end of x
	if ( lm_end > m)
	{
		printf("length(x)(%d) > length(y)(%d)\n", lm_end, m);
		return false;
	}

	xx.resize(lm_end);
	distVec.clear();
	slopeMat.clear();
	R2Mat.clear();
	for (dist=0; dist<maxDist; dist+=10)
	{
		// count new x
		sumTmp = 0;
		for (i=lm_start; i<lm_end; i++)
		{
			if (i < dist)
				xx[i] = x[i];
			else
				xx[i] = x[i] + x[i-dist];
			sumTmp += xx[i];
		}

		// normalize new x
		for (i=0; i<lm_end; i++)
			xx[i] /= sumTmp;

		// begin regression
		slopeVec.clear();
		R2Vec.clear();

		Sxx = 0;
		for (i=lm_start; i<lm_end; i++)
			Sxx += xx[i] * xx[i];

		if (strand == '+' || strand == 'f')
		{
			for (j=(m-1); j>=(lm_end-1); j--)
			{
				Syy = 0;
				Sxy = 0;
				for (i=lm_start; i<lm_end; i++)
				{
					Syy += y[j-i] * y[j-i];
					Sxy += xx[i] * y[j-i];
				}
				tmp1 = Sxy / Sxx;
				tmp2 = tmp1 * Sxy / Syy;
				slopeVec.push_back(tmp1);
				R2Vec.push_back(tmp2);
			}
		}
		else if (strand == '-' || strand == 'r')
		{
			for (j=0; j<=(m-lm_end); j++)
			{
				Syy = 0;
				Sxy = 0;
				for (i=lm_start; i<lm_end; i++)
				{
					Syy += y[j+i] * y[j+i];
					Sxy += xx[i] * y[j+i];
				}
				tmp1 = Sxy / Sxx;
				tmp2 = tmp1 * Sxy / Syy;
				slopeVec.push_back(tmp1);
				R2Vec.push_back(tmp2);
			}
		}
		else
			return false;
		
		distVec.push_back(dist);
		slopeMat.push_back(slopeVec);
		R2Mat.push_back(R2Vec);
	}

#ifdef DEBUG_LM1
	for (i=0; i<=(m - lm_end); i++)
	{
		printf("Pos:%d\tslope: %f\tR2: %f\n", i, slopeVec[i], R2Vec[i]);
	}
#endif
	return true;
}
