#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

FILE *mustOpen(char *fileName, char *mode);
char* readLine(FILE* fh);
int bedColNum(char* line);
void chopLine(char* line, int colNum, char** wordArray);
bool isInt(char* str);
bool isChr(char* str);
bool isStrand(char *str);

template <class T> int printVec(vector<T> &vec) // 1: success; 0: failure
{
	int len = vec.size();
	if (len <= 0)
	{
		cout<<"empty vector"<<endl;
		return 0;
	}
	int i;
	for (i=0; i<len; i++)
	{
		cout<<vec[i]<<"\t";
	}
	cout<<endl;
	return 1;
}

template <class T> void clearVector(vector<T> &vt)
{
	vector<T> vtTemp;
	vtTemp.swap(vt);
}

template <class T> T sum(vector<T> &vt)
{
	int vecSize = vt.size();
	int i;
	T vecSum = 0;
	for (i=0; i<vecSize; i++)
		vecSum += vt[i];
	return vecSum;
}

template <class T> double weightAvg(vector<T> &vt)
{
	int vecSize = vt.size();
	int i;
	T vecWSum = 0;
	T vecSum = 0;
	for (i=0; i<vecSize; i++)
	{
		vecWSum += vt[i] * i;
		vecSum += vt[i];
	}
	return (double) vecWSum / (double) vecSum ;
}

template <class T> double weightAvg(vector<T> &vt, int thk_s, int thk_e)
{
	int vecSize = vt.size();
  if( thk_s < 0 || thk_e >= vecSize) 
  {
    fprintf(stderr, "ERROR: start or end exceeds the vector boundary.\n");
    exit (0);
  }
	int i;
	T vecWSum = 0;
	T vecSum = 0;
	for (i=thk_s; i<=thk_e; i++)
	{
		vecWSum += vt[i] * i;
		vecSum += vt[i];
	}
	return (double) vecWSum / (double) vecSum ;
}

template <class T> double avg(vector<T> &vt)
{
	int vecSize = vt.size();
	int i;
	T vecSum = 0;
	for (i=0; i<vecSize; i++)
		vecSum += vt[i] ;
	return (double) vecSum / (double) vecSize ;
}

template <class T> double var(vector<T> &vt)
{
	int vecSize = vt.size();
	int i;
	T vecSum = 0;
	T vecSum2 = 0;
	double S2;
	for (i=0; i<vecSize; i++)
	{
		vecSum += vt[i];
		vecSum2 += vt[i] * vt[i];
	}
	S2 = vecSum2 - vecSum * vecSum / (double) vecSize;
	return S2 / (double) (vecSize - 1);
}

unsigned int minCntGivenFDR(int readCnt, double genomeSize, double FDR, int regionLength);

// count cdf for Poisson 
double upper_pois(int n, double lambda); // the cdf not including n
double lower_pois(int n, double lambda); // the cdf including n
unsigned int pois_cdf_inv(double lambda, double cdf); // inverse poisson distribution

// for Fisher's exact test
double logSum(int lower, int upper);
double fisherTest(int aa, int bb, int n, int m);

#endif /* GLOBAL_H */
