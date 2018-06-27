#ifndef _GETMAX_H_
#define _GETMAX_H_

#include <vector>

using namespace std;

template <class T> int getMax(vector<T> &vec, T &max, int &maxPos);

template <class T>
int getMax(vector<T> &vec, T &max, int &maxPos)
{
	int len = vec.size();
	if (len <= 0)
	{
		return 0;
	}

	int i = 0;
	typename vector<T>::const_iterator cit;
	cit = vec.begin();
	max = *cit++;
	maxPos = 0;

	for (; cit!=vec.end(); cit++)
	{
		i ++;
		if (*cit > max)
		{
			max = *cit;
			maxPos = i;
		}
	}
	return 1;
}

template <class T>
int getMatMax(vector<vector<T> > &mat, T &max, int &maxOuterPos, int &maxPos)
{
	int outerSize = mat.size();
	if (outerSize <= 0)
		return 0;
	int innerSize = mat[0].size();
	if (innerSize <= 0)
		return 0;

	int i, j;
	maxOuterPos = 0;
	maxPos = 0;
	max = mat[0][0];

	for (i=0; i<outerSize; i++)
	{
		for (j=0; j<innerSize; j++)
		{
			if (max < mat[i][j])
			{
				max = mat[i][j];
				maxPos = j;
				maxOuterPos = i;
			}
		}
	}
	return 1;
}

template <class T>
int sortIndex(vector<T> &vec, vector<int> &sortedIdx)
{
	int len = vec.size();
	int i, j, tmpIdx;
	T tmp;
	if (len <= 0)
	{
		return 0;
	}
	sortedIdx.resize((unsigned int)len);
	for (i=0; i<len; i++)
		sortedIdx[i] = i;

	for (i=0; i<len; i++)
		for (j=(i-1); j>=0; j--)
			if (vec[j+1] > vec[j])
			{
				tmp = vec[j+1];
				vec[j+1] = vec[j];
				vec[j] = tmp;

				tmpIdx = sortedIdx[j+1];
				sortedIdx[j+1] = sortedIdx[j];
				sortedIdx[j] = tmpIdx;
			}
	return 1;
}

template <class T>
int getTopN(vector<T> &vec, vector<T> &controlVec, vector<T> &maxVec, vector<int> &maxPosVec, int n, int dist)
{
	int i, j, k;
	int len = vec.size();
	int controlLen = controlVec.size();
	if (len != controlLen || len <= 0 || n < 1)
		return 0;
	
	T threshold;
	vector<int> sortedIdx;
	vector<T> sortedVal;
	sortedIdx.resize(len);
	sortedVal.resize(len);

	maxVec.clear();
	maxPosVec.clear();

	for (i=0; i<len; i++)
		sortedVal[i] = vec[i];
	k = 1;

	sortIndex(sortedVal, sortedIdx);

	maxVec.push_back(sortedVal[0]);
	maxPosVec.push_back(sortedIdx[0]);
	threshold = controlVec[sortedIdx[0]] / 3;
	bool submax;
	for (i=1;i<len; i++)
	{
		submax = true;
		for (j=0; j<i; j++)
		{
			if (sortedIdx[i]+1 == sortedIdx[j] || sortedIdx[i]-1 == sortedIdx[j])
			{
				submax = false;
				break;
			}
		}
		if (submax && sortedIdx[i]!=0 && sortedIdx[i]!=(len-1))
		{
			if (controlVec[sortedIdx[i]] > threshold && (sortedIdx[i] < (sortedIdx[0] - dist) || sortedIdx[i] > (sortedIdx[0] + dist)))
			{
				maxVec.push_back(sortedVal[i]);
				maxPosVec.push_back(sortedIdx[i]);
				k ++;
			}
			else
				;
			if (k == n)
				return 1;
		}
	}
	return 1;
}

template <class T>
int multiVec(vector<T> &vec1, vector<T> &vec2, vector<T> &vecRST)
{
	int n, m, i;
	n = vec1.size();
	m = vec2.size();
	if (n != m || n <= 0)
		return 0;

	vecRST.resize(n);
	for (i=0; i<n; i++)
		vecRST[i] = vec1[i] * vec2[i];
	return 1;
}

template <class T>
int logVec(vector<T> &vec1, vector<T> &vecRST)
{
  int n, i;
  n = vec1.size();
  if (n <= 0)
    return 0;

  vecRST.resize(n);
  for (i=0; i<n; i++)
    vecRST[i] = log(vec1[i]);
  return 1;
}

template <class T>
int get2BSTop(vector<vector<T> > &mat, vector<T> &max, vector<int> &maxPos, vector<int> &distBSIdx)
{
	if (mat.size() <= 0)
		return 0;

	max.clear();
	maxPos.clear();
	distBSIdx.clear();
	int i, maxPos0;
	T max0;
	if (!getMatMax(mat, max0, i, maxPos0))
		return 0;
	max.push_back(max0);
	maxPos.push_back(maxPos0);
	distBSIdx.push_back(i);
	return 1;
}

template <class T>
int getThickRegion(vector<T> &vec, T &max, int &summit, int &thickS, int &thickE, T &cutValue) 
{
  if(vec.size() <= 0 || summit > (int)vec.size() || cutValue > max)
  {
    //fprintf(stderr, "\tmax: %.2f\tcut: %.2f\n", max, cutValue); 
    return 0;
  }
  thickS = 0;
  thickE = (int) vec.size() - 1;
  int i;
  for (i=summit; i>=0; i--)
  {
    if (vec[i] < cutValue)
    {
      thickS = i + 1;
      break;
    }
  }
  for (i=summit; i<(int)vec.size(); i++)
  {
    if (vec[i] < cutValue)
    {
      thickE = i - 1;
      break;
    }
  }
  return 1;
}

#endif //_GETMAX_H_
