#ifndef _READS_H
#define _READS_H

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include "lm.h"
#include "getMax.h"
#include "global.h"
#include "seqsite.h"
#include "parameter.h"
extern parameter param;

using namespace std;

class reads
{
public:
	reads(void);

	int multiBSconut;
	int singleBScount;
	int bCountArray[6];
	//int ambiguousBScount;

	// the number of reads
	int readsCnt;
	// # forward reads
	int fwdReadsCnt;
	// # reverse reads
	int rvsReadsCnt;
	// # (discarded reads)
	int discarderCnt;
	// chromosome name
	vector<string> chrName;
	// forward read position
	vector<vector<int> > fpos;
	// reverse read position
	vector<vector<int> > rpos;
	// chromosome map
	map< string, int > chrIdxMap;
	// the number of chromosome
	int totalChr;

	// fwd background lambda
	double lambda_bg_f;
	// rvs background lambda
	double lambda_bg_r;

	vector<vector<vector<int> > > fisland;
	vector<vector<vector<int> > > risland;

private:
	// map iterator
	map< string, int >::iterator map_iter;

	// the maximum slope for this data, corresponding to 100 in output.
	double maxSlope;
	double maxR2;
	double maxFC;

public:
	class bingdingSite {
	public:
		bingdingSite(int chrIdx, int pos, double R2, double slope, char strand, double pval, double foldChange, double qval1, double qval2):
		  chrIdx(chrIdx), pos(pos), R2(R2), slope(slope), strand(strand), pval(pval), foldChange(foldChange),qval1(qval1), qval2(qval2){}
		//bingdingSite(int pos, double R2, char strand):pos(pos), R2(R2), strand(strand){}
		int chrIdx;
		int pos;
		double R2;
		double slope;
		char strand;
		double pval;
		double foldChange;
		double qval1;
		double qval2;
		int type; // type: 0 for one sided, 1 for two sided but not merged, 2 for merged binding sites.
		bool operator < (const bingdingSite &bs ) const 
		{
			return chrIdx < bs.chrIdx || (chrIdx == bs.chrIdx && pos < bs.pos);
		}
	};
	vector<bingdingSite> BSVec;
	vector<vector<bingdingSite> > BSClusterVec;

	class readsCluster {	// a read cluster 
	public:
		readsCluster(vector<int> read_vec, int chrIdx, char strand){
			r_count = read_vec.size();
			r_start = read_vec[0];
			r_end = read_vec[r_count - 1];
			r_length = r_end - r_start + 1;
			r_count = read_vec.size();
			r_strand = strand;
			r_chrIdx = chrIdx;
			r_dense = (double) r_count / (double) (r_length + 2*param.minDist);

			// compute the signal curve + smooth the signal curve
			int j, k, pos;
      if (param.smoothArm > 0)  // old flat window smoothing
      {
			  r_lengthEx = r_length + 2 * param.smoothArm;
			  r_vec.resize(r_lengthEx);
      }
      else // for new QuEST-like Gaussian smoothing
      {
        r_lengthEx = r_length + 2 * 3 * param.smoothBandwidth;
        r_vec.resize(r_lengthEx);
      }

			for (j=0; j<r_lengthEx; j++)
			{
				r_vec[j] = 0.0;
			}

      if (param.smoothArm > 0) // old flat window smoothing
      {
			  for (j=0; j<r_count; j++)
			  {
				  pos = read_vec[j];
				  for (k=(pos-param.smoothArm); k<=(pos+param.smoothArm); k++)
				  {
					  r_vec[ k - r_start + param.smoothArm ] ++;
				  }
			  }
        for (j=0; j<r_lengthEx; j++)
        {
          r_vec[j] = r_vec[j] / (2.0 * param.smoothArm + 1.0);
        }
      }
      else  // for new QuEST-like Gaussian smoothing
      {
        for (j=0; j<r_count; j++)
        {
          pos = read_vec[j];
          for (k=(pos - 3*param.smoothBandwidth); k<=(pos + 3*param.smoothBandwidth); k++)
          {
            r_vec[ k - r_start + 3*param.smoothBandwidth] += exp ( - pow((double)(pos-k)/(double)param.smoothBandwidth, 2) / 2);
          }
        }
        for (j=0; j<(int)r_vec.size(); j++)
        {
          r_vec[j] /= param.smoothBandwidth * sqrt (2 * 3.1415926);
        }
      }

      if( ! getMax(r_vec, r_summitVal, r_summit)) fprintf(stderr, "ERROR: failed in getting tag cluster summit!\n");
      double cutRatio = 0.1;
      double cutValue = cutRatio * r_summitVal;
      if( ! getThickRegion(r_vec, r_summitVal, r_summit, r_thickStart, r_thickEnd, cutValue)) fprintf(stderr, "ERROR: failed in getting thick region!\n");

      if (param.smoothArm > 0)  // old flat window smoothing
      {
			  r_avgPos = weightAvg(r_vec, r_thickStart, r_thickEnd) + r_start - param.smoothArm;
        r_summit += r_start - param.smoothArm;
        r_thickStart += r_start - param.smoothArm;
        r_thickEnd += r_start - param.smoothArm;
      }
      else
      {
        r_avgPos = weightAvg(r_vec, r_thickStart, r_thickEnd) + r_start - 3 * param.smoothBandwidth;
        r_summit += r_start - 3 * param.smoothBandwidth;
        r_thickStart += r_start - 3 * param.smoothBandwidth;
        r_thickEnd += r_start - 3 * param.smoothBandwidth;
      }
		}

    // for thick region
    int r_summit;
    double r_summitVal;
    int r_thickStart;
    int r_thickEnd;
    // end for thick region

		int r_start;
		int r_end;
		int r_length;
		int r_lengthEx;
		int r_count;
		char r_strand;
		int r_chrIdx;
		double r_dense;
		double r_avgPos; // weigh average position
		double pval;
		double foldChange;
		double qval1;
		double qval2;
		bool bindRegion;
		vector<double> r_vec;

		double r_maxR2;
		int r_maxR2Pos;

		bool sortByDense (const readsCluster &rc ) const 
		{
			return r_dense > rc.r_dense;
		}
		bool sortByPos (const readsCluster &rc ) const
		{
			return r_chrIdx < rc.r_chrIdx || (r_chrIdx == rc.r_chrIdx && r_avgPos < rc.r_avgPos );
		}
		bool sortByPval (const readsCluster &rc ) const
		{
			return pval < rc.pval;
		}
	};
	vector<readsCluster> readsClusterVec;
	vector<readsCluster> SWReadClusterVec;

	//struct bingdingSite 
	//{
	//	int pos;
	//	double slope;
	//	double R2;
	//	char strand;
	//};
	//vector< struct bingdingSite > bsV;
	//bool pos_compare(const struct bingdingSite &a, const struct bingdingSite &b)
	//{
	//	return a.pos < b.pos;
	//}

	// read short reads in bed file.
	void readReads(char *filename);
	// sort read position from small to big
	void sortReads();
	// cluster read positions to read clusters
	void clusterReads();
	// sort the readsClusterVec, and get the top dense ones to estimate the mean of fragment length
	void estimateFragLengthMean();

	// call peaks only considering the uniform background // added in v.1.0
	void callPeaks();

	// detect binding sites
	void detectBS();
	// detect 2 binding sites once (combination information)
	void detectBS2();
	// deal with multi-binding sites
	void analyzeMultiBS();
	// write bed formatted binding regions to file
	void writeBRBed();
	// write bed formatted binding clusters to file
	void writeBSBed();
	// write bar formatted binding sites to file
	void writeBSBar();

public:
	~reads(void);
};


#endif

