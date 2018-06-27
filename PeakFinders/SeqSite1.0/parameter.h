#ifndef _PARAM_H
#define _PARAM_H

typedef struct _parameter 
{
	double genomeSize;  // the effective genome size
	int minDist;		// the minimum gap size between two peaks
	unsigned int minReadCnt;	// the read count threshold for a peak
	int minClusterLength;	// the minimum peak size
	int bandwidth;		// the bandwidth for reads added to a peak
	double minDense;	//** not used **//
	double pval_cutoff;	// the p-value cutoff for calling peaks
	double FDR;			// false discovery rate
	double fragLengthMean;	// the average DNA fragment length
	double fragLengthVar;	//** not used **//
	double R2_cutoff;		// the R-square threshold for binding sites
	bool VERBOSE;		// TRUE for output verbose info
	char *infile;		// the name of the file to be treated
	char *controlfile;	// the name of the control file
	bool useControl;		// TRUE for using control
	double control2ChIPseqRatio;	// the normalizing factor from control data to ChIPseq data
	double control2ChIPseqRatioFwd;	// the normalizing factor for forward strand
	double control2ChIPseqRatioRvs;	// the normalizing factor for reverse strand
	char *outbedfile;	// the name for the output BED file
	bool writeBedFile;	// TRUE for output BAR file
	char *outbarfile;	// the name for the output BAR file
	bool filterDuplicateRead; // TRUE for filtering the duplicate reads
	int readLength;		//** not used **//
	int motifWidth;		// the width of motif
	int smoothArm;		// the arm length for smoothing the ChIP-seq signal
	int smoothBandwidth;		// the bandwidth for smoothing the ChIP-seq signal (adopted from QuEST)
	int lmLength;		// the curve fitting length
	int lmArm;			// the half length of lmLength
	int minBSclusterDist;	// the minimum gap size between two binding sites
	int BSextent;		// the half length of a binding region having one binding site
	bool toEstimateFragLengthMean; //if false, the average fragment length should be given by user
	bool gridSearch;	// TRUE for an alternative method (grid search) for locating binding site
	int topPercent;		// the percentage for estimating the DNA fragment length
	int lambdaRL1;		// the length of region centered at the peak summit for calculate local lambda_1
	int lambdaRL2;		// the length of region centered at the peak summit for calculate local lambda_2
	int lambdaRL3;		// the length of region centered at the peak summit for calculate local lambda_3
	bool filterSingleStrand;

}parameter;

#endif
