#include <iostream>
#include <string>
#include "global.h"
#include "reads.h"
#include "seqsite.h"
#include "lm.h"
#include "bindSig.h"
#include "parameter.h"

using namespace std;
parameter param;


int main (int argc, char **argv) 
{
	// initialize default values
	param.VERBOSE = true;
	param.infile = NULL;
	param.outbedfile = NULL;
	param.outbarfile = NULL;
	param.controlfile = NULL;
	param.genomeSize = 2400000000;
	param.bandwidth = 100;
	param.fragLengthMean = 150;
	param.lmLength = (int) (0.8 * param.fragLengthMean);
	param.lmArm = param.lmLength / 2;
	param.pval_cutoff = 0.001;
	param.FDR = 0.1;
	param.minDist = 30;
	param.minReadCnt = 10;
	param.minClusterLength = 10;
	//param.minDense = 0.08;	/* IMPORTANT: depends on FDR, genome size, # reads */
	param.R2_cutoff = 0.2; // 0.2 can pass everything, not used
	param.smoothArm = 20;
	param.smoothBandwidth = 0;
	param.motifWidth = 20;  /* supposed motif width */
	param.BSextent = 50;
	param.minBSclusterDist = param.BSextent * 2;
	param.toEstimateFragLengthMean = true;
	param.useControl = false;
	param.writeBedFile = true;
	param.gridSearch = false;
	param.lambdaRL1 = 500;
	param.lambdaRL2 = 2000;
	param.lambdaRL3 = 10000;
	param.control2ChIPseqRatio = 1.0;
	param.control2ChIPseqRatioFwd = 1.0;
	param.control2ChIPseqRatioRvs = 1.0;
	param.topPercent = 5;
	param.filterDuplicateRead = false;
  param.filterSingleStrand = false;
	// end of initialization

	int i;

	if (param.smoothArm > param.minDist)
	{
		fprintf(stderr, "Smooth arm (%d) should be less than minimum clustering distance (%d).\n", param.smoothArm, param.minDist);
		exit(1);
	}
	// end param default

	parseParam(argc, argv);

	if (param.VERBOSE)
	{
		printf("\n");
		printf("The parameters specified are:\n");
		printf(
			"  The ChIP-seq data file:  %s\n", param.infile);
		if (param.useControl)
			printf(
			"  The control data file:   %s\n", param.controlfile);
		printf(
			"  The output BAR file:     %s\n", param.outbarfile);
		if (param.writeBedFile)
			printf(
			"  The output BED file:     %s\n", param.outbedfile);
		printf(
			"  Effective genome size:   %.2g bp\n", param.genomeSize);
		printf(
			"  P-value cutoff:          %g\n", param.pval_cutoff);
		if (param.toEstimateFragLengthMean)
			printf(
			"  DNA fragment length:     estimate from the data\n"
			"  Top percentage for FLE:  %d%%\n"
			"  * FLE: fragment length estimation\n"
			, param.topPercent);
		else
			printf(
			"  DNA fragment length:     %.2f (nt)\n", param.fragLengthMean);
		printf("\n");
	}

	reads *readsI;
	readsI = new reads();

	if (param.VERBOSE)
		printf("Reading ChIP-seq file...\n");
	readsI->readReads(param.infile);
	if (param.VERBOSE)
		printf("  Total ChIP-seq reads: %d (+:%d -:%d)\n", readsI->readsCnt, readsI->fwdReadsCnt, readsI->rvsReadsCnt);
	
	// the min dense according to the given FDR
	//param.minDense = minCntGivenFDR(readsI->readsCnt, param.genomeSize, param.FDR, (int) param.fragLengthMean ) / param.fragLengthMean;

	if (readsI->readsCnt < 1000 && param.toEstimateFragLengthMean)
	{
		fprintf(stderr, "ERROR: Need more reads to estimate the fragment length!\nExiting...\n");
		//exit(1);
	}
	if (param.VERBOSE)
		printf("Sorting ChIP-seq reads...\n");
	readsI->sortReads();
	if (param.VERBOSE)
		printf("Clustering ChIP-seq reads...\n");
	readsI->clusterReads();
	//readsI->callPeaks();

	//********** reading control file ***********//
	reads *readsControl;
	readsControl = new reads();

	if (param.useControl)
	{
		if (param.VERBOSE)
			printf("Reading control file...\n");
		readsControl->readReads(param.controlfile);
		if (param.VERBOSE)
			printf("  Total control reads: %d (+:%d -:%d)\n", readsControl->readsCnt, readsControl->fwdReadsCnt, readsControl->rvsReadsCnt);

		param.control2ChIPseqRatio = (double) readsI->readsCnt / readsControl->readsCnt;
		param.control2ChIPseqRatioFwd = (double) readsI->fwdReadsCnt / readsControl->fwdReadsCnt;
		param.control2ChIPseqRatioRvs = (double) readsI->rvsReadsCnt / readsControl->rvsReadsCnt;
		//if (param.VERBOSE)
		//	printf("  The normalizing factor from control to ChIP-seq data is: %.2f\n", param.control2ChIPseqRatio);

		if (param.VERBOSE)
			printf("Sorting control reads...\n");
		readsControl->sortReads();
	}
	//**********END OF reading control file ***********//
	
	if (param.VERBOSE)
		printf("Filtering read clusters using local background...\n");
	if (!param.useControl)
	{
		readsControl->fpos.resize(readsI->totalChr);
		readsControl->rpos.resize(readsI->totalChr);
		for (i=0; i<readsI->totalChr; i++)
		{
			readsControl->fpos[i].resize(readsI->fpos[i].size());
			copy(readsI->fpos[i].begin(), readsI->fpos[i].end(), readsControl->fpos[i].begin());
			readsControl->rpos[i].resize(readsI->rpos[i].size());
			copy(readsI->rpos[i].begin(), readsI->rpos[i].end(), readsControl->rpos[i].begin());
		}
	}
	peakFilter(readsI, readsControl);

	
	if (param.toEstimateFragLengthMean)
	{
		if (param.VERBOSE)
			printf("Estimating fragment length...\n");
		readsI->estimateFragLengthMean();
		if (param.VERBOSE && param.fragLengthMean > 0)
			printf("  The estimated fragment length is: %.2f (nt)\n", param.fragLengthMean);
		if (param.fragLengthMean < 60 || param.fragLengthMean > 500)
		{
			fprintf(stderr, "WARNING: Failed to estimate fragment length from data. \n"
				"120-bp is used as the default fragment length, or specify the fragment length yourself.\n");
			param.fragLengthMean = 120;
			//exit(1);
		}
	}
	if (param.VERBOSE)
		printf("Pinpointing binding sites...\n");
	if (!param.gridSearch)
		readsI->detectBS();
	else
		readsI->detectBS2();
	readsI->analyzeMultiBS();
	//if (param.VERBOSE)
		//printf("  %d (%.2f%%) reads are enriched in binding regions\n", readsI->readsCnt - readsI->discarderCnt, 
		//(1 - readsI->discarderCnt / (double) readsI->readsCnt) * 100);
	if (param.VERBOSE)
		printf("Writing output file...\n");
	readsI->writeBSBar();
	if (param.writeBedFile)
		readsI->writeBRBed();
//	if (param.VERBOSE)
//		printf(
//		"\n"
//		"Summary of the %d detected binding regions:\n"
//		"  # Binding regions having variant number of binding sites:\n"
//		"        1      2      3      4      5     5+\n"
//		"    %5d  %5d  %5d  %5d  %5d  %5d\n"
//		"Done.\n",
//		readsI->singleBScount + readsI->multiBSconut, 
//		readsI->bCountArray[1], readsI->bCountArray[2], readsI->bCountArray[3], 
//		readsI->bCountArray[4], readsI->bCountArray[5], readsI->bCountArray[0]);
	return 1;
}

void parseParam(int argc, char **argv)
{
	int i;
	//if (argc < 3)
		//errorMsg();

	for (i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			break;

		switch(argv[i][1])
		{
		case 'l':	// average DNA fragment length
			if(++i>=argc)
				errorMsg();
			param.fragLengthMean = atof(argv[i]);
			param.toEstimateFragLengthMean = false;
			break;
		case 'p':
			if(++i>=argc)
				errorMsg();
			param.pval_cutoff = atof(argv[i]);
			break;
		case 'f':
			if(++i>=argc)
				errorMsg();
			param.FDR = atof(argv[i]);
			break;
		case 'g':
			if(++i>=argc)
				errorMsg();
			param.genomeSize = atof(argv[i]);
			break;
		//case 'e':
		//	if(++i>=argc)
		//		errorMsg();
		//	param.BSextent = atoi(argv[i]);
		//	param.minBSclusterDist = param.BSextent * 2;
		//	break;
		case 't':
			if(++i>=argc)
				errorMsg();
			param.topPercent = atoi(argv[i]);
			break;
		case 's':
			if(++i>=argc)
				errorMsg();
			param.smoothArm = atoi(argv[i]);
			break;
		case 'k':
			if(++i>=argc)
				errorMsg();
			param.smoothBandwidth = atoi(argv[i]);
      param.smoothArm = 0;
			break;
		case 'w':
			if(++i>=argc)
				errorMsg();
			param.motifWidth = atoi(argv[i]);
			break;
		case 'd':	// clustering distance
			if(++i>=argc)
				errorMsg();
			param.minDist = atoi(argv[i]);
			break;
		case 'n':	// read count threshold in a read cluster
			if(++i>=argc)
				errorMsg();
			param.minReadCnt = atoi(argv[i]);
			break;
		case 'c':
			if(++i>=argc)
				errorMsg();
			param.controlfile = argv[i];
			param.useControl = true;
			break;
//		case 'b':	// output detailed binding sites in bar file
//			if(++i>=argc)
//				errorMsg();
//			param.outbedfile = argv[i];
//			param.writeBedFile = true;
//			break;
		case 'q':
			param.VERBOSE = false;
			break;
		case 'F':
			param.filterDuplicateRead = true;
			break;
		case 'S':
			param.filterSingleStrand = true;
			break;
		//case 'i':
		//	param.gridSearch = true;
		//	break;
		case 'v': // show version
			versionMsg();
			break;
		case 'h':	// show help message
			errorMsg();
			break;
		case 'a':
			authorMsg();
			break;
		default:
			errorMsg();
			break;
		}
	}
	if(i>=argc)
		errorMsg();
	param.infile = argv[i];
	if(++i>=argc)
		errorMsg();
	param.outbarfile = argv[i];
	if(++i>=argc)
		errorMsg();
	param.outbedfile = argv[i];
	
	if (0)
	{
		fprintf(stderr, "%s\t%s\t%f\t%f\t%d\t%d\n", param.infile,param.outbedfile, param.fragLengthMean, param.fragLengthVar,param.minDist, param.minReadCnt);
		exit(0);
	}
}

void errorMsg()
{
	fprintf(stderr, 
		"\n"
		"SeqSite: ChIP-Seq binding site identification\n"
		"\n"
		"USAGE: SeqSite [options] <input.bed> <output.bar> <output.bed>\n"
		"\tinput.bed    ChIP-seq data in BED format (4 fields required: chrId, start, end, and strand)\n"
		"\toutput.bar   BAR file containing binding sites identified\n"
		"\toutput.bed   BED file containing binding regions detected\n"
		// output file is bed formatted, with single, multiple, ambiguous binding events annotated. if(detail==true), multi-binding sites are specified in the file.
		// SIGNAL_CUTOFF 0.01	// if(#reads/(dist+FRAG_LENGTH)) < SIGNAL_CUTOFF) discard;
		// MULTI_CUTOFF 0.2	// if(#reads/(#dist+FRAG_LENGTH)) > MULTI_CUTOFF)  mutil-bs analysis;

		/*** TODO: 
		two types of output:
		bar: a line for a binding site, with the signal strength
		bed: a line for a binding cluster, with each binding location
		log file? : summary of the detecting
		***/

		"Options: (* advanced)\n"
		"\t-c <string>  control data in BED format (4 fields required) (default: not use)\n"
		"\t-g <int>     effective genome size (default: 2.4e+9 for the human genome)\n"
    "\t-d <int>     * tag clustering distance (default: 30)\n"
    "\t-n <int>     * min tag count in a tag cluster (default: 10)\n"
    "\t-S           * filter single-strand tag clusters (default: not filter)\n"
		"\t-l <double>  * average DNA fragment length (default: estimate from data)\n"
		"\t-t <int>     * top <int>%% tag clusters for frag. length estimating (default: 5)\n"
		"\t-p <double>  p-value cutoff for binding region detection (default: 1e-3)\n"
		"\t-f <double>  FDR for binding region detection (default: 0.1)\n"
		"\t-s <int>     * arm length for smoothing tag signal (default: 20)\n"
		"\t-k <int>     * kernel density bandwidth for smoothing tag signal (default: use -s)\n"
		"\t-w <int>     * experimental motif width (default: 20)\n"
		"\t-F           * filter out the duplicate reads (default: FALSE)\n"
		"\t-q           quiet: no screen display (default: show progress)\n"
		"Help Options:\n"
		"\t-h           show this help message\n"
		"\t-v           show version information\n"
		"\t-a           about SeqSite\n"
		"\n"
		);
	exit(1);
}

void authorMsg()
{
	fprintf(stderr, 
		"The authors:\n"
		"\tXi Wang & Xuegong Zhang\n"
		"Contact:\n"
		"\twang-xi05@mails.tsinghua.edu.cn\n"
		"Citation:\n"
		"\tXi Wang and Xuegong Zhang. Pinpointing transcription factor binding sites from ChIP-seq data with SeqSite. 2010. Submitted.\n"
		"\n"
		);
	exit(1);
}

void versionMsg()
{
	fprintf(stderr, 
		"SeqSite version 1.1.2 (Nov 18, 2010)\n"
		);
	exit(1);
}
