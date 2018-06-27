#define DEBUG_R111
#include "reads.h"
reads::~reads(void)
{
}

reads::reads(void)
{
	readsCnt = 0;
	fwdReadsCnt = 0;
	rvsReadsCnt = 0;
	multiBSconut = 0;
	singleBScount = 0;
	discarderCnt = 0;
	//ambiguousBScount = 0;
	int i;
	for (i=0; i<6; i++)
		bCountArray[i] = 0;
	maxSlope = 0;
	maxR2 = 0;
	maxFC = 0;
}

void reads::readReads(char *filename)
{
	FILE *readF;
	char *line;
	int colNum = 0;
	int strandPos = 0;
	int i, pos;
	bool firstLine = true;
	readF = mustOpen(filename, "r"); 
	//int chrNo; //chromosome number, without "chr"
	char *chr;
	char **cols; //save temp columns
	int bufsize = 64;

	while ((line=readLine(readF)) != NULL)
	{
		readsCnt ++;
		if (line[0] == '/' || line[0] == '#')
		{
			continue;
		}

		if (firstLine)
		{
			colNum = bedColNum(line);
		}

		cols = new char*[colNum];
		for (i=0; i<colNum; i++)
			cols[i] = new char[bufsize];
		chopLine(line, colNum, cols);

		if (firstLine)
		{
			firstLine = false;
			for (i=3; i<colNum; i++)
			{
				if (isStrand(cols[i]))
				{
					strandPos = i;
					break;
				}
			}
#ifdef DEBUG_R
			printf("strand is col: %d\n", i);
#endif
		}

		// deal with the first col in bed file: chromosome
		chr = cols[0];
		int chrInd = -1;
		map_iter = chrIdxMap.find(chr);
		if (map_iter == chrIdxMap.end()) // new chromosome
		{
			chrInd = chrIdxMap.size();
			chrIdxMap[chr] = chrInd;
			fpos.push_back(vector<int>());
			rpos.push_back(vector<int>());
			chrName.push_back(chr);
#ifdef DEBUG_R
			printf("%d: %s\n", chrInd, chr);
#endif
		}
		else
			chrInd = map_iter->second;

		if (cols[strandPos][0] == '+') //forward strand, save start point to fpos
		{
			fwdReadsCnt ++;
			pos = atoi(cols[1]);
			(fpos[chrInd]).push_back(pos);
#ifdef DEBUG_R
			printf("\t%s\t%d\t+\n", chr, pos);
#endif
		}
		else  //reverse strand, save end point to rpos
		{
			rvsReadsCnt ++;
			pos = atoi(cols[2]);
			(rpos[chrInd]).push_back(pos);
#ifdef DEBUG_R
			printf("\t%s\t%d\t-\n", chr, pos);
#endif
		}
		for (i=0; i<colNum; i++)
			delete[] cols[i];
		delete[] cols;
		if (line != NULL)
			free(line);
		
		// print the reading progress
		if (param.VERBOSE && readsCnt % 1000 == 0)
			fprintf(stderr, "  %d reads done...\r", readsCnt);
	}
	totalChr = chrIdxMap.size();
	//lambda_bg_f = param.bandwidth * fwdReadsCnt / param.genomeSize;
	//lambda_bg_r = param.bandwidth * rvsReadsCnt / param.genomeSize;
}

void reads::sortReads()
{
	int chrIndex;
	for (chrIndex=0; chrIndex<totalChr; chrIndex++)
	{
		sort(fpos[chrIndex].begin(), fpos[chrIndex].end());
		sort(rpos[chrIndex].begin(), rpos[chrIndex].end());
	}
}

void reads::clusterReads()
{
	int chrIndex;
	int curPos, lastPos;
	vector<int>::iterator it;
	vector<int> tmp;

	for (map_iter=chrIdxMap.begin(); map_iter!=chrIdxMap.end(); map_iter++)
	{
		chrIndex = map_iter->second;

		/******* for forward reads ********/
		lastPos = 0;
		for (it=fpos[chrIndex].begin(); it!=fpos[chrIndex].end(); it++)
		{
			curPos = *it;
			
			if (param.filterDuplicateRead) // filter the duplicate reads
				if (curPos == lastPos)
					continue;

			if (curPos - lastPos > param.minDist)
			{
				if (tmp.size() >= param.minReadCnt && tmp.back()-tmp.front() > param.minClusterLength && tmp.size() > minCntGivenFDR(readsCnt, param.genomeSize, param.pval_cutoff, tmp.back()-tmp.front()+2*param.minDist))
				{
					// add to readsCluster vector
					readsCluster rctmp(tmp, chrIndex, '+');
					readsClusterVec.push_back(rctmp);
				}
				else
					discarderCnt += tmp.size();
				tmp.clear();
			}
			tmp.push_back(curPos);
			lastPos = curPos;
		}
		if (tmp.size() >= param.minReadCnt  && tmp.back()-tmp.front() > param.minClusterLength && tmp.size() > minCntGivenFDR(readsCnt, param.genomeSize, param.pval_cutoff, tmp.back()-tmp.front()+2*param.minDist))
		{
			// add to readsCluster vector
			readsCluster rctmp(tmp, chrIndex, '+');
			readsClusterVec.push_back(rctmp);
		}
		else
			discarderCnt += tmp.size();
		tmp.clear();
		/***** end of for forward reads ****/

		/******* for reverse reads ********/
		lastPos = 0;
		for (it=rpos[chrIndex].begin(); it!=rpos[chrIndex].end(); it++)
		{
			curPos = *it;
			
			if (param.filterDuplicateRead) // filter the duplicate reads
				if (curPos == lastPos)
					continue;

			if (curPos - lastPos > param.minDist)
			{
				if (tmp.size() >= param.minReadCnt  && tmp.back()-tmp.front() > param.minClusterLength && tmp.size() > minCntGivenFDR(readsCnt, param.genomeSize, param.pval_cutoff, tmp.back()-tmp.front()+2*param.minDist))
				{
					// add to readsCluster vector
					readsCluster rctmp(tmp, chrIndex, '-');
					readsClusterVec.push_back(rctmp);
				}
				else
					discarderCnt += tmp.size();
				tmp.clear();
			}
			tmp.push_back(curPos);
			lastPos = curPos;
		}
		if (tmp.size() >= param.minReadCnt  && tmp.back()-tmp.front() > param.minClusterLength && tmp.size() > minCntGivenFDR(readsCnt, param.genomeSize, param.pval_cutoff, tmp.back()-tmp.front()+2*param.minDist))
		{
			// add to readsCluster vector
			readsCluster rctmp(tmp, chrIndex, '-');
			readsClusterVec.push_back(rctmp);
		}
		else
			discarderCnt += tmp.size();
		tmp.clear();
		/***** end of for reverse reads ****/

#ifdef DEBUG_R11
		int i, j;
		printf("chrInd %d: +\n", chrIndex);
		for (i=0; i<fisland[chrIndex].size(); i++)
		{
			for (j=0; j<fisland[chrIndex][i].size(); j++)
			{
				printf("%d\t", fisland[chrIndex][i][j]);
			}
			printf("\n%d\n", fisland[chrIndex][i].size());
		}
		printf("chrInd %d: -\n", chrIndex);
		for (i=0; i<risland[chrIndex].size(); i++)
		{
			for (j=0; j<risland[chrIndex][i].size(); j++)
			{
				printf("%d\t", risland[chrIndex][i][j]);
			}
			printf("\n%d\n", risland[chrIndex][i].size());
		}
#endif
		//clearVector(fpos[chrIndex]);
		//clearVector(rpos[chrIndex]);
	}
	//clearVector(fpos);
	//clearVector(rpos);
}

void reads::callPeaks()
{
	int chrIndex;
	//int curPos, lastPos;
	vector<int>::iterator it, itBegin, itEnd;
	vector<int> tmp;
	int curChrSize;
	//int beginPos, endPos;

	// calculate the background model
	double unifBGFwdLambda = fwdReadsCnt * param.bandwidth / (double) param.genomeSize;
	double unifBGRvsLambda = rvsReadsCnt * param.bandwidth / (double) param.genomeSize;
	unsigned int minCntFwd = pois_cdf_inv(unifBGFwdLambda, 1.0 - param.pval_cutoff);
	unsigned int minCntRvs = pois_cdf_inv(unifBGRvsLambda, 1.0 - param.pval_cutoff);
	//printf("minCntFwd: %d\tminCntRvs: %d\n", minCntFwd, minCntRvs);

	for (map_iter=chrIdxMap.begin(); map_iter!=chrIdxMap.end(); map_iter++)
	{
		chrIndex = map_iter->second;

		/******* for forward reads ********/
		curChrSize = fpos[chrIndex].size();
		itBegin = fpos[chrIndex].begin();
		itEnd = itBegin;
		while (itEnd < fpos[chrIndex].end())
		{
			// extend itEnd
			while (distance(itBegin, itEnd) < (minCntFwd - 1) )
			{
				itEnd ++;
				if (itEnd == fpos[chrIndex].end())
					break;
				while (*itEnd - *itBegin > param.bandwidth)
				{
					//printf("%d\n", *itBegin);
					itBegin ++;
				}
			}
			if (itEnd == fpos[chrIndex].end())
				break;
			if (*itEnd - *itBegin <= param.bandwidth)
			{
				tmp.push_back(*itBegin);
				itBegin ++;
				itEnd ++;
			}
			else 
			{
				for (it=itBegin; it<itEnd; it++)
					tmp.push_back(*it);
				itBegin = itEnd;
				
				if (tmp.size() >= minCntFwd)
				{
					readsCluster rctmp(tmp, chrIndex, '+');
					readsClusterVec.push_back(rctmp);
					//printf("put: %d\n", tmp.size());
				}
				//printVec(tmp);
				tmp.clear();
			}
		}
		if (itBegin != fpos[chrIndex].end())
			for (it=itBegin; it<itEnd; it++)
				tmp.push_back(*it);
		if (tmp.size() >= minCntFwd  && tmp.size() > 0)
		{
			readsCluster rctmp(tmp, chrIndex, '+');
			readsClusterVec.push_back(rctmp);
			//printf("put: %d\n", tmp.size());
		}
		//printVec(tmp);
		tmp.clear();
		//printf("num of read cluster: %d\n", readsClusterVec.size());
		/***** end of for forward reads ****/

		/******* for reverse reads ********/
		curChrSize = rpos[chrIndex].size();
		itBegin = rpos[chrIndex].begin();
		itEnd = itBegin;
		while (itEnd < rpos[chrIndex].end())
		{
			// extend itEnd
			while (distance(itBegin, itEnd) < (minCntRvs - 1) )
			{
				itEnd ++;
				if (itEnd == rpos[chrIndex].end())
					break;
				while (*itEnd - *itBegin > param.bandwidth)
				{
					//printf("%d\n", *itBegin);
					itBegin ++;
				}
			}
			if (itEnd == rpos[chrIndex].end())
				break;
			if (*itEnd - *itBegin <= param.bandwidth)
			{
				tmp.push_back(*itBegin);
				itBegin ++;
				itEnd ++;
			}
			else 
			{
				for (it=itBegin; it<itEnd; it++)
					tmp.push_back(*it);
				itBegin = itEnd;

				if (tmp.size() >= minCntRvs)
				{
					readsCluster rctmp(tmp, chrIndex, '-');
					readsClusterVec.push_back(rctmp);
					//printf("put: %d\n", tmp.size());
				}
				//printVec(tmp);
				tmp.clear();
			}
		}
		for (it=itBegin; it<itEnd; it++)
			tmp.push_back(*it);
		if (tmp.size() >= minCntRvs && tmp.size() > 0)
		{
			readsCluster rctmp(tmp, chrIndex, '-');
			readsClusterVec.push_back(rctmp);
			//printf("put: %d\n", tmp.size());
		}
		//printVec(tmp);
		if (tmp.size() != 0)
			tmp.clear();
		//printf("num of read cluster: %d\n", readsClusterVec.size());
		/***** end of for reverse reads ****/

#ifdef DEBUG_R11
		int i, j;
		printf("chrInd %d: +\n", chrIndex);
		for (i=0; i<fisland[chrIndex].size(); i++)
		{
			for (j=0; j<fisland[chrIndex][i].size(); j++)
			{
				printf("%d\t", fisland[chrIndex][i][j]);
			}
			printf("\n%d\n", fisland[chrIndex][i].size());
		}
		printf("chrInd %d: -\n", chrIndex);
		for (i=0; i<risland[chrIndex].size(); i++)
		{
			for (j=0; j<risland[chrIndex][i].size(); j++)
			{
				printf("%d\t", risland[chrIndex][i][j]);
			}
			printf("\n%d\n", risland[chrIndex][i].size());
		}
#endif
		//clearVector(fpos[chrIndex]);
		//clearVector(rpos[chrIndex]);
	}
	//clearVector(fpos);
	//clearVector(rpos);
}

void reads::estimateFragLengthMean()
{
	int i;
	int totalRC = readsClusterVec.size();
	int topNo = (int) (param.topPercent * totalRC / 100.0);
	vector<readsCluster> sub_rcv;

	sort(readsClusterVec.begin(), readsClusterVec.end(), mem_fun_ref(&readsCluster::sortByDense));
	sort(readsClusterVec.begin(), readsClusterVec.begin()+topNo, mem_fun_ref(&readsCluster::sortByPos));
	for (i=0; i<(topNo-1); i++)
	{
		if (readsClusterVec[i].r_chrIdx == readsClusterVec[i+1].r_chrIdx && readsClusterVec[i].r_end > readsClusterVec[i+1].r_start)
		{
			sub_rcv.push_back(readsClusterVec[i++]);
			sub_rcv.push_back(readsClusterVec[i]);
		}
	}
	totalRC = sub_rcv.size();
	double avgPos, lastAvgPos = 0;
	vector<double> fragLengthSample;
	for (i=0; i<totalRC; i++)
	{
		avgPos = sub_rcv[i].r_avgPos;
		if (i%2==0)
			lastAvgPos = avgPos;
		else
			fragLengthSample.push_back(avgPos - lastAvgPos);
		//printf("%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%.1f\t%.lf\t%c\n", 
    //    sub_rcv[i].r_chrIdx, sub_rcv[i].r_start, sub_rcv[i].r_end, sub_rcv[i].r_length, sub_rcv[i].r_count, sub_rcv[i].r_dense, sub_rcv[i].r_thickStart, sub_rcv[i].r_thickEnd, sub_rcv[i].r_summit, sub_rcv[i].r_summitVal, avgPos, sub_rcv[i].r_strand);
	}
	if (fragLengthSample.size() == 0)
	{
		fprintf(stderr, "WARNING: not enough binding regions to estimate the everage fragment length.\n");
		param.fragLengthMean = -1;
		return;
	}
	//printVec(fragLengthSample);
	//printf("%f\t%f\n", avg(fragLengthSample), var(fragLengthSample));
	param.fragLengthMean = avg(fragLengthSample); // avg(fragLengthSample) under estimates the fragment length mean, because of the truncation at the tail and overlap bindings.
}

void reads::detectBS()
{
	param.fragLengthMean += param.smoothArm; // add the soomth arm to the average length
	//param.fragLengthVar  = param.fragLengthMean * param.fragLengthMean / 8.0;
	param.fragLengthVar  = param.fragLengthMean * param.fragLengthMean / 10.0;
	//param.fragLengthVar  = param.fragLengthMean * 20.0;
	param.lmLength = (int) (0.8 * param.fragLengthMean);
	param.lmArm = param.lmLength / 2;
	vector<int> f_maxPos, r_maxPos;
	vector<double> f_max, r_max, tmp, tmp2;
	//double R2_theshold = param.R2_cutoff;
	int i, j;
	int totalRC = readsClusterVec.size();
	//int arm = param.smoothArm + param.lmArm;
  int arm;
  if (param.smoothArm > 0) 
	  arm = param.smoothArm + param.lmArm + 5; //param.motifWidth / 4; // add the half motif width for shift the binding site a litter farther from the signal
  else
	  arm = 3 * param.smoothBandwidth + param.lmArm + 5; //param.motifWidth / 4;
	vector<readsCluster> sub_rcv;
	char strand;
	double pval;
	double foldChange;
	double qval1, qval2;
	int chrIdx;
	int pos = 0;
	double R2 = 0;
	double slope = 0;

	if (totalRC == 0)
		return;

	// for linear model
	lm *lm_signal;
	lm_signal = new lm();
	lm_signal->generateX();

	sort(readsClusterVec.begin(), readsClusterVec.end(), mem_fun_ref(&readsCluster::sortByPos));
	for (i=0; i<totalRC; i++)
	{
		strand = readsClusterVec[i].r_strand;
		chrIdx = readsClusterVec[i].r_chrIdx;
		pval = readsClusterVec[i].pval;
		foldChange = readsClusterVec[i].foldChange;
		qval1 = readsClusterVec[i].qval1;
		qval2 = readsClusterVec[i].qval2;
		if (strand == 'f' || strand == '+')
		{

			if (lm_signal->slipSolve(readsClusterVec[i].r_vec, 'f'))
			{
        logVec(lm_signal->slopeVec, tmp);
				multiVec(tmp, lm_signal->R2Vec, tmp2);

        //printf("fwd R2 slope log_slope R2*log_slope vecs:\n");
        //printVec(lm_signal->R2Vec);
        //printVec(lm_signal->slopeVec);
        //printVec(tmp);
        //printVec(tmp2);

				getTopN(tmp2, lm_signal->slopeVec, f_max, f_maxPos, 2, param.motifWidth);
				//getTopN(lm_signal->R2Vec, lm_signal->slopeVec, f_max, f_maxPos, 2, param.motifWidth);
				for (j=0; j<(int)f_maxPos.size(); j++)
				{
					pos = readsClusterVec[i].r_end + arm - f_maxPos[j];
					R2 = lm_signal->R2Vec[f_maxPos[j]];
					slope = lm_signal->slopeVec[f_maxPos[j]];

		      bingdingSite newBS(chrIdx, pos, R2, slope, strand, pval, foldChange, qval1, qval2);
		      BSVec.push_back(newBS);
		      if (slope > maxSlope)
			      maxSlope = slope;
		      if (R2 > maxR2)
			      maxR2 = R2;
		      if (foldChange > maxFC)
			      maxFC = foldChange;

					/*if (lm_signal->R2Vec[f_maxPos[j]] > R2_theshold)
					{
						bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_end + arm - f_maxPos[j], lm_signal->R2Vec[f_maxPos[j]], lm_signal->slopeVec[f_maxPos[j]], '+');
						//printf("%d\t%d\t%d\t%d\t%f\t%f\t+\n",readsClusterVec[i].r_end, arm, f_maxPos[j], readsClusterVec[i].r_end + arm - f_maxPos[j], lm_signal->R2Vec[f_maxPos[j]], lm_signal->slopeVec[f_maxPos[j]]);
						BSVec.push_back(newBS);
						if (lm_signal->slopeVec[f_maxPos[j]] > maxSlope)
							maxSlope = lm_signal->slopeVec[f_maxPos[j]];
					}*/
				}
			}
		}
		if (strand == 'r' || strand == '-')
		{
			if (lm_signal->slipSolve(readsClusterVec[i].r_vec, 'r'))
			{
        logVec(lm_signal->slopeVec, tmp);
				multiVec(tmp, lm_signal->R2Vec, tmp2);

        //printf("rvs R2 slope log_slope R2*log_slope vecs:\n");
        //printVec(lm_signal->R2Vec);
        //printVec(lm_signal->slopeVec);
        //printVec(tmp);
        //printVec(tmp2);

				getTopN(tmp2, lm_signal->slopeVec, r_max, r_maxPos, 2, param.motifWidth);
				//getTopN(lm_signal->R2Vec, lm_signal->slopeVec, r_max, r_maxPos, 2, param.motifWidth);
				for (j=0; j<(int)r_maxPos.size(); j++)
				{
					pos = readsClusterVec[i].r_start - arm + r_maxPos[j];
					R2 = lm_signal->R2Vec[r_maxPos[j]];
					slope = lm_signal->slopeVec[r_maxPos[j]];

		      bingdingSite newBS(chrIdx, pos, R2, slope, strand, pval, foldChange, qval1, qval2);
		      BSVec.push_back(newBS);
		      if (slope > maxSlope)
			      maxSlope = slope;
		      if (R2 > maxR2)
			      maxR2 = R2;
		      if (foldChange > maxFC)
			      maxFC = foldChange;

					/*if (lm_signal->R2Vec[r_maxPos[j]] > R2_theshold)
					{
						bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_start - arm + r_maxPos[j], lm_signal->R2Vec[r_maxPos[j]], lm_signal->slopeVec[r_maxPos[j]], '-');
						//printf("%d\t%d\t%d\t%d\t%f\t%f\t-\n",readsClusterVec[i].r_start, arm, r_maxPos[j], readsClusterVec[i].r_start - arm + r_maxPos[j], lm_signal->R2Vec[r_maxPos[j]], lm_signal->slopeVec[r_maxPos[j]]);
						BSVec.push_back(newBS);
						if (lm_signal->slopeVec[r_maxPos[j]] > maxSlope)
							maxSlope = lm_signal->slopeVec[r_maxPos[j]];
					}*/
				}
			}
		}
		//chrIdx(chrIdx), pos(pos), R2(R2), slope(slope), strand(strand), pval(pval), foldChange(foldChange),qval(qval)

		/*if (param.VERBOSE && i % 100 == 0)
			printf("  %.2f%% done...\r", 100.0 * (double) i / (double) totalRC);*/
	}
}

void reads::detectBS2()
{
	param.fragLengthMean += param.smoothArm;
	param.fragLengthVar  = param.fragLengthMean * param.fragLengthMean / 8.0;
	param.lmLength = (int) (0.8 * param.fragLengthMean);
	param.lmArm = param.lmLength / 2;
	vector<int> f_maxPos, r_maxPos, distBSIdx;
	//int f_maxPos, r_maxPos, distBS;
	vector<double> f_max, r_max, tmp;
	//double f_max, r_max;
	double R2_theshold = param.R2_cutoff;
	int i, j;
	int totalRC = readsClusterVec.size();
	//int arm = param.smoothArm + param.lmArm;
	vector<readsCluster> sub_rcv;
	char strand;

	// for linear model
	lm *lm_signal;
	lm_signal = new lm();
	lm_signal->generateX();

	sort(readsClusterVec.begin(), readsClusterVec.end(), mem_fun_ref(&readsCluster::sortByPos));
	for (i=0; i<totalRC; i++)
	{
		strand = readsClusterVec[i].r_strand;
		if (strand == 'f' || strand == '+')
		{
			if (lm_signal->slipSolve2BS(readsClusterVec[i].r_vec, 'f'))
			{
				//multiVec(lm_signal->R2Vec, lm_signal->slopeVec, tmp);
				//getTopN(tmp, f_max, f_maxPos, 2, 0.01);
				//getTopN(lm_signal->R2Vec, f_max, f_maxPos, 2, param.R2_cutoff);
				get2BSTop(lm_signal->R2Mat, f_max, f_maxPos, distBSIdx);
				for (j=0; j<(int)f_maxPos.size(); j++)
				{
					if (lm_signal->R2Mat[distBSIdx[j]][f_maxPos[j]] > R2_theshold)
					{
						//printf("%d\n", distBSIdx[j]);
//						bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_end + arm - f_maxPos[j], lm_signal->R2Mat[distBSIdx[j]][f_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]], '+');
						//printf("%d\t%f\t%f\t+\n",readsClusterVec[i].r_end + arm - f_maxPos[j], lm_signal->R2Mat[distBSIdx[j]][f_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]]);
//						BSVec.push_back(newBS);
						if (distBSIdx[j] > 0)
						{
//							bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_end + arm - f_maxPos[j] - lm_signal->distVec[distBSIdx[j]], lm_signal->R2Mat[distBSIdx[j]][f_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]], '+');
							//printf("%d\t%f\t%f\t+\n",readsClusterVec[i].r_end + arm - f_maxPos[j] - lm_signal->distVec[distBSIdx[j]], lm_signal->R2Mat[distBSIdx[j]][f_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]]);
//							BSVec.push_back(newBS);
						}
						if (lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]] > maxSlope)
							maxSlope = lm_signal->slopeMat[distBSIdx[j]][f_maxPos[j]];
					}
				}
			}
		}
		if (strand == 'r' || strand == '-')
		{
			if (lm_signal->slipSolve2BS(readsClusterVec[i].r_vec, 'r'))
			{
				//multiVec(lm_signal->R2Vec, lm_signal->slopeVec, tmp);
				//getTopN(tmp, r_max, r_maxPos, 2, 0.01);
				//getTopN(lm_signal->R2Vec, r_max, r_maxPos, 2, param.R2_cutoff);
				get2BSTop(lm_signal->R2Mat, r_max, r_maxPos, distBSIdx);
				for (j=0; j<(int)r_maxPos.size(); j++)
				{
					if (lm_signal->R2Mat[distBSIdx[j]][r_maxPos[j]] > R2_theshold)
					{
						//printf("%d\n", distBSIdx[j]);
//						bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_start - arm + r_maxPos[j], lm_signal->R2Mat[distBSIdx[j]][r_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]], '-');
						//printf("%d\t%f\t%f\t-\n",readsClusterVec[i].r_start - arm + r_maxPos[j], lm_signal->R2Mat[distBSIdx[j]][r_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]]);
//						BSVec.push_back(newBS);
						if (distBSIdx[j] > 0)
						{
//							bingdingSite newBS(readsClusterVec[i].r_chrIdx, readsClusterVec[i].r_start - arm + r_maxPos[j] + lm_signal->distVec[distBSIdx[j]], lm_signal->R2Mat[distBSIdx[j]][r_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]], '-');
							//printf("%d\t%f\t%f\t-\n",readsClusterVec[i].r_start - arm + r_maxPos[j] + lm_signal->distVec[distBSIdx[j]], lm_signal->R2Mat[distBSIdx[j]][r_maxPos[j]], lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]]);
//							BSVec.push_back(newBS);
						}
						if (lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]] > maxSlope)
							maxSlope = lm_signal->slopeMat[distBSIdx[j]][r_maxPos[j]];
					}
				}
			}
		}
	}
}

void reads::analyzeMultiBS()
{
	int totalBS, i, lastPos, curPos, midPos, curIdx, lastIdx;
	char midStrand;
	double midR2, midSlope;
	double minPval, maxFC, minQval1, minQval2;
	vector<bingdingSite> BStmp;

	if (BSVec.size() == 0)
		return;

	sort(BSVec.begin(), BSVec.end());
	totalBS = BSVec.size();

	lastPos = BSVec[0].pos;
	lastIdx = BSVec[0].chrIdx;
	BStmp.push_back(BSVec[0]);
	for (i=1; i<totalBS; i++)
	{
		//cout<<BSVec[i].chrIdx<<"\t"<<BSVec[i].pos<<"\t"<<BSVec[i].R2<<"\t"<<BSVec[i].strand<<endl;
		curPos = BSVec[i].pos;
		curIdx = BSVec[i].chrIdx;
		if (curIdx != lastIdx || curPos - lastPos > param.minBSclusterDist)
		{
			BSClusterVec.push_back(BStmp);
			BStmp.clear();
			BStmp.push_back(BSVec[i]);
		}
		else if (curPos - lastPos < param.motifWidth)
		{
			midPos = (int) ((lastPos * BSVec[i-1].slope + curPos * BSVec[i].slope) / (BSVec[i-1].slope + BSVec[i].slope));
			midR2 = (BSVec[i].R2 + BSVec[i-1].R2) / 2.0;
			midSlope = (BSVec[i].slope + BSVec[i-1].slope) / 2.0;
			minPval = min(BSVec[i].pval, BSVec[i-1].pval);
			maxFC = max(BSVec[i].foldChange, BSVec[i-1].foldChange);
			minQval1 = min(BSVec[i].qval1, BSVec[i-1].qval1);
			minQval2 = min(BSVec[i].qval2, BSVec[i-1].qval2);
			if (BSVec[i].strand == BSVec[i-1].strand)
			{
				midStrand = BSVec[i].strand;
				//printf("!! check reads.cpp Ln659 !!\n");
			}
			else
				midStrand = 'b';
			bingdingSite newMidBS(curIdx, midPos, midR2, midSlope, midStrand, minPval, maxFC, minQval1, minQval2);
			BStmp.pop_back();
			BStmp.push_back(newMidBS);
		}
		else
			BStmp.push_back(BSVec[i]);
		lastPos = curPos;
		lastIdx = curIdx;
	}
	BSClusterVec.push_back(BStmp);
	BStmp.clear();
}

void reads::writeBRBed()
{
	FILE *writeBedF = mustOpen(param.outbedfile, "w");
	int i, totalRC = (int)readsClusterVec.size();
	char strand;
	double pval;
	double foldChange;
	double qval1, qval2;
	int chrIdx;
	const char* chr;
	int start, end;
	char* name = (char*) malloc(sizeof(char)*256);
	int score;
  bool flag_out;
  int rc; // read count

	if (totalRC == 0)
		return;

	for (i=0; i<totalRC; i++)
	{
		strand = readsClusterVec[i].r_strand;
		chrIdx = readsClusterVec[i].r_chrIdx;
		chr = chrName[chrIdx].c_str();
    if(strand == 'f' | strand == '+') 
    {
		  start = readsClusterVec[i].r_thickStart + (int) (param.fragLengthMean / 2);
		  end = readsClusterVec[i].r_thickEnd + (int) (param.fragLengthMean /2);
     }
    else 
     {
		  start = readsClusterVec[i].r_thickStart - (int) (param.fragLengthMean / 2);
		  end = readsClusterVec[i].r_thickEnd - (int) (param.fragLengthMean / 2);
    }
		pval = readsClusterVec[i].pval;
		foldChange = readsClusterVec[i].foldChange;
		qval1 = readsClusterVec[i].qval1;
		qval2 = readsClusterVec[i].qval2;
    rc = readsClusterVec[i].r_count;

    if ( (strand == 'f' || strand == '+') && (i+1) < totalRC && ( readsClusterVec[i+1].r_strand == '-' || readsClusterVec[i+1].r_strand == 'r' ))
    {
      if (chrIdx == readsClusterVec[i+1].r_chrIdx)
      {
        if ( (readsClusterVec[i].r_end + (int) (param.fragLengthMean)/2) > (readsClusterVec[i+1].r_start - (int) (param.fragLengthMean)/2) )
        {
          i ++;
          flag_out = 1;
          strand = 'b';
          start = min(start, readsClusterVec[i].r_thickStart - (int) (param.fragLengthMean/2) );
          end = max(end, readsClusterVec[i].r_thickEnd - (int) (param.fragLengthMean/2) );

          if (param.filterSingleStrand)
          {
            pval = sqrt(pval * readsClusterVec[i].pval);
            foldChange = 0.5 * (foldChange +  readsClusterVec[i].foldChange); 
            qval1 = sqrt(qval1 * readsClusterVec[i].qval1);
            qval2 = sqrt(qval2 * readsClusterVec[i].qval2);
            rc = (rc + readsClusterVec[i].r_count) / 2;
          }
          else
          {
            pval = min(pval, readsClusterVec[i].pval);
            foldChange = max(foldChange, readsClusterVec[i].foldChange); 
            qval1 = min(qval1, readsClusterVec[i].qval1);
            qval2 = min(qval2, readsClusterVec[i].qval2);
            rc = (rc + readsClusterVec[i].r_count) / 2;
          }
        }
        else
          flag_out = 0;
      }
      else 
        flag_out = 0;
    }
    else flag_out = 0;

    if (param.filterSingleStrand && ! flag_out) 
    {
      continue;
    }
    if (start < 0) start = 0;
		score = (int)(foldChange / maxFC * 1000);
		//sprintf(name, "%d|%.2f|%g|%g|%g", rc, foldChange,pval,qval1,qval2);
		sprintf(name, "%d|%.2f|%g|%g", rc, foldChange,pval,qval1);
		//fprintf(writeBedF, "%s\t%d\t%d\t%s\t%d\t%c\n", chr, start, end, name, score, strand);
		fprintf(writeBedF, "%s\t%d\t%d\t%s\t%d\t+\n", chr, start, end, name, score);
	}
}		

void reads::writeBSBed()
{
	FILE *writeBedF = mustOpen(param.outbedfile, "w");

	int totalBSC, totalBS, i, j;
	const char *chr;

	char* name = (char*) malloc(sizeof(char)*32);
	int score;
	char* itemRgb = "0";
	int start;
	int end;
	int extent = param.BSextent;
	double dtmp;
	int itmp;
	char* blockSize = "1,";
	char* blockSizes = (char*) malloc(sizeof(char)*256);
	char* blockStarts = (char*) malloc(sizeof(char)*256);

	//sort(chrName.begin(), chrName.end());

	totalBSC = BSClusterVec.size();
	for (i=0; i<totalBSC; i++)
	{
		totalBS = BSClusterVec[i].size();
		if (totalBS == 0)
		{
			printf("ambiguous error! check the code\n");
			exit(1);
		}

		if (totalBS == 1)
			singleBScount ++;
		else
			multiBSconut ++;

		// # BS clusters with n bindings
		if (totalBS > 5)
			bCountArray[0] ++;
		else
			bCountArray[totalBS] ++;

		chr = chrName[BSClusterVec[i][0].chrIdx].c_str();

		start = BSClusterVec[i][0].pos - extent;
		if (start < 0)
			start = 0;
		end  = BSClusterVec[i][totalBS-1].pos + extent;
		sprintf(name, "SeqSite-%d", i+1);
		score = 0;
		blockSizes[0] = '\0';
		blockStarts[0] = '\0';
		for (j=0; j<totalBS; j++)
		{
			//dtmp = BSClusterVec[i][j].slope / maxSlope * 100;
			//itmp = (int) (dtmp * 10);
			//if (score < itmp)
			//	score = itmp;
			dtmp = BSClusterVec[i][j].foldChange * BSClusterVec[i][j].slope / maxFC * 100 / maxSlope * 100;
			itmp = (int) (dtmp * 10);
			if (score < itmp)
				score = itmp;
			blockSizes = strcat(blockSizes, blockSize);
			sprintf(blockStarts, "%s%d,", blockStarts, BSClusterVec[i][j].pos - start);
		}
		
		fprintf(writeBedF, "%s\t%d\t%d\t%s\t%d\t+\t%d\t%d\t%s\t%d\t%s\t%s\n", chr, start, end, name, score, start, end, itemRgb, totalBS, blockSizes, blockStarts);

		/*
		bar file refer to: 
			http://biogibbs.stanford.edu/%7Ejiangh/browser/README.html#formats
		bed file refer to: 
			http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED
		*/
	}
}

void reads::writeBSBar()
{
	FILE *writeBarF = mustOpen(param.outbarfile, "w");
	int totalBSC, totalBS, i, j, ii;
	const char *chr;
	char tmpStrand;
	int type; // binding site type: sigle direction bs in a cluster: 0 for forward, 1 for revers; both direction bs in a cluster: 2 for no merged bs, 4 for merged bs

	//sort(chrName.begin(), chrName.end());

	totalBSC = BSClusterVec.size();
	for (i=0; i<totalBSC; i++)
	{
		totalBS = BSClusterVec[i].size();
		chr = chrName[BSClusterVec[i][0].chrIdx].c_str();
		for (j=0; j<totalBS; j++)
		{
			type = 0;
			tmpStrand = BSClusterVec[i][0].strand;
			for (ii=1; ii<totalBS; ii++)
			{
				if (tmpStrand != BSClusterVec[i][ii].strand)
				{
					type = 2;
					break;
				}
			}
			for (ii=0; ii<totalBS; ii++)
				if (BSClusterVec[i][ii].strand == 'f'|| BSClusterVec[i][ii].strand == '+')
					BSClusterVec[i][ii].type = type + 0;
				else if (BSClusterVec[i][ii].strand == 'r'|| BSClusterVec[i][ii].strand == '-')
					BSClusterVec[i][ii].type = type + 1;

			for (ii=0; ii<totalBS; ii++)
				if (BSClusterVec[i][ii].strand == 'b')
					BSClusterVec[i][ii].type = 4;

			//fprintf(writeBarF, "%s\t%d\t%d\t%e\t%f\t%e\t%e\t%f\t%f\n", chr, BSClusterVec[i][j].pos, BSClusterVec[i][j].type, BSClusterVec[i][j].pval, BSClusterVec[i][j].foldChange, BSClusterVec[i][j].qval1, BSClusterVec[i][j].qval2, 
			fprintf(writeBarF, "%s\t%d\t%e\t%f\t%e\t%f\t%f\n", chr, BSClusterVec[i][j].pos, BSClusterVec[i][j].pval, BSClusterVec[i][j].foldChange, BSClusterVec[i][j].qval1, 
				BSClusterVec[i][j].R2, BSClusterVec[i][j].slope / maxSlope * 100);
			//printf("%s\t%d\t%f\n", chr, BSClusterVec[i][j].pos, BSClusterVec[i][j].slope / maxSlope * 100);
		}
	}
}
