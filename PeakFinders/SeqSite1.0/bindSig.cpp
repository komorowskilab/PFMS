/*
bindRegion: for detecting binding regions w/wo control data 
*/

#include "bindSig.h"

void peakFilter(reads* ChIPseq, reads* control)
{
	int i;
	int *mapChr;
	mapChr = new int[ChIPseq->totalChr];
	int filteredNum = 0;
	int filteredReadCnt = 0;

	if (param.useControl)
		chrMap(ChIPseq, control, mapChr);
	else
		for (i=0; i<ChIPseq->totalChr; i++)
			mapChr[i] = i;
	
	sort(ChIPseq->readsClusterVec.begin(), ChIPseq->readsClusterVec.end(), mem_fun_ref(&reads::readsCluster::sortByPos));
	calPval(ChIPseq, control, mapChr);
	pval2qval(ChIPseq);
	filteredNum = filterPeak(ChIPseq, &filteredReadCnt);
	//printf("  # Filtered: %d\t# Filtered Reads: %d\n", filteredNum, filteredReadCnt);
}

void chrMap(reads* ChIPseq, reads* control, int* mapChr)
{
	map< string, int >::iterator map_iter;
	map< string, int >::iterator map_iter_control;

	if (ChIPseq->totalChr > control->totalChr)
		fprintf(stderr, "Warning: the control data may be defective.\n");

	for (map_iter = ChIPseq->chrIdxMap.begin(); map_iter != ChIPseq->chrIdxMap.end(); map_iter ++)
	{
		mapChr[map_iter->second] = -1;
		for (map_iter_control = control->chrIdxMap.begin(); map_iter_control != control->chrIdxMap.end(); map_iter_control ++)
		{
			if (map_iter_control->first.compare(map_iter->first) == 0)
			{
				//printf("%s\t%d\t%s\t%d\n", map_iter->first.c_str(), map_iter->second, map_iter_control->first.c_str(), map_iter_control->second);
				mapChr[map_iter->second] = map_iter_control->second;
				break;
			}
		}
		if (mapChr[map_iter->second] == -1)
			fprintf(stderr, "  %s is not in control data.\n", map_iter->first.c_str());
	}
}

void calPval(reads* ChIPseq, reads* control, int* mapChr)
{
	int i;
	int totalRC = ChIPseq->readsClusterVec.size();
	vector<int>::iterator lowerBound, upperBound;
	int ChIPseqChr, controlChr;
	int controlReadCnt;
  int ReadCntR1;// ReadCntR2, ReadCntR3;
	double controlBgLambda = 0;
	double lambdaR1 = 0, lambdaR2 = 0, lambdaR3 = 0;
	double maxLambda;
	int rStart, rEnd, rMid, rLength;
	int lambdaRArm1 = param.lambdaRL1 / 2;
	int lambdaRArm2 = param.lambdaRL2 / 2;
	int lambdaRArm3 = param.lambdaRL3 / 2;
	int lambdaRStart, lambdaREnd;
	double pval;
	double foldChange;
	double *control2ChIPseqRatioFwd, *control2ChIPseqRatioRvs;
	control2ChIPseqRatioFwd = new double[ChIPseq->totalChr];
	control2ChIPseqRatioRvs = new double[ChIPseq->totalChr];
	double uniformBG = ChIPseq->readsCnt / (double) param.genomeSize;

	for (i=0; i<ChIPseq->totalChr; i++)
	{
		control2ChIPseqRatioFwd[i] = 1;
		control2ChIPseqRatioRvs[i] = 1;
	}
	if (param.useControl)
	{
		for (i=0; i<ChIPseq->totalChr; i++)
		{
			if ((controlChr = mapChr[i]) == -1)
				continue;
			control2ChIPseqRatioFwd[i] = (double) ( 1 + ChIPseq->fpos[i].size() ) / (double) ( 1 + control->fpos[controlChr].size() );
			control2ChIPseqRatioRvs[i] = (double) ( 1 + ChIPseq->rpos[i].size() ) / (double) ( 1 + control->rpos[controlChr].size() );
		}
	}
	
	for (i=0; i<totalRC; i++)
	{
		ChIPseqChr = ChIPseq->readsClusterVec[i].r_chrIdx;
		if ((controlChr = mapChr[ChIPseqChr]) == -1)
			continue;

		rStart = ChIPseq->readsClusterVec[i].r_start;
		rEnd = ChIPseq->readsClusterVec[i].r_end;
		//rLength = rEnd - rStart + 1;
    rLength = ChIPseq->readsClusterVec[i].r_thickEnd - ChIPseq->readsClusterVec[i].r_thickStart + 1 - param.smoothArm;
    //if(rLength > rEnd - rStart + 1) fprintf(stderr,"length not match\n");
		rMid = (rStart + rEnd) / 2;
	
    if ( rLength < param.bandwidth)
    {
      ChIPseq->readsClusterVec[i].pval = 1.0;
      ChIPseq->readsClusterVec[i].foldChange = 0.0;
      continue;
    }

		//printf("%s\t%s\n", ChIPseq->chrName[ChIPseqChr].c_str(), control->chrName[controlChr].c_str());

		if (ChIPseq->readsClusterVec[i].r_strand == 'f'|| ChIPseq->readsClusterVec[i].r_strand == '+')
		{
			if (param.useControl)
			{
				// for binding region
				lowerBound = lower_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), rStart - 10);
				upperBound = upper_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), rEnd + 10);
				controlReadCnt = distance(lowerBound, upperBound);
				controlBgLambda = controlReadCnt * control2ChIPseqRatioFwd[ChIPseqChr];
				//controlBgLambda = controlReadCnt * param.control2ChIPseqRatioFwd;

				// for lambda1
				lambdaRStart = rMid - lambdaRArm1;
				lambdaREnd = rMid + lambdaRArm1;
				lowerBound = lower_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaRStart);
				upperBound = upper_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaREnd);
				ReadCntR1 = distance(lowerBound, upperBound) / param.lambdaRL1 * rLength;
				lambdaR1 = ReadCntR1 * control2ChIPseqRatioFwd[ChIPseqChr];
				//lambdaR1 = ReadCntR1 * param.control2ChIPseqRatioFwd;
			}
			else
			{
				// for lambda2
				lambdaRStart = rMid - lambdaRArm2;
				lambdaREnd = rMid + lambdaRArm2;
				lowerBound = lower_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaRStart);
				upperBound = upper_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaREnd);
				controlReadCnt = distance(lowerBound, upperBound);
				lambdaR2 = controlReadCnt * control2ChIPseqRatioFwd[ChIPseqChr] / param.lambdaRL2 * rLength;
				//lambdaR2 = controlReadCnt * param.control2ChIPseqRatioFwd / param.lambdaRL2 * rLength;
				// for lambda3
				lambdaRStart = rMid - lambdaRArm3;
				lambdaREnd = rMid + lambdaRArm3;
				lowerBound = lower_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaRStart);
				upperBound = upper_bound(control->fpos[controlChr].begin(), control->fpos[controlChr].end(), lambdaREnd);
				controlReadCnt = distance(lowerBound, upperBound);
				lambdaR3 = controlReadCnt * control2ChIPseqRatioFwd[ChIPseqChr] / param.lambdaRL3 * rLength;
				//lambdaR3 = controlReadCnt * param.control2ChIPseqRatioFwd / param.lambdaRL3 * rLength;
			}
			
			//printf("%d\n", controlReadCnt);
		}
		else if (ChIPseq->readsClusterVec[i].r_strand == 'r'|| ChIPseq->readsClusterVec[i].r_strand == '-')
		{
			if (param.useControl)
			{
				// for binding region
				lowerBound = lower_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), rStart - 10);
				upperBound = upper_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), rEnd + 10);
				controlReadCnt = distance(lowerBound, upperBound);
				controlBgLambda = controlReadCnt * control2ChIPseqRatioRvs[ChIPseqChr];
				//controlBgLambda = controlReadCnt * param.control2ChIPseqRatioRvs;

				// for lambda1
				lambdaRStart = rMid - lambdaRArm1;
				lambdaREnd = rMid + lambdaRArm1;
				lowerBound = lower_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaRStart);
				upperBound = upper_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaREnd);
				ReadCntR1 = distance(lowerBound, upperBound) / param.lambdaRL1 * rLength;
				lambdaR1 = ReadCntR1 * control2ChIPseqRatioRvs[ChIPseqChr];
				//lambdaR1 = controlReadCnt * param.control2ChIPseqRatioRvs;
			}
      else
      {
			  // for lambda2
			  lambdaRStart = rMid - lambdaRArm2;
			  lambdaREnd = rMid + lambdaRArm2;
			  lowerBound = lower_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaRStart);
			  upperBound = upper_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaREnd);
			  controlReadCnt = distance(lowerBound, upperBound);
			  lambdaR2 = controlReadCnt * control2ChIPseqRatioRvs[ChIPseqChr] / param.lambdaRL2 * rLength;
			  //lambdaR2 = controlReadCnt * param.control2ChIPseqRatioRvs / param.lambdaRL2 * rLength;
			  // for lambda3
			  lambdaRStart = rMid - lambdaRArm3;
			  lambdaREnd = rMid + lambdaRArm3;
			  lowerBound = lower_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaRStart);
			  upperBound = upper_bound(control->rpos[controlChr].begin(), control->rpos[controlChr].end(), lambdaREnd);
			  controlReadCnt = distance(lowerBound, upperBound);
			  lambdaR3 = controlReadCnt * control2ChIPseqRatioRvs[ChIPseqChr] / param.lambdaRL3 * rLength;
			  //lambdaR3 = controlReadCnt * param.control2ChIPseqRatioRvs / param.lambdaRL3 * rLength;
      }

			//printf("%d\n", controlReadCnt);
		}
		else
			fprintf(stderr, "Error: the third direction: %c.\n", ChIPseq->readsClusterVec[i].r_strand);

		if (param.useControl)
    {
      maxLambda = max(uniformBG, max(controlBgLambda, lambdaR1));
			//maxLambda = max(controlBgLambda, lambdaR1);
      // for Fisher's exact test
		  pval = fisherTest(ChIPseq->readsClusterVec[i].r_count, max(controlReadCnt,ReadCntR1), ChIPseq->readsCnt, control->readsCnt);
      //printf("Fisher:\t%d\t%d\t%d\t%d\t%g\n", ChIPseq->readsClusterVec[i].r_count, max(controlReadCnt,ReadCntR1), ChIPseq->readsCnt, control->readsCnt, pval);
    }
    else
    {
      maxLambda = max(uniformBG, max(lambdaR2, lambdaR3));
      // for Poisson
		  pval = upper_pois(ChIPseq->readsClusterVec[i].r_count, maxLambda);
    }
		
    ChIPseq->readsClusterVec[i].pval = (pval > 1e-12) ? pval : 0.0;
		foldChange = (double) (ChIPseq->readsClusterVec[i].r_count + 1) / (maxLambda + 1);
		ChIPseq->readsClusterVec[i].foldChange = foldChange;

    // for fold change < 1.1
    if(foldChange < 1.1)
      ChIPseq->readsClusterVec[i].pval = 1.0;

    maxLambda = maxLambda/(double)rLength;
    if (param.smoothArm > 0)  // old flat window smoothing
    {
      int tmpSummit = ChIPseq->readsClusterVec[i].r_summit - (ChIPseq->readsClusterVec[i].r_start - param.smoothArm);
      int thickStart, thickEnd;
      if( getThickRegion(ChIPseq->readsClusterVec[i].r_vec, ChIPseq->readsClusterVec[i].r_summitVal, tmpSummit, 
            thickStart, thickEnd, maxLambda)) 
        //fprintf(stderr, "ERROR: failed in getting thick region!\n");
      {
        ChIPseq->readsClusterVec[i].r_thickStart = thickStart + ChIPseq->readsClusterVec[i].r_start - param.smoothArm;
        ChIPseq->readsClusterVec[i].r_thickEnd   = thickEnd   + ChIPseq->readsClusterVec[i].r_start - param.smoothArm;
      }
    }
    else
    {
      int tmpSummit = ChIPseq->readsClusterVec[i].r_summit - (ChIPseq->readsClusterVec[i].r_start + 3 * param.smoothBandwidth);
      int thickStart, thickEnd;
      if( getThickRegion(ChIPseq->readsClusterVec[i].r_vec, ChIPseq->readsClusterVec[i].r_summitVal, tmpSummit, 
            thickStart, thickEnd, maxLambda)) 
        //fprintf(stderr, "ERROR: failed in getting thick region!\n");
      {
        ChIPseq->readsClusterVec[i].r_thickStart = thickStart + ChIPseq->readsClusterVec[i].r_start - 3 * param.smoothBandwidth;
        ChIPseq->readsClusterVec[i].r_thickEnd   = thickEnd   + ChIPseq->readsClusterVec[i].r_start - 3 * param.smoothBandwidth;
      }
    }
		
		/*printf("%d\t%d\t%c\t%d\t%f\t%.3e\t%.3f\n", ChIPseq->readsClusterVec[i].r_start, 
			ChIPseq->readsClusterVec[i].r_end, ChIPseq->readsClusterVec[i].r_strand, ChIPseq->readsClusterVec[i].r_count, maxLambda, 
			ChIPseq->readsClusterVec[i].pval,  ChIPseq->readsClusterVec[i].foldChange);  */
	}
}

int filterPeak(reads* ChIPseq, int* filteredReadCnt)
{
	int filteredNum;
	vector<reads::readsCluster>::iterator tmpIterator;
	sort(ChIPseq->readsClusterVec.begin(), ChIPseq->readsClusterVec.end(), mem_fun_ref(&reads::readsCluster::sortByPval));
	
	//tmpIterator = binSearchPval(ChIPseq);
	tmpIterator = binSearchQval1(ChIPseq);
	//tmpIterator = binSearchQval2(ChIPseq);
	filteredNum = distance(tmpIterator, ChIPseq->readsClusterVec.end());
	ChIPseq->readsClusterVec.erase(tmpIterator, ChIPseq->readsClusterVec.end());
	/*
	for (tmpIterator = ChIPseq->readsClusterVec.begin(); tmpIterator < ChIPseq->readsClusterVec.end();)
	{
    //printf("param.bandwidth: %d\n", param.bandwidth);
		if ( (*tmpIterator).pval > param.pval_cutoff || ( (*tmpIterator).rEnd - (*tmpIterator).rStart) < param.bandwidth )
		{
      //printf("\tparam.bandwidth: %d\n", param.bandwidth);
			filteredNum ++;
			*filteredReadCnt += (*tmpIterator).r_count;
			// discard those regions
			ChIPseq->readsClusterVec.erase(tmpIterator);
		}
		else
			tmpIterator ++;
	}*/
	return filteredNum;
}

void pval2qval(reads* ChIPseq)
{
	int i, peakCnt = ChIPseq->readsClusterVec.size();
	sort(ChIPseq->readsClusterVec.begin(), ChIPseq->readsClusterVec.end(), mem_fun_ref(&reads::readsCluster::sortByPval));
	for (i=0; i<peakCnt; i++) 
	{
		ChIPseq->readsClusterVec[i].qval1 = ChIPseq->readsClusterVec[i].pval * peakCnt / (double) (i + 1) ;
		ChIPseq->readsClusterVec[i].qval2 = ChIPseq->readsClusterVec[i].pval * (double) (peakCnt - i) ;
    if (ChIPseq->readsClusterVec[i].qval1 > 1) ChIPseq->readsClusterVec[i].qval1 = 1;
    if (ChIPseq->readsClusterVec[i].qval2 > 1) ChIPseq->readsClusterVec[i].qval2 = 1;
	}
}

vector<reads::readsCluster>::iterator binSearchPval(reads *ChIPseq)
{
	int totalRC = ChIPseq->readsClusterVec.size();
	int leftB = 0, rightB = totalRC, mid = totalRC / 2;
	if (totalRC < 1)
		return ChIPseq->readsClusterVec.begin();
	vector<reads::readsCluster>::iterator tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
	while (1)
	{
		if ( tmpIterator->pval <= param.pval_cutoff)
			leftB = mid;
		else if ( tmpIterator->pval > param.pval_cutoff)
			rightB = mid;
		mid = (leftB + rightB) / 2;
		tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
		if ( tmpIterator->pval > param.pval_cutoff)
		{
			if (tmpIterator == ChIPseq->readsClusterVec.begin())
				return ChIPseq->readsClusterVec.begin();
			else if ( (tmpIterator-1)->pval <= param.pval_cutoff )
				return tmpIterator;
		}
		else  //if ( tmpIterator->pval <= param.pval_cutoff )
		{
			if ( (tmpIterator+1) == ChIPseq->readsClusterVec.end())
				return ChIPseq->readsClusterVec.end();
			else if ( (tmpIterator+1)->pval > param.pval_cutoff )
				return tmpIterator + 1;
		}
	}
}

vector<reads::readsCluster>::iterator binSearchQval1(reads *ChIPseq)
{
	int totalRC = ChIPseq->readsClusterVec.size();
	int leftB = 0, rightB = totalRC, mid = totalRC / 2;
	if (totalRC < 1)
		return ChIPseq->readsClusterVec.begin();
	vector<reads::readsCluster>::iterator tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
	while (1)
	{
		if ( tmpIterator->qval1 <= param.FDR)
			leftB = mid;
		else if ( tmpIterator->qval1 > param.FDR)
			rightB = mid;
		mid = (leftB + rightB) / 2;
		tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
		if ( tmpIterator->qval1 > param.FDR)
		{
			if (tmpIterator == ChIPseq->readsClusterVec.begin())
				return ChIPseq->readsClusterVec.begin();
			else if ( (tmpIterator-1)->qval1 <= param.FDR )
				return tmpIterator;
		}
		else  //if ( tmpIterator->qval1 <= param.FDR )
		{
			if ( (tmpIterator+1) == ChIPseq->readsClusterVec.end())
				return ChIPseq->readsClusterVec.end();
			else if ( (tmpIterator+1)->qval1 > param.FDR )
				return tmpIterator + 1;
		}
	}
}

vector<reads::readsCluster>::iterator binSearchQval2(reads *ChIPseq)
{
	int totalRC = ChIPseq->readsClusterVec.size();
	int leftB = 0, rightB = totalRC, mid = totalRC / 2;
	if (totalRC < 1)
		return ChIPseq->readsClusterVec.begin();
	vector<reads::readsCluster>::iterator tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
	while (1)
	{
		if ( tmpIterator->qval2 <= param.FDR)
			leftB = mid;
		else if ( tmpIterator->qval2 > param.FDR)
			rightB = mid;
		mid = (leftB + rightB) / 2;
		tmpIterator = ChIPseq->readsClusterVec.begin() + mid;
		if ( tmpIterator->qval2 > param.FDR)
		{
			if (tmpIterator == ChIPseq->readsClusterVec.begin())
				return ChIPseq->readsClusterVec.begin();
			else if ( (tmpIterator-1)->qval2 <= param.FDR )
				return tmpIterator;
		}
		else  //if ( tmpIterator->qval2 <= param.FDR )
		{
			if ( (tmpIterator+1) == ChIPseq->readsClusterVec.end())
				return ChIPseq->readsClusterVec.end();
			else if ( (tmpIterator+1)->qval2 > param.FDR )
				return tmpIterator + 1;
		}
	}
}
