#ifndef __BINDSIG_H
#define __BINDSIG_H

#include <iostream>
#include <string>
#include "global.h"
#include "reads.h"
#include "seqsite.h"
#include "lm.h"
#include "parameter.h"
extern parameter param;

using namespace std;

/*
bindSig: calculate binding regions significance using w/ control data 
*/

// map the control chrId to the ChIPseq ChrId
void chrMap(reads* ChIPseq, reads* control, int* mapChr);
// filter the non-significance binding regions based on p-value calculated
void peakFilter(reads* ChIPseq, reads* control);
// calculate p-value given control data
void calPval(reads* ChIPseq, reads* control, int* mapChr);
// summary the peak filtering
int filterPeak(reads* ChIPseq, int* filteredReadCnt);
// Benjamini multiple testing correction
void pval2qval(reads* ChIPseq);
// binary search the position just larger than the p-val cutoff
vector<reads::readsCluster>::iterator binSearchPval(reads *ChIPseq);
// binary search the position just larger than the q-val_1 cutoff
vector<reads::readsCluster>::iterator binSearchQval1(reads *ChIPseq);
// binary search the position just larger than the q-val_2 cutoff
vector<reads::readsCluster>::iterator binSearchQval2(reads *ChIPseq);

#endif
