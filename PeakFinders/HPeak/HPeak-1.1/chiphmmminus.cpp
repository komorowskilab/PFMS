// Steve Qin
// 05/28/08

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<math.h>
#include <time.h>
using namespace std;

#define MAX_LENGTH 500
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

void readData(const char dataFileName[], vector <int> &order, vector <double> &data, int *nRow);
double logPoisson(const int k, const double lambda);
void pathFinder(const int count, const vector <int> &order, const vector <double> &data, const int total, 
	int *&path, double *&proba, int *&hits, const double fgcaselambda, const double fgcontrollambda, 
	const double bglambda,const double enrichedregion, const double totalmappedreads, const double medianPeakBinCount);

double unran(int *na, int *nb, int *nc);
double gammaln(double xx);

int main(int argc, char **argv)
{
	int j,numCycles=1; 
	int count,total;
	int windowSize;
	double readcoverage;//added 04/16/08
	vector <int> order;
	vector <double> data;
	int *path;
	int *hits;
	double *proba;
	char inputname[MAX_LENGTH]="jy9.10.select.chr22.txt";//"sle.txt";
        char parasname[MAX_LENGTH]="paras.txt";
	char outputname[MAX_LENGTH]="out.txt";
	ofstream outPutFile;
	double totalmappedreads;
	double fgcaselambda, fgcontrollambda, bglambda; 
	double enrichedregion,casereads,controlreads;
        istringstream iss;
    	string lineString;   
	double temval;
	double medianPeakWidth, medianPeakBinCount;
	vector <double> datas;
	if(argc<6)
        {
          printf("5 options need to be specified:\n\tinput file name,\n\tinformation file name,\n\toutputfile name,\n\tbin size,\n\tchromosome size.\n");
            exit(0);
        }  
        for(j=0;j<MAX_LENGTH;j++)
        {
            inputname[j]=argv[1][j];
            parasname[j] = argv[2][j];
	    outputname[j] = argv[3][j];
        }       
	windowSize = atoi(argv[4]);
	total = atoi(argv[5]);
	ifstream inFile(parasname);
    	if (! inFile)
    	{
        	cout << "Error opening input parameter file" <<parasname<< endl;
        	exit(0);
    	}
	for(j=0;j<7;j++)
	{
	        getline(inFile,lineString);
        	iss.clear();
        	iss.str(lineString +" ");
        	iss >> temval ;
		datas.push_back(temval);
	}

        enrichedregion = datas[0];
	casereads = datas[1]/1000000;
	controlreads = datas[2]/1000000;
	totalmappedreads = datas[3]/1000000;
       	readcoverage = datas[5];
	medianPeakWidth = (double) datas[6];
	
	if(medianPeakWidth > 800)
		medianPeakWidth = 800;	
	medianPeakBinCount = medianPeakWidth /windowSize;
	fgcaselambda = casereads*readcoverage/enrichedregion;
        fgcontrollambda = controlreads*readcoverage/enrichedregion;
        bglambda = totalmappedreads*readcoverage/(3100*0.9);

	path = new int[total];
        proba = new double[total];   
	srand((unsigned) time(NULL));
	readData(inputname,order, data, &count);
	for(j=0;j<numCycles;j++)
	{
		pathFinder(count,order,data,total,path,proba,hits,fgcaselambda,fgcontrollambda,bglambda,enrichedregion,totalmappedreads, medianPeakBinCount);
	}
	outPutFile.open(outputname);
	if(!outPutFile)
	{
		cout << "ERROR: Unable to open file: " << outputname << endl;
	    exit(30);
	}//end of if
	for(j=0;j<total;j++)
	{
		if(proba[j] > 0.01)
		{
			outPutFile << j<<"	"<<path[j]<<"	"<<proba[j]<<"	"<<hits[j]<<endl;
		}
	}
	delete [] path;
	delete [] proba;
	delete [] hits;
	return 0;
}//end of main

void pathFinder(const int count, const vector <int> &order, const vector <double> &data, const int total,
	int *&path, double *&proba, int *&hits, const double fgcaselambda, const double fgcontrollambda,
        const double bglambda, const double enrichedregion, const double totalmappedreads, const double medianPeakBinCount)
{
	int j,k;
	double (*logfnh)[2],trp[2][2];
	double p[2],p0,p1,pp0[100],pp1[100];
	double ratio,compa;
	double dif=0,inside1,inside2,inside3,inside4;
	int na,nb,nc;
	double lambdaback, lambdafore;
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
	double sum01,sum02,tem;
	double poispfgcase[100],poispfgcontrol[100],poispbg[100];

	hits = new int[total];
	for(j=0;j<total;j++)
	{
		hits[j] = 0;
		//path[j] = 0;
	}
        for(j=0;j<count;j++)
        {                            
                if(order[j]>=total)
                {
                        cout <<order[j]<<" **"<<count<<endl;   
                        exit(0);
                }    
               hits[order[j]] = (int) data[j];
        }

	p[1] = enrichedregion/(3100*0.9);
	p[0] = 1 - p[1];
	trp[0][0] = p[0];
	trp[0][1] = p[1];
        double a = exp(1.0/medianPeakBinCount * log(0.5));
        trp[1][1] = exp(1.0/medianPeakBinCount * log(0.5)); 
        trp[1][0] = 1 - trp[1][1];             
	for(j=0;j<100;j++)
	{
		poispfgcase[j] = exp(logPoisson(j,fgcaselambda));
		poispfgcontrol[j] = exp(logPoisson(j,fgcontrollambda));
		poispbg[j] = exp(logPoisson(j,bglambda));	
	}
	for(j=0;j<100;j++)
        {
        	sum01=0;
		for(k=0;k<(100-j);k++)
                {
                        tem = poispbg[k+j] *poispbg[k];
                                sum01 = sum01 + tem;
                } 
                pp0[j] = sum01;
		sum02 = 0;
		for(k=0;k<(100-j);k++)
		{
			tem = poispfgcase[k+j] *poispfgcontrol[k];
				sum02 = sum02 + tem;
		}
		pp1[j] = sum02;
	}
	logfnh = new double[total][2];	
	logfnh[0][0] = 0;
	logfnh[0][1] = -1000;
	for(j=1;j<total;j++)
	{
		if(hits[j] ==0)
		{
			p0 = 1;
			p1 = 0;
		}
		else if(hits[j]>=100)
		{
			p0 = 0;
			p1 = 1;
		}
		else
		{
			p0 = pp0[hits[j]];
			p1 = pp1[hits[j]];
			p0 = p0/(p0+p1);
			p1 = 1 - p0;
		}
		inside3 = logfnh[j-1][1]+log(trp[1][0])-logfnh[j-1][0]-log(trp[0][0]);
		inside4 = logfnh[j-1][1]+log(trp[1][1])-logfnh[j-1][0]-log(trp[0][1]);
		if(p0 ==0)
		{
			logfnh[j][0] = -1000;
			logfnh[j][1] = 0;
		}
		else if(p1 ==0)
		{
			logfnh[j][1] = -1000;
			logfnh[j][0] = 0;
		}
		else
		{
			if(inside3 < 0)
			{
				logfnh[j][0] = log(p0)+logfnh[j-1][0]+log(trp[0][0])+log(1 + exp(inside3));
			}
			else
				logfnh[j][0] = log(p0)+logfnh[j-1][1]+log(trp[1][0])+log(1 + exp(-inside3));
			if(inside4 < 0)	
				logfnh[j][1] = log(p1)+logfnh[j-1][0]+log(trp[0][1])+log(1 + exp(inside4));
			else
				logfnh[j][1] = log(p1)+logfnh[j-1][1]+log(trp[1][1])+log(1 + exp(-inside4));
		}
	}
	ratio = 1.0/(1.0+exp(logfnh[total-1][1]-logfnh[total-1][0]));
	proba[total-1] = 1 - ratio;
	compa = unran(&na, &nb, &nc);
	if(compa <= ratio)
		path[total-1]=0;
	else
		path[total-1]=1;
	for(j = total - 2;j>=0;j--)
	{
		ratio = 1/(1 + (trp[1][path[j+1]]/trp[0][path[j+1]])*exp(logfnh[j][1]-logfnh[j][0]));
		proba[j] = 1 - ratio;
		compa = unran(&na, &nb, &nc);
		if(compa <= ratio)
			path[j]=0;
		else
			path[j]=1;
	}
}//end of pathFinder

double logPoisson(const int k, const double lambda)
{
	double x;

	x = k*log(lambda) -lambda - gammaln((double)k + 1);
	return(x);
}//end of poissonBackground

void readData(const char dataFileName[], vector <int> &order, vector <double> &data, int *nRow)
{
	int count=0;
	int temOrder;
	double temVal;
	istringstream iss;
    	string lineString;

	ifstream inFile(dataFileName);
    if (! inFile)
    {
        cout << "Error opening input file" <<dataFileName<< endl;
        exit(0);
    }
	count=0;
	while (inFile)
	{
		if(inFile)
		{
			getline(inFile,lineString);
			iss.clear();
			iss.str(lineString +" ");
			iss >> temOrder >> temVal ;
			if(iss)
			{
				order.push_back(temOrder);
				data.push_back(temVal);
			}//end of if
		}//end of if
		count ++;
	}//end of while
	*nRow=count-1;
	cout <<"There are "<<*nRow<<" nonzero counts."<<endl;
}

double unran(int *na, int *nb, int *nc)
{
	double random;
	*na=(171*(*na))%30269;
	*nb=(172*(*nb))%30307;
	*nc=(170*(*nc))%30323;
	random=(double) *na/30269.0+(double) *nb/30307.0+(double) *nc/30323.0;
	random=random-floor(random);
	return random;
}

double gammaln(double xx)
{
	double ser,stp,tmp,x,y,cof[6],gam;
	int j;
	cof[0]=76.18009172947146;
	cof[1]=-86.50532032941677;
	cof[2]=24.01409824083091;
	cof[3]=-1.231739572450155;
	cof[4]=0.1208650973866179*0.01;
	cof[5]=-0.5395239384953*0.00001;
	stp=2.5066282746310005;
	x=xx;
	y=x;
	tmp=x+5.5;
	tmp=(x+0.5)*log(tmp)-tmp;
	ser=1.000000000190015;
	for (j=0;j<6;j++) 
	{
		y=y+1.0;
		ser=ser+cof[j]/y;
	}
	gam=tmp+log(stp*ser/x);
	return gam;
}
