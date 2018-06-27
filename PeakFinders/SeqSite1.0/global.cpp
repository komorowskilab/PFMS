//#define DEBUG_G 
#include "global.h"

FILE *mustOpen(char *fileName, char *mode)	/* Open a file - or squawk and die. */                                        
{                                                                             
	FILE *f;

	if ((f = fopen(fileName, mode)) == NULL)
	{
		char *modeName = "";
		if (mode)
		{
			if (mode[0] == 'r')
				modeName = " to read";
			else if (mode[0] == 'w')
				modeName = " to write";
			else if (mode[0] == 'a')
				modeName = " to append";
		}
		fprintf(stderr, "Can't open %s%s.", fileName, modeName);
		exit(1);
	}
	return f;
}

char* readLine(FILE* fh)
/* Read a line of any size into dynamic memory, return null on EOF */

{
	int bufCapacity = 256;
	int bufSize = 0;
	char* buf = (char*) malloc(sizeof(char)*bufCapacity);
	int ch;

	/* loop until EOF of EOLN */
	while (((ch = getc(fh)) != EOF) && (ch != '\n'))
	{
		/* expand if almost full, always keep one extra char for
		 * zero termination */
		if (bufSize >= bufCapacity-2)
		{
			bufCapacity *= 2;
			buf = (char*) realloc(buf, bufCapacity);
		}
		buf[bufSize++] = ch;
	}
	/* only return EOF if no data was read */
	if ((ch == EOF) && (bufSize == 0))
	{
		free(buf);
		return NULL;
	}
	buf[bufSize] = '\0';
	return buf;
}

int bedColNum(char* line)
/* Chop the first line of bed file, return the number of cols */
{
	int colNum = 0;;
	int i = 0;
	while (line[i] != '\0')
	{
		if (line[i] != 32 && line[i] != 9)
			i ++;
		else
		{
			colNum ++;

			i ++;
			while(line[i] == 32)
				i ++;
		}
	}
	colNum ++;
	return colNum;
}

void chopLine(char* line, int colNum, char** wordArray)
/* Chop each line of bed file, return the col array*/
{
	int i, j, k;
	i = 0;
	j = 0;
	k = 0;
	while (line[i] != '\0')
	{
		if (line[i] != 32 && line[i] != 9)
		{
			wordArray[k][j] = line[i];
			i ++;
			j ++;
		}
		else
		{
			wordArray[k][j] = '\0';
			j = 0;
#ifdef DEBUG_G
			printf("%s\n", wordArray[k]);
#endif
			k ++;

			i ++;
			while(line[i] == 32)
				i ++;
		}
	}
	wordArray[k][j] = '\0';

#ifdef DEBUG_G
	printf("%s\n", wordArray[k]);
#endif
}

bool isInt(char* str)
{
	int i = 0;
	while (str[i] != '\0')
	{
		if(str[i] < 48 && str[i] > 57)
			return false;
		i ++;
	}
	return true;
}

bool isChr(char *str)
{
	if (str[0] == 'c' && str[1] == 'h' && str[2] == 'r')
		return true;
	else
		return false;
}

bool isStrand(char *str)
{
	if (str[0] == '+' || str[0] == '-')
		return true;
	else
		return false;
}

unsigned int minCntGivenFDR(int readCnt, double genomeSize, double FDR, int regionLength)
{
	int i;
	int factorial = 1;
	double lambda = readCnt / (double) genomeSize * regionLength;
	double factor1 = exp( - lambda );
	double cdf = factor1;
	for (i=1; ;i++)
	{
		factorial *= i;
		cdf += factor1 * pow(lambda, i) / factorial;
		if ((1 - cdf) < FDR)
			return i + 1;
	}
}

double upper_pois(int n, double lambda) // the cdf not including n
{
	if (n < 0)
		return 1;
	return 1 - lower_pois(n, lambda);
}

double lower_pois(int n, double lambda) // the cdf including n
{
	if (n < 0)
		return 0;

	double factor1 = exp(-lambda);
	double factor2 = 1.0;
	double sum = factor2;
	int i;
	for (i=1; i<=n; i++)
	{
		factor2 = factor2 * lambda / (double)i;
		sum += factor2;
	}
	sum *= factor1;
//	if (sum > 1)
//		sum = 1.0;
	return sum;
}

unsigned int pois_cdf_inv(double lambda, double cdf) // cdf: lower cdf
{
	if (cdf > 1 || cdf < 0)
	{
		fprintf(stderr, "cdf must be in [0,1]!\n");
		exit(1);
	}
	else if (cdf == 0)
		return 0;
	
	double factor1 = exp(-lambda);
	cdf /= factor1;

	double factor2 = 1.0;
	unsigned int i, i_max = 10;
	double sum = 1.0, sum_last;
	if (sum >= cdf)
		return 0;

	for (i=1; i<i_max; i++)
	{
		sum_last = sum;
		factor2 = factor2 * lambda / (double) i;
		sum += factor2;
		if (sum >= cdf && sum_last <= cdf)
			return i;
	}
	fprintf(stderr, "pois_cdf_inv exceeds max i 1000.\n");
	return i_max;
}

double logSum(int lower, int upper)
{
  if(lower > upper)
    return 0.0;

  int i;
  double sum = 0.0;
  for(i=lower; i<=upper; i++)
    sum += log(i);

  return sum;
}

double fisherTest(int aa, int bb, int n, int m)
{
  if(aa > n || bb > m)
    return 1.0;

  int ab = aa + bb;
  double p = 0, log_p;
  int a, b;
  for (a=aa; a<=min(n, ab); a++)
  {
    b = ab - a;
    log_p = logSum(a+1, a+b) + logSum(n-a+1, n) + logSum(m-b+1, m) - logSum(1, b) - logSum(n+m-a-b+1, n+m);
    p += exp(log_p);
  }
  return p;
}

