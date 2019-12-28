/*	------------------------------------------------------------------------------------

		                          * megaclust *
      unbiased hierarchical density based parallel clustering of large datasets


    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2008-2019 Nicolas Guex
    Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


	Code:       Nicolas Guex, 2008-2019
	Contact:    Nicolas.Guex@unil.ch
	Repository: https://github.com/sib-swiss/megaclust


	Articles:   megaclust was used here
                    https://www.ncbi.nlm.nih.gov/pubmed/29241546
                    https://www.ncbi.nlm.nih.gov/pubmed/23396282




	Machine :	Unix
	Language:	C
	Requires:	mpi, pthread

	Version information

	Version:	1.0  Dec.  2019 Public release of code under GPL2+ license




	Compiling:   (you will need mpi on your system)

	
	make all



	Testing:

    ./test/unit_test1.sh
    ./test/unit_test2.sh


	------------------------------------------------------------------------------------
*/

	/*------------------------- I N T E R F A C E ----------------------- */

/*

	dselect.c: process and select data to cluster and creates binary input files.

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "dclust.h"

#define kLoadingProgressReporting 500000

/* ------------------------------------------------------------------------------------ */

static unsigned int verbose = 0;

/* ------------------------------------------------------------------------------------ */

typedef	struct	__attribute__((__aligned__(16))) FACSDATA_struct	FACSDATA;
struct	FACSDATA_struct
{
	CELLNAMEIDX cellnameidx;
	unsigned short data[kMaxInputCol];
};


typedef	struct	FACSNAME_struct	FACSNAME;
struct	FACSNAME_struct
{
	char name[kMaxCellName];	
	unsigned  int selcnt;	
	unsigned  int leftovercnt;	
};


typedef	struct LEFTOVER_struct	LEFTOVER;
#pragma pack(1)
struct	LEFTOVER_struct
{
	CELLNAMEIDX cellnameidx;
	float val[kMaxInputCol];
};

/* ------------------------------------------------------------------------------------ */

static unsigned int ReadUnAssignedFileHeader(FILE *uf,unsigned int *rowcnt,unsigned int *colcnt)
{
	unsigned int ui;
	unsigned int maxclusterid;
	int endian;
	char hdr[kHeaderSize];

	fread(&hdr[0],sizeof(char),kHeaderSize,uf);
	if (strcmp(hdr,"dclust unassigned file v1.0   \n") != 0)
	{
		printf("Error: input is not a dclust unassigned file\n");
		return(0);
	}	

	fread(&endian,sizeof(int),1,uf);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(0);
	}
	fread(&ui,sizeof(int),1,uf);
	*rowcnt = ui;
		
	fread(&ui,sizeof(int),1,uf);
	*colcnt = ui;
	
	fread(&maxclusterid,sizeof(int),1,uf);
	return(maxclusterid);
	
} /* ReadUnAssignedFileHeader */

/* ------------------------------------------------------------------------------------ */
static unsigned int ReadAssignedFileHeader(FILE *uf,unsigned int *rowcnt,unsigned int *colcnt)
{
	unsigned int ui;
	unsigned int maxclusterid;
	int endian;
	char hdr[kHeaderSize];

	fread(&hdr[0],sizeof(char),kHeaderSize,uf);
	if (strcmp(hdr,"dclust assigned file v1.0     \n") != 0)
	{
		printf("Error: input is not a dclust assigned file\n");
		return(0);
	}	

	fread(&endian,sizeof(int),1,uf);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(0);
	}
	fread(&ui,sizeof(int),1,uf);
	*rowcnt = ui;
		
	fread(&ui,sizeof(int),1,uf);
	*colcnt = ui;
	
	fread(&maxclusterid,sizeof(int),1,uf);
	return(maxclusterid);
	
} /* ReadAssignedFileHeader */

/* ------------------------------------------------------------------------------------ */

static unsigned short processUnAssignedFile(FILE *f, FILE *af, FACSDATA *facsdata,FACSNAME *uniquecellnames,unsigned short cellnamecnt, unsigned int inputrcnt,unsigned int *selectedcnt,unsigned int colcnt,unsigned short *key)
{
	FACSDATA *facsp;
	unsigned int i;
	char linbuf[kMaxLineBuf];
	int endian;
	unsigned int rcnt;
	char hdr[kHeaderSize];
	long long sum[kMaxInputCol];
	double score[kMaxInputCol];
	double bestscore;
	unsigned int 	selcnt = 0;
	unsigned short cn;


		

	strcpy(hdr,"dclust input file v1.0        \n");
	fwrite(&hdr[0],sizeof(char),kHeaderSize,af);

	endian = 1;
	fwrite(&endian,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af); // reserve space for rowcnt
	fwrite(&colcnt,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af);//place for loadEveryNsample, which is  = 1

	/* read header */
	fseek(f,(kHeaderSize+4*sizeof(int)),SEEK_SET);
	fread(&linbuf[0],kMaxLineBuf,1,f);

	/* write header */
	fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,af);

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
		sum[cn] = 0;

	for (cn = 0; cn<cellnamecnt;cn++)
		uniquecellnames[cn].selcnt = 0;

	facsp = &facsdata[0];
	for (rcnt = 0; rcnt< inputrcnt ; rcnt++)
	{
		float	floatdata[kMaxInputCol];

		fread(&facsp->cellnameidx,sizeof(CELLNAMEIDX),1,f);
		fread(&floatdata[0],sizeof(float),colcnt,f);

		for (i=0; i<colcnt; i++)
		{
			facsp->data[i] = (unsigned short)floatdata[i];
			sum[i]+=(long long)facsp->data[i];
		}
		facsp++;
		selcnt++;

	}
	fseek(af,36L,SEEK_SET);
	fwrite(&selcnt,sizeof(int),1,af);

	printf("LOG: %10d columns in input file\n",colcnt);
	printf("LOG: %10d events in input file\n",inputrcnt);
	printf("LOG: %10d events selected\n",selcnt);

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
	{
		/* compute SDDEV for a column */
		int mean = (int)(sum[cn] / selcnt);
		long long sd = 0;
		facsp = &facsdata[0];
		for (i = 0; i<selcnt;i++)
		{
			int diff = (int)facsp->data[cn] - mean;
			sd += (long long)(diff * diff); 
			facsp++;
		}
		score[cn] = (double)sd;
		if (verbose > 0)
			printf("LOG: column %4d: mean=%5d stdev=%.0f\n",cn,mean,sqrt(score[cn]/selcnt));
	}
	/* select column with largest sddev as sorting key */
	bestscore = 0.0;
	for (cn = 0; cn<colcnt;cn++)
	{
		if (score[cn] > bestscore)
		{
			bestscore = score[cn];
			*key = cn;
		}
	}
	if (verbose > 1)
		printf("LOG: colkey = %d\n",*key);	
	
	*selectedcnt = selcnt;
	return(cellnamecnt);
	
} /* processUnAssignedFile */

/* ------------------------------------------------------------------------------------ */

static unsigned short processAssignedFile(FILE *f, FILE *af, FACSDATA *facsdata,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int cluster, unsigned int inputrcnt,unsigned int *selectedcnt,unsigned int colcnt,unsigned short *key)
{
	FACSDATA *facsp;
	unsigned int i;
	char linbuf[kMaxLineBuf];
	int endian;
	unsigned int rcnt;
	char hdr[kHeaderSize];
	long long sum[kMaxInputCol];
	double score[kMaxInputCol];
	double bestscore;
	unsigned int 	selcnt = 0;
	unsigned short cn;


		

	strcpy(hdr,"dclust input file v1.0        \n");
	fwrite(&hdr[0],sizeof(char),kHeaderSize,af);

	endian = 1;
	fwrite(&endian,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af); // reserve space for rowcnt
	fwrite(&colcnt,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af);//place for loadEveryNsample, which is  = 1

	/* read header */
	fseek(f,(kHeaderSize+4*sizeof(int)),SEEK_SET);
	fread(&linbuf[0],kMaxLineBuf,1,f);
	/* remove the ,cluster from header */
	i=strlen(linbuf);
	linbuf[i-8] = 0;  /* remove the ,cluster */				
	linbuf[i-7] = 0;  /* remove the ,cluster */				
	linbuf[i-6] = 0;  /* remove the ,cluster */				
	linbuf[i-5] = 0;  /* remove the ,cluster */				
	linbuf[i-4] = 0;  /* remove the ,cluster */				
	linbuf[i-3] = 0;  /* remove the ,cluster */				
	linbuf[i-2] = 0;  /* remove the ,cluster */				
	linbuf[i-1] = 0;  /* remove the ,cluster */				

	/* write header */
	fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,af);


	for (cn = 0; cn<(unsigned short)colcnt;cn++)
		sum[cn] = 0;

	for (cn = 0; cn<cellnamecnt;cn++)
		uniquecellnames[cn].selcnt = 0;

	facsp = &facsdata[0];
	for (rcnt = 0; rcnt< inputrcnt ; rcnt++)
	{
		float	floatdata[kMaxInputCol];
		unsigned int clusterid;

		fread(&facsp->cellnameidx,sizeof(CELLNAMEIDX),1,f);
		fread(&floatdata[0],sizeof(float),colcnt,f);
		fread(&clusterid,sizeof(int),1,f);

		if (clusterid == cluster)
		{
			for (i=0; i<colcnt; i++)
			{
				facsp->data[i] = (unsigned short)floatdata[i];
				sum[i]+=(long long)facsp->data[i];
			}
			facsp++;
			selcnt++;
		}

	}
	fseek(af,36L,SEEK_SET);
	fwrite(&selcnt,sizeof(int),1,af);

	printf("LOG: %10d columns in input file\n",colcnt);
	printf("LOG: %10d events in input file\n",inputrcnt);
	printf("LOG: %10d events selected\n",selcnt);

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
	{
		/* compute SDDEV for a column */
		int mean = (int)(sum[cn] / selcnt);
		long long sd = 0;
		facsp = &facsdata[0];
		for (i = 0; i<selcnt;i++)
		{
			int diff = (int)facsp->data[cn] - mean;
			sd += (long long)(diff * diff); 
			facsp++;
		}
		score[cn] = (double)sd;
		if (verbose > 0)
			printf("LOG: column %4d: mean=%5d stdev=%f\n",cn,mean,sqrt(score[cn]/selcnt));
	}
	/* select column with largest sddev as sorting key */
	bestscore = 0.0;
	for (cn = 0; cn<colcnt;cn++)
	{
		if (score[cn] > bestscore)
		{
			bestscore = score[cn];
			*key = cn;
		}
	}
	if (verbose > 1)
		printf("LOG: colkey = %d\n",*key);	
	
	*selectedcnt = selcnt;
	return(cellnamecnt);
	
} /* processAssignedFile */

/* ------------------------------------------------------------------------------------ */

static unsigned short processInputFile(FILE *f, FILE *af, FILE *lf,unsigned int loadEveryNsample, FACSDATA *facsdata,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int *selectedcnt,unsigned int *columns,unsigned short *key,unsigned int firstColIsSelectFlag,unsigned  short *minkeyval,unsigned short *maxkeyval)
{
	FACSDATA *facsp;
	unsigned short flt10000[64];
	unsigned short flt1000[64];
	unsigned short flt100[64];
	unsigned short flt10[64];
	unsigned int i;
	unsigned int colcnt;
	char linbuf[kMaxLineBuf];
	int endian;
	unsigned int rowcnt = 0;
	unsigned int leftovercnt = 0;
	char hdr[kHeaderSize];
	long long sum[kMaxInputCol];
	double score[kMaxInputCol];
	double  bestscore;
	unsigned int skip;
	unsigned short cn;
	unsigned int nextprintout = kLoadingProgressReporting;
	unsigned int leftoverstructsize;
	LEFTOVER lo;
	unsigned short minval = 65535;
	unsigned short maxval = 0;

	for (i=0; i<64; i++)
	{
		flt10000[i] = 0;
		flt1000[i] = 0;
		flt100[i] = 0;
		flt10[i] = 0;
	}
	for (i=49; i<=57; i++)  /* from '1' to '9' */
	{
		flt10000[i] = 10000*(i-48);
		flt1000[i] = 1000*(i-48);
		flt100[i] =  100*(i-48);
		flt10[i] =  10*(i-48);
	}
		
	facsp = &facsdata[0];
	
	if (firstColIsSelectFlag)
	{
		char c;
		/* suppress first column header */
		while((c = getc(f)) != ',') {}
	}
	fgets(linbuf,kMaxLineBuf,f);	/* retrieve number of columns from  header */
	colcnt = 0; /* start at 0 to return effective number of columns containing data, thus skipping sample name  */
	i = 0;
	do 
	{
		if (linbuf[i] == ',') colcnt++;
		if ((linbuf[i] == '\n') || (linbuf[i] == '\r')) linbuf[i] = 0;;
	} while (linbuf[i++] != 0);
	while (i++ < kMaxLineBuf) { linbuf[i] = 0; }

	strcpy(hdr,"dclust input file v1.0        \n");
	fwrite(&hdr[0],sizeof(char),kHeaderSize,af);

	endian = 1;
	fwrite(&endian,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af); // reserve space for rowcnt
	fwrite(&colcnt,sizeof(int),1,af);
	fwrite(&loadEveryNsample,sizeof(int),1,af);

	/* write header */
	fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,af);

	if (lf)
	{
		strcpy(hdr,"dclust unassigned file v1.0   \n");
		fwrite(&hdr[0],sizeof(char),kHeaderSize,lf);

		endian = 1;
		fwrite(&endian,sizeof(int),1,lf);
		fwrite(&endian,sizeof(int),1,lf); // reserve space for rowcnt
		fwrite(&colcnt,sizeof(int),1,lf);
		/* write spacer */
		endian = 0;
		fwrite(&endian,sizeof(int),1,lf);
		/* write header */
		fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,lf);
	}

	skip = 0;
	for (cn = 0; cn<(unsigned short)colcnt;cn++)
		sum[cn] = 0;

	leftoverstructsize = sizeof(CELLNAMEIDX) + colcnt*sizeof(float);
	cn = 0;
	do
	{
		int l;
		unsigned int canselect;

		/* read one line */
		fgets(linbuf,kMaxLineBuf,f);
		if (feof(f))
			break;

		canselect = 1;
		if (firstColIsSelectFlag)
			sscanf(&linbuf[0],"%u,%u,%n",&canselect,&lo.cellnameidx,&l);
		else
			sscanf(&linbuf[0],"%u,%n",&lo.cellnameidx,&l);

		if ((skip != 0) || (canselect == 0)) /* put in leftover */
		{
			unsigned short c1,c2,c3,c4,c5;
			char *cp = &linbuf[l];
			
			for (i=0; i<(colcnt-1); i++)
			{
				c1 = *cp++;
				c2 = *cp++;
				if (c2 == ',')
				{
					lo.val[i] = (float)(c1-48);
					continue;
				}
				c3 = *cp++;
				if (c3 == ',')
				{
					lo.val[i]=(float)(flt10[c1]+c2-48);
					continue;
				}
				c4 = *cp++;
				if (c4 == ',')
				{
					lo.val[i]=(float)(flt100[c1]+flt10[c2]+c3-48);
					continue;
				}
				c5 = *cp++;
				if (c5 == ',')
				{
					lo.val[i]=(float)(flt1000[c1]+flt100[c2]+flt10[c3]+c4-48);
					continue;
				}
				cp++;
				lo.val[i]=(float)(flt10000[c1]+flt1000[c2]+flt100[c3]+flt10[c4]+c5-48);
			}
			/* read last col */
			c1 = *cp++;
			c2 = *cp++;
			if (c2 == '\n')
			{
				lo.val[i] = (float)(c1-48);
				goto skipit;
			}
			c3 = *cp++;
			if (c3 == '\n')
			{
				lo.val[i]=(float)(flt10[c1]+c2-48);
				goto skipit;
			}
			c4 = *cp++;
			if (c4 == '\n')
			{
				lo.val[i]=(float)(flt100[c1]+flt10[c2]+c3-48);
				goto skipit;
			}
			c5 = *cp++;
			if (c5 == ',')
			{
				lo.val[i]=(float)(flt1000[c1]+flt100[c2]+flt10[c3]+c4-48);
				continue;
			}
			cp++;
			lo.val[i]=(float)(flt10000[c1]+flt1000[c2]+flt100[c3]+flt10[c4]+c5-48);
		skipit:
		
		
			fwrite(&lo.cellnameidx,leftoverstructsize,1,lf);
			leftovercnt++;
		}
		else /* retain */
		{
			unsigned short tot = l;
			facsp->cellnameidx = lo.cellnameidx;

			for (i=0; i<colcnt; i++)
			{
				sscanf(&linbuf[tot],"%hu,%n",&facsp->data[i],&l); tot+=l;
				sum[i]+=(long long)facsp->data[i];
			}
			facsp++;
			rowcnt++;
		}
		if (canselect) /* otherwise was a no select and does not count */
			skip++;
		if (skip == loadEveryNsample)
		{
			skip = 0;
			if (rowcnt > nextprintout)
			{
				printf("LOG: %10d events selected so far...\n",nextprintout);
				fflush(stdout);
				nextprintout += kLoadingProgressReporting;
			}
		}

	} while(1);
	fseek(af,36L,SEEK_SET);
	fwrite(&rowcnt,sizeof(int),1,af);

	if (lf)
	{
		fseek(lf,36L,SEEK_SET);
		fwrite(&leftovercnt,sizeof(int),1,lf);
	}
	printf("LOG: %10d columns in input file\n",colcnt);
	printf("LOG: %10d events in input file\n",rowcnt+leftovercnt);
	printf("LOG: %10d events selected\n",rowcnt);
	printf("LOG: %10d events leftover\n",leftovercnt);

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
	{
		/* compute SDDEV for a column */
		int mean = (int)(sum[cn] / rowcnt);
		long long sd = 0;
		facsp = &facsdata[0];
		for (i = 0; i<rowcnt;i++)
		{
			int diff = (int)facsp->data[cn] - mean;
			sd += (long long)(diff * diff);
			facsp++;
		}
		score[cn] = (double)sd;
		if (verbose > 0)
			printf("LOG: column %4d: mean=%5d stdev=%f\n",cn,mean,sqrt(score[cn]/rowcnt));
	}
	/* select column with largest sddev as sorting key */
	bestscore = 0.0;
	for (cn = 0; cn<colcnt;cn++)
	{
		if (score[cn] > bestscore)
		{
			bestscore = score[cn];
			*key = cn;
		}
	}
	facsp = &facsdata[0];
	cn = *key;
	for (i = 0; i<rowcnt;i++)
	{
		if (facsp->data[cn] > maxval)
			maxval = facsp->data[cn]; 
		if (facsp->data[cn] < minval)
			minval = facsp->data[cn]; 
		facsp++;
	}
	if (verbose > 1)
		printf("LOG: colkey = %d; (%u - %u)\n",*key,minval,maxval);	
	
	*minkeyval = minval;
	*maxkeyval = maxval;
	*selectedcnt = rowcnt;
	*columns = colcnt;
	
	cellnamecnt = 1;
	
	return(cellnamecnt);
	
} /* processInputFile */

/* ------------------------------------------------------------------------------------ */
static unsigned short SafeProcessInputFile(FILE *f, FILE *af, FILE *lf,unsigned int loadEveryNsample, FACSDATA *facsdata,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int *selectedcnt,unsigned int *columns,unsigned short *key, unsigned int firstColIsSelectFlag,unsigned  short *minkeyval,unsigned short *maxkeyval)
{
	FACSDATA *facsp;
	unsigned int i;
	unsigned int colcnt;
	char linbuf[kMaxLineBuf];
	float val[kMaxInputCol];
	int endian;
	unsigned int rowcnt = 0;
	unsigned int leftovercnt = 0;
	char hdr[kHeaderSize];
	long long sum[kMaxInputCol];
	double score[kMaxInputCol];
	unsigned int skip;
	CELLNAMEIDX cn;
	double  bestscore;
	unsigned int nextprintout = kLoadingProgressReporting;
	unsigned short minval = 65535;
	unsigned short maxval = 0;
	int maxInputVal = *maxkeyval;
		
	facsp = &facsdata[0];
	
	if (firstColIsSelectFlag)
	{
		char c;
		/* suppress first column header */
		while((c = getc(f)) != ',') {}
	}
	fgets(linbuf,kMaxLineBuf,f);	/* retrieve number of columns from  header (last colun is clusterid) */
	colcnt = 0; /* start at 0 to return effective number of columns containing data, thus skipping sample name  */
	i = 0;
	do 
	{
		if (linbuf[i] == ',') colcnt++;
		if ((linbuf[i] == '\n') || (linbuf[i] == '\r')) linbuf[i] = 0;;
	} while (linbuf[i++] != 0);
	while (i++ < kMaxLineBuf) { linbuf[i] = 0; }

	if (colcnt > kMaxInputCol)
	{
		printf("LOG:Fatal: number of requested columns exceed maximum allowed (%u > %d)\n",colcnt,kMaxInputCol);
		return(0);
	}	
	
	strcpy(hdr,"dclust input file v1.0        \n");
	fwrite(&hdr[0],sizeof(char),kHeaderSize,af);

	endian = 1;
	fwrite(&endian,sizeof(int),1,af);
	fwrite(&endian,sizeof(int),1,af); // reserve space for rowcnt
	fwrite(&colcnt,sizeof(int),1,af);
	fwrite(&loadEveryNsample,sizeof(int),1,af);

	/* write header */
	fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,af);

	if (lf)
	{
		strcpy(hdr,"dclust unassigned file v1.0   \n");
		fwrite(&hdr[0],sizeof(char),kHeaderSize,lf);

		endian = 1;
		fwrite(&endian,sizeof(int),1,lf);
		fwrite(&endian,sizeof(int),1,lf); // reserve space for rowcnt
		fwrite(&colcnt,sizeof(int),1,lf);
		/* write spacer */
		endian = 0;
		fwrite(&endian,sizeof(int),1,lf);
		/* write header */
		fwrite(&linbuf[0],sizeof(char),kMaxLineBuf,lf);
	}
	skip = 0;

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
		sum[cn] = 0;

	cn = 0;
	do
	{
		int tot,l;
		unsigned int canselect;
		int inputVal;
		
		/* read one line */
		fgets(linbuf,kMaxLineBuf,f);
		if (feof(f))
			break;

		tot = 0;

		canselect = 1;
		if (firstColIsSelectFlag)
			sscanf(&linbuf[0],"%u,%u,%n",&canselect,&cn,&l);
		else
			sscanf(&linbuf[0],"%u,%n",&cn,&l);
		tot+=l;

		if ((skip != 0) || (canselect == 0)) /* put in leftover */
		{
			for (i=0; i<colcnt; i++)
			{
				unsigned int tmp;
				sscanf(&linbuf[tot],"%u,%n",&tmp,&l); tot+=l;		
				val[i] = (float)tmp;
			}

			fwrite(&cn,sizeof(CELLNAMEIDX),1,lf);
			fwrite(&val[0],sizeof(float),colcnt,lf);
			leftovercnt++;
		}
		else /* retain */
		{
			facsp->cellnameidx = cn;
			for (i=0; i<colcnt; i++)
			{
				sscanf(&linbuf[tot],"%d,%n",&inputVal,&l); tot+=l;
				if ((inputVal < 0) || inputVal > maxInputVal)
				{
					printf("LOG:Fatal: input value %d out of valid supported unsigned range [0..%d]. Please rescale your data first\n",inputVal,maxInputVal);
					return(0);
				}
				facsp->data[i] = (unsigned short)inputVal;
				sum[i]+=(long long)facsp->data[i];
			}
			facsp++;
			rowcnt++;
		}
		if (canselect) /* otherwise was a no select and does not count */
			skip++;
		if (skip == loadEveryNsample)
		{
			skip = 0;
			if (rowcnt > nextprintout)
			{
				printf("LOG: %10d events selected so far...\n",nextprintout);
				fflush(stdout);
				nextprintout += kLoadingProgressReporting;
			}
		}

	} while(1);
	fseek(af,36L,SEEK_SET);
	fwrite(&rowcnt,sizeof(int),1,af);

	if (lf)
	{
		fseek(lf,36L,SEEK_SET);
		fwrite(&leftovercnt,sizeof(int),1,lf);
	}
	printf("LOG: %10d columns in input file\n",colcnt);
	printf("LOG: %10d events in input file\n",rowcnt+leftovercnt);
	printf("LOG: %10d events selected\n",rowcnt);
	printf("LOG: %10d events leftover\n",leftovercnt);

	for (cn = 0; cn<(unsigned short)colcnt;cn++)
	{
		/* compute SDDEV for a column */
		int mean = (int)(sum[cn] / rowcnt);
		long long sd = 0;
		facsp = &facsdata[0];
		for (i = 0; i<rowcnt;i++)
		{
			int diff = (int)facsp->data[cn] - mean;
			sd += (long long)(diff * diff); 
			facsp++;
		}
		score[cn] = (double)sd;
		if (verbose > 0)
			printf("LOG: column %4d: mean=%5d stdev=%f\n",cn,mean,sqrt(score[cn]/rowcnt));
	}
	/* select column with largest sddev as sorting key */
	bestscore = 0.0;
	for (cn = 0; cn<colcnt;cn++)
	{
		if (score[cn] > bestscore)
		{
			bestscore = score[cn];
			*key = cn;
		}
	}

	facsp = &facsdata[0];
	cn = *key;
	for (i = 0; i<rowcnt;i++)
	{
		if (facsp->data[cn] > maxval)
			maxval = facsp->data[cn]; 
		if (facsp->data[cn] < minval)
			minval = facsp->data[cn]; 
		facsp++;
	}
	if (verbose > 1)
		printf("LOG: colkey = %d; (%u - %u)\n",*key,minval,maxval);	
	
	*selectedcnt = rowcnt;
	*columns = colcnt;
	*minkeyval = minval;
	*maxkeyval = maxval;

	cellnamecnt = 1;

	return(cellnamecnt);
	
} /* SafeProcessInputFile */

/* ------------------------------------------------------------------------------------ */
static int WriteUniqueCellNames(FILE *af,FACSNAME *facsname,FACSDATA	*facsdata,unsigned short cellnamecnt,unsigned int rcnt,unsigned int colcnt,unsigned short key,unsigned short minkeyval,unsigned short maxkeyval)
{
	unsigned int i;
	unsigned short val;
	FACSDATA	*facsp;
	unsigned int cnt[65536];
	unsigned int *sortedcellnameidx=NULL;
	char emptyStr[kMaxCellName]; 
	strcpy(emptyStr,"!NO_CATEGORIES!");

		/* skip header */
		fseek(af,(kHeaderSize+4*sizeof(int)+kMaxLineBuf),SEEK_SET);

		/* write column key */
		fwrite(&key,sizeof(unsigned short),1,af);
		/* write unique cellnames count */
		fwrite(&cellnamecnt,sizeof(unsigned short),1,af);
		/* add unique cellnames */		
		fwrite(&emptyStr,sizeof(char),kMaxCellName,af);
	

		
		/* write data sorted according to key */
		printf("LOG: Writing Selected Events\n");

	
		sortedcellnameidx = malloc(rcnt*sizeof(unsigned int));
		if (sortedcellnameidx)
		{
			unsigned int start;
			unsigned int last;
			unsigned short val;
			unsigned int key2;
			unsigned int *sortedcellnameidxcopy=NULL;
			for (i = minkeyval; i <= maxkeyval; i++)
				cnt[i] = 0;

			facsp = &facsdata[0];
			for (i = 0; i<rcnt;i++)
			{
				cnt[facsp->data[key]]++;
				facsp++;
			}

			for (i = minkeyval+1; i <= maxkeyval; i++)
				cnt[i] += cnt[i-1];
			for (i = minkeyval; i <= maxkeyval; i++)
				cnt[i]--;

			for (i = rcnt-1; i>=1;i--)
			{
				facsp--;
				sortedcellnameidx[cnt[facsp->data[key]]--] = i;
			}
			facsp--;
			sortedcellnameidx[cnt[facsp->data[key]]--] = 0;


			/* sort by second key */
			sortedcellnameidxcopy = malloc(rcnt*sizeof(unsigned int));
			if (sortedcellnameidxcopy)
			{
				/* make a copy of the sorted by first key, we will read from here and modify the original as we go */
				memcpy(sortedcellnameidxcopy,sortedcellnameidx,rcnt*sizeof(unsigned int));
				if (key > 0)
					key2 = key-1;
				else
					key2 = 1;
				start = 0;
				do
				{
					
					for (i = 0; i <= 65535; i++)
						cnt[i] = 0;
					val = facsdata[sortedcellnameidxcopy[start]].data[key];
					minkeyval=65535;
					maxkeyval=0;
					last = start;
					while (facsdata[sortedcellnameidxcopy[last]].data[key] == val)
					{
						unsigned short val=facsdata[sortedcellnameidxcopy[last++]].data[key2];
						if (val > maxkeyval)
							maxkeyval=val;
						if (val < minkeyval)
							minkeyval=val;
						cnt[val]++;
						if (last >= rcnt)
							break;
					};
					if (maxkeyval > minkeyval)
					{
						unsigned short val;
						
						for (i = minkeyval+1; i <= maxkeyval; i++)
							cnt[i] += cnt[i-1];
						for (i = minkeyval; i <= maxkeyval; i++)
							cnt[i]--;


						for (i = last-1; i>=(start+1);i--)
						{
							val=facsdata[sortedcellnameidxcopy[i]].data[key2];
							sortedcellnameidx[start+ cnt[val] ] = sortedcellnameidxcopy[i] ;
							cnt[val]--;
						}
						val=facsdata[sortedcellnameidxcopy[start]].data[key2];
						sortedcellnameidx[start + cnt[val]] = sortedcellnameidxcopy[start];
					}
					start=last;
				} while(last < rcnt);
				free(sortedcellnameidxcopy);
			}
			

			/* write results */
			for (i = 0; i<rcnt;i++)
			{
				facsp = &facsdata[sortedcellnameidx[i]];
				fwrite(&facsp->cellnameidx,sizeof(CELLNAMEIDX),1,af);
				fwrite(&facsp->data,sizeof(unsigned short),colcnt,af);
			}
			free(sortedcellnameidx);
		}
		else	/* resort to slower sort without memory overhead */
		{
			for (val = minkeyval; val<= maxkeyval;val++)
			{
				facsp = &facsdata[0];
				for (i = 0; i<rcnt;i++)
				{
					if (facsp->data[key] == val)
					{
						fwrite(&facsp->cellnameidx,sizeof(CELLNAMEIDX),1,af);
						fwrite(&facsp->data,sizeof(unsigned short),colcnt,af);
					}
					facsp++;
				}
			}
		}
		return(0);
		
} /* WriteUniqueCellNames */
/* ------------------------------------------------------------------------------------ */

int main (int argc, char **argv)
{	
	FACSDATA	*facsdata;
	char	*version="VERSION 1.0; 2019-12-26";
	char ifn[kMaxFilename];
	char wfn[kMaxFilename];
	char ofn[kMaxFilename];
	char nfn[kMaxFilename];
	unsigned int quickprocess = 0;
	unsigned int totalrowcnt = 0;
	unsigned int colcnt = 0;
	FILE *f;
	FILE *wf;
	FILE *lf = NULL;
	FACSNAME *facsname = NULL;
	int c;
	int err;
	int loadEveryNsample;
	unsigned int lastclusterid = 0;
	unsigned int firstColIsSelectFlag = 0;
	unsigned short cellnamecnt = 0;
	unsigned int binary = 0;
	unsigned int readFromUnassigned = 0;
	int  keyOverride = -1;

	/* --------- process arguments */

	loadEveryNsample = 1;

	ifn[0] = 0;
	wfn[0] = 0;
	ofn[0] = 0;
	nfn[0] = 0;

	opterr = 0;
	while ((c = getopt (argc, argv, "i:b:o:s:qfv:k:u:")) != -1)
	switch (c)
	{      
	  case 'i':
			strcpy(ifn,optarg);
        break;

	  case 'b':
			binary = 1;
			readFromUnassigned = 0;
			strcpy(ifn,optarg);
        break;

	  case 'u':
			binary = 0;
			readFromUnassigned = 1;
			strcpy(ifn,optarg);
        break;

	  case 's':
			sscanf(optarg,"%d",&loadEveryNsample);
        break;

	  case 'o':
			strcpy(ofn,optarg);
        break;

	  case 'q':
			quickprocess = 1;
        break;

	  case 'f':
			firstColIsSelectFlag = 1;
        break;

	  case 'v':
			sscanf(optarg,"%d",&verbose);
        break;

	  case 'k':
			sscanf(optarg,"%d",&keyOverride);
        break;
	}

	if (ofn[0] == 0)
	{
			strcpy(ofn,ifn);
	}	
	if (ifn[0] == 0)
	{
		printf("usage:\n\n");
		printf("dselect -i|-b|-u InputFile [-o OutputFileRootName] [-s LoadEveryNsample][-q][-f][-n NamesFile][-v level]\n\n");
		printf("        -i InputFile           : Comma delimited file. A header is expected.\n");
		printf("                                 Values must be integers in the [0..1023] range.\n");
		printf("        -b InputFile           : should be a dclust .assigned output file.\n");
		printf("                                 It will produce one dclust binary input file 'OutputFileRootName.clusterID.selected'\n");
		printf("                                 for each clusterID present in the file specified. Options -s -q and -f are ignored but -n is mandatory.\n");
		printf("                                 The sole purpose of this option is to make input file for the splitclust.sh process\n");
		printf("        -u InputFile           : should be a dclust .unassigned output file.\n");
		printf("                                 It will produce one dclust binary input file 'OutputFileRootName.selected'\n");
		printf("                                 Options -s -q and -f are ignored but -n is mandatory.\n");
		printf("                                 The sole purpose of this option is to allow further clustering on unassigned events\n");
		printf("        -s LoadEveryNsample    : Select every 1/LoadEveryNsample sample.\n");
		printf("                                 A value of 2 will load one sample, skip the next one and so on.\n");
		printf("                                 Default value is 1 (all events selected).\n");
		printf("        -o OutputFileRootName  : Selected events will be written in binary format suitable for dclust\n");
		printf("                                 in the file 'OutputFileRootName.selected'\n");
		printf("                                 A second binary file 'OutputFileRootName.leftover' is also produced\n");
		printf("                                 The default OutputFileRootName is equal to InputFile\n");
		printf("        -q                     : quick parsing, no check is done, input values must be comma separated without spaces\n");
		printf("        -f                     : first column of the InputFile has to be treated as a selection flag\n");
		printf("                               : a content of 0 will place the event in the leftover bin (transparent for option -s)\n");
		printf("                               : a content of 1 will select or leave the event depending on option -s\n");
		printf("      -v level                 : specifies the verbose level; default is 0.\n\n");
		printf("VERSION\n");
		printf("\n%s\n",version);
		printf("Author:  Nicolas Guex; 2008-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
		printf("CONTACT: Nicolas.Guex@unil.ch\n");
		printf("SEE ALSO\n");
		printf("      dclust cextract\n\n");
		return(1);
	}


	sprintf(wfn,"%s.selected",ofn);
	/* --------- process */
	if (strcmp(ifn,wfn) == 0)
	{
		printf("Error:Input Filename must differ from Output Filename\n");
		return(1);
	}
	
	if (binary || readFromUnassigned)
	{
		if (nfn[0]==0)
		{
			printf("Error:option -b and -u require use of option -n\n");
			return(1);
		}
		f = fopen(ifn,"rb");
	}
	else
		f = fopen(ifn,"r");
		
	if (!f)
	{
		printf("Error:Cannot read Input File %s\n",ifn);
		return(1);
	}
	

	if (binary) /* retrieve actual count of events from input file */
	{
		loadEveryNsample = 1;
		firstColIsSelectFlag = 0;
		lastclusterid = ReadAssignedFileHeader(f,&totalrowcnt,&colcnt);
		facsdata = calloc(totalrowcnt,sizeof(FACSDATA));
	}
	else 	/* allocate memory for input data */
	{
		if (readFromUnassigned)
		{
			loadEveryNsample = 1;
			ReadUnAssignedFileHeader(f,&totalrowcnt,&colcnt);
			facsdata = calloc(totalrowcnt,sizeof(FACSDATA));
		}
		else
			facsdata = calloc(kMAXEVENTS,sizeof(FACSDATA));
	}
	if (!facsdata)
	{
		printf("Error:Cannot Allocate Memory to select up to %d events.\n",kMAXEVENTS);
		fclose(f);
		return(1);
	}



	err = 0;		
	if ((loadEveryNsample > 1) || (firstColIsSelectFlag))
	{
		sprintf(wfn,"%s.leftover",ofn);
		lf = fopen(wfn,"wb");
		if (!lf)
		{
			printf("Error:Cannot write leftover File %s\n",wfn);
			err = 1; 
		}
	}
	
	
	if (err == 0)
	{
		unsigned int rcnt = 0;
		unsigned short key = 0;
		
		if (binary)
		{
			unsigned int cluster;
			for (cluster = 1; cluster <= lastclusterid; cluster++)
			{
				printf("LOG: =====================================================================\n");
				printf("LOG: processing cluster %u\n",cluster);
				sprintf(wfn,"%s.%u.selected",ofn,cluster);
				wf = fopen(wfn,"wb+");
				if (!wf)
				{
					printf("Error:Cannot write Output File %s\n",wfn);
					err = 1;
				}
				else
				{
					cellnamecnt = processAssignedFile(f,wf,facsdata,facsname,cellnamecnt,cluster,totalrowcnt,&rcnt,colcnt,&key);
					err = WriteUniqueCellNames(wf,facsname,facsdata,cellnamecnt,rcnt,colcnt,key,0,kMAX_ALLOWED_INPUT_VALUE);
					fclose(wf);
				}
			}
		}
		else if (readFromUnassigned)
		{
				printf("LOG: =====================================================================\n");
				printf("LOG: processing unassigned file\n");
				sprintf(wfn,"%s.selected",ofn);
				wf = fopen(wfn,"wb+");
				if (!wf)
				{
					printf("Error:Cannot write Output File %s\n",wfn);
					err = 1;
				}
				else
				{
					cellnamecnt = processUnAssignedFile(f,wf,facsdata,facsname,cellnamecnt,totalrowcnt,&rcnt,colcnt,&key);
					err = WriteUniqueCellNames(wf,facsname,facsdata,cellnamecnt,rcnt,colcnt,key,0,kMAX_ALLOWED_INPUT_VALUE);
					fclose(wf);
				}
		}
		else /* regular .csv input file */
		{
			sprintf(wfn,"%s.selected",ofn);
			wf = fopen(wfn,"wb+");
			if (!wf)
			{
				printf("Error:Cannot write Output File %s\n",wfn);
				err = 1;
			}
			else
			{
				unsigned short minkeyval=0;
				unsigned short maxkeyval=(kMAX_ALLOWED_INPUT_VALUE-1);
				if (quickprocess)
					cellnamecnt = processInputFile(f,wf,lf,loadEveryNsample,facsdata,facsname,cellnamecnt,&rcnt,&colcnt,&key,firstColIsSelectFlag,&minkeyval,&maxkeyval);
				else
					cellnamecnt = SafeProcessInputFile(f,wf,lf,loadEveryNsample,facsdata,facsname,cellnamecnt,&rcnt,&colcnt,&key,firstColIsSelectFlag,&minkeyval,&maxkeyval);
				if (keyOverride != -1)
					key = keyOverride;
				if (cellnamecnt > 0)
				{
					err = WriteUniqueCellNames(wf,facsname,facsdata,cellnamecnt,rcnt,colcnt,key,minkeyval,maxkeyval);
				}
				else
				{
					err = 1;
				}
				fclose(wf);
			}
		}
		
	}


	free(facsdata);
	if (f)
		fclose(f);
	if (lf)
		fclose(lf);

	if (err)
		printf("Error\n");
	
	return(err);
	
	
} /* main */
/* ------------------------------------------------------------------------------------ */
