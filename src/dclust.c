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

	dclust.c: performs the parallel density based hirearchical clustering.

*/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <limits.h>
#include <time.h>
#include "mpi.h"
#include "dclust.h"
#include <pthread.h>

/* ------------------------------------------------------------------------------------ */

#define kMaxMergeRequests 786432	/* 6144 Kb */
#define kMaxCluster 1000000

#define kMaxCPU 129

#define kRawPrint   0x01
#define kSplitPrint 0x02
#define kStartLocalCluster 4000000  /* 4 millions in practice, allows for up to 512 processors but each cpu should not reach more than 4 millions distinct clusters */  

/* status of computations chunks */
#define kChunkStatusToDo 0
#define kChunkStatusComputing 1

/* MPI messages */
#define kWhichBlocksToCompute 0
#define kClusterMsg1 1
#define kClusterMsg2 2
#define kCPUdoneMsg 4
#define kFinalCntRequest 8
#define kFinalMergeRequestCntRequest 16
#define kSendFinalMergeRequestRequest 32
#define kJoinListRequest 64
#define kJoinListRequestDone 128
#define kRepeatWithNewDistMsg 512
#define kInitialClusterCntMsg 1024

#define kRenameClusterCntRequest 8192
#define kRenameClusterRequest 16384

#define kLeftoverData 32768
#define kLeftoverDataLength 32769
#define kLeftoverClusters 65536

/* status of cpu nodes */
#define kCPU_availaible 2147483647
#define kNoMoreBlocks   2147483647

#define kMaskedEvent 0

/* ------------------------------------------------------------------------------------ */
typedef	struct	FACSDATA_struct	FACSDATA;
struct	FACSDATA_struct
{
	unsigned short data[kMaxInputCol];
};
#ifdef COLUMNS_BY_32BLOCK
typedef	struct	FACSBLOCK32DATA_struct	FACSBLOCK32DATA;
struct	FACSBLOCK32DATA_struct
{
	unsigned short data[kMaxInputCol];
};
#endif

typedef	struct	FACSNAME_struct	FACSNAME;
struct	FACSNAME_struct
{
	CELLNAMEIDX condition;	
};


typedef	struct	CHUNK_struct	CHUNK;
struct	CHUNK_struct
{
	unsigned int ii;
	unsigned int jj;
	unsigned int iilast;
	unsigned int jjlast;
	unsigned int status;
};

typedef	struct	MERGECLUSTER_struct	MERGECLUSTER;
struct	MERGECLUSTER_struct
{
	 unsigned int cluster1;
	 unsigned int cluster2;
};

typedef	struct	MERGECLUSTERLIST_struct	MERGECLUSTERLIST;
struct	MERGECLUSTERLIST_struct
{
	 unsigned int cluster1;
	 unsigned int cluster2;
	 MERGECLUSTERLIST *next;
};

typedef	struct	MERGECLUSTERHASH_struct	MERGECLUSTERHASH;
struct	MERGECLUSTERHASH_struct
{
	unsigned int key;
	 unsigned int cluster1;
	 unsigned int cluster2;
};

typedef	struct	JOINREQUEST_struct	JOINREQUEST;
struct	JOINREQUEST_struct
{
	unsigned int getFromCPU;
	unsigned int sendToCPU;
};

typedef	struct	CPU_struct	CPU;
struct	CPU_struct
{
	unsigned int ii;
	unsigned int jj;
	unsigned int iilast;
	unsigned int jjlast;
};

typedef	struct	STATS_struct	STATS;
struct	STATS_struct
{
	float dist;
	float pctAssigned;
	int rawClustersCnt;
	int trimmedClustersCnt;
};

typedef	struct	CLUSTERHISTORY_struct	CLUSTERHISTORY;
struct	CLUSTERHISTORY_struct
{
	unsigned int pass;
	float dist;
	int cluster;
	int descendfrom;
	int mergedto;
	char retain;
};

typedef	struct EXECUTIONPLAN_struct  EXECUTIONPLAN;
struct EXECUTIONPLAN_struct
{
	FACSDATA *facsdata;
	unsigned int *clusterid;
	CPU chunk;
};


/* ------------------------------------------------------------------------------------ */

static unsigned int gTestDist;
static volatile unsigned int clustercnt;
static unsigned int mergerequestcnt;
static MERGECLUSTER mergerequest[kMaxMergeRequests];

static char header[kMaxLineBuf];
static char headerWithCluster[kMaxLineBuf];

/* those constants have to stay like that... */
#define kClusterTooSmall 2
#define kClusterEliminated 1
#define kClusterLargeEnough 0


static int *clustersnum[kMaxCPU];
static 	CLUSTERHISTORY *clusterhistory = NULL;
static unsigned int clusterhistorycnt = 0;
static int printwarnmergereq = 1;



#ifdef COLUMNS_BY_32BLOCK
static pthread_mutex_t mergingrequestmutex;
#endif

static pthread_mutex_t clustercntmutex;
static unsigned short sortkey;
static unsigned int thread_mergerequestcnt;

static MERGECLUSTER thread_mergerequest[kMaxMergeRequests];

/* ------------------------------------------------------------------------------------ */

static void DistributeLeftoverToClosestCluster(FACSDATA *facs, unsigned int *clusterid,unsigned int loaded, unsigned int colcnt,FACSDATA *leftoverfacs,unsigned int *leftoverclusterid,unsigned int leftoverloaded)
{
	unsigned int i,j,col;
	unsigned int dmin;
	unsigned int jmin;
	unsigned int ambiguousFlag;

	for (i = 0; i < leftoverloaded; i++)
	{
			FACSDATA *facspi = &leftoverfacs[i];
			dmin = 0xffffffff;
			jmin = i; // just to initialize something.
			ambiguousFlag = 0;

			// compute distance of unassigned (facspi) to each assigned (facspj) and record closest event of each cluster.
			for (j = 0; j < loaded; j++)
			{
				if (clusterid[j] > 0)
				{
					unsigned int d = 0;
					FACSDATA *facspj = &facs[j];
					for (col=0;col<colcnt;col++)
					{
						int diff = (int)facspj->data[col] - facspi->data[col];
						d += diff*diff;
					}
					// record closest distance
					if (d < dmin)
					{
						dmin = d;
						jmin = j;
						ambiguousFlag = 0;
					}
					else if ((d == dmin) && (clusterid[j] != clusterid[jmin]))
					{
						ambiguousFlag = 1;
					}
				}
			}
			// set the id of the closest cluster (but with an offset to make sure it will not affect loop above - cluster > maxclusterid for now)
			if (ambiguousFlag)
			{
				// unnecessary as we used calloc        leftoverclusterid[i] = 0;
			}
			else
			{
				if (dmin <= gTestDist) { leftoverclusterid[i] = clusterid[jmin]; }  // unnecessary as we used calloc    else { leftoverclusterid[i] = 0; } 
			}
	}

} // DistributeLeftoverToClosestCluster

/* ------------------------------------------------------------------------------------ */
static int DoProcessLeftoverbinaryFile(char *fn,FACSDATA *facs,unsigned int *clusterid,unsigned int loaded,unsigned int colcnt,int nproc,unsigned int verbose)
{
	unsigned int i,ccnt,leftoverrowcnt;
	int endian;
	char hdr[kHeaderSize];
	char lfn[kMaxFilename];
	char ofn[kMaxFilename];
	FILE *f = NULL;
	FILE *of = NULL;
	FACSDATA *leftoverfacs = NULL;
	FACSNAME *leftovername = NULL;
	unsigned int *leftoverclusterid = NULL;
	char *p;

		leftoverrowcnt = 0;
		sprintf(lfn,"%s",fn);
		p = strstr(lfn,"selected");
		if (p)
			strcpy(p,"leftover");
		f = fopen(lfn,"rb");
		if (!f)
		{
			printf("LOG: Warning: no leftover file to treat (%s)\n",lfn);
			goto bail;
		}
		
		fread(&hdr[0],sizeof(char),kHeaderSize,f);
		if (strcmp(hdr,"dclust unassigned file v1.0   \n") != 0)
		{
			printf("Error: input is not a dclust unassigned file\n");
			return(0);
		}	
		fread(&endian,sizeof(int),1,f);
		if (endian != 1)
		{
			printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
			return(0);
		}
		fread(&leftoverrowcnt,sizeof(int),1,f);
		fread(&ccnt,sizeof(int),1,f);
		if (ccnt != colcnt)
		{
			printf("Error: column count in binary file=%u, expected=%u\n",ccnt,colcnt);
			return(0);
		}
		/* read spacer */
		fread(&endian,sizeof(int),1,f);
		/* skip header */
		fseek(f,kMaxLineBuf,SEEK_CUR);

		leftoverfacs = calloc(leftoverrowcnt,sizeof(FACSDATA));
		if (!leftoverfacs)
		{
			printf("Error:Cannot Allocate Memory.\n");
			goto bail;
		}
		leftovername = calloc(leftoverrowcnt,sizeof(FACSNAME));
		if (!leftovername)
		{
			printf("Error:Cannot Allocate Memory.\n");
			goto bail;
		}
		leftoverclusterid = calloc(leftoverrowcnt,sizeof(FACSNAME));
		if (!leftoverclusterid)
		{
			printf("Error:Cannot Allocate Memory.\n");
			goto bail;
		}

		sprintf(ofn,"%s",fn);
		p = strstr(ofn,"selected");
		if (p)
			strcpy(p,"leftover.clusters");
		of = fopen(ofn,"w");
		if (!of)
		{
			printf("Error:Cannot Create Leftover Cluster results File %s\n",ofn);
			goto bail;
		}

		if (verbose > 0)
			printf("LOG: reading leftover file %s which contains %u events\n",lfn,leftoverrowcnt);
		fflush(stdout);
		for (i = 0; i < leftoverrowcnt; i++)
		{
			unsigned int j;
			float tmpfloat;
			fread(&leftovername[i],sizeof(CELLNAMEIDX),1,f);
			for (j = 0; j < ccnt; j++)
			{
				fread(&tmpfloat,sizeof(float),1,f);
				leftoverfacs[i].data[j] = (unsigned short)tmpfloat;
			}
		}
		fclose(f);

		if (verbose > 0)
			printf("LOG: processing leftover file %s which contains %u events\n",lfn,leftoverrowcnt);
		fflush(stdout);
		MPI_Bcast (&leftoverrowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (leftoverrowcnt > 0)
		{
			unsigned int ii;
			unsigned int starti,lasti;
			unsigned int datachunk = 1+(leftoverrowcnt/nproc);
			unsigned int actualChunkSize;
			MPI_Bcast (&clusterid[0], loaded, MPI_INT, 0, MPI_COMM_WORLD);

			for (ii = 1; ii<nproc; ii++)
			{
				starti = ii*datachunk;
				if (starti < leftoverrowcnt)
				{
					lasti = starti + datachunk;
					if (lasti > leftoverrowcnt)
						lasti = leftoverrowcnt;
					actualChunkSize = (lasti-starti);				
					MPI_Send(&actualChunkSize, 1, MPI_INT,  ii, kLeftoverDataLength, MPI_COMM_WORLD);
					MPI_Send(&leftoverfacs[starti].data[0], actualChunkSize*sizeof(FACSDATA),MPI_CHAR,  ii, kLeftoverData, MPI_COMM_WORLD);
				}
			}
			fflush(stdout);

			DistributeLeftoverToClosestCluster(facs,clusterid, loaded, colcnt, leftoverfacs, leftoverclusterid, datachunk);

			for (ii = 1; ii<nproc; ii++)
			{
				starti = ii*datachunk;
				if (starti < leftoverrowcnt)
				{
					lasti = starti + datachunk;
					if (lasti > leftoverrowcnt)
						lasti = leftoverrowcnt;
					fflush(stdout);
					MPI_Recv(&leftoverclusterid[starti], (lasti-starti), MPI_INT,  ii, kLeftoverClusters, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			}
			
		}
		if (verbose > 0)
			printf("LOG: writing cluster assignment of leftover events to file %s\n",ofn);
		fflush(stdout);
		// we do not want any header for easier subsequent merging... fprintf(of,"N,cluster\n"); 
		for (i = 0; i < leftoverrowcnt; i++)
		{
			fprintf(of,"%u,%u\n",leftovername[i].condition,leftoverclusterid[i]);
		}
		fclose(of);

		free(leftovername);
		free(leftoverfacs);		
		free(leftoverclusterid);	
		return(0);

bail :
		fflush(stdout);
		MPI_Bcast (&leftoverrowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (f)
			fclose(f);
		if (leftovername)
			free(leftovername);
		if (leftovername)
			free(leftoverfacs);
		if (leftoverclusterid)
			free(leftoverclusterid);
		return(1);

}

/* ------------------------------------------------------------------------------------ */

static void thread_InsertMergeRequest(unsigned int cluster1, unsigned int cluster2)  /* note: cluster2 > cluster1 is always true */
{
	unsigned int ii;
	
		if (thread_mergerequestcnt == 0)
		{
			thread_mergerequest[0].cluster1 = cluster1;
			thread_mergerequest[0].cluster2 = cluster2;
			thread_mergerequestcnt++;
			return;
		}
	
		if (thread_mergerequestcnt < kMaxMergeRequests)
		{
			/* quickly find position of cluster2 - 1 using a binary divide and conquer... */
			int left;
			unsigned int val;
			int right = thread_mergerequestcnt-1;

			left = 0;
			val=cluster2-1;
			do
			{
				ii=(left+right)>>1;
				if (thread_mergerequest[ii].cluster2 > val)
					right=ii-1;
				else
					left=ii+1;
			} while ((thread_mergerequest[ii].cluster2 != val) && (left <= right));
			
			/* move up from this point to do our real check, using cluster2... */
			/* test if and where to insert request */
			for( ;  ii< thread_mergerequestcnt; ii++)
			{
				if 	 (thread_mergerequest[ii].cluster2 > cluster2)
					break;  /* we can insert and be done */
				if (thread_mergerequest[ii].cluster2 == cluster2) 
				{
					if (thread_mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (thread_mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = thread_mergerequest[ii].cluster1;
						thread_mergerequest[ii].cluster1 = cluster1;
						thread_InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						thread_InsertMergeRequest(thread_mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&thread_mergerequest[ii+1].cluster1, &thread_mergerequest[ii].cluster1, (thread_mergerequestcnt-ii)*sizeof(MERGECLUSTER));

			thread_mergerequest[ii].cluster1 = cluster1;
			thread_mergerequest[ii].cluster2 = cluster2;
			thread_mergerequestcnt++;
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				printf("LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}


} /* thread_InsertMergeRequest */
/* ------------------------------------------------------------------------------------ */

static void *computesim_funcion(void *ep)
{
	FACSDATA *facspi;
	FACSDATA *facspj;
	unsigned int *clusterpi;
	unsigned int *clusterpj;
	unsigned int i,j;
	unsigned int starti,startj;
	unsigned int lasti,lastj;
	unsigned  int actualstartj;

	FACSDATA *facs = ((EXECUTIONPLAN*)ep)->facsdata;
	unsigned int *clusterid = ((EXECUTIONPLAN*)ep)->clusterid;
	CPU *chunk = &((EXECUTIONPLAN*)ep)->chunk;

	starti = chunk->ii;
	lasti  = chunk->iilast;
	startj = chunk->jj;
	lastj  = chunk->jjlast;


		{
				unsigned short valii,valjj;

				facspi = &facs[lasti-1];
				facspj = &facs[startj];
				valii = facspi->data[sortkey];
				valjj = facspj->data[sortkey];
				if (valjj > valii)
				{
					unsigned  int mindist = valjj-valii;
					mindist *= mindist;
					if (mindist > gTestDist)
					{
						goto skipchunk;
					}
				}

		}


		facspi = &facs[starti];
		clusterpi = &clusterid[starti];
		for(i = starti;  i< lasti; i++)
		{
			if (startj <= i)
				actualstartj = i+1;
			else 
				actualstartj = startj;
			facspj = &facs[actualstartj];
			clusterpj = &clusterid[actualstartj];
			for(j = actualstartj;  j< lastj; j++)
			{
				int diff;
				unsigned int d;

				if ((*clusterpj) && (*clusterpj == *clusterpi))
					goto straight2next;


				diff = (int)facspj->data[0] - facspi->data[0];
				d = diff*diff;
				diff = (int)facspj->data[1] - facspi->data[1];
				d += diff*diff;
				diff = (int)facspj->data[2] - facspi->data[2];
				d += diff*diff;
				diff = (int)facspj->data[3] - facspi->data[3];
				d += diff*diff;
	#if defined(COLUMNS_8) || defined(COLUMNS_12) || defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[4] - facspi->data[4];
				d += diff*diff;
				diff = (int)facspj->data[5] - facspi->data[5];
				d += diff*diff;
				diff = (int)facspj->data[6] - facspi->data[6];
				d += diff*diff;
				diff = (int)facspj->data[7] - facspi->data[7];
				d += diff*diff;
	#endif
	#if defined(COLUMNS_12) || defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[8] - facspi->data[8];
				d += diff*diff;
				diff = (int)facspj->data[9] - facspi->data[9];
				d += diff*diff;
				diff = (int)facspj->data[10] - facspi->data[10];
				d += diff*diff;
				diff = (int)facspj->data[11] - facspi->data[11];
				d += diff*diff;
	#endif
	#if defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[12] - facspi->data[12];
				d += diff*diff;
				diff = (int)facspj->data[13] - facspi->data[13];
				d += diff*diff;
				diff = (int)facspj->data[14] - facspi->data[14];
				d += diff*diff;
				diff = (int)facspj->data[15] - facspi->data[15];
				d += diff*diff;
	#endif 
	#if  defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
	#if  defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				/* computation done for the first 16 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;
	#endif
				#undef kDATAOFFSET
				#define kDATAOFFSET 16
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+4] - facspi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+5] - facspi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+6] - facspi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+7] - facspi->data[kDATAOFFSET+7];
				d += diff*diff;
	#if defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[kDATAOFFSET+8] - facspi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+9] - facspi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+10] - facspi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+11] - facspi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+12] - facspi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+13] - facspi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+14] - facspi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+15] - facspi->data[kDATAOFFSET+15];
				d += diff*diff;
	#endif
	#endif
	#if    defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				/* computation done for the first 32 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;

				#undef kDATAOFFSET
				#define kDATAOFFSET 32
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+4] - facspi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+5] - facspi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+6] - facspi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+7] - facspi->data[kDATAOFFSET+7];
				d += diff*diff;
	  #if    defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[kDATAOFFSET+8] - facspi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+9] - facspi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+10] - facspi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+11] - facspi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+12] - facspi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+13] - facspi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+14] - facspi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+15] - facspi->data[kDATAOFFSET+15];
				d += diff*diff;
	  #endif
	#endif
	#if   defined(COLUMNS_52)
			#undef kDATAOFFSET
			#define kDATAOFFSET 48
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
	#endif

				if (d <= gTestDist) 
				{
					if (*clusterpi)
					{
						if (*clusterpj == 0)
						{
							*clusterpj = *clusterpi;
						}
						else /* i=assigned, j=assigned */
						{
							unsigned int cluster1,cluster2;
							if(*clusterpj > *clusterpi)
							{
								cluster1 = *clusterpi;
								cluster2 = *clusterpj;
							}
							else
							{
								cluster1 = *clusterpj;
								cluster2 = *clusterpi;
							}
							thread_InsertMergeRequest(cluster1,cluster2);
						}
						
					}
					else
					{
						if (*clusterpj == 0)   /* i=not yes assigned, j=not yet assigned */
						{
							pthread_mutex_lock(&clustercntmutex);
							*clusterpi = ++clustercnt;
							*clusterpj = clustercnt;
							pthread_mutex_unlock(&clustercntmutex);
						}
						else /* i=not yes assigned, j=assigned */
						{
							*clusterpi = *clusterpj;
						}
					}
				}
		straight2next:		
				facspj++;
				clusterpj++;
			}
			facspi++;
			clusterpi++;
		}
skipchunk:
	pthread_exit(NULL);

} /* computesim_funcion */
/* ------------------------------------------------------------------------------------ */

static void DoMergeLists(int nproc,int verbose)
{
		JOINREQUEST joinRequest;
		int cpu,sendingCPUoffset;
		unsigned int received,joiningprocesses;
		
		sendingCPUoffset = 1;
		do
		{
			cpu = 1;

			if ((cpu + sendingCPUoffset) >= nproc)
				break;

			if (verbose > 2)
			{	
				printf("LOG:Aggregating merge requests (cpus +%d)\n",sendingCPUoffset);				
				fflush(stdout);
			}
			/* send instruction to join list */
			joiningprocesses = 0;
			do
			{
				if ((cpu + sendingCPUoffset) >= nproc)
					break;
				joinRequest.getFromCPU = cpu + sendingCPUoffset;
				joinRequest.sendToCPU = cpu;				
				MPI_Send(&joinRequest, 2, MPI_INT,  joinRequest.sendToCPU,  kJoinListRequest, MPI_COMM_WORLD);
				MPI_Send(&joinRequest, 2, MPI_INT,  joinRequest.getFromCPU, kJoinListRequest, MPI_COMM_WORLD);
				cpu += sendingCPUoffset+sendingCPUoffset;
				joiningprocesses++;
			} while (1);

			/* wait until done */
			received = 0;
			cpu = 1;
			do
			{
				unsigned int done;
				if ((cpu + sendingCPUoffset) >= nproc)
					break;
				MPI_Recv(&done, 1, MPI_INT,  cpu, kJoinListRequestDone, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
				cpu += sendingCPUoffset+sendingCPUoffset;
			} while (++received < joiningprocesses);

			sendingCPUoffset *= 2;

		} while(1);

		joinRequest.getFromCPU  = 0;
		joinRequest.sendToCPU = 0;
		for (cpu = 1; cpu < nproc; cpu++)
			MPI_Send(&joinRequest, 2, MPI_INT,  cpu,  kJoinListRequest, MPI_COMM_WORLD);

} /* DoMergeLists */

/* ------------------------------------------------------------------------------------ */

static unsigned int dclustFileReadHeader(FILE *f,int idproc,unsigned int *rcnt,unsigned int *ccnt,unsigned short *colkey,unsigned int *loadEveryNsample)
{
	unsigned int rowcnt;
	unsigned int colcnt;
	unsigned short cellnamecnt,key;
	int endian;
	unsigned int skip;
	char hdr[kHeaderSize];
	int err = 1 /* always pessimist*/ ;

	fread(&hdr[0],sizeof(char),kHeaderSize,f);
	if (strcmp(hdr,"dclust input file v1.0        \n") != 0)
	{
		if (idproc == 0)
			printf("Error: input is not a dclust input file\n");
		return(err);
	}	

	fread(&endian,sizeof(int),1,f);
	if (endian != 1)
	{
		if (idproc == 0)
			printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(err);
	}
	fread(&rowcnt,sizeof(int),1,f);
	fread(&colcnt,sizeof(float),1,f);
	if (colcnt > kMaxInputCol)
	{
		if (idproc == 0)
			printf("LOG:Fatal: number of requested columns exceed maximum allowed (%u > %d)\n",colcnt,kMaxInputCol);
		return(err);
	}	
	fread(&skip,sizeof(int),1,f);

	/* read header */
	fread(&header,sizeof(char),kMaxLineBuf,f);
	sprintf(headerWithCluster,"%s,cluster",header);

	/* read column key used to sort data */
	fread(&key,sizeof(unsigned short),1,f);
	/* read number of unique cellnames */
	fread(&cellnamecnt,sizeof(unsigned short),1,f);

	/* skip cellnames to point on actual data */
	fseek(f,((unsigned int)cellnamecnt*kMaxCellName*sizeof(char)),SEEK_CUR);

	err = 0;
	*rcnt = rowcnt;
	*ccnt = colcnt;
	*colkey = key;
	*loadEveryNsample = skip;
	return(err);
	
} /* dclustFileReadHeader */
/* ------------------------------------------------------------------------------------ */

static unsigned int loaddata(FILE *f, FACSNAME *facsname, FACSDATA *facs,unsigned int rowcnt,unsigned int colcnt,int idproc)
{
	FACSDATA *facsp;
	FACSNAME *facsnamep;
	unsigned int cnt;

	cnt = 0;

	facsp = &facs[0];	
	if (idproc == 0)
	{
		facsnamep = &facsname[0];
		for (cnt = 0; cnt < rowcnt; cnt++)
		{
			fread(facsnamep,sizeof(CELLNAMEIDX),1,f);
			fread(facsp,sizeof(short),colcnt,f);
			facsp++;
			facsnamep++;
		}
	}
	else
	{
		for (cnt = 0; cnt < rowcnt; cnt++)
		{
			/* skip nameindex */
			fseek(f,sizeof(CELLNAMEIDX),SEEK_CUR);
			fread(facsp,sizeof(short),colcnt,f);
			facsp++;
		}
	}
	
	return(cnt);
	
} /* loaddata */
/* ------------------------------------------------------------------------------------ */

static unsigned int CountAssigned(unsigned int *clusterid,unsigned int rowcnt,unsigned int maxclusterid)
{
	unsigned int *clusteridp;
	unsigned int i;
	unsigned int assigned = 0;

			clusteridp = &clusterid[0];
			for(i = 0;  i< rowcnt; i++)
			{
					if ((*clusteridp > 0) && (*clusteridp<=maxclusterid))
						assigned++;
					clusteridp++;
			}//i

		return(assigned);
} /* CountAssigned */
/* ------------------------------------------------------------------------------------ */
static void WriteClusterIndices(unsigned int *clusterid,unsigned int rowcnt,float distcutoff,char *dir)
{
	unsigned int *clusteridp;
	FILE *af;
	char fn[kMaxFilename];

		sprintf(fn,"%s-%.6f",dir,distcutoff);
		af=fopen(fn,"wb");
		if (af)
		{
			clusteridp = &clusterid[0];
			fwrite(clusteridp,sizeof(int),rowcnt,af);
			fclose(af);
		}
		else
			printf("LOG: Error writing %s\n",fn);
	
} /* WriteClusterIndices */
/* ------------------------------------------------------------------------------------ */
static unsigned int WriteSplitBinFile(FACSNAME *facsname,FACSDATA *facs,unsigned int *clusterid,unsigned int rowcnt,unsigned int colcnt,unsigned int maxclusterid,char *ofn)
{
	unsigned short *usp;
	unsigned int *clusteridp;
	unsigned int i;
	unsigned int col;
	FILE *uf;
	FILE *af;
	FACSDATA *facsp;
	FACSNAME *facsnamep;
	char fn[kMaxFilename];
	char hdr[kHeaderSize];
	float val[kMaxInputCol];
	unsigned int assigned = 0;
	int endian = 1;

		sprintf(fn,"%s.assigned",ofn);
		af=fopen(fn,"wb");
		sprintf(fn,"%s.unassigned",ofn);
		uf=fopen(fn,"wb");

		if (af && uf)
		{

			strcpy(hdr,"dclust assigned file v1.0     \n");
			fwrite(&hdr[0],sizeof(char),kHeaderSize,af);
			strcpy(hdr,"dclust unassigned file v1.0   \n");
			fwrite(&hdr[0],sizeof(char),kHeaderSize,uf);

			fwrite(&endian,sizeof(int),1,af);
			fwrite(&endian,sizeof(int),1,uf);
			endian--;
			fwrite(&endian,sizeof(int),1,af); // keep space for number of rows
			fwrite(&endian,sizeof(int),1,uf); // keep space for number of rows

			fwrite(&colcnt,sizeof(int),1,af);
			fwrite(&colcnt,sizeof(int),1,uf);

			fwrite(&maxclusterid,sizeof(int),1,af); 
			fwrite(&maxclusterid,sizeof(int),1,uf);

			fwrite(&header,sizeof(char),kMaxLineBuf,uf);
			fwrite(&headerWithCluster,sizeof(char),kMaxLineBuf,af);

			clusteridp = &clusterid[0];
			facsp = &facs[0];
			facsnamep = &facsname[0];
			for(i = 0;  i< rowcnt; i++)
			{
					usp = &facsp->data[0];
					if ((*clusteridp > 0) && (*clusteridp<=maxclusterid))
					{
						fwrite(facsnamep,sizeof(CELLNAMEIDX),1,af);
						for (col=0;col<colcnt;col++)
						{
								val[col] = (float)(*usp);
								usp++;
						}
						fwrite(val,sizeof(float),colcnt,af);
						fwrite(clusteridp,sizeof(int),1,af);
						assigned++;
					}
					else
					{
						fwrite(facsnamep,sizeof(CELLNAMEIDX),1,uf);
						for (col=0;col<colcnt;col++)
						{
								val[col] = (float)(*usp);
								usp++;
						}
						fwrite(val,sizeof(float),colcnt,uf);
					}
					clusteridp++;
					facsp++;
					facsnamep++;
			}//i
		}
		if (af)
		{
			fseek(af,36L,SEEK_SET);
			fwrite(&assigned,sizeof(int),1,af);
			fclose(af);
		}
		else
			printf("LOG: Cannot write results to %s.assigned\n",ofn);
		if (uf)
		{
			fseek(uf,36L,SEEK_SET);
			rowcnt -= assigned;
			fwrite(&rowcnt,sizeof(int),1,uf);
			fclose(uf);
		}
		else
			printf("LOG: Cannot write results to %s.unassigned\n",ofn);
		return(assigned);

} /* WriteSplitBinFile */
/* ------------------------------------------------------------------------------------ */

static void removeclustersnum(int id)
{
	int *cnp;
	int cpu =   id / kStartLocalCluster ;
	int cluster = id - cpu*kStartLocalCluster;
	cnp = clustersnum[cpu];
	cnp[cluster] = kClusterEliminated;  /* flag this specific cluster id as being removed */

} /* removeclustersnum */
/* ------------------------------------------------------------------------------------ */
static void AdjustClustersID(unsigned int *clusters, unsigned int loaded,unsigned int nproc,int verbose,int trimmedclustercnt,int *firstAvailClusterID)
{
	unsigned int i;
	int j;
	int *cnp;
	unsigned int clusterid = 1; /* first cluster will have id=1 */
	unsigned int tinyClustersId = trimmedclustercnt +1 ;  /* start to pile up number of clusters too small to pass the min size cutoff after "good" clusters */
	unsigned int ii;

	/* loop over each cpu, which has attributed its own ids */
	for (i = 1; i < nproc; i++)
	{
		cnp = clustersnum[i];
		for (j = 1; j<=cnp[0]; j++)
		{
			if (cnp[j]==kClusterLargeEnough)
			{
				/* rename clusterid (attribute clusterid as a new id) */
				if (verbose > 2)
					printf("LOG:Renaming clusterid %u to %u\n",i*kStartLocalCluster+j,clusterid);
				cnp[j]=clusterid++;
			}
			else if (cnp[j]==kClusterTooSmall)
			{
				cnp[j]=tinyClustersId++;
			}
			else
				cnp[j] = 0;
		}

		for(ii = 0;  ii< loaded; ii++)
		{
			if ((int)clusters[ii] >= i*kStartLocalCluster)
			{
				int idx = (int)clusters[ii] - i*kStartLocalCluster;
				if (idx < kStartLocalCluster)
				{
					clusters[ii] = cnp[idx];
				}
			}
		}
		
	}
	

	*firstAvailClusterID = tinyClustersId-1;
	
} /* AdjustClustersID */
/* ------------------------------------------------------------------------------------ */
static unsigned int RemoveSmallClusters(unsigned int *clusters,unsigned int loaded,unsigned int nproc,unsigned int minevents)
{
	unsigned int ii;
	int i,j;
	int *cnp;
	unsigned int retainedClusterCnt = 0;

	/* loop over each cpu, which has attributed its own ids */
	for (i = 1; i < nproc; i++)
	{
		unsigned int *cnt;
		cnp = clustersnum[i];
		cnt = calloc((cnp[0]+1),sizeof(unsigned int));

		for(ii = 0;  ii< loaded; ii++)
		{
			int idx = (int)clusters[ii] - i*kStartLocalCluster;
			if ((idx >= 0) && (idx <= cnp[0])) /* cluster belongs to proc i */
				cnt[idx]++;
		}
		for (j = 1; j<=cnp[0]; j++)
		{
			if (cnt[j] < minevents)
			{
				if (cnp[j] != kClusterEliminated)
					cnp[j]=kClusterTooSmall;  /* flag this specific cluster id as being removed */
			}
			else
			{
				retainedClusterCnt++;
			}
		}
		free(cnt);
	}
	return(retainedClusterCnt);
	
} /* RemoveSmallClusters */
/* ------------------------------------------------------------------------------------ */

static void PrintClusterStatus(CLUSTERHISTORY *clusterhistory,unsigned int clusterhistorycnt)
{
	unsigned int ii;
	for (ii = 0; ii<clusterhistorycnt; ii++)
		printf("LOG: %4u %6.2f | %3d ^ %3d -> %3d | %c\n",clusterhistory[ii].pass,clusterhistory[ii].dist,clusterhistory[ii].cluster,clusterhistory[ii].descendfrom,clusterhistory[ii].mergedto,clusterhistory[ii].retain);
		
				
} /* PrintClusterStatus */
/* ------------------------------------------------------------------------------------ */

static void DistributeUnassignedToClosestCluster(FACSDATA *facs, unsigned int *clusterid,unsigned int loaded, unsigned int colcnt,unsigned int maxclusterid,unsigned int starti,unsigned int lasti)
{
	unsigned int i,j,col;
	unsigned int dmin;
	unsigned int jmin;
	unsigned int ambiguousFlag;
	unsigned int ambiguousCnt = 0;
	unsigned int reassigned = 0;

	for (i = starti; i < lasti; i++)
	{
		if (clusterid[i] == 9999999) /* flagged for reassignment */
		{
			FACSDATA *facspi = &facs[i];
			dmin = 0xffffffff;
			jmin = i; // just to initialize something.
			ambiguousFlag = 0;

			// compute distance of unassigned (facspi) to each assigned (facspj) and record closest event of each cluster.
			for (j = 0; j < loaded; j++)
			{
				if ((clusterid[j] > 0) && (clusterid[j] <= maxclusterid)) /* assigned */
				{
					unsigned int d = 0;
					FACSDATA *facspj = &facs[j];
					for (col=0;col<colcnt;col++)
					{
						int diff = (int)facspj->data[col] - facspi->data[col];
						d += diff*diff;
					}
					// record closest distance
					if (d < dmin)
					{
						dmin = d;
						jmin = j;
						ambiguousFlag = 0;
					}
					else if ((d == dmin) && (clusterid[j] != clusterid[jmin]))
					{
						ambiguousFlag = 1;
					}
				}
			}
			// set the id of the closest cluster (but with an offset to make sure it will not affect loop above - cluster > maxclusterid for now)
			if (ambiguousFlag)
			{
				clusterid[i] = 0 + maxclusterid + 1;
				ambiguousCnt++;
			}
			else
			{
				clusterid[i] = clusterid[jmin] + maxclusterid + 1;
				reassigned++;
			}
		}
	}
	// adjust clusterid of reassigned sequences.
	for (i = starti; i < lasti; i++)
	{
		if (clusterid[i] > maxclusterid)
			clusterid[i] -= (maxclusterid + 1);
	}


} // DistributeUnassignedToClosestCluster

/* ------------------------------------------------------------------------------------ */
static unsigned int FlagSequencesToReassign(unsigned int *clusterid,unsigned int loaded,char *fn)
{
	unsigned int i;
	unsigned int cnt = 0;

	FILE *f=fopen(fn,"rb");
	if (!f)
	{
		printf("LOG: ERROR Cannot open last clustering distance status file (%s)\n",fn);
		return(0);
	}
	
	for(i = 0;  i< loaded; i++)
	{
		unsigned int cl;
		fread(&cl,sizeof(int),1,f);
		if (cl == 0) // make sure that sequence was assigned at last step.
			continue;
		if (clusterid[i]==0) // need to reassign.
		{
			clusterid[i] = 9999999;
			cnt++;
		}
	}
	fclose(f);
	return(cnt);
	
} /* FlagSequencesToReassign */
/* ------------------------------------------------------------------------------------ */
static void UpdateClusterHistory(unsigned int pass, int clustercnt,int prevclustercnt, float dist,int verbose)
{
	unsigned int ii;
	unsigned int prevcnt = clusterhistorycnt;
	
		if (verbose > 1)
			printf("LOG: clusterhistorycnt = %u  (pass %u)\n",clusterhistorycnt,pass);
		if (pass == 0)
		{
			for (ii = 1; ii<=clustercnt; ii++)
			{
				clusterhistory[clusterhistorycnt].pass = pass;
				clusterhistory[clusterhistorycnt].dist = dist;
				clusterhistory[clusterhistorycnt].cluster = ii;
				clusterhistory[clusterhistorycnt].descendfrom = 0;
				clusterhistory[clusterhistorycnt].mergedto = 0;
				clusterhistory[clusterhistorycnt].retain = ' ';
				clusterhistorycnt++;
			}
		}
		else
		{			
			int parent = 1;
			for (ii = 1; ii<=clustercnt; ii++)
			{
				unsigned int x;
				
				if (parent <= prevclustercnt)
				{
					for (x = 0; x<prevcnt; x++)
					{
						if ((clusterhistory[x].pass == (pass-1)) && (clusterhistory[x].cluster == parent))
						{
							if (clusterhistory[x].mergedto != 0)
								parent++;
						}
					}
					if (parent <= prevclustercnt)
						clusterhistory[clusterhistorycnt].descendfrom = parent;
					else
						clusterhistory[clusterhistorycnt].descendfrom = 0;
					parent++;
				}
				else
					clusterhistory[clusterhistorycnt].descendfrom = 0;
				clusterhistory[clusterhistorycnt].pass = pass;
				clusterhistory[clusterhistorycnt].dist = dist;
				clusterhistory[clusterhistorycnt].cluster = ii;
				clusterhistory[clusterhistorycnt].mergedto = 0;
				clusterhistory[clusterhistorycnt].retain = ' ';
				clusterhistorycnt++;
			}
		}
		
		
} /* UpdateClusterHistory */
/* ------------------------------------------------------------------------------------ */
static char CheckClusterNotYetRetained(int from,int cluster)
{
	int x;
	for (x = from-1; x >= 0; x--)
	{
		if (clusterhistory[x].cluster == cluster)
		{
			if (clusterhistory[x].retain != ' ')
				return('n');
			if (clusterhistory[x].descendfrom == 0)
				return('y');
			return(CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom));
		}
	}
	return('!');
	
} /* CheckClusterNotYetRetained */
/* ------------------------------------------------------------------------------------ */

static int SelectClusterHistory(unsigned int rowcnt, unsigned int *clusterid,char *dir,int verbose)
{
	unsigned int ii;
	int x;
	int cnt = 1;
	
	for (ii = 0; ii<clusterhistorycnt; ii++)
	if (clusterhistory[ii].mergedto)
	{
		if (clusterhistory[ii].descendfrom == 0)
			clusterhistory[ii].retain = 'y';
		else
			clusterhistory[ii].retain = CheckClusterNotYetRetained(ii,clusterhistory[ii].descendfrom);
			
		for (x = 0; x<ii; x++)
		{
			if ((clusterhistory[x].pass == (clusterhistory[ii].pass )) && (clusterhistory[x].cluster == clusterhistory[ii].mergedto))
			{
				if (clusterhistory[x].descendfrom == 0)
					clusterhistory[x].retain = 'y';
				else
					clusterhistory[x].retain = CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom);
				break;
			}
		}
		
	}


	for (x = clusterhistorycnt-1; x >= 0; x--)
	{
		if (clusterhistory[x].pass != clusterhistory[clusterhistorycnt-1].pass)
			break;
		if (clusterhistory[x].retain != ' ')
			continue;
		clusterhistory[x].retain = CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom);
	}
	
	
	printf("LOG:**************************************************************\n");
	for (ii = 0; ii<clusterhistorycnt; ii++)
	{
		if (clusterhistory[ii].retain == 'y')
		{
			char fn[kMaxFilename];
			FILE *f = NULL;
			unsigned int evtcnt = 0;
			
			sprintf(fn,"%s-%.6f",dir,clusterhistory[ii].dist);
			f=fopen(fn,"rb");
			if (verbose > 1)
				printf("%s.%d\n",fn,clusterhistory[ii].cluster);
			if (f)
			{
				unsigned int i;
				for(i = 0;  i< rowcnt; i++)
				{
					unsigned int cl;
					fread(&cl,sizeof(int),1,f);
					if (cl == clusterhistory[ii].cluster)
					{
						clusterid[i] = cnt;
						evtcnt++;
					}
				} 
				fclose(f);
			}
			else
				printf("LOG: Error reading %s\n",fn);

			printf("LOG: Cluster %5d: %10u events\n",cnt,evtcnt);
			cnt++;
		}
	}

	return(cnt-1);
		
} /* SelectClusterHistory */
/* ------------------------------------------------------------------------------------ */


static unsigned int ProcessMergeRequests(unsigned int *clusterid,unsigned int loaded,unsigned int mergerequestcnt,unsigned int nproc,int previouslyRetainedClusterCnt,float previousSamplingDist,unsigned int pass)
{
	unsigned int ii;
	unsigned  int isMergingPreexistingClusters = 0;
	if (previouslyRetainedClusterCnt > -1)
	{
		/* identifies which preexisting clusters might have been merged and should subsequently be split */
		for (ii = 0; ii<mergerequestcnt; ii++)
		{
			if (mergerequest[ii].cluster2 <= (1*kStartLocalCluster+previouslyRetainedClusterCnt))
			{
				unsigned int b;

				isMergingPreexistingClusters = 1;
				for (b = 0; b< clusterhistorycnt; b++)
				{
					if ((clusterhistory[b].pass == (pass-1)) && clusterhistory[b].cluster == (mergerequest[ii].cluster2-(1*kStartLocalCluster)))
						clusterhistory[b].mergedto = (mergerequest[ii].cluster1-(1*kStartLocalCluster));
				}
			}
		}
	}
	
{
	int i,j,mrg;
	// init has table: no change.
	for (i = 1; i < nproc; i++)
	{
		int *cnp = clustersnum[i];
		for (j = 1; j<=cnp[0]; j++)
		{
			cnp[j]=i*kStartLocalCluster+j;
		}
	}
	// update table according to merge request processing.
	for (mrg = mergerequestcnt-1; mrg >= 0 ; mrg--)
	{
		for (i = 1; i < nproc; i++)
		{
			int *cnp = clustersnum[i];
			for (j = 1; j<=cnp[0]; j++)
			{
				if (cnp[j]==mergerequest[mrg].cluster2)
					cnp[j] = mergerequest[mrg].cluster1;
			}
		}
	}
	// update clusterids
	for(ii = 0;  ii< loaded; ii++)
	{
		if (clusterid[ii] > 0)
		{
			int cpuid = (int)clusterid[ii] / kStartLocalCluster;
			int hashidx = (int)clusterid[ii] - cpuid * kStartLocalCluster;
			int *cnp = clustersnum[cpuid];
			
			clusterid[ii] = cnp[hashidx];
		}
	}

	// get back to calloc state
	for (i = 1; i < nproc; i++)
	{
		int *cnp = clustersnum[i];
		bzero(&cnp[1], cnp[0]*sizeof(int));
	}
}
	
	return(isMergingPreexistingClusters);
	
} /* ProcessMergeRequests */
/* ------------------------------------------------------------------------------------ */

static void InsertMergeRequest(unsigned int cluster1, unsigned int cluster2)
{
		if (mergerequestcnt == 0)
		{
			mergerequest[0].cluster1 = cluster1;
			mergerequest[0].cluster2 = cluster2;
			mergerequestcnt++;
			return;
		}
		
		
		if (mergerequestcnt < kMaxMergeRequests)
		{
			/* quickly find position of cluster2 - 1 using a binary divide and conquer... */
			unsigned int ii;
			         int left = 0;
			         int right = mergerequestcnt-1;
			unsigned int val=cluster2-1;
			do
			{
				ii=(left+right)>>1;
				if (mergerequest[ii].cluster2 > val)
					right=ii-1;
				else
					left=ii+1;
			} while ((mergerequest[ii].cluster2 != val) && (left <= right));
			
			/* move up from this point to do our real check, using cluster2... */
			/* test if and where to insert request */
			for( ;  ii< mergerequestcnt; ii++)
			{
				if 	 (mergerequest[ii].cluster2 > cluster2)
					break;  /* we can insert and be done */
				if (mergerequest[ii].cluster2 == cluster2) 
				{
					if (mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = mergerequest[ii].cluster1;
						mergerequest[ii].cluster1 = cluster1;
						InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						InsertMergeRequest(mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&mergerequest[ii+1].cluster1, &mergerequest[ii].cluster1, (mergerequestcnt-ii)*sizeof(MERGECLUSTER));

			mergerequest[ii].cluster1 = cluster1;
			mergerequest[ii].cluster2 = cluster2;
			mergerequestcnt++;
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				printf("LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}


} /* InsertMergeRequest */

/* ------------------------------------------------------------------------------------ */
static void InsertMergeRequestWhere(unsigned int cluster1, unsigned int cluster2,unsigned int *where)
{

		if (mergerequestcnt < kMaxMergeRequests)
		{
			unsigned int ii;
			
			/* test if and where to insert request (start from where we inserted previous one as lists are sorted according to cluster2  */
			for(ii = *where;  ii< mergerequestcnt; ii++)
			{
				if 	 (mergerequest[ii].cluster2 > cluster2) 
					break;  /* we can insert and be done */
				if (mergerequest[ii].cluster2 == cluster2) 
				{
					if (mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = mergerequest[ii].cluster1;
						mergerequest[ii].cluster1 = cluster1;
						InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						InsertMergeRequest(mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&mergerequest[ii+1].cluster1, &mergerequest[ii].cluster1, (mergerequestcnt-ii)*sizeof(MERGECLUSTER));
			mergerequest[ii].cluster1 = cluster1;
			mergerequest[ii].cluster2 = cluster2;
			mergerequestcnt++;
			*where = ii; /* remember where we did insert to save time */
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				printf("LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}


} /* InsertMergeRequestWhere */
/* ------------------------------------------------------------------------------------ */

#ifdef COLUMNS_BY_32BLOCK
static void computesim(FACSDATA *facs, unsigned int *clusterid, CPU *chunk)
{
	FACSDATA *facspi;
	FACSDATA *facspj;

	unsigned int *clusterpi;
	unsigned int *clusterpj;
	unsigned int i,j;
	unsigned int starti,startj;
	unsigned int lasti,lastj;
	unsigned  int actualstartj;
	
	starti = chunk->ii;
	startj = chunk->jj;
	lasti = chunk->iilast;
	lastj = chunk->jjlast;
	facspi = &facs[starti];
	clusterpi = &clusterid[starti];
	for(i = starti;  i< lasti; i++)
	{
		if (startj <= i)
			actualstartj = i+1;
		else 
			actualstartj = startj;
		facspj = &facs[actualstartj];
		clusterpj = &clusterid[actualstartj];
		for(j = actualstartj;  j< lastj; j++)
		{
			FACSBLOCK32DATA  *facspi32 = (FACSBLOCK32DATA*)facspi;
			FACSBLOCK32DATA  *facspj32 = (FACSBLOCK32DATA*)facspj;
			unsigned int b = kBLOCK32CNT;
			unsigned int dsum = 0;
			while(b-- > 0)  /* do 32 dimensions blocks */
			{ 
				int diff;
				unsigned int d;

				diff = (int)facspj32->data[0] - facspi32->data[0];
				d = diff*diff;
				diff = (int)facspj32->data[1] - facspi32->data[1];
				d += diff*diff;
				diff = (int)facspj32->data[2] - facspi32->data[2];
				d += diff*diff;
				diff = (int)facspj32->data[3] - facspi32->data[3];
				d += diff*diff;
				diff = (int)facspj32->data[4] - facspi32->data[4];
				d += diff*diff;
				diff = (int)facspj32->data[5] - facspi32->data[5];
				d += diff*diff;
				diff = (int)facspj32->data[6] - facspi32->data[6];
				d += diff*diff;
				diff = (int)facspj32->data[7] - facspi32->data[7];
				d += diff*diff;

				diff = (int)facspj32->data[8] - facspi32->data[8];
				d += diff*diff;
				diff = (int)facspj32->data[9] - facspi32->data[9];
				d += diff*diff;
				diff = (int)facspj32->data[10] - facspi32->data[10];
				d += diff*diff;
				diff = (int)facspj32->data[11] - facspi32->data[11];
				d += diff*diff;
				diff = (int)facspj32->data[12] - facspi32->data[12];
				d += diff*diff;
				diff = (int)facspj32->data[13] - facspi32->data[13];
				d += diff*diff;
				diff = (int)facspj32->data[14] - facspi32->data[14];
				d += diff*diff;
				diff = (int)facspj32->data[15] - facspi32->data[15];
				d += diff*diff;

				#undef kDATAOFFSET
				#define kDATAOFFSET 16
				diff = (int)facspj32->data[kDATAOFFSET+0] - facspi32->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+1] - facspi32->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+2] - facspi32->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+3] - facspi32->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+4] - facspi32->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+5] - facspi32->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+6] - facspi32->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+7] - facspi32->data[kDATAOFFSET+7];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+8] - facspi32->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+9] - facspi32->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+10] - facspi32->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+11] - facspi32->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+12] - facspi32->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+13] - facspi32->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+14] - facspi32->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)facspj32->data[kDATAOFFSET+15] - facspi32->data[kDATAOFFSET+15];
				d += diff*diff;

				dsum += d;
				if (dsum > gTestDist)
					goto straight2next;
				facspj32++;
				facspi32++;
			}

			/* note: 4294967295 is the max we can get to with an unsigned int32 without overflow. */
			/* using 512 columns, it limits to a range of [0..2896]; in practice, we will not want to test such a large gTestDist... */
			/* as we loop using blocks of 32, using a range of [0..10000] would max out at 3200640032, which is fine. */
			
			if (dsum <= gTestDist) 
			{
				if (*clusterpi)
				{
					if (*clusterpj == 0)
					{
						*clusterpj = *clusterpi;
					}
					else /* i=assigned, j=assigned */
					{
						if(*clusterpi != *clusterpj)
						{
							unsigned int cluster1,cluster2;
							if(*clusterpj > *clusterpi)
							{
								cluster1 = *clusterpi;
								cluster2 = *clusterpj;
							}
							else
							{
								cluster1 = *clusterpj;
								cluster2 = *clusterpi;
							}
							pthread_mutex_lock(&mergingrequestmutex);
							InsertMergeRequest(cluster1,cluster2);
							pthread_mutex_unlock(&mergingrequestmutex);
						}
					}
					
				}
				else
				{
					if (*clusterpj == 0)   /* i=not yes assigned, j=not yet assigned */
					{
						*clusterpi = ++clustercnt;
						*clusterpj = clustercnt;
					}
					else /* i=not yes assigned, j=assigned */
					{
						*clusterpi = *clusterpj;
					}
				}
			}
	straight2next:		
			facspj++;
			clusterpj++;
		}
		facspi++;
		clusterpi++;
	}
	
} /* computesim */
#else
static void computesim(FACSDATA *facs, unsigned int *clusterid, CPU *chunk)
{
	FACSDATA *facspi;
	FACSDATA *facspj;
	unsigned int *clusterpi;
	unsigned int *clusterpj;
	unsigned int i,j;
	unsigned int starti,startj;
	unsigned int lasti,lastj;
	unsigned  int actualstartj;


	starti = chunk->ii;
	startj = chunk->jj;
	lasti = chunk->iilast;
	lastj = chunk->jjlast;


		if (chunk->ii != chunk->jj)
		{
				unsigned short valii,valjj;

				facspi = &facs[lasti-1];
				facspj = &facs[startj];
				valii = facspi->data[sortkey];
				valjj = facspj->data[sortkey];
				if (valjj > valii)
				{
					unsigned  int mindist = valjj-valii;
					mindist *= mindist;
					if (mindist > gTestDist)
					{
						goto skipthischunk;
					}
				}
		}
	
		facspi = &facs[starti];
		clusterpi = &clusterid[starti];
		for(i = starti;  i< lasti; i++)
		{
			if (startj <= i)
				actualstartj = i+1;
			else 
				actualstartj = startj;
			facspj = &facs[actualstartj];
			clusterpj = &clusterid[actualstartj];
			for(j = actualstartj;  j< lastj; j++)
			{
				int diff;
				unsigned int d;

				if ((*clusterpj) && (*clusterpj == *clusterpi))
					goto straight2next;


				diff = (int)facspj->data[0] - facspi->data[0];
				d = diff*diff;
				diff = (int)facspj->data[1] - facspi->data[1];
				d += diff*diff;
				diff = (int)facspj->data[2] - facspi->data[2];
				d += diff*diff;
				diff = (int)facspj->data[3] - facspi->data[3];
				d += diff*diff;
	#if defined(COLUMNS_8) || defined(COLUMNS_12) || defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[4] - facspi->data[4];
				d += diff*diff;
				diff = (int)facspj->data[5] - facspi->data[5];
				d += diff*diff;
				diff = (int)facspj->data[6] - facspi->data[6];
				d += diff*diff;
				diff = (int)facspj->data[7] - facspi->data[7];
				d += diff*diff;
	#endif
	#if defined(COLUMNS_12) || defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[8] - facspi->data[8];
				d += diff*diff;
				diff = (int)facspj->data[9] - facspi->data[9];
				d += diff*diff;
				diff = (int)facspj->data[10] - facspi->data[10];
				d += diff*diff;
				diff = (int)facspj->data[11] - facspi->data[11];
				d += diff*diff;
	#endif
	#if defined(COLUMNS_16) || defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[12] - facspi->data[12];
				d += diff*diff;
				diff = (int)facspj->data[13] - facspi->data[13];
				d += diff*diff;
				diff = (int)facspj->data[14] - facspi->data[14];
				d += diff*diff;
				diff = (int)facspj->data[15] - facspi->data[15];
				d += diff*diff;
	#endif 
	#if  defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
	#if  defined(COLUMNS_24) || defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				/* computation done for the first 16 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;
	#endif
				#undef kDATAOFFSET
				#define kDATAOFFSET 16
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+4] - facspi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+5] - facspi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+6] - facspi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+7] - facspi->data[kDATAOFFSET+7];
				d += diff*diff;
	#if defined(COLUMNS_32) || defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[kDATAOFFSET+8] - facspi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+9] - facspi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+10] - facspi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+11] - facspi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+12] - facspi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+13] - facspi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+14] - facspi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+15] - facspi->data[kDATAOFFSET+15];
				d += diff*diff;
	#endif
	#endif
	#if   defined(COLUMNS_40) || defined(COLUMNS_48) || defined(COLUMNS_52)
				/* computation done for the first 32 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;

				#undef kDATAOFFSET
				#define kDATAOFFSET 32
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+4] - facspi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+5] - facspi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+6] - facspi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+7] - facspi->data[kDATAOFFSET+7];
				d += diff*diff;
	  #if    defined(COLUMNS_48) || defined(COLUMNS_52)
				diff = (int)facspj->data[kDATAOFFSET+8] - facspi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+9] - facspi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+10] - facspi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+11] - facspi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+12] - facspi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+13] - facspi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+14] - facspi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+15] - facspi->data[kDATAOFFSET+15];
				d += diff*diff;
	  #endif
	#endif
	#if   defined(COLUMNS_52)
				#undef kDATAOFFSET
				#define kDATAOFFSET 48
				diff = (int)facspj->data[kDATAOFFSET+0] - facspi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+1] - facspi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+2] - facspi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)facspj->data[kDATAOFFSET+3] - facspi->data[kDATAOFFSET+3];
				d += diff*diff;
	#endif
	
				if (d <= gTestDist) 
				{
					if (*clusterpi)
					{
						if (*clusterpj == 0)
						{
							*clusterpj = *clusterpi;
						}
						else /* i=assigned, j=assigned */
						{
							unsigned int cluster1,cluster2;
							if(*clusterpj > *clusterpi)
							{
								cluster1 = *clusterpi;
								cluster2 = *clusterpj;
							}
							else
							{
								cluster1 = *clusterpj;
								cluster2 = *clusterpi;
							}
							InsertMergeRequest(cluster1,cluster2);
						}
						
					}
					else
					{
						if (*clusterpj == 0)   /* i=not yes assigned, j=not yet assigned */
						{
							pthread_mutex_lock(&clustercntmutex);
							*clusterpi = ++clustercnt;
							*clusterpj = clustercnt;
							pthread_mutex_unlock(&clustercntmutex);
						}
						else /* i=not yes assigned, j=assigned */
						{
							*clusterpi = *clusterpj;
						}
					}
				}
		straight2next:		
				facspj++;
				clusterpj++;
			}
			facspi++;
			clusterpi++;
		}
skipthischunk: ;
	
} /* computesim */
#endif

/* ------------------------------------------------------------------------------------ */
static void DoComputingSlave(FACSDATA *facsdata,unsigned int *clusterid,int idproc,unsigned int initialClusterCnt)
{
			CPU	  cpudata;
			MPI_Request mpireq;

			/* determine the first value to use to start recording new clusterids for this processor */
			clustercnt = (idproc*kStartLocalCluster+initialClusterCnt);
			
			do
			{
				MPI_Recv(&cpudata, 4, MPI_INT,  0, kWhichBlocksToCompute, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (cpudata.ii != kNoMoreBlocks)
				{
					int sndcnt,rcvcnt;

					/* ------ update clusterid for each data of that might be touched by the process */

					rcvcnt =  (cpudata.iilast-cpudata.ii);
					MPI_Recv(&clusterid[cpudata.ii],rcvcnt, MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if (cpudata.jj != cpudata.ii)
					{
						rcvcnt =  (cpudata.jjlast-cpudata.jj);
						MPI_Recv(&clusterid[cpudata.jj],rcvcnt, MPI_INT,  0, kClusterMsg2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}

					/* ------ do the heavy computation */

					{
						
						if ((cpudata.jj != cpudata.ii))
						{
							pthread_t thread;
							EXECUTIONPLAN ep;
							CPU			cpudatasection;
							unsigned int tmpl;
							unsigned int where;

							ep.facsdata = facsdata;
							ep.clusterid = clusterid;
							
							// block 1 against 1
							ep.chunk.ii = ep.chunk.iilast = cpudata.ii;
							ep.chunk.iilast += ((cpudata.iilast - cpudata.ii) >> 1);
							ep.chunk.jj = ep.chunk.jjlast = cpudata.jj;
							ep.chunk.jjlast += ((cpudata.jjlast - cpudata.jj) >> 1);

							thread_mergerequestcnt = 0;
							thread_mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test thread_mergerequestcnt == 0 in thread_InsertMergeRequest */
							if (pthread_create (&thread, NULL, &computesim_funcion, &ep)) 
								printf("Error: Failed creating thread\n");

							// block 2 against 2
							cpudatasection.ii = cpudata.ii + ((cpudata.iilast - cpudata.ii) >> 1) ;
							cpudatasection.iilast = cpudata.iilast;
							cpudatasection.jj = cpudata.jj + ((cpudata.jjlast - cpudata.jj) >> 1) ;
							cpudatasection.jjlast = cpudata.jjlast;
							computesim(facsdata,clusterid,&cpudatasection);

							if (pthread_join (thread, NULL))
								printf("Error: Failed pthread_join\n");
							where = 0;
							for (tmpl = 0; tmpl < thread_mergerequestcnt; tmpl++)
							{
								InsertMergeRequestWhere(thread_mergerequest[tmpl].cluster1,thread_mergerequest[tmpl].cluster2,&where);
							}

							// block 1 against 2
							ep.chunk.jj = ep.chunk.jjlast;
							ep.chunk.jjlast = cpudata.jjlast;

							thread_mergerequestcnt = 0;
							if (pthread_create (&thread, NULL, &computesim_funcion, &ep)) 
								printf("Error: Failed creating thread\n");

							// block 2 against 1
							cpudatasection.jjlast = cpudatasection.jj;
							cpudatasection.jj = cpudata.jj;
							computesim(facsdata,clusterid,&cpudatasection);

							if (pthread_join (thread, NULL))
								printf("Error: Failed pthread_join\n");

							where = 0;
							for (tmpl = 0; tmpl < thread_mergerequestcnt; tmpl++)
							{
								InsertMergeRequestWhere(thread_mergerequest[tmpl].cluster1,thread_mergerequest[tmpl].cluster2,&where);
							}

						}
						else 
							computesim(facsdata,clusterid,&cpudata);

					}


					/* ------ send back updated data */

					MPI_Isend(&idproc, 1, MPI_INT,  0, kCPUdoneMsg, MPI_COMM_WORLD,&mpireq);
					sndcnt = cpudata.iilast-cpudata.ii;
					MPI_Send(&clusterid[cpudata.ii], sndcnt, MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD);
					if (cpudata.jj != cpudata.ii)
					{
						sndcnt = cpudata.jjlast-cpudata.jj;
						MPI_Send(&clusterid[cpudata.jj], sndcnt, MPI_INT,  0, kClusterMsg2, MPI_COMM_WORLD);
					}
				}
				else /* send the final count of "new clusters" allocated by this proc. */
				{
					if (cpudata.jj != kNoMoreBlocks)
					{
						MPI_Send((void *)&clustercnt, 1, MPI_INT,  0, kFinalCntRequest, MPI_COMM_WORLD);
					}
				}
			} while (cpudata.ii != kNoMoreBlocks);

			/* ------- JOIN mergerequest list (cpus obtain the list from an other cpu until everything is merged) ------- */

			do
			{
				JOINREQUEST joinRequest;

				MPI_Recv(&joinRequest, 2, MPI_INT,  0, kJoinListRequest, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
			if (joinRequest.getFromCPU == 0)
				break;
				
				if (joinRequest.getFromCPU == idproc) /*  we are  the one to send */
				{
					MPI_Send(&mergerequestcnt, 1, MPI_INT,  joinRequest.sendToCPU, kFinalMergeRequestCntRequest, MPI_COMM_WORLD);
					if (mergerequestcnt > 0)
						MPI_Send(&mergerequest, mergerequestcnt*2, MPI_INT,  joinRequest.sendToCPU, kSendFinalMergeRequestRequest, MPI_COMM_WORLD);
				}
				else /* we are the one receiving */
				{
					MERGECLUSTER cpumergerequest[kMaxMergeRequests];
					unsigned int cpumergerequestcnt;
					unsigned int where;
					int tmpl;
					MPI_Recv(&cpumergerequestcnt, 1, MPI_INT,  joinRequest.getFromCPU, kFinalMergeRequestCntRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					if (cpumergerequestcnt > 0)
						MPI_Recv(&cpumergerequest, cpumergerequestcnt*2, MPI_INT,  joinRequest.getFromCPU, kSendFinalMergeRequestRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					where = 0;
					for (tmpl = 0; tmpl < cpumergerequestcnt; tmpl++)
					{
						InsertMergeRequestWhere(cpumergerequest[tmpl].cluster1,cpumergerequest[tmpl].cluster2,&where);
					}
					MPI_Send(&mergerequestcnt,1, MPI_INT,  0, kJoinListRequestDone, MPI_COMM_WORLD);
				}
			} while(1);
			
			if (idproc == 1) /* send final list to master */
			{
				MPI_Send(&mergerequestcnt, 1, MPI_INT,  0, kFinalMergeRequestCntRequest, MPI_COMM_WORLD);
				if (mergerequestcnt > 0)
					MPI_Send(&mergerequest, mergerequestcnt*2, MPI_INT,  0, kSendFinalMergeRequestRequest, MPI_COMM_WORLD);
			}

	
} /* DoComputingSlave */
/* ------------------------------------------------------------------------------------ */


/* function repeatadly called only by the master to identify a suitable computing chunk to asign to an available slave */
static unsigned int AssignChunk(unsigned int *clusterid,CHUNK *chunk, CPU *cpu, unsigned int nproc)
{
	unsigned int i;
	int sndcnt;
	MPI_Request mpireq;
	unsigned int candidate = 0;

	/* test if computing in already in progress somewhere for one of those blocks */
	// start from last to first proc, as according to lsf, first proc will have lowest load and we submit to the last identified avail proc.
	for (i = (nproc-1); i>0; i--)
	{
		if ((cpu[i].ii == chunk->ii) || (cpu[i].jj == chunk->jj) || (cpu[i].ii == chunk->jj) || (cpu[i].jj == chunk->ii))
			return(0); 

		if (cpu[i].ii == kCPU_availaible)
			candidate = i;
	}

	cpu[candidate].ii = chunk->ii;
	cpu[candidate].jj = chunk->jj;
	cpu[candidate].iilast = chunk->iilast;
	cpu[candidate].jjlast = chunk->jjlast;
	MPI_Isend(&cpu[candidate], 4, MPI_INT,  candidate, kWhichBlocksToCompute, MPI_COMM_WORLD,&mpireq);

	sndcnt = cpu[candidate].iilast-cpu[candidate].ii;
	MPI_Send(&clusterid[cpu[candidate].ii], sndcnt, MPI_INT,  candidate, kClusterMsg1, MPI_COMM_WORLD);
	if (cpu[candidate].ii != cpu[candidate].jj)
	{
		sndcnt = cpu[candidate].jjlast-cpu[candidate].jj;
		MPI_Send(&clusterid[cpu[candidate].jj], sndcnt, MPI_INT,  candidate, kClusterMsg2, MPI_COMM_WORLD);
	}
	chunk->status = kChunkStatusComputing;

	return(candidate);
	
} /* AssignChunk */
/* ------------------------------------------------------------------------------------ */
int main (int argc, char **argv)
{	

	char	fn[kMaxFilename];
	char	ofn[kMaxFilename];
	unsigned int colcnt;
	char	*version="VERSION 1.0; 2019-12-26";
	FACSNAME *facsname = NULL;
	FACSDATA *facsdata = NULL;
	FILE	 *f=NULL;
	unsigned int loaded;
	float distcutoff;
	float distcutoffincreasestep;
	float lastdistcutoff;
	float bestdistcutoff;

	int idproc, nproc;
	int verbose;
	int goOnEvenIfClusterCntDecreases;
	unsigned int rowcnt;
	float pctEventsToKeepCluster;
	unsigned int cntcutoff = 0;
	int c;
	unsigned int desiredBlockSize = 0;
	float stopWhenPctAssigned;
	struct timeval ts;
	unsigned int printClusterStatus = 0;
	unsigned int assignUnassigned = 0;
	unsigned int assignLeftover = 0;
	
	/* must be first instruction */
    if (MPI_Init(&argc, &argv))
		return(1);

	gettimeofday(&ts, NULL); 

	/* --------- process arguments */

	colcnt = 0;
	rowcnt = 0;

	fn[0] = 0;
	ofn[0] = 0;
	distcutoff = 0.0;
	distcutoffincreasestep = 0.5;
	lastdistcutoff = -1.0;
	bestdistcutoff = 0.0;
	pctEventsToKeepCluster = 0.5;
	goOnEvenIfClusterCntDecreases = 0;
	verbose = 0;
	stopWhenPctAssigned = 95.0;
	opterr = 0;
	while ((c = getopt (argc, argv, "i:o:f:l:s:k:n:p:b:v:gMUL")) != -1)
	switch (c)
	{
      case 'i':
			strcpy(fn,optarg);
        break;
      
	  case 'o':
			strcpy(ofn,optarg);
        break;
      
	  case 'f':
			sscanf(optarg,"%f",&distcutoff);
        break;

	  case 'l':
			sscanf(optarg,"%f",&lastdistcutoff);
        break;

	  case 's':
			sscanf(optarg,"%f",&distcutoffincreasestep);
			if (distcutoffincreasestep < 0.0)
				distcutoffincreasestep = -distcutoffincreasestep;
        break;

	  case 'g':
			goOnEvenIfClusterCntDecreases = 1;
        break;

	  case 'k':
			sscanf(optarg,"%f",&pctEventsToKeepCluster);
        break;

	  case 'n':
			sscanf(optarg,"%u",&cntcutoff);
        break;

	  case 'p':
			sscanf(optarg,"%f",&stopWhenPctAssigned);
        break;
		
	  case 'b':
			sscanf(optarg,"%u",&desiredBlockSize);
		break;
			
	  case 'M':
			printClusterStatus = 1;
		break;

	  case 'U':
			assignUnassigned = 1;
		break;

	  case 'L':
			assignLeftover = 1;
		break;

	  case 'v':
			sscanf(optarg,"%d",&verbose);
        break;

	}
	if (ofn[0] == 0)
		strcpy(ofn,fn);

	if ((fn[0] == 0) || (distcutoff < 0.00001))
	{
		printf("usage:\n\n");
		printf("dclust -i InputFile -f FirstDistanceCutoff [-l LastDistanceCutoff [-s Step] [-g]] [-o OutputFile] [-k PctEventsToKeepCluster | -n numEventsToKeepCluster] [-p pctAssigned] [ -v level]\n\n");
		printf("       -i InputFile              : dselect binary output file.\n");
		printf("       -f FirstDistanceCutoff    : First Floating point cutoff value used to place events in the same cluster.\n");
		printf("       -l LastDistanceCutoff     : Last Distance cutoff to test. Defaults is the same as DistanceCutoff.\n");
		printf("       -s Step                   : Floating point increment of DistanceCutoff to test. Default is %f\n",distcutoffincreasestep);
		printf("                                   Must be positive. To scan decreasing distances, use a DistanceCutoff > LastDistanceCutoff\n");
		printf("       -g                        : go on with scanning even if the number of clusters retained decreases. Default is to stop.\n");
		printf("                                   Must be positive. To scan decreasing distances, use a DistanceCutoff > LastDistanceCutoff\n");
		printf("       -p pctAssigned            : Stop sampling as soon as pctAssigned events have been assigned. Defaults to %f %%\n",stopWhenPctAssigned);
		printf("       -o OutputFile             : Rootname for the output files. Various extensions will be added.\n");
		printf("                                   Default is same name as InputFile\n");
		printf("       -k PctEventsToKeepCluster : Keep only clusters with more than Percent Events. Default is 0.5\n");
		printf("       -n numEventsToKeepCluster : Keep only clusters with at least this number of events.\n");
		printf("                                   Default is computed from the -k option (0.5 percent of input)\n");
		printf("       -M                        : Report cluster Merging history\n");
		printf("       -U                        : assign Unassigned to discovered clusters\n");
		printf("       -L                        : assign Leftover (see dselect) to discovered clusters\n");
		printf("       -v level                  : specifies the verbose level; default is 0.\n\n");
		printf("VERSION\n");
		printf("\n%s\n",version);
		printf("Author:  Nicolas Guex; 2008-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
		printf("CONTACT: Nicolas.Guex@unil.ch\n");
		printf("SEE ALSO\n");
		printf("      dselect cextract\n\n");
		return(1);
	}

	/* --------- initialize MPI */
	nproc=0; 
	
	if (MPI_Comm_rank(MPI_COMM_WORLD, &idproc))
		return(1);

	if (MPI_Comm_size(MPI_COMM_WORLD, &nproc))
		return(1);

	if (nproc >= kMaxCPU)
	{
		printf("Software has been compiled to run at most on %d cpus\n",(kMaxCPU-1));
		return(1);
	}
	

	/* --------- process */
	if (lastdistcutoff < 0.0)
		lastdistcutoff = distcutoff; 

	if (idproc == 0)
	{
		printf("LOG: %s; max_rows: %d; max_columns: %d; max_input_value: %d\n",version,kMAXEVENTS,kMaxInputCol,kMAX_ALLOWED_INPUT_VALUE);
		printf("LOG: first_dist: %f; last_dist: %f; dist_step: %f\n",distcutoff,lastdistcutoff,distcutoffincreasestep);
		if (verbose > 0)
		{
			printf("LOG:InputFile=%s\n",fn);
			printf("LOG:OutputFileRootname=%s\n",ofn);
		}
	}

	if (idproc == 0)
	{
		f = fopen(fn,"r");
		if (!f)
		{
			printf("LOG:Cannot Open InputFile %s\n",fn);
			rowcnt = 0;	/* signal slaves that they can bail */
			MPI_Bcast (&rowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&rowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&rowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}		

	{
		unsigned int ii;
		unsigned int jj;
		unsigned int *clusterid = NULL;
		
		int mrg;
		int *cnp;
		unsigned int loadEveryNsample;
		struct timeval te;
		unsigned short key;

		if (idproc == 0)
		{
			if (dclustFileReadHeader(f,idproc,&rowcnt,&colcnt,&key,&loadEveryNsample) != 0)
				goto abort;
			MPI_Bcast (&rowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&colcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&key, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
		else
		{
			loadEveryNsample = 0; /* unused for slave, stop compiler warnings */
			MPI_Bcast (&rowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&colcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&key, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (rowcnt == 0)
				goto abort;
		}

		sortkey = key;

		/* --------- allocate memory */
		facsdata = calloc(rowcnt,sizeof(FACSDATA));
		if (!facsdata)
			goto abort;

		if (idproc == 0)
		{
			facsname = calloc(rowcnt,sizeof(FACSNAME));
			if (!facsname)
				goto abort;
		}
		clusterid = calloc(rowcnt,sizeof(MPI_INT));

		/* should test if calloc worked */
		if (idproc == 0)
		{
			unsigned int tosend = rowcnt;
			unsigned int sendfrom = 0;
			loaded = loaddata(f,facsname,facsdata,rowcnt,colcnt,idproc);
			if (loaded == 0)
			{
				printf("LOG:error (no input data)\n");
				goto abort;
			}
			while(tosend > 0)
			{
	                    unsigned int smallchunk;
			    if (tosend <= 1000000)
			    {
					smallchunk = tosend;
					tosend = 0;
			    }
			    else
			    {
					smallchunk = 1000000;
					tosend -= 1000000;
 			    }
			    MPI_Bcast (&facsdata[sendfrom], smallchunk*sizeof(FACSDATA), MPI_CHAR, 0, MPI_COMM_WORLD);
			    sendfrom += 1000000;
			}
			MPI_Bcast(&clusterid[0], rowcnt, MPI_INT,  0, MPI_COMM_WORLD);
		}
		else
		{
			unsigned int tosend = rowcnt;
			unsigned int sendfrom = 0;

			loaded=rowcnt;  // ???????? really useful for slave ?????
			while(tosend > 0)
			{
				unsigned int smallchunk;    
				if (tosend <= 1000000)
				{   
					smallchunk = tosend;
					tosend = 0;
				}   
				else
				{   
				   smallchunk = 1000000;
					tosend -= 1000000;
				}
				MPI_Bcast (&facsdata[sendfrom], smallchunk*sizeof(FACSDATA), MPI_CHAR, 0, MPI_COMM_WORLD);
				sendfrom += 1000000;
			}
			MPI_Bcast(&clusterid[0], rowcnt, MPI_INT,  0, MPI_COMM_WORLD);
		}
	
		if (cntcutoff > 0)
		{
			if (cntcutoff < colcnt)
			{
				printf("LOG:Warning:pctEventsToKeepCluster lower than column count (%u < %u)\n",cntcutoff,colcnt);
			}
			pctEventsToKeepCluster = (float)cntcutoff / (float)(rowcnt)*100.0;
		}
		else /* compute cntctoff */
		{
			cntcutoff = (unsigned int)(((float)(rowcnt) / 100.0 * pctEventsToKeepCluster ) );
			if (cntcutoff < colcnt)
			{
				printf("LOG:Warning:pctEventsToKeepCluster lower than column count (%u < %u)\n",cntcutoff,colcnt);
			}
		}
		
		
		if (idproc == 0)  /* ---------------- master node ------------- */
		{
			printf("LOG:Loaded %u rows of %u columns from %s\n",rowcnt,colcnt,fn);
			printf("LOG:PctEventsToKeepCluster=%.3f (%u events)\n",pctEventsToKeepCluster,cntcutoff);
			printf("DBH:totseconds,cpus,loadEveryNsample,distcutoff,loaded,assigned,unassigned,pctassigned,pctunassigned,clustercnt,trimclustercnt\n");
			fclose(f);
		}
		gTestDist = (unsigned int)(distcutoff*distcutoff*colcnt);

		if (idproc == 0)  /* ---------------- master node ------------- */
		{
			CHUNK *chunk;
			CPU	  cpu[kMaxCPU];
			unsigned int whichcpu;
			unsigned int alldone = 1;  // will be initialized, ignore compiler whining.
			unsigned int submitted;
			unsigned int chunkcnt;
			unsigned int chunckcnt;
			unsigned  int unassigned;
			unsigned int processingBlockSize = 131072;
			int trimmedclustercnt;
			int initialClusterCnt;
			int highesttrimmedclustercnt = -1;
			unsigned int passcnt = 0;
			unsigned int  isMergingPreexistingClusters;
			float distOfLastClusterIndices = 0.0;
			STATS stats[2];

			clusterhistory = malloc(kMaxCluster*sizeof(CLUSTERHISTORY));
			if (!clusterhistory)
			{
				fprintf(stderr,"LOG: ERROR: not enough memory\n");
				goto abort;
			}

			if (desiredBlockSize == 0)
			{
				
				/* arrange to keep each slave node busy with at least about 100 computations, but do not go below blocksize of 256 events */
				do 
				{
					chunkcnt = ((loaded/processingBlockSize+2)*(loaded/processingBlockSize+2))/2;
					processingBlockSize >>= 1;
				} while( ((chunkcnt / (nproc-1)) < 100) && (processingBlockSize > 128) );
				processingBlockSize <<= 1;
			}
			else
			{
				processingBlockSize = desiredBlockSize;
				chunkcnt = ((loaded/processingBlockSize+2)*(loaded/processingBlockSize+2))/2;
			}

			chunk = calloc(chunkcnt,sizeof(CHUNK));
			if (!chunk)
			{
				printf("LOG:Not enough memory to allocate %u chunks; recompile with larger processingBlockSize\n",chunkcnt);	
				cpu[0].ii = kNoMoreBlocks;
				cpu[0].jj = kNoMoreBlocks;
				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&cpu[0], 4, MPI_INT,  ii, kWhichBlocksToCompute, MPI_COMM_WORLD);
				goto abort;
			}				


			
			stats[0].dist = 0.0;
			stats[0].rawClustersCnt = -1;
			stats[0].trimmedClustersCnt = -1;
			stats[0].pctAssigned = 0.0;
			
repeatWithNewDist:
			stats[1].dist = distcutoff;
			stats[1].rawClustersCnt = -1;
			stats[1].trimmedClustersCnt = -1;
			stats[1].pctAssigned = 0.0;
			isMergingPreexistingClusters = 0;

			printf("LOG:DistanceCutoff=%.3f\n",distcutoff);
			mergerequestcnt = 0;
			mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test mergerequestcnt == 0 in InsertMergeRequest */
			trimmedclustercnt = -1;
			for (ii = 1; ii<nproc; ii++)
			{
				cpu[ii].ii = kCPU_availaible;
				cpu[ii].jj = kCPU_availaible;
			}
			chunckcnt = 0;
			for (ii = 0; ii < loaded; ii += processingBlockSize)
			for (jj = ii; jj < loaded; jj += processingBlockSize)
			{
				FACSDATA *facspi;
				FACSDATA *facspj;
				unsigned short valii,valjj;
				
				chunk[chunckcnt].ii = ii;
				chunk[chunckcnt].jj = jj;
				if ((ii+processingBlockSize) < loaded)
					chunk[chunckcnt].iilast = ii+processingBlockSize; 
				else
					chunk[chunckcnt].iilast  = loaded;
				if ((jj+processingBlockSize) < loaded)
					chunk[chunckcnt].jjlast = jj+processingBlockSize;
				else
					chunk[chunckcnt].jjlast  = loaded;
					
				if (ii != jj)
				{
					facspi = &facsdata[chunk[chunckcnt].iilast-1];
					facspj = &facsdata[jj];
					valii = facspi->data[key];
					valjj = facspj->data[key];
					if (valjj > valii)
					{
						unsigned  int mindist = valjj-valii;
						mindist *= mindist;
						if (mindist > gTestDist)
						{
							continue;
						}
					}
				}
				chunk[chunckcnt].status = kChunkStatusToDo;
				chunckcnt++;
			}

			submitted = 0;
			fflush(stdout);
			do
			{				
				/* if at least one cpu idle, try to submit as many jobs as possible */
				if (submitted < (nproc-1))
				{
					alldone = 1;
					for (ii = 0; ii < chunckcnt; ii++)
					{
						if (chunk[ii].status == kChunkStatusToDo)
						{
							alldone = 0;
							if (AssignChunk(clusterid,&chunk[ii],cpu,nproc))
							{
									if (verbose > 1)
									{
										printf("LOG:Sent computation Request for chunk %8u / %8u\n",ii,chunckcnt);
										fflush(stdout);
									}
									submitted++;
									if (submitted == (nproc-1))
										break;
							}
						}
					}
				}
				/* wait until one of the slave cpu is done and update available list accordingly */
				if (submitted)
				{
					int rcvcnt;

					MPI_Recv(&whichcpu, 1, MPI_INT,  MPI_ANY_SOURCE, kCPUdoneMsg, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					rcvcnt =  (cpu[whichcpu].iilast-cpu[whichcpu].ii);
					MPI_Recv(&clusterid[cpu[whichcpu].ii],rcvcnt, MPI_INT,  whichcpu, kClusterMsg1, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					if (cpu[whichcpu].jj != cpu[whichcpu].ii)
					{
						rcvcnt =  (cpu[whichcpu].jjlast-cpu[whichcpu].jj);
						MPI_Recv(&clusterid[cpu[whichcpu].jj],rcvcnt, MPI_INT,  whichcpu, kClusterMsg2, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					}
					cpu[whichcpu].ii = cpu[whichcpu].jj = kCPU_availaible;
					submitted--;
				}
			} while ((submitted > 0) || (!alldone));

			if (verbose > 1)
			{
				printf("LOG:Master Starts Collecting Results\n");
				fflush(stdout);
			}


			
			/* signal to all nodes that they should clean up */
			/* collect cluster number assigned by each proc and adjust clusters from 1..clustercnt */
			cpu[0].ii = kNoMoreBlocks;
			cpu[0].jj = 0;
			for (ii = 1; ii<nproc; ii++)
			{
				int finalCPUcnt;

				MPI_Send(&cpu[0], 4, MPI_INT,  ii, kWhichBlocksToCompute, MPI_COMM_WORLD);
				MPI_Recv(&finalCPUcnt, 1, MPI_INT,  ii, kFinalCntRequest, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
				finalCPUcnt -= (ii*kStartLocalCluster);
				if (verbose > 2)
				{
					printf("LOG:Final Cluster Cnt For CPU %3u = %6d\n",ii,finalCPUcnt);
					fflush(stdout);
				}
				clustersnum[ii] = calloc((finalCPUcnt+1),sizeof(int));
				if (!clustersnum[ii])
				{
					printf("LOG:Not enough memory to allocate cluster ID %u\n",ii);	
					goto abort;
				}				
				cnp = clustersnum[ii];
				*cnp = finalCPUcnt;

				if (ii == 1)
				{
					clustercnt = finalCPUcnt;
				}
				else
				{
					clustercnt += finalCPUcnt;
				}
			}

			DoMergeLists(nproc,verbose);
			MPI_Recv(&mergerequestcnt, 1, MPI_INT,  1, kFinalMergeRequestCntRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if (mergerequestcnt > 0)
				MPI_Recv(&mergerequest, mergerequestcnt*2, MPI_INT,  1, kSendFinalMergeRequestRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			if (verbose > 1)
			{
				printf("LOG:Master received %8u mergerequests resulting in %u clusters\n",mergerequestcnt,(clustercnt-mergerequestcnt));
				fflush(stdout);
			}

			gettimeofday(&te, NULL);
			isMergingPreexistingClusters = ProcessMergeRequests(clusterid,loaded,mergerequestcnt,nproc,stats[0].trimmedClustersCnt,stats[0].dist,passcnt);
			gettimeofday(&te, NULL);


			for (mrg = mergerequestcnt-1; mrg >= 0 ; mrg--)
			{
				removeclustersnum(mergerequest[mrg].cluster2);
			}

			/* Discarding clusters with too few events */
			trimmedclustercnt = RemoveSmallClusters(clusterid,loaded,nproc,cntcutoff);
			stats[1].trimmedClustersCnt = trimmedclustercnt;
			printf("LOG: %12d Clusters retained with at least %9.3f %% events (at least %u events)\n",trimmedclustercnt,pctEventsToKeepCluster,cntcutoff);
			fflush(stdout);

			UpdateClusterHistory(passcnt,trimmedclustercnt,stats[0].trimmedClustersCnt,distcutoff,verbose);
			AdjustClustersID(clusterid,loaded,nproc,verbose,trimmedclustercnt,&initialClusterCnt);
			unassigned = loaded;
			if (trimmedclustercnt == 0)
			{
				printf("LOG: %12d Assigned    (%5.1f %%)\n",0,0.0);	
				printf("LOG: %12u Unassigned  (%5.1f %%)\n",unassigned,100.0);	
				fflush(stdout);
			}
			else  /* write binary results file */
			{
				if (trimmedclustercnt > highesttrimmedclustercnt)
				{
					highesttrimmedclustercnt = trimmedclustercnt;
					bestdistcutoff = distcutoff;
				}
				else if ((trimmedclustercnt == highesttrimmedclustercnt) && (lastdistcutoff >= distcutoff)) /* if we go up, and cluster still as good, overwrite */
				{
					bestdistcutoff = distcutoff;
				}
				
				{

					if (verbose > 1)
					{
						printf("LOG:Adjusting Cluster Numbers\n");
						fflush(stdout);
					}
					if (isMergingPreexistingClusters == 1)  //here should write under correct name only if merging, otherwise, save status unde "last dist" just in case we were indeed last dist...
					{
						char oldfn[kMaxFilename];
						char newfn[kMaxFilename];
						sprintf(oldfn,"%s-%.6f",ofn,0.0);
						sprintf(newfn,"%s-%.6f",ofn,distOfLastClusterIndices);
						rename(oldfn,newfn);
						WriteClusterIndices(clusterid,loaded,distcutoff,ofn);
					}
					else
					{
						WriteClusterIndices(clusterid,loaded,0.0,ofn);
					}
					distOfLastClusterIndices = distcutoff;


					{
						if (verbose > 2)
						{
							printf("LOG:Not Writing Results\n");
							fflush(stdout);
						}
						unassigned = loaded-CountAssigned(clusterid,loaded,trimmedclustercnt);
					}
					stats[1].pctAssigned = 100.0*(loaded-unassigned)/loaded;
					printf("LOG: %12u TotalEvents\n",loaded);	
					printf("LOG: %12u Assigned    (%5.1f %%)\n",(loaded-unassigned),stats[1].pctAssigned);	
					printf("LOG: %12u Unassigned  (%5.1f %%)\n",unassigned,100.0-stats[1].pctAssigned);	
					fflush(stdout);
				}

			}
						
			for (ii = 1; ii<nproc; ii++)
			{
				if (clustersnum[ii])
					free(clustersnum[ii]);
			}

			gettimeofday(&te, NULL);
			printf("DBV:%d,%d,%u,%.3f,%u,%u,%u,%.3f,%.3f,%u,%d\n",((int)te.tv_sec-(int)ts.tv_sec),nproc,loadEveryNsample,distcutoff,loaded,loaded-(unassigned),(unassigned),stats[1].pctAssigned,(100.0 - stats[1].pctAssigned),(clustercnt-mergerequestcnt),trimmedclustercnt);
			fflush(stdout);

			/* test if should keep scanning */
			ii = 0;
			if (lastdistcutoff == distcutoff) 
				ii = 1;
			else
			{
				if (lastdistcutoff < distcutoff)  /* we go down */
				{
					if (trimmedclustercnt == 0)  /* no need to go further down!!! */
						ii = 1;
					else
					{
						distcutoff -= distcutoffincreasestep;
						if (distcutoff < lastdistcutoff)  // all done.
							ii = 1;
						if ((goOnEvenIfClusterCntDecreases == 0) && (trimmedclustercnt < highesttrimmedclustercnt))
						{
							printf("LOG: Stopping. (number of clusters decreases)\n");	
							ii = 1;
						}
					}
				}
				else  /* we go up */
				{
					/*we stop if number of untrimmed clusters becomes 1 and at least one trimmed cluster exists,
					 or at least stopWhenPctAssigned % of events have been assigned to a retained cluster */
					if (  (((clustercnt-mergerequestcnt) == 1) && (trimmedclustercnt >= 1) && (passcnt > 0)) || (stats[1].pctAssigned >= stopWhenPctAssigned)) 
						ii = 1;
					else
					{
						/* if less than 0.1% change between two tested distances and only 1 cluster left, with over 50% of events assigned, double sampling step. */
						if ((  (trimmedclustercnt == 1) && (stats[1].pctAssigned > 50.0) && ((stats[1].pctAssigned - stats[0].pctAssigned) <= 0.1) ) 
						|| ((stats[1].pctAssigned >= 99.0) && (stats[1].trimmedClustersCnt == stats[0].trimmedClustersCnt) && ((stats[1].pctAssigned - stats[0].pctAssigned) <= 0.1)))
						{
							distcutoffincreasestep *= 2.0;
							printf("LOG: Increasing distance cutoff sampling step to %.3f\n",distcutoffincreasestep);	
						}
						distcutoff += distcutoffincreasestep;
						if (((clustercnt-mergerequestcnt) == 0) && (stats[0].rawClustersCnt == 0))
						{
							distcutoff += distcutoffincreasestep;
						}
						stats[1].rawClustersCnt = (clustercnt-mergerequestcnt);
						
						if (distcutoff > lastdistcutoff)  // all done.
							ii = 1;
						if ((goOnEvenIfClusterCntDecreases == 0) && (trimmedclustercnt < highesttrimmedclustercnt))
						{
							printf("LOG: Stopping. (number of clusters decreases)\n");	
							ii = 1;
						}
					}
				}
			}
			
			if (ii == 1)  // all done.
			{
				char oldfn[kMaxFilename];
				char newfn[kMaxFilename];
				sprintf(oldfn,"%s-%.6f",ofn,0.0);
				sprintf(newfn,"%s-%.6f",ofn,distOfLastClusterIndices);
				rename(oldfn,newfn);

				memset(clusterid,0,rowcnt*sizeof(MPI_INT));  // reset all clusterid
				trimmedclustercnt = SelectClusterHistory(loaded,clusterid,ofn,verbose);
				if (printClusterStatus)
					PrintClusterStatus(clusterhistory,clusterhistorycnt);
				unassigned = loaded-WriteSplitBinFile(facsname,facsdata,clusterid,loaded,colcnt,trimmedclustercnt,ofn);
				printf("LOG: %12u TotalEvents\n",loaded);	
				printf("LOG: %12u Assigned    (%5.1f %%)\n",(loaded-unassigned),100.0*(loaded-unassigned)/loaded);	
				printf("LOG: %12u Unassigned  (%5.1f %%)\n",unassigned,100.0*unassigned/loaded);	
				fflush(stdout);



				free(chunk); 
				gTestDist = 0;
				printf("LOG: Master is all done and identified a max of %d clusters at distance %.3f; notifying slaves.\n",highesttrimmedclustercnt,bestdistcutoff);	
				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&gTestDist, 1, MPI_INT,  ii,  kRepeatWithNewDistMsg, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD); 

				/* assign leftover */
				{
					unsigned int starti,lasti;
					unsigned int datachunk = 1+(loaded/nproc);
					distcutoff += distcutoffincreasestep;
					gTestDist = (unsigned int)(distcutoff*distcutoff*colcnt);

					gettimeofday(&te, NULL);
					if (assignUnassigned)
					{
						unsigned int cnt2reassign;
						cnt2reassign = FlagSequencesToReassign(clusterid,loaded,newfn);
						printf("LOG:%d sec; Distributing %u events to the %d discovered clusters.\n",((int)te.tv_sec-(int)ts.tv_sec),cnt2reassign,trimmedclustercnt);
					}
					MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD); // in fact won't be used during DistributeUnassignedToClosestCluster.
					MPI_Bcast (&trimmedclustercnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast (&clusterid[0], loaded, MPI_INT, 0, MPI_COMM_WORLD);
					lasti = 0 + datachunk;
					if (lasti > loaded)
						lasti = loaded;
					DistributeUnassignedToClosestCluster(facsdata,clusterid,loaded,colcnt,trimmedclustercnt,0,lasti);
					gettimeofday(&te, NULL);
					printf("LOG:%d sec; Collecting Results\n",((int)te.tv_sec-(int)ts.tv_sec));	
					for (ii = 1; ii<nproc; ii++)
					{
						starti = ii*datachunk;
						if (starti < loaded)
						{
							lasti = starti + datachunk;
							if (lasti > loaded)
								lasti = loaded;
							MPI_Recv(&clusterid[starti], (lasti-starti), MPI_INT,  ii, kClusterMsg1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD); 

					gettimeofday(&te, NULL);
					if (assignLeftover)
					{
						printf("LOG:%d sec; Processing leftover file\n",((int)te.tv_sec-(int)ts.tv_sec));
						DoProcessLeftoverbinaryFile(fn,facsdata,clusterid,loaded,colcnt,nproc,verbose);
					}
					else
					{
						unsigned int leftoverrowcnt = 0;
						MPI_Bcast (&leftoverrowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
					}
					MPI_Barrier(MPI_COMM_WORLD); 

					gettimeofday(&te, NULL);
					printf("LOG:%d sec; Writing Clustering Results\n",((int)te.tv_sec-(int)ts.tv_sec));	

					unassigned = loaded-WriteSplitBinFile(facsname,facsdata,clusterid,loaded,colcnt,trimmedclustercnt,ofn);

					printf("LOG: %12u TotalEvents\n",rowcnt);
					printf("LOG: %12u Assigned    (%5.1f %%)\n",(rowcnt-unassigned),100.0*(rowcnt-unassigned)/rowcnt);	
					printf("LOG: %12u Unassigned  (%5.1f %%)\n",unassigned,100.0*unassigned/rowcnt);	
				}


			}
			else
			{
				gTestDist = (unsigned int)(distcutoff*distcutoff*colcnt);

				/* if gDist increases, do as if computing node #1 had discovered valid clusters so keep resuls already valid and reduce the number of merging events */
				for(ii = 0;  ii< rowcnt; ii++)
				{
					if ((clusterid[ii] > 0))
					{
						clusterid[ii] += 1*kStartLocalCluster;
					}
					else
						clusterid[ii] = 0;
					
				}

				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&gTestDist, 1, MPI_INT,  ii,  kRepeatWithNewDistMsg, MPI_COMM_WORLD);
				MPI_Send(&initialClusterCnt, 1, MPI_INT,  1,  kInitialClusterCntMsg, MPI_COMM_WORLD); // send starting clustercount to slave node #1
				MPI_Barrier(MPI_COMM_WORLD);
				stats[0] = stats[1];
				passcnt++;

				printf("LOG:**************************************************************\n");
				goto repeatWithNewDist;
			}
		}
		else /* ---------------------- slave node  -------------------- */
		{
			unsigned int initialClusterCnt = 0;

			do
			{
				mergerequestcnt = 0;
				mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test mergerequestcnt == 0 in InsertMergeRequest */

				DoComputingSlave(facsdata,clusterid,idproc,initialClusterCnt);
				memset(clusterid,0,rowcnt*sizeof(MPI_INT));  // reset clusterid

				MPI_Recv(&gTestDist, 1, MPI_INT,0, kRepeatWithNewDistMsg,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				if ((idproc == 1) && (gTestDist > 0))
				{
					MPI_Recv(&initialClusterCnt, 1, MPI_INT,0, kInitialClusterCntMsg,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				fflush(stdout);			

				MPI_Barrier(MPI_COMM_WORLD);
				if (gTestDist == 0)
				{
					unsigned int leftoverrowcnt;
					unsigned int trimmedclustercnt;
					unsigned int starti,lasti;
					unsigned int datachunk = 1+(loaded/nproc);
					MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast (&trimmedclustercnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast (&clusterid[0], loaded, MPI_INT, 0, MPI_COMM_WORLD);
					starti = idproc*datachunk;
					if (starti < loaded)
					{
						lasti = starti + datachunk;
						if (lasti > loaded)
							lasti = loaded;
						DistributeUnassignedToClosestCluster(facsdata,clusterid,loaded,colcnt,trimmedclustercnt,starti,lasti);
						MPI_Send(&clusterid[starti], (lasti-starti), MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD);
					}
					MPI_Barrier(MPI_COMM_WORLD); 

					MPI_Bcast (&leftoverrowcnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
					fflush(stdout);
					if (leftoverrowcnt > 0)
					{
						FACSDATA *leftoverfacs=NULL;
						unsigned int *leftoverclusterid=NULL;
						
						MPI_Bcast (&clusterid[0], loaded, MPI_INT, 0, MPI_COMM_WORLD);
						fflush(stdout);
						MPI_Recv(&leftoverrowcnt, 1, MPI_INT,  0, kLeftoverDataLength, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						fflush(stdout);
						leftoverfacs = calloc(leftoverrowcnt,sizeof(FACSDATA));
						if (!leftoverfacs)
						{
							printf("LOG:CPU %d Error:Cannot Allocate Memory.\n",idproc);
							goto bail;
						}
						leftoverclusterid = calloc(leftoverrowcnt,sizeof(FACSNAME));
						if (!leftoverclusterid)
						{
							printf("LOG:CPU %d Error:Cannot Allocate Memory.\n",idproc);
							goto bail;
						}
						MPI_Recv(leftoverfacs, (int)(leftoverrowcnt*sizeof(FACSDATA)),MPI_CHAR,0, kLeftoverData,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						fflush(stdout);
						DistributeLeftoverToClosestCluster(facsdata,clusterid, loaded, colcnt, leftoverfacs, leftoverclusterid, leftoverrowcnt);
						fflush(stdout);
						MPI_Send(&leftoverclusterid[0], leftoverrowcnt, MPI_INT,  0, kLeftoverClusters, MPI_COMM_WORLD);
					bail:
						if(leftoverfacs)
							free(leftoverfacs);
						if(leftoverclusterid)
							free(leftoverclusterid);
							
					}
					MPI_Barrier(MPI_COMM_WORLD); 
					break;
				}
			} while(1);
		}
abort:
		if (facsdata)
			free(facsdata);
		if (facsname)
			free(facsname);
		if (clusterid)
			free(clusterid);
		if (clusterhistory)
			free(clusterhistory);
		MPI_Finalize();
		return(0);
	} // f

		MPI_Finalize();
		return(1);

} /* main */
/* ------------------------------------------------------------------------------------ */
