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

	cextract.c: Summarizes results from the clustering process and extract clusters from
				binary files.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <float.h>
#include "dclust.h"

/* ------------------------------------------------------------------------------------ */

typedef	struct	FACSNAME_struct	FACSNAME;
struct	FACSNAME_struct
{
	char name[kMaxCellName];	
};


typedef	struct	FACSDATA_struct	FACSDATA;
struct	FACSDATA_struct
{
	float	data[kMaxInputCol];
};


static int firsttime = 1;
static unsigned int extractedcnt = 0;
/* ------------------------------------------------------------------------------------ */

static void PrintSummaryTable(FILE *f,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int *summary,unsigned int maxclusterid)
{
	unsigned short cn;
	unsigned int clusterid;
	
		for (cn = 0; cn<cellnamecnt;cn++)
		{
			fprintf(f,"%s",uniquecellnames[cn].name);
			for (clusterid = 0; clusterid<=maxclusterid;clusterid++)
			{
				fprintf(f,",%d",summary[(cn*(maxclusterid+1))+clusterid]);
			}
			fprintf(f,"\n");
		}

} /* PrintSummaryTable */

/* ------------------------------------------------------------------------------------ */

static unsigned int GetMaxClusterId(FILE *f,unsigned short cellnamecnt,unsigned int **summary)
{
	unsigned int rowcnt;
	unsigned int colcnt;
	char	hdr[kHeaderSize];
	unsigned int maxclusterid;
	int endian;
	
	fread(&hdr,sizeof(char),kHeaderSize,f);
	if (strcmp(hdr,"dclust assigned file v1.0     \n") != 0)
	{
		printf("Error: input is not a dclust assigned file\n");
		return(0);
	}	
	fread(&endian,sizeof(int),1,f);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(0);
	}
	fread(&rowcnt,sizeof(int),1,f);
	fread(&colcnt,sizeof(int),1,f);
	if (colcnt > kMaxInputCol)
	{
		printf("Error: number of requested data columns exceed maximum allowed (%d > %d)\n",colcnt,kMaxInputCol);
		return(0);
	}	

	fread(&maxclusterid,sizeof(int),1,f);


	if (!*summary)
	{
		*summary = calloc(cellnamecnt*(maxclusterid+1),sizeof(unsigned int));
		if (!*summary)
		{
			printf("Error:Cannot Allocate Memory for Summary Table.\n");
			return(0);
		}
	}

	return(maxclusterid);
	
} /* GetMaxClusterId */

/* ------------------------------------------------------------------------------------ */
static unsigned short ExtractDataFromUnassigned(FILE *f, FILE *of,unsigned int printCID)
{
	unsigned int rowcnt,rcnt;
	unsigned int colcnt;
	char	header[kMaxLineBuf];
	char	hdr[kHeaderSize];
	unsigned int mcid;
	int endian;


	fread(&hdr,sizeof(char),kHeaderSize,f);
	if (strcmp(hdr,"dclust unassigned file v1.0   \n") != 0)
	{
		printf("Error: input is not a dclust unassigned file\n");
		return(1);
	}	
	fread(&endian,sizeof(int),1,f);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(1);
	}
	fread(&rowcnt,sizeof(int),1,f);
	fread(&colcnt,sizeof(int),1,f);
	if (colcnt > kMaxInputCol)
	{
		printf("Error: number of requested data columns exceed maximum allowed (%d > %d)\n",colcnt,kMaxInputCol);
		return(1);
	}	

	fread(&mcid,sizeof(int),1,f);

	/* skip header */
	fread(&header,sizeof(char),kMaxLineBuf,f);

	if (of)
	{
		for (rcnt = 0; rcnt < rowcnt; rcnt++)
		{
			CELLNAMEIDX cellnameidx;
			unsigned short col;
			float value[kMaxInputCol]; 

			/* skip cellname */
			fread(&cellnameidx,sizeof(CELLNAMEIDX),1,f);

			fread(&value[0],sizeof(float),colcnt,f);

			fprintf(of,"%u",cellnameidx);

			for(col=0;col<colcnt;col++)
				fprintf(of,",%d",(int)value[col]);
			if (printCID)
				fprintf(of,",0\n");
			else
				fprintf(of,"\n");
			extractedcnt++;

		} 
	}	
	return(0);
	
	
} /* ExtractDataFromUnassigned */

/* ------------------------------------------------------------------------------------ */

static unsigned int loaddatabinary(FILE *f,FILE *of,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int cid)
{
	unsigned int rowcnt,rcnt;
	unsigned int colcnt;
	char	header[kMaxLineBuf];
	char	hdr[kHeaderSize];
	unsigned int mcid;
	int endian;


	fread(&hdr,sizeof(char),kHeaderSize,f);
	if (strcmp(hdr,"dclust assigned file v1.0     \n") != 0)
	{
		printf("Error: input is not a dclust assigned file\n");
		return(1);
	}	
	fwrite(&hdr,sizeof(char),kHeaderSize,of);
	
	fread(&endian,sizeof(int),1,f);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(1);
	}
	fwrite(&endian,sizeof(int),1,of);
	
	fread(&rowcnt,sizeof(int),1,f);
	fwrite(&rowcnt,sizeof(int),1,of);

	fread(&colcnt,sizeof(int),1,f);
	if (colcnt > kMaxInputCol)
	{
		printf("Error: number of requested data columns exceed maximum allowed (%d > %d)\n",colcnt,kMaxInputCol);
		return(1);
	}	
	fwrite(&colcnt,sizeof(int),1,of);

	fread(&mcid,sizeof(int),1,f);
	mcid=1; /* we extract only  one cluster */
	fwrite(&mcid,sizeof(int),1,of);


	/* skip header */
	fread(&header,sizeof(char),kMaxLineBuf,f);
	fwrite(&header,sizeof(char),kMaxLineBuf,of);

	for (rcnt = 0; rcnt < rowcnt; rcnt++)
	{
		CELLNAMEIDX cellnameidx;
		float value[kMaxInputCol]; 
		unsigned int clusterid;

		/* skip cellname */
		fread(&cellnameidx,sizeof(CELLNAMEIDX),1,f);
		fread(&value[0],sizeof(float),colcnt,f);
		fread(&clusterid,sizeof(unsigned int),1,f);
		if (clusterid == cid)
		{
			fwrite(&cellnameidx,sizeof(CELLNAMEIDX),1,of);
			fwrite(&value[0],sizeof(float),colcnt,of);
			clusterid=1;
			fwrite(&clusterid,sizeof(unsigned int),1,of);
			extractedcnt++;
		}
	} 

	/* adjust number of events assigned  in cluster */
	fseek(of,(kHeaderSize+sizeof(int)),SEEK_SET);
	fwrite(&extractedcnt,sizeof(unsigned int),1,of);

	return(0);
	
} /* loaddatabinary */

/* ------------------------------------------------------------------------------------ */
static unsigned int loaddata(FILE *f,FILE *of,FACSNAME *uniquecellnames,unsigned short cellnamecnt,unsigned int *summary,unsigned int maxclusterid,unsigned int cid,unsigned int nid,unsigned int printCID)
{
	unsigned int rowcnt,rcnt;
	unsigned int colcnt;
	char	header[kMaxLineBuf];
	char	hdr[kHeaderSize];
	unsigned int mcid;
	int endian;


	fread(&hdr,sizeof(char),kHeaderSize,f);
	if (strcmp(hdr,"dclust assigned file v1.0     \n") != 0)
	{
		printf("Error: input is not a dclust assigned file\n");
		return(1);
	}	
	fread(&endian,sizeof(int),1,f);
	if (endian != 1)
	{
		printf("Error: File not supported: Wrong platform (little/Big endian incompatibility)\n");
		return(1);
	}
	fread(&rowcnt,sizeof(int),1,f);
	fread(&colcnt,sizeof(int),1,f);
	if (colcnt > kMaxInputCol)
	{
		printf("Error: number of requested data columns exceed maximum allowed (%d > %d)\n",colcnt,kMaxInputCol);
		return(1);
	}	

	fread(&mcid,sizeof(int),1,f);
	if ((mcid > 0) && (mcid != maxclusterid))
	{
		printf("Error: Maximum ClusterID=%d differs from Maximum ClusterID=%d from .assigned file\n",mcid,maxclusterid);
		return(1);
	}	

	/* skip header */
	fread(&header,sizeof(char),kMaxLineBuf,f);

	if (of)
	{
		if (firsttime) /* make sure header is printed just once */
		{
			if (!printCID)
				header[strlen(header)-8] = 0;  /* remove the ,cluster */				
			fprintf(of,"%s\n",header);
			firsttime = 0;
		}
		for (rcnt = 0; rcnt < rowcnt; rcnt++)
		{
			CELLNAMEIDX cellnameidx;
			unsigned short col;
			float value[kMaxInputCol]; 
			unsigned int clusterid;

			/* skip cellname */
			fread(&cellnameidx,sizeof(CELLNAMEIDX),1,f);
			fread(&value[0],sizeof(float),colcnt,f);
			fread(&clusterid,sizeof(unsigned int),1,f);

			if ((cid == -1) || (clusterid == cid))
			{
				fprintf(of,"%u",cellnameidx);
				for(col=0;col<colcnt;col++)
					fprintf(of,",%d",(int)value[col]);
				if (printCID)
					fprintf(of,",%d\n",clusterid);
				else
					fprintf(of,"\n");
				extractedcnt++;
			}
		} 
	}
	else
	{
		for (rcnt = 0; rcnt < rowcnt; rcnt++)
		{
			CELLNAMEIDX cellnameidx;
			unsigned int clusterid;

			/* skip cellname */
			fread(&cellnameidx,sizeof(CELLNAMEIDX),1,f);
			fseek(f,colcnt*sizeof(float),SEEK_CUR);
			fread(&clusterid,sizeof(unsigned int),1,f);
		}
	}
	
	return(0);
	
} /* loaddata */

/* ------------------------------------------------------------------------------------ */

int main (int argc, char **argv)
{	
	int clusterid,cellnameid;
	unsigned int processAssigned = 0;
	char	*version="VERSION 1.0; 2019-12-26";
	char fn[kMaxFilename];
	char rfn[kMaxFilename];
	char sfn[kMaxFilename];
	char ofn[kMaxFilename];
	char extractcellname[kMaxCellName];
	unsigned int maxclusterid = 0;
	unsigned short cellnamecnt;
	unsigned int binary = 0;
	unsigned int ExtractFromBinaryUnassigned = 0;
	FACSNAME *uniquecellnames = NULL;
	unsigned int *summary = NULL;
	FILE *f = NULL;
	FILE *of = NULL;
	int c;
	int verbose;
	unsigned int printCID;

	/* --------- process arguments */

	fn[0] = 0;
	rfn[0] = 0;
	sfn[0] = 0;
	ofn[0] = 0;
	extractcellname[0] = 0;
	clusterid = -1;
	cellnameid = -1;
	verbose = 0;
	printCID = 1;

	opterr = 0;
	while ((c = getopt (argc, argv, "f:as:o:c:n:i:v:0bU:")) != -1)
	switch (c)
	{
      case 'f':
			strcpy(rfn,optarg);
        break;
			
	  case 'a':
			processAssigned = 1;
        break;

	  case 's':
			strcpy(sfn,optarg);
        break;

	  case 'o':
			strcpy(ofn,optarg);
        break;

	  case 'c':
			sscanf(optarg,"%d",&clusterid);
        break;

	  case 'n':
			strcpy(extractcellname,optarg);
        break;

	  case 'i':
			sscanf(optarg,"%d",&cellnameid);
        break;

	  case '0':
			printCID = 0;
        break;

	  case 'b':
			binary = 1;
        break;

	  case 'U':
			strcpy(rfn,optarg);
			ExtractFromBinaryUnassigned = 1;
        break;

	  case 'v':
			sscanf(optarg,"%d",&verbose);
        break;
	}
	
	if ( (rfn[0] == 0) || ((clusterid >= 0) && (ofn[0] == 0)) )
	{
		printf("usage\n\n");
		printf("cextract -f RootName [-a] [-s SummaryOutputFile] [-o OutputFile [-c ClusterID [-0|-b]] [-i CellNameID] [-n CellName]] [-v]  [-U UnassignedFile ]\n");
		printf("         -f RootName         : use dclust input file (RootName.selected) that contains the events selected for the clustering\n");
		printf("         -a                  : process dclust output file (RootName.selected.assigned) containing events with assigned clusterid\n");
		printf("         -s SummaryOutputFile: Write a comma separated file with the number of events for each cluster and cellname\n");
		printf("         -c ClusterID        : specifies which ClusterID to extract in OutputFile. Default is -1 for all.\n");
		printf("                               with option -0 is specified the ClusterID is not written to the last column of the .csv output file\n");
		printf("                               with option -b a dclust .assigned (with -a) or .unassigned (with -u) file will be produced instead of a .csv file\n");
		printf("                               option -b requires a unique ClusterID and is incompatible with most other options\n");
		printf("                               and must be used like this: cextract -f RootName -a -o OutputFile -c ClusterID -b\n");
		printf("                               and must be used like this: cextract -f RootName -u -o OutputFile -c 0         -b\n");
		printf("         -n CellName         : specifies which CellName to extract in OutputFile. Default is empty (all).\n");
		printf("         -i CellNameID       : specifies which CellNameID to extract in OutputFile. Default is -1 for all.\n");
		printf("                             : set verbose level (-v) to 1 to list CellNameIDs\n");
		printf("         -o OutputFile       : Name of the comma separated file with clusterids\n");
		printf("         -v level            : specifies the verbose level; default is 0.\n");
		printf("                               level 1 lists the available CellNameIDs\n");
		printf("         -U UnassignedFile   : extract events from a dclust .unassigned output file. Clusterid will be 0 for every event.\n");
		printf("                               This option is incompatible with any other option.\n\n");
		printf("VERSION\n");
		printf("\n%s\n",version);
		printf("Author:  Nicolas Guex; 2008-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
		printf("CONTACT: Nicolas.Guex@unil.ch\n");
		printf("SEE ALSO\n");
		return(1);
	}


	if ((sfn[0] == 0) && (ofn[0]==0))
	{
		printf("Warning: no SummaryOutputFile or OutputFile specified.\n");		
	}
	else
	{
		if (ExtractFromBinaryUnassigned == 0)
		if (processAssigned == 0)
			printf("Warning: option -a is not selected. Output Files will be empty.\n");
	}
	if ((printCID == 0) && (clusterid == -1))
	{
		printf("Error: attempt to output a file containing several clusters without clusterID.\n");		
		printf("       option -0 can only be used in combination with option -c\n");	
		return(1);
	}	
	if ((binary == 1) && ((clusterid == -1) || (sfn[0] != 0)))
	{
		printf("Error: attempt to output a binary file containing several clusters.\n");		
		printf("       option -b can only be used in combination with option -c followed by a unique cluster id\n");	
		printf("       option -s are incompatible with this option.\n");
		return(1);
	}	
	if (ExtractFromBinaryUnassigned && (ofn[0]==0))
	{
		printf("Error: no OutputFile specified.\n");		
		return(1);
	}

	{
	}

	/* --------- process */

	if (ofn[0])
	{
		if (binary)
			of = fopen(ofn,"wb+");
		else
			of = fopen(ofn,"w");
		if (!of)
		{
			printf("Error:Cannot Write File %s\n",ofn);
			goto bail;
		}
	}

	if (ExtractFromBinaryUnassigned)
	{
		f = fopen(rfn,"r");
		if (!f)
		{
			printf("Error:Cannot Open Input File %s\n",rfn);
			goto bail;
		}
		ExtractDataFromUnassigned(f,of,printCID);
		fclose(f);
		cellnamecnt = 0;
	}
	else
	{	
		sprintf(fn,"%s.selected",rfn);
		f = fopen(fn,"r");
		if (!f)
		{
			printf("Error:Cannot Open Input File %s\n",fn);
			goto bail;
		}
		cellnamecnt = 1;
		fclose(f);
	}

	if (cellnamecnt > 0)
	{
		/* convert cellname into cellname id */
		if (extractcellname[0])
		{
			unsigned int i;
			for (i=0; i< cellnamecnt; i++)
			{
				if (strcmp(uniquecellnames[i].name,extractcellname) == 0)
					break;
			}
			if (i < cellnamecnt)
				cellnameid = i;
			else
			{
				printf("Warning: '%s' is not a valid complete cellname\n",extractcellname);
			}
		}
		
		
		sprintf(fn,"%s.selected.assigned",rfn);
		f = fopen(fn,"r");
		if (!f)
		{
			printf("Error:Cannot Open Input File %s\n",fn);
			goto bail;
		}
		maxclusterid = GetMaxClusterId(f,cellnamecnt,&summary);

		if (processAssigned)
		{
			unsigned int err;
			printf("Processing %s  (cluster id = %d)\n",fn,clusterid);
			rewind(f);
			if (binary)
				err = loaddatabinary(f,of,uniquecellnames,cellnamecnt,clusterid);
			else
				err = loaddata(f,of,uniquecellnames,cellnamecnt,summary,maxclusterid,clusterid,cellnameid,printCID);
			fclose(f);
			if (err)
				goto bail;
		}
		else
			fclose(f);
			
		if (ofn[0])
			printf("Extracted %12d events to %s\n",extractedcnt,ofn);
		
		if (summary && sfn[0])
		{
			f = fopen(sfn,"w");
			if (f)
			{
				printf("Writing SummaryTable to %s\n",sfn);
				PrintSummaryTable(f,uniquecellnames,cellnamecnt,summary,maxclusterid);
				fclose(f);
			}
			else
				printf("Error:Cannot Write Output File %s\n",sfn);
		}		
	}
	
	if (of)
		fclose(of);
	if (summary)
		free(summary);
	free(uniquecellnames);
	return(0);
	
	/* abort */
bail:
	if (uniquecellnames)
		free(uniquecellnames);
	return(1);
	
} /* main */
/* ------------------------------------------------------------------------------------ */
