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


/*

	dclust.h: unique set of constants for the various components of megaclust.

*/


#define kMaxLineBuf 2048 
#define kMaxFilename 1024
#define kMaxCellName 32
#define kHeaderSize 32
#define kMAXEVENTS 15000000
#define kMAX_ALLOWED_INPUT_VALUE 16384

/* ------------------------------------------------------------------------------------ */

#define CELLNAMEIDX  unsigned int


#ifdef COLUMNS_BY_32BLOCK

#ifdef COLUMNS_64
#define kBLOCK32CNT 2
#else
#define kBLOCK32CNT 4
#endif
#define kMaxInputCol (32*kBLOCK32CNT)      


#else


#ifdef COLUMNS_52
#define kMaxInputCol 52      
#else
#ifdef COLUMNS_48
#define kMaxInputCol 48      
#else
#ifdef COLUMNS_32
#define kMaxInputCol 32      
#else
#ifdef COLUMNS_24
#define kMaxInputCol 24      
#else
#ifdef COLUMNS_16
#define kMaxInputCol 16      
#else
#ifdef COLUMNS_12
#define kMaxInputCol 12      
#else
#ifdef COLUMNS_8
#define kMaxInputCol 8      
#else
#define kMaxInputCol 4   
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#endif  /* !COLUMNS_BY_32BLOCK */


/* ------------------------------------------------------------------------------------ */
