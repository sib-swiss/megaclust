#!/bin/bash


#	------------------------------------------------------------------------------------
#
#                                 * megaclust *
#     unbiased hierarchical density based parallel clustering of large datasets
#
#
#   Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2008-2019 Nicolas Guex
#   Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex
#
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#
#	Code:       Nicolas Guex, 2008-2019
#	Contact:    Nicolas.Guex@unil.ch
#	Repository: https://github.com/sib-swiss/megaclust
#
#
#	------------------------------------------------------------------------------------



#	This is a wrapper to simplify the use of megaclust components.
#   i.e. launching the clustering and collecting the results.


MEGACLUSTDIR=./


#------------------------------------------------------


usage()
{
cat << EOF
usage: $0 -N cpu -C columns_count -d directory -f fisrtdist -l lastdist -s stepincrement [-k PctEventsToKeepCluster | -n NumberOfEventsToKeepCluster] [-p pctAssigned] [-U] [-L] [-M] [-v level]

PURPOSE:
This wrapper script submits the recursive density clustering process to dclust.
It works with binary files created by the companion software dselect.

OPTIONS:
-N      number of CPUs to use with mpirun (defaults to 16)
-C      number of effective data columns present in input file (e.g. not taking into account the selection or row_number columns).
-d      full path (with directory) of file to process
-f      first distance to test
-l      last distance to be tested
-s      step increment for the distance test
-k      minimum percent of events needed to retain a cluster; defaults is 0.5
-n      minimum number of events needed to retain a cluster; defaults is undefined as option -k is taken
-p      pctAssigned            : Stop sampling as soon as pctAssigned events have been assigned. Defaults to 95.0 pct
-U      assign Unassigned to discovered clusters
-L      assign Leftover (see dselect) to discovered clusters
-M      report cluster Merging history
-v      verbose level (default is 0)

SEE ALSO  dselect dclust cextract

Author:  Nicolas Guex 2009-2019
Contact: Nicolas.Guex@unil.ch
EOF
}

#------------------------------------------------------

CPUs=16
DIR=
FIRSTDIST=
LASTDIST=
K=0.5
N=
VERBOSE=0
COLUMNS=4
PCTASSIGNED=95.0
ASSIGN_LEFTOVER=
ASSIGN_UNASSIGNED=
REPORT_MERGING_HIST=

while getopts "N:C:d:f:l:s:k:n:p:v:ULM" OPTION
do
     case $OPTION in
         N)
             CPUs=$OPTARG
             ;;
         C)
             COLCNT=$OPTARG
             ;;
         d)
             DIR=$OPTARG
             ;;
         f)
             FIRSTDIST=$OPTARG
             ;;
         l)
             LASTDIST=$OPTARG
             ;;
         s)
             STEP=$OPTARG
             ;;
         k)
             K=$OPTARG
             ;;
         n)
             N=$OPTARG
             ;;
         p)
             PCTASSIGNED=$OPTARG
             ;;
         v)
             VERBOSE=$OPTARG
             ;;
         U)
             ASSIGN_UNASSIGNED='-U'
             ;;
         L)
             ASSIGN_LEFTOVER='-L'
             ;;
         M)
             REPORT_MERGING_HIST='-M'
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $DIR ]] || [[ -z $FIRSTDIST ]] || [[ -z $LASTDIST ]] || [[ -z $STEP ]]
then
     usage
     exit 1
fi

if [ $COLCNT -gt 128 ]; then
 usage
 exit 1
fi

if [ $COLCNT -gt 4 ]; then
 COLUMNS=8
fi
if [ $COLCNT -gt 8 ]; then
 COLUMNS=12
fi
if [ $COLCNT -gt 12 ]; then
 COLUMNS=16
fi
if [ $COLCNT -gt 16 ]; then
 COLUMNS=24
fi
if [ $COLCNT -gt 24 ]; then
 COLUMNS=32
fi
if [ $COLCNT -gt 32 ]; then
 COLUMNS=48
fi
if [ $COLCNT -gt 48 ]; then
 COLUMNS=52
fi
if [ $COLCNT -gt 52 ]; then
 COLUMNS=64
fi
if [ $COLCNT -gt 64 ]; then
 COLUMNS=128
fi

#------------------------------------------------------
MPIRUN="mpirun -n $CPUs"
BINDIR=$MEGACLUSTDIR/bin

DCLUST="$MPIRUN $BINDIR/dclust$COLUMNS"
CEXTRACT=$BINDIR/cextract$COLUMNS


if [ ! -f $DIR.selected ]; then
 echo "dclust input file $DIR.selected is missing."
 echo "use dselect to create a suitable file."
 exit 1
fi


FINALCLUSTERCNT=0
# FIRSTDIST=`printf "%0.3f" $FIRSTDIST`
# LASTDIST=`printf "%0.3f" $LASTDIST`
# STEP=`printf "%0.3f" $STEP`

if [ -z $N ]; then
 echo "megaclust.sh: Processing $DIR with distance $FIRSTDIST and $K PctEventsToKeepCluster"
 COUNTOPTION="-k $K"
else
 echo "megaclust.sh: Processing $DIR with distance $FIRSTDIST and $N EventsToKeepCluster"
 COUNTOPTION="-n $N"
fi


echo "megaclust.sh:     $DCLUST -i $DIR.selected -f $FIRSTDIST -l $LASTDIST -s $STEP -p $PCTASSIGNED $COUNTOPTION -v $VERBOSE $ASSIGN_UNASSIGNED $ASSIGN_LEFTOVER $REPORT_MERGING_HIST -g > $DIR.dclust"
                        $DCLUST -i $DIR.selected -f $FIRSTDIST -l $LASTDIST -s $STEP -p $PCTASSIGNED $COUNTOPTION -v $VERBOSE $ASSIGN_UNASSIGNED $ASSIGN_LEFTOVER $REPORT_MERGING_HIST -g > $DIR.dclust



CMD="$CEXTRACT -f $DIR -a -o $DIR.clusters"
if [ $VERBOSE -gt 0 ]; then echo "megaclust.sh: $CMD"; fi
$CMD > $DIR.cextract.log

CMD="$CEXTRACT -U $DIR.selected.unassigned -o $DIR.unassigned"
if [ $VERBOSE -gt 0 ]; then echo "megaclust.sh: $CMD"; fi
$CMD >> $DIR.cextract.log

CMD="cat $DIR.clusters $DIR.unassigned"
if [ $VERBOSE -gt 0 ]; then echo "$CMD | sed -e s/,.*,/,/g | sort -g > $DIR.clusters.sort"; fi
$CMD | sed -e s/,.*,/,/g | sort -g > $DIR.clusters.sort

if [ -f $DIR.leftover.clusters ]; then
 mv $DIR.clusters.sort $DIR.clusters.sort.noleftover
 cat $DIR.clusters.sort.noleftover $DIR.leftover.clusters  | sort -g > $DIR.clusters.sort
fi

