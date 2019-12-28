#!/bin/sh

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


if [ $# -eq 0 ]; then
	CPUs=8
else
	CPUs=$1
fi

OK2RUN=$(which mpirun)
if [ ! -x "$OK2RUN" ]; then 
  echo $OK2RUN
  echo "Error: mpirun is not installed on your system"
  echo
  exit 1
fi

if [ ! -f "./test/shapes.csv" ]; then
  echo "Error: test input data not found"
  echo
  exit 1
fi

if [ ! -f "./test/unit_test1.expected" ]; then
  echo "Error: test validation data not found"
  echo
  exit 1
fi


### warning, must be an absolute path

DIR=/tmp/megaclust_test.$$

#### remove any previous test ###

if [ -e $DIR ] ; then rm -r $DIR ; fi

#### run Megaclust test ###

echo "Creating result directory: $DIR"
mkdir $DIR

echo "Preparing clustering"
cut -d, -f1-4 ./test/shapes.csv > $DIR/shapes.csv
./bin/dselect4 -i $DIR/shapes.csv -o $DIR/shapes > $DIR/shapes.dselect.log

ERR=`grep ^Error $DIR/shapes.dselect.log | wc -l`
if [ $ERR != 0 ]; then
  echo "FAILED: dselect error"
  exit 1
fi

echo "Running using $CPUs processors (this will take several minutes, you can monitor progress in file $DIR/shapes.dclust)"
./megaclust.sh -N $CPUs -C 3 -d $DIR/shapes -f 1 -l 10 -s 0.1 -k 0.5 -p 99.0 -v 1 > $DIR/shapes.megaclust.log

echo "done"
echo

CLUSTER0cnt=`grep ",0$" $DIR/shapes.clusters.sort | wc -l`
if [ $CLUSTER0cnt != 43695 ]; then
	echo "FAILED: Validation failed: unassigned=$CLUSTER0cnt but 43695 was expected."
fi

grep "^LOG: Cluster" $DIR/shapes.dclust | cut -c20- | sort -gr > $DIR/unit_test1.obtained

ERR=`diff ./test/unit_test1.expected $DIR/unit_test1.obtained | wc -l`
if [ $ERR == 0 ]; then
	echo "PASSED: Clustering Process was Successfull"
else
	echo "FAILED: Validation failed; see differences with expected results:"
	sdiff ./test/unit_test1.expected $DIR/unit_test1.obtained
fi

echo ""
echo "results are not erased, you can do it yourself with the following command"
echo "rm -r $DIR"
echo

