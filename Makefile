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

ARCH=$(shell uname -s)
SRC=src
include config/$(ARCH).def

# the default target
default: usage

4col: prep
	$(CC) $(CFLAGS) -lm -o bin/dselect4     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -o bin/dclust4      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -o bin/cextract4    $(SRC)/cextract.c

8col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_8 -o bin/dselect8     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_8 -o bin/dclust8      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_8 -o bin/cextract8    $(SRC)/cextract.c

12col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_12 -o bin/dselect12     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_12 -o bin/dclust12      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_12 -o bin/cextract12    $(SRC)/cextract.c

16col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_16 -o bin/dselect16     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_16 -o bin/dclust16      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_16 -o bin/cextract16    $(SRC)/cextract.c

24col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_24 -o bin/dselect24     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_24 -o bin/dclust24      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_24 -o bin/cextract24    $(SRC)/cextract.c

32col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_32 -o bin/dselect32     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_32 -o bin/dclust32      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_32 -o bin/cextract32    $(SRC)/cextract.c

48col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_48 -o bin/dselect48     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_48 -o bin/dclust48      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_48 -o bin/cextract48    $(SRC)/cextract.c

52col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_52 -o bin/dselect52     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS) -DCOLUMNS_52 -o bin/dclust52      $(SRC)/dclust.c
	$(CC) $(CFLAGS)    -DCOLUMNS_52 -o bin/cextract52    $(SRC)/cextract.c

64col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_64 -DCOLUMNS_BY_32BLOCK -o bin/dselect64     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS)  -DCOLUMNS_64 -DCOLUMNS_BY_32BLOCK -o bin/dclust64      $(SRC)/dclust.c
	$(CC) $(CFLAGS)     -DCOLUMNS_64 -DCOLUMNS_BY_32BLOCK -o bin/cextract64    $(SRC)/cextract.c

128col: prep
	$(CC) $(CFLAGS) -lm -DCOLUMNS_BY_32BLOCK -o bin/dselect128     $(SRC)/dselect.c
	$(MPICC) $(CFLAGS)  -DCOLUMNS_BY_32BLOCK -o bin/dclust128      $(SRC)/dclust.c
	$(CC) $(CFLAGS)     -DCOLUMNS_BY_32BLOCK -o bin/cextract128    $(SRC)/cextract.c



all: 4col 8col 12col 16col 24col 32col 48col 52col 64col 128col


prep:
	@type $(MPICC) >/dev/null 2>&1 || (echo "ERROR: $(MPICC) is not installed"; exit 1)
	-mkdir -p bin

clean:
	-rm -f *.o  bin/dselect*  bin/dclust*  bin/cextract*

usage:
	@echo "Possible usage:"
	@echo "make clean"
	@echo "make 4col"
	@echo "make 8col"
	@echo "make 12col"
	@echo "make 16col"
	@echo "make 24col"
	@echo "make 32col"
	@echo "make 48col"
	@echo "make 52col"
	@echo "make 64col"
	@echo "make 128col"
	@echo "make all"

