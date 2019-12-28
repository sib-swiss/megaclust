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
