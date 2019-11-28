/*  Example file to use MOSES on a data file and an arbitrary number of
    template files. Call the programme as
    ./main datafile tempfile1 tempfile2 ... tempfileN shotnoisevalue
    Copyright (C) 2018  Benedict Kalus

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include "densityfield.h"


int main(int argc, const char *argv[]){
		if (argc<3){
			std::cerr<<"Error: Too few arguments: datafile tempfile1 tempfile2 ... tempfileN shotnoisevalue"<<std::endl;
			return -4;
		}
		PowerSpec P("./P_Patchy.txt");										// read in average power spectrum from Patchy mocks
		int tempnum=argc-3;																// determin template number by number of template files in call
		densfield F(tempnum,argv);												// read in density field and templates

	  F.setSN(atof(argv[tempnum+2]));										// set shot noise obtained from calc_F.cpp
		PowerSpec PFKP=F.PFKP();													// calculate power without error mitigation
		PowerSpec PepsilonBF=F.PepsilonBF(P);							// apply debiased mode subtraction with Patchy prior
		PowerSpec Pmod=F.Pmodsample(P);										// interpolate Patchy mocks for the same k-bins as the measured powers
		// print power spectrum results
		for (int ibin=0; ibin<binnum; ibin++) {
			std::cout<<PFKP.getk(ibin)<<'\t'<<F.getNk(ibin)<<'\t'<<PFKP.getP(ibin)<<'\t'<<PepsilonBF.getP(ibin)<<'\t'<<Pmod.getP(ibin)<<std::endl;
		}
		std::cout<<"Cleaned density field:"<<std::endl;
		F.get_cleaned(PepsilonBF).print();
    return 0;
}
