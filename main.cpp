#include <iostream>
#include "densityfield.h"

// Example file to use MOSES on a data file and an arbitrary number of template files. Call the programme as ./main datafile tempfile1 tempfile2 ... tempfileN shotnoisevalue

int main(int argc, const char *argv[]){
		if (argc<3){
			std::cerr<<"Error: Too few arguments"<<std::endl;
			return -4;
		}
		PowerSpec P("./P_Patchy.txt");										// read in average power spectrum from Patchy mocks
		int tempnum=argc-3;																// determin template number by number of template files in call
		densfield F(tempnum,argv);												// read in density field and templates

	  F.setSN(atof(argv[tempnum+2]));										// set shot noise obtained from calc_F.cpp
		PowerSpec PFKP=F.PFKP();													// calculate power without error mitigation
		PowerSpec PepsilonBF=F.PepsilonBF(P);							// apply debiased mode subtraction with Patchy prior
		PowerSpec Pmod=F.Pmodsample(P);										// interpolate Patchy mocks for the same k-bins as the measured powers
		// print power spectrum results, mode projection only printed for single template
		for (int ibin=0; ibin<binnum; ibin++) {
			std::cout<<PFKP.getk(ibin)<<'\t'<<F.getNk(ibin)<<'\t'<<PFKP.getP(ibin)<<'\t'<<PepsilonBF.getP(ibin)<<'\t'<<Pmod.getP(ibin)<<std::endl;
		}
    return 0;
}
