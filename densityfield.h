/*  File defining the densfield class and its routines, including debiased
    mode subtraction
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

#ifndef DENSFIELD
#define DENSFIELD

#include <cmath>
#include <complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "countlines.h"
#include "powerspec.h"

// set the binning (hard wired)
const int binnum=29;
const double bins[binnum+1]={0,0.006,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.2};

// define density field class
class densfield {
	int size;										// number of modes
	bool init;									// flag for whether effective k-value and number of modes for bins have been initialised
	double SN;									// shot noise
	int tempnum;								// number of templates
	double RP;									// parameter as defined in Eq. (15) of Kalus, Percival, Bacon & Samushia (arXiv:1607.02417)
	double SP;									// Eq. (19)
	double complex *F;					// array of field values
	double complex *Ff;					// array of template values
	// k-vector entries
	double *k1;
	double *k2;
	double *k3;
	int *binid;									// id of bins into which modes falls
	double *keff;								// effective k-values of the bins
	int *Nk;										// number of modes in bins
	gsl_matrix_complex * R;			// multi-template extension of RP
	gsl_vector_complex * S;			// multi-template extension of SP
	gsl_matrix_complex * Rinv;	// inverse of R

	public:

		// constructor of an uninitialised density field container
		densfield(int size): size(size){
				k1 = (double*)malloc(sizeof(double)*size);
				k2 = (double*)malloc(sizeof(double)*size);
				k3 = (double*)malloc(sizeof(double)*size);
				binid=(int*)malloc(sizeof(int)*size);
				keff=(double*)malloc(sizeof(double)*size);
				Nk=(int*)malloc(sizeof(int)*size);
				Ff = (double complex*)malloc(sizeof(double complex)*size);
				F = (double complex*)malloc(sizeof(double complex)*size);
				for (int ik=0; ik<size; ik++){
						binid[ik]=0;
						k1[ik]=0;
						k2[ik]=0;
						k3[ik]=0;
						F[ik]=0;
						Ff[ik]=0;
				}
				for (int ibin=0; ibin<binnum; ibin++) {
						keff[ibin]=0;
						Nk[ibin]=0;
				}
				SN=0;
				init=false;
				tempnum=1;
		}

	// constructor for a density field reading in a single template field and a density field
	densfield(const char *tempfile, const char *densfile){
		FILE *file;
		if ((file=fopen(tempfile,"r"))==NULL){
			std::cerr<<"error: "<<tempfile<<" cannot be opened"<<std::endl;
			exit(-1);
		}
		const int linecount=countlines(file);
		size=linecount;
		double Ffr, Ffi;																							// dummy variables for real and imaginary parts of template
		k1 = (double*)malloc(sizeof(double)*linecount);
		k2 = (double*)malloc(sizeof(double)*linecount);
		k3 = (double*)malloc(sizeof(double)*linecount);
		binid=(int*)malloc(sizeof(int)*linecount);
		keff=(double*)malloc(sizeof(double)*linecount);
		Nk=(int*)malloc(sizeof(int)*linecount);
		Ff = (double complex*)malloc(sizeof(double complex)*linecount);
		F = (double complex*)malloc(sizeof(double complex)*linecount);
		for (int ik=0; ik<linecount; ik++){
			binid[ik]=0;
		}
		for (int ibin=0; ibin<binnum; ibin++) {
			keff[ibin]=0;
			Nk[ibin]=0;
		}
		const int linelength=100;
		char readline[linelength];
		int ik=0;
		while (fgets(readline, linelength, file)){
			sscanf(readline, "%lf%lf%lf%lf%lf", &k1[ik], &k2[ik], &k3[ik], &Ffr, &Ffi);
			Ff[ik]=Ffr+I*Ffi;
			++ik;
		}
		fclose(file);
		if ((file=fopen(densfile,"r"))==NULL){
			std::cerr<<"error: "<<densfile<<" cannot be opened"<<std::endl;
			exit(-1);
		}
		double Fr, Fi;																							// dummy variables for real and imaginary parts of field values
		ik=0;
		while (fgets(readline, linelength, file) && ik<linecount){
			double k1b=k1[ik];
			double k2b=k2[ik];
			double k3b=k3[ik];
			double k=sqrt(k1b*k1b+k2b*k2b+k3b*k3b);
			sscanf(readline, "%le%le%le%le%le", &k1[ik], &k2[ik], &k3[ik], &Fr, &Fi);
			// check whether k-values in density and template files agree
			if ((fabs(k1b-k1[ik])>0.00001)||(fabs(k2b-k2[ik])>0.00001)||(fabs(k3b-k3[ik])>0.00001)) {
				std::cout<<ik<<'\t'<<k1b<<"="<<k1[ik]<<'\t'<<k2b<<"="<<k2[ik]<<'\t'<<k3b<<"="<<k3[ik]<<std::endl;
				std::cerr<<"k-values of field and template do not agree\n"<<readline<<std::endl;
				exit(-2);
			}
			F[ik]=Fr+I*Fi;
			for(int bid=0; bid<binnum; bid++){
				if (bins[bid]<k && bins[bid+1]>k) {
					keff[bid]+=k;																				// sum up k values to obtain effective k of bin
					++Nk[bid];																					// count modes in bin
					binid[ik]=bid;																			// assign bin to mode
					break;
				}
			}
			++ik;
		}
		for (int ibin=0; ibin<binnum; ibin++) {
				keff[ibin]/=Nk[ibin];																	// devide sum of k-values by number of bins
		}
		fclose(file);
		SN=0;																											// set shot noise to zero by default
		init=true;
		tempnum=1;																								// this constructor only works for a single template
	}

	// constructor for a density field with an arbitrary number of templates tempnum
	densfield(const int tempnum, const char *tempfilelist[]):tempnum(tempnum){
		FILE *file;
		double Ffr, Ffi;
		if ((file=fopen(tempfilelist[1],"r"))==NULL){
			std::cerr<<"error: "<<tempfilelist[1]<<" cannot be opened"<<std::endl;
			exit(-1);
		}
		const int linecount=countlines(file);
		size=linecount;
		k1 = (double*)malloc(sizeof(double)*linecount);
		k2 = (double*)malloc(sizeof(double)*linecount);
		k3 = (double*)malloc(sizeof(double)*linecount);
		binid=(int*)malloc(sizeof(int)*linecount);
		keff=(double*)malloc(sizeof(double)*linecount);
		Nk=(int*)malloc(sizeof(int)*linecount);
		Ff = (double complex*)malloc(sizeof(double complex)*linecount*tempnum);				// array comprises all templates
		F = (double complex*)malloc(sizeof(double complex)*linecount);
		double Fr, Fi;
		const int linelength=100;
		char readline[linelength];
		int ik=0;
		while (fgets(readline, linelength, file) && ik<linecount){
			sscanf(readline, "%lf%lf%lf%lf%lf", &k1[ik], &k2[ik], &k3[ik], &Fr, &Fi);
			F[ik]=Fr+I*Fi;
			++ik;
		}
		fclose(file);
		// iterate over list of template files
		for(int temp=0; temp<tempnum; temp++){
			if ((file=fopen(tempfilelist[temp+2],"r"))==NULL){
				std::cerr<<"error: "<<tempfilelist[temp+2]<<" cannot be opened"<<std::endl;
				exit(-1);
			}
			ik=0;
			while (fgets(readline, linelength, file)){
				double k1b=k1[ik];
				double k2b=k2[ik];
				double k3b=k3[ik];
				double k=sqrt(k1b*k1b+k2b*k2b+k3b*k3b);
				sscanf(readline, "%lf%lf%lf%lf%lf", &k1[ik], &k2[ik], &k3[ik], &Ffr, &Ffi);
				Ff[ik+temp*size]=Ffr+I*Ffi;
				if ((fabs(k1b-k1[ik])>0.00001)||(fabs(k2b-k2[ik])>0.00001)||(fabs(k3b-k3[ik])>0.00001)) {
					std::cout<<ik<<'\t'<<k1b<<"="<<k1[ik]<<'\t'<<k2b<<"="<<k2[ik]<<'\t'<<k3b<<"="<<k3[ik]<<std::endl;
					std::cerr<<"k-values of field and template do not agree\n"<<readline<<std::endl;
					exit(-2);
				}
				++ik;
			}
			fclose(file);
		}
		SN=0;
		initialise();
	}

	// copy constructor
	densfield(const densfield& other){
		size=other.size;
		k1 = (double*)malloc(sizeof(double)*size);
		k2 = (double*)malloc(sizeof(double)*size);
		k3 = (double*)malloc(sizeof(double)*size);
		Ff = (double complex*)malloc(sizeof(double complex)*size);
		F = (double complex*)malloc(sizeof(double complex)*size);
		for (int ik=0; ik<size; ik++) {
			k1[ik]=other.k1[ik];
			k2[ik]=other.k2[ik];
			k3[ik]=other.k3[ik];
			Ff[ik]=other.Ff[ik];
			F[ik]=other.F[ik];
		}
		init=other.init;
		SN=other.SN;
		tempnum=other.tempnum;
	}

	// destructor
	~densfield(){
		if (k1) {
			free(k1);
		}
		if (k2) {
			free(k2);
		}
		if (k3) {
			free(k3);
		}
		if (Ff) {
			free(Ff);
		}
		if (F) {
			free(F);
		}
	}

	// overload assignment operator almost like copy constructor
	densfield& operator=(const densfield& other){
		if (this != &other) {
			size=other.size;
			k1 = (double*)malloc(sizeof(double)*size);
			k2 = (double*)malloc(sizeof(double)*size);
			k3 = (double*)malloc(sizeof(double)*size);
			Ff = (double complex*)malloc(sizeof(double complex)*size);
			F = (double complex*)malloc(sizeof(double complex)*size);
			for (int ik=0; ik<size; ik++) {
				k1[ik]=other.k1[ik];
				k2[ik]=other.k2[ik];
				k3[ik]=other.k3[ik];
				Ff[ik]=other.Ff[ik];
				F[ik]=other.F[ik];
			}
		}
		init=other.init;
		SN=other.SN;
		tempnum=other.tempnum;
		return *this;
	}

	// procedure to set single entries
	void set_entry(int ik, double k1s, double k2s, double k3s, double Fr, double Fi, double Ffr, double Ffi){
		k1[ik]=k1s;
		k2[ik]=k2s;
		k3[ik]=k3s;
		Ff[ik]=Ffr+I*Ffi;
		F[ik]=Fr+I*Fi;
		init=false;
	}

	// various get functions
	double get_k1(int ik){
		return k1[ik];
	}

	double get_k2(int ik){
		return k2[ik];
	}

	double get_k3(int ik){
		return k3[ik];
	}

	double get_k(int ik){
		if (ik<size){
			return sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
		}else{
			std::cout<<ik<<" is out of range"<<std::endl;
			exit(-5);
		}
	}

	double get_Fr(int ik){
		if (ik<size){
			return creal(F[ik]);
		}else{
			std::cout<<ik<<" is out of range"<<std::endl;
			exit(-5);
		}
	}

	double get_Fi(int ik){
		return cimag(F[ik]);
	}

	double get_Ffr(int ik){
		return creal(Ff[ik]);
	}

	double get_Ffi(int ik){
		return cimag(Ff[ik]);
	}

	// set the shot noise value, computed using external code
	void setSN(double SNnew){
			SN=SNnew;
	}

	// compute effective k-values and count modes in bins, set RP and SP to zero
	void initialise(){
		for (int ibin=0; ibin<binnum; ibin++) {
				keff[ibin]=0;
				Nk[ibin]=0;
		}
		for(int ik=0; ik<size; ik++){
				double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
				for(int bid=0; bid<binnum; bid++){
					if (bins[bid]<k && bins[bid+1]>k) {
								keff[bid]+=k;
								if (k!=0.) ++Nk[bid];
								binid[ik]=bid;
								break;
					}
				}
		}
		for (int ibin=0; ibin<binnum; ibin++) {
				keff[ibin]/=Nk[ibin];
		}
		RP=0;
		SP=0;
		init=true;
	}

	// compute RP and SP as defined in Kalus, Percival, Bacon & Samushia (2016)
	void calc_RP_SP(PowerSpec Pmod){
		if(!init) initialise();
		RP=0;
		SP=0;
		for (int ik=0; ik<size; ik++){
				double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
				if (Pmod.getP(binid[ik])>0){
					RP+=2*cabs(Ff[ik])*cabs(Ff[ik])/Pmod.getP(binid[ik]);
					SP+=2*creal(conj(F[ik])*Ff[ik])/Pmod.getP(binid[ik]);
				}
		}
		std::cout<<"RP="<<RP<<"\nSP="<<SP<<std::endl;
		std::cout<<"epsilonBF="<<SP/RP<<"+-"<<sqrt(1./RP)<<std::endl<<std::endl;
	}

	// compute R and S for more than one template
	void calc_R_S(PowerSpec Pmodel){
		if(!init) initialise();
		R = gsl_matrix_complex_alloc (tempnum, tempnum);
		S = gsl_vector_complex_alloc (tempnum);
		gsl_matrix_complex_set_zero (R);
		gsl_vector_complex_set_zero (S);
		for (int ik=0; ik<size; ik++) {
			double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
			for (int A=0; A<tempnum; A++){
				double complex Sadd=Ff[ik+size*A]*conj(F[ik])/Pmodel.getP(binid[ik]);
				gsl_vector_complex_set (S, A, gsl_complex_add(gsl_vector_complex_get (S, A), gsl_complex_rect(creal(Sadd),cimag(Sadd))));
				for (int B=0; B<tempnum; B++){
					double complex Radd=Ff[ik+size*A]*conj(Ff[ik+size*B])/Pmodel.getP(binid[ik]);
					gsl_matrix_complex_set (R, A, B, gsl_complex_add(gsl_matrix_complex_get (R, A, B), gsl_complex_rect(creal(Radd),cimag(Radd))));
				}
			}
		}
		// also invert R to get Rinv to find best fitting epsilon
		int s;
		Rinv = gsl_matrix_complex_alloc (tempnum, tempnum);
		gsl_permutation * p = gsl_permutation_alloc (tempnum);
		gsl_linalg_complex_LU_decomp (R, p, &s);
		gsl_linalg_complex_LU_invert (R, p, Rinv);
		gsl_permutation_free (p);
	}

	// get best fitting epsilon for single template as the ratio of SP and RP
	double get_epsilon(PowerSpec Pmod){
		calc_RP_SP(Pmod);
		return SP/RP;
	}

	// get the number of modes in a k-bin
	int getNk(int bin){
			if(!init) initialise();
			if (bin<binnum) {
					return Nk[bin];
			}
			return 0;
	}

	// get the size of the density field array
	int getsize(){
			return size;
	}

	// print a list of the modes and their associated field and template values (only implemented for single template)
	void print(){
			for (int ik=0; ik<size; ik++) {
					std::cout<<k1[ik]<<'\t'<<k2[ik]<<'\t'<<k3[ik]<<'\t'<<creal(F[ik])<<((cimag(F[ik])>0)?'+':'-')<<abs(cimag(F[ik]))<<"i\t"<<creal(Ff[ik])<<((cimag(Ff[ik])>0)?'+':'-')<<abs(cimag(Ff[ik]))<<'i'<<std::endl;
			}
	}

	// sample the model power spectrum by a weighted average of model power spectrum values for k modes inside each bin
	PowerSpec Pmodsample(PowerSpec Pmod){
			if(!init) initialise();
			double Phat[binnum+1];
			for (int bin=0; bin<binnum+1; bin++) {
					Phat[bin]=0;
			}
			for (int ik=0; ik<size; ik++) {
					double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
					Phat[binid[ik]]+=Pmod.Pm(k)/Nk[binid[ik]];
			}
			return PowerSpec(binnum, keff, Phat);
	}

	// estimate the power spectrum as the average absolute density field value squared in a k-bin without accounting for the systematic templates
	PowerSpec PFKP(){
		if(!init) initialise();
		double Phat[binnum+1]={0};
		for (int bin=0; bin<binnum+1; bin++) {
			Phat[bin]=-SN;
		}
		for (int ik=0; ik<size; ik++) {
			Phat[binid[ik]]+=cabs(F[ik])*cabs(F[ik])/Nk[binid[ik]];
		}
		return PowerSpec(binnum, keff, Phat);
	}

	// debiased mode subtraction
	PowerSpec PepsilonBF(PowerSpec Pmod){
			if(!init) initialise();
			double Phat[binnum]={0};
			for (int bin=0; bin<binnum; bin++) {
					Phat[bin]=0;
			}
			// implementation for single template case
			if (tempnum==1){
				calc_RP_SP(Pmod);
				double epsilonbf=SP/RP;
				for (int ik=0; ik<size; ik++) {
					double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
				 	Phat[binid[ik]]+=pow(cabs(F[ik]-epsilonbf*Ff[ik]),2)/(1-1./RP*cabs(Ff[ik])*cabs(Ff[ik])/Pmod.Pm(k));
				}
			// implementation for an arbitrary number of templates
			}else{
				calc_R_S(Pmod);
				double epsilonBF[tempnum];											// array of best fitting epsilon for each template
				for (int A=0; A<tempnum; A++){
					epsilonBF[A]=0;
					for (int B=0; B<tempnum; B++){
						epsilonBF[A]+=GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get(Rinv, A, B), gsl_vector_complex_get(S, B)));
					}
					std::cout<<"epsilonBF"<<A<<"="<<epsilonBF[A]<<std::endl;
				}
				for (int ik=0; ik<size; ik++) {
					double complex D=F[ik];												// corrected density field is measured density field...
					for (int A=0; A<tempnum; A++){
						D-=epsilonBF[A]*Ff[ik+A*size];							// ... minus templates with best fitting amplitudes
					}
					double corrf=0;																// correction factor for debiasing
					double k=sqrt(k1[ik]*k1[ik]+k2[ik]*k2[ik]+k3[ik]*k3[ik]);
					for (int A=0; A<tempnum; A++){
						for (int B=0; B<tempnum; B++){
							double complex RinvAB=GSL_REAL(gsl_matrix_complex_get(Rinv, A, B))+I*GSL_IMAG(gsl_matrix_complex_get(Rinv, A, B));
							corrf+=creal(Ff[ik+A*size]*RinvAB*Ff[ik+B*size]);
						}
					}
					// summing up absolute corrected field values squared with bias correction
					Phat[binid[ik]]+=pow(cabs(D),2)/(1-corrf/Pmod.Pm(k));
				}
			}
			for (int bin=0; bin<binnum; bin++){
					Phat[bin]/=Nk[bin];								// devide by mode number to get average
					Phat[bin]-=SN;										// subtract shot noise
			}
			return PowerSpec(binnum, keff, Phat);
	}
};

#endif
