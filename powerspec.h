/*  File defining the powerspec class and its routines
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

#ifndef POWERSPEC
#define POWERSPEC

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <time.h>

#include "countlines.h"

class PowerSpec{
	double *k;								// list of k-values
	double *P;								// list of power spectrum values
	unsigned ksteps;					// size of k and P arrays
	gsl_interp_accel *acc;		// used in spline interpolations
	gsl_spline *Pmodelspline;	// cubic spline based on k and P arrays
	bool splineinit;					// flag whether spline has been initialised

	public:
	// empty container constructor
  PowerSpec(){
		ksteps=0;
		k=0;
		P=0;
		splineinit=false;
	}

	// copy constructor
  PowerSpec(const PowerSpec& other){
        k=new double[other.ksteps];
        P=new double[other.ksteps];
        for (unsigned i=0; i<other.ksteps; i++) {
            k[i]=other.k[i];
            P[i]=other.P[i];
        }
				ksteps=other.ksteps;
				splineinit=false;							// would be better to copy the spline, but it didn't work naively
  }

	// constructor reading in power spectrum file Pfile
	PowerSpec(const char *Pfile){
		FILE *Pin;
		ksteps=0;
		const unsigned rllen=100;
		char readline[rllen];
		if (!(Pin=fopen(Pfile,"r"))){
			std::cerr<<"File "<<Pfile<<" cannot be opened."<<std::endl;
			exit(-10);
		}
		const unsigned linecount=countlines(Pin);
		k=new double[linecount];
		P=new double[linecount];
		while(fgets(readline,rllen,Pin)){
			if(ksteps==linecount){
				std::cerr<<"error: Power spectrum file is too long. Increase linecount in powerspec.h.\n";
				exit(-1);
			}
			sscanf(readline,"%le\t%le",&k[ksteps],&P[ksteps]);
			++ksteps;
		}
		// initialise cubic spline
		makespline();
	}

	// constructor defining PowerSpec object by k and P arrays with length
	PowerSpec(const int length, const double knew[], const double Pnew[]){
		k=new double[length];
		P=new double[length];
		for (int i=0; i<length; i++) {
			k[i]=knew[i];
			P[i]=Pnew[i];
		}
		ksteps=length;
		splineinit=false;
	}

	// destructor
	~PowerSpec(){
		if (k) delete[] k;
		if (P) delete[] P;
		if (splineinit) {
			gsl_spline_free(Pmodelspline);
			gsl_interp_accel_free(acc);
		}
		splineinit=false;
	}

	// overload assignment operator
	PowerSpec& operator=(const PowerSpec& other){
		if (this != &other) {
			k=new double[other.ksteps];
			P=new double[other.ksteps];
			for (unsigned i=0; i<other.ksteps; i++) {
				k[i]=other.k[i];
				P[i]=other.P[i];
			}
			ksteps=other.ksteps;
			splineinit=false;
		}
		return *this;
	}

	// print a list of k and P values
	void print(){
		for (unsigned i=0; i<ksteps; i++) {
			std::cout<<k[i]<<'\t'<<getP(i)<<std::endl;
		}
		std::cout<<std::endl;
	}

	// returns an interpolation value of the power spectrum at km if positive, zero otherwise
  inline double Pm(double km){
		if (!splineinit) makespline();
		double retvar=gsl_spline_eval (Pmodelspline, km, acc);
		return (retvar>0)?retvar:0;
	}

	// resample the power spectrum for a new k array knew using cubic spline interpolation
	void resample(double knew[], unsigned newkbins){
		if (ksteps<newkbins) {
			double *Pnew=new double[newkbins];
			double *kn=new double[newkbins];
			for (unsigned ik=0; ik<newkbins; ik++) {
				Pnew[ik]=Pm(knew[ik]);
				kn[ik]=knew[ik];
			}
			delete[] P;
			delete[] k;
			P=Pnew;
			k=kn;
		} else {
			for (unsigned ik=0; ik<newkbins; ik++) {
				P[ik]=Pm(knew[ik]);
				k[ik]=knew[ik];
			}
		}
		ksteps=newkbins;
		gsl_spline_free (Pmodelspline);
		makespline();
	}

	void makespline(){
		acc = gsl_interp_accel_alloc ();
		Pmodelspline = gsl_spline_alloc (gsl_interp_cspline, ksteps);
		gsl_spline_init (Pmodelspline, k, P, ksteps);
		splineinit=true;
	}

	// set new k-value
	void setk(int ik, double knew){
		k[ik]=knew;
		splineinit=false;
	}

	// set new power spectrum value
	void setP(int ik, double Pnew){
		P[ik]=Pnew;
		if (splineinit) {
			gsl_spline_free(Pmodelspline);
			gsl_interp_accel_free(acc);
		}
		acc = gsl_interp_accel_alloc ();
		Pmodelspline = gsl_spline_alloc (gsl_interp_cspline, ksteps);
		gsl_spline_init (Pmodelspline, k, P, ksteps);
		splineinit=true;
	}

	// get the size of the arrays
	unsigned getksteps(){
		return ksteps;
	}

	// overload the subscript operator returning the ith element of the power spectrum array if i is smaller than the size of the array
	double operator[](int i){
		return (i<(int)ksteps)?P[i]:0;
	}

	// get ith k value
	double getk(int i){
		return (i<(int)ksteps)?k[i]:9e99;
	}

	// get power spectrum value, just references operator[]
	inline double getP(int i){
		return P[i];
	}

};

#endif
