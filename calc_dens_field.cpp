/*  Programme to calculate the density field from a galaxy catalogue
		Call as
		.calc_dens_field randomcatalogue galaxycatalogue outputfile
    Copyright (C) 2018  Benedict Kalus based on calc_pow.c by Will Percival
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

#include "cosmoparams.h"
#include "header.h"
#include <fftw3.h>

const double Om_m = cosmoparam::OmegaM;
const double Om_v = cosmoparam::OmegaLambda;
double calc_dp(double);																	// distance to redshift

int main(int argc, char *argv[]) {
	if (argc!=4){
		std::cerr<<"Error: Too few arguments. Call as .calc_dens_field randomcatalogue galaxycatalogue outputfile"<<std::endl;
		return -4;
	}
	const double area=6851.4, REDMIN=0.43, REDMAX=0.70;		// specs for BOSS CMASS

	// *********************************************************
	// memory requirements

	long NX,NY,NZ; // size of FFT input in x, y, z
	NX=NY=NZ=128;
	const long NTOT = NX*NY*(NZ/2+1);	// size of FFT output

	double *ddg = (double*)malloc(sizeof(double)*NX*NY*NZ);												// real space galaxy number
	double *ddr = (double*)malloc(sizeof(double)*NX*NY*NZ);												// real space random number
	fftw_complex *ddgFourier = (fftw_complex*)malloc(sizeof(fftw_complex)*NTOT);	// Fourier space over-density

	// *********************************************************

	float XMIN = -1800;
	float XMAX = 1700;
	float YMIN = -1800;
	float YMAX = 1700;
	float ZMIN = -1800;
	float ZMAX = 1700;

	float dx  = (XMAX-XMIN)/(float)NX;
	float dy  = (YMAX-YMIN)/(float)NY;
	float dz  = (ZMAX-ZMIN)/(float)NZ;

	// set up fast Fourier transform

	fftw_plan dp_r2c;
	dp_r2c = fftw_plan_dft_r2c_3d(NX,NY,NZ,ddg,ddgFourier,FFTW_ESTIMATE);

	// allocate memory for galaxies and randoms
	const long MAX_GAL = 100000000;
	struct use_gal *gal;
	if(!(gal = (struct use_gal*)malloc(MAX_GAL*sizeof(struct use_gal))-1)) {
		std::cerr<<"memory allocation problem for galaxies"<<std::endl;
		exit(-42);
	}

	const long MAX_RAN = 100000000;
	struct use_gal *ran;
	if(!(ran = (struct use_gal*)malloc(MAX_RAN*sizeof(struct use_gal))-1)) {
		std::cerr<<"memory allocation problem for randoms"<<std::endl;
		exit(-43);
	}

	FILE *fout;

	// *********************************************************
	// read in random data

	long nran = MAX_RAN;
	read_gal_file(ran,argv[1],&nran,stderr,REDMIN,REDMAX,12);

	// *********************************************************
	// read in galaxy data

	long ngal = MAX_GAL;
	read_gal_file(gal,argv[2],&ngal,stderr,REDMIN,REDMAX,12);

	// *********************************************************
	// apply nbar and weights to galaxies & randoms
	// assumes read in weights are systematic only (no FKP)


	const double areaRad=area/(180.*180.)*pi*pi;
	const int numrbins=1000;													// number of redshift bins in nbar calculation
	const double rmin=calc_dp(REDMIN);								// mininum distance in survey
	const double rmax=calc_dp(REDMAX);								// mininum distance in survey
	const double dr=(rmax-rmin)/numrbins;							// width of bins
	double numgalrbin[numrbins]={0};									// number of galaxies in each bin
	double numranrbin[numrbins]={0};									// number of randoms in each bin
	double Volrbin[numrbins];													// volume of each bin
	for (int ir=0; ir<numrbins; ir++) {
		double rlow=rmin+ir*dr;
		double rup=rmin+(ir+1)*dr;
		Volrbin[ir]=areaRad/3.*(rup*rup*rup-rlow*rlow*rlow);
	}
  for(long ig=1;ig<=ngal;ig++) {
		double wfkp = 1.0/(1.0+Pfkp*gal[ig].nbar);
		gal[ig].wght *= wfkp;
		numgalrbin[(int)((gal[ig].dist()-rmin)/dr)]+=gal[ig].wght;
	}

  for(long ir=1;ir<=nran;ir++) {
		double wfkp = 1.0/(1.0+Pfkp*ran[ir].nbar);
		ran[ir].wght *= wfkp;
		numranrbin[(int)((ran[ir].dist()-rmin)/dr)]+=ran[ir].wght;
  }

	// *********************************************************
	// some integrals used in the code
  double gal_nbw=0.0, gal_nbwsq=0.0, gal_nbsqwsq=0.0;
  double ran_nbw=0.0, ran_nbwsq=0.0, ran_nbsqwsq=0.0;
  for(long ig=1;ig<=ngal;ig++) {
		gal_nbw      +=              gal[ig].wght;
		gal_nbwsq    +=              gal[ig].wght*gal[ig].wght;
  }

  for(long ir=1;ir<=nran;ir++) {
		ran_nbw     +=              ran[ir].wght;
		ran_nbwsq   +=              ran[ir].wght*ran[ir].wght;
  }

	for (int ib=0; ib<numrbins; ib++) {
		gal_nbsqwsq += (rmin+(ib+0.5)*dr)*(rmin+(ib+0.5)*dr)*numgalrbin[ib]*numgalrbin[ib]/(Volrbin[ib]*Volrbin[ib]);
		ran_nbsqwsq += (rmin+(ib+0.5)*dr)*(rmin+(ib+0.5)*dr)*numranrbin[ib]*numranrbin[ib]/(Volrbin[ib]*Volrbin[ib]);
	}
	gal_nbsqwsq*=areaRad*dr;
	ran_nbsqwsq*=areaRad*dr;
	double SN = gal_nbwsq/gal_nbsqwsq;					// shot noise
  double alpha = gal_nbw / ran_nbw;

	SN          *= 1.+alpha;
	std::cout<<"shot noise: "<<SN<<std::endl;
	ran_nbwsq   *= alpha*alpha;
	ran_nbsqwsq *= alpha;

	// calculate Nyquist frequency
	double ny_x = (pi/dx);
	double ny_y = (pi/dy);
	double ny_z = (pi/dz);
	double min_nyquist=ny_x;
	if(ny_y<min_nyquist) min_nyquist=ny_y;
	if(ny_z<min_nyquist) min_nyquist=ny_z;
	std::cout<<"Nyquist: "<<min_nyquist<<std::endl;

	// big file for density
	for(long ind=0;ind<NX*NY*NZ;ind++){
		ddg[ind]=0.0;
		ddr[ind]=0.0;
	}

	// add galaxies to grid
	for(long ig=1;ig<=ngal;ig++) {

		long ix    = (long)( (double)(gal[ig].x-XMIN)/dx );
		long iy    = (long)( (double)(gal[ig].y-YMIN)/dy );
		long iz    = (long)( (double)(gal[ig].z-ZMIN)/dz );

		long ipos = iz+NZ*(iy+NY*ix);
		ddg[ipos]+=gal[ig].wght;
	}

	// subtract randoms from grid
	for(long ir=1;ir<=nran;ir++) {

		long ix    = (long)( (double)(ran[ir].x-XMIN)/dx );
		long iy    = (long)( (double)(ran[ir].y-YMIN)/dy );
		long iz    = (long)( (double)(ran[ir].z-ZMIN)/dz );

		long ipos = iz+NZ*(iy+NY*ix);
		ddr[ipos]+=ran[ir].wght*alpha;
	}
	for (long ipos=0; ipos<NX*NY*NZ; ipos++){
		ddg[ipos]-=ddr[ipos];
	}


	// Fourier transform density field

	fftw_execute(dp_r2c);


	// Work out |k| for each frequency component
	double fx, fy, fz;
	fout=fopen(argv[3],"w");
	for(int i=0;i<NX;i++) {

		// frequency in x
		if(i<=NX/2) fx = (float)i/((float)NX*dx);
		else  fx = ((float)i-(float)NX)/((float)NX*dx);

		for(int j=0;j<NY;j++) {

			// frequency in y
			if(j<=NY/2) fy = (float)j/((float)NY*dy);
			else  fy = ((float)j-(float)NY)/((float)NY*dy);

			for(int k=0;k<=NZ/2;k++) {

				// frequency in z
				fz=(float)k/((float)NZ*dz);

				// length of k vector and bin
				double fktot=2.*pi*sqrt(fx*fx+fy*fy+fz*fz);

				if (fktot < (0.5*min_nyquist) && fktot > 0.0){
					// set up correction for gridding - in effect we're
					// convolving the density field with a top-hat function in
					// each direction, so we're multiplying each Fourier mode by
					// a sinc function. To correct this, we therefore divide by
					// the sinc functions.
					double sinc_x = 1.0, sinc_y=1.0, sinc_z=1.0;
					double ax = pi*fx*dx;
					double ay = pi*fy*dy;
					double az = pi*fz*dz;
					if(fx!=0.0) sinc_x = sin(ax)/ax;
					if(fy!=0.0) sinc_y = sin(ay)/ay;
					if(fz!=0.0) sinc_z = sin(az)/az;
					double grid_cor = 1.0/(sinc_x*sinc_y*sinc_z*sqrt(alpha*ran_nbsqwsq));

					long ipos = k+(NZ/2+1)*(j+NY*i);
					double dkr = ddgFourier[ipos][0];
					double dki = ddgFourier[ipos][1];
					fprintf(fout,"%f\t%f\t%f\t%f\t%f\n",2.*pi*fx,2.*pi*fy,2.*pi*fz,dkr*grid_cor,dki*grid_cor);
				}
			}
		}
	}
	fclose(fout);
	return 0;
} // end main

// compute distance to redshift red
double calc_dp(double red) {
  double qsimp(double (*func)(double), double, double);
  double dpbit(double);
  return qsimp(dpbit,0.,red);
}

double dpbit(double z) {
	return 2997.92458/sqrt((1.+z)*(1.+z)*(1.+z)*Om_m+Om_v);
}
