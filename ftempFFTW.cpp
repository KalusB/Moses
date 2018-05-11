
/*  Programme to Fourier transform a template in configuration space
    to be used for Mode Subtraction
    Call as
    ./ftempFFTW infile outfile
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

#include "cosmoparams.h"
#include "header.h"
#include <fftw3.h>

const double Om_m = cosmoparam::OmegaM;
const double Om_v = cosmoparam::OmegaLambda;
double calc_dp(double);

int main(int argc, char *argv[]) {
	if (argc!=3){
		std::cerr<<"Error: ftempFFTW infile outfile"<<std::endl;
		return -1;
	}

	// *********************************************************
	// memory requirements

	long NX,NY,NZ; // size of FFT input in x, y, z
	NX=NY=NZ=128;
	const long NTOT = NX*NY*(NZ/2+1);	// size of FFT output

	double *ddg = (double*)malloc(sizeof(double)*NX*NY*NZ);
	fftw_complex *ddgFourier = (fftw_complex*)malloc(sizeof(fftw_complex)*NTOT);

	fftw_plan dp_r2c;
	dp_r2c = fftw_plan_dft_r2c_3d(NX,NY,NZ,ddg,ddgFourier,FFTW_ESTIMATE);

	float XMIN = -1800;
	float XMAX = 1700;
	float YMIN = -1800;
	float YMAX = 1700;
	float ZMIN = -1800;
	float ZMAX = 1700;

	float dx  = (XMAX-XMIN)/(float)NX;
	float dy  = (YMAX-YMIN)/(float)NY;
	float dz  = (ZMAX-ZMIN)/(float)NZ;

	// calculate Nyquist frequency
	double ny_x = (pi/dx);
	double ny_y = (pi/dy);
	double ny_z = (pi/dz);
	double min_nyquist=ny_x;
	if(ny_y<min_nyquist) min_nyquist=ny_y;
	if(ny_z<min_nyquist) min_nyquist=ny_z;

	FILE* infile;
	if ((infile=fopen(argv[1],"r"))==NULL){
		std::cerr<<"error: "<<argv[1]<<" cannot be opened"<<std::endl;
		exit(-1);
	}
	// big file for density
	for(long ind=0;ind<NX*NY*NZ;ind++){
		double buf;
		char buff[100];
		fgets(buff,100,infile);
		sscanf(buff,"%le",&ddg[ind]);
	}

	puts("big file for density");

	// Fourier transform density field
	fftw_execute(dp_r2c);

	// Work out |k| for each frequency component
	double fx, fy, fz;
	FILE *fout;
	fout=fopen(argv[2],"w");
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
					double grid_cor = 1.0/(sinc_x*sinc_y*sinc_z);

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

double calc_dp(double red) {
  double qsimp(double (*func)(double), double, double);
  double dpbit(double);
  return qsimp(dpbit,0.,red);
}

double dpbit(double z) {
	return 2997.92458/sqrt((1.+z)*(1.+z)*(1.+z)*Om_m+Om_v);
}
