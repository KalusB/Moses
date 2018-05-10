/*  Header file for .calc_dens_field
    Copyright (C) 2018  Benedict Kalus
    Modified version of original file by Will Percival
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
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// structures for various tables
struct use_gal {
	double x,y,z,wght,alpha,nbar;

	void print(){
		std::cout<<x<<'\t'<<y<<'\t'<<z<<'\t'<<wght<<'\t'<<alpha<<'\t'<<nbar<<std::endl;
	}

	double dist(){
		return sqrt(x*x+y*y+z*z);
	}
};

const double pi    = 3.1415926536; // Pi

// Estimate of P(k) for FKP weight
const float Pfkp=10000;


// integration functions: integration.c
double qsimp(double (*func)(double),double,double);

// SDSS table input: table_input.c
void read_gal_file(struct use_gal*,char*,long*,FILE*,double,double,int);
