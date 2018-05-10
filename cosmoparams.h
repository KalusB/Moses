/*  Setting cosmological parameters used in various codes
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

#ifndef COSMOPARAMS
#define COSMOPARAMS

#include <cmath>

namespace cosmoparam{
	const double H0(70);													// Hubble parameter at current epoch
	const double h(H0/100.);											// reduced Hubble parameter
	const double OmegaB(0.0226/h/h);							// present day baryon density parameter
	const double OmegaC(0.112/h/h);								// present day cold dark matter density parameter
	const double OmegaM(OmegaB+OmegaC);						// present day matter density parameter
	const double ns(0.96);												// spectral tilt of primordial power spectrum
	const double BaryonsInMatter(OmegaB/OmegaM);	// baryon to total matter ratio
	const double sig8(0.8);												// sigma8
	const double om_m(OmegaM);										// alternative name for OmegaM
	const double OmegaLambda(1-OmegaM);						// present day cosmological constant energy density parameter
	const double deltac=1.686;										// critical spherical over-density for dark matter halos to collapse at redshift 0

// matter density parameter at redshift z
	double Omega_M(double z) {
  		double aa=1./(1.+z);
  		return OmegaM/(OmegaM+OmegaLambda*aa*aa*aa);
	}

// cosmological constant energy density parameter at redshift z
	double Omega_Lambda(double z){
  		double aa=1./(1.+z);
		return OmegaLambda/(OmegaLambda+OmegaM/aa/aa/aa);;
	}
}

#endif
