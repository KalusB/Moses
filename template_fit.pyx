'''
    Programme to fit polynomials to the ratio of observed galaxies to
    expected galaxies as a function of a contaminant, and to generate
    templates in configuration space based on that information
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
'''

import numpy as np
cimport numpy as np
import healpy as hp
import fitsio
from scipy.integrate import quad
import GalaxyvsContaminant

cdef int NSIDE=GalaxyvsContaminant.NSIDE
cdef int NX, NY, NZ
NX=NY=NZ=128
cdef double XMIN = -1800.
cdef double XMAX = 1700.
cdef double YMIN = -1800.
cdef double YMAX = 1700.
cdef double ZMIN = -1800.
cdef double ZMAX = 1700.

cdef double dx  = float(XMAX-XMIN)/NX
cdef double dy  = float(YMAX-YMIN)/NY
cdef double dz  = float(ZMAX-ZMIN)/NZ

def loadsyst(char* systfilename):
	cdef np.ndarray fsys=np.loadtxt(systfilename)
	cdef np.ndarray mask=np.zeros(len(fsys))
	cdef int cell, neighbour
	for cell in np.arange(len(fsys)):
		if fsys[cell]==0:
			mask[cell]=0
		else:
			mask[cell]=1
			for neighbour in hp.get_interp_weights(NSIDE,cell,nest=True)[0]:
				if fsys[neighbour]==0:
					mask[cell]-=0.25
	fsys/=mask
	return fsys, mask

def polynomial2invweight(np.ndarray p1, np.ndarray fsys, np.ndarray Wtot):
	cdef np.ndarray W=np.zeros(NX*NY*NZ)
	cdef np.ndarray f1=np.poly1d(p1)(fsys)
	cdef double xmin, xmax, ymin, ymax, zmin, zmax, dpmin, dpmax, dp, x, y, z, dec, ra
	for ix in np.arange(NX):
		for iy in np.arange(NY):
			for iz in np.arange(NZ):
				ipos = iz+NZ*(iy+NY*ix)
				xmin=xmax=ymin=ymax=zmin=zmax=0
				if ix*dx+XMIN<0:
					xmin=(ix+1)*dx+XMIN
					xmax=ix*dx+XMIN
				else:
					xmin=ix*dx+XMIN
					xmax=(ix+1)*dx+XMIN
				if iy*dy+YMIN<0:
					ymin=(iy+1)*dy+YMIN
					ymax=iy*dy+YMIN
				else:
					ymin=iy*dy+YMIN
					ymax=(iy+1)*dy+YMIN
				if iz*dz+ZMIN<0:
					zmin=(iz+1)*dz+ZMIN
					zmax=iz*dz+ZMIN
				else:
					zmin=iz*dz+ZMIN
					zmax=(iz+1)*dz+ZMIN
				dpmin=np.sqrt(xmin*xmin+ymin*ymin+zmin*zmin)
				dpmax=np.sqrt(xmax*xmax+ymax*ymax+zmax*zmax)
				dp=(dpmin+dpmax)/2.
				x=(xmin+xmax)/2.
				y=(ymin+ymax)/2.
				z=(zmin+zmax)/2.
				dec=np.arcsin(z/dp)
				ra=np.arctan2(y,x)
				pix=hp.pixelfunc.ang2pix(NSIDE,-dec+0.5*np.pi,ra)
				W[ipos]=f1[pix]-Wtot[ipos]
	Wtot+=W
	return W

def main(char* systfilename, char* label, char* rhorfilename, double minfsys, double maxfsys, char[:,:] files, char[:,:] randoms, double redshiftmin=0.43, double redshiftmax=0.7, int totdeg=3, char* GalaxyvsContaminantfile=None):
	cdef np.ndarray fsys, mask
	fsys, mask = loadsyst(systfilename)
	cdef int deg, rhobar, ix, iy, iz, ipos, pix
	cdef np.ndarray rho_r=np.loadtxt(rhorfilename)
	cdef np.ndarray data1
	if GalaxyvsContaminantfile:
		data1=np.loadtxt(GalaxyvsContaminantfile)
	else:
		data1=GalaxyvsContaminant.main(fsys, mask, minfsys, maxfsys, files, randoms, redshiftmin, redshiftmax)
	cdef np.ndarray Wtot=np.ones(NX*NY*NZ)
	cdef np.ndarray f=np.zeros((NX*NY*NZ,totdeg))
	cdef np.ndarray W, fwbar, pos, p1
	cdef double Wrhobar_over_rhobar
	for deg in np.arange(1,totdeg+1):
		fwbar=np.zeros(NX*NY*NZ)
		pos=np.zeros(NX*NY*NZ)
		p1=np.polyfit(data1[:,0],data1[:,1],deg,w=1./data1[:,2])
		W=polynomial2invweight(p1, Wtot)
		Wrhobar_over_rhobar=np.dot(rho_r*W)/np.sum(rho_r)
		for ipos in np.arange(NX*NY*NZ):
			f[ipos][deg-1]=(W[ipos]-Wrhobar_over_rhobar)*rho_r[ipos]
	cdef np.ndarray fout=np.zeros(NX*NY*NZ)
	for deg in np.arange(1,totdeg+1):
		for ipos in np.arange(NX*NY*NZ):
			fout[ipos]=f[ipos][deg-1]*Wtot[ipos]
		np.savetxt("f"+label+str(deg)+"_NGC_BOX128.txt",f[:,deg-1])
