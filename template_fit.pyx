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
import time
import GalaxyvsContaminant

cdef int NSIDE=64
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

def loadsyst(char* systfilename,char* oi='NESTED'):
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
	cdef int NSIDEin=hp.pixelfunc.get_nside(fsys)
	mask=hp.ud_grade(mask,nside_out=NSIDE,order_in=oi, order_out='RING')
	fsys=hp.ud_grade(fsys,nside_out=NSIDE,order_in=oi, order_out='RING')
	fsys*=float(hp.nside2npix(NSIDEin))/hp.nside2npix(NSIDE)
	return fsys, mask

def polynomial2invweight(np.ndarray p1, np.ndarray fsys, np.ndarray Wtot, char* ifib2filename=""):
	cdef np.ndarray W=np.zeros(NX*NY*NZ)
	cdef bint ifib2flag=(len(ifib2filename)>0)
	cdef np.ndarray f1=np.empty([len(p1),len(fsys)]), ifib2av
	if ifib2flag:
		ifib2av=np.loadtxt(ifib2filename)
		for i in np.arange(len(p1)):
			f1[i,:]=np.poly1d(p1[i])(fsys)
	else:
		f1=np.poly1d(p1)(fsys)
	cdef double dp, x, y, z, dec, ra, ifib2
	cdef int ifib2ind
	for ix in np.arange(NX):
		for iy in np.arange(NY):
			for iz in np.arange(NZ):
				ipos = iz+NZ*(iy+NY*ix)
				x=(ix+0.5)*dx+XMIN
				y=(iy+0.5)*dy+YMIN
				z=(iz+0.5)*dz+ZMIN
				dp=np.sqrt(x*x+y*y+z*z)
				dec=np.arcsin(z/dp)
				ra=np.arctan2(y,x)
				pix=hp.pixelfunc.ang2pix(NSIDE,-dec+0.5*np.pi,ra)
				if ifib2flag:
					ifib2=ifib2av[np.searchsorted(ifib2av[1,:],dp),2]
					ifib2ind=np.searchsorted(np.arange(20.3,21.8,0.3),ifib2)
					W[ipos]=f1[ifib2ind,pix]-Wtot[ipos]
				else:
					W[ipos]=f1[pix]-Wtot[ipos]
	Wtot+=W
	return W

def main_ifib2(char* systfilename, char* label, char* rhorfilename, double minfsys, double maxfsys, char[:,:] files, char[:,:] randoms, char* GalaxyvsContaminantfile, char* ifib2filename, double redshiftmin=0.43, double redshiftmax=0.7, int totdeg=3):
	cdef np.ndarray fsys, mask
	fsys, mask = loadsyst(systfilename)
	cdef int deg, rhobar, ix, iy, iz, ipos, pix
	cdef np.ndarray rho_r=np.loadtxt(rhorfilename)
	cdef np.ndarray data1, data2, data3, data4, data5
	data1=np.loadtxt(GalaxyvsContaminantfile[0])
	data2=np.loadtxt(GalaxyvsContaminantfile[1])
	data3=np.loadtxt(GalaxyvsContaminantfile[2])
	data4=np.loadtxt(GalaxyvsContaminantfile[3])
	data5=np.loadtxt(GalaxyvsContaminantfile[4])
	cdef np.ndarray Wtot=np.ones(NX*NY*NZ)
	cdef np.ndarray f=np.zeros((NX*NY*NZ,totdeg))
	cdef np.ndarray W, fwbar, pos, p1
	cdef double Wrhobar_over_rhobar
	for deg in np.arange(1,totdeg+1):
		fwbar=np.zeros(NX*NY*NZ)
		pos=np.zeros(NX*NY*NZ)
		p1=np.array([np.polyfit(data1[:,0],data1[:,1],deg,w=1./data1[:,2]), np.polyfit(data2[:,0],data2[:,1],deg,w=1./data2[:,2]), np.polyfit(data3[:,0],data3[:,1],deg,w=1./data3[:,2]), np.polyfit(data4[:,0],data4[:,1],deg,w=1./data4[:,2]), np.polyfit(data5[:,0],data5[:,1],deg,w=1./data5[:,2])])
		W=polynomial2invweight(p1, Wtot, ifib2filename)
		Wrhobar_over_rhobar=np.dot(rho_r,W)/np.sum(rho_r)
		for ix in np.arange(NX):
			for iy in np.arange(NY):
				for iz in np.arange(NZ):
					ipos = iz+NZ*(iy+NY*ix)
					f[ipos][deg-1]=(W[ipos]-Wrhobar_over_rhobar)*rho_r[ipos]
	cdef np.ndarray fout=np.zeros(NX*NY*NZ)
	for deg in np.arange(1,totdeg+1):
		np.savetxt("f"+label+str(deg)+"_NGC_BOX128.txt",f[:,deg-1])

def main(char* systfilename, char* label, char* rhorfilename, double minfsys, double maxfsys, char[:,:] files, char[:,:] randoms, double redshiftmin=0.43, double redshiftmax=0.7, int totdeg=3, GalaxyvsContaminantfile=None):
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
		Wrhobar_over_rhobar=np.dot(rho_r,W)/np.sum(rho_r)
		for ipos in np.arange(NX*NY*NZ):
			f[ipos][deg-1]=(W[ipos]-Wrhobar_over_rhobar)*rho_r[ipos]
	cdef np.ndarray fout=np.zeros(NX*NY*NZ)
	for deg in np.arange(1,totdeg+1):
		np.savetxt("f"+label+str(deg)+"_NGC_BOX128.txt",f[:,deg-1])
