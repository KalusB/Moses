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

def main(char* systfilename, char* label, char* rhorfilename, double minfsys, double maxfsys, char[:,:] files, char[:,:] randoms, double redshiftmin=0.43, double redshiftmax=0.7, int totdeg=3):
	cdef np.ndarray fsys=np.loadtxt(systfilename)
	cdef np.ndarray mask=np.zeros(len(fsys))
	cdef int cell, neighbour, deg, rhobar, ix, iy, iz, ipos, pix
	for cell in np.arange(len(fsys)):
		if fsys[cell]==0:
			mask[cell]=0
		else:
			mask[cell]=1
			for neighbour in hp.get_interp_weights(NSIDE,cell,nest=True)[0]:
				if fsys[neighbour]==0:
					mask[cell]-=0.25
	fsys/=mask
	cdef np.ndarray rho_r=np.loadtxt(rhorfilename)
	cdef np.ndarray data1=GalaxyvsContaminant.main(fsys, mask, minfsys, maxfsys, files, randoms, redshiftmin, redshiftmax)
	f1=[1]
	cdef np.ndarray Wtot=np.ones(NX*NY*NZ)
	cdef np.ndarray f=np.zeros((NX*NY*NZ,totdeg))
	cdef np.ndarray W, fwbar, pos, p1
	cdef double xmin, xmax, ymin, ymax, zmin, zmax, dpmin, dpmax, dp, x, y, z, dec, ra, Wrhobar_over_rhobar
	for deg in np.arange(1,totdeg+1):
		W=np.zeros(NX*NY*NZ)
		fwbar=np.zeros(NX*NY*NZ)
		pos=np.zeros(NX*NY*NZ)
		Wrhobar=0
		rhobar=0
		p1=np.polyfit(data1[:,0],data1[:,1],deg,w=1./data1[:,2])
		f1.append(np.poly1d(p1)(fsys))
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
					W[ipos]=f1[deg][pix]-Wtot[ipos]
					Wtot[ipos]+=W[ipos]
					Wrhobar+=rho_r[ipos]*W[ipos]
					rhobar+=rho_r[ipos]
		Wrhobar_over_rhobar=Wrhobar/rhobar
		for ipos in np.arange(NX*NY*NZ):
			f[ipos][deg-1]=(W[ipos]-Wrhobar_over_rhobar)*rho_r[ipos]
		print "order "+str(deg)+"/"+str(totdeg)+" complete"
	cdef np.ndarray fout=np.zeros(NX*NY*NZ)
	for deg in np.arange(1,totdeg+1):
		for ipos in np.arange(NX*NY*NZ):
			fout[ipos]=f[ipos][deg-1]*Wtot[ipos]
		np.savetxt("f"+label+str(deg)+"_NGC_BOX128.txt",f[:,deg-1])
