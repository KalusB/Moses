'''
    Programme to compute the ratio of observed galaxies to expected 
    galaxies in contaminant bins
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

DTYPE = np.int
ctypedef np.int_t DTYPE_t

# Find HealPix index corresponding to Equatorial coordinates
def DeclRaToIndex(double decl,double RA):
	return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))

cdef int NSIDE = 64
cdef int numbins=10	# number of contaminant bins

# Contaminant bin corresponding to contamination level f
def binid(double f, double minfsys, double maxfsys):
	cdef int bid=int(np.nan_to_num(float(f-minfsys)/(maxfsys-minfsys))*numbins)
	if bid<0 or bid>numbins:
		return numbins
	else:
		return bid

#main routine, fsys is a healpix map of the contaminant, mask is an auxiliary map providing the survey footprint, minfsys and maxfsys correspond to the contaminant values of the minimum and maximum bin, Galaxyfiles is a list containing the filenames of galaxy catalogue FITS files (can be more than one to allow e.g. for data on different hemispheres), randoms provides the filenames for the corresponding random catalogues, redshiftmin and redshiftmax are the minimal and maximal redshift of the survey
def main(np.ndarray fsys, np.ndarray mask, double minfsys, double maxfsys, np.ndarray Galaxyfiles, np.ndarray randoms, double redshiftmin=0.43, double redshiftmax=0.7):
	cdef int cell
	cdef np.ndarray nsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray nbarbynsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray sigmanbarbynsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray ngal=np.zeros(len(Galaxyfiles))

	cdef int i, nonzerocells, bin
	cdef double ng, nr, nbarsum, nbarrsum, alpha
	cdef np.ndarray Galaxydata, Randomdata, nbar, nbarr, nbarb, nbarrb, numfsys

	for i in np.arange(len(Galaxyfiles)):
		ng=0
		Galaxydata = fitsio.read(Galaxyfiles[i])
		print "opened "+Galaxyfiles[i]
		nbar=np.zeros(hp.nside2npix(NSIDE))
		nbarr=np.zeros(hp.nside2npix(NSIDE))
		nbarb=np.zeros(numbins)
		nbarrb=np.zeros(numbins)
		nonzerocells=0
		nr=0
		for gal in Galaxydata:
			if gal['Z']>redshiftmin and gal['Z']<redshiftmax:
				nbar[int(DeclRaToIndex(gal[1],gal[0]))]+=(gal['WEIGHT_NOZ']+gal['WEIGHT_CP']-1.)*gal['WEIGHT_FKP']*gal['WEIGHT_SEEING']
				ng+=(gal['WEIGHT_NOZ']+gal['WEIGHT_CP']-1.)*gal['WEIGHT_FKP']*gal['WEIGHT_SEEING']
		nbarsum=0
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			nbarsum+=nbar[cell]
		print "ng:",ng,"sum of galaxies in cells:",nbarsum
		print "read "+Galaxyfiles[i]
		ngal[i]=Galaxydata.shape[0]
		del Galaxydata
		Randomdata = fitsio.read(randoms[i])
		print "opened "+randoms[i]
		for ran in Randomdata:
			if ran['Z']>redshiftmin and ran['Z']<redshiftmax:
				nbarr[int(DeclRaToIndex(ran[1],ran[0]))]+=ran['WEIGHT_FKP']
				nr+=ran['WEIGHT_FKP']
		nbarrsum=0
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			nbarrsum+=nbarr[cell]
		print "nr:",nr,"sum of randoms in cells:",nbarrsum
		print "read "+randoms[i]
		numfsys=np.zeros(numbins)
		alpha=nr/ng
		del Randomdata
		print "ratio of randoms per galaxy:",alpha
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			if binid(fsys[cell])<numbins:
				nbarb[binid(fsys[cell])]+=nbar[cell]
				nbarrb[binid(fsys[cell])]+=nbarr[cell]
				nsys[binid(fsys[cell]),i]+=fsys[cell]
				numfsys[binid(fsys[cell])]+=1
		for bin in range(0,numbins):
			if numfsys[bin]>0:
				nbarbynsys[bin,i]=alpha*nbarb[bin]/nbarrb[bin]
				sigmanbarbynsys[bin,i]=np.sqrt(nbarb[bin])/nbarrb[bin]*alpha
				nsys[bin,i]/=numfsys[bin]
				print nsys[bin,:], nbarbynsys[bin,:], sigmanbarbynsys[bin,:],numfsys[bin]
	cdef np.ndarray nsystot=np.zeros(numbins)
	cdef np.ndarray nbarbynsystot=np.zeros(numbins)
	cdef np.ndarray sigmanbarbynsystot=np.zeros(numbins)
	cdef double N=0
	for i in np.arange(len(Galaxyfiles)):
		nsystot+=nsys[:,i]*ngal[i]
		nbarbynsystot+=nbarbynsys[:,i]*ngal[i]
		sigmanbarbynsystot+=sigmanbarbynsys[:,i]*ngal[i]
		N+=ngal[i]
		print ngal[i],nsys[1,i],nsystot[1]/N,nbarbynsys[1,i],nbarbynsystot[1]/N,sigmanbarbynsys[1,i],sigmanbarbynsystot[1]/N
	nsystot/=N
	nbarbynsystot/=N
	sigmanbarbynsystot/=N
	return np.transpose([nsystot,nbarbynsystot,sigmanbarbynsystot])
