/*  Routine to read in galaxy catalogue files 
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

#include "header.h"

double calc_dp(double);

// ************************************************************************
// READ IN CATALOGUES. Various formats for columns!
//
void read_gal_file(struct use_gal *gal,char *fname,long *nmax,FILE *fout,double zmin,double zmax,int fflag) {

	printf("%s\n",fname);

	// open file and go through lines
	FILE *fp;
	if((fp=fopen(fname,"r"))==NULL) err_handler("cannot open galaxy file");
	const int bsz=400; char buf[bsz];

	long count_gal=0;
	while(fgets(buf,bsz,fp)) if(strncmp(buf,"#",1)!=0 && strncmp(buf+1,"#",1)!=0) {

		double ra,dec,red=-1.,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
		int tflag=0, nread=0;

		// catalog format: "ra,dec,z,wght_fkp,wght_cp,wght_noz,wght_sdc"
		if(fflag==0) {
			nread = sscanf(buf,"%le %le %le %le %le %le %le",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);
		}

		// catalog format: "ra,dec,z"
		if(fflag==1) nread = sscanf(buf,"%lf %lf %lf",&ra,&dec,&red);

		// catalog format: "ra,dec,z,wght_fkp,z,wght_fkp", with the second z chosen
		if(fflag==2) nread = sscanf(buf,"%lf %lf %*f %*f %lf %lf",&ra,&dec,&red,&tmp1);

		// catalog format: "ra,dec,zflag,z,wght_cp,wght_noz,wght_star"
		if(fflag==3) nread = sscanf(buf,"%lf %lf %d %lf %lf %lf %lf",&ra,&dec,&tflag,&red,&tmp1,&tmp2,&tmp3);

		// catalog format: "ra,dec,z,ipoly,wght_fkp,wght_cp,wght_noz,wght_veto"
		if(fflag==4) nread = sscanf(buf,"%lf %lf %lf %*d %lf %lf %lf %lf",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);

		// catalog format: "ra,dec,z,wght", where wght is the combined weight
		if(fflag==5) nread = sscanf(buf,"%lf %lf %lf %lf",&ra,&dec,&red,&tmp1);

		// catalog format: "ra,dec,z,wght", where wght is the combined weight
		if(fflag==6) nread = sscanf(buf,"%lf %lf %lf %lf %lf",&ra,&dec,&red,&tmp1,&tmp2);

		// catalog format: "ra,dec,z,ipoly,wght_fkp,wght_cp,wght_noz,wght_veto"
		if(fflag==7) nread = sscanf(buf,"%lf %lf %lf %*d %lf %lf %lf %*f %*d %*f %*f %lf",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);

		// catalog format: "ra,dec,z,wght,wght_fkp,nbar"
		if(fflag==8) {
			nread = sscanf(buf,"%lf %lf %lf %lf %lf %lf",&ra,&dec,&red,&tmp1,&tmp2,&tmp3);
			if(nread==5) {
				tmp1=1.0;
				tmp3=tmp2;
			}
		}

		// catalogue format: "ra,dec,z,weight_fkp,weight_cp,weight_noz,weight_star,weight_seeing,weight_systot"
		if (fflag==10) {
			nread=sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&tmp6);
		}

		// catalogue format: "ra,dec,z,nbar,bias,vetoflag,fibercollision"
		if (fflag==11) {
			nread=sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);
		}

		// catalogue format: "ra,dec,z,weight_fkp"
		if (fflag==12) {
			nread=sscanf(buf,"%le %le %le %le",&ra,&dec,&red,&tmp1);
		}

    // catalog format: "ra,dec,z,wght_fkp,wght_cp,wght_noz,wght_sdc" but ignoring wght_sdc
    if(fflag==13){
      nread = sscanf(buf,"%le %le %le %le %le %le %le",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);
		  tmp4=1;
    }

		// Patchy_V6 data format
		if(fflag==14){
			nread = sscanf(buf,"%le %le %le %le %le %le %le %le",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
		}

		// Patchy_V6 random format
		if(fflag==15){
			nread = sscanf(buf,"%le %le %le %le %le %le %le",&ra,&dec,&red,&tmp1,&tmp2,&tmp3,&tmp4);
		}

		if(red>=zmin && red<=zmax && (fflag!=3 || (tflag==1 || tflag==2)) && (fflag!=4 || (tmp1>0.0 && tmp2>0.0 && tmp3>0.0 && tmp4>0.1)) &&(fflag!=7 || (tmp1>0.0 && tmp2>0.0 && tmp3>0.0 && tmp4>0.1))) {

			if(++count_gal>(*nmax)) {
				count_gal--;
				break;
			}

			if(fflag!=8) {
				ra   = ra * pi/180.;
				dec  = dec * pi/180.;
			}
			double dp   = calc_dp(red);
			gal[count_gal].x  = dp*cos(dec)*cos(ra);
			gal[count_gal].y  = dp*cos(dec)*sin(ra);
			gal[count_gal].z  = dp*sin(dec);

			gal[count_gal].wght = 1.0;
			if(fflag==0 || fflag==13) {
				gal[count_gal].wght = tmp4*(tmp2+tmp3-1.0);
				gal[count_gal].nbar = (1./tmp1-1.)/Pfkp;
			}
			if(fflag==4 || fflag==7) {
				gal[count_gal].wght = (tmp2+tmp3-1.0);
			}
			if(fflag==2) {
				gal[count_gal].wght = tmp1;
				gal[count_gal].nbar = (1./tmp1-1.)/Pfkp;
			}
			if(fflag==3) {
				gal[count_gal].wght = tmp3*(tmp1+tmp2-1.0);
			}
			if(fflag==5) {
				gal[count_gal].wght = tmp1;
			}
			if(fflag==6) {
				gal[count_gal].wght = tmp1;
				gal[count_gal].nbar = tmp2;
			}
			if(fflag==8) {
				gal[count_gal].wght = tmp1;
				gal[count_gal].nbar = tmp3;
			}
			if (fflag==10) {
				gal[count_gal].wght=tmp6*(tmp2+tmp3-1.);
				gal[count_gal].nbar=(1./tmp1-1.)/Pfkp;
			}
			if (fflag==11) {
				gal[count_gal].wght=tmp3*tmp4;
				gal[count_gal].nbar=tmp1;
			}
			if (fflag==12) {
				gal[count_gal].nbar=(1./tmp1-1.)/Pfkp;
			}
			if (fflag==11) {
				gal[count_gal].wght=tmp3*tmp4;
				gal[count_gal].nbar=tmp1;
			}
			if (fflag==14) {
				gal[count_gal].wght=tmp4*tmp5;
				gal[count_gal].nbar=tmp2;
			}
			if (fflag==15) {
				gal[count_gal].wght=tmp3*tmp4;
				gal[count_gal].nbar=tmp1;
			}

			// catalog format: "x,y,z,wght"
			if(fflag==9) {
				nread = sscanf(buf,"%lf %lf %lf %lf",&tmp1,&tmp2,&tmp3,&tmp4);
				if(++count_gal>(*nmax)) {
					count_gal--;
					break;
				}
				gal[count_gal].x  = tmp1;
				gal[count_gal].y  = tmp2;
				gal[count_gal].z  = tmp3;
				gal[count_gal].wght = tmp4;
			}
		}
	}
	fprintf(stderr,"read in %ld galaxies, maximum %ld\n",count_gal,(*nmax));
	(*nmax) = count_gal;
}
