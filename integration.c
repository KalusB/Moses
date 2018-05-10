/*  Integration routine 
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

#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>

#define QSIMP_EPS     1.0e-7       // accuracy for numerical integration
#define QMIDINF_EPS   1.0e-5       // accuracy for numerical integration
#define QMIDINF2_EPS  1.0e-3       // accuracy for numerical integration
#define QSIMPMID_EPS  1.0e-3       // accuracy for numerical integration

#define QSIMP_MAX     40           // max loops for numerical integration
#define QMIDINF_MAX   40           // max loops for numerical integration
#define QMIDINF2_MAX  20           // max loops for numerical integration
#define QSIMPMID_MAX  20           // max loops for numerical integration

#define FUNC_SIMP(x) ((*func)(x))
#define FUNC_MIDINF(a,b) ((*func)(1.0/(a),(b))/((a)*(a)))
#define FUNC_MIDINF2(a,b,c) ((*func)(1.0/(a),(b),(c))/((a)*(a)))
#define FUNC_SIMPMID(v,w,x) ((*func)(v,w,x))
#define FUNC_SIMPMID2(v,w) ((*func)(v,w))

void err_handler(char*);

double qsimp(double (*func)(double), double a, double b)
{
  int j;
  double s,st,ost,os;
  double trapzd(double (*func)(double), double, double, int);

  ost = os = -1.0e30;
  for (j=1;j<=QSIMP_MAX;j++) {
    st=trapzd(func,a,b,j);
    s=(4.0*st-ost)/3.0;
    if (fabs(s-os) < QSIMP_EPS*fabs(os)) return s;
    if (s == 0.0 && os == 0.0 && j > 6) return s;
    os=s;
    ost=st;
  }
  err_handler("Too many steps in routine qsimp");
  return 0.0;
}

double trapzd(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC_SIMP(a)+FUNC_SIMP(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC_SIMP(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
