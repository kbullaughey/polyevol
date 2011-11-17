#include <valarray>
#include "sim_rand.h"

using std::valarray;

double
ran1() {
  return drand48();
}

#define PI 3.141592654
int
poidev(double xm) {
   float gammln(float);
   static float sq,alxm,g,oldm=(-1.0);
   float em,t,y;
   
   if (xm < 12.0) {
      if (xm != oldm) {
         oldm=xm;
         g=exp(-xm);
      }
      em = -1;
      t=1.0;
      do {
         ++em;
         t *= ran1();
      } while (t > g);
   } else {
      if (xm != oldm) {
         oldm=xm;
         sq=sqrt(2.0*xm);
         alxm=log(xm);
         g=xm*alxm-gammln(xm+1.0);
      }
      do {
         do {
            y=tan(PI*ran1());
            em=sq*y+xm;
         } while (em < 0.0);
         em=floor(em);
         t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (ran1() > t);
   }
   return (int)em;
}
#undef PI

void
ranint(int n, valarray<int> &ranout) {
  int i, spot;
  valarray<int> buf(n);
  for (i=0; i<n; i++) buf[i] = i;

  for (i=0 ; i<(n-1); i++){
    spot = (int)floor( (n-i)*ran1() ) ; 
    ranout[i] = buf[spot];
    buf[spot] = buf[ n-i-1];
  }
  ranout[n-1] = buf[0];
}

float
gammln(float xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

