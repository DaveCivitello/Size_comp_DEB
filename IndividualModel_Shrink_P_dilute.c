/* file IndividualModel_Shrink.c */
#include <R.h>
#include <math.h>
static double parms[25];

#define iM parms[0]
#define k parms[1]
#define M parms[2]
#define EM parms[3]
#define Fh parms[4]
#define muD parms[5]
#define DR parms[6]
#define fe parms[7]
#define yRP parms[8]
#define ph parms[9]
#define yPE parms[10]
#define iPM parms[11]
#define eh parms[12]
#define mP parms[13]
#define alpha parms[14]
#define yEF parms[15]
#define LM parms[16]
#define kR parms[17]
#define d0 parms[18]
#define kk parms[19]
#define hb parms[20]
#define theta parms[21]
#define mR parms[22]
#define yVE parms[23]
#define startage parms[24]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=25;
odeparms(&N, parms);
}

/* Derivatives and 2 output variables */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");

double Chi = M/(1 + EM);
double fH = y[0]/(y[0]+Fh);
double fP = y[2]/(y[2] + eh);

double L = y[1];
double LG = fmax(y[1], yout[1]);

double GVOL = pow(LG, 3);
double VOL = pow(L,3);
double SA = pow(L,2);
double fbf = 0;

if( t[0] < startage){fH = fe;}else{fbf = 1;};
double Dens = y[5]/(Chi*GVOL);
double kstar = fmin(k + y[5]*alpha, 1);
double g = 1/(yVE*kstar*EM);
double aM = iM*yEF;
double mV = aM*k/(LM*Chi);
double mD = muD*mV;
double rp = Dens*Dens/(ph*ph + Dens*Dens);
double Jec = y[2]*g/(g + y[2])*(aM*SA + yVE*EM*(mV+mR*EM*y[7])*Chi*VOL);

ydot[0] = -iM*SA*fH*fbf;
ydot[1] = yVE/(3*Chi*SA)*(kstar*Jec - (mV+mR*EM*y[7])*Chi*VOL);
ydot[2] = aM/(Chi*EM*L)*(fH - y[2]) - iPM*y[5]*fP/(EM*Chi*VOL);
ydot[5] = yPE*iPM*fP*(1 - rp)*y[5] - mP*y[5];
ydot[6] = fmax(yRP*yPE*iPM*fP*rp*y[5],0);
if(y[3] < DR){
  ydot[3] = (1 - kstar)*Jec - mD*y[3];
  ydot[4] = 0;}else{
  ydot[3] = fmin(0, (1 - kstar)*Jec - mD*DR);
  ydot[4] = fmax((1 - kstar)*Jec - mD*DR, 0);}
ydot[7] = theta/(Chi*VOL)*ydot[6] + kR*(1-y[2]) - kR*y[7] - 3*y[7]*ydot[1]/L;
ydot[8] = kk*fmax(y[7] - d0, 0) + hb;
if(y[2] <= 0){
  ydot[0] = 0;
  ydot[1] = 0;
  ydot[2] = 0;
  ydot[3] = 0;
  ydot[4] = 0;
  ydot[5] = 0;
  ydot[6] = 0;
  ydot[7] = 0;
  ydot[8] = 100;}

  yout[0] = exp(-y[8]);
  yout[1] = LG;

}

/* END file IndividualModel_Shrink.c */ 
