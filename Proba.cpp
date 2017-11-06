#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
//#include <Values.h>
#include <string.h>
#define F_real double_
#include "doub.h"
#include "stype.h"
#include "gspuser.h"
int FORSAG(F_real *P,F_real *SD,F_real *XLFX,F_real *GG,F_real *T,F_real *TF,F_real *ALFA,F_real *PDF,F_real *GTS,F_real *PFX,F_real *FFX,F_real *PF,F_real *GTF,F_real *QTF,F_real *ALFAF);
int IT(F_real *T,
       F_real *ITResult);

int AIT(F_real *TA,
	F_real *AL,
	F_real *Z0,
	F_real *AITResult);

int RCKM(F_real *A,
	    F_real *T,
	    F_real *R,
	    F_real *XK,
	    F_real *XM);

int QLAM(F_real *XL,
	 F_real *XK,
	 F_real *RES);

int vid(VECT *X,VECT *L,double *a,double *d);


int ogr_raz(double, double* ,int );
int t =15,r=0;
    double_ operator +(double_ c1, double_ c2)
     {
       double d1,d2,ma;int i;
       ogr_raz(c1.a,&d1,t); ma=fabs(d1);
       ogr_raz(c2.a,&d2,t);
       if (fabs(d2)>ma) {ma=fabs(d2);}
       for (i=0;i<t-r;i++) ma=ma/10.0;
       d1=d1+d2;
       if (ma>fabs(d1)) {return double_(0.0);}
       ogr_raz(d1,&d2,t);
       return double_(d2);
    }

     double_ operator -(double_ c1, double_ c2)
     {
       double d1,d2,ma;int i;
       ogr_raz(c1.a,&d1,t);   ma=fabs(d1);
       ogr_raz(c2.a,&d2,t);
       if (fabs(d2)>ma) {ma=fabs(d2);}
       for (i=0;i<t-r;i++) ma=ma/10.0;
       d1=d1-d2;
       if (ma>fabs(d1)) {return double_(0.0);}
       ogr_raz(d1,&d2,t);
       return double_(d2);
     }


double_ operator /(double_ c1, double_ c2)
      {
       double d1,d2;
       ogr_raz(c1.a,&d1,t);
       ogr_raz(c2.a,&d2,t);   d1=d1/d2;
       ogr_raz(d1,&d2,t);
       return double_(d2);
      }

 double_ operator *(double_ c1, double_ c2)
      {
       double d1,d2;
       ogr_raz(c1.a,&d1,t);
       ogr_raz(c2.a,&d2,t);   d1=d1*d2;
       ogr_raz(d1,&d2,t);
       return double_(d2);
      }

int ogr_raz(double c1, double* c2,int t)
   {
       char *endptr;
       char *str, str1[25],str2[5];
       int dec1, sign1=0, ndig =15;
 strcpy(str1,"0.");
 str=ecvt(c1, ndig, &dec1, &sign1);
 strncat(str1,str,t);  strncat(str1,"e",1);
 itoa(dec1, str2,10);
 strcat(str1,str2);
 *c2 = strtod(str1, &endptr);
 if (sign1 !=0) *c2=-*c2;
 return 1;
 }

 int bool2(double *Xc0, double *Yc0)
{
F_real P,SD,XLFX,GG,T,TF,ALFA,PDF,GTS,PFX,FFX,PF,GTF,QTF,ALFAF;

P=150.0;
SD=0.95;
XLFX=0.1;
GG=72.0;
T=800.0;
TF=1500.0;
ALFA=1.0;
PDF=0.98;
GTS=50.0;
TF=*Xc0; T=*Yc0;

FORSAG(&P,&SD,&XLFX,&GG,&T,&TF,&ALFA,&PDF,&GTS,&PFX,&FFX,&PF,&GTF,&QTF,&ALFAF);
return 1;
 }
//F_real AL,Z0,R,XM,XK,HU,AIH,AIK,QTSM,ZFX,AIKG,RF,XKF,XMF,ZLF,XLF,MM,QL;
/*   основная функция   */
int FORSAG(F_real *P,
	    F_real *SD,
	    F_real *XLFX,
	    F_real *GG,
	    F_real *T,
	    F_real *TF,
	    F_real *ALFA,
	    F_real *PDF,
	    F_real *GTS,
	    F_real *PFX,
	    F_real *FFX,
	    F_real *PF,
	    F_real *GTF,
	    F_real *QTF,
	    F_real *ALFAF)

{
F_real AL,Z0,R,XM,XK,HU,AIH,AIK,QTSM,ZFX,AIKG,RF,XKF,XMF,ZLF,XLF,MM,QL;
F_real ZLF1,YF,XLFX1,sqr1;
//FILE  *out;

 AL=1.0;
  (*PFX)=(*P)* (*SD);
  Z0=14.78;
  if ((*ALFA).a==0)
  {
    R=29.27;
    XM=0.6847;
    XK=1.4;
    goto VOZD;
  }
  RCKM(ALFA,T,&R,&XK,&XM);
 VOZD:
  ZFX=0.5*((*XLFX)+1/ (*XLFX));
  QLAM(XLFX,&XK,&QL);

  if(R.a * (*T).a/9.81 > 0)
	(*FFX)=(*GG)*sqrt(R.a * (*T).a/9.81)/XM/ (*PFX)/QL;
  else
	(*FFX)=(*GG)/XM/ (*PFX)/QL;

  if ((*TF).a>0)
  {
    IT(TF,&AIK);
    IT(T,&AIH);
    AIT(T,&AL,&Z0,&AIKG);
    AIT(TF,&AL,&Z0,&MM);
    HU=10250.0+14.78*(AIH-70.05)-15.78*(AIKG-72.35)+15.0;
    *QTF=(AIK-AIH)/(HU* (*PDF)-15.78*(MM-AIKG)+14.78*(AIK-AIH));
    if ((*ALFA).a==0)
      QTSM=0;
    else
    {
      QTSM=1/ (*ALFA)/14.78;
      *GTF=((*GG)- (*GTS))* (*QTF);
      *ALFAF=((*GG)- (*GTS))/((*GTS)+ (*GTF))/14.78;
      RCKM(ALFAF,TF,&RF,&XKF,&XMF);
	  if((XK.a+1)/XK.a*XKF.a/(XKF.a+1)*R.a/RF.a* (*T).a/(*TF).a > 0)
		ZLF=(1+QTSM)/(1+QTSM+ *QTF)*ZFX*sqrt((XK.a+1)/XK.a*XKF.a/(XKF.a+1)*R.a/RF.a* (*T).a/(*TF).a);
	  else
		ZLF=(1+QTSM)/(1+QTSM+ *QTF)*ZFX;

      if (ZLF.a<1.0)
		{
		  ZLF=1.0;
	
		}
	  if(ZLF.a*ZLF.a-1 > 0)
		XLF=ZLF-sqrt(ZLF.a*ZLF.a-1);
	  else
		XLF=ZLF;
      QLAM(&XLF,&XKF,&QL);
	  if((*TF).a*RF.a/9.81 > 0)
		*PF=((*GG)+ (*GTF))*sqrt((*TF).a*RF.a/9.81)/XMF/ (*FFX)/QL;
	  else
		*PF=((*GG)+ (*GTF))/XMF/ (*FFX)/QL;
    }
  }
  
  return 1;
}

int IT(F_real *T,
       F_real *ITResult)
{ F_real T1,R;
 T1=(*T)/100.0;
 R=(0.24242+T1*(-0.24746e-2+T1*(0.52835e-3+T1*(-0.24661e-4+T1*0.3869e-6))))/(1.0-exp((-1)*(T1.a-10.101)*(T1.a-10.101)));
  *ITResult=(*T)*R;
 return 1;
}


int AIT(F_real *TA,
	F_real *AL,
	F_real *Z0,
	F_real *AITResult)
{
F_real T100,AB,A1,R,B;
T100=(*TA)/100.0;
IT(TA,&AB);
     R=.24303+T100*(.41667E-3+T100*(.36978E-3+T100*(-.1887E-4+T100*.29886E-6)));
     A1=(*TA)*R;
     B=(1.0+(*Z0))/(1.0+(*AL)*(*Z0));
     *AITResult=(1-B)*AB+B*A1;
return 1;
}
int RCKM(F_real *A,
	    F_real *T,
	    F_real *R,
	    F_real *XK,
	    F_real *XM)
{F_real Z0,AITResult1,AITResult2;
 F_real T10;
 T10=(*T)+10.0;
 *R=( (1.0-1.0/(*A)) * 29.27 + (1.0/(*A)/14.78+1.0/(*A))*29.35 ) / (1.0+1.0/(*A)/14.78);
 Z0=14.78;
 AIT(&T10,A,&Z0,&AITResult1);
 AIT(T,   A,&Z0,&AITResult2);
 (*XK)=1.0/(1.0-(*R)/42.7/(AITResult1-AITResult2) );

 if((2./ ((*XK).a+1.0)>0) && ((*XK).a* pow( 2./ ((*XK).a+1.0),((*XK).a+1.0)/((*XK).a-1.0) ) > 0))
	(*XM)=sqrt( (*XK).a* pow( 2./ ((*XK).a+1.0),((*XK).a+1.0)/((*XK).a-1.0) ) );
 else
	(*XM)=1;
 return 1;
}

int QLAM(F_real *XL,
	 F_real *XK,
	 F_real *RES)
{
  if((1- ((*XK).a-1)/((*XK).a+1)*(*XL).a*(*XL).a )>0)
	*RES= (*XL)*pow( fabs(((*XK).a+1)/2),1/((*XK).a-1) ) * pow( (1- ((*XK).a-1)/((*XK).a+1)*(*XL).a*(*XL).a ) , 1/((*XK).a-1) );
  else
	*RES= (*XL)*pow( fabs(((*XK).a+1)/2),1/((*XK).a-1) );
  return 1;
}

/*

   void main()
   {
   int i=0;
   //double_ b(12318), c(12312);
   double a=1600,d=1000;
   //long double exp(long double (x));
   //long double logl(long double (x));
   bool2(&a,&d);

   }

*/