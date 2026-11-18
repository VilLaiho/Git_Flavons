#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */
extern int access(const char *pathname, int mode);

int nModelParticles=44;
static ModelPrtclsStr ModelPrtcls_[44]=
{
  {"a","a",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","WZ",2,1,3,0}
, {"W+","W-",0, 24, "MW","WW",2,1,3,3}
, {"g","g",1, 21, "0","0",2,8,16,0}
, {"ve","ve~",0, 12, "0","0",1,1,2,0}
, {"vm","vm~",0, 14, "0","0",1,1,2,0}
, {"vt","vt~",0, 16, "0","0",1,1,2,0}
, {"e-","e+",0, 11, "Me","0",1,1,2,-3}
, {"mu-","mu+",0, 13, "MMU","0",1,1,2,-3}
, {"ta-","ta+",0, 15, "MTA","0",1,1,2,-3}
, {"u","u~",0, 2, "MU","0",1,3,6,2}
, {"c","c~",0, 4, "MC","0",1,3,6,2}
, {"t","t~",0, 6, "MT","WT",1,3,6,2}
, {"d","d~",0, 1, "MD","0",1,3,6,-1}
, {"s","s~",0, 3, "MS","0",1,3,6,-1}
, {"b","b~",0, 5, "MB","0",1,3,6,-1}
, {"H","H",1, 25, "MH","WH",0,1,1,0}
, {"X11","X11",1, 10137, "mX11","WX11",0,1,1,0}
, {"X12","X12",1, 10012, "mX12","WX12",0,1,1,0}
, {"X13","X13",1, 10013, "mX13","WX13",0,1,1,0}
, {"X21","X21",1, 10021, "mX21","WX21",0,1,1,0}
, {"X22","X22",1, 10022, "mX22","WX22",0,1,1,0}
, {"X23","X23",1, 10023, "mX23","WX23",0,1,1,0}
, {"X31","X31",1, 10031, "mX31","WX31",0,1,1,0}
, {"X32","X32",1, 10032, "mX32","WX32",0,1,1,0}
, {"X33","X33",1, 10033, "mX33","WX33",0,1,1,0}
, {"~ZL11","~ZL11",1, 10000, "mZL11","WZL11",0,1,1,0}
, {"~ZL12","~ZL12",1, 10001, "mZL12","WZL12",0,1,1,0}
, {"~ZL13","~ZL13",1, 10002, "mZL13","WZL13",0,1,1,0}
, {"~ZL21","~ZL21",1, 10003, "mZL21","WZL21",0,1,1,0}
, {"~ZL22","~ZL22",1, 10004, "mZL22","WZL22",0,1,1,0}
, {"~ZL23","~ZL23",1, 10005, "mZL23","WZL23",0,1,1,0}
, {"~ZL31","~ZL31",1, 10006, "mZL31","WZL31",0,1,1,0}
, {"~ZL32","~ZL32",1, 10007, "mZL32","WZL32",0,1,1,0}
, {"~ZL33","~ZL33",1, 10008, "mZL33","WZL33",0,1,1,0}
, {"~ZR11","~ZR11",1, 10009, "mZR11","WZR11",0,1,1,0}
, {"~ZR12","~ZR12",1, 10010, "mZR12","WZR12",0,1,1,0}
, {"~ZR13","~ZR13",1, 10011, "mZR13","WZR13",0,1,1,0}
, {"~ZR21","~ZR21",1, 10087, "mZR21","WZR21",0,1,1,0}
, {"~ZR22","~ZR22",1, 10078, "mZR22","WZR22",0,1,1,0}
, {"~ZR23","~ZR23",1, 10014, "mZR23","WZR23",0,1,1,0}
, {"~ZR31","~ZR31",1, 10015, "mZR31","WZR31",0,1,1,0}
, {"~ZR32","~ZR32",1, 10016, "mZR32","WZR32",0,1,1,0}
, {"~ZR33","~ZR33",1, 10017, "mZR33","WZR33",0,1,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=69;
int nModelFunc=2;
static int nCurrentVars=68;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[71]={
 "aEWM1","Gf","aS","ymdo","ymup","yms","ymc","ymb","ymt","yme"
,"ymm","ymtau","cabi","lP","lL1","lR1","gLR","gLL","gRR","gLP"
,"gRP","K1","K2","K3","K4","K5","K6","K7","K8","MZ"
,"Me","MMU","MTA","MU","MC","MT","MD","MS","MB","MH"
,"mX11","mX12","mX13","mX21","mX22","mX23","mX31","mX32","mX33","mZL11"
,"mZL12","mZL13","mZL21","mZL22","mZL23","mZL31","mZL32","mZL33","mZR11","mZR12"
,"mZR13","mZR21","mZR22","mZR23","mZR31","mZR32","mZR33","E","Pi","aEW"
,"MW"};
char**varNames=varNames_;
static REAL varValues_[71]={
   1.279000E+02,  1.166370E-05,  1.184000E-01,  5.040000E-03,  2.550000E-03,  1.010000E-01,  1.270000E+00,  4.700000E+00,  1.720000E+02,  5.110000E-04
,  1.056600E-01,  1.777000E+00,  2.277360E-01,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00
,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  9.118760E+01
,  5.110000E-04,  1.056600E-01,  1.777000E+00,  2.550000E-03,  1.270000E+00,  1.720000E+02,  5.040000E-03,  1.010000E-01,  4.700000E+00,  1.250000E+02
,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02
,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02
,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.550000E+02,  2.718282E+00,  3.141593E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   nCurrentVars=69;
   V[69]=Pow(V[0],-1);
   if(!isfinite(V[69]) || FError) return 69;
   nCurrentVars=70;
   V[70]=Pow(Pow(V[29],2)/(2.)+Pow(Pow(V[29],4)/(4.)-V[69]*V[68]*Pow(2,-0.5)*Pow(V[1],-1)*Pow(V[29],2),0.5),0.5);
   if(!isfinite(V[70]) || FError) return 70;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
