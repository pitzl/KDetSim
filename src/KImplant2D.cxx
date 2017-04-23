#include "KImplant2D.h"
#include "TMath.h"

ClassImp(KImplant2D)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KImplant2D                                                           //
// A mesh generator for the non-equividistant bins                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


Double_t KImplant2D::ImplEdge(Double_t *x,Double_t *par)
{
  // x  = coordinate x
  // par[0]= curvature in X;
  // par[1]= curvature in Y;

  if(x[0]>par[0]) return -1;  
  return par[1]*TMath::Power(1-TMath::Power(TMath::Abs(x[0]/par[0]),par[2]), 1/par[2] );  
}

Double_t KImplant2D::Derivative(Double_t *x,Double_t *par)
{
  // x  = coordinate x
  // par[0]= width in X;
  // par[1]= width in Y;
  // par[2]= curvature
  if(x[0]>par[0]) return 1;  
  return (-par[1]/TMath::Power(par[0],par[2])*TMath::Power(1-TMath::Power(x[0]/par[0],par[2]), 1/par[2]-1 )*TMath::Power(x[0],par[2]-1));  
}

Double_t KImplant2D::Distance(Double_t *x,Double_t *par, Double_t *X)
{
  // par[0]= width in X;
  // par[1]= width in Y;
  // par[2]= curvature in XY;


 Double_t yl,yr,ym;
 Double_t xl=0,xr=par[0],xm;
 Double_t y0,dist;
 Int_t d=0;


 xm=(xr+xl)/2.;
 yl=TMath::Tan(TMath::ATan(Derivative(&xl,par))+TMath::Pi()/2)*(x[0]-xl)+ImplEdge(&xl,par)-x[1];
 yr=TMath::Tan(TMath::ATan(Derivative(&xr,par))+TMath::Pi()/2)*(x[0]-xr)+ImplEdge(&xr,par)-x[1];
 ym=TMath::Tan(TMath::ATan(Derivative(&xm,par))+TMath::Pi()/2)*(x[0]-xm)+ImplEdge(&xm,par)-x[1];

 while ((TMath::Abs(ym)>1e-5 && TMath::Abs(xl-xr)>1e-15) && d<50 )
   {
     if(yl*ym<0) xr=xm; else xl=xm;
     xm=(xl+xr)/2;
     
     yl=TMath::Tan(TMath::ATan(Derivative(&xl,par))+TMath::Pi()/2)*(x[0]-xl)+ImplEdge(&xl,par)-x[1];
     yr=TMath::Tan(TMath::ATan(Derivative(&xr,par))+TMath::Pi()/2)*(x[0]-xr)+ImplEdge(&xr,par)-x[1];
     ym=TMath::Tan(TMath::ATan(Derivative(&xm,par))+TMath::Pi()/2)*(x[0]-xm)+ImplEdge(&xm,par)-x[1];

     //     printf("%d::  %f %f (%f %f, %f %f)\n",d,xm,ym, xl,yl,xr,yr);
     d++;
   }

 y0=ImplEdge(&xm,par);
 dist=TMath::Sqrt(TMath::Power(xm-x[0],2)+TMath::Power(y0-x[1],2)); 

 if(X!=NULL)
   {
 X[0]=xm;
 X[1]=y0;
   }

 if(xm>x[0]) dist=-dist;
 return dist;
}

Double_t KImplant2D::Distance1(Double_t *x,Double_t *par, Double_t *X)
{
  // par[0]= curvature in Y;
  // par[1]= curvature in X;
  // par[2]= Point to which the distance is measured X1;
  // par[3]= Point to which the distance is measured Y1;

 Double_t yl[2],yr[2];
 Double_t vl,vr;
 Double_t dist;
 Int_t d=0;

 yl[0]=0; yr[0]=par[0];
 yl[1]=ImplEdge(yl[0],par); yr[1]=ImplEdge(yr[0],par); 
 vl=PDistance(x,yl); vr=PDistance(x,yr); 
 printf("%d :: %f %f %f %f , %f %f\n",d,yl[0],yr[0],yl[1],yr[1],vl,vr);
 while ((TMath::Abs(vl-vr)>1e-10 && TMath::Abs(yl[0]-yr[0])>1e-15) && d<50 )
   {
     if(vl<vr) 
       { yr[0]=(yl[0]+yr[0])*0.5;  yr[1]=ImplEdge(yr[0],par); vr=PDistance(x,yr); }
     else 
       { yl[0]=(yl[0]+yr[0])*0.5;  yl[1]=ImplEdge(yl[0],par); vl=PDistance(x,yl); }
 printf("%d :: %f %f %f %f , %f %f\n",d,yl[0],yr[0],yl[1],yr[1],vl,vr);

     d++;
   }

 dist=PDistance(x,yl);

 if(X!=NULL)
   {
 X[0]=yl[0];
 X[1]=yl[1];
   }

 if(yl[0]>x[0]) dist=-dist;
 return dist;
}



Double_t KImplant2D::Conc(Double_t *x, Double_t Thresh)
{
  // x[0]  = coordinate x
  // x[1]  = coordinate y 

  Double_t dist=Distance(x,Dim); 
  if(dist<Thresh) return 0; else
  return fConc->Eval(dist);
  
}


KImplant2D::KImplant2D(Double_t *x, Double_t Sigma, Double_t Nimpl)
{
  // x[0]  = coordinate x
  // x[1]  = coordinate y 

  for(Int_t i=0;i<3;i++) Dim[i]=x[i];
  fConc=new TF1("fConc","TMath::Erfc((x-[0])/[1])*[2]+[3]",-20,20);
  fConc->SetParameter(3,0);
  fConc->SetParameter(2,Nimpl);
  fConc->SetParameter(1,Sigma);
  fConc->SetParameter(0,0);

}
