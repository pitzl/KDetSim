#include "KImplant3D.h"
#include "KImplant2D.h"
#include "TMath.h"

ClassImp(KImplant3D)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KImplant3D                                                           //
// A mesh generator for the non-equividistant bins                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


Double_t KImplant3D::ImplEdge(Double_t *x,Double_t *par)
{
  // par[0]= width in X 
  // par[1]= width in Y
  // par[2]= width in Z
  // par[3]= curvature in X;
  // par[4]= curvature in XZ;
  // par[5]= curvature in YZ;
  Double_t par2D[3];
  Double_t y;
  if(x[1]>par[2]) return -1;
  par2D[0]=par[2];  par2D[1]=par[0];  par2D[2]=par[4];
  Double_t ap=KImplant2D::ImplEdge(x[1],par2D);
  par2D[0]=par[2];  par2D[1]=par[1];  par2D[2]=par[5];
  Double_t bp=KImplant2D::ImplEdge(x[1],par2D);
  //  Double_t ap=par[0]*TMath::Power(1-TMath::Power(x[1]/par[2],par[4]), 1/par[4] );  
  //  Double_t bp=par[1]*TMath::Power(1-TMath::Power(x[1]/par[2],par[5]), 1/par[5] );     
  //  printf("ap=%f, bp=%f \n",ap,bp);
  
  
  if(x[0]>ap) return -2;
  
  if(ap!=0)
    {
      par2D[0]=ap;  par2D[1]=bp;  par2D[2]=par[3];
      y=KImplant2D::ImplEdge(x[0],par2D);
      //      y=bp*TMath::Power(1-TMath::Power(x[0]/ap,par[3]), 1/par[3] );  
    }
     else
       y=0;

  return y;
}

Double_t KImplant3D::Distance(Double_t *R,Double_t *par, Double_t *RR)
{
  // par[0]= width in X 
  // par[1]= width in Y
  // par[2]= width in Z
  // par[3]= curvature in XY;
  // par[4]= curvature in XZ;
  // par[5]= curvature in YZ;
  
  Double_t par2D[3];
  Double_t SS[3];  
  Double_t XX[3],ap,bp,dist,mindist=1e6;
  Int_t k;
  Double_t z=R[2]>par[2]?par[2]/2:0;
  while(z<=par[2])
    {
       par2D[0]=par[2];  par2D[1]=par[0];  par2D[2]=par[4];
       ap=KImplant2D::ImplEdge(z,par2D);
       par2D[0]=par[2];  par2D[1]=par[1];  par2D[2]=par[5];
       bp=KImplant2D::ImplEdge(z,par2D); 
 
       par2D[0]=ap;  par2D[1]=bp;  par2D[2]=par[3];
       KImplant2D::Distance(R,par2D,XX);      
       XX[2]=z;

       dist=0;
       for(k=0;k<3;k++) dist+=TMath::Power(R[k]-XX[k],2);
       dist=TMath::Sqrt(dist);

       if(dist<mindist) 
	 {
	     mindist=dist; 
	     for(k=0;k<3;k++) SS[k]=XX[k];
	 }
       else break;

        if(z<0.9*par[2]) z+=0.1; else 
	 z+=0.02;
    }

  //if((Float_t)R[0]<=(Float_t)SS[0] && (Float_t)R[2]<=(Float_t)SS[2] && (Float_t)R[1]<=(Float_t)SS[1]) mindist=-mindist;
    XX[0]=R[0];
    XX[1]=R[2];
    XX[2]=KImplant3D::ImplEdge(XX,par);
    if(XX[2]>0 && XX[2]>=R[1]) mindist=-mindist;
  
  if(RR!=NULL) for(k=0;k<3;k++) RR[k]=SS[k];
  return mindist;
}


Double_t KImplant3D::Conc(Double_t *x, Double_t Thresh)
{
  // x[0]  = coordinate x
  // x[1]  = coordinate y
  // x[2]  = coordinate z

  Double_t dist=Distance(x,Dim);
  if(dist<Thresh) return 0; else
  return fConc->Eval(dist);
}

KImplant3D::KImplant3D(Double_t *x, Double_t Sigma, Double_t Nimpl)
{
  // x[0]  = coordinate x
  // x[1]  = coordinate y 

  for(Int_t i=0;i<6;i++) Dim[i]=x[i];
  fConc=new TF1("fConc","TMath::Erfc((x-[0])/[1])*[2]+[3]",-20,20);
  fConc->SetParameter(3,0);
  fConc->SetParameter(2,Nimpl);
  fConc->SetParameter(1,Sigma);
  fConc->SetParameter(0,0);

}
