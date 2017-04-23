#include "KStrip.h"
#include "KMesh.h"
ClassImp(KStrip)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KStrip                                                               //
//                                                                      //
// Class for description of silicon microstrip detector                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


KStrip::KStrip(Float_t x1, Float_t x2, Float_t x3, Int_t x4, Float_t x5)
{
   //Constructor:K
   // Float_t x1;  Strip Pitch
   // Float_t x2;  Strip Width
   // Int_t x3;    Number of strips
   // Float_t x4;  Detector thickness

  Pitch=x1;
  Width=x2;
  Depth=x3;
  NoStrips=x4;
  CellY=x5;
}




void KStrip::SetUpVolume(Int_t numy,Float_t *Ybins, Int_t numx, Float_t *Xbins)
{
  Float_t Zbins=0;
  nx=numx;
  ny=numy;

  CellX=Pitch*NoStrips;

  EG=new TH3I("EG","EG",nx,Xbins,ny,Ybins,1,&Zbins);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");

  DM=new TH3I("DM","DM",nx,Xbins,ny,Ybins,1,&Zbins);
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");
}

void KStrip::SetUpVolume(Float_t StS,Float_t StE,Int_t Num)
{
KMesh Ym(CellY);
Float_t *Ybins=new Float_t [1000];
ny=Ym.GetBins(Num,StS,StE,Ybins);

Int_t i,j;
CellX=Pitch*NoStrips;
nx=(Int_t)(CellX/StS);

 Float_t *Xbins=new Float_t[nx+1];
Float_t Zbins;
 Zbins=1;
 Ybins[0]=0;
 Xbins[0]=0;

 // for(Int_t j=1;j<=ny;j++) 
 // if(j*St1<=Yb) Ybins[j]=Ybins[j-1]+St1; else Ybins[j]=Ybins[j-1]+St2;
 for(Int_t i=1;i<=nx;i++) Xbins[i]= Xbins[i-1]+StS;

  EG=new TH3I("EG","EG",nx,Xbins,ny,Ybins,1,&Zbins);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");

  DM=new TH3I("DM","DM",nx,Xbins,ny,Ybins,1,&Zbins);
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");

  delete Xbins;
  delete Ybins;


}

void KStrip::SetUpVolume(Float_t St1,Float_t St2)
{
CellX=Pitch*NoStrips;
nx=(Int_t)(CellX/St1);
ny=(Int_t)(CellY/St2);
  EG=new TH3I("EG","EG",nx,0,CellX,ny,0,CellY,1,0,1);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");

  DM=new TH3I("DM","DM",nx,0,CellX,ny,0,CellY,1,0,1);
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");

}


void KStrip::SetUpMaterial(Int_t mat)
{
 for(int i=1;i<=nx;i++)
   for(int j=1;j<=ny;j++) DM->SetBinContent(i,j,1,mat);
}


void KStrip::SetUpElectrodes(Int_t rs)
{
Float_t sc,sr,sl;
Int_t i,j,k;
 // Setup default ramo strip
 if(rs>=0) RamoStrip=rs; else
 if(rs==-1) RamoStrip=NoStrips/2; else RamoStrip=-2;
//Back plane metalization
 for(i=1;i<=nx;i++) EG->SetBinContent(i,1,1,2);
//DC Strips
 for(i=0;i<NoStrips;i++)
   {
     sc=(i+1)*(CellX/NoStrips)-(Pitch/2.);
     sl=sc-(Width/2.);
     sr=sc+(Width/2.);
    
     for(k=EG->GetXaxis()->FindBin(sl);k<=EG->GetXaxis()->FindBin(sr);k++) 
       {
	 for(j=EG->GetYaxis()->FindBin(CellY-Depth);j<=ny;j++)
	 if(i==RamoStrip || RamoStrip==-2) 
	   EG->SetBinContent(k,j,1,16385); 
	   else
	   EG->SetBinContent(k,j,1,1); 
       }
   }

 //Setup default track
  enp[0]=((Float_t)NoStrips)/2*Pitch; enp[1]=1;     enp[2]=0.5;
  exp[0]=((Float_t)NoStrips)/2*Pitch; exp[1]=CellY; exp[2]=0.5;
  
}

