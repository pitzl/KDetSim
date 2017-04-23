#include "K3D.h"


ClassImp(K3D)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// K3D                                                                  //
//                                                                      //
// Description of the 3D detector.                                      //
//////////////////////////////////////////////////////////////////////////

K3D::~K3D()
{
    delete PosD; 
    delete PosX; 
    delete PosY;
    delete PosR;
    delete PosW;
}

K3D::K3D(Int_t x1, Float_t x, Float_t y, Float_t z)
{
  Col=x1;
  PosD=new Float_t [Col];
  PosX=new Float_t [Col];
  PosY=new Float_t [Col];
  PosR=new Float_t [Col];
  PosW=new Short_t [Col];
  PosM=new Short_t [Col];

  for(Int_t i=0;i<Col;i++) 
    {PosD[i]=0; PosX[i]=0; PosY[i]=0; PosR[i]=0; PosW[i]=0; PosM[i]=0;}

  CellZ=z;
  CellX=x;
  CellY=y;
}

void K3D::SetUpColumn(Int_t n, Float_t posX, Float_t posY, Float_t R, Float_t Depth, Short_t Wei,Short_t Mat)
{
  PosD[n]=Depth;
  PosR[n]=R;
  PosX[n]=posX;
  PosY[n]=posY;
  PosW[n]=Wei;
  PosM[n]=Mat;
}


void K3D::SetUpVolume(Float_t St1, Float_t St2)
{
  nx=(int)(CellX/St1);
  ny=(int)(CellY/St1);
  nz=(int)(CellZ/St2);
  
 EG=new TH3I("EG","EG",nx,0,CellX,ny,0,CellY,nz,0,CellZ);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");
  EG->GetZaxis()->SetTitle("z [#mum]");

  GetGrid(EG,1);
}

void K3D::SetUpElectrodes(Int_t back)
{
  Float_t Pos[3],L=0;
  Int_t i,j,k;
  if(back)
    {
    for(k=1;k<=nz;k++)
     for(j=1;j<=ny;j++)
      for(i=1;i<=nx;i++)
	if(k==1) EG->SetBinContent(i,j,k,2); 	
	      else 
	         EG->SetBinContent(i,j,k,0); 
    }

  for(Int_t i=0;i<Col;i++)
    {
      Pos[0]=PosX[i];
      Pos[1]=PosY[i];
      Pos[2]=PosD[i]>=0?Pos[2]=PosD[i]/2.:(2*CellZ+PosD[i])/2.; 
      L=TMath::Abs(PosD[i]/2.); 
      ElCylinder(Pos,PosR[i],L,3,PosW[i],PosM[i]);
    }

}

void K3D::SetUpMaterial(Int_t mat)
{
 for(int i=1;i<=nx;i++)
   for(int j=1;j<=ny;j++)
     for(int k=1;k<=nz;k++) DM->SetBinContent(i,j,k,mat);
}


