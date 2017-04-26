
// Pad detector = diode

#include "Rtypes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "KPad.h"
#include "KMaterial.h"

#define JMAX 40

ClassImp(KPad)


KPad::KPad(Float_t x,Float_t y)
{
  //Constructor of the class KPad:
  //		Float_t x ; width of the diode in um. For the simulation it is not neccessary to put real dimensions
  //		Float_t y ; thickness in um
  Neff=NULL;
  CellX=x;
  CellY=y;
}
 

void KPad::SetUpVolume(Float_t St1)
{
  nx=(Int_t) (CellX/St1);
  ny=(Int_t) (CellY/St1);
  nz=1;
  //Set the boundary condition matrix //
  EG=new TH3I("EG","EG",nx,0,CellX,ny,0,CellY,1,0,1);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");

  DM=new TH3I("DM","DM",nx,0,CellX,ny,0,CellY,1,0,1);
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");
}

void KPad::SetUpElectrodes()
{
 
 for(int i=1;i<=nx;i++){ EG->SetBinContent(i,1,1,2);  EG->SetBinContent(i,ny,1,16385);} 
 KMaterial::Mat=0;
  //Default track
 enp[0]=CellX/2;  exp[0]=enp[0];
 enp[1]=1;        exp[1]=CellY;
 if(Neff!=NULL ) 
   {CalField(0); CalField(1);} 
 else 
   printf("Please define space charge function Neff before field calculation\n");
 
}


//_________________________________________________________________________________
Float_t  KPad::rtbis(float x1, float x2, float xacc)
{
  void nrerror( std::string error_text );
  int j;
  float dx,f,fmid,xmid,rtb;

  f=PoEqSolve(x1);
  fmid=PoEqSolve(x2);
  if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=JMAX;j++) {
    fmid=PoEqSolve(xmid=rtb+(dx *= 0.5));
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}

//___________________________________________________________________________________
void KPad::rk4(float y[], float dydx[], int n, float x, float h, float yout[])
{
  float *vector(long,long);
  void free_vector(float*,long,long);
  int i;
  float xh,hh,h6,*dym,*dyt,*yt;

  dym=vector(1,n);
  dyt=vector(1,n);
  yt=vector(1,n);
  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
  Derivs(xh,yt,dyt);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  Derivs(xh,yt,dym);
  for (i=1;i<=n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  Derivs(x+h,yt,dyt);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  free_vector(yt,1,n);
  free_vector(dyt,1,n);
  free_vector(dym,1,n);
}

//------------------------------------------------------------------------------
KPad::~KPad()
{
}


void KPad::Derivs(float x,float y[],float dydx[])
{
  //Double_t permMat=(material==0)?perm:permDi;
  Double_t permMat=Perm();
  Double_t permV  = perm0*1e-6;;
  dydx[1]=y[2];
  dydx[2]=-Neff->Eval(x)*e_0/(permMat*permV); 
  // if(x>150) dydx[2]*=permMat;
}


Float_t KPad::PoEqSolve(Float_t der)
{
  //TArrayF PhyField(ny+2);
  //TArrayF PhyPot(ny+2);
  //Float_t Step=1;

  Int_t i;
  Float_t h,x=0;

  Float_t y[3],dydx[3],yout[3];

  y[1]=Voltage;
  y[2]=der;
  PhyField[0]=y[2]; PhyPot[0]=y[1];   
  Derivs(x,y,dydx);
  for (i=1;i<=ny;i++) 
    {
      h=GetStepSize(1,i);
      rk4(y,dydx,2,x,h,yout);
      //printf("%f %f %f\n",x+h,yout[1],yout[2]);
      y[1]=yout[1];
      y[2]=yout[2]; 
      PhyField[i]=y[2]; PhyPot[i]=y[1]; 
      x=x+h;
      Derivs(x,y,dydx);
    }
     
  //     printf("y[1]=%f\n",xp1);
  return y[1];
}


void KPad::GetRamoField()
{
  Int_t i,j;
  Double_t *x=new Double_t [nx*ny+1];
  for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
      x[i+(j-1)*nx]=(Double_t)(j-1)/(Double_t)(ny-1); 
  Ramo[0].U=MapToGeometry(x);
  Ramo[0].CalField();
  delete [] x;

}


void KPad::GetRamoField(TH1F *rf)
{

  Int_t i,j;
  Double_t *x=new Double_t [nx*ny+1];

  for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
      x[i+(j-1)*nx]=rf->GetBinContent(j); 

  Ramo[0].U=MapToGeometry(x);
  Ramo[0].CalField();
  delete [] x;
}


void KPad::GetField()
{
  Int_t i,j;
  //Float_t Step=1;
  Float_t aa;
  PhyPot=TArrayF(ny+2);
  PhyField=TArrayF(ny+2);
  Double_t *x=new Double_t [nx*ny+1];
  //TArrayD PhyPot2D(nx*ny+1);
  //TArrayI StripPosition=TArrayI(2); StripPosition[0]=1; StripPosition[1]=nx; 

  aa=rtbis(-100,100,0.000001);
  //if(!Invert && GetDepletionVoltage()>Voltage) aa=rtbis(1,300,0.000001); else aa=rtbis(-10,10,0.000001);

  for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
      {
	x[i+(j-1)*nx]=(Double_t)PhyPot[j-1];
      }
  //Real=EField(PhyPot2D,nx,ny);
  //Real.CalField(Step,1);

  // for(i=0;i<nx*ny+1;i++)     {printf("%f ",PhyPot2D[i]); if(i%nx==0) printf("\n");}

  Real.U=MapToGeometry(x);
  Real.CalField();
  delete [] x;
}


void KPad::GetField(TH1F *rf)
{
  Int_t i,j;
  Double_t *x=new Double_t [nx*ny+1];
  PhyPot=TArrayF(ny+2);
  for(j=1;j<=ny;j++)

    for(i=1;i<=nx;i++)
      x[i+(j-1)*nx]=rf->GetBinContent(j); 
  PhyPot[j-1]=rf->GetBinContent(j); 
  Real.U=MapToGeometry(x);
  Real.CalField();
  delete [] x;
}

void KPad::GetField(TF1 *rf)
{
  Int_t i,j;
  Double_t *x=new Double_t [nx*ny+1];
  PhyPot=TArrayF(ny+2); 
  for(j=1;j<=ny;j++)

    for(i=1;i<=nx;i++)
      x[i+(j-1)*nx]=rf->Eval(EG->GetXaxis()->GetBinCenter(j));
  PhyPot[j-1]=rf->Eval(EG->GetXaxis()->GetBinCenter(j));
  Real.U=MapToGeometry(x);
  Real.CalField();
  delete [] x;
}


// void KPad::GetField(TF1 *potential)
//  {
//  //Set the field as defined in the potential function.
//  Int_t i,j,index;
//  Float_t aa,val;
//  Double_t abspot=0;
//  PhyPot=TArrayF(ny+2);
//  PhyField=TArrayF(ny+2);
//  TArrayD PhyPot2D(nx*ny+1);

//  Real=EField(PhyPot2D,nx,ny);
//   abspot=potential->Integral(0,300);// printf("napetost=%e\n",abspot);
//  for(j=1;j<=ny;j++)
//    {
//    for(i=1;i<=nx;i++)
//      {
//        index=i+(j-1)*nx;
//        Real.SetEfx(index,0);
//        Real.SetEfy(index,-potential->Eval((Double_t)(ny-j)));
//        Real.SetEf (index,-potential->Eval((Double_t)(ny-j)));	  
//      }  
//    PhyField[j]=-potential->Eval((Double_t)(ny-j));
//    if(j==1) PhyPot[j]=abspot; else
//     PhyPot[j]=potential->Integral((Double_t)(ny-j),(Double_t)(ny-j-1))+PhyPot[j-1]; 
//   }
// //Real=EField(PhyPot2D,nx,ny,StripPosition);
// //Real.CalField(Step,1);

// // for(i=0;i<nx*ny+1;i++)     {printf("%f ",PhyPot2D[i]); if(i%nx==0) printf("\n");}
// }



TGraph *KPad::DrawPad(char *option)
{
  //Draws potential "p" or  electric field "f"
  Char_t Opt[10];
  TGraph *gr;
  TArrayF xx=TArrayF(ny+1);
  for(Int_t i=0;i<=ny;i++) xx[i]=(Float_t)i*GetStepSize(1,i);

  if(strchr(option,'p')!=NULL) gr=new TGraph(ny+1,xx.GetArray(),PhyPot.GetArray());
  if(strchr(option,'f')!=NULL) gr=new TGraph(ny+1,xx.GetArray(),PhyField.GetArray());
  if(strchr(option,'s')!=NULL) sprintf(Opt,"CP"); else  sprintf(Opt,"ACP"); 
  //if(! strcmp(option,"p")) gr=new TGraph(ny,xx.GetArray(),PhyPot.GetArray());
  //if(! strcmp(option,"f")) gr=new TGraph(ny,xx.GetArray(),PhyField.GetArray());
  gr->Draw(Opt);
  gr->GetHistogram()->Draw();
  gr->Draw(Opt);
  return gr;
}

Double_t laser(Double_t *x, Double_t *par)
{
  Double_t dvg=0;
  dvg=par[0]*TMath::Exp(- TMath::Power((*x-par[1]),2)/par[2])+par[3]+par[4]*(*x)+par[5]*TMath::Power(*x,2);
  return dvg;
}

// void KPad::LaserV(float *enp, float *exip, int div,Float_t lambda)
// {
//   // The simulation of the drift for the red laser illumination! 
//   // A penetration profile  is devided into Int_ div buckets. Each bucket is drifted in the field. The
//   // induced currents for each carrier is calculated as the sum  all buckets. 

// TH1F trn[4],trp[4];
// Float_t sp[3];;
// DStruct seg;
// Double_t MaxPudu=2.7,pudu=0.;
// Int_t PuduStep=27;
// int i,ii=0,j,k=0,z,cutoff=27,shiftn,shiftp;
// //void printsegdrift(struct segdrift *);
// Double_t *cc= new Double_t[cutoff]; 
// Double_t *tt= new Double_t[cutoff]; 
// Double_t cd;

// float xe,xs;//= new Double_t[cutoff]; 
//  TH1F *histop  = new TH1F("chl+","charge+",200,0,STEPH);
//  TH1F *histon  = new TH1F("chl-","charge-",200,0,STEPH);
//  TH1F *outhisp,*outhisn;

//  Double_t par[6]={-0.038553,1.07634,1.14664,1.54445e-2,-6.08484e-3,6.08802e-4};
//  TF1 *las=new TF1("laser",laser,0,5,6);
//  las->SetParameters(par);

// while(pudu<MaxPudu) 
//   {
//    cc[k]=las->Integral(pudu,pudu+MaxPudu/PuduStep)/las->Integral(MaxPudu/PuduStep,MaxPudu);
//    tt[k]=pudu*1e-9;
//    //  printf("%e %e\n",cc[k],tt[k]);
//    pudu+=MaxPudu/PuduStep; k++; 
//   } 

// sum->Reset(); pos->Reset(); neg->Reset();
//   for(i=0;i<div;i++) 
//   {
//      for(j=0;j<3;j++) sp[j]=((exip[j]-enp[j])/div)*i+enp[j]+(exip[j]-enp[j])/(2*div);
//      xs=((exip[1]-enp[1])/div)*i; 
//      xe=xs+(exip[1]-enp[1])/div;
//      cd=(TMath::Exp(-fabs(xs)/lambda)-TMath::Exp(-fabs(xe)/lambda));
 
//      //      printf("%f %f %f %f %f\n",sp[0],sp[1],fabs(xs),fabs(xe),cd); 

//      for(ii=0;ii<average;ii++)
//        {
//      Drift(sp[0],sp[1],1,&seg,MobMod);
//      seg.GetCH(histop,1); 
//      Drift(sp[0],sp[1],-1,&seg,MobMod);
//      seg.GetCH(histon,1); 
//        }
//      histop->Scale(1/(Float_t)average); histon->Scale(1/(Float_t)average); 
//       // pad2->cd(); histon->Draw(); c1->Update();   printf("%d %e\n",i,histon->Integral());

//      if(trapping) {
//        ht->Trapping(histop,trp);  outhisp=&trp[3];
//        et->Trapping(histon,trn); outhisn=&trn[3]; 
//        //pad2->cd(); histon->Draw(); c1->Update();  
//                   } 
//      else {outhisp=histop; outhisn=histon;}

//      shiftp=(Int_t) (MaxPudu*1e-9/outhisp->GetBinWidth(1))+1;
//      shiftn=(Int_t) (MaxPudu*1e-9/outhisn->GetBinWidth(1))+1;
//      //printf("%d %d\n",shiftp,shiftn);
//      for(k=0;k<PuduStep;k++)
//      {
//      for(j=1;j<=outhisn->GetNbinsX()-shiftn;j++) {z=(Int_t)(tt[k]/histon->GetBinWidth(1)); neg->AddBinContent(j+z,outhisn->GetBinContent(j)*cd*fabs(cc[k]));}
//      for(j=1;j<=outhisp->GetNbinsX()-shiftp;j++) {z=(Int_t)(tt[k]/histop->GetBinWidth(1)); pos->AddBinContent(j+z,outhisp->GetBinContent(j)*cd*fabs(cc[k]));}
//      }
//      histon->Reset(); histop->Reset();
//   }
     
// sum->Add(neg);
// sum->Add(pos);   
// delete histon;
// delete histop;
// delete cc;
// delete tt;
// delete las;
// //shapper(taush,cht);
// }









// void KPad::CalM(DStruct *seg, Double_t *data)
//     {
//       Int_t i,j,numreg=0;
//       Double_t *M=new Double_t [seg->Steps+1];
//       Double_t *DIF=new Double_t [seg->Steps+1];
//       Double_t dx,dy,sum=1.,sum2=0,sum1=0,dif,dif1,dif2;
//       printf("Number of Steps=%d\n",seg->Steps);




    
//       for(i=1; i<seg->Steps;i++)
// 	{
// 	  //	  if(i==seg->Steps) DIF[i]=0;
// 	  dx=TMath::Sqrt(TMath::Power((seg->Xtrack[i+1]-seg->Xtrack[i]),2)+TMath::Power((seg->Ytrack[i+1]-seg->Ytrack[i]),2));
// 	  dif=Real.alpha(0.5*(seg->Efield[i+1]+seg->Efield[i]));
// 	  sum*=(1+dif*dx);
// 	  DIF[i]=sum;
// 	}

//       for(i=1;i<seg->Steps;i++) { printf("Step=%d [%f %f], E=%f , a[i]=%f M(i)=%f\n",i,seg->Xtrack[i],seg->Ytrack[i],seg->Efield[i],Real.alpha(seg->Efield[i]),DIF[i]); }


//       data[0]=DIF[seg->Steps-1]; data[1]=0; data[2]=0; data[3]=0;
//       printf("KKK %f\n",data[0]);
//       for(i=1;i<seg->Steps;i++) 
// 	{
// 	  if(DIF[i]>data[0]*0.2) 
// 	    {
// 	      data[1]+=seg->Xtrack[i]; data[2]+=seg->Ytrack[i]; numreg++;
// 	      data[3]+=seg->Time[i];
// 	       printf("::::: %e %e %e %d\n",data[1],data[2],data[3],numreg); 	      
// 	    }
// 	}
//    data[1]/=numreg; data[2]/=numreg; data[3]/=numreg; 


        
//     }

