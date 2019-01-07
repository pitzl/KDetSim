
#include "KMaterial.h"

Double_t KM( TH1D *his, Float_t Start, Short_t Rev )
{
  // The derivation of the calculation can be seen in textbooks or e.g. here:
  // http://www.iue.tuwien.ac.at/phd/park/node36.html
  // The function calculates Gain factor in given electric field
  // TH1D *his;  - electric field profile
  // Float_t Start; - creation point of the e-h pair
  //                  for multiplication junction Start=0;
  //                  the opposite side of the junction Start=thickness;
  // Short_t Rev; - in the calculation it is assumed that electrons drift in 
  //                the high field region (Rev=0). If the holes drift then 
  //                Rev=1;

  Int_t NumBin=his->GetNbinsX();
  Int_t sb=his->GetXaxis()->FindBin(Start);
      
  // printf("Start Bin=%d, NumBins=%d\n",sb,NumBin);
       
  Double_t dx,EF;
  Double_t I1=0,I2=0,It=0;

  //  Calculate the integral in numerator

  for( Int_t i = NumBin; i > sb; i-- ) {
    dx = his->GetXaxis()->GetBinCenter(i) - his->GetXaxis()->GetBinCenter(i-1);
    EF = 0.5 * ( his->GetBinContent(i) + his->GetBinContent(i-1) );
    I1 += ( KAlpha(EF, -Rev, KMaterial::ImpactIonization ) -
	    KAlpha(EF,  Rev, KMaterial::ImpactIonization ) ) * dx;
  }

  //   Calculate the integral in denumerator

  It = I1;

  for( Int_t i = sb; i > 1; i-- ) {

    dx = his->GetXaxis()->GetBinCenter(i) - his->GetXaxis()->GetBinCenter(i-1);
    EF = 0.5 * ( his->GetBinContent(i) + his->GetBinContent(i-1) );
    It += ( KAlpha( EF,-Rev, KMaterial::ImpactIonization ) -
	    KAlpha( EF, Rev, KMaterial::ImpactIonization ) ) * dx;
    I2 += KAlpha( EF, Rev,KMaterial::ImpactIonization ) * TMath::Exp(It) * dx; 
  }

  printf("I1=%f, I2=%f It=%f (dx=%e, expI1=%e)\n",I1,I2,It, dx,TMath::Exp(I1));
  return (TMath::Exp(I1)/(1-I2));

} // KM


Double_t KAlpha( Double_t E, Short_t Charg, Int_t which )
{
  // Function calculates impact ionization coefficientf
  // for a given E [V/um]. 
  // Short_t Charg;  ---> Charg=1; holes
  //                 ---> Charg=-1; electrons
  // Int_t which;    ---> 0 -> silicon
  //                 ---> 10 -> diamond Trew parametrization
  //                 ---> 11 -> diamond Watanabe parametrization
  //                 ---> 12 -> diamond Hiraiwa parametrization

  Double_t  alp,A,B,a=TMath::Sqrt(10),b=TMath::Sqrt(10);

  switch(which)
    {
    case 0:  // silicon
      if(Charg>0)
	alp=1.3e-3*TMath::Exp(-13.2*(2e7/(E*1e6)-1));
      else
	alp=2.3e-1*TMath::Exp(-6.78*(2e7/(E*1e6)-1));
      break;
    case 10:
      // Trew parametrization
      A=1.935e4;
      B=7.749e2;      
      alp=A*TMath::Exp(-B/E);      
      break;
    case 11:
      // Watanabe parametrization
      if(Charg>0)
	{A=19.3; B=4.41e2;}
      else
	{A=46.2; B=7.59e2;}
      alp=A*TMath::Exp(-B/E);    
      break;
    case 12:
      // Hiraiwa parametrization
      if(Charg>0)
	{A=19.3/a; B=4.41e2*b;}
      else
	{A=46.2/a; B=7.59e2*b;}
      alp=A*TMath::Exp(-B/E);      
      break;
    }
  return alp;
	
}


ClassImp(KMaterial)

Int_t KMaterial::Mat = 1; // silicon
Float_t KMaterial::Temperature=293;
Int_t KMaterial::Mobility = 1; // Canali
Int_t KMaterial::ImpactIonization=0;

Float_t KMaterial::Perm(Int_t Material)
{
  Float_t perm;
  switch(Material)
    {
    case 0: perm=11.7; break; //silicon 
    case 1: perm=11.7; break; //poly silicon 
    case 2: perm=3.9;  break; //silicon oxide 2.648 
    case 10: perm=5.7;  break; //diamond
    case 20: perm=1; break;   //air
    case 100: perm=1; break;   //aluminium
    default: perm=1; break;
    }
  return perm;
}


Int_t KMaterial::MobMod()
{
  Int_t ret;
  switch(Mat)
    {
    case 0: ret=Mobility; break; //if(Mobility==1) ret=1; else ret=0; break; //silicon 
    case 1: ret=8; break;   //poly silicon 
    case 2: ret=9;  break;  //silicon oxide 2.648 
    case 10: ret=10; break; //diamond
    }
  return ret;
}


Float_t KMaterial::dEX( Double_t E, Double_t *x, Double_t *y, Double_t eps )
{
  Int_t k=0;
  Float_t E0=E,p=0;
  Double_t xx=0,dE;
  eps=eps*1e-4;
  KMaterial::Mat=2; // SiO2
  while(E>0.6)
    { 
      dE=dEdx(E)*eps;
      if(dE<0)
	E+=dE; else {dE=-E; E=0;}
      x[k]=xx*1e4;
      xx+=eps;
      y[k]=TMath::Abs(dE/E0);
      //      p+=y[k];
      //      printf("E=%f , dE=%f , x=%f\n",E,y[k],x[k]); 
      k++;
    } 
  y[k]=TMath::Abs(E/E0);
  x[k]=xx*1e4; 
  printf("E=%f , dE=%f , x=%f\n",0.0,y[k],x[k]); 
  //  printf("p=%f",p+y[k]);
  return (Float_t) x[k];
}


Double_t KMaterial::dEdx(Double_t E)
{
  Double_t konst=0.1535; 
  Double_t A=28.086; //atomic mass Si
  Double_t Z=14; //atomic number Si
  Double_t rho=2.33; //density of silicon
  Double_t z=2;  //alpha particles;
  Double_t mass=3727; //alpha particles [MeV]
  Double_t me=0.511; //alpha particles;

  Double_t C0=-4.44,a=0.1492,m=3.25;
  Double_t X0=0.2014,X1=2.87;
  
  Double_t C=0;
  Double_t gamma=E/mass+1;
  Double_t eta=TMath::Sqrt(gamma*gamma-1);
  Double_t beta=TMath::Sqrt(1-1/(gamma*gamma));
  Double_t Wmax=2*me*eta*eta; //me<<mass
  
  Double_t X=TMath::Log10(eta);
  Double_t delta=0;
  Double_t dE;
  Double_t I = (9.76 * Z + 58.8 * TMath::Power(Z,-0.19)) * 1e-6; // Ionisation energy [MeV]

  if(X<X0) delta=0; 
  if(X0<X && X<X1) delta=4.6052*X+C0+a*TMath::Power((X1-X),m); 
  if(X>X1) delta=4.6052*X+C0;

  Double_t logarg=2*me*eta*eta*Wmax/(I*I);
  //  printf("logarg=%e %e %e :::::::: ",logarg,TMath::Log(logarg) , I);
  dE=konst*Z/A*TMath::Power(z/beta,2)*rho*(TMath::Log(logarg)-2*beta*beta-delta-2*C/Z); // Bethe-Bloch

  return(-dE);

}
