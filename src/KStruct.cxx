
#include "KStruct.h"

//------------------------------------------------------------------------------
void KStruct::Clear()
{
  for( Int_t i = 0; i < MAXPOINT; i++ ) {
    Xtrack[i]=0;
    Ytrack[i]=0;
    Ztrack[i]=0; 
    Charge[i]=0;
    Time[i]=0; 
    Efield[i]=0;
    MulCar[i]=0;
  }
  Ylength=0;
  Xlength=0;
  TTime=0;
  Steps=0;
  DStrip=0;
  for( int ipx = 0; ipx < 99; ++ipx )
    TCharge[ipx] = 0;
}

//------------------------------------------------------------------------------
KStruct::KStruct()
{
  //std::cout << " KStruct" << std::endl;
  Clear();
}

//------------------------------------------------------------------------------
KStruct:: ~KStruct()
{
  //std::cout << "~KStruct" << std::endl;
}

//------------------------------------------------------------------------------
void KStruct::Info()
{
  printf( "\nParticle Charge= %d\n", PCharge );
  printf( "Number of Steps= %d\n", Steps );
  printf( "Drift Strip    = %d\n", DStrip );
  printf( "Pathlength X   = %f\n", Xlength );
  printf( "Pathlength Y   = %f\n", Ylength );
  printf( "Pathlength Z   = %f\n", Zlength );
  printf( "Total Time     = %f\n", TTime*1e9 );
  printf( "Total Charge   = %f\n", TCharge[1] );
}

//------------------------------------------------------------------------------
void KStruct::Draw(Char_t *option)
{
  Float_t *x,*y;
  if(! strcmp(option,"xy") || ! strcmp(option,"yx")) {x=Xtrack; y=Ytrack;}
  if(! strcmp(option,"xt") || ! strcmp(option,"tx")) {x=Time; y=Xtrack;}
  if(! strcmp(option,"xc") || ! strcmp(option,"cx")) {x=Xtrack; y=Charge;}
  if(! strcmp(option,"yt") || ! strcmp(option,"ty")) {x=Time; y=Ytrack;}
  if(! strcmp(option,"yc") || ! strcmp(option,"cy")) {x=Ytrack; y=Charge;}
  if(! strcmp(option,"tc") || ! strcmp(option,"ct")) {x=Time; y=Charge;}
  
  TGraph *gr = new TGraph(Steps-1,&x[3],&y[3]);
  gr->SetTitle( "Vector Plot" );
  gr->Draw("AL"); gr->GetHistogram()->Draw(); gr->Draw("AL*");

}

//------------------------------------------------------------------------------
void KStruct::GetCH( TH1F *histo, Int_t Update, Float_t Mult, Float_t tau )
{
  // Float_t tau;  trapping time [ns]

  TH1F * his = new TH1F( "ch", "Charge Histogram",
			 histo->GetNbinsX(),
			 histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax() );

  Double_t * ch = new Double_t[Steps+1];

  Axis_t * ti = new Axis_t[Steps+1];

  //std::cout << "KStruct::GetCH at " << ch << std::endl;

  for( Int_t i = 1; i < Steps+1; ++i ) {
    ch[i] = Charge[i];
    ti[i] = Time[i];
  }

  his->FillN( Steps, &ti[1], &ch[1] );

  // Trappping is included if tau>0   // added 20.1.2010
  // exponential damping with time

  if( tau > 0 )
    for( Int_t i = 1; i < his->GetNbinsX(); i++ ) 
      his->SetBinContent( i, his->GetBinContent(i) *
			  TMath::Exp( -( his->GetBinCenter(i) - Time[0] ) / tau ) );

  if( Update ) {
    his->Scale(Mult);
    histo->Add(his);
  }
  else
    his->Copy(*histo);

  delete his;
  delete [] ti;
  delete [] ch;

} // GetCH

//------------------------------------------------------------------------------
Float_t KStruct::GetCHMult(TH1F *histo,Int_t Update,Float_t Mult, Float_t tau)
{
  //   Float_t tau;  trapping time [s]

  Double_t dx,dy,dif,dift,sum=1,summ=1.,mulf,traf;
  Int_t i;

  TH1F *his=new TH1F("ch","Charge Histogram",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());

  Double_t *ch=new Double_t[Steps+1];

  Double_t *Multi=new Double_t[Steps+1];

  Axis_t *ti=new Axis_t[Steps+1];

  for(i=1;i<Steps+1;i++) {
    if( PCharge < 0 )
      dif = KAlpha( 0.5*(Efield[i+1]+Efield[i] ), -1, KMaterial::ImpactIonization );
    else
      dif = KAlpha( 0.5*(Efield[i+1]+Efield[i] ),  1, KMaterial::ImpactIonization );

    dift=Time[i+1]-Time[i];
    dx = TMath::Sqrt( TMath::Power((Xtrack[i+1]-Xtrack[i]),2) + TMath::Power((Ytrack[i+1]-Ytrack[i]),2) );

    // multiplication and trapping factors:

    mulf = (1+dif*dx);
    traf = (1-dift/tau);
    MulCar[i] = (mulf-1)*sum*traf; 
    summ *= mulf;
    sum *= mulf;
    if( tau > 0 ) sum *= traf;

    ch[i] = Charge[i]*sum;
    ti[i] = Time[i];

    // printf("%d :: X=%4.1f , Y=%4.1f :: E=%4.2e ::  Time=%4.1e ; Charge=%4.1e ; dif=%5.2e ; MultT=%5.4e Mult=%5.4f hole=%5.3e\n",i,Xtrack[i],Ytrack[i],Efield[i],ti[i],ch[i],dif,sum,summ,MulCar[i]);

  } 

  his->FillN( Steps, &ti[1], &ch[1] );

  if( Update ) {
    his->Scale(Mult);
    histo->Add(his);
  }
  else
    his->Copy(*histo);

  delete his;
  delete [] ti;
  delete [] ch;
  delete [] Multi;
  return summ;

} // GetCHMult

//------------------------------------------------------------------------------
TH1D * KStruct::GetElFieldAlongTheDrift()
{
  Double_t Step = TMath::Sqrt( TMath::Power( Xtrack[1] - Xtrack[0], 2 ) +
			       TMath::Power( Ytrack[1] - Ytrack[0], 2 ) +
			       TMath::Power( Ztrack[1] - Ztrack[0], 2 ) );

  TH1D *his = new TH1D( "E-Track", "E-Track", Steps, 0, Steps*Step );

  for( Int_t i = 0; i < Steps; i++ ) // bin 0 is underfow ?
    his->SetBinContent( i, Efield[i] );

  return his;
}
