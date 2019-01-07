
#include "KField.h"
#include "KGeometry.h"
#include <iostream> // cout

//------------------------------------------------------------------------------
Float_t KInterpolate2D( TH3F *his, Float_t x, Float_t y, Int_t dir, Int_t bin )
{
  Int_t EX1,EX2,EY1,EY2;
  Float_t t,u,ret;
  TAxis *ax1,*ax2;
  Float_t v11,v21,v22,v12;

  ret=0;

  switch(dir)
    {
    case 3: ax1=his->GetXaxis(); ax2=his->GetYaxis(); break;
    case 2: ax1=his->GetXaxis(); ax2=his->GetZaxis(); break;
    case 1: ax1=his->GetYaxis(); ax2=his->GetZaxis(); break;
    }

  EX1=ax1->FindBin(x);
  if(ax1->GetBinCenter(EX1)<=x)
    EX2 = EX1+1;
  else {
    EX2=EX1;
    EX1--;
  }
  EY1 = ax2->FindBin(y);
  if(ax2->GetBinCenter(EY1)<=y)
    EY2=EY1+1;
  else {
    EY2=EY1;
    EY1--;
  }

  if( EY2 > ax2->GetNbins() ) {
    u=0;
    EY2=ax2->GetNbins();
  }
  else
    if( EY1 < 1 ) {
      u=0;
      EY1=1;
    }
    else
      u=(y-ax2->GetBinCenter(EY1))/(ax2->GetBinCenter(EY2)-ax2->GetBinCenter(EY1));

  if(EX2>ax1->GetNbins()) {
    t=0;
    EX2=ax1->GetNbins();
  }
  else
    if(EX1<1) {
      t=0;
      EX1=1;
    }
    else
      t=(x-ax1->GetBinCenter(EX1))/(ax1->GetBinCenter(EX2)-ax1->GetBinCenter(EX1));

  //     printf("Points are:: %d %d %d %d [dir=%d, bin=%d] (t=%f u=%f)\n",EX1,EX2,EY1,EY2, dir, bin, t,u);

  switch(dir)
    {
    case 3:
      v11=his->GetBinContent(EX1,EY1,bin); v21=his->GetBinContent(EX2,EY1,bin);
      v22=his->GetBinContent(EX2,EY2,bin); v12=his->GetBinContent(EX1,EY2,bin);
      break;
    case 2:
      v11=his->GetBinContent(EX1,bin,EY1); v21=his->GetBinContent(EX2,bin,EY1);
      v22=his->GetBinContent(EX2,bin,EY2); v12=his->GetBinContent(EX1,bin,EY2);
      break;
    case 1:
      v11=his->GetBinContent(bin,EX1,EY1); v21=his->GetBinContent(bin,EX2,EY1);
      v22=his->GetBinContent(bin,EX2,EY2); v12=his->GetBinContent(bin,EX1,EY2);
      break;
    }

  ret=(1-t)*(1-u)*v11;
  ret+=t*(1-u)*v21;
  ret+=t*u*v22;
  ret+=(1-t)*u*v12;
  return ret;
}

//------------------------------------------------------------------------------
// E = grad U
Float_t KField::GetFieldPoint( Float_t *X, Float_t *Y )
{
  Float_t a,b,k12,k23;
  k12 = ( Y[0] - Y[1] ) / ( X[0] - X[1] );
  k23 = ( Y[1] - Y[2] ) / ( X[1] - X[2] );
  a = (k23-k12) / ( X[2] - X[0] );
  b = k12 - a * ( X[0] + X[1] );
  return 2*a*X[1] + b;
}

//------------------------------------------------------------------------------
Int_t KField::CalField() // E = -grd U
{
  Float_t X[3], Y[3], EE;

  if( U == NULL ) {
    printf( "Cannot calculate field - no potential array!\n" );
    return -1;
  };

  Int_t Nx = U->GetNbinsX();
  Int_t Ny = U->GetNbinsY();
  Int_t Nz = U->GetNbinsZ();

  if( Nz == 1 ) {
    printf("2D field!\n");
    dim = 2;
  }
  else
    dim = 3;

  Ex = new TH3F();
  U->Copy(*Ex);
  Ex->Reset();
  //Ex->SetTitle( Form( "%s_%i", U->GetName(), GetNhs() ) );

  Ey = new TH3F();
  U->Copy(*Ey);
  Ey->Reset();
  //Ey->SetTitle( Form( "%s_%i", U->GetName(), GetNhs() ) );

  Ez = new TH3F();
  U->Copy(*Ez);
  Ez->Reset();
  //Ez->SetTitle( Form( "%s_%i", U->GetName(), GetNhs() ) );

  E = new TH3F();
  U->Copy(*E);
  E->Reset();
  //E->SetTitle( Form( "%s_%i", U->GetName(), GetNhs() ) );

  for( int k = 1; k <= Nz; ++k )
    for( int j = 1; j <= Ny; ++j )
      for( int i = 1; i <= Nx; ++i ) {

	if( i==1 || i==Nx )
	  Ex->SetBinContent( i, j, k, 0 );
	else {
	  for( int q = 0; q <= 2; ++q ) {
	    X[q] = U->GetXaxis()->GetBinCenter(i+q-1); // left, center, right
	    Y[q] = U->GetBinContent(i+q-1,j,k);
	  }
	  Ex->SetBinContent( i, j, k, GetFieldPoint(X,Y) );
	}

	if( j==1 || j==Ny )
	  Ey->SetBinContent( i, j, k, 0 ); // periodic boundary conditions?
	else {
	  for( int q = 0; q <= 2; ++q ) {
	    X[q] = U->GetYaxis()->GetBinCenter(j+q-1);
	    Y[q] = U->GetBinContent(i,j+q-1,k);
	  }
	  Ey->SetBinContent( i, j, k, GetFieldPoint(X,Y) );
	}

	if( k==1 || k==Nz )
	  Ez->SetBinContent( i, j, k, 0 );
	else {
	  for( int q = 0; q <= 2; ++q ) {
	    X[q] = U->GetZaxis()->GetBinCenter(k+q-1);
	    Y[q] = U->GetBinContent(i,j,k+q-1);
	  }
	  Ez->SetBinContent( i, j, k, GetFieldPoint(X,Y) );
	}
	
	EE = TMath::Sqrt( TMath::Power( Ex->GetBinContent(i,j,k), 2 ) +
			  TMath::Power( Ey->GetBinContent(i,j,k), 2 ) +
			  TMath::Power( Ez->GetBinContent(i,j,k), 2 ) );
	E->SetBinContent(i,j,k,EE);

      }
  /*
  std::cout << "  made Ex named " << Ex->GetName()
	    << " titled " << Ex->GetTitle()
	    << std::endl;
  std::cout << "  made Ey named " << Ey->GetName()
	    << " titled " << Ey->GetTitle()
	    << std::endl;
  std::cout << "  made Ez named " << Ez->GetName()
	    << " titled " << Ez->GetTitle()
	    << std::endl;
  std::cout << "  made E named " << E->GetName()
	    << " titled " << E->GetTitle()
	    << std::endl;
  */
  return 0;

} // CalField

//------------------------------------------------------------------------------
void  KField::CalFieldXYZ( Float_t x, Float_t y, Float_t z, Float_t *E )
{
  if( dim == 2 ) {
    E[1] = KInterpolate2D( Ex, x, y );
    E[2] = KInterpolate2D( Ey, x, y );
    E[3] = 0;
  }
  else {

    Int_t Nx = Ez->GetXaxis()->GetNbins(); 
    Int_t Ny = Ez->GetYaxis()->GetNbins(); 
    Int_t Nz = Ez->GetZaxis()->GetNbins(); 

    if( Ez->GetZaxis()->FindBin(z) <= 1 ) {
      E[1] = KInterpolate2D( Ex, x, y, 3, 1 );
      E[2] = KInterpolate2D( Ey, x, y, 3, 1 );
      E[3] = KInterpolate2D( Ez, x, y, 3, 1 );
    }

    else if( Ez->GetZaxis()->FindBin(z) >= Nz ) {
      E[1] = KInterpolate2D( Ex, x, y, 3, Nz );
      E[2] = KInterpolate2D( Ey, x, y, 3, Nz );
      E[3] = KInterpolate2D( Ez, x, y, 3, Nz );
    }

    else if( Ez->GetYaxis()->FindBin(y) <= 1 ) {
      E[1] = KInterpolate2D( Ex, x, z, 2, 1 );
      E[2] = KInterpolate2D( Ey, x, z, 2, 1 );
      E[3] = KInterpolate2D( Ez, x, z, 2, 1 );
    }

    else if( Ez->GetYaxis()->FindBin(y) >= Ny ) {
      E[1] = KInterpolate2D( Ex, x, z, 2, Ny );
      E[2] = KInterpolate2D( Ey, x, z, 2, Ny );
      E[3] = KInterpolate2D( Ez, x, z, 2, Ny );
    }

    else if( Ez->GetXaxis()->FindBin(x) <= 1 ) {
      E[1] = KInterpolate2D( Ex, y, z, 1, 1 );
      E[2] = KInterpolate2D( Ey, y, z, 1, 1 );
      E[3] = KInterpolate2D( Ez, y, z, 1, 1 );
    }

    else if( Ez->GetXaxis()->FindBin(x) >= Nx ) {
      E[1] = KInterpolate2D( Ex, y, z, 1, Nx );
      E[2] = KInterpolate2D( Ey, y, z, 1, Nx );
      E[3] = KInterpolate2D( Ez, y, z, 1, Nx );
    }

    else {
      E[1] = Ex->Interpolate( x, y, z );
      E[2] = Ey->Interpolate( x, y, z );
      E[3] = Ez->Interpolate( x, y, z );
    }

  } // 3 D

  E[0] = TMath::Sqrt( E[1]*E[1] + E[2]*E[2] + E[3]*E[3] ); // magnitude

} // CalFieldXYZ

//------------------------------------------------------------------------------
TVector3 * KField::CalFieldXYZ( Float_t x, Float_t y, Float_t z )
{
  Float_t E[4];
  CalFieldXYZ( x, y, z, E ); 
  TVector3 * vec = new TVector3( E[1], E[2], E[3] );
  return vec;
  delete vec;
}

//------------------------------------------------------------------------------
Float_t KField::CalPotXYZ(Float_t x, Float_t y, Float_t z)
{
  Float_t ret = 0;

  if( dim == 2 )
    ret = KInterpolate2D(U,x,y);

  else {

    Int_t Nz = U->GetZaxis()->GetNbins(); 
    Int_t Ny = U->GetYaxis()->GetNbins(); 
    Int_t Nx = U->GetXaxis()->GetNbins(); 

    if( z >= U->GetZaxis()->GetBinCenter(Nz) )
      ret = KInterpolate2D( U, x, y, 3, Nz );

    else if( y >= U->GetYaxis()->GetBinCenter(Ny) )
      ret = KInterpolate2D( U, x, z, 2, Ny );

    else if( x >= U->GetXaxis()->GetBinCenter(Nx) )
      ret = KInterpolate2D( U, y, z, 1, Nx );

    else if( z <= U->GetZaxis()->GetBinCenter(1)  )
      ret = KInterpolate2D( U, x, y, 3, 1 );

    else if( y <= U->GetYaxis()->GetBinCenter(1)  )
      ret = KInterpolate2D( U, x, z, 2, 1 );

    else if( x <= U->GetXaxis()->GetBinCenter(1)  )
      ret = KInterpolate2D( U, y, z, 1, 1 );

    else 
      ret = U->Interpolate(x,y,z);

  } // 3D

  return ret;

} // CalPotXYZ

//------------------------------------------------------------------------------
Float_t KField::DriftVelocity( Float_t E, Float_t Charg, Float_t T, Double_t Neff, Int_t which )
{
  E *= 1e4; // [V/cm]
  return( (float) ( Mobility( E, T, Charg, Neff, which ) * E ) );
}

//------------------------------------------------------------------------------
Float_t KField::DriftVelocity( Float_t cx, Float_t cy, Float_t cz, Float_t Charg,
			       Float_t T, Double_t Neff, Int_t which )
{
  Float_t E[4];
  CalFieldXYZ( cx, cy, cz, E);
  E[0] *= 1e4; // [V/cm]

  return((float) ( Mobility( E[0], T, Charg, Neff, which ) * (Double_t) E[0]) ); 
}

//------------------------------------------------------------------------------
Double_t KField::Mobility( Float_t cx, Float_t cy, Float_t cz, Float_t T,
			   Float_t Charg, Double_t Neff, Int_t which )
{
  Float_t E[4]; // DP: 3 -> 4
  CalFieldXYZ( cx, cy, cz, E );
  E[0] *= 1e4; 
  return( Mobility( E[0], T, Charg, Neff, which ) );
}

//------------------------------------------------------------------------------
Double_t KField::Mobility( Float_t E, Float_t T, Float_t Charg, Double_t Neff, Int_t which )
{
  Double_t lfm=0,hfm=0;
  Double_t vsatn,vsatp,vsat;
  Double_t betap,betan;
  Double_t alpha;

  switch( which ) // 1 is default: Canali
    {
    case 0:
      alpha=0.72*TMath::Power(T/300,0.065);
      if( Charg > 0 ) {
	Double_t ulp=460*TMath::Power(T/300,-2.18);
	Double_t uminp=45*TMath::Power(T/300,-0.45);
	Double_t Crefp=2.23e17*TMath::Power(T/300,3.2);
	betap=1;
	vsatp=9.05e6*TMath::Sqrt(TMath::TanH(312/T));
	lfm=uminp+(ulp-uminp)/(1+TMath::Power(Neff/Crefp,alpha));
	hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatp,betap),1/betap));
      }
      else {
	Double_t uln=1430*TMath::Power(T/300,-2);
	Double_t uminn=80*TMath::Power(T/300,-0.45);
	Double_t Crefn=1.12e17*TMath::Power(T/300,3.2);
	betan=2;
	vsatn=1.45e7*TMath::Sqrt(TMath::TanH(155/T));
	lfm=uminn+(uln-uminn)/(1+TMath::Power(Neff/Crefn,alpha));
	hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatn,betan),1/betan));
      }
      break;
    case 1:
      //printf("%e ",par[0]);
      if( Charg > 0 ) {
	lfm=8.54e5*TMath::Power(T,-1.075)*TMath::Exp(1-T/124.);
	vsatp=1.445e7*TMath::Exp(-T/435.9);
	betap=2.49*TMath::Exp(-T/270.3);
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
      }
      else {
	lfm=2.712e8*TMath::Power(T,-2.133);
	vsatn=1.586e7*TMath::Exp(-T/723.6);
	betan=-8.262e-8*TMath::Power(T,3)+6.817e-5*TMath::Power(T,2)-1.847e-2*T+2.429;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
      }
      break;
    case 2:   // WF2
      if( Charg > 0 ) {
	lfm=480;
	vsatp=9.5e6;
	betap=1;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
      }
      else {
	lfm=1350;
	vsatn=1.1e7;
	betan=0.5;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
      }
      break;
    case 3: // Klanner Scharf
      Double_t bb,cc,E0;
      if( Charg > 0 ) {
	E0=2970*TMath::Power(T/300,5.63);
	bb=9.57e-8*TMath::Power(T/300,-0.155);
	cc=-3.24e-13;
	lfm=457*TMath::Power(T/300,-2.80);
	if(E>E0) hfm=1./(1/lfm+bb*(E-E0)+cc*TMath::Power(E-E0,2)); else hfm=lfm;
      }
      else {
	E0=2970*TMath::Power(T/300,5.63);
	lfm=1430*TMath::Power(T/300,-1.99);
	vsatn=1.05e7*TMath::Power(T/300,-3.02);
	if(E>E0) hfm=1./(1/lfm+1/vsatn*(E-E0)); else hfm=lfm;
      }
      break;
    case 4:   // Jacoboni
      if( Charg > 0 ) {
	lfm = 474 * TMath::Power(T/300., -2.619);
	vsatp = 0.940e7  * TMath::Power(T/300., -0.226);
	betap = 1.181 * TMath::Power(T/300., 0.633 ); // <100> orientation
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,betap),1/betap);
      }
      else {
	lfm = 1440*  TMath::Power(T/300., -2.260);
	vsatn = 1.054e7  *  TMath::Power(T/300., -0.602);
	betan = 0.992 *  TMath::Power(T/300., 0.572); // <100> orientation
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,betan),1/betan);
      }
      break;
    case 9:
      if( Charg > 0 )
	hfm=0;
      else
	hfm=0;
      break;
    case 10: //Diamond parametrization
      if( Charg > 0 ) {
	lfm=2064;
	vsat=14.1e6;
      }
      else {
	lfm=1714;
	vsat=9.6e6;
      };
      hfm = lfm / ( 1+(lfm*E)/vsat );
      break;
    }
  return hfm;

} // Mobility

//------------------------------------------------------------------------------
TH2F * KField::Draw( Char_t *opt, Int_t b1, Int_t b2 )
{
  TH2F * ret = NULL;
  if( !strcmp( "U", opt ) ) ret = KHisProject( U,  b1, b2 );
  if( !strcmp( "E", opt ) ) ret = KHisProject( E,  b1, b2 );
  if( !strcmp( "X", opt ) ) ret = KHisProject( Ex, b1, b2 );
  if( !strcmp( "Y", opt ) ) ret = KHisProject( Ey, b1, b2 );
  if( !strcmp( "Z", opt ) ) ret = KHisProject( Ez, b1, b2 );
  return ret;
}

//------------------------------------------------------------------------------
// multiplication
Float_t KField::M( Int_t dir, Float_t a1, Float_t a2, Float_t a3 )
{
  //http://www.iue.tuwien.ac.at/phd/park/node36.html
  Int_t NumBin;
  Int_t p1 = E->GetXaxis()->FindBin(a1); NumBin = E->GetNbinsX();
  Int_t p2 = E->GetYaxis()->FindBin(a2); NumBin = E->GetNbinsY();
  Int_t p3 = E->GetZaxis()->FindBin(a3); NumBin = E->GetNbinsZ();

  Float_t dx,EF;
  Double_t I1 = 0, I2 = 0, It = 0;

  for(Int_t i = E->GetNbinsY(); i >= p2; i-- ) {
    dx = E->GetYaxis()->GetBinCenter(i)-E->GetYaxis()->GetBinCenter(i-1);
    //      dx*=1e-6;
    EF = E->GetBinContent(p1,i,p3);
    I1+= (alpha(EF)-beta(EF))*dx;
    if(i%10==0 || i>2940) {printf("i=%d, x=%f, EF=%e alpha,beta (%e %e) I1=%e\n",i,E->GetYaxis()->GetBinCenter(i),EF,alpha(EF),beta(EF),I1);}
  }
  printf("I1=%f (dx=%e, EF=%e)\n",I1,dx,EF);

  for( Int_t i = p2; i >= 1; i-- ) {
    It = 0;
    for( Int_t j = E->GetNbinsY();j >= i; j-- ) { 
      dx = E->GetYaxis()->GetBinCenter(j)-E->GetYaxis()->GetBinCenter(j-1);
      EF = E->GetBinContent(p1,j,p3);
      It += ( alpha(EF) - beta(EF) ) * dx;
    }

    dx = E->GetYaxis()->GetBinCenter(i)-E->GetYaxis()->GetBinCenter(i-1);
    EF = E->GetBinContent(p1,i,p3);
    I2 += beta(EF) * TMath::Exp(It)*dx; 
    //if(i%10==0 || i>2940) printf("I2=%e It=%e ::x=%e (dx=%e, EF=%e)\n",I2,It,E->GetYaxis()->GetBinCenter(i),dx,EF);
  }

  printf("I2=%f (dx=%e, EF=%e)\n",I2,dx,EF);
  return (TMath::Exp(I1)/(1-I2));

} // M

//------------------------------------------------------------------------------
KField::~KField()
{
  if( U  != NULL ) delete U; 
  if( Ex != NULL ) delete Ex; 
  if( Ey != NULL ) delete Ey; 
  if( Ez != NULL ) delete Ez;
}
