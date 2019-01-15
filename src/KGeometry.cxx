
#include "KGeometry.h"
#include "TMath.h"
#include "math.h"
#include <iostream> // cout

ClassImp(KGeometry)

//------------------------------------------------------------------------------
KGeometry::KGeometry()
{
  EG = NULL;
  DM = NULL;
  nx = 1;
  ny = 1;
  nz = 1;
}

//------------------------------------------------------------------------------
KGeometry::~KGeometry()
{
  if( EG != NULL ) delete EG;
  if( DM != NULL ) delete DM;
}

//------------------------------------------------------------------------------
void KGeometry::GetGrid( TH3I *x, Short_t which )
{
  switch(which)
    {
    case 0: 
      if( EG != NULL ) delete EG;
      EG = new TH3I();
      x->Copy(*EG);
      nx = EG->GetNbinsX();
      ny = EG->GetNbinsY();
      nz = EG->GetNbinsZ();
      break;
    case 1:
      if( DM != NULL ) delete DM;
      DM = new TH3I();
      x->Copy(*DM); 
      break;
    }

  if( DM != NULL ) if( DM->GetNbinsX() != nx ) printf( "Warning: dimenssions mismatch - X !\n");
  if( DM != NULL ) if( DM->GetNbinsY() != ny ) printf( "Warning: dimenssions mismatch - Y !\n");
  if( DM != NULL ) if( DM->GetNbinsZ() != nz ) printf( "Warning: dimenssions mismatch - Z !\n");

} // GetGrid

//------------------------------------------------------------------------------
void KGeometry::Reset( Int_t Which, Int_t What )
{
  for( int k = 1; k <= nz; ++k )
    for( int j = 1; j <= ny; ++j )
      for( int i = 1; i <= nx; ++i )
	if( Which )
	  if( EG != NULL )
	    EG->SetBinContent( i, j, k, What );
	  else 
	    if( DM != NULL )
	      DM->SetBinContent( i, j, k, What );
}

//------------------------------------------------------------------------------
Int_t KGeometry::SetBoundaryConditions()
{
  if( EG == NULL ) {
    printf( "Please set the geometry first ! \n" );
    return -1;
  }

  std::cout << "KGeometry::SetBoundaryConditions:" << std::endl;

  nx = EG->GetNbinsX();
  ny = EG->GetNbinsY();
  nz = EG->GetNbinsZ();

  //Bit 1 = 1 -> GND - 0 V bias
  //Bit 2 = 2 -> Voltage (usual bias voltage)
  //Bit 15-32 = 32768 -> Additional Voltages 

  // Bits determining boundary conditions:
  // bit 2 = 4  -> down val
  // bit 3 = 8  -> up val 
  // bit 4 = 16 -> left val
  // bit 5 = 32 -> right val
  // bit 6 = 64  -> down der
  // bit 7 = 128  -> up der 
  // bit 8 = 256 -> left der
  // bit 9 = 512 -> right der
  // the 3D section is separated:
  // bit 10 = 1024 -> out val
  // bit 11 = 2048 -> in val
  // bit 12 = 4096 -> out der
  // bit 13 = 8192 -> in der
  // bit 14 = 16384 -> read out node

  for( int k = 1; k <= nz; ++k )

    for( int j = 1; j <= ny; ++j )

      for( int i = 1; i <= nx; ++i ) {

	int cval = EG->GetBinContent( i, j, k );

	if( !( cval & 1 || cval & 2 || cval >= 32768 ) ) {

	  int nval = 0;

	  if( i+1 <= nx ) {
	    int val = EG->GetBinContent( i+1, j, k );
	    if( val & 1 || val & 2 || val >= 32768 ) nval |= 32;
	  }
	  else
	    nval |= 512;

	  if( i-1 > 0 ) {
	    int val = EG->GetBinContent( i-1, j, k );
	    if( val & 1 || val & 2 || val >= 32768 ) nval |= 16;
	  }
	  else
	    nval |= 256;

	  if( j+1 <= ny ) {
	    int val = EG->GetBinContent( i, j+1, k );
	    if( val & 1 || val & 2 || val >= 32768 ) nval |= 8;
	  }
	  else
	    nval |= 128;

	  if( j-1 > 0 ) {
	    int val = EG->GetBinContent( i, j-1, k ); 
	    if( val & 1 || val & 2 || val >= 32768 )  nval |= 4; 
	  }
	  else
	    nval |= 64;

	  if( k+1 <= nz ) {
	    int val = EG->GetBinContent( i, j, k+1 );
	    if( val & 1 || val & 2 || val >= 32768 ) nval |= 2048;
	  }
	  else
	    if( nz != 1 ) nval |= 8192;

	  if( k > 1 ) {
	    int val = EG->GetBinContent( i, j, k-1 );
	    if( val & 1 || val & 2 || val >= 32768 ) { // k is next to external electrode
	      nval |= 1024; // GK
	      //nval |= 4096; // DP, gives E spike
	      //nval |= 2048; // DP, no convergence
	      //nval |= 8192; // DP, no convergence
	      //printf( "%5d %3d %3d : %5d\n", i, j, k, nval );
	    }
	  }
	  else
	    if( nz != 1 ) nval |= 1024; // DP, no change
	  //if( nz != 1 ) nval |= 4096; // GK
	  //if( nz != 1 ) nval |= 8192; // DP for n in p, no change

	  EG->SetBinContent( i, j, k, nval );
	  // if(k==nz) printf("%d %d %d :: %d \n",i,j,k,nval);

	} // cval

      } // loops i,j,k

  return 0;

} // SetBoundaryConditions

//------------------------------------------------------------------------------
// map NR vector x into 3D field
TH3D * KGeometry::MapToGeometry( Double_t *x, Double_t Scale )
{
  TH3D *fhis = new TH3D();
  EG->Copy(*fhis);
  fhis->Reset();

  // Map the array of values: E, U, W ... to the geometry.
  //   Double_t xb=EG->GetXaxis()->GetBinUpEdge(nx);
  //   Double_t yb=EG->GetYaxis()->GetBinUpEdge(ny);
  //   Double_t zb=EG->GetZaxis()->GetBinUpEdge(nz);
  //   Double_t bx=EG->GetXaxis()->GetBinLowEdge(1);
  //   Double_t by=EG->GetYaxis()->GetBinLowEdge(1);
  //   Double_t bz=EG->GetZaxis()->GetBinLowEdge(1);

  //   TH3D *fhis=new TH3D("Pot3d","Pot3d",nx,bx,xb,ny,by,yb,nz,bz,zb);

  for( int k = 1; k <= nz; ++k )
    for( int j = 1; j <= ny; ++j )
      for( int i = 1; i <= nx; ++i ) {
	int n = (k-1)*nx*ny + (j-1)*nx + i; 
	fhis->SetBinContent( i, j, k, x[n] * Scale );
      }
  return fhis;
}

//------------------------------------------------------------------------------
Double_t KGeometry::GetStepSize(Int_t dir, Int_t i)
{
  Double_t Lo,Hi;
  Double_t ret;
  switch(dir)
    {
    case 0:
      Hi=fabs(EG->GetXaxis()->GetBinCenter(i+1)-EG->GetXaxis()->GetBinCenter(i));
      Lo=fabs(EG->GetXaxis()->GetBinCenter(i)-EG->GetXaxis()->GetBinCenter(i-1));
      break;
    case 1:
      Hi=fabs(EG->GetYaxis()->GetBinCenter(i+1)-EG->GetYaxis()->GetBinCenter(i));
      Lo=fabs(EG->GetYaxis()->GetBinCenter(i)-EG->GetYaxis()->GetBinCenter(i-1));
      break;
    case 2:
      Hi=fabs(EG->GetZaxis()->GetBinCenter(i+1)-EG->GetZaxis()->GetBinCenter(i));
      Lo=fabs(EG->GetZaxis()->GetBinCenter(i)-EG->GetZaxis()->GetBinCenter(i-1));
      break;
    default:
      Hi=-1; Lo=-1;
      break;
    }
  ret=0.5*Hi+0.5*Lo;
  return ret;
}

//------------------------------------------------------------------------------
Double_t KGeometry::GetStepSize(Int_t dir, Float_t x)
{
  Int_t bin;
  switch(dir)
    {
    case 0:
      bin=EG->GetXaxis()->FindBin(x);
      break;
    case 1:
      bin=EG->GetYaxis()->FindBin(x);
      break;
    case 2:
      bin=EG->GetZaxis()->FindBin(x);
      break;
    default:
      bin=0;
      break;
    }
  return GetStepSize(dir, bin);
}

//------------------------------------------------------------------------------
Float_t KGeometry::GetLowEdge(Int_t dir)
{
  Float_t ret=0;
  switch(dir)
    {
    case 0: ret=EG->GetXaxis()->GetBinLowEdge(1); break;
    case 1: ret=EG->GetYaxis()->GetBinLowEdge(1); break;
    case 2: ret=EG->GetZaxis()->GetBinLowEdge(1); break;
    default: printf("Index %i out of scope!\n",dir); ret=0; break;
    }
  return ret;
}

//------------------------------------------------------------------------------
Float_t KGeometry::GetUpEdge(Int_t dir)
{
  Float_t ret=0;
  switch(dir)
    {
    case 0: ret=EG->GetXaxis()->GetBinUpEdge(nx); break;
    case 1: ret=EG->GetYaxis()->GetBinUpEdge(ny); break;
    case 2: ret=EG->GetZaxis()->GetBinUpEdge(nz); break;
    default: printf("Index %i out of scope!\n", dir); ret=0; break;
    }
  return ret;
}

//------------------------------------------------------------------------------
void KGeometry::ElLine(Float_t *r0,Float_t *r1, Float_t *W, Int_t Wei, Int_t Mat )
{
  Int_t i,j,q,k,Bx,By,Bz;
  Float_t t;
  Float_t p[3],r[3];

  for( i = 0; i < 3; ++i )
    p[i] = r1[i] - r0[i];

  for( t = 0; t <= 1; t += 0.01 ) {

    for( i = 0; i < 3; ++i )
      r[i] = p[i]*t + r0[i];

    Bz = EG->GetZaxis()->FindBin(r[2]);
    By = EG->GetYaxis()->FindBin(r[1]);      
    Bx = EG->GetXaxis()->FindBin(r[0]);

    for( k = -(Int_t)(W[2]/GetStepSize(2,Bz))+Bz; k <= (Int_t)(W[2]/GetStepSize(2,Bz))+Bz; ++k )
      if( k <= nz && k >= 1 ) {
	for( j = -(Int_t)(W[1]/GetStepSize(1,By))+By; j <= (Int_t)(W[1]/GetStepSize(1,By))+By; ++j )
	  if( j <= ny && j >= 1 )
	    for( q = -(Int_t)(W[0]/GetStepSize(0,Bx))+Bx;q<=(Int_t)(W[0]/GetStepSize(0,Bx))+Bx;++q )
	      if( j <= nx && j >= 1 ) {
		if( EG != NULL ) EG->SetBinContent( q, j, k ,Wei );
		if( DM != NULL ) DM->SetBinContent( q, j, k ,Mat );
	      }
      } // k, j, q

  } // t

}

//------------------------------------------------------------------------------
void KGeometry::ElRectangle(Float_t *Pos, Float_t *Size, Int_t Wei, Int_t Mat)
{
  // Sets Up 
  Int_t i,j,k,q;
  Int_t xpl,ypl,zpl,xpr,ypr,zpr;
  //Set up left edge of the box
  xpl=EG->GetXaxis()->FindBin(Pos[0]-Size[0]);
  ypl=EG->GetYaxis()->FindBin(Pos[1]-Size[1]);
  zpl=EG->GetZaxis()->FindBin(Pos[2]-Size[2]);
  //Set up rigth edge of the box	   
  xpr=EG->GetXaxis()->FindBin(Pos[0]+Size[0]);
  ypr=EG->GetYaxis()->FindBin(Pos[1]+Size[1]);
  zpr=EG->GetZaxis()->FindBin(Pos[2]+Size[2]);
  //	   printf("(%d %d),(%d %d),(%d %d)\n",xpl,xpr,ypl,ypr,zpl,zpr);
  // Fill the geometry histogram
  for( k = zpl; k <= zpr; ++k )	   
    for( j = ypl; j <= ypr; ++j )
      for( i = xpl; i <= xpr; ++i ) {
	//	  if(x0+i>1 && x0+i<=nx && y0+i>1 && y0+i<=ny)
	//	  printf("%d %d %d %d\n",i,j,k,Wei);
	if( EG != NULL) EG->SetBinContent( i, j, k, Wei );
	if( DM != NULL) DM->SetBinContent( i, j, k, Mat );
      }
}

//------------------------------------------------------------------------------
void KGeometry::ElCylinder(Float_t *Pos,Float_t R, Float_t L,Int_t O, Int_t Wei, Int_t Mat)
{
  // Cylindrical electrode 
  // Float_t *Pos;  - postion of the cone center 
  Float_t Dist,D,x,y,z,Bu,Bd;
  Int_t i,j,k,q;

  for( k = 1; k <= nz; ++k )

    for( j = 1; j <= ny; ++j )

      for( i = 1; i <= nx; ++i ) {

	x = EG->GetXaxis()->GetBinCenter(i);
	y = EG->GetYaxis()->GetBinCenter(j);
	z = EG->GetZaxis()->GetBinCenter(k);

	switch(O)
	  {
	  case 3:   D = TMath::Sqrt( TMath::Power(x-Pos[0],2) + TMath::Power(y-Pos[1],2) ) - R; break;
	  case 2:   D = TMath::Sqrt( TMath::Power(x-Pos[0],2) + TMath::Power(z-Pos[2],2) ) - R; break;
	  case 1:   D = TMath::Sqrt( TMath::Power(y-Pos[1],2) + TMath::Power(z-Pos[2],2) ) - R; break;
	  }

	if( D <= 0 ) {
	  switch(O)
	    {
	    case 3:
	      Dist = EG->GetZaxis()->GetBinCenter(k); 
	      Bu=Pos[2]+L;
	      Bd=Pos[2]-L;
	      break;
	    case 2:
	      Dist = EG->GetYaxis()->GetBinCenter(j); 
	      Bu=Pos[1]+L;
	      Bd=Pos[1]-L;
	      break;
	    case 1:
	      Dist = EG->GetXaxis()->GetBinCenter(i); 
	      Bu=Pos[0]+L;
	      Bd=Pos[0]-L;
	      break;
	    }

	  if( Dist <= Bu && Dist >= Bd) {
	    if( EG != NULL ) EG->SetBinContent( i, j, k, Wei ); 
	    if( DM != NULL ) DM->SetBinContent( i, j, k, Mat );
	  }

	} // D < 0

      } // loops i,j,k
}

//------------------------------------------------------------------------------
TH3D * KGeometry::GetGeom()
{
  TH3D * dhis = new TH3D();

  EG->Copy(*dhis);
  dhis->SetName( Form( "%s_%i", EG->GetName(), GetNhs() ) );
  dhis->SetTitle( "electrode geometry" );

  for( int k = 1; k <= nz; ++k )
    for( int j = 1; j <= ny; ++j )
      for( int i = 1; i <= nx; ++i ) {
	int bin = dhis->GetBinContent( i, j, k );
	int col = 0;
	if( bin > 32768 ) col = 1; 
	if( bin == 1 || bin == 2 ) col = 2;
	if( bin == 16385) col = 3;	 
	dhis->SetBinContent( i, j, k, col );
      }
  //  dhis->Draw("glbox");
  return dhis;
}

//------------------------------------------------------------------------------
TH2D * KHisProject( void * hisIn, Int_t axis, Int_t Bin1 )
{
  //Projects any quantity maped to geometry in different views

  TH3D * his = (TH3D *) hisIn;

  Int_t Nx = his->GetNbinsX();
  Int_t Ny = his->GetNbinsY();
  Int_t Nz = his->GetNbinsZ();

  Double_t *Xbins=new Double_t [Nx+1];
  Double_t *Ybins=new Double_t [Ny+1];;
  Double_t *Zbins=new Double_t [Nz+1];;

  //  printf("%d %d %d\n", Nx,Ny,Nz);
  for( int i = 0; i <= Nx; ++i ) {
    Xbins[i] = his->GetXaxis()->GetBinLowEdge(i+1);
    //   printf("x:: %d %f\n",i,Xbins[i]);
  }
  for( int i = 0; i <= Ny; ++i ) {
    Ybins[i] = his->GetYaxis()->GetBinLowEdge(i+1);
    //   printf("x:: %d %f\n",i,Ybins[i]);
  }

  for( int i = 0; i <= Nz; ++i ) {
    Zbins[i] = his->GetZaxis()->GetBinLowEdge(i+1);
    //    printf("x:: %d %f\n",i,Zbins[i]);
  }

  TH2D *his2D;
  switch(axis)
    {
    case 1:
      his2D = new TH2D( Form( "YZ_%i", GetNhs() ),
			Form( "YZ plane at x %f", his->GetXaxis()->GetBinCenter(Bin1) ),
			Ny, Ybins, Nz, Zbins );
      for( int i = 1; i <= Ny; ++i )
	for( int j = 1; j <= Nz; ++j )
	  his2D->SetBinContent( i, j, his->GetBinContent(Bin1,i,j) );
      his2D->GetXaxis()->SetTitle("y [#mum]");
      his2D->GetYaxis()->SetTitle("z [#mum]");
      break;
    case 2:
      his2D = new TH2D( Form( "XZ_%i", GetNhs() ),
			Form( "XZ plane at y %f", his->GetYaxis()->GetBinCenter(Bin1) ),
			Nx, Xbins, Nz, Zbins );
      for( int i = 1; i <= Nx; ++i )
	for( int j = 1;j <= Nz; ++j )
	  his2D->SetBinContent( i, j, his->GetBinContent(i,Bin1,j) );
      his2D->GetXaxis()->SetTitle("x [#mum]");
      his2D->GetYaxis()->SetTitle("z [#mum]");
      break;
    case 3:
      his2D = new TH2D( Form( "XY_%i", GetNhs() ),
			Form( "XY plane at z %f", his->GetZaxis()->GetBinCenter(Bin1) ),
			Nx, Xbins, Ny, Ybins );
      for( int i = 1; i <= Nx; ++i )
	for( int j = 1; j <= Ny; ++j )
	  his2D->SetBinContent( i, j, his->GetBinContent(i,j,Bin1) );
      his2D->GetXaxis()->SetTitle("x [#mum]");
      his2D->GetYaxis()->SetTitle("y [#mum]");
      break;
    }

  return his2D;

} // KHisProject

//------------------------------------------------------------------------------
int GetNhs()
{
  static int nhs{0};
  ++nhs;
  //std::cout << "  nhs " << nhs << std::endl;
  return nhs;
};
