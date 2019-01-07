
#include "KPixel.h"

ClassImp(KPixel)

//------------------------------------------------------------------------------
KPixel::KPixel( Int_t n, Float_t x, Float_t y, Float_t z )
{
  nPix = n;
  PSx = new Float_t [nPix];
  PSy = new Float_t [nPix];
  PSd = new Float_t [nPix];
  PSWx = new Float_t [nPix];
  PSWy = new Float_t [nPix];
  PSW = new Short_t [nPix];

  for( Int_t i = 0; i < nPix; i++ ) {
    PSx[i] = 0;
    PSy[i] = 0;
    PSd[i] = 0;
    PSWx[i] = 0;
    PSWy[i] = 0;
    PSW[i] = 0;
  }

  CellX = x; // [um]
  CellY = y;
  CellZ = z;

}

//------------------------------------------------------------------------------
KPixel::~KPixel()
{
  delete PSx; 
  delete PSy; 
  delete PSd;
  delete PSWx;
  delete PSWy;
  delete PSW;
}

//------------------------------------------------------------------------------
void KPixel::SetUpVolume( Float_t St1, Float_t St2, Float_t St3 )
{
  nx = (int)(CellX/St1);
  ny = (int)(CellY/St2);
  nz = (int)(CellZ/St3);

  std::cout << "KPixel::SetUpVolume nx " << nx
	    << ", ny " << ny
	    << ", nz " << nz
	    << std::endl;

  EG = new TH3I( "EG", "EG_0", nx, 0, CellX, ny, 0, CellY, nz, 0, CellZ );
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");
  EG->GetZaxis()->SetTitle("z [#mum]");

  DM = new TH3I( "DM", "DM_0", nx, 0, CellX, ny, 0, CellY, nz, 0, CellZ );
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");
  DM->GetZaxis()->SetTitle("z [#mum]");

}

//------------------------------------------------------------------------------
void KPixel::SetUpPixel( Int_t i, Float_t posX, Float_t posY, Float_t WX, Float_t WY,
			 Float_t Depth, Short_t Weight )
{
  if( i < 0 || i >= nPix ) {
    std::cout << "KPixel::SetUpPixel wrong index " << i
	      << " (must be 0 < " << nPix << ")"
	      << std::endl;
    return;
  }
  PSx[i] = posX;
  PSy[i] = posY;
  PSWx[i] = WX;
  PSWy[i] = WY;
  PSd[i] = Depth;
  PSW[i] = Weight;
}

//------------------------------------------------------------------------------
void KPixel::SetPixelW( Int_t i, Short_t Weight )
{
  if( i < 0 || i >= nPix ) {
    std::cout << "KPixel::SetPixelW wrong index " << i
	      << " (must be 0 < " << nPix << ")"
	      << std::endl;
    return;
  }
  PSW[i] = Weight;
}

//------------------------------------------------------------------------------
void KPixel::SetUpElectrodes( Int_t Material )
{
  for( Int_t k = 1;k <= nz; ++k )
    for( Int_t j = 1;j <= ny; ++j )
      for( Int_t i = 1;i <= nx; ++i ) {
	if( k == 1 ) { // surface
	  EG->SetBinContent( i, j, k, 2 ); // bit 2 = bias voltage here
	  DM->SetBinContent( i, j, k, 100 ); // KMaterial: 100 = Alu, overwritten below
	}
	else { // DP: { E-field does not converge
	  EG->SetBinContent( i, j, k, 0 );
	  //DM->SetBinContent( i, j, k, Material );
	} // DP
	DM->SetBinContent( i, j, k, Material ); // intentional?
      } // loops k,j,i

  for( Int_t p = 0; p < nPix; ++p ) {

    // find left and right pixel edges:

    Int_t xpl = EG->GetXaxis()->FindBin( PSx[p] - PSWx[p] );
    Int_t ypl = EG->GetYaxis()->FindBin( PSy[p] - PSWy[p] );
    if( xpl < 1 ) xpl = 1;
    if( ypl < 1 ) ypl = 1;
    Int_t zpl = nz; // top bin = pixel electrode

    Int_t xpr = EG->GetXaxis()->FindBin( PSWx[p] + PSx[p] );
    Int_t ypr = EG->GetYaxis()->FindBin( PSWy[p] + PSy[p] );
    if( xpl > nx ) xpr = nx;
    if( ypl > ny ) ypr = ny;
    Int_t zpr = EG->GetZaxis()->FindBin( CellZ - PSd[p] ); // implant thickness

    // printf("Hole %d:: Bins:  X(%d %d) Y(%d %d) Z(%d %d)\n",p,xpl,xpr,ypl,ypr,zpl,zpr);

    for( Int_t k = zpl; k >= zpr; --k ) // top layer(s)
      for( Int_t j = ypl; j <= ypr; ++j )
	for( Int_t i = xpl; i <= xpr; ++i ) {
	  //	  if(x0+i>1 && x0+i<=nx && y0+i>1 && y0+i<=ny)
	  EG->SetBinContent( i, j, k, PSW[p] ); // weight = electrode flag
	}

  } // px

  enp[0] = CellX/2;  exp[0] = enp[0];
  enp[1] = 1;        exp[1] = CellY;
  enp[2] = CellZ/2;  exp[2] = CellZ/2;

}
