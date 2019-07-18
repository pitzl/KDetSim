
#include "TPolyLine3D.h"
#include "KDetector.h"
#include "TFile.h"

#define ABS(x) x>0?x:-x
#define PREDZNAK(x) x>0?1:-1

// Declaration of macros used to setup the solver of partial differential equations

// #define C1(x,y,z) y3[n]=-2./(x*x)-2./(y*y)-2./(z*z);

// #define L1(x) y2[n]=1./(x*x);
// #define R1(x) y4[n]=1./(x*x);
// #define U1(x) y5[n]=1./(x*x);
// #define D1(x) y6[n]=1./(x*x);
// #define I1(x) y7[n]=1./(x*x);
// #define O1(x) y8[n]=1./(x*x);

// #define L2(x) y2[n]=2./(x*x);
// #define R2(x) y4[n]=2./(x*x);
// #define U2(x) y5[n]=2./(x*x);
// #define D2(x) y6[n]=2./(x*x);
// #define I2(x) y7[n]=2./(x*x);
// #define O2(x) y8[n]=2./(x*x);

#define C1(x,y,z) y3[n]=x+y+z;

#define L1(x) y2[n]=x;
#define R1(x) y4[n]=x;
#define U1(x) y5[n]=x;
#define D1(x) y6[n]=x;
#define I1(x) y7[n]=x;
#define O1(x) y8[n]=x;

#define L2(x) y2[n]=2*x;
#define R2(x) y4[n]=2*x;
#define U2(x) y5[n]=2*x;
#define D2(x) y6[n]=2*x;
#define I2(x) y7[n]=2*x;
#define O2(x) y8[n]=2*x;

#define C0 y3[n]=1.;
#define U0 y5[n]=0.;
#define D0 y6[n]=0.;
#define R0 y4[n]=0.;
#define L0 y2[n]=0.;
#define I0 y7[n]=0.;
#define O0 y8[n]=0.;

#define  PI  3.1415927
#define EPS 1.0e-14
#define STEP_DET_CUR 25e-9

ClassImp(KDetector)

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KDetector                                                            //
//                                                                      //
// The base class for all detectors. It incorporates the calculation of //
// electric and weithing field as well as simulation of drift.          //
//                                                                      //
// Calculation of the drift:                                            //
// BEGIN_LATEX                                                          
//  #frac{d^{2}f(x)}{dx^{2}}=#frac{f(x+h1)-2*f(x+0.5(h1-h2))+f(x-h2)}{(0.5 h1+ 0.5 h2)^{2}}
//  h1=h2=h  
//  #frac{d^{2}f(x)}{dx^{2}}=#frac{f(x+h)-2*f(x)+f(x-h)}{h^{2}}
// END_LATEX                                                            
//////////////////////////////////////////////////////////////////////////

double **a, *b, *y2, *y3, *y4, *y5, *y6, *y7, *y8;
double *dvector(long, long);
void free_dvector(double*,long,long);

//------------------------------------------------------------------------------
bool isOK( double x ) // check against NAN
{
  if( x >= 0 ) return 1;
  if( x <= 0 ) return 1;
  return 0;
}

//------------------------------------------------------------------------------
/**********************************************************
 snrm: Calculates the norm of vector 
 unsigned long n - dimension of the vector
 double sx[]     - components
 itol            - <=3 real norm , otherwise max component.
************************************************************/
double snrm( unsigned long n, double sx[], int itol )
{
  unsigned long i,isamax;
  double ans;

  if( itol <= 3 ) {
    ans = 0.0;
    for( i = 1; i <= n; ++i )
      ans += sx[i]*sx[i];
    return sqrt(ans);
  }
  else {
    isamax = 1;
    for( i = 1; i <= n; ++i )
      if( fabs(sx[i]) > fabs(sx[isamax]) )
	isamax = i;
    return fabs( sx[isamax] );
  }
}

//------------------------------------------------------------------------------
void atimes( unsigned long n, int dim[], double x[], double r[], int itrnsp )
{
  // multiply the vectors with matrices!
  // Used for calculation of the electric and ramo field

  int nx = dim[0];
  int ny = dim[1];
  int nz = dim[2];

  int q = 0;

  for( int k = 1; k <= nz; ++k )
    for( int j = 1; j <= ny; ++j )          /*mnozenje po stolpcu*/
      for( int i = 1; i <= nx; ++i ) {      /*mnozenje po vrstici*/ 
	++q;
	double C,L,D,O,R,U,I;
	C=y3[q]*x[q];
	if(q-1>1)        L=y2[q]*x[q-1]; else L=0;
	if(q-nx>1)       D=y6[q]*x[q-nx]; else D=0;
	if((q-nx*ny)>1)  I=y7[q]*x[q-ny*nx]; else I=0;

	if(q+1<=n)       R=y4[q]*x[q+1]; else R=0;
	if(q+nx<=n)      U=y5[q]*x[q+nx]; else U=0;
	if((q+nx*ny)<=n) O=y8[q]*x[q+ny*nx]; else O=0;
	 
	r[q]=C+L+D+O+R+U+I;
      }

  if( n != q ) printf( "\n Error in matrix solving!" );

  return;
}

//------------------------------------------------------------------------------
void asolve(unsigned long n, double b[], double x[], int itrnsp)
{
  //Solve a very simple system of n equations
  //with the diagonal elements to be the only ones!

  for( unsigned long i = 1; i <= n; ++i )
    x[i] = ( y3[i] != 0.0 ? b[i]/y3[i] : b[i] );
}

//------------------------------------------------------------------------------
// The main function for electric field calculation:
/* (C) Copr. 1986-92 Numerical Recipes Software &c&):)+!. */
/**********************************************************
 linbcg: Solves linear eqution sparse system Ax=b
 unsigned long n - number of equations
 int nx = # of points in x direction
 int ny = # of points in y direction
 double b[] = right side vector
 double x[] = solution of the system
  int itol  = a way to calculte norm
 double tol
 int itmax  = max # of iterations
 int *iter =
 double *err =
************************************************************/
void linbcg( unsigned long n, int dim[], double b[], double x[], int itol, double tol,
	     int itmax, int *iter, double *err)
{
  void asolve( unsigned long n, double b[], double x[], int itrnsp);
  void atimes( unsigned long n, int dim[],double x[], double r[], int itrnsp);
  double snrm( unsigned long n, double sx[], int itol);
  double *dvector(long, long);
  void free_dvector(double *, long, long);
  void nrerror( std::string error_text );
  unsigned long j;
  double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
  double *p,*pp,*r,*rr,*z,*zz;

  bool ldb = 0;

  p  = dvector(1,n);
  pp = dvector(1,n);
  r  = dvector(1,n);
  rr = dvector(1,n);
  z  = dvector(1,n);
  zz = dvector(1,n);

  atimes( n, dim, x, r, 0 );

  for( j = 1; j <= n; ++j ) {
    r[j] = b[j]-r[j];
    rr[j] = r[j];
  }

  atimes( n, dim, r, rr, 0 ); // minimal residual invariant

  znrm = 1.0;
  if( itol == 1 )
    bnrm = snrm( n, b, itol );
  else if( itol == 2 ) {
    asolve( n, b, z, 0 );
    bnrm = snrm( n, z, itol );
  }
  else if( itol == 3 || itol == 4 ) {
    asolve(n,b,z,0);
    bnrm = snrm(n,z,itol);
    asolve(n,r,z,0);
    znrm = snrm(n,z,itol);
  }
  else
    nrerror( "illegal itol in linbcg" );

  if( ldb )
    std::cout << "linbcg bnrm " << bnrm << std::endl;

  asolve( n, r, z, 0 );

  *iter = 0;

  if( isOK( bnrm ) && fabs( bnrm ) > tol ) {

    while( *iter < itmax ) {

      ++(*iter);

      zm1nrm = znrm;
      asolve(n,rr,zz,1);
      for( bknum = 0.0, j = 1; j <= n; ++j )
	bknum += z[j]*rr[j];
      if( *iter == 1 ) {
	for( j = 1; j <= n; ++j ) {
	  p[j] = z[j];
	  pp[j] = zz[j];
	}
      }
      else {
	bk = bknum / bkden;
	for( j = 1; j <= n; ++j ) {
	  p[j] = bk*p[j] + z[j];
	  pp[j] = bk*pp[j] + zz[j];
	}
      }
      bkden = bknum;
      atimes( n, dim, p, z, 0 );
      for( akden = 0.0, j = 1; j <= n; ++j )
	akden += z[j]*pp[j];
      ak = bknum / akden;
      atimes( n, dim, pp, zz, 1 );
      for( j = 1; j <= n; ++j ) {
	x[j] += ak*p[j];
	r[j] -= ak*z[j];
	rr[j] -= ak*zz[j];
      }

      asolve( n, r, z, 0 );

      if( itol == 1 || itol == 2) {
	znrm=1.0;
	*err = snrm( n, r, itol ) / bnrm;
      }
      else if( itol == 3 || itol == 4) {
	znrm=snrm(n,z,itol);
	if( fabs(zm1nrm-znrm) > EPS*znrm) {
	  dxnrm=fabs(ak)*snrm(n,p,itol);
	  *err=znrm/fabs(zm1nrm-znrm)*dxnrm;
	}
	else {
	  *err=znrm/bnrm;
	  continue;
	}
	xnrm = snrm( n, x, itol );
	if( *err <= 0.5*xnrm)
	  *err /= xnrm;
	else {
	  *err=znrm/bnrm;
	  continue;
	}
      }

      if( ldb )
	std::cout << "  " << *iter << "  " << *err << std::endl << std::flush;

      if( *iter == itmax )
	std::cout
	  << "iterations limit " << *iter
	  << ", error " << *err
	  << std::endl;

      if( *err <= tol ) break;

    } // iter

  } // OK
  else
    std::cout << "ERROR in linbcg: bnrm " << bnrm << std::endl;

  free_dvector(p,1,n);
  free_dvector(pp,1,n);
  free_dvector(r,1,n);
  free_dvector(rr,1,n);
  free_dvector(z,1,n);
  free_dvector(zz,1,n);

} // linbcg

//------------------------------------------------------------------------------
KDetector::KDetector()
{
  // Author: Gregor Kramberger
  // Default constructor for KDetector class
  // Default = no space charge -> default space charge is function

  NeffF = new TF3( "Profile", "x[0]*x[1]*x[2]*0+[0]", 0, 1000, 0, 1000, 0, 1000 );
  NeffF->SetParameter(0,0); // no doping
  NeffH = NULL; // no Neff histogram

  for( Int_t i = 0; i < 3; ++i )
    B[i] = 0; // magnetic field

  // setting up default random generator for diffusion

  ran = new TRandom3(4357); // DP

  // Calculation parameters

  CalErr = 1e-6; // tolerance for linbcg
  MaxIter = 2000; // GK

  // histograms for storing the drift
  pos = NULL;
  neg = NULL;
  sum = NULL;
  SetDriftHisto( 25.0, 500 ); // [ns]
  //SetDriftHisto( 10, 500 ); // [ns] DP for edge9 at 400V, tauh=2/F
  qnode[0] = 0;

  // setting up general variables

  Landau = 1;
  taue = -1;   // no electron trapping 
  tauh = -1;   // no hole trapping 
  MTresh = -1; // no multiplication 
  BDTresh = -1;

  // drift:
  Deps = 1e-5; // precision of tracking
  // MobMod=1;  // Mobility parametrization
  average = 1; // average over waveforms
  diff = 0;    // diffusion
  SStep = 1;   // Step size of simulation [um]
  Temperature=263; // temperature
  BreakDown=0; //  no breakdown
  Debug=0;     //  bo printing of debug information
  Voltage2=0;

} // constructor

//------------------------------------------------------------------------------
double KDetector::V( int val, int dowhat )
{
  double voltage;
  int k = 0;
  if( dowhat == 0 ) { // E-field
    if( val&1 ) voltage = 0;  
    if( val&2 ) voltage = Voltage;
    if( val & 32768)
      //if bit 15 is on - many voltages
      voltage = Voltages[val>>16];
  }
  else { // Ramo
    // numerical calculation converges faster if 10'000 is used instead of 1
    // therefore the potential is scaled after calculation to 1
    if( val & 16384 ) // readout electrode
      voltage = 10000; // 1E4
    else
      voltage = 0;
  }
  return voltage;

} // V

//------------------------------------------------------------------------------
Double_t KDetector::kappa(int i,int j, int k,  int dowhat )
{
  //Sets the effective space charge values for given point in the mesh

  Double_t x,y,z,ret;
 
  //  if(NeffF!=NULL && NeffH!=NULL) printf( "Warning:: Histogram values will be taken for Neff!\n" );

  //Position in space
  x = EG->GetXaxis()->GetBinCenter(i);
  y = EG->GetYaxis()->GetBinCenter(j);
  z = EG->GetZaxis()->GetBinCenter(k);

  if( DM != NULL )
    KMaterial::Mat = DM->GetBinContent(i,j,k);
  else
    KMaterial::Mat = 0;

  if( dowhat == 0 ) { // E-field

    if( NeffF != NULL )  // Neff=v enotah [um-3]
      //            ret=(NeffF->Eval(x,y,z)*1e6*e_0)/(KMaterial::Perm()*perm0); /*printf( "i=%d,j=%d,y=%e\n",i,j,y);*/
      ret = ( NeffF->Eval(x,y,z) *1e6 * e_0 ) / (perm0); /*printf( "i=%d,j=%d,y=%e\n",i,j,y);*/

    if( NeffH != NULL )
      //   ret=(NeffH->GetBinContent(i,j,k)*1e6*e_0)/(KMaterial::Perm()*perm0);
      ret = ( NeffH->GetBinContent(i,j,k) * 1e6 * e_0 ) / (perm0);

  }

  //if( dowhat==0) if(j>nc) y=(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0); else y=-(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0);

  else // Ramo
    ret = 0;

  return ret;

}

//------------------------------------------------------------------------------
void KDetector::Declaration( Int_t dowhat )
{
  for( int k = 1; k <= nz; ++k )

    for( int j = 1; j <= ny; ++j )

      for( int i = 1; i <= nx; ++i ) {

	long n = (k-1)*nx*ny + (j-1)*nx + i; // make index of the matrix element

	Int_t ii,jj,kk;
	if( j-1 < 1 ) jj = 1; else jj = j-1;
	if( i-1 < 1 ) ii = 1; else ii = i-1; 
	if( k-1 < 1 ) kk = 1; else kk = k-1; 

	/////////// DEFINE STEPS IN X //////////////////////////////////////

	Double_t Rd = fabs( EG->GetXaxis()->GetBinCenter(i+1) - EG->GetXaxis()->GetBinCenter(i) );
	Double_t Ld = fabs( EG->GetXaxis()->GetBinCenter(i) - EG->GetXaxis()->GetBinCenter(i-1) );
	if( i+1 > nx ) Rd=Ld;
	if( i-1 <  1 ) Ld=Rd;

	////////// DEFINE PERMITIVITY IN X - normal surface ////////////////////////////

	Double_t PRd = Perm( DM->GetBinContent(i,j,k) ) + Perm( DM->GetBinContent(i,jj,k) );

	if( nz != 1 ) {
	  PRd += Perm( DM->GetBinContent(i,j,kk) ) + Perm( DM->GetBinContent(i,jj,kk) );
	  PRd /= 4;
	}
	else
	  PRd /= 2;

	Double_t PLd = Perm( DM->GetBinContent(ii,j,k) ) + Perm( DM->GetBinContent(ii,jj,k) );
	if( nz != 1 ) {
	  PLd += Perm( DM->GetBinContent(ii,j,kk) ) + Perm( DM->GetBinContent(ii,jj,kk) );
	  PLd /= 4;
	}
	else
	  PLd /= 2;

	/////////// DEFINE STEPS IN Y //////////////////////////////////////

	Double_t Ud = fabs( EG->GetYaxis()->GetBinCenter(j+1) -
			    EG->GetYaxis()->GetBinCenter(j) );
	Double_t Dd = fabs( EG->GetYaxis()->GetBinCenter(j) -
			    EG->GetYaxis()->GetBinCenter(j-1) );
	if( j+1 > ny ) Ud = Dd;
	if( j-1 <  1 ) Dd = Ud;

	////////// DEFINE PERMITIVITY IN Y ////////////////////////////

	Double_t PUd = Perm( DM->GetBinContent(i,j,k) ) + Perm( DM->GetBinContent(ii,j,k) );
	if( nz != 1 ) {
	  PUd += Perm( DM->GetBinContent(i,j,kk) ) + Perm( DM->GetBinContent(ii,j,kk) );
	  PUd /= 4;
	}
	else
	  PUd /= 2;

	Double_t PDd = Perm( DM->GetBinContent(i,jj,k) ) + Perm( DM->GetBinContent(ii,jj,k) );
	if( nz != 1 ) {
	  PDd += Perm( DM->GetBinContent(i,jj,kk) ) + Perm( DM->GetBinContent(ii,jj,kk) );
	  PDd /= 4;
	}
	else
	  PDd /= 2;

	/////////// DEFINE STEPS IN Z //////////////////////////////////////

	Double_t Od = fabs( EG->GetZaxis()->GetBinCenter(k+1) -
			    EG->GetZaxis()->GetBinCenter(k) );
	Double_t Id = fabs( EG->GetZaxis()->GetBinCenter(k) -
			    EG->GetZaxis()->GetBinCenter(k-1) );
	if( k+1 > nz ) Od=Id;
	if( k-1 <  1 ) Id=Od;

	//////////DEFINE PERMITIVITY IN Z ////////////////////////////

	Double_t POd = 0;
	Double_t PId = 0;
	if( nz != 1 ) {
	  POd = Perm( DM->GetBinContent(i,jj,k) ) + Perm( DM->GetBinContent(i,j,k) ) +
	    Perm( DM->GetBinContent(ii,j,k) ) + Perm( DM->GetBinContent(ii,jj,k) );

	  PId = Perm( DM->GetBinContent(i,jj,kk) ) + Perm( DM->GetBinContent(i,j,kk) ) +
	    Perm( DM->GetBinContent(ii,j,kk) ) + Perm( DM->GetBinContent(ii,jj,kk) );

	  POd /= 4;
	  PId /= 4;
	}

	if( dowhat > 0 ) { // Ramo
	  PRd=1;
	  PLd=1;
	  PUd=1;
	  PDd=1;
	  POd=1;
	  PId=1;
	}
	 
	Double_t Xr = PRd / (0.5*Rd*(Rd+Ld));
	Double_t Xl = PLd / (0.5*Ld*(Rd+Ld));
	Double_t Xc = -(Xr+Xl);
	Double_t Yr = PUd / (0.5*Ud*(Ud+Dd));
	Double_t Yl = PDd / (0.5*Dd*(Ud+Dd));
	Double_t Yc = -(Yr+Yl);

	Double_t Zr = 0;
	Double_t Zl = 0;
	Double_t Zc = 0;
	if( nz != 1 ) {
	  Zr = POd/(0.5*Od*(Od+Id));
	  Zl = PId/(0.5*Id*(Od+Id));
	  Zc = -(Zr+Zl);
	}

	b[n] = 0;

	int val = EG->GetBinContent(i,j,k);

	if( nz == 1 ) {
	  C1(Xc,Yc,0)
	    I0
	    O0
	    }
	else {
	  C1(Xc,Yc,Zc)
	    I1(Zl)
	    O1(Zr)
	    }

	R1(Xr) // #define above
	  U1(Yr)
	  L1(Xl)
	  D1(Yl)
	  ;

	if(val&4)    { D0 b[n] -= V( EG->GetBinContent(i,j-1,k), dowhat ) * Yl; }
	if(val&8)    { U0 b[n] -= V( EG->GetBinContent(i,j+1,k), dowhat ) * Yr; }
	if(val&16)   { L0 b[n] -= V( EG->GetBinContent(i-1,j,k), dowhat ) * Xl; }
	if(val&32)   { R0 b[n] -= V( EG->GetBinContent(i+1,j,k), dowhat ) * Xr; }
	if(val&1024) { I0 b[n] -= V( EG->GetBinContent(i,j,k-1), dowhat ) * Zl; }
	if(val&2048) { O0 b[n] -= V( EG->GetBinContent(i,j,k+1), dowhat ) * Zr; }

	if(val&64)   { U2(Yr) D0 if(val&8)    { U0 b[n] -= V( EG->GetBinContent(i,j+1,k), dowhat ) * Yr; } }
	if(val&128)  { D2(Yl) U0 if(val&4)    { D0 b[n] -= V( EG->GetBinContent(i,j-1,k), dowhat ) * Yl; } }
	if(val&256)  { R2(Xr) L0 if(val&32)   { R0 b[n] -= V( EG->GetBinContent(i+1,j,k), dowhat ) * Xr; } }
	if(val&512)  { L2(Xl) R0 if(val&16)   { L0 b[n] -= V( EG->GetBinContent(i-1,j,k), dowhat ) * Xl; } }
	if(val&4096) { O2(Zr) I0 if(val&2048) { O0 b[n] -= V( EG->GetBinContent(i,j,k+1), dowhat ) * Zr; } }
	if(val&8192) { I2(Zl) O0 if(val&1024) { I0 b[n] -= V( EG->GetBinContent(i,j,k-1), dowhat ) * Zl; } }

	b[n] -= kappa( i, j, k, dowhat );

	if(val&1 || val&2 || val >= 32768 ) {
	  U0
	    D0
	    L0
	    R0
	    C0
	    O0
	    I0
	    b[n] = V( val, dowhat );
	}

	//if(j<=2 && i<=2) printf( "stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f :: %d\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n],Mat);
	//      if(k==nz && (j==2 || j==ny-1)) printf( "stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n]);

      } // i,j,k

} // Declaration

//------------------------------------------------------------------------------
void KDetector::CalField( Int_t what )
{
  //booking memory:
  int num = nx*ny*nz;
  b  = dvector( 1, num );
  y6 = dvector( 1, num );
  y2 = dvector( 1, num ); 
  y3 = dvector( 1, num );
  y4 = dvector( 1, num );
  y5 = dvector( 1, num );
  y7 = dvector( 1, num );
  y8 = dvector( 1, num );
  Double_t *x = dvector( 1, num );    

  // Setting up the boundary conditions
  std::cout << "KDetector::CalField setting up matrix " << what << " ... " << std::endl << std::flush;
  Declaration( what );

  std::cout << "KDetector::CalField solving matrix with " << num << " cells... " << std::endl << std::flush;

  // matrix solving

  for( int i = 1; i <=num; ++i )
    x[i] = 1.0;

  int dim[3];
  dim[0] = nx;
  dim[1] = ny;
  dim[2] = nz;
  double err;
  int iter;

  linbcg( num, dim, b, x, 1, CalErr, MaxIter, &iter, &err );

  std::cout << "linbcg done after " << iter << " iterations" << std::endl << std::flush;

  if( iter >= MaxIter ) return;

  // Calculating the field
  if( what == 0 ) {
    Real.U = MapToGeometry(x);
    Real.U->SetName( "U" );
    Real.U->SetTitle( "E_pot" );

    Real.CalField(); // E = -grd U
    Real.Ex->SetName( "E_x" );
    Real.Ey->SetName( "E_y" );
    Real.Ez->SetName( "E_z" );
    Real.E->SetName( "E3" );
    Real.Ex->SetTitle( "E_{x}" );
    Real.Ey->SetTitle( "E_{y}" );
    Real.Ez->SetTitle( "E_{z}" );
    Real.E->SetTitle( "|E|" );
  }
  else {
    // Scale back to 1 from 10'000 when storing the potential
    Ramo[what-1].U = MapToGeometry( x, 1e-4 );
    Ramo[what-1].U->SetName( Form( "Ramo_%i", what-1 ) );
    Ramo[what-1].U->SetTitle( Form( "Ramo_%i", what-1 ) );
    //Ramo[what-1].CalField(); // not needed
  }
  // Freeing the memory
  free_dvector(x, 1,num);   free_dvector(b, 1,num);   free_dvector(y2, 1,num); 
  free_dvector(y3, 1,num);  free_dvector(y4, 1,num);  free_dvector(y5, 1,num); 
  free_dvector(y6, 1,num);  free_dvector(y7, 1,num);  free_dvector(y8, 1,num);

} // CalField

//------------------------------------------------------------------------------
void KDetector::Drift( Double_t sx, Double_t sy, Double_t sz, Float_t qsign,
		       KStruct *seg, Double_t t0 )
{
  // Drift simulation for a point charge (Float_t qsign;)
  // starting from ( sx,sy, sz)
  // KStruct *seg stores the drift paths, drift times and induced charges

  // Inclusion of Magnetic field 28.8.2001 - revised 15.10.2012

  TVector3 BB(B);     // Create a magnetic field vector
  Float_t muhe=1650;  // parametrization of the magnetic field 
  Float_t muhh=310;   // is based on simply this mobility parameters

  // Start time in the  absolute domain (used for delayed charge generation in multiplication 
  Double_t t = t0;    // drift time [s]

  //std::cout << "Drift seg at " << seg << std::endl; // always same mem loc

  seg->Clear();
  seg->PCharge = (Int_t) qsign;

  // start drift:

  Double_t cx = sx; // set current coordinates
  Double_t cy = sy;
  Double_t cz = sz;

  Int_t st = 0;

  seg->Xtrack[st] = cx; // put the first point in the KStruct 
  seg->Ytrack[st] = cy;
  seg->Ztrack[st] = cz;
  seg->Time[st] = t*1E9; // [ns]
  seg->Charge[st] = 0;

  Float_t EE[4];
  Real.CalFieldXYZ( cx, cy, cz, EE );

  seg->Efield[st] = EE[0]; // magnitude

  Float_t pathlen = 0;
  Double_t sumc[99];
  for( int ipx = 0; ipx < 99; ++ipx )
    sumc[ipx] = 0;

  Int_t is_hit = 0;

  bool ldb = 0;
  if( ldb )
    std::cout << "Drift from " << sx
	      << " " << sy
	      << " " << sz
	      << std::endl << std::flush;

  do {

    ++st;

    TVector3 FF; // Lorentz drift
    TVector3 EV( EE[1], EE[2], EE[3] );
    if( qsign > 0 )
      FF = EV + muhh * EV.Cross(BB);
    else
      FF = EV - muhe * EV.Cross(BB);

    //printf( "Field : %f %f %f (%f %f %f)  ---- ",FF[0],FF[1],FF[2],(*EE)[0],(*EE)[1],(*EE)[2]);

    Double_t deltacx = 0;
    Double_t deltacy = 0;
    Double_t deltacz = 0;

    if( FF.Mag() != 0 ) {
      deltacy = -SStep*qsign * FF.y() / FF.Mag();
      deltacx = -SStep*qsign * FF.x() / FF.Mag();
      deltacz = -SStep*qsign * FF.z() / FF.Mag();
    }

    if( ldb )
      std::cout << "KDetector::Drift stp " << st
		<< ": FF " << FF.Mag()
		<< std::flush;

    if( DM != NULL )
      KMaterial::Mat = DM->GetBinContent( DM->FindBin( cx, cy, cz ) );
    else
      KMaterial::Mat=0;

    if( ldb )
      std::cout << ", mat " << KMaterial::Mat
		<< std::endl << std::flush;

    // check boundaries:

    Double_t ncx = cx + deltacx;
    if( ncx < GetLowEdge(0) ) {
      ncx = GetLowEdge(0) + Deps;
      is_hit = 3;
    }
    if( ncx > GetUpEdge(0) ) {
      ncx = GetUpEdge(0) - Deps;
      is_hit = 4;
    }

    Double_t ncy = cy + deltacy;
    if( ncy < GetLowEdge(1) ) {
      ncy = GetLowEdge(1) + Deps;
      is_hit = 5;
    }
    if( ncy > GetUpEdge(1) ) {
      ncy = GetUpEdge(1) - Deps;
      is_hit = 6;
    }

    Double_t ncz = cz + deltacz;
    if( ncz < GetLowEdge(2) ) {
      ncz = GetLowEdge(2) + Deps;
      is_hit = 7;
    }
    if( ncz > GetUpEdge(2) ) {
      ncz = GetUpEdge(2) - Deps;
      is_hit = 8;
    }

    if( ldb )
      std::cout << "KDetector::Drift stp " << st
		<< " to " << ncx
		<< ", " << ncy
		<< ", " << ncz
		<< std::flush;

    Float_t EEN[4];
    Real.CalFieldXYZ( ncx, ncy, ncz, EEN );

    if( ldb )
      std::cout << " um, EEN = " << EEN[0]
		<< std::flush;

    Double_t vel =
      Real.DriftVelocity( 0.5 * ( EEN[0] + EE[0] ),
			  qsign, Temperature,
			  TMath::Abs( NeffF->Eval( cx, cy, cz ) ),
			  MobMod() );
    if( ldb )
      std::cout << ", v = " << vel
		<< std::endl << std::flush;

    Double_t difx = 0;
    Double_t dify = 0;
    Double_t difz = 0;

    if( vel < Deps ) {
      deltacx = 0;
      deltacy = 0;
      deltacz = 0;
      is_hit = 10;
    }
    else if( diff ) { // is diffusion ON
      Double_t Stime = SStep*1e-4/vel; // calculate step time [s]
      Double_t sigma =
	TMath::Sqrt( 2 * Kboltz * Real.Mobility( EE[0],
						 Temperature, qsign,
						 TMath::Abs( NeffF->Eval(cx,cy,cz) ),
						 MobMod() ) * Temperature * Stime );
      dify = ran->Gaus(0,sigma)*1e4; // [um]
      difx = ran->Gaus(0,sigma)*1e4;
      if( nz != 1 )
	difz = ran->Gaus(0,sigma)*1e4;

      ncx = cx + deltacx + difx;
      if( ncx < GetLowEdge(0) ) {
	ncx = GetLowEdge(0) + Deps;
	is_hit = 3;
      }
      if( ncx > GetUpEdge(0) ) {
	ncx = GetUpEdge(0) - Deps;
	is_hit = 4;
      }

      ncy = cy + deltacy + dify;
      if( ncy < GetLowEdge(1) ) {
	ncy = GetLowEdge(1) + Deps;
	is_hit = 5;
      }
      if( ncy > GetUpEdge(1) ) {
	ncy = GetUpEdge(1) - Deps;
	is_hit = 6;
      }

      ncz = cz + deltacz + difz;
      if( ncz < GetLowEdge(2) ) {
	ncz = GetLowEdge(2) + Deps;
	is_hit = 7;
      }
      if( ncz > GetUpEdge(2) ) {
	ncz = GetUpEdge(2) - Deps;
	is_hit = 8;
      }

    } // diff

    if( Debug )
      printf( "%d %f : x:%f->%f y:%f->%f z:%f->%f (%f %f %f): Mat=%d :: ",
	      st,qsign,cx,ncx,cy,ncy,cz,ncz,deltacx,deltacy,deltacz, KMaterial::Mat );

    if( ldb )
      std::cout << "KDetector::Drift stp " << st
		<< " at " << ncx
		<< ", " << ncy
		<< ", " << ncz
		<< std::flush;

    // induced signal:

    for( int ipx = 0; ipx < 99; ++ipx ) {
      if( Ramo[ipx].U == NULL ) break;
      Double_t qstp = qsign *
	( Ramo[ipx].CalPotXYZ( ncx, ncy, ncz ) - Ramo[ipx].CalPotXYZ( cx, cy, cz ) );
      seg->Charge[st] += qstp;
      sumc[ipx] += qstp;
    }

    if( vel > Deps )  {
      t = t + SStep*1e-4/vel; //else t+=Stime;
      pathlen += SStep;
    }

    if( ldb )
      std::cout << " um: Ramo signal " << seg->Charge[st]
		<< ", t " << t
		<< ", l " << pathlen
		<< std::endl << std::flush;

    cx = ncx;
    cy = ncy;
    cz = ncz;

    seg->Xtrack[st] = cx;
    seg->Ytrack[st] = cy;
    seg->Ztrack[st] = cz;
    seg->Time[st] = t*1E9; // [ns]

    if( ldb )
      std::cout << "KDetector::Drift stp " << st
		<< " at " << cx
		<< ", " << cy
		<< ", " << cz
		<< std::flush;

    Real.CalFieldXYZ( cx, cy, cz, EE );

    seg->Efield[st] = EE[0];

    if( ldb )
      std::cout << " um, EE = " << EE[0]
		<< std::endl << std::flush;

    // termination of the drift:

    if( t > 25E-9 ) // [s] LHC
      //if( t > 10E-9 ) // [s] DP for edge9 at 400V, tauh=2/F
      is_hit = 9;
    /*
    Float_t WPot = Ramo[0].CalPotXYZ(cx,cy,cz); // should be all nodes
    if( WPot > ( 1 - Deps ) ) is_hit = 1;
    if( TMath::Abs(WPot) < Deps ) is_hit = 2;
    if( ldb )
      std::cout << ", WPot " << WPot
		<< std::endl << std::flush;
    */
    if( st >= MAXPOINT-1 )
      is_hit = 11;

    if( Debug )
      printf( "(t=%e, vel=%e) [Ch=%f ChInt=%f] Ishit=%d \n",
	      t, vel, seg->Charge[st], sumc[0], is_hit );

  }
  while( !is_hit ); // Do until the end of drift

  (*seg).Xlength = pathlen;
  (*seg).Ylength = pathlen;
  (*seg).TTime = t*1E9; // [ns]
  for( int ipx = 0; ipx < 99; ++ipx )
    (*seg).TCharge[ipx] = sumc[ipx];
  (*seg).Steps = st;

  //delete EE;

  return;

} // Drift

//------------------------------------------------------------------------------
void KDetector::MipIR( Int_t ndiv )
{
  // a track is divided into ndiv pieces
  // pos and neg charges are drifted in 1 um steps to the electrodes
  // induced currents are calculated and stored per step

  double dist = TMath::Sqrt( TMath::Power( enp[0] - exp[0], 2 ) +
			     TMath::Power( enp[1] - exp[1], 2 ) +
			     TMath::Power( enp[2] - exp[2], 2 ) );

  bool ldb = 0;
  if( ldb )
    std::cout << "KDetector::MipIR"
	      << ", x " << enp[0] << " to " << exp[0]
	      << ", y " << enp[1] << " to " << exp[1]
	      << ", z " << enp[2] << " to " << exp[2]
	      << ", l " << dist
	      << ", ndiv " << ndiv
	      << ", MTresh " << MTresh
	      << ", average " << average
	      << std::endl << std::flush;

  if( ldb )
    std::cout
      << "drift bins " << pos->GetNbinsX()
      << ", " << pos->GetXaxis()->GetXmin()
      << ", " << pos->GetXaxis()->GetXmax()
      << std::endl << std::flush;

  TH1D * histop = new
    TH1D( "chp", "charge+ vs time",
	  pos->GetNbinsX(), pos->GetXaxis()->GetXmin(), pos->GetXaxis()->GetXmax() );
  TH1D * histon  = new
    TH1D( "chm", "charge- vs time",
	  neg->GetNbinsX(), neg->GetXaxis()->GetXmin(), neg->GetXaxis()->GetXmax() ); 

  sum->Reset();
  pos->Reset();
  neg->Reset();
  for( int ipx = 0; ipx < 99; ++ipx )
    qnode[ipx] = 0;

  double Q = dist*0.1; // 0.1 ke/um mean

  // energy loss distribution:
  // the sum of many Landaus becomes Gaussian (central limit theorem)
  // sample from overall Landau:

  TF1 * lan = new TF1( "lan"," TMath::Landau(x,[0],[1])", 0, 10*Q );
  lan->SetParameter( 0, Q ); // peak from mean
  lan->SetParameter( 1, Q/10/sqrt(dist/150) ); // width, adjusted DP 10/2017

  if( Landau )
    Q = lan->GetRandom();

  // distribute over steps:

  double qstp = Q/ndiv; // [ke] mean

  // re-sample with step-wise Landau:

  TF1 * lan2 = new TF1( "lan2"," TMath::Landau(x,[0],[1])", 0, 10*qstp );
  lan2->SetParameter( 0, qstp*11.25/15.5 ); // peak from mean
  lan2->SetParameter( 1, qstp/10 ); // width, adjusted

  for( int i = 0; i < ndiv; ++i ) { 

    Float_t sp[3];
    for( int j = 0; j < 3; ++j )
      sp[j] = enp[j] + ( i + 0.5 ) * ( exp[j] - enp[j] ) / ndiv; // point on track

    // printf( "#i=%d ndiv=%d, pointx=%f, pointy=%f pointz=%f \n",i,ndiv,sp[0],sp[1],sp[2]); 

    for( int j = 0; j < average; ++j ) { 

      KStruct seg;

      // holes:

      Drift( sp[0], sp[1], sp[2], 1, &seg );

      //double q = qstp; // DP 10/2017
      //double q = 1.0 / ndiv; // DP Dec 2018
      double q = 1.0; // DP Dec 2018
      if( Landau )
	q = lan2->GetRandom();
      //std::cout << "  " << q;

      //if( ldb )
	std::cout << "h drift from z " << sp[2]
		  << ", y " << sp[1]
		  << ": t " << seg.TTime
		  << ", S0 " << seg.TCharge[0]
		  << ", S1 " << seg.TCharge[1]
		  << ", S2 " << seg.TCharge[2]
	  //<< ", q " << seg.TCharge[0]*q
	  //<< ", N " << seg.Steps
		  << ", dist " << seg.Xlength
		  << std::endl;

      for( int ipx = 0; ipx < 99; ++ipx )
	qnode[ipx] += q*seg.TCharge[ipx]; // without trapping

      seg.GetCH( histop, 1, q, tauh ); // all nodes vs time, with trapping (KStruct)

      // e:

      Drift( sp[0], sp[1], sp[2], -1, &seg );

      if( ldb )
	std::cout << "e drift from z " << sp[2]
		  << ", y " << sp[1]
		  << ": t " << seg.TTime
		  << ", S " << seg.TCharge[0]
	  //<< ", q " << seg.TCharge[0]*q
	  //<< ", N " << seg.Steps
		  << ", d " << seg.Xlength
		  << std::endl;

      for( int ipx = 0; ipx < 99; ++ipx )
	qnode[ipx] += q*seg.TCharge[ipx]; // without trapping

      Float_t mule = 0;

      if( MTresh > 1 ) {
	mule = seg.GetCHMult( histon, 1, q, taue );  // performing multiplication
	// if(Debug) printf( ":: Mstep = %f ::",mule);
	std::cout << "e drift from z " << sp[2]
		  << ", y " << sp[1]
		  << ": t " << seg.TTime
		  << ", S0 " << seg.TCharge[0]
		  << ", S1 " << seg.TCharge[1]
		  << ", S2 " << seg.TCharge[2]
	  //<< ", q " << seg.TCharge[0]*q
	  //<< ", N " << seg.Steps
		  << ", dist " << seg.Xlength
		  << ", mul " << mule
		  << std::endl;
      }
      else
	seg.GetCH( histon, 1, q, taue ); // with trapping

      //if the multiplication is large enough then do the hole tracking:

      if( MTresh > 1 && mule > MTresh ) {
	if( ldb )
	  std::cout << "drift extra holes" << std::endl;
	for( int e = 1; e < seg.Steps+1; ++e ) {
	  KStruct segmul;
	  Drift( seg.Xtrack[e], seg.Ytrack[e], seg.Ztrack[e], 1, &segmul, seg.Time[e]*1E-9 );
	  Float_t mulh = segmul.GetCHMult( histop, 1, seg.MulCar[e], tauh );
	  if( mulh > BDTresh ) {
	    printf( "HOLE MULTIPLICATION - BREAKDOWN\n" );
	    BreakDown = 1;
	  }
	} // Steps
      } // MTresh

    } // average

    histop->Scale( 1 / ( (Float_t) average ) );
    histon->Scale( 1 / ( (Float_t) average ) ); 

    pos->Add(histop);
    neg->Add(histon);

    histop->Reset();
    histon->Reset();

  } // ndiv

  //std::cout << std::endl << std::flush;

  for( int ipx = 0; ipx < 99; ++ipx )
    qnode[ipx] /= (Float_t) average;

  sum->Add(neg);
  sum->Add(pos);

  delete histop;
  delete histon;

} // MipIR

//------------------------------------------------------------------------------
void KDetector::ShowMipIR( Int_t ndiv, Int_t color, Int_t how )
{
  // The simulation of the drift for the minimum ionizing particles. 
  // A track is devided into Int_ ndiv buckets. Each bucket is drifted in the field. The
  // induced currents for each carrier is calculated as the sum of all buckets. 

  TGraph *gr;
  TPolyLine3D *gr3D;
  Float_t sp[3];
  KStruct seg;
  TLine *line;

  std::cout << "KDetector::ShowMipIR"
	    << ", x " << enp[0] << " to " << exp[0]
	    << ", y " << enp[1] << " to " << exp[1]
	    << ", z " << enp[2] << " to " << exp[2]
	    << ", ndiv " << ndiv
	    << ", MTresh " << MTresh
	    << ", average " << average
	    << std::endl;

  // Draw histograms 

  if( EG != NULL ) { 

    if( nz == 1 )
      KHisProject( EG, 3, how )->Draw( "COL" );
    else {
      TH3D * hh = GetGeom();
      hh->SetFillColor(color);
      hh->SetTitle( "drift display" );
      hh->Draw( "iso" );
    }

  }

  // Draw drift paths

  for( int i = 0; i < ndiv; ++i ) {

    for( int j = 0; j < 3; ++j )
      sp[j] = ( ( exp[j] - enp[j] ) / ndiv ) *  i + enp[j] + ( exp[j] - enp[j] ) / (2*ndiv);

    if( Debug )
      printf( "drift start point %f %f %f \n", sp[0], sp[1], sp[2] );

    // hole drift

    Drift( sp[0],sp[1],sp[2],1, &seg );

    if( nz == 1 )
      gr = new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
    else
      gr3D = new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	

    if( nz == 1 ) {
      gr->SetLineColor(4);
      gr->SetLineStyle(1);
      gr->Draw( "L" );
    }
    else {
      gr3D->SetLineColor(4); 
      gr3D->SetLineStyle(1);  
      gr3D->Draw( "SAME" ); 
    }

    // electron drift:

    Drift( sp[0], sp[1], sp[2], -1, &seg );

    if(nz==1)
      gr = new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
    else
      gr3D = new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	

    if( nz == 1 ) {
      gr->SetLineColor(2);
      gr->SetLineStyle(3); 
      gr->Draw( "L" );
    }
    else {
      gr3D->SetLineColor(2); 
      gr3D->SetLineStyle(1);  
      gr3D->Draw( "SAME" ); 
    }

  }

} // ShowMipIR

//------------------------------------------------------------------------------
TH2D * KDetector::Draw( std::string option, Float_t proj )
{
  //The function draws weighting and electric field and also the event display
  // Char_t *option:
  // Which field:
  // 	  W  weighting  
  // 	  E  electric 
  //      G  geometry
  //      M  material  
  //      N  Neff = space charge
  // What component:
  // 	  P  potential 
  // 	  F  |field|
  //      X  x component of the field
  //      Y  y component of the field
  //      Z  z component of the field
  // In case of 3D which plane?
  //     yz  the cross section is along the x value of proj
  //     xz  the cross section is along the y value of proj  
  //     xy  the cross section is along the z value of proj
  // Float_t proj; position along the axis of projection

  KField *cf;
  TH3D   *ch;
  TH2D * histo = new TH2D();
  TString opt(option);

  Int_t i = 0;
  Int_t iproj = 0;
  Int_t What = 1;      // default is E
  Int_t Which = 0;     // default is potential
  Int_t Which3D = 3;   // default is XY

  if( opt.Contains( "N" ) ) What=5; // Neff
  if( opt.Contains( "M" ) ) What=4; // mat
  if( opt.Contains( "G" ) ) What=3; // Geo
  if( opt.Contains( "W" ) ) What=2; // weighting  
  if( opt.Contains( "E" ) ) What=1; // E

  if( nz > 1 ) {
    if( opt.Contains( "xy" ) ) Which3D = 3;
    if( opt.Contains( "xz" ) ) Which3D = 2;
    if( opt.Contains( "yz" ) ) Which3D = 1;
  }

  if( opt.Contains( "P" ) ) Which = 0; // potential
  if( opt.Contains( "F" ) ) Which = 1; // |E|
  if( opt.Contains( "X" ) ) Which = 2; // Ex
  if( opt.Contains( "Y" ) ) Which = 3;
  if( opt.Contains( "Z" ) ) Which = 4; // Ez

  if( What <= 2 ) {
    if( What == 2 )
      cf = &Ramo[0];
    else
      cf = &Real;
    switch(Which)
      {
      case 0: ch = cf->U;  break;
      case 1: ch = cf->E;  break;
      case 2: ch = cf->Ex; break;
      case 3: ch = cf->Ey; break;
      case 4: ch = cf->Ez; break;
      default: ch = cf->U;  break;
      }
  } 
  else {
    if( What == 3 ) ch = (TH3D *)EG;  
    if( What == 4 ) ch = (TH3D *)DM;  
    if( What == 5 ) ch = (TH3D *)NeffH;  
  }

  if( nz == 1 )
    histo = KHisProject( ch, 3, 1 ); 
  else {
    switch(Which3D)
      {
      case 3: iproj = ch->GetZaxis()->FindBin(proj); break;
      case 2: iproj = ch->GetYaxis()->FindBin(proj); break;
      case 1: iproj = ch->GetXaxis()->FindBin(proj); break;
      }
    histo = KHisProject( ch, Which3D, iproj );
    histo->SetStats(0);
    if( What == 1 ) {
      if( Which == 2 )
	histo->SetTitle( Form( "E_{X} at 3rd = %i", int(proj) ) );
      if( Which == 3 )
	histo->SetTitle( Form( "E_{Y} at 3rd = %i", int(proj) ) );
      if( Which == 4 )
	histo->SetTitle( Form( "E_{Z} at 3rd = %i", int(proj) ) );
    }
  }
  //histo->Draw( "COLZ" );

  return histo;

} // Draw

//------------------------------------------------------------------------------
TH1D * KDetector::Draw1D( std::string option, Float_t proj, Int_t axis, Float_t pos )
{
  // Draws the 1D projection of the Field
  // Char_t option;  see ::Draw()
  // Char_t proj;    see ::Draw() 
  // Int_t axis;     0=x, 1=y, 2=z;
  // Float_t pos;    position along the missing coordinate if 2D z=0.5;

  TH2D * h2 = Draw( option, proj );

  Int_t iter;
  Float_t low,up;

  if( axis == 0 ) {
    iter = h2->GetNbinsX(); 
    up = h2->GetXaxis()->GetBinUpEdge(iter);
    low = h2->GetXaxis()->GetBinLowEdge(1);
  }
  else {
    iter = h2->GetNbinsY();
    up = h2->GetYaxis()->GetBinUpEdge(iter);
    low = h2->GetYaxis()->GetBinLowEdge(1);
  }

  TH1D * h1 = new
    TH1D( Form( "field_%i", GetNhs() ),
	  Form( "field along dimension %i at %i and %i",
		axis+1, int(pos), int(proj) ),
	  iter, low, up );
  
  for( Int_t i = 1; i <= iter; ++i ) {
    if( axis == 0 )
      h1->SetBinContent( i, h2->GetBinContent(i, h2->GetYaxis()->FindBin(pos) ) );
    else
      h1->SetBinContent( i, h2->GetBinContent( h2->GetXaxis()->FindBin(pos), i ) );
  }

  return h1;

} // Draw1D

//------------------------------------------------------------------------------
KDetector::~KDetector()
{
  //Destructor of the detector class

  if(NeffF!=NULL) delete NeffF;
  if(NeffH!=NULL) delete NeffH;
  if(ran!=NULL) delete ran;
  if(pos!=NULL) delete pos;
  if(neg!=NULL) delete neg;
  if(sum!=NULL) delete sum;

}

//------------------------------------------------------------------------------
void KDetector::CalM( KStruct *seg, Double_t *data, Int_t CarrierType )
{
  Int_t i,j,numreg=0;
  Double_t *M=new Double_t [seg->Steps+1];
  Double_t *DIF=new Double_t [seg->Steps+1];
  Double_t dx,dy,sum=1.,sum2=0,sum1=0,dif,dif1,dif2;
  //      printf( "Number of Steps=%d\n",seg->Steps);
    
  for( i = 1; i < seg->Steps; ++i ) {
    //	  if(i==seg->Steps) DIF[i]=0;
    dx = TMath::Sqrt( TMath::Power( (seg->Xtrack[i+1] - seg->Xtrack[i] ), 2 ) +
		      TMath::Power( (seg->Ytrack[i+1] - seg->Ytrack[i] ), 2 ) );

    if( CarrierType < 0 )
      dif = Real.alpha( 0.5 * ( seg->Efield[i+1] + seg->Efield[i] ) );
    else
      dif = Real.beta( 0.5 * ( seg->Efield[i+1] + seg->Efield[i] ) );

    sum *= (1+dif*dx);
    DIF[i] = sum;
  }

  //      for(i=1;i<seg->Steps;i++) { printf( "Step=%d [%f %f], E=%f , a[i]=%f M(i)=%f\n",i,seg->Xtrack[i],seg->Ytrack[i],seg->Efield[i],Real.alpha(seg->Efield[i]),DIF[i]); }

  data[0] = DIF[seg->Steps-1];
  data[1]=0;
  data[2]=0;
  data[3]=0;
  // printf( "KKK %f\n",data[0]);
  for( i = 1; i < seg->Steps; ++i ) {
    if( (DIF[i]-1) / (data[0]-1) > 0.8 && DIF[i] > 1.02 ) {
      data[1] += seg->Xtrack[i];
      data[2] += seg->Ytrack[i]; 
      data[3] += seg->Time[i]*1E-9;
      ++numreg;
      //printf( "MUL=%f (X=%3.1f Y=%3.1f)::::: %e %e %e %d\n",DIF[i],seg->Xtrack[i],seg->Ytrack[i],data[1],data[2],data[3],numreg);
    }
  }
  if( numreg != 0 ) {
    data[1] /= numreg;
    data[2] /= numreg;
    data[3] /= numreg;
  }

} // CalM

//------------------------------------------------------------------------------
void KDetector::SetDriftHisto( Float_t x, Int_t numbins )
{
  if( pos != NULL ) delete pos;
  pos = new TH1D( "qpos", "hole signal", numbins, 0, x );
  if( neg != NULL ) delete neg;
  neg = new TH1D( "qneg", "electron signal", numbins, 0, x ); 	
  if( sum != NULL ) delete sum;
  sum = new TH1D( "charge", "sum signal", numbins, 0, x ); 

  sum->SetXTitle( "t [ns]" );
  neg->SetXTitle( "t [ns]" );
  pos->SetXTitle( "t [ns]" );
  sum->SetYTitle( "I [arb.]" );
  neg->SetYTitle( "I [arb.]" );
  pos->SetYTitle( "I [arb.]" );
  sum->GetYaxis()->SetTitleOffset(1.4);
  neg->GetYaxis()->SetTitleOffset(1.4);
  pos->GetYaxis()->SetTitleOffset(1.4);
  pos->SetLineColor(4);
  neg->SetLineColor(2);
  sum->SetLineColor(1);	

}

//------------------------------------------------------------------------------
void  KDetector::Save( std::string name, std::string file, bool all=1 )
{
  TFile * fn = new TFile( file.c_str(), "UPDATE" );

  std::cout << "Writing fields to " << file << std::endl;

  Char_t str[100];
  sprintf( str, "E_%s", name.c_str() );
  if( Real.U != NULL ) {
    if( all )
      Real.U->Write(str);
  }
  else {
    printf( "Histogram U (Electric) does not exist!\n" );
    return;
  }

  for( int ipx = 0; ipx < 99; ++ipx ) {
    sprintf( str, "W_%i_%s", ipx, name.c_str() );
    if( Ramo[ipx].U != NULL )
      Ramo[ipx].U->Write(str);
  }

  sprintf( str, "G_%s", name.c_str() );
  if( EG != NULL ) {
    if( all )
      EG->Write(str);
  }
  else {
    printf( "Histogram Geometry does not exist!\n" );
    return;
  }

  sprintf( str, "M_%s", name.c_str() );
  if( DM != NULL ) {
    if( all )
      DM->Write(str);
  }
  else {
    printf( "Histogram Material does not exist!\n" );
    return;
  }

  fn->Close();

} // Save

//------------------------------------------------------------------------------
TFile * KDetector::Read( std::string name, std::string file, bool all=1 )
{
  TFile * fn = new TFile( file.c_str() );

  Char_t str[100];

  if( Real.U == NULL )
    Real.U = new TH3D();
  if( all ) {
    sprintf( str, "E_%s", name.c_str() );
    std::cout << "read " << str << std::endl;
    Real.U->Read(str);
  }

  for( int ipx = 0; ipx < 99; ++ipx ) {
    sprintf( str,"W_%i_%s", ipx, name.c_str() );
    if( fn->GetKey(str) ) {
      if( Ramo[ipx].U == NULL )
	Ramo[ipx].U = new TH3D();
      std::cout << "read " << str << std::endl;
      Ramo[ipx].U->Read(str);
    }
    else
      Ramo[ipx].U = NULL;
  }

  if( EG == NULL )
    EG = new TH3I();
  if( all ) {
    sprintf( str, "G_%s", name.c_str() );
    std::cout << "read " << str << std::endl;
    EG->Read(str);
  }

  if( DM == NULL )
    DM = new TH3I();
  if( all ) {
    sprintf( str, "M_%s", name.c_str() );
    std::cout << "read " << str << std::endl;
    DM->Read(str);
  }

  nx = Ramo[0].U->GetNbinsX();
  ny = Ramo[0].U->GetNbinsY();
  nz = Ramo[0].U->GetNbinsZ();

  if( all )
    Real.CalField(); // E = -grd U

  //Ramo[0].CalField(); DP: not needed for KPixel

  return fn;
}
