
// Daniel Pitzl, Dec 2018
// edge-on, irradiated, 3 pix

// root -l edg3F6.C+

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg3F6( double Vbias = 450 )
{
  gStyle->SetTitleFont( 62, "XYZ" ); // 62 = Helvetica bold
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.2, "y" );
  gStyle->SetTitleOffset( 1.0, "z" );

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetHistFillStyle(0); // 1001 = solid

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  double dlt = 0.0009 * Vbias/250; // [cm]
  double N0 = 1620 * 250/Vbias; // [1/um^3] adjust at 800V, scaled 1/dlt^2 at fixed Vbias mdy 1.96

  if( Vbias < 620 )
    N0 = 1870 * 250/Vbias; // [1/um^3] adjust at 600V, scaled 1/dlt^2 at fixed Vbias, mdy 1.54

  if( Vbias < 520 )
    N0 =  900 * 250/Vbias; // [1/um^3] adjust at 500V, scaled 1/dlt^2 at fixed Vbias, mdy 1.58

  if( Vbias < 470 )
    N0 = 1930 * 250/Vbias; // [1/um^3] adjust at 450V, scaled 1/dlt^2 at fixed Vbias, mdy 2.5

  if( Vbias < 420 )
    N0 = 1940 * 250/Vbias; // [1/um^3] adjust at 400V, scaled 1/dlt^2 at fixed Vbias

  if( Vbias < 370 )
    N0 = 1950 * 250/Vbias; // [1/um^3] adjust at 350V, scaled 1/dlt^2 at fixed Vbias mdy 3.2

  if( Vbias < 320 )
    N0 = 1970 * 250/Vbias; // [1/um^3] adjust at 300V, scaled 1/dlt^2 at fixed Vbias mdy 1.43

  if( Vbias < 270 )
    N0 = 1990 * 250/Vbias; // [1/um^3] adjust at 250V, scaled 1/dlt^2 at fixed Vbias mdy 1.93

  if( Vbias < 220 )
    N0 = 2020 * 250/Vbias; // [1/um^3] adjust at 200V, scaled 1/dlt^2 at fixed Vbias mdy 1.28

  if( Vbias < 170 )
    N0 = 2050 * 250/Vbias; // [1/um^3] adjust at 150V, scaled 1/dlt^2 at fixed Vbias mdy 2.12

  double px =  50; // [um] double pixel for good weightfield
  double py =  50; // [um]
  //double pz = 158; // [um] even, adjusted 193

  double pz = 156; // [um] even, adjusted 174, 800 V
  if( Vbias < 620 )
    pz = 154; // [um] even, adjusted 174, 600 V

  // doping (space charge):

  double dj = 1.0; // 600V

  if( Vbias < 520 )
    dj = 0.0; // 500 V

  if( Vbias < 470 )
    dj = 0.6; // 450 V

  if( Vbias < 420 )
    dj = 0.5; // 350 V

  if( Vbias < 320 )
    dj = 0.2; // 250 V

  if( Vbias < 270 )
    dj = 0.0; // 250 V

  TF3 * f3;
  f3 = new
    TF3( "f3",
	 "x[0]*x[1]*0+[0]*pow(1-x[2]/[2],2)+[1]*exp(-([2]-x[2])/[3])", // exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );
  f3->SetParameter( 0, dj * pow( Vbias/150, 3 ) ); // [1/um^3] for double junction
  f3->SetParameter( 1, -N0 ); // [1/um^3]
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt*1E4 ); // range [um]

  // Geo:

  KPixel * det = new KPixel( 3, 3*px, py, pz ); // pixel, [um]

  det->SetUpVolume( 1.0, 2.0, 1.0 ); // um cubes (speed vs granularity)

  //                  center             half-implant    depth [um]
  det->SetUpPixel( 0, 0.5*px,  0.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 1, 1.5*px,  0.5*py, 0.45*px, 0.45*py, 1, 16385 ); // collecting electrode at Grnd
  det->SetUpPixel( 2, 2.5*px,  0.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd

  // bias applied at bot, negative to collect e on top
  det->Voltage = -Vbias;

  det->NeffF = f3;

  det->Temperature = 248;
  det->diff = 0;
  det->Landau = 0;
  det->MTresh = 1.01; // multiplication
  det->BDTresh = 1.5; // breakdown
  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  double F = 6.6; // [E15 p/cm2] fluence

  // e-trapping [ns]:

  //det->taue = 7.5/F; // [ns] from slope 193 at 800V
  det->taue = 6.9/F; // [ns] from slope 174 at 800V, 450 and below
  det->taue = 3.3/F; // [ns] from slope 174 at 500V

  // hole trapping: induces anti-charge

  det->tauh = 14/F; // [ns] tune 800V
  det->tauh = 18/F; // [ns] tune 800V

  TFile * fn = new
    TFile( Form( "edg3_px%d_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( px + 0.2 ), int( F + 0.1 ),
		 int( dlt*250/Vbias*1e4 + 0.2 ),
		 int(det->taue*F+0.01), int(det->tauh*F+0.01), int(-det->Voltage) ),
	   "RECREATE" );

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  TCanvas cEZ;
  cEZ.SetTitle("Ez map");
  det->Draw( "EZxz", 0.5*py )->Draw("COLZ");
  cEZ.Update();

  TH1D * field = det->Draw1D( "EZxz1", 0.5*py, 2, 1.5*px );

  cout << "E-field points " << field->GetNbinsX() << endl;
  double Temp = det->Temperature;
  cout << "Temp " << Temp << endl;

  TProfile * Nvsy = new TProfile( "nvsy", "space charge;height [#mum];space charge [1/#mum^{3}]",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );
  TProfile * Evsy = new TProfile( "evsy", "field;height [#mum];E field [V/#mum]",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );
  TProfile * Vvsy = new TProfile( "vvsy", "e drift velocity;height [#mum];e drift velocity [#mum/ns]",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );
  TProfile * Rvsy = new TProfile( "rvsy", "e range;height [#mum];e range [#mum]",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );

  for( int ii = 1; ii <= field->GetNbinsX(); ++ii ) {

    double z = field->GetBinCenter(ii);
    double E = field->GetBinContent(ii);
    double N = f3->Eval( 1.5*px, 0.5*py, z );
    double v = det->Real.DriftVelocity( fabs(E), -1, Temp, fabs(N), det->Mobility );
    cout << ii
      	 << "  " << z
	 << "  " << N
	 << "  " << E
	 << "  " << v*1E-5 // [um/ns]
	 << endl;

    double y = z - int(0.5*pz); // like data
    Nvsy->Fill( y, -N );
    Evsy->Fill( y, E );
    Vvsy->Fill( y, v*1E-5 );
    Rvsy->Fill( y, fabs(v)*1E-5*det->taue );

  }

  TCanvas cny;
  cny.SetTitle("N vs y");
  Nvsy->Draw();
  cny.Update();

  TCanvas cey;
  cey.SetTitle("E vs y");
  Evsy->Draw(); // evsy->Fit("expo","w")
  cey.Update();

  TCanvas cvy;
  cvy.SetTitle("V vs y");
  Vvsy->Draw();
  cvy.Update();

  TCanvas cry;
  cry.SetTitle("R vs y");
  Rvsy->Draw();
  cry.Update();

  TCanvas cWP;
  cWP.SetTitle("W map");
  det->Draw( "WPxz", 0.5*py )->Draw("COLZ");
  cWP.Update();

  TCanvas cWP1;
  cWP1.SetTitle("W vs z");
  TH1D * wpot = det->Draw1D( "WPxz1", 0.5*py, 2, 1.5*px );
  wpot->Draw();  
  cWP1.Update();

  TProfile * Wvsy = new TProfile( "wvsy", "weighting potential;height [#mum];weighting potential",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );

  for( int ii = 1; ii <= wpot->GetNbinsX(); ++ii ) {
    double z = wpot->GetBinCenter(ii);
    double w = wpot->GetBinContent(ii);
    double y = z - int(0.5*pz); // like data
    Wvsy->Fill( y, w );
  }
  TCanvas cwy;
  cwy.SetTitle("W vs y");
  Wvsy->Draw();
  cwy.Update();

  // track: along y

  det->enp[0] = 1.5*px; // mid
  det->enp[1] = 0.0*py; // entry point y
  det->enp[2] = 0.5*pz; // entry point z

  det->exp[0] = 1.5*px;
  det->exp[1] = 1.0*py;
  det->exp[2] = 0.5*pz;

  int ndiv = 10;
  TCanvas c1;
  c1.SetTitle("drift");
  det->ShowMipIR(ndiv);
  c1.Update();

  gStyle->SetOptStat(111111);

  TH1I * hq = new TH1I( "q", "charge;charge [ke];tracks", 400, 0, 40 );

  TProfile * Qvsy = new
    TProfile( "qvsy", "signal vs height;height [#mum];<readout signal> [relative]",
	      int(pz+0.5), -0.5*pz, 0.5*pz );

  TProfile * Pvsy = new
    TProfile( "pvsy", "signal vs height;height [mm];<readout signal> [relative]",
	      2*int(pz+0.5), -pz*1e-3, pz*1e-3 ); // like data

  TProfile * Svsy = new
    TProfile( "svsy", "smeared signal vs height;height [mm];<readout signal> [smeared]",
	      2*int(pz+0.5), -pz*1e-3, pz*1e-3 ); // like data

  int nev = 1;
  if( det->diff )
    nev = 9;
  if( det->Landau )
    nev = 25;

  TCanvas csteps;
  csteps.SetTitle("pulse");

  for( double z = 0.5; z < pz; z += 1 ) {

    det->enp[2] = z; // height
    det->exp[2] = z;

    double y = z - int(0.5*pz) - 2; // align 174 800 V
    if( Vbias < 620 )
      y = z - int(0.5*pz) - 0; // align 174 200 V

    double sumq = 0;

    for( int iev = 0; iev < nev; ++iev ) {

      det->MipIR(ndiv);
      double q = -det->sum->Integral(); // time integrated
      q /= ndiv; // norm to 1
      hq->Fill(q);
      Qvsy->Fill( y, q ); // Profile
      Pvsy->Fill( y*1E-3, q ); // [mm]
      sumq += q;

    } // events

    cout << z << "  " << sumq/nev << endl << flush;
    det->sum->Draw();
    csteps.Update();

  } // z

  Double_t smr = 0.010; // [mm] 

  double invsq2pi = 0.3989422804014; // 1/sqrt(2pi) for Gaussian norm
  double stp = Pvsy->GetBinWidth(1); // for normalization
  double norm = stp * invsq2pi / smr;

  // Convolution of spectrum with Gaussian:

  for( int ii = 1; ii <= Pvsy->GetNbinsX(); ++ii ) {

    Double_t yy = Pvsy->GetBinCenter(ii); // [mm]

    double sum = 0;

    for( int jj = 1; jj <= Pvsy->GetNbinsX(); ++jj ) {

      Double_t xx = Pvsy->GetBinCenter(jj);
      double qq = Pvsy->GetBinContent(jj);
      sum += qq * TMath::Gaus( yy, xx, smr );

    } // jj

    Svsy->Fill( yy, norm * sum ); // [mm]

  } // ii

  TCanvas cq;
  cq.SetTitle("q");
  hq->Draw();
  cq.Update();

  TCanvas cqy;
  cqy.SetTitle("Q vs y");
  Qvsy->Draw();
  cqy.Update();

  fn->Write(); // histos
  cout << fn->GetName() << endl;

  delete det;

  cout << "N0 " << N0*Vbias/250 << endl;
  cout << "Emax = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;

}
