
// Daniel Pitzl, Jan 2019
// edge-on, irradiated, 3x3 pix phase I

// root -l edg9F2.C

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg9F2( double Vbias = 800 )
{
  gStyle->SetTitleFont( 62, "XYZ" ); // 62 = Helvetica bold
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.2, "y" );
  gStyle->SetTitleOffset( 1.0, "z" );

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetHistFillStyle(0); // 1001 = solid

  gStyle->SetOptStat(1);
  gStyle->SetStatX(0.4); // right edge of Stat box
  gStyle->SetStatY(0.9); // top   edge of Stat box
  gROOT->ForceStyle();

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  // Geo: doRamo if you change!

  double px = 450; // [um] // double px
  double py = 100; // [um] shoot along 100
  //double pz = 288; // [um] double even, doRamo !
  double pz = 284; // [um] double even, doRamo !

  KPixel * det = new KPixel( 9, 3*px, 3*py, pz ); // pixel, [um]

  det->SetUpVolume( 6, 4, 2 ); // um cubes (speed vs granularity)

  //                  center             half-implant    depth [um]
  det->SetUpPixel( 0, 0.5*px,  0.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 1, 1.5*px,  0.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 2, 2.5*px,  0.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd

  det->SetUpPixel( 3, 0.5*px,  1.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 4, 1.5*px,  1.5*py, 0.45*px, 0.45*py, 1, 16385 ); // collecting electrode at Grnd
  det->SetUpPixel( 5, 2.5*px,  1.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd

  det->SetUpPixel( 6, 0.5*px,  2.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 7, 1.5*px,  2.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd
  det->SetUpPixel( 8, 2.5*px,  2.5*py, 0.45*px, 0.45*py, 1,     1 ); // Grnd

  // N0 exp(-z/dlt):

  double dlt = 11 * Vbias/250; // [um] at 275 V
  //double dlt = 16 * Vbias/250; // [um] at 325 V

  //double N0 = 2240; // [1/um^3] at 225V, dtl 9, pow 4, scale 1/dlt^2 at fixed Vbias
  //double N0 = 2290; // [1/um^3] at 225V, dlt 9, pow 2, scale 1/dlt^2 at fixed Vbias
  //double N0 = 1855; // [1/um^3] at 225V, dlt 10, pow 2, scale 1/dlt^2 at fixed Vbias
  //double N0 = 1505; // [1/um^3] at 225V, dlt 11, pow 3, scale 1/dlt^2 at fixed Vbias
  double N0 = 1520; // [1/um^3] at 225V, dlt 11, pow 2, scale 1/dlt^2 at fixed Vbias: mdy 0.254

  if( Vbias > 270 )
    N0 = 1216; // dlt 11, pow 3: mdy 0.217 ke

  if( Vbias > 320 )
    N0 =  996; // dlt 11, pwr 4, dj 5.3: mdy 0.195
  //N0 = 1010; // dlt 11, pwr 3: mdy 0.194

  if( Vbias > 390 )
    N0 = 741; // dlt 11, pwr 5: mdy 0.155

  if( Vbias > 490 )
    N0 = 492; // dlt 11, pwr 6: mdy 0.109 ke

  if( Vbias > 590 )
    N0 = 249; // dlt 11, pwr 6: mdy 0.079

  if( Vbias > 790 )
    //N0 = 109; // dlt 15, taue 10, tauh 28: mdy 0.134
    N0 = 127; // pwr 6: mdy 0.133

  // double junction amplitude:

  double pwr = 2.0; // position of minimum
  if( Vbias > 270 )
    pwr = 3.0;
  if( Vbias > 320 )
    pwr = 4.0;
  if( Vbias > 390 )
    pwr = 5.0;
  if( Vbias > 490 )
    pwr = 6.0;

  //double dj = 2.4; // pow 4 at 225V
  double dj = 1.0; // pow 3 at 225V: mdy 0.254
  //double dj = 0.8; // pow 2 at 225V: mdy 0.337

  if( Vbias > 270 )
    //dj = 3.6; // pow 4: mdy 0.218 ke
    dj = 2.2; // pow 3: mdy 0.217 ke
    //dj = 1.4; // pow 2

  if( Vbias > 320 )
    dj = 5.3; // pwr 4: mdy 0.195
  //dj = 3.0; // pwr 3

  if( Vbias > 390 )
    dj = 20; // pwr 5

  if( Vbias > 490 )
    dj =  70; // pwr 6

  if( Vbias > 590 )
    dj =  168; // pwr 6: mdy 0.

  if( Vbias > 790 )
    dj = 286; // pwr 6
    //dj = 170; // dlt 15, taue 11, tauh 28: mdy 0.122

  // doping (space charge):

  TF3 * f3 = new
    TF3( "f3",
	 "x[0]*x[1]*0+[0]*pow(1-x[2]/[2],[4])+[1]*exp(-([2]-x[2])/[3])", // dj + exp
	 //"x[0]*x[1]*0+[0]*exp(-x[2]/[4])+[1]*exp(-([2]-x[2])/[3])", // dj + exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );

  f3->SetParameter( 0, dj ); // [1/um^3] for double junction for 500V
  f3->SetParameter( 1, -N0 ); // [1/um^3]
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt ); // width [um] n-side
  f3->SetParameter( 4, pwr ); // width [um] p-side

  double F = 1.9; // [E15 p/cm2] fluence

  // e-trapping [ns]:

  det->taue = 10.1/F; // [ns] from left peak at 600V
  //det->taue =  7.5/F; // [ns] at 275V, left peak not smaller, right peak narrower
  //det->taue = 11.1/F; // [ns] slope at 800 V

  // hole trapping [ns]:

  det->tauh = 28/F; // [ns] tune 500 V

  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  // bias applied at bot, negative to collect e on top
  det->Voltage = -Vbias;

  det->NeffF = f3;

  det->Temperature = 263; // Mar 2017
  det->diff = 0;
  det->Landau = 0;
  det->MTresh = 1.01; // multiplication
  det->BDTresh = 1.5; // breakdown

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  det->CalField(0); // E-field

  bool doRamo = 0; // do this once for each geometry

  if( doRamo ) {
    det->CalField(1); // Ramo weight field

    det->Save( "edg9", // model
	       Form( "edg9_px%d_pz%d_Ramo.root",
		     int( px + 0.2 ),int( pz + 0.2 ) ), // file
	       0 ); // only Ramo
  }
  else // save time
    det->Read( "edg9", // model
	       Form( "edg9_px%d_pz%d_Ramo.root",
		     int( px + 0.2 ),int( pz + 0.2 ) ), // file
	       0 ); // only Ramo

  // det->Read opens a root file
  // new file must be defined afterwards...

  TFile * fn = new
    TFile( Form( "edg9_px%d_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( px + 0.2 ), int( F + 0.1 ),
		 int( dlt*250/Vbias + 0.2 ),
		 int( det->taue*F + 0.01 ), int( det->tauh*F + 0.01 ),
		 int( -det->Voltage+0.5 ) ),
	   "RECREATE" );

  TCanvas cEZ;
  cEZ.SetTitle("Ez map");
  det->Draw( "EZxz", 1.5*py )->Draw("COLZ");
  cEZ.Update();

  TH1D * field = det->Draw1D( "EZxz1", 1.5*py, 2, 1.5*px );

  cout << "E-field points " << field->GetNbinsX() << endl;
  double Temp = det->Temperature;
  cout << "Temp " << Temp << endl;

  TProfile * Nvsy = new TProfile( "nvsy", "space charge;height [#mum];space charge [1/#mum^{3}]",
				  int(pz+0.1)/2, -0.5*pz, 0.5*pz );
  TProfile * Evsy = new TProfile( "evsy", "field;height [#mum];E field [V/#mum]",
				  int(pz+0.1)/2, -0.5*pz, 0.5*pz );
  TProfile * Vvsy = new TProfile( "vvsy", "e drift velocity;height [#mum];e drift velocity [#mum/ns]",
				  int(pz+0.1)/2, -0.5*pz, 0.5*pz );
  TProfile * Rvsy = new TProfile( "rvsy", "e range;height [#mum];e range [#mum]",
				  int(pz+0.1)/2, -0.5*pz, 0.5*pz );

  for( int ii = 1; ii <= field->GetNbinsX(); ++ii ) {

    double z = field->GetBinCenter(ii);
    double E = field->GetBinContent(ii);
    double N = f3->Eval(1.5*px,1.5*py,z);
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
  det->Draw( "WPxz", 1.5*py )->Draw("COLZ");
  cWP.Update();

  TCanvas cWP1;
  cWP1.SetTitle("W vs z");
  TH1D * wpot = det->Draw1D( "WPxz1", 1.5*py, 2, 1.5*px );
  wpot->Draw();  
  cWP1.Update();

  TProfile * Wvsy = new TProfile( "wvsy", "weighting potential;height [#mum];weighting potential",
				  int(pz+0.1)/2, -0.5*pz, 0.5*pz );

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
  det->enp[2] = 0.8*pz; // entry point z

  det->exp[0] = 1.5*px;
  det->exp[1] = 3.0*py;
  det->exp[2] = 0.8*pz;

  int ndiv = 30;
  TCanvas c1;
  c1.SetTitle("drift");
  det->ShowMipIR(ndiv);
  c1.Update();

  TH1I * hq = new TH1I( "q", "charge;charge [ke];tracks", 400, 0, 40 );

  TProfile * Qvsy = new
    TProfile( "qvsy", "signal vs height;height [#mum];<readout signal> [relative]",
	      int(pz+0.1)/2, -0.5*pz, 0.5*pz );

  TProfile * Pvsy = new
    TProfile( "pvsy", "signal vs height;height [mm];<readout signal> [relative]",
	      2*int(pz+0.1)/2, -pz*1e-3, pz*1e-3 ); // like data

  TProfile * Svsy = new
    TProfile( "svsy", "signal vs smeared height;smeared height [mm];<readout signal> [relative]",
	      2*int(pz+0.1)/2, -pz*1e-3, pz*1e-3 ); // like data

  int nev = 1;
  if( det->diff )
    nev = 9;
  if( det->Landau )
    nev = 25;

  TCanvas csteps;
  csteps.SetTitle("pulse");

  for( double z = 1; z < pz; z += 2 ) {

    det->enp[2] = z; // height
    det->exp[2] = z;

    double y = z - int(0.5*pz); // shift (don't invert)

    double sumq = 0;

    for( int iev = 0; iev < nev; ++iev ) {

      det->MipIR(ndiv);
      double q = -det->sum->Integral(); // time integrated
      q *= 3; // norm to middle pixel
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

  Double_t smr = 0.020; // [mm] width at 325V

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

  TCanvas csy;
  csy.SetTitle("S vs y");
  Svsy->Draw();
  csy.Update();

  fn->Write(); // histos
  cout << fn->GetName() << endl;

  delete det;

  cout << "N0 " << N0 << endl;
  cout << "Emax = " << dlt*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;

}
