
// Daniel Pitzl, Jan 2019
// edge-on, irradiated, 3x3 pix

// root -l edg9F6.C

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg9F6( double Vbias = 800 )
{
  gStyle->SetTitleFont( 62, "XYZ" ); // 62 = Helvetica bold
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.2, "y" );
  gStyle->SetTitleOffset( 1.0, "z" );

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetHistFillStyle(0); // 1001 = solid

  gStyle->SetStatX(0.5); // right edge of Stat box
  gStyle->SetStatY(0.9); // top   edge of Stat box

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  // Geo: doRamo if you change!

  double px =  50; // [um] combine two, for shallow charge at 800V
  double py = 100; // [um]
  double pz = 154; // [um] doRamo !

  KPixel * det = new KPixel( 9, 3*px, 3*py, pz ); // pixel, [um]

  det->SetUpVolume( 1.0, 2.0, 1.0 ); // um cubes (speed vs granularity)

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

  double dlt = 0.0009 * Vbias/250; // [cm]

  double N0 = 2010 * 250/Vbias; // [1/um^3] tune at 150V: dj 0, mob4: mdy 1.24

  if( Vbias > 180 )
    N0 = 1982 * 250/Vbias; // [1/um^3] tune at 200V: dj 0, mob4: mdy 1.06

  if( Vbias > 230 )
    N0 = 1960 * 250/Vbias; // [1/um^3] tune at 250V: dj 0, mob4: mdy 1.86

  if( Vbias > 280 )
    N0 = 1955 * 250/Vbias; // [1/um^3] tune at 300V: dj 0, mob4: mdy 1.71

  if( Vbias > 330 )
    N0 = 1945 * 250/Vbias; // [1/um^3] tune at 350V: dj 10, mob4: mdy 3.39, 3.47

  if( Vbias > 380 )
    //N0 = 1930 * 250/Vbias; // [1/um^3] tune at 400V with dj 10, mob4:  mdy 1.96
    N0 = 1910 * 250/Vbias; // [1/um^3] tune at 400V with dj 10, mob3:  mdy 2.16

  if( Vbias > 430 )
    //N0 = 1920 * 250/Vbias; // [1/um^3] tune at 450V with dj 15, mob4: mdy 2.00
    N0 = 1910 * 250/Vbias; // [1/um^3] tune at 450V with dj 15, mob3: mdy 2.11

  if( Vbias > 480 )
    //N0 = 1880 * 250/Vbias; // [1/um^3] tune at 500V, dj 30, tauh 28, mob3: mdy 2.74
    //N0 = 1850 * 250/Vbias; // [1/um^3] tune at 500V, dj 40, tauh 28, mob3: mdy 2.64
    N0 = 1820 * 250/Vbias; // [1/um^3] tune at 500V, dj 50, tauh 28, mob3: mdy 2.54
    //N0 = 1780 * 250/Vbias; // [1/um^3] tune at 500V, dj 60, tauh 28, mob3: mdy 2.66
    //N0 = 1870 * 250/Vbias; // [1/um^3] tune at 500V, dj 50, pow3, tauh 28, mob4: mdy 2.77
    //N0 = 1870 * 250/Vbias; // [1/um^3] tune at 500V, dj 50, pow4, tauh 28, mob4: mdy 3.11
    //N0 = 1790 * 250/Vbias; // [1/um^3] tune at 500V, dj 70, pow3, tauh 28, mob4: mdy 2.87
    //N0 = 1830 * 250/Vbias; // [1/um^3] tune at 500V, dj 70, pow4, tauh 28, mob4: mdy 2.90

  if( Vbias > 580 )
    //N0 = 1910 * 250/Vbias; // [1/um^3] tune at 600V, dj 50, tauh 21: mdy 2.21
    //N0 = 1390 * 250/Vbias; // [1/um^3] tune at 600V, dj 150, tauh 14: mdy 1.43
    //N0 = 1320 * 250/Vbias; // [1/um^3] tune at 600V, dj 200, tauh 28, mob4: mdy 1.04
    N0 = 1305 * 250/Vbias; // [1/um^3] tune at 600V, dj 200, tauh 28, mob3: mdy 1.13

  if( Vbias > 680 )
    //N0 = 1960 * 250/Vbias; // [1/um^3] tune at 700V, dj 50, tauh 21: mdy 2.73
    //N0 = 1250 * 250/Vbias; // [1/um^3] tune at 700V, dj 250, tauh 21: mdy 2.06
    N0 = 1050 * 250/Vbias; // [1/um^3] tune at 700V, dj 290, tauh 21: mdy 1.94

  if( Vbias > 780 )
    //N0 = 1300 * 250/Vbias; // [1/um^3] tune at 800V, dj 290, tauh 21, mob3: mdy 2.10
    N0 =  900 * 250/Vbias; // [1/um^3] tune at 800V, dj 350, tauh 21, mob4: mdy 2.05
    //N0 =  300 * 250/Vbias; // [1/um^3] tune at 800V, dj 450, tauh 28, mob4: mdy bad
    //N0 = 1900 * 250/Vbias; // [1/um^3] tune at 800V, dj 50, tauh 21, mob4: mdy 2.61

  // double junction amplitude:

  double dj = 0; // 200V

  if( Vbias > 230 )
    dj = 0.6;
  if( Vbias > 280 )
    dj = 2;
  if( Vbias > 330 )
    dj = 5; // mdy 3.47
  if( Vbias > 380 )
    dj = 10;
  if( Vbias > 430 )
    //dj = 18; // mdy 2.10
    dj = 20;
  if( Vbias > 480 )
    //dj = 30; // mdy 2.74
    //dj = 40; // mdy 2.64
    dj = 50; // pow3, mdy 2.54
  //dj = 60; // pow3, mdy 2.66
  //dj = 50; // pow4, mdy 3.11
  //dj = 70; // pow4, mdy 2.87
  if( Vbias > 580 )
    dj = 200;
  if( Vbias > 680 )
    //dj = 50;
    //dj = 250; // mdy 2.06
    dj = 290; // mdy 1.94
  if( Vbias > 780 )
    //dj = 50;
    dj = 350;
  //dj = 450;

  // doping (space charge):

  TF3 * f3 = new
    TF3( "f3",
	 "x[0]*x[1]*0+[0]*pow(1-x[2]/[2],3)+[1]*exp(-([2]-x[2])/[3])", // dj + exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );

  f3->SetParameter( 0, dj ); // [1/um^3] for double junction for 500V
  f3->SetParameter( 1, -N0 ); // [1/um^3]
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt*1E4 ); // range [um]

  double F = 6.6; // [E15 p/cm2] fluence

  // e-trapping [ns]:

  det->taue = 8.8/F; // [ns] from slope at 800 V, mob3: mdy 1.82
  //det->taue = 11.9/F; // [ns] from slope at 800 V, tauh 28, dj 350, mob4: mdy 2.05
  //det->taue = 11.4/F; // [ns] from slope at 800 V, tauh 21, dj 350, mob4: mdy 2.10
  //det->taue = 12.8/F; // [ns] from slope at 800 V, dj 50, mob4: mdy 2.61

  // hole trapping [ns]:

  //det->tauh =  7/F; // [ns] at 500V
  //det->tauh = 14/F; // [ns] at 600 V
  //det->tauh = 21/F; // [ns] at 800 V
  det->tauh = 28/F; // [ns] at 600 V

  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  // bias applied at bot, negative to collect e on top
  det->Voltage = -Vbias;

  det->NeffF = f3;

  det->Temperature = 258; // Sep 2018
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
	       Form( "edg9_px%d_Ramo.root", int( px + 0.2 ) ), // file
	       0 ); // only Ramo
  }
  else // save time
    det->Read( "edg9", // model
	       Form( "edg9_px%d_Ramo.root", int( px + 0.2 ) ), // file
	       0 ); // only Ramo

  // det->Read opens a root file
  // new file must be defined afterwards...

  TFile * fn = new
    TFile( Form( "edg9_px%d_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( px + 0.2 ), int( F + 0.1 ),
		 int( dlt*250/Vbias*1e4 + 0.2 ),
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
				  int(pz+0.1), -0.5*pz, 0.5*pz );
  TProfile * Evsy = new TProfile( "evsy", "field;height [#mum];E field [V/#mum]",
				  int(pz+0.1), -0.5*pz, 0.5*pz );
  TProfile * Vvsy = new TProfile( "vvsy", "e drift velocity;height [#mum];e drift velocity [#mum/ns]",
				  int(pz+0.1), -0.5*pz, 0.5*pz );
  TProfile * Rvsy = new TProfile( "rvsy", "e range;height [#mum];e range [#mum]",
				  int(pz+0.1), -0.5*pz, 0.5*pz );

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
				  int(pz+0.1), -0.5*pz, 0.5*pz );

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
  det->exp[1] = 3.0*py;
  det->exp[2] = 0.5*pz;

  int ndiv = 30;
  TCanvas c1;
  c1.SetTitle("drift");
  det->ShowMipIR(ndiv);
  c1.Update();

  gStyle->SetOptStat(111111);

  TH1I * hq = new TH1I( "q", "charge;charge [ke];tracks", 400, 0, 40 );

  TProfile * Qvsy = new
    TProfile( "qvsy", "signal vs height;height [#mum];<readout signal> [relative]",
	      int(pz+0.1), -0.5*pz, 0.5*pz );

  TProfile * Pvsy = new
    TProfile( "pvsy", "signal vs height;height [mm];<readout signal> [relative]",
	      2*int(pz+0.1), -pz*1e-3, pz*1e-3 ); // like data

  TProfile * Svsy = new
    TProfile( "svsy", "signal vs smeared height;smeared height [mm];<readout signal> [relative]",
	      2*int(pz+0.1), -pz*1e-3, pz*1e-3 ); // like data

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
    if( Vbias < 720 )
      y = z - int(0.5*pz) - 0; // align 174 200 V

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

  Double_t smr = 0.010; // [mm] tuned to shallow charge at 800V

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

  cout << "N0 " << N0*Vbias/250 << endl;
  cout << "Emax = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;

}
