
// Daniel Pitzl, Dec 2018
// edge-on, irradiated, 3 pix

// root -l edg3t100.C+

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );


#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg3t100( double Vbias = 800 )
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

  double dlt = 0.0010*Vbias/250; // [cm]
  double N0 = 1647*250/Vbias; // [1/um^3] adjust at 400V, scaled 1/dlt^2 at fixed Vbias

  double px =  50; // [um] influences shallow region at 150V: weightfield of trapped holes
  double py =  50; // [um]
  double pz = 102; // [um] compensate implants

  // doping (space charge):

  TF3 * f3 = new
    TF3( "f3",
	 "x[0]*x[1]*0+[0]+[1]*exp(-([2]-x[2])/[3])", // exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );
  f3->SetParameter( 0, 1*pow(Vbias/500,3) ); // [1/um^3] for double junction for 500V
  f3->SetParameter( 1, -N0 ); // [1/um^3]
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt*1E4 ); // range [um]

  // Geo:

  KPixel * det = new KPixel( 3, 3*px, py, pz ); // pixel, [um]

  det->SetUpVolume( 1.0, 2.0, 1.0 ); // um cubes (speed vs granularity)

  //                  center             half-implant    [um]  pot
  det->SetUpPixel( 0, 0.5*px,  0.5*py, 0.45*px, 0.45*py, 0.5,     1 ); // Grnd
  det->SetUpPixel( 1, 1.5*px,  0.5*py, 0.45*px, 0.45*py, 0.5, 16385 ); // collecting electrode at Grnd
  det->SetUpPixel( 2, 2.5*px,  0.5*py, 0.45*px, 0.45*py, 0.5,     1 ); // Grnd

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

  double F = 8; // [E15 n/cm2] fluence
  F = 16;
  F = 32;

  // e-trapping [ns]:

  det->taue = 4/F; // [ns] from slope at 800 V

  // hole trapping: induces anti-charge

  det->tauh =  4/F; // [ns]

  TFile * fn = new
    TFile( Form( "edg3_t100_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( F + 0.1 ), int( dlt*250/Vbias*1e4 + 0.2 ),
		 int(det->taue*F+0.01), int(det->tauh*F+0.01), int(-det->Voltage) ),
	   "RECREATE" );

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  TCanvas cEZ;
  cEZ.SetTitle("Ez map");
  det->Draw( "EZxz", 0.5*py )->Draw("COLZ");
  cEZ.Update();

  TH1F * field = det->Draw1D( "EZxz1", 0.5*py, 2, 1.5*px );

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

  for( int ii = 2; ii < field->GetNbinsX(); ++ii ) {

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

    double y = int(0.5*pz) - z; // like data
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
  TH1F * wpot = det->Draw1D( "WPxz1", 0.5*py, 2, 1.5*px );
  wpot->Draw();  
  cWP1.Update();

  TProfile * Wvsy = new TProfile( "wvsy", "weighting potential;height [#mum];weighting potential",
				  int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );

  for( int ii = 1; ii <= wpot->GetNbinsX(); ++ii ) {
    double z = wpot->GetBinCenter(ii);
    double w = wpot->GetBinContent(ii);
    double y = int(0.5*pz) - z;
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
	      int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );
  TProfile * Pvsy = new
    TProfile( "pvsy", "signal vs height;height [mm];<readout signal> [relative]",
	      2*int(pz+0.5)+1, -(pz+0.5)*1e-3, (pz+0.5)*1e-3 ); // like data
  TProfile * Svsy = new
    TProfile( "svsy", "signal vs smeared height;smeared height [mm];<readout signal> [relative]",
	      2*int(pz+0.5)+1, -(pz+0.5)*1e-3, (pz+0.5)*1e-3 ); // like data

  int nev = 1;
  if( det->diff )
    nev = 9;
  if( det->Landau )
    nev = 25;

  TCanvas csteps;
  csteps.SetTitle("pulse");

  for( double z = 1; z < pz; z += 1 ) {

    det->enp[2] = z; // height
    det->exp[2] = z;

    double y = int(0.5*pz) - z; // like data

    double sumq = 0;

    for( int iev = 0; iev < nev; ++iev ) {

      det->MipIR(ndiv);
      double q = -det->sum->Integral(); // time integrated
      hq->Fill(q);
      Qvsy->Fill( y, q ); // Profile
      Pvsy->Fill( y*1E-3, q ); // [mm]
      sumq += q;

    } // events

    cout << z << "  " << sumq/nev << endl << flush;
    det->sum->Draw();
    csteps.Update();

  } // z

  Double_t smr = 0.011; // [mm] tuned at 250V

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

  cout << "Emax = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;
}
