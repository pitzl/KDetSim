
// Daniel Pitzl, Dec 2018
// edge-on, irradiated, 3 pix

// root -l edg3F4.C+

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg3F4( double Vbias = 500 )
{
  gStyle->SetTitleFont( 62, "XYZ" ); // 62 = Helvetica bold
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.2, "y" );
  gStyle->SetTitleOffset( 1.0, "z" );

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetHistFillStyle(0); // 1001 = solid

  //double Vbias = 1200;
  //double Vbias = 800;
  //double Vbias = 600;
  //double Vbias = 500; // full depletion at 8E15 n/cm2
  //double Vbias = 400; // convex at taue=4/F
  //double Vbias = 300;
  //double Vbias = 250; // half depletion at 8E15 n/cm2
  //double Vbias = 200;
  //double Vbias = 150;
  //double Vbias = 100; not sensitive to taue
  //double Vbias =  50;

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  int nmode = 3;

  double dlt;
  double N0;

  if( nmode == 2 ) { // 1/sqrt space charge, d^1.5 potential, E-field not dropping fast enough
    dlt = 20;
    N0 = 260/sqrt(dlt); // full depletion at 500V for F=4
  }

  else if( nmode == 4 ) { // 1/z^2 space charge, field spike, mult, breakdown
    N0 = 160; // full depletion at 500V for F=4
    dlt = 20; // [um]
  }

  else if( nmode == 3 ) { // exp(-z/l)

    if( Vbias < 120 ) {
      dlt = 0.0016*Vbias/250; // [cm]
      N0 = 650*250/Vbias; // [1/um^3] tune at 100V
    }
    else if( Vbias < 170 ) {
      dlt = 0.0015*Vbias/250; // [cm]
      N0 = 734*250/Vbias; // [1/um^3] tune at 150V
    }
    else if( Vbias < 320 ) {
      dlt = 0.00121*Vbias/250; // [cm]
      N0 = 1121*250/Vbias; // [1/um^3] adjust at 300V, scaled 1/dlt^2 at fixed Vbias
    }
    else if( Vbias < 420 ) {
      dlt = 0.0011*Vbias/250; // [cm]
      N0 = 1356*250/Vbias; // [1/um^3] adjust at 400V, scaled 1/dlt^2 at fixed Vbias
    }
    else if( Vbias < 820 ) {
      dlt = 0.0010*Vbias/250; // [cm]
      N0 = 1640*250/Vbias; // [1/um^3] adjust at 500V, scaled 1/dlt^2 at fixed Vbias
    }
    else {
      dlt = 0.0009*Vbias/250; // [cm]
      N0 = 2025*250/Vbias; // [1/um^3] adjust at 600V, scaled 1/dlt^2 at fixed Vbias
    }

  } // mode

  double px =  50; // [um] double pixel for good weightfield
  double py =  50; // [um]
  double pz = 150; // [um] adjusted

  // doping (space charge):

  TF3 * f3;
  if( nmode == 2 ) { // 1/sqrt
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*pow(([3]-x[2])/[2],-0.5)", // 1/sqrt
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 0 ); // no offset, no  double junction
    f3->SetParameter( 1, -N0 ); // sqrt(cm) -> sqrt(um), negative for n in p
    f3->SetParameter( 2, dlt );
    f3->SetParameter( 3, pz ); // end of depletion
  }

  else if( nmode == 4 ) { // 1/z^2
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*pow(([3]-x[2])/[2],-2)", // 1/z^2
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 0 ); // no offset, no  double junction
    f3->SetParameter( 1, -N0 ); // sqrt(cm) -> sqrt(um), negative for n in p
    f3->SetParameter( 2, dlt );
    f3->SetParameter( 3, pz ); // end of depletion
  }

  else if( nmode == 3 ) { // exp
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*exp(-([2]-x[2])/[3])", // exp
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    //f3->SetParameter( 0, 1*pow(Vbias/400,3) ); // [1/um^3] for double junction for 400V
    f3->SetParameter( 0, 0 ); // 
    f3->SetParameter( 1, -N0 ); // [1/um^3]
    f3->SetParameter( 2, pz ); // thickness
    f3->SetParameter( 3, dlt*1E4 ); // range [um]

  }

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

  double F = 4; // [E15 n/cm2] fluence

  // e-trapping [ns]:

  det->taue = 4.3/F; // [ns] from slope, sensitive
  if( Vbias > 720 )
    det->taue = 5.1/F; // [ns] from slope at 800 V with tauh 6

  // hole trapping: induces anti-charge

  //det->tauh = 4/F; // [ns] tune 200V
  //det->tauh = 6/F; // [ns] tune 800V less curvature
  det->tauh = 8/F; // [ns] tune 800V shifted right

  TFile * fn;
  if( nmode == 3 )
    fn = new
      TFile( Form( "edg3_px%d_F%d_dlt%d_e%d_h%d_V%d.root",
		   int( px + 0.2 ), int( F + 0.1 ),
		   int( dlt*250/Vbias*1e4 + 0.2 ),
		   int(det->taue*F+0.01), int(det->tauh*F+0.01), int(-det->Voltage) ),
	     "RECREATE" );
  else
    fn = new
      TFile( Form( "edg3_px%d_F%d_mode%d_e%d_h%d_V%d.root",
		   int(px+0.2), int(F+0.1), nmode,
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
  TH1D * wpot = det->Draw1D( "WPxz1", 0.5*py, 2, 1.5*px );
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

    double y = int(0.5*pz) - z; // like data

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

  cout << "N0 " << N0 << endl;
  if( nmode == 3 )
    cout << "Emax = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;
}
