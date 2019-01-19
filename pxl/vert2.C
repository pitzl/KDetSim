
// Daniel Pitzl, Jan 2019
// model for thicker phase I sensors, vertical incidence
// 3x3 px

// root -l vert2.C // 1.9E15p/cm2

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

/* interactive
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void vert2( double Vbias = 800 )
*/
{
  double Vbias = 800;

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
  double N0 = 127; // [1/um^3] at 800V

  // double junction amplitude:

  double pwr = 6.0; // 800V
  double dj = 286; // pwr 6

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

  // bias applied at bot, negative to collect e on top
  det->Voltage = -Vbias;

  det->NeffF = f3;

  det->Temperature = 263;
  det->diff = 0;
  det->Landau = 0;
  det->MTresh = 1.01; // multiplication
  det->BDTresh = 1.5; // breakdown
  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  TFile * fn = new
    TFile( Form( "vert2_pz%d_dlt%d_V%d.root",
		 int( pz+0.2 ),
		 int( dlt*250/Vbias + 0.2 ), int(-det->Voltage+0.1) ),
	   "RECREATE" );

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

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

  for( int ii = 1; ii < field->GetNbinsX(); ++ii ) {

    double z = field->GetBinCenter(ii);
    double E = field->GetBinContent(ii);
    double N = f3->Eval( 1.5*px, 1.5*py, z );
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
    double y = int(0.5*pz) - z;
    Wvsy->Fill( y, w );
  }
  TCanvas cwy;
  cwy.SetTitle("W vs y");
  Wvsy->Draw();
  cwy.Update();

  // track: along z

  det->enp[0] = 1.5*px; // mid
  det->enp[1] = 1.5*py; // entry point y
  det->enp[2] = 0.0*pz; // entry point z

  det->exp[0] = 1.5*px;
  det->exp[1] = 1.5*py;
  det->exp[2] = 1.0*pz;

  int ndiv = 30;
  TCanvas c1;
  c1.SetTitle("drift");
  det->ShowMipIR(ndiv);
  c1.Update();

  int nev = 1;
  if( det->diff )
    nev = 9;
  if( det->Landau )
    nev = 25;

  TCanvas csteps;
  csteps.SetTitle("pulse");

  TProfile * CvsF = new
    TProfile( "cvsf", "charge collection;fluence [10^{15 }p/cm^{2}];<readout signal> [relative]",
	      43, -0.25, 21.25 );
  TProfile * QvsF = new
    TProfile( "qvsf", "charge collection;fluence [10^{15 }p/cm^{2}];<charge in 285 #mum> [ke]",
	      43, -0.25, 21.25 );

  // trapping: fluence scan (same field)

  for( double F = 0; F < 21.25; F += 0.5 ) {// [E15 p/cm2] fluence

    det->taue = 10/F; // [ns] from edg9F2
    det->tauh = 28/F; // [ns]

    double sumq = 0;

    for( int iev = 0; iev < nev; ++iev ) {

      det->MipIR(ndiv);
      det->sum->Draw();
      csteps.Update();
      double q = -det->sum->Integral(); // time integrated
      q /= ndiv; // norm to 1
      CvsF->Fill( F, q );
      QvsF->Fill( F, q*21.5 ); // ke in 285 um
      sumq += q;

    } // events

    cout << F << "  " << sumq/nev << endl << flush;

  } // F

  TCanvas cqf;
  cqf.SetTitle("Q vs F");
  QvsF->Draw();
  cqf.Update();

  fn->Write(); // histos
  cout << fn->GetName() << endl;

  delete det;

  cout << "Emax = " << dlt*q0*N0/eps*1E4 << " V/um" << endl;

}
