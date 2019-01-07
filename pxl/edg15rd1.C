
// Daniel Pitzl, Dec 2018
// edge-on, irradiated, 5x3 pix
// read one px, scan dx

// root -l edg15r1.C

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
/*
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg15F8( double Vbias = 300 ) */
{
  gStyle->SetTitleFont( 62, "XYZ" ); // 62 = Helvetica bold
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.2, "y" );
  gStyle->SetTitleOffset( 1.0, "z" );

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetHistFillStyle(0); // 1001 = solid

  bool doField = 0;

  double Vbias = 800;

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  double dlt = 0.0011 * Vbias/250; // [cm] for 250 V
  double N0 = 1338 * 250/Vbias; // [1/um^3] adjust at 250V, scaled 1/dlt^2 at fixed Vbias

  double px =  25; // [um] true
  double py = 100; // [um]
  double pz = 152; // [um] compensate implants

  double F = 8; // [E15 n/cm2] fluence

  // e-trapping [ns]:

  taue = 2.7/F; // [ns] F8 at 800V

  // hole trapping: induces anti-charge

  tauh = 1.9/F; // [ns] F8 at 800V

  // doping (space charge):

  double dj = 0.6; // [1/um^3] double junction for 500V
  if( Vbias < 270 ) dj = 0;

  TF3 * f3 = new
    TF3( "f3",
	 "x[0]*x[1]*0+[0]+[1]*exp(-([2]-x[2])/[3])", // exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );
  //f3->SetParameter( 0, dj*pow( Vbias/500, 3 ) ); // [1/um^3] for double junction
  f3->SetParameter( 0, dj*pow( Vbias/500, 5 ) ); // [1/um^3] for double junction
  f3->SetParameter( 1, -N0 ); // [1/um^3] negative for n in p
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt*1E4 ); // range [um]

  // Geo:

  KPixel * det = new KPixel( 15, 5*px, 3*py, pz ); // pixel, [um] inherits KDetector

  det->SetUpVolume( 1.0, 2.0, 1.0 ); // um cubes (speed vs granularity) KPixel: EG and DM

  //                  center           half-implant   depth pot [um]
  det->SetUpPixel(  0, 0.5*px, 0.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  1, 1.5*px, 0.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  2, 2.5*px, 0.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  3, 3.5*px, 0.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  4, 4.5*px, 0.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd

  det->SetUpPixel(  5, 0.5*px, 1.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  6, 1.5*px, 1.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd, will be readout 0
  det->SetUpPixel(  7, 2.5*px, 1.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  8, 3.5*px, 1.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel(  9, 4.5*px, 1.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd

  det->SetUpPixel( 10, 0.5*px, 2.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel( 11, 1.5*px, 2.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel( 12, 2.5*px, 2.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel( 13, 3.5*px, 2.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd
  det->SetUpPixel( 14, 4.5*px, 2.5*py, 0.45*px, 0.45*py, 1, 1 ); // Grnd

  // bias applied at bot, negative to collect e on top
  det->Voltage = -Vbias;

  det->NeffF = f3;
  det->taue = taue;
  det->tauh = tauh;

  det->Temperature = 248;
  det->diff = 0;
  det->Landau = 0;
  det->MTresh = 1.01; // multiplication
  det->BDTresh = 1.5; // breakdown
  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  if( doField ) {
    det->SetUpElectrodes(); // KPixel, surface contacts
    det->SetBoundaryConditions(); // KGeometry, fills EG
    det->CalField(0); // KDetector E-field
  }
  else {
    det->Read( "edg15F8", // model
	       Form( "edg15_px%d_rd1_F%d_dlt%d_V%d_fields.root",
		     int( px + 0.2 ), int( F + 0.1 ),
		     int( dlt*250/Vbias*1e4 + 0.2 ),
		     int(Vbias+0.1) ) ); // model and file
  }

  // det->Read opens a root file
  // new file must be defined afterwards...

  TFile * fn = new
    TFile( Form( "edg15_px%d_rd1_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( px + 0.2 ), int( F + 0.1 ),
		 int( dlt*250/Vbias*1e4 + 0.2 ),
		 int(taue*F+0.01), int(tauh*F+0.01), int(Vbias+0.1) ),
	   "RECREATE" );

  TCanvas cEZ;
  cEZ.SetTitle("Ez map");
  det->Draw( "EZxz", 1.5*py )->Draw("COLZ");
  cEZ.Update();

  TH1F * field = det->Draw1D( "EZxz1", 1.5*py, 2, 2.5*px );

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
    double N = f3->Eval( 2.5*px, 1.5*py, z );
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

  // Ramo weighting field, one readout pixel at a time:

  if( doField ) {

    int iRamo = 0;
    for( int ipx = 6; ipx <= 6; ++ipx ) {
      det->SetPixelW( ipx, 16385 ); // bits 14 and 1: readout node at Grnd
      det->SetUpElectrodes(); // in KPixel
      det->SetBoundaryConditions();
      det->CalField(1+iRamo); // Ramo weighting potential
      det->SetPixelW( ipx, 1 ); // back to Grnd
      ++iRamo;
    }

    det->Save( "edg15F8", // model
	       Form( "edg15_px%d_rd1_F%d_dlt%d_V%d_fields.root",
		     int( px + 0.2 ), int( F + 0.1 ),
		     int( dlt*250/Vbias*1e4 + 0.2 ),
		     int(Vbias+0.1) ) ); // model and file
    fn->cd();

  } // doField

  TCanvas cWP;
  cWP.SetTitle("W map");
  det->Draw( "WPxz", 1.5*py )->Draw("COLZ"); // for Ramo[0]
  cWP.Update();

  TCanvas cWP1;
  cWP1.SetTitle("W vs z");
  TH1F * wpot1 = det->Draw1D( "WPxz1", 1.5*py, 2, 1.5*px ); // Ramo[0] left
  wpot1->Draw();  
  TH1F * wpot2 = det->Draw1D( "WPxz1", 1.5*py, 2, 2.5*px ); // Ramo[0] mid
  wpot2->Draw("same");  
  cWP1.Update();

  TProfile * W1vsy = new
    TProfile( "w1vsy", "this weighting potential;height [#mum];this weighting potential",
	      int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );
  TProfile * W2vsy = new
    TProfile( "w2vsy", "next weighting potential;height [#mum];next weighting potential",
	      int(pz+0.5)+1, -0.5*(pz+1), 0.5*(pz+1) );

  for( int ii = 1; ii <= wpot1->GetNbinsX(); ++ii ) {
    double z = wpot1->GetBinCenter(ii);
    double y = int(0.5*pz) - z;
    W1vsy->Fill( y, wpot1->GetBinContent(ii) );
    W2vsy->Fill( y, wpot2->GetBinContent(ii) );
  }
  TCanvas cwy;
  cwy.SetTitle("W vs y");
  W1vsy->SetMaximum(1);
  W1vsy->Draw();
  W2vsy->Draw("same");
  cwy.Update();

  // track: along y

  det->enp[0] = 1.5*px; // under readout
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

  gStyle->SetOptStat(111111);

  TProfile2D * Qvsdxz = new
    TProfile2D( "qvsdxz", "signal vs dx and z;dx [#mum];z [#mum];<readout signal> [relative]",
		151, -75.5, 75.5, int(pz+0.1)+60, -30, pz+30 ); // space for smearing

  TProfile2D * Svsdxz = new
    TProfile2D( "svsdxz",
	      "signal vs smeared dx and z;smeared dx [#mum];smeared z [#mum];<readout signal> [relative]",
	      151, -75.5, 75.5, int(pz+0.1)+60, -30, pz+30 );

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

    for( double dx = 0; dx < 75.5; dx += 1 ) {

      det->enp[0] = dx + 1.5*px;
      det->exp[0] = dx + 1.5*px;

      double sumq = 0;

      for( int iev = 0; iev < nev; ++iev ) {

	det->MipIR(ndiv);
	double q = -det->sum->Integral(); // time integrated, summed readout nodes, after trapping
	q *= 3./ndiv; // norm to middle readout row
	Qvsdxz->Fill( dx, z, q ); // Profile
	Qvsdxz->Fill(-dx, z, q ); // symmetric
	sumq += q;

      } // events

      cout << dx << "  " << sumq/nev << endl << flush;
      det->sum->Draw();
      csteps.Update();

    } // dx

  } // z

  TCanvas cqm;
  cqm.SetTitle("Q map");
  Qvsdxz->Draw("colz");
  cqm.Update();

  double smr = 10; // [um]
  double smr2 = sqrt(2)*smr; // 2D
  double invsq2pi = 0.3989422804014; // 1/sqrt(2pi) for Gaussian norm
  double stpx = Qvsdxz->GetXaxis()->GetBinWidth(1); // for normalization
  double stpz = Qvsdxz->GetYaxis()->GetBinWidth(1); // for normalization
  double norm = stpx*stpz * invsq2pi / smr2/smr2; // 2D Gaussian

  // Convolution of map with 2-D Gaussian:

  for( int ix = 1; ix <= Qvsdxz->GetNbinsX(); ++ix ) {

    Double_t dx = Qvsdxz->GetXaxis()->GetBinCenter(ix); // [mm]

    for( int iz = 1; iz <= Qvsdxz->GetNbinsY(); ++iz ) {

      Double_t z = Qvsdxz->GetYaxis()->GetBinCenter(iz); // [mm]

      double sum = 0;

      for( int jx = 1; jx <= Qvsdxz->GetNbinsX(); ++jx ) {

	Double_t xx = Qvsdxz->GetXaxis()->GetBinCenter(jx);

	for( int jz = 1; jz <= Qvsdxz->GetNbinsY(); ++jz ) {

	  Double_t zz = Qvsdxz->GetYaxis()->GetBinCenter(jz);

	  double qq = Qvsdxz->GetBinContent(jx,jz);
	  sum += qq * exp( -0.5 * ( pow( (xx-dx)/smr2, 2 ) + pow( (zz-z)/smr2, 2 ) ) );

	} // jz

      } // jx

      Svsdxz->Fill( dx, z, norm * sum ); // [mm]

    } // iz

  } // ix

  TCanvas csm;
  csm.SetTitle("S map");
  Svsdxz->Draw("colz");
  csm.Update();

  cout << "dlt   = " << dlt*250/Vbias*1e4 << " um at 250 V" << endl;
  cout << "N0    = " << N0*Vbias/250 << " /um^3 at 250 V" << endl;
  cout << "Emax  = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "taue  = " << taue*F << " ns/cm^2" << endl;
  cout << "tauh  = " << tauh*F << " ns/cm^2" << endl;
  cout << "track = " << smr*1e3 << " um smearing" << endl;

  fn->Write(); // histos
  cout << fn->GetName() << endl;

  delete det;

}
