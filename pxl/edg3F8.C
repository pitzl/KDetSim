
// Daniel Pitzl, Dec 2018
// edge-on, irradiated, 3 pix

// root -l edg3F8.C+

// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "../inc/KPixel.h"

void edg3F8( double Vbias = 800 )
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

  double eps_0 = 8.854187817e-14; // [F/cm]
  double eps_Si = 11.7;
  double eps = eps_0 * eps_Si; // 1.0359e-12 F/cm
  double q0 = 1.60217733e-19; // [C]

  int nmode = 3;

  double N0;
  double depl;
  double dlt = 0;

  if( nmode == 0 ) {  // constant space charge: quadratic potential, poor model for qvsy

    N0 = 80; // [1/um^3] = [1E12/cm3] tuned at 250V
    if( Vbias > 470 )
      N0 = 60; // [1/um^3] = [1E12/cm3] tuned at 500V
    if( Vbias > 570 )
      N0 = 48; // [1/um^3] = [1E12/cm3] tuned at 600V

    depl = sqrt( 2*eps*Vbias / ( N0*1e12 * q0 ) )*1E4; // [um]
  }

  else if( nmode == 1 ) { // linear space charge: cubic potential, too deep depletion at 150V

    if(      Vbias > 770 )
      N0 =   1/0.015; // [1e12/cm3/cm] tuned at 800V mdy 4.69
    else if( Vbias > 570 )
      N0 = 110/0.015; // [1e12/cm3/cm] tuned at 600V
    else if( Vbias > 470 )
      N0 =  220/0.015; // [1e12/cm3/cm] tuned at 500V
    else if( Vbias > 370 )
      N0 = 300/0.015; // [1e12/cm3/cm] tune  at 400V
    else if( Vbias > 220 )
      N0 = 600/0.015; // [1e12/cm3/cm] tuned at 250V
    else
      N0 = 350/0.015; // [1e12/cm3/cm] tune  at 150V

    depl = pow( 6*eps*Vbias / ( N0*1e12 * q0 ), 1.0/3.0 )*1E4; // [um]
  }

  else if( nmode == 2 ) { // 1/sqrt space charge, d^1.5 potential, E-field not dropping fast enough
    N0 = 1.32; // full depletion at 500V for F=8
    depl = pow( 3.0/4.0*eps*Vbias / ( N0*1e12 * q0 ), 2.0/3.0 )*1E4; // [um]
  }

  else if( nmode == 3 ) { // exp(-z/l)

    //dlt = 0.0012 * Vbias/250; // [cm]
    //N0 = 1148 * 250/Vbias; // [1/um^3] adjust at 400V
    //N0 = 1150 * 250/Vbias; // [1/um^3] adjust at 300V
    //N0 = 1155 * 250/Vbias; // [1/um^3] tune at 250V
    //N0 = 1163 * 250/Vbias; // [1/um^3] tune at 200V
    //N0 = 1177 * 250/Vbias; // [1/um^3] tune at 150V

    dlt = 0.0011 * Vbias/250; // [cm] best mdy at 250V
    N0 = 1371 * 250/Vbias; // [1/um^3] adjust at 250V, scaled 1/dlt^2 at fixed Vbias

    if( Vbias > 270 ) {
      dlt = 0.0010 * Vbias/250; // [cm]
      //N0 = 1664 * 250/Vbias; // [1/um^3] adjust at 250V, scaled 1/dlt^2 at fixed Vbias
      N0 = 1647 * 250/Vbias; // [1/um^3] adjust at 400V, scaled 1/dlt^2 at fixed Vbias
    }

    depl = dlt*log(eps*Vbias/(dlt*dlt*q0*N0))*1e4; // [um]

  } // mode

  cout << "depletion depth " << depl << " um" << endl;

  //double px =  25; // [um] see weightfield for trapped holes at 150V
  double px =  50; // [um] influences shallow region at 150V: weightfield of trapped holes
  double py =  50; // [um]
  double pz = 152; // [um] compensate implants

  double z0 = pz - depl; // [um] begin of depleted region for pure n
  if( nmode == 3 )
    z0 = 0;

  // doping (space charge):

  TF3 * f3;

  if( nmode == 0 ) { // const
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*fmax(0,(x[2]-[2])/fabs(x[2]-[2]))", // truncated const
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 0.5 ); // trapped holes 250V
    if( Vbias > 470 )
      f3->SetParameter( 0, 3 ); // trapped holes 500V
    f3->SetParameter( 1,-N0 ); // n in p [1/um^3], 150 um at 600 V
    f3->SetParameter( 2, z0 ); // end of depletion
  }

  else if( nmode == 1 ) { // linear: from e trapping
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*fmax(0,x[2]-[2])", // linear
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 0.6 ); // trapped holes 400V
    if( Vbias > 470 )
      f3->SetParameter( 0, 2 ); // trapped holes 500V
    //f3->SetParameter( 0, 1*fmax( 0, depl-pz ) ); // trapped holes, scale with overdepletion
    f3->SetParameter( 1, -N0*1e-4 ); // /cm -> /um, negative for n in p
    f3->SetParameter( 2, z0 ); // end of depletion
  }

  else if( nmode == 2 ) { // 1/sqrt
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*((x[2]-[2]>0.1)?1/sqrt(152-x[2]):0)", // 1/sqrt
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 0 ); // no offset, no  double junction
    f3->SetParameter( 1, -N0*100 ); // sqrt(cm) -> sqrt(um), negative for n in p
    f3->SetParameter( 2, z0 ); // end of depletion
  }

  else if( nmode == 3 ) { // exp
    f3 = new
      TF3( "f3",
	   "x[0]*x[1]*0+[0]+[1]*exp(-([2]-x[2])/[3])", // exp
	   -1000, 1000, -1000, 1000, -1000, 1000 );
    f3->SetParameter( 0, 1*pow( Vbias/500, 3 ) ); // [1/um^3] for double junction for 500V
    //f3->SetParameter( 0, 0*pow( Vbias/500, 3 ) ); // [1/um^3] for double junction for 500V
    f3->SetParameter( 1, -N0 ); // [1/um^3] negative for n in p
    f3->SetParameter( 2, pz ); // thickness
    f3->SetParameter( 3, dlt*1E4 ); // range [um]
  }

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
  //F = 16;
  //F = 32;

  // e-trapping [ns]:

  //det->taue = 8/F; // [ns]
  //det->taue = 5/F; // [ns] from slope at 800 V, sensitive
  det->taue = 4/F; // [ns] from slope at 800 V
  if( Vbias > 770 )
    det->taue = 4.5/F; // [ns] from slope at 800 V, sensitive
  //det->taue = 2/F; // [ns]

  // hole trapping: induces anti-charge

  //det->tauh = 200/F; // [ns] shallow charge 100%
  //det->tauh = 20/F; // [ns]
  //det->tauh =  5/F; // [ns]
  det->tauh =  4/F; // [ns]
  //det->tauh =  3/F; // [ns]
  //det->tauh =  2.5/F; // [ns]
  //det->tauh =  2/F; // [ns]
  //det->tauh =  1/F; // [ns] requires drift time truncation at 10 ns

  if( nmode == 3 && det->Mobility == 4 ) { // expo with Jacoboni
    //det->taue = 6.9/F; // [ns] 800V mdy 2.22
    det->taue = 6.8/F; // [ns] 800V mdy 2.00
    //det->taue = 6.6/F; // [ns] 800V mdy 2.06
    //det->taue = 6.4/F; // [ns] 800V mdy 2.53

    //det->tauh = 7/F; // [ns] 2.22
    det->tauh = 6/F; // [ns] 2.00
    //det->tauh = 5/F; // [ns] 2.06
    //det->tauh = 4/F; // [ns]  mdy 2.53
  }

  if( nmode == 0 ) { // const
    det->taue = 2.4/F; // [ns] mdy 

    det->tauh = 2/F; // [ns] 3.66
    if( det->Mobility == 4 )
      det->tauh = 0.5/F; // [ns] 
  }

  if( nmode == 1 ) { // linear

    //det->taue = 2/F; // [ns] mdy 7.41
    det->taue = 3/F; // [ns] mdy 4.44
    //det->taue = 4/F; // [ns] mdy 4.69
    //det->taue = 5/F; // [ns] mdy 5.80

    det->tauh = 3/F; // [ns] 4.44
    //det->tauh = 4/F; // [ns] 4.47
    //det->tauh = 6/F; // [ns] 4.88
  }

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

  // truncate negative E-field
  /*
  for( int i = 1; i <= det->nx; ++i )
    for( int j = 1; j <= det->ny; ++j )
      for( int k = 2; k < det->nz; ++k ) { // first and last are electrodes
	double Ek = det->Real.Ez->GetBinContent( i, j, k );
	if( Ek < 0 )
	  det->Real.Ez->SetBinContent( i, j, k, 0.02 ); // [V/um] low field region
	  //det->Real.Ez->SetBinContent( i, j, k, 0.005 ); // low field region
	double EE = sqrt( pow( det->Real.Ex->GetBinContent( i, j, k ), 2 ) +
			  pow( det->Real.Ey->GetBinContent( i, j, k ), 2 ) +
			  pow( det->Real.Ez->GetBinContent( i, j, k ), 2 ) );
	det->Real.E->SetBinContent( i, j, k, EE );
      }
  */
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
    TProfile( "svsy", "signal vs smeared height;smeared height [mm];<readout signal> [relative]",
	      2*int(pz+0.5), -pz*1e-3, pz*1e-3 ); // like data

  int nev = 1;
  if( det->diff )
    nev = 9;
  if( det->Landau )
    nev = 25;

  TCanvas csteps;
  csteps.SetTitle("pulse");

  //for( double z = fmax( 0.5, z0 ); z < pz; z += 1 ) {
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

  //Double_t smr = 0.012; // [mm] tuned at 250V
  Double_t smr = 0.011; // [mm] tuned at 250V
  //Double_t smr = 0.010; // [mm] tuned to shallow charge at 800V
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

  cout
    << "N0 " << N0
    << ", depletion depth " << depl << " um"
    << ", z0 " << z0
    << endl;
  if( nmode == 3 )
    cout << "Emax = " << dlt*1E4*q0*N0/eps*1E4 << " V/um" << endl;
  cout << "track smearing " << smr*1e3 << " um" << endl;

  cout << "enter any character to end..." << endl;
  string ans;
  cin >> ans;
}
