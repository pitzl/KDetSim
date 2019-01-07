
// Daniel Pitzl, Dec 2018
// edge-on, irradiated

// root -l edg9.C

{
  double Vbias = 800;
  double dlt = 0.0010*Vbias/250; // [cm]
  double Neff = 1643*250/Vbias; // [1/um^3] adjust at 250V

  double px =  50; // [um] combine two, for shallow charge at 800V
  double py = 100; // [um]
  double pz = 150; // [um]

  TF3 * f3 = new  // doping (space charge):
    TF3( "f3",
	 "x[0]*x[1]*0+[0]+[1]*exp(-([2]-x[2])/[3])", // exp
	 -1000, 1000, -1000, 1000, -1000, 1000 );

  f3->SetParameter( 0, 1*pow(Vbias/500,3) ); // [1/um^3] for double junction for 500V
  f3->SetParameter( 1, -Neff ); // [1/um^3] negative for n in p
  f3->SetParameter( 2, pz ); // thickness
  f3->SetParameter( 3, dlt*1E4 ); // range [um]

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

  det->Voltage = -Vbias;  // bias applied at bot, negative to collect e on top
  det->NeffF = f3;

  det->Temperature = 248;
  det->diff = 0;
  det->Landau = 0;
  det->MTresh = 1.01; // multiplication
  det->BDTresh = 1.5; // breakdown
  det->Mobility = 3; // mobility model: 1=Canali, 3=Scharf, 4=Jacoboni (KField)

  double F = 8; // [E15 n/cm2] fluence
  det->taue = 4/F; // [ns] e-trapping
  det->tauh = 3/F; // [ns] hole trapping

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  // track: along y

  det->enp[0] = 1.5*px; // mid
  det->enp[1] = 0.0*py; // entry point y
  det->enp[2] = 0.5*pz; // entry point z

  det->exp[0] = 1.5*px;
  det->exp[1] = 3.0*py;
  det->exp[2] = 0.5*pz;

  int ndiv = 30;
  int nev = 1;

  TFile * fn = new
    TFile( Form( "edg9_px%d_F%d_dlt%d_e%d_h%d_V%d.root",
		 int( px + 0.2 ), int( F + 0.1 ),
		 int( dlt*250/Vbias*1e4 + 0.2 ),
		 int( det->taue*F + 0.01 ), int( det->tauh*F + 0.01 ),
		 int( -det->Voltage+0.5 ) ),
	   "RECREATE" );

  TProfile * Pvsy = new
    TProfile( "pvsy", "signal vs height;height [mm];<readout signal> [relative]",
	      2*int(pz+0.5)+1, -(pz+0.5)*1e-3, (pz+0.5)*1e-3 ); // like data
  TProfile * Svsy = new
    TProfile( "svsy", "signal vs smeared height;smeared height [mm];<readout signal> [relative]",
	      2*int(pz+0.5)+1, -(pz+0.5)*1e-3, (pz+0.5)*1e-3 ); // like data

  for( double z = 1; z < pz-1; z += 1 ) {

    det->enp[2] = z; // height
    det->exp[2] = z;

    double y = 0.5*pz - z; // like data

    for( int iev = 0; iev < nev; ++iev ) {

      det->MipIR(ndiv); // track

      double q = -det->sum->Integral(); // time integrated
      q *= 3; // norm to middle pixel
      Pvsy->Fill( y*1E-3, q ); // [mm]

    } // events

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

  fn->Write(); // histos
  cout << fn->GetName() << endl;

  delete det;

}
