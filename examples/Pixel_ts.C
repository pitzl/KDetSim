
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .L fitep0.C
// .x Pixel_ts.C

// tilt scan

{
  TFile * fn = new TFile( "tilt_scan.root", "RECREATE" );

  double px = 100; // [um]
  double py =  60; // [um]
  double pz = 285; // [um]
  KPixel * det = new KPixel( 5, 5*px, 1*py, pz ); // 3 pixel, [um]

  det->SetUpVolume( 5, 5, 3 ); // um cubes (speed vs granularity)

  det->SetUpPixel( 0, 0.5*px,  0.5*py, 0.4*px, 0.5*py-1, 1, 1 );
  det->SetUpPixel( 1, 1.5*px,  0.5*py, 0.4*px, 0.5*py-1, 1, 1 );
  det->SetUpPixel( 2, 2.5*px,  0.5*py, 0.4*px, 0.5*py-1, 1, 1 );
  det->SetUpPixel( 3, 3.5*px,  0.5*py, 0.4*px, 0.5*py-1, 1, 1 );
  det->SetUpPixel( 4, 4.5*px,  0.5*py, 0.4*px, 0.5*py-1, 1, 1 );

  // bias applied at bot, negative to collect e on top
  //det->Voltage = -67; // just full depletion for Neff = 1, d = 300
  det->Voltage = -150; // above full depletion for Neff = 1, d = 300

  // E-field unit is V/um

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  // doping (space charge):

  TF3 * f3 = new TF3( "f3", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 300, 0, 500 );
  f3->SetParameter( 0, 1 ); // n in n [1/um^3]
  det->NeffF = f3;

  det->CalField(0); // E-field

  TCanvas cEZ2;
  det->Draw( "EZxz2", 50 )->Draw("COLZ");
  cEZ2.Update();

  TCanvas cEZ1;
  det->Draw1D( "EZxz1", 50, 2, 150 )->Draw();
  //Projection->SetTitle("E_{z} vs z");
  cEZ1.Update();

  // 1st readout node:
  det->SetPixelW( 1, 16385 ); // readout node at Grnd

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  det->CalField(1); // Ramo weight field

  TCanvas cWP;
  det->Draw( "WPxz", 50 )->Draw("COLZ");
  cWP.Update();

  // 2nd readout node:
  det->SetPixelW( 1,     1 );
  det->SetPixelW( 2, 16385 ); // 2nd readout node

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  det->CalField(2); // 2nd Ramo weighting potential

  // 3rd readout node:
  det->SetPixelW( 2,     1 );
  det->SetPixelW( 3, 16385 ); // 3rd readout node

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  det->CalField(3); // 3rd Ramo weighting potential

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // track:

  double pi = 4*atan(1);
  double wt = 180/pi;

  det->enp[0] = 2.5*px; // entry point x [um]
  det->enp[1] = 0.5*py; // entry point y
  det->enp[2] =  2; // entry point z

  det->exp[0] = 2.5*px;
  det->exp[1] = 0.5*py;
  det->exp[2] = pz-2;

  det->Temperature = 293; // [K]
  det->diff = 1;
  det->SStep = 5; // [um] drift step
  det->Landau = 1;

  double thr = 1.5; // [ke]

  int nev = 16;

  TProfile * sxvst = new
    TProfile( "sxvst",
	      "#sigma x vs tilt;tilt angle [deg];resolution [um]",
	      41, -0.5, 40.5, 0, 999 );

  for( int tilt = 0; tilt < 40; ++tilt ) {

    double ang = tilt/wt; // flat, optimum, at geometric minimum

    int ndiv = pz/cos(ang)/5+1; // steps

    cout << "track angle " << ang*wt << endl;
    cout << "track divisions " << ndiv << endl;

    TH1I * hqpx = new TH1I( Form( "qpx_%i", tilt ),
			    Form( "pixel charge tilt %i;pixel charge [ke];pixel", tilt ),
			    100, 0, 50 );
    TH1I * hq = new TH1I( Form( "q_%i", tilt ),
			  Form( "charge tilt %i;charge [ke];tracks", tilt ),
			  100, 0, 100 );
    TH1I * hdx = new TH1I( Form( "dx_%i", tilt ),
			   Form( "residual tilt %i;#Deltax [um];tracks", tilt ),
			   100, -50, 50 );
    TProfile * dxvsx = new
      TProfile( Form( "dxvsx_%i", tilt ),
		Form( "#Deltax vs x tilt %i;x [um];<#Deltax [um]", tilt ),
		101, 1.5*px-0.5, 2.5*px+0.5, -999, 999 );

  // step across one pixel in x:

    for( double x = 1.5*px; x <= 2.5*px; x += 1 ) {

      det->enp[0] = x - 0.5*pz*tan(ang);
      det->exp[0] = x + 0.5*pz*tan(ang);

      for( int iev = 0; iev < nev; ++iev ) {

	det->MipIR(ndiv);

	double q0 = -det->qnode[0];
	double q1 = -det->qnode[1];
	double q2 = -det->qnode[2];

	double npx = 0;

	if( q0 > thr ) {
	  ++npx;
	  hqpx->Fill( q0 );
	}
	else
	  q0 = 0;

	if( q1 > thr ) {
	  ++npx;
	  hqpx->Fill( q1 );
	}
	else
	  q1 = 0;

	if( q2 > thr ) {
	  ++npx;
	  hqpx->Fill( q2 );
	}
	else
	  q2 = 0;

	double q = q0+q1+q2;
	hq->Fill( q );

	double xcog = ( q0*1.5*px + q1*2.5*px + q2*3.5*px ) / q;
	double dx = xcog - x;
	hdx->Fill( dx );
	dxvsx->Fill( x, dx );

	if( iev == 0 )
	  cout << "x " << x
	       << ": npx " << npx
	       << ", q " << q
	       << ", dx " << dx
	       << endl;

      } // ev

    } // x

    TCanvas c0;
    hqpx->Draw();
    c0.Update();

    TCanvas c1;
    hq->Draw();
    c1.Update();

    TCanvas c2;
    hdx->Draw();
    fitep0( hdx->GetName() );
    TF1 * f1( (TF1*) hdx->GetListOfFunctions()->Last() );
    double sx = f1->GetParameter(1);
    cout << "tilt " << tilt << ", sx " << sx << endl;
    sxvst->Fill( tilt, sx );
    c2.Update();

    TCanvas c3;
    dxvsx->Draw();
    c3.Update();

  } // angles

  TCanvas c4;
  sxvst->SetMarkerStyle(20);
  sxvst->SetMarkerSize(1);
  sxvst->SetMarkerColor(6);
  sxvst->Draw("P");

  fn->Write();

}
