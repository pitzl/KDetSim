
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .L fitep0.C
// TMemStat mm("gnubuiltin"); // memstat.root
// .x Pixel_t.C

// .x root/tutorials/memstat/memstatExample.C ?

// tilt

{
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

  //double tilt = 18.3/wt; // pos slope dxvsx
  double tilt = 19.3/wt; // flat, optimum, at geometric minimum
  //double tilt = 20.3/wt; // neg slope dxvsx

  int ndiv = pz/cos(tilt)/5+1; // steps

  cout << "track angle " << tilt*wt << endl;
  cout << "track divisions " << ndiv << endl;

  det->enp[0] = 2.5*px - 0.5*pz*tan(tilt);; // entry point x [um]
  det->enp[1] = 0.5*py; // entry point y
  det->enp[2] =  2; // entry point z

  det->exp[0] = 2.5*px + 0.5*pz*tan(tilt);
  det->exp[1] = 0.5*py;
  det->exp[2] = pz-2;

  det->Temperature = 293; // [K]
  det->diff = 1; // 0: still mem leak
  det->SStep = 5; // [um] drift step [1: use 42.8% mem, sys: 9.5%]
  det->Landau = 1;

  double thr = 1.5; // [ke]

  TCanvas cM;
  det->ShowMipIR(ndiv); // divisions
  cM.Update();

  // step across one pixel in x:

  int nev = 16;

  gStyle->SetOptStat(111111);
  TH1I * hqpx = new TH1I( "qpx", "pixel charge;pixel charge [ke];pixel",  100, 0, 50 );
  TH1I * hq = new TH1I( "q", "charge;charge [ke];tracks",  100, 0, 100 );
  TH1I * hdx = new TH1I( "dx", "residual;#Deltax [um];tracks",  100, -50, 50 );
  TProfile * dxvsx = new
    TProfile( "dxvsx", "#Deltax vs x;x [um];<#Deltax [um]", 101, 1.5*px-0.5, 2.5*px+0.5,
	      -999, 999 );

  for( double x = 1.5*px; x <= 2.5*px; x += 1 ) {

    det->enp[0] = x - 0.5*pz*tan(tilt);
    det->exp[0] = x + 0.5*pz*tan(tilt);

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

  TCanvas c1;
  hq->Draw();

  TCanvas c2;
  hdx->Draw();
  fitep0( hdx->GetName() );
  TF1 * f1( (TF1*) hdx->GetListOfFunctions()->Last() );
  double rx = f1->GetParameter(1);
  cout << rx << endl;

  TCanvas c3;
  dxvsx->Draw();

}
