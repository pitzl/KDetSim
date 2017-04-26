
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .L fitlang.C
// .x Pixel_1l.C
// .ls
// one pixel, with Landau flucts

{
  KPixel * det = new KPixel( 1, 1*150, 1*100, 300 ); // 1 pixel, [um]

  det->Voltage = -150; // applied at back, collect e on top

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 1 ); // n-in-n doping [1/um^3]
  det->NeffF = f2;

  det->SetUpVolume( 5, 5, 3 ); // micron cubes
  det->SetUpPixel( 0, 0.5*150,  0.5*100, 65, 40, 1, 16385 ); // collecting electrode

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  // track entry and exit point:

  det->Temperature = 293;
  det->diff = 1;
  det->Landau = 0;

  det->enp[0] = 0.5*150; // entry point x [um]
  det->enp[1] = 0.5*100; // entry point y
  det->enp[2] =  3; // entry point z

  det->exp[0] = 0.5*150;
  det->exp[1] = 0.5*100;
  det->exp[2] = 297; // [um]

  double dist = TMath::Sqrt( TMath::Power( det->enp[0] - det->exp[0], 2 ) +
			     TMath::Power( det->enp[1] - det->exp[1], 2 ) +
			     TMath::Power( det->enp[2] - det->exp[2], 2 ) );

  // Landau path: 294 um

  const int nseg = 21;

  // segment vector:
  Double_t vs[3] = { ( (double) det->exp[0]-det->enp[0] ) / nseg,
		     ( (double) det->exp[1]-det->enp[1] ) / nseg,
		     ( (double) det->exp[2]-det->enp[2] ) / nseg };

  std::cout << "segment vector: " << vs[0] << ", " << vs[1] << ", " << vs[2] << std::endl;

  // track segments:
  Double_t * entx = new Double_t [nseg];
  Double_t * extx = new Double_t [nseg];
  Double_t * enty = new Double_t [nseg];
  Double_t * exty = new Double_t [nseg];
  Double_t * entz = new Double_t [nseg];
  Double_t * extz = new Double_t [nseg];

  TH1F * hisp[nseg];
  TH1F * hisn[nseg];
  TH1F * hiss[nseg]; 

  std::cout << "segments:" << std::endl;

  for( int i = 0; i < nseg; ++i ) {

    entx[i] = vs[0]*i + det->enp[0];
    enty[i] = vs[1]*i + det->enp[1];
    entz[i] = vs[2]*i + det->enp[2];
    extx[i] = vs[0]*(i+1) + det->enp[0];
    exty[i] = vs[1]*(i+1) + det->enp[1];
    extz[i] = vs[2]*(i+1) + det->enp[2];
    hisp[i] = new TH1F(); 
    hisn[i] = new TH1F();
    hiss[i] = new TH1F();  
    std::cout << i
	      << ":  " << entx[i]
	      << " - " << extx[i]
	      << ",  " << enty[i]
	      << " - " << exty[i]
	      << ",  " << entz[i]
	      << " - " << extz[i]
	      << std::endl;

  } // seg

  // induced currents per segment without fluctuations:

  std::cout << "segment currents:" << std::endl;

  for( int i = 0; i < nseg; ++i ) { // must be new loop: det->enp is overwritten

    std::cout << i;

    det->sum->Reset();
    det->pos->Reset();
    det->neg->Reset();

    det->SetEntryPoint( entx[i], enty[i], entz[i] );
    det->SetExitPoint( extx[i], exty[i], extz[i] ); // bug fix DP

    det->MipIR( dist / nseg ); // track segment [um]

    det->sum->Copy( *hiss[i] );
    det->pos->Copy( *hisp[i] );
    det->neg->Copy( *hisn[i] );

    std::cout << ": " << hiss[i]->Integral()
	      << std::endl;

  } // seg

  // event loop:

  // calculate most probable charge per segement in ke:
  double lanseg = dist/nseg*75E-3; // 0.075 ke/um

  // energy loss distribution:
  TF1 * lan = new TF1( "lan"," TMath::Landau(x,[0],[1])", 0, 20*lanseg );
  lan->SetParameter( 0, lanseg ); // peak
  lan->SetParameter( 1, lanseg/8. ); // width

  TCanvas c1;

  const int pev = 10;

  TH1F * hev[pev];

  TH1I hqs( "hqs", "segment charge;charge [ke];segments", 100, 0, 10*lanseg );
  TH1I hqt( "hqt", "track charge;track charge [ke];events", 100, 0, 100 );

  bool ldb = 0;

  for( int j = 0; j < 10000; ++j ) {

    if( j < pev )
      hev[j] = new TH1F( Form( "landau_%i", j ), "current vs time;time [ns];current",
			 det->pos->GetNbinsX(),
			 det->pos->GetXaxis()->GetXmin(),
			 det->pos->GetXaxis()->GetXmax() );

    // apply Landau randoms:

    if( ldb ) std::cout << j << ":";

    TH1F hp( "pulse", "signal pulse",
	     det->pos->GetNbinsX(),
	     det->pos->GetXaxis()->GetXmin(),
	     det->pos->GetXaxis()->GetXmax() );

    for( int i = 0; i < nseg; ++i ) {
      double q = lan->GetRandom();
      hqs.Fill( q );
      if( ldb ) std::cout << " " << q;
      if( j < pev )
	hev[j]->Add( hiss[i], q/nseg );
      hp.Add( hiss[i], q/nseg );
    }

    // pulse vs time for m events:

    if( j == 0 )
	hev[j]->Draw();
    else if( j < pev ) {
      hev[j]->SetLineColor(j+1);
      hev[j]->Draw("SAME");
    }

    if( ldb )
      std::cout << " = " << hp.GetEntries()
		<< ", " << hp.GetSumOfWeights()
		<< ", " << hp.Integral()
		<< std::endl;

    hqt.Fill( hp.Integral() );

  } // events

  TCanvas c2;
  hqs.Draw();

  TCanvas c3;
  hqt.Draw();

  fitlang( hqt.GetName() );
}
