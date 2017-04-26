
// root -l
// .x Lan.C

// Landau additive? no

{
  // calculate most probable charge per segement in ke:

  double LQ = 22; // [ke]

  int nseg = 22;

  TF1 * lan = new TF1( "lan"," TMath::Landau(x,[0],[1])", 0, 20*LQ/nseg );
  lan->SetParameter( 0, LQ/nseg ); // peak
  lan->SetParameter( 1, LQ/nseg/15 ); // width

  TF1 * lan2 = new TF1( "lan2"," TMath::Landau(x,[0],[1])", 0, 20*LQ );
  lan2->SetParameter( 0, LQ ); // peak
  lan2->SetParameter( 1, LQ/15 ); // width

  TH1I hqs( "hqs", "segment charge;charge [ke];segments", 100, 0, 10*LQ/nseg );
  TH1I hqt( "hqt", "track charge;track charge [ke];events", 100, 0, 100 );
  TH1I hq2( "hq2", "track charge;track charge [ke];events", 100, 0, 100 );

  for( int j = 0; j < 10000; ++j ) {

    double Q = 0;
    for( int i = 0; i < nseg; ++i ) {
      double q = lan->GetRandom();
      hqs.Fill( q );
      Q += q;
    }

    hqt.Fill( Q ); // summed
    hq2.Fill( lan2->GetRandom() ); // full

  } // events

  TCanvas c1;
  hqs.Draw();

  TCanvas c2;
  hqt.Draw();

  TCanvas c3;
  hq2.Draw();

}
