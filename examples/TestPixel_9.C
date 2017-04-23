
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x TestPixel_9.C
// .ls

{
  KPixel * det = new KPixel( 9, 3*150, 3*100, 300 ); // 9 pixels, [um]

  Char_t str[100];
  int vbias = 200;
  sprintf( str, "Px9_V%i", vbias );
  det->Read( str, "px9-grd3.root" );

  // track:

  det->enp[0] = 1.5*150; // entry point x [um]
  det->enp[1] = 1.5*100; // entry point y
  det->enp[2] =  3; // entry point z

  det->exp[0] = 1.5*160;
  det->exp[1] = 1.5*100;
  det->exp[2] = 297;

  int ndiv = 49; // 6*49 = 294
  
  det->Temperature = 293;
  det->diff = 1;

  TCanvas c1;
  det->ShowMipIR(ndiv); // divisions
  c1.Update();

  // scan over central pixel:

  TCanvas c3;
  c3.cd();
  TH1I hq( "q", "charge;charge [ke];events", 100, 0, 100 );
  TProfile2D npxvsxy( "npxvsxy", "Npix vs x-y;x [um];y [um];<Npix>",
		      75, 0, 150, 50, 0, 100, -0.5, 99.5 );

  det->SetAverage( 1 );

  for( double x = 151.0; x < 300; x += 2 ) {
    det->enp[0] = x;
    det->exp[0] = x;
    for( double y = 101.0; y < 200; y += 2 ) {
      det->enp[1] = y;
      det->exp[1] = y;

      det->MipIR(ndiv); // one e per div

      int npx = 0;
      double q = 0;
      for( int ipx = 0; ipx < 9; ++ipx ) {
	double qpx = -det->qnode[ipx] * 22.0 / ndiv; // scale to 22 ke
	q += qpx;
	if( qpx > 1.7 ) // threshold [ke]
	  ++npx;
      }
      cout << x
	   << "  " << y
	   << ": " << q
	   << ", " << npx
	   << endl;
      hq.Fill( q );
      npxvsxy.Fill( x-150, y-100, npx );
    }
  }

  npxvsxy.Draw("colz");

}
