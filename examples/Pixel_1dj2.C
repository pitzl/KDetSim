
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pixel_1dj2.C
// .ls

// one pixel, double junction, z scan
{
  //int abias =  40; // |V|, too much double peak, negative field
  int abias =  80; // |V|, double peak
  //int abias = 400; // |V|, merged peaks

  TFile * fn = new TFile( Form( "dj01_v%i.root", abias ), "RECREATE" );

  KPixel * det = new KPixel( 1, 1*150, 1*100, 300 ); // 1 pixel, [um]

  det->SetUpVolume( 5, 5, 2 ); // um cubes (speed vs granularity)

  det->SetUpPixel( 0, 0.5*150,  0.5*100, 65, 40, 1, 16385 ); // collecting electrode at Grnd

  // bias applied at bot, negative to collect e on top
  det->Voltage = -abias; // double junction

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  // doping (space charge): irradiated: quadratic

  TF3 * f3 = new TF3( "f3",
		      "[0]+[2]*pow((x[2]-[1])/300,2)",
		      -9000, 9000, -9000, 9000, 0, 300 );
  f3->SetParameter( 0,  -4 ); // right doping (at pix) [1/um^3]
  f3->SetParameter( 1, 300 ); // x mid
  f3->SetParameter( 2,  12 ); // X**2
  det->NeffF = f3;

  det->taue = 1; // trapping [ns]
  det->tauh = 5;

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  TCanvas cWZ;
  det->Draw1D( "WPxz", 50, 2, 75 )->Draw();
  cWZ.Update();

  TCanvas cEZ;
  //det->Draw( "EZxz", 50 )->Draw("COLZ");
  det->Draw1D( "EZxz", 50, 2, 75 )->Draw();
  cEZ.Update();

  // track:

  det->enp[0] = 0.0*150; // left
  det->enp[1] = 0.5*100; // entry point y
  det->enp[2] = 0.5*300; // height

  det->exp[0] = 1.0*150; // right
  det->exp[1] = 0.5*100;
  det->exp[2] = 0.5*300; // height

  det->Temperature = 293;
  det->diff = 1;

  //TCanvas c1; det->ShowMipIR(30); // divisions, without trapping

  // pulse shape:
  /*
  TCanvas c2;
  det->MipIR(30);
  det->sum->Draw();
  det->pos->Draw("SAME");
  det->neg->Draw("SAME");
  c2.Update();
  cout
    << "q " << det->sum->Integral()
    << ", " << det->sum->GetSumOfWeights()
    << endl;
  */
  // z scan:

  TCanvas c3;
  c3.cd();
  TProfile qvsz( "qvsz", "Q vs z;z [um];Q", 50, 0, 300, -9E99, 9E99 );
  det->SetAverage( 9 );
  for( int iz = 3; iz < 300; iz += 6 ) {
    det->enp[2] = iz;
    det->exp[2] = iz;
    det->MipIR(30);
    double q = -det->sum->Integral(); // with trapping [e]
    cout << "z " << iz
	 << ": q " << q
	 << endl;
    qvsz.Fill( iz, q );
  }
  qvsz.SetLineColor(2);
  qvsz.Draw( );

  fn->Write();

}
