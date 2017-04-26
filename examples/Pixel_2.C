
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pixel_2.C

// study charge sharing between 2 pixels
{
  double px = 100; // [um]
  double py =  40; // [um]
  double pz = 300; // [um]
  KPixel * det = new KPixel( 2, 2*px, py, pz ); // 2 pixels, scaled down

  det->SetUpVolume( 2, 2, 2 ); // um cubes (speed vs granularity)
  //                                                16385 = bit 14 and 1
  det->SetUpPixel( 0, 0.5*px, 0.5*py, 0.4*px, 0.4*py, 2, 16385 ); // readout node at ground
  det->SetUpPixel( 1, 1.5*px, 0.5*py, 0.4*px, 0.4*py, 2,     1 ); // 1 = Grnd

  // bias applied at bot, negative to collect e on top
  //det->Voltage = -67; // just full depletion for Neff = 1, d = 300
  det->Voltage = -150; // above full depletion for Neff  =1, d = 300

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  TF3 * f2 = new TF3( "f2", "[0]+x[2]*0", 0, 2*px, 0, py, 0, pz );
  f2->SetParameter( 0, 1 ); // n in n [1/um^3]
  //f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  det->CalField(0); // E-field, in KDetector
  det->CalField(1); // Ramo weighting potential

  // 2nd readout node:
  det->SetPixelW( 0,     1 );
  det->SetPixelW( 1, 16385 ); // 2nd readout node

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  det->CalField(2); // 2nd Ramo weighting potential

  TCanvas cWP;
  det->Draw( "WPxz", 25 )->Draw("COLZ");

  // track:

  double pi = 4*atan(1);
  double wt = 180/pi;
  double bet = 0;
  //double bet = atan( px/pz); // incident angle
  cout << "track angle " << bet*wt << endl;

  det->enp[0] = 0.5*px; // entry point x [um]
  det->enp[1] = 0.5*py; // entry point y
  det->enp[2] =  2; // entry point z

  det->exp[0] = 0.5*px + pz*tan(bet);
  det->exp[1] = 0.5*py;
  det->exp[2] = pz-2;

  det->Temperature = 293;
  det->diff = 1;

  // show drift paths:

  TCanvas c1;
  det->ShowMipIR(37); // z divisions

  cout
    << "q " << det->sum->Integral()
    << ", " << det->sum->GetSumOfWeights()
    << endl;

  // current vs time:

  TCanvas c2;
  det->MipIR(37);
  det->sum->Draw();
  det->pos->Draw("SAME");
  det->neg->Draw("SAME");
  c2.Update();
  cout
    << "q " << det->sum->Integral()
    << ", " << det->sum->GetSumOfWeights()
    << endl;

  // step across gap in x:

  TCanvas c3;
  c3.cd();
  TProfile qvsx( "qvsx", "Q vs x;x [um];Q",     int(px)+1, 0.5*px-0.5, 1.5*px+0.5, -9E9, 9E9 );
  TProfile q0vsx( "q0vsx", "QL vs x;x [um];QL", int(px)+1, 0.5*px-0.5, 1.5*px+0.5, -9E9, 9E9 );
  TProfile q1vsx( "q1vsx", "QR vs x;x [um];QR", int(px)+1, 0.5*px-0.5, 1.5*px+0.5, -9E9, 9E9 );
  det->SetAverage( 9 );
  for( double x = 0.5*px; x <= 1.5*px; x += 1 ) {
    det->enp[0] = x;
    det->exp[0] = x + pz*tan(bet);
    det->MipIR(37);
    double q = det->sum->Integral();
    double q0 = det->qnode[0];
    double q1 = det->qnode[1];
    cout << "x " << x
	 << ": q " << q
	 << ", q0 " << q0
	 << " + q1 " << q1
	 << " = " << q0+q1
	 << endl;
    qvsx.Fill( x, q );
    q0vsx.Fill( x, q0 );
    q1vsx.Fill( x, q1 );
  }
  q0vsx.SetLineColor(2);
  q0vsx.Draw(   );
  q1vsx.SetLineColor(8);
  q1vsx.Draw("same");
  qvsx.Draw("same");

}
