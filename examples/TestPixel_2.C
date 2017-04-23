
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x TestPixel_2.C

// study crack between pixels
{
  KPixel * det = new KPixel( 2, 2*60, 60,  60 ); // 2 pixel, [um]

  det->SetUpVolume( 3, 3, 1 ); // um cubes (speed vs granularity)
  //                                                16385 = bit 14 and 1 = readout node at ground
  det->SetUpPixel( 0, 0.5*60,  0.5*60, 20, 30, 1, 16385 );
  det->SetUpPixel( 1, 1.5*60,  0.5*60, 20, 30, 1,     1 ); // 1 = Grnd

  // bias applied at bot, negative to collect e on top
  //det->Voltage = -75; // just full depletion for Neff = 1, d = 300
  //det->Voltage = -200;
  det->Voltage = -15;

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 1 ); // n in n
  //f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  det->CalField(0); // E-field, in KDetector
  det->CalField(1); // Ramo weighting potential

  // 2nd readout node:
  det->SetUpPixel( 0, 0.5*60,  0.5*60, 20, 30, 1,     1 );
  det->SetUpPixel( 1, 1.5*60,  0.5*60, 20, 30, 1, 16385 ); // 2nd readout node

  det->SetUpElectrodes(); // KPixel, sets EG and DM
  det->SetBoundaryConditions(); // KGeometry

  det->CalField(2); // 2nd Ramo weighting potential

  // track:

  det->enp[0] = 0.5*60; // entry point x [um]
  det->enp[1] = 0.5*60; // entry point y
  det->enp[2] =  1; // entry point z

  det->exp[0] = 0.5*60;
  det->exp[1] = 0.5*60;
  det->exp[2] = 59;

  det->Temperature = 293;
  det->diff = 1;

  // show drift paths:

  TCanvas c1;
  det->ShowMipIR(58); // z divisions

  cout
    << "q " << det->sum->Integral()
    << ", " << det->sum->GetSumOfWeights()
    << endl;

  // current vs time:

  TCanvas c2;
  det->MipIR(58);
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
  TProfile qvsx( "qvsx", "Q vs x;x [um];Q", 41, 39.5, 80.5, -9E99, 9E99 );
  TProfile q0vsx( "q0vsx", "Q0 vs x;x [um];Q0", 41, 39.5, 80.5, -9E99, 9E99 );
  TProfile q1vsx( "q1vsx", "Q1 vs x;x [um];Q1", 41, 39.5, 80.5, -9E99, 9E99 );
  det->SetAverage( 9 );
  for( int ix = 60-20; ix <= 60+20; ++ix ) {
    det->enp[0] = ix;
    det->exp[0] = ix;
    det->MipIR(58);
    double q = det->sum->Integral();
    double q0 = det->qnode[0];
    double q1 = det->qnode[1];
    cout << "x " << ix
	 << ": q " << q
	 << ", q0 " << q0
	 << " + q1 " << q1
	 << " = " << q0+q1
	 << endl;
    qvsx.Fill( ix, q );
    q0vsx.Fill( ix, q0 );
    q1vsx.Fill( ix, q1 );
  }
  q0vsx.SetLineColor(2);
  q0vsx.Draw( );
  q1vsx.SetLineColor(8);
  q1vsx.Draw("same");
  qvsx.Draw("same");

}
