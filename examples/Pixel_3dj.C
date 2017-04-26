
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pixel_8.C
// .ls

// shoot quer through 3 pixels, double junction
{
  KPixel * det = new KPixel( 3, 3*100, 1*150, 300 ); // 3 pixel, [um]

  det->SetUpVolume( 5, 5, 5 ); // 5x5x5 um cubes (speed vs granularity)

  det->SetUpPixel( 0, 0.5*100,  0.5*150, 40, 65, 1, 1     ); // Grnd
  det->SetUpPixel( 1, 1.5*100,  0.5*150, 40, 65, 1, 16385 ); // collecting electrode at Grnd
  det->SetUpPixel( 2, 2.5*100,  0.5*150, 40, 65, 1, 1     ); // Grnd

  // bias applied at bot, negative to collect e on top
  det->Voltage = -20; // zero mid
  //det->Voltage = -30; // dip mid
  //det->Voltage = -50; // raised mid
  //det->Voltage = -99; // raised mid

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  // doping (space charge):

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]+[1]*x[2]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0,  1 ); // start [1/um^3]
  f2->SetParameter( 1, -2.0/300 ); // slope for double junction
  det->NeffF = f2;

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  // track:

  det->enp[0] = 0.0*100; // left
  det->enp[1] = 0.5*150; // entry point y
  det->enp[2] = 0.0*300; // entry point z

  det->exp[0] = 3.0*100;
  det->exp[1] = 0.5*150;
  det->exp[2] = 1.0*300;

  det->Temperature = 293;
  det->diff = 1;

  TCanvas c1;
  det->ShowMipIR(30); // divisions

  TCanvas cEZ;
  det->Draw( "EZxz", 50 )->Draw("COLZ");

  TCanvas c2;
  det->MipIR(30);
  det->sum->Draw();
  det->pos->Draw("SAME");
  det->neg->Draw("SAME");

  cout
    << "q " << det->sum->Integral()
    << endl;

}
