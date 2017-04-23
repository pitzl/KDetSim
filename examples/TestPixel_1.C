
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x TestPixel_8.C
// .ls

{
  KPixel * det = new KPixel( 1, 1*150, 1*100, 300 ); // 1 pixel, [um]

  det->SetUpVolume( 5, 5, 5 ); // 5x5x5 um cubes (much faster)
  det->SetUpPixel( 0, 0.5*150,  0.5*100, 65, 40, 1, 16385 ); // collecting electrode at Grnd

  // bias applied at bot, negative to collect e on top
  //det->Voltage = -50; // below full depletion for Neff = 1, d = 300
  //det->Voltage = -70; // below full depletion for Neff = 1
  //det->Voltage = -75; // just full depletion for Neff = 1, d = 300
  //det->Voltage = -80; // full depletion for Neff = 1
  det->Voltage = -99; // above full depletion for Neff = 1

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  // doping (space charge):

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 1 ); // n in n
  //f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  det->enp[0] = 0.5*150; // entry point x [um]
  det->enp[1] = 0.5*100; // entry point y
  det->enp[2] =  0; // entry point z

  det->exp[0] = 0.5*150;
  det->exp[1] = 0.5*100;
  det->exp[2] = 300;

  det->Temperature = 293;
  det->diff = 1;

  TCanvas c1;
  det->ShowMipIR(30); // divisions

  TCanvas c2;
  det->MipIR(30);
  det->sum->Draw();
  det->pos->Draw("SAME");
  det->neg->Draw("SAME");

  cout
    << "q " << det->sum->Integral()
    << ", " << det->sum->GetSumOfWeights()
    << endl;

}
