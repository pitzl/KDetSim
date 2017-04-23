
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x TestPixel_8.C
// .ls

// shoot quer through 3 pixels
{
  KPixel * det = new KPixel( 3, 3*100, 1*150, 300 ); // 3 pixel, [um]

  det->SetUpVolume( 5, 5, 5 ); // 5x5x5 um cubes (speed vs granularity)

  det->SetUpPixel( 0, 0.5*100,  0.5*150, 40, 65, 1, 1     ); // Grnd
  det->SetUpPixel( 1, 1.5*100,  0.5*150, 40, 65, 1, 16385 ); // collecting electrode at Grnd
  det->SetUpPixel( 2, 2.5*100,  0.5*150, 40, 65, 1, 1     ); // Grnd

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

}
