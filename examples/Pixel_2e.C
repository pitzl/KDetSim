
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pixel_2e.C

// draw E fields
{
  KPixel * det = new KPixel( 2, 2*60, 60,  60 ); // 2 pixel, [um]

  det->SetUpVolume( 3, 3, 1 ); // um cubes (speed vs granularity)
  //                                                16385 = bit 14 and 1 = readout node at Grnd
  det->SetUpPixel( 0, 0.5*60,  0.5*60, 20, 30, 1, 16385 );
  det->SetUpPixel( 1, 1.5*60,  0.5*60, 20, 30, 1,     1 );

  // bias applied at bot, negative to collect e on top
  //det->Voltage = -75; // just full depletion for Neff = 1, d = 300
  //det->Voltage = -200;
  det->Voltage = -5;

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 1 ); // n in n [1/um^3]
  //f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  det->CalField(0); // E-field
  det->CalField(1); // Ramo weight field

  TCanvas cGeo;
  det->GetGeom()->Draw("iso");

  TCanvas cEZ;
  det->Draw("EZxz",50)->Draw("COLZ");

  TCanvas cEX;
  det->Draw("EXxz",50)->Draw("COLZ");
}
