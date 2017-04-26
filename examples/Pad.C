
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pad.C
// .ls

{
  TF1 *neff = new TF1( "neff", "[0]+x[0]*0", 0, 1000 );
  neff->SetParameter( 0, 1 ); // n-in-n doping [1/um^3]

  KPad det( 50, 300 ); // simple diode
  det.Neff = neff;
  det.Voltage = -200; // applied at back, collect at top

  det.SetUpVolume(1);
  det.SetUpElectrodes();

  TCanvas c1;
  det.SetEntryPoint( 25, 299.9, 0.5 );
  det.SetExitPoint( 25, 1., 0.5);
  det.Temperature = 293;
  det.diff = 1;
  det.ShowMipIR(200);

  TCanvas c2;
  det.MipIR(200);
  det.sum->Draw();
  det.pos->Draw("SAME");
  det.neg->Draw("SAME");

}
