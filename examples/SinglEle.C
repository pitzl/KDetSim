
// define a 3D detector with 5 electrodes

// root -l
// gSystem->Load( "/home/YOU/KDetSim/lib/KDetSim.sl" );
// .x Test3D.C

{
  gStyle->SetCanvasPreferGL(kTRUE);

  // x=100 , y is 50 and thickness 120
  K3D * det = new K3D( 5, 100, 25, 120 ); // units: um

  // define the voltage
  det->Voltage = -40;

  // define the drift mesh size and simulation mesh size in microns

  det->SetUpVolume( 1, 1 );

  // define  columns #, postions, weighting factor 2=0 , material Al=1

  det->SetUpColumn( 0,   0,  0, 3, -99,     2, 1 );
  det->SetUpColumn( 1, 100,  0, 3, -99,     2, 1 );
  det->SetUpColumn( 2,   0, 25, 3, -99,     2, 1 );
  det->SetUpColumn( 3, 100, 25, 3, -99,     2, 1 );
  det->SetUpColumn( 4,  50, 12.5, 3,-99, 16385, 1 );

  Float_t Pos[3]  = { 100, 25, 1 };
  Float_t Size[3] = { 100, 25, 2 };
  det->ElRectangle( Pos, Size, 0, 20 );

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  //define the space charge:

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 2 );//charge density e.c./um^3
  det->NeffF = f2;

  // calculate weigting field
  // calculate electric field
  det->CalField(0);
  det->CalField(1);

  // Show electric field:

  TCanvas c1;
  c1.cd();
  c1.SetTitle( "field" );
  det->Draw( "EFxy", 60 )->Draw("COLZ"); // |E| vs xy at z = 60 um
  
  TCanvas c11;
  c11.cd();
  c11.SetTitle( "EX" );
  det->Draw( "EXxy", 60 )->Draw("COLZ"); // Ex vs xy at z = 60 um

  TCanvas c12;
  c12.cd();
  c12.SetTitle( "EY" );
  det->Draw( "EYxy", 60 )->Draw("COLZ"); // Ey vs xy at z = 60 um


  // Show electric potential:

  TCanvas c2;
  c2.cd();
  c2.SetTitle( "potential" );
  det->Draw( "EPxy", 60 )->Draw("COLZ");

  // set entry points of the track

  det->enp[0] =  30;
  det->enp[1] =  7;
  det->enp[2] =   0;

  det->exp[0] =  30;
  det->exp[1] =  7;
  det->exp[2] = 120;

  // switch on the diffusion

  det->diff = 1;

  // Show mip track

  TCanvas c3;
  c3.SetTitle( "drift" );
  c3.cd();
  det->ShowMipIR(30);

  // calculate induced current

  TCanvas c4;
  c4.SetTitle( "signal" );
  c4.cd();
  det->MipIR(100);
  det->sum->Draw();
  det->neg->Draw("SAME");
  det->pos->Draw("SAME");

}
