
// root -l
// gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
// .x Pixelee_9.C
// .ls

// calculate E-field for 3x3 pixels (slow)
{
  KPixel * det = new KPixel( 9, 3*150, 3*100, 300 ); // 9 pixels, [um]

  det->SetUpVolume( 5, 5, 5 ); // um cubes (speed vs granularity)
  //det->SetUpVolume( 3, 3, 3 ); // um cubes (speed vs granularity)

  //Bit 1 = 1 -> GND - 0 V bias
  //Bit 2 = 2 -> Voltage (usual bias voltage)
  //Bit 15-32 = 32768 -> Additional Voltages 

  //Bits determining boundary conditions:
  //bit 2 = 4  -> down val
  //bit 3 = 8  -> up val 
  //bit 4 = 16 -> left val
  //bit 5 = 32 -> rigth val
  //bit 6 = 64  -> down der
  //bit 7 = 128  -> up der 
  //bit 8 = 256 -> left der
  //bit 9 = 512 -> rigth der
  //the 3D section is separated
  //bit 10= 1024 -> out val
  //bit 11= 2048 -> in val
  //bit 12= 4096 -> out der
  //bit 13= 8192 -> in der
  //bit 14= 16384-> read out node

  det->SetUpPixel( 0, 0.5*150,  0.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 1, 1.5*150,  0.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 2, 2.5*150,  0.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 3, 0.5*150,  1.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 4, 1.5*150,  1.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 5, 2.5*150,  1.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 6, 0.5*150,  2.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 7, 1.5*150,  2.5*100, 65, 40, 1, 1 );
  det->SetUpPixel( 8, 2.5*150,  2.5*100, 65, 40, 1, 1 );

  int abias = 200;
  det->Voltage = -abias;

  det->SetUpElectrodes(); // in KPixel
  det->SetBoundaryConditions(); // in KGeometry

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 1 ); // n-in-n doping [1/um^3]
  //f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  det->CalField(0); // E-field, in KDetector

  // Ramo weighting field, one pixel at a time:

  for( int ipx = 0; ipx < 9; ++ipx ) {
    det->SetPixelW( ipx, 16385 ); // bits 14 and 1: readout node at Grdn
    det->SetUpElectrodes(); // in KPixel
    det->SetBoundaryConditions();
    det->CalField(1+ipx); // Ramo weighting potential
    det->SetPixelW( ipx, 1 ); // back to Grnd
  }

  Char_t str[100];
  sprintf( str, "Px9_V%i", abias );
  det->Save( str, "px9.root" );

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

}
