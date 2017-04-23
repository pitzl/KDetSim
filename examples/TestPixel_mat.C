{
  Int_t i,j,k;
  KPixel * det = new KPixel( 1, 120, 120, 100 ); // one pixel

  det->Voltage = -200;

  det->SetUpVolume( 2, 2, 2 );
  det->SetUpPixel( 0, 60, 60, 20, 20, 2,16385 );

  det->SetUpElectrodes();
  det->SetBoundaryConditions();

  TF3 * f2 = new TF3( "f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000 );
  f2->SetParameter( 0, 0 ); // no doping
  det->NeffF = f2;

  for( k=1; k<=det->nz; k++ )
    for( j=1; j<=det->ny; j++ )
      for( i=1; i<=det->nx; i++ )
      	if( (j<10 || j>50 || i<10 || i>50 ) && k>45 )
	  det->DM->SetBinContent( i, j, k, 2 ); // material
	else
	  det->DM->SetBinContent( i, j, k, 0 );

  det->CalField(0);
  det->CalField(1);
   
  det->enp[0] = 60; // entry point x
  det->enp[1] = 60; // entry point y
  det->enp[2] =  1; // entry point z

  det->exp[0] = 60;
  det->exp[1] = 60;
  det->exp[2] = 99;

}
