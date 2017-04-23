
{
  gStyle->SetCanvasPreferGL(kTRUE);
  Int_t i,j,k;
  KPixel *det=new KPixel(1,50,50,100);
  det->Voltage=200;

  det->SetUpVolume(1,1,1);
  det->SetUpPixel(0,25,25,25,25,1,16385);
  det->SetUpElectrodes();


  for(k=1;k<=det->nz;k++)
    for(j=1;j<=det->ny;j++)
      for(i=1;i<=det->nx;i++)
	if(k>=30) det->DM->SetBinContent(i,j,k,1); else det->DM->SetBinContent(i,j,k,1);

  det->SetBoundaryConditions();

  TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,0);
  det->NeffF=f2;

  det->CalField(0);
  //   det->CalField(1);

  det->enp[0]=100;
  det->enp[1]=100;
  det->enp[2]=1;
  det->exp[0]=100;
  det->exp[1]=100;
  det->exp[2]=100;
  
  //   det->Save("DD3D","tempp.root");
   
}
