
{
  gStyle->SetCanvasPreferGL(kTRUE);
  KDetector det;
  det.nx=120;
  det.ny=100;
  det.nz=50;
  det.EG=new TH3I("EG","EG",det.nx,0,240,det.ny,0,500,det.nz,0,100);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");

  det.DM=new TH3I("DM","DM",det.nx,0,240,det.ny,0,500,det.nz,0,100);
  det.DM->GetXaxis()->SetTitle("x [#mum]");
  det.DM->GetYaxis()->SetTitle("y [#mum]");
  det.DM->GetZaxis()->SetTitle("z [#mum]");
  
  det.Voltage=500;
  TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,-2);
  det.NeffF=f2;

  //BackPlane
  Float_t BackPos[3]={120,250,0.5};
  Float_t BackSiz[3]={119.9,249.9,0.1};
  det.ElRectangle(BackPos,BackSiz,2,0);

  //Strips
  Float_t BackPos[3]={40,350,99.5};
  Float_t BackSiz[3]={9.9,149.,0.1}; 
  for(Int_t i=0;i<3;i++)
    { 
      BackPos[0]=i*80+40;
      if(i==1) 
	det.ElRectangle(BackPos,BackSiz,16385,0); 
      else  
	det.ElRectangle(BackPos,BackSiz,1,0);
    }
  
  //Guard
  Float_t BackPos[3]={120,100,99.5};
  Float_t BackSiz[3]={119.1,30,0.1};
  det.ElRectangle(BackPos,BackSiz,1,0);

  det->SetBoundaryConditions();
  det->CalField(0);
  det->CalField(1);

  det->SetEntryPoint(110,250,1);
  det->SetExitPoint(120,250,49);
  det->ShowMipIR(20);

}
