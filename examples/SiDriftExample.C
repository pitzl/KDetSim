{
  KDetector det;
  det.SStep=3;
  det.nx=1000;
  det.ny=100;
  det.nz=1;
  //init geometry
  det.EG=new TH3I("EG","EG",det.nx,0,10000,det.ny,0,200,det.nz,0,1);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");
  //init material
  det.DM=new TH3I("DM","DM",det.nx,0,10000,det.ny,0,200,det.nz,0,1);
  det.DM->GetXaxis()->SetTitle("x [#mum]");
  det.DM->GetYaxis()->SetTitle("y [#mum]");
  det.DM->GetZaxis()->SetTitle("z [#mum]");
  //init space charge histo
  det.NeffH=new TH3F("Neff","Neff",det.nx,0,10000,det.ny,0,200,det.nz,0,1);
  det.NeffH->GetXaxis()->SetTitle("x [#mum]");
  det.NeffH->GetYaxis()->SetTitle("y [#mum]");
  det.NeffH->GetZaxis()->SetTitle("z [#mum]");

  // Collection electrode
  Float_t BackPos[3]={9700,199,0.5};
  Float_t BackSiz[3]={5,0.5,0.5};
  det.ElRectangle(BackPos,BackSiz,16385,0);  
 
  // Strip definition
  Float_t StripPos[3]={50,200,0};
  Float_t StripSiz[3]={250,0.5,0};

  // definition of votlages to be applied to the field strips
  for(int i=0; i<10; i++) det.Voltages[i]=-10*(9-i)-10;
  for(int i=0; i<10; i++)
    {
    StripPos[0]=i*1000; 
 
     StripPos[1]=199;
     det.ElRectangle(StripPos,StripSiz,(1<<(15+i)),0);  
    
     StripPos[1]=1;
     det.ElRectangle(StripPos,StripSiz,(1<<(15+i)),0);  
    }

  //  SetUpMaterial
  for(int j=0;j<=det.ny;j++)
    for(int i=0;i<=det.nx;i++) 
          det->DM->SetBinContent(i,j,1,0); 
	  det->NeffH->SetBinContent(i,j,1,0.1);


  det.Voltage=-1000;
  det.Voltage2=0;


  det->SetBoundaryConditions();

    det->CalField(0);
    det->CalField(1);

    det->SetEntryPoint(7000,200,0.5);
    det->SetExitPoint(7000,1,0.5);

    det.SetDriftHisto(10000e-9);
    det.SetPrecision(1e-9);
    det->diff=0;
    det->ShowMipIR(200);
}
