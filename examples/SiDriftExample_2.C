{
  KDetector det;
  det.SStep=3;
  det.nx=3500;
  det.ny=100;
  det.nz=1;
  Float_t dimX=35000;
  Float_t dimY=300;
  //init geometry
  det.EG=new TH3I("EG","EG",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,1);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");
  //init material
  det.DM=new TH3I("DM","DM",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,1);
  det.DM->GetXaxis()->SetTitle("x [#mum]");
  det.DM->GetYaxis()->SetTitle("y [#mum]");
  det.DM->GetZaxis()->SetTitle("z [#mum]");
  //init space charge histo
  det.NeffH=new TH3F("Neff","Neff",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,1);
  det.NeffH->GetXaxis()->SetTitle("x [#mum]");
  det.NeffH->GetYaxis()->SetTitle("y [#mum]");
  det.NeffH->GetZaxis()->SetTitle("z [#mum]");

  // Collection electrode
  Float_t BackPos[3]={34930,299,0.5};
  Float_t BackSiz[3]={50,0.5,0.5};
  det.ElRectangle(BackPos,BackSiz,16385,0);  
 
  // Strip definition
  Float_t StripPos[3]={50,300,0};
  Float_t StripSiz[3]={30,0.5,0};

  // definition of votlages to be applied to the field strips

  det->Voltages->Set(291);
  for(int i=0; i<=290; i++)
    {
      StripPos[0]=i*120; 
     det->Voltages[i]=-1800+i*6;
     StripPos[1]=299;
     det.ElRectangle(StripPos,StripSiz,det->SetElecVolt(i),0);  
    
     StripPos[1]=1;
     det.ElRectangle(StripPos,StripSiz,det->SetElecVolt(i),0);  
    }

  //  SetUpMaterial
  for(int j=0;j<=det.ny;j++)
    for(int i=0;i<=det.nx;i++) 
      {
          det->DM->SetBinContent(i,j,1,0); 
          det->NeffH->SetBinContent(i,j,1,0.5);  // very high resistivity Neff=1e11 cm-3
      }

  det->SetBoundaryConditions();

    det->CalField(0);
    det->CalField(1);

    det->SetEntryPoint(27000,300,0.5);
    det->SetExitPoint(27000,1,0.5);

    det.SetDriftHisto(2000e-9);
    det.SetPrecision(1e-9);
    det->diff=1;
    det->ShowMipIR(300);
}
