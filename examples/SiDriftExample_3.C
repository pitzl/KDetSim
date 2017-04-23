{
  KDetector det;
  det.SStep=3;
  det.nx=100;
  det.ny=100;
  det.nz=150;
  Float_t dimX=1000;
  Float_t dimY=300;
  Float_t dimZ=900;
  //init geometry
  det.EG=new TH3I("EG","EG",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,dimZ);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");
  //init material
  det.DM=new TH3I("DM","DM",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,dimZ);
  det.DM->GetXaxis()->SetTitle("x [#mum]");
  det.DM->GetYaxis()->SetTitle("y [#mum]");
  det.DM->GetZaxis()->SetTitle("z [#mum]");
  //init space charge histo
  det.NeffH=new TH3F("Neff","Neff",det.nx,0,dimX,det.ny,0,dimY,det.nz,0,dimZ);
  det.NeffH->GetXaxis()->SetTitle("x [#mum]");
  det.NeffH->GetYaxis()->SetTitle("y [#mum]");
  det.NeffH->GetZaxis()->SetTitle("z [#mum]");

  // Collection electrodes
  Float_t BackPos[3]={950,299,450};
  Float_t BackSiz[3]={30,0.5,100};
  for(int el=0;el<3;el++) {
    BackPos[2]=el*300+150;; 
    if(el==1)
      det.ElRectangle(BackPos,BackSiz,16385,0); // readout node
    else
      det.ElRectangle(BackPos,BackSiz,1,0);  
  }
 
  // Field Strips - definitions
  Float_t StripPos[3]={50,300,dimZ/2};
  Float_t StripSiz[3]={30,0.5,dimZ/2};

  // definition of votlages to be applied to the field strips

  det->Voltages->Set(291);

  for(int i=0; i<=7; i++) {
    StripPos[0]=i*120; 
    det->Voltages[i]=-70+i*6;
    StripPos[1]=299;
    det.ElRectangle(StripPos,StripSiz,det->SetElecVolt(i),0);  

    StripPos[1]=1;
    det.ElRectangle(StripPos,StripSiz,det->SetElecVolt(i),0);  
  }

  //  SetUpMaterial
  for(int k=0;k<=det.nz;k++) 
    for(int j=0;j<=det.ny;j++)
      for(int i=0;i<=det.nx;i++) {
	det->DM->SetBinContent(i,j,k,0); 
	det->NeffH->SetBinContent(i,j,k,0.5);  // very high resistivity Neff=1e11 cm-3
      }

  det->SetBoundaryConditions();

  det->CalField(0);
  det->CalField(1);

  det->SetEntryPoint(200,300,450);
  det->SetExitPoint(200,1,450);

  det.SetDriftHisto(500e-9);
  det.SetPrecision(1e-9);
  det->diff=1;
  det->ShowMipIR(300);

}
