{
  gStyle->SetCanvasPreferGL(kTRUE);

  K3D *det=new K3D(1,200,200,100);
  det->Voltage=200;
  det->SetUpVolume(2,2);
  
  Float_t Pos[3]={100,100,98};
  Float_t Size[3]={100,100,2};
  det->ElRectangle(Pos,Size,0,20);
  
  det->SetUpColumn( 0, 100, 100, 20, -40, 16385, 1 ); // readout node at Grnd

  det->SetUpElectrodes(1);
  det->SetBoundaryConditions();

  TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,-2);
  det->NeffF=f2;

  det->CalField(0); // E field
  det->CalField(1); // Ramo weighting potential

  det->enp[0]=50;
  det->enp[1]=50;
  det->enp[2]=1;
  det->exp[0]=50;
  det->exp[1]=50;
  det->exp[2]=100;

  det->diff=1;
  
  //   det->Save("DD3D","tempp.root");
  TCanvas c1;
  c1.cd();
  det->ShowMipIR(30);

  TCanvas c2;
  c2.cd();
  det->Draw("EPxy",99)->Draw("COLZ");

  TCanvas c3;
  c3.cd();
  det.MipIR(100);
  det->sum.Draw();
  det->neg.Draw("SAME");
  det->pos.Draw("SAME");
   
}
