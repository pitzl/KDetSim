
{
  gStyle->SetCanvasPreferGL(kTRUE);

  KDetector det;
  det.nx=100;
  det.ny=88;
  det.nz=50;
  det.EG=new TH3I("EG","EG",det.nx,0,200,det.ny,0,176,det.nz,0,100);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");
  det->GetGrid(det.EG,1);
  
  det.Voltage=100;
  TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,-2);
  det.NeffF=f2;
  
  Float_t W[3]={3,3,50};
  
  Float_t P1[3]={53,3,50};
  Float_t P2[3]={147,3,50};

  Float_t P3[3]={3,86,50};
  Float_t P4[3]={197,86,50};

  Float_t P5[3]={53,172,50};
  Float_t P6[3]={147,172,50};

  Float_t P7[3]={100,86,50};
  Float_t P8[3]={100,86,50};

  det.ElLine(P1,P2,W,2,1);
  det.ElLine(P2,P4,W,2,1);
  det.ElLine(P1,P3,W,2,1);
  det.ElLine(P3,P5,W,2,1);
  det.ElLine(P4,P6,W,2,1);
  det.ElLine(P5,P6,W,2,1);

  det->ElCylinder(P7,10,50,3,16385,2);
  det->ElCylinder(P7,8,50,3,16385,1);

  det->SetBoundaryConditions();

  det.CalField(0);
  det.CalField(1);

  det->enp[0]=40;
  det->enp[1]=100;
  det->enp[2]=1;
  det->exp[0]=40;
  det->exp[1]=90;
  det->exp[2]=100;

}
