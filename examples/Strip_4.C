{
  // Different mesh structures Y[]= denotes the depth in um 
  // and correspoding step in X[]. The depts follow from rigth to left and they shoul
  // obviousy sum to detector thickness.

  //   Float_t Y[]={90,75,70,30,35};
  //   Float_t X[]={3,2.5,2,1.5,1};

  //   Float_t Y[]={150,150}; 
  //   Float_t X[]={2,1};

  //   Float_t Y[]={60,55,50,45,30,35,15,10};
  //   Float_t X[]={3,2.75,2.5,2.25,2,1.75,1.5,1.25};

  Float_t Y[]={100,75,  69.5,37.5,18};
  Float_t X[]={  2, 1,0.50,0.375,0.25};
  Int_t i,j;

  KMesh Mesh(300);
  Float_t Bins[1000];
 
  // Generation of bins for Y direction
  Int_t N=Mesh.GetBins(5,Y,X,Bins);

  // Generation of bins for X direction 
  Float_t *Xbins=new Float_t[321];
  Xbins[0]=0;
  for(Int_t i=1;i<=140;i++) Xbins[i]= Xbins[i-1]+2;;

  // Uniform space charge of -4e12 cm-3
  TF3 *f3=new TF3("f2","x*y*z*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,0);
 
  // Define strip detector

  KStrip *det=new KStrip(80,40,1,1,300);
  det->Voltage=-500;

  det->SetUpVolume(N,Bins,40,Xbins);
  det->SetUpElectrodes();

 for(j=1;j<=det->ny;j++)
      for(i=1;i<=det->nx;i++)
	if(j>=70) det->DM->SetBinContent(i,j,1,20); else det->DM->SetBinContent(i,j,1,1);

  det->SetBoundaryConditions();


  det->NeffF=f3;

  det->SetCalculationParameters(1e-6,25000);
  det.CalField(0);
  det->SetCalculationParameters(1e-6,2000);
  det.CalField(1);
}
