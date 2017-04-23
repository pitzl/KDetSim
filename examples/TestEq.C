
{
  gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );

  Int_t i=30,j=210,k=1;
  TH3F *test=det.Real.U;
  Float_t Xl=test->GetBinContent(i-1,j,k);
  Float_t Xr=test->GetBinContent(i+1,j,k);
  Float_t Xc=test->GetBinContent(i,j,k);
  Float_t bXl=test->GetXaxis()->GetBinUpEdge(i)-test->GetXaxis()->GetBinUpEdge(i-1);
  Float_t bXr=test->GetXaxis()->GetBinUpEdge(i+1)-test->GetXaxis()->GetBinUpEdge(i);

  Float_t Yl=test->GetBinContent(i,j-1,k);
  Float_t Yr=test->GetBinContent(i,j+1,k);
  Float_t Yc=test->GetBinContent(i,j,k);
  Float_t bYl=test->GetYaxis()->GetBinUpEdge(j)-test->GetYaxis()->GetBinUpEdge(j-1);
  Float_t bYr=test->GetYaxis()->GetBinUpEdge(j+1)-test->GetYaxis()->GetBinUpEdge(j);

  Float_t Zl=test->GetBinContent(i,j,k-1);
  Float_t Zr=test->GetBinContent(i,j,k+1);
  Float_t Zc=test->GetBinContent(i,j,k);
  Float_t bZl=test->GetZaxis()->GetBinUpEdge(k)-test->GetZaxis()->GetBinUpEdge(k-1);
  Float_t bZr=test->GetZaxis()->GetBinUpEdge(k+1)-test->GetZaxis()->GetBinUpEdge(k);

  
  printf("%f \n",(Xr-2*Xc+Xl)/(bXl*bXl)+(Yr-2*Yc+Yl)/(bYl*bYl));
  
}
