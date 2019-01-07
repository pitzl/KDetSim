
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l s4-100.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  gStyle->SetHistFillStyle(0);

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();
 
  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //               topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "s4-100.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.65, 0.70, 0.95, 0.9 );

  TFile * fdata = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34451.root" );
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "4 #upoint10^{15 }n/cm^{2} edge-on 100 V" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " Data 100 V ", "l" );
  int i0 = colpvsy->FindBin(-0.0575);
  double ph0 = colpvsy->GetBinContent(i0);
  cout << "ph0 " << ph0 << endl;

  TH1F * data = (TH1F*)colpvsy->Clone();

  // smeared sim:

  //TFile * fsim = TFile::Open( "edg3_px50_F4_lam12_e4_h4_V100.root" ); // too narrow
  TFile * fsim = TFile::Open( "edg3_px50_F4_lam16_e4_h4_V100.root" ); // better

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(-0.0575);
  double qh0 = svsy->GetBinContent(k0);
  cout << "qh0 " << qh0 << endl;
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, " KDetSim", "l" );

  lgnd->Draw( "same" );

  c1.Print( "s4-100.ps" );

  double sumdy = 0;
  double chisq = 0;
  int npnt = 0;

  for( int ii = 1; ii <= data->GetNbinsX(); ++ii ) {

    double y = data->GetBinCenter(ii);
    if( y < -0.08 ) continue;
    if( y >  0.03 ) break; // truncated

    double q = data->GetBinContent(ii);
    double e = data->GetBinError(ii);

    int jj = svsy->FindBin(y);
    double s = svsy->GetBinContent(jj);

    cout << ii
	 << "  " << y*1E3
	 << "  " << svsy->GetBinCenter(jj)*1E3
	 << "  " << q
	 << "  " << s
	 << "  " << e
	 << "  " << s-q
	 << endl;

    sumdy += fabs(s-q);
    chisq += pow( s-q, 2 );
    ++npnt;

  }

  cout
    << "mean dy " << sumdy/npnt << endl
    << "mean sy " << sqrt(chisq/npnt) << endl
    << "norm ph " << ph0/qh0
    << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s4-100.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s4-100.ps" );
  ierr = system( "rm -f  s4-100.ps" );
  cout << "evince s4-100.pdf" << endl;

}
