
// Daniel Pitzl, Jan 2019
// KDetSim edg

// root -l s2-325.C

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

  c1.Print( "s2-325.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.20, 0.75, 0.45, 0.9 );

  TFile * fdata = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28605.root");

  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetTitle( "1.9 #upoint10^{15 }p/cm^{2} shallow 325 V" );
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->GetXaxis()->SetRangeUser(-0.18,0.18);
  cmsqyvsd->Draw();
  lgnd->AddEntry( cmsqyvsd, " Data 325 V ", "l" );
  int i0 = cmsqyvsd->FindBin(0.085);
  double q0 = cmsqyvsd->GetBinContent(i0);
  cout << "q0 " << q0 << endl;

  TH1F * data = (TH1F*)cmsqyvsd->Clone();

  // smeared sim:

  TFile * fsim = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V325.root" ); // mdy 

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(0.085);
  double s0 = svsy->GetBinContent(k0);
  cout << "s0 " << s0 << endl;
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, " KDetSim", "l" );

  lgnd->Draw( "same" );

  c1.Print( "s2-325.ps" );

  double sumdy = 0;
  double chisq = 0;
  int npnt = 0;

  for( int ii = 1; ii <= data->GetNbinsX(); ++ii ) {

    double y = data->GetBinCenter(ii);
    if( y < -0.16 ) continue;
    if( y >  0.16 ) break;

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
    chisq += pow( (s-q)/e, 2 );
    ++npnt;

  }

  cout
    << "mean dy " << sumdy/npnt << endl
    << "mean sy " << sqrt(chisq/npnt) << endl
    << "norm ke " << q0/s0
    << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s2-325.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s2-325.ps" );
  ierr = system( "rm -f  s2-325.ps" );
  cout << "evince s2-325.pdf" << endl;

}
