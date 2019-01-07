
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l s8.C

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

  c1.Print( "s8.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.80, 0.50, 0.95, 0.9 );

  TFile * fd800 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34437.root" );
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "8 #upoint10^{15 }n/cm^{2} edge-on" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " 800 V", "l" );
  int i0 = colpvsy->FindBin(-0.0525);
  double ph0 = colpvsy->GetBinContent(i0);
  cout << "ph0 " << ph0 << endl;

  // smeared sim:

  //TFile * fs800 = TFile::Open( "edg3_px50_F8_lam10_e4p5_h4_V800_smr11.root" ); // best
  TFile * fs800 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V800.root" ); // mdy 1.06
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(-0.0525);
  double qh0 = svsy->GetBinContent(k0);
  cout << "qh0 " << qh0 << endl;
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 600 V:

  TFile * fd600 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34430.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 600 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs600 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V600.root" );
  TFile * fs600 = TFile::Open( "edg15_px25_F8_dlt11_e2_h2_V600.root" ); // mdy 1.59
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 500 V:

  TFile * fd500 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34431.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 500 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs500 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V500.root" );
  TFile * fs500 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V500.root" ); // mean dy 1.09
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 400 V:

  TFile * fd400 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34432.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 400 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs400 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400.root" );
  TFile * fs400 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V400.root" ); // mdy 1.41
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 300 V:

  TFile * fd300 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34433.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 300 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs300 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V300.root" );
  TFile * fs300 = TFile::Open( "edg15_px25_F8_dlt11_e2_h3_V300.root" ); // mdy 1.18
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 250 V:

  TFile * fd250 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34434.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 250 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs250 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V250.root" );
  TFile * fs250 = TFile::Open( "edg15_px25_F8_dlt11_e3_h3_V250.root" ); // mean dy 1.52
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 200 V:

  TFile * fd200 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34435.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 200 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs200 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V200.root" );
  TFile * fs200 = TFile::Open( "edg15_px25_F8_dlt11_e3_h4_V200.root" ); // mdy 1.02
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 150 V:

  TFile * fd150 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34436.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 150 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  //TFile * fs150 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V150.root" );
  TFile * fs150 = TFile::Open( "edg15_px25_F8_dlt11_e4_h3_V150.root" ); // mdy 1.10
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  lgnd->Draw( "same" );

  c1.Print( "s8.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s8.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s8.ps" );
  ierr = system( "rm -f  s8.ps" );
  cout << "evince s8.pdf" << endl;

}
