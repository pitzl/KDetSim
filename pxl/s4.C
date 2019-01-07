
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l s4.C

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

  c1.Print( "s4.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.80, 0.50, 0.95, 0.9 );

  TFile * fd800 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34453.root" );
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "4 #upoint10^{15 }n/cm^{2} edge-on" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " 800 V", "l" );
  int i0 = colpvsy->FindBin(-0.0425);
  double ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs800 = TFile::Open( "edg3_px50_F4_lam10_e5_h6_V800.root" ); // best
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(-0.0425);
  double qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 600 V:

  TFile * fd600 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34443.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 600 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs600 = TFile::Open( "edg3_px50_F4_lam10_e4_h4_V600.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 500 V:

  TFile * fd500 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34444.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 500 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs500 = TFile::Open( "edg3_px50_F4_lam10_e4_h4_V500.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 400 V:

  TFile * fd400 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34445.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 400 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs400 = TFile::Open( "edg3_px50_F4_lam11_e4_h4_V400.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 300 V:

  TFile * fd300 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34446.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 300 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs300 = TFile::Open( "edg3_px50_F4_lam12_e4_h4_V300.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 250 V:

  TFile * fd250 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34447.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 250 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs250 = TFile::Open( "edg3_px50_F4_lam12_e4_h4_V250.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 200 V:

  TFile * fd200 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34449.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 200 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs200 = TFile::Open( "edg3_px50_F4_lam12_e4_h4_V200.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 150 V:

  TFile * fd150 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34450.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 150 V", "l" );
  i0 = colpvsy->FindBin(-0.0525);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs150 = TFile::Open( "edg3_px50_F4_lam15_e4_h4_V150.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0525);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  // 100 V:

  TFile * fd100 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34451.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 100 V", "l" );
  i0 = colpvsy->FindBin(-0.0575);
  ph0 = colpvsy->GetBinContent(i0);

  TFile * fs100 = TFile::Open( "edg3_px50_F4_lam16_e4_h4_V100.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(-0.0575);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );

  lgnd->Draw( "same" );

  c1.Print( "s4.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s4.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s4.ps" );
  ierr = system( "rm -f  s4.ps" );
  cout << "evince s4.pdf" << endl;

}
