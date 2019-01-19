
// Daniel Pitzl, Jan 2019
// KDetSim edg

// root -l s2.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.4, "y" );

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

  c1.Print( "s2-800.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.87, 0.5, 0.99, 0.9 );

  TFile * fd800 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28600.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetTitle( "1.9 #upoint10^{15 }p/cm^{2} shallow" );
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->GetXaxis()->SetRangeUser( -0.18, 0.18 );
  cmsqyvsd->Draw();
  lgnd->AddEntry( cmsqyvsd, " 800 V", "l" );
  int i0 = cmsqyvsd->FindBin(0.075);
  double q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs800 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V800.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(0.075);
  double s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 600:

  TFile * fd600 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28601.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 600 V", "l" );
  i0 = cmsqyvsd->FindBin(0.085);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs600 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V600.root" ); // mdy 
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.085);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 500:

  TFile * fd500 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28602.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 500 V", "l" );
  i0 = cmsqyvsd->FindBin(0.085);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs500 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V500.root" ); // mdy 
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.085);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 400:

  TFile * fd400 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28610.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 400 V", "l" );
  i0 = cmsqyvsd->FindBin(0.095);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs400 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V400.root" ); // mdy 
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.095);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 325:

  TFile * fd325 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28605.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 325 V", "l" );
  i0 = cmsqyvsd->FindBin(0.085);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs325 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V325.root" ); // mdy 
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.085);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 275:

  TFile * fd275 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28606.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 275 V", "l" );
  i0 = cmsqyvsd->FindBin(0.095);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs275 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V275.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.095);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  // 225:

  TFile * fd225 = TFile::Open( "/home/pitzl/eudaq/tele-scope/shallow28607.root");
  cmsqyvsd->SetStats(0);
  cmsqyvsd->SetLineColor(1);
  cmsqyvsd->SetLineWidth(3);
  cmsqyvsd->Draw("same");
  lgnd->AddEntry( cmsqyvsd, " 225 V", "l" );
  i0 = cmsqyvsd->FindBin(0.095);
  q0 = cmsqyvsd->GetBinContent(i0);

  // smeared sim:

  TFile * fs225 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V225.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.095);
  s0 = svsy->GetBinContent(k0);
  svsy->Scale(q0/s0);
  svsy->Draw( "same" );
  cout << "norm q " << q0/s0 << endl;

  lgnd->Draw( "same" );

  c1.Print( "s2.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s2.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s2.ps" );
  ierr = system( "rm -f  s2.ps" );
  cout << "evince s2.pdf" << endl;

}
