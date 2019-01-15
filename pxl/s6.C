
// Daniel Pitzl, Jan 2019
// KDetSim edg

// root -l s6.C

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

  c1.Print( "s6-800.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.85, 0.5, 0.99, 0.9 );

  TFile * fd800 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33592.root" ); // 174i long
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "6.6 #upoint10^{15 }p/cm^{2} edge-on" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " 800 V", "l" );
  int i0 = colpvsy->FindBin(0.0475);
  double ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs800 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V800.root" ); // mob4: mean dy 2.05
  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(0.0475);
  double qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 700:

  TFile * fd700 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33581.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 700 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs700 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V700.root" ); // mdy 2.07

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 600:

  TFile * fd600 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33582.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 600 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs600 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V600.root" ); // mdy 1.15

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 500:

  TFile * fd500 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33583.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 500 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs500 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V500.root" ); // mob4 mdy 

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 450:

  TFile * fd450 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33584.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 450 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs450 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V450.root" ); // mean dy 2.00

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 400:

  TFile * fd400 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33585.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 400 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs400 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V400.root" ); // mdy 1.93

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 350:

  TFile * fd350 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33586.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 350 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs350 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V350.root" ); // mean dy

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 300:

  TFile * fd300 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33587.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 300 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs300 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V300.root" ); // mean dy 1.71

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 250:

  TFile * fd250 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33588.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 250 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs250 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V250.root" ); // mean dy 1.97

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 200:

  TFile * fd200 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33589.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 200 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs200 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V200.root" ); // mean dy 1.07

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  // 150:

  TFile * fd150 = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33590.root" );
  colpvsy->SetStats(0);
  colpvsy->SetLineColor(1);
  colpvsy->Draw("same");
  lgnd->AddEntry( colpvsy, " 150 V", "l" );
  i0 = colpvsy->FindBin(0.0475);
  ph0 = colpvsy->GetBinContent(i0);

  // smeared sim:

  TFile * fs150 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V150.root" ); // mean dy 1.24

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  k0 = svsy->FindBin(0.0475);
  qh0 = svsy->GetBinContent(k0);
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  cout << "norm ph " << ph0/qh0 << endl;

  lgnd->Draw( "same" );

  c1.Print( "s6.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s6.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s6.ps" );
  ierr = system( "rm -f  s6.ps" );
  cout << "evince s6.pdf" << endl;

}
