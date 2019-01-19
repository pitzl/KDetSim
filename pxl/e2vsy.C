
// Daniel Pitzl, Jan 2019
// KDetSim edg

// root -l e2vsy.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

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

  c1.Print( "evsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.6, 0.4, 0.8, 0.9 );

  // 225 V:

  TFile * fs225 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V225.root" );
  evsy->SetTitle( "E-field" );
  evsy->SetStats(0);
  evsy->GetXaxis()->SetRangeUser( -140, 140 );
  evsy->SetMinimum(0);
  evsy->SetLineColor(11); // grey
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 225 V", "l" );

  // 275 V:

  TFile * fs275 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V275.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(93); // gold
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 250 V", "l" );

  // 325 V:

  TFile * fs325 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V325.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(51); // purple
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 300 V", "l" );

  // 400 V:

  TFile * fs400 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V400.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(6);
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 400 V", "l" );

  // 500 V:

  TFile * fs500 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V500.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(4);
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 500 V", "l" );

  // 600 V:

  TFile * fs600 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V600.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(2);
  evsy->Draw( "samel" );
  lgnd->AddEntry( evsy, " 600 V", "l" );

  // 800 V:

  TFile * fs800 = TFile::Open( "edg9_px450_F2_dlt11_e10_h28_V800.root" );
  evsy->SetStats(0);
  evsy->SetLineColor(1);
  evsy->Draw("samel");
  lgnd->AddEntry( evsy, " 800 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "evsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "evsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf evsy.ps" );
  ierr = system( "rm -f  evsy.ps" );
  cout << "evince evsy.pdf" << endl;

}
