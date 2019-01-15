
// Daniel Pitzl, Jan 2019
// KDetSim edg

// root -l v6vsy.C

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

  c1.Print( "vvsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.88, 0.16, 0.99, 0.67 );

  // 150 V:

  TFile * fs150 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V150.root" ); // mean dy 1.24
  vvsy->SetTitle( "e velocity" );
  vvsy->SetStats(0);
  vvsy->GetXaxis()->SetRangeUser(-75,75);
  vvsy->SetNdivisions(-506);
  //vvsy->SetMaximum(25);
  vvsy->SetLineColor(95); // orange
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 150 V", "l" );

  // 200 V:

  TFile * fs200 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V200.root" ); // mean dy 1.07
  vvsy->SetStats(0);
  vvsy->SetLineColor(28); // brown
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 200 V", "l" );

  // 250 V:

  TFile * fs250 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V250.root" ); // mean dy 1.97
  vvsy->SetStats(0);
  vvsy->SetLineColor(51); // purple
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 250 V", "l" );

  // 300 V:

  TFile * fs300 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V300.root" ); // mean dy 1.71
  vvsy->SetStats(0);
  vvsy->SetLineColor(11); // grey
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 300 V", "l" );

  // 350 V:

  TFile * fs350 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V350.root" ); // mean dy
  vvsy->SetStats(0);
  vvsy->SetLineColor(6); // magenta
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 350 V", "l" );

  // 400 V:

  TFile * fs400 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V400.root" ); // mdy 1.93
  vvsy->SetStats(0);
  vvsy->SetLineColor(7);
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 400 V", "l" );

  // 450 V:

  TFile * fs450 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V450.root" ); // mean dy 2.00
  vvsy->SetStats(0);
  vvsy->SetLineColor(91); // gold
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 450 V", "l" );

  // 500 V:

  TFile * fs500 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V500.root" ); // mob4 mdy 
  vvsy->SetStats(0);
  vvsy->SetLineColor(8);
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 500 V", "l" );

  // 600 V:

  TFile * fs600 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V600.root" ); // mdy 1.15
  vvsy->SetStats(0);
  vvsy->SetLineColor(4);
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 600 V", "l" );

  // 700 V:

  TFile * fs700 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V700.root" ); // mdy 2.07
  vvsy->SetStats(0);
  vvsy->SetLineColor(2);
  vvsy->Draw( "samel" );
  lgnd->AddEntry( vvsy, " 700 V", "l" );

  // 800 V:

  TFile * fs800 = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V800.root" ); // mob4: mean dy 2.05
  vvsy->SetStats(0);
  vvsy->SetLineColor(1);
  vvsy->Draw("samel");
  lgnd->AddEntry( vvsy, " 800 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "vvsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "vvsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf vvsy.ps" );
  ierr = system( "rm -f  vvsy.ps" );
  cout << "evince vvsy.pdf" << endl;

}
