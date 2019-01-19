
// Daniel Pitzl, Jan 2019
// KDetSim vert

// root -l qvsf.C

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

  c1.Print( "qvsf.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.6, 0.7, 0.93, 0.88 );

  // 285 um:

  TFile * fs285 = TFile::Open( "vert2_pz284_dlt11_V800.root" );
  qvsf->SetTitle( "charge vs p fluence" );
  qvsf->GetYaxis()->SetTitle( "collected charge [ke]" );
  qvsf->SetStats(0);
  qvsf->SetMinimum(0);
  qvsf->SetLineColor(2);
  qvsf->Draw( "samel" );
  lgnd->AddEntry( qvsf, " 285 #mum, 800 V", "l" );

  // 275 V:

  TFile * fs150 = TFile::Open( "vert6_pz154_dlt9_V800.root" );
  qvsf->SetStats(0);
  qvsf->SetLineColor(4);
  qvsf->Draw( "samel" );
  lgnd->AddEntry( qvsf, " 150 #mum, 800 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "qvsf.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "qvsf.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf qvsf.ps" );
  ierr = system( "rm -f  qvsf.ps" );
  cout << "evince qvsf.pdf" << endl;

}
