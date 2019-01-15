
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l s6-600.C

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

  c1.Print( "s6-600.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.20, 0.75, 0.45, 0.9 );

  TFile * fdata = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33582.root" ); // mdy 1.54
  //TFile * fdata = TFile::Open( "/home/pitzl/eudaq/tele-scope/scoped33799.root" );
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "6.6 #upoint10^{15 }p/cm^{2} edge-on 600 V" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " Data 600 V ", "l" );
  int i0 = colpvsy->FindBin(0.0475);
  double ph0 = colpvsy->GetBinContent(i0);
  cout << "ph0 " << ph0 << endl;

  TH1F * data = (TH1F*)colpvsy->Clone();

  // smeared sim:

  //TFile * fsim = TFile::Open( "edg3_px50_F6_dlt9_e6_h14_V600.root" ); // mean dy 2.5

  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e7_h14_V600.root" ); // mdy 1.32

  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e11_h21_V600.root" ); // mdy 1.29
  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e12_h21_V600.root" ); // mdy 1.22
  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e12_h14_V600.root" ); // mdy 1.43
  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e12_h28_V600.root" ); // mdy 1.15
  TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e11_h28_V600.root" ); // mob4, mdy 1.04
  //TFile * fsim = TFile::Open( "edg9_px50_F6_dlt9_e8_h28_V600.root" ); // mob3, mdy 1.15

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(0.0475);
  double qh0 = svsy->GetBinContent(k0);
  cout << "qh0 " << qh0 << endl;
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, " KDetSim", "l" );

  lgnd->Draw( "same" );

  c1.Print( "s6-600.ps" );

  double sumdy = 0;
  double chisq = 0;
  int npnt = 0;

  for( int ii = 1; ii <= data->GetNbinsX(); ++ii ) {

    double y = data->GetBinCenter(ii);
    if( y < -0.08 ) continue;
    if( y >  0.08 ) break;

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
    << "norm ph " << ph0/qh0
    << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s6-600.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s6-600.ps" );
  ierr = system( "rm -f  s6-600.ps" );
  cout << "evince s6-600.pdf" << endl;

}
