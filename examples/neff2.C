
// doping profile: 3rd order polynomial
// root -l cv.C neff2.C
// f2->GetParameter(1)
// f2->SetParameter(1,-2)
// f2->Draw()

{
  TF1 * f2 = new TF1( "f2",
		      "[0]+[2]*pow((x[0]-[1])/300,2)",
		      0, 300 );
  f2->SetParameter( 0,  -2 ); // left
  f2->SetParameter( 1, 300 ); // mid
  f2->SetParameter( 2,   6 ); // para
  f2->Draw("l");
}

