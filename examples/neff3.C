
// doping profile: 3rd order polynomial
// root -l cv.C neff3.C
// f3->GetParameter(1)
// f3->SetParameter(1,-2)
// f3->Draw()

{
  TF1 * f3 = new TF1( "f3",
		      "[0]+[1]*(x[0]-[2])/150+[3]*pow((x[0]-[2])/150,3)",
		      0, 300 );
  f3->SetParameter( 0,   0 ); // mid doping
  f3->SetParameter( 1,  -8 ); // slope
  f3->SetParameter( 2, 135 ); // x mid
  f3->SetParameter( 3,  -4 ); // X**3
  f3->Draw("l");
}

