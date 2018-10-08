int PlotQ2X(){

  TChain * ResultTree=new TChain("ResultTree");
  TString pattern;
  // pattern = "Results/Pieces/Breit_pythia.ep.20x250.1Mevents.RadCor=0.Q2=10*";
  pattern = "Results/Pieces/Breit_pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0*";
  
  ResultTree->Add( pattern );
  ResultTree->Draw("Q2:TMath::Log10(X)","","colz");

  // gPad->SaveAs("plots/Q2X10-100.png");
  gPad->SaveAs("plots/Q2X1-10.png");

  return 0;
}
