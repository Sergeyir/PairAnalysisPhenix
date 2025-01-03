struct
{
	std::string runName = "Run14HeAu200";
	TFile *dataFile = TFile::Open(("data/Real/" + runName + "/SingleTrack/sum.root").c_str());
	TFile *simFile = TFile::Open(("data/PostSim/" + runName + "/Heatmaps/all.root").c_str());
} Par;

void ELossTOFe()
{
   TH2F *eloss = (TH2F *) Par.dataFile->Get("ELoss: TOFe");
   eloss->Draw("COLZ");

   gPad->SetLogz();

   TF1 cutFunc("ELoss cut", "[0]*pow(x, [1])");
   cutFunc.SetParameters(0.0005, -2.5);
   cutFunc.SetLineWidth(2);
   cutFunc.DrawClone("SAME");
}
