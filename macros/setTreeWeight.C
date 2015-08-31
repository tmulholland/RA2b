// use this if skimming
// check to make sure WeightHist properly
// stored the samples cross-section 
setTreeWeight(char* file,  float L) {

  TFile f(file,"update");
  tree = (TTree*)f.Get("TreeMaker2/PreSelection");
  
  TH1F *WeightHist = new TH1F("WeightHist","WeightHist",1,-10000,10000);
  tree->Project("WeightHist","Weight","HT>-10");

  // weight gives 4/fb 
  float weight = WeightHist->GetMean();

  cout << "sample lumi = " << 4./weight << endl;
  cout << "new tree lumi = " << L << endl;

  weight *= L/4.;

  tree->SetWeight(weight);
  tree->AutoSave(); 
  
}

// use this if you didn't run over the entire miniAOD
// or if WeightHist doesn't have the correct
// weight applied
setTreeWeight(char* file,  float L, float sigma) {

  TFile f(file,"update");
  tree = (TTree*)f.Get("TreeMaker2/PreSelection");
  
  gDirectory->cd("TreeMaker2");

  // cross section should be given in pb
  // 1 pb = 1000 fb
  sigma *= 1000.;
  float SampleLumi = (PreSelection->GetEntries()/sigma);
  cout << "sample lumi = " << SampleLumi << endl;
  cout << "new tree lumi = " << L << endl;

  float weight = L/SampleLumi;

  tree->SetWeight(weight);
  tree->AutoSave(); 
  
}
