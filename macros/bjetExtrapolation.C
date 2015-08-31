void bjetExtrapolation() {

  gROOT->Reset();
  
  int lowRange = 500;
  int highRange = 5500;
  int lowRangeMHT = 200;
  int highRangeMHT = 5200;
  int nBins = 1;

  TH1F *HThistData = new TH1F("HThistData","HThistData",nBins,lowRange,highRange);
  TH1F *HThistDataObs = new TH1F("HThistDataObs","HThistDataObs",nBins,lowRange,highRange);

  TH1F *HTJet = new TH1F("HTJet","HTJet",nBins,lowRange,highRange);
  
  TH1F *h_obs = new TH1F("h_obs","h_obs",12,1,13);
  TH1F *h_sf = new TH1F("h_sf","h_sf",12,1,13);
  TH1F *h_sfBinomCorr = new TH1F("h_sfBinomCorr","h_sfBinomCorr",12,1,13);

  // for parallel processing with TProof, set proof = true
  // *****Warning, TProof seems to give inconsistent results*****
  bool proof = false;
  int binIter = 0;
  // to check closure, set doLumi = effective luminosity of sample(s)
  float doLumi = 10.;

  // choose sample:
  // 0: signal selection
  // 1: Z->mumu
  // 2: Z->ee
  // 3: Z->ll (mumu or ee)
  int sample = 3;

  float prob = 0.045 ;

  if(proof)
    TProof *p = TProof::Open("");
  
  TChain dataChain("TreeMaker2/PreSelection");
  
  if(sample == 1 || sample == 2 || sample == 3) {

    setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-100to200-RunIISpring15.root",  doLumi, 139.4);
    setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-200to400-RunIISpring15.root",  doLumi, 42.75);
    setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-400to600-RunIISpring15.root",  doLumi, 5.497);
    setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-600toInf-RunIISpring15.root",  doLumi, 2.21);
    dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-100to200-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-200to400-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-400to600-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-600toInf-RunIISpring15.root") ;
    
  } else if(sample == 0) {

    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-100to200-RunIISpring15.root",  doLumi, 280.47);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-200to400-RunIISpring15.root",  doLumi, 78.36);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-400to600-RunIISpring15.root",  doLumi, 10.94);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-600toInf-RunIISpring15.root",  doLumi, 4.20);
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-100to200-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-200to400-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-400to600-RunIISpring15.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-600toInf-RunIISpring15.root") ;

  } else {
    cout << "Picked unknown sample" << endl;
    break;
  }

  
  if(proof) 
    dataChain.SetProof();

  gROOT->ProcessLine(".L /home/ww/work/RA2b/macros/bjetExtrapolationAndRatio.C");    

  int nJetBin = 1;
  int bJetBin = -1;
  float SimpleSF[4];
  float SimpleSFerror[4];

  // statistical error on correction
  // obtained from running toys
  // small contribution to error
  float binomStatError[3][4] = {
    {0.0, 0.0, 0.0, 0.0}, // njet 4-6 no correction error
    {0.0, 0.01, 0.022, 0.04}, // njet 7-8
    {0.0, 0.027, 0.062, 0.11} // njet 9+
  };
  
  HThistData->Sumw2();  
  HThistDataObs->Sumw2();  
  gStyle->SetOptStat(0);

  float njets[13];

  // need to get the yields in each njets bin
  // for the binomial correction term inclusive
  // in other variables
  for(int i = 4; i<13; i++) {
    TCut Cuts = getCuts(i, -1, -1, sample);
    dataChain.Project("HTJet", "HT", Cuts);
    njets[i] = HTJet->Integral();
  }
  
  cout << "begin loop" << endl;
  for(int nJetBin=1; nJetBin<=3; nJetBin++) {
    for(int bJetBin=0; bJetBin<=3; bJetBin++) {

	TCut ZlepCuts = getCuts(nJetBin, bJetBin, -1, sample);

	dataChain.Project("HThistData", "HT", ZlepCuts);
	
	if(nJetBin==1) {
	  SimpleSF[bJetBin] = HThistData->Integral();
	  SimpleSFerror[bJetBin] = HThistData->GetBinError(1);
	}

	int bin = 4*(nJetBin-1)+bJetBin+1;
	h_obs->Fill(bin,HThistData->Integral());
	h_obs->SetBinError(bin,HThistData->GetBinError(1));

	float scaledPred = h_obs->GetBinContent(4*(nJetBin-1)+1)*SimpleSF[bJetBin]/SimpleSF[0];

	float scaledPredBinomCorr = scaledPred*getCorrection(bJetBin, njets, prob, nJetBin);
	
	float scaledPredError;

	scaledPredError = scaledPred*TMath::Sqrt( 1/SimpleSF[bJetBin]+1/SimpleSF[0]
						  + binomStatError[nJetBin-1][bJetBin]**2);
	  
	float scaledPredBinomCorrError = 0;

	scaledPredBinomCorrError =  scaledPredBinomCorr*TMath::Sqrt(1/SimpleSF[bJetBin]+1/SimpleSF[0]);
	
	h_sf->Fill(bin,scaledPred);
	h_sf->SetBinError(bin,scaledPredError);

	h_sfBinomCorr->Fill(bin,scaledPredBinomCorr);
	h_sfBinomCorr->SetBinError(bin,scaledPredBinomCorrError);

	// monte carlo corrected sf
	float scaledPredMCCorr = h_obs->GetBinContent(bin);

	float scaledPredMCCorrError = HThistData->GetBinError(1);

	cout << "nJet "
	     << nJetBin
	     << " btags "
	     << bJetBin << endl;
	cout << "Extrapolation Factor "
	     << scaledPredBinomCorr/h_sfBinomCorr->GetBinContent(1+(nJetBin-1)*4) << endl;
	cout << "Extrapolation Factor Error "
	     << scaledPredBinomCorrError/h_sfBinomCorr->GetBinContent(1+(nJetBin-1)*4)  << endl;
	cout << "Extrapolation Factor MC "
	     << scaledPredMCCorr/h_obs->GetBinContent(4*(nJetBin-1)+1) << endl;
	cout << "Extrapolation Factor MC Error "
	     << scaledPredMCCorrError/h_obs->GetBinContent(4*(nJetBin-1)+1) << endl;
	
    }
  }

  std::ostringstream ss;
  ss << doLumi;
  std::string L(ss.str());
  
  if (sample == 0) {
    TString title = "CMS Simulation, Z#rightarrow#nu#bar{#nu}, #sqrt{s} = 13TeV, L = "
      +L+" fb^{-1}";
    TString fileName = "Znn_L"+L+"_ExtrapolationPlot.pdf";
  } else if (sample == 1 || sample == 2 || sample == 3) {
    TString title = "CMS Simulation, Z#rightarrow ll, #sqrt{s} = 13TeV, L = "
      +L+" fb^{-1}";
    TString fileName = "Zll_L"+L+"_ExtrapolationPlot.pdf";
  } else {
    TString title = "CMS Simulation, #sqrt{s} = 13TeV, L = "
      +L+" fb^{-1}";    
  }
  
  TCanvas *canvas = bjetExtrapolationAndRatio(h_sf, h_sfBinomCorr, h_obs,title,"Bin","Events", true);
  //  canvas->SaveAs(fileName);
  
}

float getMax(float a, float b) {
  return (a<b)?b:a;   
}


float getCorrection(int b, float NJets[], float p, int bin){

  float NJets5o4 = NJets[5]/NJets[4];
  float NJets6o4 = NJets[6]/NJets[4];
  float NJets8o7 = NJets[8]/NJets[7];
  float NJets10o9 = NJets[10]/NJets[9];
  float NJets11o9 = NJets[11]/NJets[9];
  float NJets12o9 = NJets[12]/NJets[9];
  
  // normTerm same for all bins
  float normTerm = (1+NJets5o4*(1-p)+NJets6o4*(1-p)**2)/
    (TMath::Binomial(4,b)+NJets5o4*TMath::Binomial(5,b)*(1-p)+NJets6o4*TMath::Binomial(6,b)*(1-p)**2);
  
  if(bin==2)
    float corrTerm = (TMath::Binomial(7,b)+NJets8o7*TMath::Binomial(8,b)*(1-p))*normTerm/
      (1+NJets8o7*(1-p));
  else if(bin==3)
    float corrTerm = (TMath::Binomial(9,b)+NJets10o9*TMath::Binomial(10,b)*(1-p)
		      +NJets11o9*TMath::Binomial(11,b)*(1-p)**2
		      +NJets12o9*TMath::Binomial(12,b)*(1-p)**3)*normTerm/
      (1+NJets10o9*(1-p)+NJets11o9*(1-p)**2+NJets12o9*(1-p)**2);
  else
    corrTerm = 1;
  
  return corrTerm;

}

// baseline: nJetBin = -1, bJetBin = -1, kinBin = -1, Zcode = 0
// dilepton baseline: nJetBin = -1, bJetBin = -1, kinBin = -1, Zcode = 3

// 3 njet bins: nJetBin = 1, 2, or 3, (-1 for baseline) otherwise
// njet = i will cut on the ith jet multiplicity
// e.g. njet = 7 returns "NJets==7"

// 4 btag bins: bJetBin = 0, 1, 2, or 3, (-1 for baseline)

// kinematic binning: kinBin = 1, 2, 3, 4, 5, or 6, (-1 for baseline)

// sample 
TCut getCuts(int nJetBin, int bJetBin, int kinBin, int Zcode){

  // if dilepton selection, we need to use cleaned vars
  if(Zcode == 1 || Zcode == 2 || Zcode == 3) {
    char *htg = "HTclean>=%i"; char *htl = "HTclean<%i";
    char *mhtg = "MHTclean>=%i"; char *mhtl = "MHTclean<%i";
    char *njg = "NJetsclean>=%i"; char *nje = "NJetsclean==%i"; char *njl = "NJetsclean<=%i";
    char *bjg = "BTagsclean>=%i"; char *bje = "BTagsclean==%i"; char *bjl = "BTagsclean<=%i";
    char *dp1g = "DeltaPhi1clean>%i"; char *dp2g = "DeltaPhi2clean>%i"; char *dp3g = "DeltaPhi3clean>%i";
    char *dp1l = "DeltaPhi1clean<=%i"; char *dp2l = "DeltaPhi2clean<=%i"; char *dp3l = "DeltaPhi3clean<=%i";
  } else {
    char *htg = "HT>=%i"; char *htl = "HT<%i";
    char *mhtg = "MHT>=%i"; char *mhtl = "MHT<%i";
    char *njg = "NJets>=%i"; char *nje = "NJets==%i"; char *njl = "NJets<=%i";
    char *bjg = "BTags>=%i"; char *bje = "BTags==%i"; char *bjl = "BTags<=%i";
    char *dp1g = "DeltaPhi1>%i"; char *dp2g = "DeltaPhi2>%i"; char *dp3g = "DeltaPhi3>%i";
    char *dp1l = "DeltaPhi1<=%i"; char *dp2l = "DeltaPhi2<=%i"; char *dp3l = "DeltaPhi3<=%i";
  }

  TCut cuts = Form(htg,10); cuts += Form(mhtg,10);

  switch (kinBin) {
  case 1:
    cuts += Form(htl,800); cuts += Form(mhtl,500);
    break;
  case 2:
    cuts += Form(htg,800); cuts += Form(htl,1200); cuts += Form(mhtl,500);
    break;
  case 3:
    cuts += Form(htg,1200); cuts += Form(mhtl,500);
    break;
  case 4:
    cuts += Form(htl,1200); cuts += Form(mhtg,500); cuts += Form(mhtl,750);
    break;
  case 5:
    cuts += Form(htg,1200); cuts += Form(mhtg,500); cuts += Form(mhtl,750);
    break;
  case 6:    
    cuts += Form(htg,800); cuts+= Form(mhtg,750);
    break;
  default:
    cuts += Form(htg,500); cuts += Form(mhtg,200);
    break;
  }

  switch (nJetBin) {
  case -1:
    cuts += Form(njg,4);
    break;
  case 1:
    cuts += Form(njg,4); cuts += Form(njl,6);
    break;
  case 2:
    cuts += Form(njg,7); cuts += Form(njl,8);
   break;
  case 3:
    cuts += Form(njg,9);
    break;
  default:
    cuts += Form(nje,nJetBin);
    break;
 }
  
  switch (bJetBin) {
  case -1:
    cuts += Form(bjg,0);
    break;
  case 3:
    cuts += Form(bjg,3);
  default:
    cuts += Form(bje,bJetBin);
    break;
  }

  switch (Zcode) {
  case 0: // SIG
    cuts += "Leptons==0";
    cuts += "isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0";
    cuts += "DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3";
    break;
  case 1: // Z->mumu
    cuts += "(@Muons.size()==2&&@Electrons.size()==0)&&Zp4.M()>=76.188&&Zp4.M()<=106.188";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  case 2: // Z->ee
    cuts += "(@Muons.size()==0&&@Electrons.size()==2)&&Zp4.M()>=76.188&&Zp4.M()<=106.188";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  case 3: // Z->ll
    cuts += "(@Muons.size()==2&&@Electrons.size()==0)||(@Muons.size()==0&&@Electrons.size()==2)";
    cuts += "Zp4.M()>=76.188&&Zp4.M()<=106.188";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  default: // picked unknown sample
    cuts += "BTags>=99999";
    cout << "******WARNING PICKED UNKNOWN SAMPLE******" << endl;
    cout << "******Plots will be empty******" << endl;
    break;
  }
  return cuts; 
}

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
