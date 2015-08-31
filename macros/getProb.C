// Algebraic way of determining the probability parameter in the correction.
// Gives consistent results as when doing fits to the binomial.

void getProb() {

  gROOT->Reset();
  
  int lowRange = 500;
  int highRange = 5500;
  int lowRangeMHT = 200;
  int highRangeMHT = 5200;
  int nBins = 1;


  TH1F *HTJet = new TH1F("HTJet","HTJet",nBins,lowRange,highRange);

  TH1F *probPlot = new TH1F("probPlot", " ", 12, 1., 13.);
  
  bool proof = false; 
  int binIter = 0;

  // set doLumi to the lumiosity that you'd like
  // to know the statistical error on p at 
  float doLumi = 10.;
  //  float doLumi = 1573.74;// TTZJets
  //  float doLumi = 382.; //DYJets
  //  float doLumi = 713.75; //Phys14 Zinv
  int sample = 3;

  gROOT->ProcessLine(".L getCuts.C");    
  gROOT->ProcessLine(".L setTreeWeight.C");    

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

    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-100to200-Phys14.root",  doLumi, 473.2);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-200to400-Phys14.root",  doLumi, 128.0);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-400to600-Phys14.root",  doLumi, 15.23);
    setTreeWeight("/home/ww/work/data/RA2bTrees/ZJets-600toInf-Phys14.root",  doLumi, 5.224);
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-100to200-Phys14.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-200to400-Phys14.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-400to600-Phys14.root") ;
    dataChain.Add("/home/ww/work/data/RA2bTrees/ZJets-600toInf-Phys14.root") ;

  } else {
    cout << "Picked unknown sample" << endl;
    break;
  }

  if(proof) 
    dataChain.SetProof();

  gStyle->SetOptStat(0);

  float nb0[13];
  //  float nb0error[13];
  float nbi[13];
  //  float nbierror[13];

  for(int i = 4; i<13; i++) {
    TCut Cuts = getCuts(i, 0, -1, sample);
    dataChain.Project("HTJet", "HT", Cuts);
    nb0[i] = HTJet->Integral();
  }
  
  for(int i = 4; i<13; i++) {
    TCut Cuts = getCuts(i, -1, -1, sample);
    dataChain.Project("HTJet", "HT", Cuts);
    nbi[i] = HTJet->Integral();
  }

  float totalProbNum = 0.;
  float totalProbDenom = 0.;
  for(int i = 4; i<13; i++) {

    if ( nb0[i] == 0 || nbi[i] == 0)
      break;
    
    float nj = i; //to avoid integer crap 
    cout << "NJet " << nj << endl;
    cout << "fraction of 0b events " << nb0[i]/nbi[i] << " +- " << TMath::Sqrt(1/nb0[i]+1/nbi[i]) << endl;

    float ProbError = TMath::Sqrt((1/nbi[i])**(2/nj)*(1/nj)**2*nb0[i]**(2*(1/nj-1))*nb0[i]
				  +nb0[i]**(2/nj)*1/nj**2*nbi[i]**(-2*(1/nj+1))*nbi[i]);
    float Prob = 1-(nb0[i]/nbi[i])**(1/nj);
    
    cout << "probability of tagging a bjet " << Prob << " +- " << ProbError << endl;
    probPlot->Fill(i, Prob);
    probPlot->SetBinError(i,ProbError);

    float weight = 1/ProbError**2;
    totalProbNum += weight*Prob;
    totalProbDenom += weight;

  }
  probPlot->SetMinimum(0.);

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,1000,600);
  
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.125);
  c1->SetRightMargin(0.125);

  probPlot->GetXaxis()->SetTitle("NJets");
  probPlot->GetYaxis()->SetTitle("Probability");
  probPlot->GetXaxis()->SetTitleSize(0.06);
  probPlot->GetXaxis()->SetLabelSize(0.05);
  probPlot->GetXaxis()->SetTitleOffset(0.80);
  probPlot->GetYaxis()->SetLabelSize(0.05);
  probPlot->GetYaxis()->SetTitleSize(0.06);
  probPlot->GetYaxis()->SetTitleOffset(0.80);

  probPlot->Draw();

  cout <<"**********************"<<endl;
  cout <<"**********************"<<endl;
  cout << "Weight average probability = " << totalProbNum/totalProbDenom << " +- " << 1/TMath::Sqrt(totalProbDenom) << endl;
  
}
